/*
 * atlas_ribbon_fit_tool.c
 *
 * Fit a complete scroll ribbon to a frozen atlas solution.  The input U,V
 * coordinates are observations, not optimization variables: this first pass
 * reconstructs geometry at the discrete ordering already accepted by the
 * atlas solver.  F0 exports exact section evidence, F1 fits positions, and F2
 * adds increasing-U tangent integration.
 */

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../unroll/piece_set.h"
#include "../whole/atlas_overlap_audit.h"
#include "../whole/atlas_ribbon_texture.h"
#include "../unroll/scaffold.h"
#include "../whole/atlas_ribbon_fit.h"
#include "../whole/atlas_solution.h"

#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ARFT_PATH_CAP 2048

static const char *mode_name(AtlasRibbonFitMode mode)
{
    switch (mode) {
    case ATLAS_RIBBON_OBSERVATIONS: return "f0_observations";
    case ATLAS_RIBBON_POSITION:     return "f1_position";
    case ATLAS_RIBBON_TANGENT:      return "f2_tangent";
    case ATLAS_RIBBON_REGISTER_ONLY: return "f3_u_register";
    case ATLAS_RIBBON_REGISTER_FIELD_ONLY: return "f4_u_field";
    case ATLAS_RIBBON_REGISTER_FIELD_SMOOTH_ONLY: return "f5_u_field_smooth";
    default:                        return "unknown";
    case ATLAS_RIBBON_REGISTER_COLLISION_ONLY: return "f6_u_collision";
    case ATLAS_RIBBON_REGISTERED_RIBBON: return "f7_registered_ribbon";
    }
}

static FILE *open_out(const char *dir, const char *name)
{
    char path[ARFT_PATH_CAP];
    int n = snprintf(path, sizeof path, "%s/%s", dir, name);
    if (n < 0 || (size_t)n >= sizeof path) return NULL;
    if (ves_ensure_parent_dir(path) != 0) return NULL;
    return fopen(path, "wb");
}

typedef struct {
    int32_t face0;
    int32_t face1;
    int32_t allowed;
    int32_t reason;
} ArftResidualRef;


typedef struct {
    AtlasOverlapAudit initial;
    AtlasOverlapAudit registered;
    double *registered_min_xyz;
    uint8_t *registered_allowed;
    int32_t *registered_reason;
    int32_t *face_source;
    ArftResidualRef *residual;
    size_t nkept_faces;
    size_t nresiduals;
    size_t persisted_pairs;
    size_t new_pairs;
    size_t resolved_pairs;
    size_t registered_known_pairs;
    size_t registered_allowed_pairs;
    size_t registered_direct_pairs;
    size_t registered_same_winding_local_pairs;
    size_t registered_self_pairs;
    size_t registered_hard_pairs;
    size_t registered_unclassified_pairs;
    AtlasRibbonCollisionBound *collision_bound;
    size_t ncollision_bounds;
    size_t collision_pairs_without_row;
    size_t collision_bound_face_pairs;
    size_t collision_interval_failures;
} ArftUvAudit;

typedef struct {
    int round;
    size_t target_pairs;
    size_t audit_checks;
    size_t hard_pairs_peak;
    size_t ledger_updates;
    double residual_rms_before;
    double residual_p95_before;
    double residual_max_before;
    double residual_rms_after;
    double residual_p95_after;
    double residual_max_after;
    double step_rms;
    double step_max;
    double cumulative_rms;
    double cumulative_max;
} ArftTextureAlignRound;

typedef struct {
    int rounds_requested;
    int rounds_completed;
    int sweeps;
    double pair_weight;
    double anchor_weight;
} ArftTextureAlignStats;

static int compare_residual_ref(const void *pa, const void *pb)
{
    const ArftResidualRef *a = (const ArftResidualRef *)pa;
    const ArftResidualRef *b = (const ArftResidualRef *)pb;
    if (a->face0 != b->face0) return a->face0 < b->face0 ? -1 : 1;
    return a->face1 < b->face1 ? -1 : (a->face1 > b->face1 ? 1 : 0);
}

static int compare_overlap_pair_key(const AtlasOverlapPair *a,
                                    const AtlasOverlapPair *b)
{
    if (a->component0 != b->component0)
        return a->component0 < b->component0 ? -1 : 1;
    if (a->component1 != b->component1)
        return a->component1 < b->component1 ? -1 : 1;
    if (a->turn != b->turn) return a->turn < b->turn ? -1 : 1;
    if (a->face0 != b->face0) return a->face0 < b->face0 ? -1 : 1;
    return a->face1 < b->face1 ? -1 : (a->face1 > b->face1 ? 1 : 0);
}

static const ArftResidualRef *find_residual_ref(
    const ArftUvAudit *audit, int32_t face0, int32_t face1)
{
    if (face0 > face1) {
        int32_t temporary = face0; face0 = face1; face1 = temporary;
    }
    ArftResidualRef key;
    memset(&key, 0, sizeof key);
    key.face0 = face0;
    key.face1 = face1;
    return (const ArftResidualRef *)bsearch(
        &key, audit->residual, audit->nresiduals,
        sizeof(*audit->residual), compare_residual_ref);
}
static double face_min_xyz(const PieceSet *ps, const int32_t *faces,
                           size_t face0, size_t face1)
{
    double best2 = DBL_MAX;
    for (int a = 0; a < 3; a++) {
        int32_t va = faces[face0 * 3 + (size_t)a];
        for (int b = 0; b < 3; b++) {
            int32_t vb = faces[face1 * 3 + (size_t)b];
            double distance2 = 0.0;
            for (int d = 0; d < 3; d++) {
                double delta = (double)ps->verts[(size_t)va * 3 + (size_t)d] -
                               (double)ps->verts[(size_t)vb * 3 + (size_t)d];
                distance2 += delta * delta;
            }
            if (distance2 < best2) best2 = distance2;
        }
    }
    return sqrt(best2);
}

static double layer_distance(const AtlasRibbonLayerSample *a,
                             const AtlasRibbonLayerSample *b)
{
    return hypot(a->p[0] - b->p[0], a->p[1] - b->p[1]);
}


static int compare_collision_bound(const void *pa, const void *pb)
{
    const AtlasRibbonCollisionBound *a = (const AtlasRibbonCollisionBound *)pa;
    const AtlasRibbonCollisionBound *b = (const AtlasRibbonCollisionBound *)pb;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->rank_lo != b->rank_lo) return a->rank_lo < b->rank_lo ? -1 : 1;
    return a->rank_hi < b->rank_hi ? -1 : (a->rank_hi > b->rank_hi ? 1 : 0);
}

static const AtlasRibbonCollisionBound *find_collision_bound(
    const AtlasRibbonCollisionBound *ledger, size_t nledger,
    const AtlasRibbonCollisionBound *key)
{
    return (const AtlasRibbonCollisionBound *)bsearch(
        key, ledger, nledger, sizeof(*ledger), compare_collision_bound);
}

static int merge_collision_ledger(
    AtlasRibbonCollisionBound **io_ledger,
    size_t *io_count,
    size_t *io_capacity,
    const AtlasRibbonCollisionBound *proposal,
    size_t nproposal,
    double relaxation,
    size_t *out_updates)
{
    if (io_ledger == NULL || io_count == NULL || io_capacity == NULL ||
        out_updates == NULL || !isfinite(relaxation) || relaxation <= 0.0 ||
        relaxation > 1.0 || (nproposal > 0 && proposal == NULL) ||
        *io_count > *io_capacity || nproposal > SIZE_MAX - *io_count)
        return -1;
    AtlasRibbonCollisionBound *ledger = *io_ledger;
    size_t count = *io_count, updates = 0;
    for (size_t i = 0; i < nproposal; i++) {
        AtlasRibbonCollisionBound candidate = proposal[i];
        candidate.lower_shift_delta = candidate.current_delta +
                                      relaxation * candidate.escape_delta;
        if (!isfinite(candidate.lower_shift_delta)) return -1;
        const AtlasRibbonCollisionBound *known =
            find_collision_bound(ledger, count, &candidate);
        if (known == NULL ||
            candidate.lower_shift_delta > known->lower_shift_delta + 1.0e-9)
            updates++;
    }
    size_t needed = count + nproposal;
    if (needed > *io_capacity) {
        size_t capacity = *io_capacity > 0 ? *io_capacity : 256;
        while (capacity < needed) {
            if (capacity > SIZE_MAX / 2) { capacity = needed; break; }
            capacity *= 2;
        }
        if (capacity > SIZE_MAX / sizeof(*ledger)) return -1;
        AtlasRibbonCollisionBound *resized =
            (AtlasRibbonCollisionBound *)realloc(
                ledger, capacity * sizeof(*ledger));
        if (resized == NULL) return -1;
        ledger = resized;
        *io_capacity = capacity;
    }
    for (size_t i = 0; i < nproposal; i++) {
        ledger[count] = proposal[i];
        ledger[count].lower_shift_delta = ledger[count].current_delta +
            relaxation * ledger[count].escape_delta;
        count++;
    }
    qsort(ledger, count, sizeof(*ledger), compare_collision_bound);
    size_t write = 0;
    for (size_t i = 0; i < count; i++) {
        AtlasRibbonCollisionBound *src = &ledger[i];
        if (write == 0 || compare_collision_bound(&ledger[write - 1], src)) {
            ledger[write++] = *src;
            continue;
        }
        AtlasRibbonCollisionBound *dst = &ledger[write - 1];
        uint32_t face_pairs = dst->face_pairs > src->face_pairs
                            ? dst->face_pairs : src->face_pairs;
        double min_xyz = fmin(dst->min_xyz, src->min_xyz);
        if (src->lower_shift_delta > dst->lower_shift_delta) *dst = *src;
        dst->face_pairs = face_pairs;
        dst->min_xyz = min_xyz;
    }
    *io_ledger = ledger;
    *io_count = write;
    *out_updates = updates;
    return 0;
}

static int triangle_u_slice(const int32_t *faces, const double *u,
                            const double *v, size_t face, double axial,
                            double *out_min, double *out_max)
{
    double minimum = DBL_MAX, maximum = -DBL_MAX;
    int intersections = 0;
    double scale = fmax(1.0, fabs(axial));
    for (int edge = 0; edge < 3; edge++) {
        int32_t va = faces[face * 3 + (size_t)edge];
        int32_t vb = faces[face * 3 + (size_t)((edge + 1) % 3)];
        double aa = v[va], ab = v[vb];
        double ua = u[va], ub = u[vb];
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
        double value = ua + t * (ub - ua);
        if (value < minimum) minimum = value;
        if (value > maximum) maximum = value;
        intersections++;
    }
    if (intersections == 0 || !isfinite(minimum) || !isfinite(maximum) ||
        minimum > maximum)
        return -1;
    *out_min = minimum;
    *out_max = maximum;
    return 0;
}

static int face_pair_translation_interval(
    const int32_t *faces, const double *u, const double *v,
    size_t face0, size_t face1, double *out_left, double *out_right,
    double *out_v)
{
    double min_v0 = DBL_MAX, max_v0 = -DBL_MAX;
    double min_v1 = DBL_MAX, max_v1 = -DBL_MAX;
    double candidate[8];
    size_t ncandidate = 0;
    for (int k = 0; k < 3; k++) {
        double value0 = v[faces[face0 * 3 + (size_t)k]];
        double value1 = v[faces[face1 * 3 + (size_t)k]];
        if (value0 < min_v0) min_v0 = value0;
        if (value0 > max_v0) max_v0 = value0;
        if (value1 < min_v1) min_v1 = value1;
        if (value1 > max_v1) max_v1 = value1;
    }
    double common_min = fmax(min_v0, min_v1);
    double common_max = fmin(max_v0, max_v1);
    double scale = fmax(1.0, fmax(fabs(common_min), fabs(common_max)));
    if (common_min > common_max + 1.0e-12 * scale) return -1;
    candidate[ncandidate++] = common_min;
    candidate[ncandidate++] = common_max;
    for (int k = 0; k < 3; k++) {
        double value0 = v[faces[face0 * 3 + (size_t)k]];
        double value1 = v[faces[face1 * 3 + (size_t)k]];
        if (value0 > common_min && value0 < common_max)
            candidate[ncandidate++] = value0;
        if (value1 > common_min && value1 < common_max)
            candidate[ncandidate++] = value1;
    }
    double left = DBL_MAX, right = -DBL_MAX;
    for (size_t i = 0; i < ncandidate; i++) {
        double min0, max0, min1, max1;
        if (triangle_u_slice(faces, u, v, face0, candidate[i],
                             &min0, &max0) != 0 ||
            triangle_u_slice(faces, u, v, face1, candidate[i],
                             &min1, &max1) != 0)
            return -1;
        double local_left = min0 - max1;
        double local_right = max0 - min1;
        if (local_left < left) left = local_left;
        if (local_right > right) right = local_right;
    }
    if (!isfinite(left) || !isfinite(right) || !(left < right)) return -1;
    *out_left = left;
    *out_right = right;
    *out_v = 0.5 * (common_min + common_max);
    return 0;
}

static int build_dense_chart_shift(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    double **out_dense)
{
    if (fit->row_chart_shift == NULL || fit->nrows != set->nrows ||
        fit->nrows > SIZE_MAX / solution->ncharts)
        return -1;
    size_t n = fit->nrows * solution->ncharts;
    uint8_t *present = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    double *dense = (double *)ARENA_CALLOC(arena, n, sizeof(*dense));
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        if (s->row < 0 || (size_t)s->row >= fit->nrows || s->chart < 0 ||
            (size_t)s->chart >= solution->ncharts)
            return -1;
        present[(size_t)s->row * solution->ncharts + (size_t)s->chart] = 1;
    }
    for (size_t chart = 0; chart < solution->ncharts; chart++) {
        size_t first = 0;
        while (first < fit->nrows &&
               !present[first * solution->ncharts + chart])
            first++;
        if (first == fit->nrows) continue;
        double first_shift =
            fit->row_chart_shift[first * solution->ncharts + chart];
        for (size_t row = 0; row <= first; row++)
            dense[row * solution->ncharts + chart] = first_shift;
        size_t previous = first;
        for (size_t next = first + 1; next < fit->nrows; next++) {
            if (!present[next * solution->ncharts + chart]) continue;
            double a = fit->row_chart_shift[
                previous * solution->ncharts + chart];
            double b = fit->row_chart_shift[next * solution->ncharts + chart];
            double span = (double)(next - previous);
            for (size_t row = previous; row <= next; row++) {
                double t = (double)(row - previous) / span;
                dense[row * solution->ncharts + chart] =
                    (1.0 - t) * a + t * b;
            }
            previous = next;
        }
        double last_shift =
            fit->row_chart_shift[previous * solution->ncharts + chart];
        for (size_t row = previous; row < fit->nrows; row++)
            dense[row * solution->ncharts + chart] = last_shift;
    }
    *out_dense = dense;
    return 0;
}

static int nearest_common_chart_row(
    const uint8_t *present, size_t nrows, size_t ncharts,
    int32_t chart0, int32_t chart1, double rowf, int32_t *out_row)
{
    if (present == NULL || out_row == NULL || nrows == 0 || chart0 < 0 ||
        chart1 < 0 || (size_t)chart0 >= ncharts || (size_t)chart1 >= ncharts)
        return -1;
    size_t center;
    if (rowf <= 0.0) center = 0;
    else if (rowf >= (double)(nrows - 1)) center = nrows - 1;
    else center = (size_t)floor(rowf + 0.5);
    for (size_t offset = 0; offset < nrows; offset++) {
        if (offset <= center) {
            size_t row = center - offset;
            if (present[row * ncharts + (size_t)chart0] &&
                present[row * ncharts + (size_t)chart1]) {
                *out_row = (int32_t)row;
                return 0;
            }
        }
        if (offset > 0 && center <= SIZE_MAX - offset) {
            size_t row = center + offset;
            if (row < nrows &&
                present[row * ncharts + (size_t)chart0] &&
                present[row * ncharts + (size_t)chart1]) {
                *out_row = (int32_t)row;
                return 0;
            }
        }
    }
    return 1;
}

static int build_collision_bounds(
    Arena_T arena,
    const int32_t *faces,
    const double *registered_parameter,
    const double *dense_shift,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonFitResult *fit,
    ArftUvAudit *out)
{
    size_t ncharts = solution->ncharts;
    if (faces == NULL || registered_parameter == NULL || dense_shift == NULL ||
        ncharts == 0 || set->nrows != fit->nrows ||
        set->nrows > SIZE_MAX / ncharts)
        return -1;
    size_t field_size = set->nrows * ncharts;
    uint8_t *present = (uint8_t *)ARENA_CALLOC(arena, field_size, 1);
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *sample = &set->layer_sample[i];
        if (sample->row < 0 || (size_t)sample->row >= set->nrows ||
            sample->chart < 0 || (size_t)sample->chart >= ncharts)
            return -1;
        present[(size_t)sample->row * ncharts + (size_t)sample->chart] = 1;
    }
    size_t capacity = out->registered_hard_pairs;
    out->collision_bound = (AtlasRibbonCollisionBound *)ARENA_ALLOC(
        arena, (capacity ? capacity : 1) * sizeof(*out->collision_bound));
    size_t raw_count = 0;
    double clearance = opts->observation_u_spacing > 0.0
                     ? opts->observation_u_spacing : 1.0e-6;
    for (size_t i = 0; i < out->registered.npairs; i++) {
        if (out->registered_allowed[i]) continue;
        const AtlasOverlapPair *pair = &out->registered.pairs[i];
        int32_t chart0 = pair->component0, chart1 = pair->component1;
        if (chart0 < 0 || chart1 < 0 || (size_t)chart0 >= ncharts ||
            (size_t)chart1 >= ncharts || chart0 == chart1)
            return -1;
        uint64_t rank0 = solution->chart[chart0].rank;
        uint64_t rank1 = solution->chart[chart1].rank;
        if (rank0 >= ncharts || rank1 >= ncharts || rank0 == rank1)
            return -1;
        int32_t chart_lo = rank0 < rank1 ? chart0 : chart1;
        int32_t chart_hi = rank0 < rank1 ? chart1 : chart0;
        size_t face_lo = rank0 < rank1 ? (size_t)pair->face0
                                      : (size_t)pair->face1;
        size_t face_hi = rank0 < rank1 ? (size_t)pair->face1
                                      : (size_t)pair->face0;
        double left, right, v;
        if (face_pair_translation_interval(
                faces, registered_parameter, solution->v,
                face_lo, face_hi, &left, &right, &v) != 0) {
            out->collision_interval_failures++;
            continue;
        }
        int32_t row;
        int row_rc = nearest_common_chart_row(
            present, set->nrows, ncharts, chart_lo, chart_hi,
            (v - fit->v0) / fit->dv, &row);
        if (row_rc < 0) return -1;
        if (row_rc > 0) {
            out->collision_pairs_without_row++;
            continue;
        }
        double current_delta =
            dense_shift[(size_t)row * ncharts + (size_t)chart_hi] -
            dense_shift[(size_t)row * ncharts + (size_t)chart_lo];
        double escape_delta = fmax(0.0, right) + clearance;
        if (!isfinite(current_delta) || !isfinite(escape_delta)) return -1;
        AtlasRibbonCollisionBound *bound = &out->collision_bound[raw_count++];
        memset(bound, 0, sizeof *bound);
        bound->row = row;
        bound->chart_lo = chart_lo;
        bound->chart_hi = chart_hi;
        bound->rank_lo = solution->chart[chart_lo].rank;
        bound->rank_hi = solution->chart[chart_hi].rank;
        bound->face_lo = out->face_source[face_lo];
        bound->face_hi = out->face_source[face_hi];
        bound->winding_lo = solution->chart[chart_lo].winding;
        bound->winding_hi = solution->chart[chart_hi].winding;
        bound->face_pairs = 1;
        bound->v = v;
        bound->min_xyz = out->registered_min_xyz[i];
        bound->current_delta = current_delta;
        bound->escape_delta = escape_delta;
        bound->lower_shift_delta = current_delta + escape_delta;
    }
    out->collision_bound_face_pairs = raw_count;
    qsort(out->collision_bound, raw_count, sizeof(*out->collision_bound),
          compare_collision_bound);
    size_t write = 0;
    for (size_t i = 0; i < raw_count; i++) {
        AtlasRibbonCollisionBound *src = &out->collision_bound[i];
        if (write == 0 || compare_collision_bound(
                &out->collision_bound[write - 1], src) != 0) {
            out->collision_bound[write++] = *src;
            continue;
        }
        AtlasRibbonCollisionBound *dst = &out->collision_bound[write - 1];
        uint32_t count = dst->face_pairs == UINT32_MAX
                       ? UINT32_MAX : dst->face_pairs + 1;
        double min_xyz = fmin(dst->min_xyz, src->min_xyz);
        if (src->lower_shift_delta > dst->lower_shift_delta) *dst = *src;
        dst->face_pairs = count;
        dst->min_xyz = min_xyz;
    }
    out->ncollision_bounds = write;
    return 0;
}

static int build_collision_block_bounds(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    const ArftUvAudit *audit,
    AtlasRibbonCollisionBound **out_bound,
    size_t *out_count)
{
    if (arena == NULL || solution == NULL || set == NULL || fit == NULL ||
        audit == NULL || out_bound == NULL || out_count == NULL ||
        fit->registered_u == NULL || fit->row_chart_shift == NULL ||
        fit->nrows != set->nrows || solution->ncharts == 0 ||
        set->nrows > SIZE_MAX / solution->ncharts ||
        solution->ncharts > SIZE_MAX / solution->ncharts)
        return -1;
    size_t ncharts = solution->ncharts;
    size_t matrix_size = ncharts * ncharts;
    size_t field_size = set->nrows * ncharts;
    uint8_t *pair = (uint8_t *)ARENA_CALLOC(arena, matrix_size, 1);
    size_t unique_pairs = 0;
    for (size_t i = 0; i < audit->ncollision_bounds; i++) {
        const AtlasRibbonCollisionBound *source = &audit->collision_bound[i];
        if (source->chart_lo < 0 || source->chart_hi < 0 ||
            (size_t)source->chart_lo >= ncharts ||
            (size_t)source->chart_hi >= ncharts ||
            source->rank_lo >= source->rank_hi ||
            solution->chart[source->chart_lo].rank != source->rank_lo ||
            solution->chart[source->chart_hi].rank != source->rank_hi)
            return -1;
        size_t key = (size_t)source->chart_lo * ncharts +
                     (size_t)source->chart_hi;
        if (!pair[key]) {
            pair[key] = 1;
            unique_pairs++;
        }
    }
    if (unique_pairs > SIZE_MAX / set->nrows) return -1;
    size_t capacity = unique_pairs * set->nrows;
    AtlasRibbonCollisionBound *bound =
        (AtlasRibbonCollisionBound *)ARENA_ALLOC(
            arena, (capacity ? capacity : 1) * sizeof(*bound));
    double *minimum = (double *)ARENA_ALLOC(
        arena, field_size * sizeof(*minimum));
    double *maximum = (double *)ARENA_ALLOC(
        arena, field_size * sizeof(*maximum));
    for (size_t i = 0; i < field_size; i++) {
        minimum[i] = DBL_MAX;
        maximum[i] = -DBL_MAX;
    }
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *sample = &set->layer_sample[i];
        if (sample->row < 0 || (size_t)sample->row >= set->nrows ||
            sample->chart < 0 || (size_t)sample->chart >= ncharts ||
            !isfinite(fit->registered_u[i]))
            return -1;
        size_t at = (size_t)sample->row * ncharts + (size_t)sample->chart;
        if (fit->registered_u[i] < minimum[at]) minimum[at] = fit->registered_u[i];
        if (fit->registered_u[i] > maximum[at]) maximum[at] = fit->registered_u[i];
    }
    double clearance = set->observation_du > 0.0
                     ? set->observation_du : 1.0e-6;
    size_t count = 0;
    for (size_t chart_lo = 0; chart_lo < ncharts; chart_lo++)
    for (size_t chart_hi = 0; chart_hi < ncharts; chart_hi++) {
        if (!pair[chart_lo * ncharts + chart_hi]) continue;
        for (size_t row = 0; row < set->nrows; row++) {
            size_t lo_at = row * ncharts + chart_lo;
            size_t hi_at = row * ncharts + chart_hi;
            if (minimum[lo_at] == DBL_MAX || maximum[lo_at] == -DBL_MAX ||
                minimum[hi_at] == DBL_MAX || maximum[hi_at] == -DBL_MAX)
                continue;
            double escape = maximum[lo_at] - minimum[hi_at] + clearance;
            if (!(escape > 1.0e-9)) continue;
            AtlasRibbonCollisionBound *dst = &bound[count++];
            memset(dst, 0, sizeof *dst);
            dst->row = (int32_t)row;
            dst->chart_lo = (int32_t)chart_lo;
            dst->chart_hi = (int32_t)chart_hi;
            dst->rank_lo = solution->chart[chart_lo].rank;
            dst->rank_hi = solution->chart[chart_hi].rank;
            dst->face_lo = dst->face_hi = -1;
            dst->winding_lo = solution->chart[chart_lo].winding;
            dst->winding_hi = solution->chart[chart_hi].winding;
            dst->v = fit->v0 + (double)row * fit->dv;
            dst->min_xyz = NAN;
            dst->current_delta =
                fit->row_chart_shift[hi_at] - fit->row_chart_shift[lo_at];
            dst->escape_delta = escape;
            dst->lower_shift_delta = dst->current_delta + escape;
        }
    }
    qsort(bound, count, sizeof(*bound), compare_collision_bound);
    *out_bound = bound;
    *out_count = count;
    return 0;
}
static int build_uv_audit(

    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonFitResult *fit,
    ArftUvAudit *out)
{
    if (arena == NULL || ps == NULL || solution == NULL || set == NULL ||
        opts == NULL || fit == NULL || out == NULL || ps->nf != solution->nfaces ||
        ps->nv != solution->nvertices || ps->uv == NULL || ps->phi == NULL ||
        solution->ncharts == 0 || fit->nrows == 0 || !(fit->dv > 0.0))
        return -1;
    memset(out, 0, sizeof *out);

    size_t nkept = 0;
    for (size_t face = 0; face < ps->nf; face++)
        if (solution->face_keep[face]) nkept++;
    if (nkept == 0 || nkept > (size_t)INT32_MAX) return -1;
    int32_t *faces = (int32_t *)ARENA_ALLOC(
        arena, nkept * 3 * sizeof(*faces));
    out->face_source = (int32_t *)ARENA_ALLOC(
        arena, nkept * sizeof(*out->face_source));
    size_t at = 0;
    for (size_t source = 0; source < ps->nf; source++) {
        if (!solution->face_keep[source]) continue;
        out->face_source[at] = (int32_t)source;
        for (int corner = 0; corner < 3; corner++) {
            int32_t vertex = ps->faces[source * 3 + (size_t)corner];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            faces[at * 3 + (size_t)corner] = vertex;
        }
        at++;
    }
    if (at != nkept) return -1;
    out->nkept_faces = nkept;

    double *dense = NULL;
    if (build_dense_chart_shift(arena, solution, set, fit, &dense) != 0)
        return -1;
    double *registered_parameter = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(*registered_parameter));
    float *source_registered = (float *)ARENA_ALLOC(
        arena, ps->nv * sizeof(*source_registered));
    for (size_t vertex = 0; vertex < ps->nv; vertex++) {
        source_registered[vertex] = ps->uv[vertex * 2];
        int32_t chart = solution->vertex_chart[vertex];
        double shift = 0.0;
        if (chart >= 0 && (size_t)chart < solution->ncharts) {
            double rowf = (solution->v[vertex] - fit->v0) / fit->dv;
            if (rowf <= 0.0) {
                shift = dense[(size_t)chart];
            } else if (rowf >= (double)(fit->nrows - 1)) {
                shift = dense[(fit->nrows - 1) * solution->ncharts +
                              (size_t)chart];
            } else {
                size_t row0 = (size_t)floor(rowf);
                size_t row1 = row0 + 1;
                double t = rowf - (double)row0;
                double a = dense[row0 * solution->ncharts + (size_t)chart];
                double b = dense[row1 * solution->ncharts + (size_t)chart];
                shift = (1.0 - t) * a + t * b;
            }
        }
        registered_parameter[vertex] = solution->u[vertex] + shift;
        if (!isfinite(registered_parameter[vertex])) return -1;
    }

    if (AtlasOverlapAudit_build(
            arena, faces, nkept, ps->nv, solution->u, solution->v,
            source_registered, ps->phi, solution->vertex_chart,
            solution->ncharts, &out->initial) != 0 ||
        !out->initial.broad_phase_complete ||
        out->initial.indexed_faces != nkept)
        return -1;
    if (AtlasOverlapAudit_build(
            arena, faces, nkept, ps->nv,
            registered_parameter, solution->v,
            source_registered, ps->phi, solution->vertex_chart,
            solution->ncharts, &out->registered) != 0 ||
        !out->registered.broad_phase_complete ||
        out->registered.indexed_faces != nkept)
        return -1;

    size_t initial_at = 0, registered_at = 0;
    while (initial_at < out->initial.npairs &&
           registered_at < out->registered.npairs) {
        int comparison = compare_overlap_pair_key(
            &out->initial.pairs[initial_at],
            &out->registered.pairs[registered_at]);
        if (comparison < 0) {
            out->resolved_pairs++;
            initial_at++;
        } else if (comparison > 0) {
            out->new_pairs++;
            registered_at++;
        } else {
            out->persisted_pairs++;
            initial_at++;
            registered_at++;
        }
    }
    out->resolved_pairs += out->initial.npairs - initial_at;
    out->new_pairs += out->registered.npairs - registered_at;

    out->residual = (ArftResidualRef *)ARENA_ALLOC(
        arena, (solution->nresiduals ? solution->nresiduals : 1) *
               sizeof(*out->residual));
    for (size_t i = 0; i < solution->nresiduals; i++) {
        int32_t face0 = solution->residual[i].face0;
        int32_t face1 = solution->residual[i].face1;
        if (face0 > face1) {
            int32_t temporary = face0; face0 = face1; face1 = temporary;
        }
        out->residual[i].face0 = face0;
        out->residual[i].face1 = face1;
        out->residual[i].allowed = solution->residual[i].allowed;
        out->residual[i].reason = solution->residual[i].reason;
    }
    qsort(out->residual, solution->nresiduals,
          sizeof(*out->residual), compare_residual_ref);
    for (size_t i = 0; i < solution->nresiduals; i++) {
        if (out->nresiduals > 0 &&
            out->residual[out->nresiduals - 1].face0 ==
                out->residual[i].face0 &&
            out->residual[out->nresiduals - 1].face1 ==
                out->residual[i].face1) {
            if (out->residual[i].allowed)
                out->residual[out->nresiduals - 1].allowed = 1;
            continue;
        }
        out->residual[out->nresiduals++] = out->residual[i];
    }

    if (solution->ncharts > SIZE_MAX / solution->ncharts) return -1;
    size_t chart_matrix = solution->ncharts * solution->ncharts;
    uint8_t *direct_pair = (uint8_t *)ARENA_CALLOC(arena, chart_matrix, 1);
    for (size_t i = 0; i < solution->nresiduals; i++) {
        const AtlasSolutionResidual *residual = &solution->residual[i];
        if (!residual->allowed || residual->chart0 < 0 ||
            residual->chart1 < 0 ||
            (size_t)residual->chart0 >= solution->ncharts ||
            (size_t)residual->chart1 >= solution->ncharts)
            continue;
        direct_pair[(size_t)residual->chart0 * solution->ncharts +
                    (size_t)residual->chart1] = 1;
        direct_pair[(size_t)residual->chart1 * solution->ncharts +
                    (size_t)residual->chart0] = 1;
    }
    size_t previous = SIZE_MAX;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *current = &set->layer_sample[i];
        if (previous != SIZE_MAX &&
            set->layer_sample[previous].row != current->row)
            previous = SIZE_MAX;
        if (current->cluster_count != 1) continue;
        if (previous != SIZE_MAX) {
            const AtlasRibbonLayerSample *prior = &set->layer_sample[previous];
            double tangent_dot = prior->tangent[0] * current->tangent[0] +
                                 prior->tangent[1] * current->tangent[1];
            if (prior->chart != current->chart && prior->chart >= 0 &&
                current->chart >= 0 &&
                layer_distance(prior, current) <= opts->local_xyz_tolerance &&
                tangent_dot >= opts->tangent_dot_min) {
                direct_pair[(size_t)prior->chart * solution->ncharts +
                            (size_t)current->chart] = 1;
                direct_pair[(size_t)current->chart * solution->ncharts +
                            (size_t)prior->chart] = 1;
            }
        }
        previous = i;
    }

    size_t npairs = out->registered.npairs;
    out->registered_min_xyz = (double *)ARENA_ALLOC(
        arena, (npairs ? npairs : 1) * sizeof(*out->registered_min_xyz));
    out->registered_allowed = (uint8_t *)ARENA_CALLOC(
        arena, (npairs ? npairs : 1), 1);
    out->registered_reason = (int32_t *)ARENA_CALLOC(
        arena, (npairs ? npairs : 1), sizeof(*out->registered_reason));
    for (size_t i = 0; i < out->registered.npairs; i++) {
        const AtlasOverlapPair *pair = &out->registered.pairs[i];
        int32_t face0 = out->face_source[pair->face0];
        int32_t face1 = out->face_source[pair->face1];
        const ArftResidualRef *known =
            find_residual_ref(out, face0, face1);
        if (known != NULL) out->registered_known_pairs++;
        double min_xyz = face_min_xyz(
            ps, faces, (size_t)pair->face0, (size_t)pair->face1);
        out->registered_min_xyz[i] = min_xyz;
        int reason = 0;
        int32_t chart0 = pair->component0, chart1 = pair->component1;
        if (chart0 == chart1) {
            reason = 3;
            out->registered_self_pairs++;
        } else if (known != NULL && known->allowed) {
            reason = known->reason > 0 ? known->reason : 1;
            out->registered_direct_pairs++;
        } else if (direct_pair[(size_t)chart0 * solution->ncharts +
                               (size_t)chart1] &&
                   (opts->local_xyz_tolerance <= 0.0 ||
                    min_xyz <= opts->local_xyz_tolerance)) {
            reason = 1;
            out->registered_direct_pairs++;
        } else if (min_xyz <= opts->local_xyz_tolerance &&
                   solution->chart[chart0].winding ==
                       solution->chart[chart1].winding) {
            reason = 2;
            out->registered_same_winding_local_pairs++;
        }
        out->registered_reason[i] = reason;
        out->registered_allowed[i] = reason != 0;
        if (reason != 0) {
            out->registered_allowed_pairs++;
        } else {
            out->registered_hard_pairs++;
            if (known == NULL) out->registered_unclassified_pairs++;
        }
    }
    if (build_collision_bounds(
            arena, faces, registered_parameter, dense, solution, set, opts,
            fit, out) != 0)
        return -1;
    return 0;
}

static void world_point(const AtlasRibbonObservationSet *set,
                        double v, const double p[2], double xyz[3])
{
    for (int d = 0; d < 3; d++)
        xyz[d] = set->axis_point[d] + v * set->axis[d] +
                 p[0] * set->basis0[d] + p[1] * set->basis1[d];
}

static int same_key(const AtlasRibbonObservation *ob,
                    const AtlasRibbonTarget *target)
{
    return ob->row == target->row && ob->column == target->column;
}

static void target_color(const AtlasRibbonTarget *target, double rgb[3])
{
    if (target->accepted && target->charts > 1) {
        rgb[0] = 0.10; rgb[1] = 0.90; rgb[2] = 0.90;
    } else if (target->accepted) {
        rgb[0] = 0.15; rgb[1] = 0.85; rgb[2] = 0.25;
    } else if (target->known_reason == 4) {
        rgb[0] = 1.00; rgb[1] = 0.55; rgb[2] = 0.05;
    } else {
        rgb[0] = 0.95; rgb[1] = 0.10; rgb[2] = 0.20;
    }
}

static int write_observations_obj(
    const char *dir,
    const AtlasRibbonObservationSet *set)
{
    FILE *fp = open_out(dir, "observations_world.obj");
    if (fp == NULL) return -1;
    fprintf(fp,
            "# Exact constant-v triangle sections at common fixed-U samples.\n"
            "# green=single observation, cyan=local duplicate, "
            "orange=known topology residual, red=unclassified remote conflict\n");
    size_t observation_at = 0;
    uint64_t vertex = 1;
    for (size_t i = 0; i < set->ntarget; i++) {
        const AtlasRibbonTarget *target = &set->target[i];
        while (observation_at < set->nobservation &&
               !same_key(&set->observation[observation_at], target))
            observation_at++;
        size_t first = observation_at;
        while (observation_at < set->nobservation &&
               same_key(&set->observation[observation_at], target))
            observation_at++;
        double rgb[3];
        target_color(target, rgb);
        for (size_t j = first; j < observation_at; j++) {
            const AtlasRibbonObservation *ob = &set->observation[j];
            double a[3], b[3];
            world_point(set, ob->v, ob->p, a);
            double q[2] = {
                ob->p[0] + 2.0 * ob->tangent[0],
                ob->p[1] + 2.0 * ob->tangent[1]
            };
            world_point(set, ob->v, q, b);
            fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f %.3f\n",
                    a[0], a[1], a[2], rgb[0], rgb[1], rgb[2]);
            fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f %.3f\n",
                    b[0], b[1], b[2], rgb[0], rgb[1], rgb[2]);
            fprintf(fp, "l %" PRIu64 " %" PRIu64 "\n",
                    vertex, vertex + 1);
            vertex += 2;
        }
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_targets_csv(
    const char *dir,
    const AtlasRibbonObservationSet *set)
{
    FILE *fp = open_out(dir, "observation_targets.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,column,u,v,accepted,observations,charts,clusters,"
                "chart0,chart1,face0,face1,known_reason,"
                "max_cluster_separation,p0,p1,tangent0,tangent1\n");
    for (size_t i = 0; i < set->ntarget; i++) {
        const AtlasRibbonTarget *t = &set->target[i];
        fprintf(fp, "%d,%d,%.17g,%.17g,%d,%u,%u,%u,%d,%d,%d,%d,%d,"
                    "%.17g,%.17g,%.17g,%.17g,%.17g\n",
                t->row, t->column, t->u, t->v, t->accepted,
                t->observations, t->charts, t->clusters,
                t->chart0, t->chart1, t->face0, t->face1,
                t->known_reason, t->max_cluster_separation,
                t->p[0], t->p[1], t->tangent[0], t->tangent[1]);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_conflicts_csv(
    const char *dir,
    const AtlasRibbonObservationSet *set)
{
    FILE *fp = open_out(dir, "remote_collisions.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,column,u,v,observations,charts,clusters,chart0,chart1,"
                "face0,face1,known_reason,max_cluster_separation\n");
    for (size_t i = 0; i < set->ntarget; i++) {
        const AtlasRibbonTarget *t = &set->target[i];
        if (t->accepted) continue;
        fprintf(fp, "%d,%d,%.17g,%.17g,%u,%u,%u,%d,%d,%d,%d,%d,%.17g\n",
                t->row, t->column, t->u, t->v, t->observations,
                t->charts, t->clusters, t->chart0, t->chart1,
                t->face0, t->face1, t->known_reason,
                t->max_cluster_separation);
    }
    return fclose(fp) == 0 ? 0 : -1;
}
static int write_bridge_cuts_csv(
    const char *dir,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts)
{
    FILE *fp = open_out(dir, "bridge_cuts.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,v,u0,u1,du,distance,metric_limit,reason,"
                "chart0,chart1,face0,face1\n");
    const AtlasRibbonTarget *previous = NULL;
    for (size_t i = 0; i < set->ntarget; i++) {
        const AtlasRibbonTarget *current = &set->target[i];
        if (!current->accepted) continue;
        if (previous != NULL && previous->row == current->row) {
            double du = current->u - previous->u;
            double dx = current->p[0] - previous->p[0];
            double dy = current->p[1] - previous->p[1];
            double distance = sqrt(dx * dx + dy * dy);
            double limit = opts->max_bridge_stretch * du +
                           opts->bridge_slack;
            int cut_u = opts->max_fill_u > 0.0 &&
                        du > opts->max_fill_u;
            int cut_metric = opts->max_bridge_stretch > 0.0 &&
                             distance > limit;
            if (cut_u || cut_metric) {
                fprintf(fp,
                    "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,"
                    "%d,%d,%d,%d\n",
                    current->row, current->v, previous->u, current->u,
                    du, distance, limit, cut_u ? "u_gap" : "metric_jump",
                    previous->chart0, current->chart0,
                    previous->face0, current->face0);
            }
        }
        previous = current;
    }
    return fclose(fp) == 0 ? 0 : -1;
}


static const char *register_kind_name(int kind)
{
    switch ((AtlasRibbonRegisterKind)kind) {
    case ATLAS_RIBBON_REGISTER_REMOTE: return "remote_layer";
    case ATLAS_RIBBON_REGISTER_METRIC: return "metric_jump";
    case ATLAS_RIBBON_REGISTER_ORDER:  return "rank_inversion";
    default:                           return "unknown";
    }
}

static int write_registered_samples_csv(
    const char *dir,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit)
{
    FILE *fp = open_out(dir, "registered_layer_samples.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,source_column,source_u,registered_u,shift_u,v,"
                "chart,rank,group,winding,charts,cluster_index,cluster_count,"
                "face,p0,p1,tangent0,tangent1\n");
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        double registered = fit->registered_u[i];
        fprintf(fp,
                "%d,%d,%.17g,%.17g,%.17g,%.17g,%d,%" PRIu64
                ",%d,%d,%u,%u,%u,%d,%.17g,%.17g,%.17g,%.17g\n",
                s->row, s->source_column, s->source_u, registered,
                registered - s->source_u, s->v, s->chart, s->rank,
                s->group, s->winding, s->charts, s->cluster_index,
                s->cluster_count, s->face, s->p[0], s->p[1],
                s->tangent[0], s->tangent[1]);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_register_bounds_csv(
    const char *dir,
    const AtlasRibbonFitResult *fit)
{
    FILE *fp = open_out(dir, "u_registration_bounds.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "kind,row,chart_lo,rank_lo,chart_hi,rank_hi,sample_lo,"
                "sample_hi,source_u_delta,endpoint_distance,required_arc,"
                "lower_shift_delta,final_shift_delta,slack,updates\n");
    for (size_t i = 0; i < fit->nregister_bounds; i++) {
        const AtlasRibbonRegisterBound *b = &fit->register_bound[i];
        fprintf(fp,
                "%s,%d,%d,%" PRIu64 ",%d,%" PRIu64 ",%" PRIu64
                ",%" PRIu64 ",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%u\n",
                register_kind_name(b->kind), b->row,
                b->chart_lo, b->rank_lo, b->chart_hi, b->rank_hi,
                b->sample_lo, b->sample_hi, b->source_u_delta,
                b->endpoint_distance, b->required_arc,
                b->lower_shift_delta, b->final_shift_delta,
                b->final_shift_delta - b->lower_shift_delta, b->updates);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int overlap_pair_present(const AtlasOverlapAudit *audit,
                                const AtlasOverlapPair *key)
{
    size_t lo = 0, hi = audit->npairs;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int comparison = compare_overlap_pair_key(&audit->pairs[mid], key);
        if (comparison < 0) lo = mid + 1;
        else hi = mid;
    }
    return lo < audit->npairs &&
           compare_overlap_pair_key(&audit->pairs[lo], key) == 0;
}

static int write_uv_audit(const char *dir, const ArftUvAudit *audit)
{
    FILE *fp = open_out(dir, "u_overlap_pairs.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "initial,known,ledger_allowed,allowed,reason,min_xyz,chart0,chart1,face0,face1,"
                "compact_face0,compact_face1,turn,phase_residual,"
                "parameter_du,registered_du,registered_parameter_shift\n");
    for (size_t i = 0; i < audit->registered.npairs; i++) {
        const AtlasOverlapPair *pair = &audit->registered.pairs[i];
        int32_t source0 = audit->face_source[pair->face0];
        int32_t source1 = audit->face_source[pair->face1];
        const ArftResidualRef *known =
            find_residual_ref(audit, source0, source1);
        fprintf(fp,
                "%d,%d,%d,%d,%d,%.17g,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g\n",
                overlap_pair_present(&audit->initial, pair), known != NULL,
                known != NULL ? known->allowed : 0,
                audit->registered_allowed[i], audit->registered_reason[i],
                audit->registered_min_xyz[i],
                pair->component0, pair->component1, source0, source1,
                pair->face0, pair->face1, pair->turn,
                pair->phase_residual, pair->parameter_du,
                pair->registered_du, pair->registered_parameter_shift);
    }
    if (fclose(fp) != 0) return -1;

    fp = open_out(dir, "u_collision_bounds.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,v,chart_lo,rank_lo,winding_lo,chart_hi,rank_hi,"
                "winding_hi,face_lo,face_hi,face_pairs,min_xyz,"
                "current_shift_delta,escape_delta,lower_shift_delta\n");
    for (size_t i = 0; i < audit->ncollision_bounds; i++) {
        const AtlasRibbonCollisionBound *bound = &audit->collision_bound[i];
        fprintf(fp,
                "%d,%.17g,%d,%" PRIu64 ",%d,%d,%" PRIu64
                ",%d,%d,%d,%u,%.17g,%.17g,%.17g,%.17g\n",
                bound->row, bound->v,
                bound->chart_lo, bound->rank_lo, bound->winding_lo,
                bound->chart_hi, bound->rank_hi, bound->winding_hi,
                bound->face_lo, bound->face_hi, bound->face_pairs,
                bound->min_xyz, bound->current_delta, bound->escape_delta,
                bound->lower_shift_delta);
    }
    if (fclose(fp) != 0) return -1;

    fp = open_out(dir, "u_overlap_audit.json");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"kept_faces\": %zu,\n"
        "  \"checkpoint_exact_pairs\": %zu,\n"
        "  \"checkpoint_cross_chart_pairs\": %zu,\n"
        "  \"checkpoint_same_chart_pairs\": %zu,\n"
        "  \"registered_exact_pairs\": %zu,\n"
        "  \"registered_cross_chart_pairs\": %zu,\n"
        "  \"registered_same_chart_pairs\": %zu,\n"
        "  \"persisted_pairs\": %zu,\n"
        "  \"new_pairs\": %zu,\n"
        "  \"resolved_pairs\": %zu,\n"
        "  \"registered_known_pairs\": %zu,\n"
        "  \"registered_allowed_pairs\": %zu,\n"
        "  \"registered_direct_pairs\": %zu,\n"
        "  \"registered_same_winding_local_pairs\": %zu,\n"
        "  \"registered_self_pairs\": %zu,\n"
        "  \"registered_hard_pairs\": %zu,\n"
        "  \"registered_unclassified_pairs\": %zu,\n"
        "  \"collision_bound_rows\": %zu,\n"
        "  \"collision_bound_face_pairs\": %zu,\n"
        "  \"collision_pairs_without_common_row\": %zu,\n"
        "  \"collision_interval_failures\": %zu,\n"
        "  \"checkpoint_broad_phase_candidates\": %zu,\n"
        "  \"registered_broad_phase_candidates\": %zu,\n"
        "  \"checkpoint_broad_phase_cell_size\": %.17g,\n"
        "  \"registered_broad_phase_cell_size\": %.17g,\n"
        "  \"broad_phase_complete\": true\n"
        "}\n",
        audit->nkept_faces,
        audit->initial.npairs, audit->initial.cross_component_pairs,
        audit->initial.same_component_pairs,
        audit->registered.npairs, audit->registered.cross_component_pairs,
        audit->registered.same_component_pairs,
        audit->persisted_pairs, audit->new_pairs, audit->resolved_pairs,
        audit->registered_known_pairs, audit->registered_allowed_pairs,
        audit->registered_direct_pairs,
        audit->registered_same_winding_local_pairs,
        audit->registered_self_pairs,
        audit->registered_hard_pairs,
        audit->registered_unclassified_pairs,
        audit->ncollision_bounds,
        audit->collision_bound_face_pairs,
        audit->collision_pairs_without_row,
        audit->collision_interval_failures,
        audit->initial.broad_phase_candidate_pairs,
        audit->registered.broad_phase_candidate_pairs,
        audit->initial.broad_phase_cell_size,
        audit->registered.broad_phase_cell_size);
    return fclose(fp) == 0 ? 0 : -1;
}
static int write_collision_ledger(
    const char *dir,
    const AtlasRibbonCollisionBound *ledger,
    size_t nledger)
{
    FILE *fp = open_out(dir, "u_collision_ledger.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,v,chart_lo,rank_lo,winding_lo,chart_hi,rank_hi,"
                "winding_hi,face_lo,face_hi,face_pairs,min_xyz,"
                "proposal_current_delta,proposal_escape_delta,"
                "applied_lower_shift_delta\n");
    for (size_t i = 0; i < nledger; i++) {
        const AtlasRibbonCollisionBound *bound = &ledger[i];
        fprintf(fp,
                "%d,%.17g,%d,%" PRIu64 ",%d,%d,%" PRIu64
                ",%d,%d,%d,%u,%.17g,%.17g,%.17g,%.17g\n",
                bound->row, bound->v,
                bound->chart_lo, bound->rank_lo, bound->winding_lo,
                bound->chart_hi, bound->rank_hi, bound->winding_hi,
                bound->face_lo, bound->face_hi, bound->face_pairs,
                bound->min_xyz, bound->current_delta, bound->escape_delta,
                bound->lower_shift_delta);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_rows_csv(
    const char *dir,
    const AtlasRibbonFitResult *fit)
{
    FILE *fp = open_out(dir, "row_metrics.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row,v,accepted_targets,fitted_nodes,supported_nodes,"
                "curve_segments,crossings,crossing_cuts,u_gap_cuts,metric_jump_cuts,"
                "direct_spans,spiral_spans,topology_cuts,"
                "u_min,u_max,speed_mean,"
                "speed_rms_error,speed_max_error,phase_backstep_fraction,"
                "phase_max_backstep,winding_span,radius_pitch_median,"
                "radius_order_tests,radius_order_violations\n");
    for (size_t row = 0; row < fit->nrows; row++) {
        const AtlasRibbonRowMetric *m = &fit->row_metric[row];
        fprintf(fp, "%d,%.17g,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu,%zu\n",
                m->row, fit->v0 + (double)row * fit->dv,
                m->accepted_targets, m->fitted_nodes, m->supported_nodes,
                m->curve_segments, m->crossings, m->crossing_cuts,
                m->u_gap_cuts, m->metric_jump_cuts,
                m->direct_spans, m->spiral_spans, m->topology_cuts, m->u_min, m->u_max,
                m->speed_mean, m->speed_rms_error, m->speed_max_error,
                m->phase_backstep_fraction, m->phase_max_backstep,
                m->winding_span, m->radius_pitch_median,
                m->radius_order_tests, m->radius_order_violations);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_ribbon_obj(
    Arena_T arena,
    const char *dir,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    size_t *out_faces)
{
    FILE *fp = open_out(dir, "ribbon.obj");
    if (fp == NULL) return -1;
    size_t nnode = fit->nrows * fit->ncolumns;
    int64_t *index = (int64_t *)ARENA_CALLOC(
        arena, nnode, sizeof(int64_t));
    double umin = DBL_MAX, vmin = DBL_MAX;
    int64_t next = 1;
    fprintf(fp, "# Fitted complete ribbon; source-space vertices are z,y,x.\n");
    fprintf(fp, "# provenance colors: observed=green chart_interp=teal "
                "seam_interp=orange spiral_fill=blue\n");
    for (size_t row = 0; row < fit->nrows; row++) {
        double v = fit->v0 + (double)row * fit->dv;
        for (size_t column = 0; column < fit->ncolumns; column++) {
            size_t node = row * fit->ncolumns + column;
            if (!fit->valid[node]) continue;
            double u = fit->u0 + (double)column * fit->du;
            double xyz[3];
            world_point(set, v, &fit->p[node * 2], xyz);
            double rgb[3];
            uint8_t provenance = fit->provenance != NULL
                               ? fit->provenance[node]
                               : (fit->support[node]
                                  ? ATLAS_RIBBON_NODE_CHART_INTERP
                                  : ATLAS_RIBBON_NODE_SPIRAL_FILL);
            if (provenance == ATLAS_RIBBON_NODE_OBSERVED) {
                rgb[0] = 0.20; rgb[1] = 0.85; rgb[2] = 0.35;
            } else if (provenance == ATLAS_RIBBON_NODE_CHART_INTERP) {
                rgb[0] = 0.15; rgb[1] = 0.72; rgb[2] = 0.62;
            } else if (provenance == ATLAS_RIBBON_NODE_SEAM_INTERP) {
                rgb[0] = 0.95; rgb[1] = 0.55; rgb[2] = 0.20;
            } else if (provenance == ATLAS_RIBBON_NODE_SPIRAL_FILL) {
                rgb[0] = 0.45; rgb[1] = 0.65; rgb[2] = 0.95;
            } else {
                rgb[0] = 0.80; rgb[1] = 0.15; rgb[2] = 0.65;
            }
            fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f %.3f\n",
                    xyz[0], xyz[1], xyz[2], rgb[0], rgb[1], rgb[2]);
            index[node] = next++;
            if (u < umin) umin = u;
            if (v < vmin) vmin = v;
        }
    }
    if (umin == DBL_MAX) {
        fclose(fp);
        return -1;
    }
    for (size_t row = 0; row < fit->nrows; row++) {
        double v = fit->v0 + (double)row * fit->dv;
        for (size_t column = 0; column < fit->ncolumns; column++) {
            size_t node = row * fit->ncolumns + column;
            if (!fit->valid[node]) continue;
            double u = fit->u0 + (double)column * fit->du;
            fprintf(fp, "vt %.17g %.17g\n", u - umin, v - vmin);
        }
    }
    size_t faces = 0;
    for (size_t row = 0; row + 1 < fit->nrows; row++) {
        for (size_t column = 0; column + 1 < fit->ncolumns; column++) {
            size_t a = row * fit->ncolumns + column;
            size_t b = a + 1;
            size_t c = a + fit->ncolumns;
            size_t d = c + 1;
            if (index[a] == 0 || index[b] == 0 ||
                index[c] == 0 || index[d] == 0 ||
                (fit->u_edge != NULL &&
                 (!fit->u_edge[a] || !fit->u_edge[c])) ||
                (fit->v_edge != NULL &&
                 (!fit->v_edge[a] || !fit->v_edge[b])))
                continue;
            fprintf(fp, "f %" PRId64 "/%" PRId64 " %" PRId64 "/%" PRId64
                        " %" PRId64 "/%" PRId64 "\n",
                    index[a], index[a], index[b], index[b],
                    index[d], index[d]);
            fprintf(fp, "f %" PRId64 "/%" PRId64 " %" PRId64 "/%" PRId64
                        " %" PRId64 "/%" PRId64 "\n",
                    index[a], index[a], index[d], index[d],
                    index[c], index[c]);
            faces += 2;
        }
    }
    *out_faces = faces;
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_ribbon_rows_obj(
    const char *dir,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit)
{
    FILE *fp = open_out(dir, "ribbon_rows_world.obj");
    if (fp == NULL) return -1;
    uint64_t next = 1;
    for (size_t row = 0; row < fit->nrows; row++) {
        double v = fit->v0 + (double)row * fit->dv;
        uint64_t run_first = 0, run_last = 0;
        for (size_t column = 0; column <= fit->ncolumns; column++) {
            size_t node = row * fit->ncolumns + column;
            int valid = column < fit->ncolumns && fit->valid[node];
            if (valid && run_first != 0 && column > 0 &&
                fit->u_edge != NULL && !fit->u_edge[node - 1]) {
                if (run_last > run_first) {
                    fprintf(fp, "l");
                    for (uint64_t i = run_first; i <= run_last; i++)
                        fprintf(fp, " %" PRIu64, i);
                    fprintf(fp, "\n");
                }
                run_first = run_last = 0;
            }
            if (!valid) {
                if (run_last > run_first) {
                    fprintf(fp, "l");
                    for (uint64_t i = run_first; i <= run_last; i++)
                        fprintf(fp, " %" PRIu64, i);
                    fprintf(fp, "\n");
                }
                run_first = run_last = 0;
                continue;
            }
            double xyz[3];
            world_point(set, v, &fit->p[node * 2], xyz);
            fprintf(fp, "v %.9g %.9g %.9g\n", xyz[0], xyz[1], xyz[2]);
            if (run_first == 0) run_first = next;
            run_last = next++;
        }
    }
    return fclose(fp) == 0 ? 0 : -1;
}

typedef enum {
    ARFT_COVERAGE_BACKGROUND = 0,
    ARFT_COVERAGE_RIBBON = 1,
    ARFT_COVERAGE_PREFIT_CULL = 2,
    ARFT_COVERAGE_INVALID_CHART = 3,
    ARFT_COVERAGE_OUTSIDE_V = 4,
    ARFT_COVERAGE_OUTSIDE_U = 5,
    ARFT_COVERAGE_NO_U_SUPPORT = 6,
    ARFT_COVERAGE_SUBGRID_U = 7,
    ARFT_COVERAGE_U_GAP = 8,
    ARFT_COVERAGE_METRIC_U = 9,
    ARFT_COVERAGE_TOPOLOGY_U = 10,
    ARFT_COVERAGE_CROSSING_U = 11,
    ARFT_COVERAGE_VERTICAL_METRIC = 12,
    ARFT_COVERAGE_GEOMETRY_MISMATCH = 13,
    ARFT_COVERAGE_REASON_COUNT = 14
} ArftCoverageReason;

static const char *coverage_reason_name(int reason)
{
    static const char *name[ARFT_COVERAGE_REASON_COUNT] = {
        "background",
        "ribboned",
        "prefit_cull",
        "invalid_chart",
        "outside_v_lattice",
        "outside_u_lattice",
        "no_u_support",
        "subgrid_u_island",
        "u_gap_cut",
        "metric_u_cut",
        "topology_u_cut",
        "crossing_u_cut",
        "vertical_metric_cut",
        "geometry_mismatch"
    };
    return reason >= 0 && reason < ARFT_COVERAGE_REASON_COUNT
         ? name[reason] : "invalid";
}

/* Priority used only when several source layers paint the same atlas pixel.
 * Every layer is still counted separately in the JSON/CSV audit. */
static int coverage_reason_priority(int reason)
{
    static const uint8_t priority[ARFT_COVERAGE_REASON_COUNT] = {
        0, 1, 4, 5, 2, 2, 3, 6, 7, 9, 10, 12, 11, 13
    };
    return reason >= 0 && reason < ARFT_COVERAGE_REASON_COUNT
         ? priority[reason] : 0;
}

static int coverage_reason_from_u_reject(uint8_t reject)
{
    switch ((AtlasRibbonEdgeReject)reject) {
    case ATLAS_RIBBON_EDGE_SUBGRID: return ARFT_COVERAGE_SUBGRID_U;
    case ATLAS_RIBBON_EDGE_U_GAP: return ARFT_COVERAGE_U_GAP;
    case ATLAS_RIBBON_EDGE_METRIC: return ARFT_COVERAGE_METRIC_U;
    case ATLAS_RIBBON_EDGE_TOPOLOGY: return ARFT_COVERAGE_TOPOLOGY_U;
    case ATLAS_RIBBON_EDGE_CROSSING: return ARFT_COVERAGE_CROSSING_U;
    default: return ARFT_COVERAGE_NO_U_SUPPORT;
    }
}

static double coverage_dense_shift(const double *dense, size_t nrows,
                                   size_t ncharts, int32_t chart,
                                   double rowf)
{
    if (rowf <= 0.0) return dense[(size_t)chart];
    if (rowf >= (double)(nrows - 1))
        return dense[(nrows - 1) * ncharts + (size_t)chart];
    size_t row0 = (size_t)floor(rowf);
    size_t row1 = row0 + 1;
    double t = rowf - (double)row0;
    return (1.0 - t) * dense[row0 * ncharts + (size_t)chart] +
           t * dense[row1 * ncharts + (size_t)chart];
}

static int coverage_classify(const AtlasRibbonObservationSet *set,
                             const AtlasRibbonFitOptions *opts,
                             const AtlasRibbonFitResult *fit,
                             double registered_u, double v,
                             const double source_p[2],
                             double *out_residual)
{
    *out_residual = NAN;
    double row_coordinate = (v - fit->v0) / fit->dv;
    int64_t row0 = (int64_t)floor(row_coordinate);
    if (row0 < 0 || row0 + 1 >= (int64_t)fit->nrows)
        return ARFT_COVERAGE_OUTSIDE_V;
    double column_coordinate = (registered_u - fit->u0) / fit->du;
    int64_t column = (int64_t)floor(column_coordinate);
    if (column < 0 || column + 1 >= (int64_t)fit->ncolumns)
        return ARFT_COVERAGE_OUTSIDE_U;

    size_t a = (size_t)row0 * fit->ncolumns + (size_t)column;
    size_t c = a + fit->ncolumns;
    int lower_u = fit->u_edge != NULL && fit->u_edge[a];
    int upper_u = fit->u_edge != NULL && fit->u_edge[c];
    if (!lower_u || !upper_u) {
        int lower_reason = fit->u_reject != NULL
                         ? coverage_reason_from_u_reject(fit->u_reject[a])
                         : ARFT_COVERAGE_NO_U_SUPPORT;
        int upper_reason = fit->u_reject != NULL
                         ? coverage_reason_from_u_reject(fit->u_reject[c])
                         : ARFT_COVERAGE_NO_U_SUPPORT;
        if (lower_u) return upper_reason;
        if (upper_u) return lower_reason;
        return coverage_reason_priority(lower_reason) >=
               coverage_reason_priority(upper_reason)
             ? lower_reason : upper_reason;
    }
    if (fit->v_edge == NULL || !fit->v_edge[a] || !fit->v_edge[a + 1])
        return ARFT_COVERAGE_VERTICAL_METRIC;
    if (!fit->valid[a] || !fit->valid[a + 1] ||
        !fit->valid[c] || !fit->valid[c + 1])
        return ARFT_COVERAGE_NO_U_SUPPORT;

    double tx = column_coordinate - (double)column;
    double ty = row_coordinate - (double)row0;
    double lower[2] = {
        (1.0 - tx) * fit->p[a * 2] + tx * fit->p[(a + 1) * 2],
        (1.0 - tx) * fit->p[a * 2 + 1] + tx * fit->p[(a + 1) * 2 + 1]
    };
    double upper[2] = {
        (1.0 - tx) * fit->p[c * 2] + tx * fit->p[(c + 1) * 2],
        (1.0 - tx) * fit->p[c * 2 + 1] + tx * fit->p[(c + 1) * 2 + 1]
    };
    double predicted[2] = {
        (1.0 - ty) * lower[0] + ty * upper[0],
        (1.0 - ty) * lower[1] + ty * upper[1]
    };
    *out_residual = hypot(source_p[0] - predicted[0],
                          source_p[1] - predicted[1]);
    double tolerance = opts->local_xyz_tolerance > 0.0
                     ? opts->local_xyz_tolerance : set->observation_du;
    return *out_residual <= tolerance
         ? ARFT_COVERAGE_RIBBON : ARFT_COVERAGE_GEOMETRY_MISMATCH;
}

static double coverage_orient(double ax, double ay, double bx, double by,
                              double px, double py)
{
    return (bx - ax) * (py - ay) - (by - ay) * (px - ax);
}

static int write_coverage_audit(
    Arena_T arena,
    const char *dir,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonFitResult *fit)
{
    if (arena == NULL || dir == NULL || ps == NULL || solution == NULL ||
        set == NULL || opts == NULL || fit == NULL || fit->u_edge == NULL ||
        fit->v_edge == NULL || fit->row_chart_shift == NULL ||
        solution->ncharts == 0 || ps->nf != solution->nfaces ||
        ps->nv != solution->nvertices)
        return -1;

    double umin = DBL_MAX, umax = -DBL_MAX;
    double vmin = DBL_MAX, vmax = -DBL_MAX;
    for (size_t face = 0; face < ps->nf; face++)
    for (int corner = 0; corner < 3; corner++) {
        int32_t vertex = ps->faces[face * 3 + (size_t)corner];
        if (vertex < 0 || (size_t)vertex >= ps->nv) continue;
        double u = solution->u[vertex], v = solution->v[vertex];
        if (!isfinite(u) || !isfinite(v)) continue;
        if (u < umin) umin = u;
        if (u > umax) umax = u;
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    if (!(umin < umax) || !(vmin < vmax)) return -1;
    double origin_u = floor(umin), origin_v = floor(vmin);
    double end_u = ceil(umax), end_v = ceil(vmax);
    if (!(end_u > origin_u) || !(end_v > origin_v) ||
        end_u - origin_u > (double)INT32_MAX ||
        end_v - origin_v > (double)INT32_MAX)
        return -1;
    int width = (int)(end_u - origin_u);
    int height = (int)(end_v - origin_v);
    if (width <= 0 || height <= 0 ||
        (size_t)width > SIZE_MAX / (size_t)height)
        return -1;
    size_t pixels = (size_t)width * (size_t)height;
    uint8_t *reason = (uint8_t *)ARENA_CALLOC(arena, pixels, 1);
    uint8_t *layers = (uint8_t *)ARENA_CALLOC(arena, pixels, 1);
    double *dense = NULL;
    if (build_dense_chart_shift(arena, solution, set, fit, &dense) != 0)
        return -1;
    if (solution->ncharts > SIZE_MAX / ARFT_COVERAGE_REASON_COUNT)
        return -1;
    uint64_t *chart_reason = (uint64_t *)ARENA_CALLOC(
        arena, solution->ncharts * ARFT_COVERAGE_REASON_COUNT,
        sizeof(*chart_reason));
    uint64_t layer_reason[ARFT_COVERAGE_REASON_COUNT] = {0};
    uint64_t display_reason[ARFT_COVERAGE_REASON_COUNT] = {0};
    uint64_t raster_faces = 0, degenerate_faces = 0;
    uint64_t geometry_tests = 0;
    double geometry_sum2 = 0.0, geometry_max = 0.0;

    for (size_t face = 0; face < ps->nf; face++) {
        int32_t vi[3] = {
            ps->faces[face * 3], ps->faces[face * 3 + 1],
            ps->faces[face * 3 + 2]
        };
        int valid_vertex = 1;
        for (int k = 0; k < 3; k++)
            if (vi[k] < 0 || (size_t)vi[k] >= ps->nv ||
                !isfinite(solution->u[vi[k]]) ||
                !isfinite(solution->v[vi[k]]))
                valid_vertex = 0;
        if (!valid_vertex) continue;
        double u[3], v[3];
        for (int k = 0; k < 3; k++) {
            u[k] = solution->u[vi[k]];
            v[k] = solution->v[vi[k]];
        }
        double area = coverage_orient(u[0], v[0], u[1], v[1], u[2], v[2]);
        if (!(fabs(area) > 1.0e-12)) {
            degenerate_faces++;
            continue;
        }
        raster_faces++;
        int x0 = (int)floor(fmin(u[0], fmin(u[1], u[2])) - origin_u);
        int x1 = (int)ceil(fmax(u[0], fmax(u[1], u[2])) - origin_u) - 1;
        int y0 = (int)floor(fmin(v[0], fmin(v[1], v[2])) - origin_v);
        int y1 = (int)ceil(fmax(v[0], fmax(v[1], v[2])) - origin_v) - 1;
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= width) x1 = width - 1;
        if (y1 >= height) y1 = height - 1;
        int32_t chart = solution->vertex_chart[vi[0]];
        int chart_valid = chart >= 0 && (size_t)chart < solution->ncharts &&
            solution->vertex_chart[vi[1]] == chart &&
            solution->vertex_chart[vi[2]] == chart;
        for (int y = y0; y <= y1; y++) {
            double pv = origin_v + (double)y + 0.5;
            for (int x = x0; x <= x1; x++) {
                double pu = origin_u + (double)x + 0.5;
                double w0 = coverage_orient(
                    u[1], v[1], u[2], v[2], pu, pv) / area;
                double w1 = coverage_orient(
                    u[2], v[2], u[0], v[0], pu, pv) / area;
                double w2 = 1.0 - w0 - w1;
                if (w0 < -1.0e-10 || w1 < -1.0e-10 || w2 < -1.0e-10)
                    continue;
                int why;
                double residual = NAN;
                if (!solution->face_keep[face]) {
                    why = ARFT_COVERAGE_PREFIT_CULL;
                } else if (!chart_valid) {
                    why = ARFT_COVERAGE_INVALID_CHART;
                } else {
                    double rowf = (pv - fit->v0) / fit->dv;
                    double shift = coverage_dense_shift(
                        dense, fit->nrows, solution->ncharts, chart, rowf);
                    double source_p[2] = {0.0, 0.0};
                    for (int d = 0; d < 3; d++) {
                        double xyz = w0 * (double)ps->verts[(size_t)vi[0] * 3 + d] +
                                     w1 * (double)ps->verts[(size_t)vi[1] * 3 + d] +
                                     w2 * (double)ps->verts[(size_t)vi[2] * 3 + d];
                        double delta = xyz - set->axis_point[d];
                        source_p[0] += delta * set->basis0[d];
                        source_p[1] += delta * set->basis1[d];
                    }
                    why = coverage_classify(
                        set, opts, fit, pu + shift, pv, source_p, &residual);
                    if (isfinite(residual)) {
                        geometry_tests++;
                        geometry_sum2 += residual * residual;
                        if (residual > geometry_max) geometry_max = residual;
                    }
                    chart_reason[(size_t)chart * ARFT_COVERAGE_REASON_COUNT +
                                 (size_t)why]++;
                }
                layer_reason[why]++;
                size_t pixel = (size_t)y * (size_t)width + (size_t)x;
                if (layers[pixel] != UINT8_MAX) layers[pixel]++;
                if (coverage_reason_priority(why) >
                    coverage_reason_priority(reason[pixel]))
                    reason[pixel] = (uint8_t)why;
            }
        }
    }

    uint64_t observed_pixels = 0, multi_pixels = 0;
    for (size_t i = 0; i < pixels; i++) {
        if (layers[i] > 0) observed_pixels++;
        if (layers[i] > 1) multi_pixels++;
        display_reason[reason[i]]++;
    }
    char path[ARFT_PATH_CAP];
    int n = snprintf(path, sizeof path, "%s/ribbon_coverage_reason.tif", dir);
    if (n < 0 || (size_t)n >= sizeof path ||
        TiffIO_save(path, reason, 1, height, width) != 0)
        return -1;
    n = snprintf(path, sizeof path, "%s/ribbon_coverage_layers.tif", dir);
    if (n < 0 || (size_t)n >= sizeof path ||
        TiffIO_save(path, layers, 1, height, width) != 0)
        return -1;

    FILE *fp = open_out(dir, "ribbon_coverage_legend.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "code,reason\n");
    for (int i = 0; i < ARFT_COVERAGE_REASON_COUNT; i++)
        fprintf(fp, "%d,%s\n", i, coverage_reason_name(i));
    if (fclose(fp) != 0) return -1;

    fp = open_out(dir, "ribbon_coverage_by_chart.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "chart");
    for (int i = 1; i < ARFT_COVERAGE_REASON_COUNT; i++)
        fprintf(fp, ",%s", coverage_reason_name(i));
    fprintf(fp, "\n");
    for (size_t chart = 0; chart < solution->ncharts; chart++) {
        fprintf(fp, "%zu", chart);
        for (int i = 1; i < ARFT_COVERAGE_REASON_COUNT; i++)
            fprintf(fp, ",%" PRIu64,
                    chart_reason[chart * ARFT_COVERAGE_REASON_COUNT +
                                 (size_t)i]);
        fprintf(fp, "\n");
    }
    if (fclose(fp) != 0) return -1;

    uint64_t layer_total = 0, kept_layer_total = 0;
    for (int i = 1; i < ARFT_COVERAGE_REASON_COUNT; i++)
        layer_total += layer_reason[i];
    for (int i = ARFT_COVERAGE_INVALID_CHART;
         i < ARFT_COVERAGE_REASON_COUNT; i++)
        kept_layer_total += layer_reason[i];
    kept_layer_total += layer_reason[ARFT_COVERAGE_RIBBON];
    fp = open_out(dir, "ribbon_coverage_stats.json");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"domain\": \"discrete_source_uv\",\n"
        "  \"origin_u\": %.17g,\n"
        "  \"origin_v\": %.17g,\n"
        "  \"width\": %d,\n"
        "  \"height\": %d,\n"
        "  \"raster_faces\": %" PRIu64 ",\n"
        "  \"degenerate_faces\": %" PRIu64 ",\n"
        "  \"observed_union_pixels\": %" PRIu64 ",\n"
        "  \"multi_layer_union_pixels\": %" PRIu64 ",\n"
        "  \"observation_layer_pixels\": %" PRIu64 ",\n"
        "  \"kept_observation_layer_pixels\": %" PRIu64 ",\n"
        "  \"ribboned_layer_pixels\": %" PRIu64 ",\n"
        "  \"ribboned_fraction_all\": %.9g,\n"
        "  \"ribboned_fraction_kept\": %.9g,\n"
        "  \"geometry_tests\": %" PRIu64 ",\n"
        "  \"geometry_residual_rms\": %.9g,\n"
        "  \"geometry_residual_max\": %.9g,\n"
        "  \"layer_reason_counts\": {\n",
        origin_u, origin_v, width, height, raster_faces, degenerate_faces,
        observed_pixels, multi_pixels, layer_total, kept_layer_total,
        layer_reason[ARFT_COVERAGE_RIBBON],
        layer_total > 0 ? (double)layer_reason[ARFT_COVERAGE_RIBBON] /
                          (double)layer_total : 0.0,
        kept_layer_total > 0 ? (double)layer_reason[ARFT_COVERAGE_RIBBON] /
                               (double)kept_layer_total : 0.0,
        geometry_tests,
        geometry_tests > 0 ? sqrt(geometry_sum2 / (double)geometry_tests) : 0.0,
        geometry_max);
    for (int i = 1; i < ARFT_COVERAGE_REASON_COUNT; i++)
        fprintf(fp, "    \"%s\": %" PRIu64 "%s\n",
                coverage_reason_name(i), layer_reason[i],
                i + 1 < ARFT_COVERAGE_REASON_COUNT ? "," : "");
    fprintf(fp, "  },\n  \"display_pixel_reason_counts\": {\n");
    for (int i = 0; i < ARFT_COVERAGE_REASON_COUNT; i++)
        fprintf(fp, "    \"%s\": %" PRIu64 "%s\n",
                coverage_reason_name(i), display_reason[i],
                i + 1 < ARFT_COVERAGE_REASON_COUNT ? "," : "");
    fprintf(fp, "  }\n}\n");
    return fclose(fp) == 0 ? 0 : -1;
}


static int build_texture_delta_targets(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonTextureResult *phase,
    double pair_weight,
    double **out_delta,
    double **out_weight,
    size_t *out_count)
{
    if (arena == NULL || solution == NULL || set == NULL || phase == NULL ||
        out_delta == NULL || out_weight == NULL || out_count == NULL ||
        set->nrows < 2 || solution->ncharts == 0 ||
        set->nrows - 1 > SIZE_MAX / solution->ncharts ||
        !isfinite(pair_weight) || !(pair_weight > 0.0))
        return -1;
    size_t ncharts = solution->ncharts;
    size_t nedges = (set->nrows - 1) * ncharts;
    double *delta = (double *)ARENA_CALLOC(
        arena, nedges, sizeof(*delta));
    double *weight = (double *)ARENA_CALLOC(
        arena, nedges, sizeof(*weight));
    size_t count = 0;
    for (size_t i = 0; i < phase->npairs; i++) {
        const AtlasRibbonTexturePair *pair = &phase->pair[i];
        if (!pair->accepted) continue;
        if (pair->row0 < 0 || pair->row1 != pair->row0 + 1 ||
            (size_t)pair->row1 >= set->nrows || pair->chart < 0 ||
            (size_t)pair->chart >= ncharts ||
            !isfinite(pair->desired_shift_delta))
            return -1;
        size_t at = (size_t)pair->row0 * ncharts + (size_t)pair->chart;
        if (weight[at] > 0.0) return -1;
        delta[at] = pair->desired_shift_delta;
        weight[at] = pair_weight;
        count++;
    }
    *out_delta = delta;
    *out_weight = weight;
    *out_count = count;
    return count > 0 ? 0 : -1;
}

static int texture_field_motion(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const double *before,
    const double *after,
    double *out_rms,
    double *out_max)
{
    if (scratch == NULL || solution == NULL || set == NULL ||
        before == NULL || after == NULL || out_rms == NULL ||
        out_max == NULL || solution->ncharts == 0 ||
        set->nrows > SIZE_MAX / solution->ncharts)
        return -1;
    size_t nfield = set->nrows * solution->ncharts;
    uint8_t *seen = (uint8_t *)ARENA_CALLOC(scratch, nfield, 1);
    double sum2 = 0.0, maximum = 0.0;
    size_t count = 0;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *sample = &set->layer_sample[i];
        if (sample->row < 0 || (size_t)sample->row >= set->nrows ||
            sample->chart < 0 ||
            (size_t)sample->chart >= solution->ncharts)
            return -1;
        size_t at = (size_t)sample->row * solution->ncharts +
                    (size_t)sample->chart;
        if (seen[at]) continue;
        seen[at] = 1;
        double delta = after[at] - before[at];
        if (!isfinite(delta)) return -1;
        sum2 += delta * delta;
        if (fabs(delta) > maximum) maximum = fabs(delta);
        count++;
    }
    *out_rms = count > 0 ? sqrt(sum2 / (double)count) : 0.0;
    *out_max = maximum;
    return count > 0 ? 0 : -1;
}

static int run_texture_alignment(
    Arena_T texture_arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonTextureResult *phase,
    int rounds,
    int sweeps,
    double pair_weight,
    double anchor_weight,
    AtlasRibbonCollisionBound **collision_ledger,
    size_t *ncollision_ledger,
    size_t *collision_ledger_capacity,
    ArftUvAudit *uv_audit,
    Arena_T *uv_audit_arena,
    int *have_uv_audit,
    AtlasRibbonFitResult *fit,
    ArftTextureAlignStats *stats,
    ArftTextureAlignRound **out_round)
{
    if (texture_arena == NULL || ps == NULL || solution == NULL ||
        set == NULL || opts == NULL || phase == NULL || rounds < 1 ||
        sweeps < 1 || !(pair_weight > 0.0) || anchor_weight < 0.0 ||
        collision_ledger == NULL || ncollision_ledger == NULL ||
        collision_ledger_capacity == NULL || uv_audit == NULL ||
        uv_audit_arena == NULL || have_uv_audit == NULL || fit == NULL ||
        stats == NULL || out_round == NULL || !*have_uv_audit ||
        uv_audit->registered_hard_pairs != 0 || solution->ncharts == 0 ||
        set->nrows > SIZE_MAX / solution->ncharts)
        return -1;
    size_t nfield = set->nrows * solution->ncharts;
    double *edge_delta = NULL, *edge_weight = NULL;
    size_t target_pairs = 0;
    if (build_texture_delta_targets(
            texture_arena, solution, set, phase, pair_weight,
            &edge_delta, &edge_weight, &target_pairs) != 0)
        return -1;
    double *initial = (double *)ARENA_ALLOC(
        texture_arena, nfield * sizeof(*initial));
    double *anchor = (double *)ARENA_ALLOC(
        texture_arena, nfield * sizeof(*anchor));
    double *candidate = (double *)ARENA_ALLOC(
        texture_arena, nfield * sizeof(*candidate));
    double *row_trust = (double *)ARENA_ALLOC(
        texture_arena, set->nrows * sizeof(*row_trust));
    uint8_t *row_mark = (uint8_t *)ARENA_ALLOC(
        texture_arena, set->nrows * sizeof(*row_mark));
    memcpy(initial, fit->row_chart_shift, nfield * sizeof(*initial));
    ArftTextureAlignRound *round_stats =
        (ArftTextureAlignRound *)ARENA_CALLOC(
            texture_arena, (size_t)rounds, sizeof(*round_stats));
    memset(stats, 0, sizeof *stats);
    stats->rounds_requested = rounds;
    stats->sweeps = sweeps;
    stats->pair_weight = pair_weight;
    stats->anchor_weight = anchor_weight;

    for (int round = 0; round < rounds; round++) {
        ArftTextureAlignRound *r = &round_stats[round];
        r->round = round + 1;
        r->target_pairs = target_pairs;
        r->residual_rms_before = phase->residual_rms;
        r->residual_p95_before = phase->residual_p95;
        r->residual_max_before = phase->residual_max;
        memcpy(anchor, fit->row_chart_shift, nfield * sizeof(*anchor));
        if (*uv_audit_arena != NULL) Arena_dispose(uv_audit_arena);
        *have_uv_audit = 0;
        memset(uv_audit, 0, sizeof *uv_audit);

        Arena_T phase_arena = Arena_new();
        int phase_rc = phase_arena == NULL ? -1 :
            AtlasRibbonFit_refine_u_deltas(
                phase_arena, solution, set, opts,
                *collision_ledger, *ncollision_ledger,
                edge_delta, edge_weight, anchor, anchor_weight,
                sweeps, fit);
        if (phase_arena != NULL) Arena_dispose(&phase_arena);
        if (phase_rc != 0) {
            fprintf(stderr,
                    "[atlas_ribbon_fit] texture U step failed round=%d\n",
                    round + 1);
            return -1;
        }
        memcpy(candidate, fit->row_chart_shift,
               nfield * sizeof(*candidate));
        for (size_t row = 0; row < set->nrows; row++) row_trust[row] = 1.0;

        for (int attempt = 0;; attempt++) {
            Arena_T audit_arena = Arena_new();
            ArftUvAudit audit;
            memset(&audit, 0, sizeof audit);
            if (audit_arena == NULL ||
                build_uv_audit(
                    audit_arena, ps, solution, set, opts, fit, &audit) != 0) {
                if (audit_arena != NULL) Arena_dispose(&audit_arena);
                return -1;
            }
            r->audit_checks++;
            if (audit.registered_hard_pairs > r->hard_pairs_peak)
                r->hard_pairs_peak = audit.registered_hard_pairs;
            fprintf(stderr,
                    "[atlas_ribbon_fit] texture round=%d audit=%d exact=%zu "
                    "hard=%zu proposals=%zu/%zu ledger=%zu\n",
                    round + 1, attempt, audit.registered.npairs,
                    audit.registered_hard_pairs, audit.ncollision_bounds,
                    audit.collision_bound_face_pairs, *ncollision_ledger);
            if (audit.registered_hard_pairs == 0) {
                *uv_audit = audit;
                *uv_audit_arena = audit_arena;
                *have_uv_audit = 1;
                break;
            }
            if (attempt >= opts->collision_rounds ||
                audit.ncollision_bounds == 0) {
                fprintf(stderr,
                        "[atlas_ribbon_fit] texture alignment could not "
                        "restore an exact collision-free field\n");
                Arena_dispose(&audit_arena);
                return -1;
            }

            double relaxation = opts->collision_relaxation;
            if (opts->collision_polish_rounds > 0 &&
                attempt >= opts->collision_rounds -
                           opts->collision_polish_rounds)
                relaxation = 1.0;
            memset(row_mark, 0, set->nrows * sizeof(*row_mark));
            size_t rollback_radius = 2 + 2 * (size_t)attempt;
            if (attempt + 1 >= opts->collision_rounds)
                rollback_radius = set->nrows;
            for (size_t i = 0; i < audit.ncollision_bounds; i++) {
                int32_t center = audit.collision_bound[i].row;
                if (center < 0 || (size_t)center >= set->nrows) {
                    Arena_dispose(&audit_arena);
                    return -1;
                }
                size_t first = (size_t)center > rollback_radius
                             ? (size_t)center - rollback_radius : 0;
                size_t last = (size_t)center + rollback_radius < set->nrows
                            ? (size_t)center + rollback_radius
                            : set->nrows - 1;
                for (size_t row = first; row <= last; row++) row_mark[row] = 1;
            }
            size_t updates = 0;
            double minimum_trust = 1.0;
            for (size_t row = 0; row < set->nrows; row++) {
                if (row_mark[row]) {
                    double next = row_trust[row] * (1.0 - relaxation);
                    if (next < 1.0e-6) next = 0.0;
                    if (next < row_trust[row] - 1.0e-12) {
                        row_trust[row] = next;
                        updates++;
                    }
                }
                if (row_trust[row] < minimum_trust)
                    minimum_trust = row_trust[row];
            }
            fprintf(stderr,
                    "[atlas_ribbon_fit] texture rollback=%d relax=%.3g "
                    "rows=%zu min_trust=%.4g\n",
                    attempt, relaxation, updates, minimum_trust);
            r->ledger_updates = updates > SIZE_MAX - r->ledger_updates
                              ? SIZE_MAX : r->ledger_updates + updates;
            Arena_dispose(&audit_arena);
            if (updates == 0) {
                fprintf(stderr,
                        "[atlas_ribbon_fit] texture rollback stalled\n");
                return -1;
            }
            for (size_t row = 0; row < set->nrows; row++)
            for (size_t chart = 0; chart < solution->ncharts; chart++) {
                size_t at = row * solution->ncharts + chart;
                fit->row_chart_shift[at] =
                    anchor[at] + row_trust[row] *
                    (candidate[at] - anchor[at]);
            }
            Arena_T refresh_arena = Arena_new();
            int refresh_rc = refresh_arena == NULL ? -1 :
                AtlasRibbonFit_refresh_u_field(
                    refresh_arena, solution, set, opts, fit);
            if (refresh_arena != NULL) Arena_dispose(&refresh_arena);
            if (refresh_rc != 0) {
                fprintf(stderr,
                        "[atlas_ribbon_fit] texture rollback refresh failed "
                        "round=%d attempt=%d\n", round + 1, attempt);
                return -1;
            }
        }

        Arena_T metric_arena = Arena_new();
        if (metric_arena == NULL ||
            AtlasRibbonTexture_refresh(
                metric_arena, solution, set, fit, phase) != 0 ||
            texture_field_motion(
                metric_arena, solution, set, anchor,
                fit->row_chart_shift, &r->step_rms, &r->step_max) != 0 ||
            texture_field_motion(
                metric_arena, solution, set, initial,
                fit->row_chart_shift,
                &r->cumulative_rms, &r->cumulative_max) != 0) {
            if (metric_arena != NULL) Arena_dispose(&metric_arena);
            return -1;
        }
        Arena_dispose(&metric_arena);
        r->residual_rms_after = phase->residual_rms;
        r->residual_p95_after = phase->residual_p95;
        r->residual_max_after = phase->residual_max;
        stats->rounds_completed = round + 1;
        fprintf(stderr,
                "[atlas_ribbon_fit] texture round=%d residual=%.4g->%.4g "
                "p95=%.4g->%.4g step=%.4g/%.4g cumulative=%.4g/%.4g "
                "hard_peak=%zu ledger_updates=%zu\n",
                round + 1, r->residual_rms_before, r->residual_rms_after,
                r->residual_p95_before, r->residual_p95_after,
                r->step_rms, r->step_max,
                r->cumulative_rms, r->cumulative_max,
                r->hard_pairs_peak, r->ledger_updates);
    }
    *out_round = round_stats;
    return 0;
}
static void json_double(FILE *fp, double value)
{
    if (isfinite(value)) fprintf(fp, "%.17g", value);
    else fputs("null", fp);
}

static void json_metric(FILE *fp, const char *name, double value, int comma)
{
    fprintf(fp, "  \"%s\": ", name);
    json_double(fp, value);
    fprintf(fp, "%s\n", comma ? "," : "");
}

static int write_texture_phase(
    const char *dir,
    const AtlasRibbonTextureOptions *opts,
    const AtlasRibbonTextureResult *phase)
{
    FILE *fp = open_out(dir, "u_texture_phase_pairs.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "row0,row1,chart,group,winding,lag_columns,samples,accepted,"
                "correlation,zero_correlation,margin,desired_shift_delta,"
                "current_shift_delta,residual\n");
    for (size_t i = 0; i < phase->npairs; i++) {
        const AtlasRibbonTexturePair *p = &phase->pair[i];
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%u,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                p->row0, p->row1, p->chart, p->group, p->winding,
                p->lag_columns, p->samples, p->accepted,
                p->correlation, p->zero_correlation, p->margin,
                p->desired_shift_delta, p->current_shift_delta, p->residual);
    }
    if (fclose(fp) != 0) return -1;
    fp = open_out(dir, "u_texture_phase_stats.json");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"search_columns\": %d,\n"
        "  \"min_samples\": %zu,\n"
        "  \"min_correlation\": %.17g,\n"
        "  \"min_margin\": %.17g,\n"
        "  \"pairs\": %zu,\n"
        "  \"accepted_pairs\": %zu,\n"
        "  \"sampled_layer_points\": %zu,\n"
        "  \"missing_layer_points\": %zu,\n"
        "  \"cubes_loaded\": %d,\n"
        "  \"cubes_missing\": %d,\n",
        opts->search_columns, opts->min_samples,
        opts->min_correlation, opts->min_margin,
        phase->npairs, phase->accepted_pairs,
        phase->sampled_layer_points, phase->missing_layer_points,
        phase->cubes_loaded, phase->cubes_missing);
    json_metric(fp, "accepted_correlation_mean",
                phase->accepted_correlation_mean, 1);
    json_metric(fp, "accepted_margin_mean", phase->accepted_margin_mean, 1);
    json_metric(fp, "current_delta_rms", phase->current_delta_rms, 1);
    json_metric(fp, "desired_delta_rms", phase->desired_delta_rms, 1);
    json_metric(fp, "residual_rms", phase->residual_rms, 1);
    json_metric(fp, "residual_p95", phase->residual_p95, 1);
    json_metric(fp, "residual_max", phase->residual_max, 0);
    fputs("}\n", fp);
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_texture_alignment(
    const char *dir,
    const ArftTextureAlignStats *stats,
    const ArftTextureAlignRound *round)
{
    if (stats == NULL || stats->rounds_completed < 1 || round == NULL)
        return -1;
    FILE *fp = open_out(dir, "u_texture_alignment_rounds.csv");
    if (fp == NULL) return -1;
    fprintf(fp,
        "round,target_pairs,audit_checks,hard_pairs_peak,rollback_row_updates,"
        "residual_rms_before,residual_p95_before,residual_max_before,"
        "residual_rms_after,residual_p95_after,residual_max_after,"
        "step_rms,step_max,cumulative_rms,cumulative_max\n");
    size_t total_updates = 0;
    for (int i = 0; i < stats->rounds_completed; i++) {
        const ArftTextureAlignRound *r = &round[i];
        total_updates = r->ledger_updates > SIZE_MAX - total_updates
            ? SIZE_MAX
            : total_updates + r->ledger_updates;
        fprintf(fp,
            "%d,%zu,%zu,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,%.17g\n",
            r->round, r->target_pairs, r->audit_checks,
            r->hard_pairs_peak, r->ledger_updates,
            r->residual_rms_before, r->residual_p95_before,
            r->residual_max_before, r->residual_rms_after,
            r->residual_p95_after, r->residual_max_after,
            r->step_rms, r->step_max,
            r->cumulative_rms, r->cumulative_max);
    }
    if (fclose(fp) != 0) return -1;
    const ArftTextureAlignRound *last =
        &round[stats->rounds_completed - 1];
    fp = open_out(dir, "u_texture_alignment_stats.json");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"rounds_requested\": %d,\n"
        "  \"rounds_completed\": %d,\n"
        "  \"sweeps_per_round\": %d,\n"
        "  \"pair_weight\": %.17g,\n"
        "  \"anchor_weight\": %.17g,\n"
        "  \"target_pairs\": %zu,\n"
        "  \"total_rollback_row_updates\": %zu,\n",
        stats->rounds_requested, stats->rounds_completed, stats->sweeps,
        stats->pair_weight, stats->anchor_weight, last->target_pairs,
        total_updates);
    json_metric(fp, "initial_residual_rms", round[0].residual_rms_before, 1);
    json_metric(fp, "initial_residual_p95", round[0].residual_p95_before, 1);
    json_metric(fp, "final_residual_rms", last->residual_rms_after, 1);
    json_metric(fp, "final_residual_p95", last->residual_p95_after, 1);
    json_metric(fp, "final_residual_max", last->residual_max_after, 1);
    json_metric(fp, "cumulative_motion_rms", last->cumulative_rms, 1);
    json_metric(fp, "cumulative_motion_max", last->cumulative_max, 0);
    fputs("}\n", fp);
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_stats_json(
    const char *dir,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    size_t ribbon_faces,
    double seconds)
{
    FILE *fp = open_out(dir, "ribbon_fit_stats.json");
    if (fp == NULL) return -1;
    size_t kept_faces = 0;
    for (size_t f = 0; f < solution->nfaces; f++)
        if (solution->face_keep[f]) kept_faces++;
    fprintf(fp,
        "{\n"
        "  \"mode\": \"%s\",\n"
        "  \"piece_vertices\": %zu,\n"
        "  \"piece_faces\": %zu,\n"
        "  \"kept_faces\": %zu,\n"
        "  \"charts\": %zu,\n"
        "  \"checkpoint_fingerprint\": \"0x%016" PRIx64 "\",\n"
        "  \"slice_spacing\": %.17g,\n"
        "  \"observation_u_spacing\": %.17g,\n"
        "  \"fit_u_spacing\": %.17g,\n"
        "  \"local_xyz_tolerance\": %.17g,\n"
        "  \"tangent_dot_min\": %.17g,\n"
        "  \"max_fill_u\": %.17g,\n"
        "  \"max_bridge_stretch\": %.17g,\n"
        "  \"bridge_slack\": %.17g,\n"
        "  \"lambda_position\": %.17g,\n"
        "  \"lambda_smooth\": %.17g,\n"
        "  \"lambda_tangent\": %.17g,\n"
        "  \"input_faces\": %zu,\n"
        "  \"sliced_faces\": %zu,\n"
        "  \"slice_segments\": %zu,\n"
        "  \"degenerate_segments\": %zu,\n"
        "  \"invalid_chart_faces\": %zu,\n"
        "  \"observations\": %zu,\n"
        "  \"targets\": %zu,\n"
        "  \"accepted_targets\": %zu,\n"
        "  \"duplicate_targets\": %zu,\n"
        "  \"conflict_targets\": %zu,\n"
        "  \"conflict_clusters\": %zu,\n"
        "  \"known_topology_conflicts\": %zu,\n"
        "  \"fit_rows\": %zu,\n"
        "  \"fit_columns\": %zu,\n"
        "  \"fitted_nodes\": %zu,\n"
        "  \"supported_nodes\": %zu,\n"
        "  \"observed_nodes\": %zu,\n"
        "  \"chart_interp_nodes\": %zu,\n"
        "  \"seam_interp_nodes\": %zu,\n"
        "  \"spiral_fill_nodes\": %zu,\n"
        "  \"curve_segments\": %zu,\n"
        "  \"ribbon_faces\": %zu,\n"
        "  \"row_crossings\": %zu,\n"
        "  \"crossing_cuts\": %zu,\n"
        "  \"u_gap_cuts\": %zu,\n"
        "  \"metric_jump_cuts\": %zu,\n"
        "  \"direct_spans\": %zu,\n"
        "  \"spiral_spans\": %zu,\n"
        "  \"topology_cuts\": %zu,\n"
        "  \"winding_rows\": %zu,\n"
        "  \"winding_phase_backsteps\": %zu,\n"
        "  \"winding_radius_order_tests\": %zu,\n"
        "  \"winding_radius_order_violations\": %zu,\n"
        "  \"parameter_samples\": %zu,\n"
        "  \"chamfer_samples\": %zu,\n"
        "  \"speed_samples\": %zu,\n",
        mode_name(opts->mode), ps->nv, ps->nf,
        kept_faces, solution->ncharts,
        solution->piece_fingerprint,
        opts->slice_spacing, opts->observation_u_spacing,
        opts->fit_u_spacing, opts->local_xyz_tolerance,
        opts->tangent_dot_min, opts->max_fill_u,
        opts->max_bridge_stretch, opts->bridge_slack,
        opts->lambda_position, opts->lambda_smooth,
        opts->lambda_tangent,
        set->input_faces, set->sliced_faces, set->slice_segments,
        set->degenerate_segments, set->invalid_chart_faces,
        set->nobservation, set->ntarget, set->accepted_targets,
        set->duplicate_targets, set->conflict_targets,
        set->conflict_clusters, set->known_topology_conflicts,
        fit->nrows, fit->ncolumns, fit->fitted_nodes,
        fit->supported_nodes, fit->observed_nodes,
        fit->chart_interp_nodes, fit->seam_interp_nodes,
        fit->spiral_fill_nodes, fit->curve_segments, ribbon_faces,
        fit->row_crossings, fit->crossing_cuts,
        fit->u_gap_cuts, fit->metric_jump_cuts,
        fit->direct_spans, fit->spiral_spans, fit->topology_cuts,
        fit->winding_rows,
        fit->winding_phase_backsteps,
        fit->winding_radius_order_tests,
        fit->winding_radius_order_violations,
        fit->parameter_samples, fit->chamfer_samples,
        fit->speed_samples);
    fprintf(fp,
        "  \"lambda_register_vertical\": %.17g,\n"
        "  \"register_sweeps_requested\": %d,\n"
        "  \"collision_relaxation\": %.17g,\n"
        "  \"collision_rounds_requested\": %d,\n"
        "  \"collision_sweeps_requested\": %d,\n"
        "  \"collision_polish_rounds\": %d,\n",
        opts->lambda_register_vertical, opts->register_sweeps,
        opts->collision_relaxation, opts->collision_rounds,
        opts->collision_sweeps, opts->collision_polish_rounds);
    fprintf(fp,
        "  \"layer_samples\": %zu,\n"
        "  \"remote_layer_samples\": %zu,\n"
        "  \"register_winding_direction\": %d,\n"
        "  \"register_iterations\": %zu,\n"
        "  \"register_constraints\": %zu,\n"
        "  \"register_remote_constraints\": %zu,\n"
        "  \"register_metric_constraints\": %zu,\n"
        "  \"register_order_constraints\": %zu,\n"
        "  \"register_unresolved\": %zu,\n"
        "  \"register_inversions_before\": %zu,\n"
        "  \"register_inversions_after\": %zu,\n"
        "  \"register_local_inversions\": %zu,\n",
        set->nlayer_samples, set->remote_layer_samples,
        fit->register_winding_direction, fit->register_iterations,
        fit->register_constraints,
        fit->register_constraints - fit->register_metric_constraints -
            fit->register_order_constraints,
        fit->register_metric_constraints, fit->register_order_constraints,
        fit->register_unresolved, fit->register_inversions_before,
        fit->register_inversions_after, fit->register_local_inversions);
    fprintf(fp,
        "  \"register_qp_failures\": %zu,\n"
        "  \"register_present_chart_rows\": %zu,\n"
        "  \"register_smoothing_sweeps\": %zu,\n"
        "  \"register_collision_rounds\": %zu,\n"
        "  \"register_collision_constraints\": %zu,\n"
        "  \"register_collision_face_pairs\": %zu,\n",
        fit->register_qp_failures, fit->register_present_chart_rows,
        fit->register_smoothing_sweeps,
        fit->register_collision_rounds,
        fit->register_collision_constraints,
        fit->register_collision_face_pairs);
    json_metric(fp, "register_smoothing_max_delta",
                fit->register_smoothing_max_delta, 1);
    json_metric(fp, "register_width_before", fit->register_width_before, 1);
    json_metric(fp, "register_width_after", fit->register_width_after, 1);
    json_metric(fp, "register_shift_rms", fit->register_shift_rms, 1);
    json_metric(fp, "register_shift_max", fit->register_shift_max, 1);
    json_metric(fp, "register_max_required_arc",
                fit->register_max_required_arc, 1);
    json_metric(fp, "register_vertical_shift_rms_before",
                fit->register_vertical_shift_rms_before, 1);
    json_metric(fp, "register_vertical_shift_max_before",
                fit->register_vertical_shift_max_before, 1);
    json_metric(fp, "register_vertical_shift_rms",
                fit->register_vertical_shift_rms, 1);
    json_metric(fp, "register_vertical_shift_max",
                fit->register_vertical_shift_max, 1);
    json_metric(fp, "register_local_gap_mean",
                fit->register_local_gap_mean, 1);
    json_metric(fp, "register_local_gap_max", fit->register_local_gap_max, 1);
    json_metric(fp, "register_row_width_ratio_mean",
                fit->register_row_width_ratio_mean, 1);
    json_metric(fp, "parameter_mean", fit->parameter_mean, 1);
    json_metric(fp, "parameter_rms", fit->parameter_rms, 1);
    json_metric(fp, "parameter_p50", fit->parameter_p50, 1);
    json_metric(fp, "parameter_p90", fit->parameter_p90, 1);
    json_metric(fp, "parameter_p95", fit->parameter_p95, 1);
    json_metric(fp, "parameter_p99", fit->parameter_p99, 1);
    json_metric(fp, "parameter_max", fit->parameter_max, 1);
    json_metric(fp, "directed_chamfer_mean", fit->chamfer_mean, 1);
    json_metric(fp, "directed_chamfer_rms", fit->chamfer_rms, 1);
    json_metric(fp, "directed_chamfer_p50", fit->chamfer_p50, 1);
    json_metric(fp, "directed_chamfer_p90", fit->chamfer_p90, 1);
    json_metric(fp, "directed_chamfer_p95", fit->chamfer_p95, 1);
    json_metric(fp, "directed_chamfer_p99", fit->chamfer_p99, 1);
    json_metric(fp, "directed_hausdorff", fit->directed_hausdorff, 1);
    json_metric(fp, "speed_mean", fit->speed_mean, 1);
    json_metric(fp, "speed_rms_error", fit->speed_rms_error, 1);
    json_metric(fp, "speed_p95_error", fit->speed_p95_error, 1);
    json_metric(fp, "speed_max_error", fit->speed_max_error, 1);
    json_metric(fp, "seconds", seconds, 0);
    fputs("}\n", fp);
    return fclose(fp) == 0 ? 0 : -1;
}

static void usage(const char *program)
{
    fprintf(stderr,
        "usage: %s <placed_dir> <atlas_solution.bin> <out_dir> [options]\n"
        "       %s --selftest\n"
        "options:\n"
        "  --mode f0|f1|f2|register|field|smooth|collision|ribbon\n"
        "                            stage to stop after (default ribbon =\n"
        "                            the whole path; earlier stages ablate)\n"
        "  --slice-spacing F        constant-v row spacing (default 4)\n"
        "  --observation-u F         exact section sample spacing (default 2)\n"
        "  --fit-u F                 ribbon control spacing (default 16)\n"
        "  --local-xyz F             local duplicate radius (default 8)\n"
        "  --tangent-dot F           local tangent compatibility (default .5)\n"
        "  --max-fill-u F            largest bridged U gap (default 512)\n"
        "  --max-bridge-stretch F    reject metrically impossible chords (1.25)\n"
        "  --bridge-slack F          bridge metric allowance (default 4)\n"
        "  --lambda-position F       fixed-parameter observation weight (1)\n"
        "  --lambda-smooth F         zero-curvature weight (8)\n"
        "  --lambda-tangent F        increasing-U tangent weight (4)\n"
        "  --lambda-register-v F     vertical chart-shift weight (16)\n"
        "  --register-sweeps N       bidirectional feasible sweeps (6)\n"
        "  --collision-rounds N      exact-audit continuation rounds (8)\n"
        "  --collision-polish N      final full-escape rounds (0)\n"
        "  --collision-sweeps N      sweeps after each collision round (4)\n"
        "  --collision-relaxation F  fraction of escape applied per round (.35)\n"
        "  --raw DIR                 audit adjacent-row RAW texture phase\n"
        "  --texture-search N        phase search in observation columns (8)\n"
        "  --texture-min-samples N   minimum matched samples (24)\n"
        "  --texture-min-corr F      accepted Pearson correlation (.20)\n"
        "  --texture-min-margin F    best-vs-other peak margin (.02)\n"
        "  --texture-normal-range F  RAW normal-max reach (2)\n"
        "  --texture-normal-steps N  RAW normal-max taps (5)\n"
        "  --texture-rounds N        soft collision-locked U rounds (0)\n"
        "  --texture-sweeps N        row sweeps per U round (4)\n"
        "  --texture-weight F        accepted phase-pair weight (1024)\n"
        "  --texture-anchor F        per-round proximal weight (4096)\n",
        program, program);
}

static int parse_mode(const char *value, AtlasRibbonFitMode *mode)
{
    if (strcmp(value, "f0") == 0 || strcmp(value, "observations") == 0)
        *mode = ATLAS_RIBBON_OBSERVATIONS;
    else if (strcmp(value, "f1") == 0 || strcmp(value, "position") == 0)
        *mode = ATLAS_RIBBON_POSITION;
    else if (strcmp(value, "f2") == 0 || strcmp(value, "tangent") == 0)
        *mode = ATLAS_RIBBON_TANGENT;
    else if (strcmp(value, "register") == 0 || strcmp(value, "f3") == 0)
        *mode = ATLAS_RIBBON_REGISTER_ONLY;
    else if (strcmp(value, "field") == 0 || strcmp(value, "f4") == 0)
        *mode = ATLAS_RIBBON_REGISTER_FIELD_ONLY;
    else if (strcmp(value, "smooth") == 0 || strcmp(value, "f5") == 0)
        *mode = ATLAS_RIBBON_REGISTER_FIELD_SMOOTH_ONLY;
    else if (strcmp(value, "collision") == 0 || strcmp(value, "f6") == 0)
        *mode = ATLAS_RIBBON_REGISTER_COLLISION_ONLY;
    else if (strcmp(value, "ribbon") == 0 || strcmp(value, "f7") == 0)
        *mode = ATLAS_RIBBON_REGISTERED_RIBBON;
    else
        return -1;
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        const char *path =
            "output/_selftest_atlas_ribbon/atlas_solution.bin";
        Arena_T arena = Arena_new();
        int failures = arena == NULL ? 1 : AtlasSolution_selftest(arena, path);
        if (arena != NULL) Arena_dispose(&arena);
        failures += AtlasRibbonFit_selftest();
        failures += AtlasRibbonTexture_selftest();
        failures += AtlasOverlapAudit_selftest();
        fprintf(stderr, "[atlas_ribbon_fit selftest] %s (%d failures)\n",
                failures == 0 ? "ALL PASS" : "FAILURES", failures);
        return failures == 0 ? 0 : 1;
    }
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }
    const char *placed_dir = argv[1];
    const char *solution_path = argv[2];
    const char *out_dir = argv[3];
    AtlasRibbonTextureOptions texture_opts;
    AtlasRibbonTextureOptions_default(&texture_opts);
    const char *raw_dir = NULL;
    int texture_rounds = 0;
    int texture_sweeps = 4;
    double texture_weight = 1024.0;
    double texture_anchor_weight = 4096.0;
    AtlasRibbonFitOptions opts;
    AtlasRibbonFitOptions_default(&opts);
    /* The library default stops at the tangent fit, which is the right base
     * for a caller assembling stages itself. Running THIS tool means asking
     * for a ribbon, so the front end defaults to the whole path; --mode still
     * selects any earlier stage for ablation. */
    opts.mode = ATLAS_RIBBON_REGISTERED_RIBBON;
    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--raw") == 0 && i + 1 < argc) {
            raw_dir = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "--texture-min-samples") == 0 && i + 1 < argc) {
            char *end = NULL;
            unsigned long value = strtoul(argv[++i], &end, 10);
            if (end == argv[i] || *end != '\0' || value < 3 || value > 100000) {
                fprintf(stderr, "%s: invalid --texture-min-samples\n", argv[0]);
                return 1;
            }
            texture_opts.min_samples = (size_t)value;
            continue;
        }
        if (strcmp(argv[i], "--mode") == 0 && i + 1 < argc) {
            if (parse_mode(argv[++i], &opts.mode) != 0) {
                fprintf(stderr, "%s: invalid mode %s\n", argv[0], argv[i]);
                return 1;
            }
            continue;
        }
        int *integer_destination = NULL;
        if (strcmp(argv[i], "--register-sweeps") == 0)
            integer_destination = &opts.register_sweeps;
        else if (strcmp(argv[i], "--collision-rounds") == 0)
            integer_destination = &opts.collision_rounds;
        else if (strcmp(argv[i], "--collision-sweeps") == 0)
            integer_destination = &opts.collision_sweeps;
        else if (strcmp(argv[i], "--collision-polish") == 0)
            integer_destination = &opts.collision_polish_rounds;
        else if (strcmp(argv[i], "--texture-search") == 0)
            integer_destination = &texture_opts.search_columns;
        else if (strcmp(argv[i], "--texture-normal-steps") == 0)
            integer_destination = &texture_opts.normal_steps;
        else if (strcmp(argv[i], "--texture-rounds") == 0)
            integer_destination = &texture_rounds;
        else if (strcmp(argv[i], "--texture-sweeps") == 0)
            integer_destination = &texture_sweeps;
        if (integer_destination != NULL) {
            const char *option = argv[i];
            if (i + 1 >= argc) {
                fprintf(stderr, "%s: incomplete option %s\n", argv[0], argv[i]);
                return 1;
            }
            char *end = NULL;
            long value = strtol(argv[++i], &end, 10);
            if (end == argv[i] || *end != '\0' || value < 0 || value > 100) {
                fprintf(stderr, "%s: invalid value for %s\n", argv[0], option);
                return 1;
            }
            *integer_destination = (int)value;
            continue;
        }
        double *destination = NULL;
        if (strcmp(argv[i], "--slice-spacing") == 0)
            destination = &opts.slice_spacing;
        else if (strcmp(argv[i], "--observation-u") == 0)
            destination = &opts.observation_u_spacing;
        else if (strcmp(argv[i], "--fit-u") == 0)
            destination = &opts.fit_u_spacing;
        else if (strcmp(argv[i], "--local-xyz") == 0)
            destination = &opts.local_xyz_tolerance;
        else if (strcmp(argv[i], "--tangent-dot") == 0)
            destination = &opts.tangent_dot_min;
        else if (strcmp(argv[i], "--max-fill-u") == 0)
            destination = &opts.max_fill_u;
        else if (strcmp(argv[i], "--lambda-position") == 0)
            destination = &opts.lambda_position;
        else if (strcmp(argv[i], "--lambda-smooth") == 0)
            destination = &opts.lambda_smooth;
        else if (strcmp(argv[i], "--max-bridge-stretch") == 0)
            destination = &opts.max_bridge_stretch;
        else if (strcmp(argv[i], "--bridge-slack") == 0)
            destination = &opts.bridge_slack;
        else if (strcmp(argv[i], "--lambda-register-v") == 0)
            destination = &opts.lambda_register_vertical;
        else if (strcmp(argv[i], "--lambda-tangent") == 0)
            destination = &opts.lambda_tangent;
        else if (strcmp(argv[i], "--collision-relaxation") == 0)
            destination = &opts.collision_relaxation;
        else if (strcmp(argv[i], "--texture-min-corr") == 0)
            destination = &texture_opts.min_correlation;
        else if (strcmp(argv[i], "--texture-min-margin") == 0)
            destination = &texture_opts.min_margin;
        else if (strcmp(argv[i], "--texture-normal-range") == 0)
            destination = &texture_opts.normal_range;
        else if (strcmp(argv[i], "--texture-weight") == 0)
            destination = &texture_weight;
        else if (strcmp(argv[i], "--texture-anchor") == 0)
            destination = &texture_anchor_weight;
        if (destination == NULL || i + 1 >= argc) {
            fprintf(stderr, "%s: unknown or incomplete option %s\n",
                    argv[0], argv[i]);
            return 1;
        }
        char *end = NULL;
        *destination = strtod(argv[++i], &end);
        if (end == argv[i] || *end != '\0' || !isfinite(*destination)) {
            fprintf(stderr, "%s: invalid value for %s\n",
                    argv[0], argv[i - 1]);
            return 1;
        }
    }

    char probe[ARFT_PATH_CAP];

    if (texture_rounds > 0 &&
        (raw_dir == NULL || opts.mode != ATLAS_RIBBON_REGISTERED_RIBBON ||
         texture_sweeps < 1 || !(texture_weight > 0.0) ||
         texture_anchor_weight < 0.0)) {
        fprintf(stderr,
                "%s: texture alignment requires ribbon mode, --raw, positive "
                "rounds/sweeps/weight, and a nonnegative anchor\n", argv[0]);
        return 1;
    }
    int probe_n = snprintf(probe, sizeof probe, "%s/.", out_dir);
    if (probe_n < 0 || (size_t)probe_n >= sizeof probe ||
        ves_ensure_parent_dir(probe) != 0) {
        fprintf(stderr, "%s: cannot create %s\n", argv[0], out_dir);
        return 1;
    }

    double start = ves_clock_sec();
    Arena_T arena = Arena_new();
    if (arena == NULL) return 1;
    PieceSet ps;
    if (PieceSet_build(arena, placed_dir, &ps) != 0) {
        fprintf(stderr, "%s: cannot load placed dir %s\n",
                argv[0], placed_dir);
        Arena_dispose(&arena);
        return 1;
    }
    ScaffoldCalib cal;
    if (Scaffold_read_calib(placed_dir, &cal) != 0) {
        fprintf(stderr, "%s: missing placed calibration\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    AtlasSolution solution;
    if (AtlasSolution_read(arena, solution_path, &ps, &solution) != 0) {
        fprintf(stderr, "%s: checkpoint does not match placed input: %s\n",
                argv[0], solution_path);
        Arena_dispose(&arena);
        return 1;
    }
    /* The tabu atlas is overlap-clean but globally wind-scrambled; ranking
     * charts by their continuous wind makes the register's one-sided rank
     * bounds agree with physical order instead of fighting it. */
    int wind_ranked = AtlasSolution_rank_by_wind(arena, &solution);
    if (wind_ranked < 0) {
        fprintf(stderr, "%s: checkpoint wind field invalid: %s\n",
                argv[0], solution_path);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr, "[atlas_ribbon_fit] chart ranks: %s\n",
            wind_ranked ? "wind order" : "producer relayout order");
    size_t kept_faces = 0;
    for (size_t f = 0; f < solution.nfaces; f++)
        if (solution.face_keep[f]) kept_faces++;
    fprintf(stderr,
            "[atlas_ribbon_fit] mode=%s cubes=%zu vertices=%zu faces=%zu "
            "kept=%zu charts=%zu\n",
            mode_name(opts.mode), ps.n_cubes, ps.nv, ps.nf,
            kept_faces, solution.ncharts);

    AtlasRibbonObservationSet observations;
    if (AtlasRibbonFit_build_observations(
            arena, &ps, &solution, &cal, &opts, &observations) != 0) {
        fprintf(stderr, "%s: observation extraction failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr,
            "[atlas_ribbon_fit] observations=%zu targets=%zu layers=%zu "
            "remote_layers=%zu accepted=%zu duplicates=%zu conflicts=%zu "
            "known_topology=%zu\n",
            observations.nobservation, observations.ntarget,
            observations.nlayer_samples, observations.remote_layer_samples,
            observations.accepted_targets, observations.duplicate_targets,
            observations.conflict_targets,
            observations.known_topology_conflicts);

    AtlasRibbonFitResult fit;
    if (AtlasRibbonFit_solve(
            arena, &ps, &solution, &observations, &opts, &fit) != 0) {
        fprintf(stderr, "%s: ribbon solve failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    ArftUvAudit uv_audit;
    int have_uv_audit = 0;
    Arena_T uv_audit_arena = NULL;
    AtlasRibbonCollisionBound *collision_ledger = NULL;
    size_t ncollision_ledger = 0, collision_ledger_capacity = 0;
    memset(&uv_audit, 0, sizeof uv_audit);
    if (opts.mode == ATLAS_RIBBON_REGISTER_FIELD_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTER_FIELD_SMOOTH_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTER_COLLISION_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON) {
        for (int round = 0;; round++) {
            Arena_T round_arena = Arena_new();
            ArftUvAudit round_audit;
            memset(&round_audit, 0, sizeof round_audit);
            if (round_arena == NULL ||
                build_uv_audit(round_arena, &ps, &solution, &observations,
                               &opts, &fit, &round_audit) != 0) {
                fprintf(stderr, "%s: exact registered U(v) audit failed\n",
                        argv[0]);
                if (round_arena != NULL) Arena_dispose(&round_arena);
                free(collision_ledger);
                Arena_dispose(&arena);
                return 1;
            }
            int stop = (opts.mode != ATLAS_RIBBON_REGISTER_COLLISION_ONLY &&
                        opts.mode != ATLAS_RIBBON_REGISTERED_RIBBON) ||
                       round_audit.registered_hard_pairs == 0 ||
                       round_audit.ncollision_bounds == 0 ||
                       round >= opts.collision_rounds;
            size_t ledger_updates = 0;
            double round_relaxation = opts.collision_relaxation;
            if (opts.collision_polish_rounds > 0 &&
                round >= opts.collision_rounds - opts.collision_polish_rounds)
                round_relaxation = 1.0;
            if (!stop && merge_collision_ledger(
                    &collision_ledger, &ncollision_ledger,
                    &collision_ledger_capacity,
                    round_audit.collision_bound,
                    round_audit.ncollision_bounds,
                    round_relaxation, &ledger_updates) != 0) {
                fprintf(stderr, "%s: collision ledger merge failed\n", argv[0]);
                Arena_dispose(&round_arena);
                free(collision_ledger);
                Arena_dispose(&arena);
                return 1;
            }
            size_t nblock_bounds = 0;
            if (!stop && round_relaxation >= 1.0) {
                AtlasRibbonCollisionBound *block_bound = NULL;
                size_t block_updates = 0;
                if (build_collision_block_bounds(
                        round_arena, &solution, &observations, &fit,
                        &round_audit, &block_bound, &nblock_bounds) != 0 ||
                    merge_collision_ledger(
                        &collision_ledger, &ncollision_ledger,
                        &collision_ledger_capacity,
                        block_bound, nblock_bounds, 1.0,
                        &block_updates) != 0) {
                    fprintf(stderr,
                            "%s: collision block-ledger merge failed\n",
                            argv[0]);
                    Arena_dispose(&round_arena);
                    free(collision_ledger);
                    Arena_dispose(&arena);
                    return 1;
                }
                ledger_updates = block_updates > SIZE_MAX - ledger_updates
                               ? SIZE_MAX : ledger_updates + block_updates;
            }
            fprintf(stderr,
                    "[atlas_ribbon_fit] collision round=%d exact=%zu "
                    "hard=%zu proposals=%zu/%zu blocks=%zu ledger=%zu updates=%zu relax=%.3g\n",
                    round, round_audit.registered.npairs,
                    round_audit.registered_hard_pairs,
                    round_audit.ncollision_bounds,
                    round_audit.collision_bound_face_pairs,
                    nblock_bounds, ncollision_ledger, ledger_updates, round_relaxation);
            if (!stop && ledger_updates == 0) stop = 1;
            if (stop) {
                uv_audit = round_audit;
                uv_audit_arena = round_arena;
                have_uv_audit = 1;
                break;
            }
            Arena_dispose(&round_arena);
            Arena_T refine_arena = Arena_new();
            int refine_rc = refine_arena == NULL ? -1 :
                AtlasRibbonFit_refine_u_collisions(
                    refine_arena, &solution, &observations, &opts,
                    collision_ledger, ncollision_ledger, &fit);
            if (refine_arena != NULL) Arena_dispose(&refine_arena);
            if (refine_rc != 0) {
                fprintf(stderr, "%s: collision U(v) refinement failed\n",
                        argv[0]);
                free(collision_ledger);
                Arena_dispose(&arena);
                return 1;
            }
        }
    }
    Arena_T texture_arena = NULL;
    AtlasRibbonTextureResult texture_phase;
    int have_texture_phase = 0;
    ArftTextureAlignStats texture_align_stats;
    ArftTextureAlignRound *texture_align_round = NULL;
    int have_texture_alignment = 0;
    memset(&texture_align_stats, 0, sizeof texture_align_stats);
    memset(&texture_phase, 0, sizeof texture_phase);
    if (raw_dir != NULL && opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON) {
        texture_arena = Arena_new();
        if (texture_arena == NULL ||
            AtlasRibbonTexture_measure(
                texture_arena, &ps, &solution, &observations, &fit,
                raw_dir, &texture_opts, &texture_phase) != 0) {
            fprintf(stderr, "%s: RAW texture-phase audit failed\n", argv[0]);
            if (texture_arena != NULL) Arena_dispose(&texture_arena);
            if (uv_audit_arena != NULL) Arena_dispose(&uv_audit_arena);
            free(collision_ledger);
            Arena_dispose(&arena);
            return 1;
        }
        have_texture_phase = 1;
        fprintf(stderr,
                "[atlas_ribbon_fit] RAW phase pairs=%zu accepted=%zu "
                "corr=%.3f margin=%.3f delta_rms=%.3f desired=%.3f "
                "residual=%.3f p95=%.3f max=%.3f samples=%zu/%zu\n",
                texture_phase.npairs, texture_phase.accepted_pairs,
                texture_phase.accepted_correlation_mean,
                texture_phase.accepted_margin_mean,
                texture_phase.current_delta_rms,
                texture_phase.desired_delta_rms,
                texture_phase.residual_rms, texture_phase.residual_p95,
                texture_phase.residual_max,
                texture_phase.sampled_layer_points,
                texture_phase.sampled_layer_points + texture_phase.missing_layer_points);
        if (texture_rounds > 0) {
            if (run_texture_alignment(
                    texture_arena, &ps, &solution, &observations, &opts,
                    &texture_phase, texture_rounds, texture_sweeps,
                    texture_weight, texture_anchor_weight,
                    &collision_ledger, &ncollision_ledger,
                    &collision_ledger_capacity, &uv_audit,
                    &uv_audit_arena, &have_uv_audit, &fit,
                    &texture_align_stats, &texture_align_round) != 0) {
                fprintf(stderr, "%s: texture U alignment failed\n", argv[0]);
                if (texture_arena != NULL) Arena_dispose(&texture_arena);
                if (uv_audit_arena != NULL)
                    Arena_dispose(&uv_audit_arena);
                free(collision_ledger);
                Arena_dispose(&arena);
                return 1;
            }
            have_texture_alignment = 1;
        }
    }
    int io = 0;
    if (opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON &&
        AtlasRibbonFit_build_registered_ribbon(
            arena, &ps, &solution, &observations, &opts, &fit) != 0) {
        fprintf(stderr, "%s: registered ribbon construction failed\n", argv[0]);
        if (texture_arena != NULL) Arena_dispose(&texture_arena);
        if (uv_audit_arena != NULL) Arena_dispose(&uv_audit_arena);
        free(collision_ledger);
        Arena_dispose(&arena);
        return 1;
    }

    io |= write_observations_obj(out_dir, &observations);
    io |= write_targets_csv(out_dir, &observations);
    io |= write_conflicts_csv(out_dir, &observations);
    io |= write_bridge_cuts_csv(out_dir, &observations, &opts);
    size_t ribbon_faces = 0;
    if (opts.mode != ATLAS_RIBBON_OBSERVATIONS) {
        io |= write_registered_samples_csv(out_dir, &observations, &fit);
        io |= write_register_bounds_csv(out_dir, &fit);
    }
    if (have_uv_audit)
        io |= write_uv_audit(out_dir, &uv_audit);
    if (ncollision_ledger > 0)
        io |= write_collision_ledger(out_dir, collision_ledger, ncollision_ledger);
    if (have_texture_phase)
        io |= write_texture_phase(out_dir, &texture_opts, &texture_phase);
    if (have_texture_alignment)
        io |= write_texture_alignment(
            out_dir, &texture_align_stats, texture_align_round);
    if (opts.mode == ATLAS_RIBBON_POSITION ||
        opts.mode == ATLAS_RIBBON_TANGENT ||
        opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON) {
        io |= write_ribbon_obj(
            arena, out_dir, &observations, &fit, &ribbon_faces);
        io |= write_ribbon_rows_obj(out_dir, &observations, &fit);
        io |= write_rows_csv(out_dir, &fit);
        if (opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON)
            io |= write_coverage_audit(
                arena, out_dir, &ps, &solution, &observations, &opts, &fit);
    }
    double seconds = ves_clock_sec() - start;
    io |= write_stats_json(out_dir, &ps, &solution, &opts,
                           &observations, &fit, ribbon_faces, seconds);
    if (io != 0) {
        fprintf(stderr, "%s: output write failed\n", argv[0]);
        if (texture_arena != NULL) Arena_dispose(&texture_arena);
        Arena_dispose(&arena);
        if (uv_audit_arena != NULL) Arena_dispose(&uv_audit_arena);
        free(collision_ledger);
        return 1;
    }
    if (opts.mode == ATLAS_RIBBON_REGISTER_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTER_FIELD_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTER_FIELD_SMOOTH_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTER_COLLISION_ONLY ||
        opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON) {
        fprintf(stderr,
                "[atlas_ribbon_fit] %s dir=%+d iterations=%zu "
                "bounds=%zu inversions=%zu->%zu unresolved=%zu "
                "width=%.4g->%.4g (%+.2f%%) shift_rms=%.4g max=%.4g "
                "max_arc=%.4g local_gap=%.4g/%.4g vertical=%.4g/%.4g "
                "row_width=%.4g qp_fail=%zu collision=%zu/%zu "
                "seconds=%.3f\n",
                mode_name(opts.mode), fit.register_winding_direction,
                fit.register_iterations, fit.register_constraints,
                fit.register_inversions_before, fit.register_inversions_after,
                fit.register_unresolved, fit.register_width_before,
                fit.register_width_after,
                100.0 * (fit.register_width_after /
                         fit.register_width_before - 1.0),
                fit.register_shift_rms, fit.register_shift_max,
                fit.register_max_required_arc,
                fit.register_local_gap_mean, fit.register_local_gap_max,
                fit.register_vertical_shift_rms,
                fit.register_vertical_shift_max,
                fit.register_row_width_ratio_mean,
                fit.register_qp_failures,
                fit.register_collision_rounds,
                fit.register_collision_constraints, seconds);
        if (have_uv_audit)
            fprintf(stderr,
                    "[atlas_ribbon_fit] exact UV checkpoint=%zu registered=%zu "
                    "persisted=%zu new=%zu resolved=%zu known=%zu allowed=%zu "
                    "hard=%zu unclassified=%zu collision_bounds=%zu/%zu "
                    "no_row=%zu interval_fail=%zu\n",
                    uv_audit.initial.npairs, uv_audit.registered.npairs,
                    uv_audit.persisted_pairs, uv_audit.new_pairs,
                    uv_audit.resolved_pairs, uv_audit.registered_known_pairs,
                    uv_audit.registered_allowed_pairs,
                    uv_audit.registered_hard_pairs,
                    uv_audit.registered_unclassified_pairs,
                    uv_audit.ncollision_bounds,
                    uv_audit.collision_bound_face_pairs,
                    uv_audit.collision_pairs_without_row,
                    uv_audit.collision_interval_failures);
    } else if (opts.mode != ATLAS_RIBBON_OBSERVATIONS) {
        fprintf(stderr,
                "[atlas_ribbon_fit] nodes=%zu faces=%zu crossings=%zu "
                "locked_rms=%.4g chamfer_mean=%.4g p95=%.4g "
                "hausdorff=%.4g backsteps=%zu radius_order=%zu/%zu "
                "seconds=%.3f\n",
                fit.fitted_nodes, ribbon_faces, fit.row_crossings,
                fit.parameter_rms, fit.chamfer_mean, fit.chamfer_p95,
                fit.directed_hausdorff, fit.winding_phase_backsteps,
                fit.winding_radius_order_violations,
                fit.winding_radius_order_tests, seconds);
    }
    if (opts.mode == ATLAS_RIBBON_REGISTERED_RIBBON)
        fprintf(stderr,
                "[atlas_ribbon_fit] ribbon nodes=%zu faces=%zu direct=%zu "
                "spiral=%zu cuts=%zu/%zu/%zu crossings=%zu crossing_cuts=%zu "
                "locked_rms=%.4g chamfer_mean=%.4g p95=%.4g "
                "hausdorff=%.4g backsteps=%zu radius_order=%zu/%zu\n",
                fit.fitted_nodes, ribbon_faces, fit.direct_spans,
                fit.spiral_spans, fit.u_gap_cuts, fit.metric_jump_cuts,
                fit.topology_cuts, fit.row_crossings, fit.crossing_cuts,
                fit.parameter_rms,
                fit.chamfer_mean, fit.chamfer_p95, fit.directed_hausdorff,
                fit.winding_phase_backsteps,
                fit.winding_radius_order_violations,
                fit.winding_radius_order_tests);
    if (texture_arena != NULL) Arena_dispose(&texture_arena);
    Arena_dispose(&arena);
    if (uv_audit_arena != NULL) Arena_dispose(&uv_audit_arena);
    free(collision_ledger);
    return 0;
}
