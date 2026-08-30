#include "atlas_track_grow.h"

#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t a, b, c, d;
    float fit_rms;
    float fit_max;
    uint8_t inferred;
    uint8_t observed_edges;
} AtgQuad;

enum {
    ATG_INFER_NONE = 0,
    ATG_INFER_FRONTIER = 1,
    ATG_INFER_UV_BANK = 2
};

enum {
    ATG_GAP_EDGE_NONE = 0,
    ATG_GAP_EDGE_TOP = 1,
    ATG_GAP_EDGE_BOTTOM = 2,
    ATG_GAP_EDGE_LEFT = 3,
    ATG_GAP_EDGE_RIGHT = 4,
    ATG_GAP_EDGE_HORIZONTAL_PAIR = 5,
    ATG_GAP_EDGE_VERTICAL_PAIR = 6
};

typedef struct {
    AtgQuad quad;
    double score;
    uint8_t missing_edge;
} AtgGapProposal;

enum {
    ATG_BRIDGE_HORIZONTAL = 0,
    ATG_BRIDGE_VERTICAL = 1,
    ATG_BRIDGE_POINT = 2
};

typedef struct {
    int32_t a;
    int32_t b;
    double delta;       /* desired global_u[b] - global_u[a] */
    uint8_t kind;
} AtgParameterBridge;

typedef struct {
    double cell;
    double p0_min;
    double p1_min;
    int32_t nx;
    int32_t ny;
    int32_t *head;
    int32_t *next;
} AtgSpatialGrid;

static int atg_intervals_overlap(double a0, double a1,
                                 double b0, double b1)
{
    return a0 < b1 && b0 < a1;
}

static double atg_dot2(const double a[2], const double b[2])
{
    return a[0] * b[0] + a[1] * b[1];
}

static double atg_distance2(const double a[2], const double b[2])
{
    double x = b[0] - a[0];
    double y = b[1] - a[1];
    return sqrt(x * x + y * y);
}

static double atg_dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void atg_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static void atg_world_point(
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonTarget *target,
    double xyz[3])
{
    for (int d = 0; d < 3; d++)
        xyz[d] = set->axis_point[d] + target->v * set->axis[d] +
                 target->p[0] * set->basis0[d] +
                 target->p[1] * set->basis1[d];
}

static int atg_direct_relations_valid(const AtlasTrackGrowOptions *o)
{
    if (o == NULL) return 0;
    if (o->direct_relation_count == 0)
        return o->direct_relation_a == NULL && o->direct_relation_b == NULL;
    if (o->direct_relation_a == NULL || o->direct_relation_b == NULL)
        return 0;
    for (size_t i = 0; i < o->direct_relation_count; i++) {
        int32_t a = o->direct_relation_a[i];
        int32_t b = o->direct_relation_b[i];
        if (a < 0 || b <= a) return 0;
        if (i > 0) {
            int32_t pa = o->direct_relation_a[i - 1];
            int32_t pb = o->direct_relation_b[i - 1];
            if (a < pa || (a == pa && b <= pb)) return 0;
        }
    }
    return 1;
}

static int atg_options_valid(const AtlasTrackGrowOptions *o)
{
    return o != NULL &&
           atg_direct_relations_valid(o) &&
           isfinite(o->horizontal_max_gap) &&
           o->horizontal_max_gap > 0.0 &&
           isfinite(o->horizontal_stretch) &&
           o->horizontal_stretch >= 1.0 &&
           isfinite(o->horizontal_slack) &&
           o->horizontal_slack >= 0.0 &&
           isfinite(o->horizontal_min_fraction) &&
           o->horizontal_min_fraction >= 0.0 &&
           o->horizontal_min_fraction <= 1.0 &&
           isfinite(o->vertical_search_u) &&
           o->vertical_search_u >= 0.0 &&
           isfinite(o->vertical_radius) &&
           o->vertical_radius >= 0.0 &&
           isfinite(o->vertical_chain_residual_max) &&
           o->vertical_chain_residual_max >= 0.0 &&
           o->vertical_winding_delta_max >= 0 &&
           isfinite(o->tangent_dot_min) &&
           o->tangent_dot_min >= -1.0 &&
           o->tangent_dot_min <= 1.0 &&
           isfinite(o->quad_normal_dot_min) &&
           o->quad_normal_dot_min >= -1.0 &&
           o->quad_normal_dot_min <= 1.0 &&
           o->min_component_quads >= 1 &&
           o->min_component_quads <= INT32_MAX &&
           o->gap_fill_rounds >= 0 &&
           o->gap_fill_rounds <= 16 &&
           isfinite(o->gap_fit_rms_max) &&
           o->gap_fit_rms_max > 0.0 &&
           isfinite(o->gap_fit_max_max) &&
           o->gap_fit_max_max >= o->gap_fit_rms_max &&
           isfinite(o->uv_bridge_max_u) &&
           o->uv_bridge_max_u >= 0.0 &&
           isfinite(o->uv_bridge_fit_rms_max) &&
           o->uv_bridge_fit_rms_max > 0.0 &&
           isfinite(o->uv_bridge_fit_max_max) &&
           o->uv_bridge_fit_max_max >= o->uv_bridge_fit_rms_max &&
           isfinite(o->parameter_stay_weight) &&
           o->parameter_stay_weight >= 0.0 &&
           isfinite(o->parameter_bridge_weight) &&
           o->parameter_bridge_weight > 0.0 &&
           o->parameter_iterations >= 0 &&
           o->parameter_iterations <= 10000;
}

void AtlasTrackGrowOptions_default(AtlasTrackGrowOptions *opts)
{
    if (opts == NULL) return;
    memset(opts, 0, sizeof *opts);
    opts->horizontal_max_gap = 24.0;
    opts->horizontal_stretch = 1.6;
    opts->horizontal_slack = 2.5;
    opts->horizontal_min_fraction = 0.35;
    opts->vertical_search_u = 12.0;
    opts->vertical_radius = 8.0;
    opts->tangent_dot_min = 0.35;
    opts->physical_neighbours = 0;
    opts->preserve_layout_parameterization = 0;
    opts->vertical_winding_delta_max = INT_MAX;
    opts->vertical_chain_residual_max = 0.0;
    opts->quad_normal_dot_min = 0.0;
    opts->min_component_quads = 64;
    opts->gap_fill_rounds = 2;
    opts->gap_fit_rms_max = 1.5;
    opts->gap_fit_max_max = 4.0;
    opts->uv_bridge_max_u = 0.0;
    opts->uv_bridge_fit_rms_max = 0.75;
    opts->uv_bridge_fit_max_max = 2.0;
    opts->parameter_stay_weight = 0.01;
    opts->parameter_bridge_weight = 0.25;
    opts->parameter_iterations = 160;
}

static int atg_targets_sorted(const AtlasRibbonObservationSet *set)
{
    for (size_t i = 1; i < set->ntarget; i++) {
        const AtlasRibbonTarget *a = &set->target[i - 1];
        const AtlasRibbonTarget *b = &set->target[i];
        if (a->row > b->row ||
            (a->row == b->row && a->column > b->column))
            return 0;
    }
    return 1;
}

static int atg_key_compare(int32_t row0, int32_t column0,
                           int32_t row1, int32_t column1)
{
    if (row0 != row1) return row0 < row1 ? -1 : 1;
    return column0 < column1 ? -1 : (column0 > column1 ? 1 : 0);
}

/* Targets are the geometry consumed by the grower; layer samples carry the
 * authoritative chart rank and per-sample snapped winding.  Join the two
 * sorted streams once so every published vertex keeps that identity. */
static int atg_map_target_layers(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    int32_t **out_layer)
{
    int32_t *map = (int32_t *)ARENA_ALLOC(
        arena, set->ntarget * sizeof(*map));
    for (size_t i = 0; i < set->ntarget; i++) map[i] = -1;
    size_t layer = 0;
    for (size_t target = 0; target < set->ntarget; target++) {
        const AtlasRibbonTarget *t = &set->target[target];
        while (layer < set->nlayer_samples) {
            const AtlasRibbonLayerSample *s = &set->layer_sample[layer];
            int comparison = atg_key_compare(
                s->row, s->source_column, t->row, t->column);
            if (comparison >= 0) break;
            layer++;
        }
        size_t end = layer;
        while (end < set->nlayer_samples) {
            const AtlasRibbonLayerSample *s = &set->layer_sample[end];
            if (atg_key_compare(
                    s->row, s->source_column, t->row, t->column) != 0)
                break;
            if (map[target] < 0 &&
                (s->chart == t->chart0 || t->chart0 < 0))
                map[target] = (int32_t)end;
            end++;
        }
        if (t->accepted && map[target] < 0) return -1;
    }
    *out_layer = map;
    return 0;
}

static int atg_winding_direction(const AtlasRibbonObservationSet *set)
{
    size_t positive = 0, negative = 0;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        double radius = hypot(s->p[0], s->p[1]);
        if (!(radius > 1.0e-6)) continue;
        double cross = s->p[0] * s->tangent[1] -
                       s->p[1] * s->tangent[0];
        if (cross > 0.05 * radius) positive++;
        else if (cross < -0.05 * radius) negative++;
    }
    return positive >= negative ? 1 : -1;
}

static int atg_horizontal_compatible_max(
    const AtlasRibbonTarget *a,
    const AtlasRibbonTarget *b,
    double gap,
    double max_gap,
    const AtlasTrackGrowOptions *opts,
    double *out_score)
{
    if (!(gap > 0.0) || gap > max_gap)
        return 0;
    double delta[2] = {
        b->p[0] - a->p[0],
        b->p[1] - a->p[1]
    };
    double distance = sqrt(atg_dot2(delta, delta));
    if (!(distance > 1.0e-9) ||
        distance > opts->horizontal_stretch * gap +
                   opts->horizontal_slack ||
        distance + opts->horizontal_slack <
            opts->horizontal_min_fraction * gap)
        return 0;
    delta[0] /= distance;
    delta[1] /= distance;
    double along_a = atg_dot2(delta, a->tangent);
    double along_b = atg_dot2(delta, b->tangent);
    double tangent_dot = atg_dot2(a->tangent, b->tangent);
    if (along_a < opts->tangent_dot_min ||
        along_b < opts->tangent_dot_min ||
        tangent_dot < opts->tangent_dot_min)
        return 0;
    if (out_score != NULL) {
        double metric = fabs(distance - gap) /
                        (gap + opts->horizontal_slack + 1.0e-12);
        *out_score = metric +
                     0.25 * (2.0 - along_a - along_b) +
                     0.10 * (1.0 - tangent_dot) +
                     1.0e-4 * gap;
    }
    return 1;
}

static int atg_horizontal_compatible(
    const AtlasRibbonTarget *a,
    const AtlasRibbonTarget *b,
    double gap,
    const AtlasTrackGrowOptions *opts,
    double *out_score)
{
    return atg_horizontal_compatible_max(
        a, b, gap, opts->horizontal_max_gap, opts, out_score);
}

static int atg_vertical_compatible(
    const AtlasRibbonTarget *a,
    const AtlasRibbonTarget *b,
    double du,
    const AtlasTrackGrowOptions *opts,
    double *out_score)
{
    if (du > opts->vertical_search_u) return 0;
    double distance = atg_distance2(a->p, b->p);
    double tangent_dot = atg_dot2(a->tangent, b->tangent);
    if (distance > opts->vertical_radius ||
        tangent_dot < opts->tangent_dot_min)
        return 0;
    if (out_score != NULL) {
        double chart_factor =
            a->chart0 >= 0 && a->chart0 == b->chart0 ? 0.85 : 1.0;
        *out_score = chart_factor *
                     (distance + 0.20 * du +
                      0.50 * (1.0 - tangent_dot));
    }
    return 1;
}

static int atg_same_layer_group(
    const AtlasRibbonObservationSet *set,
    const int32_t *target_layer,
    size_t a,
    size_t b)
{
    int32_t la = target_layer[a];
    int32_t lb = target_layer[b];
    if (la < 0 || lb < 0 ||
        (size_t)la >= set->nlayer_samples ||
        (size_t)lb >= set->nlayer_samples)
        return 0;
    return set->layer_sample[la].group ==
           set->layer_sample[lb].group;
}

static int atg_direct_relation_contains(
    const AtlasTrackGrowOptions *opts,
    int32_t a,
    int32_t b)
{
    if (a == b) return a >= 0;
    if (a < 0 || b < 0 || opts == NULL ||
        opts->direct_relation_count == 0)
        return 0;
    if (b < a) {
        int32_t swap = a;
        a = b;
        b = swap;
    }
    size_t lo = 0, hi = opts->direct_relation_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int32_t ea = opts->direct_relation_a[mid];
        int32_t eb = opts->direct_relation_b[mid];
        if (ea < a || (ea == a && eb < b)) lo = mid + 1;
        else hi = mid;
    }
    return lo < opts->direct_relation_count &&
           opts->direct_relation_a[lo] == a &&
           opts->direct_relation_b[lo] == b;
}

static int atg_chart_join_allowed(
    const AtlasRibbonObservationSet *set,
    const int32_t *target_layer,
    size_t a,
    size_t b,
    const AtlasTrackGrowOptions *opts,
    int legacy_group_gate)
{
    int32_t la = target_layer[a];
    int32_t lb = target_layer[b];
    if (la < 0 || lb < 0 ||
        (size_t)la >= set->nlayer_samples ||
        (size_t)lb >= set->nlayer_samples)
        return 0;
    int32_t ca = set->layer_sample[la].chart;
    int32_t cb = set->layer_sample[lb].chart;
    if (ca < 0 || cb < 0) return 0;
    if (ca == cb) return 1;
    if (opts->direct_relation_count > 0)
        return atg_direct_relation_contains(opts, ca, cb);
    return !legacy_group_gate || atg_same_layer_group(
        set, target_layer, a, b);
}

static int atg_vertical_winding_compatible(
    const AtlasRibbonObservationSet *set,
    const int32_t *target_layer,
    size_t a,
    size_t b,
    int maximum_delta)
{
    int32_t la = target_layer[a];
    int32_t lb = target_layer[b];
    if (la < 0 || lb < 0 ||
        (size_t)la >= set->nlayer_samples ||
        (size_t)lb >= set->nlayer_samples)
        return 0;
    int64_t delta =
        (int64_t)set->layer_sample[la].winding -
        (int64_t)set->layer_sample[lb].winding;
    if (delta < 0) delta = -delta;
    return delta <= maximum_delta;
}

static int atg_physical_horizontal_compatible(
    const AtlasRibbonTarget *a,
    const AtlasRibbonTarget *b,
    double nominal_step,
    const AtlasTrackGrowOptions *opts,
    double *out_score)
{
    double delta[2] = {
        b->p[0] - a->p[0],
        b->p[1] - a->p[1]
    };
    double distance = sqrt(atg_dot2(delta, delta));
    double minimum =
        opts->horizontal_min_fraction * nominal_step;
    if (!(distance > 1.0e-9) || distance < minimum ||
        distance > opts->horizontal_max_gap)
        return 0;
    double tangent_dot = atg_dot2(a->tangent, b->tangent);
    if (tangent_dot < opts->tangent_dot_min) return 0;
    double along_a = atg_dot2(delta, a->tangent) / distance;
    double along_b = atg_dot2(delta, b->tangent) / distance;
    if (along_a < opts->tangent_dot_min ||
        along_b < opts->tangent_dot_min)
        return 0;
    double cross_a = fabs(delta[0] * a->tangent[1] -
                          delta[1] * a->tangent[0]);
    double cross_b = fabs(delta[0] * b->tangent[1] -
                          delta[1] * b->tangent[0]);
    if (cross_a > opts->horizontal_slack ||
        cross_b > opts->horizontal_slack)
        return 0;
    if (out_score != NULL) {
        double metric = fabs(distance - nominal_step) /
                        (nominal_step + opts->horizontal_slack + 1.0e-12);
        *out_score = metric +
                     0.50 * (cross_a + cross_b) +
                     0.25 * (2.0 - along_a - along_b) +
                     0.10 * (1.0 - tangent_dot) +
                     1.0e-4 * distance;
    }
    return 1;
}

static int atg_physical_vertical_compatible(
    const AtlasRibbonTarget *a,
    const AtlasRibbonTarget *b,
    const AtlasTrackGrowOptions *opts,
    double *out_score)
{
    double distance = atg_distance2(a->p, b->p);
    double tangent_dot = atg_dot2(a->tangent, b->tangent);
    if (distance > opts->vertical_radius ||
        tangent_dot < opts->tangent_dot_min)
        return 0;
    if (out_score != NULL) {
        double chart_factor =
            a->chart0 >= 0 && a->chart0 == b->chart0 ? 0.85 : 1.0;
        *out_score = chart_factor *
                     (distance + 0.50 * (1.0 - tangent_dot));
    }
    return 1;
}

static int atg_spatial_grid_init(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    double cell,
    AtgSpatialGrid *grid)
{
    if (arena == NULL || set == NULL || grid == NULL ||
        !(cell > 0.0) || !isfinite(cell))
        return -1;
    double p0_min = DBL_MAX, p0_max = -DBL_MAX;
    double p1_min = DBL_MAX, p1_max = -DBL_MAX;
    for (size_t i = 0; i < set->ntarget; i++) {
        const AtlasRibbonTarget *target = &set->target[i];
        if (!target->accepted) continue;
        if (!isfinite(target->p[0]) || !isfinite(target->p[1]))
            return -1;
        if (target->p[0] < p0_min) p0_min = target->p[0];
        if (target->p[0] > p0_max) p0_max = target->p[0];
        if (target->p[1] < p1_min) p1_min = target->p[1];
        if (target->p[1] > p1_max) p1_max = target->p[1];
    }
    if (p0_min == DBL_MAX || p1_min == DBL_MAX) return -1;
    double nx_real = floor((p0_max - p0_min) / cell) + 1.0;
    double ny_real = floor((p1_max - p1_min) / cell) + 1.0;
    if (!(nx_real >= 1.0) || !(ny_real >= 1.0) ||
        nx_real > INT32_MAX || ny_real > INT32_MAX)
        return -1;
    int32_t nx = (int32_t)nx_real;
    int32_t ny = (int32_t)ny_real;
    if ((size_t)nx > SIZE_MAX / (size_t)ny ||
        (size_t)nx * (size_t)ny > (size_t)INT32_MAX)
        return -1;
    grid->cell = cell;
    grid->p0_min = p0_min;
    grid->p1_min = p1_min;
    grid->nx = nx;
    grid->ny = ny;
    grid->head = (int32_t *)ARENA_ALLOC(
        arena, (size_t)nx * (size_t)ny * sizeof(*grid->head));
    grid->next = (int32_t *)ARENA_ALLOC(
        arena, set->ntarget * sizeof(*grid->next));
    return 0;
}

static int32_t atg_spatial_x(
    const AtgSpatialGrid *grid,
    double p0)
{
    int64_t value = (int64_t)floor((p0 - grid->p0_min) / grid->cell);
    if (value < 0) value = 0;
    if (value >= grid->nx) value = grid->nx - 1;
    return (int32_t)value;
}

static int32_t atg_spatial_y(
    const AtgSpatialGrid *grid,
    double p1)
{
    int64_t value = (int64_t)floor((p1 - grid->p1_min) / grid->cell);
    if (value < 0) value = 0;
    if (value >= grid->ny) value = grid->ny - 1;
    return (int32_t)value;
}

static void atg_spatial_grid_build(
    const AtlasRibbonObservationSet *set,
    size_t first,
    size_t last,
    AtgSpatialGrid *grid)
{
    size_t cells = (size_t)grid->nx * (size_t)grid->ny;
    for (size_t i = 0; i < cells; i++) grid->head[i] = -1;
    for (size_t i = first; i < last; i++) {
        if (!set->target[i].accepted) continue;
        int32_t x = atg_spatial_x(grid, set->target[i].p[0]);
        int32_t y = atg_spatial_y(grid, set->target[i].p[1]);
        size_t bucket = (size_t)y * (size_t)grid->nx + (size_t)x;
        grid->next[i] = grid->head[bucket];
        grid->head[bucket] = (int32_t)i;
    }
}

static size_t atg_lower_bound_column(
    const AtlasRibbonTarget *target,
    size_t first,
    size_t last,
    int32_t column)
{
    while (first < last) {
        size_t middle = first + (last - first) / 2;
        if (target[middle].column < column)
            first = middle + 1;
        else
            last = middle;
    }
    return first;
}

static int atg_quad_valid(
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonTarget *target,
    int32_t a,
    int32_t b,
    int32_t c,
    int32_t d,
    double normal_dot_min)
{
    double pa[3], pb[3], pc[3], pd[3];
    atg_world_point(set, &target[a], pa);
    atg_world_point(set, &target[b], pb);
    atg_world_point(set, &target[c], pc);
    atg_world_point(set, &target[d], pd);
    double ab[3], ac[3], dc[3], db[3];
    for (int k = 0; k < 3; k++) {
        ab[k] = pb[k] - pa[k];
        ac[k] = pc[k] - pa[k];
        dc[k] = pc[k] - pd[k];
        db[k] = pb[k] - pd[k];
    }
    double n0[3], n1[3];
    atg_cross3(ab, ac, n0);
    atg_cross3(dc, db, n1);
    double l0 = sqrt(atg_dot3(n0, n0));
    double l1 = sqrt(atg_dot3(n1, n1));
    if (!(l0 > 1.0e-9) || !(l1 > 1.0e-9)) return 0;
    return atg_dot3(n0, n1) >= normal_dot_min * l0 * l1;
}

static int atg_solve3(
    const double matrix[3][3],
    const double rhs[3],
    double solution[3])
{
    double a[3][4];
    for (int row = 0; row < 3; row++) {
        for (int column = 0; column < 3; column++)
            a[row][column] = matrix[row][column];
        a[row][3] = rhs[row];
    }
    for (int pivot = 0; pivot < 3; pivot++) {
        int best = pivot;
        for (int row = pivot + 1; row < 3; row++)
            if (fabs(a[row][pivot]) > fabs(a[best][pivot])) best = row;
        if (!(fabs(a[best][pivot]) > 1.0e-12)) return -1;
        if (best != pivot)
            for (int column = pivot; column < 4; column++) {
                double temporary = a[pivot][column];
                a[pivot][column] = a[best][column];
                a[best][column] = temporary;
            }
        double scale = a[pivot][pivot];
        for (int column = pivot; column < 4; column++)
            a[pivot][column] /= scale;
        for (int row = 0; row < 3; row++) {
            if (row == pivot) continue;
            double factor = a[row][pivot];
            for (int column = pivot; column < 4; column++)
                a[row][column] -= factor * a[pivot][column];
        }
    }
    for (int row = 0; row < 3; row++) solution[row] = a[row][3];
    return 0;
}

/* Best local affine map XYZ = c0 + c1*U + c2*V over observed UV samples.
 * Centering and scaling make the small normal equations insensitive to the
 * absolute winding coordinate (tens of thousands of voxels on real data). */
static int atg_uvxyz_fit(
    const AtlasRibbonObservationSet *set,
    const double *target_u,
    const int32_t *index,
    size_t count,
    double *out_rms,
    double *out_max)
{
    enum { ATG_FIT_MAX_POINTS = 16 };
    if (index == NULL || count < 4 || count > ATG_FIT_MAX_POINTS)
        return -1;
    double xyz[ATG_FIT_MAX_POINTS][3];
    double u[ATG_FIT_MAX_POINTS], v[ATG_FIT_MAX_POINTS];
    double umean = 0.0, vmean = 0.0;
    for (size_t i = 0; i < count; i++) {
        u[i] = target_u[index[i]];
        v[i] = set->target[index[i]].v;
        if (!isfinite(u[i]) || !isfinite(v[i])) return -1;
        atg_world_point(set, &set->target[index[i]], xyz[i]);
        umean += u[i];
        vmean += v[i];
    }
    umean /= (double)count;
    vmean /= (double)count;
    double uscale2 = 0.0, vscale2 = 0.0;
    for (size_t i = 0; i < count; i++) {
        uscale2 += (u[i] - umean) * (u[i] - umean);
        vscale2 += (v[i] - vmean) * (v[i] - vmean);
    }
    double uscale = sqrt(uscale2 / (double)count);
    double vscale = sqrt(vscale2 / (double)count);
    if (!(uscale > 1.0e-9) || !(vscale > 1.0e-9)) return -1;

    double design[ATG_FIT_MAX_POINTS][3];
    double normal[3][3] = {{0}};
    double rhs[3][3] = {{0}};
    for (size_t i = 0; i < count; i++) {
        design[i][0] = 1.0;
        design[i][1] = (u[i] - umean) / uscale;
        design[i][2] = (v[i] - vmean) / vscale;
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < 3; column++)
                normal[row][column] += design[i][row] * design[i][column];
            for (int axis = 0; axis < 3; axis++)
                rhs[axis][row] += design[i][row] * xyz[i][axis];
        }
    }
    double coefficient[3][3];
    for (int axis = 0; axis < 3; axis++)
        if (atg_solve3(normal, rhs[axis], coefficient[axis]) != 0)
            return -1;

    double sum2 = 0.0, maximum = 0.0;
    for (size_t i = 0; i < count; i++) {
        double residual2 = 0.0;
        for (int axis = 0; axis < 3; axis++) {
            double predicted = 0.0;
            for (int term = 0; term < 3; term++)
                predicted += coefficient[axis][term] * design[i][term];
            double residual = predicted - xyz[i][axis];
            residual2 += residual * residual;
        }
        double residual = sqrt(residual2);
        sum2 += residual2;
        if (residual > maximum) maximum = residual;
    }
    *out_rms = sqrt(sum2 / (double)count);
    *out_max = maximum;
    return 0;
}

static int atg_uvxyz_fit4(
    const AtlasRibbonObservationSet *set,
    const double *target_u,
    int32_t ia,
    int32_t ib,
    int32_t ic,
    int32_t id,
    double *out_rms,
    double *out_max)
{
    int32_t index[4] = {ia, ib, ic, id};
    return atg_uvxyz_fit(
        set, target_u, index, 4, out_rms, out_max);
}

enum {
    ATG_GAP_VALID = 0,
    ATG_GAP_REJECT_EVIDENCE = 1,
    ATG_GAP_REJECT_TOPOLOGY = 2,
    ATG_GAP_REJECT_FIT = 3
};

static int atg_make_gap_proposal(
    const AtlasRibbonObservationSet *set,
    const double *target_u,
    const AtlasTrackGrowOptions *opts,
    int32_t a,
    int32_t b,
    int32_t c,
    int32_t d,
    uint8_t missing_edge,
    uint8_t observed_edges,
    AtgGapProposal *proposal)
{
    if (a < 0 || b < 0 || c < 0 || d < 0 ||
        a == b || a == c || a == d || b == c || b == d || c == d)
        return ATG_GAP_REJECT_TOPOLOGY;
    const AtlasRibbonTarget *ta = &set->target[a];
    const AtlasRibbonTarget *tb = &set->target[b];
    const AtlasRibbonTarget *tc = &set->target[c];
    const AtlasRibbonTarget *td = &set->target[d];
    if (!ta->accepted || !tb->accepted || !tc->accepted || !td->accepted)
        return ATG_GAP_REJECT_EVIDENCE;
    if (ta->row != tb->row || tc->row != td->row ||
        tc->row != ta->row + 1 ||
        !(tb->column > ta->column) || !(td->column > tc->column))
        return ATG_GAP_REJECT_TOPOLOGY;

    double edge_score = 0.0;
    int compatible = 0;
    if (missing_edge == ATG_GAP_EDGE_TOP) {
        double gap = (double)(tb->column - ta->column) * set->observation_du;
        compatible = atg_horizontal_compatible(
            ta, tb, gap, opts, &edge_score);
    } else if (missing_edge == ATG_GAP_EDGE_BOTTOM) {
        double gap = (double)(td->column - tc->column) * set->observation_du;
        compatible = atg_horizontal_compatible(
            tc, td, gap, opts, &edge_score);
    } else if (missing_edge == ATG_GAP_EDGE_LEFT) {
        double du = fabs((double)(tc->column - ta->column) *
                         set->observation_du);
        compatible = atg_vertical_compatible(
            ta, tc, du, opts, &edge_score);
    } else if (missing_edge == ATG_GAP_EDGE_RIGHT) {
        double du = fabs((double)(td->column - tb->column) *
                         set->observation_du);
        compatible = atg_vertical_compatible(
            tb, td, du, opts, &edge_score);
    } else if (missing_edge == ATG_GAP_EDGE_HORIZONTAL_PAIR) {
        double top_gap =
            (double)(tb->column - ta->column) * set->observation_du;
        double bottom_gap =
            (double)(td->column - tc->column) * set->observation_du;
        double top_score = 0.0, bottom_score = 0.0;
        compatible = atg_horizontal_compatible(
                         ta, tb, top_gap, opts, &top_score) &&
                     atg_horizontal_compatible(
                         tc, td, bottom_gap, opts, &bottom_score);
        edge_score = 0.5 * (top_score + bottom_score);
    } else if (missing_edge == ATG_GAP_EDGE_VERTICAL_PAIR) {
        double left_du = fabs((double)(tc->column - ta->column) *
                              set->observation_du);
        double right_du = fabs((double)(td->column - tb->column) *
                               set->observation_du);
        double left_score = 0.0, right_score = 0.0;
        compatible = atg_vertical_compatible(
                         ta, tc, left_du, opts, &left_score) &&
                     atg_vertical_compatible(
                         tb, td, right_du, opts, &right_score);
        edge_score = 0.5 * (left_score + right_score);
    } else {
        compatible = 1;
    }
    if (!compatible) return ATG_GAP_REJECT_EVIDENCE;
    if (!atg_quad_valid(set, set->target, a, b, c, d,
                        opts->quad_normal_dot_min))
        return ATG_GAP_REJECT_TOPOLOGY;

    double fit_rms = 0.0, fit_max = 0.0;
    if (atg_uvxyz_fit4(
            set, target_u, a, b, c, d, &fit_rms, &fit_max) != 0 ||
        fit_rms > opts->gap_fit_rms_max ||
        fit_max > opts->gap_fit_max_max)
        return ATG_GAP_REJECT_FIT;
    memset(proposal, 0, sizeof *proposal);
    proposal->quad.a = a;
    proposal->quad.b = b;
    proposal->quad.c = c;
    proposal->quad.d = d;
    proposal->quad.fit_rms = (float)fit_rms;
    proposal->quad.fit_max = (float)fit_max;
    proposal->quad.inferred = 1;
    proposal->quad.observed_edges = observed_edges;
    proposal->missing_edge = missing_edge;
    proposal->score = fit_rms + 0.25 * fit_max + 0.10 * edge_score -
                      0.01 * (double)observed_edges;
    if (missing_edge == ATG_GAP_EDGE_HORIZONTAL_PAIR ||
        missing_edge == ATG_GAP_EDGE_VERTICAL_PAIR)
        proposal->score += 0.25;
    return ATG_GAP_VALID;
}

static int atg_gap_proposal_compare(const void *left, const void *right)
{
    const AtgGapProposal *a = (const AtgGapProposal *)left;
    const AtgGapProposal *b = (const AtgGapProposal *)right;
    if (a->score != b->score) return a->score < b->score ? -1 : 1;
    if (a->quad.a != b->quad.a) return a->quad.a < b->quad.a ? -1 : 1;
    if (a->quad.b != b->quad.b) return a->quad.b < b->quad.b ? -1 : 1;
    if (a->quad.c != b->quad.c) return a->quad.c < b->quad.c ? -1 : 1;
    if (a->quad.d != b->quad.d) return a->quad.d < b->quad.d ? -1 : 1;
    return 0;
}

static void atg_consider_gap_proposal(
    const AtlasRibbonObservationSet *set,
    const double *target_u,
    const AtlasTrackGrowOptions *opts,
    int32_t a,
    int32_t b,
    int32_t c,
    int32_t d,
    uint8_t missing_edge,
    uint8_t observed_edges,
    AtgGapProposal *best,
    int *have_best,
    AtlasTrackGrowResult *out)
{
    AtgGapProposal candidate;
    out->gap_fill_candidates++;
    int status = atg_make_gap_proposal(
        set, target_u, opts, a, b, c, d, missing_edge,
        observed_edges, &candidate);
    if (status == ATG_GAP_REJECT_EVIDENCE) {
        out->gap_fill_evidence_rejects++;
        return;
    }
    if (status == ATG_GAP_REJECT_TOPOLOGY) {
        out->gap_fill_topology_rejects++;
        return;
    }
    if (status == ATG_GAP_REJECT_FIT) {
        out->gap_fill_fit_rejects++;
        return;
    }
    if (!*have_best || atg_gap_proposal_compare(&candidate, best) < 0) {
        *best = candidate;
        *have_best = 1;
    }
}

static int atg_set_edge(
    int32_t *forward,
    int32_t *backward,
    int32_t a,
    int32_t b);

static int atg_apply_gap_proposal(
    const AtgGapProposal *proposal,
    int32_t *right,
    int32_t *left,
    int32_t *down,
    int32_t *up,
    uint8_t *right_inferred,
    uint8_t *down_inferred,
    AtlasTrackGrowResult *out)
{
    int32_t a = proposal->quad.a;
    int32_t b = proposal->quad.b;
    int32_t c = proposal->quad.c;
    int32_t d = proposal->quad.d;
    switch (proposal->missing_edge) {
        case ATG_GAP_EDGE_NONE:
            return right[a] == b && right[c] == d &&
                   down[a] == c && down[b] == d;
        case ATG_GAP_EDGE_TOP:
            if (down[a] != c || right[c] != d || down[b] != d ||
                right[a] >= 0 || left[b] >= 0)
                return 0;
            if (atg_set_edge(right, left, a, b) != 0) return 0;
            right_inferred[a] = 1;
            out->gap_fill_horizontal_edges++;
            return 1;
        case ATG_GAP_EDGE_BOTTOM:
            if (right[a] != b || down[a] != c || down[b] != d ||
                right[c] >= 0 || left[d] >= 0)
                return 0;
            if (atg_set_edge(right, left, c, d) != 0) return 0;
            right_inferred[c] = 1;
            out->gap_fill_horizontal_edges++;
            return 1;
        case ATG_GAP_EDGE_LEFT:
            if (right[a] != b || right[c] != d || down[b] != d ||
                down[a] >= 0 || up[c] >= 0)
                return 0;
            if (atg_set_edge(down, up, a, c) != 0) return 0;
            down_inferred[a] = 1;
            out->gap_fill_vertical_edges++;
            return 1;
        case ATG_GAP_EDGE_RIGHT:
            if (right[a] != b || down[a] != c || right[c] != d ||
                down[b] >= 0 || up[d] >= 0)
                return 0;
            if (atg_set_edge(down, up, b, d) != 0) return 0;
            down_inferred[b] = 1;
            out->gap_fill_vertical_edges++;
            return 1;
        case ATG_GAP_EDGE_HORIZONTAL_PAIR:
            if (down[a] != c || down[b] != d ||
                right[a] >= 0 || left[b] >= 0 ||
                right[c] >= 0 || left[d] >= 0)
                return 0;
            if (atg_set_edge(right, left, a, b) != 0 ||
                atg_set_edge(right, left, c, d) != 0)
                return 0;
            right_inferred[a] = right_inferred[c] = 1;
            out->gap_fill_horizontal_edges += 2;
            out->gap_fill_pair_quads++;
            return 1;
        case ATG_GAP_EDGE_VERTICAL_PAIR:
            if (right[a] != b || right[c] != d ||
                down[a] >= 0 || up[c] >= 0 ||
                down[b] >= 0 || up[d] >= 0)
                return 0;
            if (atg_set_edge(down, up, a, c) != 0 ||
                atg_set_edge(down, up, b, d) != 0)
                return 0;
            down_inferred[a] = down_inferred[b] = 1;
            out->gap_fill_vertical_edges += 2;
            out->gap_fill_pair_quads++;
            return 1;
        default:
            return 0;
    }
}

static void atg_claim_edge(
    int32_t *owner,
    int32_t edge,
    int32_t quad,
    UnionFind *uf)
{
    if (owner[edge] >= 0)
        uf_union(uf, owner[edge], quad);
    else
        owner[edge] = quad;
}

static double atg_edge_distance(
    const AtlasRibbonTarget *target,
    int32_t a,
    int32_t b)
{
    return atg_distance2(target[a].p, target[b].p);
}

static double atg_lattice_edge_span(
    const AtlasRibbonObservationSet *set,
    int32_t a,
    int32_t b)
{
    double steps = nearbyint(
        atg_edge_distance(set->target, a, b) / set->observation_du);
    if (steps < 1.0) steps = 1.0;
    return steps * set->observation_du;
}

static int atg_set_edge(
    int32_t *forward,
    int32_t *backward,
    int32_t a,
    int32_t b)
{
    if ((forward[a] >= 0 && forward[a] != b) ||
        (backward[b] >= 0 && backward[b] != a))
        return -1;
    forward[a] = b;
    backward[b] = a;
    return 0;
}

static int atg_uv_bank_winding_compatible(
    const AtlasRibbonObservationSet *set,
    const int32_t *target_layer,
    int32_t a,
    int32_t b,
    int32_t c,
    int32_t d)
{
    int32_t layer[4] = {
        target_layer[a], target_layer[b], target_layer[c], target_layer[d]
    };
    for (int i = 0; i < 4; i++)
        if (layer[i] < 0 || (size_t)layer[i] >= set->nlayer_samples)
            return 0;
    const AtlasRibbonLayerSample *sa = &set->layer_sample[layer[0]];
    const AtlasRibbonLayerSample *sb = &set->layer_sample[layer[1]];
    const AtlasRibbonLayerSample *sc = &set->layer_sample[layer[2]];
    const AtlasRibbonLayerSample *sd = &set->layer_sample[layer[3]];
    int64_t top_delta = (int64_t)sb->winding - (int64_t)sa->winding;
    int64_t bottom_delta = (int64_t)sd->winding - (int64_t)sc->winding;
    return sa->winding == sc->winding &&
           sb->winding == sd->winding &&
           top_delta == bottom_delta &&
           top_delta >= -1 && top_delta <= 1;
}

static int atg_row_has_face_node_between(
    const AtlasRibbonObservationSet *set,
    const size_t *row_start,
    const uint8_t *face_node,
    int32_t row,
    int32_t column_lo,
    int32_t column_hi)
{
    if (row < 0 || (size_t)row >= set->nrows ||
        column_hi <= column_lo + 1)
        return 0;
    size_t first = row_start[row];
    size_t last = row_start[(size_t)row + 1];
    size_t at = atg_lower_bound_column(
        set->target, first, last, column_lo + 1);
    for (; at < last && set->target[at].column < column_hi; at++)
        if (face_node[at]) return 1;
    return 0;
}

/* Bridge a horizontal ribbon hole from two complete core bank cells.  The
 * outer cells contribute eight UV->XYZ observations; neither the four bridge
 * corners alone nor a free extrapolation can admit a face.  Processing rows
 * in order lets one accepted edge become the top of the next row's bridge,
 * while every row still needs fresh original core support on both banks. */
static int atg_add_uv_bank_bridges(
    const AtlasRibbonObservationSet *set,
    const double *target_u,
    const int32_t *target_layer,
    const AtlasTrackGrowOptions *opts,
    const size_t *row_start,
    int32_t *right,
    int32_t *left,
    int32_t *down,
    int32_t *up,
    uint8_t *right_inferred,
    uint8_t *down_inferred,
    AtgQuad *quad,
    size_t quad_capacity,
    size_t *io_nquad,
    int32_t *quad_at_a,
    uint8_t *blocked_a,
    uint8_t *horizontal_face,
    uint8_t *vertical_face,
    AtlasTrackGrowResult *out)
{
    if (!(opts->uv_bridge_max_u > 0.0)) return 0;
    size_t n = set->ntarget;
    uint8_t *face_node = (uint8_t *)calloc(n, 1);
    if (face_node == NULL) return -1;
    for (size_t q = 0; q < *io_nquad; q++) {
        face_node[quad[q].a] = face_node[quad[q].b] = 1;
        face_node[quad[q].c] = face_node[quad[q].d] = 1;
    }

    double fit_sum2 = 0.0, fit_maximum = 0.0;
    size_t fit_count = 0;
    for (size_t ia = 0; ia < n; ia++) {
        int32_t a = (int32_t)ia;
        if (!set->target[a].accepted || blocked_a[a] ||
            quad_at_a[a] >= 0)
            continue;
        int32_t c = down[a];
        int32_t la = left[a];
        if (c < 0 || la < 0 || down_inferred[a]) continue;
        int32_t lc = down[la];
        int32_t qleft = quad_at_a[la];
        if (lc < 0 || qleft < 0 || quad[qleft].inferred != ATG_INFER_NONE ||
            right[la] != a || right[lc] != c || left[c] != lc ||
            down[la] != lc)
            continue;

        int32_t b = right[a];
        if (b < 0) {
            size_t last = row_start[(size_t)set->target[a].row + 1];
            for (size_t ib = ia + 1; ib < last; ib++) {
                double gap = target_u[ib] - target_u[a];
                if (gap > opts->uv_bridge_max_u) break;
                if (!(gap > 0.0) || !set->target[ib].accepted ||
                    left[ib] >= 0)
                    continue;
                int32_t candidate = (int32_t)ib;
                int32_t candidate_down = down[candidate];
                int32_t candidate_right = right[candidate];
                if (candidate_down < 0 || candidate_right < 0 ||
                    down_inferred[candidate] ||
                    quad_at_a[candidate] < 0 ||
                    quad[quad_at_a[candidate]].inferred != ATG_INFER_NONE)
                    continue;
                b = candidate;
                break;
            }
        }
        if (b < 0 || set->target[b].row != set->target[a].row)
            continue;
        int32_t d = down[b];
        int32_t rb = right[b];
        if (d < 0 || rb < 0 || down_inferred[b]) continue;
        int32_t rd = down[rb];
        int32_t qright = quad_at_a[b];
        if (rd < 0 || qright < 0 ||
            quad[qright].inferred != ATG_INFER_NONE ||
            right[d] != rd || up[d] != b || up[rd] != rb ||
            set->target[c].row != set->target[d].row)
            continue;

        int top_open = right[a] < 0 && left[b] < 0;
        int top_same = right[a] == b && left[b] == a;
        int bottom_open = right[c] < 0 && left[d] < 0;
        int bottom_same = right[c] == d && left[d] == c;
        if ((!top_open && !top_same) || (!bottom_open && !bottom_same) ||
            (!top_open && !bottom_open))
            continue;

        out->uv_bridge_candidates++;
        double top_gap = target_u[b] - target_u[a];
        double bottom_gap = target_u[d] - target_u[c];
        double top_score = 0.0, bottom_score = 0.0;
        if (!atg_uv_bank_winding_compatible(
                set, target_layer, a, b, c, d) ||
            !atg_horizontal_compatible_max(
                &set->target[a], &set->target[b], top_gap,
                opts->uv_bridge_max_u, opts, &top_score) ||
            !atg_horizontal_compatible_max(
                &set->target[c], &set->target[d], bottom_gap,
                opts->uv_bridge_max_u, opts, &bottom_score)) {
            out->uv_bridge_evidence_rejects++;
            continue;
        }
        if (atg_row_has_face_node_between(
                set, row_start, face_node, set->target[a].row,
                set->target[a].column, set->target[b].column) ||
            atg_row_has_face_node_between(
                set, row_start, face_node, set->target[c].row,
                set->target[c].column, set->target[d].column) ||
            !atg_quad_valid(
                set, set->target, a, b, c, d,
                opts->quad_normal_dot_min)) {
            out->uv_bridge_topology_rejects++;
            continue;
        }

        int32_t support[8] = {la, a, lc, c, b, rb, d, rd};
        double fit_rms = 0.0, fit_max = 0.0;
        if (atg_uvxyz_fit(
                set, target_u, support, 8, &fit_rms, &fit_max) != 0 ||
            fit_rms > opts->uv_bridge_fit_rms_max ||
            fit_max > opts->uv_bridge_fit_max_max) {
            out->uv_bridge_fit_rejects++;
            continue;
        }
        if (*io_nquad >= quad_capacity || *io_nquad > (size_t)INT32_MAX) {
            free(face_node);
            return -1;
        }
        if (top_open) {
            if (atg_set_edge(right, left, a, b) != 0) {
                free(face_node);
                return -1;
            }
            right_inferred[a] = 1;
            out->uv_bridge_horizontal_edges++;
        }
        if (bottom_open) {
            if (atg_set_edge(right, left, c, d) != 0) {
                free(face_node);
                return -1;
            }
            right_inferred[c] = 1;
            out->uv_bridge_horizontal_edges++;
        }
        AtgQuad *cell = &quad[*io_nquad];
        memset(cell, 0, sizeof *cell);
        cell->a = a;
        cell->b = b;
        cell->c = c;
        cell->d = d;
        cell->fit_rms = (float)fit_rms;
        cell->fit_max = (float)fit_max;
        cell->inferred = ATG_INFER_UV_BANK;
        cell->observed_edges = 2;
        quad_at_a[a] = (int32_t)*io_nquad;
        horizontal_face[a] = horizontal_face[c] = 1;
        vertical_face[a] = vertical_face[b] = 1;
        face_node[a] = face_node[b] = face_node[c] = face_node[d] = 1;
        (*io_nquad)++;
        out->uv_bridge_quads++;
        fit_sum2 += fit_rms * fit_rms;
        if (fit_max > fit_maximum) fit_maximum = fit_max;
        fit_count++;
        double span = top_gap > bottom_gap ? top_gap : bottom_gap;
        if (span > out->uv_bridge_u_span_max)
            out->uv_bridge_u_span_max = span;
    }
    free(face_node);
    out->uv_bridge_fit_rms = fit_count > 0
        ? sqrt(fit_sum2 / (double)fit_count) : 0.0;
    out->uv_bridge_fit_max = fit_maximum;
    return 0;
}

static int atg_mask_needs_split(uint8_t quadrant_mask)
{
    return quadrant_mask == UINT8_C(0x09) ||
           quadrant_mask == UINT8_C(0x06);
}

static int atg_cross_output_pair(
    const int32_t *source_head,
    const int32_t *source_next,
    const int32_t *vertex_component,
    int32_t source_a,
    int32_t source_b,
    int32_t *out_a,
    int32_t *out_b)
{
    int32_t best_a = -1, best_b = -1;
    for (int32_t a = source_head[source_a]; a >= 0; a = source_next[a])
        for (int32_t b = source_head[source_b]; b >= 0; b = source_next[b]) {
            if (vertex_component[a] == vertex_component[b]) continue;
            if (best_a < 0 || vertex_component[a] < vertex_component[best_a] ||
                (vertex_component[a] == vertex_component[best_a] &&
                 vertex_component[b] < vertex_component[best_b])) {
                best_a = a;
                best_b = b;
            }
        }
    if (best_a < 0) return 0;
    *out_a = best_a;
    *out_b = best_b;
    return 1;
}

static size_t atg_collect_parameter_bridges(
    const AtlasRibbonObservationSet *set,
    const int32_t *source_right,
    const int32_t *source_down,
    const int32_t *source_head,
    const int32_t *source_next,
    const AtlasTrackGrowResult *out,
    AtgParameterBridge *bridge,
    size_t *horizontal,
    size_t *vertical,
    size_t *point)
{
    size_t count = 0, nh = 0, nv = 0, np = 0;

    /* A source sample duplicated to split a point contact or disconnected fan
     * remains one parameter-space observation.  This couples its gauges while
     * deliberately leaving the published face topology split. */
    for (size_t source = 0; source < set->ntarget; source++) {
        int32_t a = source_head[source];
        if (a < 0) continue;
        for (int32_t b = source_next[a]; b >= 0; b = source_next[b]) {
            if (bridge != NULL) {
                bridge[count].a = a;
                bridge[count].b = b;
                bridge[count].delta = 0.0;
                bridge[count].kind = ATG_BRIDGE_POINT;
            }
            count++;
            np++;
            a = b;
        }
    }

    for (size_t source = 0; source < set->ntarget; source++) {
        int32_t target = source_right[source];
        int32_t a = -1, b = -1;
        if (target >= 0 && atg_cross_output_pair(
                source_head, source_next, out->vertex_component,
                (int32_t)source, target, &a, &b)) {
            if (bridge != NULL) {
                bridge[count].a = a;
                bridge[count].b = b;
                bridge[count].delta =
                    set->target[target].u - set->target[source].u;
                bridge[count].kind = ATG_BRIDGE_HORIZONTAL;
            }
            count++;
            nh++;
        }
        target = source_down[source];
        a = b = -1;
        if (target >= 0 && atg_cross_output_pair(
                source_head, source_next, out->vertex_component,
                (int32_t)source, target, &a, &b)) {
            if (bridge != NULL) {
                bridge[count].a = a;
                bridge[count].b = b;
                bridge[count].delta =
                    set->target[target].u - set->target[source].u;
                bridge[count].kind = ATG_BRIDGE_VERTICAL;
            }
            count++;
            nv++;
        }
    }
    if (horizontal != NULL) *horizontal = nh;
    if (vertical != NULL) *vertical = nv;
    if (point != NULL) *point = np;
    return count;
}

typedef struct {
    size_t component;
    size_t quads;
    double width;
    double height;
} AtgPackItem;

typedef struct {
    int32_t a;
    int32_t b;
    uint32_t weight;
    double desired;
} AtgChainConstraint;

typedef struct {
    int32_t chain_a;
    int32_t chain_b;
    int32_t source;
    double desired;
} AtgVerticalChainEdge;

static int atg_compare_pack_item(const void *left, const void *right)
{
    const AtgPackItem *a = (const AtgPackItem *)left;
    const AtgPackItem *b = (const AtgPackItem *)right;
    if (a->quads != b->quads) return a->quads > b->quads ? -1 : 1;
    return a->component < b->component ? -1 :
           (a->component > b->component ? 1 : 0);
}

static int atg_compare_chain_constraint(const void *left, const void *right)
{
    const AtgChainConstraint *a = (const AtgChainConstraint *)left;
    const AtgChainConstraint *b = (const AtgChainConstraint *)right;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b ? 1 : 0);
}

static int atg_compare_vertical_chain_edge(const void *left, const void *right)
{
    const AtgVerticalChainEdge *a =
        (const AtgVerticalChainEdge *)left;
    const AtgVerticalChainEdge *b =
        (const AtgVerticalChainEdge *)right;
    if (a->chain_a != b->chain_a)
        return a->chain_a < b->chain_a ? -1 : 1;
    if (a->chain_b != b->chain_b)
        return a->chain_b < b->chain_b ? -1 : 1;
    if (a->desired != b->desired)
        return a->desired < b->desired ? -1 : 1;
    return a->source < b->source ? -1 :
           (a->source > b->source ? 1 : 0);
}

/* A nearest physical point in the next row is not enough to establish sheet
 * continuation where two turns pass close together.  First integrate exact
 * arclength along every mutual horizontal chain.  Cross-row matches between
 * the same pair of chains must then imply one common relative U shift.  Cut
 * individual matches that disagree with the pair median before any quads or
 * connected components are constructed. */
static int atg_filter_vertical_chain_edges(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    const AtlasTrackGrowOptions *opts,
    const int32_t *right,
    const int32_t *left,
    int32_t *down,
    int32_t *up,
    int32_t *down_best,
    int32_t *up_best,
    AtlasTrackGrowResult *out)
{
    double residual_max = opts->vertical_chain_residual_max;
    if (!(residual_max > 0.0) || set->ntarget == 0) return 0;
    size_t n = set->ntarget;
    int32_t *chain = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*chain));
    double *chain_u = (double *)ARENA_ALLOC(
        arena, n * sizeof(*chain_u));
    for (size_t i = 0; i < n; i++) {
        chain[i] = -1;
        chain_u[i] = 0.0;
    }

    int32_t nchains = 0;
    for (int pass = 0; pass < 2; pass++) {
        for (size_t seed = 0; seed < n; seed++) {
            if (chain[seed] >= 0 ||
                (pass == 0 && left[seed] >= 0))
                continue;
            int32_t current = (int32_t)seed;
            double u = 0.0;
            while (current >= 0 && chain[(size_t)current] < 0) {
                chain[(size_t)current] = nchains;
                chain_u[(size_t)current] = u;
                int32_t next = right[(size_t)current];
                if (next >= 0)
                    u += atg_lattice_edge_span(set, current, next);
                current = next;
            }
            nchains++;
        }
    }
    if (nchains <= 0) return -1;

    size_t nedges = 0;
    for (size_t i = 0; i < n; i++)
        if (down[i] >= 0) nedges++;
    AtgVerticalChainEdge *edge = (AtgVerticalChainEdge *)ARENA_ALLOC(
        arena, (nedges ? nedges : 1) * sizeof(*edge));
    size_t at = 0;
    for (size_t i = 0; i < n; i++) {
        int32_t b = down[i];
        if (b < 0) continue;
        if (chain[i] < 0 || chain[(size_t)b] < 0) return -1;
        edge[at].chain_a = chain[i];
        edge[at].chain_b = chain[(size_t)b];
        edge[at].source = (int32_t)i;
        edge[at].desired = chain_u[i] - chain_u[(size_t)b];
        at++;
    }
    if (at != nedges) return -1;
    qsort(edge, nedges, sizeof(*edge),
          atg_compare_vertical_chain_edge);

    double sum2 = 0.0;
    double maximum = 0.0;
    for (size_t first = 0; first < nedges;) {
        size_t last = first + 1;
        while (last < nedges &&
               edge[last].chain_a == edge[first].chain_a &&
               edge[last].chain_b == edge[first].chain_b)
            last++;
        size_t count = last - first;
        size_t middle = first + count / 2;
        double median = edge[middle].desired;
        if ((count & 1) == 0)
            median = 0.5 * (edge[middle - 1].desired + median);
        out->vertical_chain_pairs++;
        for (size_t i = first; i < last; i++) {
            double residual = fabs(edge[i].desired - median);
            sum2 += residual * residual;
            if (residual > maximum) maximum = residual;
            out->vertical_chain_checks++;
            if (residual <= residual_max) continue;
            int32_t a = edge[i].source;
            int32_t b = down[(size_t)a];
            if (b < 0 || up[(size_t)b] != a) return -1;
            down[(size_t)a] = -1;
            up[(size_t)b] = -1;
            if (down_best[(size_t)a] == b)
                down_best[(size_t)a] = -1;
            if (up_best[(size_t)b] == a)
                up_best[(size_t)b] = -1;
            if (out->vertical_edges == 0) return -1;
            out->vertical_edges--;
            out->vertical_chain_cuts++;
        }
        first = last;
    }
    out->vertical_chain_residual_rms = nedges > 0
        ? sqrt(sum2 / (double)nedges) : 0.0;
    out->vertical_chain_residual_max = maximum;

    /* Pair-local agreement does not expose a wrong link that is internally
     * coherent but closes an impossible cycle through several row chains.
     * Solve one shift per chain, remove individual links whose residual cannot
     * fit that common layout, and repeat before constructing any cells. */
    size_t max_constraints = out->vertical_chain_pairs;
    AtgChainConstraint *constraint = (AtgChainConstraint *)ARENA_ALLOC(
        arena, (max_constraints ? max_constraints : 1) *
               sizeof(*constraint));
    size_t *offset = (size_t *)ARENA_ALLOC(
        arena, ((size_t)nchains + 1) * sizeof(*offset));
    int32_t *incident = (int32_t *)ARENA_ALLOC(
        arena, (max_constraints ? max_constraints * 2 : 1) *
               sizeof(*incident));
    size_t *cursor = (size_t *)ARENA_ALLOC(
        arena, (size_t)nchains * sizeof(*cursor));
    double *shift = (double *)ARENA_ALLOC(
        arena, (size_t)nchains * sizeof(*shift));
    double *initial_shift = (double *)ARENA_ALLOC(
        arena, (size_t)nchains * sizeof(*initial_shift));
    uint8_t *seen = (uint8_t *)ARENA_ALLOC(
        arena, (size_t)nchains);
    uint8_t *anchor = (uint8_t *)ARENA_ALLOC(
        arena, (size_t)nchains);
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, (size_t)nchains * sizeof(*queue));
    int filter_iterations = opts->parameter_iterations;
    if (filter_iterations < 160) filter_iterations = 160;
    if (filter_iterations > 2000) filter_iterations = 2000;

    for (int round = 0; round < 16; round++) {
        size_t nconstraints = 0;
        for (size_t first = 0; first < nedges;) {
            size_t last = first + 1;
            while (last < nedges &&
                   edge[last].chain_a == edge[first].chain_a &&
                   edge[last].chain_b == edge[first].chain_b)
                last++;
            double desired_sum = 0.0;
            uint32_t weight = 0;
            for (size_t i = first; i < last; i++) {
                if (down[(size_t)edge[i].source] < 0) continue;
                desired_sum += edge[i].desired;
                weight++;
            }
            if (weight > 0) {
                if (nconstraints >= max_constraints) return -1;
                constraint[nconstraints].a = edge[first].chain_a;
                constraint[nconstraints].b = edge[first].chain_b;
                constraint[nconstraints].weight = weight;
                constraint[nconstraints].desired =
                    desired_sum / (double)weight;
                nconstraints++;
            }
            first = last;
        }

        memset(offset, 0, ((size_t)nchains + 1) * sizeof(*offset));
        for (size_t e = 0; e < nconstraints; e++) {
            offset[(size_t)constraint[e].a + 1]++;
            offset[(size_t)constraint[e].b + 1]++;
        }
        for (int32_t c = 1; c <= nchains; c++)
            offset[(size_t)c] += offset[(size_t)c - 1];
        memcpy(cursor, offset, (size_t)nchains * sizeof(*cursor));
        for (size_t e = 0; e < nconstraints; e++) {
            incident[cursor[(size_t)constraint[e].a]++] = (int32_t)e;
            incident[cursor[(size_t)constraint[e].b]++] = (int32_t)e;
        }

        memset(shift, 0, (size_t)nchains * sizeof(*shift));
        memset(initial_shift, 0,
               (size_t)nchains * sizeof(*initial_shift));
        memset(seen, 0, (size_t)nchains);
        memset(anchor, 0, (size_t)nchains);
        for (int32_t seed = 0; seed < nchains; seed++) {
            if (seen[(size_t)seed]) continue;
            seen[(size_t)seed] = 1;
            anchor[(size_t)seed] = 1;
            size_t head = 0, tail = 0;
            queue[tail++] = seed;
            while (head < tail) {
                int32_t current = queue[head++];
                for (size_t pos = offset[(size_t)current];
                     pos < offset[(size_t)current + 1]; pos++) {
                    AtgChainConstraint *link =
                        &constraint[(size_t)incident[pos]];
                    int32_t next;
                    double value;
                    if (link->a == current) {
                        next = link->b;
                        value = shift[(size_t)current] + link->desired;
                    } else {
                        next = link->a;
                        value = shift[(size_t)current] - link->desired;
                    }
                    if (seen[(size_t)next]) continue;
                    seen[(size_t)next] = 1;
                    shift[(size_t)next] = value;
                    queue[tail++] = next;
                }
            }
        }
        memcpy(initial_shift, shift,
               (size_t)nchains * sizeof(*shift));

        for (int iteration = 0;
             iteration < filter_iterations; iteration++) {
            double max_change = 0.0;
            for (int32_t c = 0; c < nchains; c++) {
                double numerator = anchor[(size_t)c]
                                 ? initial_shift[(size_t)c] : 0.0;
                double denominator = anchor[(size_t)c] ? 1.0 : 0.0;
                for (size_t pos = offset[(size_t)c];
                     pos < offset[(size_t)c + 1]; pos++) {
                    AtgChainConstraint *link =
                        &constraint[(size_t)incident[pos]];
                    double weight = (double)link->weight;
                    if (link->a == c)
                        numerator += weight *
                            (shift[(size_t)link->b] - link->desired);
                    else
                        numerator += weight *
                            (shift[(size_t)link->a] + link->desired);
                    denominator += weight;
                }
                if (!(denominator > 0.0)) continue;
                double next = numerator / denominator;
                double change = fabs(next - shift[(size_t)c]);
                if (change > max_change) max_change = change;
                shift[(size_t)c] = next;
            }
            if (max_change < 1.0e-9) break;
        }

        size_t checked = 0;
        size_t cut_this_round = 0;
        double global_sum2 = 0.0;
        double global_max = 0.0;
        for (size_t i = 0; i < nedges; i++) {
            int32_t a = edge[i].source;
            int32_t b = down[(size_t)a];
            if (b < 0) continue;
            double residual = fabs(
                (shift[(size_t)edge[i].chain_b] -
                 shift[(size_t)edge[i].chain_a]) - edge[i].desired);
            global_sum2 += residual * residual;
            if (residual > global_max) global_max = residual;
            checked++;
            out->vertical_global_checks++;
            if (residual <= residual_max) continue;
            if (up[(size_t)b] != a) return -1;
            down[(size_t)a] = -1;
            up[(size_t)b] = -1;
            if (down_best[(size_t)a] == b)
                down_best[(size_t)a] = -1;
            if (up_best[(size_t)b] == a)
                up_best[(size_t)b] = -1;
            if (out->vertical_edges == 0) return -1;
            out->vertical_edges--;
            out->vertical_global_cuts++;
            cut_this_round++;
        }
        out->vertical_global_rounds = round + 1;
        out->vertical_global_residual_rms = checked > 0
            ? sqrt(global_sum2 / (double)checked) : 0.0;
        out->vertical_global_residual_max = global_max;
        if (cut_this_round == 0) break;
    }
    return 0;
}

/* Once topology has been recovered from local XYZ evidence, the checkpoint U
 * values are no longer a valid parameterization: they are exactly the bad
 * layout that fragmented the sheet.  The fitting domain is the discrete quad
 * lattice, not a separate continuous ruler for every row: quantize each
 * measured H span to whole observation cells so skipped samples remain empty
 * columns for the later fit, keep V edges isoparametric, solve the row-chain
 * offsets, and shelf-pack whole connected sheets. */
static int atg_build_intrinsic_parameterization(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    const AtlasTrackGrowOptions *opts,
    const int32_t *mesh_right,
    const int32_t *mesh_left,
    const int32_t *mesh_down,
    const int32_t *mesh_up,
    AtlasTrackGrowResult *out)
{
    if (out->nv == 0) return 0;
    int32_t *chain = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*chain));
    for (size_t i = 0; i < out->nv; i++) chain[i] = -1;
    for (size_t i = 0; i < out->nv; i++) {
        const AtlasRibbonTarget *target =
            &set->target[out->source_target[i]];
        out->global_uv[i * 2] = 0.0;
        out->global_uv[i * 2 + 1] = target->v;
    }

    /* Each mutual H graph is a simple row chain on the quad lattice.  Long
     * evidence links retain an integer number of missing lattice cells.  The
     * only remaining unknown is one scalar lattice shift per chain. */
    int32_t nchains = 0;
    for (size_t pass = 0; pass < 2; pass++) {
        for (size_t seed = 0; seed < out->nv; seed++) {
            if (chain[seed] >= 0 ||
                (pass == 0 && mesh_left[seed] >= 0))
                continue;
            int32_t current = (int32_t)seed;
            double u = 0.0;
            while (current >= 0 && chain[(size_t)current] < 0) {
                chain[(size_t)current] = nchains;
                out->global_uv[(size_t)current * 2] = u;
                int32_t next = mesh_right[(size_t)current];
                if (next >= 0)
                    u += atg_lattice_edge_span(
                        set,
                        out->source_target[(size_t)current],
                        out->source_target[(size_t)next]);
                current = next;
            }
            nchains++;
        }
    }
    if (nchains <= 0) return -1;

    size_t raw_count = 0;
    for (size_t i = 0; i < out->nv; i++)
        if (mesh_down[i] >= 0) raw_count++;
    AtgChainConstraint *constraint = (AtgChainConstraint *)ARENA_ALLOC(
        arena, (raw_count ? raw_count : 1) * sizeof(*constraint));
    size_t raw_at = 0;
    for (size_t i = 0; i < out->nv; i++) {
        int32_t down = mesh_down[i];
        if (down < 0) continue;
        int32_t a = chain[i], b = chain[(size_t)down];
        if (a == b) continue;
        AtgChainConstraint *edge = &constraint[raw_at++];
        edge->a = a;
        edge->b = b;
        edge->weight = 1;
        edge->desired =
            out->global_uv[i * 2] -
            out->global_uv[(size_t)down * 2];
    }
    qsort(constraint, raw_at, sizeof(*constraint),
          atg_compare_chain_constraint);
    size_t nconstraints = 0;
    for (size_t i = 0; i < raw_at; i++) {
        AtgChainConstraint *source = &constraint[i];
        if (nconstraints > 0 &&
            constraint[nconstraints - 1].a == source->a &&
            constraint[nconstraints - 1].b == source->b) {
            AtgChainConstraint *target = &constraint[nconstraints - 1];
            target->desired += source->desired;
            target->weight++;
        } else {
            constraint[nconstraints++] = *source;
        }
    }
    for (size_t i = 0; i < nconstraints; i++)
        constraint[i].desired /= (double)constraint[i].weight;

    size_t *offset = (size_t *)ARENA_CALLOC(
        arena, (size_t)nchains + 1, sizeof(*offset));
    for (size_t e = 0; e < nconstraints; e++) {
        offset[(size_t)constraint[e].a + 1]++;
        offset[(size_t)constraint[e].b + 1]++;
    }
    for (int32_t c = 1; c <= nchains; c++)
        offset[(size_t)c] += offset[(size_t)c - 1];
    int32_t *incident = (int32_t *)ARENA_ALLOC(
        arena, offset[(size_t)nchains] * sizeof(*incident));
    size_t *cursor = (size_t *)ARENA_ALLOC(
        arena, (size_t)nchains * sizeof(*cursor));
    memcpy(cursor, offset, (size_t)nchains * sizeof(*cursor));
    for (size_t e = 0; e < nconstraints; e++) {
        incident[cursor[(size_t)constraint[e].a]++] = (int32_t)e;
        incident[cursor[(size_t)constraint[e].b]++] = (int32_t)e;
    }

    double *shift = (double *)ARENA_CALLOC(
        arena, (size_t)nchains, sizeof(*shift));
    double *initial_shift = (double *)ARENA_CALLOC(
        arena, (size_t)nchains, sizeof(*initial_shift));
    size_t *chain_observations = (size_t *)ARENA_CALLOC(
        arena, (size_t)nchains, sizeof(*chain_observations));
    uint8_t *chain_seen = (uint8_t *)ARENA_CALLOC(
        arena, (size_t)nchains, 1);
    uint8_t *chain_anchor = (uint8_t *)ARENA_CALLOC(
        arena, (size_t)nchains, 1);
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, (size_t)nchains * sizeof(*queue));
    /* The cube-sheet layout fixes the otherwise free integer offset of every
     * horizontal chain.  Average its per-sample U offsets first; vertical
     * lattice constraints may polish them, but they may no longer place a
     * disconnected row fragment on an arbitrary coincident wrap. */
    for (size_t i = 0; i < out->nv; i++) {
        int32_t c = chain[i];
        const AtlasRibbonTarget *target =
            &set->target[out->source_target[i]];
        if (c < 0 || !isfinite(target->u)) return -1;
        initial_shift[(size_t)c] +=
            target->u - out->global_uv[i * 2];
        chain_observations[(size_t)c]++;
    }
    for (int32_t c = 0; c < nchains; c++) {
        if (chain_observations[(size_t)c] == 0) return -1;
        initial_shift[(size_t)c] /=
            (double)chain_observations[(size_t)c];
        shift[(size_t)c] = initial_shift[(size_t)c];
    }
    for (int32_t seed = 0; seed < nchains; seed++) {
        if (chain_seen[(size_t)seed]) continue;
        out->parameter_islands++;
        chain_seen[(size_t)seed] = 1;
        chain_anchor[(size_t)seed] = 1;
        size_t head = 0, tail = 0;
        queue[tail++] = seed;
        while (head < tail) {
            int32_t current = queue[head++];
            for (size_t at = offset[(size_t)current];
                 at < offset[(size_t)current + 1]; at++) {
                AtgChainConstraint *edge =
                    &constraint[(size_t)incident[at]];
                int32_t next;
                if (edge->a == current) {
                    next = edge->b;
                } else {
                    next = edge->a;
                }
                if (chain_seen[(size_t)next]) continue;
                chain_seen[(size_t)next] = 1;
                queue[tail++] = next;
            }
        }
    }
    double max_change = 0.0;
    int iterations = 0;
    for (int iteration = 0;
         iteration < opts->parameter_iterations;
         iteration++) {
        max_change = 0.0;
        for (int32_t c = 0; c < nchains; c++) {
            double anchor_weight = opts->parameter_stay_weight *
                (double)chain_observations[(size_t)c];
            if (!(anchor_weight > 0.0) && chain_anchor[(size_t)c])
                anchor_weight = 1.0;
            double numerator = anchor_weight * initial_shift[(size_t)c];
            double denominator = anchor_weight;
            for (size_t at = offset[(size_t)c];
                 at < offset[(size_t)c + 1]; at++) {
                AtgChainConstraint *edge =
                    &constraint[(size_t)incident[at]];
                double weight = (double)edge->weight;
                if (edge->a == c)
                    numerator += weight *
                        (shift[(size_t)edge->b] - edge->desired);
                else
                    numerator += weight *
                        (shift[(size_t)edge->a] + edge->desired);
                denominator += weight;
            }
            if (!(denominator > 0.0)) continue;
            double next = numerator / denominator;
            double change = fabs(next - shift[(size_t)c]);
            if (change > max_change) max_change = change;
            shift[(size_t)c] = next;
        }
        iterations = iteration + 1;
        if (max_change < 1.0e-9) break;
    }
    out->parameter_iterations = iterations;
    out->parameter_max_change = max_change;
    for (int32_t c = 0; c < nchains; c++)
        shift[(size_t)c] = set->observation_du *
            nearbyint(shift[(size_t)c] / set->observation_du);
    for (size_t i = 0; i < out->nv; i++)
        out->global_uv[i * 2] += shift[(size_t)chain[i]];

    double hsum2 = 0.0, hmax = 0.0;
    double vsum2 = 0.0, vmax = 0.0;
    size_t hn = 0, vn = 0;
    double anchor_sum2 = 0.0, anchor_max = 0.0;
    for (size_t i = 0; i < out->nv; i++) {
        int32_t neighbor = mesh_right[i];
        if (neighbor >= 0) {
            double distance = atg_edge_distance(
                set->target, out->source_target[i],
                out->source_target[(size_t)neighbor]);
            double residual = fabs(
                (out->global_uv[(size_t)neighbor * 2] -
                 out->global_uv[i * 2]) - distance);
            hsum2 += residual * residual;
            if (residual > hmax) hmax = residual;
            hn++;
        }
        neighbor = mesh_down[i];
        if (neighbor >= 0) {
            double residual = fabs(
                (out->global_uv[(size_t)neighbor * 2] -
                 out->global_uv[i * 2]));
            static const double limit[9] = {
                0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0
            };
            size_t bin = 0;
            while (bin < 9 && residual > limit[bin]) bin++;
            out->vertical_residual_hist[bin]++;
            vsum2 += residual * residual;
            if (residual > vmax) vmax = residual;
            vn++;
        }
        double residual = fabs(
            shift[(size_t)chain[i]] -
            initial_shift[(size_t)chain[i]]);
        anchor_sum2 += residual * residual;
        if (residual > anchor_max) anchor_max = residual;
    }
    out->horizontal_residual_rms =
        hn > 0 ? sqrt(hsum2 / (double)hn) : 0.0;
    out->horizontal_residual_max = hmax;
    out->vertical_residual_rms =
        vn > 0 ? sqrt(vsum2 / (double)vn) : 0.0;
    out->vertical_residual_max = vmax;
    out->parameter_anchor_rms =
        sqrt(anchor_sum2 / (double)out->nv);
    out->parameter_anchor_max = anchor_max;

    size_t selected = 0;
    for (size_t ci = 0; ci < out->ncomponents; ci++)
        if (out->component[ci].selected) selected++;
    AtgPackItem *item = (AtgPackItem *)ARENA_ALLOC(
        arena, (selected ? selected : 1) * sizeof(*item));
    size_t nitem = 0;
    double max_width = 0.0;
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        AtlasTrackGrowComponent *component = &out->component[ci];
        if (!component->selected) continue;
        size_t begin = component->first_vertex;
        size_t end = begin + component->vertices;
        double umin = DBL_MAX, umax = -DBL_MAX;
        double vmin = DBL_MAX, vlocal_max = -DBL_MAX;
        component->winding_min = INT32_MAX;
        component->winding_max = INT32_MIN;
        component->rank_min = UINT64_MAX;
        component->rank_max = 0;
        for (size_t i = begin; i < end; i++) {
            double u = out->global_uv[i * 2];
            double v = out->global_uv[i * 2 + 1];
            if (u < umin) umin = u;
            if (u > umax) umax = u;
            if (v < vmin) vmin = v;
            if (v > vlocal_max) vlocal_max = v;
            if (out->vertex_winding[i] < component->winding_min)
                component->winding_min = out->vertex_winding[i];
            if (out->vertex_winding[i] > component->winding_max)
                component->winding_max = out->vertex_winding[i];
            if (out->vertex_rank[i] < component->rank_min)
                component->rank_min = out->vertex_rank[i];
            if (out->vertex_rank[i] > component->rank_max)
                component->rank_max = out->vertex_rank[i];
        }
        if (umin == DBL_MAX || vmin == DBL_MAX) return -1;
        component->global_u0 = umin;
        component->global_v0 = vmin;
        component->local_u_span = umax - umin;
        component->local_v_span = vlocal_max - vmin;
        if (component->local_u_span > out->local_u_span_max)
            out->local_u_span_max = component->local_u_span;
        for (size_t i = begin; i < end; i++) {
            out->local_uv[i * 2] = out->global_uv[i * 2] - umin;
            out->local_uv[i * 2 + 1] = out->global_uv[i * 2 + 1] - vmin;
        }
        item[nitem].component = ci;
        item[nitem].quads = component->quads;
        item[nitem].width = component->local_u_span;
        item[nitem].height = component->local_v_span;
        if (item[nitem].width > max_width)
            max_width = item[nitem].width;
        nitem++;
    }
    if (nitem != selected) return -1;
    qsort(item, nitem, sizeof(*item), atg_compare_pack_item);

    const double margin = 16.0;
    double target_width = max_width > 65536.0 ? max_width : 65536.0;
    double x = 0.0, y = 0.0, shelf_height = 0.0;
    int band = 0;
    out->atlas_u_span = 0.0;
    out->atlas_v_span = 0.0;
    for (size_t at = 0; at < nitem; at++) {
        AtgPackItem *place = &item[at];
        AtlasTrackGrowComponent *component =
            &out->component[place->component];
        if (x > 0.0 && x + place->width > target_width) {
            y += shelf_height + margin;
            x = 0.0;
            shelf_height = 0.0;
            band++;
        }
        component->atlas_u0 = x;
        component->atlas_v0 = y;
        component->atlas_band = band;
        size_t begin = component->first_vertex;
        size_t end = begin + component->vertices;
        for (size_t i = begin; i < end; i++) {
            out->uv[i * 2] = out->local_uv[i * 2] + x;
            out->uv[i * 2 + 1] = out->local_uv[i * 2 + 1] + y;
        }
        double right = x + place->width;
        double bottom = y + place->height;
        if (right > out->atlas_u_span) out->atlas_u_span = right;
        if (bottom > out->atlas_v_span) out->atlas_v_span = bottom;
        x = right + margin;
        if (place->height > shelf_height) shelf_height = place->height;
    }
    out->atlas_bands = nitem > 0 ? band + 1 : 0;
    return 0;
}

static int atg_build_parameterization(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    const ScaffoldCalib *calibration,
    const AtlasTrackGrowOptions *opts,
    const int32_t *source_right,
    const int32_t *source_down,
    const int32_t *mesh_right,
    const int32_t *mesh_left,
    const int32_t *mesh_down,
    const int32_t *mesh_up,
    AtlasTrackGrowResult *out)
{
    if (out->nv == 0) return 0;
    if (calibration == NULL) return -1;
    if (opts->physical_neighbours &&
        !opts->preserve_layout_parameterization)
        return atg_build_intrinsic_parameterization(
            arena, set, opts, mesh_right, mesh_left,
            mesh_down, mesh_up, out);

    int32_t *source_head = (int32_t *)ARENA_ALLOC(
        arena, set->ntarget * sizeof(*source_head));
    int32_t *source_next = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*source_next));
    for (size_t i = 0; i < set->ntarget; i++) source_head[i] = -1;
    for (size_t i = 0; i < out->nv; i++) {
        int32_t source = out->source_target[i];
        if (source < 0 || (size_t)source >= set->ntarget) return -1;
        source_next[i] = source_head[source];
        source_head[source] = (int32_t)i;
    }

    size_t bridge_count = atg_collect_parameter_bridges(
        set, source_right, source_down, source_head, source_next, out,
        NULL, NULL, NULL, NULL);
    AtgParameterBridge *bridge = bridge_count > 0
        ? (AtgParameterBridge *)ARENA_ALLOC(
            arena, bridge_count * sizeof(*bridge))
        : NULL;
    size_t collected = atg_collect_parameter_bridges(
        set, source_right, source_down, source_head, source_next, out,
        bridge, &out->parameter_horizontal_bridges,
        &out->parameter_vertical_bridges,
        &out->parameter_point_bridges);
    if (collected != bridge_count) return -1;
    out->parameter_bridge_edges = bridge_count;

    /* Index every soft bridge at both endpoints so it participates in the
     * same normal-equation sweep as the cycle edges.  Every edge preserves its
     * difference in the checkpoint's packed U gauge.  Consequently the exact
     * checkpoint coordinates are a zero-residual global solution, including
     * across split point contacts and evidence-supported component bridges. */
    size_t *bridge_offset = (size_t *)ARENA_CALLOC(
        arena, out->nv + 1, sizeof(*bridge_offset));
    for (size_t e = 0; e < bridge_count; e++) {
        bridge_offset[(size_t)bridge[e].a + 1]++;
        bridge_offset[(size_t)bridge[e].b + 1]++;
    }
    for (size_t i = 1; i <= out->nv; i++)
        bridge_offset[i] += bridge_offset[i - 1];
    if (bridge_count > SIZE_MAX / 2 ||
        bridge_offset[out->nv] != bridge_count * 2)
        return -1;
    int32_t *bridge_incident = bridge_count > 0
        ? (int32_t *)ARENA_ALLOC(
            arena, bridge_count * 2 * sizeof(*bridge_incident))
        : NULL;
    size_t *bridge_cursor = (size_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*bridge_cursor));
    memcpy(bridge_cursor, bridge_offset,
           out->nv * sizeof(*bridge_cursor));
    for (size_t e = 0; e < bridge_count; e++) {
        bridge_incident[bridge_cursor[(size_t)bridge[e].a]++] = (int32_t)e;
        bridge_incident[bridge_cursor[(size_t)bridge[e].b]++] = (int32_t)e;
    }

    UnionFind component_uf = UF_new(arena, (int32_t)out->ncomponents);
    for (size_t e = 0; e < bridge_count; e++) {
        int32_t ca = out->vertex_component[bridge[e].a];
        int32_t cb = out->vertex_component[bridge[e].b];
        if (ca != cb) uf_union(&component_uf, ca, cb);
    }
    uint8_t *island_seen = (uint8_t *)ARENA_CALLOC(
        arena, out->ncomponents, 1);
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        if (!out->component[ci].selected) continue;
        int32_t root = uf_find(&component_uf, (int32_t)ci);
        if (!island_seen[root]) {
            island_seen[root] = 1;
            out->parameter_islands++;
        }
    }

    double *initial = (double *)ARENA_ALLOC(
        arena, out->nv * sizeof(*initial));
    for (size_t i = 0; i < out->nv; i++) {
        const AtlasRibbonTarget *target =
            &set->target[out->source_target[i]];
        initial[i] = target->u;
        if (!isfinite(initial[i]) || !isfinite(target->v)) return -1;
        out->global_uv[i * 2] = initial[i];
        out->global_uv[i * 2 + 1] = target->v;
    }

    double max_change = 0.0;
    int iterations = 0;
    for (int iteration = 0;
         iteration < opts->parameter_iterations;
         iteration++) {
        max_change = 0.0;
        for (size_t oi = 0; oi < out->nv; oi++) {
            double numerator =
                opts->parameter_stay_weight * initial[oi];
            double denominator = opts->parameter_stay_weight;

            int32_t neighbor = mesh_left[oi];
            if (neighbor >= 0) {
                double delta = initial[oi] - initial[(size_t)neighbor];
                numerator += out->global_uv[(size_t)neighbor * 2] + delta;
                denominator += 1.0;
            }
            neighbor = mesh_right[oi];
            if (neighbor >= 0) {
                double delta = initial[oi] - initial[(size_t)neighbor];
                numerator += out->global_uv[(size_t)neighbor * 2] + delta;
                denominator += 1.0;
            }
            neighbor = mesh_up[oi];
            if (neighbor >= 0) {
                double delta = initial[oi] - initial[(size_t)neighbor];
                numerator += out->global_uv[(size_t)neighbor * 2] + delta;
                denominator += 1.0;
            }
            neighbor = mesh_down[oi];
            if (neighbor >= 0) {
                double delta = initial[oi] - initial[(size_t)neighbor];
                numerator += out->global_uv[(size_t)neighbor * 2] + delta;
                denominator += 1.0;
            }
            for (size_t at = bridge_offset[oi];
                 at < bridge_offset[oi + 1]; at++) {
                const AtgParameterBridge *edge =
                    &bridge[(size_t)bridge_incident[at]];
                double desired;
                if ((size_t)edge->a == oi)
                    desired = out->global_uv[(size_t)edge->b * 2] -
                              edge->delta;
                else
                    desired = out->global_uv[(size_t)edge->a * 2] +
                              edge->delta;
                numerator += opts->parameter_bridge_weight * desired;
                denominator += opts->parameter_bridge_weight;
            }
            if (!(denominator > 0.0)) continue;
            double next = numerator / denominator;
            double old = out->global_uv[oi * 2];
            double change = fabs(next - old);
            if (change > max_change) max_change = change;
            out->global_uv[oi * 2] = next;
        }
        iterations = iteration + 1;
        if (max_change < 1.0e-7) break;
    }
    out->parameter_iterations = iterations;
    out->parameter_max_change = max_change;

    double hsum2 = 0.0, hmax = 0.0;
    double vsum2 = 0.0, vmax = 0.0;
    size_t hn = 0, vn = 0;
    for (size_t oi = 0; oi < out->nv; oi++) {
        int32_t source = out->source_target[oi];
        int32_t neighbor = mesh_right[oi];
        if (neighbor >= 0) {
            double step = atg_edge_distance(
                set->target, source, out->source_target[neighbor]);
            double residual = fabs(
                (out->global_uv[(size_t)neighbor * 2] -
                 out->global_uv[oi * 2]) - step);
            hsum2 += residual * residual;
            if (residual > hmax) hmax = residual;
            hn++;
        }
        neighbor = mesh_down[oi];
        if (neighbor >= 0) {
            double expected = opts->preserve_layout_parameterization
                ? initial[(size_t)neighbor] - initial[oi] : 0.0;
            double residual = fabs(
                out->global_uv[(size_t)neighbor * 2] -
                out->global_uv[oi * 2] - expected);
            vsum2 += residual * residual;
            if (residual > vmax) vmax = residual;
            vn++;
        }
    }
    double bridge_sum2 = 0.0, bridge_max = 0.0;
    for (size_t e = 0; e < bridge_count; e++) {
        double residual = fabs(
            (out->global_uv[(size_t)bridge[e].b * 2] -
             out->global_uv[(size_t)bridge[e].a * 2]) - bridge[e].delta);
        bridge_sum2 += residual * residual;
        if (residual > bridge_max) bridge_max = residual;
        if (bridge[e].kind == ATG_BRIDGE_HORIZONTAL) {
            hsum2 += residual * residual;
            if (residual > hmax) hmax = residual;
            hn++;
        } else if (bridge[e].kind == ATG_BRIDGE_VERTICAL) {
            vsum2 += residual * residual;
            if (residual > vmax) vmax = residual;
            vn++;
        }
    }
    out->horizontal_residual_rms =
        hn > 0 ? sqrt(hsum2 / (double)hn) : 0.0;
    out->horizontal_residual_max = hmax;
    out->vertical_residual_rms =
        vn > 0 ? sqrt(vsum2 / (double)vn) : 0.0;
    out->vertical_residual_max = vmax;
    out->parameter_bridge_residual_rms = bridge_count > 0
        ? sqrt(bridge_sum2 / (double)bridge_count) : 0.0;
    out->parameter_bridge_residual_max = bridge_max;

    double anchor_sum2 = 0.0, anchor_max = 0.0;
    for (size_t i = 0; i < out->nv; i++) {
        double residual = fabs(out->global_uv[i * 2] - initial[i]);
        anchor_sum2 += residual * residual;
        if (residual > anchor_max) anchor_max = residual;
    }
    out->parameter_anchor_rms = sqrt(anchor_sum2 / (double)out->nv);
    out->parameter_anchor_max = anchor_max;

    double global_umin = DBL_MAX, global_umax = -DBL_MAX;
    double global_vmin = DBL_MAX, global_vmax = -DBL_MAX;
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        AtlasTrackGrowComponent *component = &out->component[ci];
        if (!component->selected) continue;
        size_t begin = component->first_vertex;
        size_t end = begin + component->vertices;
        double umin = DBL_MAX, umax = -DBL_MAX;
        double vmin = DBL_MAX, vlocal_max = -DBL_MAX;
        component->winding_min = INT32_MAX;
        component->winding_max = INT32_MIN;
        component->rank_min = UINT64_MAX;
        component->rank_max = 0;
        for (size_t i = begin; i < end; i++) {
            double u = out->global_uv[i * 2];
            double v = out->global_uv[i * 2 + 1];
            if (u < umin) umin = u;
            if (u > umax) umax = u;
            if (v < vmin) vmin = v;
            if (v > vlocal_max) vlocal_max = v;
            if (out->vertex_winding[i] < component->winding_min)
                component->winding_min = out->vertex_winding[i];
            if (out->vertex_winding[i] > component->winding_max)
                component->winding_max = out->vertex_winding[i];
            if (out->vertex_rank[i] < component->rank_min)
                component->rank_min = out->vertex_rank[i];
            if (out->vertex_rank[i] > component->rank_max)
                component->rank_max = out->vertex_rank[i];
        }
        if (umin == DBL_MAX || vmin == DBL_MAX) return -1;
        component->global_u0 = umin;
        component->global_v0 = vmin;
        component->local_u_span = umax - umin;
        component->local_v_span = vlocal_max - vmin;
        if (component->local_u_span > out->local_u_span_max)
            out->local_u_span_max = component->local_u_span;
        for (size_t i = begin; i < end; i++) {
            out->local_uv[i * 2] = out->global_uv[i * 2] - umin;
            out->local_uv[i * 2 + 1] = out->global_uv[i * 2 + 1] - vmin;
        }
        if (umin < global_umin) global_umin = umin;
        if (umax > global_umax) global_umax = umax;
        if (vmin < global_vmin) global_vmin = vmin;
        if (vlocal_max > global_vmax) global_vmax = vlocal_max;
    }
    if (global_umin == DBL_MAX || global_vmin == DBL_MAX) return -1;
    out->atlas_u_span = global_umax - global_umin;
    out->atlas_v_span = global_vmax - global_vmin;
    out->atlas_bands = 1;
    for (size_t i = 0; i < out->nv; i++) {
        out->uv[i * 2] = out->global_uv[i * 2] - global_umin;
        out->uv[i * 2 + 1] = out->global_uv[i * 2 + 1] - global_vmin;
    }
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        AtlasTrackGrowComponent *component = &out->component[ci];
        if (!component->selected) continue;
        component->atlas_u0 = component->global_u0 - global_umin;
        component->atlas_v0 = component->global_v0 - global_vmin;
        component->atlas_band = 0;
    }

    for (size_t a = 0; a < out->ncomponents; a++) {
        const AtlasTrackGrowComponent *ca = &out->component[a];
        if (!ca->selected) continue;
        for (size_t b = a + 1; b < out->ncomponents; b++) {
            const AtlasTrackGrowComponent *cb = &out->component[b];
            if (!cb->selected) continue;
            if (atg_intervals_overlap(
                    ca->global_u0,
                    ca->global_u0 + ca->local_u_span,
                    cb->global_u0,
                    cb->global_u0 + cb->local_u_span) &&
                atg_intervals_overlap(
                    ca->global_v0,
                    ca->global_v0 + ca->local_v_span,
                    cb->global_v0,
                    cb->global_v0 + cb->local_v_span)) {
                out->atlas_placed_overlaps++;
                if (ca->winding_min <= cb->winding_max &&
                    cb->winding_min <= ca->winding_max)
                    out->atlas_winding_overlaps++;
            }
        }
    }
    return 0;
}

static double atg_triangle_area(
    const double *xyz,
    int32_t ia,
    int32_t ib,
    int32_t ic)
{
    const double *a = &xyz[(size_t)ia * 3];
    const double *b = &xyz[(size_t)ib * 3];
    const double *c = &xyz[(size_t)ic * 3];
    double ab[3], ac[3], cross[3];
    for (int d = 0; d < 3; d++) {
        ab[d] = b[d] - a[d];
        ac[d] = c[d] - a[d];
    }
    atg_cross3(ab, ac, cross);
    return 0.5 * sqrt(atg_dot3(cross, cross));
}

int AtlasTrackGrow_build(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    const ScaffoldCalib *calibration,
    const AtlasTrackGrowOptions *opts,
    AtlasTrackGrowResult *out)
{
    if (arena == NULL || set == NULL || calibration == NULL || out == NULL ||
        !atg_options_valid(opts) ||
        !isfinite(calibration->spiral_a) ||
        !isfinite(calibration->spiral_b) ||
        (calibration->sense != -1 && calibration->sense != 1) ||
        (set->ntarget > 0 && set->target == NULL) ||
        (set->nlayer_samples > 0 && set->layer_sample == NULL) ||
        !isfinite(set->observation_du) ||
        !(set->observation_du > 0.0) ||
        set->ntarget > (size_t)INT32_MAX ||
        set->nlayer_samples > (size_t)INT32_MAX ||
        set->nrows > (size_t)INT32_MAX ||
        !atg_targets_sorted(set))
        return -1;
    memset(out, 0, sizeof *out);
    out->targets = set->ntarget;
    out->selected_row_min = INT32_MAX;
    out->selected_row_max = INT32_MIN;
    out->selected_column_min = INT32_MAX;
    out->selected_column_max = INT32_MIN;
    out->winding_min = INT32_MAX;
    out->winding_max = INT32_MIN;
    out->winding_direction = atg_winding_direction(set);
    if (set->ntarget == 0 || set->nrows == 0) {
        out->selected_row_min = out->selected_row_max = 0;
        out->selected_column_min = out->selected_column_max = 0;
        out->winding_min = out->winding_max = 0;
        return 0;
    }

    size_t n = set->ntarget;
    int32_t *target_layer = NULL;
    if (atg_map_target_layers(arena, set, &target_layer) != 0) return -1;
    double *target_u = (double *)ARENA_ALLOC(
        arena, n * sizeof(*target_u));
    for (size_t i = 0; i < n; i++) {
        target_u[i] = NAN;
        if (set->target[i].accepted) {
            out->accepted_targets++;
            int32_t layer = target_layer[i];
            if (layer < 0 || (size_t)layer >= set->nlayer_samples)
                return -1;
            target_u[i] = set->target[i].u;
            if (!isfinite(target_u[i])) return -1;
        } else {
            out->conflict_targets++;
        }
    }

    size_t *row_start = (size_t *)ARENA_ALLOC(
        arena, (set->nrows + 1) * sizeof(*row_start));
    size_t cursor = 0;
    for (size_t row = 0; row < set->nrows; row++) {
        row_start[row] = cursor;
        while (cursor < n &&
               set->target[cursor].row == (int32_t)row)
            cursor++;
    }
    row_start[set->nrows] = cursor;
    if (cursor != n) return -1;

    int32_t *right_best = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*right_best));
    int32_t *left_best = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*left_best));
    int32_t *down_best = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*down_best));
    int32_t *up_best = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*up_best));
    int32_t *right = (int32_t *)ARENA_ALLOC(arena, n * sizeof(*right));
    int32_t *left = (int32_t *)ARENA_ALLOC(arena, n * sizeof(*left));
    int32_t *down = (int32_t *)ARENA_ALLOC(arena, n * sizeof(*down));
    int32_t *up = (int32_t *)ARENA_ALLOC(arena, n * sizeof(*up));
    double *right_score = (double *)ARENA_ALLOC(
        arena, n * sizeof(*right_score));
    double *left_score = (double *)ARENA_ALLOC(
        arena, n * sizeof(*left_score));
    double *down_score = (double *)ARENA_ALLOC(
        arena, n * sizeof(*down_score));
    double *up_score = (double *)ARENA_ALLOC(
        arena, n * sizeof(*up_score));
    for (size_t i = 0; i < n; i++) {
        right_best[i] = left_best[i] = -1;
        down_best[i] = up_best[i] = -1;
        right[i] = left[i] = down[i] = up[i] = -1;
        right_score[i] = left_score[i] = DBL_MAX;
        down_score[i] = up_score[i] = DBL_MAX;
    }

    AtgSpatialGrid physical_grid;
    memset(&physical_grid, 0, sizeof physical_grid);
    if (opts->physical_neighbours) {
        double grid_cell = opts->vertical_radius > 0.0
                         ? opts->vertical_radius
                         : opts->horizontal_max_gap;
        if (grid_cell > opts->horizontal_max_gap)
            grid_cell = opts->horizontal_max_gap;
        if (atg_spatial_grid_init(
                arena, set, grid_cell, &physical_grid) != 0)
            return -1;
        int cell_radius = (int)ceil(
            opts->horizontal_max_gap / physical_grid.cell);
        for (size_t row = 0; row < set->nrows; row++) {
            size_t first = row_start[row];
            size_t last = row_start[row + 1];
            atg_spatial_grid_build(
                set, first, last, &physical_grid);
            for (size_t ia = first; ia < last; ia++) {
                const AtlasRibbonTarget *a = &set->target[ia];
                if (!a->accepted) continue;
                int32_t cx = atg_spatial_x(&physical_grid, a->p[0]);
                int32_t cy = atg_spatial_y(&physical_grid, a->p[1]);
                int32_t y0 = cy - cell_radius;
                int32_t y1 = cy + cell_radius;
                int32_t x0 = cx - cell_radius;
                int32_t x1 = cx + cell_radius;
                if (y0 < 0) y0 = 0;
                if (x0 < 0) x0 = 0;
                if (y1 >= physical_grid.ny) y1 = physical_grid.ny - 1;
                if (x1 >= physical_grid.nx) x1 = physical_grid.nx - 1;
                for (int32_t y = y0; y <= y1; y++) {
                    for (int32_t x = x0; x <= x1; x++) {
                        size_t bucket =
                            (size_t)y * (size_t)physical_grid.nx +
                            (size_t)x;
                        for (int32_t ib = physical_grid.head[bucket];
                             ib >= 0;
                             ib = physical_grid.next[(size_t)ib]) {
                            if ((size_t)ib == ia)
                                continue;
                            if (!atg_chart_join_allowed(
                                    set, target_layer, ia, (size_t)ib,
                                    opts, 1)) {
                                out->horizontal_relation_rejects++;
                                continue;
                            }
                            const AtlasRibbonTarget *b =
                                &set->target[(size_t)ib];
                            double score = 0.0;
                            if (!atg_physical_horizontal_compatible(
                                    a, b, set->observation_du,
                                    opts, &score))
                                continue;
                            out->horizontal_candidates++;
                            if (score < right_score[ia]) {
                                right_score[ia] = score;
                                right_best[ia] = ib;
                            }
                            if (score < left_score[(size_t)ib]) {
                                left_score[(size_t)ib] = score;
                                left_best[(size_t)ib] = (int32_t)ia;
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (size_t row = 0; row < set->nrows; row++) {
            size_t first = row_start[row];
            size_t last = row_start[row + 1];
            for (size_t ia = first; ia < last; ia++) {
                const AtlasRibbonTarget *a = &set->target[ia];
                if (!a->accepted) continue;
                for (size_t ib = ia + 1; ib < last; ib++) {
                    const AtlasRibbonTarget *b = &set->target[ib];
                    double gap =
                        ((double)b->column - (double)a->column) *
                        set->observation_du;
                    if (gap > opts->horizontal_max_gap) break;
                    if (!b->accepted || !(gap > 0.0)) continue;
                    if (!atg_chart_join_allowed(
                            set, target_layer, ia, ib, opts, 0)) {
                        out->horizontal_relation_rejects++;
                        continue;
                    }
                    double score = 0.0;
                    if (!atg_horizontal_compatible(
                            a, b, gap, opts, &score))
                        continue;
                    out->horizontal_candidates++;
                    if (score < right_score[ia]) {
                        right_score[ia] = score;
                        right_best[ia] = (int32_t)ib;
                    }
                    if (score < left_score[ib]) {
                        left_score[ib] = score;
                        left_best[ib] = (int32_t)ia;
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < n; i++) {
        int32_t j = right_best[i];
        if (j >= 0 && left_best[j] == (int32_t)i) {
            right[i] = j;
            left[j] = (int32_t)i;
            out->horizontal_edges++;
        }
    }

    if (opts->physical_neighbours) {
        int cell_radius = (int)ceil(
            opts->vertical_radius / physical_grid.cell);
        for (size_t row = 0; row + 1 < set->nrows; row++) {
            size_t first = row_start[row];
            size_t last = row_start[row + 1];
            size_t next_first = row_start[row + 1];
            size_t next_last = row_start[row + 2];
            atg_spatial_grid_build(
                set, next_first, next_last, &physical_grid);
            for (size_t ia = first; ia < last; ia++) {
                const AtlasRibbonTarget *a = &set->target[ia];
                if (!a->accepted) continue;
                int32_t cx = atg_spatial_x(&physical_grid, a->p[0]);
                int32_t cy = atg_spatial_y(&physical_grid, a->p[1]);
                int32_t y0 = cy - cell_radius;
                int32_t y1 = cy + cell_radius;
                int32_t x0 = cx - cell_radius;
                int32_t x1 = cx + cell_radius;
                if (y0 < 0) y0 = 0;
                if (x0 < 0) x0 = 0;
                if (y1 >= physical_grid.ny) y1 = physical_grid.ny - 1;
                if (x1 >= physical_grid.nx) x1 = physical_grid.nx - 1;
                for (int32_t y = y0; y <= y1; y++) {
                    for (int32_t x = x0; x <= x1; x++) {
                        size_t bucket =
                            (size_t)y * (size_t)physical_grid.nx +
                            (size_t)x;
                        for (int32_t ib = physical_grid.head[bucket];
                             ib >= 0;
                             ib = physical_grid.next[(size_t)ib]) {
                            if (!atg_chart_join_allowed(
                                    set, target_layer, ia, (size_t)ib,
                                    opts, 1)) {
                                out->vertical_relation_rejects++;
                                continue;
                            }
                            if (!atg_vertical_winding_compatible(
                                    set, target_layer, ia, (size_t)ib,
                                    opts->vertical_winding_delta_max)) {
                                out->vertical_winding_rejects++;
                                continue;
                            }
                            const AtlasRibbonTarget *b =
                                &set->target[(size_t)ib];
                            double score = 0.0;
                            if (!atg_physical_vertical_compatible(
                                    a, b, opts, &score))
                                continue;
                            out->vertical_candidates++;
                            if (score < down_score[ia]) {
                                down_score[ia] = score;
                                down_best[ia] = ib;
                            }
                            if (score < up_score[(size_t)ib]) {
                                up_score[(size_t)ib] = score;
                                up_best[(size_t)ib] = (int32_t)ia;
                            }
                        }
                    }
                }
            }
        }
    } else {
        int64_t search_columns = (int64_t)ceil(
            opts->vertical_search_u / set->observation_du);
        for (size_t row = 0; row + 1 < set->nrows; row++) {
            size_t first = row_start[row];
            size_t last = row_start[row + 1];
            size_t next_first = row_start[row + 1];
            size_t next_last = row_start[row + 2];
            for (size_t ia = first; ia < last; ia++) {
                const AtlasRibbonTarget *a = &set->target[ia];
                if (!a->accepted) continue;
                int64_t lo64 = (int64_t)a->column - search_columns;
                int64_t hi64 = (int64_t)a->column + search_columns;
                int32_t lo = lo64 < INT32_MIN ? INT32_MIN : (int32_t)lo64;
                int32_t hi = hi64 > INT32_MAX ? INT32_MAX : (int32_t)hi64;
                size_t begin = atg_lower_bound_column(
                    set->target, next_first, next_last, lo);
                for (size_t ib = begin;
                     ib < next_last && set->target[ib].column <= hi;
                     ib++) {
                    const AtlasRibbonTarget *b = &set->target[ib];
                    if (!b->accepted) continue;
                    if (!atg_chart_join_allowed(
                            set, target_layer, ia, ib, opts, 0)) {
                        out->vertical_relation_rejects++;
                        continue;
                    }
                    double du = fabs(
                        ((double)b->column - (double)a->column) *
                        set->observation_du);
                    double score = 0.0;
                    if (!atg_vertical_compatible(a, b, du, opts, &score))
                        continue;
                    out->vertical_candidates++;
                    if (score < down_score[ia]) {
                        down_score[ia] = score;
                        down_best[ia] = (int32_t)ib;
                    }
                    if (score < up_score[ib]) {
                        up_score[ib] = score;
                        up_best[ib] = (int32_t)ia;
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < n; i++) {
        int32_t j = down_best[i];
        if (j >= 0 && up_best[j] == (int32_t)i) {
            down[i] = j;
            up[j] = (int32_t)i;
            out->vertical_edges++;
        }
    }
    if (opts->physical_neighbours &&
        atg_filter_vertical_chain_edges(
            arena, set, opts,
            right, left, down, up, down_best, up_best, out) != 0)
        return -1;

    if (n > (size_t)INT32_MAX / 2) return -1;
    size_t quad_capacity = n * 2;
    AtgQuad *quad = (AtgQuad *)ARENA_ALLOC(
        arena, quad_capacity * sizeof(*quad));
    int32_t *quad_at_a = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*quad_at_a));
    uint8_t *blocked_a = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    uint8_t *horizontal_face = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    uint8_t *vertical_face = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    uint8_t *right_inferred = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    uint8_t *down_inferred = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    for (size_t i = 0; i < n; i++) quad_at_a[i] = -1;
    size_t nquad = 0;
    double core_fit_sum2 = 0.0, core_fit_max = 0.0;
    size_t core_fit_count = 0;
    for (size_t ia = 0; ia < n; ia++) {
        int32_t a = (int32_t)ia;
        int32_t b = right[a];
        int32_t c = down[a];
        if (b < 0 || c < 0) continue;
        int32_t d = down[b];
        if (d < 0 || right[c] != d) continue;
        out->quad_candidates++;
        if (!atg_quad_valid(
                set, set->target, a, b, c, d,
                opts->quad_normal_dot_min)) {
            out->quad_twist_rejects++;
            blocked_a[a] = 1;
            continue;
        }
        AtgQuad *cell = &quad[nquad];
        memset(cell, 0, sizeof *cell);
        cell->a = a;
        cell->b = b;
        cell->c = c;
        cell->d = d;
        cell->observed_edges = 4;
        double fit_rms = 0.0, fit_max = 0.0;
        if (atg_uvxyz_fit4(
                set, target_u, a, b, c, d, &fit_rms, &fit_max) == 0) {
            cell->fit_rms = (float)fit_rms;
            cell->fit_max = (float)fit_max;
            core_fit_sum2 += fit_rms * fit_rms;
            if (fit_max > core_fit_max) core_fit_max = fit_max;
            core_fit_count++;
        }
        quad_at_a[a] = (int32_t)nquad;
        horizontal_face[a] = horizontal_face[c] = 1;
        vertical_face[a] = vertical_face[b] = 1;
        nquad++;
    }
    out->cycle_quads = nquad;
    out->core_fit_rms = core_fit_count > 0
        ? sqrt(core_fit_sum2 / (double)core_fit_count) : 0.0;
    out->core_fit_max = core_fit_max;
    if (nquad == 0) {
        out->selected_row_min = out->selected_row_max = 0;
        out->selected_column_min = out->selected_column_max = 0;
        return 0;
    }

    size_t core_quads = nquad;
    AtgGapProposal *proposal = opts->gap_fill_rounds > 0
        ? (AtgGapProposal *)ARENA_ALLOC(
            arena, n * sizeof(*proposal))
        : NULL;
    double gap_fit_sum2 = 0.0, gap_fit_max = 0.0;
    size_t gap_fit_count = 0;
    for (int round = 0; round < opts->gap_fill_rounds; round++) {
        size_t nproposal = 0;
        uint8_t minimum_observed = round == 0 ? 3 : 2;
        for (size_t ia = 0; ia < n; ia++) {
            int32_t a = (int32_t)ia;
            if (quad_at_a[a] >= 0 || blocked_a[a] ||
                !set->target[a].accepted)
                continue;
            AtgGapProposal best = {0};
            int have_best = 0;
            int32_t b = right[a];
            int32_t c = down[a];
            int32_t d = -1;

            /* A closure made exact by a fill from the preceding round. */
            if (b >= 0 && c >= 0 && (d = down[b]) >= 0 && right[c] == d) {
                uint8_t observed = (uint8_t)(
                    !right_inferred[a] + !right_inferred[c] +
                    !down_inferred[a] + !down_inferred[b]);
                int frontier = horizontal_face[a] || horizontal_face[c] ||
                               vertical_face[a] || vertical_face[b];
                if (frontier && observed >= minimum_observed)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_NONE, observed,
                        &best, &have_best, out);
            }

            /* Missing bottom edge: top + both verticals are trusted. */
            b = right[a];
            c = down[a];
            d = b >= 0 ? down[b] : -1;
            if (b >= 0 && c >= 0 && d >= 0 &&
                right[c] < 0 && left[d] < 0 &&
                (right_best[c] == d || left_best[d] == c)) {
                uint8_t observed = (uint8_t)(
                    !right_inferred[a] + !down_inferred[a] +
                    !down_inferred[b]);
                int frontier = horizontal_face[a] || vertical_face[a] ||
                               vertical_face[b];
                if (frontier && observed >= minimum_observed)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_BOTTOM, observed,
                        &best, &have_best, out);
            }

            /* Missing right edge: both horizontals + left are trusted. */
            b = right[a];
            c = down[a];
            d = c >= 0 ? right[c] : -1;
            if (b >= 0 && c >= 0 && d >= 0 &&
                down[b] < 0 && up[d] < 0 &&
                (down_best[b] == d || up_best[d] == b)) {
                uint8_t observed = (uint8_t)(
                    !right_inferred[a] + !down_inferred[a] +
                    !right_inferred[c]);
                int frontier = horizontal_face[a] || vertical_face[a] ||
                               horizontal_face[c];
                if (frontier && observed >= minimum_observed)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_RIGHT, observed,
                        &best, &have_best, out);
            }

            /* Missing left edge: both horizontals + right are trusted. */
            b = right[a];
            d = b >= 0 ? down[b] : -1;
            c = d >= 0 ? left[d] : -1;
            if (b >= 0 && c >= 0 && d >= 0 &&
                down[a] < 0 && up[c] < 0 &&
                (down_best[a] == c || up_best[c] == a)) {
                uint8_t observed = (uint8_t)(
                    !right_inferred[a] + !down_inferred[b] +
                    !right_inferred[c]);
                int frontier = horizontal_face[a] || vertical_face[b] ||
                               horizontal_face[c];
                if (frontier && observed >= minimum_observed)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_LEFT, observed,
                        &best, &have_best, out);
            }

            /* Missing top edge: bottom + both verticals are trusted. */
            c = down[a];
            d = c >= 0 ? right[c] : -1;
            b = d >= 0 ? up[d] : -1;
            if (b >= 0 && c >= 0 && d >= 0 &&
                right[a] < 0 && left[b] < 0 &&
                (right_best[a] == b || left_best[b] == a)) {
                uint8_t observed = (uint8_t)(
                    !down_inferred[a] + !down_inferred[b] +
                    !right_inferred[c]);
                int frontier = vertical_face[a] || vertical_face[b] ||
                               horizontal_face[c];
                if (frontier && observed >= minimum_observed)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_TOP, observed,
                        &best, &have_best, out);
            }

            /* Missing both verticals: two original horizontal frontier edges
             * bound the cell, and each discarded one-way vertical choice
             * supplies an endpoint pairing.  This is admitted only from the
             * second round because it has two rather than three observed
             * sides; the affine UV->XYZ gate remains decisive. */
            b = right[a];
            c = down_best[a];
            d = b >= 0 ? down_best[b] : -1;
            if (b >= 0 && c >= 0 && d >= 0 && right[c] == d &&
                !right_inferred[a] && !right_inferred[c] &&
                down[a] < 0 && up[c] < 0 &&
                down[b] < 0 && up[d] < 0) {
                int frontier = horizontal_face[a] || horizontal_face[c];
                if (frontier && minimum_observed <= 2)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_VERTICAL_PAIR, 2,
                        &best, &have_best, out);
            }

            /* Missing both horizontals: symmetrically, retain two original
             * vertical frontier edges and require one-way horizontal support
             * at both rows before asking the UV->XYZ fit to close the cell. */
            c = down[a];
            b = right_best[a];
            d = b >= 0 ? down[b] : -1;
            if (b >= 0 && c >= 0 && d >= 0 &&
                (right_best[c] == d || left_best[d] == c) &&
                !down_inferred[a] && !down_inferred[b] &&
                right[a] < 0 && left[b] < 0 &&
                right[c] < 0 && left[d] < 0) {
                int frontier = vertical_face[a] || vertical_face[b];
                if (frontier && minimum_observed <= 2)
                    atg_consider_gap_proposal(
                        set, target_u, opts, a, b, c, d,
                        ATG_GAP_EDGE_HORIZONTAL_PAIR, 2,
                        &best, &have_best, out);
            }
            if (have_best) proposal[nproposal++] = best;
        }
        if (nproposal == 0) break;
        qsort(proposal, nproposal, sizeof(*proposal),
              atg_gap_proposal_compare);
        size_t added = 0;
        for (size_t p = 0; p < nproposal; p++) {
            AtgGapProposal *candidate = &proposal[p];
            int32_t a = candidate->quad.a;
            if (quad_at_a[a] >= 0 ||
                !atg_apply_gap_proposal(
                    candidate, right, left, down, up,
                    right_inferred, down_inferred, out)) {
                out->gap_fill_topology_rejects++;
                continue;
            }
            if (nquad >= quad_capacity || nquad > (size_t)INT32_MAX)
                return -1;
            quad[nquad] = candidate->quad;
            quad_at_a[a] = (int32_t)nquad;
            horizontal_face[candidate->quad.a] = 1;
            horizontal_face[candidate->quad.c] = 1;
            vertical_face[candidate->quad.a] = 1;
            vertical_face[candidate->quad.b] = 1;
            double fit_rms = candidate->quad.fit_rms;
            double fit_maximum = candidate->quad.fit_max;
            gap_fit_sum2 += fit_rms * fit_rms;
            if (fit_maximum > gap_fit_max) gap_fit_max = fit_maximum;
            gap_fit_count++;
            nquad++;
            added++;
        }
        if (added == 0) break;
    }
    out->gap_fill_quads = nquad - core_quads;
    out->gap_fill_fit_rms = gap_fit_count > 0
        ? sqrt(gap_fit_sum2 / (double)gap_fit_count) : 0.0;
    out->gap_fill_fit_max = gap_fit_max;
    if (atg_add_uv_bank_bridges(
            set, target_u, target_layer, opts, row_start,
            right, left, down, up, right_inferred, down_inferred,
            quad, quad_capacity, &nquad, quad_at_a, blocked_a,
            horizontal_face, vertical_face, out) != 0)
        return -1;
    if (nquad > (size_t)INT32_MAX) return -1;

    UnionFind quad_uf = UF_new(arena, (int32_t)nquad);
    int32_t *horizontal_owner = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*horizontal_owner));
    int32_t *vertical_owner = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*vertical_owner));
    for (size_t i = 0; i < n; i++)
        horizontal_owner[i] = vertical_owner[i] = -1;
    for (size_t q = 0; q < nquad; q++) {
        atg_claim_edge(
            horizontal_owner, quad[q].a, (int32_t)q, &quad_uf);
        atg_claim_edge(
            horizontal_owner, quad[q].c, (int32_t)q, &quad_uf);
        atg_claim_edge(
            vertical_owner, quad[q].a, (int32_t)q, &quad_uf);
        atg_claim_edge(
            vertical_owner, quad[q].b, (int32_t)q, &quad_uf);
    }

    size_t *root_quads = (size_t *)ARENA_CALLOC(
        arena, nquad, sizeof(*root_quads));
    size_t *root_core_quads = (size_t *)ARENA_CALLOC(
        arena, nquad, sizeof(*root_core_quads));
    for (size_t q = 0; q < nquad; q++) {
        int32_t root = uf_find(&quad_uf, (int32_t)q);
        root_quads[root]++;
        if (!quad[q].inferred) root_core_quads[root]++;
    }
    for (size_t q = 0; q < nquad; q++) {
        if (root_quads[q] == 0) continue;
        out->ncomponents++;
    }
    out->component = (AtlasTrackGrowComponent *)ARENA_CALLOC(
        arena, out->ncomponents, sizeof(*out->component));
    int32_t *root_component = (int32_t *)ARENA_ALLOC(
        arena, nquad * sizeof(*root_component));
    for (size_t q = 0; q < nquad; q++) root_component[q] = -1;
    size_t component = 0;
    for (size_t q = 0; q < nquad; q++) {
        if (root_quads[q] == 0) continue;
        out->component[component].id = (int32_t)component;
        out->component[component].quads = root_quads[q];
        out->component[component].core_quads = root_core_quads[q];
        out->component[component].inferred_quads =
            root_quads[q] - root_core_quads[q];
        out->component[component].selected =
            root_core_quads[q] >= (size_t)opts->min_component_quads;
        root_component[q] = (int32_t)component;
        if (out->component[component].selected) {
            out->selected_components++;
            out->selected_quads += root_quads[q];
            out->selected_core_quads += root_core_quads[q];
            out->selected_inferred_quads +=
                root_quads[q] - root_core_quads[q];
        }
        component++;
    }
    if (component != out->ncomponents) return -1;

    uint8_t *cycle_node = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    uint8_t *selected_node = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    for (size_t q = 0; q < nquad; q++) {
        const AtgQuad *cell = &quad[q];
        if (!cell->inferred) {
            cycle_node[cell->a] = cycle_node[cell->b] = 1;
            cycle_node[cell->c] = cycle_node[cell->d] = 1;
        }
        int32_t root = uf_find(&quad_uf, (int32_t)q);
        int32_t ci = root_component[root];
        if (ci < 0 || (size_t)ci >= out->ncomponents) return -1;
        if (!out->component[ci].selected) continue;
        selected_node[cell->a] = selected_node[cell->b] = 1;
        selected_node[cell->c] = selected_node[cell->d] = 1;
    }
    for (size_t i = 0; i < n; i++) {
        if (cycle_node[i]) out->cycle_nodes++;
        if (selected_node[i]) out->selected_cycle_nodes++;
    }
    if (out->selected_quads == 0) {
        out->selected_row_min = out->selected_row_max = 0;
        out->selected_column_min = out->selected_column_max = 0;
        return 0;
    }
    if (out->selected_quads > SIZE_MAX / 2) return -1;

    size_t *component_quad_start = (size_t *)ARENA_ALLOC(
        arena, (out->ncomponents + 1) * sizeof(*component_quad_start));
    size_t selected_quad_cursor = 0;
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        component_quad_start[ci] = selected_quad_cursor;
        if (out->component[ci].selected)
            selected_quad_cursor += out->component[ci].quads;
    }
    component_quad_start[out->ncomponents] = selected_quad_cursor;
    if (selected_quad_cursor != out->selected_quads) return -1;
    int32_t *selected_quad = (int32_t *)ARENA_ALLOC(
        arena, out->selected_quads * sizeof(*selected_quad));
    size_t *component_quad_fill = (size_t *)ARENA_ALLOC(
        arena, out->ncomponents * sizeof(*component_quad_fill));
    memcpy(component_quad_fill, component_quad_start,
           out->ncomponents * sizeof(*component_quad_fill));
    for (size_t q = 0; q < nquad; q++) {
        int32_t root = uf_find(&quad_uf, (int32_t)q);
        int32_t ci = root_component[root];
        if (ci < 0 || !out->component[ci].selected) continue;
        size_t at = component_quad_fill[ci]++;
        if (at >= out->selected_quads) return -1;
        selected_quad[at] = (int32_t)q;
    }
    for (size_t ci = 0; ci < out->ncomponents; ci++)
        if (out->component[ci].selected &&
            component_quad_fill[ci] != component_quad_start[ci] +
                                       out->component[ci].quads)
            return -1;

    int32_t *node_stamp = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*node_stamp));
    uint8_t *quadrant_mask = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    int32_t *touched_source = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*touched_source));
    for (size_t i = 0; i < n; i++) node_stamp[i] = -1;
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        AtlasTrackGrowComponent *dst = &out->component[ci];
        if (!dst->selected) continue;
        dst->first_vertex = out->nv;
        dst->first_face = out->nf;
        dst->faces = dst->quads * 2;
        size_t ntouched = 0;
        size_t qend = component_quad_start[ci] + dst->quads;
        for (size_t at = component_quad_start[ci]; at < qend; at++) {
            const AtgQuad *cell = &quad[selected_quad[at]];
            int32_t corner[4] = {cell->a, cell->b, cell->c, cell->d};
            for (int k = 0; k < 4; k++) {
                int32_t source = corner[k];
                if (node_stamp[source] != (int32_t)ci) {
                    node_stamp[source] = (int32_t)ci;
                    quadrant_mask[source] = 0;
                    touched_source[ntouched++] = source;
                }
                quadrant_mask[source] |= (uint8_t)(UINT8_C(1) << k);
            }
        }
        for (size_t i = 0; i < ntouched; i++) {
            dst->vertices++;
            if (atg_mask_needs_split(quadrant_mask[touched_source[i]]))
                dst->vertices++;
        }
        if (out->nv > SIZE_MAX - dst->vertices ||
            out->nf > SIZE_MAX - dst->faces)
            return -1;
        out->nv += dst->vertices;
        out->nf += dst->faces;
    }
    if (out->nv > (size_t)INT32_MAX ||
        out->nf != out->selected_quads * 2)
        return -1;
    out->xyz = (double *)ARENA_ALLOC(
        arena, out->nv * 3 * sizeof(*out->xyz));
    out->local_uv = (double *)ARENA_ALLOC(
        arena, out->nv * 2 * sizeof(*out->local_uv));
    out->global_uv = (double *)ARENA_ALLOC(
        arena, out->nv * 2 * sizeof(*out->global_uv));
    out->uv = (double *)ARENA_ALLOC(
        arena, out->nv * 2 * sizeof(*out->uv));
    out->faces = (int32_t *)ARENA_ALLOC(
        arena, out->nf * 3 * sizeof(*out->faces));
    out->source_target = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*out->source_target));
    out->vertex_component = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*out->vertex_component));
    out->face_component = (int32_t *)ARENA_ALLOC(
        arena, out->nf * sizeof(*out->face_component));
    out->face_inferred = (uint8_t *)ARENA_ALLOC(
        arena, out->nf * sizeof(*out->face_inferred));
    out->vertex_winding = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*out->vertex_winding));
    out->vertex_rank = (uint64_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*out->vertex_rank));

    int32_t *output_right = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*output_right));
    int32_t *output_left = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*output_left));
    int32_t *output_down = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*output_down));
    int32_t *output_up = (int32_t *)ARENA_ALLOC(
        arena, out->nv * sizeof(*output_up));
    for (size_t i = 0; i < out->nv; i++)
        output_right[i] = output_left[i] =
        output_down[i] = output_up[i] = -1;
    int32_t *source_to_output = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*source_to_output));
    int32_t *source_to_output_second = (int32_t *)ARENA_ALLOC(
        arena, n * sizeof(*source_to_output_second));
    for (size_t i = 0; i < n; i++) node_stamp[i] = -1;

    int32_t max_chart = -1;
    size_t vertex_cursor = 0, face_cursor = 0;
    for (size_t ci = 0; ci < out->ncomponents; ci++) {
        AtlasTrackGrowComponent *dst = &out->component[ci];
        if (!dst->selected) continue;
        if (vertex_cursor != dst->first_vertex ||
            face_cursor != dst->first_face)
            return -1;
        size_t qend = component_quad_start[ci] + dst->quads;
        size_t ntouched = 0;
        for (size_t at = component_quad_start[ci]; at < qend; at++) {
            const AtgQuad *cell = &quad[selected_quad[at]];
            int32_t source[4] = {cell->a, cell->b, cell->c, cell->d};
            for (int k = 0; k < 4; k++) {
                int32_t si = source[k];
                if (node_stamp[si] != (int32_t)ci) {
                    node_stamp[si] = (int32_t)ci;
                    quadrant_mask[si] = 0;
                    touched_source[ntouched++] = si;
                }
                quadrant_mask[si] |= (uint8_t)(UINT8_C(1) << k);
            }
        }
        for (size_t i = 0; i < ntouched; i++) {
            int32_t si = touched_source[i];
            int copies = atg_mask_needs_split(quadrant_mask[si]) ? 2 : 1;
            for (int copy = 0; copy < copies; copy++) {
                if (vertex_cursor > (size_t)INT32_MAX) return -1;
                int32_t output = (int32_t)vertex_cursor;
                if (copy == 0)
                    source_to_output[si] = output;
                else
                    source_to_output_second[si] = output;
                out->source_target[vertex_cursor] = si;
                out->vertex_component[vertex_cursor] = (int32_t)ci;
                int32_t layer = target_layer[si];
                if (layer < 0 || (size_t)layer >= set->nlayer_samples)
                    return -1;
                out->vertex_winding[vertex_cursor] =
                    set->layer_sample[layer].winding;
                out->vertex_rank[vertex_cursor] =
                    set->layer_sample[layer].rank;
                if (out->vertex_winding[vertex_cursor] < out->winding_min)
                    out->winding_min = out->vertex_winding[vertex_cursor];
                if (out->vertex_winding[vertex_cursor] > out->winding_max)
                    out->winding_max = out->vertex_winding[vertex_cursor];
                atg_world_point(
                    set, &set->target[si],
                    &out->xyz[vertex_cursor * 3]);
                const AtlasRibbonTarget *t = &set->target[si];
                if (t->row < out->selected_row_min)
                    out->selected_row_min = t->row;
                if (t->row > out->selected_row_max)
                    out->selected_row_max = t->row;
                if (t->column < out->selected_column_min)
                    out->selected_column_min = t->column;
                    if (t->column > out->selected_column_max)
                        out->selected_column_max = t->column;
                    if (t->chart0 > max_chart) max_chart = t->chart0;
                    if (t->chart1 > max_chart) max_chart = t->chart1;
                    vertex_cursor++;
            }
            if (copies == 1)
                source_to_output_second[si] = source_to_output[si];
        }
        for (size_t at = component_quad_start[ci]; at < qend; at++) {
            const AtgQuad *cell = &quad[selected_quad[at]];
            int32_t source[4] = {cell->a, cell->b, cell->c, cell->d};
            int32_t output[4];
            for (int k = 0; k < 4; k++) {
                int32_t si = source[k];
                uint8_t mask = quadrant_mask[si];
                int second_fan =
                    (mask == UINT8_C(0x09) && k == 3) ||
                    (mask == UINT8_C(0x06) && k == 2);
                output[k] = second_fan
                    ? source_to_output_second[si]
                    : source_to_output[si];
                if (output[k] < 0) return -1;
            }
            int32_t a = output[0], b = output[1];
            int32_t c = output[2], d = output[3];
            if (atg_set_edge(output_right, output_left, a, b) != 0 ||
                atg_set_edge(output_right, output_left, c, d) != 0 ||
                atg_set_edge(output_down, output_up, a, c) != 0 ||
                atg_set_edge(output_down, output_up, b, d) != 0)
                return -1;
            out->faces[face_cursor * 3] = a;
            out->faces[face_cursor * 3 + 1] = b;
            out->faces[face_cursor * 3 + 2] = c;
            out->face_component[face_cursor] = (int32_t)ci;
            out->face_inferred[face_cursor++] = cell->inferred;
            out->faces[face_cursor * 3] = b;
            out->faces[face_cursor * 3 + 1] = d;
            out->faces[face_cursor * 3 + 2] = c;
            out->face_component[face_cursor] = (int32_t)ci;
            out->face_inferred[face_cursor++] = cell->inferred;
        }
        if (vertex_cursor != dst->first_vertex + dst->vertices ||
            face_cursor != dst->first_face + dst->faces)
            return -1;
    }
    if (vertex_cursor != out->nv || face_cursor != out->nf) return -1;
    out->selected_rows =
        (size_t)((int64_t)out->selected_row_max -
                 (int64_t)out->selected_row_min + 1);
    out->source_u_span =
        ((double)out->selected_column_max -
         (double)out->selected_column_min) *
        set->observation_du;

    if (max_chart >= 0) {
        uint8_t *chart_seen = (uint8_t *)ARENA_CALLOC(
            arena, (size_t)max_chart + 1, 1);
        for (size_t i = 0; i < out->nv; i++) {
            const AtlasRibbonTarget *target =
                &set->target[out->source_target[i]];
            int32_t chart[2] = {target->chart0, target->chart1};
            for (int k = 0; k < 2; k++) {
                if (chart[k] >= 0 && !chart_seen[chart[k]]) {
                    chart_seen[chart[k]] = 1;
                    out->selected_charts++;
                }
            }
        }
    }

    if (atg_build_parameterization(
            arena, set, calibration, opts, right, down,
            output_right, output_left, output_down, output_up, out) != 0)
        return -1;
    for (size_t f = 0; f < out->nf; f++)
        out->surface_area += atg_triangle_area(
            out->xyz,
            out->faces[f * 3],
            out->faces[f * 3 + 1],
            out->faces[f * 3 + 2]);
    return 0;
}

static void atg_selftest_set(
    AtlasRibbonObservationSet *set,
    AtlasRibbonTarget target[16],
    AtlasRibbonLayerSample layer[16])
{
    memset(set, 0, sizeof *set);
    memset(target, 0, 16 * sizeof(*target));
    memset(layer, 0, 16 * sizeof(*layer));
    set->target = target;
    set->ntarget = 16;
    set->layer_sample = layer;
    set->nlayer_samples = 16;
    set->nrows = 4;
    set->observation_du = 2.0;
    set->axis[0] = 1.0;
    set->basis0[1] = 1.0;
    set->basis1[2] = 1.0;
    for (int row = 0; row < 4; row++)
        for (int column = 0; column < 4; column++) {
            AtlasRibbonTarget *t = &target[row * 4 + column];
            t->row = row;
            t->column = column;
            t->accepted = 1;
            t->chart0 = column;
            t->u = 2.0 * (double)column;
            t->v = 4.0 * (double)row;
            t->p[0] = 2.0 * (double)column;
            t->p[1] = 0.0;
            t->tangent[0] = 1.0;
            layer[row * 4 + column].row = row;
            layer[row * 4 + column].source_column = column;
            layer[row * 4 + column].chart = column;
            layer[row * 4 + column].rank = (uint64_t)column;
            layer[row * 4 + column].winding = 0;
            layer[row * 4 + column].p[0] = t->p[0];
            layer[row * 4 + column].tangent[0] = 1.0;
        }
}

static void atg_selftest_point_contact_set(
    AtlasRibbonObservationSet *set,
    AtlasRibbonTarget target[7],
    AtlasRibbonLayerSample layer[7])
{
    static const int32_t coordinate[7][2] = {
        {0, 0}, {0, 1},
        {1, 0}, {1, 1}, {1, 2},
        {2, 1}, {2, 2}
    };
    memset(set, 0, sizeof *set);
    memset(target, 0, 7 * sizeof(*target));
    memset(layer, 0, 7 * sizeof(*layer));
    set->target = target;
    set->ntarget = 7;
    set->layer_sample = layer;
    set->nlayer_samples = 7;
    set->nrows = 3;
    set->observation_du = 2.0;
    set->axis[0] = 1.0;
    set->basis0[1] = 1.0;
    set->basis1[2] = 1.0;
    for (int i = 0; i < 7; i++) {
        AtlasRibbonTarget *t = &target[i];
        t->row = coordinate[i][0];
        t->column = coordinate[i][1];
        t->accepted = 1;
        t->chart0 = i < 4 ? 0 : 1;
        t->u = 2.0 * (double)t->column;
        t->v = 4.0 * (double)t->row;
        t->p[0] = 2.0 * (double)t->column;
        t->tangent[0] = 1.0;
        layer[i].row = t->row;
        layer[i].source_column = t->column;
        layer[i].chart = t->chart0;
        layer[i].rank = (uint64_t)t->chart0;
        layer[i].winding = 0;
        layer[i].p[0] = t->p[0];
        layer[i].tangent[0] = 1.0;
    }
}

static void atg_selftest_gap_fill_set(
    AtlasRibbonObservationSet *set,
    AtlasRibbonTarget target[8],
    AtlasRibbonLayerSample layer[8])
{
    static const int32_t coordinate[8][2] = {
        {0, 0}, {0, 2},
        {1, 0}, {1, 2},
        {2, 0}, {2, 1}, {2, 2}, {2, 3}
    };
    /* On the final row c->d is the correct surface edge, but d's one-way
     * left choice is the interloper e and e's right choice is f.  Mutual-best
     * therefore drops c->d.  The trusted cell above supplies a frontier, and
     * the four intended corners lie exactly on one affine UV->XYZ patch. */
    static const double p0[8] = {
        0.0, 3.0,
        0.0, 3.0,
        0.0, 0.5, 3.0, 4.5
    };
    memset(set, 0, sizeof *set);
    memset(target, 0, 8 * sizeof(*target));
    memset(layer, 0, 8 * sizeof(*layer));
    set->target = target;
    set->ntarget = 8;
    set->layer_sample = layer;
    set->nlayer_samples = 8;
    set->nrows = 3;
    set->observation_du = 2.0;
    set->axis[0] = 1.0;
    set->basis0[1] = 1.0;
    set->basis1[2] = 1.0;
    for (int i = 0; i < 8; i++) {
        AtlasRibbonTarget *t = &target[i];
        t->row = coordinate[i][0];
        t->column = coordinate[i][1];
        t->accepted = 1;
        t->chart0 = i;
        t->u = 2.0 * (double)t->column;
        t->v = 4.0 * (double)t->row;
        t->p[0] = p0[i];
        t->tangent[0] = 1.0;
        layer[i].row = t->row;
        layer[i].source_column = t->column;
        layer[i].chart = i;
        layer[i].rank = (uint64_t)t->column;
        layer[i].winding = t->column;
        layer[i].p[0] = t->p[0];
        layer[i].tangent[0] = 1.0;
    }
}

static void atg_selftest_uv_bank_set(
    AtlasRibbonObservationSet *set,
    AtlasRibbonTarget target[12],
    AtlasRibbonLayerSample layer[12],
    double right_offset)
{
    memset(set, 0, sizeof *set);
    memset(target, 0, 12 * sizeof(*target));
    memset(layer, 0, 12 * sizeof(*layer));
    set->target = target;
    set->ntarget = 12;
    set->layer_sample = layer;
    set->nlayer_samples = 12;
    set->nrows = 3;
    set->observation_du = 2.0;
    set->axis[0] = 1.0;
    set->basis0[1] = 1.0;
    set->basis1[2] = 1.0;
    int at = 0;
    for (int row = 0; row < 3; row++) {
        int column[4] = {0, 1, 4, 5};
        for (int j = 0; j < 4; j++, at++) {
            AtlasRibbonTarget *t = &target[at];
            t->row = row;
            t->column = column[j];
            t->accepted = 1;
            t->chart0 = j < 2 ? 0 : 1;
            t->u = 2.0 * (double)t->column;
            t->v = 4.0 * (double)row;
            t->p[0] = t->u + (j < 2 ? 0.0 : right_offset);
            t->tangent[0] = 1.0;
            layer[at].row = t->row;
            layer[at].source_column = t->column;
            layer[at].chart = t->chart0;
            layer[at].rank = (uint64_t)t->chart0;
            layer[at].winding = 0;
            layer[at].p[0] = t->p[0];
            layer[at].tangent[0] = 1.0;
        }
    }
}

int AtlasTrackGrow_selftest(void)
{
    int failures = 0;
    {
        int32_t right[4] = {1, -1, 3, -1};
        int32_t left[4] = {-1, 0, -1, 2};
        int32_t down[4] = {-1, -1, -1, -1};
        int32_t up[4] = {-1, -1, -1, -1};
        uint8_t right_inferred[4] = {0};
        uint8_t down_inferred[4] = {0};
        AtgGapProposal proposal = {0};
        AtlasTrackGrowResult pair_result = {0};
        proposal.quad.a = 0;
        proposal.quad.b = 1;
        proposal.quad.c = 2;
        proposal.quad.d = 3;
        proposal.missing_edge = ATG_GAP_EDGE_VERTICAL_PAIR;
        if (!atg_apply_gap_proposal(
                &proposal, right, left, down, up,
                right_inferred, down_inferred, &pair_result) ||
            down[0] != 2 || down[1] != 3 ||
            pair_result.gap_fill_pair_quads != 1 ||
            pair_result.gap_fill_vertical_edges != 2) {
            fprintf(stderr,
                    "[atlas_track_grow selftest] FAIL vertical pair fill\n");
            failures++;
        }
    }
    {
        int32_t right[4] = {-1, -1, -1, -1};
        int32_t left[4] = {-1, -1, -1, -1};
        int32_t down[4] = {2, 3, -1, -1};
        int32_t up[4] = {-1, -1, 0, 1};
        uint8_t right_inferred[4] = {0};
        uint8_t down_inferred[4] = {0};
        AtgGapProposal proposal = {0};
        AtlasTrackGrowResult pair_result = {0};
        proposal.quad.a = 0;
        proposal.quad.b = 1;
        proposal.quad.c = 2;
        proposal.quad.d = 3;
        proposal.missing_edge = ATG_GAP_EDGE_HORIZONTAL_PAIR;
        if (!atg_apply_gap_proposal(
                &proposal, right, left, down, up,
                right_inferred, down_inferred, &pair_result) ||
            right[0] != 1 || right[2] != 3 ||
            pair_result.gap_fill_pair_quads != 1 ||
            pair_result.gap_fill_horizontal_edges != 2) {
            fprintf(stderr,
                    "[atlas_track_grow selftest] FAIL horizontal pair fill\n");
            failures++;
        }
    }
    Arena_T arena = Arena_new();
    if (arena == NULL) return 1;
    AtlasRibbonObservationSet set;
    AtlasRibbonTarget target[16];
    AtlasRibbonLayerSample layer[16];
    atg_selftest_set(&set, target, layer);
    AtlasTrackGrowOptions opts;
    AtlasTrackGrowOptions_default(&opts);
    opts.min_component_quads = 1;
    opts.parameter_stay_weight = 0.0;
    ScaffoldCalib calibration;
    Scaffold_calib_default(&calibration);
    AtlasTrackGrowResult result;
    if (AtlasTrackGrow_build(
            arena, &set, &calibration, &opts, &result) != 0 ||
        result.horizontal_edges != 12 ||
        result.vertical_edges != 12 ||
        result.cycle_quads != 9 ||
        result.ncomponents != 1 ||
        result.selected_components != 1 ||
        result.selected_quads != 9 ||
        result.nv != 16 ||
        result.nf != 18 ||
        result.parameter_islands != 1 ||
        result.horizontal_residual_rms > 1.0e-6 ||
        result.vertical_residual_rms > 1.0e-6) {
        fprintf(stderr,
                "[atlas_track_grow selftest] FAIL regular lattice: "
                "h=%zu v=%zu q=%zu comp=%zu selected=%zu nv=%zu nf=%zu "
                "hrms=%.6g vrms=%.6g\n",
                result.horizontal_edges, result.vertical_edges,
                result.cycle_quads, result.ncomponents,
                result.selected_quads, result.nv, result.nf,
                result.horizontal_residual_rms,
                result.vertical_residual_rms);
        failures++;
    }
    Arena_dispose(&arena);

    AtlasTrackGrowOptions_default(&opts);
    opts.min_component_quads = 1;
    opts.parameter_stay_weight = 0.0;
    arena = Arena_new();
    if (arena == NULL) return failures + 1;
    AtlasRibbonTarget bank_target[12];
    AtlasRibbonLayerSample bank_layer[12];
    atg_selftest_uv_bank_set(&set, bank_target, bank_layer, 0.0);
    AtlasTrackGrowOptions_default(&opts);
    opts.min_component_quads = 1;
    opts.gap_fill_rounds = 0;
    opts.horizontal_max_gap = 3.0;
    opts.uv_bridge_max_u = 8.0;
    opts.parameter_stay_weight = 0.0;
    if (AtlasTrackGrow_build(
            arena, &set, &calibration, &opts, &result) != 0 ||
        result.cycle_quads != 4 ||
        result.uv_bridge_quads != 2 ||
        result.uv_bridge_horizontal_edges != 3 ||
        result.selected_core_quads != 4 ||
        result.selected_inferred_quads != 2 ||
        result.selected_quads != 6 ||
        result.ncomponents != 1 ||
        result.nv != 12 || result.nf != 12 ||
        result.uv_bridge_fit_rms > 1.0e-6) {
        fprintf(stderr,
                "[atlas_track_grow selftest] FAIL UV bank bridge: "
                "core=%zu bank=%zu edges=%zu selected=%zu+%zu "
                "comp=%zu nv=%zu nf=%zu fit=%.6g\n",
                result.cycle_quads, result.uv_bridge_quads,
                result.uv_bridge_horizontal_edges,
                result.selected_core_quads,
                result.selected_inferred_quads,
                result.ncomponents, result.nv, result.nf,
                result.uv_bridge_fit_rms);
        failures++;
    }
    Arena_dispose(&arena);

    arena = Arena_new();
    if (arena == NULL) return failures + 1;
    atg_selftest_uv_bank_set(&set, bank_target, bank_layer, -2.0);
    opts.uv_bridge_fit_rms_max = 1.0e-4;
    opts.uv_bridge_fit_max_max = 1.0e-3;
    if (AtlasTrackGrow_build(
            arena, &set, &calibration, &opts, &result) != 0 ||
        result.cycle_quads != 4 ||
        result.uv_bridge_quads != 0 ||
        result.uv_bridge_fit_rejects == 0 ||
        result.ncomponents != 2 || result.selected_quads != 4) {
        fprintf(stderr,
                "[atlas_track_grow selftest] FAIL UV bank fit reject: "
                "core=%zu bank=%zu fit_reject=%zu comp=%zu selected=%zu\n",
                result.cycle_quads, result.uv_bridge_quads,
                result.uv_bridge_fit_rejects,
                result.ncomponents, result.selected_quads);
        failures++;
    }
    Arena_dispose(&arena);

    AtlasTrackGrowOptions_default(&opts);
    opts.min_component_quads = 1;
    opts.parameter_stay_weight = 0.0;
    arena = Arena_new();
    if (arena == NULL) return failures + 1;
    AtlasRibbonTarget point_target[7];
    AtlasRibbonLayerSample point_layer[7];
    atg_selftest_point_contact_set(&set, point_target, point_layer);
    if (AtlasTrackGrow_build(
            arena, &set, &calibration, &opts, &result) != 0 ||
        result.cycle_quads != 2 ||
        result.ncomponents != 2 ||
        result.selected_components != 2 ||
        result.selected_quads != 2 ||
        result.selected_cycle_nodes != 7 ||
        result.nv != 8 ||
        result.nf != 4 ||
        result.nv - result.selected_cycle_nodes != 1 ||
        result.parameter_point_bridges < 1 ||
        result.parameter_islands != 1) {
        fprintf(stderr,
                "[atlas_track_grow selftest] FAIL point contact: "
                "q=%zu comp=%zu selected_comp=%zu selected_q=%zu "
                "nodes=%zu nv=%zu nf=%zu point_bridges=%zu islands=%zu\n",
                result.cycle_quads, result.ncomponents,
                result.selected_components, result.selected_quads,
                result.selected_cycle_nodes, result.nv, result.nf,
                result.parameter_point_bridges, result.parameter_islands);
        failures++;
    }
    Arena_dispose(&arena);

    arena = Arena_new();
    if (arena == NULL) return failures + 1;
    AtlasRibbonTarget gap_target[8];
    AtlasRibbonLayerSample gap_layer[8];
    atg_selftest_gap_fill_set(&set, gap_target, gap_layer);
    if (AtlasTrackGrow_build(
            arena, &set, &calibration, &opts, &result) != 0) {
        fprintf(stderr,
                "[atlas_track_grow selftest] FAIL UV->XYZ gap build\n");
        failures++;
    } else {
        size_t inferred_faces = 0;
        for (size_t f = 0; f < result.nf; f++)
            if (result.face_inferred[f]) inferred_faces++;
        if (result.cycle_quads != 1 ||
            result.gap_fill_quads != 1 ||
            result.gap_fill_horizontal_edges != 1 ||
            result.selected_core_quads != 1 ||
            result.selected_inferred_quads != 1 ||
            result.selected_quads != 2 ||
            result.ncomponents != 1 ||
            result.selected_components != 1 ||
            result.nv != 6 || result.nf != 4 ||
            inferred_faces != 2 ||
            result.gap_fill_fit_rms > 1.0e-5) {
            fprintf(stderr,
                    "[atlas_track_grow selftest] FAIL UV->XYZ gap: "
                    "core=%zu fill=%zu fill_h=%zu selected=%zu+%zu "
                    "comp=%zu nv=%zu nf=%zu inferred_faces=%zu fit=%.6g\n",
                    result.cycle_quads, result.gap_fill_quads,
                    result.gap_fill_horizontal_edges,
                    result.selected_core_quads,
                    result.selected_inferred_quads,
                    result.ncomponents, result.nv, result.nf,
                    inferred_faces, result.gap_fill_fit_rms);
            failures++;
        }
    }
    Arena_dispose(&arena);

    arena = Arena_new();
    if (arena == NULL) return failures + 1;
    memset(&set, 0, sizeof set);
    set.observation_du = 2.0;
    if (AtlasTrackGrow_build(
            arena, &set, &calibration, &opts, &result) != 0 ||
        result.nv != 0 || result.nf != 0) {
        fprintf(stderr,
                "[atlas_track_grow selftest] FAIL empty evidence\n");
        failures++;
    }
    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_track_grow selftest] %s (%d failures)\n",
            failures == 0 ? "ALL PASS" : "FAILURES", failures);
    return failures;
}
