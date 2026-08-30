#include "atlas_developable_complete.h"

#include "../common/union_find.h"

#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ADC_TAU 6.283185307179586476925286766559

typedef struct {
    int32_t row;
    int32_t column;
    int32_t target;
    int32_t vertex;
    int32_t component;
} AdcSelected;

typedef struct {
    int direction;
    int32_t row;
    int32_t column_left;
    int32_t column_right;
    int32_t target_left;
    int32_t target_right;
    int32_t vertex_left;
    int32_t vertex_right;
    int32_t component_left;
    int32_t component_right;
    double gap_u;
    double chord;
    double score;
} AdcRowGap;

typedef struct {
    size_t first;
    size_t count;
    int32_t component_left;
    int32_t component_right;
} AdcRun;

typedef struct {
    double primal;
    double dual;
    double tail_before;
    double tail_after;
    int pcg_iterations;
    int iterations;
} AdcAdmmStats;

typedef struct {
    uint64_t key;
    double xyz_sum[3];
    double weight;
    int32_t component;
    int32_t winding;
    uint64_t rank;
    int32_t vertex;
    uint8_t conflict;
} AdcLatticeNode;

typedef struct {
    Arena_T arena;
    AdcLatticeNode *node;
    size_t node_capacity;
    size_t node_count;
    uint8_t *cell;
    size_t rows;
    size_t columns;
    size_t cell_count;
    size_t conflicts;
    size_t rejected_cells;
    double merge_tolerance;
} AdcLattice;

typedef struct {
    int32_t *head;
    int32_t *next;
    size_t rows;
    size_t columns;
} AdcSeedGrid;

typedef struct {
    int32_t row;
    int32_t column;
    int32_t vertex[4];
    int direction;
    double score;
} AdcDirectCell;

typedef struct {
    int32_t gap_a;
    int32_t gap_b;
    int32_t span;
    int direction;
    double score;
} AdcShortStrip;

static uint64_t adc_hash64(uint64_t value)
{
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return value;
}

static uint64_t adc_lattice_key(size_t row, size_t column)
{
    return ((uint64_t)(uint32_t)row << 32) | (uint64_t)(uint32_t)column;
}

static int adc_lattice_rehash(AdcLattice *lattice, size_t capacity)
{
    if (lattice == NULL || capacity < 16 ||
        (capacity & (capacity - 1)) != 0)
        return -1;
    AdcLatticeNode *next = (AdcLatticeNode *)ARENA_ALLOC(
        lattice->arena, capacity * sizeof(*next));
    for (size_t i = 0; i < capacity; i++) {
        memset(&next[i], 0, sizeof(next[i]));
        next[i].key = UINT64_MAX;
        next[i].vertex = -1;
    }
    for (size_t i = 0; i < lattice->node_capacity; i++) {
        if (lattice->node[i].key == UINT64_MAX) continue;
        size_t at = (size_t)adc_hash64(lattice->node[i].key) &
                    (capacity - 1);
        while (next[at].key != UINT64_MAX)
            at = (at + 1) & (capacity - 1);
        next[at] = lattice->node[i];
    }
    lattice->node = next;
    lattice->node_capacity = capacity;
    return 0;
}

static int adc_lattice_reserve_additions(
    AdcLattice *lattice, size_t additions)
{
    if (lattice == NULL || additions == 0) return 0;
    if (additions > SIZE_MAX - lattice->node_count) return -1;
    size_t desired = lattice->node_count + additions;
    size_t capacity = lattice->node_capacity;
    if (capacity == 0) capacity = 4096;
    /*
     * Grow before a patch saves its temporary arena mark.  A rehash performed
     * while adc_lattice_add_patch() is consuming temporary XYZ samples would
     * otherwise be rolled back with those samples, leaving lattice->node
     * pointing at reclaimed arena storage.
     */
    while (desired >= (capacity / 10) * 7) {
        if (capacity > SIZE_MAX / 2) return -1;
        capacity *= 2;
    }
    if (capacity == lattice->node_capacity) return 0;
    return adc_lattice_rehash(lattice, capacity);
}

static AdcLatticeNode *adc_lattice_find(
    AdcLattice *lattice, uint64_t key, int create)
{
    if (lattice == NULL) return NULL;
    if (lattice->node_capacity == 0) {
        if (!create || adc_lattice_rehash(lattice, 4096) != 0) return NULL;
    }
    if (create &&
        (lattice->node_count + 1) * 10 >= lattice->node_capacity * 7) {
        if (lattice->node_capacity > SIZE_MAX / 2 ||
            adc_lattice_rehash(
                lattice, lattice->node_capacity * 2) != 0)
            return NULL;
    }
    size_t at = (size_t)adc_hash64(key) & (lattice->node_capacity - 1);
    while (lattice->node[at].key != UINT64_MAX &&
           lattice->node[at].key != key)
        at = (at + 1) & (lattice->node_capacity - 1);
    if (lattice->node[at].key == UINT64_MAX) {
        if (!create) return NULL;
        memset(&lattice->node[at], 0, sizeof(lattice->node[at]));
        lattice->node[at].key = key;
        lattice->node[at].vertex = -1;
        lattice->node_count++;
    }
    return &lattice->node[at];
}

static int adc_lattice_init(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasDevelopableCompleteOptions *opts,
    AdcLattice *lattice)
{
    if (arena == NULL || evidence == NULL || opts == NULL || lattice == NULL)
        return -1;
    memset(lattice, 0, sizeof *lattice);
    lattice->arena = arena;
    lattice->rows = evidence->nrows;
    lattice->merge_tolerance = opts->bank_fit_max_max;
    int32_t maximum_column = -1;
    for (size_t i = 0; i < evidence->ntarget; i++) {
        if (evidence->target[i].column > maximum_column)
            maximum_column = evidence->target[i].column;
    }
    if (lattice->rows == 0 || maximum_column < 0) return 0;
    lattice->columns = (size_t)maximum_column + 1;
    if (lattice->rows > SIZE_MAX / lattice->columns) return -1;
    lattice->cell = (uint8_t *)ARENA_CALLOC(
        arena, lattice->rows * lattice->columns, sizeof(*lattice->cell));
    return adc_lattice_rehash(lattice, 4096);
}

static int adc_lattice_add_node(
    AdcLattice *lattice,
    int32_t row,
    int32_t column,
    const double xyz[3],
    int32_t component,
    int32_t winding,
    uint64_t rank)
{
    if (lattice == NULL || row < 0 || column < 0 ||
        (size_t)row >= lattice->rows ||
        (size_t)column >= lattice->columns)
        return -1;
    AdcLatticeNode *node = adc_lattice_find(
        lattice, adc_lattice_key((size_t)row, (size_t)column), 1);
    if (node == NULL) return -1;
    if (node->conflict) return 0;
    if (node->weight > 0.0) {
        double distance2 = 0.0;
        for (int d = 0; d < 3; d++) {
            double delta = xyz[d] - node->xyz_sum[d] / node->weight;
            distance2 += delta * delta;
        }
        if (sqrt(distance2) > lattice->merge_tolerance) {
            node->conflict = 1;
            lattice->conflicts++;
            return 0;
        }
    } else {
        node->component = component;
        node->winding = winding;
        node->rank = rank;
    }
    for (int d = 0; d < 3; d++) node->xyz_sum[d] += xyz[d];
    node->weight += 1.0;
    return 0;
}

static void adc_lattice_mark_cell(
    AdcLattice *lattice, int32_t row, int32_t column)
{
    if (lattice == NULL || row < 0 || column < 0 ||
        (size_t)row + 1 >= lattice->rows ||
        (size_t)column + 1 >= lattice->columns)
        return;
    size_t index = (size_t)row * lattice->columns + (size_t)column;
    if (!lattice->cell[index]) {
        lattice->cell[index] = 1;
        lattice->cell_count++;
    }
}

static double adc_clamp(double value, double low, double high)
{
    if (value < low) return low;
    if (value > high) return high;
    return value;
}

static double adc_dot2(const double a[2], const double b[2])
{
    return a[0] * b[0] + a[1] * b[1];
}

static double adc_distance2(const double a[2], const double b[2])
{
    double x = b[0] - a[0];
    double y = b[1] - a[1];
    return sqrt(x * x + y * y);
}

static double adc_distance3(const double *a, const double *b)
{
    double x = b[0] - a[0];
    double y = b[1] - a[1];
    double z = b[2] - a[2];
    return sqrt(x * x + y * y + z * z);
}

static double adc_dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void adc_cross3(
    const double a[3], const double b[3], double output[3])
{
    output[0] = a[1] * b[2] - a[2] * b[1];
    output[1] = a[2] * b[0] - a[0] * b[2];
    output[2] = a[0] * b[1] - a[1] * b[0];
}

static int adc_normalize3(double value[3])
{
    double length = sqrt(adc_dot3(value, value));
    if (!(length > 1.0e-12) || !isfinite(length)) return -1;
    for (int d = 0; d < 3; d++) value[d] /= length;
    return 0;
}

static int adc_compare_selected(const void *left, const void *right)
{
    const AdcSelected *a = (const AdcSelected *)left;
    const AdcSelected *b = (const AdcSelected *)right;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->column != b->column) return a->column < b->column ? -1 : 1;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->target != b->target) return a->target < b->target ? -1 : 1;
    return 0;
}

static int adc_compare_selected_column(const void *left, const void *right)
{
    const AdcSelected *a = (const AdcSelected *)left;
    const AdcSelected *b = (const AdcSelected *)right;
    if (a->column != b->column)
        return a->column < b->column ? -1 : 1;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->target != b->target) return a->target < b->target ? -1 : 1;
    return 0;
}

static int adc_compare_gap(const void *left, const void *right)
{
    const AdcRowGap *a = (const AdcRowGap *)left;
    const AdcRowGap *b = (const AdcRowGap *)right;
    if (a->direction != b->direction)
        return a->direction < b->direction ? -1 : 1;
    if (a->component_left != b->component_left)
        return a->component_left < b->component_left ? -1 : 1;
    if (a->component_right != b->component_right)
        return a->component_right < b->component_right ? -1 : 1;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->score < b->score) return -1;
    if (a->score > b->score) return 1;
    if (a->column_left != b->column_left)
        return a->column_left < b->column_left ? -1 : 1;
    if (a->column_right != b->column_right)
        return a->column_right < b->column_right ? -1 : 1;
    return 0;
}

static int adc_compare_gap_spatial(const void *left, const void *right)
{
    const AdcRowGap *a = (const AdcRowGap *)left;
    const AdcRowGap *b = (const AdcRowGap *)right;
    if (a->direction != b->direction)
        return a->direction < b->direction ? -1 : 1;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->column_left != b->column_left)
        return a->column_left < b->column_left ? -1 : 1;
    if (a->column_right != b->column_right)
        return a->column_right < b->column_right ? -1 : 1;
    if (a->score < b->score) return -1;
    if (a->score > b->score) return 1;
    if (a->component_left != b->component_left)
        return a->component_left < b->component_left ? -1 : 1;
    if (a->component_right != b->component_right)
        return a->component_right < b->component_right ? -1 : 1;
    return 0;
}

static int adc_compare_short_strip(const void *left, const void *right)
{
    const AdcShortStrip *a = (const AdcShortStrip *)left;
    const AdcShortStrip *b = (const AdcShortStrip *)right;
    if (a->span != b->span) return a->span < b->span ? -1 : 1;
    if (a->score < b->score) return -1;
    if (a->score > b->score) return 1;
    if (a->direction != b->direction)
        return a->direction < b->direction ? -1 : 1;
    if (a->gap_a != b->gap_a) return a->gap_a < b->gap_a ? -1 : 1;
    if (a->gap_b != b->gap_b) return a->gap_b < b->gap_b ? -1 : 1;
    return 0;
}

void AtlasDevelopableCompleteOptions_default(
    AtlasDevelopableCompleteOptions *opts)
{
    if (opts == NULL) return;
    memset(opts, 0, sizeof *opts);
    opts->min_gap_u = 4.0;
    opts->max_gap_u = 96.0;
    opts->min_gap_v = 4.0;
    opts->max_gap_v = 96.0;
    opts->min_patch_rows = 4;
    opts->max_bank_shift_columns = 6;
    opts->sample_columns = 0;
    opts->fill_internal_gaps = 0;
    opts->fill_vertical_gaps = 0;
    opts->allow_component_drift = 0;
    opts->shared_lattice_fill = 0;
    opts->fill_aligned_gaps = 0;
    opts->short_strip_max_span_u = 8;
    opts->short_strip_max_span_v = 8;
    opts->tangent_dot_min = 0.2;
    opts->normal_cone_dot_min = 0.25;
    opts->chord_min_fraction = 0.15;
    opts->chord_max_stretch = 1.8;
    opts->chord_slack = 4.0;
    opts->vertical_stretch = 2.0;
    opts->vertical_slack = 4.0;
    opts->bank_fit_rms_max = 2.0;
    opts->bank_fit_max_max = 6.0;
    opts->hermite_prior_weight = 0.2;
    opts->boundary_weight = 64.0;
    opts->nuclear_lambda = 0.05;
    opts->rank_tail_lambda = 0.2;
    opts->admm_rho = 1.0;
    opts->nuclear_iterations = 8;
    opts->rank_tail_iterations = 24;
    opts->pcg_iterations = 160;
    opts->pcg_tolerance = 1.0e-7;
    opts->max_patches = 4096;
}

static int adc_options_valid(const AtlasDevelopableCompleteOptions *opts)
{
    return opts != NULL &&
           isfinite(opts->min_gap_u) && opts->min_gap_u > 0.0 &&
           isfinite(opts->max_gap_u) &&
           opts->max_gap_u >= opts->min_gap_u &&
           isfinite(opts->min_gap_v) && opts->min_gap_v > 0.0 &&
           isfinite(opts->max_gap_v) &&
           opts->max_gap_v >= opts->min_gap_v &&
           opts->min_patch_rows >= 2 &&
           opts->max_bank_shift_columns >= 0 &&
           opts->sample_columns >= 0 &&
           (opts->fill_internal_gaps == 0 ||
            opts->fill_internal_gaps == 1) &&
           (opts->fill_vertical_gaps == 0 ||
            opts->fill_vertical_gaps == 1) &&
           (opts->allow_component_drift == 0 ||
             opts->allow_component_drift == 1) &&
           (opts->shared_lattice_fill == 0 ||
             opts->shared_lattice_fill == 1) &&
           (opts->fill_aligned_gaps == 0 ||
            opts->fill_aligned_gaps == 1) &&
           (!opts->fill_aligned_gaps || opts->shared_lattice_fill) &&
           opts->short_strip_max_span_u >= 2 &&
           opts->short_strip_max_span_u <= 1024 &&
           opts->short_strip_max_span_v >= 2 &&
           opts->short_strip_max_span_v <= 1024 &&
           isfinite(opts->tangent_dot_min) &&
           opts->tangent_dot_min >= -1.0 &&
           opts->tangent_dot_min <= 1.0 &&
           isfinite(opts->normal_cone_dot_min) &&
           opts->normal_cone_dot_min >= 0.0 &&
           opts->normal_cone_dot_min <= 1.0 &&
           isfinite(opts->chord_min_fraction) &&
           opts->chord_min_fraction >= 0.0 &&
           isfinite(opts->chord_max_stretch) &&
           opts->chord_max_stretch >= 1.0 &&
           isfinite(opts->chord_slack) && opts->chord_slack >= 0.0 &&
           isfinite(opts->vertical_stretch) &&
           opts->vertical_stretch >= 1.0 &&
           isfinite(opts->vertical_slack) &&
           opts->vertical_slack >= 0.0 &&
           isfinite(opts->bank_fit_rms_max) &&
           opts->bank_fit_rms_max > 0.0 &&
           isfinite(opts->bank_fit_max_max) &&
           opts->bank_fit_max_max >= opts->bank_fit_rms_max &&
           isfinite(opts->hermite_prior_weight) &&
           opts->hermite_prior_weight > 0.0 &&
           isfinite(opts->boundary_weight) &&
           opts->boundary_weight >= opts->hermite_prior_weight &&
           isfinite(opts->nuclear_lambda) &&
           opts->nuclear_lambda >= 0.0 &&
           isfinite(opts->rank_tail_lambda) &&
           opts->rank_tail_lambda >= 0.0 &&
           isfinite(opts->admm_rho) && opts->admm_rho > 0.0 &&
           opts->nuclear_iterations >= 0 &&
           opts->rank_tail_iterations >= 0 &&
           opts->pcg_iterations >= 1 &&
           isfinite(opts->pcg_tolerance) &&
           opts->pcg_tolerance > 0.0 &&
           opts->max_patches >= 1;
}

static size_t adc_hessian_count(size_t rows, size_t columns)
{
    if (rows < 3 || columns < 3) return 0;
    return (rows - 2) * (columns - 2);
}

static void adc_hessian_apply(
    const double *x,
    size_t rows,
    size_t columns,
    double inv_du2,
    double inv_dv2,
    double inv_4du_dv,
    double *out)
{
    size_t at = 0;
    for (size_t row = 1; row + 1 < rows; row++) {
        for (size_t column = 1; column + 1 < columns; column++) {
            size_t center = row * columns + column;
            double xx = (x[center - 1] - 2.0 * x[center] +
                         x[center + 1]) * inv_du2;
            double yy = (x[center - columns] - 2.0 * x[center] +
                         x[center + columns]) * inv_dv2;
            double xy = (x[center + columns + 1] -
                         x[center + columns - 1] -
                         x[center - columns + 1] +
                         x[center - columns - 1]) * inv_4du_dv;
            out[at * 3] = xx;
            out[at * 3 + 1] = xy;
            out[at * 3 + 2] = yy;
            at++;
        }
    }
}

static void adc_hessian_adjoint(
    const double *value,
    size_t rows,
    size_t columns,
    double inv_du2,
    double inv_dv2,
    double inv_4du_dv,
    double *out)
{
    size_t nodes = rows * columns;
    memset(out, 0, nodes * sizeof(*out));
    size_t at = 0;
    for (size_t row = 1; row + 1 < rows; row++) {
        for (size_t column = 1; column + 1 < columns; column++) {
            size_t center = row * columns + column;
            double xx = value[at * 3];
            double xy = value[at * 3 + 1];
            double yy = value[at * 3 + 2];
            out[center - 1] += inv_du2 * xx;
            out[center] -= 2.0 * inv_du2 * xx;
            out[center + 1] += inv_du2 * xx;
            out[center - columns] += inv_dv2 * yy;
            out[center] -= 2.0 * inv_dv2 * yy;
            out[center + columns] += inv_dv2 * yy;
            out[center + columns + 1] += inv_4du_dv * xy;
            out[center + columns - 1] -= inv_4du_dv * xy;
            out[center - columns + 1] -= inv_4du_dv * xy;
            out[center - columns - 1] += inv_4du_dv * xy;
            at++;
        }
    }
}

static void adc_build_diagonal(
    const double *weight,
    size_t rows,
    size_t columns,
    double rho,
    double inv_du2,
    double inv_dv2,
    double inv_4du_dv,
    double *diagonal)
{
    size_t nodes = rows * columns;
    memcpy(diagonal, weight, nodes * sizeof(*diagonal));
    double xx = rho * inv_du2 * inv_du2;
    double yy = rho * inv_dv2 * inv_dv2;
    double xy = rho * inv_4du_dv * inv_4du_dv;
    for (size_t row = 1; row + 1 < rows; row++) {
        for (size_t column = 1; column + 1 < columns; column++) {
            size_t center = row * columns + column;
            diagonal[center - 1] += xx;
            diagonal[center] += 4.0 * xx;
            diagonal[center + 1] += xx;
            diagonal[center - columns] += yy;
            diagonal[center] += 4.0 * yy;
            diagonal[center + columns] += yy;
            diagonal[center + columns + 1] += xy;
            diagonal[center + columns - 1] += xy;
            diagonal[center - columns + 1] += xy;
            diagonal[center - columns - 1] += xy;
        }
    }
    for (size_t i = 0; i < nodes; i++)
        if (!(diagonal[i] > DBL_MIN) || !isfinite(diagonal[i]))
            diagonal[i] = 1.0;
}

static void adc_system_apply(
    const double *x,
    const double *weight,
    size_t rows,
    size_t columns,
    double rho,
    double inv_du2,
    double inv_dv2,
    double inv_4du_dv,
    double *hessian,
    double *adjoint,
    double *out)
{
    size_t nodes = rows * columns;
    adc_hessian_apply(
        x, rows, columns, inv_du2, inv_dv2, inv_4du_dv, hessian);
    adc_hessian_adjoint(
        hessian, rows, columns, inv_du2, inv_dv2, inv_4du_dv,
        adjoint);
    for (size_t i = 0; i < nodes; i++)
        out[i] = weight[i] * x[i] + rho * adjoint[i];
}

static double adc_dot(const double *a, const double *b, size_t count)
{
    double value = 0.0;
    for (size_t i = 0; i < count; i++) value += a[i] * b[i];
    return value;
}

static int adc_pcg(
    const double *weight,
    const double *diagonal,
    const double *rhs,
    size_t rows,
    size_t columns,
    double rho,
    double inv_du2,
    double inv_dv2,
    double inv_4du_dv,
    int max_iterations,
    double tolerance,
    double *hessian,
    double *adjoint,
    double *residual,
    double *preconditioned,
    double *direction,
    double *product,
    double *x)
{
    size_t nodes = rows * columns;
    adc_system_apply(
        x, weight, rows, columns, rho,
        inv_du2, inv_dv2, inv_4du_dv,
        hessian, adjoint, product);
    for (size_t i = 0; i < nodes; i++) {
        residual[i] = rhs[i] - product[i];
        preconditioned[i] = residual[i] / diagonal[i];
        direction[i] = preconditioned[i];
    }
    double rz = adc_dot(residual, preconditioned, nodes);
    double rhs_norm = sqrt(adc_dot(rhs, rhs, nodes));
    double target = tolerance * (rhs_norm > 1.0 ? rhs_norm : 1.0);
    double residual_norm = sqrt(adc_dot(residual, residual, nodes));
    if (residual_norm <= target) return 0;
    if (!(rz > 0.0) || !isfinite(rz)) return -1;

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        adc_system_apply(
            direction, weight, rows, columns, rho,
            inv_du2, inv_dv2, inv_4du_dv,
            hessian, adjoint, product);
        double denominator = adc_dot(direction, product, nodes);
        if (!(denominator > DBL_MIN) || !isfinite(denominator)) return -1;
        double alpha = rz / denominator;
        for (size_t i = 0; i < nodes; i++) {
            x[i] += alpha * direction[i];
            residual[i] -= alpha * product[i];
        }
        residual_norm = sqrt(adc_dot(residual, residual, nodes));
        if (residual_norm <= target) return iteration + 1;
        for (size_t i = 0; i < nodes; i++)
            preconditioned[i] = residual[i] / diagonal[i];
        double next_rz = adc_dot(residual, preconditioned, nodes);
        if (!(next_rz > 0.0) || !isfinite(next_rz)) return -1;
        double beta = next_rz / rz;
        for (size_t i = 0; i < nodes; i++)
            direction[i] = preconditioned[i] + beta * direction[i];
        rz = next_rz;
    }
    return max_iterations;
}

static void adc_eigenvalues2(
    double a,
    double b,
    double c,
    double *lambda0,
    double *lambda1,
    double *cosine,
    double *sine)
{
    double angle = 0.5 * atan2(2.0 * b, a - c);
    double cs = cos(angle);
    double sn = sin(angle);
    double trace = 0.5 * (a + c);
    double radius = hypot(0.5 * (a - c), b);
    *lambda0 = trace + radius;
    *lambda1 = trace - radius;
    *cosine = cs;
    *sine = sn;
}

static double adc_soft_threshold(double value, double threshold)
{
    double magnitude = fabs(value);
    if (magnitude <= threshold) return 0.0;
    return copysign(magnitude - threshold, value);
}

static void adc_prox_hessian(
    const double input[3],
    double threshold,
    int rank_tail,
    double output[3])
{
    double lambda0 = 0.0, lambda1 = 0.0;
    double cosine = 1.0, sine = 0.0;
    adc_eigenvalues2(
        input[0], input[1], input[2],
        &lambda0, &lambda1, &cosine, &sine);
    if (rank_tail) {
        if (fabs(lambda0) >= fabs(lambda1))
            lambda1 = adc_soft_threshold(lambda1, threshold);
        else
            lambda0 = adc_soft_threshold(lambda0, threshold);
    } else {
        lambda0 = adc_soft_threshold(lambda0, threshold);
        lambda1 = adc_soft_threshold(lambda1, threshold);
    }
    double cc = cosine * cosine;
    double ss = sine * sine;
    double cs = cosine * sine;
    output[0] = lambda0 * cc + lambda1 * ss;
    output[1] = (lambda0 - lambda1) * cs;
    output[2] = lambda0 * ss + lambda1 * cc;
}

static double adc_rank_tail_measure(const double *hessian, size_t count)
{
    if (count == 0) return 0.0;
    double total = 0.0;
    for (size_t i = 0; i < count; i++) {
        double lambda0 = 0.0, lambda1 = 0.0;
        double cosine = 1.0, sine = 0.0;
        adc_eigenvalues2(
            hessian[i * 3], hessian[i * 3 + 1],
            hessian[i * 3 + 2],
            &lambda0, &lambda1, &cosine, &sine);
        total += fmin(fabs(lambda0), fabs(lambda1));
    }
    return total / (double)count;
}

static int adc_admm_solve(
    Arena_T arena,
    const double *prior,
    const double *weight,
    size_t rows,
    size_t columns,
    double du,
    double dv,
    const AtlasDevelopableCompleteOptions *opts,
    double *value,
    AdcAdmmStats *stats)
{
    memset(stats, 0, sizeof *stats);
    size_t nodes = rows * columns;
    size_t hessian_count = adc_hessian_count(rows, columns);
    if (nodes == 0 || prior == NULL || weight == NULL || value == NULL)
        return -1;
    memcpy(value, prior, nodes * sizeof(*value));
    if (hessian_count == 0) return 0;
    if (!(du > 0.0) || !(dv > 0.0)) return -1;

    Arena_Mark mark = Arena_save(arena);
    size_t coefficients = hessian_count * 3;
    double *hx = (double *)ARENA_ALLOC(
        arena, coefficients * sizeof(*hx));
    double *z = (double *)ARENA_ALLOC(
        arena, coefficients * sizeof(*z));
    double *z_previous = (double *)ARENA_ALLOC(
        arena, coefficients * sizeof(*z_previous));
    double *dual = (double *)ARENA_CALLOC(
        arena, coefficients, sizeof(*dual));
    double *prox_input = (double *)ARENA_ALLOC(
        arena, coefficients * sizeof(*prox_input));
    double *adjoint = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*adjoint));
    double *rhs = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*rhs));
    double *diagonal = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*diagonal));
    double *residual = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*residual));
    double *preconditioned = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*preconditioned));
    double *direction = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*direction));
    double *product = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*product));

    double inv_du2 = 1.0 / (du * du);
    double inv_dv2 = 1.0 / (dv * dv);
    double inv_4du_dv = 1.0 / (4.0 * du * dv);
    adc_hessian_apply(
        value, rows, columns, inv_du2, inv_dv2, inv_4du_dv, hx);
    memcpy(z, hx, coefficients * sizeof(*z));
    stats->tail_before = adc_rank_tail_measure(hx, hessian_count);
    adc_build_diagonal(
        weight, rows, columns, opts->admm_rho,
        inv_du2, inv_dv2, inv_4du_dv, diagonal);

    int total_iterations =
        opts->nuclear_iterations + opts->rank_tail_iterations;
    int current_phase = -1;
    for (int iteration = 0; iteration < total_iterations; iteration++) {
        int rank_tail = iteration >= opts->nuclear_iterations;
        int phase = rank_tail ? 1 : 0;
        if (phase != current_phase) {
            memset(dual, 0, coefficients * sizeof(*dual));
            adc_hessian_apply(
                value, rows, columns,
                inv_du2, inv_dv2, inv_4du_dv, z);
            current_phase = phase;
        }
        for (size_t i = 0; i < coefficients; i++)
            prox_input[i] = z[i] - dual[i];
        adc_hessian_adjoint(
            prox_input, rows, columns,
            inv_du2, inv_dv2, inv_4du_dv, adjoint);
        for (size_t i = 0; i < nodes; i++)
            rhs[i] = weight[i] * prior[i] +
                     opts->admm_rho * adjoint[i];

        int pcg = adc_pcg(
            weight, diagonal, rhs, rows, columns, opts->admm_rho,
            inv_du2, inv_dv2, inv_4du_dv,
            opts->pcg_iterations, opts->pcg_tolerance,
            hx, adjoint, residual, preconditioned,
            direction, product, value);
        if (pcg < 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        stats->pcg_iterations += pcg;
        adc_hessian_apply(
            value, rows, columns,
            inv_du2, inv_dv2, inv_4du_dv, hx);
        memcpy(z_previous, z, coefficients * sizeof(*z_previous));
        double lambda = rank_tail
            ? opts->rank_tail_lambda : opts->nuclear_lambda;
        double threshold = lambda / opts->admm_rho;
        double primal2 = 0.0, dual2 = 0.0;
        for (size_t i = 0; i < hessian_count; i++) {
            double input[3] = {
                hx[i * 3] + dual[i * 3],
                hx[i * 3 + 1] + dual[i * 3 + 1],
                hx[i * 3 + 2] + dual[i * 3 + 2]
            };
            adc_prox_hessian(input, threshold, rank_tail, &z[i * 3]);
            for (int d = 0; d < 3; d++) {
                size_t index = i * 3 + (size_t)d;
                double primal = hx[index] - z[index];
                double dual_change = z[index] - z_previous[index];
                dual[index] += primal;
                primal2 += primal * primal;
                dual2 += dual_change * dual_change;
            }
        }
        stats->primal = sqrt(primal2 / (double)coefficients);
        stats->dual =
            opts->admm_rho * sqrt(dual2 / (double)coefficients);
        stats->iterations = iteration + 1;
        if (iteration + 1 >= opts->nuclear_iterations &&
            stats->primal < 1.0e-7 && stats->dual < 1.0e-7)
            break;
    }
    adc_hessian_apply(
        value, rows, columns, inv_du2, inv_dv2, inv_4du_dv, hx);
    stats->tail_after = adc_rank_tail_measure(hx, hessian_count);
    Arena_restore(arena, mark);
    return 0;
}

static int adc_copy_seed(
    Arena_T arena,
    const AtlasTrackGrowResult *seed,
    size_t vertex_capacity,
    size_t face_capacity,
    AtlasTrackGrowResult *mesh)
{
    if (seed == NULL || mesh == NULL || vertex_capacity < seed->nv ||
        face_capacity < seed->nf)
        return -1;
    *mesh = *seed;
    mesh->xyz = (double *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) * 3 *
               sizeof(*mesh->xyz));
    mesh->local_uv = (double *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) * 2 *
               sizeof(*mesh->local_uv));
    mesh->global_uv = (double *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) * 2 *
               sizeof(*mesh->global_uv));
    mesh->uv = (double *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) * 2 *
               sizeof(*mesh->uv));
    mesh->source_target = (int32_t *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) *
               sizeof(*mesh->source_target));
    mesh->vertex_component = (int32_t *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) *
               sizeof(*mesh->vertex_component));
    mesh->vertex_winding = (int32_t *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) *
               sizeof(*mesh->vertex_winding));
    mesh->vertex_rank = (uint64_t *)ARENA_ALLOC(
        arena, (vertex_capacity ? vertex_capacity : 1) *
               sizeof(*mesh->vertex_rank));
    mesh->faces = (int32_t *)ARENA_ALLOC(
        arena, (face_capacity ? face_capacity : 1) * 3 *
               sizeof(*mesh->faces));
    mesh->face_component = (int32_t *)ARENA_ALLOC(
        arena, (face_capacity ? face_capacity : 1) *
               sizeof(*mesh->face_component));
    mesh->face_inferred = (uint8_t *)ARENA_ALLOC(
        arena, (face_capacity ? face_capacity : 1) *
               sizeof(*mesh->face_inferred));
    if (seed->nv > 0) {
        memcpy(mesh->xyz, seed->xyz,
               seed->nv * 3 * sizeof(*mesh->xyz));
        memcpy(mesh->local_uv, seed->local_uv,
               seed->nv * 2 * sizeof(*mesh->local_uv));
        memcpy(mesh->global_uv, seed->global_uv,
               seed->nv * 2 * sizeof(*mesh->global_uv));
        memcpy(mesh->uv, seed->uv,
               seed->nv * 2 * sizeof(*mesh->uv));
        memcpy(mesh->source_target, seed->source_target,
               seed->nv * sizeof(*mesh->source_target));
        memcpy(mesh->vertex_component, seed->vertex_component,
               seed->nv * sizeof(*mesh->vertex_component));
        memcpy(mesh->vertex_winding, seed->vertex_winding,
               seed->nv * sizeof(*mesh->vertex_winding));
        memcpy(mesh->vertex_rank, seed->vertex_rank,
               seed->nv * sizeof(*mesh->vertex_rank));
    }
    if (seed->nf > 0) {
        memcpy(mesh->faces, seed->faces,
               seed->nf * 3 * sizeof(*mesh->faces));
        memcpy(mesh->face_component, seed->face_component,
               seed->nf * sizeof(*mesh->face_component));
        memcpy(mesh->face_inferred, seed->face_inferred,
               seed->nf * sizeof(*mesh->face_inferred));
    }
    return 0;
}

static int adc_seed_existing_cells(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    uint8_t **out_cell,
    size_t *out_rows,
    size_t *out_columns)
{
    *out_cell = NULL;
    *out_rows = evidence->nrows;
    *out_columns = 0;
    int32_t maximum_column = -1;
    for (size_t i = 0; i < evidence->ntarget; i++) {
        if (evidence->target[i].column > maximum_column)
            maximum_column = evidence->target[i].column;
    }
    if (*out_rows == 0 || maximum_column < 0) return 0;
    *out_columns = (size_t)maximum_column + 1;
    if (*out_rows > SIZE_MAX / *out_columns) return -1;
    size_t count = *out_rows * *out_columns;
    uint8_t *cell = (uint8_t *)ARENA_CALLOC(
        arena, count, sizeof(*cell));
    for (size_t face = 0; face < seed->nf; face++) {
        int32_t row_min = INT32_MAX, row_max = INT32_MIN;
        int32_t column_min = INT32_MAX, column_max = INT32_MIN;
        int valid = 1;
        for (int corner = 0; corner < 3; corner++) {
            int32_t vertex = seed->faces[face * 3 + (size_t)corner];
            if (vertex < 0 || (size_t)vertex >= seed->nv) {
                valid = 0;
                break;
            }
            int32_t target_index = seed->source_target[(size_t)vertex];
            if (target_index < 0 ||
                (size_t)target_index >= evidence->ntarget) {
                valid = 0;
                break;
            }
            const AtlasRibbonTarget *target = &evidence->target[target_index];
            if (target->row < row_min) row_min = target->row;
            if (target->row > row_max) row_max = target->row;
            if (target->column < column_min) column_min = target->column;
            if (target->column > column_max) column_max = target->column;
        }
        if (valid && row_max == row_min + 1 &&
            column_max == column_min + 1 && row_min >= 0 &&
            column_min >= 0 && (size_t)row_max < *out_rows &&
            (size_t)column_max < *out_columns)
            cell[(size_t)row_min * *out_columns +
                 (size_t)column_min] = 1;
    }
    *out_cell = cell;
    return 0;
}

static int adc_existing_u_edge(
    const uint8_t *cell,
    size_t rows,
    size_t columns,
    int32_t row,
    int32_t column)
{
    if (cell == NULL || row < 0 || column < 0 ||
        (size_t)row >= rows || (size_t)column + 1 >= columns)
        return 0;
    if ((size_t)row + 1 < rows &&
        cell[(size_t)row * columns + (size_t)column])
        return 1;
    return row > 0 &&
        cell[((size_t)row - 1) * columns + (size_t)column];
}

static int adc_existing_v_edge(
    const uint8_t *cell,
    size_t rows,
    size_t columns,
    int32_t row,
    int32_t column)
{
    if (cell == NULL || row < 0 || column < 0 ||
        (size_t)row + 1 >= rows || (size_t)column >= columns)
        return 0;
    if ((size_t)column + 1 < columns &&
        cell[(size_t)row * columns + (size_t)column])
        return 1;
    return column > 0 &&
        cell[(size_t)row * columns + (size_t)column - 1];
}

static int adc_gap_from_pair(
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    const uint8_t *existing_cell,
    size_t cell_rows,
    size_t cell_columns,
    const AdcSelected *left,
    const AdcSelected *right,
    AdcRowGap *out)
{
    if ((!opts->fill_internal_gaps &&
         left->component == right->component) ||
        right->column <= left->column)
        return 0;
    if (right->column == left->column + 1 &&
        adc_existing_u_edge(
            existing_cell, cell_rows, cell_columns,
            left->row, left->column))
        return 0;
    const AtlasRibbonTarget *a = &evidence->target[left->target];
    const AtlasRibbonTarget *b = &evidence->target[right->target];
    double gap_u = b->u - a->u;
    if (!(gap_u >= opts->min_gap_u && gap_u <= opts->max_gap_u))
        return 0;
    double chord = adc_distance2(a->p, b->p);
    if (chord + opts->chord_slack <
            opts->chord_min_fraction * gap_u ||
        chord > opts->chord_max_stretch * gap_u + opts->chord_slack)
        return 0;
    double tangent_dot = adc_dot2(a->tangent, b->tangent);
    if (tangent_dot < opts->tangent_dot_min) return 0;
    int64_t winding_delta =
        (int64_t)seed->vertex_winding[right->vertex] -
        (int64_t)seed->vertex_winding[left->vertex];
    if (winding_delta < -2 || winding_delta > 2) return 0;
    double chord_direction[2] = {
        b->p[0] - a->p[0], b->p[1] - a->p[1]
    };
    if (chord > 1.0e-12) {
        chord_direction[0] /= chord;
        chord_direction[1] /= chord;
        if (adc_dot2(a->tangent, chord_direction) < -0.35 ||
            adc_dot2(b->tangent, chord_direction) < -0.35)
            return 0;
    }
    memset(out, 0, sizeof *out);
    out->direction = 0;
    out->row = left->row;
    out->column_left = left->column;
    out->column_right = right->column;
    out->target_left = left->target;
    out->target_right = right->target;
    out->vertex_left = left->vertex;
    out->vertex_right = right->vertex;
    out->component_left = left->component;
    out->component_right = right->component;
    out->gap_u = gap_u;
    out->chord = chord;
    out->score = fabs(chord / gap_u - 1.0) + 0.5 * (1.0 - tangent_dot);
    return 1;
}

static int adc_gap_v_from_pair(
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    const uint8_t *existing_cell,
    size_t cell_rows,
    size_t cell_columns,
    const AdcSelected *top,
    const AdcSelected *bottom,
    AdcRowGap *out)
{
    if ((!opts->fill_internal_gaps &&
         top->component == bottom->component) ||
        bottom->row <= top->row)
        return 0;
    if (bottom->row == top->row + 1 &&
        adc_existing_v_edge(
            existing_cell, cell_rows, cell_columns,
            top->row, top->column))
        return 0;
    const AtlasRibbonTarget *a = &evidence->target[top->target];
    const AtlasRibbonTarget *b = &evidence->target[bottom->target];
    double gap_v = b->v - a->v;
    if (!(gap_v >= opts->min_gap_v && gap_v <= opts->max_gap_v))
        return 0;
    double support_slack = evidence->observation_du > 0.0
        ? 2.0 * evidence->observation_du : 4.0;
    if (fabs(b->u - a->u) > support_slack) return 0;
    const double *xyz_top = &seed->xyz[(size_t)top->vertex * 3];
    const double *xyz_bottom = &seed->xyz[(size_t)bottom->vertex * 3];
    double chord = adc_distance3(xyz_top, xyz_bottom);
    if (chord + opts->chord_slack <
            opts->chord_min_fraction * gap_v ||
        chord > opts->chord_max_stretch * gap_v + opts->chord_slack)
        return 0;
    double tangent_dot = adc_dot2(a->tangent, b->tangent);
    if (tangent_dot < opts->tangent_dot_min) return 0;
    int64_t winding_delta =
        (int64_t)seed->vertex_winding[bottom->vertex] -
        (int64_t)seed->vertex_winding[top->vertex];
    if (winding_delta < -1 || winding_delta > 1) return 0;

    memset(out, 0, sizeof *out);
    out->direction = 1;
    out->row = top->column;
    out->column_left = top->row;
    out->column_right = bottom->row;
    out->target_left = top->target;
    out->target_right = bottom->target;
    out->vertex_left = top->vertex;
    out->vertex_right = bottom->vertex;
    out->component_left = top->component;
    out->component_right = bottom->component;
    out->gap_u = gap_v;
    out->chord = chord;
    out->score = fabs(chord / gap_v - 1.0) +
                 0.5 * (1.0 - tangent_dot);
    return 1;
}

static int adc_build_gaps(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    AdcRowGap **out_gap,
    size_t *out_count)
{
    *out_gap = NULL;
    *out_count = 0;
    if (seed->nv == 0 || evidence->ntarget == 0) return 0;
    AdcSelected *selected = (AdcSelected *)ARENA_ALLOC(
        arena, seed->nv * sizeof(*selected));
    size_t selected_count = 0;
    for (size_t vertex = 0; vertex < seed->nv; vertex++) {
        int32_t target = seed->source_target[vertex];
        if (target < 0 || (size_t)target >= evidence->ntarget) continue;
        const AtlasRibbonTarget *sample = &evidence->target[target];
        if (!sample->accepted) continue;
        AdcSelected *entry = &selected[selected_count++];
        entry->row = sample->row;
        entry->column = sample->column;
        entry->target = target;
        entry->vertex = (int32_t)vertex;
        entry->component = seed->vertex_component[vertex];
    }
    if (selected_count < 2) return 0;
    qsort(selected, selected_count, sizeof(*selected), adc_compare_selected);

    uint8_t *existing_cell = NULL;
    size_t cell_rows = 0, cell_columns = 0;
    if (adc_seed_existing_cells(
            arena, evidence, seed, &existing_cell,
            &cell_rows, &cell_columns) != 0)
        return -1;

    size_t unique_count = 0;
    for (size_t i = 0; i < selected_count; i++) {
        if (unique_count > 0 &&
            selected[i].target == selected[unique_count - 1].target &&
            selected[i].component == selected[unique_count - 1].component)
            continue;
        selected[unique_count++] = selected[i];
    }
    if (unique_count > SIZE_MAX / (2 * sizeof(AdcRowGap))) return -1;
    AdcRowGap *gap = (AdcRowGap *)ARENA_ALLOC(
        arena, unique_count * 2 * sizeof(*gap));
    size_t gap_count = 0;
    size_t first = 0;
    while (first < unique_count) {
        size_t last = first + 1;
        while (last < unique_count &&
               selected[last].row == selected[first].row)
            last++;
        const AdcSelected *previous = &selected[first];
        for (size_t i = first + 1; i < last; i++) {
            const AdcSelected *current = &selected[i];
            if (current->column == previous->column) {
                if (current->component == previous->component)
                    previous = current;
                continue;
            }
            if (adc_gap_from_pair(
                    evidence, seed, opts, existing_cell,
                    cell_rows, cell_columns, previous, current,
                    &gap[gap_count]))
                gap_count++;
            previous = current;
        }
        first = last;
    }
    if (opts->fill_vertical_gaps) {
        qsort(selected, unique_count, sizeof(*selected),
              adc_compare_selected_column);
        first = 0;
        while (first < unique_count) {
            size_t last = first + 1;
            while (last < unique_count &&
                   selected[last].column == selected[first].column)
                last++;
            const AdcSelected *previous = &selected[first];
            for (size_t i = first + 1; i < last; i++) {
                const AdcSelected *current = &selected[i];
                if (current->row == previous->row) {
                    if (current->component == previous->component)
                        previous = current;
                    continue;
                }
                if (adc_gap_v_from_pair(
                        evidence, seed, opts, existing_cell,
                        cell_rows, cell_columns, previous, current,
                        &gap[gap_count]))
                    gap_count++;
                previous = current;
            }
            first = last;
        }
    }
    if (gap_count == 0) return 0;
    qsort(gap, gap_count, sizeof(*gap),
          opts->allow_component_drift
              ? adc_compare_gap_spatial : adc_compare_gap);

    size_t compact = 0;
    for (size_t i = 0; i < gap_count;) {
        size_t last = i + 1;
        if (opts->allow_component_drift) {
            while (last < gap_count &&
                   gap[last].direction == gap[i].direction &&
                   gap[last].row == gap[i].row &&
                   gap[last].column_left == gap[i].column_left &&
                   gap[last].column_right == gap[i].column_right)
                last++;
        } else {
            while (last < gap_count &&
                    gap[last].direction == gap[i].direction &&
                    gap[last].component_left == gap[i].component_left &&
                    gap[last].component_right == gap[i].component_right &&
                    gap[last].row == gap[i].row &&
                    gap[last].column_left == gap[i].column_left &&
                    gap[last].column_right == gap[i].column_right)
                last++;
        }
        gap[compact++] = gap[i];
        i = last;
    }
    *out_gap = gap;
    *out_count = compact;
    return 0;
}

static int adc_rows_connect(
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *previous,
    const AdcRowGap *current)
{
    if (current->direction != previous->direction) return 0;
    if (current->row != previous->row + 1) return 0;
    if (abs(current->column_left - previous->column_left) >
            opts->max_bank_shift_columns ||
        abs(current->column_right - previous->column_right) >
            opts->max_bank_shift_columns)
        return 0;
    double support_spacing = previous->direction == 0
        ? evidence->dv : evidence->observation_du;
    if (!(support_spacing > 0.0)) support_spacing = 1.0;
    double maximum = opts->vertical_stretch * support_spacing +
                     opts->vertical_slack;
    const double *left0 = &seed->xyz[(size_t)previous->vertex_left * 3];
    const double *left1 = &seed->xyz[(size_t)current->vertex_left * 3];
    const double *right0 = &seed->xyz[(size_t)previous->vertex_right * 3];
    const double *right1 = &seed->xyz[(size_t)current->vertex_right * 3];
    return adc_distance3(left0, left1) <= maximum &&
           adc_distance3(right0, right1) <= maximum;
}

static int adc_build_runs(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *gap,
    size_t gap_count,
    AdcRowGap **out_ordered_gap,
    AdcRun **out_run,
    size_t *out_count,
    size_t *out_pairs)
{
    *out_ordered_gap = NULL;
    *out_run = NULL;
    *out_count = 0;
    *out_pairs = 0;
    if (gap_count == 0) return 0;
    AdcRowGap *ordered = (AdcRowGap *)ARENA_ALLOC(
        arena, gap_count * sizeof(*ordered));
    AdcRun *run = (AdcRun *)ARENA_ALLOC(
        arena, gap_count * sizeof(*run));
    size_t ordered_count = 0, run_count = 0;
    size_t pair_count = 0;

    if (!opts->allow_component_drift) {
        for (size_t i = 0; i < gap_count; i++) {
            if (i == 0 || gap[i].direction != gap[i - 1].direction ||
                gap[i].component_left != gap[i - 1].component_left ||
                gap[i].component_right != gap[i - 1].component_right)
                pair_count++;
        }
    }

    /* A component pair does not identify a unique hole: the same two sheets
     * can bound several disjoint intervals in one support column.  Track each
     * spatial interval independently.  In exact mode, continuations must keep
     * the component pair; drift mode may change it but pays a small cost. */
    AdcRowGap *spatial = (AdcRowGap *)ARENA_ALLOC(
        arena, gap_count * sizeof(*spatial));
    memcpy(spatial, gap, gap_count * sizeof(*spatial));
    qsort(spatial, gap_count, sizeof(*spatial), adc_compare_gap_spatial);
    uint8_t *used = (uint8_t *)ARENA_CALLOC(
        arena, gap_count, sizeof(*used));
    size_t support_rows = 0;
    for (size_t i = 0; i < gap_count; i++) {
        if (spatial[i].row < 0) return -1;
        size_t extent = (size_t)spatial[i].row + 1;
        if (extent > support_rows) support_rows = extent;
    }
    if (support_rows == 0 || support_rows > SIZE_MAX / 2) return -1;
    size_t buckets = support_rows * 2;
    size_t *row_begin = (size_t *)ARENA_CALLOC(
        arena, buckets + 1, sizeof(*row_begin));
    for (size_t i = 0; i < gap_count; i++) {
        size_t key = (size_t)spatial[i].direction * support_rows +
                     (size_t)spatial[i].row;
        if (key >= buckets) return -1;
        row_begin[key + 1]++;
    }
    for (size_t key = 0; key < buckets; key++)
        row_begin[key + 1] += row_begin[key];

    for (size_t start = 0; start < gap_count; start++) {
        if (used[start]) continue;
        size_t current = start;
        run[run_count].first = ordered_count;
        run[run_count].component_left = spatial[start].component_left;
        run[run_count].component_right = spatial[start].component_right;
        do {
            used[current] = 1;
            ordered[ordered_count++] = spatial[current];
            int32_t next_row = spatial[current].row + 1;
            size_t best = SIZE_MAX;
            double best_cost = DBL_MAX;
            if (next_row >= 0 && (size_t)next_row < support_rows) {
                size_t key =
                    (size_t)spatial[current].direction * support_rows +
                    (size_t)next_row;
                size_t first = row_begin[key];
                size_t last = row_begin[key + 1];
                for (size_t candidate = first; candidate < last;
                     candidate++) {
                    if (used[candidate]) continue;
                    if (!opts->allow_component_drift &&
                        (spatial[candidate].component_left !=
                             spatial[current].component_left ||
                         spatial[candidate].component_right !=
                             spatial[current].component_right))
                        continue;
                    if (!adc_rows_connect(
                            evidence, seed, opts,
                            &spatial[current], &spatial[candidate]))
                        continue;
                    double cost = (double)(
                        abs(spatial[candidate].column_left -
                            spatial[current].column_left) +
                        abs(spatial[candidate].column_right -
                            spatial[current].column_right));
                    cost += 0.01 * spatial[candidate].score;
                    if (opts->allow_component_drift &&
                        spatial[candidate].component_left ==
                            spatial[current].component_left &&
                        spatial[candidate].component_right ==
                            spatial[current].component_right)
                        cost -= 0.25;
                    if (cost < best_cost) {
                        best = candidate;
                        best_cost = cost;
                    }
                }
            }
            current = best;
        } while (current != SIZE_MAX);
        run[run_count].count = ordered_count - run[run_count].first;
        run_count++;
    }
    if (opts->allow_component_drift) pair_count = run_count;
    if (ordered_count != gap_count) return -1;
    *out_ordered_gap = ordered;
    *out_run = run;
    *out_count = run_count;
    *out_pairs = pair_count;
    return 0;
}

static size_t adc_patch_columns(
    const AtlasRibbonObservationSet *evidence,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *gap,
    size_t count)
{
    if (opts->sample_columns >= 3)
        return (size_t)opts->sample_columns;
    double maximum = 0.0;
    for (size_t i = 0; i < count; i++)
        if (gap[i].gap_u > maximum) maximum = gap[i].gap_u;
    double spacing = gap[0].direction == 0
        ? evidence->observation_du : evidence->dv;
    if (!(spacing > 0.0)) spacing = 2.0;
    double raw = ceil(maximum / spacing) + 1.0;
    if (raw < 3.0) raw = 3.0;
    if (raw > 1025.0) raw = 1025.0;
    return (size_t)raw;
}

static int adc_solve3(double matrix[3][3], double rhs[3], double out[3])
{
    double augmented[3][4];
    for (int row = 0; row < 3; row++) {
        for (int column = 0; column < 3; column++)
            augmented[row][column] = matrix[row][column];
        augmented[row][3] = rhs[row];
    }
    for (int pivot = 0; pivot < 3; pivot++) {
        int best = pivot;
        for (int row = pivot + 1; row < 3; row++)
            if (fabs(augmented[row][pivot]) >
                fabs(augmented[best][pivot]))
                best = row;
        if (fabs(augmented[best][pivot]) < 1.0e-12) return -1;
        if (best != pivot) {
            for (int column = pivot; column < 4; column++) {
                double temporary = augmented[pivot][column];
                augmented[pivot][column] = augmented[best][column];
                augmented[best][column] = temporary;
            }
        }
        double inverse = 1.0 / augmented[pivot][pivot];
        for (int column = pivot; column < 4; column++)
            augmented[pivot][column] *= inverse;
        for (int row = 0; row < 3; row++) {
            if (row == pivot) continue;
            double factor = augmented[row][pivot];
            for (int column = pivot; column < 4; column++)
                augmented[row][column] -=
                    factor * augmented[pivot][column];
        }
    }
    for (int row = 0; row < 3; row++) out[row] = augmented[row][3];
    return 0;
}

static int adc_bank_fit(
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *gap,
    size_t count,
    double coefficient[3][3],
    double frame_normal[3],
    double *out_rms,
    double *out_maximum)
{
    if (count < 2) return -1;
    double gram[3][3] = {{0.0}};
    double rhs[3][3] = {{0.0}};
    for (size_t row = 0; row < count; row++) {
        double v = count > 1
            ? (double)row / (double)(count - 1) : 0.0;
        int32_t vertex[2] = {
            gap[row].vertex_left, gap[row].vertex_right
        };
        for (int side = 0; side < 2; side++) {
            double q[3] = {1.0, (double)side, v};
            const double *point = &seed->xyz[(size_t)vertex[side] * 3];
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++)
                    gram[a][b] += q[a] * q[b];
                for (int d = 0; d < 3; d++)
                    rhs[d][a] += q[a] * point[d];
            }
        }
    }
    for (int d = 0; d < 3; d++) {
        double matrix[3][3];
        memcpy(matrix, gram, sizeof matrix);
        if (adc_solve3(matrix, rhs[d], coefficient[d]) != 0) return -1;
    }
    double direction_u[3] = {
        coefficient[0][1], coefficient[1][1], coefficient[2][1]
    };
    double direction_v[3] = {
        coefficient[0][2], coefficient[1][2], coefficient[2][2]
    };
    adc_cross3(direction_u, direction_v, frame_normal);
    if (adc_normalize3(frame_normal) != 0) return -1;
    double sum2 = 0.0, maximum = 0.0;
    size_t samples = 0;
    for (size_t row = 0; row < count; row++) {
        double v = count > 1
            ? (double)row / (double)(count - 1) : 0.0;
        int32_t vertex[2] = {
            gap[row].vertex_left, gap[row].vertex_right
        };
        for (int side = 0; side < 2; side++) {
            double q[3] = {1.0, (double)side, v};
            const double *point = &seed->xyz[(size_t)vertex[side] * 3];
            double residual_vector[3];
            double residual2 = 0.0;
            for (int d = 0; d < 3; d++) {
                double predicted = coefficient[d][0] * q[0] +
                                   coefficient[d][1] * q[1] +
                                   coefficient[d][2] * q[2];
                residual_vector[d] = point[d] - predicted;
                residual2 += residual_vector[d] * residual_vector[d];
            }
            double height = adc_dot3(residual_vector, frame_normal);
            double tangent2 = residual2 - height * height;
            if (tangent2 < 0.0) tangent2 = 0.0;
            double tangent_residual = sqrt(tangent2);
            sum2 += tangent2;
            if (tangent_residual > maximum) maximum = tangent_residual;
            samples++;
        }
    }
    *out_rms = samples > 0 ? sqrt(sum2 / (double)samples) : 0.0;
    *out_maximum = maximum;
    return 0;
}

static int adc_patch_overlaps(
    const AtlasDevelopablePatch *patch,
    size_t count,
    int32_t row_first,
    int32_t row_last,
    int32_t column_first,
    int32_t column_last)
{
    for (size_t i = 0; i < count; i++) {
        const AtlasDevelopablePatch *other = &patch[i];
        if (!other->accepted) continue;
        int rows_overlap = row_first <= other->row_last &&
                           other->row_first <= row_last;
        int columns_overlap = column_first < other->column_right_max &&
                              other->column_left_min < column_last;
        if (rows_overlap && columns_overlap) return 1;
    }
    return 0;
}

static void adc_patch_grid_bounds(
    const AdcRowGap *gap,
    size_t count,
    int32_t *row_first,
    int32_t *row_last,
    int32_t *column_first,
    int32_t *column_last)
{
    if (gap[0].direction == 0) {
        *row_first = gap[0].row;
        *row_last = gap[count - 1].row;
        *column_first = gap[0].column_left;
        *column_last = gap[0].column_right;
        for (size_t i = 1; i < count; i++) {
            if (gap[i].column_left < *column_first)
                *column_first = gap[i].column_left;
            if (gap[i].column_right > *column_last)
                *column_last = gap[i].column_right;
        }
    } else {
        *column_first = gap[0].row;
        *column_last = gap[count - 1].row;
        *row_first = gap[0].column_left;
        *row_last = gap[0].column_right;
        for (size_t i = 1; i < count; i++) {
            if (gap[i].column_left < *row_first)
                *row_first = gap[i].column_left;
            if (gap[i].column_right > *row_last)
                *row_last = gap[i].column_right;
        }
    }
}

static double adc_phase_fraction(int direction, const double p[2])
{
    double angle = atan2(p[1], p[0]) * (double)direction;
    double wrapped = fmod(angle, ADC_TAU);
    if (wrapped < 0.0) wrapped += ADC_TAU;
    return wrapped / ADC_TAU;
}

static int32_t adc_interpolated_winding(
    const AtlasTrackGrowResult *seed,
    int32_t left,
    int32_t right,
    const AtlasRibbonTarget *target_left,
    const AtlasRibbonTarget *target_right,
    const double p[2],
    double t)
{
    double left_turn = (double)seed->vertex_winding[left] +
        adc_phase_fraction(seed->winding_direction, target_left->p);
    double right_turn = (double)seed->vertex_winding[right] +
        adc_phase_fraction(seed->winding_direction, target_right->p);
    double turn = (1.0 - t) * left_turn + t * right_turn;
    double phase = adc_phase_fraction(seed->winding_direction, p);
    double winding = nearbyint(turn - phase);
    if (winding < (double)INT32_MIN) return INT32_MIN;
    if (winding > (double)INT32_MAX) return INT32_MAX;
    return (int32_t)winding;
}

static void adc_hermite_prior(
    const AtlasRibbonTarget *left,
    const AtlasRibbonTarget *right,
    double t,
    double output[2])
{
    double t2 = t * t;
    double t3 = t2 * t;
    double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    double h10 = t3 - 2.0 * t2 + t;
    double h01 = -2.0 * t3 + 3.0 * t2;
    double h11 = t3 - t2;
    double gap = right->u - left->u;
    for (int d = 0; d < 2; d++)
        output[d] = h00 * left->p[d] +
                    h10 * gap * left->tangent[d] +
                    h01 * right->p[d] +
                    h11 * gap * right->tangent[d];
}

static int adc_tangent_normal_valid(
    const AtlasRibbonObservationSet *evidence,
    const double tangent[2],
    const double frame_normal[3],
    double minimum_dot,
    int direction)
{
    double world_tangent[3];
    for (int d = 0; d < 3; d++)
        world_tangent[d] = tangent[0] * evidence->basis0[d] +
                           tangent[1] * evidence->basis1[d];
    if (adc_normalize3(world_tangent) != 0) return 0;
    double normal[3];
    if (direction == 0)
        adc_cross3(world_tangent, evidence->axis, normal);
    else
        adc_cross3(evidence->axis, world_tangent, normal);
    if (adc_normalize3(normal) != 0) return 0;
    return adc_dot3(normal, frame_normal) >= minimum_dot;
}

static int adc_normal_cone_valid(
    const AtlasRibbonObservationSet *evidence,
    const AdcRowGap *gap,
    size_t rows,
    const double frame_normal[3],
    double minimum_dot)
{
    static const double sample_t[3] = {0.25, 0.5, 0.75};
    for (size_t row = 0; row < rows; row++) {
        const AtlasRibbonTarget *left =
            &evidence->target[gap[row].target_left];
        const AtlasRibbonTarget *right =
            &evidence->target[gap[row].target_right];
        if (!adc_tangent_normal_valid(
                evidence, left->tangent, frame_normal, minimum_dot,
                gap[row].direction) ||
            !adc_tangent_normal_valid(
                evidence, right->tangent, frame_normal, minimum_dot,
                gap[row].direction))
            return 0;
        double span = gap[row].direction == 0
            ? right->u - left->u : right->v - left->v;
        if (!(span > 1.0e-12)) return 0;
        for (size_t sample = 0; sample < 3; sample++) {
            double t = sample_t[sample];
            double t2 = t * t;
            double dh00 = 6.0 * t2 - 6.0 * t;
            double dh10 = 3.0 * t2 - 4.0 * t + 1.0;
            double dh01 = -dh00;
            double dh11 = 3.0 * t2 - 2.0 * t;
            double tangent[2];
            if (gap[row].direction == 0) {
                for (int d = 0; d < 2; d++)
                    tangent[d] = (dh00 * left->p[d] +
                                  dh10 * span * left->tangent[d] +
                                  dh01 * right->p[d] +
                                  dh11 * span * right->tangent[d]) / span;
            } else {
                for (int d = 0; d < 2; d++)
                    tangent[d] = (1.0 - t) * left->tangent[d] +
                                 t * right->tangent[d];
            }
            if (!adc_tangent_normal_valid(
                    evidence, tangent, frame_normal, minimum_dot,
                    gap[row].direction))
                return 0;
        }
    }
    return 1;
}

static int adc_patch_geometry_valid(
    const double *xyz,
    const AdcRowGap *gap,
    size_t rows,
    size_t columns,
    double dv,
    const double frame_normal[3],
    const AtlasDevelopableCompleteOptions *opts)
{
    for (size_t i = 0; i < rows * columns * 3; i++)
        if (!isfinite(xyz[i])) return 0;
    for (size_t row = 0; row < rows; row++) {
        double expected = gap[row].gap_u / (double)(columns - 1);
        double maximum = opts->chord_max_stretch * expected +
                         opts->chord_slack;
        for (size_t column = 0; column + 1 < columns; column++) {
            const double *a = &xyz[(row * columns + column) * 3];
            const double *b = &xyz[(row * columns + column + 1) * 3];
            double edge = adc_distance3(a, b);
            if (!(edge > 1.0e-8) || edge > maximum) return 0;
        }
    }
    double vertical_maximum = opts->vertical_stretch * dv +
                              opts->vertical_slack;
    for (size_t row = 0; row + 1 < rows; row++) {
        for (size_t column = 0; column < columns; column++) {
            const double *a = &xyz[(row * columns + column) * 3];
            const double *b = &xyz[((row + 1) * columns + column) * 3];
            double edge = adc_distance3(a, b);
            if (!(edge > 1.0e-8) || edge > vertical_maximum) return 0;
        }
    }
    for (size_t row = 0; row + 1 < rows; row++) {
        for (size_t column = 0; column + 1 < columns; column++) {
            const double *a = &xyz[(row * columns + column) * 3];
            const double *b = &xyz[(row * columns + column + 1) * 3];
            const double *c = &xyz[((row + 1) * columns + column) * 3];
            const double *d = &xyz[((row + 1) * columns + column + 1) * 3];
            double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            double dc[3] = {c[0] - d[0], c[1] - d[1], c[2] - d[2]};
            double db[3] = {b[0] - d[0], b[1] - d[1], b[2] - d[2]};
            double n0[3] = {
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0]
            };
            double n1[3] = {
                dc[1] * db[2] - dc[2] * db[1],
                dc[2] * db[0] - dc[0] * db[2],
                dc[0] * db[1] - dc[1] * db[0]
            };
            double length0 = sqrt(
                n0[0] * n0[0] + n0[1] * n0[1] + n0[2] * n0[2]);
            double length1 = sqrt(
                n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
            if (!(length0 > 1.0e-10) || !(length1 > 1.0e-10)) return 0;
            double frame_dot0 = adc_dot3(n0, frame_normal) / length0;
            double frame_dot1 = adc_dot3(n1, frame_normal) / length1;
            if (frame_dot0 < opts->normal_cone_dot_min ||
                frame_dot1 < opts->normal_cone_dot_min)
                return 0;
            double dot = (n0[0] * n1[0] + n0[1] * n1[1] +
                          n0[2] * n1[2]) / (length0 * length1);
            if (dot < -0.25) return 0;
        }
    }
    return 1;
}

static int adc_lattice_add_patch(
    AdcLattice *lattice,
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *gap,
    size_t rows,
    size_t columns,
    const double *xyz)
{
    if (lattice == NULL) return 0;
    if (seed == NULL || gap == NULL || xyz == NULL || rows < 2 || columns < 3)
        return -1;
    for (size_t support = 0; support < rows; support++) {
        int32_t first = gap[support].column_left;
        int32_t last = gap[support].column_right;
        if (first < 0 || last <= first) return -1;
        for (int32_t coordinate = first; coordinate <= last; coordinate++) {
            double t = (double)(coordinate - first) /
                       (double)(last - first);
            double sample = t * (double)(columns - 1);
            size_t a = (size_t)floor(sample);
            if (a + 1 >= columns) a = columns - 2;
            double fraction = sample - (double)a;
            if (t >= 1.0) fraction = 1.0;
            const double *p0 = &xyz[(support * columns + a) * 3];
            const double *p1 = &xyz[(support * columns + a + 1) * 3];
            double point[3];
            for (int d = 0; d < 3; d++)
                point[d] = (1.0 - fraction) * p0[d] + fraction * p1[d];
            int32_t source = t < 0.5
                ? gap[support].vertex_left : gap[support].vertex_right;
            if (source < 0 || (size_t)source >= seed->nv) return -1;
            int32_t row = gap[support].direction == 0
                ? gap[support].row : coordinate;
            int32_t column = gap[support].direction == 0
                ? coordinate : gap[support].row;
            if (adc_lattice_add_node(
                    lattice, row, column, point,
                    gap[support].component_left,
                    seed->vertex_winding[(size_t)source],
                    seed->vertex_rank[(size_t)source]) != 0)
                return -1;
        }
    }
    for (size_t support = 0; support + 1 < rows; support++) {
        int32_t first = gap[support].column_left >
                                gap[support + 1].column_left
            ? gap[support].column_left : gap[support + 1].column_left;
        int32_t last = gap[support].column_right <
                               gap[support + 1].column_right
            ? gap[support].column_right : gap[support + 1].column_right;
        for (int32_t coordinate = first; coordinate < last; coordinate++) {
            int32_t row = gap[support].direction == 0
                ? gap[support].row : coordinate;
            int32_t column = gap[support].direction == 0
                ? coordinate : gap[support].row;
            adc_lattice_mark_cell(lattice, row, column);
        }
    }
    return 0;
}

static int adc_complete_run(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *gap,
    size_t rows,
    size_t vertex_capacity,
    size_t face_capacity,
    AdcLattice *lattice,
    AtlasDevelopablePatch *patch,
    AtlasTrackGrowResult *mesh)
{
    patch->direction = gap[0].direction;
    patch->rows = rows;
    patch->component_left = gap[0].component_left;
    patch->component_right = gap[0].component_right;
    int32_t column_first, column_last;
    adc_patch_grid_bounds(
        gap, rows, &patch->row_first, &patch->row_last,
        &column_first, &column_last);
    if (gap[0].direction == 0) {
        patch->column_left_min = patch->column_left_max =
            gap[0].column_left;
        patch->column_right_min = patch->column_right_max =
            gap[0].column_right;
    } else {
        patch->column_left_min = patch->column_left_max = column_first;
        patch->column_right_min = patch->column_right_max = column_last;
    }
    double gap_sum = 0.0, chord_ratio_sum = 0.0;
    for (size_t row = 0; row < rows; row++) {
        if (row > 0 &&
            (gap[row].component_left != gap[row - 1].component_left ||
             gap[row].component_right != gap[row - 1].component_right))
            patch->component_pair_changes++;
        if (gap[row].direction == 0) {
            if (gap[row].column_left < patch->column_left_min)
                patch->column_left_min = gap[row].column_left;
            if (gap[row].column_left > patch->column_left_max)
                patch->column_left_max = gap[row].column_left;
            if (gap[row].column_right < patch->column_right_min)
                patch->column_right_min = gap[row].column_right;
            if (gap[row].column_right > patch->column_right_max)
                patch->column_right_max = gap[row].column_right;
        }
        gap_sum += gap[row].gap_u;
        if (gap[row].gap_u > patch->gap_u_max)
            patch->gap_u_max = gap[row].gap_u;
        chord_ratio_sum += gap[row].chord / gap[row].gap_u;
    }
    patch->gap_u_mean = gap_sum / (double)rows;
    patch->chord_ratio_mean = chord_ratio_sum / (double)rows;
    if (rows < (size_t)opts->min_patch_rows) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_SHORT;
        return 0;
    }
    double coefficient[3][3];
    double frame_normal[3];
    if (adc_bank_fit(
            seed, gap, rows, coefficient, frame_normal,
            &patch->bank_fit_rms, &patch->bank_fit_max) != 0 ||
        patch->bank_fit_rms > opts->bank_fit_rms_max ||
        patch->bank_fit_max > opts->bank_fit_max_max) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_BANK_METRIC;
        return 0;
    }
    if (!adc_normal_cone_valid(
            evidence, gap, rows, frame_normal,
            opts->normal_cone_dot_min)) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
        return 0;
    }

    size_t columns = adc_patch_columns(evidence, opts, gap, rows);
    patch->columns = columns;
    size_t added_vertices = rows * (columns - 2);
    size_t added_faces = (rows - 1) * (columns - 1) * 2;
    if (added_vertices > vertex_capacity ||
        added_faces > face_capacity ||
        mesh->nv > vertex_capacity - added_vertices ||
        mesh->nf > face_capacity - added_faces ||
        mesh->nv + added_vertices > (size_t)INT32_MAX) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
        return 0;
    }

    if (lattice != NULL) {
        size_t additions = 0;
        for (size_t row = 0; row < rows; row++) {
            if (gap[row].column_right < gap[row].column_left) {
                patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
                return 0;
            }
            size_t count =
                (size_t)(gap[row].column_right - gap[row].column_left) + 1;
            if (count > SIZE_MAX - additions) {
                patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
                return 0;
            }
            additions += count;
        }
        if (adc_lattice_reserve_additions(lattice, additions) != 0) {
            patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
            return 0;
        }
    }

    Arena_Mark mark = Arena_save(arena);
    size_t nodes = rows * columns;
    double *prior = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*prior));
    double *value = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*value));
    double *weight = (double *)ARENA_ALLOC(
        arena, nodes * sizeof(*weight));
    double *xyz = (double *)ARENA_ALLOC(
        arena, nodes * 3 * sizeof(*xyz));
    for (size_t row = 0; row < rows; row++) {
        const AtlasRibbonTarget *left =
            &evidence->target[gap[row].target_left];
        const AtlasRibbonTarget *right =
            &evidence->target[gap[row].target_right];
        double row_coordinate = rows > 1
            ? (double)row / (double)(rows - 1) : 0.0;
        for (size_t column = 0; column < columns; column++) {
            double t = (double)column / (double)(columns - 1);
            double p[2];
            if (gap[row].direction == 0) {
                adc_hermite_prior(left, right, t, p);
            } else {
                for (int d = 0; d < 2; d++)
                    p[d] = (1.0 - t) * left->p[d] +
                           t * right->p[d];
            }
            size_t index = row * columns + column;
            double axial = (1.0 - t) * left->v + t * right->v;
            double observed[3];
            for (int d = 0; d < 3; d++)
                observed[d] = evidence->axis_point[d] +
                    axial * evidence->axis[d] +
                    p[0] * evidence->basis0[d] +
                    p[1] * evidence->basis1[d];
            if (column == 0)
                memcpy(observed,
                       &seed->xyz[(size_t)gap[row].vertex_left * 3],
                       sizeof observed);
            else if (column + 1 == columns)
                memcpy(observed,
                       &seed->xyz[(size_t)gap[row].vertex_right * 3],
                       sizeof observed);
            double residual[3];
            for (int d = 0; d < 3; d++) {
                double base = coefficient[d][0] +
                    coefficient[d][1] * t +
                    coefficient[d][2] * row_coordinate;
                residual[d] = observed[d] - base;
            }
            prior[index] = adc_dot3(residual, frame_normal);
            weight[index] = column == 0 || column + 1 == columns
                ? opts->boundary_weight : opts->hermite_prior_weight;
        }
    }
    double direction_u[3] = {
        coefficient[0][1], coefficient[1][1], coefficient[2][1]
    };
    double direction_v[3] = {
        coefficient[0][2], coefficient[1][2], coefficient[2][2]
    };
    double du = sqrt(adc_dot3(direction_u, direction_u)) /
                (double)(columns - 1);
    double dv = sqrt(adc_dot3(direction_v, direction_v)) /
                (double)(rows - 1);
    AdcAdmmStats stats;
    if (adc_admm_solve(
            arena, prior, weight, rows, columns, du, dv,
            opts, value, &stats) != 0) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
        Arena_restore(arena, mark);
        return 0;
    }
    for (size_t row = 0; row < rows; row++) {
        double row_coordinate = rows > 1
            ? (double)row / (double)(rows - 1) : 0.0;
        for (size_t column = 0; column < columns; column++) {
            size_t index = row * columns + column;
            double t = (double)column / (double)(columns - 1);
            double *point = &xyz[index * 3];
            for (int d = 0; d < 3; d++)
                point[d] = coefficient[d][0] +
                    coefficient[d][1] * t +
                    coefficient[d][2] * row_coordinate +
                    frame_normal[d] * value[index];
            if (column == 0)
                memcpy(point,
                       &seed->xyz[(size_t)gap[row].vertex_left * 3],
                       3 * sizeof(*point));
            else if (column + 1 == columns)
                memcpy(point,
                       &seed->xyz[(size_t)gap[row].vertex_right * 3],
                       3 * sizeof(*point));
        }
    }
    if (!adc_patch_geometry_valid(
            xyz, gap, rows, columns, dv, frame_normal, opts)) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_GEOMETRY;
        Arena_restore(arena, mark);
        return 0;
    }
    if (adc_lattice_add_patch(
            lattice, seed, gap, rows, columns, xyz) != 0) {
        patch->reject_reason = ATLAS_DEVELOPABLE_REJECT_PARAMETER;
        Arena_restore(arena, mark);
        return 0;
    }

    size_t vertex_first = mesh->nv;
    for (size_t row = 0; row < rows; row++) {
        int32_t left_vertex = gap[row].vertex_left;
        int32_t right_vertex = gap[row].vertex_right;
        const AtlasRibbonTarget *left =
            &evidence->target[gap[row].target_left];
        const AtlasRibbonTarget *right =
            &evidence->target[gap[row].target_right];
        for (size_t column = 1; column + 1 < columns; column++) {
            double t = (double)column / (double)(columns - 1);
            size_t grid = row * columns + column;
            size_t vertex = mesh->nv++;
            memcpy(&mesh->xyz[vertex * 3], &xyz[grid * 3],
                   3 * sizeof(*mesh->xyz));
            for (int d = 0; d < 2; d++) {
                mesh->local_uv[vertex * 2 + (size_t)d] =
                    (1.0 - t) *
                        seed->local_uv[(size_t)left_vertex * 2 + (size_t)d] +
                    t * seed->local_uv[(size_t)right_vertex * 2 + (size_t)d];
                mesh->global_uv[vertex * 2 + (size_t)d] =
                    (1.0 - t) *
                        seed->global_uv[(size_t)left_vertex * 2 + (size_t)d] +
                    t * seed->global_uv[(size_t)right_vertex * 2 + (size_t)d];
                mesh->uv[vertex * 2 + (size_t)d] =
                    (1.0 - t) *
                        seed->uv[(size_t)left_vertex * 2 + (size_t)d] +
                    t * seed->uv[(size_t)right_vertex * 2 + (size_t)d];
            }
            mesh->source_target[vertex] = -1;
            mesh->vertex_component[vertex] = gap[row].component_left;
            double relative[3] = {
                mesh->xyz[vertex * 3] - evidence->axis_point[0],
                mesh->xyz[vertex * 3 + 1] - evidence->axis_point[1],
                mesh->xyz[vertex * 3 + 2] - evidence->axis_point[2]
            };
            double p[2] = {
                adc_dot3(relative, evidence->basis0),
                adc_dot3(relative, evidence->basis1)
            };
            mesh->vertex_winding[vertex] = adc_interpolated_winding(
                seed, left_vertex, right_vertex, left, right, p, t);
            mesh->vertex_rank[vertex] = t < 0.5
                ? seed->vertex_rank[left_vertex]
                : seed->vertex_rank[right_vertex];
        }
    }
    for (size_t row = 0; row + 1 < rows; row++) {
        for (size_t column = 0; column + 1 < columns; column++) {
            int32_t corner[4];
            size_t grid_row[2] = {row, row + 1};
            size_t grid_column[2] = {column, column + 1};
            for (int vr = 0; vr < 2; vr++) {
                size_t source_row = grid_row[vr];
                for (int uc = 0; uc < 2; uc++) {
                    size_t source_column = grid_column[uc];
                    int32_t vertex = -1;
                    if (source_column == 0)
                        vertex = gap[source_row].vertex_left;
                    else if (source_column + 1 == columns)
                        vertex = gap[source_row].vertex_right;
                    else
                        vertex = (int32_t)(vertex_first +
                            source_row * (columns - 2) +
                            source_column - 1);
                    corner[vr * 2 + uc] = vertex;
                }
            }
            mesh->faces[mesh->nf * 3] = corner[0];
            mesh->faces[mesh->nf * 3 + 1] = corner[1];
            mesh->faces[mesh->nf * 3 + 2] = corner[2];
            mesh->face_component[mesh->nf] = gap[0].component_left;
            mesh->face_inferred[mesh->nf++] = UINT8_C(1);
            mesh->faces[mesh->nf * 3] = corner[1];
            mesh->faces[mesh->nf * 3 + 1] = corner[3];
            mesh->faces[mesh->nf * 3 + 2] = corner[2];
            mesh->face_component[mesh->nf] = gap[0].component_left;
            mesh->face_inferred[mesh->nf++] = UINT8_C(1);
        }
    }
    patch->added_vertices = added_vertices;
    patch->added_faces = added_faces;
    patch->tail_before = stats.tail_before;
    patch->tail_after = stats.tail_after;
    patch->primal_residual = stats.primal;
    patch->dual_residual = stats.dual;
    patch->pcg_iterations = stats.pcg_iterations;
    patch->accepted = 1;
    patch->reject_reason = ATLAS_DEVELOPABLE_ACCEPTED;
    Arena_restore(arena, mark);
    return 1;
}

static int32_t adc_lattice_existing_vertex(
    const AtlasTrackGrowResult *seed,
    const int32_t *head,
    const int32_t *next,
    size_t grid_index,
    int32_t preferred_component)
{
    int32_t fallback = head[grid_index];
    for (int32_t vertex = fallback; vertex >= 0;
         vertex = next[(size_t)vertex]) {
        if (preferred_component >= 0 &&
            seed->vertex_component[(size_t)vertex] == preferred_component)
            return vertex;
    }
    return fallback;
}

static int32_t adc_lattice_corner_vertex(
    AdcLattice *lattice,
    const AtlasTrackGrowResult *seed,
    const int32_t *head,
    const int32_t *next,
    size_t row,
    size_t column)
{
    size_t grid = row * lattice->columns + column;
    AdcLatticeNode *node = adc_lattice_find(
        lattice, adc_lattice_key(row, column), 0);
    int32_t preferred = node != NULL ? node->component : -1;
    int32_t existing = adc_lattice_existing_vertex(
        seed, head, next, grid, preferred);
    if (existing >= 0) return existing;
    if (node == NULL || node->conflict || !(node->weight > 0.0)) return -1;
    return node->vertex;
}

static int adc_seed_grid_init(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AdcLattice *lattice,
    AdcSeedGrid *grid)
{
    if (arena == NULL || evidence == NULL || seed == NULL ||
        lattice == NULL || grid == NULL || lattice->columns == 0 ||
        lattice->rows > SIZE_MAX / lattice->columns)
        return -1;
    memset(grid, 0, sizeof *grid);
    grid->rows = lattice->rows;
    grid->columns = lattice->columns;
    size_t count = grid->rows * grid->columns;
    grid->head = (int32_t *)ARENA_ALLOC(
        arena, (count ? count : 1) * sizeof(*grid->head));
    grid->next = (int32_t *)ARENA_ALLOC(
        arena, (seed->nv ? seed->nv : 1) * sizeof(*grid->next));
    for (size_t i = 0; i < count; i++) grid->head[i] = -1;
    for (size_t vertex = 0; vertex < seed->nv; vertex++) {
        grid->next[vertex] = -1;
        int32_t target_index = seed->source_target[vertex];
        if (target_index < 0 || (size_t)target_index >= evidence->ntarget)
            continue;
        const AtlasRibbonTarget *target = &evidence->target[target_index];
        if (target->row < 0 || target->column < 0 ||
            (size_t)target->row >= grid->rows ||
            (size_t)target->column >= grid->columns)
            continue;
        size_t at = (size_t)target->row * grid->columns +
                    (size_t)target->column;
        grid->next[vertex] = grid->head[at];
        grid->head[at] = (int32_t)vertex;
    }
    return 0;
}

static int32_t adc_seed_grid_vertex(
    const AdcSeedGrid *grid,
    const AtlasTrackGrowResult *seed,
    int32_t row,
    int32_t column,
    int32_t preferred_component)
{
    if (grid == NULL || seed == NULL || row < 0 || column < 0 ||
        (size_t)row >= grid->rows || (size_t)column >= grid->columns)
        return -1;
    size_t at = (size_t)row * grid->columns + (size_t)column;
    return adc_lattice_existing_vertex(
        seed, grid->head, grid->next, at, preferred_component);
}

static int adc_lattice_point_xyz(
    AdcLattice *lattice,
    const AdcSeedGrid *grid,
    const AtlasTrackGrowResult *seed,
    int32_t row,
    int32_t column,
    int32_t preferred_component,
    double xyz[3])
{
    if (lattice == NULL || grid == NULL || seed == NULL || xyz == NULL ||
        row < 0 || column < 0 || (size_t)row >= lattice->rows ||
        (size_t)column >= lattice->columns)
        return 0;
    AdcLatticeNode *node = adc_lattice_find(
        lattice, adc_lattice_key((size_t)row, (size_t)column), 0);
    int32_t preferred = preferred_component;
    if (preferred < 0 && node != NULL) preferred = node->component;
    int32_t vertex = adc_seed_grid_vertex(
        grid, seed, row, column, preferred);
    if (vertex >= 0) {
        memcpy(xyz, &seed->xyz[(size_t)vertex * 3], 3 * sizeof(*xyz));
        return 1;
    }
    if (node == NULL || node->conflict || !(node->weight > 0.0)) return 0;
    for (int d = 0; d < 3; d++) xyz[d] = node->xyz_sum[d] / node->weight;
    return 1;
}

/* Lift one bracketed line from the closest aligned evidence line.  This is a
 * one-dimensional Coons correction: preserve the neighbor's measured curve,
 * then interpolate the two endpoint offsets so the new line lands exactly on
 * both observed banks.  With neighbors on both sides, average their lifts. */
static int adc_aligned_neighbor_prediction(
    AdcLattice *lattice,
    const AdcSeedGrid *grid,
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *gap,
    int32_t coordinate,
    int offset,
    double xyz[3])
{
    int32_t first = gap->column_left;
    int32_t last = gap->column_right;
    if (last <= first || coordinate < first || coordinate > last)
        return 0;
    int32_t row0, column0, row1, column1, row, column;
    if (gap->direction == 0) {
        row0 = row1 = row = gap->row + offset;
        column0 = first;
        column1 = last;
        column = coordinate;
    } else {
        row0 = first;
        row1 = last;
        row = coordinate;
        column0 = column1 = column = gap->row + offset;
    }
    double q0[3], q1[3], q[3];
    if (!adc_lattice_point_xyz(
            lattice, grid, seed, row0, column0,
            gap->component_left, q0) ||
        !adc_lattice_point_xyz(
            lattice, grid, seed, row1, column1,
            gap->component_right, q1) ||
        !adc_lattice_point_xyz(
            lattice, grid, seed, row, column,
            coordinate - first < last - coordinate
                ? gap->component_left : gap->component_right,
            q))
        return 0;
    double t = (double)(coordinate - first) / (double)(last - first);
    const double *a = &seed->xyz[(size_t)gap->vertex_left * 3];
    const double *b = &seed->xyz[(size_t)gap->vertex_right * 3];
    for (int d = 0; d < 3; d++)
        xyz[d] = q[d] + (1.0 - t) * (a[d] - q0[d]) +
                          t * (b[d] - q1[d]);
    return 1;
}

static int adc_aligned_neighbor_segment(
    AdcLattice *lattice,
    const AdcSeedGrid *grid,
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *gap,
    int32_t coordinate,
    int offset)
{
    int32_t row0, column0, row1, column1;
    if (gap->direction == 0) {
        row0 = row1 = gap->row + offset;
        column0 = coordinate;
        column1 = coordinate + 1;
    } else {
        row0 = coordinate;
        row1 = coordinate + 1;
        column0 = column1 = gap->row + offset;
    }
    double xyz[3];
    return adc_lattice_point_xyz(
               lattice, grid, seed, row0, column0, -1, xyz) &&
           adc_lattice_point_xyz(
               lattice, grid, seed, row1, column1, -1, xyz);
}

static void adc_aligned_mark_cell(
    AdcLattice *lattice,
    const AdcRowGap *gap,
    int32_t coordinate,
    int offset)
{
    if (gap->direction == 0)
        adc_lattice_mark_cell(
            lattice, gap->row + (offset < 0 ? -1 : 0), coordinate);
    else
        adc_lattice_mark_cell(
            lattice, coordinate,
            gap->row + (offset < 0 ? -1 : 0));
}

static int adc_lattice_add_aligned_gap(
    AdcLattice *lattice,
    const AdcSeedGrid *grid,
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *gap,
    size_t *out_nodes,
    size_t *out_cells)
{
    *out_nodes = 0;
    *out_cells = 0;
    if (lattice == NULL || grid == NULL || seed == NULL || gap == NULL ||
        gap->column_right <= gap->column_left)
        return 0;
    int64_t span64 = (int64_t)gap->column_right - gap->column_left;
    if (span64 <= 0 || (uint64_t)span64 >= SIZE_MAX) return -1;
    size_t count = (size_t)span64 + 1;
    uint8_t *needed = (uint8_t *)ARENA_CALLOC(
        lattice->arena, count, sizeof(*needed));
    size_t touched = 0;
    for (int offset = -1; offset <= 1; offset += 2) {
        for (int32_t coordinate = gap->column_left;
             coordinate < gap->column_right; coordinate++) {
            if (!adc_aligned_neighbor_segment(
                    lattice, grid, seed, gap, coordinate, offset))
                continue;
            size_t index = (size_t)(coordinate - gap->column_left);
            needed[index] = needed[index + 1] = 1;
            touched++;
        }
    }
    if (touched == 0) return 0;
    size_t additions = 0;
    for (size_t i = 0; i < count; i++) additions += needed[i] != 0;
    if (adc_lattice_reserve_additions(lattice, additions) != 0) return -1;

    size_t proposals = 0;
    for (size_t index = 0; index < count; index++) {
        if (!needed[index]) continue;
        int32_t coordinate = gap->column_left + (int32_t)index;
        int32_t row = gap->direction == 0 ? gap->row : coordinate;
        int32_t column = gap->direction == 0 ? coordinate : gap->row;
        int32_t preferred = index * 2 < count
            ? gap->component_left : gap->component_right;
        if (adc_seed_grid_vertex(
                grid, seed, row, column, preferred) >= 0)
            continue;

        double t = (double)index / (double)(count - 1);
        const double *a = &seed->xyz[(size_t)gap->vertex_left * 3];
        const double *b = &seed->xyz[(size_t)gap->vertex_right * 3];
        double prediction[3] = {
            (1.0 - t) * a[0] + t * b[0],
            (1.0 - t) * a[1] + t * b[1],
            (1.0 - t) * a[2] + t * b[2]
        };
        double lifted[3] = {0.0, 0.0, 0.0};
        int neighbors = 0;
        for (int offset = -1; offset <= 1; offset += 2) {
            double candidate[3];
            if (!adc_aligned_neighbor_prediction(
                    lattice, grid, seed, gap, coordinate,
                    offset, candidate))
                continue;
            for (int d = 0; d < 3; d++) lifted[d] += candidate[d];
            neighbors++;
        }
        if (neighbors > 0) {
            for (int d = 0; d < 3; d++)
                prediction[d] = lifted[d] / (double)neighbors;
        }
        int32_t source = index * 2 < count
            ? gap->vertex_left : gap->vertex_right;
        AdcLatticeNode *before = adc_lattice_find(
            lattice, adc_lattice_key((size_t)row, (size_t)column), 0);
        if (adc_lattice_add_node(
                lattice, row, column, prediction,
                seed->vertex_component[(size_t)source],
                seed->vertex_winding[(size_t)source],
                seed->vertex_rank[(size_t)source]) != 0)
            return -1;
        if (before == NULL) (*out_nodes)++;
        proposals++;
    }
    size_t cells_before = lattice->cell_count;
    for (int offset = -1; offset <= 1; offset += 2) {
        for (int32_t coordinate = gap->column_left;
             coordinate < gap->column_right; coordinate++) {
            if (adc_aligned_neighbor_segment(
                    lattice, grid, seed, gap, coordinate, offset))
                adc_aligned_mark_cell(lattice, gap, coordinate, offset);
        }
    }
    *out_cells = lattice->cell_count - cells_before;
    return proposals > 0 || *out_cells > 0;
}

static int adc_lattice_cell_valid(
    const AtlasTrackGrowResult *mesh,
    const int32_t corner[4],
    const AtlasRibbonObservationSet *evidence,
    const AtlasDevelopableCompleteOptions *opts)
{
    const double *a = &mesh->xyz[(size_t)corner[0] * 3];
    const double *b = &mesh->xyz[(size_t)corner[1] * 3];
    const double *c = &mesh->xyz[(size_t)corner[2] * 3];
    const double *d = &mesh->xyz[(size_t)corner[3] * 3];
    double maximum = opts->chord_max_stretch *
        (evidence->observation_du > evidence->dv
             ? evidence->observation_du : evidence->dv) +
        opts->chord_slack;
    if (adc_distance3(a, b) > maximum || adc_distance3(a, c) > maximum ||
        adc_distance3(b, d) > maximum || adc_distance3(c, d) > maximum)
        return 0;
    double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    double dc[3] = {c[0] - d[0], c[1] - d[1], c[2] - d[2]};
    double db[3] = {b[0] - d[0], b[1] - d[1], b[2] - d[2]};
    double n0[3], n1[3];
    adc_cross3(ab, ac, n0);
    adc_cross3(dc, db, n1);
    double length0 = sqrt(adc_dot3(n0, n0));
    double length1 = sqrt(adc_dot3(n1, n1));
    if (!(length0 > 1.0e-10) || !(length1 > 1.0e-10)) return 0;
    return adc_dot3(n0, n1) / (length0 * length1) >= -0.25;
}

static int adc_direct_cell_valid(
    const AtlasTrackGrowResult *seed,
    const int32_t corner[4],
    const AtlasRibbonObservationSet *evidence,
    const AtlasDevelopableCompleteOptions *opts,
    double *out_score)
{
    if (seed == NULL || corner == NULL || evidence == NULL || opts == NULL)
        return 0;
    for (int i = 0; i < 4; i++)
        if (corner[i] < 0 || (size_t)corner[i] >= seed->nv)
            return 0;
    static const int edge[4][2] = {{0, 1}, {2, 3}, {0, 2}, {1, 3}};
    double du = evidence->observation_du > 0.0
        ? evidence->observation_du : 1.0;
    double dv = evidence->dv > 0.0 ? evidence->dv : 1.0;
    double maximum[4] = {
        opts->chord_max_stretch * du + opts->chord_slack,
        opts->chord_max_stretch * du + opts->chord_slack,
        opts->vertical_stretch * dv + opts->vertical_slack,
        opts->vertical_stretch * dv + opts->vertical_slack
    };
    double score = 0.0;
    for (int i = 0; i < 4; i++) {
        int32_t a = corner[edge[i][0]];
        int32_t b = corner[edge[i][1]];
        double length = adc_distance3(
            &seed->xyz[(size_t)a * 3], &seed->xyz[(size_t)b * 3]);
        if (!(length > 1.0e-8) || length > maximum[i]) return 0;
        double expected = i < 2 ? du : dv;
        score += fabs(length / expected - 1.0);
        if (seed->vertex_winding != NULL) {
            int64_t delta =
                (int64_t)seed->vertex_winding[(size_t)b] -
                (int64_t)seed->vertex_winding[(size_t)a];
            int64_t limit = i < 2 ? 2 : 1;
            if (delta < -limit || delta > limit) return 0;
        }
    }
    if (!adc_lattice_cell_valid(seed, corner, evidence, opts)) return 0;
    if (out_score != NULL) *out_score = score;
    return 1;
}

static void adc_direct_cell_consider(
    const AtlasTrackGrowResult *seed,
    const AtlasRibbonObservationSet *evidence,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *a,
    const AdcRowGap *b,
    int direction,
    int32_t row,
    int32_t column,
    size_t grid,
    int32_t *chosen,
    AdcDirectCell *cell,
    size_t *cell_count)
{
    int32_t corner[4];
    if (direction == 0) {
        corner[0] = a->vertex_left;
        corner[1] = a->vertex_right;
        corner[2] = b->vertex_left;
        corner[3] = b->vertex_right;
    } else {
        corner[0] = a->vertex_left;
        corner[1] = b->vertex_left;
        corner[2] = a->vertex_right;
        corner[3] = b->vertex_right;
    }
    double geometry_score = 0.0;
    if (!adc_direct_cell_valid(
            seed, corner, evidence, opts, &geometry_score))
        return;
    double score = a->score + b->score + geometry_score;
    int32_t at = chosen[grid];
    if (at < 0) {
        if (*cell_count > (size_t)INT32_MAX) return;
        at = (int32_t)(*cell_count);
        chosen[grid] = at;
        (*cell_count)++;
    } else if (!(score < cell[(size_t)at].score)) {
        return;
    }
    AdcDirectCell *out = &cell[(size_t)at];
    out->row = row;
    out->column = column;
    memcpy(out->vertex, corner, sizeof out->vertex);
    out->direction = direction;
    out->score = score;
}

static void adc_short_strip_point(
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *gap,
    double t,
    double point[3])
{
    const double *a = &seed->xyz[(size_t)gap->vertex_left * 3];
    const double *b = &seed->xyz[(size_t)gap->vertex_right * 3];
    for (int d = 0; d < 3; d++)
        point[d] = (1.0 - t) * a[d] + t * b[d];
}

static int adc_short_strip_valid(
    const AtlasTrackGrowResult *seed,
    const AtlasRibbonObservationSet *evidence,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *a,
    const AdcRowGap *b,
    double *out_score)
{
    if (a->direction != b->direction || b->row != a->row + 1 ||
        a->column_left != b->column_left ||
        a->column_right != b->column_right)
        return 0;
    int32_t span = a->column_right - a->column_left;
    int span_limit = a->direction == 0
        ? opts->short_strip_max_span_u : opts->short_strip_max_span_v;
    if (span <= 1 || span > span_limit) return 0;
    if (seed->vertex_winding != NULL) {
        int64_t limit = a->direction == 0 ? 1 : 2;
        int64_t delta0 =
            (int64_t)seed->vertex_winding[(size_t)b->vertex_left] -
            (int64_t)seed->vertex_winding[(size_t)a->vertex_left];
        int64_t delta1 =
            (int64_t)seed->vertex_winding[(size_t)b->vertex_right] -
            (int64_t)seed->vertex_winding[(size_t)a->vertex_right];
        if (delta0 < -limit || delta0 > limit ||
            delta1 < -limit || delta1 > limit)
            return 0;
    }
    double score = a->score + b->score;
    double du = evidence->observation_du > 0.0
        ? evidence->observation_du : 1.0;
    double dv = evidence->dv > 0.0 ? evidence->dv : 1.0;
    double horizontal_maximum =
        opts->chord_max_stretch * du + opts->chord_slack;
    double vertical_maximum =
        opts->vertical_stretch * dv + opts->vertical_slack;
    static const int edge[4][2] = {{0, 1}, {2, 3}, {0, 2}, {1, 3}};
    for (int32_t step = 0; step < span; step++) {
        double t0 = (double)step / (double)span;
        double t1 = (double)(step + 1) / (double)span;
        double xyz[12];
        if (a->direction == 0) {
            adc_short_strip_point(seed, a, t0, &xyz[0]);
            adc_short_strip_point(seed, a, t1, &xyz[3]);
            adc_short_strip_point(seed, b, t0, &xyz[6]);
            adc_short_strip_point(seed, b, t1, &xyz[9]);
        } else {
            adc_short_strip_point(seed, a, t0, &xyz[0]);
            adc_short_strip_point(seed, b, t0, &xyz[3]);
            adc_short_strip_point(seed, a, t1, &xyz[6]);
            adc_short_strip_point(seed, b, t1, &xyz[9]);
        }
        for (int e = 0; e < 4; e++) {
            const double *p = &xyz[(size_t)edge[e][0] * 3];
            const double *q = &xyz[(size_t)edge[e][1] * 3];
            double length = adc_distance3(p, q);
            double maximum = e < 2
                ? horizontal_maximum : vertical_maximum;
            double expected = e < 2 ? du : dv;
            if (!(length > 1.0e-8) || length > maximum) return 0;
            score += fabs(length / expected - 1.0);
        }
        AtlasTrackGrowResult proxy;
        memset(&proxy, 0, sizeof proxy);
        proxy.xyz = xyz;
        static const int32_t corner[4] = {0, 1, 2, 3};
        if (!adc_lattice_cell_valid(
                &proxy, corner, evidence, opts))
            return 0;
    }
    if (out_score != NULL) *out_score = score;
    return 1;
}

static int adc_lattice_add_short_strip(
    AdcLattice *lattice,
    const AtlasTrackGrowResult *seed,
    const AdcRowGap *a,
    const AdcRowGap *b,
    size_t *out_cells)
{
    *out_cells = 0;
    int32_t span = a->column_right - a->column_left;
    if (lattice == NULL || span <= 1 ||
        b->row != a->row + 1 || a->direction != b->direction)
        return 0;
    int needed = 0;
    for (int32_t coordinate = a->column_left;
         coordinate < a->column_right; coordinate++) {
        int32_t row = a->direction == 0 ? a->row : coordinate;
        int32_t column = a->direction == 0 ? coordinate : a->row;
        if (row >= 0 && column >= 0 &&
            (size_t)row + 1 < lattice->rows &&
            (size_t)column + 1 < lattice->columns &&
            !lattice->cell[(size_t)row * lattice->columns +
                           (size_t)column]) {
            needed = 1;
            break;
        }
    }
    if (!needed) return 0;
    size_t additions = ((size_t)span + 1) * 2;
    if (adc_lattice_reserve_additions(lattice, additions) != 0) return -1;
    const AdcRowGap *line[2] = {a, b};
    for (int support = 0; support < 2; support++) {
        const AdcRowGap *g = line[support];
        int32_t winding0 =
            seed->vertex_winding[(size_t)g->vertex_left];
        int32_t winding1 =
            seed->vertex_winding[(size_t)g->vertex_right];
        for (int32_t coordinate = g->column_left;
             coordinate <= g->column_right; coordinate++) {
            double t = (double)(coordinate - g->column_left) /
                       (double)span;
            double point[3];
            adc_short_strip_point(seed, g, t, point);
            int32_t source = t < 0.5
                ? g->vertex_left : g->vertex_right;
            int32_t row = g->direction == 0 ? g->row : coordinate;
            int32_t column = g->direction == 0 ? coordinate : g->row;
            int32_t winding = (int32_t)llround(
                (1.0 - t) * (double)winding0 + t * (double)winding1);
            int32_t component = t < 0.5
                ? g->component_left : g->component_right;
            if (adc_lattice_add_node(
                    lattice, row, column, point, component, winding,
                    seed->vertex_rank[(size_t)source]) != 0)
                return -1;
        }
    }
    size_t before = lattice->cell_count;
    for (int32_t coordinate = a->column_left;
         coordinate < a->column_right; coordinate++) {
        int32_t row = a->direction == 0 ? a->row : coordinate;
        int32_t column = a->direction == 0 ? coordinate : a->row;
        adc_lattice_mark_cell(lattice, row, column);
    }
    *out_cells = lattice->cell_count - before;
    return 1;
}

/*
 * A one-cell gap has no unknown XYZ samples: its four corners are already
 * observed.  Pair such edges by atlas location, not by source component IDs.
 * Component IDs routinely change at cube boundaries, and requiring a stable
 * pair there turns a geometrically trivial square into thousands of isolated
 * one-support runs.  We still retain the endpoint gates from adc_build_gaps()
 * and validate the complete quad in 3D before admitting it.
 */
static int adc_build_direct_cells(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    const AdcRowGap *gap,
    size_t gap_count,
    size_t rows,
    size_t columns,
    AdcLattice *lattice,
    AdcDirectCell **out_cell,
    size_t *out_count,
    size_t *out_u_count,
    size_t *out_v_count,
    size_t *out_short_patches,
    size_t *out_short_cells,
    size_t *out_short_u_cells,
    size_t *out_short_v_cells)
{
    *out_cell = NULL;
    *out_count = *out_u_count = *out_v_count = 0;
    *out_short_patches = *out_short_cells = 0;
    *out_short_u_cells = *out_short_v_cells = 0;
    if (gap_count == 0 || rows < 2 || columns < 2) return 0;
    if (gap_count > (size_t)INT32_MAX || rows > SIZE_MAX / columns)
        return -1;
    size_t grid_size = rows * columns;
    size_t capacity = gap_count < grid_size ? gap_count : grid_size;
    AdcDirectCell *cell = (AdcDirectCell *)ARENA_ALLOC(
        arena, (capacity ? capacity : 1) * sizeof(*cell));
    if (gap_count > SIZE_MAX / (2 * sizeof(AdcShortStrip))) return -1;
    AdcShortStrip *short_strip = (AdcShortStrip *)ARENA_ALLOC(
        arena, (gap_count ? gap_count * 2 : 1) * sizeof(*short_strip));
    Arena_Mark mark = Arena_save(arena);
    int32_t *head = (int32_t *)ARENA_ALLOC(
        arena, grid_size * sizeof(*head));
    int32_t *next = (int32_t *)ARENA_ALLOC(
        arena, gap_count * sizeof(*next));
    int32_t *chosen = (int32_t *)ARENA_ALLOC(
        arena, grid_size * sizeof(*chosen));
    for (size_t i = 0; i < grid_size; i++) {
        head[i] = -1;
        chosen[i] = -1;
    }
    for (size_t i = 0; i < gap_count; i++) next[i] = -1;

    for (size_t i = 0; i < gap_count; i++) {
        if (gap[i].direction != 0 ||
            gap[i].column_right != gap[i].column_left + 1 ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].row >= rows ||
            (size_t)gap[i].column_right >= columns)
            continue;
        size_t key = (size_t)gap[i].row * columns +
                     (size_t)gap[i].column_left;
        next[i] = head[key];
        head[key] = (int32_t)i;
    }
    size_t count = 0;
    for (size_t i = 0; i < gap_count; i++) {
        if (gap[i].direction != 0 ||
            gap[i].column_right != gap[i].column_left + 1 ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].row + 1 >= rows ||
            (size_t)gap[i].column_right >= columns)
            continue;
        size_t grid = (size_t)gap[i].row * columns +
                      (size_t)gap[i].column_left;
        for (int32_t other = head[grid + columns]; other >= 0;
             other = next[(size_t)other])
            adc_direct_cell_consider(
                seed, evidence, opts, &gap[i], &gap[(size_t)other], 0,
                gap[i].row, gap[i].column_left, grid,
                chosen, cell, &count);
    }

    for (size_t i = 0; i < grid_size; i++) head[i] = -1;
    for (size_t i = 0; i < gap_count; i++) next[i] = -1;
    for (size_t i = 0; i < gap_count; i++) {
        if (gap[i].direction != 1 ||
            gap[i].column_right != gap[i].column_left + 1 ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].column_right >= rows ||
            (size_t)gap[i].row >= columns)
            continue;
        size_t key = (size_t)gap[i].column_left * columns +
                     (size_t)gap[i].row;
        next[i] = head[key];
        head[key] = (int32_t)i;
    }
    for (size_t i = 0; i < gap_count; i++) {
        if (gap[i].direction != 1 ||
            gap[i].column_right != gap[i].column_left + 1 ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].column_right >= rows ||
            (size_t)gap[i].row + 1 >= columns)
            continue;
        size_t grid = (size_t)gap[i].column_left * columns +
                      (size_t)gap[i].row;
        for (int32_t other = head[grid + 1]; other >= 0;
             other = next[(size_t)other])
            adc_direct_cell_consider(
                seed, evidence, opts, &gap[i], &gap[(size_t)other], 1,
                gap[i].column_left, gap[i].row, grid,
                chosen, cell, &count);
    }
    for (size_t i = 0; i < count; i++) {
        if (cell[i].direction == 0) (*out_u_count)++;
        else (*out_v_count)++;
    }

    size_t short_count = 0;
    for (size_t i = 0; i < grid_size; i++) {
        head[i] = -1;
        chosen[i] = -1;
    }
    for (size_t i = 0; i < gap_count; i++) next[i] = -1;
    for (size_t i = 0; i < gap_count; i++) {
        int32_t span = gap[i].column_right - gap[i].column_left;
        if (gap[i].direction != 0 || span <= 1 ||
            span > opts->short_strip_max_span_u ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].row >= rows ||
            (size_t)gap[i].column_right >= columns)
            continue;
        size_t key = (size_t)gap[i].row * columns +
                     (size_t)gap[i].column_left;
        next[i] = head[key];
        head[key] = (int32_t)i;
    }
    for (size_t i = 0; i < gap_count; i++) {
        int32_t span = gap[i].column_right - gap[i].column_left;
        if (gap[i].direction != 0 || span <= 1 ||
            span > opts->short_strip_max_span_u ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].row + 1 >= rows ||
            (size_t)gap[i].column_right >= columns)
            continue;
        size_t grid = (size_t)gap[i].row * columns +
                      (size_t)gap[i].column_left;
        for (int32_t other = head[grid + columns]; other >= 0;
             other = next[(size_t)other]) {
            double score = 0.0;
            if (!adc_short_strip_valid(
                    seed, evidence, opts,
                    &gap[i], &gap[(size_t)other], &score))
                continue;
            int32_t at = chosen[grid];
            if (at < 0) {
                if (short_count >= gap_count * 2) return -1;
                at = (int32_t)short_count++;
                chosen[grid] = at;
            } else if (!(score < short_strip[(size_t)at].score)) {
                continue;
            }
            short_strip[(size_t)at].gap_a = (int32_t)i;
            short_strip[(size_t)at].gap_b = other;
            short_strip[(size_t)at].span = span;
            short_strip[(size_t)at].direction = 0;
            short_strip[(size_t)at].score = score;
        }
    }

    for (size_t i = 0; i < grid_size; i++) {
        head[i] = -1;
        chosen[i] = -1;
    }
    for (size_t i = 0; i < gap_count; i++) next[i] = -1;
    for (size_t i = 0; i < gap_count; i++) {
        int32_t span = gap[i].column_right - gap[i].column_left;
        if (gap[i].direction != 1 || span <= 1 ||
            span > opts->short_strip_max_span_v ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].column_right >= rows ||
            (size_t)gap[i].row >= columns)
            continue;
        size_t key = (size_t)gap[i].column_left * columns +
                     (size_t)gap[i].row;
        next[i] = head[key];
        head[key] = (int32_t)i;
    }
    for (size_t i = 0; i < gap_count; i++) {
        int32_t span = gap[i].column_right - gap[i].column_left;
        if (gap[i].direction != 1 || span <= 1 ||
            span > opts->short_strip_max_span_v ||
            gap[i].row < 0 || gap[i].column_left < 0 ||
            (size_t)gap[i].column_right >= rows ||
            (size_t)gap[i].row + 1 >= columns)
            continue;
        size_t grid = (size_t)gap[i].column_left * columns +
                      (size_t)gap[i].row;
        for (int32_t other = head[grid + 1]; other >= 0;
             other = next[(size_t)other]) {
            double score = 0.0;
            if (!adc_short_strip_valid(
                    seed, evidence, opts,
                    &gap[i], &gap[(size_t)other], &score))
                continue;
            int32_t at = chosen[grid];
            if (at < 0) {
                if (short_count >= gap_count * 2) return -1;
                at = (int32_t)short_count++;
                chosen[grid] = at;
            } else if (!(score < short_strip[(size_t)at].score)) {
                continue;
            }
            short_strip[(size_t)at].gap_a = (int32_t)i;
            short_strip[(size_t)at].gap_b = other;
            short_strip[(size_t)at].span = span;
            short_strip[(size_t)at].direction = 1;
            short_strip[(size_t)at].score = score;
        }
    }
    Arena_restore(arena, mark);
    qsort(short_strip, short_count, sizeof(*short_strip),
          adc_compare_short_strip);
    for (size_t i = 0; i < short_count; i++) {
        size_t cells = 0;
        int accepted = adc_lattice_add_short_strip(
            lattice, seed,
            &gap[(size_t)short_strip[i].gap_a],
            &gap[(size_t)short_strip[i].gap_b], &cells);
        if (accepted < 0) return -1;
        if (!accepted || cells == 0) continue;
        (*out_short_patches)++;
        *out_short_cells += cells;
        if (short_strip[i].direction == 0) *out_short_u_cells += cells;
        else *out_short_v_cells += cells;
    }
    *out_cell = cell;
    *out_count = count;
    return 0;
}

static int adc_lattice_rebuild(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    AdcLattice *lattice,
    const AdcDirectCell *direct_cell,
    size_t direct_count,
    AtlasTrackGrowResult *mesh,
    size_t *out_vertices,
    size_t *out_faces,
    size_t *out_direct_cells)
{
    *out_vertices = 0;
    *out_faces = 0;
    *out_direct_cells = 0;
    if (lattice == NULL || lattice->rows == 0 || lattice->columns == 0) {
        return adc_copy_seed(arena, seed, seed->nv, seed->nf, mesh);
    }
    if (lattice->rows > SIZE_MAX / lattice->columns) return -1;
    size_t grid_size = lattice->rows * lattice->columns;
    int32_t *head = (int32_t *)ARENA_ALLOC(
        arena, grid_size * sizeof(*head));
    for (size_t i = 0; i < grid_size; i++) head[i] = -1;
    int32_t *next = (int32_t *)ARENA_ALLOC(
        arena, (seed->nv ? seed->nv : 1) * sizeof(*next));
    for (size_t vertex = 0; vertex < seed->nv; vertex++) {
        next[vertex] = -1;
        int32_t target_index = seed->source_target[vertex];
        if (target_index < 0 || (size_t)target_index >= evidence->ntarget)
            continue;
        const AtlasRibbonTarget *target = &evidence->target[target_index];
        if (target->row < 0 || target->column < 0 ||
            (size_t)target->row >= lattice->rows ||
            (size_t)target->column >= lattice->columns)
            continue;
        size_t grid = (size_t)target->row * lattice->columns +
                      (size_t)target->column;
        next[vertex] = head[grid];
        head[grid] = (int32_t)vertex;
    }

    uint8_t *existing_cell = (uint8_t *)ARENA_CALLOC(
        arena, grid_size, sizeof(*existing_cell));
    for (size_t face = 0; face < seed->nf; face++) {
        int32_t row_min = INT32_MAX, row_max = INT32_MIN;
        int32_t col_min = INT32_MAX, col_max = INT32_MIN;
        int valid = 1;
        for (int corner = 0; corner < 3; corner++) {
            int32_t vertex = seed->faces[face * 3 + (size_t)corner];
            if (vertex < 0 || (size_t)vertex >= seed->nv) {
                valid = 0;
                break;
            }
            int32_t target_index = seed->source_target[(size_t)vertex];
            if (target_index < 0 ||
                (size_t)target_index >= evidence->ntarget) {
                valid = 0;
                break;
            }
            const AtlasRibbonTarget *target = &evidence->target[target_index];
            if (target->row < row_min) row_min = target->row;
            if (target->row > row_max) row_max = target->row;
            if (target->column < col_min) col_min = target->column;
            if (target->column > col_max) col_max = target->column;
        }
        if (valid && row_max == row_min + 1 && col_max == col_min + 1 &&
            row_min >= 0 && col_min >= 0 &&
            (size_t)row_max < lattice->rows &&
            (size_t)col_max < lattice->columns)
            existing_cell[(size_t)row_min * lattice->columns +
                          (size_t)col_min] = 1;
    }

    uint8_t *direct_mask = (uint8_t *)ARENA_CALLOC(
        arena, grid_size, sizeof(*direct_mask));
    for (size_t i = 0; i < direct_count; i++) {
        if (direct_cell[i].row < 0 || direct_cell[i].column < 0 ||
            (size_t)direct_cell[i].row + 1 >= lattice->rows ||
            (size_t)direct_cell[i].column + 1 >= lattice->columns)
            return -1;
        size_t grid = (size_t)direct_cell[i].row * lattice->columns +
                      (size_t)direct_cell[i].column;
        direct_mask[grid] = 1;
    }

    if (lattice->node_count > SIZE_MAX - seed->nv ||
        lattice->cell_count > (SIZE_MAX - seed->nf) / 2 ||
        direct_count > (SIZE_MAX - seed->nf - lattice->cell_count * 2) / 2)
        return -1;
    size_t vertex_capacity = seed->nv + lattice->node_count;
    size_t face_capacity = seed->nf + lattice->cell_count * 2 +
                           direct_count * 2;
    if (vertex_capacity > (size_t)INT32_MAX ||
        adc_copy_seed(
            arena, seed, vertex_capacity, face_capacity, mesh) != 0)
        return -1;
    double offset[2] = {0.0, 0.0};
    if (seed->nv > 0) {
        offset[0] = seed->global_uv[0] - seed->uv[0];
        offset[1] = seed->global_uv[1] - seed->uv[1];
    }
    for (size_t i = 0; i < lattice->node_capacity; i++) {
        AdcLatticeNode *node = &lattice->node[i];
        if (node->key == UINT64_MAX || node->conflict ||
            !(node->weight > 0.0))
            continue;
        size_t row = (size_t)(uint32_t)(node->key >> 32);
        size_t column = (size_t)(uint32_t)node->key;
        size_t grid = row * lattice->columns + column;
        int32_t existing = adc_lattice_existing_vertex(
            seed, head, next, grid, node->component);
        if (existing >= 0) {
            node->vertex = existing;
            continue;
        }
        if (mesh->nv >= vertex_capacity || mesh->nv > (size_t)INT32_MAX)
            return -1;
        size_t vertex = mesh->nv++;
        node->vertex = (int32_t)vertex;
        for (int d = 0; d < 3; d++)
            mesh->xyz[vertex * 3 + (size_t)d] =
                node->xyz_sum[d] / node->weight;
        mesh->global_uv[vertex * 2] = evidence->observation_u0 +
            (double)column * evidence->observation_du;
        mesh->global_uv[vertex * 2 + 1] = evidence->v0 +
            (double)row * evidence->dv;
        for (int d = 0; d < 2; d++) {
            mesh->uv[vertex * 2 + (size_t)d] =
                mesh->global_uv[vertex * 2 + (size_t)d] - offset[d];
            mesh->local_uv[vertex * 2 + (size_t)d] =
                mesh->uv[vertex * 2 + (size_t)d];
        }
        mesh->source_target[vertex] = -1;
        mesh->vertex_component[vertex] = node->component;
        mesh->vertex_winding[vertex] = node->winding;
        mesh->vertex_rank[vertex] = node->rank;
    }

    size_t rejected = 0;
    for (size_t row = 0; row + 1 < lattice->rows; row++) {
        for (size_t column = 0; column + 1 < lattice->columns; column++) {
            size_t grid = row * lattice->columns + column;
            if (!lattice->cell[grid] || existing_cell[grid] ||
                direct_mask[grid])
                continue;
            int32_t corner[4] = {
                adc_lattice_corner_vertex(
                    lattice, seed, head, next, row, column),
                adc_lattice_corner_vertex(
                    lattice, seed, head, next, row, column + 1),
                adc_lattice_corner_vertex(
                    lattice, seed, head, next, row + 1, column),
                adc_lattice_corner_vertex(
                    lattice, seed, head, next, row + 1, column + 1)
            };
            if (corner[0] < 0 || corner[1] < 0 ||
                corner[2] < 0 || corner[3] < 0 ||
                !adc_lattice_cell_valid(mesh, corner, evidence, opts)) {
                rejected++;
                continue;
            }
            if (mesh->nf + 2 > face_capacity) return -1;
            mesh->faces[mesh->nf * 3] = corner[0];
            mesh->faces[mesh->nf * 3 + 1] = corner[1];
            mesh->faces[mesh->nf * 3 + 2] = corner[2];
            mesh->face_component[mesh->nf] =
                mesh->vertex_component[(size_t)corner[0]];
            mesh->face_inferred[mesh->nf++] = UINT8_C(1);
            mesh->faces[mesh->nf * 3] = corner[1];
            mesh->faces[mesh->nf * 3 + 1] = corner[3];
            mesh->faces[mesh->nf * 3 + 2] = corner[2];
            mesh->face_component[mesh->nf] =
                mesh->vertex_component[(size_t)corner[0]];
            mesh->face_inferred[mesh->nf++] = UINT8_C(1);
        }
    }
    size_t emitted_direct = 0;
    for (size_t i = 0; i < direct_count; i++) {
        size_t grid = (size_t)direct_cell[i].row * lattice->columns +
                      (size_t)direct_cell[i].column;
        if (existing_cell[grid] ||
            !adc_direct_cell_valid(
                mesh, direct_cell[i].vertex, evidence, opts, NULL))
            continue;
        const int32_t *corner = direct_cell[i].vertex;
        if (mesh->nf + 2 > face_capacity) return -1;
        mesh->faces[mesh->nf * 3] = corner[0];
        mesh->faces[mesh->nf * 3 + 1] = corner[1];
        mesh->faces[mesh->nf * 3 + 2] = corner[2];
        mesh->face_component[mesh->nf] =
            mesh->vertex_component[(size_t)corner[0]];
        mesh->face_inferred[mesh->nf++] = UINT8_C(1);
        mesh->faces[mesh->nf * 3] = corner[1];
        mesh->faces[mesh->nf * 3 + 1] = corner[3];
        mesh->faces[mesh->nf * 3 + 2] = corner[2];
        mesh->face_component[mesh->nf] =
            mesh->vertex_component[(size_t)corner[0]];
        mesh->face_inferred[mesh->nf++] = UINT8_C(1);
        emitted_direct++;
    }
    lattice->rejected_cells = rejected;
    *out_vertices = mesh->nv - seed->nv;
    *out_faces = mesh->nf - seed->nf;
    *out_direct_cells = emitted_direct;
    return 0;
}

static int adc_remap_components(
    Arena_T arena,
    const AtlasTrackGrowResult *seed,
    AtlasTrackGrowResult *mesh,
    size_t *out_components)
{
    *out_components = 0;
    if (mesh->nv == 0 || mesh->nf == 0) {
        mesh->selected_components = 0;
        return 0;
    }
    if (seed->ncomponents == 0 ||
        seed->ncomponents > (size_t)INT32_MAX ||
        mesh->nv > (size_t)INT32_MAX)
        return -1;

    /*
     * Relabel from the emitted face graph, not inherited source labels.
     * A synthetic lattice vertex records which source bank proposed it, but a
     * rejected neighboring cell can leave that vertex disconnected from the
     * source component.  Unioning those labels therefore invents connectivity
     * that is absent from the OBJ and confuses the winding organizer.
     */
    UnionFind topology = UF_new(arena, (int32_t)mesh->nv);
    uint8_t *used = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nv, sizeof(*used));
    for (size_t face = 0; face < mesh->nf; face++) {
        int32_t vertex[3];
        for (int corner = 0; corner < 3; corner++) {
            vertex[corner] = mesh->faces[face * 3 + (size_t)corner];
            if (vertex[corner] < 0 ||
                (size_t)vertex[corner] >= mesh->nv)
                return -1;
            used[(size_t)vertex[corner]] = 1;
        }
        uf_union(&topology, vertex[0], vertex[1]);
        uf_union(&topology, vertex[0], vertex[2]);
    }
    int32_t *root_id = (int32_t *)ARENA_ALLOC(
        arena, mesh->nv * sizeof(*root_id));
    int32_t *source_id = (int32_t *)ARENA_ALLOC(
        arena, seed->ncomponents * sizeof(*source_id));
    for (size_t i = 0; i < mesh->nv; i++) root_id[i] = -1;
    for (size_t i = 0; i < seed->ncomponents; i++) source_id[i] = -1;
    size_t next = 0;
    for (size_t i = 0; i < mesh->nv; i++) {
        if (!used[i]) continue;
        int32_t source = mesh->vertex_component[i];
        if (source < 0 || (size_t)source >= seed->ncomponents)
            return -1;
        int32_t root = uf_find(&topology, (int32_t)i);
        if (root_id[root] < 0) {
            if (next > (size_t)INT32_MAX) return -1;
            root_id[root] = (int32_t)next++;
        }
        int32_t component = root_id[root];
        if (source_id[source] < 0) source_id[source] = component;
        mesh->vertex_component[i] = component;
    }
    for (size_t i = 0; i < mesh->nf; i++) {
        int32_t root = uf_find(&topology, mesh->faces[i * 3]);
        if (root_id[root] < 0) return -1;
        mesh->face_component[i] = root_id[root];
    }
    /* Unreferenced lattice nodes do not define mesh components.  Retain a
     * valid gid for exporters by attaching them to the first emitted piece
     * carrying the same source label (or piece zero as a final fallback). */
    for (size_t i = 0; i < mesh->nv; i++) {
        if (used[i]) continue;
        int32_t source = mesh->vertex_component[i];
        if (source < 0 || (size_t)source >= seed->ncomponents)
            return -1;
        mesh->vertex_component[i] = source_id[source] >= 0
            ? source_id[source] : 0;
    }
    mesh->selected_components = next;
    *out_components = next;
    return 0;
}

int AtlasDevelopableComplete_build(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    AtlasDevelopableCompleteResult *out)
{
    if (arena == NULL || evidence == NULL || seed == NULL || out == NULL ||
        !adc_options_valid(opts))
        return -1;
    memset(out, 0, sizeof *out);
    out->seed_components = seed->selected_components;
    out->seed_vertices = seed->nv;
    out->seed_faces = seed->nf;

    AdcRowGap *gap = NULL;
    size_t gap_count = 0;
    if (adc_build_gaps(
            arena, evidence, seed, opts, &gap, &gap_count) != 0)
        return -1;
    out->candidate_rows = gap_count;
    AdcRowGap *ordered_gap = NULL;
    AdcRun *run = NULL;
    size_t run_count = 0, pair_count = 0;
    if (adc_build_runs(
            arena, evidence, seed, opts, gap, gap_count,
            &ordered_gap, &run, &run_count, &pair_count) != 0)
        return -1;
    AdcRowGap *candidate_gap = gap;
    size_t candidate_gap_count = gap_count;
    gap = ordered_gap;
    out->candidate_pairs = pair_count;
    size_t attempt_count = run_count;
    if (attempt_count > (size_t)opts->max_patches)
        attempt_count = (size_t)opts->max_patches;

    size_t vertex_capacity = seed->nv;
    size_t face_capacity = seed->nf;
    for (size_t i = 0; i < attempt_count; i++) {
        if (run[i].count < (size_t)opts->min_patch_rows) continue;
        size_t columns = adc_patch_columns(
            evidence, opts, &gap[run[i].first], run[i].count);
        size_t added_vertices = run[i].count * (columns - 2);
        size_t added_faces = (run[i].count - 1) * (columns - 1) * 2;
        if (added_vertices > SIZE_MAX - vertex_capacity ||
            added_faces > SIZE_MAX - face_capacity)
            return -1;
        vertex_capacity += added_vertices;
        face_capacity += added_faces;
    }
    if (vertex_capacity > (size_t)INT32_MAX) return -1;
    if (adc_copy_seed(
            arena, seed, vertex_capacity, face_capacity, &out->mesh) != 0)
        return -1;
    out->patch = (AtlasDevelopablePatch *)ARENA_CALLOC(
        arena, (attempt_count ? attempt_count : 1), sizeof(*out->patch));
    out->npatches = attempt_count;
    AdcLattice lattice;
    memset(&lattice, 0, sizeof lattice);
    AdcLattice *lattice_ptr = NULL;
    if (opts->shared_lattice_fill) {
        if (adc_lattice_init(arena, evidence, opts, &lattice) != 0)
            return -1;
        lattice_ptr = &lattice;
    }
    AdcDirectCell *direct_cell = NULL;
    size_t direct_count = 0, direct_u_count = 0, direct_v_count = 0;
    size_t short_patches = 0, short_cells = 0;
    size_t short_u_cells = 0, short_v_cells = 0;

    for (size_t i = 0; i < attempt_count; i++) {
        AtlasDevelopablePatch *record = &out->patch[i];
        record->id = (int32_t)i;
        const AdcRowGap *bank = &gap[run[i].first];
        size_t rows = run[i].count;
        int32_t row_first, row_last, column_first, column_last;
        adc_patch_grid_bounds(
            bank, rows, &row_first, &row_last,
            &column_first, &column_last);
        if (!opts->shared_lattice_fill &&
            rows >= (size_t)opts->min_patch_rows &&
            adc_patch_overlaps(
                out->patch, i, row_first, row_last,
                column_first, column_last)) {
            record->direction = bank[0].direction;
            record->component_left = bank[0].component_left;
            record->component_right = bank[0].component_right;
            record->row_first = row_first;
            record->row_last = row_last;
            record->column_left_min = column_first;
            record->column_right_max = column_last;
            record->rows = rows;
            record->reject_reason = ATLAS_DEVELOPABLE_REJECT_OVERLAP;
            out->rejected_patches++;
            continue;
        }
        int accepted = adc_complete_run(
            arena, evidence, seed, opts, bank, rows,
            vertex_capacity, face_capacity, lattice_ptr,
            record, &out->mesh);
        if (accepted < 0) return -1;
        if (accepted) {
            out->accepted_patches++;
            out->added_vertices += record->added_vertices;
            out->added_faces += record->added_faces;
            out->admm_iterations += (size_t)(
                opts->nuclear_iterations + opts->rank_tail_iterations);
            out->pcg_iterations += (size_t)record->pcg_iterations;
            out->tail_before += record->tail_before;
            out->tail_after += record->tail_after;
            if (record->primal_residual > out->primal_residual_max)
                out->primal_residual_max = record->primal_residual;
            if (record->dual_residual > out->dual_residual_max)
                out->dual_residual_max = record->dual_residual;
        } else {
            out->rejected_patches++;
        }
    }
    size_t fitted_accepted = out->accepted_patches;
    if (opts->fill_aligned_gaps) {
        AdcSeedGrid seed_grid;
        if (lattice_ptr == NULL ||
            adc_seed_grid_init(
                arena, evidence, seed, lattice_ptr, &seed_grid) != 0)
            return -1;
        /* A line can become usable after its neighbor is filled, so allow a
         * few deterministic closure waves.  Every admitted line remains
         * directly bracketed by evidence at both ends. */
        for (int round = 0; round < 8; round++) {
            size_t progress = 0;
            for (size_t i = 0; i < attempt_count; i++) {
                AtlasDevelopablePatch *record = &out->patch[i];
                if (record->reject_reason != ATLAS_DEVELOPABLE_REJECT_SHORT ||
                    run[i].count != 1)
                    continue;
                const AdcRowGap *line = &gap[run[i].first];
                size_t nodes = 0, cells = 0;
                int accepted = adc_lattice_add_aligned_gap(
                    lattice_ptr, &seed_grid, seed, line, &nodes, &cells);
                if (accepted < 0) return -1;
                if (!accepted) continue;
                record->accepted = 1;
                record->reject_reason = ATLAS_DEVELOPABLE_ACCEPTED;
                record->columns =
                    (size_t)(line->column_right - line->column_left) + 1;
                record->added_vertices = nodes;
                record->added_faces = cells * 2;
                out->accepted_patches++;
                if (out->rejected_patches == 0) return -1;
                out->rejected_patches--;
                out->aligned_gap_lines++;
                out->aligned_gap_cells += cells;
                out->added_vertices += nodes;
                out->added_faces += cells * 2;
                progress++;
            }
            if (progress == 0) break;
        }
        if (adc_build_direct_cells(
                arena, evidence, seed, opts,
                candidate_gap, candidate_gap_count,
                lattice.rows, lattice.columns, &lattice,
                &direct_cell, &direct_count,
                &direct_u_count, &direct_v_count,
                &short_patches, &short_cells,
                &short_u_cells, &short_v_cells) != 0)
            return -1;
    }
    if (fitted_accepted > 0) {
        out->tail_before /= (double)fitted_accepted;
        out->tail_after /= (double)fitted_accepted;
    }
    if (opts->shared_lattice_fill) {
        AtlasTrackGrowResult shared_mesh;
        size_t shared_vertices = 0, shared_faces = 0, direct_emitted = 0;
        if (adc_lattice_rebuild(
                arena, evidence, seed, opts, &lattice,
                direct_cell, direct_count,
                &shared_mesh, &shared_vertices, &shared_faces,
                &direct_emitted) != 0)
            return -1;
        if (direct_emitted != direct_count) return -1;
        out->mesh = shared_mesh;
        out->added_vertices = shared_vertices;
        out->added_faces = shared_faces;
        out->lattice_nodes = shared_vertices;
        out->lattice_cells = shared_faces / 2;
        out->lattice_conflicts = lattice.conflicts;
        out->lattice_rejected_cells = lattice.rejected_cells;
        out->adjacent_gap_cells = direct_emitted;
        out->adjacent_gap_u_cells = direct_u_count;
        out->adjacent_gap_v_cells = direct_v_count;
        out->short_strip_patches = short_patches;
        out->short_strip_cells = short_cells;
        out->short_strip_u_cells = short_u_cells;
        out->short_strip_v_cells = short_v_cells;
    }
    if (out->mesh.nv != seed->nv + out->added_vertices ||
        out->mesh.nf != seed->nf + out->added_faces)
        return -1;
    if (adc_remap_components(
            arena, seed, &out->mesh, &out->completed_components) != 0)
        return -1;
    return 0;
}

static int adc_selftest_adjoint(Arena_T arena)
{
    size_t rows = 7, columns = 8;
    size_t nodes = rows * columns;
    size_t coefficients = adc_hessian_count(rows, columns) * 3;
    double *x = (double *)ARENA_ALLOC(arena, nodes * sizeof(*x));
    double *y = (double *)ARENA_ALLOC(
        arena, coefficients * sizeof(*y));
    double *hx = (double *)ARENA_ALLOC(
        arena, coefficients * sizeof(*hx));
    double *hty = (double *)ARENA_ALLOC(arena, nodes * sizeof(*hty));
    for (size_t i = 0; i < nodes; i++)
        x[i] = sin(0.37 * (double)i) + 0.1 * (double)i;
    for (size_t i = 0; i < coefficients; i++)
        y[i] = cos(0.19 * (double)i) - 0.03 * (double)i;
    adc_hessian_apply(x, rows, columns, 0.25, 0.0625, 0.03125, hx);
    adc_hessian_adjoint(
        y, rows, columns, 0.25, 0.0625, 0.03125, hty);
    double left = adc_dot(hx, y, coefficients);
    double right = adc_dot(x, hty, nodes);
    double error = fabs(left - right) /
        (1.0 + fabs(left) + fabs(right));
    if (error > 1.0e-12) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL adjoint "
                "left=%.17g right=%.17g error=%.3g\n",
                left, right, error);
        return 1;
    }
    return 0;
}

static int adc_selftest_admm(Arena_T arena)
{
    size_t rows = 11, columns = 13;
    size_t nodes = rows * columns;
    double *prior = (double *)ARENA_ALLOC(arena, nodes * sizeof(*prior));
    double *weight = (double *)ARENA_ALLOC(arena, nodes * sizeof(*weight));
    double *value = (double *)ARENA_ALLOC(arena, nodes * sizeof(*value));
    for (size_t row = 0; row < rows; row++) {
        double y = (double)row - 5.0;
        for (size_t column = 0; column < columns; column++) {
            double x = (double)column - 6.0;
            size_t index = row * columns + column;
            prior[index] = 0.08 * x * x + 0.025 * y * y;
            weight[index] = row == 0 || row + 1 == rows ||
                            column == 0 || column + 1 == columns
                ? 64.0 : 0.2;
        }
    }
    AtlasDevelopableCompleteOptions opts;
    AtlasDevelopableCompleteOptions_default(&opts);
    opts.nuclear_iterations = 4;
    opts.rank_tail_iterations = 12;
    opts.pcg_iterations = 200;
    AdcAdmmStats stats;
    if (adc_admm_solve(
            arena, prior, weight, rows, columns, 1.0, 1.0,
            &opts, value, &stats) != 0 ||
        !isfinite(stats.tail_after) ||
        stats.tail_after > stats.tail_before + 1.0e-6) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL ADMM "
                "tail %.6g -> %.6g primal=%.3g dual=%.3g\n",
                stats.tail_before, stats.tail_after,
                stats.primal, stats.dual);
        return 1;
    }
    return 0;
}

static int adc_selftest_bridge(Arena_T arena)
{
    enum { ROWS = 6, PER_ROW = 4, NV = ROWS * PER_ROW };
    AtlasRibbonObservationSet evidence;
    memset(&evidence, 0, sizeof evidence);
    evidence.ntarget = NV;
    evidence.nrows = ROWS;
    evidence.observation_du = 2.0;
    evidence.dv = 2.0;
    evidence.axis[0] = 1.0;
    evidence.basis0[1] = 1.0;
    evidence.basis1[2] = 1.0;
    evidence.target = (AtlasRibbonTarget *)ARENA_CALLOC(
        arena, NV, sizeof(*evidence.target));

    AtlasTrackGrowResult seed;
    memset(&seed, 0, sizeof seed);
    seed.nv = NV;
    seed.nf = (ROWS - 1) * 4;
    seed.ncomponents = 2;
    seed.selected_components = 2;
    seed.winding_direction = 1;
    seed.xyz = (double *)ARENA_ALLOC(
        arena, NV * 3 * sizeof(*seed.xyz));
    seed.local_uv = (double *)ARENA_ALLOC(
        arena, NV * 2 * sizeof(*seed.local_uv));
    seed.global_uv = (double *)ARENA_ALLOC(
        arena, NV * 2 * sizeof(*seed.global_uv));
    seed.uv = (double *)ARENA_ALLOC(
        arena, NV * 2 * sizeof(*seed.uv));
    seed.source_target = (int32_t *)ARENA_ALLOC(
        arena, NV * sizeof(*seed.source_target));
    seed.vertex_component = (int32_t *)ARENA_ALLOC(
        arena, NV * sizeof(*seed.vertex_component));
    seed.vertex_winding = (int32_t *)ARENA_CALLOC(
        arena, NV, sizeof(*seed.vertex_winding));
    seed.vertex_rank = (uint64_t *)ARENA_ALLOC(
        arena, NV * sizeof(*seed.vertex_rank));
    seed.faces = (int32_t *)ARENA_ALLOC(
        arena, seed.nf * 3 * sizeof(*seed.faces));
    seed.face_component = (int32_t *)ARENA_ALLOC(
        arena, seed.nf * sizeof(*seed.face_component));
    seed.face_inferred = (uint8_t *)ARENA_CALLOC(
        arena, seed.nf, sizeof(*seed.face_inferred));
    seed.component = (AtlasTrackGrowComponent *)ARENA_CALLOC(
        arena, 2, sizeof(*seed.component));

    static const int column[PER_ROW] = {0, 1, 5, 6};
    for (int row = 0; row < ROWS; row++) {
        for (int local = 0; local < PER_ROW; local++) {
            int index = row * PER_ROW + local;
            int component = local < 2 ? 0 : 1;
            double u = 2.0 * (double)column[local];
            double v = 2.0 * (double)row;
            AtlasRibbonTarget *target = &evidence.target[index];
            target->row = row;
            target->column = column[local];
            target->chart0 = component;
            target->chart1 = -1;
            target->accepted = 1;
            target->u = u;
            target->v = v;
            target->p[0] = u;
            target->tangent[0] = 1.0;
            seed.xyz[(size_t)index * 3] = v;
            seed.xyz[(size_t)index * 3 + 1] = u;
            seed.xyz[(size_t)index * 3 + 2] = 0.0;
            seed.local_uv[(size_t)index * 2] = u;
            seed.local_uv[(size_t)index * 2 + 1] = v;
            seed.global_uv[(size_t)index * 2] = u;
            seed.global_uv[(size_t)index * 2 + 1] = v;
            seed.uv[(size_t)index * 2] = u;
            seed.uv[(size_t)index * 2 + 1] = v;
            seed.source_target[index] = index;
            seed.vertex_component[index] = component;
            seed.vertex_rank[index] = (uint64_t)component;
        }
    }
    size_t face = 0;
    for (int component = 0; component < 2; component++) {
        int offset = component * 2;
        for (int row = 0; row + 1 < ROWS; row++) {
            int32_t a = row * PER_ROW + offset;
            int32_t b = a + 1;
            int32_t c = (row + 1) * PER_ROW + offset;
            int32_t d = c + 1;
            seed.faces[face * 3] = a;
            seed.faces[face * 3 + 1] = b;
            seed.faces[face * 3 + 2] = c;
            seed.face_component[face++] = component;
            seed.faces[face * 3] = b;
            seed.faces[face * 3 + 1] = d;
            seed.faces[face * 3 + 2] = c;
            seed.face_component[face++] = component;
        }
    }
    AtlasDevelopableCompleteOptions opts;
    AtlasDevelopableCompleteOptions_default(&opts);
    opts.max_gap_u = 20.0;
    opts.sample_columns = 5;
    opts.nuclear_iterations = 2;
    opts.rank_tail_iterations = 4;
    AtlasDevelopableCompleteResult result;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &result) != 0 ||
        result.accepted_patches != 1 || result.added_vertices != 18 ||
        result.added_faces != 40 || result.completed_components != 1 ||
        result.mesh.nv != NV + 18 || result.mesh.nf != seed.nf + 40) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL bridge "
                "patch=%zu vertices=%zu faces=%zu components=%zu "
                "mesh=%zu/%zu\n",
                result.accepted_patches, result.added_vertices,
                result.added_faces, result.completed_components,
                result.mesh.nv, result.mesh.nf);
        return 1;
    }

    for (size_t vertex = 0; vertex < seed.nv; vertex++)
        seed.vertex_component[vertex] = 0;
    for (size_t triangle = 0; triangle < seed.nf; triangle++)
        seed.face_component[triangle] = 0;
    seed.ncomponents = 1;
    seed.selected_components = 1;
    AtlasDevelopableCompleteResult internal_off;
        if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &internal_off) != 0 ||
        internal_off.accepted_patches != 0 ||
        internal_off.added_vertices != 0 ||
        internal_off.completed_components != 2) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL internal-off "
                "patch=%zu vertices=%zu components=%zu\n",
                internal_off.accepted_patches,
                internal_off.added_vertices,
                internal_off.completed_components);
        return 1;
    }
    opts.fill_internal_gaps = 1;
    AtlasDevelopableCompleteResult internal_on;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &internal_on) != 0 ||
        internal_on.accepted_patches != 1 ||
        internal_on.added_vertices != 18 ||
        internal_on.added_faces != 40 ||
        internal_on.completed_components != 1) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL internal-on "
                "patch=%zu vertices=%zu faces=%zu components=%zu\n",
                internal_on.accepted_patches,
                internal_on.added_vertices,
                internal_on.added_faces,
                internal_on.completed_components);
        return 1;
    }

    for (int row = 0; row < ROWS; row++) {
        for (int local = 0; local < PER_ROW; local++) {
            int index = row * PER_ROW + local;
            seed.vertex_component[index] =
                local < 2 ? row * 2 : row * 2 + 1;
        }
    }
    for (size_t triangle = 0; triangle < seed.nf; triangle++)
        seed.face_component[triangle] = 0;
    seed.ncomponents = ROWS * 2;
    seed.selected_components = ROWS * 2;
    opts.fill_internal_gaps = 0;
    opts.allow_component_drift = 0;
    AtlasDevelopableCompleteResult drift_off;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &drift_off) != 0 ||
        drift_off.accepted_patches != 0 ||
        drift_off.added_vertices != 0) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL drift-off "
                "patch=%zu vertices=%zu\n",
                drift_off.accepted_patches, drift_off.added_vertices);
        return 1;
    }
    opts.allow_component_drift = 1;
    AtlasDevelopableCompleteResult drift_on;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &drift_on) != 0 ||
        drift_on.accepted_patches != 1 ||
        drift_on.added_vertices != 18 ||
        drift_on.added_faces != 40 ||
        drift_on.completed_components != 1 ||
        drift_on.patch[0].component_pair_changes != ROWS - 1) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL drift-on "
                "patch=%zu vertices=%zu faces=%zu components=%zu "
                "changes=%zu\n",
                drift_on.accepted_patches,
                drift_on.added_vertices,
                drift_on.added_faces,
                drift_on.completed_components,
                drift_on.npatches > 0
                    ? drift_on.patch[0].component_pair_changes : 0);
        return 1;
    }

    evidence.nrows = 7;
    for (int support = 0; support < ROWS; support++) {
        for (int local = 0; local < PER_ROW; local++) {
            int index = support * PER_ROW + local;
            int component = local < 2 ? 0 : 1;
            double u = 2.0 * (double)support;
            double v = 2.0 * (double)column[local];
            AtlasRibbonTarget *target = &evidence.target[index];
            target->row = column[local];
            target->column = support;
            target->chart0 = component;
            target->u = u;
            target->v = v;
            target->p[0] = u;
            target->p[1] = 0.0;
            target->tangent[0] = 1.0;
            target->tangent[1] = 0.0;
            seed.xyz[(size_t)index * 3] = v;
            seed.xyz[(size_t)index * 3 + 1] = u;
            seed.xyz[(size_t)index * 3 + 2] = 0.0;
            seed.local_uv[(size_t)index * 2] = u;
            seed.local_uv[(size_t)index * 2 + 1] = v;
            seed.global_uv[(size_t)index * 2] = u;
            seed.global_uv[(size_t)index * 2 + 1] = v;
            seed.uv[(size_t)index * 2] = u;
            seed.uv[(size_t)index * 2 + 1] = v;
            seed.vertex_component[index] = component;
            seed.vertex_rank[index] = (uint64_t)component;
        }
    }
    for (size_t triangle = 0; triangle < seed.nf; triangle++)
        seed.face_component[triangle] =
            triangle < seed.nf / 2 ? 0 : 1;
    seed.ncomponents = 2;
    seed.selected_components = 2;
    AtlasDevelopableCompleteOptions_default(&opts);
    opts.max_gap_v = 20.0;
    opts.sample_columns = 5;
    opts.nuclear_iterations = 2;
    opts.rank_tail_iterations = 4;
    opts.fill_vertical_gaps = 1;
    AtlasDevelopableCompleteResult vertical;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &vertical) != 0 ||
        vertical.accepted_patches != 1 ||
        vertical.added_vertices != 18 ||
        vertical.added_faces != 40 ||
        vertical.completed_components != 1 ||
        vertical.patch[0].direction != 1) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL vertical "
                "patch=%zu vertices=%zu faces=%zu components=%zu "
                "direction=%d\n",
                vertical.accepted_patches,
                vertical.added_vertices,
                vertical.added_faces,
                vertical.completed_components,
                vertical.npatches > 0 ? vertical.patch[0].direction : -1);
        return 1;
    }

    /* Two spatially disjoint V holes may have the same component pair at
     * every U support sample.  They must remain two independent runs; a
     * component-pair-only uniqueness key used to discard one of them. */
    static const int repeated_v_row[PER_ROW] = {0, 2, 3, 5};
    for (int support = 0; support < ROWS; support++) {
        for (int local = 0; local < PER_ROW; local++) {
            int index = support * PER_ROW + local;
            double u = 2.0 * (double)support;
            double v = 2.0 * (double)repeated_v_row[local];
            AtlasRibbonTarget *target = &evidence.target[index];
            target->row = repeated_v_row[local];
            target->column = support;
            target->chart0 = 0;
            target->u = u;
            target->v = v;
            target->p[0] = u;
            target->p[1] = 0.0;
            target->tangent[0] = 1.0;
            target->tangent[1] = 0.0;
            seed.xyz[(size_t)index * 3] = v;
            seed.xyz[(size_t)index * 3 + 1] = u;
            seed.xyz[(size_t)index * 3 + 2] = 0.0;
            seed.local_uv[(size_t)index * 2] = u;
            seed.local_uv[(size_t)index * 2 + 1] = v;
            seed.global_uv[(size_t)index * 2] = u;
            seed.global_uv[(size_t)index * 2 + 1] = v;
            seed.uv[(size_t)index * 2] = u;
            seed.uv[(size_t)index * 2 + 1] = v;
            seed.vertex_component[index] = 0;
            seed.vertex_rank[index] = 0;
        }
    }
    seed.nf = 0;
    seed.ncomponents = 1;
    seed.selected_components = 1;
    AtlasDevelopableCompleteOptions_default(&opts);
    opts.max_gap_v = 20.0;
    opts.sample_columns = 5;
    opts.nuclear_iterations = 2;
    opts.rank_tail_iterations = 4;
    opts.fill_vertical_gaps = 1;
    opts.fill_internal_gaps = 1;
    AtlasDevelopableCompleteResult repeated_vertical;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &repeated_vertical) != 0 ||
        repeated_vertical.accepted_patches != 2 ||
        repeated_vertical.added_vertices != 36 ||
        repeated_vertical.added_faces != 80) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL repeated-v "
                "patch=%zu vertices=%zu faces=%zu\n",
                repeated_vertical.accepted_patches,
                repeated_vertical.added_vertices,
                repeated_vertical.added_faces);
        return 1;
    }
    opts.shared_lattice_fill = 1;
    AtlasDevelopableCompleteResult repeated_vertical_lattice;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts,
            &repeated_vertical_lattice) != 0 ||
        repeated_vertical_lattice.accepted_patches != 2 ||
        repeated_vertical_lattice.added_vertices != 12 ||
        repeated_vertical_lattice.added_faces != 40 ||
        repeated_vertical_lattice.lattice_cells != 20 ||
        repeated_vertical_lattice.lattice_conflicts != 0) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL repeated-v-grid "
                "patch=%zu vertices=%zu faces=%zu cells=%zu conflicts=%zu\n",
                repeated_vertical_lattice.accepted_patches,
                repeated_vertical_lattice.added_vertices,
                repeated_vertical_lattice.added_faces,
                repeated_vertical_lattice.lattice_cells,
                repeated_vertical_lattice.lattice_conflicts);
        return 1;
    }

    /* The identical ambiguity exists in U: one component pair may bound more
     * than one disjoint horizontal hole in every V evidence row. */
    static const int repeated_u_column[PER_ROW] = {0, 2, 3, 5};
    for (int row = 0; row < ROWS; row++) {
        for (int local = 0; local < PER_ROW; local++) {
            int index = row * PER_ROW + local;
            double u = 2.0 * (double)repeated_u_column[local];
            double v = 2.0 * (double)row;
            AtlasRibbonTarget *target = &evidence.target[index];
            target->row = row;
            target->column = repeated_u_column[local];
            target->chart0 = 0;
            target->u = u;
            target->v = v;
            target->p[0] = u;
            target->p[1] = 0.0;
            target->tangent[0] = 1.0;
            target->tangent[1] = 0.0;
            seed.xyz[(size_t)index * 3] = v;
            seed.xyz[(size_t)index * 3 + 1] = u;
            seed.xyz[(size_t)index * 3 + 2] = 0.0;
            seed.local_uv[(size_t)index * 2] = u;
            seed.local_uv[(size_t)index * 2 + 1] = v;
            seed.global_uv[(size_t)index * 2] = u;
            seed.global_uv[(size_t)index * 2 + 1] = v;
            seed.uv[(size_t)index * 2] = u;
            seed.uv[(size_t)index * 2 + 1] = v;
            seed.vertex_component[index] = 0;
            seed.vertex_rank[index] = 0;
        }
    }
    AtlasDevelopableCompleteOptions_default(&opts);
    opts.max_gap_u = 20.0;
    opts.sample_columns = 5;
    opts.nuclear_iterations = 2;
    opts.rank_tail_iterations = 4;
    opts.fill_internal_gaps = 1;
    AtlasDevelopableCompleteResult repeated_horizontal;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts, &repeated_horizontal) != 0 ||
        repeated_horizontal.accepted_patches != 2 ||
        repeated_horizontal.added_vertices != 36 ||
        repeated_horizontal.added_faces != 80) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL repeated-u "
                "patch=%zu vertices=%zu faces=%zu\n",
                repeated_horizontal.accepted_patches,
                repeated_horizontal.added_vertices,
                repeated_horizontal.added_faces);
        return 1;
    }
    opts.shared_lattice_fill = 1;
    AtlasDevelopableCompleteResult repeated_horizontal_lattice;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &opts,
            &repeated_horizontal_lattice) != 0 ||
        repeated_horizontal_lattice.accepted_patches != 2 ||
        repeated_horizontal_lattice.added_vertices != 12 ||
        repeated_horizontal_lattice.added_faces != 40 ||
        repeated_horizontal_lattice.lattice_cells != 20 ||
        repeated_horizontal_lattice.lattice_conflicts != 0) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL repeated-u-grid "
                "patch=%zu vertices=%zu faces=%zu cells=%zu conflicts=%zu\n",
                repeated_horizontal_lattice.accepted_patches,
                repeated_horizontal_lattice.added_vertices,
                repeated_horizontal_lattice.added_faces,
                repeated_horizontal_lattice.lattice_cells,
                repeated_horizontal_lattice.lattice_conflicts);
        return 1;
    }
    return 0;
}

static int adc_selftest_lattice_cross(Arena_T arena)
{
    AtlasRibbonObservationSet evidence;
    memset(&evidence, 0, sizeof evidence);
    evidence.nrows = 6;
    evidence.observation_du = 2.0;
    evidence.dv = 2.0;
    evidence.ntarget = 1;
    evidence.target = (AtlasRibbonTarget *)ARENA_CALLOC(
        arena, 1, sizeof(*evidence.target));
    evidence.target[0].column = 5;

    AtlasDevelopableCompleteOptions opts;
    AtlasDevelopableCompleteOptions_default(&opts);
    AdcLattice lattice;
    if (adc_lattice_init(arena, &evidence, &opts, &lattice) != 0)
        return 1;
    for (int32_t row = 2; row <= 3; row++) {
        for (int32_t column = 0; column < 6; column++) {
            double xyz[3] = {
                2.0 * (double)row, 2.0 * (double)column, 0.0
            };
            if (adc_lattice_add_node(
                    &lattice, row, column, xyz, 0, 0, 0) != 0)
                return 1;
        }
    }
    for (int32_t column = 2; column <= 3; column++) {
        for (int32_t row = 0; row < 6; row++) {
            double xyz[3] = {
                2.0 * (double)row, 2.0 * (double)column, 0.0
            };
            if (adc_lattice_add_node(
                    &lattice, row, column, xyz, 0, 0, 0) != 0)
                return 1;
        }
    }
    for (int32_t column = 0; column < 5; column++)
        adc_lattice_mark_cell(&lattice, 2, column);
    for (int32_t row = 0; row < 5; row++)
        adc_lattice_mark_cell(&lattice, row, 2);

    AtlasTrackGrowResult seed;
    AtlasTrackGrowResult mesh;
    memset(&seed, 0, sizeof seed);
    seed.ncomponents = 1;
    seed.selected_components = 1;
    size_t added_vertices = 0, added_faces = 0, direct_cells = 0;
    if (adc_lattice_rebuild(
            arena, &evidence, &seed, &opts, &lattice,
            NULL, 0, &mesh, &added_vertices, &added_faces,
            &direct_cells) != 0 ||
        lattice.node_count != 20 || lattice.cell_count != 9 ||
        lattice.conflicts != 0 || lattice.rejected_cells != 0 ||
        direct_cells != 0 || added_vertices != 20 || added_faces != 18 ||
        mesh.nv != 20 || mesh.nf != 18) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL lattice-cross "
                "nodes=%zu cells=%zu conflicts=%zu rejected=%zu "
                "mesh=%zu/%zu\n",
                lattice.node_count, lattice.cell_count, lattice.conflicts,
                lattice.rejected_cells, mesh.nv, mesh.nf);
        return 1;
    }
    return 0;
}

static int adc_selftest_aligned_gap_direction(
    Arena_T arena, int direction)
{
    enum { N = 7 };
    AtlasRibbonObservationSet evidence;
    AtlasTrackGrowResult seed;
    AtlasDevelopableCompleteOptions opts;
    memset(&evidence, 0, sizeof evidence);
    memset(&seed, 0, sizeof seed);
    AtlasDevelopableCompleteOptions_default(&opts);
    evidence.nrows = direction == 0 ? 2 : 5;
    evidence.observation_du = 2.0;
    evidence.dv = 2.0;
    evidence.ntarget = N;
    evidence.target = (AtlasRibbonTarget *)ARENA_CALLOC(
        arena, N, sizeof(*evidence.target));
    seed.nv = N;
    seed.ncomponents = seed.selected_components = 1;
    seed.xyz = (double *)ARENA_ALLOC(
        arena, N * 3 * sizeof(*seed.xyz));
    seed.source_target = (int32_t *)ARENA_ALLOC(
        arena, N * sizeof(*seed.source_target));
    seed.vertex_component = (int32_t *)ARENA_CALLOC(
        arena, N, sizeof(*seed.vertex_component));
    seed.vertex_winding = (int32_t *)ARENA_CALLOC(
        arena, N, sizeof(*seed.vertex_winding));
    seed.vertex_rank = (uint64_t *)ARENA_CALLOC(
        arena, N, sizeof(*seed.vertex_rank));
    for (int i = 0; i < 5; i++) {
        evidence.target[i].row = direction == 0 ? 0 : i;
        evidence.target[i].column = direction == 0 ? i : 0;
    }
    evidence.target[5].row = 0;
    evidence.target[5].column = direction == 0 ? 0 : 1;
    evidence.target[6].row = direction == 0 ? 1 : 4;
    evidence.target[6].column = direction == 0 ? 4 : 1;
    if (direction == 0) evidence.target[5].row = 1;
    for (int i = 0; i < N; i++) {
        seed.source_target[i] = i;
        seed.xyz[(size_t)i * 3] =
            2.0 * (double)evidence.target[i].row;
        seed.xyz[(size_t)i * 3 + 1] =
            2.0 * (double)evidence.target[i].column;
        seed.xyz[(size_t)i * 3 + 2] = 0.0;
    }

    AdcLattice lattice;
    AdcSeedGrid grid;
    if (adc_lattice_init(arena, &evidence, &opts, &lattice) != 0 ||
        adc_seed_grid_init(arena, &evidence, &seed, &lattice, &grid) != 0)
        return 1;
    AdcRowGap gap;
    memset(&gap, 0, sizeof gap);
    gap.direction = direction;
    gap.row = 1;
    gap.column_left = 0;
    gap.column_right = 4;
    gap.target_left = 5;
    gap.target_right = 6;
    gap.vertex_left = 5;
    gap.vertex_right = 6;
    size_t nodes = 0, cells = 0;
    int accepted = adc_lattice_add_aligned_gap(
        &lattice, &grid, &seed, &gap, &nodes, &cells);
    int failed = accepted != 1 || nodes != 3 || cells != 4 ||
                 lattice.node_count != 3 || lattice.cell_count != 4;
    for (int coordinate = 1; coordinate < 4 && !failed; coordinate++) {
        int32_t row = direction == 0 ? 1 : coordinate;
        int32_t column = direction == 0 ? coordinate : 1;
        AdcLatticeNode *node = adc_lattice_find(
            &lattice, adc_lattice_key((size_t)row, (size_t)column), 0);
        if (node == NULL || node->conflict || !(node->weight > 0.0) ||
            fabs(node->xyz_sum[0] / node->weight - 2.0 * row) > 1.0e-12 ||
            fabs(node->xyz_sum[1] / node->weight - 2.0 * column) > 1.0e-12)
            failed = 1;
    }
    if (failed) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL aligned-%c "
                "accepted=%d nodes=%zu/%zu cells=%zu/%zu\n",
                direction == 0 ? 'u' : 'v', accepted, nodes,
                lattice.node_count, cells, lattice.cell_count);
        return 1;
    }
    return 0;
}

static int adc_selftest_aligned_gaps(Arena_T arena)
{
    if (adc_selftest_aligned_gap_direction(arena, 0) != 0) return 1;
    return adc_selftest_aligned_gap_direction(arena, 1);
}

static int adc_selftest_existing_edge_filter(Arena_T arena)
{
    AtlasRibbonObservationSet evidence;
    AtlasTrackGrowResult seed;
    AtlasDevelopableCompleteOptions opts;
    memset(&evidence, 0, sizeof evidence);
    memset(&seed, 0, sizeof seed);
    AtlasDevelopableCompleteOptions_default(&opts);
    evidence.nrows = 2;
    evidence.observation_du = evidence.dv = 2.0;
    evidence.ntarget = 4;
    evidence.target = (AtlasRibbonTarget *)ARENA_CALLOC(
        arena, 4, sizeof(*evidence.target));
    seed.nv = 4;
    seed.nf = 2;
    seed.xyz = (double *)ARENA_CALLOC(
        arena, 12, sizeof(*seed.xyz));
    seed.source_target = (int32_t *)ARENA_ALLOC(
        arena, 4 * sizeof(*seed.source_target));
    seed.vertex_component = (int32_t *)ARENA_CALLOC(
        arena, 4, sizeof(*seed.vertex_component));
    seed.vertex_winding = (int32_t *)ARENA_CALLOC(
        arena, 4, sizeof(*seed.vertex_winding));
    seed.faces = (int32_t *)ARENA_ALLOC(
        arena, 6 * sizeof(*seed.faces));
    static const int32_t faces[6] = {0, 1, 2, 1, 3, 2};
    memcpy(seed.faces, faces, sizeof faces);
    for (int row = 0; row < 2; row++) {
        for (int column = 0; column < 2; column++) {
            int index = row * 2 + column;
            AtlasRibbonTarget *target = &evidence.target[index];
            target->accepted = 1;
            target->row = row;
            target->column = column;
            target->u = 2.0 * column;
            target->v = 2.0 * row;
            target->p[0] = target->u;
            target->tangent[0] = 1.0;
            seed.source_target[index] = index;
            seed.xyz[(size_t)index * 3] = target->v;
            seed.xyz[(size_t)index * 3 + 1] = target->u;
        }
    }
    uint8_t *cell = NULL;
    size_t rows = 0, columns = 0;
    AdcRowGap *gap = NULL;
    size_t gap_count = 0;
    opts.min_gap_u = opts.min_gap_v = 2.0;
    opts.fill_internal_gaps = opts.fill_vertical_gaps = 1;
    int failed =
        adc_seed_existing_cells(
            arena, &evidence, &seed, &cell, &rows, &columns) != 0 ||
        rows != 2 || columns != 2 || cell == NULL || !cell[0] ||
        !adc_existing_u_edge(cell, rows, columns, 0, 0) ||
        !adc_existing_u_edge(cell, rows, columns, 1, 0) ||
        !adc_existing_v_edge(cell, rows, columns, 0, 0) ||
        !adc_existing_v_edge(cell, rows, columns, 0, 1) ||
        adc_build_gaps(
            arena, &evidence, &seed, &opts, &gap, &gap_count) != 0 ||
        gap_count != 0;
    if (failed) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL occupied-edge "
                "grid=%zux%zu gaps=%zu\n",
                rows, columns, gap_count);
        return 1;
    }
    return 0;
}

static int adc_selftest_direct_cell_direction(Arena_T arena, int direction)
{
    AtlasRibbonObservationSet evidence;
    AtlasTrackGrowResult seed;
    AtlasDevelopableCompleteOptions opts;
    memset(&evidence, 0, sizeof evidence);
    memset(&seed, 0, sizeof seed);
    AtlasDevelopableCompleteOptions_default(&opts);
    evidence.nrows = 2;
    evidence.observation_du = evidence.dv = 2.0;
    evidence.ntarget = 4;
    evidence.target = (AtlasRibbonTarget *)ARENA_CALLOC(
        arena, 4, sizeof(*evidence.target));
    seed.nv = 4;
    seed.ncomponents = seed.selected_components = 4;
    seed.xyz = (double *)ARENA_CALLOC(
        arena, 12, sizeof(*seed.xyz));
    seed.local_uv = (double *)ARENA_CALLOC(
        arena, 8, sizeof(*seed.local_uv));
    seed.global_uv = (double *)ARENA_CALLOC(
        arena, 8, sizeof(*seed.global_uv));
    seed.uv = (double *)ARENA_CALLOC(
        arena, 8, sizeof(*seed.uv));
    seed.source_target = (int32_t *)ARENA_ALLOC(
        arena, 4 * sizeof(*seed.source_target));
    seed.vertex_component = (int32_t *)ARENA_ALLOC(
        arena, 4 * sizeof(*seed.vertex_component));
    seed.vertex_winding = (int32_t *)ARENA_CALLOC(
        arena, 4, sizeof(*seed.vertex_winding));
    seed.vertex_rank = (uint64_t *)ARENA_CALLOC(
        arena, 4, sizeof(*seed.vertex_rank));
    for (int row = 0; row < 2; row++) {
        for (int column = 0; column < 2; column++) {
            int index = row * 2 + column;
            AtlasRibbonTarget *target = &evidence.target[index];
            target->accepted = 1;
            target->row = row;
            target->column = column;
            target->u = 2.0 * column;
            target->v = 2.0 * row;
            target->p[0] = target->u;
            target->tangent[0] = 1.0;
            seed.source_target[index] = index;
            seed.vertex_component[index] = index;
            seed.xyz[(size_t)index * 3] = target->v;
            seed.xyz[(size_t)index * 3 + 1] = target->u;
            seed.local_uv[(size_t)index * 2] = target->u;
            seed.local_uv[(size_t)index * 2 + 1] = target->v;
            seed.global_uv[(size_t)index * 2] = target->u;
            seed.global_uv[(size_t)index * 2 + 1] = target->v;
            seed.uv[(size_t)index * 2] = target->u;
            seed.uv[(size_t)index * 2 + 1] = target->v;
        }
    }
    AdcRowGap gap[2];
    memset(gap, 0, sizeof gap);
    for (int support = 0; support < 2; support++) {
        gap[support].direction = direction;
        gap[support].row = support;
        gap[support].column_left = 0;
        gap[support].column_right = 1;
        if (direction == 0) {
            gap[support].vertex_left = support * 2;
            gap[support].vertex_right = support * 2 + 1;
        } else {
            gap[support].vertex_left = support;
            gap[support].vertex_right = support + 2;
        }
    }
    AdcDirectCell *direct = NULL;
    size_t direct_count = 0, u_count = 0, v_count = 0;
    AdcLattice lattice;
    AtlasTrackGrowResult mesh;
    memset(&lattice, 0, sizeof lattice);
    memset(&mesh, 0, sizeof mesh);
    size_t vertices = 0, faces_added = 0, emitted = 0, components = 0;
    size_t short_patches = 0, short_cells = 0;
    size_t short_u = 0, short_v = 0;
    int failed = adc_lattice_init(
        arena, &evidence, &opts, &lattice) != 0;
    if (!failed)
        failed = adc_build_direct_cells(
            arena, &evidence, &seed, &opts, gap, 2, 2, 2, &lattice,
            &direct, &direct_count, &u_count, &v_count,
            &short_patches, &short_cells, &short_u, &short_v) != 0;
    if (!failed)
        failed = direct_count != 1 ||
            u_count != (size_t)(direction == 0) ||
            v_count != (size_t)(direction == 1) ||
            short_patches != 0 || short_cells != 0;
    if (!failed)
        failed = adc_lattice_rebuild(
            arena, &evidence, &seed, &opts, &lattice,
            direct, direct_count, &mesh, &vertices, &faces_added,
            &emitted) != 0;
    if (!failed)
        failed = vertices != 0 || faces_added != 2 || emitted != 1 ||
            mesh.nv != 4 || mesh.nf != 2 ||
            adc_remap_components(
                arena, &seed, &mesh, &components) != 0 ||
            components != 1;
    if (failed) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL adjacent-%c "
                "cells=%zu u/v=%zu/%zu emitted=%zu mesh=%zu/%zu "
                "components=%zu\n",
                direction == 0 ? 'u' : 'v', direct_count,
                u_count, v_count, emitted, mesh.nv, mesh.nf, components);
        return 1;
    }
    return 0;
}

static int adc_selftest_direct_cells(Arena_T arena)
{
    if (adc_selftest_direct_cell_direction(arena, 0) != 0) return 1;
    return adc_selftest_direct_cell_direction(arena, 1);
}

static int adc_selftest_short_strip_direction(Arena_T arena, int direction)
{
    AtlasRibbonObservationSet evidence;
    AtlasTrackGrowResult seed, mesh;
    AtlasDevelopableCompleteOptions opts;
    memset(&evidence, 0, sizeof evidence);
    memset(&seed, 0, sizeof seed);
    memset(&mesh, 0, sizeof mesh);
    AtlasDevelopableCompleteOptions_default(&opts);
    evidence.nrows = direction == 0 ? 2 : 4;
    evidence.observation_du = evidence.dv = 2.0;
    evidence.ntarget = 4;
    evidence.target = (AtlasRibbonTarget *)ARENA_CALLOC(
        arena, 4, sizeof(*evidence.target));
    seed.nv = 4;
    seed.ncomponents = seed.selected_components = 4;
    seed.xyz = (double *)ARENA_CALLOC(
        arena, 12, sizeof(*seed.xyz));
    seed.local_uv = (double *)ARENA_CALLOC(
        arena, 8, sizeof(*seed.local_uv));
    seed.global_uv = (double *)ARENA_CALLOC(
        arena, 8, sizeof(*seed.global_uv));
    seed.uv = (double *)ARENA_CALLOC(
        arena, 8, sizeof(*seed.uv));
    seed.source_target = (int32_t *)ARENA_ALLOC(
        arena, 4 * sizeof(*seed.source_target));
    seed.vertex_component = (int32_t *)ARENA_ALLOC(
        arena, 4 * sizeof(*seed.vertex_component));
    seed.vertex_winding = (int32_t *)ARENA_CALLOC(
        arena, 4, sizeof(*seed.vertex_winding));
    seed.vertex_rank = (uint64_t *)ARENA_CALLOC(
        arena, 4, sizeof(*seed.vertex_rank));
    static const int corner_row[4] = {0, 0, 1, 1};
    static const int corner_column[4] = {0, 1, 0, 1};
    for (int i = 0; i < 4; i++) {
        int row = corner_row[i] * (direction == 0 ? 1 : 3);
        int column = corner_column[i] * (direction == 0 ? 3 : 1);
        AtlasRibbonTarget *target = &evidence.target[i];
        target->accepted = 1;
        target->row = row;
        target->column = column;
        target->u = 2.0 * column;
        target->v = 2.0 * row;
        target->p[0] = target->u;
        target->tangent[0] = 1.0;
        seed.source_target[i] = i;
        seed.vertex_component[i] = i;
        seed.xyz[(size_t)i * 3] = target->v;
        seed.xyz[(size_t)i * 3 + 1] = target->u;
        seed.local_uv[(size_t)i * 2] = target->u;
        seed.local_uv[(size_t)i * 2 + 1] = target->v;
        seed.global_uv[(size_t)i * 2] = target->u;
        seed.global_uv[(size_t)i * 2 + 1] = target->v;
        seed.uv[(size_t)i * 2] = target->u;
        seed.uv[(size_t)i * 2 + 1] = target->v;
    }
    AdcRowGap gap[2];
    memset(gap, 0, sizeof gap);
    for (int support = 0; support < 2; support++) {
        gap[support].direction = direction;
        gap[support].row = support;
        gap[support].column_left = 0;
        gap[support].column_right = 3;
        if (direction == 0) {
            gap[support].vertex_left = support * 2;
            gap[support].vertex_right = support * 2 + 1;
        } else {
            gap[support].vertex_left = support;
            gap[support].vertex_right = support + 2;
        }
        gap[support].component_left =
            seed.vertex_component[(size_t)gap[support].vertex_left];
        gap[support].component_right =
            seed.vertex_component[(size_t)gap[support].vertex_right];
    }
    AdcLattice lattice;
    memset(&lattice, 0, sizeof lattice);
    AdcDirectCell *direct = NULL;
    size_t direct_count = 0, direct_u = 0, direct_v = 0;
    size_t short_patches = 0, short_cells = 0;
    size_t short_u = 0, short_v = 0;
    size_t vertices = 0, faces = 0, emitted = 0, components = 0;
    int failed = adc_lattice_init(
        arena, &evidence, &opts, &lattice) != 0;
    if (!failed)
        failed = adc_build_direct_cells(
            arena, &evidence, &seed, &opts, gap, 2,
            evidence.nrows, direction == 0 ? 4 : 2, &lattice,
            &direct, &direct_count, &direct_u, &direct_v,
            &short_patches, &short_cells, &short_u, &short_v) != 0;
    if (!failed)
        failed = direct_count != 0 || short_patches != 1 ||
            short_cells != 3 ||
            short_u != (size_t)(direction == 0 ? 3 : 0) ||
            short_v != (size_t)(direction == 1 ? 3 : 0);
    if (!failed)
        failed = adc_lattice_rebuild(
            arena, &evidence, &seed, &opts, &lattice,
            direct, direct_count, &mesh, &vertices, &faces,
            &emitted) != 0;
    if (!failed)
        failed = vertices != 4 || faces != 6 || emitted != 0 ||
            mesh.nv != 8 || mesh.nf != 6 ||
            adc_remap_components(
                arena, &seed, &mesh, &components) != 0 ||
            components != 1;
    if (failed) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL short-strip-%c "
                "direct=%zu short=%zu/%zu u/v=%zu/%zu mesh=%zu/%zu "
                "components=%zu\n",
                direction == 0 ? 'u' : 'v', direct_count,
                short_patches, short_cells, short_u, short_v,
                mesh.nv, mesh.nf, components);
        return 1;
    }
    return 0;
}

static int adc_selftest_short_strips(Arena_T arena)
{
    if (adc_selftest_short_strip_direction(arena, 0) != 0) return 1;
    return adc_selftest_short_strip_direction(arena, 1);
}

static int adc_selftest_topology_remap(Arena_T arena)
{
    AtlasTrackGrowResult seed;
    AtlasTrackGrowResult mesh;
    memset(&seed, 0, sizeof seed);
    memset(&mesh, 0, sizeof mesh);
    seed.ncomponents = 2;
    mesh.nv = 6;
    mesh.nf = 2;
    mesh.faces = (int32_t *)ARENA_ALLOC(
        arena, mesh.nf * 3 * sizeof(*mesh.faces));
    mesh.vertex_component = (int32_t *)ARENA_ALLOC(
        arena, mesh.nv * sizeof(*mesh.vertex_component));
    mesh.face_component = (int32_t *)ARENA_ALLOC(
        arena, mesh.nf * sizeof(*mesh.face_component));
    static const int32_t faces[6] = {0, 1, 2, 3, 4, 5};
    static const int32_t source[6] = {0, 0, 0, 0, 1, 1};
    memcpy(mesh.faces, faces, sizeof faces);
    memcpy(mesh.vertex_component, source, sizeof source);
    mesh.face_component[0] = 0;
    mesh.face_component[1] = 0;

    size_t components = 0;
    if (adc_remap_components(
            arena, &seed, &mesh, &components) != 0 ||
        components != 2 || mesh.selected_components != 2 ||
        mesh.face_component[0] == mesh.face_component[1] ||
        mesh.vertex_component[0] != mesh.vertex_component[2] ||
        mesh.vertex_component[3] != mesh.vertex_component[5] ||
        mesh.vertex_component[0] == mesh.vertex_component[3]) {
        fprintf(stderr,
                "[atlas_developable_complete selftest] FAIL topology-remap "
                "components=%zu faces=%d/%d vertices=%d/%d\n",
                components, mesh.face_component[0], mesh.face_component[1],
                mesh.vertex_component[0], mesh.vertex_component[3]);
        return 1;
    }
    return 0;
}

int AtlasDevelopableComplete_selftest(void)
{
    int failures = 0;
    Arena_T arena = Arena_new();
    if (arena == NULL) return 1;
    failures += adc_selftest_adjoint(arena);
    Arena_free(arena);
    failures += adc_selftest_admm(arena);
    Arena_free(arena);
    failures += adc_selftest_bridge(arena);
    Arena_free(arena);
    failures += adc_selftest_lattice_cross(arena);
    Arena_free(arena);
    failures += adc_selftest_aligned_gaps(arena);
    Arena_free(arena);
    failures += adc_selftest_existing_edge_filter(arena);
    Arena_free(arena);
    failures += adc_selftest_direct_cells(arena);
    Arena_free(arena);
    failures += adc_selftest_short_strips(arena);
    Arena_free(arena);
    failures += adc_selftest_topology_remap(arena);
    Arena_free(arena);
    {
        AtlasRibbonObservationSet evidence;
        AtlasTrackGrowResult seed;
        AtlasDevelopableCompleteOptions opts;
        AtlasDevelopableCompleteResult result;
        memset(&evidence, 0, sizeof evidence);
        memset(&seed, 0, sizeof seed);
        AtlasDevelopableCompleteOptions_default(&opts);
        if (AtlasDevelopableComplete_build(
                arena, &evidence, &seed, &opts, &result) != 0 ||
            result.mesh.nv != 0 || result.mesh.nf != 0) {
            fprintf(stderr,
                    "[atlas_developable_complete selftest] FAIL empty\n");
            failures++;
        }
    }
    Arena_dispose(&arena);
    fprintf(stderr,
            "[atlas_developable_complete selftest] %s (%d failures)\n",
            failures == 0 ? "ALL PASS" : "FAILURES", failures);
    return failures;
}
