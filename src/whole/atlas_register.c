#include "atlas_register.h"

#include "../common/union_find.h"
#include "../flatten/sparse_solve.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t island;
    double value;
    double weight;
} ArPrior;

static int ar_prior_compare(const void *left, const void *right)
{
    const ArPrior *a = (const ArPrior *)left;
    const ArPrior *b = (const ArPrior *)right;
    if (a->island < b->island) return -1;
    if (a->island > b->island) return 1;
    if (a->value < b->value) return -1;
    if (a->value > b->value) return 1;
    return 0;
}

void AtlasRegisterOptions_default(AtlasRegisterOptions *options)
{
    if (options == NULL) return;
    memset(options, 0, sizeof(*options));
    options->lambda_weld = 1.0;
    options->lambda_radial = 0.25;
    options->lambda_local = 1.0;
    options->irls_rounds = 4;
    options->eigen_power_iterations = 20;
    options->eigen_inverse_iterations = 6;
    options->max_condition = 1.0e9;
    options->energy_relative_tolerance = 1.0e-9;
    options->convergence_tolerance = 1.0e-8;
}

static double ar_kind_lambda(const AtlasRegisterOptions *options,
                             AtlasRegisterEdgeKind kind)
{
    switch (kind) {
    case ATLAS_REGISTER_WELD: return options->lambda_weld;
    case ATLAS_REGISTER_RADIAL: return options->lambda_radial;
    case ATLAS_REGISTER_LOCAL: return options->lambda_local;
    default: return 0.0;
    }
}

static int ar_edge_usable(const AtlasRegisterProblem *problem,
                          const AtlasRegisterOptions *options, size_t index)
{
    const AtlasRegisterEdge *edge = &problem->edges[index];
    return edge->component0 >= 0 && edge->component1 >= 0 &&
           (size_t)edge->component0 < problem->ncomponents &&
           (size_t)edge->component1 < problem->ncomponents &&
           edge->component0 != edge->component1 && isfinite(edge->target) &&
           isfinite(edge->weight) && edge->weight > 0.0 &&
           isfinite(edge->robust_scale) &&
           (edge->kind == ATLAS_REGISTER_WELD ||
            edge->kind == ATLAS_REGISTER_RADIAL ||
            edge->kind == ATLAS_REGISTER_LOCAL) &&
           isfinite(ar_kind_lambda(options, edge->kind)) &&
           ar_kind_lambda(options, edge->kind) > 0.0;
}

static double ar_residual(const AtlasRegisterEdge *edge,
                          const double *relative_shift)
{
    return relative_shift[edge->component1] -
           relative_shift[edge->component0] - edge->target;
}

static double ar_huber_weight(double residual, double scale)
{
    double magnitude = fabs(residual);
    if (!(scale > 0.0) || magnitude <= scale || magnitude == 0.0) return 1.0;
    return scale / magnitude;
}

static double ar_huber_loss(double residual, double scale)
{
    double magnitude = fabs(residual);
    if (!(scale > 0.0) || magnitude <= scale)
        return 0.5 * residual * residual;
    return scale * (magnitude - 0.5 * scale);
}

static void ar_energies(const AtlasRegisterProblem *problem,
                        const AtlasRegisterOptions *options,
                        const unsigned char *usable,
                        const double *frozen_weight,
                        const double *relative_shift,
                        double *frozen_energy, double *robust_energy)
{
    double quadratic = 0.0;
    double robust = 0.0;
    for (size_t i = 0; i < problem->nedges; i++) {
        if (!usable[i]) continue;
        const AtlasRegisterEdge *edge = &problem->edges[i];
        double residual = ar_residual(edge, relative_shift);
        double base = edge->weight * ar_kind_lambda(options, edge->kind);
        quadratic += 0.5 * frozen_weight[i] * residual * residual;
        robust += base * ar_huber_loss(residual, edge->robust_scale);
    }
    *frozen_energy = quadratic;
    *robust_energy = robust;
}

static void ar_residual_rms(const AtlasRegisterProblem *problem,
                            const AtlasRegisterOptions *options,
                            const unsigned char *usable,
                            const double *relative_shift,
                            AtlasRegisterEdgeKind kind, double *out_rms)
{
    double weighted_square = 0.0;
    double weight_sum = 0.0;
    for (size_t i = 0; i < problem->nedges; i++) {
        if (!usable[i] || problem->edges[i].kind != kind) continue;
        const AtlasRegisterEdge *edge = &problem->edges[i];
        double weight = edge->weight * ar_kind_lambda(options, kind);
        double residual = ar_residual(edge, relative_shift);
        weighted_square += weight * residual * residual;
        weight_sum += weight;
    }
    *out_rms = weight_sum > 0.0 ? sqrt(weighted_square / weight_sum) : 0.0;
}

static void ar_matvec(size_t n, size_t nt, const int *row, const int *column,
                      const double *value, const double *x, double *result)
{
    memset(result, 0, n * sizeof(*result));
    for (size_t k = 0; k < nt; k++) {
        int r = row[k];
        int c = column[k];
        result[r] += value[k] * x[c];
        if (r != c) result[c] += value[k] * x[r];
    }
}

static double ar_dot(size_t n, const double *a, const double *b)
{
    double result = 0.0;
    for (size_t i = 0; i < n; i++) result += a[i] * b[i];
    return result;
}

static int ar_normalize(size_t n, double *x)
{
    double norm2 = ar_dot(n, x, x);
    if (!(norm2 > 0.0) || !isfinite(norm2)) return -1;
    double inverse_norm = 1.0 / sqrt(norm2);
    for (size_t i = 0; i < n; i++) x[i] *= inverse_norm;
    return 0;
}

/* Power and inverse-power Rayleigh estimates on the actual equilibrated
 * matrix.  Inverse iteration deliberately uses the same SPD Cholesky path as
 * the production solve; a factorization failure is therefore caught before
 * accepting an IRLS step. */
static int ar_condition_estimate(Arena_T arena, size_t n, size_t nt,
                                 const int *row, const int *column,
                                 const double *value,
                                 const AtlasRegisterOptions *options,
                                 double *out_minimum, double *out_maximum)
{
    if (n == 0) {
        *out_minimum = 1.0;
        *out_maximum = 1.0;
        return 0;
    }
    if (n > (size_t)INT_MAX || nt > (size_t)INT_MAX) return -1;
    Arena_Mark mark = Arena_save(arena);
    double *vector = (double *)ARENA_ALLOC(arena, n * sizeof(*vector));
    double *product = (double *)ARENA_ALLOC(arena, n * sizeof(*product));
    double *solution = (double *)ARENA_ALLOC(arena, n * sizeof(*solution));

    for (size_t i = 0; i < n; i++)
        vector[i] = 1.0 + (double)((i * 37u) % 19u) / 19.0;
    if (ar_normalize(n, vector) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    int power_iterations = options->eigen_power_iterations;
    if (power_iterations < 2) power_iterations = 2;
    for (int iteration = 0; iteration < power_iterations; iteration++) {
        ar_matvec(n, nt, row, column, value, vector, product);
        if (ar_normalize(n, product) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        memcpy(vector, product, n * sizeof(*vector));
    }
    ar_matvec(n, nt, row, column, value, vector, product);
    double maximum = ar_dot(n, vector, product);

    for (size_t i = 0; i < n; i++)
        vector[i] = 0.5 + (double)((i * 53u + 7u) % 23u) / 23.0;
    if (ar_normalize(n, vector) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    int inverse_iterations = options->eigen_inverse_iterations;
    if (inverse_iterations < 2) inverse_iterations = 2;
    for (int iteration = 0; iteration < inverse_iterations; iteration++) {
        memset(solution, 0, n * sizeof(*solution));
        if (Sparse_solve_sym((int)n, (int)nt, row, column, value,
                             vector, solution, SPARSE_SPD) != 0 ||
            ar_normalize(n, solution) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        memcpy(vector, solution, n * sizeof(*vector));
    }
    ar_matvec(n, nt, row, column, value, vector, product);
    double minimum = ar_dot(n, vector, product);
    Arena_restore(arena, mark);
    if (!(minimum > 0.0) || !(maximum >= minimum) ||
        !isfinite(minimum) || !isfinite(maximum))
        return -1;
    *out_minimum = minimum;
    *out_maximum = maximum;
    return 0;
}

static size_t ar_apply_island_gauge(Arena_T arena,
                                    const AtlasRegisterProblem *problem,
                                    const int32_t *island,
                                    size_t nislands,
                                    const double *relative_shift,
                                    double *absolute_shift)
{
    Arena_Mark mark = Arena_save(arena);
    ArPrior *prior = (ArPrior *)ARENA_ALLOC(
        arena, (problem->ncomponents ? problem->ncomponents : 1) *
                   sizeof(*prior));
    size_t nprior = 0;
    if (problem->prior_shift != NULL) {
        for (size_t component = 0; component < problem->ncomponents;
             component++) {
            double value = problem->prior_shift[component];
            double weight = problem->prior_weight != NULL
                                ? problem->prior_weight[component] : 1.0;
            if (!isfinite(value) || !isfinite(weight) || !(weight > 0.0))
                continue;
            prior[nprior].island = island[component];
            prior[nprior].value = value - relative_shift[component];
            prior[nprior].weight = weight;
            nprior++;
        }
    }
    qsort(prior, nprior, sizeof(*prior), ar_prior_compare);
    double *gauge = (double *)ARENA_CALLOC(
        arena, nislands ? nislands : 1, sizeof(*gauge));
    for (size_t first = 0; first < nprior; ) {
        size_t last = first + 1;
        while (last < nprior && prior[last].island == prior[first].island)
            last++;
        double total = 0.0;
        for (size_t i = first; i < last; i++) total += prior[i].weight;
        double cumulative = 0.0;
        double median = prior[first].value;
        for (size_t i = first; i < last; i++) {
            cumulative += prior[i].weight;
            median = prior[i].value;
            if (2.0 * cumulative >= total) break;
        }
        gauge[prior[first].island] = median;
        first = last;
    }
    for (size_t component = 0; component < problem->ncomponents; component++)
        absolute_shift[component] =
            relative_shift[component] + gauge[island[component]];
    Arena_restore(arena, mark);
    return nprior;
}

int AtlasRegister_solve(Arena_T arena,
                        const AtlasRegisterProblem *problem,
                        const AtlasRegisterOptions *input_options,
                        double *out_shift,
                        AtlasRegisterStats *stats)
{
    if (arena == NULL || problem == NULL || out_shift == NULL ||
        problem->ncomponents == 0 || problem->ncomponents > (size_t)INT32_MAX ||
        (problem->nedges > 0 && problem->edges == NULL))
        return -1;

    AtlasRegisterOptions defaults;
    AtlasRegisterOptions_default(&defaults);
    const AtlasRegisterOptions *options =
        input_options != NULL ? input_options : &defaults;
    if (!(options->lambda_weld > 0.0) ||
        !(options->lambda_radial > 0.0) ||
        !(options->lambda_local > 0.0) || options->irls_rounds < 1 ||
        options->irls_rounds > ATLAS_REGISTER_MAX_IRLS ||
        !(options->max_condition > 1.0) ||
        !(options->energy_relative_tolerance >= 0.0) ||
        !(options->convergence_tolerance >= 0.0))
        return -1;

    AtlasRegisterStats report;
    memset(&report, 0, sizeof(report));
    report.components = problem->ncomponents;
    report.condition_min = DBL_MAX;
    report.raw_diagonal_min = DBL_MAX;

    unsigned char *usable = (unsigned char *)ARENA_CALLOC(
        arena, problem->nedges ? problem->nedges : 1, sizeof(*usable));
    UnionFind components = UF_new(arena, (int32_t)problem->ncomponents);
    double *degree = (double *)ARENA_CALLOC(
        arena, problem->ncomponents, sizeof(*degree));
    for (size_t i = 0; i < problem->nedges; i++) {
        if (!ar_edge_usable(problem, options, i)) {
            report.rejected_edges++;
            continue;
        }
        usable[i] = 1;
        report.usable_edges++;
        if (problem->edges[i].kind == ATLAS_REGISTER_WELD)
            report.weld_edges++;
        else if (problem->edges[i].kind == ATLAS_REGISTER_RADIAL)
            report.radial_edges++;
        else
            report.local_edges++;
        int32_t c0 = problem->edges[i].component0;
        int32_t c1 = problem->edges[i].component1;
        uf_union(&components, c0, c1);
        double weight = problem->edges[i].weight *
                        ar_kind_lambda(options, problem->edges[i].kind);
        degree[c0] += weight;
        degree[c1] += weight;
    }

    int32_t *root_island = (int32_t *)ARENA_ALLOC(
        arena, problem->ncomponents * sizeof(*root_island));
    int32_t *island = (int32_t *)ARENA_ALLOC(
        arena, problem->ncomponents * sizeof(*island));
    int32_t *anchor = (int32_t *)ARENA_ALLOC(
        arena, problem->ncomponents * sizeof(*anchor));
    for (size_t i = 0; i < problem->ncomponents; i++) {
        root_island[i] = -1;
        anchor[i] = -1;
    }
    size_t nislands = 0;
    for (size_t component = 0; component < problem->ncomponents; component++) {
        int32_t root = uf_find(&components, (int32_t)component);
        if (root_island[root] < 0) root_island[root] = (int32_t)nislands++;
        island[component] = root_island[root];
        int32_t *choice = &anchor[island[component]];
        if (*choice < 0 || degree[component] > degree[*choice])
            *choice = (int32_t)component;
    }
    report.evidence_islands = nislands;

    int32_t *variable = (int32_t *)ARENA_ALLOC(
        arena, problem->ncomponents * sizeof(*variable));
    size_t nvariables = 0;
    for (size_t component = 0; component < problem->ncomponents; component++) {
        if (anchor[island[component]] == (int32_t)component)
            variable[component] = -1;
        else
            variable[component] = (int32_t)nvariables++;
    }
    report.variables = nvariables;

    double *relative = (double *)ARENA_CALLOC(
        arena, problem->ncomponents, sizeof(*relative));
    double *candidate = (double *)ARENA_CALLOC(
        arena, problem->ncomponents, sizeof(*candidate));
    double *frozen_weight = (double *)ARENA_CALLOC(
        arena, problem->nedges ? problem->nedges : 1,
        sizeof(*frozen_weight));

    ar_residual_rms(problem, options, usable, relative,
                    ATLAS_REGISTER_WELD, &report.weld_residual_rms_before);
    ar_residual_rms(problem, options, usable, relative,
                    ATLAS_REGISTER_RADIAL, &report.radial_residual_rms_before);
    ar_residual_rms(problem, options, usable, relative,
                    ATLAS_REGISTER_LOCAL, &report.local_residual_rms_before);

    if (nvariables == 0) {
        report.trusted_priors = ar_apply_island_gauge(
            arena, problem, island, nislands, relative, out_shift);
        report.condition_min = 1.0;
        report.condition_max = 1.0;
        report.raw_diagonal_min = 0.0;
        if (stats != NULL) *stats = report;
        return 0;
    }
    if (nvariables > (size_t)INT_MAX ||
        report.usable_edges > (size_t)INT_MAX - nvariables) {
        report.solve_failed = 1;
        if (stats != NULL) *stats = report;
        return -1;
    }

    size_t triplet_capacity = nvariables + report.usable_edges;
    int *row = (int *)ARENA_ALLOC(arena, triplet_capacity * sizeof(*row));
    int *column =
        (int *)ARENA_ALLOC(arena, triplet_capacity * sizeof(*column));
    double *value =
        (double *)ARENA_ALLOC(arena, triplet_capacity * sizeof(*value));
    double *diagonal =
        (double *)ARENA_ALLOC(arena, nvariables * sizeof(*diagonal));
    double *scale = (double *)ARENA_ALLOC(arena, nvariables * sizeof(*scale));
    double *rhs = (double *)ARENA_ALLOC(arena, nvariables * sizeof(*rhs));
    double *balanced_rhs =
        (double *)ARENA_ALLOC(arena, nvariables * sizeof(*balanced_rhs));
    double *solution =
        (double *)ARENA_ALLOC(arena, nvariables * sizeof(*solution));

    for (int round = 0; round < options->irls_rounds; round++) {
        size_t downweighted = 0;
        for (size_t edge_index = 0; edge_index < problem->nedges;
             edge_index++) {
            if (!usable[edge_index]) continue;
            const AtlasRegisterEdge *edge = &problem->edges[edge_index];
            double robust_weight = ar_huber_weight(
                ar_residual(edge, relative), edge->robust_scale);
            if (robust_weight < 1.0 - 1.0e-12) downweighted++;
            frozen_weight[edge_index] = edge->weight *
                ar_kind_lambda(options, edge->kind) * robust_weight;
        }

        memset(diagonal, 0, nvariables * sizeof(*diagonal));
        memset(rhs, 0, nvariables * sizeof(*rhs));
        size_t nt = nvariables;
        for (size_t i = 0; i < problem->nedges; i++) {
            if (!usable[i]) continue;
            const AtlasRegisterEdge *edge = &problem->edges[i];
            int32_t v0 = variable[edge->component0];
            int32_t v1 = variable[edge->component1];
            double weight = frozen_weight[i];
            if (!(weight > 0.0) || !isfinite(weight)) {
                report.solve_failed = 1;
                if (stats != NULL) *stats = report;
                return -1;
            }
            if (v0 >= 0) {
                diagonal[v0] += weight;
                rhs[v0] -= weight * edge->target;
            }
            if (v1 >= 0) {
                diagonal[v1] += weight;
                rhs[v1] += weight * edge->target;
            }
            if (v0 >= 0 && v1 >= 0) {
                int r = v0 > v1 ? v0 : v1;
                int c = v0 > v1 ? v1 : v0;
                row[nt] = r;
                column[nt] = c;
                value[nt] = -weight;
                nt++;
            }
        }

        double diagonal_min = DBL_MAX;
        double diagonal_max = 0.0;
        for (size_t i = 0; i < nvariables; i++) {
            if (!(diagonal[i] > 0.0) || !isfinite(diagonal[i])) {
                report.solve_failed = 1;
                if (stats != NULL) *stats = report;
                return -1;
            }
            if (diagonal[i] < diagonal_min) diagonal_min = diagonal[i];
            if (diagonal[i] > diagonal_max) diagonal_max = diagonal[i];
            scale[i] = 1.0 / sqrt(diagonal[i]);
            row[i] = (int)i;
            column[i] = (int)i;
            value[i] = 1.0;
            balanced_rhs[i] = scale[i] * rhs[i];
        }
        for (size_t i = nvariables; i < nt; i++)
            value[i] *= scale[row[i]] * scale[column[i]];

        AtlasRegisterIterationStats iteration;
        memset(&iteration, 0, sizeof(iteration));
        iteration.iteration = round + 1;
        iteration.downweighted_edges = downweighted;
        iteration.raw_diagonal_min = diagonal_min;
        iteration.raw_diagonal_max = diagonal_max;
        ar_energies(problem, options, usable, frozen_weight, relative,
                    &iteration.frozen_energy_before,
                    &iteration.robust_energy_before);
        if (round == 0)
            report.robust_energy_initial = iteration.robust_energy_before;

        if (ar_condition_estimate(arena, nvariables, nt, row, column, value,
                                  options, &iteration.lambda_min,
                                  &iteration.lambda_max) != 0) {
            report.solve_failed = 1;
            if (stats != NULL) *stats = report;
            return -1;
        }
        iteration.condition_estimate =
            iteration.lambda_max / iteration.lambda_min;
        if (!isfinite(iteration.condition_estimate) ||
            iteration.condition_estimate > options->max_condition) {
            report.condition_failed = 1;
            if (stats != NULL) *stats = report;
            return -1;
        }

        memset(solution, 0, nvariables * sizeof(*solution));
        if (Sparse_solve_sym((int)nvariables, (int)nt, row, column, value,
                             balanced_rhs, solution, SPARSE_SPD) != 0) {
            report.solve_failed = 1;
            if (stats != NULL) *stats = report;
            return -1;
        }
        memcpy(candidate, relative,
               problem->ncomponents * sizeof(*candidate));
        for (size_t component = 0; component < problem->ncomponents;
             component++) {
            int32_t v = variable[component];
            if (v >= 0) candidate[component] = scale[v] * solution[v];
            else candidate[component] = 0.0;
            if (!isfinite(candidate[component])) {
                report.solve_failed = 1;
                if (stats != NULL) *stats = report;
                return -1;
            }
        }
        for (size_t component = 0; component < problem->ncomponents;
             component++) {
            double step = fabs(candidate[component] - relative[component]);
            if (step > iteration.max_relative_step)
                iteration.max_relative_step = step;
        }
        ar_energies(problem, options, usable, frozen_weight, candidate,
                    &iteration.frozen_energy_after,
                    &iteration.robust_energy_after);
        double frozen_tolerance = options->energy_relative_tolerance *
            fmax(1.0, fabs(iteration.frozen_energy_before));
        double robust_tolerance = options->energy_relative_tolerance *
            fmax(1.0, fabs(iteration.robust_energy_before));
        if (iteration.frozen_energy_after >
                iteration.frozen_energy_before + frozen_tolerance ||
            iteration.robust_energy_after >
                iteration.robust_energy_before + robust_tolerance) {
            report.energy_failed = 1;
            if (stats != NULL) *stats = report;
            return -1;
        }

        memcpy(relative, candidate,
               problem->ncomponents * sizeof(*relative));
        ar_residual_rms(problem, options, usable, relative,
                        ATLAS_REGISTER_WELD,
                        &iteration.weld_residual_rms);
        ar_residual_rms(problem, options, usable, relative,
                        ATLAS_REGISTER_RADIAL,
                        &iteration.radial_residual_rms);
        ar_residual_rms(problem, options, usable, relative,
                        ATLAS_REGISTER_LOCAL,
                        &iteration.local_residual_rms);
        report.trusted_priors = ar_apply_island_gauge(
            arena, problem, island, nislands, relative, out_shift);

        report.iteration[round] = iteration;
        report.irls_rounds_run = round + 1;
        report.downweighted_edges = downweighted;
        report.robust_energy_final = iteration.robust_energy_after;
        if (iteration.condition_estimate < report.condition_min)
            report.condition_min = iteration.condition_estimate;
        if (iteration.condition_estimate > report.condition_max)
            report.condition_max = iteration.condition_estimate;
        if (diagonal_min < report.raw_diagonal_min)
            report.raw_diagonal_min = diagonal_min;
        if (diagonal_max > report.raw_diagonal_max)
            report.raw_diagonal_max = diagonal_max;

        if (options->iteration_fn != NULL &&
            options->iteration_fn(options->iteration_context, problem,
                                  out_shift, &report.iteration[round]) != 0) {
            report.trace_failed = 1;
            if (stats != NULL) *stats = report;
            return -1;
        }
        if (iteration.max_relative_step <= options->convergence_tolerance)
            break;
    }

    ar_residual_rms(problem, options, usable, relative,
                    ATLAS_REGISTER_WELD, &report.weld_residual_rms_after);
    ar_residual_rms(problem, options, usable, relative,
                    ATLAS_REGISTER_RADIAL, &report.radial_residual_rms_after);
    ar_residual_rms(problem, options, usable, relative,
                    ATLAS_REGISTER_LOCAL, &report.local_residual_rms_after);
    double shift_square = 0.0;
    for (size_t component = 0; component < problem->ncomponents; component++) {
        double magnitude = fabs(out_shift[component]);
        shift_square += out_shift[component] * out_shift[component];
        if (magnitude > report.shift_max) report.shift_max = magnitude;
    }
    report.shift_rms = sqrt(shift_square / (double)problem->ncomponents);
    if (report.condition_min == DBL_MAX) report.condition_min = 1.0;
    if (report.raw_diagonal_min == DBL_MAX) report.raw_diagonal_min = 0.0;
    if (stats != NULL) *stats = report;
    return 0;
}

typedef struct {
    int calls;
    int failed;
} ArSelftestTrace;

static int ar_selftest_trace(void *context,
                             const AtlasRegisterProblem *problem,
                             const double *shift,
                             const AtlasRegisterIterationStats *iteration)
{
    ArSelftestTrace *trace = (ArSelftestTrace *)context;
    trace->calls++;
    if (problem == NULL || shift == NULL || iteration == NULL ||
        iteration->iteration != trace->calls ||
        !isfinite(iteration->condition_estimate) ||
        iteration->frozen_energy_after > iteration->frozen_energy_before +
                                             1.0e-7 ||
        iteration->robust_energy_after > iteration->robust_energy_before +
                                             1.0e-7)
        trace->failed = 1;
    return trace->failed ? -1 : 0;
}

int AtlasRegister_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;
    AtlasRegisterEdge edge[] = {
        {0, 1, 10.0, 50.0, 2.0, ATLAS_REGISTER_WELD},
        {1, 2, 10.0, 50.0, 2.0, ATLAS_REGISTER_WELD},
        {2, 3, 10.0, 50.0, 2.0, ATLAS_REGISTER_WELD},
        {0, 2, 20.0, 10.0, 5.0, ATLAS_REGISTER_RADIAL},
        {1, 3, 20.0, 10.0, 5.0, ATLAS_REGISTER_RADIAL},
        /* Contradictory radial evidence must be downweighted, not unioned. */
        {0, 3, 80.0, 1.0, 5.0, ATLAS_REGISTER_RADIAL}
    };
    double prior[] = {100.0, 110.0, 120.0, 130.0, 7.0};
    double prior_weight[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    AtlasRegisterProblem problem;
    problem.ncomponents = 5;
    problem.edges = edge;
    problem.nedges = sizeof(edge) / sizeof(edge[0]);
    problem.prior_shift = prior;
    problem.prior_weight = prior_weight;

    AtlasRegisterOptions options;
    AtlasRegisterOptions_default(&options);
    ArSelftestTrace trace;
    memset(&trace, 0, sizeof(trace));
    options.iteration_fn = ar_selftest_trace;
    options.iteration_context = &trace;
    double shift[5];
    AtlasRegisterStats stats;
    int rc = AtlasRegister_solve(arena, &problem, &options, shift, &stats);
    int failed = rc != 0 || trace.failed || trace.calls != stats.irls_rounds_run ||
                 stats.evidence_islands != 2 || stats.variables != 3 ||
                 stats.downweighted_edges == 0 || stats.condition_failed ||
                 stats.energy_failed || stats.solve_failed ||
                 !(stats.condition_max < options.max_condition) ||
                 fabs(shift[4] - 7.0) > 1.0e-9;
    for (int i = 0; i < 4; i++) {
        double expected = 100.0 + 10.0 * (double)i;
        if (fabs(shift[i] - expected) > 0.25) failed = 1;
    }
    if (!failed)
        fprintf(stderr,
                "[atlas_register selftest] PASS: islands=%zu vars=%zu "
                "rounds=%d cond_max=%.6g energy=%.6g->%.6g\n",
                stats.evidence_islands, stats.variables,
                stats.irls_rounds_run, stats.condition_max,
                stats.robust_energy_initial, stats.robust_energy_final);
    else
        fprintf(stderr,
                "[atlas_register selftest] FAIL: rc=%d trace=%d/%d "
                "islands=%zu vars=%zu downweighted=%zu cond=%.6g "
                "shift=[%.6g %.6g %.6g %.6g %.6g]\n",
                rc, trace.calls, trace.failed, stats.evidence_islands,
                stats.variables, stats.downweighted_edges,
                stats.condition_max, shift[0], shift[1], shift[2], shift[3],
                shift[4]);
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}
