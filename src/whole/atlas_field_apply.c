#include "atlas_field_apply.h"

#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t root;
    double value;
} AfaGaugeResidual;

typedef struct {
    size_t nvertices;
    size_t nsamples;
    size_t nsupport;
    const float *raw_u;
    const double *sample_target;
    const int32_t *sample_component;
    const int32_t *vertex_mesh;
    double *corrected_target;
    double *field_u;
    AtlasFieldApplyIterationFn iteration_fn;
    void *iteration_context;
} AfaIterationContext;

static double afa_dist3(const float *a, const float *b)
{
    double x = (double)a[0] - (double)b[0];
    double y = (double)a[1] - (double)b[1];
    double z = (double)a[2] - (double)b[2];
    return sqrt(x * x + y * y + z * z);
}

static int afa_compare_gauge_residual(const void *pa, const void *pb)
{
    const AfaGaugeResidual *a = (const AfaGaugeResidual *)pa;
    const AfaGaugeResidual *b = (const AfaGaugeResidual *)pb;
    if (a->root != b->root) return a->root < b->root ? -1 : 1;
    return a->value < b->value ? -1 : (a->value > b->value ? 1 : 0);
}

static void afa_apply_shifts(const AfaIterationContext *context,
                             const double *shift)
{
    for (size_t i = 0; i < context->nsamples; i++) {
        int32_t support = context->sample_component[i];
        context->corrected_target[i] =
            context->sample_target[i] + shift[support];
    }
    for (size_t i = 0; i < context->nvertices; i++) {
        int32_t mesh = context->vertex_mesh[i];
        context->field_u[i] = (double)context->raw_u[i] +
            shift[context->nsupport + (size_t)mesh];
    }
}

static int afa_registration_iteration(
    void *opaque, const AtlasRegisterProblem *problem,
    const double *absolute_shift,
    const AtlasRegisterIterationStats *iteration)
{
    AfaIterationContext *context = (AfaIterationContext *)opaque;
    if (context == NULL || problem == NULL || absolute_shift == NULL ||
        iteration == NULL ||
        problem->ncomponents < context->nsupport)
        return -1;
    afa_apply_shifts(context, absolute_shift);
    return context->iteration_fn == NULL
        ? 0
        : context->iteration_fn(context->iteration_context,
                                context->corrected_target, context->field_u,
                                iteration);
}

int AtlasFieldApply_solve(Arena_T arena,
                          const float *vertices,
                          size_t nvertices,
                          const int32_t *triangles,
                          size_t ntriangles,
                          const float *raw_u,
                          const float *reference_u,
                          const AtlasCandidateSampleRef *sample_ref,
                          const double *sample_target,
                          const int32_t *sample_component,
                          size_t nsamples,
                          size_t nsupport_components,
                          AtlasFieldApplyIterationFn iteration_fn,
                          void *iteration_context,
                          double *out_corrected_target,
                          double *out_u,
                          AtlasFieldApplyStats *stats)
{
    if (arena == NULL || vertices == NULL || nvertices == 0 ||
        nvertices > (size_t)INT32_MAX || triangles == NULL ||
        ntriangles == 0 || raw_u == NULL || sample_ref == NULL ||
        sample_target == NULL || nsamples == 0 || sample_component == NULL ||
        nsupport_components == 0 ||
        nsupport_components > (size_t)INT32_MAX ||
        out_corrected_target == NULL || out_u == NULL || stats == NULL)
        return -1;
    memset(stats, 0, sizeof(*stats));

    for (size_t i = 0; i < nvertices; i++) {
        if (!isfinite(raw_u[i])) return -1;
        for (int d = 0; d < 3; d++)
            if (!isfinite(vertices[3 * i + (size_t)d])) return -1;
    }

    UnionFind mesh_graph = UF_new(arena, (int32_t)nvertices);
    for (size_t f = 0; f < ntriangles; f++) {
        int32_t a = triangles[3 * f];
        int32_t b = triangles[3 * f + 1];
        int32_t c = triangles[3 * f + 2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nvertices || (size_t)b >= nvertices ||
            (size_t)c >= nvertices)
            return -1;
        uf_union(&mesh_graph, a, b);
        uf_union(&mesh_graph, a, c);
    }

    int32_t *mesh_index = (int32_t *)ARENA_ALLOC(
        arena, nvertices * sizeof(*mesh_index));
    int32_t *vertex_mesh = (int32_t *)ARENA_ALLOC(
        arena, nvertices * sizeof(*vertex_mesh));
    for (size_t i = 0; i < nvertices; i++) mesh_index[i] = -1;
    size_t nmesh = 0;
    for (size_t i = 0; i < nvertices; i++) {
        int32_t root = uf_find(&mesh_graph, (int32_t)i);
        if (root == (int32_t)i) {
            if (nmesh >= (size_t)INT32_MAX) return -1;
            mesh_index[root] = (int32_t)nmesh++;
        }
    }
    for (size_t i = 0; i < nvertices; i++) {
        int32_t root = uf_find(&mesh_graph, (int32_t)i);
        if (mesh_index[root] < 0) return -1;
        vertex_mesh[i] = mesh_index[root];
    }
    if (nsupport_components > (size_t)INT32_MAX - nmesh) return -1;
    size_t nnode = nsupport_components + nmesh;

    AtlasRegisterEdge *edge = (AtlasRegisterEdge *)ARENA_ALLOC(
        arena, nsamples * sizeof(*edge));
    size_t *support_count = (size_t *)ARENA_CALLOC(
        arena, nsupport_components, sizeof(*support_count));
    unsigned char *observed_mesh = (unsigned char *)ARENA_CALLOC(
        arena, nmesh, sizeof(*observed_mesh));
    for (size_t i = 0; i < nsamples; i++) {
        int32_t a = sample_ref[i].mesh_vertex[0];
        int32_t b = sample_ref[i].mesh_vertex[1];
        int32_t support = sample_component[i];
        double t = sample_ref[i].mesh_t;
        if (a < 0 || b < 0 || (size_t)a >= nvertices ||
            (size_t)b >= nvertices || support < 0 ||
            (size_t)support >= nsupport_components || !isfinite(t) ||
            t < 0.0 || t > 1.0 || !isfinite(sample_target[i]) ||
            vertex_mesh[a] != vertex_mesh[b])
            return -1;
        int32_t mesh = vertex_mesh[a];
        double raw_value =
            (1.0 - t) * (double)raw_u[a] + t * (double)raw_u[b];
        double local_length = afa_dist3(&vertices[(size_t)a * 3],
                                        &vertices[(size_t)b * 3]);
        edge[i].component0 =
            (int32_t)(nsupport_components + (size_t)mesh);
        edge[i].component1 = support;
        edge[i].target = raw_value - sample_target[i];
        edge[i].weight = 1.0;
        /* Residuals larger than two local mesh-edge lengths are gauge/model
         * outliers, not evidence that a triangle should be stretched. */
        edge[i].robust_scale = fmax(1.0, 2.0 * local_length);
        edge[i].kind = ATLAS_REGISTER_LOCAL;
        support_count[support]++;
        observed_mesh[mesh] = 1;
    }

    double *prior_shift = (double *)ARENA_ALLOC(
        arena, nnode * sizeof(*prior_shift));
    double *prior_weight = (double *)ARENA_CALLOC(
        arena, nnode, sizeof(*prior_weight));
    for (size_t i = 0; i < nnode; i++) prior_shift[i] = (double)NAN;
    if (reference_u != NULL) {
        AfaGaugeResidual *residual = (AfaGaugeResidual *)ARENA_ALLOC(
            arena, nsamples * sizeof(*residual));
        for (size_t i = 0; i < nsamples; i++) {
            int32_t a = sample_ref[i].mesh_vertex[0];
            int32_t b = sample_ref[i].mesh_vertex[1];
            double t = sample_ref[i].mesh_t;
            double reference_value =
                (1.0 - t) * (double)reference_u[a] +
                t * (double)reference_u[b];
            if (!isfinite(reference_value)) return -1;
            residual[i].root = sample_component[i];
            residual[i].value = reference_value - sample_target[i];
        }
        qsort(residual, nsamples, sizeof(*residual),
              afa_compare_gauge_residual);
        for (size_t first = 0; first < nsamples; ) {
            size_t last = first + 1;
            while (last < nsamples &&
                   residual[last].root == residual[first].root)
                last++;
            size_t count = last - first;
            double median = residual[first + count / 2].value;
            if ((count & 1u) == 0)
                median = 0.5 *
                    (median + residual[first + count / 2 - 1].value);
            int32_t support = residual[first].root;
            prior_shift[support] = median;
            prior_weight[support] = (double)count;
            first = last;
        }
        for (size_t vertex = 0; vertex < nvertices; vertex++) {
            int32_t root = uf_find(&mesh_graph, (int32_t)vertex);
            if (root != (int32_t)vertex) continue;
            int32_t mesh = mesh_index[root];
            size_t node = nsupport_components + (size_t)mesh;
            if (!isfinite(reference_u[vertex])) return -1;
            prior_shift[node] =
                (double)reference_u[vertex] - (double)raw_u[vertex];
            /* This is decisive for an unobserved singleton and deliberately
             * negligible beside sample-count-weighted support priors. */
            prior_weight[node] = 1.0;
        }
    }

    AtlasRegisterProblem registration;
    registration.ncomponents = nnode;
    registration.edges = edge;
    registration.nedges = nsamples;
    registration.prior_shift = prior_shift;
    registration.prior_weight = prior_weight;

    AfaIterationContext trace;
    trace.nvertices = nvertices;
    trace.nsamples = nsamples;
    trace.nsupport = nsupport_components;
    trace.raw_u = raw_u;
    trace.sample_target = sample_target;
    trace.sample_component = sample_component;
    trace.vertex_mesh = vertex_mesh;
    trace.corrected_target = out_corrected_target;
    trace.field_u = out_u;
    trace.iteration_fn = iteration_fn;
    trace.iteration_context = iteration_context;

    AtlasRegisterOptions options;
    AtlasRegisterOptions_default(&options);
    options.lambda_local = 1.0;
    options.iteration_fn =
        iteration_fn != NULL ? afa_registration_iteration : NULL;
    options.iteration_context = &trace;
    double *shift = (double *)ARENA_ALLOC(arena, nnode * sizeof(*shift));
    if (AtlasRegister_solve(arena, &registration, &options, shift,
                            &stats->gauge_solver) != 0)
        return -1;
    afa_apply_shifts(&trace, shift);

    stats->mesh_components = nmesh;
    stats->observations = nsamples;
    stats->gauge_graph_components =
        stats->gauge_solver.evidence_islands;
    for (size_t mesh = 0; mesh < nmesh; mesh++) {
        if (observed_mesh[mesh])
            stats->observed_components++;
        else
            stats->anchored_unobserved_components++;
    }

    double correction2 = 0.0;
    double consistency2 = 0.0;
    double reference_shift2 = 0.0;
    double reference_residual2 = 0.0;
    for (size_t i = 0; i < nsamples; i++) {
        int32_t a = sample_ref[i].mesh_vertex[0];
        int32_t b = sample_ref[i].mesh_vertex[1];
        int32_t support = sample_component[i];
        int32_t mesh = vertex_mesh[a];
        double t = sample_ref[i].mesh_t;
        double correction = shift[support];
        double mesh_shift = shift[nsupport_components + (size_t)mesh];
        double field_value = (1.0 - t) * out_u[a] + t * out_u[b];
        double residual = out_corrected_target[i] - field_value;
        correction2 += correction * correction;
        consistency2 += residual * residual;
        reference_shift2 += mesh_shift * mesh_shift;
        if (fabs(correction) > stats->gauge_correction_max)
            stats->gauge_correction_max = fabs(correction);
        if (fabs(mesh_shift) > stats->reference_gauge_shift_max)
            stats->reference_gauge_shift_max = fabs(mesh_shift);
        if (fabs(residual) > stats->observation_max)
            stats->observation_max = fabs(residual);
        if (reference_u != NULL) {
            double reference_value =
                (1.0 - t) * (double)reference_u[a] +
                t * (double)reference_u[b];
            double reference_residual =
                out_corrected_target[i] - reference_value;
            reference_residual2 +=
                reference_residual * reference_residual;
        }
    }
    stats->gauge_correction_rms =
        sqrt(correction2 / (double)nsamples);
    stats->gauge_consistency_rms =
        sqrt(consistency2 / (double)nsamples);
    stats->observation_rms = stats->gauge_consistency_rms;
    stats->reference_gauge_shift_rms =
        sqrt(reference_shift2 / (double)nsamples);
    stats->reference_residual_rms = reference_u != NULL
        ? sqrt(reference_residual2 / (double)nsamples) : 0.0;

    double gradient_error2 = 0.0;
    size_t ngradient = 0;
    for (size_t f = 0; f < ntriangles; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = triangles[3 * f + (size_t)k];
            int32_t b = triangles[3 * f + (size_t)((k + 1) % 3)];
            double length = afa_dist3(&vertices[(size_t)a * 3],
                                      &vertices[(size_t)b * 3]);
            if (length <= 1e-12) continue;
            double before = (double)raw_u[b] - (double)raw_u[a];
            double after = out_u[b] - out_u[a];
            double error = (after - before) / length;
            gradient_error2 += error * error;
            ngradient++;
        }
    }
    stats->edge_gradient_delta_rms = ngradient
        ? sqrt(gradient_error2 / (double)ngradient) : 0.0;
    if (!isfinite(stats->edge_gradient_delta_rms) ||
        stats->edge_gradient_delta_rms > 1.0e-9)
        return -1;
    return 0;
}

typedef struct {
    const float *raw_u;
    int calls;
    int failed;
} AfaSelftestTrace;

static int afa_selftest_iteration(
    void *opaque, const double *corrected_sample, const double *field_u,
    const AtlasRegisterIterationStats *iteration)
{
    AfaSelftestTrace *trace = (AfaSelftestTrace *)opaque;
    trace->calls++;
    if (corrected_sample == NULL || field_u == NULL || iteration == NULL ||
        iteration->iteration != trace->calls ||
        iteration->frozen_energy_after >
            iteration->frozen_energy_before + 1.0e-8 ||
        iteration->robust_energy_after >
            iteration->robust_energy_before + 1.0e-8 ||
        fabs((field_u[1] - field_u[0]) -
             ((double)trace->raw_u[1] - (double)trace->raw_u[0])) > 1.0e-10)
        trace->failed = 1;
    return trace->failed ? -1 : 0;
}

int AtlasFieldApply_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;
    const float vertex[18] = {
        0, 0, 0,  1, 0, 0,  0, 1, 0,
        3, 0, 0,  4, 0, 0,  3, 1, 0
    };
    const int32_t triangle[6] = {0, 1, 2, 3, 4, 5};
    const float raw[6] = {0, 1, 0.25f, 10, 11, 10.25f};
    const float reference[6] = {
        100, 101, 100.25f, -40, -39, -39.75f
    };
    AtlasCandidateSampleRef sample[4];
    memset(sample, 0, sizeof(sample));
    sample[0].mesh_vertex[0] = 0;
    sample[0].mesh_vertex[1] = 1;
    sample[0].mesh_t = 0.25;
    sample[1].mesh_vertex[0] = 0;
    sample[1].mesh_vertex[1] = 2;
    sample[1].mesh_t = 0.5;
    sample[2].mesh_vertex[0] = 3;
    sample[2].mesh_vertex[1] = 4;
    sample[2].mesh_t = 0.75;
    sample[3].mesh_vertex[0] = 3;
    sample[3].mesh_vertex[1] = 5;
    sample[3].mesh_t = 0.5;
    const int32_t support[4] = {0, 1, 0, 1};
    /*
     * Deliberately make the support/chart K2,2 cycle inconsistent by 500
     * u-units.  The old per-vertex field extension could pay this
     * contradiction by shearing a few triangles.  Registration may leave
     * or downweight a residual, but it must preserve each chart's raw edge
     * differences exactly.
     */
    const double target[4] = {100.25, 100.125, -39.25, 460.125};
    double corrected[4], output[6];
    AtlasFieldApplyStats stats;
    AfaSelftestTrace trace = {raw, 0, 0};
    int rc = AtlasFieldApply_solve(
        arena, vertex, 6, triangle, 2, raw, reference, sample, target,
        support, 4, 2, afa_selftest_iteration, &trace, corrected, output,
        &stats);
    int failed = rc != 0 || trace.calls == 0 || trace.failed ||
        stats.gauge_solver.downweighted_edges == 0 ||
        stats.edge_gradient_delta_rms > 1.0e-10 ||
        stats.gauge_solver.condition_max > 1.0e9 ||
        fabs((output[1] - raw[1]) - (output[0] - raw[0])) > 1.0e-9 ||
        fabs((output[2] - raw[2]) - (output[0] - raw[0])) > 1.0e-9 ||
        fabs((output[4] - raw[4]) - (output[3] - raw[3])) > 1.0e-9 ||
        fabs((output[5] - raw[5]) - (output[3] - raw[3])) > 1.0e-9;
    fprintf(stderr,
        "[atlas_field_apply selftest] %s: rounds=%d condition=%.6g "
        "gradient_delta=%.3g\n",
        failed ? "FAILED" : "PASSED", stats.gauge_solver.irls_rounds_run,
        stats.gauge_solver.condition_max, stats.edge_gradient_delta_rms);
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}
