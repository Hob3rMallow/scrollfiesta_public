#include "atlas_xyz_weld_audit.h"

#include "atlas_overlap_audit.h"
#include "../remesh/seam_planes.h"
#include "../remesh/seam_refine.h"
#include "../remesh/seam_weld.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t vertex;
    int32_t chart;
    double u;
} AxwaPosition;

typedef struct {
    int32_t lo;
    int32_t hi;
} AxwaEdge;

typedef struct {
    int32_t vertex;
    double u;
} AxwaEmbeddingPosition;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    double target;
    double length;
} AxwaShiftObservation;

static int axwa_compare_position(const void *pa, const void *pb)
{
    const AxwaPosition *a = (const AxwaPosition *)pa;
    const AxwaPosition *b = (const AxwaPosition *)pb;
    if (a->vertex != b->vertex) return a->vertex < b->vertex ? -1 : 1;
    return a->chart < b->chart ? -1 : (a->chart > b->chart ? 1 : 0);
}

static int axwa_compare_edge(const void *pa, const void *pb)
{
    const AxwaEdge *a = (const AxwaEdge *)pa;
    const AxwaEdge *b = (const AxwaEdge *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    return a->hi < b->hi ? -1 : (a->hi > b->hi ? 1 : 0);
}

static int axwa_compare_embedding(const void *pa, const void *pb)
{
    const AxwaEmbeddingPosition *a = (const AxwaEmbeddingPosition *)pa;
    const AxwaEmbeddingPosition *b = (const AxwaEmbeddingPosition *)pb;
    if (a->vertex != b->vertex) return a->vertex < b->vertex ? -1 : 1;
    return a->u < b->u ? -1 : (a->u > b->u ? 1 : 0);
}

static int axwa_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int axwa_compare_shift_observation(const void *pa, const void *pb)
{
    const AxwaShiftObservation *a = (const AxwaShiftObservation *)pa;
    const AxwaShiftObservation *b = (const AxwaShiftObservation *)pb;
    if (a->chart0 != b->chart0)
        return a->chart0 < b->chart0 ? -1 : 1;
    if (a->chart1 != b->chart1)
        return a->chart1 < b->chart1 ? -1 : 1;
    return a->target < b->target ? -1 : (a->target > b->target ? 1 : 0);
}

static double axwa_quantile(double *value, size_t count, double q)
{
    if (count == 0) return 0.0;
    qsort(value, count, sizeof(double), axwa_compare_double);
    double x = q * (double)(count - 1);
    size_t lo = (size_t)floor(x), hi = (size_t)ceil(x);
    double t = x - (double)lo;
    return value[lo] * (1.0 - t) + value[hi] * t;
}

static size_t axwa_vertex_cube(const PieceSet *ps, size_t vertex)
{
    size_t lo = 0, hi = ps->n_cubes;
    while (lo + 1 < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (ps->cube_voff[mid] <= vertex) lo = mid;
        else hi = mid;
    }
    return lo;
}

static double axwa_xyz_edge(const float *vertices, int32_t a, int32_t b)
{
    double d0 = (double)vertices[(size_t)b * 3] -
                (double)vertices[(size_t)a * 3];
    double d1 = (double)vertices[(size_t)b * 3 + 1] -
                (double)vertices[(size_t)a * 3 + 1];
    double d2 = (double)vertices[(size_t)b * 3 + 2] -
                (double)vertices[(size_t)a * 3 + 2];
    return sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

static double axwa_uv_edge(double ua, double va, double ub, double vb)
{
    double du = ub - ua, dv = vb - va;
    return sqrt(du * du + dv * dv);
}

static double axwa_symmetric_stretch(double uv, double xyz)
{
    if (!(uv > 1.0e-15) || !(xyz > 1.0e-15)) return DBL_MAX;
    double r = uv / xyz;
    return r >= 1.0 ? r : 1.0 / r;
}

int AtlasXyzWeldTopology_build(
    Arena_T arena, const PieceSet *ps, float cube_size, float rho,
    float rho_max, float band, const BpaBridgeGate *gate,
    AtlasXyzWeldTopology *out)
{
    if (arena == NULL || ps == NULL || ps->verts == NULL ||
        ps->faces == NULL || ps->nv == 0 || ps->nf == 0 ||
        ps->n_cubes < 2 || ps->cube_voff == NULL || out == NULL ||
        !(cube_size > 0.0f) || !(rho > 0.0f) || !(band > 0.0f))
        return -1;
    memset(out, 0, sizeof(*out));
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, ps->nv, sizeof(uint8_t));
    for (size_t f = 0; f < ps->nf; f++)
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            used[vertex] = 1;
        }
    SeamPlane planes[64];
    size_t nplanes = SeamPlanes_detect(
        ps->verts, ps->nv, used, (double)cube_size, (double)band,
        planes, sizeof(planes) / sizeof(planes[0]));
    SeamRefineParams refine_params;
    SeamRefine_default_params(&refine_params);
    refine_params.band = band;
    float *refined_vertices = NULL;
    int32_t *refined_faces = NULL;
    int32_t *parent0 = NULL, *parent1 = NULL;
    size_t refined_nv = 0, refined_nf = 0, nnew = 0;
    SeamRefineStats refine_stats;
    if (SeamRefine_process_with_parents(
            arena, ps->verts, ps->nv, ps->faces, ps->nf,
            planes, nplanes, &refine_params,
            &refined_vertices, &refined_nv, &refined_faces, &refined_nf,
            &parent0, &parent1, &nnew, &refine_stats) != 0 ||
        refined_vertices == NULL || refined_faces == NULL ||
        refined_nv < ps->nv || refined_nv - ps->nv != nnew ||
        refined_nf < ps->nf || (nnew > 0 &&
        (parent0 == NULL || parent1 == NULL)))
        return -1;
    int32_t *combined = NULL;
    size_t combined_faces = 0, bridge_faces = 0;
    if (SeamWeld_bridge(
            arena, refined_vertices, refined_nv, refined_faces, refined_nf,
            cube_size, rho, rho_max, band, NULL, gate,
            &combined, &combined_faces, &bridge_faces) != 0 ||
        combined == NULL || bridge_faces > combined_faces)
        return -1;
    out->vertices = refined_vertices;
    out->nvertices = refined_nv;
    out->source_vertices = ps->nv;
    out->parent0 = parent0;
    out->parent1 = parent1;
    out->faces = combined + (combined_faces - bridge_faces) * 3;
    out->nfaces = bridge_faces;
    out->refined_faces = refined_nf;
    out->refined_vertices_added = nnew;
    out->combined_faces = combined_faces;
    return 0;
}

int AtlasXyzWeldTopology_collect_connections(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const int32_t *face_chart, size_t ncharts,
    AtlasXyzWeldConnection **out_connection, size_t *out_nconnection,
    AtlasXyzWeldConnectionStats *stats)
{
    if (arena == NULL || ps == NULL || topology == NULL ||
        topology->vertices == NULL || topology->faces == NULL ||
        topology->nfaces == 0 || topology->source_vertices != ps->nv ||
        topology->nvertices < ps->nv || face_chart == NULL || ncharts == 0 ||
        ncharts > (size_t)INT32_MAX || out_connection == NULL ||
        out_nconnection == NULL || stats == NULL || ps->faces == NULL ||
        ps->cube_voff == NULL || ps->n_cubes < 2 ||
        (topology->nvertices > ps->nv &&
         (topology->parent0 == NULL || topology->parent1 == NULL)) ||
        ncharts > SIZE_MAX / ncharts)
        return -1;
    memset(stats, 0, sizeof(*stats));
    *out_connection = NULL;
    *out_nconnection = 0;

    int32_t *vertex_chart = (int32_t *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(int32_t));
    int32_t *vertex_cube = (int32_t *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(int32_t));
    for (size_t i = 0; i < topology->nvertices; i++) {
        vertex_chart[i] = -1;
        vertex_cube[i] = -1;
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_chart[f];
        if (chart < 0 || (size_t)chart >= ncharts) return -1;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            if (vertex_chart[vertex] == -1) vertex_chart[vertex] = chart;
            else if (vertex_chart[vertex] != chart) vertex_chart[vertex] = -2;
        }
    }
    for (size_t vertex = 0; vertex < ps->nv; vertex++) {
        size_t cube = axwa_vertex_cube(ps, vertex);
        if (cube > (size_t)INT32_MAX) return -1;
        vertex_cube[vertex] = (int32_t)cube;
    }
    for (size_t vertex = ps->nv; vertex < topology->nvertices; vertex++) {
        size_t ni = vertex - ps->nv;
        int32_t parent0 = topology->parent0[ni];
        int32_t parent1 = topology->parent1[ni];
        if (parent0 < 0 || parent1 < 0 ||
            (size_t)parent0 >= vertex || (size_t)parent1 >= vertex ||
            vertex_cube[parent0] != vertex_cube[parent1])
            return -1;
        vertex_cube[vertex] = vertex_cube[parent0];
        vertex_chart[vertex] =
            vertex_chart[parent0] >= 0 &&
            vertex_chart[parent0] == vertex_chart[parent1]
                ? vertex_chart[parent0] : -2;
    }

    if (topology->nfaces > SIZE_MAX / 3) return -1;
    size_t edge_capacity = 3 * topology->nfaces;
    AxwaEdge *edge = (AxwaEdge *)ARENA_ALLOC(
        arena, (edge_capacity ? edge_capacity : 1) * sizeof(AxwaEdge));
    size_t nedge = 0;
    for (size_t f = 0; f < topology->nfaces; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = topology->faces[f * 3 + (size_t)k];
            int32_t b = topology->faces[f * 3 + (size_t)((k + 1) % 3)];
            if (a < 0 || b < 0 || (size_t)a >= topology->nvertices ||
                (size_t)b >= topology->nvertices || a == b)
                return -1;
            edge[nedge].lo = a < b ? a : b;
            edge[nedge].hi = a < b ? b : a;
            nedge++;
        }
    }
    qsort(edge, nedge, sizeof(AxwaEdge), axwa_compare_edge);
    size_t unique_edge = 0;
    for (size_t i = 0; i < nedge; i++) {
        if (i > 0 && edge[i].lo == edge[i - 1].lo &&
            edge[i].hi == edge[i - 1].hi)
            continue;
        edge[unique_edge++] = edge[i];
    }
    nedge = unique_edge;

    size_t pair_cells = ncharts * ncharts;
    size_t *pair_edges = (size_t *)ARENA_CALLOC(
        arena, pair_cells, sizeof(size_t));
    double *pair_length = (double *)ARENA_CALLOC(
        arena, pair_cells, sizeof(double));
    for (size_t i = 0; i < nedge; i++) {
        int32_t a = edge[i].lo, b = edge[i].hi;
        if (vertex_cube[a] == vertex_cube[b]) continue;
        stats->cross_cube_edges++;
        int32_t chart0 = vertex_chart[a], chart1 = vertex_chart[b];
        if (chart0 < 0 || chart1 < 0) {
            stats->ambiguous_chart_edges++;
            continue;
        }
        if (chart0 == chart1) {
            stats->same_chart_edges++;
            continue;
        }
        if (chart0 > chart1) {
            int32_t temporary = chart0;
            chart0 = chart1;
            chart1 = temporary;
        }
        double length = axwa_xyz_edge(topology->vertices, a, b);
        if (!(length > 0.0) || !isfinite(length)) return -1;
        size_t cell = (size_t)chart0 * ncharts + (size_t)chart1;
        pair_edges[cell]++;
        pair_length[cell] += length;
        stats->cross_chart_edges++;
    }
    size_t nconnection = 0;
    for (size_t chart0 = 0; chart0 < ncharts; chart0++)
        for (size_t chart1 = chart0 + 1; chart1 < ncharts; chart1++)
            if (pair_edges[chart0 * ncharts + chart1] != 0) nconnection++;
    AtlasXyzWeldConnection *connection =
        (AtlasXyzWeldConnection *)ARENA_ALLOC(
            arena, (nconnection ? nconnection : 1) *
                       sizeof(AtlasXyzWeldConnection));
    size_t output = 0;
    for (size_t chart0 = 0; chart0 < ncharts; chart0++) {
        for (size_t chart1 = chart0 + 1; chart1 < ncharts; chart1++) {
            size_t cell = chart0 * ncharts + chart1;
            if (pair_edges[cell] == 0) continue;
            connection[output].chart0 = (int32_t)chart0;
            connection[output].chart1 = (int32_t)chart1;
            connection[output].cross_cube_edges = pair_edges[cell];
            connection[output].total_xyz_edge_length = pair_length[cell];
            output++;
        }
    }
    if (output != nconnection) return -1;
    stats->relations = nconnection;
    *out_connection = connection;
    *out_nconnection = nconnection;
    return 0;
}

int AtlasXyzWeldTopology_collect_shift_constraints(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const double *base_u, const int32_t *face_chart, size_t ncharts,
    AtlasXyzWeldShiftConstraint **out_constraint, size_t *out_nconstraint,
    AtlasXyzWeldConnectionStats *stats)
{
    if (arena == NULL || ps == NULL || topology == NULL || base_u == NULL ||
        topology->vertices == NULL || topology->faces == NULL ||
        topology->nfaces == 0 || topology->source_vertices != ps->nv ||
        topology->nvertices < ps->nv || face_chart == NULL || ncharts == 0 ||
        ncharts > (size_t)INT32_MAX || out_constraint == NULL ||
        out_nconstraint == NULL || stats == NULL || ps->faces == NULL ||
        ps->cube_voff == NULL || ps->n_cubes < 2 ||
        (topology->nvertices > ps->nv &&
         (topology->parent0 == NULL || topology->parent1 == NULL)))
        return -1;
    memset(stats, 0, sizeof(*stats));
    *out_constraint = NULL;
    *out_nconstraint = 0;

    int32_t *vertex_chart = (int32_t *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(*vertex_chart));
    int32_t *vertex_cube = (int32_t *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(*vertex_cube));
    double *vertex_u = (double *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(*vertex_u));
    for (size_t i = 0; i < topology->nvertices; i++) {
        vertex_chart[i] = -1;
        vertex_cube[i] = -1;
        vertex_u[i] = 0.0;
    }
    for (size_t i = 0; i < ps->nv; i++) {
        if (!isfinite(base_u[i])) return -1;
        vertex_u[i] = base_u[i];
        size_t cube = axwa_vertex_cube(ps, i);
        if (cube > (size_t)INT32_MAX) return -1;
        vertex_cube[i] = (int32_t)cube;
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_chart[f];
        if (chart < 0 || (size_t)chart >= ncharts) return -1;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            if (vertex_chart[vertex] == -1) vertex_chart[vertex] = chart;
            else if (vertex_chart[vertex] != chart) vertex_chart[vertex] = -2;
        }
    }
    for (size_t vertex = ps->nv; vertex < topology->nvertices; vertex++) {
        size_t offset = vertex - ps->nv;
        int32_t parent0 = topology->parent0[offset];
        int32_t parent1 = topology->parent1[offset];
        if (parent0 < 0 || parent1 < 0 ||
            (size_t)parent0 >= vertex || (size_t)parent1 >= vertex ||
            vertex_cube[parent0] != vertex_cube[parent1])
            return -1;
        vertex_cube[vertex] = vertex_cube[parent0];
        vertex_chart[vertex] =
            vertex_chart[parent0] >= 0 &&
            vertex_chart[parent0] == vertex_chart[parent1]
                ? vertex_chart[parent0] : -2;
        vertex_u[vertex] = 0.5 * (vertex_u[parent0] + vertex_u[parent1]);
    }

    if (topology->nfaces > SIZE_MAX / 3) return -1;
    size_t edge_capacity = 3 * topology->nfaces;
    AxwaEdge *edge = (AxwaEdge *)ARENA_ALLOC(
        arena, (edge_capacity ? edge_capacity : 1) * sizeof(*edge));
    size_t nedge = 0;
    for (size_t f = 0; f < topology->nfaces; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = topology->faces[f * 3 + (size_t)k];
            int32_t b = topology->faces[f * 3 + (size_t)((k + 1) % 3)];
            if (a < 0 || b < 0 || (size_t)a >= topology->nvertices ||
                (size_t)b >= topology->nvertices || a == b)
                return -1;
            edge[nedge].lo = a < b ? a : b;
            edge[nedge].hi = a < b ? b : a;
            nedge++;
        }
    }
    qsort(edge, nedge, sizeof(*edge), axwa_compare_edge);
    size_t unique_edge = 0;
    for (size_t i = 0; i < nedge; i++) {
        if (i > 0 && edge[i].lo == edge[i - 1].lo &&
            edge[i].hi == edge[i - 1].hi)
            continue;
        edge[unique_edge++] = edge[i];
    }
    nedge = unique_edge;

    AxwaShiftObservation *observation =
        (AxwaShiftObservation *)ARENA_ALLOC(
            arena, (nedge ? nedge : 1) * sizeof(*observation));
    size_t nobservation = 0;
    for (size_t i = 0; i < nedge; i++) {
        int32_t a = edge[i].lo;
        int32_t b = edge[i].hi;
        if (vertex_cube[a] == vertex_cube[b]) continue;
        stats->cross_cube_edges++;
        int32_t chart0 = vertex_chart[a];
        int32_t chart1 = vertex_chart[b];
        if (chart0 < 0 || chart1 < 0) {
            stats->ambiguous_chart_edges++;
            continue;
        }
        if (chart0 == chart1) {
            stats->same_chart_edges++;
            continue;
        }
        int32_t vertex0 = a;
        int32_t vertex1 = b;
        if (chart0 > chart1) {
            int32_t temporary = chart0;
            chart0 = chart1;
            chart1 = temporary;
            temporary = vertex0;
            vertex0 = vertex1;
            vertex1 = temporary;
        }
        double length = axwa_xyz_edge(topology->vertices, vertex0, vertex1);
        if (!(length > 0.0) || !isfinite(length)) return -1;
        observation[nobservation].chart0 = chart0;
        observation[nobservation].chart1 = chart1;
        observation[nobservation].target =
            vertex_u[vertex0] - vertex_u[vertex1];
        observation[nobservation].length = length;
        nobservation++;
        stats->cross_chart_edges++;
    }
    qsort(observation, nobservation, sizeof(*observation),
          axwa_compare_shift_observation);

    size_t nconstraint = 0;
    for (size_t first = 0; first < nobservation; ) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].chart0 == observation[first].chart0 &&
               observation[last].chart1 == observation[first].chart1)
            last++;
        nconstraint++;
        first = last;
    }
    AtlasXyzWeldShiftConstraint *constraint =
        (AtlasXyzWeldShiftConstraint *)ARENA_ALLOC(
            arena, (nconstraint ? nconstraint : 1) * sizeof(*constraint));
    double *deviation = (double *)ARENA_ALLOC(
        arena, (nobservation ? nobservation : 1) * sizeof(*deviation));
    size_t output = 0;
    for (size_t first = 0; first < nobservation; ) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].chart0 == observation[first].chart0 &&
               observation[last].chart1 == observation[first].chart1)
            last++;
        size_t count = last - first;
        double median = observation[first + count / 2].target;
        if ((count & 1u) == 0)
            median = 0.5 * (median +
                            observation[first + count / 2 - 1].target);
        double total_length = 0.0;
        for (size_t i = first; i < last; i++) {
            deviation[i - first] = fabs(observation[i].target - median);
            total_length += observation[i].length;
        }
        qsort(deviation, count, sizeof(*deviation), axwa_compare_double);
        double mad = deviation[count / 2];
        if ((count & 1u) == 0)
            mad = 0.5 * (mad + deviation[count / 2 - 1]);
        constraint[output].chart0 = observation[first].chart0;
        constraint[output].chart1 = observation[first].chart1;
        constraint[output].cross_cube_edges = count;
        constraint[output].total_xyz_edge_length = total_length;
        constraint[output].target_shift_median = median;
        constraint[output].target_shift_mad = mad;
        output++;
        first = last;
    }
    if (output != nconstraint) return -1;
    stats->relations = nconstraint;
    *out_constraint = constraint;
    *out_nconstraint = nconstraint;
    return 0;
}

int AtlasXyzWeldAudit_evaluate_with_mesh(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const double *base_u, const double *v, const int32_t *face_chart,
    size_t ncharts, const double *chart_shift,
    AtlasXyzWeldAuditStats *stats, AtlasXyzWeldAuditMesh *mesh)
{
    if (arena == NULL || ps == NULL || topology == NULL ||
        topology->vertices == NULL || topology->faces == NULL ||
        topology->nfaces == 0 || topology->nvertices < ps->nv ||
        topology->source_vertices != ps->nv ||
        (topology->nvertices > ps->nv &&
         (topology->parent0 == NULL || topology->parent1 == NULL)) ||
        base_u == NULL || v == NULL || face_chart == NULL ||
        ncharts == 0 || ncharts > (size_t)INT32_MAX ||
        chart_shift == NULL || stats == NULL || ps->faces == NULL ||
        ps->verts == NULL || ps->uv == NULL || ps->phi == NULL ||
        ps->cube_voff == NULL || ps->n_cubes < 2 ||
        ps->nf > (size_t)INT32_MAX || topology->nfaces > (size_t)INT32_MAX ||
        ps->nf > SIZE_MAX - topology->nfaces)
        return -1;
    memset(stats, 0, sizeof(*stats));
    if (mesh != NULL) memset(mesh, 0, sizeof(*mesh));
    stats->bridge_faces = topology->nfaces;

    if (ps->nf > SIZE_MAX / 3) return -1;
    size_t source_position_capacity = 3 * ps->nf;
    AxwaPosition *source_position = (AxwaPosition *)ARENA_ALLOC(
        arena, (source_position_capacity ? source_position_capacity : 1) *
                   sizeof(AxwaPosition));
    size_t nsource_position = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_chart[f];
        if (chart < 0 || (size_t)chart >= ncharts) return -1;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv ||
                !isfinite(base_u[vertex]) || !isfinite(v[vertex]) ||
                !isfinite(chart_shift[chart]))
                return -1;
            source_position[nsource_position].vertex = vertex;
            source_position[nsource_position].chart = chart;
            source_position[nsource_position].u =
                base_u[vertex] + chart_shift[chart];
            nsource_position++;
        }
    }
    qsort(source_position, nsource_position, sizeof(AxwaPosition),
          axwa_compare_position);
    size_t unique_position = 0;
    for (size_t i = 0; i < nsource_position; i++) {
        if (i > 0 &&
            source_position[i].vertex == source_position[i - 1].vertex &&
            source_position[i].chart == source_position[i - 1].chart)
            continue;
        source_position[unique_position++] = source_position[i];
    }
    nsource_position = unique_position;
    size_t *source_position_offset = (size_t *)ARENA_CALLOC(
        arena, ps->nv + 1, sizeof(size_t));
    for (size_t i = 0; i < nsource_position; i++)
        source_position_offset[(size_t)source_position[i].vertex + 1]++;
    for (size_t i = 0; i < ps->nv; i++)
        source_position_offset[i + 1] += source_position_offset[i];

    size_t nnew = topology->nvertices - ps->nv;
    if (nnew > (SIZE_MAX - nsource_position - 1) / 4) return -1;
    size_t position_capacity = nsource_position + 4 * nnew + 1;
    AxwaPosition *position = (AxwaPosition *)ARENA_ALLOC(
        arena, position_capacity * sizeof(AxwaPosition));
    size_t *position_offset = (size_t *)ARENA_ALLOC(
        arena, (topology->nvertices + 1) * sizeof(size_t));
    double *topology_base_u = (double *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(double));
    double *topology_v = (double *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(double));
    float *topology_registered = (float *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(float));
    float *topology_phi = (float *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(float));
    int32_t *topology_cube = (int32_t *)ARENA_ALLOC(
        arena, topology->nvertices * sizeof(int32_t));
    size_t nposition = 0;
    for (size_t vertex = 0; vertex < ps->nv; vertex++) {
        position_offset[vertex] = nposition;
        for (size_t i = source_position_offset[vertex];
             i < source_position_offset[vertex + 1]; i++)
            position[nposition++] = source_position[i];
        position_offset[vertex + 1] = nposition;
        topology_base_u[vertex] = base_u[vertex];
        topology_v[vertex] = v[vertex];
        topology_registered[vertex] = ps->uv[vertex * 2];
        topology_phi[vertex] = ps->phi[vertex];
        size_t cube = axwa_vertex_cube(ps, vertex);
        if (cube > (size_t)INT32_MAX) return -1;
        topology_cube[vertex] = (int32_t)cube;
    }
    const double two_pi = 6.283185307179586476925286766559;
    for (size_t vertex = ps->nv; vertex < topology->nvertices; vertex++) {
        size_t ni = vertex - ps->nv;
        int32_t parent0 = topology->parent0[ni];
        int32_t parent1 = topology->parent1[ni];
        if (parent0 < 0 || parent1 < 0 ||
            (size_t)parent0 >= vertex || (size_t)parent1 >= vertex ||
            topology_cube[parent0] != topology_cube[parent1])
            return -1;
        position_offset[vertex] = nposition;
        topology_base_u[vertex] = 0.5 *
            (topology_base_u[parent0] + topology_base_u[parent1]);
        size_t i0 = position_offset[parent0];
        size_t i1 = position_offset[parent1];
        size_t end0 = position_offset[(size_t)parent0 + 1];
        size_t end1 = position_offset[(size_t)parent1 + 1];
        while (i0 < end0 || i1 < end1) {
            int32_t chart;
            if (i1 >= end1 ||
                (i0 < end0 && position[i0].chart < position[i1].chart)) {
                chart = position[i0++].chart;
            } else if (i0 >= end0 || position[i1].chart < position[i0].chart) {
                chart = position[i1++].chart;
            } else {
                chart = position[i0].chart;
                i0++; i1++;
            }
            if (nposition == position_capacity) {
                if (position_capacity > SIZE_MAX / 2) return -1;
                size_t next_capacity = 2 * position_capacity;
                AxwaPosition *next = (AxwaPosition *)ARENA_ALLOC(
                    arena, next_capacity * sizeof(AxwaPosition));
                memcpy(next, position, nposition * sizeof(AxwaPosition));
                position = next;
                position_capacity = next_capacity;
            }
            position[nposition].vertex = (int32_t)vertex;
            position[nposition].chart = chart;
            position[nposition].u =
                topology_base_u[vertex] + chart_shift[chart];
            nposition++;
        }
        position_offset[vertex + 1] = nposition;
        if (position_offset[vertex] == nposition) return -1;
        topology_v[vertex] =
            0.5 * (topology_v[parent0] + topology_v[parent1]);
        topology_registered[vertex] = 0.5f *
            (topology_registered[parent0] + topology_registered[parent1]);
        double dphi = remainder(
            (double)topology_phi[parent1] - topology_phi[parent0], two_pi);
        topology_phi[vertex] = (float)(topology_phi[parent0] + 0.5 * dphi);
        topology_cube[vertex] = topology_cube[parent0];
    }
    stats->incident_chart_positions = nposition;

    uint8_t *bridge_vertex = (uint8_t *)ARENA_CALLOC(
        arena, topology->nvertices, sizeof(uint8_t));
    if (topology->nfaces > SIZE_MAX / 3) return -1;
    size_t edge_capacity = 3 * topology->nfaces;
    AxwaEdge *edge = (AxwaEdge *)ARENA_ALLOC(
        arena, (edge_capacity ? edge_capacity : 1) * sizeof(AxwaEdge));
    size_t nedge = 0;
    for (size_t f = 0; f < topology->nfaces; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = topology->faces[f * 3 + (size_t)k];
            int32_t b = topology->faces[f * 3 + (size_t)((k + 1) % 3)];
            if (a < 0 || b < 0 || (size_t)a >= topology->nvertices ||
                (size_t)b >= topology->nvertices || a == b)
                return -1;
            if (position_offset[(size_t)a] == position_offset[(size_t)a + 1] ||
                position_offset[(size_t)b] == position_offset[(size_t)b + 1])
                return -1;
            bridge_vertex[a] = bridge_vertex[b] = 1;
            edge[nedge].lo = a < b ? a : b;
            edge[nedge].hi = a < b ? b : a;
            nedge++;
        }
    }
    qsort(edge, nedge, sizeof(AxwaEdge), axwa_compare_edge);
    size_t unique_edge = 0;
    for (size_t i = 0; i < nedge; i++) {
        if (i > 0 && edge[i].lo == edge[i - 1].lo &&
            edge[i].hi == edge[i - 1].hi)
            continue;
        edge[unique_edge++] = edge[i];
    }
    nedge = unique_edge;
    stats->bridge_edges = nedge;

    double *vertex_spread = (double *)ARENA_ALLOC(
        arena, (topology->nvertices ? topology->nvertices : 1) *
                   sizeof(double));
    size_t nvertex_spread = 0;
    for (size_t vertex = 0; vertex < topology->nvertices; vertex++) {
        if (!bridge_vertex[vertex]) continue;
        stats->bridge_vertices++;
        size_t first = position_offset[vertex], last = position_offset[vertex + 1];
        double lo = HUGE_VAL, hi = -HUGE_VAL;
        for (size_t i = first; i < last; i++) {
            if (position[i].u < lo) lo = position[i].u;
            if (position[i].u > hi) hi = position[i].u;
        }
        double spread = hi - lo;
        vertex_spread[nvertex_spread++] = spread;
        if (last - first > 1) stats->split_bridge_vertices++;
    }
    stats->incident_u_spread_median =
        axwa_quantile(vertex_spread, nvertex_spread, 0.50);
    stats->incident_u_spread_p95 =
        axwa_quantile(vertex_spread, nvertex_spread, 0.95);
    stats->incident_u_spread_max =
        axwa_quantile(vertex_spread, nvertex_spread, 1.00);

    double *edge_stretch = (double *)ARENA_ALLOC(
        arena, (nedge ? nedge : 1) * sizeof(double));
    double *edge_error = (double *)ARENA_ALLOC(
        arena, (nedge ? nedge : 1) * sizeof(double));
    size_t ncross = 0;
    for (size_t e = 0; e < nedge; e++) {
        int32_t a = edge[e].lo, b = edge[e].hi;
        if (topology_cube[a] == topology_cube[b])
            continue;
        double xyz = axwa_xyz_edge(topology->vertices, a, b);
        if (!(xyz > 1.0e-12) || !isfinite(xyz)) return -1;
        double best_score = HUGE_VAL, best_uv = HUGE_VAL;
        for (size_t ia = position_offset[a]; ia < position_offset[(size_t)a + 1]; ia++)
            for (size_t ib = position_offset[b]; ib < position_offset[(size_t)b + 1]; ib++) {
                double uv = axwa_uv_edge(
                    position[ia].u, topology_v[a],
                    position[ib].u, topology_v[b]);
                double score = axwa_symmetric_stretch(uv, xyz);
                if (score < best_score ||
                    (score == best_score && uv < best_uv)) {
                    best_score = score;
                    best_uv = uv;
                }
            }
        if (!isfinite(best_score) || !isfinite(best_uv)) return -1;
        double stretch = best_uv / xyz;
        edge_stretch[ncross] = stretch;
        edge_error[ncross] = fabs(best_uv - xyz);
        if (best_uv > fmax(4.0 * xyz, 25.0))
            stats->exploded_cross_cube_edges++;
        if (best_uv < 0.25 * xyz)
            stats->compressed_cross_cube_edges++;
        ncross++;
    }
    stats->cross_cube_bridge_edges = ncross;
    stats->best_cross_cube_edge_stretch_median =
        axwa_quantile(edge_stretch, ncross, 0.50);
    stats->best_cross_cube_edge_stretch_p95 =
        axwa_quantile(edge_stretch, ncross, 0.95);
    stats->best_cross_cube_edge_stretch_max =
        axwa_quantile(edge_stretch, ncross, 1.00);
    stats->best_cross_cube_edge_abs_error_median =
        axwa_quantile(edge_error, ncross, 0.50);
    stats->best_cross_cube_edge_abs_error_p95 =
        axwa_quantile(edge_error, ncross, 0.95);
    stats->best_cross_cube_edge_abs_error_max =
        axwa_quantile(edge_error, ncross, 1.00);

    size_t total_faces = ps->nf + topology->nfaces;
    if (total_faces > (size_t)INT32_MAX / 3 ||
        ncharts > (size_t)INT32_MAX - topology->nfaces)
        return -1;
    size_t total_corners = 3 * total_faces;
    int32_t *audit_faces = (int32_t *)ARENA_ALLOC(
        arena, total_corners * sizeof(int32_t));
    double *audit_u = (double *)ARENA_ALLOC(
        arena, total_corners * sizeof(double));
    double *audit_v = (double *)ARENA_ALLOC(
        arena, total_corners * sizeof(double));
    float *audit_registered = (float *)ARENA_ALLOC(
        arena, total_corners * sizeof(float));
    float *audit_phi = (float *)ARENA_ALLOC(
        arena, total_corners * sizeof(float));
    int32_t *audit_component = (int32_t *)ARENA_ALLOC(
        arena, total_corners * sizeof(int32_t));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_chart[f];
        for (int k = 0; k < 3; k++) {
            size_t corner = f * 3 + (size_t)k;
            int32_t source = ps->faces[corner];
            audit_faces[corner] = (int32_t)corner;
            audit_u[corner] = base_u[source] + chart_shift[chart];
            audit_v[corner] = v[source];
            audit_registered[corner] = ps->uv[(size_t)source * 2];
            audit_phi[corner] = ps->phi[source];
            audit_component[corner] = chart;
        }
    }

    double *face_stretch = (double *)ARENA_ALLOC(
        arena, topology->nfaces * sizeof(double));
    AxwaEmbeddingPosition *embedded = (AxwaEmbeddingPosition *)ARENA_ALLOC(
        arena, 3 * topology->nfaces * sizeof(AxwaEmbeddingPosition));
    size_t nembedded = 0;
    for (size_t f = 0; f < topology->nfaces; f++) {
        int32_t source[3] = {
            topology->faces[f * 3], topology->faces[f * 3 + 1],
            topology->faces[f * 3 + 2]
        };
        double xyz_edge[3] = {
            axwa_xyz_edge(topology->vertices, source[0], source[1]),
            axwa_xyz_edge(topology->vertices, source[1], source[2]),
            axwa_xyz_edge(topology->vertices, source[2], source[0])
        };
        if (!(xyz_edge[0] > 1.0e-12) || !(xyz_edge[1] > 1.0e-12) ||
            !(xyz_edge[2] > 1.0e-12))
            return -1;
        double best_score = HUGE_VAL, best_sum = HUGE_VAL;
        size_t best_i[3] = {position_offset[source[0]],
                            position_offset[source[1]],
                            position_offset[source[2]]};
        for (size_t i0 = position_offset[source[0]];
             i0 < position_offset[(size_t)source[0] + 1]; i0++)
            for (size_t i1 = position_offset[source[1]];
                 i1 < position_offset[(size_t)source[1] + 1]; i1++)
                for (size_t i2 = position_offset[source[2]];
                     i2 < position_offset[(size_t)source[2] + 1]; i2++) {
                    double uv_edge[3] = {
                        axwa_uv_edge(position[i0].u, topology_v[source[0]],
                                     position[i1].u, topology_v[source[1]]),
                        axwa_uv_edge(position[i1].u, topology_v[source[1]],
                                     position[i2].u, topology_v[source[2]]),
                        axwa_uv_edge(position[i2].u, topology_v[source[2]],
                                     position[i0].u, topology_v[source[0]])
                    };
                    double score = 1.0, sum = 0.0;
                    for (int e = 0; e < 3; e++) {
                        double s = axwa_symmetric_stretch(
                            uv_edge[e], xyz_edge[e]);
                        if (s > score) score = s;
                        double l = isfinite(s) ? log(s) : HUGE_VAL;
                        sum += l * l;
                    }
                    if (score < best_score ||
                        (score == best_score && sum < best_sum)) {
                        best_score = score;
                        best_sum = sum;
                        best_i[0] = i0; best_i[1] = i1; best_i[2] = i2;
                    }
                }
        if (!isfinite(best_score)) return -1;
        face_stretch[f] = best_score;
        size_t output_face = ps->nf + f;
        double uv_edge2[3];
        for (int k = 0; k < 3; k++) {
            size_t corner = output_face * 3 + (size_t)k;
            int32_t vertex = source[k];
            double selected_u = position[best_i[k]].u;
            audit_faces[corner] = (int32_t)corner;
            audit_u[corner] = selected_u;
            audit_v[corner] = topology_v[vertex];
            audit_registered[corner] = topology_registered[vertex];
            audit_phi[corner] = topology_phi[vertex];
            audit_component[corner] = (int32_t)(ncharts + f);
            embedded[nembedded].vertex = vertex;
            embedded[nembedded].u = selected_u;
            nembedded++;
        }
        int exploded = 0;
        for (int e = 0; e < 3; e++) {
            int next = (e + 1) % 3;
            double uv = axwa_uv_edge(
                position[best_i[e]].u, topology_v[source[e]],
                position[best_i[next]].u, topology_v[source[next]]);
            uv_edge2[e] = uv * uv;
            if (uv > fmax(4.0 * xyz_edge[e], 25.0)) exploded = 1;
        }
        if (exploded) stats->exploded_bridge_faces++;
        double area2 = fabs(
            (audit_u[output_face * 3 + 1] - audit_u[output_face * 3]) *
            (audit_v[output_face * 3 + 2] - audit_v[output_face * 3]) -
            (audit_u[output_face * 3 + 2] - audit_u[output_face * 3]) *
            (audit_v[output_face * 3 + 1] - audit_v[output_face * 3]));
        double aspect = area2 > 2.0e-12
            ? (uv_edge2[0] + uv_edge2[1] + uv_edge2[2]) /
                  (2.0 * sqrt(3.0) * area2)
            : DBL_MAX;
        if (!isfinite(aspect) || aspect > 100.0)
            stats->degenerate_bridge_faces++;
    }
    stats->best_bridge_face_symmetric_stretch_median =
        axwa_quantile(face_stretch, topology->nfaces, 0.50);
    stats->best_bridge_face_symmetric_stretch_p95 =
        axwa_quantile(face_stretch, topology->nfaces, 0.95);
    stats->best_bridge_face_symmetric_stretch_max =
        axwa_quantile(face_stretch, topology->nfaces, 1.00);

    qsort(embedded, nembedded, sizeof(AxwaEmbeddingPosition),
          axwa_compare_embedding);
    double *embedded_spread = (double *)ARENA_ALLOC(
        arena, (nembedded ? nembedded : 1) * sizeof(double));
    size_t nembedded_spread = 0;
    double continuity_epsilon = 1.0e-6;
    for (size_t first = 0; first < nembedded;) {
        size_t last = first + 1;
        while (last < nembedded &&
               embedded[last].vertex == embedded[first].vertex)
            last++;
        double spread = embedded[last - 1].u - embedded[first].u;
        embedded_spread[nembedded_spread++] = spread;
        if (spread > continuity_epsilon)
            stats->discontinuous_embedded_bridge_vertices++;
        first = last;
    }
    stats->embedded_bridge_vertices = nembedded_spread;
    stats->embedded_u_spread_median =
        axwa_quantile(embedded_spread, nembedded_spread, 0.50);
    stats->embedded_u_spread_p95 =
        axwa_quantile(embedded_spread, nembedded_spread, 0.95);
    stats->embedded_u_spread_max =
        axwa_quantile(embedded_spread, nembedded_spread, 1.00);

    AtlasOverlapAudit overlap;
    size_t ncomponents = ncharts + topology->nfaces;
    if (AtlasOverlapAudit_build(
            arena, audit_faces, total_faces, total_corners,
            audit_u, audit_v, audit_registered, audit_phi,
            audit_component, ncomponents, &overlap) != 0 ||
        !overlap.broad_phase_complete ||
        overlap.indexed_faces != total_faces)
        return -1;
    stats->exact_overlap_pairs = overlap.exact_face_pairs;
    stats->broad_phase_records = overlap.broad_phase_records;
    stats->broad_phase_cells = overlap.broad_phase_cells;
    stats->broad_phase_candidate_pairs =
        overlap.broad_phase_candidate_pairs;
    stats->broad_phase_cell_size = overlap.broad_phase_cell_size;
    stats->broad_phase_complete = overlap.broad_phase_complete;
    for (size_t i = 0; i < overlap.npairs; i++) {
        int base0 = (size_t)overlap.pairs[i].face0 < ps->nf;
        int base1 = (size_t)overlap.pairs[i].face1 < ps->nf;
        if (base0 && base1) stats->base_base_overlap_pairs++;
        else if (base0 || base1) stats->bridge_base_overlap_pairs++;
        else stats->bridge_bridge_overlap_pairs++;
    }
    if (mesh != NULL) {
        mesh->faces = audit_faces;
        mesh->u = audit_u;
        mesh->v = audit_v;
        mesh->component = audit_component;
        mesh->nvertices = total_corners;
        mesh->nfaces = total_faces;
        mesh->base_faces = ps->nf;
        mesh->bridge_faces = topology->nfaces;
    }
    return 0;
}

int AtlasXyzWeldAudit_evaluate(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const double *base_u, const double *v, const int32_t *face_chart,
    size_t ncharts, const double *chart_shift,
    AtlasXyzWeldAuditStats *stats)
{
    return AtlasXyzWeldAudit_evaluate_with_mesh(
        arena, ps, topology, base_u, v, face_chart, ncharts, chart_shift,
        stats, NULL);
}
