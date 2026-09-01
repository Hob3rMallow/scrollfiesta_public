#define _USE_MATH_DEFINES
#include "winding_field.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define WF_INV_FOUR_PI (1.0 / (4.0 * M_PI))
#define WF_DEFAULT_LEAF_FACES 32u
#define WF_DEFAULT_BETA 2.0
#define WF_DEFAULT_EXACT_WORK UINT64_C(50000000)
#define WF_MORTON_BITS 10u
#define WF_CHILDREN 8u

typedef struct {
    uint32_t code;
    int32_t face;
} WfMortonFace;

typedef struct {
    uint64_t key;
    int8_t direction;
} WfHalfEdge;

typedef struct {
    int32_t lo, hi;
    int32_t multiplicity; /* signed net lo->hi incidence in the source */
} WfBoundaryEdge;

typedef struct {
    double center[3];
    double radius;
    double e0[3];
    double e1[9];
    double e2[27];
    float bmin[3], bmax[3];
    size_t face_first, face_count;
    size_t child_first;
    uint8_t child_count;
} WfNode;

struct WindingField_T {
    struct WindingField_T *self;
    const float *vertices;
    const int32_t *faces;
    size_t nvertices, nfaces;
    int32_t *face_order;
    WfNode *node;
    size_t nnodes, nleaves, root;
    WfBoundaryEdge *boundary;
    size_t nboundary, boundary_multiplicity;
    WindingFieldOptions options;
    double bbox_min[3], bbox_max[3], bbox_diagonal;
};

typedef struct {
    size_t exact_triangles;
    size_t approximate_nodes;
    int degenerate_ray;
} WfQueryStats;

#ifdef _OPENMP
static int wf_thread_count(void)
{
    int processors = omp_get_num_procs();
    return processors > 1 ? processors / 2 : 1;
}
#endif

static int wf_valid(WindingField_T field)
{
    return field != NULL && field->self == field &&
           field->vertices != NULL && field->faces != NULL &&
           field->nvertices > 0 && field->nfaces > 0 &&
           field->node != NULL && field->nnodes > 0;
}

void WindingField_default_options(WindingFieldOptions *options)
{
    if (options == NULL) return;
    options->backend = WINDING_FIELD_AUTO;
    options->leaf_faces = WF_DEFAULT_LEAF_FACES;
    options->beta = WF_DEFAULT_BETA;
    options->exact_work_limit = WF_DEFAULT_EXACT_WORK;
}

static int wf_compare_morton(const void *pa, const void *pb)
{
    const WfMortonFace *a = (const WfMortonFace *)pa;
    const WfMortonFace *b = (const WfMortonFace *)pb;
    if (a->code != b->code) return a->code < b->code ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

static int wf_compare_half_edge(const void *pa, const void *pb)
{
    const WfHalfEdge *a = (const WfHalfEdge *)pa;
    const WfHalfEdge *b = (const WfHalfEdge *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (int)a->direction - (int)b->direction;
}

static uint32_t wf_morton_spread(uint32_t value)
{
    value &= UINT32_C(0x000003ff);
    value = (value | (value << 16)) & UINT32_C(0x030000ff);
    value = (value | (value << 8)) & UINT32_C(0x0300f00f);
    value = (value | (value << 4)) & UINT32_C(0x030c30c3);
    value = (value | (value << 2)) & UINT32_C(0x09249249);
    return value;
}

static uint32_t wf_morton_code(const double p[3],
                               const double bmin[3],
                               const double bmax[3])
{
    uint32_t q[3] = { 0, 0, 0 };
    const double scale = (double)((UINT32_C(1) << WF_MORTON_BITS) - 1u);
    for (int k = 0; k < 3; k++) {
        double extent = bmax[k] - bmin[k];
        double t = extent > 0.0 ? (p[k] - bmin[k]) / extent : 0.5;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        q[k] = (uint32_t)floor(t * scale + 0.5);
    }
    return wf_morton_spread(q[0]) |
           (wf_morton_spread(q[1]) << 1) |
           (wf_morton_spread(q[2]) << 2);
}

static void wf_face_geometry(
    WindingField_T field, size_t face_index,
    double centroid[3], double area_normal[3], double *area)
{
    const int32_t *face = &field->faces[face_index * 3];
    const float *a = &field->vertices[(size_t)face[0] * 3];
    const float *b = &field->vertices[(size_t)face[1] * 3];
    const float *c = &field->vertices[(size_t)face[2] * 3];
    double e1[3] = { (double)b[0] - a[0], (double)b[1] - a[1],
                     (double)b[2] - a[2] };
    double e2[3] = { (double)c[0] - a[0], (double)c[1] - a[1],
                     (double)c[2] - a[2] };
    area_normal[0] = 0.5 * (e1[1] * e2[2] - e1[2] * e2[1]);
    area_normal[1] = 0.5 * (e1[2] * e2[0] - e1[0] * e2[2]);
    area_normal[2] = 0.5 * (e1[0] * e2[1] - e1[1] * e2[0]);
    *area = sqrt(area_normal[0] * area_normal[0] +
                 area_normal[1] * area_normal[1] +
                 area_normal[2] * area_normal[2]);
    for (int k = 0; k < 3; k++)
        centroid[k] = ((double)a[k] + b[k] + c[k]) / 3.0;
}

static int wf_validate_faces(
    const float *vertices, size_t nvertices,
    const int32_t *faces, size_t nfaces,
    double bmin[3], double bmax[3])
{
    for (int k = 0; k < 3; k++) {
        bmin[k] = DBL_MAX;
        bmax[k] = -DBL_MAX;
    }
    for (size_t i = 0; i < nvertices; i++) {
        for (int k = 0; k < 3; k++) {
            double value = vertices[i * 3 + (size_t)k];
            if (!isfinite(value)) return -1;
            if (value < bmin[k]) bmin[k] = value;
            if (value > bmax[k]) bmax[k] = value;
        }
    }
    for (size_t i = 0; i < nfaces * 3; i++) {
        if (faces[i] < 0 || (size_t)faces[i] >= nvertices) return -1;
    }
    return 0;
}

static size_t wf_node_capacity(size_t nfaces, size_t leaf_faces,
                               size_t *out_leaves)
{
    size_t level = (nfaces + leaf_faces - 1) / leaf_faces;
    size_t total = level;
    *out_leaves = level;
    while (level > 1) {
        level = (level + WF_CHILDREN - 1) / WF_CHILDREN;
        if (total > SIZE_MAX - level) return 0;
        total += level;
    }
    return total;
}

static void wf_compute_node(WindingField_T field, WfNode *node)
{
    double area_total = 0.0;
    double centroid_sum[3] = { 0.0, 0.0, 0.0 };
    double fallback_sum[3] = { 0.0, 0.0, 0.0 };
    memset(node->e0, 0, sizeof node->e0);
    memset(node->e1, 0, sizeof node->e1);
    memset(node->e2, 0, sizeof node->e2);
    for (int k = 0; k < 3; k++) {
        node->bmin[k] = FLT_MAX;
        node->bmax[k] = -FLT_MAX;
    }
    for (size_t at = node->face_first;
         at < node->face_first + node->face_count; at++) {
        size_t f = (size_t)field->face_order[at];
        double centroid[3], area_normal[3], area = 0.0;
        wf_face_geometry(field, f, centroid, area_normal, &area);
        area_total += area;
        for (int k = 0; k < 3; k++) {
            centroid_sum[k] += area * centroid[k];
            fallback_sum[k] += centroid[k];
        }
        for (int corner = 0; corner < 3; corner++) {
            const float *p = &field->vertices[
                (size_t)field->faces[f * 3 + (size_t)corner] * 3];
            for (int k = 0; k < 3; k++) {
                if (p[k] < node->bmin[k]) node->bmin[k] = p[k];
                if (p[k] > node->bmax[k]) node->bmax[k] = p[k];
            }
        }
    }
    for (int k = 0; k < 3; k++) {
        node->center[k] = area_total > 0.0
            ? centroid_sum[k] / area_total
            : fallback_sum[k] / (double)node->face_count;
    }
    node->radius = 0.0;
    for (size_t at = node->face_first;
         at < node->face_first + node->face_count; at++) {
        size_t f = (size_t)field->face_order[at];
        double centroid[3], area_normal[3], area = 0.0;
        double d[3];
        wf_face_geometry(field, f, centroid, area_normal, &area);
        (void)area;
        for (int k = 0; k < 3; k++) {
            d[k] = centroid[k] - node->center[k];
            node->e0[k] += area_normal[k];
        }
        for (int u = 0; u < 3; u++) {
            for (int v = 0; v < 3; v++) {
                node->e1[u * 3 + v] += d[u] * area_normal[v];
                for (int k = 0; k < 3; k++)
                    node->e2[k * 9 + u * 3 + v] +=
                        0.5 * d[k] * d[u] * area_normal[v];
            }
        }
        for (int corner = 0; corner < 3; corner++) {
            const float *p = &field->vertices[
                (size_t)field->faces[f * 3 + (size_t)corner] * 3];
            double r2 = 0.0;
            for (int k = 0; k < 3; k++) {
                double delta = (double)p[k] - node->center[k];
                r2 += delta * delta;
            }
            if (r2 > node->radius * node->radius)
                node->radius = sqrt(r2);
        }
    }
}

static int wf_build_bvh(Arena_T arena, WindingField_T field,
                        size_t leaf_faces)
{
    size_t nleaves = 0;
    size_t capacity = wf_node_capacity(field->nfaces, leaf_faces, &nleaves);
    if (capacity == 0 || capacity > (size_t)LONG_MAX / sizeof(WfNode) ||
        field->nfaces > (size_t)LONG_MAX / sizeof(int32_t) ||
        field->nfaces > (size_t)LONG_MAX / sizeof(WfMortonFace))
        return -1;
    field->face_order = (int32_t *)ARENA_ALLOC(
        arena, (long)(field->nfaces * sizeof(int32_t)));
    field->node = (WfNode *)ARENA_CALLOC(
        arena, (long)capacity, (long)sizeof(WfNode));
    field->nleaves = nleaves;
    Arena_Mark mark = Arena_save(arena);
    WfMortonFace *order = (WfMortonFace *)ARENA_ALLOC(
        arena, (long)(field->nfaces * sizeof(WfMortonFace)));
    for (size_t f = 0; f < field->nfaces; f++) {
        double centroid[3], area_normal[3], area = 0.0;
        wf_face_geometry(field, f, centroid, area_normal, &area);
        (void)area_normal;
        (void)area;
        order[f].code = wf_morton_code(
            centroid, field->bbox_min, field->bbox_max);
        order[f].face = (int32_t)f;
    }
    qsort(order, field->nfaces, sizeof(WfMortonFace), wf_compare_morton);
    for (size_t f = 0; f < field->nfaces; f++)
        field->face_order[f] = order[f].face;
    Arena_restore(arena, mark);

    size_t nnodes = 0;
    for (size_t leaf = 0; leaf < nleaves; leaf++) {
        WfNode *node = &field->node[nnodes++];
        node->face_first = leaf * leaf_faces;
        node->face_count = field->nfaces - node->face_first;
        if (node->face_count > leaf_faces) node->face_count = leaf_faces;
        node->child_count = 0;
        wf_compute_node(field, node);
    }
    size_t level_first = 0;
    size_t level_count = nleaves;
    while (level_count > 1) {
        size_t next_first = nnodes;
        for (size_t first = 0; first < level_count; first += WF_CHILDREN) {
            size_t count = level_count - first;
            if (count > WF_CHILDREN) count = WF_CHILDREN;
            WfNode *node = &field->node[nnodes++];
            WfNode *first_child = &field->node[level_first + first];
            WfNode *last_child = &field->node[level_first + first + count - 1];
            node->child_first = level_first + first;
            node->child_count = (uint8_t)count;
            node->face_first = first_child->face_first;
            node->face_count = last_child->face_first +
                               last_child->face_count - node->face_first;
            wf_compute_node(field, node);
        }
        level_first = next_first;
        level_count = nnodes - next_first;
    }
    if (nnodes != capacity) return -1;
    field->nnodes = nnodes;
    field->root = nnodes - 1;
    return 0;
}

static void wf_make_half_edges(WindingField_T field, WfHalfEdge *edge)
{
    size_t at = 0;
    for (size_t f = 0; f < field->nfaces; f++) {
        const int32_t *face = &field->faces[f * 3];
        for (int k = 0; k < 3; k++) {
            int32_t a = face[k];
            int32_t b = face[(k + 1) % 3];
            uint32_t lo = (uint32_t)(a < b ? a : b);
            uint32_t hi = (uint32_t)(a < b ? b : a);
            edge[at].key = ((uint64_t)lo << 32) | (uint64_t)hi;
            edge[at].direction = a < b ? 1 : -1;
            at++;
        }
    }
}

static int wf_count_boundaries(
    const WfHalfEdge *edge, size_t nedges,
    size_t *out_count, size_t *out_multiplicity)
{
    size_t count = 0, multiplicity = 0;
    for (size_t first = 0; first < nedges;) {
        size_t last = first + 1;
        int64_t signed_count = edge[first].direction;
        while (last < nedges && edge[last].key == edge[first].key) {
            signed_count += edge[last].direction;
            last++;
        }
        if (signed_count != 0) {
            uint64_t magnitude = signed_count < 0
                ? (uint64_t)(-signed_count) : (uint64_t)signed_count;
            if (magnitude > (uint64_t)INT32_MAX ||
                multiplicity > SIZE_MAX - (size_t)magnitude)
                return -1;
            count++;
            multiplicity += (size_t)magnitude;
        }
        first = last;
    }
    *out_count = count;
    *out_multiplicity = multiplicity;
    return 0;
}

static int wf_build_boundaries(Arena_T arena, WindingField_T field)
{
    if (field->nfaces > SIZE_MAX / 3) return -1;
    size_t nedges = field->nfaces * 3;
    if (nedges > (size_t)LONG_MAX / sizeof(WfHalfEdge)) return -1;
    size_t count = 0, multiplicity = 0;
    {
        Arena_Mark mark = Arena_save(arena);
        WfHalfEdge *edge = (WfHalfEdge *)ARENA_ALLOC(
            arena, (long)(nedges * sizeof(WfHalfEdge)));
        wf_make_half_edges(field, edge);
        qsort(edge, nedges, sizeof(WfHalfEdge), wf_compare_half_edge);
        if (wf_count_boundaries(edge, nedges, &count, &multiplicity) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        Arena_restore(arena, mark);
    }
    if (count > (size_t)LONG_MAX / sizeof(WfBoundaryEdge)) return -1;
    field->boundary = count > 0 ? (WfBoundaryEdge *)ARENA_ALLOC(
        arena, (long)(count * sizeof(WfBoundaryEdge))) : NULL;
    field->nboundary = count;
    field->boundary_multiplicity = multiplicity;
    if (count == 0) return 0;
    {
        Arena_Mark mark = Arena_save(arena);
        WfHalfEdge *edge = (WfHalfEdge *)ARENA_ALLOC(
            arena, (long)(nedges * sizeof(WfHalfEdge)));
        size_t out = 0;
        wf_make_half_edges(field, edge);
        qsort(edge, nedges, sizeof(WfHalfEdge), wf_compare_half_edge);
        for (size_t first = 0; first < nedges;) {
            size_t last = first + 1;
            int64_t signed_count = edge[first].direction;
            while (last < nedges && edge[last].key == edge[first].key) {
                signed_count += edge[last].direction;
                last++;
            }
            if (signed_count != 0) {
                field->boundary[out].lo =
                    (int32_t)(uint32_t)(edge[first].key >> 32);
                field->boundary[out].hi =
                    (int32_t)(uint32_t)edge[first].key;
                field->boundary[out].multiplicity = (int32_t)signed_count;
                out++;
            }
            first = last;
        }
        if (out != count) {
            Arena_restore(arena, mark);
            return -1;
        }
        Arena_restore(arena, mark);
    }
    return 0;
}

int WindingField_build(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const int32_t *faces, size_t nfaces,
    const WindingFieldOptions *options,
    WindingField_T *out_field,
    WindingFieldStats *stats)
{
    WindingFieldOptions selected;
    WindingField_T field = NULL;
    double diag2 = 0.0;
    if (arena == NULL || vertices == NULL || nvertices == 0 ||
        faces == NULL || nfaces == 0 || out_field == NULL)
        return -1;
    WindingField_default_options(&selected);
    if (options != NULL) selected = *options;
    if (selected.backend < WINDING_FIELD_AUTO ||
        selected.backend > WINDING_FIELD_FAST ||
        selected.leaf_faces == 0 || selected.leaf_faces > SIZE_MAX / 2 ||
        !(selected.beta > 0.0) || !isfinite(selected.beta))
        return -1;
    ARENA_NEW(arena, field);
    memset(field, 0, sizeof *field);
    field->self = field;
    field->vertices = vertices;
    field->faces = faces;
    field->nvertices = nvertices;
    field->nfaces = nfaces;
    field->options = selected;
    if (wf_validate_faces(vertices, nvertices, faces, nfaces,
                          field->bbox_min, field->bbox_max) != 0)
        return -1;
    for (int k = 0; k < 3; k++) {
        double d = field->bbox_max[k] - field->bbox_min[k];
        diag2 += d * d;
    }
    field->bbox_diagonal = sqrt(diag2);
    if (wf_build_bvh(arena, field, selected.leaf_faces) != 0 ||
        wf_build_boundaries(arena, field) != 0)
        return -1;
    if (stats != NULL) {
        memset(stats, 0, sizeof *stats);
        stats->vertices = nvertices;
        stats->faces = nfaces;
        stats->bvh_nodes = field->nnodes;
        stats->bvh_leaves = field->nleaves;
        stats->exterior_edges = field->nboundary;
        stats->exterior_multiplicity = field->boundary_multiplicity;
    }
    *out_field = field;
    return 0;
}

static double wf_solid_angle_points(
    const double a[3], const double b[3], const double c[3],
    const double q[3])
{
    double x[3] = { a[0] - q[0], a[1] - q[1], a[2] - q[2] };
    double y[3] = { b[0] - q[0], b[1] - q[1], b[2] - q[2] };
    double z[3] = { c[0] - q[0], c[1] - q[1], c[2] - q[2] };
    double lx = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double ly = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    double lz = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    double det = x[0]*(y[1]*z[2] - y[2]*z[1]) -
                 x[1]*(y[0]*z[2] - y[2]*z[0]) +
                 x[2]*(y[0]*z[1] - y[1]*z[0]);
    double denom = lx*ly*lz +
        (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])*lz +
        (y[0]*z[0] + y[1]*z[1] + y[2]*z[2])*lx +
        (z[0]*x[0] + z[1]*x[1] + z[2]*x[2])*ly;
    if (lx == 0.0 || ly == 0.0 || lz == 0.0) return NAN;
    return 2.0 * atan2(det, denom);
}

static double wf_face_solid_angle(
    WindingField_T field, size_t f, const double q[3])
{
    const int32_t *face = &field->faces[f * 3];
    double p[3][3];
    for (int corner = 0; corner < 3; corner++)
        for (int k = 0; k < 3; k++)
            p[corner][k] = field->vertices[
                (size_t)face[corner] * 3 + (size_t)k];
    return wf_solid_angle_points(p[0], p[1], p[2], q);
}

static double wf_eval_direct(
    WindingField_T field, const double q[3], WfQueryStats *stats)
{
    double total = 0.0;
    for (size_t f = 0; f < field->nfaces; f++) {
        double value = wf_face_solid_angle(field, f, q);
        if (!isfinite(value)) return NAN;
        total += value;
    }
    stats->exact_triangles += field->nfaces;
    return total * WF_INV_FOUR_PI;
}

static double wf_eval_expansion(const WfNode *node, const double q[3])
{
    double x[3] = { node->center[0] - q[0],
                    node->center[1] - q[1],
                    node->center[2] - q[2] };
    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    if (!(r2 > 0.0) || !isfinite(r2)) return NAN;
    double r = sqrt(r2);
    double inv_r3 = WF_INV_FOUR_PI / (r2 * r);
    double inv_r5 = inv_r3 / r2;
    double inv_r7 = inv_r5 / r2;
    double value = (x[0]*node->e0[0] + x[1]*node->e0[1] +
                    x[2]*node->e0[2]) * inv_r3;
    for (int u = 0; u < 3; u++) {
        for (int v = 0; v < 3; v++) {
            double h = -3.0 * x[u] * x[v] * inv_r5;
            if (u == v) h += inv_r3;
            value += h * node->e1[u * 3 + v];
        }
    }
    for (int k = 0; k < 3; k++) {
        for (int u = 0; u < 3; u++) {
            for (int v = 0; v < 3; v++) {
                double t = 15.0 * x[k] * x[u] * x[v] * inv_r7;
                double diagonal = 0.0;
                if (u == v) diagonal += x[k];
                if (u == k) diagonal += x[v];
                if (v == k) diagonal += x[u];
                t -= 3.0 * diagonal * inv_r5;
                value += t * node->e2[k * 9 + u * 3 + v];
            }
        }
    }
    return value;
}

static double wf_eval_fast_node(
    WindingField_T field, size_t node_index, const double q[3],
    double beta, WfQueryStats *stats)
{
    const WfNode *node = &field->node[node_index];
    double d0 = node->center[0] - q[0];
    double d1 = node->center[1] - q[1];
    double d2 = node->center[2] - q[2];
    double distance = sqrt(d0*d0 + d1*d1 + d2*d2);
    if (distance > beta * node->radius && distance > 0.0) {
        double value = wf_eval_expansion(node, q);
        if (isfinite(value)) {
            stats->approximate_nodes++;
            return value;
        }
    }
    if (node->child_count == 0) {
        double total = 0.0;
        for (size_t at = node->face_first;
             at < node->face_first + node->face_count; at++) {
            double value = wf_face_solid_angle(
                field, (size_t)field->face_order[at], q);
            if (!isfinite(value)) return NAN;
            total += value * WF_INV_FOUR_PI;
        }
        stats->exact_triangles += node->face_count;
        return total;
    }
    double total = 0.0;
    for (size_t child = 0; child < node->child_count; child++) {
        double value = wf_eval_fast_node(
            field, node->child_first + child, q, beta, stats);
        if (!isfinite(value)) return NAN;
        total += value;
    }
    return total;
}

static int wf_ray_triangle(
    WindingField_T field, size_t f, const double q[3], int *ambiguous)
{
    const int32_t *face = &field->faces[f * 3];
    const float *a = &field->vertices[(size_t)face[0] * 3];
    const float *b = &field->vertices[(size_t)face[1] * 3];
    const float *c = &field->vertices[(size_t)face[2] * 3];
    double denom = ((double)b[1] - c[1]) * ((double)a[0] - c[0]) +
                   ((double)c[0] - b[0]) * ((double)a[1] - c[1]);
    double scale = fabs((double)a[0]) + fabs((double)a[1]) +
                   fabs((double)b[0]) + fabs((double)b[1]) +
                   fabs((double)c[0]) + fabs((double)c[1]) + 1.0;
    double tol = 64.0 * DBL_EPSILON * scale;
    if (fabs(denom) <= tol) {
        double min0 = fmin((double)a[0], fmin((double)b[0], c[0]));
        double max0 = fmax((double)a[0], fmax((double)b[0], c[0]));
        double min1 = fmin((double)a[1], fmin((double)b[1], c[1]));
        double max1 = fmax((double)a[1], fmax((double)b[1], c[1]));
        if (q[0] >= min0 - tol && q[0] <= max0 + tol &&
            q[1] >= min1 - tol && q[1] <= max1 + tol)
            *ambiguous = 1;
        return 0;
    }
    double l0 = (((double)b[1] - c[1]) * (q[0] - c[0]) +
                 ((double)c[0] - b[0]) * (q[1] - c[1])) / denom;
    double l1 = (((double)c[1] - a[1]) * (q[0] - c[0]) +
                 ((double)a[0] - c[0]) * (q[1] - c[1])) / denom;
    double l2 = 1.0 - l0 - l1;
    double bary_tol = 128.0 * DBL_EPSILON;
    if (l0 < -bary_tol || l1 < -bary_tol || l2 < -bary_tol) return 0;
    if (l0 <= bary_tol || l1 <= bary_tol || l2 <= bary_tol) {
        *ambiguous = 1;
        return 0;
    }
    double hit = l0 * a[2] + l1 * b[2] + l2 * c[2];
    double hit_tol = 128.0 * DBL_EPSILON *
                     (fabs(hit) + fabs(q[2]) + 1.0);
    if (fabs(hit - q[2]) <= hit_tol) {
        *ambiguous = 1;
        return 0;
    }
    if (hit < q[2]) return 0;
    return denom > 0.0 ? 1 : -1;
}

static int64_t wf_ray_count_node(
    WindingField_T field, size_t node_index, const double q[3],
    int *ambiguous)
{
    const WfNode *node = &field->node[node_index];
    double tol = 64.0 * DBL_EPSILON *
        (fabs(q[0]) + fabs(q[1]) + fabs(q[2]) +
         field->bbox_diagonal + 1.0);
    if (q[0] < (double)node->bmin[0] - tol ||
        q[0] > (double)node->bmax[0] + tol ||
        q[1] < (double)node->bmin[1] - tol ||
        q[1] > (double)node->bmax[1] + tol ||
        q[2] > (double)node->bmax[2] + tol)
        return 0;
    if (node->child_count == 0) {
        int64_t count = 0;
        for (size_t at = node->face_first;
             at < node->face_first + node->face_count; at++) {
            count += wf_ray_triangle(
                field, (size_t)field->face_order[at], q, ambiguous);
            if (*ambiguous) return 0;
        }
        return count;
    }
    int64_t count = 0;
    for (size_t child = 0; child < node->child_count; child++) {
        count += wf_ray_count_node(
            field, node->child_first + child, q, ambiguous);
        if (*ambiguous) return 0;
    }
    return count;
}

static double wf_eval_boundary_exact(
    WindingField_T field, const double q[3], WfQueryStats *stats)
{
    int ambiguous = 0;
    int64_t ray_count = wf_ray_count_node(field, field->root, q, &ambiguous);
    if (ambiguous) {
        stats->degenerate_ray = 1;
        return wf_eval_direct(field, q, stats);
    }
    double distance = field->bbox_diagonal > 0.0
        ? 2.0 * field->bbox_diagonal + 1.0 : 1.0;
    double apex[3] = { q[0], q[1], q[2] - distance };
    double cone = 0.0;
    for (size_t i = 0; i < field->nboundary; i++) {
        const WfBoundaryEdge *edge = &field->boundary[i];
        double hi[3], lo[3];
        for (int k = 0; k < 3; k++) {
            hi[k] = field->vertices[(size_t)edge->hi * 3 + (size_t)k];
            lo[k] = field->vertices[(size_t)edge->lo * 3 + (size_t)k];
        }
        double angle = wf_solid_angle_points(hi, lo, apex, q);
        if (!isfinite(angle)) {
            stats->degenerate_ray = 1;
            return wf_eval_direct(field, q, stats);
        }
        cone += (double)edge->multiplicity * angle;
    }
    stats->exact_triangles += field->boundary_multiplicity;
    return (double)ray_count - cone * WF_INV_FOUR_PI;
}

static WindingFieldBackend wf_select_backend(
    WindingField_T field, WindingFieldBackend requested, size_t nqueries)
{
    WindingFieldBackend backend = requested;
    if (backend == WINDING_FIELD_AUTO) backend = field->options.backend;
    if (backend != WINDING_FIELD_AUTO) return backend;
    if (field->boundary_multiplicity == 0)
        return WINDING_FIELD_BOUNDARY_EXACT;
    if (nqueries == 0) return WINDING_FIELD_BOUNDARY_EXACT;
    if (field->boundary_multiplicity <=
        field->options.exact_work_limit / nqueries)
        return WINDING_FIELD_BOUNDARY_EXACT;
    return WINDING_FIELD_FAST;
}

static double wf_eval_one(
    WindingField_T field, const double q[3], WindingFieldBackend backend,
    WfQueryStats *stats)
{
    if (!isfinite(q[0]) || !isfinite(q[1]) || !isfinite(q[2])) return NAN;
    if (backend == WINDING_FIELD_DIRECT)
        return wf_eval_direct(field, q, stats);
    if (backend == WINDING_FIELD_BOUNDARY_EXACT)
        return wf_eval_boundary_exact(field, q, stats);
    if (backend == WINDING_FIELD_FAST)
        return wf_eval_fast_node(
            field, field->root, q, field->options.beta, stats);
    return NAN;
}

int WindingField_evaluate_many(
    WindingField_T field,
    const double *queries, size_t nqueries,
    WindingFieldBackend backend,
    double *out_winding,
    WindingFieldEvalStats *stats)
{
    if (!wf_valid(field) || (nqueries > 0 &&
        (queries == NULL || out_winding == NULL)))
        return -1;
    WindingFieldBackend selected = wf_select_backend(field, backend, nqueries);
    if (selected < WINDING_FIELD_DIRECT || selected > WINDING_FIELD_FAST)
        return -1;
    unsigned long long exact_triangles = 0;
    unsigned long long approximate_nodes = 0;
    unsigned long long degenerate_rays = 0;
    int64_t count = (int64_t)nqueries;
    int64_t qi = 0;
#ifdef _OPENMP
    int thread_count = wf_thread_count();
#pragma omp parallel for schedule(dynamic, 64) \
    num_threads(thread_count) \
    reduction(+:exact_triangles, approximate_nodes, degenerate_rays)
#endif
    for (qi = 0; qi < count; qi++) {
        WfQueryStats local;
        memset(&local, 0, sizeof local);
        out_winding[(size_t)qi] = wf_eval_one(
            field, &queries[(size_t)qi * 3], selected, &local);
        exact_triangles += (unsigned long long)local.exact_triangles;
        approximate_nodes += (unsigned long long)local.approximate_nodes;
        degenerate_rays += (unsigned long long)local.degenerate_ray;
    }
    if (stats != NULL) {
        memset(stats, 0, sizeof *stats);
        stats->backend_used = selected;
        stats->queries = nqueries;
        stats->direct_queries = selected == WINDING_FIELD_DIRECT ? nqueries : 0;
        stats->boundary_queries =
            selected == WINDING_FIELD_BOUNDARY_EXACT ? nqueries : 0;
        stats->fast_queries = selected == WINDING_FIELD_FAST ? nqueries : 0;
        stats->degenerate_ray_fallbacks = (size_t)degenerate_rays;
        stats->exact_triangle_evaluations = (size_t)exact_triangles;
        stats->approximate_node_evaluations = (size_t)approximate_nodes;
    }
    return 0;
}

int WindingField_evaluate_sides(
    WindingField_T field,
    const float *points, const float *normals, size_t npoints,
    double eps, WindingFieldBackend backend,
    float *out_mean, float *out_jump,
    WindingFieldEvalStats *stats)
{
    if (!wf_valid(field) || (npoints > 0 &&
        (points == NULL || normals == NULL || out_mean == NULL ||
         out_jump == NULL)) || !(eps > 0.0) || !isfinite(eps) ||
        npoints > (size_t)INT64_MAX)
        return -1;
    size_t query_count = npoints <= SIZE_MAX / 2 ? 2 * npoints : SIZE_MAX;
    WindingFieldBackend selected = wf_select_backend(
        field, backend, query_count);
    if (selected < WINDING_FIELD_DIRECT || selected > WINDING_FIELD_FAST)
        return -1;
    unsigned long long exact_triangles = 0;
    unsigned long long approximate_nodes = 0;
    unsigned long long degenerate_rays = 0;
    unsigned long long invalid_normals = 0;
    int64_t count = (int64_t)npoints;
    int64_t pi = 0;
#ifdef _OPENMP
    int thread_count = wf_thread_count();
#pragma omp parallel for schedule(dynamic, 64) \
    num_threads(thread_count) \
    reduction(+:exact_triangles, approximate_nodes, degenerate_rays, invalid_normals)
#endif
    for (pi = 0; pi < count; pi++) {
        size_t i = (size_t)pi;
        double n[3] = { normals[i*3], normals[i*3+1], normals[i*3+2] };
        double length = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        WfQueryStats plus_stats, minus_stats;
        memset(&plus_stats, 0, sizeof plus_stats);
        memset(&minus_stats, 0, sizeof minus_stats);
        if (!(length > 0.0) || !isfinite(length)) {
            out_mean[i] = NAN;
            out_jump[i] = NAN;
            invalid_normals++;
            continue;
        }
        double plus[3], minus[3];
        for (int k = 0; k < 3; k++) {
            double unit = n[k] / length;
            plus[k] = (double)points[i*3 + (size_t)k] + eps * unit;
            minus[k] = (double)points[i*3 + (size_t)k] - eps * unit;
        }
        double wp = wf_eval_one(field, plus, selected, &plus_stats);
        double wm = wf_eval_one(field, minus, selected, &minus_stats);
        out_mean[i] = isfinite(wp) && isfinite(wm)
                    ? (float)(0.5 * (wp + wm)) : NAN;
        out_jump[i] = isfinite(wp) && isfinite(wm)
                    ? (float)(wm - wp) : NAN;
        exact_triangles += (unsigned long long)(
            plus_stats.exact_triangles + minus_stats.exact_triangles);
        approximate_nodes += (unsigned long long)(
            plus_stats.approximate_nodes + minus_stats.approximate_nodes);
        degenerate_rays += (unsigned long long)(
            plus_stats.degenerate_ray + minus_stats.degenerate_ray);
    }
    if (stats != NULL) {
        memset(stats, 0, sizeof *stats);
        stats->backend_used = selected;
        stats->queries = query_count;
        stats->direct_queries = selected == WINDING_FIELD_DIRECT
                              ? query_count : 0;
        stats->boundary_queries = selected == WINDING_FIELD_BOUNDARY_EXACT
                                ? query_count : 0;
        stats->fast_queries = selected == WINDING_FIELD_FAST ? query_count : 0;
        stats->degenerate_ray_fallbacks = (size_t)degenerate_rays;
        stats->exact_triangle_evaluations = (size_t)exact_triangles;
        stats->approximate_node_evaluations = (size_t)approximate_nodes;
        stats->invalid_normals = (size_t)invalid_normals;
    }
    return 0;
}

static void wf_selftest_check(int condition, const char *message, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[winding field selftest] FAIL: %s\n", message);
    (*fails)++;
}

int WindingField_selftest(void)
{
    static const float cube_vertices[8*3] = {
        0,0,0, 0,0,1, 0,1,0, 0,1,1,
        1,0,0, 1,0,1, 1,1,0, 1,1,1
    };
    static const int32_t cube_faces[12*3] = {
        0,1,3, 0,3,2,  4,6,7, 4,7,5,
        0,4,5, 0,5,1,  2,3,7, 2,7,6,
        0,2,6, 0,6,4,  1,5,7, 1,7,3
    };
    static const double queries[4*3] = {
        0.5,0.5,0.5,  2.0,2.0,2.0,
        0.25,0.25,0.25, -1.0,0.5,0.5
    };
    Arena_T arena = Arena_new();
    WindingField_T closed = NULL, open = NULL;
    WindingFieldOptions options;
    WindingFieldStats build_stats;
    double direct[4], boundary[4], fast[4];
    int fails = 0;
    if (arena == NULL) return 1;
    WindingField_default_options(&options);
    options.leaf_faces = 2;
    wf_selftest_check(
        WindingField_build(arena, cube_vertices, 8, cube_faces, 12,
                           &options, &closed, &build_stats) == 0,
        "closed cube builds", &fails);
    if (closed != NULL) {
        wf_selftest_check(build_stats.exterior_edges == 0,
                          "closed cube has no exterior edges", &fails);
        wf_selftest_check(
            WindingField_evaluate_many(closed, queries, 4,
                WINDING_FIELD_DIRECT, direct, NULL) == 0 &&
            WindingField_evaluate_many(closed, queries, 4,
                WINDING_FIELD_BOUNDARY_EXACT, boundary, NULL) == 0 &&
            WindingField_evaluate_many(closed, queries, 4,
                WINDING_FIELD_FAST, fast, NULL) == 0,
            "all closed-cube backends evaluate", &fails);
        for (int i = 0; i < 4; i++) {
            wf_selftest_check(fabs(direct[i] - boundary[i]) < 1e-10,
                              "closed boundary-cone agrees with direct", &fails);
            wf_selftest_check(fabs(direct[i] - fast[i]) < 2e-3,
                              "closed fast field agrees with direct", &fails);
        }
        wf_selftest_check(fabs(fabs(direct[0]) - 1.0) < 1e-10 &&
                          fabs(direct[1]) < 1e-10,
                          "cube inside/outside winding", &fails);
    }

    Arena_Mark mark = Arena_save(arena);
    wf_selftest_check(
        WindingField_build(arena, cube_vertices, 8, cube_faces, 10,
                           &options, &open, &build_stats) == 0,
        "open cube builds", &fails);
    if (open != NULL) {
        wf_selftest_check(build_stats.exterior_edges == 4,
                          "open cube exposes four oriented edges", &fails);
        wf_selftest_check(
            WindingField_evaluate_many(open, queries, 4,
                WINDING_FIELD_DIRECT, direct, NULL) == 0 &&
            WindingField_evaluate_many(open, queries, 4,
                WINDING_FIELD_BOUNDARY_EXACT, boundary, NULL) == 0,
            "open direct/exact backends evaluate", &fails);
        for (int i = 0; i < 4; i++)
            wf_selftest_check(fabs(direct[i] - boundary[i]) < 1e-10,
                              "open boundary-cone agrees with direct", &fails);
    }
    Arena_restore(arena, mark);
    Arena_dispose(&arena);
    fprintf(stderr, "[winding field selftest] %s (%d failure(s))\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}
