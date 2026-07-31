#include "overlap_quality.h"

#include "../common/csr.h"
#include "../common/mesh_manifold.h"
#include "../common/obj_io.h"
#include "../common/pca.h"
#include "../remesh/cvt_remesh.h"
#include "../remesh/orient_mesh.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ANSI_DECLARATORS
#define REAL double
#define VOID int
#include "triangle.h"

extern jmp_buf triangle_jmpbuf;
extern int triangle_jmpbuf_set;

typedef struct {
    uint64_t key;
    int32_t  face;
} QualityEdge;

typedef struct {
    size_t faces;
    size_t degenerate;
    size_t slivers_10;
    double mean_min_angle;
    double worst_min_angle;
} QualityMeasure;

static uint64_t quality_edge_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static int compare_quality_edge(const void *pa, const void *pb)
{
    const QualityEdge *a = (const QualityEdge *)pa;
    const QualityEdge *b = (const QualityEdge *)pb;
    if (a->key < b->key) return -1;
    if (a->key > b->key) return 1;
    return 0;
}

static void free_quality_triangle_output(struct triangulateio *out)
{
    if (out->pointlist) trifree((VOID *)out->pointlist);
    if (out->pointattributelist) trifree((VOID *)out->pointattributelist);
    if (out->pointmarkerlist) trifree((VOID *)out->pointmarkerlist);
    if (out->trianglelist) trifree((VOID *)out->trianglelist);
    if (out->triangleattributelist)
        trifree((VOID *)out->triangleattributelist);
    if (out->neighborlist) trifree((VOID *)out->neighborlist);
    if (out->segmentlist) trifree((VOID *)out->segmentlist);
    if (out->segmentmarkerlist) trifree((VOID *)out->segmentmarkerlist);
    if (out->edgelist) trifree((VOID *)out->edgelist);
    if (out->edgemarkerlist) trifree((VOID *)out->edgemarkerlist);
    if (out->normlist) trifree((VOID *)out->normlist);
}

static int constrained_cvt_dual(Arena_T arena,
                                const double *uv, size_t nv,
                                const int32_t *segments, size_t n_segments,
                                int32_t **out_faces, size_t *out_nf)
{
    *out_faces = NULL;
    *out_nf = 0;
    if (nv < 3 || nv > (size_t)INT_MAX || n_segments < 3 ||
        n_segments > (size_t)INT_MAX)
        return -1;

    struct triangulateio input;
    struct triangulateio output;
    memset(&input, 0, sizeof(input));
    memset(&output, 0, sizeof(output));
    input.numberofpoints = (int)nv;
    input.pointlist = (REAL *)malloc(nv * 2 * sizeof(REAL));
    input.numberofsegments = (int)n_segments;
    input.segmentlist = (int *)malloc(n_segments * 2 * sizeof(int));
    if (!input.pointlist || !input.segmentlist) {
        free(input.pointlist);
        free(input.segmentlist);
        return -1;
    }
    memcpy(input.pointlist, uv, nv * 2 * sizeof(REAL));
    for (size_t i = 0; i < n_segments * 2; i++)
        input.segmentlist[i] = (int)segments[i];

    triangle_jmpbuf_set = 1;
    if (setjmp(triangle_jmpbuf) != 0) {
        triangle_jmpbuf_set = 0;
        free(input.pointlist);
        free(input.segmentlist);
        free_quality_triangle_output(&output);
        return -1;
    }
    /* The CVT already selected and relaxed the sites.  Triangle supplies only
     * the exact PSLG dual here: no point insertion, no boundary splitting. */
    triangulate("pYzQ", &input, &output, NULL);
    triangle_jmpbuf_set = 0;
    free(input.pointlist);
    free(input.segmentlist);

    if (output.numberofpoints != (int)nv || !output.pointlist ||
        output.numberoftriangles <= 0 || !output.trianglelist) {
        free_quality_triangle_output(&output);
        return -1;
    }
    for (size_t vi = 0; vi < nv; vi++) {
        if (fabs((double)output.pointlist[vi * 2] - uv[vi * 2]) > 1.0e-9 ||
            fabs((double)output.pointlist[vi * 2 + 1] -
                 uv[vi * 2 + 1]) > 1.0e-9) {
            free_quality_triangle_output(&output);
            return -1;
        }
    }
    size_t nf = (size_t)output.numberoftriangles;
    int32_t *faces = (int32_t *)ARENA_ALLOC(
        arena, (long)(nf * 3 * sizeof(*faces)));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t a = (int32_t)output.trianglelist[fi * 3];
        int32_t b = (int32_t)output.trianglelist[fi * 3 + 1];
        int32_t c = (int32_t)output.trianglelist[fi * 3 + 2];
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv ||
            (size_t)b >= nv || (size_t)c >= nv) {
            free_quality_triangle_output(&output);
            return -1;
        }
        double area2 = (uv[(size_t)b * 2] - uv[(size_t)a * 2]) *
                       (uv[(size_t)c * 2 + 1] - uv[(size_t)a * 2 + 1]) -
                       (uv[(size_t)b * 2 + 1] - uv[(size_t)a * 2 + 1]) *
                       (uv[(size_t)c * 2] - uv[(size_t)a * 2]);
        if (fabs(area2) <= 1.0e-14) {
            free_quality_triangle_output(&output);
            return -1;
        }
        faces[fi * 3] = a;
        faces[fi * 3 + 1] = area2 > 0.0 ? b : c;
        faces[fi * 3 + 2] = area2 > 0.0 ? c : b;
    }
    free_quality_triangle_output(&output);
    *out_faces = faces;
    *out_nf = nf;
    return 0;
}

static double clamp_unit(double x)
{
    if (x < -1.0) return -1.0;
    if (x > 1.0) return 1.0;
    return x;
}

static double face_min_angle(const float *verts, const int32_t tri[3],
                             double *out_area)
{
    double side[3];
    for (int edge = 0; edge < 3; edge++) {
        const float *a = &verts[(size_t)tri[(edge + 1) % 3] * 3];
        const float *b = &verts[(size_t)tri[(edge + 2) % 3] * 3];
        double d0 = (double)b[0] - a[0];
        double d1 = (double)b[1] - a[1];
        double d2 = (double)b[2] - a[2];
        side[edge] = sqrt(d0 * d0 + d1 * d1 + d2 * d2);
    }

    const float *a = &verts[(size_t)tri[0] * 3];
    const float *b = &verts[(size_t)tri[1] * 3];
    const float *c = &verts[(size_t)tri[2] * 3];
    double ab[3] = {
        (double)b[0] - a[0], (double)b[1] - a[1],
        (double)b[2] - a[2]
    };
    double ac[3] = {
        (double)c[0] - a[0], (double)c[1] - a[1],
        (double)c[2] - a[2]
    };
    double cross[3] = {
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0]
    };
    double area = 0.5 * sqrt(cross[0] * cross[0] +
                             cross[1] * cross[1] +
                             cross[2] * cross[2]);
    if (out_area) *out_area = area;
    if (side[0] <= 1.0e-12 || side[1] <= 1.0e-12 ||
        side[2] <= 1.0e-12 || area <= 1.0e-14)
        return 0.0;

    const double radians_to_degrees =
        57.295779513082320876798154814105;
    double minimum = 180.0;
    for (int corner = 0; corner < 3; corner++) {
        double opposite = side[corner];
        double first = side[(corner + 1) % 3];
        double second = side[(corner + 2) % 3];
        double cosine = (first * first + second * second -
                         opposite * opposite) / (2.0 * first * second);
        double angle = acos(clamp_unit(cosine)) * radians_to_degrees;
        if (angle < minimum) minimum = angle;
    }
    return minimum;
}

static QualityMeasure measure_quality(const float *verts,
                                      const int32_t *faces, size_t nf)
{
    QualityMeasure result;
    memset(&result, 0, sizeof(result));
    result.faces = nf;
    result.worst_min_angle = DBL_MAX;
    for (size_t fi = 0; fi < nf; fi++) {
        double area = 0.0;
        double angle = face_min_angle(verts, &faces[fi * 3], &area);
        if (area <= 1.0e-12) result.degenerate++;
        if (angle < 10.0) result.slivers_10++;
        if (angle < result.worst_min_angle) result.worst_min_angle = angle;
        result.mean_min_angle += angle;
    }
    if (nf > 0) result.mean_min_angle /= (double)nf;
    else result.worst_min_angle = 0.0;
    return result;
}

static void add_csr_weight(const int32_t *offset, const int32_t *target,
                           double *weight, int32_t a, int32_t b,
                           double value)
{
    for (int32_t edge = offset[a]; edge < offset[a + 1]; edge++) {
        if (target[edge] == b) {
            weight[edge] += value;
            return;
        }
    }
}

static double orient_chart(const float *uv, int32_t a, int32_t b, int32_t c)
{
    double ab0 = (double)uv[(size_t)b * 2] - uv[(size_t)a * 2];
    double ab1 = (double)uv[(size_t)b * 2 + 1] - uv[(size_t)a * 2 + 1];
    double ac0 = (double)uv[(size_t)c * 2] - uv[(size_t)a * 2];
    double ac1 = (double)uv[(size_t)c * 2 + 1] - uv[(size_t)a * 2 + 1];
    return ab0 * ac1 - ab1 * ac0;
}

static double cotangent_chart(const float *uv,
                              int32_t a, int32_t b, int32_t opposite)
{
    double ax = (double)uv[(size_t)a * 2] -
                uv[(size_t)opposite * 2];
    double ay = (double)uv[(size_t)a * 2 + 1] -
                uv[(size_t)opposite * 2 + 1];
    double bx = (double)uv[(size_t)b * 2] -
                uv[(size_t)opposite * 2];
    double by = (double)uv[(size_t)b * 2 + 1] -
                uv[(size_t)opposite * 2 + 1];
    double cross = fabs(ax * by - ay * bx);
    if (cross <= 1.0e-14) return 0.0;
    return (ax * bx + ay * by) / cross;
}

static int harmonic_height_patch(Arena_T arena,
                                 const float *source_verts,
                                 const int32_t *faces, size_t nf, size_t nv,
                                 const float *chart_uv,
                                 const uint8_t *boundary,
                                 const float center[3],
                                 const float normal_in[3],
                                 float *fair_verts)
{
    float normal[3] = {normal_in[0], normal_in[1], normal_in[2]};
    double length = sqrt((double)normal[0] * normal[0] +
                         (double)normal[1] * normal[1] +
                         (double)normal[2] * normal[2]);
    if (length <= 1.0e-12) return -1;
    for (int d = 0; d < 3; d++) normal[d] = (float)(normal[d] / length);

    float axis_u[3], axis_v[3];
    PCA_orthonormal_basis(normal, axis_u, axis_v);
    CSR_T adjacency = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *offset = CSR_offset(adjacency);
    const int32_t *target = CSR_target(adjacency);
    size_t directed_edges = (size_t)offset[nv];
    double *weight = (double *)ARENA_CALLOC(
        arena, (long)directed_edges, (long)sizeof(*weight));

    for (size_t fi = 0; fi < nf; fi++) {
        int32_t a = faces[fi * 3];
        int32_t b = faces[fi * 3 + 1];
        int32_t c = faces[fi * 3 + 2];
        if (fabs(orient_chart(chart_uv, a, b, c)) <= 1.0e-14)
            return -1;
        double cot_c = cotangent_chart(chart_uv, a, b, c);
        double cot_a = cotangent_chart(chart_uv, b, c, a);
        double cot_b = cotangent_chart(chart_uv, c, a, b);
        add_csr_weight(offset, target, weight, a, b, cot_c);
        add_csr_weight(offset, target, weight, b, a, cot_c);
        add_csr_weight(offset, target, weight, b, c, cot_a);
        add_csr_weight(offset, target, weight, c, b, cot_a);
        add_csr_weight(offset, target, weight, c, a, cot_b);
        add_csr_weight(offset, target, weight, a, c, cot_b);
    }

    float *height = (float *)ARENA_ALLOC(
        arena, (long)(nv * sizeof(*height)));
    float *next = (float *)ARENA_ALLOC(
        arena, (long)(nv * sizeof(*next)));
    double boundary_sum = 0.0;
    size_t boundary_count = 0;
    for (size_t vi = 0; vi < nv; vi++) {
        const float *p = &source_verts[vi * 3];
        height[vi] = (p[0] - center[0]) * normal[0] +
                     (p[1] - center[1]) * normal[1] +
                     (p[2] - center[2]) * normal[2];
        if (boundary[vi]) {
            boundary_sum += height[vi];
            boundary_count++;
        }
    }
    if (boundary_count < 3) return -1;
    float initial = (float)(boundary_sum / (double)boundary_count);
    for (size_t vi = 0; vi < nv; vi++)
        if (!boundary[vi]) height[vi] = initial;

    for (int iteration = 0; iteration < 3000; iteration++) {
        memcpy(next, height, nv * sizeof(*next));
        double max_delta = 0.0;
        for (size_t vi = 0; vi < nv; vi++) {
            if (boundary[vi]) continue;
            double sum = 0.0, sum_weight = 0.0;
            for (int32_t edge = offset[vi]; edge < offset[vi + 1]; edge++) {
                int32_t vj = target[edge];
                double w = weight[edge];
                if (w <= 1.0e-8) {
                    double du = (double)chart_uv[vi * 2] -
                                chart_uv[(size_t)vj * 2];
                    double dv = (double)chart_uv[vi * 2 + 1] -
                                chart_uv[(size_t)vj * 2 + 1];
                    w = 1.0 / fmax(1.0e-6, sqrt(du * du + dv * dv));
                }
                sum += w * height[vj];
                sum_weight += w;
            }
            if (sum_weight <= 0.0) continue;
            float value = (float)(sum / sum_weight);
            double delta = fabs((double)value - height[vi]);
            if (delta > max_delta) max_delta = delta;
            next[vi] = value;
        }
        float *swap = height;
        height = next;
        next = swap;
        if (max_delta < 1.0e-6) break;
    }

    for (size_t vi = 0; vi < nv; vi++) {
        if (boundary[vi]) {
            memcpy(&fair_verts[vi * 3], &source_verts[vi * 3],
                   3 * sizeof(float));
            continue;
        }
        double u = chart_uv[vi * 2];
        double v = chart_uv[vi * 2 + 1];
        for (int d = 0; d < 3; d++) {
            fair_verts[vi * 3 + (size_t)d] =
                center[d] + (float)(u * axis_u[d] + v * axis_v[d]) +
                height[vi] * normal[d];
        }
    }
    return 0;
}

static int32_t quality_find(int32_t *parent, int32_t x)
{
    int32_t root = x;
    while (parent[root] != root) root = parent[root];
    while (parent[x] != x) {
        int32_t next = parent[x];
        parent[x] = root;
        x = next;
    }
    return root;
}

static void quality_union(int32_t *parent, uint8_t *rank,
                          int32_t a, int32_t b)
{
    a = quality_find(parent, a);
    b = quality_find(parent, b);
    if (a == b) return;
    if (rank[a] < rank[b]) parent[a] = b;
    else if (rank[a] > rank[b]) parent[b] = a;
    else {
        parent[b] = a;
        rank[a]++;
    }
}

static size_t edge_group_count(const QualityEdge *edges, size_t n_edges,
                               uint64_t key)
{
    size_t lo = 0, hi = n_edges;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (edges[mid].key < key) lo = mid + 1;
        else hi = mid;
    }
    if (lo >= n_edges || edges[lo].key != key) return 0;
    size_t end = lo + 1;
    while (end < n_edges && edges[end].key == key) end++;
    return end - lo;
}

static int cvt_patch_conforming(Arena_T arena,
                                const int32_t *input_faces, size_t input_nf,
                                size_t input_nv,
                                const uint8_t *input_boundary,
                                const int32_t *cvt_faces, size_t cvt_nf,
                                size_t cvt_nv,
                                const int32_t *site_src, size_t n_pin,
                                size_t expected_boundary_edges)
{
    int32_t *pin_for_input = (int32_t *)ARENA_ALLOC(
        arena, (long)(input_nv * sizeof(*pin_for_input)));
    for (size_t vi = 0; vi < input_nv; vi++) pin_for_input[vi] = -1;
    for (size_t pin = 0; pin < n_pin; pin++) {
        int32_t src = site_src[pin];
        if (src < 0 || (size_t)src >= input_nv || !input_boundary[src] ||
            pin_for_input[src] >= 0)
            return 0;
        pin_for_input[src] = (int32_t)pin;
    }

    QualityEdge *input_edges = (QualityEdge *)ARENA_ALLOC(
        arena, (long)(input_nf * 3 * sizeof(*input_edges)));
    for (size_t fi = 0; fi < input_nf; fi++) {
        for (int edge = 0; edge < 3; edge++) {
            int32_t a = input_faces[fi * 3 + (size_t)edge];
            int32_t b = input_faces[fi * 3 + (size_t)((edge + 1) % 3)];
            input_edges[fi * 3 + (size_t)edge].key =
                quality_edge_key(a, b);
            input_edges[fi * 3 + (size_t)edge].face = (int32_t)fi;
        }
    }
    qsort(input_edges, input_nf * 3, sizeof(*input_edges),
          compare_quality_edge);

    QualityEdge *output_edges = (QualityEdge *)ARENA_ALLOC(
        arena, (long)(cvt_nf * 3 * sizeof(*output_edges)));
    for (size_t fi = 0; fi < cvt_nf; fi++) {
        for (int edge = 0; edge < 3; edge++) {
            int32_t a = cvt_faces[fi * 3 + (size_t)edge];
            int32_t b = cvt_faces[fi * 3 + (size_t)((edge + 1) % 3)];
            if (a < 0 || b < 0 || (size_t)a >= cvt_nv ||
                (size_t)b >= cvt_nv || a == b)
                return 0;
            output_edges[fi * 3 + (size_t)edge].key =
                quality_edge_key(a, b);
            output_edges[fi * 3 + (size_t)edge].face = (int32_t)fi;
        }
    }
    qsort(output_edges, cvt_nf * 3, sizeof(*output_edges),
          compare_quality_edge);

    size_t input_boundary_edges = 0;
    for (size_t begin = 0; begin < input_nf * 3;) {
        size_t end = begin + 1;
        while (end < input_nf * 3 &&
               input_edges[end].key == input_edges[begin].key) end++;
        if (end - begin == 1) {
            input_boundary_edges++;
            int32_t a = (int32_t)(input_edges[begin].key >> 32);
            int32_t b = (int32_t)(input_edges[begin].key &
                                  UINT32_C(0xffffffff));
            if (pin_for_input[a] < 0 || pin_for_input[b] < 0 ||
                edge_group_count(output_edges, cvt_nf * 3,
                                 quality_edge_key(pin_for_input[a],
                                                  pin_for_input[b])) != 1)
                return 0;
        }
        begin = end;
    }
    if (input_boundary_edges != expected_boundary_edges) return 0;

    size_t output_boundary_edges = 0;
    int32_t *parent = (int32_t *)ARENA_ALLOC(
        arena, (long)(cvt_nf * sizeof(*parent)));
    uint8_t *rank = (uint8_t *)ARENA_CALLOC(
        arena, (long)cvt_nf, (long)sizeof(*rank));
    for (size_t fi = 0; fi < cvt_nf; fi++) parent[fi] = (int32_t)fi;
    for (size_t begin = 0; begin < cvt_nf * 3;) {
        size_t end = begin + 1;
        while (end < cvt_nf * 3 &&
               output_edges[end].key == output_edges[begin].key) end++;
        size_t count = end - begin;
        if (count == 1) {
            int32_t a = (int32_t)(output_edges[begin].key >> 32);
            int32_t b = (int32_t)(output_edges[begin].key &
                                  UINT32_C(0xffffffff));
            if ((size_t)a >= n_pin || (size_t)b >= n_pin)
                return 0;
            output_boundary_edges++;
        } else if (count == 2) {
            quality_union(parent, rank, output_edges[begin].face,
                          output_edges[begin + 1].face);
        } else {
            return 0;
        }
        begin = end;
    }
    if (output_boundary_edges != expected_boundary_edges) return 0;

    size_t components = 0;
    for (size_t fi = 0; fi < cvt_nf; fi++)
        if (quality_find(parent, (int32_t)fi) == (int32_t)fi) components++;
    if (components != 1) return 0;

    MeshManifoldStats manifold =
        MeshManifold_audit(arena, cvt_nv, cvt_faces, cvt_nf);
    if (!MeshManifold_ok(&manifold) ||
        manifold.boundary_edges != expected_boundary_edges)
        return 0;
    if (cvt_nf != 2 * cvt_nv - expected_boundary_edges - 2)
        return 0;
    return 1;
}

static double env_quality_ratio(void)
{
    const char *text = getenv("VES_OVERLAP_MERGE_QUALITY_RATIO");
    if (!text || !*text) return 1.0;
    char *end = NULL;
    double value = strtod(text, &end);
    if (end == text || !isfinite(value)) return 1.0;
    if (value < 0.75) value = 0.75;
    if (value > 2.0) value = 2.0;
    return value;
}

int OverlapQuality_improve(Arena_T arena,
                           ComponentMesh *candidate,
                           OverlapMergeStats *merge,
                           const char *debug_dir,
                           OverlapQualityStats *stats)
{
    assert(arena && candidate && merge && stats);
    memset(stats, 0, sizeof(*stats));
    stats->reason = "invalid-input";
    if (!candidate->verts || !candidate->faces || !merge->vertex_role ||
        !merge->chart_uv || merge->patch_face_start >= candidate->nf)
        return -1;

    const size_t old_nv = candidate->nv;
    const size_t old_nf = candidate->nf;
    const size_t patch_start = merge->patch_face_start;
    const size_t patch_nf = old_nf - patch_start;
    int32_t *old_to_local = (int32_t *)ARENA_ALLOC(
        arena, (long)(old_nv * sizeof(*old_to_local)));
    for (size_t vi = 0; vi < old_nv; vi++) old_to_local[vi] = -1;
    size_t patch_nv = 0;
    for (size_t fi = patch_start; fi < old_nf; fi++) {
        for (int k = 0; k < 3; k++) {
            int32_t vi = candidate->faces[fi * 3 + (size_t)k];
            if (vi < 0 || (size_t)vi >= old_nv) return -1;
            if (old_to_local[vi] < 0)
                old_to_local[vi] = (int32_t)patch_nv++;
        }
    }
    if (patch_nv < 4 || patch_nf < 2) return -1;

    int32_t *local_to_old = (int32_t *)ARENA_ALLOC(
        arena, (long)(patch_nv * sizeof(*local_to_old)));
    for (size_t old = 0; old < old_nv; old++)
        if (old_to_local[old] >= 0)
            local_to_old[(size_t)old_to_local[old]] = (int32_t)old;

    float *patch_verts = (float *)ARENA_ALLOC(
        arena, (long)(patch_nv * 3 * sizeof(*patch_verts)));
    float *patch_chart = (float *)ARENA_ALLOC(
        arena, (long)(patch_nv * 2 * sizeof(*patch_chart)));
    for (size_t local = 0; local < patch_nv; local++) {
        size_t old = (size_t)local_to_old[local];
        memcpy(&patch_verts[local * 3], &candidate->verts[old * 3],
               3 * sizeof(float));
        memcpy(&patch_chart[local * 2], &merge->chart_uv[old * 2],
               2 * sizeof(float));
    }
    int32_t *patch_faces = (int32_t *)ARENA_ALLOC(
        arena, (long)(patch_nf * 3 * sizeof(*patch_faces)));
    for (size_t local_fi = 0; local_fi < patch_nf; local_fi++) {
        size_t fi = patch_start + local_fi;
        for (int k = 0; k < 3; k++)
            patch_faces[local_fi * 3 + (size_t)k] =
                old_to_local[candidate->faces[fi * 3 + (size_t)k]];
    }

    QualityEdge *patch_edges = (QualityEdge *)ARENA_ALLOC(
        arena, (long)(patch_nf * 3 * sizeof(*patch_edges)));
    for (size_t fi = 0; fi < patch_nf; fi++) {
        for (int edge = 0; edge < 3; edge++) {
            int32_t a = patch_faces[fi * 3 + (size_t)edge];
            int32_t b = patch_faces[fi * 3 + (size_t)((edge + 1) % 3)];
            patch_edges[fi * 3 + (size_t)edge].key =
                quality_edge_key(a, b);
            patch_edges[fi * 3 + (size_t)edge].face = (int32_t)fi;
        }
    }
    qsort(patch_edges, patch_nf * 3, sizeof(*patch_edges),
          compare_quality_edge);
    uint8_t *boundary = (uint8_t *)ARENA_CALLOC(
        arena, (long)patch_nv, (long)sizeof(*boundary));
    size_t boundary_edges = 0;
    for (size_t begin = 0; begin < patch_nf * 3;) {
        size_t end = begin + 1;
        while (end < patch_nf * 3 &&
               patch_edges[end].key == patch_edges[begin].key) end++;
        if (end - begin == 1) {
            int32_t a = (int32_t)(patch_edges[begin].key >> 32);
            int32_t b = (int32_t)(patch_edges[begin].key &
                                  UINT32_C(0xffffffff));
            boundary[a] = boundary[b] = 1;
            boundary_edges++;
        } else if (end - begin != 2) {
            stats->reason = "input-patch-nonmanifold";
            return -1;
        }
        begin = end;
    }
    size_t boundary_vertices = 0;
    for (size_t local = 0; local < patch_nv; local++) {
        if (!boundary[local]) continue;
        boundary_vertices++;
    }
    if (boundary_vertices != boundary_edges || boundary_vertices < 3 ||
        boundary_vertices >= patch_nv) {
        stats->reason = "input-patch-not-one-disk";
        return -1;
    }

    QualityMeasure input_quality =
        measure_quality(patch_verts, patch_faces, patch_nf);
    float *fair_verts = (float *)ARENA_ALLOC(
        arena, (long)(patch_nv * 3 * sizeof(*fair_verts)));
    if (harmonic_height_patch(
            arena, patch_verts, patch_faces, patch_nf, patch_nv,
            patch_chart, boundary, merge->plane_center,
            merge->plane_normal, fair_verts) != 0) {
        stats->reason = "harmonic-fairing-failed";
        return -1;
    }
    QualityMeasure fair_quality =
        measure_quality(fair_verts, patch_faces, patch_nf);
    if (fair_quality.degenerate != 0) {
        stats->reason = "harmonic-fairing-degenerate";
        return -1;
    }

    double area = 0.0;
    for (size_t fi = 0; fi < patch_nf; fi++) {
        double face_area = 0.0;
        (void)face_min_angle(fair_verts, &patch_faces[fi * 3], &face_area);
        area += face_area;
    }
    size_t old_interior = patch_nv - boundary_vertices;
    double ratio = env_quality_ratio();
    size_t desired_interior =
        (size_t)floor((double)old_interior * ratio + 0.5);
    if (desired_interior < 1) desired_interior = 1;
    double target_h = sqrt(2.0 * area / (double)desired_interior);
    if (!isfinite(target_h) || target_h <= 1.0e-6) {
        stats->reason = "invalid-quality-spacing";
        return -1;
    }

    CvtOpts options;
    CVT_default_opts(&options);
    options.n_iters = 10;
    options.seed = UINT32_C(0x4d455247);
    options.verbose = 0;
    float *cvt_verts = NULL;
    int32_t *rvd_faces = NULL;
    int32_t *site_src = NULL;
    size_t cvt_nv = 0, rvd_nf = 0, n_pin = 0;
    int cvt_rc = CVT_remesh_pinned(
        arena, fair_verts, patch_nv, patch_faces, patch_nf, target_h,
        &options, &cvt_verts, &cvt_nv, &rvd_faces, &rvd_nf,
        &site_src, &n_pin);
    if (cvt_rc != 0 || !cvt_verts || !rvd_faces || !site_src ||
        n_pin != boundary_vertices || cvt_nv <= n_pin || rvd_nf == 0) {
        stats->reason = "pinned-cvt-failed";
        return -1;
    }

    /*
     * RVD is the right relaxation operator, but its restricted dual is not a
     * conforming patch boundary guarantee: a nearby interior site can steal a
     * pinned-pinned boundary adjacency.  Keep the relaxed sites, then construct
     * their exact planar constrained-Delaunay dual using every original
     * frontier edge as a PSLG segment.
     */
    int32_t *pin_for_input = (int32_t *)ARENA_ALLOC(
        arena, (long)(patch_nv * sizeof(*pin_for_input)));
    for (size_t local = 0; local < patch_nv; local++)
        pin_for_input[local] = -1;
    for (size_t pin = 0; pin < n_pin; pin++) {
        int32_t src = site_src[pin];
        if (src < 0 || (size_t)src >= patch_nv || !boundary[src] ||
            pin_for_input[src] >= 0) {
            stats->reason = "pinned-cvt-site-map-failed";
            return -1;
        }
        pin_for_input[src] = (int32_t)pin;
    }

    int32_t *cvt_segments = (int32_t *)ARENA_ALLOC(
        arena, (long)(boundary_edges * 2 * sizeof(*cvt_segments)));
    size_t segment_count = 0;
    for (size_t begin = 0; begin < patch_nf * 3;) {
        size_t end = begin + 1;
        while (end < patch_nf * 3 &&
               patch_edges[end].key == patch_edges[begin].key) end++;
        if (end - begin == 1) {
            int32_t a = (int32_t)(patch_edges[begin].key >> 32);
            int32_t b = (int32_t)(patch_edges[begin].key &
                                  UINT32_C(0xffffffff));
            if (pin_for_input[a] < 0 || pin_for_input[b] < 0 ||
                segment_count >= boundary_edges) {
                stats->reason = "pinned-cvt-segment-map-failed";
                return -1;
            }
            cvt_segments[segment_count * 2] = pin_for_input[a];
            cvt_segments[segment_count * 2 + 1] = pin_for_input[b];
            segment_count++;
        }
        begin = end;
    }
    if (segment_count != boundary_edges) {
        stats->reason = "pinned-cvt-segment-count-failed";
        return -1;
    }

    float normal[3] = {
        merge->plane_normal[0], merge->plane_normal[1],
        merge->plane_normal[2]
    };
    float axis_u[3], axis_v[3];
    PCA_orthonormal_basis(normal, axis_u, axis_v);
    double *cvt_uv = (double *)ARENA_ALLOC(
        arena, (long)(cvt_nv * 2 * sizeof(*cvt_uv)));
    for (size_t site = 0; site < n_pin; site++) {
        size_t src = (size_t)site_src[site];
        cvt_uv[site * 2] = patch_chart[src * 2];
        cvt_uv[site * 2 + 1] = patch_chart[src * 2 + 1];
    }
    for (size_t site = n_pin; site < cvt_nv; site++) {
        double d0 = (double)cvt_verts[site * 3] - merge->plane_center[0];
        double d1 = (double)cvt_verts[site * 3 + 1] -
                    merge->plane_center[1];
        double d2 = (double)cvt_verts[site * 3 + 2] -
                    merge->plane_center[2];
        cvt_uv[site * 2] =
            d0 * axis_u[0] + d1 * axis_u[1] + d2 * axis_u[2];
        cvt_uv[site * 2 + 1] =
            d0 * axis_v[0] + d1 * axis_v[1] + d2 * axis_v[2];
    }

    int32_t *cvt_faces = NULL;
    size_t cvt_nf = 0;
    if (constrained_cvt_dual(
            arena, cvt_uv, cvt_nv, cvt_segments, boundary_edges,
            &cvt_faces, &cvt_nf) != 0 ||
        !cvt_faces || cvt_nf == 0) {
        stats->reason = "pinned-cdt-failed";
        return -1;
    }
    if (!cvt_patch_conforming(
            arena, patch_faces, patch_nf, patch_nv, boundary,
            cvt_faces, cvt_nf, cvt_nv, site_src, n_pin, boundary_edges)) {
        stats->reason = "pinned-cdt-frontier-mismatch";
        return -1;
    }

    QualityMeasure output_quality =
        measure_quality(cvt_verts, cvt_faces, cvt_nf);
    fprintf(stderr,
            "    delamination mesh quality candidate: "
            "V/F %zu/%zu -> %zu/%zu, boundary=%zu, "
            "min-angle mean %.2f -> %.2f harmonic -> %.2f CVT, "
            "worst %.2f -> %.2f -> %.2f, <10deg %zu -> %zu -> %zu\n",
            patch_nv, patch_nf, cvt_nv, cvt_nf, boundary_vertices,
            input_quality.mean_min_angle, fair_quality.mean_min_angle,
            output_quality.mean_min_angle, input_quality.worst_min_angle,
            fair_quality.worst_min_angle, output_quality.worst_min_angle,
            input_quality.slivers_10, fair_quality.slivers_10,
            output_quality.slivers_10);

    if (debug_dir && *debug_dir) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/03a_harmonic_pre_cvt_patch.obj",
                 debug_dir);
        ObjIO_write(path, fair_verts, patch_nv, patch_faces, patch_nf);
        snprintf(path, sizeof(path), "%s/03b_pinned_cvt_patch.obj",
                 debug_dir);
        ObjIO_write(path, cvt_verts, cvt_nv, cvt_faces, cvt_nf);
    }

    if (output_quality.degenerate != 0 ||
        output_quality.mean_min_angle + 0.5 <
            fair_quality.mean_min_angle ||
        output_quality.slivers_10 > fair_quality.slivers_10 ||
        output_quality.worst_min_angle + 0.25 <
            fair_quality.worst_min_angle) {
        /*
         * The planar heal is already a constrained Delaunay triangulation.
         * Its harmonic height solve can therefore be the better quality
         * result: RVD relaxation is optional, not grounds to discard a valid
         * smooth heal.  Select the harmonic candidate only when it materially
         * improves the actual lifted input patch.
         */
        if (fair_quality.degenerate != 0 ||
            fair_quality.slivers_10 > input_quality.slivers_10 ||
            fair_quality.worst_min_angle + 0.25 <
                input_quality.worst_min_angle) {
            stats->reason = "harmonic-and-cvt-quality-regression";
            return -1;
        }
        for (size_t local = 0; local < patch_nv; local++) {
            size_t old = (size_t)local_to_old[local];
            memcpy(&candidate->verts[old * 3], &fair_verts[local * 3],
                   3 * sizeof(float));
        }
        PCA_normal(candidate->verts, candidate->nv,
                   candidate->pca_normal, candidate->centroid);
        candidate->self = candidate;

        stats->reason = "applied-harmonic";
        stats->applied = 1;
        stats->input_vertices = patch_nv;
        stats->input_faces = patch_nf;
        stats->output_vertices = patch_nv;
        stats->output_faces = patch_nf;
        stats->boundary_vertices = boundary_vertices;
        stats->target_h = target_h;
        stats->input_slivers_10 = input_quality.slivers_10;
        stats->fair_slivers_10 = fair_quality.slivers_10;
        stats->output_slivers_10 = fair_quality.slivers_10;
        stats->input_mean_min_angle = input_quality.mean_min_angle;
        stats->fair_mean_min_angle = fair_quality.mean_min_angle;
        stats->output_mean_min_angle = fair_quality.mean_min_angle;
        stats->input_worst_min_angle = input_quality.worst_min_angle;
        stats->fair_worst_min_angle = fair_quality.worst_min_angle;
        stats->output_worst_min_angle = fair_quality.worst_min_angle;
        fprintf(stderr,
                "    delamination mesh quality: selected harmonic patch "
                "(CVT regressed worst/slivers); V/F %zu/%zu, "
                "mean min-angle %.2f -> %.2f, worst %.2f -> %.2f, "
                "<10deg %zu -> %zu\n",
                patch_nv, patch_nf, input_quality.mean_min_angle,
                fair_quality.mean_min_angle, input_quality.worst_min_angle,
                fair_quality.worst_min_angle, input_quality.slivers_10,
                fair_quality.slivers_10);
        return 0;
    }

    uint8_t *old_kept = (uint8_t *)ARENA_CALLOC(
        arena, (long)old_nv, (long)sizeof(*old_kept));
    for (size_t fi = 0; fi < patch_start; fi++)
        for (int k = 0; k < 3; k++)
            old_kept[candidate->faces[fi * 3 + (size_t)k]] = 1;
    /* The disk frontier consists of both the retained collar interface and
     * any genuine open boundary of the source sheet.  The latter has no
     * exterior incident face, but it is still an exact positional/topological
     * constraint and must survive as a pinned vertex. */
    for (size_t local = 0; local < patch_nv; local++) {
        if (!boundary[local]) continue;
        old_kept[(size_t)local_to_old[local]] = 1;
    }

    int32_t *old_remap = (int32_t *)ARENA_ALLOC(
        arena, (long)(old_nv * sizeof(*old_remap)));
    size_t kept_nv = 0;
    for (size_t old = 0; old < old_nv; old++) {
        if (old_kept[old]) old_remap[old] = (int32_t)kept_nv++;
        else old_remap[old] = -1;
    }
    size_t new_interior = cvt_nv - n_pin;
    size_t new_nv = kept_nv + new_interior;
    size_t new_nf = patch_start + cvt_nf;
    float *new_verts = (float *)ARENA_ALLOC(
        arena, (long)(new_nv * 3 * sizeof(*new_verts)));
    uint8_t *new_role = (uint8_t *)ARENA_CALLOC(
        arena, (long)new_nv, (long)sizeof(*new_role));
    float *new_chart = (float *)ARENA_ALLOC(
        arena, (long)(new_nv * 2 * sizeof(*new_chart)));
    for (size_t old = 0; old < old_nv; old++) {
        if (old_remap[old] < 0) continue;
        size_t next = (size_t)old_remap[old];
        memcpy(&new_verts[next * 3], &candidate->verts[old * 3],
               3 * sizeof(float));
        memcpy(&new_chart[next * 2], &merge->chart_uv[old * 2],
               2 * sizeof(float));
        new_role[next] = merge->vertex_role[old];
    }

    for (size_t site = n_pin; site < cvt_nv; site++) {
        size_t next = kept_nv + site - n_pin;
        memcpy(&new_verts[next * 3], &cvt_verts[site * 3],
               3 * sizeof(float));
        double d0 = (double)cvt_verts[site * 3] - merge->plane_center[0];
        double d1 = (double)cvt_verts[site * 3 + 1] - merge->plane_center[1];
        double d2 = (double)cvt_verts[site * 3 + 2] - merge->plane_center[2];
        new_chart[next * 2] = (float)(
            d0 * axis_u[0] + d1 * axis_u[1] + d2 * axis_u[2]);
        new_chart[next * 2 + 1] = (float)(
            d0 * axis_v[0] + d1 * axis_v[1] + d2 * axis_v[2]);
        new_role[next] = OVERLAP_MERGE_VERTEX_REPAIR;
    }

    int32_t *new_faces = (int32_t *)ARENA_ALLOC(
        arena, (long)(new_nf * 3 * sizeof(*new_faces)));
    for (size_t fi = 0; fi < patch_start; fi++) {
        for (int k = 0; k < 3; k++) {
            int32_t old = candidate->faces[fi * 3 + (size_t)k];
            new_faces[fi * 3 + (size_t)k] = old_remap[old];
            if (old_remap[old] < 0) {
                stats->reason = "retained-remap-failed";
                return -1;
            }
        }
    }
    int32_t *site_global = (int32_t *)ARENA_ALLOC(
        arena, (long)(cvt_nv * sizeof(*site_global)));
    for (size_t site = 0; site < n_pin; site++) {
        int32_t local = site_src[site];
        int32_t old = local_to_old[local];
        site_global[site] = old_remap[old];
        if (site_global[site] < 0) {
            stats->reason = "pinned-remap-failed";
            return -1;
        }
        new_role[site_global[site]] = OVERLAP_MERGE_VERTEX_ANCHOR;
    }
    for (size_t site = n_pin; site < cvt_nv; site++)
        site_global[site] = (int32_t)(kept_nv + site - n_pin);
    for (size_t fi = 0; fi < cvt_nf; fi++)
        for (int k = 0; k < 3; k++)
            new_faces[(patch_start + fi) * 3 + (size_t)k] =
                site_global[cvt_faces[fi * 3 + (size_t)k]];

    size_t flips = 0, components = 0, residual = 0;
    OrientMesh_consistent(arena, new_verts, new_nv, NULL,
                          new_faces, new_nf, &flips, &components, &residual);
    MeshManifoldStats before = MeshManifold_audit(
        arena, candidate->nv, candidate->faces, candidate->nf);
    MeshManifoldStats after =
        MeshManifold_audit(arena, new_nv, new_faces, new_nf);
    if (residual != 0 || components != 1 || !MeshManifold_ok(&after) ||
        after.same_dir_edges != 0 ||
        after.boundary_edges != before.boundary_edges) {
        stats->reason = "improved-candidate-topology-failed";
        return -1;
    }

    ComponentMesh improved = *candidate;
    improved.verts = new_verts;
    improved.faces = new_faces;
    improved.nv = new_nv;
    improved.nf = new_nf;
    improved.pin_mask = NULL;
    improved.vert_normals = NULL;
    PCA_normal(improved.verts, improved.nv,
               improved.pca_normal, improved.centroid);
    improved.self = &improved;

    *candidate = improved;
    candidate->self = candidate;
    merge->vertex_role = new_role;
    merge->chart_uv = new_chart;
    merge->patch_faces_added = cvt_nf;

    stats->reason = "applied";
    stats->applied = 1;
    stats->input_vertices = patch_nv;
    stats->input_faces = patch_nf;
    stats->output_vertices = cvt_nv;
    stats->output_faces = cvt_nf;
    stats->boundary_vertices = boundary_vertices;
    stats->target_h = target_h;
    stats->input_slivers_10 = input_quality.slivers_10;
    stats->fair_slivers_10 = fair_quality.slivers_10;
    stats->output_slivers_10 = output_quality.slivers_10;
    stats->input_mean_min_angle = input_quality.mean_min_angle;
    stats->fair_mean_min_angle = fair_quality.mean_min_angle;
    stats->output_mean_min_angle = output_quality.mean_min_angle;
    stats->input_worst_min_angle = input_quality.worst_min_angle;
    stats->fair_worst_min_angle = fair_quality.worst_min_angle;
    stats->output_worst_min_angle = output_quality.worst_min_angle;
    fprintf(stderr,
            "    delamination mesh quality: pinned CVT h=%.3f "
            "V/F %zu/%zu -> %zu/%zu, boundary=%zu, "
            "min-angle mean %.2f -> %.2f harmonic -> %.2f CVT, "
            "worst %.2f -> %.2f -> %.2f, <10deg %zu -> %zu -> %zu\n",
            target_h, patch_nv, patch_nf, cvt_nv, cvt_nf,
            boundary_vertices, input_quality.mean_min_angle,
            fair_quality.mean_min_angle, output_quality.mean_min_angle,
            input_quality.worst_min_angle, fair_quality.worst_min_angle,
            output_quality.worst_min_angle, input_quality.slivers_10,
            fair_quality.slivers_10, output_quality.slivers_10);
    (void)flips;
    return 0;
}
