/* slim_refine.c -- metric refinement after an injective disk initialization.
 *
 * This is a C local/global SLIM-style solver for symmetric Dirichlet energy.
 * The local step builds a weighted rotation proxy per triangle.  The global
 * step solves the coupled 2-D block system with PCG.  A whole-map line search
 * accepts only lower-energy candidates with positive triangle Jacobians and
 * simple component boundaries. */
#include "slim_refine.h"
#include "../common/uv_guard.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include "../common/ves_omp.h"
#endif

#define SR_AREA_EPS 1.0e-12

typedef struct SRFace {
    double gx[3];
    double gy[3];
    double area;
} SRFace;

typedef struct SRBoundaryEdge {
    int32_t a;
    int32_t b;
    int32_t component;
} SRBoundaryEdge;

typedef struct SREdgeBox {
    int32_t a;
    int32_t b;
    int32_t component;
    double min_u;
    double max_u;
    double min_v;
    double max_v;
} SREdgeBox;

typedef struct SRSystem {
    size_t nv;
    size_t nf;
    size_t ne;
    size_t nb;
    size_t ncomponents;
    size_t metric_uv_faces;
    size_t metric_fallback_faces;
    const int32_t *faces;
    SRFace *face;
    int32_t *edge_a;
    int32_t *edge_b;
    int32_t *face_edge;
    size_t *offset;
    int32_t *neighbor;
    size_t *adj_edge;
    int32_t *vertex_component;
    size_t *component_faces;
    size_t *component_vertex_offset;
    int32_t *component_vertex;
    size_t *component_face_offset;
    size_t *component_face;
    size_t *component_boundary_offset;
    size_t *component_boundary;
    SRBoundaryEdge *boundary;
    double *diagonal;              /* [nv*3]: 00,01,11 */
    double *edge_block;            /* [ne*3]: 00,01,11 */
    double *rhs;                   /* [nv*2] */
    double proximal_min;
    double proximal_mean;
    double proximal_max;
    int threads;
} SRSystem;

#define SR_AMG_MAX_LEVELS 24
#define SR_AMG_COARSE_VERTICES 64
#define SR_AMG_AGGREGATE_SIZE 4
#define SR_AMG_SMOOTH_SWEEPS 2
#define SR_AMG_OMEGA 0.5

typedef struct SRAMGLevel {
    size_t n;
    size_t ne;
    int32_t *edge_a;
    int32_t *edge_b;
    size_t *offset;
    int32_t *neighbor;
    size_t *adj_edge;
    double *diagonal;
    double *edge_block;
    int32_t *aggregate;
    int32_t *edge_to_coarse;
    size_t ncoarse;
    double *cholesky;
} SRAMGLevel;

typedef struct SRAMG {
    int levels;
    int threads;
    SRAMGLevel level[SR_AMG_MAX_LEVELS];
} SRAMG;

typedef struct SRAMGWork {
    double *x[SR_AMG_MAX_LEVELS];
    double *r[SR_AMG_MAX_LEVELS];
    double *temporary[SR_AMG_MAX_LEVELS];
} SRAMGWork;

typedef struct SRCoarseRef {
    uint64_t key;
    size_t fine_edge;
} SRCoarseRef;

typedef struct SRMeasure {
    double energy;
    double qc_mean;
    double qc_max;
    double min_det;
    size_t flips;
} SRMeasure;

typedef struct SRPack {
    int32_t component;
    size_t faces;
    double min_u;
    double max_u;
    double min_v;
    double max_v;
} SRPack;

static int sr_threads(int requested)
{
#ifdef _OPENMP
    return requested > 0 ? requested : omp_get_max_threads();
#else
    (void)requested;
    return 1;
#endif
}

static int sr_compare_u64(const void *left, const void *right)
{
    uint64_t a = *(const uint64_t *)left;
    uint64_t b = *(const uint64_t *)right;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static uint64_t sr_edge_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static size_t sr_find_key(const uint64_t *key, size_t count, uint64_t target)
{
    size_t lo = 0;
    size_t hi = count;
    while (lo < hi) {
        size_t middle = lo + (hi - lo) / 2;
        if (key[middle] < target) lo = middle + 1;
        else hi = middle;
    }
    return lo < count && key[lo] == target ? lo : SIZE_MAX;
}

static int sr_face_rest(const float *p0, const float *p1, const float *p2,
                        SRFace *face)
{
    double e1[3] = { 0.0, 0.0, 0.0 };
    double e2[3] = { 0.0, 0.0, 0.0 };
    double unit[3] = { 0.0, 0.0, 0.0 };
    double length = 0.0;
    double projection = 0.0;
    double height = 0.0;
    size_t axis = 0;
    for (axis = 0; axis < 3; axis++) {
        e1[axis] = (double)p1[axis] - (double)p0[axis];
        e2[axis] = (double)p2[axis] - (double)p0[axis];
        length += e1[axis] * e1[axis];
    }
    length = sqrt(length);
    if (!(length > 1.0e-12) || !isfinite(length)) return -1;
    for (axis = 0; axis < 3; axis++) {
        unit[axis] = e1[axis] / length;
        projection += e2[axis] * unit[axis];
    }
    /* Squared-length subtraction loses the height of a thin source face. */
    double cross[3] = {e1[1]*e2[2]-e1[2]*e2[1],
                       e1[2]*e2[0]-e1[0]*e2[2],
                       e1[0]*e2[1]-e1[1]*e2[0]};
    height = hypot(hypot(cross[0],cross[1]),cross[2])/length;
    if (!(height > 1.0e-12) || !isfinite(height)) return -1;

    face->gx[0] = -1.0 / length;
    face->gx[1] = 1.0 / length;
    face->gx[2] = 0.0;
    face->gy[0] = projection / (length * height) - 1.0 / height;
    face->gy[1] = -projection / (length * height);
    face->gy[2] = 1.0 / height;
    face->area = 0.5 * length * height;
    return 0;
}

/* Build the same rotation-invariant triangle rest frame from a trusted 2-D
 * winding/developable parameterization.  Only its local metric is transferred:
 * global translations, rotations, reflections, and overlaps in metric_uv have
 * no effect on the objective. */
static int sr_face_rest_uv(const float *p0, const float *p1, const float *p2,
                           SRFace *face)
{
    double e10 = (double)p1[0] - (double)p0[0];
    double e11 = (double)p1[1] - (double)p0[1];
    double e20 = (double)p2[0] - (double)p0[0];
    double e21 = (double)p2[1] - (double)p0[1];
    double length = hypot(e10, e11);
    double projection = 0.0;
    double height = 0.0;

    if (!(length > 1.0e-12) || !isfinite(length)) return -1;
    projection = (e20 * e10 + e21 * e11) / length;
    height = fabs(e10*e21-e11*e20)/length;
    if (!(height > 1.0e-12) || !isfinite(height)) return -1;

    face->gx[0] = -1.0 / length;
    face->gx[1] = 1.0 / length;
    face->gx[2] = 0.0;
    face->gy[0] = projection / (length * height) - 1.0 / height;
    face->gy[1] = -projection / (length * height);
    face->gy[2] = 1.0 / height;
    face->area = 0.5 * length * height;
    return 0;
}

static int sr_face_rest_uv_double(const double *p0, const double *p1,
                                  const double *p2, SRFace *face)
{
    double e10 = p1[0] - p0[0];
    double e11 = p1[1] - p0[1];
    double e20 = p2[0] - p0[0];
    double e21 = p2[1] - p0[1];
    double length = hypot(e10, e11);
    double projection = 0.0;
    double height = 0.0;
    if (!(length > 1.0e-12) || !isfinite(length)) return -1;
    projection = (e20 * e10 + e21 * e11) / length;
    /* Cross-product height avoids cancellation for long, thin, valid
     * scaffold triangles. Squared-length subtraction can round to zero. */
    height = fabs(e10 * e21 - e11 * e20) / length;
    if (!(height > 1.0e-12) || !isfinite(height)) return -1;
    face->gx[0] = -1.0 / length;
    face->gx[1] = 1.0 / length;
    face->gx[2] = 0.0;
    face->gy[0] = projection / (length * height) - 1.0 / height;
    face->gy[1] = -projection / (length * height);
    face->gy[2] = 1.0 / height;
    face->area = 0.5 * length * height;
    return 0;
}

/* Interpolate the upper-triangular, rotation-invariant rest frame rather than
 * its inverse gradients.  Positive edge length and height stay positive for
 * the entire path, so every continuation stage remains a valid metric. */
static int sr_blend_rest_faces(SRSystem *system, const SRFace *initial,
                               const SRFace *target, double t)
{
    if (system == NULL || initial == NULL || target == NULL ||
        !(t >= 0.0 && t <= 1.0))
        return -1;
    for (size_t f = 0; f < system->nf; f++) {
        double l0 = 1.0 / initial[f].gx[1];
        double h0 = 1.0 / initial[f].gy[2];
        double p0 = -initial[f].gy[1] * l0 * h0;
        double l1 = 1.0 / target[f].gx[1];
        double h1 = 1.0 / target[f].gy[2];
        double p1 = -target[f].gy[1] * l1 * h1;
        double length = (1.0 - t) * l0 + t * l1;
        double height = (1.0 - t) * h0 + t * h1;
        double projection = (1.0 - t) * p0 + t * p1;
        SRFace *face = &system->face[f];
        if (!(length > 1.0e-12) || !(height > 1.0e-12) ||
            !isfinite(length) || !isfinite(height) ||
            !isfinite(projection))
            return -1;
        face->gx[0] = -1.0 / length;
        face->gx[1] = 1.0 / length;
        face->gx[2] = 0.0;
        face->gy[0] = projection / (length * height) - 1.0 / height;
        face->gy[1] = -projection / (length * height);
        face->gy[2] = 1.0 / height;
        face->area = 0.5 * length * height;
    }
    return 0;
}

static int sr_build_system(Arena_T arena, const float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const SlimRefineOpts *opts, int threads,
                           SRSystem *system)
{
    uint64_t *keys = NULL;
    uint8_t *is_boundary = NULL;
    size_t key_count = 0;
    size_t unique_count = 0;
    size_t boundary_count = 0;
    size_t face_index = 0;
    size_t edge_index = 0;
    size_t vertex = 0;
    size_t *cursor = NULL;
    int32_t *queue = NULL;
    size_t queue_head = 0;
    size_t queue_tail = 0;
    int32_t component = 0;

    if (arena == NULL || verts == NULL || faces == NULL || system == NULL ||
        nv == 0 || nf == 0 || nv > (size_t)INT32_MAX ||
        nf > SIZE_MAX / 3 || nf * 3 > SIZE_MAX / sizeof(*keys)) return -1;
    memset(system, 0, sizeof(*system));
    system->nv = nv;
    system->nf = nf;
    system->faces = faces;
    system->threads = sr_threads(threads);
    key_count = nf * 3;
    keys = (uint64_t *)ARENA_ALLOC(arena, key_count * sizeof(*keys));
    is_boundary = (uint8_t *)ARENA_CALLOC(arena, key_count, sizeof(*is_boundary));

    for (face_index = 0; face_index < nf; face_index++) {
        size_t corner = 0;
        for (corner = 0; corner < 3; corner++) {
            int32_t a = faces[face_index * 3 + corner];
            int32_t b = faces[face_index * 3 + (corner + 1) % 3];
            if (a < 0 || b < 0 || (size_t)a >= nv || (size_t)b >= nv || a == b)
                return -1;
            keys[face_index * 3 + corner] = sr_edge_key(a, b);
        }
    }
    qsort(keys, key_count, sizeof(*keys), sr_compare_u64);
    for (edge_index = 0; edge_index < key_count;) {
        size_t end = edge_index + 1;
        while (end < key_count && keys[end] == keys[edge_index]) end++;
        if (end - edge_index > 2) {
            fprintf(stderr, "  [slim] non-manifold edge multiplicity %zu\n",
                    end - edge_index);
            return -1;
        }
        keys[unique_count] = keys[edge_index];
        is_boundary[unique_count] = (uint8_t)(end - edge_index == 1);
        if (is_boundary[unique_count]) boundary_count++;
        unique_count++;
        edge_index = end;
    }
    system->ne = unique_count;
    system->nb = boundary_count;
    system->edge_a = (int32_t *)ARENA_ALLOC(
        arena, unique_count * sizeof(*system->edge_a));
    system->edge_b = (int32_t *)ARENA_ALLOC(
        arena, unique_count * sizeof(*system->edge_b));
    for (edge_index = 0; edge_index < unique_count; edge_index++) {
        system->edge_a[edge_index] = (int32_t)(keys[edge_index] >> 32);
        system->edge_b[edge_index] =
            (int32_t)(keys[edge_index] & UINT64_C(0xffffffff));
    }

    system->face_edge = (int32_t *)ARENA_ALLOC(
        arena, key_count * sizeof(*system->face_edge));
    system->face = (SRFace *)ARENA_ALLOC(arena, nf * sizeof(*system->face));
    for (face_index = 0; face_index < nf; face_index++) {
        size_t corner = 0;
        int32_t a = faces[face_index * 3];
        int32_t b = faces[face_index * 3 + 1];
        int32_t c = faces[face_index * 3 + 2];
        int rest_status = -1;
        const float *metric_uv = opts->metric_uv;
        if (opts->scaffold_uv != NULL && face_index >= opts->scaffold_first_face) {
            rest_status = sr_face_rest_uv_double(
                &opts->scaffold_uv[(size_t)a * 2],
                &opts->scaffold_uv[(size_t)b * 2],
                &opts->scaffold_uv[(size_t)c * 2], &system->face[face_index]);
            if (rest_status == 0)
                system->face[face_index].area = opts->scaffold_weight;
        } else if (metric_uv != NULL) {
            rest_status = sr_face_rest_uv(&metric_uv[(size_t)a * 2],
                                          &metric_uv[(size_t)b * 2],
                                          &metric_uv[(size_t)c * 2],
                                          &system->face[face_index]);
            if (rest_status == 0) {
                system->metric_uv_faces++;
            } else {
                /* A winding field loaded from float OBJ can quantize an
                 * isolated sub-pixel triangle to a line.  Fall back atomically
                 * for that face; never discard the good guide metric globally. */
                rest_status = sr_face_rest(
                    &verts[(size_t)a * 3], &verts[(size_t)b * 3],
                    &verts[(size_t)c * 3], &system->face[face_index]);
                if (rest_status == 0) system->metric_fallback_faces++;
            }
        } else {
            rest_status = sr_face_rest(
                &verts[(size_t)a * 3], &verts[(size_t)b * 3],
                &verts[(size_t)c * 3], &system->face[face_index]);
        }
        if (rest_status != 0) {
            fprintf(stderr, "  [slim] degenerate rest face %zu\n", face_index);
            return -1;
        }
        for (corner = 0; corner < 3; corner++) {
            int32_t p = faces[face_index * 3 + corner];
            int32_t q = faces[face_index * 3 + (corner + 1) % 3];
            size_t found = sr_find_key(keys, unique_count, sr_edge_key(p, q));
            if (found == SIZE_MAX || found > (size_t)INT32_MAX) return -1;
            system->face_edge[face_index * 3 + corner] = (int32_t)found;
        }
    }

    system->offset = (size_t *)ARENA_CALLOC(
        arena, nv + 1, sizeof(*system->offset));
    for (edge_index = 0; edge_index < unique_count; edge_index++) {
        system->offset[(size_t)system->edge_a[edge_index] + 1]++;
        system->offset[(size_t)system->edge_b[edge_index] + 1]++;
    }
    for (vertex = 0; vertex < nv; vertex++)
        system->offset[vertex + 1] += system->offset[vertex];
    system->neighbor = (int32_t *)ARENA_ALLOC(
        arena, system->offset[nv] * sizeof(*system->neighbor));
    system->adj_edge = (size_t *)ARENA_ALLOC(
        arena, system->offset[nv] * sizeof(*system->adj_edge));
    cursor = (size_t *)ARENA_ALLOC(arena, nv * sizeof(*cursor));
    memcpy(cursor, system->offset, nv * sizeof(*cursor));
    for (edge_index = 0; edge_index < unique_count; edge_index++) {
        int32_t a = system->edge_a[edge_index];
        int32_t b = system->edge_b[edge_index];
        size_t pa = cursor[(size_t)a]++;
        size_t pb = cursor[(size_t)b]++;
        system->neighbor[pa] = b;
        system->adj_edge[pa] = edge_index;
        system->neighbor[pb] = a;
        system->adj_edge[pb] = edge_index;
    }

    system->vertex_component = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*system->vertex_component));
    queue = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(*queue));
    for (vertex = 0; vertex < nv; vertex++) system->vertex_component[vertex] = -1;
    for (vertex = 0; vertex < nv; vertex++) {
        if (system->vertex_component[vertex] >= 0) continue;
        queue_head = 0;
        queue_tail = 0;
        queue[queue_tail++] = (int32_t)vertex;
        system->vertex_component[vertex] = component;
        while (queue_head < queue_tail) {
            int32_t current = queue[queue_head++];
            size_t position = 0;
            for (position = system->offset[(size_t)current];
                 position < system->offset[(size_t)current + 1]; position++) {
                int32_t next = system->neighbor[position];
                if (system->vertex_component[(size_t)next] < 0) {
                    system->vertex_component[(size_t)next] = component;
                    queue[queue_tail++] = next;
                }
            }
        }
        component++;
    }
    system->ncomponents = (size_t)component;
    system->component_faces = (size_t *)ARENA_CALLOC(
        arena, system->ncomponents, sizeof(*system->component_faces));
    for (face_index = 0; face_index < nf; face_index++) {
        int32_t a = faces[face_index * 3];
        int32_t b = faces[face_index * 3 + 1];
        int32_t c = faces[face_index * 3 + 2];
        int32_t label = system->vertex_component[(size_t)a];
        if (label < 0 || system->vertex_component[(size_t)b] != label ||
            system->vertex_component[(size_t)c] != label) return -1;
        system->component_faces[(size_t)label]++;
    }

    system->boundary = (SRBoundaryEdge *)ARENA_ALLOC(
        arena, boundary_count * sizeof(*system->boundary));
    boundary_count = 0;
    for (edge_index = 0; edge_index < unique_count; edge_index++) {
        if (is_boundary[edge_index]) {
            int32_t a = system->edge_a[edge_index];
            int32_t b = system->edge_b[edge_index];
            int32_t label = system->vertex_component[(size_t)a];
            if (label < 0 || system->vertex_component[(size_t)b] != label)
                return -1;
            system->boundary[boundary_count].a = a;
            system->boundary[boundary_count].b = b;
            system->boundary[boundary_count].component = label;
            boundary_count++;
        }
    }
    if (boundary_count != system->nb) return -1;

    system->component_vertex_offset = (size_t *)ARENA_CALLOC(
        arena, system->ncomponents + 1,
        sizeof(*system->component_vertex_offset));
    system->component_face_offset = (size_t *)ARENA_CALLOC(
        arena, system->ncomponents + 1,
        sizeof(*system->component_face_offset));
    system->component_boundary_offset = (size_t *)ARENA_CALLOC(
        arena, system->ncomponents + 1,
        sizeof(*system->component_boundary_offset));
    system->component_vertex = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*system->component_vertex));
    system->component_face = (size_t *)ARENA_ALLOC(
        arena, nf * sizeof(*system->component_face));
    system->component_boundary = (size_t *)ARENA_ALLOC(
        arena, system->nb * sizeof(*system->component_boundary));
    for (vertex = 0; vertex < nv; vertex++)
        system->component_vertex_offset[
            (size_t)system->vertex_component[vertex] + 1]++;
    for (face_index = 0; face_index < nf; face_index++)
        system->component_face_offset[
            (size_t)system->vertex_component[
                (size_t)faces[face_index * 3]] + 1]++;
    for (edge_index = 0; edge_index < system->nb; edge_index++)
        system->component_boundary_offset[
            (size_t)system->boundary[edge_index].component + 1]++;
    for (edge_index = 0; edge_index < system->ncomponents; edge_index++) {
        system->component_vertex_offset[edge_index + 1] +=
            system->component_vertex_offset[edge_index];
        system->component_face_offset[edge_index + 1] +=
            system->component_face_offset[edge_index];
        system->component_boundary_offset[edge_index + 1] +=
            system->component_boundary_offset[edge_index];
    }
    memcpy(cursor, system->component_vertex_offset,
           system->ncomponents * sizeof(*cursor));
    for (vertex = 0; vertex < nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        system->component_vertex[cursor[label]++] = (int32_t)vertex;
    }
    memcpy(cursor, system->component_face_offset,
           system->ncomponents * sizeof(*cursor));
    for (face_index = 0; face_index < nf; face_index++) {
        size_t label = (size_t)system->vertex_component[
            (size_t)faces[face_index * 3]];
        system->component_face[cursor[label]++] = face_index;
    }
    memcpy(cursor, system->component_boundary_offset,
           system->ncomponents * sizeof(*cursor));
    for (edge_index = 0; edge_index < system->nb; edge_index++) {
        size_t label = (size_t)system->boundary[edge_index].component;
        system->component_boundary[cursor[label]++] = edge_index;
    }

    system->diagonal = (double *)ARENA_ALLOC(
        arena, nv * 3 * sizeof(*system->diagonal));
    system->edge_block = (double *)ARENA_ALLOC(
        arena, unique_count * 3 * sizeof(*system->edge_block));
    system->rhs = (double *)ARENA_ALLOC(arena, nv * 2 * sizeof(*system->rhs));
    return 0;
}

static double sr_orient2(const double *a, const double *b, const double *c)
{
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

static int sr_orientation(const double *a, const double *b, const double *c)
{
    double value = sr_orient2(a, b, c);
    double scale = fabs((b[0] - a[0]) * (c[1] - a[1])) +
                   fabs((b[1] - a[1]) * (c[0] - a[0])) + 1.0;
    double epsilon = 1.0e-12 * scale;
    return value > epsilon ? 1 : (value < -epsilon ? -1 : 0);
}

static int sr_on_segment(const double *a, const double *b, const double *p)
{
    double epsilon = 1.0e-10 * (hypot(b[0] - a[0], b[1] - a[1]) + 1.0);
    return sr_orientation(a, b, p) == 0 &&
           p[0] >= fmin(a[0], b[0]) - epsilon &&
           p[0] <= fmax(a[0], b[0]) + epsilon &&
           p[1] >= fmin(a[1], b[1]) - epsilon &&
           p[1] <= fmax(a[1], b[1]) + epsilon;
}

static int sr_segments_intersect(const double *a, const double *b,
                                 const double *c, const double *d)
{
    int o0 = sr_orientation(a, b, c);
    int o1 = sr_orientation(a, b, d);
    int o2 = sr_orientation(c, d, a);
    int o3 = sr_orientation(c, d, b);
    if (o0 * o1 < 0 && o2 * o3 < 0) return 1;
    if (o0 == 0 && sr_on_segment(a, b, c)) return 1;
    if (o1 == 0 && sr_on_segment(a, b, d)) return 1;
    if (o2 == 0 && sr_on_segment(c, d, a)) return 1;
    if (o3 == 0 && sr_on_segment(c, d, b)) return 1;
    return 0;
}

static int sr_compare_box(const void *left, const void *right)
{
    const SREdgeBox *a = (const SREdgeBox *)left;
    const SREdgeBox *b = (const SREdgeBox *)right;
    if (a->component != b->component)
        return a->component < b->component ? -1 : 1;
    if (a->min_u != b->min_u) return a->min_u < b->min_u ? -1 : 1;
    if (a->max_u != b->max_u) return a->max_u < b->max_u ? -1 : 1;
    if (a->min_v != b->min_v) return a->min_v < b->min_v ? -1 : 1;
    return a->max_v < b->max_v ? -1 : (a->max_v > b->max_v ? 1 : 0);
}

static size_t sr_boundary_intersections(const SRSystem *system,
                                        const double *uv, SREdgeBox *box)
{
    size_t edge = 0;
    size_t intersections = 0;
    for (edge = 0; edge < system->nb; edge++) {
        int32_t a = system->boundary[edge].a;
        int32_t b = system->boundary[edge].b;
        const double *pa = &uv[(size_t)a * 2];
        const double *pb = &uv[(size_t)b * 2];
        box[edge].a = a;
        box[edge].b = b;
        box[edge].component = system->boundary[edge].component;
        box[edge].min_u = fmin(pa[0], pb[0]);
        box[edge].max_u = fmax(pa[0], pb[0]);
        box[edge].min_v = fmin(pa[1], pb[1]);
        box[edge].max_v = fmax(pa[1], pb[1]);
    }
    qsort(box, system->nb, sizeof(*box), sr_compare_box);
    for (edge = 0; edge < system->nb; edge++) {
        size_t other = edge + 1;
        const SREdgeBox *a = &box[edge];
        double epsilon = 1.0e-10 *
            (hypot(a->max_u - a->min_u, a->max_v - a->min_v) + 1.0);
        for (; other < system->nb; other++) {
            const SREdgeBox *b = &box[other];
            if (b->component != a->component || b->min_u > a->max_u + epsilon)
                break;
            if (b->max_v < a->min_v - epsilon ||
                b->min_v > a->max_v + epsilon) continue;
            if (a->a == b->a || a->a == b->b ||
                a->b == b->a || a->b == b->b) continue;
            if (sr_segments_intersect(&uv[(size_t)a->a * 2],
                                      &uv[(size_t)a->b * 2],
                                      &uv[(size_t)b->a * 2],
                                      &uv[(size_t)b->b * 2])) {
                intersections++;
                if (intersections >= 16) return intersections;
            }
        }
    }
    return intersections;
}

static size_t sr_component_boundary_intersections(const SRSystem *system,
                                                  size_t component,
                                                  const double *uv,
                                                  SREdgeBox *box)
{
    size_t begin = system->component_boundary_offset[component];
    size_t end = system->component_boundary_offset[component + 1];
    size_t count = end - begin;
    size_t local = 0;
    size_t intersections = 0;
    for (local = 0; local < count; local++) {
        size_t edge = system->component_boundary[begin + local];
        int32_t a = system->boundary[edge].a;
        int32_t b = system->boundary[edge].b;
        const double *pa = &uv[(size_t)a * 2];
        const double *pb = &uv[(size_t)b * 2];
        box[local].a = a;
        box[local].b = b;
        box[local].component = (int32_t)component;
        box[local].min_u = fmin(pa[0], pb[0]);
        box[local].max_u = fmax(pa[0], pb[0]);
        box[local].min_v = fmin(pa[1], pb[1]);
        box[local].max_v = fmax(pa[1], pb[1]);
    }
    qsort(box, count, sizeof(*box), sr_compare_box);
    for (local = 0; local < count; local++) {
        size_t other = local + 1;
        const SREdgeBox *a = &box[local];
        double epsilon = 1.0e-10 *
            (hypot(a->max_u - a->min_u, a->max_v - a->min_v) + 1.0);
        for (; other < count; other++) {
            const SREdgeBox *b = &box[other];
            if (b->min_u > a->max_u + epsilon) break;
            if (b->max_v < a->min_v - epsilon ||
                b->min_v > a->max_v + epsilon) continue;
            if (a->a == b->a || a->a == b->b ||
                a->b == b->a || a->b == b->b) continue;
            if (sr_segments_intersect(&uv[(size_t)a->a * 2],
                                      &uv[(size_t)a->b * 2],
                                      &uv[(size_t)b->a * 2],
                                      &uv[(size_t)b->b * 2])) {
                intersections++;
                if (intersections >= 16) return intersections;
            }
        }
    }
    return intersections;
}

static void sr_face_jacobian(const SRFace *face, const int32_t *indices,
                             const double *uv, double jacobian[4])
{
    size_t corner = 0;
    jacobian[0] = 0.0;
    jacobian[1] = 0.0;
    jacobian[2] = 0.0;
    jacobian[3] = 0.0;
    for (corner = 0; corner < 3; corner++) {
        size_t vertex = (size_t)indices[corner];
        jacobian[0] += face->gx[corner] * uv[vertex * 2];
        jacobian[1] += face->gy[corner] * uv[vertex * 2];
        jacobian[2] += face->gx[corner] * uv[vertex * 2 + 1];
        jacobian[3] += face->gy[corner] * uv[vertex * 2 + 1];
    }
}

static int sr_measure(const SRSystem *system, const double *uv,
                      SRMeasure *measure)
{
    double energy = 0.0;
    double qc_sum = 0.0;
    double area_sum = 0.0;
    double qc_max = 1.0;
    double min_det = DBL_MAX;
    size_t flips = 0;
    size_t face_index = 0;
    for (face_index = 0; face_index < system->nf; face_index++) {
        double jacobian[4] = { 0.0, 0.0, 0.0, 0.0 };
        double det = 0.0;
        double frobenius = 0.0;
        double p = 0.0;
        double q = 0.0;
        double r = 0.0;
        double disc = 0.0;
        double lambda_max = 0.0;
        double qc = 0.0;
        const SRFace *face = &system->face[face_index];
        sr_face_jacobian(face, &system->faces[face_index * 3], uv, jacobian);
        det = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
        if (!(det > SR_AREA_EPS) || !isfinite(det)) flips++;
        if (det < min_det) min_det = det;
        frobenius = jacobian[0] * jacobian[0] +
                    jacobian[1] * jacobian[1] +
                    jacobian[2] * jacobian[2] +
                    jacobian[3] * jacobian[3];
        if (!(det > 0.0)) return -1;
        energy += face->area * frobenius * (1.0 + 1.0 / (det * det));
        p = jacobian[0] * jacobian[0] + jacobian[1] * jacobian[1];
        q = jacobian[0] * jacobian[2] + jacobian[1] * jacobian[3];
        r = jacobian[2] * jacobian[2] + jacobian[3] * jacobian[3];
        disc = hypot(p - r, 2.0 * q);
        lambda_max = fmax(0.0, 0.5 * (p + r + disc));
        /* For a 2x2 Jacobian, sigma_max/sigma_min =
         * lambda_max / |det(J)|.  Recovering lambda_min by subtracting two
         * nearly equal values catastrophically loses precision on exactly the
         * skinny triangles this diagnostic is meant to rank. */
        qc = det > DBL_MIN && isfinite(lambda_max)
           ? lambda_max / det : 1.0e15;
        if (!(qc >= 1.0) || !isfinite(qc)) qc = 1.0e15;
        if (qc > 1.0e15) qc = 1.0e15;
        qc_sum += face->area * qc;
        area_sum += face->area;
        if (qc > qc_max) qc_max = qc;
    }
    if (!isfinite(energy) || !isfinite(qc_sum) || !(area_sum > 0.0)) return -1;
    measure->energy = energy;
    measure->qc_mean = qc_sum / area_sum;
    measure->qc_max = qc_max;
    measure->min_det = min_det;
    measure->flips = flips;
    return 0;
}

static int sr_measure_component(const SRSystem *system, const double *uv,
                                size_t component, SRMeasure *measure)
{
    size_t begin = system->component_face_offset[component];
    size_t end = system->component_face_offset[component + 1];
    double energy = 0.0;
    double qc_sum = 0.0;
    double area_sum = 0.0;
    double qc_max = 1.0;
    double min_det = DBL_MAX;
    size_t position = 0;
    size_t flips = 0;
    if (begin == end) return -1;
    for (position = begin; position < end; position++) {
        size_t face_index = system->component_face[position];
        double jacobian[4] = { 0.0, 0.0, 0.0, 0.0 };
        double det = 0.0;
        double frobenius = 0.0;
        double p = 0.0;
        double q = 0.0;
        double r = 0.0;
        double disc = 0.0;
        double lambda_max = 0.0;
        double qc = 0.0;
        const SRFace *face = &system->face[face_index];
        sr_face_jacobian(face, &system->faces[face_index * 3], uv, jacobian);
        det = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
        if (!(det > SR_AREA_EPS) || !isfinite(det)) flips++;
        if (det < min_det) min_det = det;
        if (!(det > 0.0)) return -1;
        frobenius = jacobian[0] * jacobian[0] +
                    jacobian[1] * jacobian[1] +
                    jacobian[2] * jacobian[2] +
                    jacobian[3] * jacobian[3];
        energy += face->area * frobenius * (1.0 + 1.0 / (det * det));
        p = jacobian[0] * jacobian[0] + jacobian[1] * jacobian[1];
        q = jacobian[0] * jacobian[2] + jacobian[1] * jacobian[3];
        r = jacobian[2] * jacobian[2] + jacobian[3] * jacobian[3];
        disc = hypot(p - r, 2.0 * q);
        lambda_max = fmax(0.0, 0.5 * (p + r + disc));
        qc = det > DBL_MIN && isfinite(lambda_max)
           ? lambda_max / det : 1.0e15;
        if (!(qc >= 1.0) || !isfinite(qc)) qc = 1.0e15;
        if (qc > 1.0e15) qc = 1.0e15;
        qc_sum += face->area * qc;
        area_sum += face->area;
        if (qc > qc_max) qc_max = qc;
    }
    if (!isfinite(energy) || !isfinite(qc_sum) || !(area_sum > 0.0)) return -1;
    measure->energy = energy;
    measure->qc_mean = qc_sum / area_sum;
    measure->qc_max = qc_max;
    measure->min_det = min_det;
    measure->flips = flips;
    return 0;
}

static double sr_smallest_positive_root(double a, double b, double c)
{
    double scale = fabs(a) + fabs(b) + fabs(c) + 1.0;
    double root = DBL_MAX;
    if (fabs(a) <= 1.0e-14 * scale) {
        if (b < 0.0) {
            double linear = -c / b;
            if (linear > 0.0 && isfinite(linear)) root = linear;
        }
    } else {
        double discriminant = b * b - 4.0 * a * c;
        if (discriminant >= 0.0 && isfinite(discriminant)) {
            double square_root = sqrt(discriminant);
            double q = -0.5 * (b + copysign(square_root, b));
            double root0 = q / a;
            double root1 = fabs(q) > DBL_MIN ? c / q : DBL_MAX;
            if (root0 > 0.0 && isfinite(root0)) root = root0;
            if (root1 > 0.0 && isfinite(root1) && root1 < root) root = root1;
        }
    }
    return root;
}

/* Largest component-wise interpolation step before any triangle reaches zero
 * signed area.  This is the analytic feasibility bound used by the original
 * SLIM flip-avoiding line search; starting there avoids dozens of blind
 * halvings on severely distorted initial maps. */
static double sr_face_flip_step(const SRSystem *system,
                                const double *current,
                                const double *destination,
                                size_t face_index)
{
    int32_t i0 = system->faces[face_index * 3];
    int32_t i1 = system->faces[face_index * 3 + 1];
    int32_t i2 = system->faces[face_index * 3 + 2];
    double e10[2];
    double e20[2];
    double de1[2];
    double de2[2];
    double delta0[2];
    double delta1[2];
    double delta2[2];
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double root = 0.0;
    for (int coordinate = 0; coordinate < 2; coordinate++) {
        delta0[coordinate] = destination[(size_t)i0 * 2 + coordinate] -
                             current[(size_t)i0 * 2 + coordinate];
        delta1[coordinate] = destination[(size_t)i1 * 2 + coordinate] -
                             current[(size_t)i1 * 2 + coordinate];
        delta2[coordinate] = destination[(size_t)i2 * 2 + coordinate] -
                             current[(size_t)i2 * 2 + coordinate];
        e10[coordinate] = current[(size_t)i1 * 2 + coordinate] -
                          current[(size_t)i0 * 2 + coordinate];
        e20[coordinate] = current[(size_t)i2 * 2 + coordinate] -
                          current[(size_t)i0 * 2 + coordinate];
        de1[coordinate] = delta1[coordinate] - delta0[coordinate];
        de2[coordinate] = delta2[coordinate] - delta0[coordinate];
    }
    c = e10[0] * e20[1] - e10[1] * e20[0];
    b = de1[0] * e20[1] - de1[1] * e20[0] +
        e10[0] * de2[1] - e10[1] * de2[0];
    a = de1[0] * de2[1] - de1[1] * de2[0];
    if (!(c > 0.0) || !isfinite(c)) return 0.0;
    root = sr_smallest_positive_root(a, b, c);
    if (root == DBL_MAX || root >= 1.25) return 1.0;
    return fmax(0.0, 0.8 * root);
}

static double sr_component_flip_step(const SRSystem *system,
                                     const double *current,
                                     const double *destination,
                                     size_t component)
{
    size_t begin = system->component_face_offset[component];
    size_t end = system->component_face_offset[component + 1];
    double alpha = 1.0;
    size_t position = 0;
    for (position = begin; position < end; position++) {
        size_t face_index = system->component_face[position];
        double face_alpha = sr_face_flip_step(
            system, current, destination, face_index);
        if (face_alpha < alpha) alpha = face_alpha;
    }
    return alpha;
}

/* A single nearly collapsed triangle must not throttle an entire very large
 * disk.  Start from the global proxy destination, damp only vertices incident
 * on flip-limiting faces, and iteratively repair any nonuniform-step failures.
 * The result is still accepted only after the exact component energy and
 * simple-boundary gates pass. */
static int sr_component_local_step(const SRSystem *system,
                                   const double *current,
                                   const double *destination,
                                   size_t component,
                                   const SlimRefineOpts *opts,
                                   double old_energy,
                                   double *trial,
                                   double *vertex_alpha,
                                   uint8_t *damp_vertex,
                                   SREdgeBox *boxes,
                                   SRMeasure *candidate,
                                   size_t *out_trials,
                                   double *out_alpha_min,
                                   double *out_alpha_mean,
                                   double *out_alpha_max)
{
    size_t vertex_begin = system->component_vertex_offset[component];
    size_t vertex_end = system->component_vertex_offset[component + 1];
    size_t face_begin = system->component_face_offset[component];
    size_t face_end = system->component_face_offset[component + 1];
    size_t boundary_begin = system->component_boundary_offset[component];
    size_t boundary_end = system->component_boundary_offset[component + 1];
    size_t trials = 0;

    for (size_t position = vertex_begin; position < vertex_end; position++) {
        size_t vertex = (size_t)system->component_vertex[position];
        vertex_alpha[vertex] = 1.0;
    }
    for (size_t position = face_begin; position < face_end; position++) {
        size_t face_index = system->component_face[position];
        double alpha = sr_face_flip_step(
            system, current, destination, face_index);
        for (size_t corner = 0; corner < 3; corner++) {
            size_t vertex =
                (size_t)system->faces[face_index * 3 + corner];
            if (alpha < vertex_alpha[vertex])
                vertex_alpha[vertex] = alpha;
        }
    }

    for (size_t line = 0; line < (size_t)opts->line_search; line++) {
        size_t invalid_faces = 0;
        size_t boundary_intersections = 0;
        for (size_t position = vertex_begin; position < vertex_end;
             position++) {
            size_t vertex = (size_t)system->component_vertex[position];
            double alpha = vertex_alpha[vertex];
            damp_vertex[vertex] = 0;
            trial[vertex * 2] = current[vertex * 2] + alpha *
                (destination[vertex * 2] - current[vertex * 2]);
            trial[vertex * 2 + 1] = current[vertex * 2 + 1] + alpha *
                (destination[vertex * 2 + 1] - current[vertex * 2 + 1]);
        }
        trials++;
        for (size_t position = face_begin; position < face_end; position++) {
            size_t face_index = system->component_face[position];
            const int32_t *indices = &system->faces[face_index * 3];
            double jacobian[4];
            double determinant;
            sr_face_jacobian(&system->face[face_index], indices,
                             trial, jacobian);
            determinant = jacobian[0] * jacobian[3] -
                          jacobian[1] * jacobian[2];
            /* Local damping changes the triangle's path after the original
             * uniform flip bound. A positive endpoint can lie beyond two
             * roots; certify this actual nonuniform proposal's whole path. */
            if (!(determinant > SR_AREA_EPS) || !isfinite(determinant) ||
                !UvGuard_interval(current,trial,indices,1.0)) {
                invalid_faces++;
                for (size_t corner = 0; corner < 3; corner++)
                    damp_vertex[(size_t)indices[corner]] = 1;
            }
        }
        if (invalid_faces != 0) {
            if (opts->verbose &&
                (line < 3 || line + 1 == (size_t)opts->line_search))
                fprintf(stderr,
                        "      [slim local] c%zu trial=%zu invalid-faces=%zu\n",
                        component, line + 1, invalid_faces);
            for (size_t position = vertex_begin; position < vertex_end;
                 position++) {
                size_t vertex = (size_t)system->component_vertex[position];
                if (damp_vertex[vertex]) vertex_alpha[vertex] *= 0.5;
            }
            continue;
        }
        if (opts->guard_boundary)
            boundary_intersections = sr_component_boundary_intersections(
                system, component, trial, boxes);
        if (boundary_intersections != 0) {
            if (opts->verbose &&
                (line < 3 || line + 1 == (size_t)opts->line_search))
                fprintf(stderr,
                        "      [slim local] c%zu trial=%zu boundary-x=%zu\n",
                        component, line + 1, boundary_intersections);
            for (size_t position = boundary_begin; position < boundary_end;
                 position++) {
                size_t edge = system->component_boundary[position];
                damp_vertex[(size_t)system->boundary[edge].a] = 1;
                damp_vertex[(size_t)system->boundary[edge].b] = 1;
            }
            for (size_t position = vertex_begin; position < vertex_end;
                 position++) {
                size_t vertex = (size_t)system->component_vertex[position];
                if (damp_vertex[vertex]) vertex_alpha[vertex] *= 0.5;
            }
            continue;
        }
        if (sr_measure_component(system, trial, component, candidate) == 0 &&
            candidate->flips == 0 && candidate->energy < old_energy) {
            double alpha_min = DBL_MAX;
            double alpha_max = 0.0;
            double alpha_sum = 0.0;
            for (size_t position = vertex_begin; position < vertex_end;
                 position++) {
                size_t vertex = (size_t)system->component_vertex[position];
                double alpha = vertex_alpha[vertex];
                if (alpha < alpha_min) alpha_min = alpha;
                if (alpha > alpha_max) alpha_max = alpha;
                alpha_sum += alpha;
            }
            if (out_trials != NULL) *out_trials = trials;
            if (out_alpha_min != NULL) *out_alpha_min = alpha_min;
            if (out_alpha_mean != NULL)
                *out_alpha_mean =
                    alpha_sum / (double)(vertex_end - vertex_begin);
            if (out_alpha_max != NULL) *out_alpha_max = alpha_max;
            return 1;
        }
        if (opts->verbose &&
            (line < 3 || line + 1 == (size_t)opts->line_search))
            fprintf(stderr,
                    "      [slim local] c%zu trial=%zu non-descent "
                    "E=%.9g candidate=%.9g\n",
                    component, line + 1, old_energy,
                    candidate->energy);
        for (size_t position = vertex_begin; position < vertex_end;
             position++) {
            size_t vertex = (size_t)system->component_vertex[position];
            vertex_alpha[vertex] *= 0.5;
        }
    }
    if (out_trials != NULL) *out_trials = trials;
    return 0;
}

/* Express each disconnected winding/developable guide chart in the similarity
 * frame of the safe initializer.  The guide supplies shape, not an absolute
 * gauge.  A globally reversed chart is reflected once; a chart with mixed face
 * orientation is marked unusable and remains at the safe map. */
static int sr_align_guide(Arena_T arena, const SRSystem *system,
                          const double *safe, const float *guide,
                          double *aligned, uint8_t *usable)
{
    size_t ncomponents = system->ncomponents;
    double *count = NULL;
    double *guide_u = NULL;
    double *guide_v = NULL;
    double *safe_u = NULL;
    double *safe_v = NULL;
    double *coefficient_a = NULL;
    double *coefficient_b = NULL;
    double *norm = NULL;
    size_t *positive = NULL;
    size_t *negative = NULL;
    uint8_t *reflect = NULL;
    size_t vertex = 0;
    size_t face_index = 0;
    size_t component = 0;

    if (arena == NULL || system == NULL || safe == NULL || guide == NULL ||
        aligned == NULL || usable == NULL) return -1;
    count = (double *)ARENA_CALLOC(arena, ncomponents, sizeof(*count));
    guide_u = (double *)ARENA_CALLOC(arena, ncomponents, sizeof(*guide_u));
    guide_v = (double *)ARENA_CALLOC(arena, ncomponents, sizeof(*guide_v));
    safe_u = (double *)ARENA_CALLOC(arena, ncomponents, sizeof(*safe_u));
    safe_v = (double *)ARENA_CALLOC(arena, ncomponents, sizeof(*safe_v));
    coefficient_a = (double *)ARENA_CALLOC(
        arena, ncomponents, sizeof(*coefficient_a));
    coefficient_b = (double *)ARENA_CALLOC(
        arena, ncomponents, sizeof(*coefficient_b));
    norm = (double *)ARENA_CALLOC(arena, ncomponents, sizeof(*norm));
    positive = (size_t *)ARENA_CALLOC(
        arena, ncomponents, sizeof(*positive));
    negative = (size_t *)ARENA_CALLOC(
        arena, ncomponents, sizeof(*negative));
    reflect = (uint8_t *)ARENA_CALLOC(
        arena, ncomponents, sizeof(*reflect));
    memcpy(aligned, safe, system->nv * 2 * sizeof(*aligned));

    for (vertex = 0; vertex < system->nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        double gu = (double)guide[vertex * 2];
        double gv = (double)guide[vertex * 2 + 1];
        if (!isfinite(gu) || !isfinite(gv)) return -1;
        count[label] += 1.0;
        guide_u[label] += gu;
        guide_v[label] += gv;
        safe_u[label] += safe[vertex * 2];
        safe_v[label] += safe[vertex * 2 + 1];
    }
    for (component = 0; component < ncomponents; component++) {
        if (!(count[component] > 0.0)) return -1;
        guide_u[component] /= count[component];
        guide_v[component] /= count[component];
        safe_u[component] /= count[component];
        safe_v[component] /= count[component];
    }
    for (face_index = 0; face_index < system->nf; face_index++) {
        int32_t ia = system->faces[face_index * 3];
        int32_t ib = system->faces[face_index * 3 + 1];
        int32_t ic = system->faces[face_index * 3 + 2];
        size_t label = (size_t)system->vertex_component[(size_t)ia];
        double ab_u = (double)guide[(size_t)ib * 2] -
                      (double)guide[(size_t)ia * 2];
        double ab_v = (double)guide[(size_t)ib * 2 + 1] -
                      (double)guide[(size_t)ia * 2 + 1];
        double ac_u = (double)guide[(size_t)ic * 2] -
                      (double)guide[(size_t)ia * 2];
        double ac_v = (double)guide[(size_t)ic * 2 + 1] -
                      (double)guide[(size_t)ia * 2 + 1];
        double area2 = ab_u * ac_v - ab_v * ac_u;
        if (area2 > 0.0) positive[label]++;
        else if (area2 < 0.0) negative[label]++;
    }
    for (component = 0; component < ncomponents; component++) {
        usable[component] =
            (positive[component] == 0) != (negative[component] == 0);
        reflect[component] =
            usable[component] && negative[component] != 0;
    }
    for (vertex = 0; vertex < system->nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        double gx = (double)guide[vertex * 2] - guide_u[label];
        double gy = (double)guide[vertex * 2 + 1] - guide_v[label];
        double sx = safe[vertex * 2] - safe_u[label];
        double sy = safe[vertex * 2 + 1] - safe_v[label];
        if (reflect[label]) gy = -gy;
        norm[label] += gx * gx + gy * gy;
        coefficient_a[label] += gx * sx + gy * sy;
        coefficient_b[label] += gx * sy - gy * sx;
    }
    for (component = 0; component < ncomponents; component++) {
        if (!usable[component] || !(norm[component] > 1.0e-20) ||
            !isfinite(norm[component])) {
            usable[component] = 0;
            continue;
        }
        coefficient_a[component] /= norm[component];
        coefficient_b[component] /= norm[component];
        if (!isfinite(coefficient_a[component]) ||
            !isfinite(coefficient_b[component]) ||
            !(coefficient_a[component] * coefficient_a[component] +
              coefficient_b[component] * coefficient_b[component] >
              1.0e-24))
            usable[component] = 0;
    }
    for (vertex = 0; vertex < system->nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        double gx = (double)guide[vertex * 2] - guide_u[label];
        double gy = (double)guide[vertex * 2 + 1] - guide_v[label];
        double a = coefficient_a[label];
        double b = coefficient_b[label];
        if (!usable[label]) continue;
        if (reflect[label]) gy = -gy;
        aligned[vertex * 2] = safe_u[label] + a * gx - b * gy;
        aligned[vertex * 2 + 1] = safe_v[label] + b * gx + a * gy;
        if (!isfinite(aligned[vertex * 2]) ||
            !isfinite(aligned[vertex * 2 + 1]))
            return -1;
    }
    return 0;
}

static int sr_normalize_orientation(Arena_T arena, const SRSystem *system,
                                    double *uv)
{
    size_t *positive = (size_t *)ARENA_CALLOC(
        arena, system->ncomponents, sizeof(*positive));
    size_t *negative = (size_t *)ARENA_CALLOC(
        arena, system->ncomponents, sizeof(*negative));
    uint8_t *reflect = (uint8_t *)ARENA_CALLOC(
        arena, system->ncomponents, sizeof(*reflect));
    size_t face_index = 0;
    size_t vertex = 0;
    size_t component = 0;
    for (face_index = 0; face_index < system->nf; face_index++) {
        int32_t a = system->faces[face_index * 3];
        int32_t b = system->faces[face_index * 3 + 1];
        int32_t c = system->faces[face_index * 3 + 2];
        double area2 = sr_orient2(&uv[(size_t)a * 2],
                                  &uv[(size_t)b * 2],
                                  &uv[(size_t)c * 2]);
        size_t label = (size_t)system->vertex_component[(size_t)a];
        if (area2 > SR_AREA_EPS) positive[label]++;
        else if (area2 < -SR_AREA_EPS) negative[label]++;
        else {
            fprintf(stderr, "  [slim] collapsed input face %zu\n", face_index);
            return -1;
        }
    }
    for (component = 0; component < system->ncomponents; component++) {
        if (positive[component] != 0 && negative[component] != 0) {
            fprintf(stderr,
                    "  [slim] component %zu has mixed input orientation\n",
                    component);
            return -1;
        }
        reflect[component] = (uint8_t)(negative[component] != 0);
    }
    for (vertex = 0; vertex < system->nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        if (reflect[label]) uv[vertex * 2 + 1] = -uv[vertex * 2 + 1];
    }
    return 0;
}

/* 2x2 polar/SVD local step. W is symmetric and R is the closest
 * orientation-preserving rotation. */
static int sr_local_proxy(const double jacobian[4], const SlimRefineOpts *opts,
                          double weight[3], double rotation[4])
{
    double a = jacobian[0];
    double b = jacobian[1];
    double c = jacobian[2];
    double d = jacobian[3];
    double det = a * d - b * c;
    double p = a * a + b * b;
    double q = a * c + b * d;
    double r = c * c + d * d;
    double disc = hypot(p - r, 2.0 * q);
    double lambda0 = fmax(0.0, 0.5 * (p + r + disc));
    double sigma0 = sqrt(lambda0);
    /* det(J) is the product of the singular values. Subtracting nearly
     * equal eigenvalues of JJ^T erases the smaller, valid stretch. */
    double sigma1 = sigma0 > 0 ? fabs(det)/sigma0 : 0;
    double theta = 0.5 * atan2(2.0 * q, p - r);
    double cosine = cos(theta);
    double sine = sin(theta);
    double w0 = 0.0;
    double w1 = 0.0;
    double polar_norm = 0.0;
    if (!(det > SR_AREA_EPS) || !isfinite(det)) return -1;
    sigma0 = fmax(sigma0, opts->sigma_min);
    sigma1 = fmax(sigma1, opts->sigma_min);
    w0 = sqrt(1.0 + 1.0 / sigma0 + 1.0 / (sigma0 * sigma0) +
              1.0 / (sigma0 * sigma0 * sigma0));
    w1 = sqrt(1.0 + 1.0 / sigma1 + 1.0 / (sigma1 * sigma1) +
              1.0 / (sigma1 * sigma1 * sigma1));
    if (w0 > opts->weight_max) w0 = opts->weight_max;
    if (w1 > opts->weight_max) w1 = opts->weight_max;
    weight[0] = w0 * cosine * cosine + w1 * sine * sine;
    weight[1] = (w0 - w1) * cosine * sine;
    weight[2] = w0 * sine * sine + w1 * cosine * cosine;

    polar_norm = hypot(a + d, c - b);
    if (!(polar_norm > 1.0e-20) || !isfinite(polar_norm)) return -1;
    rotation[0] = (a + d) / polar_norm;
    rotation[1] = (b - c) / polar_norm;
    rotation[2] = (c - b) / polar_norm;
    rotation[3] = (a + d) / polar_norm;
    return 0;
}

static int sr_assemble(SRSystem *system, const double *uv,
                       const SlimRefineOpts *opts)
{
    size_t face_index = 0;
    memset(system->diagonal, 0, system->nv * 3 * sizeof(*system->diagonal));
    memset(system->edge_block, 0,
           system->ne * 3 * sizeof(*system->edge_block));
    memset(system->rhs, 0, system->nv * 2 * sizeof(*system->rhs));

    for (face_index = 0; face_index < system->nf; face_index++) {
        const SRFace *face = &system->face[face_index];
        const int32_t *indices = &system->faces[face_index * 3];
        double jacobian[4] = { 0.0, 0.0, 0.0, 0.0 };
        double weight[3] = { 0.0, 0.0, 0.0 };
        double rotation[4] = { 0.0, 0.0, 0.0, 0.0 };
        double h00 = 0.0;
        double h01 = 0.0;
        double h11 = 0.0;
        double target_x[2] = { 0.0, 0.0 };
        double target_y[2] = { 0.0, 0.0 };
        size_t corner = 0;
        sr_face_jacobian(face, indices, uv, jacobian);
        if (sr_local_proxy(jacobian, opts, weight, rotation) != 0)
            return -1;
        h00 = face->area *
            (weight[0] * weight[0] + weight[1] * weight[1]);
        h01 = face->area *
            (weight[0] * weight[1] + weight[1] * weight[2]);
        h11 = face->area *
            (weight[1] * weight[1] + weight[2] * weight[2]);
        target_x[0] = h00 * rotation[0] + h01 * rotation[2];
        target_x[1] = h01 * rotation[0] + h11 * rotation[2];
        target_y[0] = h00 * rotation[1] + h01 * rotation[3];
        target_y[1] = h01 * rotation[1] + h11 * rotation[3];

        for (corner = 0; corner < 3; corner++) {
            size_t v = (size_t)indices[corner];
            double coefficient = face->gx[corner] * face->gx[corner] +
                                 face->gy[corner] * face->gy[corner];
            system->diagonal[v * 3] += coefficient * h00;
            system->diagonal[v * 3 + 1] += coefficient * h01;
            system->diagonal[v * 3 + 2] += coefficient * h11;
            system->rhs[v * 2] += face->gx[corner] * target_x[0] +
                                  face->gy[corner] * target_y[0];
            system->rhs[v * 2 + 1] += face->gx[corner] * target_x[1] +
                                      face->gy[corner] * target_y[1];
        }
        for (corner = 0; corner < 3; corner++) {
            size_t next = (corner + 1) % 3;
            size_t edge = (size_t)system->face_edge[face_index * 3 + corner];
            double coefficient = face->gx[corner] * face->gx[next] +
                                 face->gy[corner] * face->gy[next];
            system->edge_block[edge * 3] += coefficient * h00;
            system->edge_block[edge * 3 + 1] += coefficient * h01;
            system->edge_block[edge * 3 + 2] += coefficient * h11;
        }
    }
    system->proximal_min = DBL_MAX;
    system->proximal_mean = 0.0;
    system->proximal_max = 0.0;
    for (size_t component = 0; component < system->ncomponents; component++) {
        size_t begin = system->component_vertex_offset[component];
        size_t end = system->component_vertex_offset[component + 1];
        double trace_sum = 0.0;
        double proximal = 0.0;
        size_t position = 0;
        for (position = begin; position < end; position++) {
            size_t v = (size_t)system->component_vertex[position];
            trace_sum += system->diagonal[v * 3] +
                         system->diagonal[v * 3 + 2];
        }
        proximal = opts->proximal * trace_sum /
                   fmax(2.0 * (double)(end - begin), 1.0);
        if (!(proximal > 0.0) || !isfinite(proximal)) return -1;
        if (proximal < system->proximal_min)
            system->proximal_min = proximal;
        if (proximal > system->proximal_max)
            system->proximal_max = proximal;
        system->proximal_mean += proximal * (double)(end - begin);
        for (position = begin; position < end; position++) {
            size_t v = (size_t)system->component_vertex[position];
            system->diagonal[v * 3] += proximal;
            system->diagonal[v * 3 + 2] += proximal;
            system->rhs[v * 2] += proximal * uv[v * 2];
            system->rhs[v * 2 + 1] += proximal * uv[v * 2 + 1];
        }
    }
    system->proximal_mean /= (double)system->nv;
    /* Exact symmetric Dirichlet elimination. Removing a fixed neighbour's
     * block moves its contribution to the free row RHS. The pinned rows are
     * identities; AMG/PCG still see a symmetric positive definite system. */
    if (opts->fixed_vertices != NULL) {
        for (size_t e = 0; e < system->ne; e++) {
            size_t a = (size_t)system->edge_a[e];
            size_t b = (size_t)system->edge_b[e];
            double *block = &system->edge_block[e * 3];
            int fa = opts->fixed_vertices[a] != 0;
            int fb = opts->fixed_vertices[b] != 0;
            if (!fa && fb) {
                system->rhs[a * 2] -= block[0] * uv[b * 2] + block[1] * uv[b * 2 + 1];
                system->rhs[a * 2 + 1] -= block[1] * uv[b * 2] + block[2] * uv[b * 2 + 1];
            } else if (fa && !fb) {
                system->rhs[b * 2] -= block[0] * uv[a * 2] + block[1] * uv[a * 2 + 1];
                system->rhs[b * 2 + 1] -= block[1] * uv[a * 2] + block[2] * uv[a * 2 + 1];
            }
            if (fa || fb) block[0] = block[1] = block[2] = 0.0;
        }
        for (size_t v = 0; v < system->nv; v++) {
            if (!opts->fixed_vertices[v]) continue;
            system->diagonal[v * 3] = system->diagonal[v * 3 + 2] = 1.0;
            system->diagonal[v * 3 + 1] = 0.0;
            system->rhs[v * 2] = uv[v * 2];
            system->rhs[v * 2 + 1] = uv[v * 2 + 1];
        }
    }
    return 0;
}

static int sr_compare_coarse_ref(const void *left, const void *right)
{
    const SRCoarseRef *a = (const SRCoarseRef *)left;
    const SRCoarseRef *b = (const SRCoarseRef *)right;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->fine_edge < b->fine_edge ? -1 :
          (a->fine_edge > b->fine_edge ? 1 : 0);
}

static double sr_block_norm(const double *block)
{
    return sqrt(block[0] * block[0] +
                2.0 * block[1] * block[1] +
                block[2] * block[2]);
}

static int sr_amg_build_rows(Arena_T arena, SRAMGLevel *level)
{
    size_t edge = 0;
    size_t vertex = 0;
    size_t *cursor = NULL;
    level->offset = (size_t *)ARENA_CALLOC(
        arena, level->n + 1, sizeof(*level->offset));
    for (edge = 0; edge < level->ne; edge++) {
        int32_t a = level->edge_a[edge];
        int32_t b = level->edge_b[edge];
        if (a < 0 || b < 0 || (size_t)a >= level->n ||
            (size_t)b >= level->n || a == b) return -1;
        level->offset[(size_t)a + 1]++;
        level->offset[(size_t)b + 1]++;
    }
    for (vertex = 0; vertex < level->n; vertex++)
        level->offset[vertex + 1] += level->offset[vertex];
    level->neighbor = (int32_t *)ARENA_ALLOC(
        arena, level->offset[level->n] * sizeof(*level->neighbor));
    level->adj_edge = (size_t *)ARENA_ALLOC(
        arena, level->offset[level->n] * sizeof(*level->adj_edge));
    cursor = (size_t *)ARENA_ALLOC(arena, level->n * sizeof(*cursor));
    memcpy(cursor, level->offset, level->n * sizeof(*cursor));
    for (edge = 0; edge < level->ne; edge++) {
        int32_t a = level->edge_a[edge];
        int32_t b = level->edge_b[edge];
        size_t pa = cursor[(size_t)a]++;
        size_t pb = cursor[(size_t)b]++;
        level->neighbor[pa] = b;
        level->adj_edge[pa] = edge;
        level->neighbor[pb] = a;
        level->adj_edge[pb] = edge;
    }
    return 0;
}

static int sr_amg_restrict_matrix(const SRAMGLevel *fine,
                                  SRAMGLevel *coarse)
{
    size_t vertex = 0;
    size_t edge = 0;
    memset(coarse->diagonal, 0,
           coarse->n * 3 * sizeof(*coarse->diagonal));
    memset(coarse->edge_block, 0,
           coarse->ne * 3 * sizeof(*coarse->edge_block));
    for (vertex = 0; vertex < fine->n; vertex++) {
        size_t aggregate = (size_t)fine->aggregate[vertex];
        size_t entry = 0;
        if (aggregate >= coarse->n) return -1;
        for (entry = 0; entry < 3; entry++)
            coarse->diagonal[aggregate * 3 + entry] +=
                fine->diagonal[vertex * 3 + entry];
    }
    for (edge = 0; edge < fine->ne; edge++) {
        size_t a = (size_t)fine->aggregate[(size_t)fine->edge_a[edge]];
        size_t b = (size_t)fine->aggregate[(size_t)fine->edge_b[edge]];
        size_t entry = 0;
        if (a == b) {
            for (entry = 0; entry < 3; entry++)
                coarse->diagonal[a * 3 + entry] +=
                    2.0 * fine->edge_block[edge * 3 + entry];
        } else {
            int32_t mapped = fine->edge_to_coarse[edge];
            if (mapped < 0 || (size_t)mapped >= coarse->ne) return -1;
            for (entry = 0; entry < 3; entry++)
                coarse->edge_block[(size_t)mapped * 3 + entry] +=
                    fine->edge_block[edge * 3 + entry];
        }
    }
    return 0;
}

static int sr_amg_build_next(Arena_T arena, SRAMGLevel *fine,
                             SRAMGLevel *coarse)
{
    int32_t *aggregate = NULL;
    SRCoarseRef *reference = NULL;
    size_t reference_count = 0;
    size_t vertex = 0;
    size_t edge = 0;
    size_t coarse_count = 0;
    size_t coarse_edges = 0;
    if (fine->n == 0 || fine->n > (size_t)INT32_MAX) return -1;
    aggregate = (int32_t *)ARENA_ALLOC(
        arena, fine->n * sizeof(*aggregate));
    for (vertex = 0; vertex < fine->n; vertex++) aggregate[vertex] = -1;
    for (vertex = 0; vertex < fine->n; vertex++) {
        int32_t member[SR_AMG_AGGREGATE_SIZE];
        int member_count = 1;
        int member_index = 0;
        if (aggregate[vertex] >= 0) continue;
        if (coarse_count > (size_t)INT32_MAX) return -1;
        aggregate[vertex] = (int32_t)coarse_count;
        member[0] = (int32_t)vertex;
        while (member_count < SR_AMG_AGGREGATE_SIZE) {
            int32_t best = -1;
            double best_score = -1.0;
            for (member_index = 0; member_index < member_count;
                 member_index++) {
                size_t source = (size_t)member[member_index];
                size_t position = 0;
                for (position = fine->offset[source];
                     position < fine->offset[source + 1]; position++) {
                    int32_t candidate = fine->neighbor[position];
                    size_t candidate_vertex = (size_t)candidate;
                    size_t candidate_edge = fine->adj_edge[position];
                    double denominator = 0.0;
                    double score = 0.0;
                    if (candidate < 0 || aggregate[candidate_vertex] >= 0)
                        continue;
                    denominator = sqrt(fmax(
                        sr_block_norm(&fine->diagonal[source * 3]) *
                        sr_block_norm(&fine->diagonal[candidate_vertex * 3]),
                        1.0e-30));
                    score = sr_block_norm(
                        &fine->edge_block[candidate_edge * 3]) / denominator;
                    if (score > best_score) {
                        best_score = score;
                        best = candidate;
                    }
                }
            }
            if (best < 0) break;
            aggregate[(size_t)best] = (int32_t)coarse_count;
            member[member_count++] = best;
        }
        coarse_count++;
    }
    fine->aggregate = aggregate;
    fine->ncoarse = coarse_count;
    fine->edge_to_coarse = (int32_t *)ARENA_ALLOC(
        arena, fine->ne * sizeof(*fine->edge_to_coarse));
    for (edge = 0; edge < fine->ne; edge++) fine->edge_to_coarse[edge] = -1;
    reference = (SRCoarseRef *)malloc(
        (fine->ne ? fine->ne : 1) * sizeof(*reference));
    if (reference == NULL) return -1;
    for (edge = 0; edge < fine->ne; edge++) {
        int32_t a = aggregate[(size_t)fine->edge_a[edge]];
        int32_t b = aggregate[(size_t)fine->edge_b[edge]];
        if (a != b) {
            reference[reference_count].key = sr_edge_key(a, b);
            reference[reference_count].fine_edge = edge;
            reference_count++;
        }
    }
    qsort(reference, reference_count, sizeof(*reference),
          sr_compare_coarse_ref);
    for (edge = 0; edge < reference_count;) {
        size_t end = edge + 1;
        while (end < reference_count &&
               reference[end].key == reference[edge].key) end++;
        coarse_edges++;
        edge = end;
    }
    memset(coarse, 0, sizeof(*coarse));
    coarse->n = coarse_count;
    coarse->ne = coarse_edges;
    coarse->edge_a = (int32_t *)ARENA_ALLOC(
        arena, coarse_edges * sizeof(*coarse->edge_a));
    coarse->edge_b = (int32_t *)ARENA_ALLOC(
        arena, coarse_edges * sizeof(*coarse->edge_b));
    coarse->diagonal = (double *)ARENA_ALLOC(
        arena, coarse_count * 3 * sizeof(*coarse->diagonal));
    coarse->edge_block = (double *)ARENA_ALLOC(
        arena, coarse_edges * 3 * sizeof(*coarse->edge_block));
    coarse_edges = 0;
    for (edge = 0; edge < reference_count;) {
        size_t end = edge + 1;
        while (end < reference_count &&
               reference[end].key == reference[edge].key) end++;
        coarse->edge_a[coarse_edges] =
            (int32_t)(reference[edge].key >> 32);
        coarse->edge_b[coarse_edges] =
            (int32_t)(reference[edge].key & UINT64_C(0xffffffff));
        for (size_t at = edge; at < end; at++)
            fine->edge_to_coarse[reference[at].fine_edge] =
                (int32_t)coarse_edges;
        coarse_edges++;
        edge = end;
    }
    free(reference);
    if (sr_amg_build_rows(arena, coarse) != 0) return -1;
    return sr_amg_restrict_matrix(fine, coarse);
}

static int sr_amg_factor_coarsest(SRAMGLevel *level)
{
    size_t dimension = level->n * 2;
    size_t row = 0;
    size_t edge = 0;
    double *matrix = level->cholesky;
    if (dimension == 0 || dimension > SIZE_MAX / dimension ||
        matrix == NULL) return -1;
    memset(matrix, 0, dimension * dimension * sizeof(*matrix));
    for (row = 0; row < level->n; row++) {
        matrix[(row * 2) * dimension + row * 2] =
            level->diagonal[row * 3];
        matrix[(row * 2) * dimension + row * 2 + 1] =
            level->diagonal[row * 3 + 1];
        matrix[(row * 2 + 1) * dimension + row * 2] =
            level->diagonal[row * 3 + 1];
        matrix[(row * 2 + 1) * dimension + row * 2 + 1] =
            level->diagonal[row * 3 + 2];
    }
    for (edge = 0; edge < level->ne; edge++) {
        size_t a = (size_t)level->edge_a[edge];
        size_t b = (size_t)level->edge_b[edge];
        const double *block = &level->edge_block[edge * 3];
        matrix[(a * 2) * dimension + b * 2] += block[0];
        matrix[(a * 2) * dimension + b * 2 + 1] += block[1];
        matrix[(a * 2 + 1) * dimension + b * 2] += block[1];
        matrix[(a * 2 + 1) * dimension + b * 2 + 1] += block[2];
        matrix[(b * 2) * dimension + a * 2] += block[0];
        matrix[(b * 2) * dimension + a * 2 + 1] += block[1];
        matrix[(b * 2 + 1) * dimension + a * 2] += block[1];
        matrix[(b * 2 + 1) * dimension + a * 2 + 1] += block[2];
    }
    for (row = 0; row < dimension; row++) {
        size_t column = 0;
        for (column = 0; column <= row; column++) {
            double value = matrix[row * dimension + column];
            size_t k = 0;
            for (k = 0; k < column; k++)
                value -= matrix[row * dimension + k] *
                         matrix[column * dimension + k];
            if (row == column) {
                double floor_value = 1.0e-14 *
                    (fabs(matrix[row * dimension + row]) + 1.0);
                if (value < floor_value) value = floor_value;
                matrix[row * dimension + column] = sqrt(value);
            } else {
                matrix[row * dimension + column] =
                    value / matrix[column * dimension + column];
            }
        }
        for (column = row + 1; column < dimension; column++)
            matrix[row * dimension + column] = 0.0;
    }
    return 0;
}

static int sr_amg_update(SRAMG *amg)
{
    int level = 0;
    for (level = 0; level + 1 < amg->levels; level++)
        if (sr_amg_restrict_matrix(&amg->level[level],
                                   &amg->level[level + 1]) != 0) return -1;
    return sr_amg_factor_coarsest(&amg->level[amg->levels - 1]);
}

static int sr_amg_build(Arena_T arena, const SRSystem *system, SRAMG *amg)
{
    SRAMGLevel *fine = NULL;
    size_t dimension = 0;
    memset(amg, 0, sizeof(*amg));
    amg->threads = system->threads;
    fine = &amg->level[0];
    fine->n = system->nv;
    fine->ne = system->ne;
    fine->edge_a = system->edge_a;
    fine->edge_b = system->edge_b;
    fine->offset = system->offset;
    fine->neighbor = system->neighbor;
    fine->adj_edge = system->adj_edge;
    fine->diagonal = system->diagonal;
    fine->edge_block = system->edge_block;
    amg->levels = 1;
    while (amg->levels < SR_AMG_MAX_LEVELS &&
           fine->n > SR_AMG_COARSE_VERTICES) {
        SRAMGLevel *coarse = &amg->level[amg->levels];
        if (sr_amg_build_next(arena, fine, coarse) != 0) return -1;
        amg->levels++;
        if (coarse->n >= fine->n ||
            coarse->n > (size_t)(0.95 * (double)fine->n)) break;
        fine = coarse;
    }
    fine = &amg->level[amg->levels - 1];
    dimension = fine->n * 2;
    if (dimension == 0 || dimension > SIZE_MAX / dimension ||
        dimension * dimension > SIZE_MAX / sizeof(*fine->cholesky))
        return -1;
    fine->cholesky = (double *)ARENA_ALLOC(
        arena, dimension * dimension * sizeof(*fine->cholesky));
    return sr_amg_update(amg);
}

static void sr_block_solve(const double *block, double r0, double r1,
                           double *x0, double *x1)
{
    double a = block[0];
    double b = block[1];
    double d = block[2];
    double determinant = a * d - b * b;
    double floor_value = 1.0e-18 * (fabs(a * d) + b * b + 1.0);
    if (!(determinant > floor_value) || !isfinite(determinant)) {
        *x0 = r0 / fmax(a, 1.0e-12);
        *x1 = r1 / fmax(d, 1.0e-12);
    } else {
        *x0 = (d * r0 - b * r1) / determinant;
        *x1 = (-b * r0 + a * r1) / determinant;
    }
}

static void sr_amg_apply_level(const SRAMG *amg, int level_index,
                               const double *x, double *result)
{
    const SRAMGLevel *level = &amg->level[level_index];
    int ii = 0;
    int count = (int)level->n;
    int threads = amg->threads;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1 && count >= 1000)
#endif
    for (ii = 0; ii < count; ii++) {
        size_t vertex = (size_t)ii;
        const double *diagonal = &level->diagonal[vertex * 3];
        double y0 = diagonal[0] * x[vertex * 2] +
                    diagonal[1] * x[vertex * 2 + 1];
        double y1 = diagonal[1] * x[vertex * 2] +
                    diagonal[2] * x[vertex * 2 + 1];
        size_t position = 0;
        for (position = level->offset[vertex];
             position < level->offset[vertex + 1]; position++) {
            size_t neighbor = (size_t)level->neighbor[position];
            size_t edge = level->adj_edge[position];
            const double *block = &level->edge_block[edge * 3];
            y0 += block[0] * x[neighbor * 2] +
                  block[1] * x[neighbor * 2 + 1];
            y1 += block[1] * x[neighbor * 2] +
                  block[2] * x[neighbor * 2 + 1];
        }
        result[vertex * 2] = y0;
        result[vertex * 2 + 1] = y1;
    }
}

static void sr_amg_smooth(const SRAMG *amg, int level_index,
                          const double *rhs, double *x, double *temporary)
{
    const SRAMGLevel *level = &amg->level[level_index];
    int sweep = 0;
    for (sweep = 0; sweep < SR_AMG_SMOOTH_SWEEPS; sweep++) {
        int ii = 0;
        int count = (int)level->n;
        int threads = amg->threads;
        sr_amg_apply_level(amg, level_index, x, temporary);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1 && count >= 1000)
#endif
        for (ii = 0; ii < count; ii++) {
            size_t vertex = (size_t)ii;
            double correction0 = 0.0;
            double correction1 = 0.0;
            sr_block_solve(&level->diagonal[vertex * 3],
                           rhs[vertex * 2] - temporary[vertex * 2],
                           rhs[vertex * 2 + 1] - temporary[vertex * 2 + 1],
                           &correction0, &correction1);
            x[vertex * 2] += SR_AMG_OMEGA * correction0;
            x[vertex * 2 + 1] += SR_AMG_OMEGA * correction1;
        }
    }
}

static int sr_amg_solve_coarsest(const SRAMGLevel *level,
                                 const double *rhs, double *x,
                                 double *temporary)
{
    size_t dimension = level->n * 2;
    size_t row = 0;
    const double *factor = level->cholesky;
    if (factor == NULL || dimension == 0) return -1;
    for (row = 0; row < dimension; row++) {
        size_t column = 0;
        double value = rhs[row];
        double diagonal = factor[row * dimension + row];
        for (column = 0; column < row; column++)
            value -= factor[row * dimension + column] * temporary[column];
        if (!(diagonal > 0.0) || !isfinite(diagonal)) return -1;
        temporary[row] = value / diagonal;
    }
    for (row = dimension; row > 0; row--) {
        size_t index = row - 1;
        size_t column = 0;
        double value = temporary[index];
        double diagonal = factor[index * dimension + index];
        for (column = index + 1; column < dimension; column++)
            value -= factor[column * dimension + index] * x[column];
        x[index] = value / diagonal;
    }
    return 0;
}

static int sr_amg_vcycle(const SRAMG *amg, SRAMGWork *work,
                         int level_index)
{
    const SRAMGLevel *level = &amg->level[level_index];
    double *x = work->x[level_index];
    double *rhs = work->r[level_index];
    double *temporary = work->temporary[level_index];
    if (level_index + 1 == amg->levels)
        return sr_amg_solve_coarsest(level, rhs, x, temporary);

    sr_amg_smooth(amg, level_index, rhs, x, temporary);
    sr_amg_apply_level(amg, level_index, x, temporary);
    {
        const SRAMGLevel *coarse = &amg->level[level_index + 1];
        double *coarse_rhs = work->r[level_index + 1];
        double *coarse_x = work->x[level_index + 1];
        size_t vertex = 0;
        memset(coarse_rhs, 0, coarse->n * 2 * sizeof(*coarse_rhs));
        for (vertex = 0; vertex < level->n; vertex++) {
            size_t aggregate = (size_t)level->aggregate[vertex];
            coarse_rhs[aggregate * 2] +=
                rhs[vertex * 2] - temporary[vertex * 2];
            coarse_rhs[aggregate * 2 + 1] +=
                rhs[vertex * 2 + 1] - temporary[vertex * 2 + 1];
        }
        memset(coarse_x, 0, coarse->n * 2 * sizeof(*coarse_x));
        if (sr_amg_vcycle(amg, work, level_index + 1) != 0) return -1;
        for (vertex = 0; vertex < level->n; vertex++) {
            size_t aggregate = (size_t)level->aggregate[vertex];
            x[vertex * 2] += coarse_x[aggregate * 2];
            x[vertex * 2 + 1] += coarse_x[aggregate * 2 + 1];
        }
    }
    sr_amg_smooth(amg, level_index, rhs, x, temporary);
    return 0;
}

static int sr_amg_work_init(Arena_T arena, const SRAMG *amg,
                            SRAMGWork *work)
{
    int level_index = 0;
    memset(work, 0, sizeof(*work));
    if (amg == NULL || amg->levels <= 0) return -1;
    for (level_index = 0; level_index < amg->levels; level_index++) {
        size_t count = amg->level[level_index].n * 2;
        work->temporary[level_index] = (double *)ARENA_ALLOC(
            arena, count * sizeof(*work->temporary[level_index]));
        if (level_index > 0) {
            work->x[level_index] = (double *)ARENA_ALLOC(
                arena, count * sizeof(*work->x[level_index]));
            work->r[level_index] = (double *)ARENA_ALLOC(
                arena, count * sizeof(*work->r[level_index]));
        }
    }
    return 0;
}

static int sr_amg_precondition(const SRAMG *amg, SRAMGWork *work,
                               const double *residual, double *result)
{
    size_t count = amg->level[0].n * 2;
    work->x[0] = result;
    work->r[0] = (double *)residual;
    memset(result, 0, count * sizeof(*result));
    return sr_amg_vcycle(amg, work, 0);
}

static void sr_apply(const SRSystem *system, const double *x, double *result)
{
    int ii = 0;
    int count = (int)system->nv;
    int threads = system->threads;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
    for (ii = 0; ii < count; ii++) {
        size_t vertex = (size_t)ii;
        double d00 = system->diagonal[vertex * 3];
        double d01 = system->diagonal[vertex * 3 + 1];
        double d11 = system->diagonal[vertex * 3 + 2];
        double y0 = d00 * x[vertex * 2] + d01 * x[vertex * 2 + 1];
        double y1 = d01 * x[vertex * 2] + d11 * x[vertex * 2 + 1];
        size_t position = 0;
        for (position = system->offset[vertex];
             position < system->offset[vertex + 1]; position++) {
            size_t neighbor = (size_t)system->neighbor[position];
            size_t edge = system->adj_edge[position];
            double b00 = system->edge_block[edge * 3];
            double b01 = system->edge_block[edge * 3 + 1];
            double b11 = system->edge_block[edge * 3 + 2];
            y0 += b00 * x[neighbor * 2] + b01 * x[neighbor * 2 + 1];
            y1 += b01 * x[neighbor * 2] + b11 * x[neighbor * 2 + 1];
        }
        result[vertex * 2] = y0;
        result[vertex * 2 + 1] = y1;
    }
}

static double sr_dot(const double *a, const double *b, size_t count,
                     int threads)
{
    double result = 0.0;
    int ii = 0;
    int n = (int)count;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:result) num_threads(threads) if(threads > 1)
#endif
    for (ii = 0; ii < n; ii++) result += a[(size_t)ii] * b[(size_t)ii];
    return result;
}

static void sr_block_precondition(const SRSystem *system,
                                  const double *residual, double *result)
{
    int ii = 0;
    int count = (int)system->nv;
    int threads = system->threads;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
    for (ii = 0; ii < count; ii++) {
        size_t vertex = (size_t)ii;
        sr_block_solve(&system->diagonal[vertex * 3],
                       residual[vertex * 2], residual[vertex * 2 + 1],
                       &result[vertex * 2], &result[vertex * 2 + 1]);
    }
}

static int sr_pcg(Arena_T arena, const SRSystem *system,
                   const SRAMG *amg, const SlimRefineOpts *opts, double *x,
                   int *out_iterations, double *out_relative_residual)
{
    Arena_Mark mark = Arena_save(arena);
    size_t count = system->nv * 2;
    double *residual = (double *)ARENA_ALLOC(arena, count * sizeof(*residual));
    double *preconditioned = (double *)ARENA_ALLOC(
        arena, count * sizeof(*preconditioned));
    double *direction = (double *)ARENA_ALLOC(arena, count * sizeof(*direction));
    double *product = (double *)ARENA_ALLOC(arena, count * sizeof(*product));
    SRAMGWork amg_work;
    double correction_norm = 0.0;
    double rz = 0.0;
    double relative = 0.0;
    int iteration = 0;
    int threads = system->threads;
    int ii = 0;
    int n = (int)count;
    int status = 0;

    memset(&amg_work, 0, sizeof(amg_work));
    if (amg != NULL && sr_amg_work_init(arena, amg, &amg_work) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }

    sr_apply(system, x, product);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
    for (ii = 0; ii < n; ii++)
        residual[(size_t)ii] = system->rhs[(size_t)ii] - product[(size_t)ii];
    /* Solve the correction to the supplied iterate. The absolute RHS can
     * contain an arbitrarily large gauge/pin coordinate, or be zero while
     * the iterate still has a nonzero residual. Neither means convergence. */
    correction_norm = sqrt(sr_dot(residual,residual,count,threads));
    if (!isfinite(correction_norm)) { Arena_restore(arena,mark); return -1; }
    if (correction_norm == 0) {
        if (out_iterations != NULL) *out_iterations = 0;
        if (out_relative_residual != NULL) *out_relative_residual = 0.0;
        Arena_restore(arena, mark);
        return 0;
    }
    relative = 1;
    if (amg != NULL) status = sr_amg_precondition(
        amg, &amg_work, residual, preconditioned);
    else sr_block_precondition(system, residual, preconditioned);
    rz = sr_dot(residual, preconditioned, count, threads);
    if (!(rz > 0.0) || !isfinite(rz)) status = -1;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
    for (ii = 0; ii < n; ii++) direction[(size_t)ii] = preconditioned[(size_t)ii];

    for (iteration = 0; status == 0 && iteration < opts->pcg_iterations;
         iteration++) {
        double p_ap = 0.0;
        double alpha = 0.0;
        double rz_next = 0.0;
        double beta = 0.0;
        sr_apply(system, direction, product);
        p_ap = sr_dot(direction, product, count, threads);
        if (!(p_ap > 0.0) || !isfinite(p_ap)) {
            status = -1;
            break;
        }
        alpha = rz / p_ap;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
        for (ii = 0; ii < n; ii++) {
            x[(size_t)ii] += alpha * direction[(size_t)ii];
            residual[(size_t)ii] -= alpha * product[(size_t)ii];
        }
        relative = sqrt(sr_dot(residual, residual, count, threads)) / correction_norm;
        if (!isfinite(relative)) {
            status = -1;
            break;
        }
        if (relative <= opts->pcg_tolerance) {
            iteration++;
            break;
        }
        if (amg != NULL) status = sr_amg_precondition(
            amg, &amg_work, residual, preconditioned);
        else sr_block_precondition(system, residual, preconditioned);
        if (status != 0) break;
        rz_next = sr_dot(residual, preconditioned, count, threads);
        if (!(rz_next > 0.0) || !isfinite(rz_next)) {
            status = -1;
            break;
        }
        beta = rz_next / rz;
        rz = rz_next;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads > 1)
#endif
        for (ii = 0; ii < n; ii++)
            direction[(size_t)ii] = preconditioned[(size_t)ii] +
                                    beta * direction[(size_t)ii];
    }
    if (out_iterations != NULL) *out_iterations = iteration;
    if (out_relative_residual != NULL) *out_relative_residual = relative;
    Arena_restore(arena, mark);
    return status;
}

static int sr_compare_pack(const void *left, const void *right)
{
    const SRPack *a = (const SRPack *)left;
    const SRPack *b = (const SRPack *)right;
    if (a->faces != b->faces) return a->faces > b->faces ? -1 : 1;
    return a->component < b->component ? -1 :
          (a->component > b->component ? 1 : 0);
}

static int sr_pack_strip(Arena_T arena, const SRSystem *system, double padding,
                         double *uv, double *out_width, double *out_height)
{
    SRPack *pack = (SRPack *)ARENA_ALLOC(
        arena, system->ncomponents * sizeof(*pack));
    double *shift_u = (double *)ARENA_ALLOC(
        arena, system->ncomponents * sizeof(*shift_u));
    double *shift_v = (double *)ARENA_ALLOC(
        arena, system->ncomponents * sizeof(*shift_v));
    double cursor = 0.0;
    double height = 0.0;
    size_t component = 0;
    size_t vertex = 0;
    if (padding < 0.0) padding = 0.0;
    for (component = 0; component < system->ncomponents; component++) {
        pack[component].component = (int32_t)component;
        pack[component].faces = system->component_faces[component];
        pack[component].min_u = DBL_MAX;
        pack[component].max_u = -DBL_MAX;
        pack[component].min_v = DBL_MAX;
        pack[component].max_v = -DBL_MAX;
    }
    for (vertex = 0; vertex < system->nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        double u = uv[vertex * 2];
        double v = uv[vertex * 2 + 1];
        if (u < pack[label].min_u) pack[label].min_u = u;
        if (u > pack[label].max_u) pack[label].max_u = u;
        if (v < pack[label].min_v) pack[label].min_v = v;
        if (v > pack[label].max_v) pack[label].max_v = v;
    }
    qsort(pack, system->ncomponents, sizeof(*pack), sr_compare_pack);
    for (component = 0; component < system->ncomponents; component++) {
        size_t label = (size_t)pack[component].component;
        double width = pack[component].max_u - pack[component].min_u;
        double local_height = pack[component].max_v - pack[component].min_v;
        if (!(width >= 0.0) || !(local_height >= 0.0) ||
            !isfinite(width) || !isfinite(local_height)) return -1;
        shift_u[label] = cursor - pack[component].min_u;
        shift_v[label] = -pack[component].min_v;
        cursor += width;
        if (component + 1 < system->ncomponents) cursor += padding;
        if (local_height > height) height = local_height;
    }
    for (vertex = 0; vertex < system->nv; vertex++) {
        size_t label = (size_t)system->vertex_component[vertex];
        uv[vertex * 2] += shift_u[label];
        uv[vertex * 2 + 1] += shift_v[label];
    }
    if (out_width != NULL) *out_width = cursor;
    if (out_height != NULL) *out_height = height;
    return 0;
}

static void sr_bounds(const double *uv, size_t nv,
                      double *out_width, double *out_height)
{
    double min_u = DBL_MAX;
    double max_u = -DBL_MAX;
    double min_v = DBL_MAX;
    double max_v = -DBL_MAX;
    size_t vertex = 0;
    for (vertex = 0; vertex < nv; vertex++) {
        double u = uv[vertex * 2];
        double v = uv[vertex * 2 + 1];
        if (u < min_u) min_u = u;
        if (u > max_u) max_u = u;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    if (out_width != NULL) *out_width = max_u - min_u;
    if (out_height != NULL) *out_height = max_v - min_v;
}

void SlimRefine_defaults(SlimRefineOpts *opts)
{
    if (opts == NULL) return;
    opts->iterations = 20;
    opts->homotopy_stages = 1;
    opts->final_iterations = 0;
    opts->pcg_iterations = 120;
    opts->pcg_tolerance = 1.0e-5;
    opts->energy_tolerance = 1.0e-6;
    opts->sigma_min = 1.0e-8;
    opts->weight_max = 1.0e8;
    opts->proximal = 1.0e-5;
    opts->line_search = 24;
    opts->use_amg = 1;
    opts->local_steps = 0;
    opts->local_step_min_faces = 100000;
    opts->guide_uv = NULL;
    opts->metric_uv = NULL;
    opts->guide_warm_start = 0;
    opts->guard_boundary = 1;
    opts->strip_pack = 1;
    opts->padding = 20.0;
    opts->threads = 0;
    opts->verbose = 0;
    opts->scaffold_uv = NULL;
    opts->scaffold_first_face = 0;
    opts->scaffold_weight = 1.0;
    opts->fixed_vertices = NULL;
}

int SlimRefine_run_double(Arena_T arena,
                          const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const double *uv_in,
                          const SlimRefineOpts *options,
                          double *uv_out,
                          SlimRefineStats *stats)
{
    SlimRefineOpts opts;
    SRSystem system;
    SRAMG amg;
    SRMeasure before;
    SRMeasure current_measure;
    SRMeasure candidate_measure;
    SRMeasure after;
    SRFace *initial_face = NULL;
    SRFace *target_face = NULL;
    double *current = NULL;
    double *destination = NULL;
    double *trial = NULL;
    double *vertex_alpha = NULL;
    uint8_t *damp_vertex = NULL;
    SRMeasure *component_measure = NULL;
    SREdgeBox *boxes = NULL;
    size_t count = 0;
    size_t iteration = 0;
    size_t continuation_iterations = 0;
    size_t total_iterations = 0;
    size_t iterations_executed = 0;
    size_t boundary_before = 0;
    size_t boundary_after = 0;
    int status = -1;
    int amg_built = 0;
    int solver_failed = 0;
    double accepted_alpha_sum = 0.0;
    double accepted_alpha_min = DBL_MAX;
    double accepted_alpha_max = 0.0;
    double guide_alpha_sum = 0.0;
    double guide_alpha_min = DBL_MAX;
    double guide_alpha_max = 0.0;
    double local_alpha_sum = 0.0;
    double local_alpha_min = DBL_MAX;
    double local_alpha_max = 0.0;
    double total_start = ves_clock_sec();
    double section_start = total_start;

    if (arena == NULL || verts == NULL || faces == NULL || uv_in == NULL ||
        uv_out == NULL || nv == 0 || nf == 0 || nv > (size_t)INT_MAX / 2)
        return -1;
    if (options != NULL) opts = *options;
    else SlimRefine_defaults(&opts);
    if (opts.iterations < 0 || opts.homotopy_stages <= 0 ||
        opts.final_iterations < 0 ||
        opts.pcg_iterations <= 0 ||
        !(opts.pcg_tolerance > 0.0) || !(opts.energy_tolerance >= 0.0) ||
        !(opts.sigma_min > 0.0) || !(opts.weight_max >= 1.0) ||
        !(opts.proximal > 0.0) || opts.line_search <= 0 ||
        (opts.local_steps && opts.local_step_min_faces == 0)) return -1;
    if (opts.scaffold_uv != NULL &&
        (opts.scaffold_first_face == 0 || opts.scaffold_first_face >= nf ||
         !(opts.scaffold_weight > 0.0) || !isfinite(opts.scaffold_weight) ||
         opts.metric_uv != NULL || opts.homotopy_stages != 1)) return -1;
    if (opts.fixed_vertices != NULL &&
        (opts.strip_pack || opts.guide_warm_start || opts.homotopy_stages != 1)) return -1;
    if (opts.homotopy_stages > 1 &&
        (opts.iterations <= 0 ||
         (opts.guide_warm_start && opts.guide_uv != NULL))) {
        fprintf(stderr,
                "  [slim] metric homotopy requires positive iterations and "
                "cannot be combined with guide warm start\n");
        return -1;
    }
    if ((size_t)opts.iterations >
        SIZE_MAX / (size_t)opts.homotopy_stages) return -1;
    continuation_iterations = (size_t)opts.iterations *
                              (size_t)opts.homotopy_stages;
    if ((size_t)opts.final_iterations >
        SIZE_MAX - continuation_iterations) return -1;
    total_iterations = continuation_iterations +
                       (size_t)opts.final_iterations;
    if (stats != NULL) memset(stats, 0, sizeof(*stats));
    memset(&system, 0, sizeof(system));
    memset(&amg, 0, sizeof(amg));
    memset(&before, 0, sizeof(before));
    memset(&current_measure, 0, sizeof(current_measure));
    memset(&candidate_measure, 0, sizeof(candidate_measure));
    memset(&after, 0, sizeof(after));

    if (sr_build_system(arena, verts, nv, faces, nf, &opts,
                        opts.threads, &system) != 0)
        return -1;
    if (opts.homotopy_stages > 1) {
        initial_face = (SRFace *)ARENA_ALLOC(
            arena, nf * sizeof(*initial_face));
        target_face = (SRFace *)ARENA_ALLOC(
            arena, nf * sizeof(*target_face));
        memcpy(target_face, system.face, nf * sizeof(*target_face));
        for (size_t face_index = 0; face_index < nf; face_index++) {
            int32_t a = faces[face_index * 3];
            int32_t b = faces[face_index * 3 + 1];
            int32_t c = faces[face_index * 3 + 2];
            if (sr_face_rest_uv_double(&uv_in[(size_t)a * 2],
                                       &uv_in[(size_t)b * 2],
                                       &uv_in[(size_t)c * 2],
                                       &initial_face[face_index]) != 0) {
                fprintf(stderr,
                        "  [slim] injective scaffold has a degenerate face "
                        "at %zu; cannot construct metric homotopy\n",
                        face_index);
                return -1;
            }
        }
    }
    if (opts.metric_uv != NULL)
        fprintf(stderr,
                "  [slim] winding rest metric: %zu face(s), "
                "%zu isolated 3-D fallback(s)\n",
                system.metric_uv_faces, system.metric_fallback_faces);
    count = nv * 2;
    current = (double *)ARENA_ALLOC(arena, count * sizeof(*current));
    destination = (double *)ARENA_ALLOC(arena, count * sizeof(*destination));
    trial = (double *)ARENA_ALLOC(arena, count * sizeof(*trial));
    vertex_alpha = (double *)ARENA_ALLOC(
        arena, nv * sizeof(*vertex_alpha));
    damp_vertex = (uint8_t *)ARENA_ALLOC(
        arena, nv * sizeof(*damp_vertex));
    component_measure = (SRMeasure *)ARENA_ALLOC(
        arena, system.ncomponents * sizeof(*component_measure));
    boxes = (SREdgeBox *)ARENA_ALLOC(
        arena, (system.nb ? system.nb : 1) * sizeof(*boxes));
    memcpy(current, uv_in, count * sizeof(*current));
    if (sr_normalize_orientation(arena, &system, current) != 0) return -1;
    if (sr_measure(&system, current, &before) != 0 || before.flips != 0)
        return -1;
    current_measure = before;
    for (size_t component = 0; component < system.ncomponents; component++)
        if (sr_measure_component(&system, current, component,
                                 &component_measure[component]) != 0)
            return -1;
    if (opts.guard_boundary) {
        boundary_before = sr_boundary_intersections(&system, current, boxes);
        if (boundary_before != 0) {
            fprintf(stderr,
                    "  [slim] input boundary has %zu non-adjacent intersections\n",
                    boundary_before);
            return -1;
        }
    }
    if (stats != NULL) {
        stats->components = system.ncomponents;
        stats->edges = system.ne;
        stats->boundary_edges = system.nb;
        stats->metric_uv_faces = system.metric_uv_faces;
        stats->metric_fallback_faces = system.metric_fallback_faces;
        stats->energy_before = before.energy;
        stats->qc_mean_before = before.qc_mean;
        stats->qc_max_before = before.qc_max;
        stats->min_det_before = before.min_det;
        stats->flips_before = before.flips;
        stats->boundary_intersections_before = boundary_before;
        stats->guide_energy_after = before.energy;
        stats->homotopy_stages_run = opts.homotopy_stages == 1 ? 1 : 0;
        stats->setup_seconds = ves_clock_sec() - section_start;
    }

    /* Tutte/harmonic is only the injective certificate.  Before metric
     * optimization, move each chart toward the winding-informed shape as far
     * as the exact face-flip and simple-boundary guards permit. */
    if (opts.guide_warm_start && opts.guide_uv != NULL) {
        section_start = ves_clock_sec();
        uint8_t *guide_usable = (uint8_t *)ARENA_CALLOC(
            arena, system.ncomponents, sizeof(*guide_usable));
        size_t accepted_components = 0;
        size_t rejected_components = 0;
        size_t guide_trials = 0;
        if (sr_align_guide(arena, &system, current, opts.guide_uv,
                           destination, guide_usable) != 0)
            return -1;
        memcpy(trial, current, count * sizeof(*trial));
        for (size_t component = 0; component < system.ncomponents;
             component++) {
            size_t vertex_begin = system.component_vertex_offset[component];
            size_t vertex_end = system.component_vertex_offset[component + 1];
            double alpha = 0.0;
            int accepted = 0;
            if (guide_usable[component])
                alpha = sr_component_flip_step(
                    &system, current, destination, component);
            for (size_t line = 0;
                 guide_usable[component] && alpha > 0.0 &&
                 line < (size_t)opts.line_search; line++) {
                for (size_t position = vertex_begin; position < vertex_end;
                     position++) {
                    size_t vertex =
                        (size_t)system.component_vertex[position];
                    trial[vertex * 2] = current[vertex * 2] + alpha *
                        (destination[vertex * 2] - current[vertex * 2]);
                    trial[vertex * 2 + 1] = current[vertex * 2 + 1] + alpha *
                        (destination[vertex * 2 + 1] -
                         current[vertex * 2 + 1]);
                }
                guide_trials++;
                if (sr_measure_component(&system, trial, component,
                                         &candidate_measure) == 0 &&
                    candidate_measure.flips == 0 &&
                    (!opts.guard_boundary ||
                     sr_component_boundary_intersections(
                         &system, component, trial, boxes) == 0) &&
                    candidate_measure.energy <
                        component_measure[component].energy) {
                    accepted = 1;
                    break;
                }
                alpha *= 0.5;
            }
            if (accepted) {
                if (opts.verbose)
                    fprintf(stderr,
                            "    [slim guide] c%zu F=%zu alpha=%.3e "
                            "E=%.6g->%.6g qcmax=%.4g->%.4g\n",
                            component, system.component_faces[component],
                            alpha, component_measure[component].energy,
                            candidate_measure.energy,
                            component_measure[component].qc_max,
                            candidate_measure.qc_max);
                component_measure[component] = candidate_measure;
                accepted_components++;
                guide_alpha_sum += alpha;
                if (alpha < guide_alpha_min) guide_alpha_min = alpha;
                if (alpha > guide_alpha_max) guide_alpha_max = alpha;
            } else {
                rejected_components++;
                for (size_t position = vertex_begin; position < vertex_end;
                     position++) {
                    size_t vertex =
                        (size_t)system.component_vertex[position];
                    trial[vertex * 2] = current[vertex * 2];
                    trial[vertex * 2 + 1] = current[vertex * 2 + 1];
                }
            }
        }
        if (accepted_components != 0) {
            memcpy(current, trial, count * sizeof(*current));
            if (sr_measure(&system, current, &current_measure) != 0 ||
                current_measure.flips != 0 ||
                !(current_measure.energy < before.energy))
                return -1;
        }
        if (opts.verbose)
            fprintf(stderr,
                    "  [slim guide] accepted=%zu/%zu trials=%zu "
                    "alpha=%.3g/%.3g/%.3g E=%.9g->%.9g\n",
                    accepted_components, system.ncomponents, guide_trials,
                    accepted_components ? guide_alpha_min : 0.0,
                    accepted_components
                        ? guide_alpha_sum / (double)accepted_components : 0.0,
                    guide_alpha_max, before.energy, current_measure.energy);
        if (stats != NULL) {
            stats->guide_components_accepted = accepted_components;
            stats->guide_components_rejected = rejected_components;
            stats->guide_line_search_trials = guide_trials;
            stats->guide_alpha_min =
                accepted_components ? guide_alpha_min : 0.0;
            stats->guide_alpha_mean = accepted_components
                ? guide_alpha_sum / (double)accepted_components : 0.0;
            stats->guide_alpha_max = guide_alpha_max;
            stats->guide_energy_after = current_measure.energy;
            stats->guide_seconds += ves_clock_sec() - section_start;
        }
    }

    for (iteration = 0; iteration < total_iterations; iteration++) {
        size_t stage = iteration < continuation_iterations &&
                       opts.iterations > 0
            ? iteration / (size_t)opts.iterations
            : (size_t)opts.homotopy_stages - 1;
        size_t stage_iteration = iteration < continuation_iterations &&
                                 opts.iterations > 0
            ? iteration % (size_t)opts.iterations
            : (size_t)opts.iterations + iteration - continuation_iterations;
        size_t stage_iteration_limit = (size_t)opts.iterations +
            (stage + 1 == (size_t)opts.homotopy_stages
                ? (size_t)opts.final_iterations : 0);
        int final_stage = stage + 1 == (size_t)opts.homotopy_stages;
        int pcg_iterations = 0;
        double pcg_relative = 0.0;
        size_t accepted_components = 0;
        size_t rejected_components = 0;
        double iteration_alpha_sum = 0.0;
        double iteration_alpha_min = DBL_MAX;
        double iteration_alpha_max = 0.0;
        double old_energy = 0.0;
        if (initial_face != NULL && stage_iteration == 0) {
            double t = (double)(stage + 1) /
                       (double)opts.homotopy_stages;
            if (sr_blend_rest_faces(&system, initial_face, target_face, t) != 0 ||
                sr_measure(&system, current, &current_measure) != 0 ||
                current_measure.flips != 0)
                return -1;
            for (size_t component = 0;
                 component < system.ncomponents; component++)
                if (sr_measure_component(&system, current, component,
                                         &component_measure[component]) != 0)
                    return -1;
            if (stats != NULL) stats->homotopy_stages_run = stage + 1;
            if (opts.verbose)
                fprintf(stderr,
                        "  [slim homotopy] stage %zu/%d t=%.6f "
                        "E=%.9g qc=%.4g/%.4g\n",
                        stage + 1, opts.homotopy_stages, t,
                        current_measure.energy, current_measure.qc_mean,
                        current_measure.qc_max);
        }
        old_energy = current_measure.energy;
        iterations_executed++;
        if (stats != NULL) stats->iterations_run = iterations_executed;
        section_start = ves_clock_sec();
        if (sr_assemble(&system, current, &opts) != 0) {
            solver_failed = 1;
            break;
        }
        if (stats != NULL)
            stats->assembly_seconds += ves_clock_sec() - section_start;
        if (opts.verbose && stage_iteration == 0) {
            double trace_min = DBL_MAX;
            double trace_max = 0.0;
            double trace_sum = 0.0;
            for (size_t vertex = 0; vertex < system.nv; vertex++) {
                double trace = system.diagonal[vertex * 3] +
                               system.diagonal[vertex * 3 + 2];
                if (trace < trace_min) trace_min = trace;
                if (trace > trace_max) trace_max = trace;
                trace_sum += trace;
            }
            fprintf(stderr,
                    "  [slim] proxy diagonal trace min/mean/max="
                    "%.3e/%.3e/%.3e proximal(relative)=%.3e "
                    "actual=%.3e/%.3e/%.3e\n",
                    trace_min, trace_sum / (double)system.nv, trace_max,
                    opts.proximal, system.proximal_min,
                    system.proximal_mean, system.proximal_max);
        }
        if (opts.use_amg) {
            section_start = ves_clock_sec();
            int amg_status = amg_built
                ? sr_amg_update(&amg)
                : sr_amg_build(arena, &system, &amg);
            if (amg_status != 0) {
                solver_failed = 1;
                if (opts.verbose)
                    fprintf(stderr,
                            "  [slim] AMG construction/update failed at "
                            "iteration %zu\n", iteration + 1);
                break;
            }
            amg_built = 1;
            if (stats != NULL) {
                stats->amg_levels = amg.levels;
                stats->amg_coarse_vertices =
                    amg.level[amg.levels - 1].n;
            }
            if (opts.verbose && iteration == 0)
            {
                int level_index = 0;
                fprintf(stderr, "  [slim] AMG hierarchy:");
                for (level_index = 0; level_index < amg.levels;
                     level_index++)
                    fprintf(stderr, " %zu",
                            amg.level[level_index].n);
                fprintf(stderr, " vertices\n");
            }
            if (stats != NULL)
                stats->amg_seconds += ves_clock_sec() - section_start;
        }
        memcpy(destination, current, count * sizeof(*destination));
        section_start = ves_clock_sec();
        if (sr_pcg(arena, &system, opts.use_amg ? &amg : NULL,
                   &opts, destination,
                   &pcg_iterations, &pcg_relative) != 0) {
            solver_failed = 1;
            if (opts.verbose)
                fprintf(stderr, "  [slim] PCG breakdown at iteration %zu\n",
                        iteration + 1);
            break;
        }
        if (stats != NULL)
            stats->pcg_seconds += ves_clock_sec() - section_start;
        if (opts.fixed_vertices != NULL) {
            for (size_t v = 0; v < nv; v++) {
                if (!opts.fixed_vertices[v]) continue;
                destination[v * 2] = current[v * 2];
                destination[v * 2 + 1] = current[v * 2 + 1];
            }
        }
        if (stats != NULL) {
            stats->pcg_iterations_total += (size_t)pcg_iterations;
            if (pcg_iterations > stats->pcg_iterations_max)
                stats->pcg_iterations_max = pcg_iterations;
            if (pcg_relative > stats->pcg_relative_residual_max)
                stats->pcg_relative_residual_max = pcg_relative;
        }
        section_start = ves_clock_sec();
        memcpy(trial, current, count * sizeof(*trial));
        for (size_t component = 0; component < system.ncomponents;
             component++) {
            size_t vertex_begin = system.component_vertex_offset[component];
            size_t vertex_end = system.component_vertex_offset[component + 1];
            double alpha = sr_component_flip_step(
                &system, current, destination, component);
            double component_old_energy = component_measure[component].energy;
            double component_old_qc = component_measure[component].qc_max;
            size_t line = 0;
            int accepted = 0;
            int used_local = opts.local_steps &&
                system.component_faces[component] >=
                    opts.local_step_min_faces;
            double local_step_min = 0.0;
            double local_step_mean = 0.0;
            double local_step_max = 0.0;
            if (used_local) {
                size_t local_trials = 0;
                accepted = sr_component_local_step(
                    &system, current, destination, component, &opts,
                    component_measure[component].energy, trial,
                    vertex_alpha, damp_vertex, boxes, &candidate_measure,
                    &local_trials, &local_step_min, &local_step_mean,
                    &local_step_max);
                alpha = local_step_mean;
                if (stats != NULL) stats->line_search_trials += local_trials;
            }
            /* Local damping is only a proposal. It can exhaust its search
             * while the uniform, analytically bounded step still descends.
             * Try that step before declaring this component stalled. */
            if (!accepted) {
                used_local = 0;
                alpha = sr_component_flip_step(
                    &system, current, destination, component);
                for (line = 0; line < (size_t)opts.line_search; line++) {
                    size_t position = 0;
                    for (position = vertex_begin; position < vertex_end;
                         position++) {
                        size_t vertex =
                            (size_t)system.component_vertex[position];
                        trial[vertex * 2] = current[vertex * 2] + alpha *
                            (destination[vertex * 2] - current[vertex * 2]);
                        trial[vertex * 2 + 1] =
                            current[vertex * 2 + 1] + alpha *
                            (destination[vertex * 2 + 1] -
                             current[vertex * 2 + 1]);
                    }
                    if (stats != NULL) stats->line_search_trials++;
                    if (sr_measure_component(&system, trial, component,
                                             &candidate_measure) == 0 &&
                        candidate_measure.flips == 0 &&
                        (!opts.guard_boundary ||
                         sr_component_boundary_intersections(
                             &system, component, trial, boxes) == 0) &&
                        candidate_measure.energy <
                            component_measure[component].energy) {
                        accepted = 1;
                        break;
                    }
                    alpha *= 0.5;
                }
            }
            if (accepted) {
                component_measure[component] = candidate_measure;
                accepted_components++;
                iteration_alpha_sum += alpha;
                if (alpha < iteration_alpha_min) iteration_alpha_min = alpha;
                if (alpha > iteration_alpha_max) iteration_alpha_max = alpha;
                accepted_alpha_sum += alpha;
                if (alpha < accepted_alpha_min) accepted_alpha_min = alpha;
                if (alpha > accepted_alpha_max) accepted_alpha_max = alpha;
                if (used_local) {
                    if (stats != NULL)
                        stats->local_component_steps_accepted++;
                    local_alpha_sum += local_step_mean;
                    if (local_step_min < local_alpha_min)
                        local_alpha_min = local_step_min;
                    if (local_step_max > local_alpha_max)
                        local_alpha_max = local_step_max;
                }
                if (opts.verbose &&
                    (alpha < 1.0e-3 ||
                     system.component_faces[component] >= 10000)) {
                    if (used_local)
                        fprintf(stderr,
                                "    [slim] c%zu F=%zu local-alpha="
                                "%.3e/%.3e/%.3e E=%.6g->%.6g "
                                "qcmax=%.4g->%.4g\n",
                                component, system.component_faces[component],
                                local_step_min, local_step_mean,
                                local_step_max, component_old_energy,
                                candidate_measure.energy, component_old_qc,
                                candidate_measure.qc_max);
                    else
                        fprintf(stderr,
                                "    [slim] c%zu F=%zu alpha=%.3e "
                                "E=%.6g->%.6g qcmax=%.4g->%.4g\n",
                                component, system.component_faces[component],
                                alpha, component_old_energy,
                                candidate_measure.energy, component_old_qc,
                                candidate_measure.qc_max);
                }
            } else {
                size_t position = 0;
                rejected_components++;
                for (position = vertex_begin; position < vertex_end;
                     position++) {
                    size_t vertex = (size_t)system.component_vertex[position];
                    trial[vertex * 2] = current[vertex * 2];
                    trial[vertex * 2 + 1] = current[vertex * 2 + 1];
                }
                if (opts.verbose)
                    fprintf(stderr,
                            "    [slim] c%zu F=%zu rejected all %d "
                            "line-search trials\n",
                            component, system.component_faces[component],
                            opts.line_search);
            }
        }
        if (stats != NULL) {
            stats->component_steps_accepted += accepted_components;
            stats->component_steps_rejected += rejected_components;
            stats->line_search_seconds += ves_clock_sec() - section_start;
        }
        if (accepted_components == 0) {
            if (!final_stage) {
                iteration = (stage + 1) * (size_t)opts.iterations - 1;
                continue;
            }
            /* Exhausting a line search does not establish stationarity. */
            break;
        }
        memcpy(current, trial, count * sizeof(*current));
        if (sr_measure(&system, current, &current_measure) != 0 ||
            !(current_measure.energy < old_energy)) {
            solver_failed = 1;
            break;
        }
        if (stats != NULL) {
            stats->accepted_steps++;
        }
        if (opts.verbose) {
            fprintf(stderr,
                    "  [slim] stage %zu/%d iter %zu/%d: E %.9g -> %.9g "
                    "components=%zu/%zu alpha=%.3g/%.3g/%.3g "
                    "pcg=%d rel=%.3e qc=%.4g/%.4g\n",
                    stage + 1, opts.homotopy_stages,
                    stage_iteration + 1, (int)stage_iteration_limit,
                    old_energy, current_measure.energy,
                    accepted_components, system.ncomponents,
                    iteration_alpha_min,
                    iteration_alpha_sum / (double)accepted_components,
                    iteration_alpha_max,
                    pcg_iterations, pcg_relative, current_measure.qc_mean,
                    current_measure.qc_max);
        }
        if ((old_energy - current_measure.energy) /
                fmax(fabs(old_energy), 1.0) <= opts.energy_tolerance) {
            if (!final_stage) {
                iteration = (stage + 1) * (size_t)opts.iterations - 1;
                continue;
            }
            if (stats != NULL) stats->converged = 1;
            break;
        }
    }

    if (solver_failed) return -1;

    section_start = ves_clock_sec();
    if (target_face != NULL)
        memcpy(system.face, target_face, nf * sizeof(*target_face));
    if (sr_measure(&system, current, &after) != 0 || after.flips != 0 ||
        after.energy > before.energy * (1.0 + 1.0e-10)) {
        memcpy(current, uv_in, count * sizeof(*current));
        if (sr_normalize_orientation(arena, &system, current) != 0) return -1;
        after = before;
        if (stats != NULL) stats->reverted = 1;
    }
    if (opts.strip_pack) {
        if (sr_pack_strip(arena, &system, opts.padding, current,
                          stats != NULL ? &stats->atlas_width : NULL,
                          stats != NULL ? &stats->atlas_height : NULL) != 0)
            return -1;
    } else if (stats != NULL) {
        sr_bounds(current, nv, &stats->atlas_width, &stats->atlas_height);
    }
    boundary_after = opts.guard_boundary
        ? sr_boundary_intersections(&system, current, boxes) : 0;
    if (boundary_after != 0) {
        fprintf(stderr,
                "  [slim] post-refinement boundary audit found %zu intersections\n",
                boundary_after);
        return -1;
    }
    if (sr_measure(&system, current, &after) != 0 || after.flips != 0)
        return -1;
    memcpy(uv_out, current, count * sizeof(*uv_out));
    if (stats != NULL) {
        stats->accepted_alpha_min = stats->component_steps_accepted > 0
            ? accepted_alpha_min : 0.0;
        stats->accepted_alpha_mean = stats->component_steps_accepted > 0
            ? accepted_alpha_sum / (double)stats->component_steps_accepted
            : 0.0;
        stats->accepted_alpha_max = accepted_alpha_max;
        stats->local_vertex_alpha_min =
            stats->local_component_steps_accepted > 0
                ? local_alpha_min : 0.0;
        stats->local_vertex_alpha_mean =
            stats->local_component_steps_accepted > 0
                ? local_alpha_sum /
                  (double)stats->local_component_steps_accepted : 0.0;
        stats->local_vertex_alpha_max = local_alpha_max;
        stats->energy_after = after.energy;
        stats->qc_mean_after = after.qc_mean;
        stats->qc_max_after = after.qc_max;
        stats->min_det_after = after.min_det;
        stats->flips_after = after.flips;
        stats->boundary_intersections_after = boundary_after;
        stats->finalize_seconds = ves_clock_sec() - section_start;
        stats->total_seconds = ves_clock_sec() - total_start;
    }
    status = 0;
    return status;
}

static int sr_local_interval_selftest(void)
{
    int failures = 0;
    for (int shifted = 0; shifted < 2; shifted++) {
        Arena_T arena = Arena_new(); SRSystem system;
        SlimRefineOpts opts; SlimRefine_defaults(&opts);
        const float xyz[12] = {0,0,0, 1,0,0, 1,1,0, 0,1,0};
        const int32_t faces[6] = {0,1,2, 0,2,3};
        double current[8] = {0,0, 20,0, 20,.05, 0,.05};
        double destination[8] = {
            -6.248038358190012,-7.73428763658029,
            16.146750749448817,6.256435730826154,
            -.8912790027333588,5.531694344158108,
            -4.570491549248102,-3.7242425104600594};
        double control[8] = {0,0, 1,0, 1,1, 0,1}, trial[8], alpha[4];
        uint8_t damp[4]; SREdgeBox boxes[4]; SRMeasure before = {0}, after = {0};
        size_t trials; double amin,amean,amax;
        opts.guard_boundary = 1;
        for (int v = 0; v < 4; v++) {
            double u = shifted ? 1e6 : 0, w = shifted ? -2e6 : 0;
            current[2*v] += u; current[2*v+1] += w;
            destination[2*v] += u; destination[2*v+1] += w;
            control[2*v] += u; control[2*v+1] += w;
        }
        int ready = !sr_build_system(arena,xyz,4,faces,2,&opts,1,&system) &&
                 !sr_measure_component(&system,current,0,&before);
        int ok = ready;
        /* The old locally damped endpoint has positive areas 117.99 and
         * 12.78 and lower energy, but the first triangle flips twice en
         * route (minimum signed area -1.35 at t=.123). The initial uniform
         * flip bounds do not certify the subsequent nonuniform movement. */
        int accepted = ok && sr_component_local_step(&system,current,destination,0,&opts,
            before.energy,trial,alpha,damp,boxes,&after,&trials,&amin,&amean,&amax);
        if (accepted) ok = !after.flips && after.energy < before.energy &&
            UvGuard_interval(current,trial,faces,1) && UvGuard_interval(current,trial,faces+3,1);
        int safe = ready && sr_component_local_step(&system,current,control,0,&opts,
            before.energy,trial,alpha,damp,boxes,&after,&trials,&amin,&amean,&amax);
        ok = ok && safe && !after.flips && after.energy < before.energy &&
            UvGuard_interval(current,trial,faces,1) && UvGuard_interval(current,trial,faces+3,1);
        if (!ok) {
            fprintf(stderr,"[slim selftest] locally damped path crossed a fold (shifted %d, accepted %d, control %d)\n",shifted,accepted,safe);
            failures++;
        }
        Arena_dispose(&arena);
    }
    return failures;
}

int SlimRefine_selftest(void)
{
    enum { NX = 13, NY = 10 };
    const size_t nv = (size_t)NX * NY;
    const size_t nf = (size_t)(NX - 1) * (NY - 1) * 2;
    float *verts = (float *)malloc(nv * 3 * sizeof(*verts));
    int32_t *faces = (int32_t *)malloc(nf * 3 * sizeof(*faces));
    double *uv = (double *)malloc(nv * 2 * sizeof(*uv));
    double *output = (double *)malloc(nv * 2 * sizeof(*output));
    float *guide = (float *)malloc(nv * 2 * sizeof(*guide));
    Arena_T arena = Arena_new();
    Arena_T homotopy_arena = NULL;
    SlimRefineOpts opts;
    SlimRefineStats stats;
    size_t x = 0;
    size_t y = 0;
    size_t face = 0;
    int failures = sr_local_interval_selftest();
    {
        /* Original 3-D and float UV rest frames need the same thin-face
         * stability as the existing double scaffold frame. */
        const float a[3] = {0,0,0}, b[3] = {10000,0,0}, c[3] = {10000,.00001f,0};
        SRFace rest, uv_rest;
        double expected = .5*(double)b[0]*c[1];
        if (sr_face_rest(a,b,c,&rest) || sr_face_rest_uv(a,b,c,&uv_rest) ||
            fabs(rest.area-expected)>1e-12 || fabs(uv_rest.area-expected)>1e-12) {
            fprintf(stderr,"[slim selftest] thin original/float rest frame lost\n"); failures++;
        }
    }
    {
        /* A large gauge value must not hide an unsolved correction. A
         * zero right-hand side also still needs to correct nonzero x. */
        SRSystem system = {0}; SlimRefineOpts options;
        size_t offset[] = {0,0}; double diagonal[] = {1,0,1};
        double rhs[2], x[2]; int iterations; double relative;
        system.nv = 1; system.threads = 1; system.offset = offset;
        system.diagonal = diagonal; system.rhs = rhs;
        SlimRefine_defaults(&options); options.pcg_iterations = 4;
        for (int translated = 0; translated < 2; translated++) {
            rhs[0] = translated ? 1e8 : 0; rhs[1] = translated ? -2e8 : 0;
            x[0] = rhs[0]+1; x[1] = rhs[1]-2;
            if (sr_pcg(arena,&system,NULL,&options,x,&iterations,&relative) ||
                iterations == 0 || x[0] != rhs[0] || x[1] != rhs[1]) {
                fprintf(stderr,"[slim selftest] PCG correction skipped (translated %d)\n",translated); failures++;
            }
        }
    }
    {
        /* Subtracting eigenvalues of JJ^T erases a valid minor stretch. */
        const double j[] = {1e9,0,0,1}; double weight[3], rotation[4];
        SlimRefineOpts options; SlimRefine_defaults(&options);
        if (sr_local_proxy(j,&options,weight,rotation) ||
            fabs(weight[2]-2)>1e-12 || fabs(rotation[0]-1)>1e-12 || fabs(rotation[3]-1)>1e-12) {
            fprintf(stderr,"[slim selftest] anisotropic minor stretch lost\n"); failures++;
        }
    }
    {
        const double a[2] = {0.0, 0.0};
        const double b[2] = {10000.0, 0.0};
        const double c[2] = {10000.0, 0.00001};
        SRFace thin;
        if (sr_face_rest_uv_double(a, b, c, &thin) != 0 ||
            fabs(thin.area - 0.05) > 1e-12) {
            fprintf(stderr, "[slim selftest] valid thin UV rest triangle rejected\n");
            failures++;
        }
    }
    {
        /* A distorted 3x3 planar grid has a descending uniform step, but
         * nonuniform damping makes its first local proposal invalid. With
         * one search trial, that failed proposal used to terminate the solve
         * and incorrectly mark the unchanged map converged. */
        float x[27];
        const int32_t f[24] = {0,1,4, 0,4,3, 1,2,5, 1,5,4,
                               3,4,7, 3,7,6, 4,5,8, 4,8,7};
        const double seed[18] = {
            -.27098144260107127, -.0025446770787633686,
            3.678257275853249, .07525575256522189,
            5.240087834561052, -.11383228570384421,
            -.6988021247795364, .2482505595255208,
            1.9256461089142474, .1753067953715581,
            4.889264325276298, .2417827913140163,
            1.1397988445138199, .6897234438162543,
            1.7980608256682153, .655216866090183,
            5.207217304297799, .6967589494232659};
        double result[18];
        SlimRefineOpts options;
        SlimRefineStats measured;
        Arena_T scratch = Arena_new();
        for (size_t v = 0; v < 9; v++) {
            x[3*v] = 0; x[3*v+1] = (float)(v%3); x[3*v+2] = (float)(v/3);
        }
        SlimRefine_defaults(&options);
        options.iterations = 1; options.pcg_iterations = 500;
        options.pcg_tolerance = 1e-10; options.threads = 1;
        options.strip_pack = 0; options.local_steps = 1;
        options.local_step_min_faces = 1; options.line_search = 1;
        options.energy_tolerance = 0;
        int rc = SlimRefine_run_double(scratch,x,9,f,8,seed,&options,result,&measured);
        if (rc || measured.accepted_steps != 1 || measured.local_component_steps_accepted ||
            measured.flips_after || measured.boundary_intersections_after || measured.converged ||
            !(measured.energy_after < .5*measured.energy_before)) {
            fprintf(stderr,"[slim selftest] local-to-global fallback failed\n");
            failures++;
        }
        Arena_dispose(&scratch);
    }
    if (verts == NULL || faces == NULL || uv == NULL || output == NULL ||
        guide == NULL) {
        free(verts); free(faces); free(uv); free(output); free(guide);
        Arena_dispose(&arena);
        return 1;
    }
    for (y = 0; y < NY; y++) {
        for (x = 0; x < NX; x++) {
            size_t v = y * NX + x;
            verts[v * 3] = 0.0f;
            verts[v * 3 + 1] = (float)y;
            verts[v * 3 + 2] = (float)x;
            uv[v * 2] = 3.0 * (double)x;
            uv[v * 2 + 1] = 0.5 * (double)y;
            guide[v * 2] = (float)x;
            guide[v * 2 + 1] = -(float)y;
        }
    }
    for (y = 0; y + 1 < NY; y++) {
        for (x = 0; x + 1 < NX; x++) {
            int32_t v00 = (int32_t)(y * NX + x);
            int32_t v01 = v00 + 1;
            int32_t v10 = v00 + NX;
            int32_t v11 = v10 + 1;
            faces[face * 3] = v00;
            faces[face * 3 + 1] = v01;
            faces[face * 3 + 2] = v11;
            face++;
            faces[face * 3] = v00;
            faces[face * 3 + 1] = v11;
            faces[face * 3 + 2] = v10;
            face++;
        }
    }
    SlimRefine_defaults(&opts);
    opts.iterations = 10;
    opts.pcg_iterations = 300;
    opts.pcg_tolerance = 1.0e-9;
    opts.energy_tolerance = 1.0e-9;
    opts.strip_pack = 0;
    opts.guide_uv = guide;
    opts.metric_uv = guide;
    opts.guide_warm_start = 1;
    opts.local_steps = 1;
    opts.local_step_min_faces = 1;
    if (SlimRefine_run_double(arena, verts, nv, faces, nf, uv, &opts,
                              output, &stats) != 0) {
        fprintf(stderr, "[slim selftest] solver failed\n");
        failures++;
    } else {
        if (!(stats.energy_after < 0.8 * stats.energy_before)) {
            fprintf(stderr,
                    "[slim selftest] energy did not fall enough: %.9g -> %.9g\n",
                    stats.energy_before, stats.energy_after);
            failures++;
        }
        if (!(stats.qc_max_after < stats.qc_max_before) ||
            stats.flips_after != 0 ||
            stats.boundary_intersections_after != 0) {
            fprintf(stderr,
                    "[slim selftest] invariant/quality failure qc %.6g -> %.6g "
                    "flips=%zu boundary=%zu\n",
                    stats.qc_max_before, stats.qc_max_after, stats.flips_after,
                    stats.boundary_intersections_after);
            failures++;
        }
        if (stats.amg_levels < 2 || stats.amg_coarse_vertices == 0) {
            fprintf(stderr,
                    "[slim selftest] AMG hierarchy missing: levels=%d "
                    "coarsest=%zu\n",
                    stats.amg_levels, stats.amg_coarse_vertices);
            failures++;
        }
        if (stats.guide_components_accepted != 1 ||
            stats.guide_components_rejected != 0 ||
            !(stats.guide_energy_after < stats.energy_before)) {
            fprintf(stderr,
                    "[slim selftest] winding warm start failed: "
                    "accepted/rejected=%zu/%zu E=%.9g->%.9g\n",
                    stats.guide_components_accepted,
                    stats.guide_components_rejected,
                    stats.energy_before, stats.guide_energy_after);
            failures++;
        }
        if (stats.local_component_steps_accepted == 0) {
            fprintf(stderr,
                    "[slim selftest] guarded local step was not exercised\n");
            failures++;
        }
    }
    homotopy_arena = Arena_new();
    SlimRefine_defaults(&opts);
    opts.iterations = 2;
    opts.homotopy_stages = 8;
    opts.pcg_iterations = 300;
    opts.pcg_tolerance = 1.0e-9;
    opts.energy_tolerance = 1.0e-9;
    opts.strip_pack = 0;
    opts.metric_uv = guide;
    if (SlimRefine_run_double(homotopy_arena, verts, nv, faces, nf, uv,
                              &opts, output, &stats) != 0) {
        fprintf(stderr, "[slim selftest] metric homotopy solver failed\n");
        failures++;
    } else if (stats.homotopy_stages_run != 8 ||
               stats.flips_after != 0 ||
               stats.boundary_intersections_after != 0 ||
               stats.reverted ||
               !(stats.energy_after < stats.energy_before)) {
        fprintf(stderr,
                "[slim selftest] metric homotopy invariant failure: "
                "stages=%zu flips=%zu boundary=%zu reverted=%d "
                "E=%.9g->%.9g\n",
                stats.homotopy_stages_run, stats.flips_after,
                stats.boundary_intersections_after, stats.reverted,
                stats.energy_before, stats.energy_after);
        failures++;
    }
    Arena_dispose(&homotopy_arena);
    free(verts); free(faces); free(uv); free(output); free(guide);
    Arena_dispose(&arena);
    fprintf(stderr, "[slim_refine selftest] %s\n",
            failures == 0 ? "ok" : "FAILED");
    return failures;
}
