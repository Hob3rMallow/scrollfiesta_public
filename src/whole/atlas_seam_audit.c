#include "atlas_seam_audit.h"

#include "../common/kdtree.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ASA_R3D 3.5
#define ASA_SKIN 4.0
#define ASA_CHUNK 128.0
#define ASA_PHASE_GATE 0.15

static int asa_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int asa_compare_observation(const void *pa, const void *pb)
{
    const AtlasSeamPair *a = (const AtlasSeamPair *)pa;
    const AtlasSeamPair *b = (const AtlasSeamPair *)pb;
    if (a->component0 != b->component0)
        return a->component0 < b->component0 ? -1 : 1;
    if (a->component1 != b->component1)
        return a->component1 < b->component1 ? -1 : 1;
    return a->target_shift < b->target_shift ? -1 :
           (a->target_shift > b->target_shift ? 1 : 0);
}

static double asa_median(double *value, size_t count)
{
    if (count == 0) return NAN;
    qsort(value, count, sizeof(double), asa_compare_double);
    return count & 1u ? value[count / 2]
                      : 0.5 * (value[count / 2 - 1] + value[count / 2]);
}

static size_t asa_vertex_cube(const PieceSet *ps, size_t vertex)
{
    size_t lo = 0, hi = ps->n_cubes;
    while (lo + 1 < hi) {
        size_t mid = (lo + hi) / 2;
        if (ps->cube_voff[mid] <= vertex) lo = mid;
        else hi = mid;
    }
    return lo;
}

int AtlasSeamAudit_build(Arena_T arena, const PieceSet *ps,
                         const double *parameter_u,
                         const int32_t *vertex_component,
                         size_t ncomponents,
                         AtlasSeamAudit *out)
{
    if (arena == NULL || ps == NULL || parameter_u == NULL ||
        vertex_component == NULL || ncomponents == 0 ||
        ncomponents > (size_t)INT32_MAX || out == NULL || ps->nv == 0 ||
        ps->nf == 0 || ps->n_cubes < 2 || ps->verts == NULL ||
        ps->uv == NULL || ps->phi == NULL || ps->gid == NULL ||
        ps->faces == NULL || ps->cube_voff == NULL || ps->cube_org == NULL)
        return -1;
    memset(out, 0, sizeof *out);
    uint8_t *referenced = (uint8_t *)ARENA_CALLOC(
        arena, ps->nv, sizeof(uint8_t));
    for (size_t f = 0; f < ps->nf; f++)
        for (int k = 0; k < 3; k++) {
            int32_t v = ps->faces[f * 3 + (size_t)k];
            if (v < 0 || (size_t)v >= ps->nv) return -1;
            referenced[v] = 1;
        }
    int32_t *boundary = (int32_t *)ARENA_ALLOC(
        arena, ps->nv * sizeof(int32_t));
    size_t nb = 0;
    for (size_t c = 0; c < ps->n_cubes; c++) {
        for (size_t i = ps->cube_voff[c]; i < ps->cube_voff[c + 1]; i++) {
            if (!referenced[i] || ps->gid[i] < 0 ||
                vertex_component[i] < 0 ||
                (size_t)vertex_component[i] >= ncomponents ||
                !isfinite(parameter_u[i]))
                continue;
            double dmin = HUGE_VAL;
            for (int axis = 0; axis < 3; axis++) {
                double t = (double)ps->verts[i * 3 + (size_t)axis] -
                           (double)ps->cube_org[c][axis];
                double d = fmin(fabs(t), fabs(ASA_CHUNK - t));
                if (d < dmin) dmin = d;
            }
            if (dmin < ASA_SKIN) boundary[nb++] = (int32_t)i;
        }
    }
    out->boundary_vertices = nb;
    if (nb < 2) return 0;
    float *point = (float *)ARENA_ALLOC(arena, nb * 3 * sizeof(float));
    for (size_t i = 0; i < nb; i++)
        for (int axis = 0; axis < 3; axis++)
            point[i * 3 + (size_t)axis] =
                ps->verts[(size_t)boundary[i] * 3 + (size_t)axis];
    KDTree_T kd = KDTree_new(arena, point, nb);
    int32_t *best = (int32_t *)ARENA_ALLOC(
        arena, ps->nv * sizeof(int32_t));
    float *best_distance = (float *)ARENA_ALLOC(
        arena, ps->nv * sizeof(float));
    for (size_t i = 0; i < nb; i++) {
        best[boundary[i]] = -1;
        best_distance[boundary[i]] = HUGE_VALF;
    }
    const double two_pi = 6.283185307179586476925286766559;
    int32_t hit[64];
    for (size_t i = 0; i < nb; i++) {
        int32_t gi = boundary[i];
        size_t ci = asa_vertex_cube(ps, (size_t)gi);
        float query[3] = {point[i * 3], point[i * 3 + 1],
                          point[i * 3 + 2]};
        size_t nh = KDTree_ball_query(kd, query, ASA_R3D * ASA_R3D,
                                      hit, 64);
        for (size_t h = 0; h < nh; h++) {
            int32_t gj = boundary[hit[h]];
            if (gj == gi || asa_vertex_cube(ps, (size_t)gj) == ci) continue;
            double phase = remainder((double)ps->phi[gj] -
                                     (double)ps->phi[gi], two_pi);
            if (fabs(phase) > ASA_PHASE_GATE) continue;
            double d2 = 0.0;
            for (int axis = 0; axis < 3; axis++) {
                double d = (double)ps->verts[(size_t)gj * 3 + (size_t)axis] -
                           (double)ps->verts[(size_t)gi * 3 + (size_t)axis];
                d2 += d * d;
            }
            if ((float)d2 < best_distance[gi]) {
                best_distance[gi] = (float)d2;
                best[gi] = gj;
            }
        }
    }
    AtlasSeamPair *observation = (AtlasSeamPair *)ARENA_ALLOC(
        arena, nb * sizeof(AtlasSeamPair));
    size_t no = 0;
    for (size_t i = 0; i < nb; i++) {
        int32_t gi = boundary[i], gj = best[gi];
        if (gj < 0 || gj < gi || best[gj] != gi) continue;
        int32_t v0 = gi, v1 = gj;
        int32_t c0 = vertex_component[gi], c1 = vertex_component[gj];
        double q0 = parameter_u[gi], q1 = parameter_u[gj];
        double r0 = ps->uv[(size_t)gi * 2], r1 = ps->uv[(size_t)gj * 2];
        double p0 = ps->phi[gi], p1 = ps->phi[gj];
        if (c0 > c1) {
            int32_t ct = c0; c0 = c1; c1 = ct;
            int32_t vt = v0; v0 = v1; v1 = vt;
            double t = q0; q0 = q1; q1 = t;
            t = r0; r0 = r1; r1 = t;
            t = p0; p0 = p1; p1 = t;
        }
        AtlasSeamPair *o = &observation[no++];
        o->vertex0 = v0;
        o->vertex1 = v1;
        o->component0 = c0;
        o->component1 = c1;
        o->bundle = -1;
        o->target_shift = q0 - q1;
        o->distance = sqrt((double)best_distance[gi]);
        o->phase_residual = remainder(p1 - p0, two_pi);
        o->registered_du = r1 - r0;
        out->mutual_pairs++;
        if (c0 == c1) out->same_component_pairs++;
        else out->cross_component_pairs++;
    }
    if (no == 0) return 0;
    qsort(observation, no, sizeof(AtlasSeamPair), asa_compare_observation);
    AtlasSeamBundle *bundle = (AtlasSeamBundle *)ARENA_ALLOC(
        arena, no * sizeof(AtlasSeamBundle));
    double *scratch = (double *)ARENA_ALLOC(arena, no * sizeof(double));
    size_t nbundle = 0;
    for (size_t first = 0; first < no;) {
        size_t last = first + 1;
        while (last < no &&
               observation[last].component0 == observation[first].component0 &&
               observation[last].component1 == observation[first].component1)
            last++;
        size_t bundle_index = nbundle++;
        AtlasSeamBundle *b = &bundle[bundle_index];
        memset(b, 0, sizeof *b);
        b->component0 = observation[first].component0;
        b->component1 = observation[first].component1;
        b->pairs = last - first;
        for (size_t j = first; j < last; j++)
            observation[j].bundle = (int32_t)bundle_index;
#define ASA_AGGREGATE(field, median_name, mad_name) do {                    \
            size_t count = 0;                                               \
            for (size_t j = first; j < last; j++)                           \
                scratch[count++] = observation[j].field;                    \
            b->median_name = asa_median(scratch, count);                    \
            count = 0;                                                      \
            for (size_t j = first; j < last; j++)                           \
                scratch[count++] = fabs(observation[j].field -              \
                                        b->median_name);                    \
            b->mad_name = asa_median(scratch, count);                       \
        } while (0)
        ASA_AGGREGATE(target_shift, target_shift_median, target_shift_mad);
        ASA_AGGREGATE(phase_residual, phase_residual_median,
                      phase_residual_mad);
        ASA_AGGREGATE(registered_du, registered_du_median, registered_du_mad);
#undef ASA_AGGREGATE
        size_t count = 0;
        for (size_t j = first; j < last; j++)
            scratch[count++] = observation[j].distance;
        b->distance_median = asa_median(scratch, count);
        first = last;
    }
    out->bundles = bundle;
    out->nbundles = nbundle;
    out->pairs = observation;
    out->npairs = no;
    return 0;
}

static void asa_check(int condition, const char *what, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[atlas_seam_audit selftest] FAIL: %s\n", what);
    (*fails)++;
}

int AtlasSeamAudit_selftest(void)
{
    Arena_T arena = Arena_new();
    PieceSet ps;
    memset(&ps, 0, sizeof ps);
    float verts[18] = {0, 0, 128, 1, 0, 128, 0, 1, 128,
                       0, 0, 128, 1, 0, 128, 0, 1, 128};
    float uv[12] = {0, 0, 1, 0, 0, 1, 20, 0, 21, 0, 20, 1};
    float phi[6] = {0, 0, 0, 0.01f, 0.01f, 0.01f};
    int32_t gid[6] = {0, 0, 0, 0, 0, 0};
    int32_t faces[6] = {0, 1, 2, 3, 4, 5};
    size_t offset[3] = {0, 3, 6};
    long origin[2][3] = {{0, 0, 0}, {0, 0, 128}};
    double parameter[6] = {0, 1, 0, 20, 21, 20};
    int32_t component[6] = {0, 0, 0, 1, 1, 1};
    ps.verts = verts; ps.uv = uv; ps.phi = phi; ps.gid = gid;
    ps.faces = faces; ps.nv = 6; ps.nf = 2; ps.n_cubes = 2;
    ps.cube_voff = offset; ps.cube_org = origin;
    AtlasSeamAudit audit;
    int fails = 0;
    int rc = AtlasSeamAudit_build(
        arena, &ps, parameter, component, 2, &audit);
    asa_check(rc == 0, "build", &fails);
    asa_check(audit.mutual_pairs == 3 && audit.nbundles == 1,
              "mutual seam bundle", &fails);
    asa_check(audit.npairs == 3 && audit.pairs != NULL &&
              audit.pairs[0].bundle == 0 &&
              audit.pairs[0].vertex0 >= 0 && audit.pairs[0].vertex0 < 3 &&
              audit.pairs[0].vertex1 >= 3 && audit.pairs[0].vertex1 < 6,
              "pointwise seam rows retained", &fails);
    asa_check(audit.nbundles == 1 &&
              fabs(audit.bundles[0].target_shift_median + 20.0) < 1e-8 &&
              audit.bundles[0].target_shift_mad < 1e-8,
              "constant closure target", &fails);
    fprintf(stderr, "[atlas_seam_audit selftest] %s\n",
            fails == 0 ? "PASSED" : "FAILED");
    Arena_dispose(&arena);
    return fails;
}
