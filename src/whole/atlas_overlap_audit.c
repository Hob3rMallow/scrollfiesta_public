#include "atlas_overlap_audit.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AOA_EDGE_EPS 1.0e-8
#define AOA_GRID_SCALE 2.0
#define AOA_RECORDS_PER_FACE 32u
#define AOA_GRID_RETRIES 64

typedef struct {
    int64_t u0, u1;
    int64_t v0, v1;
} AoaCellBounds;

typedef struct {
    int64_t cell_u, cell_v;
    int32_t face;
} AoaCellRecord;

typedef struct {
    size_t records;
    size_t cells;
    size_t candidates;
    double cell_size;
} AoaBroadPhaseStats;

static int aoa_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int aoa_compare_i32(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa, b = *(const int32_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int aoa_compare_cell_record(const void *pa, const void *pb)
{
    const AoaCellRecord *a = (const AoaCellRecord *)pa;
    const AoaCellRecord *b = (const AoaCellRecord *)pb;
    if (a->cell_v != b->cell_v)
        return a->cell_v < b->cell_v ? -1 : 1;
    if (a->cell_u != b->cell_u)
        return a->cell_u < b->cell_u ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

static int aoa_compare_observation(const void *pa, const void *pb)
{
    const AtlasOverlapPair *a = (const AtlasOverlapPair *)pa;
    const AtlasOverlapPair *b = (const AtlasOverlapPair *)pb;
    if (a->component0 != b->component0)
        return a->component0 < b->component0 ? -1 : 1;
    if (a->component1 != b->component1)
        return a->component1 < b->component1 ? -1 : 1;
    if (a->turn != b->turn) return a->turn < b->turn ? -1 : 1;
    if (a->face0 != b->face0) return a->face0 < b->face0 ? -1 : 1;
    return a->face1 < b->face1 ? -1 : (a->face1 > b->face1 ? 1 : 0);
}

static double aoa_median(double *value, size_t count)
{
    if (count == 0) return NAN;
    qsort(value, count, sizeof(double), aoa_compare_double);
    return count & 1 ? value[count / 2]
                     : 0.5 * (value[count / 2 - 1] + value[count / 2]);
}

static double aoa_triangle_scale(const double x[3], const double y[3])
{
    double scale = 0.0;
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) % 3;
        double length = hypot(x[j] - x[i], y[j] - y[i]);
        if (length > scale) scale = length;
    }
    return scale;
}

/* Strict convex-triangle intersection.  Boundary contact has zero area and
 * is deliberately not reported.  Testing all six edge normals also catches
 * folded one-rings whose faces share a mesh vertex or edge; topology is not
 * a license for their interiors to overlap.  Coordinates are projected from
 * a local origin so the predicate is stable under a large global U shift. */
static int aoa_triangles_overlap(const double *u, const double *v,
                                 const int32_t *faces,
                                 int32_t fa, int32_t fb)
{
    double ax[3], ay[3], bx[3], by[3];
    for (int i = 0; i < 3; i++) {
        int32_t va = faces[(size_t)fa * 3 + (size_t)i];
        int32_t vb = faces[(size_t)fb * 3 + (size_t)i];
        ax[i] = u[va]; ay[i] = v[va];
        bx[i] = u[vb]; by[i] = v[vb];
    }
    double scale_a = aoa_triangle_scale(ax, ay);
    double scale_b = aoa_triangle_scale(bx, by);
    double scale = scale_a > scale_b ? scale_a : scale_b;
    double overlap_scale = scale_a < scale_b ? scale_a : scale_b;
    if (!(scale > 0.0) || !isfinite(scale)) return 0;
    double area_a = fabs((ax[1] - ax[0]) * (ay[2] - ay[0]) -
                         (ax[2] - ax[0]) * (ay[1] - ay[0]));
    double area_b = fabs((bx[1] - bx[0]) * (by[2] - by[0]) -
                         (bx[2] - bx[0]) * (by[1] - by[0]));
    if (area_a <= AOA_EDGE_EPS * scale_a * scale_a ||
        area_b <= AOA_EDGE_EPS * scale_b * scale_b)
        return 0;

    double origin_x = ax[0], origin_y = ay[0];
    for (int source = 0; source < 2; source++) {
        const double *x = source == 0 ? ax : bx;
        const double *y = source == 0 ? ay : by;
        for (int edge = 0; edge < 3; edge++) {
            int next = (edge + 1) % 3;
            double dx = x[next] - x[edge];
            double dy = y[next] - y[edge];
            double nx = -dy, ny = dx;
            double axis_length = hypot(nx, ny);
            if (!(axis_length > 0.0)) return 0;
            double amin = HUGE_VAL, amax = -HUGE_VAL;
            double bmin = HUGE_VAL, bmax = -HUGE_VAL;
            for (int i = 0; i < 3; i++) {
                double pa = (ax[i] - origin_x) * nx +
                            (ay[i] - origin_y) * ny;
                double pb = (bx[i] - origin_x) * nx +
                            (by[i] - origin_y) * ny;
                if (pa < amin) amin = pa;
                if (pa > amax) amax = pa;
                if (pb < bmin) bmin = pb;
                if (pb > bmax) bmax = pb;
            }
            double lo = amin > bmin ? amin : bmin;
            double hi = amax < bmax ? amax : bmax;
            double tolerance = AOA_EDGE_EPS * axis_length * overlap_scale;
            if (hi - lo <= tolerance) return 0;
        }
    }
    return 1;
}

static int aoa_cell_index(double value, double origin, double cell_size,
                          int64_t *out)
{
    long double scaled = ((long double)value - (long double)origin) /
                         (long double)cell_size;
    if (!isfinite(scaled) || scaled < 0.0L ||
        scaled > (long double)INT64_MAX)
        return -1;
    *out = (int64_t)floorl(scaled);
    return 0;
}

static int aoa_pair_pass(const AoaCellRecord *record, size_t nrecord,
                         const AoaCellBounds *bounds,
                         const int32_t *faces, const double *u,
                         const double *v, int32_t *out_f0,
                         int32_t *out_f1, size_t output_capacity,
                         size_t *out_exact, size_t *out_candidates,
                         size_t *out_cells)
{
    size_t exact = 0, candidates = 0, cells = 0;
    for (size_t first = 0; first < nrecord;) {
        size_t last = first + 1;
        while (last < nrecord &&
               record[last].cell_u == record[first].cell_u &&
               record[last].cell_v == record[first].cell_v)
            last++;
        cells++;
        for (size_t ia = first; ia < last; ia++) {
            for (size_t ib = ia + 1; ib < last; ib++) {
                int32_t fa = record[ia].face;
                int32_t fb = record[ib].face;
                if (fa == fb) continue;
                int64_t canonical_u = bounds[fa].u0 > bounds[fb].u0 ?
                                      bounds[fa].u0 : bounds[fb].u0;
                int64_t canonical_v = bounds[fa].v0 > bounds[fb].v0 ?
                                      bounds[fa].v0 : bounds[fb].v0;
                if (record[first].cell_u != canonical_u ||
                    record[first].cell_v != canonical_v)
                    continue;
                if (candidates == SIZE_MAX) return -1;
                candidates++;
                if (!aoa_triangles_overlap(u, v, faces, fa, fb)) continue;
                if (exact == SIZE_MAX) return -1;
                if (out_f0 != NULL) {
                    if (exact >= output_capacity) return -1;
                    if (fa > fb) {
                        int32_t temporary = fa; fa = fb; fb = temporary;
                    }
                    out_f0[exact] = fa;
                    out_f1[exact] = fb;
                }
                exact++;
            }
        }
        first = last;
    }
    *out_exact = exact;
    if (out_candidates != NULL) *out_candidates = candidates;
    if (out_cells != NULL) *out_cells = cells;
    return 0;
}

static int aoa_detect_pairs(Arena_T arena, const int32_t *faces,
                            size_t nf, size_t nv,
                            const double *u, const double *v,
                            int32_t **out_f0, int32_t **out_f1,
                            size_t *out_count, int *out_truncated,
                            AoaBroadPhaseStats *stats)
{
    double *area = (double *)ARENA_ALLOC(arena, nf * sizeof(double));
    double umin = HUGE_VAL, umax = -HUGE_VAL;
    double vmin = HUGE_VAL, vmax = -HUGE_VAL;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f * 3], b = faces[f * 3 + 1];
        int32_t c = faces[f * 3 + 2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv)
            return -1;
        double ax = u[a], ay = v[a];
        double bx = u[b], by = v[b];
        double cx = u[c], cy = v[c];
        area[f] = 0.5 * fabs((bx - ax) * (cy - ay) -
                             (cx - ax) * (by - ay));
        for (int k = 0; k < 3; k++) {
            int32_t vertex = faces[f * 3 + (size_t)k];
            if (u[vertex] < umin) umin = u[vertex];
            if (u[vertex] > umax) umax = u[vertex];
            if (v[vertex] < vmin) vmin = v[vertex];
            if (v[vertex] > vmax) vmax = v[vertex];
        }
    }
    qsort(area, nf, sizeof(double), aoa_compare_double);
    double median_area = area[nf / 2];
    if (!(median_area >= 1.0e-12)) median_area = 1.0e-6;
    double cell_size = AOA_GRID_SCALE * sqrt(median_area);
    if (!(cell_size >= 1.0e-6)) cell_size = 1.0e-6;

    if (nf > SIZE_MAX / AOA_RECORDS_PER_FACE) return -1;
    size_t record_target = nf * AOA_RECORDS_PER_FACE;
    AoaCellBounds *bounds = (AoaCellBounds *)ARENA_ALLOC(
        arena, nf * sizeof(AoaCellBounds));
    size_t nrecord = 0;
    int grid_ready = 0;
    for (int attempt = 0; attempt < AOA_GRID_RETRIES; attempt++) {
        long double record_estimate = 0.0L;
        int index_overflow = 0;
        for (size_t f = 0; f < nf; f++) {
            double face_umin = HUGE_VAL, face_umax = -HUGE_VAL;
            double face_vmin = HUGE_VAL, face_vmax = -HUGE_VAL;
            for (int k = 0; k < 3; k++) {
                int32_t vertex = faces[f * 3 + (size_t)k];
                if (u[vertex] < face_umin) face_umin = u[vertex];
                if (u[vertex] > face_umax) face_umax = u[vertex];
                if (v[vertex] < face_vmin) face_vmin = v[vertex];
                if (v[vertex] > face_vmax) face_vmax = v[vertex];
            }
            if (aoa_cell_index(face_umin, umin, cell_size,
                               &bounds[f].u0) != 0 ||
                aoa_cell_index(face_umax, umin, cell_size,
                               &bounds[f].u1) != 0 ||
                aoa_cell_index(face_vmin, vmin, cell_size,
                               &bounds[f].v0) != 0 ||
                aoa_cell_index(face_vmax, vmin, cell_size,
                               &bounds[f].v1) != 0) {
                index_overflow = 1;
                break;
            }
            long double width = (long double)bounds[f].u1 -
                                (long double)bounds[f].u0 + 1.0L;
            long double height = (long double)bounds[f].v1 -
                                 (long double)bounds[f].v0 + 1.0L;
            record_estimate += width * height;
        }
        if (!index_overflow &&
            record_estimate <= (long double)record_target) {
            nrecord = (size_t)record_estimate;
            grid_ready = 1;
            break;
        }

        long double growth = 2.0L;
        if (index_overflow) {
            long double span_u = (long double)umax - (long double)umin;
            long double span_v = (long double)vmax - (long double)vmin;
            long double span = span_u > span_v ? span_u : span_v;
            long double required = span / ((long double)INT64_MAX / 4.0L);
            if (required > (long double)cell_size)
                growth = 2.0L * required / (long double)cell_size;
        } else if (record_estimate > 0.0L) {
            growth = sqrtl(record_estimate /
                           (long double)record_target) * 1.05L;
            if (growth < 2.0L) growth = 2.0L;
        }
        long double next_size = (long double)cell_size * growth;
        if (!isfinite(next_size) || next_size > (long double)HUGE_VAL)
            return -1;
        cell_size = (double)next_size;
        if (!isfinite(cell_size) || !(cell_size > 0.0)) return -1;
    }
    if (!grid_ready || nrecord > SIZE_MAX / sizeof(AoaCellRecord)) return -1;

    AoaCellRecord *record = (AoaCellRecord *)ARENA_ALLOC(
        arena, (nrecord ? nrecord : 1) * sizeof(AoaCellRecord));
    size_t cursor = 0;
    for (size_t f = 0; f < nf; f++) {
        int64_t cell_v = bounds[f].v0;
        for (;;) {
            int64_t cell_u = bounds[f].u0;
            for (;;) {
                if (cursor >= nrecord) return -1;
                record[cursor].cell_u = cell_u;
                record[cursor].cell_v = cell_v;
                record[cursor].face = (int32_t)f;
                cursor++;
                if (cell_u == bounds[f].u1) break;
                cell_u++;
            }
            if (cell_v == bounds[f].v1) break;
            cell_v++;
        }
    }
    if (cursor != nrecord) return -1;
    qsort(record, nrecord, sizeof(AoaCellRecord), aoa_compare_cell_record);

    size_t exact = 0, candidates = 0, occupied_cells = 0;
    if (aoa_pair_pass(record, nrecord, bounds, faces, u, v,
                      NULL, NULL, 0, &exact, &candidates,
                      &occupied_cells) != 0)
        return -1;
    if (exact > SIZE_MAX / sizeof(int32_t)) return -1;
    int32_t *f0 = NULL, *f1 = NULL;
    if (exact > 0) {
        f0 = (int32_t *)ARENA_ALLOC(arena, exact * sizeof(int32_t));
        f1 = (int32_t *)ARENA_ALLOC(arena, exact * sizeof(int32_t));
        size_t verified_exact = 0;
        if (aoa_pair_pass(record, nrecord, bounds, faces, u, v,
                          f0, f1, exact, &verified_exact, NULL, NULL) != 0 ||
            verified_exact != exact)
            return -1;
    }
    *out_f0 = f0;
    *out_f1 = f1;
    *out_count = exact;
    *out_truncated = 0;
    stats->records = nrecord;
    stats->cells = occupied_cells;
    stats->candidates = candidates;
    stats->cell_size = cell_size;
    return 0;
}

static int32_t aoa_face_component(const int32_t *faces, size_t f,
                                  const int32_t *vertex_component,
                                  size_t ncomponents)
{
    int32_t c = vertex_component[faces[f * 3]];
    if (c < 0 || (size_t)c >= ncomponents) return -1;
    return vertex_component[faces[f * 3 + 1]] == c &&
           vertex_component[faces[f * 3 + 2]] == c ? c : -1;
}

static double aoa_face_mean_float(const int32_t *faces, size_t f,
                                  const float *value)
{
    return ((double)value[faces[f * 3]] +
            (double)value[faces[f * 3 + 1]] +
            (double)value[faces[f * 3 + 2]]) / 3.0;
}

static double aoa_face_mean_double(const int32_t *faces, size_t f,
                                   const double *value)
{
    return (value[faces[f * 3]] + value[faces[f * 3 + 1]] +
            value[faces[f * 3 + 2]]) / 3.0;
}

int AtlasOverlapAudit_build(Arena_T arena,
                            const int32_t *faces, size_t nfaces,
                            size_t nvertices,
                            const double *u, const double *v,
                            const float *registered_u,
                            const float *phi,
                            const int32_t *vertex_component,
                            size_t ncomponents,
                            AtlasOverlapAudit *out)
{
    if (arena == NULL || faces == NULL || nfaces == 0 ||
        nfaces > (size_t)INT32_MAX || nvertices < 3 || u == NULL ||
        v == NULL || registered_u == NULL || phi == NULL ||
        vertex_component == NULL || ncomponents == 0 || out == NULL)
        return -1;
    memset(out, 0, sizeof *out);
    for (size_t i = 0; i < nvertices; i++) {
        if (!isfinite(u[i]) || !isfinite(v[i]) ||
            !isfinite((double)registered_u[i]) ||
            !isfinite((double)phi[i]))
            return -1;
    }
    size_t *component_faces = (size_t *)ARENA_CALLOC(
        arena, ncomponents, sizeof(size_t));
    for (size_t f = 0; f < nfaces; f++) {
        int32_t c = aoa_face_component(
            faces, f, vertex_component, ncomponents);
        if (c < 0) return -1;
        component_faces[c]++;
    }
    int32_t *pair0 = NULL, *pair1 = NULL;
    size_t npair = 0;
    AoaBroadPhaseStats broad_phase;
    memset(&broad_phase, 0, sizeof broad_phase);
    if (aoa_detect_pairs(arena, faces, nfaces, nvertices, u, v,
                         &pair0, &pair1, &npair,
                         &out->pair_buffer_truncated,
                         &broad_phase) != 0)
        return -1;
    out->indexed_faces = nfaces;
    out->broad_phase_records = broad_phase.records;
    out->broad_phase_cells = broad_phase.cells;
    out->broad_phase_candidate_pairs = broad_phase.candidates;
    out->broad_phase_cell_size = broad_phase.cell_size;
    out->broad_phase_complete = 1;
    out->exact_face_pairs = npair;
    AtlasOverlapPair *obs = (AtlasOverlapPair *)ARENA_ALLOC(
        arena, (npair ? npair : 1) * sizeof(AtlasOverlapPair));
    const double two_pi = 6.283185307179586476925286766559;
    for (size_t i = 0; i < npair; i++) {
        int32_t f0 = pair0[i], f1 = pair1[i];
        int32_t c0 = aoa_face_component(
            faces, (size_t)f0, vertex_component, ncomponents);
        int32_t c1 = aoa_face_component(
            faces, (size_t)f1, vertex_component, ncomponents);
        double p0 = aoa_face_mean_float(faces, (size_t)f0, phi);
        double p1 = aoa_face_mean_float(faces, (size_t)f1, phi);
        double q0 = aoa_face_mean_double(faces, (size_t)f0, u);
        double q1 = aoa_face_mean_double(faces, (size_t)f1, u);
        double r0 = aoa_face_mean_float(faces, (size_t)f0, registered_u);
        double r1 = aoa_face_mean_float(faces, (size_t)f1, registered_u);
        if (c0 > c1) {
            int32_t ct = c0; c0 = c1; c1 = ct;
            double t = p0; p0 = p1; p1 = t;
            t = q0; q0 = q1; q1 = t;
            t = r0; r0 = r1; r1 = t;
            int32_t ft = f0; f0 = f1; f1 = ft;
        }
        double turns = (p1 - p0) / two_pi;
        long turn = lround(turns);
        if (turn < INT32_MIN || turn > INT32_MAX) return -1;
        obs[i].component0 = c0;
        obs[i].component1 = c1;
        obs[i].face0 = f0;
        obs[i].face1 = f1;
        obs[i].bundle = -1;
        obs[i].turn = (int32_t)turn;
        obs[i].phase_residual = p1 - p0 - two_pi * (double)turn;
        obs[i].parameter_du = q1 - q0;
        obs[i].registered_du = r1 - r0;
        obs[i].registered_parameter_shift =
            (r1 - q1) - (r0 - q0);
        if (c0 == c1) out->same_component_pairs++;
        else out->cross_component_pairs++;
    }
    qsort(obs, npair, sizeof(AtlasOverlapPair), aoa_compare_observation);
    AtlasOverlapBundle *bundle = (AtlasOverlapBundle *)ARENA_ALLOC(
        arena, (npair ? npair : 1) * sizeof(AtlasOverlapBundle));
    double *scratch = (double *)ARENA_ALLOC(
        arena, (npair ? npair : 1) * sizeof(double));
    int32_t *face_scratch = (int32_t *)ARENA_ALLOC(
        arena, (npair ? npair : 1) * sizeof(int32_t));
    size_t nbundle = 0;
    for (size_t first = 0; first < npair;) {
        size_t last = first + 1;
        while (last < npair &&
               obs[last].component0 == obs[first].component0 &&
               obs[last].component1 == obs[first].component1)
            last++;
        size_t bundle_index = nbundle++;
        AtlasOverlapBundle *b = &bundle[bundle_index];
        memset(b, 0, sizeof *b);
        b->component0 = obs[first].component0;
        b->component1 = obs[first].component1;
        b->face_pairs = last - first;
        b->component_faces0 = component_faces[b->component0];
        b->component_faces1 = component_faces[b->component1];
        for (size_t i = first; i < last; i++)
            obs[i].bundle = (int32_t)bundle_index;
        size_t mode_count = 0;
        for (size_t r0 = first; r0 < last;) {
            size_t r1 = r0 + 1;
            while (r1 < last && obs[r1].turn == obs[r0].turn) r1++;
            if (r1 - r0 > mode_count) {
                mode_count = r1 - r0;
                b->turn_mode = obs[r0].turn;
            }
            r0 = r1;
        }
        b->turn_agreement = (double)mode_count / (double)b->face_pairs;
        size_t count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = obs[i].phase_residual;
        }
        b->phase_residual_median = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = fabs(obs[i].phase_residual -
                                      b->phase_residual_median);
        }
        b->phase_residual_mad = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = obs[i].parameter_du;
        }
        b->parameter_du_median = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = fabs(obs[i].parameter_du -
                                      b->parameter_du_median);
        }
        b->parameter_du_mad = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = obs[i].registered_du;
        }
        b->registered_du_median = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = fabs(obs[i].registered_du -
                                      b->registered_du_median);
        }
        b->registered_du_mad = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = obs[i].registered_parameter_shift;
        }
        b->registered_parameter_shift_median = aoa_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++) {
            if (obs[i].turn != b->turn_mode) continue;
            scratch[count++] = fabs(
                obs[i].registered_parameter_shift -
                b->registered_parameter_shift_median);
        }
        b->registered_parameter_shift_mad = aoa_median(scratch, count);

        count = 0;
        for (size_t i = first; i < last; i++)
            if (obs[i].turn == b->turn_mode)
                face_scratch[count++] = obs[i].face0;
        qsort(face_scratch, count, sizeof(int32_t), aoa_compare_i32);
        for (size_t i = 0; i < count; i++)
            if (i == 0 || face_scratch[i] != face_scratch[i - 1])
                b->unique_faces0++;
        count = 0;
        for (size_t i = first; i < last; i++)
            if (obs[i].turn == b->turn_mode)
                face_scratch[count++] = obs[i].face1;
        qsort(face_scratch, count, sizeof(int32_t), aoa_compare_i32);
        for (size_t i = 0; i < count; i++)
            if (i == 0 || face_scratch[i] != face_scratch[i - 1])
                b->unique_faces1++;
        first = last;
    }
    out->bundles = bundle;
    out->nbundles = nbundle;
    out->pairs = obs;
    out->npairs = npair;
    return 0;
}

static void aoa_check(int condition, const char *what, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[atlas_overlap_audit selftest] FAIL: %s\n", what);
    (*fails)++;
}

int AtlasOverlapAudit_selftest(void)
{
    Arena_T arena = Arena_new();
    int fails = 0;
    double u[8] = {0, 4, 4, 0, 0, 4, 4, 0};
    double v[8] = {0, 0, 4, 4, 0, 0, 4, 4};
    float reg[8] = {0, 4, 4, 0, 100, 104, 104, 100};
    float phi[8] = {0, 0, 0, 0, 6.283185307f, 6.283185307f,
                    6.283185307f, 6.283185307f};
    int32_t component[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    int32_t face[12] = {0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7};
    AtlasOverlapAudit audit;
    int rc = AtlasOverlapAudit_build(
        arena, face, 4, 8, u, v, reg, phi, component, 2, &audit);
    aoa_check(rc == 0, "build", &fails);
    aoa_check(audit.cross_component_pairs > 0 &&
              audit.same_component_pairs == 0,
              "cross-component overlap only", &fails);
    aoa_check(audit.nbundles == 1 &&
              audit.bundles[0].turn_mode == 1 &&
              audit.bundles[0].turn_agreement > 0.999 &&
              fabs(audit.bundles[0].registered_du_median - 100.0) < 1e-6 &&
              fabs(audit.bundles[0].registered_parameter_shift_median -
                   100.0) < 1e-6,
              "coherent winding and registered displacement", &fails);
    aoa_check(audit.npairs == audit.exact_face_pairs &&
              audit.npairs > 0 && audit.pairs != NULL &&
              audit.pairs[0].bundle == 0 &&
              audit.pairs[0].face0 >= 0 && audit.pairs[0].face1 >= 0,
              "exact face pairs retained", &fails);
    aoa_check(audit.broad_phase_complete &&
              audit.indexed_faces == 4 &&
              audit.broad_phase_records >= audit.indexed_faces &&
              !audit.pair_buffer_truncated,
              "broad phase reports complete coverage", &fails);

    {
        double wide_u[9] = {
            0, 1, 0,
            100000, 100004, 100000,
            100000, 100004, 100000
        };
        double shifted_u[9];
        double wide_v[9] = {
            20, 20, 21,
            0, 0, 4,
            0, 0, 4
        };
        float wide_reg[9] = {0, 1, 0, 0, 4, 0, 100, 104, 100};
        float wide_phi[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        int32_t wide_component[9] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
        int32_t wide_face[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        AtlasOverlapAudit wide, shifted;
        int wide_rc = AtlasOverlapAudit_build(
            arena, wide_face, 3, 9, wide_u, wide_v, wide_reg, wide_phi,
            wide_component, 3, &wide);
        for (int i = 0; i < 9; i++) shifted_u[i] = wide_u[i] + 1.0e9;
        int shifted_rc = AtlasOverlapAudit_build(
            arena, wide_face, 3, 9, shifted_u, wide_v, wide_reg, wide_phi,
            wide_component, 3, &shifted);
        aoa_check(wide_rc == 0 && shifted_rc == 0 &&
                  wide.broad_phase_complete && shifted.broad_phase_complete &&
                  wide.indexed_faces == 3 && shifted.indexed_faces == 3,
                  "large-span audits cover every face", &fails);
        aoa_check(wide.npairs == 1 && shifted.npairs == wide.npairs,
                  "large-span and translated overlap invariance", &fails);
    }

    {
        double fold_u[6] = {0, 4, 0, 3, 0.5, 0};
        double fold_v[6] = {0, 0, 4, 0.5, 3, -4};
        float fold_reg[6] = {0, 4, 0, 3, 0.5f, 0};
        float fold_phi[6] = {0, 0, 0, 0, 0, 0};
        int32_t fold_component[6] = {0, 0, 0, 0, 0, 0};
        int32_t fold_face[9] = {0, 1, 2, 0, 3, 4, 1, 0, 5};
        AtlasOverlapAudit fold;
        int fold_rc = AtlasOverlapAudit_build(
            arena, fold_face, 3, 6, fold_u, fold_v, fold_reg, fold_phi,
            fold_component, 1, &fold);
        aoa_check(fold_rc == 0 && fold.npairs == 1 &&
                  fold.same_component_pairs == 1 &&
                  fold.pairs[0].face0 == 0 && fold.pairs[0].face1 == 1,
                  "shared-vertex fold counts but boundary contact does not",
                  &fails);
    }

    {
        enum { dense_faces = 32, dense_vertices = dense_faces * 3 };
        double dense_u[dense_vertices], dense_v[dense_vertices];
        float dense_reg[dense_vertices], dense_phi[dense_vertices];
        int32_t dense_component[dense_vertices];
        int32_t dense_face[dense_faces * 3];
        for (int f = 0; f < dense_faces; f++) {
            const double triangle_u[3] = {0, 4, 0};
            const double triangle_v[3] = {0, 0, 4};
            for (int k = 0; k < 3; k++) {
                int vertex = f * 3 + k;
                dense_u[vertex] = triangle_u[k];
                dense_v[vertex] = triangle_v[k];
                dense_reg[vertex] = (float)triangle_u[k];
                dense_phi[vertex] = 0.0f;
                dense_component[vertex] = f;
                dense_face[vertex] = vertex;
            }
        }
        AtlasOverlapAudit dense;
        int dense_rc = AtlasOverlapAudit_build(
            arena, dense_face, dense_faces, dense_vertices,
            dense_u, dense_v, dense_reg, dense_phi, dense_component,
            dense_faces, &dense);
        size_t expected = (size_t)dense_faces * (dense_faces - 1) / 2;
        aoa_check(dense_rc == 0 && dense.npairs == expected &&
                  dense.npairs > (size_t)dense_faces * 10 &&
                  !dense.pair_buffer_truncated &&
                  dense.broad_phase_candidate_pairs == expected,
                  "overlap output is not capped at ten pairs per face",
                  &fails);
    }
    fprintf(stderr, "[atlas_overlap_audit selftest] %s\n",
            fails == 0 ? "PASSED" : "FAILED");
    Arena_dispose(&arena);
    return fails;
}
