#include "atlas_ribbon_fit.h"

#include "monotone_qp.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

typedef struct {
    int32_t chart;
    int32_t face;
    double p[2];
    double tangent[2];
    uint32_t count;
} ArfChartSample;

typedef struct {
    double p[2];
    double tangent[2];
    uint32_t charts;
    int32_t chart;
    int32_t face;
    uint64_t rank;
} ArfCluster;

typedef struct {
    int n;
    int bandwidth;
    double *lower;
} ArfBand;

typedef struct {
    int32_t index;
    double min_x, max_x;
    double min_y, max_y;
    double a[2], b[2];
} ArfSegmentBox;

typedef struct {
    size_t sample;
    uint64_t rank;
    double u;
} ArfRegisterRef;

typedef struct {
    double *lower;
    uint8_t *kind;
    uint32_t *updates;
    uint64_t *sample_lo;
    uint64_t *sample_hi;
    double *required_arc;
    double *endpoint_distance;
    size_t ncharts;
} ArfRegisterBounds;

static double arf_dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void arf_cross3(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

static double arf_normalize2(double p[2])
{
    double n = sqrt(p[0] * p[0] + p[1] * p[1]);
    if (n > 1.0e-30) {
        p[0] /= n;
        p[1] /= n;
    }
    return n;
}

static double arf_normalize3(double p[3])
{
    double n = sqrt(arf_dot3(p, p));
    if (n > 1.0e-30) {
        p[0] /= n;
        p[1] /= n;
        p[2] /= n;
    }
    return n;
}

static double arf_distance2(const double a[2], const double b[2])
{
    double x = a[0] - b[0];
    double y = a[1] - b[1];
    return sqrt(x * x + y * y);
}

static double arf_mod_tau(double angle);
static int arf_registration_winding_direction(
    const AtlasRibbonObservationSet *set);

static void arf_make_frame(const ScaffoldCalib *cal,
                           double axis[3], double basis0[3], double basis1[3])
{
    axis[0] = (double)cal->axis_dir[0];
    axis[1] = (double)cal->axis_dir[1];
    axis[2] = (double)cal->axis_dir[2];
    if (arf_normalize3(axis) < 1.0e-12) {
        axis[0] = 1.0;
        axis[1] = axis[2] = 0.0;
    }
    double reference[3] = {1.0, 0.0, 0.0};
    if (fabs(axis[0]) > 0.8) {
        reference[0] = 0.0;
        reference[1] = 1.0;
    }
    arf_cross3(axis, reference, basis0);
    arf_normalize3(basis0);
    arf_cross3(axis, basis0, basis1);
    arf_normalize3(basis1);
}

static void arf_project_point(const AtlasRibbonObservationSet *set,
                              const float p[3], double out[2])
{
    double d[3] = {
        (double)p[0] - set->axis_point[0],
        (double)p[1] - set->axis_point[1],
        (double)p[2] - set->axis_point[2]
    };
    out[0] = arf_dot3(d, set->basis0);
    out[1] = arf_dot3(d, set->basis1);
}

void AtlasRibbonFitOptions_default(AtlasRibbonFitOptions *opts)
{
    if (opts == NULL) return;
    memset(opts, 0, sizeof *opts);
    opts->slice_spacing = 4.0;
    opts->observation_u_spacing = 2.0;
    opts->fit_u_spacing = 16.0;
    opts->local_xyz_tolerance = 8.0;
    opts->tangent_dot_min = 0.5;
    opts->max_fill_u = 512.0;
    opts->max_bridge_stretch = 1.25;
    opts->bridge_slack = 4.0;
    opts->lambda_position = 1.0;
    opts->lambda_smooth = 8.0;
    opts->lambda_tangent = 4.0;
    opts->lambda_register_vertical = 16.0;
    opts->register_sweeps = 6;
    opts->collision_relaxation = 0.35;
    opts->collision_rounds = 8;
    opts->collision_sweeps = 4;
    opts->collision_polish_rounds = 0;
    opts->mode = ATLAS_RIBBON_TANGENT;
}

static int arf_options_valid(const AtlasRibbonFitOptions *o)
{
    return o != NULL &&
           isfinite(o->slice_spacing) && o->slice_spacing > 0.0 &&
           isfinite(o->observation_u_spacing) &&
           o->observation_u_spacing > 0.0 &&
           isfinite(o->fit_u_spacing) && o->fit_u_spacing > 0.0 &&
           isfinite(o->local_xyz_tolerance) &&
           o->local_xyz_tolerance >= 0.0 &&
           isfinite(o->tangent_dot_min) &&
           o->tangent_dot_min >= -1.0 && o->tangent_dot_min <= 1.0 &&
           isfinite(o->max_fill_u) && o->max_fill_u >= 0.0 &&
           isfinite(o->max_bridge_stretch) &&
           o->max_bridge_stretch >= 0.0 &&
           isfinite(o->bridge_slack) && o->bridge_slack >= 0.0 &&
           isfinite(o->lambda_position) && o->lambda_position > 0.0 &&
           isfinite(o->lambda_smooth) && o->lambda_smooth >= 0.0 &&
           isfinite(o->lambda_tangent) && o->lambda_tangent >= 0.0 &&
           isfinite(o->lambda_register_vertical) &&
           o->lambda_register_vertical >= 0.0 &&
           o->register_sweeps >= 0 && o->register_sweeps <= 100 &&
           isfinite(o->collision_relaxation) &&
           o->collision_relaxation > 0.0 && o->collision_relaxation <= 1.0 &&
           o->collision_rounds >= 0 && o->collision_rounds <= 100 &&
           o->collision_sweeps >= 0 && o->collision_sweeps <= 100 &&
           o->collision_polish_rounds >= 0 &&
           o->collision_polish_rounds <= o->collision_rounds &&
           o->mode >= ATLAS_RIBBON_OBSERVATIONS &&
           o->mode <= ATLAS_RIBBON_REGISTERED_RIBBON;
}

static int arf_compare_observation(const void *pa, const void *pb)
{
    const AtlasRibbonObservation *a =
        (const AtlasRibbonObservation *)pa;
    const AtlasRibbonObservation *b =
        (const AtlasRibbonObservation *)pb;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->column != b->column)
        return a->column < b->column ? -1 : 1;
    if (a->chart != b->chart) return a->chart < b->chart ? -1 : 1;
    if (a->face != b->face) return a->face < b->face ? -1 : 1;
    return 0;
}
static int arf_compare_cluster_rank(const void *pa, const void *pb)
{
    const ArfCluster *a = (const ArfCluster *)pa;
    const ArfCluster *b = (const ArfCluster *)pb;
    if (a->rank != b->rank) return a->rank < b->rank ? -1 : 1;
    if (a->chart != b->chart) return a->chart < b->chart ? -1 : 1;
    return 0;
}
static int arf_compare_register_ref(const void *pa, const void *pb)
{
    const ArfRegisterRef *a = (const ArfRegisterRef *)pa;
    const ArfRegisterRef *b = (const ArfRegisterRef *)pb;
    if (a->u < b->u) return -1;
    if (a->u > b->u) return 1;
    if (a->rank != b->rank) return a->rank < b->rank ? -1 : 1;
    return a->sample < b->sample ? -1 : a->sample > b->sample;
}



static int arf_known_reason(const AtlasSolution *solution,
                            int32_t chart0, int32_t chart1)
{
    if (chart0 > chart1) {
        int32_t t = chart0;
        chart0 = chart1;
        chart1 = t;
    }
    int answer = -1;
    for (size_t i = 0; i < solution->nresiduals; i++) {
        int32_t a = solution->residual[i].chart0;
        int32_t b = solution->residual[i].chart1;
        if (a > b) {
            int32_t t = a;
            a = b;
            b = t;
        }
        if (a != chart0 || b != chart1) continue;
        int reason = solution->residual[i].reason;
        if (reason == 4) return reason;
        if (reason > answer) answer = reason;
    }
    return answer;
}

static int arf_face_chart(const AtlasSolution *solution,
                          const int32_t vertex[3])
{
    int32_t c = solution->vertex_chart[vertex[0]];
    if (c < 0 || solution->vertex_chart[vertex[1]] != c ||
        solution->vertex_chart[vertex[2]] != c ||
        (size_t)c >= solution->ncharts)
        return -1;
    return c;
}

typedef struct {
    double u;
    double p[2];
} ArfEndpoint;

static int arf_face_plane_segment(
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const int32_t vertex[3],
    double plane_v,
    ArfEndpoint endpoint[2])
{
    int count = 0;
    for (int e = 0; e < 3; e++) {
        int32_t a = vertex[e];
        int32_t b = vertex[(e + 1) % 3];
        double va = solution->v[a];
        double vb = solution->v[b];
        double lo = va < vb ? va : vb;
        double hi = va < vb ? vb : va;
        if (!(lo < plane_v && plane_v < hi)) continue;
        if (count >= 2) return 0;
        double t = (plane_v - va) / (vb - va);
        endpoint[count].u =
            (1.0 - t) * solution->u[a] + t * solution->u[b];
        float p[3];
        for (int d = 0; d < 3; d++)
            p[d] = (float)((1.0 - t) *
                           (double)ps->verts[(size_t)a * 3 + (size_t)d] +
                           t * (double)ps->verts[(size_t)b * 3 + (size_t)d]);
        arf_project_point(set, p, endpoint[count].p);
        count++;
    }
    return count == 2 ? 1 : 0;
}

static int arf_emit_observations(
    const PieceSet *ps,
    const AtlasSolution *solution,
    AtlasRibbonObservationSet *set,
    AtlasRibbonObservation *out,
    size_t capacity,
    size_t *out_count,
    int collect_stats)
{
    size_t count = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!solution->face_keep[f]) continue;
        int32_t vertex[3] = {
            ps->faces[f * 3],
            ps->faces[f * 3 + 1],
            ps->faces[f * 3 + 2]
        };
        int chart = arf_face_chart(solution, vertex);
        if (chart < 0) {
            if (collect_stats) set->invalid_chart_faces++;
            continue;
        }
        double vlo = solution->v[vertex[0]];
        double vhi = vlo;
        for (int k = 1; k < 3; k++) {
            double vv = solution->v[vertex[k]];
            if (vv < vlo) vlo = vv;
            if (vv > vhi) vhi = vv;
        }
        int64_t row0 = (int64_t)ceil((vlo - set->v0) / set->dv);
        int64_t row1 = (int64_t)floor((vhi - set->v0) / set->dv);
        if (row0 < 0) row0 = 0;
        if (row1 >= (int64_t)set->nrows) row1 = (int64_t)set->nrows - 1;
        int face_sliced = 0;
        for (int64_t row = row0; row <= row1; row++) {
            double plane_v = set->v0 + (double)row * set->dv;
            ArfEndpoint endpoint[2];
            if (!arf_face_plane_segment(ps, solution, set, vertex,
                                        plane_v, endpoint))
                continue;
            face_sliced = 1;
            if (collect_stats) set->slice_segments++;
            double du = endpoint[1].u - endpoint[0].u;
            double tangent[2] = {
                endpoint[1].p[0] - endpoint[0].p[0],
                endpoint[1].p[1] - endpoint[0].p[1]
            };
            if (du < 0.0) {
                du = -du;
                tangent[0] = -tangent[0];
                tangent[1] = -tangent[1];
                ArfEndpoint swap = endpoint[0];
                endpoint[0] = endpoint[1];
                endpoint[1] = swap;
            }
            if (!(du > 1.0e-10) || arf_normalize2(tangent) < 1.0e-10) {
                if (collect_stats) set->degenerate_segments++;
                continue;
            }
            int64_t col0 = (int64_t)ceil(
                (endpoint[0].u - set->observation_u0) /
                set->observation_du);
            int64_t col1 = (int64_t)floor(
                (endpoint[1].u - set->observation_u0) /
                set->observation_du);
            for (int64_t column = col0; column <= col1; column++) {
                if (column < INT32_MIN || column > INT32_MAX) return -1;
                double u = set->observation_u0 +
                           (double)column * set->observation_du;
                double eps = 1.0e-10 *
                             (1.0 + fabs(endpoint[0].u) +
                              fabs(endpoint[1].u));
                if (!(u > endpoint[0].u + eps &&
                      u < endpoint[1].u - eps))
                    continue;
                if (out != NULL) {
                    if (count >= capacity) return -1;
                    double t = (u - endpoint[0].u) /
                               (endpoint[1].u - endpoint[0].u);
                    AtlasRibbonObservation *dst = &out[count];
                    dst->row = (int32_t)row;
                    dst->column = (int32_t)column;
                    dst->chart = chart;
                    dst->face = (int32_t)f;
                    dst->u = u;
                    dst->v = plane_v;
                    dst->p[0] = (1.0 - t) * endpoint[0].p[0] +
                                t * endpoint[1].p[0];
                    dst->p[1] = (1.0 - t) * endpoint[0].p[1] +
                                t * endpoint[1].p[1];
                    dst->tangent[0] = tangent[0];
                    dst->tangent[1] = tangent[1];
                }
                count++;
            }
        }
        if (collect_stats && face_sliced) set->sliced_faces++;
    }
    *out_count = count;
    return 0;
}

int AtlasRibbonFit_build_observations(
    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const ScaffoldCalib *cal,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonObservationSet *out)
{
    if (arena == NULL || ps == NULL || solution == NULL || cal == NULL ||
        out == NULL || !arf_options_valid(opts) ||
        solution->nvertices != ps->nv || solution->nfaces != ps->nf)
        return -1;
    memset(out, 0, sizeof *out);
    out->input_faces = ps->nf;
    arf_make_frame(cal, out->axis, out->basis0, out->basis1);
    for (int d = 0; d < 3; d++)
        out->axis_point[d] = (double)cal->axis_point[d];

    double umin = DBL_MAX, umax = -DBL_MAX;
    double vmin = DBL_MAX, vmax = -DBL_MAX;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!solution->face_keep[f]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            double u = solution->u[vertex];
            double v = solution->v[vertex];
            if (!isfinite(u) || !isfinite(v)) return -1;
            if (u < umin) umin = u;
            if (u > umax) umax = u;
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }
    }
    if (!(umin < umax) || !(vmin < vmax)) return -1;
    out->u_min = umin;
    out->u_max = umax;
    out->v_min = vmin;
    out->v_max = vmax;
    out->observation_du = opts->observation_u_spacing;
    out->observation_u0 =
        floor(umin / out->observation_du) * out->observation_du +
        0.5 * out->observation_du;
    out->dv = opts->slice_spacing;
    out->v0 = floor(vmin / out->dv) * out->dv + 0.5 * out->dv;
    int64_t nrows = (int64_t)floor((vmax - out->v0) / out->dv) + 1;
    if (nrows <= 0 || (uint64_t)nrows > (uint64_t)INT32_MAX)
        return -1;
    out->nrows = (size_t)nrows;

    size_t count = 0;
    if (arf_emit_observations(ps, solution, out, NULL, 0,
                              &count, 0) != 0 ||
        count == 0)
        return -1;
    out->observation = (AtlasRibbonObservation *)ARENA_ALLOC(
        arena, count * sizeof(AtlasRibbonObservation));
    out->nobservation = count;
    size_t emitted = 0;
    if (arf_emit_observations(ps, solution, out, out->observation,
                              count, &emitted, 1) != 0 ||
        emitted != count)
        return -1;
    qsort(out->observation, out->nobservation,
          sizeof out->observation[0], arf_compare_observation);

    size_t max_group = 0, nkey = 0;
    for (size_t first = 0; first < out->nobservation;) {
        size_t last = first + 1;
        while (last < out->nobservation &&
               out->observation[last].row ==
                   out->observation[first].row &&
               out->observation[last].column ==
                   out->observation[first].column)
            last++;
        size_t size = last - first;
        if (size > max_group) max_group = size;
        nkey++;
        first = last;
    }
    out->target = (AtlasRibbonTarget *)ARENA_ALLOC(
        arena, nkey * sizeof(AtlasRibbonTarget));
    out->layer_sample = (AtlasRibbonLayerSample *)ARENA_ALLOC(
        arena, out->nobservation * sizeof(AtlasRibbonLayerSample));
    ArfChartSample *chart_sample = (ArfChartSample *)ARENA_ALLOC(
        arena, max_group * sizeof(ArfChartSample));
    ArfCluster *cluster = (ArfCluster *)ARENA_ALLOC(
        arena, max_group * sizeof(ArfCluster));

    size_t target_count = 0, layer_count = 0;
    for (size_t first = 0; first < out->nobservation;) {
        size_t last = first + 1;
        while (last < out->nobservation &&
               out->observation[last].row ==
                   out->observation[first].row &&
               out->observation[last].column ==
                   out->observation[first].column)
            last++;

        size_t nsample = 0;
        for (size_t begin = first; begin < last;) {
            size_t end = begin + 1;
            while (end < last &&
                   out->observation[end].chart ==
                       out->observation[begin].chart)
                end++;
            ArfChartSample *sample = &chart_sample[nsample++];
            memset(sample, 0, sizeof *sample);
            sample->chart = out->observation[begin].chart;
            sample->face = out->observation[begin].face;
            sample->count = (uint32_t)(end - begin);
            for (size_t i = begin; i < end; i++) {
                sample->p[0] += out->observation[i].p[0];
                sample->p[1] += out->observation[i].p[1];
                sample->tangent[0] += out->observation[i].tangent[0];
                sample->tangent[1] += out->observation[i].tangent[1];
            }
            double inv = 1.0 / (double)(end - begin);
            sample->p[0] *= inv;
            sample->p[1] *= inv;
            arf_normalize2(sample->tangent);
            begin = end;
        }

        size_t ncluster = 0;
        for (size_t i = 0; i < nsample; i++) {
            size_t best = SIZE_MAX;
            double best_distance = DBL_MAX;
            for (size_t c = 0; c < ncluster; c++) {
                double distance =
                    arf_distance2(chart_sample[i].p, cluster[c].p);
                double tangent_dot =
                    chart_sample[i].tangent[0] * cluster[c].tangent[0] +
                    chart_sample[i].tangent[1] * cluster[c].tangent[1];
                if (distance <= opts->local_xyz_tolerance &&
                    tangent_dot >= opts->tangent_dot_min &&
                    distance < best_distance) {
                    best = c;
                    best_distance = distance;
                }
            }
            if (best == SIZE_MAX) {
                ArfCluster *dst = &cluster[ncluster++];
                memset(dst, 0, sizeof *dst);
                dst->p[0] = chart_sample[i].p[0];
                dst->p[1] = chart_sample[i].p[1];
                dst->tangent[0] = chart_sample[i].tangent[0];
                dst->tangent[1] = chart_sample[i].tangent[1];
                dst->charts = 1;
                dst->chart = chart_sample[i].chart;
                dst->face = chart_sample[i].face;
                dst->rank = solution->chart[dst->chart].rank;
            } else {
                ArfCluster *dst = &cluster[best];
                double old = (double)dst->charts;
                double den = old + 1.0;
                dst->p[0] = (old * dst->p[0] +
                             chart_sample[i].p[0]) / den;
                dst->p[1] = (old * dst->p[1] +
                             chart_sample[i].p[1]) / den;
                dst->tangent[0] += chart_sample[i].tangent[0];
                dst->tangent[1] += chart_sample[i].tangent[1];
                arf_normalize2(dst->tangent);
                dst->charts++;
                uint64_t rank = solution->chart[chart_sample[i].chart].rank;
                if (rank < dst->rank) {
                    dst->rank = rank;
                    dst->chart = chart_sample[i].chart;
                    dst->face = chart_sample[i].face;
                }
            }
        }
        qsort(cluster, ncluster, sizeof cluster[0],
              arf_compare_cluster_rank);
        for (size_t c = 0; c < ncluster; c++) {
            if (layer_count >= out->nobservation) return -1;
            const ArfCluster *src = &cluster[c];
            int32_t chart = src->chart;
            if (chart < 0 || (size_t)chart >= solution->ncharts ||
                solution->chart[chart].chart != chart)
                return -1;
            const AtlasSolutionChart *chart_info = &solution->chart[chart];
            AtlasRibbonLayerSample *dst = &out->layer_sample[layer_count++];
            memset(dst, 0, sizeof *dst);
            dst->row = out->observation[first].row;
            dst->source_column = out->observation[first].column;
            dst->chart = chart;
            dst->face = src->face;
            dst->group = chart_info->group;
            /* dst->winding is assigned by the calibrated per-sample snap
             * at the end of this function, never copied per chart. */
            dst->rank = src->rank;
            dst->charts = src->charts;
            dst->cluster_index = (uint32_t)c;
            dst->cluster_count = (uint32_t)ncluster;
            dst->source_u = out->observation[first].u;
            dst->v = out->observation[first].v;
            dst->p[0] = src->p[0];
            dst->p[1] = src->p[1];
            dst->tangent[0] = src->tangent[0];
            dst->tangent[1] = src->tangent[1];
            if (ncluster > 1) out->remote_layer_samples++;
        }


        AtlasRibbonTarget *target = &out->target[target_count++];
        memset(target, 0, sizeof *target);
        target->row = out->observation[first].row;
        target->column = out->observation[first].column;
        target->u = out->observation[first].u;
        target->v = out->observation[first].v;
        target->observations = (uint32_t)(last - first);
        target->charts = (uint32_t)nsample;
        target->clusters = (uint32_t)ncluster;
        target->chart0 = ncluster > 0 ? cluster[0].chart : -1;
        target->face0 = ncluster > 0 ? cluster[0].face : -1;
        target->chart1 = ncluster > 1 ? cluster[1].chart : -1;
        target->face1 = ncluster > 1 ? cluster[1].face : -1;
        target->known_reason = ncluster > 1
            ? arf_known_reason(solution, target->chart0, target->chart1)
            : -1;
        for (size_t a = 0; a < ncluster; a++)
            for (size_t b = a + 1; b < ncluster; b++) {
                double distance = arf_distance2(cluster[a].p, cluster[b].p);
                if (distance > target->max_cluster_separation)
                    target->max_cluster_separation = distance;
            }
        if (ncluster == 1) {
            target->accepted = 1;
            target->p[0] = cluster[0].p[0];
            target->p[1] = cluster[0].p[1];
            target->tangent[0] = cluster[0].tangent[0];
            target->tangent[1] = cluster[0].tangent[1];
            out->accepted_targets++;
            if (nsample > 1) out->duplicate_targets++;
        } else {
            target->accepted = 0;
            out->conflict_targets++;
            out->conflict_clusters += ncluster;
            if (target->known_reason == 4)
                out->known_topology_conflicts++;
        }
        first = last;
    }
    out->nlayer_samples = layer_count;
    out->ntarget = target_count;

    /* Per-sample integer winding snap.  The checkpoint's continuous wind
     * is in the producer's sense*phi convention; the span math downstream
     * counts wraps of THIS side's angle alpha = mod_tau(direction*theta).
     * The two conventions relate by an unknown fixed rotation and possible
     * reflection, so the sign is calibrated from the data itself: the
     * circular resultant of (sigma*wind - alpha/2pi) is near 1 for the
     * sign whose convention matches and near 0 for the other.  Snapping
     * each sample against its own alpha makes the integer increment at
     * the consumer's own theta cut, tolerating half a turn of continuous
     * error and charts that straddle the cut. */
    if (layer_count > 0) {
        int direction = arf_registration_winding_direction(out);
        double cos_sum[2] = {0.0, 0.0};
        double sin_sum[2] = {0.0, 0.0};
        for (size_t i = 0; i < layer_count; i++) {
            const AtlasRibbonLayerSample *s = &out->layer_sample[i];
            double wind = solution->chart[s->chart].wind;
            double alpha = arf_mod_tau(
                (double)direction * atan2(s->p[1], s->p[0])) /
                (2.0 * M_PI);
            for (int sign = 0; sign < 2; sign++) {
                double sigma = sign == 0 ? 1.0 : -1.0;
                double turns = sigma * wind - alpha;
                cos_sum[sign] += cos(2.0 * M_PI * turns);
                sin_sum[sign] += sin(2.0 * M_PI * turns);
            }
        }
        double r_plus = hypot(cos_sum[0], sin_sum[0]) /
                        (double)layer_count;
        double r_minus = hypot(cos_sum[1], sin_sum[1]) /
                         (double)layer_count;
        double sigma = r_plus >= r_minus ? 1.0 : -1.0;
        double resultant = r_plus >= r_minus ? r_plus : r_minus;
        if (resultant >= 0.25) {
            for (size_t i = 0; i < layer_count; i++) {
                AtlasRibbonLayerSample *s = &out->layer_sample[i];
                double wind = solution->chart[s->chart].wind;
                double alpha = arf_mod_tau(
                    (double)direction * atan2(s->p[1], s->p[0])) /
                    (2.0 * M_PI);
                s->winding = (int32_t)lround(sigma * wind - alpha);
            }
            fprintf(stderr, "[atlas_ribbon_fit] wind snap: sign=%+d "
                    "resultant=%.3f (other %.3f) direction=%+d\n",
                    (int)sigma, resultant,
                    sigma > 0.0 ? r_minus : r_plus, direction);
        } else {
            for (size_t i = 0; i < layer_count; i++) {
                AtlasRibbonLayerSample *s = &out->layer_sample[i];
                s->winding = (int32_t)lround(
                    solution->chart[s->chart].wind);
            }
            fprintf(stderr, "[atlas_ribbon_fit] wind snap: frames "
                    "incoherent (resultant %.3f/%.3f) -- per-chart "
                    "rounding only\n", r_plus, r_minus);
        }
    }
    return 0;
}

static double arf_mod_tau(double angle)
{
    double value = fmod(angle, 2.0 * M_PI);
    if (value < 0.0) value += 2.0 * M_PI;
    return value;
}

static int arf_registration_winding_direction(
    const AtlasRibbonObservationSet *set)
{
    size_t positive = 0, negative = 0;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        double radius = sqrt(s->p[0] * s->p[0] + s->p[1] * s->p[1]);
        if (!(radius > 1.0e-6)) continue;
        double cross = s->p[0] * s->tangent[1] -
                       s->p[1] * s->tangent[0];
        if (cross > 0.05 * radius) positive++;
        else if (cross < -0.05 * radius) negative++;
    }
    return positive >= negative ? 1 : -1;
}

static int arf_layer_samples_local(const AtlasRibbonLayerSample *a,
                                   const AtlasRibbonLayerSample *b,
                                   const AtlasRibbonFitOptions *opts)
{
    double tangent_dot = a->tangent[0] * b->tangent[0] +
                         a->tangent[1] * b->tangent[1];
    return arf_distance2(a->p, b->p) <= opts->local_xyz_tolerance &&
           tangent_dot >= opts->tangent_dot_min;
}

static double arf_spiral_arc(double radius0, double radius1,
                             double angle_delta)
{
    enum { INTERVALS = 32 };
    double dr = radius1 - radius0;
    double sum = 0.0;
    for (int i = 0; i <= INTERVALS; i++) {
        double t = (double)i / (double)INTERVALS;
        double radius = radius0 + t * dr;
        double speed = sqrt(dr * dr +
                            radius * radius * angle_delta * angle_delta);
        double weight = i == 0 || i == INTERVALS ? 1.0
                      : (i & 1) ? 4.0 : 2.0;
        sum += weight * speed;
    }
    return sum / (3.0 * (double)INTERVALS);
}

/* Estimate the amount of intrinsic U needed to continue from the lower-rank
 * observation to the higher-rank observation without cutting a Cartesian
 * chord through the scroll.  The chart winding disambiguates a positive turn
 * when it is coherent; rank plus measured tangents supplies the direction. */
static double arf_required_registration_arc(
    const AtlasRibbonLayerSample *lo,
    const AtlasRibbonLayerSample *hi,
    int direction)
{
    double radius0 = sqrt(lo->p[0] * lo->p[0] + lo->p[1] * lo->p[1]);
    double radius1 = sqrt(hi->p[0] * hi->p[0] + hi->p[1] * hi->p[1]);
    double theta0 = atan2(lo->p[1], lo->p[0]);
    double theta1 = atan2(hi->p[1], hi->p[0]);
    double alpha0 = arf_mod_tau((double)direction * theta0);
    double alpha1 = arf_mod_tau((double)direction * theta1);
    int64_t winding_delta = (int64_t)hi->winding - (int64_t)lo->winding;
    double angle_delta;
    if (winding_delta > 0 && winding_delta < INT32_MAX) {
        angle_delta = alpha1 - alpha0 +
                      2.0 * M_PI * (double)winding_delta;
        while (angle_delta <= 1.0e-12) angle_delta += 2.0 * M_PI;
    } else {
        /* A tiny same-winding phase reversal is registration noise, not
         * evidence that the missing ribbon traverses an extra full turn. */
        angle_delta = fabs(atan2(sin(alpha1 - alpha0),
                                 cos(alpha1 - alpha0)));
    }
    double arc = arf_spiral_arc(radius0, radius1, angle_delta);
    double chord = arf_distance2(lo->p, hi->p);
    if (!isfinite(arc) || arc < chord) arc = chord;
    return arc;
}

static int arf_register_add_bound(
    ArfRegisterBounds *bounds,
    const AtlasRibbonObservationSet *set,
    size_t sample0,
    size_t sample1,
    AtlasRibbonRegisterKind kind,
    int winding_direction,
    int keep_nonpositive,
    double *max_required_arc)
{
    if (sample0 >= set->nlayer_samples || sample1 >= set->nlayer_samples)
        return -1;
    const AtlasRibbonLayerSample *lo = &set->layer_sample[sample0];
    const AtlasRibbonLayerSample *hi = &set->layer_sample[sample1];
    if (lo->rank == hi->rank) return 0;
    if (lo->rank > hi->rank) {
        const AtlasRibbonLayerSample *swap = lo;
        lo = hi;
        hi = swap;
        size_t index_swap = sample0;
        sample0 = sample1;
        sample1 = index_swap;
    }
    if (hi->rank >= bounds->ncharts) return -1;
    double required = arf_required_registration_arc(
        lo, hi, winding_direction);
    double source_delta = hi->source_u - lo->source_u;
    double lower = required - source_delta;
    if (!isfinite(required) || !isfinite(lower)) return -1;
    if (required > *max_required_arc) *max_required_arc = required;
    if (!keep_nonpositive && !(lower > 1.0e-9)) return 0;
    size_t at = (size_t)lo->rank * bounds->ncharts + (size_t)hi->rank;
    bounds->updates[at]++;
    if (lower <= bounds->lower[at] + 1.0e-9) return 0;
    bounds->lower[at] = lower;
    bounds->kind[at] = (uint8_t)kind;
    bounds->sample_lo[at] = (uint64_t)sample0;
    bounds->sample_hi[at] = (uint64_t)sample1;
    bounds->required_arc[at] = required;
    bounds->endpoint_distance[at] = arf_distance2(lo->p, hi->p);
    return 1;
}

static void arf_register_solve_bounds(const ArfRegisterBounds *bounds,
                                      double *rank_shift)
{
    for (size_t hi = 0; hi < bounds->ncharts; hi++) {
        double value = hi > 0 ? rank_shift[hi - 1] : 0.0;
        for (size_t lo = 0; lo < hi; lo++) {
            double lower = bounds->lower[lo * bounds->ncharts + hi];
            if (!(lower > -DBL_MAX)) continue;
            double candidate = rank_shift[lo] + lower;
            if (candidate > value) value = candidate;
        }
        rank_shift[hi] = value;
    }
}

static size_t arf_register_scan_rows(
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const double *rank_shift,
    ArfRegisterBounds *bounds,
    ArfRegisterRef *ref,
    int winding_direction,
    int add_bounds,
    size_t *out_added,
    size_t *out_local_inversions,
    size_t *out_same_rank_unresolved,
    double *max_required_arc)
{
    size_t inversions = 0, added = 0, local_inversions = 0;
    size_t same_rank_unresolved = 0;
    for (size_t first = 0; first < set->nlayer_samples;) {
        size_t last = first + 1;
        while (last < set->nlayer_samples &&
               set->layer_sample[last].row == set->layer_sample[first].row)
            last++;
        size_t count = last - first;
        for (size_t j = 0; j < count; j++) {
            size_t sample = first + j;
            ref[j].sample = sample;
            ref[j].rank = set->layer_sample[sample].rank;
            ref[j].u = set->layer_sample[sample].source_u +
                         rank_shift[ref[j].rank];
        }
        qsort(ref, count, sizeof ref[0], arf_compare_register_ref);
        for (size_t j = 1; j < count; j++) {
            const AtlasRibbonLayerSample *a =
                &set->layer_sample[ref[j - 1].sample];
            const AtlasRibbonLayerSample *b =
                &set->layer_sample[ref[j].sample];
            if (a->rank == b->rank) {
                double gap = ref[j].u - ref[j - 1].u;
                double distance = arf_distance2(a->p, b->p);
                if (distance > opts->max_bridge_stretch * gap +
                               opts->bridge_slack)
                    same_rank_unresolved++;
                continue;
            }
            int local = arf_layer_samples_local(a, b, opts);
            if (a->rank > b->rank) {
                if (local) {
                    local_inversions++;
                    continue;
                }
                inversions++;
                if (add_bounds) {
                    int rc = arf_register_add_bound(
                        bounds, set, ref[j - 1].sample, ref[j].sample,
                        ATLAS_RIBBON_REGISTER_ORDER, winding_direction, 0,
                        max_required_arc);
                    if (rc < 0) return SIZE_MAX;
                    added += (size_t)rc;
                }
                continue;
            }
            if (local) continue;
            double gap = ref[j].u - ref[j - 1].u;
            double distance = arf_distance2(a->p, b->p);
            if (distance <= opts->max_bridge_stretch * gap +
                            opts->bridge_slack)
                continue;
            if (add_bounds) {
                int rc = arf_register_add_bound(
                    bounds, set, ref[j - 1].sample, ref[j].sample,
                    ATLAS_RIBBON_REGISTER_METRIC, winding_direction, 0,
                    max_required_arc);
                if (rc < 0) return SIZE_MAX;
                added += (size_t)rc;
            }
        }
        first = last;
    }
    if (out_added != NULL) *out_added = added;
    if (out_local_inversions != NULL) *out_local_inversions = local_inversions;
    if (out_same_rank_unresolved != NULL)
        *out_same_rank_unresolved = same_rank_unresolved;
    return inversions;
}

static int arf_register_u(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonFitResult *out)
{
    if (set->nlayer_samples == 0 || solution->ncharts == 0 ||
        solution->ncharts > SIZE_MAX / solution->ncharts)
        return -1;
    size_t ncharts = solution->ncharts;
    size_t matrix_size = ncharts * ncharts;
    if (matrix_size > SIZE_MAX / sizeof(double)) return -1;

    size_t *rank_to_chart = (size_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(*rank_to_chart));
    uint8_t *rank_seen = (uint8_t *)ARENA_CALLOC(arena, ncharts, 1);
    for (size_t c = 0; c < ncharts; c++) {
        if (solution->chart[c].chart != (int32_t)c ||
            solution->chart[c].rank >= ncharts ||
            rank_seen[solution->chart[c].rank])
            return -1;
        rank_seen[solution->chart[c].rank] = 1;
        rank_to_chart[solution->chart[c].rank] = c;
    }
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        if (s->chart < 0 || (size_t)s->chart >= ncharts ||
            s->rank != solution->chart[s->chart].rank)
            return -1;
    }

    ArfRegisterBounds bounds;
    memset(&bounds, 0, sizeof bounds);
    bounds.ncharts = ncharts;
    bounds.lower = (double *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.lower));
    bounds.kind = (uint8_t *)ARENA_CALLOC(arena, matrix_size, 1);
    bounds.updates = (uint32_t *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*bounds.updates));
    bounds.sample_lo = (uint64_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.sample_lo));
    bounds.sample_hi = (uint64_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.sample_hi));
    bounds.required_arc = (double *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*bounds.required_arc));
    bounds.endpoint_distance = (double *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*bounds.endpoint_distance));
    for (size_t i = 0; i < matrix_size; i++) {
        bounds.lower[i] = -DBL_MAX;
        bounds.sample_lo[i] = bounds.sample_hi[i] = UINT64_MAX;
    }

    out->register_winding_direction =
        arf_registration_winding_direction(set);
    double max_required_arc = 0.0;

    /* Every remote cluster at an otherwise identical source UV coordinate is
     * a retained physical layer.  Consecutive immutable ranks must receive
     * enough parameter length to travel between them. */
    for (size_t first = 0; first < set->nlayer_samples;) {
        const AtlasRibbonLayerSample *head = &set->layer_sample[first];
        size_t count = head->cluster_count;
        if (count == 0 || first + count > set->nlayer_samples) return -1;
        for (size_t j = 0; j < count; j++) {
            const AtlasRibbonLayerSample *s = &set->layer_sample[first + j];
            if (s->row != head->row ||
                s->source_column != head->source_column ||
                s->cluster_count != count || s->cluster_index != j)
                return -1;
        }
        for (size_t j = 1; j < count; j++) {
            int rc = arf_register_add_bound(
                &bounds, set, first + j - 1, first + j,
                ATLAS_RIBBON_REGISTER_REMOTE,
                out->register_winding_direction, 0, &max_required_arc);
            if (rc < 0) return -1;
        }
        first += count;
    }

    /* Seed the exact failures which exposed the old Cartesian formulation:
     * consecutive single-layer observations whose chord cannot fit in their
     * source-U separation. */
    size_t previous = SIZE_MAX;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *current = &set->layer_sample[i];
        if (current->cluster_count != 1) continue;
        if (previous != SIZE_MAX &&
            set->layer_sample[previous].row == current->row) {
            const AtlasRibbonLayerSample *prior = &set->layer_sample[previous];
            double gap = current->source_u - prior->source_u;
            double distance = arf_distance2(prior->p, current->p);
            int within_fill = opts->max_fill_u <= 0.0 ||
                              gap <= opts->max_fill_u;
            if (within_fill && gap >= 0.0 &&
                distance > opts->max_bridge_stretch * gap +
                           opts->bridge_slack &&
                prior->rank != current->rank) {
                int rc = arf_register_add_bound(
                    &bounds, set, previous, i,
                    ATLAS_RIBBON_REGISTER_METRIC,
                    out->register_winding_direction, 0, &max_required_arc);
                if (rc < 0) return -1;
            }
        }
        previous = i;
    }

    size_t max_row_samples = 0;
    for (size_t first = 0; first < set->nlayer_samples;) {
        size_t last = first + 1;
        while (last < set->nlayer_samples &&
               set->layer_sample[last].row == set->layer_sample[first].row)
            last++;
        if (last - first > max_row_samples) max_row_samples = last - first;
        first = last;
    }
    ArfRegisterRef *ref = (ArfRegisterRef *)ARENA_ALLOC(
        arena, max_row_samples * sizeof(*ref));
    double *rank_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(*rank_shift));

    size_t initial_local = 0, initial_same_rank = 0;
    out->register_inversions_before = arf_register_scan_rows(
        set, opts, rank_shift, &bounds, ref,
        out->register_winding_direction, 0, NULL,
        &initial_local, &initial_same_rank, &max_required_arc);
    if (out->register_inversions_before == SIZE_MAX) return -1;

    enum { MAX_REGISTER_ITERATIONS = 64 };
    int converged = 0;
    for (int iteration = 0; iteration < MAX_REGISTER_ITERATIONS; iteration++) {
        arf_register_solve_bounds(&bounds, rank_shift);
        size_t added = 0;
        size_t inversions = arf_register_scan_rows(
            set, opts, rank_shift, &bounds, ref,
            out->register_winding_direction, 1, &added,
            NULL, NULL, &max_required_arc);
        if (inversions == SIZE_MAX) return -1;
        out->register_iterations = (size_t)iteration + 1;
        if (added == 0) {
            converged = 1;
            break;
        }
    }
    arf_register_solve_bounds(&bounds, rank_shift);

    size_t final_local = 0, same_rank_unresolved = 0;
    out->register_inversions_after = arf_register_scan_rows(
        set, opts, rank_shift, &bounds, ref,
        out->register_winding_direction, 0, NULL,
        &final_local, &same_rank_unresolved, &max_required_arc);
    if (out->register_inversions_after == SIZE_MAX) return -1;
    out->register_local_inversions = final_local;
    out->register_unresolved = out->register_inversions_after +
                               same_rank_unresolved + (converged ? 0 : 1);

    double weighted_center = 0.0;
    for (size_t i = 0; i < set->nlayer_samples; i++)
        weighted_center += rank_shift[set->layer_sample[i].rank];
    weighted_center /= (double)set->nlayer_samples;
    for (size_t r = 0; r < ncharts; r++) rank_shift[r] -= weighted_center;

    out->registered_u = (double *)ARENA_ALLOC(
        arena, set->nlayer_samples * sizeof(*out->registered_u));
    out->chart_shift = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(*out->chart_shift));
    double source_min = DBL_MAX, source_max = -DBL_MAX;
    double registered_min = DBL_MAX, registered_max = -DBL_MAX;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        double u = s->source_u + rank_shift[s->rank];
        out->registered_u[i] = u;
        if (s->source_u < source_min) source_min = s->source_u;
        if (s->source_u > source_max) source_max = s->source_u;
        if (u < registered_min) registered_min = u;
        if (u > registered_max) registered_max = u;
    }
    double shift_sum2 = 0.0;
    for (size_t r = 0; r < ncharts; r++) {
        size_t chart = rank_to_chart[r];
        out->chart_shift[chart] = rank_shift[r];
        shift_sum2 += rank_shift[r] * rank_shift[r];
        double magnitude = fabs(rank_shift[r]);
        if (magnitude > out->register_shift_max)
            out->register_shift_max = magnitude;
    }
    out->register_shift_rms = sqrt(shift_sum2 / (double)ncharts);
    out->register_width_before = source_max - source_min;
    out->register_width_after = registered_max - registered_min;
    out->register_max_required_arc = max_required_arc;

    size_t nbound = 0;
    for (size_t lo = 0; lo < ncharts; lo++)
        for (size_t hi = lo + 1; hi < ncharts; hi++)
            if (bounds.lower[lo * ncharts + hi] > -DBL_MAX) nbound++;
    out->register_bound = nbound > 0
        ? (AtlasRibbonRegisterBound *)ARENA_ALLOC(
            arena, nbound * sizeof(*out->register_bound))
        : NULL;
    out->nregister_bounds = nbound;
    size_t at = 0;
    for (size_t lo = 0; lo < ncharts; lo++)
    for (size_t hi = lo + 1; hi < ncharts; hi++) {
        size_t bi = lo * ncharts + hi;
        if (!(bounds.lower[bi] > -DBL_MAX)) continue;
        AtlasRibbonRegisterBound *dst = &out->register_bound[at++];
        memset(dst, 0, sizeof *dst);
        dst->chart_lo = (int32_t)rank_to_chart[lo];
        dst->chart_hi = (int32_t)rank_to_chart[hi];
        dst->rank_lo = lo;
        dst->rank_hi = hi;
        dst->sample_lo = bounds.sample_lo[bi];
        dst->sample_hi = bounds.sample_hi[bi];
        dst->updates = bounds.updates[bi];
        dst->kind = bounds.kind[bi];
        dst->lower_shift_delta = bounds.lower[bi];
        dst->final_shift_delta = rank_shift[hi] - rank_shift[lo];
        dst->required_arc = bounds.required_arc[bi];
        dst->endpoint_distance = bounds.endpoint_distance[bi];
        if (dst->sample_lo < set->nlayer_samples &&
            dst->sample_hi < set->nlayer_samples) {
            const AtlasRibbonLayerSample *a =
                &set->layer_sample[dst->sample_lo];
            const AtlasRibbonLayerSample *b =
                &set->layer_sample[dst->sample_hi];
            dst->row = a->row == b->row ? a->row : -1;
            dst->source_u_delta = b->source_u - a->source_u;
        } else {
            dst->row = -1;
            dst->source_u_delta = NAN;
        }
        out->register_constraints++;
        if (dst->kind == ATLAS_RIBBON_REGISTER_METRIC)
            out->register_metric_constraints++;
        else if (dst->kind == ATLAS_RIBBON_REGISTER_ORDER)
            out->register_order_constraints++;
    }
    if (at != nbound) return -1;
    return 0;
}
static size_t arf_register_scan_one_row(
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    size_t first,
    size_t last,
    const double *rank_shift,
    ArfRegisterBounds *bounds,
    ArfRegisterRef *ref,
    int winding_direction,
    int add_bounds,
    size_t *out_added,
    size_t *out_local_inversions,
    size_t *out_same_rank_unresolved,
    double *max_required_arc)
{
    size_t count = last - first;
    for (size_t j = 0; j < count; j++) {
        size_t sample = first + j;
        ref[j].sample = sample;
        ref[j].rank = set->layer_sample[sample].rank;
        ref[j].u = set->layer_sample[sample].source_u +
                     rank_shift[ref[j].rank];
    }
    qsort(ref, count, sizeof ref[0], arf_compare_register_ref);
    size_t inversions = 0, added = 0, local_inversions = 0;
    size_t same_rank_unresolved = 0;
    for (size_t j = 1; j < count; j++) {
        const AtlasRibbonLayerSample *a =
            &set->layer_sample[ref[j - 1].sample];
        const AtlasRibbonLayerSample *b =
            &set->layer_sample[ref[j].sample];
        double gap = ref[j].u - ref[j - 1].u;
        double distance = arf_distance2(a->p, b->p);
        if (a->rank == b->rank) {
            if (distance > opts->max_bridge_stretch * gap +
                           opts->bridge_slack)
                same_rank_unresolved++;
            continue;
        }
        int local = arf_layer_samples_local(a, b, opts);
        if (a->rank > b->rank) {
            if (local) {
                local_inversions++;
                continue;
            }
            inversions++;
            if (add_bounds) {
                int rc = arf_register_add_bound(
                    bounds, set, ref[j - 1].sample, ref[j].sample,
                    ATLAS_RIBBON_REGISTER_ORDER, winding_direction, 1,
                    max_required_arc);
                if (rc < 0) return SIZE_MAX;
                added += (size_t)rc;
            }
        } else if (!local &&
                   distance > opts->max_bridge_stretch * gap +
                              opts->bridge_slack && add_bounds) {
            int rc = arf_register_add_bound(
                bounds, set, ref[j - 1].sample, ref[j].sample,
                ATLAS_RIBBON_REGISTER_METRIC, winding_direction, 1,
                max_required_arc);
            if (rc < 0) return SIZE_MAX;
            added += (size_t)rc;
        }
    }
    if (out_added != NULL) *out_added = added;
    if (out_local_inversions != NULL) *out_local_inversions = local_inversions;
    if (out_same_rank_unresolved != NULL)
        *out_same_rank_unresolved = same_rank_unresolved;
    return inversions;
}

static int arf_register_row_qp(
    Arena_T arena,
    const ArfRegisterBounds *bounds,
    const uint32_t *present,
    const uint32_t *local_pair,
    const double *vertical_target_sum,
    const double *vertical_weight,
    double *rank_shift)
{
    size_t ncharts = bounds->ncharts;
    size_t nlocal = 0, nbound = 0;
    for (size_t lo = 0; lo < ncharts; lo++)
    for (size_t hi = lo + 1; hi < ncharts; hi++) {
        size_t at = lo * ncharts + hi;
        if (local_pair[at] > 0) nlocal++;
        if (bounds->lower[at] > -DBL_MAX) nbound++;
    }
    if (nbound == 0 && vertical_target_sum == NULL) {
        memset(rank_shift, 0, ncharts * sizeof(*rank_shift));
        return 0;
    }
    if (ncharts + nlocal > SIZE_MAX / sizeof(MonotoneQpRow) ||
        ncharts + 2 * nlocal > SIZE_MAX / sizeof(MonotoneQpCoeff))
        return -1;
    Arena_Mark mark = Arena_save(arena);
    size_t nrows = ncharts + nlocal;
    size_t ncoeff = ncharts + 2 * nlocal;
    MonotoneQpRow *row = (MonotoneQpRow *)ARENA_ALLOC(
        arena, nrows * sizeof(*row));
    MonotoneQpCoeff *coeff = (MonotoneQpCoeff *)ARENA_ALLOC(
        arena, ncoeff * sizeof(*coeff));
    MonotoneQpBound *bound = (MonotoneQpBound *)ARENA_ALLOC(
        arena, (nbound ? nbound : 1) * sizeof(*bound));
    size_t cursor = 0;
    for (size_t r = 0; r < ncharts; r++) {
        double base_weight = present[r] > 0
                           ? 1.0 + 0.02 * (double)present[r] : 1.0e-4;
        double v_weight = vertical_weight != NULL ? vertical_weight[r] : 0.0;
        row[r].first = cursor;
        row[r].count = 1;
        row[r].target = v_weight > 0.0
                      ? vertical_target_sum[r] / (base_weight + v_weight)
                      : 0.0;
        row[r].weight = base_weight + v_weight;
        row[r].kind = 0;
        row[r].owner = (int32_t)r;
        coeff[cursor].var = (int32_t)r;
        coeff[cursor].value = 1.0;
        cursor++;
    }
    size_t ri = ncharts, bi = 0;
    for (size_t lo = 0; lo < ncharts; lo++)
    for (size_t hi = lo + 1; hi < ncharts; hi++) {
        size_t at = lo * ncharts + hi;
        if (local_pair[at] > 0) {
            row[ri].first = cursor;
            row[ri].count = 2;
            row[ri].target = 0.0;
            row[ri].weight = 25.0 * (double)local_pair[at];
            row[ri].kind = 1;
            row[ri].owner = ri <= (size_t)INT32_MAX ? (int32_t)ri : INT32_MAX;
            coeff[cursor].var = (int32_t)lo;
            coeff[cursor].value = -1.0;
            cursor++;
            coeff[cursor].var = (int32_t)hi;
            coeff[cursor].value = 1.0;
            cursor++;
            ri++;
        }
        if (bounds->lower[at] > -DBL_MAX) {
            bound[bi].lo = (int32_t)lo;
            bound[bi].hi = (int32_t)hi;
            bound[bi].lower = bounds->lower[at];
            bound[bi].stroke = -1;
            bound[bi].edge = bi <= (size_t)INT32_MAX
                           ? (int32_t)bi : INT32_MAX;
            bi++;
        }
    }
    if (ri != nrows || cursor != ncoeff || bi != nbound) {
        Arena_restore(arena, mark);
        return -1;
    }

    arf_register_solve_bounds(bounds, rank_shift);
    double center_sum = 0.0, center_weight = 0.0;
    for (size_t r = 0; r < ncharts; r++)
        if (present[r] > 0) {
            center_sum += (double)present[r] * rank_shift[r];
            center_weight += (double)present[r];
        }
    double center = center_weight > 0.0 ? center_sum / center_weight : 0.0;
    for (size_t r = 0; r < ncharts; r++) rank_shift[r] -= center;
    double *feasible = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(*feasible));
    memcpy(feasible, rank_shift, ncharts * sizeof(*feasible));

    MonotoneQpProblem problem;
    memset(&problem, 0, sizeof problem);
    problem.nvar = ncharts;
    problem.rows = row;
    problem.nrows = nrows;
    problem.coeff = coeff;
    problem.ncoeff = ncoeff;
    problem.bounds = bound;
    problem.nbounds = nbound;
    MonotoneQpOptions options;
    MonotoneQpOptions_default(&options);
    size_t limit = 8 * (ncharts + nbound) + 128;
    options.max_active_iterations = limit > (size_t)INT_MAX
                                  ? INT_MAX : (int)limit;
    MonotoneQpStats stats;
    int rc = MonotoneQp_solve(
        arena, &problem, &options, rank_shift, NULL, &stats);
    if (rc != 0 || stats.max_violation > 1.0e-6) {
        memcpy(rank_shift, feasible, ncharts * sizeof(*rank_shift));
        rc = 1;
    } else {
        rc = 0;
    }
    Arena_restore(arena, mark);
    return rc;
}

static void arf_register_vertical_stats(
    size_t nrows,
    size_t ncharts,
    const uint32_t *row_present,
    const double *field,
    double *out_rms,
    double *out_max)
{
    double sum2 = 0.0, maximum = 0.0;
    size_t count = 0;
    for (size_t row = 1; row < nrows; row++)
    for (size_t chart = 0; chart < ncharts; chart++) {
        if (row_present[(row - 1) * ncharts + chart] == 0 ||
            row_present[row * ncharts + chart] == 0)
            continue;
        double delta = field[row * ncharts + chart] -
                       field[(row - 1) * ncharts + chart];
        sum2 += delta * delta;
        count++;
        if (fabs(delta) > maximum) maximum = fabs(delta);
    }
    *out_rms = count > 0 ? sqrt(sum2 / (double)count) : 0.0;
    *out_max = maximum;
}

static int arf_register_smooth_u_field(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const size_t *rank_to_chart,
    const uint32_t *row_present,
    ArfRegisterBounds *bounds,
    const AtlasRibbonCollisionBound *collision_bound,
    size_t ncollision_bounds,
    const double *phase_delta,
    const double *phase_weight,
    const double *anchor_target,
    double anchor_weight,
    int sweep_limit,
    uint32_t *present,
    uint32_t *local_pair,
    double *rank_shift,
    AtlasRibbonFitResult *out)
{
    size_t nrows = set->nrows, ncharts = solution->ncharts;
    if (nrows > SIZE_MAX / ncharts || sweep_limit < 0 ||
        (ncollision_bounds > 0 && collision_bound == NULL) ||
        ((phase_delta == NULL) != (phase_weight == NULL)) ||
        !isfinite(anchor_weight) || anchor_weight < 0.0 ||
        (anchor_weight > 0.0 && anchor_target == NULL))
        return -1;
    size_t field_size = nrows * ncharts;
    int32_t *previous_row = (int32_t *)ARENA_ALLOC(
        arena, field_size * sizeof(*previous_row));
    int32_t *next_row = (int32_t *)ARENA_ALLOC(
        arena, field_size * sizeof(*next_row));
    for (size_t chart = 0; chart < ncharts; chart++) {
        int32_t previous = -1;
        for (size_t row = 0; row < nrows; row++) {
            size_t at = row * ncharts + chart;
            previous_row[at] = previous;
            if (row_present[at] > 0) previous = (int32_t)row;
        }
        int32_t next = -1;
        for (size_t row = nrows; row-- > 0;) {
            size_t at = row * ncharts + chart;
            next_row[at] = next;
            if (row_present[at] > 0) next = (int32_t)row;
        }
    }

    size_t *row_start = (size_t *)ARENA_ALLOC(
        arena, (nrows + 1) * sizeof(*row_start));
    size_t sample_at = 0;
    for (size_t row = 0; row < nrows; row++) {
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row < (int32_t)row)
            sample_at++;
        row_start[row] = sample_at;
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row == (int32_t)row)
            sample_at++;
    }
    row_start[nrows] = sample_at;
    if (sample_at != set->nlayer_samples) return -1;

    size_t *bound_start = (size_t *)ARENA_ALLOC(
        arena, (nrows + 1) * sizeof(*bound_start));
    size_t bound_at = 0;
    for (size_t row = 0; row < nrows; row++) {
        while (bound_at < out->nregister_bounds &&
               out->register_bound[bound_at].row < (int32_t)row)
            bound_at++;
        bound_start[row] = bound_at;
        while (bound_at < out->nregister_bounds &&
               out->register_bound[bound_at].row == (int32_t)row)
            bound_at++;
    }
    bound_start[nrows] = bound_at;
    if (bound_at != out->nregister_bounds) return -1;

    size_t *collision_start = (size_t *)ARENA_ALLOC(
        arena, (nrows + 1) * sizeof(*collision_start));
    size_t collision_at = 0;
    for (size_t row = 0; row < nrows; row++) {
        while (collision_at < ncollision_bounds &&
               collision_bound[collision_at].row < (int32_t)row)
            collision_at++;
        collision_start[row] = collision_at;
        while (collision_at < ncollision_bounds &&
               collision_bound[collision_at].row == (int32_t)row) {
            const AtlasRibbonCollisionBound *b = &collision_bound[collision_at];
            if (b->chart_lo < 0 || b->chart_hi < 0 ||
                (size_t)b->chart_lo >= ncharts ||
                (size_t)b->chart_hi >= ncharts ||
                b->rank_lo >= ncharts || b->rank_hi >= ncharts ||
                b->rank_lo >= b->rank_hi || !isfinite(b->lower_shift_delta) ||
                solution->chart[b->chart_lo].rank != b->rank_lo ||
                solution->chart[b->chart_hi].rank != b->rank_hi)
                return -1;
            collision_at++;
        }
    }
    collision_start[nrows] = collision_at;
    if (collision_at != ncollision_bounds) return -1;

    double *target_sum = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(*target_sum));
    double *target_weight = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(*target_weight));
    double *before = (double *)ARENA_ALLOC(
        arena, field_size * sizeof(*before));
    size_t sweep_base = out->register_smoothing_sweeps;
    for (int sweep = 0; sweep < sweep_limit; sweep++) {
        memcpy(before, out->row_chart_shift,
               field_size * sizeof(*before));
        for (int direction = 0; direction < 2; direction++) {
            for (size_t step = 0; step < nrows; step++) {
                size_t row = direction == 0 ? step : nrows - 1 - step;
                size_t first = row_start[row], last = row_start[row + 1];
                if (first == last) continue;
                for (size_t i = 0; i < ncharts * ncharts; i++) {
                    bounds->lower[i] = -DBL_MAX;
                    local_pair[i] = 0;
                }
                memset(present, 0, ncharts * sizeof(*present));
                memset(target_sum, 0, ncharts * sizeof(*target_sum));
                memset(target_weight, 0, ncharts * sizeof(*target_weight));
                for (size_t rank = 0; rank < ncharts; rank++) {
                    size_t chart = rank_to_chart[rank];
                    size_t field_at = row * ncharts + chart;
                    present[rank] = row_present[field_at];
                    rank_shift[rank] = out->row_chart_shift[field_at];
                    if (present[rank] == 0) continue;
                    int32_t adjacent[2] = {
                        previous_row[field_at], next_row[field_at]
                    };
                    for (int side = 0; side < 2; side++) {
                        if (adjacent[side] < 0) continue;
                        size_t other = (size_t)adjacent[side];
                        double gap = fabs((double)other - (double)row);
                        double weight = opts->lambda_register_vertical / gap;
                        target_sum[rank] += weight *
                            out->row_chart_shift[other * ncharts + chart];
                        target_weight[rank] += weight;
                    }
                    if (phase_weight != NULL && row > 0) {
                        size_t edge = (row - 1) * ncharts + chart;
                        double weight = phase_weight[edge];
                        double delta = phase_delta[edge];
                        if (!isfinite(weight) || weight < 0.0 ||
                            (weight > 0.0 && !isfinite(delta)))
                            return -1;
                        if (weight > 0.0 &&
                            row_present[(row - 1) * ncharts + chart] > 0) {
                            target_sum[rank] += weight *
                                (out->row_chart_shift[
                                     (row - 1) * ncharts + chart] + delta);
                            target_weight[rank] += weight;
                        }
                    }
                    if (phase_weight != NULL && row + 1 < nrows) {
                        size_t edge = row * ncharts + chart;
                        double weight = phase_weight[edge];
                        double delta = phase_delta[edge];
                        if (!isfinite(weight) || weight < 0.0 ||
                            (weight > 0.0 && !isfinite(delta)))
                            return -1;
                        if (weight > 0.0 &&
                            row_present[(row + 1) * ncharts + chart] > 0) {
                            target_sum[rank] += weight *
                                (out->row_chart_shift[
                                     (row + 1) * ncharts + chart] - delta);
                            target_weight[rank] += weight;
                        }
                    }
                    if (anchor_weight > 0.0) {
                        double target = anchor_target[field_at];
                        if (!isfinite(target)) return -1;
                        target_sum[rank] += anchor_weight * target;
                        target_weight[rank] += anchor_weight;
                    }
                }
                for (size_t i = bound_start[row];
                     i < bound_start[row + 1]; i++) {
                    const AtlasRibbonRegisterBound *b =
                        &out->register_bound[i];
                    if (b->rank_lo >= ncharts || b->rank_hi >= ncharts ||
                        b->rank_lo >= b->rank_hi)
                        return -1;
                    bounds->lower[b->rank_lo * ncharts + b->rank_hi] =
                        b->lower_shift_delta;
                }
                size_t previous = SIZE_MAX;
                for (size_t i = collision_start[row];
                     i < collision_start[row + 1]; i++) {
                    const AtlasRibbonCollisionBound *b = &collision_bound[i];
                    size_t key = (size_t)b->rank_lo * ncharts +
                                 (size_t)b->rank_hi;
                    if (b->lower_shift_delta > bounds->lower[key])
                        bounds->lower[key] = b->lower_shift_delta;
                }
                for (size_t i = first; i < last; i++) {
                    const AtlasRibbonLayerSample *current =
                        &set->layer_sample[i];
                    if (current->cluster_count != 1) continue;
                    if (previous != SIZE_MAX) {
                        const AtlasRibbonLayerSample *prior =
                            &set->layer_sample[previous];
                        if (prior->rank != current->rank &&
                            arf_layer_samples_local(prior, current, opts)) {
                            size_t lo = prior->rank < current->rank
                                      ? prior->rank : current->rank;
                            size_t hi = prior->rank < current->rank
                                      ? current->rank : prior->rank;
                            size_t pair = lo * ncharts + hi;
                            if (local_pair[pair] != UINT32_MAX)
                                local_pair[pair]++;
                        }
                    }
                    previous = i;
                }
                int qrc = arf_register_row_qp(
                    arena, bounds, present, local_pair,
                    target_sum, target_weight, rank_shift);
                if (qrc < 0) return -1;
                if (qrc > 0) out->register_qp_failures++;
                for (size_t rank = 0; rank < ncharts; rank++) {
                    size_t chart = rank_to_chart[rank];
                    if (present[rank] > 0)
                        out->row_chart_shift[row * ncharts + chart] =
                            rank_shift[rank];
                }
            }
        }
        double max_delta = 0.0;
        for (size_t i = 0; i < field_size; i++)
            if (row_present[i] > 0 &&
                fabs(out->row_chart_shift[i] - before[i]) > max_delta)
                max_delta = fabs(out->row_chart_shift[i] - before[i]);
        out->register_smoothing_sweeps = sweep_base + (size_t)sweep + 1;
        out->register_smoothing_max_delta = max_delta;
        out->register_iterations += 2;
        if (max_delta < 1.0e-3) break;
    }
    return 0;
}

static int arf_register_refresh_u_field(
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const size_t *rank_to_chart,
    const uint32_t *row_present,
    double *rank_shift,
    ArfRegisterRef *ref,
    double *chart_weight,
    AtlasRibbonFitResult *out)
{
    size_t ncharts = solution->ncharts;
    out->register_inversions_after = 0;
    out->register_local_inversions = 0;
    out->register_unresolved = 0;
    out->register_present_chart_rows = 0;
    out->register_vertical_shift_rms = 0.0;
    out->register_vertical_shift_max = 0.0;
    out->register_local_gap_mean = 0.0;
    out->register_local_gap_max = 0.0;
    out->register_row_width_ratio_mean = 0.0;
    out->register_shift_rms = 0.0;
    out->register_shift_max = 0.0;
    memset(out->chart_shift, 0, ncharts * sizeof(*out->chart_shift));
    memset(chart_weight, 0, ncharts * sizeof(*chart_weight));

    double source_min = DBL_MAX, source_max = -DBL_MAX;
    double registered_min = DBL_MAX, registered_max = -DBL_MAX;
    double shift_sum2 = 0.0, local_gap_sum = 0.0, row_ratio_sum = 0.0;
    size_t shift_count = 0, local_gap_count = 0, row_ratio_count = 0;
    size_t sample_at = 0;
    for (size_t row = 0; row < set->nrows; row++) {
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row < (int32_t)row)
            sample_at++;
        size_t first = sample_at;
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row == (int32_t)row)
            sample_at++;
        size_t last = sample_at;
        if (first == last) continue;
        for (size_t rank = 0; rank < ncharts; rank++) {
            size_t chart = rank_to_chart[rank];
            rank_shift[rank] =
                out->row_chart_shift[row * ncharts + chart];
        }
        size_t final_local = 0, final_same = 0;
        size_t final_inversions = arf_register_scan_one_row(
            set, opts, first, last, rank_shift, NULL, ref,
            out->register_winding_direction, 0, NULL,
            &final_local, &final_same, NULL);
        if (final_inversions == SIZE_MAX) return -1;
        out->register_inversions_after += final_inversions;
        out->register_local_inversions += final_local;
        out->register_unresolved += final_same;

        double row_source_min = DBL_MAX, row_source_max = -DBL_MAX;
        double row_registered_min = DBL_MAX, row_registered_max = -DBL_MAX;
        for (size_t i = first; i < last; i++) {
            const AtlasRibbonLayerSample *sample = &set->layer_sample[i];
            double registered = sample->source_u + rank_shift[sample->rank];
            out->registered_u[i] = registered;
            if (sample->source_u < source_min) source_min = sample->source_u;
            if (sample->source_u > source_max) source_max = sample->source_u;
            if (registered < registered_min) registered_min = registered;
            if (registered > registered_max) registered_max = registered;
            if (sample->source_u < row_source_min)
                row_source_min = sample->source_u;
            if (sample->source_u > row_source_max)
                row_source_max = sample->source_u;
            if (registered < row_registered_min)
                row_registered_min = registered;
            if (registered > row_registered_max)
                row_registered_max = registered;
        }
        if (row_source_max > row_source_min) {
            row_ratio_sum += (row_registered_max - row_registered_min) /
                             (row_source_max - row_source_min);
            row_ratio_count++;
        }
        for (size_t rank = 0; rank < ncharts; rank++) {
            size_t chart = rank_to_chart[rank];
            uint32_t samples = row_present[row * ncharts + chart];
            if (samples == 0) continue;
            double shift = rank_shift[rank];
            out->register_present_chart_rows++;
            shift_sum2 += shift * shift;
            shift_count++;
            if (fabs(shift) > out->register_shift_max)
                out->register_shift_max = fabs(shift);
            out->chart_shift[chart] += shift * (double)samples;
            chart_weight[chart] += (double)samples;
        }

        for (size_t j = 0; j < last - first; j++) {
            size_t sample = first + j;
            ref[j].sample = sample;
            ref[j].rank = set->layer_sample[sample].rank;
            ref[j].u = out->registered_u[sample];
        }
        qsort(ref, last - first, sizeof ref[0], arf_compare_register_ref);
        for (size_t j = 1; j < last - first; j++) {
            const AtlasRibbonLayerSample *a =
                &set->layer_sample[ref[j - 1].sample];
            const AtlasRibbonLayerSample *b =
                &set->layer_sample[ref[j].sample];
            if (a->rank == b->rank || !arf_layer_samples_local(a, b, opts))
                continue;
            double gap = ref[j].u - ref[j - 1].u;
            local_gap_sum += gap;
            local_gap_count++;
            if (gap > out->register_local_gap_max)
                out->register_local_gap_max = gap;
        }
    }
    if (sample_at != set->nlayer_samples) return -1;
    out->register_width_before = source_max - source_min;
    out->register_width_after = registered_max - registered_min;
    out->register_shift_rms = shift_count > 0
        ? sqrt(shift_sum2 / (double)shift_count) : 0.0;
    out->register_local_gap_mean = local_gap_count > 0
        ? local_gap_sum / (double)local_gap_count : 0.0;
    out->register_row_width_ratio_mean = row_ratio_count > 0
        ? row_ratio_sum / (double)row_ratio_count : 0.0;
    for (size_t chart = 0; chart < ncharts; chart++)
        if (chart_weight[chart] > 0.0)
            out->chart_shift[chart] /= chart_weight[chart];
    for (size_t i = 0; i < out->nregister_bounds; i++) {
        AtlasRibbonRegisterBound *bound = &out->register_bound[i];
        if (bound->row < 0 || (size_t)bound->row >= set->nrows ||
            bound->chart_lo < 0 || (size_t)bound->chart_lo >= ncharts ||
            bound->chart_hi < 0 || (size_t)bound->chart_hi >= ncharts)
            return -1;
        size_t row = (size_t)bound->row;
        bound->final_shift_delta =
            out->row_chart_shift[row * ncharts + (size_t)bound->chart_hi] -
            out->row_chart_shift[row * ncharts + (size_t)bound->chart_lo];
    }
    arf_register_vertical_stats(
        set->nrows, ncharts, row_present, out->row_chart_shift,
        &out->register_vertical_shift_rms,
        &out->register_vertical_shift_max);
    return 0;
}

static int arf_register_u_field(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonFitResult *out)
{
    if (set->nlayer_samples == 0 || solution->ncharts == 0 ||
        solution->ncharts > SIZE_MAX / solution->ncharts ||
        set->nrows > SIZE_MAX / solution->ncharts)
        return -1;
    size_t ncharts = solution->ncharts;
    size_t matrix_size = ncharts * ncharts;
    size_t field_size = set->nrows * ncharts;
    size_t *rank_to_chart = (size_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(*rank_to_chart));
    uint8_t *rank_seen = (uint8_t *)ARENA_CALLOC(arena, ncharts, 1);
    for (size_t c = 0; c < ncharts; c++) {
        if (solution->chart[c].chart != (int32_t)c ||
            solution->chart[c].rank >= ncharts ||
            rank_seen[solution->chart[c].rank])
            return -1;
        rank_seen[solution->chart[c].rank] = 1;
        rank_to_chart[solution->chart[c].rank] = c;
    }
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *s = &set->layer_sample[i];
        if (s->chart < 0 || (size_t)s->chart >= ncharts ||
            s->rank != solution->chart[s->chart].rank)
            return -1;
    }

    out->nrows = set->nrows;
    out->v0 = set->v0;
    out->dv = set->dv;
    out->register_winding_direction =
        arf_registration_winding_direction(set);
    out->registered_u = (double *)ARENA_ALLOC(
        arena, set->nlayer_samples * sizeof(*out->registered_u));
    out->chart_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(*out->chart_shift));
    out->row_chart_shift = (double *)ARENA_CALLOC(
        arena, field_size, sizeof(*out->row_chart_shift));
    uint32_t *row_present = (uint32_t *)ARENA_CALLOC(
        arena, field_size, sizeof(*row_present));
    double *chart_weight = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(*chart_weight));

    ArfRegisterBounds bounds;
    memset(&bounds, 0, sizeof bounds);
    bounds.ncharts = ncharts;
    bounds.lower = (double *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.lower));
    bounds.kind = (uint8_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.kind));
    bounds.updates = (uint32_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.updates));
    bounds.sample_lo = (uint64_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.sample_lo));
    bounds.sample_hi = (uint64_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.sample_hi));
    bounds.required_arc = (double *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.required_arc));
    bounds.endpoint_distance = (double *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bounds.endpoint_distance));
    uint32_t *present = (uint32_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(*present));
    uint32_t *local_pair = (uint32_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*local_pair));
    double *rank_shift = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(*rank_shift));

    size_t max_row_samples = 0;
    for (size_t first = 0; first < set->nlayer_samples;) {
        size_t last = first + 1;
        while (last < set->nlayer_samples &&
               set->layer_sample[last].row == set->layer_sample[first].row)
            last++;
        if (last - first > max_row_samples) max_row_samples = last - first;
        first = last;
    }
    ArfRegisterRef *ref = (ArfRegisterRef *)ARENA_ALLOC(
        arena, max_row_samples * sizeof(*ref));
    size_t bound_capacity = set->nlayer_samples;
    out->register_bound = (AtlasRibbonRegisterBound *)ARENA_ALLOC(
        arena, bound_capacity * sizeof(*out->register_bound));
    size_t nstored_bounds = 0;

    double source_min = DBL_MAX, source_max = -DBL_MAX;
    double registered_min = DBL_MAX, registered_max = -DBL_MAX;
    double shift_sum2 = 0.0, local_gap_sum = 0.0;
    double row_ratio_sum = 0.0;
    size_t shift_count = 0, local_gap_count = 0, row_ratio_count = 0;
    double max_required_arc = 0.0;
    size_t sample_at = 0;
    for (size_t row_index = 0; row_index < set->nrows; row_index++) {
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row < (int32_t)row_index)
            sample_at++;
        size_t first = sample_at;
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row == (int32_t)row_index)
            sample_at++;
        size_t last = sample_at;
        if (first == last) continue;

        for (size_t i = 0; i < matrix_size; i++) {
            bounds.lower[i] = -DBL_MAX;
            bounds.kind[i] = 0;
            bounds.updates[i] = 0;
            bounds.sample_lo[i] = bounds.sample_hi[i] = UINT64_MAX;
            bounds.required_arc[i] = 0.0;
            bounds.endpoint_distance[i] = 0.0;
            local_pair[i] = 0;
        }
        memset(present, 0, ncharts * sizeof(*present));
        memset(rank_shift, 0, ncharts * sizeof(*rank_shift));
        for (size_t i = first; i < last; i++)
            present[set->layer_sample[i].rank]++;

        for (size_t key = first; key < last;) {
            size_t count = set->layer_sample[key].cluster_count;
            if (count == 0 || key + count > last) return -1;
            for (size_t j = 1; j < count; j++) {
                int rc = arf_register_add_bound(
                    &bounds, set, key + j - 1, key + j,
                    ATLAS_RIBBON_REGISTER_REMOTE,
                    out->register_winding_direction, 1, &max_required_arc);
                if (rc < 0) return -1;
            }
            key += count;
        }

        size_t previous = SIZE_MAX;
        for (size_t i = first; i < last; i++) {
            const AtlasRibbonLayerSample *current = &set->layer_sample[i];
            if (current->cluster_count != 1) continue;
            if (previous != SIZE_MAX) {
                const AtlasRibbonLayerSample *prior =
                    &set->layer_sample[previous];
                double gap = current->source_u - prior->source_u;
                double distance = arf_distance2(prior->p, current->p);
                int within_fill = opts->max_fill_u <= 0.0 ||
                                  gap <= opts->max_fill_u;
                if (prior->rank != current->rank &&
                    arf_layer_samples_local(prior, current, opts)) {
                    uint64_t lo = prior->rank < current->rank
                                ? prior->rank : current->rank;
                    uint64_t hi = prior->rank < current->rank
                                ? current->rank : prior->rank;
                    size_t at = (size_t)lo * ncharts + (size_t)hi;
                    if (local_pair[at] != UINT32_MAX) local_pair[at]++;
                } else if (within_fill && gap >= 0.0 &&
                           distance > opts->max_bridge_stretch * gap +
                                      opts->bridge_slack &&
                           prior->rank != current->rank) {
                    int rc = arf_register_add_bound(
                        &bounds, set, previous, i,
                        ATLAS_RIBBON_REGISTER_METRIC,
                        out->register_winding_direction, 1, &max_required_arc);
                    if (rc < 0) return -1;
                }
            }
            previous = i;
        }

        size_t initial_local = 0, initial_same = 0;
        size_t initial_inversions = arf_register_scan_one_row(
            set, opts, first, last, rank_shift, &bounds, ref,
            out->register_winding_direction, 0, NULL,
            &initial_local, &initial_same, &max_required_arc);
        if (initial_inversions == SIZE_MAX) return -1;
        out->register_inversions_before += initial_inversions;

        int qrc = arf_register_row_qp(
            arena, &bounds, present, local_pair, NULL, NULL, rank_shift);
        if (qrc < 0) return -1;
        if (qrc > 0) out->register_qp_failures++;
        if (out->register_iterations < 1)
            out->register_iterations = 1;

        size_t final_local = 0, final_same = 0;
        size_t final_inversions = arf_register_scan_one_row(
            set, opts, first, last, rank_shift, &bounds, ref,
            out->register_winding_direction, 0, NULL,
            &final_local, &final_same, &max_required_arc);
        if (final_inversions == SIZE_MAX) return -1;
        out->register_inversions_after += final_inversions;
        out->register_local_inversions += final_local;
        out->register_unresolved += final_same + (qrc > 0 ? 1 : 0);

        double row_source_min = DBL_MAX, row_source_max = -DBL_MAX;
        double row_registered_min = DBL_MAX, row_registered_max = -DBL_MAX;
        for (size_t i = first; i < last; i++) {
            const AtlasRibbonLayerSample *s = &set->layer_sample[i];
            double u = s->source_u + rank_shift[s->rank];
            out->registered_u[i] = u;
            if (s->source_u < source_min) source_min = s->source_u;
            if (s->source_u > source_max) source_max = s->source_u;
            if (u < registered_min) registered_min = u;
            if (u > registered_max) registered_max = u;
            if (s->source_u < row_source_min) row_source_min = s->source_u;
            if (s->source_u > row_source_max) row_source_max = s->source_u;
            if (u < row_registered_min) row_registered_min = u;
            if (u > row_registered_max) row_registered_max = u;
        }
        if (row_source_max > row_source_min) {
            row_ratio_sum += (row_registered_max - row_registered_min) /
                             (row_source_max - row_source_min);
            row_ratio_count++;
        }
        for (size_t r = 0; r < ncharts; r++) {
            size_t chart = rank_to_chart[r];
            out->row_chart_shift[row_index * ncharts + chart] = rank_shift[r];
            row_present[row_index * ncharts + chart] = present[r];
            if (present[r] == 0) continue;
            out->register_present_chart_rows++;
            shift_sum2 += rank_shift[r] * rank_shift[r];
            shift_count++;
            double magnitude = fabs(rank_shift[r]);
            if (magnitude > out->register_shift_max)
                out->register_shift_max = magnitude;
            out->chart_shift[chart] += rank_shift[r] * (double)present[r];
            chart_weight[chart] += (double)present[r];
        }

        for (size_t j = 0; j < last - first; j++) {
            size_t sample = first + j;
            ref[j].sample = sample;
            ref[j].rank = set->layer_sample[sample].rank;
            ref[j].u = out->registered_u[sample];
        }
        qsort(ref, last - first, sizeof ref[0], arf_compare_register_ref);
        for (size_t j = 1; j < last - first; j++) {
            const AtlasRibbonLayerSample *a =
                &set->layer_sample[ref[j - 1].sample];
            const AtlasRibbonLayerSample *b =
                &set->layer_sample[ref[j].sample];
            if (a->rank == b->rank ||
                !arf_layer_samples_local(a, b, opts))
                continue;
            double gap = ref[j].u - ref[j - 1].u;
            local_gap_sum += gap;
            local_gap_count++;
            if (gap > out->register_local_gap_max)
                out->register_local_gap_max = gap;
        }

        for (size_t lo = 0; lo < ncharts; lo++)
        for (size_t hi = lo + 1; hi < ncharts; hi++) {
            size_t bi = lo * ncharts + hi;
            if (!(bounds.lower[bi] > -DBL_MAX)) continue;
            if (nstored_bounds >= bound_capacity) return -1;
            AtlasRibbonRegisterBound *dst =
                &out->register_bound[nstored_bounds++];
            memset(dst, 0, sizeof *dst);
            dst->chart_lo = (int32_t)rank_to_chart[lo];
            dst->chart_hi = (int32_t)rank_to_chart[hi];
            dst->row = (int32_t)row_index;
            dst->rank_lo = lo;
            dst->rank_hi = hi;
            dst->sample_lo = bounds.sample_lo[bi];
            dst->sample_hi = bounds.sample_hi[bi];
            dst->updates = bounds.updates[bi];
            dst->kind = bounds.kind[bi];
            dst->lower_shift_delta = bounds.lower[bi];
            dst->final_shift_delta = rank_shift[hi] - rank_shift[lo];
            dst->required_arc = bounds.required_arc[bi];
            dst->endpoint_distance = bounds.endpoint_distance[bi];
            if (dst->sample_lo < set->nlayer_samples &&
                dst->sample_hi < set->nlayer_samples)
                dst->source_u_delta =
                    set->layer_sample[dst->sample_hi].source_u -
                    set->layer_sample[dst->sample_lo].source_u;
            else
                dst->source_u_delta = NAN;
            out->register_constraints++;
            if (dst->kind == ATLAS_RIBBON_REGISTER_METRIC)
                out->register_metric_constraints++;
            else if (dst->kind == ATLAS_RIBBON_REGISTER_ORDER)
                out->register_order_constraints++;
        }
    }
    out->nregister_bounds = nstored_bounds;
    arf_register_vertical_stats(
        set->nrows, ncharts, row_present, out->row_chart_shift,
        &out->register_vertical_shift_rms_before,
        &out->register_vertical_shift_max_before);
    if ((opts->mode == ATLAS_RIBBON_REGISTER_FIELD_SMOOTH_ONLY ||
         opts->mode == ATLAS_RIBBON_REGISTER_COLLISION_ONLY ||
         opts->mode == ATLAS_RIBBON_REGISTERED_RIBBON) &&
        arf_register_smooth_u_field(
            arena, solution, set, opts, rank_to_chart, row_present,
            &bounds, NULL, 0, NULL, NULL, NULL, 0.0,
            opts->register_sweeps,
            present, local_pair, rank_shift, out) != 0)
        return -1;
    out->register_width_before = source_max - source_min;
    out->register_width_after = registered_max - registered_min;
    out->register_shift_rms = shift_count > 0
        ? sqrt(shift_sum2 / (double)shift_count) : 0.0;
    out->register_local_gap_mean = local_gap_count > 0
        ? local_gap_sum / (double)local_gap_count : 0.0;
    out->register_row_width_ratio_mean = row_ratio_count > 0
        ? row_ratio_sum / (double)row_ratio_count : 0.0;
    out->register_max_required_arc = max_required_arc;
    for (size_t c = 0; c < ncharts; c++)
        if (chart_weight[c] > 0.0)
            out->chart_shift[c] /= chart_weight[c];

    double vertical_sum2 = 0.0;
    size_t vertical_count = 0;
    for (size_t row = 1; row < set->nrows; row++)
    for (size_t c = 0; c < ncharts; c++) {
        if (row_present[(row - 1) * ncharts + c] == 0 ||
            row_present[row * ncharts + c] == 0)
            continue;
        double delta = out->row_chart_shift[row * ncharts + c] -
                       out->row_chart_shift[(row - 1) * ncharts + c];
        vertical_sum2 += delta * delta;
        vertical_count++;
        if (fabs(delta) > out->register_vertical_shift_max)
            out->register_vertical_shift_max = fabs(delta);
    }
    out->register_vertical_shift_rms = vertical_count > 0
        ? sqrt(vertical_sum2 / (double)vertical_count) : 0.0;
    if (arf_register_refresh_u_field(
            solution, set, opts, rank_to_chart, row_present,
            rank_shift, ref, chart_weight, out) != 0)
        return -1;
    return 0;
}

static int arf_refine_u_field(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonCollisionBound *collision_bound,
    size_t ncollision_bounds,
    const double *phase_delta,
    const double *phase_weight,
    const double *anchor_target,
    double anchor_weight,
    int sweeps,
    AtlasRibbonFitResult *out)
{
    if (scratch == NULL || solution == NULL || set == NULL || out == NULL ||
        !arf_options_valid(opts) ||
        (ncollision_bounds > 0 && collision_bound == NULL) ||
        solution->ncharts == 0 || set->nlayer_samples == 0 ||
        set->nrows != out->nrows || out->row_chart_shift == NULL ||
        out->registered_u == NULL || out->chart_shift == NULL ||
        set->nrows > SIZE_MAX / solution->ncharts ||
        solution->ncharts > SIZE_MAX / solution->ncharts)
        return -1;
    size_t ncharts = solution->ncharts;
    size_t field_size = set->nrows * ncharts;
    size_t matrix_size = ncharts * ncharts;
    size_t *rank_to_chart = (size_t *)ARENA_ALLOC(
        scratch, ncharts * sizeof(*rank_to_chart));
    uint8_t *rank_seen = (uint8_t *)ARENA_CALLOC(scratch, ncharts, 1);
    for (size_t chart = 0; chart < ncharts; chart++) {
        uint64_t rank = solution->chart[chart].rank;
        if (solution->chart[chart].chart != (int32_t)chart ||
            rank >= ncharts || rank_seen[rank])
            return -1;
        rank_seen[rank] = 1;
        rank_to_chart[rank] = chart;
    }

    uint32_t *row_present = (uint32_t *)ARENA_CALLOC(
        scratch, field_size, sizeof(*row_present));
    size_t max_row_samples = 0, current_row_samples = 0;
    int32_t previous_row = -1;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *sample = &set->layer_sample[i];
        if (sample->row < 0 || (size_t)sample->row >= set->nrows ||
            sample->chart < 0 || (size_t)sample->chart >= ncharts ||
            sample->rank != solution->chart[sample->chart].rank)
            return -1;
        size_t at = (size_t)sample->row * ncharts + (size_t)sample->chart;
        if (row_present[at] != UINT32_MAX) row_present[at]++;
        if (sample->row != previous_row) {
            if (current_row_samples > max_row_samples)
                max_row_samples = current_row_samples;
            current_row_samples = 0;
            previous_row = sample->row;
        }
        current_row_samples++;
    }
    if (current_row_samples > max_row_samples)
        max_row_samples = current_row_samples;
    if (max_row_samples == 0) return -1;

    ArfRegisterBounds register_bounds;
    memset(&register_bounds, 0, sizeof register_bounds);
    register_bounds.ncharts = ncharts;
    register_bounds.lower = (double *)ARENA_ALLOC(
        scratch, matrix_size * sizeof(*register_bounds.lower));
    uint32_t *present = (uint32_t *)ARENA_ALLOC(
        scratch, ncharts * sizeof(*present));
    uint32_t *local_pair = (uint32_t *)ARENA_ALLOC(
        scratch, matrix_size * sizeof(*local_pair));
    double *rank_shift = (double *)ARENA_ALLOC(
        scratch, ncharts * sizeof(*rank_shift));
    double *chart_weight = (double *)ARENA_ALLOC(
        scratch, ncharts * sizeof(*chart_weight));
    ArfRegisterRef *ref = (ArfRegisterRef *)ARENA_ALLOC(
        scratch, max_row_samples * sizeof(*ref));

    if (arf_register_smooth_u_field(
            scratch, solution, set, opts, rank_to_chart, row_present,
            &register_bounds, collision_bound, ncollision_bounds,
            phase_delta, phase_weight, anchor_target, anchor_weight, sweeps,
            present, local_pair, rank_shift, out) != 0)
        return -1;
    if (arf_register_refresh_u_field(
            solution, set, opts, rank_to_chart, row_present,
            rank_shift, ref, chart_weight, out) != 0)
        return -1;
    return 0;
}

int AtlasRibbonFit_refine_u_collisions(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonCollisionBound *collision_bound,
    size_t ncollision_bounds,
    AtlasRibbonFitResult *out)
{
    if (opts == NULL ||
        arf_refine_u_field(
            scratch, solution, set, opts, collision_bound, ncollision_bounds,
            NULL, NULL, NULL, 0.0, opts->collision_sweeps, out) != 0)
        return -1;
    size_t face_pairs = 0;
    for (size_t i = 0; i < ncollision_bounds; i++) {
        size_t count = collision_bound[i].face_pairs;
        face_pairs = count > SIZE_MAX - face_pairs ? SIZE_MAX
                                                   : face_pairs + count;
    }
    out->register_collision_rounds++;
    out->register_collision_constraints = ncollision_bounds;
    out->register_collision_face_pairs = face_pairs;
    return 0;
}

int AtlasRibbonFit_refine_u_deltas(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const AtlasRibbonCollisionBound *collision_bound,
    size_t ncollision_bounds,
    const double *edge_delta,
    const double *edge_weight,
    const double *anchor_target,
    double anchor_weight,
    int sweeps,
    AtlasRibbonFitResult *out)
{
    if (edge_delta == NULL || edge_weight == NULL || sweeps < 1)
        return -1;
    return arf_refine_u_field(
        scratch, solution, set, opts, collision_bound, ncollision_bounds,
        edge_delta, edge_weight, anchor_target, anchor_weight, sweeps, out);
}

int AtlasRibbonFit_refresh_u_field(
    Arena_T scratch,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonFitResult *out)
{
    return arf_refine_u_field(
        scratch, solution, set, opts, NULL, 0,
        NULL, NULL, NULL, 0.0, 0, out);
}

static double *arf_band_at(ArfBand *band, int row, int column)
{
    if (row < column) {
        int swap = row;
        row = column;
        column = swap;
    }
    int difference = row - column;
    if (row < 0 || row >= band->n || column < 0 ||
        difference > band->bandwidth)
        return NULL;
    return &band->lower[
        (size_t)row * (size_t)(band->bandwidth + 1) +
        (size_t)difference];
}

static int arf_add_row(ArfBand *band,
                       double *rhs0,
                       double *rhs1,
                       const int *variable,
                       const double *coefficient,
                       int count,
                       double target0,
                       double target1,
                       double weight)
{
    if (!(weight >= 0.0) || !isfinite(weight)) return -1;
    if (weight == 0.0) return 0;
    for (int i = 0; i < count; i++) {
        int vi = variable[i];
        if (vi < 0 || vi >= band->n || !isfinite(coefficient[i]))
            return -1;
        rhs0[vi] += weight * coefficient[i] * target0;
        rhs1[vi] += weight * coefficient[i] * target1;
        for (int j = 0; j <= i; j++) {
            double *entry = arf_band_at(band, vi, variable[j]);
            if (entry == NULL) return -1;
            *entry += weight * coefficient[i] * coefficient[j];
        }
    }
    return 0;
}

static int arf_band_factor(ArfBand *band)
{
    int n = band->n;
    int bw = band->bandwidth;
    for (int i = 0; i < n; i++) {
        int j0 = i - bw;
        if (j0 < 0) j0 = 0;
        for (int j = j0; j <= i; j++) {
            double *aij = arf_band_at(band, i, j);
            double sum = *aij;
            int k0 = i - bw;
            if (k0 < j - bw) k0 = j - bw;
            if (k0 < 0) k0 = 0;
            for (int k = k0; k < j; k++) {
                double *lik = arf_band_at(band, i, k);
                double *ljk = arf_band_at(band, j, k);
                if (lik != NULL && ljk != NULL) sum -= *lik * *ljk;
            }
            if (i == j) {
                if (!(sum > 1.0e-18) || !isfinite(sum)) return -1;
                *aij = sqrt(sum);
            } else {
                double *ljj = arf_band_at(band, j, j);
                *aij = sum / *ljj;
            }
        }
    }
    return 0;
}

static void arf_band_solve(const ArfBand *band, double *rhs)
{
    int n = band->n;
    int bw = band->bandwidth;
    for (int i = 0; i < n; i++) {
        int j0 = i - bw;
        if (j0 < 0) j0 = 0;
        double sum = rhs[i];
        for (int j = j0; j < i; j++)
            sum -= band->lower[
                (size_t)i * (size_t)(bw + 1) + (size_t)(i - j)] * rhs[j];
        rhs[i] = sum / band->lower[(size_t)i * (size_t)(bw + 1)];
    }
    for (int i = n - 1; i >= 0; i--) {
        int j1 = i + bw;
        if (j1 >= n) j1 = n - 1;
        double sum = rhs[i];
        for (int j = i + 1; j <= j1; j++)
            sum -= band->lower[
                (size_t)j * (size_t)(bw + 1) + (size_t)(j - i)] * rhs[j];
        rhs[i] = sum / band->lower[(size_t)i * (size_t)(bw + 1)];
    }
}

static int arf_solve_segment(
    Arena_T arena,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    const size_t *target_index,
    size_t first,
    size_t last,
    AtlasRibbonFitResult *out,
    size_t row)
{
    double u_first = set->target[target_index[first]].u;
    double u_last = set->target[target_index[last - 1]].u;
    int64_t global0 = (int64_t)floor((u_first - out->u0) / out->du) - 1;
    int64_t global1 = (int64_t)ceil((u_last - out->u0) / out->du) + 1;
    if (global0 < 0) global0 = 0;
    if (global1 >= (int64_t)out->ncolumns)
        global1 = (int64_t)out->ncolumns - 1;
    int n = (int)(global1 - global0 + 1);
    if (n < 2 || n > INT_MAX / 3) return -1;

    Arena_Mark mark = Arena_save(arena);
    ArfBand band;
    band.n = n;
    band.bandwidth = 2;
    band.lower = (double *)ARENA_CALLOC(
        arena, (size_t)n * 3, sizeof(double));
    double *rhs0 = (double *)ARENA_CALLOC(arena, n, sizeof(double));
    double *rhs1 = (double *)ARENA_CALLOC(arena, n, sizeof(double));

    for (size_t q = first; q < last; q++) {
        const AtlasRibbonTarget *target =
            &set->target[target_index[q]];
        double coordinate = (target->u - out->u0) / out->du;
        int64_t left_global = (int64_t)floor(coordinate);
        double t = coordinate - (double)left_global;
        if (left_global < global0) {
            left_global = global0;
            t = 0.0;
        }
        if (left_global >= global1) {
            left_global = global1 - 1;
            t = 1.0;
        }
        int variable[2] = {
            (int)(left_global - global0),
            (int)(left_global + 1 - global0)
        };
        double coefficient[2] = {1.0 - t, t};
        if (arf_add_row(&band, rhs0, rhs1, variable, coefficient, 2,
                        target->p[0], target->p[1],
                        opts->lambda_position) != 0) {
            Arena_restore(arena, mark);
            return -1;
        }
        if (opts->mode >= ATLAS_RIBBON_TANGENT &&
            opts->lambda_tangent > 0.0) {
            double deriv_coeff[2] = {-1.0, 1.0};
            if (arf_add_row(
                    &band, rhs0, rhs1, variable, deriv_coeff, 2,
                    out->du * target->tangent[0],
                    out->du * target->tangent[1],
                    opts->lambda_tangent) != 0) {
                Arena_restore(arena, mark);
                return -1;
            }
        }
    }
    if (opts->lambda_smooth > 0.0) {
        for (int i = 1; i + 1 < n; i++) {
            int variable[3] = {i - 1, i, i + 1};
            double coefficient[3] = {1.0, -2.0, 1.0};
            if (arf_add_row(&band, rhs0, rhs1, variable, coefficient, 3,
                            0.0, 0.0, opts->lambda_smooth) != 0) {
                Arena_restore(arena, mark);
                return -1;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double *diagonal = arf_band_at(&band, i, i);
        *diagonal += 1.0e-9;
    }
    if (arf_band_factor(&band) != 0) {
        Arena_restore(arena, mark);
        return -1;
    }
    arf_band_solve(&band, rhs0);
    arf_band_solve(&band, rhs1);
    for (int i = 0; i < n; i++) {
        size_t column = (size_t)(global0 + i);
        size_t index = (row * out->ncolumns + column);
        out->p[index * 2] = rhs0[i];
        out->p[index * 2 + 1] = rhs1[i];
        if (!out->valid[index]) {
            out->valid[index] = 1;
            out->fitted_nodes++;
        }
    }
    size_t nearest = first;
    for (int i = 0; i < n; i++) {
        double u = out->u0 + (double)(global0 + i) * out->du;
        while (nearest + 1 < last &&
               fabs(set->target[target_index[nearest + 1]].u - u) <
               fabs(set->target[target_index[nearest]].u - u))
            nearest++;
        if (fabs(set->target[target_index[nearest]].u - u) <= out->du) {
            size_t column = (size_t)(global0 + i);
            size_t index = row * out->ncolumns + column;
            if (!out->support[index]) {
                out->support[index] = 1;
                out->supported_nodes++;
            }
        }
    }
    Arena_restore(arena, mark);
    return 0;
}

static int arf_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa;
    double b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double arf_quantile_sorted(const double *x, size_t n, double q)
{
    if (n == 0) return NAN;
    double position = q * (double)(n - 1);
    size_t lo = (size_t)floor(position);
    size_t hi = (size_t)ceil(position);
    double t = position - (double)lo;
    return (1.0 - t) * x[lo] + t * x[hi];
}

static void arf_distribution(double *value, size_t count,
                             double *mean, double *rms,
                             double *p50, double *p90, double *p95,
                             double *p99, double *maximum)
{
    if (count == 0) {
        *mean = *rms = *p50 = *p90 = *p95 = *p99 = *maximum = NAN;
        return;
    }
    double sum = 0.0, sum2 = 0.0, maxv = 0.0;
    for (size_t i = 0; i < count; i++) {
        sum += value[i];
        sum2 += value[i] * value[i];
        if (value[i] > maxv) maxv = value[i];
    }
    qsort(value, count, sizeof(double), arf_compare_double);
    *mean = sum / (double)count;
    *rms = sqrt(sum2 / (double)count);
    *p50 = arf_quantile_sorted(value, count, 0.50);
    *p90 = arf_quantile_sorted(value, count, 0.90);
    *p95 = arf_quantile_sorted(value, count, 0.95);
    *p99 = arf_quantile_sorted(value, count, 0.99);
    *maximum = maxv;
}

static int arf_eval_row(const AtlasRibbonFitResult *fit,
                        size_t row, double u, double p[2])
{
    double coordinate = (u - fit->u0) / fit->du;
    int64_t left = (int64_t)floor(coordinate);
    double t = coordinate - (double)left;
    if (left < 0 || left + 1 >= (int64_t)fit->ncolumns) return -1;
    size_t a = row * fit->ncolumns + (size_t)left;
    size_t b = a + 1;
    if (!fit->valid[a] || !fit->valid[b]) return -1;
    if (fit->u_edge != NULL && !fit->u_edge[a]) return -1;
    p[0] = (1.0 - t) * fit->p[a * 2] + t * fit->p[b * 2];
    p[1] = (1.0 - t) * fit->p[a * 2 + 1] + t * fit->p[b * 2 + 1];
    return 0;
}

static double arf_orient2(const double a[2],
                          const double b[2],
                          const double c[2])
{
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

static int arf_segments_cross(const ArfSegmentBox *a,
                              const ArfSegmentBox *b)
{
    double aa = arf_orient2(a->a, a->b, b->a);
    double ab = arf_orient2(a->a, a->b, b->b);
    double ba = arf_orient2(b->a, b->b, a->a);
    double bb = arf_orient2(b->a, b->b, a->b);
    double eps = 1.0e-10;
    return ((aa > eps && ab < -eps) || (aa < -eps && ab > eps)) &&
           ((ba > eps && bb < -eps) || (ba < -eps && bb > eps));
}

static int arf_compare_segment_box(const void *pa, const void *pb)
{
    const ArfSegmentBox *a = (const ArfSegmentBox *)pa;
    const ArfSegmentBox *b = (const ArfSegmentBox *)pb;
    if (a->min_x != b->min_x) return a->min_x < b->min_x ? -1 : 1;
    return a->index < b->index ? -1 : (a->index > b->index ? 1 : 0);
}

static size_t arf_row_crossings(Arena_T arena,
                                const AtlasRibbonFitResult *fit,
                                size_t row)
{
    size_t count = 0;
    for (size_t i = 0; i + 1 < fit->ncolumns; i++) {
        size_t a = row * fit->ncolumns + i;
        if (fit->valid[a] && fit->valid[a + 1] &&
            (fit->u_edge == NULL || fit->u_edge[a])) count++;
    }
    if (count < 2) return 0;
    Arena_Mark mark = Arena_save(arena);
    ArfSegmentBox *segment = (ArfSegmentBox *)ARENA_ALLOC(
        arena, count * sizeof(ArfSegmentBox));
    size_t n = 0;
    double max_extent = 0.0;
    for (size_t i = 0; i + 1 < fit->ncolumns; i++) {
        size_t a = row * fit->ncolumns + i;
        if (!fit->valid[a] || !fit->valid[a + 1] ||
            (fit->u_edge != NULL && !fit->u_edge[a])) continue;
        ArfSegmentBox *s = &segment[n++];
        s->index = (int32_t)i;
        s->a[0] = fit->p[a * 2];
        s->a[1] = fit->p[a * 2 + 1];
        s->b[0] = fit->p[(a + 1) * 2];
        s->b[1] = fit->p[(a + 1) * 2 + 1];
        s->min_x = s->a[0] < s->b[0] ? s->a[0] : s->b[0];
        s->max_x = s->a[0] > s->b[0] ? s->a[0] : s->b[0];
        s->min_y = s->a[1] < s->b[1] ? s->a[1] : s->b[1];
        s->max_y = s->a[1] > s->b[1] ? s->a[1] : s->b[1];
        double extent = s->max_x - s->min_x;
        if (extent > max_extent) max_extent = extent;
    }
    qsort(segment, n, sizeof segment[0], arf_compare_segment_box);
    size_t crossings = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t jj = i; jj-- > 0;) {
            if (segment[i].min_x - segment[jj].min_x > max_extent)
                break;
            int64_t delta = (int64_t)segment[i].index -
                            (int64_t)segment[jj].index;
            if (delta >= -1 && delta <= 1) continue;
            if (segment[jj].max_x < segment[i].min_x ||
                segment[jj].max_y < segment[i].min_y ||
                segment[i].max_y < segment[jj].min_y)
                continue;
            if (arf_segments_cross(&segment[i], &segment[jj]))
                crossings++;
        }
    }
    Arena_restore(arena, mark);
    return crossings;
}

/* A full ribbon must not preserve an invented connection through a crossing.
 * Remove the less-supported (then less-direct, then more distorted) segment of
 * each strict row intersection.  The node positions stay locked; this merely
 * turns ambiguous synthesized geometry into an explicit seam. */
static size_t arf_cut_row_crossings(Arena_T arena,
                                    AtlasRibbonFitResult *fit,
                                    size_t row)
{
    if (fit->u_edge == NULL) return 0;
    size_t cuts = 0;
    for (;;) {
        size_t count = 0;
        for (size_t i = 0; i + 1 < fit->ncolumns; i++) {
            size_t node = row * fit->ncolumns + i;
            if (fit->valid[node] && fit->valid[node + 1] &&
                fit->u_edge[node])
                count++;
        }
        if (count < 2) break;
        Arena_Mark mark = Arena_save(arena);
        ArfSegmentBox *segment = (ArfSegmentBox *)ARENA_ALLOC(
            arena, count * sizeof(*segment));
        size_t n = 0;
        double max_extent = 0.0;
        for (size_t i = 0; i + 1 < fit->ncolumns; i++) {
            size_t node = row * fit->ncolumns + i;
            if (!fit->valid[node] || !fit->valid[node + 1] ||
                !fit->u_edge[node])
                continue;
            ArfSegmentBox *s = &segment[n++];
            s->index = (int32_t)i;
            s->a[0] = fit->p[node * 2];
            s->a[1] = fit->p[node * 2 + 1];
            s->b[0] = fit->p[(node + 1) * 2];
            s->b[1] = fit->p[(node + 1) * 2 + 1];
            s->min_x = fmin(s->a[0], s->b[0]);
            s->max_x = fmax(s->a[0], s->b[0]);
            s->min_y = fmin(s->a[1], s->b[1]);
            s->max_y = fmax(s->a[1], s->b[1]);
            double extent = s->max_x - s->min_x;
            if (extent > max_extent) max_extent = extent;
        }
        qsort(segment, n, sizeof(*segment), arf_compare_segment_box);
        int32_t index0 = -1, index1 = -1;
        for (size_t i = 0; i < n && index0 < 0; i++) {
            for (size_t jj = i; jj-- > 0;) {
                if (segment[i].min_x - segment[jj].min_x > max_extent)
                    break;
                int64_t delta = (int64_t)segment[i].index -
                                (int64_t)segment[jj].index;
                if (delta >= -1 && delta <= 1) continue;
                if (segment[jj].max_x < segment[i].min_x ||
                    segment[jj].max_y < segment[i].min_y ||
                    segment[i].max_y < segment[jj].min_y)
                    continue;
                if (arf_segments_cross(&segment[i], &segment[jj])) {
                    index0 = segment[i].index;
                    index1 = segment[jj].index;
                    break;
                }
            }
        }
        if (index0 < 0 || index1 < 0) {
            Arena_restore(arena, mark);
            break;
        }
        size_t node0 = row * fit->ncolumns + (size_t)index0;
        size_t node1 = row * fit->ncolumns + (size_t)index1;
        int support0 = (int)fit->support[node0] + (int)fit->support[node0 + 1];
        int support1 = (int)fit->support[node1] + (int)fit->support[node1 + 1];
        int provenance0 = fit->provenance[node0] > fit->provenance[node0 + 1]
                        ? fit->provenance[node0] : fit->provenance[node0 + 1];
        int provenance1 = fit->provenance[node1] > fit->provenance[node1 + 1]
                        ? fit->provenance[node1] : fit->provenance[node1 + 1];
        double length0 = hypot(fit->p[(node0 + 1) * 2] - fit->p[node0 * 2],
                               fit->p[(node0 + 1) * 2 + 1] - fit->p[node0 * 2 + 1]);
        double length1 = hypot(fit->p[(node1 + 1) * 2] - fit->p[node1 * 2],
                               fit->p[(node1 + 1) * 2 + 1] - fit->p[node1 * 2 + 1]);
        double error0 = fabs(length0 / fit->du - 1.0);
        double error1 = fabs(length1 / fit->du - 1.0);
        size_t cut = support0 != support1 ? (support0 < support1 ? node0 : node1)
                   : provenance0 != provenance1
                         ? (provenance0 > provenance1 ? node0 : node1)
                   : error0 != error1 ? (error0 > error1 ? node0 : node1)
                   : (index0 > index1 ? node0 : node1);
        Arena_restore(arena, mark);
        fit->u_edge[cut] = 0;
        if (fit->u_reject != NULL)
            fit->u_reject[cut] = ATLAS_RIBBON_EDGE_CROSSING;
        cuts++;
    }
    return cuts;
}

static void arf_measure_row_winding(Arena_T arena,
                                    const AtlasRibbonFitResult *fit,
                                    size_t row,
                                    AtlasRibbonRowMetric *metric)
{
    size_t best_first = 0, best_count = 0;
    for (size_t first = 0; first < fit->ncolumns;) {
        while (first < fit->ncolumns &&
               !fit->valid[row * fit->ncolumns + first]) first++;
        if (first >= fit->ncolumns) break;
        size_t last = first + 1;
        while (last < fit->ncolumns &&
               fit->valid[row * fit->ncolumns + last] &&
               (fit->u_edge == NULL ||
                fit->u_edge[row * fit->ncolumns + last - 1]))
            last++;
        if (last - first > best_count) {
            best_first = first;
            best_count = last - first;
        }
        first = last;
    }
    if (best_count < 4) return;

    Arena_Mark mark = Arena_save(arena);
    double *wind = (double *)ARENA_ALLOC(arena, best_count * sizeof(double));
    double *radius = (double *)ARENA_ALLOC(
        arena, best_count * sizeof(double));
    double previous_angle = 0.0, unwrapped = 0.0;
    for (size_t i = 0; i < best_count; i++) {
        size_t index = row * fit->ncolumns + best_first + i;
        double x = fit->p[index * 2];
        double y = fit->p[index * 2 + 1];
        radius[i] = sqrt(x * x + y * y);
        double angle = atan2(y, x);
        if (i == 0) previous_angle = angle;
        else {
            double delta = angle - previous_angle;
            while (delta > M_PI) delta -= 2.0 * M_PI;
            while (delta < -M_PI) delta += 2.0 * M_PI;
            unwrapped += delta;
            previous_angle = angle;
        }
        wind[i] = unwrapped / (2.0 * M_PI);
    }
    double direction = wind[best_count - 1] >= 0.0 ? 1.0 : -1.0;
    size_t backsteps = 0;
    double max_backstep = 0.0;
    for (size_t i = 0; i < best_count; i++) wind[i] *= direction;
    for (size_t i = 1; i < best_count; i++) {
        double delta = wind[i] - wind[i - 1];
        if (delta < -1.0e-5) {
            backsteps++;
            if (-delta > max_backstep) max_backstep = -delta;
        }
    }
    metric->phase_backstep_fraction =
        (double)backsteps / (double)(best_count - 1);
    metric->phase_max_backstep = max_backstep;
    metric->winding_span = wind[best_count - 1] - wind[0];
    if (metric->winding_span < 0.0) metric->winding_span = 0.0;

    int nturn = (int)floor(metric->winding_span);
    if (nturn >= 2) {
        double *median = (double *)ARENA_ALLOC(
            arena, (size_t)(nturn + 1) * sizeof(double));
        int valid_turns = 0;
        for (int turn = 0; turn <= nturn; turn++) {
            size_t count = 0;
            for (size_t i = 0; i < best_count; i++)
                if ((int)floor(wind[i]) == turn) count++;
            if (count == 0) {
                median[turn] = NAN;
                continue;
            }
            double *sample = (double *)ARENA_ALLOC(
                arena, count * sizeof(double));
            size_t at = 0;
            for (size_t i = 0; i < best_count; i++)
                if ((int)floor(wind[i]) == turn) sample[at++] = radius[i];
            qsort(sample, count, sizeof(double), arf_compare_double);
            median[turn] = arf_quantile_sorted(sample, count, 0.5);
            valid_turns++;
        }
        if (valid_turns >= 2) {
            double trend = median[nturn] - median[0];
            double sign = trend >= 0.0 ? 1.0 : -1.0;
            double *pitch = (double *)ARENA_ALLOC(
                arena, (size_t)nturn * sizeof(double));
            size_t npitch = 0;
            for (int turn = 0; turn < nturn; turn++) {
                if (!isfinite(median[turn]) ||
                    !isfinite(median[turn + 1]))
                    continue;
                double delta = median[turn + 1] - median[turn];
                metric->radius_order_tests++;
                if (sign * delta <= 0.0)
                    metric->radius_order_violations++;
                pitch[npitch++] = fabs(delta);
            }
            if (npitch > 0) {
                qsort(pitch, npitch, sizeof(double), arf_compare_double);
                metric->radius_pitch_median =
                    arf_quantile_sorted(pitch, npitch, 0.5);
            }
        }
    }
    Arena_restore(arena, mark);
}

static double arf_point_segment_distance2(const double p[2],
                                          const double a[2],
                                          const double b[2])
{
    double dx = b[0] - a[0];
    double dy = b[1] - a[1];
    double den = dx * dx + dy * dy;
    double t = den > 1.0e-30
             ? ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / den
             : 0.0;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double x = a[0] + t * dx - p[0];
    double y = a[1] + t * dy - p[1];
    return sqrt(x * x + y * y);
}

static int arf_interpolated_node(const AtlasRibbonFitResult *fit,
                                 int64_t row0, int64_t row1,
                                 double t, size_t column, double p[2])
{
    if (row0 < 0 || row1 < 0 ||
        row0 >= (int64_t)fit->nrows ||
        row1 >= (int64_t)fit->nrows)
        return -1;
    size_t a = (size_t)row0 * fit->ncolumns + column;
    size_t b = (size_t)row1 * fit->ncolumns + column;
    if (!fit->valid[a] || !fit->valid[b]) return -1;
    if (fit->v_edge != NULL && row0 != row1) {
        int64_t lower = row0 < row1 ? row0 : row1;
        int64_t upper = row0 < row1 ? row1 : row0;
        if (upper != lower + 1) return -1;
        size_t vertical = (size_t)lower * fit->ncolumns + column;
        if (!fit->v_edge[vertical]) return -1;
    }
    p[0] = (1.0 - t) * fit->p[a * 2] + t * fit->p[b * 2];
    p[1] = (1.0 - t) * fit->p[a * 2 + 1] + t * fit->p[b * 2 + 1];
    return 0;
}

static void arf_measure_vertex_distances(
    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    AtlasRibbonFitResult *fit,
    const double *dense_chart_shift, size_t ncharts)
{
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, ps->nv, 1);
    for (size_t f = 0; f < ps->nf; f++) {
        if (!solution->face_keep[f]) continue;
        used[ps->faces[f * 3]] = 1;
        used[ps->faces[f * 3 + 1]] = 1;
        used[ps->faces[f * 3 + 2]] = 1;
    }
    double *chamfer = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    double *parameter = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    size_t nchamfer = 0, nparameter = 0;
    for (size_t vertex = 0; vertex < ps->nv; vertex++) {
        if (!used[vertex]) continue;
        double point[2];
        arf_project_point(set, &ps->verts[vertex * 3], point);
        double row_coordinate =
            (solution->v[vertex] - fit->v0) / fit->dv;
        int64_t row0 = (int64_t)floor(row_coordinate);
        int64_t row1 = row0 + 1;
        double row_t = row_coordinate - (double)row0;
        if (row0 < 0) {
            row0 = row1 = 0;
            row_t = 0.0;
        } else if (row1 >= (int64_t)fit->nrows) {
            row0 = row1 = (int64_t)fit->nrows - 1;
            row_t = 0.0;
        }
        double best = DBL_MAX;
        for (size_t column = 0; column + 1 < fit->ncolumns; column++) {
            if (fit->u_edge != NULL &&
                (!fit->u_edge[(size_t)row0 * fit->ncolumns + column] ||
                 !fit->u_edge[(size_t)row1 * fit->ncolumns + column]))
                continue;
            double a[2], b[2];
            if (arf_interpolated_node(fit, row0, row1, row_t,
                                      column, a) != 0 ||
                arf_interpolated_node(fit, row0, row1, row_t,
                                      column + 1, b) != 0)
                continue;
            double distance = arf_point_segment_distance2(point, a, b);
            if (distance < best) best = distance;
        }
        if (best < DBL_MAX) chamfer[nchamfer++] = best;

        double sample_u = solution->u[vertex];
        if (dense_chart_shift != NULL) {
            int32_t chart = solution->vertex_chart[vertex];
            if (chart < 0 || (size_t)chart >= ncharts) continue;
            double shift0 = dense_chart_shift[
                (size_t)row0 * ncharts + (size_t)chart];
            double shift1 = dense_chart_shift[
                (size_t)row1 * ncharts + (size_t)chart];
            sample_u += (1.0 - row_t) * shift0 + row_t * shift1;
        }
        double u_coordinate = (sample_u - fit->u0) / fit->du;
        int64_t column0 = (int64_t)floor(u_coordinate);
        double column_t = u_coordinate - (double)column0;
        int horizontal_ok = fit->u_edge == NULL ||
            (column0 >= 0 && column0 + 1 < (int64_t)fit->ncolumns &&
             fit->u_edge[(size_t)row0 * fit->ncolumns + (size_t)column0] &&
             fit->u_edge[(size_t)row1 * fit->ncolumns + (size_t)column0]);
        if (column0 >= 0 &&
            column0 + 1 < (int64_t)fit->ncolumns && horizontal_ok) {
            double a0[2], a1[2], b0[2], b1[2];
            if (arf_interpolated_node(fit, row0, row0, 0.0,
                                      (size_t)column0, a0) == 0 &&
                arf_interpolated_node(fit, row0, row0, 0.0,
                                      (size_t)column0 + 1, a1) == 0 &&
                arf_interpolated_node(fit, row1, row1, 0.0,
                                      (size_t)column0, b0) == 0 &&
                arf_interpolated_node(fit, row1, row1, 0.0,
                                      (size_t)column0 + 1, b1) == 0) {
                double pa[2] = {
                    (1.0 - column_t) * a0[0] + column_t * a1[0],
                    (1.0 - column_t) * a0[1] + column_t * a1[1]
                };
                double pb[2] = {
                    (1.0 - column_t) * b0[0] + column_t * b1[0],
                    (1.0 - column_t) * b0[1] + column_t * b1[1]
                };
                double predicted[2] = {
                    (1.0 - row_t) * pa[0] + row_t * pb[0],
                    (1.0 - row_t) * pa[1] + row_t * pb[1]
                };
                parameter[nparameter++] =
                    arf_distance2(point, predicted);
            }
        }
    }
    fit->chamfer_samples = nchamfer;
    arf_distribution(chamfer, nchamfer,
                     &fit->chamfer_mean, &fit->chamfer_rms,
                     &fit->chamfer_p50, &fit->chamfer_p90,
                     &fit->chamfer_p95, &fit->chamfer_p99,
                     &fit->directed_hausdorff);
    fit->parameter_samples = nparameter;
    arf_distribution(parameter, nparameter,
                     &fit->parameter_mean, &fit->parameter_rms,
                     &fit->parameter_p50, &fit->parameter_p90,
                     &fit->parameter_p95, &fit->parameter_p99,
                     &fit->parameter_max);
}

static int arf_build_dense_chart_shift(
    Arena_T arena,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitResult *fit,
    double **out_dense)
{
    if (arena == NULL || solution == NULL || set == NULL || fit == NULL ||
        out_dense == NULL || solution->ncharts == 0 ||
        fit->row_chart_shift == NULL || fit->nrows != set->nrows ||
        fit->nrows > SIZE_MAX / solution->ncharts)
        return -1;
    size_t ncharts = solution->ncharts;
    size_t n = fit->nrows * ncharts;
    uint8_t *present = (uint8_t *)ARENA_CALLOC(arena, n, 1);
    double *dense = (double *)ARENA_CALLOC(arena, n, sizeof(*dense));
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *sample = &set->layer_sample[i];
        if (sample->row < 0 || (size_t)sample->row >= fit->nrows ||
            sample->chart < 0 || (size_t)sample->chart >= ncharts)
            return -1;
        present[(size_t)sample->row * ncharts + (size_t)sample->chart] = 1;
    }
    for (size_t chart = 0; chart < ncharts; chart++) {
        size_t first = 0;
        while (first < fit->nrows && !present[first * ncharts + chart])
            first++;
        if (first == fit->nrows) continue;
        double first_shift = fit->row_chart_shift[first * ncharts + chart];
        for (size_t row = 0; row <= first; row++)
            dense[row * ncharts + chart] = first_shift;
        size_t previous = first;
        for (size_t next = first + 1; next < fit->nrows; next++) {
            if (!present[next * ncharts + chart]) continue;
            double a = fit->row_chart_shift[previous * ncharts + chart];
            double b = fit->row_chart_shift[next * ncharts + chart];
            double span = (double)(next - previous);
            for (size_t row = previous; row <= next; row++) {
                double t = (double)(row - previous) / span;
                dense[row * ncharts + chart] = (1.0 - t) * a + t * b;
            }
            previous = next;
        }
        double last_shift = fit->row_chart_shift[previous * ncharts + chart];
        for (size_t row = previous; row < fit->nrows; row++)
            dense[row * ncharts + chart] = last_shift;
    }
    *out_dense = dense;
    return 0;
}

static int arf_registered_spiral_span(
    const AtlasRibbonLayerSample *a,
    const AtlasRibbonLayerSample *b,
    int direction,
    double *out_angle_delta,
    double *out_arc)
{
    if (a == NULL || b == NULL || out_angle_delta == NULL || out_arc == NULL ||
        (direction != -1 && direction != 1))
        return -1;
    double radius0 = sqrt(a->p[0] * a->p[0] + a->p[1] * a->p[1]);
    double radius1 = sqrt(b->p[0] * b->p[0] + b->p[1] * b->p[1]);
    if (!(radius0 > 1.0e-6) || !(radius1 > 1.0e-6)) return -1;
    double theta0 = atan2(a->p[1], a->p[0]);
    double theta1 = atan2(b->p[1], b->p[0]);
    double alpha0 = arf_mod_tau((double)direction * theta0);
    double alpha1 = arf_mod_tau((double)direction * theta1);
    int64_t winding_delta = (int64_t)b->winding - (int64_t)a->winding;
    if (winding_delta < 0 || winding_delta >= INT32_MAX) return -1;
    double angle_delta;
    if (winding_delta > 0) {
        angle_delta = alpha1 - alpha0 +
                      2.0 * M_PI * (double)winding_delta;
        while (angle_delta <= 1.0e-12) angle_delta += 2.0 * M_PI;
    } else {
        angle_delta = atan2(sin(alpha1 - alpha0), cos(alpha1 - alpha0));
        if (angle_delta < -1.0e-8) return -1;
        if (angle_delta < 0.0) angle_delta = 0.0;
    }
    double arc = arf_spiral_arc(radius0, radius1, angle_delta);
    if (!isfinite(arc) || !isfinite(angle_delta)) return -1;
    *out_angle_delta = angle_delta;
    *out_arc = arc;
    return 0;
}

static void arf_registered_span_point(
    const AtlasRibbonLayerSample *a,
    const AtlasRibbonLayerSample *b,
    int direction,
    int kind,
    double angle_delta,
    double t,
    double p[2])
{
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    if (kind == 1) {
        p[0] = (1.0 - t) * a->p[0] + t * b->p[0];
        p[1] = (1.0 - t) * a->p[1] + t * b->p[1];
        return;
    }
    double radius0 = sqrt(a->p[0] * a->p[0] + a->p[1] * a->p[1]);
    double radius1 = sqrt(b->p[0] * b->p[0] + b->p[1] * b->p[1]);
    double radius = (1.0 - t) * radius0 + t * radius1;
    double theta = atan2(a->p[1], a->p[0]) +
                   (double)direction * angle_delta * t;
    p[0] = radius * cos(theta);
    p[1] = radius * sin(theta);
}

/* Mark every fitted U cell whose open interval intersects [u0,u1].  A valid
 * edge clears this mark later; until then the numerically larger reason is
 * the more specific/decisive rejection provenance. */
static void arf_mark_u_reject(AtlasRibbonFitResult *fit, size_t row,
                              double u0, double u1,
                              AtlasRibbonEdgeReject reason)
{
    if (fit == NULL || fit->u_reject == NULL || row >= fit->nrows ||
        fit->ncolumns < 2 || !(fit->du > 0.0) ||
        reason == ATLAS_RIBBON_EDGE_NONE)
        return;
    if (u1 < u0) {
        double swap = u0;
        u0 = u1;
        u1 = swap;
    }
    if (!(u1 > u0)) return;
    int64_t first = (int64_t)floor((u0 - fit->u0) / fit->du);
    int64_t last = (int64_t)ceil((u1 - fit->u0) / fit->du) - 1;
    if (last < 0 || first >= (int64_t)fit->ncolumns - 1) return;
    if (first < 0) first = 0;
    if (last >= (int64_t)fit->ncolumns - 1)
        last = (int64_t)fit->ncolumns - 2;
    for (int64_t column = first; column <= last; column++) {
        size_t node = row * fit->ncolumns + (size_t)column;
        if (fit->u_reject[node] < (uint8_t)reason)
            fit->u_reject[node] = (uint8_t)reason;
    }
}

int AtlasRibbonFit_build_registered_ribbon(
    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonFitResult *fit)
{
    if (arena == NULL || ps == NULL || solution == NULL || set == NULL ||
        fit == NULL || !arf_options_valid(opts) ||
        solution->nvertices != ps->nv || solution->nfaces != ps->nf ||
        fit->registered_u == NULL || fit->row_chart_shift == NULL ||
        fit->nrows != set->nrows || fit->nrows == 0 ||
        set->nlayer_samples == 0 || !(opts->fit_u_spacing > 0.0))
        return -1;

    double registered_min = DBL_MAX, registered_max = -DBL_MAX;
    for (size_t i = 0; i < set->nlayer_samples; i++) {
        double u = fit->registered_u[i];
        if (!isfinite(u)) return -1;
        if (u < registered_min) registered_min = u;
        if (u > registered_max) registered_max = u;
    }
    if (!(registered_min < registered_max)) return -1;
    fit->du = opts->fit_u_spacing;
    fit->u0 = floor(registered_min / fit->du) * fit->du;
    int64_t ncolumns =
        (int64_t)ceil((registered_max - fit->u0) / fit->du) + 1;
    if (ncolumns < 2 || (uint64_t)ncolumns > (uint64_t)INT32_MAX ||
        fit->nrows > SIZE_MAX / (size_t)ncolumns)
        return -1;
    fit->ncolumns = (size_t)ncolumns;
    size_t nnode = fit->nrows * fit->ncolumns;
    fit->p = (double *)ARENA_ALLOC(arena, nnode * 2 * sizeof(*fit->p));
    fit->valid = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    fit->provenance = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    fit->support = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    fit->u_edge = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    fit->v_edge = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    fit->u_reject = (uint8_t *)ARENA_ALLOC(arena, nnode);
    memset(fit->u_reject, ATLAS_RIBBON_EDGE_NO_SUPPORT, nnode);
    fit->row_metric = (AtlasRibbonRowMetric *)ARENA_CALLOC(
        arena, fit->nrows, sizeof(*fit->row_metric));
    for (size_t i = 0; i < nnode * 2; i++) fit->p[i] = NAN;

    fit->fitted_nodes = fit->supported_nodes = fit->curve_segments = 0;
    fit->observed_nodes = fit->chart_interp_nodes = 0;
    fit->seam_interp_nodes = fit->spiral_fill_nodes = 0;
    fit->row_crossings = fit->u_gap_cuts = fit->metric_jump_cuts = 0;
    fit->crossing_cuts = 0;
    fit->direct_spans = fit->spiral_spans = fit->topology_cuts = 0;
    fit->winding_rows = fit->winding_phase_backsteps = 0;
    fit->winding_radius_order_tests = 0;
    fit->winding_radius_order_violations = 0;
    fit->parameter_samples = fit->chamfer_samples = fit->speed_samples = 0;
    fit->parameter_mean = fit->parameter_rms = NAN;
    fit->parameter_p50 = fit->parameter_p90 = NAN;
    fit->parameter_p95 = fit->parameter_p99 = fit->parameter_max = NAN;
    fit->chamfer_mean = fit->chamfer_rms = NAN;
    fit->chamfer_p50 = fit->chamfer_p90 = NAN;
    fit->chamfer_p95 = fit->chamfer_p99 = fit->directed_hausdorff = NAN;
    fit->speed_mean = fit->speed_rms_error = NAN;
    fit->speed_p95_error = fit->speed_max_error = NAN;

    Arena_Mark scratch_mark = Arena_save(arena);
    size_t max_row_samples = 0;
    for (size_t first = 0; first < set->nlayer_samples;) {
        size_t last = first + 1;
        while (last < set->nlayer_samples &&
               set->layer_sample[last].row == set->layer_sample[first].row)
            last++;
        if (last - first > max_row_samples) max_row_samples = last - first;
        first = last;
    }
    if (max_row_samples < 2) {
        Arena_restore(arena, scratch_mark);
        return -1;
    }
    ArfRegisterRef *ref = (ArfRegisterRef *)ARENA_ALLOC(
        arena, max_row_samples * sizeof(*ref));
    uint8_t *kind = (uint8_t *)ARENA_ALLOC(arena, max_row_samples);
    uint8_t *reject = (uint8_t *)ARENA_ALLOC(arena, max_row_samples);
    double *angle_delta = (double *)ARENA_ALLOC(
        arena, max_row_samples * sizeof(*angle_delta));

    size_t sample_at = 0;
    for (size_t row = 0; row < fit->nrows; row++) {
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row < (int32_t)row)
            sample_at++;
        size_t first = sample_at;
        while (sample_at < set->nlayer_samples &&
               set->layer_sample[sample_at].row == (int32_t)row)
            sample_at++;
        size_t count = sample_at - first;
        AtlasRibbonRowMetric *metric = &fit->row_metric[row];
        metric->row = (int32_t)row;
        if (count < 2) continue;
        for (size_t j = 0; j < count; j++) {
            ref[j].sample = first + j;
            ref[j].rank = set->layer_sample[first + j].rank;
            ref[j].u = fit->registered_u[first + j];
        }
        qsort(ref, count, sizeof(*ref), arf_compare_register_ref);
        size_t compact = 0;
        for (size_t j = 0; j < count; j++) {
            if (compact > 0 && ref[j].u <= ref[compact - 1].u + 1.0e-9) {
                const AtlasRibbonLayerSample *a =
                    &set->layer_sample[ref[compact - 1].sample];
                const AtlasRibbonLayerSample *b =
                    &set->layer_sample[ref[j].sample];
                if (!arf_layer_samples_local(a, b, opts)) {
                    metric->topology_cuts++;
                    fit->topology_cuts++;
                }
                continue;
            }
            ref[compact++] = ref[j];
        }
        count = compact;
        metric->accepted_targets = count;
        if (count < 2) continue;

        for (size_t j = 0; j + 1 < count; j++) {
            const AtlasRibbonLayerSample *a =
                &set->layer_sample[ref[j].sample];
            const AtlasRibbonLayerSample *b =
                &set->layer_sample[ref[j + 1].sample];
            double gap = ref[j + 1].u - ref[j].u;
            double limit = opts->max_bridge_stretch > 0.0
                ? opts->max_bridge_stretch * gap + opts->bridge_slack
                : DBL_MAX;
            kind[j] = 0;
            reject[j] = ATLAS_RIBBON_EDGE_TOPOLOGY;
            angle_delta[j] = 0.0;
            if (!(gap > 1.0e-9)) {
                metric->topology_cuts++;
                fit->topology_cuts++;
            } else if (opts->max_fill_u > 0.0 && gap > opts->max_fill_u) {
                reject[j] = ATLAS_RIBBON_EDGE_U_GAP;
                metric->u_gap_cuts++;
                fit->u_gap_cuts++;
            } else if (arf_layer_samples_local(a, b, opts)) {
                if (arf_distance2(a->p, b->p) <= limit) {
                    kind[j] = 1;
                    reject[j] = ATLAS_RIBBON_EDGE_SUBGRID;
                    metric->direct_spans++;
                    fit->direct_spans++;
                } else {
                    reject[j] = ATLAS_RIBBON_EDGE_METRIC;
                    metric->metric_jump_cuts++;
                    fit->metric_jump_cuts++;
                }
            } else if (a->rank > b->rank) {
                metric->topology_cuts++;
                fit->topology_cuts++;
            } else {
                double arc = 0.0;
                if (arf_registered_spiral_span(
                        a, b, fit->register_winding_direction,
                        &angle_delta[j], &arc) != 0) {
                    metric->topology_cuts++;
                    fit->topology_cuts++;
                } else if (arc > limit) {
                    reject[j] = ATLAS_RIBBON_EDGE_METRIC;
                    metric->metric_jump_cuts++;
                    fit->metric_jump_cuts++;
                } else {
                    kind[j] = 2;
                    reject[j] = ATLAS_RIBBON_EDGE_SUBGRID;
                    metric->spiral_spans++;
                    fit->spiral_spans++;
                }
            }
            arf_mark_u_reject(fit, row, ref[j].u, ref[j + 1].u,
                              (AtlasRibbonEdgeReject)reject[j]);
        }

        size_t chain = 0;
        while (chain + 1 < count) {
            while (chain + 1 < count && kind[chain] == 0) chain++;
            if (chain + 1 >= count) break;
            size_t chain_last = chain + 1;
            while (chain_last + 1 < count && kind[chain_last] != 0)
                chain_last++;
            double ua = ref[chain].u;
            double ub = ref[chain_last].u;
            int64_t column0 = (int64_t)ceil(
                (ua - fit->u0) / fit->du - 1.0e-10);
            int64_t column1 = (int64_t)floor(
                (ub - fit->u0) / fit->du + 1.0e-10);
            if (column0 < 0) column0 = 0;
            if (column1 >= (int64_t)fit->ncolumns)
                column1 = (int64_t)fit->ncolumns - 1;
            size_t link = chain;
            for (int64_t column = column0; column <= column1; column++) {
                double u = fit->u0 + (double)column * fit->du;
                while (link + 1 < chain_last &&
                       ref[link + 1].u < u - 1.0e-9)
                    link++;
                double den = ref[link + 1].u - ref[link].u;
                if (!(den > 1.0e-12)) continue;
                double t = (u - ref[link].u) / den;
                const AtlasRibbonLayerSample *a =
                    &set->layer_sample[ref[link].sample];
                const AtlasRibbonLayerSample *b =
                    &set->layer_sample[ref[link + 1].sample];
                double p[2];
                arf_registered_span_point(
                    a, b, fit->register_winding_direction, kind[link],
                    angle_delta[link], t, p);
                uint8_t provenance = ATLAS_RIBBON_NODE_SPIRAL_FILL;
                if (kind[link] == 1) {
                    if (fabs(t) <= 1.0e-9 || fabs(t - 1.0) <= 1.0e-9)
                        provenance = ATLAS_RIBBON_NODE_OBSERVED;
                    else if (a->chart == b->chart)
                        provenance = ATLAS_RIBBON_NODE_CHART_INTERP;
                    else
                        provenance = ATLAS_RIBBON_NODE_SEAM_INTERP;
                }
                size_t node = row * fit->ncolumns + (size_t)column;
                if (!fit->valid[node]) {
                    fit->p[node * 2] = p[0];
                    fit->p[node * 2 + 1] = p[1];
                    fit->valid[node] = 1;
                    fit->provenance[node] = provenance;
                } else if (fit->provenance[node] == ATLAS_RIBBON_NODE_INVALID ||
                           provenance < fit->provenance[node]) {
                    fit->provenance[node] = provenance;
                }
                if (provenance != ATLAS_RIBBON_NODE_SPIRAL_FILL)
                    fit->support[node] = 1;
            }
            if (column0 <= column1)
                for (int64_t column = column0; column < column1; column++) {
                    size_t node = row * fit->ncolumns + (size_t)column;
                    if (fit->valid[node] && fit->valid[node + 1]) {
                        fit->u_edge[node] = 1;
                        fit->u_reject[node] = ATLAS_RIBBON_EDGE_NONE;
                    }
                }
            chain = chain_last;
        }
    }
    if (sample_at != set->nlayer_samples) {
        Arena_restore(arena, scratch_mark);
        return -1;
    }

    for (size_t row = 0; row < fit->nrows; row++) {
        size_t cuts = arf_cut_row_crossings(arena, fit, row);
        fit->crossing_cuts += cuts;
        fit->topology_cuts += cuts;
        fit->row_metric[row].crossing_cuts += cuts;
        fit->row_metric[row].topology_cuts += cuts;
    }

    double vertical_limit = opts->max_bridge_stretch > 0.0
        ? opts->max_bridge_stretch * fit->dv + opts->bridge_slack
        : DBL_MAX;
    for (size_t row = 0; row + 1 < fit->nrows; row++)
        for (size_t column = 0; column < fit->ncolumns; column++) {
            size_t a = row * fit->ncolumns + column;
            size_t b = a + fit->ncolumns;
            if (!fit->valid[a] || !fit->valid[b]) continue;
            double dp0 = fit->p[b * 2] - fit->p[a * 2];
            double dp1 = fit->p[b * 2 + 1] - fit->p[a * 2 + 1];
            double distance = sqrt(dp0 * dp0 + dp1 * dp1 + fit->dv * fit->dv);
            if (distance <= vertical_limit) fit->v_edge[a] = 1;
        }

    double *speed_error = (double *)ARENA_ALLOC(
        arena, nnode * sizeof(*speed_error));
    size_t nspeed = 0;
    for (size_t row = 0; row < fit->nrows; row++) {
        AtlasRibbonRowMetric *metric = &fit->row_metric[row];
        double speed_sum = 0.0, speed_error2 = 0.0, speed_max = 0.0;
        int have_extent = 0;
        for (size_t column = 0; column < fit->ncolumns; column++) {
            size_t node = row * fit->ncolumns + column;
            if (fit->valid[node]) {
                metric->fitted_nodes++;
                fit->fitted_nodes++;
                switch (fit->provenance[node]) {
                case ATLAS_RIBBON_NODE_OBSERVED:
                    fit->observed_nodes++;
                    break;
                case ATLAS_RIBBON_NODE_CHART_INTERP:
                    fit->chart_interp_nodes++;
                    break;
                case ATLAS_RIBBON_NODE_SEAM_INTERP:
                    fit->seam_interp_nodes++;
                    break;
                case ATLAS_RIBBON_NODE_SPIRAL_FILL:
                    fit->spiral_fill_nodes++;
                    break;
                default:
                    break;
                }
                if (fit->support[node]) {
                    metric->supported_nodes++;
                    fit->supported_nodes++;
                }
                if (!have_extent) {
                    metric->u_min = fit->u0 + (double)column * fit->du;
                    have_extent = 1;
                }
                metric->u_max = fit->u0 + (double)column * fit->du;
            }
            if (column + 1 >= fit->ncolumns || !fit->u_edge[node]) continue;
            double dx = fit->p[(node + 1) * 2] - fit->p[node * 2];
            double dy = fit->p[(node + 1) * 2 + 1] - fit->p[node * 2 + 1];
            double speed = sqrt(dx * dx + dy * dy) / fit->du;
            double error = fabs(speed - 1.0);
            speed_error[nspeed++] = error;
            speed_sum += speed;
            speed_error2 += error * error;
            if (error > speed_max) speed_max = error;
            metric->curve_segments++;
            fit->curve_segments++;
        }
        if (metric->curve_segments > 0) {
            metric->speed_mean = speed_sum / (double)metric->curve_segments;
            metric->speed_rms_error =
                sqrt(speed_error2 / (double)metric->curve_segments);
            metric->speed_max_error = speed_max;
            metric->crossings = arf_row_crossings(arena, fit, row);
            fit->row_crossings += metric->crossings;
            arf_measure_row_winding(arena, fit, row, metric);
            if (metric->winding_span > 0.0) {
                fit->winding_rows++;
                fit->winding_phase_backsteps += (size_t)llround(
                    metric->phase_backstep_fraction *
                    (double)metric->curve_segments);
                fit->winding_radius_order_tests += metric->radius_order_tests;
                fit->winding_radius_order_violations +=
                    metric->radius_order_violations;
            }
        }
    }
    fit->speed_samples = nspeed;
    if (nspeed > 0) {
        double dummy50, dummy90, dummy99;
        arf_distribution(speed_error, nspeed,
                         &fit->speed_mean, &fit->speed_rms_error,
                         &dummy50, &dummy90, &fit->speed_p95_error,
                         &dummy99, &fit->speed_max_error);
        double speed_sum = 0.0;
        for (size_t row = 0; row < fit->nrows; row++)
            speed_sum += fit->row_metric[row].speed_mean *
                         (double)fit->row_metric[row].curve_segments;
        fit->speed_mean = speed_sum / (double)nspeed;
    }
    double *dense_chart_shift = NULL;
    if (arf_build_dense_chart_shift(
            arena, solution, set, fit, &dense_chart_shift) != 0) {
        Arena_restore(arena, scratch_mark);
        return -1;
    }
    arf_measure_vertex_distances(
        arena, ps, solution, set, fit, dense_chart_shift, solution->ncharts);
    int ok = fit->fitted_nodes > 0 && fit->curve_segments > 0;
    Arena_restore(arena, scratch_mark);
    return ok ? 0 : -1;
}

int AtlasRibbonFit_solve(
    Arena_T arena,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const AtlasRibbonObservationSet *set,
    const AtlasRibbonFitOptions *opts,
    AtlasRibbonFitResult *out)
{
    if (arena == NULL || ps == NULL || solution == NULL || set == NULL ||
        out == NULL || !arf_options_valid(opts))
        return -1;
    memset(out, 0, sizeof *out);
    out->parameter_mean = out->parameter_rms = NAN;
    out->parameter_p50 = out->parameter_p90 = NAN;
    out->parameter_p95 = out->parameter_p99 = out->parameter_max = NAN;
    out->chamfer_mean = out->chamfer_rms = NAN;
    out->chamfer_p50 = out->chamfer_p90 = NAN;
    out->chamfer_p95 = out->chamfer_p99 = NAN;
    out->directed_hausdorff = NAN;
    out->speed_mean = out->speed_rms_error = NAN;
    out->speed_p95_error = out->speed_max_error = NAN;

    if (opts->mode == ATLAS_RIBBON_OBSERVATIONS) return 0;
    if (opts->mode == ATLAS_RIBBON_REGISTER_FIELD_ONLY ||
        opts->mode == ATLAS_RIBBON_REGISTER_FIELD_SMOOTH_ONLY ||
        opts->mode == ATLAS_RIBBON_REGISTER_COLLISION_ONLY ||
        opts->mode == ATLAS_RIBBON_REGISTERED_RIBBON) {
        if (arf_register_u_field(arena, solution, set, opts, out) != 0)
            return -1;
        return 0;
    }
    if (arf_register_u(arena, solution, set, opts, out) != 0)
        return -1;
    if (opts->mode == ATLAS_RIBBON_REGISTER_ONLY) return 0;


    out->nrows = set->nrows;
    out->v0 = set->v0;
    out->dv = set->dv;
    out->du = opts->fit_u_spacing;
    out->u0 = floor(set->u_min / out->du) * out->du;
    int64_t ncolumns =
        (int64_t)ceil((set->u_max - out->u0) / out->du) + 1;
    if (ncolumns < 2 || (uint64_t)ncolumns > (uint64_t)INT32_MAX ||
        out->nrows > SIZE_MAX / (size_t)ncolumns)
        return -1;
    out->ncolumns = (size_t)ncolumns;
    size_t nnode = out->nrows * out->ncolumns;
    out->p = (double *)ARENA_ALLOC(arena, nnode * 2 * sizeof(double));
    out->valid = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    out->support = (uint8_t *)ARENA_CALLOC(arena, nnode, 1);
    out->row_metric = (AtlasRibbonRowMetric *)ARENA_CALLOC(
        arena, out->nrows, sizeof(AtlasRibbonRowMetric));
    for (size_t i = 0; i < nnode * 2; i++) out->p[i] = NAN;

    size_t target_at = 0;
    for (size_t row = 0; row < out->nrows; row++) {
        while (target_at < set->ntarget &&
               set->target[target_at].row < (int32_t)row)
            target_at++;
        size_t row_first = target_at;
        while (target_at < set->ntarget &&
               set->target[target_at].row == (int32_t)row)
            target_at++;
        size_t row_last = target_at;
        size_t accepted_count = 0;
        for (size_t i = row_first; i < row_last; i++)
            if (set->target[i].accepted) accepted_count++;
        AtlasRibbonRowMetric *metric = &out->row_metric[row];
        metric->row = (int32_t)row;
        metric->accepted_targets = accepted_count;
        if (accepted_count < 2) continue;
        Arena_Mark mark = Arena_save(arena);
        size_t *accepted = (size_t *)ARENA_ALLOC(
            arena, accepted_count * sizeof(size_t));
        size_t at = 0;
        for (size_t i = row_first; i < row_last; i++)
            if (set->target[i].accepted) accepted[at++] = i;
        size_t first = 0;
        while (first < accepted_count) {
            size_t last = first + 1;
            while (last < accepted_count) {
                double gap = set->target[accepted[last]].u -
                             set->target[accepted[last - 1]].u;
                int cut_u = opts->max_fill_u > 0.0 &&
                            gap > opts->max_fill_u;
                const AtlasRibbonTarget *ta =
                    &set->target[accepted[last - 1]];
                const AtlasRibbonTarget *tb =
                    &set->target[accepted[last]];
                double bridge_distance = arf_distance2(ta->p, tb->p);
                int cut_metric = opts->max_bridge_stretch > 0.0 &&
                    bridge_distance >
                        opts->max_bridge_stretch * gap + opts->bridge_slack;
                if (cut_u || cut_metric) {
                    if (cut_u) {
                        metric->u_gap_cuts++;
                        out->u_gap_cuts++;
                    } else {
                        metric->metric_jump_cuts++;
                        out->metric_jump_cuts++;
                    }
                    break;
                }
                last++;
            }
            if (last - first >= 2 &&
                arf_solve_segment(arena, set, opts, accepted,
                                  first, last, out, row) != 0) {
                Arena_restore(arena, mark);
                return -1;
            }
            first = last;
        }
        Arena_restore(arena, mark);
    }

    double *speed_error = (double *)ARENA_ALLOC(
        arena, nnode * sizeof(double));
    size_t nspeed = 0;
    for (size_t row = 0; row < out->nrows; row++) {
        AtlasRibbonRowMetric *metric = &out->row_metric[row];
        double speed_sum = 0.0, speed_error2 = 0.0, speed_max = 0.0;
        int have_extent = 0;
        for (size_t column = 0; column < out->ncolumns; column++) {
            size_t index = row * out->ncolumns + column;
            if (out->valid[index]) {
                metric->fitted_nodes++;
                if (out->support[index]) metric->supported_nodes++;
                if (!have_extent) {
                    metric->u_min = out->u0 + (double)column * out->du;
                    have_extent = 1;
                }
                metric->u_max = out->u0 + (double)column * out->du;
            }
            if (column + 1 >= out->ncolumns ||
                !out->valid[index] || !out->valid[index + 1])
                continue;
            double dx = out->p[(index + 1) * 2] - out->p[index * 2];
            double dy = out->p[(index + 1) * 2 + 1] -
                        out->p[index * 2 + 1];
            double speed = sqrt(dx * dx + dy * dy) / out->du;
            double error = fabs(speed - 1.0);
            speed_error[nspeed++] = error;
            speed_sum += speed;
            speed_error2 += error * error;
            if (error > speed_max) speed_max = error;
            metric->curve_segments++;
            out->curve_segments++;
        }
        if (metric->curve_segments > 0) {
            metric->speed_mean =
                speed_sum / (double)metric->curve_segments;
            metric->speed_rms_error =
                sqrt(speed_error2 / (double)metric->curve_segments);
            metric->speed_max_error = speed_max;
            metric->crossings = arf_row_crossings(arena, out, row);
            out->row_crossings += metric->crossings;
            arf_measure_row_winding(arena, out, row, metric);
            if (metric->winding_span > 0.0) {
                out->winding_rows++;
                out->winding_phase_backsteps +=
                    (size_t)llround(metric->phase_backstep_fraction *
                                    (double)metric->curve_segments);
                out->winding_radius_order_tests +=
                    metric->radius_order_tests;
                out->winding_radius_order_violations +=
                    metric->radius_order_violations;
            }
        }
    }
    out->speed_samples = nspeed;
    if (nspeed > 0) {
        double dummy50, dummy90, dummy99;
        arf_distribution(speed_error, nspeed,
                         &out->speed_mean, &out->speed_rms_error,
                         &dummy50, &dummy90, &out->speed_p95_error,
                         &dummy99, &out->speed_max_error);
        double speed_sum = 0.0;
        for (size_t row = 0; row < out->nrows; row++)
            speed_sum += out->row_metric[row].speed_mean *
                         (double)out->row_metric[row].curve_segments;
        out->speed_mean = speed_sum / (double)nspeed;
    }
    arf_measure_vertex_distances(arena, ps, solution, set, out, NULL, 0);
    return 0;
}

int AtlasRibbonFit_selftest(void)
{
    Arena_T arena = Arena_new();
    if (arena == NULL) return 1;
    enum { LAYERS = 3, VERTS_PER_LAYER = 4, FACES_PER_LAYER = 2 };
    float vertices[LAYERS * VERTS_PER_LAYER * 3];
    double u[LAYERS * VERTS_PER_LAYER];
    double v[LAYERS * VERTS_PER_LAYER];
    int32_t vertex_chart[LAYERS * VERTS_PER_LAYER];
    int32_t faces[LAYERS * FACES_PER_LAYER * 3];
    int32_t face_cube[LAYERS * FACES_PER_LAYER];
    uint8_t keep[LAYERS * FACES_PER_LAYER];
    memset(vertices, 0, sizeof vertices);
    for (int layer = 0; layer < LAYERS; layer++) {
        double y = layer < 2 ? (double)layer : 100.0;
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
                int index = layer * 4 + j * 2 + i;
                vertices[index * 3] = (float)j;
                vertices[index * 3 + 1] = (float)y;
                vertices[index * 3 + 2] = (float)(10 * i);
                u[index] = (double)(10 * i);
                v[index] = (double)j;
                vertex_chart[index] = layer;
            }
        int vo = layer * 4;
        int fo = layer * 6;
        faces[fo] = vo;
        faces[fo + 1] = vo + 1;
        faces[fo + 2] = vo + 3;
        faces[fo + 3] = vo;
        faces[fo + 4] = vo + 3;
        faces[fo + 5] = vo + 2;
        face_cube[layer * 2] = face_cube[layer * 2 + 1] = 0;
        keep[layer * 2] = keep[layer * 2 + 1] = 1;
    }
    char ids[1][48];
    memset(ids, 0, sizeof ids);
    size_t cube_voff[2] = {0, LAYERS * VERTS_PER_LAYER};
    PieceSet ps;
    memset(&ps, 0, sizeof ps);
    ps.verts = vertices;
    ps.faces = faces;
    ps.face_cube = face_cube;
    ps.ids = ids;
    ps.cube_voff = cube_voff;
    ps.nv = LAYERS * VERTS_PER_LAYER;
    ps.nf = LAYERS * FACES_PER_LAYER;
    ps.n_cubes = 1;
    AtlasSolution solution;
    memset(&solution, 0, sizeof solution);
    solution.nvertices = ps.nv;
    solution.nfaces = ps.nf;
    solution.ncharts = LAYERS;
    solution.u = u;
    solution.v = v;
    solution.face_keep = keep;
    solution.vertex_chart = vertex_chart;
    AtlasSolutionChart chart[LAYERS];
    memset(chart, 0, sizeof chart);
    for (int i = 0; i < LAYERS; i++) {
        chart[i].chart = i;
        chart[i].rank = (uint64_t)i;
    }
    solution.chart = chart;
    ScaffoldCalib cal;
    memset(&cal, 0, sizeof cal);
    cal.axis_dir[0] = 1.0f;
    cal.pitch = 9.5;
    AtlasRibbonFitOptions opts;
    AtlasRibbonFitOptions_default(&opts);
    opts.slice_spacing = 1.0;
    opts.observation_u_spacing = 2.0;
    opts.fit_u_spacing = 2.0;
    opts.local_xyz_tolerance = 8.0;
    AtlasRibbonObservationSet conflict;
    int failures = 0;
    if (AtlasRibbonFit_build_observations(
            arena, &ps, &solution, &cal, &opts, &conflict) != 0 ||
        conflict.conflict_targets == 0 ||
        conflict.accepted_targets != 0) {
        fprintf(stderr, "[atlas_ribbon_fit selftest] conflict fixture failed\n");
        failures++;
    }

    keep[4] = keep[5] = 0;
    Arena_Mark mark = Arena_save(arena);
    AtlasRibbonObservationSet local;
    if (AtlasRibbonFit_build_observations(
            arena, &ps, &solution, &cal, &opts, &local) != 0 ||
        local.accepted_targets == 0 ||
        local.duplicate_targets == 0 ||
        local.conflict_targets != 0) {
        fprintf(stderr, "[atlas_ribbon_fit selftest] duplicate fixture failed\n");
        failures++;
    } else {
        AtlasRibbonFitResult fit;
        opts.mode = ATLAS_RIBBON_TANGENT;
        opts.lambda_smooth = 1.0;
        opts.lambda_tangent = 8.0;
        if (AtlasRibbonFit_solve(
                arena, &ps, &solution, &local, &opts, &fit) != 0 ||
            fit.fitted_nodes == 0 || fit.row_crossings != 0 ||
            !(fit.speed_rms_error < 0.1) ||
            !(fit.chamfer_p95 < 2.0)) {
            fprintf(stderr,
                    "[atlas_ribbon_fit selftest] row fit failed "
                    "(nodes=%zu crossings=%zu speed=%.6g chamfer=%.6g)\n",
                    fit.fitted_nodes, fit.row_crossings,
                    fit.speed_rms_error, fit.chamfer_p95);
            failures++;
        }
    }
    Arena_restore(arena, mark);
    Arena_dispose(&arena);
    if (failures == 0)
        fprintf(stderr, "[atlas_ribbon_fit selftest] PASS\n");
    return failures;
}
