#include "atlas_collision_register.h"

#include "atlas_overlap_audit.h"
#include "monotone_qp.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t *faces;
    size_t nfaces;
    size_t nvertices;
    double *base_u;
    double *v;
    float *dummy;
    int32_t *component;
    int32_t *source_face;
    float *xyz;
} AcrSoup;

typedef struct {
    size_t at;
    int32_t lo;
    int32_t hi;
    int32_t face_lo;
    int32_t face_hi;
    double lower;
    double priority;
} AcrCollisionCandidate;

static int acr_candidate_compare(const void *left, const void *right)
{
    const AcrCollisionCandidate *a = (const AcrCollisionCandidate *)left;
    const AcrCollisionCandidate *b = (const AcrCollisionCandidate *)right;
    if (a->priority < b->priority) return -1;
    if (a->priority > b->priority) return 1;
    return a->at < b->at ? -1 : a->at > b->at;
}

static double acr_kind_scale(AtlasRegisterEdgeKind kind)
{
    return kind == ATLAS_REGISTER_RADIAL ? 0.25 : 1.0;
}

static int acr_edge_usable(const AtlasCollisionRegisterProblem *p, size_t i)
{
    const AtlasRegisterEdge *e = &p->edges[i];
    return e->component0 >= 0 && e->component1 >= 0 &&
           (size_t)e->component0 < p->ncharts &&
           (size_t)e->component1 < p->ncharts &&
           e->component0 != e->component1 && isfinite(e->target) &&
           isfinite(e->weight) && e->weight > 0.0 &&
           isfinite(e->robust_scale) &&
           (e->kind == ATLAS_REGISTER_WELD ||
            e->kind == ATLAS_REGISTER_RADIAL ||
            e->kind == ATLAS_REGISTER_LOCAL);
}


/* Freeze the Huber IRLS weight at the unconstrained metric solution.  The
 * collision QP is then the local quadratic model of the same robust objective,
 * rather than allowing one bad metric edge to dominate once a bound activates. */
static double acr_edge_frozen_weight(const AtlasCollisionRegisterProblem *p,
                                     const AtlasRegisterEdge *e)
{
    double residual = -e->target;
    if (p->desired_shift != NULL)
        residual += p->desired_shift[e->component1] -
                    p->desired_shift[e->component0];
    double robust = 1.0;
    if (e->robust_scale > 0.0 && fabs(residual) > e->robust_scale)
        robust = e->robust_scale / fabs(residual);
    return e->weight * acr_kind_scale(e->kind) * robust;
}
static int acr_build_soup(Arena_T arena,
                          const AtlasCollisionRegisterProblem *p,
                          AcrSoup *soup)
{
    memset(soup, 0, sizeof(*soup));
    size_t kept = 0;
    for (size_t f = 0; f < p->nfaces; f++) {
        if (p->face_keep != NULL && !p->face_keep[f]) continue;
        int32_t a = p->faces[f * 3];
        int32_t b = p->faces[f * 3 + 1];
        int32_t c = p->faces[f * 3 + 2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= p->nvertices || (size_t)b >= p->nvertices ||
            (size_t)c >= p->nvertices)
            return -1;
        int32_t chart = p->vertex_chart[a];
        if (chart < 0 || (size_t)chart >= p->ncharts ||
            p->vertex_chart[b] != chart || p->vertex_chart[c] != chart)
            return -1;
        kept++;
    }
    if (kept == 0 || kept > SIZE_MAX / 3) return -1;
    size_t corners = kept * 3;
    soup->faces = (int32_t *)ARENA_ALLOC(
        arena, corners * sizeof(*soup->faces));
    soup->base_u = (double *)ARENA_ALLOC(
        arena, corners * sizeof(*soup->base_u));
    soup->v = (double *)ARENA_ALLOC(arena, corners * sizeof(*soup->v));
    soup->dummy = (float *)ARENA_CALLOC(
        arena, corners, sizeof(*soup->dummy));
    soup->component = (int32_t *)ARENA_ALLOC(
        arena, corners * sizeof(*soup->component));
    soup->source_face = (int32_t *)ARENA_ALLOC(
        arena, kept * sizeof(*soup->source_face));
    if (p->xyz != NULL)
        soup->xyz = (float *)ARENA_ALLOC(
            arena, corners * 3 * sizeof(*soup->xyz));
    size_t out = 0;
    size_t face_out = 0;
    for (size_t f = 0; f < p->nfaces; f++) {
        if (p->face_keep != NULL && !p->face_keep[f]) continue;
        int32_t chart = p->vertex_chart[p->faces[f * 3]];
        soup->source_face[face_out++] = (int32_t)f;
        for (int k = 0; k < 3; k++) {
            int32_t source = p->faces[f * 3 + (size_t)k];
            soup->faces[out] = (int32_t)out;
            soup->base_u[out] = p->base_u[source];
            soup->v[out] = p->v[source];
            soup->component[out] = chart;
            if (soup->xyz != NULL) {
                soup->xyz[out * 3] = p->xyz[(size_t)source * 3];
                soup->xyz[out * 3 + 1] = p->xyz[(size_t)source * 3 + 1];
                soup->xyz[out * 3 + 2] = p->xyz[(size_t)source * 3 + 2];
            }
            out++;
        }
    }
    if (out != corners || face_out != kept) return -1;
    soup->nfaces = kept;
    soup->nvertices = corners;
    return 0;
}

static int acr_triangle_u_slice(const AcrSoup *soup, size_t face,
                                double axial, double *out_min,
                                double *out_max)
{
    double minimum = DBL_MAX, maximum = -DBL_MAX;
    int intersections = 0;
    double scale = fmax(1.0, fabs(axial));
    for (int edge = 0; edge < 3; edge++) {
        int32_t va = soup->faces[face * 3 + (size_t)edge];
        int32_t vb = soup->faces[face * 3 + (size_t)((edge + 1) % 3)];
        double aa = soup->v[va], ab = soup->v[vb];
        double ua = soup->base_u[va], ub = soup->base_u[vb];
        scale = fmax(scale, fmax(fabs(aa), fabs(ab)));
        double tolerance = 1.0e-12 * scale;
        double dv = ab - aa;
        if (fabs(dv) <= tolerance) {
            if (fabs(axial - aa) <= tolerance) {
                if (ua < minimum) minimum = ua;
                if (ua > maximum) maximum = ua;
                if (ub < minimum) minimum = ub;
                if (ub > maximum) maximum = ub;
                intersections += 2;
            }
            continue;
        }
        double lo = aa < ab ? aa : ab;
        double hi = aa < ab ? ab : aa;
        if (axial < lo - tolerance || axial > hi + tolerance) continue;
        double t = (axial - aa) / dv;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        double u = ua + t * (ub - ua);
        if (u < minimum) minimum = u;
        if (u > maximum) maximum = u;
        intersections++;
    }
    if (intersections == 0 || !isfinite(minimum) || !isfinite(maximum) ||
        minimum > maximum)
        return -1;
    *out_min = minimum;
    *out_max = maximum;
    return 0;
}

/* Complete open interval of relative translations delta = shift[face1] -
 * shift[face0] for which two triangles have positive-area intersection. */
static int acr_face_pair_interval(const AcrSoup *soup, size_t face0,
                                  size_t face1, double *out_left,
                                  double *out_right)
{
    double min_v0 = DBL_MAX, max_v0 = -DBL_MAX;
    double min_v1 = DBL_MAX, max_v1 = -DBL_MAX;
    double candidate[8];
    size_t ncandidate = 0;
    for (int k = 0; k < 3; k++) {
        double value0 = soup->v[soup->faces[face0 * 3 + (size_t)k]];
        double value1 = soup->v[soup->faces[face1 * 3 + (size_t)k]];
        if (value0 < min_v0) min_v0 = value0;
        if (value0 > max_v0) max_v0 = value0;
        if (value1 < min_v1) min_v1 = value1;
        if (value1 > max_v1) max_v1 = value1;
    }
    double common_min = fmax(min_v0, min_v1);
    double common_max = fmin(max_v0, max_v1);
    double vscale = fmax(1.0, fmax(fabs(common_min), fabs(common_max)));
    if (common_min > common_max + 1.0e-12 * vscale) return -1;
    candidate[ncandidate++] = common_min;
    candidate[ncandidate++] = common_max;
    for (int k = 0; k < 3; k++) {
        double value0 = soup->v[soup->faces[face0 * 3 + (size_t)k]];
        double value1 = soup->v[soup->faces[face1 * 3 + (size_t)k]];
        if (value0 > common_min && value0 < common_max)
            candidate[ncandidate++] = value0;
        if (value1 > common_min && value1 < common_max)
            candidate[ncandidate++] = value1;
    }
    double left = DBL_MAX, right = -DBL_MAX;
    for (size_t i = 0; i < ncandidate; i++) {
        double min0, max0, min1, max1;
        if (acr_triangle_u_slice(soup, face0, candidate[i],
                                 &min0, &max0) != 0 ||
            acr_triangle_u_slice(soup, face1, candidate[i],
                                 &min1, &max1) != 0)
            return -1;
        double local_left = min0 - max1;
        double local_right = max0 - min1;
        if (local_left < left) left = local_left;
        if (local_right > right) right = local_right;
    }
    if (!isfinite(left) || !isfinite(right) || !(left < right)) return -1;
    *out_left = left;
    *out_right = right;
    return 0;
}

static double acr_face_min_xyz(const AcrSoup *soup,
                               size_t face0, size_t face1);

static int acr_pair_allowed(const AtlasCollisionRegisterProblem *p,
                            const AcrSoup *soup, size_t face0, size_t face1,
                            int *out_reason)
{
    int32_t chart0 = soup->component[soup->faces[face0 * 3]];
    int32_t chart1 = soup->component[soup->faces[face1 * 3]];
    double min_xyz = acr_face_min_xyz(soup, face0, face1);
    int direct = p->allowed_chart_pair != NULL &&
        (p->allowed_chart_pair[(size_t)chart0 * p->ncharts +
                               (size_t)chart1] ||
         p->allowed_chart_pair[(size_t)chart1 * p->ncharts +
                               (size_t)chart0]);
    int reason = 0;
    if (chart0 == chart1) {
        reason = 3;
    } else if (direct && (min_xyz < 0.0 ||
               p->same_sheet_xyz_tolerance <= 0.0 ||
               min_xyz <= p->same_sheet_xyz_tolerance)) {
        reason = 1;
    } else if (min_xyz >= 0.0 && p->chart_radius != NULL &&
               p->same_sheet_radius_tolerance > 0.0 &&
               p->same_sheet_xyz_tolerance > 0.0 &&
               min_xyz <= p->same_sheet_xyz_tolerance &&
               isfinite(p->chart_radius[chart0]) &&
               isfinite(p->chart_radius[chart1]) &&
               fabs(p->chart_radius[chart0] - p->chart_radius[chart1]) <=
                   p->same_sheet_radius_tolerance) {
        reason = 2;
    }
    if (out_reason != NULL) *out_reason = reason;
    return reason != 0;
}

static double acr_face_min_xyz(const AcrSoup *soup,
                               size_t face0, size_t face1)
{
    if (soup->xyz == NULL) return -1.0;
    double best2 = DBL_MAX;
    for (int a = 0; a < 3; a++) {
        int32_t va = soup->faces[face0 * 3 + (size_t)a];
        for (int b = 0; b < 3; b++) {
            int32_t vb = soup->faces[face1 * 3 + (size_t)b];
            double dx = (double)soup->xyz[(size_t)va * 3] -
                        (double)soup->xyz[(size_t)vb * 3];
            double dy = (double)soup->xyz[(size_t)va * 3 + 1] -
                        (double)soup->xyz[(size_t)vb * 3 + 1];
            double dz = (double)soup->xyz[(size_t)va * 3 + 2] -
                        (double)soup->xyz[(size_t)vb * 3 + 2];
            double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < best2) best2 = d2;
        }
    }
    return best2 < DBL_MAX ? sqrt(best2) : -1.0;
}


/* Existing difference constraints are feasible.  Tightening lo->hi by
 * `lower` creates a positive cycle exactly when the longest existing path
 * hi->lo closes above zero. */
static int acr_can_tighten(size_t ncharts,
                           const int32_t *edge_lo, const int32_t *edge_hi,
                           const double *edge_lower, size_t nedges,
                           int32_t lo, int32_t hi, double lower,
                           double *distance)
{
    for (size_t c = 0; c < ncharts; c++) distance[c] = -DBL_MAX;
    distance[hi] = 0.0;
    for (size_t pass = 0; pass + 1 < ncharts; pass++) {
        int changed = 0;
        for (size_t e = 0; e < nedges; e++) {
            int32_t a = edge_lo[e], b = edge_hi[e];
            if (!(distance[a] > -DBL_MAX)) continue;
            double candidate = distance[a] + edge_lower[e];
            double tolerance = 1.0e-10 *
                fmax(1.0, fmax(fabs(candidate), fabs(distance[b])));
            if (candidate > distance[b] + tolerance) {
                distance[b] = candidate;
                changed = 1;
            }
        }
        if (!changed) break;
    }
    if (!(distance[lo] > -DBL_MAX)) return 1;
    double cycle = distance[lo] + lower;
    double tolerance = 1.0e-9 *
        fmax(1.0, fmax(fabs(distance[lo]), fabs(lower)));
    return cycle <= tolerance;
}

static int acr_make_feasible(size_t ncharts,
                             const MonotoneQpBound *bound, size_t nbounds,
                             double *shift)
{
    for (size_t pass = 0; pass < ncharts; pass++) {
        int changed = 0;
        for (size_t e = 0; e < nbounds; e++) {
            int32_t lo = bound[e].lo, hi = bound[e].hi;
            double candidate = shift[lo] + bound[e].lower;
            double tolerance = 1.0e-10 *
                fmax(1.0, fmax(fabs(candidate), fabs(shift[hi])));
            if (candidate > shift[hi] + tolerance) {
                shift[hi] = candidate;
                changed = 1;
            }
        }
        if (!changed) return 0;
    }
    return -1;
}
static int acr_collect_bounds(Arena_T work,
                              const AtlasCollisionRegisterProblem *p,
                              const AcrSoup *soup, const size_t *rank,
                              const double *shift, double margin,
                              double *lower, int32_t *bound_face_lo,
                              int32_t *bound_face_hi,
                              uint32_t *bound_updates,
                              uint8_t *bound_collision,
                              uint8_t *blocked_pair,
                              size_t *out_updates, size_t *out_rejected,
                              size_t *out_exact, size_t *out_cross,
                              size_t *out_allowed, size_t *out_hard,
                              size_t *out_topology_limited,
                              AtlasOverlapAudit *out_audit)
{
    size_t ncharts = p->ncharts;
    size_t matrix_size = ncharts * ncharts;
    double *u = (double *)ARENA_ALLOC(
        work, soup->nvertices * sizeof(*u));
    for (size_t i = 0; i < soup->nvertices; i++)
        u[i] = soup->base_u[i] + shift[soup->component[i]];
    AtlasOverlapAudit audit;
    if (AtlasOverlapAudit_build(
            work, soup->faces, soup->nfaces, soup->nvertices, u, soup->v,
            soup->dummy, soup->dummy, soup->component, ncharts, &audit) != 0 ||
        !audit.broad_phase_complete || audit.indexed_faces != soup->nfaces)
        return -1;

    int32_t *candidate_index = (int32_t *)ARENA_ALLOC(
        work, matrix_size * sizeof(*candidate_index));
    for (size_t i = 0; i < matrix_size; i++) candidate_index[i] = -1;
    AcrCollisionCandidate *candidate = (AcrCollisionCandidate *)ARENA_ALLOC(
        work, (audit.npairs ? audit.npairs : 1) * sizeof(*candidate));
    size_t ncandidate = 0, allowed = 0, hard = 0;
    for (size_t i = 0; i < audit.npairs; i++) {
        const AtlasOverlapPair *pair = &audit.pairs[i];
        int32_t c0 = pair->component0, c1 = pair->component1;
        if (c0 == c1) continue;
        if (pair->face0 < 0 || pair->face1 < 0 ||
            (size_t)pair->face0 >= soup->nfaces ||
            (size_t)pair->face1 >= soup->nfaces)
            return -1;
        size_t f0 = (size_t)pair->face0, f1 = (size_t)pair->face1;
        int reason = 0;
        if (acr_pair_allowed(p, soup, f0, f1, &reason)) {
            allowed++;
            continue;
        }
        hard++;
        int32_t lo = c0, hi = c1;
        size_t flo = f0, fhi = f1;
        if (rank[c0] > rank[c1]) {
            lo = c1; hi = c0;
            flo = f1; fhi = f0;
        }
        double left, right;
        if (acr_face_pair_interval(soup, flo, fhi, &left, &right) != 0)
            return -1;
        double bound = right + margin;
        size_t at = (size_t)lo * ncharts + (size_t)hi;
        int32_t ci = candidate_index[at];
        if (ci < 0) {
            if (ncandidate > (size_t)INT32_MAX) return -1;
            ci = (int32_t)ncandidate++;
            candidate_index[at] = ci;
            candidate[ci].at = at;
            candidate[ci].lo = lo;
            candidate[ci].hi = hi;
            candidate[ci].face_lo = soup->source_face[flo];
            candidate[ci].face_hi = soup->source_face[fhi];
            candidate[ci].lower = bound;
        } else if (bound > candidate[ci].lower) {
            candidate[ci].face_lo = soup->source_face[flo];
            candidate[ci].face_hi = soup->source_face[fhi];
            candidate[ci].lower = bound;
        }
        (void)left;
    }

    for (size_t i = 0; i < ncandidate; i++) {
        int32_t lo = candidate[i].lo, hi = candidate[i].hi;
        double desired_lo = p->desired_shift != NULL
                          ? p->desired_shift[lo] : 0.0;
        double desired_hi = p->desired_shift != NULL
                          ? p->desired_shift[hi] : 0.0;
        candidate[i].priority = candidate[i].lower -
                                (desired_hi - desired_lo);
    }
    qsort(candidate, ncandidate, sizeof(*candidate), acr_candidate_compare);

    int32_t *edge_lo = (int32_t *)ARENA_ALLOC(
        work, matrix_size * sizeof(*edge_lo));
    int32_t *edge_hi = (int32_t *)ARENA_ALLOC(
        work, matrix_size * sizeof(*edge_hi));
    double *edge_lower = (double *)ARENA_ALLOC(
        work, matrix_size * sizeof(*edge_lower));
    size_t *edge_index = (size_t *)ARENA_ALLOC(
        work, matrix_size * sizeof(*edge_index));
    for (size_t i = 0; i < matrix_size; i++) edge_index[i] = SIZE_MAX;
    size_t nconstraint = 0;
    for (size_t lo = 0; lo < ncharts; lo++)
    for (size_t hi = 0; hi < ncharts; hi++) {
        size_t at = lo * ncharts + hi;
        if (!(lower[at] > -DBL_MAX)) continue;
        edge_index[at] = nconstraint;
        edge_lo[nconstraint] = (int32_t)lo;
        edge_hi[nconstraint] = (int32_t)hi;
        edge_lower[nconstraint] = lower[at];
        nconstraint++;
    }
    double *distance = (double *)ARENA_ALLOC(
        work, ncharts * sizeof(*distance));

    size_t updates = 0, rejected = 0;
    for (size_t i = 0; i < ncandidate; i++) {
        AcrCollisionCandidate *q = &candidate[i];
        double tolerance = 1.0e-10 * fmax(1.0, fabs(q->lower));
        if (q->lower <= lower[q->at] + tolerance) {
            bound_collision[q->at] = 1;
            blocked_pair[q->at] = 0;
            if (bound_face_lo[q->at] < 0) {
                bound_face_lo[q->at] = q->face_lo;
                bound_face_hi[q->at] = q->face_hi;
            }
            continue;
        }
        int compatible = p->topology_deviation_limit <= 0.0 ||
            acr_can_tighten(ncharts, edge_lo, edge_hi, edge_lower,
                            nconstraint, q->lo, q->hi, q->lower,
                            distance);
        if (!compatible) {
            blocked_pair[q->at] = 1;
            rejected++;
            continue;
        }
        lower[q->at] = q->lower;
        bound_face_lo[q->at] = q->face_lo;
        bound_face_hi[q->at] = q->face_hi;
        bound_collision[q->at] = 1;
        blocked_pair[q->at] = 0;
        if (bound_updates[q->at] < UINT32_MAX) bound_updates[q->at]++;
        size_t ei = edge_index[q->at];
        if (ei == SIZE_MAX) {
            ei = nconstraint++;
            edge_index[q->at] = ei;
            edge_lo[ei] = q->lo;
            edge_hi[ei] = q->hi;
        }
        edge_lower[ei] = q->lower;
        updates++;
    }

    size_t topology_limited = 0;
    for (size_t i = 0; i < audit.npairs; i++) {
        const AtlasOverlapPair *pair = &audit.pairs[i];
        int32_t c0 = pair->component0, c1 = pair->component1;
        if (c0 == c1) continue;
        size_t f0 = (size_t)pair->face0, f1 = (size_t)pair->face1;
        if (acr_pair_allowed(p, soup, f0, f1, NULL)) continue;
        int32_t lo = c0, hi = c1;
        if (rank[lo] > rank[hi]) { int32_t t = lo; lo = hi; hi = t; }
        if (blocked_pair[(size_t)lo * ncharts + (size_t)hi])
            topology_limited++;
    }

    *out_updates = updates;
    *out_rejected = rejected;
    *out_exact = audit.exact_face_pairs;
    *out_cross = audit.cross_component_pairs;
    *out_allowed = allowed;
    *out_hard = hard;
    *out_topology_limited = topology_limited;
    if (out_audit != NULL) *out_audit = audit;
    return 0;
}

static size_t acr_count_bounds(const double *lower, size_t ncharts)
{
    size_t count = 0;
    for (size_t a = 0; a < ncharts; a++)
    for (size_t b = 0; b < ncharts; b++)
        if (lower[a * ncharts + b] > -DBL_MAX) count++;
    return count;
}

static int acr_store_reports(Arena_T arena,
                             const AtlasCollisionRegisterProblem *p,
                             const AcrSoup *soup,
                             const double *lower,
                             const int32_t *bound_face_lo,
                             const int32_t *bound_face_hi,
                             const uint32_t *bound_updates,
                             const uint8_t *bound_topology,
                             const uint8_t *bound_collision,
                             const uint8_t *initial_bound,
                             const uint8_t *blocked_pair,
                             const size_t *rank,
                             const double *shift,
                             const AtlasOverlapAudit *audit,
                             AtlasCollisionRegisterStats *stats)
{
    size_t nbounds = acr_count_bounds(lower, p->ncharts);
    AtlasCollisionRegisterBound *bounds =
        (AtlasCollisionRegisterBound *)ARENA_ALLOC(
            arena, (nbounds ? nbounds : 1) * sizeof(*bounds));
    size_t bi = 0;
    for (size_t lo = 0; lo < p->ncharts; lo++)
    for (size_t hi = 0; hi < p->ncharts; hi++) {
        size_t at = lo * p->ncharts + hi;
        if (!(lower[at] > -DBL_MAX)) continue;
        if (bi >= nbounds || (!bound_topology[at] && !bound_collision[at]) ||
            (bound_collision[at] &&
             (bound_face_lo[at] < 0 || bound_face_hi[at] < 0)))
            return -1;
        double delta = shift[hi] - shift[lo];
        double slack = delta - lower[at];
        double tolerance = 1.0e-7 * fmax(1.0, fabs(lower[at]));
        bounds[bi].chart_lo = (int32_t)lo;
        bounds[bi].chart_hi = (int32_t)hi;
        bounds[bi].face_lo = bound_face_lo[at];
        bounds[bi].face_hi = bound_face_hi[at];
        bounds[bi].lower = lower[at];
        bounds[bi].final_delta = delta;
        bounds[bi].slack = slack;
        bounds[bi].updates = bound_updates[at];
        bounds[bi].topology = bound_topology[at] != 0;
        bounds[bi].collision = bound_collision[at] != 0;
        bounds[bi].initial = initial_bound[at] != 0;
        bounds[bi].active = slack <= tolerance;
        bi++;
    }
    if (bi != nbounds) return -1;

    AtlasCollisionRegisterResidual *residual =
        (AtlasCollisionRegisterResidual *)ARENA_ALLOC(
            arena, (audit->npairs ? audit->npairs : 1) * sizeof(*residual));
    for (size_t i = 0; i < audit->npairs; i++) {
        const AtlasOverlapPair *pair = &audit->pairs[i];
        if (pair->face0 < 0 || pair->face1 < 0 ||
            (size_t)pair->face0 >= soup->nfaces ||
            (size_t)pair->face1 >= soup->nfaces)
            return -1;
        residual[i].chart0 = pair->component0;
        residual[i].chart1 = pair->component1;
        residual[i].face0 = soup->source_face[pair->face0];
        residual[i].face1 = soup->source_face[pair->face1];
        residual[i].min_xyz = acr_face_min_xyz(
            soup, (size_t)pair->face0, (size_t)pair->face1);
        residual[i].radius_delta =
            p->chart_radius != NULL &&
            isfinite(p->chart_radius[pair->component0]) &&
            isfinite(p->chart_radius[pair->component1])
            ? fabs(p->chart_radius[pair->component0] -
                   p->chart_radius[pair->component1])
            : -1.0;
        int reason = 0;
        residual[i].allowed = acr_pair_allowed(
            p, soup, (size_t)pair->face0, (size_t)pair->face1, &reason);
        if (!residual[i].allowed && pair->component0 != pair->component1) {
            int32_t lo = pair->component0, hi = pair->component1;
            if (rank[lo] > rank[hi]) { int32_t t = lo; lo = hi; hi = t; }
            if (blocked_pair[(size_t)lo * p->ncharts + (size_t)hi])
                reason = 4;
        }
        residual[i].reason = reason;
    }

    stats->bound = bounds;
    stats->nbound = nbounds;
    stats->residual = residual;
    stats->nresidual = audit->npairs;
    return 0;
}

int AtlasCollisionRegister_solve(
    Arena_T arena, const AtlasCollisionRegisterProblem *p,
    double *out_shift, AtlasCollisionRegisterStats *stats)
{
    if (arena == NULL || p == NULL || out_shift == NULL || stats == NULL ||
        p->faces == NULL || p->nfaces == 0 ||
        p->nfaces > (size_t)INT32_MAX || p->nvertices == 0 ||
        p->base_u == NULL || p->v == NULL || p->vertex_chart == NULL ||
        p->ncharts == 0 || p->ncharts > (size_t)INT32_MAX ||
        p->chart_rank == NULL ||
        (p->nedges > 0 && p->edges == NULL) ||
        !isfinite(p->topology_deviation_limit) ||
        p->topology_deviation_limit < 0.0 ||
        !isfinite(p->collision_margin) || p->collision_margin < 0.0 ||
        !isfinite(p->same_sheet_radius_tolerance) ||
        p->same_sheet_radius_tolerance < 0.0 ||
        !isfinite(p->same_sheet_xyz_tolerance) ||
        p->same_sheet_xyz_tolerance < 0.0)
        return -1;
    memset(stats, 0, sizeof(*stats));
    memset(out_shift, 0, p->ncharts * sizeof(*out_shift));
    if (p->ncharts > SIZE_MAX / p->ncharts ||
        p->ncharts * p->ncharts > SIZE_MAX / sizeof(double))
        return -1;

    uint8_t *rank_seen = (uint8_t *)ARENA_CALLOC(
        arena, p->ncharts, sizeof(*rank_seen));
    for (size_t c = 0; c < p->ncharts; c++) {
        size_t r = p->chart_rank[c];
        if (r >= p->ncharts || rank_seen[r]) return -1;
        rank_seen[r] = 1;
    }

    AcrSoup soup;
    if (acr_build_soup(arena, p, &soup) != 0) return -1;
    stats->kept_faces = soup.nfaces;

    size_t nmetric = 0;
    for (size_t i = 0; i < p->nedges; i++) {
        if (acr_edge_usable(p, i)) nmetric++;
        else stats->rejected_edges++;
    }
    stats->metric_edges = nmetric;
    if (nmetric > SIZE_MAX / 2) return -1;
    if (p->ncharts > SIZE_MAX - nmetric ||
        p->ncharts > SIZE_MAX - 2 * nmetric)
        return -1;
    size_t nrows = p->ncharts + nmetric;
    size_t ncoeff = p->ncharts + 2 * nmetric;
    MonotoneQpRow *row = (MonotoneQpRow *)ARENA_ALLOC(
        arena, nrows * sizeof(*row));
    MonotoneQpCoeff *coeff = (MonotoneQpCoeff *)ARENA_ALLOC(
        arena, ncoeff * sizeof(*coeff));
    size_t cursor = 0;
    for (size_t c = 0; c < p->ncharts; c++) {
        double target = p->desired_shift != NULL ? p->desired_shift[c] : 0.0;
        double weight = p->desired_weight != NULL ? p->desired_weight[c] : 1.0;
        if (!isfinite(target) || !isfinite(weight) || !(weight > 0.0))
            return -1;
        row[c].first = cursor;
        row[c].count = 1;
        row[c].target = target;
        row[c].weight = weight;
        row[c].kind = 0;
        row[c].owner = (int32_t)c;
        coeff[cursor].var = (int32_t)c;
        coeff[cursor].value = 1.0;
        cursor++;
    }
    size_t ri = p->ncharts;
    for (size_t i = 0; i < p->nedges; i++) {
        if (!acr_edge_usable(p, i)) continue;
        const AtlasRegisterEdge *e = &p->edges[i];
        row[ri].first = cursor;
        row[ri].count = 2;
        row[ri].target = e->target;
        row[ri].weight = acr_edge_frozen_weight(p, e);
        row[ri].kind = 1 + (int32_t)e->kind;
        row[ri].owner = i <= (size_t)INT32_MAX ? (int32_t)i : INT32_MAX;
        coeff[cursor].var = e->component0;
        coeff[cursor].value = -1.0;
        cursor++;
        coeff[cursor].var = e->component1;
        coeff[cursor].value = 1.0;
        cursor++;
        ri++;
    }
    if (ri != nrows || cursor != ncoeff) return -1;

    size_t matrix_size = p->ncharts * p->ncharts;
    double *lower = (double *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*lower));
    int32_t *bound_face_lo = (int32_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bound_face_lo));
    int32_t *bound_face_hi = (int32_t *)ARENA_ALLOC(
        arena, matrix_size * sizeof(*bound_face_hi));
    uint32_t *bound_updates = (uint32_t *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*bound_updates));
    uint8_t *bound_topology = (uint8_t *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*bound_topology));
    uint8_t *bound_collision = (uint8_t *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*bound_collision));
    uint8_t *blocked_pair = (uint8_t *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*blocked_pair));
    uint8_t *initial_bound = (uint8_t *)ARENA_CALLOC(
        arena, matrix_size, sizeof(*initial_bound));
    for (size_t i = 0; i < matrix_size; i++) {
        lower[i] = -DBL_MAX;
        bound_face_lo[i] = bound_face_hi[i] = -1;
    }

    if (p->topology_deviation_limit > 0.0) {
        for (size_t i = 0; i < p->nedges; i++) {
            if (!acr_edge_usable(p, i)) continue;
            const AtlasRegisterEdge *e = &p->edges[i];
            if (e->kind == ATLAS_REGISTER_RADIAL) continue;
            int32_t a = e->component0, b = e->component1;
            double da = p->desired_shift != NULL ? p->desired_shift[a] : 0.0;
            double db = p->desired_shift != NULL ? p->desired_shift[b] : 0.0;
            double center = db - da;
            size_t ab = (size_t)a * p->ncharts + (size_t)b;
            size_t ba = (size_t)b * p->ncharts + (size_t)a;
            double lower_ab = center - p->topology_deviation_limit;
            double lower_ba = -center - p->topology_deviation_limit;
            if (lower_ab > lower[ab]) lower[ab] = lower_ab;
            if (lower_ba > lower[ba]) lower[ba] = lower_ba;
            bound_topology[ab] = bound_topology[ba] = 1;
        }
        for (size_t i = 0; i < matrix_size; i++)
            if (bound_topology[i]) stats->topology_bounds++;
    }

    Arena_T initial_work = Arena_new();
    if (initial_work == NULL) return -1;
    size_t updates = 0, rejected = 0;
    if (acr_collect_bounds(initial_work, p, &soup, p->chart_rank,
                           out_shift, p->collision_margin, lower,
                           bound_face_lo, bound_face_hi, bound_updates,
                           bound_collision, blocked_pair,
                           &updates, &rejected, &stats->exact_pairs_before,
                           &stats->exact_cross_before,
                           &stats->exact_allowed_before,
                           &stats->exact_hard_before,
                           &stats->exact_topology_limited_before, NULL) != 0) {
        Arena_dispose(&initial_work);
        return -1;
    }
    Arena_dispose(&initial_work);
    stats->collision_bounds_added += updates;
    stats->collision_bounds_rejected += rejected;
    memcpy(initial_bound, bound_collision,
           matrix_size * sizeof(*initial_bound));

    int max_outer = p->max_outer_iterations > 0
                  ? p->max_outer_iterations : 64;
    for (int outer = 0; outer < max_outer; outer++) {
        Arena_T work = Arena_new();
        if (work == NULL) return -1;
        size_t nbounds = acr_count_bounds(lower, p->ncharts);
        MonotoneQpBound *bound = (MonotoneQpBound *)ARENA_ALLOC(
            work, (nbounds ? nbounds : 1) * sizeof(*bound));
        size_t bi = 0;
        for (size_t a = 0; a < p->ncharts; a++)
        for (size_t b = 0; b < p->ncharts; b++) {
            double value = lower[a * p->ncharts + b];
            if (!(value > -DBL_MAX)) continue;
            bound[bi].lo = (int32_t)a;
            bound[bi].hi = (int32_t)b;
            bound[bi].lower = value;
            bound[bi].stroke = -1;
            bound[bi].edge = bi <= (size_t)INT32_MAX ? (int32_t)bi : INT32_MAX;
            bi++;
        }
        if (bi != nbounds) {
            Arena_dispose(&work);
            return -1;
        }
        stats->total_bounds = nbounds;
        stats->collision_bounds = 0;
        for (size_t i = 0; i < matrix_size; i++)
            if (bound_collision[i]) stats->collision_bounds++;

        for (size_t c = 0; c < p->ncharts; c++)
            out_shift[c] = p->desired_shift != NULL ? p->desired_shift[c] : 0.0;
        if (acr_make_feasible(p->ncharts, bound, nbounds, out_shift) != 0) {
            stats->topology_infeasible = 1;
            stats->collision_failed = 1;
            Arena_dispose(&work);
            return 0;
        }

        MonotoneQpProblem qp;
        memset(&qp, 0, sizeof(qp));
        qp.nvar = p->ncharts;
        qp.rows = row;
        qp.nrows = nrows;
        qp.coeff = coeff;
        qp.ncoeff = ncoeff;
        qp.bounds = bound;
        qp.nbounds = nbounds;
        MonotoneQpOptions options;
        MonotoneQpOptions_default(&options);
        size_t limit = 8 * (p->ncharts + nbounds) + 128;
        options.max_active_iterations = limit > (size_t)INT_MAX
                                      ? INT_MAX : (int)limit;
        MonotoneQpStats qstats;
        int qrc = MonotoneQp_solve(
            work, &qp, &options, out_shift, NULL, &qstats);
        if (qrc != 0) {
            stats->solve_failed = 1;
            Arena_dispose(&work);
            return 0;
        }
        if (qstats.iterations > INT_MAX - stats->qp_iterations)
            stats->qp_iterations = INT_MAX;
        else
            stats->qp_iterations += qstats.iterations;
        stats->qp_active_bounds = qstats.active_final;
        stats->qp_objective = qstats.objective_final;
        stats->qp_stationarity = qstats.stationarity_residual;
        stats->qp_max_violation = qstats.max_violation;
        stats->outer_iterations = (size_t)outer + 1;

        size_t exact = 0, cross = 0, allowed = 0, hard = 0;
        size_t topology_limited = 0;
        AtlasOverlapAudit final_audit;
        memset(&final_audit, 0, sizeof(final_audit));
        updates = rejected = 0;
        if (acr_collect_bounds(work, p, &soup, p->chart_rank,
                               out_shift, p->collision_margin, lower,
                               bound_face_lo, bound_face_hi, bound_updates,
                               bound_collision, blocked_pair,
                               &updates, &rejected, &exact, &cross,
                               &allowed, &hard, &topology_limited,
                               &final_audit) != 0) {
            Arena_dispose(&work);
            return -1;
        }
        stats->exact_pairs_after = exact;
        stats->exact_cross_after = cross;
        stats->exact_allowed_after = allowed;
        stats->exact_hard_after = hard;
        stats->exact_topology_limited_after = topology_limited;
        stats->collision_bounds_added += updates;
        stats->collision_bounds_rejected += rejected;
        if (hard == 0 ||
            (updates == 0 && hard == topology_limited)) {
            if (acr_store_reports(arena, p, &soup, lower, bound_face_lo,
                                  bound_face_hi, bound_updates,
                                  bound_topology, bound_collision,
                                  initial_bound, blocked_pair, p->chart_rank,
                                  out_shift, &final_audit, stats) != 0) {
                Arena_dispose(&work);
                return -1;
            }
            double sum2 = 0.0;
            for (size_t c = 0; c < p->ncharts; c++) {
                double a = fabs(out_shift[c]);
                sum2 += out_shift[c] * out_shift[c];
                if (a > 1.0e-10) stats->moved_charts++;
                if (a > stats->shift_max) stats->shift_max = a;
            }
            stats->shift_rms = sqrt(sum2 / (double)p->ncharts);
            Arena_dispose(&work);
            return 0;
        }
        if (updates == 0) {
            stats->collision_failed = 1;
            Arena_dispose(&work);
            return 0;
        }
        Arena_dispose(&work);
    }
    stats->collision_failed = 1;
    return 0;
}

int AtlasCollisionRegister_selftest(void)
{
    const int32_t faces[6] = {0, 1, 2, 3, 4, 5};
    const double u[6] = {0, 2, 0, 0, 2, 0};
    const double v[6] = {0, 0, 2, 0, 0, 2};
    const int32_t chart[6] = {0, 0, 0, 1, 1, 1};
    const size_t rank[2] = {0, 1};
    const double desired[2] = {0, 0};
    AtlasRegisterEdge edge;
    memset(&edge, 0, sizeof(edge));
    edge.component0 = 0;
    edge.component1 = 1;
    edge.target = 0.0;
    edge.weight = 10.0;
    edge.kind = ATLAS_REGISTER_WELD;
    AtlasCollisionRegisterProblem p;
    memset(&p, 0, sizeof(p));
    p.faces = faces;
    p.nfaces = 2;
    p.nvertices = 6;
    p.base_u = u;
    p.v = v;
    p.vertex_chart = chart;
    p.ncharts = 2;
    p.chart_rank = rank;
    p.desired_shift = desired;
    p.edges = &edge;
    p.nedges = 1;
    p.collision_margin = 1.0e-5;
    p.max_outer_iterations = 8;
    Arena_T arena = Arena_new();
    if (arena == NULL) return -1;
    double shift[2];
    AtlasCollisionRegisterStats stats;
    int rc = AtlasCollisionRegister_solve(arena, &p, shift, &stats);
    int failed = rc != 0 || stats.solve_failed || stats.collision_failed ||
                 stats.exact_cross_before != 1 ||
                 stats.exact_hard_before != 1 ||
                 stats.exact_cross_after != 0 ||
                 stats.exact_hard_after != 0 ||
                 shift[1] - shift[0] < 2.0;
    Arena_dispose(&arena);

    if (failed) return -1;

    /* The identical geometric situation is deliberately retained when the
     * charts are known observations of one sheet.  It must be reported, but
     * must not create a separating bound or a metric jump. */
    const uint8_t allowed_pair[4] = {0, 1, 1, 0};
    const float near_xyz[18] = {
        0, 0, 0,  2, 0, 0,  0, 2, 0,
        0, 0, 0,  2, 0, 0,  0, 2, 0
    };
    p.xyz = near_xyz;
    p.same_sheet_xyz_tolerance = 8.0;
    p.allowed_chart_pair = allowed_pair;
    arena = Arena_new();
    if (arena == NULL) return -1;
    rc = AtlasCollisionRegister_solve(arena, &p, shift, &stats);
    failed = rc != 0 || stats.solve_failed || stats.collision_failed ||
             stats.exact_cross_before != 1 ||
             stats.exact_allowed_before != 1 ||
             stats.exact_hard_before != 0 ||
             stats.exact_cross_after != 1 ||
             stats.exact_allowed_after != 1 ||
             stats.exact_hard_after != 0 ||
             stats.collision_bounds != 0 ||
             stats.nresidual != 1 ||
             !stats.residual[0].allowed ||
             stats.residual[0].reason != 1 ||
             fabs(shift[0]) > 1.0e-9 || fabs(shift[1]) > 1.0e-9;
    Arena_dispose(&arena);

    if (failed) return -1;

    /* A direct topology edge does not exempt remote colliding faces.  Topology
     * still has priority when separating them would exceed the seam budget. */
    const float remote_xyz[18] = {
        0, 0, 0,    2, 0, 0,    0, 2, 0,
        0, 0, 100,  2, 0, 100,  0, 2, 100
    };
    p.xyz = remote_xyz;
    p.topology_deviation_limit = 0.25;
    arena = Arena_new();
    if (arena == NULL) return -1;
    rc = AtlasCollisionRegister_solve(arena, &p, shift, &stats);
    failed = rc != 0 || stats.solve_failed || stats.collision_failed ||
             stats.topology_infeasible || stats.topology_bounds != 2 ||
             stats.collision_bounds != 0 ||
             stats.collision_bounds_rejected == 0 ||
             stats.exact_hard_before != 1 ||
             stats.exact_topology_limited_before != 1 ||
             stats.exact_hard_after != 1 ||
             stats.exact_topology_limited_after != 1 ||
             stats.nresidual != 1 ||
             stats.residual[0].allowed ||
             stats.residual[0].reason != 4 ||
             stats.residual[0].min_xyz < 99.0 ||
             fabs(shift[0]) > 1.0e-9 || fabs(shift[1]) > 1.0e-9;
    Arena_dispose(&arena);
    return failed ? -1 : 0;
}
