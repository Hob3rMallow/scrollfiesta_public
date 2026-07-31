#include "atlas_candidates.h"

#include "../common/kdtree.h"
#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define AC_MAX_CANDIDATES 16
#define AC_BALL_CAP 256

typedef struct {
    int32_t edge[2];
    int32_t plane;
    int32_t cube;
    double mesh_t;
    double p[3];
    double raw_u;
    int32_t adjacent[2];
    int32_t degree;
} SliceNode;

typedef struct {
    int32_t node[2];
    int32_t face;
    uint8_t visited;
} SliceSegment;

typedef struct {
    int32_t edge0;
    int32_t edge1;
    int32_t plane;
    int32_t node;
} EdgePlaneSlot;

typedef struct {
    int32_t a;
    int32_t b;
    int32_t stroke;
    int32_t plane;
} CoarseEdge;

typedef struct {
    int32_t edge;
    int32_t stroke;
    double t;
    double lateral2;
    double score;
} CandidateHit;

static double ac_dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double ac_normalize3(double v[3])
{
    double n = sqrt(ac_dot3(v, v));
    if (n > 1e-15) {
        v[0] /= n;
        v[1] /= n;
        v[2] /= n;
    }
    return n;
}

static double ac_dist3(const double a[3], const double b[3])
{
    double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return sqrt(ac_dot3(d, d));
}

/*
 * Pinned spiral geometry.  The in-plane frame is byte-identical to the
 * driver's reference maps (deterministic tangent seed, f1 = axis x f0), so
 * u_geo here equals the driver's radius-anchored expected_geo field.
 */
typedef struct {
    int enabled;
    double pitch;
    double spiral_a;
    double spiral_b;
    double sense;
    double tol;
    double axis[3];
    double f0[3];
    double f1[3];
    double origin[3];
} AcSpiral;

static void ac_spiral_init(AcSpiral *sp, const AtlasCandidateOptions *o,
                           const double axis[3], const float axis_point[3])
{
    memset(sp, 0, sizeof *sp);
    if (o->pitch <= 0.0) return;
    sp->enabled = 1;
    sp->pitch = o->pitch;
    sp->spiral_a = o->spiral_a;
    sp->spiral_b = o->spiral_b;
    sp->sense = o->sense >= 0 ? 1.0 : -1.0;
    sp->tol = o->wind_tol;
    memcpy(sp->axis, axis, sizeof sp->axis);
    for (int d = 0; d < 3; d++) sp->origin[d] = (double)axis_point[d];
    double t[3] = {1.0, 0.0, 0.0};
    if (fabs(axis[0]) > 0.9) { t[0] = 0.0; t[1] = 1.0; }
    double proj = ac_dot3(t, axis);
    for (int d = 0; d < 3; d++) sp->f0[d] = t[d] - proj * axis[d];
    ac_normalize3(sp->f0);
    sp->f1[0] = axis[1] * sp->f0[2] - axis[2] * sp->f0[1];
    sp->f1[1] = axis[2] * sp->f0[0] - axis[0] * sp->f0[2];
    sp->f1[2] = axis[0] * sp->f0[1] - axis[1] * sp->f0[0];
}

static void ac_spiral_polar(const AcSpiral *sp, const double p[3],
                            double *out_r, double *out_theta)
{
    double q[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < 3; d++) q[d] = p[d] - sp->origin[d];
    double along = ac_dot3(q, sp->axis);
    double pa = 0.0, pb = 0.0;
    for (int d = 0; d < 3; d++) {
        double w = q[d] - along * sp->axis[d];
        pa += w * sp->f0[d];
        pb += w * sp->f1[d];
    }
    *out_r = sqrt(pa * pa + pb * pb);
    *out_theta = atan2(pb, pa);
}

static double ac_spiral_wind(const AcSpiral *sp, double r, double theta)
{
    return sp->sense * r / sp->pitch - theta / (2.0 * M_PI);
}

static double ac_spiral_u_geo(const AcSpiral *sp, double r, double theta)
{
    double w = ac_spiral_wind(sp, r, theta);
    double phi = theta + 2.0 * M_PI * floor(w + 0.5);
    return sp->spiral_a * phi +
           sp->spiral_b * phi * phi / (4.0 * M_PI);
}

/*
 * Branch-cut-free pairwise winding delta -- the seam gate's own test.  Local
 * differentials only: axis wander is common-mode, pitch error is second
 * order at relation scales, so this stays trustworthy where the absolute
 * winding does not.  ~0 for a same-wrap pair by the spiral property, ~+-1
 * for a cross-wrap short-circuit.
 */
static double ac_spiral_dwind(const AcSpiral *sp,
                              double r0, double th0,
                              double r1, double th1)
{
    double dth = th1 - th0;
    while (dth > M_PI) dth -= 2.0 * M_PI;
    while (dth < -M_PI) dth += 2.0 * M_PI;
    return sp->sense * (r1 - r0) / sp->pitch - dth / (2.0 * M_PI);
}

static int ac_wind_bin(double dw)
{
    double a = fabs(dw);
    if (a < 0.1) return 0;
    if (a < 0.35) return 1;
    if (a < 0.65) return 2;
    if (a < 1.5) return 3;
    return 4;
}

static int ac_wind_fusing(double dw)
{
    double m = floor(dw + 0.5);
    return m != 0.0 && fabs(dw - m) <= 0.25;
}

void AtlasCandidateOptions_default(AtlasCandidateOptions *opts)
{
    if (opts == NULL) return;
    opts->slice_spacing = 8.0;
    opts->min_stroke_length = 6.0;
    opts->sample_spacing = 4.0;
    opts->match_radius = 4.0;
    opts->match_angle_deg = 25.0;
    opts->max_slice_stride = 3;
    opts->max_candidates = 4;
    opts->continuation_radius = 4.0;
    opts->continuation_angle_deg = 30.0;
    opts->pitch = 0.0;          /* <= 0: every geometric mechanism off */
    opts->spiral_a = 0.0;
    opts->spiral_b = 0.0;
    opts->sense = -1;
    opts->wind_gate = 1;
    opts->wind_tol = 0.35;
    opts->geo_gauge = 1;
    opts->wind_cut = 1;
}

void AtlasCandidateRefineOptions_default(AtlasCandidateRefineOptions *opts)
{
    if (opts == NULL) return;
    memset(opts, 0, sizeof *opts);
    opts->support_floor = 1e-4;
    opts->gross_initial_residual = 512.0;
    opts->initial_coherence_mad_limit = 16.0;
    opts->overlap_residual_limit = 4.0;
    opts->overlap_coherence_mad_limit = 4.0;
    opts->persistent_contact_fraction = 0.75;
    opts->minimum_topology_alternative_fraction = 0.50;
    opts->ambiguity_ratio_limit = 0.65;
    opts->minimum_prior = 0.05;
    opts->minimum_top_prior_fraction = 0.60;
    opts->minimum_arc_correlation = 0.90;
    opts->minimum_links = 8;
    opts->minimum_planes = 3;
    opts->front_minimum_links = 8;
    opts->front_minimum_planes = 3;
    opts->front_maximum_extension = 2;
    opts->front_close_gap = 1.5;
    opts->front_minimum_contact_overlap = 0.50;
}

static int ac_validate_options(const AtlasCandidateOptions *o)
{
    if (o == NULL) return 0;
    if (o->pitch > 0.0 &&
        !(isfinite(o->pitch) &&
          isfinite(o->spiral_a) && isfinite(o->spiral_b) &&
          (o->sense == 1 || o->sense == -1) &&
          isfinite(o->wind_tol) &&
          o->wind_tol > 0.0 && o->wind_tol <= 0.5))
        return 0;
    return isfinite(o->slice_spacing) && o->slice_spacing > 0.0 &&
           isfinite(o->min_stroke_length) && o->min_stroke_length >= 0.0 &&
           isfinite(o->sample_spacing) && o->sample_spacing > 0.0 &&
           isfinite(o->match_radius) && o->match_radius > 0.0 &&
           isfinite(o->match_angle_deg) &&
           o->match_angle_deg >= 0.0 && o->match_angle_deg < 90.0 &&
           o->max_slice_stride >= 1 && o->max_slice_stride <= 16 &&
           o->max_candidates >= 1 &&
           o->max_candidates <= AC_MAX_CANDIDATES &&
           isfinite(o->continuation_radius) &&
           o->continuation_radius > 0.0 &&
           isfinite(o->continuation_angle_deg) &&
           o->continuation_angle_deg >= 0.0 &&
           o->continuation_angle_deg < 90.0;
}

/*
 * Global planes are t_k=(k+1/2)h.  Strict-interior intersection makes a
 * vertex exactly on a plane consistently belong to neither incident edge,
 * avoiding zero-length duplicate segments.
 */
static void ac_plane_range(double tlo, double thi, double h,
                           int32_t *out_k0, int32_t *out_k1)
{
    double x0 = tlo / h - 0.5;
    double x1 = thi / h - 0.5;
    double a = floor(x0) + 1.0;
    double b = ceil(x1) - 1.0;
    if (a < (double)INT32_MIN) a = (double)INT32_MIN;
    if (a > (double)INT32_MAX) a = (double)INT32_MAX;
    if (b < (double)INT32_MIN) b = (double)INT32_MIN;
    if (b > (double)INT32_MAX) b = (double)INT32_MAX;
    *out_k0 = (int32_t)a;
    *out_k1 = (int32_t)b;
}

static size_t ac_next_pow2(size_t value)
{
    size_t p = 1;
    while (p < value && p <= SIZE_MAX / 2) p <<= 1;
    return p;
}

static uint64_t ac_hash_edge_plane(int32_t a, int32_t b, int32_t plane)
{
    uint64_t x = (uint64_t)(uint32_t)a * 0x9e3779b185ebca87ull;
    x ^= (uint64_t)(uint32_t)b * 0xc2b2ae3d27d4eb4full;
    x ^= (uint64_t)(uint32_t)plane * 0x165667b19e3779f9ull;
    x ^= x >> 29;
    x *= 0x9fb21c651e98df25ull;
    return x ^ (x >> 32);
}

static int32_t ac_get_node(EdgePlaneSlot *table, size_t mask,
                           SliceNode *nodes, size_t *nnode,
                           int32_t va, int32_t vb, int32_t plane,
                           int32_t cube, const float *vertices,
                           const float *vertex_u, const double *axial,
                           double spacing)
{
    int32_t lo = va < vb ? va : vb;
    int32_t hi = va < vb ? vb : va;
    size_t slot = (size_t)ac_hash_edge_plane(lo, hi, plane) & mask;
    while (table[slot].node >= 0) {
        if (table[slot].edge0 == lo && table[slot].edge1 == hi &&
            table[slot].plane == plane)
            return table[slot].node;
        slot = (slot + 1) & mask;
    }

    double tp = ((double)plane + 0.5) * spacing;
    double den = axial[hi] - axial[lo];
    double t = fabs(den) > 1e-15 ? (tp - axial[lo]) / den : 0.5;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    if (*nnode > (size_t)INT32_MAX) return -1;
    int32_t id = (int32_t)(*nnode);
    SliceNode *node = &nodes[(*nnode)++];
    memset(node, 0, sizeof *node);
    node->edge[0] = lo;
    node->edge[1] = hi;
    node->plane = plane;
    node->cube = cube;
    node->mesh_t = t;
    node->adjacent[0] = node->adjacent[1] = -1;
    for (int d = 0; d < 3; d++)
        node->p[d] = (1.0 - t) * (double)vertices[(size_t)lo * 3 + (size_t)d] +
                     t * (double)vertices[(size_t)hi * 3 + (size_t)d];
    node->raw_u = (1.0 - t) * (double)vertex_u[lo] +
                  t * (double)vertex_u[hi];

    table[slot].edge0 = lo;
    table[slot].edge1 = hi;
    table[slot].plane = plane;
    table[slot].node = id;
    return id;
}

static size_t ac_walk_chain(SliceSegment *segments, const SliceNode *nodes,
                            int32_t start_node, int32_t start_segment,
                            int32_t *sequence, size_t capacity,
                            int *out_closed)
{
    size_t n = 0;
    int32_t node = start_node;
    int32_t segment = start_segment;
    *out_closed = 0;

    while (segment >= 0 && !segments[segment].visited && n < capacity) {
        sequence[n++] = node;
        segments[segment].visited = 1;
        int32_t next = segments[segment].node[0] == node
                     ? segments[segment].node[1]
                     : segments[segment].node[0];
        int32_t next_segment = -1;
        int ncheck = nodes[next].degree < 2 ? nodes[next].degree : 2;
        for (int j = 0; j < ncheck; j++) {
            int32_t candidate = nodes[next].adjacent[j];
            if (candidate >= 0 && candidate != segment &&
                !segments[candidate].visited) {
                next_segment = candidate;
                break;
            }
        }
        if (next_segment < 0) {
            if (next == start_node)
                *out_closed = 1;
            else if (n < capacity)
                sequence[n++] = next;
            break;
        }
        node = next;
        segment = next_segment;
    }
    return n;
}

static double ac_sample_dual_width(const AtlasStripSample *samples,
                                   const AtlasStripStroke *strokes,
                                   int32_t sample)
{
    const AtlasStripSample *s = &samples[sample];
    const AtlasStripStroke *st = &strokes[s->stroke];
    double before = 0.0, after = 0.0;
    if (s->ordinal > 0)
        before = s->s - samples[sample - 1].s;
    if (s->ordinal + 1 < st->count)
        after = samples[sample + 1].s - s->s;
    double width = 0.5 * (before + after);
    if (before == 0.0) width = 0.5 * after;
    if (after == 0.0) width = 0.5 * before;
    return width > 1e-6 ? width : 1.0;
}

static void ac_member_direct(AtlasStripMember *m,
                             const AtlasStripSample *samples,
                             const AtlasStripStroke *strokes,
                             int32_t sample, double membership,
                             int32_t observation)
{
    memset(m, 0, sizeof *m);
    const AtlasStripSample *s = &samples[sample];
    const AtlasStripStroke *st = &strokes[s->stroke];
    int32_t lo, hi;
    if (s->ordinal + 1 < st->count) {
        lo = sample;
        hi = sample + 1;
    } else {
        lo = sample - 1;
        hi = sample;
    }
    m->value0 = sample;
    m->value1 = -1;
    m->deriv_lo = lo;
    m->deriv_hi = hi;
    m->deriv_length = samples[hi].s - samples[lo].s;
    memcpy(m->p, s->p, sizeof m->p);
    memcpy(m->tangent, s->tangent, sizeof m->tangent);
    m->dual_width = ac_sample_dual_width(samples, strokes, sample);
    m->membership = membership;
    m->base_weight = 1.0;
    m->observation = observation;
}

static void ac_member_edge(AtlasStripMember *m,
                           const AtlasStripSample *samples,
                           const CoarseEdge *edge, double t,
                           double membership, int32_t observation)
{
    memset(m, 0, sizeof *m);
    m->value0 = edge->a;
    m->value1 = edge->b;
    m->value_t = t;
    m->deriv_lo = edge->a;
    m->deriv_hi = edge->b;
    m->deriv_length = samples[edge->b].s - samples[edge->a].s;
    for (int d = 0; d < 3; d++) {
        m->p[d] = (1.0 - t) * samples[edge->a].p[d] +
                  t * samples[edge->b].p[d];
        m->tangent[d] = (1.0 - t) * samples[edge->a].tangent[d] +
                        t * samples[edge->b].tangent[d];
    }
    if (ac_normalize3(m->tangent) < 1e-12)
        memcpy(m->tangent, samples[edge->a].tangent, sizeof m->tangent);
    m->dual_width = m->deriv_length;
    m->membership = membership;
    m->base_weight = 1.0;
    m->observation = observation;
}

static void ac_keep_candidate(CandidateHit *hits, int *nhit, int capacity,
                              const CandidateHit *candidate)
{
    for (int i = 0; i < *nhit; i++) {
        if (hits[i].stroke != candidate->stroke) continue;
        if (candidate->score > hits[i].score) hits[i] = *candidate;
        return;
    }
    if (*nhit < capacity) {
        hits[(*nhit)++] = *candidate;
        return;
    }
    int worst = 0;
    for (int i = 1; i < *nhit; i++)
        if (hits[i].score < hits[worst].score) worst = i;
    if (candidate->score > hits[worst].score) hits[worst] = *candidate;
}

static int ac_compare_candidate(const void *pa, const void *pb)
{
    const CandidateHit *a = (const CandidateHit *)pa;
    const CandidateHit *b = (const CandidateHit *)pb;
    return a->score > b->score ? -1 : (a->score < b->score ? 1 : 0);
}

int AtlasCandidates_build(Arena_T arena,
                          const float *vertices,
                          size_t nvertices,
                          const int32_t *triangles,
                          const int32_t *face_cube,
                          size_t ntriangles,
                          const float *vertex_u,
                          const float axis_point[3],
                          const float axis_dir[3],
                          const AtlasCandidateOptions *input_opts,
                          AtlasCandidateSet *out)
{
    if (arena == NULL || vertices == NULL || nvertices < 3 ||
        nvertices > (size_t)INT32_MAX || triangles == NULL ||
        ntriangles == 0 || vertex_u == NULL || axis_point == NULL ||
        axis_dir == NULL || out == NULL)
        return -1;

    AtlasCandidateOptions defaults;
    AtlasCandidateOptions_default(&defaults);
    const AtlasCandidateOptions *o = input_opts != NULL ? input_opts : &defaults;
    if (!ac_validate_options(o)) return -1;
    memset(out, 0, sizeof *out);
    out->stats.input_faces = ntriangles;

    double axis[3] = {
        (double)axis_dir[0], (double)axis_dir[1], (double)axis_dir[2]
    };
    if (ac_normalize3(axis) < 1e-12) return -1;
    for (int d = 0; d < 3; d++)
        if (!isfinite((double)axis_point[d]) ||
            !isfinite((double)axis_dir[d])) return -1;

    double *axial = (double *)ARENA_ALLOC(
        arena, nvertices * sizeof(double));
    for (size_t i = 0; i < nvertices; i++) {
        double q[3];
        for (int d = 0; d < 3; d++) {
            double v = (double)vertices[i * 3 + (size_t)d];
            if (!isfinite(v) || !isfinite((double)vertex_u[i])) return -1;
            q[d] = v - (double)axis_point[d];
        }
        axial[i] = ac_dot3(q, axis);
    }

    AcSpiral sp;
    ac_spiral_init(&sp, o, axis, axis_point);

    /* Retained-triangle connectivity is the hard evidence used by the first
     * pass.  Cube identity is not enough: one cube can contain several real
     * sheets, while a retained component is already a certified surface. */
    UnionFind mesh_topology = UF_new(arena, (int32_t)nvertices);
    uint8_t *vertex_used = (uint8_t *)ARENA_CALLOC(
        arena, nvertices, sizeof(uint8_t));
    size_t segment_capacity = 0;
    for (size_t f = 0; f < ntriangles; f++) {
        int32_t v[3] = {
            triangles[f * 3 + 0],
            triangles[f * 3 + 1],
            triangles[f * 3 + 2]
        };
        for (int j = 0; j < 3; j++)
            if (v[j] < 0 || (size_t)v[j] >= nvertices) return -1;
        vertex_used[v[0]] = vertex_used[v[1]] = vertex_used[v[2]] = 1;
        uf_union(&mesh_topology, v[0], v[1]);
        uf_union(&mesh_topology, v[1], v[2]);
        double lo = axial[v[0]], hi = axial[v[0]];
        for (int j = 1; j < 3; j++) {
            if (axial[v[j]] < lo) lo = axial[v[j]];
            if (axial[v[j]] > hi) hi = axial[v[j]];
        }
        int32_t k0, k1;
        ac_plane_range(lo, hi, o->slice_spacing, &k0, &k1);
        if (k1 >= k0) {
            size_t add = (size_t)((int64_t)k1 - (int64_t)k0 + 1);
            if (segment_capacity > SIZE_MAX - add) return -1;
            segment_capacity += add;
        }
    }
    if (segment_capacity == 0 ||
        segment_capacity > (size_t)INT32_MAX / 2) return -1;

    int32_t *mesh_component = (int32_t *)ARENA_ALLOC(
        arena, nvertices * sizeof(int32_t));
    int32_t *mesh_root_label = (int32_t *)ARENA_ALLOC(
        arena, nvertices * sizeof(int32_t));
    for (size_t i = 0; i < nvertices; i++) {
        mesh_component[i] = -1;
        mesh_root_label[i] = -1;
    }
    size_t nmesh_component = 0;
    for (size_t i = 0; i < nvertices; i++) {
        if (!vertex_used[i]) continue;
        int32_t root = uf_find(&mesh_topology, (int32_t)i);
        if (mesh_root_label[root] < 0) {
            if (nmesh_component > (size_t)INT32_MAX) return -1;
            mesh_root_label[root] = (int32_t)nmesh_component++;
        }
        mesh_component[i] = mesh_root_label[root];
    }

    SliceNode *nodes = (SliceNode *)ARENA_ALLOC(
        arena, 2 * segment_capacity * sizeof(SliceNode));
    SliceSegment *segments = (SliceSegment *)ARENA_ALLOC(
        arena, segment_capacity * sizeof(SliceSegment));
    size_t table_need = segment_capacity > SIZE_MAX / 4
                      ? SIZE_MAX : 4 * segment_capacity;
    size_t table_cap = ac_next_pow2(table_need < 16 ? 16 : table_need);
    if (table_cap < table_need) return -1;
    EdgePlaneSlot *table = (EdgePlaneSlot *)ARENA_ALLOC(
        arena, table_cap * sizeof(EdgePlaneSlot));
    for (size_t i = 0; i < table_cap; i++) table[i].node = -1;

    size_t nnode = 0, nsegment = 0;
    int32_t plane_min = INT32_MAX, plane_max = INT32_MIN;
    for (size_t f = 0; f < ntriangles; f++) {
        int32_t fv[3] = {
            triangles[f * 3 + 0],
            triangles[f * 3 + 1],
            triangles[f * 3 + 2]
        };
        double lo = axial[fv[0]], hi = axial[fv[0]];
        for (int j = 1; j < 3; j++) {
            if (axial[fv[j]] < lo) lo = axial[fv[j]];
            if (axial[fv[j]] > hi) hi = axial[fv[j]];
        }
        int32_t k0, k1;
        ac_plane_range(lo, hi, o->slice_spacing, &k0, &k1);
        for (int64_t kk = (int64_t)k0; kk <= (int64_t)k1; kk++) {
            int32_t k = (int32_t)kk;
            double tp = ((double)k + 0.5) * o->slice_spacing;
            int32_t crossing[3];
            int ncross = 0;
            for (int e = 0; e < 3; e++) {
                int32_t a = fv[e], b = fv[(e + 1) % 3];
                double ea = axial[a], eb = axial[b];
                double elo = ea < eb ? ea : eb;
                double ehi = ea < eb ? eb : ea;
                if (!(elo < tp && tp < ehi)) continue;
                int32_t node = ac_get_node(
                    table, table_cap - 1, nodes, &nnode, a, b, k,
                    face_cube != NULL ? face_cube[f] : 0,
                    vertices, vertex_u, axial, o->slice_spacing);
                if (node < 0) return -1;
                if (ncross < 3) crossing[ncross] = node;
                ncross++;
            }
            if (ncross != 2 || crossing[0] == crossing[1]) continue;
            if (nsegment >= segment_capacity) return -1;
            segments[nsegment].node[0] = crossing[0];
            segments[nsegment].node[1] = crossing[1];
            segments[nsegment].face = (int32_t)f;
            segments[nsegment].visited = 0;
            nsegment++;
            if (k < plane_min) plane_min = k;
            if (k > plane_max) plane_max = k;
        }
    }
    if (nsegment == 0 || nnode < 2) return -1;

    for (size_t s = 0; s < nsegment; s++) {
        for (int endpoint = 0; endpoint < 2; endpoint++) {
            SliceNode *node = &nodes[segments[s].node[endpoint]];
            if (node->degree < 2)
                node->adjacent[node->degree] = (int32_t)s;
            node->degree++;
        }
    }
    for (size_t i = 0; i < nnode; i++)
        if (nodes[i].degree > 2) out->stats.branch_nodes++;

    size_t sample_capacity = 2 * nsegment + 2;
    AtlasStripSample *samples = (AtlasStripSample *)ARENA_ALLOC(
        arena, sample_capacity * sizeof(AtlasStripSample));
    AtlasCandidateSampleRef *sample_ref =
        (AtlasCandidateSampleRef *)ARENA_ALLOC(
            arena, sample_capacity * sizeof(AtlasCandidateSampleRef));
    double *initial_u = (double *)ARENA_ALLOC(
        arena, sample_capacity * sizeof(double));
    AtlasStripStroke *strokes = (AtlasStripStroke *)ARENA_ALLOC(
        arena, (nsegment + 1) * sizeof(AtlasStripStroke));
    int32_t *sequence = (int32_t *)ARENA_ALLOC(
        arena, (nsegment + 2) * sizeof(int32_t));
    int32_t *compact = (int32_t *)ARENA_ALLOC(
        arena, (nsegment + 2) * sizeof(int32_t));
    double *selected_s = (double *)ARENA_ALLOC(
        arena, (nsegment + 2) * sizeof(double));
    int32_t *selected_k = (int32_t *)ARENA_ALLOC(
        arena, (nsegment + 2) * sizeof(int32_t));
    int32_t *seg_off = (int32_t *)ARENA_ALLOC(
        arena, (nsegment + 3) * sizeof(int32_t));
    double *chain_r = NULL, *chain_th = NULL, *chain_ugeo = NULL;
    double *geo_u = NULL;
    if (sp.enabled) {
        chain_r = (double *)ARENA_ALLOC(
            arena, (nsegment + 2) * sizeof(double));
        chain_th = (double *)ARENA_ALLOC(
            arena, (nsegment + 2) * sizeof(double));
        chain_ugeo = (double *)ARENA_ALLOC(
            arena, (nsegment + 2) * sizeof(double));
        geo_u = (double *)ARENA_ALLOC(
            arena, sample_capacity * sizeof(double));
    }

    size_t nsample = 0, nstroke = 0;
    for (int pass = 0; pass < 2; pass++) {
        for (size_t ni = 0; ni < nnode; ni++) {
            int ncheck = nodes[ni].degree < 2 ? nodes[ni].degree : 2;
            if (ncheck == 0) continue;
            if (pass == 0 && nodes[ni].degree != 1) continue;
            for (int aj = 0; aj < ncheck; aj++) {
                int32_t start_segment = nodes[ni].adjacent[aj];
                if (start_segment < 0 || segments[start_segment].visited)
                    continue;
                int closed = 0;
                size_t nseq = ac_walk_chain(
                    segments, nodes, (int32_t)ni, start_segment,
                    sequence, nsegment + 2, &closed);
                if (nseq < 2) continue;

                size_t nc = 0;
                for (size_t j = 0; j < nseq; j++) {
                    if (nc > 0 &&
                        ac_dist3(nodes[compact[nc - 1]].p,
                                 nodes[sequence[j]].p) <= 1e-8)
                        continue;
                    compact[nc++] = sequence[j];
                }
                if (nc < 2) continue;

                /*
                 * Wrap-crossing cuts.  On a coarse mesh whose interior edges
                 * exceed the pitch, a plane-slice chain can snake across
                 * wraps through in-mesh fusions; one such stroke then owes
                 * whole circumferences while carrying only its own
                 * arclength -- the collapse channel no relation gate can
                 * reach.  Segment the chain wherever the PAIRWISE winding
                 * accumulated from the segment anchor exceeds the
                 * tolerance; the continuation gate below refuses to re-join
                 * the halves across the cut (dwind ~ 1 there), while an
                 * accidental same-wrap cut is harmlessly re-glued
                 * (dwind ~ 0).
                 */
                size_t nseg_here = 1;
                seg_off[0] = 0;
                seg_off[1] = (int32_t)nc;
                if (sp.enabled) {
                    double wc = 0.0, w_anchor = 0.0;
                    for (size_t j = 0; j < nc; j++) {
                        ac_spiral_polar(&sp, nodes[compact[j]].p,
                                        &chain_r[j], &chain_th[j]);
                        chain_ugeo[j] = ac_spiral_u_geo(
                            &sp, chain_r[j], chain_th[j]);
                        if (j > 0) {
                            wc += ac_spiral_dwind(
                                &sp, chain_r[j - 1], chain_th[j - 1],
                                chain_r[j], chain_th[j]);
                            if (o->wind_cut &&
                                fabs(wc - w_anchor) > sp.tol) {
                                seg_off[nseg_here++] = (int32_t)j;
                                w_anchor = wc;
                            }
                        }
                    }
                    seg_off[nseg_here] = (int32_t)nc;
                    if (nseg_here > 1) {
                        out->stats.chains_wind_split++;
                        out->stats.chain_wind_cuts += nseg_here - 1;
                    }
                }

                for (size_t segi = 0; segi < nseg_here; segi++) {
                size_t c0 = (size_t)seg_off[segi];
                size_t ncs = (size_t)seg_off[segi + 1] - c0;
                int32_t *cs = compact + c0;
                double *sr = sp.enabled ? chain_r + c0 : NULL;
                double *sth = sp.enabled ? chain_th + c0 : NULL;
                double *sug = sp.enabled ? chain_ugeo + c0 : NULL;
                if (ncs < 2) continue;

                double total_length = 0.0;
                for (size_t j = 1; j < ncs; j++)
                    total_length += ac_dist3(nodes[cs[j - 1]].p,
                                             nodes[cs[j]].p);
                if (total_length < o->min_stroke_length) {
                    out->stats.short_strokes_dropped++;
                    continue;
                }

                /* Winding span of the EMITTED stroke; with the cut enabled
                 * this is the verification that segmentation worked. */
                int trusted = 0;
                double mean_r = 0.0;
                if (sp.enabled) {
                    double wc = 0.0, wc_min = 0.0, wc_max = 0.0;
                    for (size_t j = 0; j < ncs; j++) {
                        mean_r += sr[j];
                        if (j > 0) {
                            wc += ac_spiral_dwind(
                                &sp, sr[j - 1], sth[j - 1],
                                sr[j], sth[j]);
                            if (wc < wc_min) wc_min = wc;
                            if (wc > wc_max) wc_max = wc;
                        }
                    }
                    mean_r /= (double)ncs;
                    double span = wc_max - wc_min;
                    if (span > out->stats.stroke_wind_span_max)
                        out->stats.stroke_wind_span_max = span;
                    if (span >= 0.5) out->stats.strokes_wind_bridge++;
                }

                double mean_s = 0.5 * total_length, mean_u = 0.0;
                int use_geo = sp.enabled && o->geo_gauge;
                for (size_t j = 0; j < ncs; j++)
                    mean_u += use_geo ? sug[j] : nodes[cs[j]].raw_u;
                mean_u /= (double)ncs;
                double covariance = 0.0, cumulative = 0.0;
                for (size_t j = 0; j < ncs; j++) {
                    if (j > 0)
                        cumulative += ac_dist3(nodes[cs[j - 1]].p,
                                               nodes[cs[j]].p);
                    double uref = use_geo ? sug[j] : nodes[cs[j]].raw_u;
                    covariance += (cumulative - mean_s) * (uref - mean_u);
                }
                if (covariance < 0.0) {
                    for (size_t j = 0; j < ncs / 2; j++) {
                        int32_t tmp = cs[j];
                        cs[j] = cs[ncs - 1 - j];
                        cs[ncs - 1 - j] = tmp;
                        if (sp.enabled) {
                            double swap = sr[j];
                            sr[j] = sr[ncs - 1 - j];
                            sr[ncs - 1 - j] = swap;
                            swap = sth[j];
                            sth[j] = sth[ncs - 1 - j];
                            sth[ncs - 1 - j] = swap;
                            swap = sug[j];
                            sug[j] = sug[ncs - 1 - j];
                            sug[ncs - 1 - j] = swap;
                        }
                    }
                }

                /*
                 * Prior trust: the ABSOLUTE winding is only usable where
                 * the local spiral model demonstrably fits, i.e. u_geo
                 * tracks arclength along the oriented chain.  Round-flips
                 * (residual bridges, eccentric lobes, mis-centring) blow
                 * the deviation up by a circumference and silence the
                 * prior for the whole stroke; such fragments are placed by
                 * evidence or reported ambiguous, never guessed.
                 */
                if (sp.enabled) {
                    double dev_max = 0.0, cum = 0.0;
                    for (size_t j = 0; j < ncs; j++) {
                        if (j > 0)
                            cum += ac_dist3(nodes[cs[j - 1]].p,
                                            nodes[cs[j]].p);
                        double dev = fabs(sug[j] - sug[0] - cum);
                        if (dev > dev_max) dev_max = dev;
                    }
                    double circumference = 2.0 * M_PI * mean_r;
                    trusted = circumference > 0.0 &&
                              dev_max < sp.tol * circumference;
                    if (!trusted) out->stats.strokes_prior_untrusted++;
                }

                /* StrokeStrip is a coarse solve.  Keep exact mesh-edge
                 * samples at approximately uniform arclength, always
                 * retaining both endpoints; the full triangle mesh returns
                 * later through the FEM extension. */
                size_t nselected = 1;
                sequence[0] = cs[0];
                selected_s[0] = 0.0;
                selected_k[0] = 0;
                cumulative = 0.0;
                double last_selected_s = 0.0;
                for (size_t j = 1; j < ncs; j++) {
                    cumulative += ac_dist3(nodes[cs[j - 1]].p,
                                           nodes[cs[j]].p);
                    if (j + 1 == ncs ||
                        cumulative - last_selected_s >= o->sample_spacing) {
                        sequence[nselected] = cs[j];
                        selected_s[nselected] = cumulative;
                        selected_k[nselected] = (int32_t)j;
                        last_selected_s = cumulative;
                        nselected++;
                    }
                }
                if (nselected < 2) continue;

                if (nsample + nselected > sample_capacity ||
                    nstroke > (size_t)INT32_MAX) return -1;
                AtlasStripStroke *stroke = &strokes[nstroke];
                memset(stroke, 0, sizeof *stroke);
                stroke->first = (int32_t)nsample;
                stroke->count = (int32_t)nselected;
                stroke->component = nodes[sequence[0]].cube;
                stroke->resolved = 1;

                double gauge = use_geo ? sug[0]
                                       : nodes[sequence[0]].raw_u;
                for (size_t j = 0; j < nselected; j++) {
                    cumulative = selected_s[j];
                    int32_t node_id = sequence[j];
                    AtlasStripSample *sample = &samples[nsample + j];
                    memset(sample, 0, sizeof *sample);
                    memcpy(sample->p, nodes[node_id].p, sizeof sample->p);
                    sample->s = cumulative;
                    sample->stroke = (int32_t)nstroke;
                    sample->ordinal = (int32_t)j;
                    sample->confidence =
                        (j == 0 || j + 1 == nselected) ? 0.75 : 1.0;

                    sample_ref[nsample + j].mesh_vertex[0] =
                        nodes[node_id].edge[0];
                    sample_ref[nsample + j].mesh_vertex[1] =
                        nodes[node_id].edge[1];
                    sample_ref[nsample + j].mesh_t = nodes[node_id].mesh_t;
                    sample_ref[nsample + j].cube = nodes[node_id].cube;
                    sample_ref[nsample + j].plane = nodes[node_id].plane;
                    sample_ref[nsample + j].mesh_component =
                        mesh_component[nodes[node_id].edge[0]];
                    sample_ref[nsample + j].raw_u = nodes[node_id].raw_u;
                    sample_ref[nsample + j].axial =
                        ((double)nodes[node_id].plane + 0.5) *
                        o->slice_spacing;
                    sample_ref[nsample + j].radius = 0.0;
                    sample_ref[nsample + j].theta = 0.0;
                    sample_ref[nsample + j].wind = 0.0;
                    sample_ref[nsample + j].u_geo = 0.0;
                    if (sp.enabled) {
                        size_t ck = (size_t)selected_k[j];
                        sample_ref[nsample + j].radius = sr[ck];
                        sample_ref[nsample + j].theta = sth[ck];
                        sample_ref[nsample + j].wind = ac_spiral_wind(
                            &sp, sr[ck], sth[ck]);
                        sample_ref[nsample + j].u_geo = sug[ck];
                        /* Offset-only prior: one row at the stroke head.
                         * A per-sample prior would carry u_geo's wobble
                         * into the shape and fight the monotone cone; the
                         * head row decides only which wind the fragment
                         * sits on. */
                        geo_u[nsample + j] = (trusted && j == 0)
                                           ? sug[ck] : (double)NAN;
                    }
                    /* Unit-speed intrinsic initialization is strictly inside
                     * the beta=1/2 feasible cone and retains only this
                     * fragment's raw gauge. */
                    initial_u[nsample + j] = gauge + cumulative;
                }
                for (size_t j = 0; j < nselected; j++) {
                    size_t ja = j > 0 ? j - 1 : j;
                    size_t jb = j + 1 < nselected ? j + 1 : j;
                    double tangent[3];
                    for (int d = 0; d < 3; d++)
                        tangent[d] =
                            samples[nsample + jb].p[d] -
                            samples[nsample + ja].p[d];
                    if (ac_normalize3(tangent) < 1e-12) {
                        tangent[0] = 0.0;
                        tangent[1] = 1.0;
                        tangent[2] = 0.0;
                    }
                    memcpy(samples[nsample + j].tangent, tangent,
                           sizeof tangent);
                }
                nsample += nselected;
                nstroke++;
                }
                (void)closed;
            }
        }
    }
    if (nstroke == 0 || nsample < 2) return -1;

    size_t nplane = (size_t)((int64_t)plane_max - (int64_t)plane_min + 1);
    if (nplane == 0 || nplane > (size_t)INT32_MAX) return -1;

    size_t ncoarse_edge = nsample - nstroke;
    CoarseEdge *coarse_edge = (CoarseEdge *)ARENA_ALLOC(
        arena, ncoarse_edge * sizeof(CoarseEdge));
    size_t *plane_edge_count = (size_t *)ARENA_CALLOC(
        arena, nplane, sizeof(size_t));
    size_t ne = 0;
    for (size_t s = 0; s < nstroke; s++) {
        const AtlasStripStroke *stroke = &strokes[s];
        int32_t plane = sample_ref[stroke->first].plane;
        for (int32_t j = 0; j + 1 < stroke->count; j++) {
            coarse_edge[ne].a = stroke->first + j;
            coarse_edge[ne].b = stroke->first + j + 1;
            coarse_edge[ne].stroke = (int32_t)s;
            coarse_edge[ne].plane = plane;
            plane_edge_count[(size_t)(plane - plane_min)]++;
            ne++;
        }
    }
    if (ne != ncoarse_edge) return -1;

    size_t *plane_edge_off = (size_t *)ARENA_ALLOC(
        arena, (nplane + 1) * sizeof(size_t));
    plane_edge_off[0] = 0;
    for (size_t p = 0; p < nplane; p++)
        plane_edge_off[p + 1] = plane_edge_off[p] + plane_edge_count[p];
    size_t *plane_cursor = (size_t *)ARENA_ALLOC(
        arena, nplane * sizeof(size_t));
    memcpy(plane_cursor, plane_edge_off, nplane * sizeof(size_t));
    int32_t *plane_edge_ids = (int32_t *)ARENA_ALLOC(
        arena, ncoarse_edge * sizeof(int32_t));
    for (size_t e = 0; e < ncoarse_edge; e++) {
        size_t p = (size_t)(coarse_edge[e].plane - plane_min);
        plane_edge_ids[plane_cursor[p]++] = (int32_t)e;
    }

    float *edge_midpoint = (float *)ARENA_ALLOC(
        arena, ncoarse_edge * 3 * sizeof(float));
    double *plane_max_half = (double *)ARENA_CALLOC(
        arena, nplane, sizeof(double));
    for (size_t p = 0; p < nplane; p++) {
        for (size_t q = plane_edge_off[p]; q < plane_edge_off[p + 1]; q++) {
            const CoarseEdge *edge = &coarse_edge[plane_edge_ids[q]];
            for (int d = 0; d < 3; d++)
                edge_midpoint[q * 3 + (size_t)d] =
                    (float)(0.5 * (samples[edge->a].p[d] +
                                   samples[edge->b].p[d]));
            double half = 0.5 * ac_dist3(samples[edge->a].p,
                                         samples[edge->b].p);
            if (half > plane_max_half[p]) plane_max_half[p] = half;
        }
    }
    KDTree_T *edge_tree = (KDTree_T *)ARENA_CALLOC(
        arena, nplane, sizeof(KDTree_T));
    for (size_t p = 0; p < nplane; p++) {
        size_t count = plane_edge_off[p + 1] - plane_edge_off[p];
        if (count > 0)
            edge_tree[p] = KDTree_new(
                arena, &edge_midpoint[plane_edge_off[p] * 3], count);
    }

    if (nsample > SIZE_MAX / (size_t)(o->max_candidates + 1))
        return -1;
    size_t member_capacity =
        nsample * (size_t)(o->max_candidates + 1);
    AtlasStripCrossSection *cross_sections =
        (AtlasStripCrossSection *)ARENA_ALLOC(
            arena, nsample * sizeof(AtlasStripCrossSection));
    AtlasStripMember *members = (AtlasStripMember *)ARENA_ALLOC(
        arena, member_capacity * sizeof(AtlasStripMember));
    size_t ncross = 0, nmember = 0;
    double cos_match = cos(o->match_angle_deg * M_PI / 180.0);
    const double match2 = o->match_radius * o->match_radius;
    int32_t ball[AC_BALL_CAP];
    CandidateHit candidate[AC_MAX_CANDIDATES];

    for (size_t source = 0; source < nsample; source++) {
        int32_t source_plane = sample_ref[source].plane;
        int ncandidate = 0;
        int used_stride = 0;
        for (int stride = 1;
             stride <= o->max_slice_stride && ncandidate == 0;
             stride++) {
            int64_t target_plane64 = (int64_t)source_plane + stride;
            if (target_plane64 > plane_max) break;
            size_t pslot =
                (size_t)((int32_t)target_plane64 - plane_min);
            if (pslot >= nplane || edge_tree[pslot] == NULL) continue;

            double axial_gap = (double)stride * o->slice_spacing;
            double reach = o->match_radius + plane_max_half[pslot];
            double query_r2 = axial_gap * axial_gap + reach * reach + 1e-6;
            float query[3] = {
                (float)samples[source].p[0],
                (float)samples[source].p[1],
                (float)samples[source].p[2]
            };
            size_t nh = KDTree_ball_query(
                edge_tree[pslot], query, (float)query_r2,
                ball, AC_BALL_CAP);
            if (nh == AC_BALL_CAP) out->stats.truncated_ball_queries++;
            for (size_t h = 0; h < nh; h++) {
                size_t local = (size_t)ball[h];
                size_t ordered = plane_edge_off[pslot] + local;
                if (ordered >= plane_edge_off[pslot + 1]) continue;
                int32_t edge_id = plane_edge_ids[ordered];
                const CoarseEdge *edge = &coarse_edge[edge_id];
                double ab[3], aq[3];
                for (int d = 0; d < 3; d++) {
                    ab[d] = samples[edge->b].p[d] -
                            samples[edge->a].p[d];
                    aq[d] = samples[source].p[d] -
                            samples[edge->a].p[d];
                }
                double den = ac_dot3(ab, ab);
                if (den <= 1e-12) continue;
                double t = ac_dot3(aq, ab) / den;
                if (t < 0.0) t = 0.0;
                if (t > 1.0) t = 1.0;
                double delta[3], tangent[3], cp[3];
                for (int d = 0; d < 3; d++) {
                    cp[d] = samples[edge->a].p[d] + t * ab[d];
                    delta[d] = cp[d] - samples[source].p[d];
                    tangent[d] = ab[d];
                }
                double axial_delta = ac_dot3(delta, axis);
                double lateral2 = ac_dot3(delta, delta) -
                                  axial_delta * axial_delta;
                if (lateral2 < 0.0) lateral2 = 0.0;
                if (lateral2 > match2) continue;
                if (ac_normalize3(tangent) < 1e-12) continue;
                double tangent_dot =
                    ac_dot3(samples[source].tangent, tangent);
                if (tangent_dot < cos_match) continue;
                /*
                 * Winding gate: a cross-section asserts "same isovalue",
                 * which is only true of same-wrap pairs.  atlas_strip.h
                 * documents the sub-half-wrap contract; this enforces it
                 * where the relation is born, with the pairwise delta.
                 */
                if (sp.enabled) {
                    double cr = 0.0, cth = 0.0;
                    ac_spiral_polar(&sp, cp, &cr, &cth);
                    double dw = ac_spiral_dwind(
                        &sp, sample_ref[source].radius,
                        sample_ref[source].theta, cr, cth);
                    out->stats.member_wind_hist[ac_wind_bin(dw)]++;
                    if (ac_wind_fusing(dw))
                        out->stats.member_wind_fusing++;
                    if (o->wind_gate && fabs(dw) > sp.tol) {
                        out->stats.member_wind_rejected++;
                        continue;
                    }
                }
                double sigma = 0.5 * o->match_radius;
                double score = exp(-0.5 * lateral2 / (sigma * sigma)) *
                               tangent_dot * tangent_dot *
                               tangent_dot * tangent_dot /
                               (double)stride;
                CandidateHit hit;
                hit.edge = edge_id;
                hit.stroke = edge->stroke;
                hit.t = t;
                hit.lateral2 = lateral2;
                hit.score = score;
                ac_keep_candidate(candidate, &ncandidate,
                                  o->max_candidates, &hit);
            }
            if (ncandidate > 0) used_stride = stride;
        }
        if (ncandidate == 0) continue;
        qsort(candidate, (size_t)ncandidate, sizeof(candidate[0]),
              ac_compare_candidate);

        if (nmember + (size_t)ncandidate + 1 > member_capacity)
            return -1;
        AtlasStripCrossSection *cs = &cross_sections[ncross];
        memset(cs, 0, sizeof *cs);
        cs->first = nmember;
        cs->count = ncandidate + 1;
        cs->weight = o->slice_spacing * (double)used_stride;
        cs->id = (int32_t)ncross;

        ac_member_direct(&members[nmember++], samples, strokes,
                         (int32_t)source, 1.0, (int32_t)source);
        double tangent_sum[3] = {
            samples[source].tangent[0],
            samples[source].tangent[1],
            samples[source].tangent[2]
        };
        for (int j = 0; j < ncandidate; j++) {
            const CoarseEdge *edge = &coarse_edge[candidate[j].edge];
            ac_member_edge(&members[nmember++], samples, edge,
                           candidate[j].t, candidate[j].score,
                           candidate[j].edge);
            for (int d = 0; d < 3; d++)
                tangent_sum[d] += candidate[j].score *
                                  members[nmember - 1].tangent[d];
        }
        if (ac_normalize3(tangent_sum) < 1e-12)
            memcpy(tangent_sum, samples[source].tangent,
                   sizeof tangent_sum);
        memcpy(cs->tangent, tangent_sum, sizeof cs->tangent);
        if (ncandidate > 1) out->stats.ambiguous_cross_sections++;
        out->stats.matched_source_samples++;
        ncross++;
    }
    if (ncross == 0 || nmember == 0) return -1;

    /* Same-plane endpoint continuation candidates. */
    size_t *plane_stroke_count = (size_t *)ARENA_CALLOC(
        arena, nplane, sizeof(size_t));
    for (size_t s = 0; s < nstroke; s++) {
        int32_t plane = sample_ref[strokes[s].first].plane;
        plane_stroke_count[(size_t)(plane - plane_min)]++;
    }
    size_t *plane_stroke_off = (size_t *)ARENA_ALLOC(
        arena, (nplane + 1) * sizeof(size_t));
    plane_stroke_off[0] = 0;
    for (size_t p = 0; p < nplane; p++)
        plane_stroke_off[p + 1] =
            plane_stroke_off[p] + plane_stroke_count[p];
    memcpy(plane_cursor, plane_stroke_off, nplane * sizeof(size_t));
    int32_t *plane_stroke_ids = (int32_t *)ARENA_ALLOC(
        arena, nstroke * sizeof(int32_t));
    for (size_t s = 0; s < nstroke; s++) {
        int32_t plane = sample_ref[strokes[s].first].plane;
        size_t p = (size_t)(plane - plane_min);
        plane_stroke_ids[plane_cursor[p]++] = (int32_t)s;
    }
    float *stroke_start = (float *)ARENA_ALLOC(
        arena, nstroke * 3 * sizeof(float));
    for (size_t p = 0; p < nplane; p++)
        for (size_t q = plane_stroke_off[p];
             q < plane_stroke_off[p + 1]; q++) {
            int32_t s = plane_stroke_ids[q];
            int32_t first = strokes[s].first;
            for (int d = 0; d < 3; d++)
                stroke_start[q * 3 + (size_t)d] =
                    (float)samples[first].p[d];
        }
    KDTree_T *start_tree = (KDTree_T *)ARENA_CALLOC(
        arena, nplane, sizeof(KDTree_T));
    for (size_t p = 0; p < nplane; p++) {
        size_t count = plane_stroke_off[p + 1] - plane_stroke_off[p];
        if (count > 0)
            start_tree[p] = KDTree_new(
                arena, &stroke_start[plane_stroke_off[p] * 3], count);
    }
    int32_t *best_start_for_end = (int32_t *)ARENA_ALLOC(
        arena, nstroke * sizeof(int32_t));
    int32_t *best_end_for_start = (int32_t *)ARENA_ALLOC(
        arena, nstroke * sizeof(int32_t));
    double *best_end_cost = (double *)ARENA_ALLOC(
        arena, nstroke * sizeof(double));
    double *best_start_cost = (double *)ARENA_ALLOC(
        arena, nstroke * sizeof(double));
    for (size_t s = 0; s < nstroke; s++) {
        best_start_for_end[s] = -1;
        best_end_for_start[s] = -1;
        best_end_cost[s] = DBL_MAX;
        best_start_cost[s] = DBL_MAX;
    }
    double cos_cont =
        cos(o->continuation_angle_deg * M_PI / 180.0);
    double cont2 = o->continuation_radius * o->continuation_radius;
    for (size_t s = 0; s < nstroke; s++) {
        int32_t last = strokes[s].first + strokes[s].count - 1;
        int32_t plane = sample_ref[last].plane;
        size_t pslot = (size_t)(plane - plane_min);
        if (start_tree[pslot] == NULL) continue;
        float query[3] = {
            (float)samples[last].p[0],
            (float)samples[last].p[1],
            (float)samples[last].p[2]
        };
        size_t nh = KDTree_ball_query(
            start_tree[pslot], query, (float)cont2,
            ball, AC_BALL_CAP);
        if (nh == AC_BALL_CAP) out->stats.truncated_ball_queries++;
        for (size_t h = 0; h < nh; h++) {
            size_t ordered = plane_stroke_off[pslot] + (size_t)ball[h];
            if (ordered >= plane_stroke_off[pslot + 1]) continue;
            int32_t target_stroke = plane_stroke_ids[ordered];
            if ((size_t)target_stroke == s) continue;
            int32_t first = strokes[target_stroke].first;
            double tangent_dot =
                ac_dot3(samples[last].tangent, samples[first].tangent);
            if (tangent_dot < cos_cont) continue;
            double gap[3], tangent[3];
            for (int d = 0; d < 3; d++) {
                gap[d] = samples[first].p[d] - samples[last].p[d];
                tangent[d] = samples[last].tangent[d] +
                             samples[first].tangent[d];
            }
            double distance2 = ac_dot3(gap, gap);
            if (distance2 > cont2) continue;
            if (ac_normalize3(tangent) < 1e-12) continue;
            double forward = ac_dot3(tangent, gap);
            double lateral2 = distance2 - forward * forward;
            if (lateral2 < 0.0) lateral2 = 0.0;
            if (forward < -0.25 * o->continuation_radius ||
                lateral2 > 0.75 * 0.75 * cont2)
                continue;
            /* Winding gate: an end-to-start join across wraps would carry
             * a gap-sized target where a circumference is owed.  The
             * tolerance here is TIGHTER than the cross-section gate: a
             * genuine continuation along one wrap has dwind ~ 0.0x by the
             * spiral property, while a tight-wound core neighbour at
             * sub-half-pitch spacing shows ~0.3 -- inside the relation
             * gate's blind zone but far outside a continuation's. */
            if (sp.enabled) {
                double cont_tol = sp.tol < 0.15 ? sp.tol : 0.15;
                double dw = ac_spiral_dwind(
                    &sp, sample_ref[last].radius, sample_ref[last].theta,
                    sample_ref[first].radius, sample_ref[first].theta);
                out->stats.continuation_wind_hist[ac_wind_bin(dw)]++;
                if (ac_wind_fusing(dw))
                    out->stats.continuation_wind_fusing++;
                if (o->wind_gate && fabs(dw) > cont_tol) {
                    out->stats.continuation_wind_rejected++;
                    continue;
                }
            }
            double cost = sqrt(distance2) /
                          (0.25 + 0.75 * tangent_dot);
            if (cost < best_end_cost[s]) {
                best_end_cost[s] = cost;
                best_start_for_end[s] = target_stroke;
            }
            if (cost < best_start_cost[target_stroke]) {
                best_start_cost[target_stroke] = cost;
                best_end_for_start[target_stroke] = (int32_t)s;
            }
        }
    }

    AtlasStripContinuation *continuations =
        (AtlasStripContinuation *)ARENA_ALLOC(
            arena, nstroke * sizeof(AtlasStripContinuation));
    size_t ncontinuation = 0;
    for (size_t s = 0; s < nstroke; s++) {
        int32_t target_stroke = best_start_for_end[s];
        if (target_stroke < 0 ||
            best_end_for_start[target_stroke] != (int32_t)s)
            continue;
        int32_t a = strokes[s].first + strokes[s].count - 1;
        int32_t b = strokes[target_stroke].first;
        double tangent[3], delta[3];
        for (int d = 0; d < 3; d++) {
            tangent[d] = samples[a].tangent[d] + samples[b].tangent[d];
            delta[d] = samples[a].p[d] - samples[b].p[d];
        }
        if (ac_normalize3(tangent) < 1e-12) continue;
        double distance2 = ac_dot3(delta, delta);
        double tangent_dot =
            ac_dot3(samples[a].tangent, samples[b].tangent);
        double sigma = 0.5 * o->continuation_radius;
        double score = exp(-0.5 * distance2 / (sigma * sigma)) *
                       tangent_dot * tangent_dot;
        AtlasStripContinuation *cn = &continuations[ncontinuation];
        cn->a = a;
        cn->b = b;
        cn->target = ac_dot3(tangent, delta);
        cn->weight = 1.0;
        cn->membership = score > 1e-6 ? score : 1e-6;
        cn->observation = (int32_t)ncontinuation;
        ncontinuation++;
    }

    UnionFind support = UF_new(arena, (int32_t)nsample);
    for (size_t s = 0; s < nstroke; s++) {
        int32_t first = strokes[s].first;
        for (int32_t j = 1; j < strokes[s].count; j++)
            uf_union(&support, first, first + j);
    }
    for (size_t i = 0; i < ncontinuation; i++)
        uf_union(&support, continuations[i].a, continuations[i].b);
    for (size_t c = 0; c < ncross; c++) {
        const AtlasStripCrossSection *cs = &cross_sections[c];
        int32_t root_sample = members[cs->first].value0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m = &members[cs->first + (size_t)j];
            uf_union(&support, root_sample, m->value0);
            if (m->value1 >= 0)
                uf_union(&support, root_sample, m->value1);
        }
    }

    size_t *root_size = (size_t *)ARENA_CALLOC(
        arena, nsample, sizeof(size_t));
    int32_t *root_best = (int32_t *)ARENA_ALLOC(
        arena, nsample * sizeof(int32_t));
    int32_t *root_component = (int32_t *)ARENA_ALLOC(
        arena, nsample * sizeof(int32_t));
    for (size_t i = 0; i < nsample; i++) {
        root_best[i] = -1;
        root_component[i] = -1;
        root_size[uf_find(&support, (int32_t)i)]++;
    }
    /* Anchor stroke preference: a PRIOR-TRUSTED stroke (finite geo_u at its
     * head) beats a longer untrusted one -- an exact anchor pinned a wrap
     * off is the one mistake the soft prior can never out-vote. */
    for (size_t s = 0; s < nstroke; s++) {
        int32_t root = uf_find(&support, strokes[s].first);
        int32_t best = root_best[root];
        int st = geo_u != NULL && isfinite(geo_u[strokes[s].first]);
        int bt = best >= 0 && geo_u != NULL &&
                 isfinite(geo_u[strokes[best].first]);
        if (best < 0 || (st && !bt) ||
            (st == bt && strokes[s].count > strokes[best].count))
            root_best[root] = (int32_t)s;
    }

    MonotoneQpAnchor *anchors = (MonotoneQpAnchor *)ARENA_ALLOC(
        arena, nstroke * sizeof(MonotoneQpAnchor));
    size_t nanchor = 0, largest = 0;
    for (size_t i = 0; i < nsample; i++) {
        if (uf_find(&support, (int32_t)i) != (int32_t)i) continue;
        int32_t best = root_best[i];
        if (best < 0) continue;
        int32_t sample = strokes[best].first;
        root_component[i] = (int32_t)nanchor;
        anchors[nanchor].var = sample;
        anchors[nanchor].value = initial_u[sample];
        anchors[nanchor].component = (int32_t)nanchor;
        if (root_size[i] > largest) largest = root_size[i];
        nanchor++;
    }
    for (size_t s = 0; s < nstroke; s++) {
        int32_t root = uf_find(&support, strokes[s].first);
        strokes[s].component = root_component[root];
        strokes[s].resolved = root_size[root] > 1;
    }

    out->samples = samples;
    out->strokes = strokes;
    out->members = members;
    out->cross_sections = cross_sections;
    out->continuations = continuations;
    out->anchors = anchors;
    out->sample_ref = sample_ref;
    out->vertex_mesh_component = mesh_component;
    out->initial_u = initial_u;
    out->geo_u = geo_u;
    out->member_prior = (double *)ARENA_ALLOC(
        arena, nmember * sizeof(double));
    out->continuation_prior = (double *)ARENA_ALLOC(
        arena, (ncontinuation ? ncontinuation : 1) * sizeof(double));
    out->member_state = (uint8_t *)ARENA_CALLOC(
        arena, nmember, sizeof(uint8_t));
    out->continuation_state = (uint8_t *)ARENA_CALLOC(
        arena, (ncontinuation ? ncontinuation : 1), sizeof(uint8_t));
    out->member_bundle = (int32_t *)ARENA_ALLOC(
        arena, nmember * sizeof(int32_t));
    out->continuation_bundle = (int32_t *)ARENA_ALLOC(
        arena, (ncontinuation ? ncontinuation : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nmember; i++)
        out->member_prior[i] = members[i].membership;
    for (size_t i = 0; i < ncontinuation; i++)
        out->continuation_prior[i] = continuations[i].membership;
    for (size_t i = 0; i < nmember; i++) out->member_bundle[i] = -1;
    for (size_t i = 0; i < ncontinuation; i++)
        out->continuation_bundle[i] = -1;
    out->plane_min = plane_min;
    out->plane_max = plane_max;
    out->mesh_components = nmesh_component;

    out->problem.samples = samples;
    out->problem.nsamples = nsample;
    out->problem.strokes = strokes;
    out->problem.nstrokes = nstroke;
    out->problem.members = members;
    out->problem.nmembers = nmember;
    out->problem.cross_sections = cross_sections;
    out->problem.ncross_sections = ncross;
    out->problem.continuations = continuations;
    out->problem.ncontinuations = ncontinuation;
    out->problem.anchors = anchors;
    out->problem.nanchors = nanchor;
    out->problem.prior_u = geo_u;

    out->stats.intersection_nodes = nnode;
    out->stats.intersection_segments = nsegment;
    out->stats.strokes = nstroke;
    out->stats.samples = nsample;
    out->stats.continuations = ncontinuation;
    out->stats.cross_sections = ncross;
    out->stats.members = nmember;
    out->stats.support_components = nanchor;
    out->stats.largest_component_samples = largest;
    out->stats.source_match_fraction =
        nsample > 0 ? (double)out->stats.matched_source_samples /
                      (double)nsample : 0.0;
    /* Pass one is deliberately the exact original all-admitted graph.  Its
     * overlaps are evidence for the later classifier, not something to erase
     * before the global map exists. */
    return AtlasCandidates_select_baseline(arena, 1e-4, out);
}

static int ac_compare_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa;
    double b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int32_t ac_member_mesh_component(const AtlasCandidateSet *set,
                                        const AtlasStripMember *member)
{
    int32_t component = set->sample_ref[member->value0].mesh_component;
    if (member->value1 >= 0 &&
        set->sample_ref[member->value1].mesh_component != component)
        return -1;
    return component;
}

static int ac_rebuild_support(Arena_T arena, AtlasCandidateSet *set,
                              double support_floor)
{
    AtlasStripProblem *p = &set->problem;
    if (p->nsamples == 0 || p->nstrokes == 0 || set->anchors == NULL)
        return -1;
    UnionFind support = UF_new(arena, (int32_t)p->nsamples);
    for (size_t s = 0; s < p->nstrokes; s++) {
        int32_t first = set->strokes[s].first;
        for (int32_t j = 1; j < set->strokes[s].count; j++)
            uf_union(&support, first, first + j);
    }
    for (size_t i = 0; i < p->ncontinuations; i++) {
        const AtlasStripContinuation *cn = &set->continuations[i];
        int active = set->continuation_state != NULL
                   ? set->continuation_state[i] != ATLAS_CANDIDATE_INACTIVE
                   : cn->membership * cn->weight > support_floor;
        if (active)
            uf_union(&support, cn->a, cn->b);
    }
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        int32_t root_sample = set->members[cs->first].value0;
        for (int32_t j = 1; j < cs->count; j++) {
            const AtlasStripMember *m =
                &set->members[cs->first + (size_t)j];
            size_t mi = cs->first + (size_t)j;
            int active = set->member_state != NULL
                       ? set->member_state[mi] != ATLAS_CANDIDATE_INACTIVE
                       : m->membership * m->base_weight > support_floor;
            if (!active) continue;
            uf_union(&support, root_sample, m->value0);
            if (m->value1 >= 0)
                uf_union(&support, root_sample, m->value1);
        }
    }

    size_t *root_size = (size_t *)ARENA_CALLOC(
        arena, p->nsamples, sizeof(size_t));
    int32_t *root_best = (int32_t *)ARENA_ALLOC(
        arena, p->nsamples * sizeof(int32_t));
    int32_t *root_component = (int32_t *)ARENA_ALLOC(
        arena, p->nsamples * sizeof(int32_t));
    for (size_t i = 0; i < p->nsamples; i++) {
        root_best[i] = -1;
        root_component[i] = -1;
        root_size[uf_find(&support, (int32_t)i)]++;
    }
    for (size_t s = 0; s < p->nstrokes; s++) {
        int32_t root = uf_find(&support, set->strokes[s].first);
        int32_t best = root_best[root];
        int st = set->geo_u != NULL &&
                 isfinite(set->geo_u[set->strokes[s].first]);
        int bt = best >= 0 && set->geo_u != NULL &&
                 isfinite(set->geo_u[set->strokes[best].first]);
        if (best < 0 || (st && !bt) ||
            (st == bt && set->strokes[s].count > set->strokes[best].count))
            root_best[root] = (int32_t)s;
    }

    size_t nanchor = 0, largest = 0;
    for (size_t i = 0; i < p->nsamples; i++) {
        if (uf_find(&support, (int32_t)i) != (int32_t)i) continue;
        int32_t best = root_best[i];
        if (best < 0) continue;
        int32_t sample = set->strokes[best].first;
        root_component[i] = (int32_t)nanchor;
        set->anchors[nanchor].var = sample;
        set->anchors[nanchor].value = set->initial_u[sample];
        set->anchors[nanchor].component = (int32_t)nanchor;
        if (root_size[i] > largest) largest = root_size[i];
        nanchor++;
    }
    for (size_t s = 0; s < p->nstrokes; s++) {
        int32_t root = uf_find(&support, set->strokes[s].first);
        set->strokes[s].component = root_component[root];
        set->strokes[s].resolved = root_size[root] > 1;
    }
    p->anchors = set->anchors;
    p->nanchors = nanchor;
    set->stats.support_components = nanchor;
    set->stats.largest_component_samples = largest;
    return 0;
}

int AtlasCandidates_reanchor(Arena_T arena, AtlasCandidateSet *set,
                             double support_floor)
{
    if (arena == NULL || set == NULL) return -1;
    return ac_rebuild_support(arena, set, support_floor);
}

int AtlasCandidates_select_baseline(Arena_T arena,
                                    double support_floor,
                                    AtlasCandidateSet *set)
{
    if (arena == NULL || set == NULL || !isfinite(support_floor) ||
        support_floor < 0.0 || set->member_prior == NULL ||
        set->continuation_prior == NULL || set->member_state == NULL ||
        set->continuation_state == NULL)
        return -1;
    AtlasStripProblem *p = &set->problem;
    if (p->nsamples == 0 || set->members == NULL ||
        set->continuations == NULL)
        return -1;

    memset(&set->selection, 0, sizeof set->selection);
    set->selection.mesh_components = set->mesh_components;
    set->bundles = NULL;
    set->nbundles = 0;
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        for (int32_t j = 0; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            set->members[mi].membership = set->member_prior[mi];
            set->member_state[mi] = j == 0
                ? ATLAS_CANDIDATE_SOURCE : ATLAS_CANDIDATE_BASELINE;
            set->member_bundle[mi] = -1;
            if (j > 0) set->selection.candidate_links++;
        }
    }
    for (size_t i = 0; i < p->ncontinuations; i++) {
        set->continuations[i].membership = set->continuation_prior[i];
        set->continuation_state[i] = ATLAS_CANDIDATE_BASELINE;
        set->continuation_bundle[i] = -1;
    }
    return ac_rebuild_support(arena, set, support_floor);
}

int AtlasCandidates_select_conservative(Arena_T arena,
                                        double support_floor,
                                        AtlasCandidateSet *set)
{
    if (arena == NULL || set == NULL || !isfinite(support_floor) ||
        support_floor < 0.0 || set->member_prior == NULL ||
        set->continuation_prior == NULL || set->member_state == NULL ||
        set->continuation_state == NULL)
        return -1;
    AtlasStripProblem *p = &set->problem;
    if (p->nsamples == 0 || set->sample_ref == NULL ||
        set->members == NULL || set->continuations == NULL)
        return -1;

    memset(&set->selection, 0, sizeof set->selection);
    set->selection.mesh_components = set->mesh_components;
    set->bundles = NULL;
    set->nbundles = 0;
    for (size_t i = 0; i < p->nmembers; i++) {
        set->members[i].membership = 0.0;
        set->member_state[i] = ATLAS_CANDIDATE_INACTIVE;
        set->member_bundle[i] = -1;
    }
    for (size_t i = 0; i < p->ncontinuations; i++) {
        set->continuations[i].membership = 0.0;
        set->continuation_state[i] = ATLAS_CANDIDATE_INACTIVE;
        set->continuation_bundle[i] = -1;
    }

    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        size_t source_index = cs->first;
        AtlasStripMember *source = &set->members[source_index];
        int32_t source_component = ac_member_mesh_component(set, source);
        source->membership = set->member_prior[source_index];
        set->member_state[source_index] = ATLAS_CANDIDATE_SOURCE;

        size_t best = SIZE_MAX;
        double best_prior = -DBL_MAX;
        size_t topology_count = 0;
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            int32_t target_component =
                ac_member_mesh_component(set, &set->members[mi]);
            set->selection.candidate_links++;
            if (target_component == source_component) {
                topology_count++;
                set->selection.topology_candidate_links++;
                if (set->member_prior[mi] > best_prior) {
                    best_prior = set->member_prior[mi];
                    best = mi;
                }
            } else {
                set->selection.inactive_cross_component_links++;
            }
        }
        if (topology_count > 1)
            set->selection.topology_ambiguous_cross_sections++;
        if (best != SIZE_MAX) {
            set->members[best].membership = set->member_prior[best];
            set->member_state[best] = ATLAS_CANDIDATE_TOPOLOGY;
            set->selection.topology_selected_cross_sections++;
        } else {
            set->selection.unmatched_cross_sections++;
        }
    }

    for (size_t i = 0; i < p->ncontinuations; i++) {
        AtlasStripContinuation *cn = &set->continuations[i];
        int32_t ca = set->sample_ref[cn->a].mesh_component;
        int32_t cb = set->sample_ref[cn->b].mesh_component;
        if (ca == cb) {
            cn->membership = set->continuation_prior[i];
            set->continuation_state[i] = ATLAS_CANDIDATE_TOPOLOGY;
            set->selection.topology_selected_continuations++;
        } else {
            set->selection.inactive_cross_component_continuations++;
        }
    }
    return ac_rebuild_support(arena, set, support_floor);
}

typedef struct {
    int32_t component0;
    int32_t component1;
    int32_t plane;
    size_t member;
    size_t continuation;
    double initial_residual;
    double residual;
    double prior;
    double arc0;
    double arc1;
    double gap;
    int top_prior;
    int has_topology_alternative;
} AcRefinementLink;

static int ac_compare_refinement_link(const void *pa, const void *pb)
{
    const AcRefinementLink *a = (const AcRefinementLink *)pa;
    const AcRefinementLink *b = (const AcRefinementLink *)pb;
    if (a->component0 != b->component0)
        return a->component0 < b->component0 ? -1 : 1;
    if (a->component1 != b->component1)
        return a->component1 < b->component1 ? -1 : 1;
    if (a->plane != b->plane) return a->plane < b->plane ? -1 : 1;
    if (a->member != b->member) return a->member < b->member ? -1 : 1;
    return a->continuation < b->continuation ? -1 :
           (a->continuation > b->continuation ? 1 : 0);
}

static double ac_median(double *value, size_t n)
{
    if (n == 0) return NAN;
    qsort(value, n, sizeof(double), ac_compare_double);
    if (n & 1u) return value[n / 2];
    return 0.5 * (value[n / 2 - 1] + value[n / 2]);
}

static double ac_member_arc(const AtlasCandidateSet *set,
                            const AtlasStripMember *m)
{
    double a = set->samples[m->value0].s;
    if (m->value1 < 0) return a;
    return (1.0 - m->value_t) * a +
           m->value_t * set->samples[m->value1].s;
}

static double ac_member_axial(const AtlasCandidateSet *set,
                              const AtlasStripMember *m)
{
    double a = set->sample_ref[m->value0].axial;
    if (m->value1 < 0) return a;
    return (1.0 - m->value_t) * a +
           m->value_t * set->sample_ref[m->value1].axial;
}

static double ac_cross_lateral_gap(const AtlasCandidateSet *set,
                                   const AtlasStripMember *a,
                                   const AtlasStripMember *b)
{
    double d2 = 0.0;
    for (int d = 0; d < 3; d++) {
        double q = b->p[d] - a->p[d];
        d2 += q * q;
    }
    double da = ac_member_axial(set, b) - ac_member_axial(set, a);
    double lateral2 = d2 - da * da;
    if (lateral2 < 0.0) lateral2 = 0.0;
    return sqrt(lateral2);
}

static int ac_validate_refine_options(const AtlasCandidateRefineOptions *o)
{
    return o != NULL && isfinite(o->support_floor) &&
           o->support_floor >= 0.0 &&
           isfinite(o->gross_initial_residual) &&
           o->gross_initial_residual > 0.0 &&
           isfinite(o->initial_coherence_mad_limit) &&
           o->initial_coherence_mad_limit > 0.0 &&
           isfinite(o->overlap_residual_limit) &&
           o->overlap_residual_limit > 0.0 &&
           isfinite(o->overlap_coherence_mad_limit) &&
           o->overlap_coherence_mad_limit > 0.0 &&
           isfinite(o->persistent_contact_fraction) &&
           o->persistent_contact_fraction > 0.0 &&
           o->persistent_contact_fraction <= 1.0 &&
           isfinite(o->minimum_topology_alternative_fraction) &&
           o->minimum_topology_alternative_fraction >= 0.0 &&
           o->minimum_topology_alternative_fraction <= 1.0 &&
           isfinite(o->ambiguity_ratio_limit) &&
           o->ambiguity_ratio_limit >= 0.0 &&
           o->ambiguity_ratio_limit < 1.0 &&
           isfinite(o->minimum_prior) && o->minimum_prior >= 0.0 &&
           isfinite(o->minimum_top_prior_fraction) &&
           o->minimum_top_prior_fraction >= 0.0 &&
           o->minimum_top_prior_fraction <= 1.0 &&
           isfinite(o->minimum_arc_correlation) &&
           o->minimum_arc_correlation >= -1.0 &&
           o->minimum_arc_correlation <= 1.0 &&
           o->minimum_links >= 1 && o->minimum_planes >= 1 &&
           o->front_minimum_links >= 1 &&
           o->front_minimum_planes >= 1 &&
           o->front_maximum_extension >= 0 &&
           isfinite(o->front_close_gap) && o->front_close_gap > 0.0 &&
           isfinite(o->front_minimum_contact_overlap) &&
           o->front_minimum_contact_overlap >= 0.0 &&
           o->front_minimum_contact_overlap <= 1.0;
}

static void ac_best_bundle_update(int32_t component, int32_t bundle,
                                  double score, int32_t *best_bundle,
                                  double *best_score, double *second_score)
{
    if (score > best_score[component]) {
        second_score[component] = best_score[component];
        best_score[component] = score;
        best_bundle[component] = bundle;
    } else if (score > second_score[component]) {
        second_score[component] = score;
    }
}

int AtlasCandidates_select_refined(Arena_T arena,
                                   const double *frozen_u,
                                   const AtlasCandidateRefineOptions *input,
                                   AtlasCandidateSet *set)
{
    AtlasCandidateRefineOptions defaults;
    AtlasCandidateRefineOptions_default(&defaults);
    const AtlasCandidateRefineOptions *o = input != NULL ? input : &defaults;
    if (arena == NULL || frozen_u == NULL || set == NULL ||
        !ac_validate_refine_options(o))
        return -1;
    AtlasStripProblem *p = &set->problem;
    for (size_t i = 0; i < p->nsamples; i++)
        if (!isfinite(frozen_u[i])) return -1;

    /* Classification starts from the exact HEAD assignment graph.  The
     * frozen solution is the observed overlap map and is never allowed to
     * validate a correspondence by itself. */
    if (AtlasCandidates_select_baseline(
            arena, o->support_floor, set) != 0)
        return -1;

    size_t capacity = p->nmembers - p->ncross_sections +
                      p->ncontinuations;
    AcRefinementLink *link = (AcRefinementLink *)ARENA_ALLOC(
        arena, (capacity ? capacity : 1) * sizeof(AcRefinementLink));
    size_t nlink = 0;

    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        const AtlasStripMember *source = &set->members[cs->first];
        int32_t source_component = ac_member_mesh_component(set, source);
        double source_u = AtlasStrip_member_value(source, frozen_u);
        double source_initial =
            AtlasStrip_member_value(source, set->initial_u);
        double source_arc = ac_member_arc(set, source);
        int32_t plane = set->sample_ref[source->value0].plane;
        int has_topology_alternative = 0;
        for (int32_t j = 1; j < cs->count; j++) {
            const AtlasStripMember *candidate =
                &set->members[cs->first + (size_t)j];
            if (ac_member_mesh_component(set, candidate) == source_component) {
                has_topology_alternative = 1;
                break;
            }
        }
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            const AtlasStripMember *target = &set->members[mi];
            int32_t target_component = ac_member_mesh_component(set, target);
            if (target_component < 0 || target_component == source_component)
                continue;
            AcRefinementLink *q = &link[nlink++];
            memset(q, 0, sizeof *q);
            q->component0 = source_component < target_component
                          ? source_component : target_component;
            q->component1 = source_component < target_component
                          ? target_component : source_component;
            q->plane = plane;
            q->member = mi;
            q->continuation = SIZE_MAX;
            double residual = AtlasStrip_member_value(target, frozen_u) -
                              source_u;
            double initial_residual =
                AtlasStrip_member_value(target, set->initial_u) -
                source_initial;
            q->residual = source_component == q->component0
                        ? residual : -residual;
            q->initial_residual = source_component == q->component0
                                ? initial_residual : -initial_residual;
            q->prior = set->member_prior[mi];
            double target_arc = ac_member_arc(set, target);
            if (source_component == q->component0) {
                q->arc0 = source_arc;
                q->arc1 = target_arc;
            } else {
                q->arc0 = target_arc;
                q->arc1 = source_arc;
            }
            q->gap = ac_cross_lateral_gap(set, source, target);
            q->top_prior = j == 1;
            q->has_topology_alternative = has_topology_alternative;
        }
    }

    for (size_t i = 0; i < p->ncontinuations; i++) {
        const AtlasStripContinuation *cn = &set->continuations[i];
        int32_t ca = set->sample_ref[cn->a].mesh_component;
        int32_t cb = set->sample_ref[cn->b].mesh_component;
        if (ca < 0 || cb < 0 || ca == cb) continue;
        AcRefinementLink *q = &link[nlink++];
        memset(q, 0, sizeof *q);
        q->component0 = ca < cb ? ca : cb;
        q->component1 = ca < cb ? cb : ca;
        q->plane = set->sample_ref[cn->a].plane;
        q->member = SIZE_MAX;
        q->continuation = i;
        double row = frozen_u[cn->a] - frozen_u[cn->b] - cn->target;
        q->residual = ca == q->component0 ? -row : row;
        double initial_row = set->initial_u[cn->a] -
                             set->initial_u[cn->b] - cn->target;
        q->initial_residual = ca == q->component0
                            ? -initial_row : initial_row;
        q->prior = set->continuation_prior[i];
        q->arc0 = q->arc1 = NAN;
        q->gap = ac_dist3(set->samples[cn->a].p,
                          set->samples[cn->b].p);
        q->top_prior = 1;
    }
    if (nlink == 0) {
        set->selection.refinement_remaining_active_cross_sections =
            p->ncross_sections;
        return 0;
    }
    qsort(link, nlink, sizeof(AcRefinementLink),
          ac_compare_refinement_link);

    int32_t *component_plane_min = (int32_t *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(int32_t));
    int32_t *component_plane_max = (int32_t *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(int32_t));
    for (size_t k = 0; k < set->mesh_components; k++) {
        component_plane_min[k] = INT32_MAX;
        component_plane_max[k] = INT32_MIN;
    }
    for (size_t i = 0; i < p->nsamples; i++) {
        int32_t component = set->sample_ref[i].mesh_component;
        int32_t plane = set->sample_ref[i].plane;
        if (plane < component_plane_min[component])
            component_plane_min[component] = plane;
        if (plane > component_plane_max[component])
            component_plane_max[component] = plane;
    }

    set->bundles = (AtlasCandidateBundle *)ARENA_CALLOC(
        arena, nlink, sizeof(AtlasCandidateBundle));
    double *scratch = (double *)ARENA_ALLOC(
        arena, nlink * sizeof(double));
    size_t first = 0;
    while (first < nlink) {
        size_t last = first + 1;
        while (last < nlink &&
               link[last].component0 == link[first].component0 &&
               link[last].component1 == link[first].component1)
            last++;
        size_t count = last - first;
        AtlasCandidateBundle *b = &set->bundles[set->nbundles];
        b->component0 = link[first].component0;
        b->component1 = link[first].component1;
        b->links = count;
        b->plane_min = INT32_MAX;
        b->plane_max = INT32_MIN;

        size_t unique_planes = 0, cross_links = 0, top_links = 0;
        size_t topology_alternative_links = 0;
        int32_t previous_plane = INT32_MIN;
        double sum0 = 0.0, sum1 = 0.0, sum00 = 0.0;
        double sum11 = 0.0, sum01 = 0.0;
        for (size_t i = first; i < last; i++) {
            if (link[i].member == SIZE_MAX) continue;
            if (link[i].plane != previous_plane) {
                unique_planes++;
                previous_plane = link[i].plane;
            }
            if (link[i].plane < b->plane_min) b->plane_min = link[i].plane;
            if (link[i].plane > b->plane_max) b->plane_max = link[i].plane;
            cross_links++;
            if (link[i].top_prior) top_links++;
            if (link[i].has_topology_alternative)
                topology_alternative_links++;
            sum0 += link[i].arc0;
            sum1 += link[i].arc1;
            sum00 += link[i].arc0 * link[i].arc0;
            sum11 += link[i].arc1 * link[i].arc1;
            sum01 += link[i].arc0 * link[i].arc1;
        }
        b->planes = unique_planes;
        b->cross_links = cross_links;
        b->continuation_links = count - cross_links;
        if (cross_links == 0) {
            b->plane_min = link[first].plane;
            b->plane_max = link[last - 1].plane;
            b->initial_residual_median = NAN;
            b->initial_residual_mad = NAN;
            b->frozen_residual_median = NAN;
            b->frozen_residual_mad = NAN;
            b->prior_median = NAN;
        } else {
            size_t k = 0;
            for (size_t i = first; i < last; i++)
                if (link[i].member != SIZE_MAX)
                    scratch[k++] = link[i].initial_residual;
            b->initial_residual_median = ac_median(scratch, k);
            k = 0;
            for (size_t i = first; i < last; i++)
                if (link[i].member != SIZE_MAX)
                    scratch[k++] = fabs(link[i].initial_residual -
                                        b->initial_residual_median);
            b->initial_residual_mad = ac_median(scratch, k);
            k = 0;
            for (size_t i = first; i < last; i++)
                if (link[i].member != SIZE_MAX)
                    scratch[k++] = link[i].residual;
            b->frozen_residual_median = ac_median(scratch, k);
            k = 0;
            for (size_t i = first; i < last; i++)
                if (link[i].member != SIZE_MAX)
                    scratch[k++] = fabs(link[i].residual -
                                        b->frozen_residual_median);
            b->frozen_residual_mad = ac_median(scratch, k);
            k = 0;
            for (size_t i = first; i < last; i++)
                if (link[i].member != SIZE_MAX)
                    scratch[k++] = link[i].prior;
            b->prior_median = ac_median(scratch, k);
        }
        b->top_prior_fraction = cross_links
            ? (double)top_links / (double)cross_links : 1.0;
        b->topology_alternative_fraction = cross_links
            ? (double)topology_alternative_links / (double)cross_links : 0.0;

        double sigma2 = o->overlap_coherence_mad_limit *
                        o->overlap_coherence_mad_limit;
        for (size_t i = first; i < last; i++) {
            if (link[i].member == SIZE_MAX || cross_links == 0) continue;
            double d = link[i].residual - b->frozen_residual_median;
            b->evidence += link[i].prior * exp(-0.5 * d * d / sigma2);
        }
        if (cross_links >= 2) {
            double n = (double)cross_links;
            double covariance = n * sum01 - sum0 * sum1;
            double variance0 = n * sum00 - sum0 * sum0;
            double variance1 = n * sum11 - sum1 * sum1;
            double den = sqrt(fmax(0.0, variance0 * variance1));
            b->arc_correlation = den > 1e-12 ? covariance / den : 0.0;
        } else {
            b->arc_correlation = NAN;
        }

        size_t ngap = 0;
        for (size_t i = first; i < last; i++)
            if (link[i].member != SIZE_MAX &&
                link[i].plane == b->plane_min)
                scratch[ngap++] = link[i].gap;
        b->first_plane_gap = ac_median(scratch, ngap);
        ngap = 0;
        for (size_t i = first; i < last; i++)
            if (link[i].member != SIZE_MAX &&
                link[i].plane == b->plane_max)
                scratch[ngap++] = link[i].gap;
        b->last_plane_gap = ac_median(scratch, ngap);

        int32_t c0min = component_plane_min[b->component0];
        int32_t c0max = component_plane_max[b->component0];
        int32_t c1min = component_plane_min[b->component1];
        int32_t c1max = component_plane_max[b->component1];
        b->component0_lifespan = c0max - c0min + 1;
        b->component1_lifespan = c1max - c1min + 1;
        b->component0_extension =
            (b->plane_min > c0min ? b->plane_min - c0min : 0) +
            (c0max > b->plane_max ? c0max - b->plane_max : 0);
        b->component1_extension =
            (b->plane_min > c1min ? b->plane_min - c1min : 0) +
            (c1max > b->plane_max ? c1max - b->plane_max : 0);
        int32_t overlap_lo = c0min > c1min ? c0min : c1min;
        int32_t overlap_hi = c0max < c1max ? c0max : c1max;
        int32_t overlap_span = overlap_hi >= overlap_lo
                             ? overlap_hi - overlap_lo + 1 : 0;
        int32_t contact_span = b->plane_max - b->plane_min + 1;
        b->contact_overlap_fraction = overlap_span > 0
            ? fmin(1.0, (double)contact_span / (double)overlap_span) : 0.0;

        int32_t bundle_index = (int32_t)set->nbundles;
        for (size_t i = first; i < last; i++) {
            if (link[i].member != SIZE_MAX)
                set->member_bundle[link[i].member] = bundle_index;
            else
                set->continuation_bundle[link[i].continuation] = bundle_index;
        }
        set->nbundles++;
        first = last;
    }

    int32_t *best_bundle = (int32_t *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(int32_t));
    double *best_score = (double *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(double));
    double *second_score = (double *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(double));
    for (size_t i = 0; i < set->mesh_components; i++) {
        best_bundle[i] = -1;
        best_score[i] = second_score[i] = -DBL_MAX;
    }
    for (size_t i = 0; i < set->nbundles; i++) {
        AtlasCandidateBundle *b = &set->bundles[i];
        ac_best_bundle_update(b->component0, (int32_t)i, b->evidence,
                              best_bundle, best_score, second_score);
        ac_best_bundle_update(b->component1, (int32_t)i, b->evidence,
                              best_bundle, best_score, second_score);
    }

    for (size_t i = 0; i < set->nbundles; i++) {
        AtlasCandidateBundle *b = &set->bundles[i];
        b->mutual_best = best_bundle[b->component0] == (int32_t)i &&
                         best_bundle[b->component1] == (int32_t)i;
        b->component0_runner_up_ratio =
            best_bundle[b->component0] == (int32_t)i
            ? (second_score[b->component0] > 0.0
               ? second_score[b->component0] / best_score[b->component0]
               : 0.0)
            : (b->evidence > 0.0
               ? best_score[b->component0] / b->evidence : DBL_MAX);
        b->component1_runner_up_ratio =
            best_bundle[b->component1] == (int32_t)i
            ? (second_score[b->component1] > 0.0
               ? second_score[b->component1] / best_score[b->component1]
               : 0.0)
            : (b->evidence > 0.0
               ? best_score[b->component1] / b->evidence : DBL_MAX);
        if (b->mutual_best) set->selection.refinement_mutual_bundles++;

        int coherent_overlap =
            b->cross_links >= (size_t)o->minimum_links &&
            b->planes >= (size_t)o->minimum_planes &&
            isfinite(b->initial_residual_median) &&
            isfinite(b->initial_residual_mad) &&
            isfinite(b->frozen_residual_median) &&
            isfinite(b->frozen_residual_mad) &&
            fabs(b->initial_residual_median) >=
                o->gross_initial_residual &&
            b->initial_residual_mad <=
                o->initial_coherence_mad_limit &&
            fabs(b->frozen_residual_median) <=
                o->overlap_residual_limit &&
            b->frozen_residual_mad <=
                o->overlap_coherence_mad_limit;
        int bounded_front =
            (b->component0_extension <= o->front_maximum_extension ||
             b->component1_extension <= o->front_maximum_extension) &&
            (b->first_plane_gap <= o->front_close_gap ||
             b->last_plane_gap <= o->front_close_gap) &&
            b->contact_overlap_fraction >=
                o->front_minimum_contact_overlap;
        int persistent_parallel = !bounded_front &&
            b->contact_overlap_fraction >=
                o->persistent_contact_fraction;
        int same_sheet_seam = persistent_parallel &&
            b->continuation_links > 0;
        int missed_sheet = persistent_parallel &&
            b->continuation_links == 0 &&
            b->topology_alternative_fraction >=
                o->minimum_topology_alternative_fraction;

        if (coherent_overlap && bounded_front &&
            b->cross_links >= (size_t)o->front_minimum_links &&
            b->planes >= (size_t)o->front_minimum_planes) {
            b->decision = ATLAS_BUNDLE_DELAMINATION;
            set->selection.refinement_delamination_bundles++;
        } else if (coherent_overlap && same_sheet_seam) {
            b->decision = ATLAS_BUNDLE_SAME_SHEET_SEAM;
            set->selection.refinement_sheet_seam_bundles++;
        } else if (coherent_overlap && missed_sheet) {
            b->decision = ATLAS_BUNDLE_SEPARATE_SHEET;
            set->selection.refinement_separate_sheet_bundles++;
        } else {
            b->decision = ATLAS_BUNDLE_UNCLASSIFIED;
            set->selection.refinement_inconclusive_bundles++;
        }
    }
    set->selection.refinement_bundles = set->nbundles;

    size_t remaining = 0;
    for (size_t c = 0; c < p->ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        int active = 0;
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            int32_t bi = set->member_bundle[mi];
            AtlasCandidateBundleDecision decision = bi >= 0
                ? set->bundles[bi].decision : ATLAS_BUNDLE_UNCLASSIFIED;
            if (decision == ATLAS_BUNDLE_SEPARATE_SHEET) {
                set->members[mi].membership = 0.0;
                set->member_state[mi] = ATLAS_CANDIDATE_INACTIVE;
                set->selection.refinement_pruned_candidate_links++;
                continue;
            }
            if (decision == ATLAS_BUNDLE_DELAMINATION)
                set->member_state[mi] =
                    ATLAS_CANDIDATE_REFINED_DELAMINATION;
            else if (decision == ATLAS_BUNDLE_SAME_SHEET_SEAM)
                set->member_state[mi] =
                    ATLAS_CANDIDATE_REFINED_SHEET_SEAM;
            active = 1;
        }
        if (active) remaining++;
    }
    for (size_t i = 0; i < p->ncontinuations; i++) {
        int32_t bi = set->continuation_bundle[i];
        if (bi < 0) continue;
        AtlasCandidateBundleDecision decision = set->bundles[bi].decision;
        int cross_cube = set->sample_ref[set->continuations[i].a].cube !=
                         set->sample_ref[set->continuations[i].b].cube;
        if (decision == ATLAS_BUNDLE_SEPARATE_SHEET && !cross_cube) {
            set->continuations[i].membership = 0.0;
            set->continuation_state[i] = ATLAS_CANDIDATE_INACTIVE;
            set->selection.refinement_pruned_continuations++;
        } else if (decision == ATLAS_BUNDLE_DELAMINATION) {
            set->continuation_state[i] =
                ATLAS_CANDIDATE_REFINED_DELAMINATION;
        } else if (decision == ATLAS_BUNDLE_SAME_SHEET_SEAM) {
            set->continuation_state[i] =
                ATLAS_CANDIDATE_REFINED_SHEET_SEAM;
        }
    }
    set->selection.refinement_remaining_active_cross_sections = remaining;
    return ac_rebuild_support(arena, set, o->support_floor);
}

static void ac_selftest_check(int condition, const char *what, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[atlas_candidates selftest] FAIL: %s\n", what);
    (*fails)++;
}

int AtlasCandidates_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* Two disjoint z-by-u rectangles with a 0.25-voxel same-slice gap. */
    float vertices[8 * 3] = {
         0,  0.00f, 0,  16,  0.00f, 0,
         0,  5.00f, 0,  16,  5.00f, 0,
         0,  5.25f, 0,  16,  5.25f, 0,
         0, 10.25f, 0,  16, 10.25f, 0
    };
    int32_t triangles[4 * 3] = {
        0, 1, 2, 1, 3, 2,
        4, 5, 6, 5, 7, 6
    };
    int32_t face_cube[4] = {0, 0, 1, 1};
    float raw_u[8] = {
        0, 0, 5, 5,
        25.25f, 25.25f, 30.25f, 30.25f
    };
    float axis_point[3] = {0, 0, 0};
    float axis_dir[3] = {1, 0, 0};

    AtlasCandidateOptions co;
    AtlasCandidateOptions_default(&co);
    co.slice_spacing = 4.0;
    co.min_stroke_length = 2.0;
    co.match_radius = 1.0;
    co.continuation_radius = 1.0;
    co.max_candidates = 3;

    AtlasCandidateSet set;
    int build_rc = AtlasCandidates_build(
        arena, vertices, 8, triangles, face_cube, 4, raw_u,
        axis_point, axis_dir, &co, &set);
    ac_selftest_check(build_rc == 0, "build", &fails);
    if (build_rc == 0) {
        ac_selftest_check(set.problem.nstrokes == 8,
                          "two strokes on four planes", &fails);
        ac_selftest_check(set.problem.ncontinuations == 4,
                          "one mutual continuation per plane", &fails);
        ac_selftest_check(set.mesh_components == 2,
                          "two retained mesh components", &fails);
        int baseline_rc = AtlasCandidates_select_baseline(
            arena, 1e-4, &set);
        ac_selftest_check(baseline_rc == 0,
                          "HEAD baseline selection", &fails);
        ac_selftest_check(set.problem.nanchors == 1,
                          "HEAD baseline retains the overlap assignment", &fails);
        ac_selftest_check(
            set.selection.candidate_links > 0,
            "HEAD baseline retains candidate links", &fails);

        double *frozen = (double *)ARENA_ALLOC(
            arena, set.problem.nsamples * sizeof(double));
        for (size_t i = 0; i < set.problem.nsamples; i++)
            frozen[i] = set.samples[i].p[1];
        AtlasCandidateRefineOptions refine_opts;
        AtlasCandidateRefineOptions_default(&refine_opts);
        refine_opts.gross_initial_residual = 1.0;
        refine_opts.initial_coherence_mad_limit = 64.0;
        refine_opts.minimum_links = 1;
        refine_opts.minimum_planes = 1;
        refine_opts.front_minimum_links = 1;
        refine_opts.front_minimum_planes = 1;
        int refine_rc = AtlasCandidates_select_refined(
            arena, frozen, &refine_opts, &set);
        ac_selftest_check(refine_rc == 0, "overlap classification", &fails);
        ac_selftest_check(set.selection.refinement_delamination_bundles == 1,
                          "bounded coherent overlap is delamination", &fails);
        ac_selftest_check(
            set.selection.refinement_pruned_candidate_links == 0 &&
            set.selection.refinement_pruned_continuations == 0,
            "delamination retains every HEAD assignment", &fails);
        ac_selftest_check(set.problem.nanchors == 1,
                          "delamination keeps the HEAD gauge", &fails);

        AtlasCandidateRefineOptions seam_opts = refine_opts;
        seam_opts.front_close_gap = 0.1;
        int seam_rc = AtlasCandidates_select_refined(
            arena, frozen, &seam_opts, &set);
        ac_selftest_check(seam_rc == 0,
                          "persistent seam classification", &fails);
        ac_selftest_check(
            set.selection.refinement_sheet_seam_bundles == 1,
            "continuation-supported overlap is a same-sheet seam", &fails);
        ac_selftest_check(
            set.selection.refinement_pruned_candidate_links == 0,
            "same-sheet seam retains HEAD overlap assignments", &fails);

        size_t saved_continuations = set.problem.ncontinuations;
        set.problem.ncontinuations = 0;
        AtlasCandidateRefineOptions separate_opts = seam_opts;
        separate_opts.minimum_topology_alternative_fraction = 0.0;
        int separate_rc = AtlasCandidates_select_refined(
            arena, frozen, &separate_opts, &set);
        ac_selftest_check(separate_rc == 0,
                          "missed-sheet classification", &fails);
        ac_selftest_check(
            set.selection.refinement_separate_sheet_bundles == 1 &&
            set.selection.refinement_pruned_candidate_links > 0,
            "unsupported persistent overlap is removed as a missed sheet",
            &fails);
        set.problem.ncontinuations = saved_continuations;

        refine_rc = AtlasCandidates_select_refined(
            arena, frozen, &refine_opts, &set);
        ac_selftest_check(refine_rc == 0,
                          "restore delamination fixture", &fails);

        size_t nvar = set.problem.nsamples +
                      set.problem.ncross_sections;
        double *relaxed = (double *)ARENA_ALLOC(
            arena, nvar * sizeof(double));
        memcpy(relaxed, set.initial_u,
               set.problem.nsamples * sizeof(double));

        AtlasStripOptions so;
        AtlasStripOptions_default(&so);
        so.lambda_align = 4.0;
        so.lambda_continuation = 4.0;
        AtlasStrip_initialize_intercepts(&set.problem, &so, relaxed);
        AtlasStripRobustOptions ro;
        AtlasStripRobustOptions_default(&ro);
        ro.l1_iterations = 2;
        ro.likelihood_iterations = 4;
        ro.likelihood_sigma = 0.5;
        double *membership = (double *)ARENA_ALLOC(
            arena, set.problem.nmembers * sizeof(double));
        MonotoneQpOptions qo;
        MonotoneQpOptions_default(&qo);
        int solve_rc = AtlasStrip_solve_robust(
            arena, &set.problem, &so, &ro, &qo, relaxed,
            membership, NULL, NULL);
        ac_selftest_check(solve_rc == 0, "robust relaxed solve", &fails);

        AtlasStripMember *final_member = (AtlasStripMember *)ARENA_ALLOC(
            arena, set.problem.nmembers * sizeof(AtlasStripMember));
        memcpy(final_member, set.problem.members,
               set.problem.nmembers * sizeof(AtlasStripMember));
        for (size_t i = 0; i < set.problem.nmembers; i++)
            final_member[i].membership = membership[i];
        AtlasStripProblem final_problem = set.problem;
        final_problem.members = final_member;
        AtlasStripOptions final_opts = so;
        final_opts.mode = ATLAS_STRIP_FINAL;
        double *final = (double *)ARENA_ALLOC(
            arena, nvar * sizeof(double));
        memcpy(final, relaxed, nvar * sizeof(double));
        AtlasStrip_initialize_intercepts(
            &final_problem, &final_opts, final);
        AtlasStripSystem final_system;
        int final_build = AtlasStrip_build(
            arena, &final_problem, &final_opts, &final_system);
        MonotoneQpStats final_stats;
        int final_rc = final_build == 0
                     ? MonotoneQp_solve(arena, &final_system.qp, &qo,
                                        final, NULL, &final_stats)
                     : -1;
        ac_selftest_check(final_rc == 0, "final solve", &fails);

        double max_continuation = 0.0;
        for (size_t i = 0; i < final_problem.ncontinuations; i++) {
            const AtlasStripContinuation *cn =
                &final_problem.continuations[i];
            double residual =
                fabs(final[cn->a] - final[cn->b] - cn->target);
            if (residual > max_continuation)
                max_continuation = residual;
        }
        ac_selftest_check(max_continuation < 1e-8,
                          "continuation recovers the shifted fragment",
                          &fails);
        fprintf(stderr,
                "[atlas_candidates selftest] strokes=%zu samples=%zu "
                "cross=%zu cont=%zu components=%zu max_cont=%.3e\n",
                set.problem.nstrokes, set.problem.nsamples,
                set.problem.ncross_sections,
                set.problem.ncontinuations,
                set.problem.nanchors, max_continuation);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_candidates selftest] %s\n",
            fails ? "FAILED" : "PASSED");
    return fails;
}
