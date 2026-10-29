/* marble_strip_uv.c -- StrokeStrip variational U repair for a quad ribbon. */

#include "marble_strip_uv.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct CoarseStrip {
    int H, W;
    int *row, *col;
    AtlasStripSample *sample;
    AtlasStripStroke *stroke;
    AtlasStripMember *member;
    AtlasStripCrossSection *cross;
    MonotoneQpAnchor *anchor;
    double *prior_u;
    double *prior_weight;
    double *sample_trust;
    double *x;
    double *membership;
    double *sample_likelihood;
    size_t nsample, nmember, ncross, nanchor, x_capacity;
} CoarseStrip;

enum {
    FINAL_CROSS_WINDOW_ROWS = 17,
    FINAL_CROSS_STEP_ROWS = 8
};

static double clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

static double norm3(const double v[3])
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static int normalize3(double v[3])
{
    double n = norm3(v);
    if (!isfinite(n) || n <= 1e-12) return -1;
    v[0] /= n; v[1] /= n; v[2] /= n;
    return 0;
}

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int compare_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

static double median_in_place(double *value, size_t n)
{
    if (n == 0) return NAN;
    qsort(value, n, sizeof *value, compare_double);
    return n & 1u ? value[n / 2] :
           0.5 * (value[n / 2 - 1] + value[n / 2]);
}

static int build_indices(int fine_n, int stride, int **out, int *out_n)
{
    int capacity, n = 0;
    int *index;
    if (fine_n < 2 || stride < 1 || out == NULL || out_n == NULL) return -1;
    capacity = (fine_n - 1 + stride - 1) / stride + 2;
    index = (int *)malloc((size_t)capacity * sizeof *index);
    if (index == NULL) return -1;
    for (int i = 0; i < fine_n - 1; i += stride) index[n++] = i;
    if (n == 0 || index[n - 1] != fine_n - 1) index[n++] = fine_n - 1;
    *out = index; *out_n = n;
    return 0;
}

static int interval_for(const int *index, int n, int x)
{
    int lo = 0, hi = n - 1;
    if (x <= index[0]) return 0;
    if (x >= index[n - 1]) return n - 2;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (index[mid] <= x) lo = mid; else hi = mid;
    }
    return lo;
}

static int build_cross_windows(int H, int **out_start, int *out_n)
{
    int *start, n = 0, capacity;
    if (H < 1 || out_start == NULL || out_n == NULL) return -1;
    capacity = (H + FINAL_CROSS_STEP_ROWS - 1) / FINAL_CROSS_STEP_ROWS + 2;
    start = (int *)malloc((size_t)capacity * sizeof *start);
    if (start == NULL) return -1;
    if (H <= FINAL_CROSS_WINDOW_ROWS) {
        start[n++] = 0;
    } else {
        int last = H - FINAL_CROSS_WINDOW_ROWS;
        for (int s = 0; s < last; s += FINAL_CROSS_STEP_ROWS) start[n++] = s;
        if (n == 0 || start[n - 1] != last) start[n++] = last;
    }
    *out_start = start; *out_n = n;
    return 0;
}

static void coarse_dispose(CoarseStrip *c)
{
    if (c == NULL) return;
    free(c->row); free(c->col); free(c->sample); free(c->stroke);
    free(c->member); free(c->cross); free(c->anchor);
    free(c->prior_u); free(c->prior_weight); free(c->sample_trust);
    free(c->x); free(c->membership); free(c->sample_likelihood);
    memset(c, 0, sizeof *c);
}

void MarbleStripUvOptions_default(MarbleStripUvOptions *o)
{
    if (o == NULL) return;
    memset(o, 0, sizeof *o);
    o->sample_stride_u = 4;
    o->sample_stride_v = 4;
    o->l1_iterations = 8;
    o->max_irls_iterations = 1000;
    o->final_iterations = 0;
    o->final_damping = 0.1;
    o->movement_tolerance = 1e-5;
    o->likelihood_sigma = 0.0;
    o->lambda_length = 1.0;
    o->lambda_align = 1.0;
    o->lambda_local = 0.05;
    o->lambda_prior = 100.0;
    o->monotone_fraction = 0.05;
    o->membership_floor = 0.05;
    o->match_radius = 12.0;
    o->match_angle_deg = 35.0;
    o->match_search_columns = 24;
    o->topology_fallback = 0.02;
    o->length_winsor = 3.0;
    o->max_displacement = 64.0;
    o->verbose = 0;
}

/* Robust latent length per topology-column edge, followed by a bounded use of
 * each row's own XYZ arclength.  The per-row normalization keeps the two clean
 * crop seams compatible while redistributing U toward locally compressed
 * regions instead of spreading their excess length over every row. */
static int effective_lengths(const float *verts, const double *uv,
                             const double *trust, int H, int W,
                             double winsor, double **out_length)
{
    size_t ne = (size_t)H * (size_t)(W - 1);
    double *length = (double *)malloc(ne * sizeof *length);
    double *latent = (double *)malloc((size_t)(W - 1) * sizeof *latent);
    double *tmp = (double *)malloc((size_t)H * sizeof *tmp);
    if (length == NULL || latent == NULL || tmp == NULL) {
        free(length); free(latent); free(tmp); return -1;
    }
    for (int c = 0; c + 1 < W; c++) {
        size_t n = 0;
        for (int r = 0; r < H; r++) {
            size_t a = (size_t)r * (size_t)W + (size_t)c, b = a + 1;
            double t = fmin(clamp01(trust[a]), clamp01(trust[b]));
            double dz = (double)verts[b * 3] - (double)verts[a * 3];
            double dy = (double)verts[b * 3 + 1] - (double)verts[a * 3 + 1];
            double dx = (double)verts[b * 3 + 2] - (double)verts[a * 3 + 2];
            double d = sqrt(dz * dz + dy * dy + dx * dx);
            if (t >= 0.5 && isfinite(d) && d > 1e-6) tmp[n++] = d;
        }
        if (n < 3) {
            n = 0;
            for (int r = 0; r < H; r++) {
                size_t a = (size_t)r * (size_t)W + (size_t)c, b = a + 1;
                double d = fabs(uv[b * 2] - uv[a * 2]);
                if (isfinite(d) && d > 1e-6) tmp[n++] = d;
            }
        }
        latent[c] = median_in_place(tmp, n);
        if (!isfinite(latent[c]) || latent[c] <= 1e-6) latent[c] = 1.0;
    }
    for (int r = 0; r < H; r++) {
        double sum = 0.0;
        for (int c = 0; c + 1 < W; c++) {
            size_t a = (size_t)r * (size_t)W + (size_t)c, b = a + 1;
            double dz = (double)verts[b * 3] - (double)verts[a * 3];
            double dy = (double)verts[b * 3 + 1] - (double)verts[a * 3 + 1];
            double dx = (double)verts[b * 3 + 2] - (double)verts[a * 3 + 2];
            double raw = sqrt(dz * dz + dy * dy + dx * dx);
            double t = fmin(clamp01(trust[a]), clamp01(trust[b]));
            double cap = 1.5 + (winsor - 1.5) * (1.0 - t);
            if (cap < 1.05) cap = 1.05;
            if (!isfinite(raw) || raw <= 1e-6) raw = latent[c];
            if (raw < latent[c] / cap) raw = latent[c] / cap;
            if (raw > latent[c] * cap) raw = latent[c] * cap;
            length[(size_t)r * (size_t)(W - 1) + (size_t)c] = raw;
            sum += raw;
        }
        double span = uv[((size_t)r * (size_t)W + (size_t)(W - 1)) * 2] -
                      uv[((size_t)r * (size_t)W) * 2];
        if (!isfinite(span) || span <= 1e-6 || sum <= 1e-12) {
            free(length); free(latent); free(tmp); return -1;
        }
        double scale = span / sum;
        for (int c = 0; c + 1 < W; c++)
            length[(size_t)r * (size_t)(W - 1) + (size_t)c] *= scale;
    }
    free(latent); free(tmp);
    *out_length = length;
    return 0;
}

static double sample_base_membership(double trust, double floor_value)
{
    trust = clamp01(trust);
    return floor_value + (1.0 - floor_value) * trust;
}

static int coarse_build(const float *verts, const double *uv,
                        const double *trust, const uint8_t *pin,
                        int H, int W, const MarbleStripUvOptions *o,
                        CoarseStrip *c)
{
    double *edge_length = NULL, *raw_tangent = NULL;
    size_t cross_capacity, member_capacity;
    memset(c, 0, sizeof *c);
    if (build_indices(H, o->sample_stride_v, &c->row, &c->H) != 0 ||
        build_indices(W, o->sample_stride_u, &c->col, &c->W) != 0 ||
        effective_lengths(verts, uv, trust, H, W, o->length_winsor,
                          &edge_length) != 0)
        goto fail;
    c->nsample = (size_t)c->H * (size_t)c->W;
    if ((size_t)(c->H - 1) > SIZE_MAX / (size_t)c->W) goto fail;
    cross_capacity = (size_t)(c->H - 1) * (size_t)c->W;
    if (cross_capacity > SIZE_MAX / 3u) goto fail;
    member_capacity = cross_capacity * 3u;
    c->sample = (AtlasStripSample *)calloc(c->nsample, sizeof *c->sample);
    c->stroke = (AtlasStripStroke *)calloc((size_t)c->H, sizeof *c->stroke);
    c->member = (AtlasStripMember *)calloc(member_capacity, sizeof *c->member);
    c->cross = (AtlasStripCrossSection *)calloc(cross_capacity, sizeof *c->cross);
    c->anchor = (MonotoneQpAnchor *)malloc(c->nsample * sizeof *c->anchor);
    c->prior_u = (double *)malloc(c->nsample * sizeof *c->prior_u);
    c->prior_weight = (double *)malloc(c->nsample * sizeof *c->prior_weight);
    c->sample_trust = (double *)malloc(c->nsample * sizeof *c->sample_trust);
    c->x_capacity = c->nsample + cross_capacity;
    c->x = (double *)calloc(c->x_capacity, sizeof *c->x);
    c->membership = (double *)malloc(member_capacity * sizeof *c->membership);
    c->sample_likelihood = (double *)malloc(c->nsample * sizeof *c->sample_likelihood);
    raw_tangent = (double *)malloc(c->nsample * 3 * sizeof *raw_tangent);
    if (c->sample == NULL || c->stroke == NULL || c->member == NULL ||
        c->cross == NULL || c->anchor == NULL || c->prior_u == NULL ||
        c->prior_weight == NULL || c->sample_trust == NULL || c->x == NULL ||
        c->membership == NULL || c->sample_likelihood == NULL ||
        raw_tangent == NULL) goto fail;

    for (int rr = 0; rr < c->H; rr++) {
        int fr = c->row[rr];
        c->stroke[rr].first = rr * c->W;
        c->stroke[rr].count = c->W;
        c->stroke[rr].component = 0;
        c->stroke[rr].resolved = 1;
        double s = 0.0;
        for (int cc = 0; cc < c->W; cc++) {
            int fc = c->col[cc];
            size_t fine = (size_t)fr * (size_t)W + (size_t)fc;
            size_t i = (size_t)rr * (size_t)c->W + (size_t)cc;
            if (cc > 0) {
                int prev = c->col[cc - 1];
                for (int x = prev; x < fc; x++)
                    s += edge_length[(size_t)fr * (size_t)(W - 1) + (size_t)x];
            }
            for (int k = 0; k < 3; k++) c->sample[i].p[k] = verts[fine * 3 + (size_t)k];
            c->sample[i].s = s;
            c->sample[i].stroke = rr;
            c->sample[i].ordinal = cc;
            c->sample_trust[i] = clamp01(trust[fine]);
            c->sample[i].confidence = 0.25 + 0.75 * c->sample_trust[i];
            c->prior_u[i] = uv[fine * 2];
            c->prior_weight[i] = c->sample_trust[i] * c->sample_trust[i];
            c->x[i] = c->prior_u[i];
            c->sample_likelihood[i] = 1.0;
            if ((pin != NULL && pin[fine] != 0) || cc == 0 || cc + 1 == c->W) {
                c->anchor[c->nanchor].var = (int32_t)i;
                c->anchor[c->nanchor].value = c->prior_u[i];
                c->anchor[c->nanchor].component = 0;
                c->nanchor++;
            }
        }
    }
    free(edge_length); edge_length = NULL;

    /* Raw row tangents, followed by a clean-row column consensus.  The latter
     * is the better angle constraint in a marble core: a crumpled observation
     * cannot vote its own tangent back into the latent strip. */
    for (int rr = 0; rr < c->H; rr++) for (int cc = 0; cc < c->W; cc++) {
        int lo = cc > 0 ? cc - 1 : cc;
        int hi = cc + 1 < c->W ? cc + 1 : cc;
        size_t i = (size_t)rr * (size_t)c->W + (size_t)cc;
        size_t a = (size_t)rr * (size_t)c->W + (size_t)lo;
        size_t b = (size_t)rr * (size_t)c->W + (size_t)hi;
        double *t = raw_tangent + i * 3;
        for (int k = 0; k < 3; k++) t[k] = c->sample[b].p[k] - c->sample[a].p[k];
        if (normalize3(t) != 0) { t[0] = 0.0; t[1] = 0.0; t[2] = 1.0; }
    }
    for (int cc = 0; cc < c->W; cc++) {
        double consensus[3] = {0,0,0}, weight_sum = 0.0;
        for (int rr = 0; rr < c->H; rr++) {
            size_t i = (size_t)rr * (size_t)c->W + (size_t)cc;
            double w = c->sample_trust[i] * c->sample_trust[i];
            for (int axis = 0; axis < 3; axis++)
                consensus[axis] += w * raw_tangent[i * 3 + (size_t)axis];
            weight_sum += w;
        }
        if (weight_sum <= 1e-12 || normalize3(consensus) != 0) {
            memcpy(consensus, raw_tangent + (size_t)cc * 3, sizeof consensus);
            if (normalize3(consensus) != 0) {
                consensus[0] = 0.0; consensus[1] = 0.0; consensus[2] = 1.0;
            }
        }
        for (int rr = 0; rr < c->H; rr++) {
            size_t i = (size_t)rr * (size_t)c->W + (size_t)cc;
            double q = c->sample_trust[i];
            for (int axis = 0; axis < 3; axis++)
                c->sample[i].tangent[axis] =
                    q * raw_tangent[i * 3 + (size_t)axis] +
                    (1.0 - q) * consensus[axis];
            if (normalize3(c->sample[i].tangent) != 0)
                memcpy(c->sample[i].tangent, consensus, sizeof consensus);
        }
    }
    free(raw_tangent); raw_tangent = NULL;

    /* StrokeStrip Sec. 5.1: a source point connects to a geometrically
     * adjacent point on the next stroke, not blindly to the same topology
     * column.  Search target-row segments in the plane normal to axial Z and
     * retain a low-prior topology alternative when geometry proposes a shift. */
    const double cosine_gate = cos(o->match_angle_deg *
                                   0.01745329251994329577);
    const double sigma = 0.5 * o->match_radius;
    const double radius2 = o->match_radius * o->match_radius;
    for (int rr = 0; rr + 1 < c->H; rr++) for (int cc = 0; cc < c->W; cc++) {
        size_t source = (size_t)rr * (size_t)c->W + (size_t)cc;
        int k0 = cc - o->match_search_columns;
        int k1 = cc + o->match_search_columns;
        int best_k = -1;
        double best_t = 0.0, best_score = -1.0;
        if (k0 < 0) k0 = 0;
        if (k1 > c->W - 2) k1 = c->W - 2;
        for (int k = k0; k <= k1; k++) {
            size_t a = (size_t)(rr + 1) * (size_t)c->W + (size_t)k;
            size_t b = a + 1;
            double ab1 = c->sample[b].p[1] - c->sample[a].p[1];
            double ab2 = c->sample[b].p[2] - c->sample[a].p[2];
            double aq1 = c->sample[source].p[1] - c->sample[a].p[1];
            double aq2 = c->sample[source].p[2] - c->sample[a].p[2];
            double denom = ab1 * ab1 + ab2 * ab2;
            if (denom <= 1e-12) continue;
            double t = (aq1 * ab1 + aq2 * ab2) / denom;
            if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
            double d1 = c->sample[a].p[1] + t * ab1 - c->sample[source].p[1];
            double d2 = c->sample[a].p[2] + t * ab2 - c->sample[source].p[2];
            double lateral2 = d1 * d1 + d2 * d2;
            if (lateral2 > radius2) continue;
            double target_tangent[3];
            for (int axis = 0; axis < 3; axis++)
                target_tangent[axis] =
                    (1.0 - t) * c->sample[a].tangent[axis] +
                    t * c->sample[b].tangent[axis];
            if (normalize3(target_tangent) != 0) continue;
            double tangent_dot = dot3(c->sample[source].tangent, target_tangent);
            if (tangent_dot < cosine_gate) continue;
            double score = exp(-0.5 * lateral2 / (sigma * sigma)) *
                           tangent_dot * tangent_dot * tangent_dot * tangent_dot;
            if (score > best_score) {
                best_score = score; best_k = k; best_t = t;
            }
        }
        int topology_k = cc < c->W - 1 ? cc : c->W - 2;
        double topology_t = cc < c->W - 1 ? 0.0 : 1.0;
        int add_topology = best_k < 0 ||
            fabs(((double)best_k + best_t) - (double)cc) > 0.5;
        int candidate_count = 1 + (best_k >= 0 && add_topology ? 1 : 0);
        if (c->ncross >= cross_capacity ||
            c->nmember + (size_t)(1 + candidate_count) > member_capacity)
            goto fail;
        AtlasStripCrossSection *cs = &c->cross[c->ncross];
        memset(cs, 0, sizeof *cs);
        cs->first = c->nmember;
        cs->count = 1 + candidate_count;
        cs->weight = 1.0;
        cs->id = (int32_t)c->ncross;
        double tangent_sum[3];
        memcpy(tangent_sum, c->sample[source].tangent, sizeof tangent_sum);

        AtlasStripMember *m = &c->member[c->nmember++];
        int sc0 = cc > 0 ? cc - 1 : cc;
        int sc1 = cc + 1 < c->W ? cc + 1 : cc;
        size_t sdlo = (size_t)rr * (size_t)c->W + (size_t)sc0;
        size_t sdhi = (size_t)rr * (size_t)c->W + (size_t)sc1;
        m->value0 = (int32_t)source; m->value1 = -1; m->value_t = 0.0;
        m->deriv_lo = (int32_t)sdlo; m->deriv_hi = (int32_t)sdhi;
        m->deriv_length = c->sample[sdhi].s - c->sample[sdlo].s;
        memcpy(m->p, c->sample[source].p, sizeof m->p);
        memcpy(m->tangent, c->sample[source].tangent, sizeof m->tangent);
        m->dual_width = 1.0; m->membership = 1.0; m->base_weight = 1.0;
        m->observation = (int32_t)(c->nmember - 1);

        for (int which = 0; which < candidate_count; which++) {
            int k = which == 0 && best_k >= 0 ? best_k : topology_k;
            double t = which == 0 && best_k >= 0 ? best_t : topology_t;
            double geometry = which == 0 && best_k >= 0 ? best_score :
                              o->topology_fallback;
            size_t a = (size_t)(rr + 1) * (size_t)c->W + (size_t)k;
            size_t b = a + 1;
            m = &c->member[c->nmember++];
            m->value0 = (int32_t)a; m->value1 = (int32_t)b; m->value_t = t;
            m->deriv_lo = (int32_t)a; m->deriv_hi = (int32_t)b;
            m->deriv_length = c->sample[b].s - c->sample[a].s;
            double target_trust = (1.0 - t) * c->sample_trust[a] +
                                  t * c->sample_trust[b];
            double pair_trust = fmin(c->sample_trust[source], target_trust);
            for (int axis = 0; axis < 3; axis++) {
                m->p[axis] = (1.0 - t) * c->sample[a].p[axis] +
                             t * c->sample[b].p[axis];
                m->tangent[axis] =
                    (1.0 - t) * c->sample[a].tangent[axis] +
                    t * c->sample[b].tangent[axis];
            }
            if (normalize3(m->tangent) != 0)
                memcpy(m->tangent, c->sample[a].tangent, sizeof m->tangent);
            double confidence = 0.25 + 0.75 * pair_trust;
            m->membership = fmax(o->membership_floor,
                                 fmin(1.0, geometry * confidence));
            m->dual_width = 1.0; m->base_weight = 1.0;
            m->observation = (int32_t)(c->nmember - 1);
            for (int axis = 0; axis < 3; axis++)
                tangent_sum[axis] += m->membership * m->tangent[axis];
        }
        if (normalize3(tangent_sum) != 0)
            memcpy(tangent_sum, c->sample[source].tangent, sizeof tangent_sum);
        memcpy(cs->tangent, tangent_sum, sizeof cs->tangent);
        c->ncross++;
    }
    if (c->ncross == 0 || c->nmember == 0) goto fail;
    return 0;
fail:
    free(edge_length); free(raw_tangent); coarse_dispose(c); return -1;
}

static void problem_from_coarse(const CoarseStrip *c, AtlasStripProblem *p)
{
    memset(p, 0, sizeof *p);
    p->samples = c->sample; p->nsamples = c->nsample;
    p->strokes = c->stroke; p->nstrokes = (size_t)c->H;
    p->members = c->member; p->nmembers = c->nmember;
    p->cross_sections = c->cross; p->ncross_sections = c->ncross;
    p->anchors = c->anchor; p->nanchors = c->nanchor;
    p->prior_u = c->prior_u; p->prior_weight = c->prior_weight;
}

static int bracket_u(const double *u, int first, int count, double value)
{
    int lo = 0, hi = count - 1;
    if (value <= u[first]) return 0;
    if (value >= u[first + count - 1]) return count - 2;
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        if (u[first + mid] <= value) lo = mid; else hi = mid;
    }
    return lo;
}

/* Rebuild C(u) from the current relaxed isovalues, as in StrokeStrip Sec. 5.3. */
static int rebuild_isolines(CoarseStrip *c, const MarbleStripUvOptions *o,
                            const int *window_start, int nwindow,
                            AtlasStripMember *member, size_t nmember,
                            AtlasStripCrossSection *cross, size_t ncross)
{
    double u0 = -DBL_MAX, u1 = DBL_MAX;
    size_t member_offset = 0, cross_offset = 0;
    for (int rr = 0; rr < c->H; rr++) {
        size_t first = (size_t)rr * (size_t)c->W;
        if (c->x[first] > u0) u0 = c->x[first];
        if (c->x[first + (size_t)c->W - 1] < u1)
            u1 = c->x[first + (size_t)c->W - 1];
    }
    if (!isfinite(u0) || !isfinite(u1) || u1 <= u0) return -1;
    for (int cc = 0; cc < c->W; cc++) {
        double iso = u0 + (u1 - u0) * (double)cc / (double)(c->W - 1);
        for (int ww = 0; ww < nwindow; ww++) {
            int first_row = window_start[ww];
            int count = c->H - first_row;
            double consensus[3] = {0,0,0}, consensus_weight = 0.0;
            double ct[3] = {0,0,0};
            if (count > FINAL_CROSS_WINDOW_ROWS) count = FINAL_CROSS_WINDOW_ROWS;
            if (count < 1 || cross_offset >= ncross ||
                member_offset + (size_t)count > nmember) return -1;
            AtlasStripCrossSection *cs = &cross[cross_offset];
            memset(cs, 0, sizeof *cs);
            cs->first = member_offset;
            cs->count = count;
            cs->weight = fmin(1.0, (double)FINAL_CROSS_STEP_ROWS /
                                   (double)count);
            cs->id = (int32_t)cross_offset;

            /* A surface ruling may curve over the full scroll height.  Fit the
             * paper's tangent on a narrow overlapping axial neighborhood. */
            for (int local = 0; local < count; local++) {
                int rr = first_row + local, first = rr * c->W;
                int k = bracket_u(c->x, first, c->W, iso);
                size_t lo = (size_t)first + (size_t)k, hi = lo + 1;
                double denom = c->x[hi] - c->x[lo];
                double t = denom > 1e-12 ? (iso - c->x[lo]) / denom : 0.0;
                double raw[3];
                t = clamp01(t);
                double trust = (1.0 - t) * c->sample_trust[lo] +
                               t * c->sample_trust[hi];
                double w = trust * trust;
                for (int axis = 0; axis < 3; axis++)
                    raw[axis] = (1.0 - t) * c->sample[lo].tangent[axis] +
                                t * c->sample[hi].tangent[axis];
                if (normalize3(raw) != 0) continue;
                for (int axis = 0; axis < 3; axis++) consensus[axis] += w * raw[axis];
                consensus_weight += w;
            }
            if (consensus_weight <= 1e-12 || normalize3(consensus) != 0) {
                int rr = first_row, first = rr * c->W;
                int k = bracket_u(c->x, first, c->W, iso);
                memcpy(consensus, c->sample[(size_t)first + (size_t)k].tangent,
                       sizeof consensus);
                if (normalize3(consensus) != 0) {
                    consensus[0] = 0.0; consensus[1] = 0.0; consensus[2] = 1.0;
                }
            }
            for (int local = 0; local < count; local++) {
                int rr = first_row + local, first = rr * c->W;
                int k = bracket_u(c->x, first, c->W, iso);
                size_t lo = (size_t)first + (size_t)k, hi = lo + 1;
                double denom = c->x[hi] - c->x[lo];
                double t = denom > 1e-12 ? (iso - c->x[lo]) / denom : 0.0;
                double raw[3];
                t = clamp01(t);
                AtlasStripMember *m = &member[member_offset + (size_t)local];
                memset(m, 0, sizeof *m);
                m->value0 = (int32_t)lo; m->value1 = (int32_t)hi; m->value_t = t;
                m->deriv_lo = (int32_t)lo; m->deriv_hi = (int32_t)hi;
                m->deriv_length = c->sample[hi].s - c->sample[lo].s;
                double trust = (1.0 - t) * c->sample_trust[lo] +
                               t * c->sample_trust[hi];
                double like = (1.0 - t) * c->sample_likelihood[lo] +
                              t * c->sample_likelihood[hi];
                for (int axis = 0; axis < 3; axis++) {
                    m->p[axis] = (1.0 - t) * c->sample[lo].p[axis] +
                                 t * c->sample[hi].p[axis];
                    raw[axis] = (1.0 - t) * c->sample[lo].tangent[axis] +
                                t * c->sample[hi].tangent[axis];
                }
                if (normalize3(raw) != 0) memcpy(raw, consensus, sizeof raw);
                for (int axis = 0; axis < 3; axis++)
                    m->tangent[axis] = trust * raw[axis] +
                                       (1.0 - trust) * consensus[axis];
                if (normalize3(m->tangent) != 0)
                    memcpy(m->tangent, consensus, sizeof m->tangent);
                m->dual_width = 1.0;
                m->membership = sample_base_membership(trust, o->membership_floor) *
                                fmax(like, 1e-7);
                m->base_weight = 1.0;
                m->observation = (int32_t)(member_offset + (size_t)local);
                for (int axis = 0; axis < 3; axis++)
                    ct[axis] += m->membership * m->tangent[axis];
            }
            if (normalize3(ct) != 0) memcpy(ct, consensus, sizeof ct);
            memcpy(cs->tangent, ct, sizeof cs->tangent);
            member_offset += (size_t)count;
            cross_offset++;
        }
    }
    return member_offset == nmember && cross_offset == ncross ? 0 : -1;
}

static int final_solve(Arena_T arena, CoarseStrip *c,
                       const MarbleStripUvOptions *o,
                       AtlasStripOptions *so,
                       AtlasStripMetrics *out_metrics,
                       int *out_solves, double *out_last_change)
{
    int *window_start = NULL, nwindow = 0;
    size_t members_per_column = 0, nmember, ncross;
    if (build_cross_windows(c->H, &window_start, &nwindow) != 0) return -1;
    for (int ww = 0; ww < nwindow; ww++) {
        int count = c->H - window_start[ww];
        if (count > FINAL_CROSS_WINDOW_ROWS) count = FINAL_CROSS_WINDOW_ROWS;
        members_per_column += (size_t)count;
    }
    if ((size_t)c->W > SIZE_MAX / (size_t)nwindow ||
        (size_t)c->W > SIZE_MAX / members_per_column) {
        free(window_start); return -1;
    }
    ncross = (size_t)c->W * (size_t)nwindow;
    nmember = (size_t)c->W * members_per_column;
    if (c->nsample + ncross > c->x_capacity || ncross > INT32_MAX ||
        nmember > INT32_MAX) {
        free(window_start); return -1;
    }
    AtlasStripMember *member = (AtlasStripMember *)malloc(nmember * sizeof *member);
    AtlasStripCrossSection *cross = (AtlasStripCrossSection *)malloc(
        ncross * sizeof *cross);
    double *previous = (double *)malloc(c->nsample * sizeof *previous);
    if (member == NULL || cross == NULL || previous == NULL) {
        free(window_start); free(member); free(cross); free(previous); return -1;
    }
    int solves = 0;
    double last_change = DBL_MAX;
    AtlasStripProblem p;
    problem_from_coarse(c, &p);
    so->mode = ATLAS_STRIP_FINAL;
    for (int iteration = 0; iteration < o->final_iterations; iteration++) {
        if (rebuild_isolines(c, o, window_start, nwindow,
                             member, nmember, cross, ncross) != 0) goto fail;
        p.members = member; p.nmembers = nmember;
        p.cross_sections = cross; p.ncross_sections = ncross;
        memcpy(previous, c->x, c->nsample * sizeof *previous);
        AtlasStrip_initialize_intercepts(&p, so, c->x);
        Arena_Mark mark = Arena_save(arena);
        AtlasStripSystem system;
        double residual = 0.0;
        if (AtlasStrip_build(arena, &p, so, &system) != 0 ||
            MonotoneQp_solve_ls(arena, &system.qp, c->x, &residual) != 0) {
            Arena_restore(arena, mark); goto fail;
        }
        AtlasStrip_project_monotone(arena, &p, so, c->x);
        Arena_restore(arena, mark);
        double raw_change = 0.0;
        last_change = 0.0;
        for (size_t i = 0; i < c->nsample; i++) {
            double delta = c->x[i] - previous[i];
            if (fabs(delta) > raw_change) raw_change = fabs(delta);
            c->x[i] = previous[i] + o->final_damping * delta;
            if (fabs(o->final_damping * delta) > last_change)
                last_change = fabs(o->final_damping * delta);
        }
        solves++;
        if (o->verbose)
            fprintf(stderr,
                    "  [marble-strip final] iteration=%d raw/accepted-dU=%.9g/%.9g\n",
                    iteration + 1, raw_change, last_change);
        if (last_change <= o->movement_tolerance) break;
    }
    if (rebuild_isolines(c, o, window_start, nwindow,
                         member, nmember, cross, ncross) != 0) goto fail;
    p.members = member; p.nmembers = nmember;
    p.cross_sections = cross; p.ncross_sections = ncross;
    AtlasStrip_initialize_intercepts(&p, so, c->x);
    AtlasStrip_measure(&p, so, c->x, out_metrics);
    *out_solves = solves; *out_last_change = last_change;
    free(window_start); free(member); free(cross); free(previous); return 0;
fail:
    free(window_start); free(member); free(cross); free(previous); return -1;
}

static void prolongate_u(const CoarseStrip *c, const double *reference_uv,
                          const uint8_t *pin,
                          int H, int W, double *candidate)
{
    for (int r = 0; r < H; r++) {
        int ri = interval_for(c->row, c->H, r);
        double tr = (double)(r - c->row[ri]) /
                    (double)(c->row[ri + 1] - c->row[ri]);
        for (int col = 0; col < W; col++) {
            int ci = interval_for(c->col, c->W, col);
            double tc = (double)(col - c->col[ci]) /
                        (double)(c->col[ci + 1] - c->col[ci]);
            size_t a = (size_t)ri * (size_t)c->W + (size_t)ci;
            size_t b = a + 1;
            size_t q = (size_t)(ri + 1) * (size_t)c->W + (size_t)ci;
            size_t d = q + 1;
            /* Multigrid correction transfer: prolongate the solved *change*
             * (c->x - prior_u), not absolute U, so trusted fine-scale UV from
             * an earlier pass survives.  The correction goes to zero on its own
             * at the strip edges and the pinned halo, which the coarse solve
             * anchors to the reference U -- so NO per-vertex trust fade is
             * applied.  A fade scales a large interior correction non-uniformly
             * and manufactures a U fold that the monotonicity guard rejects. */
            double da = c->x[a] - c->prior_u[a];
            double db = c->x[b] - c->prior_u[b];
            double dq = c->x[q] - c->prior_u[q];
            double dd = c->x[d] - c->prior_u[d];
            double top = (1.0 - tc) * da + tc * db;
            double bottom = (1.0 - tc) * dq + tc * dd;
            size_t i = (size_t)r * (size_t)W + (size_t)col;
            double delta = (1.0 - tr) * top + tr * bottom;
            candidate[i] = reference_uv[i * 2] + delta;
            if (pin != NULL && pin[i]) candidate[i] = reference_uv[i * 2];
        }
    }
}

static double safe_commit_alpha(const double *reference_uv,
                                 const double *candidate,
                                 int H, int W, double fraction,
                                 double max_displacement,
                                 double *out_maxmove,
                                 double *out_worst_ratio,
                                 int *out_worst_row,
                                 int *out_worst_col)
{
    /* Largest uniform alpha in [0,1] such that committing du = alpha*(candidate
     * - reference) keeps every fine U edge monotone at >= the required minimum
     * speed and no vertex past the displacement trust radius.  A single global
     * scale preserves monotonicity where a per-vertex taper cannot: the
     * committed gap is a convex blend (1-alpha)*before + alpha*after that stays
     * >= req*before for alpha <= the per-edge bound.  Returning the graded alpha
     * (instead of 0/1) commits the largest safe fraction of the correction
     * rather than discarding a 96%-good solve over one slightly tight edge. */
    double maxmove = 0.0;
    double worst_ratio = DBL_MAX;
    int worst_row = -1, worst_col = -1;
    double req = fmin(0.5, fmax(1e-4, fraction));
    double alpha = 1.0;
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t i = (size_t)r * (size_t)W + (size_t)c;
        double move = fabs(candidate[i] - reference_uv[i * 2]);
        if (move > maxmove) maxmove = move;
        if (max_displacement > 0.0 && move > max_displacement) {
            double a = max_displacement / move;
            if (a < alpha) alpha = a;
        }
        if (c + 1 < W) {
            size_t j = i + 1;
            double before = reference_uv[j * 2] - reference_uv[i * 2];
            double after = candidate[j] - candidate[i];
            double ratio = before > 1e-9 && isfinite(after)
                         ? after / before : -DBL_MAX;
            if (ratio < worst_ratio) {
                worst_ratio = ratio; worst_row = r; worst_col = c;
            }
            if (!isfinite(after)) {
                alpha = 0.0;
            } else if (before > 1e-9) {
                if (after < before) {
                    double a = (1.0 - req) * before / (before - after);
                    if (a < alpha) alpha = a;
                }
            } else if (after < 0.0) {
                /* Degenerate reference edge: forbid only a backward step. */
                alpha = 0.0;
            }
        }
    }
    if (out_maxmove != NULL) *out_maxmove = maxmove;
    if (out_worst_ratio != NULL) *out_worst_ratio = worst_ratio;
    if (out_worst_row != NULL) *out_worst_row = worst_row;
    if (out_worst_col != NULL) *out_worst_col = worst_col;
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;
    return alpha;
}

int MarbleStripUv_solve(Arena_T arena, const float *verts,
                        const double *reference_uv, const double *trust,
                        const uint8_t *pin, int H, int W,
                        const MarbleStripUvOptions *input_options,
                        double *output_uv, MarbleStripUvStats *stats)
{
    MarbleStripUvOptions defaults;
    MarbleStripUvOptions_default(&defaults);
    const MarbleStripUvOptions *o = input_options != NULL ? input_options : &defaults;
    if (arena == NULL || verts == NULL || reference_uv == NULL || trust == NULL ||
        output_uv == NULL || H < 2 || W < 2 || o->sample_stride_u < 1 ||
        o->sample_stride_v < 1 || o->l1_iterations < 0 ||
        o->max_irls_iterations < 1 + o->l1_iterations ||
        o->max_irls_iterations > 1000 ||
        o->final_iterations < 0 || o->final_iterations > 1000 ||
        !(o->final_damping > 0.0) || o->final_damping > 1.0 ||
        !(o->movement_tolerance > 0.0) || o->lambda_length < 0.0 ||
        !(o->lambda_align > 0.0) || o->lambda_local < 0.0 ||
        o->lambda_prior < 0.0 || o->monotone_fraction < 0.0 ||
        o->membership_floor <= 0.0 || o->membership_floor > 1.0 ||
        !(o->match_radius > 0.0) || o->match_angle_deg <= 0.0 ||
        o->match_angle_deg >= 90.0 || o->match_search_columns < 1 ||
        o->match_search_columns > 1024 || o->topology_fallback <= 0.0 ||
        o->topology_fallback > 1.0 ||
        o->length_winsor < 1.05 || o->max_displacement < 0.0)
        return -1;
    MarbleStripUvStats local;
    memset(&local, 0, sizeof local);
    CoarseStrip c;
    if (coarse_build(verts, reference_uv, trust, pin, H, W, o, &c) != 0) {
        if (o->verbose)
            fprintf(stderr, "[marble_strip_uv] coarse correspondence build failed\n");
        return -1;
    }
    if (o->verbose)
        fprintf(stderr,
                "[marble_strip_uv] coarse=%dx%d samples=%zu cross=%zu "
                "members=%zu anchors=%zu\n",
                c.H, c.W, c.nsample, c.ncross, c.nmember, c.nanchor);
    local.coarse_h = c.H; local.coarse_w = c.W;
    local.coarse_samples = c.nsample;
    local.coarse_cross_sections = c.ncross;
    local.coarse_members = c.nmember;
    local.exact_anchors = c.nanchor;

    AtlasStripProblem p;
    problem_from_coarse(&c, &p);
    AtlasStripOptions so;
    AtlasStripOptions_default(&so);
    so.mode = ATLAS_STRIP_RELAXED;
    so.lambda_length = o->lambda_length;
    so.lambda_align = o->lambda_align;
    so.lambda_local = o->lambda_local;
    so.lambda_prior = o->lambda_prior;
    so.monotone_fraction = o->monotone_fraction;
    so.membership_floor = 1e-10;
    AtlasStrip_initialize_intercepts(&p, &so, c.x);
    AtlasStripMetrics before;
    AtlasStrip_measure(&p, &so, c.x, &before);

    AtlasStripRobustOptions ro;
    AtlasStripRobustOptions_default(&ro);
    ro.l1_iterations = o->l1_iterations;
    /* The public cap is a cap on global variational solves, not merely on the
     * last phase: initial + L1 + Gaussian rounds never exceed it. */
    ro.likelihood_iterations = o->max_irls_iterations -
                               o->l1_iterations - 1;
    ro.likelihood_sigma = o->likelihood_sigma;
    ro.convergence_tolerance = o->movement_tolerance;
    ro.use_active_set = 0;
    MonotoneQpOptions qo;
    MonotoneQpOptions_default(&qo);
    qo.verbose = o->verbose;
    AtlasStripRobustStats robust;
    memset(&robust, 0, sizeof robust);
    if (AtlasStrip_solve_robust(arena, &p, &so, &ro, &qo, c.x,
                               c.membership, NULL, &robust) != 0) {
        if (o->verbose)
            fprintf(stderr, "[marble_strip_uv] robust variational solve failed\n");
        coarse_dispose(&c); return -1;
    }
    local.robust_solves = robust.total_qp_solves;
    local.likelihood_sigma = robust.likelihood_sigma;
    local.downweighted_members = robust.downweighted_members;
    local.last_max_u_change = robust.max_u_change;

    local.membership_min = DBL_MAX;
    double *like_sum = (double *)calloc(c.nsample, sizeof *like_sum);
    double *like_weight = (double *)calloc(c.nsample, sizeof *like_weight);
    if (like_sum == NULL || like_weight == NULL) {
        free(like_sum); free(like_weight); coarse_dispose(&c); return -1;
    }
    for (size_t mi = 0; mi < c.nmember; mi++) {
        double prior = c.member[mi].membership;
        double like = prior > 0.0 ? c.membership[mi] / prior : 0.0;
        like = clamp01(like);
        if (like < local.membership_min) local.membership_min = like;
        if (like > local.membership_max) local.membership_max = like;
        local.membership_mean += like;
        int32_t a = c.member[mi].value0, b = c.member[mi].value1;
        if (a >= 0 && (size_t)a < c.nsample) {
            double w = b >= 0 ? 1.0 - c.member[mi].value_t : 1.0;
            like_sum[a] += w * like; like_weight[a] += w;
        }
        if (b >= 0 && (size_t)b < c.nsample) {
            double w = c.member[mi].value_t;
            like_sum[b] += w * like; like_weight[b] += w;
        }
    }
    local.membership_mean /= (double)c.nmember;
    if (local.membership_min == DBL_MAX) local.membership_min = 0.0;
    for (size_t si = 0; si < c.nsample; si++)
        c.sample_likelihood[si] = like_weight[si] > 0.0
                                ? like_sum[si] / like_weight[si] : 1.0;
    free(like_sum); free(like_weight);
    for (size_t mi = 0; mi < c.nmember; mi++)
        c.member[mi].membership = c.membership[mi];

    AtlasStripMetrics after;
    int final_solves = 0;
    double final_change = robust.max_u_change;
    if (o->final_iterations > 0) {
        if (final_solve(arena, &c, o, &so, &after,
                        &final_solves, &final_change) != 0) {
            coarse_dispose(&c); return -1;
        }
        local.converged = final_change <= o->movement_tolerance;
    } else {
        so.mode = ATLAS_STRIP_RELAXED;
        AtlasStrip_initialize_intercepts(&p, &so, c.x);
        AtlasStrip_measure(&p, &so, c.x, &after);
        local.converged = robust.max_u_change <= o->movement_tolerance &&
                          robust.max_membership_change <= o->movement_tolerance;
    }
    local.final_solves = final_solves;
    local.last_max_u_change = final_change;
    local.length_rms_before = before.rms_length_row;
    local.length_rms_after = after.rms_length_row;
    local.align_rms_before = before.rms_align_row;
    local.align_rms_after = after.rms_align_row;
    local.local_speed_rms_before = before.rms_local_speed_error;
    local.local_speed_rms_after = after.rms_local_speed_error;
    local.min_monotone_ratio = after.min_monotone_ratio;

    size_t n = (size_t)H * (size_t)W;
    double *candidate = (double *)malloc(n * sizeof *candidate);
    if (candidate == NULL) { coarse_dispose(&c); return -1; }
    prolongate_u(&c, reference_uv, pin, H, W, candidate);
    double candidate_maxmove = 0.0, candidate_worst_ratio = DBL_MAX;
    int candidate_worst_row = -1, candidate_worst_col = -1;
    double alpha = safe_commit_alpha(reference_uv, candidate, H, W,
                                      o->monotone_fraction,
                                      o->max_displacement,
                                      &candidate_maxmove,
                                      &candidate_worst_ratio,
                                      &candidate_worst_row,
                                      &candidate_worst_col);
    if (!(alpha > 0.0)) {
        if (o->verbose)
            fprintf(stderr,
                    "[marble_strip_uv] full candidate rejected by monotonicity/"
                    "displacement guard: maxmove=%.9g limit=%.6g "
                    "worst-ratio=%.9g at (%d,%d)\n",
                    candidate_maxmove, o->max_displacement,
                    candidate_worst_ratio,
                    candidate_worst_row, candidate_worst_col);
        free(candidate); coarse_dispose(&c); return -1;
    }
    local.safe_commit_alpha = alpha;
    if (o->verbose)
        fprintf(stderr,
                "[marble_strip_uv] committed alpha=%.6g maxmove=%.6g "
                "pre-scale worst-ratio=%.6g at (%d,%d)\n",
                alpha, candidate_maxmove, candidate_worst_ratio,
                candidate_worst_row, candidate_worst_col);
    double sumsq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double du = alpha * (candidate[i] - reference_uv[i * 2]);
        output_uv[i * 2] = reference_uv[i * 2] + du;
        output_uv[i * 2 + 1] = reference_uv[i * 2 + 1];
        if (fabs(du) > 1e-12) local.moved_vertices++;
        sumsq += du * du;
        if (fabs(du) > local.max_u_displacement)
            local.max_u_displacement = fabs(du);
    }
    local.rms_u_displacement = sqrt(sumsq / (double)n);
    free(candidate); coarse_dispose(&c);
    if (stats != NULL) *stats = local;
    return 0;
}

int MarbleStripUv_selftest(void)
{
    enum { H = 17, W = 33 };
    size_t n = (size_t)H * (size_t)W;
    float *verts = (float *)malloc(n * 3 * sizeof *verts);
    double *uv = (double *)malloc(n * 2 * sizeof *uv);
    double *trust = (double *)malloc(n * sizeof *trust);
    uint8_t *pin = (uint8_t *)calloc(n, 1);
    double *out = (double *)malloc(n * 2 * sizeof *out);
    if (verts == NULL || uv == NULL || trust == NULL || pin == NULL || out == NULL) {
        free(verts); free(uv); free(trust); free(pin); free(out); return 1;
    }
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t i = (size_t)r * W + (size_t)c;
        verts[i * 3] = (float)r;
        verts[i * 3 + 1] = 0.0f;
        verts[i * 3 + 2] = (float)(2 * c);
        uv[i * 2] = (double)(2 * c); uv[i * 2 + 1] = (double)r;
        trust[i] = 1.0;
        if (c == 0 || c + 1 == W) pin[i] = 1;
        if (r >= 5 && r <= 11 && c >= 9 && c <= 23) {
            double phase = 3.14159265358979323846 * (double)(c - 9) / 14.0;
            verts[i * 3 + 1] = (float)(10.0 * sin(phase));
            trust[i] = 0.0;
        }
    }
    Arena_T arena = Arena_new();
    MarbleStripUvOptions o; MarbleStripUvOptions_default(&o);
    o.sample_stride_u = 2; o.sample_stride_v = 2;
    o.l1_iterations = 2; o.max_irls_iterations = 20; o.final_iterations = 0;
    o.lambda_prior = 30.0; o.max_displacement = 32.0;
    MarbleStripUvStats s;
    int fails = 0;
    if (arena == NULL || MarbleStripUv_solve(arena, verts, uv, trust, pin,
            H, W, &o, out, &s) != 0) {
        fprintf(stderr, "[marble_strip_uv selftest] FAIL solve\n"); fails++;
    } else {
        double min_du = DBL_MAX, max_v_error = 0.0;
        for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
            size_t i = (size_t)r * W + (size_t)c;
            if (c + 1 < W) {
                double d = out[(i + 1) * 2] - out[i * 2];
                if (d < min_du) min_du = d;
            }
            double e = fabs(out[i * 2 + 1] - uv[i * 2 + 1]);
            if (e > max_v_error) max_v_error = e;
        }
        if (!(min_du > 0.0)) {
            fprintf(stderr, "[marble_strip_uv selftest] FAIL monotone %.6g\n", min_du);
            fails++;
        }
        if (max_v_error != 0.0) {
            fprintf(stderr, "[marble_strip_uv selftest] FAIL V changed %.6g\n", max_v_error);
            fails++;
        }
        if (s.moved_vertices == 0 || !(s.safe_commit_alpha > 0.0)) {
            fprintf(stderr, "[marble_strip_uv selftest] FAIL no variational move\n");
            fails++;
        }
    }
    Arena_dispose(&arena);
    free(verts); free(uv); free(trust); free(pin); free(out);
    fprintf(stderr, "[marble_strip_uv selftest] %s\n", fails ? "FAILED" : "ok");
    return fails;
}
