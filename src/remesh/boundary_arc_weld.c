/*
 * boundary_arc_weld.c -- topology/provenance-aware exact boundary zipper.
 *
 * A connected component is not a layer label on a folded scroll.  Consequently
 * the late exact-position cleanup must not merge every coincident pair in one
 * component.  This module first recovers the directed manifold boundary, then
 * accepts only reciprocal pairs whose face normals agree and whose directed
 * boundary tangents oppose.  Isolated contacts are rejected unless they close
 * an existing one-vertex link.  Non-local pairs need repeated, monotone support
 * along both boundary chains.  The final Weld_verts_filtered call merely
 * executes this precomputed matching; it makes no geometric decisions itself.
 */
#include "boundary_arc_weld.h"

#include "../common/vert_weld.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int32_t lo, hi;
    int32_t va, vb;
} BawHalfEdge;

typedef struct {
    int64_t z, y, x;
    int32_t v;
} BawCell;

typedef struct {
    int32_t a, b;
    uint8_t local_link;
} BawPair;

typedef struct {
    int32_t a, b;
    double length;
    uint8_t edge_adjacent;
} BawArcLink;

typedef struct {
    int32_t root, chart_a, chart_b;
    size_t nodes;
    double support;
} BawChartCandidate;

typedef struct {
    size_t nv;
    int32_t *allowed_partner;
    BoundaryArcWeldStats *stats;
} BawFilter;

static int baw_halfedge_cmp(const void *pa, const void *pb)
{
    const BawHalfEdge *a = (const BawHalfEdge *)pa;
    const BawHalfEdge *b = (const BawHalfEdge *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}

static int baw_cell_cmp(const void *pa, const void *pb)
{
    const BawCell *a = (const BawCell *)pa;
    const BawCell *b = (const BawCell *)pb;
    if (a->z != b->z) return a->z < b->z ? -1 : 1;
    if (a->y != b->y) return a->y < b->y ? -1 : 1;
    if (a->x != b->x) return a->x < b->x ? -1 : 1;
    if (a->v != b->v) return a->v < b->v ? -1 : 1;
    return 0;
}

static int baw_link_cmp(const void *pa, const void *pb)
{
    const BawArcLink *a = (const BawArcLink *)pa;
    const BawArcLink *b = (const BawArcLink *)pb;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    return 0;
}

static int baw_chart_candidate_cmp(const void *pa, const void *pb)
{
    const BawChartCandidate *a = (const BawChartCandidate *)pa;
    const BawChartCandidate *b = (const BawChartCandidate *)pb;
    if (a->support != b->support) return a->support > b->support ? -1 : 1;
    if (a->nodes != b->nodes) return a->nodes > b->nodes ? -1 : 1;
    return a->root < b->root ? -1 : (a->root > b->root ? 1 : 0);
}

static int baw_cell_coord_cmp(const BawCell *a, int64_t z, int64_t y, int64_t x)
{
    if (a->z != z) return a->z < z ? -1 : 1;
    if (a->y != y) return a->y < y ? -1 : 1;
    if (a->x != x) return a->x < x ? -1 : 1;
    return 0;
}

static size_t baw_cell_lower(const BawCell *cells, size_t n,
                             int64_t z, int64_t y, int64_t x)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t m = lo + (hi - lo) / 2;
        if (baw_cell_coord_cmp(&cells[m], z, y, x) < 0) lo = m + 1;
        else hi = m;
    }
    return lo;
}

static double baw_dist(const float *verts, int32_t a, int32_t b)
{
    double dz = (double)verts[(size_t)a*3+0] - verts[(size_t)b*3+0];
    double dy = (double)verts[(size_t)a*3+1] - verts[(size_t)b*3+1];
    double dx = (double)verts[(size_t)a*3+2] - verts[(size_t)b*3+2];
    return sqrt(dz*dz + dy*dy + dx*dx);
}

static double baw_dot3(const float *a, const float *b)
{
    return (double)a[0]*b[0] + (double)a[1]*b[1] + (double)a[2]*b[2];
}

static int32_t baw_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void baw_union(int32_t *parent, int32_t a, int32_t b)
{
    a = baw_find(parent, a); b = baw_find(parent, b);
    if (a == b) return;
    if (a < b) parent[b] = a;
    else parent[a] = b;
}

void BoundaryArcWeld_default_params(BoundaryArcWeldParams *p)
{
    if (!p) return;
    p->normal_dot_min = 0.50;
    p->tangent_dot_max = -0.50;
    p->ambiguity_margin = 0.10;
    p->max_anchor_gap = 32.0;
    p->min_support_length = 12.0;
    p->min_anchors = 2;
}

/* Walk one directed boundary chain to its first other reciprocal anchor. */
static int baw_next_anchor(const float *verts,
                           const int32_t *prev, const int32_t *next,
                           const int32_t *pair_id,
                           int32_t start, int forward, int32_t self_pair,
                           double max_length,
                           int32_t *out_pair, int32_t *out_vert,
                           double *out_length, int *out_steps)
{
    int32_t v = start;
    double length = 0.0;
    /* The length bound is the real termination condition.  The step cap is
     * only corruption insurance for zero-length/duplicate boundary edges. */
    for (int step = 0; step < 512; step++) {
        int32_t w = forward ? next[v] : prev[v];
        if (w < 0 || w == v) return 0;
        length += baw_dist(verts, v, w);
        if (length > max_length) return 0;
        v = w;
        if (pair_id[v] >= 0 && pair_id[v] != self_pair) {
            *out_pair = pair_id[v];
            *out_vert = v;
            *out_length = length;
            if (out_steps) *out_steps = step + 1;
            return 1;
        }
    }
    return 0;
}

static bool baw_pair_filter(size_t i, size_t j, void *context)
{
    BawFilter *f = (BawFilter *)context;
    if (f->stats) f->stats->filter_tests++;
    if (i >= f->nv || j >= f->nv) return false;
    if (f->allowed_partner[i] != (int32_t)j ||
        f->allowed_partner[j] != (int32_t)i) return false;
    if (f->stats) f->stats->filter_accepts++;
    return true;
}

int BoundaryArcWeld_process_with_policy(Arena_T arena,
                             const float *verts, size_t nv,
                            int32_t *faces, size_t nf,
                            float eps,
                             const int32_t *merge_group,
                             const BoundaryArcWeldParams *params_in,
                             const BoundaryArcWeldPolicy *policy,
                             float **out_verts, size_t *out_nv,
                            size_t *out_nf,
                            BoundaryArcWeldStats *stats)
{
    BoundaryArcWeldParams defaults;
    const BoundaryArcWeldParams *p = params_in;
    BawHalfEdge *he = NULL;
    int32_t *prev = NULL, *next = NULL;
    uint8_t *bdeg = NULL;
    float *normal = NULL, *tangent = NULL;
    BawCell *cells = NULL;
    int32_t *best = NULL, *pair_id = NULL, *partner = NULL;
    int32_t *allowed_partner = NULL, *pair_parent = NULL;
    double *best_score = NULL, *second_score = NULL;
    BawPair *pairs = NULL;
    BawArcLink *links = NULL;
    size_t *cluster_nodes = NULL;
    double *cluster_support = NULL;
    size_t *cluster_links = NULL;
    uint32_t *pair_degree = NULL;
    uint8_t *cluster_edge_zipper = NULL, *cluster_accept = NULL;
    int32_t *cluster_chart_a = NULL, *cluster_chart_b = NULL;
    uint8_t *cluster_chart_bad = NULL, *cluster_self = NULL;
    uint8_t *cluster_nondisk = NULL;
    BawChartCandidate *chart_candidate = NULL;
    int32_t *chart_parent = NULL;
    size_t n_cells = 0, n_pairs = 0, n_links = 0;
    int rc = -1;

    if (!arena || !verts || !faces || !out_verts || !out_nv || !out_nf ||
        eps <= 0.0f) return -1;
    if (!p) { BoundaryArcWeld_default_params(&defaults); p = &defaults; }
    if (stats) memset(stats, 0, sizeof(*stats));
    if (nv == 0 || nf == 0) {
        Weld_verts_filtered(arena, verts, nv, NULL, faces, nf, out_nf,
                            eps, true, merge_group, NULL, NULL,
                            out_verts, out_nv, NULL);
        return 0;
    }

    if (nv > (size_t)INT32_MAX || nf > SIZE_MAX/3) goto cleanup;

    /* Directed run-1 halfedges recover both the boundary topology and its
     * orientation.  Opposite sides of one valid seam traverse in opposite
     * directions when the incident faces have compatible winding. */
    {
        size_t hn = nf * 3;
        he = (BawHalfEdge *)malloc(hn * sizeof(*he));
        prev = (int32_t *)malloc(nv * sizeof(*prev));
        next = (int32_t *)malloc(nv * sizeof(*next));
        bdeg = (uint8_t *)calloc(nv, 1);
        normal = (float *)calloc(nv * 3, sizeof(*normal));
        tangent = (float *)calloc(nv * 3, sizeof(*tangent));
        if (!he || !prev || !next || !bdeg || !normal || !tangent) goto cleanup;
        for (size_t v = 0; v < nv; v++) prev[v] = next[v] = -1;

        for (size_t f = 0; f < nf; f++) {
            int32_t a = faces[f*3+0], b = faces[f*3+1], c = faces[f*3+2];
            if (a < 0 || b < 0 || c < 0 ||
                (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv)
                goto cleanup;
            int32_t tri[3] = {a,b,c};
            for (int e = 0; e < 3; e++) {
                int32_t x = tri[e], y = tri[(e+1)%3];
                BawHalfEdge *h = &he[f*3+(size_t)e];
                h->lo = x < y ? x : y; h->hi = x < y ? y : x;
                h->va = x; h->vb = y;
            }
            const float *pa = &verts[(size_t)a*3];
            const float *pb = &verts[(size_t)b*3];
            const float *pc = &verts[(size_t)c*3];
            double u0=pb[0]-pa[0], u1=pb[1]-pa[1], u2=pb[2]-pa[2];
            double v0=pc[0]-pa[0], v1=pc[1]-pa[1], v2=pc[2]-pa[2];
            float n0=(float)(u1*v2-u2*v1);
            float n1=(float)(u2*v0-u0*v2);
            float n2=(float)(u0*v1-u1*v0);
            int32_t ids[3] = {a,b,c};
            for (int k = 0; k < 3; k++) {
                normal[(size_t)ids[k]*3+0] += n0;
                normal[(size_t)ids[k]*3+1] += n1;
                normal[(size_t)ids[k]*3+2] += n2;
            }
        }
        qsort(he, hn, sizeof(*he), baw_halfedge_cmp);
        for (size_t i = 0; i < hn; ) {
            size_t j = i + 1;
            while (j < hn && he[j].lo == he[i].lo && he[j].hi == he[i].hi) j++;
            if (j - i == 1) {
                int32_t a = he[i].va, b = he[i].vb;
                if (bdeg[a] < UINT8_MAX) bdeg[a]++;
                if (bdeg[b] < UINT8_MAX) bdeg[b]++;
                if (next[a] == -1) next[a] = b; else next[a] = -2;
                if (prev[b] == -1) prev[b] = a; else prev[b] = -2;
            }
            i = j;
        }
        for (size_t v = 0; v < nv; v++) {
            double nl = sqrt(baw_dot3(&normal[v*3], &normal[v*3]));
            if (nl > 1e-20) {
                normal[v*3+0] = (float)(normal[v*3+0]/nl);
                normal[v*3+1] = (float)(normal[v*3+1]/nl);
                normal[v*3+2] = (float)(normal[v*3+2]/nl);
            }
            if (bdeg[v] == 2 && prev[v] >= 0 && next[v] >= 0 &&
                prev[v] != next[v]) {
                double t0 = (double)verts[(size_t)next[v]*3+0] -
                            verts[(size_t)prev[v]*3+0];
                double t1 = (double)verts[(size_t)next[v]*3+1] -
                            verts[(size_t)prev[v]*3+1];
                double t2 = (double)verts[(size_t)next[v]*3+2] -
                            verts[(size_t)prev[v]*3+2];
                double tl = sqrt(t0*t0+t1*t1+t2*t2);
                if (tl > 1e-20 && nl > 1e-20) {
                    tangent[v*3+0] = (float)(t0/tl);
                    tangent[v*3+1] = (float)(t1/tl);
                    tangent[v*3+2] = (float)(t2/tl);
                    n_cells++;
                }
            }
        }
    }
    if (stats) stats->boundary_vertices = n_cells;

    cells = (BawCell *)malloc((n_cells ? n_cells : 1) * sizeof(*cells));
    best = (int32_t *)malloc(nv * sizeof(*best));
    best_score = (double *)malloc(nv * sizeof(*best_score));
    second_score = (double *)malloc(nv * sizeof(*second_score));
    if (!cells || !best || !best_score || !second_score) goto cleanup;
    {
        double inv = 1.0 / (double)eps;
        size_t k = 0;
        for (size_t v = 0; v < nv; v++) {
            best[v] = -1; best_score[v] = -DBL_MAX; second_score[v] = -DBL_MAX;
            if (bdeg[v] != 2 || prev[v] < 0 || next[v] < 0 ||
                baw_dot3(&tangent[v*3], &tangent[v*3]) < 0.5) continue;
            cells[k].z = (int64_t)floor((double)verts[v*3+0] * inv);
            cells[k].y = (int64_t)floor((double)verts[v*3+1] * inv);
            cells[k].x = (int64_t)floor((double)verts[v*3+2] * inv);
            cells[k].v = (int32_t)v;
            k++;
        }
        qsort(cells, n_cells, sizeof(*cells), baw_cell_cmp);

        double eps2 = (double)eps * eps;
        for (size_t ci = 0; ci < n_cells; ci++) {
            int32_t v = cells[ci].v;
            for (int dz = -1; dz <= 1; dz++)
            for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                int64_t z = cells[ci].z + dz;
                int64_t y = cells[ci].y + dy;
                int64_t x = cells[ci].x + dx;
                size_t q = baw_cell_lower(cells, n_cells, z, y, x);
                for (; q < n_cells &&
                       baw_cell_coord_cmp(&cells[q], z, y, x) == 0; q++) {
                    int32_t w = cells[q].v;
                    if (w == v) continue;
                    if (merge_group && merge_group[v] != merge_group[w]) continue;
                    double d0=(double)verts[(size_t)w*3+0]-verts[(size_t)v*3+0];
                    double d1=(double)verts[(size_t)w*3+1]-verts[(size_t)v*3+1];
                    double d2=(double)verts[(size_t)w*3+2]-verts[(size_t)v*3+2];
                    double dd=d0*d0+d1*d1+d2*d2;
                    if (dd > eps2) continue;
                    if (stats && v < w) stats->spatial_candidates++;
                    double nd = baw_dot3(&normal[(size_t)v*3],
                                         &normal[(size_t)w*3]);
                    double td = baw_dot3(&tangent[(size_t)v*3],
                                         &tangent[(size_t)w*3]);
                    if (nd < p->normal_dot_min || td > p->tangent_dot_max)
                        continue;
                    /* Differential compatibility dominates; distance only
                     * breaks sub-epsilon ties. */
                    double score = nd - td - 0.25*sqrt(dd)/(double)eps;
                    if (score > best_score[v]) {
                        second_score[v] = best_score[v];
                        best_score[v] = score; best[v] = w;
                    } else if (score > second_score[v]) {
                        second_score[v] = score;
                    }
                }
            }
        }
    }

    /* A local nearest neighbour is evidence only when it is unambiguous in
     * both directions.  This is the one-to-one part of the correspondence. */
    for (size_t v = 0; v < nv; v++) {
        if (best[v] >= 0 && second_score[v] > -DBL_MAX/2 &&
            best_score[v] - second_score[v] < p->ambiguity_margin) {
            best[v] = -1;
            if (stats) stats->ambiguous_vertices++;
        }
    }
    pairs = (BawPair *)malloc((n_cells/2 + 1) * sizeof(*pairs));
    pair_id = (int32_t *)malloc(nv * sizeof(*pair_id));
    partner = (int32_t *)malloc(nv * sizeof(*partner));
    allowed_partner = (int32_t *)malloc(nv * sizeof(*allowed_partner));
    if (!pairs || !pair_id || !partner || !allowed_partner) goto cleanup;
    for (size_t v = 0; v < nv; v++) {
        pair_id[v] = partner[v] = allowed_partner[v] = -1;
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t w = best[v];
        if (w < 0 || (size_t)w <= v || best[w] != (int32_t)v) continue;
        BawPair *pair = &pairs[n_pairs];
        pair->a = (int32_t)v; pair->b = w;
        pair->local_link =
            (next[v] == prev[w] || prev[v] == next[w]) ? 1u : 0u;
        pair_id[v] = pair_id[w] = (int32_t)n_pairs;
        partner[v] = w; partner[w] = (int32_t)v;
        n_pairs++;
    }
    if (stats) stats->reciprocal_pairs = n_pairs;

    pair_parent = (int32_t *)malloc((n_pairs ? n_pairs : 1) * sizeof(*pair_parent));
    links = (BawArcLink *)malloc((n_pairs*2 + 1) * sizeof(*links));
    cluster_nodes = (size_t *)calloc((n_pairs ? n_pairs : 1), sizeof(*cluster_nodes));
    cluster_support = (double *)calloc((n_pairs ? n_pairs : 1), sizeof(*cluster_support));
    if (!pair_parent || !links || !cluster_nodes || !cluster_support) goto cleanup;
    for (size_t i = 0; i < n_pairs; i++) pair_parent[i] = (int32_t)i;

    /* Monotone support: advancing along A must retreat along B.  The first
     * anchor in each direction must be the same reciprocal pair, and the two
     * arclengths must agree up to sampling jitter. */
    for (size_t i = 0; i < n_pairs; i++) {
        int32_t a = pairs[i].a, b = pairs[i].b;
        for (int side = 0; side < 2; side++) {
            int32_t pa=-1, pb=-1, va=-1, vb=-1;
            int sa=0, sb=0;
            double la=0.0, lb=0.0;
            int oka = baw_next_anchor(verts, prev, next, pair_id, a,
                                      side == 0, (int32_t)i,
                                      p->max_anchor_gap, &pa, &va, &la, &sa);
            int okb = baw_next_anchor(verts, prev, next, pair_id, b,
                                      side != 0, (int32_t)i,
                                      p->max_anchor_gap, &pb, &vb, &lb, &sb);
            if (!oka || !okb || pa != pb || pa == (int32_t)i ||
                partner[va] != vb) continue;
            double tol = 0.35 * (la > lb ? la : lb);
            if (tol < 3.0) tol = 3.0;
            if (fabs(la-lb) > tol) continue;
            if ((int32_t)i < pa) {
                links[n_links].a = (int32_t)i;
                links[n_links].b = pa;
                links[n_links].length = 0.5*(la+lb);
                links[n_links].edge_adjacent = (sa == 1 && sb == 1) ? 1u : 0u;
                n_links++;
            }
        }
    }
    /* The same neighbouring anchor is normally discovered from both ends.
     * Collapse that duplicate before measuring path topology/support. */
    if (n_links > 1) {
        size_t write = 0;
        qsort(links, n_links, sizeof(*links), baw_link_cmp);
        for (size_t i = 0; i < n_links; ) {
            size_t j = i + 1;
            BawArcLink link = links[i];
            while (j < n_links && links[j].a == link.a && links[j].b == link.b) {
                if (links[j].length < link.length) link.length = links[j].length;
                link.edge_adjacent =
                    (uint8_t)(link.edge_adjacent && links[j].edge_adjacent);
                j++;
            }
            links[write++] = link;
            i = j;
        }
        n_links = write;
    }
    cluster_links = (size_t *)calloc((n_pairs ? n_pairs : 1),
                                     sizeof(*cluster_links));
    pair_degree = (uint32_t *)calloc((n_pairs ? n_pairs : 1),
                                    sizeof(*pair_degree));
    cluster_edge_zipper = (uint8_t *)malloc(n_pairs ? n_pairs : 1);
    cluster_accept = (uint8_t *)calloc((n_pairs ? n_pairs : 1), 1);
    if (!cluster_links || !pair_degree || !cluster_edge_zipper ||
        !cluster_accept) goto cleanup;
    memset(cluster_edge_zipper, 1, n_pairs ? n_pairs : 1);
    for (size_t i = 0; i < n_links; i++)
        baw_union(pair_parent, links[i].a, links[i].b);
    for (size_t i = 0; i < n_pairs; i++) {
        int32_t r = baw_find(pair_parent, (int32_t)i);
        cluster_nodes[r]++;
    }
    for (size_t i = 0; i < n_links; i++) {
        int32_t r = baw_find(pair_parent, links[i].a);
        cluster_support[r] += links[i].length;
        cluster_links[r]++;
        pair_degree[links[i].a]++;
        pair_degree[links[i].b]++;
        if (!links[i].edge_adjacent) cluster_edge_zipper[r] = 0;
    }

    for (size_t i = 0; i < n_pairs; i++) {
        if (baw_find(pair_parent, (int32_t)i) != (int32_t)i) continue;
        if (stats) stats->arc_clusters++;
        if ((!policy || !policy->cross_chart_forest) &&
            cluster_nodes[i] >= p->min_anchors &&
            cluster_support[i] >= p->min_support_length) {
            if (stats) stats->accepted_clusters++;
        }
    }

    if (policy && policy->cross_chart_forest) {
        size_t ncandidate = 0;
        if (!policy->chart_component || !policy->chart_is_disk ||
            policy->n_chart_components == 0 ||
            policy->n_chart_components > (size_t)INT32_MAX)
            goto cleanup;
        cluster_chart_a = (int32_t *)malloc((n_pairs ? n_pairs : 1) * sizeof(int32_t));
        cluster_chart_b = (int32_t *)malloc((n_pairs ? n_pairs : 1) * sizeof(int32_t));
        cluster_chart_bad = (uint8_t *)calloc((n_pairs ? n_pairs : 1), 1);
        cluster_self = (uint8_t *)calloc((n_pairs ? n_pairs : 1), 1);
        cluster_nondisk = (uint8_t *)calloc((n_pairs ? n_pairs : 1), 1);
        chart_candidate = (BawChartCandidate *)malloc(
            (n_pairs ? n_pairs : 1) * sizeof(*chart_candidate));
        chart_parent = (int32_t *)malloc(policy->n_chart_components * sizeof(int32_t));
        if (!cluster_chart_a || !cluster_chart_b || !cluster_chart_bad ||
            !cluster_self || !cluster_nondisk || !chart_candidate ||
            !chart_parent) goto cleanup;
        for (size_t i = 0; i < n_pairs; i++)
            cluster_chart_a[i] = cluster_chart_b[i] = -1;
        for (size_t c = 0; c < policy->n_chart_components; c++)
            chart_parent[c] = (int32_t)c;
        for (size_t i = 0; i < n_pairs; i++) {
            int32_t r = baw_find(pair_parent, (int32_t)i);
            int32_t ca = policy->chart_component[pairs[i].a];
            int32_t cb = policy->chart_component[pairs[i].b];
            if (ca < 0 || cb < 0 ||
                (size_t)ca >= policy->n_chart_components ||
                (size_t)cb >= policy->n_chart_components) {
                cluster_chart_bad[r] = 1;
                continue;
            }
            if (ca == cb) cluster_self[r] = 1;
            if (!policy->chart_is_disk[ca] || !policy->chart_is_disk[cb])
                cluster_nondisk[r] = 1;
            if (ca > cb) { int32_t t = ca; ca = cb; cb = t; }
            if (cluster_chart_a[r] < 0) {
                cluster_chart_a[r] = ca;
                cluster_chart_b[r] = cb;
            } else if (cluster_chart_a[r] != ca || cluster_chart_b[r] != cb) {
                cluster_chart_bad[r] = 1;
            }
        }
        for (size_t i = 0; i < n_pairs; i++) {
            size_t endpoints = 0;
            int path_ok = 1;
            if (baw_find(pair_parent, (int32_t)i) != (int32_t)i) continue;
            if (cluster_nodes[i] < p->min_anchors ||
                cluster_support[i] < p->min_support_length)
                continue;
            if (stats) stats->chart_clusters_considered++;
            if (cluster_self[i]) {
                if (stats) stats->chart_clusters_rejected_self++;
                continue;
            }
            if (cluster_nondisk[i]) {
                if (stats) stats->chart_clusters_rejected_nondisk++;
                continue;
            }
            if (cluster_chart_bad[i] || cluster_chart_a[i] < 0) {
                if (stats) stats->chart_clusters_rejected_nonpath++;
                continue;
            }
            for (size_t j = 0; j < n_pairs; j++) {
                if (baw_find(pair_parent, (int32_t)j) != (int32_t)i) continue;
                if (pair_degree[j] == 1) endpoints++;
                else if (pair_degree[j] != 2) path_ok = 0;
            }
            if (!path_ok || endpoints != 2 ||
                cluster_links[i] + 1 != cluster_nodes[i]) {
                if (stats) stats->chart_clusters_rejected_nonpath++;
                continue;
            }
            if (policy->require_edge_zipper && !cluster_edge_zipper[i]) {
                if (stats) stats->chart_clusters_rejected_sparse++;
                continue;
            }
            chart_candidate[ncandidate].root = (int32_t)i;
            chart_candidate[ncandidate].chart_a = cluster_chart_a[i];
            chart_candidate[ncandidate].chart_b = cluster_chart_b[i];
            chart_candidate[ncandidate].nodes = cluster_nodes[i];
            chart_candidate[ncandidate].support = cluster_support[i];
            ncandidate++;
        }
        qsort(chart_candidate, ncandidate, sizeof(*chart_candidate),
              baw_chart_candidate_cmp);
        for (size_t i = 0; i < ncandidate; i++) {
            int32_t ca = baw_find(chart_parent, chart_candidate[i].chart_a);
            int32_t cb = baw_find(chart_parent, chart_candidate[i].chart_b);
            if (ca == cb) {
                if (stats) stats->chart_clusters_rejected_cycle++;
                continue;
            }
            baw_union(chart_parent, ca, cb);
            cluster_accept[chart_candidate[i].root] = 1;
            if (stats) {
                stats->chart_clusters_accepted++;
                stats->accepted_clusters++;
            }
        }
    }
    for (size_t i = 0; i < n_pairs; i++) {
        int32_t r = baw_find(pair_parent, (int32_t)i);
        int accept_arc = cluster_nodes[r] >= p->min_anchors &&
                         cluster_support[r] >= p->min_support_length;
        if (pairs[i].local_link && stats) stats->local_link_pairs++;
        if (policy && policy->cross_chart_forest)
            accept_arc = cluster_accept[r] ? 1 : 0;
        if ((policy && policy->cross_chart_forest && accept_arc) ||
            ((!policy || !policy->cross_chart_forest) &&
             (pairs[i].local_link || accept_arc))) {
            int32_t a = pairs[i].a, b = pairs[i].b;
            allowed_partner[a] = b; allowed_partner[b] = a;
            if (stats) stats->accepted_pairs++;
        }
    }
    if (stats) stats->arc_links = n_links;

    {
        BawFilter filter;
        filter.nv = nv; filter.allowed_partner = allowed_partner;
        filter.stats = stats;
        Weld_verts_filtered(arena, verts, nv, NULL, faces, nf, out_nf,
                            eps, true, merge_group,
                            baw_pair_filter, &filter,
                            out_verts, out_nv, NULL);
    }
    rc = 0;

cleanup:
    free(chart_parent); free(chart_candidate);
    free(cluster_nondisk); free(cluster_self); free(cluster_chart_bad);
    free(cluster_chart_b); free(cluster_chart_a);
    free(cluster_accept); free(cluster_edge_zipper);
    free(pair_degree); free(cluster_links);
    free(cluster_support); free(cluster_nodes); free(links); free(pair_parent);
    free(allowed_partner); free(partner); free(pair_id); free(pairs);
    free(second_score); free(best_score); free(best); free(cells);
    free(tangent); free(normal); free(bdeg); free(next); free(prev); free(he);
    return rc;
}

int BoundaryArcWeld_process(Arena_T arena,
                            const float *verts, size_t nv,
                            int32_t *faces, size_t nf,
                            float eps,
                            const int32_t *merge_group,
                            const BoundaryArcWeldParams *params,
                            float **out_verts, size_t *out_nv,
                            size_t *out_nf,
                            BoundaryArcWeldStats *stats)
{
    return BoundaryArcWeld_process_with_policy(
        arena, verts, nv, faces, nf, eps, merge_group, params, NULL,
        out_verts, out_nv, out_nf, stats);
}

static int baw_test_strips(int same_side, int split_groups,
                           size_t *out_nv, BoundaryArcWeldStats *out_stats)
{
    enum { N = 6, NV = N*4, NF = (N-1)*4 };
    float verts[NV*3];
    int32_t faces[NF*3], group[NV];
    size_t nf = 0;
    for (int i = 0; i < N; i++) {
        float x = (float)(i*5);
        float row[4][3] = {
            {0, 0, x}, {0, -1, x}, {0, 0, x},
            {0, same_side ? -2.0f : 1.0f, x}
        };
        for (int r = 0; r < 4; r++) {
            int v = r*N+i;
            memcpy(&verts[v*3], row[r], 3*sizeof(float));
            group[v] = (split_groups && r >= 2) ? 1 : 0;
        }
    }
    for (int i = 0; i < N-1; i++) {
        int32_t a=i, an=i+1, c=N+i, cn=N+i+1;
        faces[nf*3+0]=c; faces[nf*3+1]=a; faces[nf*3+2]=an; nf++;
        faces[nf*3+0]=c; faces[nf*3+1]=an; faces[nf*3+2]=cn; nf++;
        int32_t b=2*N+i, bn=2*N+i+1, d=3*N+i, dn=3*N+i+1;
        if (!same_side) {
            faces[nf*3+0]=b; faces[nf*3+1]=d; faces[nf*3+2]=dn; nf++;
            faces[nf*3+0]=b; faces[nf*3+1]=dn; faces[nf*3+2]=bn; nf++;
        } else {
            faces[nf*3+0]=d; faces[nf*3+1]=b; faces[nf*3+2]=bn; nf++;
            faces[nf*3+0]=d; faces[nf*3+1]=bn; faces[nf*3+2]=dn; nf++;
        }
    }
    Arena_T arena = Arena_new();
    float *welded = NULL; size_t welded_nv = 0, welded_nf = 0;
    BoundaryArcWeldStats st;
    int rc = BoundaryArcWeld_process(arena, verts, NV, faces, nf, 1e-3f,
                                     group, NULL, &welded, &welded_nv,
                                     &welded_nf, &st);
    if (out_nv) *out_nv = welded_nv;
    if (out_stats) *out_stats = st;
    Arena_dispose(&arena);
    return rc;
}

int BoundaryArcWeld_selftest(void)
{
    int fails = 0;
    {
        size_t nv = 0; BoundaryArcWeldStats s;
        if (baw_test_strips(0, 0, &nv, &s) != 0 ||
            s.accepted_pairs < 4 || nv > 20) {
            fprintf(stderr, "[selftest] boundary arc positive: pairs=%zu "
                    "nv=%zu (want >=4, <=20) -> FAIL\n",
                    s.accepted_pairs, nv);
            fails++;
        }
    }
    {
        size_t nv = 0; BoundaryArcWeldStats s;
        if (baw_test_strips(1, 0, &nv, &s) != 0 ||
            s.accepted_pairs != 0 || nv != 24) {
            fprintf(stderr, "[selftest] same-facing boundary contact: pairs=%zu "
                    "nv=%zu (want 0, 24) -> FAIL\n",
                    s.accepted_pairs, nv);
            fails++;
        }
    }
    {
        size_t nv = 0; BoundaryArcWeldStats s;
        if (baw_test_strips(0, 1, &nv, &s) != 0 ||
            s.accepted_pairs != 0 || nv != 24) {
            fprintf(stderr, "[selftest] cross-component boundary contact: "
                    "pairs=%zu nv=%zu (want 0, 24) -> FAIL\n",
                    s.accepted_pairs, nv);
            fails++;
        }
    }
    {
        /* One perfectly oriented point contact has no reciprocal arc support
         * and no shared boundary-link endpoint: it must remain split. */
        const float verts[18] = {
            0,-1,-1,  0,0,0,  0,-1,1,
            0, 1, 1,  0,0,0,  0, 1,-1
        };
        int32_t faces[6] = {0,1,2, 3,4,5};
        int32_t group[6] = {0,0,0,0,0,0};
        Arena_T arena = Arena_new();
        float *w = NULL; size_t nv = 0, nf = 0;
        BoundaryArcWeldStats s;
        if (BoundaryArcWeld_process(arena, verts, 6, faces, 2, 1e-3f,
                                    group, NULL, &w, &nv, &nf, &s) != 0 ||
            s.accepted_pairs != 0 || nv != 6) {
            fprintf(stderr, "[selftest] isolated opposite contact: pairs=%zu "
                    "nv=%zu (want 0, 6) -> FAIL\n", s.accepted_pairs, nv);
            fails++;
        }
        Arena_dispose(&arena);
    }
    if (!fails)
        fprintf(stderr, "[selftest] boundary-arc weld -> ok "
                "(supported seam accepted; ambiguous/isolated/layer contacts rejected)\n");
    return fails;
}
