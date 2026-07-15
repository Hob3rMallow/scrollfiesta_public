#include "ves_platform.h"
#include "qem.h"
#include "kdtree.h"
#include "csr.h"
#include "pipeline_constants.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * Quadric: symmetric 4x4 stored as 10 doubles (upper triangle)
 *
 *   [ q[0] q[1] q[2] q[3] ]
 *   [      q[4] q[5] q[6] ]
 *   [           q[7] q[8] ]
 *   [                q[9] ]
 * ================================================================ */

typedef struct {
    double q[10];
} Quadric;

static void quadric_zero(Quadric *out)
{
    memset(out->q, 0, sizeof(out->q));
}

static void quadric_add(Quadric *dst, const Quadric *src)
{
    size_t i = 0;
    for (i = 0; i < 10; i++) {
        dst->q[i] += src->q[i];
    }
}

/* Build quadric from plane equation ax + by + cz + d = 0 */
static void quadric_from_plane(Quadric *out, double a, double b, double c, double d)
{
    out->q[0] = a * a;  out->q[1] = a * b;  out->q[2] = a * c;  out->q[3] = a * d;
    out->q[4] = b * b;  out->q[5] = b * c;  out->q[6] = b * d;
    out->q[7] = c * c;  out->q[8] = c * d;
    out->q[9] = d * d;
}

/* Evaluate v^T Q v for point (x, y, z) */
static double quadric_eval(const Quadric *Q, double x, double y, double z)
{
    return Q->q[0]*x*x + 2.0*Q->q[1]*x*y + 2.0*Q->q[2]*x*z + 2.0*Q->q[3]*x
         + Q->q[4]*y*y + 2.0*Q->q[5]*y*z + 2.0*Q->q[6]*y
         + Q->q[7]*z*z + 2.0*Q->q[8]*z
         + Q->q[9];
}

/* Solve 3x3 upper-left of Q for optimal position via Cramer's rule.
 * Returns 0 on success, -1 if singular (caller uses midpoint fallback). */
static int quadric_solve(const Quadric *Q, float out[3])
{
    /* System: A * p = -b where A is upper-left 3x3, b is [q3, q6, q8] */
    double a00 = Q->q[0], a01 = Q->q[1], a02 = Q->q[2];
    double a11 = Q->q[4], a12 = Q->q[5];
    double a22 = Q->q[7];
    double b0 = -Q->q[3], b1 = -Q->q[6], b2 = -Q->q[8];

    /* Determinant by cofactor expansion */
    double c00 = a11 * a22 - a12 * a12;
    double c01 = a01 * a22 - a02 * a12;
    double c02 = a01 * a12 - a02 * a11;

    double det = a00 * c00 - a01 * c01 + a02 * c02;

    if (fabs(det) < QEM_DET_THRESHOLD) return -1;

    double inv_det = 1.0 / det;

    /* Cramer's rule — replace each column of A with b */
    double dx = b0 * c00 - a01 * (b1 * a22 - b2 * a12) + a02 * (b1 * a12 - b2 * a11);
    double dy = a00 * (b1 * a22 - b2 * a12) - a01 * (b0 * a22 - b2 * a02) + a02 * (b0 * a12 - b1 * a02);
    double dz = a00 * (a11 * b2 - a12 * b1) - a01 * (a01 * b2 - a12 * b0) + a02 * (a01 * b1 - a11 * b0);

    out[0] = (float)(dx * inv_det);
    out[1] = (float)(dy * inv_det);
    out[2] = (float)(dz * inv_det);
    return 0;
}

/* ================================================================
 * Per-vertex state
 * ================================================================ */

typedef struct {
    float    pos[3];
    Quadric  Q;
    int32_t  remap;     /* collapse chain: who this vertex was merged into */
    uint8_t  alive;
} QEM_Vertex;

/* ================================================================
 * Edge with collapse cost
 * ================================================================ */

typedef struct {
    int32_t  v0, v1;       /* v0 < v1 always */
    double   cost;
    float    opt_pos[3];   /* optimal collapse position */
    uint8_t  face_count;   /* 1 = boundary, 2 = interior */
    uint8_t  dead;
} QEM_Edge;

/* ================================================================
 * Adjacency: vertex -> list of incident edge indices
 * ================================================================ */

typedef struct AdjNode {
    int32_t edge_idx;
    int32_t next;       /* index into adj_pool, -1 = end */
} AdjNode;

typedef struct {
    int32_t  *head;     /* [nv], head index per vertex, -1 = empty */
    AdjNode  *pool;
    int32_t   pool_size;
    int32_t   pool_cap;
} AdjList;

static void adj_add(AdjList *adj, int32_t vertex, int32_t edge_idx)
{
    assert(adj->pool_size < adj->pool_cap);
    int32_t ni = adj->pool_size++;
    adj->pool[ni].edge_idx = edge_idx;
    adj->pool[ni].next = adj->head[vertex];
    adj->head[vertex] = ni;
}

/* ================================================================
 * Face adjacency: vertex -> list of incident face indices
 * Used for normal-flip check during edge collapse.
 * ================================================================ */

typedef struct {
    int32_t face_idx;
    int32_t next;       /* index into pool, -1 = end */
} FaceAdjNode;

typedef struct {
    int32_t      *head;     /* [nv], head index per vertex, -1 = empty */
    FaceAdjNode  *pool;
    int32_t       pool_size;
    int32_t       pool_cap;
} FaceAdj;

static void fadj_add(FaceAdj *fa, int32_t vertex, int32_t face_idx)
{
    assert(fa->pool_size < fa->pool_cap);
    int32_t ni = fa->pool_size++;
    fa->pool[ni].face_idx = face_idx;
    fa->pool[ni].next = fa->head[vertex];
    fa->head[vertex] = ni;
}

/* ================================================================
 * Half-edge for dedup sort
 * ================================================================ */

typedef struct {
    int32_t v0, v1;     /* v0 < v1 */
    int32_t face_idx;
} HalfEdge;

static int halfedge_cmp(const void *a, const void *b)
{
    const HalfEdge *ha = (const HalfEdge *)a;
    const HalfEdge *hb = (const HalfEdge *)b;
    if (ha->v0 != hb->v0) return (ha->v0 < hb->v0) ? -1 : 1;
    if (ha->v1 != hb->v1) return (ha->v1 < hb->v1) ? -1 : 1;
    return 0;
}

/* ================================================================
 * Compute edge cost and optimal position
 * ================================================================ */

static void compute_edge_cost(QEM_Edge *e, const QEM_Vertex *verts,
                              const uint8_t *vert_pinned)
{
    Quadric Qsum;
    quadric_zero(&Qsum);
    quadric_add(&Qsum, &verts[e->v0].Q);
    quadric_add(&Qsum, &verts[e->v1].Q);

    /* Pinned-vertex handling: if both endpoints are pinned the edge is
     * un-collapsible (would shrink a boundary loop or move a Dirichlet
     * vertex). If exactly one endpoint is pinned, force opt_pos to the
     * pinned vertex so the survivor of the collapse is bit-identical to
     * its pre-collapse position. */
    int p0 = vert_pinned ? (vert_pinned[e->v0] != 0) : 0;
    int p1 = vert_pinned ? (vert_pinned[e->v1] != 0) : 0;

    if (p0 && p1) {
        e->dead = 1;
        e->opt_pos[0] = verts[e->v0].pos[0];
        e->opt_pos[1] = verts[e->v0].pos[1];
        e->opt_pos[2] = verts[e->v0].pos[2];
        e->cost = 0.0;  /* never popped; dead is the gate */
        return;
    }
    if (p0) {
        e->opt_pos[0] = verts[e->v0].pos[0];
        e->opt_pos[1] = verts[e->v0].pos[1];
        e->opt_pos[2] = verts[e->v0].pos[2];
    } else if (p1) {
        e->opt_pos[0] = verts[e->v1].pos[0];
        e->opt_pos[1] = verts[e->v1].pos[1];
        e->opt_pos[2] = verts[e->v1].pos[2];
    } else if (quadric_solve(&Qsum, e->opt_pos) == 0) {
        /* Try optimal solve. Clamp: if optimal point is too far from
         * midpoint, use midpoint. */
        float mx = 0.5f * (verts[e->v0].pos[0] + verts[e->v1].pos[0]);
        float my = 0.5f * (verts[e->v0].pos[1] + verts[e->v1].pos[1]);
        float mz = 0.5f * (verts[e->v0].pos[2] + verts[e->v1].pos[2]);
        float dx = e->opt_pos[0] - mx;
        float dy = e->opt_pos[1] - my;
        float dz = e->opt_pos[2] - mz;
        float ex = verts[e->v0].pos[0] - verts[e->v1].pos[0];
        float ey = verts[e->v0].pos[1] - verts[e->v1].pos[1];
        float ez = verts[e->v0].pos[2] - verts[e->v1].pos[2];
        float edge_len_sq = ex*ex + ey*ey + ez*ez;
        float opt_dist_sq = dx*dx + dy*dy + dz*dz;
        if (opt_dist_sq > QEM_DISPLACEMENT_CLAMP * QEM_DISPLACEMENT_CLAMP * edge_len_sq) {
            e->opt_pos[0] = mx;
            e->opt_pos[1] = my;
            e->opt_pos[2] = mz;
        }
    } else {
        /* Fallback: midpoint */
        e->opt_pos[0] = 0.5f * (verts[e->v0].pos[0] + verts[e->v1].pos[0]);
        e->opt_pos[1] = 0.5f * (verts[e->v0].pos[1] + verts[e->v1].pos[1]);
        e->opt_pos[2] = 0.5f * (verts[e->v0].pos[2] + verts[e->v1].pos[2]);
    }

    e->cost = quadric_eval(&Qsum, (double)e->opt_pos[0],
                                  (double)e->opt_pos[1],
                                  (double)e->opt_pos[2]);
    /* Ensure non-negative cost */
    if (e->cost < 0.0) e->cost = 0.0;
}

/* ================================================================
 * Find root of remap chain with path compression
 * ================================================================ */

static int32_t find_root(QEM_Vertex *verts, int32_t v)
{
    while (verts[v].remap != v) {
        verts[v].remap = verts[verts[v].remap].remap;  /* path compression */
        v = verts[v].remap;
    }
    return v;
}

/* ================================================================
 * Normal-flip check for edge collapse
 *
 * Before accepting a collapse of rv1 into rv0 at opt_pos, verify
 * that no incident triangle would have its normal reversed.
 * Flipped normals cause self-intersecting geometry.
 *
 * Returns 1 if any flip detected (reject collapse), 0 if safe.
 * ================================================================ */

static int collapse_would_flip(
    QEM_Vertex      *verts,
    const int32_t   *in_faces,
    const uint8_t   *face_alive,
    const FaceAdj   *fadj,
    int32_t          rv0,
    int32_t          rv1,
    const float      opt_pos[3])
{
    int32_t check_verts[2];
    int cvi = 0;
    check_verts[0] = rv0;
    check_verts[1] = rv1;

    for (cvi = 0; cvi < 2; cvi++) {
        int32_t ni = fadj->head[check_verts[cvi]];
        while (ni >= 0) {
            int32_t fi = fadj->pool[ni].face_idx;
            ni = fadj->pool[ni].next;

            if (!face_alive[fi]) continue;

            /* Resolve face vertices to current roots */
            int32_t ri0 = find_root(verts, in_faces[fi * 3 + 0]);
            int32_t ri1 = find_root(verts, in_faces[fi * 3 + 1]);
            int32_t ri2 = find_root(verts, in_faces[fi * 3 + 2]);

            /* Skip already-degenerate faces */
            if (ri0 == ri1 || ri1 == ri2 || ri0 == ri2) continue;

            /* Count how many vertices are rv0 or rv1 */
            int touches = 0;
            if (ri0 == rv0 || ri0 == rv1) touches++;
            if (ri1 == rv0 || ri1 == rv1) touches++;
            if (ri2 == rv0 || ri2 == rv1) touches++;
            if (touches == 0) continue;  /* not incident */
            if (touches >= 2) continue;  /* will become degenerate — skip */

            /* Get current vertex positions */
            float p0[3], p1[3], p2[3];
            p0[0] = verts[ri0].pos[0]; p0[1] = verts[ri0].pos[1]; p0[2] = verts[ri0].pos[2];
            p1[0] = verts[ri1].pos[0]; p1[1] = verts[ri1].pos[1]; p1[2] = verts[ri1].pos[2];
            p2[0] = verts[ri2].pos[0]; p2[1] = verts[ri2].pos[1]; p2[2] = verts[ri2].pos[2];

            /* Compute old normal */
            float e1x = p1[0] - p0[0], e1y = p1[1] - p0[1], e1z = p1[2] - p0[2];
            float e2x = p2[0] - p0[0], e2y = p2[1] - p0[1], e2z = p2[2] - p0[2];
            float nox = e1y * e2z - e1z * e2y;
            float noy = e1z * e2x - e1x * e2z;
            float noz = e1x * e2y - e1y * e2x;
            float old_len_sq = nox * nox + noy * noy + noz * noz;
            if (old_len_sq < 1e-20f) continue;  /* zero-area triangle */

            /* Replace the affected vertex with opt_pos */
            if (ri0 == rv0 || ri0 == rv1) {
                p0[0] = opt_pos[0]; p0[1] = opt_pos[1]; p0[2] = opt_pos[2];
            }
            if (ri1 == rv0 || ri1 == rv1) {
                p1[0] = opt_pos[0]; p1[1] = opt_pos[1]; p1[2] = opt_pos[2];
            }
            if (ri2 == rv0 || ri2 == rv1) {
                p2[0] = opt_pos[0]; p2[1] = opt_pos[1]; p2[2] = opt_pos[2];
            }

            /* Compute new normal */
            float ne1x = p1[0] - p0[0], ne1y = p1[1] - p0[1], ne1z = p1[2] - p0[2];
            float ne2x = p2[0] - p0[0], ne2y = p2[1] - p0[1], ne2z = p2[2] - p0[2];
            float nnx = ne1y * ne2z - ne1z * ne2y;
            float nny = ne1z * ne2x - ne1x * ne2z;
            float nnz = ne1x * ne2y - ne1y * ne2x;
            float new_len_sq = nnx * nnx + nny * nny + nnz * nnz;

            /* Reject if new face collapses to zero area */
            if (new_len_sq < 1e-20f) return 1;

            /* Reject if the normal swings past the fold threshold. Normalize the
             * dot of the (unnormalized) old/new normals to a cosine so the test
             * is on the ANGLE, not area-weighted -- this catches sub-90deg
             * fold-starts that accumulate into a self-intersection across LME
             * rounds, not just a full (>90deg) reversal. */
            float dot = nox * nnx + noy * nny + noz * nnz;
            float denom = sqrtf(old_len_sq * new_len_sq);
            if (denom < 1e-20f) return 1;
            if (dot / denom < QEM_FLIP_COS_THRESHOLD) return 1;
        }
    }
    return 0;
}

/* ================================================================
 * Link condition for edge collapse (Dey, Edelsbrunner, Guha 1999)
 *
 * Collapsing edge (rv0, rv1) is manifold-safe iff the only vertices
 * adjacent to BOTH endpoints are the apex vertices of the faces that
 * share the edge. A common neighbour w that is NOT such an apex means
 * edges (rv0,w) and (rv1,w) exist without (rv0,rv1,w) being a face;
 * after the merge, edge (survivor,w) inherits faces from both sides and
 * carries a THIRD face — a non-manifold edge. The normal-flip and
 * proximity guards do not catch this purely topological case (it was the
 * source of the raw_snap non-manifold edges).
 *
 * Returns 1 if the collapse would violate the link condition (reject),
 * 0 if safe. Same fadj / find_root machinery as collapse_would_flip.
 * Conservatively rejects if either neighbourhood exceeds LINK_MAXDEG, so
 * manifoldness is never risked for a degree we cannot buffer.
 * ================================================================ */

#define LINK_MAXDEG 64

static int collapse_violates_link(
    QEM_Vertex      *verts,
    const int32_t   *in_faces,
    const uint8_t   *face_alive,
    const FaceAdj   *fadj,
    int32_t          rv0,
    int32_t          rv1)
{
    int32_t nbr0[LINK_MAXDEG];
    int     n0 = 0;
    int     shared_faces = 0;

    /* Unique neighbour roots of rv0 from its alive incident faces, plus a
     * count of faces containing BOTH endpoints (the shared/edge faces). */
    for (int32_t ni = fadj->head[rv0]; ni >= 0; ni = fadj->pool[ni].next) {
        int32_t fi = fadj->pool[ni].face_idx;
        if (!face_alive[fi]) continue;
        int32_t r0 = find_root(verts, in_faces[fi * 3 + 0]);
        int32_t r1 = find_root(verts, in_faces[fi * 3 + 1]);
        int32_t r2 = find_root(verts, in_faces[fi * 3 + 2]);
        if (r0 == r1 || r1 == r2 || r0 == r2) continue;  /* degenerate */
        int has0 = (r0 == rv0) || (r1 == rv0) || (r2 == rv0);
        if (!has0) continue;
        int has1 = (r0 == rv1) || (r1 == rv1) || (r2 == rv1);
        if (has1) shared_faces++;
        int32_t fr[3] = { r0, r1, r2 };
        for (int k = 0; k < 3; k++) {
            int32_t w = fr[k];
            if (w == rv0 || w == rv1) continue;
            int dup = 0;
            for (int t = 0; t < n0; t++) { if (nbr0[t] == w) { dup = 1; break; } }
            if (!dup) {
                if (n0 >= LINK_MAXDEG) return 1;   /* too complex -> reject (safe) */
                nbr0[n0++] = w;
            }
        }
    }

    /* Count UNIQUE common neighbours: rv1 neighbours that appear in nbr0. */
    int32_t seen1[LINK_MAXDEG];
    int     s1 = 0;
    int     common = 0;
    for (int32_t ni = fadj->head[rv1]; ni >= 0; ni = fadj->pool[ni].next) {
        int32_t fi = fadj->pool[ni].face_idx;
        if (!face_alive[fi]) continue;
        int32_t r0 = find_root(verts, in_faces[fi * 3 + 0]);
        int32_t r1 = find_root(verts, in_faces[fi * 3 + 1]);
        int32_t r2 = find_root(verts, in_faces[fi * 3 + 2]);
        if (r0 == r1 || r1 == r2 || r0 == r2) continue;
        int has1 = (r0 == rv1) || (r1 == rv1) || (r2 == rv1);
        if (!has1) continue;
        int32_t fr[3] = { r0, r1, r2 };
        for (int k = 0; k < 3; k++) {
            int32_t w = fr[k];
            if (w == rv0 || w == rv1) continue;
            int seen = 0;
            for (int t = 0; t < s1; t++) { if (seen1[t] == w) { seen = 1; break; } }
            if (seen) continue;
            if (s1 >= LINK_MAXDEG) return 1;
            seen1[s1++] = w;
            for (int t = 0; t < n0; t++) {
                if (nbr0[t] == w) { common++; break; }
            }
        }
    }

    /* Safe iff every common neighbour is an apex of a shared face. A common
     * neighbour beyond the shared-face apexes would become a non-manifold edge. */
    return (common > shared_faces) ? 1 : 0;
}

/* ================================================================
 * Normal consistency check: reject if vertex normals oppose
 *
 * Returns 1 if normals oppose (reject collapse), 0 if compatible.
 * ================================================================ */

static int normals_oppose(const float *vert_normals, int32_t v0, int32_t v1)
{
    float dot = vert_normals[v0 * 3 + 0] * vert_normals[v1 * 3 + 0]
              + vert_normals[v0 * 3 + 1] * vert_normals[v1 * 3 + 1]
              + vert_normals[v0 * 3 + 2] * vert_normals[v1 * 3 + 2];
    return dot < QEM_NORMAL_DOT_THRESHOLD;
}

/* ================================================================
 * Float comparison for qsort (median edge length computation)
 * ================================================================ */

static int float_cmp(const void *a, const void *b)
{
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

/* ================================================================
 * Build KD-tree from alive vertices for proximity checking.
 * Returns NULL if no alive vertices.
 * ================================================================ */

static KDTree_T build_prox_tree(Arena_T arena, const QEM_Vertex *verts,
                                 size_t in_nv,
                                 float **kd_points_out,
                                 int32_t **kd_to_vert_out,
                                 size_t *kd_count_out)
{
    /* Count alive vertices */
    size_t alive = 0;
    size_t i = 0;
    for (i = 0; i < in_nv; i++) {
        if (verts[i].alive) alive++;
    }
    if (alive == 0) {
        *kd_count_out = 0;
        return NULL;
    }

    float *pts = (float *)ARENA_ALLOC(arena, (long)(alive * 3 * sizeof(float)));
    int32_t *map = (int32_t *)ARENA_ALLOC(arena, (long)(alive * sizeof(int32_t)));
    size_t ki = 0;
    for (i = 0; i < in_nv; i++) {
        if (verts[i].alive) {
            pts[ki * 3 + 0] = verts[i].pos[0];
            pts[ki * 3 + 1] = verts[i].pos[1];
            pts[ki * 3 + 2] = verts[i].pos[2];
            map[ki] = (int32_t)i;
            ki++;
        }
    }

    *kd_points_out = pts;
    *kd_to_vert_out = map;
    *kd_count_out = alive;
    return KDTree_new(arena, pts, alive);
}

/* ================================================================
 * Spatial proximity check: reject if opt_pos is too close to
 * a non-neighbor vertex (indicates cross-sheet collision).
 *
 * Returns 1 if too close (reject collapse), 0 if safe.
 * ================================================================ */

static int collapse_too_close(KDTree_T tree, const int32_t *kd_to_vert,
                               const float opt_pos[3], float safe_radius_sq,
                               QEM_Vertex *verts, const FaceAdj *fadj,
                               const int32_t *in_faces, const uint8_t *face_alive,
                               int32_t rv0, int32_t rv1)
{
    if (!tree) return 0;

    int32_t results[64];
    size_t n_results = KDTree_ball_query(tree, opt_pos, safe_radius_sq,
                                          results, 64);
    if (n_results == 0) return 0;

    /* Build exclusion set: all vertex roots from faces incident to rv0 and rv1.
     * Typically ~12-18 unique vertices. Linear scan in fixed-size array. */
    int32_t excl[128];
    int excl_count = 0;

    int32_t check_v[2];
    int cvi = 0;
    check_v[0] = rv0;
    check_v[1] = rv1;
    for (cvi = 0; cvi < 2; cvi++) {
        int32_t ni = fadj->head[check_v[cvi]];
        while (ni >= 0) {
            int32_t fi = fadj->pool[ni].face_idx;
            ni = fadj->pool[ni].next;
            if (!face_alive[fi]) continue;
            int k = 0;
            for (k = 0; k < 3; k++) {
                int32_t r = find_root(verts, in_faces[fi * 3 + k]);
                /* Check if already in excl */
                int found = 0;
                int j = 0;
                for (j = 0; j < excl_count; j++) {
                    if (excl[j] == r) { found = 1; break; }
                }
                if (!found && excl_count < 128) {
                    excl[excl_count++] = r;
                }
            }
        }
    }

    /* Check each ball query result against exclusion set */
    size_t ri = 0;
    for (ri = 0; ri < n_results; ri++) {
        int32_t kd_idx = results[ri];
        int32_t vert_idx = kd_to_vert[kd_idx];
        int32_t root = find_root(verts, vert_idx);
        if (!verts[root].alive) continue;

        /* Check exclusion */
        int excluded = 0;
        int j = 0;
        for (j = 0; j < excl_count; j++) {
            if (excl[j] == root) { excluded = 1; break; }
        }
        if (!excluded) return 1; /* too close to a non-neighbor */
    }

    return 0;
}

/* ================================================================
 * Maintenance helpers: edge flips + tangential smoothing
 * between collapse passes for triangle quality improvement.
 * Adapted from tri_quality.c but self-contained.
 * ================================================================ */

/* Compute unnormalized face normal */
static void qem_maint_face_normal(const float *verts,
                                   int32_t a, int32_t b, int32_t c,
                                   float out[3])
{
    float e1[3], e2[3];
    int k = 0;
    for (k = 0; k < 3; k++) {
        e1[k] = verts[b * 3 + k] - verts[a * 3 + k];
        e2[k] = verts[c * 3 + k] - verts[a * 3 + k];
    }
    out[0] = e1[1] * e2[2] - e1[2] * e2[1];
    out[1] = e1[2] * e2[0] - e1[0] * e2[2];
    out[2] = e1[0] * e2[1] - e1[1] * e2[0];
}

/* Compute area-weighted vertex normals, normalized */
static void qem_maint_vertex_normals(const float *verts, size_t nv,
                                      const int32_t *faces, size_t nf,
                                      float *normals)
{
    size_t i = 0;
    memset(normals, 0, nv * 3 * sizeof(float));

    for (i = 0; i < nf; i++) {
        int32_t a = faces[i * 3 + 0];
        int32_t b = faces[i * 3 + 1];
        int32_t c = faces[i * 3 + 2];
        float fn[3];
        int k = 0;
        qem_maint_face_normal(verts, a, b, c, fn);
        for (k = 0; k < 3; k++) {
            normals[a * 3 + k] += fn[k];
            normals[b * 3 + k] += fn[k];
            normals[c * 3 + k] += fn[k];
        }
    }

    for (i = 0; i < nv; i++) {
        float nx = normals[i * 3 + 0];
        float ny = normals[i * 3 + 1];
        float nz = normals[i * 3 + 2];
        float len = sqrtf(nx * nx + ny * ny + nz * nz);
        if (len > 1e-12f) {
            float inv = 1.0f / len;
            normals[i * 3 + 0] *= inv;
            normals[i * 3 + 1] *= inv;
            normals[i * 3 + 2] *= inv;
        }
    }
}

/* Angle at vertex a in triangle (a, b, c) */
static float qem_maint_angle_at(const float *verts,
                                 int32_t a, int32_t b, int32_t c)
{
    float ab[3], ac[3];
    int k = 0;
    float dot = 0.0f, la = 0.0f, lb = 0.0f, denom = 0.0f, cos_a = 0.0f;
    for (k = 0; k < 3; k++) {
        ab[k] = verts[b * 3 + k] - verts[a * 3 + k];
        ac[k] = verts[c * 3 + k] - verts[a * 3 + k];
    }
    dot = ab[0]*ac[0] + ab[1]*ac[1] + ab[2]*ac[2];
    la = sqrtf(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    lb = sqrtf(ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2]);
    denom = la * lb;
    if (denom < 1e-12f) return 0.0f;
    cos_a = dot / denom;
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;
    return acosf(cos_a);
}

/* Minimum angle of triangle */
static float qem_maint_min_angle(const float *verts,
                                  int32_t a, int32_t b, int32_t c)
{
    float a1 = qem_maint_angle_at(verts, a, b, c);
    float a2 = qem_maint_angle_at(verts, b, c, a);
    float a3 = qem_maint_angle_at(verts, c, a, b);
    float m = a1;
    if (a2 < m) m = a2;
    if (a3 < m) m = a3;
    return m;
}

/* Half-edge entry for maintenance edge flips (needs opposite vertex) */
typedef struct {
    int32_t v0, v1;     /* sorted: v0 < v1 */
    int32_t face;
    int32_t opposite;
    int8_t  fwd;        /* 1 if original edge was v0→v1 in face winding, 0 if v1→v0 */
} MaintHalfEdge;

static int maint_he_cmp(const void *a, const void *b)
{
    const MaintHalfEdge *ea = (const MaintHalfEdge *)a;
    const MaintHalfEdge *eb = (const MaintHalfEdge *)b;
    if (ea->v0 != eb->v0) return (ea->v0 < eb->v0) ? -1 : 1;
    if (ea->v1 != eb->v1) return (ea->v1 < eb->v1) ? -1 : 1;
    return 0;
}

/* True if edge (u,w) already exists in the sorted-by-(v0,v1) half-edge list.
 * An edge flip (a,b)->(c,d) is only manifold-safe if (c,d) is NOT already an
 * edge: otherwise the two flipped faces give edge (c,d) a third/fourth incident
 * face -> a non-manifold edge. (The classic Delaunay-flip precondition; missing
 * here, it was the source of qslim's non-manifold decimation output.) */
static int maint_edge_exists(const MaintHalfEdge *mhe, size_t n_he,
                             int32_t u, int32_t w)
{
    int32_t v0 = (u < w) ? u : w;
    int32_t v1 = (u < w) ? w : u;
    size_t lo = 0, hi = n_he;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (mhe[mid].v0 < v0 || (mhe[mid].v0 == v0 && mhe[mid].v1 < v1))
            lo = mid + 1;
        else
            hi = mid;
    }
    return (lo < n_he && mhe[lo].v0 == v0 && mhe[lo].v1 == v1) ? 1 : 0;
}

/* Detect boundary vertices from face list.
 * Uses existing HalfEdge type and halfedge_cmp from above. */
static void qem_detect_boundary(Arena_T arena,
                                 const int32_t *faces, size_t nf,
                                 size_t nv, uint8_t *is_boundary)
{
    Arena_Mark mark = Arena_save(arena);
    size_t n_he = nf * 3;
    size_t i = 0;

    HalfEdge *he = (HalfEdge *)ARENA_ALLOC(arena,
        (long)(n_he * sizeof(HalfEdge)));

    for (i = 0; i < nf; i++) {
        int32_t fi[3];
        size_t k = 0;
        fi[0] = faces[i * 3 + 0];
        fi[1] = faces[i * 3 + 1];
        fi[2] = faces[i * 3 + 2];
        for (k = 0; k < 3; k++) {
            int32_t a = fi[k];
            int32_t b = fi[(k + 1) % 3];
            if (a > b) { int32_t t = a; a = b; b = t; }
            he[i * 3 + k].v0 = a;
            he[i * 3 + k].v1 = b;
            he[i * 3 + k].face_idx = (int32_t)i;
        }
    }

    qsort(he, n_he, sizeof(HalfEdge), halfedge_cmp);

    memset(is_boundary, 0, nv * sizeof(uint8_t));

    i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && he[j].v0 == he[i].v0
                        && he[j].v1 == he[i].v1) {
            j++;
        }
        if (j - i == 1) {
            /* Boundary edge: appears in only one face */
            is_boundary[he[i].v0] = 1;
            is_boundary[he[i].v1] = 1;
        }
        i = j;
    }

    Arena_restore(arena, mark);
}

/* One pass of Delaunay edge flips.
 * Returns number of flips performed. */
static size_t qem_edge_flip_pass(Arena_T arena,
                                  const float *verts,
                                  int32_t *faces, size_t nf,
                                  size_t nv,
                                  const uint8_t *is_boundary)
{
    Arena_Mark mark = Arena_save(arena);
    size_t n_he = nf * 3;
    size_t n_flipped = 0;
    size_t f = 0;
    size_t i = 0;

    MaintHalfEdge *mhe = (MaintHalfEdge *)ARENA_ALLOC(arena,
        (long)(n_he * sizeof(MaintHalfEdge)));

    for (f = 0; f < nf; f++) {
        int32_t v[3];
        int e = 0;
        v[0] = faces[f * 3 + 0];
        v[1] = faces[f * 3 + 1];
        v[2] = faces[f * 3 + 2];
        for (e = 0; e < 3; e++) {
            int32_t a = v[e];
            int32_t b = v[(e + 1) % 3];
            int32_t opp = v[(e + 2) % 3];
            size_t idx = f * 3 + (size_t)e;
            mhe[idx].v0 = (a < b) ? a : b;
            mhe[idx].v1 = (a < b) ? b : a;
            mhe[idx].face = (int32_t)f;
            mhe[idx].opposite = opp;
            mhe[idx].fwd = (a < b) ? 1 : 0;
        }
    }

    qsort(mhe, n_he, sizeof(MaintHalfEdge), maint_he_cmp);

    uint8_t *vert_used = (uint8_t *)ARENA_ALLOC(arena,
        (long)(nv * sizeof(uint8_t)));
    uint8_t *face_used = (uint8_t *)ARENA_ALLOC(arena,
        (long)(nf * sizeof(uint8_t)));
    memset(vert_used, 0, nv * sizeof(uint8_t));
    memset(face_used, 0, nf * sizeof(uint8_t));

    i = 0;
    while (i + 1 < n_he) {
        if (mhe[i].v0 == mhe[i + 1].v0 && mhe[i].v1 == mhe[i + 1].v1) {
            int32_t a = mhe[i].v0;    /* sorted: a < b */
            int32_t b = mhe[i].v1;
            int32_t c = 0, d = 0;     /* c = fwd-face opposite, d = bwd-face opposite */
            int32_t fc = 0, fd = 0;   /* face indices for c-side and d-side */
            float n_orig[3], n1[3], n2[3];
            float dot1 = 0.0f, dot2 = 0.0f;
            float cur_min = 0.0f, cur_min2 = 0.0f;
            float flip_min = 0.0f, flip_min2 = 0.0f;

            /* Determine winding: fwd half-edge (a→b in face) owns "c",
             * bwd half-edge (b→a in face) owns "d".
             * For consistent orientation, one is fwd and the other is bwd. */
            if (mhe[i].fwd && !mhe[i + 1].fwd) {
                c = mhe[i].opposite;      fc = mhe[i].face;
                d = mhe[i + 1].opposite;  fd = mhe[i + 1].face;
            } else if (!mhe[i].fwd && mhe[i + 1].fwd) {
                c = mhe[i + 1].opposite;  fc = mhe[i + 1].face;
                d = mhe[i].opposite;      fd = mhe[i].face;
            } else {
                /* Both same direction — inconsistent winding, skip */
                i += 2;
                continue;
            }

            /* Skip boundary edges */
            if (is_boundary[a] && is_boundary[b]) {
                i += 2;
                continue;
            }

            /* Skip if already used */
            if (vert_used[a] || vert_used[b] || vert_used[c] || vert_used[d]
                || face_used[fc] || face_used[fd]) {
                i += 2;
                continue;
            }

            /* Topological guard: never flip onto an edge that already exists.
             * If (c,d) is already an edge elsewhere, the flip gives it a third
             * incident face -> non-manifold. The vert_used independence guard
             * above guarantees no flip performed earlier in THIS pass touched
             * c or d, so the original sorted half-edge list `mhe` is the
             * correct edge set to query. */
            if (maint_edge_exists(mhe, n_he, c, d)) {
                i += 2;
                continue;
            }

            /* Check angle improvement.
             * Original faces: (a,b,c) and (b,a,d) ≡ (d,b,a).
             * min_angle only cares about the triangle, not winding. */
            cur_min = qem_maint_min_angle(verts, a, b, c);
            cur_min2 = qem_maint_min_angle(verts, b, a, d);
            if (cur_min2 < cur_min) cur_min = cur_min2;

            flip_min = qem_maint_min_angle(verts, c, a, d);
            flip_min2 = qem_maint_min_angle(verts, d, b, c);
            if (flip_min2 < flip_min) flip_min = flip_min2;

            if (flip_min <= cur_min + 0.001f) {
                i += 2;
                continue;
            }

            /* Convexity: reference normal from fwd face (a,b,c).
             * The orientation-preserving flip of the quad a->d->b->c->a across
             * the new diagonal c-d yields faces (c,a,d) and (d,b,c) -- both keep
             * the quad's boundary-edge directions, so their normals must agree
             * with n_orig. (A reversed write here would negate them and inject a
             * same_dir edge on every quad boundary edge.) */
            qem_maint_face_normal(verts, a, b, c, n_orig);
            qem_maint_face_normal(verts, c, a, d, n1);
            qem_maint_face_normal(verts, d, b, c, n2);
            dot1 = n_orig[0]*n1[0] + n_orig[1]*n1[1] + n_orig[2]*n1[2];
            dot2 = n_orig[0]*n2[0] + n_orig[1]*n2[1] + n_orig[2]*n2[2];
            if (dot1 <= 0.0f || dot2 <= 0.0f) {
                i += 2;
                continue;
            }

            /* Perform flip: (a,b,c) + (b,a,d) -> (c,a,d) + (d,b,c)
             * fc is the face that had (a,b,c), fd had (b,a,d). Winding is
             * preserved: both new tris keep the quad boundary a->d->b->c. */
            faces[fc * 3 + 0] = c;
            faces[fc * 3 + 1] = a;
            faces[fc * 3 + 2] = d;

            faces[fd * 3 + 0] = d;
            faces[fd * 3 + 1] = b;
            faces[fd * 3 + 2] = c;

            vert_used[a] = vert_used[b] = vert_used[c] = vert_used[d] = 1;
            face_used[fc] = face_used[fd] = 1;
            n_flipped++;

            i += 2;
        } else {
            i++;
        }
    }

    Arena_restore(arena, mark);
    return n_flipped;
}

/* One iteration of tangential Laplacian smoothing.
 * Displacement projected onto surface tangent plane (removes normal component).
 * Uses CSR_from_faces for adjacency.
 *
 * Skips boundary verts AND pinned verts (if pin_mask is non-NULL). Without
 * the pin guard the maintenance smooth would displace interior pinned
 * verts, breaking the Dirichlet contract that QEM_simplify_pinned promises
 * to its callers. */
static void qem_smooth_pass(Arena_T arena,
                              float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const float *vert_normals,
                              const uint8_t *is_boundary,
                              const uint8_t *pin_mask,
                              float lambda)
{
    Arena_Mark mark = Arena_save(arena);
    size_t i = 0;

    CSR_T adj = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *off = CSR_offset(adj);
    const int32_t *tgt = CSR_target(adj);

    float *disp = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    memset(disp, 0, nv * 3 * sizeof(float));

    for (i = 0; i < nv; i++) {
        int32_t degree = 0;
        float mean[3] = {0.0f, 0.0f, 0.0f};
        float d[3] = {0.0f, 0.0f, 0.0f};
        float nx = 0.0f, ny = 0.0f, nz = 0.0f, dot = 0.0f;
        int32_t j = 0;
        float inv_d = 0.0f;

        if (is_boundary[i]) continue;
        if (pin_mask && pin_mask[i]) continue;

        degree = off[(int32_t)i + 1] - off[(int32_t)i];
        if (degree == 0) continue;

        for (j = off[(int32_t)i]; j < off[(int32_t)i + 1]; j++) {
            int32_t nb = tgt[j];
            mean[0] += verts[nb * 3 + 0];
            mean[1] += verts[nb * 3 + 1];
            mean[2] += verts[nb * 3 + 2];
        }
        inv_d = 1.0f / (float)degree;
        mean[0] *= inv_d;
        mean[1] *= inv_d;
        mean[2] *= inv_d;

        d[0] = mean[0] - verts[i * 3 + 0];
        d[1] = mean[1] - verts[i * 3 + 1];
        d[2] = mean[2] - verts[i * 3 + 2];

        /* Project out normal component: d -= dot(d, n) * n */
        nx = vert_normals[i * 3 + 0];
        ny = vert_normals[i * 3 + 1];
        nz = vert_normals[i * 3 + 2];
        dot = d[0]*nx + d[1]*ny + d[2]*nz;
        d[0] -= dot * nx;
        d[1] -= dot * ny;
        d[2] -= dot * nz;

        disp[i * 3 + 0] = lambda * d[0];
        disp[i * 3 + 1] = lambda * d[1];
        disp[i * 3 + 2] = lambda * d[2];
    }

    /* Apply displacement */
    for (i = 0; i < nv; i++) {
        if (is_boundary[i]) continue;
        if (pin_mask && pin_mask[i]) continue;
        verts[i * 3 + 0] += disp[i * 3 + 0];
        verts[i * 3 + 1] += disp[i * 3 + 1];
        verts[i * 3 + 2] += disp[i * 3 + 2];
    }

    Arena_restore(arena, mark);
}

/* Full maintenance pass: detect boundary, 1 edge flip, 2 smooth iters.
 * If pin_mask is non-NULL, interior pinned verts are also held fixed (the
 * smooth pass skips them). Edge flips do not move vertices and remain
 * pin-safe by construction. */
static void qem_maintenance_pass(Arena_T arena,
                                  float *verts, size_t nv,
                                  int32_t *faces, size_t nf,
                                  const uint8_t *pin_mask)
{
    Arena_Mark mark = Arena_save(arena);
    int si = 0;
    size_t n_flips = 0;

    uint8_t *is_boundary = (uint8_t *)ARENA_ALLOC(arena,
        (long)(nv * sizeof(uint8_t)));
    qem_detect_boundary(arena, faces, nf, nv, is_boundary);

    float *vnormals = (float *)ARENA_ALLOC(arena,
        (long)(nv * 3 * sizeof(float)));
    qem_maint_vertex_normals(verts, nv, faces, nf, vnormals);

    /* Edge-flip to convergence (Surazhsky-Gotsman min-angle improvement, normal-
     * guarded so geometry never folds). A single pass locks each vert to one flip
     * and so only flips non-adjacent edges; iterate until no edge improves so the
     * heavier QEM's slivers are actually cleared. is_boundary is stable across
     * flips (interior edges flip to interior edges; flips never move verts). */
    for (int fr = 0; fr < QEM_FLIP_MAX_ROUNDS; fr++) {
        size_t k = qem_edge_flip_pass(arena, verts, faces, nf, nv, is_boundary);
        n_flips += k;
        if (k == 0) break;
    }

    /* 2 iterations of gentle tangential smoothing */
    for (si = 0; si < 2; si++) {
        qem_maint_vertex_normals(verts, nv, faces, nf, vnormals);
        qem_smooth_pass(arena, verts, nv, faces, nf, vnormals,
                         is_boundary, pin_mask, QEM_MAINT_SMOOTH_LAMBDA);
    }

#ifdef VESUVIUS_DEBUG
    fprintf(stderr, "    QEM maintenance: %zu flips, 2 smooth iters\n",
            n_flips);
#endif

    Arena_restore(arena, mark);
}

/* ================================================================
 * Inner collapse pass (single QEM pass)
 * ================================================================ */

static int qem_collapse_pass(Arena_T arena,
                              const float *in_verts, size_t in_nv,
                              const int32_t *in_faces, size_t in_nf,
                              const uint8_t *caller_pin_mask,
                              size_t target_nf,
                              float **out_verts, size_t *out_nv,
                              int32_t **out_faces, size_t *out_nf,
                              uint8_t **out_pin_mask)
{
    assert(arena);
    assert(out_verts && out_nv && out_faces && out_nf);

    /* Handle trivial cases */
    if (in_nv == 0 || in_nf == 0) {
        *out_verts = NULL;
        *out_nv = 0;
        *out_faces = NULL;
        *out_nf = 0;
        if (out_pin_mask) *out_pin_mask = NULL;
        return 0;
    }

    if (in_nf <= target_nf) {
        /* Copy input unchanged */
        float *vout = (float *)ARENA_ALLOC(arena, (long)(in_nv * 3 * sizeof(float)));
        memcpy(vout, in_verts, in_nv * 3 * sizeof(float));
        int32_t *fout = (int32_t *)ARENA_ALLOC(arena, (long)(in_nf * 3 * sizeof(int32_t)));
        memcpy(fout, in_faces, in_nf * 3 * sizeof(int32_t));
        *out_verts = vout;
        *out_nv = in_nv;
        *out_faces = fout;
        *out_nf = in_nf;
        if (out_pin_mask) {
            if (caller_pin_mask) {
                uint8_t *pout = (uint8_t *)ARENA_ALLOC(arena,
                    (long)(in_nv * sizeof(uint8_t)));
                memcpy(pout, caller_pin_mask, in_nv * sizeof(uint8_t));
                *out_pin_mask = pout;
            } else {
                *out_pin_mask = NULL;
            }
        }
        return 0;
    }

    /* ---- Phase 1: Init vertices ---- */
    QEM_Vertex *verts = (QEM_Vertex *)ARENA_ALLOC(arena, (long)(in_nv * sizeof(QEM_Vertex)));
    {
        size_t i = 0;
        for (i = 0; i < in_nv; i++) {
            verts[i].pos[0] = in_verts[i * 3 + 0];
            verts[i].pos[1] = in_verts[i * 3 + 1];
            verts[i].pos[2] = in_verts[i * 3 + 2];
            quadric_zero(&verts[i].Q);
            verts[i].remap = (int32_t)i;
            verts[i].alive = 1;
        }
    }

    /* ---- Phase 2: Face quadrics ---- */
    /* Per-face alive tracking */
    uint8_t *face_alive = (uint8_t *)ARENA_ALLOC(arena, (long)(in_nf * sizeof(uint8_t)));
    memset(face_alive, 1, in_nf * sizeof(uint8_t));

    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            int32_t i0 = in_faces[f * 3 + 0];
            int32_t i1 = in_faces[f * 3 + 1];
            int32_t i2 = in_faces[f * 3 + 2];

            double v0x = (double)in_verts[i0 * 3 + 0];
            double v0y = (double)in_verts[i0 * 3 + 1];
            double v0z = (double)in_verts[i0 * 3 + 2];
            double v1x = (double)in_verts[i1 * 3 + 0];
            double v1y = (double)in_verts[i1 * 3 + 1];
            double v1z = (double)in_verts[i1 * 3 + 2];
            double v2x = (double)in_verts[i2 * 3 + 0];
            double v2y = (double)in_verts[i2 * 3 + 1];
            double v2z = (double)in_verts[i2 * 3 + 2];

            /* Cross product for normal */
            double e1x = v1x - v0x, e1y = v1y - v0y, e1z = v1z - v0z;
            double e2x = v2x - v0x, e2y = v2y - v0y, e2z = v2z - v0z;
            double nx = e1y * e2z - e1z * e2y;
            double ny = e1z * e2x - e1x * e2z;
            double nz = e1x * e2y - e1y * e2x;
            double len = sqrt(nx * nx + ny * ny + nz * nz);
            if (len < 1e-15) continue;  /* degenerate face */
            nx /= len; ny /= len; nz /= len;

            double d = -(nx * v0x + ny * v0y + nz * v0z);

            Quadric Qf;
            quadric_from_plane(&Qf, nx, ny, nz, d);

            quadric_add(&verts[i0].Q, &Qf);
            quadric_add(&verts[i1].Q, &Qf);
            quadric_add(&verts[i2].Q, &Qf);
        }
    }

    /* ---- Phase 2a: Probabilistic quadric noise (Trettner & Kobbelt 2020) ----
     * Adds isotropic regularization sigma^2 * ||x - q||^2 to each vertex
     * quadric, making A = nn^T + sigma^2*I always full-rank. This biases
     * optimal positions toward original vertex positions, producing more
     * uniform triangulations. sigma is much smaller than boundary weight
     * so boundary preservation is unaffected. */
    {
        Arena_Mark sigma_mark = Arena_save(arena);
        size_t sample_faces = in_nf < 1000 ? in_nf : 1000;
        size_t sample_edges = sample_faces * 3;
        float *prob_edge_lens = (float *)ARENA_ALLOC(arena,
            (long)(sample_edges * sizeof(float)));
        size_t pstep = in_nf / sample_faces;
        size_t psi = 0;
        size_t pfi = 0;
        double sigma_n_sq = 0.0;

        if (pstep < 1) pstep = 1;

        for (pfi = 0; psi < sample_edges && pfi < in_nf; pfi += pstep) {
            int32_t i0 = in_faces[pfi * 3 + 0];
            int32_t i1 = in_faces[pfi * 3 + 1];
            int32_t i2 = in_faces[pfi * 3 + 2];
            float dx, dy, dz;
            /* Edge 0-1 */
            dx = in_verts[i1*3+0] - in_verts[i0*3+0];
            dy = in_verts[i1*3+1] - in_verts[i0*3+1];
            dz = in_verts[i1*3+2] - in_verts[i0*3+2];
            prob_edge_lens[psi++] = sqrtf(dx*dx + dy*dy + dz*dz);
            /* Edge 1-2 */
            dx = in_verts[i2*3+0] - in_verts[i1*3+0];
            dy = in_verts[i2*3+1] - in_verts[i1*3+1];
            dz = in_verts[i2*3+2] - in_verts[i1*3+2];
            prob_edge_lens[psi++] = sqrtf(dx*dx + dy*dy + dz*dz);
            /* Edge 2-0 */
            dx = in_verts[i0*3+0] - in_verts[i2*3+0];
            dy = in_verts[i0*3+1] - in_verts[i2*3+1];
            dz = in_verts[i0*3+2] - in_verts[i2*3+2];
            prob_edge_lens[psi++] = sqrtf(dx*dx + dy*dy + dz*dz);
        }

        if (psi > 0) {
            qsort(prob_edge_lens, psi, sizeof(float), float_cmp);
            double sigma_n = QEM_PROB_SIGMA_FACTOR * (double)prob_edge_lens[psi / 2];
            sigma_n_sq = sigma_n * sigma_n;
        }
        Arena_restore(arena, sigma_mark);

        /* Add sigma^2 * ||x - q||^2 to each vertex quadric */
        if (sigma_n_sq > 0.0) {
            size_t vi = 0;
            for (vi = 0; vi < in_nv; vi++) {
                double qx = (double)verts[vi].pos[0];
                double qy = (double)verts[vi].pos[1];
                double qz = (double)verts[vi].pos[2];
                /* A diagonal += sigma^2 */
                verts[vi].Q.q[0] += sigma_n_sq;
                verts[vi].Q.q[4] += sigma_n_sq;
                verts[vi].Q.q[7] += sigma_n_sq;
                /* b terms -= sigma^2 * q */
                verts[vi].Q.q[3] -= sigma_n_sq * qx;
                verts[vi].Q.q[6] -= sigma_n_sq * qy;
                verts[vi].Q.q[8] -= sigma_n_sq * qz;
                /* c term += sigma^2 * ||q||^2 */
                verts[vi].Q.q[9] += sigma_n_sq * (qx*qx + qy*qy + qz*qz);
            }
        }
    }

    /* ---- Phase 2b: Compute area-weighted vertex normals ---- */
    float *vert_normals = (float *)ARENA_ALLOC(arena, (long)(in_nv * 3 * sizeof(float)));
    memset(vert_normals, 0, in_nv * 3 * sizeof(float));
    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            int32_t i0 = in_faces[f * 3 + 0];
            int32_t i1 = in_faces[f * 3 + 1];
            int32_t i2 = in_faces[f * 3 + 2];
            float e1x = in_verts[i1*3+0] - in_verts[i0*3+0];
            float e1y = in_verts[i1*3+1] - in_verts[i0*3+1];
            float e1z = in_verts[i1*3+2] - in_verts[i0*3+2];
            float e2x = in_verts[i2*3+0] - in_verts[i0*3+0];
            float e2y = in_verts[i2*3+1] - in_verts[i0*3+1];
            float e2z = in_verts[i2*3+2] - in_verts[i0*3+2];
            float nx = e1y*e2z - e1z*e2y;
            float ny = e1z*e2x - e1x*e2z;
            float nz = e1x*e2y - e1y*e2x;
            /* Area-weighted: unnormalized cross product = 2x triangle area */
            vert_normals[i0*3+0] += nx; vert_normals[i0*3+1] += ny; vert_normals[i0*3+2] += nz;
            vert_normals[i1*3+0] += nx; vert_normals[i1*3+1] += ny; vert_normals[i1*3+2] += nz;
            vert_normals[i2*3+0] += nx; vert_normals[i2*3+1] += ny; vert_normals[i2*3+2] += nz;
        }
        /* Normalize to unit length */
        size_t i = 0;
        for (i = 0; i < in_nv; i++) {
            float nx = vert_normals[i*3+0];
            float ny = vert_normals[i*3+1];
            float nz = vert_normals[i*3+2];
            float len = sqrtf(nx*nx + ny*ny + nz*nz);
            if (len > 1e-10f) {
                vert_normals[i*3+0] = nx / len;
                vert_normals[i*3+1] = ny / len;
                vert_normals[i*3+2] = nz / len;
            }
        }
    }

    /* ---- Phase 3: Build edges via half-edge sort+dedup ---- */
    size_t n_halfedges = in_nf * 3;
    HalfEdge *halfedges = (HalfEdge *)ARENA_ALLOC(arena, (long)(n_halfedges * sizeof(HalfEdge)));
    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            int32_t fi[3];
            fi[0] = in_faces[f * 3 + 0];
            fi[1] = in_faces[f * 3 + 1];
            fi[2] = in_faces[f * 3 + 2];
            size_t k = 0;
            for (k = 0; k < 3; k++) {
                int32_t a = fi[k];
                int32_t b = fi[(k + 1) % 3];
                if (a > b) { int32_t t = a; a = b; b = t; }
                halfedges[f * 3 + k].v0 = a;
                halfedges[f * 3 + k].v1 = b;
                halfedges[f * 3 + k].face_idx = (int32_t)f;
            }
        }
    }

    qsort(halfedges, n_halfedges, sizeof(HalfEdge), halfedge_cmp);

    /* Count unique edges */
    size_t n_edges = 0;
    {
        size_t i = 0;
        while (i < n_halfedges) {
            size_t j = i + 1;
            while (j < n_halfedges && halfedges[j].v0 == halfedges[i].v0
                                   && halfedges[j].v1 == halfedges[i].v1) {
                j++;
            }
            n_edges++;
            i = j;
        }
    }

    QEM_Edge *edges = (QEM_Edge *)ARENA_ALLOC(arena, (long)(n_edges * sizeof(QEM_Edge)));
    {
        size_t ei = 0;
        size_t i = 0;
        while (i < n_halfedges) {
            int32_t v0 = halfedges[i].v0;
            int32_t v1 = halfedges[i].v1;
            size_t j = i + 1;
            while (j < n_halfedges && halfedges[j].v0 == v0 && halfedges[j].v1 == v1) {
                j++;
            }
            uint8_t fc = (uint8_t)(j - i);
            if (fc > 2) fc = 2;  /* cap at 2 (non-manifold edges) */

            edges[ei].v0 = v0;
            edges[ei].v1 = v1;
            edges[ei].face_count = fc;
            edges[ei].dead = 0;
            ei++;
            i = j;
        }
        assert(ei == n_edges);
    }

    /* ---- Phase 4: Boundary penalty ---- */
    {
        size_t i = 0;
        for (i = 0; i < n_edges; i++) {
            if (edges[i].face_count != 1) continue;  /* only boundary edges */

            int32_t v0 = edges[i].v0;
            int32_t v1 = edges[i].v1;

            /* Edge direction */
            double ex = (double)verts[v1].pos[0] - (double)verts[v0].pos[0];
            double ey = (double)verts[v1].pos[1] - (double)verts[v0].pos[1];
            double ez = (double)verts[v1].pos[2] - (double)verts[v0].pos[2];
            double elen = sqrt(ex * ex + ey * ey + ez * ez);
            if (elen < 1e-15) continue;
            ex /= elen; ey /= elen; ez /= elen;

            /* Average vertex quadric normal approximation:
             * Use the gradient of the summed quadric as a normal proxy.
             * For boundary preservation we need a plane perpendicular to
             * the boundary edge. We use two perpendicular planes. */

            /* Plane 1: perpendicular in XY plane (arbitrary but stable) */
            /* Use cross product of edge direction with each axis to get
             * two orthogonal planes containing the edge */
            double n1x = ey, n1y = -ex, n1z = 0.0;
            double n1len = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
            if (n1len < 1e-15) {
                /* Edge is along Z axis, use different perpendicular */
                n1x = ez; n1y = 0.0; n1z = -ex;
                n1len = sqrt(n1x * n1x + n1y * n1y + n1z * n1z);
            }
            if (n1len < 1e-15) continue;
            n1x /= n1len; n1y /= n1len; n1z /= n1len;

            /* Second perpendicular: cross of edge and n1 */
            double n2x = ey * n1z - ez * n1y;
            double n2y = ez * n1x - ex * n1z;
            double n2z = ex * n1y - ey * n1x;
            double n2len = sqrt(n2x * n2x + n2y * n2y + n2z * n2z);
            if (n2len < 1e-15) continue;
            n2x /= n2len; n2y /= n2len; n2z /= n2len;

            /* Midpoint of boundary edge */
            double mx = 0.5 * ((double)verts[v0].pos[0] + (double)verts[v1].pos[0]);
            double my = 0.5 * ((double)verts[v0].pos[1] + (double)verts[v1].pos[1]);
            double mz = 0.5 * ((double)verts[v0].pos[2] + (double)verts[v1].pos[2]);

            /* Plane 1 */
            double d1 = -(n1x * mx + n1y * my + n1z * mz);
            Quadric Qb1;
            quadric_from_plane(&Qb1, n1x, n1y, n1z, d1);
            {
                size_t qi = 0;
                for (qi = 0; qi < 10; qi++) Qb1.q[qi] *= QEM_BOUNDARY_WEIGHT;
            }
            quadric_add(&verts[v0].Q, &Qb1);
            quadric_add(&verts[v1].Q, &Qb1);

            /* Plane 2 */
            double d2 = -(n2x * mx + n2y * my + n2z * mz);
            Quadric Qb2;
            quadric_from_plane(&Qb2, n2x, n2y, n2z, d2);
            {
                size_t qi = 0;
                for (qi = 0; qi < 10; qi++) Qb2.q[qi] *= QEM_BOUNDARY_WEIGHT;
            }
            quadric_add(&verts[v0].Q, &Qb2);
            quadric_add(&verts[v1].Q, &Qb2);
        }
    }

    /* ---- Phase 4b: Auto-pin boundary vertices + caller pins ----
     * A vertex incident on any boundary edge (face_count == 1) is part
     * of an open boundary loop. Collapsing such loops sealed apertures
     * shut and created spurious k2 voids in the diagnostic run; the
     * boundary-penalty quadrics from Phase 4 prevent the boundary from
     * MOVING but do not prevent loop shrinkage via successive collapse.
     * Auto-pin makes boundary verts un-collapsible into each other and
     * forces interior collapses to keep the boundary vertex as the
     * survivor (see compute_edge_cost). The caller-supplied pin_mask
     * (halo seam verts, Dirichlet pins) is unioned in. */
    uint8_t *vert_pinned = (uint8_t *)ARENA_ALLOC(arena,
        (long)(in_nv * sizeof(uint8_t)));
    if (caller_pin_mask) {
        memcpy(vert_pinned, caller_pin_mask, in_nv * sizeof(uint8_t));
    } else {
        memset(vert_pinned, 0, in_nv * sizeof(uint8_t));
    }
    size_t n_caller_pins = 0;
    size_t n_boundary_verts = 0;
    {
        size_t i = 0;
        for (i = 0; i < in_nv; i++) {
            if (vert_pinned[i]) n_caller_pins++;
        }
        for (i = 0; i < n_edges; i++) {
            if (edges[i].face_count != 1) continue;
            if (!vert_pinned[edges[i].v0]) {
                vert_pinned[edges[i].v0] = 1;
                n_boundary_verts++;
            }
            if (!vert_pinned[edges[i].v1]) {
                vert_pinned[edges[i].v1] = 1;
                n_boundary_verts++;
            }
        }
    }
    (void)n_caller_pins;     /* available for debug logging */
    (void)n_boundary_verts;

    /* ---- Phase 5: Build adjacency + heap ---- */
    AdjList adj;
    adj.head = (int32_t *)ARENA_ALLOC(arena, (long)(in_nv * sizeof(int32_t)));
    {
        size_t i = 0;
        for (i = 0; i < in_nv; i++) adj.head[i] = -1;
    }
    /* Each edge appears in two adjacency lists, so pool needs 2*n_edges entries.
     * During collapse, we add more entries (merge adjacency). Budget extra. */
    adj.pool_cap = (int32_t)(n_edges * 4);
    adj.pool = (AdjNode *)ARENA_ALLOC(arena, (long)((size_t)adj.pool_cap * sizeof(AdjNode)));
    adj.pool_size = 0;

    {
        size_t i = 0;
        for (i = 0; i < n_edges; i++) {
            compute_edge_cost(&edges[i], verts, vert_pinned);
            adj_add(&adj, edges[i].v0, (int32_t)i);
            adj_add(&adj, edges[i].v1, (int32_t)i);
        }
    }

    /* Face adjacency: vertex -> incident faces (for normal-flip check) */
    FaceAdj fadj;
    fadj.head = (int32_t *)ARENA_ALLOC(arena, (long)(in_nv * sizeof(int32_t)));
    {
        size_t i = 0;
        for (i = 0; i < in_nv; i++) fadj.head[i] = -1;
    }
    fadj.pool_cap = (int32_t)(in_nf * 4);
    fadj.pool = (FaceAdjNode *)ARENA_ALLOC(arena,
        (long)((size_t)fadj.pool_cap * sizeof(FaceAdjNode)));
    fadj.pool_size = 0;
    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            fadj_add(&fadj, in_faces[f * 3 + 0], (int32_t)f);
            fadj_add(&fadj, in_faces[f * 3 + 1], (int32_t)f);
            fadj_add(&fadj, in_faces[f * 3 + 2], (int32_t)f);
        }
    }

    /* LME-rounds scratch arrays (see Phase 6 below). One allocation per
     * pass, reused round-to-round. */
    int32_t *best_edge_per_vert = (int32_t *)ARENA_ALLOC(arena,
        (long)(in_nv * sizeof(int32_t)));
    double  *best_cost_per_vert = (double *)ARENA_ALLOC(arena,
        (long)(in_nv * sizeof(double)));
    int32_t *lme_round = (int32_t *)ARENA_ALLOC(arena,
        (long)(n_edges * sizeof(int32_t)));
    uint8_t *vert_touched = (uint8_t *)ARENA_ALLOC(arena,
        (long)(in_nv * sizeof(uint8_t)));
    uint8_t *edge_recosted = (uint8_t *)ARENA_ALLOC(arena,
        (long)(n_edges * sizeof(uint8_t)));

    /* ---- Phase 5b: Compute median edge length for proximity radius ---- */
    float median_edge_len = 0.0f;
    {
        Arena_Mark scratch_mark = Arena_save(arena);
        size_t sample_count = n_edges < 1000 ? n_edges : 1000;
        float *edge_lens = (float *)ARENA_ALLOC(arena,
            (long)(sample_count * sizeof(float)));
        size_t step = n_edges / sample_count;
        if (step < 1) step = 1;
        size_t si = 0;
        size_t i = 0;
        for (i = 0; si < sample_count && i < n_edges; i += step) {
            float dx = verts[edges[i].v0].pos[0] - verts[edges[i].v1].pos[0];
            float dy = verts[edges[i].v0].pos[1] - verts[edges[i].v1].pos[1];
            float dz = verts[edges[i].v0].pos[2] - verts[edges[i].v1].pos[2];
            edge_lens[si] = sqrtf(dx*dx + dy*dy + dz*dz);
            si++;
        }
        if (si > 0) {
            qsort(edge_lens, si, sizeof(float), float_cmp);
            median_edge_len = edge_lens[si / 2];
        }
        Arena_restore(arena, scratch_mark);
    }
    float safe_radius = QEM_SAFE_RADIUS_FACTOR * median_edge_len;
    float safe_radius_sq = safe_radius * safe_radius;

    /* ---- Phase 5c: Build KD-tree for spatial proximity check ---- */
    Arena_Mark kd_mark = Arena_save(arena);
    float *kd_points = NULL;
    int32_t *kd_to_vert = NULL;
    KDTree_T prox_tree = NULL;
    size_t kd_alive_count = 0;
    prox_tree = build_prox_tree(arena, verts, in_nv,
                                 &kd_points, &kd_to_vert, &kd_alive_count);

    /* ---- Phase 6: LME-rounds collapse loop ----
     *
     * Each round:
     *   (1) Compute per-vertex "best (lowest-cost) incident edge".
     *   (2) Mark every edge that is the best edge for BOTH endpoints
     *       (reciprocal best) as an LME (Local Minimum Edge).
     *   (3) Collapse every LME in this round; LMEs by construction
     *       share no vertices, so they are non-conflicting.
     *   (4) Re-cost every edge incident to a touched vertex.
     *
     * The set of LMEs is uniquely determined by the local mesh — it
     * does not depend on edge-traversal order, no priority queue, no
     * tie-break bias. Two patches that share the 2-ring neighbourhood
     * of any seam edge therefore pick the same LMEs and collapse them
     * the same way (Ozaki & Kanai EGPGV 2015; Grund et al. VMV 2011).
     * That replaces our prior cross-cube quadric augmentation hack.
     */
    size_t current_nf = in_nf;
    size_t collapses = 0;
    size_t rejected = 0;
    size_t rejected_normal = 0;
    size_t rejected_prox = 0;
    size_t rejected_link = 0;
    int    round_idx = 0;

    while (current_nf > target_nf && round_idx < QEM_LME_MAX_ROUNDS) {

        /* (1) Per-vertex best (lowest-cost) incident edge. */
        {
            size_t v = 0;
            for (v = 0; v < in_nv; v++) {
                best_edge_per_vert[v] = -1;
                best_cost_per_vert[v] = DBL_MAX;
            }
        }
        {
            size_t ei = 0;
            for (ei = 0; ei < n_edges; ei++) {
                if (edges[ei].dead) continue;
                int32_t v0 = find_root(verts, edges[ei].v0);
                int32_t v1 = find_root(verts, edges[ei].v1);
                if (v0 == v1) { edges[ei].dead = 1; continue; }
                /* Staleness check: if our union-find chain has migrated
                 * either endpoint into a vertex whose pin state (or
                 * absorbed quadric) differs from the last time this
                 * edge was costed, the stored cost+opt_pos are wrong.
                 * Cheapest fix: update edges[].v0/v1 to current roots
                 * and re-call compute_edge_cost. */
                int32_t a = (v0 < v1) ? v0 : v1;
                int32_t b = (v0 < v1) ? v1 : v0;
                if (a != edges[ei].v0 || b != edges[ei].v1) {
                    edges[ei].v0 = a;
                    edges[ei].v1 = b;
                    compute_edge_cost(&edges[ei], verts, vert_pinned);
                    if (edges[ei].dead) continue;  /* now both-pinned -> dead */
                }
                double c = edges[ei].cost;
                if (c < best_cost_per_vert[v0]) {
                    best_cost_per_vert[v0] = c;
                    best_edge_per_vert[v0] = (int32_t)ei;
                }
                if (c < best_cost_per_vert[v1]) {
                    best_cost_per_vert[v1] = c;
                    best_edge_per_vert[v1] = (int32_t)ei;
                }
            }
        }

        /* (2) Identify LMEs (edges reciprocally best for both endpoints). */
        size_t n_lmes = 0;
        {
            size_t ei = 0;
            for (ei = 0; ei < n_edges; ei++) {
                if (edges[ei].dead) continue;
                int32_t v0 = find_root(verts, edges[ei].v0);
                int32_t v1 = find_root(verts, edges[ei].v1);
                if (best_edge_per_vert[v0] == (int32_t)ei
                    && best_edge_per_vert[v1] == (int32_t)ei) {
                    lme_round[n_lmes++] = (int32_t)ei;
                }
            }
        }
        if (n_lmes == 0) break;  /* no progress possible */

        /* (3) Collapse each LME. They are non-conflicting by construction,
         *     but successive collapses may invalidate later LMEs (e.g. via
         *     normals-oppose / flip / proximity guards), so the per-LME
         *     filters still apply per collapse. */
        memset(vert_touched, 0, in_nv * sizeof(uint8_t));
        size_t round_collapses = 0;
        {
            size_t k = 0;
            for (k = 0; k < n_lmes; k++) {
                int32_t ei = lme_round[k];
                if (edges[ei].dead) continue;

                int32_t v0 = find_root(verts, edges[ei].v0);
                int32_t v1 = find_root(verts, edges[ei].v1);
                if (v0 == v1 || !verts[v0].alive || !verts[v1].alive) {
                    edges[ei].dead = 1;
                    continue;
                }
                /* Belt-and-braces: LMEs share no vertex by construction, but
                 * find_root may have shifted that after earlier collapses
                 * this round. Skip rather than risk double-modifying a vert. */
                if (vert_touched[v0] || vert_touched[v1]) continue;

                if (normals_oppose(vert_normals, v0, v1)) {
                    edges[ei].dead = 1; rejected_normal++; continue;
                }
                if (collapse_would_flip(verts, in_faces, face_alive, &fadj,
                                        v0, v1, edges[ei].opt_pos)) {
                    edges[ei].dead = 1; rejected++; continue;
                }
                if (collapse_too_close(prox_tree, kd_to_vert,
                                        edges[ei].opt_pos, safe_radius_sq,
                                        verts, &fadj, in_faces, face_alive,
                                        v0, v1)) {
                    edges[ei].dead = 1; rejected_prox++; continue;
                }
                /* Topological guard: never collapse across a shared non-triangle
                 * neighbour (would create a non-manifold edge). See
                 * collapse_violates_link. */
                if (collapse_violates_link(verts, in_faces, face_alive, &fadj,
                                           v0, v1)) {
                    edges[ei].dead = 1; rejected_link++; continue;
                }

                /* Apply collapse v1 -> v0 (identical mechanics to prior
                 * heap-based handler, minus the heap_update/remove calls). */
                verts[v0].pos[0] = edges[ei].opt_pos[0];
                verts[v0].pos[1] = edges[ei].opt_pos[1];
                verts[v0].pos[2] = edges[ei].opt_pos[2];
                quadric_add(&verts[v0].Q, &verts[v1].Q);
                if (vert_pinned[v1]) vert_pinned[v0] = 1;
                verts[v1].alive = 0;
                verts[v1].remap = v0;
                edges[ei].dead  = 1;
                vert_touched[v0] = 1;
                vert_touched[v1] = 1;

                /* Update vertex normal for v0 (average + renormalize). */
                {
                    float nnx = vert_normals[v0*3+0] + vert_normals[v1*3+0];
                    float nny = vert_normals[v0*3+1] + vert_normals[v1*3+1];
                    float nnz = vert_normals[v0*3+2] + vert_normals[v1*3+2];
                    float nlen = sqrtf(nnx*nnx + nny*nny + nnz*nnz);
                    if (nlen > 1e-10f) {
                        vert_normals[v0*3+0] = nnx / nlen;
                        vert_normals[v0*3+1] = nny / nlen;
                        vert_normals[v0*3+2] = nnz / nlen;
                    }
                }

                /* Mark degenerate faces dead, decrement current_nf. */
                {
                    int32_t check_v[2];
                    int cvi = 0;
                    check_v[0] = v0;
                    check_v[1] = v1;
                    for (cvi = 0; cvi < 2; cvi++) {
                        int32_t fni = fadj.head[check_v[cvi]];
                        while (fni >= 0) {
                            int32_t fi = fadj.pool[fni].face_idx;
                            fni = fadj.pool[fni].next;
                            if (!face_alive[fi]) continue;
                            int32_t ri0 = find_root(verts, in_faces[fi * 3 + 0]);
                            int32_t ri1 = find_root(verts, in_faces[fi * 3 + 1]);
                            int32_t ri2 = find_root(verts, in_faces[fi * 3 + 2]);
                            if (ri0 == ri1 || ri1 == ri2 || ri0 == ri2) {
                                face_alive[fi] = 0;
                                current_nf--;
                            }
                        }
                    }
                }

                /* Splice fadj: v1's incident-face list onto v0. */
                if (fadj.head[v1] >= 0) {
                    int32_t tail = fadj.head[v1];
                    while (fadj.pool[tail].next >= 0) {
                        tail = fadj.pool[tail].next;
                    }
                    fadj.pool[tail].next = fadj.head[v0];
                    fadj.head[v0] = fadj.head[v1];
                    fadj.head[v1] = -1;
                }

                /* Splice adj (edge incidence): v1's edges onto v0. Critical
                 * for the end-of-round recost — without this splice, an
                 * edge whose original endpoints were (X, Y) both later
                 * merging into different chains is never visited via
                 * adj.head[current_root_of_X] (the edge is still hanging
                 * off adj.head[X], not adj.head[current_root]), so its
                 * cost/opt_pos go stale. Heap-based code papered over
                 * this with per-edge adj_add inside its merge loop; we
                 * just splice. */
                if (adj.head[v1] >= 0) {
                    int32_t tail = adj.head[v1];
                    while (adj.pool[tail].next >= 0) {
                        tail = adj.pool[tail].next;
                    }
                    adj.pool[tail].next = adj.head[v0];
                    adj.head[v0] = adj.head[v1];
                    adj.head[v1] = -1;
                }

                round_collapses++;
                collapses++;
            }
        }

        /* (4) Re-cost every edge incident to a touched vertex. Without a
         *     priority queue we just re-read the cost from edges[].cost in
         *     the next round; compute_edge_cost mutates that field in place. */
        memset(edge_recosted, 0, n_edges * sizeof(uint8_t));
        {
            size_t v = 0;
            for (v = 0; v < in_nv; v++) {
                if (!vert_touched[v]) continue;
                /* Walk both surviving v0's adj list AND dead v1's adj list.
                 * We DO NOT skip dead verts here — their adj entries still
                 * reference real edges, just with v1 now resolving via
                 * find_root to v0. */
                int32_t ni = adj.head[v];
                while (ni >= 0) {
                    int32_t edge_idx = adj.pool[ni].edge_idx;
                    ni = adj.pool[ni].next;
                    if (edges[edge_idx].dead) continue;
                    if (edge_recosted[edge_idx]) continue;
                    edge_recosted[edge_idx] = 1;

                    int32_t ea = find_root(verts, edges[edge_idx].v0);
                    int32_t eb = find_root(verts, edges[edge_idx].v1);
                    if (ea == eb) {
                        edges[edge_idx].dead = 1;
                        continue;
                    }
                    if (ea > eb) { int32_t t = ea; ea = eb; eb = t; }
                    edges[edge_idx].v0 = ea;
                    edges[edge_idx].v1 = eb;
                    compute_edge_cost(&edges[edge_idx], verts, vert_pinned);
                }
            }
        }

        /* Periodic KD-tree rebuild. Trigger when total collapses cross a
         * QEM_KDTREE_REBUILD_INTERVAL boundary this round. */
        if (prox_tree && round_collapses > 0) {
            size_t prev_total = collapses - round_collapses;
            if ((collapses / QEM_KDTREE_REBUILD_INTERVAL)
                > (prev_total / QEM_KDTREE_REBUILD_INTERVAL)) {
                Arena_restore(arena, kd_mark);
                kd_mark = Arena_save(arena);
                prox_tree = build_prox_tree(arena, verts, in_nv,
                                             &kd_points, &kd_to_vert,
                                             &kd_alive_count);
            }
        }

        round_idx++;

        /* Convergence: if this round collapsed fewer than the minimum
         * progress fraction of remaining faces, no point continuing. */
        if (round_collapses
            < (size_t)((double)current_nf * QEM_LME_MIN_PROGRESS_FRAC + 0.5)) {
            break;
        }
    }

    /* Free KD-tree before compaction (no longer needed) */
    Arena_restore(arena, kd_mark);

    /* ---- Phase 6b: Compute area-weighted reference normal for winding repair ---- */
    /* Cumulative float precision during collapses can produce a small number of
     * inverted faces (dot ≈ 0 passes but later flips). We compute a majority-vote
     * reference normal and fix any faces whose normal opposes it. */
    double ref_nx = 0.0, ref_ny = 0.0, ref_nz = 0.0;
    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            if (!face_alive[f]) continue;
            int32_t ri0 = find_root(verts, in_faces[f * 3 + 0]);
            int32_t ri1 = find_root(verts, in_faces[f * 3 + 1]);
            int32_t ri2 = find_root(verts, in_faces[f * 3 + 2]);
            if (ri0 == ri1 || ri1 == ri2 || ri0 == ri2) continue;
            /* Unnormalized cross product = area-weighted normal */
            float e1[3], e2[3];
            int k = 0;
            for (k = 0; k < 3; k++) {
                e1[k] = verts[ri1].pos[k] - verts[ri0].pos[k];
                e2[k] = verts[ri2].pos[k] - verts[ri0].pos[k];
            }
            ref_nx += (double)(e1[1]*e2[2] - e1[2]*e2[1]);
            ref_ny += (double)(e1[2]*e2[0] - e1[0]*e2[2]);
            ref_nz += (double)(e1[0]*e2[1] - e1[1]*e2[0]);
        }
    }

    /* ---- Phase 7: Compact ---- */

    /* Mark verts that appear in at least one non-degenerate live face.
     * Without this filter, boundary-pinned verts whose incident triangles
     * have all collapsed away survive into the output as dangling verts
     * (alive=1 but degree=0), which crashes downstream snap_cg via a
     * singular ILU diagonal. */
    uint8_t *vert_used = (uint8_t *)ARENA_ALLOC(arena,
        (long)(in_nv * sizeof(uint8_t)));
    memset(vert_used, 0, in_nv * sizeof(uint8_t));
    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            int32_t ri0 = find_root(verts, in_faces[f * 3 + 0]);
            int32_t ri1 = find_root(verts, in_faces[f * 3 + 1]);
            int32_t ri2 = find_root(verts, in_faces[f * 3 + 2]);
            if (ri0 == ri1 || ri1 == ri2 || ri0 == ri2) continue;
            vert_used[ri0] = 1;
            vert_used[ri1] = 1;
            vert_used[ri2] = 1;
        }
    }

    /* Build old-to-new vertex index map */
    int32_t *vert_map = (int32_t *)ARENA_ALLOC(arena, (long)(in_nv * sizeof(int32_t)));
    memset(vert_map, -1, in_nv * sizeof(int32_t));

    size_t new_nv = 0;
    {
        size_t i = 0;
        for (i = 0; i < in_nv; i++) {
            if (verts[i].alive && (vert_used[i] || vert_pinned[i])) {
                /* Keep a pinned vertex even if it became dangling after
                 * LME-rounds collapse — its position is the Dirichlet
                 * contract and must survive the simplification. The
                 * caller (raw_snap, grid_weld) treats isolated pinned
                 * verts as boundary samples. */
                vert_map[i] = (int32_t)new_nv;
                new_nv++;
            }
        }
    }

    /* Emit live, used vertices */
    float *new_verts = (float *)ARENA_ALLOC(arena, (long)(new_nv * 3 * sizeof(float)));
    {
        size_t i = 0;
        size_t vi = 0;
        for (i = 0; i < in_nv; i++) {
            if (verts[i].alive && (vert_used[i] || vert_pinned[i])) {
                new_verts[vi * 3 + 0] = verts[i].pos[0];
                new_verts[vi * 3 + 1] = verts[i].pos[1];
                new_verts[vi * 3 + 2] = verts[i].pos[2];
                vi++;
            }
        }
    }

    /* Emit non-degenerate faces */
    /* First pass: count */
    size_t new_nf = 0;
    {
        size_t f = 0;
        for (f = 0; f < in_nf; f++) {
            int32_t i0 = find_root(verts, in_faces[f * 3 + 0]);
            int32_t i1 = find_root(verts, in_faces[f * 3 + 1]);
            int32_t i2 = find_root(verts, in_faces[f * 3 + 2]);
            if (i0 == i1 || i1 == i2 || i0 == i2) continue;  /* degenerate */
            new_nf++;
        }
    }

    int32_t *new_faces = (int32_t *)ARENA_ALLOC(arena, (long)(new_nf * 3 * sizeof(int32_t)));
    size_t winding_fixes = 0;
    {
        size_t f = 0;
        size_t fi = 0;
        for (f = 0; f < in_nf; f++) {
            int32_t i0 = find_root(verts, in_faces[f * 3 + 0]);
            int32_t i1 = find_root(verts, in_faces[f * 3 + 1]);
            int32_t i2 = find_root(verts, in_faces[f * 3 + 2]);
            if (i0 == i1 || i1 == i2 || i0 == i2) continue;

            /* Check winding against reference normal; swap if inverted */
            {
                float fe1[3], fe2[3];
                int k = 0;
                for (k = 0; k < 3; k++) {
                    fe1[k] = verts[i1].pos[k] - verts[i0].pos[k];
                    fe2[k] = verts[i2].pos[k] - verts[i0].pos[k];
                }
                float fnx = fe1[1]*fe2[2] - fe1[2]*fe2[1];
                float fny = fe1[2]*fe2[0] - fe1[0]*fe2[2];
                float fnz = fe1[0]*fe2[1] - fe1[1]*fe2[0];
                double dot = (double)fnx * ref_nx + (double)fny * ref_ny
                           + (double)fnz * ref_nz;
                if (dot < 0.0) {
                    /* Swap i1 and i2 to reverse winding */
                    int32_t tmp = i1; i1 = i2; i2 = tmp;
                    winding_fixes++;
                }
            }

            assert(vert_map[i0] >= 0 && vert_map[i1] >= 0 && vert_map[i2] >= 0);
            new_faces[fi * 3 + 0] = vert_map[i0];
            new_faces[fi * 3 + 1] = vert_map[i1];
            new_faces[fi * 3 + 2] = vert_map[i2];
            fi++;
        }
        assert(fi == new_nf);
    }

#ifdef VESUVIUS_DEBUG
    fprintf(stderr, "    QEM: %zu -> %zu faces (%zu -> %zu verts, %zu collapses, "
            "rej: %zu flip, %zu normal, %zu prox, %zu link, %zu winding fixes)\n",
            in_nf, new_nf, in_nv, new_nv, collapses,
            rejected, rejected_normal, rejected_prox, rejected_link, winding_fixes);
#endif

    *out_verts = new_verts;
    *out_nv = new_nv;
    *out_faces = new_faces;
    *out_nf = new_nf;

    /* Emit compacted pin_mask if requested. vert_pinned has been updated
     * during the collapse loop to propagate pins through one-pinned
     * collapses (see line 1494: `if (vert_pinned[v1]) vert_pinned[v0]=1`).
     * The compaction map vert_map[] turns old vert indices into new ones. */
    if (out_pin_mask) {
        uint8_t *pout = (uint8_t *)ARENA_CALLOC(arena,
            (long)new_nv, 1L);
        for (size_t i = 0; i < in_nv; i++) {
            if (vert_map[i] >= 0 && vert_pinned[i]) {
                pout[vert_map[i]] = 1;
            }
        }
        *out_pin_mask = pout;
    }

    return 0;
}

/* ================================================================
 * Public QEM simplification.
 *
 * Single collapse pass to target, then one edge-flip + tangential-
 * smooth maintenance pass on the result to improve triangle quality.
 * ================================================================ */

int QEM_simplify(Arena_T arena,
                 const float *in_verts, size_t in_nv,
                 const int32_t *in_faces, size_t in_nf,
                 size_t target_nf,
                 float **out_verts, size_t *out_nv,
                 int32_t **out_faces, size_t *out_nf)
{
    return QEM_simplify_pinned(arena, in_verts, in_nv, in_faces, in_nf,
                               NULL, NULL, target_nf,
                               out_verts, out_nv, out_faces, out_nf,
                               NULL);
}

int QEM_simplify_pinned(Arena_T arena,
                        const float *in_verts, size_t in_nv,
                        const int32_t *in_faces, size_t in_nf,
                        const uint8_t *pin_mask,
                        const float *ref_normal,
                        size_t target_nf,
                        float **out_verts, size_t *out_nv,
                        int32_t **out_faces, size_t *out_nf,
                        uint8_t **out_pin_mask)
{
    int rc = 0;

    assert(arena);
    assert(out_verts && out_nv && out_faces && out_nf);

    /* Handle trivial cases */
    if (in_nv == 0 || in_nf == 0) {
        *out_verts = NULL;
        *out_nv = 0;
        *out_faces = NULL;
        *out_nf = 0;
        if (out_pin_mask) *out_pin_mask = NULL;
        return 0;
    }

    if (in_nf <= target_nf) {
        float *vout = (float *)ARENA_ALLOC(arena,
            (long)(in_nv * 3 * sizeof(float)));
        memcpy(vout, in_verts, in_nv * 3 * sizeof(float));
        int32_t *fout = (int32_t *)ARENA_ALLOC(arena,
            (long)(in_nf * 3 * sizeof(int32_t)));
        memcpy(fout, in_faces, in_nf * 3 * sizeof(int32_t));
        *out_verts = vout;
        *out_nv = in_nv;
        *out_faces = fout;
        *out_nf = in_nf;
        if (out_pin_mask) {
            if (pin_mask) {
                uint8_t *pout = (uint8_t *)ARENA_ALLOC(arena,
                    (long)(in_nv * sizeof(uint8_t)));
                memcpy(pout, pin_mask, in_nv * sizeof(uint8_t));
                *out_pin_mask = pout;
            } else {
                *out_pin_mask = NULL;
            }
        }
        return 0;
    }

    /* Single collapse pass all the way to target */
    rc = qem_collapse_pass(arena, in_verts, in_nv, in_faces, in_nf,
                            pin_mask, target_nf,
                            out_verts, out_nv, out_faces, out_nf,
                            out_pin_mask);
    if (rc != 0) return rc;

    /* One maintenance pass (edge flips + smooth) on the final mesh.
     * Pin-aware: qem_smooth_pass skips pinned verts so the Dirichlet
     * contract is preserved through maintenance. */
    if (*out_nf > 100) {
        const uint8_t *post_pins = (out_pin_mask && *out_pin_mask)
                                       ? *out_pin_mask : NULL;
        qem_maintenance_pass(arena, *out_verts, *out_nv,
                              *out_faces, *out_nf,
                              post_pins);
    }

    /* ---- Final winding repair ----
     * Both the collapse pass and maintenance pass can leave residual
     * inverted faces from float precision edge cases. Fix them by
     * computing an area-weighted reference normal and flipping any
     * face whose normal opposes it. Iterate until convergence since
     * fixing faces shifts the reference normal slightly. O(nf) per
     * iteration, typically converges in 1-2 passes. */
    if (*out_nf > 0) {
        float *fv = *out_verts;
        int32_t *ff = *out_faces;
        size_t fnf = *out_nf;
        size_t total_fixed = 0;
        int iter = 0;
        int max_iter = 5;  /* safety limit */

        for (iter = 0; iter < max_iter; iter++) {
            double rnx = 0.0, rny = 0.0, rnz = 0.0;
            size_t fi = 0;
            size_t n_fixed = 0;

            if (ref_normal) {
                /* Caller supplied reference (e.g. cm->pca_normal which
                 * is (1,1,1)-hemisphere oriented per pca.c so adjacent
                 * cubes share the same sign). Use it directly so the
                 * face flips are cross-cube deterministic. */
                rnx = (double)ref_normal[0];
                rny = (double)ref_normal[1];
                rnz = (double)ref_normal[2];
            } else {
                /* Accumulate area-weighted reference normal from the
                 * post-collapse mesh. This is per-component and so can
                 * point in opposite directions across adjacent cubes
                 * whose meshes are different subsets of the same
                 * underlying surface — see ref_normal docs in qem.h. */
                for (fi = 0; fi < fnf; fi++) {
                    int32_t a = ff[fi * 3 + 0];
                    int32_t b = ff[fi * 3 + 1];
                    int32_t c = ff[fi * 3 + 2];
                    float e1[3], e2[3];
                    int k = 0;
                    for (k = 0; k < 3; k++) {
                        e1[k] = fv[b * 3 + k] - fv[a * 3 + k];
                        e2[k] = fv[c * 3 + k] - fv[a * 3 + k];
                    }
                    rnx += (double)(e1[1]*e2[2] - e1[2]*e2[1]);
                    rny += (double)(e1[2]*e2[0] - e1[0]*e2[2]);
                    rnz += (double)(e1[0]*e2[1] - e1[1]*e2[0]);
                }
            }

            /* Fix any face opposing the reference */
            for (fi = 0; fi < fnf; fi++) {
                int32_t a = ff[fi * 3 + 0];
                int32_t b = ff[fi * 3 + 1];
                int32_t c = ff[fi * 3 + 2];
                float e1[3], e2[3];
                int k = 0;
                for (k = 0; k < 3; k++) {
                    e1[k] = fv[b * 3 + k] - fv[a * 3 + k];
                    e2[k] = fv[c * 3 + k] - fv[a * 3 + k];
                }
                float fnx = e1[1]*e2[2] - e1[2]*e2[1];
                float fny = e1[2]*e2[0] - e1[0]*e2[2];
                float fnz = e1[0]*e2[1] - e1[1]*e2[0];
                double dot = (double)fnx * rnx + (double)fny * rny
                           + (double)fnz * rnz;
                if (dot < 0.0) {
                    ff[fi * 3 + 1] = c;
                    ff[fi * 3 + 2] = b;
                    n_fixed++;
                }
            }

            total_fixed += n_fixed;
            if (n_fixed == 0) break;  /* converged */
        }

#ifdef VESUVIUS_DEBUG
        if (total_fixed > 0) {
            fprintf(stderr, "    QEM final winding repair: %zu faces fixed "
                    "(%d iterations)\n", total_fixed, iter);
        }
#endif
    }

    return 0;
}
