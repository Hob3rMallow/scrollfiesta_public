#include "patch_repair.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Local tangent-plane Delaunay re-triangulation of BPA tear patches.
 * See patch_repair.h for the contract and invariants.
 *
 * Implementation notes:
 *  - All scratch is arena-allocated; per the Weld_verts convention the arena
 *    is NOT restored (out_faces lives on it; caller's lifetime subsumes us).
 *  - Deterministic: every loop iterates in ascending index order; the ear
 *    clipper takes the first valid ear in walk order; interior points are
 *    inserted in ascending vertex index; Lawson runs full passes to a
 *    fixpoint.
 *  - int32_t face/vert indices internally (a per-component mesh is far below
 *    2^31); size_t only at the API boundary.
 *  - Safety posture: every ambiguous situation skips the patch (mesh is left
 *    untouched there). A missed repair is a hole downstream hole-fill may
 *    still take; a wrong repair is a topology error. We only take wins.
 */

/* Tuning constants — mirrored in pipeline_constants.h (PR_*). Local copies
 * per the ball_pivot.c convention. */
#define PR_TEAR_MAX_EDGES     64    /* boundary cluster larger than this is not a tear */
#define PR_PATCH_RINGS         2    /* face rings grown around a tear                   */
#define PR_PATCH_MAX_FACES   256    /* patch bigger than this is skipped                */
#define PR_PATCH_MAX_VERTS   512
#define PR_PLANAR_MAX_THICK 0.75    /* sqrt(smallest PCA eigenvalue) cap, vox           */
#define PR_EPS_2D          1e-5     /* 2D degeneracy epsilon (orient / on-edge)         */
#define PR_MAX_PATCHES      4096    /* safety cap on tear clusters processed            */
#define PR_LAWSON_MAX_PASSES 64     /* full no-flip-pass convergence cap                */

/* ------------------------------------------------------------------ */
/* small helpers                                                       */
/* ------------------------------------------------------------------ */

static int32_t uf_find(int32_t *p, int32_t x)
{
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}
static void uf_union(int32_t *p, int32_t a, int32_t b)
{
    int32_t ra = uf_find(p, a), rb = uf_find(p, b);
    if (ra == rb) return;
    if (ra < rb) p[rb] = ra; else p[ra] = rb;   /* smallest index is root */
}

/* 3x3 symmetric Jacobi eigensolver (ascending eigenvalues). Local mirror of
 * src/common/mls_project.c:jacobi_3x3 — kept private per repo convention. */
static void jacobi_3x3(double a[3][3], double eigenvalues[3],
                       double eigenvectors[3][3])
{
    int i = 0, j = 0;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
    double m[3][3];
    memcpy(m, a, sizeof(m));
    for (int iter = 0; iter < 100; iter++) {
        double max_off = 0.0;
        int p = 0, q = 1;
        for (i = 0; i < 3; i++)
            for (j = i + 1; j < 3; j++)
                if (fabs(m[i][j]) > max_off) { max_off = fabs(m[i][j]); p = i; q = j; }
        if (max_off < 1e-15) break;
        double theta = 0.0;
        double diff = m[q][q] - m[p][p];
        if (fabs(diff) < 1e-30) theta = 3.14159265358979323846 / 4.0;
        else theta = 0.5 * atan2(2.0 * m[p][q], diff);
        double c = cos(theta), s = sin(theta);
        double app = m[p][p], aqq = m[q][q], apq = m[p][q];
        m[p][p] = c*c*app - 2.0*s*c*apq + s*s*aqq;
        m[q][q] = s*s*app + 2.0*s*c*apq + c*c*aqq;
        m[p][q] = 0.0; m[q][p] = 0.0;
        for (int r = 0; r < 3; r++) {
            if (r == p || r == q) continue;
            double rp = m[r][p], rq = m[r][q];
            m[r][p] = c*rp - s*rq; m[p][r] = m[r][p];
            m[r][q] = s*rp + c*rq; m[q][r] = m[r][q];
        }
        for (int r = 0; r < 3; r++) {
            double vp = eigenvectors[r][p], vq = eigenvectors[r][q];
            eigenvectors[r][p] = c*vp - s*vq;
            eigenvectors[r][q] = s*vp + c*vq;
        }
    }
    double evals[3] = { m[0][0], m[1][1], m[2][2] };
    int order[3] = { 0, 1, 2 };
    for (i = 1; i < 3; i++) {
        int k = i;
        while (k > 0 && evals[order[k-1]] > evals[order[k]]) {
            int t = order[k-1]; order[k-1] = order[k]; order[k] = t; k--;
        }
    }
    double sorted[3][3];
    for (i = 0; i < 3; i++) {
        eigenvalues[i] = evals[order[i]];
        for (j = 0; j < 3; j++) sorted[j][i] = eigenvectors[j][order[i]];
    }
    memcpy(eigenvectors, sorted, sizeof(sorted));
}

/* ------------------------------------------------------------------ */
/* unique-edge table                                                   */
/* ------------------------------------------------------------------ */

typedef struct { int32_t u, v; int32_t f; } ERec;          /* one face-edge  */
typedef struct { int32_t u, v; int32_t f0, f1; int32_t mult; } UEdge;

static int cmp_erec(const void *a, const void *b)
{
    const ERec *ea = (const ERec *)a, *eb = (const ERec *)b;
    if (ea->u != eb->u) return (ea->u < eb->u) ? -1 : 1;
    if (ea->v != eb->v) return (ea->v < eb->v) ? -1 : 1;
    if (ea->f != eb->f) return (ea->f < eb->f) ? -1 : 1;
    return 0;
}

/* binary search UE for (u,v), u<v. Returns index or -1. */
static int32_t ue_find(const UEdge *ue, int32_t n_ue, int32_t u, int32_t v)
{
    int32_t lo = 0, hi = n_ue;
    while (lo < hi) {
        int32_t mid = lo + (hi - lo) / 2;
        if (ue[mid].u < u || (ue[mid].u == u && ue[mid].v < v)) lo = mid + 1;
        else hi = mid;
    }
    if (lo < n_ue && ue[lo].u == u && ue[lo].v == v) return lo;
    return -1;
}

/* ------------------------------------------------------------------ */
/* 2D predicates                                                       */
/* ------------------------------------------------------------------ */

static double orient2d(const double *a, const double *b, const double *c)
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}

/* strict segment crossing (shared endpoints excluded by the caller) */
static int segs_cross(const double *a, const double *b,
                      const double *c, const double *d)
{
    double d1 = orient2d(c, d, a), d2 = orient2d(c, d, b);
    double d3 = orient2d(a, b, c), d4 = orient2d(a, b, d);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
        ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) return 1;
    return 0;
}

/* point strictly inside CCW triangle when eps > 0; inclusive when eps < 0 */
static int pt_in_tri(const double *p, const double *a, const double *b,
                     const double *c, double eps)
{
    double o1 = orient2d(a, b, p);
    double o2 = orient2d(b, c, p);
    double o3 = orient2d(c, a, p);
    return (o1 > eps && o2 > eps && o3 > eps);
}

/* squared distance from p to segment ab */
static double pt_seg_d2(const double *p, const double *a, const double *b)
{
    double abx = b[0]-a[0], aby = b[1]-a[1];
    double apx = p[0]-a[0], apy = p[1]-a[1];
    double L2 = abx*abx + aby*aby;
    double t = (L2 > 1e-30) ? (apx*abx + apy*aby) / L2 : 0.0;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    double dx = apx - t*abx, dy = apy - t*aby;
    return dx*dx + dy*dy;
}

/* incircle for Lawson: (a,b,c) CCW; > 0 <=> d strictly inside circumcircle.
 * Plain double is ample at patch scale (coords O(10) vox, spacing ~1 vox). */
static int incircle(const double *a, const double *b, const double *c,
                    const double *d)
{
    double ax = a[0]-d[0], ay = a[1]-d[1];
    double bx = b[0]-d[0], by = b[1]-d[1];
    double cx = c[0]-d[0], cy = c[1]-d[1];
    double det = (ax*ax + ay*ay) * (bx*cy - cx*by)
               - (bx*bx + by*by) * (ax*cy - cx*ay)
               + (cx*cx + cy*cy) * (ax*by - bx*ay);
    return det > 1e-12;
}

/* ------------------------------------------------------------------ */
/* per-patch triangulation (local indices 0..npv-1)                    */
/* ------------------------------------------------------------------ */

typedef struct { int32_t v[3]; } LTri;

/* Ear-clip a CCW simple polygon ring[0..nr-1] (local indices into P[2*i]).
 * Appends nr-2 CCW triangles to T. `work` is caller scratch >= nr ints. */
static int ear_clip(const double *P, const int32_t *ring_in, int32_t nr,
                    LTri *T, int32_t *n_tri, int32_t *work)
{
    if (nr < 3) return -1;
    int32_t *ring = work;
    memcpy(ring, ring_in, (size_t)nr * sizeof(int32_t));
    int32_t n = nr;
    int32_t guard = nr * nr + 8;
    while (n > 3 && guard-- > 0) {
        int clipped = 0;
        for (int32_t i = 0; i < n; i++) {
            int32_t ip = (i + n - 1) % n, in = (i + 1) % n;
            const double *a = &P[2*ring[ip]], *b = &P[2*ring[i]], *c = &P[2*ring[in]];
            if (orient2d(a, b, c) <= PR_EPS_2D) continue;     /* reflex/flat */
            int blocked = 0;
            for (int32_t k = 0; k < n && !blocked; k++) {
                if (k == ip || k == i || k == in) continue;
                if (pt_in_tri(&P[2*ring[k]], a, b, c, -PR_EPS_2D)) blocked = 1;
            }
            if (blocked) continue;
            T[*n_tri].v[0] = ring[ip];
            T[*n_tri].v[1] = ring[i];
            T[*n_tri].v[2] = ring[in];
            (*n_tri)++;
            for (int32_t k = i; k < n - 1; k++) ring[k] = ring[k+1];
            n--;
            clipped = 1;
            break;
        }
        if (!clipped) return -1;                    /* stuck (degenerate)    */
    }
    if (n == 3) {
        if (orient2d(&P[2*ring[0]], &P[2*ring[1]], &P[2*ring[2]]) <= PR_EPS_2D)
            return -1;
        T[*n_tri].v[0] = ring[0];
        T[*n_tri].v[1] = ring[1];
        T[*n_tri].v[2] = ring[2];
        (*n_tri)++;
    }
    return (guard > 0) ? 0 : -1;
}

/* Insert point p into the triangulation. Strictly inside a triangle =>
 * 1->3 split. Exactly on an INTERIOR edge (both adjacent triangles present)
 * => 2->4 split — required for lattice-regular data where points fall
 * exactly on diagonals. On a ring edge / not found => -1 (caller skips;
 * ring-edge proximity was already excluded by the pt_seg_d2 gate).
 * Every successful insertion adds exactly 2 triangles (Euler holds). */
static int insert_point(const double *P, int32_t p, LTri *T, int32_t *n_tri)
{
    /* pass 1: strict containment */
    for (int32_t t = 0; t < *n_tri; t++) {
        const double *a = &P[2*T[t].v[0]], *b = &P[2*T[t].v[1]], *c = &P[2*T[t].v[2]];
        if (!pt_in_tri(&P[2*p], a, b, c, PR_EPS_2D)) continue;
        int32_t va = T[t].v[0], vb = T[t].v[1], vc = T[t].v[2];
        T[t].v[0] = va; T[t].v[1] = vb; T[t].v[2] = p;
        T[*n_tri].v[0] = vb; T[*n_tri].v[1] = vc; T[*n_tri].v[2] = p; (*n_tri)++;
        T[*n_tri].v[0] = vc; T[*n_tri].v[1] = va; T[*n_tri].v[2] = p; (*n_tri)++;
        return 0;
    }
    /* pass 2: on an interior edge — locate the edge and both triangles */
    for (int32_t t = 0; t < *n_tri; t++) {
        if (!pt_in_tri(&P[2*p], &P[2*T[t].v[0]], &P[2*T[t].v[1]], &P[2*T[t].v[2]],
                       -PR_EPS_2D)) continue;              /* inclusive miss  */
        for (int e = 0; e < 3; e++) {
            int32_t va = T[t].v[e], vb = T[t].v[(e+1)%3], vc = T[t].v[(e+2)%3];
            if (fabs(orient2d(&P[2*va], &P[2*vb], &P[2*p])) > PR_EPS_2D) continue;
            /* p collinear with (va,vb); must lie between them */
            if (pt_seg_d2(&P[2*p], &P[2*va], &P[2*vb]) > PR_EPS_2D*PR_EPS_2D)
                continue;
            /* find the other triangle carrying (vb,va) */
            int32_t t2 = -1, vd = -1;
            for (int32_t k = 0; k < *n_tri && t2 < 0; k++) {
                if (k == t) continue;
                for (int e2 = 0; e2 < 3; e2++) {
                    int32_t x = T[k].v[e2], y = T[k].v[(e2+1)%3];
                    if ((x == vb && y == va) || (x == va && y == vb)) {
                        t2 = k; vd = T[k].v[(e2+2)%3]; break;
                    }
                }
            }
            if (t2 < 0) return -1;              /* ring edge: refuse         */
            /* split t: (va,vb,vc) -> (va,p,vc) + (p,vb,vc)                  */
            T[t].v[0] = va; T[t].v[1] = p; T[t].v[2] = vc;
            T[*n_tri].v[0] = p; T[*n_tri].v[1] = vb; T[*n_tri].v[2] = vc; (*n_tri)++;
            /* split t2: contains (vb,va) or (va,vb); rewrite as fan around p */
            /* determine t2's winding of the shared edge */
            int rev = 0;                        /* t2 uses (vb,va)?          */
            for (int e2 = 0; e2 < 3; e2++) {
                int32_t x = T[t2].v[e2], y = T[t2].v[(e2+1)%3];
                if (x == vb && y == va) { rev = 1; break; }
                if (x == va && y == vb) { rev = 0; break; }
            }
            if (rev) {                          /* (vb,va,vd) CCW            */
                T[t2].v[0] = vb; T[t2].v[1] = p; T[t2].v[2] = vd;
                T[*n_tri].v[0] = p; T[*n_tri].v[1] = va; T[*n_tri].v[2] = vd; (*n_tri)++;
            } else {                            /* (va,vb,vd) CCW            */
                T[t2].v[0] = va; T[t2].v[1] = p; T[t2].v[2] = vd;
                T[*n_tri].v[0] = p; T[*n_tri].v[1] = vb; T[*n_tri].v[2] = vd; (*n_tri)++;
            }
            return 0;
        }
    }
    return -1;
}

/* Lawson flips, full passes to a fixpoint (no restart-per-flip). Ring edges
 * have only one triangle in this triangulation, so they can never flip.
 * Returns 0 on convergence, -1 on pass-cap (treated as failure). */
static int lawson_flips(const double *P, LTri *T, int32_t n_tri)
{
    for (int pass = 0; pass < PR_LAWSON_MAX_PASSES; pass++) {
        int32_t flips = 0;
        for (int32_t t1 = 0; t1 < n_tri; t1++) {
            for (int e1 = 0; e1 < 3; e1++) {
                int32_t a = T[t1].v[e1], b = T[t1].v[(e1+1)%3];
                int32_t c = T[t1].v[(e1+2)%3];
                int32_t t2f = -1, d = -1;
                for (int32_t t2 = t1 + 1; t2 < n_tri && t2f < 0; t2++) {
                    for (int e2 = 0; e2 < 3; e2++) {
                        int32_t x = T[t2].v[e2], y = T[t2].v[(e2+1)%3];
                        if ((x == b && y == a) || (x == a && y == b)) {
                            t2f = t2;
                            d = T[t2].v[(e2+2)%3];
                            break;
                        }
                    }
                }
                if (t2f < 0) continue;                  /* boundary edge      */
                if (!incircle(&P[2*a], &P[2*b], &P[2*c], &P[2*d])) continue;
                /* flip (a,b) -> (c,d), keep both CCW; skip degenerate flips */
                if (orient2d(&P[2*a], &P[2*d], &P[2*c]) <= PR_EPS_2D) continue;
                if (orient2d(&P[2*d], &P[2*b], &P[2*c]) <= PR_EPS_2D) continue;
                T[t1].v[0] = a; T[t1].v[1] = d; T[t1].v[2] = c;
                T[t2f].v[0] = d; T[t2f].v[1] = b; T[t2f].v[2] = c;
                flips++;
            }
        }
        if (flips == 0) return 0;
    }
    return -1;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int PatchRepair_repair(Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       int32_t **out_faces, size_t *out_nf,
                       PatchRepairStats *stats)
{
    assert(arena);
    assert(out_faces);
    assert(out_nf);

    PatchRepairStats st;
    memset(&st, 0, sizeof st);

    *out_faces = NULL;
    *out_nf = 0;

    /* trivial inputs: hand back a copy (contract: out is arena-owned) */
    if (nf == 0 || nv < 3 || !verts || !faces) {
        int32_t *cp = (int32_t *)ARENA_ALLOC(arena, (long)((nf ? nf : 1) * 3 * sizeof(int32_t)));
        if (nf) memcpy(cp, faces, nf * 3 * sizeof(int32_t));
        *out_faces = cp;
        *out_nf = nf;
        if (stats) *stats = st;
        return 0;
    }
    assert(nv <= 0x7FFFFFFF && nf <= 0x7FFFFFFF);
    int32_t NV = (int32_t)nv, NF = (int32_t)nf;

    /* ---- 1. unique-edge table ---- */
    ERec *er = (ERec *)ARENA_ALLOC(arena, (long)NF * 3L * (long)sizeof(ERec));
    for (int32_t f = 0; f < NF; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t u = faces[(size_t)f*3 + (size_t)k];
            int32_t v = faces[(size_t)f*3 + (size_t)((k+1)%3)];
            if (u > v) { int32_t t = u; u = v; v = t; }
            er[f*3+k].u = u; er[f*3+k].v = v; er[f*3+k].f = f;
        }
    }
    qsort(er, (size_t)NF*3, sizeof(ERec), cmp_erec);
    UEdge *ue = (UEdge *)ARENA_ALLOC(arena, (long)NF * 3L * (long)sizeof(UEdge));
    int32_t n_ue = 0;
    for (int32_t i = 0; i < NF*3; ) {
        int32_t j = i;
        while (j < NF*3 && er[j].u == er[i].u && er[j].v == er[i].v) j++;
        ue[n_ue].u = er[i].u; ue[n_ue].v = er[i].v;
        ue[n_ue].mult = j - i;
        ue[n_ue].f0 = er[i].f;
        ue[n_ue].f1 = (j - i > 1) ? er[i+1].f : -1;
        n_ue++;
        i = j;
    }

    /* ---- 2. face-connectivity components (over shared edges) ---- */
    int32_t *vcomp = (int32_t *)ARENA_ALLOC(arena, (long)NV * (long)sizeof(int32_t));
    for (int32_t i = 0; i < NV; i++) vcomp[i] = i;
    for (int32_t e = 0; e < n_ue; e++) uf_union(vcomp, ue[e].u, ue[e].v);
    for (int32_t i = 0; i < NV; i++) vcomp[i] = uf_find(vcomp, i);
    {
        uint8_t *seen = (uint8_t *)ARENA_CALLOC(arena, (long)NV, 1L);
        for (int32_t f = 0; f < NF; f++) {
            int32_t r = vcomp[faces[(size_t)f*3]];
            if (!seen[r]) { seen[r] = 1; st.n_components++; }
        }
    }

    /* ---- 3. boundary clusters (UF over boundary-edge endpoints) ---- */
    int32_t *bcl = (int32_t *)ARENA_ALLOC(arena, (long)NV * (long)sizeof(int32_t));
    for (int32_t i = 0; i < NV; i++) bcl[i] = i;
    for (int32_t e = 0; e < n_ue; e++)
        if (ue[e].mult == 1) uf_union(bcl, ue[e].u, ue[e].v);
    int32_t *cl_edges = (int32_t *)ARENA_CALLOC(arena, (long)NV, (long)sizeof(int32_t));
    for (int32_t e = 0; e < n_ue; e++)
        if (ue[e].mult == 1) cl_edges[uf_find(bcl, ue[e].u)]++;
    /* per component: perimeter cluster = the one with the most edges
     * (ties: smallest root, for determinism) */
    int32_t *comp_perim = (int32_t *)ARENA_ALLOC(arena, (long)NV * (long)sizeof(int32_t));
    int32_t *comp_perim_n = (int32_t *)ARENA_CALLOC(arena, (long)NV, (long)sizeof(int32_t));
    for (int32_t i = 0; i < NV; i++) comp_perim[i] = -1;
    for (int32_t e = 0; e < n_ue; e++) {
        if (ue[e].mult != 1) continue;
        int32_t cl = uf_find(bcl, ue[e].u);
        int32_t cp = vcomp[ue[e].u];
        if (cl_edges[cl] > comp_perim_n[cp] ||
            (cl_edges[cl] == comp_perim_n[cp] && (comp_perim[cp] < 0 || cl < comp_perim[cp]))) {
            comp_perim[cp] = cl;
            comp_perim_n[cp] = cl_edges[cl];
        }
    }

    /* ---- 4. tear clusters ---- */
    int32_t *tear_root = (int32_t *)ARENA_ALLOC(arena, (long)PR_MAX_PATCHES * (long)sizeof(int32_t));
    int32_t n_tear = 0;
    for (int32_t r = 0; r < NV && n_tear < PR_MAX_PATCHES; r++) {
        if (cl_edges[r] == 0) continue;
        if (uf_find(bcl, r) != r) continue;
        if (cl_edges[r] > PR_TEAR_MAX_EDGES) continue;
        if (comp_perim[vcomp[r]] == r) continue;        /* the perimeter      */
        tear_root[n_tear++] = r;
    }
    st.n_tear_clusters = (size_t)n_tear;

    if (n_tear == 0) {
        int32_t *cp = (int32_t *)ARENA_ALLOC(arena, (long)NF * 3L * (long)sizeof(int32_t));
        memcpy(cp, faces, (size_t)NF * 3 * sizeof(int32_t));
        *out_faces = cp;
        *out_nf = (size_t)NF;
        if (stats) *stats = st;
        return 0;
    }

    /* tear-cluster membership: is boundary-cluster root r a tear root?     */
    uint8_t *is_tear_root = (uint8_t *)ARENA_CALLOC(arena, (long)NV, 1L);
    for (int32_t t = 0; t < n_tear; t++) is_tear_root[tear_root[t]] = 1;

    /* ---- 5. vert -> face incidence (CSR) ---- */
    int32_t *vf_cnt = (int32_t *)ARENA_CALLOC(arena, (long)(NV + 1), (long)sizeof(int32_t));
    for (int32_t f = 0; f < NF; f++)
        for (int k = 0; k < 3; k++) vf_cnt[faces[(size_t)f*3+(size_t)k] + 1]++;
    for (int32_t i = 0; i < NV; i++) vf_cnt[i+1] += vf_cnt[i];
    int32_t *vf = (int32_t *)ARENA_ALLOC(arena, (long)NF * 3L * (long)sizeof(int32_t));
    {
        int32_t *cur = (int32_t *)ARENA_ALLOC(arena, (long)NV * (long)sizeof(int32_t));
        memcpy(cur, vf_cnt, (size_t)NV * sizeof(int32_t));
        for (int32_t f = 0; f < NF; f++)
            for (int k = 0; k < 3; k++) vf[cur[faces[(size_t)f*3+(size_t)k]]++] = f;
    }

    /* ---- 6. grow K-ring patches per tear; merge overlaps ---- */
    int32_t *fpatch = (int32_t *)ARENA_ALLOC(arena, (long)NF * (long)sizeof(int32_t));
    for (int32_t f = 0; f < NF; f++) fpatch[f] = -1;
    int32_t *puf = (int32_t *)ARENA_ALLOC(arena, (long)n_tear * (long)sizeof(int32_t));
    for (int32_t t = 0; t < n_tear; t++) puf[t] = t;
    uint8_t *vmark = (uint8_t *)ARENA_CALLOC(arena, (long)NV, 1L);
    int32_t *vlist = (int32_t *)ARENA_ALLOC(arena, (long)NV * (long)sizeof(int32_t));
    int32_t *pf_all = (int32_t *)ARENA_ALLOC(arena, (long)NF * (long)sizeof(int32_t));
    int32_t pf_used = 0;
    int32_t *pf_beg = (int32_t *)ARENA_ALLOC(arena, (long)(n_tear + 1) * (long)sizeof(int32_t));

    for (int32_t t = 0; t < n_tear; t++) {
        pf_beg[t] = pf_used;
        int32_t nvl = 0;
        for (int32_t e = 0; e < n_ue; e++) {
            if (ue[e].mult != 1) continue;
            if (uf_find(bcl, ue[e].u) != tear_root[t]) continue;
            if (!vmark[ue[e].u]) { vmark[ue[e].u] = 1; vlist[nvl++] = ue[e].u; }
            if (!vmark[ue[e].v]) { vmark[ue[e].v] = 1; vlist[nvl++] = ue[e].v; }
        }
        int32_t ring_lo = 0;
        int overflow = 0;
        for (int ring = 0; ring < PR_PATCH_RINGS && !overflow; ring++) {
            int32_t ring_hi = nvl;
            for (int32_t vi = ring_lo; vi < ring_hi && !overflow; vi++) {
                int32_t wv2 = vlist[vi];
                for (int32_t k = vf_cnt[wv2]; k < vf_cnt[wv2+1]; k++) {
                    int32_t f = vf[k];
                    if (fpatch[f] == t) continue;
                    if (fpatch[f] >= 0) {               /* overlap: merge     */
                        uf_union(puf, fpatch[f], t);
                    } else {
                        fpatch[f] = t;
                        if (pf_used >= NF) { overflow = 1; break; }
                        pf_all[pf_used++] = f;
                        for (int c = 0; c < 3; c++) {
                            int32_t w2 = faces[(size_t)f*3+(size_t)c];
                            if (!vmark[w2]) { vmark[w2] = 1; vlist[nvl++] = w2; }
                        }
                    }
                }
            }
            ring_lo = ring_hi;
        }
        for (int32_t vi = 0; vi < nvl; vi++) vmark[vlist[vi]] = 0;
    }
    pf_beg[n_tear] = pf_used;

    for (int32_t t = 0; t < n_tear; t++)
        if (uf_find(puf, t) == t) st.n_patches++;

    /* ---- 7. per-patch scratch (arena, allocated once) ---- */
    enum { MAXPE = PR_PATCH_MAX_FACES * 3 + 8 };
    int32_t *P_faces  = (int32_t *)ARENA_ALLOC(arena, (long)(PR_PATCH_MAX_FACES * 4) * (long)sizeof(int32_t));
    int32_t *P_verts  = (int32_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * (long)sizeof(int32_t));
    int32_t *l_of_g   = (int32_t *)ARENA_ALLOC(arena, (long)NV * (long)sizeof(int32_t));
    int32_t *pe_u     = (int32_t *)ARENA_ALLOC(arena, (long)MAXPE * (long)sizeof(int32_t));
    int32_t *pe_v     = (int32_t *)ARENA_ALLOC(arena, (long)MAXPE * (long)sizeof(int32_t));
    int32_t *pe_n     = (int32_t *)ARENA_ALLOC(arena, (long)MAXPE * (long)sizeof(int32_t));
    int32_t *pe_ue    = (int32_t *)ARENA_ALLOC(arena, (long)MAXPE * (long)sizeof(int32_t));
    int32_t *ring_eu  = (int32_t *)ARENA_ALLOC(arena, (long)MAXPE * (long)sizeof(int32_t));
    int32_t *ring_ev  = (int32_t *)ARENA_ALLOC(arena, (long)MAXPE * (long)sizeof(int32_t));
    uint8_t *on_ring  = (uint8_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS);
    int32_t *ring_deg = (int32_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * (long)sizeof(int32_t));
    int32_t *ring_nb  = (int32_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * 2L * (long)sizeof(int32_t));
    int32_t *ring_ord = (int32_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * (long)sizeof(int32_t));
    int32_t *ec_work  = (int32_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * (long)sizeof(int32_t));
    int32_t *interior = (int32_t *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * (long)sizeof(int32_t));
    double  *P2       = (double  *)ARENA_ALLOC(arena, (long)PR_PATCH_MAX_VERTS * 2L * (long)sizeof(double));
    LTri    *T        = (LTri    *)ARENA_ALLOC(arena, (long)(2 * PR_PATCH_MAX_VERTS + 8) * (long)sizeof(LTri));
    for (int32_t i = 0; i < NV; i++) l_of_g[i] = -1;

    uint8_t *fdel = (uint8_t *)ARENA_CALLOC(arena, (long)NF, 1L);
    int32_t new_cap = 4096;
    int32_t *new_tris = (int32_t *)ARENA_ALLOC(arena, (long)new_cap * 3L * (long)sizeof(int32_t));
    int32_t n_new = 0;

    /* ---- 8. repair each merged patch group ---- */
    for (int32_t g = 0; g < n_tear; g++) {
        if (uf_find(puf, g) != g) continue;

        /* collect group faces */
        int32_t npf = 0;
        int too_big = 0;
        for (int32_t t = 0; t < n_tear && !too_big; t++) {
            if (uf_find(puf, t) != g) continue;
            for (int32_t k = pf_beg[t]; k < pf_beg[t+1]; k++) {
                if (npf >= PR_PATCH_MAX_FACES * 4) { too_big = 1; break; }
                P_faces[npf++] = pf_all[k];
            }
        }
        if (too_big || npf > PR_PATCH_MAX_FACES || npf == 0) { st.n_skip_gate++; continue; }
        for (int32_t i = 1; i < npf; i++) {             /* sort (tiny)        */
            int32_t x = P_faces[i], j = i - 1;
            while (j >= 0 && P_faces[j] > x) { P_faces[j+1] = P_faces[j]; j--; }
            P_faces[j+1] = x;
        }

        /* patch verts + local map */
        int32_t npv = 0;
        int bad = 0;
        for (int32_t i = 0; i < npf && !bad; i++) {
            for (int c = 0; c < 3; c++) {
                int32_t w2 = faces[(size_t)P_faces[i]*3+(size_t)c];
                if (l_of_g[w2] >= 0) continue;
                if (npv >= PR_PATCH_MAX_VERTS) { bad = 1; break; }
                l_of_g[w2] = 0;
                P_verts[npv++] = w2;
            }
        }
        if (bad) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_gate++;
            continue;
        }
        for (int32_t i = 1; i < npv; i++) {
            int32_t x = P_verts[i], j = i - 1;
            while (j >= 0 && P_verts[j] > x) { P_verts[j+1] = P_verts[j]; j--; }
            P_verts[j+1] = x;
        }
        for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = i;

        /* classify patch edges.
         * ring edge:  patch-mult 1, mesh-mult 2 (frontier to kept mesh)
         * tear edge:  patch-mult 1, mesh-mult 1, cluster IS one of this
         *             group's tear roots (dissolves in the refill)
         * any other boundary edge in the patch (perimeter, foreign cluster)
         * or a non-manifold edge => skip the patch. */
        int32_t npe = 0;
        int skip = 0;
        for (int32_t i = 0; i < npf && !skip; i++) {
            int32_t f = P_faces[i];
            for (int c = 0; c < 3; c++) {
                int32_t u = faces[(size_t)f*3+(size_t)c];
                int32_t v = faces[(size_t)f*3+(size_t)((c+1)%3)];
                if (u > v) { int32_t tmp = u; u = v; v = tmp; }
                int32_t found = -1;
                for (int32_t k = 0; k < npe; k++)
                    if (pe_u[k] == u && pe_v[k] == v) { found = k; break; }
                if (found >= 0) { pe_n[found]++; continue; }
                if (npe >= MAXPE) { skip = 1; break; }
                pe_u[npe] = u; pe_v[npe] = v; pe_n[npe] = 1;
                pe_ue[npe] = ue_find(ue, n_ue, u, v);
                npe++;
            }
        }
        int32_t n_ring_e = 0;
        for (int32_t k = 0; k < npe && !skip; k++) {
            int32_t uei = pe_ue[k];
            assert(uei >= 0);
            if (ue[uei].mult > 2) { skip = 1; break; }
            if (pe_n[k] == 1 && ue[uei].mult == 2) {
                ring_eu[n_ring_e] = pe_u[k];
                ring_ev[n_ring_e] = pe_v[k];
                n_ring_e++;
            } else if (ue[uei].mult == 1) {
                if (!is_tear_root[uf_find(bcl, pe_u[k])]) { skip = 1; break; }
            }
        }
        if (skip || n_ring_e < 3) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_gate++;
            continue;
        }

        /* vertex-frontier gate: any patch vert used by an outside face must
         * lie on a ring edge (else deleting the patch pinches it). Also, a
         * ring edge's outside face must not already be deleted by an earlier
         * committed patch (stale orientation reference). */
        memset(on_ring, 0, (size_t)npv);
        for (int32_t k = 0; k < n_ring_e; k++) {
            on_ring[l_of_g[ring_eu[k]]] = 1;
            on_ring[l_of_g[ring_ev[k]]] = 1;
        }
        for (int32_t i = 0; i < npv && !skip; i++) {
            if (on_ring[i]) continue;
            int32_t w2 = P_verts[i];
            for (int32_t k = vf_cnt[w2]; k < vf_cnt[w2+1]; k++) {
                int32_t f = vf[k];
                if (fpatch[f] < 0 || uf_find(puf, fpatch[f]) != g) { skip = 1; break; }
            }
        }
        if (skip) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_gate++;
            continue;
        }

        /* ring must be ONE simple closed cycle */
        memset(ring_deg, 0, (size_t)npv * sizeof(int32_t));
        int ring_bad = 0;
        for (int32_t k = 0; k < n_ring_e && !ring_bad; k++) {
            int32_t a = l_of_g[ring_eu[k]], b = l_of_g[ring_ev[k]];
            if (ring_deg[a] >= 2 || ring_deg[b] >= 2) { ring_bad = 1; break; }
            ring_nb[a*2 + ring_deg[a]] = b; ring_deg[a]++;
            ring_nb[b*2 + ring_deg[b]] = a; ring_deg[b]++;
        }
        int32_t nr = 0;
        if (!ring_bad) {
            int32_t start = -1;
            for (int32_t i = 0; i < npv; i++) {
                if (ring_deg[i] == 1) { ring_bad = 1; break; }
                if (ring_deg[i] == 2 && start < 0) start = i;
            }
            if (start < 0) ring_bad = 1;
            if (!ring_bad) {
                int32_t cur = start, prev = -1;
                for (int32_t s = 0; s <= n_ring_e && cur >= 0; s++) {
                    ring_ord[nr++] = cur;
                    int32_t nx = (ring_nb[cur*2+0] != prev) ? ring_nb[cur*2+0] : ring_nb[cur*2+1];
                    prev = cur; cur = nx;
                    if (cur == start) break;
                }
                if (nr != n_ring_e) ring_bad = 1;       /* multiple cycles    */
            }
        }
        if (ring_bad) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_ring++;
            continue;
        }

        /* planarity gate: PCA of patch verts */
        double C[3] = {0,0,0};
        for (int32_t i = 0; i < npv; i++)
            for (int c = 0; c < 3; c++) C[c] += (double)verts[(size_t)P_verts[i]*3+(size_t)c];
        for (int c = 0; c < 3; c++) C[c] /= (double)npv;
        double cov[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for (int32_t i = 0; i < npv; i++) {
            double d[3];
            for (int c = 0; c < 3; c++) d[c] = (double)verts[(size_t)P_verts[i]*3+(size_t)c] - C[c];
            for (int a = 0; a < 3; a++)
                for (int b = a; b < 3; b++) cov[a][b] += d[a]*d[b];
        }
        cov[1][0]=cov[0][1]; cov[2][0]=cov[0][2]; cov[2][1]=cov[1][2];
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) cov[a][b] /= (double)npv;
        double evals[3], evecs[3][3];
        jacobi_3x3(cov, evals, evecs);
        if (sqrt(fmax(evals[0], 0.0)) > PR_PLANAR_MAX_THICK || evals[1] < 1e-12) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_gate++;
            continue;
        }
        double e1[3] = { evecs[0][2], evecs[1][2], evecs[2][2] };
        double e2[3] = { evecs[0][1], evecs[1][1], evecs[2][1] };

        /* 2D projection */
        for (int32_t i = 0; i < npv; i++) {
            double d[3];
            for (int c = 0; c < 3; c++) d[c] = (double)verts[(size_t)P_verts[i]*3+(size_t)c] - C[c];
            P2[2*i+0] = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2];
            P2[2*i+1] = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
        }

        /* ring area + CCW; simplicity; interior containment */
        double area2 = 0.0;
        for (int32_t i = 0; i < nr; i++) {
            const double *a = &P2[2*ring_ord[i]], *b = &P2[2*ring_ord[(i+1)%nr]];
            area2 += a[0]*b[1] - b[0]*a[1];
        }
        if (fabs(area2) < 10.0 * PR_EPS_2D) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_project++;
            continue;
        }
        if (area2 < 0.0) {
            for (int32_t i = 0; i < nr/2; i++) {
                int32_t t2 = ring_ord[i];
                ring_ord[i] = ring_ord[nr-1-i];
                ring_ord[nr-1-i] = t2;
            }
        }
        int proj_bad = 0;
        for (int32_t i = 0; i < nr && !proj_bad; i++) {
            int32_t a1 = ring_ord[i], b1 = ring_ord[(i+1)%nr];
            for (int32_t j = i+1; j < nr && !proj_bad; j++) {
                int32_t a2 = ring_ord[j], b2 = ring_ord[(j+1)%nr];
                if (a1==a2||a1==b2||b1==a2||b1==b2) continue;
                if (segs_cross(&P2[2*a1], &P2[2*b1], &P2[2*a2], &P2[2*b2])) proj_bad = 1;
            }
        }
        for (int32_t i = 0; i < npv && !proj_bad; i++)
            for (int32_t j = i+1; j < npv; j++) {
                double dx = P2[2*i]-P2[2*j], dy = P2[2*i+1]-P2[2*j+1];
                if (dx*dx + dy*dy < PR_EPS_2D*PR_EPS_2D) { proj_bad = 1; break; }
            }
        int32_t n_int = 0;
        for (int32_t i = 0; i < npv && !proj_bad; i++) {
            if (on_ring[i]) continue;
            int cn = 0;
            const double *p = &P2[2*i];
            for (int32_t k = 0; k < nr; k++) {
                const double *a = &P2[2*ring_ord[k]], *b = &P2[2*ring_ord[(k+1)%nr]];
                if (((a[1] <= p[1]) && (b[1] > p[1])) || ((a[1] > p[1]) && (b[1] <= p[1]))) {
                    double xt = a[0] + (p[1]-a[1]) / (b[1]-a[1]) * (b[0]-a[0]);
                    if (p[0] < xt) cn ^= 1;
                }
                if (pt_seg_d2(p, a, b) < PR_EPS_2D*PR_EPS_2D) { proj_bad = 1; break; }
            }
            if (proj_bad) break;
            if (!cn) proj_bad = 1;
            else interior[n_int++] = i;
        }
        if (proj_bad) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_project++;
            continue;
        }

        /* ---- triangulate ---- */
        int32_t n_tri = 0;
        int tri_bad = (ear_clip(P2, ring_ord, nr, T, &n_tri, ec_work) != 0);
        for (int32_t i = 0; i < n_int && !tri_bad; i++)
            tri_bad = (insert_point(P2, interior[i], T, &n_tri) != 0);
        if (!tri_bad) tri_bad = (lawson_flips(P2, T, n_tri) != 0);
        if (!tri_bad && n_tri != 2*n_int + nr - 2) tri_bad = 1;   /* Euler    */
        if (tri_bad) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_tri++;
            continue;
        }

        /* ---- orientation: every ring edge must be used exactly once by the
         * new triangles, REVERSED vs its kept outside face. Try forward and
         * globally-flipped; if neither satisfies all ring edges, skip. ---- */
        int ok_fwd = 1, ok_rev = 1;
        for (int32_t k = 0; k < n_ring_e && (ok_fwd || ok_rev); k++) {
            int32_t gu = ring_eu[k], gv = ring_ev[k];
            int32_t uei = ue_find(ue, n_ue, (gu<gv)?gu:gv, (gu<gv)?gv:gu);
            int32_t f0 = ue[uei].f0, f1 = ue[uei].f1;
            int f0_in = (fpatch[f0] >= 0 && uf_find(puf, fpatch[f0]) == g);
            int32_t fo = f0_in ? f1 : f0;
            if (fo < 0 || fdel[fo]) { ok_fwd = ok_rev = 0; break; }  /* stale */
            int32_t du = -1, dv = -1;
            for (int c = 0; c < 3; c++) {
                int32_t x = faces[(size_t)fo*3+(size_t)c];
                int32_t y = faces[(size_t)fo*3+(size_t)((c+1)%3)];
                if ((x == gu && y == gv) || (x == gv && y == gu)) { du = x; dv = y; break; }
            }
            if (du < 0) { ok_fwd = ok_rev = 0; break; }
            int32_t lu = l_of_g[du], lv = l_of_g[dv];
            int n_use = 0, dir_same = 0;
            for (int32_t t2 = 0; t2 < n_tri; t2++) {
                for (int c = 0; c < 3; c++) {
                    int32_t x = T[t2].v[c], y = T[t2].v[(c+1)%3];
                    if (x == lu && y == lv) { n_use++; dir_same = 1; }
                    else if (x == lv && y == lu) { n_use++; dir_same = 0; }
                }
            }
            if (n_use != 1) { ok_fwd = ok_rev = 0; break; }
            if (dir_same) ok_fwd = 0; else ok_rev = 0;
        }
        if (!ok_fwd && !ok_rev) {
            for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
            st.n_skip_tri++;
            continue;
        }

        /* ---- commit ---- */
        for (int32_t i = 0; i < npf; i++) fdel[P_faces[i]] = 1;
        st.faces_deleted += (size_t)npf;
        if (n_new + n_tri > new_cap) {
            int32_t cap2 = new_cap * 2 + n_tri;
            int32_t *nt2 = (int32_t *)ARENA_ALLOC(arena, (long)cap2 * 3L * (long)sizeof(int32_t));
            memcpy(nt2, new_tris, (size_t)n_new * 3 * sizeof(int32_t));
            new_tris = nt2;
            new_cap = cap2;
        }
        for (int32_t t2 = 0; t2 < n_tri; t2++) {
            int32_t a = P_verts[T[t2].v[0]];
            int32_t b = P_verts[T[t2].v[1]];
            int32_t c = P_verts[T[t2].v[2]];
            new_tris[n_new*3+0] = a;
            if (ok_fwd) { new_tris[n_new*3+1] = b; new_tris[n_new*3+2] = c; }
            else        { new_tris[n_new*3+1] = c; new_tris[n_new*3+2] = b; }
            n_new++;
        }
        st.faces_added += (size_t)n_tri;
        st.n_repaired++;

        for (int32_t i = 0; i < npv; i++) l_of_g[P_verts[i]] = -1;
    }

    /* ---- 9. assemble output ---- */
    size_t kept = 0;
    for (int32_t f = 0; f < NF; f++) if (!fdel[f]) kept++;
    size_t total = kept + (size_t)n_new;
    int32_t *out = (int32_t *)ARENA_ALLOC(arena, (long)((total ? total : 1) * 3 * sizeof(int32_t)));
    size_t w = 0;
    for (int32_t f = 0; f < NF; f++) {
        if (fdel[f]) continue;
        out[w*3+0] = faces[(size_t)f*3+0];
        out[w*3+1] = faces[(size_t)f*3+1];
        out[w*3+2] = faces[(size_t)f*3+2];
        w++;
    }
    memcpy(out + w*3, new_tris, (size_t)n_new * 3 * sizeof(int32_t));
    *out_faces = out;
    *out_nf = total;
    if (stats) *stats = st;
    return 0;
}
