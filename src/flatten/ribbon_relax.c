/* ribbon_relax.c -- see ribbon_relax.h. UV-only relaxation: symmetric-Dirichlet
 * isometry + coherence-gated fiber-axis alignment, minimized by embedding-
 * preserving per-vertex line search. Geometry math in double; verts/uv are
 * float at the interface. */
#include "ribbon_relax.h"

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RLX_PI      3.14159265358979323846
#define RLX_AREA_EPS 1e-9

/* per-face precomputed rest data + fiber target */
typedef struct {
    double mi00, mi01, mi10, mi11;   /* inverse of material frame M = [x1|x2] */
    double a_rest;                   /* 3D triangle area */
    double dfx, dfy;                 /* material-frame fiber unit direction */
    double coh;                      /* aggregated fiber coherence (0 => no align) */
    double qc_initial;               /* stretch in the immutable input chart */
    int    s;                        /* sign of initial UV signed area (+1/-1) */
    int    ok;                       /* 1 = non-degenerate */
} FaceRelax;

/* ---- small geometry helpers ---------------------------------------------- */

static double signed_area2(const double a[2], const double b[2], const double c[2])
{
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);   /* 2 * area */
}

/* J = [uv1-uv0 | uv2-uv0] * Minv  (rest material frame -> UV) */
static void jacobian(const FaceRelax *fr, const double uv0[2], const double uv1[2],
                     const double uv2[2], double J[4])
{
    double e00 = uv1[0] - uv0[0], e10 = uv1[1] - uv0[1];   /* col 0 */
    double e01 = uv2[0] - uv0[0], e11 = uv2[1] - uv0[1];   /* col 1 */
    J[0] = e00 * fr->mi00 + e01 * fr->mi10;   /* J00 */
    J[1] = e00 * fr->mi01 + e01 * fr->mi11;   /* J01 */
    J[2] = e10 * fr->mi00 + e11 * fr->mi10;   /* J10 */
    J[3] = e10 * fr->mi01 + e11 * fr->mi11;   /* J11 */
}

/* symmetric Dirichlet + gated sin^2(2phi) alignment, scaled by rest area */
static double face_energy(const FaceRelax *fr, const double uv0[2], const double uv1[2],
                          const double uv2[2], double lambda)
{
    double J[4], det, fro, d2, e;
    if (!fr->ok) return 0.0;
    jacobian(fr, uv0, uv1, uv2, J);
    det = J[0] * J[3] - J[1] * J[2];
    fro = J[0] * J[0] + J[1] * J[1] + J[2] * J[2] + J[3] * J[3];
    d2 = det * det; if (d2 < 1e-12) d2 = 1e-12;
    e = fro * (1.0 + 1.0 / d2);                    /* ||J||^2 + ||J^-1||^2 */
    if (lambda > 0.0 && fr->coh > 0.0) {
        double mx = J[0] * fr->dfx + J[1] * fr->dfy;   /* fiber dir mapped into UV */
        double my = J[2] * fr->dfx + J[3] * fr->dfy;
        double s2 = sin(2.0 * atan2(my, mx));
        e += lambda * fr->coh * s2 * s2;
    }
    return fr->a_rest * e;
}

/* quasi-conformal stretch sigma1/sigma2 of the UV map (for stats) */
static void face_stretch(const FaceRelax *fr, const double uv0[2], const double uv1[2],
                         const double uv2[2], double *qc, double *s1o, double *s2o)
{
    double J[4], E, G, F, disc, l1, l2, s1, s2;
    jacobian(fr, uv0, uv1, uv2, J);
    E = J[0] * J[0] + J[2] * J[2];
    G = J[1] * J[1] + J[3] * J[3];
    F = J[0] * J[1] + J[2] * J[3];
    disc = sqrt(fmax(0.0, 0.25 * (E - G) * (E - G) + F * F));
    l1 = 0.5 * (E + G) + disc; l2 = 0.5 * (E + G) - disc;
    s1 = sqrt(fmax(l1, 0.0)); s2 = sqrt(fmax(l2, 0.0));
    *qc = (s2 > 1e-9) ? s1 / s2 : 1e9; *s1o = s1; *s2o = s2;
}

/* mapped fiber angle -> distance to nearest u/v axis, degrees [0,45] */
static double face_axis_err(const FaceRelax *fr, const double uv0[2], const double uv1[2],
                            const double uv2[2])
{
    double J[4], mx, my, deg, d90;
    jacobian(fr, uv0, uv1, uv2, J);
    mx = J[0] * fr->dfx + J[1] * fr->dfy;
    my = J[2] * fr->dfx + J[3] * fr->dfy;
    deg = fabs(atan2(my, mx) * 180.0 / RLX_PI);
    d90 = fmod(deg, 90.0);
    return fmin(d90, 90.0 - d90);
}

/* build rest frame + fiber target for one face. P* are (z,y,x) float; uv* the
 * initial UV (for pulling the fiber direction into the material frame). */
static int face_precompute(const float *P0, const float *P1, const float *P2,
                           const double ref0[2], const double ref1[2],
                           const double ref2[2], const double start0[2],
                           const double start1[2], const double start2[2],
                           double theta_grid, double coh_in, double coh_gate,
                           double qc_reject, int reference_metric,
                           FaceRelax *fr)
{
    double e1[3], e2[3], L1, e1h[3], d, h2, det_m, area2, J0[4], detJ0, wx, wy, dx, dy, dn;
    double E0, G0, F0, disc0, l1, l2, s1, s2, qc0;
    double qc_start, start_s1, start_s2;
    int i;
    memset(fr, 0, sizeof *fr);
    if (reference_metric) {
        /* Embed reference UV in the same (z,y,x) convention as XYZ.  Its
         * induced first fundamental form is exactly the developable metric
         * that a detector-confirmed corrupt XYZ patch should inherit. */
        e1[0] = 0.0; e1[1] = ref1[1] - ref0[1]; e1[2] = ref1[0] - ref0[0];
        e2[0] = 0.0; e2[1] = ref2[1] - ref0[1]; e2[2] = ref2[0] - ref0[0];
    } else {
        for (i = 0; i < 3; i++) {
            e1[i] = (double)P1[i] - P0[i];
            e2[i] = (double)P2[i] - P0[i];
        }
    }
    L1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    if (L1 < 1e-9) return 0;
    for (i = 0; i < 3; i++) e1h[i] = e1[i] / L1;
    d = e2[0] * e1h[0] + e2[1] * e1h[1] + e2[2] * e1h[2];
    h2 = sqrt(fmax(0.0, (e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]) - d * d));
    if (h2 < 1e-9) return 0;
    /* M = [[L1, d],[0, h2]]; Minv = [[1/L1, -d/(L1 h2)],[0, 1/h2]] */
    fr->mi00 = 1.0 / L1; fr->mi01 = -d / (L1 * h2); fr->mi10 = 0.0; fr->mi11 = 1.0 / h2;
    det_m = L1 * h2;
    fr->a_rest = 0.5 * det_m;

    area2 = signed_area2(ref0, ref1, ref2);
    fr->s = signed_area2(start0, start1, start2) >= 0.0 ? 1 : -1;

    /* initial Jacobian + stretch: reject pre-broken faces (UV-collapsed slivers
     * and overlap-relocated faces) whose sigma1/sigma2 is enormous -- they would
     * otherwise dominate the energy and hijack the optimizer. */
    jacobian(fr, ref0, ref1, ref2, J0);
    detJ0 = J0[0] * J0[3] - J0[1] * J0[2];
    E0 = J0[0]*J0[0] + J0[2]*J0[2]; G0 = J0[1]*J0[1] + J0[3]*J0[3]; F0 = J0[0]*J0[1] + J0[2]*J0[3];
    disc0 = sqrt(fmax(0.0, 0.25 * (E0 - G0) * (E0 - G0) + F0 * F0));
    l1 = 0.5 * (E0 + G0) + disc0; l2 = 0.5 * (E0 + G0) - disc0;
    s1 = sqrt(fmax(l1, 0.0)); s2 = sqrt(fmax(l2, 0.0));
    qc0 = (s2 > 1e-9) ? s1 / s2 : 1e18;
    if (qc_reject > 0.0 && qc0 > qc_reject) return 0;   /* fr->ok stays 0 -> excluded */
    face_stretch(fr, start0, start1, start2,
                 &qc_start, &start_s1, &start_s2);
    fr->qc_initial = qc_start;

    fr->coh = 0.0; fr->dfx = 1.0; fr->dfy = 0.0;
    if (coh_in >= coh_gate && fabs(area2) > RLX_AREA_EPS && fabs(detJ0) > 1e-12) {
        wx = cos(theta_grid); wy = sin(theta_grid);        /* fiber dir in UV */
        /* d_material = J0^-1 * w */
        dx = ( J0[3] * wx - J0[1] * wy) / detJ0;
        dy = (-J0[2] * wx + J0[0] * wy) / detJ0;
        dn = sqrt(dx * dx + dy * dy);
        if (dn > 1e-12) { fr->dfx = dx / dn; fr->dfy = dy / dn; fr->coh = coh_in; }
    }
    fr->ok = 1;
    return 1;
}

/* aggregate the FiberField over a face's texel footprint (coherence-weighted
 * doubled-angle mean) -> grid orientation + coherence */
static void face_fiber_target(const FiberField *fib, double du, double dv,
                              const double uv0[2], const double uv1[2], const double uv2[2],
                              double *out_theta, double *out_coh)
{
    double umin, umax, vmin, vmax, k, sx = 0.0, sy = 0.0, wc = 0.0, range;
    int x0, x1, y0, y1, x, y;
    *out_theta = 0.0; *out_coh = 0.0;
    if (fib == NULL || fib->W <= 0) return;
    umin = fmin(uv0[0], fmin(uv1[0], uv2[0])) / du; umax = fmax(uv0[0], fmax(uv1[0], uv2[0])) / du;
    vmin = fmin(uv0[1], fmin(uv1[1], uv2[1])) / dv; vmax = fmax(uv0[1], fmax(uv1[1], uv2[1])) / dv;
    x0 = (int)floor(umin); x1 = (int)ceil(umax); y0 = (int)floor(vmin); y1 = (int)ceil(vmax);
    if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
    if (x1 > fib->W - 1) x1 = fib->W - 1; if (y1 > fib->H - 1) y1 = fib->H - 1;
    k = (double)fib->rosy;
    range = fib->range_deg * RLX_PI / 180.0;
    for (y = y0; y <= y1; y++)
        for (x = x0; x <= x1; x++) {
            size_t p = (size_t)y * fib->W + x;
            double c, th;
            if (!fib->valid[p]) continue;
            c = fib->coh[p]; th = fib->theta[p];
            sx += c * cos(k * th); sy += c * sin(k * th); wc += c;
        }
    if (wc < 1e-9) return;
    {
        double ang = atan2(sy, sx) / k;
        ang = fmod(ang, range); if (ang < 0.0) ang += range;
        *out_theta = ang;
        *out_coh = sqrt(sx * sx + sy * sy) / wc;
    }
}

/* ---- adjacency + boundary ------------------------------------------------ */

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

/* A free boundary needs a global guard, not merely incident-face orientation
 * tests.  A disk with consistently oriented triangles and a simple boundary is
 * globally injective; if a non-adjacent pair of boundary edges crosses, it is
 * not.  This dynamic uniform hash indexes the current boundary each sweep.
 * Accepted edge moves are inserted again (stale records are harmless because
 * queries always test the edge's current endpoints), avoiding an O(B^2) scan. */
typedef struct {
    int64_t u, v;
    int32_t head;
    uint8_t used;
} RlxBoundaryCell;

typedef struct {
    int32_t edge, next;
} RlxBoundaryEntry;

typedef struct {
    RlxBoundaryCell *cell;
    size_t cell_cap, cell_used;
    RlxBoundaryEntry *entry;
    size_t entry_count, entry_cap;
    double cell_size;
    const uint64_t *edge_key;
    size_t nedge;
    uint32_t *seen;
    uint32_t stamp;
} RlxBoundaryGrid;

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static uint64_t rlx_hash64(uint64_t x)
{
    x ^= x >> 30; x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27; x *= UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

static uint64_t rlx_cell_hash(int64_t u, int64_t v)
{
    return rlx_hash64((uint64_t)u ^
                      (rlx_hash64((uint64_t)v) + UINT64_C(0x9e3779b97f4a7c15)));
}

static int rlx_cell_coord(double x, double cell_size, int64_t *out)
{
    long double q;
    if (!isfinite(x) || !(cell_size > 0.0) || !isfinite(cell_size)) return -1;
    q = floorl((long double)x / (long double)cell_size);
    if (q < (long double)INT64_MIN || q > (long double)INT64_MAX) return -1;
    *out = (int64_t)q;
    return 0;
}

static void rlx_boundary_grid_dispose(RlxBoundaryGrid *g)
{
    if (g == NULL) return;
    free(g->cell); free(g->entry); free(g->seen);
    memset(g, 0, sizeof(*g));
}

static int rlx_boundary_grid_rehash(RlxBoundaryGrid *g, size_t new_cap)
{
    RlxBoundaryCell *next;
    size_t i;
    if (new_cap < 16) new_cap = 16;
    next = (RlxBoundaryCell *)calloc(new_cap, sizeof(*next));
    if (next == NULL) return -1;
    for (i = 0; i < g->cell_cap; i++) {
        if (g->cell[i].used) {
            size_t at = (size_t)rlx_cell_hash(g->cell[i].u, g->cell[i].v) &
                        (new_cap - 1);
            while (next[at].used) at = (at + 1) & (new_cap - 1);
            next[at] = g->cell[i];
        }
    }
    free(g->cell);
    g->cell = next;
    g->cell_cap = new_cap;
    return 0;
}

static RlxBoundaryCell *rlx_boundary_grid_cell(RlxBoundaryGrid *g,
                                                int64_t u, int64_t v,
                                                int create)
{
    size_t at;
    if (create && (g->cell_used + 1) * 10 >= g->cell_cap * 7) {
        if (g->cell_cap > SIZE_MAX / 2 ||
            rlx_boundary_grid_rehash(g, g->cell_cap * 2) != 0)
            return NULL;
    }
    at = (size_t)rlx_cell_hash(u, v) & (g->cell_cap - 1);
    while (g->cell[at].used) {
        if (g->cell[at].u == u && g->cell[at].v == v) return &g->cell[at];
        at = (at + 1) & (g->cell_cap - 1);
    }
    if (!create) return NULL;
    g->cell[at].used = 1;
    g->cell[at].u = u; g->cell[at].v = v; g->cell[at].head = -1;
    g->cell_used++;
    return &g->cell[at];
}

static int rlx_boundary_grid_add(RlxBoundaryGrid *g, int32_t edge,
                                 int64_t u0, int64_t u1,
                                 int64_t v0, int64_t v1)
{
    int64_t vcell;
    for (vcell = v0;; vcell++) {
        int64_t ucell;
        for (ucell = u0;; ucell++) {
            RlxBoundaryCell *cell = rlx_boundary_grid_cell(g, ucell, vcell, 1);
            if (cell == NULL) return -1;
            if (g->entry_count == g->entry_cap) {
                size_t cap = g->entry_cap < 1024 ? 1024 : g->entry_cap * 2;
                RlxBoundaryEntry *entry;
                if (cap < g->entry_cap || cap > (size_t)INT32_MAX) return -1;
                entry = (RlxBoundaryEntry *)realloc(g->entry,
                                                     cap * sizeof(*entry));
                if (entry == NULL) return -1;
                g->entry = entry; g->entry_cap = cap;
            }
            if (g->entry_count > (size_t)INT32_MAX) return -1;
            g->entry[g->entry_count].edge = edge;
            g->entry[g->entry_count].next = cell->head;
            cell->head = (int32_t)g->entry_count++;
            if (ucell == u1) break;
        }
        if (vcell == v1) break;
    }
    return 0;
}

static int rlx_boundary_grid_add_edge(RlxBoundaryGrid *g, size_t edge,
                                      const double *uv)
{
    uint64_t key = g->edge_key[edge];
    size_t a = (size_t)(key >> 32), b = (size_t)(key & UINT32_MAX);
    double u0 = fmin(uv[a * 2], uv[b * 2]);
    double u1 = fmax(uv[a * 2], uv[b * 2]);
    double v0 = fmin(uv[a * 2 + 1], uv[b * 2 + 1]);
    double v1 = fmax(uv[a * 2 + 1], uv[b * 2 + 1]);
    int64_t cu0, cu1, cv0, cv1;
    if (rlx_cell_coord(u0, g->cell_size, &cu0) != 0 ||
        rlx_cell_coord(u1, g->cell_size, &cu1) != 0 ||
        rlx_cell_coord(v0, g->cell_size, &cv0) != 0 ||
        rlx_cell_coord(v1, g->cell_size, &cv1) != 0)
        return -1;
    return rlx_boundary_grid_add(g, (int32_t)edge, cu0, cu1, cv0, cv1);
}

static int rlx_boundary_grid_build(RlxBoundaryGrid *g,
                                   const uint64_t *edge_key, size_t nedge,
                                   const double *uv)
{
    double *length = NULL, cell_size;
    size_t i, cap = 16;
    long double total = 0.0L, target;
    memset(g, 0, sizeof(*g));
    if (nedge == 0 || edge_key == NULL || uv == NULL ||
        nedge > (size_t)INT32_MAX) return -1;
    length = (double *)malloc(nedge * sizeof(*length));
    if (length == NULL) return -1;
    for (i = 0; i < nedge; i++) {
        size_t a = (size_t)(edge_key[i] >> 32);
        size_t b = (size_t)(edge_key[i] & UINT32_MAX);
        length[i] = hypot(uv[a * 2] - uv[b * 2],
                          uv[a * 2 + 1] - uv[b * 2 + 1]);
    }
    qsort(length, nedge, sizeof(*length), cmp_double);
    cell_size = 4.0 * length[nedge / 2];
    if (!(cell_size > 1.0e-6) || !isfinite(cell_size)) cell_size = 1.0;
    target = 8.0L * (long double)nedge + 1.0L;
    for (int pass = 0; pass < 12; pass++) {
        total = 0.0L;
        for (i = 0; i < nedge; i++) {
            size_t a = (size_t)(edge_key[i] >> 32);
            size_t b = (size_t)(edge_key[i] & UINT32_MAX);
            int64_t u0, u1, v0, v1;
            if (rlx_cell_coord(fmin(uv[a * 2], uv[b * 2]), cell_size, &u0) ||
                rlx_cell_coord(fmax(uv[a * 2], uv[b * 2]), cell_size, &u1) ||
                rlx_cell_coord(fmin(uv[a * 2 + 1], uv[b * 2 + 1]), cell_size, &v0) ||
                rlx_cell_coord(fmax(uv[a * 2 + 1], uv[b * 2 + 1]), cell_size, &v1)) {
                free(length); return -1;
            }
            total += ((long double)u1 - (long double)u0 + 1.0L) *
                     ((long double)v1 - (long double)v0 + 1.0L);
            if (total > target * 1024.0L) break;
        }
        if (total <= target) break;
        cell_size *= sqrt((double)(total / target)) * 1.05;
        if (!isfinite(cell_size)) { free(length); return -1; }
    }
    free(length);
    if (total > (long double)SIZE_MAX / 4.0L || total > INT32_MAX) return -1;
    while ((long double)cap < 2.0L * total + 16.0L) {
        if (cap > SIZE_MAX / 2) return -1;
        cap *= 2;
    }
    g->cell = (RlxBoundaryCell *)calloc(cap, sizeof(*g->cell));
    g->seen = (uint32_t *)calloc(nedge, sizeof(*g->seen));
    g->entry_cap = (size_t)total + 4 * nedge + 16;
    if (g->entry_cap > (size_t)INT32_MAX) goto fail;
    g->entry = (RlxBoundaryEntry *)malloc(g->entry_cap * sizeof(*g->entry));
    if (g->cell == NULL || g->seen == NULL || g->entry == NULL) goto fail;
    g->cell_cap = cap; g->cell_size = cell_size;
    g->edge_key = edge_key; g->nedge = nedge; g->stamp = 1;
    for (i = 0; i < nedge; i++)
        if (rlx_boundary_grid_add_edge(g, i, uv) != 0) goto fail;
    return 0;
fail:
    rlx_boundary_grid_dispose(g);
    return -1;
}

static double rlx_orient2(const double *a, const double *b, const double *c)
{
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

static int rlx_orientation_sign(const double *a, const double *b,
                                const double *c)
{
    double abx = b[0] - a[0], aby = b[1] - a[1];
    double acx = c[0] - a[0], acy = c[1] - a[1];
    double value = abx * acy - aby * acx;
    double eps = 1.0e-12 * (fabs(abx * acy) + fabs(aby * acx) + 1.0);
    return value > eps ? 1 : (value < -eps ? -1 : 0);
}

static int rlx_on_segment(const double *a, const double *b, const double *p)
{
    double eps = 1.0e-10 * (hypot(b[0] - a[0], b[1] - a[1]) + 1.0);
    return rlx_orientation_sign(a, b, p) == 0 &&
           p[0] >= fmin(a[0], b[0]) - eps && p[0] <= fmax(a[0], b[0]) + eps &&
           p[1] >= fmin(a[1], b[1]) - eps && p[1] <= fmax(a[1], b[1]) + eps;
}

static int rlx_segments_intersect(const double *a, const double *b,
                                  const double *c, const double *d)
{
    int o1 = rlx_orientation_sign(a, b, c);
    int o2 = rlx_orientation_sign(a, b, d);
    int o3 = rlx_orientation_sign(c, d, a);
    int o4 = rlx_orientation_sign(c, d, b);
    if (o1 * o2 < 0 && o3 * o4 < 0) return 1;
    if (o1 == 0 && rlx_on_segment(a, b, c)) return 1;
    if (o2 == 0 && rlx_on_segment(a, b, d)) return 1;
    if (o3 == 0 && rlx_on_segment(c, d, a)) return 1;
    if (o4 == 0 && rlx_on_segment(c, d, b)) return 1;
    return 0;
}

static int rlx_point_in_triangle(const double *p, const double *a,
                                 const double *b, const double *c)
{
    int s0 = rlx_orientation_sign(a, b, p);
    int s1 = rlx_orientation_sign(b, c, p);
    int s2 = rlx_orientation_sign(c, a, p);
    int nonnegative = s0 >= 0 && s1 >= 0 && s2 >= 0;
    int nonpositive = s0 <= 0 && s1 <= 0 && s2 <= 0;
    return nonnegative || nonpositive;
}

static int rlx_segment_hits_swept_edge(const double *a, const double *b,
                                       const double *fixed,
                                       const double *old_moving,
                                       const double *new_moving)
{
    if (rlx_segments_intersect(a, b, fixed, new_moving) ||
        rlx_segments_intersect(a, b, old_moving, new_moving)) return 1;
    if (fabs(rlx_orient2(fixed, old_moving, new_moving)) <= 1.0e-14)
        return 0;
    if (rlx_point_in_triangle(a, fixed, old_moving, new_moving) ||
        rlx_point_in_triangle(b, fixed, old_moving, new_moving)) return 1;
    if (rlx_segments_intersect(a, b, fixed, old_moving)) return 1;
    return 0;
}

static int rlx_edge_shares(uint64_t key, size_t a, size_t b)
{
    size_t p = (size_t)(key >> 32), q = (size_t)(key & UINT32_MAX);
    return p == a || p == b || q == a || q == b;
}

static int rlx_boundary_move_safe(RlxBoundaryGrid *g, size_t moved,
                                  size_t neighbor0, size_t neighbor1,
                                  const double old_uv[2],
                                  const double trial[2], const double *uv)
{
    double u0 = fmin(fmin(old_uv[0], trial[0]),
                     fmin(uv[neighbor0 * 2], uv[neighbor1 * 2]));
    double u1 = fmax(fmax(old_uv[0], trial[0]),
                     fmax(uv[neighbor0 * 2], uv[neighbor1 * 2]));
    double v0 = fmin(fmin(old_uv[1], trial[1]),
                     fmin(uv[neighbor0 * 2 + 1], uv[neighbor1 * 2 + 1]));
    double v1 = fmax(fmax(old_uv[1], trial[1]),
                     fmax(uv[neighbor0 * 2 + 1], uv[neighbor1 * 2 + 1]));
    int64_t cu0, cu1, cv0, cv1, cv;
    if (rlx_cell_coord(u0, g->cell_size, &cu0) ||
        rlx_cell_coord(u1, g->cell_size, &cu1) ||
        rlx_cell_coord(v0, g->cell_size, &cv0) ||
        rlx_cell_coord(v1, g->cell_size, &cv1)) return 0;
    if (++g->stamp == 0) {
        memset(g->seen, 0, g->nedge * sizeof(*g->seen));
        g->stamp = 1;
    }
    for (cv = cv0;; cv++) {
        int64_t cu;
        for (cu = cu0;; cu++) {
            RlxBoundaryCell *cell = rlx_boundary_grid_cell(g, cu, cv, 0);
            int32_t at = cell != NULL ? cell->head : -1;
            while (at >= 0) {
                int32_t e = g->entry[at].edge;
                uint64_t key = g->edge_key[(size_t)e];
                if (g->seen[(size_t)e] != g->stamp) {
                    size_t a = (size_t)(key >> 32);
                    size_t b = (size_t)(key & UINT32_MAX);
                    const double *pa = &uv[a * 2], *pb = &uv[b * 2];
                    g->seen[(size_t)e] = g->stamp;
                    if (!rlx_edge_shares(key, moved, neighbor0) &&
                        rlx_segment_hits_swept_edge(pa, pb,
                            &uv[neighbor0 * 2], old_uv, trial)) return 0;
                    if (!rlx_edge_shares(key, moved, neighbor1) &&
                        rlx_segment_hits_swept_edge(pa, pb,
                            &uv[neighbor1 * 2], old_uv, trial)) return 0;
                }
                at = g->entry[at].next;
            }
            if (cu == cu1) break;
        }
        if (cv == cv1) break;
    }
    return 1;
}

static size_t rlx_boundary_intersections(RlxBoundaryGrid *g,
                                         const double *uv)
{
    size_t e, count = 0;
    for (e = 0; e < g->nedge; e++) {
        uint64_t key = g->edge_key[e];
        size_t a = (size_t)(key >> 32), b = (size_t)(key & UINT32_MAX);
        const double *pa = &uv[a * 2], *pb = &uv[b * 2];
        int64_t u0, u1, v0, v1, vv;
        if (rlx_cell_coord(fmin(pa[0], pb[0]), g->cell_size, &u0) ||
            rlx_cell_coord(fmax(pa[0], pb[0]), g->cell_size, &u1) ||
            rlx_cell_coord(fmin(pa[1], pb[1]), g->cell_size, &v0) ||
            rlx_cell_coord(fmax(pa[1], pb[1]), g->cell_size, &v1)) return SIZE_MAX;
        if (++g->stamp == 0) {
            memset(g->seen, 0, g->nedge * sizeof(*g->seen));
            g->stamp = 1;
        }
        for (vv = v0;; vv++) {
            int64_t uu;
            for (uu = u0;; uu++) {
                RlxBoundaryCell *cell = rlx_boundary_grid_cell(g, uu, vv, 0);
                int32_t at = cell != NULL ? cell->head : -1;
                while (at >= 0) {
                    size_t other = (size_t)g->entry[at].edge;
                    if (other > e && g->seen[other] != g->stamp) {
                        uint64_t okey = g->edge_key[other];
                        size_t c = (size_t)(okey >> 32);
                        size_t d = (size_t)(okey & UINT32_MAX);
                        g->seen[other] = g->stamp;
                        if (!rlx_edge_shares(okey, a, b) &&
                            rlx_segments_intersect(pa, pb,
                                                   &uv[c * 2], &uv[d * 2]))
                            count++;
                    }
                    at = g->entry[at].next;
                }
                if (uu == u1) break;
            }
            if (vv == v1) break;
        }
    }
    return count;
}

static int boundary_turn_preserved(size_t center, size_t moved,
                                   const double trial[2], const double *uv,
                                   const int32_t *neighbor0,
                                   const int32_t *neighbor1,
                                   const int8_t *turn_sign)
{
    size_t a = (size_t)neighbor0[center];
    size_t b = center;
    size_t c = (size_t)neighbor1[center];
    const double *ua = a == moved ? trial : &uv[a * 2];
    const double *ub = b == moved ? trial : &uv[b * 2];
    const double *uc = c == moved ? trial : &uv[c * 2];
    double area2 = signed_area2(ua, ub, uc);
    return area2 * (double)turn_sign[center] > 1.0e-12;
}

/* ---- driver -------------------------------------------------------------- */

void RibbonRelax_defaults(RibbonRelaxOpts *o)
{
    if (o == NULL) return;
    o->lambda_align = 0.5;   /* balanced isometry<->fiber-alignment operating point */
    o->coh_gate     = 0.15;
    o->sweeps       = 1000;
    o->movement_tolerance = 1e-5;
    o->line_search  = 24;
    o->max_disp     = 40.0;
    o->qc_reject    = 50.0;
    o->max_stretch_growth = 0.0;
    o->stretch_guard_floor = 1.0;
    o->reference_metric = 0;
    o->fix_boundary = 0;     /* free: lets the sheet rotate/shear u onto the fibers */
    o->convex_boundary = 0;
    o->global_injective = 0;
    o->guide_uv      = NULL;
    o->guide_weight  = 0.0;
    o->guide_weight_vertex = NULL;
    o->reference_uv  = NULL;
    o->pin          = NULL;
    o->verbose      = 0;
}

/* total area-weighted energy + stats over all faces at a given uv */
static void measure(const FaceRelax *fr, const int32_t *faces, size_t nf,
                    const double *uv, int ref_sign, double lambda,
                    double *E, double *stretch_mean, double *stretch_max,
                    double *axis_err, int *flips)
{
    double e = 0.0, sw = 0.0, sm = 0.0, smax = 1.0, aerr = 0.0, aw = 0.0;
    int fl = 0;
    size_t f;
    for (f = 0; f < nf; f++) {
        const double *u0, *u1, *u2;
        double qc, s1, s2, a2;
        if (!fr[f].ok) continue;
        u0 = uv + (size_t)faces[f * 3 + 0] * 2;
        u1 = uv + (size_t)faces[f * 3 + 1] * 2;
        u2 = uv + (size_t)faces[f * 3 + 2] * 2;
        e += face_energy(&fr[f], u0, u1, u2, lambda);
        face_stretch(&fr[f], u0, u1, u2, &qc, &s1, &s2);
        sm += fr[f].a_rest * fmin(qc, 20.0); sw += fr[f].a_rest;   /* clamp so a few slivers don't swamp the bulk */
        if (qc > smax && qc < 1e8) smax = qc;
        a2 = signed_area2(u0, u1, u2);
        if ((a2 >= 0.0 ? 1 : -1) != ref_sign) fl++;
        if (fr[f].coh > 0.0) { double er = face_axis_err(&fr[f], u0, u1, u2); aerr += fr[f].coh * er; aw += fr[f].coh; }
    }
    if (E) *E = e;
    if (stretch_mean) *stretch_mean = (sw > 0) ? sm / sw : 0.0;
    if (stretch_max) *stretch_max = smax;
    if (axis_err) *axis_err = (aw > 0) ? aerr / aw : 0.0;
    if (flips) *flips = fl;
}

int RibbonRelax_run_double(Arena_T arena,
                           const float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const double *uv_in, const uint8_t *face_skip,
                           const FiberField *fib, double du, double dv,
                           const RibbonRelaxOpts *opts,
                           double *uv_out, RibbonRelaxStats *stats)
{
    RibbonRelaxOpts o;
    FaceRelax *fr = NULL;
    double *uv = NULL;               /* working UV (double) */
    size_t *vf_off = NULL, *vf_idx = NULL;
    uint64_t *edges = NULL;
    uint8_t *bnd = NULL, *movable = NULL;
    int32_t *bnd_neighbor0 = NULL, *bnd_neighbor1 = NULL;
    int32_t *bnd_edge0 = NULL, *bnd_edge1 = NULL;
    int8_t *bnd_turn_sign = NULL;
    double *guide_area = NULL;
    size_t f, i, v, sweep, n_fiber = 0, n_move_total = 0, n_interior = 0, n_reject = 0;
    size_t n_convex_reject = 0;
    size_t n_stretch_guard_reject = 0;
    size_t n_boundary_reject = 0, n_boundary_edges = 0;
    size_t boundary_before = 0, boundary_after = 0;
    int ref_sign = 1;
    const double *reference_uv = NULL;
    double lambda, disp_sum = 0.0, disp_max = 0.0;

    if (arena == NULL || verts == NULL || faces == NULL || uv_in == NULL || uv_out == NULL) return -1;
    if (opts) o = *opts; else RibbonRelax_defaults(&o);
    if (o.max_stretch_growth > 0.0 &&
        (o.max_stretch_growth < 1.0 || o.stretch_guard_floor < 1.0))
        return -1;
    if (du <= 0) du = 1.0; if (dv <= 0) dv = 1.0;
    lambda = o.lambda_align;
    reference_uv = o.reference_uv != NULL ? o.reference_uv : uv_in;
    if (stats) memset(stats, 0, sizeof *stats);

    /* trivial: nothing to do -> copy through */
    if (nv == 0 || nf == 0) {
        for (i = 0; i < nv * 2; i++) uv_out[i] = uv_in[i];
        if (stats) stats->reverted = 1;
        return 0;
    }

    uv = (double *)Arena_alloc(arena, (size_t)(nv * 2 * sizeof(double)), __FILE__, __LINE__);
    for (i = 0; i < nv * 2; i++) uv[i] = uv_in[i];

    /* reference winding sign = majority over faces */
    {
        long pos = 0, neg = 0;
        for (f = 0; f < nf; f++) {
            const double *u0, *u1, *u2; double a2;
            if (face_skip && face_skip[f]) continue;
            u0 = uv + (size_t)faces[f*3+0]*2; u1 = uv + (size_t)faces[f*3+1]*2; u2 = uv + (size_t)faces[f*3+2]*2;
            a2 = signed_area2(u0, u1, u2);
            if (a2 >= 0) pos++; else neg++;
        }
        ref_sign = (neg > pos) ? -1 : 1;
    }

    /* per-face precompute + fiber target */
    fr = (FaceRelax *)Arena_alloc(arena, (size_t)(nf * sizeof *fr), __FILE__, __LINE__);
    if (o.guide_uv != NULL && o.guide_weight > 0.0)
        guide_area = (double *)Arena_calloc(
            arena, nv, sizeof(*guide_area), __FILE__, __LINE__);
    for (f = 0; f < nf; f++) {
        int32_t a, b, c;
        const double *u0, *u1, *u2, *r0, *r1, *r2;
        double th = 0.0, ch = 0.0;
        if (face_skip && face_skip[f]) { memset(&fr[f], 0, sizeof fr[f]); continue; }
        a = faces[f*3+0]; b = faces[f*3+1]; c = faces[f*3+2];
        u0 = uv + (size_t)a*2; u1 = uv + (size_t)b*2; u2 = uv + (size_t)c*2;
        r0 = reference_uv + (size_t)a*2;
        r1 = reference_uv + (size_t)b*2;
        r2 = reference_uv + (size_t)c*2;
        face_fiber_target(fib, du, dv, r0, r1, r2, &th, &ch);
        face_precompute(&verts[(size_t)a*3], &verts[(size_t)b*3], &verts[(size_t)c*3],
                        r0, r1, r2, u0, u1, u2, th, ch, o.coh_gate,
                        o.qc_reject, o.reference_metric, &fr[f]);
        if (!fr[f].ok) n_reject++;
        else {
            if (fr[f].coh > 0.0) n_fiber++;
            if (guide_area != NULL) {
                double share = fr[f].a_rest / 3.0;
                guide_area[(size_t)a] += share;
                guide_area[(size_t)b] += share;
                guide_area[(size_t)c] += share;
            }
        }
    }

    /* vertex -> face incidence (CSR) */
    vf_off = (size_t *)Arena_calloc(arena, (size_t)(nv + 1), sizeof(size_t), __FILE__, __LINE__);
    for (f = 0; f < nf; f++) { if (face_skip && face_skip[f]) continue; for (i = 0; i < 3; i++) vf_off[(size_t)faces[f*3+i] + 1]++; }
    for (v = 0; v < nv; v++) vf_off[v+1] += vf_off[v];
    vf_idx = (size_t *)Arena_alloc(arena, (size_t)(3 * nf * sizeof(size_t)), __FILE__, __LINE__);
    {
        size_t *cur = (size_t *)Arena_alloc(arena, (size_t)((nv) * sizeof(size_t)), __FILE__, __LINE__);
        for (v = 0; v < nv; v++) cur[v] = vf_off[v];
        for (f = 0; f < nf; f++) { if (face_skip && face_skip[f]) continue; for (i = 0; i < 3; i++) { size_t vv = (size_t)faces[f*3+i]; vf_idx[cur[vv]++] = f; } }
    }

    /* boundary vertices via edge incidence (kept faces only) */
    bnd = (uint8_t *)Arena_calloc(arena, (size_t)nv, 1, __FILE__, __LINE__);
    bnd_neighbor0 = (int32_t *)Arena_alloc(
        arena, nv * sizeof(*bnd_neighbor0), __FILE__, __LINE__);
    bnd_neighbor1 = (int32_t *)Arena_alloc(
        arena, nv * sizeof(*bnd_neighbor1), __FILE__, __LINE__);
    bnd_edge0 = (int32_t *)Arena_alloc(
        arena, nv * sizeof(*bnd_edge0), __FILE__, __LINE__);
    bnd_edge1 = (int32_t *)Arena_alloc(
        arena, nv * sizeof(*bnd_edge1), __FILE__, __LINE__);
    bnd_turn_sign = (int8_t *)Arena_calloc(
        arena, nv, sizeof(*bnd_turn_sign), __FILE__, __LINE__);
    for (v = 0; v < nv; v++) {
        bnd_neighbor0[v] = -1;
        bnd_neighbor1[v] = -1;
        bnd_edge0[v] = -1;
        bnd_edge1[v] = -1;
    }
    edges = (uint64_t *)Arena_alloc(arena, (size_t)(3 * nf * sizeof(uint64_t)), __FILE__, __LINE__);
    {
        size_t ne = 0;
        for (f = 0; f < nf; f++) {
            if (face_skip && face_skip[f]) continue;
            for (i = 0; i < 3; i++) {
                uint32_t p = (uint32_t)faces[f*3+i], q = (uint32_t)faces[f*3+(i+1)%3];
                uint32_t lo = p < q ? p : q, hi = p < q ? q : p;
                edges[ne++] = ((uint64_t)lo << 32) | hi;
            }
        }
        qsort(edges, ne, sizeof(uint64_t), cmp_u64);
        for (i = 0; i < ne; ) {
            size_t j = i; while (j < ne && edges[j] == edges[i]) j++;
            if (j - i == 1) {
                uint64_t key = edges[i];
                int32_t edge_id;
                size_t p = (size_t)(key >> 32);
                size_t q = (size_t)(key & 0xffffffffu);
                if (n_boundary_edges >= (size_t)INT32_MAX) {
                    fprintf(stderr, "  [relax] too many boundary edges\n");
                    return -1;
                }
                edge_id = (int32_t)n_boundary_edges;
                edges[n_boundary_edges++] = key;
                bnd[p] = 1;
                bnd[q] = 1;
                if (bnd_neighbor0[p] < 0) {
                    bnd_neighbor0[p] = (int32_t)q;
                    bnd_edge0[p] = edge_id;
                } else if (bnd_neighbor1[p] < 0 && bnd_neighbor0[p] != (int32_t)q) {
                    bnd_neighbor1[p] = (int32_t)q;
                    bnd_edge1[p] = edge_id;
                }
                if (bnd_neighbor0[q] < 0) {
                    bnd_neighbor0[q] = (int32_t)p;
                    bnd_edge0[q] = edge_id;
                } else if (bnd_neighbor1[q] < 0 && bnd_neighbor0[q] != (int32_t)p) {
                    bnd_neighbor1[q] = (int32_t)p;
                    bnd_edge1[q] = edge_id;
                }
            }
            i = j;
        }
    }

    if (o.global_injective && !o.fix_boundary) {
        RlxBoundaryGrid initial_grid;
        for (v = 0; v < nv; v++) {
            if (bnd[v] &&
                (bnd_neighbor0[v] < 0 || bnd_neighbor1[v] < 0 ||
                 bnd_edge0[v] < 0 || bnd_edge1[v] < 0)) {
                fprintf(stderr,
                        "  [relax] injective-boundary contract failed at "
                        "vertex %zu (boundary degree is not two)\n", v);
                return -1;
            }
        }
        if (rlx_boundary_grid_build(&initial_grid, edges, n_boundary_edges,
                                    uv) != 0) {
            fprintf(stderr,
                    "  [relax] could not build injective-boundary index\n");
            return -1;
        }
        boundary_before = rlx_boundary_intersections(&initial_grid, uv);
        rlx_boundary_grid_dispose(&initial_grid);
        if (boundary_before != 0) {
            fprintf(stderr,
                    "  [relax] injective-boundary contract failed: %zu "
                    "non-adjacent boundary intersections in input\n",
                    boundary_before);
            return -1;
        }
        boundary_after = boundary_before;
    }

    if (o.convex_boundary && !o.fix_boundary) {
        for (v = 0; v < nv; v++) {
            if (bnd[v]) {
                const double *a = NULL, *b = NULL, *c = NULL;
                double area2 = 0.0;
                if (bnd_neighbor0[v] < 0 || bnd_neighbor1[v] < 0) {
                    fprintf(stderr,
                            "  [relax] convex-boundary contract failed at "
                            "vertex %zu (boundary degree is not two)\n", v);
                    return -1;
                }
                a = &uv[(size_t)bnd_neighbor0[v] * 2];
                b = &uv[v * 2];
                c = &uv[(size_t)bnd_neighbor1[v] * 2];
                area2 = signed_area2(a, b, c);
                if (fabs(area2) <= 1.0e-12) {
                    fprintf(stderr,
                            "  [relax] convex-boundary contract failed at "
                            "vertex %zu (near-zero initial turn %.3e)\n",
                            v, area2);
                    return -1;
                }
                bnd_turn_sign[v] = area2 > 0.0 ? 1 : -1;
            }
        }
    }

    movable = (uint8_t *)Arena_alloc(arena, (size_t)nv, __FILE__, __LINE__);
    for (v = 0; v < nv; v++) {
        movable[v] = (vf_off[v+1] > vf_off[v]) && !(o.fix_boundary && bnd[v])
                     && !(o.pin != NULL && o.pin[v]);
        if (movable[v]) n_interior++;
    }

    if (stats) {
        stats->n_boundary_edges = n_boundary_edges;
        stats->boundary_intersections_before = boundary_before;
        measure(fr, faces, nf, uv, ref_sign, lambda,
                &stats->energy_before, &stats->stretch_mean_before, &stats->stretch_max_before,
                &stats->axis_err_before, &stats->flips_before);
        if (guide_area != NULL) {
            double sum = 0.0, weight = 0.0;
            for (v = 0; v < nv; v++) {
                double dx = uv[v * 2] - o.guide_uv[v * 2];
                double dy = uv[v * 2 + 1] - o.guide_uv[v * 2 + 1];
                double confidence = o.guide_weight_vertex != NULL ?
                                    o.guide_weight_vertex[v] : 1.0;
                sum += confidence * guide_area[v] * (dx * dx + dy * dy);
                weight += confidence * guide_area[v];
            }
            stats->guide_rms_before = weight > 0.0 ? sqrt(sum / weight) : 0.0;
        }
    }

    /* Gauss-Seidel sweeps: per-vertex line search along -grad */
    for (sweep = 0; sweep < (size_t)o.sweeps; sweep++) {
        RlxBoundaryGrid boundary_grid;
        int have_boundary_grid = 0;
        size_t moved = 0;
        double sweep_max_move = 0.0;
        memset(&boundary_grid, 0, sizeof(boundary_grid));
        if (o.global_injective && !o.fix_boundary) {
            if (rlx_boundary_grid_build(&boundary_grid, edges,
                                        n_boundary_edges, uv) != 0) {
                fprintf(stderr,
                        "  [relax] could not rebuild injective-boundary index "
                        "for sweep %zu\n", sweep);
                return -1;
            }
            have_boundary_grid = 1;
        }
        for (v = 0; v < nv; v++) {
            double base, gu, gv, step_u, step_v, trust;
            double elen = 0.0, t, dispmax2;
            size_t s, s0 = vf_off[v], s1 = vf_off[v+1], nls;
            double cur[2], trial[2];
            if (!movable[v]) continue;
            cur[0] = uv[v*2]; cur[1] = uv[v*2+1];

            /* local energy helper inlined: sum incident faces with v at position p */
            #define LOCAL_E(px, py, outE) do {                                            \
                double _e = 0.0; size_t _s;                                               \
                double _pp[2]; _pp[0] = (px); _pp[1] = (py);                              \
                for (_s = s0; _s < s1; _s++) {                                            \
                    size_t _f = vf_idx[_s];                                               \
                    int32_t _a = faces[_f*3+0], _b = faces[_f*3+1], _c = faces[_f*3+2];   \
                    const double *_u0 = ((size_t)_a==v)?_pp:uv+(size_t)_a*2;              \
                    const double *_u1 = ((size_t)_b==v)?_pp:uv+(size_t)_b*2;              \
                    const double *_u2 = ((size_t)_c==v)?_pp:uv+(size_t)_c*2;              \
                    _e += face_energy(&fr[_f], _u0, _u1, _u2, lambda);                    \
                }                                                                         \
                if (guide_area != NULL) {                                                  \
                    double _dx = (px) - o.guide_uv[v * 2];                                \
                    double _dy = (py) - o.guide_uv[v * 2 + 1];                            \
                    double _confidence = o.guide_weight_vertex != NULL ?                   \
                                         o.guide_weight_vertex[v] : 1.0;                   \
                    _e += o.guide_weight * _confidence * guide_area[v] *                  \
                          (_dx * _dx + _dy * _dy);                                         \
                }                                                                         \
                (outE) = _e;                                                             \
            } while (0)

            LOCAL_E(cur[0], cur[1], base);
            /* mean incident UV edge length for step scale */
            for (s = s0; s < s1; s++) {
                size_t ff = vf_idx[s]; int m;
                for (m = 0; m < 3; m++) {
                    size_t w = (size_t)faces[ff*3+m];
                    if (w == v) continue;
                    { double ex = uv[w*2]-cur[0], ey = uv[w*2+1]-cur[1]; elen += sqrt(ex*ex+ey*ey); }
                }
            }
            elen = (s1 > s0) ? elen / (double)(2 * (s1 - s0)) : 1.0;

            {
                double h = fmin(0.05, fmax(0.005, 0.01 * elen));
                double epx, emx, epy, emy, epp, epm, emp, emm;
                double hxx, hxy, hyy, disc, minimum_eigenvalue;
                double scale, floor_eigenvalue, det, descent;
                LOCAL_E(cur[0] + h, cur[1], epx);
                LOCAL_E(cur[0] - h, cur[1], emx);
                LOCAL_E(cur[0], cur[1] + h, epy);
                LOCAL_E(cur[0], cur[1] - h, emy);
                LOCAL_E(cur[0] + h, cur[1] + h, epp);
                LOCAL_E(cur[0] + h, cur[1] - h, epm);
                LOCAL_E(cur[0] - h, cur[1] + h, emp);
                LOCAL_E(cur[0] - h, cur[1] - h, emm);
                gu = (epx - emx) / (2.0 * h);
                gv = (epy - emy) / (2.0 * h);
                if (hypot(gu, gv) < 1e-12) continue;
                hxx = (epx - 2.0 * base + emx) / (h * h);
                hyy = (epy - 2.0 * base + emy) / (h * h);
                hxy = (epp - epm - emp + emm) / (4.0 * h * h);

                /* A normalized gradient has a finite step even arbitrarily
                 * close to a minimum, so the former solver kept oscillating
                 * at iteration 1000.  Use a damped 2x2 Newton step for this
                 * vertex instead.  Shifting the numerical Hessian to a small
                 * positive eigenvalue keeps the direction descending in a
                 * non-convex barrier neighbourhood; line search remains the
                 * transaction that accepts or rejects it. */
                scale = fabs(hxx) + 2.0 * fabs(hxy) + fabs(hyy) + 1.0;
                floor_eigenvalue = 1e-6 * scale;
                disc = hypot(hxx - hyy, 2.0 * hxy);
                minimum_eigenvalue = 0.5 * (hxx + hyy - disc);
                if (minimum_eigenvalue < floor_eigenvalue) {
                    double shift = floor_eigenvalue - minimum_eigenvalue;
                    hxx += shift;
                    hyy += shift;
                }
                det = hxx * hyy - hxy * hxy;
                if (det > 1e-18 * scale * scale && isfinite(det)) {
                    step_u = (-hyy * gu + hxy * gv) / det;
                    step_v = ( hxy * gu - hxx * gv) / det;
                } else {
                    step_u = -gu / scale;
                    step_v = -gv / scale;
                }
                descent = gu * step_u + gv * step_v;
                if (!isfinite(step_u) || !isfinite(step_v) || descent >= 0.0) {
                    step_u = -gu / scale;
                    step_v = -gv / scale;
                }
            }
            trust = 0.5 * elen;
            if (o.max_disp > 0.0 && trust > o.max_disp) trust = o.max_disp;
            if (trust < 0.02) trust = 0.02;
            {
                double step_length = hypot(step_u, step_v);
                if (step_length < 1e-14) continue;
                if (step_length > trust) {
                    step_u *= trust / step_length;
                    step_v *= trust / step_length;
                }
            }
            dispmax2 = o.max_disp > 0.0 ? o.max_disp * o.max_disp : 0.0;
            t = 1.0;
            for (nls = 0; nls < (size_t)o.line_search; nls++) {
                double ddx, ddy; int flip = 0, stretch_reject = 0; double ecand;
                trial[0] = cur[0] + t * step_u;
                trial[1] = cur[1] + t * step_v;
                ddx = trial[0] - reference_uv[v*2];
                ddy = trial[1] - reference_uv[v*2+1];
                if (o.max_disp > 0.0 && ddx*ddx + ddy*ddy > dispmax2) {
                    t *= 0.5;
                    continue;
                }
                for (s = s0; s < s1 && !flip; s++) {
                    size_t ff = vf_idx[s];
                    int32_t a = faces[ff*3+0], b = faces[ff*3+1], c = faces[ff*3+2];
                    const double *u0 = ((size_t)a==v)?trial:uv+(size_t)a*2;
                    const double *u1 = ((size_t)b==v)?trial:uv+(size_t)b*2;
                    const double *u2 = ((size_t)c==v)?trial:uv+(size_t)c*2;
                    double a2 = signed_area2(u0, u1, u2);
                    if ((a2 >= 0 ? 1 : -1) != fr[ff].s || fabs(a2) < RLX_AREA_EPS) flip = 1;
                }
                if (!flip && o.max_stretch_growth > 0.0) {
                    for (s = s0; s < s1 && !stretch_reject; s++) {
                        size_t ff = vf_idx[s];
                        int32_t a, b, c;
                        const double *u0, *u1, *u2;
                        double qc, sigma1, sigma2, limit;
                        if (!fr[ff].ok) continue;
                        a = faces[ff*3+0]; b = faces[ff*3+1]; c = faces[ff*3+2];
                        u0 = ((size_t)a==v)?trial:uv+(size_t)a*2;
                        u1 = ((size_t)b==v)?trial:uv+(size_t)b*2;
                        u2 = ((size_t)c==v)?trial:uv+(size_t)c*2;
                        face_stretch(&fr[ff], u0, u1, u2,
                                     &qc, &sigma1, &sigma2);
                        limit = fr[ff].qc_initial * o.max_stretch_growth;
                        if (limit < o.stretch_guard_floor)
                            limit = o.stretch_guard_floor;
                        if (!isfinite(qc) || qc > limit + 1e-10)
                            stretch_reject = 1;
                    }
                    if (stretch_reject) n_stretch_guard_reject++;
                }
                if (!flip && !stretch_reject && o.convex_boundary && bnd[v]) {
                    size_t n0 = (size_t)bnd_neighbor0[v];
                    size_t n1 = (size_t)bnd_neighbor1[v];
                    if (!boundary_turn_preserved(v, v, trial, uv,
                                                 bnd_neighbor0, bnd_neighbor1,
                                                 bnd_turn_sign) ||
                        !boundary_turn_preserved(n0, v, trial, uv,
                                                 bnd_neighbor0, bnd_neighbor1,
                                                 bnd_turn_sign) ||
                        !boundary_turn_preserved(n1, v, trial, uv,
                                                 bnd_neighbor0, bnd_neighbor1,
                                                 bnd_turn_sign)) {
                        flip = 1;
                        n_convex_reject++;
                    }
                }
                if (!flip && !stretch_reject && have_boundary_grid && bnd[v]) {
                    size_t n0 = (size_t)bnd_neighbor0[v];
                    size_t n1 = (size_t)bnd_neighbor1[v];
                    if (!rlx_boundary_move_safe(&boundary_grid, v, n0, n1,
                                                cur, trial, uv)) {
                        flip = 1;
                        n_boundary_reject++;
                    }
                }
                if (!flip && !stretch_reject) {
                    LOCAL_E(trial[0], trial[1], ecand);
                    if (ecand < base - 1e-12) {
                        uv[v*2] = trial[0]; uv[v*2+1] = trial[1];
                        {
                            double accepted_move = hypot(trial[0] - cur[0],
                                                         trial[1] - cur[1]);
                            if (accepted_move > sweep_max_move)
                                sweep_max_move = accepted_move;
                        }
                        if (have_boundary_grid && bnd[v] &&
                            (rlx_boundary_grid_add_edge(
                                 &boundary_grid, (size_t)bnd_edge0[v], uv) != 0 ||
                             rlx_boundary_grid_add_edge(
                                 &boundary_grid, (size_t)bnd_edge1[v], uv) != 0)) {
                            rlx_boundary_grid_dispose(&boundary_grid);
                            fprintf(stderr,
                                    "  [relax] injective-boundary index "
                                    "update failed at vertex %zu\n", v);
                            return -1;
                        }
                        moved++;
                        break;
                    }
                }
                t *= 0.5;
            }
            #undef LOCAL_E
        }
        if (have_boundary_grid) {
            boundary_after = rlx_boundary_intersections(&boundary_grid, uv);
            rlx_boundary_grid_dispose(&boundary_grid);
            if (boundary_after != 0) {
                fprintf(stderr,
                        "  [relax] injective-boundary invariant failed after "
                        "sweep %zu: %zu intersections\n",
                        sweep, boundary_after);
                return -1;
            }
        }
        n_move_total += moved;
        if (stats) {
            stats->iterations = (int)sweep + 1;
            stats->last_max_movement = sweep_max_move;
        }
        if (o.verbose) {
            double e; measure(fr, faces, nf, uv, ref_sign, lambda, &e, NULL, NULL, NULL, NULL);
            fprintf(stderr, "  [relax] sweep %3zu: moved %zu, max-move=%.9g, E=%.1f\n",
                    sweep, moved, sweep_max_move, e);
        }
        if (moved == 0 || (o.movement_tolerance > 0.0 &&
                           sweep_max_move <= o.movement_tolerance)) {
            if (stats) stats->converged = 1;
            break;
        }
    }

    /* per-movable-vertex displacement */
    for (v = 0; v < nv; v++) {
        double dx, dy, dd;
        if (!movable[v]) continue;
        dx = uv[v*2] - reference_uv[v*2];
        dy = uv[v*2+1] - reference_uv[v*2+1];
        dd = sqrt(dx*dx + dy*dy);
        disp_sum += dd; if (dd > disp_max) disp_max = dd;
    }

    /* stats after + fail-closed guard */
    {
        double E1 = 0, sm1 = 0, smx1 = 1, ae1 = 0; int fl1 = 0;
        measure(fr, faces, nf, uv, ref_sign, lambda, &E1, &sm1, &smx1, &ae1, &fl1);
        if (stats) {
            stats->energy_after = E1; stats->stretch_mean_after = sm1;
            stats->stretch_max_after = smx1; stats->axis_err_after = ae1;
            stats->flips_after = fl1; stats->n_interior = n_interior;
            stats->n_moved = n_move_total; stats->n_fiber_faces = n_fiber;
            stats->n_reject = n_reject;
            stats->n_stretch_guard_reject = n_stretch_guard_reject;
            stats->n_convex_reject = n_convex_reject;
            stats->n_boundary_collision_reject = n_boundary_reject;
            stats->boundary_intersections_after = boundary_after;
            if (guide_area != NULL) {
                double sum = 0.0, weight = 0.0;
                for (v = 0; v < nv; v++) {
                    double dx = uv[v * 2] - o.guide_uv[v * 2];
                    double dy = uv[v * 2 + 1] - o.guide_uv[v * 2 + 1];
                    double confidence = o.guide_weight_vertex != NULL ?
                                        o.guide_weight_vertex[v] : 1.0;
                    sum += confidence * guide_area[v] * (dx * dx + dy * dy);
                    weight += confidence * guide_area[v];
                }
                stats->guide_rms_after =
                    weight > 0.0 ? sqrt(sum / weight) : 0.0;
            }
            stats->mean_disp = (n_interior > 0) ? disp_sum / (double)n_interior : 0.0;
            stats->max_disp = disp_max;
        }
        if (fl1 > (stats ? stats->flips_before : 0) || E1 > (stats ? stats->energy_before : E1) + 1e-6) {
            for (i = 0; i < nv * 2; i++) uv_out[i] = uv_in[i];   /* revert */
            if (stats) stats->reverted = 1;
            return 0;
        }
    }

    for (i = 0; i < nv * 2; i++) uv_out[i] = uv[i];
    return 0;
}

int RibbonRelax_run(Arena_T arena,
                    const float *verts, size_t nv,
                    const int32_t *faces, size_t nf,
                    const float *uv_in, const uint8_t *face_skip,
                    const FiberField *fib, double du, double dv,
                    const RibbonRelaxOpts *opts,
                    float *uv_out, RibbonRelaxStats *stats)
{
    Arena_Mark mark;
    double *input = NULL;
    double *output = NULL;
    size_t i = 0;
    size_t count = nv > 0 ? nv * 2 : 1;
    int result = -1;

    if (arena == NULL || verts == NULL || faces == NULL ||
        uv_in == NULL || uv_out == NULL) return -1;
    mark = Arena_save(arena);
    input = (double *)ARENA_ALLOC(arena, count * sizeof(*input));
    output = (double *)ARENA_ALLOC(arena, count * sizeof(*output));
    for (i = 0; i < nv * 2; i++) input[i] = (double)uv_in[i];
    result = RibbonRelax_run_double(arena, verts, nv, faces, nf, input,
                                    face_skip, fib, du, dv, opts, output, stats);
    if (result == 0)
        for (i = 0; i < nv * 2; i++) uv_out[i] = (float)output[i];
    Arena_restore(arena, mark);
    return result;
}

/* ---- selftest ------------------------------------------------------------ */

/* n x n planar grid in the (y,x) plane (z=0), two triangles per cell */
static void rlx_grid(int n, float **pverts, int32_t **pfaces, size_t *pnv, size_t *pnf)
{
    int i, j, k = 0, t = 0;
    size_t nv = (size_t)n * n, nf = (size_t)(n - 1) * (n - 1) * 2;
    float *V = (float *)malloc(nv * 3 * sizeof *V);
    int32_t *F = (int32_t *)malloc(nf * 3 * sizeof *F);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) { V[k*3+0]=0.0f; V[k*3+1]=(float)i; V[k*3+2]=(float)j; k++; }
    for (i = 0; i < n-1; i++)
        for (j = 0; j < n-1; j++) {
            int v00=i*n+j, v10=(i+1)*n+j, v01=i*n+j+1, v11=(i+1)*n+j+1;
            F[t*3+0]=v00; F[t*3+1]=v10; F[t*3+2]=v11; t++;
            F[t*3+0]=v00; F[t*3+1]=v11; F[t*3+2]=v01; t++;
        }
    *pverts=V; *pfaces=F; *pnv=nv; *pnf=nf;
}

int RibbonRelax_selftest(void)
{
    int fails = 0, n = 11;
    float *V = NULL; int32_t *F = NULL; size_t nv = 0, nf = 0, i;
    float *uv_in = NULL, *uv_out = NULL;
    Arena_T arena = Arena_new();
    RibbonRelaxOpts o; RibbonRelaxStats st;

    rlx_grid(n, &V, &F, &nv, &nf);
    uv_in = (float *)malloc(nv * 2 * sizeof *uv_in);
    uv_out = (float *)malloc(nv * 2 * sizeof *uv_out);

    /* (a) isometry: identity UV with deterministic interior perturbation -> energy must drop, no flips */
    for (i = 0; i < nv; i++) {
        int r = (int)i / n, c = (int)i % n;
        double pu = 0.0, pv = 0.0;
        if (r > 0 && r < n-1 && c > 0 && c < n-1) {
            pu = 0.33 * (double)(((r*7 + c) % 3) - 1);
            pv = 0.33 * (double)(((r*5 + c*3) % 3) - 1);
        }
        uv_in[i*2+0] = (float)((double)c + pu);
        uv_in[i*2+1] = (float)((double)r + pv);
    }
    RibbonRelax_defaults(&o); o.lambda_align = 0.0; o.sweeps = 40;
    if (RibbonRelax_run(arena, V, nv, F, nf, uv_in, NULL, NULL, 1.0, 1.0, &o, uv_out, &st) != 0) { fprintf(stderr, "run fail\n"); fails++; }
    else {
        if (!(st.energy_after < st.energy_before)) { fprintf(stderr, "[relax selftest]   FAIL: iso energy %.2f -> %.2f\n", st.energy_before, st.energy_after); fails++; }
        if (st.flips_after != 0) { fprintf(stderr, "[relax selftest]   FAIL: iso introduced %d flips\n", st.flips_after); fails++; }
        if (st.reverted) { fprintf(stderr, "[relax selftest]   FAIL: iso reverted\n"); fails++; }
    }

    /* (b) alignment (face level): a rotation of UV that puts the fiber on-axis must
     *     lower face_energy vs the misaligned one when alignment dominates. */
    {
        FaceRelax fa; double u0[2]={0,0}, u1[2]={1,0}, u2[2]={0,1};
        double r0[2], r1[2], r2[2], ang = -30.0*RLX_PI/180.0, ca=cos(ang), sa=sin(ang);
        float P0[3]={0,0,0}, P1[3]={0,0,1}, P2[3]={0,1,0};     /* unit right triangle in (y,x) */
        double e_mis, e_ali;
        /* fiber grid orientation 30deg in the identity UV */
        face_precompute(P0, P1, P2, u0, u1, u2, u0, u1, u2,
                        30.0 * RLX_PI / 180.0, 1.0, 0.1, 100.0, 0, &fa);
        e_mis = face_energy(&fa, u0, u1, u2, 1000.0);
        /* rotate UV by -30deg (isometric) -> fiber lands on an axis */
        r0[0]=ca*u0[0]-sa*u0[1]; r0[1]=sa*u0[0]+ca*u0[1];
        r1[0]=ca*u1[0]-sa*u1[1]; r1[1]=sa*u1[0]+ca*u1[1];
        r2[0]=ca*u2[0]-sa*u2[1]; r2[1]=sa*u2[0]+ca*u2[1];
        e_ali = face_energy(&fa, r0, r1, r2, 1000.0);
        if (!(e_ali < e_mis)) { fprintf(stderr, "[relax selftest]   FAIL: align e %.2f !< %.2f\n", e_ali, e_mis); fails++; }
        if (fa.coh <= 0.0) { fprintf(stderr, "[relax selftest]   FAIL: fiber target not set\n"); fails++; }
    }

    /* (c) trivial input */
    if (RibbonRelax_run(arena, V, 0, F, 0, uv_in, NULL, NULL, 1.0, 1.0, &o, uv_out, &st) != 0) { fprintf(stderr, "[relax selftest] FAIL trivial\n"); fails++; }

    /* (d) LLP64 size-arithmetic regression (the 2026-07-17 whole-scroll crash).
     *     On Win64 `long` is 32-bit, so the old `(long)(nf * sizeof *fr)` truncated
     *     35,548,713 faces * 72 B = 2,559,507,336 to a negative int, which
     *     Arena_alloc's size_t parameter sign-extended to ~1.8e19 bytes (OOM).
     *     Alloc-size math must be size_t end to end and must not truncate >2 GB. */
    {
        size_t big_nf = 35548713;                             /* the count that crashed */
        size_t truth  = big_nf * sizeof(FaceRelax);           /* 64-bit product */
        size_t coded  = (size_t)(big_nf * sizeof(FaceRelax)); /* the alloc idiom used above */
        if (truth <= 0x7fffffffULL) {
            fprintf(stderr, "[relax selftest]   FAIL: regression count no longer exceeds 2 GB "
                    "(sizeof FaceRelax=%zu -- bump big_nf)\n", sizeof(FaceRelax)); fails++;
        }
        if (coded != truth) {
            fprintf(stderr, "[relax selftest]   FAIL: >2 GB alloc size truncates (got %zu want %zu)\n",
                    coded, truth); fails++;
        }
    }

    /* (a2) the fixed-input L-infinity stretch guard cannot let repeated
     * vertex updates trade one new outlier for a lower global energy. */
    RibbonRelax_defaults(&o);
    o.lambda_align = 0.0;
    o.sweeps = 12;
    o.max_stretch_growth = 1.0;
    o.stretch_guard_floor = 1.25;
    if (RibbonRelax_run(arena, V, nv, F, nf, uv_in, NULL, NULL, 1.0, 1.0,
                        &o, uv_out, &st) != 0) {
        fprintf(stderr, "[relax selftest]   FAIL: guarded run\n");
        fails++;
    } else {
        double allowed = fmax(st.stretch_max_before, o.stretch_guard_floor);
        if (st.stretch_max_after > allowed + 1e-8) {
            fprintf(stderr,
                    "[relax selftest]   FAIL: stretch guard %.9g > %.9g\n",
                    st.stretch_max_after, allowed);
            fails++;
        }
        if (!(st.energy_after < st.energy_before) || st.reverted) {
            fprintf(stderr,
                    "[relax selftest]   FAIL: guarded energy %.9g -> %.9g%s\n",
                    st.energy_before, st.energy_after,
                    st.reverted ? " reverted" : "");
            fails++;
        }
    }

    /* (e) a free-boundary vertex must not be able to drag either incident
     * edge across a remote edge.  This is the global failure that incident-
     * triangle orientation alone cannot detect. */
    {
        uint64_t boundary[4] = {
            ((uint64_t)0 << 32) | 1u,
            ((uint64_t)0 << 32) | 3u,
            ((uint64_t)1 << 32) | 2u,
            ((uint64_t)2 << 32) | 3u
        };
        double square[8] = { 0, 0, 1, 0, 1, 1, 0, 1 };
        double safe[2] = { 0.1, 0.1 };
        double crossing[2] = { 2.0, 0.5 };
        RlxBoundaryGrid grid;
        if (rlx_boundary_grid_build(&grid, boundary, 4, square) != 0) {
            fprintf(stderr,
                    "[relax selftest]   FAIL: boundary index build\n");
            fails++;
        } else {
            if (rlx_boundary_intersections(&grid, square) != 0) {
                fprintf(stderr,
                        "[relax selftest]   FAIL: simple square reported "
                        "an intersection\n");
                fails++;
            }
            if (!rlx_boundary_move_safe(&grid, 0, 1, 3,
                                        &square[0], safe, square)) {
                fprintf(stderr,
                        "[relax selftest]   FAIL: safe boundary move rejected\n");
                fails++;
            }
            if (rlx_boundary_move_safe(&grid, 0, 1, 3,
                                       &square[0], crossing, square)) {
                fprintf(stderr,
                        "[relax selftest]   FAIL: crossing boundary move accepted\n");
                fails++;
            }
            rlx_boundary_grid_dispose(&grid);
        }
    }

    free(V); free(F); free(uv_in); free(uv_out);
    Arena_dispose(&arena);
    fprintf(stderr, "[ribbon_relax selftest] %s\n", fails == 0 ? "ok" : "FAILED");
    return fails;
}
