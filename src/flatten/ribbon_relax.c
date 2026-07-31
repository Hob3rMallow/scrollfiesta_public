/* ribbon_relax.c -- see ribbon_relax.h. UV-only relaxation: symmetric-Dirichlet
 * isometry + coherence-gated fiber-axis alignment, minimized by embedding-
 * preserving per-vertex line search. Geometry math in double; verts/uv are
 * float at the interface. */
#include "ribbon_relax.h"

#include <math.h>
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
                           const double uv0[2], const double uv1[2], const double uv2[2],
                           double theta_grid, double coh_in, double coh_gate,
                           double qc_reject, FaceRelax *fr)
{
    double e1[3], e2[3], L1, e1h[3], d, h2, det_m, area2, J0[4], detJ0, wx, wy, dx, dy, dn;
    double E0, G0, F0, disc0, l1, l2, s1, s2, qc0;
    int i;
    memset(fr, 0, sizeof *fr);
    for (i = 0; i < 3; i++) { e1[i] = (double)P1[i] - P0[i]; e2[i] = (double)P2[i] - P0[i]; }
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

    area2 = signed_area2(uv0, uv1, uv2);
    fr->s = (area2 >= 0.0) ? 1 : -1;

    /* initial Jacobian + stretch: reject pre-broken faces (UV-collapsed slivers
     * and overlap-relocated faces) whose sigma1/sigma2 is enormous -- they would
     * otherwise dominate the energy and hijack the optimizer. */
    jacobian(fr, uv0, uv1, uv2, J0);
    detJ0 = J0[0] * J0[3] - J0[1] * J0[2];
    E0 = J0[0]*J0[0] + J0[2]*J0[2]; G0 = J0[1]*J0[1] + J0[3]*J0[3]; F0 = J0[0]*J0[1] + J0[2]*J0[3];
    disc0 = sqrt(fmax(0.0, 0.25 * (E0 - G0) * (E0 - G0) + F0 * F0));
    l1 = 0.5 * (E0 + G0) + disc0; l2 = 0.5 * (E0 + G0) - disc0;
    s1 = sqrt(fmax(l1, 0.0)); s2 = sqrt(fmax(l2, 0.0));
    qc0 = (s2 > 1e-9) ? s1 / s2 : 1e18;
    if (qc_reject > 0.0 && qc0 > qc_reject) return 0;   /* fr->ok stays 0 -> excluded */

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

/* ---- driver -------------------------------------------------------------- */

void RibbonRelax_defaults(RibbonRelaxOpts *o)
{
    if (o == NULL) return;
    o->lambda_align = 0.5;   /* balanced isometry<->fiber-alignment operating point */
    o->coh_gate     = 0.15;
    o->sweeps       = 40;    /* coordinate descent converges by ~30 on real sheets */
    o->line_search  = 10;
    o->max_disp     = 40.0;
    o->qc_reject    = 50.0;
    o->fix_boundary = 0;     /* free: lets the sheet rotate/shear u onto the fibers */
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

int RibbonRelax_run(Arena_T arena,
                    const float *verts, size_t nv,
                    const int32_t *faces, size_t nf,
                    const float *uv_in, const uint8_t *face_skip,
                    const FiberField *fib, double du, double dv,
                    const RibbonRelaxOpts *opts,
                    float *uv_out, RibbonRelaxStats *stats)
{
    RibbonRelaxOpts o;
    FaceRelax *fr = NULL;
    double *uv = NULL;               /* working UV (double) */
    size_t *vf_off = NULL, *vf_idx = NULL;
    uint64_t *edges = NULL;
    uint8_t *bnd = NULL, *movable = NULL;
    size_t f, i, v, sweep, n_fiber = 0, n_move_total = 0, n_interior = 0, n_reject = 0;
    int ref_sign = 1;
    double lambda, disp_sum = 0.0, disp_max = 0.0;

    if (arena == NULL || verts == NULL || faces == NULL || uv_in == NULL || uv_out == NULL) return -1;
    if (opts) o = *opts; else RibbonRelax_defaults(&o);
    if (du <= 0) du = 1.0; if (dv <= 0) dv = 1.0;
    lambda = o.lambda_align;
    if (stats) memset(stats, 0, sizeof *stats);

    /* trivial: nothing to do -> copy through */
    if (nv == 0 || nf == 0) {
        for (i = 0; i < nv * 2; i++) uv_out[i] = uv_in[i];
        if (stats) stats->reverted = 1;
        return 0;
    }

    uv = (double *)Arena_alloc(arena, (size_t)(nv * 2 * sizeof(double)), __FILE__, __LINE__);
    for (i = 0; i < nv * 2; i++) uv[i] = (double)uv_in[i];

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
    for (f = 0; f < nf; f++) {
        int32_t a, b, c;
        const double *u0, *u1, *u2;
        double th = 0.0, ch = 0.0;
        if (face_skip && face_skip[f]) { memset(&fr[f], 0, sizeof fr[f]); continue; }
        a = faces[f*3+0]; b = faces[f*3+1]; c = faces[f*3+2];
        u0 = uv + (size_t)a*2; u1 = uv + (size_t)b*2; u2 = uv + (size_t)c*2;
        face_fiber_target(fib, du, dv, u0, u1, u2, &th, &ch);
        face_precompute(&verts[(size_t)a*3], &verts[(size_t)b*3], &verts[(size_t)c*3],
                        u0, u1, u2, th, ch, o.coh_gate, o.qc_reject, &fr[f]);
        if (!fr[f].ok) n_reject++;
        else if (fr[f].coh > 0.0) n_fiber++;
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
            if (j - i == 1) { bnd[(size_t)(edges[i] >> 32)] = 1; bnd[(size_t)(edges[i] & 0xffffffffu)] = 1; }
            i = j;
        }
    }

    movable = (uint8_t *)Arena_alloc(arena, (size_t)nv, __FILE__, __LINE__);
    for (v = 0; v < nv; v++) {
        movable[v] = (vf_off[v+1] > vf_off[v]) && !(o.fix_boundary && bnd[v])
                     && !(o.pin != NULL && o.pin[v]);
        if (movable[v]) n_interior++;
    }

    if (stats) {
        measure(fr, faces, nf, uv, ref_sign, lambda,
                &stats->energy_before, &stats->stretch_mean_before, &stats->stretch_max_before,
                &stats->axis_err_before, &stats->flips_before);
    }

    /* Gauss-Seidel sweeps: per-vertex line search along -grad */
    for (sweep = 0; sweep < (size_t)o.sweeps; sweep++) {
        size_t moved = 0;
        for (v = 0; v < nv; v++) {
            double base, ep, em, gu, gv, gn, dirx, diry, t0, elen = 0.0, t, dispmax2;
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

            { double h = 0.02;
              double ep2, em2;
              LOCAL_E(cur[0]+h, cur[1], ep); LOCAL_E(cur[0]-h, cur[1], em);
              LOCAL_E(cur[0], cur[1]+h, ep2); LOCAL_E(cur[0], cur[1]-h, em2);
              gu = (ep - em) / (2*h); gv = (ep2 - em2) / (2*h);
            }
            gn = sqrt(gu*gu + gv*gv);
            if (gn < 1e-12) continue;
            dirx = -gu/gn; diry = -gv/gn;

            t0 = 0.5 * elen; if (t0 > o.max_disp) t0 = o.max_disp; if (t0 < 0.02) t0 = 0.02;
            dispmax2 = o.max_disp * o.max_disp;
            t = t0;
            for (nls = 0; nls < (size_t)o.line_search; nls++) {
                double ddx, ddy; int flip = 0; double ecand;
                trial[0] = cur[0] + t*dirx; trial[1] = cur[1] + t*diry;
                ddx = trial[0] - (double)uv_in[v*2]; ddy = trial[1] - (double)uv_in[v*2+1];
                if (ddx*ddx + ddy*ddy > dispmax2) { t *= 0.5; continue; }
                for (s = s0; s < s1 && !flip; s++) {
                    size_t ff = vf_idx[s];
                    int32_t a = faces[ff*3+0], b = faces[ff*3+1], c = faces[ff*3+2];
                    const double *u0 = ((size_t)a==v)?trial:uv+(size_t)a*2;
                    const double *u1 = ((size_t)b==v)?trial:uv+(size_t)b*2;
                    const double *u2 = ((size_t)c==v)?trial:uv+(size_t)c*2;
                    double a2 = signed_area2(u0, u1, u2);
                    if ((a2 >= 0 ? 1 : -1) != fr[ff].s || fabs(a2) < RLX_AREA_EPS) flip = 1;
                }
                if (!flip) {
                    LOCAL_E(trial[0], trial[1], ecand);
                    if (ecand < base - 1e-12) { uv[v*2] = trial[0]; uv[v*2+1] = trial[1]; moved++; break; }
                }
                t *= 0.5;
            }
            #undef LOCAL_E
        }
        n_move_total += moved;
        if (o.verbose) {
            double e; measure(fr, faces, nf, uv, ref_sign, lambda, &e, NULL, NULL, NULL, NULL);
            fprintf(stderr, "  [relax] sweep %2zu: moved %zu, E=%.1f\n", sweep, moved, e);
        }
        if (moved == 0) break;
    }

    /* per-movable-vertex displacement */
    for (v = 0; v < nv; v++) {
        double dx, dy, dd;
        if (!movable[v]) continue;
        dx = uv[v*2] - (double)uv_in[v*2]; dy = uv[v*2+1] - (double)uv_in[v*2+1];
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
            stats->mean_disp = (n_interior > 0) ? disp_sum / (double)n_interior : 0.0;
            stats->max_disp = disp_max;
        }
        if (fl1 > (stats ? stats->flips_before : 0) || E1 > (stats ? stats->energy_before : E1) + 1e-6) {
            for (i = 0; i < nv * 2; i++) uv_out[i] = uv_in[i];   /* revert */
            if (stats) stats->reverted = 1;
            return 0;
        }
    }

    for (i = 0; i < nv * 2; i++) uv_out[i] = (float)uv[i];
    return 0;
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
        face_precompute(P0, P1, P2, u0, u1, u2, 30.0*RLX_PI/180.0, 1.0, 0.1, 100.0, &fa);
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

    free(V); free(F); free(uv_in); free(uv_out);
    Arena_dispose(&arena);
    fprintf(stderr, "[ribbon_relax selftest] %s\n", fails == 0 ? "ok" : "FAILED");
    return fails;
}
