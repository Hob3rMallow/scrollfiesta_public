// attene_predicates_wrap.cpp — exact Voronoi "side" predicates for restricted
// power-diagram / RVD clipping, C-callable, backed by the vendored Attene machinery.
//
// A Voronoi bisector plane is defined by two SITE points, not three mesh points, so
// it does not fit the vendored AtteneLPI3D/AtteneTPI3D (which model planes as a mesh
// triangle RST). We therefore evaluate the equidistance functional directly. For a
// point v with exact homogeneous coordinates (L/d):
//
//   |v-p|^2 - |v-q|^2 = ( 2*L.(q-p) + (p.p - q.q)*d ) / d
//
// so  side(v ; p,q) = sign(g)*sign(d)  with  g = 2*L.(q-p) + (p.p - q.q)*d.
// sign < 0  =>  v closer to p (inside p's half-space).
//
// (L,d) for v:
//   explicit           : L = v,            d = 1
//   edge x bisector(i,j): segment (a,b) meeting plane n.x = e, n = 2(sj-si),
//                         e = |sj|^2-|si|^2 ; d = n.(b-a), L = a*d + (e - n.a)*(b-a)
//   facet x bis(i,j),bis(i,k): three-plane meet (facet RST plane + two bisectors),
//                         Cramer:  d = m.(n1 x n2),
//                         L = f*(n1 x n2) + e1*(n2 x m) + e2*(m x n1)
//
// The functional is evaluated once as an attene_interval (fast, directed-rounding
// filter) and, only when the sign is unreliable, again as an exact BSNumber. A single
// templated core computes (L,d,g) so the two tiers can never diverge.

#include "attene_predicates_wrap.h"
#include "attene_predicates.h"   // attene_interval, attene_setFPUModeToRoundUP/NEAR
#include "bsnumber.h"            // BSNumber (exact)

#include <vector>
#include <cstdio>
#include <cmath>

namespace {

enum { PT_EXPLICIT = 0, PT_EDGE_BISECTOR = 1, PT_FACET_BISECTOR2 = 2 };

struct PtRec {
    int     kind;
    double  x, y, z;   // explicit coordinates (PT_EXPLICIT)
    int32_t id[6];     // defining EXPLICIT ids:
                       //   edge  : {a, b, si, sj}          -> id[0..3]
                       //   facet : {r, s, t, si, sj, sk}   -> id[0..5]
};

// ---- generic 3-vector over T in {attene_interval, BSNumber, double} ------------

template <class T> struct V3 { T x, y, z; };

template <class T> static inline V3<T> vmk(const double* p) {
    return V3<T>{ T(p[0]), T(p[1]), T(p[2]) };
}
template <class T> static inline V3<T> vsub(const V3<T>& a, const V3<T>& b) {
    return V3<T>{ a.x - b.x, a.y - b.y, a.z - b.z };
}
template <class T> static inline V3<T> vadd(const V3<T>& a, const V3<T>& b) {
    return V3<T>{ a.x + b.x, a.y + b.y, a.z + b.z };
}
template <class T> static inline T vdot(const V3<T>& a, const V3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <class T> static inline V3<T> vcross(const V3<T>& a, const V3<T>& b) {
    return V3<T>{ a.y * b.z - a.z * b.y,
                  a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x };
}
template <class T> static inline V3<T> vscale(const T& s, const V3<T>& a) {
    return V3<T>{ s * a.x, s * a.y, s * a.z };
}

static inline void fetch(const PtRec* pts, int32_t id, double out[3]) {
    const PtRec& r = pts[id];
    out[0] = r.x; out[1] = r.y; out[2] = r.z;
}

// Homogeneous (L, d) of a point, in arithmetic type T.
template <class T>
static void point_LD(const PtRec* pts, const PtRec& v, V3<T>& L, T& d) {
    if (v.kind == PT_EXPLICIT) {
        L = V3<T>{ T(v.x), T(v.y), T(v.z) };
        d = T(1.0);
        return;
    }

    const T two(2.0);

    if (v.kind == PT_EDGE_BISECTOR) {
        double A[3], B[3], SI[3], SJ[3];
        fetch(pts, v.id[0], A); fetch(pts, v.id[1], B);
        fetch(pts, v.id[2], SI); fetch(pts, v.id[3], SJ);
        V3<T> a = vmk<T>(A), b = vmk<T>(B), si = vmk<T>(SI), sj = vmk<T>(SJ);

        V3<T> n   = vscale(two, vsub(sj, si));       // n = 2(sj - si)
        T     e   = vdot(sj, sj) - vdot(si, si);     // e = |sj|^2 - |si|^2
        V3<T> bma = vsub(b, a);
        d = vdot(n, bma);                            // d = n.(b - a)
        T ena = e - vdot(n, a);                      // e - n.a
        L = vadd(vscale(d, a), vscale(ena, bma));    // L = a*d + (e - n.a)*(b - a)
        return;
    }

    // PT_FACET_BISECTOR2: facet plane (r,s,t) meets bisector(i,j) and bisector(i,k)
    double R[3], S[3], TT[3], SI[3], SJ[3], SK[3];
    fetch(pts, v.id[0], R); fetch(pts, v.id[1], S); fetch(pts, v.id[2], TT);
    fetch(pts, v.id[3], SI); fetch(pts, v.id[4], SJ); fetch(pts, v.id[5], SK);
    V3<T> r = vmk<T>(R), s = vmk<T>(S), t = vmk<T>(TT);
    V3<T> si = vmk<T>(SI), sj = vmk<T>(SJ), sk = vmk<T>(SK);

    V3<T> m  = vcross(vsub(s, r), vsub(t, r));       // facet normal
    T     f  = vdot(m, r);                           // facet offset (m.x = f)
    V3<T> n1 = vscale(two, vsub(sj, si)); T e1 = vdot(sj, sj) - vdot(si, si);
    V3<T> n2 = vscale(two, vsub(sk, si)); T e2 = vdot(sk, sk) - vdot(si, si);

    V3<T> n1xn2 = vcross(n1, n2);
    V3<T> n2xm  = vcross(n2, m);
    V3<T> mxn1  = vcross(m, n1);
    d = vdot(m, n1xn2);
    L = vadd(vadd(vscale(f, n1xn2), vscale(e1, n2xm)), vscale(e2, mxn1));
}

// g = 2*L.(q - p) + (p.p - q.q)*d , the numerator of |v-p|^2 - |v-q|^2.
template <class T>
static void side_gd(const PtRec* pts, const PtRec& v,
                    const double P[3], const double Q[3], T& g, T& d) {
    V3<T> L;
    point_LD<T>(pts, v, L, d);
    V3<T> p = vmk<T>(P), q = vmk<T>(Q);
    const T two(2.0);
    g = two * vdot(L, vsub(q, p)) + (vdot(p, p) - vdot(q, q)) * d;
}

} // namespace

struct AttenePredCtx {
    std::vector<PtRec> pts;
};

extern "C" {

AttenePredCtx *AttenePred_new(void) { return new AttenePredCtx(); }
void           AttenePred_free(AttenePredCtx *ctx) { delete ctx; }
void           AttenePred_reset(AttenePredCtx *ctx) { ctx->pts.clear(); }

int32_t AttenePred_count(AttenePredCtx *ctx) { return (int32_t)ctx->pts.size(); }
void    AttenePred_truncate(AttenePredCtx *ctx, int32_t n) {
    if (n >= 0 && (size_t)n < ctx->pts.size()) ctx->pts.resize((size_t)n);
}

AttenePt AttenePred_explicit(AttenePredCtx *ctx, double x, double y, double z) {
    PtRec r{};
    r.kind = PT_EXPLICIT; r.x = x; r.y = y; r.z = z;
    ctx->pts.push_back(r);
    return (AttenePt)(ctx->pts.size() - 1);
}

AttenePt AttenePred_edge_bisector(AttenePredCtx *ctx,
                                  AttenePt a, AttenePt b, AttenePt si, AttenePt sj) {
    PtRec r{};
    r.kind = PT_EDGE_BISECTOR;
    r.id[0] = a; r.id[1] = b; r.id[2] = si; r.id[3] = sj;
    ctx->pts.push_back(r);
    return (AttenePt)(ctx->pts.size() - 1);
}

AttenePt AttenePred_facet_bisector2(AttenePredCtx *ctx,
                                    AttenePt r_, AttenePt s_, AttenePt t_,
                                    AttenePt si, AttenePt sj, AttenePt sk) {
    PtRec r{};
    r.kind = PT_FACET_BISECTOR2;
    r.id[0] = r_; r.id[1] = s_; r.id[2] = t_; r.id[3] = si; r.id[4] = sj; r.id[5] = sk;
    ctx->pts.push_back(r);
    return (AttenePt)(ctx->pts.size() - 1);
}

int AttenePred_side_bisector(AttenePredCtx *ctx, AttenePt v, AttenePt si, AttenePt sj) {
    const PtRec *pts = ctx->pts.data();
    const PtRec &vr = pts[v];
    double P[3], Q[3];
    fetch(pts, si, P);
    fetch(pts, sj, Q);

    // Fast path: interval filter under directed (round-up) arithmetic.
    attene_setFPUModeToRoundUP();
    attene_interval gi, di;
    side_gd<attene_interval>(pts, vr, P, Q, gi, di);
    bool reliable = gi.signIsReliable() && di.signIsReliable();
    int gs = gi.sign(), ds = di.sign();
    attene_setFPUModeToRoundNEAR();
    if (reliable) return gs * ds;

    // Exact fallback: BSNumber (arbitrary precision, never rounds).
    BSNumber gb, db;
    side_gd<BSNumber>(pts, vr, P, Q, gb, db);
    int dbs = db.sign();
    if (dbs == 0) return 0;          // degenerate construction (no unique point)
    return gb.sign() * dbs;
}

int AttenePred_xyz(AttenePredCtx *ctx, AttenePt v, double *x, double *y, double *z) {
    const PtRec *pts = ctx->pts.data();
    const PtRec &vr = pts[v];
    if (vr.kind == PT_EXPLICIT) { *x = vr.x; *y = vr.y; *z = vr.z; return 0; }

    V3<double> L;
    double d;
    point_LD<double>(pts, vr, L, d);
    if (d == 0.0 || !std::isfinite(d)) return -1;
    *x = L.x / d; *y = L.y / d; *z = L.z / d;
    return 0;
}

// ---- self-test -----------------------------------------------------------------

int AttenePred_selftest(void) {
    AttenePredCtx *c = AttenePred_new();
    int rc = 0;
    #define CHECK(cond, msg) do { if (!(cond)) { \
        fprintf(stderr, "AttenePred_selftest FAIL: %s\n", (msg)); rc = -1; goto done; } } while (0)

    {
        // sites si=(0,0,0), sj=(2,0,0) -> bisector plane x = 1.
        AttenePt si = AttenePred_explicit(c, 0, 0, 0);
        AttenePt sj = AttenePred_explicit(c, 2, 0, 0);
        AttenePt sk = AttenePred_explicit(c, 0, 2, 0);   // for a second bisector

        AttenePt near_i = AttenePred_explicit(c, 0.5, 0.3, -0.2);
        AttenePt near_j = AttenePred_explicit(c, 1.5, -0.4, 0.1);
        AttenePt on_bis = AttenePred_explicit(c, 1.0, 5.0, -3.0);   // x==1 -> equidistant

        CHECK(AttenePred_side_bisector(c, near_i, si, sj) < 0, "explicit near si");
        CHECK(AttenePred_side_bisector(c, near_j, si, sj) > 0, "explicit near sj");
        CHECK(AttenePred_side_bisector(c, on_bis, si, sj) == 0, "explicit on bisector");
        // swapping the sites flips the sign
        CHECK(AttenePred_side_bisector(c, near_i, sj, si) > 0, "explicit near si (swapped)");

        // edge (a=si .. b=sj) crosses bisector(si,sj) exactly at (1,0,0):
        AttenePt ev = AttenePred_edge_bisector(c, si, sj, si, sj);
        CHECK(AttenePred_side_bisector(c, ev, si, sj) == 0, "edge pt on its own bisector");
        // (1,0,0): |.-si|^2 = 1, |.-sk|^2 = 1+4 = 5 -> closer to si
        CHECK(AttenePred_side_bisector(c, ev, si, sk) < 0, "edge pt vs bisector(si,sk)");
        double ex, ey, ez;
        CHECK(AttenePred_xyz(c, ev, &ex, &ey, &ez) == 0 &&
              std::fabs(ex - 1.0) < 1e-9 && std::fabs(ey) < 1e-9 && std::fabs(ez) < 1e-9,
              "edge pt coords");

        // facet z=0 plane (si,sj,sk all have z=0) meets bisector(si,sj) [x=1] and
        // bisector(si,sk) [y=1] at (1,1,0).
        AttenePt fv = AttenePred_facet_bisector2(c, si, sj, sk, si, sj, sk);
        double fx, fy, fz;
        CHECK(AttenePred_xyz(c, fv, &fx, &fy, &fz) == 0 &&
              std::fabs(fx - 1.0) < 1e-9 && std::fabs(fy - 1.0) < 1e-9 && std::fabs(fz) < 1e-9,
              "facet pt coords");
        CHECK(AttenePred_side_bisector(c, fv, si, sj) == 0, "facet pt on bisector(si,sj)");
        CHECK(AttenePred_side_bisector(c, fv, si, sk) == 0, "facet pt on bisector(si,sk)");
    }

    {
        // Large-coordinate stress: sites far from origin, a point one unit inside.
        AttenePt si = AttenePred_explicit(c, 1e7, 1e7, 1e7);
        AttenePt sj = AttenePred_explicit(c, 1e7 + 4.0, 1e7, 1e7);   // bisector x = 1e7+2
        AttenePt vi = AttenePred_explicit(c, 1e7 + 1.9, 1e7, 1e7);   // closer to si
        AttenePt vj = AttenePred_explicit(c, 1e7 + 2.1, 1e7, 1e7);   // closer to sj
        CHECK(AttenePred_side_bisector(c, vi, si, sj) < 0, "large-coord near si");
        CHECK(AttenePred_side_bisector(c, vj, si, sj) > 0, "large-coord near sj");
    }

done:
    #undef CHECK
    AttenePred_free(c);
    if (rc == 0) fprintf(stderr, "AttenePred_selftest: OK\n");
    return rc;
}

} // extern "C"
