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
// The functional is evaluated in plain double with a running forward error bound
// (2026-10-06; decides nearly every test without the rounding-mode window), then,
// only if that cannot certify the sign, as an attene_interval (directed-rounding
// filter), and only when the interval is unreliable, as an exact BSNumber. A single
// templated core computes (L,d,g) so the tiers can never diverge.

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
    int     ld_ok;     // interval (L,d) below is cached
    attene_interval Lx, Ly, Lz, d;   // interval homogeneous coordinates
    int     ed_ok;     // double (L,d) and error bounds below are cached
    double  edL[3], edLe[3], edD, edDe;
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

// ---- error-bounded double (the 2026-10-06 first tier) ---------------------------
// ED carries a round-to-nearest value v and a bound e >= |v - x|, x being the exact
// real value of the same expression over the same (exact, double) inputs.  Per
// operation (u = 2^-53): a sum adds the operands' bounds plus its own rounding,
// |fl(s) - s| <= u|s| <= 2u|fl(s)|; a product of a = av - da, b = bv - db adds
// |av| eb + |bv| ea + ea eb plus its rounding.  ED_TINY covers absolute rounding
// near underflow.  The bound arithmetic itself rounds (relative error ~ depth * u,
// depth < 100), which ED_SAFE absorbs ~10^5-fold.  A sign is certified only when
// |v| > e * ED_SAFE; NaN or Inf never certifies.
struct ED {
    double v, e;
    ED() : v(0.0), e(0.0) {}
    ED(double x) : v(x), e(0.0) {}
    ED(double x, double err) : v(x), e(err) {}
};
static const double ED_U    = 2.220446049250313e-16;   // 2^-52 = 2u
static const double ED_TINY = 1e-300;
static const double ED_SAFE = 1.0 + 1e-9;
static inline ED operator+(const ED& a, const ED& b) {
    double s = a.v + b.v;
    return ED(s, a.e + b.e + (ED_U * std::fabs(s) + ED_TINY));
}
static inline ED operator-(const ED& a, const ED& b) {
    double s = a.v - b.v;
    return ED(s, a.e + b.e + (ED_U * std::fabs(s) + ED_TINY));
}
static inline ED operator*(const ED& a, const ED& b) {
    double p = a.v * b.v;
    return ED(p, std::fabs(a.v) * b.e + std::fabs(b.v) * a.e + a.e * b.e +
                 (ED_U * std::fabs(p) + ED_TINY));
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

// The same g, from a point's (L,d) and the plane's precomputed
// qmp = q - p and pq = p.p - q.q: identical operations on identical operands.
template <class T>
static void side_g(const V3<T>& L, const T& d, const V3<T>& qmp, const T& pq, T& g) {
    const T two(2.0);
    g = two * vdot(L, qmp) + pq * d;
}

// Interval (L,d) of a point, computed by point_LD on its first side test and
// cached in the record (2026-10-04).  Called only inside the FE_UPWARD window,
// as the uncached evaluation was, so the cached bits are exactly what every
// later recomputation would produce: a constructed point is tested against
// several bisectors, and recomputing a facet point cost ~50 interval products.
static inline void cached_LD(const PtRec* pts, PtRec& v,
                             V3<attene_interval>& L, attene_interval& d) {
    if (!v.ld_ok) {
        point_LD<attene_interval>(pts, v, L, d);
        v.Lx = L.x; v.Ly = L.y; v.Lz = L.z; v.d = d;
        v.ld_ok = 1;
        return;
    }
    L = V3<attene_interval>{ v.Lx, v.Ly, v.Lz };
    d = v.d;
}

// Double (L,d) with error bounds, cached per constructed point like cached_LD.
// Runs under round-to-nearest (outside the interval window).
static inline void cached_LD_ed(const PtRec* pts, PtRec& v, V3<ED>& L, ED& d) {
    if (v.kind == PT_EXPLICIT) {
        L = V3<ED>{ ED(v.x), ED(v.y), ED(v.z) };
        d = ED(1.0);
        return;
    }
    if (!v.ed_ok) {
        point_LD<ED>(pts, v, L, d);
        v.edL[0] = L.x.v; v.edL[1] = L.y.v; v.edL[2] = L.z.v;
        v.edLe[0] = L.x.e; v.edLe[1] = L.y.e; v.edLe[2] = L.z.e;
        v.edD = d.v; v.edDe = d.e;
        v.ed_ok = 1;
        return;
    }
    L = V3<ED>{ ED(v.edL[0], v.edLe[0]), ED(v.edL[1], v.edLe[1]), ED(v.edL[2], v.edLe[2]) };
    d = ED(v.edD, v.edDe);
}

// First tier: the exact side sign when the double evaluation certifies both
// sign(g) and sign(d), else 2 (undecided).
static inline int side_double(const PtRec* pts, PtRec& v, const V3<ED>& qmp, const ED& pq) {
    V3<ED> L;
    ED d, g;
    cached_LD_ed(pts, v, L, d);
    side_g<ED>(L, d, qmp, pq, g);
    if (std::fabs(g.v) > g.e * ED_SAFE && std::fabs(d.v) > d.e * ED_SAFE)
        return ((g.v > 0.0) - (g.v < 0.0)) * ((d.v > 0.0) - (d.v < 0.0));
    return 2;
}

// Exact fallback (BSNumber, never rounds); run under round-to-nearest.
static int side_exact(const PtRec* pts, const PtRec& v, const double P[3], const double Q[3]) {
    BSNumber gb, db;
    side_gd<BSNumber>(pts, v, P, Q, gb, db);
    int dbs = db.sign();
    if (dbs == 0) return 0;          // degenerate construction (no unique point)
    return gb.sign() * dbs;
}

// The pre-2026-10-04 evaluation, kept verbatim as the selftest reference.
static int side_reference(const PtRec* pts, const PtRec& vr, const double P[3], const double Q[3]) {
    attene_setFPUModeToRoundUP();
    attene_interval gi, di;
    side_gd<attene_interval>(pts, vr, P, Q, gi, di);
    bool reliable = gi.signIsReliable() && di.signIsReliable();
    int gs = gi.sign(), ds = di.sign();
    attene_setFPUModeToRoundNEAR();
    if (reliable) return gs * ds;
    return side_exact(pts, vr, P, Q);
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
    int out = 0;
    AttenePred_side_bisector_batch(ctx, &v, 1, si, sj, &out);
    return out;
}

void AttenePred_side_bisector_batch(AttenePredCtx *ctx, const AttenePt *v, int n,
                                    AttenePt si, AttenePt sj, int *out) {
    PtRec *pts = ctx->pts.data();
    double P[3], Q[3];
    fetch(pts, si, P);
    fetch(pts, sj, Q);
    int need_interval = 0, need_exact = 0;

    // First tier (2026-10-06): error-bounded double, no rounding-mode switch.
    // Certified signs are exact, so the result equals the old two-tier path's.
    {
        V3<ED> p = vmk<ED>(P), q = vmk<ED>(Q);
        V3<ED> qmp = vsub(q, p);
        ED pq = vdot(p, p) - vdot(q, q);
        for (int i = 0; i < n; i++) {
            out[i] = side_double(pts, pts[v[i]], qmp, pq);
            if (out[i] == 2) need_interval = 1;
        }
    }
    if (!need_interval) return;

    // Second tier: interval filter under directed (round-up) arithmetic, one
    // rounding-mode window and one set of plane constants for the undecided.
    attene_setFPUModeToRoundUP();
    {
        V3<attene_interval> p = vmk<attene_interval>(P), q = vmk<attene_interval>(Q);
        V3<attene_interval> qmp = vsub(q, p);
        attene_interval pq = vdot(p, p) - vdot(q, q);
        for (int i = 0; i < n; i++) {
            if (out[i] != 2) continue;
            V3<attene_interval> L;
            attene_interval di, gi;
            cached_LD(pts, pts[v[i]], L, di);
            side_g<attene_interval>(L, di, qmp, pq, gi);
            if (gi.signIsReliable() && di.signIsReliable()) {
                out[i] = gi.sign() * di.sign();
            } else {
                out[i] = 3;                  // decided exactly below
                need_exact = 1;
            }
        }
    }
    attene_setFPUModeToRoundNEAR();
    if (!need_exact) return;

    // Exact fallback: BSNumber (arbitrary precision, never rounds).
    for (int i = 0; i < n; i++)
        if (out[i] == 3) out[i] = side_exact(pts, pts[v[i]], P, Q);
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

    {
        // Cached + batched evaluation against the verbatim reference, on random
        // sites with constructed edge/facet points, exact on-bisector points and
        // repeated tests of the same point (cache reuse).
        AttenePred_reset(c);
        uint32_t st = 0x9E3779B9u;
        auto rnd = [&st](double lo, double hi) {
            st ^= st << 13; st ^= st >> 17; st ^= st << 5;
            return lo + (hi - lo) * ((double)(st >> 8) * (1.0 / 16777216.0));
        };
        const int NS = 40;
        for (int i = 0; i < NS; i++)
            AttenePred_explicit(c, 1.0e4 + rnd(-6, 6), 2.5e3 + rnd(-6, 6), 3.1e3 + rnd(-6, 6));
        AttenePt tri[3];
        for (int k = 0; k < 3; k++)
            tri[k] = AttenePred_explicit(c, 1.0e4 + rnd(-6, 6), 2.5e3 + rnd(-6, 6), 3.1e3 + rnd(-1, 1));
        std::vector<AttenePt> probe;
        for (int i = 0; i < 300; i++) {
            int a = (int)rnd(0, NS - 0.01), b = (int)rnd(0, NS - 0.01), k = (int)rnd(0, NS - 0.01);
            if (a == b || a == k || b == k) continue;
            probe.push_back(AttenePred_edge_bisector(c, tri[i % 3], tri[(i + 1) % 3], a, b));
            probe.push_back(AttenePred_facet_bisector2(c, tri[0], tri[1], tri[2], a, b, k));
            probe.push_back(AttenePred_explicit(c, 1.0e4 + rnd(-6, 6), 2.5e3 + rnd(-6, 6), 3.1e3 + rnd(-6, 6)));
        }
        // constructed points that lie exactly on the bisector they are tested against
        for (int i = 0; i + 1 < NS; i += 2) {
            probe.push_back(AttenePred_edge_bisector(c, tri[0], tri[1], i, i + 1));
            probe.push_back(AttenePred_facet_bisector2(c, tri[0], tri[1], tri[2], i, i + 1, (i + 2) % NS));
        }
        long n_checks = 0, n_bad = 0, n_zero = 0;
        std::vector<int> got(probe.size());
        for (int rep2 = 0; rep2 < 3; rep2++) {            // rep 1,2 hit the cache
            for (int i = 0; i + 1 < NS; i++) {
                AttenePt a = (AttenePt)i, b = (AttenePt)(i + 1);
                AttenePred_side_bisector_batch(c, probe.data(), (int)probe.size(), a, b, got.data());
                double P[3], Q[3];
                fetch(c->pts.data(), a, P); fetch(c->pts.data(), b, Q);
                for (size_t j = 0; j < probe.size(); j++) {
                    int ref = side_reference(c->pts.data(), c->pts[(size_t)probe[j]], P, Q);
                    int one = AttenePred_side_bisector(c, probe[j], a, b);
                    n_checks++;
                    if (ref == 0) n_zero++;
                    if (got[j] != ref || one != ref) n_bad++;
                }
            }
        }
        CHECK(n_bad == 0, "cached/batched side != reference");
        CHECK(n_zero > 0, "no exactly-on-bisector case exercised");
        fprintf(stderr, "AttenePred_selftest: %ld cached/batched side tests == reference (%ld on-bisector)\n",
                n_checks, n_zero);

        // The double tier on its own: every certified sign must equal the
        // reference, and it must certify nearly every off-bisector test.
        long n_dec = 0, n_off = 0, n_dbad = 0;
        for (int i = 0; i + 1 < NS; i++) {
            double P[3], Q[3];
            fetch(c->pts.data(), (AttenePt)i, P); fetch(c->pts.data(), (AttenePt)(i + 1), Q);
            V3<ED> p = vmk<ED>(P), q = vmk<ED>(Q);
            V3<ED> qmp = vsub(q, p);
            ED pq = vdot(p, p) - vdot(q, q);
            for (size_t j = 0; j < probe.size(); j++) {
                int ref = side_reference(c->pts.data(), c->pts[(size_t)probe[j]], P, Q);
                int dbl = side_double(c->pts.data(), c->pts[(size_t)probe[j]], qmp, pq);
                if (ref != 0) n_off++;
                if (dbl != 2) { n_dec++; if (dbl != ref) n_dbad++; }
            }
        }
        CHECK(n_dbad == 0, "double tier certified a wrong sign");
        CHECK(n_dec * 100 >= n_off * 95, "double tier certifies under 95% of off-bisector tests");
        fprintf(stderr, "AttenePred_selftest: double tier certified %ld of %ld off-bisector tests, 0 wrong\n",
                n_dec, n_off);
    }

    {
        // Explicit points one ulp either side of an exact bisector, at world-sized
        // coordinates: the double tier must certify only correct signs.
        AttenePred_reset(c);
        long n_dec = 0, n_dbad = 0, n_tests = 0;
        uint32_t st = 0x2545F491u;
        auto rnd = [&st](double lo, double hi) {
            st ^= st << 13; st ^= st >> 17; st ^= st << 5;
            return lo + (hi - lo) * ((double)(st >> 8) * (1.0 / 16777216.0));
        };
        for (int t = 0; t < 2000; t++) {
            double cx = 1.0e4 + rnd(-50, 50), cy = 2.5e3 + rnd(-50, 50), cz = 3.1e3 + rnd(-50, 50);
            double h = (double)(1 + (int)rnd(0, 7.99)) * 0.25;   // exact half-gap
            AttenePt a = AttenePred_explicit(c, cx - h, cy, cz);
            AttenePt b = AttenePred_explicit(c, cx + h, cy, cz);  // bisector: x == cx exactly
            double P[3], Q[3];
            fetch(c->pts.data(), a, P); fetch(c->pts.data(), b, Q);
            V3<ED> p = vmk<ED>(P), q = vmk<ED>(Q);
            V3<ED> qmp = vsub(q, p);
            ED pq = vdot(p, p) - vdot(q, q);
            double xs[5] = { cx, std::nextafter(cx, 0.0), std::nextafter(cx, 1e9),
                             cx - 1e-9, cx + 1e-9 };
            for (int k = 0; k < 5; k++) {
                AttenePt pt = AttenePred_explicit(c, xs[k], cy + rnd(-9, 9), cz + rnd(-9, 9));
                int ref = side_reference(c->pts.data(), c->pts[(size_t)pt], P, Q);
                int dbl = side_double(c->pts.data(), c->pts[(size_t)pt], qmp, pq);
                n_tests++;
                if (dbl != 2) { n_dec++; if (dbl != ref) n_dbad++; }
            }
        }
        CHECK(n_dbad == 0, "double tier wrong one ulp from a bisector");
        fprintf(stderr, "AttenePred_selftest: %ld near-bisector tests, double tier certified %ld, 0 wrong\n",
                n_tests, n_dec);
    }

done:
    #undef CHECK
    AttenePred_free(c);
    if (rc == 0) fprintf(stderr, "AttenePred_selftest: OK\n");
    return rc;
}

} // extern "C"
