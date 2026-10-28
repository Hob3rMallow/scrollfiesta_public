/* Exact generalized-metric signs, exposed to the native C optimizer.
 * If X and U are the physical and UV edge Gram matrices, respectively,
 * sigma_min >= l iff U-l*l*X is PSD; sigma_max <= h iff h*h*X-U is PSD.
 * The three principal minors of a symmetric 2x2 matrix decide PSD, without
 * a tangent frame, division, square root, or subtraction of absolute UVs
 * against rounded shape derivatives. The source Gram must be positive definite.
 *
 * The filter and exact fallback evaluate the same polynomial. This wrapper
 * owns any temporary allocations made inside the vendored exact library.
 * Build with IEEE arithmetic, without fast-math reassociation.
 */
#include "metric_bounds.h"
#include "bsnumber.h"
#include <cfloat>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <new>
#if defined(__SSE2__) || defined(_M_X64)
#include <xmmintrin.h>
#define MB_HAS_SSE 1
#else
#define MB_HAS_SSE 0
#endif

namespace {
struct MbFilter {
    double value, error;
    explicit MbFilter(float v) : value(v), error(0) {
        uint32_t bits; std::memcpy(&bits, &v, sizeof bits);
        /* Detect binary32 subnormals before promotion can flush them. */
        if ((bits & UINT32_C(0x7f800000)) == 0 &&
            (bits & UINT32_C(0x007fffff)) != 0)
            error = std::numeric_limits<double>::infinity();
    }
    explicit MbFilter(double v = 0) : value(v), error(0) {
        uint64_t bits; std::memcpy(&bits, &v, sizeof bits);
        /* Subnormal inputs go directly to exact arithmetic. This also keeps
         * the filter conservative when a caller has enabled DAZ/FTZ. */
        if ((bits & UINT64_C(0x7ff0000000000000)) == 0 &&
            (bits & UINT64_C(0x000fffffffffffff)) != 0)
            error = std::numeric_limits<double>::infinity();
    }
    MbFilter(double v, double e) : value(v), error(e) {}
};

static double mb_round(double v)
{
    /* Covers rounding in any IEEE mode, gradual or flushed underflow.
     * The factors below also bound rounding while accumulating the errors. */
    return 2*DBL_EPSILON*std::fabs(v) + 8*DBL_MIN;
}
static MbFilter operator+(const MbFilter& a, const MbFilter& b)
{
    double v = a.value+b.value;
    return MbFilter(v, 4*(a.error+b.error+mb_round(v)));
}
static MbFilter operator-(const MbFilter& a, const MbFilter& b)
{
    double v = a.value-b.value;
    return MbFilter(v, 4*(a.error+b.error+mb_round(v)));
}
static MbFilter operator*(const MbFilter& a, const MbFilter& b)
{
    double v = a.value*b.value;
    double e = std::fabs(a.value)*b.error + std::fabs(b.value)*a.error + a.error*b.error;
    return MbFilter(v, 4*(e+mb_round(v)));
}

/* Return 2 for an uncertain filtered sign; the exact tier returns -1,0,1. */
static int mb_sign(const MbFilter& a)
{
    if (!std::isfinite(a.value) || !std::isfinite(a.error)) return 2;
    return a.value > a.error ? 1 : a.value < -a.error ? -1 : 2;
}
static int mb_sign(const BSNumber& a) { return a.sign(); }

template<class T> struct MbGram { T a, b, c; };
template<class T, class Input> static MbGram<T> mb_gram(const Input *p, int dimension)
{
    MbGram<T> g;
    for (int k = 0; k < dimension; k++) {
        T a = T(p[dimension+k])-T(p[k]);
        T b = T(p[2*dimension+k])-T(p[k]);
        g.a = g.a+a*a; g.b = g.b+a*b; g.c = g.c+b*b;
    }
    return g;
}

template<class T> static int mb_classify(const float xyz[9], const double uv[6],
    const double bounds[4], unsigned *violated)
{
    MbGram<T> x = mb_gram<T>(xyz,3), u = mb_gram<T>(uv,2);
    int source = mb_sign(x.a*x.c-x.b*x.b);
    if (source == 2) return 1;
    if (source <= 0) return -1;
    unsigned bad = 0;
    for (int k = 0; k < 4; k++) {
        T scale = T(bounds[k])*T(bounds[k]);
        T a = u.a-scale*x.a, b = u.b-scale*x.b, c = u.c-scale*x.c;
        int sa = mb_sign(a), sc = mb_sign(c), sd = mb_sign(a*c-b*b);
        /* Upper bounds negate both diagonals; the determinant is unchanged. */
        if ((k & 1) && sa != 2) sa = -sa;
        if ((k & 1) && sc != 2) sc = -sc;
        if (sa == -1 || sc == -1 || sd == -1) bad |= 1u << k;
        else if (sa == 2 || sc == 2 || sd == 2) return 1;
    }
    *violated = bad;
    return 0;
}

static bool mb_valid(const float xyz[9], const double uv[6], const double bounds[4])
{
    if (!xyz || !uv || !bounds) return false;
    for (int k = 0; k < 9; k++) if (!std::isfinite(xyz[k])) return false;
    for (int k = 0; k < 6; k++) if (!std::isfinite(uv[k])) return false;
    uint64_t ordered[4]; std::memcpy(ordered,bounds,sizeof ordered);
    /* Positive finite binary64 encodings have the same integer order, even
     * when DAZ would make a comparison treat a subnormal bound as zero. */
    for (int k = 0; k < 4; k++)
        if (!ordered[k] || ordered[k] >= UINT64_C(0x7ff0000000000000)) return false;
    return ordered[0] <= ordered[1] && ordered[2] <= ordered[3];
}
}

extern "C" int MetricBounds_classify(const float xyz[9], const double uv[6],
    const double bounds[4], MetricBoundsResult *result)
{
    if (!result || !mb_valid(xyz,uv,bounds)) return -1;
    try {
        MetricBoundsResult r = {0,0};
        int rc = mb_classify<MbFilter>(xyz,uv,bounds,&r.violated);
        if (rc > 0) {
            r.exact_fallback = 1;
            rc = mb_classify<BSNumber>(xyz,uv,bounds,&r.violated);
        }
        if (rc) return -1;
        *result = r;
        return 0;
    } catch (const std::bad_alloc&) { return -1; }
}

#include "metric_bounds_test.inc"
