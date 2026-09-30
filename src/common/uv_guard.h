#ifndef UV_GUARD_INCLUDED
#define UV_GUARD_INCLUDED

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

/* Strict orientation in the face's own translated frame. The arithmetic
 * margin is relative to products, not chart position or physical scale. */
static inline int UvGuard_positive(const double *uv, const int32_t *f)
{
    double ax = uv[2*(size_t)f[1]]-uv[2*(size_t)f[0]];
    double ay = uv[2*(size_t)f[1]+1]-uv[2*(size_t)f[0]+1];
    double bx = uv[2*(size_t)f[2]]-uv[2*(size_t)f[0]];
    double by = uv[2*(size_t)f[2]+1]-uv[2*(size_t)f[0]+1];
    double p = ax*by, q = ay*bx;
    return isfinite(p) && isfinite(q) && p-q > 32*DBL_EPSILON*(fabs(p)+fabs(q));
}

/* Protect the whole straight path, not just the endpoint: an orientation-
 * preserving endpoint can lie beyond two determinant roots. The quadratic's
 * minimum on [0,alpha] is evaluated explicitly. */
static inline int UvGuard_interval(const double *uv, const double *target,
                                  const int32_t *f, double alpha)
{
    double e[2], g[2], de[2], dg[2];
    if (!UvGuard_positive(uv,f) || !(alpha >= 0 && alpha <= 1)) return 0;
    for (size_t k = 0; k < 2; k++) {
        size_t a = 2*(size_t)f[0]+k, b = 2*(size_t)f[1]+k, c = 2*(size_t)f[2]+k;
        e[k] = uv[b]-uv[a]; g[k] = uv[c]-uv[a];
        de[k] = (target[b]-target[a])-e[k];
        dg[k] = (target[c]-target[a])-g[k];
    }
    double a = de[0]*dg[1]-de[1]*dg[0];
    double b = de[0]*g[1]-de[1]*g[0]+e[0]*dg[1]-e[1]*dg[0];
    double c = e[0]*g[1]-e[1]*g[0];
    double minimum = fmin(c,(a*alpha+b)*alpha+c);
    if (a > 0) {
        double t = -.5*b/a;
        if (t > 0 && t < alpha) minimum = fmin(minimum,(a*t+b)*t+c);
    }
    double scale = fabs(c)+fabs(b)*alpha+fabs(a)*alpha*alpha;
    return isfinite(scale) && minimum > 64*DBL_EPSILON*scale;
}

static inline int UvGuard_mesh(const double *uv, const int32_t *faces, size_t nf)
{
    for (size_t f = 0; f < nf; f++) if (!UvGuard_positive(uv,faces+3*f)) return 0;
    return 1;
}

#endif
