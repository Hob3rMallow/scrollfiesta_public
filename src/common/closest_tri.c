#include "closest_tri.h"
#include <string.h>

static double v3dot(const double a[3], const double b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

/* closest point on triangle abc to p (Ericson, Real-Time Collision Detection) */
void ClosestTri_point(const double p[3], const double a[3],
                      const double b[3], const double c[3],
                      double out[3], double *out_u, double *out_v)
{
    double ab[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
    double ac[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
    double ap[3] = { p[0]-a[0], p[1]-a[1], p[2]-a[2] };
    double d1 = v3dot(ab, ap), d2 = v3dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) { memcpy(out, a, 3*sizeof(double)); *out_u = 0; *out_v = 0; return; }
    double bp[3] = { p[0]-b[0], p[1]-b[1], p[2]-b[2] };
    double d3 = v3dot(ab, bp), d4 = v3dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) { memcpy(out, b, 3*sizeof(double)); *out_u = 1; *out_v = 0; return; }
    double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double w = d1 / (d1 - d3);
        for (int i = 0; i < 3; i++) out[i] = a[i] + w * ab[i];
        *out_u = w; *out_v = 0; return;
    }
    double cp[3] = { p[0]-c[0], p[1]-c[1], p[2]-c[2] };
    double d5 = v3dot(ab, cp), d6 = v3dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) { memcpy(out, c, 3*sizeof(double)); *out_u = 0; *out_v = 1; return; }
    double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        for (int i = 0; i < 3; i++) out[i] = a[i] + w * ac[i];
        *out_u = 0; *out_v = w; return;
    }
    double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        for (int i = 0; i < 3; i++) out[i] = b[i] + w * (c[i] - b[i]);
        *out_u = 1.0 - w; *out_v = w; return;
    }
    double denom = 1.0 / (va + vb + vc);
    double v_ = vb * denom, w_ = vc * denom;
    for (int i = 0; i < 3; i++) out[i] = a[i] + ab[i]*v_ + ac[i]*w_;
    *out_u = v_; *out_v = w_;
}

double ClosestTri_dist2(const double p[3], const double a[3],
                        const double b[3], const double c[3])
{
    double q[3], u, v;
    ClosestTri_point(p, a, b, c, q, &u, &v);
    double d0 = p[0]-q[0], d1 = p[1]-q[1], d2 = p[2]-q[2];
    return d0*d0 + d1*d1 + d2*d2;
}
