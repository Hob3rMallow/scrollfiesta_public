/* Necessary corner-angle bounds for a distortion-limited planar map.
 *
 * A linear triangle map whose principal stretches have ratio at most K sends
 * a corner theta into [2 atan(tan(theta/2)/K), 2 atan(K tan(theta/2))].
 * For a locally injective interior vertex, the mapped corners sum to 2*pi.
 * Failing this interval test proves that the original fan needs a cut or a
 * geometry repair; passing it is only a necessary condition, not a UV audit.
 */
#ifndef INTRINSIC_ANGLE_INCLUDED
#define INTRINSIC_ANGLE_INCLUDED

#include <math.h>

#define INTRINSIC_ANGLE_TAU 6.28318530717958647693
#define INTRINSIC_ANGLE_TOL 1e-9

static inline int IntrinsicAngle_corner_bounds(const float *a, const float *b,
    const float *c, double condition, double *lower, double *upper)
{
    double u[3],v[3],cross[3],dot=0.0;
    for(int d=0;d<3;d++) {
        u[d]=(double)b[d]-a[d]; v[d]=(double)c[d]-a[d];
        dot+=u[d]*v[d];
    }
    cross[0]=u[1]*v[2]-u[2]*v[1];
    cross[1]=u[2]*v[0]-u[0]*v[2];
    cross[2]=u[0]*v[1]-u[1]*v[0];
    double area2=hypot(hypot(cross[0],cross[1]),cross[2]);
    if(!(condition>=1.0)||!isfinite(condition)||!(area2>0.0)||!isfinite(dot)||!isfinite(area2))
        return 0;
    double half=0.5*atan2(area2,dot),s=sin(half),co=cos(half);
    *lower=2.0*atan2(s,condition*co);
    *upper=2.0*atan2(condition*s,co);
    return isfinite(*lower)&&isfinite(*upper);
}

static inline int IntrinsicAngle_closed_fan_possible(double lower,double upper)
{
    return isfinite(lower)&&isfinite(upper)&&
        lower<=INTRINSIC_ANGLE_TAU+INTRINSIC_ANGLE_TOL&&
        upper>=INTRINSIC_ANGLE_TAU-INTRINSIC_ANGLE_TOL;
}
#endif
