#include "closeup_fiber_volume.h"
#include "fiber_frame.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double cfv_dot(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

static void cfv_cross(const double a[3], const double b[3], double out[3])
{
    for (int k = 0; k < 3; k++) out[k] = a[(k+1)%3]*b[(k+2)%3]-a[(k+2)%3]*b[(k+1)%3];
}

static int cfv_unit(double a[3])
{
    double length = sqrt(cfv_dot(a,a));
    if (!(length > 1e-12) || !isfinite(length)) return -1;
    for (int k = 0; k < 3; k++) a[k] /= length;
    return 0;
}

static int cfv_face(const float *xyz, const float *normal, const int32_t *face,
                    double center[3], double first[3], double second[3],
                    double depth[3], double *area)
{
    double diagonal[3] = {0}, cross[3] = {0};
    size_t a = (size_t)face[0], b = (size_t)face[1], c = (size_t)face[2];
    for (int k = 0; k < 3; k++) {
        size_t q = (size_t)k;
        first[k] = (double)xyz[3*b+q]-xyz[3*a+q];
        diagonal[k] = (double)xyz[3*c+q]-xyz[3*a+q];
        center[k] = ((double)xyz[3*a+q]+xyz[3*b+q]+xyz[3*c+q])/3;
        depth[k] = normal ? (double)normal[3*a+q]+normal[3*b+q]+normal[3*c+q] : 0;
    }
    cfv_cross(first,diagonal,cross);
    *area = .5*sqrt(cfv_dot(cross,cross));
    if (!(*area > 0) || !isfinite(*area) || cfv_unit(first) || cfv_unit(cross)) return -1;
    cfv_cross(cross,first,second);
    if (cfv_unit(depth)) memcpy(depth,cross,3*sizeof(double));
    return 0;
}

int CloseupFiberVolume_measure(const float *xyz, const float *normal, size_t nv,
                              const int32_t *faces, size_t nf, const char *field,
                              SlimRefineDirection *out, CloseupFiberStats *stats)
{
    double started = ves_clock_sec();
    CloseupFiberStats result = {0};
    FiberVolume_T volume = NULL;
    int rc = -1;
    const int depth_radius = FIBER_FLATTEN_HALF_DEPTH;
    if (stats) memset(stats,0,sizeof *stats);
    if (!nf) return 0;
    if (!xyz || !faces || !field || !out || !nv || nv > INT32_MAX ||
            nv > SIZE_MAX/(3*sizeof(float)) || nf > SIZE_MAX/(3*sizeof(int32_t)) ||
            nf > SIZE_MAX/sizeof *out) return -1;
    memset(out,0,nf*sizeof *out);
    for (size_t k = 0; k < 3*nv; k++)
        if (!isfinite(xyz[k]) || (normal && !isfinite(normal[k]))) goto done;
    for (size_t k = 0; k < 3*nf; k++)
        if (faces[k] < 0 || (size_t)faces[k] >= nv) goto done;
    if (FiberVolume_open(field,&volume)) goto done;
    result.patches = nf;
    for (size_t f = 0; f < nf; f++) {
        double center[3] = {0}, first[3] = {0}, second[3] = {0}, depth[3] = {0};
        double area = 0, real = 0, imaginary = 0, evidence = 0;
        if (cfv_face(xyz,normal,faces+3*f,center,first,second,depth,&area)) goto done;
        result.area += area;
        for (int offset = -depth_radius; offset <= depth_radius; offset++) {
            double position[3] = {0}, decoded[3] = {0}, direction[3] = {0};
            for (int k = 0; k < 3; k++) {
                position[k] = center[k]+(double)offset*depth[k];
                /* FiberVolume_sample rounds to C long before the bounds check. */
                if (!isfinite(position[k]) || position[k] < (double)LONG_MIN+1 ||
                        position[k] > (double)LONG_MAX-1) goto done;
            }
            double presence = FiberVolume_sample(volume,position[2],position[1],position[0],decoded);
            if (!isfinite(presence) || presence < 0 || presence > 1) goto done;
            if (!(presence > 0)) continue;
            for (int k = 0; k < 3; k++) direction[k] = decoded[2-k];
            if (cfv_unit(direction)) continue;
            double a = cfv_dot(direction,first), b = cfv_dot(direction,second);
            double projected2 = a*a+b*b;
            if (!(projected2 > 1e-20) || !isfinite(projected2)) continue;
            double weight = presence*fmin(1.0,projected2);
            double denominator = projected2*projected2;
            /* exp(4 i psi), independent of sign and of the two depth families. */
            real += weight*(a*a*a*a-6*a*a*b*b+b*b*b*b)/denominator;
            imaginary += weight*4*a*b*(a*a-b*b)/denominator;
            evidence += weight;
            result.samples++;
        }
        if (!(evidence > 0)) continue;
        double coherence = fmin(1.0,hypot(real,imaginary)/evidence);
        if (!isfinite(coherence) || coherence < ASM_FIBRE_MIN_COHERENCE) continue;
        double angle = .25*atan2(imaginary,real), cosine = cos(angle), sine = sin(angle);
        for (int k = 0; k < 3; k++) {
            out[f].h[k] = cosine*first[k]+sine*second[k];
            out[f].v[k] = -sine*first[k]+cosine*second[k];
        }
        out[f].weight = coherence*evidence/(double)(2*depth_radius+1);
        result.faces++; result.measured_patches++; result.measured_area += area;
        if ((f+1)%262144 == 0 || f+1 == nf) {
            fprintf(stderr,"[closeup fiber volume] faces %zu/%zu measured=%zu area=%.2f%% samples=%zu\n",
                    f+1,nf,result.faces,100*result.measured_area/result.area,result.samples);
            fflush(stderr);
        }
    }
    rc = 0;
done:
    result.seconds = ves_clock_sec()-started;
    if (stats) *stats = result;
    FiberVolume_close(&volume);
    if (rc) fprintf(stderr,"ERROR: direct physical fiber-volume measurement failed\n");
    return rc;
}

#include "closeup_fiber_volume_test.inc"
