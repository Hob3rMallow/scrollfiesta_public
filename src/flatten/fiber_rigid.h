#ifndef FIBER_RIGID_INCLUDED
#define FIBER_RIGID_INCLUDED

#include "slim_refine.h"

typedef struct FiberRigidStats {
    size_t faces, supported_faces;
    double area, supported_area, evidence_mass;
    double radians, center[2];
    double rms_before_degrees, rms_after_degrees;
    double signed_mean_before_degrees, signed_mean_after_degrees;
} FiberRigidStats;

/* One common SO(2) turn for the entire input, including disconnected pieces,
 * point fans and isolated vertices. Minimize the physical-area/evidence weighted
 * squared angular distance of both mapped fibre arms to complementary axes.
 * The exact periodic scalar least-squares fit is a sorted quadratic sweep;
 * it has no iterative linear solve or topology-dependent admission rule.
 * Missing evidence is excluded and reported. A successful fit certifies only
 * its angular objective, not the input's flattening or material validity.
 * Return -1 for invalid geometry/directions or no usable evidence. */
int FiberRigid_fit(Arena_T arena, const float *xyz, size_t nv,
                   const int32_t *faces, size_t nf, const double *uv,
                   const SlimRefineDirection *directions, FiberRigidStats *stats);

/* Rotate every UV around the reported physical-area centroid. Input/output may
 * alias. No scaling, reflection, piece packing, mesh repair or XYZ change. */
int FiberRigid_apply(const double *uv, size_t nv, const FiberRigidStats *stats,
                     double *out);

int FiberRigid_selftest(void);
#endif
