#ifndef VESUVIUS_SHAFT_WARP_H
#define VESUVIUS_SHAFT_WARP_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32) && defined(SHAFT_WARP_BUILD_DLL)
#define SHAFT_WARP_API __declspec(dllexport)
#else
#define SHAFT_WARP_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* An ordered natural cubic in physical micrometres. Scanner Z is never sorted.
 * Its chord parameter, physical arc and transported normal frame are distinct.
 * This immutable object may be queried concurrently after successful creation. */
typedef struct ShaftWarp ShaftWarp;

enum {
    SHAFT_WARP_ENDPOINT = 1u,             /* finite endpoint foot, diagnostic */
    SHAFT_WARP_AMBIGUOUS = 2u,            /* distant feet tied in distance */
    SHAFT_WARP_OUTSIDE_SUPPORT = 4u,      /* endpoint extrapolation requested */
    SHAFT_WARP_FOLDED = 8u,               /* nonpositive local Jacobian */
    SHAFT_WARP_OUTSIDE_RADIUS = 16u,
    SHAFT_WARP_BELOW_JACOBIAN = 32u,
    SHAFT_WARP_DIFFERENT_BRANCH = 64u
};

typedef struct {
    double parameter_um, s_um;
    double foot_um_zyx[3], tangent_zyx[3];
    double normal1_zyx[3], normal2_zyx[3], curvature_per_um_zyx[3];
    double distance_um, normal_residual_um, local_jacobian;
    uint32_t flags;
} ShaftWarpProjection;

SHAFT_WARP_API uint32_t ShaftWarp_abi(void);
SHAFT_WARP_API size_t ShaftWarp_projection_size(void);

/* out must point to NULL. initial_normal may be NULL for a deterministic frame.
 * Returns 0, -1 (invalid input), -2 (allocation), or -3 (unresolved geometry).
 * Failed creation leaves *out unchanged. Duplicate/stationary curves fail. */
SHAFT_WARP_API int ShaftWarp_create(ShaftWarp **out, const double *points_um_zyx,
                                   size_t count, const double initial_normal_zyx[3]);
SHAFT_WARP_API void ShaftWarp_free(ShaftWarp *warp);
SHAFT_WARP_API double ShaftWarp_length_um(const ShaftWarp *warp);
SHAFT_WARP_API double ShaftWarp_parameter_length_um(const ShaftWarp *warp);
SHAFT_WARP_API int ShaftWarp_eval_parameter(const ShaftWarp *warp, double parameter_um,
                                           ShaftWarpProjection *out);
SHAFT_WARP_API int ShaftWarp_eval_arc(const ShaftWarp *warp, double s_um,
                                     ShaftWarpProjection *out);
/* Conservative physical AABB of every normal disk of radius radius_um over
 * the closed arc interval. Uses cubic coordinate extrema, expanded by radius;
 * it can include points outside the tube. Bounds do not certify invertibility.
 * Require 0 <= begin <= end <= length and finite nonnegative radius. Return
 * 0, -1 (invalid input), or -3 (unresolved); failure leaves both outputs intact. */
SHAFT_WARP_API int ShaftWarp_bounds_arc(const ShaftWarp *warp, double begin_um,
    double end_um, double radius_um, double lower_um_zyx[3], double upper_um_zyx[3]);
/* Search all potentially nearest cubic spans. ambiguity_um is a numerical
 * distance-tie tolerance, not a shaft-accuracy bound. Flags are not suppressed. */
SHAFT_WARP_API int ShaftWarp_project(const ShaftWarp *warp, const double point_um_zyx[3],
                                    double ambiguity_um, ShaftWarpProjection *out);

/* Checked world <-> (arc,normal1,normal2) transforms, all in physical um.
 * Return 1 for refused coordinates and negative for errors; outputs remain
 * unchanged on failure. Exact endpoint disks are allowed, extrapolation is not.
 * The inverse also verifies the global nearest branch. Geometry checks do not
 * establish that the supplied shaft is accurate or supported by scroll RAW. */
SHAFT_WARP_API int ShaftWarp_to_metric(const ShaftWarp *warp, const double point_um_zyx[3],
    double maximum_radius_um, double minimum_jacobian, double ambiguity_um,
    double out_suv_um[3], uint32_t *flags);
SHAFT_WARP_API int ShaftWarp_from_metric(const ShaftWarp *warp, const double suv_um[3],
    double maximum_radius_um, double minimum_jacobian, double ambiguity_um,
    double out_point_um_zyx[3], uint32_t *flags);

#ifdef __cplusplus
}
#endif
#endif
