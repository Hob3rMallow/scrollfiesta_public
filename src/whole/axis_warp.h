#ifndef VESUVIUS_AXIS_WARP_H
#define VESUVIUS_AXIS_WARP_H

#include <stddef.h>
#include <stdint.h>
#include "shaft_warp.h"

typedef struct {
    double source_voxel_um_zyx[3], source_origin_um_zyx[3];
    double metric_voxel_um;
    double maximum_radius_um, minimum_jacobian, ambiguity_um;
    double initial_normal_zyx[3];
    int has_initial_normal;
} AxisWarpPhysical;

/* A legacy slice table or an ordered physical shaft used to straighten metric
 * coordinates. Source/world vertices are retained separately for CT sampling
 * and export. The physical mode replaces scanner Z with shaft arc length. */
typedef struct {
    double *z;
    double *y;
    double *x;
    size_t n;
    double reference_y;
    double reference_x;
    /* When physical is set, z/y/x store ordered physical um, not scanner-Z
     * lookup samples. source voxels <-> physical um is explicit. */
    ShaftWarp *physical;
    AxisWarpPhysical physical_config;
} AxisWarp;

void AxisWarp_init(AxisWarp *warp);
void AxisWarp_dispose(AxisWarp *warp);

/* CSV rows are z,y,x. Non-numeric headers/comments are skipped; retained rows
 * must have strictly increasing, finite z. The straight reference center is
 * the interpolated curve midpoint. */
int AxisWarp_load_csv(AxisWarp *warp, const char *path);

/* Build an ordered physical shaft. On failure warp is unchanged. Metric
 * coordinates are (s,u,v)/metric_voxel_um + (0,reference_y,reference_x).
 * The reference Y/X is an output-coordinate origin, not an estimated axis. */
int AxisWarp_create_physical(AxisWarp *warp, const double *points_um_zyx,
                            size_t count, const AxisWarpPhysical *config,
                            double reference_y, double reference_x);
int AxisWarp_to_metric(const AxisWarp *warp, const double source_zyx[3],
                       double metric_suv[3], uint32_t *flags);
int AxisWarp_to_world(const AxisWarp *warp, const double metric_suv[3],
                      double source_zyx[3], uint32_t *flags);
/* Stable semantic identity, including mode, units, geometry and domain.
 * Revision must change when coordinate behavior changes. */
uint64_t AxisWarp_fingerprint(const AxisWarp *warp);

int AxisWarp_valid(const AxisWarp *warp);
/* Legacy scanner-Z lookup only. Physical shafts return NaN: a world point
 * or physical arc is required to identify a position through a bend. */
void AxisWarp_eval(const AxisWarp *warp, double z, double *y, double *x);

/* Legacy mode preserves Z and translates horizontal slices. Physical mode
 * uses the checked transform above and writes NaN on refusal; consumers that
 * must handle refusals should call AxisWarp_to_metric. Input/output may alias. */
void AxisWarp_straighten_point(const AxisWarp *warp,
                               const float in_zyx[3], float out_zyx[3]);
void AxisWarp_straighten_vertices(const AxisWarp *warp,
                                  const float *in_zyx, float *out_zyx,
                                  size_t nvertices);

int AxisWarp_selftest(void);

#endif
