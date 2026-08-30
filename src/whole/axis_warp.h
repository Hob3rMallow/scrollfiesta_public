#ifndef VESUVIUS_AXIS_WARP_H
#define VESUVIUS_AXIS_WARP_H

#include <stddef.h>

/* A sampled z,y,x umbilicus used to straighten only the metric coordinates
 * seen by the whole-scroll parameterizer.  Source/world vertices are retained
 * separately for CT sampling and export. */
typedef struct {
    double *z;
    double *y;
    double *x;
    size_t n;
    double reference_y;
    double reference_x;
} AxisWarp;

void AxisWarp_init(AxisWarp *warp);
void AxisWarp_dispose(AxisWarp *warp);

/* CSV rows are z,y,x. Non-numeric headers/comments are skipped; retained rows
 * must have strictly increasing, finite z. The straight reference center is
 * the interpolated curve midpoint. */
int AxisWarp_load_csv(AxisWarp *warp, const char *path);

int AxisWarp_valid(const AxisWarp *warp);
void AxisWarp_eval(const AxisWarp *warp, double z, double *y, double *x);

/* Preserve z and translate each horizontal slice so the sampled curve maps to
 * (reference_y, reference_x). Input and output may alias. */
void AxisWarp_straighten_point(const AxisWarp *warp,
                               const float in_zyx[3], float out_zyx[3]);
void AxisWarp_straighten_vertices(const AxisWarp *warp,
                                  const float *in_zyx, float *out_zyx,
                                  size_t nvertices);

int AxisWarp_selftest(void);

#endif
