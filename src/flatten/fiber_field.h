#ifndef FIBER_FIELD_INCLUDED
#define FIBER_FIELD_INCLUDED

#include <stdint.h>
#include "../common/arena.h"

/* ============================================================================
 * fiber_field.h -- papyrus fiber orientation field via the structure tensor.
 *
 * Given a baked rawtex texture (the (u,v) raster of RAW CT intensity over the
 * unrolled sheet), estimate at every texel the dominant local fiber
 * orientation and how strongly oriented it is (coherence). Papyrus is a
 * near-orthogonal cross-hatch, so a well-unrolled sheet has fibers running
 * along the u and v axes. This field is the target the ribbon-relax stage
 * aligns to (horizontal fibers -> u, vertical -> v), gated by coherence so
 * cracks/holes/damage (low coherence) never drive the alignment.
 *
 * Two symmetry modes (`rosy`):
 *   rosy = 2  LINE orientation (mod 180 deg). Plain Bigun structure tensor:
 *             gradient orientation phi = 0.5*atan2(2 Sxy, Sxx-Syy); a fiber
 *             runs perpendicular, theta = phi + pi/2 (mod pi). Picks whichever
 *             ONE family locally dominates -> the two cross-hatch families
 *             compete and the field speckles.
 *   rosy = 4  CROSS orientation (mod 90 deg): square the complex trace-free
 *             Hessian eigendirection and spatially average the resulting
 *             fourth-order representation E*(cos(4 theta),sin(4 theta)).
 *             Swapping two orthogonal ridge families leaves that value
 *             unchanged, yielding one unordered cross {+/-d,+/-perp(d)}
 *             rather than pretending to recover a directed fiber. theta is
 *             in [0,pi/2), with 0 = axes aligned.
 *
 * Method: (1) Gaussian pre-smooth I at grad_sigma; central-difference gradient
 * (gx,gy), x = column = u (horizontal), y = row = v (vertical). (2) accumulate
 * the second gradient moment for rosy=2 or squared trace-free Hessian moment
 * for rosy=4.
 * (3) Gaussian-smooth that moment at tensor_sigma (integration scale).
 * (4) extract the line/cross orientation and coherence.
 *
 * Convention: theta = 0 -> fiber/grid runs along +u (horizontal); pi/2 ->
 * along +v (vertical). theta is an axis, not a vector.
 * ==========================================================================*/

#define FIBER_NBINS 36   /* orientation histogram bins over [0, range_deg) */

typedef struct FiberField {
    float   *theta;   /* [W*H] orientation, radians in [0, range); 0 where !valid */
    float   *coh;     /* [W*H] coherence in [0,1]; 0 where !valid */
    uint8_t *valid;   /* [W*H] 1 = texel had coverage with full tensor support */
    int W, H;
    int rosy;                    /* 2 or 4 */
    double range_deg;            /* orientation domain: 180 (rosy2) or 90 (rosy4) */

    /* aggregate diagnostics (over valid texels; coherence-weighted where noted) */
    double  mean_coh;            /* mean coherence over valid texels */
    double  frac_valid;          /* valid / (W*H) */
    double  hist[FIBER_NBINS];   /* coherence-weighted orientation histogram, sums to 1 */
    double  peak_deg;            /* dominant histogram peak orientation (deg, [0,range)) */
    double  axis_frac;           /* coh-weighted fraction within axis_tol of an axis */
    double  mean_axis_err_deg;   /* coh-weighted mean |angle to nearest axis| (deg, [0,45]) */
} FiberField;

/* Compute the field from a row-major uint8 image [W*H]. Texels with img==0 are
 * background; the valid mask also drops texels whose tensor-scale neighbourhood
 * is not fully covered, so edges/hole rims don't masquerade as fibers.
 * grad_sigma  <= 0 -> 1.0 ; tensor_sigma <= 0 -> 4.0 ; axis_tol_deg <= 0 -> 15.0.
 * rosy != 4 is treated as 2. Arena-allocates out->theta/coh/valid.
 * Returns 0 on success, -1 on bad args. */
int FiberField_compute(Arena_T arena, const uint8_t *img, int W, int H,
                       double grad_sigma, double tensor_sigma, double axis_tol_deg,
                       int rosy, FiberField *out);

/* In-process unit test (synthetic oriented gratings + a cross-hatch). 0 = pass. */
int FiberField_selftest(void);

#endif
