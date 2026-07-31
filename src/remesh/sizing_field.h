#ifndef SIZING_FIELD_INCLUDED
#define SIZING_FIELD_INCLUDED

/*
 * CvtField — a spatially-varying target-edge-length (sizing) field for the CVT
 * remesher, used to grade density: fine (weldable) near the cube's owned-box
 * planes, coarse in the interior.
 *
 * The field is a smoothstep ramp in DEPTH into the owned box: at a plane the
 * target edge length is h_seam (derived from the seam-weld bridge cap so seams
 * always close); band_lo voxels in it is still h_seam (flat pad); by band_hi it
 * has reached h_interior and stays there. Halo points (depth < 0) clamp to the
 * densest end — correct, they hold the near-plane verts the trim keeps.
 *
 * Density enters the CVT ONLY as a scalar weight on triangle area inside
 * Rvd_accumulate; it touches no geometric predicate, so a NULL field leaves the
 * exact clip/dual byte-identical. Seeding uses point density 1/h^2 (spacing ∝ h);
 * the CWF energy uses 1/h^energy_exp with energy_exp = 4 so a restricted 2-D CVT's
 * fixed point (spacing ∝ rho_E^-1/4) equals the seed — the gradient is preserved,
 * not relaxed away, over the short iteration budget.
 */

#include <math.h>

typedef struct {
    double lo, hi;          /* owned-box planes (scalar, all axes), local frame  */
    double h_seam;          /* target edge length at the owned-box planes (vox)  */
    double h_interior;      /* target edge length deep in the interior (vox)     */
    double band_lo;         /* flat full-density pad depth from each plane (vox)  */
    double band_hi;         /* depth at which h has reached h_interior (vox)      */
    int    energy_exp;      /* CWF area-weight exponent (4 stable; 2 legacy)      */
} CvtField;

/* Target edge length h(x) at a point (local z,y,x). */
static inline double sizing_h(const CvtField *f, double c0, double c1, double c2) {
    double d0 = c0 - f->lo, d1 = f->hi - c0;
    double d2 = c1 - f->lo, d3 = f->hi - c1;
    double d4 = c2 - f->lo, d5 = f->hi - c2;
    double d = d0;
    if (d1 < d) d = d1;
    if (d2 < d) d = d2;
    if (d3 < d) d = d3;
    if (d4 < d) d = d4;
    if (d5 < d) d = d5;                       /* signed depth into the owned box */
    double w = f->band_hi - f->band_lo;
    double u = (w > 1e-9) ? (d - f->band_lo) / w : (d >= f->band_lo ? 1.0 : 0.0);
    if (u < 0.0) u = 0.0; else if (u > 1.0) u = 1.0;
    double t = u * u * (3.0 - 2.0 * u);        /* smoothstep (C1: zero end slopes) */
    return f->h_seam + (f->h_interior - f->h_seam) * t;
}

/* CWF energy density weight ρ_E(x) = 1 / h(x)^energy_exp. */
static inline double Rvd_density_rho(const CvtField *f, double c0, double c1, double c2) {
    double h = sizing_h(f, c0, c1, c2);
    if (h < 1e-9) h = 1e-9;
    double h2 = h * h;
    return (f->energy_exp == 2) ? (1.0 / h2) : (1.0 / (h2 * h2));  /* default exp=4 */
}

/* Seeding point-density weight = 1 / h(x)^2 (spacing ∝ h). */
static inline double sizing_seed_rho(const CvtField *f, double c0, double c1, double c2) {
    double h = sizing_h(f, c0, c1, c2);
    if (h < 1e-9) h = 1e-9;
    return 1.0 / (h * h);
}

#endif /* SIZING_FIELD_INCLUDED */
