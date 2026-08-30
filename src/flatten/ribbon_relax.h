#ifndef RIBBON_RELAX_INCLUDED
#define RIBBON_RELAX_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "fiber_field.h"

/* ============================================================================
 * ribbon_relax.h -- step 4 of the ribbon pipeline: relax the UV parameterization
 * to minimize parametric distortion (keep it a flat sheet) while aligning the
 * papyrus fibers to the texture axes -- WITHOUT introducing foldovers or
 * self-intersections. Runs after snap; verts are frozen, only UV moves.
 *
 * Per-face energy over the 2x2 Jacobian J (rest isometric frame -> UV):
 *   E_iso(J)   = ||J||_F^2 + ||J^-1||_F^2   (symmetric Dirichlet; = 4 at J=rotation,
 *                -> +inf as det J -> 0, so it is also the foldover barrier)
 *   E_align(J) = coh * sin^2(2*phi),  phi = angle of J * d_fiber in UV
 *                (0 when the fiber direction lands on a u/v axis, mod 90)
 *   E_face     = A_rest * (E_iso + lambda * E_align)
 * The fiber target d_fiber (material direction) and its coherence come from a
 * FiberField (rosy=4 grid orientation) aggregated over each face's texel
 * footprint; low-coherence faces (cracks/holes) get no alignment term.
 *
 * Optimizer: Gauss-Seidel over interior vertices (boundary pinned). Each vertex
 * takes a backtracking line search along -grad of its incident-face energy;
 * any step that would flip an incident UV triangle (sign of signed area) or
 * exceed max_disp is rejected -> the embedding is preserved by construction.
 * Fail-closed: if the result raises total distortion or the flip count, the
 * input UV is returned unchanged.
 * ==========================================================================*/

typedef struct RibbonRelaxOpts {
    double lambda_align;   /* weight of the fiber-alignment term; balanced operating point (def 0.5) */
    double coh_gate;       /* min aggregated coherence to apply alignment (def 0.15) */
    int    sweeps;         /* Gauss-Seidel cap (default 1000) */
    double movement_tolerance; /* stop when the largest accepted UV move in a
                                  sweep is <= this (default 1e-5) */
    int    line_search;    /* backtracking steps per vertex (def 24; enough to
                              resolve the default movement tolerance) */
    double max_disp;       /* cap on per-vertex |uv - uv_in| in vox; <=0 disables
                            * the cap (the injectivity guards still apply) (def 40) */
    double qc_reject;      /* drop faces whose INITIAL stretch sigma1/sigma2 exceeds this
                            * (pre-broken UV-collapsed slivers + overlap-relocated faces);
                            * they otherwise dominate the energy. <=0 disables. (def 50) */
    int    reference_metric; /* 1 = use reference_uv's planar metric as the
                            * rest metric instead of frozen XYZ (default 0).
                            * Local robust repair uses this when the detector
                            * says XYZ itself is not a trustworthy metric
                            * observation inside the selected patch. */
    int    fix_boundary;   /* 1 = pin boundary-loop UVs; 0 = free (lets the sheet globally
                            * rotate/shear u onto the fibers + flatten) (def 0) */
    int    convex_boundary;/* when boundary is free, preserve every boundary
                            * turn of an initially convex disk (def 0) */
    int    global_injective;/* when boundary is free, preserve its ambient
                             * isotopy: no non-adjacent boundary edges may
                             * touch or cross during a vertex move (def 0) */
    const double *guide_uv;/* optional [nv*2] globally aligned winding guide */
    double guide_weight;   /* area-weighted soft guide strength (def 0) */
    const double *guide_weight_vertex; /* optional [nv] multiplier for the
                            * guide term. Local repair uses a smooth 0..1
                            * confidence ramp instead of imprinting a binary
                            * detector boundary into UV. */
    const double *reference_uv; /* optional [nv*2] fixed material/raster
                            * coordinates for nested coarse-to-fine starts.
                            * Fiber lookup, material pullback, QC, and max_disp
                            * use this reference; uv_in remains the numerical
                            * starting iterate. NULL means uv_in. */
    const uint8_t *pin;    /* [nv] nullable: nonzero = vertex immovable (the
                            * whole-grid band tiling pins halo verts with
                            * this). def NULL = none */
    int    verbose;        /* per-sweep energy logging */
} RibbonRelaxOpts;

/* Fill with defaults. */
void RibbonRelax_defaults(RibbonRelaxOpts *o);

typedef struct RibbonRelaxStats {
    int    iterations, converged;
    double last_max_movement;
    double stretch_mean_before, stretch_mean_after;  /* area-wtd mean quasi-conformal sigma1/sigma2 */
    double stretch_max_before,  stretch_max_after;
    double axis_err_before, axis_err_after;          /* coh-wtd mean fiber axis error, deg */
    double energy_before, energy_after;              /* total E (area-wtd) */
    double mean_disp, max_disp;                      /* per-movable-vertex |uv - uv_in| in vox */
    double guide_rms_before, guide_rms_after;
    int    flips_before, flips_after;                /* UV signed-area sign flips vs first face */
    size_t n_interior, n_moved, n_fiber_faces, n_reject;
    size_t n_convex_reject;
    size_t n_boundary_edges, n_boundary_collision_reject;
    size_t boundary_intersections_before, boundary_intersections_after;
    int    reverted;                                 /* 1 = guard tripped, uv_out == uv_in */
} RibbonRelaxStats;

/* verts[nv*3] (z,y,x, frozen), faces[nf*3], uv_in[nv*2]. `fib` is a rosy=4
 * FiberField computed on the current baked rawtex; du/dv are that raster's
 * pixel pitch in vox (texel = uv/du). face_skip (nullable) omits faces from the
 * energy entirely (e.g. the overlap-relocated/dropped faces that are not baked),
 * so they never tug the sheet; vertices left with no kept face are pinned.
 * Writes uv_out[nv*2] (may be uv_in on a revert). opts NULL -> defaults.
 * Returns 0 on success, -1 on bad args. */
int RibbonRelax_run(Arena_T arena,
                    const float *verts, size_t nv,
                    const int32_t *faces, size_t nf,
                    const float *uv_in, const uint8_t *face_skip,
                    const FiberField *fib, double du, double dv,
                    const RibbonRelaxOpts *opts,
                    float *uv_out, RibbonRelaxStats *stats);

/* Double-precision UV entry point. Large packed disk atlases can have
 * sub-voxel triangles tens of thousands of voxels from the origin; converting
 * those coordinates to float before the barrier solve can collapse an
 * otherwise valid face. The optimizer and all guards are identical to
 * RibbonRelax_run, but UV stays double from input through output. */
int RibbonRelax_run_double(Arena_T arena,
                           const float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const double *uv_in, const uint8_t *face_skip,
                           const FiberField *fib, double du, double dv,
                           const RibbonRelaxOpts *opts,
                           double *uv_out, RibbonRelaxStats *stats);

/* In-process unit test. 0 = pass. */
int RibbonRelax_selftest(void);

#endif
