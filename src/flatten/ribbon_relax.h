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
    int    sweeps;         /* Gauss-Seidel sweeps (def 40; converges by ~30) */
    int    line_search;    /* backtracking steps per vertex (def 10) */
    double max_disp;       /* cap on per-vertex |uv - uv_in| in vox (def 40) */
    double qc_reject;      /* drop faces whose INITIAL stretch sigma1/sigma2 exceeds this
                            * (pre-broken UV-collapsed slivers + overlap-relocated faces);
                            * they otherwise dominate the energy. <=0 disables. (def 50) */
    int    fix_boundary;   /* 1 = pin boundary-loop UVs; 0 = free (lets the sheet globally
                            * rotate/shear u onto the fibers + flatten) (def 0) */
    const uint8_t *pin;    /* [nv] nullable: nonzero = vertex immovable (the
                            * whole-grid band tiling pins halo verts with
                            * this). def NULL = none */
    int    verbose;        /* per-sweep energy logging */
} RibbonRelaxOpts;

/* Fill with defaults. */
void RibbonRelax_defaults(RibbonRelaxOpts *o);

typedef struct RibbonRelaxStats {
    double stretch_mean_before, stretch_mean_after;  /* area-wtd mean quasi-conformal sigma1/sigma2 */
    double stretch_max_before,  stretch_max_after;
    double axis_err_before, axis_err_after;          /* coh-wtd mean fiber axis error, deg */
    double energy_before, energy_after;              /* total E (area-wtd) */
    double mean_disp, max_disp;                      /* per-movable-vertex |uv - uv_in| in vox */
    int    flips_before, flips_after;                /* UV signed-area sign flips vs first face */
    size_t n_interior, n_moved, n_fiber_faces, n_reject;
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

/* In-process unit test. 0 = pass. */
int RibbonRelax_selftest(void);

#endif
