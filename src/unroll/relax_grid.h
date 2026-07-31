#ifndef RELAX_GRID_INCLUDED
#define RELAX_GRID_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "../flatten/ribbon_relax.h"
#include "piece_set.h"

/* ============================================================================
 * relax_grid.h -- banded whole-grid UV relax (step4 back half). A one-shot
 * RibbonRelax over 30M verts is a 1-1.7h serial Gauss-Seidel; instead the
 * strip is tiled into u-bands relaxed in parallel RED-BLACK waves (even bands
 * then odd), each band pinning its halo verts (the new RibbonRelaxOpts.pin)
 * so concurrent writers never touch a shared vertex. Runs AFTER UVWeld_run:
 * the welded seams couple bands' interiors to their (pinned) halos, and the
 * alternating waves let seam strain diffuse across band boundaries.
 *
 * Each band works in a BAND-LOCAL uv frame (shifted by the band window
 * origin): that one shift both indexes the FiberField texels correctly
 * (RibbonRelax reads texel = uv/du absolute) and rescues float32 precision
 * at |u| ~ 1.2e6. The fiber field itself comes from the just-baked step3
 * rawtex, streamed per band via TiffIO_rows_*.
 * ==========================================================================*/

typedef struct RelaxGridOpts {
    RibbonRelaxOpts base;   /* sweeps default 16 here (E-curve gains <= ~12) */
    size_t band_px;         /* band width, strip px (def 32768) */
    size_t halo_px;         /* pinned halo, px (def 512) */
    int    rounds;          /* red-black rounds (def 1 = even wave + odd wave) */
    int    threads;         /* concurrent bands (def 16; fiber scratch cap) */
    double fib_grad_sigma, fib_tensor_sigma, fib_axis_tol;  /* 1.5/6/15 */
    int    verbose;
} RelaxGridOpts;

void RelaxGridOpts_default(RelaxGridOpts *o);

typedef struct RelaxGridStats {
    size_t n_bands, bands_run, bands_reverted, bands_no_fiber;
    size_t n_moved;
    double energy_before, energy_after;        /* summed over run bands */
    double stretch_before, stretch_after;      /* band-mean averages */
    double mean_disp, max_disp;
    double seconds;
} RelaxGridStats;

/* Mutates ps->uv in place. rawtex_tif = step3 bake (NULL or unreadable ->
 * isometry-only relax, counted in bands_no_fiber). u0/v0/du/dv = that bake's
 * canvas mapping (RasterStats). Returns 0, -1 on hard failure. */
int RelaxGrid_run(Arena_T arena, PieceSet *ps, const char *rawtex_tif,
                  double u0, double v0, double du, double dv,
                  const RelaxGridOpts *opts, RelaxGridStats *out);

int RelaxGrid_selftest(void);

#endif
