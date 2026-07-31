#ifndef SNAP_SOLVE_INCLUDED
#define SNAP_SOLVE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "surface_snap.h"

/* ============================================================================
 * snap_solve.h -- SOLVER stage: displace FIXABLE-ISLAND vertices onto the bright
 * papyrus target the detector found (target = vert + voff*vdir), keeping the
 * moved patch smooth and continuous with the surrounding GOOD surface, and NEVER
 * crossing into another wrap.
 *
 * One data-weighted Laplacian solve over the whole mesh (`SnapCG_solve`): FIXABLE-
 * island verts are pulled toward their target; every other vertex is PINNED in
 * place (the Dirichlet anchor that gives seam continuity). After the solve, each
 * moved vertex is re-checked against the ORIGINAL mesh's occupancy KD-tree and
 * REVERTED if it entered another sheet. Content is already baked into the target
 * by detection, so the solve does no RAW sampling.
 * ==========================================================================*/

typedef struct SnapSolveOpts {
    double alpha;        /* pull blend for fixable verts: w = alpha/(1-alpha)*deg
                          * (RawSnap convention); def 0.85: land close to the
                          * measured bright target while retaining 1-ring
                          * continuity */
    double w_anchor;     /* data weight pinning every non-fixable vertex (def 1e4) */
    int    max_iter;     /* CG iterations (def 96) */
    double tol;          /* CG tolerance (def 1e-4) */
    double occ_thresh;   /* occupancy radius for the revert guard, vox (def 0.9) */
    double local_r;      /* own-patch exclusion radius, vox (def 2.0) */
    double max_disp;     /* hard cap on per-vertex displacement, vox (def 8.0) */
    int    verbose;
} SnapSolveOpts;

void SnapSolveOpts_default(SnapSolveOpts *o);

typedef struct SnapSolveResult {
    size_t n_moved;         /* fixable-island verts displaced */
    size_t n_reverted;      /* moves rolled back by the occupancy guard */
    double mean_disp;       /* mean displacement of moved verts (vox) */
    double max_disp;        /* max displacement kept (vox) */
} SnapSolveResult;

/* Displace `verts` in place. faces/nf define the 1-ring; det = detection result
 * (reads vclass, vregion, rclass, voff, vdir). Returns 0 on success, -1 on
 * degenerate input. */
int SnapSolve_run(Arena_T arena, float *verts, size_t nv,
                  const int32_t *faces, size_t nf,
                  const SnapResult *det, const SnapSolveOpts *opts,
                  SnapSolveResult *out);

/* In-process unit tests. Returns 0 if all pass, else #failures. */
int SnapSolve_selftest(void);

#endif
