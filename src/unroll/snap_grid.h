#ifndef SNAP_GRID_INCLUDED
#define SNAP_GRID_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "../common/raw_sample.h"
#include "piece_set.h"

/* ============================================================================
 * snap_grid.h -- whole-grid two-pass snap (step3_snap): dark-whiff quilting
 * repair followed by one iterative global recto placement, reformulated for
 * the 30M-vert piece set.
 *
 * Scale decomposition -- what stays GLOBAL and what factorizes:
 *   GLOBAL: the occupancy KD-tree over ALL verts (a band-local tree would
 *   miss another wrap at the same 3D spot -- probes would sail through other
 *   sheets, violating the no-merger north star), the contrast window, and
 *   the per-vertex cv sampling.
 *   PER-CUBE (exact only pre-weld): an unwelded piece set shares no vertices
 *   across cubes, so the mesh 1-ring -- and with it the GCO dark cut, quilting
 *   target MRF, dark-region grouping, and pinned Laplacian solve -- factorizes into
 *   1702 independent per-cube problems, run in parallel. AFTER the
 *   step2_join weld cross-cube edges exist: set global_mode (one whole-mesh
 *   problem, exact at 4x5x5 scale). Per-cube mode guards itself by skipping
 *   faces that reference verts outside their cube's range.
 *
 * Mutates ps->verts in place (uv untouched); caller refreshes normals and
 * re-bakes. On-disk placed records never change.
 * ==========================================================================*/

typedef struct SnapGridOpts {
    /* detection -- defaults mirror SnapOpts_default (surface_snap.h) */
    float  axis_point[3], axis_dir[3];
    double reach, occ_thresh, local_r, gain_bias, step;
    double tensor_weight, tensor_radius; /* legacy, ignored by pass 1 */
    double band_frac, min_gain, anchor_frac;
    int    region_cap, min_region;
    double pct_lo, pct_hi;
    int    dark_thresh;
    double dark_lambda, dark_sigma;
    double normal_range;
    int    normal_samples;
    /* pass-1 coherent repair target */
    int    target_mode, depth_bins;
    double w_match, w_close, smooth_mu, smooth_tau;
    /* solve -- defaults mirror SnapSolveOpts_default (snap_solve.h) */
    double alpha, w_anchor;
    int    cg_max_iter;
    double cg_tol, max_disp;
    /* pass-2 oriented recto refinement; recto_iters=0 disables it */
    int    recto_iters, recto_inner_iters;
    double recto_range, recto_max_total;
    /* scale */
    int    threads;          /* OpenMP workers over cubes (def 32) */
    int    global_mode;      /* treat the whole mesh as ONE cube (exact when
                                cross-cube edges exist, i.e. after the
                                step2_join weld; falls back to per-cube with
                                a face-range guard above SG_GLOBAL_NV_CAP) */
    int    verbose;
} SnapGridOpts;

void SnapGridOpts_default(SnapGridOpts *o);

typedef struct SnapGridStats {
    size_t n_sampled, n_dark, n_fixable, n_crack, n_rejected;
    size_t nreg, n_reg_fixable, n_reg_crack, n_reg_anchorless, n_reg_mixed,
           n_reg_reject;
    size_t max_region_size;
    size_t n_moved, n_reverted;
    double mean_disp, max_disp;
    double mean_target_tensor; /* legacy diagnostic, always zero */
    double mean_quilt_cost, mean_target_dist;
    double win_lo, win_hi;
    size_t n_cubes_run, n_cubes_gco_fallback;
    size_t n_quilt_fallback;
    size_t n_recto_supported, n_recto_moved, n_recto_reverted,
           n_recto_slope_limited;
    int recto_iterations;
    double recto_mean_disp, recto_max_disp, recto_mean_gradient;
    double sec_occ, sec_cv, sec_cubes, sec_recto;
} SnapGridStats;

/* ct MUST be pre-warmed (cubetable_prewarm_all) -- workers read it
 * concurrently. Returns 0 (empty ok), -1 on hard failure. */
int SnapGrid_run(Arena_T arena, PieceSet *ps, CubeTable *ct,
                 const SnapGridOpts *opts, SnapGridStats *out);

int SnapGrid_selftest(void);

#endif
