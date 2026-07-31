#ifndef ARC_REG_INCLUDED
#define ARC_REG_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "placed_cube.h"
#include "weld_track.h"

/* ============================================================================
 * arc_reg.h -- ArcReg: the L2 quilting fix. Generalizes group_graph's per-
 * (cube,gid) scalar `du` into a smooth piecewise-linear u-WARP w(phi),
 * constrained by the WTRK cross-cube correspondences (individual same-papyrus
 * springs, conf-weighted). This closes the arc-length continuity gap that a
 * single constant du per group leaves (the per-cube tiling's stretch floor
 * 1.35 vs the welded 1.115): where two cubes' copies of one wrap drift in
 * arc length across a seam, a constant shift can't correct a phi-DEPENDENT
 * mismatch, but a warp can.
 *
 * The winding integers k (group_graph) are UNTOUCHED -- winding is solved.
 * ArcReg only refines du -> du + delta(phi), delta a small zero-ish shape.
 *
 * Applied at PlacedCube_finalize as
 *     u_reg = u_raw + DeltaU(phi_raw, k) + du + delta_g(phi_raw)
 * so it composes with the existing per-group registration and everything
 * downstream (bake, ownership, atlas) is unchanged.
 *
 * A correspondence c (a in node (cube_a,gid_a), b in (cube_b,gid_b)) asserts
 * the two boundary verts are the same physical papyrus, so their registered u
 * should match. With r0 = u_reg_a - u_reg_b (the current mismatch),
 *   minimize  sum_c conf_c*(r0_c + d_a(phi_a) - d_b(phi_b))^2
 *           + lambda_smooth * sum_node sum_knots (d[i+1]-d[i])^2
 *           + lambda_ridge  * sum d^2          (removes the global-u gauge)
 * solved matrix-free CG on the normal equations; d is magnitude-bounded.
 * ==========================================================================*/

typedef struct {
    int    n_knots;         /* piecewise-linear warp knots per node (def 5) */
    double lambda_smooth;   /* inter-knot smoothness weight (def 4.0) */
    double lambda_ridge;    /* ridge toward delta=0; pins the gauge (def 0.02) */
    int    cg_iters;        /* CG iterations on the normal equations (def 200) */
    int    relax_iters;     /* broad StrokeStrip likelihood rounds (def 5) */
    int    final_iters;     /* tighter final likelihood rounds (def 2) */
    double sigma_floor;     /* residual-likelihood sigma floor, vox (def 2) */
    double max_warp;        /* clamp |delta| per knot, vox (def 40) */
    int    min_node_corr;   /* nodes with fewer springs get a scalar warp
                             * (n_knots forced to 1 -> a constant du refine)
                             * (def 8) */
    int    conf_min;        /* use only corrs with conf >= this (def 1) */
    double max_mismatch;    /* drop springs whose current |u_a - u_b| exceeds
                             * this (vox): a whole-circumference mismatch is a
                             * winding turn-off (two cubes disagree on the
                             * wrap's k -- group_graph's job), NOT arc-length
                             * drift; keeping it saturates the warp clamp.
                             * (def 60 = a few pitch, << a mid-wrap circumference) */
    int    verbose;
} ArcRegOpts;

void ArcRegOpts_default(ArcRegOpts *o);

/* One node's warp: piecewise-linear over [phi0, phi0 + (n_knots-1)*dphi],
 * clamped-constant outside. delta[] are the knot heights (vox). */
typedef struct {
    int32_t cube, gid;
    int     n_knots;
    double  phi0, dphi;     /* knot 0 phi, knot spacing (rad) */
    double  delta[8];       /* <= ARC_REG_MAX_KNOTS heights */
} ArcRegWarp;

#define ARC_REG_MAX_KNOTS 8

typedef struct {
    size_t n_corr_used;     /* springs after the conf + mismatch gate */
    size_t n_gated;         /* springs dropped as winding turn-off (|r0|>max) */
    size_t n_downweighted;  /* final likelihood < 0.1 (WTS outlier candidates) */
    size_t n_nodes;         /* distinct (cube,gid) nodes with springs */
    size_t n_scalar_nodes;  /* nodes forced to a scalar (1-knot) warp */
    double med_before;      /* MEDIAN corr u-mismatch before the warp, vox
                             * (robust; RMS is dominated by the turn-off tail) */
    double med_after;       /* median after (the quilting improvement) */
    double rms_before, rms_after;   /* over the kept (gated) springs */
    double p95_before, p95_after;
    size_t monotone_clamps; /* legacy name: knots clipped to max_warp */
    double max_abs_warp;    /* largest |delta| applied, vox */
    int    reweight_rounds; /* completed likelihood/solve alternations */
    double seconds;
} ArcRegStats;

/* Solve the warp over all cubes. `wt_to_local[wt->n_cubes]` maps a WTRK cube
 * index to the caller's cube index (skins/nskin/regs order); -1 = absent.
 * `skins[c]/nskin[c]` are the RAW boundary skins (SkinVert, indexed like the
 * WTRK corr's local vert index); `regs[c]` the group_graph registration.
 * On return, `out_warps` (nullable) receives up to *out_n_warps ArcRegWarp
 * records (caller sizes the buffer; pass NULL to only measure). Returns 0. */
int ArcReg_solve(Arena_T arena,
                 const WeldTrack *wt, const int32_t *wt_to_local,
                 SkinVert *const *skins, const size_t *nskin,
                 const PlacedReg *regs, size_t n_cubes,
                 double spiral_a, double spiral_b,
                 const ArcRegOpts *opts,
                 ArcRegWarp *out_warps, size_t *out_n_warps,
                 ArcRegStats *stats);

/* Evaluate a node's warp at raw phi (clamped-constant outside the knot span).
 * NULL / 0-knot warp -> 0. Used by finalize to apply the warp. */
double ArcReg_eval(const ArcRegWarp *w, double phi);

int ArcReg_selftest(void);

#endif
