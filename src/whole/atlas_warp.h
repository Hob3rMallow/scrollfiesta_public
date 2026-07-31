#ifndef ATLAS_WARP_INCLUDED
#define ATLAS_WARP_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_radial_order.h"
#include "atlas_strip.h"
#include "monotone_qp.h"

/* ============================================================================
 * atlas_warp.h -- re-solve the coarse atlas under radial ordering.
 *
 * The predecessor stages moved whole charts by a constant in u.  That cannot
 * express what a scroll actually needs: one turn of correction shifts u by
 * 2*pi*a*k + b*k*phi + pi*b*k^2, which varies along the wrap, so collapsing it
 * to a chart mean converts a shear into a teleport (measured: 21k-vox gauge
 * shifts, cross-cube weld stretch p95 3.2 -> 1663).
 *
 * Here every coarse sample keeps its own unknown and the StrokeStrip rows that
 * produced the HEAD atlas are re-solved verbatim, with two additions:
 *
 *   proximity   u_i ~= base_i,  weight lambda_prior
 *   ordering    the hard bounds and soft spacing rows from atlas_radial_order
 *
 * The proximity prior is what keeps the good 95% of the atlas still.  It is
 * also what lets the ANCHORS GO.  AtlasStrip pins one exact gauge per support
 * component; a ray bound between two anchored components is then unsatisfiable
 * and would be silently dropped as infeasible -- discarding exactly the
 * cross-component evidence the stage exists to apply.  A prior on every sample
 * makes the Hessian positive definite on its own, so gauge is determined
 * without pinning anything, and a component free-floating by a turn can be
 * carried to where the ordering evidence says it belongs.
 *
 * Division of labour between the two ordering strengths:
 *
 *   The hard bound forbids local inversion -- the outer wrap dipping below the
 *   inner one at the same angle, which is what shreds a region into interleaved
 *   overlap.  It is deliberately small, so a misgated pair costs voxels.
 *
 *   Wrap-to-wrap SEPARATION (the outer wrap's u interval starting after the
 *   inner wrap's ends) is carried by the soft spacing target plus the strip's
 *   own unit-speed length rows.  Where a sheet is connected through the turn,
 *   arc length already supplies it for free and the spacing row agrees; where a
 *   piece floats free, the spacing row is the only evidence there is, so it
 *   must outweigh the proximity prior -- hence lambda_spacing >> lambda_prior.
 * ==========================================================================*/

/*
 * COST NOTE, measured on the 4x5x5 fixture.  The hard ordering bounds are what
 * force an active-set method, and MonotoneQp admits ONE bound per iteration
 * with a full refactorization each time.  With 104k ray bounds of which half
 * start violated -- because half the wrap pairs in the incoming atlas are
 * coincident -- that is tens of thousands of factorizations of a 259k-variable
 * system.  It does not finish.
 *
 * The bounds are also the least valuable constraint: the warp fixture shows an
 * 8-vox ordering bound cannot even detect a whole-turn error, which is the
 * error that actually matters.  Separation is carried entirely by the spacing
 * rows.  So hard bounds are OFF by default and the warp is a plain weighted
 * least-squares solve -- one factorization per IRLS round.  Turn them on only
 * to clean up residual local inversions, ideally on a lazily-generated subset.
 */
typedef struct {
    double lambda_prior;    /* proximity to the incoming field (0.05) */
    double lambda_spacing;  /* scales the radial spacing rows (0.25) */
    int    irls_rounds;     /* robust reweighting of the spacing evidence (4) */
    double irls_scale;      /* tolerated error as a fraction of one turn (0.35) */
    int    hard_bounds;     /* 1 enables the ordering inequalities (default 0) */
    int    keep_anchors;    /* 1 restores AtlasStrip's exact gauges (default 0) */
} AtlasWarpOptions;

void AtlasWarpOptions_default(AtlasWarpOptions *opts);

typedef struct {
    size_t nvar;
    size_t nrows;
    size_t nbounds;
    size_t ray_bounds;
    size_t ray_bounds_violated_before; /* disorder present in the input */
    size_t ray_bounds_disabled;        /* dropped by the bootstrap */
    size_t spacing_rows;

    AtlasRadialOrderBootstrapStats bootstrap;

    /* 0 once the QP has been reached; negative names the setup step that
     * failed, so an early return is still diagnosable from the stats. */
    int    stage_rc;
    int    qp_rc;
    int    qp_iterations;
    int    qp_active_final;
    double objective_initial;
    double objective_final;
    double min_slack;
    double stationarity;

    /* Displacement from the incoming field, the teleportation detector.
     * Reported relative to gauge_shift, since a constant added to every u is
     * not a displacement at all -- it is the same atlas. */
    double gauge_shift;
    double shift_rms;
    double shift_p95;
    double shift_max;

    size_t ray_bounds_violated_after;

    /*
     * Spacing residual: how far the field is from putting radial neighbours
     * one circumference apart.  This, not the bound-violation count, is the
     * measure of wrap-level disorder -- a component a whole turn out of place
     * still satisfies an 8-vox ordering bound comfortably, and only the
     * spacing evidence can see it.
     */
    double spacing_rms_before;
    double spacing_rms_after;
    double spacing_p95_before;
    double spacing_p95_after;

    /* Robust reweighting: how much of the radial evidence survived, and how
     * much of it the solve had to reject as contradictory. */
    int    irls_rounds_run;
    double irls_sigma;
    size_t spacing_rows_downweighted;
    double spacing_weight_retained;

    double span_before;
    double span_after;
} AtlasWarpStats;

/*
 * problem      the HEAD candidate graph with frozen memberships (FINAL mode)
 * order        ray constraints built against the same sample indexing
 * base_sample  incoming coarse field, [problem->nsamples]
 * out_sample   solved coarse field, [problem->nsamples]; may alias base_sample
 */
int AtlasWarp_solve(Arena_T arena,
                    const AtlasStripProblem *problem,
                    const AtlasStripOptions *strip_opts,
                    const MonotoneQpOptions *qp_opts,
                    const AtlasRadialOrderSet *order,
                    const double *base_sample,
                    const AtlasWarpOptions *opts,
                    double *out_sample,
                    AtlasWarpStats *stats);

/* ==========================================================================
 * Per-component gauge solve.
 *
 * Cutting the wrap-fusing links does not, by itself, push anything apart -- it
 * only converts wraps that were welded into one rigid body into separate
 * bodies, each with one free constant.  THIS is what then places them: the
 * radial evidence says component A and component B must differ by one local
 * circumference, and with the gauges finally free that is a solvable system.
 *
 * It is also the cheap form.  The full per-sample warp adds 100k long-range
 * rows and wrecks the elimination tree (a single factorization ran past ten
 * minutes).  Here the unknowns are one shift per support component -- 558 on
 * the 4x5x5 fixture after pruning, against ~100k pair equations -- so the
 * normal equations are a few hundred square and solve instantly.
 *
 * Pairs INSIDE one component are skipped: they carry no information about a
 * variable that exists.  If most pairs are intra-component, that is the signal
 * that the graph still needs cutting, and it is reported.
 * ========================================================================== */
typedef struct {
    size_t components;
    size_t equations;          /* cross-component pairs used */
    size_t intra_skipped;      /* pairs with nothing to say -- graph not cut */
    int    irls_rounds_run;
    size_t downweighted;
    double residual_rms_before;
    double residual_rms_after;
    double shift_rms;
    double shift_max;
    double gauge_shift;
    size_t not_converged;      /* IRLS rounds that hit the sweep cap */
    double final_delta;
    size_t equation_islands;   /* disconnected pieces of the equation graph */
    int    solve_failed;
} AtlasWarpGaugeStats;

/*
 * sample_component[i]  support component of sample i
 * base_sample[i]       incoming coarse field
 * out_shift[c]         solved additive shift per component
 */
int AtlasWarp_solve_gauges(Arena_T arena,
                           const AtlasRadialOrderSet *order,
                           const int32_t *sample_component,
                           size_t nsamples,
                           size_t ncomponents,
                           const double *base_sample,
                           const AtlasWarpOptions *opts,
                           double *out_shift,
                           AtlasWarpGaugeStats *stats);

int AtlasWarp_selftest(void);

enum { ATLAS_WARP_ROW_PRIOR = 12 };

#endif
