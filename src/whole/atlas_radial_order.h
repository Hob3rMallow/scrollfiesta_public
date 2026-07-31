#ifndef ATLAS_RADIAL_ORDER_INCLUDED
#define ATLAS_RADIAL_ORDER_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../unroll/scaffold.h"
#include "atlas_strip.h"
#include "cube_register.h"
#include "monotone_qp.h"

/* ============================================================================
 * atlas_radial_order.h -- the one fact that makes a scroll a scroll.
 *
 * Along a ray leaving the umbilicus, wrap k+1 is farther out in u than wrap k.
 * Nothing in the atlas solve ever asserted that, so overlap removal had to be
 * bought by translating charts until they stopped intersecting -- which tears
 * the physical weld (edge stretch p95 3.2 -> 1663 on the 4x5x5 fixture).
 *
 * The order is MEASURED, not optimized.  Coarse samples already sit on axial
 * slice planes, so a ray is just a narrow angular window on one plane.  Every
 * sample is joined to the nearest sample outside it in that window -- its
 * "next wrap out" -- and that single pair becomes a difference constraint on u.
 * Radius is used ONLY to rank and to count wraps; no registered quantity
 * enters, so the constraints cannot inherit the error they are meant to fix.
 *
 * The window is an ARC length, not a fixed angle: its half-width is
 * arc_window / r, so it spans the same few samples whether it sits near the
 * core or out at the rim.  Fixed angular bins were tried first and were wrong
 * for two reasons -- a ray straddling a bin edge silently loses its
 * constraints, and the sample spacing per bin varies by an order of magnitude
 * across the scroll.
 *
 * Because a pair exists only when the outer radius exceeds the inner one by a
 * real separation, the ray graph alone is ACYCLIC by construction.  A positive
 * cycle can therefore only appear once stroke monotonicity is added, and that
 * is precisely the signature of a component that physically fuses two wraps.
 *
 * TWO SIGNALS, TWO DOMAINS.  Radius FINDS the neighbour.  Inside one upstream
 * registration domain (normally one cube), phi CLASSIFIES it: the unknown
 * integer gauge is common-mode and cancels in a difference, so turn zero still
 * identifies a delamination.  Across domains, absolute phi is invalid -- its
 * integer part was rounded independently -- and using it creates exactly the
 * one-circumference seam error this module is supposed to repair.  There the
 * classifier is the branch-safe local geometric winding differential
 *
 *     dw = sense * (r_outer-r_inner) / pitch - wrap(dtheta) / (2*pi).
 *
 * The signed spacing target never uses absolute phi.  For m turns it is the
 * physical Archimedean-spiral arclength at the pair midpoint,
 *
 *     du = sense * m * 2*pi*r_mid.
 *
 * This is exact for a calibrated Archimedean spiral and cannot inherit an
 * upstream per-cube turn rounding.
 *
 * Two strengths per admitted pair:
 *
 *   hard   u[outer] - u[inner] >= gap_min          (ordering, ~8 vox)
 *   soft   u[outer] - u[inner] ~= m * 2*pi*r_mid   (spacing, m turns)
 *
 * The hard bound is deliberately tiny.  It states only "outer is outside",
 * which is what a crossing violates; the spacing target carries the magnitude
 * and yields to evidence.  A misgated pair therefore costs ~8 vox of local
 * distortion instead of a half-circumference tear.
 *
 * Delamination is exempt BY CONSTRUCTION.  Two plies of one physical wrap sit
 * at nearly the same radius, so their separation falls under the cross-wrap
 * threshold and no constraint is emitted -- the legitimate case never needs a
 * whitelist lookup.  The optional whitelist covers only eccentric geometry
 * where a classified delamination is radially wide enough to look cross-wrap.
 *
 * ==========================================================================*/

typedef struct {
    double arc_window;           /* angular half-window as an arc length (6.0) */
    double separation_min_ratio; /* |dr| >= this * pitch is cross-wrap (0.6) */
    double turn_tolerance;       /* |dphi/2pi - m| gate for admission (0.30) */
    double core_radius_pitches;  /* skip pairs inside this many pitches (2.0) */
    double gap_min;              /* hard ordering bound, vox (8.0) */
    double lambda_spacing;       /* soft spacing row weight (0.10) */
    int    max_turns;            /* drop pairs spanning more wraps (4) */
} AtlasRadialOrderOptions;

void AtlasRadialOrderOptions_default(AtlasRadialOrderOptions *opts);

/* One admitted radial-adjacency observation.  inner/outer are sample indices
 * ordered by 3D radius; the u-space direction is applied when the bound is
 * emitted, so bound.hi is whichever of the two carries the larger u. */
typedef struct {
    int32_t inner;
    int32_t outer;
    int32_t plane;
    int32_t turns;          /* wraps crossed, from phi; always >= 1 */
    double  radius_inner;
    double  radius_gap;
    double  phi_inner;       /* diagnostic only; never used for the target */
    double  spacing_target; /* SIGNED u target for (outer - inner) */
    double  weight;
} AtlasRadialOrderPair;

typedef struct {
    size_t samples_considered;
    size_t rejected_core;         /* inside the collapsed core */
    size_t rejected_no_neighbour; /* nothing outside it in the window */
    size_t rejected_delamination; /* whitelisted component pair */
    size_t rejected_same_wrap;    /* phi says turn 0: delamination, exempt */
    size_t rejected_ambiguous;    /* phi not near a whole turn */
    size_t rejected_wide;         /* more than max_turns out */
    size_t admitted;
    size_t planes_used;
    size_t window_visits;         /* neighbour scans, for cost accounting */

    /* Independent cross-check: how often the radial gap implies the same wrap
     * count phi does.  Low agreement means the axis wanders enough locally
     * that radius is only good for finding a neighbour, not for counting. */
    size_t turn_agree_radius;
    size_t turn_disagree_radius;

    /* Disorder present in the incoming field, in units of one circumference. */
    size_t coincident_pairs;      /* |du| < half the target: an overlap */
    size_t satisfied_pairs;       /* du within half a turn of the target */
    size_t reversed_pairs;        /* du of the wrong sign entirely */
    double median_radius_gap;
} AtlasRadialOrderStats;

typedef struct {
    AtlasRadialOrderPair *pairs;
    size_t npairs;

    MonotoneQpBound *bounds;   /* npairs entries; stroke = -1 marks a ray */
    size_t nbounds;

    /* Soft spacing rows in MonotoneQp form, ready to append to a strip system.
     * Each row is (x[hi] - x[lo]) ~= target, so coeff holds 2 entries. */
    MonotoneQpRow   *rows;
    size_t           nrows;
    MonotoneQpCoeff *coeff;
    size_t           ncoeff;

    AtlasRadialOrderStats stats;
} AtlasRadialOrderSet;

enum { ATLAS_RADIAL_ROW_SPACING = 11 };

/*
 * Build the ray constraints.
 *
 * sample_plane[i]     axial slice index of sample i (AtlasCandidateSampleRef.plane)
 * sample_component[i] retained-triangle component, or NULL to skip the whitelist
 * sample_phi[i]       multi-turn winding angle, radians; trusted only when
 *                     both samples have the same phase_domain
 * phase_domain[i]     common integer-gauge domain (normally cube id), or NULL
 *                     when sample_phi is globally registered
 * base_u[i]           incoming atlas value; used ONLY to report disorder
 * delam_keys          sorted ((uint64_t)lo << 32) | hi component pairs to exempt
 */
int AtlasRadialOrder_build(Arena_T arena,
                           const AtlasStripSample *samples,
                           const int32_t *sample_plane,
                           const int32_t *sample_component,
                           const double *sample_phi,
                           const int32_t *phase_domain,
                           size_t nsamples,
                           const double *base_u,
                           const ScaffoldCalib *cal,
                           const uint64_t *delam_keys,
                           size_t ndelam_keys,
                           const AtlasRadialOrderOptions *opts,
                           AtlasRadialOrderSet *out);

typedef struct {
    size_t rounds;
    size_t relaxations;
    size_t moved_variables;
    size_t anchor_blocked;   /* bounds dropped: they would move an anchor */
    size_t cycle_broken;     /* ray bounds dropped to break a positive cycle */
    double max_shift;
    double total_shift;
} AtlasRadialOrderBootstrapStats;

/*
 * Make x feasible for a general set of difference bounds x[hi] - x[lo] >= lower.
 *
 * MonotoneQp requires a feasible start and will not repair one.  This is the
 * usual longest-path push (SPFA), raising x[hi] until every bound holds, which
 * terminates exactly when the constraint graph has no positive cycle.
 *
 * It is deliberately fail-soft.  A positive cycle means contradictory ordering
 * evidence -- in practice one connected component that physically fuses two
 * wraps, which no injective warp can fix -- so instead of aborting, the ray
 * bound that closed the cycle is disabled and counted.  Stroke monotonicity is
 * never disabled: it is upstream truth.  A large cycle_broken count is the
 * signal to cut the fused component, not to loosen the solve.
 *
 * disabled[] (optional, nbounds entries) receives 1 for each dropped bound.
 */
int AtlasRadialOrder_bootstrap(Arena_T arena,
                               const MonotoneQpBound *bounds, size_t nbounds,
                               const MonotoneQpAnchor *anchors, size_t nanchors,
                               size_t nvar,
                               double tolerance,
                               double *x,
                               uint8_t *disabled,
                               AtlasRadialOrderBootstrapStats *stats);

int AtlasRadialOrder_selftest(void);

#endif
