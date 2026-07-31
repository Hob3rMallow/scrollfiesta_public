#ifndef ATLAS_TURN_ADVANCE_INCLUDED
#define ATLAS_TURN_ADVANCE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* ============================================================================
 * atlas_turn_advance.h -- does u actually advance along the scroll?
 *
 * Everything upstream asks where a piece SITS.  This asks whether the field
 * has the right SIZE, which turns out to be the prior question: on the 4x5x5
 * fixture, support component 3 holds 57,714 vox of papyrus along its own
 * arclength and the atlas gives it 1,719 -- 34x compressed.  No amount of
 * moving pieces around fixes a field that never grew.
 *
 * The pinned map makes "owed" exact rather than estimated.  With
 * u(phi) = a*phi + b*phi^2/(4pi), the derivative is du/dphi = a + b*phi/(2pi)
 * = r(phi), so u IS arclength and u(phi_b) - u(phi_a) is literally the length
 * of papyrus between two angles.  Every quantity here is a comparison of
 * delivered u against that.
 *
 * TWO VIEWS, because they localise different things.
 *
 *   Bins  group a piece's samples by integer turn and take the median u of
 *         each.  The step between adjacent turns should be one local
 *         circumference.  This says HOW MUCH advance is missing.
 *
 *   Edges score the candidate graph one relation at a time.  The three kinds
 *         make different promises and must be judged against different
 *         targets: a stroke's unit-speed rows should deliver full arclength; a
 *         continuation should deliver its own measured physical gap; a
 *         cross-section is an ALIGNMENT row and deliberately delivers ZERO,
 *         which is correct only if its members really are one ruling.  A
 *         cross-section whose members span a real angle is therefore not a
 *         mild error -- it is an active instruction to collapse a turn, and
 *         the collapsed-edge count is what finds those.
 *
 * Nothing here modifies a field.  It is a measurement, and its numbers are
 * meant to be quoted.
 * ==========================================================================*/

typedef struct {
    int32_t piece;
    int32_t turn;             /* floor(sense * phi / 2pi) */
    size_t  samples;
    double  u_median;
    double  radius_median;
} AtlasTurnAdvanceBin;

typedef struct {
    size_t edges;
    size_t scored;            /* edges with a meaningful target */
    size_t collapsed;         /* real angle spanned, no u delivered */
    size_t reversed;          /* u moved against the target */
    double dphi_turns_median;
    double du_median;
    double du_expected_median;
    double ratio_median;      /* delivered / owed, per edge */
    double owed_total;        /* sum |owed| -- the advance this kind carries */
    double delivered_total;   /* sum |delivered| */
} AtlasTurnAdvanceEdgeStats;

typedef struct {
    AtlasTurnAdvanceBin *bins;
    size_t nbins;
    size_t steps;             /* adjacent-turn pairs measured */
    size_t steps_backward;    /* steps whose u moved the wrong way */
    double step_median;       /* observed |u| between adjacent turns */
    double circumference_median;   /* what one turn owes there */
    double advance_fraction;  /* the headline: step_median / circumference */
} AtlasTurnAdvanceResult;

/*
 * u, phi, radius, piece   [nsamples]; piece < 0 excludes the sample
 * npieces                 upper bound on the piece id
 * sense                   winding sense from the pinned calibration (+/-1)
 */
int AtlasTurnAdvance_bins(Arena_T arena,
                          const double *u, const double *phi,
                          const double *radius, const int32_t *piece,
                          size_t nsamples, size_t npieces, double sense,
                          AtlasTurnAdvanceResult *out);

/*
 * Score one kind of relation.
 *
 * target[i] (optional) is the u difference the relation itself asserts, used
 * for continuations, which carry a measured physical gap rather than a spiral
 * arc.  Pass NULL to score against the spiral, or a zero array to score an
 * alignment row -- with `expect_zero` set, an edge is judged by how much angle
 * it spans while asserting no u change, which is the collapse signature.
 */
int AtlasTurnAdvance_edges(Arena_T arena,
                           const double *u, const double *phi,
                           const int32_t *edge_a, const int32_t *edge_b,
                           const double *target, size_t nedges,
                           size_t nsamples,
                           double spiral_a, double spiral_b,
                           int expect_zero,
                           AtlasTurnAdvanceEdgeStats *out);

int AtlasTurnAdvance_selftest(void);

#endif
