#ifndef LANE_TURN_MRF_INCLUDED
#define LANE_TURN_MRF_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Integer turn offsets for the ribbon's reconstruction lanes.
 *
 * WHY A SECOND MRF (2026-09-01).  The winding registration's MRF has one site
 * per mesh connected component and is handed the same relation set that built
 * its own spanning forest, so on the 4x5x5 it converges in one round with zero
 * label changes -- a tautology.  Meanwhile the ribbon fit measures 12,452
 * places where two claims land in the SAME (row, u-cell) while sitting a
 * median 15.9 vox apart in 3-D, i.e. 1.7 wraps.  That is direct, independent
 * evidence the certificate is wrong, and it was being printed and discarded.
 *
 * This module reuses the winding MRF verbatim -- ordinal integer labels, Huber
 * unary, sum of weight * |(k_b - k_a) - target| pairwise, GCO alpha-beta swap,
 * Gibbs abstention -- at LANE granularity, where that evidence lives:
 *
 *   sites    the reconstruction lanes (branch-aware output identities)
 *   labels   integer turn offsets, bounded, centred on the carried certificate
 *   unary    Huber trust on 0, so a lane only moves when evidence outweighs it
 *   edges    ORDER  two lanes co-claiming a cell at a radial separation that
 *                   snaps to an integer number of pitches: target = that
 *                   integer, so the pair is pushed onto adjacent turns
 *            SAME   a certified chain continuation crossing a lane boundary:
 *                   target = 0, priced like the register's continuations
 *
 * The caller supplies the aggregated edges and applies the resulting shift; a
 * turn is worth a different arclength at every radius, so converting a label
 * into a u displacement is deliberately NOT this module's business. */

typedef struct {
    int32_t a, b;      /* lane ids, a != b */
    int32_t target;    /* preferred shift[b] - shift[a] */
    double  weight;    /* non-negative evidence mass */
} LaneTurnEdge;

typedef struct {
    size_t  sites, edges, labels, rounds;
    size_t  changed, abstained, moved;
    int32_t shift_min, shift_max;
    double  energy_before, energy_after;
    double  mean_confidence;
} LaneTurnStats;

/* Bound on how far one lane may move.  The evidence is a radial ORDER
 * observation, which is only trustworthy for a small integer number of wraps:
 * WR_ORDER_MULTI in the register is compiled out at k >= 2 precisely because a
 * ray's radial neighbour is often a different sheet.  Two turns is already
 * generous for a pair the lattice says are colliding. */
#define LANE_TURN_MAX_SHIFT 2

/* Solve for one integer turn offset per lane.
 *
 * `lane_support` (nlanes entries, may be NULL) is the evidence mass behind
 * each lane; the largest-supported lane is pinned so the global gauge cannot
 * drift.  `trust_weight` prices the Huber unary that holds every lane at 0.
 *
 * out_shift receives an arena-owned nlanes array.  Returns 0 on success. */
int LaneTurn_solve(Arena_T arena,
                   const LaneTurnEdge *edges, size_t nedges,
                   const double *lane_support, size_t nlanes,
                   double trust_weight,
                   int32_t **out_shift,
                   LaneTurnStats *stats);

int LaneTurn_selftest(void);

#endif
