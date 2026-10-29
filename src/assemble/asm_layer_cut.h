#ifndef ASM_LAYER_CUT_INCLUDED
#define ASM_LAYER_CUT_INCLUDED

#include <stddef.h>
#include <stdio.h>

#include "asm_types.h"
#include "asm_axis.h"

/* ============================================================================
 * asm_layer_cut.h -- drop the joins that fuse wraps, from 3-D evidence alone.
 *
 * relate's layer pairs are charts stacked one or more wraps apart: a ray
 * along one chart's normal crosses the other.  On a single sheet two such
 * charts are joined only by going around the axis (a full turn per wrap).  A
 * join path between them that winds ZERO times around the axis must cross a
 * bridge between wraps, and the region of the join graph holding it cannot
 * be flattened without laying one wrap over the other.
 *
 * Each pass searches, for every layer pair inside one join component, the
 * shortest join path whose lifted angle about the axis returns to its own
 * branch (winding 0) -- a BFS on the lift of the join graph, so a path that
 * winds once is a different state and never a witness, at any radius.  Every
 * witness path then nominates its most used join (the bridge is the
 * bottleneck of all paths through it); the nominated joins are dropped
 * (ASM_REL_DROPPED) and the witnesses through them searched again, until none
 * remains.  Nothing depends on the pose layout: the 21x21x21's pose component
 * of 45k charts overlaps itself everywhere (1.78 M contradicting pairs of
 * ~40 vox^2 each), and there clean2's layout witnesses cannot localize a
 * bridge.
 * ==========================================================================*/

typedef struct AsmLayerCutOpts {
    const AsmAxis *axis;   /* chart angles about the scroll axis; NULL disables the pass */
    int    max_hops;       /* longest witness path searched, joins */
    double min_radius;     /* charts nearer the axis have no angle and are not traversed, vox */
    int    min_hits;       /* layer pairs with fewer ray hits are not evidence */
    int    max_passes;     /* nomination passes; witnesses left after the last are reported */
} AsmLayerCutOpts;

typedef struct AsmLayerCutStats {
    size_t pairs;              /* layer pairs used as evidence (both charts traversable) */
    size_t pairs_joined;       /* of those, inside one join component at the start */
    size_t witnesses_first;    /* zero-winding join paths within max_hops before any drop */
    size_t witnesses_last;     /* ... after the last pass */
    size_t joins_first;        /* joins among traversable charts at the start */
    size_t dropped;            /* joins dropped */
    int    passes;
    double sec;
} AsmLayerCutStats;

void AsmLayerCut_defaults(AsmLayerCutOpts *o);

/* Drops fusing joins in place.  ledger (nullable) gets one CSV row per
 * dropped join.  0 = success (including "nothing to do"). */
int AsmLayerCut_run(AsmRun *run, const AsmLayerCutOpts *o, AsmLayerCutStats *st, FILE *ledger);

int AsmLayerCut_selftest(void);

#endif
