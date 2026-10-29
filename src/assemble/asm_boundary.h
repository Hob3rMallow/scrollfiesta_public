#ifndef ASM_BOUNDARY_INCLUDED
#define ASM_BOUNDARY_INCLUDED
#include "asm_types.h"

typedef struct AsmBoundaryStats {
    size_t runs, supported_runs, ordered_samples, short_run_samples, duplicate_samples;
    double support_length;
} AsmBoundaryStats;

/* Reorder one relation's original correspondences by the actual mesh boundary
 * walks on BOTH charts. No UVs or chart IDs influence the walks. Scratch is
 * self-restoring and bounded by the two incident charts and the sample count.
 * Returns -1 on invalid/nonmanifold topology or conflicting vertex matches;
 * the correspondence array is then unchanged. Exact duplicate observations
 * are retained with run=-1 so they cannot supply additional votes.
 *
 * A complete run needs the declared minimum sample count and source arc length.
 * If one exists, every other unique valid observation remains in its own run,
 * including short runs and isolated points; ordering must not erase a difficult
 * residual. If no complete run exists, all runs remain unresolved (-1).
 * Only correspondence order and run labels change. Source geometry, UVs,
 * physical trim offsets and valid bits are preserved exactly. */
int AsmBoundary_order(Arena_T arena, const AsmChart *a, const AsmChart *b,
                      AsmCorr *samples, size_t count, AsmBoundaryStats *stats);
int AsmBoundary_selftest(void);
/* Physical quadrature on the same original boundary walks. Isolated unique
 * observations retain the mean of their two original boundary Voronoi cells.
 * No UV distances, uniform point counts, or doubled directional mass. Invalid
 * or non-unique observations return -1; they remain unresolved obligations. */
int AsmBoundary_weights(Arena_T arena, const AsmChart *a, const AsmChart *b,
                        const AsmCorr *samples, size_t count, double *weights);
/* A chart's original boundary loops, independent of any seam: built once per
 * chart (in `arena`, which must outlive its use) and shared by every seam of
 * that chart. Never NULL; an ill-formed chart gives loops that fail every
 * seam exactly as AsmBoundary_weights would. */
typedef struct AsmBoundaryLoops AsmBoundaryLoops;
AsmBoundaryLoops *AsmBoundary_loops(Arena_T arena, const AsmChart *c);
/* AsmBoundary_weights with both charts' loops already built: identical results. */
int AsmBoundary_weights_loops(Arena_T arena, const AsmChart *a, const AsmBoundaryLoops *la,
                              const AsmChart *b, const AsmBoundaryLoops *lb,
                              const AsmCorr *samples, size_t count, double *weights);
#endif
