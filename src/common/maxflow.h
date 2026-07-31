#ifndef MAXFLOW_INCLUDED
#define MAXFLOW_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

/* ============================================================================
 * maxflow.h -- s-t max-flow / min-cut (Dinic) on a small directed graph, for
 * exact submodular binary-energy moves (Kolmogorov-Zabih reduction). Built
 * for the group-graph collective-shift solver: ~1e4 nodes / ~1e5 arcs solve
 * in milliseconds; this is NOT a BK-tree implementation for megapixel grids.
 *
 * Node ids 0..n-1 are the problem nodes; source/sink are implicit (managed
 * internally). Capacities are int64 -- callers scale float energies up
 * (e.g. x4096) and never worry about overflow.
 * ==========================================================================*/

typedef struct Maxflow_T *Maxflow_T;

/* Graph over n problem nodes, expecting about edge_hint pairwise arcs
 * (grows as needed). Arena-allocated. */
Maxflow_T Maxflow_new(Arena_T arena, int32_t n, size_t edge_hint);

/* Directed pairwise capacity u->v of cap, plus the reverse arc v->u of rcap
 * (rcap may be 0). Both must be >= 0. */
void Maxflow_add_edge(Maxflow_T mf, int32_t u, int32_t v,
                      int64_t cap, int64_t rcap);

/* Terminal capacities: source->v of cap_s and v->sink of cap_t (either may
 * be 0). Accumulates across calls. */
void Maxflow_add_terminal(Maxflow_T mf, int32_t v, int64_t cap_s,
                          int64_t cap_t);

/* Run Dinic. Returns the max-flow value (== min-cut cost). */
int64_t Maxflow_solve(Maxflow_T mf);

/* After solve: 1 if v is on the SOURCE side of the min cut (reachable in the
 * residual graph), else 0. */
int Maxflow_in_source_side(const Maxflow_T mf, int32_t v);

/* In-process unit tests: hand-built graphs with known cuts, plus exhaustive
 * 2^n brute-force equality on random submodular binary energies (the
 * correctness gate for the KZ move construction). Returns 0 if all pass. */
int Maxflow_selftest(void);

#endif
