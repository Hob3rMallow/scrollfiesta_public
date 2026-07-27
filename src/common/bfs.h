#ifndef BFS_INCLUDED
#define BFS_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "csr.h"

/* k-ring BFS: expand from seed by k hops on CSR graph.
 * Writes visited vertex indices into out_vertices.
 * Returns count of visited vertices.
 * dist[nrows]: caller-provided scratch, set to UINT32_MAX = unvisited. */
size_t BFS_kring(const CSR_T graph, int32_t seed, int32_t k,
                 int32_t *out_vertices, size_t max_out,
                 uint32_t *dist);

/* Multi-source BFS: expand from all seeds simultaneously.
 * seeds[0..n_seeds-1] are starting vertices. dist[seed] = 0.
 * Fills dist[] for all reachable vertices.
 * dist must be pre-initialized to UINT32_MAX by caller. */
void BFS_multi_source(const CSR_T graph, const int32_t *seeds,
                      size_t n_seeds, uint32_t *dist);

/* Edmonds-Karp BFS: find shortest path from source to sink in residual
 * graph. Returns 1 if path found, 0 if sink unreachable.
 * parent_edge[v] = edge index used to reach v (-1 = unvisited).
 * queue[n_nodes]: caller-provided scratch. */
int BFS_augmenting_path(const int32_t *offset, const int32_t *target,
                        const int32_t *cap, const int32_t *flow,
                        int32_t n_nodes, int32_t source, int32_t sink,
                        int32_t *parent_edge,
                        int32_t *queue);

#endif
