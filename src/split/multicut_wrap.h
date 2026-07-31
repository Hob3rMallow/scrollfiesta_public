#ifndef MULTICUT_WRAP_INCLUDED
#define MULTICUT_WRAP_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Solve lifted multicut via greedy additive edge contraction
 * followed by Kernighan-Lin refinement (Keuper et al. 2015).
 *
 * Two edge sets define two graphs over the same num_nodes vertices:
 *   - original graph (adj_*): mesh face adjacency edges only
 *   - lifted graph (lifted_*): adjacency + overlap edges
 *
 * Weight convention:
 *   positive = prefer same cluster (merge / attractive)
 *   negative = prefer different cluster (cut / repulsive)
 *
 * The solver guarantees each output cluster is connected in the original graph.
 *
 * out_labels: pre-allocated [num_nodes], filled with 0-based cluster IDs.
 * out_num_clusters: number of distinct clusters.
 * Returns 0 on success, -1 on error.
 */
int LiftedMulticut_kernighan_lin(
    int32_t        num_nodes,
    /* Original graph (adjacency only) */
    int32_t        num_adj_edges,
    const int32_t *adj_from,
    const int32_t *adj_to,
    /* Lifted graph (adjacency + overlap); weights for ALL lifted edges */
    int32_t        num_lifted_edges,
    const int32_t *lifted_from,
    const int32_t *lifted_to,
    const double  *lifted_weights,
    /* Oracle hint: minimum number of clusters (0 = no constraint) */
    int32_t        min_clusters,
    /* Output */
    int32_t       *out_labels,
    int32_t       *out_num_clusters);

/* Regression for coalesced parallel lifted edges: their weights must be
 * accumulated by graph edge id rather than shifted by input position. */
int LiftedMulticut_selftest(void);

#ifdef __cplusplus
}
#endif

#endif
