/*
 * multicut_wrap.cpp -- C wrapper for Andres lifted multicut solver.
 *
 * Uses greedy additive edge contraction for initial clustering,
 * then Kernighan-Lin local search for refinement.
 *
 * The Andres library (deps/src/graph/) is header-only C++.
 */
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "multicut_wrap.h"

#include <vector>
#include <cstddef>
#include <cstdint>

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"

extern "C" {

int LiftedMulticut_kernighan_lin(
    int32_t        num_nodes,
    int32_t        num_adj_edges,
    const int32_t *adj_from,
    const int32_t *adj_to,
    int32_t        num_lifted_edges,
    const int32_t *lifted_from,
    const int32_t *lifted_to,
    const double  *lifted_weights,
    int32_t        min_clusters,
    int32_t       *out_labels,
    int32_t       *out_num_clusters)
{
    if (num_nodes <= 0 || num_lifted_edges <= 0) {
        *out_num_clusters = 1;
        for (int32_t i = 0; i < num_nodes; i++)
            out_labels[i] = 0;
        return 0;
    }

    try {
        /* Build original graph (adjacency only - defines connectivity) */
        andres::graph::Graph<> original(static_cast<size_t>(num_nodes));
        for (int32_t i = 0; i < num_adj_edges; i++) {
            original.insertEdge(
                static_cast<size_t>(adj_from[i]),
                static_cast<size_t>(adj_to[i]));
        }

        /* Build lifted graph (adjacency + overlap edges) */
        andres::graph::Graph<> lifted(static_cast<size_t>(num_nodes));
        std::vector<double> costs;
        costs.reserve(static_cast<size_t>(num_lifted_edges));
        for (int32_t i = 0; i < num_lifted_edges; i++) {
            const size_t edge = lifted.insertEdge(
                static_cast<size_t>(lifted_from[i]),
                static_cast<size_t>(lifted_to[i]));
            /* Graph<> coalesces parallel edges.  The caller deliberately can
             * supply the same pair once as an attractive original edge and
             * again as a repulsive lifted overlap.  Indexing costs by input
             * position silently shifted every later weight after the first
             * duplicate.  Accumulate by the graph's returned edge id so the
             * combined objective is represented exactly. */
            if (edge >= costs.size()) costs.resize(edge + 1, 0.0);
            costs[edge] += lifted_weights[i];
        }
        if (costs.size() != lifted.numberOfEdges()) return -1;

        /* Phase 1: Greedy additive edge contraction */
        size_t min_k = (min_clusters > 1)
            ? static_cast<size_t>(min_clusters) : 1;
        std::vector<size_t> labels =
            andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(
                original, lifted, costs, min_k);

        /* Phase 2: Kernighan-Lin refinement */
        andres::graph::multicut_lifted::KernighanLinSettings settings;
        settings.numberOfOuterIterations = 100;
        settings.epsilon = 1e-6;
        settings.introduce_new_sets = true;

        labels = andres::graph::multicut_lifted::kernighanLin(
            original, lifted, costs, labels, settings);

        /* Copy results and compute cluster count */
        int32_t max_label = 0;
        for (size_t i = 0; i < static_cast<size_t>(num_nodes); i++) {
            out_labels[i] = static_cast<int32_t>(labels[i]);
            if (out_labels[i] > max_label)
                max_label = out_labels[i];
        }
        *out_num_clusters = max_label + 1;
        return 0;
    }
    catch (...) {
        *out_num_clusters = 1;
        for (int32_t i = 0; i < num_nodes; i++)
            out_labels[i] = 0;
        return -1;
    }
}

int LiftedMulticut_selftest(void)
{
    const int32_t adj_from[2] = {0, 1};
    const int32_t adj_to[2] = {1, 2};
    const int32_t lifted_from[3] = {0, 0, 1};
    const int32_t lifted_to[3] = {1, 1, 2};
    const double lifted_weight[3] = {10.0, -30.0, 10.0};
    int32_t label[3] = {-1, -1, -1};
    int32_t clusters = 0;
    int rc = LiftedMulticut_kernighan_lin(
        3, 2, adj_from, adj_to, 3, lifted_from, lifted_to,
        lifted_weight, 0, label, &clusters);
    if (rc != 0 || clusters != 2) return 1;
    if (label[0] == label[1] || label[1] != label[2]) return 1;
    return 0;
}

} /* extern "C" */
