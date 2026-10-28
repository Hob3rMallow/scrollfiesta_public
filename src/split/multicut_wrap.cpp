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
#include "../common/pipeline_constants.h"

#include <algorithm>
#include <vector>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include "andres/graph/graph.hxx"
#include "andres/graph/multicut-lifted/greedy-additive.hxx"
#include "andres/graph/multicut-lifted/kernighan-lin.hxx"
#include "multicut_kl_fast.hxx"
#include "multicut_gaec_fast.hxx"

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
        /* Phase timing (stderr only): graph build, GAEC, Kernighan-Lin. */
        typedef std::chrono::steady_clock Clock;
        const Clock::time_point t_build = Clock::now();
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

        const Clock::time_point t_gaec = Clock::now();
        /* Phase 1: Greedy additive edge contraction */
        size_t min_k = (min_clusters > 1)
            ? static_cast<size_t>(min_clusters) : 1;
        /* identical to andres::...::greedyAdditiveEdgeContraction, with pooled
         * map nodes (multicut_gaec_fast.hxx; LiftedMulticut_kl_selftest) */
        std::vector<size_t> labels =
            ves_gaec::greedyAdditiveEdgeContraction(
                original, lifted, costs, min_k);

        const Clock::time_point t_kl = Clock::now();
        /* Phase 2: Kernighan-Lin refinement */
        andres::graph::multicut_lifted::KernighanLinSettings settings;
        settings.numberOfOuterIterations = 100;
        settings.epsilon = 1e-6;
        settings.introduce_new_sets = true;
        if (MULTICUT_KL_MAX_INNER > 0)   /* A/B lever; 0 = vendored (unlimited) */
            settings.numberOfInnerIterations = (size_t)MULTICUT_KL_MAX_INNER;

        /* identical result to andres::...::kernighanLin, with an O(log n)
         * move selection (multicut_kl_fast.hxx; LiftedMulticut_kl_selftest) */
        ves_kl::KlStats kls;
        labels = ves_kl::kernighanLin(
            original, lifted, costs, labels, settings, &kls);
        const Clock::time_point t_end = Clock::now();
        std::fprintf(stderr,
            "      kl detail: %ld updates, %ld moves, %ld outer: diff=%.2f setup=%.2f "
            "moves=%.2f finish=%.2f outer=%.2f s\n",
            kls.updates, kls.moves, kls.outer_iters, kls.t_diff, kls.t_setup,
            kls.t_moves, kls.t_finish, kls.t_outer);
        {
            /* A cap of K moves per update reproduces this result exactly
             * when every update's best prefix ends within K moves. */
            size_t max_best = 0;
            for (const auto &bm : kls.best_moves)
                if (bm.first > max_best) max_best = bm.first;
            std::fprintf(stderr, "      kl best: max_best=%zu |", max_best);
            const size_t caps[5] = { 1000, 4000, 16000, 64000, 256000 };
            for (size_t c : caps) {
                size_t kept = 0, over = 0;
                for (const auto &bm : kls.best_moves) {
                    kept += bm.second < c ? bm.second : c;
                    if (bm.first > c) over++;
                }
                std::fprintf(stderr, " cap%zu: moves %zu, over %zu;", c, kept, over);
            }
            std::fprintf(stderr, "\n");
        }
        std::fprintf(stderr,
            "      multicut timing: build=%.3f gaec=%.3f kl=%.3f s "
            "(nodes %d, adj %d, lifted %d)\n",
            std::chrono::duration<double>(t_gaec - t_build).count(),
            std::chrono::duration<double>(t_kl - t_gaec).count(),
            std::chrono::duration<double>(t_end - t_kl).count(),
            (int)num_nodes, (int)num_adj_edges, (int)num_lifted_edges);

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

/* ---- ves_kl::kernighanLin must reproduce the vendored Kernighan-Lin label
 * for label on mesh-like instances: a grid adjacency with attractive costs
 * from a small set (ties), lifted repulsive "overlap" edges, GAEC starts and
 * random starts. ---- */
static uint32_t kl_st_x = 0x2545F491u;
static uint32_t kl_st_next(void) {
    uint32_t x = kl_st_x; x ^= x << 13; x ^= x >> 17; x ^= x << 5; kl_st_x = x; return x;
}

int LiftedMulticut_kl_selftest(void)
{
    int fails = 0, cases = 0;
    long moved_total = 0;
    kl_st_x = 0x2545F491u;
    for (int inst = 0; inst < 12 && !fails; inst++) {
        const int W = 24 + (inst % 4) * 12, H = 20 + (inst % 3) * 10;
        const size_t n = (size_t)W * (size_t)H;
        andres::graph::Graph<> original(n), lifted(n);
        std::vector<double> costs;
        auto add = [&](size_t a, size_t b, double c) {
            const size_t e = lifted.insertEdge(a, b);
            if (e >= costs.size()) costs.resize(e + 1, 0.0);
            costs[e] += c;
        };
        static const double att[4] = { 1.0, 0.5, 2.0, 1.0 };
        for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++) {
                size_t v = (size_t)y * (size_t)W + (size_t)x;
                if (x + 1 < W) {
                    original.insertEdge(v, v + 1);
                    add(v, v + 1, att[kl_st_next() % 4]);
                }
                if (y + 1 < H) {
                    original.insertEdge(v, v + (size_t)W);
                    add(v, v + (size_t)W, att[kl_st_next() % 4]);
                }
            }
        /* repulsive lifted edges between nearby cells (overlaps), two strengths */
        const size_t nrep = n / 3;
        for (size_t r = 0; r < nrep; r++) {
            size_t a = kl_st_next() % n;
            int dx = (int)(kl_st_next() % 7) - 3, dy = (int)(kl_st_next() % 7) - 3;
            int ax = (int)(a % (size_t)W), ay = (int)(a / (size_t)W);
            int bx = ax + dx, by = ay + dy;
            if (bx < 0 || by < 0 || bx >= W || by >= H) continue;
            size_t b = (size_t)by * (size_t)W + (size_t)bx;
            if (a == b) continue;
            add(a, b, (kl_st_next() & 1) ? -1.0e6 : -3.0);
        }
        andres::graph::multicut_lifted::KernighanLinSettings settings;
        settings.numberOfOuterIterations = 100;
        settings.epsilon = 1e-6;
        settings.introduce_new_sets = true;
        for (int start = 0; start < 2 && !fails; start++) {
            std::vector<size_t> init;
            if (start == 0) {
                init = andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(
                    original, lifted, costs, (size_t)1);
                for (size_t mink : { (size_t)1, (size_t)3 }) {
                    std::vector<size_t> g_ref = andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(
                        original, lifted, costs, mink);
                    std::vector<size_t> g_got = ves_gaec::greedyAdditiveEdgeContraction(
                        original, lifted, costs, mink);
                    if (g_ref != g_got) {
                        std::fprintf(stderr, "LiftedMulticut_kl_selftest FAIL: GAEC instance %d min_k %zu\n",
                                     inst, mink);
                        fails++;
                    }
                }
            } else {
                init.resize(n);
                for (size_t v = 0; v < n; v++) init[v] = kl_st_next() % 5;
            }
            std::vector<size_t> ref = andres::graph::multicut_lifted::kernighanLin(
                original, lifted, costs, init, settings);
            std::vector<size_t> got = ves_kl::kernighanLin(
                original, lifted, costs, init, settings);
            cases++;
            for (size_t v = 0; v < n; v++) if (ref[v] != init[v]) moved_total++;
            if (ref != got) {
                std::fprintf(stderr, "LiftedMulticut_kl_selftest FAIL: instance %d start %d\n",
                             inst, start);
                fails++;
            }
        }
    }
    std::fprintf(stderr, "LiftedMulticut_kl_selftest: %s (%d instances == vendored KL and GAEC, %ld labels changed by KL)\n",
                 fails ? "FAIL" : "OK", cases, moved_total);
    return fails ? 1 : 0;
}

/* ---- ves_gaec::greedyAdditiveEdgeContraction (hash-map adjacency) must
 * reproduce the vendored GAEC label for label: grids up to 30k vertices with
 * weights drawn from small sets (many queue ties), repulsive and attractive
 * lifted edges of several ranges, several cluster lower bounds, and a hub
 * vertex lifted to every other vertex (its map outgrows the 1 MiB slab blocks). ---- */
int LiftedMulticut_gaec_selftest(void)
{
    struct Spec { int W, H, rep_div, range, hub; };
    static const Spec specs[] = {
        { 12, 10, 3, 3, 0 }, { 40, 30, 2, 3, 0 }, { 64, 48, 1, 5, 0 },
        { 120, 90, 3, 9, 0 }, { 200, 150, 2, 4, 0 }, { 250, 250, 4, 3, 1 },
    };
    static const double att[4] = { 1.0, 0.5, 2.0, 1.0 };
    static const double rep[4] = { -3.0, -1.0e6, -0.5, -3.0 };
    int fails = 0, cases = 0, kl_cases = 0;
    size_t clusters_total = 0, kl_moved = 0;
    kl_st_x = 0x9E3779B9u;
    for (const Spec &sp : specs) {
        const size_t n = (size_t)sp.W * (size_t)sp.H;
        andres::graph::Graph<> original(n), lifted(n);
        std::vector<double> costs;
        auto add = [&](size_t a, size_t b, double c) {
            const size_t e = lifted.insertEdge(a, b);
            if (e >= costs.size()) costs.resize(e + 1, 0.0);
            costs[e] += c;
        };
        for (int y = 0; y < sp.H; y++)
            for (int x = 0; x < sp.W; x++) {
                const size_t v = (size_t)y * (size_t)sp.W + (size_t)x;
                if (x + 1 < sp.W) { original.insertEdge(v, v + 1); add(v, v + 1, att[kl_st_next() % 4]); }
                if (y + 1 < sp.H) { original.insertEdge(v, v + (size_t)sp.W); add(v, v + (size_t)sp.W, att[kl_st_next() % 4]); }
            }
        const size_t nlift = n / (size_t)sp.rep_div;
        for (size_t r = 0; r < nlift; r++) {
            const size_t a = kl_st_next() % n;
            const int span = 2 * sp.range + 1;
            const int dx = (int)(kl_st_next() % (uint32_t)span) - sp.range;
            const int dy = (int)(kl_st_next() % (uint32_t)span) - sp.range;
            const int bx = (int)(a % (size_t)sp.W) + dx, by = (int)(a / (size_t)sp.W) + dy;
            if (bx < 0 || by < 0 || bx >= sp.W || by >= sp.H) continue;
            const size_t b = (size_t)by * (size_t)sp.W + (size_t)bx;
            if (a == b) continue;
            /* mostly repulsive; one in eight attractive (merges can then raise weights) */
            add(a, b, (kl_st_next() % 8 == 0) ? att[kl_st_next() % 4] : rep[kl_st_next() % 4]);
        }
        if (sp.hub)
            for (size_t v = 1; v < n; v++)
                add(0, v, (v % 3 == 0) ? -0.25 : 0.125);
        for (size_t mink : { (size_t)1, (size_t)2, (size_t)7, n / 4 }) {
            const std::vector<size_t> ref = andres::graph::multicut_lifted::greedyAdditiveEdgeContraction(
                original, lifted, costs, mink);
            const std::vector<size_t> got = ves_gaec::greedyAdditiveEdgeContraction(
                original, lifted, costs, mink);
            cases++;
            clusters_total += ref.empty() ? 0 : *std::max_element(ref.begin(), ref.end()) + 1;
            if (ref != got) {
                std::fprintf(stderr, "LiftedMulticut_gaec_selftest FAIL: %dx%d min_k %zu\n", sp.W, sp.H, mink);
                fails++;
            }
            /* the fast KL on the same instance from this GAEC start (deep 8-ary
             * trees: border capacities up to 62.5k positions) */
            if (mink <= 2) {
                andres::graph::multicut_lifted::KernighanLinSettings settings;
                settings.numberOfOuterIterations = 100;
                settings.epsilon = 1e-6;
                settings.introduce_new_sets = true;
                const std::vector<size_t> kref = andres::graph::multicut_lifted::kernighanLin(
                    original, lifted, costs, ref, settings);
                const std::vector<size_t> kgot = ves_kl::kernighanLin(original, lifted, costs, ref, settings);
                kl_cases++;
                for (size_t v = 0; v < n; v++) if (kref[v] != ref[v]) kl_moved++;
                if (kref != kgot) {
                    std::fprintf(stderr, "LiftedMulticut_gaec_selftest FAIL: KL %dx%d from GAEC min_k %zu\n",
                                 sp.W, sp.H, mink);
                    fails++;
                }
            }
        }
    }
    std::fprintf(stderr, "LiftedMulticut_gaec_selftest: %s (%d cases == vendored GAEC, %zu clusters; "
                 "%d KL cases == vendored KL, %zu labels changed by KL)\n",
                 fails ? "FAIL" : "OK", cases, clusters_total, kl_cases, kl_moved);
    return fails ? 1 : 0;
}

} /* extern "C" */
