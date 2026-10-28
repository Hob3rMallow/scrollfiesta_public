#pragma once
#ifndef VES_MULTICUT_KL_FAST_HXX
#define VES_MULTICUT_KL_FAST_HXX

/*
 * multicut_kl_fast.hxx -- the vendored Andres lifted-multicut Kernighan-Lin
 * (deps/src/graph/include/andres/graph/multicut-lifted/kernighan-lin.hxx)
 * with the same control flow and arithmetic, faster data structures.
 *
 * Move selection (2026-10-04).  The vendored inner loop rescans the whole
 * border for every move (compacting zero-reference entries by swap-with-last,
 * then taking the FIRST maximal difference in array order), and it keeps
 * moving until the border is empty, so one bipartition update costs
 * O(moves x border).  Here the border array evolves identically:
 *   - compaction touches only the entries whose reference count dropped to
 *     zero since the last scan, in ascending position order, with the same
 *     swap-with-last-and-recheck rule, so the array ends up identical;
 *   - a max tree over array positions returns the maximal difference with ties
 *     to the lowest position -- the first maximum of the strict '>' scan.
 *
 * Memory layout (2026-10-05; same results, fewer cache misses):
 *   - the graphs are copied once into CSR arrays in the andres adjacency order,
 *     the lifted entries carrying their cost inline, so every sum runs over the
 *     same terms in the same order;
 *   - a vertex's difference, border position, label, reference count and
 *     moved flag share one 32-byte record (the vendored code keeps five
 *     arrays); the record's label mirrors vertex_labels, which stays the
 *     authoritative copy for the outer loop and is written in lockstep;
 *   - the max tree is 8-ary (a node's children are one cache line of doubles,
 *     7 levels for 512k positions instead of 19 binary ones).
 * Every difference, reference count, move, label and objective is computed by
 * the same operations in the same order, so the labels are identical
 * (LiftedMulticut_kl_selftest compares against the vendored function).
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <stack>
#include <utility>
#include <vector>

#include "andres/graph/multicut-lifted/kernighan-lin.hxx"

namespace ves_kl {

using andres::graph::multicut_lifted::KernighanLinSettings;

/* Phase timing of one kernighanLin call (diagnostics only). */
struct KlStats {
    double t_diff = 0, t_setup = 0, t_moves = 0, t_finish = 0, t_outer = 0, t_total = 0;
    long   updates = 0, moves = 0, outer_iters = 0;
    /* per bipartition update: (moves up to the best prefix, moves made) */
    std::vector<std::pair<size_t, size_t>> best_moves;
};
using KlClock = std::chrono::steady_clock;
inline double kl_secs(KlClock::time_point a, KlClock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

/* 8-ary max tree over border positions; ties -> lowest position.  Level 0 holds
 * the leaves; node j of level k+1 covers children 8j..8j+7 of level k and
 * stores the first maximum among them (strict '>' in child order), i.e. the
 * lowest position attaining its subtree's maximum.  The tree is a pure
 * function of the leaf values, so the order of updates does not matter. */
class BorderMaxTree8 {
public:
    void reset(size_t capacity) {
        nlev_ = 0;
        size_t s = round8(capacity + 1);
        for (;;) {
            if (nlev_ == lv_.size()) lv_.emplace_back();
            Level &L = lv_[nlev_++];
            L.n = s;
            if (L.val.size() < s) { L.val.resize(s); L.pos.resize(s); }
            std::fill(L.val.begin(), L.val.begin() + (std::ptrdiff_t)s, kEmpty);
            if (s == 1) break;
            const size_t parents = (s + 7) / 8;
            s = parents == 1 ? 1 : round8(parents);
        }
        for (size_t i = 0; i < lv_[0].n; ++i) lv_[0].pos[i] = (int64_t)i;
        for (size_t k = 1; k < nlev_; ++k)
            for (size_t j = 0; j < lv_[k].n; ++j)
                lv_[k].pos[j] = 8 * j < lv_[k - 1].n ? lv_[k - 1].pos[8 * j] : -1;
    }
    /* leaf write without propagation; call build() before any query */
    void put_leaf(size_t p, double v) { lv_[0].val[p] = v; }
    void build() {
        for (size_t k = 1; k < nlev_; ++k) {
            const size_t nparents = (lv_[k - 1].n + 7) / 8;
            for (size_t j = 0; j < nparents; ++j) pull(k, j);
        }
    }
    void set(size_t p, double v) {
        lv_[0].val[p] = v;
        size_t j = p;
        for (size_t k = 1; k < nlev_; ++k) {
            j >>= 3;
            /* stop once a node's (max, position) is unchanged: everything above
             * it is then unchanged too, so the tree equals full propagation */
            if (!pull(k, j)) break;
        }
    }
    void clear(size_t p) { set(p, kEmpty); }
    /* position of the max (or -1 when every leaf is empty) */
    int64_t argmax() const { const Level &T = lv_[nlev_ - 1]; return T.val[0] == kEmpty ? -1 : T.pos[0]; }
    double max() const { return lv_[nlev_ - 1].val[0]; }
private:
    static constexpr double kEmpty = -std::numeric_limits<double>::infinity();
    struct Level { size_t n = 0; std::vector<double> val; std::vector<int64_t> pos; };
    static size_t round8(size_t s) { return (s + 7) & ~(size_t)7; }
    /* recompute node j of level k from its 8 children; false when unchanged */
    bool pull(size_t k, size_t j) {
        const Level &C = lv_[k - 1];
        const double *cv = &C.val[8 * j];
        size_t b = 0;
        for (size_t c = 1; c < 8; ++c)
            if (cv[c] > cv[b]) b = c;
        Level &P = lv_[k];
        const double nv = cv[b];
        const int64_t np = C.pos[8 * j + b];
        if (P.val[j] == nv && P.pos[j] == np) return false;
        P.val[j] = nv;
        P.pos[j] = np;
        return true;
    }
    std::vector<Level> lv_;
    size_t nlev_ = 0;
};

template<typename ORIGINAL_GRAPH, typename LIFTED_GRAPH, typename ECA, typename VERTEXLABELS>
inline
auto kernighanLin(
    ORIGINAL_GRAPH const& original_graph,
    LIFTED_GRAPH const& lifted_graph,
    ECA const& edge_costs,
    VERTEXLABELS const& input_vertex_labels,
    const KernighanLinSettings settings = KernighanLinSettings(),
    KlStats *stats = nullptr) -> VERTEXLABELS
{
    KlStats st_local;
    KlStats &st = stats ? *stats : st_local;
    const KlClock::time_point t_call = KlClock::now();
    const size_t n = original_graph.numberOfVertices();

    /* CSR copies in the andres adjacency order (costs inline for the lifted graph) */
    struct LAdj { size_t v; double c; };
    std::vector<size_t> l_off(n + 1, 0), o_off(n + 1, 0);
    std::vector<LAdj> l_adj;
    std::vector<size_t> o_adj;
    l_adj.reserve(2 * lifted_graph.numberOfEdges());
    o_adj.reserve(2 * original_graph.numberOfEdges());
    for (size_t v = 0; v < n; ++v) {
        for (auto it = lifted_graph.adjacenciesFromVertexBegin(v); it != lifted_graph.adjacenciesFromVertexEnd(v); ++it)
            l_adj.push_back(LAdj{ it->vertex(), edge_costs[it->edge()] });
        l_off[v + 1] = l_adj.size();
        for (auto it = original_graph.adjacenciesFromVertexBegin(v); it != original_graph.adjacenciesFromVertexEnd(v); ++it)
            o_adj.push_back(it->vertex());
        o_off[v + 1] = o_adj.size();
    }

    /* per-vertex hot record: differences, border_pos, (mirrored) label,
     * referenced_by, is_moved of the vendored Buffers */
    struct alignas(32) VState {
        double   diff = 0.0;
        int64_t  bpos = -1;
        size_t   label = 0;
        uint32_t ref = 0;
        uint32_t moved = 0;
    };
    std::vector<VState> vs(n);
    VERTEXLABELS vertex_labels;                 // authoritative labels (buffer.vertex_labels)
    size_t max_not_used_label = 0;
    auto set_label = [&vs, &vertex_labels](size_t v, size_t l) { vertex_labels[v] = l; vs[v].label = l; };
    auto sync_labels = [&vs, &vertex_labels, n]() { for (size_t v = 0; v < n; ++v) vs[v].label = vertex_labels[v]; };

    BorderMaxTree8 tree;
    std::vector<size_t> pending_zero;      // border positions whose ref count hit 0
    std::vector<size_t> border;

    auto update_bipartition = [&](std::vector<size_t>& A, std::vector<size_t>& B)
    {
        st.updates++;
        const KlClock::time_point t0 = KlClock::now();
        struct Move
        {
            int v { -1 };
            double difference { std::numeric_limits<double>::lowest() };
            size_t new_label;
        };

        auto gain_from_merging = .0;

        auto compute_differences = [&](const std::vector<size_t>& X, size_t label_X, size_t label_Y)
        {
            for (long int i = 0; i < static_cast<long int>(X.size()); ++i)
            {
                const size_t a = X[static_cast<size_t>(i)];
                double diffExt = .0;
                double diffInt = .0;
                size_t ref_cnt = 0;

                for (size_t q = l_off[a]; q < l_off[a + 1]; ++q)
                {
                    const auto lbl = vs[l_adj[q].v].label;

                    if (lbl == label_X)
                        diffInt += l_adj[q].c;
                    else if (lbl == label_Y)
                        diffExt += l_adj[q].c;
                }

                for (size_t q = o_off[a]; q < o_off[a + 1]; ++q)
                    if (vs[o_adj[q]].label == label_Y)
                        ++ref_cnt;

                vs[a].diff = diffExt - diffInt;
                vs[a].ref = static_cast<uint32_t>(ref_cnt);
                vs[a].moved = 0;

                gain_from_merging += diffExt;
            }
        };


        if (A.empty())
            return .0;

        auto label_A = vertex_labels[A[0]];
        auto label_B = (!B.empty()) ? vertex_labels[B[0]] : max_not_used_label;

        compute_differences(A, label_A, label_B);
        compute_differences(B, label_B, label_A);
        const KlClock::time_point t1 = KlClock::now();
        st.t_diff += kl_secs(t0, t1);

        gain_from_merging /= 2.0;

        border.clear();

        for (auto a : A)
            if (vs[a].ref > 0)
                border.push_back(a);

        for (auto b : B)
            if (vs[b].ref > 0)
                border.push_back(b);

        /* Border bookkeeping: positions, the max tree, and pending compaction. */
        tree.reset(A.size() + B.size());
        pending_zero.clear();
        for (size_t p = 0; p < border.size(); ++p) {
            vs[border[p]].bpos = static_cast<int64_t>(p);
            tree.put_leaf(p, vs[border[p]].diff);
        }
        tree.build();
        auto border_move = [&](size_t from, size_t to) {   // border[to] = border[from]
            const size_t v = border[from];
            border[to] = v;
            vs[v].bpos = static_cast<int64_t>(to);
            tree.set(to, vs[v].diff);
        };

        std::vector<Move> moves;
        double cumulative_diff = .0;
        std::pair<double, size_t> max_move { std::numeric_limits<double>::lowest(), 0 };
        const KlClock::time_point t2 = KlClock::now();
        st.t_setup += kl_secs(t1, t2);

        for (size_t k = 0; k < settings.numberOfInnerIterations; ++k)
        {
            Move m;

            if (B.empty() && k == 0)
            {
                for (auto a : A)
                    if (vs[a].diff > m.difference)
                    {
                        m.v = static_cast<int>(a);
                        m.difference = vs[a].diff;
                    }
            }
            else
            {
                /* The vendored scan: for i ascending, a zero-reference entry is
                 * swapped with the last one and position i is re-examined.  Only
                 * zero entries trigger swaps, and they are visited in ascending
                 * position order, so processing just those positions reproduces
                 * the array exactly. */
                size_t size = border.size();
                std::sort(pending_zero.begin(), pending_zero.end());
                for (size_t z : pending_zero)
                {
                    while (z < size && vs[border[z]].ref == 0)
                    {
                        const size_t gone = border[z];
                        --size;
                        if (z != size) {
                            border_move(size, z);
                            border[size] = gone;
                        }
                        vs[gone].bpos = -1;
                        tree.clear(size);
                        if (z == size) break;
                    }
                }
                pending_zero.clear();
                border.erase(border.begin() + static_cast<std::ptrdiff_t>(size), border.end());

                /* first maximal difference in array order (strict '>' from lowest) */
                const int64_t p = tree.argmax();
                if (p >= 0 && tree.max() > m.difference)
                {
                    m.v = static_cast<int>(border[static_cast<size_t>(p)]);
                    m.difference = vs[static_cast<size_t>(m.v)].diff;
                }
            }

            if (m.v == -1)
                break;

            const size_t mv = static_cast<size_t>(m.v);
            const auto old_label = vs[mv].label;

            if (old_label == label_A)
                m.new_label = label_B;
            else
                m.new_label = label_A;

            // update differences and references
            for (size_t q = l_off[mv]; q < l_off[mv + 1]; ++q)
            {
                VState &w = vs[l_adj[q].v];
                if (w.moved)
                    continue;

                const auto lbl = w.label;

                // edge to an element of the new set
                if (lbl == m.new_label)
                    w.diff -= 2.0*l_adj[q].c;
                // edge to an element of the old set
                else if (lbl == old_label)
                    w.diff += 2.0*l_adj[q].c;
                else
                    continue;

                if (w.bpos >= 0)
                    tree.set(static_cast<size_t>(w.bpos), w.diff);
            }

            for (size_t q = o_off[mv]; q < o_off[mv + 1]; ++q)
            {
                const size_t wv = o_adj[q];
                VState &w = vs[wv];
                if (w.moved)
                    continue;

                const auto lbl = w.label;

                // edge to an element of the new set
                if (lbl == m.new_label)
                {
                    --w.ref;

                    if (w.ref == 0)
                    {
                        if (w.bpos >= 0)
                            pending_zero.push_back(static_cast<size_t>(w.bpos));
                    }
                }
                // edge to an element of the old set
                else if (lbl == old_label)
                {
                    ++w.ref;

                    if (w.ref == 1)
                    {
                        const size_t p = border.size();
                        border.push_back(wv);
                        w.bpos = static_cast<int64_t>(p);
                        tree.set(p, w.diff);
                    }
                }
            }

            set_label(mv, m.new_label);
            vs[mv].ref = 0;
            vs[mv].diff = std::numeric_limits<double>::lowest();
            vs[mv].moved = 1;
            {
                const int64_t bp = vs[mv].bpos;
                if (bp >= 0) {
                    tree.set(static_cast<size_t>(bp), vs[mv].diff);
                    pending_zero.push_back(static_cast<size_t>(bp));
                }
            }
            moves.push_back(m);

            cumulative_diff += m.difference;

            if (cumulative_diff > max_move.first)
                max_move = std::make_pair(cumulative_diff, moves.size());
        }

        for (size_t v : border)
            vs[v].bpos = -1;
        const KlClock::time_point t3 = KlClock::now();
        st.t_moves += kl_secs(t2, t3);
        st.moves += static_cast<long>(moves.size());
        st.best_moves.emplace_back(max_move.second, moves.size());
        struct FinishTimer {
            KlStats &s; KlClock::time_point t;
            ~FinishTimer() { s.t_finish += kl_secs(t, KlClock::now()); }
        } finish_timer{ st, t3 };


        if (gain_from_merging > max_move.first && gain_from_merging > settings.epsilon)
        {
            A.insert(A.end(), B.begin(), B.end());

            for (auto a : A)
                set_label(a, label_A);

            for (auto b : B)
                set_label(b, label_A);

            B.clear();

            return gain_from_merging;
        }
        else if (max_move.first > settings.epsilon)
        {
            // revert some changes
            for (size_t i = max_move.second; i < moves.size(); ++i)
            {
                const size_t v = static_cast<size_t>(moves[i].v);
                vs[v].moved = 0;

                if (moves[i].new_label == label_B)
                    set_label(v, label_A);
                else
                    set_label(v, label_B);
            }

            // make sure that this is unique label
            if (B.empty())
                ++max_not_used_label;

            A.erase(std::partition(A.begin(), A.end(), [&](size_t a) { return !vs[a].moved; }), A.end());
            B.erase(std::partition(B.begin(), B.end(), [&](size_t b) { return !vs[b].moved; }), B.end());

            for (size_t i = 0; i < max_move.second; ++i)
                // move vertex to the other set
                if (moves[i].new_label == label_B)
                    B.push_back(static_cast<size_t>(moves[i].v));
                else
                    A.push_back(static_cast<size_t>(moves[i].v));

            return max_move.first;
        }
        else
            for (size_t i = 0; i < moves.size(); ++i)
                if (moves[i].new_label == label_B)
                    set_label(static_cast<size_t>(moves[i].v), label_A);
                else
                    set_label(static_cast<size_t>(moves[i].v), label_B);

        return .0;
    };

    auto compute_obj_value = [&lifted_graph, &edge_costs] (std::vector<size_t> const& labels)
    {
        double obj_value = .0;

        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        {
            auto const v0 = lifted_graph.vertexOfEdge(i, 0);
            auto const v1 = lifted_graph.vertexOfEdge(i, 1);

            if (labels[v0] != labels[v1])
                obj_value += edge_costs[i];
        }

        return obj_value;
    };

    auto form_partitions = [] (std::vector<size_t> const& labels)
    {
        auto const number_of_components = *std::max_element(labels.begin(), labels.end()) + 1;

        std::vector<std::vector<size_t>> partitions(number_of_components);

        for (size_t i = 0; i < labels.size(); ++i)
            partitions[labels[i]].push_back(i);

        return partitions;
    };

    /* connected components of equal labels on the original graph; a vertex's
     * new label is its component's rank by lowest vertex, so the DFS order
     * inside a component does not matter */
    auto project_onto_feasable_set = [&o_off, &o_adj, n] (std::vector<size_t> const& labels)
    {
        std::stack<size_t> S;
        std::vector<char> visited(n);

        std::vector<size_t> feasable_vertex_labels(labels.size());

        for (size_t i = 0, new_label = 0; i < n; ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto const label = labels[i];

                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    feasable_vertex_labels[v] = new_label;

                    for (size_t q = o_off[v]; q < o_off[v + 1]; ++q)
                        if (labels[o_adj[q]] == label && !visited[o_adj[q]])
                        {
                            S.push(o_adj[q]);
                            visited[o_adj[q]] = 1;
                        }
                }

                ++new_label;
            }

        return feasable_vertex_labels;
    };

    /* flags only (set-to-1), so the visiting order does not matter */
    auto mark_partitions_that_changed_shape = [&o_off, &o_adj, n] (std::vector<size_t> const& previous_vertex_labels, std::vector<size_t> const& current_vertex_labels)
    {
        auto const number_of_components = *std::max_element(current_vertex_labels.begin(), current_vertex_labels.end()) + 1;

        std::vector<char> changed(number_of_components);

        std::stack<size_t> S;
        std::vector<char> visited(n);

        for (size_t i = 0; i < n; ++i)
            if (!visited[i])
            {
                S.push(i);
                visited[i] = 1;

                auto const label_new = current_vertex_labels[i];
                auto const label_old = previous_vertex_labels[i];

                while (!S.empty())
                {
                    auto const v = S.top();
                    S.pop();

                    for (size_t q = o_off[v]; q < o_off[v + 1]; ++q)
                    {
                        const size_t w = o_adj[q];
                        if (
                            (previous_vertex_labels[w] == label_old && current_vertex_labels[w] != label_new)
                            ||
                            (previous_vertex_labels[w] != label_old && current_vertex_labels[w] == label_new)
                            )
                            changed[label_new] = 1;

                        if (visited[w])
                            continue;

                        if (current_vertex_labels[w] == label_new)
                        {
                            S.push(w);
                            visited[w] = 1;

                            if (previous_vertex_labels[w] != label_old)
                                changed[label_new] = 1;
                        }
                    }
                }
            }

        return changed;
    };



    vertex_labels = project_onto_feasable_set(input_vertex_labels);
    sync_labels();

    auto current_obj_value = compute_obj_value(vertex_labels);

    auto partitions = form_partitions(vertex_labels);
    max_not_used_label = partitions.size();

    auto last_good_vertex_labels = vertex_labels;
    auto vertex_labels_for_changed_thing = vertex_labels;

    std::vector<std::pair<size_t, size_t>> cut_pairs;   // (min, max) partition pairs on original edges
    for (size_t k = 0; k < settings.numberOfOuterIterations; ++k)
    {
        st.outer_iters++;
        // 1 if i-th partitioned changed since last iteration, 0 otherwise
        std::vector<char> changed(partitions.size(), 1);

        if (k > 0)
            changed = mark_partitions_that_changed_shape(vertex_labels_for_changed_thing, vertex_labels);

        /* the vendored per-partition std::set of neighbouring partitions, as one
         * sorted, deduplicated pair list: same pairs, same ascending order */
        cut_pairs.clear();
        for (size_t e = 0; e < original_graph.numberOfEdges(); ++e)
        {
            auto const v0 = vertex_labels[original_graph.vertexOfEdge(e, 0)];
            auto const v1 = vertex_labels[original_graph.vertexOfEdge(e, 1)];

            if (v0 != v1)
                cut_pairs.emplace_back(std::min(v0, v1), std::max(v0, v1));
        }
        std::sort(cut_pairs.begin(), cut_pairs.end());
        cut_pairs.erase(std::unique(cut_pairs.begin(), cut_pairs.end()), cut_pairs.end());

        for (size_t q = 0; q < cut_pairs.size(); )
        {
            const size_t i = cut_pairs[q].first;
            size_t q_end = q;
            while (q_end < cut_pairs.size() && cut_pairs[q_end].first == i)
                ++q_end;

            if (!partitions[i].empty())
                for (size_t r = q; r < q_end; ++r)
                {
                    const size_t j = cut_pairs[r].second;
                    if (!partitions[j].empty() && (changed[j] || changed[i]))
                    {
                        update_bipartition(partitions[i], partitions[j]);

                        if (partitions[i].size() == 0)
                            break;
                    }
                }

            q = q_end;
        }


        vertex_labels = project_onto_feasable_set(vertex_labels);
        sync_labels();

        partitions = form_partitions(vertex_labels);
        max_not_used_label = partitions.size();

        auto pair_updates_decrease = current_obj_value - compute_obj_value(vertex_labels);
        current_obj_value -= pair_updates_decrease;



        if (k > 0)
            changed = mark_partitions_that_changed_shape(vertex_labels_for_changed_thing, vertex_labels);
        else
            changed.resize(partitions.size(), 1);

        // try to intoduce new partitions
        if (settings.introduce_new_sets)
            for (size_t i = 0; i < partitions.size(); ++i)
            {
                if (!changed[i])
                    continue;

                while (1)
                {
                    std::vector<size_t> new_set;
                    update_bipartition(partitions[i], new_set);

                    if (new_set.empty())
                        break;
                }
            }


        vertex_labels = project_onto_feasable_set(vertex_labels);
        sync_labels();

        partitions = form_partitions(vertex_labels);
        max_not_used_label = partitions.size();

        auto new_sets_decrease = current_obj_value - compute_obj_value(vertex_labels);
        current_obj_value -= new_sets_decrease;


        if (new_sets_decrease + pair_updates_decrease < std::numeric_limits<double>::epsilon())
            break;

        bool didnt_change = true;
        for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
        {
            auto v0 = lifted_graph.vertexOfEdge(i, 0);
            auto v1 = lifted_graph.vertexOfEdge(i, 1);

            auto edge_label = vertex_labels[v0] == vertex_labels[v1] ? 0 : 1;

            if (static_cast<bool>(edge_label) != (last_good_vertex_labels[v0] != last_good_vertex_labels[v1]))
                didnt_change = false;
        }

        if (didnt_change)
            break;

        vertex_labels_for_changed_thing = last_good_vertex_labels;

        last_good_vertex_labels = vertex_labels;
    }

    st.t_total = kl_secs(t_call, KlClock::now());
    st.t_outer = st.t_total - st.t_diff - st.t_setup - st.t_moves - st.t_finish;
    return last_good_vertex_labels;
}

} // namespace ves_kl

#endif
