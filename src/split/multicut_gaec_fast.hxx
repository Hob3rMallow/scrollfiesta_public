#pragma once
#ifndef VES_MULTICUT_GAEC_FAST_HXX
#define VES_MULTICUT_GAEC_FAST_HXX

/*
 * multicut_gaec_fast.hxx -- the vendored Andres greedy additive edge
 * contraction (deps/src/graph/include/andres/graph/multicut-lifted/
 * greedy-additive.hxx) with its per-vertex std::map adjacency replaced by
 * open-addressing hash maps.  Labels are identical:
 *   - the priority queue is the same std::priority_queue<Edge>, fed the same
 *     edges with the same weights and editions in the same order, so its pops
 *     (including the order among equal weights) are the same;
 *   - the only adjacency iteration whose ORDER reaches the queue is the loop
 *     over the merged vertex's lifted neighbours, which std::map visits in
 *     ascending key order -- here that vertex's entries are copied and sorted
 *     by key first;
 *   - every other use is membership, lookup, size() or set-like insertion and
 *     erasure, none of which depends on iteration order;
 *   - every floating-point operation (p.second + t) is unchanged.
 * 2026-10-05: the map operations (O(log d) pointer chasing per lookup) were
 * the cost -- GAEC 33 s on a 374k-node component.  The previous pooled-node
 * std::map version is superseded; the vendored function stays the reference
 * in LiftedMulticut_kl_selftest.
 */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <new>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "andres/partition.hxx"

namespace ves_gaec {

/* Power-of-two blocks from 4 MiB chunks with a free list per size class;
 * blocks above 1 MiB go to operator new.  Chunks are freed with the pool. */
class SlabPool {
public:
    SlabPool() { for (auto &f : free_) f = nullptr; }
    ~SlabPool() { for (void *c : chunks_) ::operator delete(c); }
    SlabPool(const SlabPool&) = delete;
    SlabPool& operator=(const SlabPool&) = delete;
    void *get(size_t bytes) {
        if (bytes > kBig) return ::operator new(bytes);
        const int k = lg(bytes);
        if (free_[k]) { void *p = free_[k]; free_[k] = *static_cast<void**>(p); return p; }
        if (static_cast<size_t>(end_ - next_) < bytes) {
            char *c = static_cast<char*>(::operator new(kChunk));
            chunks_.push_back(c);
            next_ = c; end_ = c + kChunk;
        }
        void *p = next_;
        next_ += bytes;
        return p;
    }
    void put(void *p, size_t bytes) {
        if (bytes > kBig) { ::operator delete(p); return; }
        const int k = lg(bytes);
        *static_cast<void**>(p) = free_[k];
        free_[k] = p;
    }
private:
    static constexpr size_t kChunk = (size_t)4 << 20;
    static constexpr size_t kBig = (size_t)1 << 20;
    static int lg(size_t bytes) { int k = 0; while (((size_t)1 << k) < bytes) ++k; return k; }
    void *free_[64];
    std::vector<void*> chunks_;
    char *next_ = nullptr, *end_ = nullptr;
};

/* uint32 key -> V, linear probing, capacity a power of two (>= 4), load <= 3/4,
 * backward-shift deletion (no tombstones).  Key 0xFFFFFFFF marks an empty slot,
 * so vertex ids must be below it (checked by the caller). */
template <class V>
class HMap {
public:
    struct Slot { uint32_t key; V val; };
    static constexpr uint32_t kEmpty = 0xFFFFFFFFu;

    uint32_t size() const { return size_; }

    const V *find(uint32_t k) const {
        if (size_ == 0) return nullptr;
        for (uint32_t i = home(k);; i = (i + 1) & mask_) {
            const Slot &s = slots_[i];
            if (s.key == k) return &s.val;
            if (s.key == kEmpty) return nullptr;
        }
    }
    bool contains(uint32_t k) const { return find(k) != nullptr; }

    /* the value of k, inserted as `init` when absent; valid until the next insert */
    V &get_or_insert(uint32_t k, V init, SlabPool &pool) {
        if ((size_ + 1) * 4 > cap_ * 3) rehash(cap_ ? cap_ * 2 : 4, pool);
        for (uint32_t i = home(k);; i = (i + 1) & mask_) {
            Slot &s = slots_[i];
            if (s.key == k) return s.val;
            if (s.key == kEmpty) { s.key = k; s.val = init; ++size_; return s.val; }
        }
    }
    void assign(uint32_t k, V v, SlabPool &pool) { get_or_insert(k, v, pool) = v; }

    void erase(uint32_t k) {
        if (size_ == 0) return;
        uint32_t i = home(k);
        for (;; i = (i + 1) & mask_) {
            if (slots_[i].key == kEmpty) return;
            if (slots_[i].key == k) break;
        }
        for (uint32_t j = i;;) {
            j = (j + 1) & mask_;
            if (slots_[j].key == kEmpty) break;
            const uint32_t h = home(slots_[j].key);
            /* the entry at j may fill the hole at i unless its home lies
             * cyclically in (i, j] */
            const bool stays = (i <= j) ? (h > i && h <= j) : (h > i || h <= j);
            if (!stays) { slots_[i] = slots_[j]; i = j; }
        }
        slots_[i].key = kEmpty;
        --size_;
    }

    template <class F> void for_each(F f) const {
        for (uint32_t i = 0; i < cap_; ++i)
            if (slots_[i].key != kEmpty) f(slots_[i].key, slots_[i].val);
    }

    void reserve(uint32_t n, SlabPool &pool) {
        uint32_t c = 4;
        while (c * 3 < n * 4) c <<= 1;
        if (c > cap_) rehash(c, pool);
    }
    void release(SlabPool &pool) {
        if (slots_) pool.put(slots_, (size_t)cap_ * sizeof(Slot));
        slots_ = nullptr; cap_ = size_ = mask_ = 0; shift_ = 32;
    }

private:
    uint32_t home(uint32_t k) const { return (uint32_t)(k * 0x9E3779B1u) >> shift_; }
    void rehash(uint32_t ncap, SlabPool &pool) {
        Slot *old = slots_;
        const uint32_t ocap = cap_;
        slots_ = static_cast<Slot*>(pool.get((size_t)ncap * sizeof(Slot)));
        for (uint32_t i = 0; i < ncap; ++i) slots_[i].key = kEmpty;
        cap_ = ncap; mask_ = ncap - 1;
        uint32_t l = 0;
        while ((1u << l) < ncap) ++l;
        shift_ = 32 - l;
        for (uint32_t i = 0; i < ocap; ++i) {
            if (old[i].key == kEmpty) continue;
            uint32_t j = home(old[i].key);
            while (slots_[j].key != kEmpty) j = (j + 1) & mask_;
            slots_[j] = old[i];
        }
        if (old) pool.put(old, (size_t)ocap * sizeof(Slot));
    }
    Slot *slots_ = nullptr;
    uint32_t cap_ = 0, size_ = 0, mask_ = 0, shift_ = 32;
};

/// Graph decomposition by greedy additive edge contraction (vendored semantics).
template<typename ORIGGRAPH, typename LIFTGRAPH, typename ECA>
std::vector<size_t> greedyAdditiveEdgeContraction(
    ORIGGRAPH const& original_graph,
    LIFTGRAPH const& lifted_graph,
    ECA const& edge_costs,
    size_t number_of_clusters_lower_bound = 1,
    size_t number_of_clusters_upper_bound = std::numeric_limits<size_t>::max()
)
{
    using W = typename ECA::value_type;
    const size_t n = original_graph.numberOfVertices();
    if (n >= (size_t)0xFFFFFFFFu)
        throw std::length_error("ves_gaec: vertex ids must fit below 0xFFFFFFFF");

    struct Edge
    {
        Edge(size_t _a, size_t _b, W _w)
        {
            if (_a > _b)
                std::swap(_a, _b);

            a = _a;
            b = _b;

            w = _w;
        }

        size_t a;
        size_t b;
        size_t edition;
        W w;

        bool operator <(Edge const& other) const
        {
            return w < other.w;
        }
    };

    SlabPool pool;
    std::vector<HMap<uint8_t>> orig(n);      // original-graph copy: membership only
    std::vector<HMap<W>> lift(n);            // lifted-graph copy: weights
    std::vector<HMap<size_t>> editions(n);   // edge_editions[a][b], a < b; never erased
    struct Release {                         // hand big blocks back before the pool goes
        SlabPool &p; std::vector<HMap<uint8_t>> &o; std::vector<HMap<W>> &l; std::vector<HMap<size_t>> &e;
        ~Release() {
            for (auto &m : o) m.release(p);
            for (auto &m : l) m.release(p);
            for (auto &m : e) m.release(p);
        }
    } release{ pool, orig, lift, editions };
    for (size_t v = 0; v < n; ++v) {
        orig[v].reserve((uint32_t)original_graph.numberOfEdgesFromVertex(v), pool);
        lift[v].reserve((uint32_t)lifted_graph.numberOfEdgesFromVertex(v), pool);
    }
    auto orig_exists = [&orig](size_t a, size_t b) { return orig[a].contains((uint32_t)b); };

    std::priority_queue<Edge> Q;

    for (size_t i = 0; i < original_graph.numberOfEdges(); ++i)
    {
        auto const a = original_graph.vertexOfEdge(i, 0);
        auto const b = original_graph.vertexOfEdge(i, 1);

        orig[a].assign((uint32_t)b, 1, pool);
        orig[b].assign((uint32_t)a, 1, pool);
    }

    for (size_t i = 0; i < lifted_graph.numberOfEdges(); ++i)
    {
        auto const a = lifted_graph.vertexOfEdge(i, 0);
        auto const b = lifted_graph.vertexOfEdge(i, 1);

        lift[a].assign((uint32_t)b, edge_costs[i], pool);
        lift[b].assign((uint32_t)a, edge_costs[i], pool);

        if (orig_exists(a, b))
        {
            auto e = Edge(a, b, edge_costs[i]);
            e.edition = ++editions[e.a].get_or_insert((uint32_t)e.b, 0, pool);

            Q.push(e);
        }
    }

    andres::Partition<size_t> partition(n);
    std::vector<std::pair<uint32_t, W>> nb;   // the merged vertex's lifted neighbours, ascending
    while (!Q.empty() && partition.numberOfSets() > number_of_clusters_lower_bound)
    {
        auto edge = Q.top();
        Q.pop();

        if (!orig_exists(edge.a, edge.b) || edge.edition < editions[edge.a].get_or_insert((uint32_t)edge.b, 0, pool))
            continue;

        if (partition.numberOfSets() <= number_of_clusters_upper_bound && edge.w < W())
            break;

        auto stable_vertex = edge.a;
        auto merge_vertex = edge.b;

        if (lift[stable_vertex].size() < lift[merge_vertex].size())
            std::swap(stable_vertex, merge_vertex);

        partition.merge(stable_vertex, merge_vertex);

        const uint32_t s32 = (uint32_t)stable_vertex, m32 = (uint32_t)merge_vertex;
        orig[merge_vertex].for_each([&](uint32_t p, uint8_t) {
            if (p == s32)
                return;
            orig[stable_vertex].assign(p, 1, pool);
            orig[p].assign(s32, 1, pool);
        });

        orig[merge_vertex].for_each([&](uint32_t p, uint8_t) { orig[p].erase(m32); });
        orig[merge_vertex].release(pool);

        nb.clear();
        lift[merge_vertex].for_each([&nb](uint32_t p, W w) { nb.emplace_back(p, w); });
        std::sort(nb.begin(), nb.end(),
                  [](const std::pair<uint32_t, W> &x, const std::pair<uint32_t, W> &y) { return x.first < y.first; });

        for (auto const& p : nb)
        {
            if (p.first == s32)
                continue;

            auto t = W();

            if (const W *ws = lift[stable_vertex].find(p.first))
                t = *ws;

            lift[stable_vertex].assign(p.first, p.second + t, pool);
            lift[p.first].assign(s32, p.second + t, pool);

            if (orig_exists(stable_vertex, p.first))
            {
                auto e = Edge(stable_vertex, p.first, p.second + t);
                e.edition = ++editions[e.a].get_or_insert((uint32_t)e.b, 0, pool);

                Q.push(e);
            }
        }

        for (auto const& p : nb)
            lift[p.first].erase(m32);
        lift[merge_vertex].release(pool);
    }

    /* labels by ascending representative, as the vendored std::map `indices` */
    std::vector<char> is_rep(n, 0);
    for (size_t v = 0; v < n; ++v)
        is_rep[partition.find(v)] = 1;

    std::vector<size_t> rep_label(n, 0);
    size_t cnt = 0;
    for (size_t r = 0; r < n; ++r)
        if (is_rep[r])
            rep_label[r] = cnt++;

    std::vector<size_t> vertex_labels(n);
    for (size_t v = 0; v < n; ++v)
        vertex_labels[v] = rep_label[partition.find(v)];

    return vertex_labels;
}

} // namespace ves_gaec

#endif
