#include "kdtree.h"
#include <assert.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ves_omp.h"

/* 16-byte node: point + original index.
 * Split axis = depth % 3. 4 nodes per 64-byte cacheline. */
typedef struct {
    float   point[3];      /* copied point coordinates */
    int32_t orig_index;    /* original index for result mapping */
} KDNode;

struct KDTree_T {
    KDNode *nodes;
    size_t  n;
    void   *self;
};

/* Swap two index entries */
static void swap_idx(int32_t *indices, size_t a, size_t b)
{
    int32_t t = indices[a];
    indices[a] = indices[b];
    indices[b] = t;
}

/* Size of the left subtree of the root of a COMPLETE binary tree with n nodes
 * (last level filled left-to-right). The heap-layout build below MUST split the
 * range at this size, NOT at n/2: n/2 leaves interior heap slots unwritten for
 * many n (e.g. n=5), and those slots keep the memset-zeroed (0,0,0)/orig_index 0
 * value, which the nearest/ball searches then treat as a real point at the origin.
 * Splitting by the complete-tree left size fills slots 0..n-1 with no gaps. */
static size_t complete_left_size(size_t n)
{
    if (n <= 1) return 0;
    size_t h = 0;                              /* floor(log2(n)) */
    while (((size_t)1 << (h + 1)) <= n) h++;
    size_t upper     = ((size_t)1 << h) - 1;   /* nodes in the full levels above the last */
    size_t last      = n - upper;              /* nodes present in the last level */
    size_t last_left = (size_t)1 << (h - 1);   /* capacity of the last level's left half */
    return (((size_t)1 << (h - 1)) - 1) + (last < last_left ? last : last_left);
}

/* Partition indices by median of axis. Uses nth_element-style partition.
 * Returns median index position. */
static size_t partition_median(const float *points, int32_t *indices,
                                size_t lo, size_t hi, int axis)
{
    size_t mid = lo + complete_left_size(hi - lo);

    /* Simple median-of-three pivot selection, then partition */
    while (lo < hi) {
        /* Choose pivot as median of lo, mid, hi-1 */
        size_t a = lo, b = mid, c = hi - 1;
        float va = points[indices[a] * 3 + axis];
        float vb = points[indices[b] * 3 + axis];
        float vc = points[indices[c] * 3 + axis];

        /* Move median to lo position as pivot */
        size_t pivot_idx;
        if (va <= vb) {
            pivot_idx = (vb <= vc) ? b : ((va <= vc) ? c : a);
        } else {
            pivot_idx = (va <= vc) ? a : ((vb <= vc) ? c : b);
        }
        swap_idx(indices, lo, pivot_idx);

        float pivot_val = points[indices[lo] * 3 + axis];
        size_t i = lo + 1;
        size_t j = hi - 1;

        while (i <= j) {
            while (i <= j && points[indices[i] * 3 + axis] <= pivot_val) {
                i++;
            }
            while (i <= j && points[indices[j] * 3 + axis] > pivot_val) {
                if (j == 0) { break; }
                j--;
            }
            if (i < j) {
                swap_idx(indices, i, j);
                i++;
                if (j > 0) { j--; }
            } else {
                break;
            }
        }
        /* j is now the last position <= pivot */
        swap_idx(indices, lo, j);

        if (j == mid) {
            return mid;
        } else if (j < mid) {
            lo = j + 1;
        } else {
            hi = j;
        }
    }
    return lo;
}

/* A subtree left for the parallel phase of the build: its heap slot, index range and depth. */
typedef struct {
    size_t node, lo, hi;
    int depth;
} KDPending;

#define KD_PARALLEL_MIN 65536   /* points from which the subtrees below KD_PARALLEL_DEPTH build in parallel */
#define KD_PARALLEL_DEPTH 6     /* 64 independent subtrees */

/* Build tree recursively into heap-layout nodes array.  With pend != NULL the
 * recursion stops at depth stop_depth and records each subtree there instead:
 * a subtree's partitions touch only its own index range and heap slots, so the
 * recorded subtrees can be built in any order (or at once) into the same tree. */
static void build_tree(const float *points, int32_t *indices,
                       KDNode *nodes, size_t node_idx,
                       size_t lo, size_t hi, int depth, size_t n_total,
                       KDPending *pend, size_t *npend, int stop_depth)
{
    if (lo >= hi || node_idx >= n_total) {
        return;
    }
    if (pend != NULL && depth == stop_depth) {
        pend[*npend].node = node_idx; pend[*npend].lo = lo; pend[*npend].hi = hi; pend[*npend].depth = depth;
        (*npend)++;
        return;
    }

    int axis = depth % 3;

    if (hi - lo == 1) {
        /* Leaf node */
        int32_t idx = indices[lo];
        nodes[node_idx].point[0] = points[idx * 3 + 0];
        nodes[node_idx].point[1] = points[idx * 3 + 1];
        nodes[node_idx].point[2] = points[idx * 3 + 2];
        nodes[node_idx].orig_index = idx;
        return;
    }

    size_t mid = partition_median(points, indices, lo, hi, axis);

    int32_t idx = indices[mid];
    nodes[node_idx].point[0] = points[idx * 3 + 0];
    nodes[node_idx].point[1] = points[idx * 3 + 1];
    nodes[node_idx].point[2] = points[idx * 3 + 2];
    nodes[node_idx].orig_index = idx;

    /* Left child: [lo, mid) */
    if (mid > lo && 2 * node_idx + 1 < n_total) {
        build_tree(points, indices, nodes, 2 * node_idx + 1,
                   lo, mid, depth + 1, n_total, pend, npend, stop_depth);
    }

    /* Right child: [mid+1, hi) */
    if (mid + 1 < hi && 2 * node_idx + 2 < n_total) {
        build_tree(points, indices, nodes, 2 * node_idx + 2,
                   mid + 1, hi, depth + 1, n_total, pend, npend, stop_depth);
    }
}

/* partition_median over elements that carry their point (a KDNode: point +
 * original index): the same comparisons and swaps, so the same permutation,
 * read from consecutive memory instead of gathered through the index array. */
static size_t partition_median_packed(KDNode *e, size_t lo, size_t hi, int axis)
{
    size_t mid = lo + complete_left_size(hi - lo);
    while (lo < hi) {
        size_t a = lo, b = mid, c = hi - 1;
        float va = e[a].point[axis], vb = e[b].point[axis], vc = e[c].point[axis];
        size_t pivot_idx;
        KDNode t;
        if (va <= vb) {
            pivot_idx = (vb <= vc) ? b : ((va <= vc) ? c : a);
        } else {
            pivot_idx = (va <= vc) ? a : ((vb <= vc) ? c : b);
        }
        t = e[lo]; e[lo] = e[pivot_idx]; e[pivot_idx] = t;
        float pivot_val = e[lo].point[axis];
        size_t i = lo + 1;
        size_t j = hi - 1;
        while (i <= j) {
            while (i <= j && e[i].point[axis] <= pivot_val) {
                i++;
            }
            while (i <= j && e[j].point[axis] > pivot_val) {
                if (j == 0) { break; }
                j--;
            }
            if (i < j) {
                t = e[i]; e[i] = e[j]; e[j] = t;
                i++;
                if (j > 0) { j--; }
            } else {
                break;
            }
        }
        t = e[lo]; e[lo] = e[j]; e[j] = t;
        if (j == mid) {
            return mid;
        } else if (j < mid) {
            lo = j + 1;
        } else {
            hi = j;
        }
    }
    return lo;
}

/* One node of build_tree over packed elements: its slot, and the subtrees
 * below it (out[0..*nout)) in build_tree's order. */
static void build_node_packed(KDNode *e, KDNode *nodes, KDPending in, size_t n_total, KDPending *out, size_t *nout)
{
    size_t mid = 0;
    *nout = 0;
    if (in.lo >= in.hi || in.node >= n_total) return;
    if (in.hi - in.lo == 1) { nodes[in.node] = e[in.lo]; return; }
    mid = partition_median_packed(e, in.lo, in.hi, in.depth % 3);
    nodes[in.node] = e[mid];
    if (mid > in.lo && 2 * in.node + 1 < n_total) {
        out[*nout].node = 2 * in.node + 1; out[*nout].lo = in.lo; out[*nout].hi = mid;
        out[*nout].depth = in.depth + 1; (*nout)++;
    }
    if (mid + 1 < in.hi && 2 * in.node + 2 < n_total) {
        out[*nout].node = 2 * in.node + 2; out[*nout].lo = mid + 1; out[*nout].hi = in.hi;
        out[*nout].depth = in.depth + 1; (*nout)++;
    }
}

static void build_tree_packed(KDNode *e, KDNode *nodes, KDPending in, size_t n_total)
{
    KDPending kids[2];
    size_t nk = 0, k = 0;
    build_node_packed(e, nodes, in, n_total, kids, &nk);
    for (k = 0; k < nk; k++) build_tree_packed(e, nodes, kids[k], n_total);
}

/* The whole tree.  Large ones over packed elements: the top KD_PARALLEL_DEPTH
 * levels one level at a time (a level's subtrees partition at once), then the
 * subtrees below in parallel (2026-10-30: two 6-9 M point trees per close-up
 * adopt round took 2.7 s serial).  Every node takes the slot and the point a
 * serial build_tree over the index array gives it. */
static void build_all(const float *points, int32_t *indices, KDNode *nodes, size_t n)
{
#ifdef _OPENMP
    if (n >= KD_PARALLEL_MIN && !omp_in_parallel()) {
        KDNode *e = (KDNode *)malloc(n * sizeof *e);
        if (e != NULL) {
            KDPending cur[(size_t)1 << KD_PARALLEL_DEPTH], kids[(size_t)2 << KD_PARALLEL_DEPTH];
            size_t ncur = 1, nk[(size_t)1 << KD_PARALLEL_DEPTH];
            long p = 0;
            int level = 0;
#pragma omp parallel for schedule(static)
            for (p = 0; p < (long)n; p++) {
                e[p].point[0] = points[(size_t)p * 3]; e[p].point[1] = points[(size_t)p * 3 + 1];
                e[p].point[2] = points[(size_t)p * 3 + 2]; e[p].orig_index = indices[p];
            }
            cur[0].node = 0; cur[0].lo = 0; cur[0].hi = n; cur[0].depth = 0;
            for (level = 0; level < KD_PARALLEL_DEPTH && ncur > 0; level++) {
                size_t k = 0, next = 0;
#pragma omp parallel for schedule(dynamic, 1)
                for (p = 0; p < (long)ncur; p++) build_node_packed(e, nodes, cur[p], n, &kids[2 * (size_t)p], &nk[p]);
                for (k = 0; k < ncur; k++) {
                    size_t c = 0;
                    for (c = 0; c < nk[k]; c++) cur[next++] = kids[2 * k + c];
                }
                ncur = next;
            }
#pragma omp parallel for schedule(dynamic, 1)
            for (p = 0; p < (long)ncur; p++) build_tree_packed(e, nodes, cur[p], n);
            free(e);
            return;
        }
    }
#endif
    build_tree(points, indices, nodes, 0, 0, n, 0, n, NULL, NULL, 0);
}

KDTree_T KDTree_new(Arena_T arena, const float *points, size_t n)
{
    assert(arena);

    KDTree_T tree;
    ARENA_NEW(arena, tree);
    tree->n = n;
    tree->self = tree;

    if (n == 0) {
        tree->nodes = NULL;
        return tree;
    }

    assert(points);

    tree->nodes = (KDNode *)ARENA_ALLOC(arena, (size_t)(n * sizeof(KDNode)));

    /* Initialize all nodes */
    memset(tree->nodes, 0, n * sizeof(KDNode));

    /* Build index array */
    int32_t *indices = (int32_t *)ARENA_ALLOC(arena,
                                               (size_t)(n * sizeof(int32_t)));
    for (size_t i = 0; i < n; i++) {
        indices[i] = (int32_t)i;
    }

    build_all(points, indices, tree->nodes, n);

    return tree;
}

/* Squared distance between two 3D points */
static inline float dist_sq(const float a[3], const float b[3])
{
    float d0 = a[0] - b[0];
    float d1 = a[1] - b[1];
    float d2 = a[2] - b[2];
    return d0 * d0 + d1 * d1 + d2 * d2;
}

/* Recursive 1-NN search (pre-2026-10-06 form, kept as the selftest reference
 * for nearest_rec below). */
static void ref_nearest_rec(const KDNode *nodes, size_t n,
                            size_t node_idx, int depth,
                            const float query[3],
                            size_t *best_idx, float *best_dist_sq)
{
    if (node_idx >= n) {
        return;
    }

    const KDNode *nd = &nodes[node_idx];
    float d = dist_sq(query, nd->point);
    if (d < *best_dist_sq) {
        *best_dist_sq = d;
        *best_idx = node_idx;
    }

    int axis = depth % 3;
    float diff = query[axis] - nd->point[axis];
    float diff_sq = diff * diff;

    /* Visit closer child first */
    size_t first  = (diff <= 0.0f) ? 2 * node_idx + 1 : 2 * node_idx + 2;
    size_t second = (diff <= 0.0f) ? 2 * node_idx + 2 : 2 * node_idx + 1;

    ref_nearest_rec(nodes, n, first, depth + 1, query, best_idx, best_dist_sq);

    /* Prune: only visit other side if splitting plane is closer than best */
    if (diff_sq < *best_dist_sq) {
        ref_nearest_rec(nodes, n, second, depth + 1, query, best_idx, best_dist_sq);
    }
}

/* One query's constant state, passed by pointer: the recursions used to carry
 * up to nine arguments, which dominated their per-node cost (2026-10-06). */
typedef struct {
    const KDNode *nodes;
    size_t        n;
    float         q[3];
    float         r2;          /* ball: radius^2; nearest: best distance^2 */
    size_t        best;        /* nearest: best node */
    int32_t      *out;         /* ball: output buffer */
    size_t        max;         /* ball: capacity */
    size_t        count;       /* ball: written */
} KDWalk;

/* 1-NN, the same traversal as ref_nearest_rec: pre-order, the closer child
 * (diff <= 0 -> left) first, the other only while diff^2 < best; the second
 * child is the loop's next iteration and the axis cycles instead of depth % 3.
 * The per-axis difference is the one dist_sq already computed (same operands,
 * same operation).  Same nodes in the same order, so the same best. */
static void nearest_rec(KDWalk *w, size_t i, int axis)
{
    const size_t n = w->n;
    while (i < n) {
        const KDNode *nd = &w->nodes[i];
        float d0 = w->q[0] - nd->point[0];
        float d1 = w->q[1] - nd->point[1];
        float d2 = w->q[2] - nd->point[2];
        float d = d0 * d0 + d1 * d1 + d2 * d2;
        if (d < w->r2) {
            w->r2 = d;
            w->best = i;
        }
        size_t l = 2 * i + 1;
        if (l >= n) return;                 /* a leaf: both children absent */
        float diff = axis == 0 ? d0 : (axis == 1 ? d1 : d2);
        int na = axis == 2 ? 0 : axis + 1;
        size_t first = (diff <= 0.0f) ? l : l + 1;
        if (first < n) nearest_rec(w, first, na);
        if (!(diff * diff < w->r2)) return;
        i = (diff <= 0.0f) ? l + 1 : l;
        axis = na;
    }
}

size_t KDTree_nearest(const KDTree_T tree, const float query[3],
                      float *out_dist_sq)
{
    assert(tree && tree->self == tree);
    assert(tree->n > 0);
    assert(query);

    KDWalk w;
    w.nodes = tree->nodes; w.n = tree->n;
    w.q[0] = query[0]; w.q[1] = query[1]; w.q[2] = query[2];
    w.r2 = FLT_MAX; w.best = 0;
    w.out = NULL; w.max = 0; w.count = 0;
    nearest_rec(&w, 0, 0);

    if (out_dist_sq) {
        *out_dist_sq = w.r2;
    }
    return (size_t)tree->nodes[w.best].orig_index;
}

/* Recursive 1-NN over points whose label differs from `exclude` (labels are
 * indexed by original point index).  Same traversal and pruning as
 * nearest_rec, so the distance found is the exact minimum of the computed
 * dist_sq over the non-excluded points that lie below the initial bound. */
static void nearest_ex_rec(const KDNode *nodes, size_t n,
                           size_t node_idx, int depth,
                           const float query[3],
                           const int32_t *labels, int32_t exclude,
                           size_t *best_idx, float *best_dist_sq)
{
    if (node_idx >= n) {
        return;
    }

    const KDNode *nd = &nodes[node_idx];
    if (labels[nd->orig_index] != exclude) {
        float d = dist_sq(query, nd->point);
        if (d < *best_dist_sq) {
            *best_dist_sq = d;
            *best_idx = node_idx;
        }
    }

    int axis = depth % 3;
    float diff = query[axis] - nd->point[axis];
    float diff_sq = diff * diff;

    size_t first  = (diff <= 0.0f) ? 2 * node_idx + 1 : 2 * node_idx + 2;
    size_t second = (diff <= 0.0f) ? 2 * node_idx + 2 : 2 * node_idx + 1;

    nearest_ex_rec(nodes, n, first, depth + 1, query, labels, exclude,
                   best_idx, best_dist_sq);
    if (diff_sq < *best_dist_sq) {
        nearest_ex_rec(nodes, n, second, depth + 1, query, labels, exclude,
                       best_idx, best_dist_sq);
    }
}

size_t KDTree_nearest_bounded(const KDTree_T tree, const float query[3],
                              float bound_sq, float *out_dist_sq)
{
    assert(tree && tree->self == tree);
    assert(query);

    KDWalk w;
    w.nodes = tree->nodes; w.n = tree->n;
    w.q[0] = query[0]; w.q[1] = query[1]; w.q[2] = query[2];
    w.r2 = bound_sq; w.best = (size_t)-1;
    w.out = NULL; w.max = 0; w.count = 0;
    if (tree->n > 0) nearest_rec(&w, 0, 0);
    if (out_dist_sq) *out_dist_sq = w.r2;
    return w.best == (size_t)-1 ? (size_t)-1 : (size_t)tree->nodes[w.best].orig_index;
}

size_t KDTree_nearest_excluding(const KDTree_T tree, const float query[3],
                                const int32_t *labels, int32_t exclude,
                                float bound_dist_sq, float *out_dist_sq)
{
    assert(tree && tree->self == tree);
    assert(query && labels);

    size_t best_idx = (size_t)-1;
    float best_dsq = bound_dist_sq;
    if (tree->n > 0)
        nearest_ex_rec(tree->nodes, tree->n, 0, 0, query, labels, exclude,
                       &best_idx, &best_dsq);
    if (out_dist_sq) {
        *out_dist_sq = best_dsq;
    }
    return best_idx == (size_t)-1 ? (size_t)-1
                                  : (size_t)tree->nodes[best_idx].orig_index;
}

/* Recursive ball query.  It stops as soon as max_results hits are written:
 * the buffer then holds the first max_results hits of this depth-first order,
 * which is exactly what the old full traversal wrote, and callers only ever
 * saw min(found, max_results), so nothing observable changes.  Traversing the
 * rest of a saturated ball was a measured CVT hot spot (2026-10-04). */
/* 2026-10-06: the context struct and the second child as the loop's next
 * iteration (as nearest_rec); the axis cycles instead of depth % 3 and the
 * per-axis difference is the one the distance already computed.  Same nodes,
 * same order, same writes as before (selftest: ref_ball_rec). */
static void ball_rec(KDWalk *w, size_t i, int axis)
{
    const size_t n = w->n;
    while (i < n && w->count < w->max) {
        const KDNode *nd = &w->nodes[i];
        float d0 = w->q[0] - nd->point[0];
        float d1 = w->q[1] - nd->point[1];
        float d2 = w->q[2] - nd->point[2];
        float d = d0 * d0 + d1 * d1 + d2 * d2;
        if (d <= w->r2) {
            w->out[w->count] = nd->orig_index;    /* count < max here */
            w->count++;
        }
        size_t l = 2 * i + 1;
        if (l >= n) return;                 /* a leaf: both children absent */
        float diff = axis == 0 ? d0 : (axis == 1 ? d1 : d2);
        int na = axis == 2 ? 0 : axis + 1;
        size_t first = (diff <= 0.0f) ? l : l + 1;
        if (first < n) ball_rec(w, first, na);
        /* Prune: only visit the other side if the splitting plane is within radius */
        if (!(diff * diff <= w->r2)) return;
        i = (diff <= 0.0f) ? l + 1 : l;
        axis = na;
    }
}

size_t KDTree_ball_query(const KDTree_T tree, const float center[3],
                         float radius_sq,
                         int32_t *out_indices, size_t max_results)
{
    assert(tree && tree->self == tree);
    assert(center);

    if (tree->n == 0) {
        return 0;
    }

    KDWalk w;
    w.nodes = tree->nodes; w.n = tree->n;
    w.q[0] = center[0]; w.q[1] = center[1]; w.q[2] = center[2];
    w.r2 = radius_sq; w.best = 0;
    w.out = out_indices; w.max = max_results; w.count = 0;
    ball_rec(&w, 0, 0);

    /* count <= max_results: the number of results written */
    return w.count;
}

/* ---- selftest: the early-exit traversal against the pre-2026-10-04 full
 * traversal, kept verbatim here as the reference. ---- */
static void ref_ball_rec(const KDNode *nodes, size_t n,
                         size_t node_idx, int depth,
                         const float center[3], float radius_sq,
                         int32_t *out_indices, size_t max_results,
                         size_t *count)
{
    if (node_idx >= n) {
        return;
    }

    const KDNode *nd = &nodes[node_idx];
    float d = dist_sq(center, nd->point);
    if (d <= radius_sq) {
        if (*count < max_results) {
            out_indices[*count] = nd->orig_index;
        }
        (*count)++;
    }

    int axis = depth % 3;
    float diff = center[axis] - nd->point[axis];
    float diff_sq = diff * diff;

    size_t first  = (diff <= 0.0f) ? 2 * node_idx + 1 : 2 * node_idx + 2;
    size_t second = (diff <= 0.0f) ? 2 * node_idx + 2 : 2 * node_idx + 1;

    ref_ball_rec(nodes, n, first, depth + 1, center, radius_sq,
                 out_indices, max_results, count);
    if (diff_sq <= radius_sq) {
        ref_ball_rec(nodes, n, second, depth + 1, center, radius_sq,
                     out_indices, max_results, count);
    }
}

static uint32_t kd_st_rng(uint32_t *s)
{
    uint32_t x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *s = x;
    return x;
}

static float kd_st_unit(uint32_t *s)
{
    return (float)(kd_st_rng(s) >> 8) * (1.0f / 16777216.0f);
}

int KDTree_selftest(void)
{
    static const size_t sizes[] = { 0, 1, 2, 3, 5, 17, 100, 1000, 5000 };
    Arena_T a = Arena_new();
    uint32_t rng = 0x5EEDu;
    int fails = 0;
    size_t checks = 0;
    {   /* the parallel build writes the serial tree slot for slot (with exact duplicates) */
        static const size_t big[] = { KD_PARALLEL_MIN, KD_PARALLEL_MIN + 12345, 4 * KD_PARALLEL_MIN + 7 };
        for (size_t bi = 0; bi < sizeof big / sizeof big[0]; bi++) {
            size_t n = big[bi];
            float *pts = (float *)ARENA_ALLOC(a, n * 3 * sizeof(float));
            KDNode *serial = (KDNode *)ARENA_CALLOC(a, n, sizeof(KDNode));
            KDNode *parallel = (KDNode *)ARENA_CALLOC(a, n, sizeof(KDNode));
            int32_t *ia = (int32_t *)ARENA_ALLOC(a, n * sizeof(int32_t));
            int32_t *ib = (int32_t *)ARENA_ALLOC(a, n * sizeof(int32_t));
            for (size_t i = 0; i < n * 3; i++) {
                pts[i] = (i % 3 == 0) ? 2.0f * kd_st_unit(&rng) : 50.0f * kd_st_unit(&rng);
                if (i >= 30 && (i / 3) % 11 == 0) pts[i] = pts[i - 30];
            }
            for (size_t i = 0; i < n; i++) ia[i] = ib[i] = (int32_t)i;
            build_tree(pts, ia, serial, 0, 0, n, 0, n, NULL, NULL, 0);
            build_all(pts, ib, parallel, n);
            if (memcmp(serial, parallel, n * sizeof(KDNode)) != 0) {
                fprintf(stderr, "  kdtree selftest FAIL: parallel build of %zu points differs from the serial build\n", n);
                fails++;
            }
        }
    }
    for (size_t si = 0; si < sizeof(sizes) / sizeof(sizes[0]) && !fails; si++) {
        size_t n = sizes[si];
        float *pts = (float *)ARENA_ALLOC(a, (size_t)((n ? n : 1) * 3 * sizeof(float)));
        for (size_t i = 0; i < n * 3; i++) {
            /* a sheet-like slab plus exact duplicates every 7th point */
            pts[i] = (i % 3 == 0) ? 2.0f * kd_st_unit(&rng) : 50.0f * kd_st_unit(&rng);
            if (i >= 21 && (i / 3) % 7 == 0) pts[i] = pts[i - 21];
        }
        KDTree_T t = KDTree_new(a, pts, n);
        size_t caps[] = { 0, 1, 2, 7, 64, 1024, n, n + 1 };
        size_t bufn = n + 16;
        int32_t *got = (int32_t *)ARENA_ALLOC(a, (size_t)(bufn * sizeof(int32_t)));
        int32_t *ref = (int32_t *)ARENA_ALLOC(a, (size_t)(bufn * sizeof(int32_t)));
        for (int q = 0; q < 60 && !fails; q++) {
            float c[3] = { 2.0f * kd_st_unit(&rng), 50.0f * kd_st_unit(&rng),
                           50.0f * kd_st_unit(&rng) };
            if (q % 5 == 0 && n > 0) {          /* query exactly at a data point */
                size_t k = (size_t)kd_st_rng(&rng) % n;
                c[0] = pts[k * 3]; c[1] = pts[k * 3 + 1]; c[2] = pts[k * 3 + 2];
            }
            float r = (q % 4 == 0) ? 1.0e4f : 30.0f * kd_st_unit(&rng);
            float r2 = r * r;
            for (size_t ci = 0; ci < sizeof(caps) / sizeof(caps[0]); ci++) {
                size_t cap = caps[ci];
                if (cap > bufn) cap = bufn;
                for (size_t i = 0; i < bufn; i++) { got[i] = -7; ref[i] = -7; }
                size_t ng = KDTree_ball_query(t, c, r2, got, cap);
                size_t rc = 0;
                if (n > 0)
                    ref_ball_rec(t->nodes, t->n, 0, 0, c, r2, ref, cap, &rc);
                size_t nr = rc < cap ? rc : cap;
                checks++;
                if (ng != nr || memcmp(got, ref, bufn * sizeof(int32_t)) != 0) {
                    fprintf(stderr, "KDTree_selftest FAIL: n=%zu cap=%zu got=%zu ref=%zu\n",
                            n, cap, ng, nr);
                    fails++;
                    break;
                }
            }
            if (n > 0) {                      /* 1-NN against the recursive reference */
                float dg = -1.0f, dr = FLT_MAX;
                size_t ig = KDTree_nearest(t, c, &dg), ir = 0;
                ref_nearest_rec(t->nodes, t->n, 0, 0, c, &ir, &dr);
                checks++;
                if (ig != (size_t)t->nodes[ir].orig_index ||
                    memcmp(&dg, &dr, sizeof(float)) != 0) {
                    fprintf(stderr, "KDTree_selftest FAIL: nearest n=%zu got %zu ref %d\n",
                            n, ig, (int)t->nodes[ir].orig_index);
                    fails++;
                }
                /* bounded by one ulp above a random point's distance (the
                 * duplicates every 7th point make exact ties) */
                size_t k = (size_t)kd_st_rng(&rng) % n;
                float e0 = c[0] - pts[k * 3], e1 = c[1] - pts[k * 3 + 1], e2 = c[2] - pts[k * 3 + 2];
                float bd = e0 * e0 + e1 * e1 + e2 * e2, db = -1.0f;
                size_t ib = KDTree_nearest_bounded(t, c, nextafterf(bd, FLT_MAX), &db);
                checks++;
                if (ib != ig || memcmp(&db, &dg, sizeof(float)) != 0) {
                    fprintf(stderr, "KDTree_selftest FAIL: bounded nearest n=%zu got %zu want %zu\n",
                            n, ib, ig);
                    fails++;
                }
            }
        }
    }
    Arena_dispose(&a);
    fprintf(stderr, "KDTree_selftest: %s (%zu capped ball and 1-NN queries vs the recursive references)\n",
            fails ? "FAIL" : "OK", checks);
    return fails ? -1 : 0;
}
