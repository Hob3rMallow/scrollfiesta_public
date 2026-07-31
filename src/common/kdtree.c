#include "kdtree.h"
#include <assert.h>
#include <string.h>
#include <float.h>

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

/* Build tree recursively into heap-layout nodes array */
static void build_tree(const float *points, int32_t *indices,
                       KDNode *nodes, size_t node_idx,
                       size_t lo, size_t hi, int depth, size_t n_total)
{
    if (lo >= hi || node_idx >= n_total) {
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
                   lo, mid, depth + 1, n_total);
    }

    /* Right child: [mid+1, hi) */
    if (mid + 1 < hi && 2 * node_idx + 2 < n_total) {
        build_tree(points, indices, nodes, 2 * node_idx + 2,
                   mid + 1, hi, depth + 1, n_total);
    }
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

    tree->nodes = (KDNode *)ARENA_ALLOC(arena, (long)(n * sizeof(KDNode)));

    /* Initialize all nodes */
    memset(tree->nodes, 0, n * sizeof(KDNode));

    /* Build index array */
    int32_t *indices = (int32_t *)ARENA_ALLOC(arena,
                                               (long)(n * sizeof(int32_t)));
    for (size_t i = 0; i < n; i++) {
        indices[i] = (int32_t)i;
    }

    build_tree(points, indices, tree->nodes, 0, 0, n, 0, n);

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

/* Recursive 1-NN search */
static void nearest_rec(const KDNode *nodes, size_t n,
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

    nearest_rec(nodes, n, first, depth + 1, query, best_idx, best_dist_sq);

    /* Prune: only visit other side if splitting plane is closer than best */
    if (diff_sq < *best_dist_sq) {
        nearest_rec(nodes, n, second, depth + 1, query, best_idx, best_dist_sq);
    }
}

size_t KDTree_nearest(const KDTree_T tree, const float query[3],
                      float *out_dist_sq)
{
    assert(tree && tree->self == tree);
    assert(tree->n > 0);
    assert(query);

    size_t best_idx = 0;
    float best_dsq = FLT_MAX;

    nearest_rec(tree->nodes, tree->n, 0, 0, query, &best_idx, &best_dsq);

    if (out_dist_sq) {
        *out_dist_sq = best_dsq;
    }
    return (size_t)tree->nodes[best_idx].orig_index;
}

/* Recursive ball query */
static void ball_rec(const KDNode *nodes, size_t n,
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

    /* Visit closer child first */
    size_t first  = (diff <= 0.0f) ? 2 * node_idx + 1 : 2 * node_idx + 2;
    size_t second = (diff <= 0.0f) ? 2 * node_idx + 2 : 2 * node_idx + 1;

    ball_rec(nodes, n, first, depth + 1, center, radius_sq,
             out_indices, max_results, count);

    /* Prune: only visit other side if splitting plane is within radius */
    if (diff_sq <= radius_sq) {
        ball_rec(nodes, n, second, depth + 1, center, radius_sq,
                 out_indices, max_results, count);
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

    size_t count = 0;
    ball_rec(tree->nodes, tree->n, 0, 0, center, radius_sq,
             out_indices, max_results, &count);

    /* Return min(count, max_results) as the number of results written */
    return (count < max_results) ? count : max_results;
}
