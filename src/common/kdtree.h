#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

typedef struct KDTree_T *KDTree_T;

/* Build balanced KD-tree from N points in R^3.
 * Points are copied into internal storage (sorted by tree order).
 * Construction: O(N log N). Arena-allocated. */
KDTree_T KDTree_new(Arena_T arena, const float *points, size_t n);

/* 1-NN query. Returns index of nearest point.
 * *out_dist_sq receives squared distance. */
size_t KDTree_nearest(const KDTree_T tree, const float query[3],
                      float *out_dist_sq);

/* 1-NN with an initial bound: the nearest point among those at squared
 * distance < bound_sq, or (size_t)-1 when none; *out_dist_sq gets its squared
 * distance (or the bound).  The traversal is KDTree_nearest's with the bound as
 * the starting best, so every point below the bound is visited in the same
 * order: with bound_sq above the distance of SOME point (e.g. nextafterf of it),
 * the result -- index and distance, ties included -- equals KDTree_nearest's. */
size_t KDTree_nearest_bounded(const KDTree_T tree, const float query[3],
                              float bound_sq, float *out_dist_sq);

/* Ball query. Writes the indices of points within squared radius into
 * out_indices, at most max_results of them (the first ones in the tree's
 * depth-first order), and returns the number written -- min(found,
 * max_results), never the total found.  A return of max_results means the
 * ball may hold more points. */
size_t KDTree_ball_query(const KDTree_T tree, const float center[3],
                         float radius_sq,
                         int32_t *out_indices, size_t max_results);

/* 1-NN among the points whose labels[orig_index] != exclude, restricted to
 * squared distance < bound_dist_sq.  Returns the original index, or
 * (size_t)-1 when no such point lies below the bound; *out_dist_sq gets the
 * squared distance found (or the bound).  Exact: same pruning as
 * KDTree_nearest, so the distance is the minimum over the eligible points. */
size_t KDTree_nearest_excluding(const KDTree_T tree, const float query[3],
                                const int32_t *labels, int32_t exclude,
                                float bound_dist_sq, float *out_dist_sq);

/* Selftest: capped ball queries against the full reference traversal.
 * Returns 0 on pass. */
int KDTree_selftest(void);

#endif
