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

/* Ball query. Returns count of points within squared radius.
 * Writes up to max_results indices into out_indices. */
size_t KDTree_ball_query(const KDTree_T tree, const float center[3],
                         float radius_sq,
                         int32_t *out_indices, size_t max_results);

#endif
