#ifndef UNION_FIND_INCLUDED
#define UNION_FIND_INCLUDED

#include <stdint.h>
#include "arena.h"

typedef struct {
    int32_t *parent;   /* parent[i] = parent of i, or i if root */
    int32_t *rank;     /* rank[i] = upper bound on height */
    int32_t  count;    /* current number of distinct components */
    int32_t  n;        /* total elements */
} UnionFind;

/* Create UF for n elements. Arena-allocated. */
UnionFind UF_new(Arena_T arena, int32_t n);

/* Find with path splitting (iterative, no recursion). */
int32_t uf_find(UnionFind *uf, int32_t x);

/* Union by rank. Decrements uf->count on merge. */
void uf_union(UnionFind *uf, int32_t a, int32_t b);

#endif
