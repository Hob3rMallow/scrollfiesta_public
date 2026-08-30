#include "union_find.h"
#include <assert.h>

UnionFind UF_new(Arena_T arena, int32_t n)
{
    UnionFind uf;
    assert(arena);
    assert(n >= 0);

    uf.n     = n;
    uf.count = n;
    uf.parent = (int32_t *)ARENA_ALLOC(arena, (long)n * (long)sizeof(int32_t));
    uf.rank   = (int32_t *)ARENA_CALLOC(arena, (long)n, (long)sizeof(int32_t));

    for (int32_t i = 0; i < n; i++) {
        uf.parent[i] = i;
    }
    return uf;
}

int32_t uf_find(UnionFind *uf, int32_t x)
{
    assert(uf);
    assert(x >= 0 && x < uf->n);

    /* Path splitting: make every node point to its grandparent */
    while (uf->parent[x] != x) {
        int32_t next = uf->parent[uf->parent[x]];
        uf->parent[x] = next;
        x = next;
    }
    return x;
}

void uf_union(UnionFind *uf, int32_t a, int32_t b)
{
    assert(uf);
    a = uf_find(uf, a);
    b = uf_find(uf, b);
    if (a == b) {
        return;
    }
    /* Union by rank */
    if (uf->rank[a] < uf->rank[b]) {
        int32_t t = a;
        a = b;
        b = t;
    }
    uf->parent[b] = a;
    if (uf->rank[a] == uf->rank[b]) {
        uf->rank[a]++;
    }
    uf->count--;
}
