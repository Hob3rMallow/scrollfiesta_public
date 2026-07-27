#include "bfs.h"
#include <assert.h>
#include <string.h>
#include <limits.h>

size_t BFS_kring(const CSR_T graph, int32_t seed, int32_t k,
                 int32_t *out_vertices, size_t max_out,
                 uint32_t *dist)
{
    assert(graph);
    assert(dist);
    assert(seed >= 0 && seed < CSR_nrows(graph));
    assert(k >= 0);

    int32_t nrows = CSR_nrows(graph);
    const int32_t *offset = CSR_offset(graph);
    const int32_t *target = CSR_target(graph);

    /* Reset dist */
    memset(dist, 0xFF, (size_t)nrows * sizeof(uint32_t));

    /* Use out_vertices as the queue (up to max_out) */
    size_t count = 0;
    size_t head = 0;

    dist[seed] = 0;
    if (count < max_out) {
        out_vertices[count] = seed;
    }
    count++;

    while (head < count) {
        size_t idx = head;
        if (idx >= max_out) {
            break;
        }
        int32_t v = out_vertices[idx];
        head++;

        uint32_t d = dist[v];
        if (d >= (uint32_t)k) {
            continue;
        }

        for (int32_t j = offset[v]; j < offset[v + 1]; j++) {
            int32_t u = target[j];
            if (dist[u] == UINT32_MAX) {
                dist[u] = d + 1;
                if (count < max_out) {
                    out_vertices[count] = u;
                }
                count++;
            }
        }
    }

    return (count < max_out) ? count : max_out;
}

void BFS_multi_source(const CSR_T graph, const int32_t *seeds,
                      size_t n_seeds, uint32_t *dist)
{
    assert(graph);
    assert(dist);

    if (n_seeds == 0) {
        return;
    }
    assert(seeds);

    int32_t nrows = CSR_nrows(graph);
    const int32_t *offset = CSR_offset(graph);
    const int32_t *target = CSR_target(graph);

    /* Use a flat queue - allocate from dist area temporarily or use
     * a separate buffer. We'll use a simple approach: store queue
     * in a portion of unvisited dist entries is not possible, so
     * we need an explicit queue. Use the fact that nrows is bounded. */
    /* For simplicity, we do BFS inline with dist as visited marker. */
    /* We need a queue - use a circular buffer approach with two passes */

    /* Simple BFS: use an allocated queue on the stack for small sizes,
     * or we accept that the caller should provide workspace.
     * For this implementation, we iterate level by level using dist. */

    /* Initialize seeds */
    for (size_t i = 0; i < n_seeds; i++) {
        assert(seeds[i] >= 0 && seeds[i] < nrows);
        dist[seeds[i]] = 0;
    }

    /* Level-synchronous BFS. Less efficient than queue-based but
     * avoids extra allocation. For pipeline use (small seed sets),
     * this is adequate. */
    int changed = 1;
    uint32_t current_dist = 0;

    while (changed) {
        changed = 0;
        for (int32_t v = 0; v < nrows; v++) {
            if (dist[v] != current_dist) {
                continue;
            }
            for (int32_t j = offset[v]; j < offset[v + 1]; j++) {
                int32_t u = target[j];
                if (dist[u] == UINT32_MAX) {
                    dist[u] = current_dist + 1;
                    changed = 1;
                }
            }
        }
        current_dist++;
    }
}

int BFS_augmenting_path(const int32_t *offset, const int32_t *target,
                        const int32_t *cap, const int32_t *flow,
                        int32_t n_nodes, int32_t source, int32_t sink,
                        int32_t *parent_edge,
                        int32_t *queue)
{
    assert(offset && target && cap && flow);
    assert(parent_edge && queue);
    assert(source >= 0 && source < n_nodes);
    assert(sink >= 0 && sink < n_nodes);

    /* Initialize parent_edge to -1 (unvisited) */
    for (int32_t i = 0; i < n_nodes; i++) {
        parent_edge[i] = -1;
    }

    parent_edge[source] = -2; /* visited but no parent */
    int32_t head = 0, tail = 0;
    queue[tail++] = source;

    while (head < tail) {
        int32_t v = queue[head++];

        for (int32_t e = offset[v]; e < offset[v + 1]; e++) {
            int32_t u = target[e];
            /* Check residual capacity: cap[e] - flow[e] > 0 */
            if (parent_edge[u] == -1 && cap[e] - flow[e] > 0) {
                parent_edge[u] = e;
                if (u == sink) {
                    return 1;
                }
                queue[tail++] = u;
            }
        }
    }

    return 0;
}
