#include "../common/ves_platform.h"
#include "../common/run_ctx.h"

#include "bridge_cut.h"
#include "../common/csr.h"
#include "../common/bfs.h"
#include "../common/kdtree.h"
#include "../common/pca.h"
#include "../common/union_find.h"
#include "../common/pipeline_constants.h"
#include "oracle.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef VESUVIUS_DEBUG
#include <stdio.h>
#endif

/* ---------- Local constants ---------- */

#define MIN_COMPONENT_VERTS  200
#define MIN_SPLIT_FACES       10
#define MAX_RESULT_PIECES    128

/* Vertex-split indexing: v_in = 2v, v_out = 2v + 1 */
#define V_IN(v)  ((int32_t)(v) * 2)
#define V_OUT(v) ((int32_t)(v) * 2 + 1)

/* ---------- Types ---------- */

typedef struct {
    int32_t  n_nodes;
    int32_t  n_edges;
    int32_t *offset;    /* [n_nodes + 1] — CSR row offsets */
    int32_t *target;    /* [n_edges] — edge target node */
    int32_t *cap;       /* [n_edges] — edge capacity */
    int32_t *flow_arr;  /* [n_edges] — current flow (init 0) */
    int32_t *rev;       /* [n_edges] — reverse edge index */
    int32_t  S;         /* super-source node index */
    int32_t  T;         /* super-sink node index */
} FlowNet;

typedef struct {
    ComponentMesh *meshes;
    size_t         count;
    size_t         capacity;
} ResultCollector;

/* ---------- is_bridge ---------- */

static int is_bridge(int32_t flow_val, float height, float split_ratio)
{
    if (split_ratio < 0.20f) return 0;

    float fh = (height > 0.0f) ? (float)flow_val / height : 1e30f;

    if (fh <= 0.36f) return 1;

    if (fh <= 1.3f && split_ratio >= 0.40f && height >= 20.0f) return 1;

    return 0;
}

/* ---------- build_flow_net ---------- */

static FlowNet build_flow_net(Arena_T arena,
                               const int32_t *mesh_offset,
                               const int32_t *mesh_target,
                               int32_t n_verts,
                               const uint8_t *is_protected,
                               const uint8_t *is_source,
                               const uint8_t *is_sink)
{
    FlowNet fn;
    int32_t N = n_verts;
    fn.n_nodes = 2 * N + 2;
    fn.S = 2 * N;
    fn.T = 2 * N + 1;

    /* Pass 1: count out-degree for each flow node */
    int32_t *deg = (int32_t *)ARENA_CALLOC(arena, (long)fn.n_nodes,
                                            (long)sizeof(int32_t));

    /* Internal edges: v_in->v_out (fwd), v_out->v_in (rev) */
    for (int32_t v = 0; v < N; v++) {
        deg[V_IN(v)]++;
        deg[V_OUT(v)]++;
    }

    /* Cross edges: for each CSR directed edge u->v, add
     * V_OUT(u)->V_IN(v) (fwd) and V_IN(v)->V_OUT(u) (rev) */
    for (int32_t u = 0; u < N; u++) {
        for (int32_t e = mesh_offset[u]; e < mesh_offset[u + 1]; e++) {
            int32_t v = mesh_target[e];
            deg[V_OUT(u)]++;
            deg[V_IN(v)]++;
        }
    }

    /* Super-source/sink edges */
    for (int32_t v = 0; v < N; v++) {
        if (is_source[v]) {
            deg[fn.S]++;
            deg[V_IN(v)]++;
        }
        if (is_sink[v]) {
            deg[V_OUT(v)]++;
            deg[fn.T]++;
        }
    }

    /* Pass 2: prefix sum -> offsets */
    fn.offset = (int32_t *)ARENA_ALLOC(arena,
                    (long)((size_t)fn.n_nodes + 1) * (long)sizeof(int32_t));
    fn.offset[0] = 0;
    for (int32_t i = 0; i < fn.n_nodes; i++) {
        fn.offset[i + 1] = fn.offset[i] + deg[i];
    }
    fn.n_edges = fn.offset[fn.n_nodes];

    /* Allocate edge arrays */
    fn.target   = (int32_t *)ARENA_ALLOC(arena,
                      (long)fn.n_edges * (long)sizeof(int32_t));
    fn.cap      = (int32_t *)ARENA_ALLOC(arena,
                      (long)fn.n_edges * (long)sizeof(int32_t));
    fn.flow_arr = (int32_t *)ARENA_CALLOC(arena,
                      (long)fn.n_edges, (long)sizeof(int32_t));
    fn.rev      = (int32_t *)ARENA_ALLOC(arena,
                      (long)fn.n_edges * (long)sizeof(int32_t));

    /* Pass 3: fill edges using write cursors */
    int32_t *cursor = (int32_t *)ARENA_ALLOC(arena,
                          (long)fn.n_nodes * (long)sizeof(int32_t));
    memcpy(cursor, fn.offset, (size_t)fn.n_nodes * sizeof(int32_t));

    /* Internal edges */
    for (int32_t v = 0; v < N; v++) {
        int32_t cap_val = is_protected[v] ? FLOW_INF : 1;
        int32_t fi = cursor[V_IN(v)]++;
        int32_t ri = cursor[V_OUT(v)]++;
        fn.target[fi] = V_OUT(v);  fn.cap[fi] = cap_val;
        fn.target[ri] = V_IN(v);   fn.cap[ri] = 0;
        fn.rev[fi] = ri;
        fn.rev[ri] = fi;
    }

    /* Cross edges */
    for (int32_t u = 0; u < N; u++) {
        for (int32_t e = mesh_offset[u]; e < mesh_offset[u + 1]; e++) {
            int32_t v = mesh_target[e];
            int32_t fi = cursor[V_OUT(u)]++;
            int32_t ri = cursor[V_IN(v)]++;
            fn.target[fi] = V_IN(v);   fn.cap[fi] = FLOW_INF;
            fn.target[ri] = V_OUT(u);  fn.cap[ri] = 0;
            fn.rev[fi] = ri;
            fn.rev[ri] = fi;
        }
    }

    /* Super-source edges */
    for (int32_t v = 0; v < N; v++) {
        if (is_source[v]) {
            int32_t fi = cursor[fn.S]++;
            int32_t ri = cursor[V_IN(v)]++;
            fn.target[fi] = V_IN(v);  fn.cap[fi] = FLOW_INF;
            fn.target[ri] = fn.S;     fn.cap[ri] = 0;
            fn.rev[fi] = ri;
            fn.rev[ri] = fi;
        }
    }

    /* Super-sink edges */
    for (int32_t v = 0; v < N; v++) {
        if (is_sink[v]) {
            int32_t fi = cursor[V_OUT(v)]++;
            int32_t ri = cursor[fn.T]++;
            fn.target[fi] = fn.T;       fn.cap[fi] = FLOW_INF;
            fn.target[ri] = V_OUT(v);   fn.cap[ri] = 0;
            fn.rev[fi] = ri;
            fn.rev[ri] = fi;
        }
    }

    return fn;
}

/* ---------- edmonds_karp ---------- */

static int32_t edmonds_karp(FlowNet *fn, int32_t *parent_edge,
                             int32_t *queue)
{
    int32_t n_nodes = fn->n_nodes;
    int32_t total_flow = 0;

    while (total_flow < MAX_FLOW_LIMIT) {
        /* BFS for augmenting path in residual graph */
        int found = BFS_augmenting_path(fn->offset, fn->target,
                                         fn->cap, fn->flow_arr,
                                         n_nodes, fn->S, fn->T,
                                         parent_edge, queue);
        if (!found) break;

        /* Find bottleneck along the path */
        int32_t bottleneck = FLOW_INF;
        int32_t u = fn->T;
        while (u != fn->S) {
            int32_t e = parent_edge[u];
            int32_t residual = fn->cap[e] - fn->flow_arr[e];
            if (residual < bottleneck) bottleneck = residual;
            u = fn->target[fn->rev[e]]; /* source of edge e */
        }

        /* Augment flow along the path */
        u = fn->T;
        while (u != fn->S) {
            int32_t e = parent_edge[u];
            fn->flow_arr[e] += bottleneck;
            fn->flow_arr[fn->rev[e]] -= bottleneck;
            u = fn->target[fn->rev[e]];
        }

        total_flow += bottleneck;
    }

    return total_flow;
}

/* ---------- extract_reachable ---------- */

/* Full BFS on residual graph from S (no early exit at sink).
 * Sets reachable[v] = 1 if V_OUT(v) is reachable from S. */
static void extract_reachable(const FlowNet *fn, int32_t n_verts,
                               int32_t *bfs_visited, int32_t *bfs_queue,
                               uint8_t *reachable)
{
    int32_t n_nodes = fn->n_nodes;

    /* Reset visited */
    memset(bfs_visited, 0, (size_t)n_nodes * sizeof(int32_t));

    bfs_visited[fn->S] = 1;
    int32_t head = 0, tail = 0;
    bfs_queue[tail++] = fn->S;

    while (head < tail) {
        int32_t v = bfs_queue[head++];
        for (int32_t e = fn->offset[v]; e < fn->offset[v + 1]; e++) {
            int32_t nb = fn->target[e];
            if (!bfs_visited[nb] && fn->cap[e] - fn->flow_arr[e] > 0) {
                bfs_visited[nb] = 1;
                bfs_queue[tail++] = nb;
            }
        }
    }

    /* Map to original vertices: check V_OUT (matching Python) */
    for (int32_t v = 0; v < n_verts; v++) {
        reachable[v] = bfs_visited[V_OUT(v)] ? (uint8_t)1 : (uint8_t)0;
    }
}

/* ---------- detect_bridge ---------- */

/* Detect bridge via Edmonds-Karp max-flow with vertex splitting.
 * All scratch is arena-allocated; caller must save/restore.
 * Returns 0 on success with out_ parameters filled.
 * Returns -1 if mesh is too small or degenerate. */
static int detect_bridge(Arena_T arena, const ComponentMesh *mesh,
                          int32_t *out_flow, float *out_height,
                          float *out_split_ratio, uint8_t **out_reachable)
{
    int32_t N = (int32_t)mesh->nv;
    size_t nf = mesh->nf;

    if (mesh->nv < MIN_COMPONENT_VERTS || nf == 0) return -1;

#ifdef VESUVIUS_DEBUG
    double t0, t1;
    t0 = ves_clock_sec();
#endif

    /* 1. PCA normal + projection */
    float normal[3] = {0.0f, 0.0f, 0.0f};
    float centroid[3] = {0.0f, 0.0f, 0.0f};
    PCA_normal(mesh->verts, mesh->nv, normal, centroid);

    float *proj = (float *)ARENA_ALLOC(arena,
                      (long)(mesh->nv * sizeof(float)));
    float proj_min = 0.0f, proj_max = 0.0f;
    PCA_project(mesh->verts, mesh->nv, normal, proj,
                &proj_min, &proj_max);

    float height = proj_max - proj_min;
    *out_height = height;

    /* 2. Source/sink seed selection */
    int32_t source_seed = 0, sink_seed = 0;
    for (int32_t i = 1; i < N; i++) {
        if (proj[i] < proj[source_seed]) source_seed = i;
        if (proj[i] > proj[sink_seed])   sink_seed = i;
    }

    /* 3. Build mesh adjacency */
    CSR_T adj = CSR_from_faces(arena, mesh->faces, nf, mesh->nv);
    const int32_t *mesh_off = CSR_offset(adj);
    const int32_t *mesh_tgt = CSR_target(adj);

    /* 4. k-ring BFS expansion */
    uint32_t *dist = (uint32_t *)ARENA_ALLOC(arena,
                         (long)(mesh->nv * sizeof(uint32_t)));
    int32_t *ring_buf = (int32_t *)ARENA_ALLOC(arena,
                         (long)(mesh->nv * sizeof(int32_t)));

    size_t n_src = BFS_kring(adj, source_seed, SEED_RING,
                              ring_buf, mesh->nv, dist);

    uint8_t *is_source = (uint8_t *)ARENA_CALLOC(arena, (long)mesh->nv, 1L);
    for (size_t i = 0; i < n_src; i++) {
        is_source[ring_buf[i]] = 1;
    }

    size_t n_snk = BFS_kring(adj, sink_seed, SEED_RING,
                              ring_buf, mesh->nv, dist);

    uint8_t *is_sink = (uint8_t *)ARENA_CALLOC(arena, (long)mesh->nv, 1L);
    for (size_t i = 0; i < n_snk; i++) {
        is_sink[ring_buf[i]] = 1;
    }

    /* Remove overlap */
    size_t actual_src = 0, actual_snk = 0;
    for (int32_t v = 0; v < N; v++) {
        if (is_source[v] && is_sink[v]) {
            is_source[v] = 0;
            is_sink[v] = 0;
        }
        if (is_source[v]) actual_src++;
        if (is_sink[v])   actual_snk++;
    }

    if (actual_src == 0 || actual_snk == 0) return -1;

    /* Protected = source OR sink */
    uint8_t *is_protected = (uint8_t *)ARENA_CALLOC(arena,
                                (long)mesh->nv, 1L);
    for (int32_t v = 0; v < N; v++) {
        is_protected[v] = (is_source[v] || is_sink[v]) ? (uint8_t)1
                                                        : (uint8_t)0;
    }

#ifdef VESUVIUS_DEBUG
    t1 = ves_clock_sec();
    fprintf(stderr, "      detect: adj+seeds %.3f s (%zu src, %zu snk)\n",
            t1 - t0,
            actual_src, actual_snk);

    t0 = ves_clock_sec();
#endif

    /* 5. Build vertex-split flow network */
    FlowNet fn = build_flow_net(arena, mesh_off, mesh_tgt, N,
                                 is_protected, is_source, is_sink);

    /* BFS workspace (reused across all augmenting path iterations) */
    int32_t *parent_edge = (int32_t *)ARENA_ALLOC(arena,
                               (long)fn.n_nodes * (long)sizeof(int32_t));
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena,
                         (long)fn.n_nodes * (long)sizeof(int32_t));

    /* 6. Edmonds-Karp max-flow */
    int32_t flow_val = edmonds_karp(&fn, parent_edge, queue);
    *out_flow = flow_val;

#ifdef VESUVIUS_DEBUG
    t1 = ves_clock_sec();
    fprintf(stderr, "      detect: flow=%d (%d nodes, %d edges) %.3f s\n",
            (int)flow_val, (int)fn.n_nodes, (int)fn.n_edges,
            t1 - t0);
    t0 = ves_clock_sec();
#endif

    /* 7. Extract reachable set via full BFS on residual graph */
    /* Reuse parent_edge as visited array, queue as BFS queue */
    uint8_t *reachable = (uint8_t *)ARENA_ALLOC(arena, (long)mesh->nv);
    extract_reachable(&fn, N, parent_edge, queue, reachable);

    /* 8. Compute split_ratio */
    int32_t n_reach = 0;
    for (int32_t v = 0; v < N; v++) {
        if (reachable[v]) n_reach++;
    }
    int32_t n_unreach = N - n_reach;
    int32_t smaller = (n_reach < n_unreach) ? n_reach : n_unreach;
    float split_ratio = (float)smaller / (float)N;
    *out_split_ratio = split_ratio;
    *out_reachable = reachable;

#ifdef VESUVIUS_DEBUG
    t1 = ves_clock_sec();
    fprintf(stderr, "      detect: reachable %d/%d (%.3f s)\n",
            (int)n_reach, (int)N,
            t1 - t0);
#endif

    return 0;
}

/* ---------- split_mesh ---------- */

/* Split mesh at boundary with CUT_GAP_DEPTH exclusion zone.
 * Assigns faces by first vertex's reachability (matching Python). */
static int split_mesh(Arena_T arena, const ComponentMesh *mesh,
                       const uint8_t *reachable,
                       ComponentMesh *out_a, ComponentMesh *out_b)
{
    size_t N = mesh->nv;
    size_t nf = mesh->nf;

    Arena_Mark scratch = Arena_save(arena);

    /* 1. Build adjacency for boundary detection */
    CSR_T adj = CSR_from_faces(arena, mesh->faces, nf, N);
    const int32_t *off = CSR_offset(adj);
    const int32_t *tgt = CSR_target(adj);

    /* 2. Find boundary vertices */
    size_t n_boundary = 0;
    uint8_t *is_boundary = (uint8_t *)ARENA_CALLOC(arena, (long)N, 1L);

    for (size_t v = 0; v < N; v++) {
        for (int32_t e = off[v]; e < off[v + 1]; e++) {
            if (reachable[v] != reachable[(size_t)tgt[e]]) {
                is_boundary[v] = 1;
                n_boundary++;
                break;
            }
        }
    }

    /* 3. Build KD-tree on boundary, compute exclusion zone */
    uint8_t *excluded = (uint8_t *)ARENA_CALLOC(arena, (long)N, 1L);

    if (n_boundary > 0) {
        float *boundary_pts = (float *)ARENA_ALLOC(arena,
                                  (long)(n_boundary * 3 * sizeof(float)));
        size_t bi = 0;
        for (size_t v = 0; v < N; v++) {
            if (is_boundary[v]) {
                boundary_pts[bi * 3 + 0] = mesh->verts[v * 3 + 0];
                boundary_pts[bi * 3 + 1] = mesh->verts[v * 3 + 1];
                boundary_pts[bi * 3 + 2] = mesh->verts[v * 3 + 2];
                bi++;
            }
        }

        KDTree_T tree = KDTree_new(arena, boundary_pts, n_boundary);
        float gap_sq = CUT_GAP_DEPTH * CUT_GAP_DEPTH;

        for (size_t v = 0; v < N; v++) {
            float dist_sq = 0.0f;
            KDTree_nearest(tree, &mesh->verts[v * 3], &dist_sq);
            if (dist_sq < gap_sq) {
                excluded[v] = 1;
            }
        }
    }

    /* 4. Compute per-vertex side assignment:
     *    0 = excluded, 1 = source (reachable), 2 = sink */
    uint8_t *vert_side = (uint8_t *)malloc(N);
    assert(vert_side);
    for (size_t v = 0; v < N; v++) {
        if (excluded[v]) {
            vert_side[v] = 0;
        } else if (reachable[v]) {
            vert_side[v] = 1;
        } else {
            vert_side[v] = 2;
        }
    }

    /* Restore scratch (frees CSR, KD-tree, boundary arrays) */
    Arena_restore(arena, scratch);

    /* 5. Partition faces and count vertices per side */
    size_t n_src_f = 0, n_snk_f = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t v0 = mesh->faces[fi * 3 + 0];
        int32_t v1 = mesh->faces[fi * 3 + 1];
        int32_t v2 = mesh->faces[fi * 3 + 2];

        if (vert_side[v0] == 0 || vert_side[v1] == 0 ||
            vert_side[v2] == 0) continue;

        if (vert_side[v0] == 1) n_src_f++;
        else                    n_snk_f++;
    }

    /* 6. Extract side A (source/reachable) */
    {
        int32_t *old_to_new = (int32_t *)ARENA_ALLOC(arena,
                                   (long)(N * sizeof(int32_t)));
        memset(old_to_new, 0xFF, N * sizeof(int32_t)); /* -1 */

        /* Mark vertices used by source faces */
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t v0 = mesh->faces[fi * 3 + 0];
            int32_t v1 = mesh->faces[fi * 3 + 1];
            int32_t v2 = mesh->faces[fi * 3 + 2];
            if (vert_side[v0] == 0 || vert_side[v1] == 0 ||
                vert_side[v2] == 0) continue;
            if (vert_side[v0] != 1) continue;
            old_to_new[v0] = 0;
            old_to_new[v1] = 0;
            old_to_new[v2] = 0;
        }

        int32_t new_nv = 0;
        for (size_t v = 0; v < N; v++) {
            if (old_to_new[v] != -1) {
                old_to_new[v] = new_nv++;
            }
        }

        out_a->nv = (size_t)new_nv;
        out_a->nf = n_src_f;

        if (new_nv > 0) {
            out_a->verts = (float *)ARENA_ALLOC(arena,
                               (long)new_nv * 3L * (long)sizeof(float));
        } else {
            out_a->verts = NULL;
        }
        if (n_src_f > 0) {
            out_a->faces = (int32_t *)ARENA_ALLOC(arena,
                               (long)n_src_f * 3L * (long)sizeof(int32_t));
        } else {
            out_a->faces = NULL;
        }

        for (size_t v = 0; v < N; v++) {
            if (old_to_new[v] >= 0) {
                size_t ni = (size_t)old_to_new[v];
                out_a->verts[ni * 3 + 0] = mesh->verts[v * 3 + 0];
                out_a->verts[ni * 3 + 1] = mesh->verts[v * 3 + 1];
                out_a->verts[ni * 3 + 2] = mesh->verts[v * 3 + 2];
            }
        }

        size_t fi_out = 0;
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t v0 = mesh->faces[fi * 3 + 0];
            int32_t v1 = mesh->faces[fi * 3 + 1];
            int32_t v2 = mesh->faces[fi * 3 + 2];
            if (vert_side[v0] == 0 || vert_side[v1] == 0 ||
                vert_side[v2] == 0) continue;
            if (vert_side[v0] != 1) continue;
            out_a->faces[fi_out * 3 + 0] = old_to_new[v0];
            out_a->faces[fi_out * 3 + 1] = old_to_new[v1];
            out_a->faces[fi_out * 3 + 2] = old_to_new[v2];
            fi_out++;
        }

        if (new_nv >= 3) {
            PCA_normal(out_a->verts, out_a->nv,
                       out_a->pca_normal, out_a->centroid);
        } else {
            memset(out_a->pca_normal, 0, sizeof(out_a->pca_normal));
            memset(out_a->centroid, 0, sizeof(out_a->centroid));
        }
        out_a->comp_id = mesh->comp_id;
        out_a->vert_normals = NULL; /* dropped on split; pca_normal still set */
        out_a->self = out_a;
    }

    /* 7. Extract side B (sink/unreachable) */
    {
        int32_t *old_to_new = (int32_t *)ARENA_ALLOC(arena,
                                   (long)(N * sizeof(int32_t)));
        memset(old_to_new, 0xFF, N * sizeof(int32_t));

        for (size_t fi = 0; fi < nf; fi++) {
            int32_t v0 = mesh->faces[fi * 3 + 0];
            int32_t v1 = mesh->faces[fi * 3 + 1];
            int32_t v2 = mesh->faces[fi * 3 + 2];
            if (vert_side[v0] == 0 || vert_side[v1] == 0 ||
                vert_side[v2] == 0) continue;
            if (vert_side[v0] != 2) continue;
            old_to_new[v0] = 0;
            old_to_new[v1] = 0;
            old_to_new[v2] = 0;
        }

        int32_t new_nv = 0;
        for (size_t v = 0; v < N; v++) {
            if (old_to_new[v] != -1) {
                old_to_new[v] = new_nv++;
            }
        }

        out_b->nv = (size_t)new_nv;
        out_b->nf = n_snk_f;

        if (new_nv > 0) {
            out_b->verts = (float *)ARENA_ALLOC(arena,
                               (long)new_nv * 3L * (long)sizeof(float));
        } else {
            out_b->verts = NULL;
        }
        if (n_snk_f > 0) {
            out_b->faces = (int32_t *)ARENA_ALLOC(arena,
                               (long)n_snk_f * 3L * (long)sizeof(int32_t));
        } else {
            out_b->faces = NULL;
        }

        for (size_t v = 0; v < N; v++) {
            if (old_to_new[v] >= 0) {
                size_t ni = (size_t)old_to_new[v];
                out_b->verts[ni * 3 + 0] = mesh->verts[v * 3 + 0];
                out_b->verts[ni * 3 + 1] = mesh->verts[v * 3 + 1];
                out_b->verts[ni * 3 + 2] = mesh->verts[v * 3 + 2];
            }
        }

        size_t fi_out = 0;
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t v0 = mesh->faces[fi * 3 + 0];
            int32_t v1 = mesh->faces[fi * 3 + 1];
            int32_t v2 = mesh->faces[fi * 3 + 2];
            if (vert_side[v0] == 0 || vert_side[v1] == 0 ||
                vert_side[v2] == 0) continue;
            if (vert_side[v0] != 2) continue;
            out_b->faces[fi_out * 3 + 0] = old_to_new[v0];
            out_b->faces[fi_out * 3 + 1] = old_to_new[v1];
            out_b->faces[fi_out * 3 + 2] = old_to_new[v2];
            fi_out++;
        }

        if (new_nv >= 3) {
            PCA_normal(out_b->verts, out_b->nv,
                       out_b->pca_normal, out_b->centroid);
        } else {
            memset(out_b->pca_normal, 0, sizeof(out_b->pca_normal));
            memset(out_b->centroid, 0, sizeof(out_b->centroid));
        }
        out_b->comp_id = mesh->comp_id;
        out_b->vert_normals = NULL; /* dropped on split; pca_normal still set */
        out_b->self = out_b;
    }

    free(vert_side);
    return 0;
}

/* ---------- add_to_results ---------- */

static void add_to_results(ResultCollector *results,
                            const ComponentMesh *mesh)
{
    if (results->count >= results->capacity) return;
    results->meshes[results->count] = *mesh;
    results->meshes[results->count].self =
        &results->meshes[results->count];
    results->count++;
}

/* ---------- extract_mesh_components ---------- */

/* Split mesh into connected components via union-find on faces.
 * Returns number of components. If mesh is connected (1 component),
 * *out_comps is set to NULL and the caller should use the original mesh.
 * Otherwise *out_comps is an arena-allocated array of ComponentMeshes. */
static int32_t extract_mesh_components(Arena_T arena,
                                        const ComponentMesh *mesh,
                                        ComponentMesh **out_comps)
{
    int32_t N = (int32_t)mesh->nv;
    size_t nf = mesh->nf;

    if (N == 0 || nf == 0) {
        *out_comps = NULL;
        return 0;
    }

    /* Union-find on face edges */
    UnionFind uf = UF_new(arena, N);
    for (size_t fi = 0; fi < nf; fi++) {
        uf_union(&uf, mesh->faces[fi * 3 + 0], mesh->faces[fi * 3 + 1]);
        uf_union(&uf, mesh->faces[fi * 3 + 0], mesh->faces[fi * 3 + 2]);
    }

    if (uf.count <= 1) {
        *out_comps = NULL;
        return 1;
    }

    int32_t n_comps = uf.count;

    /* Map root IDs to sequential 0-based component IDs */
    int32_t *root_to_comp = (int32_t *)ARENA_ALLOC(arena,
                                 (long)N * (long)sizeof(int32_t));
    memset(root_to_comp, 0xFF, (size_t)N * sizeof(int32_t)); /* -1 */

    int32_t *vert_comp = (int32_t *)ARENA_ALLOC(arena,
                              (long)N * (long)sizeof(int32_t));
    int32_t comp_id = 0;
    for (int32_t v = 0; v < N; v++) {
        int32_t r = uf_find(&uf, v);
        if (root_to_comp[r] == -1) {
            root_to_comp[r] = comp_id++;
        }
        vert_comp[v] = root_to_comp[r];
    }

    /* Count vertices and faces per component */
    int32_t *comp_nv = (int32_t *)ARENA_CALLOC(arena,
                            (long)n_comps, (long)sizeof(int32_t));
    int32_t *comp_nf = (int32_t *)ARENA_CALLOC(arena,
                            (long)n_comps, (long)sizeof(int32_t));

    for (int32_t v = 0; v < N; v++) {
        comp_nv[vert_comp[v]]++;
    }
    for (size_t fi = 0; fi < nf; fi++) {
        comp_nf[vert_comp[mesh->faces[fi * 3]]]++;
    }

    /* Allocate component meshes (zero-init so pin_mask = NULL unless
     * explicitly set; downstream raw_snap reads pin_mask defensively). */
    ComponentMesh *comps = (ComponentMesh *)ARENA_CALLOC(arena,
                               (long)n_comps, (long)sizeof(ComponentMesh));

    /* Vertex reindex: old vertex -> new index within its component */
    int32_t *old_to_new = (int32_t *)ARENA_ALLOC(arena,
                               (long)N * (long)sizeof(int32_t));
    int32_t *vert_counter = (int32_t *)ARENA_CALLOC(arena,
                                 (long)n_comps, (long)sizeof(int32_t));

    for (int32_t v = 0; v < N; v++) {
        old_to_new[v] = vert_counter[vert_comp[v]]++;
    }

    for (int32_t c = 0; c < n_comps; c++) {
        comps[c].nv = (size_t)comp_nv[c];
        comps[c].nf = (size_t)comp_nf[c];
        if (comp_nv[c] > 0) {
            comps[c].verts = (float *)ARENA_ALLOC(arena,
                                 (long)comp_nv[c] * 3L * (long)sizeof(float));
        } else {
            comps[c].verts = NULL;
        }
        if (comp_nf[c] > 0) {
            comps[c].faces = (int32_t *)ARENA_ALLOC(arena,
                                 (long)comp_nf[c] * 3L * (long)sizeof(int32_t));
        } else {
            comps[c].faces = NULL;
        }
        comps[c].comp_id = mesh->comp_id;
        comps[c].self = &comps[c];
    }

    /* Fill vertices */
    for (int32_t v = 0; v < N; v++) {
        int32_t c = vert_comp[v];
        size_t ni = (size_t)old_to_new[v];
        comps[c].verts[ni * 3 + 0] = mesh->verts[v * 3 + 0];
        comps[c].verts[ni * 3 + 1] = mesh->verts[v * 3 + 1];
        comps[c].verts[ni * 3 + 2] = mesh->verts[v * 3 + 2];
    }

    /* Fill faces */
    int32_t *face_counter = (int32_t *)ARENA_CALLOC(arena,
                                 (long)n_comps, (long)sizeof(int32_t));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t v0 = mesh->faces[fi * 3 + 0];
        int32_t c = vert_comp[v0];
        int32_t fo = face_counter[c]++;
        comps[c].faces[fo * 3 + 0] = old_to_new[mesh->faces[fi * 3 + 0]];
        comps[c].faces[fo * 3 + 1] = old_to_new[mesh->faces[fi * 3 + 1]];
        comps[c].faces[fo * 3 + 2] = old_to_new[mesh->faces[fi * 3 + 2]];
    }

    /* Compute PCA for each component */
    for (int32_t c = 0; c < n_comps; c++) {
        if (comps[c].nv >= 3) {
            PCA_normal(comps[c].verts, comps[c].nv,
                       comps[c].pca_normal, comps[c].centroid);
        } else {
            memset(comps[c].pca_normal, 0, sizeof(comps[c].pca_normal));
            memset(comps[c].centroid, 0, sizeof(comps[c].centroid));
        }
    }

    *out_comps = comps;
    return n_comps;
}

/* ---------- process_recursive ---------- */

static void process_recursive(Arena_T arena, const ComponentMesh *mesh,
                               int depth, int grid_size,
                               double *timeout_at_sec,
                               ResultCollector *results)
{
    /* Depth limit */
    if (depth > BRIDGE_MAX_DEPTH) {
        add_to_results(results, mesh);
        return;
    }

    /* Result capacity */
    if (results->count >= results->capacity) return;

    /* Mesh too small */
    if (mesh->nv < MIN_COMPONENT_VERTS || mesh->nf < 20) {
        add_to_results(results, mesh);
        return;
    }

    /* Cancellation (API callers) piggybacks on the timeout poll site. */
    RunCtx_check();

    /* Timeout check */
    if (timeout_at_sec) {
        if (ves_clock_sec() >= *timeout_at_sec) {
            add_to_results(results, mesh);
            return;
        }
    }

    /* Split disconnected mesh components before bridge detection */
    {
        Arena_Mark conn_mark = Arena_save(arena);
        ComponentMesh *comps = NULL;
        int32_t n_comps = extract_mesh_components(arena, mesh, &comps);

        if (n_comps > 1) {
#ifdef VESUVIUS_DEBUG
            fprintf(stderr, "    split[d=%d]: mesh has %d disconnected "
                    "components\n", depth, (int)n_comps);
#endif
            for (int32_t c = 0; c < n_comps; c++) {
                process_recursive(arena, &comps[c], depth, grid_size,
                                  timeout_at_sec, results);
            }
            return;
        }
        /* Single component — free UF scratch and proceed */
        Arena_restore(arena, conn_mark);
    }

    /* Oracle look-ahead */
    int sheets = Oracle_count_sheets(arena, mesh, grid_size);
    if (sheets <= 1) {
        add_to_results(results, mesh);
        return;
    }

    /* Detect bridge */
    Arena_Mark flow_mark = Arena_save(arena);

    int32_t flow_val = 0;
    float height = 0.0f, split_ratio = 0.0f;
    uint8_t *reachable = NULL;

    int rc = detect_bridge(arena, mesh, &flow_val, &height,
                            &split_ratio, &reachable);

    if (rc != 0 || !is_bridge(flow_val, height, split_ratio)) {
#ifdef VESUVIUS_DEBUG
        if (rc == 0) {
            float fh = (height > 0.0f) ? (float)flow_val / height : 1e30f;
            fprintf(stderr, "    bridge[d=%d]: NO flow=%d h=%.1f "
                    "f/h=%.3f split=%.1f%%\n",
                    depth, (int)flow_val, (double)height,
                    (double)fh, (double)(split_ratio * 100.0f));
        }
#endif
        Arena_restore(arena, flow_mark);
        add_to_results(results, mesh);
        return;
    }

#ifdef VESUVIUS_DEBUG
    {
        float fh = (height > 0.0f) ? (float)flow_val / height : 1e30f;
        fprintf(stderr, "    bridge[d=%d]: BRIDGE flow=%d h=%.1f "
                "f/h=%.3f split=%.1f%%\n",
                depth, (int)flow_val, (double)height,
                (double)fh, (double)(split_ratio * 100.0f));
    }
#endif

    /* Copy reachable via malloc, then free flow network */
    uint8_t *reach_copy = (uint8_t *)malloc(mesh->nv);
    if (!reach_copy) {
        /* OOM fallback: keep original mesh */
        Arena_restore(arena, flow_mark);
        add_to_results(results, mesh);
        return;
    }
    memcpy(reach_copy, reachable, mesh->nv);

    Arena_restore(arena, flow_mark); /* frees flow network (~44MB) */

    /* Split mesh */
    ComponentMesh side_a, side_b;
    memset(&side_a, 0, sizeof(side_a));
    memset(&side_b, 0, sizeof(side_b));

    split_mesh(arena, mesh, reach_copy, &side_a, &side_b);
    free(reach_copy);

    /* Guard: degenerate split */
    if (side_a.nf < MIN_SPLIT_FACES || side_b.nf < MIN_SPLIT_FACES) {
#ifdef VESUVIUS_DEBUG
        fprintf(stderr, "    bridge[d=%d]: degenerate split "
                "(%zu / %zu faces)\n",
                depth, side_a.nf, side_b.nf);
#endif
        add_to_results(results, mesh);
        return;
    }

    /* Recurse on both pieces (larger first) */
    if (side_a.nf >= side_b.nf) {
        process_recursive(arena, &side_a, depth + 1, grid_size,
                          timeout_at_sec, results);
        process_recursive(arena, &side_b, depth + 1, grid_size,
                          timeout_at_sec, results);
    } else {
        process_recursive(arena, &side_b, depth + 1, grid_size,
                          timeout_at_sec, results);
        process_recursive(arena, &side_a, depth + 1, grid_size,
                          timeout_at_sec, results);
    }
}

/* ---------- BridgeCut_process ---------- */

int BridgeCut_process(Arena_T               arena,
                      const ComponentMesh   *mesh,
                      int                    grid_size,
                      double                 timeout_sec,
                      ComponentMesh        **out_meshes,
                      size_t                *out_count)
{
    assert(arena);
    assert(out_meshes);
    assert(out_count);

    /* Degenerate input */
    if (!mesh || mesh->nv == 0 || mesh->nf == 0) {
        *out_meshes = NULL;
        *out_count = 0;
        return -1;
    }

    /* Set up timeout */
    double timeout_at_sec = 0.0;
    int use_timeout = (timeout_sec > 0.0);
    if (use_timeout) {
        timeout_at_sec = ves_clock_sec() + timeout_sec;
    }

    /* Allocate result collector */
    ResultCollector results;
    results.capacity = MAX_RESULT_PIECES;
    results.count = 0;
    results.meshes = (ComponentMesh *)ARENA_ALLOC(arena,
                         (long)(results.capacity *
                                (long)sizeof(ComponentMesh)));

    /* Run recursive bridge detection */
    process_recursive(arena, mesh, 0, grid_size,
                      use_timeout ? &timeout_at_sec : NULL,
                      &results);

    /* Ensure at least 1 result */
    if (results.count == 0) {
        add_to_results(&results, mesh);
    }

    /* Filter out degenerate pieces (gap exclusion can create empty meshes) */
    {
        size_t dst = 0;
        for (size_t i = 0; i < results.count; i++) {
            if (results.meshes[i].nf >= MIN_SPLIT_FACES &&
                results.meshes[i].nv >= 3) {
                if (dst != i) {
                    results.meshes[dst] = results.meshes[i];
                }
                dst++;
            }
        }
        /* If filtering removed everything, keep the original */
        if (dst == 0) {
            add_to_results(&results, mesh);
            dst = results.count;
        }
        results.count = dst;
    }

    *out_meshes = results.meshes;
    *out_count = results.count;

    return 0;
}
