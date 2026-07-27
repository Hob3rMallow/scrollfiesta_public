/*
 * hole_fill.c — Interior hole detection and filling.
 *
 * 1. Find boundary edges (undirected edge count == 1)
 * 2. Chain into closed loops
 * 3. Classify interior vs boundary holes
 * 4. For each interior hole:
 *    a. PCA project to 2D
 *    b. CDT via Triangle library
 *    c. Embed Steiner points in 3D
 *    d. Smooth with cotangent Laplacian CG solve
 *    e. Orient fill faces to match mesh
 * 5. Stitch fills into mesh
 */
#include "../common/ves_platform.h"
#include "../common/run_ctx.h"

#include "../common/arena.h"
#include "../common/csr.h"
#include "../common/mesh_types.h"
#include "../common/pca.h"
#include "../common/pipeline_constants.h"
#include "clipper2_wrap.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Forward decl: used before its static definition (clang errors otherwise). */
static int loop_diag_in_mesh(int32_t va, int32_t vb,
                             const int32_t *mesh_faces, size_t mesh_nf);

/* Triangle library */
#define ANSI_DECLARATORS
#define REAL double
#define VOID int
#include "triangle.h"

#include <setjmp.h>

/* From triangle.c: longjmp guard for triexit() */
extern jmp_buf triangle_jmpbuf;

/* Per-hole debug tracing is gated behind HOLEFILL_DEBUG (off by default). On the
 * full-grid weld it fired thousands of times -- each with an fflush(stderr) --
 * and dominated the hole-fill runtime. Summary/error lines stay unconditional. */
static int hf_log_on(void) {
    static int d = -1;
    if (d < 0) d = sf_env("HOLEFILL_DEBUG") ? 1 : 0;
    return d;
}
#define HFLOG(...)   do { if (hf_log_on()) fprintf(stderr, __VA_ARGS__); } while (0)
#define HFFLUSH()    do { if (hf_log_on()) fflush(stderr); } while (0)
extern int triangle_jmpbuf_set;

#ifndef _MSC_VER
#include <signal.h>
/* Linux: sigsetjmp-based SIGSEGV guard for Triangle crashes */
static sigjmp_buf triangle_segv_jmpbuf;
static volatile sig_atomic_t triangle_segv_guard = 0;

static void triangle_segv_handler(int sig) {
    if (triangle_segv_guard) {
        siglongjmp(triangle_segv_jmpbuf, 1);
    }
    /* Not our crash — re-raise with default handler */
    signal(sig, SIG_DFL);
    raise(sig);
}
#endif

/* Constants */
#define HOLE_MARGIN        2.0f    /* (unused) interior/exterior gate removed */
#define MIN_LOOP_VERTS     3       /* minimum hole size */
#define MAX_LOOP_VERTS     500     /* maximum hole size (vertex-count cap) */
#define CG_MAX_ITER        500     /* CG solver max iterations */
#define CG_TOL             1e-6    /* CG relative tolerance */
#define MICRO_HOLE_MAX     6       /* (legacy) */
/* Degenerate-fill prune (prune_degenerate_fill): the minimal-surface smoother
 * collapses interior Steiner points of a thin-slit fill into a sub-voxel cluster,
 * leaving needle triangles. Merge Steiner endpoints of edges shorter than this and
 * drop faces below the area floor. Well under a real BPA fill edge (~0.5-1.5 vox),
 * comfortably above the observed 0.05-0.1 vox cluster gaps. */
#define FILL_DEGEN_EPS_VOX 0.25f
#define FILL_DEGEN_AREA    1e-4f
/* No-merger backstop: a 4+ closed loop wider than this (voxels) is skipped.
 * Per-component separation is the primary wrap guard (each sheet is its own
 * component post-split); this catches a pathological loop bridging a sheet
 * fold. Generous so real interior holes / welded-shut bays still fill. */
#define HOLEFILL_MAX_DIAM_VOX 48.0f

/* ------------------------------------------------------------------ */
/* Boundary edge detection                                             */
/* ------------------------------------------------------------------ */
typedef struct {
    int32_t v0, v1;  /* undirected edge (v0 < v1) */
} UEdge;

static int cmp_uedge(const void *a, const void *b)
{
    const UEdge *ea = (const UEdge *)a;
    const UEdge *eb = (const UEdge *)b;
    if (ea->v0 != eb->v0) return (ea->v0 < eb->v0) ? -1 : 1;
    return (ea->v1 < eb->v1) ? -1 : (ea->v1 > eb->v1) ? 1 : 0;
}

/* Compare face keys: 4 x int32 {min_v, mid_v, max_v, orig_idx}.
 * Sort by the first three fields (canonical vertex triple). */
static int cmp_face_key(const void *a, const void *b)
{
    const int32_t *ka = (const int32_t *)a;
    const int32_t *kb = (const int32_t *)b;
    if (ka[0] != kb[0]) return (ka[0] < kb[0]) ? -1 : 1;
    if (ka[1] != kb[1]) return (ka[1] < kb[1]) ? -1 : 1;
    if (ka[2] != kb[2]) return (ka[2] < kb[2]) ? -1 : 1;
    return 0;
}

static size_t find_boundary_edges(Arena_T arena,
                                  const int32_t *faces, size_t nf,
                                  UEdge **out_edges)
{
    size_t n_he = nf * 3;
    UEdge *all_edges = (UEdge *)ARENA_ALLOC(arena,
                         (long)n_he * (long)sizeof(UEdge));

    for (size_t f = 0; f < nf; f++) {
        int32_t v[3];
        v[0] = faces[f * 3 + 0];
        v[1] = faces[f * 3 + 1];
        v[2] = faces[f * 3 + 2];
        for (int e = 0; e < 3; e++) {
            int32_t a = v[e];
            int32_t b = v[(e + 1) % 3];
            size_t idx = f * 3 + (size_t)e;
            all_edges[idx].v0 = (a < b) ? a : b;
            all_edges[idx].v1 = (a < b) ? b : a;
        }
    }

    qsort(all_edges, n_he, sizeof(UEdge), cmp_uedge);

    /* Count boundary edges (appear exactly once) */
    size_t n_bdry = 0;
    size_t i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && all_edges[j].v0 == all_edges[i].v0 &&
               all_edges[j].v1 == all_edges[i].v1) {
            j++;
        }
        if (j - i == 1) n_bdry++;
        i = j;
    }

    UEdge *bdry = (UEdge *)ARENA_ALLOC(arena,
                     (long)(n_bdry + 1) * (long)sizeof(UEdge));
    size_t bi = 0;
    i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && all_edges[j].v0 == all_edges[i].v0 &&
               all_edges[j].v1 == all_edges[i].v1) {
            j++;
        }
        if (j - i == 1) {
            bdry[bi++] = all_edges[i];
        }
        i = j;
    }

    *out_edges = bdry;
    return n_bdry;
}

/* ------------------------------------------------------------------ */
/* chain_into_loops — walk boundary edges into closed loops            */
/*                                                                     */
/* 1. Build boundary adjacency                                         */
/* 2. Walk maximal chains (extend in both directions from seed edge)   */
/* 3. Closed chains become loops immediately                           */
/* 4. Open chains: try to close by matching endpoints at same position */
/*    (handles duplicate vertices from marching cubes / Poisson)       */
/* ------------------------------------------------------------------ */
typedef struct {
    int32_t *verts;    /* vertex indices, last == first (closed) */
    size_t   len;      /* including the closing vertex */
} BoundaryLoop;

/* Helper: find unused edge (cur,nb) in sorted edge list, mark used.
 * Returns edge index or (size_t)-1 if not found/already used. */
static size_t find_and_mark_edge(int32_t cur, int32_t nb,
                                 const UEdge *edges, size_t n_edges,
                                 uint8_t *edge_used)
{
    UEdge key;
    key.v0 = (cur < nb) ? cur : nb;
    key.v1 = (cur < nb) ? nb : cur;
    UEdge *found = (UEdge *)bsearch(&key, edges, n_edges,
                                    sizeof(UEdge), cmp_uedge);
    if (found) {
        size_t fidx = (size_t)(found - edges);
        if (!edge_used[fidx]) {
            edge_used[fidx] = 1;
            return fidx;
        }
    }
    return (size_t)-1;
}

/* Walk forward from cur (coming from prev) along unused boundary edges.
 * Appends vertices to buf starting at buf[*len].
 * Returns the last vertex reached (== dead-end or == some earlier vertex). */
static int32_t walk_chain_forward(int32_t cur, int32_t prev,
                                  const int32_t *adj_off,
                                  const int32_t *adj_list,
                                  const UEdge *edges, size_t n_edges,
                                  uint8_t *edge_used,
                                  int32_t *buf, size_t *len,
                                  size_t max_len)
{
    while (*len < max_len) {
        buf[(*len)++] = cur;

        int32_t next = -1;
        for (int32_t j = adj_off[cur]; j < adj_off[cur + 1]; j++) {
            int32_t nb = adj_list[j];
            if (nb == prev) continue;
            if (find_and_mark_edge(cur, nb, edges, n_edges, edge_used)
                != (size_t)-1) {
                next = nb;
                break;
            }
        }

        if (next < 0) return cur; /* dead end */
        prev = cur;
        cur = next;
    }
    return cur;
}

static size_t chain_into_loops(Arena_T arena,
                               const UEdge *edges, size_t n_edges,
                               size_t nv, const float *verts,
                               BoundaryLoop **out_loops)
{
    /* Build adjacency: boundary_adj[v] = list of boundary neighbors */
    int32_t *degree = (int32_t *)ARENA_CALLOC(arena, (long)nv,
                                               (long)sizeof(int32_t));
    for (size_t i = 0; i < n_edges; i++) {
        degree[edges[i].v0]++;
        degree[edges[i].v1]++;
    }

    int32_t *adj_off = (int32_t *)ARENA_ALLOC(arena,
                         (long)(nv + 1) * (long)sizeof(int32_t));
    adj_off[0] = 0;
    for (size_t v = 0; v < nv; v++) {
        adj_off[v + 1] = adj_off[v] + degree[v];
    }

    int32_t *adj_list = (int32_t *)ARENA_ALLOC(arena,
                          (long)(n_edges * 2) * (long)sizeof(int32_t));
    int32_t *cursor = (int32_t *)ARENA_ALLOC(arena,
                        (long)nv * (long)sizeof(int32_t));
    memcpy(cursor, adj_off, nv * sizeof(int32_t));

    for (size_t i = 0; i < n_edges; i++) {
        int32_t a = edges[i].v0;
        int32_t b = edges[i].v1;
        adj_list[cursor[a]++] = b;
        adj_list[cursor[b]++] = a;
    }

    uint8_t *edge_used = (uint8_t *)ARENA_CALLOC(arena, (long)n_edges,
                                                   (long)sizeof(uint8_t));

    /* Phase 1: walk maximal chains (both directions from each seed edge).
     * Collect closed loops and open chains separately. */
    size_t max_chains = n_edges;
    BoundaryLoop *loops = (BoundaryLoop *)ARENA_ALLOC(arena,
                            (long)max_chains * (long)sizeof(BoundaryLoop));
    size_t n_loops = 0;

    /* Open chains stored temporarily — head_v, tail_v, buf, len */
    typedef struct {
        int32_t  head;   /* first vertex in chain */
        int32_t  tail;   /* last vertex in chain */
        int32_t *buf;
        size_t   len;
    } OpenChain;
    OpenChain *open_chains = (OpenChain *)ARENA_ALLOC(arena,
                               (long)max_chains * (long)sizeof(OpenChain));
    size_t n_open = 0;

    /* Chain-walk scratch: allocate ONCE and reuse every iteration. These are
     * pure scratch (each chain is copied into its own loop_buf/chain_buf below),
     * so a single pair of buffers is reused. Allocating fwd+bwd (each n_edges
     * long) PER seed edge leaked O(n_edges) per chain -> tens of GB and an OOM
     * (hole_fill.c:293) on the full-grid weld's millions of boundary edges. */
    int32_t *fwd = (int32_t *)ARENA_ALLOC(arena,
                      (long)(n_edges + 2) * (long)sizeof(int32_t));
    int32_t *bwd = (int32_t *)ARENA_ALLOC(arena,
                      (long)(n_edges + 2) * (long)sizeof(int32_t));

    for (size_t ei = 0; ei < n_edges; ei++) {
        if (edge_used[ei]) continue;

        int32_t v0 = edges[ei].v0;
        int32_t v1 = edges[ei].v1;
        edge_used[ei] = 1;

        /* Walk forward from v1 (coming from v0). fwd is reused scratch. */
        size_t fwd_len = 0;
        int32_t fwd_end = walk_chain_forward(v1, v0, adj_off, adj_list,
                                              edges, n_edges, edge_used,
                                              fwd, &fwd_len, n_edges);

        if (fwd_end == v0) {
            /* Closed loop. walk_chain_forward already appended the closing v0 as
             * fwd's LAST element, so [v0, fwd...] = [v0, v1, ..., v0] is the
             * closed cycle (first == last). Do NOT append v0 a second time:
             * that put v0 at positions 0, fwd_len AND fwd_len+1, so the
             * duplicate-vertex guard saw verts[0]==verts[loop_n-1] and wrongly
             * rejected every clean Phase-1-closed loop (the BPA path's holes). */
            int32_t *loop_buf = (int32_t *)ARENA_ALLOC(arena,
                                  (long)(fwd_len + 1) * (long)sizeof(int32_t));
            loop_buf[0] = v0;
            memcpy(loop_buf + 1, fwd, fwd_len * sizeof(int32_t));
            loops[n_loops].verts = loop_buf;
            loops[n_loops].len = fwd_len + 1;
            n_loops++;
            continue;
        }

        /* Walk backward from v0 (coming from v1). bwd is reused scratch. */
        size_t bwd_len = 0;
        walk_chain_forward(v0, v1, adj_off, adj_list,
                           edges, n_edges, edge_used,
                           bwd, &bwd_len, n_edges);

        /* Combine: reverse(bwd) + v0's edge neighbor(v1) side(fwd)
         * Full chain: bwd[bwd_len-1], ..., bwd[0], v0_seed_edge_v1_side...
         * Actually bwd starts from v0 and walks away from v1.
         * So bwd = [v0's-other-neighbor, ..., dead-end]
         * Wait — walk_chain_forward starts by appending cur to buf.
         * bwd[0] = v0 (the start vertex passed as cur).
         * No — v0 is 'cur' param and gets appended first. Then it walks
         * to v0's other neighbor (not v1).
         *
         * Actually: walk_chain_forward(v0, v1, ...) starts cur=v0, prev=v1.
         * buf[0] = v0, then finds next != v1, walks that way.
         * So bwd = [v0, nb_of_v0, ..., dead-end].
         *
         * Full chain = reverse(bwd) + fwd
         *   = [dead-end, ..., nb_of_v0, v0] + [v1, ..., fwd_end]
         *   = [dead-end, ..., v0, v1, ..., fwd_end]
         */
        size_t total = bwd_len + fwd_len;
        int32_t *chain_buf = (int32_t *)ARENA_ALLOC(arena,
                                (long)(total + 2) * (long)sizeof(int32_t));
        /* Reverse bwd into chain_buf */
        for (size_t i = 0; i < bwd_len; i++) {
            chain_buf[i] = bwd[bwd_len - 1 - i];
        }
        /* Append fwd */
        memcpy(chain_buf + bwd_len, fwd, fwd_len * sizeof(int32_t));

        int32_t head = chain_buf[0];
        int32_t tail = chain_buf[total - 1];

        /* Check if head == tail (shouldn't happen but be safe) */
        if (head == tail && total >= 3) {
            chain_buf[total] = head;
            loops[n_loops].verts = chain_buf;
            loops[n_loops].len = total + 1;
            n_loops++;
            continue;
        }

        open_chains[n_open].head = head;
        open_chains[n_open].tail = tail;
        open_chains[n_open].buf = chain_buf;
        open_chains[n_open].len = total;
        n_open++;
    }

    /* Phase 2: try to close open chains.
     * A chain can be closed if head and tail are at nearly the same 3D
     * position.  Tolerance is generous (2 voxels) because non-manifold
     * edges from Poisson / density trim can split a single boundary at
     * a junction, leaving endpoints ~1 voxel apart. */
    const float CLOSE_DIST_SQ = 4.0f; /* 2.0 voxels squared */
    for (size_t ci = 0; ci < n_open; ci++) {
        if (open_chains[ci].len == 0) continue; /* already merged */

        int32_t head = open_chains[ci].head;
        int32_t tail = open_chains[ci].tail;
        size_t  clen = open_chains[ci].len;

        /* Check self-close: head and tail at same position */
        float dz = verts[head * 3 + 0] - verts[tail * 3 + 0];
        float dy = verts[head * 3 + 1] - verts[tail * 3 + 1];
        float dx = verts[head * 3 + 2] - verts[tail * 3 + 2];
        if (dz * dz + dy * dy + dx * dx < CLOSE_DIST_SQ && clen >= 3) {
            open_chains[ci].buf[clen] = head; /* close */
            loops[n_loops].verts = open_chains[ci].buf;
            loops[n_loops].len = clen + 1;
            n_loops++;
            open_chains[ci].len = 0; /* consumed */
            continue;
        }

        /* Try to merge with another open chain at matching endpoints */
        for (size_t cj = ci + 1; cj < n_open; cj++) {
            if (open_chains[cj].len == 0) continue;

            int32_t h2 = open_chains[cj].head;
            int32_t t2 = open_chains[cj].tail;
            size_t  l2 = open_chains[cj].len;

            /* Four possible joins: tail-i to head-j, tail-i to tail-j,
             * head-i to head-j, head-i to tail-j.
             * Match by position (epsilon). */
            int join = 0; /* 0=none, 1=tail-head, 2=tail-tail,
                             3=head-head, 4=head-tail */
            float d2;

            /* tail_i == head_j ? */
            dz = verts[tail*3+0] - verts[h2*3+0];
            dy = verts[tail*3+1] - verts[h2*3+1];
            dx = verts[tail*3+2] - verts[h2*3+2];
            d2 = dz*dz + dy*dy + dx*dx;
            if (d2 < CLOSE_DIST_SQ) { join = 1; }

            if (!join) {
                /* tail_i == tail_j ? (need to reverse j) */
                dz = verts[tail*3+0] - verts[t2*3+0];
                dy = verts[tail*3+1] - verts[t2*3+1];
                dx = verts[tail*3+2] - verts[t2*3+2];
                d2 = dz*dz + dy*dy + dx*dx;
                if (d2 < CLOSE_DIST_SQ) { join = 2; }
            }
            if (!join) {
                /* head_i == head_j ? (need to reverse i) */
                dz = verts[head*3+0] - verts[h2*3+0];
                dy = verts[head*3+1] - verts[h2*3+1];
                dx = verts[head*3+2] - verts[h2*3+2];
                d2 = dz*dz + dy*dy + dx*dx;
                if (d2 < CLOSE_DIST_SQ) { join = 3; }
            }
            if (!join) {
                /* head_i == tail_j ? */
                dz = verts[head*3+0] - verts[t2*3+0];
                dy = verts[head*3+1] - verts[t2*3+1];
                dx = verts[head*3+2] - verts[t2*3+2];
                d2 = dz*dz + dy*dy + dx*dx;
                if (d2 < CLOSE_DIST_SQ) { join = 4; }
            }

            if (!join) continue;

            /* Merge chains i and j into one chain */
            size_t merged_len = clen + l2;
            int32_t *merged = (int32_t *)ARENA_ALLOC(arena,
                                (long)(merged_len + 2) * (long)sizeof(int32_t));

            if (join == 1) {
                /* chain_i + chain_j */
                memcpy(merged, open_chains[ci].buf, clen * sizeof(int32_t));
                memcpy(merged + clen, open_chains[cj].buf,
                       l2 * sizeof(int32_t));
            } else if (join == 2) {
                /* chain_i + reverse(chain_j) */
                memcpy(merged, open_chains[ci].buf, clen * sizeof(int32_t));
                for (size_t k = 0; k < l2; k++)
                    merged[clen + k] = open_chains[cj].buf[l2 - 1 - k];
            } else if (join == 3) {
                /* reverse(chain_i) + chain_j */
                for (size_t k = 0; k < clen; k++)
                    merged[k] = open_chains[ci].buf[clen - 1 - k];
                memcpy(merged + clen, open_chains[cj].buf,
                       l2 * sizeof(int32_t));
            } else { /* join == 4 */
                /* chain_j + chain_i */
                memcpy(merged, open_chains[cj].buf, l2 * sizeof(int32_t));
                memcpy(merged + l2, open_chains[ci].buf,
                       clen * sizeof(int32_t));
            }

            /* Update chain i with merged result, consume chain j */
            open_chains[ci].buf = merged;
            open_chains[ci].len = merged_len;
            open_chains[ci].head = merged[0];
            open_chains[ci].tail = merged[merged_len - 1];
            open_chains[cj].len = 0; /* consumed */

            /* Update local vars for continued matching */
            head = open_chains[ci].head;
            tail = open_chains[ci].tail;
            clen = merged_len;

            /* Check if the merged chain is now self-closing */
            dz = verts[head*3+0] - verts[tail*3+0];
            dy = verts[head*3+1] - verts[tail*3+1];
            dx = verts[head*3+2] - verts[tail*3+2];
            if (dz*dz + dy*dy + dx*dx < CLOSE_DIST_SQ && clen >= 3) {
                merged[clen] = head;
                loops[n_loops].verts = merged;
                loops[n_loops].len = clen + 1;
                n_loops++;
                open_chains[ci].len = 0;
                break;
            }
        }
    }

    /* Report any remaining open chains */
    {
        size_t n_remaining = 0;
        for (size_t ci = 0; ci < n_open; ci++) {
            if (open_chains[ci].len > 0) {
                n_remaining++;
                if (n_remaining <= 3) {
                    fprintf(stderr, "  [hole_fill] open chain %zu: %zu verts, "
                            "head=v%d (%.1f,%.1f,%.1f) tail=v%d (%.1f,%.1f,%.1f)\n",
                            ci, open_chains[ci].len,
                            open_chains[ci].head,
                            (double)verts[open_chains[ci].head * 3 + 0],
                            (double)verts[open_chains[ci].head * 3 + 1],
                            (double)verts[open_chains[ci].head * 3 + 2],
                            open_chains[ci].tail,
                            (double)verts[open_chains[ci].tail * 3 + 0],
                            (double)verts[open_chains[ci].tail * 3 + 1],
                            (double)verts[open_chains[ci].tail * 3 + 2]);
                }
            }
        }
        if (n_remaining > 0) {
            fprintf(stderr, "  [hole_fill] %zu open chains could not be closed\n",
                    n_remaining);
        }
    }

    *out_loops = loops;
    return n_loops;
}

/* ------------------------------------------------------------------ */
/* is_interior_hole — check if all verts are far from cube boundary    */
/* ------------------------------------------------------------------ */
static int is_interior_hole(const int32_t *loop_verts, size_t n_loop_verts,
                            const float *verts,
                            const int cube_shape[3], float margin)
{
    for (size_t i = 0; i < n_loop_verts; i++) {
        int32_t vi = loop_verts[i];
        for (int k = 0; k < 3; k++) {
            float coord = verts[vi * 3 + k];
            if (coord <= margin || coord >= (float)cube_shape[k] - 1.0f - margin) {
                return 0;
            }
        }
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* is_simple_polygon — O(n²) segment intersection check               */
/* ------------------------------------------------------------------ */
static int segments_intersect(double ax, double ay, double bx, double by,
                              double cx, double cy, double dx, double dy)
{
    double d1 = (dx - cx) * (ay - cy) - (dy - cy) * (ax - cx);
    double d2 = (dx - cx) * (by - cy) - (dy - cy) * (bx - cx);
    double d3 = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    double d4 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax);

    if (((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) &&
        ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))) {
        return 1;
    }

    return 0;
}

static int is_simple_polygon(const double *pts_2d, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        size_t i2 = (i + 1) % n;
        for (size_t j = i + 2; j < n; j++) {
            if (i == 0 && j == n - 1) continue; /* adjacent */
            size_t j2 = (j + 1) % n;
            if (segments_intersect(
                    pts_2d[i * 2], pts_2d[i * 2 + 1],
                    pts_2d[i2 * 2], pts_2d[i2 * 2 + 1],
                    pts_2d[j * 2], pts_2d[j * 2 + 1],
                    pts_2d[j2 * 2], pts_2d[j2 * 2 + 1])) {
                return 0;
            }
        }
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* pca_project_2d — project boundary loop to 2D via PCA               */
/* Returns u_axis, v_axis, centroid for back-projection.               */
/* ------------------------------------------------------------------ */
static int pca_project_2d(Arena_T arena,
                          const float *verts, const int32_t *loop_verts,
                          size_t n_loop,
                          double *out_2d,
                          float out_u[3], float out_v[3],
                          float out_centroid[3])
{
    /* Gather 3D points for this loop */
    float *pts = (float *)ARENA_ALLOC(arena,
                   (long)n_loop * 3 * (long)sizeof(float));
    for (size_t i = 0; i < n_loop; i++) {
        pts[i * 3 + 0] = verts[loop_verts[i] * 3 + 0];
        pts[i * 3 + 1] = verts[loop_verts[i] * 3 + 1];
        pts[i * 3 + 2] = verts[loop_verts[i] * 3 + 2];
    }

    /* PCA normal + orthonormal basis (float API — outputs still float) */
    float normal[3] = {0.0f, 0.0f, 0.0f};
    int rc = PCA_normal(pts, n_loop, normal, out_centroid);
    if (rc != 0) return -1;

    PCA_orthonormal_basis(normal, out_u, out_v);

    /* Compute centroid in double — PCA_normal truncates to float internally,
     * so we recompute here to avoid losing precision before the projection. */
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (size_t i = 0; i < n_loop; i++) {
        cx += (double)pts[i * 3 + 0];
        cy += (double)pts[i * 3 + 1];
        cz += (double)pts[i * 3 + 2];
    }
    double inv_n = 1.0 / (double)n_loop;
    cx *= inv_n;
    cy *= inv_n;
    cz *= inv_n;

    /* Promote axes to double once */
    double du[3] = { (double)out_u[0], (double)out_u[1], (double)out_u[2] };
    double dv[3] = { (double)out_v[0], (double)out_v[1], (double)out_v[2] };

    /* Project to 2D — all arithmetic in double to avoid false
     * self-intersections from float rounding in the dot products. */
    for (size_t i = 0; i < n_loop; i++) {
        double dx = (double)pts[i * 3 + 0] - cx;
        double dy = (double)pts[i * 3 + 1] - cy;
        double dz = (double)pts[i * 3 + 2] - cz;
        out_2d[i * 2 + 0] = dx * du[0] + dy * du[1] + dz * du[2];
        out_2d[i * 2 + 1] = dx * dv[0] + dy * dv[1] + dz * dv[2];
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/* is_valid_for_cdt — check polygon has no zero-length edges or        */
/* near-duplicate non-adjacent vertices that crash Triangle.           */
/* ------------------------------------------------------------------ */
static int is_valid_for_cdt(const double *pts_2d, size_t n)
{
    double min_dist = 1e-6;
    double min_dist_sq = min_dist * min_dist;

    /* Check for zero-length edges (adjacent vertices at same position).
     * This WILL crash Triangle's CDT — must reject. */
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        double dx = pts_2d[j * 2 + 0] - pts_2d[i * 2 + 0];
        double dy = pts_2d[j * 2 + 1] - pts_2d[i * 2 + 1];
        double d2 = dx * dx + dy * dy;
        if (d2 < min_dist_sq) return 0;
    }

    /* Near-duplicate non-adjacent vertices are OK for CDT as long as
     * the polygon is simple (no edge crossings).  The simplicity check
     * (is_simple_polygon) is the real guard.  Rejecting here was overly
     * conservative and forced polygons through the lossy Clipper2 path,
     * which loses boundary vertices. */

    return 1;
}

/* ------------------------------------------------------------------ */
/* safe_triangulate — crash-guarded wrapper around Triangle library     */
/*                                                                      */
/* Triangle can crash (SIGSEGV / access violation) or call exit() via  */
/* triexit() on degenerate inputs.  This wrapper catches both:          */
/*   - Windows: SEH __try/__except for access violations                */
/*   - Linux:   sigsetjmp + SIGSEGV handler                             */
/*   - Both:    triexit() longjmps to triangle_jmpbuf instead of exit() */
/* Returns 0 on success, -1 on any crash/error.                         */
/* ------------------------------------------------------------------ */

/* Per-hole triangulation timeout (ms). Triangle's exact-arithmetic robustness
 * can fail NON-deterministically on a near-degenerate hole polygon and spin
 * forever -- the 100-cube weld hung here, and the hang even MOVED between runs
 * as malloc addresses shifted (so no flag combination is a guaranteed fix, and
 * -S does NOT bound the segment-recovery loop). We therefore BOUND the call:
 * if triangulate() runs longer than this, abandon the hole (it stays unfilled,
 * which is fine for a tiny interior hole) and carry on. Override via env. */
static int tri_timeout_ms(void)
{
    const char *e = sf_env("SEAM_TRI_TIMEOUT_MS");
    int ms = e ? atoi(e) : 10000;   /* 10 s: a legit large hole CDT is < ~3.5 s */
    if (ms < 100) ms = 100;
    return ms;
}

#ifdef _MSC_VER
#include <process.h>   /* _beginthreadex */

/* SEH must be in its own function (MSVC forbids mixing __try with setjmp). */
static int triangulate_seh(char *flags, struct triangulateio *in,
                           struct triangulateio *out)
{
    __try {
        triangulate(flags, in, out, NULL);
    } __except(1 /* EXCEPTION_EXECUTE_HANDLER */) {
        HFLOG("      [safe_triangulate] SEH exception caught\n");
        return -1;
    }
    return 0;
}

typedef struct {
    char *flags;
    struct triangulateio *in;
    struct triangulateio *out;
    int rc;
} TriJob;

static unsigned __stdcall tri_thread_fn(void *p)
{
    TriJob *j = (TriJob *)p;
    /* triexit() longjmp guard, armed in THIS worker thread (triangle_jmpbuf is a
     * global; holefill runs one triangulate at a time so single-arming is safe).
     * setjmp here + __try inside triangulate_seh are in separate functions, as
     * MSVC requires. */
    triangle_jmpbuf_set = 1;
    if (setjmp(triangle_jmpbuf) != 0) {
        triangle_jmpbuf_set = 0;
        j->rc = -1;
        return 0;
    }
    j->rc = triangulate_seh(j->flags, j->in, j->out);
    triangle_jmpbuf_set = 0;
    return 0;
}

/* Windows watchdog: run triangulate() on a worker thread, wait with a timeout.
 * On timeout we do NOT TerminateThread (Triangle may hold the CRT heap lock mid-
 * malloc -> deadlock). Instead drop the worker to IDLE priority and detach it:
 * it can no longer starve the main weld thread, the OS reclaims it at process
 * exit, and we report the hole unfilled. */
static int safe_triangulate(char *flags, struct triangulateio *in,
                            struct triangulateio *out)
{
    TriJob job;
    job.flags = flags; job.in = in; job.out = out; job.rc = -1;
    uintptr_t h = _beginthreadex(NULL, 0, tri_thread_fn, &job, 0, NULL);
    if (h == 0) {
        /* Spawn failed: run inline (no timeout protection, but better than not
         * triangulating at all). */
        tri_thread_fn(&job);
        return job.rc;
    }
    DWORD ms = (DWORD)tri_timeout_ms();
    DWORD w = WaitForSingleObject((HANDLE)h, ms);
    if (w != WAIT_OBJECT_0) {
        /* Kill the worker. It must NOT survive: Triangle's globals
         * (triangle_jmpbuf, randomseed) are unsafe for two concurrent
         * triangulations, so a detached-but-still-spinning zombie races the
         * next hole's call and crashes the run (it did, at loop 4350).
         * TerminateThread is safe here: the runaway is a geometric flip /
         * segment-recovery cycle that does not grow memory (no OOM across 200 s
         * of spinning), so it is not killed mid-malloc holding the CRT heap
         * lock. Its in-progress Triangle pools leak (a few KB) -- fine for the
         * rare abandoned hole. */
        TerminateThread((HANDLE)h, 1);
        WaitForSingleObject((HANDLE)h, 2000);   /* reap before next call */
        CloseHandle((HANDLE)h);
        triangle_jmpbuf_set = 0;
        fprintf(stderr,
            "      [safe_triangulate] *** TIMEOUT after %lu ms -- killed worker, abandoning hole ***\n",
            (unsigned long)ms);
        HFFLUSH();
        return -1;
    }
    CloseHandle((HANDLE)h);
    return job.rc;
}
#else
/* Linux: SIGSEGV (crash) AND SIGALRM (watchdog timeout) both siglongjmp out of
 * triangle_segv_handler. The setjmp(triangle_jmpbuf) catches triexit(). */
static int safe_triangulate(char *flags, struct triangulateio *in,
                            struct triangulateio *out)
{
    triangle_jmpbuf_set = 1;
    int jrc = setjmp(triangle_jmpbuf);
    if (jrc != 0) {
        triangle_jmpbuf_set = 0;
        HFLOG("      [safe_triangulate] triexit caught (code=%d)\n", jrc);
        return -1;
    }

    struct sigaction sa_new, sa_old_segv, sa_old_alrm;
    memset(&sa_new, 0, sizeof(sa_new));
    sa_new.sa_handler = triangle_segv_handler;
    sigemptyset(&sa_new.sa_mask);
    sa_new.sa_flags = 0;
    sigaction(SIGSEGV, &sa_new, &sa_old_segv);
    sigaction(SIGALRM, &sa_new, &sa_old_alrm);

    triangle_segv_guard = 1;
    int rc;
    if (sigsetjmp(triangle_segv_jmpbuf, 1) != 0) {
        /* Returned from SIGSEGV or SIGALRM (timeout) handler */
        triangle_segv_guard = 0;
        triangle_jmpbuf_set = 0;
        alarm(0);
        sigaction(SIGSEGV, &sa_old_segv, NULL);
        sigaction(SIGALRM, &sa_old_alrm, NULL);
        HFLOG("      [safe_triangulate] *** SIGSEGV/TIMEOUT caught -- abandoning hole ***\n");
        return -1;
    }

    {
        int ms = tri_timeout_ms();
        alarm((unsigned)((ms + 999) / 1000));   /* seconds, rounded up */
    }
    triangulate(flags, in, out, NULL);
    rc = 0;
    alarm(0);

    triangle_segv_guard = 0;
    sigaction(SIGSEGV, &sa_old_segv, NULL);
    sigaction(SIGALRM, &sa_old_alrm, NULL);

    triangle_jmpbuf_set = 0;
    return rc;
}
#endif

/* ------------------------------------------------------------------ */
/* triangulate_polygon — CDT via Triangle library                      */
/* Returns triangulated vertices and faces.                            */
/* ------------------------------------------------------------------ */
static int triangulate_polygon(const double *pts_2d, size_t n,
                               double **out_verts_2d, size_t *out_nv,
                               int32_t **out_faces, size_t *out_nf,
                               int skip_quality)
{
    if (n < 3) return -1;

    /* Validate input before calling Triangle */
    if (!is_valid_for_cdt(pts_2d, n)) {
        return -1;
    }

    struct triangulateio in_tri;
    struct triangulateio out_tri;

    memset(&in_tri, 0, sizeof(in_tri));
    memset(&out_tri, 0, sizeof(out_tri));

    /* Input: vertices */
    in_tri.numberofpoints = (int)n;
    in_tri.pointlist = (REAL *)malloc(n * 2 * sizeof(REAL));
    if (!in_tri.pointlist) return -1;
    memcpy(in_tri.pointlist, pts_2d, n * 2 * sizeof(REAL));

    /* Input: segments (closed polygon boundary) */
    in_tri.numberofsegments = (int)n;
    in_tri.segmentlist = (int *)malloc(n * 2 * sizeof(int));
    if (!in_tri.segmentlist) {
        free(in_tri.pointlist);
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        in_tri.segmentlist[i * 2 + 0] = (int)i;
        in_tri.segmentlist[i * 2 + 1] = (int)((i + 1) % n);
    }

    /* Compute edge length statistics */
    double total_len = 0.0;
    double min_edge_len = 1e30;
    double max_edge_len = 0.0;
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        double dx = pts_2d[j * 2 + 0] - pts_2d[i * 2 + 0];
        double dy = pts_2d[j * 2 + 1] - pts_2d[i * 2 + 1];
        double len = sqrt(dx * dx + dy * dy);
        total_len += len;
        if (len > 0 && len < min_edge_len) min_edge_len = len;
        if (len > max_edge_len) max_edge_len = len;
    }
    double avg_edge = total_len / (double)n;

    /* Auto-detect degenerate aspect ratio: if max/min edge ratio is extreme,
     * quality refinement will try to insert Steiner points forever or crash */
    double edge_ratio = (min_edge_len > 0) ? max_edge_len / min_edge_len : 1e6;
    if (edge_ratio > 50.0) {
        HFLOG("      [triangulate_polygon] extreme edge ratio %.1f "
                "(min=%.2e max=%.2e), forcing skip_quality\n",
                edge_ratio, min_edge_len, max_edge_len);
        skip_quality = 1;
    }

    /* Flags: p=PSLG, Y=no Steiner on boundary, D=Delaunay, z=zero-indexed,
     * Q=quiet.  When skip_quality is off, also add q20=min angle 20° and
     * a=area constraint to get nice triangles for CG smoothing. */
    char flags[64] = {0};
    if (skip_quality) {
        /* No quality refinement -> constrained Delaunay only. NOTE: this still is
         * not hang-proof -- Triangle's segment recovery can loop on degenerate
         * input regardless of flags (it ignores -S). The real guard is the
         * watchdog timeout in safe_triangulate(); the flags only reduce how often
         * it trips. */
        snprintf(flags, sizeof(flags), "pYDzQ");
    } else {
        /* Quality refinement (q20 + area). Triangle's refinement can loop FOREVER
         * inserting Steiner points at a sub-20-degree input corner -- this is what
         * hung the 100-cube weld (a Clipper2-output polygon, loop 2560). The -S
         * switch caps added Steiner points so refinement ALWAYS terminates: on
         * hitting the cap Triangle returns the best mesh so far -- still a valid,
         * conforming triangulation, just with some triangles missing the angle
         * bound, which is fine for a hole fill we CG-smooth afterwards. Cap scales
         * with loop size and is clamped to a sane ceiling (50k Steiner < ~50 ms;
         * the runaway case would otherwise insert millions). */
        double max_area = 0.5 * avg_edge * avg_edge;
        int steiner_cap = (int)(200 * n + 2000);
        if (steiner_cap > 50000) steiner_cap = 50000;
        snprintf(flags, sizeof(flags), "pq20YDa%.6fzQS%d", max_area, steiner_cap);
    }

    HFLOG("      [triangulate_polygon] n=%zu avg_edge=%.4f skip_quality=%d flags='%s'\n",
            n, avg_edge, skip_quality, flags);
    HFFLUSH();

    /* Run triangulation (crash-guarded) */
    int tri_rc = safe_triangulate(flags, &in_tri, &out_tri);
    if (tri_rc != 0) {
        HFLOG("      [triangulate_polygon] Triangle crashed — skipping hole\n");
        free(in_tri.pointlist);
        free(in_tri.segmentlist);
        return -1;
    }

    HFLOG("      [triangulate_polygon] done: %d pts, %d tris\n",
            out_tri.numberofpoints, out_tri.numberoftriangles);
    HFFLUSH();

    if (out_tri.numberoftriangles == 0) {
        free(in_tri.pointlist);
        free(in_tri.segmentlist);
        if (out_tri.pointlist) trifree((VOID *)out_tri.pointlist);
        if (out_tri.pointmarkerlist) trifree((VOID *)out_tri.pointmarkerlist);
        if (out_tri.segmentlist) trifree((VOID *)out_tri.segmentlist);
        if (out_tri.segmentmarkerlist) trifree((VOID *)out_tri.segmentmarkerlist);
        return -1;
    }

    /* Copy output */
    size_t out_n = (size_t)out_tri.numberofpoints;
    size_t out_f = (size_t)out_tri.numberoftriangles;

    *out_verts_2d = (double *)malloc(out_n * 2 * sizeof(double));
    *out_faces = (int32_t *)malloc(out_f * 3 * sizeof(int32_t));

    if (!*out_verts_2d || !*out_faces) {
        free(*out_verts_2d);
        free(*out_faces);
        free(in_tri.pointlist);
        free(in_tri.segmentlist);
        trifree((VOID *)out_tri.pointlist);
        trifree((VOID *)out_tri.trianglelist);
        return -1;
    }

    memcpy(*out_verts_2d, out_tri.pointlist, out_n * 2 * sizeof(double));
    for (size_t f = 0; f < out_f; f++) {
        (*out_faces)[f * 3 + 0] = (int32_t)out_tri.trianglelist[f * 3 + 0];
        (*out_faces)[f * 3 + 1] = (int32_t)out_tri.trianglelist[f * 3 + 1];
        (*out_faces)[f * 3 + 2] = (int32_t)out_tri.trianglelist[f * 3 + 2];
    }

    *out_nv = out_n;
    *out_nf = out_f;

    /* Free Triangle's memory */
    free(in_tri.pointlist);
    free(in_tri.segmentlist);
    trifree((VOID *)out_tri.pointlist);
    trifree((VOID *)out_tri.trianglelist);
    if (out_tri.pointmarkerlist) trifree((VOID *)out_tri.pointmarkerlist);
    if (out_tri.segmentlist) trifree((VOID *)out_tri.segmentlist);
    if (out_tri.segmentmarkerlist) trifree((VOID *)out_tri.segmentmarkerlist);

    return 0;
}

/* ------------------------------------------------------------------ */
/* embed_steiner_3d — place Steiner points on the best-fit plane       */
/* ------------------------------------------------------------------ */
static void embed_steiner_3d(const double *pts_2d_tri, size_t total_nv,
                             const float *boundary_3d, size_t n_boundary,
                             const float u_axis[3], const float v_axis[3],
                             const float centroid[3],
                             float *out_3d)
{
    /* Boundary vertices: copy from boundary_3d */
    memcpy(out_3d, boundary_3d, n_boundary * 3 * sizeof(float));

    /* Steiner points: project from 2D to 3D via PCA basis.
     * Compute in double, truncate to float only at the end. */
    for (size_t i = n_boundary; i < total_nv; i++) {
        double u = pts_2d_tri[i * 2 + 0];
        double v = pts_2d_tri[i * 2 + 1];
        out_3d[i * 3 + 0] = (float)((double)centroid[0] + u * (double)u_axis[0] + v * (double)v_axis[0]);
        out_3d[i * 3 + 1] = (float)((double)centroid[1] + u * (double)u_axis[1] + v * (double)v_axis[1]);
        out_3d[i * 3 + 2] = (float)((double)centroid[2] + u * (double)u_axis[2] + v * (double)v_axis[2]);
    }
}

/* ------------------------------------------------------------------ */
/* CG solver for cotangent Laplacian — all double internals            */
/* Solve L @ x = b with Jacobi preconditioner.                         */
/* ------------------------------------------------------------------ */
static void csr_matvec_double(const int32_t *off, const int32_t *tgt,
                              const float *weight, int32_t nrows,
                              const double *x, double *y)
{
    for (int32_t i = 0; i < nrows; i++) {
        double sum = 0.0;
        for (int32_t j = off[i]; j < off[i + 1]; j++) {
            sum += (double)weight[j] * x[tgt[j]];
        }
        y[i] = sum;
    }
}

static int pcg_solve(const int32_t *off, const int32_t *tgt,
                     const float *weight, int32_t nrows,
                     const double *diag, const double *b, double *x)
{
    /* Jacobi preconditioner: M^{-1} = 1/diag */
    double *r  = (double *)calloc((size_t)nrows, sizeof(double));
    double *z  = (double *)calloc((size_t)nrows, sizeof(double));
    double *p  = (double *)calloc((size_t)nrows, sizeof(double));
    double *Ap = (double *)calloc((size_t)nrows, sizeof(double));

    if (!r || !z || !p || !Ap) {
        free(r); free(z); free(p); free(Ap);
        return -1;
    }

    /* r = b - A @ x */
    csr_matvec_double(off, tgt, weight, nrows, x, Ap);
    double bnorm = 0.0;
    for (int32_t i = 0; i < nrows; i++) {
        r[i] = b[i] - Ap[i];
        bnorm += b[i] * b[i];
    }
    bnorm = sqrt(bnorm);
    if (bnorm < 1e-15) {
        /* b is zero, x=0 is the solution */
        free(r); free(z); free(p); free(Ap);
        return 0;
    }

    /* z = M^{-1} r, p = z */
    double rz = 0.0;
    for (int32_t i = 0; i < nrows; i++) {
        double d = diag[i];
        z[i] = (fabs(d) > 1e-15) ? r[i] / d : r[i];
        p[i] = z[i];
        rz += r[i] * z[i];
    }

    int converged = 0;
    for (int iter = 0; iter < CG_MAX_ITER; iter++) {
        /* Ap = A @ p */
        csr_matvec_double(off, tgt, weight, nrows, p, Ap);

        /* alpha = rz / (p . Ap) */
        double pAp = 0.0;
        for (int32_t i = 0; i < nrows; i++) pAp += p[i] * Ap[i];
        if (fabs(pAp) < 1e-30) break;
        double alpha = rz / pAp;

        /* x += alpha * p, r -= alpha * Ap */
        double rnorm = 0.0;
        for (int32_t i = 0; i < nrows; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            rnorm += r[i] * r[i];
        }
        rnorm = sqrt(rnorm);

        if (rnorm / bnorm < CG_TOL) {
            converged = 1;
            break;
        }

        /* z = M^{-1} r */
        double rz_new = 0.0;
        for (int32_t i = 0; i < nrows; i++) {
            double d = diag[i];
            z[i] = (fabs(d) > 1e-15) ? r[i] / d : r[i];
            rz_new += r[i] * z[i];
        }

        /* p = z + beta * p */
        double beta = rz_new / (fabs(rz) > 1e-30 ? rz : 1e-30);
        for (int32_t i = 0; i < nrows; i++) {
            p[i] = z[i] + beta * p[i];
        }
        rz = rz_new;
    }

    free(r);
    free(z);
    free(p);
    free(Ap);
    return converged ? 0 : 1; /* 1 = did not converge, but not fatal */
}

/* ------------------------------------------------------------------ */
/* smooth_fill_cotangent — cotangent Laplacian CG solve on fill patch  */
/* Boundary verts = Dirichlet BC. Interior verts are smoothed.         */
/* ------------------------------------------------------------------ */
static int smooth_fill_cotangent(Arena_T arena,
                                 float *fill_verts,
                                 const int32_t *fill_faces, size_t fill_nf,
                                 size_t fill_nv, size_t n_boundary)
{
    if (fill_nv <= n_boundary) return 0; /* no interior points */

    Arena_Mark mark = Arena_save(arena);

    /* Build cotangent weight matrix (full NxN, boundary rows = identity).
     * Use COO triples accumulation. */
    size_t max_nnz = fill_nf * 6 + fill_nv; /* edges + diagonal */

    /* We'll use malloc for COO since the CG solver uses malloc too */
    int32_t *coo_r = (int32_t *)malloc(max_nnz * 2 * sizeof(int32_t));
    int32_t *coo_c = (int32_t *)malloc(max_nnz * 2 * sizeof(int32_t));
    float *coo_v = (float *)malloc(max_nnz * 2 * sizeof(float));
    double *diag_val = (double *)calloc(fill_nv, sizeof(double));

    if (!coo_r || !coo_c || !coo_v || !diag_val) {
        free(coo_r); free(coo_c); free(coo_v); free(diag_val);
        Arena_restore(arena, mark);
        return -1;
    }

    size_t nnz = 0;

    /* For each face, accumulate cotangent weights for each edge */
    for (size_t f = 0; f < fill_nf; f++) {
        int32_t vi[3];
        vi[0] = fill_faces[f * 3 + 0];
        vi[1] = fill_faces[f * 3 + 1];
        vi[2] = fill_faces[f * 3 + 2];

        for (int e = 0; e < 3; e++) {
            int32_t i = vi[e];
            int32_t j = vi[(e + 1) % 3];
            int32_t k = vi[(e + 2) % 3]; /* opposite vertex */

            /* cot(angle at k) for edge (i,j) */
            float ei[3], ej[3];
            for (int d = 0; d < 3; d++) {
                ei[d] = fill_verts[i * 3 + d] - fill_verts[k * 3 + d];
                ej[d] = fill_verts[j * 3 + d] - fill_verts[k * 3 + d];
            }
            float dot = ei[0] * ej[0] + ei[1] * ej[1] + ei[2] * ej[2];
            float cx = ei[1] * ej[2] - ei[2] * ej[1];
            float cy = ei[2] * ej[0] - ei[0] * ej[2];
            float cz = ei[0] * ej[1] - ei[1] * ej[0];
            float sin_val = sqrtf(cx * cx + cy * cy + cz * cz);

            float cot_w = 0.0f;
            if (sin_val > COT_SIN_GUARD) {
                cot_w = dot / sin_val;
            }
            float w = 0.5f * cot_w;
            if (w < COT_WEIGHT_MIN) w = COT_WEIGHT_MIN;

            /* Add off-diagonal entries (i,j) and (j,i) */
            /* Only for interior rows */
            if ((size_t)i >= n_boundary) {
                coo_r[nnz] = i; coo_c[nnz] = j; coo_v[nnz] = w; nnz++;
                diag_val[i] -= (double)w;
            }
            if ((size_t)j >= n_boundary) {
                coo_r[nnz] = j; coo_c[nnz] = i; coo_v[nnz] = w; nnz++;
                diag_val[j] -= (double)w;
            }
        }
    }

    /* Add diagonal entries */
    for (size_t i = 0; i < fill_nv; i++) {
        if (i < n_boundary) {
            /* Boundary row: identity */
            coo_r[nnz] = (int32_t)i;
            coo_c[nnz] = (int32_t)i;
            coo_v[nnz] = 1.0f;
            diag_val[i] = 1.0;
            nnz++;
        } else {
            /* Interior row: negative sum of off-diag */
            coo_r[nnz] = (int32_t)i;
            coo_c[nnz] = (int32_t)i;
            coo_v[nnz] = (float)diag_val[i];
            nnz++;
        }
    }

    /* Build CSR from COO */
    CSR_T L = CSR_from_coo(arena, coo_r, coo_c, coo_v, nnz, fill_nv);

    const int32_t *off = CSR_offset(L);
    const int32_t *tgt = CSR_target(L);
    const float *wt = CSR_weight(L);
    int32_t nrows = CSR_nrows(L);

    /* Solve for each coordinate independently */
    for (size_t coord = 0; coord < 3; coord++) {
        double *b = (double *)calloc(fill_nv, sizeof(double));
        double *x = (double *)calloc(fill_nv, sizeof(double));
        if (!b || !x) {
            free(b); free(x);
            free(coo_r); free(coo_c); free(coo_v); free(diag_val);
            Arena_restore(arena, mark);
            return -1;
        }

        /* RHS: boundary = boundary coord, interior = 0 */
        for (size_t i = 0; i < fill_nv; i++) {
            if (i < n_boundary) {
                b[i] = (double)fill_verts[i * 3 + coord];
                x[i] = b[i];
            } else {
                x[i] = (double)fill_verts[i * 3 + coord];
            }
        }

        pcg_solve(off, tgt, wt, nrows, diag_val, b, x);

        /* Write back (interior only — boundary stays fixed) */
        for (size_t i = n_boundary; i < fill_nv; i++) {
            fill_verts[i * 3 + coord] = (float)x[i];
        }

        free(b);
        free(x);
    }

    free(coo_r);
    free(coo_c);
    free(coo_v);
    free(diag_val);
    Arena_restore(arena, mark);
    return 0;
}

/* ------------------------------------------------------------------ */
/* orient_fill — ensure fill boundary edges oppose mesh boundary edges */
/* The fill's boundary edge (bdry_map[0]→bdry_map[1]) must be the     */
/* OPPOSITE of the mesh's half-edge direction for proper stitching.    */
/* If the mesh has half-edge a→b, the fill must have b→a.              */
/* ------------------------------------------------------------------ */
static void orient_fill(int32_t *fill_faces, size_t fill_nf,
                        const int32_t *mesh_faces, size_t mesh_nf,
                        const int32_t *bdry_map, size_t n_boundary)
{
    if (n_boundary < 2 || fill_nf == 0) return;

    /* Find the actual half-edge direction that fill faces produce for
     * boundary edges, then compare against mesh.
     *
     * We sample multiple boundary edges to use majority vote, since
     * individual edges may have local defects from vertex merging. */
    int n_consistent = 0;  /* fill opposes mesh = correct */
    int n_inconsist = 0;   /* fill matches mesh = wrong, need flip */

    for (size_t bi = 0; bi < n_boundary; bi++) {
        size_t bj = (bi + 1) % n_boundary;
        int32_t mesh_a = bdry_map[bi];
        int32_t mesh_b = bdry_map[bj];

        /* Find the fill face that has local boundary edge {bi, bj}.
         * Determine which direction the fill has: bi→bj or bj→bi. */
        int fill_has_ij = 0;
        int fill_has_ji = 0;
        for (size_t f = 0; f < fill_nf; f++) {
            int32_t v0 = fill_faces[f*3+0];
            int32_t v1 = fill_faces[f*3+1];
            int32_t v2 = fill_faces[f*3+2];
            if ((v0==(int32_t)bi && v1==(int32_t)bj) ||
                (v1==(int32_t)bi && v2==(int32_t)bj) ||
                (v2==(int32_t)bi && v0==(int32_t)bj)) {
                fill_has_ij = 1;
                break;
            }
            if ((v0==(int32_t)bj && v1==(int32_t)bi) ||
                (v1==(int32_t)bj && v2==(int32_t)bi) ||
                (v2==(int32_t)bj && v0==(int32_t)bi)) {
                fill_has_ji = 1;
                break;
            }
        }
        if (!fill_has_ij && !fill_has_ji) continue;

        /* The fill face, after stitch, will map local bi→bj to mesh_a→mesh_b,
         * and local bj→bi to mesh_b→mesh_a.
         * Now check what direction the MESH has for edge {mesh_a, mesh_b}. */
        int mesh_has_ab = 0;
        int mesh_has_ba = 0;
        for (size_t f = 0; f < mesh_nf; f++) {
            int32_t v0 = mesh_faces[f*3+0];
            int32_t v1 = mesh_faces[f*3+1];
            int32_t v2 = mesh_faces[f*3+2];
            if ((v0==mesh_a && v1==mesh_b) ||
                (v1==mesh_a && v2==mesh_b) ||
                (v2==mesh_a && v0==mesh_b)) {
                mesh_has_ab = 1;
                break;
            }
            if ((v0==mesh_b && v1==mesh_a) ||
                (v1==mesh_b && v2==mesh_a) ||
                (v2==mesh_b && v0==mesh_a)) {
                mesh_has_ba = 1;
                break;
            }
        }
        if (!mesh_has_ab && !mesh_has_ba) continue;

        /* After stitch, fill will produce:
         *   fill_has_ij → mesh_a→mesh_b
         *   fill_has_ji → mesh_b→mesh_a
         * For correct manifold, fill must OPPOSE mesh at boundary.
         * So if mesh has a→b, fill should produce b→a (fill_has_ji).
         *    if mesh has b→a, fill should produce a→b (fill_has_ij). */
        if ((mesh_has_ab && fill_has_ji) || (mesh_has_ba && fill_has_ij)) {
            n_consistent++;
        } else {
            n_inconsist++;
        }
    }

    /* Flip all fill faces if majority are inconsistent */
    if (n_inconsist > n_consistent) {
        for (size_t f = 0; f < fill_nf; f++) {
            int32_t tmp = fill_faces[f * 3 + 1];
            fill_faces[f * 3 + 1] = fill_faces[f * 3 + 2];
            fill_faces[f * 3 + 2] = tmp;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Fill result for one hole                                            */
/* ------------------------------------------------------------------ */
typedef struct {
    float   *verts;         /* [nv * 3] */
    int32_t *faces;         /* [nf * 3] */
    size_t   nv;
    size_t   nf;
    size_t   n_boundary;    /* first n_boundary verts are boundary */
    int32_t *boundary_map;  /* [n_boundary]: fill boundary idx -> mesh vertex idx */
    /* Vertex merges from 2D edge collapse: merge_src[i] -> merge_dst[i] */
    int32_t *merge_src;     /* mesh vertex indices to be replaced */
    int32_t *merge_dst;     /* mesh vertex indices to replace with */
    size_t   n_merge;       /* number of merge pairs */
} HoleFillResult;

/* ------------------------------------------------------------------ */
/* prune_degenerate_fill — drop sliver/duplicate triangles from a fill */
/* ------------------------------------------------------------------ */
/* The minimal-surface smoother pins the boundary and relaxes each interior
 * (Steiner) point to the cotangent-weighted average of its neighbours. For a
 * thin slit hole that legitimately drives the interior points into a sub-voxel
 * cluster, so the fill emerges with near-coincident Steiner vertices and needle
 * triangles between them (observed on the Z=4480 weld seam: 8 tris under 0.02
 * vox area, a 0.06-vox edge, an exact-duplicate Steiner pair).
 *
 * Steiner vertices are FILL-PRIVATE until stitch_fills runs, so collapsing them
 * changes nothing outside this patch -- topology, the seam, the boundary loop and
 * bdry_map are all untouched. We union the Steiner endpoints of every triangle
 * edge shorter than eps_vox (boundary verts [0,nb) are never merged or moved --
 * they are shared mesh verts), drop faces that collapse to < 3 distinct verts or
 * below area_eps, and compact the surviving Steiner block. The whole thing is
 * built into scratch and committed ONLY if it stays manifold -- a residual sliver
 * is a lesser evil than a non-manifold edge. Returns 1 if it changed the fill, 0
 * otherwise. Exposed (non-static) for the unit test; no public prototype. */
static int32_t pf_find(int32_t *uf, int32_t x)
{ while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

static int pf_u64cmp(const void *a, const void *b)
{ uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return x < y ? -1 : (x > y ? 1 : 0); }

static double pf_tri_area(const float *V, int32_t a, int32_t b, int32_t c)
{
    double e1x = (double)V[b*3+0]-V[a*3+0], e1y = (double)V[b*3+1]-V[a*3+1], e1z = (double)V[b*3+2]-V[a*3+2];
    double e2x = (double)V[c*3+0]-V[a*3+0], e2y = (double)V[c*3+1]-V[a*3+1], e2z = (double)V[c*3+2]-V[a*3+2];
    double cx = e1y*e2z-e1z*e2y, cy = e1z*e2x-e1x*e2z, cz = e1x*e2y-e1y*e2x;
    return 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

/* every undirected edge incident to <= 2 faces, no repeated-vertex face */
static int pf_is_manifold(const int32_t *faces, size_t nf)
{
    if (nf == 0) return 1;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3+0], b = faces[f*3+1], c = faces[f*3+2];
        if (a == b || b == c || a == c) return 0;
    }
    uint64_t *ed = (uint64_t *)malloc(nf * 3 * sizeof(uint64_t));
    if (!ed) return 1;                       /* can't check -> don't block commit */
    size_t n = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t v[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            uint32_t a = (uint32_t)v[e], b = (uint32_t)v[(e+1)%3];
            uint32_t lo = a < b ? a : b, hi = a < b ? b : a;
            ed[n++] = ((uint64_t)lo << 32) | hi;
        }
    }
    qsort(ed, n, sizeof(uint64_t), pf_u64cmp);
    int ok = 1; size_t i = 0;
    while (i < n) { size_t j = i + 1; while (j < n && ed[j] == ed[i]) j++;
                    if (j - i > 2) { ok = 0; break; } i = j; }
    free(ed);
    return ok;
}

int prune_degenerate_fill(float *verts, int32_t *faces,
                          size_t *p_nv, size_t *p_nf,
                          size_t n_boundary,
                          float eps_vox, float area_eps)
{
    size_t nv = *p_nv, nf = *p_nf;
    if (nv <= n_boundary || nf == 0) return 0;        /* no Steiner to prune */
    const double eps2 = (double)eps_vox * (double)eps_vox;

    int32_t *uf = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!uf) return 0;
    for (size_t i = 0; i < nv; i++) uf[i] = (int32_t)i;

    /* union Steiner endpoints of every short triangle edge (original coords) */
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f*3 + e], b = faces[f*3 + (e+1)%3];
            int32_t ra = pf_find(uf, a), rb = pf_find(uf, b);
            if (ra == rb) continue;
            if ((size_t)ra < n_boundary && (size_t)rb < n_boundary) continue; /* never merge two boundary */
            double dx = (double)verts[a*3+0]-verts[b*3+0];
            double dy = (double)verts[a*3+1]-verts[b*3+1];
            double dz = (double)verts[a*3+2]-verts[b*3+2];
            if (dx*dx + dy*dy + dz*dz >= eps2) continue;
            if      ((size_t)ra < n_boundary) uf[rb] = ra;   /* boundary wins */
            else if ((size_t)rb < n_boundary) uf[ra] = rb;
            else if (ra < rb)                 uf[rb] = ra;
            else                              uf[ra] = rb;
        }
    }

    /* compacted index per old vertex: boundary identity, Steiner reps repacked */
    int32_t *nidx = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!nidx) { free(uf); return 0; }
    size_t next = n_boundary;
    int any = 0;
    for (size_t i = 0; i < n_boundary; i++) nidx[i] = (int32_t)i;
    for (size_t i = n_boundary; i < nv; i++) {
        int32_t r = pf_find(uf, (int32_t)i);
        if (r == (int32_t)i) nidx[i] = (int32_t)next++;      /* surviving rep */
        else { nidx[i] = -1; any = 1; }                      /* merged away */
    }
    if (!any) { free(uf); free(nidx); return 0; }            /* nothing collapsed */
    for (size_t i = n_boundary; i < nv; i++) {
        int32_t r = pf_find(uf, (int32_t)i);
        if (r != (int32_t)i) nidx[i] = ((size_t)r < n_boundary) ? r : nidx[r];
    }
    size_t new_nv = next;

    /* scratch verts (boundary verbatim, Steiner reps compacted) */
    float *sv = (float *)malloc(new_nv * 3 * sizeof(float));
    int32_t *sf = (int32_t *)malloc(nf * 3 * sizeof(int32_t));
    if (!sv || !sf) { free(uf); free(nidx); free(sv); free(sf); return 0; }
    for (size_t i = 0; i < n_boundary; i++)
        for (int k = 0; k < 3; k++) sv[i*3+k] = verts[i*3+k];
    for (size_t i = n_boundary; i < nv; i++)
        if (pf_find(uf, (int32_t)i) == (int32_t)i)
            for (int k = 0; k < 3; k++) sv[nidx[i]*3+k] = verts[i*3+k];

    /* surviving faces, dropping collapsed / zero-area */
    size_t out = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = nidx[faces[f*3+0]], b = nidx[faces[f*3+1]], c = nidx[faces[f*3+2]];
        if (a == b || b == c || a == c) continue;
        if (pf_tri_area(sv, a, b, c) < (double)area_eps) continue;
        sf[out*3+0] = a; sf[out*3+1] = b; sf[out*3+2] = c; out++;
    }

    /* commit only if the pruned fill is still manifold and non-empty */
    int changed = 0;
    if (out > 0 && pf_is_manifold(sf, out)) {
        for (size_t i = 0; i < new_nv * 3; i++) verts[i] = sv[i];
        for (size_t i = 0; i < out * 3; i++)    faces[i] = sf[i];
        *p_nv = new_nv; *p_nf = out;
        changed = 1;
    }
    free(uf); free(nidx); free(sv); free(sf);
    return changed;
}

/* fill_one_hole's manifold pre-check calls this before its definition below */
static int loop_diag_in_mesh(int32_t va, int32_t vb,
                             const int32_t *mesh_faces, size_t mesh_nf);

/* ------------------------------------------------------------------ */
/* fill_one_hole — fill a single boundary loop                         */
/* ------------------------------------------------------------------ */
static int fill_one_hole(Arena_T arena,
                         const float *mesh_verts,
                         const int32_t *mesh_faces, size_t mesh_nf,
                         const int32_t *loop_verts, size_t n_loop,
                         HoleFillResult *result)
{
    Arena_Mark mark = Arena_save(arena);

    /* Gather 3D boundary coordinates and build boundary->mesh vertex map */
    float *bdry_3d = (float *)ARENA_ALLOC(arena,
                       (long)n_loop * 3 * (long)sizeof(float));
    int32_t *bdry_map = (int32_t *)ARENA_ALLOC(arena,
                          (long)n_loop * (long)sizeof(int32_t));
    for (size_t i = 0; i < n_loop; i++) {
        bdry_3d[i * 3 + 0] = mesh_verts[loop_verts[i] * 3 + 0];
        bdry_3d[i * 3 + 1] = mesh_verts[loop_verts[i] * 3 + 1];
        bdry_3d[i * 3 + 2] = mesh_verts[loop_verts[i] * 3 + 2];
        bdry_map[i] = loop_verts[i];
    }

    /* PCA project to 2D */
    double *pts_2d = (double *)ARENA_ALLOC(arena,
                       (long)n_loop * 2 * (long)sizeof(double));
    float u_axis[3] = {0.0f, 0.0f, 0.0f};
    float v_axis[3] = {0.0f, 0.0f, 0.0f};
    float centroid[3] = {0.0f, 0.0f, 0.0f};

    int rc = pca_project_2d(arena, mesh_verts, loop_verts, n_loop,
                            pts_2d, u_axis, v_axis, centroid);
    if (rc != 0) {
        Arena_restore(arena, mark);
        return -1;
    }

    /* Collapse short 2D edges: edges much shorter than the average are
     * degenerate in the 2D projection (vertices differ in 3D but project
     * nearly to the same 2D point).  Remove the shorter vertex from the
     * loop so CDT sees no degenerate edges.  Use a relative threshold
     * (0.1% of avg edge) rather than absolute to handle varying scales. */

    /* Compute average 2D edge length for relative threshold */
    double total_edge_len = 0.0;
    for (size_t i = 0; i < n_loop; i++) {
        size_t j = (i + 1) % n_loop;
        double dx = pts_2d[j * 2] - pts_2d[i * 2];
        double dy = pts_2d[j * 2 + 1] - pts_2d[i * 2 + 1];
        total_edge_len += sqrt(dx * dx + dy * dy);
    }
    double avg_2d_edge = total_edge_len / (double)n_loop;

    HFLOG("    [fill_one_hole] n_loop=%zu avg_2d_edge=%.4f, checking short edges...\n",
            n_loop, avg_2d_edge);
    /* Track vertex merges from collapse: src[i] -> dst[i] in mesh indices.
     * Allocate worst case (n_loop entries). */
    int32_t *collapse_src = (int32_t *)ARENA_ALLOC(arena,
                              (long)n_loop * (long)sizeof(int32_t));
    int32_t *collapse_dst = (int32_t *)ARENA_ALLOC(arena,
                              (long)n_loop * (long)sizeof(int32_t));
    size_t n_collapse = 0;
    {
        /* Relative threshold: collapse edges shorter than 0.1% of avg edge.
         * Also floor at 1e-12 to avoid false positives on tiny polygons. */
        double collapse_thr = avg_2d_edge * 1e-3;
        if (collapse_thr < 1e-6) collapse_thr = 1e-6;
        double min_dist_sq = collapse_thr * collapse_thr;
        size_t dst = 0;
        for (size_t i = 0; i < n_loop; i++) {
            size_t j = (i + 1) % n_loop;
            double dx = pts_2d[j * 2 + 0] - pts_2d[i * 2 + 0];
            double dy = pts_2d[j * 2 + 1] - pts_2d[i * 2 + 1];
            double d2 = dx * dx + dy * dy;
            if (d2 >= min_dist_sq || j == 0) {
                /* Keep vertex i (good edge, or last→first wrap) */
                pts_2d[dst * 2 + 0] = pts_2d[i * 2 + 0];
                pts_2d[dst * 2 + 1] = pts_2d[i * 2 + 1];
                bdry_3d[dst * 3 + 0] = bdry_3d[i * 3 + 0];
                bdry_3d[dst * 3 + 1] = bdry_3d[i * 3 + 1];
                bdry_3d[dst * 3 + 2] = bdry_3d[i * 3 + 2];
                bdry_map[dst] = bdry_map[i];
                dst++;
            } else {
                /* Zero-length edge i→j: drop vertex i, merge into j.
                 * Record the merge (mesh vertex i → mesh vertex j) only
                 * when they are distinct mesh vertices. */
                HFLOG("    [fill_one_hole] collapsing zero-length edge: "
                        "v%d -> v%d (2D dist=%.2e)\n",
                        bdry_map[i], bdry_map[j], sqrt(d2));
                if (bdry_map[i] != bdry_map[j]) {
                    collapse_src[n_collapse] = bdry_map[i];
                    collapse_dst[n_collapse] = bdry_map[j];
                    n_collapse++;
                }
            }
        }
        if (dst < n_loop) {
            HFLOG("    [fill_one_hole] collapsed %zu -> %zu loop verts\n",
                    n_loop, dst);
            n_loop = dst;
        }
        if (n_loop < 3) {
            HFLOG("    [fill_one_hole] loop too small after collapse (%zu)\n",
                    n_loop);
            Arena_restore(arena, mark);
            return -1;
        }
    }
    /* Also check wrap-around edge (last→first) with same threshold */
    {
        double collapse_thr_wrap = avg_2d_edge * 1e-3;
        if (collapse_thr_wrap < 1e-6) collapse_thr_wrap = 1e-6;
        double min_dist_sq_wrap = collapse_thr_wrap * collapse_thr_wrap;
        double dx = pts_2d[0] - pts_2d[(n_loop - 1) * 2 + 0];
        double dy = pts_2d[1] - pts_2d[(n_loop - 1) * 2 + 1];
        if (dx * dx + dy * dy < min_dist_sq_wrap) {
            HFLOG("    [fill_one_hole] collapsing wrap-around zero-length edge: "
                    "v%d -> v%d\n", bdry_map[n_loop - 1], bdry_map[0]);
            if (bdry_map[n_loop - 1] != bdry_map[0]) {
                collapse_src[n_collapse] = bdry_map[n_loop - 1];
                collapse_dst[n_collapse] = bdry_map[0];
                n_collapse++;
            }
            n_loop--;
            if (n_loop < 3) {
                HFLOG("    [fill_one_hole] loop too small after wrap collapse (%zu)\n",
                        n_loop);
                Arena_restore(arena, mark);
                return -1;
            }
        }
    }

    /* Check if simple polygon and valid for CDT */
    HFLOG("    [fill_one_hole] checking is_simple_polygon (n=%zu)...\n", n_loop);
    HFFLUSH();
    int is_simple = is_simple_polygon(pts_2d, n_loop);
    HFLOG("    [fill_one_hole] is_simple=%d\n", is_simple);
    HFFLUSH();
    int is_valid = is_valid_for_cdt(pts_2d, n_loop);
    HFLOG("    [fill_one_hole] is_valid_for_cdt=%d\n", is_valid);
    HFFLUSH();
    int good_for_cdt = is_simple && is_valid;

    /* If PCA projection fails, try axis-aligned projections (XY, XZ, YZ) */
    if (!good_for_cdt) {
        HFLOG("    [fill_one_hole] PCA projection not good for CDT, trying axis-aligned...\n");
        HFFLUSH();
        for (size_t ap = 0; ap < 3; ap++) {
            size_t d0 = ap < 2 ? 0 : 1;  /* XY->0, XZ->0, YZ->1 */
            size_t d1 = ap == 0 ? 1 : 2;  /* XY->1, XZ->2, YZ->2 */

            double *test_2d = (double *)ARENA_ALLOC(arena,
                                (long)n_loop * 2 * (long)sizeof(double));
            double mu = 0.0, mv = 0.0;
            for (size_t i = 0; i < n_loop; i++) {
                test_2d[i * 2 + 0] = (double)bdry_3d[i * 3 + d0];
                test_2d[i * 2 + 1] = (double)bdry_3d[i * 3 + d1];
                mu += test_2d[i * 2 + 0];
                mv += test_2d[i * 2 + 1];
            }
            mu /= (double)n_loop;
            mv /= (double)n_loop;
            for (size_t i = 0; i < n_loop; i++) {
                test_2d[i * 2 + 0] -= mu;
                test_2d[i * 2 + 1] -= mv;
            }

            HFLOG("    [fill_one_hole] axis %zu: checking is_simple (n=%zu)...\n",
                    ap, n_loop);
            HFFLUSH();
            int ax_simple = is_simple_polygon(test_2d, n_loop);
            HFLOG("    [fill_one_hole] axis %zu: is_simple=%d\n", ap, ax_simple);
            HFFLUSH();
            if (ax_simple &&
                is_valid_for_cdt(test_2d, n_loop)) {
                memcpy(pts_2d, test_2d, n_loop * 2 * sizeof(double));
                u_axis[0] = 0.0f; u_axis[1] = 0.0f; u_axis[2] = 0.0f;
                v_axis[0] = 0.0f; v_axis[1] = 0.0f; v_axis[2] = 0.0f;
                u_axis[d0] = 1.0f;
                v_axis[d1] = 1.0f;
                good_for_cdt = 1;
                break;
            }
        }
        /* If axis-aligned also failed, will use Clipper2 path below */
    }

    /* Enforce CCW winding before CDT */
    if (good_for_cdt) {
        double signed_area2 = 0.0;
        for (size_t i = 0; i < n_loop; i++) {
            size_t j = (i + 1) % n_loop;
            signed_area2 += pts_2d[i * 2] * pts_2d[j * 2 + 1]
                          - pts_2d[j * 2] * pts_2d[i * 2 + 1];
        }
        if (signed_area2 < 0.0) {
            for (size_t i = 0; i < n_loop / 2; i++) {
                size_t j = n_loop - 1 - i;
                double t0 = pts_2d[i * 2];
                double t1 = pts_2d[i * 2 + 1];
                pts_2d[i * 2]     = pts_2d[j * 2];
                pts_2d[i * 2 + 1] = pts_2d[j * 2 + 1];
                pts_2d[j * 2]     = t0;
                pts_2d[j * 2 + 1] = t1;

                float s0 = bdry_3d[i * 3];
                float s1 = bdry_3d[i * 3 + 1];
                float s2 = bdry_3d[i * 3 + 2];
                bdry_3d[i * 3]     = bdry_3d[j * 3];
                bdry_3d[i * 3 + 1] = bdry_3d[j * 3 + 1];
                bdry_3d[i * 3 + 2] = bdry_3d[j * 3 + 2];
                bdry_3d[j * 3]     = s0;
                bdry_3d[j * 3 + 1] = s1;
                bdry_3d[j * 3 + 2] = s2;

                int32_t m_tmp = bdry_map[i];
                bdry_map[i] = bdry_map[j];
                bdry_map[j] = m_tmp;
            }
        }
    }

    double *tri_pts_2d = NULL;
    int32_t *tri_faces_raw = NULL;
    size_t tri_nv = 0, tri_nf = 0;

    if (good_for_cdt) {
        /* Direct CDT — skip quality refinement if edges were collapsed.
         * Collapsed polygons have degenerate geometry that can crash Triangle's
         * quality refinement non-deterministically. */
        int skip_qual = (n_collapse > 0) ? 1 : 0;
        HFLOG("    [fill_one_hole] calling triangulate_polygon (direct CDT, n=%zu, skip_quality=%d)...\n",
                n_loop, skip_qual);
        HFFLUSH();
        rc = triangulate_polygon(pts_2d, n_loop,
                                 &tri_pts_2d, &tri_nv,
                                 &tri_faces_raw, &tri_nf,
                                 skip_qual);
        HFLOG("    [fill_one_hole] triangulate_polygon returned rc=%d nv=%zu nf=%zu\n",
                rc, tri_nv, tri_nf);
        HFFLUSH();
    } else {
        /* Self-intersecting boundary: use Clipper2 to resolve, then CDT.
         * After CDT, remap boundary vertices back to original loop via
         * nearest-neighbor in 2D so stitch can share mesh vertices. */
        double *clip_pts = NULL;
        int *clip_counts = NULL;
        int clip_n_polys = 0;

        HFLOG("    [fill_one_hole] calling Clipper2_union (n=%zu)...\n", n_loop);
        HFFLUSH();
        rc = Clipper2_union(pts_2d, n_loop, &clip_pts, &clip_counts,
                            &clip_n_polys);
        HFLOG("    [fill_one_hole] Clipper2_union returned rc=%d, n_polys=%d\n",
                rc, clip_n_polys);
        HFFLUSH();

        if (rc != 0 || clip_n_polys <= 0 || !clip_pts || !clip_counts) {
            free(clip_pts);
            free(clip_counts);
            Arena_restore(arena, mark);
            return -1;
        }

        /* Take largest Clipper2 sub-polygon */
        int best_poly = 0;
        size_t best_n = (size_t)clip_counts[0];
        size_t poly_offset = 0;
        for (int p = 0; p < clip_n_polys; p++) {
            if ((size_t)clip_counts[p] > best_n) {
                best_n = (size_t)clip_counts[p];
                best_poly = p;
            }
        }
        poly_offset = 0;
        for (int p = 0; p < best_poly; p++) {
            poly_offset += (size_t)clip_counts[p];
        }
        double *clip_poly = &clip_pts[poly_offset * 2];
        size_t clip_n = best_n;

        if (clip_n < 3) {
            free(clip_pts);
            free(clip_counts);
            Arena_restore(arena, mark);
            return -1;
        }

        /* Enforce CCW winding on Clipper2 sub-polygon */
        {
            double sa2 = 0.0;
            for (size_t i = 0; i < clip_n; i++) {
                size_t j = (i + 1) % clip_n;
                sa2 += clip_poly[i * 2] * clip_poly[j * 2 + 1]
                     - clip_poly[j * 2] * clip_poly[i * 2 + 1];
            }
            if (sa2 < 0.0) {
                for (size_t i = 0; i < clip_n / 2; i++) {
                    size_t j = clip_n - 1 - i;
                    double t0 = clip_poly[i * 2];
                    double t1 = clip_poly[i * 2 + 1];
                    clip_poly[i * 2]     = clip_poly[j * 2];
                    clip_poly[i * 2 + 1] = clip_poly[j * 2 + 1];
                    clip_poly[j * 2]     = t0;
                    clip_poly[j * 2 + 1] = t1;
                }
            }
        }

        /* CDT the cleaned polygon */
        double *clip_tri_pts = NULL;
        int32_t *clip_tri_faces = NULL;
        size_t clip_tri_nv = 0, clip_tri_nf = 0;

        HFLOG("    [fill_one_hole] CDT on Clipper2 output (clip_n=%zu)...\n",
                clip_n);
        HFFLUSH();
        /* skip_quality=1: the Clipper2 output is a post-boolean fallback polygon
         * (the input loop self-intersected). It commonly has sharp / near-
         * degenerate corners that send Triangle's q20 refinement into the
         * infinite Steiner loop that hung the 100-cube weld. We only need a VALID
         * fill here, not pretty triangles (the interior is CG-smoothed after), so
         * use a plain constrained Delaunay -- which always terminates. */
        rc = triangulate_polygon(clip_poly, clip_n,
                                 &clip_tri_pts, &clip_tri_nv,
                                 &clip_tri_faces, &clip_tri_nf,
                                 1);
        HFLOG("    [fill_one_hole] Clipper2 CDT returned rc=%d nv=%zu nf=%zu\n",
                rc, clip_tri_nv, clip_tri_nf);
        HFFLUSH();
        free(clip_pts);
        free(clip_counts);

        if (rc != 0 || clip_tri_nv == 0 || clip_tri_nf == 0) {
            free(clip_tri_pts);
            free(clip_tri_faces);
            Arena_restore(arena, mark);
            return -1;
        }

        /* Map each CDT vertex to either an original loop vertex (by 2D
         * nearest-neighbor) or mark as a new Steiner point.
         * vertex_map[i] = original loop index (0..n_loop-1) if matched,
         *                 or -1 if unmatched (new interior/intersection pt). */
        int32_t *vertex_map = (int32_t *)ARENA_ALLOC(arena,
                                (long)clip_tri_nv * (long)sizeof(int32_t));
        uint8_t *loop_used = (uint8_t *)ARENA_CALLOC(arena,
                                (long)n_loop, (long)sizeof(uint8_t));

        double match_tol_sq = 1e-6;  /* 2D distance² threshold */

        /* Compute a reasonable tolerance from boundary edge lengths */
        {
            double min_edge_sq = 1e30;
            for (size_t i = 0; i < n_loop; i++) {
                size_t j = (i + 1) % n_loop;
                double dx = pts_2d[j * 2] - pts_2d[i * 2];
                double dy = pts_2d[j * 2 + 1] - pts_2d[i * 2 + 1];
                double d2 = dx * dx + dy * dy;
                if (d2 < min_edge_sq && d2 > 0.0) min_edge_sq = d2;
            }
            /* Tolerance = (min_edge_len * 0.1)² */
            match_tol_sq = min_edge_sq * 0.01;
            if (match_tol_sq < 1e-12) match_tol_sq = 1e-12;
        }

        for (size_t i = 0; i < clip_tri_nv; i++) {
            vertex_map[i] = -1;
            double best_d2 = match_tol_sq;
            int32_t best_j = -1;
            for (size_t j = 0; j < n_loop; j++) {
                if (loop_used[j]) continue;
                double dx = clip_tri_pts[i * 2] - pts_2d[j * 2];
                double dy = clip_tri_pts[i * 2 + 1] - pts_2d[j * 2 + 1];
                double d2 = dx * dx + dy * dy;
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_j = (int32_t)j;
                }
            }
            if (best_j >= 0) {
                vertex_map[i] = best_j;
                loop_used[best_j] = 1;
            }
        }

        /* Count unmatched (extra/Steiner) verts */
        size_t n_extra = 0;
        size_t n_matched = 0;
        for (size_t i = 0; i < clip_tri_nv; i++) {
            if (vertex_map[i] < 0) n_extra++;
            else n_matched++;
        }
        /* Count unmatched original loop verts */
        size_t n_loop_unmatched = 0;
        for (size_t j = 0; j < n_loop; j++) {
            if (!loop_used[j]) n_loop_unmatched++;
        }
        if (n_loop_unmatched > 0) {
            fprintf(stderr, "  [fill_one_hole] Clipper2 path: %zu/%zu loop verts "
                    "unmatched\n", n_loop_unmatched, n_loop);
        }

        /* Build remapped output: first n_loop boundary verts, then extras.
         * extra_map[i] = new index for unmatched CDT vertex i */
        size_t remap_nv = n_loop + n_extra;
        int32_t *remap = (int32_t *)ARENA_ALLOC(arena,
                           (long)clip_tri_nv * (long)sizeof(int32_t));
        size_t extra_idx = n_loop;
        for (size_t i = 0; i < clip_tri_nv; i++) {
            if (vertex_map[i] >= 0) {
                remap[i] = vertex_map[i]; /* maps to original loop index */
            } else {
                remap[i] = (int32_t)extra_idx;
                extra_idx++;
            }
        }

        /* Build remapped face array */
        tri_nf = clip_tri_nf;
        tri_nv = remap_nv;
        tri_pts_2d = (double *)malloc(remap_nv * 2 * sizeof(double));
        tri_faces_raw = (int32_t *)malloc(tri_nf * 3 * sizeof(int32_t));
        if (!tri_pts_2d || !tri_faces_raw) {
            free(tri_pts_2d);
            free(tri_faces_raw);
            free(clip_tri_pts);
            free(clip_tri_faces);
            Arena_restore(arena, mark);
            return -1;
        }

        /* Fill 2D coords: boundary from original, extras from CDT */
        for (size_t i = 0; i < n_loop; i++) {
            tri_pts_2d[i * 2 + 0] = pts_2d[i * 2 + 0];
            tri_pts_2d[i * 2 + 1] = pts_2d[i * 2 + 1];
        }
        for (size_t i = 0; i < clip_tri_nv; i++) {
            if (vertex_map[i] < 0) {
                size_t dst = (size_t)remap[i];
                tri_pts_2d[dst * 2 + 0] = clip_tri_pts[i * 2 + 0];
                tri_pts_2d[dst * 2 + 1] = clip_tri_pts[i * 2 + 1];
            }
        }

        /* Remap face indices */
        for (size_t f = 0; f < tri_nf; f++) {
            for (size_t k = 0; k < 3; k++) {
                int32_t old_idx = clip_tri_faces[f * 3 + k];
                tri_faces_raw[f * 3 + k] = remap[old_idx];
            }
        }

        free(clip_tri_pts);
        free(clip_tri_faces);
    }

    if (rc != 0 || tri_nv == 0 || tri_nf == 0) {
        free(tri_pts_2d);
        free(tri_faces_raw);
        Arena_restore(arena, mark);
        return -1;
    }

    /* Embed in 3D: boundary verts from mesh, Steiner from PCA plane */
    HFLOG("    [fill_one_hole] embedding %zu verts in 3D (%zu boundary, %zu Steiner)...\n",
            tri_nv, n_loop, tri_nv - n_loop);
    HFFLUSH();
    float *fill_3d = (float *)ARENA_ALLOC(arena,
                       (long)tri_nv * 3 * (long)sizeof(float));
    embed_steiner_3d(tri_pts_2d, tri_nv, bdry_3d, n_loop,
                     u_axis, v_axis, centroid, fill_3d);

    /* Copy faces to arena */
    int32_t *fill_faces = (int32_t *)ARENA_ALLOC(arena,
                            (long)tri_nf * 3 * (long)sizeof(int32_t));
    memcpy(fill_faces, tri_faces_raw, tri_nf * 3 * sizeof(int32_t));

    free(tri_pts_2d);
    free(tri_faces_raw);

    /* Manifold pre-check: CDT may create diagonals between two
     * boundary verts (no Steiner involvement). If any such diagonal
     * already exists in the mesh as a 2-face interior edge, stitching
     * the fill creates a 3-fan non-manifold edge. Skip the fill rather
     * than emit non-manifold geometry. Steiner-involving diagonals
     * can't collide (Steiner is a fresh vert). */
    {
        int found_bad = 0;
        for (size_t f = 0; f < tri_nf && !found_bad; f++) {
            int32_t local_idx[3] = {
                fill_faces[f*3+0], fill_faces[f*3+1], fill_faces[f*3+2] };
            for (int e = 0; e < 3 && !found_bad; e++) {
                int32_t la = local_idx[e];
                int32_t lb = local_idx[(e+1)%3];
                if (la >= (int32_t)n_loop || lb >= (int32_t)n_loop) continue;
                /* Diagonal between two boundary verts — check if the
                 * mesh edge (bdry_map[la], bdry_map[lb]) is already
                 * incident to two faces. If so, abort the fill. */
                /* Only care about non-adjacent loop verts; consecutive
                 * loop verts are boundary edges (count=1 in mesh)
                 * which become manifold (count=2) when the fill
                 * stitches over them. */
                size_t diff = (size_t)((la > lb) ? (la - lb) : (lb - la));
                if (diff == 1 || diff == n_loop - 1) continue;
                int32_t mesh_a = bdry_map[la];
                int32_t mesh_b = bdry_map[lb];
                if (loop_diag_in_mesh(mesh_a, mesh_b, mesh_faces, mesh_nf)) {
                    found_bad = 1;
                }
            }
        }
        if (found_bad) {
            HFLOG("    [fill_one_hole] manifold pre-check: boundary-to-"
                  "boundary diagonal collides with mesh edge — skip fill\n");
            Arena_restore(arena, mark);
            return -1;
        }
    }

    /* Smooth Steiner points with cotangent Laplacian CG solve */
    if (tri_nv > n_loop) {
        HFLOG("    [fill_one_hole] smoothing %zu Steiner points...\n",
                tri_nv - n_loop);
        HFFLUSH();
        smooth_fill_cotangent(arena, fill_3d, fill_faces, tri_nf,
                              tri_nv, n_loop);
        HFLOG("    [fill_one_hole] smoothing done\n");
        HFFLUSH();

        /* Smoothing can collapse interior Steiner points of a thin-slit fill
         * into a sub-voxel cluster (needle triangles). Merge the coincident
         * Steiner verts and drop the slivers -- safe because Steiner verts are
         * fill-private until stitch_fills runs. Set HOLEFILL_NO_PRUNE=1 to
         * disable (A/B diagnostics). */
        static int prune_off = -1;
        if (prune_off < 0) { const char *e = sf_env("HOLEFILL_NO_PRUNE");
                             prune_off = (e && e[0] && e[0] != '0') ? 1 : 0; }
        size_t pre_nv = tri_nv, pre_nf = tri_nf;
        if (!prune_off &&
            prune_degenerate_fill(fill_3d, fill_faces, &tri_nv, &tri_nf, n_loop,
                                  FILL_DEGEN_EPS_VOX, FILL_DEGEN_AREA)) {
            HFLOG("    [fill_one_hole] pruned degenerate fill: "
                    "%zu->%zu verts, %zu->%zu faces\n",
                    pre_nv, tri_nv, pre_nf, tri_nf);
        }
    }

    /* Orient fill boundary to oppose mesh boundary half-edges */
    HFLOG("    [fill_one_hole] orienting fill faces...\n");
    HFFLUSH();
    orient_fill(fill_faces, tri_nf,
                mesh_faces, mesh_nf,
                bdry_map, n_loop);

    HFLOG("    [fill_one_hole] done: %zu verts, %zu faces\n", tri_nv, tri_nf);
    HFFLUSH();

    result->verts = fill_3d;
    result->faces = fill_faces;
    result->nv = tri_nv;
    result->nf = tri_nf;
    result->n_boundary = n_loop;
    result->boundary_map = bdry_map;
    result->merge_src = collapse_src;
    result->merge_dst = collapse_dst;
    result->n_merge = n_collapse;
    return 0;
}

/* ------------------------------------------------------------------ */
/* fill_micro_hole — fan triangulation for tiny holes (<=6 verts)      */
/*                                                                     */
/* Manifold pre-check: the naive fan-from-vertex-0 adds diagonal edges */
/* (loop[0], loop[2]), (loop[0], loop[3]), etc. If one of those edges  */
/* already exists in the mesh as a 2-face interior edge, the fan       */
/* triangles would make it the third reference, creating a non-        */
/* manifold edge. To detect this we look at the loop vertices'         */
/* incident faces and check: does any pair (loop[0], loop[k]) for      */
/* k = 2..n_loop-2 already appear as a mesh edge? If yes, try fanning  */
/* from a different apex. If no apex is safe, return error and let the */
/* hole stay open (BPA output's "phantom" boundaries occur on already- */
/* closed regions; leaving them open is harmless, filling them creates */
/* non-manifold output).                                               */
/* ------------------------------------------------------------------ */
static int loop_diag_in_mesh(int32_t va, int32_t vb,
                             const int32_t *mesh_faces, size_t mesh_nf)
{
    if (va == vb) return 0;
    for (size_t f = 0; f < mesh_nf; f++) {
        int32_t a = mesh_faces[f*3+0];
        int32_t b = mesh_faces[f*3+1];
        int32_t c = mesh_faces[f*3+2];
        if ((a == va && b == vb) || (b == va && a == vb)) return 1;
        if ((b == va && c == vb) || (c == va && b == vb)) return 1;
        if ((a == va && c == vb) || (c == va && a == vb)) return 1;
    }
    return 0;
}

static int fan_safe(int apex_idx, const int32_t *loop_verts, size_t n_loop,
                    const int32_t *mesh_faces, size_t mesh_nf)
{
    int32_t apex_v = loop_verts[apex_idx];
    for (size_t k = 0; k < n_loop; k++) {
        if ((int)k == apex_idx) continue;
        size_t k_next = (k + 1) % n_loop;
        size_t k_prev = (k + n_loop - 1) % n_loop;
        if ((int)k_next == apex_idx) continue;
        if ((int)k_prev == apex_idx) continue;
        /* k is not the apex and not a boundary-edge neighbour of the
         * apex; the edge (apex, k) is a diagonal that the fan would
         * introduce. Reject if already in the mesh. */
        if (loop_diag_in_mesh(apex_v, loop_verts[k], mesh_faces, mesh_nf))
            return 0;
    }
    return 1;
}

static int fill_micro_hole(Arena_T arena,
                           const float *mesh_verts,
                           const int32_t *mesh_faces, size_t mesh_nf,
                           const int32_t *loop_verts, size_t n_loop,
                           HoleFillResult *result)
{
    if (n_loop < 3) return -1;

    /* Manifold-preserving apex selection: try each loop vertex as the
     * fan apex; the first one whose diagonals don't already exist in
     * the mesh wins. */
    int apex = -1;
    for (size_t i = 0; i < n_loop; i++) {
        if (fan_safe((int)i, loop_verts, n_loop, mesh_faces, mesh_nf)) {
            apex = (int)i;
            break;
        }
    }
    if (apex < 0) return -2; /* every fan would create non-manifold edges */

    /* Fan triangulation from `apex` vertex. The output uses local
     * indices 0..n_loop-1 in the fill_faces; the apex gets local
     * index 0 by reordering loop into [apex, apex+1, ..., apex-1]. */
    size_t n_tri = n_loop - 2;
    float *fill_verts = (float *)ARENA_ALLOC(arena,
                          (long)n_loop * 3 * (long)sizeof(float));
    int32_t *fill_faces = (int32_t *)ARENA_ALLOC(arena,
                            (long)n_tri * 3 * (long)sizeof(int32_t));

    /* Rotate loop so apex is at index 0. */
    int32_t *rotated = (int32_t *)ARENA_ALLOC(arena,
                          (long)n_loop * (long)sizeof(int32_t));
    for (size_t i = 0; i < n_loop; i++) {
        rotated[i] = loop_verts[(apex + i) % n_loop];
    }

    for (size_t i = 0; i < n_loop; i++) {
        fill_verts[i * 3 + 0] = mesh_verts[rotated[i] * 3 + 0];
        fill_verts[i * 3 + 1] = mesh_verts[rotated[i] * 3 + 1];
        fill_verts[i * 3 + 2] = mesh_verts[rotated[i] * 3 + 2];
    }

    for (size_t i = 0; i < n_tri; i++) {
        fill_faces[i * 3 + 0] = 0;
        fill_faces[i * 3 + 1] = (int32_t)(i + 1);
        fill_faces[i * 3 + 2] = (int32_t)(i + 2);
    }

    /* Micro holes have no reversal — map is identity to rotated loop */
    int32_t *bmap = (int32_t *)ARENA_ALLOC(arena,
                      (long)n_loop * (long)sizeof(int32_t));
    for (size_t i = 0; i < n_loop; i++) {
        bmap[i] = rotated[i];
    }

    /* Orient fill boundary to oppose mesh boundary half-edges */
    orient_fill(fill_faces, n_tri, mesh_faces, mesh_nf, bmap, n_loop);

    result->verts = fill_verts;
    result->faces = fill_faces;
    result->nv = n_loop;
    result->nf = n_tri;
    result->n_boundary = n_loop;
    result->boundary_map = bmap;
    result->merge_src = NULL;
    result->merge_dst = NULL;
    result->n_merge = 0;
    return 0;
}

/* ------------------------------------------------------------------ */
/* stitch_fills — merge fill patches into the mesh                     */
/* ------------------------------------------------------------------ */
static int stitch_fills(Arena_T arena,
                        float **verts, int32_t **faces,
                        size_t *nv, size_t *nf,
                        const HoleFillResult *fills, size_t n_fills)
{
    if (n_fills == 0) return 0;

    /* Count total new Steiner verts and new faces */
    size_t total_steiner = 0;
    size_t total_new_faces = 0;
    for (size_t fi = 0; fi < n_fills; fi++) {
        if (fills[fi].nv < fills[fi].n_boundary) return -1;
        total_steiner += fills[fi].nv - fills[fi].n_boundary;
        total_new_faces += fills[fi].nf;
    }

    size_t old_nv = *nv;
    size_t old_nf = *nf;
    size_t new_nv = old_nv + total_steiner;
    size_t new_nf = old_nf + total_new_faces;

    float *new_verts = (float *)ARENA_ALLOC(arena,
                         (long)new_nv * 3 * (long)sizeof(float));
    int32_t *new_faces = (int32_t *)ARENA_ALLOC(arena,
                           (long)new_nf * 3 * (long)sizeof(int32_t));

    /* Copy original mesh */
    memcpy(new_verts, *verts, old_nv * 3 * sizeof(float));
    memcpy(new_faces, *faces, old_nf * 3 * sizeof(int32_t));

    size_t steiner_offset = old_nv;
    size_t face_offset = old_nf;

    for (size_t fi = 0; fi < n_fills; fi++) {
        const HoleFillResult *fill = &fills[fi];

        /* Copy Steiner vertices */
        for (size_t i = fill->n_boundary; i < fill->nv; i++) {
            size_t dst = steiner_offset + (i - fill->n_boundary);
            new_verts[dst * 3 + 0] = fill->verts[i * 3 + 0];
            new_verts[dst * 3 + 1] = fill->verts[i * 3 + 1];
            new_verts[dst * 3 + 2] = fill->verts[i * 3 + 2];
        }

        /* Copy and remap faces */
        for (size_t f = 0; f < fill->nf; f++) {
            for (size_t k = 0; k < 3; k++) {
                int32_t old_idx = fill->faces[f * 3 + k];
                int32_t new_idx;
                if ((size_t)old_idx < fill->n_boundary) {
                    /* Map to original mesh vertex via boundary_map
                     * (accounts for CCW reversal in fill_one_hole) */
                    new_idx = fill->boundary_map[old_idx];
                } else {
                    /* Map to new Steiner vertex */
                    new_idx = (int32_t)(steiner_offset +
                              (size_t)old_idx - fill->n_boundary);
                }
                new_faces[(face_offset + f) * 3 + k] = new_idx;
            }
        }

        steiner_offset += fill->nv - fill->n_boundary;
        face_offset += fill->nf;
    }

    /* Apply vertex merges from 2D edge collapses.
     * When fill_one_hole collapses a zero-length 2D edge (vertex A and B
     * differ in 3D but project to the same 2D point), vertex A is removed
     * from the fill boundary.  After stitching, A remains a mesh boundary
     * vertex with no fill connection.  Fix by merging A → B in the face
     * array so A is no longer referenced. */
    size_t total_merges_applied = 0;
    for (size_t fi = 0; fi < n_fills; fi++) {
        const HoleFillResult *fill = &fills[fi];
        if (fill->n_merge == 0) continue;
        for (size_t mi = 0; mi < fill->n_merge; mi++) {
            int32_t src = fill->merge_src[mi];
            int32_t dst_v = fill->merge_dst[mi];
            /* Replace all occurrences of src with dst_v in face array */
            for (size_t f = 0; f < new_nf * 3; f++) {
                if (new_faces[f] == src)
                    new_faces[f] = dst_v;
            }
            total_merges_applied++;
        }
    }
    if (total_merges_applied > 0) {
        /* Remove degenerate faces created by merges */
        size_t dst_f = 0;
        for (size_t f = 0; f < new_nf; f++) {
            int32_t a = new_faces[f * 3 + 0];
            int32_t b = new_faces[f * 3 + 1];
            int32_t c = new_faces[f * 3 + 2];
            if (a == b || b == c || a == c) continue;
            new_faces[dst_f * 3 + 0] = a;
            new_faces[dst_f * 3 + 1] = b;
            new_faces[dst_f * 3 + 2] = c;
            dst_f++;
        }
        fprintf(stderr, "  [stitch] applied %zu vertex merges from edge collapse, "
                "removed %zu degenerate faces\n",
                total_merges_applied, new_nf - dst_f);
        new_nf = dst_f;
    }

    *verts = new_verts;
    *faces = new_faces;
    *nv = new_nv;
    *nf = new_nf;
    return 0;
}

/* ------------------------------------------------------------------ */
/* repair_mesh_orientation — BFS propagation to fix flipped faces       */
/* Returns the number of faces flipped.                                */
/* ------------------------------------------------------------------ */

/* Comparator for half-edge entries: sort by (lo, hi) undirected edge */
static int cmp_he_entry(const void *a, const void *b)
{
    const int32_t *ea = (const int32_t *)a;
    const int32_t *eb = (const int32_t *)b;
    if (ea[0] != eb[0]) return (ea[0] < eb[0]) ? -1 : 1;
    if (ea[1] != eb[1]) return (ea[1] < eb[1]) ? -1 : 1;
    return 0; /* stable for same undirected edge */
}

static size_t repair_mesh_orientation(Arena_T arena,
                                       int32_t *faces, size_t nf)
{
    if (nf < 2) return 0;

    Arena_Mark mark = Arena_save(arena);

    /* Build half-edge table: for each half-edge (a→b), store:
     *   [0] lo = min(a,b)
     *   [1] hi = max(a,b)
     *   [2] face_idx
     *   [3] is_forward (1 if a < b, i.e. half-edge goes lo→hi) */
    size_t n_he = nf * 3;
    int32_t *he = (int32_t *)ARENA_ALLOC(arena,
                     (long)n_he * 4 * (long)sizeof(int32_t));

    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f * 3 + e];
            int32_t b = faces[f * 3 + ((e + 1) % 3)];
            size_t idx = f * 3 + (size_t)e;
            he[idx * 4 + 0] = (a < b) ? a : b;
            he[idx * 4 + 1] = (a < b) ? b : a;
            he[idx * 4 + 2] = (int32_t)f;
            he[idx * 4 + 3] = (a < b) ? 1 : 0;
        }
    }

    /* Sort by (lo, hi) to group shared edges */
    qsort(he, n_he, 4 * sizeof(int32_t), cmp_he_entry);

    /* Build face adjacency: for each manifold edge (shared by exactly 2
     * faces), record (face_a, face_b, same_direction).
     * same_direction = 1 means both half-edges go the same way → need flip. */
    size_t max_adj = nf * 3;
    int32_t *adj_f = (int32_t *)ARENA_ALLOC(arena,
                        (long)max_adj * 2 * (long)sizeof(int32_t));
    int8_t *adj_same = (int8_t *)ARENA_ALLOC(arena,
                          (long)max_adj * (long)sizeof(int8_t));
    size_t n_adj = 0;

    size_t i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && he[j * 4] == he[i * 4] && he[j * 4 + 1] == he[i * 4 + 1])
            j++;
        if (j - i == 2) {
            /* Manifold edge — two faces share it */
            int32_t fa = he[i * 4 + 2];
            int32_t fb = he[(i + 1) * 4 + 2];
            int fwd_a = he[i * 4 + 3];
            int fwd_b = he[(i + 1) * 4 + 3];
            /* If both go same direction → inconsistent → need flip */
            adj_f[n_adj * 2 + 0] = fa;
            adj_f[n_adj * 2 + 1] = fb;
            adj_same[n_adj] = (int8_t)(fwd_a == fwd_b ? 1 : 0);
            n_adj++;
        }
        i = j;
    }

    /* Build face→adj_list using CSR-like structure */
    int32_t *f_degree = (int32_t *)ARENA_CALLOC(arena, (long)nf,
                                                   (long)sizeof(int32_t));
    for (size_t ei = 0; ei < n_adj; ei++) {
        f_degree[adj_f[ei * 2 + 0]]++;
        f_degree[adj_f[ei * 2 + 1]]++;
    }
    int32_t *f_off = (int32_t *)ARENA_ALLOC(arena,
                        (long)(nf + 1) * (long)sizeof(int32_t));
    f_off[0] = 0;
    for (size_t f = 0; f < nf; f++)
        f_off[f + 1] = f_off[f] + f_degree[f];

    /* adj_data[i] = (neighbor_face, same_direction_flag) */
    int32_t *adj_nb = (int32_t *)ARENA_ALLOC(arena,
                         (long)(n_adj * 2) * (long)sizeof(int32_t));
    int8_t *adj_sd = (int8_t *)ARENA_ALLOC(arena,
                        (long)(n_adj * 2) * (long)sizeof(int8_t));
    int32_t *cursor = (int32_t *)ARENA_ALLOC(arena,
                         (long)nf * (long)sizeof(int32_t));
    memcpy(cursor, f_off, nf * sizeof(int32_t));

    for (size_t ei = 0; ei < n_adj; ei++) {
        int32_t fa = adj_f[ei * 2 + 0];
        int32_t fb = adj_f[ei * 2 + 1];
        int8_t sd = adj_same[ei];
        adj_nb[cursor[fa]] = fb;  adj_sd[cursor[fa]] = sd;  cursor[fa]++;
        adj_nb[cursor[fb]] = fa;  adj_sd[cursor[fb]] = sd;  cursor[fb]++;
    }

    /* BFS from face 0.  Track whether each face should be flipped
     * relative to face 0 (which we define as "correct"). */
    uint8_t *visited = (uint8_t *)ARENA_CALLOC(arena, (long)nf,
                                                  (long)sizeof(uint8_t));
    uint8_t *should_flip = (uint8_t *)ARENA_CALLOC(arena, (long)nf,
                                                      (long)sizeof(uint8_t));
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena,
                        (long)nf * (long)sizeof(int32_t));
    size_t q_head = 0, q_tail = 0;

    visited[0] = 1;
    should_flip[0] = 0;
    queue[q_tail++] = 0;

    while (q_head < q_tail) {
        int32_t cf = queue[q_head++];
        for (int32_t ai = f_off[cf]; ai < f_off[cf + 1]; ai++) {
            int32_t nb = adj_nb[ai];
            if (visited[nb]) continue;
            visited[nb] = 1;
            /* If the shared edge has same_direction (inconsistent),
             * the neighbor needs to be flipped relative to current face.
             * XOR with current face's flip state. */
            int8_t sd = adj_sd[ai];
            should_flip[nb] = (uint8_t)(should_flip[cf] ^ (uint8_t)sd);
            queue[q_tail++] = nb;
        }
    }

    /* Apply flips */
    size_t n_flipped = 0;
    for (size_t f = 0; f < nf; f++) {
        if (should_flip[f]) {
            /* Swap v1 and v2 to reverse winding */
            int32_t tmp = faces[f * 3 + 1];
            faces[f * 3 + 1] = faces[f * 3 + 2];
            faces[f * 3 + 2] = tmp;
            n_flipped++;
        }
    }

    if (n_flipped > 0) {
        fprintf(stderr, "  [hole_fill] orientation repair: flipped %zu/%zu faces\n",
                n_flipped, nf);
    }

    Arena_restore(arena, mark);
    return n_flipped;
}

/* ------------------------------------------------------------------ */
/* Directed boundary half-edges + geometric interior/exterior test     */
/* ------------------------------------------------------------------ */
typedef struct { int32_t lo, hi; int32_t a, b; int32_t face; } DirBHE;

static int cmp_dirbhe(const void *pa, const void *pb)
{
    const DirBHE *x = (const DirBHE *)pa, *y = (const DirBHE *)pb;
    if (x->lo != y->lo) return x->lo < y->lo ? -1 : 1;
    if (x->hi != y->hi) return x->hi < y->hi ? -1 : 1;
    return 0;
}

/* Directed boundary half-edges (edges in exactly one face), each carrying its
 * in-face direction (a->b, surface on its left) and face index, sorted by
 * undirected key for bsearch. */
static size_t build_dir_boundary(Arena_T arena, const int32_t *faces, size_t nf,
                                 DirBHE **out)
{
    size_t n_he = nf * 3;
    DirBHE *all = (DirBHE *)ARENA_ALLOC(arena,
                    (long)(n_he ? n_he : 1) * (long)sizeof(DirBHE));
    for (size_t f = 0; f < nf; f++) {
        int32_t tri[3] = { faces[f*3+0], faces[f*3+1], faces[f*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = tri[e], b = tri[(e+1)%3];
            DirBHE *h = &all[f*3 + (size_t)e];
            h->a = a; h->b = b; h->face = (int32_t)f;
            h->lo = a < b ? a : b; h->hi = a < b ? b : a;
        }
    }
    qsort(all, n_he, sizeof(DirBHE), cmp_dirbhe);
    DirBHE *bhe = (DirBHE *)ARENA_ALLOC(arena,
                    (long)(n_he ? n_he : 1) * (long)sizeof(DirBHE));
    size_t nb = 0, i = 0;
    while (i < n_he) {
        size_t j = i + 1;
        while (j < n_he && all[j].lo == all[i].lo && all[j].hi == all[i].hi) j++;
        if (j - i == 1) bhe[nb++] = all[i];
        i = j;
    }
    *out = bhe;
    return nb;
}

static const DirBHE *find_dirbhe(const DirBHE *bhe, size_t n, int32_t u, int32_t w)
{
    DirBHE key; key.lo = u < w ? u : w; key.hi = u < w ? w : u;
    return (const DirBHE *)bsearch(&key, bhe, n, sizeof(DirBHE), cmp_dirbhe);
}

/* Geometric interior/exterior classification of a boundary loop. INTERIOR (the
 * surface surrounds the loop) iff its signed area -- taken in the plane of the
 * averaged incident-face normal, with the traversal aligned to the boundary
 * half-edges (surface on the left) -- is NEGATIVE, i.e. it winds opposite the
 * outer perimeter. World-coord agnostic; relies only on consistent winding
 * (repair_mesh_orientation ran above). Returns 1 interior, 0 exterior or
 * indeterminate (fail-safe: indeterminate => do not fill). */
static int loop_is_interior(const float *verts, const int32_t *faces,
                            const DirBHE *bhe, size_t n_bhe,
                            const int32_t *loop_verts, size_t loop_n)
{
    if (loop_n < 3) return 0;
    double C[3] = {0,0,0};
    for (size_t k = 0; k < loop_n; k++) {
        const float *p = &verts[(size_t)loop_verts[k]*3];
        C[0]+=p[0]; C[1]+=p[1]; C[2]+=p[2];
    }
    C[0]/=(double)loop_n; C[1]/=(double)loop_n; C[2]/=(double)loop_n;

    double N[3] = {0,0,0};
    int dir = 0;
    for (size_t k = 0; k < loop_n; k++) {
        int32_t a = loop_verts[k], b = loop_verts[(k+1)%loop_n];
        const DirBHE *h = find_dirbhe(bhe, n_bhe, a, b);
        if (!h) continue;
        const float *p0 = &verts[(size_t)faces[(size_t)h->face*3+0]*3];
        const float *p1 = &verts[(size_t)faces[(size_t)h->face*3+1]*3];
        const float *p2 = &verts[(size_t)faces[(size_t)h->face*3+2]*3];
        double e1[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
        double e2[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
        N[0] += e1[1]*e2[2]-e1[2]*e2[1];
        N[1] += e1[2]*e2[0]-e1[0]*e2[2];
        N[2] += e1[0]*e2[1]-e1[1]*e2[0];
        if (dir == 0) dir = (h->a == a && h->b == b) ? 1 : -1;
    }
    if (dir == 0) return 0;
    if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0) return 0;

    double S = 0.0;
    for (size_t k = 0; k < loop_n; k++) {
        const float *pa = &verts[(size_t)loop_verts[k]*3];
        const float *pb = &verts[(size_t)loop_verts[(k+1)%loop_n]*3];
        double qa[3] = { pa[0]-C[0], pa[1]-C[1], pa[2]-C[2] };
        double qb[3] = { pb[0]-C[0], pb[1]-C[1], pb[2]-C[2] };
        S += (qa[1]*qb[2]-qa[2]*qb[1])*N[0]
           + (qa[2]*qb[0]-qa[0]*qb[2])*N[1]
           + (qa[0]*qb[1]-qa[1]*qb[0])*N[2];
    }
    if ((dir * S) < 0.0) return 1;

    /* Degenerate-projection tie-break. A tiny TWISTED loop (a bowtie quad --
     * observed as 4-edge slots the seam bridge leaves at a grazing seam) has
     * its two projected triangles cancel to S ~ 0, landing in the fail-safe
     * "indeterminate => exterior" arm and never filling. A genuine exterior
     * loop is a perimeter with hundreds of edges, so for small loops the
     * fail-safe points the wrong way: fill. eps is relative to the loop's own
     * scale (S ~ |N| * diam^2). */
    if (loop_n <= 8) {
        double nmag = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
        double diam2 = 1.0;
        for (size_t k = 0; k < loop_n; k++) {
            const float *pa = &verts[(size_t)loop_verts[k]*3];
            double dz = pa[0]-C[0], dy = pa[1]-C[1], dx = pa[2]-C[2];
            double d2 = dz*dz + dy*dy + dx*dx;
            if (d2 > diam2) diam2 = d2;
        }
        if (fabs(S) < 1e-6 * nmag * diam2) return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* dump_hole_obj — debug: write one hole's geometry to its own OBJ.    */
/*                                                                     */
/* PRE  (res == NULL): the input boundary loop as a tri-fan, RED.      */
/*   Written BEFORE fill_one_hole, so if Triangle hangs on this hole   */
/*   its polygon is the last *_in.obj on disk -- the exact culprit.    */
/*   World coords, so it overlays the --dump-stages OBJs in MeshLab.   */
/* POST (res != NULL): the filled patch (res->verts/faces), GREEN.     */
/*   Written after a successful fill.                                  */
/* Enabled only when dir != NULL (SEAM_HOLE_DUMP_DIR); no-op otherwise.*/
/* ------------------------------------------------------------------ */
static void dump_hole_obj(const char *dir, const char *prefix, size_t seq,
                          size_t loop_idx, size_t loop_n,
                          const float *mesh_verts, const int32_t *loop_verts,
                          const HoleFillResult *res)
{
    if (!dir) return;
    char path[1024];
    if (res == NULL) {
        snprintf(path, sizeof(path), "%s/%s_hole%04zu_loop%zu_%zuv_in.obj",
                 dir, prefix, seq, loop_idx, loop_n);
    } else {
        snprintf(path, sizeof(path), "%s/%s_hole%04zu_loop%zu_filled.obj",
                 dir, prefix, seq, loop_idx);
    }
    ves_ensure_parent_dir(path);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "  [hole_fill] WARNING: cannot write hole dump %s\n", path);
        return;
    }
    if (res == NULL) {
        fprintf(fp, "# hole %zu boundary loop %zu (%zu verts) -- PRE-fill input\n",
                seq, loop_idx, loop_n);
        for (size_t i = 0; i < loop_n; i++) {
            int32_t v = loop_verts[i];
            fprintf(fp, "v %.6f %.6f %.6f 1.0 0.0 0.0\n",
                    (double)mesh_verts[v * 3 + 0],
                    (double)mesh_verts[v * 3 + 1],
                    (double)mesh_verts[v * 3 + 2]);
        }
        /* tri-fan from vertex 0 so it renders as a (possibly non-simple) patch */
        for (size_t i = 1; i + 1 < loop_n; i++) {
            fprintf(fp, "f 1 %zu %zu\n", i + 1, i + 2);
        }
    } else {
        fprintf(fp, "# hole %zu loop %zu -- filled patch (%zu verts, %zu faces)\n",
                seq, loop_idx, res->nv, res->nf);
        for (size_t i = 0; i < res->nv; i++) {
            fprintf(fp, "v %.6f %.6f %.6f 0.0 1.0 0.0\n",
                    (double)res->verts[i * 3 + 0],
                    (double)res->verts[i * 3 + 1],
                    (double)res->verts[i * 3 + 2]);
        }
        for (size_t f = 0; f < res->nf; f++) {
            fprintf(fp, "f %d %d %d\n",
                    res->faces[f * 3 + 0] + 1,
                    res->faces[f * 3 + 1] + 1,
                    res->faces[f * 3 + 2] + 1);
        }
    }
    fclose(fp);
}

/* ------------------------------------------------------------------ */
/* log_loop_geom — degeneracy stats of a boundary loop (mesh coords).  */
/* Logged per hole so a hole that fails / times out can be checked for */
/* the near-degenerate edges (tiny edge, extreme ratio, sub-degree     */
/* corner, near-coincident non-adjacent verts) that send Triangle's    */
/* CDT non-terminating. Helps decide how aggressively to pre-clean.    */
/* ------------------------------------------------------------------ */
static void log_loop_geom(const char *tag, size_t loop_idx,
                          const float *mesh_verts, const int32_t *loop_verts,
                          size_t loop_n)
{
    if (loop_n < 3) return;
    double min_e = 1e30, max_e = 0.0, min_ang = 360.0, min_pair = 1e30;
    for (size_t i = 0; i < loop_n; i++) {
        const float *p = &mesh_verts[(size_t)loop_verts[(i + loop_n - 1) % loop_n] * 3];
        const float *c = &mesh_verts[(size_t)loop_verts[i] * 3];
        const float *q = &mesh_verts[(size_t)loop_verts[(i + 1) % loop_n] * 3];
        double ex = c[0]-q[0], ey = c[1]-q[1], ez = c[2]-q[2];
        double e = sqrt(ex*ex + ey*ey + ez*ez);
        if (e < min_e) min_e = e;
        if (e > max_e) max_e = e;
        /* interior angle at c, between edges c->p and c->q */
        double ux = p[0]-c[0], uy = p[1]-c[1], uz = p[2]-c[2];
        double vx = q[0]-c[0], vy = q[1]-c[1], vz = q[2]-c[2];
        double nu = sqrt(ux*ux+uy*uy+uz*uz), nv = sqrt(vx*vx+vy*vy+vz*vz);
        if (nu > 1e-12 && nv > 1e-12) {
            double cosang = (ux*vx + uy*vy + uz*vz) / (nu*nv);
            if (cosang > 1.0) cosang = 1.0; else if (cosang < -1.0) cosang = -1.0;
            double ang = acos(cosang) * 57.29577951308232;
            if (ang < min_ang) min_ang = ang;
        }
    }
    /* nearest non-adjacent vertex pair (pinch / near-duplicate detector) */
    for (size_t i = 0; i < loop_n; i++) {
        const float *a = &mesh_verts[(size_t)loop_verts[i] * 3];
        for (size_t j = i + 2; j < loop_n; j++) {
            if (i == 0 && j == loop_n - 1) continue;   /* wrap-adjacent */
            const float *b = &mesh_verts[(size_t)loop_verts[j] * 3];
            double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d < min_pair) min_pair = d;
        }
    }
    fprintf(stderr, "  [hole_fill] %s loop %zu geom: min_edge=%.4g max_edge=%.4g "
            "ratio=%.1f min_angle=%.2f deg min_nonadj=%.4g\n",
            tag, loop_idx, min_e, max_e,
            (min_e > 1e-12 ? max_e / min_e : 0.0), min_ang, min_pair);
    HFFLUSH();
}

/* ------------------------------------------------------------------ */
/* split_pinched_loop — retopologize a self-pinching loop into simple  */
/* sub-cycles.                                                         */
/*                                                                     */
/* A weld can fuse two cubes' boundaries into a figure-8 "<><>": the   */
/* boundary passes through ONE position 2+ times -- either the same    */
/* mesh vertex index repeated, or two coincident-but-distinct indices  */
/* (the two cubes each contributed a vertex at the seam point). That   */
/* polygon is non-simple; feeding it to CDT hangs Triangle. The closed */
/* walk is the concatenation of simple cycles joined at the pinch       */
/* point(s) -- we recover them with a stack: walking the loop, whenever */
/* we re-reach a position already on the stack we pop the cycle in      */
/* between and emit it (Hopcroft-Tarjan style cycle peeling). Each      */
/* emitted cycle is a simple sub-loop the CDT path can fill.            */
/*                                                                     */
/* loop_verts[0..n-1]  mesh vertex indices, cyclic (no closing repeat). */
/* Coincidence is tested by position within `eps`. If the loop is       */
/* already simple, returns 1 with out_sub[0] ALIASING loop_verts (no    */
/* copy). Sub-cycles with < 3 verts are dropped. Caps at max_subs.      */
/* ------------------------------------------------------------------ */
#define MAX_SUBLOOPS 32
#define PINCH_EPS    0.05f   /* verts within this (vox) are the same pinch pt */

/* Non-static so the unit test (holefill_interior_test, gate G4) can drive it
 * directly with hand-built figure-8 loops -- there is no public prototype. */
size_t split_pinched_loop(Arena_T arena,
                                 int32_t *loop_verts, size_t n,
                                 const float *verts, float eps,
                                 int32_t **out_sub, size_t *out_len,
                                 size_t max_subs)
{
    if (n < 3 || max_subs == 0) return 0;

    /* Canonical position id per slot: smallest j<=i at the same position. */
    int32_t *cid = (int32_t *)ARENA_ALLOC(arena, (long)n * (long)sizeof(int32_t));
    float eps2 = eps * eps;
    int any_pinch = 0;
    for (size_t i = 0; i < n; i++) {
        cid[i] = (int32_t)i;
        const float *pi = &verts[(size_t)loop_verts[i] * 3];
        for (size_t j = 0; j < i; j++) {
            const float *pj = &verts[(size_t)loop_verts[j] * 3];
            float dz = pi[0]-pj[0], dy = pi[1]-pj[1], dx = pi[2]-pj[2];
            if (dz*dz + dy*dy + dx*dx <= eps2) { cid[i] = cid[j]; any_pinch = 1; break; }
        }
    }
    if (!any_pinch) {              /* already simple -- no copy */
        out_sub[0] = loop_verts;
        out_len[0] = n;
        return 1;
    }

    /* Peel simple cycles with a stack of slot indices; stack_pos[cid] is the
     * slot's stack position or -1. cid in [0,n), so size stack_pos by n. */
    int32_t *stack     = (int32_t *)ARENA_ALLOC(arena, (long)n * (long)sizeof(int32_t));
    int32_t *stack_pos = (int32_t *)ARENA_ALLOC(arena, (long)n * (long)sizeof(int32_t));
    for (size_t i = 0; i < n; i++) stack_pos[i] = -1;
    size_t top = 0, n_sub = 0;

    for (size_t i = 0; i < n; i++) {
        int32_t c = cid[i];
        if (stack_pos[c] != -1) {
            /* Re-reached position c: slots [k..top) form a closed cycle. */
            size_t k = (size_t)stack_pos[c];
            size_t clen = top - k;
            if (clen >= 3 && n_sub < max_subs) {
                int32_t *sub = (int32_t *)ARENA_ALLOC(arena,
                                  (long)clen * (long)sizeof(int32_t));
                for (size_t t = 0; t < clen; t++) sub[t] = loop_verts[stack[k + t]];
                out_sub[n_sub] = sub; out_len[n_sub] = clen; n_sub++;
            }
            /* Pop the cycle's interior slots (keep slot k as the join point). */
            for (size_t t = k + 1; t < top; t++) stack_pos[cid[stack[t]]] = -1;
            top = k + 1;
            /* slot i is the same position as slot k -> absorbed, not pushed. */
        } else {
            stack[top] = (int32_t)i;
            stack_pos[c] = (int32_t)top;
            top++;
        }
    }
    /* The walk is closed, so the residual stack is the outermost cycle. */
    if (top >= 3 && n_sub < max_subs) {
        int32_t *sub = (int32_t *)ARENA_ALLOC(arena, (long)top * (long)sizeof(int32_t));
        for (size_t t = 0; t < top; t++) sub[t] = loop_verts[stack[t]];
        out_sub[n_sub] = sub; out_len[n_sub] = top; n_sub++;
    }
    return n_sub;
}

/* ------------------------------------------------------------------ */
/* HoleFill_process[_ex] — public entry point                          */
/*                                                                     */
/* interior_only != 0: fill ONLY geometrically-interior holes (the     */
/* surface surrounds them), leaving outer perimeters and open bays.    */
/* interior_only == 0: fill every cleanly-fillable closed 4+ loop (the */
/* per-cube behaviour; HoleFill_process wraps this).                   */
/* ------------------------------------------------------------------ */
int HoleFill_process_ex(Arena_T arena,
                     float **verts, int32_t **faces,
                     size_t *nv, size_t *nf,
                     const int cube_shape[3],
                     int interior_only,
                     size_t *out_n_loops,
                     size_t *out_n_interior,
                     size_t *out_n_filled)
{
    assert(arena && verts && *verts && faces && *faces && nv && nf);
    (void)cube_shape;   /* interior test is geometric (below); cube_shape unused */

    if (out_n_loops)    *out_n_loops = 0;
    if (out_n_interior) *out_n_interior = 0;
    if (out_n_filled)   *out_n_filled = 0;
    if (*nv < 3 || *nf == 0) return 0;

    Arena_Mark mark = Arena_save(arena);

    /* 0. Collapse duplicate vertices (same 3D position).
     * Poisson output can have co-located vertices that produce zero-length
     * boundary edges, crashing CDT.  We merge them in the face array and
     * remove any resulting degenerate and duplicate faces. */
    size_t n_loops_pre_dedup = 0;
    {
        const float *v = *verts;
        int32_t *f = *faces;
        size_t fnv = *nv;
        size_t fnf = *nf;

        /* Build a merge map: merge[i] = canonical vertex for vertex i.
         * Scan all pairs of boundary-adjacent vertices.  For speed, only
         * check boundary edges (non-boundary vertices can't create the
         * degenerate loops that hurt CDT). */
        UEdge *pre_edges = NULL;
        size_t pre_n = find_boundary_edges(arena, f, fnf, &pre_edges);

        /* Count loops before dedup for consistency check */
        {
            BoundaryLoop *pre_loops = NULL;
            n_loops_pre_dedup = chain_into_loops(arena, pre_edges, pre_n,
                                                  fnv, v, &pre_loops);
        }
        int32_t *merge = (int32_t *)ARENA_ALLOC(arena,
                           (long)fnv * (long)sizeof(int32_t));
        for (size_t i = 0; i < fnv; i++) merge[i] = (int32_t)i;

        size_t n_merged = 0;
        for (size_t ei = 0; ei < pre_n; ei++) {
            int32_t a = pre_edges[ei].v0;
            int32_t b = pre_edges[ei].v1;
            if (v[a*3+0] == v[b*3+0] &&
                v[a*3+1] == v[b*3+1] &&
                v[a*3+2] == v[b*3+2]) {
                /* Merge b into a (a < b since UEdge stores sorted) */
                merge[b] = a;
                n_merged++;
            }
        }

        if (n_merged > 0) {
            /* Chase merge chains: merge[merge[b]] ... */
            for (size_t i = 0; i < fnv; i++) {
                while (merge[merge[i]] != merge[i])
                    merge[i] = merge[merge[i]];
            }

            /* Rewrite face indices */
            for (size_t fi = 0; fi < fnf * 3; fi++) {
                f[fi] = merge[f[fi]];
            }

            /* Remove degenerate faces (two or more identical vertices) */
            size_t dst = 0;
            size_t n_degen = 0;
            for (size_t fi = 0; fi < fnf; fi++) {
                int32_t a = f[fi * 3 + 0];
                int32_t b = f[fi * 3 + 1];
                int32_t c = f[fi * 3 + 2];
                if (a == b || b == c || a == c) { n_degen++; continue; }
                f[dst * 3 + 0] = a;
                f[dst * 3 + 1] = b;
                f[dst * 3 + 2] = c;
                dst++;
            }

            /* Remove duplicate faces created by the merge.
             * After merging B→A, face (B,C,D) becomes (A,C,D) which
             * may already exist.  Sort faces by canonical vertex triple
             * (min, mid, max) and remove exact duplicates. */
            size_t n_dup = 0;
            {
                /* Build sortable keys: for each face, store directed canonical
                 * form (min_v, next_in_cycle, prev_in_cycle, orig_idx).
                 * Rotate so minimum vertex is first but PRESERVE cyclic
                 * winding order.  This way, faces with opposite winding
                 * (e.g. (A,B,C) vs (A,C,B)) are NOT treated as duplicates
                 * — only true duplicates (same cyclic order) are removed. */
                int32_t *keys = (int32_t *)ARENA_ALLOC(arena,
                                  (long)dst * 4 * (long)sizeof(int32_t));
                for (size_t fi = 0; fi < dst; fi++) {
                    int32_t a = f[fi * 3 + 0];
                    int32_t b = f[fi * 3 + 1];
                    int32_t c = f[fi * 3 + 2];
                    /* Rotate so minimum vertex is first */
                    if (b <= a && b <= c) {
                        int32_t t = a; a = b; b = c; c = t;
                    } else if (c <= a && c <= b) {
                        int32_t t = c; c = b; b = a; a = t;
                    }
                    keys[fi * 4 + 0] = a;
                    keys[fi * 4 + 1] = b;
                    keys[fi * 4 + 2] = c;
                    keys[fi * 4 + 3] = (int32_t)fi;
                }
                /* Sort by (a, b, c) using inline comparator on 4-int32 keys.
                 * cmp_uedge compares first two int32s — extend for third. */
                qsort(keys, dst, 4 * sizeof(int32_t), cmp_face_key);

                /* Mark duplicates in a flag array */
                uint8_t *keep = (uint8_t *)ARENA_ALLOC(arena,
                                  (long)dst * (long)sizeof(uint8_t));
                memset(keep, 1, dst);
                for (size_t fi = 1; fi < dst; fi++) {
                    if (keys[fi * 4 + 0] == keys[(fi-1) * 4 + 0] &&
                        keys[fi * 4 + 1] == keys[(fi-1) * 4 + 1] &&
                        keys[fi * 4 + 2] == keys[(fi-1) * 4 + 2]) {
                        /* Duplicate — mark for removal */
                        keep[keys[fi * 4 + 3]] = 0;
                        n_dup++;
                    }
                }

                if (n_dup > 0) {
                    size_t dst2 = 0;
                    for (size_t fi = 0; fi < dst; fi++) {
                        if (!keep[fi]) continue;
                        f[dst2 * 3 + 0] = f[fi * 3 + 0];
                        f[dst2 * 3 + 1] = f[fi * 3 + 1];
                        f[dst2 * 3 + 2] = f[fi * 3 + 2];
                        dst2++;
                    }
                    dst = dst2;
                }
            }

            fprintf(stderr, "  [hole_fill] collapsed %zu duplicate vertex pairs, "
                    "removed %zu degenerate + %zu duplicate faces (%zu -> %zu)\n",
                    n_merged, n_degen, n_dup, fnf, dst);
            *nf = dst;
        }
    }

    /* 0b. Repair orientation: BFS propagation to fix any flipped faces
     * from Poisson reconstruction output. */
    repair_mesh_orientation(arena, *faces, *nf);

    /* 1. Find boundary edges */
    UEdge *bdry_edges = NULL;
    size_t n_bdry = find_boundary_edges(arena, *faces, *nf, &bdry_edges);

    if (n_bdry == 0) {
        Arena_restore(arena, mark);
        return 0; /* watertight mesh */
    }

    /* 2. Chain into loops */
    BoundaryLoop *loops = NULL;
    size_t n_loops = chain_into_loops(arena, bdry_edges, n_bdry, *nv, *verts, &loops);

    /* Verify dedup preserved loop count (topology invariant) */
    if (n_loops_pre_dedup > 0 && n_loops != n_loops_pre_dedup) {
        fprintf(stderr, "  [hole_fill] WARNING: dedup changed loop count "
                "%zu -> %zu (topology corruption!)\n",
                n_loops_pre_dedup, n_loops);
    }

    if (n_loops == 0) {
        Arena_restore(arena, mark);
        return 0;
    }

    /* Find largest loop (outer boundary) — skip it */
    size_t largest_idx = 0;
    size_t largest_len = 0;
    for (size_t i = 0; i < n_loops; i++) {
        if (loops[i].len > largest_len) {
            largest_len = loops[i].len;
            largest_idx = i;
        }
    }

    /* 3. Classify and fill interior holes. A pinched loop splits into several
     * sub-loops, each its own fill, so size for more than n_loops (guarded
     * below as a hard backstop against overflow). */
    size_t max_fills = n_loops * 4 + 64;
    HoleFillResult *fills = (HoleFillResult *)ARENA_ALLOC(arena,
                              (long)max_fills * (long)sizeof(HoleFillResult));
    size_t n_fills = 0;
    size_t n_interior = 0;
    size_t n_micro = 0;
    double t_hf0 = ves_clock_sec();      /* per-hole timing: flag slow / hung holes */
    /* Per-hole OBJ dumps (debug). grid_weld --dump-stages sets these so each
     * hole's input boundary is dumped before the fill and its patch after --
     * a hang leaves the culprit's polygon as the last *_in.obj on disk. */
    const char *hole_dump_dir = sf_env("SEAM_HOLE_DUMP_DIR");
    const char *hole_dump_prefix = sf_env("SEAM_HOLE_DUMP_PREFIX");
    if (!hole_dump_prefix) hole_dump_prefix = "weld";

    /* interior_only: directed boundary half-edges for the signed-area test. */
    DirBHE *bhe = NULL; size_t n_bhe = 0;
    if (interior_only) n_bhe = build_dir_boundary(arena, *faces, *nf, &bhe);

    fprintf(stderr, "  [hole_fill] %zu boundary loops (largest=%zu with %zu verts)\n",
            n_loops, largest_idx, largest_len);
    if (hole_dump_dir)
        fprintf(stderr, "  [hole_fill] per-hole OBJ dumps -> %s\n", hole_dump_dir);
    HFFLUSH();

    for (size_t i = 0; i < n_loops; i++) {
        RunCtx_check();   /* one poll per hole (CDT fills can be slow) */
        /* Never fill the single largest loop: it is an outer perimeter, not a
         * hole (the giant loop is "most likely the outer one"). This applies in
         * BOTH modes. In interior_only mode the signed-area winding test (below)
         * does the heavy lifting -- it rejects EVERY component's perimeter and
         * every open bay, which is essential for a multi-component mesh where
         * only ONE of the N perimeters is the global largest -- and this
         * single-largest skip is a cheap belt-and-suspenders guard on top. */
        if (i == largest_idx) continue;

        size_t parent_n = loops[i].len - 1; /* exclude closing vertex */
        /* 3-loops are PinholeFill's job (the exact single-triangle path); this
         * CDT/Liepa path owns 4+ closed loops only. */
        if (parent_n < 4) continue;

        /* Retopologize a self-pinching "<><>" loop into simple sub-cycles. A
         * weld can fuse two cubes' boundaries into a figure-8 that passes
         * through one POSITION twice (same vertex index, or two coincident
         * distinct indices). That non-simple polygon hangs Triangle's CDT (the
         * 100-cube weld hung exactly here). Split at each pinch and fill each
         * simple sub-loop. A non-pinched loop yields one sub-loop aliasing the
         * original (no copy), so the common case is unchanged. */
        int32_t *subs[MAX_SUBLOOPS];
        size_t   sublen[MAX_SUBLOOPS];
        size_t   n_subs = split_pinched_loop(arena, loops[i].verts, parent_n,
                                             *verts, PINCH_EPS,
                                             subs, sublen, MAX_SUBLOOPS);
        if (n_subs > 1)
            HFLOG("  [hole_fill] loop %zu: %zu verts -> PINCHED, "
                    "retopologized into %zu sub-loops\n", i, parent_n, n_subs);

        for (size_t s = 0; s < n_subs; s++) {
            int32_t *lv = subs[s];
            size_t loop_n = sublen[s];
            /* A sub-loop of < 4 verts (e.g. a 3-cycle peeled at a pinch) is the
             * pinhole/fan path's job, not CDT. */
            if (loop_n < 4) continue;

            /* Size backstops apply ONLY to the per-cube fill-everything path
             * (interior_only == 0). That path has no geometric interior test,
             * so a vertex-count cap and a no-merger diameter cap stop it from
             * capping a wrap-spanning loop. The component/grid path
             * (interior_only != 0) does NOT size-gate: the single-largest skip
             * (above) drops the outer perimeter and the signed-area winding test
             * (below) drops every other perimeter / open bay, so a genuine LARGE
             * interior hole -- the giant comp-007 hole the size cap used to leave
             * gaping -- is filled rather than skipped. The diameter test is also
             * O(loop_n^2), so skipping it here keeps big interior holes cheap. */
            if (!interior_only) {
                /* Skip large holes (vertex-count cap) */
                if (loop_n > (size_t)MAX_LOOP_VERTS) {
                    HFLOG("  [hole_fill] loop %zu.%zu: %zu verts -> SKIP (too large, max=%d)\n",
                            i, s, loop_n, MAX_LOOP_VERTS);
                    continue;
                }

                /* No-merger diameter backstop. */
                {
                    const float *vp = *verts;
                    float dmax2 = 0.0f;
                    for (size_t a = 0; a < loop_n; a++) {
                        const float *pa = &vp[(size_t)lv[a] * 3];
                        for (size_t b = a + 1; b < loop_n; b++) {
                            const float *pb = &vp[(size_t)lv[b] * 3];
                            float dz = pa[0]-pb[0], dy = pa[1]-pb[1], dx = pa[2]-pb[2];
                            float d2 = dz*dz + dy*dy + dx*dx;
                            if (d2 > dmax2) dmax2 = d2;
                        }
                    }
                    if (dmax2 > HOLEFILL_MAX_DIAM_VOX * HOLEFILL_MAX_DIAM_VOX) {
                        HFLOG("  [hole_fill] loop %zu.%zu: %zu verts -> SKIP "
                                "(diameter %.1f > %.1f)\n", i, s, loop_n,
                                (double)sqrtf(dmax2), (double)HOLEFILL_MAX_DIAM_VOX);
                        continue;
                    }
                }
            }

            /* interior_only: fill only loops the surface surrounds. Leaves outer
             * perimeters and still-open bays. Geometric, world-coord agnostic. */
            if (interior_only &&
                !loop_is_interior(*verts, *faces, bhe, n_bhe, lv, loop_n)) {
                HFLOG("  [hole_fill] loop %zu.%zu: %zu verts -> SKIP "
                        "(exterior: perimeter / open bay)\n", i, s, loop_n);
                continue;
            }

            if (n_fills >= max_fills) {     /* hard backstop -- never expected */
                fprintf(stderr, "  [hole_fill] WARNING: fill cap %zu reached, "
                        "stopping\n", max_fills);
                break;
            }
            n_interior++;

            /* Fill this (now simple) sub-loop. The "-> filling" line AND the
             * PRE-fill OBJ dump are flushed BEFORE the call, so if fill_one_hole
             * still hangs (the watchdog in safe_triangulate bounds it) this log
             * line + the hole's *_in.obj are the LAST artifacts -- the culprit. */
            HFLOG("  [hole_fill] loop %zu.%zu: %zu verts -> filling (interior #%zu)...\n",
                    i, s, loop_n, n_interior);
            HFFLUSH();
            dump_hole_obj(hole_dump_dir, hole_dump_prefix, n_interior, i, loop_n,
                          *verts, lv, NULL /* PRE: input boundary */);
            log_loop_geom("geom", i, *verts, lv, loop_n);
            double t_hole0 = ves_clock_sec();
            HoleFillResult res;
            memset(&res, 0, sizeof(res));
            int rc = fill_one_hole(arena, *verts, *faces, *nf, lv, loop_n, &res);
            double t_hole_ms = (ves_clock_sec() - t_hole0) * 1000.0;
            if (rc == 0) {
                HFLOG("  [hole_fill] loop %zu.%zu: filled OK (%zu verts, %zu faces) [%.1f ms]\n",
                        i, s, res.nv, res.nf, t_hole_ms);
                HFFLUSH();
                dump_hole_obj(hole_dump_dir, hole_dump_prefix, n_interior, i, loop_n,
                              *verts, lv, &res /* POST: filled patch */);
                fills[n_fills] = res;
                n_fills++;
            } else {
                HFLOG("  [hole_fill] loop %zu.%zu: fill FAILED (%zu verts) [%.1f ms]\n",
                        i, s, loop_n, t_hole_ms);
                HFFLUSH();
            }
            if (t_hole_ms > 250.0) {
                HFLOG("  [hole_fill] *** SLOW HOLE: loop %zu.%zu (%zu verts) took %.1f ms ***\n",
                        i, s, loop_n, t_hole_ms);
                HFFLUSH();
            }
        }
    }

    fprintf(stderr, "  [hole_fill] %zu loops, %zu interior, %zu micro, "
            "%zu filled in %.2f s\n", n_loops, n_interior, n_micro, n_fills,
            ves_clock_sec() - t_hf0);
    HFFLUSH();

    /* 4. Stitch fills */
    if (n_fills > 0) {
        int rc = stitch_fills(arena, verts, faces, nv, nf,
                              fills, n_fills);
        if (rc != 0) {
            fprintf(stderr, "  [hole_fill] Stitch failed\n");
            Arena_restore(arena, mark);
            return -1;
        }
    }

    if (out_n_loops)    *out_n_loops    = n_loops;
    if (out_n_interior) *out_n_interior = n_interior;
    if (out_n_filled)   *out_n_filled   = n_fills;

    return 0;
}

/* Fill EVERY cleanly-fillable closed 4+ loop (per-cube behaviour). */
int HoleFill_process(Arena_T arena,
                     float **verts, int32_t **faces,
                     size_t *nv, size_t *nf,
                     const int cube_shape[3],
                     size_t *out_n_loops,
                     size_t *out_n_interior,
                     size_t *out_n_filled)
{
    return HoleFill_process_ex(arena, verts, faces, nv, nf, cube_shape,
                               0 /* interior_only */,
                               out_n_loops, out_n_interior, out_n_filled);
}
