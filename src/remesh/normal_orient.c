/*
 * normal_orient.c — Hoppe et al. 1992 §3.3 consistent tangent-plane
 * orientation. See normal_orient.h.
 *
 * Reuses the project KD-tree (neighbour queries) and union-find (Kruskal MST).
 * Transient graph/MST/BFS buffers use malloc/free (local to one call, same
 * convention as ball_pivot.c); the arena holds only the KD-tree + UF and is
 * restored on exit.
 */
#include "normal_orient.h"
#include "../common/run_ctx.h"

#include "../common/kdtree.h"
#include "../common/union_find.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int32_t a, b; float w; } Edge;

/* Sort by weight ascending; tie-break on (a,b) so the MST — and hence the
 * propagated orientation — is independent of KD-tree query order. */
static int edge_cmp(const void *pa, const void *pb)
{
    const Edge *ea = (const Edge *)pa, *eb = (const Edge *)pb;
    if (ea->w < eb->w) return -1;
    if (ea->w > eb->w) return 1;
    if (ea->a != eb->a) return (ea->a < eb->a) ? -1 : 1;
    if (ea->b != eb->b) return (ea->b < eb->b) ? -1 : 1;
    return 0;
}

/* --- Diagnostics (NORMAL_ORIENT_DEBUG) ----------------------------------- *
 * agree fraction = #{edges with n_i.n_j > 0} / #edges. Measures how sign-
 * consistent the field is over a given edge set. ~1.0 means "already
 * consistent"; BFS must drive the MST edge-set to ~1.0 if it is correct. */
static double agree_frac_edges(const Edge *edges, size_t ne, const float *normals)
{
    if (ne == 0) return 1.0;
    size_t pos = 0;
    for (size_t e = 0; e < ne; e++) {
        const float *na = &normals[(size_t)edges[e].a * 3];
        const float *nb = &normals[(size_t)edges[e].b * 3];
        float d = na[0]*nb[0] + na[1]*nb[1] + na[2]*nb[2];
        if (d > 0.0f) pos++;
    }
    return (double)pos / (double)ne;
}

static double agree_frac_csr(const int32_t *off, const int32_t *adj, int n,
                             const float *normals)
{
    size_t pos = 0, tot = 0;
    for (int p = 0; p < n; p++) {
        const float *np = &normals[(size_t)p * 3];
        for (int32_t e = off[p]; e < off[p+1]; e++) {
            int32_t ch = adj[e];
            if (ch <= p) continue;                  /* each undirected edge once */
            const float *nc = &normals[(size_t)ch * 3];
            float d = np[0]*nc[0] + np[1]*nc[1] + np[2]*nc[2];
            tot++;
            if (d > 0.0f) pos++;
        }
    }
    return tot ? (double)pos / (double)tot : 1.0;
}

int NormalOrient_consistent(Arena_T arena,
                            const float *verts,
                            float *normals,
                            size_t nv,
                            float radius,
                            size_t *out_flips)
{
    if (out_flips) *out_flips = 0;
    if (nv < 2) return 0;
    int n = (int)nv;
    int dbg = sf_env("NORMAL_ORIENT_DEBUG") != NULL;

    Arena_Mark mark = Arena_save(arena);
    KDTree_T tree = KDTree_new(arena, verts, nv);

    float r2 = radius * radius;
    enum { MAXNBR = 512 };
    int32_t *nbr = (int32_t *)malloc(MAXNBR * sizeof(int32_t));

    /* --- Riemannian graph: edges to neighbours within `radius`, j > i. --- */
    size_t ecap = (size_t)n * 16 + 16, ne = 0;
    Edge *edges = (Edge *)malloc(ecap * sizeof(Edge));
    size_t dbg_clamp = 0, dbg_iso = 0, dbg_degsum = 0;
    size_t dbg_degmin = (size_t)-1, dbg_degmax = 0;
    for (int i = 0; i < n; i++) {
        float c[3] = { verts[i*3+0], verts[i*3+1], verts[i*3+2] };
        size_t m = KDTree_ball_query(tree, c, r2, nbr, MAXNBR);
        if (dbg) {
            size_t raw = m;                    /* true count, may exceed buffer */
            if (raw > MAXNBR) dbg_clamp++;
            if (raw <= 1) dbg_iso++;           /* only itself within radius */
            dbg_degsum += raw;
            if (raw < dbg_degmin) dbg_degmin = raw;
            if (raw > dbg_degmax) dbg_degmax = raw;
        }
        if (m > MAXNBR) m = MAXNBR;            /* count may exceed buffer */
        const float *ni = &normals[i*3];
        for (size_t t = 0; t < m; t++) {
            int32_t j = nbr[t];
            if (j == i) continue;
            /* Store the edge canonically (lo,hi) and keep BOTH directions rather
             * than filtering j>i: KDTree_ball_query can return an INCOMPLETE
             * neighbour set (observed min-degree 3 on a dense grid where ~12 is
             * expected), so query(i) finding j does not imply query(j) finds i.
             * A j>i filter then silently drops that undirected edge and spuriously
             * fragments the graph (a lone vertex becomes its own MST component and
             * keeps a wrong sign). Taking the union of both directions' results
             * connects the edge whenever EITHER query sees it; Kruskal dedups the
             * (lo,hi) duplicates. */
            int32_t lo = (j < i) ? j : i, hi = (j < i) ? i : j;
            const float *nj = &normals[j*3];
            float d = ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2];
            float ad = (d < 0.0f) ? -d : d;
            if (ne >= ecap) { ecap *= 2; edges = (Edge *)realloc(edges, ecap * sizeof(Edge)); }
            edges[ne].a = lo; edges[ne].b = hi; edges[ne].w = 1.0f - ad; ne++;
        }
    }
    free(nbr);

    double dbg_agree_before = 0.0;
    if (dbg) {
        dbg_agree_before = agree_frac_edges(edges, ne, normals);
        fprintf(stderr,
            "[norient nv=%d r=%.3f] graph: edges=%zu  degree(incl.self) min=%zu max=%zu mean=%.1f"
            "  isolated=%zu  clamp(>%d)=%zu  AGREE_BEFORE=%.4f\n",
            n, (double)radius, ne,
            (dbg_degmin == (size_t)-1 ? 0 : dbg_degmin), dbg_degmax,
            n ? (double)dbg_degsum / (double)n : 0.0,
            dbg_iso, (int)MAXNBR, dbg_clamp, dbg_agree_before);
    }

    qsort(edges, ne, sizeof(Edge), edge_cmp);

    /* --- Kruskal MST. --- */
    UnionFind uf = UF_new(arena, n);
    size_t mcap = (size_t)(n > 1 ? n - 1 : 1);
    int32_t *ma = (int32_t *)malloc(mcap * sizeof(int32_t));
    int32_t *mb = (int32_t *)malloc(mcap * sizeof(int32_t));
    size_t nm = 0;
    double dbg_mst_wsum = 0.0; float dbg_mst_wmax = 0.0f;
    for (size_t e = 0; e < ne && (int)nm < n - 1; e++) {
        int32_t a = edges[e].a, b = edges[e].b;
        if (uf_find(&uf, a) != uf_find(&uf, b)) {
            uf_union(&uf, a, b);
            ma[nm] = a; mb[nm] = b; nm++;
            if (dbg) { dbg_mst_wsum += (double)edges[e].w;
                       if (edges[e].w > dbg_mst_wmax) dbg_mst_wmax = edges[e].w; }
        }
    }
    /* edges kept alive for AGREE_AFTER; freed at end. */

    if (dbg) {
        int comps = n - (int)nm;
        fprintf(stderr,
            "[norient nv=%d] MST: edges=%zu  components=%d  mst_w mean=%.4f max=%.4f"
            "  (graph_w min=%.4f med=%.4f max=%.4f)\n",
            n, nm, comps,
            nm ? dbg_mst_wsum / (double)nm : 0.0, (double)dbg_mst_wmax,
            ne ? (double)edges[0].w : 0.0,
            ne ? (double)edges[ne/2].w : 0.0,
            ne ? (double)edges[ne-1].w : 0.0);
        /* component-size histogram via UF roots */
        int32_t *csz = (int32_t *)calloc((size_t)n, sizeof(int32_t));
        for (int i = 0; i < n; i++) csz[uf_find(&uf, i)]++;
        size_t h1=0,h2=0,h3=0,h4=0,h5=0; int cmax=0;
        for (int i = 0; i < n; i++) {
            int s = csz[i];
            if (s == 0) continue;
            if (s > cmax) cmax = s;
            if (s == 1) h1++; else if (s <= 10) h2++; else if (s <= 100) h3++;
            else if (s <= 1000) h4++; else h5++;
        }
        fprintf(stderr,
            "[norient nv=%d] comp sizes: max=%d  buckets [1]=%zu [2-10]=%zu"
            " [11-100]=%zu [101-1000]=%zu [1000+]=%zu\n",
            n, cmax, h1, h2, h3, h4, h5);
        free(csz);
    }

    /* --- MST adjacency in CSR (undirected). --- */
    int32_t *deg = (int32_t *)calloc((size_t)n, sizeof(int32_t));
    for (size_t e = 0; e < nm; e++) { deg[ma[e]]++; deg[mb[e]]++; }
    int32_t *off = (int32_t *)malloc((size_t)(n + 1) * sizeof(int32_t));
    off[0] = 0;
    for (int i = 0; i < n; i++) off[i+1] = off[i] + deg[i];
    int32_t *adj = (int32_t *)malloc((size_t)off[n] * sizeof(int32_t));
    int32_t *cur = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    memcpy(cur, off, (size_t)n * sizeof(int32_t));
    for (size_t e = 0; e < nm; e++) {
        int32_t a = ma[e], b = mb[e];
        adj[cur[a]++] = b; adj[cur[b]++] = a;
    }
    free(cur); free(ma); free(mb); free(deg);

    /* --- Seed per component (max-z member, normal forced toward +z), then
     *     BFS the MST flipping any child whose normal opposes its parent. --- */
    int32_t *seed_of_root = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    for (int i = 0; i < n; i++) seed_of_root[i] = -1;
    for (int i = 0; i < n; i++) {
        int32_t r = uf_find(&uf, i);
        int32_t s = seed_of_root[r];
        /* coords are (z, y, x): index 0 is z */
        if (s < 0 || verts[i*3+0] > verts[s*3+0]) seed_of_root[r] = i;
    }

    /* Snapshot the original normals so each component's GLOBAL sign can be
     * re-anchored by majority vote after propagation. Without this, an already-
     * consistent component whose seed happens to point -z would be globally
     * negated (a no-op for consistency, but it perturbs BPA on the already-clean
     * components -- measured: boundary 908->1198 on one cube). */
    float *orig = (float *)malloc((size_t)n * 3 * sizeof(float));
    memcpy(orig, normals, (size_t)n * 3 * sizeof(float));

    uint8_t *vis = (uint8_t *)calloc((size_t)n, 1);
    int32_t *queue = (int32_t *)malloc((size_t)n * sizeof(int32_t));
    size_t seeds_total = 0, reanchored = 0;
    for (int i = 0; i < n; i++) {
        if (uf_find(&uf, i) != i) continue;          /* component roots only */
        int32_t seed = seed_of_root[i];
        if (seed < 0 || vis[seed]) continue;
        seeds_total++;
        /* BFS the MST, flipping any child whose normal opposes its parent. The
         * seed KEEPS its own sign here; the component's absolute sign is fixed
         * by the majority vote below, so the seed choice cannot matter. */
        size_t qh = 0, qt = 0;
        queue[qt++] = seed; vis[seed] = 1;
        while (qh < qt) {
            int32_t p = queue[qh++];
            float np0 = normals[p*3+0], np1 = normals[p*3+1], np2 = normals[p*3+2];
            for (int32_t e = off[p]; e < off[p+1]; e++) {
                int32_t ch = adj[e];
                if (vis[ch]) continue;
                vis[ch] = 1;
                float *nc = &normals[ch*3];
                float d = np0*nc[0] + np1*nc[1] + np2*nc[2];
                if (d < 0.0f) { nc[0] = -nc[0]; nc[1] = -nc[1]; nc[2] = -nc[2]; }
                queue[qt++] = ch;
            }
        }
        /* Re-anchor: the component is now internally consistent but its global
         * sign is arbitrary (seed-dependent). Pick the sign that AGREES with the
         * original MLS normals on the MAJORITY of the component's vertices, so an
         * already-consistent field is preserved (0 net flips) and a split field
         * flips only its minority region. */
        size_t agree = 0;
        for (size_t q = 0; q < qt; q++) {
            int32_t v = queue[q];
            float dd = normals[v*3+0]*orig[v*3+0] + normals[v*3+1]*orig[v*3+1]
                     + normals[v*3+2]*orig[v*3+2];
            if (dd > 0.0f) agree++;
        }
        if (agree * 2 < qt) {                        /* minority agrees -> negate */
            for (size_t q = 0; q < qt; q++) {
                int32_t v = queue[q];
                normals[v*3+0] = -normals[v*3+0];
                normals[v*3+1] = -normals[v*3+1];
                normals[v*3+2] = -normals[v*3+2];
            }
            reanchored++;
        }
    }
    /* out_flips = vertices whose FINAL normal differs from the original. */
    size_t flips = 0;
    for (int i = 0; i < n; i++) {
        float dd = normals[i*3+0]*orig[i*3+0] + normals[i*3+1]*orig[i*3+1]
                 + normals[i*3+2]*orig[i*3+2];
        if (dd < 0.0f) flips++;
    }
    free(orig);
    if (dbg) {
        double aft_all = agree_frac_edges(edges, ne, normals);
        double aft_mst = agree_frac_csr(off, adj, n, normals);
        const char *bug = (aft_mst < 0.999)
            ? " <<< MST NOT CONSISTENT -- BFS PROPAGATION BUG" : "";
        const char *verdict;
        if (aft_mst < 0.999)
            verdict = "BUG: BFS failed to make the tree consistent";
        else if (dbg_agree_before >= 0.999 && aft_all >= 0.999)
            verdict = "field was ALREADY consistent; Hoppe left it unchanged (majority re-anchor)";
        else if (aft_all > dbg_agree_before + 0.001)
            verdict = "Hoppe IMPROVED graph consistency";
        else if (aft_all < dbg_agree_before - 0.001)
            verdict = "Hoppe REDUCED graph consistency (fragmentation: independently-seeded components disagree)";
        else
            verdict = "no net change in graph consistency";
        fprintf(stderr,
            "[norient nv=%d] BFS: seeds=%zu reanchored=%zu total_flips=%zu (%.1f%%)"
            "  AGREE_AFTER all=%.4f MST=%.4f%s\n",
            n, seeds_total, reanchored, flips,
            n ? 100.0 * (double)flips / (double)n : 0.0,
            aft_all, aft_mst, bug);
        fprintf(stderr, "[norient nv=%d] VERDICT: %s\n", n, verdict);
    }
    free(edges);
    free(queue); free(vis); free(adj); free(off); free(seed_of_root);

    Arena_restore(arena, mark);
    if (out_flips) *out_flips = flips;
    return 0;
}
