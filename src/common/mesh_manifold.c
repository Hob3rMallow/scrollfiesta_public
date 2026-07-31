/* mesh_manifold.c -- combinatorial 2-manifold audit (see mesh_manifold.h).
 *
 * Algorithms ported from the validated src/tools/manifold_check.c, adapted to
 * arena scratch and the pipeline's int32_t face layout. Self-restoring: every
 * allocation is rolled back via Arena_save/Arena_restore before return. */
#include "mesh_manifold.h"

#include <stdlib.h>   /* qsort */
#include <string.h>   /* memset */
#include <assert.h>

/* ---- edge audit: directed half-edges, sorted by undirected key ---- */
typedef struct { int32_t src, dst; } DHE;

static int cmp_dhe_undirected(const void *pa, const void *pb) {
    const DHE *a = (const DHE *)pa, *b = (const DHE *)pb;
    int32_t a0 = (a->src < a->dst) ? a->src : a->dst;
    int32_t a1 = (a->src < a->dst) ? a->dst : a->src;
    int32_t b0 = (b->src < b->dst) ? b->src : b->dst;
    int32_t b1 = (b->src < b->dst) ? b->dst : b->src;
    if (a0 != b0) return (a0 < b0) ? -1 : 1;
    if (a1 != b1) return (a1 < b1) ? -1 : 1;
    return 0;
}

static void edge_audit(Arena_T arena, const int32_t *F, size_t nf,
                       MeshManifoldStats *out) {
    if (nf == 0) return;
    size_t hn = nf * 3;
    DHE *he = (DHE *)ARENA_ALLOC(arena, (long)(hn * sizeof(DHE)));
    for (size_t f = 0; f < nf; f++) {
        int32_t v0 = F[f*3+0], v1 = F[f*3+1], v2 = F[f*3+2];
        he[f*3+0].src = v0; he[f*3+0].dst = v1;
        he[f*3+1].src = v1; he[f*3+1].dst = v2;
        he[f*3+2].src = v2; he[f*3+2].dst = v0;
    }
    qsort(he, hn, sizeof(DHE), cmp_dhe_undirected);
    size_t i = 0;
    while (i < hn) {
        size_t j = i + 1;
        int32_t a0 = (he[i].src < he[i].dst) ? he[i].src : he[i].dst;
        int32_t a1 = (he[i].src < he[i].dst) ? he[i].dst : he[i].src;
        while (j < hn) {
            int32_t b0 = (he[j].src < he[j].dst) ? he[j].src : he[j].dst;
            int32_t b1 = (he[j].src < he[j].dst) ? he[j].dst : he[j].src;
            if (b0 != a0 || b1 != a1) break;
            j++;
        }
        size_t run = j - i;
        if (run == 1) { out->boundary_edges++; }
        else if (run > 2) { out->nm_edges++; }
        else {
            int dir0 = (he[i].src   < he[i].dst)   ? 0 : 1;
            int dir1 = (he[i+1].src < he[i+1].dst) ? 0 : 1;
            if (dir0 == dir1) out->same_dir_edges++;
            else out->manifold_edges++;
        }
        i = j;
    }
}

/* ---- vertex audit: per-vertex fan count via shared-neighbour union-find ---- */
typedef struct { int32_t nbr; int32_t loc; } NbrPair;
static int cmp_nbrpair(const void *a, const void *b) {
    const NbrPair *x = (const NbrPair *)a, *y = (const NbrPair *)b;
    return (x->nbr < y->nbr) ? -1 : (x->nbr > y->nbr ? 1 : 0);
}
static int32_t lf_find(int32_t *p, int32_t x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

static void vertex_audit(Arena_T arena, size_t nv, const int32_t *F, size_t nf,
                         MeshManifoldStats *out, uint8_t *out_mask) {
    if (nv == 0 || nf == 0) return;
    size_t *deg = (size_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(size_t));
    for (size_t f = 0; f < nf; f++)
        for (int k = 0; k < 3; k++) deg[F[f*3+(size_t)k]]++;

    size_t *off = (size_t *)ARENA_CALLOC(arena, (long)(nv + 1), (long)sizeof(size_t));
    size_t maxdeg = 0;
    for (size_t v = 0; v < nv; v++) {
        off[v+1] = off[v] + deg[v];
        if (deg[v] > maxdeg) maxdeg = deg[v];
    }
    size_t total = off[nv];                              /* == 3*nf */
    int32_t *cn0 = (int32_t *)ARENA_ALLOC(arena, (long)(total * sizeof(int32_t)));
    int32_t *cn1 = (int32_t *)ARENA_ALLOC(arena, (long)(total * sizeof(int32_t)));
    size_t  *cur = (size_t  *)ARENA_ALLOC(arena, (long)(nv * sizeof(size_t)));
    for (size_t v = 0; v < nv; v++) cur[v] = off[v];
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F[f*3+0], b = F[f*3+1], c = F[f*3+2];
        size_t pa = cur[a]++; cn0[pa] = b; cn1[pa] = c;
        size_t pb = cur[b]++; cn0[pb] = a; cn1[pb] = c;
        size_t pc = cur[c]++; cn0[pc] = a; cn1[pc] = b;
    }

    size_t cap = (maxdeg > 0 ? maxdeg : 1);
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena, (long)(cap * sizeof(int32_t)));
    NbrPair *pairs  = (NbrPair *)ARENA_ALLOC(arena, (long)(cap * 2 * sizeof(NbrPair)));
    for (size_t v = 0; v < nv; v++) {
        size_t s = off[v], e = off[v+1];
        size_t k = e - s;
        if (k <= 1) continue;                            /* 0/1 face -> one fan */
        for (size_t t = 0; t < k; t++) parent[t] = (int32_t)t;
        size_t np = 0;
        for (size_t t = 0; t < k; t++) {
            int32_t n0 = cn0[s+t], n1 = cn1[s+t];
            if (n0 != (int32_t)v) { pairs[np].nbr = n0; pairs[np].loc = (int32_t)t; np++; }
            if (n1 != (int32_t)v) { pairs[np].nbr = n1; pairs[np].loc = (int32_t)t; np++; }
        }
        qsort(pairs, np, sizeof(NbrPair), cmp_nbrpair);
        for (size_t t = 1; t < np; t++) {
            if (pairs[t].nbr == pairs[t-1].nbr) {
                int32_t ra = lf_find(parent, pairs[t].loc);
                int32_t rb = lf_find(parent, pairs[t-1].loc);
                if (ra != rb) parent[rb] = ra;
            }
        }
        size_t fans = 0;
        for (size_t t = 0; t < k; t++) if (lf_find(parent, (int32_t)t) == (int32_t)t) fans++;
        if (fans > 1) {
            out->nm_verts++;
            if (out_mask) out_mask[v] = 1;
        }
    }
}

MeshManifoldStats MeshManifold_audit(Arena_T arena, size_t nv,
                                     const int32_t *faces, size_t nf) {
    MeshManifoldStats s;
    memset(&s, 0, sizeof s);
    s.nv = nv; s.nf = nf;
    Arena_Mark mark = Arena_save(arena);
    edge_audit(arena, faces, nf, &s);
    vertex_audit(arena, nv, faces, nf, &s, NULL);
    Arena_restore(arena, mark);
    return s;
}

size_t MeshManifold_mark_nonmanifold_vertices(Arena_T arena, size_t nv,
                                              const int32_t *faces, size_t nf,
                                              uint8_t *out_mask) {
    MeshManifoldStats s;
    memset(&s, 0, sizeof s);
    if (!out_mask) return 0;
    memset(out_mask, 0, nv * sizeof(*out_mask));
    Arena_Mark mark = Arena_save(arena);
    vertex_audit(arena, nv, faces, nf, &s, out_mask);
    Arena_restore(arena, mark);
    return s.nm_verts;
}

/* ---- selftest: same five cases as src/tools/manifold_check.c ---- */
int MeshManifold_selftest(Arena_T arena) {
    int fails = 0;
    /* 1: closed, consistently-wound tetrahedron -> fully manifold */
    {
        int32_t F[] = { 0,2,1,  0,1,3,  0,3,2,  1,2,3 };
        MeshManifoldStats s = MeshManifold_audit(arena, 4, F, 4);
        if (!(s.boundary_edges==0 && s.manifold_edges==6 && s.same_dir_edges==0 &&
              s.nm_edges==0 && s.nm_verts==0 && MeshManifold_ok(&s))) fails++;
    }
    /* 2: two triangles sharing one edge -> manifold with boundary */
    {
        int32_t F[] = { 0,1,2,  1,3,2 };
        MeshManifoldStats s = MeshManifold_audit(arena, 4, F, 2);
        if (!(s.boundary_edges==4 && s.manifold_edges==1 && s.nm_edges==0 &&
              s.nm_verts==0 && MeshManifold_ok(&s))) fails++;
    }
    /* 3: non-manifold edge -- 3 triangles share edge {0,1} */
    {
        int32_t F[] = { 0,1,2,  0,1,3,  0,1,4 };
        MeshManifoldStats s = MeshManifold_audit(arena, 5, F, 3);
        if (!(s.nm_edges==1 && !MeshManifold_ok(&s))) fails++;
    }
    /* 4: bowtie vertex -- two tris share only vertex 0 (edge audit blind) */
    {
        int32_t F[] = { 0,1,2,  0,3,4 };
        MeshManifoldStats s = MeshManifold_audit(arena, 5, F, 2);
        if (!(s.nm_edges==0 && s.nm_verts==1 && !MeshManifold_ok(&s))) fails++;
    }
    /* 5: doubled face -> 3 same_dir edges, still edge+vertex manifold */
    {
        int32_t F[] = { 0,1,2,  0,1,2 };
        MeshManifoldStats s = MeshManifold_audit(arena, 3, F, 2);
        if (!(s.same_dir_edges==3 && s.nm_edges==0 && s.nm_verts==0)) fails++;
    }
    return fails;
}
