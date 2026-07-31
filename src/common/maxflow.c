/* maxflow.c -- Dinic max-flow. See maxflow.h.
 *
 * Adjacency as paired arcs in one flat array: arc i and arc i^1 are each
 * other's reverse (the classic trick), heads/caps in growable arrays, then
 * bucketed into CSR at solve time. Level BFS + blocking-flow DFS with the
 * current-arc optimization. Everything int64; sizes here are tiny (~1e4
 * nodes) so simplicity beats cleverness. */
#include "maxflow.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Maxflow_T {
    Arena_T arena;
    int32_t n;            /* problem nodes; source = n, sink = n+1 */
    int32_t nn;           /* n + 2 */
    /* arc store (paired: arc a's reverse is a^1) */
    int32_t *head;        /* [cap_arcs] target node of arc */
    int32_t *tail;        /* [cap_arcs] source node of arc (for CSR build) */
    int64_t *cap;         /* [cap_arcs] residual capacity */
    size_t   n_arcs, cap_arcs;
    /* CSR built at solve time */
    int32_t *aoff;        /* [nn+1] */
    int32_t *alist;       /* [n_arcs] arc ids sorted by tail */
    int32_t *level;       /* [nn] */
    int32_t *cur;         /* [nn] current-arc iterator */
    int32_t *queue;       /* [nn] BFS */
    int      solved;
};

static void mf_grow(Maxflow_T mf, size_t need)
{
    if (mf->n_arcs + need <= mf->cap_arcs)
        return;
    size_t nc = mf->cap_arcs * 2;
    while (nc < mf->n_arcs + need) nc *= 2;
    int32_t *nh = (int32_t *)ARENA_ALLOC(mf->arena, nc * sizeof(int32_t));
    int32_t *nt = (int32_t *)ARENA_ALLOC(mf->arena, nc * sizeof(int32_t));
    int64_t *ncp = (int64_t *)ARENA_ALLOC(mf->arena, nc * sizeof(int64_t));
    memcpy(nh, mf->head, mf->n_arcs * sizeof(int32_t));
    memcpy(nt, mf->tail, mf->n_arcs * sizeof(int32_t));
    memcpy(ncp, mf->cap, mf->n_arcs * sizeof(int64_t));
    mf->head = nh;
    mf->tail = nt;
    mf->cap = ncp;
    mf->cap_arcs = nc;
}

Maxflow_T Maxflow_new(Arena_T arena, int32_t n, size_t edge_hint)
{
    assert(arena);
    assert(n >= 1);
    Maxflow_T mf = (Maxflow_T)ARENA_CALLOC(arena, 1, sizeof(*mf));
    mf->arena = arena;
    mf->n = n;
    mf->nn = n + 2;
    mf->cap_arcs = 4 * (edge_hint > 16 ? edge_hint : 16) + 4 * (size_t)n;
    mf->head = (int32_t *)ARENA_ALLOC(arena, mf->cap_arcs * sizeof(int32_t));
    mf->tail = (int32_t *)ARENA_ALLOC(arena, mf->cap_arcs * sizeof(int32_t));
    mf->cap = (int64_t *)ARENA_ALLOC(arena, mf->cap_arcs * sizeof(int64_t));
    mf->n_arcs = 0;
    mf->solved = 0;
    return mf;
}

static void mf_arc_pair(Maxflow_T mf, int32_t u, int32_t v,
                        int64_t c_uv, int64_t c_vu)
{
    mf_grow(mf, 2);
    mf->tail[mf->n_arcs] = u;
    mf->head[mf->n_arcs] = v;
    mf->cap[mf->n_arcs] = c_uv;
    mf->n_arcs++;
    mf->tail[mf->n_arcs] = v;
    mf->head[mf->n_arcs] = u;
    mf->cap[mf->n_arcs] = c_vu;
    mf->n_arcs++;
}

void Maxflow_add_edge(Maxflow_T mf, int32_t u, int32_t v,
                      int64_t cap, int64_t rcap)
{
    assert(mf);
    assert(u >= 0 && u < mf->n && v >= 0 && v < mf->n && u != v);
    assert(cap >= 0 && rcap >= 0);
    assert(!mf->solved);
    mf_arc_pair(mf, u, v, cap, rcap);
}

void Maxflow_add_terminal(Maxflow_T mf, int32_t v, int64_t cap_s,
                          int64_t cap_t)
{
    assert(mf);
    assert(v >= 0 && v < mf->n);
    assert(cap_s >= 0 && cap_t >= 0);
    assert(!mf->solved);
    if (cap_s > 0)
        mf_arc_pair(mf, mf->n, v, cap_s, 0);
    if (cap_t > 0)
        mf_arc_pair(mf, v, mf->n + 1, cap_t, 0);
}

static int mf_bfs(Maxflow_T mf)
{
    for (int32_t i = 0; i < mf->nn; i++) mf->level[i] = -1;
    int32_t qh = 0, qt = 0;
    mf->queue[qt++] = mf->n;   /* source */
    mf->level[mf->n] = 0;
    while (qh < qt) {
        int32_t u = mf->queue[qh++];
        for (int32_t e = mf->aoff[u]; e < mf->aoff[u + 1]; e++) {
            int32_t a = mf->alist[e];
            int32_t v = mf->head[a];
            if (mf->cap[a] > 0 && mf->level[v] < 0) {
                mf->level[v] = mf->level[u] + 1;
                mf->queue[qt++] = v;
            }
        }
    }
    return mf->level[mf->n + 1] >= 0;
}

/* one augmenting path per call, iterative (explicit path stack -- a level
 * path can be O(n) long, and 1 MB of Windows stack is not a place to bet) */
static int64_t mf_augment(Maxflow_T mf, int32_t *path)
{
    int32_t sink = mf->n + 1;
    int32_t depth = 0;
    int32_t u = mf->n;
    for (;;) {
        if (u == sink) {
            /* bottleneck along path[0..depth-1] */
            int64_t bn = INT64_MAX;
            for (int32_t i = 0; i < depth; i++)
                if (mf->cap[path[i]] < bn) bn = mf->cap[path[i]];
            for (int32_t i = 0; i < depth; i++) {
                mf->cap[path[i]] -= bn;
                mf->cap[path[i] ^ 1] += bn;
            }
            return bn;
        }
        int advanced = 0;
        for (; mf->cur[u] < mf->aoff[u + 1]; mf->cur[u]++) {
            int32_t a = mf->alist[mf->cur[u]];
            int32_t v = mf->head[a];
            if (mf->cap[a] > 0 && mf->level[v] == mf->level[u] + 1) {
                path[depth++] = a;
                u = v;
                advanced = 1;
                break;
            }
        }
        if (advanced)
            continue;
        /* dead end: retreat (or done if at source) */
        mf->level[u] = -1;   /* prune */
        if (depth == 0)
            return 0;
        depth--;
        u = mf->tail[path[depth]];
        mf->cur[u]++;
    }
}

int64_t Maxflow_solve(Maxflow_T mf)
{
    assert(mf);
    assert(!mf->solved);
    Arena_T a = mf->arena;
    /* CSR over arcs by tail */
    mf->aoff = (int32_t *)ARENA_CALLOC(a, (size_t)mf->nn + 1,
                                       sizeof(int32_t));
    mf->alist = (int32_t *)ARENA_ALLOC(a, (mf->n_arcs + 1)
                                       * sizeof(int32_t));
    mf->level = (int32_t *)ARENA_ALLOC(a, (size_t)mf->nn * sizeof(int32_t));
    mf->cur = (int32_t *)ARENA_ALLOC(a, (size_t)mf->nn * sizeof(int32_t));
    mf->queue = (int32_t *)ARENA_ALLOC(a, (size_t)mf->nn * sizeof(int32_t));
    for (size_t i = 0; i < mf->n_arcs; i++)
        mf->aoff[mf->tail[i] + 1]++;
    for (int32_t v = 0; v < mf->nn; v++)
        mf->aoff[v + 1] += mf->aoff[v];
    {
        int32_t *fill = (int32_t *)ARENA_ALLOC(a, (size_t)mf->nn
                                               * sizeof(int32_t));
        memcpy(fill, mf->aoff, (size_t)mf->nn * sizeof(int32_t));
        for (size_t i = 0; i < mf->n_arcs; i++)
            mf->alist[fill[mf->tail[i]]++] = (int32_t)i;
    }

    int64_t flow = 0;
    int32_t *path = (int32_t *)ARENA_ALLOC(a, (size_t)mf->nn
                                           * sizeof(int32_t));
    while (mf_bfs(mf)) {
        memcpy(mf->cur, mf->aoff, (size_t)mf->nn * sizeof(int32_t));
        int64_t got = 0;
        do {
            got = mf_augment(mf, path);
            flow += got;
        } while (got > 0);
    }
    /* final BFS marks the source side (level >= 0 == residual-reachable) */
    (void)mf_bfs(mf);
    mf->solved = 1;
    return flow;
}

int Maxflow_in_source_side(const Maxflow_T mf, int32_t v)
{
    assert(mf);
    assert(mf->solved);
    assert(v >= 0 && v < mf->n);
    return mf->level[v] >= 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

static int mfst_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[maxflow selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* deterministic LCG so the fuzz cases are reproducible */
static uint32_t mfst_rng(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return *s;
}

/* Solve a general submodular binary energy exactly via KZ reduction:
 *   E(x) = sum_i U[i][x_i] + sum_e T[e][x_u][x_v],  x in {0,1}^n
 * with T[e][0][1] + T[e][1][0] >= T[e][0][0] + T[e][1][1] (submodular).
 * Label 1 = sink side. Returns min energy; labels into out_x (may be NULL). */
static int64_t mfst_solve_binary(Arena_T arena, int n,
                                 const int64_t (*U)[2],
                                 int ne, const int32_t (*euv)[2],
                                 const int64_t (*T)[2][2], int8_t *out_x)
{
    Maxflow_T mf = Maxflow_new(arena, n, (size_t)ne + 1);
    int64_t constant = 0;
    /* unary accumulation D[i] = cost(1) - cost(0), plus KZ pairwise terms */
    int64_t *D = (int64_t *)ARENA_CALLOC(arena, (size_t)n, sizeof(int64_t));
    for (int i = 0; i < n; i++) {
        constant += U[i][0];
        D[i] += U[i][1] - U[i][0];
    }
    for (int e = 0; e < ne; e++) {
        int32_t u = euv[e][0], v = euv[e][1];
        int64_t t00 = T[e][0][0], t01 = T[e][0][1];
        int64_t t10 = T[e][1][0], t11 = T[e][1][1];
        /* KZ: E = t00 + (t10-t00)[x_u] + (t11-t10)[x_v]
         *       + (t01+t10-t00-t11)[x_u=0,x_v=1] */
        constant += t00;
        D[u] += t10 - t00;
        D[v] += t11 - t10;
        int64_t beta = t01 + t10 - t00 - t11;
        assert(beta >= 0);   /* submodularity */
        if (beta > 0)
            Maxflow_add_edge(mf, u, v, beta, 0);
    }
    for (int i = 0; i < n; i++) {
        if (D[i] > 0)
            Maxflow_add_terminal(mf, i, D[i], 0);   /* pay D to take label 1 */
        else if (D[i] < 0) {
            constant += D[i];
            Maxflow_add_terminal(mf, i, 0, -D[i]);  /* pay -D to stay 0 */
        }
    }
    int64_t cut = Maxflow_solve(mf);
    if (out_x != NULL)
        for (int i = 0; i < n; i++)
            out_x[i] = (int8_t)(Maxflow_in_source_side(mf, i) ? 0 : 1);
    return constant + cut;
}

/* brute-force the same energy */
static int64_t mfst_brute(int n, const int64_t (*U)[2], int ne,
                          const int32_t (*euv)[2], const int64_t (*T)[2][2])
{
    int64_t best = INT64_MAX;
    for (uint32_t m = 0; m < (1u << n); m++) {
        int64_t e = 0;
        for (int i = 0; i < n; i++) e += U[i][(m >> i) & 1];
        for (int k = 0; k < ne; k++)
            e += T[k][(m >> euv[k][0]) & 1][(m >> euv[k][1]) & 1];
        if (e < best) best = e;
    }
    return best;
}

int Maxflow_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* t1: classic hand graph -- two disjoint 2-hop paths s->a->t (5,4) and
     * s->b->t (3,7): flow = min(5,4)+min(3,7) = 7; a on sink side of its
     * bottleneck, b on source side of its bottleneck */
    {
        Maxflow_T mf = Maxflow_new(arena, 2, 4);
        Maxflow_add_terminal(mf, 0, 5, 4);
        Maxflow_add_terminal(mf, 1, 3, 7);
        int64_t f = Maxflow_solve(mf);
        mfst_check(f == 7, "t1 flow value", &fails);
        mfst_check(Maxflow_in_source_side(mf, 0) == 1
                   && Maxflow_in_source_side(mf, 1) == 0,
                   "t1 cut sides", &fails);
    }

    /* t2: pairwise arc bottleneck -- s->0 (10), 0->1 (3), 1->t (10): flow 3,
     * cut crosses the middle arc, node 0 source side, node 1 sink side...
     * (sink side = NOT residual-reachable) */
    {
        Maxflow_T mf = Maxflow_new(arena, 2, 4);
        Maxflow_add_terminal(mf, 0, 10, 0);
        Maxflow_add_terminal(mf, 1, 0, 10);
        Maxflow_add_edge(mf, 0, 1, 3, 0);
        int64_t f = Maxflow_solve(mf);
        mfst_check(f == 3, "t2 flow value", &fails);
        mfst_check(Maxflow_in_source_side(mf, 0) == 1
                   && Maxflow_in_source_side(mf, 1) == 0,
                   "t2 cut through pairwise arc", &fails);
    }

    /* t3: fuzz -- random submodular binary energies on n<=10 nodes, KZ+Dinic
     * vs exhaustive enumeration. This is the correctness gate for the
     * collective-shift moves. */
    {
        uint32_t seed = 0xC0FFEEu;
        int bad = 0;
        for (int trial = 0; trial < 200 && !bad; trial++) {
            Arena_Mark mark = Arena_save(arena);
            int n = 2 + (int)(mfst_rng(&seed) % 9);          /* 2..10 */
            int ne_max = n * (n - 1) / 2;
            int ne = 1 + (int)(mfst_rng(&seed) % (uint32_t)ne_max);
            int64_t (*U)[2] = (int64_t (*)[2])ARENA_ALLOC(
                arena, (size_t)n * sizeof(*U));
            int32_t (*euv)[2] = (int32_t (*)[2])ARENA_ALLOC(
                arena, (size_t)ne * sizeof(*euv));
            int64_t (*T)[2][2] = (int64_t (*)[2][2])ARENA_ALLOC(
                arena, (size_t)ne * sizeof(*T));
            for (int i = 0; i < n; i++) {
                U[i][0] = (int64_t)(mfst_rng(&seed) % 1000);
                U[i][1] = (int64_t)(mfst_rng(&seed) % 1000);
            }
            for (int e = 0; e < ne; e++) {
                int32_t u = (int32_t)(mfst_rng(&seed) % (uint32_t)n);
                int32_t v = (int32_t)(mfst_rng(&seed) % (uint32_t)n);
                if (v == u) v = (v + 1) % n;
                euv[e][0] = u;
                euv[e][1] = v;
                /* random submodular table: pick t00,t11 and inflate the
                 * anti-diagonal so t01+t10 >= t00+t11 */
                int64_t t00 = (int64_t)(mfst_rng(&seed) % 500);
                int64_t t11 = (int64_t)(mfst_rng(&seed) % 500);
                int64_t extra = (int64_t)(mfst_rng(&seed) % 500);
                int64_t t01 = (t00 + t11 + extra) / 2;
                int64_t t10 = t00 + t11 + extra - t01;
                T[e][0][0] = t00;
                T[e][0][1] = t01;
                T[e][1][0] = t10;
                T[e][1][1] = t11;
            }
            int64_t got = mfst_solve_binary(arena, n, U, ne, euv, T, NULL);
            int64_t want = mfst_brute(n, U, ne, euv, T);
            if (got != want) {
                fprintf(stderr, "[maxflow selftest]   FAIL: fuzz trial %d "
                        "(n=%d ne=%d): kz=%lld brute=%lld\n", trial, n, ne,
                        (long long)got, (long long)want);
                bad = 1;
            }
            Arena_restore(arena, mark);
        }
        mfst_check(!bad, "t3 200x brute-force parity", &fails);
    }

    /* t4: the actual move shape -- residual chain: 4 nodes, edges with
     * r={0,0,-1} (a branch cut at the last edge), delta=+1: shifting the
     * whole chain but the last node fixes r3 at cost of nothing else if the
     * last node is anchored... encode theta and verify the optimal S found
     * by KZ matches brute force including labels */
    {
        enum { N = 4, NE = 3 };
        int64_t U[N][2];
        int32_t euv[NE][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 } };
        int64_t T[NE][2][2];
        int64_t r[NE] = { 0, 0, -1 };
        int64_t conf[NE] = { 100, 100, 5 };
        int delta = 1;
        for (int i = 0; i < N; i++) { U[i][0] = 0; U[i][1] = 0; }
        U[3][1] = 1000;   /* node 3 anchored at label 0 (do not shift) */
        for (int e = 0; e < NE; e++) {
            int64_t c = conf[e];
            int64_t rr = r[e];
#define MFST_ABS(x) ((x) < 0 ? -(x) : (x))
            T[e][0][0] = c * MFST_ABS(rr);
            T[e][1][1] = c * MFST_ABS(rr);
            T[e][0][1] = c * MFST_ABS(rr - delta);
            T[e][1][0] = c * MFST_ABS(rr + delta);
#undef MFST_ABS
        }
        int8_t x[N];
        int64_t got = mfst_solve_binary(arena, N, U, NE, euv, T, x);
        int64_t want = mfst_brute(N, U, NE, euv, T);
        mfst_check(got == want, "t4 move energy optimal", &fails);
        /* optimal: shift nothing (labels all 0) costs 5; shifting {0,1,2}
         * makes r3 -> -1+(1-0)... r_e + delta*(x_u - x_v): edge 2-3 with
         * x2=1,x3=0 -> r=-1+1=0 cost 0; edges 0-1,1-2 stay 0. total 0 */
        mfst_check(x[0] == 1 && x[1] == 1 && x[2] == 1 && x[3] == 0
                   && got == 0, "t4 chain flip found", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[maxflow selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
