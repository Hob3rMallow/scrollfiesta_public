#include "reeb.h"

#include "../common/csr.h"
#include "../common/bfs.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

void ReebOpts_default(ReebOpts *o)
{
    if (!o) return;
    o->neck_max_width = 8;
    o->rebound_factor = 2.5f;
    o->min_shoulder   = 12;
}

/* Single-source BFS; returns the farthest reachable vertex (graph eccentricity
 * endpoint). dist is caller scratch [nv], left holding the field on return. */
static int32_t bfs_farthest(const CSR_T g, int32_t seed, uint32_t *dist, int32_t n)
{
    int32_t seeds1[1];
    seeds1[0] = seed;
    memset(dist, 0xFF, (size_t)n * sizeof(uint32_t));   /* UINT32_MAX = unvisited */
    BFS_multi_source(g, seeds1, (size_t)1, dist);
    int32_t far = seed;
    uint32_t best = 0;
    for (int32_t v = 0; v < n; v++)
        if (dist[v] != UINT32_MAX && dist[v] > best) { best = dist[v]; far = v; }
    return far;
}

/* Is `l` (1..mmax-1) an interior valley deep enough to be a bridge neck? */
static int is_valley(size_t l, const int32_t *cut,
                     const int32_t *Lmax, const int32_t *Rmax,
                     const ReebOpts *o)
{
    int32_t w = cut[l];
    if (w > o->neck_max_width) return 0;
    int32_t lsh = Lmax[l - 1];
    int32_t rsh = Rmax[l + 1];
    if (lsh < o->min_shoulder || rsh < o->min_shoulder) return 0;
    float wv = (float)(w > 0 ? w : 1);
    if ((float)lsh < o->rebound_factor * wv) return 0;
    if ((float)rsh < o->rebound_factor * wv) return 0;
    return 1;
}

int Reeb_analyze(Arena_T arena,
                 const float *verts, size_t nv,
                 const int32_t *faces, size_t nf,
                 const ReebOpts *opts,
                 ReebResult *out)
{
    assert(arena && out);
    (void)verts;
    memset(out, 0, sizeof(*out));
    if (!faces || nv == 0 || nf == 0) return -1;

    ReebOpts o;
    ReebOpts_default(&o);
    if (opts) o = *opts;

    int32_t n = (int32_t)nv;
    CSR_T g = CSR_from_faces(arena, faces, nf, nv);
    const int32_t *off = CSR_offset(g);
    const int32_t *tgt = CSR_target(g);

    /* Always-present outputs. */
    uint32_t *morse     = (uint32_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(uint32_t));
    uint8_t  *face_flag = (uint8_t  *)ARENA_CALLOC(arena, (long)nf, 1L);
    out->morse      = morse;
    out->face_flag  = face_flag;
    out->bottleneck = -1;

    /* A referenced (positive-degree) start vertex. */
    int32_t v0 = -1;
    for (int32_t v = 0; v < n; v++)
        if (off[v + 1] > off[v]) { v0 = v; break; }
    if (v0 < 0) {                                   /* no edges -> nothing to do */
        out->cut_profile = (int32_t *)ARENA_CALLOC(arena, 1L, (long)sizeof(int32_t));
        out->morse_max = 0;
        return 0;
    }

    /* Two-sweep diameter endpoint, then the Morse field from it. */
    uint32_t *dist = (uint32_t *)ARENA_ALLOC(arena, (long)((size_t)nv * sizeof(uint32_t)));
    int32_t a = bfs_farthest(g, v0, dist, n);
    (void)bfs_farthest(g, a, dist, n);              /* dist now = field from a */

    uint32_t mmax = 0;
    for (size_t v = 0; v < nv; v++) {
        uint32_t d = (dist[v] == UINT32_MAX) ? 0u : dist[v];   /* isolated -> level 0 */
        morse[v] = d;
        if (d > mmax) mmax = d;
    }
    out->morse_max = mmax;

    size_t L = (size_t)mmax + 1;
    int32_t *cut = (int32_t *)ARENA_CALLOC(arena, (long)(L + 1), (long)sizeof(int32_t));

    /* Difference array over deduped undirected edges (count each once, u>v). */
    for (int32_t v = 0; v < n; v++) {
        uint32_t fv = morse[v];
        for (int32_t j = off[v]; j < off[v + 1]; j++) {
            int32_t u = tgt[j];
            if (u <= v) continue;
            uint32_t fu = morse[u];
            uint32_t lo = (fv < fu) ? fv : fu;
            uint32_t hi = (fv < fu) ? fu : fv;
            if (lo < hi) { cut[lo]++; cut[hi]--; }
        }
    }
    for (size_t l = 1; l < L; l++) cut[l] += cut[l - 1];   /* prefix -> C(l) */
    out->cut_profile = cut;

    if (mmax < 2) return 0;                          /* no interior level */

    /* Prefix / suffix maxima of the cut profile. */
    int32_t *Lmax = (int32_t *)ARENA_ALLOC(arena, (long)(L * sizeof(int32_t)));
    int32_t *Rmax = (int32_t *)ARENA_ALLOC(arena, (long)(L * sizeof(int32_t)));
    Lmax[0] = cut[0];
    for (size_t l = 1; l < L; l++) Lmax[l] = (cut[l] > Lmax[l - 1]) ? cut[l] : Lmax[l - 1];
    Rmax[L - 1] = cut[L - 1];
    for (size_t l = L - 1; l > 0; l--) {
        size_t li = l - 1;
        Rmax[li] = (cut[li] > Rmax[li + 1]) ? cut[li] : Rmax[li + 1];
    }

    /* Pass 1: count valley groups (contiguous candidate levels = one neck). */
    size_t n_necks = 0;
    {
        size_t l = 1;
        while (l + 1 <= (size_t)mmax) {              /* l in [1, mmax-1] */
            if (is_valley(l, cut, Lmax, Rmax, &o)) {
                n_necks++;
                while (l + 1 <= (size_t)mmax && is_valley(l, cut, Lmax, Rmax, &o)) l++;
            } else {
                l++;
            }
        }
    }
    if (n_necks == 0) return 0;

    ReebNeck *necks = (ReebNeck *)ARENA_CALLOC(arena, (long)n_necks,
                                               (long)sizeof(ReebNeck));
    /* Pass 2: emit one neck per group, taking the deepest level in the group. */
    {
        size_t k = 0, l = 1;
        while (l + 1 <= (size_t)mmax) {
            if (is_valley(l, cut, Lmax, Rmax, &o)) {
                size_t start = l, best = l;
                while (l + 1 <= (size_t)mmax && is_valley(l, cut, Lmax, Rmax, &o)) {
                    if (cut[l] < cut[best]) best = l;
                    l++;
                }
                size_t end = l - 1;
                int32_t lsh = Lmax[best - 1], rsh = Rmax[best + 1];
                necks[k].level    = (uint32_t)best;
                necks[k].lo       = (uint32_t)start;
                necks[k].hi       = (uint32_t)end;
                necks[k].width    = cut[best];
                necks[k].shoulder = (lsh < rsh) ? lsh : rsh;
                necks[k].n_faces  = 0;
                k++;
            } else {
                l++;
            }
        }
        n_necks = k;
    }
    out->necks   = necks;
    out->n_necks = n_necks;

    /* Localize: a face straddles a neck level l* if its Morse span contains l*.
     * Flag those faces and tally per-neck band size; track the bottleneck. */
    int32_t bottleneck = INT32_MAX;
    for (size_t k = 0; k < n_necks; k++)
        if (necks[k].width < bottleneck) bottleneck = necks[k].width;
    out->bottleneck = (n_necks > 0) ? bottleneck : -1;

    for (size_t fi = 0; fi < nf; fi++) {
        int32_t i0 = faces[fi * 3 + 0];
        int32_t i1 = faces[fi * 3 + 1];
        int32_t i2 = faces[fi * 3 + 2];
        if (i0 < 0 || i1 < 0 || i2 < 0 ||
            (size_t)i0 >= nv || (size_t)i1 >= nv || (size_t)i2 >= nv) continue;
        uint32_t f0 = morse[i0], f1 = morse[i1], f2 = morse[i2];
        uint32_t fmin = f0, fmax = f0;
        if (f1 < fmin) fmin = f1; if (f1 > fmax) fmax = f1;
        if (f2 < fmin) fmin = f2; if (f2 > fmax) fmax = f2;
        if (fmin == fmax) continue;                 /* face lies within one level */
        for (size_t k = 0; k < n_necks; k++) {
            /* flag faces whose Morse span overlaps the valley level band */
            if (fmin <= necks[k].hi && fmax >= necks[k].lo) {
                if (!face_flag[fi]) face_flag[fi] = 1;
                necks[k].n_faces++;
            }
        }
    }
    return 0;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

/* Triangulated grid nu x nh (wrap closes the u direction -> annulus). Verts are
 * laid out in a plane; geometry is irrelevant to the (graph-only) Morse sweep. */
static void rb_build_grid(Arena_T arena, int nu, int nh, int wrap,
                          float **out_v, int32_t **out_f, size_t *out_nv, size_t *out_nf)
{
    size_t nvv = (size_t)nu * (size_t)nh;
    int ucols = wrap ? nu : nu - 1;
    size_t nff = (size_t)ucols * (size_t)(nh - 1) * 2;
    float   *v = (float *)  ARENA_ALLOC(arena, (long)(nvv * 3 * sizeof(float)));
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nff * 3 * sizeof(int32_t)));
    for (int j = 0; j < nh; j++)
        for (int i = 0; i < nu; i++) {
            size_t idx = (size_t)j * (size_t)nu + (size_t)i;
            v[idx * 3 + 0] = (float)j;
            v[idx * 3 + 1] = (float)i;
            v[idx * 3 + 2] = 0.0f;
        }
    size_t fi = 0;
    for (int j = 0; j < nh - 1; j++)
        for (int i = 0; i < ucols; i++) {
            int i2 = wrap ? (i + 1) % nu : i + 1;
            int32_t aa = (int32_t)((size_t)j * (size_t)nu + (size_t)i);
            int32_t bb = (int32_t)((size_t)j * (size_t)nu + (size_t)i2);
            int32_t cc = (int32_t)((size_t)(j + 1) * (size_t)nu + (size_t)i);
            int32_t dd = (int32_t)((size_t)(j + 1) * (size_t)nu + (size_t)i2);
            f[fi * 3 + 0] = aa; f[fi * 3 + 1] = bb; f[fi * 3 + 2] = cc; fi++;
            f[fi * 3 + 0] = bb; f[fi * 3 + 1] = dd; f[fi * 3 + 2] = cc; fi++;
        }
    *out_v = v; *out_f = f; *out_nv = nvv; *out_nf = fi;
}

static int rb_chk(const char *name, int cond)
{
    fprintf(stderr, "  [reeb %s] %s\n", name, cond ? "ok" : "FAIL");
    return cond ? 0 : 1;
}

int Reeb_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    ReebResult r;

    /* Disk: open 20x12 grid -> monotone cut profile, no interior neck. */
    {
        float *v; int32_t *f; size_t lnv, lnf;
        rb_build_grid(arena, 20, 12, 0, &v, &f, &lnv, &lnf);
        int rc = Reeb_analyze(arena, v, lnv, f, lnf, NULL, &r);
        fprintf(stderr, "  disk: rc=%d mmax=%u necks=%zu bottleneck=%d\n",
                rc, r.morse_max, r.n_necks, r.bottleneck);
        fails += rb_chk("disk-no-neck", rc == 0 && r.n_necks == 0);
    }

    /* Annulus: 24x8 wrap grid -> uniform cross-section, no thin pinch. */
    {
        float *v; int32_t *f; size_t lnv, lnf;
        rb_build_grid(arena, 24, 8, 1, &v, &f, &lnv, &lnf);
        int rc = Reeb_analyze(arena, v, lnv, f, lnf, NULL, &r);
        fprintf(stderr, "  annulus: rc=%d mmax=%u necks=%zu bottleneck=%d\n",
                rc, r.morse_max, r.n_necks, r.bottleneck);
        fails += rb_chk("annulus-no-neck", rc == 0 && r.n_necks == 0);
    }

    /* Two 8x8 grids joined by a SINGLE triangle (vertex 63 of A to edge 64-65 of
     * B). The cut profile must pinch to ~2 at the bridge, with wide grids on both
     * sides -> exactly one neck, localized onto the bridge triangle. */
    {
        float *vA; int32_t *fA; size_t nvA, nfA;
        float *vB; int32_t *fB; size_t nvB, nfB;
        rb_build_grid(arena, 8, 8, 0, &vA, &fA, &nvA, &nfA);
        rb_build_grid(arena, 8, 8, 0, &vB, &fB, &nvB, &nfB);

        size_t nv = nvA + nvB;
        size_t nf = nfA + nfB + 1;                 /* + one bridge triangle */
        float   *v = (float *)  ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
        int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nf * 3 * sizeof(int32_t)));
        for (size_t i = 0; i < nvA * 3; i++) v[i] = vA[i];
        for (size_t i = 0; i < nvB; i++) {          /* offset B far from A in space */
            v[(nvA + i) * 3 + 0] = vB[i * 3 + 0] + 100.0f;
            v[(nvA + i) * 3 + 1] = vB[i * 3 + 1];
            v[(nvA + i) * 3 + 2] = vB[i * 3 + 2];
        }
        size_t k = 0;
        for (size_t t = 0; t < nfA * 3; t++) f[k++] = fA[t];
        for (size_t t = 0; t < nfB * 3; t++) f[k++] = (int32_t)(fB[t] + (int32_t)nvA);
        /* bridge triangle: A's corner 63, B's 0 and 1 (-> global nvA+0, nvA+1) */
        int32_t bridge = (int32_t)(nf - 1);
        f[k++] = 63;
        f[k++] = (int32_t)nvA + 0;
        f[k++] = (int32_t)nvA + 1;

        int rc = Reeb_analyze(arena, v, nv, f, nf, NULL, &r);
        int found = (rc == 0 && r.n_necks >= 1);
        int localized = (rc == 0 && r.face_flag && r.face_flag[bridge] == 1);
        int thin = (rc == 0 && r.bottleneck >= 0 && r.bottleneck <= 4);
        fprintf(stderr, "  bridge: rc=%d mmax=%u necks=%zu bottleneck=%d "
                "bridge_flagged=%d\n",
                rc, r.morse_max, r.n_necks, r.bottleneck,
                r.face_flag ? r.face_flag[bridge] : -1);
        fails += rb_chk("bridge-found", found);
        fails += rb_chk("bridge-localized", localized);
        fails += rb_chk("bridge-thin", thin);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "=== Reeb_selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails;
}
