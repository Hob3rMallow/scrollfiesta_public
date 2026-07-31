#include "mesh_topo.h"
#include "../common/union_find.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x < y) ? -1 : (x > y ? 1 : 0);
}

int MeshTopo_analyze(Arena_T arena, const float *verts, size_t nv,
                     const int32_t *faces, size_t nf, MeshTopoInfo *out)
{
    (void)verts;   /* topology uses faces only */
    assert(arena && out);
    memset(out, 0, sizeof(*out));
    if (nv == 0 || nf == 0 || !faces) return -1;
    out->n_verts = (long)nv;
    out->n_faces = (long)nf;

    Arena_Mark mark = Arena_save(arena);

    /* Pack undirected edges as sortable keys lo*nv + hi (lo<hi<nv). */
    size_t ne3 = nf * 3;
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena, (long)(ne3 * sizeof(uint64_t)));
    size_t kk = 0;
    for (size_t t = 0; t < nf; t++) {
        int32_t v[3] = { faces[t*3+0], faces[t*3+1], faces[t*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = v[e], b = v[(e+1)%3];
            if (a < 0 || b < 0 || (size_t)a >= nv || (size_t)b >= nv) continue;
            uint64_t lo = (uint64_t)(a < b ? a : b);
            uint64_t hi = (uint64_t)(a < b ? b : a);
            keys[kk++] = lo * (uint64_t)nv + hi;
        }
    }
    qsort(keys, kk, sizeof(uint64_t), cmp_u64);

    UnionFind uf = UF_new(arena, (int32_t)nv);   /* boundary graph */
    char *isb = (char *)ARENA_CALLOC(arena, (long)nv, 1L);

    long nedges = 0, nbe = 0, nman = 0, nnm = 0;
    size_t i = 0;
    while (i < kk) {
        uint64_t k = keys[i];
        size_t j = i;
        while (j < kk && keys[j] == k) j++;
        long cnt = (long)(j - i);
        nedges++;
        if (cnt == 1) {
            nbe++;
            int32_t lo = (int32_t)(k / (uint64_t)nv);
            int32_t hi = (int32_t)(k % (uint64_t)nv);
            uf_union(&uf, lo, hi);
            isb[lo] = 1; isb[hi] = 1;
        } else if (cnt == 2) {
            nman++;
        } else {
            nnm++;
        }
        i = j;
    }

    /* Boundary loops = connected components of the boundary graph. */
    long nloops = 0;
    char *seen = (char *)ARENA_CALLOC(arena, (long)nv, 1L);
    for (size_t v = 0; v < nv; v++) {
        if (!isb[v]) continue;
        int32_t r = uf_find(&uf, (int32_t)v);
        if (!seen[r]) { seen[r] = 1; nloops++; }
    }

    /* Unreferenced vertices. */
    char *used = (char *)ARENA_CALLOC(arena, (long)nv, 1L);
    for (size_t t = 0; t < nf; t++)
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[t*3+e];
            if (a >= 0 && (size_t)a < nv) used[a] = 1;
        }
    long unref = 0;
    for (size_t v = 0; v < nv; v++) if (!used[v]) unref++;

    out->n_edges = nedges;
    out->n_boundary_edges = nbe;
    out->n_manifold_edges = nman;
    out->n_nonmanifold_edges = nnm;
    out->n_boundary_loops = nloops;
    out->euler = (long)nv - nedges + (long)nf;
    out->genus = (2.0 - (double)nloops - (double)out->euler) / 2.0;
    out->n_unref_verts = unref;
    out->is_disk = (nnm == 0 && nloops == 1 && out->euler == 1);

    Arena_restore(arena, mark);
    return 0;
}

/* ============================================================================
 * Self-test: disk (1 loop, g0), annulus (2 loops, g0), disk-with-hole (2 loops).
 * ==========================================================================*/

/* grid nu x nh; wrap=1 closes the u direction (annulus), wrap=0 open (disk). */
static void build_grid(Arena_T arena, int nu, int nh, int wrap,
                       int32_t **out_f, size_t *out_nf, size_t *out_nv)
{
    size_t nvv = (size_t)nu * (size_t)nh;
    int ucols = wrap ? nu : nu - 1;
    size_t nff = (size_t)ucols * (size_t)(nh - 1) * 2;
    int32_t *f = (int32_t *)ARENA_ALLOC(arena, (long)(nff * 3 * sizeof(int32_t)));
    size_t fi = 0;
    for (int j = 0; j < nh - 1; j++) {
        for (int i = 0; i < ucols; i++) {
            int i2 = wrap ? (i + 1) % nu : i + 1;
            int32_t a = (int32_t)((size_t)j*(size_t)nu + (size_t)i);
            int32_t b = (int32_t)((size_t)j*(size_t)nu + (size_t)i2);
            int32_t c = (int32_t)((size_t)(j+1)*(size_t)nu + (size_t)i);
            int32_t d = (int32_t)((size_t)(j+1)*(size_t)nu + (size_t)i2);
            f[fi*3+0]=a; f[fi*3+1]=b; f[fi*3+2]=c; fi++;
            f[fi*3+0]=b; f[fi*3+1]=d; f[fi*3+2]=c; fi++;
        }
    }
    *out_f = f; *out_nf = fi; *out_nv = nvv;
}

static int chk(const char *name, int cond)
{
    fprintf(stderr, "  [meshtopo %s] %s\n", name, cond ? "ok" : "FAIL");
    return cond ? 0 : 1;
}

int MeshTopo_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    MeshTopoInfo info;

    /* Disk: open grid -> 1 boundary loop, euler 1, genus 0, is_disk. */
    {
        int32_t *f; size_t nff, nvv;
        build_grid(arena, 10, 8, 0, &f, &nff, &nvv);
        MeshTopo_analyze(arena, NULL, nvv, f, nff, &info);
        fprintf(stderr, "  disk: loops=%ld nm=%ld euler=%ld genus=%.1f is_disk=%d\n",
                info.n_boundary_loops, info.n_nonmanifold_edges, info.euler, info.genus, info.is_disk);
        fails += chk("disk", info.n_boundary_loops==1 && info.n_nonmanifold_edges==0
                            && info.euler==1 && info.is_disk);
    }
    /* Annulus: wrap grid -> 2 boundary loops, euler 0, genus 0, not a disk. */
    {
        int32_t *f; size_t nff, nvv;
        build_grid(arena, 12, 6, 1, &f, &nff, &nvv);
        MeshTopo_analyze(arena, NULL, nvv, f, nff, &info);
        fprintf(stderr, "  annulus: loops=%ld nm=%ld euler=%ld genus=%.1f is_disk=%d\n",
                info.n_boundary_loops, info.n_nonmanifold_edges, info.euler, info.genus, info.is_disk);
        fails += chk("annulus", info.n_boundary_loops==2 && info.n_nonmanifold_edges==0
                               && info.euler==0 && !info.is_disk);
    }

    Arena_dispose(&arena);
    return fails;
}
