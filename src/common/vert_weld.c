#include "vert_weld.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Spatial-hash welding with union-find. The hash cells are sized to the
 * weld eps, so the only neighbours a vert can possibly merge with are
 * in its 27-cell box.
 *
 * Union-find rule: parent[max(r_i, r_j)] = min(r_i, r_j). The smallest
 * index in each component is the root. Paths are compressed during find.
 *
 * The spatial-hash cell layout is built from the SORTED entry list, so
 * iteration order is stable across runs given the same input order.
 */

typedef struct {
    int32_t cell_z, cell_y, cell_x;
    int32_t vert_idx;
} CellEntry;

static inline uint64_t cell_key(int32_t cz, int32_t cy, int32_t cx)
{
    uint64_t uz = (uint64_t)(cz + (1 << 20)) & 0x1FFFFF;
    uint64_t uy = (uint64_t)(cy + (1 << 20)) & 0x1FFFFF;
    uint64_t ux = (uint64_t)(cx + (1 << 20)) & 0x1FFFFF;
    return (uz << 42) | (uy << 21) | ux;
}

static int cmp_cell_entry(const void *a, const void *b)
{
    const CellEntry *ea = (const CellEntry *)a;
    const CellEntry *eb = (const CellEntry *)b;
    uint64_t ka = cell_key(ea->cell_z, ea->cell_y, ea->cell_x);
    uint64_t kb = cell_key(eb->cell_z, eb->cell_y, eb->cell_x);
    if (ka < kb) return -1;
    if (ka > kb) return 1;
    if (ea->vert_idx < eb->vert_idx) return -1;
    if (ea->vert_idx > eb->vert_idx) return 1;
    return 0;
}

/* Stable LSD radix sort of the entries by cell key, 8 bits per pass, passes
 * whose byte is constant skipped (2026-10-06).  The entries arrive in vertex
 * order and cmp_cell_entry breaks key ties by vertex index, so this yields the
 * array qsort produced.  keys[] comes back aligned with the sorted entries. */
static void weld_sort_entries(Arena_T arena, CellEntry *entries, uint64_t *keys, size_t nv)
{
    CellEntry *te = (CellEntry *)ARENA_ALLOC(arena, nv * sizeof(CellEntry));
    uint64_t *tk = (uint64_t *)ARENA_ALLOC(arena, nv * sizeof(uint64_t));
    size_t *cnt = (size_t *)ARENA_ALLOC(arena, 256 * sizeof(size_t));
    for (int shift = 0; shift < 64; shift += 8) {
        memset(cnt, 0, 256 * sizeof(size_t));
        for (size_t i = 0; i < nv; i++) cnt[(keys[i] >> shift) & 0xFF]++;
        if (cnt[(keys[0] >> shift) & 0xFF] == nv) continue;   /* constant byte */
        size_t run = 0;
        for (int b = 0; b < 256; b++) { size_t c = cnt[b]; cnt[b] = run; run += c; }
        for (size_t i = 0; i < nv; i++) {
            size_t d = cnt[(keys[i] >> shift) & 0xFF]++;
            te[d] = entries[i];
            tk[d] = keys[i];
        }
        memcpy(entries, te, nv * sizeof(CellEntry));
        memcpy(keys, tk, nv * sizeof(uint64_t));
    }
}

/* splitmix64 finalizer (as mls_project.c) for the run table below. */
static inline uint64_t weld_hash64(uint64_t h)
{
    h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL;
    h ^= h >> 27; h *= 0x94D049BB133111EBULL;
    h ^= h >> 31;
    return h;
}

#define WELD_EMPTY UINT64_MAX

/* First index of key's run in the sorted keys, or SIZE_MAX when absent: what
 * the binary search found, without its ~log2(nv) steps per cell. */
static inline size_t weld_run_find(const uint64_t *tab_key, const size_t *tab_lo,
                                   size_t T, uint64_t key)
{
    size_t h = (size_t)(weld_hash64(key) & (uint64_t)(T - 1));
    for (;;) {
        uint64_t k = tab_key[h];
        if (k == key) return tab_lo[h];
        if (k == WELD_EMPTY) return SIZE_MAX;
        h = (h + 1) & (T - 1);
    }
}

/* Union-find with smaller-index-as-root + path compression. */
static int32_t uf_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]]; /* halving */
        x = parent[x];
    }
    return x;
}

static void uf_union(int32_t *parent, int32_t a, int32_t b)
{
    int32_t ra = uf_find(parent, a);
    int32_t rb = uf_find(parent, b);
    if (ra == rb) return;
    if (ra < rb) parent[rb] = ra;
    else parent[ra] = rb;
}

void Weld_verts_filtered(Arena_T arena,
                         const float *verts, size_t nv,
                         const float *in_normals,
                         int32_t *faces, size_t nf, size_t *out_nf,
                         float eps,
                         bool guard_orient,
                         const int32_t *merge_group,
                         WeldPairFilter pair_filter,
                         void *pair_filter_context,
                         float **out_verts, size_t *out_nv,
                         float **out_normals)
{
    assert(arena);
    assert(verts);
    assert(out_verts);
    assert(out_nv);
    assert(out_nf);
    assert(eps > 0.0f);

    if (nv == 0) {
        *out_verts = NULL;
        *out_nv = 0;
        *out_nf = 0;
        if (out_normals) *out_normals = NULL;
        return;
    }

    Arena_Mark mark = Arena_save(arena);

    double E = (double)eps;
    double E2 = E * E;
    double inv_E = 1.0 / E;

    /* -------- Build sorted cell entries -------- */
    CellEntry *entries = (CellEntry *)ARENA_ALLOC(arena,
                            (size_t)nv * sizeof(CellEntry));
    for (size_t i = 0; i < nv; i++) {
        float z = verts[i*3+0];
        float y = verts[i*3+1];
        float x = verts[i*3+2];
        entries[i].cell_z = (int32_t)floorf(z * (float)inv_E);
        entries[i].cell_y = (int32_t)floorf(y * (float)inv_E);
        entries[i].cell_x = (int32_t)floorf(x * (float)inv_E);
        entries[i].vert_idx = (int32_t)i;
    }
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena,
                        (size_t)nv * sizeof(uint64_t));
    for (size_t i = 0; i < nv; i++) {
        keys[i] = cell_key(entries[i].cell_z, entries[i].cell_y, entries[i].cell_x);
    }
    weld_sort_entries(arena, entries, keys, nv);

    /* Run table: each distinct key -> the first index of its run. */
    size_t T = 64;
    while (T < nv * 2) T <<= 1;
    uint64_t *tab_key = (uint64_t *)ARENA_ALLOC(arena, T * sizeof(uint64_t));
    size_t *tab_lo = (size_t *)ARENA_ALLOC(arena, T * sizeof(size_t));
    memset(tab_key, 0xFF, T * sizeof(uint64_t));   /* all WELD_EMPTY */
    for (size_t i = 0; i < nv; i++) {
        if (i > 0 && keys[i] == keys[i - 1]) continue;
        size_t h = (size_t)(weld_hash64(keys[i]) & (uint64_t)(T - 1));
        while (tab_key[h] != WELD_EMPTY) h = (h + 1) & (T - 1);
        tab_key[h] = keys[i];
        tab_lo[h] = i;
    }

    /* -------- Union-find merging -------- */
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena,
                         (size_t)nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) parent[i] = (int32_t)i;

    /* Orientation guard (guard_orient && nf>0): a coincident pair is fused
     * only when its area-weighted winding normals agree. Merging two
     * opposite-facing (recto vs verso) verts fuses the sheet to its own back
     * = a non-orientable knot (a `same_dir` edge no winding-repair BFS can
     * fix). The mesh is consistently wound here (the caller runs OrientMesh
     * first), so an opposing normal dot is exactly a would-be same_dir fusion.
     * Per-vertex normal = normalized sum of un-normalized incident-face cross
     * products (flip-sensitive by construction); see orient_weld.c. */
    int do_guard = (guard_orient && nf > 0) ? 1 : 0;
    double *vn = NULL;
    if (do_guard) {
        vn = (double *)ARENA_CALLOC(arena, (size_t)nv * 3L, sizeof(double));
        for (size_t f = 0; f < nf; f++) {
            int32_t i0 = faces[f*3+0], i1 = faces[f*3+1], i2 = faces[f*3+2];
            const float *p0 = &verts[(size_t)i0*3];
            const float *p1 = &verts[(size_t)i1*3];
            const float *p2 = &verts[(size_t)i2*3];
            double a0 = p1[0]-p0[0], a1 = p1[1]-p0[1], a2 = p1[2]-p0[2];
            double b0 = p2[0]-p0[0], b1 = p2[1]-p0[1], b2 = p2[2]-p0[2];
            double c0 = a1*b2 - a2*b1, c1 = a2*b0 - a0*b2, c2 = a0*b1 - a1*b0;
            int32_t idx[3] = { i0, i1, i2 };
            for (int t = 0; t < 3; t++) {
                vn[(size_t)idx[t]*3+0] += c0;
                vn[(size_t)idx[t]*3+1] += c1;
                vn[(size_t)idx[t]*3+2] += c2;
            }
        }
    }

    /* For each vert (in original index order), look at all neighbours
     * within eps in the 27-cell box and union with them. Visiting in
     * original index order makes the union pattern stable. */
    for (size_t i = 0; i < nv; i++) {
        float vz = verts[i*3+0];
        float vy = verts[i*3+1];
        float vx = verts[i*3+2];
        int32_t ciz = (int32_t)floorf(vz * (float)inv_E);
        int32_t ciy = (int32_t)floorf(vy * (float)inv_E);
        int32_t cix = (int32_t)floorf(vx * (float)inv_E);

        for (int dz = -1; dz <= 1; dz++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    uint64_t target = cell_key(ciz + dz, ciy + dy, cix + dx);
                    size_t lo = weld_run_find(tab_key, tab_lo, T, target);
                    if (lo == SIZE_MAX) continue;
                    for (size_t k = lo; k < nv && keys[k] == target; k++) {
                        int32_t j = entries[k].vert_idx;
                        /* Each unordered pair is considered exactly once. */
                        if ((size_t)j <= i) continue;
                        double dvz = (double)verts[j*3+0] - (double)vz;
                        double dvy = (double)verts[j*3+1] - (double)vy;
                        double dvx = (double)verts[j*3+2] - (double)vx;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 <= E2) {
                            /* A late cleanup may encounter bit-coincident points
                             * on two independent sheets.  Normal agreement is
                             * insufficient there: two stacked sheets can be
                             * same-facing.  Preserve the caller's connectivity
                             * partition before considering geometric guards. */
                            if (merge_group != NULL &&
                                merge_group[i] != merge_group[(size_t)j])
                                continue;
                            if (pair_filter != NULL &&
                                !pair_filter(i, (size_t)j,
                                             pair_filter_context))
                                continue;
                            if (do_guard) {
                                const double *ni = &vn[i*3];
                                const double *nj = &vn[(size_t)j*3];
                                /* Opposing winding normals => recto vs verso.
                                 * Fusing them makes a non-orientable knot; skip
                                 * this pair (leave the slit for the hole-fill). */
                                if (ni[0]*nj[0] + ni[1]*nj[1] + ni[2]*nj[2] < 0.0)
                                    continue;
                            }
                            uf_union(parent, (int32_t)i, j);
                        }
                    }
                }
            }
        }
    }

    /* -------- Compress all paths and count unique roots -------- */
    for (size_t i = 0; i < nv; i++) {
        parent[i] = uf_find(parent, (int32_t)i);
    }

    /* Map root -> new index. Walk in ascending index, assigning a new
     * index the first time each root is seen. Roots are by construction
     * the smallest index in their component, so this is identical to
     * "assign new indices in root-index order." */
    int32_t *new_idx = (int32_t *)ARENA_ALLOC(arena,
                         (size_t)nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) new_idx[i] = -1;

    size_t n_new = 0;
    for (size_t i = 0; i < nv; i++) {
        int32_t r = parent[i];
        if (new_idx[r] < 0) {
            new_idx[r] = (int32_t)n_new;
            n_new++;
        }
    }
    /* Propagate new index to every vert (so new_idx[i] is valid for all i). */
    for (size_t i = 0; i < nv; i++) {
        new_idx[i] = new_idx[parent[i]];
    }

    /* -------- Accumulate centroid (and optional normals) per group -------- */
    /* Final output arrays live on the caller's arena; allocate now. */
    float *welded = (float *)ARENA_ALLOC(arena, (size_t)(n_new * 3 * sizeof(float)));
    float *welded_n = NULL;
    if (out_normals && in_normals) {
        welded_n = (float *)ARENA_ALLOC(arena, (size_t)(n_new * 3 * sizeof(float)));
    }
    double *acc = (double *)ARENA_CALLOC(arena, (size_t)n_new * 3L, sizeof(double));
    double *acc_n = NULL;
    if (welded_n) {
        acc_n = (double *)ARENA_CALLOC(arena, (size_t)n_new * 3L, sizeof(double));
    }
    int32_t *count = (int32_t *)ARENA_CALLOC(arena, (size_t)n_new, sizeof(int32_t));

    for (size_t i = 0; i < nv; i++) {
        int32_t k = new_idx[i];
        acc[k*3+0] += (double)verts[i*3+0];
        acc[k*3+1] += (double)verts[i*3+1];
        acc[k*3+2] += (double)verts[i*3+2];
        count[k]++;
        if (acc_n) {
            acc_n[k*3+0] += (double)in_normals[i*3+0];
            acc_n[k*3+1] += (double)in_normals[i*3+1];
            acc_n[k*3+2] += (double)in_normals[i*3+2];
        }
    }
    for (size_t k = 0; k < n_new; k++) {
        double c = count[k] > 0 ? (1.0 / (double)count[k]) : 0.0;
        welded[k*3+0] = (float)(acc[k*3+0] * c);
        welded[k*3+1] = (float)(acc[k*3+1] * c);
        welded[k*3+2] = (float)(acc[k*3+2] * c);
        if (welded_n) {
            double nz = acc_n[k*3+0] * c;
            double ny = acc_n[k*3+1] * c;
            double nx = acc_n[k*3+2] * c;
            double nlen = sqrt(nz*nz + ny*ny + nx*nx);
            if (nlen > 1e-15) {
                welded_n[k*3+0] = (float)(nz / nlen);
                welded_n[k*3+1] = (float)(ny / nlen);
                welded_n[k*3+2] = (float)(nx / nlen);
            } else {
                welded_n[k*3+0] = 0.0f;
                welded_n[k*3+1] = 0.0f;
                welded_n[k*3+2] = 1.0f;
            }
        }
    }

    /* -------- Remap faces, drop degenerates -------- */
    size_t kept = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = new_idx[faces[f*3+0]];
        int32_t b = new_idx[faces[f*3+1]];
        int32_t c = new_idx[faces[f*3+2]];
        if (a == b || b == c || a == c) continue;
        faces[kept*3+0] = a;
        faces[kept*3+1] = b;
        faces[kept*3+2] = c;
        kept++;
    }

    *out_verts = welded;
    *out_nv = n_new;
    *out_nf = kept;
    if (out_normals) *out_normals = welded_n;

    /* NOTE: we do NOT Arena_restore here — the output arrays are on this
     * arena and must outlive this call. The caller's arena lifetime
     * subsumes them. The scratch (entries, keys, parent, new_idx, acc,
     * count, acc_n) is dead after this point, so a future Arena_save in
     * the caller would only reclaim the bytes if it brackets this call;
     * that's the existing convention in step0_mesh_extract.c. */
    (void)mark;
}

void Weld_verts_grouped(Arena_T arena,
                        const float *verts, size_t nv,
                        const float *in_normals,
                        int32_t *faces, size_t nf, size_t *out_nf,
                        float eps,
                        bool guard_orient,
                        const int32_t *merge_group,
                        float **out_verts, size_t *out_nv,
                        float **out_normals)
{
    Weld_verts_filtered(arena, verts, nv, in_normals, faces, nf, out_nf,
                        eps, guard_orient, merge_group, NULL, NULL,
                        out_verts, out_nv, out_normals);
}

void Weld_verts(Arena_T arena,
                const float *verts, size_t nv,
                const float *in_normals,
                int32_t *faces, size_t nf, size_t *out_nf,
                float eps,
                bool guard_orient,
                float **out_verts, size_t *out_nv,
                float **out_normals)
{
    Weld_verts_filtered(arena, verts, nv, in_normals, faces, nf, out_nf,
                        eps, guard_orient, NULL, NULL, NULL,
                        out_verts, out_nv, out_normals);
}

/* Selftest (cube_mesh --selftest): the radix sort equals qsort with
 * cmp_cell_entry, and the run table finds exactly the binary search's first
 * index, on random cells with heavy duplication and negative coordinates. */
int Weld_selftest(void)
{
    Arena_T arena = Arena_new();
    uint64_t s = 0x2545F4914F6CDD1DULL;
    int fails = 0;
    size_t checks = 0;
    for (int trial = 0; trial < 40 && !fails; trial++) {
        size_t nv = (size_t)(1 + trial * trial * 37);
        int span = 1 + trial % 9;
        Arena_Mark m = Arena_save(arena);
        CellEntry *a = (CellEntry *)ARENA_ALLOC(arena, nv * sizeof(CellEntry));
        CellEntry *b = (CellEntry *)ARENA_ALLOC(arena, nv * sizeof(CellEntry));
        uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena, nv * sizeof(uint64_t));
        for (size_t i = 0; i < nv; i++) {
            s = weld_hash64(s + 0x9E3779B97F4A7C15ULL);
            a[i].cell_z = (int32_t)(s % (uint64_t)span) - span / 2 - (trial & 1) * 70000;
            a[i].cell_y = (int32_t)((s >> 16) % (uint64_t)span) - span / 2;
            a[i].cell_x = (int32_t)((s >> 32) % (uint64_t)span) + (trial & 2) * 1000;
            a[i].vert_idx = (int32_t)i;
            b[i] = a[i];
            keys[i] = cell_key(a[i].cell_z, a[i].cell_y, a[i].cell_x);
        }
        qsort(b, nv, sizeof(CellEntry), cmp_cell_entry);
        weld_sort_entries(arena, a, keys, nv);
        checks++;
        if (memcmp(a, b, nv * sizeof(CellEntry)) != 0) { fails++; break; }
        for (size_t i = 0; i < nv; i++)
            if (keys[i] != cell_key(a[i].cell_z, a[i].cell_y, a[i].cell_x)) { fails++; break; }
        size_t T = 64;
        while (T < nv * 2) T <<= 1;
        uint64_t *tab_key = (uint64_t *)ARENA_ALLOC(arena, T * sizeof(uint64_t));
        size_t *tab_lo = (size_t *)ARENA_ALLOC(arena, T * sizeof(size_t));
        memset(tab_key, 0xFF, T * sizeof(uint64_t));
        for (size_t i = 0; i < nv; i++) {
            if (i > 0 && keys[i] == keys[i - 1]) continue;
            size_t h = (size_t)(weld_hash64(keys[i]) & (uint64_t)(T - 1));
            while (tab_key[h] != WELD_EMPTY) h = (h + 1) & (T - 1);
            tab_key[h] = keys[i];
            tab_lo[h] = i;
        }
        for (int q = 0; q < 400; q++) {
            s = weld_hash64(s + 1);
            int32_t cz = (int32_t)(s % (uint64_t)(span + 4)) - span / 2 - 2 - (trial & 1) * 70000;
            int32_t cy = (int32_t)((s >> 16) % (uint64_t)(span + 4)) - span / 2 - 2;
            int32_t cx = (int32_t)((s >> 32) % (uint64_t)(span + 4)) - 2 + (trial & 2) * 1000;
            uint64_t target = cell_key(cz, cy, cx);
            size_t lo = 0, hi = nv;
            while (lo < hi) {
                size_t mid = lo + (hi - lo) / 2;
                if (keys[mid] < target) lo = mid + 1;
                else hi = mid;
            }
            size_t want = (lo < nv && keys[lo] == target) ? lo : SIZE_MAX;
            checks++;
            if (weld_run_find(tab_key, tab_lo, T, target) != want) { fails++; break; }
        }
        Arena_restore(arena, m);
    }
    Arena_dispose(&arena);
    fprintf(stderr, "Weld_selftest: %s (%zu radix-vs-qsort sorts and run lookups)\n",
            fails ? "FAIL" : "OK", checks);
    return fails ? -1 : 0;
}
