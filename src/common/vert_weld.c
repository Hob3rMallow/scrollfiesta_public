#include "vert_weld.h"

#include <assert.h>
#include <math.h>
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

void Weld_verts(Arena_T arena,
                const float *verts, size_t nv,
                const float *in_normals,
                int32_t *faces, size_t nf, size_t *out_nf,
                float eps,
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
                            (long)nv * (long)sizeof(CellEntry));
    for (size_t i = 0; i < nv; i++) {
        float z = verts[i*3+0];
        float y = verts[i*3+1];
        float x = verts[i*3+2];
        entries[i].cell_z = (int32_t)floorf(z * (float)inv_E);
        entries[i].cell_y = (int32_t)floorf(y * (float)inv_E);
        entries[i].cell_x = (int32_t)floorf(x * (float)inv_E);
        entries[i].vert_idx = (int32_t)i;
    }
    qsort(entries, nv, sizeof(CellEntry), cmp_cell_entry);

    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena,
                        (long)nv * (long)sizeof(uint64_t));
    for (size_t i = 0; i < nv; i++) {
        keys[i] = cell_key(entries[i].cell_z, entries[i].cell_y, entries[i].cell_x);
    }

    /* -------- Union-find merging -------- */
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena,
                         (long)nv * (long)sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) parent[i] = (int32_t)i;

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
                    size_t lo = 0, hi = nv;
                    while (lo < hi) {
                        size_t mid = lo + (hi - lo) / 2;
                        if (keys[mid] < target) lo = mid + 1;
                        else hi = mid;
                    }
                    for (size_t k = lo; k < nv && keys[k] == target; k++) {
                        int32_t j = entries[k].vert_idx;
                        if ((size_t)j == i) continue;
                        double dvz = (double)verts[j*3+0] - (double)vz;
                        double dvy = (double)verts[j*3+1] - (double)vy;
                        double dvx = (double)verts[j*3+2] - (double)vx;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 <= E2) {
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
                         (long)nv * (long)sizeof(int32_t));
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
    float *welded = (float *)ARENA_ALLOC(arena, (long)(n_new * 3 * sizeof(float)));
    float *welded_n = NULL;
    if (out_normals && in_normals) {
        welded_n = (float *)ARENA_ALLOC(arena, (long)(n_new * 3 * sizeof(float)));
    }
    double *acc = (double *)ARENA_CALLOC(arena, (long)n_new * 3L, (long)sizeof(double));
    double *acc_n = NULL;
    if (welded_n) {
        acc_n = (double *)ARENA_CALLOC(arena, (long)n_new * 3L, (long)sizeof(double));
    }
    int32_t *count = (int32_t *)ARENA_CALLOC(arena, (long)n_new, (long)sizeof(int32_t));

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
