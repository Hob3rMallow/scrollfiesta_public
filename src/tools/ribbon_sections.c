/* ============================================================================
 * ribbon_sections -- carve a labeled ribbon into per-label section dirs.
 *
 *   ribbon_sections <ribbon.vmesh> <labels.i32> <out_dir>
 *                   [--top N] [--min-frac F] [--census]
 *                   [--phase <phase.f32>] [--recon <recon.i32>]
 *   ribbon_sections --selftest
 *
 * The fit's big material sections are internally grid-coherent: each one IS
 * a ribbon already.  This tool reports the label census (the kibble report)
 * and extracts the top-N labels (default 1) -- or every label holding >= F
 * of the vertices -- each into its own subdirectory:
 *
 *   <out_dir>/secNN/section.vmesh                    geometry + UV, compacted
 *   <out_dir>/secNN/ribbon_phase.f32                 sliced phase sidecar
 *   <out_dir>/secNN/ribbon_material_identity.i32     sliced labels
 *   <out_dir>/secNN/ribbon_reconstruction_component.i32   (if --recon given)
 *
 * The sidecar names are the metric re-fit's FIXED discovery names beside its
 * input: carving them is what arms island placement downstream.  Dropping
 * them was the measured logic bug -- an unplaced atlas packing fed the
 * consistency solve 4.2M cross-packing ties and the factor exploded.
 * Everything not kept is the GUTTER: counted, never emitted.  Faces are kept
 * only when all three vertices share the label, so sections never bridge
 * materials.
 * ==========================================================================*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/ves_platform.h"

typedef struct {
    int32_t label;
    size_t nv;
} RsCount;

static int rs_cmp_count_desc(const void *pa, const void *pb)
{
    const RsCount *a = (const RsCount *)pa;
    const RsCount *b = (const RsCount *)pb;
    if (a->nv != b->nv) return a->nv > b->nv ? -1 : 1;
    return a->label < b->label ? -1 : (a->label > b->label ? 1 : 0);
}

static void *rs_load_raw(Arena_T arena, const char *path, size_t elem,
                         size_t count, const char *what)
{
    FILE *f = fopen(path, "rb");
    void *data = NULL;
    size_t got = 0;
    if (f == NULL) {
        fprintf(stderr, "ribbon_sections: cannot open %s (%s)\n", path, what);
        return NULL;
    }
    data = ARENA_ALLOC(arena, (size_t)(count * elem));
    got = fread(data, elem, count, f);
    fclose(f);
    if (got != count) {
        fprintf(stderr,
                "ribbon_sections: %s holds %zu %s for %zu vertices\n",
                path, got, what, count);
        return NULL;
    }
    return data;
}

static int rs_write_raw(const char *path, const void *data, size_t elem,
                        size_t count)
{
    FILE *f = fopen(path, "wb");
    size_t put = 0;
    if (f == NULL) return -1;
    put = fwrite(data, elem, count, f);
    if (fclose(f) != 0) return -1;
    return put == count ? 0 : -1;
}

/* Slice a per-vertex 4-byte sidecar through the section remap. */
static int rs_slice_sidecar(Arena_T arena, const char *dir, const char *name,
                            const void *src, const int32_t *remap,
                            size_t nv, size_t snv)
{
    Arena_Mark mark = Arena_save(arena);
    uint32_t *out = (uint32_t *)ARENA_ALLOC(
        arena, (size_t)(snv * sizeof(uint32_t)));
    const uint32_t *in = (const uint32_t *)src;
    char path[2048];
    int rc = 0;
    for (size_t i = 0; i < nv; i++)
        if (remap[i] >= 0) out[remap[i]] = in[i];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    rc = rs_write_raw(path, out, sizeof(uint32_t), snv);
    if (rc != 0)
        fprintf(stderr, "ribbon_sections: cannot write %s\n", path);
    Arena_restore(arena, mark);
    return rc;
}

/* Coarsen one label onto a strided lattice: keep vertices sitting on every
 * (KU-th column, KV-th row) of the fit grid, rebuild quad connectivity on
 * the coarse lattice, and carry original UV coordinates so downstream u/v
 * semantics are unchanged.  A ~2-vox fit lattice at --coarsen 5 4 becomes a
 * 10x8-vox quadribbon with 5% of the vertices -- the factorization and
 * untangle costs drop ~20x while the bake still samples every raster texel
 * barycentrically.  Small holes (< one coarse cell) close, matching the
 * solid-ribbon contract; larger holes survive because their coarse corner
 * nodes are absent. */
static int rs_extract_coarse(Arena_T arena, const MeshBinData *m,
                             const int32_t *labels, const float *phase,
                             const int32_t *recon, int32_t label, int index,
                             const char *out_dir, long ku, long kv,
                             double cell_du, double cell_dv)
{
    Arena_Mark mark = Arena_save(arena);
    size_t cap = 1, nkeep = 0, snf = 0;
    int64_t *hkey = NULL;
    int32_t *hsrc = NULL, *hidx = NULL;
    char dir[1800], path[2048];
    int rc = -1;
    for (size_t i = 0; i < m->nv; i++)
        if (labels[i] == label) nkeep++;
    while (cap < 2 * nkeep + 16) cap <<= 1;
    hkey = (int64_t *)ARENA_ALLOC(arena, (size_t)(cap * sizeof(int64_t)));
    hsrc = (int32_t *)ARENA_ALLOC(arena, (size_t)(cap * sizeof(int32_t)));
    hidx = (int32_t *)ARENA_ALLOC(arena, (size_t)(cap * sizeof(int32_t)));
    for (size_t h = 0; h < cap; h++) hkey[h] = -1;
    nkeep = 0;
    for (size_t i = 0; i < m->nv; i++) {
        long col = 0, row = 0;
        int64_t key = 0;
        size_t h = 0;
        if (labels[i] != label) continue;
        col = lround((double)m->uv[i * 2 + 0] / cell_du);
        row = lround((double)m->uv[i * 2 + 1] / cell_dv);
        if (col % ku != 0 || row % kv != 0) continue;
        key = (((int64_t)(row / kv) + (1 << 20)) << 32) |
              ((int64_t)(col / ku) + (1 << 30));
        h = (size_t)(((uint64_t)key * 2654435761u) & (cap - 1));
        for (;; h = (h + 1) & (cap - 1)) {
            if (hkey[h] < 0) {
                hkey[h] = key;
                hsrc[h] = (int32_t)i;
                hidx[h] = (int32_t)nkeep++;
                break;
            }
            if (hkey[h] == key) break;   /* first vertex wins the cell */
        }
    }
    if (nkeep < 4) {
        fprintf(stderr,
                "  section %02d label %d: %zu coarse node(s); skipped\n",
                index, label, nkeep);
        Arena_restore(arena, mark);
        return -1;
    }
    {
        float *sverts = (float *)ARENA_ALLOC(
            arena, (size_t)(nkeep * 3 * sizeof(float)));
        float *suv = (float *)ARENA_ALLOC(
            arena, (size_t)(nkeep * 2 * sizeof(float)));
        int32_t *ssrc = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(nkeep * sizeof(int32_t)));
        int32_t *sfaces = NULL;
        size_t fcap = 0;
        double u_lo = 1e300, u_hi = -1e300;
        for (size_t h = 0; h < cap; h++) {
            if (hkey[h] < 0) continue;
            size_t r = (size_t)hidx[h], s = (size_t)hsrc[h];
            for (int k = 0; k < 3; k++)
                sverts[r * 3 + k] = m->verts[s * 3 + k];
            suv[r * 2 + 0] = m->uv[s * 2 + 0];
            suv[r * 2 + 1] = m->uv[s * 2 + 1];
            ssrc[r] = (int32_t)s;
            if (suv[r * 2] < u_lo) u_lo = suv[r * 2];
            if (suv[r * 2] > u_hi) u_hi = suv[r * 2];
        }
        /* coarse quads: (r,c)-(r,c+1)-(r+1,c+1)-(r+1,c) where all exist */
        fcap = nkeep * 2 + 16;
        sfaces = (int32_t *)ARENA_ALLOC(
            arena, (size_t)(fcap * 3 * sizeof(int32_t)));
        for (size_t h = 0; h < cap; h++) {
            int64_t key = hkey[h], k01, k10, k11;
            int32_t i00 = -1, i01 = -1, i10 = -1, i11 = -1;
            if (key < 0) continue;
            i00 = hidx[h];
            k01 = key + 1;
            k10 = key + ((int64_t)1 << 32);
            k11 = k10 + 1;
            {
                int64_t want[3] = { k01, k10, k11 };
                int32_t *got[3] = { &i01, &i10, &i11 };
                for (int q = 0; q < 3; q++) {
                    size_t hh = (size_t)(((uint64_t)want[q] * 2654435761u) &
                                         (cap - 1));
                    for (;; hh = (hh + 1) & (cap - 1)) {
                        if (hkey[hh] < 0) break;
                        if (hkey[hh] == want[q]) {
                            *got[q] = hidx[hh];
                            break;
                        }
                    }
                }
            }
            if (i01 < 0 || i10 < 0 || i11 < 0) continue;
            if (snf + 2 > fcap) continue;
            sfaces[snf * 3 + 0] = i00;
            sfaces[snf * 3 + 1] = i01;
            sfaces[snf * 3 + 2] = i11;
            snf++;
            sfaces[snf * 3 + 0] = i00;
            sfaces[snf * 3 + 1] = i11;
            sfaces[snf * 3 + 2] = i10;
            snf++;
        }
        snprintf(dir, sizeof dir, "%s/sec%02d", out_dir, index);
        ves_mkdir(dir);
        snprintf(path, sizeof path, "%s/section.vmesh", dir);
        rc = MeshBin_write(path, sverts, nkeep, sfaces, snf, suv);
        if (rc == 0) {
            /* sidecars: gather through ssrc (coarse node -> source vertex) */
            Arena_Mark smark = Arena_save(arena);
            uint32_t *buf = (uint32_t *)ARENA_ALLOC(
                arena, (size_t)(nkeep * sizeof(uint32_t)));
            if (phase != NULL && rc == 0) {
                const uint32_t *in = (const uint32_t *)(const void *)phase;
                for (size_t r = 0; r < nkeep; r++)
                    buf[r] = in[ssrc[r]];
                snprintf(path, sizeof path, "%s/ribbon_phase.f32", dir);
                rc = rs_write_raw(path, buf, sizeof(uint32_t), nkeep);
            }
            if (rc == 0) {
                const uint32_t *in = (const uint32_t *)(const void *)labels;
                for (size_t r = 0; r < nkeep; r++)
                    buf[r] = in[ssrc[r]];
                snprintf(path, sizeof path,
                         "%s/ribbon_material_identity.i32", dir);
                rc = rs_write_raw(path, buf, sizeof(uint32_t), nkeep);
            }
            if (recon != NULL && rc == 0) {
                const uint32_t *in = (const uint32_t *)(const void *)recon;
                for (size_t r = 0; r < nkeep; r++)
                    buf[r] = in[ssrc[r]];
                snprintf(path, sizeof path,
                         "%s/ribbon_reconstruction_component.i32", dir);
                rc = rs_write_raw(path, buf, sizeof(uint32_t), nkeep);
            }
            Arena_restore(arena, smark);
        }
        if (rc == 0)
            fprintf(stderr,
                    "  section %02d label %d: COARSE %ldx%ld -> %zu verts, "
                    "%zu faces, u=[%.1f,%.1f], sidecars(phase=%s "
                    "recon=%s) -> %s\n",
                    index, label, ku, kv, nkeep, snf, u_lo, u_hi,
                    phase != NULL ? "yes" : "NO",
                    recon != NULL ? "yes" : "no", dir);
        else
            fprintf(stderr,
                    "  section %02d label %d: COARSE WRITE FAILED under "
                    "%s\n", index, label, dir);
    }
    Arena_restore(arena, mark);
    return rc;
}

/* Extract one label into <out_dir>/secNN/.  Returns 0 on success. */
static int rs_extract(Arena_T arena, const MeshBinData *m,
                      const int32_t *labels, const float *phase,
                      const int32_t *recon, int32_t label, int index,
                      const char *out_dir)
{
    Arena_Mark mark = Arena_save(arena);
    int32_t *remap = (int32_t *)ARENA_ALLOC(
        arena, (size_t)(m->nv * sizeof(int32_t)));
    size_t snv = 0, snf = 0;
    float *sverts = NULL, *suv = NULL;
    int32_t *sfaces = NULL;
    char dir[1800], path[2048];
    double u_lo = 1e300, u_hi = -1e300;
    int rc = -1;
    for (size_t i = 0; i < m->nv; i++) remap[i] = -1;
    for (size_t f = 0; f < m->nf; f++) {
        const int32_t *fv = &m->faces[f * 3];
        if (labels[fv[0]] == label && labels[fv[1]] == label &&
            labels[fv[2]] == label)
            snf++;
    }
    if (snf == 0) {
        fprintf(stderr, "  section %02d label %d: no faces; skipped\n",
                index, label);
        Arena_restore(arena, mark);
        return -1;
    }
    sfaces = (int32_t *)ARENA_ALLOC(
        arena, (size_t)(snf * 3 * sizeof(int32_t)));
    snf = 0;
    for (size_t f = 0; f < m->nf; f++) {
        const int32_t *fv = &m->faces[f * 3];
        if (labels[fv[0]] != label || labels[fv[1]] != label ||
            labels[fv[2]] != label)
            continue;
        for (int k = 0; k < 3; k++) {
            int32_t v = fv[k];
            if (remap[v] < 0) remap[v] = (int32_t)snv++;
            sfaces[snf * 3 + k] = remap[v];
        }
        snf++;
    }
    sverts = (float *)ARENA_ALLOC(arena, (size_t)(snv * 3 * sizeof(float)));
    suv = m->uv != NULL
        ? (float *)ARENA_ALLOC(arena, (size_t)(snv * 2 * sizeof(float)))
        : NULL;
    for (size_t i = 0; i < m->nv; i++) {
        int32_t r = remap[i];
        if (r < 0) continue;
        for (int k = 0; k < 3; k++)
            sverts[(size_t)r * 3 + k] = m->verts[i * 3 + k];
        if (suv != NULL) {
            suv[(size_t)r * 2 + 0] = m->uv[i * 2 + 0];
            suv[(size_t)r * 2 + 1] = m->uv[i * 2 + 1];
            if (m->uv[i * 2] < u_lo) u_lo = m->uv[i * 2];
            if (m->uv[i * 2] > u_hi) u_hi = m->uv[i * 2];
        }
    }
    snprintf(dir, sizeof dir, "%s/sec%02d", out_dir, index);
    ves_mkdir(dir);
    snprintf(path, sizeof path, "%s/section.vmesh", dir);
    rc = MeshBin_write(path, sverts, snv, sfaces, snf, suv);
    if (rc == 0 && phase != NULL)
        rc = rs_slice_sidecar(arena, dir, "ribbon_phase.f32", phase, remap,
                              m->nv, snv);
    if (rc == 0)
        rc = rs_slice_sidecar(arena, dir, "ribbon_material_identity.i32",
                              labels, remap, m->nv, snv);
    if (rc == 0 && recon != NULL)
        rc = rs_slice_sidecar(arena, dir,
                              "ribbon_reconstruction_component.i32", recon,
                              remap, m->nv, snv);
    if (rc == 0)
        fprintf(stderr,
                "  section %02d label %d: %zu verts, %zu faces, "
                "u=[%.1f,%.1f], sidecars(phase=%s recon=%s) -> %s\n",
                index, label, snv, snf, u_lo, u_hi,
                phase != NULL ? "yes" : "NO",
                recon != NULL ? "yes" : "no", dir);
    else
        fprintf(stderr, "  section %02d label %d: WRITE FAILED under %s\n",
                index, label, dir);
    Arena_restore(arena, mark);
    return rc;
}

/* Prolong solved u from a coarse lattice back onto the fine section:
 * every fine vertex knows its (col,row); the coarse solve holds u at
 * (col/KU, row/KV) nodes; bilinear interpolation over the enclosing coarse
 * cell (falling back to the nearest populated corner) transfers the solved
 * parameterization exactly -- solve at 5%, bake at 100%. v is kept from
 * the fine mesh (rows never move). */
static int rs_prolong(Arena_T arena, const char *fine_path,
                      const char *coarse_path, const char *out_path,
                      long ku, long kv, double cell_du, double cell_dv)
{
    MeshBinData fine, coarse;
    size_t cap = 1;
    int64_t *hkey = NULL;
    int32_t *hidx = NULL;
    float *fuv = NULL;
    size_t missing = 0;
    memset(&fine, 0, sizeof fine);
    memset(&coarse, 0, sizeof coarse);
    if (MeshBin_read_arena(arena, fine_path, &fine) != 0 ||
        fine.uv == NULL) {
        fprintf(stderr, "ribbon_sections: cannot read fine %s\n", fine_path);
        return -1;
    }
    if (MeshBin_read_arena(arena, coarse_path, &coarse) != 0 ||
        coarse.uv == NULL) {
        fprintf(stderr, "ribbon_sections: cannot read coarse %s\n",
                coarse_path);
        return -1;
    }
    /* hash coarse nodes by ORIGINAL (col,row) recovered from v (unchanged)
     * and the coarse mesh's own stored ORIGINAL u?  No -- the coarse solve
     * rewrote u.  The invariant that survives the solve is v (rows) plus
     * the coarse mesh's vertex ORDER... which we did not persist.  What
     * does survive: the coarse VERTEX POSITIONS are exact fine vertices, so
     * match by 3D position hash (quantized to 1/64 vox). */
    while (cap < 2 * coarse.nv + 16) cap <<= 1;
    hkey = (int64_t *)ARENA_ALLOC(arena, (size_t)(cap * sizeof(int64_t)));
    hidx = (int32_t *)ARENA_ALLOC(arena, (size_t)(cap * sizeof(int32_t)));
    for (size_t h = 0; h < cap; h++) hkey[h] = -1;
    for (size_t i = 0; i < coarse.nv; i++) {
        int64_t kz = (int64_t)llround((double)coarse.verts[i * 3 + 0] * 64.0);
        int64_t ky = (int64_t)llround((double)coarse.verts[i * 3 + 1] * 64.0);
        int64_t kx = (int64_t)llround((double)coarse.verts[i * 3 + 2] * 64.0);
        int64_t key = (kz * 73856093) ^ (ky * 19349663) ^ (kx * 83492791);
        size_t h = (size_t)(((uint64_t)key * 2654435761u) & (cap - 1));
        for (;; h = (h + 1) & (cap - 1)) {
            if (hkey[h] < 0) {
                hkey[h] = key;
                hidx[h] = (int32_t)i;
                break;
            }
            if (hkey[h] == key) break;
        }
    }
    /* Every fine vertex: locate its coarse lattice cell by (col,row), find
     * the up-to-4 coarse corners AS FINE VERTICES (they are), read their
     * SOLVED u from the coarse mesh via the position hash, and blend. */
    fuv = (float *)ARENA_ALLOC(arena, (size_t)(fine.nv * 2 * sizeof(float)));
    {
        /* index fine vertices by (col,row) so corner lookups are exact */
        size_t fcap = 1;
        int64_t *fkey = NULL;
        int32_t *fidx = NULL;
        while (fcap < 2 * fine.nv + 16) fcap <<= 1;
        fkey = (int64_t *)ARENA_ALLOC(arena, (size_t)(fcap * sizeof(int64_t)));
        fidx = (int32_t *)ARENA_ALLOC(arena, (size_t)(fcap * sizeof(int32_t)));
        for (size_t h = 0; h < fcap; h++) fkey[h] = -1;
        for (size_t i = 0; i < fine.nv; i++) {
            long col = lround((double)fine.uv[i * 2 + 0] / cell_du);
            long row = lround((double)fine.uv[i * 2 + 1] / cell_dv);
            int64_t key = (((int64_t)row + (1 << 20)) << 32) |
                          ((int64_t)col + (1 << 30));
            size_t h = (size_t)(((uint64_t)key * 2654435761u) & (fcap - 1));
            for (;; h = (h + 1) & (fcap - 1)) {
                if (fkey[h] < 0) {
                    fkey[h] = key;
                    fidx[h] = (int32_t)i;
                    break;
                }
                if (fkey[h] == key) break;
            }
        }
        for (size_t i = 0; i < fine.nv; i++) {
            long col = lround((double)fine.uv[i * 2 + 0] / cell_du);
            long row = lround((double)fine.uv[i * 2 + 1] / cell_dv);
            long c0 = (col / ku) * ku, r0 = (row / kv) * kv;
            double fu = ku > 1 ? (double)(col - c0) / (double)ku : 0.0;
            double fv = kv > 1 ? (double)(row - r0) / (double)kv : 0.0;
            double usum = 0.0, wsum = 0.0;
            for (int dc = 0; dc <= 1; dc++)
                for (int dr = 0; dr <= 1; dr++) {
                    long cc = c0 + dc * ku, rr = r0 + dr * kv;
                    int64_t key = (((int64_t)rr + (1 << 20)) << 32) |
                                  ((int64_t)cc + (1 << 30));
                    size_t h = (size_t)(((uint64_t)key * 2654435761u) &
                                        (fcap - 1));
                    int32_t fv_idx = -1;
                    for (;; h = (h + 1) & (fcap - 1)) {
                        if (fkey[h] < 0) break;
                        if (fkey[h] == key) {
                            fv_idx = fidx[h];
                            break;
                        }
                    }
                    if (fv_idx < 0) continue;
                    {
                        /* corner fine vertex -> coarse solved u by 3D hash */
                        int64_t kz = (int64_t)llround(
                            (double)fine.verts[(size_t)fv_idx * 3 + 0] * 64.0);
                        int64_t ky = (int64_t)llround(
                            (double)fine.verts[(size_t)fv_idx * 3 + 1] * 64.0);
                        int64_t kx = (int64_t)llround(
                            (double)fine.verts[(size_t)fv_idx * 3 + 2] * 64.0);
                        int64_t pkey = (kz * 73856093) ^ (ky * 19349663) ^
                                       (kx * 83492791);
                        size_t ph = (size_t)(((uint64_t)pkey * 2654435761u) &
                                             (cap - 1));
                        int32_t ci = -1;
                        for (;; ph = (ph + 1) & (cap - 1)) {
                            if (hkey[ph] < 0) break;
                            if (hkey[ph] == pkey) {
                                ci = hidx[ph];
                                break;
                            }
                        }
                        if (ci < 0) continue;
                        {
                            double w = (dc ? fu : 1.0 - fu) *
                                       (dr ? fv : 1.0 - fv);
                            if (w <= 1e-12) w = 1e-12;
                            usum += w * (double)coarse.uv[(size_t)ci * 2];
                            wsum += w;
                        }
                    }
                }
            if (wsum > 0.0) {
                fuv[i * 2 + 0] = (float)(usum / wsum);
            } else {
                fuv[i * 2 + 0] = fine.uv[i * 2 + 0];
                missing++;
            }
            fuv[i * 2 + 1] = fine.uv[i * 2 + 1];   /* v: rows never move */
        }
    }
    if (MeshBin_write(out_path, fine.verts, fine.nv, fine.faces, fine.nf,
                      fuv) != 0) {
        fprintf(stderr, "ribbon_sections: cannot write %s\n", out_path);
        return -1;
    }
    fprintf(stderr,
            "ribbon_sections: prolonged %zu fine verts from %zu coarse "
            "nodes (%zu kept original u; %.2f%%) -> %s\n",
            fine.nv, coarse.nv, missing,
            100.0 * (double)missing / (double)fine.nv, out_path);
    return 0;
}

static int rs_selftest(void)
{
    /* two labeled patches; biggest extracts with UV and sliced sidecars */
    Arena_T arena = Arena_new();
    float verts[7 * 3] = { 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1,
                           5, 0, 0, 5, 1, 0, 5, 0, 1 };
    float uv[7 * 2] = { 0, 0, 1, 0, 0, 1, 1, 1, 10, 0, 11, 0, 10, 1 };
    int32_t faces[3 * 3] = { 0, 1, 2, 1, 3, 2, 4, 5, 6 };
    int32_t labels[7] = { 7, 7, 7, 7, 3, 3, 3 };
    float phase[7] = { 0.5f, 0.5f, 0.5f, 0.5f, 9.f, 9.f, 9.f };
    MeshBinData m;
    char tmp[1024];
    char out[1400];
    int fails = 0;
    memset(&m, 0, sizeof m);
    m.verts = verts;
    m.uv = uv;
    m.faces = faces;
    m.nv = 7;
    m.nf = 3;
    if (ves_temp_dir(tmp, sizeof tmp) != 0)
        snprintf(tmp, sizeof tmp, ".");
    snprintf(out, sizeof out, "%s/ribbon_sections_selftest", tmp);
    ves_mkdir(out);
    if (rs_extract(arena, &m, labels, phase, NULL, 7, 0, out) != 0) fails++;
    {
        char path[2048];
        MeshBinData s;
        float *ph = NULL;
        memset(&s, 0, sizeof s);
        snprintf(path, sizeof path, "%s/sec00/section.vmesh", out);
        if (MeshBin_read_arena(arena, path, &s) != 0 || s.nv != 4 ||
            s.nf != 2 || s.uv == NULL) {
            fprintf(stderr, "[sections selftest] extract mismatch\n");
            fails++;
        }
        remove(path);
        snprintf(path, sizeof path, "%s/sec00/ribbon_phase.f32", out);
        ph = (float *)rs_load_raw(arena, path, sizeof(float), 4, "phase");
        if (ph == NULL || ph[0] < 0.49f || ph[0] > 0.51f) {
            fprintf(stderr, "[sections selftest] phase sidecar mismatch\n");
            fails++;
        }
        remove(path);
        snprintf(path, sizeof path,
                 "%s/sec00/ribbon_material_identity.i32", out);
        remove(path);
    }
    Arena_dispose(&arena);
    fprintf(stderr, "ribbon_sections --selftest: %s\n",
            fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    MeshBinData m;
    const char *mesh_path = NULL, *labels_path = NULL, *out_dir = NULL;
    const char *phase_path = NULL, *recon_path = NULL;
    int32_t *labels = NULL, *recon = NULL;
    float *phase = NULL;
    long top = 1, coarsen_u = 1, coarsen_v = 1;
    double min_frac = 0.0, cell_du = 2.0, cell_dv = 2.0;
    int census_only = 0;
    setvbuf(stderr, NULL, _IONBF, 0);
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return rs_selftest();
    if (argc >= 5 && strcmp(argv[1], "--prolong") == 0) {
        long pku = 5, pkv = 4;
        double pdu = 2.0, pdv = 2.0;
        Arena_T parena = Arena_new();
        int prc = 1;
        for (int i = 5; i < argc; i++) {
            if (strcmp(argv[i], "--coarsen") == 0 && i + 2 < argc) {
                pku = strtol(argv[++i], NULL, 10);
                pkv = strtol(argv[++i], NULL, 10);
            } else if (strcmp(argv[i], "--cell") == 0 && i + 2 < argc) {
                pdu = strtod(argv[++i], NULL);
                pdv = strtod(argv[++i], NULL);
            }
        }
        if (parena != NULL)
            prc = rs_prolong(parena, argv[2], argv[3], argv[4],
                             pku, pkv, pdu, pdv) == 0 ? 0 : 1;
        if (parena != NULL) Arena_dispose(&parena);
        return prc;
    }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--top") == 0 && i + 1 < argc)
            top = strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--min-frac") == 0 && i + 1 < argc)
            min_frac = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--phase") == 0 && i + 1 < argc)
            phase_path = argv[++i];
        else if (strcmp(argv[i], "--recon") == 0 && i + 1 < argc)
            recon_path = argv[++i];
        else if (strcmp(argv[i], "--coarsen") == 0 && i + 2 < argc) {
            coarsen_u = strtol(argv[++i], NULL, 10);
            coarsen_v = strtol(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--cell") == 0 && i + 2 < argc) {
            cell_du = strtod(argv[++i], NULL);
            cell_dv = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--census") == 0)
            census_only = 1;
        else if (mesh_path == NULL)
            mesh_path = argv[i];
        else if (labels_path == NULL)
            labels_path = argv[i];
        else if (out_dir == NULL)
            out_dir = argv[i];
    }
    if (coarsen_u < 1) coarsen_u = 1;
    if (coarsen_v < 1) coarsen_v = 1;
    if (mesh_path == NULL || labels_path == NULL ||
        (out_dir == NULL && !census_only)) {
        fprintf(stderr,
                "usage: ribbon_sections <ribbon.vmesh> <labels.i32> "
                "<out_dir> [--top N] [--min-frac F] [--census]\n"
                "                       [--phase <phase.f32>] "
                "[--recon <recon.i32>]\n"
                "                       [--coarsen KU KV] [--cell DU DV]\n"
                "       ribbon_sections --prolong <fine.vmesh> "
                "<coarse_solved.vmesh> <out.vmesh>\n"
                "                       [--coarsen KU KV] [--cell DU DV]\n"
                "       ribbon_sections --selftest\n");
        return 2;
    }
    arena = Arena_new();
    if (arena == NULL) return 1;
    memset(&m, 0, sizeof m);
    if (MeshBin_read_arena(arena, mesh_path, &m) != 0) {
        fprintf(stderr, "ribbon_sections: cannot read %s\n", mesh_path);
        return 1;
    }
    fprintf(stderr, "ribbon_sections: %s: %zu verts, %zu faces, uv=%s\n",
            mesh_path, m.nv, m.nf, m.uv != NULL ? "yes" : "NO");
    labels = (int32_t *)rs_load_raw(arena, labels_path, sizeof(int32_t),
                                    m.nv, "labels");
    if (labels == NULL) return 1;
    if (phase_path != NULL) {
        phase = (float *)rs_load_raw(arena, phase_path, sizeof(float),
                                     m.nv, "phases");
        if (phase == NULL) return 1;
    }
    if (recon_path != NULL) {
        recon = (int32_t *)rs_load_raw(arena, recon_path, sizeof(int32_t),
                                       m.nv, "recon ids");
        if (recon == NULL) return 1;
    }

    /* census */
    int32_t max_label = -1;
    size_t negatives = 0;
    for (size_t i = 0; i < m.nv; i++) {
        if (labels[i] > max_label) max_label = labels[i];
        if (labels[i] < 0) negatives++;
    }
    if (max_label < 0) {
        fprintf(stderr, "ribbon_sections: no non-negative labels\n");
        return 1;
    }
    size_t nlab = (size_t)max_label + 1;
    size_t *count = (size_t *)ARENA_CALLOC(arena, nlab, sizeof(size_t));
    for (size_t i = 0; i < m.nv; i++)
        if (labels[i] >= 0) count[labels[i]]++;
    RsCount *rank = (RsCount *)ARENA_ALLOC(
        arena, (size_t)(nlab * sizeof(RsCount)));
    size_t nused = 0;
    for (size_t l = 0; l < nlab; l++) {
        if (count[l] == 0) continue;
        rank[nused].label = (int32_t)l;
        rank[nused].nv = count[l];
        nused++;
    }
    qsort(rank, nused, sizeof(RsCount), rs_cmp_count_desc);
    fprintf(stderr,
            "ribbon_sections: %zu populated label(s), %zu unlabeled "
            "vert(s); top of census:\n",
            nused, negatives);
    for (size_t k = 0; k < nused && k < 12; k++)
        fprintf(stderr, "  #%zu label %d: %zu verts (%.2f%%)\n",
                k, rank[k].label, rank[k].nv,
                100.0 * (double)rank[k].nv / (double)m.nv);

    /* selection */
    size_t keep = 0, kept_verts = 0;
    for (size_t k = 0; k < nused; k++) {
        int in = min_frac > 0.0
               ? ((double)rank[k].nv >= min_frac * (double)m.nv)
               : (k < (size_t)(top > 0 ? top : 0));
        if (!in) break;
        keep++;
        kept_verts += rank[k].nv;
    }
    fprintf(stderr,
            "ribbon_sections: keeping %zu section(s) (%zu verts, %.2f%%); "
            "gutter: %zu label(s), %zu verts (%.2f%%)\n",
            keep, kept_verts, 100.0 * (double)kept_verts / (double)m.nv,
            nused - keep + (negatives ? 1u : 0u),
            m.nv - kept_verts,
            100.0 * (double)(m.nv - kept_verts) / (double)m.nv);
    if (census_only) return 0;
    ves_mkdir(out_dir);
    for (size_t k = 0; k < keep; k++) {
        int erc = coarsen_u > 1 || coarsen_v > 1
                ? rs_extract_coarse(arena, &m, labels, phase, recon,
                                    rank[k].label, (int)k, out_dir,
                                    coarsen_u, coarsen_v, cell_du, cell_dv)
                : rs_extract(arena, &m, labels, phase, recon, rank[k].label,
                             (int)k, out_dir);
        if (erc != 0) return 1;
    }
    Arena_dispose(&arena);
    return 0;
}
