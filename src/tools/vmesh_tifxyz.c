/* vmesh_tifxyz -- export the pieces of a flattened VMESH as VC3D tifxyz segments.
 *
 *   vmesh_tifxyz <in.vmesh> <out_dir> [--prefix S] [--top K] [--min-area A] [--max-dim N]
 *   vmesh_tifxyz --selftest
 *
 * A piece is a connected component of the face graph (faces sharing a vertex), joined with any
 * component whose UV footprint comes within --uv-join voxels (charts of one assembled piece abut in
 * UV without sharing vertices; separately packed pieces keep a wide gap).
 * Each exported piece becomes <out_dir>/<prefix>_p<rank>/ holding x/y/z/mask.tif +
 * meta.json at one grid cell per voxel of its scientific UV (TifXYZ_write_double:
 * contested multi-covers are invalidated, never blended).  Pieces are ranked by UV
 * area; <out_dir>/tifxyz_index.json lists every exported piece.  Geometry and UVs
 * are copied, never modified. */
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/union_find.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../flatten/tifxyz_write.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TIFXYZ_INVALID
#define TIFXYZ_INVALID (-1.0f)   /* the writer's invalid-cell value (flatten/tifxyz_write.c) */
#endif

typedef struct {
    int32_t root;
    size_t  nf;
    double  area;
} Piece;

static int   g_seam_reach = 2;       /* --seam-fill: widest gap (cells) closed between opposite neighbours */
static float g_seam_tol = 4.5f;      /* --seam-tol: max 3-D distance (vox) between those neighbours */

static int cmp_piece(const void *a, const void *b)
{
    const Piece *p = (const Piece *)a, *q = (const Piece *)b;
    if (p->area > q->area) return -1;
    if (p->area < q->area) return 1;
    return (p->root > q->root) - (p->root < q->root);
}

/* Write s as a JSON string literal (Windows paths carry backslashes). */
static void json_str(FILE *f, const char *s)
{
    fputc('"', f);
    for (; *s; s++) {
        if (*s == '"' || *s == '\\') fputc('\\', f);
        if ((unsigned char)*s < 0x20) fprintf(f, "\\u%04x", (unsigned)(unsigned char)*s);
        else fputc(*s, f);
    }
    fputc('"', f);
}

static double uv_area(const double *uv, const int32_t *f)
{
    double ax = uv[f[1] * 2] - uv[f[0] * 2], ay = uv[f[1] * 2 + 1] - uv[f[0] * 2 + 1];
    double bx = uv[f[2] * 2] - uv[f[0] * 2], by = uv[f[2] * 2 + 1] - uv[f[0] * 2 + 1];
    return 0.5 * fabs(ax * by - ay * bx);
}

/* Close thin seam gaps: charts of neighbouring cubes abut in UV without sharing vertices, so a row
 * of pixel centres between them can stay uncovered.  An invalid cell whose nearest valid cells on
 * opposite sides (left/right or up/down) are at most `reach` cells apart and at most `tol` voxels
 * apart in 3-D takes their distance-weighted interpolation; wider or 3-D-discontinuous holes stay
 * invalid.  Reads a snapshot of the mask so fills never chain.  Returns the number of cells filled. */
static size_t seam_fill(Arena_T arena, float *x, float *y, float *z, uint8_t *mask, int W, int H, int reach,
                        float tol)
{
    size_t n = (size_t)W * (size_t)H, filled = 0;
    uint8_t *m0 = (uint8_t *)ARENA_ALLOC(arena, n);
    memcpy(m0, mask, n);
    for (int r = 0; r < H; r++)
        for (int c = 0; c < W; c++) {
            size_t i = (size_t)r * (size_t)W + (size_t)c;
            if (m0[i]) continue;
            int best_span = reach + 2;
            size_t ba = 0, bb = 0;
            int da_best = 0, db_best = 0;
            for (int dir = 0; dir < 2; dir++) {
                int dc = dir == 0 ? 1 : 0, dr = dir == 0 ? 0 : 1;
                int da = 0, db = 0;
                size_t ia = 0, ib = 0;
                for (int k = 1; k <= reach && !da; k++) {
                    int rr = r - k * dr, cc = c - k * dc;
                    if (rr < 0 || cc < 0) break;
                    size_t j = (size_t)rr * (size_t)W + (size_t)cc;
                    if (m0[j]) { da = k; ia = j; }
                }
                for (int k = 1; k <= reach && !db; k++) {
                    int rr = r + k * dr, cc = c + k * dc;
                    if (rr >= H || cc >= W) break;
                    size_t j = (size_t)rr * (size_t)W + (size_t)cc;
                    if (m0[j]) { db = k; ib = j; }
                }
                if (!da || !db || da + db - 1 > reach || da + db >= best_span) continue;
                float ex = x[ia] - x[ib], ey = y[ia] - y[ib], ez = z[ia] - z[ib];
                if (ex * ex + ey * ey + ez * ez > tol * tol) continue;
                best_span = da + db;
                ba = ia; bb = ib; da_best = da; db_best = db;
            }
            if (best_span > reach + 1) continue;
            float wa = (float)db_best / (float)(da_best + db_best), wb = 1.0f - wa;
            x[i] = wa * x[ba] + wb * x[bb];
            y[i] = wa * y[ba] + wb * y[bb];
            z[i] = wa * z[ba] + wb * z[bb];
            mask[i] = 255;
            filled++;
        }
    return filled;
}

/* Re-open a written tifxyz directory, close its seam gaps and write it back.  Returns cells filled
 * or (size_t)-1 on I/O failure. */
static size_t seam_fill_dir(const char *dir, int reach, float tol)
{
    if (reach <= 0) return 0;
    Arena_T a = Arena_new();
    char px[1100], py[1100], pz[1100], pm[1100];
    snprintf(px, sizeof px, "%s/x.tif", dir);
    snprintf(py, sizeof py, "%s/y.tif", dir);
    snprintf(pz, sizeof pz, "%s/z.tif", dir);
    snprintf(pm, sizeof pm, "%s/mask.tif", dir);
    float *x = NULL, *y = NULL, *z = NULL;
    int W = 0, H = 0, W2 = 0, H2 = 0, W3 = 0, H3 = 0;
    size_t filled = (size_t)-1;
    if (TiffIO_load_float2d(a, px, &x, &W, &H) == 0 && TiffIO_load_float2d(a, py, &y, &W2, &H2) == 0 &&
        TiffIO_load_float2d(a, pz, &z, &W3, &H3) == 0 && W == W2 && W == W3 && H == H2 && H == H3) {
        size_t n = (size_t)W * (size_t)H;
        uint8_t *mask = (uint8_t *)ARENA_ALLOC(a, n);
        for (size_t i = 0; i < n; i++) mask[i] = (x[i] != TIFXYZ_INVALID || y[i] != TIFXYZ_INVALID) ? 255 : 0;
        /* seam crossings have no valid row/column neighbour until their seams are closed: up to 3 passes */
        filled = 0;
        for (int pass = 0; pass < 3; pass++) {
            size_t f = seam_fill(a, x, y, z, mask, W, H, reach, tol);
            filled += f;
            if (f == 0) break;
        }
        if (filled > 0 && (TiffIO_save_float2d(px, x, W, H) != 0 || TiffIO_save_float2d(py, y, W, H) != 0 ||
                           TiffIO_save_float2d(pz, z, W, H) != 0 || TiffIO_save(pm, mask, 1, H, W) != 0))
            filled = (size_t)-1;
    }
    Arena_dispose(&a);
    return filled;
}

/* Export the pieces of (verts, faces, uv).  Returns the number of pieces written, -1 on error. */
/* Join components whose UV footprints come within `cell` voxels: charts of one assembled piece
 * abut in UV while separately packed pieces keep a wide gap.  Face UV boxes are marked on a
 * coarse grid; a cell claimed by two components unions them. */
static int uv_join(Arena_T arena, UnionFind *uf, const int32_t *faces, size_t nf, const double *uv, size_t nv,
                   double cell)
{
    double umin = INFINITY, vmin = INFINITY, umax = -INFINITY, vmax = -INFINITY;
    for (size_t i = 0; i < nv; i++) {
        double u = uv[i * 2], v = uv[i * 2 + 1];
        if (u < umin) umin = u;
        if (u > umax) umax = u;
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    size_t W = (size_t)((umax - umin) / cell) + 2, H = (size_t)((vmax - vmin) / cell) + 2;
    if (W * H > ((size_t)1 << 31)) return -1;
    int32_t *owner = (int32_t *)ARENA_ALLOC(arena, W * H * sizeof(int32_t));
    for (size_t i = 0; i < W * H; i++) owner[i] = -1;
    for (size_t f = 0; f < nf; f++) {
        const int32_t *t = faces + f * 3;
        double lo_u = INFINITY, hi_u = -INFINITY, lo_v = INFINITY, hi_v = -INFINITY;
        for (int c = 0; c < 3; c++) {
            double u = uv[(size_t)t[c] * 2], v = uv[(size_t)t[c] * 2 + 1];
            if (u < lo_u) lo_u = u;
            if (u > hi_u) hi_u = u;
            if (v < lo_v) lo_v = v;
            if (v > hi_v) hi_v = v;
        }
        size_t c0 = (size_t)((lo_u - umin) / cell), c1 = (size_t)((hi_u - umin) / cell);
        size_t r0 = (size_t)((lo_v - vmin) / cell), r1 = (size_t)((hi_v - vmin) / cell);
        for (size_t r = r0; r <= r1 && r < H; r++)
            for (size_t c = c0; c <= c1 && c < W; c++) {
                int32_t *o = owner + r * W + c;
                if (*o < 0) *o = t[0];
                else uf_union(uf, *o, t[0]);
            }
    }
    return 0;
}

static int export_pieces(Arena_T arena, const float *verts, size_t nv, const int32_t *faces, size_t nf,
                         const double *uv, const char *out_dir, const char *prefix, int top, double min_area,
                         int max_dim, double join_cell, FILE *index)
{
    if (nv == 0 || nf == 0 || nv > (size_t)INT32_MAX) return 0;
    UnionFind uf = UF_new(arena, (int32_t)nv);
    for (size_t f = 0; f < nf; f++) {
        const int32_t *t = faces + f * 3;
        if (t[0] < 0 || t[1] < 0 || t[2] < 0 || (size_t)t[0] >= nv || (size_t)t[1] >= nv || (size_t)t[2] >= nv)
            return -1;
        uf_union(&uf, t[0], t[1]);
        uf_union(&uf, t[0], t[2]);
    }
    if (join_cell > 0.0 && uv_join(arena, &uf, faces, nf, uv, nv, join_cell) != 0) return -1;
    /* per-root face counts and UV areas */
    size_t *slot = (size_t *)ARENA_ALLOC(arena, nv * sizeof(size_t));
    for (size_t i = 0; i < nv; i++) slot[i] = (size_t)-1;
    Piece *pieces = (Piece *)ARENA_ALLOC(arena, nv * sizeof(Piece));
    size_t np = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t r = uf_find(&uf, faces[f * 3]);
        if (slot[r] == (size_t)-1) {
            slot[r] = np;
            pieces[np].root = r;
            pieces[np].nf = 0;
            pieces[np].area = 0.0;
            np++;
        }
        pieces[slot[r]].nf++;
        pieces[slot[r]].area += uv_area(uv, faces + f * 3);
    }
    qsort(pieces, np, sizeof(Piece), cmp_piece);
    int32_t *froot = (int32_t *)ARENA_ALLOC(arena, nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) froot[f] = uf_find(&uf, faces[f * 3]);
    int32_t *vmap = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    int written = 0;
    for (size_t k = 0; k < np && (top <= 0 || written < top); k++) {
        if (pieces[k].area < min_area) break;
        Arena_T scratch = Arena_new();
        size_t pnf = pieces[k].nf, pnv = 0;
        int32_t *pf = (int32_t *)ARENA_ALLOC(scratch, pnf * 3 * sizeof(int32_t));
        for (size_t i = 0; i < nv; i++) vmap[i] = -1;
        size_t fi = 0;
        for (size_t f = 0; f < nf; f++) {
            if (froot[f] != pieces[k].root) continue;
            for (int c = 0; c < 3; c++) {
                int32_t v = faces[f * 3 + c];
                if (vmap[v] < 0) vmap[v] = (int32_t)pnv++;
                pf[fi * 3 + c] = vmap[v];
            }
            fi++;
        }
        float *pv = (float *)ARENA_ALLOC(scratch, pnv * 3 * sizeof(float));
        double *puv = (double *)ARENA_ALLOC(scratch, pnv * 2 * sizeof(double));
        for (size_t i = 0; i < nv; i++) {
            if (vmap[i] < 0) continue;
            memcpy(pv + (size_t)vmap[i] * 3, verts + i * 3, 3 * sizeof(float));
            puv[(size_t)vmap[i] * 2] = uv[i * 2];
            puv[(size_t)vmap[i] * 2 + 1] = uv[i * 2 + 1];
        }
        char name[256], dir[1024];
        snprintf(name, sizeof name, "%s_p%03d", prefix, written);
        snprintf(dir, sizeof dir, "%s/%s", out_dir, name);
        TifXYZOpts opt;
        memset(&opt, 0, sizeof opt);
        opt.px_per_vox = 1.0f;
        opt.max_dim = max_dim;
        TifXYZStats st;
        memset(&st, 0, sizeof st);
        int rc = TifXYZ_write_double(scratch, dir, name, pv, pnv, pf, pnf, puv, &opt, &st);
        if (rc != 0) {
            Arena_dispose(&scratch);
            return -1;
        }
        size_t seam = seam_fill_dir(dir, g_seam_reach, g_seam_tol);
        if (seam == (size_t)-1) {
            Arena_dispose(&scratch);
            return -1;
        }
        if (index)
            fprintf(index, "%s  {\"dir\": \"%s\", \"rank\": %d, \"faces\": %zu, \"verts\": %zu, \"uv_area\": %.3f, "
                           "\"width\": %d, \"height\": %d, \"scale\": [%.6f, %.6f], \"cells_filled\": %zu, "
                           "\"cells_contested\": %zu, \"seam_filled_cells\": %zu, "
                           "\"bbox_xyz\": [[%.1f, %.1f, %.1f], [%.1f, %.1f, %.1f]]}",
                    written ? ",\n" : "", name, written, pnf, pnv, pieces[k].area, st.width, st.height,
                    st.scale_u, st.scale_v, st.cells_filled, st.cells_contested, seam, st.bbox_lo[0], st.bbox_lo[1],
                    st.bbox_lo[2], st.bbox_hi[0], st.bbox_hi[1], st.bbox_hi[2]);
        fprintf(stderr, "[vmesh_tifxyz] %s: %zu faces, uv area %.0f vox^2, grid %dx%d, %zu cells, %zu contested\n",
                name, pnf, pieces[k].area, st.width, st.height, st.cells_filled, st.cells_contested);
        Arena_dispose(&scratch);
        written++;
    }
    fprintf(stderr, "[vmesh_tifxyz] %zu pieces, %d exported\n", np, written);
    return written;
}

static int run(const char *in, const char *out_dir, const char *prefix, int top, double min_area, int max_dim,
               double join_cell)
{
    Arena_T arena = Arena_new();
    MeshBinData m;
    memset(&m, 0, sizeof m);
    if (MeshBin_read_precise_arena(arena, in, &m) != 0 || (m.uv == NULL && m.uv64 == NULL)) {
        fprintf(stderr, "[vmesh_tifxyz] cannot read a UV mesh from %s\n", in);
        Arena_dispose(&arena);
        return 1;
    }
    double *uv = m.uv64;
    if (uv == NULL) {
        uv = (double *)ARENA_ALLOC(arena, m.nv * 2 * sizeof(double));
        for (size_t i = 0; i < m.nv * 2; i++) uv[i] = (double)m.uv[i];
    }
    char path[1024];
    snprintf(path, sizeof path, "%s/tifxyz_index.json", out_dir);
    ves_ensure_parent_dir(path);
    FILE *index = fopen(path, "w");
    if (!index) {
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(index, "{\"source\": ");
    json_str(index, in);
    fprintf(index, ", \"verts\": %zu, \"faces\": %zu, \"pieces\": [\n", m.nv, m.nf);
    int n = export_pieces(arena, m.verts, m.nv, m.faces, m.nf, uv, out_dir, prefix, top, min_area, max_dim, join_cell,
                          index);
    fprintf(index, "\n]}\n");
    fclose(index);
    Arena_dispose(&arena);
    return n < 0;
}

static int seam_fill_selftest(void)
{
    int fails = 0;
    enum { W = 12, H = 3 };
    Arena_T a = Arena_new();
    float x[W * H], y[W * H], z[W * H];
    uint8_t m[W * H];
    /* row 0: 1-cell gap at c=4 between 3-D neighbours (x = c): must close with x = 4
     * row 1: 1-cell gap at c=4 but the right side sits 10 vox away in z: must stay open
     * row 2: 4-cell gap at c=3..6: wider than reach 2, must stay open */
    for (int r = 0; r < H; r++)
        for (int c = 0; c < W; c++) {
            int i = r * W + c;
            int hole = (r < 2 && c == 4) || (r == 2 && c >= 3 && c <= 6);
            x[i] = hole ? TIFXYZ_INVALID : (float)c;
            y[i] = hole ? TIFXYZ_INVALID : (float)r * 100.0f;
            z[i] = hole ? TIFXYZ_INVALID : (r == 1 && c > 4 ? 10.0f : 0.0f);
            m[i] = hole ? 0 : 255;
        }
    size_t n = seam_fill(a, x, y, z, m, W, H, 2, 4.5f);
    if (n != 1 || m[4] != 255 || fabsf(x[4] - 4.0f) > 1e-4f || m[W + 4] != 0 || m[2 * W + 4] != 0) {
        fprintf(stderr, "selftest: seam fill filled %zu (row0 %d x=%.2f, row1 %d, row2 %d)\n", n, m[4], x[4], m[W + 4],
                m[2 * W + 4]);
        fails++;
    }
    Arena_dispose(&a);
    return fails;
}

static int selftest(void)
{
    int fails = TifXYZ_selftest() + seam_fill_selftest();
    /* two separated unit squares of 40x20 vox (2 triangles each) -> two pieces, the larger first */
    Arena_T arena = Arena_new();
    float v[12 * 3] = {
        100, 10, 10,  100, 10, 50,  120, 10, 50,  120, 10, 10,       /* piece A: 40 x 20 */
        100, 60, 10,  100, 60, 20,  110, 60, 20,  110, 60, 10,       /* piece B: 10 x 10, far away in UV */
        100, 10, 52,  100, 10, 60,  120, 10, 60,  120, 10, 52 };     /* chart A2: abuts A in UV (2 vox gap) */
    int32_t f[6 * 3] = { 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 8, 9, 10, 8, 10, 11 };
    double uv[12 * 2] = { 0, 0, 40, 0, 40, 20, 0, 20, 100, 0, 110, 0, 110, 10, 100, 10,
                          42, 0, 50, 0, 50, 20, 42, 20 };
    char dir[512];
    snprintf(dir, sizeof dir, "vmesh_tifxyz_selftest_%d", (int)ves_getpid());
    char path[640];
    snprintf(path, sizeof path, "%s/index.json", dir);
    ves_ensure_parent_dir(path);
    FILE *ix = fopen(path, "w");
    if (ix) {
        fprintf(ix, "{\"source\": ");
        json_str(ix, "C:\\dir \"quoted\"\\x.vmesh");
        fprintf(ix, ", \"pieces\": [\n");
    }
    int n = export_pieces(arena, v, 12, f, 6, uv, dir, "t", 0, 1.0, 8192, 4.0, ix);
    if (ix) {
        fprintf(ix, "\n]}\n");
        fclose(ix);
        /* the index must be JSON: the escaped path round-trips as \\ and \" */
        FILE *r = fopen(path, "rb");
        char buf[256] = { 0 };
        size_t got = r ? fread(buf, 1, sizeof buf - 1, r) : 0;
        if (r) fclose(r);
        if (got == 0 || strstr(buf, "\"C:\\\\dir \\\"quoted\\\"\\\\x.vmesh\"") == NULL) {
            fprintf(stderr, "selftest: index path not JSON-escaped: %.80s\n", buf);
            fails++;
        }
    }
    if (n != 2) { fprintf(stderr, "selftest: expected 2 joined pieces, got %d\n", n); fails++; }
    char dir2[512];
    snprintf(dir2, sizeof dir2, "vmesh_tifxyz_selftest_nojoin_%d", (int)ves_getpid());
    int n2 = export_pieces(arena, v, 12, f, 6, uv, dir2, "t", 0, 1.0, 8192, 0.0, NULL);
    if (n2 != 3) { fprintf(stderr, "selftest: expected 3 unjoined pieces, got %d\n", n2); fails++; }
    float *x = NULL;
    int W = 0, H = 0;
    snprintf(path, sizeof path, "%s/t_p000/x.tif", dir);
    Arena_T a2 = Arena_new();
    if (TiffIO_load_float2d(a2, path, &x, &W, &H) != 0 || W != 51 || H != 21) {
        fprintf(stderr, "selftest: largest piece grid %dx%d, expected 51x21\n", W, H);
        fails++;
    } else {
        /* x.tif <- vertex x (third component): the centre cell lies inside piece A (x in [10,50]) */
        float cx = x[(size_t)10 * (size_t)W + 20];
        if (!(cx > 9.0f && cx < 51.0f)) { fprintf(stderr, "selftest: x.tif centre %.2f\n", cx); fails++; }
    }
    Arena_dispose(&a2);
    Arena_dispose(&arena);
    fprintf(stderr, "vmesh_tifxyz selftest: %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) return selftest() != 0;
    if (argc < 3) goto usage;
    const char *prefix = "piece";
    int top = 16, max_dim = 32768;
    double min_area = 10000.0, join_cell = 4.0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--prefix") && i + 1 < argc) prefix = argv[++i];
        else if (!strcmp(argv[i], "--top") && i + 1 < argc) top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min-area") && i + 1 < argc) min_area = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max-dim") && i + 1 < argc) max_dim = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--uv-join") && i + 1 < argc) join_cell = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seam-fill") && i + 1 < argc) g_seam_reach = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seam-tol") && i + 1 < argc) g_seam_tol = (float)atof(argv[++i]);
        else goto usage;
    }
    if (max_dim < 2 || max_dim > 65535) goto usage;
    return run(argv[1], argv[2], prefix, top, min_area, max_dim, join_cell);
usage:
    fprintf(stderr, "usage: vmesh_tifxyz <in.vmesh> <out_dir> [--prefix S] [--top K] [--min-area A] [--max-dim N]"
                    " [--uv-join F] [--seam-fill N] [--seam-tol F]\n"
                    "       vmesh_tifxyz --selftest\n"
                    "Writes one VC3D tifxyz directory per piece (largest UV area first). A piece is a mesh\n"
                    "component joined with components whose UV footprints come within F vox (default 4;\n"
                    "0 = mesh connectivity only): charts of one packed piece abut, pieces keep a gap.\n"
                    "Seam gaps up to N cells (default 2; 0 = off) between opposite valid cells at most F vox\n"
                    "apart in 3-D (default 4.5) are closed by interpolation; wider holes stay invalid.\n");
    return 2;
}
