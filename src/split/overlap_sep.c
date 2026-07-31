/* overlap_sep.c -- Step 3: Overlap-Detection Multicut Separator
 *
 * Separates overlapping sheets by:
 *   Phase 1: Project to 2D via PCA normal
 *   Phase 2: Build face adjacency via edge hashing
 *   Phase 3: Grid-based 2D triangle overlap detection
 *   Phase 4: Lifted multicut (greedy additive + KL refinement)
 *   Phase 5: Merge small clusters, partition faces by label, extract sub-meshes
 *
 * Based on Limper et al. 2018 Section 4.1.
 */
#include "../common/ves_platform.h"

#include "overlap_sep.h"
#include "overlap_merge.h"
#include "overlap_quality.h"
#include "multicut_wrap.h"
#include "../common/pipeline_constants.h"
#include "../common/pca.h"
#include "../common/csr.h"
#include "../common/gco_wrap.h"
#include "../common/raw_sample.h"
#include "../flatten/surface_snap.h"
#include "../flatten/snap_solve.h"

#include "../common/obj_io.h"

#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/* ===== Debug OBJ output ===== */

static const char *g_debug_dir = NULL;

void OverlapSep_set_debug_dir(const char *dir)
{
    g_debug_dir = dir;
}

/* Compact and write exactly the selected source triangles.  This deliberately
 * does no ring growth, welding, projection change, or topology repair: it is a
 * literal view of the faces named by the production overlap detector. */
static int debug_write_face_subset(Arena_T arena,
                                   const char *path,
                                   const float *verts, size_t nv,
                                   const int32_t *faces, size_t nf,
                                   const uint8_t *keep,
                                   const float color[3],
                                   size_t *out_nv, size_t *out_nf)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, nv, sizeof(*used));
    size_t kept_nf = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!keep[fi]) continue;
        kept_nf++;
        used[faces[fi * 3]] = 1;
        used[faces[fi * 3 + 1]] = 1;
        used[faces[fi * 3 + 2]] = 1;
    }
    if (kept_nf == 0) {
        Arena_restore(arena, mark);
        return -1;
    }

    int32_t *remap = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*remap));
    size_t kept_nv = 0;
    for (size_t vi = 0; vi < nv; vi++)
        remap[vi] = used[vi] ? (int32_t)kept_nv++ : -1;

    float *compact_verts = (float *)ARENA_ALLOC(
        arena, kept_nv * 3 * sizeof(*compact_verts));
    for (size_t vi = 0; vi < nv; vi++) {
        if (remap[vi] < 0) continue;
        memcpy(&compact_verts[(size_t)remap[vi] * 3], &verts[vi * 3],
               3 * sizeof(*compact_verts));
    }
    int32_t *compact_faces = (int32_t *)ARENA_ALLOC(
        arena, kept_nf * 3 * sizeof(*compact_faces));
    size_t cursor = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        if (!keep[fi]) continue;
        compact_faces[cursor * 3] = remap[faces[fi * 3]];
        compact_faces[cursor * 3 + 1] = remap[faces[fi * 3 + 1]];
        compact_faces[cursor * 3 + 2] = remap[faces[fi * 3 + 2]];
        cursor++;
    }
    int rc = color ? ObjIO_write_colored(
                         path, compact_verts, kept_nv,
                         compact_faces, kept_nf, color)
                   : ObjIO_write(path, compact_verts, kept_nv,
                                 compact_faces, kept_nf);
    if (out_nv) *out_nv = kept_nv;
    if (out_nf) *out_nf = kept_nf;
    Arena_restore(arena, mark);
    return rc;
}

static int debug_write_patch_stage(Arena_T arena, const char *path,
                                   const ComponentMesh *mesh,
                                   size_t patch_face_start,
                                   const float subtract_offset[3])
{
    if (!mesh || patch_face_start >= mesh->nf) return -1;
    Arena_Mark mark = Arena_save(arena);
    uint8_t *keep = (uint8_t *)ARENA_CALLOC(
        arena, mesh->nf, sizeof(*keep));
    for (size_t fi = patch_face_start; fi < mesh->nf; fi++) keep[fi] = 1;
    const float *verts = mesh->verts;
    float *shifted = NULL;
    if (subtract_offset) {
        shifted = (float *)ARENA_ALLOC(
            arena, mesh->nv * 3 * sizeof(*shifted));
        for (size_t vi = 0; vi < mesh->nv; vi++)
            for (int d = 0; d < 3; d++)
                shifted[vi * 3 + (size_t)d] =
                    mesh->verts[vi * 3 + (size_t)d] - subtract_offset[d];
        verts = shifted;
    }
    int rc = debug_write_face_subset(
        arena, path, verts, mesh->nv, mesh->faces, mesh->nf,
        keep, NULL, NULL, NULL);
    Arena_restore(arena, mark);
    return rc;
}

/* Label-to-color palette (8 distinct colors) */
static const float g_palette[][3] = {
    {1.0f, 0.2f, 0.2f},  /* red */
    {0.2f, 0.6f, 1.0f},  /* blue */
    {0.2f, 0.9f, 0.3f},  /* green */
    {1.0f, 0.8f, 0.1f},  /* yellow */
    {0.8f, 0.3f, 0.9f},  /* purple */
    {1.0f, 0.5f, 0.1f},  /* orange */
    {0.1f, 0.9f, 0.9f},  /* cyan */
    {0.9f, 0.5f, 0.7f},  /* pink */
};
#define N_PALETTE (sizeof(g_palette)/sizeof(g_palette[0]))

static int debug_write_split_merge_side_by_side(
    Arena_T arena, const char *path, const ComponentMesh *split,
    size_t n_split, const ComponentMesh *merged)
{
    if (!split || n_split == 0 || !merged || merged->nv == 0) return -1;
    Arena_Mark mark = Arena_save(arena);
    size_t split_nv = 0, split_nf = 0;
    float min_x = FLT_MAX, max_x = -FLT_MAX;
    for (size_t i = 0; i < n_split; i++) {
        split_nv += split[i].nv;
        split_nf += split[i].nf;
        for (size_t vi = 0; vi < split[i].nv; vi++) {
            float x = split[i].verts[vi * 3 + 2];
            if (x < min_x) min_x = x;
            if (x > max_x) max_x = x;
        }
    }
    for (size_t vi = 0; vi < merged->nv; vi++) {
        float x = merged->verts[vi * 3 + 2];
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
    }
    float separation = (max_x - min_x) + 12.0f;
    size_t total_nv = split_nv + merged->nv;
    size_t total_nf = split_nf + merged->nf;
    float *verts = (float *)ARENA_ALLOC(
        arena, total_nv * 3 * sizeof(*verts));
    float *colors = (float *)ARENA_ALLOC(
        arena, total_nv * 3 * sizeof(*colors));
    int32_t *faces = (int32_t *)ARENA_ALLOC(
        arena, total_nf * 3 * sizeof(*faces));
    size_t voff = 0, foff = 0;
    for (size_t i = 0; i < n_split; i++) {
        const float *color = g_palette[i % N_PALETTE];
        for (size_t vi = 0; vi < split[i].nv; vi++) {
            memcpy(&verts[(voff + vi) * 3], &split[i].verts[vi * 3],
                   3 * sizeof(*verts));
            memcpy(&colors[(voff + vi) * 3], color, 3 * sizeof(*colors));
        }
        for (size_t fi = 0; fi < split[i].nf; fi++)
            for (int k = 0; k < 3; k++)
                faces[(foff + fi) * 3 + (size_t)k] =
                    (int32_t)voff + split[i].faces[fi * 3 + (size_t)k];
        voff += split[i].nv;
        foff += split[i].nf;
    }
    const float merged_color[3] = {0.20f, 0.90f, 0.35f};
    for (size_t vi = 0; vi < merged->nv; vi++) {
        memcpy(&verts[(voff + vi) * 3], &merged->verts[vi * 3],
               3 * sizeof(*verts));
        verts[(voff + vi) * 3 + 2] += separation;
        memcpy(&colors[(voff + vi) * 3], merged_color,
               3 * sizeof(*colors));
    }
    for (size_t fi = 0; fi < merged->nf; fi++)
        for (int k = 0; k < 3; k++)
            faces[(foff + fi) * 3 + (size_t)k] =
                (int32_t)voff + merged->faces[fi * 3 + (size_t)k];
    int rc = ObjIO_write_per_vertex_color(
        path, verts, total_nv, faces, total_nf, colors);
    Arena_restore(arena, mark);
    return rc;
}

/* ===== Timer ===== */

static double now_sec(void)
{
    return ves_clock_sec();
}

/* ===== Return original mesh unchanged ===== */

static void return_unchanged(Arena_T arena, const ComponentMesh *mesh,
                             ComponentMesh **out_meshes, size_t *out_count)
{
    ComponentMesh *out;
    ARENA_NEW(arena, out);
    *out = *mesh;
    out->self = out;
    *out_meshes = out;
    *out_count = 1;
}

/* ===== qsort comparison functions ===== */

typedef struct {
    uint64_t key;
    int32_t  face;
} EdgeRecord;

static int compare_edge_records(const void *a, const void *b)
{
    const EdgeRecord *ea = (const EdgeRecord *)a;
    const EdgeRecord *eb = (const EdgeRecord *)b;
    if (ea->key < eb->key) return -1;
    if (ea->key > eb->key) return 1;
    return 0;
}

static int compare_u64(const void *a, const void *b)
{
    uint64_t va = *(const uint64_t *)a;
    uint64_t vb = *(const uint64_t *)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

static int compare_float(const void *a, const void *b)
{
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

/* ===== Vertex sharing test ===== */

static int faces_share_vertex(const int32_t *faces, int32_t fa, int32_t fb)
{
    const int32_t *a = faces + fa * 3;
    const int32_t *b = faces + fb * 3;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (a[i] == b[j]) return 1;
    return 0;
}

/* ===== 2D geometry helpers ===== */

static int segments_intersect(double p0x, double p0y, double p1x, double p1y,
                              double q0x, double q0y, double q1x, double q1y)
{
    /* Orientation test: (q - p) x (r - p) */
    double d1x = p1x - p0x, d1y = p1y - p0y;
    double d2x = q1x - q0x, d2y = q1y - q0y;

    double o1 = d1x * (q0y - p0y) - d1y * (q0x - p0x);
    double o2 = d1x * (q1y - p0y) - d1y * (q1x - p0x);
    double o3 = d2x * (p0y - q0y) - d2y * (p0x - q0x);
    double o4 = d2x * (p1y - q0y) - d2y * (p1x - q0x);

    /* Segments straddle each other */
    if (o1 * o2 < -(double)OVERLAP_EDGE_EPS &&
        o3 * o4 < -(double)OVERLAP_EDGE_EPS)
        return 1;

    return 0;
}

static int point_in_triangle(double px, double py,
                             double ax, double ay, double bx, double by,
                             double cx, double cy)
{
    double d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by);
    double d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy);
    double d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay);

    int has_neg = (d1 < -(double)OVERLAP_EDGE_EPS) ||
                  (d2 < -(double)OVERLAP_EDGE_EPS) ||
                  (d3 < -(double)OVERLAP_EDGE_EPS);
    int has_pos = (d1 > (double)OVERLAP_EDGE_EPS) ||
                  (d2 > (double)OVERLAP_EDGE_EPS) ||
                  (d3 > (double)OVERLAP_EDGE_EPS);

    return !(has_neg && has_pos);
}

static int triangles_overlap_2d(const float *proj, const int32_t *faces,
                                int32_t fa, int32_t fb)
{
    double ax[3], ay[3], bx[3], by[3];
    for (int i = 0; i < 3; i++) {
        int32_t va = faces[fa * 3 + i];
        int32_t vb = faces[fb * 3 + i];
        ax[i] = (double)proj[va * 2 + 0]; ay[i] = (double)proj[va * 2 + 1];
        bx[i] = (double)proj[vb * 2 + 0]; by[i] = (double)proj[vb * 2 + 1];
    }

    /* Edge-edge intersection: 3x3 = 9 tests */
    for (int i = 0; i < 3; i++) {
        int i2 = (i + 1) % 3;
        for (int j = 0; j < 3; j++) {
            int j2 = (j + 1) % 3;
            if (segments_intersect(ax[i], ay[i], ax[i2], ay[i2],
                                   bx[j], by[j], bx[j2], by[j2]))
                return 1;
        }
    }

    /* Vertex-in-triangle: 3+3 = 6 tests */
    for (int i = 0; i < 3; i++) {
        if (point_in_triangle(ax[i], ay[i], bx[0], by[0], bx[1], by[1],
                              bx[2], by[2]))
            return 1;
        if (point_in_triangle(bx[i], by[i], ax[0], ay[0], ax[1], ay[1],
                              ax[2], ay[2]))
            return 1;
    }

    return 0;
}

/* ===== Face adjacency via edge hashing ===== */

static void build_face_adjacency(Arena_T arena,
                                 const int32_t *faces, size_t nf,
                                 const float *proj, size_t nv,
                                 int32_t **out_fa, int32_t **out_fb,
                                 double **out_edge_len, size_t *out_n_adj,
                                 double *out_avg_edge_len)
{
    (void)nv;
    size_t n_records = nf * 3;
    EdgeRecord *records = (EdgeRecord *)ARENA_ALLOC(arena,
        (long)(n_records * sizeof(EdgeRecord)));

    /* Build edge records: 3 edges per face */
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t v[3];
        v[0] = faces[fi * 3 + 0];
        v[1] = faces[fi * 3 + 1];
        v[2] = faces[fi * 3 + 2];
        int edges[3][2] = {{0,1}, {1,2}, {2,0}};
        for (int e = 0; e < 3; e++) {
            int32_t va = v[edges[e][0]], vb = v[edges[e][1]];
            int32_t lo = va < vb ? va : vb;
            int32_t hi = va < vb ? vb : va;
            size_t idx = fi * 3 + (size_t)e;
            records[idx].key = ((uint64_t)(uint32_t)lo << 32) |
                               (uint64_t)(uint32_t)hi;
            records[idx].face = (int32_t)fi;
        }
    }

    /* Sort by edge key */
    qsort(records, n_records, sizeof(EdgeRecord), compare_edge_records);

    /* Scan for adjacent face pairs (consecutive records with same key) */
    size_t adj_cap = nf * 3 / 2 + 1;
    int32_t *adj_fa = (int32_t *)ARENA_ALLOC(arena,
        (long)(adj_cap * sizeof(int32_t)));
    int32_t *adj_fb = (int32_t *)ARENA_ALLOC(arena,
        (long)(adj_cap * sizeof(int32_t)));
    double *adj_len = (double *)ARENA_ALLOC(arena,
        (long)(adj_cap * sizeof(double)));
    size_t n_adj = 0;
    double total_edge_len = 0.0;
    size_t total_edges = 0;

    size_t i = 0;
    while (i < n_records) {
        size_t j = i + 1;
        while (j < n_records && records[j].key == records[i].key)
            j++;

        /* Extract shared edge vertices from key */
        int32_t v0 = (int32_t)(records[i].key >> 32);
        int32_t v1 = (int32_t)(records[i].key & 0xFFFFFFFF);
        float dx = proj[v0 * 2 + 0] - proj[v1 * 2 + 0];
        float dy = proj[v0 * 2 + 1] - proj[v1 * 2 + 1];
        double elen = sqrt((double)(dx * dx + dy * dy));
        total_edge_len += elen;
        total_edges++;

        /* Each pair of faces sharing this edge = one adjacency */
        for (size_t a = i; a < j; a++) {
            for (size_t b = a + 1; b < j; b++) {
                if (n_adj < adj_cap) {
                    adj_fa[n_adj] = records[a].face;
                    adj_fb[n_adj] = records[b].face;
                    adj_len[n_adj] = elen;
                    n_adj++;
                }
            }
        }
        i = j;
    }

    double avg = (total_edges > 0) ? total_edge_len / (double)total_edges : 1.0;

    *out_fa = adj_fa;
    *out_fb = adj_fb;
    *out_edge_len = adj_len;
    *out_n_adj = n_adj;
    *out_avg_edge_len = avg;
}

/* ===== Grid-based 2D overlap detection ===== */

static void detect_overlaps(Arena_T arena,
                            const int32_t *faces, size_t nf,
                            const float *proj, size_t nv,
                            int32_t **out_fa, int32_t **out_fb,
                            size_t *out_n_ovl)
{
    (void)nv;

    /* Compute triangle areas in 2D and find median */
    float *areas = (float *)ARENA_ALLOC(arena, (long)(nf * sizeof(float)));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t va = faces[fi * 3 + 0];
        int32_t vb = faces[fi * 3 + 1];
        int32_t vc = faces[fi * 3 + 2];
        float ax = proj[va * 2], ay = proj[va * 2 + 1];
        float bx = proj[vb * 2], by = proj[vb * 2 + 1];
        float cx = proj[vc * 2], cy = proj[vc * 2 + 1];
        areas[fi] = 0.5f * fabsf((bx - ax) * (cy - ay) -
                                  (cx - ax) * (by - ay));
    }
    qsort(areas, nf, sizeof(float), compare_float);
    float median_area = areas[nf / 2];
    if (median_area < 1e-12f) median_area = 1e-6f;

    float cell_size = (float)(OVERLAP_GRID_SCALE * sqrt((double)median_area));
    if (cell_size < 1e-6f) cell_size = 1e-6f;

    /* Compute 2D bounding box of projected vertices */
    float minx = proj[0], maxx = proj[0];
    float miny = proj[1], maxy = proj[1];
    for (size_t fi = 0; fi < nf; fi++) {
        for (int j = 0; j < 3; j++) {
            int32_t vi = faces[fi * 3 + j];
            float x = proj[vi * 2 + 0], y = proj[vi * 2 + 1];
            if (x < minx) minx = x;
            if (x > maxx) maxx = x;
            if (y < miny) miny = y;
            if (y > maxy) maxy = y;
        }
    }

    int32_t nx = (int32_t)((maxx - minx) / cell_size) + 1;
    int32_t ny = (int32_t)((maxy - miny) / cell_size) + 1;
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    if (nx > 2048) nx = 2048;
    if (ny > 2048) ny = 2048;

    int32_t n_cells = nx * ny;

    /* Count faces per cell */
    int32_t *cell_count = (int32_t *)ARENA_CALLOC(arena, (long)n_cells,
                                                   (long)sizeof(int32_t));
    for (size_t fi = 0; fi < nf; fi++) {
        float fminx = 1e30f, fmaxx = -1e30f;
        float fminy = 1e30f, fmaxy = -1e30f;
        for (int j = 0; j < 3; j++) {
            int32_t vi = faces[fi * 3 + j];
            float x = proj[vi * 2 + 0], y = proj[vi * 2 + 1];
            if (x < fminx) fminx = x; if (x > fmaxx) fmaxx = x;
            if (y < fminy) fminy = y; if (y > fmaxy) fmaxy = y;
        }
        int32_t ix0 = (int32_t)((fminx - minx) / cell_size);
        int32_t ix1 = (int32_t)((fmaxx - minx) / cell_size);
        int32_t iy0 = (int32_t)((fminy - miny) / cell_size);
        int32_t iy1 = (int32_t)((fmaxy - miny) / cell_size);
        if (ix0 < 0) ix0 = 0; if (ix1 >= nx) ix1 = nx - 1;
        if (iy0 < 0) iy0 = 0; if (iy1 >= ny) iy1 = ny - 1;
        for (int32_t iy = iy0; iy <= iy1; iy++)
            for (int32_t ix = ix0; ix <= ix1; ix++)
                cell_count[iy * nx + ix]++;
    }

    /* Prefix sum -> offsets */
    int32_t *cell_off = (int32_t *)ARENA_ALLOC(arena,
        (long)((n_cells + 1) * (long)sizeof(int32_t)));
    cell_off[0] = 0;
    for (int32_t ci = 0; ci < n_cells; ci++)
        cell_off[ci + 1] = cell_off[ci] + cell_count[ci];
    int32_t total_entries = cell_off[n_cells];

    int32_t *cell_faces = (int32_t *)ARENA_ALLOC(arena,
        (long)(total_entries * (long)sizeof(int32_t)));
    int32_t *cell_cursor = (int32_t *)ARENA_CALLOC(arena, (long)n_cells,
                                                    (long)sizeof(int32_t));

    /* Fill grid */
    for (size_t fi = 0; fi < nf; fi++) {
        float fminx = 1e30f, fmaxx = -1e30f;
        float fminy = 1e30f, fmaxy = -1e30f;
        for (int j = 0; j < 3; j++) {
            int32_t vi = faces[fi * 3 + j];
            float x = proj[vi * 2 + 0], y = proj[vi * 2 + 1];
            if (x < fminx) fminx = x; if (x > fmaxx) fmaxx = x;
            if (y < fminy) fminy = y; if (y > fmaxy) fmaxy = y;
        }
        int32_t ix0 = (int32_t)((fminx - minx) / cell_size);
        int32_t ix1 = (int32_t)((fmaxx - minx) / cell_size);
        int32_t iy0 = (int32_t)((fminy - miny) / cell_size);
        int32_t iy1 = (int32_t)((fmaxy - miny) / cell_size);
        if (ix0 < 0) ix0 = 0; if (ix1 >= nx) ix1 = nx - 1;
        if (iy0 < 0) iy0 = 0; if (iy1 >= ny) iy1 = ny - 1;
        for (int32_t iy = iy0; iy <= iy1; iy++) {
            for (int32_t ix = ix0; ix <= ix1; ix++) {
                int32_t ci = iy * nx + ix;
                cell_faces[cell_off[ci] + cell_cursor[ci]] = (int32_t)fi;
                cell_cursor[ci]++;
            }
        }
    }

    /* Compute max potential overlap pairs for allocation */
    size_t max_pairs = 0;
    for (int32_t ci = 0; ci < n_cells; ci++) {
        size_t k = (size_t)(cell_off[ci + 1] - cell_off[ci]);
        max_pairs += k * (k - 1) / 2;
    }
    /* Cap to prevent excessive allocation on degenerate meshes */
    if (max_pairs > nf * 10) max_pairs = nf * 10;
    if (max_pairs == 0) {
        *out_fa = NULL;
        *out_fb = NULL;
        *out_n_ovl = 0;
        return;
    }

    /* Collect overlap pairs (as packed uint64 keys, will dedup) */
    uint64_t *ovl_keys = (uint64_t *)ARENA_ALLOC(arena,
        (long)(max_pairs * sizeof(uint64_t)));
    size_t n_raw = 0;

    for (int32_t ci = 0; ci < n_cells; ci++) {
        int32_t start = cell_off[ci];
        int32_t end = cell_off[ci + 1];
        for (int32_t a = start; a < end; a++) {
            for (int32_t b = a + 1; b < end; b++) {
                int32_t fa = cell_faces[a];
                int32_t fb = cell_faces[b];
                if (fa == fb) continue;
                if (faces_share_vertex(faces, fa, fb)) continue;
                if (triangles_overlap_2d(proj, faces, fa, fb)) {
                    int32_t lo = fa < fb ? fa : fb;
                    int32_t hi = fa < fb ? fb : fa;
                    if (n_raw < max_pairs) {
                        ovl_keys[n_raw++] =
                            ((uint64_t)(uint32_t)lo << 32) |
                            (uint64_t)(uint32_t)hi;
                    }
                }
            }
        }
    }

    if (n_raw == 0) {
        *out_fa = NULL;
        *out_fb = NULL;
        *out_n_ovl = 0;
        return;
    }

    /* Sort and dedup */
    qsort(ovl_keys, n_raw, sizeof(uint64_t), compare_u64);
    size_t n_unique = 0;
    for (size_t idx = 0; idx < n_raw; idx++) {
        if (idx == 0 || ovl_keys[idx] != ovl_keys[idx - 1]) {
            ovl_keys[n_unique++] = ovl_keys[idx];
        }
    }

    /* Extract face pairs */
    int32_t *ofa = (int32_t *)ARENA_ALLOC(arena,
        (long)(n_unique * sizeof(int32_t)));
    int32_t *ofb = (int32_t *)ARENA_ALLOC(arena,
        (long)(n_unique * sizeof(int32_t)));
    for (size_t idx = 0; idx < n_unique; idx++) {
        ofa[idx] = (int32_t)(ovl_keys[idx] >> 32);
        ofb[idx] = (int32_t)(ovl_keys[idx] & 0xFFFFFFFF);
    }

    *out_fa = ofa;
    *out_fb = ofb;
    *out_n_ovl = n_unique;
}

/* ===== Delamination-merge audit/fitting helpers ===== */

typedef struct {
    float u_axis[3];
    float v_axis[3];
    const int32_t *original_ovl_fa;
    const int32_t *original_ovl_fb;
    size_t n_original_ovl;
    const int32_t *accepted_exterior_source_faces;
    size_t accepted_patch_face_start;
} MergeAuditFrame;

typedef struct {
    size_t total;
    size_t allowed_exterior;
    size_t unexpected_exterior;
    size_t patch;
    size_t patch_patch;
    size_t patch_exterior;
} MergeOverlapAuditStats;

static int merge_original_overlap_allowed(const MergeAuditFrame *frame,
                                          int32_t fa, int32_t fb)
{
    int32_t lo = fa < fb ? fa : fb;
    int32_t hi = fa < fb ? fb : fa;
    uint64_t key = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
    size_t left = 0, right = frame->n_original_ovl;
    while (left < right) {
        size_t mid = left + (right - left) / 2;
        int32_t oa = frame->original_ovl_fa[mid];
        int32_t ob = frame->original_ovl_fb[mid];
        int32_t olo = oa < ob ? oa : ob;
        int32_t ohi = oa < ob ? ob : oa;
        uint64_t other = ((uint64_t)(uint32_t)olo << 32) | (uint32_t)ohi;
        if (other < key) left = mid + 1;
        else right = mid;
    }
    if (left >= frame->n_original_ovl) return 0;
    int32_t oa = frame->original_ovl_fa[left];
    int32_t ob = frame->original_ovl_fb[left];
    int32_t olo = oa < ob ? oa : ob;
    int32_t ohi = oa < ob ? ob : oa;
    uint64_t found = ((uint64_t)(uint32_t)olo << 32) | (uint32_t)ohi;
    return found == key;
}

static int merge_candidate_overlap_audit(
    Arena_T arena, const ComponentMesh *candidate,
    size_t patch_face_start, const int32_t *exterior_source_faces,
    MergeAuditFrame *frame, int emit_debug, MergeOverlapAuditStats *out_stats)
{
    Arena_Mark mark = Arena_save(arena);
    float *candidate_proj = (float *)ARENA_ALLOC(
        arena, candidate->nv * 2 * sizeof(float));
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        const float *p = &candidate->verts[vi * 3];
        candidate_proj[vi * 2] =
            p[0] * frame->u_axis[0] +
            p[1] * frame->u_axis[1] +
            p[2] * frame->u_axis[2];
        candidate_proj[vi * 2 + 1] =
            p[0] * frame->v_axis[0] +
            p[1] * frame->v_axis[1] +
            p[2] * frame->v_axis[2];
    }

    int32_t *fa = NULL, *fb = NULL;
    size_t n_overlap = 0;
    detect_overlaps(arena, candidate->faces, candidate->nf,
                    candidate_proj, candidate->nv,
                    &fa, &fb, &n_overlap);
    size_t allowed_exterior = 0;
    size_t unexpected_exterior = 0;
    size_t patch_overlap = 0;
    size_t patch_patch = 0;
    size_t patch_exterior = 0;
    uint8_t *overlap_face = emit_debug && g_debug_dir
        ? (uint8_t *)ARENA_CALLOC(
              arena, candidate->nf, sizeof(*overlap_face))
        : NULL;
    for (size_t i = 0; i < n_overlap; i++) {
        if ((size_t)fa[i] >= patch_face_start ||
            (size_t)fb[i] >= patch_face_start) {
            patch_overlap++;
            if ((size_t)fa[i] >= patch_face_start &&
                (size_t)fb[i] >= patch_face_start) patch_patch++;
            else patch_exterior++;
            if (overlap_face) {
                overlap_face[fa[i]] = 1;
                overlap_face[fb[i]] = 1;
            }
        } else if (merge_original_overlap_allowed(
                       frame, exterior_source_faces[fa[i]],
                       exterior_source_faces[fb[i]])) {
            allowed_exterior++;
        } else {
            unexpected_exterior++;
        }
    }
    if (out_stats) {
        out_stats->total = n_overlap;
        out_stats->allowed_exterior = allowed_exterior;
        out_stats->unexpected_exterior = unexpected_exterior;
        out_stats->patch = patch_overlap;
        out_stats->patch_patch = patch_patch;
        out_stats->patch_exterior = patch_exterior;
    }
    if (emit_debug) {
        fprintf(stderr,
                "      merge overlap audit: total=%zu allowed=%zu "
                "unexpected-exterior=%zu patch=%zu (patch-patch=%zu, "
                "patch-exterior=%zu)\n",
                n_overlap, allowed_exterior, unexpected_exterior,
                patch_overlap, patch_patch, patch_exterior);
    }
    int clean = unexpected_exterior == 0 && patch_overlap == 0;
    if (emit_debug && !clean && g_debug_dir) {
        char path[1024];
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_rejected_patchstart_%zu.obj",
                 g_debug_dir, patch_face_start);
        ObjIO_write(path, candidate->verts, candidate->nv,
                    candidate->faces, candidate->nf);
        fprintf(stderr, "      debug: wrote %s\n", path);
        snprintf(path, sizeof(path),
                 "%s/07_fitted_residual_overlap_faces.obj", g_debug_dir);
        debug_write_face_subset(
            arena, path, candidate->verts, candidate->nv,
            candidate->faces, candidate->nf, overlap_face, NULL, NULL,
            NULL);
    }
    if (clean) {
        frame->accepted_exterior_source_faces = exterior_source_faces;
        frame->accepted_patch_face_start = patch_face_start;
    }
    Arena_restore(arena, mark);
    return clean;
}

static int merge_candidate_overlap_free(Arena_T arena,
                                        const ComponentMesh *candidate,
                                        size_t patch_face_start,
                                        const int32_t *exterior_source_faces,
                                        void *user)
{
    return merge_candidate_overlap_audit(
        arena, candidate, patch_face_start, exterior_source_faces,
        (MergeAuditFrame *)user, 1, NULL);
}

static size_t merge_candidate_patch_overlap_vertices(
    Arena_T arena, const ComponentMesh *candidate, size_t patch_face_start,
    const MergeAuditFrame *frame, uint8_t *vertex_mask)
{
    Arena_Mark mark = Arena_save(arena);
    memset(vertex_mask, 0, candidate->nv * sizeof(*vertex_mask));
    float *proj = (float *)ARENA_ALLOC(
        arena, candidate->nv * 2 * sizeof(*proj));
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        const float *p = &candidate->verts[vi * 3];
        proj[vi * 2] =
            p[0] * frame->u_axis[0] + p[1] * frame->u_axis[1] +
            p[2] * frame->u_axis[2];
        proj[vi * 2 + 1] =
            p[0] * frame->v_axis[0] + p[1] * frame->v_axis[1] +
            p[2] * frame->v_axis[2];
    }
    int32_t *fa = NULL, *fb = NULL;
    size_t n_overlap = 0, patch_overlap = 0;
    detect_overlaps(arena, candidate->faces, candidate->nf, proj,
                    candidate->nv, &fa, &fb, &n_overlap);
    for (size_t i = 0; i < n_overlap; i++) {
        int patch_a = (size_t)fa[i] >= patch_face_start;
        int patch_b = (size_t)fb[i] >= patch_face_start;
        if (!patch_a && !patch_b) continue;
        patch_overlap++;
        if (patch_a) {
            const int32_t *tri = &candidate->faces[(size_t)fa[i] * 3];
            for (int k = 0; k < 3; k++) vertex_mask[tri[k]] = 1;
        }
        if (patch_b) {
            const int32_t *tri = &candidate->faces[(size_t)fb[i] * 3];
            for (int k = 0; k < 3; k++) vertex_mask[tri[k]] = 1;
        }
    }
    Arena_restore(arena, mark);
    return patch_overlap;
}

static int pred_is_surface(const uint8_t *vol,
                           size_t D, size_t H, size_t W,
                           size_t z, size_t y, size_t x)
{
    size_t hw = H * W;
    size_t index = z * hw + y * W + x;
    if (vol[index] == 0) return 0;
    if (z == 0 || vol[index - hw] == 0) return 1;
    if (z + 1 >= D || vol[index + hw] == 0) return 1;
    if (y == 0 || vol[index - W] == 0) return 1;
    if (y + 1 >= H || vol[index + W] == 0) return 1;
    if (x == 0 || vol[index - 1] == 0) return 1;
    if (x + 1 >= W || vol[index + 1] == 0) return 1;
    return 0;
}

static int merge_prediction_fit(Arena_T arena,
                                 const ComponentMesh *candidate,
                                 const uint8_t *vertex_role,
                                 const uint8_t *vol,
                                 size_t D, size_t H, size_t W,
                                 const float coordinate_offset[3],
                                 float distance,
                                 double *out_matched_fraction,
                                 double *out_mean,
                                 double *out_p95)
{
    size_t n_free = 0;
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        if (!vertex_role || vertex_role[vi] ==
            OVERLAP_MERGE_VERTEX_REPAIR) n_free++;
    }
    if (n_free == 0 || !vol || D == 0 || H == 0 || W == 0) return -1;

    float *distances = (float *)ARENA_ALLOC(
        arena, n_free * sizeof(*distances));
    int radius = (int)ceil((double)distance);
    double max_d2 = (double)distance * distance;
    size_t matched = 0;
    double sum = 0.0;
    size_t cursor = 0;

    for (size_t vi = 0; vi < candidate->nv; vi++) {
        if (vertex_role && vertex_role[vi] !=
            OVERLAP_MERGE_VERTEX_REPAIR) continue;
        const float *p = &candidate->verts[vi * 3];
        double pz = p[0] + (coordinate_offset ? coordinate_offset[0] : 0.0);
        double py = p[1] + (coordinate_offset ? coordinate_offset[1] : 0.0);
        double px = p[2] + (coordinate_offset ? coordinate_offset[2] : 0.0);
        int cz = (int)floor(pz + 0.5);
        int cy = (int)floor(py + 0.5);
        int cx = (int)floor(px + 0.5);
        double best = DBL_MAX;
        for (int dz = -radius; dz <= radius; dz++) {
            int z = cz + dz;
            if (z < 0 || (size_t)z >= D) continue;
            for (int dy = -radius; dy <= radius; dy++) {
                int y = cy + dy;
                if (y < 0 || (size_t)y >= H) continue;
                for (int dx = -radius; dx <= radius; dx++) {
                    int x = cx + dx;
                    if (x < 0 || (size_t)x >= W) continue;
                    if (!pred_is_surface(vol, D, H, W,
                                         (size_t)z, (size_t)y, (size_t)x))
                        continue;
                    double ez = pz - z;
                    double ey = py - y;
                    double ex = px - x;
                    double d2 = ez * ez + ey * ey + ex * ex;
                    if (d2 < best) best = d2;
                }
            }
        }
        if (best <= max_d2) {
            distances[cursor] = (float)sqrt(best);
            matched++;
            sum += distances[cursor];
        } else {
            distances[cursor] = distance + 1.0f;
        }
        cursor++;
    }

    qsort(distances, n_free, sizeof(*distances), compare_float);
    size_t p95_index = (size_t)floor(0.95 * (double)(n_free - 1));
    *out_matched_fraction = (double)matched / (double)n_free;
    *out_mean = (matched > 0) ? sum / (double)matched : DBL_MAX;
    *out_p95 = distances[p95_index];
    return 0;
}

static void merge_translate_for_prediction(ComponentMesh *candidate,
                                            const float offset[3],
                                            float direction)
{
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        candidate->verts[vi * 3] += direction * offset[0];
        candidate->verts[vi * 3 + 1] += direction * offset[1];
        candidate->verts[vi * 3 + 2] += direction * offset[2];
    }
}

typedef struct MergeRawFitStats {
    size_t n_free;
    size_t n_targeted;
    size_t n_supported;
    size_t n_moved;
    size_t n_reverted;
    double mean_raw_before;
    double mean_raw_after;
    double mean_raw_target;
    double mean_data_cost;
    double quilt_weight;
    double support_fraction;
    double line_search_scale;
    long long mrf_energy;
    double mean_displacement;
    double max_displacement;
    size_t harmonic_patch_overlaps;
    size_t clamped_patch_overlaps;
    size_t clamped_patch_patch;
    size_t clamped_patch_exterior;
    size_t fitted_patch_overlaps;
    size_t n_clamped;
    double clamp_scale;
    int held_clamp_boundary;
    size_t n_injective_adjusted;
    double max_injective_adjustment;
    int local_relax_rings;
    double local_relax_scale;
    size_t raw_full_patch_overlaps;
    size_t raw_local_adjusted;
    int raw_local_relax_rings;
    double raw_local_relax_scale;
    double raw_local_max_adjustment;
    size_t source_interface_pairs;
    size_t fitted_interface_pairs;
    double source_interface_support;
    double fitted_interface_support;
    double source_quilt_mean;
    double fitted_quilt_mean;
    double source_quilt_p95;
    double fitted_quilt_p95;
    double source_bend_p95;
    double fitted_bend_p95;
    double source_best_material;
    double baseline_material;
    double fitted_material;
    double source_material_support;
    double baseline_material_support;
    double fitted_material_support;
    int evidence_valid;
} MergeRawFitStats;

/* Fit only the newly triangulated repair vertices to the ORIGINAL raw CT.
 * `pin_mask` is deliberately repurposed by OverlapMerge as a Dirichlet mask:
 * every retained source vertex is pinned and only new patch vertices may move.
 *
 * SnapDetect searches raw intensity along the across-sheet direction and uses
 * the existing mesh as an occupancy guard.  SnapSolve then performs the
 * data-weighted Laplacian fit.  The prediction volume is not visible here. */
static int merge_raw_volume_fit(Arena_T arena,
                                ComponentMesh *candidate,
                                const OverlapSepOptions *options,
                                const float plane_normal[3],
                                MergeRawFitStats *stats)
{
    memset(stats, 0, sizeof(*stats));
    if (!candidate || !candidate->verts || !candidate->faces ||
        !candidate->pin_mask || !options || !options->raw_dir ||
        !*options->raw_dir) return -1;

    for (size_t vi = 0; vi < candidate->nv; vi++)
        if (!candidate->pin_mask[vi]) stats->n_free++;
    if (stats->n_free == 0) return -1;

    float axis_u[3], axis_v[3];
    PCA_orthonormal_basis(plane_normal, axis_u, axis_v);
    float *uv = (float *)ARENA_ALLOC(
        arena, candidate->nv * 2 * sizeof(*uv));
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        const float *p = &candidate->verts[vi * 3];
        uv[vi * 2] = p[0] * axis_u[0] + p[1] * axis_u[1] +
                     p[2] * axis_u[2];
        uv[vi * 2 + 1] = p[0] * axis_v[0] + p[1] * axis_v[1] +
                         p[2] * axis_v[2];
    }

    /* Raw cube lookup is in scroll/world coordinates; overlap repair runs in
     * cube-local coordinates. */
    merge_translate_for_prediction(candidate, options->world_offset, 1.0f);

    SnapOpts snap;
    SnapOpts_default(&snap);
    snap.raw_dir = options->raw_dir;
    snap.chunk = options->raw_chunk > 0 ? options->raw_chunk : 128;
    snap.reach = options->raw_snap_reach > 0.0f
        ? options->raw_snap_reach : 8.0;
    snap.min_region = 1;
    snap.region_cap = INT_MAX;
    snap.anchor_frac = 0.0;
    snap.dilate_rings = 1;
    snap.target_mode = SNAP_TARGET_QUILT_MRF;
    snap.depth_bins = 33;
    snap.w_match = 30.0;
    snap.w_close = 200.0;
    snap.smooth_mu = 2000.0;
    snap.smooth_tau = 4.0;
    snap.verbose = 1;

    SnapResult before;
    if (SnapDetect_run(arena, candidate->verts, candidate->nv,
                       candidate->faces, candidate->nf, uv,
                       &snap, &before) != 0) {
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }

    size_t *region_targets = (size_t *)ARENA_CALLOC(
        arena, before.nreg > 0 ? before.nreg : 1,
        sizeof(*region_targets));
    double raw_sum_before = 0.0;
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        if (candidate->pin_mask[vi]) continue;
        raw_sum_before += before.vcur[vi];
        int32_t r = before.vregion[vi];
        if (r >= 0 && (size_t)r < before.nreg) {
            if (before.vclass[vi] == SNAP_FIXABLE) region_targets[r]++;
        }
    }
    stats->mean_raw_before = raw_sum_before / (double)stats->n_free;

    /* Region classification was computed over the complete candidate.  For
     * this local repair, a region with a raw target on the movable patch is a
     * valid island even when unrelated pinned exterior vertices make the
     * global region look anchorless.  Pinned vertices themselves are forced
     * GOOD, so SnapSolve cannot move them. */
    for (size_t r = 0; r < before.nreg; r++) {
        if (region_targets[r] > 0) before.rclass[r] = SNAPREG_FIXABLE;
    }
    for (size_t vi = 0; vi < candidate->nv; vi++) {
        if (candidate->pin_mask[vi]) {
            before.vclass[vi] = SNAP_GOOD;
        } else if (before.vclass[vi] == SNAP_FIXABLE) {
            int32_t r = before.vregion[vi];
            if (r >= 0 && (size_t)r < before.nreg &&
                before.rclass[r] == SNAPREG_FIXABLE)
                stats->n_targeted++;
        }
    }

    if (g_debug_dir) {
        char path[1024];
        float *target = (float *)ARENA_ALLOC(
            arena, candidate->nv * 3 * sizeof(*target));
        float *colors = (float *)ARENA_ALLOC(
            arena, candidate->nv * 3 * sizeof(*colors));
        for (size_t vi = 0; vi < candidate->nv; vi++) {
            for (int d = 0; d < 3; d++) {
                target[vi * 3 + (size_t)d] =
                    candidate->verts[vi * 3 + (size_t)d] -
                    options->world_offset[d];
                colors[vi * 3 + (size_t)d] = 0.45f;
            }
            if (candidate->pin_mask[vi]) {
                colors[vi * 3] = 0.20f;
                colors[vi * 3 + 1] = 0.55f;
                colors[vi * 3 + 2] = 0.95f;
            } else if (before.vclass[vi] == SNAP_FIXABLE) {
                for (int d = 0; d < 3; d++)
                    target[vi * 3 + (size_t)d] +=
                        before.voff[vi] * before.vdir[vi * 3 + (size_t)d];
                colors[vi * 3] = 0.95f;
                colors[vi * 3 + 1] = 0.20f;
                colors[vi * 3 + 2] = 0.85f;
            } else if (before.vclass[vi] == SNAP_CRACK) {
                colors[vi * 3] = 0.95f;
                colors[vi * 3 + 1] = 0.20f;
                colors[vi * 3 + 2] = 0.15f;
            }
        }
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_raw_targets.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(path, target, candidate->nv,
                                     candidate->faces, candidate->nf,
                                     colors);
    }

    SnapSolveOpts solve;
    SnapSolveOpts_default(&solve);
    if (options->raw_snap_alpha > 0.0f &&
        options->raw_snap_alpha < 1.0f)
        solve.alpha = options->raw_snap_alpha;
    solve.max_disp = snap.reach;
    solve.verbose = 1;
    SnapSolveResult solved;
    int solve_rc = SnapSolve_run(
        arena, candidate->verts, candidate->nv,
        candidate->faces, candidate->nf, &before, &solve, &solved);
    stats->n_moved = solved.n_moved;
    stats->n_reverted = solved.n_reverted;
    stats->mean_displacement = solved.mean_disp;
    stats->max_displacement = solved.max_disp;

    SnapResult after;
    memset(&after, 0, sizeof(after));
    int after_rc = solve_rc == 0 ? SnapDetect_run(
        arena, candidate->verts, candidate->nv,
        candidate->faces, candidate->nf, uv, &snap, &after) : -1;
    if (after_rc == 0) {
        double raw_sum_after = 0.0;
        for (size_t vi = 0; vi < candidate->nv; vi++)
            if (!candidate->pin_mask[vi]) raw_sum_after += after.vcur[vi];
        stats->mean_raw_after = raw_sum_after / (double)stats->n_free;
    }

    merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
    return solve_rc == 0 && after_rc == 0 ? 0 : -1;
}

static int32_t rawfit_cost(double value)
{
    if (value < 0.0) value = 0.0;
    if (value > 9500.0) value = 9500.0;
    return (int32_t)lrint(value);
}

typedef struct RawfitWeightedValue {
    double value;
    double weight;
} RawfitWeightedValue;

typedef struct MergeRawInterfaceStats {
    size_t pairs;
    size_t supported_pairs;
    double support_fraction;
    double brightness_offset;
    double quilt_mean;
    double quilt_p95;
    double bend_p95_degrees;
} MergeRawInterfaceStats;

static int compare_rawfit_weighted_value(const void *a, const void *b)
{
    const RawfitWeightedValue *va = (const RawfitWeightedValue *)a;
    const RawfitWeightedValue *vb = (const RawfitWeightedValue *)b;
    if (va->value < vb->value) return -1;
    if (va->value > vb->value) return 1;
    return 0;
}

static double rawfit_weighted_percentile(RawfitWeightedValue *values,
                                         size_t count, double percentile)
{
    if (count == 0) return DBL_MAX;
    qsort(values, count, sizeof(*values), compare_rawfit_weighted_value);
    double total = 0.0;
    for (size_t i = 0; i < count; i++) total += values[i].weight;
    if (total <= 0.0) return values[count - 1].value;
    double target = percentile * total;
    double accum = 0.0;
    for (size_t i = 0; i < count; i++) {
        accum += values[i].weight;
        if (accum >= target) return values[i].value;
    }
    return values[count - 1].value;
}

static int rawfit_face_frame(const ComponentMesh *mesh, size_t face,
                             const float offset[3], float centroid[3],
                             float normal[3], double *area)
{
    const int32_t *tri = &mesh->faces[face * 3];
    const float *a = &mesh->verts[(size_t)tri[0] * 3];
    const float *b = &mesh->verts[(size_t)tri[1] * 3];
    const float *c = &mesh->verts[(size_t)tri[2] * 3];
    double ab[3], ac[3], cross[3];
    for (int d = 0; d < 3; d++) {
        ab[d] = (double)b[d] - a[d];
        ac[d] = (double)c[d] - a[d];
        centroid[d] = (a[d] + b[d] + c[d]) / 3.0f + offset[d];
    }
    cross[0] = ab[1] * ac[2] - ab[2] * ac[1];
    cross[1] = ab[2] * ac[0] - ab[0] * ac[2];
    cross[2] = ab[0] * ac[1] - ab[1] * ac[0];
    double length = sqrt(cross[0] * cross[0] + cross[1] * cross[1] +
                         cross[2] * cross[2]);
    if (length <= 1.0e-12) return -1;
    for (int d = 0; d < 3; d++) normal[d] = (float)(cross[d] / length);
    *area = 0.5 * length;
    return 0;
}

static double rawfit_shared_edge_length(const ComponentMesh *mesh,
                                        size_t face_a, size_t face_b)
{
    const int32_t *a = &mesh->faces[face_a * 3];
    const int32_t *b = &mesh->faces[face_b * 3];
    int32_t shared[2] = {-1, -1};
    int count = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (a[i] == b[j] && count < 2) shared[count++] = a[i];
        }
    }
    if (count != 2) return 0.0;
    const float *p = &mesh->verts[(size_t)shared[0] * 3];
    const float *q = &mesh->verts[(size_t)shared[1] * 3];
    double d0 = (double)p[0] - q[0];
    double d1 = (double)p[1] - q[1];
    double d2 = (double)p[2] - q[2];
    return sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

static int merge_raw_interface_metrics(
    Arena_T arena, CubeTable *ct, const ComponentMesh *mesh,
    const int32_t *adj_fa, const int32_t *adj_fb, size_t n_adj,
    const uint8_t *patch_face, const float offset[3],
    MergeRawInterfaceStats *out, double *face_energy)
{
    memset(out, 0, sizeof(*out));
    if (face_energy)
        for (size_t fi = 0; fi < mesh->nf; fi++) face_energy[fi] = NAN;

    double *difference = (double *)ARENA_ALLOC(
        arena, n_adj * sizeof(*difference));
    double *angle = (double *)ARENA_ALLOC(arena, n_adj * sizeof(*angle));
    double *weight = (double *)ARENA_ALLOC(arena, n_adj * sizeof(*weight));
    int32_t *pair_patch = (int32_t *)ARENA_ALLOC(
        arena, n_adj * sizeof(*pair_patch));
    int32_t *pair_anchor = (int32_t *)ARENA_ALLOC(
        arena, n_adj * sizeof(*pair_anchor));
    size_t count = 0;
    for (size_t i = 0; i < n_adj; i++) {
        int32_t fa = adj_fa[i], fb = adj_fb[i];
        if (patch_face[fa] == patch_face[fb]) continue;
        out->pairs++;
        int32_t fp = patch_face[fa] ? fa : fb;
        int32_t fe = patch_face[fa] ? fb : fa;
        float cp[3], ce[3], np[3], ne[3];
        double ap = 0.0, ae = 0.0;
        if (rawfit_face_frame(mesh, (size_t)fp, offset, cp, np, &ap) != 0 ||
            rawfit_face_frame(mesh, (size_t)fe, offset, ce, ne, &ae) != 0)
            continue;
        double vp = sample_vertex(ct, cp, np, 2.0, 5);
        double ve = sample_vertex(ct, ce, ne, 2.0, 5);
        if (vp < 0.0 || ve < 0.0) continue;
        double dot = (double)np[0] * ne[0] + (double)np[1] * ne[1] +
                     (double)np[2] * ne[2];
        if (dot < -1.0) dot = -1.0;
        if (dot > 1.0) dot = 1.0;
        double edge_weight = rawfit_shared_edge_length(
            mesh, (size_t)fp, (size_t)fe);
        if (edge_weight <= 1.0e-8) edge_weight = 0.5 * (ap + ae);
        difference[count] = vp - ve;
        angle[count] = acos(dot) * (180.0 / 3.14159265358979323846);
        weight[count] = edge_weight;
        pair_patch[count] = fp;
        pair_anchor[count] = fe;
        count++;
    }
    out->supported_pairs = count;
    out->support_fraction = out->pairs > 0
        ? (double)count / (double)out->pairs : 0.0;
    if (count == 0) return -1;

    double weight_sum = 0.0, offset_sum = 0.0;
    for (size_t i = 0; i < count; i++) {
        weight_sum += weight[i];
        offset_sum += weight[i] * difference[i];
    }
    if (weight_sum <= 0.0) return -1;
    out->brightness_offset = offset_sum / weight_sum;
    RawfitWeightedValue *quilt = (RawfitWeightedValue *)ARENA_ALLOC(
        arena, count * sizeof(*quilt));
    RawfitWeightedValue *bend = (RawfitWeightedValue *)ARENA_ALLOC(
        arena, count * sizeof(*bend));
    double error_sum = 0.0;
    double *face_sum = face_energy ? (double *)ARENA_CALLOC(
        arena, mesh->nf, sizeof(*face_sum)) : NULL;
    double *face_weight = face_energy ? (double *)ARENA_CALLOC(
        arena, mesh->nf, sizeof(*face_weight)) : NULL;
    for (size_t i = 0; i < count; i++) {
        double error = fabs(difference[i] - out->brightness_offset);
        quilt[i].value = error;
        quilt[i].weight = weight[i];
        bend[i].value = angle[i];
        bend[i].weight = weight[i];
        error_sum += error * weight[i];
        if (face_energy) {
            face_sum[pair_patch[i]] += error * weight[i];
            face_sum[pair_anchor[i]] += error * weight[i];
            face_weight[pair_patch[i]] += weight[i];
            face_weight[pair_anchor[i]] += weight[i];
        }
    }
    if (face_energy) {
        for (size_t fi = 0; fi < mesh->nf; fi++)
            if (face_weight[fi] > 0.0)
                face_energy[fi] = face_sum[fi] / face_weight[fi];
    }
    out->quilt_mean = error_sum / weight_sum;
    out->quilt_p95 = rawfit_weighted_percentile(quilt, count, 0.95);
    out->bend_p95_degrees = rawfit_weighted_percentile(bend, count, 0.95);
    return 0;
}

static int merge_raw_patch_material(
    Arena_T arena, CubeTable *ct, const ComponentMesh *source,
    const uint8_t *source_patch, const int32_t *face_labels,
    int32_t num_labels, const float source_offset[3],
    const ComponentMesh *candidate, size_t patch_face_start,
    double *out_source_best, double *out_candidate,
    double *out_source_support, double *out_candidate_support)
{
    double *label_sum = (double *)ARENA_CALLOC(
        arena, (size_t)num_labels, sizeof(*label_sum));
    double *label_weight = (double *)ARENA_CALLOC(
        arena, (size_t)num_labels, sizeof(*label_weight));
    double source_total = 0.0, source_valid = 0.0;
    for (size_t fi = 0; fi < source->nf; fi++) {
        if (!source_patch[fi]) continue;
        float c[3], n[3];
        double area = 0.0;
        if (rawfit_face_frame(source, fi, source_offset, c, n, &area) != 0)
            continue;
        source_total += area;
        double value = sample_vertex(ct, c, n, 2.0, 5);
        if (value < 0.0) continue;
        int32_t label = face_labels[fi];
        if (label < 0 || label >= num_labels) continue;
        label_sum[label] += value * area;
        label_weight[label] += area;
        source_valid += area;
    }
    *out_source_best = -1.0;
    for (int32_t label = 0; label < num_labels; label++) {
        if (label_weight[label] <= 0.0) continue;
        double mean = label_sum[label] / label_weight[label];
        if (mean > *out_source_best) *out_source_best = mean;
    }
    *out_source_support = source_total > 0.0
        ? source_valid / source_total : 0.0;

    double candidate_sum = 0.0, candidate_total = 0.0;
    double candidate_valid = 0.0;
    const float zero_offset[3] = {0.0f, 0.0f, 0.0f};
    for (size_t fi = patch_face_start; fi < candidate->nf; fi++) {
        float c[3], n[3];
        double area = 0.0;
        if (rawfit_face_frame(candidate, fi, zero_offset, c, n, &area) != 0)
            continue;
        candidate_total += area;
        double value = sample_vertex(ct, c, n, 2.0, 5);
        if (value < 0.0) continue;
        candidate_sum += value * area;
        candidate_valid += area;
    }
    *out_candidate = candidate_valid > 0.0
        ? candidate_sum / candidate_valid : -1.0;
    *out_candidate_support = candidate_total > 0.0
        ? candidate_valid / candidate_total : 0.0;
    return *out_source_best >= 0.0 && *out_candidate >= 0.0 ? 0 : -1;
}

static int rawfit_patch_no_fold(const float *verts,
                                const int32_t *faces, size_t nf,
                                size_t patch_start,
                                const double normal[3], int reference_sign)
{
    for (size_t fi = patch_start; fi < nf; fi++) {
        const int32_t *tri = &faces[fi * 3];
        const float *a = &verts[(size_t)tri[0] * 3];
        const float *b = &verts[(size_t)tri[1] * 3];
        const float *c = &verts[(size_t)tri[2] * 3];
        double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        double dot =
            (ab[1] * ac[2] - ab[2] * ac[1]) * normal[0] +
            (ab[2] * ac[0] - ab[0] * ac[2]) * normal[1] +
            (ab[0] * ac[1] - ab[1] * ac[0]) * normal[2];
        if ((double)reference_sign * dot <= 1.0e-10) return 0;
    }
    return 1;
}

static int rawfit_patch_reference_sign(const float *verts,
                                       const int32_t *faces, size_t nf,
                                       size_t patch_start,
                                       const double normal[3])
{
    double sum = 0.0;
    for (size_t fi = patch_start; fi < nf; fi++) {
        const int32_t *tri = &faces[fi * 3];
        const float *a = &verts[(size_t)tri[0] * 3];
        const float *b = &verts[(size_t)tri[1] * 3];
        const float *c = &verts[(size_t)tri[2] * 3];
        double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        sum += (ab[1] * ac[2] - ab[2] * ac[1]) * normal[0] +
               (ab[2] * ac[0] - ab[0] * ac[2]) * normal[1] +
               (ab[0] * ac[1] - ab[1] * ac[0]) * normal[2];
    }
    return sum >= 0.0 ? 1 : -1;
}

static void rawfit_apply_height_blend(float *verts,
                                      const int32_t *vertex_of,
                                      size_t n_free,
                                      const float *repair_base,
                                      const float *dirichlet_height,
                                      const float *clamped_height,
                                      const uint8_t *hold_full,
                                      const double normal[3], double scale)
{
    for (size_t site = 0; site < n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        double h = hold_full && hold_full[vi]
            ? clamped_height[vi]
            : (double)dirichlet_height[vi] + scale *
              ((double)clamped_height[vi] - dirichlet_height[vi]);
        for (int d = 0; d < 3; d++)
            verts[vi * 3 + (size_t)d] =
                repair_base[site * 3 + (size_t)d] +
                (float)(h * normal[d]);
    }
}

static void rawfit_apply_local_height_blend(
    float *verts, const int32_t *vertex_of, size_t n_free,
    const float *repair_base, const float *safe_height,
    const float *target_height, const uint8_t *relax_mask,
    const double normal[3], double scale)
{
    for (size_t site = 0; site < n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        double h = target_height[vi];
        if (relax_mask[vi])
            h = (double)safe_height[vi] + scale *
                ((double)target_height[vi] - safe_height[vi]);
        for (int d = 0; d < 3; d++)
            verts[vi * 3 + (size_t)d] =
                repair_base[site * 3 + (size_t)d] +
                (float)(h * normal[d]);
    }
}

static void rawfit_apply_local_position_blend(
    float *verts, const int32_t *vertex_of, size_t n_free,
    const float *safe_position, const float *target_position,
    const uint8_t *relax_mask, double scale)
{
    for (size_t site = 0; site < n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        const float *safe = &safe_position[vi * 3];
        const float *target = &target_position[vi * 3];
        for (int d = 0; d < 3; d++) {
            verts[vi * 3 + (size_t)d] = relax_mask[vi]
                ? safe[d] + (float)(scale * (target[d] - safe[d]))
                : target[d];
        }
    }
}

static double rawfit_projected_twice_area(const float *verts,
                                          const int32_t tri[3],
                                          const MergeAuditFrame *frame)
{
    double p[3][2];
    for (int k = 0; k < 3; k++) {
        const float *v = &verts[(size_t)tri[k] * 3];
        p[k][0] = (double)v[0] * frame->u_axis[0] +
                  (double)v[1] * frame->u_axis[1] +
                  (double)v[2] * frame->u_axis[2];
        p[k][1] = (double)v[0] * frame->v_axis[0] +
                  (double)v[1] * frame->v_axis[1] +
                  (double)v[2] * frame->v_axis[2];
    }
    return (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) -
           (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]);
}

/* Project the fair height field onto the production projection's orientation
 * half-spaces.  The common chart and exterior boundary are fixed.  Moving a
 * repair vertex along the common-plane normal changes a projected triangle's
 * signed area linearly, so cyclic orthogonal projections give a local,
 * deterministic untangle without relaxing the collar tangent. */
static int rawfit_enforce_projected_orientation(
    ComponentMesh *candidate, size_t patch_face_start,
    const double *reference_area, const int32_t *site_of,
    const int32_t *vertex_of, size_t n_free, const float *repair_base,
    float *height, const uint8_t *fixed, const double normal[3],
    const MergeAuditFrame *frame, size_t *out_adjusted,
    double *out_max_adjustment)
{
    const size_t patch_nf = candidate->nf - patch_face_start;
    const double qx = normal[0] * frame->u_axis[0] +
                      normal[1] * frame->u_axis[1] +
                      normal[2] * frame->u_axis[2];
    const double qy = normal[0] * frame->v_axis[0] +
                      normal[1] * frame->v_axis[1] +
                      normal[2] * frame->v_axis[2];
    if (qx * qx + qy * qy <= 1.0e-16) return -1;

    float *initial = (float *)malloc(n_free * sizeof(*initial));
    if (!initial) return -1;
    for (size_t site = 0; site < n_free; site++)
        initial[site] = height[(size_t)vertex_of[site]];

    int feasible = 0;
    for (int iteration = 0; iteration < 4000; iteration++) {
        size_t violated = 0;
        double worst = 0.0;
        for (size_t local_fi = 0; local_fi < patch_nf; local_fi++) {
            size_t fi = patch_face_start + local_fi;
            const int32_t *tri = &candidate->faces[fi * 3];
            double p[3][2];
            for (int k = 0; k < 3; k++) {
                const float *v = &candidate->verts[(size_t)tri[k] * 3];
                p[k][0] = (double)v[0] * frame->u_axis[0] +
                          (double)v[1] * frame->u_axis[1] +
                          (double)v[2] * frame->u_axis[2];
                p[k][1] = (double)v[0] * frame->v_axis[0] +
                          (double)v[1] * frame->v_axis[1] +
                          (double)v[2] * frame->v_axis[2];
            }
            double area = (p[1][0] - p[0][0]) *
                          (p[2][1] - p[0][1]) -
                          (p[1][1] - p[0][1]) *
                          (p[2][0] - p[0][0]);
            double sign = reference_area[local_fi] >= 0.0 ? 1.0 : -1.0;
            double target = 0.05 * fabs(reference_area[local_fi]);
            if (target < 1.0e-8) target = 1.0e-8;
            double deficit = target - sign * area;
            if (deficit <= 1.0e-10) continue;

            double d_area[3] = {
                qx * (p[1][1] - p[2][1]) -
                    qy * (p[1][0] - p[2][0]),
                qx * (p[2][1] - p[0][1]) -
                    qy * (p[2][0] - p[0][0]),
                qx * (p[0][1] - p[1][1]) -
                    qy * (p[0][0] - p[1][0])
            };
            double denom = 0.0;
            for (int k = 0; k < 3; k++) {
                size_t vi = (size_t)tri[k];
                if (site_of[vi] < 0 || (fixed && fixed[vi])) {
                    d_area[k] = 0.0;
                    continue;
                }
                d_area[k] *= sign;
                denom += d_area[k] * d_area[k];
            }
            if (denom <= 1.0e-20) {
                free(initial);
                return -1;
            }
            double step = deficit / denom;
            for (int k = 0; k < 3; k++) {
                if (d_area[k] == 0.0) continue;
                size_t vi = (size_t)tri[k];
                size_t site = (size_t)site_of[vi];
                height[vi] += (float)(step * d_area[k]);
                for (int d = 0; d < 3; d++)
                    candidate->verts[vi * 3 + (size_t)d] =
                        repair_base[site * 3 + (size_t)d] +
                        height[vi] * (float)normal[d];
            }
            violated++;
            if (deficit > worst) worst = deficit;
        }
        if (violated == 0 || worst <= 1.0e-9) {
            feasible = 1;
            break;
        }
    }

    size_t adjusted = 0;
    double max_adjustment = 0.0;
    for (size_t site = 0; site < n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        double change = fabs((double)height[vi] - initial[site]);
        if (change > 1.0e-6) adjusted++;
        if (change > max_adjustment) max_adjustment = change;
    }
    free(initial);
    if (out_adjusted) *out_adjusted = adjusted;
    if (out_max_adjustment) *out_max_adjustment = max_adjustment;
    return feasible ? 0 : -1;
}

/* All-patch quotient-surface fit.  Unlike SurfaceSnap, this is not a dark-
 * defect classifier: every replacement vertex participates.  The source
 * positions define a harmonic boundary-height field, while a multi-label GCO
 * chooses one coherent raw-supported residual depth over that field. */
static int merge_raw_volume_fit_mrf(Arena_T arena,
                                    ComponentMesh *candidate,
                                    const ComponentMesh *source,
                                    const int32_t *face_labels,
                                    int32_t num_labels,
                                    const int32_t *source_adj_fa,
                                    const int32_t *source_adj_fb,
                                    size_t source_n_adj,
                                    const OverlapSepOptions *options,
                                    const OverlapMergeStats *merge,
                                    MergeAuditFrame *audit_frame,
                                    MergeRawFitStats *stats)
{
    memset(stats, 0, sizeof(*stats));
    if (!candidate || !candidate->verts || !candidate->faces || !source ||
        !source->verts || !source->faces || !face_labels ||
        num_labels <= 0 || !source_adj_fa || !source_adj_fb ||
        !options || !options->raw_dir || !*options->raw_dir || !merge ||
        !audit_frame || !merge->vertex_role ||
        !merge->exterior_source_faces ||
        merge->patch_face_start >= candidate->nf)
        return -1;

    const size_t nv = candidate->nv;
    const int labels_n = 33;
    const int center_label = labels_n / 2;
    double reach = options->raw_snap_reach > 0.0f
        ? options->raw_snap_reach : 8.0;
    double depth_step = (2.0 * reach) / (double)(labels_n - 1);
    double normal[3] = {
        merge->plane_normal[0], merge->plane_normal[1],
        merge->plane_normal[2]
    };
    double nn = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                     normal[2] * normal[2]);
    if (nn <= 1.0e-9) return -1;
    for (int d = 0; d < 3; d++) normal[d] /= nn;

    int32_t *site_of = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*site_of));
    for (size_t vi = 0; vi < nv; vi++) site_of[vi] = -1;
    for (size_t vi = 0; vi < nv; vi++) {
        if (merge->vertex_role[vi] == OVERLAP_MERGE_VERTEX_REPAIR)
            site_of[vi] = (int32_t)stats->n_free++;
    }
    if (stats->n_free < 3 || stats->n_free > (size_t)INT_MAX) return -1;
    int32_t *vertex_of = (int32_t *)ARENA_ALLOC(
        arena, stats->n_free * sizeof(*vertex_of));
    for (size_t vi = 0; vi < nv; vi++)
        if (site_of[vi] >= 0) vertex_of[(size_t)site_of[vi]] = (int32_t)vi;

    merge_translate_for_prediction(candidate, options->world_offset, 1.0f);
    double center_world[3] = {
        merge->plane_center[0] + options->world_offset[0],
        merge->plane_center[1] + options->world_offset[1],
        merge->plane_center[2] + options->world_offset[2]
    };

    CSR_T adjacency = CSR_from_faces(
        arena, candidate->faces, candidate->nf, candidate->nv);
    const int32_t *adj_offset = CSR_offset(adjacency);
    const int32_t *adj_target = CSR_target(adjacency);
    float *height = (float *)ARENA_CALLOC(arena, nv, sizeof(*height));
    float *next_height = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*next_height));
    size_t n_anchor = 0;
    double anchor_height_sum = 0.0;
    for (size_t vi = 0; vi < nv; vi++) {
        const float *p = &candidate->verts[vi * 3];
        height[vi] = (float)(((double)p[0] - center_world[0]) * normal[0] +
                             ((double)p[1] - center_world[1]) * normal[1] +
                             ((double)p[2] - center_world[2]) * normal[2]);
        if (merge->vertex_role[vi] == OVERLAP_MERGE_VERTEX_REPAIR)
            height[vi] = 0.0f;
        else if (merge->vertex_role[vi] == OVERLAP_MERGE_VERTEX_ANCHOR) {
            n_anchor++;
            anchor_height_sum += height[vi];
        }
    }
    if (n_anchor < 3) {
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }

    float mean_anchor_height = (float)(anchor_height_sum / (double)n_anchor);
    for (size_t site = 0; site < stats->n_free; site++)
        height[(size_t)vertex_of[site]] = mean_anchor_height;

    /* Position-only Dirichlet solution is our guaranteed conservative base.
     * Keep it separately so the tangent-clamped solution can be backed off
     * continuously if the production projection would cease to be injective. */
    for (int iteration = 0; iteration < 800; iteration++) {
        memcpy(next_height, height, nv * sizeof(*next_height));
        double max_delta = 0.0;
        for (size_t site = 0; site < stats->n_free; site++) {
            size_t vi = (size_t)vertex_of[site];
            double sum = 0.0;
            int count = 0;
            for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
                sum += height[adj_target[e]];
                count++;
            }
            if (count == 0) continue;
            float value = (float)(sum / (double)count);
            double delta = fabs((double)value - height[vi]);
            if (delta > max_delta) max_delta = delta;
            next_height[vi] = value;
        }
        float *swap = height;
        height = next_height;
        next_height = swap;
        if (max_delta < 1.0e-5) break;
    }
    float *dirichlet_height = (float *)ARENA_ALLOC(
        arena, nv * sizeof(*dirichlet_height));
    memcpy(dirichlet_height, height, nv * sizeof(*dirichlet_height));

    /* Minimize squared graph-Laplacian curvature over the replacement while
     * keeping the retained sheet fixed.  The energy includes the fixed anchor
     * one-ring, so its residual carries the exterior tangent across the seam;
     * unlike a membrane/Dirichlet solve, this is a thin-plate fairing. */
    uint8_t *fair_active = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*fair_active));
    uint8_t *clamped_ring = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*clamped_ring));
    double *interface_target_sum = (double *)ARENA_CALLOC(
        arena, nv, sizeof(*interface_target_sum));
    double *interface_target_weight = (double *)ARENA_CALLOC(
        arena, nv, sizeof(*interface_target_weight));
    int32_t *fair_adj_fa = NULL, *fair_adj_fb = NULL;
    double *fair_adj_len = NULL;
    size_t fair_n_adj = 0;
    double fair_avg_edge = 0.0;
    build_face_adjacency(
        arena, candidate->faces, candidate->nf, merge->chart_uv,
        candidate->nv, &fair_adj_fa, &fair_adj_fb, &fair_adj_len,
        &fair_n_adj, &fair_avg_edge);
    (void)fair_adj_len;
    (void)fair_avg_edge;
    for (size_t edge = 0; edge < fair_n_adj; edge++) {
        int32_t fa = fair_adj_fa[edge], fb = fair_adj_fb[edge];
        int patch_a = (size_t)fa >= merge->patch_face_start;
        int patch_b = (size_t)fb >= merge->patch_face_start;
        if (patch_a == patch_b) continue;
        int32_t fp = patch_a ? fa : fb;
        int32_t fe = patch_a ? fb : fa;
        const int32_t *pt = &candidate->faces[(size_t)fp * 3];
        const int32_t *et = &candidate->faces[(size_t)fe * 3];
        int32_t interior = -1;
        for (int k = 0; k < 3; k++) {
            int shared = 0;
            for (int j = 0; j < 3; j++)
                if (pt[k] == et[j]) shared = 1;
            if (!shared && merge->vertex_role[pt[k]] ==
                               OVERLAP_MERGE_VERTEX_REPAIR)
                interior = pt[k];
        }
        if (interior < 0) continue;

        int32_t e0 = et[0], e1 = et[1], e2 = et[2];
        double u0 = merge->chart_uv[(size_t)e0 * 2];
        double v0 = merge->chart_uv[(size_t)e0 * 2 + 1];
        double du1 = (double)merge->chart_uv[(size_t)e1 * 2] - u0;
        double dv1 = (double)merge->chart_uv[(size_t)e1 * 2 + 1] - v0;
        double du2 = (double)merge->chart_uv[(size_t)e2 * 2] - u0;
        double dv2 = (double)merge->chart_uv[(size_t)e2 * 2 + 1] - v0;
        double det = du1 * dv2 - dv1 * du2;
        if (fabs(det) <= 1.0e-10) continue;
        double dh1 = (double)height[e1] - height[e0];
        double dh2 = (double)height[e2] - height[e0];
        double gu = (dh1 * dv2 - dh2 * dv1) / det;
        double gv = (du1 * dh2 - du2 * dh1) / det;
        double slope = sqrt(gu * gu + gv * gv);
        if (slope > 4.0) {
            gu *= 4.0 / slope;
            gv *= 4.0 / slope;
        }
        double du = (double)merge->chart_uv[(size_t)interior * 2] - u0;
        double dv = (double)merge->chart_uv[(size_t)interior * 2 + 1] - v0;
        double predicted = height[e0] + gu * du + gv * dv;
        double weight = rawfit_shared_edge_length(
            candidate, (size_t)fp, (size_t)fe);
        if (weight <= 1.0e-8) weight = 1.0;
        interface_target_sum[interior] += weight * predicted;
        interface_target_weight[interior] += weight;
    }

    size_t n_interface_repair = 0;
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        fair_active[vi] = 1;
        for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
            int32_t vj = adj_target[e];
            fair_active[vj] = 1;
        }
        if (interface_target_weight[vi] > 0.0) {
            height[vi] = (float)(interface_target_sum[vi] /
                                 interface_target_weight[vi]);
            clamped_ring[vi] = 1;
            n_interface_repair++;
        }
    }
    stats->n_clamped = n_interface_repair;
    float *laplacian = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*laplacian));
    const double fair_step = 0.05;
    for (int iteration = 0; iteration < 3000; iteration++) {
        for (size_t vi = 0; vi < nv; vi++) {
            if (!fair_active[vi]) continue;
            int32_t degree = adj_offset[vi + 1] - adj_offset[vi];
            if (degree <= 0) {
                laplacian[vi] = 0.0f;
                continue;
            }
            double average = 0.0;
            for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++)
                average += height[adj_target[e]];
            laplacian[vi] = (float)(height[vi] - average / degree);
        }
        memcpy(next_height, height, nv * sizeof(*next_height));
        double max_delta = 0.0;
        for (size_t site = 0; site < stats->n_free; site++) {
            size_t vi = (size_t)vertex_of[site];
            if (clamped_ring[vi]) continue;
            double gradient = laplacian[vi];
            for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
                int32_t vj = adj_target[e];
                if (!fair_active[vj]) continue;
                int32_t degree = adj_offset[vj + 1] - adj_offset[vj];
                if (degree > 0) gradient -= laplacian[vj] / degree;
            }
            double value = (double)height[vi] - fair_step * gradient;
            double delta = fabs(value - height[vi]);
            if (delta > max_delta) max_delta = delta;
            next_height[vi] = (float)value;
        }
        float *swap = height;
        height = next_height;
        next_height = swap;
        if (max_delta < 1.0e-6) break;
    }
    float *repair_base = (float *)ARENA_ALLOC(
        arena, stats->n_free * 3 * sizeof(*repair_base));
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        const float *p = &candidate->verts[vi * 3];
        double base_height =
            ((double)p[0] - center_world[0]) * normal[0] +
            ((double)p[1] - center_world[1]) * normal[1] +
            ((double)p[2] - center_world[2]) * normal[2];
        for (int d = 0; d < 3; d++)
            repair_base[site * 3 + (size_t)d] =
                p[d] - (float)(base_height * normal[d]);
    }
    int base_reference_sign = rawfit_patch_reference_sign(
        candidate->verts, candidate->faces, candidate->nf,
        merge->patch_face_start, normal);

    MergeOverlapAuditStats dirichlet_audit, clamped_audit;
    rawfit_apply_height_blend(
        candidate->verts, vertex_of, stats->n_free, repair_base,
        dirichlet_height, height, NULL, normal, 0.0);
    size_t patch_nf = candidate->nf - merge->patch_face_start;
    double *dirichlet_projected_area = (double *)ARENA_ALLOC(
        arena, patch_nf * sizeof(*dirichlet_projected_area));
    for (size_t local_fi = 0; local_fi < patch_nf; local_fi++) {
        size_t fi = merge->patch_face_start + local_fi;
        dirichlet_projected_area[local_fi] =
            rawfit_projected_twice_area(
                candidate->verts, &candidate->faces[fi * 3], audit_frame);
    }
    int dirichlet_clean = merge_candidate_overlap_audit(
        arena, candidate, merge->patch_face_start,
        merge->exterior_source_faces, audit_frame, 0, &dirichlet_audit) &&
        rawfit_patch_no_fold(candidate->verts, candidate->faces,
                             candidate->nf, merge->patch_face_start,
                             normal, base_reference_sign);
    if (!dirichlet_clean) {
        stats->harmonic_patch_overlaps = dirichlet_audit.patch;
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }

    rawfit_apply_height_blend(
        candidate->verts, vertex_of, stats->n_free, repair_base,
        dirichlet_height, height, NULL, normal, 1.0);
    int clamped_clean = merge_candidate_overlap_audit(
        arena, candidate, merge->patch_face_start,
        merge->exterior_source_faces, audit_frame, 0, &clamped_audit) &&
        rawfit_patch_no_fold(candidate->verts, candidate->faces,
                             candidate->nf, merge->patch_face_start,
                             normal, base_reference_sign);
    size_t full_clamped_overlaps = clamped_audit.patch;
    stats->clamped_patch_patch = clamped_audit.patch_patch;
    stats->clamped_patch_exterior = clamped_audit.patch_exterior;
    float *full_clamped_height = (float *)ARENA_ALLOC(
        arena, nv * sizeof(*full_clamped_height));
    memcpy(full_clamped_height, height, nv * sizeof(*full_clamped_height));
    if (g_debug_dir) {
        char path[1024];
        snprintf(path, sizeof(path),
                 "%s/04a_full_clamped_patch.obj", g_debug_dir);
        debug_write_patch_stage(
            arena, path, candidate, merge->patch_face_start,
            options->world_offset);
    }
    const uint8_t *held_boundary = NULL;
    if (!clamped_clean) {
        int constrained = rawfit_enforce_projected_orientation(
            candidate, merge->patch_face_start, dirichlet_projected_area,
            site_of, vertex_of, stats->n_free, repair_base, height,
            clamped_ring, normal, audit_frame,
            &stats->n_injective_adjusted,
            &stats->max_injective_adjustment) == 0;
        if (constrained) {
            clamped_clean = merge_candidate_overlap_audit(
                arena, candidate, merge->patch_face_start,
                merge->exterior_source_faces, audit_frame, 0,
                &clamped_audit) &&
                rawfit_patch_no_fold(
                    candidate->verts, candidate->faces, candidate->nf,
                    merge->patch_face_start, normal, base_reference_sign);
            if (clamped_clean) held_boundary = clamped_ring;
        }
        if (!clamped_clean) {
            memcpy(height, full_clamped_height,
                   nv * sizeof(*full_clamped_height));
            rawfit_apply_height_blend(
                candidate->verts, vertex_of, stats->n_free, repair_base,
                dirichlet_height, height, NULL, normal, 1.0);
        }
    }
    if (!clamped_clean) {
        uint8_t *relax_mask = (uint8_t *)ARENA_CALLOC(
            arena, nv, sizeof(*relax_mask));
        uint8_t *grown_mask = (uint8_t *)ARENA_CALLOC(
            arena, nv, sizeof(*grown_mask));
        size_t local_pairs = merge_candidate_patch_overlap_vertices(
            arena, candidate, merge->patch_face_start, audit_frame,
            relax_mask);
        for (size_t vi = 0; vi < nv; vi++)
            if (site_of[vi] < 0) relax_mask[vi] = 0;

        for (int rings = 0; local_pairs > 0 && rings <= 5; rings++) {
            rawfit_apply_local_height_blend(
                candidate->verts, vertex_of, stats->n_free, repair_base,
                dirichlet_height, full_clamped_height, relax_mask, normal,
                0.0);
            int local_base_clean = merge_candidate_overlap_audit(
                arena, candidate, merge->patch_face_start,
                merge->exterior_source_faces, audit_frame, 0,
                &clamped_audit) &&
                rawfit_patch_no_fold(
                    candidate->verts, candidate->faces, candidate->nf,
                    merge->patch_face_start, normal, base_reference_sign);
            if (local_base_clean) {
                double clean_scale = 0.0, dirty_scale = 1.0;
                for (int iteration = 0; iteration < 14; iteration++) {
                    double trial = 0.5 * (clean_scale + dirty_scale);
                    rawfit_apply_local_height_blend(
                        candidate->verts, vertex_of, stats->n_free,
                        repair_base, dirichlet_height, full_clamped_height,
                        relax_mask, normal, trial);
                    int trial_clean = merge_candidate_overlap_audit(
                        arena, candidate, merge->patch_face_start,
                        merge->exterior_source_faces, audit_frame, 0,
                        &clamped_audit) &&
                        rawfit_patch_no_fold(
                            candidate->verts, candidate->faces,
                            candidate->nf, merge->patch_face_start, normal,
                            base_reference_sign);
                    if (trial_clean) clean_scale = trial;
                    else dirty_scale = trial;
                }
                rawfit_apply_local_height_blend(
                    candidate->verts, vertex_of, stats->n_free, repair_base,
                    dirichlet_height, full_clamped_height, relax_mask,
                    normal, clean_scale);
                clamped_clean = merge_candidate_overlap_audit(
                    arena, candidate, merge->patch_face_start,
                    merge->exterior_source_faces, audit_frame, 0,
                    &clamped_audit) &&
                    rawfit_patch_no_fold(
                        candidate->verts, candidate->faces, candidate->nf,
                        merge->patch_face_start, normal,
                        base_reference_sign);
                if (clamped_clean) {
                    stats->local_relax_rings = rings;
                    stats->local_relax_scale = clean_scale;
                    stats->n_injective_adjusted = 0;
                    stats->max_injective_adjustment = 0.0;
                    for (size_t site = 0; site < stats->n_free; site++) {
                        size_t vi = (size_t)vertex_of[site];
                        height[vi] = full_clamped_height[vi];
                        if (!relax_mask[vi]) continue;
                        height[vi] = (float)(dirichlet_height[vi] +
                            clean_scale *
                            (full_clamped_height[vi] -
                             dirichlet_height[vi]));
                        double change = fabs((double)height[vi] -
                                             full_clamped_height[vi]);
                        if (change > 1.0e-6)
                            stats->n_injective_adjusted++;
                        if (change > stats->max_injective_adjustment)
                            stats->max_injective_adjustment = change;
                    }
                    break;
                }
            }

            memcpy(grown_mask, relax_mask, nv * sizeof(*grown_mask));
            for (size_t vi = 0; vi < nv; vi++) {
                if (!relax_mask[vi] || site_of[vi] < 0) continue;
                for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
                    int32_t vj = adj_target[e];
                    if (site_of[vj] >= 0) grown_mask[vj] = 1;
                }
            }
            memcpy(relax_mask, grown_mask, nv * sizeof(*relax_mask));
        }
        if (!clamped_clean) {
            memcpy(height, full_clamped_height,
                   nv * sizeof(*full_clamped_height));
            rawfit_apply_height_blend(
                candidate->verts, vertex_of, stats->n_free, repair_base,
                dirichlet_height, height, NULL, normal, 1.0);
        }
    }
    double clamp_scale = clamped_clean ? 1.0 : 0.0;
    if (!clamped_clean) {
        rawfit_apply_height_blend(
            candidate->verts, vertex_of, stats->n_free, repair_base,
            dirichlet_height, height, clamped_ring, normal, 0.0);
        int held_base_clean = merge_candidate_overlap_audit(
            arena, candidate, merge->patch_face_start,
            merge->exterior_source_faces, audit_frame, 0,
            &clamped_audit) &&
            rawfit_patch_no_fold(
                candidate->verts, candidate->faces, candidate->nf,
                merge->patch_face_start, normal, base_reference_sign);
        if (held_base_clean) held_boundary = clamped_ring;
        double clean_scale = 0.0, dirty_scale = 1.0;
        for (int iteration = 0; iteration < 12; iteration++) {
            double trial = 0.5 * (clean_scale + dirty_scale);
            rawfit_apply_height_blend(
                candidate->verts, vertex_of, stats->n_free, repair_base,
                dirichlet_height, height, held_boundary, normal, trial);
            int trial_clean = merge_candidate_overlap_audit(
                arena, candidate, merge->patch_face_start,
                merge->exterior_source_faces, audit_frame, 0,
                &clamped_audit) &&
                rawfit_patch_no_fold(
                    candidate->verts, candidate->faces, candidate->nf,
                    merge->patch_face_start, normal, base_reference_sign);
            if (trial_clean) clean_scale = trial;
            else dirty_scale = trial;
        }
        clamp_scale = clean_scale;
    }
    rawfit_apply_height_blend(
        candidate->verts, vertex_of, stats->n_free, repair_base,
        dirichlet_height, height, held_boundary, normal, clamp_scale);
    MergeOverlapAuditStats harmonic_audit;
    int harmonic_projection_clean = merge_candidate_overlap_audit(
        arena, candidate, merge->patch_face_start,
        merge->exterior_source_faces, audit_frame, 0, &harmonic_audit);
    stats->harmonic_patch_overlaps = harmonic_audit.patch;
    stats->clamped_patch_overlaps = full_clamped_overlaps;
    stats->clamp_scale = clamp_scale;
    stats->held_clamp_boundary = held_boundary != NULL;
    if (g_debug_dir) {
        char path[1024];
        snprintf(path, sizeof(path),
                 "%s/04b_projection_safe_base_patch.obj", g_debug_dir);
        debug_write_patch_stage(
            arena, path, candidate, merge->patch_face_start,
            options->world_offset);
    }

    CubeTable ct;
    if (cubetable_init(&ct, arena, options->raw_dir,
                       options->raw_chunk > 0 ? options->raw_chunk : 128,
                       candidate->verts, candidate->nv, reach + 4.0) != 0) {
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }

    float axis_u_f[3], axis_v_f[3];
    float normal_f[3] = {
        (float)normal[0], (float)normal[1], (float)normal[2]
    };
    PCA_orthonormal_basis(normal_f, axis_u_f, axis_v_f);
    double axis_u[3] = {axis_u_f[0], axis_u_f[1], axis_u_f[2]};
    double axis_v[3] = {axis_v_f[0], axis_v_f[1], axis_v_f[2]};

    double anchor_sum = 0.0;
    size_t anchor_samples = 0;
    for (size_t vi = 0; vi < nv; vi++) {
        if (merge->vertex_role[vi] != OVERLAP_MERGE_VERTEX_ANCHOR) continue;
        const float *p = &candidate->verts[vi * 3];
        double sample = sample_trilinear(&ct, p[0], p[1], p[2]);
        if (sample >= 0.0) {
            anchor_sum += sample;
            anchor_samples++;
        }
    }
    double anchor_mean = anchor_samples > 0
        ? anchor_sum / (double)anchor_samples : 128.0;

    size_t sample_count = stats->n_free * (size_t)labels_n;
    float *raw_sample = (float *)ARENA_ALLOC(
        arena, sample_count * sizeof(*raw_sample));
    float *structure = (float *)ARENA_CALLOC(
        arena, sample_count, sizeof(*structure));
    float *valid_values = (float *)ARENA_ALLOC(
        arena, sample_count * sizeof(*valid_values));
    size_t n_valid_values = 0;
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        const float *p = &candidate->verts[vi * 3];
        for (int label = 0; label < labels_n; label++) {
            double depth = ((double)label - center_label) * depth_step;
            double q[3] = {
                p[0] + depth * normal[0],
                p[1] + depth * normal[1],
                p[2] + depth * normal[2]
            };
            float qf[3] = {(float)q[0], (float)q[1], (float)q[2]};
            double sample = sample_vertex(&ct, qf, normal_f, 2.0, 5);
            size_t pos = site * (size_t)labels_n + (size_t)label;
            raw_sample[pos] = (float)sample;
            if (sample >= 0.0) valid_values[n_valid_values++] = (float)sample;
        }
    }
    if (n_valid_values < stats->n_free) {
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }
    qsort(valid_values, n_valid_values, sizeof(*valid_values), compare_float);
    double raw_lo = valid_values[(size_t)(0.05 * (double)(n_valid_values - 1))];
    double raw_hi = valid_values[(size_t)(0.95 * (double)(n_valid_values - 1))];
    if (raw_hi - raw_lo < 5.0) {
        raw_lo = 0.0;
        raw_hi = 255.0;
    }

    uint8_t *anchor_neighbors = (uint8_t *)ARENA_CALLOC(
        arena, stats->n_free, sizeof(*anchor_neighbors));
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
            int32_t vj = adj_target[e];
            if (merge->vertex_role[vj] == OVERLAP_MERGE_VERTEX_ANCHOR &&
                anchor_neighbors[site] < UINT8_MAX)
                anchor_neighbors[site]++;
        }
    }

    int32_t *data = (int32_t *)ARENA_ALLOC(
        arena, sample_count * sizeof(*data));
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        const float *p = &candidate->verts[vi * 3];
        for (int label = 0; label < labels_n; label++) {
            size_t pos = site * (size_t)labels_n + (size_t)label;
            double sample = raw_sample[pos];
            double depth = ((double)label - center_label) * depth_step;
            if (clamped_ring[vi] && label != center_label) {
                data[pos] = 9500;
                continue;
            }
            if (sample < 0.0) {
                data[pos] = 9500;
                continue;
            }
            double brightness = (sample - raw_lo) / (raw_hi - raw_lo);
            if (brightness < 0.0) brightness = 0.0;
            if (brightness > 1.0) brightness = 1.0;
            if (brightness > 0.15) {
                double q[3] = {
                    p[0] + depth * normal[0],
                    p[1] + depth * normal[1],
                    p[2] + depth * normal[2]
                };
                RawTangentTensor tensor;
                if (sample_tangent_tensor(
                        &ct, q, axis_u, axis_v, 2.0, &tensor) == 0) {
                    structure[pos] = (float)(tensor.ridge_center *
                        (0.25 + 0.75 * tensor.coherence));
                }
            }
            double material_cost = 4200.0 * (1.0 - brightness) +
                                   700.0 * (1.0 - structure[pos]);
            double match_cost = 600.0 * fabs(sample - anchor_mean) /
                                (raw_hi - raw_lo);
            double close_cost = 35.0 * fabs(depth);
            double boundary_cost = 0.0;
            if (anchor_neighbors[site] > 0) {
                int delta = abs(label - center_label);
                if (delta > 8) delta = 8;
                boundary_cost = 120.0 * anchor_neighbors[site] * delta;
            }
            if (clamped_ring[vi]) {
                int delta = abs(label - center_label);
                if (delta > 8) delta = 8;
                boundary_cost += 900.0 * delta;
            }
            data[pos] = rawfit_cost(material_cost + match_cost +
                                    close_cost + boundary_cost);
        }
    }

    int32_t *smooth = (int32_t *)ARENA_ALLOC(
        arena, (size_t)labels_n * labels_n * sizeof(*smooth));
    double quilt_weight = options->raw_quilt_weight > 0.0f
        ? options->raw_quilt_weight : 250.0;
    stats->quilt_weight = quilt_weight;
    for (int a = 0; a < labels_n; a++) {
        for (int b = 0; b < labels_n; b++) {
            int delta = abs(a - b);
            if (delta > 8) delta = 8;
            smooth[(size_t)a * labels_n + (size_t)b] =
                rawfit_cost(quilt_weight * delta);
        }
    }
    GCO_Handle gc = GCO_create((int)stats->n_free, labels_n);
    if (!gc) {
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }
    GCO_set_data_cost(gc, data);
    GCO_set_smooth_cost(gc, smooth);
    for (size_t site = 0; site < stats->n_free; site++) {
        int best_label = 0;
        int32_t best_cost = data[site * (size_t)labels_n];
        for (int label = 1; label < labels_n; label++) {
            int32_t cost = data[site * (size_t)labels_n + (size_t)label];
            if (cost < best_cost) {
                best_cost = cost;
                best_label = label;
            }
        }
        GCO_set_label(gc, (int)site, best_label);
        size_t vi = (size_t)vertex_of[site];
        for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
            int32_t vj = adj_target[e];
            int32_t sj = site_of[vj];
            if (sj > (int32_t)site) GCO_set_neighbor(gc, (int)site, sj, 1);
        }
    }
    stats->mrf_energy = GCO_expansion(gc, -1);
    if (stats->mrf_energy < 0) {
        GCO_destroy(gc);
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }
    int *labels = (int *)ARENA_ALLOC(
        arena, stats->n_free * sizeof(*labels));
    GCO_get_labels(gc, labels, (int)stats->n_free);
    GCO_destroy(gc);

    uint8_t *movable = (uint8_t *)ARENA_CALLOC(
        arena, nv, sizeof(*movable));
    float *target_depth = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*target_depth));

    double raw_sum_before = 0.0;
    double raw_target_sum = 0.0;
    size_t raw_target_count = 0;
    double selected_cost = 0.0;
    float *target_world = (float *)ARENA_ALLOC(
        arena, nv * 3 * sizeof(*target_world));
    memcpy(target_world, candidate->verts, nv * 3 * sizeof(*target_world));
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        const float *p = &candidate->verts[vi * 3];
        double current = sample_trilinear(&ct, p[0], p[1], p[2]);
        if (current >= 0.0) raw_sum_before += current;
        int label = clamped_ring[vi] ? center_label : labels[site];
        size_t pos = site * (size_t)labels_n + (size_t)label;
        double depth = ((double)label - center_label) * depth_step;
        selected_cost += data[pos];
        if (raw_sample[pos] < 0.0) continue;
        raw_target_sum += raw_sample[pos];
        raw_target_count++;
        double brightness = (raw_sample[pos] - raw_lo) / (raw_hi - raw_lo);
        if (brightness >= 0.30 || structure[pos] >= 0.15f)
            stats->n_supported++;
        if (!clamped_ring[vi]) {
            movable[vi] = 1;
            target_depth[vi] = (float)depth;
            stats->n_targeted++;
        }
        for (int d = 0; d < 3; d++) {
            target_world[vi * 3 + (size_t)d] +=
                (float)(depth * normal[d]);
        }
    }
    stats->mean_raw_before = raw_sum_before / (double)stats->n_free;
    stats->mean_raw_target = raw_target_count > 0
        ? raw_target_sum / (double)raw_target_count : 0.0;
    stats->mean_data_cost = selected_cost / (double)stats->n_free;
    stats->support_fraction = (double)stats->n_supported /
                              (double)stats->n_free;

    if (g_debug_dir) {
        char path[1024];
        float *target_local = (float *)ARENA_ALLOC(
            arena, nv * 3 * sizeof(*target_local));
        float *colors = (float *)ARENA_ALLOC(
            arena, nv * 3 * sizeof(*colors));
        for (size_t vi = 0; vi < nv; vi++) {
            for (int d = 0; d < 3; d++)
                target_local[vi * 3 + (size_t)d] =
                    target_world[vi * 3 + (size_t)d] -
                    options->world_offset[d];
            colors[vi * 3] = colors[vi * 3 + 1] =
                colors[vi * 3 + 2] = 0.35f;
            if (merge->vertex_role[vi] == OVERLAP_MERGE_VERTEX_ANCHOR) {
                colors[vi * 3] = 0.15f;
                colors[vi * 3 + 1] = 0.55f;
                colors[vi * 3 + 2] = 1.0f;
            } else if (merge->vertex_role[vi] ==
                       OVERLAP_MERGE_VERTEX_REPAIR) {
                colors[vi * 3] = 0.95f;
                colors[vi * 3 + 1] = movable[vi]
                    ? 0.20f : 0.05f;
                colors[vi * 3 + 2] = movable[vi]
                    ? 0.85f : 0.05f;
            }
        }
        snprintf(path, sizeof(path),
                 "%s/overlap_merge_raw_targets.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(path, target_local, nv,
                                     candidate->faces, candidate->nf, colors);
        snprintf(path, sizeof(path),
                 "%s/05_raw_depth_targets.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(path, target_local, nv,
                                     candidate->faces, candidate->nf, colors);
    }

    float *solve_base = (float *)ARENA_ALLOC(
        arena, nv * 3 * sizeof(*solve_base));
    memcpy(solve_base, candidate->verts, nv * 3 * sizeof(*solve_base));
    int reference_sign = rawfit_patch_reference_sign(
        solve_base, candidate->faces, candidate->nf,
        merge->patch_face_start, normal);
    /* Solve only the scalar across-sheet residual.  The common chart remains
     * fixed, the clamped first ring remains exactly on its collar tangent, and
     * no unrelated in-plane Laplacian drift can create a false fold. */
    double snap_alpha = options->raw_snap_alpha > 0.0f &&
                        options->raw_snap_alpha < 1.0f
        ? options->raw_snap_alpha : 0.88;
    double data_ratio = snap_alpha / (1.0 - snap_alpha);
    float *residual = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*residual));
    float *next_residual = (float *)ARENA_CALLOC(
        arena, nv, sizeof(*next_residual));
    for (int iteration = 0; iteration < 1000; iteration++) {
        memcpy(next_residual, residual, nv * sizeof(*next_residual));
        double max_delta = 0.0;
        for (size_t site = 0; site < stats->n_free; site++) {
            size_t vi = (size_t)vertex_of[site];
            if (!movable[vi]) {
                next_residual[vi] = 0.0f;
                continue;
            }
            int32_t degree = adj_offset[vi + 1] - adj_offset[vi];
            if (degree <= 0) continue;
            double neighbor_sum = 0.0;
            for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++)
                neighbor_sum += residual[adj_target[e]];
            double data_weight = data_ratio * degree;
            double value = (data_weight * target_depth[vi] + neighbor_sum) /
                           (data_weight + degree);
            double delta = fabs(value - residual[vi]);
            if (delta > max_delta) max_delta = delta;
            next_residual[vi] = (float)value;
        }
        float *swap = residual;
        residual = next_residual;
        next_residual = swap;
        if (max_delta < 1.0e-5) break;
    }
    double displacement_sum = 0.0;
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        double displacement = fabs(residual[vi]);
        if (displacement > 1.0e-6) {
            stats->n_moved++;
            displacement_sum += displacement;
            if (displacement > stats->max_displacement)
                stats->max_displacement = displacement;
        }
        for (int d = 0; d < 3; d++)
            candidate->verts[vi * 3 + (size_t)d] =
                solve_base[vi * 3 + (size_t)d] +
                residual[vi] * (float)normal[d];
    }
    stats->mean_displacement = stats->n_moved > 0
        ? displacement_sum / (double)stats->n_moved : 0.0;

    float *solve_full = (float *)ARENA_ALLOC(
        arena, nv * 3 * sizeof(*solve_full));
    memcpy(solve_full, candidate->verts, nv * 3 * sizeof(*solve_full));
    double accepted_scale = 1.0;
    MergeOverlapAuditStats fitted_audit;
    int projection_clean = merge_candidate_overlap_audit(
        arena, candidate, merge->patch_face_start,
        merge->exterior_source_faces, audit_frame, 0, &fitted_audit);
    int fold_clean = rawfit_patch_no_fold(
        candidate->verts, candidate->faces, candidate->nf,
        merge->patch_face_start, normal, reference_sign);
    stats->raw_full_patch_overlaps = fitted_audit.patch;

    /* A coherent raw target can still create one or two non-local projected
     * intersections.  Keep the raw fit everywhere else and relax only the
     * exact colliding repair faces (plus the minimum required rings) toward
     * the already-safe collar field. */
    if ((!projection_clean || !fold_clean) && harmonic_projection_clean) {
        uint8_t *relax_mask = (uint8_t *)ARENA_CALLOC(
            arena, nv, sizeof(*relax_mask));
        uint8_t *grown_mask = (uint8_t *)ARENA_CALLOC(
            arena, nv, sizeof(*grown_mask));
        size_t local_pairs = merge_candidate_patch_overlap_vertices(
            arena, candidate, merge->patch_face_start, audit_frame,
            relax_mask);
        for (size_t vi = 0; vi < nv; vi++)
            if (site_of[vi] < 0) relax_mask[vi] = 0;

        for (int rings = 0; local_pairs > 0 && rings <= 5; rings++) {
            rawfit_apply_local_position_blend(
                candidate->verts, vertex_of, stats->n_free, solve_base,
                solve_full, relax_mask, 0.0);
            int local_base_clean = merge_candidate_overlap_audit(
                arena, candidate, merge->patch_face_start,
                merge->exterior_source_faces, audit_frame, 0,
                &fitted_audit) &&
                rawfit_patch_no_fold(
                    candidate->verts, candidate->faces, candidate->nf,
                    merge->patch_face_start, normal, reference_sign);
            if (local_base_clean) {
                double clean_scale = 0.0, dirty_scale = 1.0;
                for (int iteration = 0; iteration < 14; iteration++) {
                    double trial = 0.5 * (clean_scale + dirty_scale);
                    rawfit_apply_local_position_blend(
                        candidate->verts, vertex_of, stats->n_free,
                        solve_base, solve_full, relax_mask, trial);
                    int trial_clean = merge_candidate_overlap_audit(
                        arena, candidate, merge->patch_face_start,
                        merge->exterior_source_faces, audit_frame, 0,
                        &fitted_audit) &&
                        rawfit_patch_no_fold(
                            candidate->verts, candidate->faces,
                            candidate->nf, merge->patch_face_start, normal,
                            reference_sign);
                    if (trial_clean) clean_scale = trial;
                    else dirty_scale = trial;
                }
                rawfit_apply_local_position_blend(
                    candidate->verts, vertex_of, stats->n_free, solve_base,
                    solve_full, relax_mask, clean_scale);
                projection_clean = merge_candidate_overlap_audit(
                    arena, candidate, merge->patch_face_start,
                    merge->exterior_source_faces, audit_frame, 0,
                    &fitted_audit);
                fold_clean = rawfit_patch_no_fold(
                    candidate->verts, candidate->faces, candidate->nf,
                    merge->patch_face_start, normal, reference_sign);
                if (projection_clean && fold_clean) {
                    stats->raw_local_relax_rings = rings;
                    stats->raw_local_relax_scale = clean_scale;
                    for (size_t site = 0; site < stats->n_free; site++) {
                        size_t vi = (size_t)vertex_of[site];
                        if (!relax_mask[vi]) continue;
                        double d2 = 0.0;
                        for (int d = 0; d < 3; d++) {
                            double delta = (1.0 - clean_scale) *
                                (solve_full[vi * 3 + (size_t)d] -
                                 solve_base[vi * 3 + (size_t)d]);
                            d2 += delta * delta;
                        }
                        double adjustment = sqrt(d2);
                        if (adjustment > 1.0e-6)
                            stats->raw_local_adjusted++;
                        if (adjustment > stats->raw_local_max_adjustment)
                            stats->raw_local_max_adjustment = adjustment;
                    }
                    break;
                }
            }

            memcpy(grown_mask, relax_mask, nv * sizeof(*grown_mask));
            for (size_t vi = 0; vi < nv; vi++) {
                if (!relax_mask[vi] || site_of[vi] < 0) continue;
                for (int32_t e = adj_offset[vi]; e < adj_offset[vi + 1]; e++) {
                    int32_t vj = adj_target[e];
                    if (site_of[vj] >= 0) grown_mask[vj] = 1;
                }
            }
            memcpy(relax_mask, grown_mask, nv * sizeof(*relax_mask));
        }
    }

    /* The raw MRF establishes the shape of one coherent depth field.  Preserve
     * that field, but choose the largest uniform amplitude whose map remains
     * injective in the production overlap projection.  A dirty harmonic base
     * is not repairable by this line search and is a hard rejection. */
    if ((!projection_clean || !fold_clean) && harmonic_projection_clean) {
        double clean_scale = 0.0;
        double dirty_scale = 1.0;
        for (int iteration = 0; iteration < 12; iteration++) {
            double trial_scale = 0.5 * (clean_scale + dirty_scale);
            for (size_t site = 0; site < stats->n_free; site++) {
                size_t vi = (size_t)vertex_of[site];
                for (int d = 0; d < 3; d++) {
                    candidate->verts[vi * 3 + (size_t)d] =
                        solve_base[vi * 3 + (size_t)d] +
                        (float)(trial_scale *
                        (solve_full[vi * 3 + (size_t)d] -
                         solve_base[vi * 3 + (size_t)d]));
                }
            }
            projection_clean = merge_candidate_overlap_audit(
                arena, candidate, merge->patch_face_start,
                merge->exterior_source_faces, audit_frame, 0,
                &fitted_audit);
            fold_clean = rawfit_patch_no_fold(
                candidate->verts, candidate->faces, candidate->nf,
                merge->patch_face_start, normal, reference_sign);
            if (projection_clean && fold_clean) clean_scale = trial_scale;
            else dirty_scale = trial_scale;
        }
        accepted_scale = clean_scale;
        for (size_t site = 0; site < stats->n_free; site++) {
            size_t vi = (size_t)vertex_of[site];
            for (int d = 0; d < 3; d++) {
                candidate->verts[vi * 3 + (size_t)d] =
                    solve_base[vi * 3 + (size_t)d] +
                    (float)(accepted_scale *
                    (solve_full[vi * 3 + (size_t)d] -
                     solve_base[vi * 3 + (size_t)d]));
            }
        }
        projection_clean = merge_candidate_overlap_audit(
            arena, candidate, merge->patch_face_start,
            merge->exterior_source_faces, audit_frame, 0, &fitted_audit);
        fold_clean = rawfit_patch_no_fold(
            candidate->verts, candidate->faces, candidate->nf,
            merge->patch_face_start, normal, reference_sign);
    }
    stats->line_search_scale = accepted_scale;
    stats->fitted_patch_overlaps = fitted_audit.patch;
    if (!harmonic_projection_clean || !projection_clean || !fold_clean) {
        merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
        return -1;
    }
    stats->n_moved = 0;
    stats->mean_displacement = 0.0;
    stats->max_displacement = 0.0;
    displacement_sum = 0.0;
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        double d2 = 0.0;
        for (int d = 0; d < 3; d++) {
            double delta = candidate->verts[vi * 3 + (size_t)d] -
                           solve_base[vi * 3 + (size_t)d];
            d2 += delta * delta;
        }
        double displacement = sqrt(d2);
        if (displacement <= 1.0e-6) continue;
        stats->n_moved++;
        displacement_sum += displacement;
        if (displacement > stats->max_displacement)
            stats->max_displacement = displacement;
    }
    stats->mean_displacement = stats->n_moved > 0
        ? displacement_sum / (double)stats->n_moved : 0.0;

    double raw_sum_after = 0.0;
    for (size_t site = 0; site < stats->n_free; site++) {
        size_t vi = (size_t)vertex_of[site];
        const float *p = &candidate->verts[vi * 3];
        double sample = sample_trilinear(&ct, p[0], p[1], p[2]);
        if (sample >= 0.0) raw_sum_after += sample;
    }
    stats->mean_raw_after = raw_sum_after / (double)stats->n_free;

    uint8_t *source_patch = (uint8_t *)ARENA_ALLOC(
        arena, source->nf * sizeof(*source_patch));
    memset(source_patch, 1, source->nf * sizeof(*source_patch));
    int source_map_ok = 1;
    for (size_t fi = 0; fi < merge->patch_face_start; fi++) {
        int32_t source_face = merge->exterior_source_faces[fi];
        if (source_face < 0 || (size_t)source_face >= source->nf) {
            source_map_ok = 0;
            break;
        }
        source_patch[source_face] = 0;
    }
    uint8_t *candidate_patch = (uint8_t *)ARENA_CALLOC(
        arena, candidate->nf, sizeof(*candidate_patch));
    for (size_t fi = merge->patch_face_start; fi < candidate->nf; fi++)
        candidate_patch[fi] = 1;

    int32_t *candidate_adj_fa = NULL, *candidate_adj_fb = NULL;
    double *candidate_adj_len = NULL;
    size_t candidate_n_adj = 0;
    double candidate_avg_edge = 0.0;
    build_face_adjacency(
        arena, candidate->faces, candidate->nf, merge->chart_uv,
        candidate->nv, &candidate_adj_fa, &candidate_adj_fb,
        &candidate_adj_len, &candidate_n_adj, &candidate_avg_edge);
    (void)candidate_adj_len;
    (void)candidate_avg_edge;

    MergeRawInterfaceStats source_interface, fitted_interface;
    memset(&source_interface, 0, sizeof(source_interface));
    memset(&fitted_interface, 0, sizeof(fitted_interface));
    const float zero_offset[3] = {0.0f, 0.0f, 0.0f};
    double *fitted_face_energy = g_debug_dir
        ? (double *)ARENA_ALLOC(
              arena, candidate->nf * sizeof(*fitted_face_energy))
        : NULL;
    int source_interface_ok = source_map_ok &&
        merge_raw_interface_metrics(
            arena, &ct, source, source_adj_fa, source_adj_fb, source_n_adj,
            source_patch, options->world_offset, &source_interface, NULL) == 0;
    int fitted_interface_ok = merge_raw_interface_metrics(
        arena, &ct, candidate, candidate_adj_fa, candidate_adj_fb,
        candidate_n_adj, candidate_patch, zero_offset,
        &fitted_interface, fitted_face_energy) == 0;
    int material_ok = source_map_ok && merge_raw_patch_material(
        arena, &ct, source, source_patch, face_labels, num_labels,
        options->world_offset, candidate, merge->patch_face_start,
        &stats->source_best_material, &stats->fitted_material,
        &stats->source_material_support,
        &stats->fitted_material_support) == 0;
    ComponentMesh baseline_candidate = *candidate;
    baseline_candidate.verts = solve_base;
    double baseline_source_ceiling = -1.0;
    double baseline_source_support = 0.0;
    int baseline_material_ok = source_map_ok && merge_raw_patch_material(
        arena, &ct, source, source_patch, face_labels, num_labels,
        options->world_offset, &baseline_candidate,
        merge->patch_face_start, &baseline_source_ceiling,
        &stats->baseline_material, &baseline_source_support,
        &stats->baseline_material_support) == 0;
    if (source_interface_ok) {
        stats->source_interface_pairs = source_interface.pairs;
        stats->source_interface_support = source_interface.support_fraction;
        stats->source_quilt_mean = source_interface.quilt_mean;
        stats->source_quilt_p95 = source_interface.quilt_p95;
        stats->source_bend_p95 = source_interface.bend_p95_degrees;
    }
    if (fitted_interface_ok) {
        stats->fitted_interface_pairs = fitted_interface.pairs;
        stats->fitted_interface_support = fitted_interface.support_fraction;
        stats->fitted_quilt_mean = fitted_interface.quilt_mean;
        stats->fitted_quilt_p95 = fitted_interface.quilt_p95;
        stats->fitted_bend_p95 = fitted_interface.bend_p95_degrees;
    }
    stats->evidence_valid = source_interface_ok && fitted_interface_ok &&
        material_ok && baseline_material_ok &&
        source_interface.support_fraction >= 0.90 &&
        fitted_interface.support_fraction >= 0.90 &&
        stats->source_material_support >= 0.90 &&
        stats->baseline_material_support >= 0.90 &&
        stats->fitted_material_support >= 0.90;
    if (g_debug_dir) {
        char path[1024];
        float *local_verts = (float *)ARENA_ALLOC(
            arena, candidate->nv * 3 * sizeof(*local_verts));
        float *colors = (float *)ARENA_ALLOC(
            arena, candidate->nv * 3 * sizeof(*colors));
        double *energy_sum = (double *)ARENA_CALLOC(
            arena, candidate->nv, sizeof(*energy_sum));
        double *energy_weight = (double *)ARENA_CALLOC(
            arena, candidate->nv, sizeof(*energy_weight));
        for (size_t vi = 0; vi < candidate->nv; vi++) {
            for (int d = 0; d < 3; d++)
                local_verts[vi * 3 + (size_t)d] =
                    candidate->verts[vi * 3 + (size_t)d] -
                    options->world_offset[d];
            colors[vi * 3] = colors[vi * 3 + 1] =
                colors[vi * 3 + 2] = 0.18f;
        }
        if (fitted_face_energy) {
            for (size_t fi = 0; fi < candidate->nf; fi++) {
                if (!isfinite(fitted_face_energy[fi])) continue;
                for (int k = 0; k < 3; k++) {
                    size_t vi = (size_t)candidate->faces[fi * 3 + (size_t)k];
                    energy_sum[vi] += fitted_face_energy[fi];
                    energy_weight[vi] += 1.0;
                }
            }
        }
        double energy_scale = fitted_interface.quilt_p95 > 1.0e-6
            ? fitted_interface.quilt_p95 : 1.0;
        for (size_t vi = 0; vi < candidate->nv; vi++) {
            if (energy_weight[vi] <= 0.0) continue;
            double t = (energy_sum[vi] / energy_weight[vi]) / energy_scale;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            colors[vi * 3] = (float)t;
            colors[vi * 3 + 1] = (float)(1.0 - fabs(2.0 * t - 1.0));
            colors[vi * 3 + 2] = (float)(1.0 - t);
        }
        snprintf(path, sizeof(path),
                 "%s/08_quilting_seam_energy.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(
            path, local_verts, candidate->nv, candidate->faces,
            candidate->nf, colors);

        for (size_t vi = 0; vi < candidate->nv; vi++) {
            colors[vi * 3] = colors[vi * 3 + 1] =
                colors[vi * 3 + 2] = 0.18f;
            if (merge->vertex_role[vi] == OVERLAP_MERGE_VERTEX_EXTERIOR)
                continue;
            const float *p = &candidate->verts[vi * 3];
            double value = sample_vertex(&ct, p, normal_f, 2.0, 5);
            if (value < 0.0) continue;
            double t = (value - raw_lo) / (raw_hi - raw_lo);
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            colors[vi * 3] = (float)t;
            colors[vi * 3 + 1] = (float)(0.15 + 0.85 * t);
            colors[vi * 3 + 2] = (float)(1.0 - t);
        }
        snprintf(path, sizeof(path),
                 "%s/09_raw_material_support.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(
            path, local_verts, candidate->nv, candidate->faces,
            candidate->nf, colors);
    }
    merge_translate_for_prediction(candidate, options->world_offset, -1.0f);
    return 0;
}

static int merge_patch_reference_sign(const ComponentMesh *candidate,
                                      size_t patch_start,
                                      const float normal[3])
{
    double sum = 0.0;
    for (size_t fi = patch_start; fi < candidate->nf; fi++) {
        const int32_t *tri = &candidate->faces[fi * 3];
        const float *a = &candidate->verts[(size_t)tri[0] * 3];
        const float *b = &candidate->verts[(size_t)tri[1] * 3];
        const float *c = &candidate->verts[(size_t)tri[2] * 3];
        double ab[3], ac[3];
        for (int d = 0; d < 3; d++) {
            ab[d] = (double)b[d] - a[d];
            ac[d] = (double)c[d] - a[d];
        }
        double cross[3] = {
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0]
        };
        sum += cross[0] * normal[0] +
               cross[1] * normal[1] +
               cross[2] * normal[2];
    }
    return (sum >= 0.0) ? 1 : -1;
}

static int merge_patch_has_no_fold(const ComponentMesh *candidate,
                                   size_t patch_start,
                                   const float normal[3],
                                   int reference_sign)
{
    for (size_t fi = patch_start; fi < candidate->nf; fi++) {
        const int32_t *tri = &candidate->faces[fi * 3];
        const float *a = &candidate->verts[(size_t)tri[0] * 3];
        const float *b = &candidate->verts[(size_t)tri[1] * 3];
        const float *c = &candidate->verts[(size_t)tri[2] * 3];
        double ab[3], ac[3];
        for (int d = 0; d < 3; d++) {
            ab[d] = (double)b[d] - a[d];
            ac[d] = (double)c[d] - a[d];
        }
        double dot =
            (ab[1] * ac[2] - ab[2] * ac[1]) * normal[0] +
            (ab[2] * ac[0] - ab[0] * ac[2]) * normal[1] +
            (ab[0] * ac[1] - ab[1] * ac[0]) * normal[2];
        if ((double)reference_sign * dot <= 1.0e-10) return 0;
    }
    return 1;
}

static double overlap_env_double(const char *name, double fallback)
{
    const char *text = getenv(name);
    if (!text || !*text) return fallback;
    char *end = NULL;
    double value = strtod(text, &end);
    return (end && end != text) ? value : fallback;
}

static size_t overlap_env_size(const char *name, size_t fallback)
{
    double value = overlap_env_double(name, (double)fallback);
    if (value < 0.0 || value > (double)SIZE_MAX) return fallback;
    return (size_t)value;
}

int OverlapSep_detect_pairs(Arena_T arena, const ComponentMesh *mesh,
                            OverlapSepPairSet *out_pairs)
{
    if (arena == NULL || mesh == NULL || out_pairs == NULL ||
        mesh->verts == NULL || mesh->faces == NULL || mesh->nv < 3 ||
        mesh->nf == 0)
        return -1;
    memset(out_pairs, 0, sizeof *out_pairs);

    float u_axis[3], v_axis[3];
    PCA_orthonormal_basis(mesh->pca_normal, u_axis, v_axis);
    float *proj = (float *)ARENA_ALLOC(
        arena, (long)(mesh->nv * 2 * sizeof(float)));
    for (size_t i = 0; i < mesh->nv; i++) {
        const float *p = &mesh->verts[i * 3];
        proj[i * 2 + 0] = p[0] * u_axis[0] + p[1] * u_axis[1] +
                          p[2] * u_axis[2];
        proj[i * 2 + 1] = p[0] * v_axis[0] + p[1] * v_axis[1] +
                          p[2] * v_axis[2];
    }
    detect_overlaps(arena, mesh->faces, mesh->nf, proj, mesh->nv,
                    &out_pairs->face0, &out_pairs->face1,
                    &out_pairs->count);
    return 0;
}

/* ================================================================
 * Entry point
 * ================================================================ */

int OverlapSep_process(Arena_T              arena,
                       const ComponentMesh  *mesh,
                       int                   sheet_count,
                       int                   n_threads,
                       double                timeout_sec,
                       ComponentMesh       **out_meshes,
                       size_t               *out_count)
{
    return OverlapSep_process_ex(arena, mesh, sheet_count, n_threads,
                                 timeout_sec, NULL, out_meshes, out_count);
}

int OverlapSep_process_ex(Arena_T                 arena,
                          const ComponentMesh     *mesh,
                          int                      sheet_count,
                          int                      n_threads,
                          double                   timeout_sec,
                          const OverlapSepOptions *options,
                          ComponentMesh          **out_meshes,
                          size_t                  *out_count)
{
    assert(arena);
    assert(mesh);
    assert(out_meshes);
    assert(out_count);
    (void)n_threads;

    size_t nv = mesh->nv;
    size_t nf = mesh->nf;

    /* Degenerate input */
    if (nf < OVERLAP_MIN_FACES || nv < 3) {
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
    }

    double t_start = now_sec();
    double t_deadline = t_start + timeout_sec;

    fprintf(stderr, "  step3-overlap: nv=%zu nf=%zu sheet_count=%d\n",
            nv, nf, sheet_count);

    /* ---- Phase 1: Project to 2D ---- */
    float u_axis[3], v_axis[3];
    PCA_orthonormal_basis(mesh->pca_normal, u_axis, v_axis);

    float *proj = (float *)ARENA_ALLOC(arena,
        (long)(nv * 2 * sizeof(float)));
    const float *verts = mesh->verts;
    for (size_t i = 0; i < nv; i++) {
        float px = verts[i * 3 + 0];
        float py = verts[i * 3 + 1];
        float pz = verts[i * 3 + 2];
        proj[i * 2 + 0] = px * u_axis[0] + py * u_axis[1] + pz * u_axis[2];
        proj[i * 2 + 1] = px * v_axis[0] + py * v_axis[1] + pz * v_axis[2];
    }

    fprintf(stderr, "    phase 1 (project): %.3fs\n", now_sec() - t_start);

    /* Debug: dump input mesh + 2D projection */
    if (g_debug_dir) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/overlap_input.obj", g_debug_dir);
        ObjIO_write(path, verts, nv, mesh->faces, nf);
        fprintf(stderr, "    debug: wrote %s\n", path);

        /* Write 2D projection as flat mesh (z=0) */
        float *proj_verts = (float *)ARENA_ALLOC(arena,
            (long)(nv * 3 * sizeof(float)));
        for (size_t vi = 0; vi < nv; vi++) {
            proj_verts[vi * 3 + 0] = proj[vi * 2 + 0];
            proj_verts[vi * 3 + 1] = proj[vi * 2 + 1];
            proj_verts[vi * 3 + 2] = 0.0f;
        }
        snprintf(path, sizeof(path), "%s/overlap_proj2d.obj", g_debug_dir);
        ObjIO_write(path, proj_verts, nv, mesh->faces, nf);
        fprintf(stderr, "    debug: wrote %s\n", path);
    }

    /* ---- Phase 2: Build face adjacency ---- */
    size_t n_adj = 0;
    int32_t *adj_fa = NULL, *adj_fb = NULL;
    double *adj_edge_len = NULL;
    double avg_edge_len = 1.0;

    build_face_adjacency(arena, mesh->faces, nf, proj, nv,
                         &adj_fa, &adj_fb, &adj_edge_len, &n_adj,
                         &avg_edge_len);

    fprintf(stderr, "    phase 2 (adjacency): %zu edges, avg_len=%.4f, %.3fs\n",
            n_adj, avg_edge_len, now_sec() - t_start);

    if (now_sec() > t_deadline) {
        fprintf(stderr, "    timeout after adjacency\n");
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
    }

    /* ---- Phase 3: Detect 2D triangle overlaps ---- */
    size_t n_ovl = 0;
    int32_t *ovl_fa = NULL, *ovl_fb = NULL;

    detect_overlaps(arena, mesh->faces, nf, proj, nv,
                    &ovl_fa, &ovl_fb, &n_ovl);

    fprintf(stderr, "    phase 3 (overlaps): %zu pairs, %.3fs\n",
            n_ovl, now_sec() - t_start);

    if (n_ovl == 0) {
        fprintf(stderr, "    no overlaps -> returning unchanged\n");
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
    }

    if (now_sec() > t_deadline) {
        fprintf(stderr, "    timeout after overlap detection\n");
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
    }

    /* ---- Phase 4: Build solver input + solve multicut ---- */
    size_t num_lifted = n_adj + n_ovl;
    int32_t *lifted_from = (int32_t *)ARENA_ALLOC(arena,
        (long)(num_lifted * sizeof(int32_t)));
    int32_t *lifted_to = (int32_t *)ARENA_ALLOC(arena,
        (long)(num_lifted * sizeof(int32_t)));
    double *lifted_weights = (double *)ARENA_ALLOC(arena,
        (long)(num_lifted * sizeof(double)));

    /* Adjacency edges: positive weight = prefer merge */
    for (size_t i = 0; i < n_adj; i++) {
        lifted_from[i] = adj_fa[i];
        lifted_to[i] = adj_fb[i];
        lifted_weights[i] = adj_edge_len[i] / avg_edge_len;
    }
    /* Overlap edges: negative weight = prefer cut */
    for (size_t i = 0; i < n_ovl; i++) {
        size_t j = n_adj + i;
        lifted_from[j] = ovl_fa[i];
        lifted_to[j] = ovl_fb[i];
        lifted_weights[j] = -(double)OVERLAP_NEG_WEIGHT;
    }

    int32_t *face_labels = (int32_t *)ARENA_ALLOC(arena,
        (long)(nf * sizeof(int32_t)));
    int32_t num_clusters = 0;
    int32_t min_k = (sheet_count > 1) ? (int32_t)sheet_count : 2;

    int rc = LiftedMulticut_kernighan_lin(
        (int32_t)nf,
        (int32_t)n_adj, adj_fa, adj_fb,
        (int32_t)num_lifted, lifted_from, lifted_to, lifted_weights,
        min_k, face_labels, &num_clusters);

    fprintf(stderr, "    phase 4 (multicut): rc=%d, clusters=%d, %.3fs\n",
            rc, num_clusters, now_sec() - t_start);

    /* Diagnostic: raw cluster-size distribution (largest first), BEFORE the
     * phase-5 sliver merge. Distinguishes "multicut separated balanced sheets"
     * (two+ comparable clusters) from "multicut made one dominant cluster +
     * slivers" (one huge, rest tiny -> multicut under-separated). */
    if (rc == 0 && num_clusters > 0) {
        int32_t *dbg_fc = (int32_t *)ARENA_CALLOC(arena,
            (long)num_clusters, (long)sizeof(int32_t));
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t l = face_labels[fi];
            if (l >= 0 && l < num_clusters) dbg_fc[l]++;
        }
        fprintf(stderr, "    cluster sizes (top): ");
        for (int t = 0; t < 8; t++) {
            int32_t bi = -1, bv = 0;
            for (int32_t l = 0; l < num_clusters; l++)
                if (dbg_fc[l] > bv) { bv = dbg_fc[l]; bi = l; }
            if (bi < 0) break;
            fprintf(stderr, "%d ", bv);
            dbg_fc[bi] = -1;
        }
        fprintf(stderr, "(of %zu faces)\n", nf);

        /* Decisive: of the repulsive (overlap, -1e6) edges, how many are CUT
         * (endpoints in different clusters = solver acted on them) vs JOINED
         * (same cluster = solver IGNORED its own -1e6 edge -> stuck/bug)? */
        size_t ov_cut = 0, ov_joined = 0;
        for (size_t i = 0; i < n_ovl; i++) {
            if (face_labels[ovl_fa[i]] != face_labels[ovl_fb[i]]) ov_cut++;
            else ov_joined++;
        }
        fprintf(stderr, "    overlap edges: %zu CUT, %zu JOINED  "
                "(JOINED = -1e6 edge left within a cluster = solver did not act)\n",
                ov_cut, ov_joined);
    }

    if (rc != 0 || num_clusters < 2) {
        fprintf(stderr, "    multicut failed or single cluster\n");
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
    }

    /* Debug: dump face-label-colored mesh */
    if (g_debug_dir) {
        float *fcolors = (float *)ARENA_ALLOC(arena,
            (long)(nv * 3 * sizeof(float)));
        /* Default: gray for unassigned */
        for (size_t vi = 0; vi < nv; vi++) {
            fcolors[vi * 3 + 0] = 0.5f;
            fcolors[vi * 3 + 1] = 0.5f;
            fcolors[vi * 3 + 2] = 0.5f;
        }
        /* Color vertices by face label (last-write-wins) */
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t fl = face_labels[fi];
            size_t ci = (size_t)fl % N_PALETTE;
            for (int j = 0; j < 3; j++) {
                size_t vi = (size_t)mesh->faces[fi * 3 + j];
                fcolors[vi * 3 + 0] = g_palette[ci][0];
                fcolors[vi * 3 + 1] = g_palette[ci][1];
                fcolors[vi * 3 + 2] = g_palette[ci][2];
            }
        }
        char path[1024];
        snprintf(path, sizeof(path), "%s/overlap_face_labels.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(path, verts, nv, mesh->faces, nf, fcolors);
        fprintf(stderr, "    debug: wrote %s (%d clusters)\n", path, num_clusters);
    }

    /* A COMPLETED multicut is ALWAYS applied. The deadline gates STARTING the
     * solve (the pre-multicut check after overlap detection, above) -- it must
     * never discard a finished cut. The old post-solve "timeout after multicut ->
     * return_unchanged" here threw away a valid separation: on a fused stacked-wrap
     * clump the lifted multicut runs ~50-70s, overran the 30s deadline, and the
     * whole stack was handed back UNCUT despite a correct ~6-cluster answer (the
     * "5-sheet monster" bug -- z04480_y02304_x02176). The KL solve is bounded by
     * problem size and terminates; once it has, its result is the cheapest correct
     * separation we will get, so we keep it regardless of the clock. */

    /* ---- Phase 5: Partition faces by label, extract sub-meshes ---- */
    /* The multicut already tells us exactly which faces belong together.
     * Just group faces by label, duplicate boundary vertices, extract. */

    /* Absorb genuinely-isolated slivers into a neighbour -- but NEVER merge two
     * clusters that have an overlap (repulsive) edge between them. Merging
     * across an overlap re-joins the faces the multicut just cut and collapses
     * the separation (the historical bug: a small cluster was merged into its
     * largest adjacent neighbour by face-count ratio alone, re-joining every
     * overlap edge). Union-find over clusters with a bitset overlap-guard;
     * num_clusters is small (tens), so the bitsets are tiny. */
    {
        int32_t nc = num_clusters;
        size_t words = ((size_t)nc + 63) / 64;
        if (words == 0) words = 1;

        int32_t *label_fcount = (int32_t *)ARENA_CALLOC(arena,
            (long)nc, (long)sizeof(int32_t));
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t l = face_labels[fi];
            if (l >= 0 && l < nc) label_fcount[l]++;
        }

        size_t n_merged = 0;
        /* Guard the bitset size; if a mesh ever yields a huge cluster count,
         * skip the sliver merge (keeping the multicut labels is always safe). */
        if (nc > 1 && (size_t)nc * words <= (size_t)(64u * 1024u * 1024u)) {
            /* members[g]: which clusters are in group g. inc[g]: which clusters
             * group g has an overlap edge to. Both at the group-ROOT index. */
            uint64_t *members = (uint64_t *)ARENA_CALLOC(arena,
                (long)((size_t)nc * words), (long)sizeof(uint64_t));
            uint64_t *inc = (uint64_t *)ARENA_CALLOC(arena,
                (long)((size_t)nc * words), (long)sizeof(uint64_t));
            for (int32_t l = 0; l < nc; l++)
                members[(size_t)l * words + (size_t)l / 64] |=
                    (uint64_t)1 << ((size_t)l % 64);
            for (size_t i = 0; i < n_ovl; i++) {
                int32_t a = face_labels[ovl_fa[i]], b = face_labels[ovl_fb[i]];
                if (a >= 0 && b >= 0 && a < nc && b < nc && a != b) {
                    inc[(size_t)a * words + (size_t)b / 64] |= (uint64_t)1 << ((size_t)b % 64);
                    inc[(size_t)b * words + (size_t)a / 64] |= (uint64_t)1 << ((size_t)a % 64);
                }
            }

            int32_t *par = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)nc * sizeof(int32_t)));
            int32_t *gsz = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)nc * sizeof(int32_t)));
            for (int32_t l = 0; l < nc; l++) { par[l] = l; gsz[l] = label_fcount[l]; }
            #define OVL_FIND(R, X) do { (R) = (X); \
                while (par[(R)] != (R)) { par[(R)] = par[par[(R)]]; (R) = par[(R)]; } } while (0)

            /* Unique cross-cluster adjacency pairs (small set). */
            uint64_t *adjpairs = (n_adj > 0) ? (uint64_t *)ARENA_ALLOC(arena,
                (long)(n_adj * sizeof(uint64_t))) : NULL;
            size_t n_ap = 0;
            for (size_t i = 0; i < n_adj; i++) {
                int32_t la = face_labels[adj_fa[i]], lb = face_labels[adj_fb[i]];
                if (la >= 0 && lb >= 0 && la < nc && lb < nc && la != lb) {
                    int32_t lo = la < lb ? la : lb, hi = la < lb ? lb : la;
                    adjpairs[n_ap++] = ((uint64_t)(uint32_t)lo << 32) | (uint64_t)(uint32_t)hi;
                }
            }
            if (n_ap > 1) {
                qsort(adjpairs, n_ap, sizeof(uint64_t), compare_u64);
                size_t u = 0;
                for (size_t i = 0; i < n_ap; i++)
                    if (i == 0 || adjpairs[i] != adjpairs[i - 1]) adjpairs[u++] = adjpairs[i];
                n_ap = u;
            }

            float ratio = (float)MIN_FRAGMENT_VERTS_S3 / 100.0f;
            int changed = 1;
            while (changed) {
                changed = 0;
                for (size_t p = 0; p < n_ap; p++) {
                    int32_t a = (int32_t)(adjpairs[p] >> 32);
                    int32_t b = (int32_t)(adjpairs[p] & 0xFFFFFFFFu);
                    int32_t ra, rb; OVL_FIND(ra, a); OVL_FIND(rb, b);
                    if (ra == rb) continue;
                    int32_t small = (gsz[ra] <= gsz[rb]) ? ra : rb;
                    int32_t big   = (small == ra) ? rb : ra;
                    if ((float)gsz[small] >= ratio * (float)gsz[big]) continue;  /* not a sliver */
                    int safe = 1;  /* overlap edge between the two groups? */
                    for (size_t w = 0; w < words; w++)
                        if (inc[(size_t)small * words + w] & members[(size_t)big * words + w]) {
                            safe = 0; break;
                        }
                    if (!safe) continue;
                    par[small] = big;
                    gsz[big] += gsz[small];
                    for (size_t w = 0; w < words; w++) {
                        members[(size_t)big * words + w] |= members[(size_t)small * words + w];
                        inc[(size_t)big * words + w]     |= inc[(size_t)small * words + w];
                    }
                    changed = 1;
                }
            }

            /* weld diagnostic: of the cross-group adjacency pairs still present,
             * how many are size-eligible to weld, and how many of THOSE are
             * overlap-safe. safe>0 => a weld was missed (bug); safe==0 => the
             * leftover slivers are genuinely overlap-trapped against every
             * neighbour (need a different reattachment than size+overlap). */
            {
                size_t se = 0, ss = 0;
                for (size_t p = 0; p < n_ap; p++) {
                    int32_t a = (int32_t)(adjpairs[p] >> 32);
                    int32_t b = (int32_t)(adjpairs[p] & 0xFFFFFFFFu);
                    int32_t ra, rb; OVL_FIND(ra, a); OVL_FIND(rb, b);
                    if (ra == rb) continue;
                    int32_t sm = (gsz[ra] <= gsz[rb]) ? ra : rb;
                    int32_t bg = (sm == ra) ? rb : ra;
                    if ((float)gsz[sm] >= ratio * (float)gsz[bg]) continue;
                    se++;
                    int safe = 1;
                    for (size_t w = 0; w < words; w++)
                        if (inc[(size_t)sm * words + w] & members[(size_t)bg * words + w]) { safe = 0; break; }
                    if (safe) ss++;
                }
                fprintf(stderr, "    weld diag: %zu size-eligible leftover pairs, %zu overlap-SAFE\n", se, ss);
            }

            for (size_t fi = 0; fi < nf; fi++) {
                int32_t old = face_labels[fi];
                if (old >= 0 && old < nc) {
                    int32_t r; OVL_FIND(r, old);
                    if (r != old) { face_labels[fi] = r; n_merged++; }
                }
            }
            #undef OVL_FIND
        }

        fprintf(stderr, "    phase 5 (merge small, overlap-safe): %zu faces relabeled\n",
                n_merged);
        size_t rejoined = 0;
        for (size_t i = 0; i < n_ovl; i++)
            if (face_labels[ovl_fa[i]] == face_labels[ovl_fb[i]]) rejoined++;
        fprintf(stderr, "    after phase-5 merge: %zu/%zu overlap edges JOINED (must be 0)\n",
                rejoined, n_ovl);
    }

    /* Count distinct labels */
    int32_t *label_fcount = (int32_t *)ARENA_CALLOC(arena,
        (long)num_clusters, (long)sizeof(int32_t));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t l = face_labels[fi];
        if (l >= 0 && l < num_clusters) label_fcount[l]++;
    }

    int32_t n_output = 0;
    for (int32_t l = 0; l < num_clusters; l++) {
        if (label_fcount[l] > 0) n_output++;
    }

    fprintf(stderr, "    phase 5: %d non-empty clusters\n", n_output);

    if (g_debug_dir) {
        uint8_t *actual = (uint8_t *)ARENA_CALLOC(
            arena, nf, sizeof(*actual));
        for (size_t i = 0; i < n_ovl; i++) {
            actual[ovl_fa[i]] = 1;
            actual[ovl_fb[i]] = 1;
        }

        char path[1024];
        size_t actual_nv = 0, actual_nf = 0;
        snprintf(path, sizeof(path),
                 "%s/overlap_actual_triangles.obj", g_debug_dir);
        if (debug_write_face_subset(arena, path, verts, nv,
                                    mesh->faces, nf, actual, NULL,
                                    &actual_nv, &actual_nf) == 0) {
            fprintf(stderr,
                    "    debug: wrote %s (exact detector faces: nv=%zu nf=%zu)\n",
                    path, actual_nv, actual_nf);
        }

        float *projected_verts = (float *)ARENA_ALLOC(
            arena, nv * 3 * sizeof(*projected_verts));
        for (size_t vi = 0; vi < nv; vi++) {
            projected_verts[vi * 3] = proj[vi * 2];
            projected_verts[vi * 3 + 1] = proj[vi * 2 + 1];
            projected_verts[vi * 3 + 2] = 0.0f;
        }
        snprintf(path, sizeof(path),
                 "%s/overlap_actual_triangles_projected.obj", g_debug_dir);
        if (debug_write_face_subset(arena, path, projected_verts, nv,
                                    mesh->faces, nf, actual, NULL,
                                    NULL, NULL) == 0)
            fprintf(stderr, "    debug: wrote %s\n", path);

        uint8_t *piece_faces = (uint8_t *)ARENA_ALLOC(
            arena, nf * sizeof(*piece_faces));
        int piece_index = 0;
        for (int32_t label = 0; label < num_clusters; label++) {
            if (label_fcount[label] == 0) continue;
            memset(piece_faces, 0, nf * sizeof(*piece_faces));
            for (size_t fi = 0; fi < nf; fi++) {
                if (actual[fi] && face_labels[fi] == label)
                    piece_faces[fi] = 1;
            }
            snprintf(path, sizeof(path),
                     "%s/overlap_actual_piece_%d_triangles.obj",
                     g_debug_dir, piece_index);
            size_t piece_nv = 0, piece_nf = 0;
            const float *color = g_palette[(size_t)piece_index % N_PALETTE];
            if (debug_write_face_subset(arena, path, verts, nv,
                                        mesh->faces, nf, piece_faces, color,
                                        &piece_nv, &piece_nf) == 0) {
                fprintf(stderr,
                        "    debug: wrote %s (label=%d nv=%zu nf=%zu)\n",
                        path, label, piece_nv, piece_nf);
            }
            piece_index++;
        }
    }

    if (n_output < 2) {
        fprintf(stderr, "    only %d cluster(s) after merge\n", n_output);
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
    }

    /* Competing small-delamination hypothesis.  The ordinary multicut split is
     * still extracted below (and dumped when debugging); only a candidate that
     * clears the bounded planar construction, ORIGINAL-RAW fit, fold check, and
     * a second run of THIS file's overlap detector replaces it at return. */
    ComponentMesh merge_candidate;
    OverlapMergeStats merge_stats;
    int merge_accepted = 0;
    memset(&merge_candidate, 0, sizeof(merge_candidate));
    memset(&merge_stats, 0, sizeof(merge_stats));
    if (options && options->enable_delamination_merge) {
        if (!options->raw_dir || !*options->raw_dir) {
            fprintf(stderr,
                    "    delamination merge: disabled for this component "
                    "(original raw volume unavailable)\n");
        } else {
            OverlapMergeConfig merge_config;
            OverlapMerge_default_config(&merge_config);
            merge_config.min_rings = (int)overlap_env_double(
                "VES_OVERLAP_MERGE_MIN_RINGS", merge_config.min_rings);
            merge_config.max_rings = (int)overlap_env_double(
                "VES_OVERLAP_MERGE_MAX_RINGS", merge_config.max_rings);
            merge_config.local_search_rings = (int)overlap_env_double(
                "VES_OVERLAP_MERGE_SEARCH_RINGS",
                merge_config.local_search_rings);
            merge_config.max_overlap_pairs = overlap_env_size(
                "VES_OVERLAP_MERGE_MAX_PAIRS",
                merge_config.max_overlap_pairs);
            merge_config.max_seed_faces = overlap_env_size(
                "VES_OVERLAP_MERGE_MAX_SEED_FACES",
                merge_config.max_seed_faces);
            merge_config.max_patch_faces = overlap_env_size(
                "VES_OVERLAP_MERGE_MAX_PATCH_FACES",
                merge_config.max_patch_faces);
            merge_config.max_patch_fraction = overlap_env_double(
                "VES_OVERLAP_MERGE_MAX_FRACTION",
                merge_config.max_patch_fraction);
            merge_config.max_seed_diameter = overlap_env_double(
                "VES_OVERLAP_MERGE_MAX_SEED_DIAMETER",
                merge_config.max_seed_diameter);
            merge_config.max_patch_diameter = overlap_env_double(
                "VES_OVERLAP_MERGE_MAX_PATCH_DIAMETER",
                merge_config.max_patch_diameter);
            merge_config.max_plane_rms = overlap_env_double(
                "VES_OVERLAP_MERGE_MAX_PLANE_RMS",
                merge_config.max_plane_rms);
            merge_config.contact_distance = overlap_env_double(
                "VES_OVERLAP_MERGE_CONTACT_DISTANCE",
                merge_config.contact_distance);
            merge_config.merge_distance = overlap_env_double(
                "VES_OVERLAP_MERGE_DISTANCE",
                merge_config.merge_distance);

            MergeAuditFrame audit_frame;
            memset(&audit_frame, 0, sizeof(audit_frame));
            memcpy(audit_frame.u_axis, u_axis, sizeof(u_axis));
            memcpy(audit_frame.v_axis, v_axis, sizeof(v_axis));
            audit_frame.original_ovl_fa = ovl_fa;
            audit_frame.original_ovl_fb = ovl_fb;
            audit_frame.n_original_ovl = n_ovl;
            int merge_rc = OverlapMerge_try(
                arena, mesh, face_labels, num_clusters,
                adj_fa, adj_fb, n_adj, ovl_fa, ovl_fb, n_ovl,
                &merge_config, &merge_candidate, &merge_stats);

            if (merge_rc == 0) {
                fprintf(stderr,
                        "    delamination merge: planar candidate rings=%d "
                        "-%zu +%zu faces, coalesced=%zu, rms=%.3f\n",
                        merge_stats.rings,
                        merge_stats.patch_faces_removed,
                        merge_stats.patch_faces_added,
                        merge_stats.patch_vertices_coalesced,
                        merge_stats.plane_rms);
                if (g_debug_dir) {
                    char path[1024];
                    snprintf(path, sizeof(path),
                             "%s/overlap_merge_planar.obj", g_debug_dir);
                    ObjIO_write(path, merge_candidate.verts,
                                merge_candidate.nv, merge_candidate.faces,
                                merge_candidate.nf);
                    snprintf(path, sizeof(path),
                             "%s/03_planar_patch.obj", g_debug_dir);
                    debug_write_patch_stage(
                        arena, path, &merge_candidate,
                        merge_stats.patch_face_start, NULL);
                }

                OverlapQualityStats quality_stats;
                int quality_ok = OverlapQuality_improve(
                    arena, &merge_candidate, &merge_stats,
                    g_debug_dir, &quality_stats) == 0;
                if (!quality_ok) {
                    fprintf(stderr,
                            "    delamination mesh quality: rejected (%s); "
                            "retaining multicut split\n",
                            quality_stats.reason ?
                                quality_stats.reason : "unknown");
                } else if (g_debug_dir) {
                    char path[1024];
                    snprintf(path, sizeof(path),
                             "%s/03c_quality_improved_candidate.obj",
                             g_debug_dir);
                    ObjIO_write(path, merge_candidate.verts,
                                merge_candidate.nv, merge_candidate.faces,
                                merge_candidate.nf);
                }

                int reference_sign = merge_patch_reference_sign(
                    &merge_candidate, merge_stats.patch_face_start,
                    merge_stats.plane_normal);
                double before_match = 0.0, before_mean = 0.0, before_p95 = 0.0;
                double after_match = 0.0, after_mean = 0.0, after_p95 = 0.0;
                float pred_safety = options->pred_safety_distance > 0.0f
                    ? options->pred_safety_distance : 8.0f;
                int have_pred = options->pred_vol && options->pred_D > 0 &&
                    options->pred_H > 0 && options->pred_W > 0;
                int fit_before_ok = 0, fit_after_ok = 0;
                if (have_pred) {
                    fit_before_ok = merge_prediction_fit(
                        arena, &merge_candidate, merge_stats.vertex_role,
                        options->pred_vol,
                        options->pred_D, options->pred_H, options->pred_W,
                        options->pred_offset, pred_safety,
                        &before_match, &before_mean,
                        &before_p95) == 0;
                }

                MergeRawFitStats raw_fit;
                memset(&raw_fit, 0, sizeof(raw_fit));
                int raw_fit_ok = quality_ok &&
                    merge_raw_volume_fit_mrf(
                        arena, &merge_candidate, mesh, face_labels,
                        num_clusters, adj_fa, adj_fb, n_adj, options,
                        &merge_stats, &audit_frame, &raw_fit) == 0;

                if (have_pred) {
                    fit_after_ok = merge_prediction_fit(
                        arena, &merge_candidate, merge_stats.vertex_role,
                        options->pred_vol,
                        options->pred_D, options->pred_H, options->pred_W,
                        options->pred_offset, pred_safety,
                        &after_match, &after_mean,
                        &after_p95) == 0;
                }
                PCA_normal(merge_candidate.verts, merge_candidate.nv,
                           merge_candidate.pca_normal,
                           merge_candidate.centroid);
                if (g_debug_dir) {
                    char path[1024];
                    snprintf(path, sizeof(path),
                             "%s/06_fitted_candidate.obj", g_debug_dir);
                    ObjIO_write(path, merge_candidate.verts,
                                merge_candidate.nv, merge_candidate.faces,
                                merge_candidate.nf);
                    snprintf(path, sizeof(path),
                             "%s/06_fitted_patch.obj", g_debug_dir);
                    debug_write_patch_stage(
                        arena, path, &merge_candidate,
                        merge_stats.patch_face_start, NULL);
                }
                int no_fold = merge_patch_has_no_fold(
                    &merge_candidate, merge_stats.patch_face_start,
                    merge_stats.plane_normal, reference_sign);
                int no_overlap = merge_candidate_overlap_free(
                    arena, &merge_candidate,
                    merge_stats.patch_face_start,
                    merge_stats.exterior_source_faces,
                    &audit_frame);
                int raw_safe = raw_fit_ok && raw_fit.n_targeted > 0 &&
                    raw_fit.support_fraction >= 0.90 &&
                    raw_fit.mean_raw_after + 2.0 >= raw_fit.mean_raw_before &&
                    raw_fit.n_reverted * 10 <= raw_fit.n_targeted &&
                    raw_fit.evidence_valid &&
                    raw_fit.fitted_quilt_mean <=
                        1.05 * raw_fit.source_quilt_mean &&
                    raw_fit.fitted_bend_p95 <=
                        raw_fit.source_bend_p95 + 10.0 &&
                    raw_fit.fitted_material >=
                        0.98 * raw_fit.baseline_material;
                /* nnU-Net is a loose guard only: raw fitting may legitimately
                 * correct its local geometry, but must not make agreement
                 * catastrophically worse. */
                int pred_safe = !have_pred ||
                    (fit_before_ok && fit_after_ok &&
                     after_match + 0.15 >= before_match &&
                     after_p95 <= pred_safety + 1.0f);

                fprintf(stderr,
                        "    delamination raw fit: free=%zu clamped=%zu "
                        "targets=%zu "
                        "supported=%zu (%.1f%%), moved=%zu reverted=%zu, "
                        "intensity %.2f -> %.2f (band-target %.2f), "
                        "data=%.1f quilt=%.0f mrf=%lld, "
                        "disp mean/max %.2f/%.2f line=%.3f, "
                        "clamp=%.3f%s (full=%zu: %zu patch-patch/%zu collar; "
                        "injective-adjust=%zu max=%.3f local=%d@%.3f), "
                        "projection patch "
                        "%zu -> full=%zu -> %zu; raw-local=%zu r%d@%.3f "
                        "max=%.3f; "
                        "fold=%s overlap=%s\n",
                        raw_fit.n_free, raw_fit.n_clamped,
                        raw_fit.n_targeted,
                        raw_fit.n_supported,
                        100.0 * raw_fit.support_fraction,
                        raw_fit.n_moved, raw_fit.n_reverted,
                        raw_fit.mean_raw_before, raw_fit.mean_raw_after,
                        raw_fit.mean_raw_target,
                        raw_fit.mean_data_cost, raw_fit.quilt_weight,
                        raw_fit.mrf_energy,
                        raw_fit.mean_displacement, raw_fit.max_displacement,
                        raw_fit.line_search_scale,
                        raw_fit.clamp_scale,
                        raw_fit.held_clamp_boundary ? "+held-boundary" : "",
                        raw_fit.clamped_patch_overlaps,
                        raw_fit.clamped_patch_patch,
                        raw_fit.clamped_patch_exterior,
                        raw_fit.n_injective_adjusted,
                        raw_fit.max_injective_adjustment,
                        raw_fit.local_relax_rings,
                        raw_fit.local_relax_scale,
                        raw_fit.harmonic_patch_overlaps,
                        raw_fit.raw_full_patch_overlaps,
                        raw_fit.fitted_patch_overlaps,
                        raw_fit.raw_local_adjusted,
                        raw_fit.raw_local_relax_rings,
                        raw_fit.raw_local_relax_scale,
                        raw_fit.raw_local_max_adjustment,
                        no_fold ? "no" : "YES",
                        no_overlap ? "zero" : "RESIDUAL");
                if (have_pred) {
                    fprintf(stderr,
                            "    delamination prediction safety only: "
                            "matched %.1f%% -> %.1f%%, mean %.3f -> %.3f, "
                            "p95 %.3f -> %.3f (radius %.1f)\n",
                            100.0 * before_match, 100.0 * after_match,
                            before_mean, after_mean, before_p95, after_p95,
                            pred_safety);
                }
                if (raw_fit.evidence_valid) {
                    fprintf(stderr,
                            "    delamination RAW evidence: seam mean/p95 "
                            "%.3f/%.3f -> %.3f/%.3f (<=5%% mean), "
                            "bend p95 %.2f -> %.2f deg (+10 max), "
                            "material source-ceiling/baseline/fitted "
                            "%.2f/%.2f/%.2f (fit >=98%% baseline); "
                            "support seam %.1f%%/%.1f%% material "
                            "%.1f%%/%.1f%%/%.1f%%\n",
                            raw_fit.source_quilt_mean,
                            raw_fit.source_quilt_p95,
                            raw_fit.fitted_quilt_mean,
                            raw_fit.fitted_quilt_p95,
                            raw_fit.source_bend_p95,
                            raw_fit.fitted_bend_p95,
                            raw_fit.source_best_material,
                            raw_fit.baseline_material,
                            raw_fit.fitted_material,
                            100.0 * raw_fit.source_interface_support,
                            100.0 * raw_fit.fitted_interface_support,
                            100.0 * raw_fit.source_material_support,
                            100.0 * raw_fit.baseline_material_support,
                            100.0 * raw_fit.fitted_material_support);
                } else {
                    fprintf(stderr,
                            "    delamination RAW evidence: unavailable or "
                            "below 90%% source coverage\n");
                }

                if (raw_safe && pred_safe && no_fold && no_overlap) {
                    merge_accepted = 1;
                    if (g_debug_dir) {
                        char path[1024];
                        snprintf(path, sizeof(path),
                                 "%s/overlap_merge_fitted.obj", g_debug_dir);
                        ObjIO_write(path, merge_candidate.verts,
                                    merge_candidate.nv,
                                    merge_candidate.faces,
                                    merge_candidate.nf);
                    }
                } else {
                    fprintf(stderr,
                            "    delamination merge: fitted candidate "
                            "rejected; retaining multicut split\n");
                }
            } else {
                fprintf(stderr,
                        "    delamination merge: ineligible (%s; seed=%zu "
                        "pairs=%zu diameter=%.2f)\n",
                        merge_stats.reason ? merge_stats.reason : "unknown",
                        merge_stats.seed_faces, merge_stats.overlap_pairs,
                        merge_stats.seed_diameter);
            }
        }
    }

    /* Debug: dump vertex-label-colored mesh (from face labels) */
    if (g_debug_dir) {
        float *vcolors = (float *)ARENA_ALLOC(arena,
            (long)(nv * 3 * sizeof(float)));
        for (size_t vi = 0; vi < nv; vi++) {
            vcolors[vi * 3 + 0] = 0.3f;
            vcolors[vi * 3 + 1] = 0.3f;
            vcolors[vi * 3 + 2] = 0.3f;
        }
        for (size_t fi = 0; fi < nf; fi++) {
            int32_t fl = face_labels[fi];
            size_t ci = (size_t)fl % N_PALETTE;
            for (int j = 0; j < 3; j++) {
                size_t vi = (size_t)mesh->faces[fi * 3 + j];
                vcolors[vi * 3 + 0] = g_palette[ci][0];
                vcolors[vi * 3 + 1] = g_palette[ci][1];
                vcolors[vi * 3 + 2] = g_palette[ci][2];
            }
        }
        char path[1024];
        snprintf(path, sizeof(path), "%s/overlap_labels_merged.obj", g_debug_dir);
        ObjIO_write_per_vertex_color(path, verts, nv, mesh->faces, nf, vcolors);
        fprintf(stderr, "    debug: wrote %s\n", path);
    }

    /* Build per-label face lists (CSR-style) */
    int32_t *lf_offset = (int32_t *)ARENA_ALLOC(arena,
        (long)((num_clusters + 1) * (long)sizeof(int32_t)));
    lf_offset[0] = 0;
    for (int32_t l = 0; l < num_clusters; l++)
        lf_offset[l + 1] = lf_offset[l] + label_fcount[l];

    int32_t total_kept = lf_offset[num_clusters];
    int32_t *lf_idx = (int32_t *)ARENA_ALLOC(arena,
        (long)(total_kept * (long)sizeof(int32_t)));
    int32_t *lf_cursor = (int32_t *)ARENA_CALLOC(arena,
        (long)num_clusters, (long)sizeof(int32_t));
    for (size_t fi = 0; fi < nf; fi++) {
        int32_t l = face_labels[fi];
        if (l >= 0 && l < num_clusters) {
            lf_idx[lf_offset[l] + lf_cursor[l]] = (int32_t)fi;
            lf_cursor[l]++;
        }
    }

    /* Extract sub-meshes: each label gets its own vertex remap.
     * Zero-init so pin_mask = NULL unless explicitly set later. */
    ComponentMesh *out = (ComponentMesh *)ARENA_CALLOC(arena,
        (long)n_output, (long)sizeof(ComponentMesh));
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena,
        (long)(nv * sizeof(int32_t)));

    int32_t out_idx = 0;
    for (int32_t l = 0; l < num_clusters; l++) {
        if (label_fcount[l] == 0) continue;

        memset(remap, 0xFF, nv * sizeof(int32_t));
        int32_t new_nv = 0;

        int32_t lf_start = lf_offset[l];
        int32_t lf_end = lf_offset[l + 1];

        /* Assign new vertex indices for this label's faces */
        const int32_t *src_faces = mesh->faces;
        for (int32_t j = lf_start; j < lf_end; j++) {
            int32_t fi = lf_idx[j];
            for (int k = 0; k < 3; k++) {
                int32_t vi = src_faces[fi * 3 + k];
                if (remap[vi] < 0)
                    remap[vi] = new_nv++;
            }
        }

        int32_t new_nf = label_fcount[l];

        float *new_verts = (float *)ARENA_ALLOC(arena,
            (long)((long)new_nv * 3L * (long)sizeof(float)));
        int32_t *new_faces = (int32_t *)ARENA_ALLOC(arena,
            (long)((long)new_nf * 3L * (long)sizeof(int32_t)));

        for (size_t v = 0; v < nv; v++) {
            if (remap[v] >= 0) {
                int32_t ni = remap[v];
                new_verts[ni * 3 + 0] = verts[v * 3 + 0];
                new_verts[ni * 3 + 1] = verts[v * 3 + 1];
                new_verts[ni * 3 + 2] = verts[v * 3 + 2];
            }
        }

        for (int32_t j = 0; j < new_nf; j++) {
            int32_t fi = lf_idx[lf_start + j];
            new_faces[j * 3 + 0] = remap[src_faces[fi * 3 + 0]];
            new_faces[j * 3 + 1] = remap[src_faces[fi * 3 + 1]];
            new_faces[j * 3 + 2] = remap[src_faces[fi * 3 + 2]];
        }

        ComponentMesh *cm = &out[out_idx++];
        cm->verts = new_verts;
        cm->faces = new_faces;
        cm->nv = (size_t)new_nv;
        cm->nf = (size_t)new_nf;
        cm->comp_id = mesh->comp_id;
        cm->nv_pre_fill = 0;
        cm->vert_normals = NULL; /* dropped on overlap split; pca_normal set */
        PCA_normal(new_verts, (size_t)new_nv, cm->pca_normal, cm->centroid);
        cm->self = cm;
    }

    *out_meshes = out;
    *out_count = (size_t)out_idx;

    double t_end = now_sec();
    if (merge_accepted) {
        fprintf(stderr,
                "  step3-overlap: split hypothesis has %zu pieces; "
                "accepted one-sheet delamination merge (%.3fs total)\n",
                *out_count, t_end - t_start);
    } else {
        fprintf(stderr, "  step3-overlap: %zu pieces, %.3fs total\n",
                *out_count, t_end - t_start);
    }

    /* Log per-piece sizes and dump debug OBJ */
    for (size_t i = 0; i < *out_count; i++) {
        fprintf(stderr, "    piece %zu: nv=%zu nf=%zu\n",
                i, (*out_meshes)[i].nv, (*out_meshes)[i].nf);
        if (g_debug_dir) {
            char path[1024];
            snprintf(path, sizeof(path), "%s/overlap_piece_%zu.obj",
                     g_debug_dir, i);
            ObjIO_write(path, (*out_meshes)[i].verts, (*out_meshes)[i].nv,
                        (*out_meshes)[i].faces, (*out_meshes)[i].nf);
            fprintf(stderr, "    debug: wrote %s\n", path);
            if (merge_accepted) {
                snprintf(path, sizeof(path),
                         "%s/review_A_split_piece_%zu.obj", g_debug_dir, i);
                ObjIO_write(path, (*out_meshes)[i].verts,
                            (*out_meshes)[i].nv, (*out_meshes)[i].faces,
                            (*out_meshes)[i].nf);
            }
        }
    }

    if (merge_accepted) {
        ComponentMesh *split_out = *out_meshes;
        size_t split_count = *out_count;
        ComponentMesh *merged_out;
        ARENA_NEW(arena, merged_out);
        *merged_out = merge_candidate;
        /* The temporary mask pinned the untouched exterior during local
         * prediction fitting.  It is not a halo ownership mask and must not
         * escape OverlapSep under that meaning. */
        merged_out->pin_mask = NULL;
        merged_out->self = merged_out;
        *out_meshes = merged_out;
        *out_count = 1;
        if (g_debug_dir) {
            char path[1024];
            snprintf(path, sizeof(path), "%s/overlap_merge_output.obj",
                     g_debug_dir);
            ObjIO_write(path, merged_out->verts, merged_out->nv,
                        merged_out->faces, merged_out->nf);
            fprintf(stderr, "    debug: wrote %s\n", path);
            snprintf(path, sizeof(path),
                     "%s/review_B_healed_component.obj", g_debug_dir);
            ObjIO_write(path, merged_out->verts, merged_out->nv,
                        merged_out->faces, merged_out->nf);
            snprintf(path, sizeof(path),
                     "%s/review_B_healed_patch.obj", g_debug_dir);
            debug_write_patch_stage(
                arena, path, merged_out, merge_stats.patch_face_start, NULL);
            snprintf(path, sizeof(path),
                     "%s/review_split_left_merge_right.obj", g_debug_dir);
            debug_write_split_merge_side_by_side(
                arena, path, split_out, split_count, merged_out);
        }
        fprintf(stderr,
                "    merge output: nv=%zu nf=%zu (split pieces retained "
                "only as debug hypotheses)\n",
                merged_out->nv, merged_out->nf);
    }

    return 0;
}
