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
#include "multicut_wrap.h"
#include "../common/pipeline_constants.h"
#include "../common/pca.h"

#include "../common/obj_io.h"

#include <math.h>
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

    if (n_output < 2) {
        fprintf(stderr, "    only %d cluster(s) after merge\n", n_output);
        return_unchanged(arena, mesh, out_meshes, out_count);
        return -1;
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
    fprintf(stderr, "  step3-overlap: %zu pieces, %.3fs total\n",
            *out_count, t_end - t_start);

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
        }
    }

    return 0;
}
