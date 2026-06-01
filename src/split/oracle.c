#include "../common/ves_platform.h"

#include "oracle.h"
#include "../common/pca.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef VESUVIUS_DEBUG
#include <stdio.h>
#endif

/* Local constants */
#define MIN_ORACLE_VERTS  100
#define MAX_SHEET_COUNT    10
#define MIN_CELL_RATIO   0.10f
#define MAX_RAST_VOXELS  500000000UL

/* ---------- helpers ---------- */

static int cmp_float(const void *a, const void *b)
{
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    if (fa < fb) return -1;
    if (fa > fb) return  1;
    return 0;
}

/* ---------- Oracle_count_sheets ---------- */

int Oracle_count_sheets(Arena_T              arena,
                        const ComponentMesh *mesh,
                        int                  grid_size)
{
    assert(arena);

    /* Degenerate input */
    if (!mesh || mesh->nv < MIN_ORACLE_VERTS || mesh->nf == 0) {
        return 1;
    }
    if (grid_size <= 0) return 1;

    Arena_Mark mark = Arena_save(arena);

    const size_t nv    = mesh->nv;
    const size_t nf    = mesh->nf;
    const float *verts = mesh->verts;
    const int32_t *faces  = mesh->faces;
    const float *normal   = mesh->pca_normal;
    const int    gs       = grid_size;        /* keep as int for arithmetic */
    const size_t n_cells  = (size_t)gs * (size_t)gs;

#ifdef VESUVIUS_DEBUG
    double t0, t1;
    t0 = ves_clock_sec();
#endif

    /* ---- Phase 1: shift vertices non-negative ---- */

    float vmin[3] = { verts[0], verts[1], verts[2] };
    for (size_t i = 1; i < nv; i++) {
        if (verts[i * 3 + 0] < vmin[0]) vmin[0] = verts[i * 3 + 0];
        if (verts[i * 3 + 1] < vmin[1]) vmin[1] = verts[i * 3 + 1];
        if (verts[i * 3 + 2] < vmin[2]) vmin[2] = verts[i * 3 + 2];
    }

    float *sv = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    float vmax[3] = { 0.0f, 0.0f, 0.0f };
    for (size_t i = 0; i < nv; i++) {
        float z = verts[i * 3 + 0] - vmin[0];
        float y = verts[i * 3 + 1] - vmin[1];
        float x = verts[i * 3 + 2] - vmin[2];
        sv[i * 3 + 0] = z;
        sv[i * 3 + 1] = y;
        sv[i * 3 + 2] = x;
        if (z > vmax[0]) vmax[0] = z;
        if (y > vmax[1]) vmax[1] = y;
        if (x > vmax[2]) vmax[2] = x;
    }

    /* Volume dimensions: ceil(max) + 2 for margin (matches Python) */
    size_t D = (size_t)(int)ceilf(vmax[0]) + 2;
    size_t H = (size_t)(int)ceilf(vmax[1]) + 2;
    size_t W = (size_t)(int)ceilf(vmax[2]) + 2;

    size_t total_vox = D * H * W;
    if (total_vox > MAX_RAST_VOXELS) {
        Arena_restore(arena, mark);
        return 1;
    }

    uint8_t *rast = (uint8_t *)ARENA_CALLOC(arena, (long)total_vox, 1L);

    /* ---- Phase 2: rasterize triangles ---- */

    for (size_t fi = 0; fi < nf; fi++) {
        size_t i0 = (size_t)faces[fi * 3 + 0];
        size_t i1 = (size_t)faces[fi * 3 + 1];
        size_t i2 = (size_t)faces[fi * 3 + 2];

        const float *v0 = &sv[i0 * 3];
        const float *v1 = &sv[i1 * 3];
        const float *v2 = &sv[i2 * 3];

        /* Edge lengths */
        float d0_sq = 0.0f, d1_sq = 0.0f, d2_sq = 0.0f;
        for (int k = 0; k < 3; k++) {
            float a = v1[k] - v0[k];
            float b = v2[k] - v0[k];
            float c = v2[k] - v1[k];
            d0_sq += a * a;
            d1_sq += b * b;
            d2_sq += c * c;
        }
        float e0 = sqrtf(d0_sq);
        float e1 = sqrtf(d1_sq);
        float e2 = sqrtf(d2_sq);

        float max_edge = e0;
        if (e1 > max_edge) max_edge = e1;
        if (e2 > max_edge) max_edge = e2;

        int subdiv = (int)ceilf(max_edge);
        if (subdiv < MIN_BARY_SUBDIV) subdiv = MIN_BARY_SUBDIV;

        float step = 1.0f / (float)subdiv;

        for (int ia = 0; ia <= subdiv; ia++) {
            float ba = (float)ia * step;
            int max_ib = subdiv - ia;
            for (int ib = 0; ib <= max_ib; ib++) {
                float bb = (float)ib * step;
                float bc = 1.0f - ba - bb;

                float pz = ba * v0[0] + bb * v1[0] + bc * v2[0];
                float py = ba * v0[1] + bb * v1[1] + bc * v2[1];
                float px = ba * v0[2] + bb * v1[2] + bc * v2[2];

                int iz = (int)roundf(pz);
                int iy = (int)roundf(py);
                int ix = (int)roundf(px);

                if (iz >= 0 && (size_t)iz < D &&
                    iy >= 0 && (size_t)iy < H &&
                    ix >= 0 && (size_t)ix < W) {
                    rast[(size_t)iz * H * W + (size_t)iy * W + (size_t)ix] = 1;
                }
            }
        }
    }

#ifdef VESUVIUS_DEBUG
    t1 = ves_clock_sec();
    fprintf(stderr, "    oracle: rasterize %zux%zux%zu (%.3f s)\n",
            D, H, W,
            t1 - t0);
    t0 = ves_clock_sec();
#endif

    /* ---- Phase 3: project occupied voxels to (u, v, w) ---- */

    float u_axis[3], v_axis[3];
    PCA_orthonormal_basis(normal, u_axis, v_axis);

    /* Count occupied voxels */
    size_t n_occ = 0;
    for (size_t i = 0; i < total_vox; i++) {
        if (rast[i]) n_occ++;
    }

    if (n_occ == 0) {
        Arena_restore(arena, mark);
        return 1;
    }

    /* Allocate projection + cell-index arrays */
    float   *proj_u = (float *)  ARENA_ALLOC(arena, (long)(n_occ * sizeof(float)));
    float   *proj_v = (float *)  ARENA_ALLOC(arena, (long)(n_occ * sizeof(float)));
    float   *proj_w = (float *)  ARENA_ALLOC(arena, (long)(n_occ * sizeof(float)));

    /* Scan volume, project occupied voxels */
    size_t pi = 0;
    for (size_t z = 0; z < D; z++) {
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < W; x++) {
                if (!rast[z * H * W + y * W + x]) continue;
                float fz = (float)z, fy = (float)y, fx = (float)x;
                proj_u[pi] = fz * u_axis[0] + fy * u_axis[1] + fx * u_axis[2];
                proj_v[pi] = fz * v_axis[0] + fy * v_axis[1] + fx * v_axis[2];
                proj_w[pi] = fz * normal[0] + fy * normal[1] + fx * normal[2];
                pi++;
            }
        }
    }

    /* UV range */
    float u_min = proj_u[0], u_max = proj_u[0];
    float v_min = proj_v[0], v_max = proj_v[0];
    for (size_t i = 1; i < n_occ; i++) {
        if (proj_u[i] < u_min) u_min = proj_u[i];
        if (proj_u[i] > u_max) u_max = proj_u[i];
        if (proj_v[i] < v_min) v_min = proj_v[i];
        if (proj_v[i] > v_max) v_max = proj_v[i];
    }

    float u_range = u_max - u_min;
    float v_range = v_max - v_min;
    if (u_range < 1e-6f) u_range = 1.0f;
    if (v_range < 1e-6f) v_range = 1.0f;

    /* Compute cell indices: cu = floor((u - u_min) / u_range * gs), clamped */
    int32_t *cell_idx = (int32_t *)ARENA_ALLOC(arena, (long)(n_occ * sizeof(int32_t)));
    float u_inv = (float)gs / u_range;
    float v_inv = (float)gs / v_range;

    for (size_t i = 0; i < n_occ; i++) {
        int cu = (int)floorf((proj_u[i] - u_min) * u_inv);
        int cv = (int)floorf((proj_v[i] - v_min) * v_inv);
        if (cu < 0) cu = 0;
        if (cu >= gs) cu = gs - 1;
        if (cv < 0) cv = 0;
        if (cv >= gs) cv = gs - 1;
        cell_idx[i] = (int32_t)(cu * gs + cv);
    }

    /* CSR construction: two-pass bin w-values by cell */
    int32_t *offsets = (int32_t *)ARENA_CALLOC(arena, (long)(n_cells + 1),
                                               (long)sizeof(int32_t));
    for (size_t i = 0; i < n_occ; i++) {
        offsets[(size_t)cell_idx[i] + 1]++;
    }
    for (size_t i = 1; i <= n_cells; i++) {
        offsets[i] += offsets[i - 1];
    }

    float   *cell_w  = (float *)  ARENA_ALLOC(arena, (long)(n_occ * sizeof(float)));
    int32_t *cursor   = (int32_t *)ARENA_CALLOC(arena, (long)n_cells,
                                                 (long)sizeof(int32_t));

    for (size_t i = 0; i < n_occ; i++) {
        size_t ci  = (size_t)cell_idx[i];
        size_t pos = (size_t)offsets[ci] + (size_t)cursor[ci];
        cell_w[pos] = proj_w[i];
        cursor[ci]++;
    }

#ifdef VESUVIUS_DEBUG
    t1 = ves_clock_sec();
    fprintf(stderr, "    oracle: raycast %zu occupied, %zu cells (%.3f s)\n",
            n_occ, n_cells,
            t1 - t0);
    t0 = ves_clock_sec();
#endif

    /* ---- Phase 4: per-cell gap counting + aggregation ---- */

    float min_gap = (float)ORACLE_MIN_GAP;
    int32_t sheet_hist[MAX_SHEET_COUNT + 1];
    memset(sheet_hist, 0, sizeof(sheet_hist));
    int32_t total_nonempty = 0;

    for (size_t ci = 0; ci < n_cells; ci++) {
        int32_t cnt = offsets[ci + 1] - offsets[ci];
        if (cnt == 0) continue;
        total_nonempty++;

        if (cnt < 2) {
            sheet_hist[1]++;
            continue;
        }

        /* Sort w-values for this cell */
        qsort(&cell_w[(size_t)offsets[ci]], (size_t)cnt,
              sizeof(float), cmp_float);

        /* Count gaps > min_gap  (strict inequality: gap of exactly 10 is NOT a gap) */
        int n_sheets = 1;
        for (int32_t j = 1; j < cnt; j++) {
            float gap = cell_w[(size_t)offsets[ci] + (size_t)j]
                      - cell_w[(size_t)offsets[ci] + (size_t)(j - 1)];
            if (gap > min_gap) {
                n_sheets++;
                if (n_sheets >= MAX_SHEET_COUNT) break;
            }
        }

        if (n_sheets > MAX_SHEET_COUNT) n_sheets = MAX_SHEET_COUNT;
        sheet_hist[n_sheets]++;
    }

    if (total_nonempty == 0) {
        Arena_restore(arena, mark);
        return 1;
    }

    /* Find highest sheet count where >= 10% of non-empty cells agree */
    int32_t min_cells = (int32_t)((float)total_nonempty * MIN_CELL_RATIO);
    int result = 1;

    for (int s = MAX_SHEET_COUNT; s >= 2; s--) {
        int32_t cells_at_least_s = 0;
        for (int t = s; t <= MAX_SHEET_COUNT; t++) {
            cells_at_least_s += sheet_hist[t];
        }
        if (cells_at_least_s >= min_cells) {
            result = s;
            break;
        }
    }

#ifdef VESUVIUS_DEBUG
    t1 = ves_clock_sec();
    fprintf(stderr, "    oracle: aggregate %d nonempty -> %d sheets (%.3f s)\n",
            (int)total_nonempty, result,
            t1 - t0);
#endif

    Arena_restore(arena, mark);
    return result;
}
