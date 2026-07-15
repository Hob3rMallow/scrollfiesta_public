#include "mls_project.h"
#include "run_ctx.h"

#include "pipeline_constants.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * Halo-deterministic MLS-midpoint projection.
 *
 * Implementation notes:
 *
 * 1. Spatial hash uses a uniform voxel-aligned grid of cell size = radius.
 *    Each cell is keyed by integer voxel coordinates so the grid is
 *    translation-invariant up to the cube origin offset (cube-local coords
 *    enter as floats; floor() of those is the cell index). Two adjacent
 *    cubes that store the same vert at the same cube-local coordinate
 *    therefore put it in the same cell, and queries return the same
 *    neighbour set in both cubes.
 *
 * 2. The 27-cell ball query iterates cells in a fixed (dz, dy, dx) order,
 *    so accumulation order is deterministic given identical inputs. No
 *    qsort on float keys, no hash-iteration order surprises.
 *
 * 3. The Jacobi 3x3 eigensolver mirrors the one in pca.c (kept inline
 *    here to avoid a public header churn; both files can be unified
 *    later if needed).
 */

/* ------------------------------------------------------------------
 * 3x3 symmetric Jacobi eigensolver (double precision, ascending evals).
 * Mirrors src/common/pca.c:jacobi_3x3 — kept private here.
 * ------------------------------------------------------------------ */
static void jacobi_3x3(double a[3][3], double eigenvalues[3],
                       double eigenvectors[3][3])
{
    int i = 0, j = 0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    double m[3][3];
    memcpy(m, a, sizeof(m));

    for (int iter = 0; iter < 100; iter++) {
        double max_off = 0.0;
        int p = 0, q = 1;
        for (i = 0; i < 3; i++) {
            for (j = i + 1; j < 3; j++) {
                double absv = fabs(m[i][j]);
                if (absv > max_off) {
                    max_off = absv;
                    p = i;
                    q = j;
                }
            }
        }
        if (max_off < 1e-15) break;

        double app = m[p][p];
        double aqq = m[q][q];
        double apq = m[p][q];

        double theta;
        double diff = aqq - app;
        if (fabs(diff) < 1e-30) {
            theta = 3.14159265358979323846 / 4.0;
        } else {
            theta = 0.5 * atan2(2.0 * apq, diff);
        }
        double c = cos(theta);
        double s = sin(theta);

        double new_m[3][3];
        memcpy(new_m, m, sizeof(new_m));
        new_m[p][p] = c*c*m[p][p] - 2.0*s*c*m[p][q] + s*s*m[q][q];
        new_m[q][q] = s*s*m[p][p] + 2.0*s*c*m[p][q] + c*c*m[q][q];
        new_m[p][q] = 0.0;
        new_m[q][p] = 0.0;
        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                new_m[r][p] = c*m[r][p] - s*m[r][q];
                new_m[p][r] = new_m[r][p];
                new_m[r][q] = s*m[r][p] + c*m[r][q];
                new_m[q][r] = new_m[r][q];
            }
        }
        memcpy(m, new_m, sizeof(m));

        for (int r = 0; r < 3; r++) {
            double vp = eigenvectors[r][p];
            double vq = eigenvectors[r][q];
            eigenvectors[r][p] = c * vp - s * vq;
            eigenvectors[r][q] = s * vp + c * vq;
        }
    }

    double evals[3] = { m[0][0], m[1][1], m[2][2] };
    int order[3] = { 0, 1, 2 };
    for (i = 1; i < 3; i++) {
        int k = i;
        while (k > 0 && evals[order[k-1]] > evals[order[k]]) {
            int t = order[k-1]; order[k-1] = order[k]; order[k] = t;
            k--;
        }
    }
    double sorted[3][3];
    for (i = 0; i < 3; i++) {
        eigenvalues[i] = evals[order[i]];
        for (j = 0; j < 3; j++) {
            sorted[j][i] = eigenvectors[j][order[i]];
        }
    }
    memcpy(eigenvectors, sorted, sizeof(sorted));
}

/* ------------------------------------------------------------------
 * Wendland C² weight, support radius R.
 *   θ(r) = (1 - r/R)^4 · (4r/R + 1)   if r < R
 *        = 0                            otherwise
 * Standard MLS kernel; smooth (C²), exactly compact, monotonic on [0,R].
 * ------------------------------------------------------------------ */
static inline double wendland(double r, double R)
{
    if (r >= R) return 0.0;
    double t = r / R;
    double s = 1.0 - t;
    double s2 = s * s;
    double s4 = s2 * s2;
    return s4 * (4.0 * t + 1.0);
}

/* ------------------------------------------------------------------
 * Spatial hash on integer cell coords. Cell size = radius (one cell-step
 * is exactly the kernel cutoff). Each cell holds a sorted linked list of
 * vert indices; iteration is by *increasing index* so the accumulation
 * order is identical across cubes when the input verts are identical.
 * ------------------------------------------------------------------ */
typedef struct {
    int32_t cell_z, cell_y, cell_x;   /* integer cell coords (signed) */
    int32_t vert_idx;                 /* index into verts[] */
    float   pos_z, pos_y, pos_x;      /* vertex position (for stable
                                       * cross-cube tiebreak — bit-exact
                                       * positions sort identically in
                                       * cubes A and B). */
} CellEntry;

/* Hash a (cz, cy, cx) triple to a uint64 key. Cantor-pair style on
 * signed-to-unsigned bits — no modulus, so we can sort by key to bucket. */
static inline uint64_t cell_key(int32_t cz, int32_t cy, int32_t cx)
{
    /* Shift signed -> unsigned, then pack into 64 bits with 21 bits per
     * axis (allows ±1e6 cells per axis — way more than we need). */
    uint64_t uz = (uint64_t)(cz + (1 << 20)) & 0x1FFFFF;
    uint64_t uy = (uint64_t)(cy + (1 << 20)) & 0x1FFFFF;
    uint64_t ux = (uint64_t)(cx + (1 << 20)) & 0x1FFFFF;
    return (uz << 42) | (uy << 21) | ux;
}

/* Stable comparator: primary key = cell_key, tiebreak by POSITION
 * lex order (z, y, x). Using positions (not vert_idx) makes the sort
 * cube-invariant: bit-identical positions in cubes A and B sort to the
 * same order regardless of how each cube's MC happened to index them.
 * Final tiebreak by vert_idx in case two verts are at the exact same
 * position (shouldn't happen for MC output, but cheap to handle). */
static int cmp_cell_entry(const void *a, const void *b)
{
    const CellEntry *ea = (const CellEntry *)a;
    const CellEntry *eb = (const CellEntry *)b;
    uint64_t ka = cell_key(ea->cell_z, ea->cell_y, ea->cell_x);
    uint64_t kb = cell_key(eb->cell_z, eb->cell_y, eb->cell_x);
    if (ka < kb) return -1;
    if (ka > kb) return 1;
    if (ea->pos_z < eb->pos_z) return -1;
    if (ea->pos_z > eb->pos_z) return 1;
    if (ea->pos_y < eb->pos_y) return -1;
    if (ea->pos_y > eb->pos_y) return 1;
    if (ea->pos_x < eb->pos_x) return -1;
    if (ea->pos_x > eb->pos_x) return 1;
    if (ea->vert_idx < eb->vert_idx) return -1;
    if (ea->vert_idx > eb->vert_idx) return 1;
    return 0;
}

void MLS_project_verts(Arena_T arena,
                       const float *verts, size_t nv,
                       float radius_vox,
                       const float cell_origin[3],
                       float *out_verts,
                       float *out_normals)
{
    assert(arena);
    assert(verts);
    assert(out_verts);
    assert(radius_vox > 0.0f);

    /* Experiment knob (MLS_RADIUS_VOX): override the LOP/MLS kernel radius at
     * runtime so a smaller radius can be swept without rebuilding. At a tight
     * scroll fold two wraps sit ~1-2 vox apart; the default R=12 spans both and
     * fuses them. Default = the caller's compiled radius. */
    {
        const char *r_env = sf_env("MLS_RADIUS_VOX");
        if (r_env && *r_env) {
            double v = atof(r_env);
            if (v > 0.0) {
                radius_vox = (float)v;
                static int mls_r_logged = 0;
                if (!mls_r_logged) {
                    mls_r_logged = 1;
                    fprintf(stderr,
                        "[MLS] kernel radius override -> %.2f vox (MLS_RADIUS_VOX)\n",
                        (double)radius_vox);
                }
            }
        }
    }

    if (nv == 0) return;

    Arena_Mark mark = Arena_save(arena);

    double R = (double)radius_vox;
    double R2 = R * R;
    double inv_R = 1.0 / R;

    /* Optional global origin offset. Cells become world-aligned (rather
     * than cube-local-aligned) when both cubes pass the same offset
     * (e.g. cube origin in world coords). */
    float ocz = 0.0f, ocy = 0.0f, ocx = 0.0f;
    if (cell_origin) {
        ocz = cell_origin[0];
        ocy = cell_origin[1];
        ocx = cell_origin[2];
    }

    /* -------- Build sorted cell entries (one per vert) -------- */
    CellEntry *entries = (CellEntry *)ARENA_ALLOC(arena,
                            (long)nv * (long)sizeof(CellEntry));
    for (size_t i = 0; i < nv; i++) {
        float z = verts[i * 3 + 0];
        float y = verts[i * 3 + 1];
        float x = verts[i * 3 + 2];
        entries[i].cell_z = (int32_t)floorf((z + ocz) * (float)inv_R);
        entries[i].cell_y = (int32_t)floorf((y + ocy) * (float)inv_R);
        entries[i].cell_x = (int32_t)floorf((x + ocx) * (float)inv_R);
        entries[i].vert_idx = (int32_t)i;
        entries[i].pos_z = z;
        entries[i].pos_y = y;
        entries[i].pos_x = x;
    }
    qsort(entries, nv, sizeof(CellEntry), cmp_cell_entry);

    /* Build a parallel array of cell keys + run-length offsets for
     * O(1) cell lookup. For each unique key, store [start, end) range. */
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena,
                        (long)nv * (long)sizeof(uint64_t));
    for (size_t i = 0; i < nv; i++) {
        keys[i] = cell_key(entries[i].cell_z,
                           entries[i].cell_y,
                           entries[i].cell_x);
    }

    /* Cache: previous-call cell_key -> [start, end) via lazy binary search.
     * For simplicity we do binary search per cell lookup; with ≤ 27 cells
     * per vertex query and log2(nv) ~ 17 comparisons, this is well under
     * neighbour-loop cost. */

    /* -------- Per-vertex MLS midpoint -------- */
    /* Each vertex is independent -- all-local accumulators, read-only cell
     * grid, one write to out_*[i] -- so this parallelizes exactly: output is
     * byte-identical to serial for any thread count/schedule (community
     * issue #3, pscamillo: serial MLS was ~83% of single-cube wall; 9.6x on
     * 16 threads). Signed int index for MSVC OpenMP 2.0 (a per-component
     * cloud is < 2^31 verts by orders of magnitude). Thread count follows
     * omp_set_num_threads(n_threads) in the driver, so grid runs with
     * threads-per-cube=1 stay serial (the fleet already fills the cores). */
    assert(nv <= (size_t)INT_MAX);
    int nv_i = (int)nv;
    int si = 0;   /* declared before the loop: MSVC OpenMP 2.0 form */
#pragma omp parallel for schedule(dynamic, 256)
    for (si = 0; si < nv_i; si++) {
        size_t i = (size_t)si;
        float vz = verts[i * 3 + 0];
        float vy = verts[i * 3 + 1];
        float vx = verts[i * 3 + 2];

        int32_t ciz = (int32_t)floorf((vz + ocz) * (float)inv_R);
        int32_t ciy = (int32_t)floorf((vy + ocy) * (float)inv_R);
        int32_t cix = (int32_t)floorf((vx + ocx) * (float)inv_R);

        /* Accumulate weighted centroid + outer-product moments in double. */
        double w_sum = 0.0;
        double cz_s = 0.0, cy_s = 0.0, cx_s = 0.0;
        double mzz = 0.0, mzy = 0.0, mzx = 0.0;
        double myy = 0.0, myx = 0.0, mxx = 0.0;
        size_t n_used = 0;

        /* Pass 1: weighted centroid */
        for (int dz = -1; dz <= 1; dz++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    uint64_t target = cell_key(ciz + dz, ciy + dy, cix + dx);

                    /* Binary search for first index with keys[idx] >= target. */
                    size_t lo = 0, hi = nv;
                    while (lo < hi) {
                        size_t mid = lo + (hi - lo) / 2;
                        if (keys[mid] < target) lo = mid + 1;
                        else hi = mid;
                    }
                    /* Walk run with matching key. */
                    for (size_t k = lo; k < nv && keys[k] == target; k++) {
                        int32_t j = entries[k].vert_idx;
                        double dvz = (double)verts[j*3+0] - (double)vz;
                        double dvy = (double)verts[j*3+1] - (double)vy;
                        double dvx = (double)verts[j*3+2] - (double)vx;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 >= R2) continue;
                        double r = sqrt(r2);
                        double w = wendland(r, R);
                        if (w <= 0.0) continue;
                        w_sum += w;
                        cz_s += w * (double)verts[j*3+0];
                        cy_s += w * (double)verts[j*3+1];
                        cx_s += w * (double)verts[j*3+2];
                        n_used++;
                    }
                }
            }
        }

        if (n_used < (size_t)MLS_MIN_NEIGHBOURS || w_sum <= 0.0) {
            /* Leave vertex in place; flag normal as undefined. */
            out_verts[i*3+0] = vz;
            out_verts[i*3+1] = vy;
            out_verts[i*3+2] = vx;
            if (out_normals) {
                out_normals[i*3+0] = 0.0f;
                out_normals[i*3+1] = 0.0f;
                out_normals[i*3+2] = 0.0f;
            }
            continue;
        }

        double inv_w = 1.0 / w_sum;
        double cz = cz_s * inv_w;
        double cy = cy_s * inv_w;
        double cx = cx_s * inv_w;

        /* Pass 2: weighted covariance around centroid (need it for normal) */
        if (out_normals) {
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
                            double dvz = (double)verts[j*3+0] - (double)vz;
                            double dvy = (double)verts[j*3+1] - (double)vy;
                            double dvx = (double)verts[j*3+2] - (double)vx;
                            double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                            if (r2 >= R2) continue;
                            double r = sqrt(r2);
                            double w = wendland(r, R);
                            if (w <= 0.0) continue;
                            double pz = (double)verts[j*3+0] - cz;
                            double py = (double)verts[j*3+1] - cy;
                            double px = (double)verts[j*3+2] - cx;
                            mzz += w * pz * pz;
                            mzy += w * pz * py;
                            mzx += w * pz * px;
                            myy += w * py * py;
                            myx += w * py * px;
                            mxx += w * px * px;
                        }
                    }
                }
            }

            double cov[3][3];
            cov[0][0] = mzz * inv_w;
            cov[0][1] = mzy * inv_w; cov[1][0] = cov[0][1];
            cov[0][2] = mzx * inv_w; cov[2][0] = cov[0][2];
            cov[1][1] = myy * inv_w;
            cov[1][2] = myx * inv_w; cov[2][1] = cov[1][2];
            cov[2][2] = mxx * inv_w;

            double eigvals[3];
            double eigvecs[3][3];
            jacobi_3x3(cov, eigvals, eigvecs);

            /* Smallest eigenvalue is in column 0 — that's the sheet normal. */
            double nz_d = eigvecs[0][0];
            double ny_d = eigvecs[1][0];
            double nx_d = eigvecs[2][0];
            double nlen = sqrt(nz_d*nz_d + ny_d*ny_d + nx_d*nx_d);
            if (nlen > 1e-15) {
                nz_d /= nlen;
                ny_d /= nlen;
                nx_d /= nlen;
            } else {
                nz_d = 0.0; ny_d = 0.0; nx_d = 1.0;
            }

            /* Sign convention: orient toward (1,1,1) hemisphere — matches
             * the convention in pca.c so any downstream code that expects
             * the pca.c sign rule still works. */
            double ref_dot = nz_d + ny_d + nx_d;
            if (ref_dot < -1e-9) {
                nz_d = -nz_d; ny_d = -ny_d; nx_d = -nx_d;
            } else if (ref_dot < 1e-9) {
                /* Near-perpendicular to (1,1,1): fall back to
                 * largest-magnitude-positive (matches pca.c). */
                int max_idx = 0;
                double max_abs = fabs(nz_d);
                if (fabs(ny_d) > max_abs) { max_abs = fabs(ny_d); max_idx = 1; }
                if (fabs(nx_d) > max_abs) { max_abs = fabs(nx_d); max_idx = 2; }
                double sgn = (max_idx == 0) ? nz_d
                           : (max_idx == 1) ? ny_d : nx_d;
                if (sgn < 0.0) {
                    nz_d = -nz_d; ny_d = -ny_d; nx_d = -nx_d;
                }
            }

            out_normals[i*3+0] = (float)nz_d;
            out_normals[i*3+1] = (float)ny_d;
            out_normals[i*3+2] = (float)nx_d;

            /* MLS tangent-plane projection: V' = V − ((V−C)·N) N
             * Only moves the vertex *perpendicular* to the local sheet
             * normal. The in-plane component of (V − C) is preserved,
             * so multi-iter LOP doesn't shrink the sheet — every iter
             * removes only the through-thickness residual without
             * smoothing the in-sheet geometry toward the cloud's center
             * of mass. Centroid replacement (the old "V = C" rule) does
             * shrink: each iter applies a smoothing kernel in all three
             * dimensions, and after ~100 iters the sheet collapses to a
             * thin curve along its principal axis (verified empirically:
             * comp004 lost 89% of its y-extent under that rule).
             *
             * Sign-flip invariant: V − d·N and V − (−d)·(−N) are equal,
             * so the (1,1,1) hemisphere bias above doesn't affect the
             * position output. */
            if (nlen > 1e-15) {
                double dvz_c = (double)vz - cz;
                double dvy_c = (double)vy - cy;
                double dvx_c = (double)vx - cx;
                double d_proj = dvz_c * nz_d + dvy_c * ny_d + dvx_c * nx_d;
                out_verts[i*3+0] = (float)((double)vz - d_proj * nz_d);
                out_verts[i*3+1] = (float)((double)vy - d_proj * ny_d);
                out_verts[i*3+2] = (float)((double)vx - d_proj * nx_d);
                continue;
            }
        }

        /* Fallback: out_normals==NULL or degenerate covariance. Use
         * the LOP-μ=0 centroid replacement. NOTE: intermediate multi-
         * iter calls pass out_normals==NULL to skip the covariance
         * pass; if you ever want tangent-plane projection on every
         * iter, plumb the normal computation through unconditionally. */
        out_verts[i*3+0] = (float)cz;
        out_verts[i*3+1] = (float)cy;
        out_verts[i*3+2] = (float)cx;
    }

    /* Validation knob (MLS_PERTURB_EPS): add a deterministic pseudo-random
     * offset in [-eps,+eps] to every output coordinate, per call -- so over
     * the 5-iteration LOP it compounds exactly like a reduced-precision
     * backend's per-pass error. Used to adjudicate the FP32 CUDA drop-in
     * from public issue #3 (measured max |delta| 1.5e-3 vox) by injecting
     * that magnitude into the reference pipeline and measuring what
     * BPA/guards/weld actually do. NOT for production runs. */
    {
        const char *pe = sf_env("MLS_PERTURB_EPS");
        if (pe && *pe) {
            double eps = atof(pe);
            for (size_t vi = 0; eps > 0.0 && vi < nv; vi++) {
                for (int c = 0; c < 3; c++) {
                    /* splitmix64 of (vi,c): deterministic across runs */
                    uint64_t h = (uint64_t)vi * 3u + (uint64_t)c
                               + 0x9E3779B97F4A7C15ULL;
                    h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL;
                    h ^= h >> 27; h *= 0x94D049BB133111EBULL;
                    h ^= h >> 31;
                    double u = ((double)(h >> 11) / 9007199254740992.0)
                             * 2.0 - 1.0;
                    out_verts[vi*3 + (size_t)c] += (float)(u * eps);
                }
            }
        }
    }

    Arena_restore(arena, mark);
}
