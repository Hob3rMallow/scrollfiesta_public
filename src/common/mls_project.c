#include "mls_project.h"

#include "pipeline_constants.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* AVX2 span kernel (see the per-vertex sweep below). Strictly optional:
 * the scalar path is the portable reference (Linux GCC without -mavx2,
 * macOS arm64, etc. all build and run it). MSVC defines __AVX2__ under
 * /arch:AVX2; GCC/Clang under -mavx2. */
#if defined(__AVX2__)
#include <immintrin.h>
#endif

/*
 * Halo-deterministic MLS-midpoint projection.
 *
 * Implementation notes (2026-07 fast engine — same math, new machinery):
 *
 * 1. Spatial index is a hash-CSR uniform grid of cell size = radius / 2
 *    (MLS_CELL_DIV). Cells are keyed by integer voxel-aligned coordinates
 *    floor((coord + cell_origin) / cell), so the grid is world-aligned when
 *    both cubes pass the same cell_origin — two adjacent cubes that see the
 *    same vert at the same world coordinate put it in the same cell. The
 *    smaller cell (R/2 instead of R) shrinks the scanned box around the
 *    query sphere from 6.5x the sphere volume to ~3x, and a conservative
 *    sphere/box cull skips corner cells entirely.
 *
 * 2. The CSR buckets are built by counting sort (count -> prefix -> stable
 *    fill), then each bucket is sorted by POSITION (z, y, x, then vert_idx)
 *    — same ordering rule as the previous qsort comparator. The ball query
 *    iterates cells in a fixed (dz, dy, dx) order, so the accumulation
 *    order is a pure function of world-aligned cells + bit-exact positions:
 *    identical across cubes whose halo neighbourhoods agree. (This is the
 *    cross-cube seam-convergence contract; see mls_project.h.)
 *
 * 3. Centroid and covariance are accumulated in ONE fused sweep as weighted
 *    moments about the query vertex v (all double):
 *        w_sum = SUM w,  s1 = SUM w*d,  M2 = SUM w*d*d^T,   d = p_j - v
 *        c   = v + s1/w_sum
 *        cov = M2/w_sum - (s1/w_sum)(s1/w_sum)^T        (covariance about c)
 *    Algebraically identical to the old two-pass version; |d| <= R keeps
 *    the moment subtraction far from catastrophic cancellation in double.
 *    Fusing halves the neighbour gather + sqrt + Wendland work (the old
 *    pass 2 recomputed everything pass 1 had already computed).
 *
 * 4. Bucket positions are snapshotted into SoA arrays at build time, so the
 *    query loop never reads `verts` for neighbours. This makes the header's
 *    "inputs/outputs may alias" contract genuinely safe (the previous
 *    implementation read verts[j] live during the parallel write-back).
 *
 * 5. The Jacobi 3x3 eigensolver mirrors the one in pca.c (kept inline
 *    here to avoid a public header churn). Measured at ~3% of per-vertex
 *    cost — not worth replacing with a closed-form solver.
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
 * Cell grid geometry.
 *
 * Cell size = radius / MLS_CELL_DIV. The ball query walks cell offsets
 * dz,dy,dx in [-MLS_CELL_DIV, +MLS_CELL_DIV]: a query point anywhere in
 * its own cell reaches at most MLS_CELL_DIV cells away for a radius-R
 * sphere (cell * DIV == R exactly). DIV=2 => 5^3 = 125 candidate cells,
 * most of which the sphere/box cull rejects.
 *
 * MLS_CULL_MARGIN_SQ: the cull may only skip a cell when even its closest
 * corner is beyond R — plus a margin generous enough that float rounding
 * in the cell ASSIGNMENT (floor((coord+origin)*inv_cell) at world scale)
 * can never orphan a true r<R neighbour into a culled cell. Assignment
 * fuzz is <= ~3e-3 vox at coord magnitude ~2e4; the margin admits a shell
 * of sqrt(R^2+0.25)-R (~0.01 vox at R=12) — points there fail the exact
 * per-point r2 test anyway, so the cull affects performance only, never
 * the neighbour set.
 * ------------------------------------------------------------------ */
enum { MLS_CELL_DIV = 2 };
#define MLS_CULL_MARGIN_SQ 0.25

/* Pack a (cz, cy, cx) triple to a uint64 key. Signed-to-unsigned shift,
 * then 21 bits per axis (allows ±1e6 cells per axis — way more than we
 * need). Top bit is always 0, so UINT64_MAX is a safe empty sentinel. */
static inline uint64_t cell_key(int32_t cz, int32_t cy, int32_t cx)
{
    uint64_t uz = (uint64_t)(cz + (1 << 20)) & 0x1FFFFF;
    uint64_t uy = (uint64_t)(cy + (1 << 20)) & 0x1FFFFF;
    uint64_t ux = (uint64_t)(cx + (1 << 20)) & 0x1FFFFF;
    return (uz << 42) | (uy << 21) | ux;
}

#define MLS_HASH_EMPTY UINT64_MAX

/* splitmix64 finalizer — same mixer as the MLS_PERTURB_EPS knob below. */
static inline uint64_t hash64(uint64_t h)
{
    h ^= h >> 30; h *= 0xBF58476D1CE4E5B9ULL;
    h ^= h >> 27; h *= 0x94D049BB133111EBULL;
    h ^= h >> 31;
    return h;
}

/* Linear-probe lookup. Returns slot index, or SIZE_MAX if absent.
 * T is a power of two. */
static inline size_t hash_find(const uint64_t *tab_key, size_t T, uint64_t key)
{
    size_t h = (size_t)(hash64(key) & (uint64_t)(T - 1));
    for (;;) {
        uint64_t kh = tab_key[h];
        if (kh == key) return h;
        if (kh == MLS_HASH_EMPTY) return SIZE_MAX;
        h = (h + 1) & (T - 1);
    }
}

#if defined(__AVX2__)
/* Lane-participation masks for the AVX2 span tail: lane l is live iff
 * l < remaining. Row 0 = full batch (remaining >= 4). */
static const uint64_t MLS_TAIL_MASK[4][4] = {
    { UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX },
    { UINT64_MAX, 0,          0,          0          },
    { UINT64_MAX, UINT64_MAX, 0,          0          },
    { UINT64_MAX, UINT64_MAX, UINT64_MAX, 0          },
};
#endif

/* Insertion sort of one CSR bucket [lo, hi) by (pos_z, pos_y, pos_x,
 * vert_idx) ascending — the same ordering rule as the old qsort
 * comparator. Buckets are small (a cell is (R/2)^3 vox of a ~1-layer
 * sheet), so insertion sort beats anything with setup cost.
 *
 * The float compares are EXACT on purpose: this is a bit-value sort for
 * cross-cube-stable ordering, not a geometric tolerance test. */
static void bucket_sort_pos(float *bz, float *by, float *bx, int32_t *bidx,
                            int32_t lo, int32_t hi)
{
    for (int32_t i = lo + 1; i < hi; i++) {
        float z = bz[i], y = by[i], x = bx[i];
        int32_t id = bidx[i];
        int32_t j = i - 1;
        while (j >= lo) {
            int gt;
            if (bz[j] != z)      gt = (bz[j] > z);
            else if (by[j] != y) gt = (by[j] > y);
            else if (bx[j] != x) gt = (bx[j] > x);
            else                 gt = (bidx[j] > id);
            if (!gt) break;
            bz[j+1] = bz[j]; by[j+1] = by[j]; bx[j+1] = bx[j];
            bidx[j+1] = bidx[j];
            j--;
        }
        bz[j+1] = z; by[j+1] = y; bx[j+1] = x; bidx[j+1] = id;
    }
}

/* One existing neighbour cell of an occupied cell: offset and CSR span. */
typedef struct {
    int32_t start, cnt;
    int8_t  dz, dy, dx, pad;
} MlsNb;

static inline int32_t cell_key_axis(uint64_t key, int shift)
{
    return (int32_t)((key >> shift) & 0x1FFFFF) - (1 << 20);
}

static const float *g_mls_point_weight = NULL;

void MLS_set_point_weights(const float *weights)
{
    g_mls_point_weight = weights;
}

/* Test hook (MLS_selftest, scripts/extract/mls_bench.c): 1 = the 2026-07
 * reference span kernel (float loads, tail mask on every batch, every sum with
 * two roundings, queries in index order).  Production never sets it. */
static int g_mls_reference_kernel = 0;

void MLS_test_reference_kernel(int on) { g_mls_reference_kernel = on; }

#if defined(__AVX2__)
/* r^2 = dz*dz + dy*dy + dx*dx of 4 lanes; with MLS_FMA_ACCUMULATE as two FMAs
 * (see pipeline_constants.h), else in the scalar order ((dz*dz + dy*dy) + dx*dx). */
static inline __m256d mls_r2_4(__m256d dz, __m256d dy, __m256d dx)
{
#if MLS_FMA_ACCUMULATE
    return _mm256_fmadd_pd(dx, dx, _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dz, dz)));
#else
    return _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dz, dz), _mm256_mul_pd(dy, dy)),
                         _mm256_mul_pd(dx, dx));
#endif
}

/* acc + a*b of 4 lanes (FMA under MLS_FMA_ACCUMULATE, else two roundings). */
static inline __m256d mls_madd4(__m256d a, __m256d b, __m256d acc)
{
#if MLS_FMA_ACCUMULATE
    return _mm256_fmadd_pd(a, b, acc);
#else
    return _mm256_add_pd(acc, _mm256_mul_pd(a, b));
#endif
}

/* Wendland weight of 4 lanes in wendland()'s op order: t = r/R, s = 1-t,
 * w = (s*s)*(s*s) * (4t + 1).  4t + 1 is ONE fmadd, which is exact: 4t is an
 * exact power-of-two scaling, so RN(4t + 1) == RN(RN(4t) + 1). */
static inline __m256d mls_wendland4(__m256d r2v, __m256d Rv)
{
    __m256d tv = _mm256_div_pd(_mm256_sqrt_pd(r2v), Rv);
    __m256d sv = _mm256_sub_pd(_mm256_set1_pd(1.0), tv);
    __m256d s2 = _mm256_mul_pd(sv, sv);
    return _mm256_mul_pd(_mm256_mul_pd(s2, s2),
                         _mm256_fmadd_pd(_mm256_set1_pd(4.0), tv, _mm256_set1_pd(1.0)));
}
#endif

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
        const char *r_env = getenv("MLS_RADIUS_VOX");
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
    assert(nv <= (size_t)INT_MAX);
    /* Scratch sizing below multiplies nv by small constants into `long`
     * ARENA_ALLOC sizes; a >1G-vert cloud is out of contract anyway. */
    assert(nv < ((size_t)1 << 30));

    Arena_Mark mark = Arena_save(arena);

    double R = (double)radius_vox;
    double R2 = R * R;

    float cell_size = radius_vox / (float)MLS_CELL_DIV;
    float inv_cell = (float)(1.0 / (double)cell_size);
    double cell_d = (double)cell_size;

    /* Optional global origin offset. Cells become world-aligned (rather
     * than cube-local-aligned) when both cubes pass the same offset
     * (e.g. cube origin in world coords). */
    float ocz = 0.0f, ocy = 0.0f, ocx = 0.0f;
    if (cell_origin) {
        ocz = cell_origin[0];
        ocy = cell_origin[1];
        ocx = cell_origin[2];
    }
    double oczd = (double)ocz, ocyd = (double)ocy, ocxd = (double)ocx;

    /* -------- Build the hash-CSR cell index --------
     * (a) hash-insert each vert's cell key, counting bucket sizes;
     * (b) prefix-sum to CSR starts;
     * (c) stable counting fill of SoA position snapshots + vert indices;
     * (d) per-bucket position sort (cross-cube-stable accumulation order).
     * All O(nv); replaces the old qsort + 54-binary-searches-per-vertex. */
    size_t T = 64;
    while (T < nv * 2) T <<= 1;

    uint64_t *tab_key  = (uint64_t *)ARENA_ALLOC(arena,
                            (size_t)T * sizeof(uint64_t));
    int32_t *tab_cnt   = (int32_t *)ARENA_ALLOC(arena,
                            (size_t)T * sizeof(int32_t));
    int32_t *tab_start = (int32_t *)ARENA_ALLOC(arena,
                            (size_t)T * sizeof(int32_t));
    uint32_t *slot_of  = (uint32_t *)ARENA_ALLOC(arena,
                            (size_t)nv * sizeof(uint32_t));
    uint32_t *uniq     = (uint32_t *)ARENA_ALLOC(arena,
                            (size_t)nv * sizeof(uint32_t));
    size_t n_uniq = 0;

    memset(tab_key, 0xFF, T * sizeof(uint64_t));   /* all MLS_HASH_EMPTY */
    memset(tab_cnt, 0, T * sizeof(int32_t));

    for (size_t i = 0; i < nv; i++) {
        float z = verts[i * 3 + 0];
        float y = verts[i * 3 + 1];
        float x = verts[i * 3 + 2];
        int32_t cz = (int32_t)floorf((z + ocz) * inv_cell);
        int32_t cy = (int32_t)floorf((y + ocy) * inv_cell);
        int32_t cx = (int32_t)floorf((x + ocx) * inv_cell);
        uint64_t key = cell_key(cz, cy, cx);
        size_t h = (size_t)(hash64(key) & (uint64_t)(T - 1));
        for (;;) {
            uint64_t kh = tab_key[h];
            if (kh == key) break;
            if (kh == MLS_HASH_EMPTY) {
                tab_key[h] = key;
                uniq[n_uniq++] = (uint32_t)h;
                break;
            }
            h = (h + 1) & (T - 1);
        }
        tab_cnt[h]++;
        slot_of[i] = (uint32_t)h;
    }

    {
        int32_t run = 0;
        for (size_t s = 0; s < T; s++) {
            tab_start[s] = run;
            run += tab_cnt[s];
        }
        assert((size_t)run == nv);
    }

    /* SoA snapshot of positions in CSR order, +3 floats of zeroed pad so
     * the AVX2 tail batch can issue a full 4-float load without reading
     * past the allocation (pad lanes are mask-excluded). The query loop
     * reads ONLY these (never `verts`), so out_verts may alias verts
     * safely. */
    float   *bz   = (float *)ARENA_ALLOC(arena,
                        (size_t)(nv + 3) * sizeof(float));
    float   *by   = (float *)ARENA_ALLOC(arena,
                        (size_t)(nv + 3) * sizeof(float));
    float   *bx   = (float *)ARENA_ALLOC(arena,
                        (size_t)(nv + 3) * sizeof(float));
    int32_t *bidx = (int32_t *)ARENA_ALLOC(arena,
                        (size_t)nv * sizeof(int32_t));
    bz[nv] = 0.0f; bz[nv+1] = 0.0f; bz[nv+2] = 0.0f;
    by[nv] = 0.0f; by[nv+1] = 0.0f; by[nv+2] = 0.0f;
    bx[nv] = 0.0f; bx[nv+1] = 0.0f; bx[nv+2] = 0.0f;
    {
        int32_t *cur = (int32_t *)ARENA_ALLOC(arena,
                            (size_t)T * sizeof(int32_t));
        memcpy(cur, tab_start, T * sizeof(int32_t));
        for (size_t i = 0; i < nv; i++) {
            size_t s = (size_t)slot_of[i];
            int32_t k = cur[s]++;
            bz[k] = verts[i * 3 + 0];
            by[k] = verts[i * 3 + 1];
            bx[k] = verts[i * 3 + 2];
            bidx[k] = (int32_t)i;
        }
    }
    for (size_t u = 0; u < n_uniq; u++) {
        size_t s = (size_t)uniq[u];
        if (tab_cnt[s] > 1) {
            bucket_sort_pos(bz, by, bx, bidx,
                            tab_start[s], tab_start[s] + tab_cnt[s]);
        }
    }

#if defined(__AVX2__)
    /* Double copies of the sorted SoA positions for the AVX2 kernel: it loads
     * them instead of converting three floats per batch (float -> double is
     * exact, so every value is unchanged).  Zero pads as for the floats. */
    double *bzd = (double *)ARENA_ALLOC(arena, (size_t)(nv + 3) * sizeof(double));
    double *byd = (double *)ARENA_ALLOC(arena, (size_t)(nv + 3) * sizeof(double));
    double *bxd = (double *)ARENA_ALLOC(arena, (size_t)(nv + 3) * sizeof(double));
    for (size_t q = 0; q < nv + 3; q++) {
        bzd[q] = (double)bz[q];
        byd[q] = (double)by[q];
        bxd[q] = (double)bx[q];
    }
#endif
    const int ref_kernel = g_mls_reference_kernel;

    /* Per occupied cell, its EXISTING neighbour cells in the query's
     * (dz, dy, dx) order, each with its CSR span.  Every vertex of a cell
     * walked the same 5^3 offsets and hashed each one; the list does those
     * lookups once per cell (a cell holds ~10^2 vertices).  The per-vertex
     * sphere/box cull below is unchanged, so the spans visited, their order
     * and every accumulation are identical (2026-10-04). */
    int32_t *uidx_of_slot = (int32_t *)ARENA_ALLOC(arena,
                                (size_t)T * sizeof(int32_t));
    for (size_t u = 0; u < n_uniq; u++) uidx_of_slot[uniq[u]] = (int32_t)u;
    int32_t *nb_off = (int32_t *)ARENA_ALLOC(arena,
                          (size_t)(n_uniq + 1) * sizeof(int32_t));
    MlsNb *nbl = NULL;
    for (int pass = 0; pass < 2; pass++) {
        size_t ne = 0;
        for (size_t u = 0; u < n_uniq; u++) {
            uint64_t key = tab_key[uniq[u]];
            int32_t ccz = cell_key_axis(key, 42);
            int32_t ccy = cell_key_axis(key, 21);
            int32_t ccx = cell_key_axis(key, 0);
            if (pass == 0) nb_off[u] = (int32_t)ne;
            for (int dz = -MLS_CELL_DIV; dz <= MLS_CELL_DIV; dz++)
                for (int dy = -MLS_CELL_DIV; dy <= MLS_CELL_DIV; dy++)
                    for (int dx = -MLS_CELL_DIV; dx <= MLS_CELL_DIV; dx++) {
                        size_t h = hash_find(tab_key, T,
                                             cell_key(ccz + dz, ccy + dy, ccx + dx));
                        if (h == SIZE_MAX) continue;
                        if (pass == 1) {
                            MlsNb *e = &nbl[ne];
                            e->start = tab_start[h];
                            e->cnt = tab_cnt[h];
                            e->dz = (int8_t)dz; e->dy = (int8_t)dy; e->dx = (int8_t)dx;
                            e->pad = 0;
                        }
                        ne++;
                    }
        }
        if (pass == 0) {
            nb_off[n_uniq] = (int32_t)ne;
            nbl = (MlsNb *)ARENA_ALLOC(arena, (size_t)((ne ? ne : 1) * sizeof(MlsNb)));
        }
    }

    /* -------- Per-vertex MLS midpoint (fused single sweep) -------- */
    /* Each vertex is independent -- all-local accumulators, read-only cell
     * index, one write to out_*[i] -- so this parallelizes exactly: output is
     * byte-identical to serial for any thread count/schedule (community
     * issue #3, pscamillo: serial MLS was ~83% of single-cube wall; 9.6x on
     * 16 threads). Signed int index for MSVC OpenMP 2.0 (a per-component
     * cloud is < 2^31 verts by orders of magnitude). Thread count follows
     * omp_set_num_threads(n_threads) in the driver. */
    int nv_i = (int)nv;
    int si = 0;   /* declared before the loop: MSVC OpenMP 2.0 form */
#pragma omp parallel for schedule(dynamic, 256)
    for (si = 0; si < nv_i; si++) {
        /* Queries in CSR (cell) order: consecutive queries share their
         * neighbour spans, which stay in cache.  Each query is independent and
         * writes only out_*[i], so the order changes no output. */
        size_t i = ref_kernel ? (size_t)si : (size_t)bidx[si];
        float vz = verts[i * 3 + 0];
        float vy = verts[i * 3 + 1];
        float vx = verts[i * 3 + 2];
        double vzd = (double)vz, vyd = (double)vy, vxd = (double)vx;

        int32_t ciz = (int32_t)floorf((vz + ocz) * inv_cell);
        int32_t ciy = (int32_t)floorf((vy + ocy) * inv_cell);
        int32_t cix = (int32_t)floorf((vx + ocx) * inv_cell);

        /* World coords for the sphere/box cull (cell boxes live in the
         * world-aligned frame the cell indices were computed in). */
        double wz = vzd + oczd;
        double wy = vyd + ocyd;
        double wx = vxd + ocxd;

        /* Weighted moments about v, accumulated in double. */
        double w_sum = 0.0;
        double sdz = 0.0, sdy = 0.0, sdx = 0.0;
        double mzz = 0.0, mzy = 0.0, mzx = 0.0;
        double myy = 0.0, myx = 0.0, mxx = 0.0;
        size_t n_used = 0;

#if defined(__AVX2__)
        /* 4-wide-double lane accumulators. Per-lane arithmetic replicates
         * the scalar op ORDER exactly (mul/add/div/sqrt only — no FMA), so
         * every element's contribution is bit-identical to the scalar
         * path's; only the accumulation GROUPING differs (4 lane partials
         * + a fixed left-assoc horizontal reduce at the end). Lane
         * assignment is a pure function of the position-sorted span
         * order, which is cross-cube stable, so the seam-convergence
         * contract is unchanged. Excluded/pad lanes are masked to w = 0
         * and contribute an exact +0.0 no-op to every accumulator. */
        __m256d aw  = _mm256_setzero_pd();
        __m256d az  = _mm256_setzero_pd();
        __m256d ay  = _mm256_setzero_pd();
        __m256d ax  = _mm256_setzero_pd();
        __m256d azz = _mm256_setzero_pd();
        __m256d azy = _mm256_setzero_pd();
        __m256d azx = _mm256_setzero_pd();
        __m256d ayy = _mm256_setzero_pd();
        __m256d ayx = _mm256_setzero_pd();
        __m256d axx = _mm256_setzero_pd();
        const __m256d vz4   = _mm256_set1_pd(vzd);
        const __m256d vy4   = _mm256_set1_pd(vyd);
        const __m256d vx4   = _mm256_set1_pd(vxd);
        const __m256d R2v   = _mm256_set1_pd(R2);
        const __m256d Rv    = _mm256_set1_pd(R);
        const __m256d one4  = _mm256_set1_pd(1.0);
        const __m256d four4 = _mm256_set1_pd(4.0);
        const __m256d zero4 = _mm256_setzero_pd();
#endif

        const int32_t uu = uidx_of_slot[slot_of[i]];
        for (int32_t e = nb_off[uu]; e < nb_off[uu + 1]; e++) {
                    const MlsNb *nb = &nbl[e];
                    const int dz = nb->dz, dy = nb->dy, dx = nb->dx;
                    double lo_z = (double)(ciz + dz) * cell_d;
                    double ddz = (wz < lo_z) ? (lo_z - wz)
                               : (wz > lo_z + cell_d) ? (wz - (lo_z + cell_d)) : 0.0;
                    double ddz2 = ddz * ddz;
                    if (ddz2 > R2 + MLS_CULL_MARGIN_SQ) continue;
                    double lo_y = (double)(ciy + dy) * cell_d;
                    double ddy = (wy < lo_y) ? (lo_y - wy)
                               : (wy > lo_y + cell_d) ? (wy - (lo_y + cell_d)) : 0.0;
                    double dzy2 = ddz2 + ddy * ddy;
                    if (dzy2 > R2 + MLS_CULL_MARGIN_SQ) continue;
                    double lo_x = (double)(cix + dx) * cell_d;
                    double ddx = (wx < lo_x) ? (lo_x - wx)
                               : (wx > lo_x + cell_d) ? (wx - (lo_x + cell_d)) : 0.0;
                    /* Conservative sphere/box cull: only skip when even the
                     * closest point of the cell box is beyond R (+ margin).
                     * The exact per-point r2 test below owns correctness. */
                    if (dzy2 + ddx * ddx > R2 + MLS_CULL_MARGIN_SQ) continue;

                    int32_t k0 = nb->start;
                    int32_t k1 = k0 + nb->cnt;
                    if (g_mls_point_weight != NULL) {
                    /* CT-weighted kernel (scalar path on every build): the
                     * Wendland weight is scaled by the neighbour's brightness
                     * weight, indexed through bidx (bucket -> original). */
                    for (int32_t k = k0; k < k1; k++) {
                        double dvz = (double)bz[k] - vzd;
                        double dvy = (double)by[k] - vyd;
                        double dvx = (double)bx[k] - vxd;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 >= R2) continue;
                        double r = sqrt(r2);
                        double w = wendland(r, R) *
                                   (double)g_mls_point_weight[bidx[k]];
                        if (w <= 0.0) continue;
                        w_sum += w;
                        double wdz = w * dvz;
                        double wdy = w * dvy;
                        double wdx = w * dvx;
                        sdz += wdz; sdy += wdy; sdx += wdx;
                        mzz += wdz * dvz;
                        mzy += wdz * dvy;
                        mzx += wdz * dvx;
                        myy += wdy * dvy;
                        myx += wdy * dvx;
                        mxx += wdx * dvx;
                        n_used++;
                    }
                    } else {

#if defined(__AVX2__)
                        if (!ref_kernel) {
                        /* 2026-10-06 kernel: double loads, the tail mask only on the
                         * span's last batch, FMA 4t+1 (exact), FMA r^2 and second
                         * moments under MLS_FMA_ACCUMULATE.  Same spans, lanes and
                         * accumulation order as the reference loop below. */
                        for (int32_t k = k0; k < k1; k += 4) {
                            __m256d dvz4 = _mm256_sub_pd(_mm256_loadu_pd(bzd + k), vz4);
                            __m256d dvy4 = _mm256_sub_pd(_mm256_loadu_pd(byd + k), vy4);
                            __m256d dvx4 = _mm256_sub_pd(_mm256_loadu_pd(bxd + k), vx4);
                            __m256d r2v = mls_r2_4(dvz4, dvy4, dvx4);
                            __m256d mdist = _mm256_cmp_pd(r2v, R2v, _CMP_LT_OQ);
                            if (k1 - k < 4)   /* the span's last batch: lanes past k1 are pads */
                                mdist = _mm256_and_pd(mdist, _mm256_castsi256_pd(
                                    _mm256_loadu_si256((const __m256i *)MLS_TAIL_MASK[k1 - k])));
                            if (_mm256_movemask_pd(mdist) == 0) continue;
                            __m256d wv = mls_wendland4(r2v, Rv);
                            __m256d mfull = _mm256_and_pd(mdist,
                                _mm256_cmp_pd(wv, zero4, _CMP_GT_OQ));
                            wv = _mm256_and_pd(wv, mfull);
                            n_used += (size_t)_mm_popcnt_u32(
                                (unsigned)_mm256_movemask_pd(mfull));
                            aw = _mm256_add_pd(aw, wv);
                            __m256d wdz = _mm256_mul_pd(wv, dvz4);
                            __m256d wdy = _mm256_mul_pd(wv, dvy4);
                            __m256d wdx = _mm256_mul_pd(wv, dvx4);
                            az = _mm256_add_pd(az, wdz);
                            ay = _mm256_add_pd(ay, wdy);
                            ax = _mm256_add_pd(ax, wdx);
                            azz = mls_madd4(wdz, dvz4, azz);
                            azy = mls_madd4(wdz, dvy4, azy);
                            azx = mls_madd4(wdz, dvx4, azx);
                            ayy = mls_madd4(wdy, dvy4, ayy);
                            ayx = mls_madd4(wdy, dvx4, ayx);
                            axx = mls_madd4(wdx, dvx4, axx);
                        }
                        } else
                    for (int32_t k = k0; k < k1; k += 4) {
                        int32_t rem = k1 - k;
                        __m256d valid = _mm256_castsi256_pd(
                            _mm256_loadu_si256((const __m256i *)
                                MLS_TAIL_MASK[(rem >= 4) ? 0 : rem]));
                        __m256d pz = _mm256_cvtps_pd(_mm_loadu_ps(bz + k));
                        __m256d py = _mm256_cvtps_pd(_mm_loadu_ps(by + k));
                        __m256d px = _mm256_cvtps_pd(_mm_loadu_ps(bx + k));
                        __m256d dvz4 = _mm256_sub_pd(pz, vz4);
                        __m256d dvy4 = _mm256_sub_pd(py, vy4);
                        __m256d dvx4 = _mm256_sub_pd(px, vx4);
                        /* r2 = ((dz*dz + dy*dy) + dx*dx) — scalar order */
                        __m256d r2v = _mm256_add_pd(
                            _mm256_add_pd(_mm256_mul_pd(dvz4, dvz4),
                                          _mm256_mul_pd(dvy4, dvy4)),
                            _mm256_mul_pd(dvx4, dvx4));
                        __m256d mdist = _mm256_and_pd(
                            _mm256_cmp_pd(r2v, R2v, _CMP_LT_OQ), valid);
                        if (_mm256_movemask_pd(mdist) == 0) continue;
                        /* Wendland, replicating wendland()'s op order:
                         * t = r/R; s = 1-t; w = (s*s)*(s*s) * (4*t + 1) */
                        __m256d rv = _mm256_sqrt_pd(r2v);
                        __m256d tv = _mm256_div_pd(rv, Rv);
                        __m256d sv = _mm256_sub_pd(one4, tv);
                        __m256d s2 = _mm256_mul_pd(sv, sv);
                        __m256d s4 = _mm256_mul_pd(s2, s2);
                        __m256d poly = _mm256_add_pd(
                            _mm256_mul_pd(four4, tv), one4);
                        __m256d wv = _mm256_mul_pd(s4, poly);
                        __m256d mfull = _mm256_and_pd(mdist,
                            _mm256_cmp_pd(wv, zero4, _CMP_GT_OQ));
                        wv = _mm256_and_pd(wv, mfull);
                        n_used += (size_t)_mm_popcnt_u32(
                            (unsigned)_mm256_movemask_pd(mfull));
                        aw = _mm256_add_pd(aw, wv);
                        __m256d wdz = _mm256_mul_pd(wv, dvz4);
                        __m256d wdy = _mm256_mul_pd(wv, dvy4);
                        __m256d wdx = _mm256_mul_pd(wv, dvx4);
                        az = _mm256_add_pd(az, wdz);
                        ay = _mm256_add_pd(ay, wdy);
                        ax = _mm256_add_pd(ax, wdx);
                        azz = _mm256_add_pd(azz, _mm256_mul_pd(wdz, dvz4));
                        azy = _mm256_add_pd(azy, _mm256_mul_pd(wdz, dvy4));
                        azx = _mm256_add_pd(azx, _mm256_mul_pd(wdz, dvx4));
                        ayy = _mm256_add_pd(ayy, _mm256_mul_pd(wdy, dvy4));
                        ayx = _mm256_add_pd(ayx, _mm256_mul_pd(wdy, dvx4));
                        axx = _mm256_add_pd(axx, _mm256_mul_pd(wdx, dvx4));
                    }
#else
                    for (int32_t k = k0; k < k1; k++) {
                        double dvz = (double)bz[k] - vzd;
                        double dvy = (double)by[k] - vyd;
                        double dvx = (double)bx[k] - vxd;
                        double r2 = dvz*dvz + dvy*dvy + dvx*dvx;
                        if (r2 >= R2) continue;
                        double r = sqrt(r2);
                        double w = wendland(r, R);
                        if (w <= 0.0) continue;
                        w_sum += w;
                        double wdz = w * dvz;
                        double wdy = w * dvy;
                        double wdx = w * dvx;
                        sdz += wdz; sdy += wdy; sdx += wdx;
                        mzz += wdz * dvz;
                        mzy += wdz * dvy;
                        mzx += wdz * dvx;
                        myy += wdy * dvy;
                        myx += wdy * dvx;
                        mxx += wdx * dvx;
                        n_used++;
                    }
#endif
                    }   /* weighted vs plain kernel */
        }

#if defined(__AVX2__)
        /* Fixed left-assoc lane reduce: ((l0+l1)+l2)+l3.  Skipped when the
         * CT-weighted scalar kernel ran: its sums live in the scalar
         * accumulators (the lanes are zero and would erase them). */
        if (g_mls_point_weight == NULL) {
            double l[4];
            _mm256_storeu_pd(l, aw);  w_sum = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, az);  sdz   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, ay);  sdy   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, ax);  sdx   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, azz); mzz   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, azy); mzy   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, azx); mzx   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, ayy); myy   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, ayx); myx   = ((l[0]+l[1])+l[2])+l[3];
            _mm256_storeu_pd(l, axx); mxx   = ((l[0]+l[1])+l[2])+l[3];
        }
#endif

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
        double pbz = sdz * inv_w;      /* centroid offset from v */
        double pby = sdy * inv_w;
        double pbx = sdx * inv_w;
        double cz = vzd + pbz;
        double cy = vyd + pby;
        double cx = vxd + pbx;

        if (out_normals) {
            /* Covariance about the centroid from the fused moments:
             * cov = M2/w_sum - pbar*pbar^T. Diagonal terms are >= 0
             * analytically; roundoff can leave a ~-1e-18 residue, which
             * the Jacobi solver is indifferent to. */
            double cov[3][3];
            cov[0][0] = mzz * inv_w - pbz * pbz;
            cov[0][1] = mzy * inv_w - pbz * pby; cov[1][0] = cov[0][1];
            cov[0][2] = mzx * inv_w - pbz * pbx; cov[2][0] = cov[0][2];
            cov[1][1] = myy * inv_w - pby * pby;
            cov[1][2] = myx * inv_w - pby * pbx; cov[2][1] = cov[1][2];
            cov[2][2] = mxx * inv_w - pbx * pbx;

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
                double dvz_c = vzd - cz;
                double dvy_c = vyd - cy;
                double dvx_c = vxd - cx;
                double d_proj = dvz_c * nz_d + dvy_c * ny_d + dvx_c * nx_d;
                out_verts[i*3+0] = (float)(vzd - d_proj * nz_d);
                out_verts[i*3+1] = (float)(vyd - d_proj * ny_d);
                out_verts[i*3+2] = (float)(vxd - d_proj * nx_d);
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
        const char *pe = getenv("MLS_PERTURB_EPS");
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

static uint64_t mls_test_rand(uint64_t *s)
{
    *s += 0x9E3779B97F4A7C15ULL;
    return hash64(*s);
}

/* Projection with the production and the reference kernel on one cloud: the
 * float positions and normals must match bit for bit (MLS_FMA_ACCUMULATE 0) or
 * within 1e-4 vox / 1e-5 (FMA: last-bit moment differences).  Returns 0 on pass. */
static int mls_selftest_cloud(Arena_T arena, const char *what, const float *v,
                              size_t nv)
{
    float *o1 = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    float *n1 = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    float *o2 = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    float *n2 = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    const float origin[3] = { 1024.0f, 2048.0f, 3072.0f };
    MLS_test_reference_kernel(0);
    MLS_project_verts(arena, v, nv, 12.0f, origin, o1, n1);
    MLS_test_reference_kernel(1);
    MLS_project_verts(arena, v, nv, 12.0f, origin, o2, n2);
    MLS_test_reference_kernel(0);
    size_t diff = 0, moved = 0;
    double dmax_v = 0.0, dmax_n = 0.0;
    for (size_t q = 0; q < nv * 3; q++) {
        moved += (o1[q] != v[q]);
        if (memcmp(&o1[q], &o2[q], sizeof(float)) != 0) {
            diff++;
            double d = fabs((double)o1[q] - (double)o2[q]);
            if (d > dmax_v) dmax_v = d;
        }
        if (memcmp(&n1[q], &n2[q], sizeof(float)) != 0) {
            diff++;
            double d = fabs((double)n1[q] - (double)n2[q]);
            if (d > dmax_n) dmax_n = d;
        }
    }
    int ok = moved > 0 && (MLS_FMA_ACCUMULATE ? (dmax_v <= 1e-4 && dmax_n <= 1e-5) : diff == 0);
    fprintf(stderr, "MLS_selftest: %s (%zu verts): kernel vs reference %zu float(s) differ "
            "(max %.3g vox, normal %.3g)%s\n", what, nv, diff, dmax_v, dmax_n,
            ok ? "" : " FAIL");
    return ok ? 0 : -1;
}

/* MLS selftest: (a) the exact 4t+1 fmadd equals the two-rounding form bit for
 * bit; (b) the production kernel against the reference kernel (cell-order
 * queries, double loads, last-batch tail mask, MLS_FMA_ACCUMULATE) on a noisy
 * curved double envelope and on a voxel-centre slab like the extract input. */
int MLS_selftest(void)
{
    int fails = 0;
#if defined(__AVX2__)
    {
        uint64_t s = 12345;
        size_t bad = 0, checked = 0;
        for (size_t n = 0; n < ((size_t)1 << 20); n++) {
            double t[4];
            for (int l = 0; l < 4; l++) {
                uint64_t h = mls_test_rand(&s);
                t[l] = (double)(h >> 11) * (1.0 / 9007199254740992.0);
                if ((h & 15) == 0) t[l] = ldexp(t[l], -(int)((h >> 4) & 63));
            }
            if (n == 0) { t[0] = 0.0; t[1] = 1.0; t[2] = 0.25; t[3] = nextafter(1.0, 0.0); }
            __m256d tv = _mm256_loadu_pd(t);
            double a[4], b[4];
            _mm256_storeu_pd(a, _mm256_fmadd_pd(_mm256_set1_pd(4.0), tv, _mm256_set1_pd(1.0)));
            _mm256_storeu_pd(b, _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(4.0), tv),
                                              _mm256_set1_pd(1.0)));
            bad += memcmp(a, b, sizeof(a)) != 0;
            checked += 4;
        }
        fprintf(stderr, "MLS_selftest: fmadd 4t+1 vs mul+add: %zu values, %zu mismatches\n",
                checked, bad);
        if (bad) fails++;
    }
#endif
    {
        Arena_T arena = Arena_new();
        uint64_t s = 777;
        const int nx = 48, ny = 40;
        size_t nv = (size_t)(2 * nx * ny), k = 0;
        float *v = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof(float));
        for (int layer = 0; layer < 2; layer++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++, k++) {
                    double jit = (double)(mls_test_rand(&s) >> 11) * (1.0 / 9007199254740992.0);
                    double x = 0.5 + (double)i, y = 0.5 + (double)j;
                    v[k * 3 + 0] = (float)(20.0 + 0.02 * (x - 24.0) * (x - 24.0) +
                                           3.0 * (double)layer + 0.4 * jit);
                    v[k * 3 + 1] = (float)y;
                    v[k * 3 + 2] = (float)x;
                }
        if (mls_selftest_cloud(arena, "noisy double envelope", v, nv) != 0) fails++;
        /* voxel centres of a 5-vox-thick curved sheet in a 40x48x48 box */
        size_t cap = 40 * 48 * 48, m = 0;
        float *w = (float *)ARENA_ALLOC(arena, cap * 3 * sizeof(float));
        for (int z = 0; z < 40; z++)
            for (int y = 0; y < 48; y++)
                for (int x = 0; x < 48; x++) {
                    double mid = 18.0 + 6.0 * sin(0.13 * (double)x) + 0.05 * (double)y;
                    if (fabs((double)z + 0.5 - mid) < 2.5) {
                        w[m * 3 + 0] = (float)z + 0.5f;
                        w[m * 3 + 1] = (float)y + 0.5f;
                        w[m * 3 + 2] = (float)x + 0.5f;
                        m++;
                    }
                }
        if (mls_selftest_cloud(arena, "voxel-centre sheet", w, m) != 0) fails++;
        Arena_dispose(&arena);
    }
    return fails ? -1 : 0;
}
