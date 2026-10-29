/* vcr_render.c -- see vcr_render.h. */
#include "flatten/vcr_render.h"

#include <math.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../common/ves_omp.h"

#include "common/arena.h"
#include "common/cuda_drv.h"
#include "common/json_read.h"
#include "common/tiff_io.h"
#include "common/ves_platform.h"
#include "common/zarr_shard.h"

#define VCR_CHUNK      128
#define VCR_CHUNK3     ((size_t)VCR_CHUNK * VCR_CHUNK * VCR_CHUNK)
#define VCR_TILE       ((size_t)VCR_CHUNK * VCR_CHUNK)      /* one output slice of a [N, 128, 128] chunk */
#define VCR_LOAD_BATCH 256                                   /* volume chunks read per upload (512 MB) */
#define VCR_BLOCK_X    128
#define VCR_BLOCK_Y    4

/* ---------------------------------------------------------------------------
 * The kernels.  Every arithmetic step is an explicit round-to-nearest
 * intrinsic in villa's order (and the module is compiled with --fmad=false),
 * so no multiply-add is ever fused except where villa calls std::fma.
 * ------------------------------------------------------------------------- */
static const char *KERNELS =
"typedef unsigned char u8;\n"
"\n"
"__device__ __forceinline__ void node_at(const float *X, const float *Y, const float *Z, int W, int r, int c,\n"
"                                        float *p)\n"
"{\n"
"    const long long i = (long long)r * W + c;\n"
"    p[0] = X[i]; p[1] = Y[i]; p[2] = Z[i];\n"
"    if (p[0] == -1.0f) { p[0] = __int_as_float(0x7fc00000); p[1] = p[0]; p[2] = p[0]; }\n"
"}\n"
"\n"
"__device__ __forceinline__ float bil(float p00, float p01, float p10, float p11)\n"
"{\n"
"    const float top = __fadd_rn(__fmul_rn(p00, 0.5f), __fmul_rn(p01, 0.5f));\n"
"    const float bot = __fadd_rn(__fmul_rn(p10, 0.5f), __fmul_rn(p11, 0.5f));\n"
"    return __fadd_rn(__fmul_rn(top, 0.5f), __fmul_rn(bot, 0.5f));\n"
"}\n"
"\n"
"/* canvas pixel (r, c): 0 = no samples; 1 = base b (maybe NaN) and unit direction d */\n"
"__device__ int pixel_geom(const float *X, const float *Y, const float *Z, int W, int H, int r, int c, int flip,\n"
"                          float *b, float *d)\n"
"{\n"
"    const int rn = r + 1, cn = c + 1;\n"
"    float q[3], l[3], rt[3], u[3], w[3], a[3], e[3], f[3], g[3];\n"
"    if (rn >= H || cn >= W) return 0;\n"
"    node_at(X, Y, Z, W, rn, cn, q);\n"
"    if ((q[0] != q[0]) | (q[1] != q[1]) | (q[2] != q[2]) | ((q[0] == -1.f) & (q[1] == -1.f) & (q[2] == -1.f)))\n"
"        return 0;\n"
"    if (rn > H - 2 || cn > W - 2) return 0;\n"
"    node_at(X, Y, Z, W, rn, cn - 1, l);\n"
"    node_at(X, Y, Z, W, rn, cn + 1, rt);\n"
"    node_at(X, Y, Z, W, rn - 1, cn, u);\n"
"    node_at(X, Y, Z, W, rn + 1, cn, w);\n"
"    const float xv0 = __fsub_rn(rt[0], l[0]), xv1 = __fsub_rn(rt[1], l[1]), xv2 = __fsub_rn(rt[2], l[2]);\n"
"    const float yv0 = __fsub_rn(w[0], u[0]), yv1 = __fsub_rn(w[1], u[1]), yv2 = __fsub_rn(w[2], u[2]);\n"
"    float n0 = __fsub_rn(__fmul_rn(xv1, yv2), __fmul_rn(xv2, yv1));\n"
"    float n1 = __fsub_rn(__fmul_rn(xv2, yv0), __fmul_rn(xv0, yv2));\n"
"    float n2 = __fsub_rn(__fmul_rn(xv0, yv1), __fmul_rn(xv1, yv0));\n"
"    const float len2 = __fadd_rn(__fadd_rn(__fmul_rn(n0, n0), __fmul_rn(n1, n1)), __fmul_rn(n2, n2));\n"
"    if (len2 == 0.0f || len2 != len2) return 0;\n"
"    const float inv = __fdiv_rn(1.0f, __fsqrt_rn(len2));\n"
"    n0 = __fmul_rn(n0, inv); n1 = __fmul_rn(n1, inv); n2 = __fmul_rn(n2, inv);\n"
"    if (flip) { n0 = -n0; n1 = -n1; n2 = -n2; }\n"
"    const float L2 = __fadd_rn(__fadd_rn(__fmul_rn(n0, n0), __fmul_rn(n1, n1)), __fmul_rn(n2, n2));\n"
"    if (L2 > 0.0f) {\n"
"        const float ia = __fdiv_rn(1.0f, __fsqrt_rn(L2));\n"
"        n0 = __fmul_rn(n0, ia); n1 = __fmul_rn(n1, ia); n2 = __fmul_rn(n2, ia);\n"
"    }\n"
"    d[0] = n0; d[1] = n1; d[2] = n2;\n"
"    node_at(X, Y, Z, W, r, c, a);\n"
"    node_at(X, Y, Z, W, r, c + 1, e);\n"
"    node_at(X, Y, Z, W, r + 1, c, f);\n"
"    node_at(X, Y, Z, W, r + 1, c + 1, g);\n"
"    for (int k = 0; k < 3; k++) b[k] = bil(a[k], e[k], f[k], g[k]);\n"
"    return 1;\n"
"}\n"
"\n"
"__device__ __forceinline__ float vox(const int *tab, const u8 *pool, int sz, int sy, int sx, int cgy, int cgx,\n"
"                                     int iz, int iy, int ix)\n"
"{\n"
"    if ((unsigned)iz >= (unsigned)sz || (unsigned)iy >= (unsigned)sy || (unsigned)ix >= (unsigned)sx) return 0.0f;\n"
"    const int slot = tab[((iz >> 7) * cgy + (iy >> 7)) * cgx + (ix >> 7)];\n"
"    if (slot < 0) return 0.0f;\n"
"    return (float)pool[(size_t)slot * 2097152u + (((size_t)(iz & 127) << 14) | ((size_t)(iy & 127) << 7)\n"
"                                                 | (size_t)(ix & 127))];\n"
"}\n"
"\n"
"__device__ float trilinear(const int *tab, const u8 *pool, int sz, int sy, int sx, int cgy, int cgx,\n"
"                           float vz, float vy, float vx)\n"
"{\n"
"    const int iz = (int)vz, iy = (int)vy, ix = (int)vx;\n"
"    const float c000 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz, iy, ix);\n"
"    const float c100 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz + 1, iy, ix);\n"
"    const float c010 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz, iy + 1, ix);\n"
"    const float c110 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz + 1, iy + 1, ix);\n"
"    const float c001 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz, iy, ix + 1);\n"
"    const float c101 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz + 1, iy, ix + 1);\n"
"    const float c011 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz, iy + 1, ix + 1);\n"
"    const float c111 = vox(tab, pool, sz, sy, sx, cgy, cgx, iz + 1, iy + 1, ix + 1);\n"
"    const float fz = __fsub_rn(vz, (float)iz), fy = __fsub_rn(vy, (float)iy), fx = __fsub_rn(vx, (float)ix);\n"
"    const float c00 = __fmaf_rn(fx, __fsub_rn(c001, c000), c000);\n"
"    const float c01 = __fmaf_rn(fx, __fsub_rn(c011, c010), c010);\n"
"    const float c10 = __fmaf_rn(fx, __fsub_rn(c101, c100), c100);\n"
"    const float c11 = __fmaf_rn(fx, __fsub_rn(c111, c110), c110);\n"
"    const float c0 = __fmaf_rn(fy, __fsub_rn(c01, c00), c00);\n"
"    const float c1 = __fmaf_rn(fy, __fsub_rn(c11, c10), c10);\n"
"    return __fmaf_rn(fz, __fsub_rn(c1, c0), c0);\n"
"}\n"
"\n"
"extern \"C\" __global__ void vcr_mark(const float *X, const float *Y, const float *Z, int W, int H, int flip,\n"
"                                    float omin, float omax, int sz, int sy, int sx, int cgy, int cgx, u8 *mask)\n"
"{\n"
"    const int c = blockIdx.x * blockDim.x + threadIdx.x, r = blockIdx.y * blockDim.y + threadIdx.y;\n"
"    if (c >= W || r >= H) return;\n"
"    float b[3], d[3];\n"
"    if (!pixel_geom(X, Y, Z, W, H, r, c, flip, b, d)) return;\n"
"    if (!((b[0] == b[0]) & (b[1] == b[1]) & (b[2] == b[2]))) return;\n"
"    const int lim[3] = {sx, sy, sz};\n"
"    int lo[3], hi[3];\n"
"    for (int k = 0; k < 3; k++) {\n"
"        const float e1 = __fadd_rn(b[k], __fmul_rn(d[k], omin)), e2 = __fadd_rn(b[k], __fmul_rn(d[k], omax));\n"
"        const float mn = fminf(e1, e2), mx = fmaxf(e1, e2), top = (float)lim[k] + 4.0f;\n"
"        if (!(mx >= -4.0f) || !(mn <= top)) return;\n"
"        lo[k] = max(0, (int)floorf(fmaxf(mn, -4.0f)) - 1);\n"
"        hi[k] = min(lim[k] - 1, (int)floorf(fminf(mx, top)) + 2);\n"
"        if (lo[k] > hi[k]) return;\n"
"    }\n"
"    for (int cz = lo[2] >> 7; cz <= hi[2] >> 7; cz++)\n"
"        for (int cy = lo[1] >> 7; cy <= hi[1] >> 7; cy++)\n"
"            for (int cx = lo[0] >> 7; cx <= hi[0] >> 7; cx++)\n"
"                mask[(cz * cgy + cy) * cgx + cx] = 1;\n"
"}\n"
"\n"
"/* one 128-row band of the layer stack, written as zarr chunks: out[tx][slice][row][col] */\n"
"extern \"C\" __global__ void vcr_layers(const float *X, const float *Y, const float *Z, int W, int H, int flip,\n"
"                                      int r0, const int *tab, const u8 *pool, int sz, int sy, int sx, int cgy,\n"
"                                      int cgx, const float *offs, int N, u8 *out)\n"
"{\n"
"    const int c = blockIdx.x * blockDim.x + threadIdx.x, ly = blockIdx.y * blockDim.y + threadIdx.y;\n"
"    const int r = r0 + ly;\n"
"    if (c >= W || ly >= 128 || r >= H) return;\n"
"    float b[3], d[3];\n"
"    const int ok = pixel_geom(X, Y, Z, W, H, r, c, flip, b, d);\n"
"    const float szf = (float)sz, syf = (float)sy, sxf = (float)sx;\n"
"    u8 *o = out + (size_t)(c >> 7) * (size_t)N * 16384u + (size_t)ly * 128u + (size_t)(c & 127);\n"
"    for (int i = 0; i < N; i++) {\n"
"        u8 val = 0;\n"
"        if (ok) {\n"
"            const float off = offs[i];\n"
"            const float px = __fadd_rn(b[0], __fmul_rn(d[0], off));\n"
"            const float py = __fadd_rn(b[1], __fmul_rn(d[1], off));\n"
"            const float pz = __fadd_rn(b[2], __fmul_rn(d[2], off));\n"
"            if (pz >= 0.f && py >= 0.f && px >= 0.f && pz < szf && py < syf && px < sxf) {\n"
"                float v = trilinear(tab, pool, sz, sy, sx, cgy, cgx, pz, py, px);\n"
"                if (v < 0.f) v = 0.f;\n"
"                if (v > 255.f) v = 255.f;\n"
"                val = (u8)(int)v;\n"
"            }\n"
"        }\n"
"        o[(size_t)i * 16384u] = val;\n"
"    }\n"
"}\n"
"\n"
"extern \"C\" __global__ void vcr_ctmax(const float *X, const float *Y, const float *Z, int W, int H, int flip,\n"
"                                     int z0, int z1, float zstep, const int *tab, const u8 *pool, int sz, int sy,\n"
"                                     int sx, int cgy, int cgx, u8 *out)\n"
"{\n"
"    const int c = blockIdx.x * blockDim.x + threadIdx.x, r = blockIdx.y * blockDim.y + threadIdx.y;\n"
"    if (c >= W || r >= H) return;\n"
"    float b[3], d[3];\n"
"    u8 val = 0;\n"
"    if (pixel_geom(X, Y, Z, W, H, r, c, flip, b, d) && isfinite(b[0])) {\n"
"        const float szf = (float)sz, syf = (float)sy, sxf = (float)sx;\n"
"        float mx = 0.0f;\n"
"        for (int k = z0; k <= z1; k++) {\n"
"            const float z = __fmul_rn((float)k, zstep);\n"
"            const float px = __fadd_rn(b[0], __fmul_rn(d[0], z));\n"
"            const float py = __fadd_rn(b[1], __fmul_rn(d[1], z));\n"
"            const float pz = __fadd_rn(b[2], __fmul_rn(d[2], z));\n"
"            float v = 0.0f;\n"
"            if (pz >= 0.f && py >= 0.f && px >= 0.f && pz < szf && py < syf && px < sxf) {\n"
"                int iz = (int)__fadd_rn(pz, 0.5f), iy = (int)__fadd_rn(py, 0.5f), ix = (int)__fadd_rn(px, 0.5f);\n"
"                if (iz >= sz) iz = sz - 1;\n"
"                if (iy >= sy) iy = sy - 1;\n"
"                if (ix >= sx) ix = sx - 1;\n"
"                v = vox(tab, pool, sz, sy, sx, cgy, cgx, iz, iy, ix);\n"
"            }\n"
"            mx = fmaxf(mx, v);\n"
"        }\n"
"        if (mx > 255.f) mx = 255.f;\n"
"        val = (u8)(int)mx;\n"
"    }\n"
"    out[(size_t)r * W + c] = val;\n"
"}\n"
#include "vcr_render_cache_kernels.inc"
;

static const char *KERNEL_NAMES[] = {"vcr_mark", "vcr_layers", "vcr_ctmax",
    "vcr_mark_rect", "vcr_layers_rect", "vcr_ctmax_rect", "vcr_patch_table"};
enum { K_MARK = 0, K_LAYERS = 1, K_CTMAX = 2, K_MARK_RECT = 3, K_LAYERS_RECT = 4, K_CTMAX_RECT = 5, K_PATCH = 6 };

/* ---------------------------------------------------------------------------
 * A render scene: the node grid on the host and a source of volume chunks.
 * ------------------------------------------------------------------------- */
typedef struct {
    int          W, H;                       /* node grid = canvas (--scale 1 of a scale-1 tifxyz) */
    const float *X, *Y, *Z;
    /* Optional bounded coordinate source; in-memory planes are test fixtures. */
    int        (*nodes)(void *ctx, int row, int col, int height, int width, float *xyz, size_t plane_stride);
    void        *node_ctx;
    long         shape[3];                   /* volume z, y, x */
    long         cg[3];                      /* chunk grid z, y, x */
    /* dst NULL: is the chunk stored?  else read it into dst (128^3).  1 stored, 0 absent, -1 error */
    int        (*load)(void *ctx, long cz, long cy, long cx, uint8_t *dst);
    void        *ctx;
} Scene;

typedef int (*BandFn)(void *ctx, int ty, int tx0, const uint8_t *band, int ntx, int N);

static double wall_s(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}

/* villa buildOffsetList: center = 0.5 * (max(1, N) - 1), offset = float((i - center) * step) */
static void layer_offsets(int n, double step, float *out)
{
    const double center = 0.5 * ((n > 1 ? n : 1) - 1.0);
    for (int i = 0; i < n; i++) out[i] = (float)(((double)i - center) * step);
}

static int gpu_render_resident(CudaDrv cu, const Scene *s, const VcrRender_job *job, BandFn band_fn, void *band_ctx,
                      uint8_t *comp_out, char *why, size_t why_len)
{
    const size_t nodes = (size_t)s->W * (size_t)s->H;
    const size_t ncg = (size_t)(s->cg[0] * s->cg[1] * s->cg[2]);
    const int N = job->composite ? 0 : job->num_slices;
    const float zstep = (float)job->slice_step;
    CudaDrv_ptr dX = 0, dY = 0, dZ = 0, dmask = 0, dtab = 0, dpool = 0, doffs = 0, dout = 0;
    uint8_t *mask = NULL, *stage = NULL, *band = NULL;
    int *tab = NULL;
    long *keys = NULL;
    float offs[512];
    int rc = -1, flip = job->flip_normals ? 1 : 0;
    int W = s->W, H = s->H, sz = (int)s->shape[0], sy = (int)s->shape[1], sx = (int)s->shape[2];
    int cgy = (int)s->cg[1], cgx = (int)s->cg[2];
    float omin = 0.0f, omax = 0.0f;
    size_t n_marked = 0, n_present = 0;
    double t0 = wall_s();

    if (!job->composite) {
        if (N < 1 || N > 512) { snprintf(why, why_len, "num-slices %d outside 1..512", N); return -1; }
        layer_offsets(N, job->slice_step, offs);
        omin = offs[0];
        omax = offs[0];
        for (int i = 1; i < N; i++) { omin = offs[i] < omin ? offs[i] : omin; omax = offs[i] > omax ? offs[i] : omax; }
    } else {
        const float a = (float)job->comp_start * zstep, b = (float)job->comp_end * zstep;
        omin = a < b ? a : b;
        omax = a < b ? b : a;
    }
    if (CudaDrv_alloc(cu, &dX, nodes * 4) || CudaDrv_alloc(cu, &dY, nodes * 4) || CudaDrv_alloc(cu, &dZ, nodes * 4)
        || CudaDrv_upload(cu, dX, s->X, nodes * 4) || CudaDrv_upload(cu, dY, s->Y, nodes * 4)
        || CudaDrv_upload(cu, dZ, s->Z, nodes * 4) || CudaDrv_alloc(cu, &dmask, ncg)
        || CudaDrv_memset8(cu, dmask, 0, ncg))
        goto gpu_fail;
    {
        void *p[] = {&dX, &dY, &dZ, &W, &H, &flip, &omin, &omax, &sz, &sy, &sx, &cgy, &cgx, &dmask};
        if (CudaDrv_launch(cu, K_MARK, (unsigned)((W + VCR_BLOCK_X - 1) / VCR_BLOCK_X),
                           (unsigned)((H + VCR_BLOCK_Y - 1) / VCR_BLOCK_Y), VCR_BLOCK_X, VCR_BLOCK_Y, p))
            goto gpu_fail;
    }
    mask = (uint8_t *)malloc(ncg);
    tab = (int *)malloc(ncg * sizeof(int));
    keys = (long *)malloc(ncg * sizeof(long));
    if (!mask || !tab || !keys) { snprintf(why, why_len, "out of host memory (chunk tables)"); goto done; }
    if (CudaDrv_download(cu, mask, dmask, ncg)) goto gpu_fail;
    for (size_t i = 0; i < ncg; i++) {
        tab[i] = -1;
        if (mask[i]) keys[n_marked++] = (long)i;
    }
    {   /* which marked chunks are stored */
        int err = 0;
        long li = 0;
        #pragma omp parallel for schedule(dynamic, 64)
        for (li = 0; li < (long)n_marked; li++) {
            const long k = keys[li], cz = k / (s->cg[1] * s->cg[2]), cy = (k / s->cg[2]) % s->cg[1], cx = k % s->cg[2];
            const int st = s->load(s->ctx, cz, cy, cx, NULL);
            if (st < 0) err = 1;
            mask[k] = (uint8_t)(st == 1 ? 2 : 1);
        }
        if (err) { snprintf(why, why_len, "volume chunk listing failed"); goto done; }
        n_marked = 0;
        for (size_t i = 0; i < ncg; i++)
            if (mask[i] == 2) { tab[i] = (int)n_present++; keys[n_marked++] = (long)i; }
    }
    printf("[vcr] marked chunks: %zu stored (%.1f GB on the card) of the %zu-chunk grid (%.1fs)\n", n_present,
           (double)n_present * VCR_CHUNK3 / 1e9, ncg, wall_s() - t0);
    fflush(stdout);
    if (CudaDrv_alloc(cu, &dpool, n_present ? n_present * VCR_CHUNK3 : 1) || CudaDrv_alloc(cu, &dtab, ncg * sizeof(int))
        || CudaDrv_upload(cu, dtab, tab, ncg * sizeof(int)))
        goto gpu_fail;
    stage = (uint8_t *)malloc((size_t)VCR_LOAD_BATCH * VCR_CHUNK3);
    if (!stage) { snprintf(why, why_len, "out of host memory (load stage)"); goto done; }
    for (size_t b0 = 0; b0 < n_present; b0 += VCR_LOAD_BATCH) {
        const size_t nb = n_present - b0 < VCR_LOAD_BATCH ? n_present - b0 : VCR_LOAD_BATCH;
        int err = 0;
        long i = 0;
        #pragma omp parallel for schedule(dynamic, 1)
        for (i = 0; i < (long)nb; i++) {
            const long k = keys[b0 + (size_t)i];
            const long cz = k / (s->cg[1] * s->cg[2]), cy = (k / s->cg[2]) % s->cg[1], cx = k % s->cg[2];
            if (s->load(s->ctx, cz, cy, cx, stage + (size_t)i * VCR_CHUNK3) != 1) err = 1;
        }
        if (err) { snprintf(why, why_len, "volume chunk read failed"); goto done; }
        if (CudaDrv_upload(cu, dpool + (CudaDrv_ptr)(b0 * VCR_CHUNK3), stage, nb * VCR_CHUNK3)) goto gpu_fail;
    }
    free(stage);
    stage = NULL;
    printf("[vcr] volume on the card (%.1fs)\n", wall_s() - t0);
    fflush(stdout);

    if (job->composite) {
        int z0 = job->comp_start, z1 = job->comp_end;
        if (CudaDrv_alloc(cu, &dout, nodes) || CudaDrv_memset8(cu, dout, 0, nodes)) goto gpu_fail;
        void *p[] = {&dX, &dY, &dZ, &W, &H, &flip, &z0, &z1, (void *)&zstep, &dtab, &dpool, &sz, &sy, &sx, &cgy, &cgx,
                     &dout};
        if (CudaDrv_launch(cu, K_CTMAX, (unsigned)((W + VCR_BLOCK_X - 1) / VCR_BLOCK_X),
                           (unsigned)((H + VCR_BLOCK_Y - 1) / VCR_BLOCK_Y), VCR_BLOCK_X, VCR_BLOCK_Y, p)
            || CudaDrv_download(cu, comp_out, dout, nodes))
            goto gpu_fail;
    } else {
        const int ntx = (W + VCR_CHUNK - 1) / VCR_CHUNK, nty = (H + VCR_CHUNK - 1) / VCR_CHUNK;
        const size_t band_bytes = (size_t)ntx * (size_t)N * VCR_TILE;
        band = (uint8_t *)malloc(band_bytes);
        if (!band) { snprintf(why, why_len, "out of host memory (band)"); goto done; }
        if (CudaDrv_alloc(cu, &doffs, (size_t)N * sizeof(float)) || CudaDrv_upload(cu, doffs, offs, (size_t)N * sizeof(float))
            || CudaDrv_alloc(cu, &dout, band_bytes))
            goto gpu_fail;
        double tk = 0.0, tw = 0.0;
        for (int ty = 0; ty < nty; ty++) {
            int r0 = ty * VCR_CHUNK;
            void *p[] = {&dX, &dY, &dZ, &W, &H, &flip, &r0, &dtab, &dpool, &sz, &sy, &sx, &cgy, &cgx, &doffs,
                         (void *)&N, &dout};
            double ta = wall_s();
            if (CudaDrv_memset8(cu, dout, 0, band_bytes)
                || CudaDrv_launch(cu, K_LAYERS, (unsigned)ntx, VCR_CHUNK / VCR_BLOCK_Y, VCR_BLOCK_X, VCR_BLOCK_Y, p)
                || CudaDrv_download(cu, band, dout, band_bytes))
                goto gpu_fail;
            double tb = wall_s();
            if (band_fn(band_ctx, ty, 0, band, ntx, N)) { snprintf(why, why_len, "writing band %d failed", ty); goto done; }
            tk += tb - ta;
            tw += wall_s() - tb;
            if (ty == nty - 1 || (ty + 1) % 8 == 0) {
                printf("[vcr] band %d/%d  render %.1fs  write %.1fs\n", ty + 1, nty, tk, tw);
                fflush(stdout);
            }
        }
    }
    rc = 0;
    goto done;
gpu_fail:
    snprintf(why, why_len, "%s", CudaDrv_error(cu));
done:
    CudaDrv_free(cu, dX);
    CudaDrv_free(cu, dY);
    CudaDrv_free(cu, dZ);
    CudaDrv_free(cu, dmask);
    CudaDrv_free(cu, dtab);
    CudaDrv_free(cu, dpool);
    CudaDrv_free(cu, doffs);
    CudaDrv_free(cu, dout);
    free(mask);
    free(tab);
    free(keys);
    free(stage);
    free(band);
    return rc;
}

#include "vcr_render_cache.inc"

/* ---------------------------------------------------------------------------
 * Files: the volume's zarr array, the tifxyz, the outputs.
 * ------------------------------------------------------------------------- */
typedef struct {
    char path[1024];
    char sep;
} VolFiles;

static int vol_chunk(void *ctx, long cz, long cy, long cx, uint8_t *dst)
{
    const VolFiles *v = (const VolFiles *)ctx;
    char p[1200];
    if (v->sep == '/') snprintf(p, sizeof p, "%s/%ld/%ld/%ld", v->path, cz, cy, cx);
    else snprintf(p, sizeof p, "%s/%ld.%ld.%ld", v->path, cz, cy, cx);
    FILE *f = fopen(p, "rb");
    if (!f) return 0;
    if (!dst) { fclose(f); return 1; }
    size_t got = fread(dst, 1, VCR_CHUNK3, f);
    fclose(f);
    return got == VCR_CHUNK3 ? 1 : -1;
}

static int open_volume(Arena_T a, const char *dir, VolFiles *vf, long shape[3], char *why, size_t n)
{
    char p[1100];
    const char *err = NULL;
    snprintf(p, sizeof p, "%s/.zarray", dir);
    const JsonValue *j = Json_parse_file(a, p, &err);
    if (!j) { snprintf(why, n, "%s: %s", p, err ? err : "unreadable"); return -1; }
    const JsonValue *sh = Json_object_get(j, "shape"), *ch = Json_object_get(j, "chunks");
    const JsonValue *comp = Json_object_get(j, "compressor"), *flt = Json_object_get(j, "filters");
    const char *dt = Json_as_string(Json_object_get(j, "dtype"));
    const char *order = Json_as_string(Json_object_get(j, "order"));
    const char *sep = Json_as_string(Json_object_get(j, "dimension_separator"));
    if (Json_array_len(sh) != 3 || Json_array_len(ch) != 3 || !dt || strcmp(dt, "|u1") != 0
        || (comp && Json_type(comp) != JSON_NULL) || (flt && Json_type(flt) != JSON_NULL)
        || (order && strcmp(order, "C") != 0)) {
        snprintf(why, n, "%s: not an uncompressed C-order uint8 3-D zarr array", p);
        return -1;
    }
    for (int k = 0; k < 3; k++) {
        shape[k] = Json_as_long(Json_array_get(sh, (size_t)k), -1);
        if (Json_as_long(Json_array_get(ch, (size_t)k), -1) != VCR_CHUNK || shape[k] <= 0 || shape[k] > 0x7fffffffL) {
            snprintf(why, n, "%s: chunks must be 128^3 (and the shape positive)", p);
            return -1;
        }
    }
    snprintf(vf->path, sizeof vf->path, "%s", dir);
    vf->sep = (sep && strcmp(sep, "/") == 0) ? '/' : '.';
    return 0;
}

/* The pipeline's TIFF cube grid as the volume (user 2026-10-17: "make it run on the tiff format"): a folder of
 * z#####_y#####_x#####.tif cubes (a grid's cubes_RAW), each 128 pages of 128 x 128 uint8, page = z.  grid_carve copies
 * the zarr's chunks into them, so cube z{128 cz}_y{128 cy}_x{128 cx} holds exactly the bytes of zarr chunk
 * (cz, cy, cx) and a render reads the same values from either.  The ids are voxel origins on the 128 grid (any other
 * id is refused); the volume's shape is the extent corner of the furthest cube; a cube that is not there reads as fill
 * 0, as a missing zarr chunk does.  The loader runs on OpenMP threads: each decodes into its own scratch arena. */
#define VCR_CUBE_SCRATCH ((size_t)4 << 20)       /* one 2 MB cube */

typedef struct {
    char     dir[1024];
    int      nthreads;
    Arena_T *scratch;                            /* [nthreads]; slot t is made and used only by thread t */
} CubeFiles;

static int thread_slot(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

static int cube_chunk(void *ctx, long cz, long cy, long cx, uint8_t *dst)
{
    CubeFiles *v = (CubeFiles *)ctx;
    char p[1200];
    snprintf(p, sizeof p, "%s/z%05ld_y%05ld_x%05ld.tif", v->dir, cz * VCR_CHUNK, cy * VCR_CHUNK, cx * VCR_CHUNK);
    FILE *f = fopen(p, "rb");
    if (!f) return 0;
    fclose(f);
    if (!dst) return 1;
    const int t = thread_slot();
    Arena_T own = NULL, a = NULL;
    if (t >= 0 && t < v->nthreads) {
        if (!v->scratch[t]) v->scratch[t] = Arena_new_sized(VCR_CUBE_SCRATCH);
        a = v->scratch[t];
    } else {
        a = own = Arena_new_sized(VCR_CUBE_SCRATCH);
    }
    const Arena_Mark m = Arena_save(a);
    uint8_t *vol = NULL;
    int D = 0, H = 0, W = 0;
    const int ok = TiffIO_load(a, p, &vol, &D, &H, &W) == 0 && D == VCR_CHUNK && H == VCR_CHUNK && W == VCR_CHUNK;
    if (ok) memcpy(dst, vol, VCR_CHUNK3);
    else fprintf(stderr, "vcr_gpu: %s: not a 128-page 128 x 128 uint8 cube\n", p);
    Arena_restore(a, m);
    if (own) Arena_dispose(&own);
    return ok ? 1 : -1;
}

/* a cube file name -> its voxel origin (z, y, x); 0 for any other name */
static int cube_name(const char *name, long o[3])
{
    long z = -1, y = -1, x = -1;
    char tail[8] = "";
    if (sscanf(name, "z%ld_y%ld_x%ld%7s", &z, &y, &x, tail) != 4 || strcmp(tail, ".tif") != 0) return 0;
    o[0] = z;
    o[1] = y;
    o[2] = x;
    return 1;
}

static void cube_note(const long o[3], long extent[3], long *count, long *bad)
{
    for (int k = 0; k < 3; k++)
        if (o[k] < 0 || o[k] % VCR_CHUNK != 0) {
            (*bad)++;
            return;
        }
    for (int k = 0; k < 3; k++)
        if (o[k] + VCR_CHUNK > extent[k]) extent[k] = o[k] + VCR_CHUNK;
    (*count)++;
}

static int open_cubes(const char *dir, CubeFiles *cf, long shape[3], char *why, size_t n)
{
    long extent[3] = {0, 0, 0}, count = 0, bad = 0, o[3] = {0, 0, 0};
#ifdef _WIN32
    char glob[1100];
    snprintf(glob, sizeof glob, "%s/*.tif", dir);
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(glob, &fd);
    if (h != INVALID_HANDLE_VALUE) {
        do {
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && cube_name(fd.cFileName, o))
                cube_note(o, extent, &count, &bad);
        } while (FindNextFileA(h, &fd));
        FindClose(h);
    }
#else
    DIR *d = opendir(dir);
    if (d) {
        const struct dirent *e = NULL;
        while ((e = readdir(d)) != NULL)
            if (cube_name(e->d_name, o)) cube_note(o, extent, &count, &bad);
        closedir(d);
    }
#endif
    if (count == 0) {
        snprintf(why, n, "%s: neither a zarr array (.zarray) nor a folder of z#####_y#####_x#####.tif cubes", dir);
        return -1;
    }
    if (bad) {
        snprintf(why, n, "%s: %ld cube id(s) off the 128-voxel grid", dir, bad);
        return -1;
    }
    for (int k = 0; k < 3; k++) shape[k] = extent[k];
    snprintf(cf->dir, sizeof cf->dir, "%s", dir);
#ifdef _OPENMP
    cf->nthreads = omp_get_max_threads();
#else
    cf->nthreads = 1;
#endif
    if (cf->nthreads < 1) cf->nthreads = 1;
    cf->scratch = (Arena_T *)calloc((size_t)cf->nthreads, sizeof(Arena_T));
    if (!cf->scratch) {
        snprintf(why, n, "out of host memory (cube scratch)");
        return -1;
    }
    printf("[vcr] %s: %ld TIFF cubes\n", dir, count);
    return 0;
}

static void close_cubes(CubeFiles *cf)
{
    if (cf->scratch)
        for (int t = 0; t < cf->nthreads; t++)
            if (cf->scratch[t]) Arena_dispose(&cf->scratch[t]);
    free(cf->scratch);
    cf->scratch = NULL;
}

typedef struct {
    TiffFloatWindow *plane[3];
} FloatFiles;

static void close_float_files(FloatFiles *f)
{
    for (int k = 0; k < 3; k++) { TiffIO_float_window_close(f->plane[k]); f->plane[k] = NULL; }
}

static int float_file_nodes(void *ctx, int row, int col, int height, int width, float *xyz, size_t plane_stride)
{
    FloatFiles *f = (FloatFiles *)ctx;
    for (int k = 0; k < 3; k++)
        if (TiffIO_float_window_read(f->plane[k], row, col, height, width, xyz + (size_t)k * plane_stride)) return -1;
    return 0;
}

static int load_tifxyz(Arena_T a, const char *dir, Scene *s, FloatFiles *files,
                        char *uuid, size_t uuid_len, char *why, size_t n)
{
    char p[1100];
    const char *err = NULL;
    snprintf(p, sizeof p, "%s/meta.json", dir);
    const JsonValue *m = Json_parse_file(a, p, &err);
    if (!m) { snprintf(why, n, "%s: %s", p, err ? err : "unreadable"); return -1; }
    const JsonValue *sc = Json_object_get(m, "scale");
    if (Json_array_len(sc) != 2 || Json_as_double(Json_array_get(sc, 0), 0) != 1.0
        || Json_as_double(Json_array_get(sc, 1), 0) != 1.0 || Json_object_get(m, "components")) {
        snprintf(why, n, "%s: only scale [1, 1] single-component tifxyz surfaces are rendered here", p);
        return -1;
    }
    const char *u = Json_as_string(Json_object_get(m, "uuid"));
    snprintf(uuid, uuid_len, "%s", u ? u : "render");
    const char *names[3] = {"x.tif", "y.tif", "z.tif"};
    int w[3] = {0, 0, 0}, h[3] = {0, 0, 0};
    for (int k = 0; k < 3; k++) {
        snprintf(p, sizeof p, "%s/%s", dir, names[k]);
        files->plane[k] = TiffIO_float_window_open(p, (size_t)64 << 20, &w[k], &h[k]);
        if (!files->plane[k]) { snprintf(why, n, "%s: cannot open bounded float TIFF reader", p); return -1; }
    }
    if (w[1] != w[0] || w[2] != w[0] || h[1] != h[0] || h[2] != h[0] || w[0] < 4 || h[0] < 4) {
        snprintf(why, n, "%s: x/y/z grids differ or are smaller than 4 x 4", dir);
        return -1;
    }
    s->W = w[0];
    s->H = h[0];
    s->nodes = float_file_nodes;
    s->node_ctx = files;
    return 0;
}

static int mkdirs(const char *path)
{
    char buf[1200];
    snprintf(buf, sizeof buf, "%s", path);
    for (char *q = buf + 1; *q; q++) {
        if (*q == '/' || *q == '\\') {
            const char keep = *q;
            if (q[-1] == ':') continue;
            *q = '\0';
            ves_mkdir(buf);
            *q = keep;
        }
    }
    ves_mkdir(buf);
    return 0;
}

static int write_text(const char *path, const char *text)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    size_t n = strlen(text);
    int ok = fwrite(text, 1, n, f) == n;
    return (fclose(f) == 0 && ok) ? 0 : -1;
}

typedef struct {
    char   arr[1100];      /* the array dir: <zarr_output>/0 */
    size_t written, skipped;
    int    failed;
    int    shard_cols;
    size_t files;
} ZarrOut;

static int write_band(void *ctx, int ty, int tx0, const uint8_t *band, int ntx, int N)
{
    ZarrOut *z = (ZarrOut *)ctx;
    char dir[1200];
    snprintf(dir, sizeof dir, z->shard_cols ? "%s/c/0/%d" : "%s/0/%d", z->arr, ty);
    mkdirs(dir);
    const size_t chunk = (size_t)N * VCR_TILE;
    size_t written = 0, skipped = 0;
    int failed = 0, tx = 0;
    if (z->shard_cols) {
        const int count = z->shard_cols, nshards = (ntx + count - 1) / count;
        if (tx0 % count) return -1;
        size_t files = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+ : written, skipped, failed, files)
        for (tx = 0; tx < nshards; tx++) {
            char p[1300];
            const int first = tx * count, supplied = ntx - first < count ? ntx - first : count;
            size_t stored = 0;
            snprintf(p, sizeof p, "%s/%d", dir, tx0 / count + tx);
            if (ZarrShard_write(p, band + (size_t)first * chunk, chunk, (size_t)supplied,
                                (size_t)count, &stored)) failed++;
            written += stored;
            skipped += (size_t)supplied - stored;
            files += stored != 0;
        }
        z->written += written;
        z->skipped += skipped;
        z->files += files;
        z->failed += failed;
        return failed ? -1 : 0;
    }
    #pragma omp parallel for schedule(dynamic, 4) reduction(+ : written, skipped, failed)
    for (tx = 0; tx < ntx; tx++) {
        const uint8_t *c = band + (size_t)tx * chunk;
        size_t i = 0;
        while (i < chunk && c[i] == 0) i++;
        if (i == chunk) { skipped++; continue; }          /* all fill: no chunk file (fill_value 0) */
        char p[1300];
        snprintf(p, sizeof p, "%s/%d", dir, tx0 + tx);
        FILE *f = fopen(p, "wb");
        if (!f || fwrite(c, 1, chunk, f) != chunk) failed++;
        if (f && fclose(f) != 0) failed++;
        written++;
    }
    z->written += written;
    z->skipped += skipped;
    z->files += written;
    z->failed += failed;
    return failed ? -1 : 0;
}

static int write_zarr_meta(const VcrRender_job *job, const Scene *s, const char *uuid)
{
    char p[1200], text[4096];
    const char *out = job->zarr_output;
    const double vs = job->voxel_size > 0 ? job->voxel_size : 1.0;
    const char *unit = job->voxel_unit ? job->voxel_unit : "micrometer";
    if (!job->zarr_shard_cols) {
        snprintf(p, sizeof p, "%s/.zgroup", out);
        if (write_text(p, "{\"zarr_format\": 2}")) return -1;
    }
    snprintf(text, sizeof text,
             "{\n \"canvas_size\": [%d, %d],\n \"chunk_size\": [%d, %d, %d],\n"
             " \"multiscales\": [{\"axes\": [{\"name\": \"z\", \"type\": \"space\", \"unit\": \"%s\"}, "
             "{\"name\": \"y\", \"type\": \"space\", \"unit\": \"%s\"}, {\"name\": \"x\", \"type\": \"space\", "
             "\"unit\": \"%s\"}],\n  \"datasets\": [{\"coordinateTransformations\": [{\"scale\": [%.6g, %.6g, %.6g], "
             "\"type\": \"scale\"}, {\"translation\": [0.0, 0.0, 0.0], \"type\": \"translation\"}], \"path\": \"0\"}],\n"
             "  \"name\": \"render\", \"version\": \"0.4\"}],\n"
             " \"note_axes_order\": \"ZYX (slice, row, col)\",\n \"num_slices\": %d,\n \"slice_step\": %.6g,\n"
             " \"source_group\": 0,\n \"surface\": \"%s\",\n"
             " \"renderer\": \"vcr_gpu: villa vc_render_tifxyz 6e53201a's surface-volume render transcribed to CUDA "
             "(src/flatten/vcr_render.c)\"\n}\n",
             s->W, s->H, job->num_slices, VCR_CHUNK, VCR_CHUNK, unit, unit, unit, vs, vs, vs, job->num_slices,
             job->slice_step, uuid);
    if (job->zarr_shard_cols) {
        char group[4600], array[1200];
        snprintf(group, sizeof group, "{\"zarr_format\": 3, \"node_type\": \"group\", \"attributes\": %s}\n", text);
        snprintf(p, sizeof p, "%s/zarr.json", out);
        snprintf(array, sizeof array, "%s/0", out);
        if (write_text(p, group)) return -1;
        return ZarrShard_array_meta(array, (size_t)job->num_slices, (size_t)s->H, (size_t)s->W,
                                     VCR_CHUNK, (size_t)job->zarr_shard_cols);
    }
    snprintf(p, sizeof p, "%s/.zattrs", out);
    if (write_text(p, text)) return -1;
    snprintf(text, sizeof text,
             "{\n \"zarr_format\": 2,\n \"shape\": [%d, %d, %d],\n \"chunks\": [%d, %d, %d],\n \"dtype\": \"|u1\",\n"
             " \"compressor\": null,\n \"fill_value\": 0,\n \"order\": \"C\",\n \"filters\": null,\n"
             " \"dimension_separator\": \"/\"\n}\n",
             job->num_slices, s->H, s->W, job->num_slices, VCR_CHUNK, VCR_CHUNK);
    snprintf(p, sizeof p, "%s/0/.zarray", out);
    return write_text(p, text);
}

int VcrRender_run(const VcrRender_job *job)
{
    Arena_T a = Arena_new();
    Scene s;
    VolFiles vf;
    CubeFiles cf;
    FloatFiles ff = {{NULL, NULL, NULL}};
    char why[2048] = "", uuid[256], zp[1100];
    uint8_t *comp = NULL;
    FILE *zf = NULL;
    int rc = 1, zarr = 0;
    double t0 = wall_s();
    memset(&s, 0, sizeof s);
    memset(&vf, 0, sizeof vf);
    memset(&cf, 0, sizeof cf);
    if (load_tifxyz(a, job->segmentation, &s, &ff, uuid, sizeof uuid, why, sizeof why) != 0)
        goto out;
    snprintf(zp, sizeof zp, "%s/.zarray", job->volume);   /* a zarr array, else the grid's TIFF cubes */
    zf = fopen(zp, "rb");
    if (zf) {
        zarr = 1;
        fclose(zf);
    }
    if (zarr ? open_volume(a, job->volume, &vf, s.shape, why, sizeof why) != 0
             : open_cubes(job->volume, &cf, s.shape, why, sizeof why) != 0)
        goto out;
    for (int k = 0; k < 3; k++) s.cg[k] = (s.shape[k] + VCR_CHUNK - 1) / VCR_CHUNK;
    s.load = zarr ? vol_chunk : cube_chunk;
    s.ctx = zarr ? (void *)&vf : (void *)&cf;
    printf("[vcr] %s: %d x %d canvas, volume %ld x %ld x %ld (%s); %s (%.1fs)\n", uuid, s.H, s.W, s.shape[0],
           s.shape[1], s.shape[2], zarr ? "zarr" : "TIFF cubes", job->composite ? "max composite" : "layer stack",
           wall_s() - t0);
    CudaDrv cu = CudaDrv_open(KERNELS, KERNEL_NAMES, 7, why, sizeof why);
    if (!cu) goto out;
    printf("[vcr] device %s\n", CudaDrv_device(cu));
    fflush(stdout);
    if (job->composite) {
        comp = (uint8_t *)malloc((size_t)s.W * (size_t)s.H);
        if (!comp) { snprintf(why, sizeof why, "out of host memory (composite)"); CudaDrv_close(cu); goto out; }
        if (gpu_render(cu, &s, job, NULL, NULL, comp, why, sizeof why) == 0) {
            char p[1400];
            mkdirs(job->tif_output);
            snprintf(p, sizeof p, "%s/%s.tif", job->tif_output, uuid);
            if (TiffIO_save(p, comp, 1, s.H, s.W) == 0) {
                printf("[vcr] wrote %s\n", p);
                rc = 0;
            } else {
                snprintf(why, sizeof why, "%s: write failed", p);
            }
        }
    } else {
        ZarrOut zo;
        memset(&zo, 0, sizeof zo);
        zo.shard_cols = job->zarr_shard_cols;
        snprintf(zo.arr, sizeof zo.arr, "%s/0", job->zarr_output);
        mkdirs(zo.arr);
        if (gpu_render(cu, &s, job, write_band, &zo, NULL, why, sizeof why) == 0 && zo.failed == 0
            && write_zarr_meta(job, &s, uuid) == 0) {
            printf("[vcr] wrote %s: [%d, %d, %d], %zu chunks (%zu all-fill chunks not stored)\n", job->zarr_output,
                   job->num_slices, s.H, s.W, zo.written, zo.skipped);
            if (zo.shard_cols) printf("[vcr] storage: %zu shard files, up to %d chunks each\n", zo.files, zo.shard_cols);
            rc = 0;
        } else if (!why[0]) {
            snprintf(why, sizeof why, "writing %s failed", job->zarr_output);
        }
    }
    CudaDrv_close(cu);
out:
    if (rc) fprintf(stderr, "vcr_gpu: %s\n", why);
    else printf("[vcr] done in %.1fs\n", wall_s() - t0);
    free(comp);
    close_cubes(&cf);
    close_float_files(&ff);
    Arena_dispose(&a);
    return rc;
}

int VcrRender_storage_probe(const char *input_array, const char *output_group, int shard_cols)
{
    char path[1600];
    Arena_T arena = Arena_new();
    uint8_t *band = NULL;
    int result = 1;
    if (input_array == NULL || output_group == NULL || shard_cols < 0 || shard_cols > ZARR_SHARD_MAX_CHUNKS) goto done;
    snprintf(path, sizeof path, "%s/.zarray", input_array);
    const JsonValue *meta = Json_parse_file(arena, path, NULL);
    const JsonValue *shape = Json_object_get(meta, "shape"), *chunks = Json_object_get(meta, "chunks");
    const char *dtype = Json_as_string(Json_object_get(meta, "dtype"));
    const char *order = Json_as_string(Json_object_get(meta, "order"));
    const char *separator = Json_as_string(Json_object_get(meta, "dimension_separator"));
    if (meta == NULL || Json_array_len(shape) != 3 || Json_array_len(chunks) != 3 ||
        dtype == NULL || strcmp(dtype, "|u1") || (order && strcmp(order, "C")) ||
        Json_member_long(meta, "zarr_format", 0) != 2 || Json_member_double(meta, "fill_value", 0) != 0 ||
        Json_type(Json_object_get(meta, "compressor")) != JSON_NULL ||
        Json_type(Json_object_get(meta, "filters")) != JSON_NULL ||
        (separator && strcmp(separator, ".") && strcmp(separator, "/"))) goto done;
    int dim[3];
    for (size_t d = 0; d < 3; d++) {
        double n = Json_as_double(Json_array_get(shape, d), 0);
        if (n < 1 || n > INT_MAX - VCR_CHUNK || n != floor(n)) goto done;
        dim[d] = (int)n;
    }
    if (dim[0] > 255 || Json_as_long(Json_array_get(chunks, 0), 0) != dim[0] ||
        Json_as_long(Json_array_get(chunks, 1), 0) != VCR_CHUNK ||
        Json_as_long(Json_array_get(chunks, 2), 0) != VCR_CHUNK) goto done;
    const int ntx = (dim[2] + VCR_CHUNK - 1) / VCR_CHUNK, nty = (dim[1] + VCR_CHUNK - 1) / VCR_CHUNK;
    const size_t chunk_bytes = (size_t)dim[0] * VCR_TILE;
    if ((size_t)ntx > SIZE_MAX / chunk_bytes) goto done;
    const size_t band_bytes = (size_t)ntx * chunk_bytes;
    band = (uint8_t *)malloc(band_bytes);
    if (band == NULL || ves_ensure_parent_dir(output_group) || ves_mkdir(output_group)) goto done;
    ZarrOut out;
    memset(&out, 0, sizeof out);
    out.shard_cols = shard_cols;
    snprintf(out.arr, sizeof out.arr, "%s/0", output_group);
    if (ves_mkdir(out.arr)) goto done;
    double read_seconds = 0, write_seconds = 0;
    size_t read_files = 0;
    for (int ty = 0; ty < nty; ty++) {
        double start = wall_s();
        memset(band, 0, band_bytes);
        int failed = 0, tx = 0;
        #pragma omp parallel for schedule(dynamic, 4) reduction(+ : failed, read_files)
        for (tx = 0; tx < ntx; tx++) {
            char source[1600];
            snprintf(source, sizeof source, separator && separator[0] == '/' ? "%s/0/%d/%d" : "%s/0.%d.%d",
                     input_array, ty, tx);
            FILE *file = fopen(source, "rb");
            if (file == NULL) { if (errno != ENOENT) failed++; continue; }
            read_files++;
            if (fread(band + (size_t)tx * chunk_bytes, 1, chunk_bytes, file) != chunk_bytes || fgetc(file) != EOF)
                failed++;
            if (fclose(file)) failed++;
        }
        read_seconds += wall_s() - start;
        if (failed) goto done;
        start = wall_s();
        if (write_band(&out, ty, 0, band, ntx, dim[0])) goto done;
        write_seconds += wall_s() - start;
        if ((ty + 1) % 16 == 0 || ty + 1 == nty) {
            printf("[vcr storage] band %d/%d read %.3fs write %.3fs\n", ty + 1, nty, read_seconds, write_seconds);
            fflush(stdout);
        }
    }
    Scene scene;
    VcrRender_job job;
    memset(&scene, 0, sizeof scene);
    memset(&job, 0, sizeof job);
    scene.H = dim[1]; scene.W = dim[2];
    job.num_slices = dim[0]; job.slice_step = 1; job.zarr_output = output_group; job.zarr_shard_cols = shard_cols;
    if (write_zarr_meta(&job, &scene, "storage-probe")) goto done;
    printf("[vcr storage] complete: %zu read files, %zu output files, %zu chunks, %zu skipped, "
           "%zu scratch bytes; read %.6f write %.6f seconds\n",
           read_files, out.files, out.written, out.skipped, band_bytes, read_seconds, write_seconds);
    result = 0;
done:
    if (result) fprintf(stderr, "vcr storage: invalid input, nonfresh output or I/O failure\n");
    free(band);
    Arena_dispose(&arena);
    return result;
}

/* ---------------------------------------------------------------------------
 * Selftest: the CPU transcription of villa's arithmetic against the kernels.
 * ------------------------------------------------------------------------- */
typedef struct {
    long     shape[3];
    long     cg[3];
    uint8_t *chunks;      /* cg0*cg1*cg2 chunks of 128^3 */
    uint8_t *stored;      /* per chunk: 1 stored */
    int      fail_read;   /* selftest I/O fault injection */
} MemVol;

static int mem_chunk(void *ctx, long cz, long cy, long cx, uint8_t *dst)
{
    const MemVol *v = (const MemVol *)ctx;
    const size_t k = (size_t)((cz * v->cg[1] + cy) * v->cg[2] + cx);
    if (!v->stored[k]) return 0;
    if (dst && v->fail_read) return -1;
    if (dst) memcpy(dst, v->chunks + k * VCR_CHUNK3, VCR_CHUNK3);
    return 1;
}

static float ref_vox(const MemVol *v, int iz, int iy, int ix)
{
    if ((unsigned)iz >= (unsigned)v->shape[0] || (unsigned)iy >= (unsigned)v->shape[1]
        || (unsigned)ix >= (unsigned)v->shape[2])
        return 0.0f;
    const size_t k = (size_t)(((long)(iz >> 7) * v->cg[1] + (iy >> 7)) * v->cg[2] + (ix >> 7));
    if (!v->stored[k]) return 0.0f;
    return (float)v->chunks[k * VCR_CHUNK3 + ((size_t)(iz & 127) << 14) + ((size_t)(iy & 127) << 7) + (size_t)(ix & 127)];
}

/* ChunkSampler::sampleTrilinear */
static float ref_trilinear(const MemVol *v, float vz, float vy, float vx)
{
    const int iz = (int)vz, iy = (int)vy, ix = (int)vx;
    const float c000 = ref_vox(v, iz, iy, ix), c100 = ref_vox(v, iz + 1, iy, ix);
    const float c010 = ref_vox(v, iz, iy + 1, ix), c110 = ref_vox(v, iz + 1, iy + 1, ix);
    const float c001 = ref_vox(v, iz, iy, ix + 1), c101 = ref_vox(v, iz + 1, iy, ix + 1);
    const float c011 = ref_vox(v, iz, iy + 1, ix + 1), c111 = ref_vox(v, iz + 1, iy + 1, ix + 1);
    const float fz = vz - (float)iz, fy = vy - (float)iy, fx = vx - (float)ix;
    const float c00 = fmaf(fx, c001 - c000, c000), c01 = fmaf(fx, c011 - c010, c010);
    const float c10 = fmaf(fx, c101 - c100, c100), c11 = fmaf(fx, c111 - c110, c110);
    const float c0 = fmaf(fy, c01 - c00, c00), c1 = fmaf(fy, c11 - c10, c10);
    return fmaf(fz, c1 - c0, c0);
}

static void ref_node(const Scene *s, int r, int c, float *p)
{
    const size_t i = (size_t)r * (size_t)s->W + (size_t)c;
    p[0] = s->X[i];
    p[1] = s->Y[i];
    p[2] = s->Z[i];
    if (p[0] == -1.0f) p[0] = p[1] = p[2] = NAN;
}

/* QuadSurface::gen + grid_normal_int + prepareBaseAndDirs for canvas pixel (r, c) */
static int ref_geom(const Scene *s, int r, int c, int flip, float *b, float *d)
{
    const int W = s->W, H = s->H, rn = r + 1, cn = c + 1;
    float q[3], l[3], rt[3], u[3], w[3], a[3], e[3], f[3], g[3];
    if (rn >= H || cn >= W) return 0;
    ref_node(s, rn, cn, q);
    if (isnan(q[0]) || isnan(q[1]) || isnan(q[2]) || (q[0] == -1.f && q[1] == -1.f && q[2] == -1.f)) return 0;
    if (rn > H - 2 || cn > W - 2) return 0;
    ref_node(s, rn, cn - 1, l);
    ref_node(s, rn, cn + 1, rt);
    ref_node(s, rn - 1, cn, u);
    ref_node(s, rn + 1, cn, w);
    const float xv0 = rt[0] - l[0], xv1 = rt[1] - l[1], xv2 = rt[2] - l[2];
    const float yv0 = w[0] - u[0], yv1 = w[1] - u[1], yv2 = w[2] - u[2];
    float n0 = xv1 * yv2 - xv2 * yv1, n1 = xv2 * yv0 - xv0 * yv2, n2 = xv0 * yv1 - xv1 * yv0;
    const float len2 = n0 * n0 + n1 * n1 + n2 * n2;
    if (len2 == 0.0f || isnan(len2)) return 0;
    const float inv = 1.0f / sqrtf(len2);
    n0 *= inv;
    n1 *= inv;
    n2 *= inv;
    if (flip) { n0 = -n0; n1 = -n1; n2 = -n2; }
    const float L2 = n0 * n0 + n1 * n1 + n2 * n2;
    if (L2 > 0.0f) {
        const float ia = 1.0f / sqrtf(L2);
        n0 *= ia;
        n1 *= ia;
        n2 *= ia;
    }
    d[0] = n0;
    d[1] = n1;
    d[2] = n2;
    ref_node(s, r, c, a);
    ref_node(s, r, c + 1, e);
    ref_node(s, r + 1, c, f);
    ref_node(s, r + 1, c + 1, g);
    for (int k = 0; k < 3; k++) {
        const float top = a[k] * 0.5f + e[k] * 0.5f, bot = f[k] * 0.5f + g[k] * 0.5f;
        b[k] = top * 0.5f + bot * 0.5f;
    }
    return 1;
}

static uint8_t ref_layer(const Scene *s, const MemVol *v, int r, int c, int flip, float off)
{
    float b[3], d[3];
    if (!ref_geom(s, r, c, flip, b, d)) return 0;
    const float px = b[0] + d[0] * off, py = b[1] + d[1] * off, pz = b[2] + d[2] * off;
    if (!(pz >= 0.f && py >= 0.f && px >= 0.f && pz < (float)v->shape[0] && py < (float)v->shape[1]
          && px < (float)v->shape[2]))
        return 0;
    float val = ref_trilinear(v, pz, py, px);
    if (val < 0.f) val = 0.f;
    if (val > 255.f) val = 255.f;
    return (uint8_t)val;
}

static uint8_t ref_ctmax(const Scene *s, const MemVol *v, int r, int c, int flip, int z0, int z1, float step)
{
    float b[3], d[3];
    if (!ref_geom(s, r, c, flip, b, d) || !isfinite(b[0])) return 0;
    float mx = 0.0f;
    for (int k = z0; k <= z1; k++) {
        const float z = (float)k * step;
        const float px = b[0] + d[0] * z, py = b[1] + d[1] * z, pz = b[2] + d[2] * z;
        float val = 0.0f;
        if (pz >= 0.f && py >= 0.f && px >= 0.f && pz < (float)v->shape[0] && py < (float)v->shape[1]
            && px < (float)v->shape[2]) {
            int iz = (int)(pz + 0.5f), iy = (int)(py + 0.5f), ix = (int)(px + 0.5f);
            if (iz >= v->shape[0]) iz = (int)v->shape[0] - 1;
            if (iy >= v->shape[1]) iy = (int)v->shape[1] - 1;
            if (ix >= v->shape[2]) ix = (int)v->shape[2] - 1;
            val = ref_vox(v, iz, iy, ix);
        }
        mx = val > mx ? val : mx;
    }
    if (mx > 255.f) mx = 255.f;
    return (uint8_t)mx;
}

typedef struct {
    uint8_t *stack;    /* [N][H][W] */
    int      W, H;
    size_t   bad_padding;
} MemOut;

static int mem_band(void *ctx, int ty, int tx0, const uint8_t *band, int ntx, int N)
{
    MemOut *m = (MemOut *)ctx;
    for (int tx = 0; tx < ntx; tx++)
        for (int i = 0; i < N; i++)
            for (int ly = 0; ly < VCR_CHUNK; ly++)
                for (int lx = 0; lx < VCR_CHUNK; lx++) {
                    const int r = ty * VCR_CHUNK + ly, c = (tx0 + tx) * VCR_CHUNK + lx;
                    if (r < m->H && c < m->W)
                        m->stack[((size_t)i * m->H + r) * m->W + c] =
                            band[(((size_t)tx * N + i) * VCR_CHUNK + ly) * VCR_CHUNK + lx];
                    else if (band[(((size_t)tx * N + i) * VCR_CHUNK + ly) * VCR_CHUNK + lx])
                        m->bad_padding++;
                }
    return 0;
}

/* The TIFF cube folder against the zarr array: the same 2 x 2 x 3 chunk volume written both ways (one chunk absent)
 * reads the same through both loaders, every chunk (staged on OpenMP threads, as the render stages them); the cube
 * folder's shape is its far corner; a folder with a cube off the 128 grid, or with no cubes, is refused.  No GPU. */
static int cube_selftest(void)
{
    char root[512], zdir[700], cdir[700], bad[700], p[1300], why[1024];
    const long cg[3] = {2, 2, 3};
    const long nch = cg[0] * cg[1] * cg[2];
    int fails = 0;
    if (ves_selftest_dir(root, sizeof root, "vcr_cubes") != 0) {
        fprintf(stderr, "vcr_render cube selftest: no scratch folder\n");
        return 1;
    }
    snprintf(zdir, sizeof zdir, "%s/raw.zarr", root);
    snprintf(cdir, sizeof cdir, "%s/cubes_RAW", root);
    snprintf(bad, sizeof bad, "%s/cubes_bad", root);
    mkdirs(zdir);
    mkdirs(cdir);
    mkdirs(bad);
    uint8_t *chunk = (uint8_t *)malloc(VCR_CHUNK3), *ref = (uint8_t *)malloc(VCR_CHUNK3);
    uint8_t *stage = (uint8_t *)malloc(VCR_CHUNK3 * (size_t)nch);
    int *st = (int *)calloc((size_t)nch, sizeof(int));
    Arena_T a = Arena_new();
    VolFiles vf;
    CubeFiles cf, cfb;
    long sz[3] = {0, 0, 0}, sc[3] = {0, 0, 0};
    memset(&vf, 0, sizeof vf);
    memset(&cf, 0, sizeof cf);
    memset(&cfb, 0, sizeof cfb);
    if (!chunk || !ref || !stage || !st) {
        fprintf(stderr, "vcr_render cube selftest: out of memory\n");
        fails++;
        goto end;
    }
    snprintf(p, sizeof p, "%s/.zarray", zdir);
    write_text(p, "{\n \"zarr_format\": 2,\n \"shape\": [256, 256, 384],\n \"chunks\": [128, 128, 128],\n"
                  " \"dtype\": \"|u1\",\n \"compressor\": null,\n \"fill_value\": 0,\n \"order\": \"C\",\n"
                  " \"filters\": null,\n \"dimension_separator\": \"/\"\n}\n");
    for (long cz = 0; cz < cg[0]; cz++)
        for (long cy = 0; cy < cg[1]; cy++)
            for (long cx = 0; cx < cg[2]; cx++) {
                if (cz == 1 && cy == 0 && cx == 2) continue;                 /* absent both ways */
                for (size_t i = 0; i < VCR_CHUNK3; i++)
                    chunk[i] = (uint8_t)(((size_t)(cz * 31 + cy * 17 + cx * 7) + i * 13 + (i >> 9)) & 0xffu);
                snprintf(p, sizeof p, "%s/%ld/%ld", zdir, cz, cy);
                mkdirs(p);
                snprintf(p, sizeof p, "%s/%ld/%ld/%ld", zdir, cz, cy, cx);
                FILE *f = fopen(p, "wb");
                if (!f || fwrite(chunk, 1, VCR_CHUNK3, f) != VCR_CHUNK3) fails++;
                if (f) fclose(f);
                snprintf(p, sizeof p, "%s/z%05ld_y%05ld_x%05ld.tif", cdir, cz * VCR_CHUNK, cy * VCR_CHUNK,
                         cx * VCR_CHUNK);
                if (TiffIO_save(p, chunk, VCR_CHUNK, VCR_CHUNK, VCR_CHUNK) != 0) fails++;
            }
    if (fails) {
        fprintf(stderr, "vcr_render cube selftest: could not write the volume\n");
        goto end;
    }
    if (open_volume(a, zdir, &vf, sz, why, sizeof why) != 0 || open_cubes(cdir, &cf, sc, why, sizeof why) != 0) {
        fprintf(stderr, "vcr_render cube selftest: %s\n", why);
        fails++;
        goto end;
    }
    if (sc[0] != 256 || sc[1] != 256 || sc[2] != 384) {
        fprintf(stderr, "vcr_render cube selftest: cube shape %ld x %ld x %ld (want 256 x 256 x 384)\n", sc[0], sc[1],
                sc[2]);
        fails++;
    }
    long li;
    #pragma omp parallel for schedule(dynamic, 1)
    for (li = 0; li < nch; li++)
        st[li] = cube_chunk(&cf, li / (cg[1] * cg[2]), (li / cg[2]) % cg[1], li % cg[2],
                            stage + (size_t)li * VCR_CHUNK3);
    for (long i = 0; i < nch; i++) {
        const long cz = i / (cg[1] * cg[2]), cy = (i / cg[2]) % cg[1], cx = i % cg[2];
        const int zs = vol_chunk(&vf, cz, cy, cx, ref);
        const int ze = vol_chunk(&vf, cz, cy, cx, NULL), ce = cube_chunk(&cf, cz, cy, cx, NULL);
        if (zs != st[i] || ze != ce || (zs == 1 && memcmp(ref, stage + (size_t)i * VCR_CHUNK3, VCR_CHUNK3) != 0)) {
            fprintf(stderr, "vcr_render cube selftest: chunk %ld %ld %ld: zarr %d/%d, cubes %d/%d%s\n", cz, cy, cx, zs,
                    ze, st[i], ce, zs == 1 && st[i] == 1 ? ", bytes differ" : "");
            fails++;
        }
    }
    if (st[1 * cg[1] * cg[2] + 0 * cg[2] + 2] != 0) {
        fprintf(stderr, "vcr_render cube selftest: the absent cube read as present\n");
        fails++;
    }
    snprintf(p, sizeof p, "%s/z00100_y00000_x00000.tif", bad);           /* an origin off the 128 grid */
    if (TiffIO_save(p, chunk, VCR_CHUNK, VCR_CHUNK, VCR_CHUNK) != 0
        || open_cubes(bad, &cfb, sz, why, sizeof why) == 0 || open_cubes(root, &cfb, sz, why, sizeof why) == 0) {
        fprintf(stderr, "vcr_render cube selftest: an off-grid or empty cube folder was accepted\n");
        fails++;
    }
end:
    close_cubes(&cf);
    close_cubes(&cfb);
    Arena_dispose(&a);
    free(chunk);
    free(ref);
    free(stage);
    free(st);
    printf("vcr_render cube selftest: %s\n", fails ? "FAIL" : "PASS");
    if (!fails) ves_selftest_clean(root);
    return fails ? 1 : 0;
}

int VcrRender_selftest(void)
{
    char why[2048] = "", scratch[1200];
    const int cube_fails = cube_selftest();
    if (ves_selftest_dir(scratch, sizeof scratch, "vcr_tiff_windows")) return 1;
    const int tiff_fails = TiffIO_selftest(scratch);
    CudaDrv cu = CudaDrv_open(KERNELS, KERNEL_NAMES, 7, why, sizeof why);
    if (!cu) {
        fprintf(stderr, "vcr_render selftest: SKIP (%s)\n", why);
        if (!tiff_fails) ves_selftest_clean(scratch);
        return cube_fails || tiff_fails;
    }
    /* a curved sheet over three chunk planes, with holes, crossing the volume's far faces */
    const int W = 1173, H = 289;
    float *X = (float *)malloc(sizeof(float) * (size_t)W * H), *Y = (float *)malloc(sizeof(float) * (size_t)W * H);
    float *Z = (float *)malloc(sizeof(float) * (size_t)W * H);
    MemVol v;
    memset(&v, 0, sizeof v);
    v.shape[0] = 256; v.shape[1] = 340; v.shape[2] = 1160;
    for (int k = 0; k < 3; k++) v.cg[k] = (v.shape[k] + VCR_CHUNK - 1) / VCR_CHUNK;
    const size_t ncg = (size_t)(v.cg[0] * v.cg[1] * v.cg[2]);
    v.chunks = (uint8_t *)calloc(ncg, VCR_CHUNK3);
    v.stored = (uint8_t *)calloc(ncg, 1);
    int fails = 0;
    if (!X || !Y || !Z || !v.chunks || !v.stored) { fprintf(stderr, "vcr_render selftest: out of memory\n"); fails++; goto end; }
    for (int r = 0; r < H; r++)
        for (int c = 0; c < W; c++) {
            const size_t i = (size_t)r * W + c;
            X[i] = 3.0f + 1.01f * (float)c + 2.5f * sinf(0.21f * (float)r);
            Y[i] = 40.0f + 0.93f * (float)r + 0.13f * (float)c;
            Z[i] = 118.0f + 7.0f * sinf(0.11f * (float)c) + 0.37f * (float)r;
            if ((r >= 20 && r < 24 && c >= 30 && c < 35) || (r == 50 && c == 60) || (c == 99 && r % 7 == 0))
                X[i] = Y[i] = Z[i] = -1.0f;
        }
    for (size_t k = 0; k < ncg; k++) {
        v.stored[k] = (uint8_t)(k != 4);        /* chunk (1, 0, 0), crossed by the sheet, absent: reads 0 */
        uint32_t h = 2166136261u ^ (uint32_t)k;
        for (size_t j = 0; j < VCR_CHUNK3; j++) {
            h ^= (uint32_t)j;
            h *= 16777619u;
            v.chunks[k * VCR_CHUNK3 + j] = (uint8_t)(h >> 13);
        }
    }
    Scene s;
    memset(&s, 0, sizeof s);
    s.W = W; s.H = H; s.X = X; s.Y = Y; s.Z = Z;
    for (int k = 0; k < 3; k++) { s.shape[k] = v.shape[k]; s.cg[k] = v.cg[k]; }
    s.load = mem_chunk;
    s.ctx = &v;
    const int cases[3][3] = {{21, 1, 0}, {24, 1, 0}, {21, 0, 0}};    /* slices, flip, - */
    for (int t = 0; t < 3; t++) {
        VcrRender_job job;
        memset(&job, 0, sizeof job);
        job.num_slices = cases[t][0];
        job.slice_step = 1.0;
        job.flip_normals = cases[t][1];
        MemOut mo = {(uint8_t *)calloc((size_t)job.num_slices * W * H, 1), W, H};
        VcrCacheStats stats;
        float offs[64];
        layer_offsets(job.num_slices, job.slice_step, offs);
        if (!mo.stack || gpu_render_cached(cu, &s, &job, mem_band, &mo, NULL, why, sizeof why, 2, 8, &stats) != 0) {
            fprintf(stderr, "vcr_render selftest: layers case %d failed: %s\n", t, why);
            fails++;
            free(mo.stack);
            continue;
        }
        size_t diff = 0, nonzero = 0;
        for (int i = 0; i < job.num_slices; i++)
            for (int r = 0; r < H; r++)
                for (int c = 0; c < W; c++) {
                    const uint8_t want = ref_layer(&s, &v, r, c, job.flip_normals, offs[i]);
                    const uint8_t got = mo.stack[((size_t)i * H + r) * W + c];
                    diff += want != got;
                    nonzero += want != 0;
                }
        if (diff || nonzero < (size_t)W * H || mo.bad_padding || !stats.evicted || !stats.splits || !stats.hits) {
            fprintf(stderr, "vcr_render selftest: layers case %d: %zu of %zu samples differ (%zu nonzero)\n", t, diff,
                    (size_t)job.num_slices * W * H, nonzero);
            fprintf(stderr, "  cache evictions=%zu splits=%zu hits=%zu bad padding=%zu\n", stats.evicted,
                    stats.splits, stats.hits, mo.bad_padding);
            fails++;
        }
        if (t == 0) {
            MemOut resident = {(uint8_t *)calloc((size_t)job.num_slices * W * H, 1), W, H};
            if (!resident.stack || gpu_render_resident(cu, &s, &job, mem_band, &resident, NULL, why, sizeof why) ||
                memcmp(resident.stack, mo.stack, (size_t)job.num_slices * W * H) || resident.bad_padding) {
                fprintf(stderr, "vcr_render selftest: resident/tiled parity failed: %s\n", why);
                fails++;
            }
            free(resident.stack);
            FloatFiles files = {{NULL, NULL, NULL}};
            const float *planes[3] = {X, Y, Z};
            int file_bad = 0;
            for (int k = 0; k < 3; k++) {
                char path[1400];
                int width = 0, height = 0;
                snprintf(path, sizeof path, "%s/coordinate_%d.tif", scratch, k);
                if (k == 0 ? TiffIO_save_float2d(path, planes[k], W, H) :
                    k == 1 ? TiffIO_save_float2d_deflate(path, planes[k], W, H) :
                             TiffIO_save_float2d_predicted(path, planes[k], W, H)) file_bad++;
                files.plane[k] = TiffIO_float_window_open(path, (size_t)2 << 20, &width, &height);
                if (!files.plane[k] || width != W || height != H) file_bad++;
            }
            Scene streamed_scene = s;
            streamed_scene.nodes = float_file_nodes;
            streamed_scene.node_ctx = &files;
            streamed_scene.X = streamed_scene.Y = streamed_scene.Z = NULL;
            MemOut streamed = {(uint8_t *)calloc((size_t)job.num_slices * W * H, 1), W, H};
            if (file_bad || !streamed.stack ||
                gpu_render_cached(cu, &streamed_scene, &job, mem_band, &streamed, NULL, why, sizeof why, 2, 8, NULL) ||
                memcmp(streamed.stack, mo.stack, (size_t)job.num_slices * W * H) || streamed.bad_padding) {
                fprintf(stderr, "vcr_render selftest: TIFF-window/resident parity failed: %s\n", why);
                fails++;
            }
            free(streamed.stack);
            close_float_files(&files);
            v.fail_read = 1;
            if (!gpu_render_cached(cu, &s, &job, mem_band, &mo, NULL, why, sizeof why, 2, 8, NULL) ||
                !strstr(why, "volume chunk read failed")) {
                fprintf(stderr, "vcr_render selftest: RAW read error was not reported: %s\n", why);
                fails++;
            }
            v.fail_read = 0;
            if (!gpu_render_cached(cu, &s, &job, mem_band, &mo, NULL, why, sizeof why, 2, 1, NULL) ||
                !strstr(why, "one sample footprint")) {
                fprintf(stderr, "vcr_render selftest: insufficient cache was not reported: %s\n", why);
                fails++;
            }
        }
        free(mo.stack);
    }
    {
        VcrRender_job job;
        memset(&job, 0, sizeof job);
        job.composite = 1;
        job.slice_step = 1.0;
        job.flip_normals = 1;
        job.comp_start = -4;
        job.comp_end = 4;
        uint8_t *comp = (uint8_t *)calloc((size_t)W * H, 1);
        if (!comp || gpu_render_cached(cu, &s, &job, NULL, NULL, comp, why, sizeof why, 2, 8, NULL) != 0) {
            fprintf(stderr, "vcr_render selftest: composite failed: %s\n", why);
            fails++;
        } else {
            size_t diff = 0, nonzero = 0;
            for (int r = 0; r < H; r++)
                for (int c = 0; c < W; c++) {
                    const uint8_t want = ref_ctmax(&s, &v, r, c, 1, -4, 4, 1.0f);
                    diff += want != comp[(size_t)r * W + c];
                    nonzero += want != 0;
                }
            if (diff || nonzero < (size_t)W * H / 2) {
                fprintf(stderr, "vcr_render selftest: composite: %zu of %d pixels differ (%zu nonzero)\n", diff, W * H,
                        nonzero);
                fails++;
            }
        }
        free(comp);
    }
end:
    printf("vcr_render selftest: %s on %s\n", fails ? "FAIL" : "PASS", CudaDrv_device(cu));
    CudaDrv_close(cu);
    free(X);
    free(Y);
    free(Z);
    free(v.chunks);
    free(v.stored);
    if (!fails && !tiff_fails) ves_selftest_clean(scratch);
    return (fails || cube_fails || tiff_fails) ? 1 : 0;
}
