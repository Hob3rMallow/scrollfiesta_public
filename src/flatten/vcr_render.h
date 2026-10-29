#ifndef VCR_RENDER_INCLUDED
#define VCR_RENDER_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* ============================================================================
 * vcr_render.h -- villa vc_render_tifxyz's surface-volume render on the GPU.
 *
 * User 2026-10-07: "can we move block rendering to the GPU please? ... this
 * would mean porting vc_render_tifxyz to the GPU ... and trying to get
 * near-identical outputs."  This is that port, for the two calls
 * scripts/ink/region_ink_all.py makes per sheet block (--scale 1
 * --group-idx 0 of a scale-1 tifxyz, slice step 1, --flip-normals):
 *
 *   layers:    --num-slices N --zarr-output O  -> O/0, uint8 [N, H, W],
 *              chunks [N, 128, 128] (uncompressed: compressor null)
 *   composite: --composite-collapse --accum-type max --composite-start/end
 *              --tif-output D                  -> D/<uuid>.tif, uint8 [H, W]
 *
 * Every float operation is transcribed from villa 6e53201a in its order
 * (MSVC /fp:precise: no contraction, so the kernels are compiled without FMA
 * contraction and use fmaf exactly where villa calls std::fma):
 *  - canvas (r, c) = grid point (c + 0.5, r + 0.5): the bilinear of nodes
 *    (r..r+1, c..c+1) as (p00*.5 + p01*.5)*.5 + (p10*.5 + p11*.5)*.5
 *    (QuadSurface::gen, warpBilinearReplicateVec3f);
 *  - valid only where the nearest node (r + 1, c + 1) is valid (lround of
 *    c + 0.5; warpNearestConstU8); tifxyz -1 nodes are NaN;
 *  - normal = grid_normal_int at node (r + 1, c + 1): cross(right - left,
 *    down - up), times 1/sqrt(len2); NaN on the grid border; then negated
 *    (--flip-normals) and renormalised (normalizeNormals: v *= 1/sqrt(L2));
 *  - sample i at b + d * off_i, off_i = float((i - (N-1)/2) * step); trilinear
 *    with villa's fma chain, clamped to [0, 255] and truncated to uint8;
 *    outside the volume or in a missing chunk: 0 (ChunkSampler);
 *  - composite: nearest samples int(v + 0.5) at z = float(k) * step for
 *    k = start..end, their max (readCompositeFast, Nearest for composites).
 *
 * The volume is an uncompressed uint8 zarr v2 array (the job's raw.zarr/0) or,
 * without a .zarray, a folder of the grid's z#####_y#####_x#####.tif cubes
 * (cubes_RAW: 128 pages of 128 x 128, the zarr's chunks byte for byte; user
 * 2026-10-17, "make it run on the tiff format"); only the chunks the block's
 * samples touch are read. Fixed-size surface tiles share a bounded GPU cache;
 * tiles whose working set exceeds the cache are subdivided without resampling.
 * ==========================================================================*/

typedef struct {
    const char *volume;         /* zarr ARRAY dir (holds .zarray), else a folder of TIFF cubes */
    const char *segmentation;   /* tifxyz dir: x.tif y.tif z.tif meta.json (scale [1, 1]) */
    int         composite;      /* 0 layers -> zarr, 1 max composite -> tif */
    int         num_slices;
    double      slice_step;
    int         flip_normals;
    int         comp_start, comp_end;
    const char *zarr_output;    /* layers */
    const char *tif_output;     /* composite */
    double      voxel_size;     /* .zattrs only */
    const char *voxel_unit;
    int         zarr_shard_cols; /* 0: legacy v2; 1..64: v3, this many adjacent chunks per shard */
    size_t      gpu_cache_mb;    /* RAW chunk pool only; 0 selects the bounded 2048 MiB default */
} VcrRender_job;

/* Render one job on the visible CUDA device; logs progress to stdout.
 * 0 on success, nonzero with a message on stderr otherwise. */
int VcrRender_run(const VcrRender_job *job);

/* GPU kernels against the CPU transcription on a synthetic surface and
 * volume (both layer stacks and the composite must agree byte for byte);
 * SKIP (0) without a CUDA device. */
int VcrRender_selftest(void);

/* Copy a native uncompressed [N,128,128]-chunked v2 array through the output
 * writer, separating read/write time for storage qualification. Fresh output
 * group only; no GPU, resampling or intensity conversion. */
int VcrRender_storage_probe(const char *input_array, const char *output_group, int shard_cols);

#endif
