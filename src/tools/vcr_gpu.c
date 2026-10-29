/* vcr_gpu -- villa vc_render_tifxyz's surface-volume render on the GPU (src/flatten/vcr_render.h).
 *
 * A drop-in for the two vc_render_tifxyz calls of scripts/ink/region_ink_all.py: it takes the same options and
 * refuses any it does not reproduce, so a block is rendered exactly as the frozen vc_render would or not at all.
 *
 *   vcr_gpu --volume <zarr array | cubes_RAW folder> --scale 1 --group-idx 0 --segmentation <tifxyz> --num-slices N
 *           --slice-step 1
 *           [--flip-normals] --zarr-output <dir> [--pyramid 0] [--cache-gb G] [--voxel-size V --voxel-unit U]
 *   vcr_gpu ... --composite-collapse --accum-type max --composite-start=A --composite-end=B --tif-output <dir>
 *   vcr_gpu --selftest
 *
 * The card: device 0 of CUDA_VISIBLE_DEVICES (CUDA_DEVICE_ORDER=PCI_BUS_ID).  NVRTC (nvrtc64_*.dll and its
 * nvrtc-builtins) must sit beside the executable.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/cuda_drv.h"
#include "common/zarr_shard.h"
#include "flatten/vcr_render.h"

static int usage(const char *why)
{
    fprintf(stderr, "vcr_gpu: %s\nusage: vcr_gpu --volume V --scale 1 --group-idx 0 --segmentation S "
                    "(--num-slices N --zarr-output O | --composite-collapse --accum-type max --composite-start=A "
                    "--composite-end=B --tif-output D) [--slice-step 1] [--flip-normals] [--pyramid 0] "
                    "[--cache-gb G] [--voxel-size V] [--voxel-unit U] [--zarr-shard-cols 1..64] [--gpu-cache-mb M]\n"
                    "       vcr_gpu --selftest | --storage-selftest\n"
                    "       vcr_gpu --storage-probe INPUT_ARRAY FRESH_OUTPUT_GROUP SHARD_COLS\n", why);
    return 2;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--storage-selftest") == 0) return ZarrShard_selftest() ? 1 : 0;
    if (argc == 5 && strcmp(argv[1], "--storage-probe") == 0) {
        char *end = NULL;
        long count = strtol(argv[4], &end, 10);
        if (end == argv[4] || *end || count < 0 || count > ZARR_SHARD_MAX_CHUNKS) return usage("invalid shard columns");
        return VcrRender_storage_probe(argv[2], argv[3], (int)count);
    }
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        int fails = ZarrShard_selftest();
        fails += CudaDrv_selftest();
        fails += VcrRender_selftest();
        return fails ? 1 : 0;
    }
    VcrRender_job job;
    memset(&job, 0, sizeof job);
    job.slice_step = 1.0;
    const char *accum = NULL;
    double scale = 1.0;
    long group = 0, pyramid = 0;
    int have_start = 0, have_end = 0;
    for (int i = 1; i < argc; i++) {
        char key[64] = "";
        const char *val = NULL, *a = argv[i], *eq = strchr(a, '=');
        if (strncmp(a, "--", 2) != 0) return usage("arguments must be options");
        if (eq) {
            size_t n = (size_t)(eq - a) < sizeof key - 1 ? (size_t)(eq - a) : sizeof key - 1;
            memcpy(key, a, n);
            key[n] = '\0';
            val = eq + 1;
        } else {
            snprintf(key, sizeof key, "%s", a);
        }
        int is_flag = !strcmp(key, "--flip-normals") || !strcmp(key, "--composite-collapse");
        if (!is_flag && !val) {
            if (i + 1 >= argc) return usage("option without a value");
            val = argv[++i];
        }
        if (!strcmp(key, "--volume")) job.volume = val;
        else if (!strcmp(key, "--segmentation")) job.segmentation = val;
        else if (!strcmp(key, "--scale")) scale = atof(val);
        else if (!strcmp(key, "--group-idx")) group = atol(val);
        else if (!strcmp(key, "--num-slices")) job.num_slices = atoi(val);
        else if (!strcmp(key, "--slice-step")) job.slice_step = atof(val);
        else if (!strcmp(key, "--flip-normals")) job.flip_normals = 1;
        else if (!strcmp(key, "--zarr-output")) job.zarr_output = val;
        else if (!strcmp(key, "--zarr-shard-cols")) {
            char *end = NULL;
            long count = strtol(val, &end, 10);
            if (end == val || *end || count < 1 || count > ZARR_SHARD_MAX_CHUNKS) return usage("invalid shard columns");
            job.zarr_shard_cols = (int)count;
        }
        else if (!strcmp(key, "--tif-output")) job.tif_output = val;
        else if (!strcmp(key, "--composite-collapse")) job.composite = 1;
        else if (!strcmp(key, "--accum-type")) accum = val;
        else if (!strcmp(key, "--composite-start")) { job.comp_start = atoi(val); have_start = 1; }
        else if (!strcmp(key, "--composite-end")) { job.comp_end = atoi(val); have_end = 1; }
        else if (!strcmp(key, "--pyramid")) pyramid = atol(val);
        else if (!strcmp(key, "--cache-gb")) { /* CPU renderer compatibility; GPU has its own bounded cache. */ }
        else if (!strcmp(key, "--gpu-cache-mb")) {
            char *end = NULL;
            long count = strtol(val, &end, 10);
            if (end == val || *end || count < 16 || count > 65536) return usage("GPU cache must be 16..65536 MiB");
            job.gpu_cache_mb = (size_t)count;
        }
        else if (!strcmp(key, "--voxel-size")) job.voxel_size = atof(val);
        else if (!strcmp(key, "--voxel-unit")) job.voxel_unit = val;
        else {
            char why[160];
            snprintf(why, sizeof why, "option %s is not reproduced on the GPU (use vc_render_tifxyz)", key);
            return usage(why);
        }
    }
    if (!job.volume || !job.segmentation) return usage("--volume and --segmentation are required");
    if (scale != 1.0 || group != 0 || pyramid != 0)
        return usage("only --scale 1 --group-idx 0 --pyramid 0 renders are reproduced");
    if (!(job.slice_step > 0.0)) return usage("--slice-step must be positive");
    if (job.composite) {
        if (!accum || strcmp(accum, "max") != 0 || !have_start || !have_end || job.comp_end < job.comp_start
            || !job.tif_output || job.zarr_output || job.zarr_shard_cols)
            return usage("a composite needs --accum-type max, --composite-start/end and --tif-output only");
    } else {
        if (job.num_slices < 1 || !job.zarr_output || job.tif_output || accum)
            return usage("a layer stack needs --num-slices and --zarr-output only");
    }
    return VcrRender_run(&job);
}
