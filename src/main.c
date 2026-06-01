/*
 * main.c -- Vesuvius C Pipeline per-cube driver.
 *
 * Usage: ./cube_mesh input.tif output.tif [--dump-obj dir] [--no-qem]
 *                                          [--no-timeout] [--halo N]
 *
 * Body of the per-cube pipeline lives in src/pipeline/pipeline_cube.c.
 * This file is just: arg parsing -> single pipeline_process_cube() call
 * -> exit status. Driven per-cube by the grid_pipeline orchestrator.
 *
 * --dump-obj dir:  write OBJ meshes under dir/<cube_id>/<cube_id>_<stage>/
 * --halo N:        load N voxels of safety boundary from neighbor cubes
 * --no-qem:        skip per-component decimation (raw MC+LOP only)
 * --no-timeout:    disable the hard wall-clock timeout
 *
 * output.tif is a positional placeholder; this binary emits OBJs only.
 */
#include "common/ves_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/arena.h"
#include "common/except.h"
#include "common/dump_obj.h"
#include "pipeline/pipeline_cube.h"

#define HARD_TIMEOUT_SEC   1e9     /* effectively infinite */

static volatile sig_atomic_t g_timeout_flag = 0;

static int get_thread_count(void)
{
    const char *env = getenv("VESUVIUS_THREADS");
    if (env) {
        int n = atoi(env);
        if (n > 0) return n;
    }
    int n = ves_cpu_count();
    return (n > 0) ? n : 2;
}

int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s input.tif output.tif [--dump-obj dir] "
                "[--no-qem] [--no-timeout] [--halo N]\n", argv[0]);
        return 1;
    }
    const char *input_path  = argv[1];
    const char *output_path = argv[2];
    const char *dump_dir    = NULL;
    int skip_qem            = 0;
    int no_timeout          = 0;
    int halo_voxels         = 0;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--dump-obj") == 0 && i + 1 < argc) {
            dump_dir = argv[++i];
        } else if (strcmp(argv[i], "--no-qem") == 0) {
            skip_qem = 1;
        } else if (strcmp(argv[i], "--qem") == 0) {
            skip_qem = 0;
        } else if (strcmp(argv[i], "--no-timeout") == 0) {
            no_timeout = 1;
        } else if (strcmp(argv[i], "--halo") == 0 && i + 1 < argc) {
            halo_voxels = atoi(argv[++i]);
            if (halo_voxels < 0) {
                fprintf(stderr, "ERROR: --halo must be >= 0\n");
                return 1;
            }
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    if (ves_ensure_parent_dir(output_path) != 0) {
        fprintf(stderr, "ERROR: cannot create output directory for %s\n",
                output_path);
        return 1;
    }

    int n_threads = get_thread_count();

    char cube_id[128] = {0};
    DumpObj_extract_cube_id(input_path, cube_id, sizeof(cube_id));

    /* Derive pred_dir = directory of input_path. */
    char pred_dir_buf[1024] = {0};
    const char *pred_dir = NULL;
    if (halo_voxels > 0) {
        strncpy(pred_dir_buf, input_path, sizeof(pred_dir_buf) - 1);
        char *last_slash = NULL;
        for (char *p = pred_dir_buf; *p; p++) {
            if (*p == '/' || *p == '\\') last_slash = p;
        }
        if (last_slash) {
            *last_slash = '\0';
            pred_dir = pred_dir_buf;
        } else {
            pred_dir = ".";
        }
    }

    fprintf(stderr, "cube_mesh: %s -> %s  (threads=%d, halo=%d)\n",
            input_path, output_path, n_threads, halo_voxels);

    double t_total = ves_clock_sec();
    double hard_timeout = no_timeout ? 1e9 : HARD_TIMEOUT_SEC;
    ves_hard_timeout_start(hard_timeout + 1.0, &g_timeout_flag);

    Arena_T arena = Arena_new();
    int pipeline_ok = 0;
    PipelineOutput out = {0};

    TRY
        PipelineInput in = {
            .tiff_path        = input_path,
            .pred_dir         = pred_dir,
            .cube_id          = cube_id,
            .halo_voxels      = halo_voxels,
            .cube_D           = 128,
            .cube_H           = 128,
            .cube_W           = 128,
            .n_threads        = n_threads,
            .qem_target_ratio = 0.0f,  /* use default */
            .dump_dir         = dump_dir,
            .skip_qem         = skip_qem,
        };
        int rc = pipeline_process_cube(arena, &in, &out);
        if (rc == 0) pipeline_ok = 1;

    EXCEPT(Timeout)
        fprintf(stderr, "TIMEOUT after %.1fs\n",
                ves_clock_sec() - t_total);
    EXCEPT(Arena_Failed)
        fprintf(stderr, "ARENA OOM\n");
    EXCEPT(IO_Failed)
        fprintf(stderr, "I/O FAILURE\n");
    END_TRY;

    ves_hard_timeout_cancel();

    double total_time = ves_clock_sec() - t_total;
    fprintf(stderr,
        "  Timings: extract=%.3f qem=%.3f trim=%.3f dump=%.3f total=%.3fs"
        "  Status: %s\n",
        out.t_extract, out.t_qem, out.t_trim, out.t_dump, total_time,
        pipeline_ok ? "OK" : "FAILED");

    Arena_dispose(&arena);
    return pipeline_ok ? 0 : 1;
}
