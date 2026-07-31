/*
 * main.c -- Vesuvius C Pipeline per-cube driver.
 *
 * Two input modes:
 *   1) TIFF:       ./cube_mesh input.tif output.tif [--dump-obj dir] [--no-qem]
 *                                                    [--no-timeout] [--halo N]
 *   2) stdin-raw:  ./cube_mesh --stdin-raw <p_size> <oz> <oy> <ox>
 *                              --dump-obj dir [--halo N] [--no-qem] [--no-timeout]
 *      Reads p_size^3 uint8 bytes (C order z,y,x) from stdin -- a padded cube
 *      whose index (0,0,0) is world voxel (oz-halo, oy-halo, ox-halo), exactly
 *      like HaloLoader_load. Lets a caller stream a cube straight from a remote
 *      zarr with no TIFF on disk. (oz,oy,ox) is the owned cube origin, so
 *      p_size must equal 128 + 2*halo. --dump-obj is required: the per-stage
 *      dumps are this mode's only output.
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
 * output.tif (TIFF mode) is a positional placeholder; this binary emits OBJs only.
 */
#include "common/ves_platform.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/arena.h"
#include "common/except.h"
#include "common/dump_obj.h"
#include "common/mls_project.h"
#include "pipeline/pipeline_cube.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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

static void usage(const char *argv0)
{
    fprintf(stderr,
        "Usage: %s input.tif output.tif [--dump-obj dir] [--dump-final-only] "
        "[--no-qem] [--simplify qem|cvt] [--no-timeout] [--halo N] [--trim-inset F]\n"
        "   or: %s --stdin-raw <p_size> <oz> <oy> <ox> --dump-obj dir "
        "[--dump-final-only] [--halo N] [--trim-inset F] [--no-qem] "
        "[--simplify qem|cvt] [--no-timeout]\n",
        argv0, argv0);
}

int main(int argc, char *argv[])
{
    int stdin_raw = (argc >= 2 && strcmp(argv[1], "--stdin-raw") == 0);

    const char *input_path  = NULL;
    const char *output_path = NULL;
    const char *dump_dir    = NULL;
    int skip_qem    = 0;
    int simplify_engine = 1;    /* 1 = CVT/RVD (default), 0 = QEM (--simplify qem) */
    int dump_final_only = 0;
    int no_timeout  = 0;
    int halo_voxels = 0;
    float trim_inset = -1.0f;   /* < 0 = BPA_OWNED_TRIM_INSET default */
    int p_size = 0, oz = 0, oy = 0, ox = 0;
    int opt_start = 0;

    if (stdin_raw) {
        if (argc < 6) { usage(argv[0]); return 1; }
        p_size = atoi(argv[2]);
        oz = atoi(argv[3]); oy = atoi(argv[4]); ox = atoi(argv[5]);
        opt_start = 6;
    } else {
        if (argc < 3) { usage(argv[0]); return 1; }
        input_path  = argv[1];
        output_path = argv[2];
        opt_start = 3;
    }

    for (int i = opt_start; i < argc; i++) {
        if (strcmp(argv[i], "--dump-obj") == 0 && i + 1 < argc) {
            dump_dir = argv[++i];
        } else if (strcmp(argv[i], "--no-qem") == 0) {
            skip_qem = 1;
        } else if (strcmp(argv[i], "--qem") == 0) {
            skip_qem = 0;
        } else if (strcmp(argv[i], "--simplify") == 0 && i + 1 < argc) {
            const char *e = argv[++i];
            if      (strcmp(e, "cvt") == 0) simplify_engine = 1;
            else if (strcmp(e, "qem") == 0) simplify_engine = 0;
            else { fprintf(stderr, "ERROR: --simplify must be qem|cvt\n"); return 1; }
        } else if (strcmp(argv[i], "--dump-final-only") == 0) {
            /* Write only step12_final (the grid_weld input); skip the
             * intermediate stage dumps. Fleet default via grid_pipeline —
             * the full stage set is ~300 MB/dense cube and IO-bounds
             * concurrent grid runs. Re-run one cube without this flag to
             * regenerate its stage dumps for seam debugging. */
            dump_final_only = 1;
        } else if (strcmp(argv[i], "--no-timeout") == 0) {
            no_timeout = 1;
        } else if (strcmp(argv[i], "--halo") == 0 && i + 1 < argc) {
            halo_voxels = atoi(argv[++i]);
            if (halo_voxels < 0) {
                fprintf(stderr, "ERROR: --halo must be >= 0\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--trim-inset") == 0 && i + 1 < argc) {
            trim_inset = (float)atof(argv[++i]);
            if (trim_inset < 0.0f || trim_inset > 8.0f) {
                fprintf(stderr, "ERROR: --trim-inset must be in [0, 8]\n");
                return 1;
            }
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    if (stdin_raw) {
        if (p_size <= 0 || oz < 0 || oy < 0 || ox < 0) {
            fprintf(stderr, "ERROR: bad --stdin-raw args\n");
            return 1;
        }
        if (p_size != 128 + 2 * halo_voxels) {
            fprintf(stderr, "ERROR: p_size (%d) must equal 128 + 2*halo (%d)\n",
                    p_size, 128 + 2 * halo_voxels);
            return 1;
        }
        if (!dump_dir) {
            /* Stage dumps are stdin mode's only output; a run without them
             * would leave nothing to inspect (house rule: always dump). */
            fprintf(stderr, "ERROR: --stdin-raw requires --dump-obj\n");
            return 1;
        }
        if (ves_stdin_set_binary() != 0) {
            fprintf(stderr, "ERROR: cannot set stdin to binary mode\n");
            return 1;
        }
    } else {
        if (ves_ensure_parent_dir(output_path) != 0) {
            fprintf(stderr, "ERROR: cannot create output directory for %s\n",
                    output_path);
            return 1;
        }
    }

    /* ALWAYS dump per-stage OBJs. A run with no --dump-obj would otherwise leave
     * no geometry to inspect (the output-path arg alone writes nothing) -- so if
     * --dump-obj was not given, default the dump dir to the output file's own
     * directory. Every cube_mesh run therefore emits its stage OBJs under
     * <out_dir>/<cube_id>/<cube_id>_<stage>/. grid_pipeline passes --dump-obj
     * explicitly, so this only changes bare cube_mesh invocations. (stdin-raw
     * mode has no output path to derive from; it requires --dump-obj above.) */
    char default_dump[1024] = {0};
    if (!dump_dir) {
        strncpy(default_dump, output_path, sizeof(default_dump) - 1);
        char *ls = NULL;
        for (char *p = default_dump; *p; p++)
            if (*p == '/' || *p == '\\') ls = p;
        if (ls) *ls = '\0';
        else strcpy(default_dump, ".");
        dump_dir = default_dump;
        fprintf(stderr,
            "  note: no --dump-obj given; defaulting per-stage dumps to %s\n",
            dump_dir);
    }

    int n_threads = get_thread_count();
#ifdef _OPENMP
    /* Bound OpenMP parallel regions (the MLS per-vertex loop) by the same
     * budget as everything else: VESUVIUS_THREADS. grid orchestrators pass
     * threads-per-cube (usually 1 -- the fleet fills the cores); a bare
     * single-cube run gets the whole machine. */
    omp_set_num_threads(n_threads);
#endif

    char cube_id[128] = {0};
    if (stdin_raw) {
        snprintf(cube_id, sizeof(cube_id), "z%05d_y%05d_x%05d", oz, oy, ox);
    } else {
        DumpObj_extract_cube_id(input_path, cube_id, sizeof(cube_id));
    }

    /* pred_dir = directory of input_path (TIFF+halo neighbor lookup only). */
    char pred_dir_buf[1024] = {0};
    const char *pred_dir = NULL;
    if (!stdin_raw && halo_voxels > 0) {
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

    fprintf(stderr, "cube_mesh: %s (threads=%d, halo=%d)\n",
            stdin_raw ? cube_id : input_path, n_threads, halo_voxels);

    if (!MLS_cubecl_preflight()) {
        fprintf(stderr, "cube_mesh: CubeCL MLS preflight failed: %s\n",
                MLS_cubecl_last_error());
        return 5;
    }
#ifdef VESUVIUS_MLS_CUBECL
    {
        const char *backend = getenv("MLS_BACKEND");
        fprintf(stderr, "cube_mesh: CubeCL MLS preflight OK (backend=%s)\n",
                (backend && backend[0]) ? backend : "auto");
    }
#endif

    double t_total = ves_clock_sec();
    double hard_timeout = no_timeout ? 1e9 : HARD_TIMEOUT_SEC;
    ves_hard_timeout_start(hard_timeout + 1.0, &g_timeout_flag);

    Arena_T arena = Arena_new();
    int pipeline_ok = 0;
    PipelineOutput out = {0};

    TRY
        uint8_t *raw_buf = NULL;
        int input_ok = 1;
        if (stdin_raw) {
            size_t n = (size_t)p_size * (size_t)p_size * (size_t)p_size;
            raw_buf = (uint8_t *)ARENA_CALLOC(arena, (long)n, 1L);
            size_t got = 0, r = 0;
            while (got < n && (r = fread(raw_buf + got, 1, n - got, stdin)) > 0)
                got += r;
            if (got != n) {
                fprintf(stderr, "stdin-raw: expected %zu bytes, got %zu\n", n, got);
                input_ok = 0;
            }
        }

        if (input_ok) {
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
                .trim_inset       = trim_inset,
                .dump_dir         = dump_dir,
                .skip_qem         = skip_qem,
                .simplify_engine  = simplify_engine,
                .dump_final_only  = dump_final_only,
                .vol_in           = raw_buf,
                .p_size_in        = stdin_raw ? p_size : 0,
                .cube_origin_zyx  = { oz, oy, ox },
            };
            int rc = pipeline_process_cube(arena, &in, &out);
            if (rc == 0) pipeline_ok = 1;
        }

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
