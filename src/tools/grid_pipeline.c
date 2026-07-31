/*
 * grid_pipeline.c -- cube-level orchestrator.
 *
 * Spawns cube_mesh.exe subprocesses, N at a time, over every cube in a
 * grid, then invokes grid_weld.exe to stitch the per-cube OBJs into one
 * welded mesh. Replaces scripts/run_grid_halo.ps1 with portable C and
 * eliminates the serial PowerShell foreach.
 *
 * Usage:
 *   grid_pipeline <grid_dir> <output_dir>
 *                 [--halo N] [--threads-per-cube N] [--max-concurrent N]
 *                 [--exe path] [--weld path]
 *                 [--max-cubes N] [--skip-weld] [--check-determinism]
 *                 [--simplify qem|cvt] [--qem | --no-qem]
 *                 [--umb-y F --umb-x F --wrap-pitch F]
 *                 [--grow-wind-tol F]
 *
 * Layout:
 *   <grid_dir>/cubes_PRED/*.tif         <- per-cube input prediction TIFFs
 *   <output_dir>/dump/<cube_id>/...     <- per-cube OBJ dumps
 *   <output_dir>/logs/<cube_id>.log     <- per-cube stderr/stdout
 *   <output_dir>/welded.obj             <- final welded mesh
 *   <output_dir>/pipeline_summary.csv   <- per-cube exit code, timing
 *
 * Concurrency: --max-concurrent caps how many cube_mesh.exe subprocesses
 * run at once. Each subprocess uses VESUVIUS_THREADS=--threads-per-cube
 * OpenMP threads. Default total = ves_cpu_count() / threads_per_cube.
 */
#include "../common/ves_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../extract/pred_reject.h"

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
  #include <sys/wait.h>
  #include <unistd.h>
  #include <fcntl.h>
#endif

/* Sentinel exit code recorded in pipeline_summary.csv for a cube skipped as a
 * garbage (solid-slab) prediction. Distinct from 0 (ok) and the cube_mesh
 * failure codes so the CSV stays self-documenting. */
#define GP_REJECT_CODE (-2)

typedef struct {
    char        cube_id[128];
    char        tiff_path[1024];
    char        log_path[1024];
    int         exit_code;
    double      wall_seconds;
} CubeJob;

/* ---- Argument plumbing ---- */

typedef struct {
    const char *grid_dir;
    const char *output_dir;
    int         halo;
    int         threads_per_cube;
    int         max_concurrent;
    const char *exe_path;
    const char *weld_path;
    int         max_cubes;
    int         skip_weld;
    int         check_determinism;
    int         skip_qem;
    int         simplify_engine; /* 1 = CVT/RVD (default), 0 = QEM (--simplify qem) */
    int         skip_existing;  /* resume: skip cubes whose step12_final OBJ is already complete */
    int         dry_run;        /* report skip/run decisions and exit; spawn nothing */
    int         reject_garbage; /* gate garbage (solid-slab) cubes pre-spawn (default 1) */
    const char *reject_list;    /* optional: skip cube_ids listed in this file (else detect inline) */
    int         full_dumps;     /* pass NO --dump-final-only to children: every cube
                                 * writes all intermediate stage OBJs (~300 MB/dense
                                 * cube -- IO-bounds concurrent runs; debug only). */
    double      umb_y, umb_x;   /* scroll umbilicus (source-space vox, y/x) */
    double      wrap_pitch;     /* wrap pitch (vox/turn). These also arm
                                 * grid_weld's seam winding/phase gates; the
                                 * Jul-12 rebuild harness dropped them and paid
                                 * 661 seam-band gaps against 16 when armed. */
    double      grow_wind_tol;  /* BPA branch growth half-width in turns; 0=off */
    int         have_umb_y, have_umb_x, have_wrap_pitch;
                                /* geometry arguments arm grid_weld's seam
                                 * winding/phase gates. Per-cube BPA gating is
                                 * separately opt-in via --grow-wind-tol. */
    int         have_grow_wind_tol;
    float       trim_inset;     /* owned-box inset passthrough. Default 0: the
                                 * whole-grid unwrap path is the documented
                                 * chain, and it pairs cross-seam skins, so
                                 * charts have to REACH the cube faces. The
                                 * 1-vox inset measurably starves that: on
                                 * PHerc0139-4x5x5 it halves the seam evidence
                                 * (2,713 pairs vs 6,441) and the registration
                                 * audit fails at 15.78% whole-turn error
                                 * instead of passing at 2.98%. Pass < 0 for
                                 * cube_mesh's own 1.0 default, which suits a
                                 * weld-only run where the bridge wants a gap. */
} GpOptions;

static void set_process_env_double(const char *name, double value)
{
    char text[64];
    snprintf(text, sizeof(text), "%.12g", value);
#ifdef _WIN32
    _putenv_s(name, text);
#else
    setenv(name, text, 1);
#endif
}

static void clear_process_env(const char *name)
{
#ifdef _WIN32
    _putenv_s(name, "");
#else
    unsetenv(name);
#endif
}


static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s <grid_dir> <output_dir> [options]\n"
        "Options:\n"
        "  --halo N                  Halo voxels (default 13 = MLS R+1)\n"
        "  --threads-per-cube N      OpenMP threads per cube (default 1)\n"
        "  --max-concurrent N        Max concurrent subprocesses (default cores/tpc)\n"
        "  --exe PATH                cube_mesh.exe path (default build/Release/cube_mesh.exe)\n"
        "  --weld PATH               grid_weld.exe path (default build/Release/grid_weld.exe)\n"
        "  --max-cubes N             Stop after N cubes (default all)\n"
        "  --skip-weld               Do not run grid_weld at end\n"
        "  --skip-existing           Resume: skip cubes whose step12_final OBJ\n"
        "                            already exists and looks complete\n"
        "  --dry-run                 With --skip-existing: print which cubes would\n"
        "                            be skipped vs run, then exit (spawns nothing)\n"
        "  --check-determinism       Run each cube twice and cmp OBJs\n"
        "  --simplify qem|cvt        Simplifier (default cvt)\n"
        "  --qem                     Alias for --simplify qem\n"
        "  --no-qem                  Pass --no-qem to cube_mesh\n"
        "  --no-reject-garbage       Process all cubes (do not skip solid-slab garbage)\n"
        "  --reject-list FILE        Skip cube_ids listed in FILE (one per line);\n"
        "                            default detects garbage inline from each TIFF\n"
        "  --trim-inset F            Owned-box inset passthrough to cube_mesh\n"
        "                            (default 0: charts reach the cube faces so\n"
        "                            the whole-grid unwrap can pair cross-seam\n"
        "                            skins; pass 1 for a weld-only run, where\n"
        "                            the bridge wants a gap)\n"
        "  --full-dumps              Children write ALL intermediate stage OBJs\n"
        "                            (default: step12_final only; the full set is\n"
        "                            ~300 MB/dense cube and IO-bounds the fleet)\n"
        "  --umb-y F --umb-x F       Scroll umbilicus (source-space voxels) and\n"
        "  --wrap-pitch F            wrap pitch (vox/turn): arms grid_weld's seam\n"
        "                            winding/phase gates, not per-cube meshing.\n"
        "  --grow-wind-tol F         EXPERIMENTAL opt-in per-cube BPA growth gate;\n"
        "                            off by default, 0 disables, must be < .5.\n"
        "                            This can fragment sheets; all three geometry\n"
        "                            arguments above are required when F > 0.\n"
        "                            PHerc0139 weld geometry:\n"
        "                            --umb-y 3405 --umb-x 2878 --wrap-pitch 9.5\n"
        "  --selftest                Run built-in unit tests and exit\n",
        prog);
}

static int parse_args(int argc, char *argv[], GpOptions *o)
{
    memset(o, 0, sizeof(*o));
    o->halo = 13;   /* >= MLS_PROJECT_RADIUS_VOX+1 so boundary LOP is fully
                     * two-sided supported and adjacent cubes agree at the seam */
    o->threads_per_cube = 1;
    o->max_concurrent = 0;  /* derived below */
    o->max_cubes = 0;
    o->simplify_engine = 1; /* CVT/RVD by default; --simplify qem selects the old QEM path */
    o->exe_path = "build/Release/cube_mesh.exe";
    o->weld_path = "build/Release/grid_weld.exe";
    o->reject_garbage = 1;  /* gate solid-slab garbage cubes by default */
    o->trim_inset = 0.0f;   /* whole-scroll path: charts reach the cube faces */
    o->grow_wind_tol = 0.0;     /* fail-safe: coherent sheets are the default */

    if (argc < 3) { usage(argv[0]); return -1; }
    o->grid_dir = argv[1];
    o->output_dir = argv[2];

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--halo") && i + 1 < argc) {
            o->halo = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--threads-per-cube") && i + 1 < argc) {
            o->threads_per_cube = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--max-concurrent") && i + 1 < argc) {
            o->max_concurrent = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--exe") && i + 1 < argc) {
            o->exe_path = argv[++i];
        } else if (!strcmp(argv[i], "--weld") && i + 1 < argc) {
            o->weld_path = argv[++i];
        } else if (!strcmp(argv[i], "--max-cubes") && i + 1 < argc) {
            o->max_cubes = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--full-dumps")) {
            o->full_dumps = 1;
        } else if (!strcmp(argv[i], "--umb-y") && i + 1 < argc) {
            o->umb_y = atof(argv[++i]);
            o->have_umb_y = 1;
        } else if (!strcmp(argv[i], "--umb-x") && i + 1 < argc) {
            o->umb_x = atof(argv[++i]);
            o->have_umb_x = 1;
        } else if (!strcmp(argv[i], "--wrap-pitch") && i + 1 < argc) {
            o->wrap_pitch = atof(argv[++i]);
            o->have_wrap_pitch = 1;
        } else if (!strcmp(argv[i], "--grow-wind-tol") && i + 1 < argc) {
            o->grow_wind_tol = atof(argv[++i]);
            o->have_grow_wind_tol = 1;
        } else if (!strcmp(argv[i], "--skip-weld")) {
            o->skip_weld = 1;
        } else if (!strcmp(argv[i], "--skip-existing") ||
                   !strcmp(argv[i], "--resume")) {
            o->skip_existing = 1;
        } else if (!strcmp(argv[i], "--dry-run")) {
            o->dry_run = 1;
        } else if (!strcmp(argv[i], "--check-determinism")) {
            o->check_determinism = 1;
        } else if (!strcmp(argv[i], "--no-qem")) {
            o->skip_qem = 1;
        } else if (!strcmp(argv[i], "--qem")) {
            o->skip_qem = 0;
            o->simplify_engine = 0;
        } else if (!strcmp(argv[i], "--simplify") && i + 1 < argc) {
            const char *e = argv[++i];
            if      (!strcmp(e, "cvt")) o->simplify_engine = 1;
            else if (!strcmp(e, "qem")) o->simplify_engine = 0;
            else { fprintf(stderr, "ERROR: --simplify must be qem|cvt\n"); return -1; }
        } else if (!strcmp(argv[i], "--no-reject-garbage")) {
            o->reject_garbage = 0;
        } else if (!strcmp(argv[i], "--reject-list") && i + 1 < argc) {
            o->reject_list = argv[++i];
        } else if (!strcmp(argv[i], "--trim-inset") && i + 1 < argc) {
            o->trim_inset = (float)atof(argv[++i]);
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return -1;
        }
    }

    if (o->grow_wind_tol < 0.0 || o->grow_wind_tol >= 0.5) {
        fprintf(stderr, "ERROR: --grow-wind-tol must be in [0, .5)\n");
        return -1;
    }
    if ((o->have_umb_y || o->have_umb_x || o->have_wrap_pitch)
        && !(o->have_umb_y && o->have_umb_x && o->have_wrap_pitch)) {
        fprintf(stderr,
                "ERROR: --umb-y, --umb-x and --wrap-pitch must be supplied together\n");
        return -1;
    }
    if (o->have_wrap_pitch && o->wrap_pitch <= 0.0) {
        fprintf(stderr, "ERROR: --wrap-pitch must be positive\n");
        return -1;
    }
    if (o->have_grow_wind_tol && o->grow_wind_tol > 0.0
        && !(o->have_umb_y && o->have_umb_x && o->have_wrap_pitch)) {
        fprintf(stderr,
                "ERROR: positive --grow-wind-tol requires --umb-y, --umb-x and --wrap-pitch\n");
        return -1;
    }

    if (o->max_concurrent == 0) {
        int cores = ves_cpu_count();
        o->max_concurrent = cores / o->threads_per_cube;
        if (o->max_concurrent < 1) o->max_concurrent = 1;
    }

    return 0;
}

/* ---- Scan <grid_dir>/cubes_PRED/ for *.tif files ---- */

static int scan_cubes(const char *grid_dir,
                      CubeJob **out_jobs, size_t *out_n)
{
    char pred_dir[1024];
    snprintf(pred_dir, sizeof(pred_dir), "%s/cubes_PRED", grid_dir);

    size_t cap = 256;
    size_t n = 0;
    CubeJob *jobs = (CubeJob *)calloc(cap, sizeof(CubeJob));
    if (!jobs) return -1;

#ifdef _WIN32
    char glob[1024];
    snprintf(glob, sizeof(glob), "%s/*.tif", pred_dir);
    WIN32_FIND_DATAA find_data;
    HANDLE hfind = FindFirstFileA(glob, &find_data);
    if (hfind == INVALID_HANDLE_VALUE) {
        free(jobs);
        return -1;
    }
    do {
        if (find_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
        const char *name = find_data.cFileName;
        size_t nlen = strlen(name);
        if (nlen < 5 || strcmp(name + nlen - 4, ".tif") != 0) continue;
        if (n >= cap) {
            cap *= 2;
            jobs = (CubeJob *)realloc(jobs, cap * sizeof(CubeJob));
            if (!jobs) { FindClose(hfind); return -1; }
        }
        memset(&jobs[n], 0, sizeof(CubeJob));
        size_t id_len = nlen - 4;
        if (id_len >= sizeof(jobs[n].cube_id)) id_len = sizeof(jobs[n].cube_id) - 1;
        memcpy(jobs[n].cube_id, name, id_len);
        snprintf(jobs[n].tiff_path, sizeof(jobs[n].tiff_path),
                 "%s/%s", pred_dir, name);
        n++;
    } while (FindNextFileA(hfind, &find_data));
    FindClose(hfind);
#else
    DIR *d = opendir(pred_dir);
    if (!d) { free(jobs); return -1; }
    struct dirent *de = NULL;
    while ((de = readdir(d)) != NULL) {
        size_t nlen = strlen(de->d_name);
        if (nlen < 5 || strcmp(de->d_name + nlen - 4, ".tif") != 0) continue;
        if (n >= cap) {
            cap *= 2;
            jobs = (CubeJob *)realloc(jobs, cap * sizeof(CubeJob));
            if (!jobs) { closedir(d); return -1; }
        }
        memset(&jobs[n], 0, sizeof(CubeJob));
        size_t id_len = nlen - 4;
        if (id_len >= sizeof(jobs[n].cube_id)) id_len = sizeof(jobs[n].cube_id) - 1;
        memcpy(jobs[n].cube_id, de->d_name, id_len);
        snprintf(jobs[n].tiff_path, sizeof(jobs[n].tiff_path),
                 "%s/%s", pred_dir, de->d_name);
        n++;
    }
    closedir(d);
#endif

    /* Sort by cube_id for deterministic ordering. */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (strcmp(jobs[i].cube_id, jobs[j].cube_id) > 0) {
                CubeJob tmp = jobs[i];
                jobs[i] = jobs[j];
                jobs[j] = tmp;
            }
        }
    }

    *out_jobs = jobs;
    *out_n = n;
    return 0;
}

/* ---- Spawn cube_mesh on one cube; redirect output to log_path ---- */

static int run_one_cube(const char *exe_path,
                        const CubeJob *job,
                        const char *output_dir,
                        const char *dump_dir,
                        int halo, int threads, int skip_qem,
                        float trim_inset, int dump_final_only,
                        int simplify_engine)
{
    char out_tif[1024];
    snprintf(out_tif, sizeof(out_tif), "%s/cubes/%s.tif",
             output_dir, job->cube_id);

    /* Compose argv. */
    char halo_str[32];
    snprintf(halo_str, sizeof(halo_str), "%d", halo);
    char trim_str[32];
    snprintf(trim_str, sizeof(trim_str), "%.3f", (double)trim_inset);

    const char *argv[20];
    int argc = 0;
    argv[argc++] = exe_path;
    argv[argc++] = job->tiff_path;
    argv[argc++] = out_tif;
    argv[argc++] = "--halo";
    argv[argc++] = halo_str;
    argv[argc++] = "--dump-obj";
    argv[argc++] = dump_dir;
    argv[argc++] = "--no-timeout";
    if (dump_final_only) argv[argc++] = "--dump-final-only";
    if (skip_qem) argv[argc++] = "--no-qem";
    if (simplify_engine == 1) { argv[argc++] = "--simplify"; argv[argc++] = "cvt"; }
    if (trim_inset >= 0.0f) {
        argv[argc++] = "--trim-inset";
        argv[argc++] = trim_str;
    }
    argv[argc] = NULL;

#ifdef _WIN32
    /* Build cmdline manually for CreateProcess with file redirection. */
    char cmdline[4096];
    size_t pos = 0;
    for (int i = 0; argv[i] != NULL; i++) {
        if (i > 0 && pos < sizeof(cmdline) - 1) cmdline[pos++] = ' ';
        if (pos < sizeof(cmdline) - 1) cmdline[pos++] = '"';
        size_t alen = strlen(argv[i]);
        if (pos + alen < sizeof(cmdline) - 2) {
            memcpy(cmdline + pos, argv[i], alen);
            pos += alen;
        }
        if (pos < sizeof(cmdline) - 1) cmdline[pos++] = '"';
    }
    cmdline[pos] = '\0';

    SECURITY_ATTRIBUTES sa;
    sa.nLength = sizeof(sa);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;
    HANDLE hLog = CreateFileA(job->log_path, GENERIC_WRITE,
                              FILE_SHARE_READ | FILE_SHARE_WRITE,
                              &sa, CREATE_ALWAYS,
                              FILE_ATTRIBUTE_NORMAL, NULL);
    if (hLog == INVALID_HANDLE_VALUE) {
        fprintf(stderr, "[%s] cannot open log %s\n",
                job->cube_id, job->log_path);
        return -1;
    }

    /* VESUVIUS_THREADS reaches the child by inheritance: the orchestrator
     * sets it ONCE (SetEnvironmentVariableA in main, before the spawn
     * threads exist) and every child gets it via lpEnvironment=NULL.
     * A per-child ANSI environment block was tried and abandoned — the
     * CreateProcessA ANSI-env conversion path rejects otherwise-valid
     * blocks with ERROR_INVALID_PARAMETER, and per-child blocks buy
     * nothing when every cube shares the same thread budget anyway. */
    (void)threads;

    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    memset(&si, 0, sizeof(si));
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = hLog;
    si.hStdError  = hLog;
    si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
    memset(&pi, 0, sizeof(pi));

    BOOL ok = CreateProcessA(exe_path, cmdline, NULL, NULL, TRUE,
                              0, NULL, NULL, &si, &pi);
    CloseHandle(hLog);
    if (!ok) {
        DWORD err = GetLastError();
        fprintf(stderr, "[%s] CreateProcess failed: %lu cmdline=%s\n",
                job->cube_id, (unsigned long)err, cmdline);
        return -1;
    }

    WaitForSingleObject(pi.hProcess, INFINITE);
    DWORD exit_code = 1;
    GetExitCodeProcess(pi.hProcess, &exit_code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    return (int)exit_code;

#else
    /* POSIX: fork + execv with stdout/stderr -> log file. */
    int fd = open(job->log_path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) return -1;

    pid_t pid = fork();
    if (pid < 0) { close(fd); return -1; }
    if (pid == 0) {
        char env_var[64];
        snprintf(env_var, sizeof(env_var), "VESUVIUS_THREADS=%d", threads);
        putenv(env_var);
        dup2(fd, STDOUT_FILENO);
        dup2(fd, STDERR_FILENO);
        close(fd);
        execv(exe_path, (char *const *)argv);
        _exit(127);
    }
    close(fd);
    int status = 0;
    waitpid(pid, &status, 0);
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return -1;
#endif
}

/* ---- Main ---- */

static int ensure_dir(const char *path)
{
    /* ves_ensure_parent_dir treats path as a file; we want the path itself. */
    char fake_child[1024];
    snprintf(fake_child, sizeof(fake_child), "%s/x", path);
    return ves_ensure_parent_dir(fake_child);
}

static int absolute_path(const char *path, char *out, size_t out_size)
{
#ifdef _WIN32
    DWORD n = GetFullPathNameA(path, (DWORD)out_size, out, NULL);
    return n > 0 && n < (DWORD)out_size ? 0 : -1;
#else
    char *resolved = realpath(path, NULL);
    if (!resolved) return -1;
    size_t n = strlen(resolved);
    if (n >= out_size) { free(resolved); return -1; }
    memcpy(out, resolved, n + 1);
    free(resolved);
    return 0;
#endif
}

static void write_json_string(FILE *f, const char *text)
{
    fputc('"', f);
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
        if (*p == '"') fputs("\\\"", f);
        else if (*p == '\\') fputs("\\\\", f);
        else if (*p == '\n') fputs("\\n", f);
        else if (*p == '\r') fputs("\\r", f);
        else if (*p == '\t') fputs("\\t", f);
        else if (*p >= 0x20) fputc(*p, f);
    }
    fputc('"', f);
}

/* Write this before any cube is launched. scrollslice walks upward from any
 * output OBJ until it finds this file, then recovers raw.zarr, pred.zarr, and
 * origin_zyx through the grid manifest. The relative default mesh becomes
 * usable when the final weld lands; intermediate meshes override it. */
static int write_scrollslice_source(const GpOptions *opts)
{
    char path[1024];
    char dataset[2048];
    int n = snprintf(path, sizeof(path),
                     "%s/scrollslice.source.json", opts->output_dir);
    if (n < 0 || (size_t)n >= sizeof(path)
        || absolute_path(opts->grid_dir, dataset, sizeof(dataset)) != 0) {
        return -1;
    }
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    fputs("{\n  \"format\": \"vesuvius-scrollslice-source-v1\",\n"
          "  \"dataset\": ", f);
    write_json_string(f, dataset);
    fputs(",\n  \"mesh_axes\": \"zyx\",\n"
          "  \"mesh\": \"welded.obj\"\n}\n", f);
    if (fclose(f) != 0) return -1;
    fprintf(stderr, "scrollslice source: %s\n", path);
    return 0;
}

/* Build the path to a cube's final per-cube mesh (the stage grid_weld reads). */
static void step12_obj_path(char *buf, size_t n,
                            const char *dump_dir, const char *cube_id)
{
    snprintf(buf, n, "%s/%s/%s_step12_final/%s_step12_final_all.obj",
             dump_dir, cube_id, cube_id, cube_id);
}

/* A cube counts as "already done" (for --skip-existing resume) iff its final
 * per-cube OBJ exists and looks complete: a non-trivial size and a clean
 * trailing newline. A cube_mesh process killed mid-dump leaves either no file
 * or a truncated tail whose last byte is not '\n', so it is re-run. The size
 * floor rejects empty/header-only stubs. This is intentionally a cheap O(1)
 * stat+tail read so a full-grid scan stays fast. */
static int obj_looks_complete(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return 0; }
    long sz = ftell(f);
    if (sz < 128) { fclose(f); return 0; }   /* header + a few verts minimum */
    if (fseek(f, -1, SEEK_END) != 0) { fclose(f); return 0; }
    int last = fgetc(f);
    fclose(f);
    return last == '\n';
}

/* ---- Optional reject-list (cube_ids to skip), loaded once before the run. ---- */
typedef struct {
    char (*ids)[128];
    size_t n;
} IdSet;

static int idset_load(const char *path, IdSet *s)
{
    s->ids = NULL; s->n = 0;
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    size_t cap = 64, n = 0;
    char (*ids)[128] = (char (*)[128])malloc(cap * sizeof(*ids));
    if (!ids) { fclose(f); return -1; }
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        size_t L = strlen(line);
        while (L > 0 && (line[L-1] == '\n' || line[L-1] == '\r' ||
                         line[L-1] == ' '  || line[L-1] == '\t')) line[--L] = '\0';
        if (L == 0) continue;
        if (n >= cap) {
            cap *= 2;
            char (*g)[128] = (char (*)[128])realloc(ids, cap * sizeof(*ids));
            if (!g) { free(ids); fclose(f); return -1; }
            ids = g;
        }
        if (L >= sizeof(ids[0])) L = sizeof(ids[0]) - 1;
        memcpy(ids[n], line, L);
        ids[n][L] = '\0';
        n++;
    }
    fclose(f);
    s->ids = ids; s->n = n;
    return 0;
}

static int idset_has(const IdSet *s, const char *id)
{
    for (size_t i = 0; i < s->n; i++)
        if (strcmp(s->ids[i], id) == 0) return 1;
    return 0;
}

static void idset_free(IdSet *s) { free(s->ids); s->ids = NULL; s->n = 0; }

/* Cheap garbage check: load a cube's OWN 128^3 prediction (no halo -- the
 * verdict is per-cube) and run the garbage detector. Per-call arena so it is
 * thread-safe inside the OpenMP loop. Returns the reject KIND:
 *   0 = keep, 1 = solid-slab garbage, 2 = empty garbage (no meshable component). */
static int cube_reject_kind(const char *tiff_path)
{
    Arena_T a = Arena_new();
    uint8_t *vol = NULL;
    int D = 0, H = 0, W = 0;
    int kind = 0;
    if (TiffIO_load(a, tiff_path, &vol, &D, &H, &W) == 0) {
        PredRejectStats st;
        if (PredReject_is_garbage(a, vol, D, H, W, &st))
            kind = st.empty ? 2 : 1;
    }
    Arena_dispose(&a);
    return kind;
}

/* Built-in unit test for the resume completeness check (--selftest). */
static int run_selftest(void)
{
    const char *p_ok    = "gp_selftest_ok.obj";
    const char *p_trunc = "gp_selftest_trunc.obj";
    const char *p_tiny  = "gp_selftest_tiny.obj";
    FILE *f = NULL;

    f = fopen(p_ok, "wb");      /* complete: verts + face + trailing newline */
    if (f) { for (int i = 0; i < 24; i++) fprintf(f, "v %d.0 %d.0 %d.0\n", i, i, i);
             fprintf(f, "f 1 2 3\n"); fclose(f); }
    f = fopen(p_trunc, "wb");   /* same bytes but killed mid-line: no final \n */
    if (f) { for (int i = 0; i < 24; i++) fprintf(f, "v %d.0 %d.0 %d.0\n", i, i, i);
             fprintf(f, "f 1 2 3"); fclose(f); }
    f = fopen(p_tiny, "wb");    /* below the size floor */
    if (f) { fprintf(f, "v 0 0 0\n"); fclose(f); }

    struct { const char *path; int expect; const char *name; } cases[] = {
        { p_ok,    1, "complete-obj" },
        { p_trunc, 0, "truncated-no-newline" },
        { p_tiny,  0, "below-size-floor" },
        { "gp_selftest_missing.obj", 0, "missing-file" },
    };
    int fails = 0;
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        int got = obj_looks_complete(cases[i].path);
        int ok  = (got == cases[i].expect);
        fprintf(stderr, "[selftest] %-22s got=%d expect=%d -> %s\n",
                cases[i].name, got, cases[i].expect, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    /* Regression: weld geometry must not silently arm per-cube BPA rejection.
     * The old coupling shattered a known PHerc0139 cube from 5 sheets into 24
     * fragments before CVT/RVD. Only an explicit positive --grow-wind-tol may
     * opt into that experimental behavior. */
    {
        char *av[] = {
            "grid_pipeline", "grid", "out",
            "--umb-y", "3405", "--umb-x", "2878",
            "--wrap-pitch", "9.5"
        };
        GpOptions o;
        int rc = parse_args((int)(sizeof(av) / sizeof(av[0])), av, &o);
        int ok = rc == 0 && o.have_umb_y && o.have_umb_x
                 && o.have_wrap_pitch && !o.have_grow_wind_tol
                 && o.grow_wind_tol == 0.0;
        fprintf(stderr,
                "[selftest] %-22s grow_explicit=%d tol=%.2f -> %s\n",
                "weld-geometry-gate-off", o.have_grow_wind_tol,
                o.grow_wind_tol, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }
    {
        char *av[] = {
            "grid_pipeline", "grid", "out",
            "--umb-y", "3405", "--umb-x", "2878",
            "--wrap-pitch", "9.5", "--grow-wind-tol", ".45"
        };
        GpOptions o;
        int rc = parse_args((int)(sizeof(av) / sizeof(av[0])), av, &o);
        int ok = rc == 0 && o.have_grow_wind_tol
                 && o.grow_wind_tol > 0.449 && o.grow_wind_tol < 0.451;
        fprintf(stderr,
                "[selftest] %-22s grow_explicit=%d tol=%.2f -> %s\n",
                "explicit-growth-opt-in", o.have_grow_wind_tol,
                o.grow_wind_tol, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    remove(p_ok); remove(p_trunc); remove(p_tiny);
    fprintf(stderr, "=== grid_pipeline selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}

int main(int argc, char *argv[])
{
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return run_selftest();

    GpOptions opts;
    if (parse_args(argc, argv, &opts) != 0) return 1;

    fprintf(stderr,
        "grid_pipeline: grid=%s output=%s halo=%d tpc=%d max_conc=%d reject_garbage=%d\n",
        opts.grid_dir, opts.output_dir, opts.halo,
        opts.threads_per_cube, opts.max_concurrent, opts.reject_garbage);

    /* Mirror the disable to cube_mesh children (which inherit our environment
     * with a NULL env block) so the in-process defensive bail in MeshExtract
     * agrees with the orchestrator's --no-reject-garbage. Set once, before any
     * subprocess is spawned, so there is no thread race. */
    if (!opts.reject_garbage) {
#ifdef _WIN32
        _putenv_s("VESUVIUS_NO_REJECT_GARBAGE", "1");
#else
        setenv("VESUVIUS_NO_REJECT_GARBAGE", "1", 1);
#endif
    }

    /* Per-cube BPA winding rejection is deliberately SEPARATE from weld
     * geometry. A 2026-07-30 4x21x21 run accidentally coupled the two: on the
     * byte-identical z04480_y03328_x02816 prediction, Resplit changed from
     * 4->5 coherent sheets to 4->24 fragments before CVT/RVD. Clear inherited
     * BPA_GROW_* variables as well, so a caller's stale environment cannot
     * silently reproduce that corruption. The experimental gate now requires
     * an explicit positive --grow-wind-tol. */
    clear_process_env("BPA_GROW_UMBILICUS_Y");
    clear_process_env("BPA_GROW_UMBILICUS_X");
    clear_process_env("BPA_GROW_WRAP_PITCH");
    clear_process_env("BPA_GROW_WIND_TOL");
    if (opts.have_grow_wind_tol && opts.grow_wind_tol > 0.0) {
        set_process_env_double("BPA_GROW_UMBILICUS_Y", opts.umb_y);
        set_process_env_double("BPA_GROW_UMBILICUS_X", opts.umb_x);
        set_process_env_double("BPA_GROW_WRAP_PITCH", opts.wrap_pitch);
        set_process_env_double("BPA_GROW_WIND_TOL", opts.grow_wind_tol);
        fprintf(stderr,
            "WARNING: experimental BPA growth gate explicitly armed: "
            "umbilicus=(%.1f,%.1f) pitch=%.3f tol=%.3f turns\n",
            opts.umb_y, opts.umb_x, opts.wrap_pitch, opts.grow_wind_tol);
    } else {
        fprintf(stderr,
            "BPA growth gate OFF (coherent-sheet default; weld geometry is scoped to grid_weld)\n");
    }

    if (opts.simplify_engine == 1) {
        const char *ratio = getenv("VES_CVT_RATIO");
        fprintf(stderr, "CVT site ratio: %s\n", ratio && ratio[0] ? ratio : "default");
    }

    CubeJob *jobs = NULL;
    size_t n_jobs = 0;
    if (scan_cubes(opts.grid_dir, &jobs, &n_jobs) != 0) {
        fprintf(stderr, "ERROR: failed to scan %s/cubes_PRED\n",
                opts.grid_dir);
        return 1;
    }
    if (opts.max_cubes > 0 && (size_t)opts.max_cubes < n_jobs) {
        n_jobs = (size_t)opts.max_cubes;
    }
    fprintf(stderr, "Processing %zu cubes\n", n_jobs);

    if (n_jobs == 0) { free(jobs); return 1; }

    /* Prepare output directory structure. */
    char dump_dir[1024], log_dir[1024], cubes_dir[1024];
    snprintf(dump_dir, sizeof(dump_dir), "%s/dump", opts.output_dir);
    snprintf(log_dir, sizeof(log_dir), "%s/logs", opts.output_dir);
    snprintf(cubes_dir, sizeof(cubes_dir), "%s/cubes", opts.output_dir);
    ensure_dir(opts.output_dir);
    ensure_dir(dump_dir);
    ensure_dir(log_dir);
    ensure_dir(cubes_dir);
    if (write_scrollslice_source(&opts) != 0) {
        fprintf(stderr, "ERROR: could not write %s/scrollslice.source.json\n",
                opts.output_dir);
        free(jobs);
        return 1;
    }

    for (size_t i = 0; i < n_jobs; i++) {
        snprintf(jobs[i].log_path, sizeof(jobs[i].log_path),
                 "%s/%s.log", log_dir, jobs[i].cube_id);
    }

    /* --dry-run: report skip/run decisions without spawning anything. */
    if (opts.dry_run) {
        size_t would_skip = 0;
        for (size_t k = 0; k < n_jobs; k++) {
            char fp[1024];
            step12_obj_path(fp, sizeof(fp), dump_dir, jobs[k].cube_id);
            if (opts.skip_existing && obj_looks_complete(fp)) {
                would_skip++;
            } else {
                fprintf(stderr, "  [dry] RUN  %s\n", jobs[k].cube_id);
            }
        }
        fprintf(stderr,
                "[dry-run] %zu cubes: %zu SKIP (already done), %zu RUN\n",
                n_jobs, would_skip, n_jobs - would_skip);
        free(jobs);
        return 0;
    }

    /* Open summary CSV. */
    char summary_path[1024];
    snprintf(summary_path, sizeof(summary_path),
             "%s/pipeline_summary.csv", opts.output_dir);
    FILE *summary = fopen(summary_path, "w");
    if (summary) {
        fprintf(summary, "cube_id,exit_code,wall_seconds\n");
    }

    /* Garbage gate: a list of rejected cube_ids to record + a log of them.
     * With --reject-list we trust a precomputed list; otherwise the detector
     * runs inline per cube (see cube_is_garbage). */
    IdSet rejset = {0};
    int have_rejlist = 0;
    if (opts.reject_garbage && opts.reject_list) {
        if (idset_load(opts.reject_list, &rejset) == 0) {
            have_rejlist = 1;
            fprintf(stderr, "reject-list: %zu cube_ids from %s\n",
                    rejset.n, opts.reject_list);
        } else {
            fprintf(stderr, "warning: cannot read --reject-list %s; "
                    "detecting garbage inline instead\n", opts.reject_list);
        }
    }
    char rej_path[1024];
    snprintf(rej_path, sizeof(rej_path), "%s/rejected_cubes.txt", opts.output_dir);
    FILE *rejlog = opts.reject_garbage ? fopen(rej_path, "w") : NULL;

    /* Run cubes in parallel via OpenMP. MSVC OpenMP 2.0 requires the
     * loop counter declared outside the for-statement. */
    double t_start = ves_clock_sec();
    omp_set_dynamic(0);
    omp_set_num_threads(opts.max_concurrent);

#ifdef _WIN32
    /* Enforce --threads-per-cube on Windows children. Set ONCE here,
     * before the spawn threads exist (so no SetEnvironmentVariable race),
     * and inherited by every CreateProcessA(lpEnvironment=NULL) child.
     * All cubes share the same thread budget, so a global parent-side set
     * is exactly equivalent to a per-child block — and sidesteps the
     * ANSI lpEnvironment conversion quirks (ERROR_INVALID_PARAMETER).
     * The POSIX branch does the equivalent putenv() in the forked child.
     * Without this, children fell back to ves_cpu_count() and every
     * concurrent cube opened a full-core OpenMP team inside MLS
     * (max_concurrent x n_cores threads of thrash). */
    {
        char thr_str[32];
        snprintf(thr_str, sizeof(thr_str), "%d", opts.threads_per_cube);
        SetEnvironmentVariableA("VESUVIUS_THREADS", thr_str);
    }
#endif

    int n_ok = 0, n_fail = 0, n_skip = 0, n_reject = 0, n_reject_empty = 0;
    int i;
    int n_jobs_i = (int)n_jobs;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:n_ok,n_fail,n_skip,n_reject,n_reject_empty)
    for (i = 0; i < n_jobs_i; i++) {
        /* Resume: a cube with an already-complete final OBJ is left untouched. */
        if (opts.skip_existing) {
            char done_path[1024];
            step12_obj_path(done_path, sizeof(done_path), dump_dir, jobs[i].cube_id);
            if (obj_looks_complete(done_path)) {
                jobs[i].exit_code = 0;
                jobs[i].wall_seconds = 0.0;
                n_ok++;
                n_skip++;
                #pragma omp critical (gp_log)
                {
                    if (summary) {
                        fprintf(summary, "%s,0,0.00\n", jobs[i].cube_id);
                        fflush(summary);
                    }
                }
                continue;
            }
        }
        /* Garbage gate: skip garbage prediction cubes entirely (no subprocess).
         * Two kinds: solid-slab (thick filled box) and empty (too little FG to
         * form one meshable component -- would otherwise mesh to a guaranteed
         * empty FAIL). Both recorded with the GP_REJECT_CODE sentinel and listed
         * in rejected_cubes.txt. A supplied --reject-list is treated as slab. */
        if (opts.reject_garbage) {
            int kind = have_rejlist
                         ? (idset_has(&rejset, jobs[i].cube_id) ? 1 : 0)
                         : cube_reject_kind(jobs[i].tiff_path);
            if (kind) {
                jobs[i].exit_code = GP_REJECT_CODE;
                jobs[i].wall_seconds = 0.0;
                n_reject++;
                if (kind == 2) n_reject_empty++;
                #pragma omp critical (gp_log)
                {
                    fprintf(stderr, "  [reject] %s (%s)\n", jobs[i].cube_id,
                            kind == 2 ? "empty prediction -- no meshable component"
                                      : "garbage solid-slab prediction");
                    if (summary) {
                        fprintf(summary, "%s,%d,0.00\n",
                                jobs[i].cube_id, GP_REJECT_CODE);
                        fflush(summary);
                    }
                    if (rejlog) {
                        fprintf(rejlog, "%s\n", jobs[i].cube_id);
                        fflush(rejlog);
                    }
                }
                continue;
            }
        }
        double cube_start = ves_clock_sec();
        int rc = run_one_cube(opts.exe_path, &jobs[i],
                              opts.output_dir, dump_dir,
                              opts.halo, opts.threads_per_cube,
                              opts.skip_qem, opts.trim_inset,
                              !opts.full_dumps, opts.simplify_engine);
        jobs[i].exit_code = rc;
        jobs[i].wall_seconds = ves_clock_sec() - cube_start;
        if (rc == 0) n_ok++; else n_fail++;
        #pragma omp critical (gp_log)
        {
            fprintf(stderr, "  [%4d/%zu] %s exit=%d t=%.1fs\n",
                    i + 1, n_jobs, jobs[i].cube_id, rc,
                    jobs[i].wall_seconds);
            if (summary) {
                fprintf(summary, "%s,%d,%.2f\n",
                        jobs[i].cube_id, rc, jobs[i].wall_seconds);
                fflush(summary);
            }
        }
    }
    if (summary) fclose(summary);
    if (rejlog) fclose(rejlog);
    idset_free(&rejset);
    double t_total = ves_clock_sec() - t_start;

    fprintf(stderr,
            "Cubes: %d ok (%d skipped already-done), %d rejected "
            "(%d solid-slab, %d empty), %d fail in %.1fs\n",
            n_ok, n_skip, n_reject, n_reject - n_reject_empty, n_reject_empty,
            n_fail, t_total);
    if (n_reject > 0)
        fprintf(stderr, "  rejected cube list -> %s\n", rej_path);

    if (opts.check_determinism && n_ok > 0) {
        char dump_dir_b[1024];
        snprintf(dump_dir_b, sizeof(dump_dir_b),
                 "%s/dump_check", opts.output_dir);
        ensure_dir(dump_dir_b);
        fprintf(stderr, "Determinism check: re-running %zu cubes\n", n_jobs);
        int n_diff = 0;
        int j;
        int n_jobs_j = (int)n_jobs;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:n_diff)
        for (j = 0; j < n_jobs_j; j++) {
            if (jobs[j].exit_code != 0) continue;
            char log_b[1024];
            snprintf(log_b, sizeof(log_b), "%s.check", jobs[j].log_path);
            CubeJob copy = jobs[j];
            strncpy(copy.log_path, log_b, sizeof(copy.log_path) - 1);
            int rc_b = run_one_cube(opts.exe_path, &copy,
                                     opts.output_dir, dump_dir_b,
                                     opts.halo, opts.threads_per_cube,
                                     opts.skip_qem, opts.trim_inset,
                                     !opts.full_dumps, opts.simplify_engine);
            if (rc_b != 0) { n_diff++; continue; }

            /* Compare the final per-cube OBJ. */
            char path_a[1024], path_b[1024];
            snprintf(path_a, sizeof(path_a),
                "%s/%s/%s_step12_final/%s_step12_final_all.obj",
                dump_dir, jobs[j].cube_id, jobs[j].cube_id, jobs[j].cube_id);
            snprintf(path_b, sizeof(path_b),
                "%s/%s/%s_step12_final/%s_step12_final_all.obj",
                dump_dir_b, jobs[j].cube_id, jobs[j].cube_id, jobs[j].cube_id);
            FILE *fa = fopen(path_a, "rb");
            FILE *fb = fopen(path_b, "rb");
            if (!fa || !fb) {
                if (fa) fclose(fa);
                if (fb) fclose(fb);
                n_diff++;
                continue;
            }
            int diff = 0;
            char ba[8192], bb[8192];
            for (;;) {
                size_t na = fread(ba, 1, sizeof(ba), fa);
                size_t nb = fread(bb, 1, sizeof(bb), fb);
                if (na != nb) { diff = 1; break; }
                if (na == 0) break;
                if (memcmp(ba, bb, na) != 0) { diff = 1; break; }
            }
            fclose(fa);
            fclose(fb);
            if (diff) {
                #pragma omp critical (gp_log)
                fprintf(stderr, "  [diff] %s OBJ differs across runs\n",
                        jobs[j].cube_id);
                n_diff++;
            }
        }
        fprintf(stderr, "Determinism check: %d/%d cubes deterministic\n",
                (int)n_jobs - n_diff, (int)n_jobs);
    }

    int weld_crashed = 0;
    int weld_audit_warn = 0;
    if (!opts.skip_weld && n_ok > 0) {
        char weld_out[1024];
        snprintf(weld_out, sizeof(weld_out), "%s/welded.obj", opts.output_dir);
        /* Arm grid_weld's seam winding + phase gates (it reads them from env:
         * SEAM_UMBILICUS_Y/X + SEAM_WRAP_PITCH). Deliberately set HERE — after
         * the cube fleet has finished, immediately before the weld spawn — so
         * per-cube meshing NEVER sees them: pinhole_fill.c also reads
         * SEAM_WRAP_PITCH, and arming it during meshing would change per-cube
         * output. Regression note (2026-07-18): the Jul-12 rebuild harness
         * stopped passing these, the weld ran gates-off, and the 4x5x5 shipped
         * 661 seam-band boundary loops instead of 16. */
        if (opts.have_umb_y && opts.have_umb_x && opts.have_wrap_pitch) {
            char envv[64];
            snprintf(envv, sizeof(envv), "%.6g", opts.umb_y);
#ifdef _WIN32
            SetEnvironmentVariableA("SEAM_UMBILICUS_Y", envv);
#else
            setenv("SEAM_UMBILICUS_Y", envv, 1);
#endif
            snprintf(envv, sizeof(envv), "%.6g", opts.umb_x);
#ifdef _WIN32
            SetEnvironmentVariableA("SEAM_UMBILICUS_X", envv);
#else
            setenv("SEAM_UMBILICUS_X", envv, 1);
#endif
            snprintf(envv, sizeof(envv), "%.6g", opts.wrap_pitch);
#ifdef _WIN32
            SetEnvironmentVariableA("SEAM_WRAP_PITCH", envv);
#else
            setenv("SEAM_WRAP_PITCH", envv, 1);
#endif
            fprintf(stderr,
                "Seam gates armed for weld: umbilicus=(%.1f,%.1f) pitch=%.2f\n",
                opts.umb_y, opts.umb_x, opts.wrap_pitch);
        } else {
            fprintf(stderr,
                "Seam gates NOT armed (pass --umb-y/--umb-x/--wrap-pitch); "
                "grid_weld's winding/phase gates run OFF\n");
        }

        fprintf(stderr, "Running grid_weld -> %s\n", weld_out);

#ifdef _WIN32
        /* Run grid_weld with stderr/stdout to a log file and capture the
         * actual exit code. ves_run_subprocess collapses all non-zero
         * exits to -1 which hides "ran but audit had warnings" — we want
         * to distinguish that from a real crash. */
        char weld_log[1024];
        snprintf(weld_log, sizeof(weld_log),
                 "%s/grid_weld.log", opts.output_dir);
        char weld_cmd[4096];
        /* The BPA seam-weld is unconditional in grid_weld now (it IS the weld;
         * the bitwise hash-join it replaced is gone). step0 trims each cube to
         * its face, so the cubes concatenate with a ~1-vox seam gap that the
         * BPA bridge closes. No flag needed. */
        snprintf(weld_cmd, sizeof(weld_cmd), "\"%s\" \"%s\" \"%s\"",
                 opts.weld_path, dump_dir, weld_out);
        SECURITY_ATTRIBUTES sa;
        sa.nLength = sizeof(sa);
        sa.lpSecurityDescriptor = NULL;
        sa.bInheritHandle = TRUE;
        HANDLE hLog = CreateFileA(weld_log, GENERIC_WRITE,
                                   FILE_SHARE_READ | FILE_SHARE_WRITE,
                                   &sa, CREATE_ALWAYS,
                                   FILE_ATTRIBUTE_NORMAL, NULL);
        STARTUPINFOA si;
        PROCESS_INFORMATION pi;
        memset(&si, 0, sizeof(si));
        si.cb = sizeof(si);
        if (hLog != INVALID_HANDLE_VALUE) {
            si.dwFlags = STARTF_USESTDHANDLES;
            si.hStdOutput = hLog;
            si.hStdError  = hLog;
            si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
        }
        memset(&pi, 0, sizeof(pi));
        BOOL ok = CreateProcessA(opts.weld_path, weld_cmd, NULL, NULL,
                                  hLog != INVALID_HANDLE_VALUE,
                                  0, NULL, NULL, &si, &pi);
        if (hLog != INVALID_HANDLE_VALUE) CloseHandle(hLog);
        if (!ok) {
            fprintf(stderr,
                "grid_weld FAILED to start (CreateProcess err=%lu)\n",
                (unsigned long)GetLastError());
            weld_crashed = 1;
        } else {
            WaitForSingleObject(pi.hProcess, INFINITE);
            DWORD exit_code = 1;
            GetExitCodeProcess(pi.hProcess, &exit_code);
            CloseHandle(pi.hProcess);
            CloseHandle(pi.hThread);
            /* Windows NTSTATUS error codes (e.g. STATUS_ACCESS_VIOLATION
             * 0xC0000005) have the top bit set. Anything in [0, 0xFF] is
             * a normal exit (success or warning). */
            if (exit_code > 0xFF) {
                fprintf(stderr,
                    "grid_weld CRASHED (exit=0x%08lX)\n",
                    (unsigned long)exit_code);
                weld_crashed = 1;
            } else if (exit_code != 0) {
                fprintf(stderr,
                    "grid_weld exit=%lu (non-zero; check %s for warnings)\n",
                    (unsigned long)exit_code, weld_log);
                weld_audit_warn = 1;
            } else {
                fprintf(stderr, "grid_weld exit=0\n");
            }
        }
#else
        const char *weld_argv[4];
        weld_argv[0] = opts.weld_path;
        weld_argv[1] = dump_dir;
        weld_argv[2] = weld_out;
        weld_argv[3] = NULL;   /* BPA seam-weld is unconditional now */
        int weld_rc = ves_run_subprocess(opts.weld_path, weld_argv, 0.0);
        /* ves_run_subprocess loses the distinction; treat -1 as crash for
         * the simple POSIX path. The Windows branch above is the
         * authoritative implementation. */
        if (weld_rc != 0) weld_audit_warn = 1;
        fprintf(stderr, "grid_weld rc=%d\n", weld_rc);
#endif

        /* Surface manifold-audit findings from weld_report.json if the
         * weld actually wrote one (succeeded enough to emit a report). */
        char report_path[1280];
        snprintf(report_path, sizeof(report_path),
                 "%s.weld_report.json", weld_out);
        FILE *fr = fopen(report_path, "rb");
        if (fr) {
            char buf[4096];
            size_t n = fread(buf, 1, sizeof(buf) - 1, fr);
            buf[n] = '\0';
            fclose(fr);
            fprintf(stderr, "weld_report.json:\n%s\n", buf);
        }
    }

    free(jobs);
    /* Cubes that failed are real errors; weld_audit_warn is just a notice;
     * weld_crashed is a real error. */
    if (n_fail > 0)     fprintf(stderr, "Cubes failed: %d\n", n_fail);
    if (weld_audit_warn) fprintf(stderr, "grid_weld reported audit warnings (non-fatal)\n");
    if (weld_crashed)   fprintf(stderr, "grid_weld crashed -- pipeline aborted\n");
    return (n_fail == 0 && !weld_crashed) ? 0 : 1;
}
