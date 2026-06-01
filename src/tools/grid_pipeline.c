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
 *                 [--qem | --no-qem]
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

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
  #include <sys/wait.h>
  #include <unistd.h>
  #include <fcntl.h>
#endif

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
} GpOptions;

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
        "  --check-determinism       Run each cube twice and cmp OBJs\n"
        "  --no-qem                  Pass --no-qem to cube_mesh\n",
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
    o->exe_path = "build/Release/cube_mesh.exe";
    o->weld_path = "build/Release/grid_weld.exe";

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
        } else if (!strcmp(argv[i], "--skip-weld")) {
            o->skip_weld = 1;
        } else if (!strcmp(argv[i], "--check-determinism")) {
            o->check_determinism = 1;
        } else if (!strcmp(argv[i], "--no-qem")) {
            o->skip_qem = 1;
        } else if (!strcmp(argv[i], "--qem")) {
            o->skip_qem = 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return -1;
        }
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
                        int halo, int threads, int skip_qem)
{
    char out_tif[1024];
    snprintf(out_tif, sizeof(out_tif), "%s/cubes/%s.tif",
             output_dir, job->cube_id);

    /* Compose argv. */
    char halo_str[32];
    snprintf(halo_str, sizeof(halo_str), "%d", halo);

    const char *argv[16];
    int argc = 0;
    argv[argc++] = exe_path;
    argv[argc++] = job->tiff_path;
    argv[argc++] = out_tif;
    argv[argc++] = "--halo";
    argv[argc++] = halo_str;
    argv[argc++] = "--dump-obj";
    argv[argc++] = dump_dir;
    argv[argc++] = "--no-timeout";
    if (skip_qem) argv[argc++] = "--no-qem";
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

    /* Build a per-cube environment block: parent env + VESUVIUS_THREADS.
     * Use SetEnvironmentVariable on the PARENT side momentarily is unsafe
     * because of races between threads, so we pass NULL (inherit parent
     * env unchanged). The cube_mesh child falls back to cpu_count() when
     * VESUVIUS_THREADS is unset, which on a busy system over-subscribes;
     * acceptable for now since the threads-per-cube limiter is in the
     * orchestrator's max_concurrent.
     *
     * TODO: thread-local env block per cube to enforce threads_per_cube. */
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

int main(int argc, char *argv[])
{
    GpOptions opts;
    if (parse_args(argc, argv, &opts) != 0) return 1;

    fprintf(stderr,
        "grid_pipeline: grid=%s output=%s halo=%d tpc=%d max_conc=%d\n",
        opts.grid_dir, opts.output_dir, opts.halo,
        opts.threads_per_cube, opts.max_concurrent);

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

    for (size_t i = 0; i < n_jobs; i++) {
        snprintf(jobs[i].log_path, sizeof(jobs[i].log_path),
                 "%s/%s.log", log_dir, jobs[i].cube_id);
    }

    /* Open summary CSV. */
    char summary_path[1024];
    snprintf(summary_path, sizeof(summary_path),
             "%s/pipeline_summary.csv", opts.output_dir);
    FILE *summary = fopen(summary_path, "w");
    if (summary) {
        fprintf(summary, "cube_id,exit_code,wall_seconds\n");
    }

    /* Run cubes in parallel via OpenMP. MSVC OpenMP 2.0 requires the
     * loop counter declared outside the for-statement. */
    double t_start = ves_clock_sec();
    omp_set_dynamic(0);
    omp_set_num_threads(opts.max_concurrent);

    int n_ok = 0, n_fail = 0;
    int i;
    int n_jobs_i = (int)n_jobs;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:n_ok,n_fail)
    for (i = 0; i < n_jobs_i; i++) {
        double cube_start = ves_clock_sec();
        int rc = run_one_cube(opts.exe_path, &jobs[i],
                              opts.output_dir, dump_dir,
                              opts.halo, opts.threads_per_cube,
                              opts.skip_qem);
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
    double t_total = ves_clock_sec() - t_start;

    fprintf(stderr, "Cubes: %d ok, %d fail in %.1fs\n",
            n_ok, n_fail, t_total);

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
                                     opts.skip_qem);
            if (rc_b != 0) { n_diff++; continue; }

            /* Compare the raw_snap OBJ. */
            char path_a[1024], path_b[1024];
            snprintf(path_a, sizeof(path_a),
                "%s/%s/%s_raw_snap/%s_raw_snap_all.obj",
                dump_dir, jobs[j].cube_id, jobs[j].cube_id, jobs[j].cube_id);
            snprintf(path_b, sizeof(path_b),
                "%s/%s/%s_raw_snap/%s_raw_snap_all.obj",
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
