/*
 * ves_platform.h — Cross-platform abstraction layer.
 *
 * Provides portable wrappers for POSIX APIs used throughout the pipeline.
 * On Linux/GCC: thin wrappers that map back to the original POSIX calls.
 * On Windows/MSVC: equivalent Win32 API calls.
 *
 * Include this header instead of directly including:
 *   <unistd.h>, <pthread.h>, <sys/stat.h>, <sys/types.h>, <sys/wait.h>,
 *   <signal.h> (for SIGALRM), <time.h> (for clock_gettime).
 *
 * Two functions are too large for inline and live in ves_platform.c:
 *   ves_run_subprocess()       — portable fork/exec or CreateProcess
 *   ves_hard_timeout_start()   — portable SIGALRM or timer thread
 *   ves_hard_timeout_cancel()
 */
#ifndef VES_PLATFORM_H
#define VES_PLATFORM_H

#include <stddef.h>
#include <stdio.h>
#include <string.h>    /* memcpy (ves_temp_dir), strcmp/strstr (ves_selftest_clean) */

/* ================================================================
 * Compiler compat
 * ================================================================ */

/* Thread-local storage for per-call solver state that concurrent callers
 * must not share (the QP deadlines of convex_qp.c and active_set_qp.c).
 * Same spelling rules as EXCEPT_THREAD_LOCAL in except.h. */
#if defined(_MSC_VER)
  #define VES_THREAD_LOCAL __declspec(thread)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_THREADS__)
  #define VES_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
  #define VES_THREAD_LOCAL __thread
#else
  #define VES_THREAD_LOCAL
#endif

#ifdef _MSC_VER
  /* MSVC does not support C99 'restrict'; use __restrict instead */
  #ifndef restrict
    #define restrict __restrict
  #endif

  /* Suppress MSVC warnings for standard C functions (fopen, snprintf, etc.) */
  #ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS
  #endif

  /* GCC __attribute__ — no-op on MSVC */
  #ifndef __attribute__
    #define __attribute__(x)
  #endif
#else
  /* GCC/Clang on POSIX — define _POSIX_C_SOURCE globally so downstream
   * headers expose clock_gettime, nanosleep, etc. */
  #ifndef _POSIX_C_SOURCE
    #define _POSIX_C_SOURCE 200809L
  #endif
#endif

/* ================================================================
 * Platform headers
 * ================================================================ */

#ifdef _WIN32
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  #include <direct.h>    /* _mkdir */
  #include <io.h>        /* _unlink, _access, _setmode */
  #include <fcntl.h>     /* _O_BINARY */
  #include <process.h>   /* _getpid */
  #include <signal.h>    /* sig_atomic_t */
  #include <time.h>      /* time_t, etc. */
#else
  #include <unistd.h>
  #include <sys/stat.h>
  #include <sys/types.h>
  #include <time.h>
  #include <signal.h>
  #include <pthread.h>
  #include <dirent.h>    /* opendir (ves_selftest_clean) */
#endif

/* ================================================================
 * High-resolution clock: ves_clock_sec()
 *
 * Returns monotonic wall-clock time in seconds as double.
 * ================================================================ */

static inline double ves_clock_sec(void)
{
#ifdef _WIN32
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (double)count.QuadPart / (double)freq.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

/* ================================================================
 * Mutex: ves_mutex_t
 *
 * Lightweight mutex for serializing libtiff calls.
 * Windows: SRWLOCK (no cleanup needed, no allocation).
 * Linux:   pthread_mutex_t.
 * ================================================================ */

#ifdef _WIN32

typedef SRWLOCK ves_mutex_t;

#define VES_MUTEX_INITIALIZER  SRWLOCK_INIT

static inline void ves_mutex_init(ves_mutex_t *m)
{
    InitializeSRWLock(m);
}

static inline void ves_mutex_lock(ves_mutex_t *m)
{
    AcquireSRWLockExclusive(m);
}

static inline void ves_mutex_unlock(ves_mutex_t *m)
{
    ReleaseSRWLockExclusive(m);
}

#else /* POSIX */

typedef pthread_mutex_t ves_mutex_t;

#define VES_MUTEX_INITIALIZER  PTHREAD_MUTEX_INITIALIZER

static inline void ves_mutex_init(ves_mutex_t *m)
{
    pthread_mutex_init(m, NULL);
}

static inline void ves_mutex_lock(ves_mutex_t *m)
{
    pthread_mutex_lock(m);
}

static inline void ves_mutex_unlock(ves_mutex_t *m)
{
    pthread_mutex_unlock(m);
}

#endif

/* ================================================================
 * CPU count: ves_cpu_count()
 * ================================================================ */

static inline int ves_cpu_count(void)
{
#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (int)si.dwNumberOfProcessors;
#else
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 2;
#endif
}

/* ================================================================
 * PID: ves_getpid()
 * ================================================================ */

static inline int ves_getpid(void)
{
#ifdef _WIN32
    return (int)_getpid();
#else
    return (int)getpid();
#endif
}

/* ================================================================
 * File operations: ves_unlink(), ves_mkdir()
 * ================================================================ */

static inline int ves_unlink(const char *path)
{
#ifdef _WIN32
    return _unlink(path);
#else
    return unlink(path);
#endif
}

static inline int ves_mkdir(const char *path)
{
#ifdef _WIN32
    return _mkdir(path);
#else
    return mkdir(path, 0755);
#endif
}

static inline int ves_rmdir(const char *path)
{
#ifdef _WIN32
    return _rmdir(path);
#else
    return rmdir(path);
#endif
}

/* ================================================================
 * Selftest scratch: VES_SELFTEST_DIR, ves_selftest_dir(), ves_selftest_clean()
 *
 * Every selftest writes its scratch under VES_SELFTEST_DIR (relative to the
 * working directory, the repo root), never under output/.  User, 2026-10-12:
 * "put all the selftests into a new selftest/ folder from now on that can be
 * blown out" and "they should clean up after a pass".  The folder may be
 * deleted at any time.
 *
 * ves_selftest_dir(buf, len, name): "selftest/<name>_<pid>" into buf, created
 *   together with selftest/.  The pid keeps two runs of the same selftest out
 *   of each other's files.  0 on success.
 * ves_selftest_clean(path): after a PASS, delete a file or a directory tree
 *   under VES_SELFTEST_DIR.  Junctions and symlinks inside are removed as
 *   links, never followed.  A path outside VES_SELFTEST_DIR, the root itself
 *   or a path containing ".." is refused (-1).  A failing selftest keeps its
 *   scratch for inspection.  0 when the path is gone.
 * ================================================================ */

#define VES_SELFTEST_DIR "selftest"

static inline int ves_selftest_dir(char *buf, size_t len, const char *name)
{
    int k = snprintf(buf, len, VES_SELFTEST_DIR "/%s_%d", name, ves_getpid());
    if (k < 0 || (size_t)k >= len) return -1;
    (void)ves_mkdir(VES_SELFTEST_DIR);
    (void)ves_mkdir(buf);
    return 0;
}

/* Recursive delete behind ves_selftest_clean(); never follows a link. */
static inline int ves_remove_tree_(const char *path)
{
    char sub[1024];
    int rc = 0;
#ifdef _WIN32
    DWORD attr = GetFileAttributesA(path);
    if (attr == INVALID_FILE_ATTRIBUTES) return 0;
    if (attr & FILE_ATTRIBUTE_READONLY)
        (void)SetFileAttributesA(path, attr & ~(DWORD)FILE_ATTRIBUTE_READONLY);
    if (!(attr & FILE_ATTRIBUTE_DIRECTORY)) return DeleteFileA(path) ? 0 : -1;
    if (attr & FILE_ATTRIBUTE_REPARSE_POINT)        /* junction / dir symlink */
        return RemoveDirectoryA(path) ? 0 : -1;
    {
        WIN32_FIND_DATAA fd;
        HANDLE h;
        int k = snprintf(sub, sizeof sub, "%s/*", path);
        if (k < 0 || (size_t)k >= sizeof sub) return -1;
        h = FindFirstFileA(sub, &fd);
        if (h != INVALID_HANDLE_VALUE) {
            do {
                if (strcmp(fd.cFileName, ".") == 0 || strcmp(fd.cFileName, "..") == 0)
                    continue;
                k = snprintf(sub, sizeof sub, "%s/%s", path, fd.cFileName);
                if (k < 0 || (size_t)k >= sizeof sub || ves_remove_tree_(sub) != 0)
                    rc = -1;
            } while (FindNextFileA(h, &fd));
            FindClose(h);
        }
    }
    if (!RemoveDirectoryA(path)) rc = -1;
#else
    struct stat st;
    DIR *d = NULL;
    if (lstat(path, &st) != 0) return 0;
    if (!S_ISDIR(st.st_mode)) return unlink(path) == 0 ? 0 : -1;  /* files, links */
    d = opendir(path);
    if (d != NULL) {
        struct dirent *e = NULL;
        while ((e = readdir(d)) != NULL) {
            int k = 0;
            if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0)
                continue;
            k = snprintf(sub, sizeof sub, "%s/%s", path, e->d_name);
            if (k < 0 || (size_t)k >= sizeof sub || ves_remove_tree_(sub) != 0)
                rc = -1;
        }
        closedir(d);
    }
    if (rmdir(path) != 0) rc = -1;
#endif
    return rc;
}

static inline int ves_selftest_clean(const char *path)
{
    size_t n = strlen(VES_SELFTEST_DIR);
    if (path == NULL || strncmp(path, VES_SELFTEST_DIR, n) != 0
        || (path[n] != '/' && path[n] != '\\') || path[n + 1] == '\0'
        || strstr(path, "..") != NULL)
        return -1;
    return ves_remove_tree_(path);
}

/* ================================================================
 * Sleep: ves_sleep_ms()
 * ================================================================ */

static inline void ves_sleep_ms(int ms)
{
#ifdef _WIN32
    Sleep((DWORD)ms);
#else
    struct timespec ts;
    ts.tv_sec  = ms / 1000;
    ts.tv_nsec = (long)(ms % 1000) * 1000000L;
    nanosleep(&ts, NULL);
#endif
}

/* ================================================================
 * Temp directory: ves_temp_dir()
 *
 * Writes the system temp directory path into buf (with trailing separator).
 * Returns 0 on success.
 * ================================================================ */

static inline int ves_temp_dir(char *buf, size_t len)
{
#ifdef _WIN32
    DWORD n = GetTempPathA((DWORD)len, buf);
    return (n > 0 && n < (DWORD)len) ? 0 : -1;
#else
    const char *tmp = "/tmp/";
    size_t slen = 5; /* strlen("/tmp/") */
    if (slen >= len) return -1;
    memcpy(buf, tmp, slen + 1);
    return 0;
#endif
}

/* ================================================================
 * Path utilities: ves_path_basename()
 *
 * Returns pointer to the filename portion of a path, handling both
 * '/' and '\' separators (Windows paths may use either).
 * ================================================================ */

static inline const char *ves_path_basename(const char *path)
{
    const char *last_sep = NULL;
    const char *p = path;
    while (*p) {
        if (*p == '/' || *p == '\\')
            last_sep = p;
        p++;
    }
    return last_sep ? last_sep + 1 : path;
}

/* ================================================================
 * posix_fadvise — no-op on Windows
 * ================================================================ */

#if defined(_WIN32) || defined(__APPLE__)
  #define VES_FADVISE_DONTNEED  0
  static inline int ves_fadvise(int fd, long offset, long len, int advice)
  {
      (void)fd; (void)offset; (void)len; (void)advice;
      return 0; /* no-op on Windows/macOS (no posix_fadvise) */
  }
#else
  #include <fcntl.h>
  #define VES_FADVISE_DONTNEED  POSIX_FADV_DONTNEED
  static inline int ves_fadvise(int fd, long offset, long len, int advice)
  {
      return posix_fadvise(fd, offset, len, advice);
  }
#endif

/* ================================================================
 * ves_stdin_set_binary() — put stdin into binary (untranslated) mode.
 *
 * On Windows, stdin defaults to text mode: CRLF pairs are collapsed and
 * a 0x1A byte reads as EOF, silently corrupting raw byte streams (e.g.
 * cube_mesh --stdin-raw). POSIX has no text/binary distinction.
 * ================================================================ */

#ifdef _WIN32
  static inline int ves_stdin_set_binary(void)
  {
      return (_setmode(_fileno(stdin), _O_BINARY) == -1) ? -1 : 0;
  }
#else
  static inline int ves_stdin_set_binary(void)
  {
      return 0; /* POSIX streams are always binary */
  }
#endif

/* ================================================================
 * Directory creation: ves_ensure_parent_dir()
 *
 * Given a file path, creates all parent directories that don't exist.
 * E.g. for "output/sub/foo.tif", creates "output/" and "output/sub/".
 * Handles both '/' and '\' separators.
 *
 * Implemented in ves_platform.c (needs string mutation).
 *
 * Returns 0 on success, -1 on failure.
 * ================================================================ */

int ves_ensure_parent_dir(const char *filepath);

/* ================================================================
 * Subprocess execution: ves_run_subprocess()
 *
 * Launches an external process with a timeout.
 * Implemented in ves_platform.c (too large for inline).
 *
 * exe:         path to executable
 * argv:        NULL-terminated argument array (argv[0] = program name)
 * timeout_sec: maximum wall time; <= 0 means no timeout
 *
 * Returns 0 on success (process exited with code 0), -1 on failure/timeout.
 * ================================================================ */

int ves_run_subprocess(const char *exe, const char *const *argv,
                       double timeout_sec);

/* Same, but the child's stdout+stderr append to log_path (NULL = suppress,
 * identical to ves_run_subprocess).  Orchestrators use this so every spawned
 * stage keeps a full log per the house logging rule. */
int ves_run_subprocess_logged(const char *exe, const char *const *argv,
                              double timeout_sec, const char *log_path);

/* ================================================================
 * Hard timeout: ves_hard_timeout_start / cancel
 *
 * Implemented in ves_platform.c.
 *
 * On Linux:  sigaction(SIGALRM) + alarm().
 * On Windows: background thread that sleeps then sets the flag.
 *
 * flag: pointer to a volatile sig_atomic_t that will be set to 1 on timeout.
 * ================================================================ */

void ves_hard_timeout_start(double seconds, volatile sig_atomic_t *flag);
void ves_hard_timeout_cancel(void);

/* ================================================================
 * Thread count
 *
 * Set the OpenMP team size for subsequent parallel regions.  n <= 0
 * means "leave the runtime default alone".  A no-op when the build
 * has OpenMP disabled, so callers need no #ifdef of their own.
 * ================================================================ */

void ves_omp_set_threads(int n);
/* Leave the process now with exit status `code`, after fflush(NULL): no
 * atexit handlers, no CRT/OpenMP/DLL teardown.  For a batch tool whose
 * outputs are all written and closed.  2026-10-05: one cube_mesh in ~2,600
 * hung at exit (output complete, 0 CPU, main + OpenMP worker + timer thread
 * alive) -- a teardown deadlock this exit path cannot enter.  Does not return. */
#if defined(_MSC_VER)
#define VES_NORETURN __declspec(noreturn)
#elif defined(__GNUC__)
#define VES_NORETURN __attribute__((noreturn))
#else
#define VES_NORETURN
#endif
VES_NORETURN void ves_fast_exit(int code);

/* Run fn(arg) on a worker thread bounded by the WORKER'S CPU time (cpu_ms) and,
 * as a backstop for a blocked or starved worker, by wall time (wall_ms).
 * Returns 0 when fn finished, 1 when its CPU budget ran out, 2 when the wall
 * backstop expired -- in both cases the worker was TerminateThread'ed, so use it
 * only for code that may be killed (no locks held while it spins) -- and -1
 * when no thread could be started (the caller should run fn inline).  A CPU
 * bound does not fire on a starved machine, so the outcome does not depend on
 * the load (2026-10-05: a 10 s WALL bound in hole_fill abandoned holes and
 * killed workers mid-flight under load).  *out_cpu_ms (optional) receives the
 * worker's CPU time.  Windows; elsewhere it returns -1. */
#ifdef _WIN32
typedef unsigned (__stdcall *VesThreadFn)(void *);
#else
typedef unsigned (*VesThreadFn)(void *);
#endif
int ves_run_thread_cpu_bounded(VesThreadFn fn, void *arg, unsigned cpu_ms,
                               unsigned wall_ms, unsigned *out_cpu_ms);

/* Selftest of ves_run_thread_cpu_bounded (quick, spinning, sleeping and
 * blocked workers).  Returns 0 on success. */
int ves_platform_selftest(void);

/* VES_PROF_OUT=<file>: start the in-process sampling profiler (~1 kHz leaf
 * samples of every thread, rewritten to <file> about once a second; resolve
 * with scripts/pipeline/ves_prof_report.py).  No-op when unset or on POSIX.
 * Diagnostics only: outputs are unaffected. */
void ves_prof_start_from_env(void);
/* Write the samples now (ves_fast_exit calls it before terminating). */
void ves_prof_flush_now(void);

#endif /* VES_PLATFORM_H */
