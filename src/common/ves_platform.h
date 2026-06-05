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

/* ================================================================
 * Compiler compat
 * ================================================================ */

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
  #include <io.h>        /* _unlink, _access */
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
  #include <string.h>    /* memcpy (ves_temp_dir); on Windows it comes via <windows.h> */
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

#ifdef _WIN32
  #define VES_FADVISE_DONTNEED  0
  static inline int ves_fadvise(int fd, long offset, long len, int advice)
  {
      (void)fd; (void)offset; (void)len; (void)advice;
      return 0; /* no-op on Windows */
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

#endif /* VES_PLATFORM_H */
