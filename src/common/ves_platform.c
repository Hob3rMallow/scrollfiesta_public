/*
 * ves_platform.c — Platform abstraction: subprocess + hard timeout.
 *
 * Contains the two functions too large for inline in ves_platform.h:
 *   1. ves_run_subprocess()      — fork/exec (Linux) or CreateProcess (Windows)
 *   2. ves_hard_timeout_start()  — SIGALRM (Linux) or timer thread (Windows)
 *      ves_hard_timeout_cancel()
 *   3. ves_fast_exit(), ves_run_thread_cpu_bounded() (2026-10-05)
 */
#include "ves_platform.h"

#include <string.h>
#include <stdlib.h>
#include <errno.h>
#ifdef _WIN32
#include <process.h>   /* _beginthreadex */
#endif

/* ================================================================
 * Directory creation: ves_ensure_parent_dir()
 * ================================================================ */

int ves_ensure_parent_dir(const char *filepath)
{
    if (!filepath) return -1;

    /* Find last separator to extract parent directory */
    const char *last_sep = NULL;
    for (const char *p = filepath; *p; p++) {
        if (*p == '/' || *p == '\\')
            last_sep = p;
    }
    if (!last_sep || last_sep == filepath) return 0; /* no parent or root */

    /* Copy parent path into mutable buffer */
    size_t parent_len = (size_t)(last_sep - filepath);
    if (parent_len >= 1024) return -1;

    char buf[1024];
    memcpy(buf, filepath, parent_len);
    buf[parent_len] = '\0';

    /* Walk path components, creating each directory */
    for (size_t i = 0; i <= parent_len; i++) {
        if (buf[i] == '/' || buf[i] == '\\' || buf[i] == '\0') {
            if (i == 0) continue; /* skip leading separator */
            char saved = buf[i];
            buf[i] = '\0';

            /* Skip drive letters on Windows (e.g., "C:") */
            if (i == 2 && buf[1] == ':') {
                buf[i] = saved;
                continue;
            }

            int rc = ves_mkdir(buf);
            if (rc != 0) {
                /* Ignore EEXIST — directory already exists */
#ifdef _WIN32
                DWORD err = GetLastError();
                if (err != ERROR_ALREADY_EXISTS) {
                    buf[i] = saved;
                    return -1;
                }
#else
                if (errno != EEXIST) {
                    buf[i] = saved;
                    return -1;
                }
#endif
            }
            buf[i] = saved;
        }
    }
    return 0;
}

/* ================================================================
 * Subprocess execution
 * ================================================================ */

#ifdef _WIN32

/* ---- Windows: CreateProcess ---- */

int ves_run_subprocess(const char *exe, const char *const *argv,
                       double timeout_sec)
{
    return ves_run_subprocess_logged(exe, argv, timeout_sec, NULL);
}

int ves_run_subprocess_logged(const char *exe, const char *const *argv,
                              double timeout_sec, const char *log_path)
{
    if (!exe || !argv) return -1;

    /* Build command line string from argv.
     * argv[0] is the program name, argv[1..] are arguments.
     * Quote each argument for safety. */
    char cmdline[4096] = {0};
    size_t pos = 0;
    for (int i = 0; argv[i] != NULL; i++) {
        if (i > 0 && pos < sizeof(cmdline) - 1) {
            cmdline[pos++] = ' ';
        }
        /* Quote the argument */
        if (pos < sizeof(cmdline) - 1) cmdline[pos++] = '"';
        size_t alen = strlen(argv[i]);
        if (pos + alen < sizeof(cmdline) - 2) {
            memcpy(cmdline + pos, argv[i], alen);
            pos += alen;
        }
        if (pos < sizeof(cmdline) - 1) cmdline[pos++] = '"';
    }
    cmdline[pos] = '\0';

    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    memset(&si, 0, sizeof(si));
    si.cb = sizeof(si);
    /* stdout/stderr go to the caller's log file (shared-read so it can be
     * tailed live), or to NUL when no log is asked for (the Poisson case) */
    HANDLE hNull = INVALID_HANDLE_VALUE;
    if (log_path != NULL) {
        SECURITY_ATTRIBUTES sa;
        memset(&sa, 0, sizeof(sa));
        sa.nLength = sizeof(sa);
        sa.bInheritHandle = TRUE;
        hNull = CreateFileA(log_path, FILE_APPEND_DATA,
                            FILE_SHARE_READ | FILE_SHARE_WRITE, &sa,
                            OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    }
    if (hNull == INVALID_HANDLE_VALUE)
        hNull = CreateFileA("NUL", GENERIC_WRITE, FILE_SHARE_WRITE,
                            NULL, OPEN_EXISTING, 0, NULL);
    if (hNull != INVALID_HANDLE_VALUE) {
        si.dwFlags = STARTF_USESTDHANDLES;
        si.hStdOutput = hNull;
        si.hStdError  = hNull;
        si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
        SetHandleInformation(hNull, HANDLE_FLAG_INHERIT, HANDLE_FLAG_INHERIT);
    }
    memset(&pi, 0, sizeof(pi));

    BOOL ok = CreateProcessA(
        exe,       /* lpApplicationName */
        cmdline,   /* lpCommandLine */
        NULL, NULL,
        (hNull != INVALID_HANDLE_VALUE) ? TRUE : FALSE, /* bInheritHandles */
        CREATE_NO_WINDOW, /* logged batch helpers do not need a console window */
        NULL, NULL,
        &si, &pi);

    if (hNull != INVALID_HANDLE_VALUE) CloseHandle(hNull);

    if (!ok) {
        fprintf(stderr, "  [ves_platform] CreateProcess failed: %lu\n",
                (unsigned long)GetLastError());
        return -1;
    }

    /* Wait with timeout */
    DWORD wait_ms = (timeout_sec > 0)
                    ? (DWORD)(timeout_sec * 1000.0)
                    : INFINITE;
    DWORD wait_result = WaitForSingleObject(pi.hProcess, wait_ms);

    if (wait_result == WAIT_TIMEOUT) {
        fprintf(stderr, "  [ves_platform] Subprocess timed out (%.1fs), terminating\n",
                timeout_sec);
        TerminateProcess(pi.hProcess, 1);
        WaitForSingleObject(pi.hProcess, 5000);
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
        return -1;
    }

    DWORD exit_code = 1;
    GetExitCodeProcess(pi.hProcess, &exit_code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    return (exit_code == 0) ? 0 : -1;
}

#else /* POSIX */

/* ---- Linux: fork + execv + waitpid ---- */

#include <sys/wait.h>
#include <spawn.h>
#include <fcntl.h>
extern char **environ;

int ves_run_subprocess(const char *exe, const char *const *argv,
                       double timeout_sec)
{
    return ves_run_subprocess_logged(exe, argv, timeout_sec, NULL);
}

int ves_run_subprocess_logged(const char *exe, const char *const *argv,
                              double timeout_sec, const char *log_path)
{
    if (!exe || !argv) return -1;

    /* Workers may launch from OpenMP threads. posix_spawn avoids calling
     * stdio/allocation routines in a forked, formerly multithreaded child. */
    pid_t pid = 0;
    posix_spawn_file_actions_t actions;
    if (posix_spawn_file_actions_init(&actions) != 0) return -1;
    const char *sink = log_path != NULL ? log_path : "/dev/null";
    int error = posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO,
                                                sink, O_WRONLY | O_CREAT | O_APPEND, 0666);
    if (!error) error = posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);
    if (!error) error = posix_spawn(&pid, exe, &actions, NULL, (char *const *)argv, environ);
    posix_spawn_file_actions_destroy(&actions);
    if (error) {
        fprintf(stderr, "  [ves_platform] posix_spawn failed: %d\n", error);
        return -1;
    }

    /* Parent: poll with timeout */
    double start = ves_clock_sec();
    int status = 0;
    for (;;) {
        pid_t w = waitpid(pid, &status, WNOHANG);
        if (w > 0) break;
        if (w < 0 && errno != EINTR) return -1;

        double elapsed = ves_clock_sec() - start;
        if (timeout_sec > 0 && elapsed > timeout_sec) {
            fprintf(stderr, "  [ves_platform] Subprocess timed out (%.1fs), killing PID %d\n",
                    elapsed, (int)pid);
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            return -1;
        }

        /* Sleep 100ms between polls */
        ves_sleep_ms(100);
    }

    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        return -1;
    }
    return 0;
}

#endif /* _WIN32 */

/* ================================================================
 * Hard timeout
 * ================================================================ */

#ifdef _WIN32

/* Windows: background thread that sets the flag after a delay */

typedef struct {
    DWORD             delay_ms;
    volatile sig_atomic_t *flag;
} TimeoutCtx;

static HANDLE g_timeout_thread = NULL;
static TimeoutCtx g_timeout_ctx = {0, NULL};

static DWORD WINAPI timeout_thread_func(LPVOID param)
{
    TimeoutCtx *ctx = (TimeoutCtx *)param;
    Sleep(ctx->delay_ms);
    if (ctx->flag) {
        *(ctx->flag) = 1;
    }
    return 0;
}

void ves_hard_timeout_start(double seconds, volatile sig_atomic_t *flag)
{
    ves_hard_timeout_cancel();

    g_timeout_ctx.delay_ms = (DWORD)(seconds * 1000.0);
    g_timeout_ctx.flag = flag;
    *flag = 0;

    g_timeout_thread = CreateThread(
        NULL, 0, timeout_thread_func, &g_timeout_ctx, 0, NULL);
}

void ves_hard_timeout_cancel(void)
{
    if (g_timeout_thread) {
        /* We can't cleanly cancel a sleeping thread, but we can
         * terminate it since it's just a simple sleeper. */
        TerminateThread(g_timeout_thread, 0);
        CloseHandle(g_timeout_thread);
        g_timeout_thread = NULL;
    }
}

VES_NORETURN void ves_fast_exit(int code)
{
    /* Lock-free flushes: a thread terminated while holding a stream lock
     * (Triangle prints to stdout) would deadlock fflush(NULL).  stderr is
     * unbuffered; every file this tool writes is closed before this point. */
    _fflush_nolock(stdout);
    _fflush_nolock(stderr);
    TerminateProcess(GetCurrentProcess(), (UINT)code);
    ExitProcess((UINT)code);            /* not reached */
}

static unsigned ves_thread_cpu_ms(HANDLE h)
{
    FILETIME c, e, k, u;
    if (!GetThreadTimes(h, &c, &e, &k, &u)) return 0;
    ULONGLONG kt = ((ULONGLONG)k.dwHighDateTime << 32) | k.dwLowDateTime;
    ULONGLONG ut = ((ULONGLONG)u.dwHighDateTime << 32) | u.dwLowDateTime;
    return (unsigned)((kt + ut) / 10000ULL);          /* 100 ns units -> ms */
}

/* The worker runs fn inside this wrapper so the watchdog knows where it is:
 * 0 = still starting (DLL_THREAD_ATTACH runs under the LOADER LOCK), 1 = in fn,
 * 2 = fn returned (thread exit and DLL_THREAD_DETACH, again under the loader
 * lock).  It is terminated only in state 1, checked again while suspended, so
 * no kill can orphan the loader lock (the PHerc0211 exit-hang stacks showed
 * exactly that: ExitProcess and a starting thread both waiting on a loader
 * lock no live thread held, after 27 wall-clock kills of starved workers). */
typedef struct {
    VesThreadFn   fn;
    void         *arg;
    volatile LONG state;
} VesBoundedJob;

static unsigned __stdcall ves_bounded_entry(void *p)
{
    VesBoundedJob *j = (VesBoundedJob *)p;
    InterlockedExchange(&j->state, 1);
    unsigned rc = j->fn(j->arg);
    InterlockedExchange(&j->state, 2);
    return rc;
}

int ves_run_thread_cpu_bounded(VesThreadFn fn, void *arg, unsigned cpu_ms,
                               unsigned wall_ms, unsigned *out_cpu_ms)
{
    VesBoundedJob job;
    job.fn = fn; job.arg = arg; job.state = 0;
    uintptr_t h = _beginthreadex(NULL, 0, ves_bounded_entry, &job, 0, NULL);
    if (h == 0) return -1;
    ULONGLONG t_run = 0;                  /* when fn was first seen running */
    for (;;) {
        DWORD w = WaitForSingleObject((HANDLE)h, 100);
        unsigned used = ves_thread_cpu_ms((HANDLE)h);
        if (out_cpu_ms) *out_cpu_ms = used;
        if (w == WAIT_OBJECT_0) {
            CloseHandle((HANDLE)h);
            return 0;
        }
        if (InterlockedCompareExchange(&job.state, 0, 0) != 1)
            continue;                     /* starting up or exiting: never kill */
        if (t_run == 0) t_run = GetTickCount64();
        int why = 0;
        if (used >= cpu_ms) why = 1;
        else if (GetTickCount64() - t_run >= (ULONGLONG)wall_ms) why = 2;
        if (!why) continue;
        if (SuspendThread((HANDLE)h) == (DWORD)-1) continue;
        if (InterlockedCompareExchange(&job.state, 0, 0) == 1) {
            TerminateThread((HANDLE)h, 1);
            WaitForSingleObject((HANDLE)h, 2000);
            CloseHandle((HANDLE)h);
            return why;
        }
        ResumeThread((HANDLE)h);          /* fn returned meanwhile: let it exit */
    }
}

static volatile LONG g_st_spin_stop = 0;
static unsigned __stdcall st_quick(void *p) { *(int *)p = 7; return 0; }
static unsigned __stdcall st_spin(void *p)
{
    volatile unsigned long long x = 0;
    (void)p;
    while (!g_st_spin_stop) x++;
    return 0;
}
static unsigned __stdcall st_sleep(void *p) { Sleep(*(DWORD *)p); return 0; }

int ves_platform_selftest(void)
{
    int fails = 0, flag = 0;
    unsigned cpu = 0;
    DWORD nap = 600, block = 5000;
    int r0 = ves_run_thread_cpu_bounded(st_quick, &flag, 300, 10000, &cpu);
    if (r0 != 0 || flag != 7) { fprintf(stderr, "ves_platform_selftest FAIL: quick %d\n", r0); fails++; }
    int r1 = ves_run_thread_cpu_bounded(st_spin, NULL, 300, 60000, &cpu);
    if (r1 != 1 || cpu < 300) { fprintf(stderr, "ves_platform_selftest FAIL: spin %d (%u ms)\n", r1, cpu); fails++; }
    /* a worker that waits 600 ms but uses no CPU is NOT killed by a 300 ms CPU bound */
    int r2 = ves_run_thread_cpu_bounded(st_sleep, &nap, 300, 60000, &cpu);
    if (r2 != 0) { fprintf(stderr, "ves_platform_selftest FAIL: sleeper killed %d\n", r2); fails++; }
    int r3 = ves_run_thread_cpu_bounded(st_sleep, &block, 300, 400, &cpu);
    if (r3 != 2) { fprintf(stderr, "ves_platform_selftest FAIL: wall backstop %d\n", r3); fails++; }
    fprintf(stderr, "ves_platform_selftest: %s (CPU-bounded worker: quick/spin/sleep/blocked)\n",
            fails ? "FAIL" : "OK");
    return fails;
}

#else /* POSIX */

/* Linux: SIGALRM + alarm() */

static volatile sig_atomic_t *g_timeout_flag_ptr = NULL;

static void ves_alarm_handler(int sig)
{
    (void)sig;
    if (g_timeout_flag_ptr) {
        *g_timeout_flag_ptr = 1;
    }
}

void ves_hard_timeout_start(double seconds, volatile sig_atomic_t *flag)
{
    g_timeout_flag_ptr = flag;
    *flag = 0;

    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = ves_alarm_handler;
    sigaction(SIGALRM, &sa, NULL);
    alarm((unsigned int)(seconds + 1.0));
}

void ves_hard_timeout_cancel(void)
{
    alarm(0);
    g_timeout_flag_ptr = NULL;
}

VES_NORETURN void ves_fast_exit(int code)
{
    fflush(stdout);
    fflush(stderr);
    _exit(code);
}

int ves_run_thread_cpu_bounded(VesThreadFn fn, void *arg, unsigned cpu_ms,
                               unsigned wall_ms, unsigned *out_cpu_ms)
{
    (void)fn; (void)arg; (void)cpu_ms; (void)wall_ms;
    if (out_cpu_ms) *out_cpu_ms = 0;
    return -1;                       /* no bounded thread here: run inline */
}

int ves_platform_selftest(void)
{
    return 0;
}

#endif /* _WIN32 */

/* ----------------------------------------------------------------
 * Thread count.  Declared in ves_platform.h; called by the public API
 * (sf_common_opts.n_threads).  Compiles to a no-op without OpenMP so
 * the caller does not need a conditional of its own.
 * ---------------------------------------------------------------- */
/* ves_omp.h declares the OpenMP entry points instead of including <omp.h>;
 * see that header for why. Without _OPENMP its omp_set_num_threads is a no-op,
 * so no conditional is needed here. */
#include "ves_omp.h"

void ves_omp_set_threads(int n)
{
    if (n <= 0) return;          /* keep the runtime's own default */
    omp_set_num_threads(n);
}
