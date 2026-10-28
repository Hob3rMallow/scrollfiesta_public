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
#include <winioctl.h>  /* FSCTL_SET_REPARSE_POINT (the scratch selftest's junction) */
#endif

static int ves_scratch_selftest(void);   /* end of file */

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

/* ---- VES_PROF_OUT: in-process sampling profiler (diagnostics only) ----
 * VES_PROF_OUT=<file> starts a background thread that, about 1000 times a
 * second, suspends every other thread of the process, records its instruction
 * pointer (GetThreadContext) and resumes it.  While a thread is suspended the
 * sampler touches only a table allocated up front with VirtualAlloc, so it can
 * never wait on a lock the suspended thread holds.  About once a second (with
 * no thread suspended) it rewrites <file> with the module map and the
 * per-address counts, so a TerminateProcess exit loses at most a second.
 * Leaf (self) samples only; scripts/pipeline/ves_prof_report.py resolves them
 * against the PDBs.  Output-neutral: nothing else reads the table. */
#include <tlhelp32.h>
#include <psapi.h>

#define VES_PROF_TAB 65536u
typedef struct { ULONG64 rip; ULONG32 count; ULONG32 pad; } VesProfSlot;
static VesProfSlot *g_prof_tab = NULL;
static ULONG64 g_prof_ticks = 0, g_prof_dropped = 0;
static volatile LONG g_prof_stop = 0;
static LONGLONG g_prof_due = -10000;             /* sample interval, 100 ns units (VES_PROF_MS, default 1 ms) */           /* selftest only: ends the sampler loop */
static HANDLE g_prof_thread = NULL;
static char g_prof_path[MAX_PATH];

static void ves_prof_add(ULONG64 rip)
{
    ULONG32 h = (ULONG32)((rip * 0x9E3779B97F4A7C15ull) >> 48);
    for (ULONG32 k = 0; k < 64; k++) {
        VesProfSlot *s = &g_prof_tab[(h + k) & (VES_PROF_TAB - 1)];
        if (s->rip == rip) { s->count++; return; }
        if (s->rip == 0) { s->rip = rip; s->count = 1; return; }
    }
    g_prof_dropped++;
}

static volatile LONG g_prof_flushing = 0;

static void ves_prof_flush(void)
{
    if (InterlockedExchange(&g_prof_flushing, 1)) return;   /* one writer at a time */
    char tmp[MAX_PATH + 8];
    snprintf(tmp, sizeof(tmp), "%s.tmp", g_prof_path);
    FILE *f = fopen(tmp, "w");
    if (!f) { InterlockedExchange(&g_prof_flushing, 0); return; }
    HMODULE mods[1024];
    DWORD need = 0;
    HANDLE proc = GetCurrentProcess();
    if (EnumProcessModules(proc, mods, (DWORD)sizeof(mods), &need)) {
        DWORD nm = need / (DWORD)sizeof(HMODULE);
        if (nm > 1024) nm = 1024;
        for (DWORD i = 0; i < nm; i++) {
            MODULEINFO mi;
            char name[MAX_PATH];
            if (!GetModuleInformation(proc, mods[i], &mi, (DWORD)sizeof(mi))) continue;
            if (!GetModuleFileNameA(mods[i], name, MAX_PATH)) continue;
            fprintf(f, "module %llx %lx %s\n", (unsigned long long)(ULONG_PTR)mi.lpBaseOfDll,
                    (unsigned long)mi.SizeOfImage, name);
        }
    }
    fprintf(f, "ticks %llu dropped %llu\n", (unsigned long long)g_prof_ticks,
            (unsigned long long)g_prof_dropped);
    for (ULONG32 k = 0; k < VES_PROF_TAB; k++)
        if (g_prof_tab[k].rip)
            fprintf(f, "%llx %lu\n", (unsigned long long)g_prof_tab[k].rip, (unsigned long)g_prof_tab[k].count);
    fclose(f);
    MoveFileExA(tmp, g_prof_path, MOVEFILE_REPLACE_EXISTING);
    InterlockedExchange(&g_prof_flushing, 0);
}

/* The thread list.  CreateToolhelp32Snapshot lists every thread on the
 * machine (19,101 threads, 0.26 s on the busy box, 2026-10-05), so a helper
 * thread refreshes the list every 2 s and the sampler never waits for it.
 * The sampler holds the lock SHARED for a whole pass -- the helper swaps (and
 * closes the old handles) only under the EXCLUSIVE lock, so no handle is
 * closed while it is being suspended -- and neither profiler thread is in the
 * list, so the sampler never suspends a thread that could hold the lock. */
#define VES_PROF_MAXTH 512
static SRWLOCK g_prof_lock = SRWLOCK_INIT;
static HANDLE g_prof_th[VES_PROF_MAXTH];
static int g_prof_nth = 0;
static volatile DWORD g_prof_sampler_tid = 0, g_prof_lister_tid = 0;
static HANDLE g_prof_lister = NULL;

static DWORD WINAPI ves_prof_lister(LPVOID unused)
{
    (void)unused;
    const DWORD pid = GetCurrentProcessId();
    g_prof_lister_tid = GetCurrentThreadId();
    while (!g_prof_stop) {
        HANDLE fresh[VES_PROF_MAXTH];
        int nf = 0;
        HANDLE snap = CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD, 0);
        if (snap != INVALID_HANDLE_VALUE) {
            THREADENTRY32 te;
            te.dwSize = (DWORD)sizeof(te);
            for (BOOL ok = Thread32First(snap, &te); ok && nf < VES_PROF_MAXTH; ok = Thread32Next(snap, &te)) {
                if (te.th32OwnerProcessID != pid || te.th32ThreadID == g_prof_sampler_tid ||
                    te.th32ThreadID == g_prof_lister_tid) continue;
                HANDLE h = OpenThread(THREAD_SUSPEND_RESUME | THREAD_GET_CONTEXT, FALSE, te.th32ThreadID);
                if (h) fresh[nf++] = h;
            }
            CloseHandle(snap);
        }
        HANDLE old[VES_PROF_MAXTH];
        AcquireSRWLockExclusive(&g_prof_lock);
        const int nold = g_prof_nth;
        memcpy(old, g_prof_th, (size_t)nold * sizeof(HANDLE));
        memcpy(g_prof_th, fresh, (size_t)nf * sizeof(HANDLE));
        g_prof_nth = nf;
        ReleaseSRWLockExclusive(&g_prof_lock);
        for (int i = 0; i < nold; i++) CloseHandle(old[i]);
        for (int w = 0; w < 20 && !g_prof_stop; w++) Sleep(100);
    }
    return 0;
}

static DWORD WINAPI ves_prof_thread(LPVOID unused)
{
    (void)unused;
    g_prof_sampler_tid = GetCurrentThreadId();
    ULONGLONG last_flush = GetTickCount64();
    HANDLE timer = CreateWaitableTimerExW(NULL, NULL, 0x00000002 /* HIGH_RESOLUTION */, TIMER_ALL_ACCESS);
    while (!g_prof_stop) {
        ULONGLONG now = GetTickCount64();
        AcquireSRWLockShared(&g_prof_lock);
        for (int i = 0; i < g_prof_nth; i++) {
            if (SuspendThread(g_prof_th[i]) == (DWORD)-1) continue;
            CONTEXT ctx;
            memset(&ctx, 0, sizeof(ctx));
            ctx.ContextFlags = CONTEXT_CONTROL;
            if (GetThreadContext(g_prof_th[i], &ctx)) ves_prof_add(ctx.Rip);
            ResumeThread(g_prof_th[i]);
        }
        ReleaseSRWLockShared(&g_prof_lock);
        g_prof_ticks++;
        if (now - last_flush >= 1000) { ves_prof_flush(); last_flush = now; }
        if (timer) {
            LARGE_INTEGER due;
            due.QuadPart = g_prof_due;        /* relative */
            if (SetWaitableTimer(timer, &due, 0, NULL, NULL, FALSE)) WaitForSingleObject(timer, INFINITE);
            else Sleep((DWORD)(-g_prof_due / 10000));
        } else {
            Sleep((DWORD)(-g_prof_due / 10000));
        }
    }
    if (timer) CloseHandle(timer);
    return 0;
}

static void ves_prof_flush_atexit(void);

void ves_prof_start_from_env(void)
{
    const char *p = getenv("VES_PROF_OUT");
    if (!p || !*p || g_prof_tab) return;
    snprintf(g_prof_path, sizeof(g_prof_path), "%s", p);
    {   /* VES_PROF_MS=<n>: one sample every n ms (default 1; 10 for long multi-threaded runs) */
        const char *ms = getenv("VES_PROF_MS");
        long v = ms ? strtol(ms, NULL, 10) : 1;
        if (v < 1) v = 1;
        if (v > 1000) v = 1000;
        g_prof_due = -(LONGLONG)v * 10000;
    }
    g_prof_tab = (VesProfSlot *)VirtualAlloc(NULL, (SIZE_T)VES_PROF_TAB * sizeof(VesProfSlot),
                                             MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
    if (!g_prof_tab) return;
    g_prof_stop = 0;
    g_prof_thread = CreateThread(NULL, 0, ves_prof_thread, NULL, 0, NULL);
    if (g_prof_thread) SetThreadPriority(g_prof_thread, THREAD_PRIORITY_HIGHEST);
    while (g_prof_thread && !g_prof_sampler_tid) Sleep(0);   /* its id must be known before the first list */
    g_prof_lister = CreateThread(NULL, 0, ves_prof_lister, NULL, 0, NULL);
    atexit(ves_prof_flush_atexit);   /* tools that return from main (sheet_assemble) */
}

void ves_prof_flush_now(void)
{
    if (g_prof_tab) ves_prof_flush();
}

static void ves_prof_flush_atexit(void)
{
    if (g_prof_tab) ves_prof_flush();
}

/* Selftest: sample this process for ~1 s while one thread spins, stop the
 * sampler, and check the file has ticks, the module map and samples. */
static unsigned __stdcall prof_st_spin(void *p)
{
    volatile double x = 1.0;
    ULONGLONG t0 = GetTickCount64();
    while (GetTickCount64() - t0 < 1000) x = x * 1.0000001 + 1e-9;
    *(double *)p = x;
    return 0;
}

static int ves_prof_selftest(void)
{
    if (g_prof_tab) return 0;                    /* a real profile is running: leave it alone */
    char dir[MAX_PATH], path[MAX_PATH];
    if (!GetTempPathA(MAX_PATH, dir)) return 1;
    snprintf(path, sizeof(path), "%sves_prof_selftest_%lu.prof", dir, (unsigned long)GetCurrentProcessId());
    if (!SetEnvironmentVariableA("VES_PROF_OUT", path)) return 1;
    _putenv_s("VES_PROF_OUT", path);
    ves_prof_start_from_env();
    double sink = 0.0;
    uintptr_t h = _beginthreadex(NULL, 0, prof_st_spin, &sink, 0, NULL);
    if (h) { WaitForSingleObject((HANDLE)h, INFINITE); CloseHandle((HANDLE)h); }
    g_prof_stop = 1;
    if (g_prof_thread) { WaitForSingleObject(g_prof_thread, 5000); CloseHandle(g_prof_thread); g_prof_thread = NULL; }
    if (g_prof_lister) { WaitForSingleObject(g_prof_lister, 5000); CloseHandle(g_prof_lister); g_prof_lister = NULL; }
    ves_prof_flush();
    for (int i = 0; i < g_prof_nth; i++) CloseHandle(g_prof_th[i]);
    g_prof_nth = 0; g_prof_sampler_tid = 0; g_prof_lister_tid = 0;
    _putenv_s("VES_PROF_OUT", "");
    SetEnvironmentVariableA("VES_PROF_OUT", NULL);
    FILE *f = fopen(path, "r");
    int modules = 0, samples = 0;
    unsigned long long ticks = 0;
    char line[1024];
    while (f && fgets(line, sizeof(line), f)) {
        if (!strncmp(line, "module ", 7)) modules++;
        else if (!strncmp(line, "ticks ", 6)) ticks = strtoull(line + 6, NULL, 10);
        else samples++;
    }
    if (f) fclose(f);
    DeleteFileA(path);
    VirtualFree(g_prof_tab, 0, MEM_RELEASE);
    g_prof_tab = NULL; g_prof_ticks = 0; g_prof_dropped = 0;
    int ok = modules > 0 && ticks > 0 && samples > 0;
    fprintf(stderr, "ves_prof_selftest: %s (%llu ticks, %d modules, %d sampled addresses)\n",
            ok ? "OK" : "FAIL", ticks, modules, samples);
    return ok ? 0 : 1;
}

VES_NORETURN void ves_fast_exit(int code)
{
    ves_prof_flush_now();               /* VES_PROF_OUT: last second of samples */
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
    if (ves_prof_selftest() != 0) fails++;
    if (ves_scratch_selftest() != 0) fails++;
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
    return ves_scratch_selftest();
}

void ves_prof_start_from_env(void) {}
void ves_prof_flush_now(void) {}

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

/* ================================================================
 * Selftest of the selftest scratch helpers (ves_selftest_dir/_clean in
 * ves_platform.h): a nested tree holding a read-only file and a junction
 * (a symlink on POSIX) to a directory OUTSIDE the tree is removed, the
 * link's target survives, and paths outside selftest/ are refused.
 * ================================================================ */

#ifdef _WIN32
/* A directory junction `link` -> `target`; junctions need no privilege,
 * unlike symlinks.  0 on success. */
static int st_junction(const char *link, const char *target)
{
    WCHAR wt[MAX_PATH], full[MAX_PATH];
    struct {                                  /* REPARSE_DATA_BUFFER, mount point */
        DWORD tag;
        WORD data_len, reserved, sub_off, sub_len, print_off, print_len;
        WCHAR buf[2 * MAX_PATH + 8];
    } rb;
    DWORD full_n = 0, bytes = 0;
    size_t sub_chars = 0;
    HANDLE h = INVALID_HANDLE_VALUE;
    int rc = -1;
    if (MultiByteToWideChar(CP_ACP, 0, target, -1, wt, MAX_PATH) == 0) return -1;
    full_n = GetFullPathNameW(wt, MAX_PATH, full, NULL);
    if (full_n == 0 || full_n >= MAX_PATH) return -1;
    if (!CreateDirectoryA(link, NULL)) return -1;
    memset(&rb, 0, sizeof rb);
    rb.tag = IO_REPARSE_TAG_MOUNT_POINT;
    memcpy(rb.buf, L"\\??\\", 4 * sizeof(WCHAR));
    memcpy(rb.buf + 4, full, (size_t)full_n * sizeof(WCHAR));
    sub_chars = 4 + (size_t)full_n;
    memcpy(rb.buf + sub_chars + 1, full, (size_t)full_n * sizeof(WCHAR));
    rb.sub_len = (WORD)(sub_chars * sizeof(WCHAR));
    rb.print_off = (WORD)((sub_chars + 1) * sizeof(WCHAR));
    rb.print_len = (WORD)((size_t)full_n * sizeof(WCHAR));
    rb.data_len = (WORD)(8 + rb.sub_len + sizeof(WCHAR) + rb.print_len + sizeof(WCHAR));
    h = CreateFileA(link, GENERIC_WRITE, 0, NULL, OPEN_EXISTING,
                    FILE_FLAG_OPEN_REPARSE_POINT | FILE_FLAG_BACKUP_SEMANTICS, NULL);
    if (h != INVALID_HANDLE_VALUE) {
        if (DeviceIoControl(h, FSCTL_SET_REPARSE_POINT, &rb, (DWORD)(8 + rb.data_len),
                            NULL, 0, &bytes, NULL))
            rc = 0;
        CloseHandle(h);
    }
    return rc;
}

static int st_exists(const char *path)
{
    return GetFileAttributesA(path) != INVALID_FILE_ATTRIBUTES;
}
#else
static int st_exists(const char *path)
{
    struct stat st;
    return lstat(path, &st) == 0;
}
#endif

static int st_write(const char *path, const char *text)
{
    FILE *f = NULL;
    if (ves_ensure_parent_dir(path) != 0) return -1;
    f = fopen(path, "wb");
    if (!f) return -1;
    fputs(text, f);
    return fclose(f) == 0 ? 0 : -1;
}

static int ves_scratch_selftest(void)
{
    int fails = 0, linked = 0;
    char dir[256], keep[256], p[600], link[600], sentinel[600];
    if (ves_selftest_dir(dir, sizeof dir, "ves_platform") != 0
        || ves_selftest_dir(keep, sizeof keep, "ves_platform_keep") != 0) {
        fprintf(stderr, "ves_platform_selftest FAIL: scratch dirs\n");
        return 1;
    }
    snprintf(p, sizeof p, "%s/a/b/c.txt", dir);
    if (st_write(p, "x") != 0) fails++;
    snprintf(p, sizeof p, "%s/ro.txt", dir);
    if (st_write(p, "x") != 0) fails++;
#ifdef _WIN32
    (void)SetFileAttributesA(p, FILE_ATTRIBUTE_READONLY);
#endif
    snprintf(sentinel, sizeof sentinel, "%s/sentinel.txt", keep);
    if (st_write(sentinel, "keep") != 0) fails++;
    snprintf(link, sizeof link, "%s/a/link", dir);
#ifdef _WIN32
    linked = st_junction(link, keep) == 0;
#else
    {
        char rel[600];
        snprintf(rel, sizeof rel, "../../%s", ves_path_basename(keep));
        linked = symlink(rel, link) == 0;
    }
#endif
    snprintf(p, sizeof p, "%s/sentinel.txt", link);
    if (!linked || !st_exists(p)) {
        fprintf(stderr, "ves_platform_selftest FAIL: link %s -> %s not made\n", link, keep);
        fails++;
    }
    if (ves_selftest_clean(dir) != 0 || st_exists(dir)) {
        fprintf(stderr, "ves_platform_selftest FAIL: %s not removed\n", dir);
        fails++;
    }
    if (!st_exists(sentinel)) {
        fprintf(stderr, "ves_platform_selftest FAIL: the clean followed the link into %s\n", keep);
        fails++;
    }
    if (ves_selftest_clean(dir) != 0) fails++;                    /* already gone: 0 */
    if (ves_selftest_clean("output/x") != -1 || ves_selftest_clean(VES_SELFTEST_DIR) != -1
        || ves_selftest_clean(VES_SELFTEST_DIR "/") != -1
        || ves_selftest_clean(VES_SELFTEST_DIR "/../output") != -1 || ves_selftest_clean(NULL) != -1) {
        fprintf(stderr, "ves_platform_selftest FAIL: a path outside %s/ was not refused\n",
                VES_SELFTEST_DIR);
        fails++;
    }
    if (ves_selftest_clean(keep) != 0 || st_exists(keep)) fails++;
    fprintf(stderr, "ves_platform_selftest: %s (selftest scratch: nested tree, read-only file, "
            "link not followed, refusals)\n", fails ? "FAIL" : "OK");
    return fails;
}
