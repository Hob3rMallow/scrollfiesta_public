/*
 * ves_platform.c — Platform abstraction: subprocess + hard timeout.
 *
 * Contains the two functions too large for inline in ves_platform.h:
 *   1. ves_run_subprocess()      — fork/exec (Linux) or CreateProcess (Windows)
 *   2. ves_hard_timeout_start()  — SIGALRM (Linux) or timer thread (Windows)
 *      ves_hard_timeout_cancel()
 */
#include "ves_platform.h"

#include <string.h>
#include <stdlib.h>
#include <errno.h>

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
    /* Redirect stdout/stderr to NUL to suppress Poisson output */
    HANDLE hNull = CreateFileA("NUL", GENERIC_WRITE, FILE_SHARE_WRITE,
                               NULL, OPEN_EXISTING, 0, NULL);
    if (hNull != INVALID_HANDLE_VALUE) {
        si.dwFlags = STARTF_USESTDHANDLES;
        si.hStdOutput = hNull;
        si.hStdError  = hNull;
        si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
    }
    memset(&pi, 0, sizeof(pi));

    BOOL ok = CreateProcessA(
        exe,       /* lpApplicationName */
        cmdline,   /* lpCommandLine */
        NULL, NULL,
        (hNull != INVALID_HANDLE_VALUE) ? TRUE : FALSE, /* bInheritHandles */
        0,         /* dwCreationFlags */
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

int ves_run_subprocess(const char *exe, const char *const *argv,
                       double timeout_sec)
{
    if (!exe || !argv) return -1;

    pid_t pid = fork();
    if (pid < 0) {
        fprintf(stderr, "  [ves_platform] fork() failed\n");
        return -1;
    }

    if (pid == 0) {
        /* Child: redirect stdout/stderr to /dev/null */
        FILE *devnull = fopen("/dev/null", "w");
        if (devnull) {
            dup2(fileno(devnull), STDOUT_FILENO);
            dup2(fileno(devnull), STDERR_FILENO);
            fclose(devnull);
        }
        /* execv expects char *const *, cast away the outer const */
        execv(exe, (char *const *)argv);
        _exit(127); /* exec failed */
    }

    /* Parent: poll with timeout */
    double start = ves_clock_sec();
    int status = 0;
    for (;;) {
        pid_t w = waitpid(pid, &status, WNOHANG);
        if (w > 0) break;

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

#endif /* _WIN32 */
