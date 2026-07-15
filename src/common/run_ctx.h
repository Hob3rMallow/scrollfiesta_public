#ifndef RUN_CTX_INCLUDED
#define RUN_CTX_INCLUDED

#include <stddef.h>

/*
 * run_ctx — per-call execution context for the public C API.
 *
 * When ScrollFiesta runs as a library (include/scrollfiesta.h), each API call
 * pushes an SfRunCtx carrying the caller's progress/cancel callback, log sink,
 * deadline, and a per-call key/value tuning table that overrides getenv for
 * the module knobs. Modules interact with it through three tiny hooks:
 *
 *   sf_env(name)              — read a tuning knob: per-call table first,
 *                               then the real environment. Drop-in for
 *                               getenv() in library code; the CLIs keep raw
 *                               getenv semantics because no context is
 *                               pushed there (table empty -> plain getenv).
 *   RunCtx_should_stop()      — cheap poll: cancel requested or deadline
 *                               passed. Safe anywhere, including OpenMP
 *                               regions and C++ TUs (never longjmps).
 *   RunCtx_check()            — poll and RAISE(Sf_Cancelled) back to the API
 *                               boundary. NEVER call inside an OpenMP
 *                               parallel region or from C++ code (longjmp
 *                               across either is undefined behaviour) —
 *                               poll with RunCtx_should_stop() there and
 *                               raise after.
 *
 * The context is thread-local: concurrent callers on different threads see
 * only their own context. With no context pushed (all CLI tools), every hook
 * degrades to a no-op / plain getenv.
 */

typedef int  (*SfRunProgressFn)(void *user, const char *stage, double fraction);
typedef void (*SfRunLogFn)(void *user, const char *line);

typedef struct SfRunKV {
    const char *key;
    const char *value;
} SfRunKV;

typedef struct SfRunCtx {
    SfRunProgressFn progress;       /* NULL = no progress reporting        */
    void           *progress_user;
    SfRunLogFn      log;            /* NULL = library keeps using stderr   */
    void           *log_user;
    volatile int    cancel_requested;
    double          deadline_sec;   /* absolute ves_clock_sec(); 0 = none  */
    int             timed_out;      /* set once the deadline fires         */
    const SfRunKV  *tuning;         /* per-call env overrides (may be NULL) */
    size_t          n_tuning;
} SfRunCtx;

/* Install / remove the calling thread's context. Push/pop must pair. */
void      RunCtx_push(SfRunCtx *ctx);
void      RunCtx_pop(void);
SfRunCtx *RunCtx_current(void);   /* NULL when no API call is active */

/* Poll for cancellation or deadline. Returns nonzero if the operation should
 * stop. Sets ctx->timed_out when the deadline (not a cancel) fired. */
int  RunCtx_should_stop(void);

/* Poll and raise Sf_Cancelled toward the API boundary when stopping.
 * longjmp-based: see the header comment for where this must NOT be called. */
void RunCtx_check(void);

/* Report progress (stage id + fraction in [0,1]) and poll the callback's
 * cancel request. Returns nonzero if the operation should stop; never raises. */
int  RunCtx_progress(const char *stage, double fraction);

/* RunCtx_progress + RAISE(Sf_Cancelled) when stopping. Same longjmp caveats
 * as RunCtx_check. */
void RunCtx_progress_check(const char *stage, double fraction);

/* printf-style line to the context's log sink, or stderr when none. */
void RunCtx_log(const char *fmt, ...);

/* getenv() with the per-call tuning table consulted first. Library modules
 * call this instead of getenv so API callers can pass knobs per call without
 * touching the process environment. */
const char *sf_env(const char *name);

#endif /* RUN_CTX_INCLUDED */
