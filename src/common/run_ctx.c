/*
 * run_ctx.c — per-call execution context (see run_ctx.h).
 */
#include "run_ctx.h"
#include "except.h"
#include "ves_platform.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Thread-local so concurrent API calls on different threads are isolated.
 * Same keyword strategy as Except_stack (see except.h). */
static EXCEPT_THREAD_LOCAL SfRunCtx *g_run_ctx = NULL;

void RunCtx_push(SfRunCtx *ctx)
{
    g_run_ctx = ctx;
}

void RunCtx_pop(void)
{
    g_run_ctx = NULL;
}

SfRunCtx *RunCtx_current(void)
{
    return g_run_ctx;
}

int RunCtx_should_stop(void)
{
    SfRunCtx *ctx = g_run_ctx;
    if (!ctx)
        return 0;
    if (ctx->cancel_requested)
        return 1;
    if (ctx->deadline_sec > 0.0 && ves_clock_sec() > ctx->deadline_sec) {
        ctx->timed_out = 1;
        return 1;
    }
    return 0;
}

void RunCtx_check(void)
{
    if (RunCtx_should_stop())
        RAISE(Sf_Cancelled);
}

int RunCtx_progress(const char *stage, double fraction)
{
    SfRunCtx *ctx = g_run_ctx;
    if (!ctx)
        return 0;
    if (ctx->progress) {
        if (fraction < 0.0) fraction = 0.0;
        if (fraction > 1.0) fraction = 1.0;
        if (ctx->progress(ctx->progress_user, stage, fraction) != 0)
            ctx->cancel_requested = 1;
    }
    return RunCtx_should_stop();
}

void RunCtx_progress_check(const char *stage, double fraction)
{
    if (RunCtx_progress(stage, fraction))
        RAISE(Sf_Cancelled);
}

void RunCtx_log(const char *fmt, ...)
{
    char    buf[1024];
    va_list ap;
    SfRunCtx *ctx = g_run_ctx;

    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    buf[sizeof buf - 1] = '\0';

    if (ctx && ctx->log)
        ctx->log(ctx->log_user, buf);
    else
        fprintf(stderr, "%s\n", buf);
}

const char *sf_env(const char *name)
{
    SfRunCtx *ctx = g_run_ctx;
    if (ctx && ctx->tuning) {
        for (size_t i = 0; i < ctx->n_tuning; i++) {
            if (ctx->tuning[i].key && strcmp(ctx->tuning[i].key, name) == 0)
                return ctx->tuning[i].value;
        }
    }
    return getenv(name);
}
