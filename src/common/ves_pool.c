#include "ves_pool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* See ves_pool.h.  One queue of runs; a run carries `tickets` helper seats.
 * A helper takes a seat (its slot), runs chunks until the range is used up,
 * and returns to the idle count; the caller cancels the seats nobody took once
 * it has run out of chunks itself, then waits for the helpers that joined. */

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <process.h>
typedef SRWLOCK vp_lock_t;
typedef CONDITION_VARIABLE vp_cond_t;
#define VP_LOCK_INIT SRWLOCK_INIT
#define VP_COND_INIT CONDITION_VARIABLE_INIT
static void vp_lock(vp_lock_t *l) { AcquireSRWLockExclusive(l); }
static void vp_unlock(vp_lock_t *l) { ReleaseSRWLockExclusive(l); }
static void vp_wait(vp_cond_t *c, vp_lock_t *l) { SleepConditionVariableSRW(c, l, INFINITE, 0); }
static void vp_wake_all(vp_cond_t *c) { WakeAllConditionVariable(c); }
static void vp_wake_one(vp_cond_t *c) { WakeConditionVariable(c); }
static long long vp_fetch_add(volatile long long *p, long long d) { return _InterlockedExchangeAdd64((volatile __int64 *)p, d); }
static long vp_add32(volatile long *p, long d) { return _InterlockedExchangeAdd(p, d) + d; }
static unsigned __stdcall vp_thread_main(void *arg);
static int vp_spawn(void) { uintptr_t h = _beginthreadex(NULL, 0, vp_thread_main, NULL, 0, NULL); if (!h) return -1; CloseHandle((HANDLE)h); return 0; }
#else
#include <pthread.h>
typedef pthread_mutex_t vp_lock_t;
typedef pthread_cond_t vp_cond_t;
#define VP_LOCK_INIT PTHREAD_MUTEX_INITIALIZER
#define VP_COND_INIT PTHREAD_COND_INITIALIZER
static void vp_lock(vp_lock_t *l) { pthread_mutex_lock(l); }
static void vp_unlock(vp_lock_t *l) { pthread_mutex_unlock(l); }
static void vp_wait(vp_cond_t *c, vp_lock_t *l) { pthread_cond_wait(c, l); }
static void vp_wake_all(vp_cond_t *c) { pthread_cond_broadcast(c); }
static void vp_wake_one(vp_cond_t *c) { pthread_cond_signal(c); }
static long long vp_fetch_add(volatile long long *p, long long d) { return __sync_fetch_and_add(p, d); }
static long vp_add32(volatile long *p, long d) { return __sync_add_and_fetch(p, d); }
static void *vp_thread_main(void *arg);
static int vp_spawn(void) { pthread_t t; if (pthread_create(&t, NULL, vp_thread_main, NULL)) return -1; pthread_detach(t); return 0; }
#endif

typedef struct VpRun {
    VesPoolFn fn;
    void *ctx;
    size_t n, chunk;
    volatile long long next;  /* start of the next chunk */
    int tickets;              /* seats not yet taken (under the lock) */
    int joined;               /* seats taken (under the lock) */
    int running;              /* helpers still running chunks (under the lock) */
    struct VpRun *link;
} VpRun;

static vp_lock_t vp_mutex = VP_LOCK_INIT;
static vp_cond_t vp_work = VP_COND_INIT, vp_done = VP_COND_INIT;
static VpRun *vp_head, *vp_tail;
static int vp_cores = 1;           /* configured cores (outer jobs + helpers) */
static int vp_started;             /* helper threads created */
static int vp_idle;                /* helpers waiting for a seat (under the lock) */
static volatile long vp_busy;      /* outer jobs running */
static size_t vp_runs, vp_borrowed;

static void vp_chunks(VpRun *r, int slot)
{
    for (;;) {
        long long s = vp_fetch_add(&r->next, (long long)r->chunk);
        if (s < 0 || (size_t)s >= r->n) break;
        size_t e = (size_t)s + r->chunk;
        if (e > r->n) e = r->n;
        r->fn(r->ctx, (size_t)s, e, slot);
    }
}

#ifdef _WIN32
static unsigned __stdcall vp_thread_main(void *arg)
#else
static void *vp_thread_main(void *arg)
#endif
{
    (void)arg;
    vp_lock(&vp_mutex);
    for (;;) {
        while (!vp_head) vp_wait(&vp_work, &vp_mutex);
        VpRun *r = vp_head;
        int slot = 1 + r->joined++;
        r->running++;
        if (--r->tickets == 0) { vp_head = r->link; if (!vp_head) vp_tail = NULL; }
        vp_unlock(&vp_mutex);
        vp_chunks(r, slot);
        vp_lock(&vp_mutex);
        vp_idle++;
        if (--r->running == 0) vp_wake_all(&vp_done);
    }
#ifndef _WIN32
    return NULL;
#endif
}

void VesPool_configure(int cores)
{
    vp_lock(&vp_mutex);
    if (cores < 1) cores = 1;
    /* Helpers once started stay; a smaller configuration only lends less. */
    vp_cores = cores;
    vp_unlock(&vp_mutex);
}

int VesPool_cores(void) { return vp_cores; }
void VesPool_job_enter(void) { vp_add32(&vp_busy, 1); }
void VesPool_job_leave(void) { vp_add32(&vp_busy, -1); }

int VesPool_share(void)
{
    int cores = vp_cores, jobs = (int)vp_busy;
    if (cores <= 1) return 0;
    if (jobs < 1) jobs = 1;
    int share = (cores - jobs) / jobs;
    return share > 0 ? share : 0;
}

int VesPool_run(size_t n, size_t chunk, VesPoolFn fn, void *ctx, int max_helpers)
{
    if (!n) return 1;
    if (chunk < 1) chunk = 1;
    int want = VesPool_share();
    if (want > max_helpers) want = max_helpers;
    size_t chunks = (n + chunk - 1) / chunk;
    if ((size_t)want > chunks - 1) want = (int)(chunks - 1);
    if (want <= 0) { fn(ctx, 0, n, 0); return 1; }
    VpRun r;
    memset(&r, 0, sizeof r);
    r.fn = fn; r.ctx = ctx; r.n = n; r.chunk = chunk;
    vp_lock(&vp_mutex);
    if (!vp_started) {
        int made = 0;
        for (int t = 0; t < vp_cores - 1; t++) if (!vp_spawn()) made++;
        vp_idle += made; vp_started = 1;
    }
    if (want > vp_idle) want = vp_idle;
    if (want <= 0) { vp_unlock(&vp_mutex); fn(ctx, 0, n, 0); return 1; }
    vp_idle -= want;
    r.tickets = want;
    if (vp_tail) vp_tail->link = &r; else vp_head = &r;
    vp_tail = &r;
    vp_runs++; vp_borrowed += (size_t)want;
    if (want == 1) vp_wake_one(&vp_work); else vp_wake_all(&vp_work);
    vp_unlock(&vp_mutex);
    vp_chunks(&r, 0);
    vp_lock(&vp_mutex);
    if (r.tickets > 0) {                      /* seats nobody took: give them back */
        VpRun **at = &vp_head, *prev = NULL;
        while (*at && *at != &r) { prev = *at; at = &(*at)->link; }
        if (*at == &r) { *at = r.link; if (vp_tail == &r) vp_tail = prev; }
        vp_idle += r.tickets; r.tickets = 0;
    }
    while (r.running > 0) vp_wait(&vp_done, &vp_mutex);
    int slots = 1 + r.joined;
    vp_unlock(&vp_mutex);
    return slots;
}

void VesPool_counts(size_t *runs, size_t *helpers)
{
    if (runs) *runs = vp_runs;
    if (helpers) *helpers = vp_borrowed;
}

/* ---------------------------------------------------------------- selftest */
typedef struct VpTest { long long sum[64]; size_t calls[64]; } VpTest;

static void vp_test_fn(void *ctx, size_t b, size_t e, int slot)
{
    VpTest *t = (VpTest *)ctx;
    for (size_t i = b; i < e; i++) t->sum[slot] += (long long)i;
    t->calls[slot]++;
}

int VesPool_selftest(void)
{
    int failures = 0, saved = vp_cores;
    VesPool_configure(8);
    /* (1) one outer job: up to 7 helpers; the slots' sums add to the total */
    VesPool_job_enter();
    VpTest t; memset(&t, 0, sizeof t);
    size_t n = 1000003;
    int slots = VesPool_run(n, 97, vp_test_fn, &t, 64);
    long long total = 0; for (int s = 0; s < 64; s++) total += t.sum[s];
    for (int s = slots; s < 64; s++) if (t.calls[s]) failures++;
    if (total != (long long)n * (long long)(n - 1) / 2 || slots < 1 || slots > 8) failures++;
    /* (2) eight outer jobs: no core idle, so no helper is borrowed */
    for (int k = 0; k < 7; k++) VesPool_job_enter();
    if (VesPool_share() != 0) failures++;
    memset(&t, 0, sizeof t);
    if (VesPool_run(n, 97, vp_test_fn, &t, 64) != 1 || t.sum[0] != (long long)n * (long long)(n - 1) / 2) failures++;
    for (int k = 0; k < 8; k++) VesPool_job_leave();
    /* (3) a range smaller than one chunk runs on the caller alone */
    memset(&t, 0, sizeof t);
    if (VesPool_run(5, 97, vp_test_fn, &t, 64) != 1 || t.sum[0] != 10) failures++;
    /* (4) repeated runs leave every helper idle again */
    for (int k = 0; k < 200; k++) { memset(&t, 0, sizeof t); VesPool_run(10000, 7, vp_test_fn, &t, 64); }
    vp_lock(&vp_mutex); if (vp_started && vp_idle != vp_cores - 1) failures++; vp_unlock(&vp_mutex);
    VesPool_configure(saved);
    fprintf(stderr, "  ves_pool lending (sums, no idle core, small range, helpers returned): %s (%d failures)\n",
            failures ? "FAIL" : "ok", failures);
    return failures;
}
