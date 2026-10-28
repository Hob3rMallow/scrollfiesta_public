#ifndef VES_POOL_H
#define VES_POOL_H
/* A process-wide pool of helper threads that lends IDLE cores to a parallel
 * loop running inside another parallel region (where OpenMP 2.0 would
 * serialize a nested team).  Outer jobs announce themselves with
 * VesPool_job_enter/leave; a loop may borrow its fair share of the cores the
 * outer jobs leave idle: (cores - jobs) / jobs helpers, so the threads running
 * never exceed `cores`.  A helper runs chunks of the caller's range; the
 * caller runs chunks too, as slot 0, and returns when every chunk is done.
 * Which thread runs which chunk is not deterministic -- callers must make
 * their result independent of it (per-slot buffers merged in a fixed order). */
#include <stddef.h>

typedef void (*VesPoolFn)(void *ctx, size_t begin, size_t end, int slot);

/* cores <= 1 disables lending (every run is serial on the caller).  Helper
 * threads (cores - 1) are created on the first run that borrows. */
void VesPool_configure(int cores);
int  VesPool_cores(void);
void VesPool_job_enter(void);
void VesPool_job_leave(void);
/* Helpers a run started now could borrow (0 when no core is idle). */
int  VesPool_share(void);
/* Run fn over [0, n) in chunks of `chunk` on the caller plus up to
 * max_helpers borrowed helpers (bounded by the share).  Returns the number of
 * slots that may have been used (1 + helpers borrowed); slots are < that. */
int  VesPool_run(size_t n, size_t chunk, VesPoolFn fn, void *ctx, int max_helpers);
/* Counters for diagnostics: runs that borrowed, helpers borrowed in total. */
void VesPool_counts(size_t *runs, size_t *helpers);
int  VesPool_selftest(void);
#endif
