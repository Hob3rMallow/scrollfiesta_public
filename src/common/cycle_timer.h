#ifndef CYCLE_TIMER_INCLUDED
#define CYCLE_TIMER_INCLUDED

/*
 * cycle_timer -- a cheap monotonic tick counter for phase instrumentation.
 *
 * cycle_now() reads the invariant TSC on x86/x64 (a few ns, no syscall), so it
 * can bracket per-face / per-site work in hot loops without perturbing them.
 * Ticks are NOT seconds: callers convert with a ratio measured over a span
 * they also time with ves_clock_sec() (e.g. a whole pipeline stage).  On other
 * targets it falls back to nanoseconds from the monotonic clock.
 *
 * Diagnostics only: nothing read from here may feed a geometric result.
 */

#include <stdint.h>

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <intrin.h>
static inline uint64_t cycle_now(void) { return (uint64_t)__rdtsc(); }
#elif defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
static inline uint64_t cycle_now(void) { return (uint64_t)__rdtsc(); }
#else
#include "ves_platform.h"
static inline uint64_t cycle_now(void)
{
    return (uint64_t)(ves_clock_sec() * 1.0e9);
}
#endif

#endif /* CYCLE_TIMER_INCLUDED */
