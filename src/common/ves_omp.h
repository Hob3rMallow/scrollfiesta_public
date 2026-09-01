#ifndef VES_OMP_H
#define VES_OMP_H

/* ================================================================
 * OpenMP runtime entry points, declared rather than included.
 *
 * No translation unit under src/ may include <omp.h> -- CI enforces
 * it with a grep. The reason is concrete: an embedding consumer can
 * compile this code against its own omp.h shim on the include path,
 * and at least one of those shims is C++-only, so pulling the real
 * header turns into a "cannot find cstddef" failure in a C build.
 *
 * `#pragma omp` directives need no declaration at all; only the
 * handful of runtime calls below do. Their signatures are fixed by
 * the OpenMP specification.
 *
 * Without _OPENMP the calls resolve to the serial answers, so call
 * sites need no #ifdef of their own.
 * ================================================================ */

#ifdef _OPENMP

extern int  omp_get_max_threads(void);
extern int  omp_get_thread_num(void);
extern int  omp_get_num_procs(void);
extern int  omp_in_parallel(void);
extern void omp_set_num_threads(int);
extern void omp_set_dynamic(int);

#else

static inline int  omp_get_max_threads(void) { return 1; }
static inline int  omp_get_thread_num(void)  { return 0; }
static inline int  omp_get_num_procs(void)   { return 1; }
static inline int  omp_in_parallel(void)     { return 0; }
static inline void omp_set_num_threads(int n) { (void)n; }
static inline void omp_set_dynamic(int n)     { (void)n; }

#endif /* _OPENMP */

#endif /* VES_OMP_H */
