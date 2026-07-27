/* C++-only poison omp.h, mirroring VC3D's core/openmp_stub/omp.h: any C
 * translation unit that includes <omp.h> by name fails to compile here.
 * ScrollFiesta library code must never include the real header (see
 * ves_platform.c for the self-declared-prototype pattern). */
#include <cstddef>

static inline int  omp_get_max_threads(void) { return 1; }
static inline int  omp_get_thread_num(void) { return 0; }
static inline int  omp_get_num_threads(void) { return 1; }
static inline void omp_set_num_threads(int n) { (void)n; }
