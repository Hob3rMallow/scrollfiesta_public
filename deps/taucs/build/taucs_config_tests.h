/* Hand-written TAUCS test config (replaces the configurator's link probes).
 * vesuvius-c flatten module. */
#ifndef TAUCS_CONFIG_TESTS_H
#define TAUCS_CONFIG_TESTS_H

/* BLAS/LAPACK symbols carry a trailing underscore (CLAPACK *_nowrap libs:
 * dgemm_, dsyrk_, dtrsm_, dpotrf_, ...). */
#define TAUCS_BLAS_UNDERSCORE

/* Platform select (the configurator would emit OSTYPE_<plat> via -D). */
#if defined(_WIN32)
#  ifndef OSTYPE_win32
#    define OSTYPE_win32
#  endif
#else
#  ifndef OSTYPE_linux
#    define OSTYPE_linux
#  endif
#endif

/* No C99 complex: we only use double-real, so taucs_complex.c takes the
 * struct-based path. (TAUCS_C99_COMPLEX intentionally left undefined.) */

#endif /* TAUCS_CONFIG_TESTS_H */
