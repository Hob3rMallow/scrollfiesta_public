/* Hand-written TAUCS build config (replaces the configurator output).
 * vesuvius-c flatten module: double-real only, no complex, no METIS.
 * See deps/taucs/configurator/taucs_structure.h for the original mechanism. */
#ifndef TAUCS_CONFIG_BUILD_H
#define TAUCS_CONFIG_BUILD_H

#define TAUCS_CONFIG_DREAL      /* enable double-precision real datatype */
/* (SREAL / DCOMPLEX / SCOMPLEX intentionally disabled) */

/* Ordering: genmmd only (minimum-degree). METIS intentionally OFF (the
 * Feb-2025 fork README notes Windows builds without METIS); AMD/COLAMD left
 * off to minimize the vendored set. genmmd is the default symmetric ordering
 * (taucs_linsolve.c) and external/src/genmmd.c is compiled in. Add
 * TAUCS_CONFIG_AMD or TAUCS_CONFIG_COLAMD (and the matching external .c files)
 * only if the linker reports missing ordering symbols. */
#define TAUCS_CONFIG_GENMMD

#endif /* TAUCS_CONFIG_BUILD_H */
