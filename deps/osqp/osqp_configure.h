#ifndef OSQP_CONFIGURE_H
#define OSQP_CONFIGURE_H
#define OSQP_ALGEBRA_BUILTIN
/* Match the reference wheel: double values and 32-bit sparse indices.
 * Dynamic setup and polishing are enabled; no Python or code generation. */
#ifdef _WIN32
#define IS_WINDOWS
#else
#define IS_LINUX
#endif
#endif
