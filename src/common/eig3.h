#ifndef EIG3_INCLUDED
#define EIG3_INCLUDED

/* ============================================================================
 * eig3.h -- shared symmetric 3x3 eigensolver (cyclic Jacobi, double precision).
 *
 * The core is the same routine that backs pca.c, lifted out so callers can read
 * the EIGENVECTOR of the smallest eigenvalue -- which pca.c hides behind its
 * PCA-normal sign convention. The developability energy (Stein/Grinspun/Crane)
 * needs lambda_min and its eigenvector of a per-vertex normal-covariance matrix,
 * so it uses this directly. All inputs symmetric; only the upper/lower triangle
 * matters (the routine treats a[][] as symmetric).
 * ==========================================================================*/

/* Full decomposition of a symmetric 3x3 matrix `a` (not modified).
 *   evals: eigenvalues in ASCENDING order.
 *   evecs: eigenvectors as COLUMNS (evecs[row][col]); evecs[*][k] is the unit
 *          eigenvector for evals[k]. */
void Eig3_sym(const double a[3][3], double evals[3], double evecs[3][3]);

/* Smallest eigenpair only.
 *   *out_eval = evals[0].
 *   out_vec   = its unit eigenvector with a DETERMINISTIC sign (the first
 *               component with magnitude > 1e-12 is made positive), so finite-
 *               difference gradient checks against it are reproducible. */
void Eig3_smallest(const double a[3][3], double *out_eval, double out_vec[3]);

/* In-process unit test. Returns 0 if all checks pass, else the failure count. */
int Eig3_selftest(void);

#endif
