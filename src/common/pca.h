#ifndef PCA_INCLUDED
#define PCA_INCLUDED

#include <stddef.h>

/* Compute PCA normal of N points.
 * out_normal: unit eigenvector of smallest eigenvalue.
 * out_centroid: mean position (may be NULL).
 * Sign convention: largest-magnitude component positive.
 * Uses double internally for covariance accumulation.
 * Returns 0 on success, -1 if N < 3 (degenerate). */
int PCA_normal(const float *points, size_t n,
               float out_normal[3], float out_centroid[3]);

/* Compute PCA principal axis of N points.
 * out_axis: unit eigenvector of the LARGEST eigenvalue (direction of greatest
 * spread) -- e.g. the long/winding axis of a scroll point cloud.
 * out_centroid: mean position (may be NULL).
 * Sign convention: oriented toward (1,1,1), matching PCA_normal.
 * Component order matches the input points (the pipeline's (z,y,x)).
 * Returns 0 on success, -1 if N < 3 (degenerate -> out_axis = (1,0,0)). */
int PCA_principal_axis(const float *points, size_t n,
                       float out_axis[3], float out_centroid[3]);

/* Project N points onto a direction vector.
 * out_proj[i] = dot(point[i], dir).
 * Returns (min_proj, max_proj) via out pointers. */
void PCA_project(const float *points, size_t n, const float dir[3],
                 float *out_proj, float *out_min, float *out_max);

/* Build orthonormal UV basis perpendicular to a given normal.
 * Uses Hughes-Moller technique for robustness. */
void PCA_orthonormal_basis(const float normal[3],
                           float out_u[3], float out_v[3]);

#endif
