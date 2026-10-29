#ifndef WINDING_FIELD_INCLUDED
#define WINDING_FIELD_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Generalized-winding-number field over an oriented triangle mesh/soup.
 *
 * DIRECT is the Van Oosterom--Strackee solid-angle sum and is the reference.
 * BOUNDARY_EXACT is Xie--Hafner--Wojtan (2026): a signed ray count minus the
 * solid angles of a generalized cone over the oriented exterior edges.
 * FAST is the second-order Barill et al. (2018) expansion, with exact triangle
 * evaluation in the near field.  AUTO selects the exact boundary method when
 * its measured boundary work is affordable, otherwise FAST. */
typedef struct WindingField_T *WindingField_T;

typedef enum {
    WINDING_FIELD_AUTO = 0,
    WINDING_FIELD_DIRECT = 1,
    WINDING_FIELD_BOUNDARY_EXACT = 2,
    WINDING_FIELD_FAST = 3
} WindingFieldBackend;

typedef struct {
    WindingFieldBackend backend;
    size_t leaf_faces;          /* Morton-BVH leaf capacity (default 32). */
    double beta;                /* far/near separation (default 2). */
    uint64_t exact_work_limit;  /* max boundary-edge evaluations in AUTO. */
} WindingFieldOptions;

typedef struct {
    size_t vertices, faces;
    size_t bvh_nodes, bvh_leaves;
    size_t exterior_edges;
    size_t exterior_multiplicity;
} WindingFieldStats;

typedef struct {
    WindingFieldBackend backend_used;
    size_t queries;
    size_t direct_queries;
    size_t boundary_queries;
    size_t fast_queries;
    size_t degenerate_ray_fallbacks;
    size_t exact_triangle_evaluations;
    size_t approximate_node_evaluations;
    size_t invalid_normals;
} WindingFieldEvalStats;

void WindingField_default_options(WindingFieldOptions *options);

/* vertices are float32 triples in the repository's (z,y,x) order.  faces are
 * oriented triangles.  Their storage must outlive the returned arena-owned
 * handle. */
int WindingField_build(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const int32_t *faces, size_t nfaces,
    const WindingFieldOptions *options,
    WindingField_T *out_field,
    WindingFieldStats *stats);

/* Evaluate double-precision query triples. */
int WindingField_evaluate_many(
    WindingField_T field,
    const double *queries, size_t nqueries,
    WindingFieldBackend backend,
    double *out_winding,
    WindingFieldEvalStats *stats);

/* Evaluate q+ = p + eps*n and q- = p - eps*n.  out_mean is
 * 0.5*(w+ + w-); out_jump is w- - w+.  A unit jump is the local clean-sheet
 * law.  Invalid/zero normals produce NAN outputs and are counted. */
int WindingField_evaluate_sides(
    WindingField_T field,
    const float *points, const float *normals, size_t npoints,
    double eps, WindingFieldBackend backend,
    float *out_mean, float *out_jump,
    WindingFieldEvalStats *stats);

/* Synthetic exact/fast/boundary-cone regression tests. */
int WindingField_selftest(void);

#endif
