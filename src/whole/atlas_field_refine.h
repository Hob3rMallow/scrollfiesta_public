#ifndef ATLAS_FIELD_REFINE_INCLUDED
#define ATLAS_FIELD_REFINE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_field.h"

typedef struct {
    size_t mesh_components;
    size_t quotient_components;
    size_t auxiliary_variables;
    size_t constraints;
    size_t anchors;
    size_t skipped_degenerate_triangles;

    size_t seam_rows;
    size_t order_rows;
    double seam_rms_before;
    double seam_rms_after;
    double seam_max_after;
    double order_rms_before;
    double order_rms_after;
    double order_max_after;

    double correction_rms;
    double correction_max;
    double edge_gradient_delta_rms;
    double edge_gradient_delta_max;
    MonotoneQpStats qp;
} AtlasFieldRefineStats;

/*
 * Solve one sparse FEM correction problem on the disjoint retained mesh.
 * Constraint variables are [0,nvertices) for mesh U values followed by
 * nauxiliary latent values. Every constraint row must be translation
 * invariant (its coefficients sum to zero); this lets the helper infer one
 * exact gauge per connected quotient component. The largest-area local mesh
 * chart in each quotient supplies that gauge at its existing base value.
 */
int AtlasFieldRefine_solve(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const int32_t *triangles, size_t ntriangles,
    const int32_t *vertex_component, size_t ncomponents,
    const double *base_u,
    const AtlasFieldConstraint *constraint, size_t nconstraint,
    const AtlasFieldConstraintCoeff *constraint_coeff,
    size_t nconstraint_coeff,
    size_t nauxiliary, const double *auxiliary_initial,
    const MonotoneQpOptions *qp_options,
    double *out_u, double *out_auxiliary,
    AtlasFieldRefineStats *stats);

#endif
