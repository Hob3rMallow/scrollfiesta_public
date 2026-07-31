#ifndef ATLAS_FIELD_INCLUDED
#define ATLAS_FIELD_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "monotone_qp.h"

/*
 * Piecewise-linear FEM extension of a resolved coarse atlas over a disjoint
 * union of retained triangle patches.  The local term preserves the raw
 * per-triangle u gradient; sparse barycentric observations determine the
 * otherwise-free chart gauges.
 */

typedef struct {
    int32_t vertex[3];
    double bary[3];
    double target;
    double weight;
    int32_t source;
} AtlasFieldObservation;

/*
 * Optional general sparse rows.  `first` indexes problem.constraint_coeff,
 * not the assembled MonotoneQp coefficient array.  These rows are used for
 * relations such as
 *
 *     U[seam_b] - U[seam_a] = 0
 *
 * and for support-averaged winding/order observations.  Keeping them in the
 * same least-squares assembly is the discrete counterpart of adding
 * A_s^T W_s A_s to the FEM stiffness matrix.
 */
typedef struct {
    size_t first;
    int32_t count;
    double target;
    double weight;
    int32_t kind;
    int32_t source;
} AtlasFieldConstraint;

typedef struct {
    int32_t variable;
    double coefficient;
} AtlasFieldConstraintCoeff;

typedef struct {
    const double *position;       /* 3*nvertex */
    const double *u0;             /* nvertex raw/local atlas */
    size_t nvertex;
    /* Optional latent support variables used to keep regional observations
     * star-sparse. They have no triangle-gradient rows of their own. */
    size_t nauxiliary;

    const int32_t *triangle;      /* 3*ntriangle */
    const double *triangle_weight;/* nullable; defaults to 1 */
    size_t ntriangle;

    const AtlasFieldObservation *observation;
    size_t nobservation;

    const AtlasFieldConstraint *constraint;
    size_t nconstraint;
    const AtlasFieldConstraintCoeff *constraint_coeff;
    size_t nconstraint_coeff;

    const MonotoneQpAnchor *anchor;
    size_t nanchor;
} AtlasFieldProblem;

enum {
    ATLAS_FIELD_ROW_GRADIENT = 20,
    ATLAS_FIELD_ROW_OBSERVATION = 21,
    ATLAS_FIELD_ROW_SEAM = 22,
    ATLAS_FIELD_ROW_ORDER = 23,
    ATLAS_FIELD_ROW_SUPPORT = 24
};

typedef struct {
    MonotoneQpProblem qp;
    size_t gradient_rows;
    size_t observation_rows;
    size_t constraint_rows;
    size_t skipped_degenerate_triangles;
} AtlasFieldSystem;

/*
 * Assemble
 *
 *   1/2 sum_t A_t w_t ||grad U - grad U0||^2
 * + 1/2 sum_o mu_o (h_o^T U - u_o)^2
 * + 1/2 sum_c w_c (a_c^T [U;q] - b_c)^2.
 *
 * Returned rows live in arena.  There are no inequality bounds at this level;
 * exact gauges, if any, are eliminated by MonotoneQp_solve.
 */
int AtlasField_build(Arena_T arena,
                     const AtlasFieldProblem *problem,
                     AtlasFieldSystem *system);

#endif
