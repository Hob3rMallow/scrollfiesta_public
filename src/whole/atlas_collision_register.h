#ifndef ATLAS_COLLISION_REGISTER_INCLUDED
#define ATLAS_COLLISION_REGISTER_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_register.h"

/*
 * Collision-aware rigid U registration with a fixed chart order.
 *
 * The metric rows and desired shifts are soft.  Exact positive-area triangle
 * intersections not explained by spatially local duplicate coverage propose
 * one-sided bounds in the supplied post-discrete order.  Bounds compatible
 * with the topology budget are hard; incompatible intersections remain in a
 * classified topology-limited ledger.  Redundant local observations are kept
 * for later ownership.  V and every chart's intrinsic UV shape are immutable.
 */
typedef struct {
    const int32_t *faces;
    size_t nfaces;
    size_t nvertices;
    const uint8_t *face_keep;       /* NULL keeps every face */
    const double *base_u;
    const double *v;
    const int32_t *vertex_chart;
    size_t ncharts;

    /* Optional evidence for classifying redundant same-sheet coverage. */
    const float *xyz;                 /* [3*nvertices], source coordinates */
    const uint8_t *allowed_chart_pair;/* [ncharts*ncharts], direct topology */
    const double *chart_radius;       /* [ncharts], optional */
    double same_sheet_radius_tolerance;
    double same_sheet_xyz_tolerance;

    /* chart_rank is a permutation inverse: chart_rank[c] is c's fixed
     * left-to-right rank in [0,ncharts). */
    const size_t *chart_rank;

    /* The unconstrained metric solution and its per-chart confidence. */
    const double *desired_shift;
    const double *desired_weight;   /* NULL means one */

    /* Difference rows shift[component1] - shift[component0] = target. */
    const AtlasRegisterEdge *edges;
    size_t nedges;
    /* Hard preservation budget around the unconstrained metric difference on
     * each non-radial topology edge.  Zero disables topology bounds. */
    double topology_deviation_limit;

    double collision_margin;        /* U clearance beyond exact contact */
    int max_outer_iterations;
} AtlasCollisionRegisterProblem;

typedef struct {
    int32_t chart_lo;
    int32_t chart_hi;
    int32_t face_lo;              /* source face which last tightened the bound */
    int32_t face_hi;
    double lower;                 /* shift[hi] - shift[lo] >= lower */
    double final_delta;
    double slack;
    uint32_t updates;
    int topology;
    int collision;
    int initial;                  /* present in the post-discrete audit */
    int active;
} AtlasCollisionRegisterBound;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    int32_t face0;                /* source face ids */
    int32_t face1;
    int allowed;
    int reason;       /* 0 hard, 1 direct, 2 XYZ-near, 3 self, 4 topology */
    double min_xyz;   /* closest source-space vertex pair on these triangles */
    double radius_delta;
} AtlasCollisionRegisterResidual;

typedef struct {
    size_t kept_faces;
    size_t metric_edges;
    size_t rejected_edges;
    size_t total_bounds;
    size_t topology_bounds;
    size_t collision_bounds_rejected;
    size_t collision_bounds;
    size_t collision_bounds_added;
    size_t outer_iterations;
    size_t exact_pairs_before;
    size_t exact_pairs_after;
    size_t exact_cross_before;
    size_t exact_cross_after;
    size_t moved_charts;
    int qp_iterations;
    int qp_active_bounds;
    double qp_objective;
    double qp_stationarity;
    size_t exact_allowed_before;
    size_t exact_allowed_after;
    size_t exact_hard_before;
    size_t exact_hard_after;
    size_t exact_topology_limited_before;
    size_t exact_topology_limited_after;
    double qp_max_violation;
    double shift_rms;
    double shift_max;
    int solve_failed;
    int collision_failed;
    int topology_infeasible;
    const AtlasCollisionRegisterBound *bound;
    size_t nbound;
    const AtlasCollisionRegisterResidual *residual;
    size_t nresidual;
} AtlasCollisionRegisterStats;

int AtlasCollisionRegister_solve(
    Arena_T arena, const AtlasCollisionRegisterProblem *problem,
    double *out_shift, AtlasCollisionRegisterStats *stats);

int AtlasCollisionRegister_selftest(void);

#endif
