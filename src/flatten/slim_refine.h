#ifndef SLIM_REFINE_INCLUDED
#define SLIM_REFINE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Global symmetric-Dirichlet refinement for an already-injective collection
 * of manifold disk charts.  The input embedding is a feasibility certificate;
 * this stage is the metric optimizer.  Every accepted iterate preserves face
 * orientation and a simple boundary for every chart. */
/* Physical material arms in the mesh's 3-D coordinate order. The solver
 * projects them into each face tangent plane and normalizes their lengths.
 * weight is relative evidence; zero leaves that face unconstrained. Evidence
 * provenance belongs to the caller and is never inferred from input UV. */
typedef struct SlimRefineDirection {
    double h[3], v[3], weight;
} SlimRefineDirection;

typedef struct SlimRefineOpts {
    int iterations;          /* local/global iterations (default 20) */
    int homotopy_stages;     /* rest-metric continuation stages (default 1) */
    int final_iterations;    /* extra iterations at the final metric */
    int pcg_iterations;      /* global proxy PCG cap (default 120) */
    double pcg_tolerance;    /* residual relative to initial correction (1e-5) */
    double energy_tolerance; /* relative accepted-energy stop (default 1e-6) */
    double sigma_min;        /* local-step singular-value floor (default 1e-8) */
    double weight_max;       /* local-step weight ceiling (default 1e8) */
    double proximal;         /* fraction of per-component mean diagonal (1e-5) */
    int line_search;         /* flip-avoiding halvings (default 24) */
    int use_amg;             /* block aggregation multigrid PCG preconditioner */
    /* Optional symmetric block-system backend. diagonal[nv*3] and blocks[ne*3]
     * store 00,01,11; rhs/result[nv*2] are interleaved. result starts at the
     * current UV and receives the proposal; return 0 on success. The regular
     * energy, face and boundary guards still apply. NULL uses AMG/PCG. */
    void *linear_context;
    int (*linear_solve)(Arena_T, void *context, size_t nv, size_t ne, const double *diagonal,
                        const int32_t *edge_a, const int32_t *edge_b,
                        const double *blocks, const double *rhs, double *result);
    int local_steps;          /* locally damp flip-limiting large charts */
    size_t local_step_min_faces; /* minimum chart size for local damping */
    int local_line_search;    /* local feasibility sweeps, independent of uniform halvings */
    int local_growth_sweeps;  /* feasible alpha growth after local damping fails; default 4, 0 disables */
    int local_relaxation;     /* coordinate proxy sweeps after feasible growth or a stalled global step */
    int local_relax_boundary; /* also relax the rim, with swept-edge collision checks */
    const SlimRefineDirection *directions; /* [nf], optional physical material arms */
    double direction_weight; /* area-weighted joint cross-axis objective; default 0 */
    int direction_rigid_gauge; /* optional one common original-frame rigid turn; default 0 */
    const float *guide_uv;   /* optional winding/developable UV, [nv*2] */
    const float *metric_uv;  /* optional locally valid winding UV whose
                              * per-face metric replaces noisy 3-D rest metric */
    int guide_warm_start;    /* guarded metric-improving morph toward guide */
    int guard_boundary;      /* reject self-intersecting chart boundaries */
    int strip_pack;          /* repack components largest-first in one row */
    double padding;          /* component bounding-box gutter (default 20) */
    int threads;             /* OpenMP threads; <=0 uses runtime default */
    int verbose;
    /* Optional solver-only scaffold. Faces before scaffold_first_face keep
     * their original 3-D rest metric. Remaining faces use scaffold_uv [nv*2]
     * as their rest shape, with equal energy weight scaffold_weight per face
     * (SCAF Eq. 3 cancels scaffold area). No source face may be
     * relabelled as scaffold. NULL disables this extension. */
    const double *scaffold_uv;
    size_t scaffold_first_face;
    double scaffold_weight;
    const uint8_t *fixed_vertices; /* [nv], hard pin to the input UV */
} SlimRefineOpts;

void SlimRefine_defaults(SlimRefineOpts *opts);

/* Open coincident UV seam banks along caller-supplied inward displacements.
 * Backtracking uses the same face and component-boundary certificates as SLIM.
 * This is seam surgery only; it does not replace metric optimization. World
 * geometry is untouched, and failure leaves uv_out unchanged. */
int SlimRefine_open_seam_seed(Arena_T arena, const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const double *uv, const double *displacement,
                              double *uv_out, double *accepted_alpha);

/* Symmetric-Dirichlet energy above its minimum, its exact Jacobian gradient,
 * and the SLIM local quadratic (2 W^T W, 00/01/11). The same proxy is used by
 * constrained seam/contact repair; its outer acceptance certificate remains
 * independent. No caller may substitute an ARAP recovery on failure. */
int SlimRefine_face_proxy(const double jacobian[4], const SlimRefineOpts *opts,
                          double *energy, double gradient[4], double block[3]);

typedef struct SlimRefineStats {
    size_t components;
    size_t edges;
    size_t boundary_edges;
    size_t metric_uv_faces;
    size_t metric_fallback_faces;
    size_t iterations_run;
    size_t homotopy_stages_run;
    size_t accepted_steps;
    size_t component_steps_accepted;
    size_t component_steps_rejected;
    size_t line_search_trials;
    double accepted_alpha_min;
    double accepted_alpha_mean;
    double accepted_alpha_max;
    size_t local_component_steps_accepted;
    size_t growth_component_steps_accepted;
    size_t growth_vertex_steps;
    double growth_vertex_alpha_mean;
    double growth_vertex_alpha_max;
    size_t relaxation_vertex_steps;
    double local_vertex_alpha_min;
    double local_vertex_alpha_mean;
    double local_vertex_alpha_max;
    size_t guide_components_accepted;
    size_t guide_components_rejected;
    size_t guide_line_search_trials;
    double guide_alpha_min;
    double guide_alpha_mean;
    double guide_alpha_max;
    double guide_energy_after;
    size_t pcg_iterations_total;
    int pcg_iterations_max;
    double pcg_relative_residual_max;
    int amg_levels;
    size_t amg_coarse_vertices;
    double energy_before;    /* pure symmetric Dirichlet; final value may not increase */
    double energy_after;
    double objective_before; /* SD plus optional direction objective */
    double objective_after;
    size_t direction_faces;
    double direction_area_fraction; /* physical area with positive evidence / all area */
    double direction_energy_before;
    double direction_energy_after;
    double direction_angle_before; /* evidence-area weighted radians^2, both arms */
    double direction_angle_after;
    double direction_gauge_radians;
    int direction_gauge_accepted;
    double qc_mean_before;
    double qc_mean_after;
    double qc_max_before;
    double qc_max_after;
    double min_det_before;
    double min_det_after;
    size_t flips_before;
    size_t flips_after;
    size_t boundary_intersections_before;
    size_t boundary_intersections_after;
    double atlas_width;
    double atlas_height;
    double setup_seconds;
    double guide_seconds;
    double assembly_seconds;
    double amg_seconds;
    double pcg_seconds;
    double line_search_seconds;
    double finalize_seconds;
    double total_seconds;
    int converged;
    int reverted;
} SlimRefineStats;

int SlimRefine_run_double(Arena_T arena,
                          const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const double *uv_in,
                          const SlimRefineOpts *opts,
                          double *uv_out,
                          SlimRefineStats *stats);

int SlimRefine_selftest(void);

#endif
