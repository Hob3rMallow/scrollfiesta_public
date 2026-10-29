#ifndef SLIM_REFINE_INCLUDED
#define SLIM_REFINE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Global symmetric-Dirichlet refinement for an already-injective collection
 * of manifold disk charts.  The input embedding is a feasibility certificate;
 * this stage is the metric optimizer.  Every accepted iterate preserves face
 * orientation and a simple boundary for every chart. */
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
    int local_steps;          /* locally damp flip-limiting large charts */
    size_t local_step_min_faces; /* minimum chart size for local damping */
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
    double energy_before;
    double energy_after;
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
