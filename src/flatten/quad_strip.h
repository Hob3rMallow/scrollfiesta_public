#ifndef QUAD_STRIP_INCLUDED
#define QUAD_STRIP_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* quad_strip -- build ONE quad strip for a single scroll component.
 *
 * Rows are v = world z and columns are u = winding arc length.  Production
 * preserves the component rectangle so missing fitted samples can be lofted
 * and corrected without silently discarding large parts of the sheet.  The
 * BOUNDED mode remains available when a deliberately sparse, alpha-backed
 * export is wanted.  Components are never bridged here.
 *
 * The z/y/x split is only internal to the PDE (SnapCG solves the three position
 * channels separately); the strip itself is one grid.
 */

typedef enum {
    QUAD_STRIP_RECT = 0,    /* production: fill the full bbox rectangle */
    QUAD_STRIP_SPAN = 1,    /* diagnostic: per-row + per-col span hull */
    QUAD_STRIP_BOUNDED = 2  /* conservative: certified domain + closed holes */
} QuadStripMode;

typedef struct QuadStripOpts {
    QuadStripMode mode;          /* default RECT */
    int           max_hole_distance; /* largest 4-neighbour distance from a
                                        closed hole cell to its fitted rim;
                                        wider holes remain alpha (default 16) */
    int           max_row_gap;   /* nearest live lattice row which may form a
                                    quad with this row (default 1).  Ragged
                                    fields commonly need 2 for alternating-z
                                    direct support. */
    double        fit_stiffness; /* anchor weight on fitted cells (default 1e4) */
    int           exact_fitted;  /* re-stamp fitted cells with the exact sample
                                    after the solve (default 1) */
    int           pde_max_iter;  /* SnapCG iterations (default 1024) */
    double        pde_tol;       /* SnapCG relative tolerance (default 1e-6) */
    int           cylindrical_fill; /* fit an affine lattice trend in
                                       (z,r,unwrapped theta), harmonically extend
                                       only its residual across holes, then map
                                       back to Cartesian (default 0) */
    double        axis_y, axis_x;/* cylindrical scroll axis point in y/x */
    int           relax_rounds;  /* ARAP developable-fill relaxation rounds after
                                    the harmonic fill (default 0 = off). Each round
                                    fits a local rotation per cell to the flat UV
                                    rest (isometric target) and re-solves positions
                                    with the fitted rim fixed, so the fill unfolds
                                    developably instead of cutting through wraps. */
    int           relax_cg_iters;/* CG iterations per global ARAP solve (default 96) */
} QuadStripOpts;

/* Full-ribbon metric ARAP.  This is deliberately separate from the historical
 * fill-only relax_rounds option above: every vertex is movable, source samples
 * are soft robust anchors, and the flat UV lattice is the intrinsic rest
 * metric.  One local frame fit plus three global Poisson solves is one outer
 * iteration. */
typedef struct QuadStripArapOpts {
    int    max_iterations;       /* outer local/global iterations (default 1000) */
    double movement_tolerance;   /* stop when max XYZ move <= this (default 1e-5) */
    double movement_tolerance_rel; /* effective tolerance = max(movement_tolerance,
                                      this * first-iteration max move).  Ties
                                      "done" to the solve's own motion scale: a
                                      full sheet starting at ~1 vox stops near
                                      1e-3 vox (sub-bake-resolution) instead of
                                      burning 1000 iterations against the absolute
                                      1e-5 (default 1e-3; 0 = legacy absolute) */
    int    movement_stall_window;  /* stop as stalled (converged=0) after this many
                                      consecutive iterations that fail to improve
                                      the running-min max move by
                                      movement_stall_fraction (default 50; 0=off) */
    double movement_stall_fraction;/* relative improvement floor (default 0.01) */
    int    cg_max_iterations;    /* iterations per scalar Poisson solve (default 256) */
    double cg_relative_tolerance;/* relative residual target (default 1e-8) */
    int    multigrid_levels;     /* aggregate hierarchy depth incl. fine (default 16) */
    int    multigrid_cycles;     /* V-cycles before PCG fallback (default 1) */
    int    multigrid_pre_sweeps; /* fine/coarse forward smoothing (default 2) */
    int    multigrid_post_sweeps;/* reverse smoothing (default 2) */
    int    multigrid_coarse_sweeps; /* coarsest symmetric sweeps (default 32) */
    int    nonlinear_multigrid_levels; /* coarse-to-fine ARAP levels incl. fine
                                          (default 4; activates above threshold) */
    int    nonlinear_coarse_iterations; /* outer ARAP rounds per coarse level
                                           (default 128) */
    int    nonlinear_min_vertices; /* minimum fine vertices before nonlinear
                                      hierarchy activates (default 4000000) */
    double source_weight;        /* measured-cell soft anchor (default 1) */
    double fill_weight;          /* introduced-cell initializer anchor (default .01) */
    double huber_delta;          /* source-anchor displacement Huber scale (default 2) */
    double minimum_source_scale; /* minimum robust source multiplier (default .05) */
    double minimum_area_ratio;   /* per-step oriented-area barrier (default .05) */
    double maximum_vertex_step;  /* local ARAP trust radius in XYZ voxels
                                    (default 2) */
    int    preserve_axial;       /* keep world z exactly at the anchor value.
                                    Cylindrical post-snap relaxation uses this
                                    because lattice v already fixes z (default 0). */
    const uint8_t *fixed_vertices; /* optional H*W exact Dirichlet mask. Fixed
                                      vertices stay at their anchor XYZ while
                                      their neighbours relax (default NULL). */
    const float *repair_weight;  /* optional H*W marble-repair confidence;
                                    0 preserves ordinary ARAP, 1 weakens the
                                    observed position and follows a fair UV
                                    tangent frame (default NULL) */
    double repair_source_scale;  /* measured-anchor multiplier at repair=1
                                    (default .0001) */
    double repair_frame_blend;   /* blend toward the fair tangent-frame target
                                    at repair=1 (default .85) */
    double repair_frame_screen;  /* screened diffusion data term for the fair
                                    frame target at repair=1 (default .002) */
    int    repair_frame_sweeps;  /* red/black screened-frame sweeps (default 128) */
    double axis_y, axis_x;       /* cylindrical scroll axis point (world y,x);
                                    read only when axis_normal_weight>0 (default 0) */
    double axis_normal_weight;   /* soft blend of each vertex frame toward the
                                    cylinder frame (V=z-hat, U=circumferential) so
                                    an axial "bucket-lid" normal is high-energy.
                                    0 = ordinary ARAP, bit-for-bit (default 0) */
    double axis_normal_radius;   /* optional ramp: fade the prior to 0 past this
                                    radius (voxels); 0 = constant weight (default 0) */
    double axis_fold_limit;      /* line-search guard: reject a step that raises a
                                    triangle |n.z|/|n| above this AND increases it.
                                    0 = barrier unchanged (default 0) */
    int    verbose;
} QuadStripArapOpts;

enum {
    QUAD_STRIP_ARAP_STOP_NONE = 0,
    QUAD_STRIP_ARAP_STOP_TOLERANCE = 1,          /* absolute movement_tolerance */
    QUAD_STRIP_ARAP_STOP_RELATIVE_TOLERANCE = 2, /* scale-aware effective tol */
    QUAD_STRIP_ARAP_STOP_STALLED = 3,            /* stall window, converged=0 */
    QUAD_STRIP_ARAP_STOP_MAX_ITERATIONS = 4,
    QUAD_STRIP_ARAP_STOP_BARRIER = 5             /* barrier retained iterate */
};

typedef struct QuadStripArapStats {
    int    iterations;
    int    converged;
    int    stop_reason;                  /* QUAD_STRIP_ARAP_STOP_* */
    double effective_movement_tolerance; /* max(abs, rel*first max move) */
    int    barrier_rejections;
    int    barrier_zero_fallbacks; /* iterations/prolongations that hard-froze
                                      implicated smooth barrier controls */
    int    barrier_stalled;        /* no valid nonzero step remained; current
                                      valid iterate was retained */
    int    maximum_cg_iterations;
    int    multigrid_levels_used;
    int    multigrid_cycles;
    int    nonlinear_coarse_levels;
    int    nonlinear_coarse_iterations;
    int    nonlinear_coarse_caps;
    double final_cg_relative_residual;
    double final_multigrid_residual_before;
    double final_multigrid_residual_after;
    double nonlinear_min_prolongation_scale;
    double maximum_vertex_movement;
    double rms_vertex_movement;
    double minimum_step_scale;
    double initial_edge_log_rms;
    double final_edge_log_rms;
    double initial_edge_log_max;
    double final_edge_log_max;
    size_t repair_vertices;
    double repair_weight_sum;
} QuadStripArapStats;

void QuadStrip_defaults(QuadStripOpts *o);
void QuadStripArap_defaults(QuadStripArapOpts *o);

/* Build the strip.  Inputs describe one component already placed on its own
 * lattice:
 *   domain[H*W]  : nonzero where the fitted field is present (support|interp)
 *   field[H*W*3] : (z,y,x) fitted position per cell (defined where domain != 0)
 *   H, W         : lattice rows (v) and cols (u)
 *   gr0, c0      : lattice origin -> uv is u=(c0+i)*grid_du, v=gr0+j
 *   grid_du      : voxels per u-column
 * Outputs are arena-allocated:
 *   *out_verts[nv*3] (z,y,x), *out_faces[nf*3] (0-based tris), *out_uv[nv*2]
 *   (u,v), *out_filled[nv] (1 = PDE-filled, 0 = sampled).  Any out_* may be NULL
 *   to skip it.  Returns 0 on success, -1 on bad input. */
int QuadStrip_build(Arena_T arena,
                    const uint8_t *domain, const float *field,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled);

/* Cylindrical production entry point.  lifted_phase is optional for Cartesian
 * builds.  With opts->cylindrical_fill it should be finite at every domain
 * cell and carry StrokeStrip's cycle-consistent universal-cover phase.  The
 * solver uses it ONLY to choose the integer 2*pi branch of atan2(field-axis):
 * the observed XYZ angle itself remains exact.  Passing NULL retains the
 * isolated legacy phase-reconstruction fallback for old callers/tests. */
int QuadStrip_build_with_phase(Arena_T arena,
                    const uint8_t *domain, const float *field,
                    const float *lifted_phase,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled);

/* Variant with topology separated from initializer evidence.  `topology`
 * says where the selected RECT/SPAN/BOUNDED policy may emit a surface;
 * `domain` says where field/phase contain Dirichlet initializer samples.
 * Cells admitted by topology but absent from domain are PDE-filled.  Passing
 * NULL topology is identical to QuadStrip_build_with_phase. */
int QuadStrip_build_topology_with_phase(Arena_T arena,
                    const uint8_t *topology,
                    const uint8_t *domain, const float *field,
                    const float *lifted_phase,
                    int H, int W, int gr0, int c0, double grid_du,
                    const QuadStripOpts *opts,
                    float **out_verts, size_t *out_nv,
                    int32_t **out_faces, size_t *out_nf,
                    float **out_uv, uint8_t **out_filled);

/* Relax one complete HxW RECT ribbon in place.  verts are z/y/x, uv are u/v,
 * and filled is 0 for measured/fitted cells and 1 for introduced cells.  The
 * optional repair_weight field is read from opts for exactly this HxW level;
 * nonlinear multigrid restricts it with the geometry. An optional fixed mask
 * supplies exact Dirichlet vertices and disables nonlinear coarsening (the
 * structured linear multigrid preconditioner remains active). The routine never
 * changes topology or UV. With preserve_axial, the z coordinate is an exact
 * per-vertex Dirichlet value while y/x remain softly anchored. Returns 0 for a safe result (including
 * a finite max-iteration stop), -1 for invalid input/allocation/solver failure,
 * and -2 if no orientation-preserving step can be found. */
int QuadStrip_metric_arap(float *verts, const float *uv,
                          const uint8_t *filled, int H, int W,
                          const QuadStripArapOpts *opts,
                          QuadStripArapStats *stats);

int QuadStrip_selftest(void);

#endif
