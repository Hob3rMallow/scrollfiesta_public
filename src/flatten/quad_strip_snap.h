#ifndef QUAD_STRIP_SNAP_INCLUDED
#define QUAD_STRIP_SNAP_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* quad_strip_snap -- PROGRESSIVE relaxation of a solid quad strip onto the real
 * papyrus surface in the RAW CT.
 *
 * One soft, iterated energy over the strip, minimized by structured red/black
 * Gauss-Seidel rounds with an optional aggregate geometric multigrid hierarchy
 * (Jacobi only for an unstructured SPAN fallback).
 * On the first proximal pass every PDE-filled cell is movable.  On later
 * active-set passes, filled cells that already sit on a credible recto
 * transition may become soft zero-depth anchors.  A fitted cell is normally
 * movable only when its signed normal profile proves BOTH that it is not already
 * on a local CT ridge and that a materially stronger ridge exists within reach.
 * One exception is a corner of an emitted triangle whose interior is bake-dark:
 * it is unpinned as a coherent seam repair, and a meaningless near-zero vertex
 * candidate is suppressed so neighbouring supported offsets can quilt across
 * it.  A bake-equivalent +/-2-voxel normal-max probe separately reports which
 * dark sites still have a substantially brighter recto transition farther
 * INWARD.  Burnt/dim fitted sites outside such raster-dark triangles remain
 * trusted.
 * Each round, for every movable cell v, minimize
 *     E = w_target * ||v - target||^2        (snap: nearest bright ridge along n,
 *                                             if a better one is within reach)
 *       + w_orig   * ||v - v_orig||^2         (stay near the harmonic/fit init)
 *       + w_smooth * sum_nbr ||v - v_nbr||^2  (slice-to-slice / circumferential
 *                                             continuity -- move as a sheet)
 *       + w_dev    * sum_nbr ||(v - v_nbr) - Rbar_nbr * (rest_v - rest_nbr)||^2
 *                                            (QUASI-developable ARAP prior: soft,
 *                                             carries the fill toward its own
 *                                             unrolled wrap, never enforces strict
 *                                             isometry -- burnt scrolls are not
 *                                             developable).
 * rest = the flat UV lattice (u, v, 0); Rbar is the average of the two cells'
 * local closest-rotations to that rest.  A movable cell without an accepted raw
 * ridge still relaxes under the sheet terms and is reported separately from
 * confidence-diffused targets; image-space darkness decides whether it is a crack.
 */

typedef struct QuadStripSnapOpts {
    double reach;         /* max march distance along +/- normal, vox (def 6) */
    double step;          /* profile step, vox (def 0.5) */
    double bright_min;    /* min ridge intensity to accept, 0..255; <= 0 derives
                           * a per-scroll threshold from fitted-cell p1/p99 (def auto) */
    double band_frac;     /* auto threshold = p1 + band_frac*(p99-p1) (def 0.30) */
    double min_gain;      /* min intensity gain over the vertex's current sample
                           * before a ridge counts as a real move (def 20) */
    double ridge_lock;    /* fitted cell is trusted/frozen when a local normal-
                           * profile maximum lies within this distance (def 0.75) */
    double dark_probe_frac; /* bake-dark diagnostic threshold within fitted
                             * p1..p99 (def .20; 0 disables). A reported target
                             * additionally requires +/-2 normal-max darkness
                             * and a stronger recto transition strictly inward. */
    double guided_reach; /* wider residual search for bake-dark/ambiguous sites
                           * and hard cap on total position displacement (def 12). */
    double winding_tube; /* minimum half-width of the adaptive same-winding
                          * corridor (def 3 vox).  On a structured cylindrical
                          * ribbon the actual limit is 45% of the measured local
                          * adjacent-wrap spacing, clamped to guided_reach.  Thus
                          * this no longer rejects valid corrections merely
                          * because the initializer missed by >3 voxels, while
                           * two neighbouring turns still cannot exchange order. */
    int    patch_recover; /* post-snap connected-dark-patch recovery.  Dark
                            * continuation cells are buffered, assigned one
                            * globally agreed +/- normal ray direction, snapped,
                            * then locally Laplacian-relaxed into their known rim
                            * (library default off; solid-ribbon production arms). */
    int    patch_buffer;  /* quad/vertex rings around each dark patch (def 2). */
    int    patch_min_cells; /* smallest connected dark-cell region (def 1). */
    double patch_ray_reach; /* one-sided first-material ray reach (def 48 vox). */
    double patch_ray_margin; /* required +/- first-hit distance margin (def 1 vox). */
    int    patch_relax_sweeps; /* red/black screened-Laplacian blend sweeps (def 64). */
    double patch_target_weight; /* direct ray-hit screen weight (def 16). */
    int    patch_fitted_dark_slack; /* fitted bake samples allowed to toggle at
                                    * moved fill seams (def 512). */
    int    local_contrast; /* allow a robust profile-local transition threshold
                            * in burnt/mineralized regions (default 1). */
    int    snap_side;     /* target policy along the oriented profile:
                           * QUAD_SNAP_SIDE_NEAREST (default) = nearest local
                           * maximum, unoriented; QUAD_SNAP_SIDE_RECTO = legacy
                           * inward bright->dark material transition;
                           * QUAD_SNAP_SIDE_MIDLINE = center of the bright
                           * papyrus band between its two dark boundaries (the
                           * sheet midline, NOT the recto face).  The solid-
                           * ribbon CLI defaults to MIDLINE with --cyl-axis. */
    int    preserve_axial;/* keep world z exactly equal to the initializer for
                           * every vertex (default 0; enabled by --cyl-axis).
                           * In cylindrical mode v already determines z, so CT
                           * snapping solves only radius/phase at that slice. */
    int    prevent_intersection_growth; /* exact whole-mesh audit after every
                           * outer round; backtrack/rollback the round until the
                           * conflict count and fold count do not grow (def 1). */
    double filled_target_anchor; /* active-set refinement: multiplier on w_orig
                                  * for PDE-filled cells already on a credible
                                  * oriented target (recto transition or band
                                  * midline, per snap_side); they also supply
                                  * zero-depth evidence but remain sheet-coupled
                                  * (def 0; solid-ribbon uses 4 after pass one) */
    double fitted_source_anchor; /* multiplier on w_orig for upstream support.
                                  * 0 preserves legacy hard freezing; a high
                                  * positive value keeps support soft but strongly
                                  * attached during every iterated CT solve. */
    const float *fitted_source_positions; /* optional persistent nv*3 source
                                           * positions for the attachment above;
                                           * when NULL, attach to this pass input. */
    double axis_y, axis_x;/* scroll-axis point used to orient inward in y/x */
    double quilt_smooth;  /* 3-D offset-field neighbour weight (def 4) */
    double quilt_seam_smooth; /* multiplier on fitted/filled neighbour edges;
                               * 1 is uniform coupling (def 1) */
    double quilt_anchor;  /* weak zero-depth anchor for no-ridge cells (def .01) */
    double quilt_huber;   /* local Huber scale for downweighting a raw target
                           * whose 3-D offset disagrees with adjacent raw targets
                           * (def 1.5 vox; 0 disables robust weighting) */
    double quilt_conf_min;/* min diffused target confidence to apply (def .02) */
    int    quilt_sweeps;  /* red/black offset/confidence sweeps per round (def 4) */
    double w_target;      /* pull toward the ridge target (def 1.0) */
    double w_orig;        /* stay near the init (def 0.25) */
    double w_smooth;      /* per-neighbour continuity (def 0.3) */
    double w_dev;         /* per-neighbour quasi-developable ARAP prior (def 0.4) */
    int    rounds;        /* iterated rounds (def 12) */
    int    solve_sweeps;  /* red/black position sweeps per outer round (def 2) */
    double solve_omega;   /* red/black under-relaxation, 0..1 (def 0.5) */
    int    multigrid_levels; /* structured hierarchy depth including the fine
                              * grid (def 1 = single-level/off; max 16) */
    int    multigrid_cycles; /* coarse corrections per scalar/position solve
                              * after the ordinary fine sweeps (def 1) */
    int    multigrid_patch;  /* coarse-cell width of independently barriered
                              * position corrections (def 16; 0 = one global
                              * transaction). */
    double bake_window_low;  /* fixed RAW window for the in-core raster-style
                              * pass gate (def 47) */
    double bake_window_high; /* fixed RAW window high endpoint (def 193) */
    int    bake_dark_u8;     /* post-window dark cutoff (def 13) */
    double remesh_motion; /* advise regridding when p95 accumulated tangential
                           * motion exceeds this many local UV cells (def .5) */
    double remesh_motion_max; /* localized tangential-motion trigger in UV cells
                               * (def 4; 0 disables) */
    double remesh_sander; /* advise regridding when bake-triangle two-sided Sander
                           * p95 reaches this value (def 2.0) */
    double remesh_over4;  /* localized trigger: fraction of emitted triangles at
                           * >=4x Sander (def .0001 = .01%; 0 disables) */
    int    require_complete_raw; /* fail before solving when any RAW chunk needed
                                  * by the padded geometry bbox is absent (def 1) */
    /* Optional caller-owned (grid_h-1)*(grid_w-1) OR-accumulators for the
     * CT-occupancy continuation trim.  Filled by the existing probes across
     * EVERY pass and round (including rejected ones, biasing toward KEEP):
     * cell_lit_accum[cell]=1 when any bake-style triangle-interior sample in
     * the cell was non-dark (CT material seen); cell_candidate_accum[cell]=1
     * when a dark sample earned a supported bright alternative within the
     * winding corridor (a recoverable "glimmer").  A filled cell with
     * neither is harmonic fill over black CT.  NULL disables (default). */
    uint8_t *cell_lit_accum;
    uint8_t *cell_candidate_accum;
    int    verbose;
} QuadStripSnapOpts;

enum {
    QUAD_STRIP_SNAP_OK = 0,
    QUAD_STRIP_SNAP_ERROR = -1,
    QUAD_STRIP_SNAP_RAW_INCOMPLETE = -2
};

enum {
    QUAD_SNAP_SIDE_NEAREST = 0,
    QUAD_SNAP_SIDE_RECTO = 1,
    QUAD_SNAP_SIDE_MIDLINE = 2
};

/* RMS physical 3-D offset jumps on axial lattice edges, partitioned by whether the
 * edge lies inside trusted fitted evidence, inside harmonic/refit fill, or
 * crosses their interface.  `all_*` is the count-weighted union of the three
 * disjoint populations (not an unweighted mean of their RMS values). */
typedef struct QuadStripSnapRegionEdgeStats {
    double all_rms;
    double fitted_fitted_rms, filled_filled_rms, fitted_filled_rms;
    size_t all_edges;
    size_t fitted_fitted_edges, filled_filled_edges, fitted_filled_edges;
} QuadStripSnapRegionEdgeStats;

/* Measurements for the mesh that leaves one proximal pass.  These are kept
 * separate from the remesh trigger so the caller can treat a pass (including a
 * preceding bounded remesh) as a transaction and reject a geometric
 * regression.  CT percentiles are direct samples at final vertex positions;
 * fitted and harmonic-fill populations are reported independently because a
 * global brightness gain can otherwise hide a worse extrapolated region. */
typedef struct QuadStripSnapStats {
    double displacement_rms;
    double tangential_p95, tangential_max;
    /* Primary distortion statistics are evaluated on the two actual triangles
     * emitted for every structured quad, using the same singular-value
     * quantization as obj_bake_raw.  A quad-averaged Jacobian can cancel an
     * alternating/checkerboard skew and materially undercount the baked tail. */
    double sander_p95, sander_p99, sander_max;
    double over2_fraction, over4_fraction;
    size_t sander_samples, quads, flipped, degenerate;

    /* Retain the former averaged-quad statistic as a diagnostic.  It is useful
     * for distinguishing broad metric drift from triangle-local skew, but is no
     * longer used by the remesh or transaction acceptance gates. */
    double quad_sander_p95, quad_sander_p99, quad_sander_max;
    double quad_over2_fraction, quad_over4_fraction;

    int snapped, nopap;
    int raw_ridge, recto_ridge, near_zero_ridge;
    int raw_no_candidate, raw_candidate_applied, quilt_only_applied;
    int guided_ridge, wide_ridge, local_contrast_ridge, winding_rejected;
    /* Midline mode only: targets that used the merged/clipped-band smoothed-
     * maximum fallback, and movable oriented-mode vertices whose band finder
     * abstained entirely. */
    int band_fallback, band_abstain;
    int dark_seed_forced_movable, dark_seed_quilt_only;
    size_t dark_face_samples, dark_face_candidates;
    size_t dark_face_outward_candidates, dark_face_no_candidate;
    size_t adaptive_winding_vertices;
    double adaptive_winding_p05, adaptive_winding_p50, adaptive_winding_p95;
    double axial_drift_max;

    /* Connected dark-patch recovery is deliberately reported separately from
     * pointwise ridge quilting. `border_consensus` means the buffered known rim
     * belongs to one connected trusted XYZ region and may select a side by a
     * strict ray majority. Ambiguous rims require every unknown-vertex ray to
     * agree on one oriented side; contradictory votes are `wobbly` and never
     * averaged. */
    size_t patch_regions, patch_regions_small;
    size_t patch_border_consensus, patch_ray_consensus;
    size_t patch_border_selected, patch_ambiguous_selected;
    size_t patch_wobbly, patch_no_consensus, patch_overlap_skipped;
    size_t patch_proposed, patch_applied, patch_attenuated, patch_geometry_rejected;
    size_t patch_vertices, patch_vertices_retained;
    size_t patch_dark_before, patch_dark_after;
    int patch_transaction_accepted;
    double patch_retained_scale;

    /* Whole-component exact triangle-conflict audit performed by the snap
     * transaction itself.  The final committed census is reusable by the
     * owning solid-ribbon pass gate; rebuilding identical pre/post BVHs there
     * would add cost without adding evidence. */
    int intersection_audit_complete;
    size_t input_intersections, output_intersections;
    size_t input_overlap, input_stab, input_fold;
    size_t output_overlap, output_stab, output_fold;
    int intersection_guard_audits, intersection_guard_backtracks;
    int intersection_guard_round_rollbacks;
    int intersection_guard_orientation_rejects;
    size_t intersection_guard_locally_rolled_back_vertices;
    int quilt_downweighted;
    double quilt_candidate_weight_mean;
    double quilt_raw_edge_rms, quilt_edge_rms;
    size_t quilt_raw_edges, quilt_edges;
    double quilt_seam_smooth;

    /* The legacy/global pair above deliberately retains its transaction-gate
     * semantics: raw requires two ridge samples and quilted requires two active
     * targets.  These three richer diagnostics expose the boundaries it omits:
     *
     * raw_regions:     both endpoints have an accepted raw ridge depth;
     * solved_regions:  every edge touching the scalar solve, with frozen fitted
     *                  endpoints represented by their zero-depth Dirichlet value;
     * applied_regions: every edge touching an active target after the confidence
     *                  cutoff, with an inactive endpoint represented as zero.
     *
     * In particular, solved/applied fitted_filled_rms measures the seam between
     * proxy-supported papyrus and PDE/bake-refit reconstruction. */
    QuadStripSnapRegionEdgeStats quilt_raw_regions;
    QuadStripSnapRegionEdgeStats quilt_solved_regions;
    QuadStripSnapRegionEdgeStats quilt_applied_regions;

    size_t fitted_ct_samples, filled_ct_samples;
    int fitted_ct_p1, fitted_ct_p25, fitted_ct_p50;
    int fitted_ct_p75, fitted_ct_p99;
    int filled_ct_p1, filled_ct_p25, filled_ct_p50;
    int filled_ct_p75, filled_ct_p99;

    /* Pass-input routing diagnostics.  "Bake-dark" uses the same +/-2 voxel,
     * five-tap normal-max support probe as obj_bake_raw, thresholded at
     * fitted_p1 + dark_probe_frac*(fitted_p99-fitted_p1).  Only the much
     * smaller recto-target populations are eligible to leave their current
     * surface; all other dark sites remain honest crack/ink/void candidates. */
    double dark_probe_threshold;
    size_t fitted_bake_dark, fitted_dark_recto_target;
    size_t filled_bake_dark, filled_dark_recto_target;

    /* Fixed-window, raster-style samples at the two real triangle interiors of
     * every structured quad.  Unlike dark_probe above, these compare the mesh
     * entering and leaving the proximal pass and therefore participate in its
     * transaction gate.  Filled/fitted classification uses the same >=2 of 4
     * filled-corner vote as audit_quad_strip_bake.py. */
    size_t input_bake_fitted_samples, input_bake_filled_samples;
    size_t input_bake_fitted_dark, input_bake_filled_dark;
    size_t input_bake_missing;
    size_t output_bake_fitted_samples, output_bake_filled_samples;
    size_t output_bake_fitted_dark, output_bake_filled_dark;
    size_t output_bake_missing;
    double input_sander_p95, input_sander_p99, input_over4_fraction;

    /* RAW coverage is measured before any geometry operation.  Production
     * callers require missing == 0; the counts remain in reports so a
     * diagnostic opt-out cannot masquerade as a complete run. */
    size_t raw_chunks_expected, raw_chunks_loaded, raw_chunks_missing;
    size_t raw_chunks_outside_volume;

    /* Last-round structured multigrid evidence.  A level count of one means
     * the historical fine-grid-only solve.  Residuals are Euclidean relative
     * residuals over all solved channels.  Position correction counts exclude
     * the ordinary barriered fine-grid red/black sweeps. */
    int quilt_mg_levels, quilt_mg_corrections;
    int position_mg_levels, position_mg_corrections;
    int position_mg_backtracks, position_mg_rejects;
    int position_mg_patches_accepted, position_mg_patches_rejected;
    double quilt_mg_residual_before, quilt_mg_residual_after;
    double position_mg_residual_before, position_mg_residual_after;
} QuadStripSnapStats;

void QuadStripSnap_defaults(QuadStripSnapOpts *o);

/* Relax the strip in place.  verts[nv*3] (z,y,x), faces[nf*3] 0-based,
 * filled[nv] (1 = a bounded-hole loft cell, 0 = certified/source geometry),
 * uv[nv*2] (the flat-UV rest for the developable term).  Upstream geometry is
 * frozen unconditionally; only bounded-hole loft vertices are movable.  With
 * snap_side RECTO, a
 * credible bright-to-dark transition along the inward-oriented normal is the
 * preferred target; with snap_side MIDLINE the target is the center of the
 * bright band between its two dark boundaries; the nearest local maximum is
 * retained only as an ambiguous-
 * evidence fallback.  With filled_target_anchor > 0, PDE-filled sites already
 * within ridge_lock of that target become stronger soft anchors, so an
 * iterative pass concentrates motion on the unresolved active set without
 * cutting recovered material out of the sheet solve.  raw_dir =
 * cubes_RAW or an uncompressed uint8 OME-Zarr,
 * chunk = cube edge (vox).  out_snapped / out_nopap (nullable) receive the last
 * round's movable-with-target / movable-without-an-accepted-raw-ridge counts.
 * The latter is measured before confidence diffusion and is not an assertion
 * that the corresponding raster must be physically empty or black. out_stats (nullable)
 * receives final geometry, quilting, and population-separated CT evidence for
 * transactional pass acceptance.  Returns QUAD_STRIP_SNAP_OK on success,
 * QUAD_STRIP_SNAP_RAW_INCOMPLETE when strict RAW coverage fails, and
 * QUAD_STRIP_SNAP_ERROR on bad input / RAW table failure. */
int QuadStripSnap_run(Arena_T arena, float *verts, size_t nv,
                      const int32_t *faces, size_t nf, const uint8_t *filled,
                      const float *uv, int grid_h, int grid_w,
                      const char *raw_dir, long chunk,
                      const QuadStripSnapOpts *opts,
                      int *out_snapped, int *out_nopap,
                      int *out_remesh_advised,
                      QuadStripSnapStats *out_stats);

/* Audit trail for the alternating structured correspondence step.  The UV
 * lattice itself stays regular; changing where its vertices sample the current
 * bilinear surface is the equivalent bounded UV/reparameterization update.
 * A sweep is applied only when exact emitted-triangle Sander statistics are
 * Pareto-nonworse and the shared map improves at least one distortion or edge-
 * distribution statistic. */
typedef struct QuadStripRemeshStats {
    int applied_sweeps;
    int used_robust_metric;
    int used_v_first;
    double alpha;
    double before_p95, before_p99, after_p95, after_p99;
    double before_edge_log_rms, after_edge_log_rms;
    size_t triangles;
    size_t before_over2, before_over4, after_over2, after_over4;
} QuadStripRemeshStats;

/* Search smooth shared row/column correspondence maps in both axis orders and
 * at bounded global step sizes.  Robust winsorized and ordinary mean edge metrics
 * are both considered.  Unlike the former per-vertex backtracking remesher, a
 * candidate is committed as one coherent map or not at all, preventing a
 * checkerboard of unrelated local step fractions.  `uv` is the unchanged flat
 * rest lattice used for the bake-faithful triangle objective.  Returns 0 for a
 * valid search (including a safe no-op), -1 on invalid input/OOM. */
int QuadStripSnap_remesh_structured(float *verts, uint8_t *filled,
                                    const float *uv, int H, int W, int sweeps,
                                    QuadStripRemeshStats *out_stats);

int QuadStripSnap_selftest(void);

#endif
