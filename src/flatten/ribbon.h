#ifndef RIBBON_INCLUDED
#define RIBBON_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* Dense storage is limited by the checked byte budget. Bounded export keeps
 * one grid resident, so collision depth need not share the old eight-grid cap. */
#define RIBBON_MAX_PEEL_LAYERS 64

/* ============================================================================
 * ribbon.h -- slice / arc-length / joint "StrokeStrip" parameterization of a
 * scroll-surface mesh, plus a clean fitted ribbon surface.
 *
 * Input: one welded recto sheet (ideally a single connected component; run
 * obj_biggest first). No winding field is required -- the winding is derived
 * in the SLICE domain, where residual inter-wrap fusion bridges (which
 * collapse any mesh-graph-integrated angle field) can be cut geometrically.
 *
 * Pipeline (StrokeStrip, Pagurek van Mossel et al. 2021, fused with the
 * diffeomorphic-spiral wrap indexing idea):
 *   A. cut the mesh with planes perpendicular to the scroll axis (thin
 *      slices); chain the per-plane triangle crossings into polylines keyed
 *      by mesh-edge identity (exact, no epsilon matching); unwind the local
 *      angle theta along each chain; CUT chains at fusion-bridge crossings
 *      (radial jumps far steeper than any spiral ramp);
 *   B. match samples across adjacent slices (WTS-adjacent pairs; distance
 *      gate < half the 7-vox inter-wrap clearance makes cross-wrap pairing
 *      impossible) + intra-slice fragment continuations;
 *   W. obtain a continuous universal-cover phase.  Preferred: interpolate the
 *      caller's gauge-synchronized reference_phi (continuation + radial layer
 *      order).  Legacy fallback: infer integer chain turns from local slice
 *      geometry.  In either case phi is an ORDER scaffold, not material U.
 *   C. jointly solve for u per sample with the StrokeStrip discretization:
 *      one average-before-square arc-length row and one latent-isovalue
 *      alignment term per maximal trusted V-connected run, weak local metric
 *      rows, physical fragment-continuation rows, monotone repair, and robust
 *      membership reweighting.  A missing V connection starts a new run;
 *      disconnected observations are never averaged together;
 *   D. transfer (u, v = axial) to every mesh vertex via the nearest slice
 *      sample from the same connected source chart and lifted-winding sheet
 *      (exact provenance + reference-phi gates);
 *   E. fit the clean ribbon: a regular (u, v) grid of 3D positions smoothed
 *      per row (tangent then position tridiagonal solves) and filled/
 *      regularized across rows.
 *
 * u is TRUE per-slice arc length (vox), unlike unwrap.c's winding*r_ref
 * scaffold.  A monotone winding integral may pin one additive constant per
 * disconnected metric solve; it never replaces the measured arclength rows.
 * When no safe cover gauge exists, detached solves are aligned to the solved
 * main sheet's empirical phi->u relation.
 * v is the axial coordinate (vox). Both start at 0.
 * ==========================================================================*/

typedef struct {
    float  axis_point[3];  /* (z,y,x) point on the scroll axis (umbilicus) */
    float  axis_dir[3];    /* (z,y,x) axis direction (normalized internally) */
    float  slice_h;        /* slicing plane spacing, vox (default 2.0) */
    float  sample_h;       /* polyline resample step, vox (default 2.0) */
    float  match_r;        /* cross-slice in-plane match radius, vox (default
                            * 3.0 -- below half the 7-vox inter-wrap clearance
                            * so a match can never join two wraps) */
    float  match_ang_deg;  /* tangent-similarity gate, deg (default 20) */
    int    scaffold_solve; /* Stage-C u is only a SCAFFOLD for a downstream
                            * re-parameterization (the quadribbon pipeline's
                            * fit stage): cap the monotone ADMM at a few
                            * iterations instead of polishing a coordinate
                            * that will be re-solved.  Claims, winding
                            * registration, and identities are unaffected. */
    int    relax_iters;    /* robust membership reweight rounds (default 5;
                             * first 3 smoothed-L1, then Gaussian likelihood) */
    int    final_iters;    /* frozen-membership tangent-align rounds (default 2) */
    int    metric_iters;   /* post-alignment continuation rounds which keep the
                             * final isovalue cross-sections fixed and ramp the
                             * physical stroke-edge unit-speed term.  Zero by
                             * default; topology-preserving parameterization
                             * enables three unless explicitly overridden. */
    float  metric_weight;  /* terminal unit-speed weight (default 4.0) */
    int    solve_threads;  /* Stage-C parallel workers; 0 = OpenMP runtime max */
    int    solve_amg;      /* geometric semi-coarsening + FMG solve (default 1;
                            * name retained for CLI compatibility) */
    const char *solve_level_prefix; /* optional human-inspection checkpoint
                             * prefix. Every robust solve round publishes its
                             * full fitted quad-ribbon immediately; the final
                             * robust solve also writes the coarse-to-fine FMG
                             * levels. NULL disables checkpoints. */
    int    solve_level_write_obj; /* also serialize text OBJ for every checkpoint;
                            * default 0 keeps only compact authoritative VMESH */
    float  grid_u;         /* ribbon grid spacing along u, vox (default 2.0) */
    int    fit_ribbon;     /* nonzero: produce the ribbon grid (default 1) */
    int    direct_ribbon;  /* nonzero: fit the quad ribbon directly from the
                              solved slice samples.  Skip the legacy source-mesh
                              UV transfer, gauge reconciliation, and atlas
                              packing stages. */
    int    stitch_solve_gauges; /* nonzero: order disconnected Stage-C solve
                              * components by the lifted-phase scaffold and
                              * concatenate only their additive U intervals.
                              * No U scaling and no V change are permitted. */
    int    fit_cover_width;/* diagnostic: size the fitted grid from pre-atlas
                             * cover U even when metric islands are compacted */
    int    verify_fit_width;/* diagnostic: fit compact and cover-padded grids
                             * from one solved SliceSet, print stage digests,
                             * and compare their shared metric-atlas prefix */
    float  wrap_spacing;   /* local radial pitch (vox/turn), used by the legacy
                             * winding fallback and diagnostics/continuation
                             * gates. It never converts radius to material U.
                             * <=0 = auto-estimate (default 0). */

    /* --- per-cube global-frame emission (scroll_whole; all default OFF, in
     * which case behavior is byte-identical to the classic single-mesh run) --- */
    int    emit_global;    /* nonzero: do NOT shift (u,v) to start at 0 -- emit
                            * u in the solve's absolute frame and v = raw axial
                            * coordinate t (world z for axis_point z=0, axis Z),
                            * fill RibbonResult.phi, and set u_origin. The
                            * ribbon grid (if fit) stays origin-based. */
    int    winding_sense;  /* 0 = auto. +1/-1 pins the outward direction of the
                             * winding scaffold so every cube agrees. */
    int    pin_orient;     /* nonzero: skip the canonical u-flip (u increases
                            * outward). The flip is per-mesh (cov(u,r) sign), so
                            * cubes of one scroll would mirror independently. */
    double spiral_a, spiral_b;  /* spiral_b != 0: explicitly pin a global spiral
                            * r ~= a + b*phi/2pi instead of fitting it, forcing
                            * every cube through the SAME arc-length init
                            * u(phi) = a*phi + b*phi^2/(4pi) (and disabling the
                            * per-component registration fallback). Calibrate
                            * (a,b) once on a seed cube, pin everywhere. */
    int reference_phi_authoritative; /* nonzero: reference_phi is the carried,
                              * gauge-synchronized registration (a sidecar from
                              * the fit), valid across disconnected components.
                              * Zero: a recomputed unwrap fallback whose phase
                              * is unanchored on detached components -- the
                              * metric projection then keeps every island's
                              * carried placement instead of phase-placing. */
    const float *reference_phi; /* optional [nv] graph-traced winding field.
                            * When supplied, slice crossings interpolate this
                            * authoritative absolute phi instead of rebuilding
                            * turns from cross-slice nearest-neighbour votes.
                            * StrokeStrip still solves true arc length; this is
                             * only its winding/turn scaffold. The array must
                             * stay alive for the duration of Ribbon_run. */
    const double *reference_u; /* optional [nv] already-solved material-U field.
                             * When supplied, slice crossings interpolate this
                             * field and Stage E fits the regular ribbon directly;
                             * the pair graph, winding assignment, metric solve,
                             * and mesh-transfer stages are skipped.  This is the
                             * fast path for a caller which has already lifted a
                             * certified chart atlas through welded topology.
                             * reference_phi should also be supplied when
                             * component_global fusion is requested: without it
                             * the fast path carries phi == u and the junction
                             * phase gate is skipped.
                               * The array must stay alive for Ribbon_run. */
    const double *reference_v; /* optional [nv] authoritative axial coordinate
                             * accompanying reference_u. projective_grid uses
                             * reference_v - dot(vertex-axis_point, axis_dir)
                             * to anchor one crop-independent slice lattice. */
    const float *reference_u_confidence; /* optional [nv] posterior confidence
                               * for the ORDER carried by reference_u, in [0,1].
                               * It weights robust per-chain gauge votes; it
                               * never clips individual samples or changes an
                               * exact within-chain metric derivative.
                               * The array must stay alive for Ribbon_run. */
    int solve_reference_u;    /* nonzero: reference_u is correspondence and
                              * initial-gauge evidence for the full Stage-C
                              * StrokeStrip solve, rather than the solved-U fast
                              * path. Geometry still determines the output U. */
    int coarse_chain_seed;    /* Opt-in carried-U initializer: first solve one
                              * offset per physical slice chain, then release
                              * every sample to the unchanged Stage-C solve.
                              * No source topology or observations are removed. */
    int projective_grid;      /* nonzero: reference_u/reference_v define an
                              * immutable signed atlas frame. Quantize the
                              * internal dense grid to that frame, restore its
                              * origins on emitted UV, retain empty U columns,
                              * and never promote peel claims into primary U. */
    int coherent_claims;      /* Carried UV lattice with independent measured
                              * cross-row continuity per peel. Preserve actual
                              * reconstruction lanes separately from lineage. */
    int metric_project_only;  /* nonzero (with reference_u + solve_reference_u):
                              * the parameterization is the direct metric
                              * projection u = s + median(reference_u - s) per
                              * slice chain -- exact XYZ arclength placed by the
                              * carried frame.  NO solver runs: no chain-gauge
                              * solve, no Stage-C StrokeStrip solve, no
                              * component interval stitching, no phi->u
                              * registration, no orientation flip.  Every gauge
                              * is a measured median of carried evidence, never
                              * a solved unknown, so the carried frame is the
                              * sole placement authority.  Pairs and the winding
                               * scaffold are still built for the correspondence
                               * audit and the transfer gates. */
    int preserve_input_rows; /* nonzero (with reference_u + solve_reference_u):
                              * use the parameterized input's own constant-V
                              * vertex rows as StrokeStrip observations and
                              * write solved U straight back to those vertices.
                              * This preserves topology but still runs the full
                              * staged global solve; unlike metric_project_only
                              * it never imposes an independent gauge per row. */
    const double *reference_atlas_shift; /* optional [nv] additive atlas shift
                             * accompanying reference_u.  Subtracting it from
                             * packed U recovers the pre-pack cover coordinate
                             * used by component-global phi consensus.  It must
                             * be constant on each connected source component;
                             * NULL means zero shift. */
    const int32_t *reference_material_island; /* optional [nv] material-domain
                             * identity accompanying reference_phi.  This must
                             * come from positive same-sheet continuation
                             * evidence only: radial-order/gauge observations
                             * must never merge these labels. Different labels
                             * are gated before claimant selection and packed
                             * into disjoint atlas intervals. */
    size_t reference_material_island_count; /* dense label upper bound for
                              * ordinary fits. Ignored for a projective fit
                              * carrying canonical owner provenance, whose
                              * reversible i31 identities are sparse. */
    const int32_t *reference_island; /* legacy alias for
                             * reference_material_island. New callers should
                             * use the explicit material field; supplying both
                             * is an error. */
    size_t reference_island_count; /* number of compact identities in
                               * reference_island. Required and nonzero when
                               * reference_island is supplied. */
    const int32_t *reference_chart; /* optional [nv] exact source-chart id.
                              * It must be constant on every connected input
                              * component.  When direct relations are supplied,
                              * cross-chart geometric correspondences are legal
                              * only for an explicitly serialized graph edge. */
    size_t reference_chart_count; /* dense label upper bound; ignored for the
                              * sparse canonical-owner contract above */
    const int32_t *reference_owner_block; /* optional [nv*3] canonical block
                              * coordinate (iz,iy,ix) for every source vertex.
                              * It must be constant on each connected input
                              * component.  In projective-grid claims mode the
                              * closest ancestor block owns an overlap before
                              * any crop-wide geometric path objective is
                              * considered; later-block alternatives remain in
                              * peel layers.  This is the finite-domain
                              * certificate which makes a nested fit immutable. */
    const int32_t *reference_relation_a; /* optional direct graph endpoints */
    const int32_t *reference_relation_b;
    const double *reference_relation_support; /* optional positive evidence */
    size_t reference_relation_count;
    int    component_global; /* nonzero: multi-chart fitted-grid ownership --
                             * slice-chain enumeration order must not choose a
                             * different wrap in each row.  HOW ownership is
                             * decided is ownership_construction below. */
    int    peel_layers; /* claims-mode collision depth.  0 selects the
                              * production default; 1 disables collision
                              * peeling; 2..8 retain that many independently
                              * selected claimants per UV cell in stacked
                              * temporary bands.  Substantial recovered runs
                              * can then be packed beside the primary atlas. */
    int    peel_layer;  /* -1 (default): retain every requested peel in one
                              * resident grid.  >=0: deterministically consume
                              * and discard every preceding claimant layer but
                              * retain only this layer.  This is the bounded-
                              * memory serialization path for projective owner
                              * certificates; it changes storage, never UV or
                              * claimant selection. */
    int    ownership_construction; /* component_global mode switch.
                             * 0 (DEFAULT): exact whole-row claimant paths.
                             * Physical continuity is lexicographic, followed
                             * by the number of source-chain/chart switches and
                             * then geometric/phase energy.  The selected path
                             * is frozen against post-fit claimant substitution.
                             * 1 (2026-08-13 experiment): ownership BY
                             * CONSTRUCTION -- certified hard chain merges plus
                             * a deterministic row-invariant fusion key; nothing
                             * is searched; losers beyond the wrap gate are
                             * reported violations; local claimant repairs are
                             * disabled and chart-run provenance is exported.
                             * Measured on the 4x5x5 0097 weld: construction
                             * skips 5x more cross-wrap triangles (58,507 vs
                             * 11,578) and visibly tears the texture. */
    int    discard_conflicting_claims; /* component_global claims-mode policy.
                             * Nonzero: a fitted-grid cell with more than one
                             * physical claimant has NO authoritative sample.
                             * All claimants are discarded, post-claim/local
                             * repair is disabled, and the cell may only be
                             * reconstructed as an unsupported variational/
                             * minimal-surface continuation from unambiguous
                             * neighbours.  Default 0 preserves legacy runs. */
    int    reference_anchor_gauges; /* nonzero: reference_u carries an already
                             * globally consistent graph registration.  Preserve
                             * exact measured arclength within each slice chain,
                             * but solve all chain gauges together with one
                             * observation-count-weighted anchor per chain.  This
                             * removes the free additive modes which let the old
                             * one-pin-per-component solve re-scramble a correct
                             * field while still allowing short, phase-certified
                             * 3-D correspondences to reconcile adjacent rows. */
    int    v_strip_planes;   /* > 0: partition the fitted grid's slice planes
                             * into axial strips of this many planes; strips
                             * fit in parallel (disjoint G rows, own scratch)
                             * and the claim selector's preceding-row term
                             * resets at strip boundaries.  < 0 = auto-size for ~2
                             * strips per OpenMP thread (clamped to [4,64]
                             * planes).  0 = one legacy strip. */
    int    radial_bridge_cut; /* -1: never split a slice chain merely because
                              * it runs radially; 0 (default): use that legacy
                              * fusion-wall heuristic only when no authoritative
                              * reference_phi exists; +1: force it on.  A wound
                              * but non-circular/deformed sheet legitimately has
                              * near-radial runs, so reference-guided flattening
                              * must preserve the mesh connectivity. */
} RibbonOpts;

/* Fill opts with the defaults above (axis = Z through the origin; callers
 * must set axis_point/axis_dir for real scrolls). */
void RibbonOpts_default(RibbonOpts *opts);

enum {
    RIB_PITCH_FALLBACK  = 0, /* estimator had insufficient support */
    RIB_PITCH_ESTIMATED = 1, /* measured from adjacent-wrap radial gaps */
    RIB_PITCH_PINNED    = 2  /* caller supplied opts->wrap_spacing */
};

/* Compressed exact claimant provenance for one fitted-grid row.  Runs never
 * span an unsupported column: a gap ends the run even when the same source
 * chart resumes on its other side. */
typedef struct {
    size_t first_col;
    size_t last_col;
    int32_t source_chart;
} RibbonChartRun;

/* Triangles ONE fitted-grid cell emits.  Corners are
 *     a = (k,j)      b = (k,j+1)
 *     c = (k+1,j)    d = (k+1,j+1)
 * A cell commits to ONE diagonal -- MAIN (b-c) or ANTI (a-d) -- never a
 * mixture, so the emitted surface stays a valid quad ribbon.
 *
 * This exists because the writer used to emit two fixed-diagonal triangles
 * whenever all four corners were PRESENT, with no other test.  On the 4x5x5
 * that drew 35,447 edges across the dark gap between physical wraps: the
 * radial spokes and staircase jogs visible in every cross-section.  The
 * parameterizer owns topology, so it now decides which triangles are legal
 * and the writer serializes that decision. */
enum RibbonCellCode {
    RIB_CELL_NONE = 0,
    RIB_CELL_ABC  = 1,   /* main diagonal, upper triangle (a,b,c) */
    RIB_CELL_BDC  = 2,   /* main diagonal, lower triangle (b,d,c) */
    RIB_CELL_MAIN = 3,   /* both -- byte-identical to the historical writer */
    RIB_CELL_ABD  = 4,   /* anti diagonal, upper triangle (a,b,d) */
    RIB_CELL_ADC  = 8,   /* anti diagonal, lower triangle (a,d,c) */
    RIB_CELL_ANTI = 12
};

typedef struct {
    /* --- per original mesh vertex (arena-allocated) --- */
    float   *uv;        /* [nv*2]: u = slice arc length (vox), v = axial (vox);
                         * both shifted to start at 0 (absolute when
                         * opts->emit_global) */
    uint8_t *uv_ok;     /* [nv]: 1 = direct slice-map transfer, 0 = neighbor-
                         * filled or unmapped. Filled values are retained for
                         * diagnostics, but are not safe face geometry. */
    float   *phi;       /* [nv] winding coordinate th + 2pi*W per vertex, from
                         * the nearest slice sample (neighbor-median filled like
                         * u). NULL unless opts->emit_global. Filled vertices
                         * remain uv_ok == 0; wholly unmapped vertices carry 0. */
    int32_t *group;     /* [nv] pair-graph winding-group id (0..w_prior_groups-1)
                         * of the nearest slice sample; -1 where unmapped. The
                         * radius anchor rounds ONE integer offset per group, so
                         * a whole-turn error is constant per group -- cross-cube
                         * registration corrects per group, not per cube. NULL
                         * unless opts->emit_global. */
    double   u_origin;  /* sample-min u (the shift emit_global skips); 0 when
                         * the run had no usable samples */

    /* --- ribbon grid (row-major: row k = slice, column j = isovalue) --- */
    float   *grid_pos;   /* [nk*nu*3] finite positions (z,y,x).  Occupancy is
                           carried only by grid_present; geometry values are
                           never validity sentinels. */
    uint8_t *grid_present;/* [nk*nu]: 1 = this ribbon vertex exists */
    uint8_t *grid_valid; /* [nk*nu]: 1 = supported by slice data (pre-fill) */
    float   *grid_phi;   /* [nk*nu]: finite gauge-synchronized lifted phase in
                          * radians. grid_valid is the measured/generated
                          * marker and preserves the authority distinction. */
    int32_t *grid_island;/* [nk*nu]: branch-aware reconstruction-component
                          * identity selected by the fitted row; -1 for
                          * unresolved derived cells.  Components are the
                          * independently reconstructable/sheet-output units. */
    int32_t *grid_material;/* [nk*nu]: original winding/material lineage for
                            * the selected claimant; -1 for unresolved derived
                            * cells.  Distinct grid_island lanes with the same
                            * lineage may be joined later, but only by the
                            * global evidence-based stitcher. */
    size_t  *grid_chart_row_offsets; /* [nk+1] into grid_chart_runs; NULL when
                                      * exact graph provenance is unavailable */
    RibbonChartRun *grid_chart_runs; /* graph-selected source-chart ownership;
                                      * unsupported/derived cells are absent */
    size_t   grid_chart_run_count;
    int32_t *source_chart_island; /* [source_chart_island_count] exact
                                   * material-domain identity for each
                                   * source chart; NULL without graph input */
    size_t   source_chart_island_count;
    uint8_t *grid_quad;  /* [(nk-1)*(nu-1)]: RIB_CELL_* triangle set this cell
                          * emits.  The writer serializes it verbatim and
                          * validates the corner ids; it applies no geometry
                          * policy of its own. */
    uint8_t *grid_cut_h; /* [nk*nu]: 1 = the lattice edge (k,j)-(k,j+1) is a
                          * topological CUT, not merely an unfilled gap.  A cut
                          * edge can never be crossed by a face, even if a
                          * later pass moves geometry back under the gate. */
    uint8_t *grid_cut_v; /* [nk*nu]: same for (k,j)-(k+1,j) */
    size_t   nu, nk;     /* grid columns (u) and rows.  nk is
                          * slices * grid_layers: the row axis carries one
                          * BAND per peel layer, so a cell that several wraps
                          * claim keeps them all instead of admitting one and
                          * discarding the rest. */
    size_t   grid_layers;/* peel layers stacked along the row axis */
    size_t   grid_layer_base;/* semantic peel index represented by row band 0.
                              * Zero for an ordinary/all-layer result; nonzero
                              * only for bounded single-peel serialization. */
    size_t   grid_source_layers;/* requested peel depth before selecting a
                                 * resident subset; zero means grid_layers */
    float    grid_du;    /* u spacing (== opts->grid_u) */
    float    grid_dv;    /* v spacing (== opts->slice_h) */
    double   grid_u_origin;/* signed U of internal grid column zero */
    double   grid_v_origin;/* signed V of internal grid row zero. In the
                            * legacy projective sampler this labels a bin's
                            * lower boundary, while XYZ is sampled at its
                            * center: source reference V = exported V + dv/2.
                            * Keep this gauge conversion explicit in audits. */
    int      grid_projective;/* writer preserves both origins and every
                              * intervening column; no crop-local compaction */

    /* --- diagnostics --- */
    int     n_slices;        /* slicing planes with any crossing */
    int     n_chains;        /* total polylines over all slices (post bridge cut) */
    int     n_closed;        /* closed loops (cut open) */
    int     n_multi_slices;  /* slices with >1 chain (fragmented) */
    size_t  n_samples;       /* resampled points (QP variables) */
    size_t  n_pairs;         /* cross-slice alignment pairs */
    size_t  n_cont_pairs;    /* intra-slice fragment continuation pairs */
    size_t  n_strip_runs;    /* maximal one-to-one V-connected StrokeStrip rows */
    size_t  n_strip_members; /* slice observations in those rows */
    size_t  n_strip_links;   /* trusted cross-slice links retained in rows */
    size_t  n_strip_links_pruned;/* non-mutual/branching links excluded */
    double  strip_length_rms;/* final RMS of averaged derivative minus one */
    double  strip_align_rms; /* final RMS latent-isovalue member residual */
    size_t  strip_gauge_components;/* solve intervals phase-stitched in U */
    double  strip_gauge_span;/* resulting phase-ordered interval span */
    double  strip_gauge_max_shift;/* largest component-constant U shift */
    size_t  bridge_cuts;     /* chain splits at fusion-bridge radial jumps */
    size_t  mono_repairs;    /* samples moved by the monotone (PAVA) repair */
    size_t  certificate_samples;/* high-confidence reference-U constraints */
    size_t  certificate_clamps;/* metric samples clipped to their trust region */
    double  certificate_drift_p95;/* pre-clip |U-U_certificate_map| p95 */
    double  certificate_drift_max;/* pre-clip maximum drift */
    double  certificate_bound_p50;/* median active trust radius (vox) */
    size_t  uv_filled;       /* no direct sample; u copied from mesh neighbours */
    size_t  uv_fallback;     /* verts still unmapped after neighbour fill */
    size_t  uv_phase_rejects;/* geometrically near sample candidates rejected
                              * because lifted phase selected another ply */
    size_t  uv_phase_blocked;/* vertex/slice queries with nearby samples but no
                              * candidate on the vertex's lifted-phase sheet */
    size_t  uv_component_rejects;/* nearby candidates from another source-mesh
                                  * chart, rejected when own-chart samples exist */
    size_t  uv_component_blocked;/* vertex/slice queries left for neighbour fill
                                  * rather than borrowing from another chart */
    size_t  uv_gauge_observations;/* same-chart continuity observations joining
                                   * disconnected Stage-C metric solves */
    size_t  uv_gauge_relations;/* robust component-pair relations retained */
    double  uv_gauge_max_shift;/* largest final additive gauge correction (vox) */
    size_t  uv_atlas_islands;/* reconstruction components with slice samples */
    size_t  uv_atlas_packed_islands;/* independent component intervals packed */
    double  uv_atlas_gutter; /* raster-safe inter-island material-U gutter (vox) */
    double  uv_atlas_pack_saved;/* cover-space span removed by metric packing */
    double  match_cover;     /* fraction of samples with a next-slice match */
    /* winding-scaffold assignment */
    int     w_groups;        /* distinct wraps found (maxW - minW + 1) */
    int     w_prior_groups;  /* chains placed by radius (== n_chains) */
    int     w_unreached;     /* chains where radius stops tracking winding (core folds) */
    size_t  w_conflicts;     /* unused (0); kept for JSON/ABI compatibility */
    double  pitch_used;      /* radial pitch used for the winding anchor (vox/turn) */
    int     pitch_source;    /* RIB_PITCH_* provenance for pitch_used */
    int     n_qp_comps;      /* connected components of the u-solve graph
                              * (chains + pairs); each is gauge-pinned and then
                              * REGISTERED onto the main chart's phi->u map */
    double  reg_max_shift;   /* largest registration shift applied (vox) */
    double  spiral_a, spiral_b, spiral_r2;  /* diagnostic fit r ~= a + b*phi/2pi */
    double  u_span, v_span;  /* extents (vox) */
    double  phi_span_turns;  /* winding extent of the samples / 2*pi */
    double  duds_err_mean;   /* mean |du/ds - 1| over chain edges (post-solve) */
    double  duds_err_max;    /* max  |du/ds - 1| */
    long    duds_hist[6];    /* |du/ds-1| in [0,1)% [1,2)% [2,5)% [5,10)% [10,20)% [20,..)% */
    size_t  grid_claim_conflicts; /* extra slice chains claiming an occupied
                                   * fitted-grid cell */
    size_t  grid_claim_conflict_cells; /* fitted-grid cells with >1 claimant */
    size_t  grid_claim_discarded_cells; /* conflicted cells deliberately made
                                         * unsupported constraints */
    size_t  grid_claim_replaced;  /* conflicts where component-global winding
                                   * consensus selected the later candidate */
    size_t  grid_claim_candidates; /* all interpolated claims, counted once */
    size_t  grid_claim_required_peels; /* maximum pre-selection cell depth */
    size_t  grid_claim_overflow;   /* claims beyond requested total capacity */
    size_t  grid_claim_stored;     /* measured claims selected for resident grids */
    size_t  grid_subcell_claim_chains;/* slice chains narrower than grid_du
                                       * retained by one nearest-column claim */
    size_t  grid_reconstruction_components;/* independently packed output
                                             * components after branch split */
    size_t  grid_branch_conflicts;/* same-row logical path pairs which overlap
                                    * in U but exceed the physical wrap gate */
    size_t  grid_branch_splits;   /* supported continuations forced to start a
                                    * new output lane by such a conflict */
    size_t  grid_branch_relation_cuts;/* robust chain-continuation relations
                                       * separated at ambiguous branches */
    double  grid_branch_relation_cut_support;/* summed support of those cuts */
    size_t  grid_long_edges_rowfit;/* defined horizontal/vertical grid edges
                                    * above the six-voxel physical gate */
    size_t  grid_claim_repairs;    /* monotone post-fit claimant substitutions */
    size_t  grid_long_edges_repaired;/* long edges after claimant repair */
    size_t  grid_row_wrap_splits;  /* adjacent-U claim pairs split by 3-D gate */
    size_t  grid_row_wrap_same_mesh;/* those pairs from one source component */
    size_t  grid_row_wrap_same_solve;/* those pairs in one metric solve */
    size_t  grid_row_wrap_same_island;/* those pairs in one winding island */
    /* hard chain merges + deterministic fusion (component_global mode) */
    size_t  grid_merge_junctions;  /* accepted certified chain continuations */
    size_t  grid_merge_ambiguous;  /* near-tied junction alternatives; nonzero
                                    * is the trigger for carrying explicit
                                    * port intervals in the layout, never a
                                    * license to guess */
    size_t  grid_merge_rejected_island;  /* junction candidates per gate */
    size_t  grid_merge_rejected_relation;
    size_t  grid_merge_rejected_du;
    size_t  grid_merge_rejected_gap;
    size_t  grid_merge_rejected_radial;
    size_t  grid_merge_rejected_phase;
    size_t  grid_fuse_same_logical;/* losing claims on the winner's own
                                    * logical chain (junction overlap) */
    size_t  grid_fuse_alternate;   /* losing claims within the wrap gate:
                                    * alternate measurements of one sheet */
    size_t  grid_fuse_violations;  /* losing claims beyond the wrap gate --
                                    * reported certificate violations */
    double  grid_fuse_violation_max;/* worst violating same-cell step (vox) */
    size_t  grid_transition_gap;   /* audited ownership transitions across an
                                    * unsupported gap (end of support) */
    size_t  grid_transition_certified;/* contiguous transitions on a direct
                                       * serialized chart relation */
    size_t  grid_transition_uncertified;/* contiguous transitions without a
                                         * relation edge -- reported */
    size_t  grid_long_edges_vfill; /* same count after bounded V hole fill */
    size_t  grid_local_outlier_repairs;/* isolated fitted samples replaced by
                                        * a certified coherent neighbour ring */
    size_t  grid_local_pair_repairs;/* adjacent two-sample outlier patches */
    size_t  grid_local_pair_slots;  /* vertices replaced by those patches */
    size_t  grid_local_run_repairs; /* one-cell-wide harmonic strip patches */
    size_t  grid_local_run_slots;   /* vertices replaced by strip patches */
    size_t  grid_local_patch_repairs;/* enclosed short-edge grid islands */
    size_t  grid_local_patch_slots;  /* vertices replaced by island patches */
    size_t  grid_local_patch_small_components;/* short-edge components at or
                                               * below the conservative cap */
    size_t  grid_local_patch_enclosed_components;/* small components whose
                                                  * complete exterior belongs
                                                  * to one large sheet */
    size_t  grid_local_patch_claim_slots;/* repaired island vertices snapped
                                          * back to measured slice claims */
    size_t  grid_local_supported_replacements;/* repaired slots originally
                                               * supported by slice data */
    size_t  grid_long_edges_local; /* count after local outlier repair */
    size_t  grid_long_edges_smooth;/* same count after guarded V smoothing */
    size_t  grid_both_diagonals_long;/* four-corner cells whose two possible
                                      * diagonals both exceed the gate */
    /* --- emitted-topology census (2026-09-01) -----------------------------
     * The counters above measure only the (k,j+1) and (k+1,j) lattice
     * neighbours.  The writer's fixed-diagonal quad also emits the diagonal
     * (k,j+1)-(k+1,j) in EVERY four-corner cell, and that edge had never been
     * measured: on the shipped 4x5x5 ribbon it carries 35,430 of the 70,631
     * long edges, so the fitted-grid line under-reported the damage by 2x. */
    size_t  grid_long_edges_diag_main;/* emitted b-c diagonal above the gate */
    size_t  grid_long_edges_diag_anti;/* the alternate a-d diagonal, measured
                                       * so an emitter that may flip a cell has
                                       * both numbers before it chooses */
    size_t  grid_cell_lane_hist[8];  /* distinct reconstruction lanes claiming
                                      * one lattice cell, bucketed 1..8+.  With
                                      * metric-island packing disabled for a
                                      * carried certificate, lanes legitimately
                                      * share cells; this is how many. */
    size_t  grid_vfill_runs;         /* bounded V continuations actually made */
    size_t  grid_vfill_rows;         /* lattice rows those continuations
                                      * generated (never measured, never
                                      * capped: the u-fill has RIB_UFILL_MAX,
                                      * the v-fill had no row bound at all) */
    /* --- lane-turn MRF (see lane_turn_mrf.h) --- */
    size_t  grid_lane_mrf_sites;     /* reconstruction lanes = MRF sites */
    size_t  grid_lane_mrf_edges;     /* accepted order + continuation factors */
    size_t  grid_lane_mrf_changed;   /* labels the solver changed */
    size_t  grid_lane_mrf_moved;     /* lanes left on a nonzero turn offset */
    size_t  grid_lane_mrf_chains;    /* chains those lanes carried */
    double  grid_lane_mrf_energy_before;
    double  grid_lane_mrf_energy_after;
    size_t  grid_lane_pair_observations;/* raw co-claim observations seen */
    size_t  grid_lane_pairs;         /* distinct lane pairs they fell into */
    size_t  grid_lane_pairs_dropped; /* observations the table could not hold */
    size_t  grid_emitted_faces;      /* faces the cell codes actually emit */
    size_t  grid_material_faces_cut; /* otherwise-legal triangles suppressed
                                      * because two known corners belong to
                                      * different material lineages */
    size_t  grid_emitted_long_edges; /* emitted edges above the wrap gate.
                                      * POSTCONDITION: must be 0.  Nonzero
                                      * means the emitter and the gate
                                      * disagree, which is a bug, not a
                                      * tuning outcome. */
    size_t  grid_cells_full;         /* four-corner cells considered */
    size_t  grid_cells_emit2;        /* ... emitting a complete quad */
    size_t  grid_cells_emit1;        /* ... emitting one triangle */
    size_t  grid_cells_emit0;        /* ... emitting nothing (the cut) */
    size_t  grid_cells_flipped;      /* ... where ANTI beat MAIN */
    size_t  grid_edge_cuts_h;        /* lattice edges cut horizontally */
    size_t  grid_edge_cuts_v;        /* ... and vertically */
    size_t  grid_ufill_reject_stretch;/* u chords refused: too long for the
                                       * arclength they replace */
    size_t  grid_ufill_reject_steps; /* ... refused: hole wider than
                                      * RIB_UFILL_MAX columns */
    size_t  grid_ufill_reject_crossing; /* invented chord crosses measured slice */
    size_t  grid_vfill_reject_rows;  /* v chords refused: run over
                                      * RIB_VFILL_MAX_ROWS */
    size_t  grid_vfill_reject_stretch;/* ... refused by the same stretch rule */
    size_t  grid_vfill_reject_crossing;/* invented chord stabs an input triangle */
    size_t  grid_vfill_run_hist[6];  /* their lengths in rows, bucketed
                                      * [1,2) [2,4) [4,8) [8,16) [16,64)
                                      * [64,inf) -- a deterministic stand-in
                                      * for a percentile */
    double  chain_gauge_pair_median_before;/* absolute cross-row constraint
                                             * residual before anchored solve */
    double  chain_gauge_pair_p95_before;
    double  chain_gauge_pair_median_after;
    double  chain_gauge_pair_p95_after;
    double  chain_gauge_shift_rms;/* sample-observation-weighted displacement
                                    * from the graph-derived chain gauge */
    double  chain_gauge_shift_p95;
    double  chain_gauge_shift_max;
} RibbonResult;

/* Run the full pipeline. verts [nv*3] (z,y,x float), faces [nf*3] 0-based.
 * All outputs are allocated from arena. Returns 0 on success, -1 on
 * degenerate input (too few verts/faces, or no plane crossings). */
int Ribbon_run(Arena_T arena,
               const float *verts, size_t nv,
               const int32_t *faces, size_t nf,
               const RibbonOpts *opts, RibbonResult *out);

typedef struct {
    size_t vertices, faces;
    size_t supported_vertices, generated_vertices;
    size_t atlas_columns, atlas_runs, empty_columns_removed;
    size_t primary_peel_vertices, peel_promotion_threshold;
    size_t promoted_peel_runs, promoted_peel_vertices;
    size_t extra_peel_runs, extra_peel_vertices;
} RibbonWriteStats;

/* Write the fitted regular grid: the authoritative UV'd VMESH companion of
 * `path` always, plus (when write_obj != 0) the compact one-vt-per-vertex text
 * OBJ at `path` itself for viewers.  The parameterizer's RIB_CELL_* decision is
 * serialized verbatim; the writer applies no geometric or identity policy.
 * Empty U-column runs are atlas-compacted.  v_rebase != 0 drops each
 * occupied column run to its own v minimum (atlas height = tallest run, not
 * the global z span -- the 21x21x21 atlas was 80% reserved-but-empty rows)
 * and records u-range -> v_offset in the <path base>_runs.json sidecar so
 * world z stays recoverable; note vt.v is then no longer globally z-zmin. */
int Ribbon_write_obj(const char *path, const RibbonResult *result,
                     int write_obj, int v_rebase, RibbonWriteStats *stats);

/* As above, but publish substantial collision-peel runs beside the primary
 * sheet in U instead of dropping every non-winning claim from the deliverable.
 * A peel run is a maximal contiguous set of occupied source-grid U columns.
 * Runs which satisfy either positive promotion criterion are promoted as
 * disjoint horizontal atlas panels: at least
 * `min_primary_share * primary_vertices`, or at least `min_vertices`.
 * Equivalently, the smaller positive threshold wins.  This lets a relative
 * rule scale down on small inputs without making a physically substantial run
 * harder to publish merely because a larger block contains more total sheet.
 * Smaller runs remain in the sibling <path base>_extras VMESH.  No geometry,
 * face, confidence, phase, or identity value is changed.  When both criteria
 * are zero this is exactly the legacy Ribbon_write_obj geometry policy.
 * Projective multi-layer output also atomically publishes <path base>_layers.json
 * with actual surface counts and explicitly non-surface/point-only layers.
 * Its previous completion record is invalidated before replacing any layer. */
int Ribbon_write_obj_promote_peels(const char *path,
                                   const RibbonResult *result,
                                   int write_obj, int v_rebase,
                                   double min_primary_share,
                                   size_t min_vertices,
                                   RibbonWriteStats *stats);

/* Deterministic in-memory coverage for the peel atlas partition/repack. */
int RibbonGridIO_selftest(void);

/* Flag faces the parameterization reveals as WRONG inter-wrap links. A genuine
 * bad weld link is PHYSICALLY LONG: it bridges the >=7-vox inter-wrap gap, so a
 * real link has a 3D edge >= len_min. A winding-collapse artifact is physically
 * SHORT (a normal ~2-vox within-sheet edge) that merely got a large |du| where
 * the core parameterization is unreliable -- those must NOT be cut, or we sever
 * real geometry and leave floaters. A face is flagged (out_bad[f]=1) only if
 * some edge satisfies ALL of:
 *   edge_3d_length >= len_min          (a real physical inter-wrap bridge)
 *   |u_a - u_b| > ratio * edge_length  (u stretched past arc length)
 *   |u_a - u_b| > floor_vox            (stretch is substantial)
 * Delamination (same sheet, |du| ~ 0) and collapsed-core short edges are kept.
 * verts [nv*3], faces [nf*3], uv [nv*2]. out_bad is caller-provided [nf].
 * Returns 0; writes the flagged count to *out_n. len_min <= 0 disables the
 * length gate (legacy du-only behavior). */
int Ribbon_flag_bad_faces(const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const float *uv, double ratio, double floor_vox,
                          double len_min, uint8_t *out_bad, size_t *out_n);

/* Map a (u,v) parameterization computed on a COARSE mesh onto the vertices of
 * the ORIGINAL mesh (simplify-first workflow, cf. successive
 * self-parameterization): closest-point projection onto the coarse surface,
 * barycentric UV interpolation, guarded by half the 7-vox inter-wrap
 * clearance so a projection can never land on a neighboring wrap. Unmapped
 * vertices (beyond the guard) are filled from mapped mesh neighbors (median);
 * *out_fallback counts vertices that stayed unmapped. Outputs arena-allocated.
 * Returns 0 on success, -1 on degenerate input. */
int Ribbon_map_uv(Arena_T arena,
                  const float *cverts, size_t cnv,
                  const int32_t *cfaces, size_t cnf,
                  const float *cuv,
                  const float *verts, size_t nv,
                  const int32_t *faces, size_t nf,
                  float **out_uv, uint8_t **out_ok, size_t *out_fallback);

/* In-process unit tests (chain builder, analytic spiral arc length, punched
 * holes, two-wrap separation, fusion wall, PAVA, coarse->fine map, degenerate
 * inputs). Returns 0 if all pass, else the number of failures. Logs stderr. */
int Ribbon_selftest(void);

#endif
