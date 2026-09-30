#ifndef QUADRIBBON_SOLIDIFY_INCLUDED
#define QUADRIBBON_SOLIDIFY_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "../whole/axis_warp.h"

/* quadribbon_solidify -- fit ONE solid quad ribbon to a fitted lattice ribbon.
 *
 * Input: any ribbon VMESH whose UV is a lattice (u = atlas column * du,
 * v = row * dv, rows = constant-z slices, optionally stacked in peel bands)
 * plus its sidecars (<stem>_support.u8 required; _phase.f32 required for the
 * winding-aware fill; _provenance.u8 / _material_identity.i32 / _lane.i32 or
 * _reconstruction_component.i32 / _stats.json optional).  The observations
 * are the vertices; cells carrying provenance QS_PROV_FILL are treated as
 * UNOBSERVED, which is what makes the same call serve as the refit after the
 * optimize stage moved the observed geometry.
 *
 * Output (<out_stem>.vmesh + sidecars): one lattice per contiguous atlas
 * column RUN (runs are separate wraps and are never bridged), every observed
 * cell exact, holes filled by quad_strip's cylindrical residual harmonic
 * solve (TAUCS), fills capped at QS_MAX_HOLE_CELLS from the nearest
 * observation and confined to the source box and outside the umbilicus core,
 * every vertex labeled by provenance.  Peel band 0 wins a cell; a band-1
 * observation is adopted only when it sits on the same wrap as its band-0
 * neighbourhood (radius within QS_BAND1_ADOPT_DR_PITCH pitches), otherwise it
 * is another wrap the certificate could not separate: counted, never planted.
 *
 * All linear algebra goes through sparse_solve.h (TAUCS); the caller's
 * executable must reserve the 512 MB stack TAUCS needs. */

enum {
    QS_PROV_OBSERVED_FIXED = 0,     /* input support 255 */
    QS_PROV_OBSERVED_GENERATED = 1, /* input support 0: the fit's own short
                                     * u/v continuation; preserved in projective
                                     * mode, but never used as direct evidence */
    QS_PROV_BAND1 = 2,              /* adopted peel-band observation */
    QS_PROV_FILL = 3,               /* harmonic fill (unobserved) */
    QS_PROV_OBSERVED_SNAPPED = 4    /* observed cell that baked dark and was moved
                                     * along its normal onto the CT ridge (<= 4 vox) */
};

typedef struct QsConfig {
    double axis_y, axis_x;      /* scroll axis point (world y, x) */
    const AxisWarp *coordinate_warp; /* shared checked metric -> source CT map */
    /* CT FRAME.  When the lane straightens the mesh about an axis table, every
     * stage after it -- this one included -- works in STRAIGHTENED coordinates,
     * while the raw CT is in world coordinates.  A probe that samples the CT at
     * a straightened position reads the wrong place (up to the axis drift, ~200
     * voxels on the 10x3x3).  These fields carry the same curve the lane
     * straightened with, so a probe can undo it:
     *     world = straightened + axis(z) - reference
     * Leave axis_curve_n at 0 when the lane is not straightened. */
    const double *axis_curve_z, *axis_curve_y, *axis_curve_x;
    size_t axis_curve_n;
    double axis_ref_y, axis_ref_x;
    double pitch;               /* wrap pitch, voxels */
    int    projective_local;    /* authoritative nested-domain mode: absolute
                                 * lattice addresses, direct observations
                                 * immutable, upstream continuation classified,
                                 * direct-only fixed-stencil non-recursive fills */
    int    projective_copy_only;/* authoritative atlas-member mode: Stage 4 is
                                 * a lossless semantic transport.  Every Stage-3
                                 * vertex/face and claimant identity is copied;
                                 * no storage-color-local fill is permitted. */
    double lattice_du, lattice_dv; /* required canonical lattice pitch in
                                    * projective mode; never inferred from a crop */
    int    local_fill_radius;   /* maximum missing sites bracketed by one
                                 * fixed-stencil fill; <=0 selects
                                 * QS_BRIDGE_MAX_CELLS */
    int    run_min_observed;    /* a column run with fewer observed cells is
                                 * confetti and is not emitted; <= 0 selects
                                 * QS_RUN_MIN_OBSERVED (the selftest lowers it
                                 * for its small synthetic lattice) */
    double good_column_frac;    /* <= 0 selects QS_GOOD_COLUMN_FRAC */
    int    fill_max_u_gap;      /* <= 0 selects QS_FILL_MAX_U_GAP */
    double skip_dr_pitch;       /* <= 0 selects QS_SKIP_DR_PITCH */
    int    run_merge_gap;       /* LAYOUT: unoccupied column gaps up to this wide stay
                                 * inside one run (0 = every gap splits; the lane
                                 * passes QS_RUN_MERGE_GAP_COLS) */
    const char *raw_zarr;       /* uncompressed u8 zarr v2 root of the RAW CT, or
                                 * NULL: with it, fill cells snap to the CT ridge
                                 * along the lattice normal (QS_FILL_SNAP_*) */
} QsConfig;

typedef struct QsReport {
    size_t nv_in, nf_in, nv_out, nf_out;
    int    n_runs, n_bands, n_slices, runs_built, runs_retried, runs_observed_only;
    double du, dv, slice_h;
    size_t observed, observed_fixed, observed_generated, input_fills;
    size_t band1_adopted, band1_rejected, band1_no_reference, duplicates;
    size_t double_cover;
    double dup_dr_pitch_p10, dup_dr_pitch_p50, dup_dr_pitch_p90;
    size_t filled, dropped_bbox, dropped_core, dropped_unreferenced;
    size_t holes, holes_refused, holes_filled;   /* hole-level fill votes */
    size_t snap_candidates, snap_moved, snap_no_ridge;   /* CT-ridge snap of fills */
    size_t fold_cells, fold_columns;            /* FOLD demotion */
    size_t dark_cells, dark_snapped, dark_unrecovered;   /* dark observed cells */
    size_t hole_cells, refused_hole_cells;      /* hole vote areas */
    size_t midline_probed, midline_outliers;    /* MIDLINE test: suspects (all) */
    size_t midline_isolated, midline_dense;     /* isolated suspects (refit) / dense (left, reported) */
    size_t scale_faces, scale_over;             /* SCALE gate: emitted triangles, and those
                                                 * whose |log2(area3d / area_uv)| exceeds
                                                 * QS_GATE_SCALE_LOG2 */
    size_t scale_over_fill;                     /* ... of which touch a fill cell */
    double scale_p50, scale_p90, scale_p99;     /* |log2| over all emitted triangles */
    double scale_obs_p90, scale_fill_p90;       /* ... split observed-only / fill-touching */
    size_t obs_smooth_cells;                     /* DEPTH SMOOTHING: cells moved */
    double obs_smooth_p50, obs_smooth_p90, obs_smooth_max;
    size_t column_probed, column_jumps;          /* COLUMN-CONTINUITY repair */
    double column_dev_p50, column_dev_p90;       /* |radius - column median|, pitches */
    size_t depth_probed, depth_outliers;         /* DEPTH-OUTLIER demotion */
    double depth_hp_p50, depth_hp_p90;
    size_t ridge_cells, ridge_moved;             /* RIDGE TRACK */
    size_t ridge_on_before, ridge_on_after;
    size_t ridge_fill_cells, ridge_fill_on_before, ridge_fill_on_after;   /* fills tracked with the wide reach */
    double ridge_shift_p50, ridge_shift_p90, ridge_shift_max;
    double ridge_ct_before, ridge_ct_after;
    size_t ramp_runs, ramp_runs_over;            /* RAMP gate: runs measured / over the limit */
    size_t ramp_cells, ramp_cells_over;          /* observed cells in them */
    double ramp_p50, ramp_p90, ramp_max;         /* |pitches of radius gained per turn| */
    size_t dropped_fill_far;                     /* FILL CAP: cells refused as unsupported */
    size_t fill_cells, fill_far;                /* FILL-SUPPORT: emitted fill cells and those
                                                 * further than QS_FILL_SUPPORT_CELLS from an
                                                 * observation */
    double fill_dist_p50, fill_dist_p90;
    size_t fill_dist_max;
    double midline_p50, midline_p90;            /* |ridge offset| of kept observed cells, vox */
    int    gates_failed;                        /* solid gates (QS_GATE_*) */
    double snap_shift_p50, snap_shift_p90;      /* |displacement| vox */
    long   cols_out;            /* emitted sheet width in columns after layout */
    size_t dropped_unbridged;   /* outside the first..last good column of the
                                 * run, or in a non-good stretch wider than
                                 * QS_FILL_MAX_U_GAP: never extrapolated in u */
    size_t dropped_skip;        /* in a gap whose two sides sit on different
                                 * wraps (radial extrapolation test) */
    size_t dropped_chord;       /* a row bridge through the cell fails the
                                 * stretch / absolute chord gate */
    size_t dropped_material;    /* a gap or a row bridge joins two lineages */
    size_t dropped_order;       /* a gap or a row bridge joins fragments placed
                                 * out of angular order (angle change vs span) */
    size_t dropped_confetti;    /* observed cells of runs below QS_RUN_MIN_OBSERVED */
    size_t faces_skipped_gate, faces_skipped_material;
    double phase_mismatch_max;
    size_t phase_mismatch_over_limit;
    double t_load, t_lattice, t_build, t_emit, t_write, t_total;
} QsReport;

/* Fit the solid ribbon of in_vmesh and write <out_stem>.vmesh + sidecars +
 * <out_stem>_provenance.png + <out_stem>_report.json.  Returns 0 on success. */
int QuadribbonSolidify_run(Arena_T arena, const QsConfig *cfg,
                           const char *in_vmesh, const char *out_stem,
                           QsReport *rep);

/* Resume predicate: the VMESH is complete AND every required sidecar has
 * exactly nv elements AND the report exists. */
int QuadribbonSolidify_complete(const char *out_stem);

/* Diagnostic access to the exact production CT-coordinate/sampling path.
 * All points must have valid inverse coordinates and available RAW samples.
 * Returns nonzero for a refused coordinate, missing chunk or volume boundary. */
int QuadribbonSolidify_probe_ct(const QsConfig *cfg, const double *metric_points,
                                size_t count, double *source_points, uint8_t *values);

int QuadribbonSolidify_selftest(void);

#endif
