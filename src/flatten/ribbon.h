#ifndef RIBBON_INCLUDED
#define RIBBON_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

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
 *   W. assign each chain an INTEGER winding index directly from radius vs. the
 *      umbilicus: W = round(median(r/pitch - theta/2pi)) per chain, so a
 *      scroll's wrap count equals its maximal winding (no vote graph to
 *      fragment and collapse). phi* = theta + 2pi*W is the diffeomorphic
 *      spiral coordinate; delaminations at a parent's radius share its W.
 *   C. jointly solve for u per sample: arc-length preservation along each
 *      polyline + isovalue alignment across slices (quadratic energy, CG in
 *      double, per-chain monotone repair, likelihood reweighting);
 *   D. transfer (u, v = axial) to every mesh vertex via the 3D-nearest
 *      slice sample (winding-safe by the distance gate);
 *   E. fit the clean ribbon: a regular (u, v) grid of 3D positions smoothed
 *      per row (tangent then position tridiagonal solves) and filled/
 *      regularized across rows.
 *
 * u is TRUE per-slice arc length (vox), unlike unwrap.c's winding*r_ref
 * approximation. v is the axial coordinate (vox). Both start at 0.
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
    int    relax_iters;    /* likelihood-reweighting solve rounds (default 5) */
    int    final_iters;    /* tightened (sigma/3) polish rounds (default 2) */
    float  grid_u;         /* ribbon grid spacing along u, vox (default 2.0) */
    int    fit_ribbon;     /* nonzero: produce the ribbon grid (default 1) */
    float  wrap_spacing;   /* radial pitch (vox/turn) for the radius-anchored
                            * winding index. <=0 = auto-estimate from per-wedge
                            * adjacent-wrap radial gaps (default 0). */

    /* --- per-cube global-frame emission (scroll_whole; all default OFF, in
     * which case behavior is byte-identical to the classic single-mesh run) --- */
    int    emit_global;    /* nonzero: do NOT shift (u,v) to start at 0 -- emit
                            * u in the solve's absolute frame and v = raw axial
                            * coordinate t (world z for axis_point z=0, axis Z),
                            * fill RibbonResult.phi, and set u_origin. The
                            * ribbon grid (if fit) stays origin-based. */
    int    winding_sense;  /* 0 = auto (cov(th,r) sign). +1/-1 pins the radius-
                            * anchor sense so every cube of one scroll agrees. */
    int    pin_orient;     /* nonzero: skip the canonical u-flip (u increases
                            * outward). The flip is per-mesh (cov(u,r) sign), so
                            * cubes of one scroll would mirror independently. */
    double spiral_a, spiral_b;  /* spiral_b != 0: pin the global spiral line
                            * r ~= a + b*phi/2pi instead of fitting it, forcing
                            * every cube through the SAME arc-length init
                            * u(phi) = a*phi + b*phi^2/(4pi) (and disabling the
                            * per-component registration fallback). Calibrate
                            * (a,b) once on a seed cube, pin everywhere. */
    const float *reference_phi; /* optional [nv] graph-traced winding field.
                            * When supplied, slice crossings interpolate this
                            * authoritative absolute phi instead of rebuilding
                            * turns from cross-slice nearest-neighbour votes.
                            * StrokeStrip still solves true arc length; this is
                            * only its winding/turn scaffold. The array must
                            * stay alive for the duration of Ribbon_run. */
} RibbonOpts;

/* Fill opts with the defaults above (axis = Z through the origin; callers
 * must set axis_point/axis_dir for real scrolls). */
void RibbonOpts_default(RibbonOpts *opts);

enum {
    RIB_PITCH_FALLBACK  = 0, /* estimator had insufficient support */
    RIB_PITCH_ESTIMATED = 1, /* measured from adjacent-wrap radial gaps */
    RIB_PITCH_PINNED    = 2  /* caller supplied opts->wrap_spacing */
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
    float   *grid_pos;   /* [nk*nu*3] positions (z,y,x); NAN where absent */
    uint8_t *grid_valid; /* [nk*nu]: 1 = supported by slice data (pre-fill) */
    size_t   nu, nk;     /* grid columns (u) and rows (v = slices) */
    float    grid_du;    /* u spacing (== opts->grid_u) */
    float    grid_dv;    /* v spacing (== opts->slice_h) */

    /* --- diagnostics --- */
    int     n_slices;        /* slicing planes with any crossing */
    int     n_chains;        /* total polylines over all slices (post bridge cut) */
    int     n_closed;        /* closed loops (cut open) */
    int     n_multi_slices;  /* slices with >1 chain (fragmented) */
    size_t  n_samples;       /* resampled points (QP variables) */
    size_t  n_pairs;         /* cross-slice alignment pairs */
    size_t  n_cont_pairs;    /* intra-slice fragment continuation pairs */
    size_t  bridge_cuts;     /* chain splits at fusion-bridge radial jumps */
    size_t  mono_repairs;    /* samples moved by the monotone (PAVA) repair */
    size_t  uv_filled;       /* no direct sample; u copied from mesh neighbours */
    size_t  uv_fallback;     /* verts still unmapped after neighbour fill */
    double  match_cover;     /* fraction of samples with a next-slice match */
    /* winding-index assignment (radius-anchored per chain) */
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
} RibbonResult;

/* Run the full pipeline. verts [nv*3] (z,y,x float), faces [nf*3] 0-based.
 * All outputs are allocated from arena. Returns 0 on success, -1 on
 * degenerate input (too few verts/faces, or no plane crossings). */
int Ribbon_run(Arena_T arena,
               const float *verts, size_t nv,
               const int32_t *faces, size_t nf,
               const RibbonOpts *opts, RibbonResult *out);

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
