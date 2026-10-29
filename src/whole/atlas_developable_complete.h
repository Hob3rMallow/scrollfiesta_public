#ifndef ATLAS_DEVELOPABLE_COMPLETE_INCLUDED
#define ATLAS_DEVELOPABLE_COMPLETE_INCLUDED

#include <stddef.h>

#include "../common/arena.h"
#include "atlas_ribbon_fit.h"
#include "atlas_track_grow.h"

/*
 * Complete bounded holes between trusted evidence-grown banks.
 *
 * This stage deliberately consumes the in-memory result built from an
 * AtlasSolution checkpoint.  It does not reconstruct state from an OBJ.
 * Existing vertices and faces are copied exactly.  A completion patch is
 * normally admitted when the same ordered bank pair supplies coherent support
 * across several consecutive atlas samples.  Horizontal patches bridge
 * missing U within evidence rows; an opt-in orthogonal family bridges missing
 * V across evidence rows.  The aligned-gap option additionally closes exact
 * cells and short, fully bracketed strips by atlas location even when source
 * component IDs change at a cube boundary.  Longer missing samples are seeded
 * from observed bank positions and tangents, then regularized as a patch-local
 * scalar height field with a Hessian rank-tail ADMM (the smaller principal
 * curvature is penalized).
 */

typedef struct {
    double min_gap_u;             /* smallest unsupported U span (4) */
    double max_gap_u;             /* largest unsupported U span (96) */
    double min_gap_v;             /* smallest unsupported V span (4) */
    double max_gap_v;             /* largest unsupported V span (96) */
    int min_patch_rows;           /* consecutive bank rows required (4) */
    int max_bank_shift_columns;   /* allowed boundary drift per row (6) */
    int sample_columns;           /* 0 = derive from observation spacing */
    int fill_internal_gaps;       /* include same-component holes (0) */
    int fill_vertical_gaps;       /* include gaps across evidence rows (0) */
    int allow_component_drift;    /* track banks across fragment IDs (0) */
    int shared_lattice_fill;      /* union validated U/V patches on grid (0) */
    int fill_aligned_gaps;        /* close aligned cells/short strips (0) */
    int short_strip_max_span_u;   /* max bracketed missing-U span (8 cells) */
    int short_strip_max_span_v;   /* max bracketed missing-V span (8 cells) */
    double tangent_dot_min;       /* endpoint tangent compatibility (.2) */
    double normal_cone_dot_min;   /* local graph normal hemisphere gate (.25) */
    double chord_min_fraction;    /* chord / U lower gate (.15) */
    double chord_max_stretch;     /* chord / U upper gate (1.8) */
    double chord_slack;           /* absolute chord allowance (4) */
    double vertical_stretch;      /* bank row-edge / dv upper gate (2) */
    double vertical_slack;        /* absolute row-edge allowance (4) */
    double bank_fit_rms_max;      /* affine bank residual gate (2) */
    double bank_fit_max_max;      /* affine bank maximum gate (6) */
    double hermite_prior_weight;  /* interior evidence-tangent prior (.2) */
    double boundary_weight;       /* exact-bank fidelity weight (64) */
    double nuclear_lambda;        /* convex Hessian shrink phase (.05) */
    double rank_tail_lambda;      /* smaller-curvature shrink phase (.2) */
    double admm_rho;              /* augmented-Lagrangian penalty (1) */
    int nuclear_iterations;       /* convex warm-start iterations (8) */
    int rank_tail_iterations;     /* rank-one refinement iterations (24) */
    int pcg_iterations;           /* h-subproblem iterations (160) */
    double pcg_tolerance;         /* relative residual tolerance (1e-7) */
    int max_patches;              /* deterministic safety cap (4096) */
} AtlasDevelopableCompleteOptions;

typedef struct {
    int32_t id;
    int direction;                /* 0 = missing U, 1 = missing V */
    int32_t component_left;
    int32_t component_right;
    int32_t row_first;
    int32_t row_last;
    int32_t column_left_min;
    int32_t column_left_max;
    int32_t column_right_min;
    int32_t column_right_max;
    size_t rows;
    size_t columns;
    size_t added_vertices;
    size_t added_faces;
    size_t component_pair_changes;
    double gap_u_mean;
    double gap_u_max;
    double chord_ratio_mean;
    double bank_fit_rms;
    double bank_fit_max;
    double tail_before;
    double tail_after;
    double primal_residual;
    double dual_residual;
    int pcg_iterations;
    int accepted;
    int reject_reason;
} AtlasDevelopablePatch;

typedef enum {
    ATLAS_DEVELOPABLE_ACCEPTED = 0,
    ATLAS_DEVELOPABLE_REJECT_SHORT = 1,
    ATLAS_DEVELOPABLE_REJECT_BANK_DRIFT = 2,
    ATLAS_DEVELOPABLE_REJECT_BANK_METRIC = 3,
    ATLAS_DEVELOPABLE_REJECT_PARAMETER = 4,
    ATLAS_DEVELOPABLE_REJECT_GEOMETRY = 5,
    ATLAS_DEVELOPABLE_REJECT_OVERLAP = 6
} AtlasDevelopableReject;

typedef struct {
    AtlasTrackGrowResult mesh;    /* seed copied, accepted patches appended */

    AtlasDevelopablePatch *patch;
    size_t npatches;
    size_t accepted_patches;
    size_t rejected_patches;
    size_t candidate_rows;
    size_t candidate_pairs;
    size_t seed_components;
    size_t completed_components;
    size_t seed_vertices;
    size_t seed_faces;
    size_t added_vertices;
    size_t added_faces;
    size_t lattice_nodes;
    size_t lattice_cells;
    size_t lattice_conflicts;
    size_t lattice_rejected_cells;
    size_t aligned_gap_lines;
    size_t aligned_gap_cells;
    size_t adjacent_gap_cells;
    size_t adjacent_gap_u_cells;
    size_t adjacent_gap_v_cells;
    size_t short_strip_patches;
    size_t short_strip_cells;
    size_t short_strip_u_cells;
    size_t short_strip_v_cells;
    size_t admm_iterations;
    size_t pcg_iterations;
    double tail_before;
    double tail_after;
    double primal_residual_max;
    double dual_residual_max;
} AtlasDevelopableCompleteResult;

void AtlasDevelopableCompleteOptions_default(
    AtlasDevelopableCompleteOptions *opts);

int AtlasDevelopableComplete_build(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteOptions *opts,
    AtlasDevelopableCompleteResult *out);

int AtlasDevelopableComplete_selftest(void);

#endif
