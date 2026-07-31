#ifndef ATLAS_SHEET_SPLIT_INCLUDED
#define ATLAS_SHEET_SPLIT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_radial_order.h"

/* ============================================================================
 * atlas_sheet_split.h -- how many sheets of papyrus is this piece of atlas,
 * and which sheet is each fragment on?
 *
 * The radial-order stage already answers the local question: for two samples
 * that sit on the same ray, phi says how many wraps apart they are.  That is a
 * DIFFERENCE.  What the atlas actually needs is a LEVEL -- an integer sheet
 * index per fragment -- because only then can a fragment be moved by a whole
 * number of turns without disturbing anything it is connected to.
 *
 * Turning differences into levels is a graph problem, and the graph is over
 * MESH COMPONENTS, not samples.  A mesh component is a physically connected
 * piece of surface, so every one of its vertices is on the same sheet by
 * definition, and a level assigned to the component moves the whole piece
 * rigidly.  Assign levels per sample instead and a shift would tear the piece
 * apart along whatever line the labelling happened to change on.
 *
 * That "by definition" is a claim about the data, not a theorem: a mesh
 * component long enough to spiral past itself would carry two sheets at once.
 * The build checks it -- `pairs_intra_component` counts radial pairs whose two
 * ends land in the SAME component -- and reports it rather than assuming.  On
 * the 4x5x5 fixture's support 3 that count is zero.
 *
 * THE ISLAND PROBLEM.  Each connected piece of the turn graph fixes its levels
 * up to one additive constant, exactly as the strip solve fixes u up to one
 * constant per support component.  Radial pairs only exist between components
 * that share a ray, so a support component typically breaks into several
 * islands with no evidence between them.  Aligning those islands needs a
 * second, independent signal, and the honest one is radius: under the pinned
 * spiral a sheet sits one pitch further out than the sheet inside it, so an
 * island's radial intercept
 *
 *     c = median over its components of (mean_radius - pitch * level)
 *
 * places it against any other island's, and the offset that aligns them is
 * lround((c - c_ref) / pitch).  This is the one step where the scroll's ~85-vox
 * axis wander can bite, so the rounding residual is reported per island in
 * PITCHES -- a residual near 0.5 means the alignment is a coin flip and should
 * be disbelieved.  phi is carried alongside as a fully independent second
 * opinion and its agreement rate is reported; it is never used to decide,
 * because phi is registered per cube and islands routinely span cubes.
 * ==========================================================================*/

typedef struct {
    /* Component-pair edges carrying fewer than this many radial pairs are
     * dropped before the forest is grown.  A single stray pair between two
     * wraps is exactly what an accidental bridge looks like. */
    int    min_edge_pairs;
    /* Rounding residual above which an island alignment is called unreliable
     * and reported as such (in pitches; 0.5 is a coin flip). */
    double align_max_residual;
} AtlasSheetSplitOptions;

void AtlasSheetSplitOptions_default(AtlasSheetSplitOptions *opts);

typedef struct {
    int32_t mesh_component;
    int32_t island;       /* connected piece of the turn graph */
    int32_t local_level;  /* turn level within the island, arbitrary zero */
    int32_t sheet;        /* globally aligned, rebased so the innermost is 0 */
    size_t  samples;
    size_t  pairs;        /* radial pairs touching this component */
    size_t  pairs_internal; /* radial pairs with BOTH ends in this piece */
    double  mean_radius;
    double  mean_phi;
    double  phi_min;
    double  phi_max;
    double  radius_min;
    double  radius_max;
    /* (phi_max - phi_min) / 2pi.  A piece with a turn span near or above 1 is
     * a SPIRAL, and no single integer sheet index describes it -- its wrap
     * number changes along its own length. */
    double  turn_span;
} AtlasSheetComponent;

typedef struct {
    int32_t island;
    size_t  components;
    size_t  samples;
    int32_t level_min;
    int32_t level_max;
    double  intercept;        /* median(mean_radius - pitch*local_level) */
    int32_t offset;           /* whole-turn shift applied to align it */
    double  residual;         /* alignment rounding residual, in pitches */
    double  mean_radius;
} AtlasSheetIsland;

typedef struct {
    size_t components;            /* mesh components in the target support */
    size_t components_with_pairs;
    size_t islands;
    size_t singleton_islands;     /* components with no radial evidence at all */

    size_t pairs_used;
    size_t pairs_intra_component; /* a component spanning two wraps: see above */
    size_t edges;                 /* distinct ordered component pairs */
    size_t edges_dropped_weak;    /* below min_edge_pairs */
    size_t edges_in_forest;
    size_t edges_consistent;      /* non-tree edges the levels already satisfy */
    size_t edges_contradictory;   /* non-tree edges they do not */

    int32_t sheet_min;
    int32_t sheet_max;
    size_t  sheets;               /* distinct occupied sheet indices */

    double align_residual_rms;    /* over islands, in pitches */
    double align_residual_max;
    size_t align_unreliable;      /* islands past align_max_residual */

    /* Independent cross-check.  Expected sheet from phi alone is
     * (phi - phi_ref)/2pi * sense; how often does it match? */
    size_t phi_agree;
    size_t phi_disagree;

    /* How badly the constant-per-piece model fits.  Once fragments are welded
     * into continuous surfaces these are the numbers that matter: a piece
     * spanning more than one turn cannot carry a single sheet index at all. */
    size_t spiral_components;      /* turn_span >= 1 */
    double turn_span_max;
    double turn_span_total;
} AtlasSheetSplitStats;

typedef struct {
    AtlasSheetComponent *components;
    size_t ncomponents;
    AtlasSheetIsland *islands;
    size_t nislands;
    /* sheet_of[c] for every mesh component id in [0, nmesh); -1 = not in the
     * target support component, or no level could be assigned. */
    int32_t *sheet_of;
    size_t nmesh;
    AtlasSheetSplitStats stats;
} AtlasSheetSplitResult;

/*
 * pairs/npairs         radial-order evidence, any support component
 * sample_component[i]  mesh component of sample i
 * sample_support[i]    support component of sample i
 * sample_radius[i]     distance from the scroll axis, vox
 * sample_phi[i]        absolute multi-turn winding angle, radians
 * target_support       restrict to this support component, or -1 for all
 * pitch, sense         from the pinned calibration
 */
int AtlasSheetSplit_solve(Arena_T arena,
                          const AtlasRadialOrderPair *pairs, size_t npairs,
                          const int32_t *sample_component,
                          const int32_t *sample_support,
                          const double *sample_radius,
                          const double *sample_phi,
                          size_t nsamples,
                          size_t nmesh_components,
                          int32_t target_support,
                          double pitch,
                          double sense,
                          const AtlasSheetSplitOptions *opts,
                          AtlasSheetSplitResult *out);

int AtlasSheetSplit_selftest(void);

#endif
