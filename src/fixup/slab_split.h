#ifndef SLAB_SPLIT_INCLUDED
#define SLAB_SPLIT_INCLUDED

#include <stdint.h>
#include "../common/arena.h"

/* slab_split — geometric mid-plane carver for pred-fused wrap slabs.
 *
 * The prediction paints radially-stacked wraps in tight contact as ONE slab
 * (measured 2026-08-26: 11,474 merged radial runs on 10x10x10, width ratio
 * p90 = 4.25 — up to four wraps as one body).  Downstream the per-cube
 * extractor lifts each 3D connected component as one point cloud, so a fused
 * slab meshes as one sheet and the ribbon fit inherits the wrong winding
 * (marble).  Raw-CT seams inside the slabs are unreliable (only ~20% show a
 * dark interface; the median interior dip is ~25%, inside fiber-texture
 * noise), so the carve is GEOMETRIC: wraps in a slab have near-equal
 * thickness, so the interface lies on the slab's medial surface.
 *
 * Per z-plane: 3-4 chamfer half-width -> primary ridge candidates (half-width
 * >= thick_r AND STRICT axis local max: >= both axis neighbours and > at
 * least one, so a constant-along-tangent interior never floods) and
 * extension candidates (same ridge test at >= ext_r, the hysteresis tier
 * that follows the interface into the fork tips where the wraps genuinely
 * diverge).  A primary survives only with z-support (a primary within
 * Chebyshev 2 on >= z_support of the +/-z_window neighbour planes) — slabs
 * are tall, plane noise is not.  Surviving primaries grow through candidates
 * (8-conn BFS, depth-capped at grow_max so the carve can bridge a fork tip
 * but can never run away down a normal sheet's medial axis) and the grown
 * set is carved to background, capped per plane at plane_cap_frac of that
 * plane's foreground (a capped plane is left UNCARVED and reported).
 * rounds > 1 re-runs the whole pass so a k-wrap slab (k >= 3) sheds one
 * interface per round.
 *
 * The carve deliberately removes interior voxels only (a ridge pixel of a
 * thick slab is >= thick_r from any background), so normal-thickness sheet
 * skins can never be breached — tangential continuity of each wrap survives
 * by construction. */

typedef struct SlabSplitParams {
    float thick_r;        /* primary ridge half-width floor (px)      [2.25] */
    float ext_r;          /* hysteresis extension floor (px)          [1.60] */
    int   grow_max;       /* hysteresis BFS depth cap (px)            [8]    */
    int   rounds;         /* full-pass repeats (k-wrap slabs)         [2]    */
    int   z_window;       /* support window half-extent (planes)      [3]    */
    int   z_support;      /* min supporting planes in the window      [2]    */
    float plane_cap_frac; /* max carved/foreground per plane          [0.10] */
} SlabSplitParams;

typedef struct SlabSplitStats {
    int64_t primary_px;   /* primary ridge candidates (pre-support), all rounds */
    int64_t supported_px; /* primaries surviving z-support, all rounds */
    int64_t carved_px;    /* voxels actually cleared */
    int64_t planes_carved;  /* planes with >= 1 carved px (round 1) */
    int64_t planes_capped;  /* planes skipped by the cap, all rounds */
    int     rounds_run;     /* rounds that carved > 0 */
} SlabSplitStats;

void SlabSplit_params_default(SlabSplitParams *p);

/* Carve vol (0/1, [D*H*W]) in place.  carve_mask (may be NULL) gets 1 at
 * every cleared voxel (caller-owned, [D*H*W], NOT cleared first).  Returns 0
 * on success. */
int SlabSplit_run(Arena_T arena, uint8_t *vol, uint8_t *carve_mask,
                  int D, int H, int W,
                  const SlabSplitParams *p, SlabSplitStats *st);

int SlabSplit_selftest(void);

#endif
