/* ridge_track.h -- put a ribbon patch ON the papyrus, coherently.
 *
 * The fit places a lattice cell where its claims say, and the solid stage
 * interpolates between them.  Measured on PHerc0139 (2026-09-03): the result
 * sits p50 3-6 voxels from the nearest CT wrap ridge in a stack whose local
 * spacing is 15-17 voxels, so a long run drifts across the gap and onto the
 * next wrap.  That drift is what bakes as the reviewer's swirl / marble: a
 * surface running on a wrap's shoulder reads a smear of two wraps.
 *
 * The wraps are SEPARATED there (measured: 0% of samples had no dark gap
 * either side), so there is a ridge to sit on and a gap to stay out of.
 *
 * This module chooses, per lattice cell, a displacement along the cell's own
 * normal that
 *   - lands on a CT ridge (bright, locally maximal),
 *   - stays close to where the fit put it (a prior, so a patch cannot walk to
 *     a different wrap wholesale),
 *   - and varies smoothly over the lattice, so neighbouring cells choose the
 *     SAME ridge.  That last term is the whole point: an independent per-cell
 *     snap was measured worse (it tore the surface at every ridge ambiguity).
 *
 * The solve is a dynamic program along u per row over a discretised offset
 * axis, then a smoothing pass across v, iterated a few times.  It is O(cells x
 * offsets) and needs no linear algebra.
 */
#ifndef RIDGE_TRACK_H
#define RIDGE_TRACK_H

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/zarr_u8.h"

typedef struct RidgeTrackOpts {
    double reach;          /* half-width of the offset search, voxels */
    double step;           /* offset discretisation, voxels */
    double prior;          /* cost per voxel of |offset| (stay near the fit) */
    double smooth_u;       /* cost per voxel of offset change between columns */
    double smooth_v;       /* cost per voxel of offset change between rows */
    double ridge_min;      /* CT value below which a sample is not papyrus */
    int    iters;          /* row-DP / column-smoothing rounds */
    int    max_step;       /* hard bound on the offset change between neighbours, in steps */
    double edge_limit;     /* no displaced lattice edge may exceed this; 0 disables */
    double flatten;        /* cost per voxel of depth deviation from the local median */
    int    flatten_win;    /* that median's window, in cells */
    int    profile_win;    /* cells over which the CT profile is averaged first */
    double fill_reach;     /* half-width of the search for WIDE cells (fills): no observation to stay near */
    double fill_prior;     /* their cost per voxel of |offset| */
} RidgeTrackOpts;

typedef struct RidgeTrackReport {
    size_t cells;          /* cells with a normal and a CT sample */
    size_t moved;          /* cells whose |offset| >= step */
    size_t on_ridge;       /* cells whose final sample is >= ridge_min */
    size_t on_ridge_before;
    double shift_p50, shift_p90, shift_max;
    size_t edge_clamped;   /* cells whose offset was shrunk to keep an edge legal */
    double edge_max;       /* longest displaced lattice edge */
    double hp_p90_before, hp_p90_after;  /* |depth - local median| over the cells */
    double neighbour_max;  /* largest |offset difference| between 4-neighbours after the
                            * Lipschitz projection: the guarantee the tangle gates need */
    double ct_before_p50, ct_after_p50;
    size_t fill_cells, fill_on_before, fill_on_after;   /* the wide (fill) cells */
} RidgeTrackReport;

void RidgeTrack_defaults(RidgeTrackOpts *opts);

/* A sampler the caller supplies so the tracker never has to know about frames:
 * it receives a position in the LATTICE's frame and returns the CT there. */
typedef double (*RidgeTrackSampler)(void *ctx, double z, double y, double x);

/* Maps a point from the lattice's frame to the frame the VERDICT measures edge
 * lengths in.  When the lane straightens about an axis curve those differ, and a
 * displacement in z changes the straightening offset, so an edge that is legal
 * in the lattice frame can be too long in world.  NULL = the frames coincide. */
typedef void (*RidgeTrackToWorld)(void *ctx, const double in_zyx[3], double out_zyx[3]);

/* pos: [W*H*3] cell positions in the LATTICE's frame, (z,y,x); modified in place.
 * valid: [W*H] nonzero where a cell exists.  Cells without a computable normal
 * or without CT are left untouched.  Returns 0 on success. */
/* `tri` (ntri triangles of 3 cell indices) is the set of faces the caller will
 * actually emit; the edge clamp uses exactly those edges, so the tangle gates
 * are guaranteed rather than approximated.  Pass NULL to fall back to the
 * lattice's own 8-neighbourhood. */
/* `depth` (one value per cell, e.g. the radius about the scroll axis) drives the
 * FLATTEN term; NULL disables it.  Moving a cell by t along its normal is taken
 * to change its depth by t, which holds while the normal is roughly radial. */
int RidgeTrack_solve(Arena_T arena, const RidgeTrackOpts *opts,
                     RidgeTrackSampler sample, RidgeTrackToWorld to_world, void *sample_ctx,
                     float *pos, const uint8_t *valid, const uint8_t *wide, int W, int H,
                     const int32_t *tri, size_t ntri, const double *depth,
                     RidgeTrackReport *rep);

int RidgeTrack_selftest(void);

#endif /* RIDGE_TRACK_H */
