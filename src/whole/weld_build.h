#ifndef WELD_BUILD_INCLUDED
#define WELD_BUILD_INCLUDED

#include <stddef.h>
#include "../common/arena.h"
#include "cube_schedule.h"   /* CubeNode */
#include "placed_cube.h"     /* SkinVert */
#include "weld_track.h"

/* ============================================================================
 * weld_build.h -- Pass 1 of the tracked-weld dataset: cross-cube boundary
 * CORRESPONDENCES from the RAW skins. Same per-pair gating the forest edge
 * builder uses (3D distance, |dr| <= pitch/2, near-integer raw dphi), but emits
 * INDIVIDUAL vertex-pair correspondences classified by seam axis (z vs y/x),
 * and flags each pair's ambiguity: conf=2 when the nearest match is the only
 * wrap within a margin, conf=1 when another wrap sits nearly as close (the
 * mis-pairing risk the winding residual traces to).
 * ==========================================================================*/

typedef struct {
    double pair_gate;      /* 3D pairing gate, vox (3.5) */
    double radius_gate;    /* |dr| gate, vox (<=0 => 0.5*pitch) */
    double frac_gate;      /* max |wrap_pi(raw dphi)|, rad (0.5) */
    double margin_frac;    /* a different-wrap skin vert within
                              (1+margin_frac)*d_nearest marks the pair ambiguous
                              (conf=1); else conf=2 (0.5) */
    int    verbose;
} WeldBuildOpts;

void WeldBuildOpts_default(WeldBuildOpts *o);

/* Build Pass-1 correspondences. cnodes/skins/nskin are indexed alike (NULL/0
 * skin = cube absent); cnodes[c].nbr gives seam adjacency; axis_point (z,y,x)
 * for radii. Fills out: ids from cnodes[].id, spiral/axis/pitch, corr[]
 * (arena). Returns 0; -1 on bad input. */
int WeldBuild_pass1(Arena_T arena, const CubeNode *cnodes, size_t n_cubes,
                    SkinVert *const *skins, const size_t *nskin,
                    const float axis_point[3],
                    double spiral_a, double spiral_b, double pitch,
                    const WeldBuildOpts *opts, WeldTrack *out);

int WeldBuild_selftest(void);

#endif
