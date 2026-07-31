#ifndef NECK_PROBE_INCLUDED
#define NECK_PROBE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * neck_probe.h -- local "side-facing triangle" detector for bridge diagnosis.
 *
 * Targets the failure the UV-difference oracle is blind to: BPA (or the seam
 * weld) joining two scroll wraps with a thin triangle whose normal faces
 * side-to-side (tangential to both sheets, perpendicular to the local stacking
 * direction). Such a triangle is a strong outlier in orientation relative to
 * its edge-adjacent neighbours.
 *
 * Purely geometric and purely local: for every triangle we compare its face
 * normal to the (area-weighted) mean normal of its edge-adjacent faces and flag
 * the triangle when the acute deviation exceeds a threshold. This is winding-
 * sign independent (folded to [0,90] deg, as in obj_components.c's crease test).
 *
 * A side-facing flag alone is NOT proof of a bridge -- a legitimate sharp crease
 * also flags. It is meant to be CROSS-CHECKED against the topological neck
 * (reeb.h) and per-component genus (mesh_topo.h): a real bridge is side-facing
 * AND sits on a thin topological neck. See the reeb_probe tool.
 * ==========================================================================*/

typedef struct {
    float deg_thresh;     /* flag when neighbourhood deviation >= this (deg).
                           * Default 60. Higher = stricter (fewer flags). */
    int   min_neighbours; /* require >= this many edge-adjacent faces before
                           * judging a fold (default 2). 0/1-neighbour tips are
                           * reported via tip_flag instead. */
} NeckProbeOpts;

typedef struct {
    size_t   nf;            /* echoes the input face count */
    size_t   n_sidefacing;  /* #triangles flagged side-facing */
    size_t   n_tips;        /* #triangles with < min_neighbours support (dangling) */
    float    max_dev_deg;   /* largest neighbourhood deviation observed */
    uint8_t *face_flag;     /* [nf]: 1 = side-facing candidate (arena-allocated) */
    float   *face_dev;      /* [nf]: per-face deviation in deg (arena-allocated) */
} NeckProbeResult;

/* Defaults: deg_thresh = 60, min_neighbours = 2. Pass NULL for opts to use them. */
void NeckProbe_default_opts(NeckProbeOpts *opts);

/* Scan a triangle mesh. verts is [nv*3] (z,y,x); faces is [nf*3] 0-based.
 * All result arrays are arena-allocated. Returns 0 on success, -1 on empty
 * input. Handles zero verts / zero faces gracefully (returns -1, out zeroed). */
int NeckProbe_scan(Arena_T arena,
                   const float *verts, size_t nv,
                   const int32_t *faces, size_t nf,
                   const NeckProbeOpts *opts,
                   NeckProbeResult *out);

/* In-process unit tests (flat sheet, single perpendicular flap, shared-edge
 * crease). Returns 0 if all pass, else the number of failures. */
int NeckProbe_selftest(void);

#endif
