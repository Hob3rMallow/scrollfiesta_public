#ifndef REEB_INCLUDED
#define REEB_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * reeb.h -- per-component topological neck finder for bridge diagnosis.
 *
 * The UV-difference oracle (oracle.c) counts sheets by gaps along ONE global
 * PCA normal and is blind to a side-facing bridge that smoothly interpolates
 * that gap. This module is global-normal-free: it sweeps a geodesic Morse
 * function across the component and looks at the CUT PROFILE
 *
 *     C(l) = number of mesh edges whose Morse span straddles level l
 *
 * which is the size of the cross-section at "depth" l. A thin bridge between
 * two wraps shows up as an INTERIOR LOCAL MINIMUM of C(l) with substantial
 * "shoulders" (wide cross-sections) on BOTH sides -- the front fills one wrap,
 * pinches through the bridge, then re-widens into the next. A clean disk's
 * profile decreases monotonically (no interior valley); a single scroll wrap
 * (annulus) has no thin pinch -- so neither is flagged. This is the false-
 * positive guard that distinguishes a bridge from honest geometry.
 *
 * The Morse function is integer geodesic hop-distance from a graph-diameter
 * endpoint (a two-sweep BFS), so the sweep runs ALONG the structure and a
 * bridge becomes an interior bottleneck. Hop distance is insensitive to the
 * sliver geometry BPA produces.
 *
 * Intended to be run per connected component. Cross-check its necks against the
 * geometric side-facing signal (neck_probe.h) and per-component genus
 * (mesh_topo.h): a real bridge is a thin neck AND side-facing.
 *
 * Allocation: ALL outputs and scratch come from `arena` with no internal
 * save/restore (the caller owns the mark, as in bridge_cut.c::detect_bridge).
 * ==========================================================================*/

typedef struct {
    int   neck_max_width; /* a valley with <= this many crossing edges is a
                           * candidate neck (default 8). */
    float rebound_factor; /* require both shoulders >= factor * valley width
                           * (default 2.5) -- a true neck pinches then re-widens. */
    int   min_shoulder;   /* ...and both shoulders >= this absolute (default 12),
                           * which kills the disk's monotone close-up + noise. */
} ReebOpts;

typedef struct {
    uint32_t level;     /* deepest Morse level of the valley */
    uint32_t lo, hi;    /* valley group level range (lo <= level <= hi) */
    int32_t  width;     /* crossing-edge count at the valley (cut size) */
    int32_t  shoulder;  /* min of the two surrounding peak widths */
    size_t   n_faces;   /* #faces straddling the valley band (the neck band) */
} ReebNeck;

typedef struct {
    uint32_t  morse_max;   /* highest Morse level */
    int32_t  *cut_profile; /* [morse_max+1] crossing edges per level (arena) */
    int32_t   bottleneck;  /* smallest reported neck width, or -1 if none */
    size_t    n_necks;     /* #necks reported */
    ReebNeck *necks;       /* [n_necks] (arena), or NULL */
    uint8_t  *face_flag;   /* [nf]: 1 = on a neck band (arena) */
    uint32_t *morse;       /* [nv]: Morse field (arena) */
    int       had_boundary;/* reserved; 0 = farthest-point seed (current) */
} ReebResult;

/* Defaults: neck_max_width 8, rebound_factor 2.5, min_shoulder 12. */
void ReebOpts_default(ReebOpts *opts);

/* Analyze one (preferably single-component) triangle mesh. verts is unused for
 * the topology (kept for signature symmetry / future geometric refinement).
 * Returns 0 on success, -1 on empty input. Graceful on zero/edge-free meshes. */
int Reeb_analyze(Arena_T arena,
                 const float *verts, size_t nv,
                 const int32_t *faces, size_t nf,
                 const ReebOpts *opts,
                 ReebResult *out);

/* In-process unit tests (disk, annulus, two grids bridged by one triangle).
 * Returns 0 if all pass, else the number of failures. */
int Reeb_selftest(void);

#endif
