#ifndef WELD_CLEANUP_INCLUDED
#define WELD_CLEANUP_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"
#include "seam_planes.h"

/*
 * Post-weld geometric cleanup of bridge slivers + T-junctions.
 *
 * SeamWeld_bridge (BPA) fuses adjacent cube meshes by rolling a ball over the
 * seam; near the seam it leaves thin "sliver" triangles, the occasional
 * zero-area face, and T-junctions (a bridge vertex landing on the interior of a
 * coarse component triangle's edge). These are NOT topological defects -- the
 * weld stays edge-manifold -- but they are geometric blemishes the per-cube QEM
 * never sees, because QEM runs BEFORE the weld and grid_weld is the terminal
 * step. This is the missing post-weld edge-collapse.
 *
 * Strategy (Surazhsky-Gotsman flip first, guarded collapse for the residue):
 *   Pass 1  min-angle edge flips (normal-guarded, vertices never moved, boundary
 *           edges never flipped) -- clears CAP slivers without removing anything.
 *   Pass 2  guarded short-edge collapse of the remaining sliver / zero-area
 *           faces -- clears NEEDLES, degenerate faces, and degenerate-cap
 *           T-junction vertices. Each collapse is gated by:
 *             - Dey-Edelsbrunner-Guha link condition (never creates a
 *               non-manifold edge),
 *             - normal-flip guard (no incident face reverses / degenerates),
 *             - max-edge-length cap (an edge longer than this is never collapsed,
 *               so even a hypothetical inter-wrap face could not merge two wraps;
 *               real mesh edges already join one sheet, this is belt-and-braces),
 *             - boundary pin (a vertex on a boundary loop is never the mover, so
 *               the outer perimeter and genuine holes are preserved exactly).
 *   Pass 3  one more flip round to tidy caps the collapse exposed.
 *
 * Vertices are never created or moved. A collapsed vertex is left ORPHANED (its
 * incident faces are remapped to the survivor and the <=2 degenerate faces are
 * dropped) so any caller-side per-vertex array (colors, cube indices, weld
 * flags) stays valid by index. cm->faces is arena-reallocated + compacted and
 * cm->nf is updated; cm->verts / cm->nv are unchanged.
 *
 * Run AFTER FoldCleanup + PinholeFill + cull, BEFORE the final winding/manifold
 * audit (the audit then double-checks manifoldness was preserved).
 *
 * Returns 0 on success. stats may be NULL.
 */

typedef struct {
    double sliver_min_alt;    /* triangle min-altitude < this (vox)  => sliver target  */
    double degen_area;        /* triangle area < this (vox^2)        => degenerate target */
    double max_collapse_len;  /* never collapse an edge longer than this (vox)         */
    int    flip_max_rounds;   /* cap on pass-1 / pass-3 flip iterations                */
    int    collapse_max_rounds;/* cap on pass-2 collapse iterations                    */
} WeldCleanupParams;

typedef struct {
    size_t n_flips;       /* total edge flips (pass 1 + pass 3)        */
    size_t n_collapses;   /* total accepted edge collapses (pass 2)    */
    size_t faces_in;      /* face count on entry                       */
    size_t faces_out;     /* face count on exit                        */
    size_t targets_in;    /* sliver + degenerate faces before          */
    size_t targets_out;   /* sliver + degenerate faces after           */
} WeldCleanupStats;

/* Fill p with defaults sourced from pipeline_constants.h. */
void WeldCleanup_default_params(WeldCleanupParams *p);

int WeldCleanup_process(Arena_T arena, ComponentMesh *cm,
                        const WeldCleanupParams *params,
                        WeldCleanupStats *stats);

/*
 * Seam-band recoarsening -- collapse the temporarily-fine seam band back toward
 * the coarse budget AFTER the bridge + hole fills have closed the seam.
 *
 * Candidate = a face whose SHORTEST edge is < collapse_below AND both endpoints
 * of that edge lie within `band` of a detected seam plane. Everything else is
 * the proven collapse machinery unchanged: only the shortest edge collapses,
 * length <= max_collapse_len (5 < 7-vox clearance: a collapse can never fuse
 * wraps -- and structurally it only contracts EXISTING edges, so it can never
 * connect components at all), Dey-Edelsbrunner-Guha link condition, normal-flip
 * guard, both-boundary edges skipped, boundary vertex always the survivor,
 * 1-ring locking, rounds to fixpoint with interleaved boundary-frozen flips.
 *
 * Boundary loops (grid edges, unbridged holes, a hierarchical level's outer
 * faces) exit with their vertex set, positions, and cyclic order BIT-IDENTICAL
 * -- simultaneously the no-reopened-seam guarantee and the requirement that
 * level L+1 of hierarchical_weld sees exactly the boundaries level L was given.
 *
 * Vertices are never created or moved (collapsed verts orphaned, caller
 * per-vertex arrays stay index-valid). cm->faces/nf updated as in
 * WeldCleanup_process. Returns 0 on success. st may be NULL.
 */

typedef struct {
    double band;              /* vox from a detected seam plane => in-band       */
    double collapse_below;    /* collapse band edges shorter than this (vox)     */
    double max_collapse_len;  /* hard cap; keep WELD_CLEANUP_MAX_COLLAPSE_LEN    */
    int    max_rounds;        /* collapse-round cap                              */
    int    flip_max_rounds;   /* interleaved / final flip iterations             */
} WeldRecoarsenParams;

typedef struct {
    size_t n_collapses;
    size_t n_flips;
    size_t faces_in;
    size_t faces_out;
} WeldRecoarsenStats;

void WeldCleanup_default_recoarsen_params(WeldRecoarsenParams *p);

int WeldCleanup_recoarsen_seam(Arena_T arena, ComponentMesh *cm,
                               const SeamPlane *planes, size_t np,
                               const WeldRecoarsenParams *params,
                               WeldRecoarsenStats *st);

/* Boundary-frozen Surazhsky-Gotsman min-angle flip rounds on a raw face
 * array (mutated in place; verts never move; boundary edges never flip;
 * normal-guarded). Thin public wrapper over the internal flip machinery so
 * the seam-band refiner can relieve its split fans without duplicating the
 * flip pass. Returns total flips. */
size_t WeldCleanup_flip_rounds(Arena_T arena, const float *verts, size_t nv,
                               int32_t *faces, size_t nf, int max_rounds);

#endif
