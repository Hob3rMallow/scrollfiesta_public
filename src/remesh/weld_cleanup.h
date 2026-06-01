#ifndef WELD_CLEANUP_INCLUDED
#define WELD_CLEANUP_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

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

#endif
