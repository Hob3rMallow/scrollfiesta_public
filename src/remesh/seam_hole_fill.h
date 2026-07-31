/* seam_hole_fill.h -- Merger-safe closure of small seam punctures left open by
 * the BPA seam bridge.
 *
 * The audit of the 4x5x5 grid (output/weld_audit_4x5x5/AUDIT.md) showed the
 * production weld leaves ~144 single-triangle + ~500 small holes straddling the
 * cube-boundary seams: the BPA bridge zips most of a grazing seam but leaves a
 * dotted line of tiny punctures, and the downstream closers do not finish them
 * (pinhole's 4.5-vox diameter gate rejects the wider ones; interior-only
 * hole-fill skips the boundary-connected ones).
 *
 * This pass closes exactly those punctures, with three gates that make an
 * inter-wrap merger impossible:
 *   1. SEAM + STRADDLE -- the loop hugs a cube-boundary plane and has verts on
 *      BOTH sides of it (a genuine weld miss between two cubes of one sheet),
 *      not a per-cube interior notch or an outer bay.
 *   2. SMALL EXTENT    -- bbox diagonal below the inter-wrap clearance, so the
 *      hole cannot span two different wraps by size alone.
 *   3. WINDING PHASE   -- (when the umbilicus + pitch are supplied, as the seam
 *      weld gate is) every loop vertex lies on the SAME wrap: the loop's phase
 *      span is under the gate tolerance. This is what protects the dense core,
 *      where wraps sit closer than the nominal clearance -- a hole spanning two
 *      fused wraps has a ~1-turn phase span and is left open.
 *
 * Only SIMPLE boundary loops (every vertex boundary-degree 2) of length
 * [3, max_loop] are considered; the fill is a fan triangulation wound to match
 * the surface, and is applied only if every fan triangle is manifold-safe
 * (no duplicate/>2-face edge), non-degenerate, and normal-consistent -- else the
 * loop is left open. Vertices are never moved or added; only faces are appended.
 *
 * Run late in grid_weld (after the bridge + interior hole-fill, before the final
 * manifold guard). Default ON; env SEAM_NO_HOLEFILL2 / --no-seamfill disables.
 */
#ifndef SEAM_HOLE_FILL_H
#define SEAM_HOLE_FILL_H

#include "../common/arena.h"
#include "../common/mesh_types.h"

typedef struct SeamHoleFillParams {
    double cube;            /* seam-plane spacing (voxels), e.g. 128 */
    double band;            /* a loop is near-seam if within this of a plane */
    double max_extent;      /* bbox-diagonal cap (vox); phase gate is the real
                             * wrap guard when armed (unarmed clamps to 4) */
    int    max_loop;        /* largest boundary loop (edges) to fill, e.g. 12 */
    double near_vox;        /* PHASE-GATED ONLY: also accept a non-straddling
                             * loop whose bbox lies entirely within this of a
                             * seam plane. The weld's trim inset + bridge pull
                             * some seam punctures fully to one side of the
                             * plane (observed 1-2 vox); the strict straddle
                             * test left those open forever. */
    /* Winding-phase gate (mirrors the seam-weld gate). pitch <= 0 disables it,
     * in which case a tighter extent cap is used instead. */
    double umb_y, umb_x, pitch;
    double wind_tol_turns;  /* max loop phase span in turns (e.g. 0.25) */
} SeamHoleFillParams;

typedef struct SeamHoleFillStats {
    size_t loops;              /* simple boundary loops examined */
    size_t seam_candidates;    /* near-seam straddling small loops */
    size_t filled;             /* loops closed */
    size_t tris_added;         /* fill triangles appended */
    size_t skip_phase;         /* rejected: verts span > tol turns (different wraps) */
    size_t skip_extent;        /* rejected: bbox diagonal too large */
    size_t skip_geom;          /* rejected: degenerate / non-convex fan */
    size_t skip_manifold;      /* rejected: fill would make a >2-face edge / dup */
    size_t skip_bubble;        /* rejected: sole boundary loop of its component */
} SeamHoleFillStats;

/* Fill p->cube/band/etc with sensible defaults (cube 128, band 6, max_extent 8,
 * max_loop 12, near_vox 3, wind_tol 0.25, phase gate disabled until
 * umbilicus/pitch set — unarmed, the extent cap clamps to 4 and the straddle
 * requirement stays strict). */
void SeamHoleFill_default_params(SeamHoleFillParams *p);

/* Close merger-safe seam punctures in `cm`. cm->faces is arena-reallocated and
 * grown; cm->verts/nv are unchanged. Returns 0 on success. */
int SeamHoleFill_process(Arena_T arena, ComponentMesh *cm,
                         const SeamHoleFillParams *p, SeamHoleFillStats *st);

#endif /* SEAM_HOLE_FILL_H */
