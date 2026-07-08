#ifndef DEPTH_PEEL_INCLUDED
#define DEPTH_PEEL_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"

/* ============================================================================
 * depth_peel.h -- Stacked-wrap separator. The FIRST separation pass in the
 * per-cube pipeline (Step 1b-peel), ahead of dev-cut, bridge-cut and overlap-sep.
 *
 * Ball-Pivoting (rho = 1.2 vox, no inter-wrap clearance) fuses papyrus wraps that
 * sit only ~2-3 vox apart into ONE connected component -- a stack of near-parallel
 * sheets welded by a few "wall" bridge triangles. Connectivity-resplit can't undo
 * it (it's one component), developability can't see it (each wrap is locally flat,
 * lambda ~ 0), and the overlap-sep lifted multicut CAN separate it but costs
 * 50-70 s on a big clump. The cheap signal is DEPTH along the surface normal: the
 * fused stack is the only place a mesh edge JUMPS across the inter-wrap gap.
 *
 * DepthPeel_process is pure (mesh in -> connected pieces out): project every vertex
 * onto the component PCA normal to get a depth w; mark a vertex as a seam vertex if
 * any incident triangle edge has |w_a - w_b| > min_gap (an inter-wrap bridge edge --
 * a long, near-normal jump that intra-sheet edges, ~0.6 vox, never make); dilate the
 * seam band by gap_depth (KD-tree, like bridge-cut's CUT_GAP_DEPTH); return the
 * surviving connected components. It mirrors dev_cut.c exactly, swapping the
 * per-vertex Crane-energy seam test for the per-edge depth-jump test.
 *
 * Why per-edge depth-jump and not UV-cell binning: it is connectivity-based, so it
 * is robust to a tilted/curved stack (depth is compared only between adjacent
 * vertices, never globally) and it SELF-GATES -- a single sheet, however folded,
 * has only short edges, hence no large depth-jump, hence one piece (the "NO
 * intra-sheet split" guarantee, not a tuned threshold). A jump only appears where
 * BPA bridged the empty inter-wrap gap. The caller (pipeline_cube) re-LOPs each new
 * piece from its OWN points; pieces that DepthPeel can't separate (a strongly curved
 * stack with no clean per-edge jump) fall through to overlap-sep unchanged.
 * ==========================================================================*/

/* Peel `mesh` into its stacked wrap layers.
 *
 *   arena          - all allocations (scratch saved/restored internally).
 *   mesh           - input component (read-only); uses mesh->pca_normal as the
 *                    depth (thickness) axis.
 *   min_gap        - depth-jump threshold (voxels) on |w_a - w_b| across a triangle
 *                    edge (DEPTH_PEEL_MIN_GAP_VOX). Above the intra-sheet edge length
 *                    (~0.6 vox), below the inter-wrap spacing (~2-3 vox).
 *   gap_depth      - exclusion radius (voxels) around seam verts (DEPTH_PEEL_GAP_DEPTH);
 *                    <= 0 removes only the seam verts themselves.
 *   min_comp_verts - drop survivor layers smaller than this (DEPTH_PEEL_MIN_COMP_VERTS);
 *                    a fold flap is small and dropped, a real wrap is large and kept.
 *   out_meshes     - arena array of sub-meshes (ORIGINAL-position space), >= 1.
 *   out_count      - number of sub-meshes. == 1 means "no separation" (self-gated).
 *
 * Returns 0 on success (incl. the no-split pass-through), -1 only on bad input. */
int DepthPeel_process(Arena_T arena, const ComponentMesh *mesh,
                      double min_gap, double gap_depth, size_t min_comp_verts,
                      ComponentMesh **out_meshes, size_t *out_count);

/* In-process unit tests. Returns 0 if all pass, else the failure count. */
int DepthPeel_selftest(void);

#endif
