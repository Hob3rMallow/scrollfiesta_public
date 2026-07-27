#ifndef DEV_CUT_INCLUDED
#define DEV_CUT_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_types.h"

/* ============================================================================
 * dev_cut.h -- Developability cut (Stein/Grinspun/Crane, "Developability of
 * Triangle Meshes", SIGGRAPH 2018, sec 4.4). The FIRST splitter in the per-cube
 * pipeline, ahead of bridge-cut and overlap-sep.
 *
 * Their cut strategy: "identify vertices with energy above a user-specified
 * tolerance eps > 0, then compute a cut passing through all such vertices." The
 * per-vertex energy lambda_v = lambda_min(A_v), A_v = sum_corner theta * N N^T
 * (Develop_vertex_energy), is ~0 on a developable (hinge/flat/cone) star and
 * large where the normals do not fit a plane -- a SEAM where the surface stops
 * being one developable piece. For a scroll that is exactly where two wraps were
 * fused, or a real fold ridge -- the RIGHT place to cut.
 *
 * DevCut_process is pure (mesh in -> connected pieces out): threshold lambda,
 * remove the seam band (+ a geometric dilation, like bridge-cut's CUT_GAP_DEPTH),
 * and return the surviving connected components. A clean developable sheet has no
 * spanning seam, so nothing disconnects and it returns ONE piece -- the
 * self-gate that enforces "NO intra-sheet split". The caller (pipeline_cube)
 * re-LOPs each new piece from its own points and keeps the split only if the
 * pieces are measurably more developable than the parent (DevCut_nondev_fraction
 * gate); otherwise it falls through to bridge-cut / overlap-sep unchanged.
 * ==========================================================================*/

/* Cut `mesh` through its non-developable seam(s).
 *
 *   arena          - all allocations (scratch saved/restored internally).
 *   mesh           - input component (read-only).
 *   eps            - seam threshold on lambda_v (DEV_CUT_EPS).
 *   gap_depth      - exclusion radius (voxels) around seam verts (DEV_CUT_GAP_DEPTH);
 *                    <= 0 removes only the seam verts themselves.
 *   min_comp_verts - drop survivor components smaller than this (DEV_CUT_MIN_COMP_VERTS).
 *   out_meshes     - arena array of sub-meshes (ORIGINAL-position space), >= 1.
 *   out_count      - number of sub-meshes. == 1 means "no separation" (the input
 *                    passed through, self-gated).
 *
 * Returns 0 on success (incl. the no-split pass-through), -1 only on bad input. */
int DevCut_process(Arena_T arena, const ComponentMesh *mesh,
                   double eps, double gap_depth, size_t min_comp_verts,
                   ComponentMesh **out_meshes, size_t *out_count);

/* Fraction of INTERIOR vertices whose Crane energy lambda_v exceeds `eps`
 * (0 = fully developable, 1 = all seam). Boundary vertices are excluded (their
 * energy is ill-defined). This is the lambda-energy gate metric: a dev-cut split
 * is accepted only when every re-LOP'd piece's fraction drops vs the parent's.
 * out_mean_lambda (optional) receives the mean interior lambda for logging.
 * Returns the fraction (0 on trivial/empty input). */
double DevCut_nondev_fraction(Arena_T arena, const ComponentMesh *mesh,
                              double eps, double *out_mean_lambda);

/* In-process unit tests. Returns 0 if all pass, else the failure count. */
int DevCut_selftest(void);

#endif
