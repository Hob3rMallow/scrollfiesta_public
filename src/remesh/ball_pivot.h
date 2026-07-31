#ifndef BALL_PIVOT_INCLUDED
#define BALL_PIVOT_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * Ball-Pivoting Algorithm (Bernardini, Mittleman, Rushmeier, Silva,
 * Taubin 1999). Triangulates an oriented point cloud by rolling a
 * ball of radius rho over the surface and emitting a triangle for
 * each (a, b, c) the ball touches simultaneously.
 *
 * This implementation re-seeds when the front empties so the whole
 * cloud is covered (the textbook BPA stops at the first front empty).
 * Manifold guarantee: each output edge is in ≤ 2 faces, enforced by
 * a pre-check that aborts any pivot whose new edges would be the
 * third reference to an interior edge. The exemplary port at
 * D:\work\bpa adds an "inner-edge already exists" rejection inside
 * the pivot which we also do.
 *
 * Determinism: vertex iteration is in index order; seed candidate
 * pairs are tried in nearest-first order; the pivot keeps the K
 * smallest-theta candidates (ties broken by vertex index) and the front
 * takes the first that passes the manifold / fold / anti-parallel guards
 * -- a deterministic choice. Halo-deterministic when verts are
 * bit-identical across cubes (LOP+pin-snap guarantees this for the halo
 * region).
 *
 * Inputs:
 *   verts[nv*3]    float (z, y, x)
 *   normals[nv*3]  float unit normals (z, y, x); LOP-PCA normals
 *                  satisfy the local-consistency requirement.
 *   nv             vertex count
 *   rho            pivot ball radius (voxels)
 *
 * Outputs (arena-allocated):
 *   *out_faces  [nf*3] int32 — 0-based vertex indices
 *   *out_nf     face count
 *
 * Returns 0 on success, non-zero on hard failure.
 */
int BallPivot_reconstruct(Arena_T arena,
                          const float *verts,
                          const float *normals,
                          size_t nv,
                          float rho,
                          int32_t **out_faces,
                          size_t *out_nf);

/*
 * Optional growth-scoped winding anchor for per-cube reconstruction.
 *
 * A local pairwise winding test is insufficient when BPA walks gradually up a
 * through-thickness wall: every individual edge can be innocent while one
 * growth front accumulates several complete turns.  The gated entry point
 * integrates the branch-cut-free local phase step
 *
 *     dw = dr / pitch - dtheta / (2*pi)
 *
 * from each seed and rejects a pivot once the accumulated phase leaves
 * [-tol,+tol].  Re-seeding then reconstructs the next wrap as a distinct mesh
 * component.  `verts` remain cube-local; origin_y/x convert them to the world
 * frame containing umb_y/x.
 *
 * Pass NULL (or pitch/tol <= 0) for the historical ungated behaviour.
 */
typedef struct {
    double umb_y, umb_x;       /* world-space scroll umbilicus */
    double origin_y, origin_x; /* cube-local -> world translation */
    double pitch;              /* radial spacing per turn (vox) */
    double tol;                /* maximum accumulated seed phase (turns) */
} BpaReconGate;

int BallPivot_reconstruct_gated(Arena_T arena,
                                const float *verts,
                                const float *normals,
                                size_t nv,
                                float rho,
                                const BpaReconGate *gate,
                                int32_t **out_faces,
                                size_t *out_nf);

/* Populate a per-cube gate from BPA_GROW_UMBILICUS_Y/X,
 * BPA_GROW_WRAP_PITCH and BPA_GROW_WIND_TOL. `origin_zyx` is the cube's world
 * origin; NULL means {0,0,0}. Returns 1 when all required values are valid and
 * the gate is armed, otherwise zeroes `out` and returns 0. */
int BpaReconGate_from_env(BpaReconGate *out, const float origin_zyx[3]);

/*
 * Seam-bridge BPA. Same rolling-ball machinery as BallPivot_reconstruct,
 * but instead of auto-seeding it is primed with an explicit initial front
 * of boundary edges and never re-seeds. Used to weld two independently
 * triangulated cube meshes across the halo seam: the open boundary loops
 * of both cubes are fed as the initial front, BPA rolls across the gap
 * over the supplied band points, and the two fronts meet in the middle —
 * producing a manifold bridge that shares the existing boundary vertices
 * (so its edges pair with the cubes' interior faces). This is the
 * "treat the halo as a new BPA weld" path: it sidesteps the fact that two
 * cubes flooding BPA from opposite sides never agree on a shared seam
 * triangulation.
 *
 * Each init edge is a boundary half-edge (va, vb) of an existing mesh
 * plus the third vertex v_opp of its single incident face. The hinge ball
 * center is reconstructed from that face (outward side, per the same
 * normal-sign convention as seeding) so the pivot advances INTO the gap,
 * away from the existing surface.
 *
 * verts/normals span ALL participating vertices: both cubes' boundary
 * verts and the band stepping-stone points, in one flat (z,y,x) array.
 * init_edges index into that array. Output faces also index into it.
 *
 * The front is grown at escalating radii rho, 1.5*rho, 2*rho (capped at
 * rho_max), so gaps the smallest ball falls through are spanned by a larger
 * one — without ever exceeding rho_max, which bounds 2*rho below the
 * inter-wrap clearance (NO inter-wrap mergers). Pass rho_max == rho for a
 * single-radius bridge.
 *
 * Emitted faces are wound consistently with the existing faces that own the
 * init edges: the directed-front glue (Bernardini §IV-D) maintains a coherent
 * orientation and rejects same-direction (non-orientable) coincidences, so the
 * bridge pairs cleanly. A single global winding pass on the combined mesh then
 * reconciles the two cubes' independent orientations.
 *
 * Returns 0 on success.
 */
typedef struct {
    int32_t va, vb, v_opp;
} BpaInitEdge;

/* Winding-gate parameters for the seam bridge. The gate rejects a bridge that
 * connects two DIFFERENT wraps, by winding phase w = r/pitch - theta/2pi about
 * the scroll umbilicus (see the gate in bpa_try_candidate). ARMED only when
 * pitch > 0 AND an umbilicus is supplied; otherwise the gate is OFF (the caller
 * uses this for a restricted-cloud permissive weld, where cloud restriction is
 * the safety, not the gate). tol = phase tolerance (turns); hard = unconditional
 * cap (turns, <=0 = none). Pass NULL to BallPivot_bridge for gate-off. */
typedef struct {
    double umb_y, umb_x;
    double pitch;
    double tol;
    double hard;
} BpaBridgeGate;

int BallPivot_bridge(Arena_T arena,
                     const float *verts,
                     const float *normals,
                     size_t nv,
                     float rho,
                     float rho_max,
                     const BpaInitEdge *init_edges,
                     size_t n_init,
                     const BpaBridgeGate *gate,
                     int32_t **out_faces,
                     size_t *out_nf);

#endif
