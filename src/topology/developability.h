#ifndef DEVELOPABILITY_INCLUDED
#define DEVELOPABILITY_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * developability.h -- Stein/Grinspun/Crane "Developability of Triangle Meshes"
 * (SIGGRAPH 2018) variational optimizer.
 *
 * Drives a mesh toward piecewise-developable by minimizing the covariance energy
 *   E_lambda = sum_{interior v} lambda_min(A_v),   A_v = sum_corner theta * N N^T
 * where the sum is over the corners of v's star, theta the corner angle and N the
 * incident face normal. A_v's smallest eigenvalue is zero exactly when the star's
 * normals are coplanar (a "hinge"), i.e. locally flattenable with a ruling line.
 * Curvature concentrates onto sparse seams; everything else flattens.
 *
 * For a scroll this is denoising toward the true (developable) papyrus surface,
 * so afterwards BFF can flatten near-isometrically. Vertices move by analytic-
 * gradient descent with Armijo backtracking, bounded by a per-vertex displacement
 * cap so ink-bearing geometry can't drift, with pinned and boundary vertices held
 * fixed. No remeshing in this version (a documented upgrade).
 * ==========================================================================*/

typedef struct {
    int    max_iters;     /* <=0 -> 100 */
    double tol_grad;      /* stop when max|grad| < this (<=0 -> 1e-6) */
    double tol_e;         /* stop on relative energy stall (<=0 -> 1e-6) */
    double max_disp_vox;  /* per-vertex displacement cap, voxels (<=0 -> 1.0) */
    int    pin_boundary;  /* 1 = hold the boundary loop fixed (default 1) */
    double c1;            /* Armijo constant (<=0 -> 1e-4) */
    double t0;            /* initial step (<=0 -> 1.0) */
    int    max_ls;        /* max line-search backtracks (<=0 -> 30) */
    int    verbose;
} DevelopOpts;

typedef struct {
    double e_initial;
    double e_final;
    int    iters;
    double max_disp;      /* largest realized per-vertex displacement (vox) */
    int    converged;     /* 1 if gradient/stall tolerance reached */
} DevelopStats;

/* Optimize `verts` (modified IN PLACE) toward developability. faces: [nf*3],
 * 0-based. pin_mask: [nv] (1 = never move) or NULL. vert_normals: [nv*3]
 * refreshed if non-NULL. opts/out may be NULL. Returns 0 on success (incl. no-op
 * on trivial input), -1 only on bad arguments. */
int Develop_optimize(Arena_T arena, float *verts, size_t nv,
                     const int32_t *faces, size_t nf,
                     const uint8_t *pin_mask, float *vert_normals,
                     const DevelopOpts *opts, DevelopStats *out);

/* Per-vertex Crane developability energy lambda_v = lambda_min(A_v),
 * A_v = sum_corner theta * N N^T, evaluated on `verts`/`faces` AS GIVEN (no
 * optimization, no vertex movement). lam_out[nv] is filled; an interior vertex's
 * energy is ~0 on a developable (hinge/flat) star and large on a seam where two
 * developable pieces meet. BOUNDARY vertices get 0 (the energy is ill-defined
 * there, so a boundary vertex is never treated as a seam). This is the read-only
 * detector the developability splitter (src/split/dev_cut.c) thresholds.
 *
 * verts: [nv*3] float (any axis order). faces: [nf*3], 0-based. lam_out: [nv]
 * caller-provided (double). Returns 0 on success (incl. all-zero on trivial
 * input), -1 only on bad arguments. arena scratch is saved/restored. */
int Develop_vertex_energy(Arena_T arena, const float *verts, size_t nv,
                          const int32_t *faces, size_t nf, double *lam_out);

int Develop_selftest(void);

#endif
