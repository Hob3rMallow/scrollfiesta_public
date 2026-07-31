#ifndef QEM_INCLUDED
#define QEM_INCLUDED

#include "arena.h"
#include <stdint.h>
#include <stddef.h>

/*
 * QEM (Quadric Error Metrics) Mesh Simplification
 * Garland-Heckbert algorithm with boundary edge preservation and
 * probabilistic-quadric (Trettner & Kobbelt 2020) regularization.
 *
 * Reduces a triangle mesh to approximately target_nf faces by
 * iteratively collapsing the lowest-cost edge. Boundary-loop verts
 * (incident on any face_count == 1 edge) are auto-pinned to prevent
 * aperture closure.
 *
 * All scratch memory is allocated from the provided arena.
 *
 * Returns 0 on success.
 * If in_nf <= target_nf, copies input unchanged.
 */
int QEM_simplify(Arena_T arena,
                 const float *in_verts, size_t in_nv,
                 const int32_t *in_faces, size_t in_nf,
                 size_t target_nf,
                 float **out_verts, size_t *out_nv,
                 int32_t **out_faces, size_t *out_nf);

/*
 * Same as QEM_simplify but accepts an optional `pin_mask` (length in_nv,
 * 1 = pinned) and an optional `ref_normal` (unit vector, used as the
 * reference normal for the final winding-repair pass instead of the
 * area-weighted average computed from the post-collapse mesh).
 *
 * The ref_normal is critical for cross-cube determinism: adjacent cubes
 * processing the same geometric surface can end up with opposite
 * area-weighted normals (their meshes are different halves of the same
 * sheet) and so flip face windings in opposite directions, producing
 * the "same_dir" winding-mismatch pairs that grid_weld counts as bad.
 * Pass cm->pca_normal (which is itself (1,1,1)-hemisphere oriented per
 * src/common/pca.c) to make the winding choice cube-shared.
 *
 * Pinned-vertex semantics:
 *   - both endpoints pinned -> edge is un-collapsible
 *   - one endpoint pinned   -> collapse keeps the pinned position exactly
 *   - pins propagate during collapse (a survivor at a pinned position
 *     becomes pinned)
 *
 * If `out_pin_mask` is non-NULL, the function emits a compacted pin_mask
 * (length *out_nv, 1 = pinned) into a fresh arena allocation and writes
 * the pointer there. This is the post-collapse pin set, indexed by the
 * NEW vertex indices in *out_verts. It also reflects propagation: a vert
 * that absorbed a pinned vertex appears pinned in the output. Callers
 * who do not need this can pass NULL.
 *
 * Pass NULL for pin_mask and/or ref_normal to fall back to the default.
 */
int QEM_simplify_pinned(Arena_T arena,
                        const float *in_verts, size_t in_nv,
                        const int32_t *in_faces, size_t in_nf,
                        const uint8_t *pin_mask,
                        const float *ref_normal,
                        size_t target_nf,
                        float **out_verts, size_t *out_nv,
                        int32_t **out_faces, size_t *out_nf,
                        uint8_t **out_pin_mask);

/*
 * Full-control variant. QEM_simplify / QEM_simplify_pinned are wrappers
 * around this with default options (auto-pin boundary, winding repair on).
 *
 *   free_boundary       -- skip the Phase-4b auto-pin of boundary-loop
 *                          vertices. The boundary QUADRIC penalty still
 *                          applies (Garland-Heckbert constraint planes), so
 *                          open edges resist drifting, but hole rims no
 *                          longer survive at full input resolution. Use for
 *                          hole-riddled welded sheets where uniform
 *                          decimation matters more than rim preservation.
 *   skip_winding_repair -- keep each surviving face's INPUT winding: skip
 *                          both global reference-normal flip repairs (the
 *                          compaction vote and the final iterative pass).
 *                          Those repairs assume a mostly-planar patch; on a
 *                          surface wrapped around an axis ~half of all face
 *                          normals oppose ANY global reference, so the
 *                          "repair" flips whole regions. A consistently
 *                          oriented input stays consistently oriented.
 *   prob_quadrics       -- use the FULL per-face anisotropic probabilistic
 *                          plane quadric (Trettner & Kobbelt 2020) instead of
 *                          the deterministic n·nᵀ quadric + per-vertex σ²I hack.
 *                          Each face plane is a noisy observation: normal
 *                          n ~ N(n̄, σ_n²I), point p ~ N(p̄, σ_p²I), and the
 *                          face quadric is E[(nᵀ(x−p))²]. The σ_n²I term makes
 *                          the 3×3 A block full-rank on flat/degenerate patches,
 *                          so the optimal-position solve is well conditioned and
 *                          places verts robustly -> fewer slivers / more uniform
 *                          triangles. Reduces to the deterministic quadric at
 *                          σ=0. Default OFF (deterministic path unchanged).
 *   prob_sigma_n/_p     -- noise magnitudes (vox). 0 => auto from the median
 *                          edge length (QEM_PROB_SIGMA_FACTOR / _POS_FACTOR).
 *                          Only consulted when prob_quadrics != 0.
 */
typedef struct {
    const uint8_t *pin_mask;      /* optional [in_nv] caller pins */
    const float   *ref_normal;    /* optional winding reference (unit) */
    int            free_boundary;
    int            skip_winding_repair;
    int            prob_quadrics;  /* full per-face probabilistic quadrics (see above) */
    float          prob_sigma_n;   /* normal-noise σ   (0 => auto from median edge) */
    float          prob_sigma_p;   /* position-noise σ (0 => auto from median edge) */
} QemOpts;

int QEM_simplify_opts(Arena_T arena,
                      const float *in_verts, size_t in_nv,
                      const int32_t *in_faces, size_t in_nf,
                      size_t target_nf,
                      const QemOpts *opts,
                      float **out_verts, size_t *out_nv,
                      int32_t **out_faces, size_t *out_nf,
                      uint8_t **out_pin_mask);

#endif
