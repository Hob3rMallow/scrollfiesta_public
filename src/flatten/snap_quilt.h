#ifndef SNAP_QUILT_INCLUDED
#define SNAP_QUILT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/raw_sample.h"

/* ============================================================================
 * snap_quilt.h -- coherent target selection for the first (repair) snap pass.
 *
 * Only vertices in `dark` are considered.  Candidate positions are sampled at
 * signed depths along the stable axis-radial across-sheet direction.  A bright
 * candidate pays:
 *
 *   quilting mismatch to the nearest clean boundary
 * + distance from the current mesh
 * + a truncated-linear 1-ring depth discontinuity
 *
 * The last term is optimized jointly with GCO alpha-expansion.  Dark vertices
 * with no reachable candidate are fixed at zero displacement, so they also act
 * as barriers against a target field jumping to the opposite side of a sheet.
 * Occupancy is tested against `occ_verts`; a ray stops before other geometry.
 * ==========================================================================*/

typedef struct SnapQuiltOpts {
    float axis_point[3];
    float axis_dir[3];
    double reach;
    double step;
    double occ_thresh;
    double local_r;
    double band;
    double min_gain;
    int target_mode;       /* SNAP_TARGET_* from surface_snap.h */
    int depth_bins;        /* odd number over [-reach,+reach] */
    double w_match;        /* |candidate RAW - propagated clean-boundary RAW| */
    double w_close;        /* |signed depth| */
    double smooth_mu;      /* truncated-linear depth difference on mesh edges */
    double smooth_tau;     /* truncation in voxels */
    int verbose;
} SnapQuiltOpts;

typedef struct SnapQuiltStats {
    size_t n_dark;
    size_t n_fixable;
    size_t n_crack;
    size_t n_bidir;
    size_t n_blocked;
    size_t n_gco_fallback;
    double mean_quilt_cost;
    double mean_distance;
} SnapQuiltStats;

/* `verts`, `normals`, cv/has/dark, and outputs are local arrays of length nv.
 * `occ` and `occ_verts` may cover a larger global mesh.  voff is non-negative;
 * vdir contains the signed winning direction.  vclass values use
 * SNAP_GOOD/SNAP_FIXABLE/SNAP_CRACK numeric parity (0/1/2).
 *
 * Returns 0 on success, -1 on invalid/unsupported input. */
int SnapQuilt_select(Arena_T arena,
                     CubeTable *ct,
                     KDTree_T occ,
                     const float *occ_verts,
                     size_t occ_index_base,
                     const float *verts,
                     const float *normals,
                     size_t nv,
                     CSR_T adj,
                     const double *cv,
                     const uint8_t *has,
                     const uint8_t *dark,
                     const SnapQuiltOpts *opts,
                     uint8_t *vclass,
                     float *voff,
                     float *vdir,
                     float *vgain,
                     float *vblock,
                     SnapQuiltStats *out);

#endif
