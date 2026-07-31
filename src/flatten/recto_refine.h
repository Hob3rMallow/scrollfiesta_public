#ifndef RECTO_REFINE_INCLUDED
#define RECTO_REFINE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/raw_sample.h"

/* ============================================================================
 * recto_refine.h -- second-pass global recto placement.
 *
 * The recto is the innermost-facing boundary of a wound papyrus sheet.  Along
 * an OUTWARD-oriented across-sheet line, that boundary is the dark-to-bright
 * RAW transition (positive derivative), not the intensity maximum.  Each outer
 * iteration finds that oriented sub-voxel edge for every supported vertex,
 * then minimizes one global scalar-displacement objective:
 *
 *   sum_i w_i (d_i - target_i)^2 + smooth_mu sum_(i,j) (d_i - d_j)^2
 *
 * Unsupported vertices are zero-displacement anchors.  The scalar field is
 * solved by matrix-free Jacobi on the mesh CSR, total-displacement
 * slope-limited on every edge, damped, and guarded against other original
 * geometry. Repeating target-find
 * + global solve updates the local normals without ever rewarding the bright
 * papyrus core or the outward (verso) edge.
 * ==========================================================================*/

typedef struct RectoRefineOpts {
    float axis_point[3];
    float axis_dir[3];
    double window_lo;
    double window_hi;
    double band_frac;       /* transition isovalue in the contrast window */
    double range;           /* local signed search radius, vox */
    double sample_step;     /* profile spacing, vox */
    double edge_half;       /* derivative half-width, vox */
    double min_edge_frac;   /* min positive edge / (window_hi-window_lo) */
    double distance_bias;   /* RAW derivative units per voxel */
    double normal_blend;    /* 0 radial, 1 outward-oriented mesh normal (def .85) */
    double data_weight;     /* supported oriented-edge target weight */
    double anchor_weight;   /* unsupported zero-displacement weight (def 24) */
    double smooth_mu;       /* global 1-ring displacement smoothness */
    double max_slope;       /* max added |di-dj| / original edge length */
    double damping;
    double max_step;        /* per outer iteration, vox */
    double max_total;       /* from second-pass input, vox */
    double occ_thresh;
    double local_r;
    int outer_iters;
    int inner_iters;
    double tol;
    int threads;
    int verbose;
} RectoRefineOpts;

void RectoRefineOpts_default(RectoRefineOpts *o);

typedef struct RectoRefineStats {
    size_t n_supported;
    size_t n_moved;
    size_t n_reverted;
    size_t n_slope_limited;
    int iterations;
    double mean_disp;
    double max_disp;
    double mean_abs_step;
    double mean_recto_gradient;
    double sec_total;
} RectoRefineStats;

/* Mutates verts in place.  The optimization is global over every face-connected
 * vertex in the supplied mesh; disconnected components naturally factorize.
 * ct must be safe for concurrent reads (the function pre-warms it serially).
 * Empty input is a clean no-op. */
int RectoRefine_run(Arena_T arena,
                    CubeTable *ct,
                    float *verts,
                    size_t nv,
                    const int32_t *faces,
                    size_t nf,
                    const RectoRefineOpts *opts,
                    RectoRefineStats *out);

/* Synthetic dark/bright slab test: converges to the inner edge, rejects the
 * outward edge, is near-idempotent, and keeps a smooth displacement field. */
int RectoRefine_selftest(void);

#endif
