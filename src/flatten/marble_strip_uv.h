#ifndef MARBLE_STRIP_UV_INCLUDED
#define MARBLE_STRIP_UV_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../whole/atlas_strip.h"

/*
 * StrokeStrip-style repair of one structured ribbon crop.
 *
 * The rows are the paper's ordered strokes, U is their common material
 * arclength parameter, and V is the already-established axial coordinate.  A
 * detector-derived trust field changes row/cross-section weights smoothly; it
 * never deletes topology.  Only U is solved, so strict row monotonicity is
 * sufficient to preserve every triangle orientation on the rectangular grid.
 */

typedef struct MarbleStripUvOptions {
    int sample_stride_u;          /* coarse StrokeStrip sample stride (default 4) */
    int sample_stride_v;          /* coarse stroke-row stride (default 4) */
    int l1_iterations;            /* robust initialization cap (default 8) */
    int max_irls_iterations;      /* total relaxed-QP cap, incl. initial/L1 (1000) */
    int final_iterations;         /* optional unknown-cross-section stage (default 0) */
    double final_damping;         /* local/global fixed-point damping (default .1) */
    double movement_tolerance;    /* max coarse U change stop (default 1e-5) */
    double likelihood_sigma;      /* <=0: StrokeStrip span/30 rule */
    double lambda_length;
    double lambda_align;
    double lambda_local;
    double lambda_prior;
    double monotone_fraction;
    double membership_floor;      /* detector-core cross-section prior */
    double match_radius;          /* lateral XYZ candidate reach (default 12) */
    double match_angle_deg;       /* oriented tangent gate (default 35 deg) */
    int match_search_columns;     /* target-row edge search radius (default 24) */
    double topology_fallback;     /* same-column alternative prior (default .02) */
    double length_winsor;         /* max XYZ/latent edge ratio (default 3) */
    double max_displacement;      /* global U trust radius (default 64) */
    int verbose;
} MarbleStripUvOptions;

typedef struct MarbleStripUvStats {
    int coarse_h, coarse_w;
    size_t coarse_samples;
    size_t coarse_cross_sections, coarse_members;
    size_t exact_anchors;
    int robust_solves;
    int final_solves;
    int converged;
    double last_max_u_change;
    double likelihood_sigma;
    size_t downweighted_members;
    double membership_min, membership_mean, membership_max;
    double length_rms_before, length_rms_after;
    double align_rms_before, align_rms_after;
    double local_speed_rms_before, local_speed_rms_after;
    double min_monotone_ratio;
    double safe_commit_alpha;
    size_t moved_vertices;
    double rms_u_displacement, max_u_displacement;
} MarbleStripUvStats;

void MarbleStripUvOptions_default(MarbleStripUvOptions *options);

/*
 * verts       frozen XYZ, row-major [H*W*3]
 * reference_uv current chart, row-major double [H*W*2]
 * trust       clean-chart confidence [0,1], row-major [H*W]
 * pin         exact-zero-displacement samples (nullable), row-major [H*W]
 * output_uv   receives [H*W*2]; V is copied bit-for-bit from reference_uv.
 *
 * The scroll adapter terminates after the paper's relaxed variational solve by
 * default.  Unlike a freehand stroke cluster, this ribbon already has a known
 * axial V foliation, so applying Eq. 7 to newly guessed 3D cross-sections adds
 * a false discrete correspondence problem rather than resolving one.
 */
int MarbleStripUv_solve(Arena_T arena,
                        const float *verts,
                        const double *reference_uv,
                        const double *trust,
                        const uint8_t *pin,
                        int H, int W,
                        const MarbleStripUvOptions *options,
                        double *output_uv,
                        MarbleStripUvStats *stats);

int MarbleStripUv_selftest(void);

#endif
