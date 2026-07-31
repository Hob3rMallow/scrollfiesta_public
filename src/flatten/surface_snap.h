#ifndef SURFACE_SNAP_INCLUDED
#define SURFACE_SNAP_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

struct CubeTable;

/* ============================================================================
 * surface_snap.h -- pass-1 DETECTION and target selection for obvious dark
 * whiffs.  Candidate brightness is only an eligibility test: the selected
 * signed depth minimizes clean-boundary quilting mismatch plus snap distance,
 * with a joint truncated-linear 1-ring depth prior.  This avoids the old
 * intensity-maximum bias toward the papyrus core or the opposite (verso) edge.
 *
 * Probes stop on OTHER existing geometry.  Only the source vertex and its
 * actual mesh 1-ring are excluded from occupancy, so a spatially close verso is
 * still a blocker.  Dark vertices with no reachable candidate stay in place as
 * cracks/anchors.  This stage only classifies and records targets; SnapSolve
 * performs the pass-1 move and RectoRefine performs the subsequent global,
 * oriented recto-edge placement.
 *
 * "No-geometry" dark (unfilled holes) is OUT OF SCOPE (a texture-coverage
 * concern, not a vertex) and deferred.
 * ==========================================================================*/

typedef struct SnapOpts {
    float  axis_point[3];   /* (z,y,x) scroll axis point (umbilicus) */
    float  axis_dir[3];     /* (z,y,x) axis direction (normalized internally) */
    const char *raw_dir;    /* cubes_RAW dir */
    long   chunk;           /* cube edge, vox (128) */
    double reach;           /* march CEILING, vox -- absolute safety cap on how
                             * far to probe; the real stop is the occupancy test,
                             * not this (def 8.0) */
    double occ_thresh;      /* a probe point within this of an existing mesh
                             * vertex is "in meshed territory", vox (def 0.9) */
    double local_r;         /* locality guard for source/actual-1-ring occupancy;
                             * non-topological nearby geometry still blocks
                             * (def 2.0 vox) */
    double gain_bias;       /* legacy target-selector knob; retained for CLI/API
                             * compatibility, ignored by the quilting pass */
    double tensor_weight;   /* legacy intensity-peak selector knob; retained for
                             * compatibility, ignored by the quilting pass */
    double tensor_radius;   /* legacy tangent-tensor radius (def 2) */
    double step;            /* probe step, vox (def 0.25) */
    double band_frac;       /* papyrus-band start as fraction of window (def 0.30) */
    double min_gain;        /* min intensity gain (raw units) to call fixable (def 25) */
    double anchor_frac;     /* min good-boundary fraction for FIXABLE-ISLAND (def 0.25) */
    int    region_cap;      /* region size above which flagged suspect (def 4000) */
    int    min_region;      /* reject dark regions smaller than this (verts) as
                             * speckle/false-positives -> back to GOOD (def 8; 0/1
                             * = keep all) */
    double pct_lo, pct_hi;  /* intensity window percentiles (def 1, 99) */
    int    dark_thresh;     /* post-stretch uint8 dark cutoff (def 64) -- the
                             * graph-cut data-term pivot dk */
    double dark_lambda;     /* dark graph-cut smoothness weight (def 20; 0 =>
                             * pure threshold, for A/B + parity) */
    double dark_sigma;      /* dark graph-cut contrast scale, raw units (def 30) */
    double normal_range;    /* legacy texture-sampling reach; geometric snap
                             * deliberately samples exactly at the vertex */
    int    normal_samples;  /* legacy texture-sampling taps (def 5) */

    /* --- region dilation (kills the halo rim around a corrected patch) --- */
    int    dilate_rings;    /* grow the dark/fixable set by this many mesh 1-rings
                             * before target-finding; a neighbour is absorbed only
                             * if it has a sample and cv < band_hi (never swallow
                             * clearly-bright papyrus). def 0 = off (legacy). */

    /* --- pass 1: quilting-energy repair target (SNAP_TARGET_*) --- */
    int    target_mode;     /* CLOSEST | QUILT_ARGMIN | QUILT_MRF (def MRF) */
    int    depth_bins;      /* odd K signed labels over [-reach,+reach]
                             * (def 33 -> 0.5 vox at reach 8) */
    double w_match;         /* quilting weight: cost += w_match*|s - mean(good nbr cv)|
                             * at a bright candidate depth (def 30) */
    double w_close;         /* small-deviation prior: cost += w_close*t (def 200) */
    double smooth_mu;       /* MRF depth smoothness weight (per 1-ring edge, QUILT_MRF
                             * only); truncated-linear on |d_u - d_v| (def 2000) */
    double smooth_tau;      /* smoothness truncation, vox (def 4) */

    int    verbose;
} SnapOpts;

void SnapOpts_default(SnapOpts *o);

/* target-selection mode: how the across-sheet stop depth is chosen.
 *   CLOSEST      -- legacy: stop at the FIRST depth clearing band + min_gain.
 *   QUILT_ARGMIN -- per-vertex argmin of a quilting/seam-continuity cost curve
 *                   along the ray (no smoothness coupling).
 *   QUILT_MRF    -- the cost curves become a GCO signed-depth labelling MRF
 *                   with clean-boundary zero pins and truncated-linear 1-ring
 *                   smoothness, preventing recto/verso jumps. This is default. */
enum { SNAP_TARGET_CLOSEST = 0, SNAP_TARGET_QUILT_ARGMIN = 1, SNAP_TARGET_QUILT_MRF = 2 };

/* per-vertex classes */
enum { SNAP_GOOD = 0, SNAP_FIXABLE = 1, SNAP_CRACK = 2 };
/* per-region classes */
enum { SNAPREG_FIXABLE = 0, SNAPREG_CRACK = 1, SNAPREG_ANCHORLESS = 2, SNAPREG_MIXED = 3,
       SNAPREG_REJECT = 4 };   /* below min_region -> verts returned to GOOD */

typedef struct SnapResult {
    /* per vertex */
    uint8_t *vclass;    /* [nv] SNAP_GOOD/FIXABLE/CRACK */
    float   *voff;      /* [nv] distance (vox) to the chosen free-space target */
    float   *vdir;      /* [nv*3] unit dir (z,y,x) of the winning ray; the snap
                         * target is verts[i] + voff[i]*vdir[i] (fixable) */
    float   *vgain;     /* [nv] intensity gain at the target (fixable) */
    float   *vtensor;   /* [nv] legacy diagnostic; zero in quilting mode */
    float   *vblock;    /* [nv] nearest distance the probe entered OTHER mesh
                         * (== reach if never blocked); small => wedged between
                         * sheets */
    float   *vcur;      /* [nv] direct raw intensity exactly at the vertex */
    float   *vkdefect;  /* [nv] angle-defect K_v (geometry-anomaly cross-signal) */
    int32_t *vregion;   /* [nv] dark-region id, or -1 */
    /* per region */
    uint8_t *rclass;    /* [nreg] SNAPREG_* */
    int32_t *rsize;     /* [nreg] vertex count */
    size_t   nreg;
    /* stats */
    size_t n_good, n_dark, n_fixable, n_crack;
    size_t n_bidir;          /* fixables where BOTH across-sheet dirs found a target */
    size_t n_blocked;        /* dark verts whose probe was blocked by other mesh */
    size_t n_reg_fixable, n_reg_crack, n_reg_anchorless, n_reg_mixed, n_reg_reject;
    size_t n_rejected;       /* dark verts returned to GOOD by the min-region filter */
    size_t n_overgrow;       /* dark verts with cv > band (cut pulled bright-ish; diag) */
    size_t max_region_size;
    double win_lo, win_hi;   /* intensity window used */
    double quilt_cost_mean;  /* pass-1 unary quilting+distance cost, fixables */
    double target_dist_mean; /* selected pass-1 target distance, vox */
    size_t n_quilt_fallback; /* GCO failure count (argmin fallback) */
    struct CubeTable *raw_table; /* arena-owned sampler reused by recto pass */
} SnapResult;

/* Run detection. verts [nv*3] (z,y,x), faces [nf*3] 0-based, uv [nv*2].
 * Outputs arena-allocated. Returns 0 on success, -1 on degenerate input. */
int SnapDetect_run(Arena_T arena,
                   const float *verts, size_t nv,
                   const int32_t *faces, size_t nf,
                   const float *uv,
                   const SnapOpts *opts, SnapResult *out);

/* In-process unit tests (off-surface->fixable, void->crack, anchorless,
 * grouping, clearance guard, degenerate). Returns 0 if all pass, else #fails. */
int Snap_selftest(void);

#endif
