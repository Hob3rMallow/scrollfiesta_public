#ifndef TIFXYZ_RECTO_INCLUDED
#define TIFXYZ_RECTO_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/pipeline_constants.h"
#include "../common/raw_sample.h"

/* ============================================================================
 * tifxyz_recto.h -- snap a tifxyz node grid onto the papyrus RECTO edge read
 * in the RAW CT (2026-10-09/10, built on PHerc0211 close-ups).
 *
 * A close-up surface is fitted to mesh observations, so where it sits INSIDE
 * the papyrus is whatever the mesher left: on cand5's 21-layer render the
 * surface lies p50 0.9 vox outward of the recto edge with a p5..p95 spread of
 * -4.0..+5.1 vox.  The ink lies on the recto.
 *
 * THE PRIOR (TifxyzRecto_prior): before anything moves, long rays are shot
 * along the node normals (every TXR_PRIOR_STRIDE-th node, +-TXR_PRIOR_RANGE_VOX)
 * and read for the wrap's structure: the band the node sits in (the crest
 * nearest to it), its inward (recto) and outward (verso) edges, the dips inside
 * it (the seam between the two fibre layers) and the gaps around it (air, the
 * darkest level), the next wraps' crests, the descent slopes and edge widths.
 * Their distributions set, PER SCROLL and per surface, what used to be global
 * constants: the search range and the move cap (from the wrap spacing and the
 * own recto's reach), the band gate (the depth that separates seam dips from
 * gaps), the edge-strength gates, the profile smoothing (the scanner's edge
 * width) and the band thickness bounds.  The TXR_* constants remain the
 * defaults and the safety bounds, and the record says which values ran.
 *
 * THE SNAP (TifxyzRecto_run): per node, the smoothed CT profile along the
 * node's INWARD normal (the direction vc_render's layer index grows along with
 * --flip-normals: normalize(dP/dcol x dP/drow) in (z, y, x) order,
 * rawtex_tifxyz.h).  The node belongs to the band whose crest is NEAREST to it
 * (a band = a local maximum standing the gate above the darkest sample); that
 * band's recto = the steepest descent to the first trough inward that stands
 * the gate below the crest (shallower dips are walked over).  A node whose
 * band is in view but whose recto is out of reach abstains and follows its
 * neighbours; it never takes another band's edge.  The per-node offsets go
 * through a point-symmetric median (the fibre weave, +-2 vox streaks 8-12 px
 * wide, is texture the detector learned, not drift), and one global scalar
 * field minimizes sum w_i (d_i - t_i)^2 + mu sum (d_i - d_j)^2 on the
 * 4-connected grid by Jacobi-preconditioned CG in double, least squares then
 * Cauchy-reweighted rounds; the move is capped per iteration and in total and
 * measured again along the moved normals.  No threshold treats cells
 * differently (the 2026-09-03 rule).
 * ==========================================================================*/

typedef struct TifxyzRectoOpts {
    int midline;    /* 0: the recto edge; 1: the band's midline */
    int inset;      /* 1: the recto edge minus TXR_RECTO_INSET_VOX (the surface stops inside the papyrus) */
    int smooth;     /* 1: the TXR_MEDIAN_HALF_SMOOTH window at stride TXR_MEDIAN_STRIDE_SMOOTH (the drift only) */
    int verbose;    /* stderr progress lines */
    size_t count_c0, count_c1;   /* the columns whose nodes the statistics count (count_c1 == 0: every column) --
                                  * a window's core: its halo is measured and moved but counted by its own window */
} TifxyzRectoOpts;

typedef struct TifxyzRectoHist {
    size_t n;
    double p5, p25, p50, p75, p95;
    double mean_abs, max_abs;
    double share_over_1p5, share_over_3;    /* |v| > 1.5 vox / > 3 vox */
} TifxyzRectoHist;

/* A measured quantity's distribution over the prior's rays. */
typedef struct TifxyzRectoQ {
    size_t n;
    double p01, p05, p25, p50, p75, p95, p99;
} TifxyzRectoQ;

typedef struct TifxyzRectoPrior {
    size_t rays, rays_with_band, rays_with_recto, rays_with_next, rays_with_seam;
    double range_vox, stride;
    TifxyzRectoQ recto_offset;      /* the own recto's offset from the node along inward, vox */
    TifxyzRectoQ spacing;           /* own crest to the neighbouring wraps' crests, vox (both sides) */
    TifxyzRectoQ thickness;         /* verso edge to recto edge, vox */
    TifxyzRectoQ gap_width;         /* recto edge to the next wrap's verso edge, vox */
    TifxyzRectoQ gap_depth;         /* the lower adjacent crest minus the gap trough, grey */
    TifxyzRectoQ seam_depth;        /* the same for dips INSIDE a band (the fibre-layer seam) */
    TifxyzRectoQ slope;             /* the recto descent's steepest slope, grey / vox */
    TifxyzRectoQ edge_width;        /* the descent's 10%..90% distance, vox */
    double gate_error;              /* share of dips the chosen gate misclassifies (seams over it + gaps under it) */
    /* the scroll's levels: a gap is a trough that falls to the AIR; the air level is the split of all troughs'
     * levels into a darker and a brighter class (Otsu), the band level the typical crest; a crest counts as a band
     * from the midpoint between them */
    TifxyzRectoQ trough_level, crest_level;
    double air_level, band_level, band_min_level, level_separability;   /* the split's between-class share of variance */
    /* every trough's width at half depth, by its level class (gap: fallen to the air level; seam: not), and the
     * width that splits all troughs into narrow and wide (Otsu): a seam between the fibre layers is ~1 vox wide, a gap
     * between wraps 3-4 (PHerc0211) -- the discriminator the levels alone cannot give on a scroll whose thin gaps do
     * not all reach the air level */
    TifxyzRectoQ gap_trough_width, seam_trough_width;
    double trough_width_split, trough_width_separability, trough_width_error;
    /* derived parameters (the snap runs on these) */
    double search_vox, max_total_vox, min_contrast, min_edge, edge_full, sigma, thick_min, thick_max;
    int    crowded;                 /* the own recto's reach exceeds what the wrap spacing allows */
} TifxyzRectoPrior;

/* Shoot the rays and derive the parameters.  pos as for the run, ct prewarmed.
 * stride <= 0 uses TXR_PRIOR_STRIDE.  Returns 0 on success; a grid with no
 * measurable band leaves the defaults in place and returns 0 with rays
 * counted (the caller reports it). */
int TifxyzRecto_prior(Arena_T arena, CubeTable *ct, const float *pos, size_t W, size_t H, int stride,
                      TifxyzRectoPrior *out);

typedef struct TifxyzRectoStats {
    size_t n_nodes, n_valid, n_dir;                 /* nodes; valid (not -1); with an inward direction */
    size_t n_supported_first, n_in_band_first, n_in_gap_first;   /* edge inward / outward of the node */
    size_t n_supported_last, n_in_band_last, n_in_gap_last;
    size_t n_beyond_first;                          /* abstained: the node's band has its recto out of reach */
    int    outer_iters;
    double mean_abs_step[TXR_OUTER_ITERS];
    int    pcg_iters[TXR_OUTER_ITERS];
    double pcg_resid[TXR_OUTER_ITERS];              /* relative residual reached */
    size_t n_capped_total;                          /* nodes that hit the total cap */
    size_t n_filtered_first, n_filtered_last;       /* nodes with a median-filtered target */
    TifxyzRectoHist target_first;                   /* recto offset along inward before any move, per node */
    TifxyzRectoHist target_last;                    /* ... measured again after the last move */
    TifxyzRectoHist target_first_med;               /* the same through the median window */
    TifxyzRectoHist target_last_med;
    TifxyzRectoHist disp;                           /* total displacement along inward */
    double slope_p99, slope_max;                    /* |disp_i - disp_j| over valid grid edges */
    double profile_first[2 * TXR_PROFILE_HALF + 1]; /* mean CT at -HALF..+HALF vox along inward */
    double profile_last[2 * TXR_PROFILE_HALF + 1];
    size_t profile_n;                               /* nodes behind each mean */
    /* the parameters the run used (the prior's, or the constants) */
    double prm_search, prm_max_total, prm_min_contrast, prm_min_edge, prm_edge_full, prm_sigma, prm_thick_min,
           prm_thick_max, prm_air_level, prm_band_min_level;
    int    prm_levels;                              /* 1: the band ends where the CT falls to the air level */
    double sec_measure, sec_solve, sec_total;
    int    windows;                                 /* TifxyzRecto_run_windows: the windows run */
} TifxyzRectoStats;

#define TXR_PROFILE_N (2 * TXR_PROFILE_HALF + 1)

/* pos: plane-major (z, y, x), W * H floats each, -1 on invalid nodes; moved
 * in place.  ct must be prewarmed (concurrent reads).  prior: the measured
 * parameters (nullable: the constants).  The five rasters (nullable, W * H
 * each) receive the diagnostics, NaN where undefined: the per-node recto
 * offset before the move and after it, the same median-filtered, and the
 * displacement.  Empty or all-invalid input is a clean no-op.  0 on success. */
int TifxyzRecto_run(Arena_T arena, CubeTable *ct, float *pos, size_t W, size_t H,
                    const TifxyzRectoOpts *opts, const TifxyzRectoPrior *prior, TifxyzRectoStats *st,
                    float *target_first, float *target_last, float *target_first_med,
                    float *target_last_med, float *disp);

/* After each window of TifxyzRecto_run_windows: its rasters (ww x H, NaN
 * undefined; the five of TifxyzRecto_run), the core columns [c0, c1) of the
 * window that the result keeps, the window's first column lc0 in the input
 * grid and its statistics (counted over the core). */
typedef void (*TifxyzRectoWindowFn)(void *ctx, size_t k, size_t n_windows, size_t lc0, size_t ww, size_t H,
                                    size_t c0, size_t c1, const TifxyzRectoStats *st, const float *target_first,
                                    const float *target_last, const float *target_first_med,
                                    const float *target_last_med, const float *disp);

/* TifxyzRecto_run over the columns [col0, col1) of a W x H grid in windows:
 * cores of window_cols columns (0: one core) tile [col0, col1); each runs
 * with `halo` columns of the input either side (clipped to the grid) and
 * keeps its core.  Every window reads pos_in, never a neighbour's result, so
 * the outcome does not depend on the order.  pos_out (3 planes of W x H)
 * receives pos_in with the cores replaced; windows whose core holds no valid
 * node are skipped.  st: the cores' statistics combined (counts summed, the
 * distributions from the core rasters, profiles weighted, solver counts
 * summed); fn (nullable) sees every window.  ct prewarmed.  0 on success. */
int TifxyzRecto_run_windows(Arena_T arena, CubeTable *ct, const float *pos_in, float *pos_out, size_t W, size_t H,
                            size_t col0, size_t col1, size_t window_cols, size_t halo,
                            const TifxyzRectoOpts *opts, const TifxyzRectoPrior *prior, TifxyzRectoStats *st,
                            TifxyzRectoWindowFn fn, void *fn_ctx);

/* Synthetic slabs (two wraps): the prior reads their spacing, thickness and
 * edges; a bent grid converges onto its band's inward edge, flat and
 * near-idempotent; a patch in the gap outside its band follows its
 * neighbourhood inward and never takes the outer wrap's edge; a hole stays
 * invalid; midline, inset and the wide window behave; windows with a halo
 * give the whole run's result and counts, a column subset leaves the rest
 * untouched; the selective prewarm loads only the reached chunk.  Returns the
 * number of failures (0 = pass). */
int TifxyzRecto_selftest(void);

#endif
