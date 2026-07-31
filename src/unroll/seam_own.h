#ifndef SEAM_OWN_INCLUDED
#define SEAM_OWN_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "../common/raw_sample.h"
#include "piece_set.h"

/* ============================================================================
 * seam_own.h -- SeamOwn: strip-space ownership for the whole-grid piece set
 * (the v2 "quilting/multicut" pass, step2_overlap). Two defects, one output:
 *
 *   TRUE overlaps -- faces from different wraps (or a fold) painting the same
 *   (u,v): radius-gated pairs (|dr| > gate), grouped into regions, split into
 *   layers (connected components; lifted multicut only for intra-cube folds),
 *   ONE layer kept per region -- anchor-radius majority vote first, Efros
 *   min-error boundary energy (Quilt_boundary_energy) as tiebreak, largest as
 *   last resort. Losers dropped (they are mis-registered or occluded copies).
 *
 *   SEAM double-paint -- same-wrap faces of ADJACENT cubes re-painting the
 *   same px near the chart boundary (|dr| <= gate, 3D-near): per seam region
 *   a min-error-boundary DP picks the cut path through the least value
 *   mismatch; each cube keeps its side. A fill-safety pass restores dropped
 *   faces wherever the keep set would lose a covered pixel.
 *
 * Output is face_keep[nf] over the CURRENT PieceSet face order -- the raster
 * bakes step2 with it (dropped sole-cover px -> diagclass 6), then the driver
 * compacts via PieceSet_apply_facekeep for the later stages. uv is IMMUTABLE
 * here (no relocation): stage TIFs stay pixel-aligned.
 * ==========================================================================*/

typedef struct SeamOwnOpts {
    float  axis_point[3];    /* (z,y,x) umbilicus (default 0,3405,2878) */
    float  axis_dir[3];      /* (z,y,x) axis (default 1,0,0) */
    const char *raw_dir;     /* cubes_RAW; NULL = no texture energy AND no
                                seam DP (vote/largest only, seams left) */
    long   raw_chunk;        /* 128 */
    CubeTable *ct;           /* shared pre-built table (NULL = build own) */
    double pitch;            /* vox/turn (9.5) */
    double radius_gate;      /* <=0 => pitch/2 */
    double cell;             /* UV cell size, vox (2.0) */
    double seam_gate3d;      /* max 3D centroid gap for a SEAM pair (8.0) --
                                guards |dr| small but 3D-far constellations */
    double spiral_a, spiral_b;  /* diagnostic only (never decides) */
    double tex_du, tex_dv;   /* local raster px (1.0) */
    double normal_range;     /* RAW normal-max reach (2.0) */
    int    normal_samples;   /* taps (5) */
    size_t region_cap;       /* faces; bigger regions get vote-only (200000) */
    int    grid_cap;         /* local raster px cap per side (3000) */
    int    seam_cut;         /* 1 = run the seam DP pass (default on) */
    int    seam_halfband;    /* DP margin around the overlap band, px (4) */
    int    min_overlap_px;   /* below this a seam region is abut-only (8) */
    int    rehome;           /* 1 = relocate vertex-disjoint losing layers to
                                their physically-correct turn instead of
                                dropping (needs ps->phi + spiral_a/b) */
    int    rehome_min_faces; /* smaller losing layers stay dropped (50) */
    double rehome_frac_tol;  /* |dr/spiral_b - k| acceptance (0.35) */
    int    verbose;
} SeamOwnOpts;

void SeamOwnOpts_default(SeamOwnOpts *o);

/* Fill axis/pitch/spiral from <placed_dir>/placed_index.json (missing keys
 * leave defaults). Returns 0 if the file was read, -1 otherwise. */
int SeamOwn_read_index(const char *placed_dir, SeamOwnOpts *o);

typedef struct SeamOwnResult {
    uint8_t *face_keep;      /* [nf] 1 = keep (arena) */
    uint8_t *face_dec;       /* [nf] 0 none | 1 ovl-kept | 3 ovl-dropped |
                                4 seam-kept | 5 seam-dropped | 6 rehomed
                                (arena) */
    int32_t *face_region;    /* [nf] TRUE-overlap region id or -1 (arena) */

    /* pair classification */
    size_t n_true_pairs, n_seam_pairs, n_mystery_pairs;
    /* TRUE-overlap resolution */
    size_t n_regions, n_layers_total, max_region_faces;
    size_t n_vote_picks, n_energy_picks, n_largest_picks;
    size_t n_multicut_fallbacks, n_regions_skipped_cap, n_regions_skipped_grid;
    size_t n_kept, n_dropped;
    /* seam pass */
    size_t n_seam_regions, n_seams_abut_only, n_seams_skipped;
    size_t n_seam_dropped, n_seam_restored;
    /* rehome pass (find homes for dropped layers instead of losing them) */
    size_t n_rehome_layers, n_rehomed_faces;
    size_t n_rehome_dup, n_rehome_blocked, n_rehome_adjacent;
    size_t n_rehome_small, n_rehome_incoherent;
    /* coverage */
    size_t multi_cells_before, multi_cells_after;
    size_t prior_agree, prior_disagree;   /* spiral-prior diagnostic */
    double energy_mean;
    int    have_texture;
    double seconds;
} SeamOwnResult;

/* Run over the piece set. uv is immutable EXCEPT for rehomed layers: a
 * vertex-disjoint losing layer whose radius says it belongs k turns away is
 * relocated by the exact DeltaU(phi,k) (before the step2 bake, so every
 * stage-2+ TIF sees the joined positions; step1 stays the pre-join record).
 * All outputs arena-allocated. Returns 0 (empty ok), -1 on hard failure. */
int SeamOwn_run(Arena_T arena, PieceSet *ps,
                const SeamOwnOpts *opts, SeamOwnResult *out);

int SeamOwn_selftest(void);

#endif
