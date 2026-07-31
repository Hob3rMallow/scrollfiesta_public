#ifndef STRIP_METRICS_INCLUDED
#define STRIP_METRICS_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "piece_set.h"

/* ============================================================================
 * strip_metrics.h -- numeric QC of a baked strip, row-streamed (never loads
 * the 610 MB TIFs whole). Three signals the v2 stages move:
 *   multi_frac -- mean-blended overlap px / covered px   (step2 target)
 *   dark_frac  -- covered px darker than dark_thresh     (step3 target)
 *   seam_ratio -- mean |img[x+1]-img[x]| over SEAM columns divided by the
 *                 same statistic over baseline columns; ~1.0 = seams as
 *                 continuous as the sheet interior         (step4 target)
 * Seam columns = the u-interval endpoints of every (cube, winding-group)
 * patch -- exactly where the quilting gutters live.
 * ==========================================================================*/

typedef struct StripMetrics {
    double fill;             /* covered px / (W*H); covered = class 1 or 4 */
    double multi_frac;       /* class-4 px / covered px */
    double dark_frac;        /* covered px with img < dark_thresh / covered */
    double seam_dcol_mean;   /* mean |dcol| across seam columns */
    double base_dcol_mean;   /* same over every 97th column */
    double seam_ratio;       /* seam_dcol_mean / base_dcol_mean (0 if no data) */
    size_t n_seam_cols;      /* seam columns with at least one covered pair */
    size_t n_pairs_used;     /* covered (x,x+1) pairs across all seam columns */
    /* v-seam (z-cube boundary) rows -- the "cut into 4 strips" artifact */
    double vseam_gap_fill;   /* covered fraction of seam rows / covered
                                fraction of reference rows (y +- 8); ~1.0 =
                                the seam rows are as covered as the sheet */
    double vseam_ratio;      /* mean |img[y+1]-img[y-1]| across seam rows
                                (both covered) / same at reference rows */
    size_t n_vseam_rows;     /* seam rows inside the canvas */
    /* synthetic fill (step6 param_extend, diagclass 7) -- kept OUT of
       `covered` so fill/multi/dark/seam stay comparable across steps 1-5 */
    size_t synth_px;         /* class-7 px */
    double synth_frac;       /* synth / (covered + synth) */
    double fill_total;       /* (covered + synth) / (W*H) */
} StripMetrics;

/* Seam columns from the piece set: per-(cube,gid) u-interval endpoints over
 * vertices referenced by kept faces, mapped to strip px (px = round(u/du)-X0
 * with X0 = round(u0/du)), sorted, deduped, canvas ends dropped. Arena-
 * allocates *out_cols. Returns 0 (possibly with *out_n == 0). */
int StripMetrics_seam_cols(Arena_T arena, const PieceSet *ps,
                           double u0, double du, size_t W,
                           int32_t **out_cols, size_t *out_n);

/* Row-streamed metrics over a rawtex + diagclass pair. seam_cols may be NULL
 * (seam stats then zero). out_col_dcol (optional) receives the per-seam-column
 * mean |dcol| (arena, [ncols]; -1.0 where fewer than H/8 covered pairs, so
 * top-K zoom picks land on substantial seams, not sliver edges). Returns 0 on
 * success, -1 on IO/shape error. */
int StripMetrics_compute(Arena_T arena, const char *rawtex_tif,
                         const char *diagclass_tif,
                         const int32_t *seam_cols, size_t ncols,
                         int dark_thresh, double **out_col_dcol,
                         StripMetrics *out);

/* v-seam row stats via targeted row reads (call after _compute; cheap).
 * Seam rows are world-z multiples of chunk mapped into the canvas:
 * y = k*chunk - v0 for integer k with 8 < y < H-9. chunk <= 0 disables. */
int StripMetrics_vseams(Arena_T arena, const char *rawtex_tif,
                        const char *diagclass_tif,
                        double v0, long chunk, StripMetrics *inout);

int StripMetrics_selftest(void);

#endif
