#ifndef PRED_REJECT_INCLUDED
#define PRED_REJECT_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * PredReject -- pre-Step-0 input validation: detect garbage nnUNet prediction
 * cubes (a big solid white rectangle / slab) so they are never meshed.
 *
 * Motivation: some cubes in a grid are not surface predictions at all -- the
 * nnUNet output is a solid filled box, which the mesh pipeline turns into piles
 * of junk fragments that pollute the weld. Gross fill fraction CANNOT separate
 * these from real cubes (a 27.7%-filled garbage slab looks identical to a
 * 26.6%-filled real dense tangle). The discriminator is SOLIDITY / THICKNESS:
 *
 *   garbage = a SOLID filled region (thick 3D interior, rectangular per frame)
 *   real    = THIN sheets (1-3 voxels thick) even when densely packed
 *
 * Two signals, AND-combined (conservative -- only reject unambiguous slabs):
 *   1. 3D erosion-survival: erode the foreground by THICK_ERODE_R voxels;
 *      thin sheets vanish, solid slabs keep a large interior.
 *   2. per-frame solid rectangle: on each axis, the largest 2D connected
 *      component of a slice fills its bounding box (a filled rectangle) over a
 *      long run of consecutive slices (a persistent slab, not a 1-3-frame
 *      flat sheet seen face-on).
 *
 * Thresholds live in pipeline_constants.h (GARBAGE_* / *RECT* group).
 */

typedef struct {
    size_t  n_vox;            /* total voxels D*H*W                          */
    size_t  fg_count;         /* foreground voxels                           */
    double  fill_frac;        /* fg_count / n_vox                            */

    /* 3D thickness (erosion-survival) signal */
    size_t  interior_count;   /* FG voxels surviving THICK_ERODE_R erosions  */
    double  interior_frac;    /* interior_count / fg_count                   */
    int     max_thickness;    /* max L1 erosion depth over FG (capped)       */

    /* per-frame solid-rectangle signal (reported for the strongest axis) */
    int     rect_axis;        /* 0=z, 1=y, 2=x                               */
    double  rect_frac;        /* fraction of slices that are rectangle frames*/
    int     max_rect_run;     /* longest run of consecutive rectangle frames */

    int     verdict;          /* 1 = reject (garbage), 0 = keep              */
    const char *reason;       /* short human-readable verdict reason         */
} PredRejectStats;

/*
 * Decide whether a binary prediction cube is garbage (solid slab) and should be
 * skipped. vol is uint8[D*H*W]; any nonzero voxel is foreground (thresholded
 * internally -- vol is never modified). Returns 1 to REJECT, 0 to KEEP. If out
 * is non-NULL it receives all computed metrics. Uses arena scratch bounded by
 * Arena_save/Arena_restore (nothing leaks into the caller's arena).
 */
int PredReject_is_garbage(Arena_T arena, const uint8_t *vol,
                          int D, int H, int W, PredRejectStats *out);

/* Self-test on synthetic volumes. Returns 0 on pass, nonzero on failure. */
int PredReject_selftest(void);

#endif
