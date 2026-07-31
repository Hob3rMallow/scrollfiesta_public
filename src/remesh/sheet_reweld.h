#ifndef SHEET_REWELD_INCLUDED
#define SHEET_REWELD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "ball_pivot.h"   /* BpaBridgeGate */

/* ============================================================================
 * sheet_reweld.h -- Phase-2 sheet-correspondence permissive re-weld.
 *
 * The primary (phase-1) seam bridge is conservative: the winding gate rejects
 * any cross-wrap connection, which correctly stops a divot from welding UP to
 * the next wrap -- but also leaves legit same-sheet closures (the divots
 * themselves, large seam-disagreement regions) as holes.
 *
 * Phase 2 recovers them SAFELY. The per-cube sheets are clean and
 * well-separated; the only ambiguity is which sheet continues into which
 * across a seam. If we confidently match sheet A_i <-> B_j, a permissive weld
 * restricted to ONLY those two sheets' vertices cannot connect to a wrong
 * sheet -- the ball physically can't reach a third. So:
 *   1. label the PRE-weld mesh into per-cube sheets (connected components),
 *   2. confirm 1:1 correspondence across each seam from BOTH the phase-1
 *      bridge's own confident connections AND geometric boundary overlap,
 *      accepting only mutually-unambiguous matches (one sheet and ONLY one),
 *   3. re-weld each confirmed pair with the winding gate OFF and a wider ball,
 *      on a cloud restricted (via SeamWeld_bridge's want_mask) to those two
 *      sheets.
 * ==========================================================================*/

typedef struct {
    double overlap_r;    /* geometric boundary-overlap match radius (vox) */
    double vote_weight;  /* weight of a phase-1 bridge vote vs an overlap match */
    double min_score;    /* min combined score to confirm a pair */
    double margin_frac;  /* runner-up score must be < margin_frac * best */
    float  rho_max;      /* wider ball cap for the permissive re-weld (vox) */
    float  band;         /* seam band, match the phase-1 bridge (vox) */
} SheetReweldParams;

typedef struct {
    size_t n_sheets;       /* connected components (per-cube sheets) */
    size_t n_candidates;   /* cross-seam sheet pairs with any evidence */
    size_t n_pairs;        /* confirmed 1:1 pairs re-welded */
    size_t n_faces_added;  /* phase-2 bridge faces merged */
} SheetReweldStats;

/* Fill defaults from pipeline_constants.h. */
void SheetReweld_default_params(SheetReweldParams *p);

/* Label the PRE-weld concatenated mesh (cubes merely concatenated, no
 * cross-cube faces yet) into per-cube sheets = connected components. Writes an
 * arena-allocated vert_sheet[nv] (sheet id >= 0, or -1 for an unused vertex)
 * and the sheet count. Call BEFORE the phase-1 seam bridge. Returns 0. */
int SheetReweld_label(Arena_T arena, const int32_t *faces, size_t nf, size_t nv,
                      int32_t **out_vert_sheet, size_t *out_n_sheets);

/* Phase-2 re-weld, run AFTER the phase-1 bridge. `n_phase1_bridge` = the number
 * of trailing faces in `faces` that are phase-1 bridge faces (used for the
 * correspondence votes; grid_weld knows this from SeamWeld_bridge). `vert_sheet`
 * / `n_sheets` come from SheetReweld_label (pre-weld). `gate` supplies the
 * umbilicus + pitch for the winding-phase geometry -- the pair welds themselves
 * run with the gate OFF (cloud restriction is the safety). On success
 * *out_faces / *out_nf hold the combined (input + phase-2 bridge) face list,
 * arena-allocated. Returns 0. */
int SheetReweld_process(Arena_T arena,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        size_t n_phase1_bridge,
                        const int32_t *vert_sheet, size_t n_sheets,
                        float cube_size,
                        const BpaBridgeGate *gate,
                        const SheetReweldParams *params,
                        int32_t **out_faces, size_t *out_nf,
                        SheetReweldStats *stats);

/* In-process unit tests. Returns 0 if all pass, else the failure count. */
int SheetReweld_selftest(void);

#endif
