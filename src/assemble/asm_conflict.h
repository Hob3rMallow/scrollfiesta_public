#ifndef ASM_CONFLICT_INCLUDED
#define ASM_CONFLICT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "asm_types.h"
#include "asm_pose.h"

/* ============================================================================
 * asm_conflict.h -- the contradiction audit of a layout, and the cleaning
 * rounds it drives.
 *
 * AUDIT: every placed chart's footprint is rasterized into cells of
 * ASM_CELL_VOX in its component's global frame.  A cell claimed by two
 * charts of the same component whose 3-D points are separated ALONG THE
 * SHEET NORMAL is a contradiction: the layout claims two different layers
 * are one place.  The weight is continuous in the local layer distance
 * (smoothstep between ASM_CONTRA_LO and ASM_CONTRA_HI of it), so no
 * threshold seam is manufactured.  The score is covered area minus
 * contradiction area (Stevens' objective with exchange rate 1).
 *
 * CLEANING: contradictions are attributed to relations by witness paths
 * (the shortest relation path between the two charts) weighted by the
 * contradiction mass; the most implicated relations are dropped in a
 * batch, the poses re-solved, and the batch kept only if the score did not
 * fall.  Charts that keep contradicting more than half of their own area
 * after their joins are gone are routed to extras.  Cuts (asm_cut) come
 * between those two steps once built.
 * ==========================================================================*/

typedef struct AsmConflictOpts {
    double cell;           /* raster cell, vox */
    double contra_lo;      /* smoothstep bounds as fractions of the local layer distance */
    double contra_hi;
    int    rounds;         /* cleaning rounds */
    int    max_path;       /* witness path length cap (relations) */
    double drop_frac;      /* fraction of flagged relations dropped per round */
    double extras_ratio;   /* chart contradiction mass / own area above this -> extras */
    const char *diag_pairs_csv;  /* nullable: the audit writes one row per contradicting chart pair here
                                  * (mass, cells, normal separation, witness path length, seam / layer evidence) */
} AsmConflictOpts;

typedef struct AsmConflictCell {
    int32_t component;
    int32_t cx, cy;
    float   w;             /* contradiction weight of the worst pair in the cell */
    float   dn;            /* its normal separation, vox */
} AsmConflictCell;

typedef struct AsmConflictReport {
    size_t cells;              /* distinct claimed cells (all components) */
    size_t cells_multi;        /* cells with >= 2 claimants of one component */
    size_t cells_contra;       /* cells with weight > 0.5 */
    double covered_area;       /* cells * cell^2 */
    double contra_area;        /* sum of cell weights * cell^2 */
    double score;              /* covered - contra */
    size_t pairs;              /* contradicting chart pairs */
    AsmConflictCell *contra_cells;   /* arena array for the report image */
    size_t n_contra_cells;
} AsmConflictReport;

typedef struct AsmConflictStats {
    int    rounds;
    size_t rels_dropped;
    size_t rels_reverted;
    size_t charts_extras;
    size_t splits_refused;        /* batches refused by the split guard */
    double split_area_refused;    /* the area those batches would have separated, vox^2 */
    double score_first, score_last;
    double contra_first, contra_last;
    double covered_first, covered_last;
    size_t pairs_first, pairs_last;
    double sec;
} AsmConflictStats;

void AsmConflict_default_opts(AsmConflictOpts *o);

/* Audit the current poses.  Contradiction masses are written into the
 * relations' flag_mass (witness attribution) and the per-chart array
 * out_chart_mass ([n_charts], nullable). */
/* INCREMENTAL GRID: the placed material's raster kept between placement
 * attempts (rebuilding it per attempt was 90% of the 21x5x5's 43 minutes).
 * Cells hash (cx, cy) -> a list of claims; a candidate's charts, in their
 * tentative poses, are rasterized and probed against it; on acceptance they
 * are inserted.  Same contradiction weight as the full audit. */
typedef struct AsmConflictGrid AsmConflictGrid;
AsmConflictGrid *AsmConflictGrid_new(Arena_T arena, const AsmConflictOpts *o, const AsmRun *run);
void AsmConflictGrid_insert(AsmConflictGrid *g, const AsmChart *c);
/* contradiction mass (vox^2) of the charts' cells against the grid; *partner = the
 * grid chart carrying the most of it, *partner_mass its share */
double AsmConflictGrid_probe(AsmConflictGrid *g, const AsmChart *const *charts, size_t n,
                             int32_t *partner, double *partner_mass);
/* Probe a relocation while ignoring only the old samples of its charts.
 * excluded is [run->n_charts], or NULL for the ordinary insertion probe. */
double AsmConflictGrid_probe_excluding(AsmConflictGrid *g, const AsmChart *const *charts, size_t n,
                                       const uint8_t *excluded, int32_t *partner, double *partner_mass);
size_t AsmConflictGrid_cells(const AsmConflictGrid *g);
/* Cells holding two claims of DIFFERENT charts a layer apart in 3-D (the same
 * weight as the audit, > 0.5): the stacked share of a placed layout, component-
 * blind.  *cells (nullable) receives the cell count. */
size_t AsmConflictGrid_stacked(const AsmConflictGrid *g, size_t *cells);
void AsmConflictGrid_dispose(AsmConflictGrid *g);   /* frees its scratch arena */

int AsmConflict_audit(AsmRun *run, const AsmConflictOpts *o,
                      double *out_chart_mass, AsmConflictReport *rep);

/* Layout-confirmed joins after placement (see ASM_JOIN_* in
 * pipeline_constants.h): every pair of placed charts of different lineages
 * whose adjacent / co-located audit cells agree in 3-D over a seam's length
 * without contradiction gets a relation flagged ASM_REL_LAYOUT taken from
 * the placed poses.  `diag_csv` (nullable) receives one row per tested pair. */
typedef struct AsmJoinStats {
    size_t pairs_tested, pairs_joined, pairs_refused_disagree, pairs_refused_short, pairs_refused_veto;
    size_t cells_agree, cells_disagree;
    size_t lineages_before, lineages_after;   /* over placed charts in the layout */
    size_t charts_shifted;      /* charts moved by the in-plane alignment (ASM_JOIN_ALIGN) */
    double shift_max;           /* the largest lineage shift applied, vox */
    double rot_max;             /* the largest lineage rotation applied, rad */
    double sec;
} AsmJoinStats;
int AsmConflict_confirm_joins(AsmRun *run, const AsmConflictOpts *o, const char *diag_csv, AsmJoinStats *st);

/* Post-placement seam re-solve (ASM_REFINE_* in pipeline_constants.h): readmits the seams the
 * placed layout confirms and re-solves the placed charts outside the largest lineage over the
 * joins with that lineage fixed; kept only if the contradiction holds.  `diag_csv` (nullable)
 * receives one row per candidate seam. */
typedef struct AsmRefineStats {
    size_t cand_join, cand_by_kind[5];        /* kinds: 1 accepted-but-dropped, 2 switched, 3 weak, 4 short */
    size_t readmitted, readmit_by_kind[5], readmit_switched, unreadmitted;
    size_t refused_gap, refused_rot, refused_parity, refused_nocorr, internal_fixed;
    size_t charts_placed, charts_fixed, charts_free, charts_pinned;
    int    iters; size_t switched; int rounds; int kept; int reverted;
    size_t moved; double moved_area, shift_p50, shift_p90, shift_max, rot_max;
    double contra_before, contra_after, covered_before, covered_after;
    size_t lineages_before, lineages_after;
    double sec;
} AsmRefineStats;
int AsmConflict_refine_seams(AsmRun *run, const AsmConflictOpts *o, const AsmPoseOpts *po, const char *diag_csv, AsmRefineStats *st);

/* Cleaning rounds: audit, drop the most implicated joins, re-solve, keep or
 * revert; then route persistently contradicting charts to extras.  The poses
 * are re-solved with `pose_opts` after every accepted change. */
int AsmConflict_clean(AsmRun *run, const AsmConflictOpts *o,
                      const AsmPoseOpts *pose_opts,
                      AsmConflictReport *final_report, AsmConflictStats *st);

int AsmConflict_selftest(void);

#endif
