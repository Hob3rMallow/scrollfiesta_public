#ifndef ASM_DISCOVER_INCLUDED
#define ASM_DISCOVER_INCLUDED

#include <stddef.h>
#include "asm_types.h"

/* ============================================================================
 * asm_discover.h -- the second pass over a layout (2026-09-16).
 *
 * Pass one is the placement we have.  Pass two looks at every chart the
 * placement left out and asks the placed charts it shares a seam with where
 * it belongs: each seam relation, composed with its placed partner's pose,
 * implies a rigid pose for the leftover.  A leftover whose placed partners
 * AGREE (every implied pose moves every vertex by at most
 * ASM_DISCOVER_AGREE_VOX of every other) on enough SIDES (distinct quadrants
 * of the leftover's own frame the partners sit on) is placed at the mean
 * pose -- even where the ribbon already carries material, PROVIDED that
 * material is the same surface (the overlapped registered vertices lie within
 * ASM_DISCOVER_SAME_SURFACE_VOX in 3-D: a delamination or a duplicate
 * observation).  Overlapped material farther away in 3-D is another wrap and
 * refuses the placement.  The leftover's rim must also be consistent with
 * the placed material around it (boundary vertices whose nearest placed
 * vertex in UV is also near in 3-D).  Rounds repeat until nothing changes,
 * so a chart admitted in one round is a partner in the next (discovery).
 *
 * At a crop face the sides beyond the face do not exist, so a chart whose
 * bounding box touches the pile's extent needs one agreeing side.
 *
 * Every examined chart is a row of <out>/stage5_discovery.csv, with the
 * reason it was or was not placed; <out>/stage5_discovered.png shows the
 * layout afterwards.  Before the rounds, the same construction is run on
 * the PLACED charts as a control: their partner-implied poses must land on
 * their actual poses (the log line reports p50/p90/max).
 * ==========================================================================*/

typedef struct AsmDiscoverOpts {
    int enabled;        /* 0 = report only, place nothing */
    int min_sides;      /* agreeing sides an interior chart needs (ASM_DISCOVER_MIN_SIDES) */
    int soft_min_sides; /* sides gate-rejected (soft) evidence alone needs (ASM_DISCOVER_SOFT_MIN_SIDES) */
    int rounds;         /* at most this many discovery rounds (ASM_DISCOVER_ROUNDS) */
    int selected_obligations; /* opt-in: require only evidence in the chosen consensus; preserve prior requirements */
} AsmDiscoverOpts;

typedef struct AsmDiscoverStats {
    size_t examined;            /* unplaced charts with a pose in the layout */
    size_t with_partner;        /* ... that share a seam with a placed chart */
    size_t admitted;            /* placed by this pass */
    double admitted_area;       /* vox^2 */
    size_t refused_disagree, refused_sides, refused_rim, refused_wrap;
    size_t admitted_overlapping; /* admitted although the ribbon carried same-surface material there */
    int    rounds;
    /* control on the placed charts: displacement of the partner-implied pose from the actual one */
    size_t control_poses, control_beyond;    /* poses tested, poses beyond ASM_DISCOVER_AGREE_VOX */
    double control_p50, control_p90, control_max;
    double sec;
} AsmDiscoverStats;

void AsmDiscover_default_opts(AsmDiscoverOpts *o);

/* Runs the pass on `run` in place (poses, flags, placement fields and the
 * continuity obligations of the admitted charts change; nothing else does).
 * out_dir may be NULL (no files).  Returns 0, or -1 on an I/O failure. */
int AsmDiscover_run(AsmRun *run, const AsmDiscoverOpts *o, const char *out_dir, AsmDiscoverStats *st);

int AsmDiscover_selftest(void);

#endif
