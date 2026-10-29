#ifndef CLOSEUP_REPAIR_INCLUDED
#define CLOSEUP_REPAIR_INCLUDED
/* ============================================================================
 * closeup_repair -- duplicate and wrong-side geometry of a close-up surface
 * (2026-10-12, built on PHerc0211 close-up surfaces).
 *
 * An ink_closeup surface mixes OBSERVED vertices (on a source chart) with FILL
 * (a quadratic extrapolated from the observed nodes to complete gaps and the
 * apron). The 10-11 revision's "unexpected pieces of letters" came from three
 * defects, measured on the two shipped surfaces by their vertex states:
 *   DUPLICATE   fill laid over papyrus the surface covers elsewhere: a vertex
 *               within CLOSEUP_REPAIR_DUP_VOX in 3-D of another that is more
 *               than CLOSEUP_REPAIR_DUP_PX away on the surface (cand5: 70% of
 *               the pairs fill over observed, 26% fill over fill, 4% observed)
 *   WRONG FILL  fill whose inward normal points away from the roll's axis
 *               (cand5: 87% of the wrong-side vertices)
 *   FLAP        a folded flap of observed papyrus read from its back: a
 *               connected set of observed vertices facing away from the axis
 *               (cand1: 71% of the large wrong-side components)
 * plus ORPHAN fill left without an observed vertex in its component.
 * Inward is the bulk's: when most observed papyrus faces away from the axis
 * the normals are read with the other sign. CloseupRepair_mark only labels;
 * the tool keeps the rest and routes flaps (observed papyrus, never deleted)
 * to an extras surface with provenance.
 * ==========================================================================*/
#include <stddef.h>
#include <stdint.h>

enum {
    CLOSEUP_REPAIR_KEEP = 0,
    CLOSEUP_REPAIR_DUPLICATE = 1,
    CLOSEUP_REPAIR_WRONG_FILL = 2,
    CLOSEUP_REPAIR_FLAP = 3,
    CLOSEUP_REPAIR_ORPHAN = 4
};

typedef struct CloseupRepairStats {
    size_t verts, observed, fill, undefined;           /* undefined: no direction to the axis */
    size_t duplicate, wrong_fill, flap, orphan;        /* vertices per reason */
    size_t pairs_fill_observed, pairs_fill_fill;       /* duplicate pairs that removed a fill vertex */
    size_t pairs_observed_observed;                    /* doubled observed papyrus (kept; the source's own) */
    size_t flap_components, flap_specks;               /* flaps routed; smaller wrong-side observed sets kept */
    size_t orphan_components;
    int inward_flipped;                                /* the normals were read with the other sign */
    double observed_facing_axis;                       /* share of observed vertices facing the axis */
    double seconds;
} CloseupRepairStats;

/* verts [nv*3] (z,y,x), uv [nv*2] (surface px), faces [nf*3], state [nv] (1 observed, 2 fill, 3 observed on an
 * adopted chart, 0 none), nrm [nv*3] (z,y,x), inward [nv*3] the unit direction from the vertex to the axis
 * (z,y,x; all 0 = undefined). reason [nv] receives CLOSEUP_REPAIR_*. Returns 0, or -1 on bad input. */
int CloseupRepair_mark(const float *verts, const double *uv, size_t nv, const int32_t *faces, size_t nf,
                       const uint8_t *state, const float *nrm, const float *inward, uint8_t *reason,
                       CloseupRepairStats *st);

int CloseupRepair_selftest(void);

#endif
