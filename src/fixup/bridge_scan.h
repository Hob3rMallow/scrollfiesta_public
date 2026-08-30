#ifndef BRIDGE_SCAN_INCLUDED
#define BRIDGE_SCAN_INCLUDED

#include <stdint.h>
#include "../common/arena.h"

/* bridge_scan — thin-neck prediction-weld detector for the 2D gap fixup.
 *
 * The inverse pathology of a gap: the model paints a 1-2px bridge fusing two
 * wraps that run 2-3px apart near the core ("lumpy joined to bumpy by one
 * pixel"). Downstream this becomes a wrap merger the mesh pipeline has to
 * fight. This pass finds them per plane:
 *
 *   thin neck   a run of skeleton pixels whose chamfer half-width is
 *               <= FIXUP_BRIDGE_MAX_WIDTH, run length <= FIXUP_BRIDGE_MAX_LEN,
 *               whose BOTH ends open into material >= FIXUP_BRIDGE_END_WIDTH
 *               thick (a neck BETWEEN bodies, not a tapering tip);
 *   certified   radial gate armed AND the neck moves radially:
 *               |r(endA)-r(endB)| >= FIXUP_BRIDGE_MIN_DR about the umbilicus
 *               and |dot(neck_dir, radial_unit)| >= FIXUP_BRIDGE_RADIAL_DOT.
 *               An along-sheet thin section is NOT certified — only a neck
 *               that steps between radii looks like a cross-wrap weld.
 *
 * Detection is always safe (report-only). CUTTING a neck breaks the fixup's
 * additive-only contract, so the driver only cuts behind --cut-bridges and
 * only necks that are certified AND persistent across >=
 * FIXUP_BRIDGE_MIN_SUPPORT planes. */

typedef struct BridgeScanHit {
    int32_t ay, ax;        /* neck run end A (skeleton px, plane coords) */
    int32_t by, bx;        /* neck run end B */
    int32_t my, mx;        /* neck midpoint */
    int32_t len;           /* run length in skeleton px */
    float   width;         /* max chamfer half-width along the run (px) */
    float   dr;            /* |r(anchorA) - r(anchorB)| (px); -1 if unarmed */
    float   radial_dot;    /* |dot(neck dir, radial unit)|; -1 if unarmed */
    uint8_t certified;     /* radial-certified cross-wrap weld */
} BridgeScanHit;

/* Scan one plane. mask is the ORIGINAL prediction (0/1); skel its skeleton
 * (from SliceTrace_run). Hits are arena-allocated. Returns 0 on success. */
int BridgeScan_run(Arena_T arena,
                   const uint8_t *mask, const uint8_t *skel, int H, int W,
                   int radial_armed, float umb_py, float umb_px,
                   BridgeScanHit **out_hits, int32_t *out_n_hits);

/* Erase one neck from the plane: clears mask pixels within ceil(width)+1 of
 * the run's skeleton pixels. Records cleared pixels in cut_mask (may be
 * NULL). Returns pixels erased. */
int32_t BridgeScan_cut(uint8_t *mask, uint8_t *cut_mask, int H, int W,
                       const uint8_t *skel, const BridgeScanHit *hit);

int BridgeScan_selftest(void);

#endif
