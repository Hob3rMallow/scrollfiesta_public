#ifndef SLICE_TRACE_INCLUDED
#define SLICE_TRACE_INCLUDED

#include <stdint.h>
#include "../common/arena.h"

/* slice_trace — per-plane skeleton tracing for the 2D prediction gap fixup.
 *
 * Skeletonizes one binary (0/1) z-plane of the prediction volume (Zhang-Suen
 * thinning + short-spur pruning + re-thin), labels the skeleton's 8-connected
 * components, and extracts curve ENDPOINTS with outward unit tangents (the
 * tangent points INTO the gap, away from the curve body). Endpoints are
 * snapped outward to the mask edge (Zhang-Suen retracts curve tips by a few
 * pixels), so gap distances measure mask-edge to mask-edge.
 *
 * Endpoints on skeleton components smaller than FIXUP_MIN_CURVE_PX are
 * dropped (spur noise). Endpoints within FIXUP_BORDER_EXCLUDE px of the plane
 * border or of a data-absent (missing-cube) region are kept but flagged
 * `excluded`: those curves are clipped by the volume, not broken papyrus,
 * and must never be matched.
 *
 * Ported from the recovered 2026-05 prototype (744d0ba^:scripts/
 * step0.5-2d-trace/gap_fill_2d.c) with arena allocation; the prototype's
 * mask-erasing small-CC removal is deliberately NOT ported — the fixup is
 * additive-only, it never deletes prediction foreground.
 *
 * Deterministic: results depend only on the plane contents (scan order). */

typedef struct SliceTraceEndpoint {
    int32_t y, x;           /* plane coords, after snap-to-mask-edge */
    int32_t skel_cc;        /* 1-based skeleton CC label */
    int32_t mask_cc;        /* 1-based mask CC label at (y,x) */
    float   tan_dy, tan_dx; /* unit outward tangent (points into the gap) */
    float   curv;           /* turning angle (radians) along curvature walk */
    uint8_t excluded;       /* 1 = border/absent-region endpoint: unmatchable */
} SliceTraceEndpoint;

/* Label 8-connected components of a binary H*W mask, 1-based (0 = bg).
 * labels is caller-allocated [H*W]; scratch comes from `arena` (mark-restored
 * internally). Returns the component count (>= 0), or -1 on bad input. */
int32_t SliceTrace_label_cc(Arena_T arena, const uint8_t *mask, int H, int W,
                            int32_t *labels);

/* Trace one plane. mask is [H*W] 0/1. absent is [H*W] (1 = data-absent
 * region) or NULL. mask_labels is the plane's SliceTrace_label_cc output.
 * Outputs are arena-allocated: *out_skel ([H*W] 0/1 skeleton, may be NULL if
 * not wanted), *out_eps / *out_n_eps the endpoint list. Empty planes yield
 * zero endpoints and succeed. Returns 0 on success. */
int SliceTrace_run(Arena_T arena,
                   const uint8_t *mask, int H, int W,
                   const uint8_t *absent,
                   const int32_t *mask_labels,
                   uint8_t **out_skel,
                   SliceTraceEndpoint **out_eps, int32_t *out_n_eps);

int SliceTrace_selftest(void);

#endif
