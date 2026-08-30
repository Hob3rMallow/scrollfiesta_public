#ifndef FIXUP_VIZ_INCLUDED
#define FIXUP_VIZ_INCLUDED

#include <stdint.h>
#include "../common/arena.h"
#include "slice_trace.h"
#include "slice_match.h"

/* fixup_viz — PNG viewers for the 2D prediction gap fixup. Both of us look
 * at these; every stage emits viewable bitmaps.
 *
 * plane overlay (one RGB PNG per plane):
 *   gray        original mask
 *   steel blue  skeleton
 *   bright green painted join pixels + accepted join curves
 *   violet      joins dropped by the cross-plane support filter
 *   red         safety rejects (anti-merger / cross-sheet / radial)
 *   orange      geometry rejects (tangent / score / arc / same-CC / evidence)
 *   teal        competition rejects (crossing ban / endpoint occupied)
 *   yellow dots endpoints (+ cyan tangent ticks); dark orange = excluded
 *
 * before/after: plain grayscale pair for blinking.
 * grid summary: max-projection over z — mask gray, painted green. */

int FixupViz_plane_overlay(const char *path, Arena_T arena,
                           const uint8_t *plane_after, /* [H*W] 0/1 post-paint */
                           const uint8_t *painted,     /* [H*W] painted px, may be NULL */
                           const uint8_t *skel,        /* [H*W] 0/1, may be NULL */
                           int H, int W,
                           const SliceTraceEndpoint *eps, int32_t n_eps,
                           const SliceMatchJoin *joins, int32_t n_joins,
                           const int32_t *join_support, /* [n_joins], may be NULL */
                           int32_t min_support,
                           const SliceMatchReject *rejects, int32_t n_rejects,
                           int scale);

int FixupViz_before_after(const char *path_before, const char *path_after,
                          Arena_T arena,
                          const uint8_t *plane_after, const uint8_t *painted,
                          int H, int W);

/* Magnified side-by-side crop around one endpoint pair: LEFT = before
 * (mask, skeleton, endpoints + tangent ticks), RIGHT = after (mask, painted
 * pixels, and the join curve if accepted, else a straight candidate line in
 * the given color). This is the viewer that actually shows a connection:
 * joins are ~5-12 px, invisible at full-plane scale. */
int FixupViz_pair_crop(const char *path, Arena_T arena,
                       const uint8_t *plane_after, const uint8_t *painted,
                       const uint8_t *skel, int H, int W,
                       const SliceTraceEndpoint *eps, int32_t n_eps,
                       int32_t a, int32_t b, int accepted,
                       uint8_t line_r, uint8_t line_g, uint8_t line_b,
                       int scale, int half);

int FixupViz_grid_summary(const char *path, Arena_T arena,
                          const uint8_t *vol, const uint8_t *painted,
                          int D, int H, int W);

int FixupViz_selftest(void);

#endif
