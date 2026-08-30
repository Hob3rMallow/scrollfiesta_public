#ifndef JOIN_PAINT_INCLUDED
#define JOIN_PAINT_INCLUDED

#include <stdint.h>
#include "../common/arena.h"
#include "slice_trace.h"

/* join_paint — cubic-Bezier join geometry, safety checks, and painting for
 * the 2D prediction gap fixup.
 *
 * A join between endpoints a and b is the cubic Bezier with control points
 *   P0 = a,  P1 = a + (|b-a|/3) * tan_a,  P2 = b + (|b-a|/3) * tan_b,  P3 = b
 * (the Hermite form in Bernstein basis; both tangents point INTO the gap).
 *
 * Safety checks here are per-join gates the matcher calls before accepting:
 *   - merger_safe: no third mask-CC within FIXUP_MERGER_MARGIN of the path
 *     (the prototype's anti-merger test);
 *   - arc_ratio_ok: arc/chord below FIXUP_MAX_ARC_RATIO (curvature sanity);
 *   - crosses_adjacent_sheet: the path corridor in planes z +/-
 *     FIXUP_CROSS_CHECK_RANGE must not contain a component that is neither
 *     endpoint's — if a different wrap occupies the gap corridor one plane
 *     up or down, this join is almost certainly crossing it.
 *
 * Painting is ADDITIVE only: it sets mask pixels to 1 and records every
 * newly painted pixel in the parallel `painted` mask (for stats and viz).
 * Ported from 744d0ba^:scripts/step0.5-2d-trace/gap_fill_2d.c. */

enum { JOINPAINT_CROSS_SAMPLES = 20 };

/* Sample the join Bezier into out_y/out_x[JOINPAINT_CROSS_SAMPLES + 1].
 * Returns the number of points written (>= 2). */
int JoinPaint_sample(const SliceTraceEndpoint *a, const SliceTraceEndpoint *b,
                     float *out_y, float *out_x);

/* 1 if the two joins' Beziers properly cross (shared endpoints excluded). */
int JoinPaint_beziers_cross(const SliceTraceEndpoint *a1,
                            const SliceTraceEndpoint *b1,
                            const SliceTraceEndpoint *a2,
                            const SliceTraceEndpoint *b2);

/* 1 if any sample of one join's Bezier comes within min_dist px of any
 * sample of the other's. Two painted joins that merely TOUCH (without
 * crossing) would fuse four fragments into one — the matcher bans them
 * with min_dist = 2 * FIXUP_PAINT_RADIUS + 1. */
int JoinPaint_beziers_too_close(const SliceTraceEndpoint *a1,
                                const SliceTraceEndpoint *b1,
                                const SliceTraceEndpoint *a2,
                                const SliceTraceEndpoint *b2,
                                float min_dist);

/* 1 if no pixel within `margin` of the path belongs to a mask CC other than
 * a's or b's. mask_labels is the plane's SliceTrace_label_cc output. */
int JoinPaint_merger_safe(const int32_t *mask_labels, int H, int W,
                          const SliceTraceEndpoint *a,
                          const SliceTraceEndpoint *b, int margin);

/* 1 if arc_length / chord_length < max_ratio. */
int JoinPaint_arc_ratio_ok(const SliceTraceEndpoint *a,
                           const SliceTraceEndpoint *b, float max_ratio);

/* 1 if, in any plane z +/- range (excluding z itself), the middle of the
 * path corridor contains foreground from a component that is neither
 * endpoint's component in that plane. all_labels is [D*H*W] of per-plane
 * SliceTrace_label_cc results. */
int JoinPaint_crosses_adjacent_sheet(const SliceTraceEndpoint *a,
                                     const SliceTraceEndpoint *b,
                                     const int32_t *all_labels,
                                     int D, int H, int W, int z, int range);

/* Paint the join into plane (sets 1) with a disk stroke of `radius`;
 * every NEWLY painted pixel is also set in painted[]. Returns px added. */
int32_t JoinPaint_draw(uint8_t *plane, uint8_t *painted, int H, int W,
                       const SliceTraceEndpoint *a,
                       const SliceTraceEndpoint *b, int radius);

int JoinPaint_selftest(void);

#endif
