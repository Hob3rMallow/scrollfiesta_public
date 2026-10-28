#ifndef METRIC_BOUNDS_INCLUDED
#define METRIC_BOUNDS_INCLUDED

/* C interface to a filtered exact predicate; this computes no new geometry.
 * XYZ and UV contain one triangle in vertex order, with 3 and 2 coordinates
 * per vertex. Bounds are {lower0, upper0, lower1, upper1}. A violation sets
 * the corresponding bit. Equality passes. Input doubles are interpreted
 * exactly, including the bounds before squaring them.
 *
 * Returns -1 for invalid input or a degenerate physical triangle, without
 * modifying result. Otherwise returns 0. The vendored bigint allocator aborts
 * on allocation exhaustion. The filter uses
 * conservative forward errors; uncertain signs use vendored BSNumber.
 */
typedef struct MetricBoundsResult {
    unsigned violated;
    int exact_fallback;
} MetricBoundsResult;

#ifdef __cplusplus
extern "C" {
#endif
int MetricBounds_classify(const float xyz[9], const double uv[6],
    const double bounds[4], MetricBoundsResult *result);
int MetricBounds_selftest(void);
#ifdef __cplusplus
}
#endif
#endif
