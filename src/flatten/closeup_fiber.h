#ifndef CLOSEUP_FIBER_INCLUDED
#define CLOSEUP_FIBER_INCLUDED

#include "slim_refine.h"
#include "../common/zarr_u8.h"

typedef struct CloseupFiberStats {
    size_t patches, measured_patches, faces, samples;
    size_t chunk_misses, chunk_missing;
    double area, measured_area, seconds;
} CloseupFiberStats;

/* Measure the CT cross on overlapping physical tangent patches of each actual
 * material component. frame_uv only indexes neighborhoods; its axes are never
 * treated as measured fibers. out[nf] receives physical arms and coherence,
 * or zero weight where the CT supplies no usable evidence. XYZ stays exact. */
int CloseupFiber_measure(const float *xyz, const float *normal,
                          const double *frame_uv, size_t nv,
                          const int32_t *faces, size_t nf, const char *raw,
                          SlimRefineDirection *out, CloseupFiberStats *stats);
/* The same measurement with a caller-owned, already-open CT volume. The caller
 * retains its cache and handle on success and failure. Use a separate handle
 * per reading worker; calls on one handle must be serial. Scratch is local to
 * this call, and chunk counters report this call's deltas. Patch raster memory
 * depends on physical bounds; frame_uv bins are not a memory limit. */
int CloseupFiber_measure_volume(ZarrU8 *volume, const float *xyz, const float *normal,
                                 const double *frame_uv, size_t nv,
                                 const int32_t *faces, size_t nf,
                                 SlimRefineDirection *out, CloseupFiberStats *stats);
int CloseupFiber_selftest(void);

#endif
