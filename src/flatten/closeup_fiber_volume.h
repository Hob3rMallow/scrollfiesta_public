#ifndef CLOSEUP_FIBER_VOLUME_INCLUDED
#define CLOSEUP_FIBER_VOLUME_INCLUDED

#include "closeup_fiber.h"

/* Read an existing inferred 3-D fiber volume at each physical face centroid
 * and at the standard depth offsets along its mean vertex normal. Mesh arrays
 * and returned arms use (z,y,x); the volume's decoded direction uses (x,y,z).
 * No UV coordinates, smoothing, geometry changes, or inferred-frame priors.
 * Unmeasured/low-coherence faces retain exactly zero arms and weight. Each
 * face is one sampling neighborhood in stats.patches; samples counts usable
 * depth samples. Cache counters are unavailable through FiberVolume's API.
 * The caller must freeze and verify the field's inference_identity.json.
 */
int CloseupFiberVolume_measure(const float *xyz, const float *normal, size_t nv,
                              const int32_t *faces, size_t nf, const char *field,
                              SlimRefineDirection *out, CloseupFiberStats *stats);
int CloseupFiberVolume_selftest(void);

#endif
