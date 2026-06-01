#ifndef ORIENT_WELD_INCLUDED
#define ORIENT_WELD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"

/*
 * Post-weld cross-component orientation.
 *
 * After the seam weld, the mesh is a set of vertex-connected components (scroll
 * wraps, fragments, bridge pieces). Each is internally winding-consistent (per-
 * cube OrientMesh + the bridge's directed glue), but a component's GLOBAL sign
 * can still be flipped relative to its neighbours -- a "backward-wound block".
 *
 * The existing OrientMesh sign-anchor cannot fix this post-weld: it anchors each
 * component against a normal field, but post-weld the only normals available are
 * recomputed FROM the winding, so a consistently-backward component agrees with
 * its own (also-backward) normals -- no disagreement to detect.
 *
 * This pass anchors each component SPATIALLY instead: two components that sit
 * within `radius` of each other should have normals pointing the same way where
 * they meet. Process components largest-first (the biggest is the trusted
 * anchor); for each next component, sum dot(n_v, n_w) over near vertex pairs to
 * an already-oriented component, and if the net is negative FLIP the whole
 * component (reverse every face winding). Components with no oriented neighbour
 * within `radius` are left unchanged.
 *
 *   verts[nv*3]  float (z,y,x), read only
 *   faces[nf*3]  int32 0-based, winding flipped IN PLACE for backward components
 *   radius       neighbour ball (voxels) for the cross-component vote
 *   out_flipped  (optional) number of components whose winding was reversed
 *
 * Uses `arena` for transient scratch (restored on return). Returns 0 on success.
 */
int OrientWeld_components(Arena_T arena,
                          const float *verts, size_t nv,
                          int32_t *faces, size_t nf,
                          float radius,
                          size_t *out_flipped);

#endif
