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

/*
 * OrientWeld_components_axis -- the same pass with a RADIAL fallback for
 * components the spatial vote cannot reach.
 *
 * A detached component (nothing from another component within `radius`
 * anywhere -- e.g. an outer-shell fragment whose connecting geometry was
 * never meshed) collects ZERO near pairs, so the spatial vote is silent and
 * the component keeps whatever global sign the per-cube (1,1,1) MLS anchor
 * gave it. For a scroll that sign is azimuth-dependent: shells on the far
 * side of the umbilicus come out backward (2026-07-10 orientation audit:
 * 13/26 components on grid_4x5x5_v5_full).
 *
 * Fallback: when a component's spatial evidence is WEAK -- fewer than 100
 * near pairs and less than one pair per component vertex (a glancing touch;
 * on the real scroll a 39k-vert shell collected 14 incidental pairs voting
 * +12 against a radial vote of -37k) -- and an axis is given (axis_point /
 * axis_dir, (z,y,x) voxel coords -- the scroll umbilicus), flip it if
 * sign( sum_v n_v . r_hat_v ) disagrees with the ANCHOR (largest)
 * component's radial sign. Radial facing is per-wrap-consistent and azimuth-
 * free, so it anchors far components without any proximity requirement.
 * Components with decisive spatial evidence are decided by the spatial vote
 * exactly as before. axis_point == NULL (or a ~zero direction) reproduces
 * OrientWeld_components bit-for-bit.
 *
 *   out_radial  (optional) components decided by the radial fallback
 *               (counted whether or not they needed a flip)
 */
int OrientWeld_components_axis(Arena_T arena,
                               const float *verts, size_t nv,
                               int32_t *faces, size_t nf,
                               float radius,
                               const float *axis_point,
                               const float *axis_dir,
                               size_t *out_flipped,
                               size_t *out_radial);

#endif
