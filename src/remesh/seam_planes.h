/*
 * seam_planes.h -- shared seam-plane detection + classification queries.
 *
 * Lifted verbatim from seam_weld.c so the weld-time seam stages (bridge,
 * band refine, band recoarsen) agree on ONE definition of "where the seams
 * are" and "does this edge/face/vert sit in a seam band". Pure queries: no
 * allocation, no state. The standalone audit tools (seam_gap, seam_audit)
 * keep their own deliberate mirrors.
 *
 * A detected plane is an interior cube-boundary plane: an axis-aligned
 * multiple of cube_size with mesh substantially on BOTH sides. A cube's
 * genuine outer boundary (grid edge, or a hierarchical block's outer face)
 * has mesh on one side only and is never detected -- which is what keeps
 * every downstream seam stage away from boundaries a later weld level owns.
 */
#ifndef SEAM_PLANES_H
#define SEAM_PLANES_H

#include <stddef.h>
#include <stdint.h>

/* A detected cube-boundary plane: mesh exists substantially on both sides. */
typedef struct { int axis; double coord; } SeamPlane;

/* Detect interior cube-boundary planes (multiples of cube_size with mesh on
 * both sides beyond the export-halo sliver). `used` is a per-vertex bitmap of
 * verts referenced by any face. Returns count, fills planes[]. */
size_t SeamPlanes_detect(const float *verts, size_t nv,
                         const uint8_t *used, double cube_size,
                         double band, SeamPlane *planes, size_t max_planes);

/* A seam-FACING boundary edge LIES IN a seam plane: both endpoints within
 * `band` of the plane AND the edge nearly parallel to it. Rejects perimeter
 * edges that merely pass near the plane while running across it. */
int SeamPlanes_edge_in(const float *verts, int32_t va, int32_t vb,
                       const SeamPlane *planes, size_t np, double band);

/* The seam plane a face CROSSES: its verts straddle the plane (some below,
 * some above) and all sit within `band` of it. Returns the plane index, or
 * -1 if the face lies entirely on one side (a fold-back, not a bridge). */
int SeamPlanes_face_straddle(const float *verts, int32_t a, int32_t b, int32_t c,
                             const SeamPlane *planes, size_t np, double band);

/* Distance of vertex v to the nearest seam plane (along that plane's axis). */
double SeamPlanes_vert_dist(const float *verts, int32_t v,
                            const SeamPlane *planes, size_t np);

#endif /* SEAM_PLANES_H */
