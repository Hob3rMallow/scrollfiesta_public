/*
 * clipper2_wrap.h — C-callable wrapper for Clipper2 polygon union.
 *
 * Resolves self-intersecting 2D polygons into simple sub-polygons
 * via EvenOdd fill rule union.  Scales double coords to int64 internally.
 *
 * All output arrays are malloc'd — caller must free() them.
 */
#ifndef CLIPPER2_WRAP_INCLUDED
#define CLIPPER2_WRAP_INCLUDED

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Clipper2_union — resolve a self-intersecting polygon into simple polygons.
 *
 * Input:
 *   pts_xy    [n * 2]  2D polygon vertices (double)
 *   n                  number of vertices
 *
 * Output (all malloc'd — caller frees):
 *   out_pts      [total_verts * 2]  vertices of all sub-polygons
 *   out_counts   [n_polys]          vertex count per sub-polygon
 *   out_n_polys                     number of resulting simple polygons
 *
 * Returns 0 on success, -1 on failure.
 */
int Clipper2_union(const double *pts_xy, size_t n,
                   double **out_pts, int **out_counts,
                   int *out_n_polys);

#ifdef __cplusplus
}
#endif

#endif /* CLIPPER2_WRAP_INCLUDED */
