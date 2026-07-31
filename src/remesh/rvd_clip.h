#ifndef RVD_CLIP_INCLUDED
#define RVD_CLIP_INCLUDED

/*
 * rvd_clip — Restricted Voronoi Diagram of a point set (sites) against a triangle
 * mesh surface. Each mesh triangle is clipped, exactly, against the Voronoi bisector
 * half-planes of the nearby sites (Yan 2009 / Levy), using the exact "side"
 * predicates in attene_predicates_wrap. The restricted cell of each site is the
 * union of the clipped polygons; we accumulate its area + centroid (for Lloyd CVT).
 *
 * Every clip-polygon vertex is exactly one of three exactly-constructible kinds —
 * mesh vertex, mesh-edge x bisector, or facet x two bisectors — so the clip is
 * crack-free and deterministic (adjacent triangles agree on shared decisions).
 */

#include "../common/arena.h"
#include "sizing_field.h"
#include <stddef.h>
#include <stdint.h>

typedef struct Rvd_T *Rvd_T;

/* Build over a triangle mesh. verts: double[nv*3] (any coord order); faces:
 * int32[nf*3], 0-based. Both are copied into the arena. */
Rvd_T Rvd_new(Arena_T arena, const double *verts, size_t nv,
              const int32_t *faces, size_t nf);

/* Compute the RVD of `sites` (double[nsites*3]) and accumulate per site the
 * restricted-cell AREA (out_area[nsites]) and area-weighted CENTROID
 * (out_centroid[nsites*3], projected to the site's own cell). Caller allocates the
 * outputs; they are overwritten. A site with an empty restricted cell gets area 0
 * and centroid = its own position.
 *
 * If out_NA/out_bNA are non-NULL (pass BOTH or neither), ALSO accumulates the CWF
 * normal-anisotropy moments per cell: out_NA[6*nsites] = symmetric N_i = integral of
 * n*n^T over the cell (entries xx,xy,xz,yy,yz,zz), and out_bNA[3*nsites] = b_i =
 * integral of (n*n^T)x, where n is the surface normal. These drive the fixed-cell
 * CWF solve (lambda_NA*N_i + lambda_CVT*m_i*I) x_i = lambda_NA*b_i + lambda_CVT*m_i*c_i.
 *
 * If out_dual_tris is non-NULL, ALSO extracts the Restricted-Delaunay dual as
 * oriented site-index triples (arena-allocated int32[*out_ndual*3]). Returns 0.
 *
 * If `dens` is non-NULL, every accumulated sub-triangle area is scaled by the
 * CWF energy weight Rvd_density_rho(dens, sub-triangle centroid) — this grades
 * the converged site density (denser where h(x) is small). `dens == NULL` is the
 * uniform path and is BYTE-IDENTICAL to the unweighted accumulation (the guard
 * applies no arithmetic when off). Density never enters a geometric predicate, so
 * the exact clip topology + dual are unchanged for a fixed site set. */
int Rvd_accumulate(Rvd_T rvd, const double *sites, size_t nsites,
                   double *out_area, double *out_centroid,
                   double *out_NA, double *out_bNA,
                   int32_t **out_dual_tris, size_t *out_ndual,
                   const CvtField *dens);

/* Total surface area of the input mesh. The sum of all restricted-cell areas from
 * Rvd_accumulate must equal this (partition property) — the core correctness check. */
double Rvd_mesh_area(Rvd_T rvd);

#endif /* RVD_CLIP_INCLUDED */
