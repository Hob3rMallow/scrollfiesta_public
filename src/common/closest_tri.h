#ifndef CLOSEST_TRI_INCLUDED
#define CLOSEST_TRI_INCLUDED

/*
 * closest_tri — closest point on a triangle to a query point (Ericson,
 * Real-Time Collision Detection). All double precision. Lifted verbatim from the
 * body that lived static in src/flatten/ribbon.c so the CVT/RVD remesher and the
 * flatten stage share one implementation.
 *
 * Coordinates are order-agnostic (the pipeline uses (z,y,x)); the routine treats
 * them as an opaque R^3 triple.
 */

/* Closest point on triangle (a,b,c) to p. Writes the point to out[3] and the
 * barycentric weights (u,v) such that out = a + u*(b-a) + v*(c-a). out_u/out_v
 * must be non-NULL. */
void ClosestTri_point(const double p[3], const double a[3],
                      const double b[3], const double c[3],
                      double out[3], double *out_u, double *out_v);

/* Squared distance from p to triangle (a,b,c). */
double ClosestTri_dist2(const double p[3], const double a[3],
                        const double b[3], const double c[3]);

#endif /* CLOSEST_TRI_INCLUDED */
