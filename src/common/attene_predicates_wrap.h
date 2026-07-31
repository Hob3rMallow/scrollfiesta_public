#ifndef ATTENE_PREDICATES_WRAP_INCLUDED
#define ATTENE_PREDICATES_WRAP_INCLUDED

/*
 * attene_predicates_wrap — C wrapper over exact Voronoi / restricted-power-diagram
 * predicates, backed by the vendored Marco-Attene indirect-predicate machinery
 * (deps/src/attene-predicates/): an interval filter (attene_interval, fast) with a
 * BSNumber arbitrary-precision exact fallback (never rounds).
 *
 * These are the Levy restricted-Voronoi "side" predicates (cf. geogram side2/side3),
 * NOT the RST-plane orient3d predicates: a Voronoi bisector plane is defined by two
 * SITE points, not three mesh points, so it does not fit AtteneLPI3D/AtteneTPI3D. We
 * evaluate the equidistance functional directly, exactly.
 *
 * Model: all points live in a per-context bump pool addressed by an int32 id.
 *   - EXPLICIT points  = mesh vertices and Voronoi sites (double coords).
 *   - CONSTRUCTED points = clip-intersection vertices, stored implicitly as the
 *     recipe that built them (edge x bisector, or facet x two bisectors). Their
 *     exact position is never materialised; predicates run on the recipe.
 *
 * Predicates are EXACT and DETERMINISTIC (same ids + coords => same sign), which is
 * what makes adjacent mesh triangles agree on shared clip decisions (no cracks) and
 * what makes adjacent cubes agree across a seam.
 *
 * Not thread-safe per-context (the interval filter flips the FPU rounding mode and
 * points cache lambdas); use one context per worker thread.
 */

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct AttenePredCtx AttenePredCtx;
typedef int32_t AttenePt;              /* id into the ctx pool; <0 = invalid */

AttenePredCtx *AttenePred_new(void);
void           AttenePred_free(AttenePredCtx *ctx);
void           AttenePred_reset(AttenePredCtx *ctx);   /* drop all points, keep capacity */

/* Pool size, and truncate back to a previous count (drops every point with id >= n,
 * keeps capacity). Lets a caller register persistent sites once (ids 0..N-1), then
 * append per-triangle mesh vertices + constructed points and truncate(N) between
 * triangles. */
int32_t AttenePred_count(AttenePredCtx *ctx);
void    AttenePred_truncate(AttenePredCtx *ctx, int32_t n);

/* Register an explicit point (mesh vertex or Voronoi site). */
AttenePt AttenePred_explicit(AttenePredCtx *ctx, double x, double y, double z);

/* Constructed: intersection of segment (a,b) with the perpendicular bisector of
 * sites (si,sj). All four ids must be EXPLICIT. */
AttenePt AttenePred_edge_bisector(AttenePredCtx *ctx,
                                  AttenePt a, AttenePt b,
                                  AttenePt si, AttenePt sj);

/* Constructed: intersection of the plane through triangle (r,s,t) with the
 * bisectors of (si,sj) and (si,sk) — the Voronoi vertex of sites i,j,k lying in the
 * facet plane. All six ids must be EXPLICIT. */
AttenePt AttenePred_facet_bisector2(AttenePredCtx *ctx,
                                    AttenePt r, AttenePt s, AttenePt t,
                                    AttenePt si, AttenePt sj, AttenePt sk);

/* Side of point v relative to the perpendicular bisector of sites (si,sj):
 *   -1 : v strictly closer to si   (inside si's Voronoi half-space)
 *   +1 : v strictly closer to sj
 *    0 : v exactly equidistant (on the bisector)
 * v may be explicit or constructed; si,sj must be explicit. Exact + deterministic. */
int AttenePred_side_bisector(AttenePredCtx *ctx, AttenePt v, AttenePt si, AttenePt sj);

/* Approximate double coordinates of any point, for energy/geometry ONLY (never for a
 * topology decision — use the side predicate for those). Returns 0 on success, <0 if
 * the construction is degenerate (parallel / no intersection). */
int AttenePred_xyz(AttenePredCtx *ctx, AttenePt v, double *x, double *y, double *z);

/* Self-test: exercises explicit/edge/facet side predicates on configurations with
 * known answers, including near-degenerate cases the plain-double path gets wrong.
 * Returns 0 on success, <0 (and logs to stderr) on the first failure. */
int AttenePred_selftest(void);

#ifdef __cplusplus
}
#endif

#endif /* ATTENE_PREDICATES_WRAP_INCLUDED */
