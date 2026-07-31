#ifndef SEAM_REFINE_INCLUDED
#define SEAM_REFINE_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "../common/arena.h"
#include "seam_planes.h"

/*
 * Weld-time seam-band refinement -- make coarse (CVT ~13-vox) cube meshes
 * bridgeable WITHOUT per-cube dense rims.
 *
 * The BPA seam bridge primes each front edge by reconstructing a hinge ball
 * from the edge's incident triangle; sphere_center fails when that triangle's
 * CIRCUMRADIUS exceeds rho (<= BRIDGE_RHO_MAX = 3.0). A ~13-vox CVT triangle
 * has circumradius ~7 -> zero bridges. Plain midpoint-splitting of boundary
 * edges alone does NOT fix this: the fan children keep the far apex, so their
 * circumradius stays ~5 at every depth. The band must become WELL-SHAPED
 * ~3.5-vox triangles (circumradius ~2.0 <= the adaptive rho) -- exactly the
 * geometry the graded CVT rim proved weldable, made here at weld time instead.
 *
 * Mechanism, per round (to a fixpoint, bounded):
 *   1. split every open-boundary edge longer than target whose endpoints touch
 *      the seam band (1 incident face -> 2, winding-preserving);
 *   2. split every interior band edge longer than target (2 faces -> 4, the
 *      independent-set face-locked pattern of isotropic_remesh::split_round);
 *   3. boundary-frozen min-angle flips (WeldCleanup_flip_rounds) restore
 *      near-Delaunay shape so the next round's splits act on good triangles.
 *
 * New vertices are chord midpoints of existing same-component edges: the mesh
 * is unchanged as a point set and homeomorphic as a complex -- refinement can
 * neither fuse nor split anything. Positions deviate from the true surface by
 * at most the parent chord's sagitta. Vertices are only appended (never moved,
 * never removed); out_new_vert_src[i] names an ENDPOINT of the edge that vert
 * (nv + i) subdivided, always an index < nv + i, so callers can propagate
 * per-vertex attributes (cube provenance) in one ordered pass.
 *
 * Outputs are arena-allocated. A no-op (np == 0 or nothing to split) returns
 * the input arrays by reference with *out_n_new = 0. Returns 0 on success.
 */

typedef struct {
    float target_len;      /* refine band edges until <= this (vox)            */
    float band;            /* an edge is refined if either endpoint lies within
                            * this of a detected seam plane                    */
    float min_parent_alt;  /* skip parents thinner than this (pre-bridge sliver
                            * cull food; splitting them just culls in pieces)  */
    int   max_rounds;      /* split-round cap (independence locks throttle a
                            * face to one split per round)                     */
    int   flip_max_rounds; /* per-round + final flip iterations               */
} SeamRefineParams;

typedef struct {
    size_t planes;         /* planes considered                                */
    size_t rounds;         /* split rounds that made progress                  */
    size_t bnd_splits;     /* boundary-edge splits (1 face -> 2)               */
    size_t int_splits;     /* interior-edge splits (2 faces -> 4)              */
    size_t flips;          /* min-angle flips across all rounds                */
    size_t verts_added;
    size_t faces_added;
} SeamRefineStats;

/* Fill p with defaults sourced from pipeline_constants.h. */
void SeamRefine_default_params(SeamRefineParams *p);

int SeamRefine_process(Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const SeamPlane *planes, size_t np,
                       const SeamRefineParams *params,
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_src, size_t *out_n_new,
                       SeamRefineStats *st);

/* Provenance-complete variant for consumers that must interpolate a field
 * onto the temporary refinement.  For new vertex (nv+i), parent0[i] and
 * parent1[i] are the two endpoints whose chord midpoint created it; both are
 * less than nv+i.  The legacy process above returns parent0 only. */
int SeamRefine_process_with_parents(
                       Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const SeamPlane *planes, size_t np,
                       const SeamRefineParams *params,
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_parent0,
                       int32_t **out_new_vert_parent1,
                       size_t *out_n_new,
                       SeamRefineStats *st);

#endif
