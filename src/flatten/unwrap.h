#ifndef UNWRAP_INCLUDED
#define UNWRAP_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * unwrap.h -- depth + winding UV parameterization of a scroll-surface mesh.
 *
 * Given the welded recto sheet, assign every vertex a 2D (u, v) coordinate:
 *   u = winding  -- cumulative angle around the scroll's central (umbilicus)
 *                   axis, integrated over the mesh graph so it grows past 2*pi
 *                   across successive wraps (radians).
 *   v = axial    -- position along the central axis (voxels).
 * Each connected component supplies a lift of polar phase to the universal
 * cover, determined up to one integer turn.  Those integer gauges are
 * synchronized from two measured invariants: close same-sheet continuation
 * correspondences and the radial ORDER of intersections on a common
 * (axial, polar-phase) ray.  Radius is deliberately not regressed into an
 * absolute turn number: a crushed/eccentric scroll is not an Archimedean
 * spiral.  Relation-graph islands with no observable registration are packed
 * after the anchored island and reported as such.
 *
 * This is the cheap "starting point" -- no linear solve. The result feeds the
 * tifxyz writer (src/flatten/tifxyz_write.c) and is the input an angle-
 * preserving optimizer (ABF++/SLIM) would later refine.
 * ==========================================================================*/

/* Which direction is the scroll wound around. */
typedef enum {
    UNWRAP_AXIS_AUTO = 0,  /* PCA principal axis (largest spread) */
    UNWRAP_AXIS_Z,         /* fixed volume Z (component 0 in (z,y,x)) */
    UNWRAP_AXIS_Y,         /* fixed volume Y (component 1) */
    UNWRAP_AXIS_X,         /* fixed volume X (component 2) */
    UNWRAP_AXIS_EXTERNAL   /* caller-supplied axis line (axis_dir + axis_point) */
} UnwrapAxisMode;

typedef struct {
    UnwrapAxisMode axis_mode;   /* default UNWRAP_AXIS_AUTO (0) */
    /* UNWRAP_AXIS_EXTERNAL: the real scroll axis as a 3D line, (z,y,x) order.
     * axis_dir need not be unit (normalized internally). axis_point is a point
     * on the umbilicus line; with EXTERNAL the per-slice centerline is pinned to
     * this point so radius/angle are measured from the true center, not the
     * mesh's local mean. */
    float  axis_dir[3];
    float  axis_point[3];
    /* Supplied local radial pitch (vox per turn, "wrap spacing").  It is used
     * only to classify nearby same-sheet continuations and to report the
     * diagnostic spiral fit; it never predicts an absolute component turn.
     * <=0 = estimate a local pitch from the largest component. */
    double wrap_spacing;
    /* Nonzero: also return the raw per-vertex winding field in out->phi and
     * its measured relation-graph island identity in out->island. */
    int    keep_phi;
} UnwrapOpts;

typedef struct {
    float  *uv;            /* [nv*2]: u = winding ARC LENGTH (vox), v = axial (vox).
                            * Both length-like (vox) so the grid is ~isotropic and
                            * matches vc_obj2tifxyz's metric UV convention. Each axis
                            * is shifted to start at 0. */
    float  *phi;           /* [nv] integrated winding angle (radians, NOT shifted;
                            * post component-anchoring). NULL unless opts->keep_phi. */
    int32_t *island;       /* [nv] winding relation-graph island identity.  The
                            * anchored island is 0; detached/unobservable islands
                            * are labeled in deterministic cover-pack order.  The
                            * label carries identity only, never a material-U
                            * distance. NULL unless opts->keep_phi. */
    int32_t *continuation_island; /* [nv] same-sheet continuation-only graph
                             * identity.  Radial-order observations never merge
                             * these labels. NULL unless opts->keep_phi. */
    int32_t *mesh_component; /* [nv] source mesh connected-component identity,
                             * before any winding observations are considered.
                             * NULL unless opts->keep_phi. */
    float   axis[3];       /* chosen winding axis, (z,y,x) order */
    float   centroid[3];   /* global centroid, (z,y,x) order */
    int     n_components;  /* connected components in the mesh graph */
    double  spiral_a;      /* spiral fit  r ~= a + b*(Phi / 2pi) */
    double  spiral_b;      /* radial gain per turn (0 if no spiral / flat) */
    double  spiral_r2;     /* fit quality in [.,1]; <0 if not computed */
    double  r_ref;         /* reference radius (median) used for u arc length */
    double  u_span;        /* winding arc-length extent (voxels) */
    double  v_span;        /* axial extent (voxels) */
    double  turns;         /* total winding / (2*pi) */

    /* Integer-gauge registration diagnostics.  A "relation component" is a
     * connected component of the measured continuation/order graph, not a mesh
     * component.  More than one means absolute inter-island U is unobservable;
     * the extra islands were explicitly packed for atlas display. */
    int     winding_sense; /* +1/-1: registered winding increases radially out */
    size_t  winding_bins;
    size_t  winding_strands;
    size_t  continuation_observations;
    size_t  order_observations;
    size_t  order_observations_suppressed;
    size_t  winding_relations;
    size_t  winding_eligible_relations;
    size_t  winding_forest_relations;
    size_t  winding_order_relations_suppressed;
    size_t  winding_continuation_components;
    size_t  winding_relation_components;
    size_t  winding_packed_relation_components;
    size_t  winding_packed_mesh_components;
    size_t  winding_relation_conflicts;
    size_t  winding_observations_dropped;
    double  continuation_satisfaction;
    double  order_satisfaction;
    int     turn_correction_min;
    int     turn_correction_max;
} UnwrapResult;

/* Compute the UV parameterization.
 * verts: [nv*3] (z,y,x). faces: [nf*3], 0-based. out->uv is arena-allocated.
 * Returns 0 on success, -1 on degenerate/empty input (nv<3 or nf<1). */
int Unwrap_run(Arena_T arena,
               const float *verts, size_t nv,
               const int32_t *faces, size_t nf,
               const UnwrapOpts *opts, UnwrapResult *out);

/* In-process unit tests (synthetic cylinder / spiral / degenerate inputs).
 * Returns 0 if all pass, else the number of failed checks. Logs to stderr. */
int Unwrap_selftest(void);

#endif
