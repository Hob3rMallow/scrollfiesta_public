#ifndef ASM_TYPES_INCLUDED
#define ASM_TYPES_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#include "../common/arena.h"
#include "../common/mesh_pile.h"

/* ============================================================================
 * asm_types.h -- shared data model of the chart-assembly unwrapper.
 *
 * The UNIT is a mesh chart: one connected component of one per-cube mesh of
 * the protected pile, kept whole.  Every chart carries its own intrinsic
 * flattening (uv in the chart frame, u x v = the face normal side).  Charts
 * are related by seam correspondences across the 2-vox cube gap and by
 * layer-neighbour probes along the normal.  A rigid pose per chart places it
 * in one global 2-D frame; contradictions of that layout drive the cleaning
 * (drop a join, cut a chart at its neck, route a piece to extras).  Nothing
 * here is an axis, a pitch constant, or an integer winding number.
 *
 * Coordinates: 3-D is world (z,y,x) like every vmesh in the tree; 2-D is
 * (u,v).  Faces are chart-local int32 triples.  All sizes are size_t.
 * ==========================================================================*/

enum {
    ASM_CHART_BLOB          = 1u << 0, /* dense prediction blob: extras, never flattened */
    ASM_CHART_SUSPECT       = 1u << 1, /* flattening leaves the stress band somewhere */
    ASM_CHART_EXTRAS        = 1u << 2, /* routed to _extras by the cleaning loop */
    ASM_CHART_CUT_CHILD     = 1u << 3, /* produced by a layout-gated cut */
    ASM_CHART_HANDLE        = 1u << 4, /* genus > 0 (opened or routed to extras) */
    ASM_CHART_NONMANIFOLD   = 1u << 5, /* manifold audit failed: extras */
    ASM_CHART_FLAT_FAILED   = 1u << 6, /* flattening failed (no disk / flips): extras */
    ASM_CHART_TINY          = 1u << 7, /* below the minimum area: extras */
    ASM_CHART_MIRROR        = 1u << 8  /* parity solve says: mirror u before placing */
};

/* Charts that never enter the layout. */
#define ASM_CHART_EXCLUDED_MASK \
    (ASM_CHART_BLOB | ASM_CHART_EXTRAS | ASM_CHART_NONMANIFOLD | \
     ASM_CHART_FLAT_FAILED | ASM_CHART_TINY)

/* Runtime evidence, deliberately absent from the clean-chart cache. Sharing
 * a coordinate frame does not establish material continuity. */
enum { ASM_PLACE_NONE, ASM_PLACE_ROOT, ASM_PLACE_SEAM, ASM_PLACE_LAYER };

typedef struct AsmChart {
    int32_t  id;          /* 0-based over the run */
    int32_t  cube;        /* index into the pile */
    int32_t  comp;        /* component ordinal within the cube mesh */
    int32_t  parent;      /* -1, or the chart id this was cut from */
    uint32_t flags;       /* ASM_CHART_* */
    int32_t  component;   /* relation-graph component after the pose solve, -1 unset */
    int32_t  placed;      /* 1 once the chart has a pose in the global frame */
    int32_t  placement_group, placement_state, placement_parent;
    size_t   placement_support;
    double   placement_rms;

    size_t   nv, nf;
    int32_t *vid;         /* [nv] local vertex index in the cube mesh */
    float   *xyz;         /* [nv*3] world (z,y,x) */
    float   *nrm;         /* [nv*3] unit vertex normals from the face winding */
    float   *uv;          /* [nv*2] intrinsic flattening, chart frame; NULL if none */
    double  *placed_uv;   /* [nv*2] repaired GLOBAL coordinates; runtime only.
                          * Overrides pose(uv) after the native repair stage. */
    int32_t *faces;       /* [nf*3] chart-local indices */
    uint8_t *boundary;    /* [nv] 1 = touches a boundary edge */
    uint8_t *repaired;    /* [nv] 1 = near fixer-added voxels, or NULL when unknown */

    double   area3d;      /* sum of 3-D triangle areas, vox^2 */
    double   area_uv;     /* sum of |2-D triangle areas| */
    double   sigma_lo;    /* smallest singular value of the map over faces */
    double   sigma_hi;    /* largest singular value */
    double   stress_frac; /* fraction of faces outside [1-band, 1+band] */
    size_t   n_flipped;   /* faces with a negative 2-D orientation */
    double   centroid[3]; /* 3-D vertex mean */
    float    bbox_lo[3], bbox_hi[3];

    /* Pose in the global 2-D frame: g = R(theta) * m(uv) + (x, y), where m
     * mirrors u when ASM_CHART_MIRROR is set. */
    double   pose_x, pose_y, pose_theta;
} AsmChart;

enum {
    ASM_REL_DROPPED   = 1u << 0, /* removed by the cleaning loop */
    ASM_REL_REPAIRED  = 1u << 1, /* runs through fixer-added voxels */
    ASM_REL_SUSPECT   = 1u << 2, /* touches a SUSPECT chart */
    ASM_REL_CONTACT   = 1u << 3, /* rejected by the tangency gate (kept for the ledger) */
    ASM_REL_PARITY    = 1u << 4, /* chart b is mirrored relative to chart a */
    ASM_REL_SWITCHED  = 1u << 5, /* switched off by the robust pose solve */
    ASM_REL_WEAK      = 1u << 6, /* gate-rejected (isometry / unimodality) but rigidly fitting: PLACEMENT evidence only,
                                  * never a pose-graph edge (the floating single-cube charts of the 10x10x10 sit beside
                                  * their holes because exactly these seams were rejected, 2026-09-09) */
    ASM_REL_LAYOUT    = 1u << 7, /* confirmed by the LAYOUT after placement: the two charts' adjacent audit cells agree
                                  * in 3-D (same layer, aligned normals, steps under the wrap gate) over a seam's length
                                  * with no contradiction; the transform is their placed relative pose.  The stitcher's
                                  * decision on geometric continuity, ledgered (2026-09-09). */
    ASM_REL_SHORT     = 1u << 8, /* rejected by the seam-length gate but rigidly fitting (rms <= 4 over >= 6 correspondences):
                                  * PLACEMENT evidence only, like WEAK (2026-09-10: 975 of the 10x10x10's 1,403 length
                                  * rejections were 15-30 vox long and two thirds of them true continuations) */
    ASM_REL_READMIT   = 1u << 9, /* readmitted as a join by the post-placement seam re-solve: the placed layout put the
                                  * pair within ASM_READMIT_VOX of the pose the seam implies (its DROPPED / SWITCHED /
                                  * WEAK / SHORT bits are cleared; this bit is the ledger's record, 2026-09-10) */
    ASM_REL_CROSSWRAP = 1u << 10, /* the UNIMODALITY gate rejected this seam: its pairwise rotations do not agree with
                                  * one rigid motion, which is what a sheet switch inside a seam looks like.  That
                                  * statistic is the cross-wrap detector (2026-09-17: disarming it for low-rms fused
                                  * seams stacked the 21x5x5's pose graph).  The bit records the verdict so the
                                  * rejection is visible after relate: WEAK alone conflates it with the isometry
                                  * rejection and with the bridge premium's demotions, and hands all three to
                                  * placement as evidence.  ASM_WRAP_VETO_ANCHOR / _JOIN / _DISCOVER decide whether it is also a veto. */
    ASM_REL_TORN      = 1u << 11  /* a recorded TEAR (user-approved policy, 2026-09-27): the layout keeps the two sides
                                  * apart and the relation is no obligation.  Tears are read with the input run; the
                                  * audit ledgers every tear (stage7_tears.csv) and measures the source regions both
                                  * with and without them. */
};

/* relations that are placement evidence only: exact anchors, never a pose-graph edge or a join */
#define ASM_REL_PLACEMENT_ONLY (ASM_REL_WEAK | ASM_REL_SHORT)



/* A seam relation: chart b's frame expressed in chart a's frame,
 *     p_a = R(theta) * p_b + (tx, ty),
 * measured on ordered boundary correspondences across one cube seam. */
typedef struct AsmRelation {
    int32_t  a, b;
    uint32_t flags;         /* ASM_REL_* */
    double   theta, tx, ty;
    double   w_xy;          /* information weight of the translation, 1/vox^2 */
    double   w_theta;       /* information weight of the rotation, 1/rad^2 */
    double   rms;           /* residual of the rigid fit, vox */
    double   seam_len;      /* matched seam length, vox */
    double   normal_gap;    /* median |displacement . n| of the pairs, vox */
    double   mean_gap;      /* median 3-D pair distance, vox */
    size_t   n_corr;
    int32_t  corr_first;    /* into AsmRun.corr */
    int32_t  corr_count;
    double   robust_w;      /* IRLS / switch weight after the pose solve (1 = trusted) */
    double   residual;      /* normalized residual after the pose solve */
    double   flag_mass;     /* accumulated contradiction witness mass */
    uint32_t continuity;    /* runtime source/geometry certificate, independent of flags */
    double   continuity_rms;
} AsmRelation;

/* Original source-face coordinates of a physical tangent displacement.
 * The inverse maps an orthonormal vertex-frame gap onto the two source
 * edges. UV deformation carries those edges; it cannot freeze the gap's
 * direction in an obsolete local coordinate system. */
typedef struct AsmTrim {
    int32_t face;
    uint32_t valid;
    double inverse[4];
} AsmTrim;

/* One boundary correspondence: chart-local vertex indices. */
typedef struct AsmCorr {
    int32_t va, vb;
    /* Measured B-A displacement in each intrinsic tangent map. Do not weld
     * trimmed source boundaries to the same UV point. */
    float gap_a[2], gap_b[2];
    int32_t run;            /* ordered boundary run, -1 = unresolved */
    uint32_t valid;         /* bit 0/1: nondegenerate source tangent on A/B */
    AsmTrim trim_a, trim_b;
} AsmCorr;

/* Layer neighbour: chart b lies `d_median` vox along chart a's normal, on
 * `side` (+1 along the normal, -1 against it).  Negative evidence and the
 * local layer distance ("pitch") ruler. */
typedef struct AsmLayerPair {
    int32_t a, b;
    int32_t count;
    int8_t  side;
    double  d_median;
    int32_t hit_first, hit_count;   /* representative hits in AsmRun.layer_hits */
    int8_t  k;                      /* median observed crossing order; holes can hide intermediate wraps */
    int8_t  nsign;                  /* majority sign of n_a . n_b over the hits (+1 same orientation, -1 opposed, 0 unknown:
                                     * snapshots written before 2026-09-17); a mixed-graph orientation bit for the wrap gate */
} AsmLayerPair;

/* One layer probe hit: from vertex va of chart a, along side * n, the ray
 * met chart b's face fb at barycentric (l0, l1, 1-l0-l1) after dist vox. */
typedef struct AsmLayerHit {
    int32_t a, va, b, fb;
    float   l0, l1;
    float   dist;
    int8_t  side;
    int8_t  order;   /* 1 = first surface crossed along the probe, 2 = second, ... */
    int8_t  nsign;   /* sign of n_a . n_face_b at the hit (+1 / -1; 0 = unknown, older snapshots) */
} AsmLayerHit;

/* Growable arena tables (doubling; old blocks are left to the arena). */
typedef struct AsmRun {
    Arena_T        arena;
    int            continuity_ready; /* source-aware placement owns pose changes */

    MeshPileEntry *pile;
    size_t         n_cubes;

    AsmChart      *charts;
    size_t         n_charts, cap_charts;

    AsmRelation   *rels;
    size_t         n_rels, cap_rels;

    AsmCorr       *corr;
    size_t         n_corr, cap_corr;

    AsmLayerPair  *layers;
    size_t         n_layers, cap_layers;

    AsmLayerHit   *layer_hits;
    size_t         n_layer_hits, cap_layer_hits;

    /* per-chart local layer distance (vox) from the layer probes; 0 = unknown */
    double        *chart_layer_d;
    double         layer_d_global;   /* median over all measured charts */
    double         cube_size;        /* pile cube edge, vox */

    /* Versioned private repair proposals, owned by arena. These bytes never
     * supply live chart coordinates; repair validates them before reuse. */
    uint8_t       *repair_resume;
    size_t         repair_resume_size;
} AsmRun;

/* Append helpers (arena doubling).  Return the index of the new element. */
size_t AsmRun_push_chart(AsmRun *run, const AsmChart *c);
size_t AsmRun_push_rel(AsmRun *run, const AsmRelation *r);
size_t AsmRun_push_corr(AsmRun *run, const AsmCorr *c);
size_t AsmRun_push_layer(AsmRun *run, const AsmLayerPair *l);
size_t AsmRun_push_layer_hit(AsmRun *run, const AsmLayerHit *h);

static inline int AsmChart_in_layout(const AsmChart *c)
{
    return (c->flags & ASM_CHART_EXCLUDED_MASK) == 0u && c->uv != NULL;
}

/* placed also marks floating pose-graph components. Only placement_state
 * certifies that a chart belongs to the shared global sheet frame. */
static inline int AsmChart_registered(const AsmChart *c)
{
    return AsmChart_in_layout(c) && c->placed && c->placement_state != ASM_PLACE_NONE;
}

/* One authoritative post-placement coordinate accessor. Never quantize a
 * repaired global field back into the float32 clean-chart cache. */
static inline void AsmChart_point(const AsmChart *c, size_t v, double p[2])
{
    if (c->placed_uv) { p[0] = c->placed_uv[2*v]; p[1] = c->placed_uv[2*v+1]; return; }
    double u = c->uv[2*v], w = c->uv[2*v+1];
    if (c->flags & ASM_CHART_MIRROR) u = -u;
    double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
    p[0] = ct*u-sn*w+c->pose_x; p[1] = sn*u+ct*w+c->pose_y;
}

/* A JOIN: a relation the gates accepted (or the layout confirmed), the cleaning
 * kept, the pose solve did not switch off, and that is not placement evidence
 * only.  Lineages (material identity) are the components of the join graph. */
static inline int AsmRel_is_join(const AsmRelation *r)
{
    return (r->flags & (ASM_REL_DROPPED | ASM_REL_CONTACT | ASM_REL_SWITCHED | ASM_REL_PLACEMENT_ONLY)) == 0u;
}

#endif
