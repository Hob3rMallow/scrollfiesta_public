#ifndef INTERSECTION_CLEANUP_INCLUDED
#define INTERSECTION_CLEANUP_INCLUDED

#include <stddef.h>
#include <stdint.h>

/*
 * Bounded triangle-conflict surgery for open surface meshes.
 *
 * The detector finds:
 *   - non-adjacent, near-coplanar triangles with positive-area overlap;
 *   - non-adjacent, non-parallel triangle interpenetrations; and
 *   - edge/vertex-adjacent fold-backs whose interiors overlap.
 *
 * Conflicting triangle pairs form a graph.  Cleanup removes a greedy vertex
 * cover of that graph, so every detected conflict is opened while deleting
 * substantially fewer faces than "drop both sides".  It never fills the
 * resulting notch: doing so without a surface-ownership decision would simply
 * recreate the intersection.
 */

enum {
    INTERSECTION_HIT_OVERLAP = 1,
    INTERSECTION_HIT_STAB = 2,
    INTERSECTION_HIT_FOLD = 3
};

enum {
    INTERSECTION_HIT_MASK_OVERLAP = 1u << (INTERSECTION_HIT_OVERLAP - 1),
    INTERSECTION_HIT_MASK_STAB = 1u << (INTERSECTION_HIT_STAB - 1),
    INTERSECTION_HIT_MASK_FOLD = 1u << (INTERSECTION_HIT_FOLD - 1),
    INTERSECTION_HIT_MASK_ALL = INTERSECTION_HIT_MASK_OVERLAP |
                                INTERSECTION_HIT_MASK_STAB |
                                INTERSECTION_HIT_MASK_FOLD
};

typedef struct {
    double gap_max;             /* near-coplanar slab thickness, vox */
    double parallel_angle_deg;  /* normals within this angle are parallel */
    double max_delete_fraction; /* fail closed above this fraction of nf */
    double max_delete_fraction_hard; /* absolute cap for small-mesh floor */
    size_t min_delete_budget_faces;  /* nominal budget floor, hard-cap limited */
    size_t max_conflicts;       /* allocation/runtime safety cap */
    int include_hinges;         /* test pairs sharing exactly one vertex */
    unsigned hit_kind_mask;     /* INTERSECTION_HIT_MASK_* kinds to report/repair */
} IntersectionCleanupParams;

typedef struct {
    size_t candidate_pairs;
    size_t overlap_pairs;
    size_t stab_pairs;
    size_t fold_pairs;
    size_t conflicts;
    size_t faces_deleted;
    size_t delete_budget;
    size_t components_touched;
    size_t max_conflict_degree;
    int budget_rejected;
} IntersectionCleanupStats;

/* Optional read-only visitor for each exact conflict discovered by the audit.
 * Return zero to continue or nonzero to abort the audit (for example, when the
 * visitor cannot grow its own component-level graph). */
typedef int (*IntersectionCleanupConflictVisitor)(size_t face_a,
                                                  size_t face_b,
                                                  int hit_kind,
                                                  void *context);

void IntersectionCleanup_default_params(IntersectionCleanupParams *params);

/*
 * Read-only form of the same detector used by IntersectionCleanup_process.
 * `face_conflict_degree`, when non-NULL, receives nf counters (the caller need
 * not initialize it). No face or vertex is changed. Chart-level callers use
 * this to reject an entire proposed geometry transaction instead of deleting
 * whichever individual triangles happen to conflict.
 */
int IntersectionCleanup_audit(const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const uint8_t *face_mask,
                              const IntersectionCleanupParams *params,
                              size_t *face_conflict_degree,
                              IntersectionCleanupStats *stats);

/* Audit variant that also reports each conflicting face pair to `visitor`.
 * Detection and statistics are otherwise identical to IntersectionCleanup_audit. */
int IntersectionCleanup_audit_visit(
                              const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const uint8_t *face_mask,
                              const IntersectionCleanupParams *params,
                              size_t *face_conflict_degree,
                              IntersectionCleanupConflictVisitor visitor,
                              void *visitor_context,
                              IntersectionCleanupStats *stats);

/* Parallel form of IntersectionCleanup_audit_visit: the serial BVH traversal
 * collects candidate pairs into bounded chunks whose exact triangle tests run
 * under OpenMP; stats, visitor calls, and degree updates fold in serially in
 * collection order.  Results are identical to the serial audit and invariant
 * to thread count.  Sole difference: a max_conflicts overflow or visitor
 * abort is detected at the next chunk flush instead of at the exact offending
 * pair (rc is -1 either way; candidate_pairs may then read larger). */
int IntersectionCleanup_audit_visit_parallel(
                              const float *verts, size_t nv,
                              const int32_t *faces, size_t nf,
                              const uint8_t *face_mask,
                              const IntersectionCleanupParams *params,
                              size_t *face_conflict_degree,
                              IntersectionCleanupConflictVisitor visitor,
                              void *visitor_context,
                              IntersectionCleanupStats *stats);

/* Apply the detector's pair predicate to two triangle index triples. */
int IntersectionCleanup_pair_test(const float *verts, size_t nv,
                                  const int32_t face_a[3],
                                  const int32_t face_b[3],
                                  const IntersectionCleanupParams *params,
                                  int *shared_out);

/* Conservative approximation budgets between two separately indexed surfaces.
 * The caller initializes every face budget in [0,max_error]. This only tightens
 * budgets, never edits geometry or interprets material labels. For every nearby
 * pair, each budget is at most one quarter of a conservative projection-gap
 * lower bound on original triangle distance. Non-separated/contact pairs get
 * zero: they must retain their original surface, not be silently removed.
 * Pairs beyond 2*max_error need no change. Exhaustive two-tree broad phase;
 * no pair-prefix cap or dense concatenation. On error discard changed budgets.
 * This protects distinct input surfaces under continuous pointwise error
 * bounds; it does not repair existing intersections or adjudicate duplicates. */
int IntersectionCleanup_cross_budget(
    const float *verts_a, size_t nv_a, const int32_t *faces_a, size_t nf_a,
    const float *verts_b, size_t nv_b, const int32_t *faces_b, size_t nf_b,
    double max_error, double *budget_a, double *budget_b,
    size_t *tested_pairs);

/* Same separation envelopes across work batches of ONE indexed source mesh.
 * owner[f] is nonnegative work ownership, not a material label. Pairs within
 * an owner are skipped because its coupled fitter retains its existing guard;
 * every spatially eligible cross-owner pair is considered once. No source
 * geometry or connectivity changes, and no new pair graph is constructed. */
int IntersectionCleanup_partition_budget(const float *verts,size_t nv,
    const int32_t *faces,size_t nf,const int32_t *owner,double max_error,
    double *budget,size_t *tested_pairs);

/*
 * Remove a bounded greedy cover of all detected conflicts.
 *
 * `face_mask`, when non-NULL, has nf entries.  Only pairs for which BOTH faces
 * are marked are considered.  Pass NULL for a whole-mesh cleanup.
 *
 * Faces are compacted in place and *pnf is updated.  Vertex positions and
 * indices are unchanged.  Deleting faces can expose a bowtie vertex; callers
 * that require vertex-manifold output should subsequently run
 * PinholeFill_split_pinches.
 *
 * Returns:
 *   0  success (including an already-clean mesh)
 *  -1  invalid input, allocation failure, conflict cap, or deletion-budget
 *      rejection.  On failure the face array and *pnf are unchanged.
 */
int IntersectionCleanup_process(const float *verts, size_t nv,
                                int32_t *faces, size_t *pnf,
                                const uint8_t *face_mask,
                                const IntersectionCleanupParams *params,
                                IntersectionCleanupStats *stats);

/* Deterministic synthetic regression tests.  Returns number of failures. */
int IntersectionCleanup_selftest(void);

#endif
