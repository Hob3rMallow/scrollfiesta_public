#ifndef WINDING_REGISTER_INCLUDED
#define WINDING_REGISTER_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/* Synchronize the one-integer-turn gauge left undetermined by independently
 * lifting wrapped polar angle on disconnected mesh components.
 *
 * q is the local lifted phase in TURNS and must increase toward larger radial
 * layer order.  Radius is used only for local continuation classification and
 * same-ray order; it is never converted to an absolute turn number. */
/* Sibling evidence from radial-layer site splitting (unwrap.c).  Each entry
 * is a pair of split sites that sat as CONSECUTIVE radial layers inside the
 * same (axial, phase) lattice keys: geometrically adjacent wraps, with the
 * integer turn delta measured per key from the lifted phases.  The register
 * treats them as order observations with the mode target repeated mode_keys
 * times, and refuses a target-0 continuation between the pair -- the split
 * has already adjudicated that a fused neck joins them. */
typedef struct {
    int32_t inner;      /* split component at the smaller radius */
    int32_t outer;      /* split component one layer out */
    int32_t target;     /* mode of lround(q_inner + 1 - q_outer) over keys */
    int32_t mode_keys;  /* keys voting for the mode */
    int32_t total_keys; /* all keys the pair shared */
} WindingSiblingPair;

typedef struct {
    size_t bins, strands;
    size_t continuation_observations, order_observations;
    size_t order_observations_suppressed;
    size_t relations, eligible_relations, forest_relations;
    size_t order_relations_suppressed;
    size_t continuation_rho_suppressed;
    size_t continuation_contacts;          /* face-to-face pairs recorded as WR_OBS_CONTACT */
    size_t continuation_contact_vetoed;    /* zero continuations vetoed by contact evidence */
    size_t continuation_tangent_rejected; /* continuation pairs meeting face-to-face (along the normal) */
    size_t continuation_front_vertices; /* original open-boundary samples in the fixed-model path */
    size_t continuation_front_rejected; /* boundary pairs not facing the same gap */
    size_t continuation_order_vetoed;   /* zero-target continuations vetoed by a stacked order relation */ /* fused-neck target-0 dropped */
    size_t overlap_observations;        /* exact cross-cube same-point votes */
    size_t order_relations_contested;/* order relations kept against a
                                      * continuation closure because
                                      * their target is nonzero: a real
                                      * contradiction for the MRF to
                                      * arbitrate, not noise to delete */
    size_t continuation_components, relation_components;
    size_t packed_relation_components, packed_mesh_components;
    /* relation_conflicts is the POST-closure residual: eligible relations
     * (including forest edges sacrificed by an accepted subtree shift) whose
     * solved difference still disagrees with their target. */
    size_t relation_conflicts, observations_dropped;
    double continuation_satisfaction, order_satisfaction;
    int correction_min, correction_max;   /* pack-inclusive; NOT drift */
    /* integer loop closure (subtree-shift repair) */
    size_t repair_closers;        /* eligible non-forest cycle relations */
    size_t repair_conflicts_pre;  /* violated eligible relations pre-repair */
    size_t repair_shifts;         /* accepted subtree shifts */
    size_t repair_capped_roots;   /* trees stopped by the work cap */
    double anchor_span_pre_turns; /* anchored island lifted-q span pre-repair */
    double anchor_span_turns;     /* same after repair, pre-pack (drift) */
    /* Small-label hierarchical MRF closure over every eligible relation. */
    size_t mrf_rounds, mrf_label_changes, mrf_abstained_sites;
    double mrf_energy_before, mrf_energy_after, mrf_mean_confidence;
    size_t mrf_field_calibrated_roots, mrf_field_calibrated_sites;
    double mrf_field_calibration_r2;
    /* Iterative UV-conflict exclusion around the relation MRF.  Each accepted
     * sparse unary forbids one source component's previously occupied integer
     * correction; graph-cut pairwise terms decide which neighbouring label it
     * should take instead. */
    size_t mrf_conflict_rounds;
    size_t mrf_conflict_bins_before, mrf_conflict_bins_after;
    size_t mrf_conflict_losing_claims;
    size_t mrf_conflict_exclusions, mrf_conflict_winner_locks;
    size_t mrf_conflict_label_changes;
    int mrf_conflict_converged;
    /* Projective parent-boundary diagnostics. */
    size_t boundary_components;
    size_t boundary_relation_cuts;
    size_t boundary_lineage_cuts;
    size_t boundary_supported_relation_components;
} WindingRegisterStats;

/* Optional generalized-winding unary, already aggregated per mesh component
 * into the anchor component's integer gauge.  center is the preferred integer
 * correction; sigma and weight express uncertainty.  All three arrays have
 * ncomponents entries and must either all be present or all be NULL. */
typedef struct {
    const double *center;
    const double *sigma;
    const double *weight;
} WindingRegisterFieldUnary;

/* Optional projective boundary condition.  The arrays are indexed by the
 * register's current component/site id and have ncomponents entries.
 *
 * correction[c] is the absolute integer which must be added to the locally
 * lifted q on a component already observed by an immutable parent domain.
 * INT32_MIN means that the component has no parent evidence.  lineage[c] is
 * the parent's stable same-sheet identity; -1 means unknown.  A continuation
 * is admissible only when it does not identify two different known lineages.
 *
 * These are algorithmic boundary data from an exact input overlap, not an
 * atlas/bake anchor.  They turn extension into a constrained graph problem:
 * established components are immovable, consistent new components propagate
 * from them, and contradictory relations are cut. */
typedef struct {
    const int32_t *correction;
    const int32_t *lineage;
} WindingRegisterBoundary;

int WindingRegister_run(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const double *axial, const double *radius, const double *theta,
    const double *q,
    const int32_t *component, int32_t ncomponents,
    const int32_t *component_size, int32_t anchor_component,
    double axial_min, double pitch, int winding_sense,
    int32_t **out_correction,
    /* Optional [ncomponents] relation-graph island labels.  The anchored
     * island is 0; unobservable islands follow in deterministic atlas-pack
     * order.  These are identities, not material-U offsets. */
    int32_t **out_relation_island,
    /* Optional [ncomponents] labels for the SAME-SHEET continuation graph
     * only.  Unlike out_relation_island these labels never union components
     * merely because their radial order constrains the integer winding gauge. */
    int32_t **out_continuation_island,
    WindingRegisterStats *stats);

/* Field-aware form used by unwrap.  out_component_confidence is optional and
 * receives the MRF's local conditional sharpness per component. */
int WindingRegister_run_with_field(
    Arena_T arena,
    const float *vertices, size_t nvertices,
    const double *axial, const double *radius, const double *theta,
    const double *q,
    const int32_t *component, int32_t ncomponents,
    const int32_t *component_size, int32_t anchor_component,
    double axial_min, double pitch, int winding_sense,
    const WindingRegisterFieldUnary *field_unary,
    const WindingRegisterBoundary *boundary,
    /* Nonzero runs the sparse winner-locked conflict-exclusion rounds after
     * the initial relation MRF.  Zero returns that initial certificate. */
    int enable_conflict_exclusion,
    const WindingSiblingPair *sibling, size_t nsibling,
    const int32_t *vertex_cube,  /* [nvertices] source-cube id per vertex, or
                                  * NULL.  Arms the OVERLAP relation family:
                                  * cross-cube vertex pairs at near-identical
                                  * positions (a halo-overlap ring) are exact
                                  * turn evidence and outrank every other
                                  * family in the forest. */
    const float *vertex_normal,  /* [nvertices*3] unit normals or NULL: the
                                  * continuation TANGENCY gate (a pair that
                                  * meets along the normal is a contact, not
                                  * a continuation) */
    int32_t **out_correction,
    int32_t **out_relation_island,
    int32_t **out_continuation_island,
    float **out_component_confidence,
    WindingRegisterStats *stats);

/* Dense, bounded evidence extraction without a local gauge solve, forest,
 * packing, or continuation closure. The scroll model maps these local chart
 * ids to stable source-chart ids and solves ALL regions' factors together.
 * In particular, geometrical overlap remains soft evidence; it is not an
 * exact source-identity constraint. No vertex stride is allowed on this API.
 * kind uses the public constants below; ineligible factors are retained for
 * audit. Every observed (pair, family, target) survives; no family wins by
 * deleting competing evidence. CONTACT is diagnostic, not an equality.
 * Unlike the legacy register this API has no environment-selected gates or
 * chart-gauge-dependent rho veto. weight is not a probability. */
enum {
    WINDING_REL_CONTINUATION = 0,
    WINDING_REL_ORDER = 1,
    WINDING_REL_SEAM_ORDER = 2,
    WINDING_REL_OVERLAP = 3,
    WINDING_REL_CONTACT = 4
};
typedef struct {
    int32_t a, b, target;
    size_t observations, mode_observations;
    double agreement, residual, weight;
    int kind, eligible;
} WindingRegisterRelation;

int WindingRegister_collect(
    Arena_T arena, const float *vertices, size_t nvertices,
    const double *axial, const double *radius, const double *theta,
    const double *q, const int32_t *component, int32_t ncomponents,
    const int32_t *component_size, double axial_origin,
    double pitch, int winding_sense, const int32_t *vertex_cube,
    const float *vertex_normal,
    /* Optional [nvertices*3] original-topology outward boundary directions.
     * Zero excludes a vertex from continuation proposals, NOT from source
     * geometry or ray-order evidence. Nonzero entries must be unit vectors.
     * Proposals require both open fronts to face the intervening gap; this
     * measures a possible continuation without stitching any mesh edges. */
    const float *boundary_direction, WindingRegisterRelation **out_relations,
    size_t *out_nrelations, WindingRegisterStats *stats);

/* Synthetic order/continuation/gauge tests.  Returns failure count. */
int WindingRegister_selftest(void);

#endif
