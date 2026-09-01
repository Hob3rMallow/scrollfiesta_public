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
typedef struct {
    size_t bins, strands;
    size_t continuation_observations, order_observations;
    size_t order_observations_suppressed;
    size_t relations, eligible_relations, forest_relations;
    size_t order_relations_suppressed;
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
    /* Nonzero runs the sparse winner-locked conflict-exclusion rounds after
     * the initial relation MRF.  Zero returns that initial certificate. */
    int enable_conflict_exclusion,
    int32_t **out_correction,
    int32_t **out_relation_island,
    int32_t **out_continuation_island,
    float **out_component_confidence,
    WindingRegisterStats *stats);

/* Synthetic order/continuation/gauge tests.  Returns failure count. */
int WindingRegister_selftest(void);

#endif
