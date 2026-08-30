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
    size_t relation_conflicts, observations_dropped;
    double continuation_satisfaction, order_satisfaction;
    int correction_min, correction_max;
} WindingRegisterStats;

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

/* Synthetic order/continuation/gauge tests.  Returns failure count. */
int WindingRegister_selftest(void);

#endif
