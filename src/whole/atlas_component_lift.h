#ifndef ATLAS_COMPONENT_LIFT_INCLUDED
#define ATLAS_COMPONENT_LIFT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/*
 * A deliberately narrow post-HEAD gauge solve.
 *
 * The input components already have useful local UV differentials.  Protected
 * edges say that two components must translate together.  Separation edges
 * observe a constant relative translation
 *
 *     shift[component1] - shift[component0] = target.
 *
 * No vertex-dependent deformation and no face deletion occurs here.
 */
typedef struct {
    int32_t component0;
    int32_t component1;
    double target;
    double confidence;
    int32_t source;
} AtlasComponentLiftEdge;

typedef struct {
    size_t ncomponents;
    const double *component_weight;
    /* Optional absolute gauge prior, used only to choose the one free
     * translation of each already-consistent separation graph. */
    const double *component_prior_shift;

    const AtlasComponentLiftEdge *protected_edges;
    size_t nprotected_edges;

    const AtlasComponentLiftEdge *separation_edges;
    size_t nseparation_edges;

    /* A contradictory cycle leaves that whole separation graph at HEAD. */
    double cycle_tolerance;
} AtlasComponentLiftProblem;

typedef struct {
    size_t protected_groups;
    size_t separation_graphs;
    size_t solved_graphs;
    size_t conflicting_graphs;
    size_t internal_conflicts;
    size_t moved_groups;
    size_t moved_components;
    double moved_weight;
    double max_abs_shift;
    double max_cycle_residual;
} AtlasComponentLiftStats;

int AtlasComponentLift_solve(Arena_T arena,
                             const AtlasComponentLiftProblem *problem,
                             double *out_component_shift,
                             int32_t *out_protected_group,
                             AtlasComponentLiftStats *stats);

int AtlasComponentLift_selftest(void);

#endif
