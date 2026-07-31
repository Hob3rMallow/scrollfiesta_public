#ifndef ATLAS_WINDING_SYNC_INCLUDED
#define ATLAS_WINDING_SYNC_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

typedef struct {
    const float *vertices;
    size_t nvertices;
    const int32_t *faces;
    size_t nfaces;
    const float *phi;
    const double *face_radius;
    const double *face_axial;
    const int32_t *face_chart;
    size_t ncharts;

    /* BoxCutter's original graph.  The first nintrinsic edges are genuine
     * shared-mesh adjacencies; the tail contains verified cube seams. */
    const int32_t *adjacency_face0;
    const int32_t *adjacency_face1;
    size_t nintrinsic_adjacency;
    size_t nadjacency;

    double spiral_a;
    double spiral_b;
    double pitch;
    int sense;
    double axial_bin_spacing;
    int phase_bins;
} AtlasWindingSyncProblem;

typedef struct {
    int32_t island0;
    int32_t island1;
    size_t observations;
    size_t cross_section_observations;
    size_t seam_observations;
    int32_t target_turn_correction;
    double mode_agreement;
    double residual_median;
    double weight;
    int eligible;
    int32_t solved_turn_correction;
    int32_t final_residual;
} AtlasWindingSyncRelation;

typedef struct {
    size_t islands;
    size_t cylindrical_bins;
    size_t strand_samples;
    size_t observations;
    size_t cross_section_observations;
    size_t seam_observations;
    size_t relations;
    size_t eligible_relations;
    size_t islands_with_relations;
    size_t relation_graph_components;
    size_t satisfied_relations;
    size_t nonzero_corrections;
    int32_t minimum_correction;
    int32_t maximum_correction;
    double relation_satisfaction_fraction;
    double prior_residual_rms;
    double final_relation_residual_rms;
    int robust_outer_iterations;
    int discrete_sweeps;
} AtlasWindingSyncStats;

/*
 * Split every BoxCutter chart into intrinsic mesh islands, then synchronize
 * the missing integer turn of each island.  Relations are discovered without
 * consulting the current absolute U: faces that coexist in cylindrical
 * (axial, phase) bins are ordered by 3-D radius, and auto-pitch predicts the
 * integer turn gap.  Verified cross-cube seams add strong equality evidence.
 *
 * The returned correction k is applied as
 *
 *     corrected_phi = phi + sense * 2*pi*k[island].
 */
int AtlasWindingSync_solve(
    Arena_T arena, const AtlasWindingSyncProblem *problem,
    int32_t **out_face_island, size_t *out_nislands,
    int32_t **out_island_turn_correction,
    double **out_island_prior, double **out_island_prior_mad,
    AtlasWindingSyncRelation **out_relation, size_t *out_nrelation,
    AtlasWindingSyncStats *stats);

#endif
