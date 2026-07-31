#ifndef ATLAS_FIELD_APPLY_INCLUDED
#define ATLAS_FIELD_APPLY_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "atlas_candidates.h"
#include "atlas_register.h"

/* Called after every accepted robust support/chart registration round.  The
 * two arrays already contain that iteration's complete reconstructed output
 * and remain valid only for the duration of the callback. */
typedef int (*AtlasFieldApplyIterationFn)(
    void *context, const double *corrected_sample, const double *field_u,
    const AtlasRegisterIterationStats *iteration);

typedef struct {
    size_t mesh_components;
    size_t observed_components;
    size_t anchored_unobserved_components;
    size_t observations;
    size_t skipped_degenerate_triangles;
    size_t gauge_graph_components;
    double gauge_correction_rms;
    double gauge_correction_max;
    double gauge_consistency_rms;
    double reference_gauge_shift_rms;
    double reference_gauge_shift_max;
    double reference_residual_rms;
    double observation_rms;
    double observation_max;
    double edge_gradient_delta_rms;
    AtlasRegisterStats gauge_solver;
} AtlasFieldApplyStats;

/*
 * Register coarse support gauges against connected raw mesh charts, then
 * apply exactly one additive shift to every vertex in each chart.  This is a
 * rigid scalar-field extension: within-chart edge differences, triangle
 * orientations, and local metric distortion are invariant by construction.
 */
int AtlasFieldApply_solve(Arena_T arena,
                          const float *vertices,
                          size_t nvertices,
                          const int32_t *triangles,
                          size_t ntriangles,
                          const float *raw_u,
                          const float *reference_u,
                          const AtlasCandidateSampleRef *sample_ref,
                          const double *sample_target,
                          const int32_t *sample_component,
                          size_t nsamples,
                          size_t nsupport_components,
                          AtlasFieldApplyIterationFn iteration_fn,
                          void *iteration_context,
                          double *out_corrected_target,
                          double *out_u,
                          AtlasFieldApplyStats *stats);

int AtlasFieldApply_selftest(void);

#endif
