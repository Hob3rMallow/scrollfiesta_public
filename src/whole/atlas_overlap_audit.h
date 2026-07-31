#ifndef ATLAS_OVERLAP_AUDIT_INCLUDED
#define ATLAS_OVERLAP_AUDIT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

typedef struct {
    int32_t component0;
    int32_t component1;
    size_t face_pairs;
    size_t unique_faces0;
    size_t unique_faces1;
    size_t component_faces0;
    size_t component_faces1;
    int32_t turn_mode;
    double turn_agreement;
    double phase_residual_median;
    double phase_residual_mad;
    double parameter_du_median;
    double parameter_du_mad;
    double registered_du_median;
    double registered_du_mad;
    /* Constant component-shift difference which would reproduce the
     * registered separation on this overlap support. */
    double registered_parameter_shift_median;
    double registered_parameter_shift_mad;
} AtlasOverlapBundle;

typedef struct {
    int32_t component0;
    int32_t component1;
    int32_t face0;
    int32_t face1;
    int32_t bundle;
    int32_t turn;
    double phase_residual;
    double parameter_du;
    double registered_du;
    double registered_parameter_shift;
} AtlasOverlapPair;

typedef struct {
    AtlasOverlapBundle *bundles;
    size_t nbundles;
    AtlasOverlapPair *pairs;
    size_t npairs;
    size_t exact_face_pairs;
    size_t same_component_pairs;
    size_t cross_component_pairs;
    /* Broad-phase completeness diagnostics.  A successful audit is only
     * authoritative when broad_phase_complete is true and indexed_faces
     * equals the input face count. */
    size_t indexed_faces;
    size_t broad_phase_records;
    size_t broad_phase_cells;
    size_t broad_phase_candidate_pairs;
    double broad_phase_cell_size;
    int broad_phase_complete;
    /* Retained for callers/log formats which predate the complete sparse
     * broad phase.  Successful builds now leave this false. */
    int pair_buffer_truncated;
} AtlasOverlapAudit;

/* Exact 2-D triangle intersections in (u,v), aggregated by retained mesh
 * component pair.  phi and registered_u are observations only; this function
 * never changes geometry or parameterization. */
int AtlasOverlapAudit_build(Arena_T arena,
                            const int32_t *faces, size_t nfaces,
                            size_t nvertices,
                            const double *u, const double *v,
                            const float *registered_u,
                            const float *phi,
                            const int32_t *vertex_component,
                            size_t ncomponents,
                            AtlasOverlapAudit *out);

int AtlasOverlapAudit_selftest(void);

#endif
