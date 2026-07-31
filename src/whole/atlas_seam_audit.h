#ifndef ATLAS_SEAM_AUDIT_INCLUDED
#define ATLAS_SEAM_AUDIT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../unroll/piece_set.h"

typedef struct {
    int32_t component0;
    int32_t component1;
    size_t pairs;
    double target_shift_median;
    double target_shift_mad;
    double distance_median;
    double phase_residual_median;
    double phase_residual_mad;
    double registered_du_median;
    double registered_du_mad;
} AtlasSeamBundle;

typedef struct {
    int32_t vertex0;
    int32_t vertex1;
    int32_t component0;
    int32_t component1;
    int32_t bundle;
    double target_shift;
    double distance;
    double phase_residual;
    double registered_du;
} AtlasSeamPair;

typedef struct {
    AtlasSeamBundle *bundles;
    size_t nbundles;
    AtlasSeamPair *pairs;
    size_t npairs;
    size_t boundary_vertices;
    size_t mutual_pairs;
    size_t same_component_pairs;
    size_t cross_component_pairs;
} AtlasSeamAudit;

/* Read-only counterpart of the later UV weld.  It detects mutual-nearest
 * cross-cube copies of the same physical boundary surface, then aggregates
 * the constant shift needed to make their current U coordinates coincide.
 * The retained `pairs` are the exact pointwise rows used by the smooth FEM
 * quotient solve. Geometry, faces, and UV are never mutated. */
int AtlasSeamAudit_build(Arena_T arena, const PieceSet *ps,
                         const double *parameter_u,
                         const int32_t *vertex_component,
                         size_t ncomponents,
                         AtlasSeamAudit *out);

int AtlasSeamAudit_selftest(void);

#endif
