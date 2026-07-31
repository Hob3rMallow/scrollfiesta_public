#ifndef OVERLAP_QUALITY_INCLUDED
#define OVERLAP_QUALITY_INCLUDED

#include <stddef.h>

#include "../common/arena.h"
#include "../common/mesh_types.h"
#include "overlap_merge.h"

typedef struct {
    const char *reason;
    int         applied;
    size_t      input_vertices;
    size_t      input_faces;
    size_t      output_vertices;
    size_t      output_faces;
    size_t      boundary_vertices;
    size_t      input_slivers_10;
    size_t      fair_slivers_10;
    size_t      output_slivers_10;
    double      target_h;
    double      input_mean_min_angle;
    double      fair_mean_min_angle;
    double      output_mean_min_angle;
    double      input_worst_min_angle;
    double      fair_worst_min_angle;
    double      output_worst_min_angle;
} OverlapQualityStats;

/*
 * Replace the provisional repair disk by a fixed-frontier CVT mesh sampled
 * from a cotangent-harmonic height field.  Exterior faces and every interface
 * vertex are preserved exactly.  The candidate is changed only after the new
 * patch passes exact-boundary, disk-topology, manifold, and quality checks.
 *
 * Returns 0 when the quality pass was applied, -1 otherwise.
 */
int OverlapQuality_improve(Arena_T arena,
                           ComponentMesh *candidate,
                           OverlapMergeStats *merge,
                           const char *debug_dir,
                           OverlapQualityStats *stats);

#endif
