#ifndef SCROLL_SOURCE_INCLUDED
#define SCROLL_SOURCE_INCLUDED

#include "../common/mesh_bin.h"
#include "../whole/axis_warp.h"

/* A bounded ORIGINAL source cube, before any geometric simplification.
 * Integer branches are integrated exactly on source edges, not by repeatedly
 * summing floating angles. Source vertex order fixes chart ids. A nonzero
 * cycle residual is explicit ambiguous topology, never silently repaired. */
typedef struct {
    double *axial, *radius, *theta, *q;
    float *straight, *normal;
    int32_t *component, *branch, *component_size;
    uint8_t *core; /* chart status: 0 usable, 1 core curl, 2 cyclic/ambiguous */
    size_t ncomponents, contradictory_edges;
} ScrollSource;

int ScrollSource_build(Arena_T arena, const MeshBinData *mesh,
                        const AxisWarp *axis, double axis_y, double axis_x,
                        int sense, double core_radius, ScrollSource *out);

typedef struct {
    size_t boundary_edges, boundary_vertices, invalid_edges, irregular_boundary_vertices;
} ScrollSourceBoundaryReport;

/* Original-topology boundary directions, in the supplied XYZ frame. Each
 * supported boundary vertex receives its normalized, edge-length-weighted
 * outward conormal. Interiors, isolated vertices and endpoints incident to
 * invalid edges receive zero. Coincident vertices with different source ids
 * remain distinct. No welding, face deletion, or coordinate solve occurs.
 * The caller must not treat a supported direction as a material-identity
 * certificate: folds, contacts and source errors still require inspection. */
int ScrollSource_boundary_directions(Arena_T arena, const MeshBinData *mesh,
                                      float **out_direction,
                                      ScrollSourceBoundaryReport *report);
int ScrollSource_selftest(void);

#endif
