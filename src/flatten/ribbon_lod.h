#ifndef RIBBON_LOD_INCLUDED
#define RIBBON_LOD_INCLUDED

#include "../common/mesh_bin.h"

/* Adaptive dyadic patches of an EXISTING regular quad ribbon. Complete
 * topological disks may be replaced; holes, cuts, ambiguous UV occupancy and
 * every physical boundary edge are retained. Output vertices are exact input
 * vertices, so all UV/provenance channels transfer via source_vertex.
 *
 * Geometry error is checked at every vertex of the overlay of the original
 * and replacement triangulations. Their difference is affine on each overlay
 * polygon, hence the maximum norm occurs at one of those tested vertices.
 * This is a continuous piecewise-linear error bound, not a sample percentile.
 * Leaves form a restricted (2:1 edge-balanced) quadtree. The 16 transition
 * masks use the primal-quad triangulation templates in Zorin et al., SIGGRAPH
 * 2000 subdivision course notes, Fig. 5.5. These are inspection triangles,
 * NOT a claim of a conforming all-quadrilateral finite-element mesh. Quad
 * leaves and hanging-edge constraints are retained separately below.
 * A crop of a fixed model must query its fixed tree, not rebuild this LOD
 * on the crop: selection depends on the supplied geometry/support. */
typedef struct {
    int32_t x, y, cells, transition_mask; /* absolute lattice; bits B,R,T,L */
    int32_t corner[4];                   /* output vertex indices, CCW */
} RibbonLodQuad;

typedef struct {
    size_t input_faces, output_faces, patches, retained_faces, isolated_vertices;
    size_t transition_cases[16], balance_splits, error_splits;
    size_t nquads;
    RibbonLodQuad *quads;
    double max_error;
} RibbonLodReport;

int RibbonLod_build(Arena_T arena, const MeshBinData *input,
                     double du, double dv, int max_cells, double max_error,
                     MeshBinData *output, int32_t **source_vertex,
                     RibbonLodReport *report);
int RibbonLod_selftest(void);
/* Published transition connectivity, corners 0..3 CCW and mids B,R,T,L
 * numbered 4..7. Shared by inspection LOD and constrained quad-field output. */
int RibbonLod_transition(int mask, int triangles[6][3], int *count);
int RibbonLod_write(const char *input_path, const char *output_path,
                     double du, double dv, int max_cells, double max_error);

#endif
