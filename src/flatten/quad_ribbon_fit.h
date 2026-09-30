#ifndef QUAD_RIBBON_FIT_INCLUDED
#define QUAD_RIBBON_FIT_INCLUDED

#include "quad_field.h"
#include "../whole/axis_warp.h"

/* Global geometry fit of an established, claimant-separated observation
 * ribbon. This stage does NOT infer material identity from a winding label.
 * Complete supported quad disks enter the coupled Q1/P1 solve; cuts, holes,
 * ambiguous UV cells and exceptional triangles remain explicit support.
 * Original physical boundaries have exact Dirichlet data. Finest cells use
 * the source triangle diagonal, so refinement can represent source creases. */
typedef struct {
    QuadFieldReport solve;
    size_t input_faces, retained_faces, output_faces, output_vertices;
    size_t iteration, orientation_failures, unresolved_cells;
    size_t clearance_faces, exact_contact_faces, clearance_failures;
    size_t source_intersections, output_intersections, intersection_failures;
    size_t intersection_passes, intersection_refinements;
    double tolerance, continuous_bound, mesh_bound, rms, elapsed_seconds;
    double assembly_seconds, solve_seconds, validation_seconds;
    double intersection_seconds;
} QuadRibbonFitReport;

typedef struct {
    MeshBinData mesh;
    int32_t *reference_vertex; /* -1 for newly evaluated material samples */
    int32_t *reference_face;   /* -1 for an original isolated vertex */
    double *reference_barycentric; /* [nv*3], original observation face */
    int32_t *face_cell;        /* -1 for retained exceptional triangles */
    int32_t *reference_triangle; /* exact original triplet AND XYZ bits, else -1 */
} QuadRibbonFitOutput;

typedef int (*QuadRibbonFitCheckpoint)(void *context, int level,
    QuadField_T field, const QuadRibbonFitOutput *output,
    const QuadRibbonFitReport *report);

/* Output storage passed to checkpoint is temporary. A nonzero callback
 * result stops immediately. No changed source geometry is written. Existing
 * proper nonadjacent triangle crossings retain their original triangles.
 * New crossings in the actual emitted mesh refine the affected finite-element
 * cells before a level is accepted. This is not a coplanar-overlap, continuous
 * Q1 embedding, or material-correspondence certificate. */
int QuadRibbonFit_run(const MeshBinData *observations, double grid_step,
    QuadRibbonFitCheckpoint checkpoint, void *context);
int QuadRibbonFit_write(const char *input, const char *directory);
/* Read every surface member declared by <input_stem>_stats.json. Straightened
 * input is converted in memory through axis (NULL means already world frame).
 * A missing member is allowed only by the writer's explicit zero-face record
 * in <input_stem>_layers.json, never by a claimant multiplicity histogram.
 * Original inter-member clearance constrains each actual coarse/fine solve.
 * Storage colors and source sidecars remain separate; no correspondence is
 * inferred and no existing source contact is repaired or suppressed. */
int QuadRibbonFit_write_atlas(const char *input_stem, const char *directory,
                               const AxisWarp *axis);
int QuadRibbonFit_selftest(void);

#endif
