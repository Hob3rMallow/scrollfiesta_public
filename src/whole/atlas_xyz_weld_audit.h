#ifndef ATLAS_XYZ_WELD_AUDIT_INCLUDED
#define ATLAS_XYZ_WELD_AUDIT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../remesh/ball_pivot.h"
#include "../unroll/piece_set.h"

/* Held-out physical-topology oracle.  Topology_build invokes the original
 * XYZ-only BPA seam weld; no UV quantity is accepted by this step. */
typedef struct {
    const float *vertices;       /* refined XYZ vertices used by the weld */
    size_t nvertices;
    size_t source_vertices;      /* first source_vertices belong to PieceSet */
    const int32_t *parent0;      /* [nvertices-source_vertices] midpoint parents */
    const int32_t *parent1;
    const int32_t *faces;       /* [nfaces*3], indices into vertices above */
    size_t nfaces;
    size_t refined_faces;
    size_t refined_vertices_added;
    size_t combined_faces;
} AtlasXyzWeldTopology;

typedef struct {
    size_t bridge_faces;
    size_t bridge_vertices;
    size_t bridge_edges;
    size_t cross_cube_bridge_edges;

    size_t incident_chart_positions;
    size_t split_bridge_vertices;
    double incident_u_spread_median;
    double incident_u_spread_p95;
    double incident_u_spread_max;

    size_t exploded_cross_cube_edges;
    size_t compressed_cross_cube_edges;
    double best_cross_cube_edge_stretch_median;
    double best_cross_cube_edge_stretch_p95;
    double best_cross_cube_edge_stretch_max;
    double best_cross_cube_edge_abs_error_median;
    double best_cross_cube_edge_abs_error_p95;
    double best_cross_cube_edge_abs_error_max;

    size_t exploded_bridge_faces;
    size_t degenerate_bridge_faces;
    double best_bridge_face_symmetric_stretch_median;
    double best_bridge_face_symmetric_stretch_p95;
    double best_bridge_face_symmetric_stretch_max;

    size_t embedded_bridge_vertices;
    size_t discontinuous_embedded_bridge_vertices;
    double embedded_u_spread_median;
    double embedded_u_spread_p95;
    double embedded_u_spread_max;

    size_t exact_overlap_pairs;
    size_t base_base_overlap_pairs;
    size_t bridge_base_overlap_pairs;
    size_t bridge_bridge_overlap_pairs;
    size_t broad_phase_records;
    size_t broad_phase_cells;
    size_t broad_phase_candidate_pairs;
    double broad_phase_cell_size;
    int broad_phase_complete;
} AtlasXyzWeldAuditStats;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    size_t cross_cube_edges;
    double total_xyz_edge_length;
} AtlasXyzWeldConnection;

typedef struct {
    size_t cross_cube_edges;
    size_t same_chart_edges;
    size_t cross_chart_edges;
    size_t ambiguous_chart_edges;
    size_t relations;
} AtlasXyzWeldConnectionStats;

/* One soft gauge equation contributed by the XYZ-only weld topology:
 *
 *     shift[chart1] - shift[chart0] = target_shift_median.
 *
 * The target is the median endpoint U mismatch over all unambiguous bridge
 * edges between the chart pair.  The equation remains soft; cross_cube_edges
 * is its evidence weight.  A single BPA edge must never hard-union two atlas
 * gauges. */
typedef struct {
    int32_t chart0;
    int32_t chart1;
    size_t cross_cube_edges;
    double total_xyz_edge_length;
    double target_shift_median;
    double target_shift_mad;
} AtlasXyzWeldShiftConstraint;

/* Triangle soup used by the overlap oracle.  The first base_faces triangles
 * are the chart-split atlas and the remainder are the XYZ/BPA bridge faces
 * embedded at their best available incident-chart positions.  Storage belongs
 * to the arena passed to AtlasXyzWeldAudit_evaluate_with_mesh. */
typedef struct {
    const int32_t *faces;
    const double *u;
    const double *v;
    const int32_t *component;
    size_t nvertices;
    size_t nfaces;
    size_t base_faces;
    size_t bridge_faces;
} AtlasXyzWeldAuditMesh;

int AtlasXyzWeldTopology_build(
    Arena_T arena, const PieceSet *ps, float cube_size, float rho,
    float rho_max, float band, const BpaBridgeGate *gate,
    AtlasXyzWeldTopology *out);

/* Evaluate one rigid chart placement against a topology already generated
 * from XYZ.  A bridge corner may use any incident chart copy of its source
 * vertex; the reported face/edge distortion is therefore a best-case bound.
 * The selected-copy continuity metric prevents that local choice from hiding
 * mutually inconsistent embeddings across adjacent bridge faces. */
int AtlasXyzWeldAudit_evaluate(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const double *base_u, const double *v, const int32_t *face_chart,
    size_t ncharts, const double *chart_shift,
    AtlasXyzWeldAuditStats *stats);

int AtlasXyzWeldAudit_evaluate_with_mesh(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const double *base_u, const double *v, const int32_t *face_chart,
    size_t ncharts, const double *chart_shift,
    AtlasXyzWeldAuditStats *stats, AtlasXyzWeldAuditMesh *mesh);

/* Aggregate the UV-independent weld topology by rigid-chart pair.  A source
 * or refined vertex incident to more than one chart is deliberately excluded:
 * such an edge cannot provide an unambiguous block-coherence observation. */
int AtlasXyzWeldTopology_collect_connections(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const int32_t *face_chart, size_t ncharts,
    AtlasXyzWeldConnection **out_connection, size_t *out_nconnection,
    AtlasXyzWeldConnectionStats *stats);

/* Aggregate the same UV-independent bridge topology into soft gauge equations
 * for a particular incoming field.  base_u is sampled at source vertices and
 * propagated through seam-refinement midpoint provenance. */
int AtlasXyzWeldTopology_collect_shift_constraints(
    Arena_T arena, const PieceSet *ps, const AtlasXyzWeldTopology *topology,
    const double *base_u, const int32_t *face_chart, size_t ncharts,
    AtlasXyzWeldShiftConstraint **out_constraint, size_t *out_nconstraint,
    AtlasXyzWeldConnectionStats *stats);

#endif
