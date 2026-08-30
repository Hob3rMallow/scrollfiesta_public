#ifndef TOPOLOGY_INVARIANTS_INCLUDED
#define TOPOLOGY_INVARIANTS_INCLUDED

/*
 * topology_invariants.h -- component-wise topological invariants for triangle
 * meshes, plus a minimum-cardinality H_1 generator basis.
 *
 * All Betti numbers are over F_2.  For an orientable surface H_1 is
 * torsion-free, so these are also the ordinary integer ranks.  A component is
 * certified homeomorphic to a disk only after the edge- and vertex-manifold,
 * boundary, orientability, Euler, and homology checks all agree.
 *
 * Generator representatives use a shortest-path-tree / maximum-dual-cotree
 * decomposition.  There are exactly beta_1 loops (the minimum possible
 * cardinality).  Within the chosen rooted shortest-path tree, the cotree is
 * selected to minimize the total rooted fundamental-loop weight.  This is not
 * advertised as a globally minimum-length homology basis over every possible
 * root and primal tree.
 */

#include <stddef.h>
#include <stdint.h>

enum {
    TOPOLOGY_DEFECT_NONMANIFOLD_EDGE       = 1u << 0,
    TOPOLOGY_DEFECT_NONMANIFOLD_VERTEX     = 1u << 1,
    TOPOLOGY_DEFECT_IRREGULAR_BOUNDARY     = 1u << 2,
    TOPOLOGY_DEFECT_NONORIENTABLE          = 1u << 3,
    TOPOLOGY_DEFECT_BOUNDARY_COUNT         = 1u << 4,
    TOPOLOGY_DEFECT_NONTRIVIAL_H1          = 1u << 5,
    TOPOLOGY_DEFECT_NONTRIVIAL_H2          = 1u << 6,
    TOPOLOGY_DEFECT_INCONSISTENT_WINDING   = 1u << 7,
    TOPOLOGY_DEFECT_GENERATOR_FAILURE      = 1u << 8,
    TOPOLOGY_DEFECT_ALGEBRA_MISMATCH       = 1u << 9
};

typedef struct {
    /* Build explicit loop representatives.  Betti numbers and the exact
     * minimum generator count are always computed, even when this is zero. */
    int emit_generators;

    /* 0 emits every generator.  A positive value keeps the shortest N
     * representatives per component while retaining the exact full rank. */
    size_t max_generators_per_component;
} TopologyAuditOptions;

typedef struct {
    size_t rank;                  /* components are sorted by face count */
    int32_t root_vertex;          /* deterministic shortest-path-tree root */

    size_t vertices, edges, faces;
    /* Geometry summary for provenance matching between pipeline stages.
     * has_geometry is zero only when TopologyAudit_analyze received verts=NULL. */
    int has_geometry;
    double bbox_min[3], bbox_max[3];
    double centroid[3];          /* arithmetic mean of referenced vertices */
    size_t boundary_edges;
    /* Connected components of the boundary graph; genuine loops exactly when
     * boundary_irregular_vertices is zero. */
    size_t boundary_loops;
    size_t nonmanifold_edges, nonmanifold_vertices;
    size_t same_direction_edges;
    size_t boundary_irregular_vertices;

    int64_t euler_characteristic;
    int orientable;
    int input_winding_consistent;
    int surface_valid;

    /* Exact F_2 Betti numbers of the triangle complex.  These remain defined
     * for pinched and non-manifold components; surface_valid is independent. */
    int64_t beta_0, beta_1, beta_2;

    /* Exactly one of orientable_genus / crosscap_number is meaningful for a
     * valid connected surface.  Unknown values are -1. */
    int64_t orientable_genus;
    int64_t crosscap_number;

    int homeomorphic_to_disk;
    uint32_t defect_mask;

    /* beta_1 is the mathematical minimum number of H_1 generators.
     * generator_basis_complete verifies the tree-cotree construction found
     * exactly that many independent representatives. */
    size_t minimal_generator_rank;
    size_t emitted_generators;
    int generator_basis_complete;
} TopologyComponentInvariant;

typedef struct {
    size_t component;             /* index into component[] */
    size_t ordinal;               /* length order within the component */
    int32_t closing_edge_a;
    int32_t closing_edge_b;

    double length;                /* reduced simple fundamental cycle */
    double rooted_length;         /* root->a + edge + b->root closed walk */

    /* Ordered simple cycle vertices.  The closing edge joins the last vertex
     * back to the first; the first vertex is not duplicated at the end. */
    size_t nvertices;
    int32_t *vertices;
} TopologyGenerator;

/* One connected component of the boundary graph.  For a valid surface every
 * entry is a simple cycle.  Keeping the irregular entries here as well is
 * useful diagnostically: a branched boundary is precisely where a nominal
 * "hole" cannot safely be handed to a polygon filler. */
typedef struct {
    size_t rank;                  /* global order: component, then size */
    size_t component;             /* index into component[] */
    size_t ordinal;               /* order within that component */
    int32_t root_vertex;          /* minimum vertex in this boundary CC */

    size_t vertices, edges;
    size_t irregular_vertices;    /* boundary degree != 2 */
    int simple_cycle;
    int component_perimeter;      /* longest boundary CC on this component */

    int has_geometry;
    double bbox_min[3], bbox_max[3];
    double length;                /* sum of boundary-edge lengths */
    double max_edge_length;
    /* Exact maximum vertex-pair distance for small boundary components.
     * Large entries carry their bbox diagonal as an upper bound instead. */
    double diameter;
    int diameter_exact;
} TopologyBoundaryInvariant;

typedef struct {
    size_t vertices, edges, faces;
    size_t face_components;
    size_t isolated_vertices;
    size_t topological_components;

    size_t disk_components;
    size_t nondisk_components;
    size_t invalid_surface_components;
    size_t nonorientable_components;
    size_t inconsistent_winding_components;

    int64_t euler_characteristic;
    int64_t beta_0, beta_1, beta_2;
    int betti_complete;
    int all_components_are_disks;

    size_t minimal_generator_rank;
    size_t emitted_generators;
    int generators_truncated;
    int generator_basis_complete;

    TopologyComponentInvariant *component;
    TopologyGenerator *generator;

    size_t boundary_components;
    TopologyBoundaryInvariant *boundary;

    /* Dense component index for every input vertex, or -1 for an isolated
     * vertex.  The index addresses component[].  Exposing this correspondence
     * lets callers make topology decisions at chart granularity instead of
     * inferring provenance again after a vertex-fan split. */
    int32_t *vertex_component;

    char error[256];
} TopologyAuditReport;

void TopologyAudit_options_default(TopologyAuditOptions *options);

/* Analyze vertices/faces.  verts may be NULL, in which case every edge has
 * unit length.  Faces must be valid non-degenerate triangles with int32
 * indices.  Returns 0 on success and -1 with report->error filled on failure. */
int TopologyAudit_analyze(const float *verts, size_t nv,
                          const int32_t *faces, size_t nf,
                          const TopologyAuditOptions *options,
                          TopologyAuditReport *report);

void TopologyAudit_dispose(TopologyAuditReport *report);

/* Native fixtures: disk, annulus, sphere, punctured torus, disconnected
 * disks, bowtie vertex, non-manifold sparse rank, and Mobius strip. */
int TopologyAudit_selftest(void);

#endif
