#ifndef QUAD_FIELD_INCLUDED
#define QUAD_FIELD_INCLUDED

#include "../common/arena.h"
#include "../common/mesh_bin.h"

/* Q1 finite elements on an explicitly supplied, edge-balanced dyadic quad
 * domain. Coordinates are absolute integer lattice addresses; chart identifies
 * an injective FITTING DOMAIN, not merely a material label (different turns
 * of the same sheet must not alias). No UV rounding or welding of different
 * chart identities occurs. The caller owns support/topology. */
typedef struct { int32_t chart, x, y, size; } QuadFieldCell;
typedef struct QuadField_T *QuadField_T;
typedef struct {
    size_t cells, nodes, unknowns, hanging_nodes, samples;
    double maximum_error, rms_error, relative_residual;
    size_t fixed_nodes;
    size_t node_storage_bytes;
    double observation_weight; /* actual quadrature weight, for bank RMS */
} QuadFieldReport;
typedef struct { int32_t chart, x, y; int constrained; } QuadFieldNodeKey;
typedef struct {
    size_t vertices, faces, excluded_vertices, excluded_faces;
    size_t reversed_faces, collapsed_faces;
    double maximum_error, rms_error;
} QuadFieldProjectionReport;

int QuadField_new(Arena_T arena, const QuadFieldCell *cells, size_t count,
                   double lattice_step, QuadField_T *out);
/* Accumulate complete source evidence in stable source order. Cell identity
 * is explicit: the caller locates the point on the immutable global tree.
 * This interface never allocates a fine raster or stores individual samples. */
int QuadField_observe(QuadField_T field, size_t cell, double u, double v,
                       const double xyz[3], double weight);
/* Global constrained least squares, factor once / solve XYZ. Coarse values
 * initialize the fine correction solve, never freeze fine unknowns. The
 * factorization includes only free masters; exact source values are eliminated
 * and scattered back without being counted against the unknown limit. */
int QuadField_solve(QuadField_T field, QuadField_T coarse,
                     QuadFieldReport *report);
/* Exact Dirichlet data on original support boundaries, NOT frozen coarse
 * iterates. Constrained nodes cannot be fixed independently: the caller must
 * refine the incident coarse leaf first. Repeated inconsistent targets fail. */
size_t QuadField_node_count(QuadField_T field);
int QuadField_node_key(QuadField_T field, size_t node, QuadFieldNodeKey *key);
int QuadField_find_node(QuadField_T field, int32_t chart, int32_t x, int32_t y,
                         size_t *node);
int QuadField_fix_node(QuadField_T field, size_t node, const double xyz[3]);
/* Finest-cell source triangles: 0 = bilinear Q1, 1 = P1 diagonal 0--2,
 * 2 = P1 diagonal 1--3 (corners CCW from lower left). Nonzero kinds require
 * size==1, hence no finer hanging edges. Set before accumulating evidence. */
int QuadField_set_element(QuadField_T field, size_t cell, int kind);
int QuadField_element(QuadField_T field, size_t cell);
/* Optional assembly cache for an unchanged observation population. Copies
 * sufficient statistics only where cell key and element basis are identical;
 * never copies constraints or freezes coefficients. Other readers need not
 * use this optimization. Call before adding observations to the new field. */
int QuadField_reuse_observations(QuadField_T field, QuadField_T previous);
size_t QuadField_observation_count(QuadField_T field, size_t cell);
int QuadField_evaluate(QuadField_T field, size_t cell, double u, double v,
                        double xyz[3]);
int QuadField_locate(QuadField_T field, int32_t chart, double u, double v,
                      size_t *cell);
int QuadField_mesh(Arena_T arena, QuadField_T field, MeshBinData *out);
/* Restriction preserves the global transition templates and coefficient
 * bits. selected==NULL exports every leaf. source_node maps output vertices
 * to the fixed field's node identities, independent of the selected area. */
int QuadField_mesh_selection(Arena_T arena, QuadField_T field,
                              const uint8_t *selected, MeshBinData *out,
                              int32_t **source_node);
/* Support-preserving DIAGNOSTIC: evaluate the fixed field at source vertices
 * and retain their exact UVs and connectivity. Negative chart identities are
 * excluded explicitly (the caller retains their geometry as extras). No face
 * is introduced across a hole, cut or separate source component. This is a
 * dense PL sampling of Q1, not a coarsened/trimmed quad tessellation, and not
 * proof that the fitted 3D embedding is physically correct. */
int QuadField_project_mesh(Arena_T arena, QuadField_T field,
                            const MeshBinData *source, const int32_t *chart,
                            MeshBinData *out, int32_t **source_vertex,
                            QuadFieldProjectionReport *report);
/* Refine marked leaves and close the 2:1 edge balance. No existing domain
 * cell is deleted; material identity and absolute lattice are unchanged. */
int QuadField_refine(Arena_T arena, QuadField_T field, const uint8_t *marked,
                      QuadField_T *out);
size_t QuadField_cell_count(QuadField_T field);
const QuadFieldCell *QuadField_cell(QuadField_T field, size_t cell);
/* Immutable, checksummed global coefficients. The source fingerprint binds
 * every saved level to its source observations; loading never solves and uses
 * the constructor's geometry/storage limits, not a free-factor unknown cap. */
int QuadField_save(QuadField_T field, const char *path, uint64_t source_fingerprint);
int QuadField_load(Arena_T arena, const char *path, uint64_t source_fingerprint,
                    QuadField_T *out);
/* Reader visits the same fixed full-source chunks on every pass. Chunk
 * storage is scratch-arena-owned. Checkpoint receives a genuinely solved
 * coarse/fine field, before the next adaptive refinement. */
typedef int (*QuadFieldReader)(void *context, Arena_T scratch, size_t chunk,
                                MeshBinData *mesh, int32_t **chart);
typedef int (*QuadFieldCheckpoint)(void *context, int level, QuadField_T field,
                                    const QuadFieldReport *report);
int QuadField_fit_stream(Arena_T arena, size_t chunks, QuadFieldReader reader,
                          void *reader_context, int root_cells,
                          double lattice_step, double tolerance,
                          QuadFieldCheckpoint checkpoint, void *checkpoint_context,
                          QuadField_T *out);
int QuadField_selftest(void);

#endif
