#ifndef SCROLL_MODEL_INCLUDED
#define SCROLL_MODEL_INCLUDED

#include <stdint.h>

#include "../common/arena.h"
#include "../common/mesh_pile.h"
#include "../common/mesh_bin.h"
#include "../whole/axis_warp.h"

typedef struct ScrollModel_T *ScrollModel_T;

typedef struct {
    const AxisWarp *axis;
    double axis_y, axis_x, pitch, core_radius;
    int sense;
    long chunk;
} ScrollModelOptions;

/* Build from the ENTIRE fixed source catalog. Area/query bounds deliberately
 * do not enter this interface. The working mesh is bounded to one evidence
 * region; the global unknowns are source charts, never fine mesh vertices.
 * Existing caches are validated against source bytes and model options.
 * Nonzero means no model is available; incomplete artifacts are diagnostic. */
int ScrollModel_build(Arena_T arena, const char *source_dir,
                       const char *model_dir, const ScrollModelOptions *options,
                       ScrollModel_T *out);

/* Read immutable global coefficients and source identities, verifying their
 * payload checksum. No graph solve, query statistics or re-gauging occurs.
 * The model may still be physically unqualified; loading is not acceptance. */
int ScrollModel_load(Arena_T arena, const char *model_dir, ScrollModel_T *out);
/* Release the native shaft owned by a loaded physical-coordinate model before
 * disposing/restoring its arena. Built models borrow their options axis. */
void ScrollModel_release(ScrollModel_T model);
uint64_t ScrollModel_coordinate_fingerprint(ScrollModel_T model);
int ScrollModel_validate_catalog(Arena_T arena, ScrollModel_T model,
                                  const char *source_dir);
/* Experimental global coarse-to-fine geometric fit. No area bounds enter
 * this API. Every level is a persisted quad coefficient field and mesh;
 * source topology/physical parameterization still require qualification. */
int ScrollModel_fit_quads(Arena_T arena, ScrollModel_T model, const char *output_dir);
/* Area diagnostics separate the untrimmed sparse control mesh from a dense
 * source-support sampling of its fixed coefficients. Both restrict without
 * solving; neither is a physically qualified, coarsened trimmed ribbon yet. */
int ScrollModel_write_quad_area(Arena_T arena, ScrollModel_T model,
                                 const char *field_path, const long *cube_bounds,
                                 const char *output_dir);

/* Source-order evaluation; no parent or query-dependent statistic. Caller
 * provides this source cube's ORIGINAL geometry and receives absolute UV,
 * winding, and material chart identity. All output arrays are arena-owned. */
int ScrollModel_evaluate_cube(Arena_T arena, ScrollModel_T model,
                              size_t cube, const MeshBinData *source_mesh,
                              float **out_uv, float **out_q,
                              int32_t **out_material);
size_t ScrollModel_cube_count(ScrollModel_T model);
const MeshPileEntry *ScrollModel_cube(ScrollModel_T model, size_t cube);

/* Diagnostic source-order query at the original sample resolution. Bounds
 * select inclusive cube origins z0,z1,y0,y1,x0,x1; NULL selects the catalog.
 * Streams original geometry and fixed UVs with stable source identities.
 * This does not fit, repair, pack or qualify an area independently. */
int ScrollModel_write_area(Arena_T arena, ScrollModel_T model,
                            const long *cube_bounds, const char *output_dir);

#endif
