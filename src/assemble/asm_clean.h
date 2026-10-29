#ifndef ASM_CLEAN_INCLUDED
#define ASM_CLEAN_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "asm_types.h"

/* ============================================================================
 * asm_clean.h -- one cube mesh -> its charts, proactively cleaned.
 *
 * A chart is one connected component of the cube mesh, kept whole.  The
 * proactive cleaning is exactly what was measured elsewhere and nothing
 * more: geometrically double-sided prediction blobs go to extras;
 * tiny confetti goes to extras; non-manifold components and
 * components with handles are flagged (extras in this version, opened by a
 * tree-cotree cut later).  Neck splitting is NOT done here: the global
 * layout is the only reliable signal for that (see asm_cut).
 * ==========================================================================*/

typedef struct AsmCleanOpts {
    double blob_area_per_face;   /* legacy config/cache field; density never excludes a chart */
    double blob_double_sided;    /* fraction of vertices with the chart's own surface 1.5-8 vox along the normal */
    double min_chart_area;       /* charts below this are confetti */
    double stress_band;          /* flattening band for the SUSPECT flag */
    int    flatten_iters;        /* ARAP iteration cap */
} AsmCleanOpts;

typedef struct AsmCleanStats {
    size_t components;
    size_t blobs, tiny, nonmanifold, handles, flat_failed, suspect, kept;
    size_t verts_in, faces_in;
    size_t verts_kept, faces_kept, verts_blob, faces_blob;
    double area_in, area_kept, area_blob;
    double flatten_sec;
    size_t degenerate_split;     /* charts flattened without their zero-area faces */
    size_t degenerate_faces;     /* zero-area faces moved into excluded remainder charts */
    size_t crumple_split;        /* unflattenable charts rescued by excising their crumpled regions */
    size_t crumple_pieces;       /* certified pieces those charts became (in the layout) */
    size_t crumple_faces;        /* faces left in their excluded FLAT_FAILED remainder charts */
    double crumple_area_excluded;
} AsmCleanStats;

void AsmClean_default_opts(AsmCleanOpts *o);

/* Split `mesh` (one cube, world frame) into charts.  Chart arrays are
 * allocated from `persist` (they outlive the call); `scratch` is saved and
 * restored.  Every chart is flattened unless excluded.  `cube` is the pile
 * index stamped into each record.  out_charts receives a persist-allocated
 * array of n charts (ids are assigned later by the run).  A component that
 * fails to flatten only because of zero-area faces (a mesher slit closed by
 * collinear triangles) is kept as the chart without them when that flattens;
 * the zero-area faces follow the components as an excluded TINY remainder
 * chart of the same component, so every source face stays accounted for.
 * A component no seed flattens is cut down to its certified pieces around
 * its crumpled regions (AsmFlatten_rescue_crumpled): the largest keeps the
 * component's slot, the others follow the components in the layout, and the
 * excised faces follow as one excluded FLAT_FAILED chart.  Returns 0. */
int AsmClean_cube(Arena_T scratch, Arena_T persist,
                  const MeshBinData *mesh, int32_t cube,
                  const AsmCleanOpts *opts,
                  AsmChart **out_charts, size_t *out_n,
                  AsmCleanStats *stats);

/* Boundary mask of a chart-local mesh: mark[nv] = 1 where the vertex touches
 * an edge with exactly one incident face.  Scratch is restored. */
void AsmClean_boundary_mask(Arena_T scratch, size_t nv,
                            const int32_t *faces, size_t nf, uint8_t *mark);

int AsmClean_selftest(void);

#endif
