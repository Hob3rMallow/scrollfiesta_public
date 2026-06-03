#ifndef MESH_RESPLIT_INCLUDED
#define MESH_RESPLIT_INCLUDED

#include <stddef.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"
#include "mesh_extract.h"   /* MeshResplitCloud */

/* Optional diagnostic dump target. When non-NULL with a stage name set, the
 * split sub-meshes (cc_stage) and/or the re-LOP'd pre-BPA point clouds
 * (mls_stage) are written into <dir>/<cube_id>_<stage>/ (world coords). NULL or
 * a NULL stage skips that dump. `dir` is the cube dump dir "<dump>/<cube_id>". */
typedef struct {
    const char *dir;
    const char *cube_id;
    const char *cc_stage;
    const char *mls_stage;
} MeshResplitDump;

/*
 * Post-extract connectivity re-split (step between Step 0 extract and Step 2).
 *
 * Extract emits one ComponentMesh per *voxel* connected-component. A voxel-CC
 * can contain several *disconnected surface sheets* surfaced together; because
 * the initial LOP runs over the whole voxel-CC cloud, those close sheets pull on
 * each other and leave cross-sheet artifacts that make the Step-2 oracle
 * false-positive and BridgeCut over-split.
 *
 * For each input mesh: split it into its actual MESH connectivity-components.
 *   - one real connectivity-component  -> pass through unchanged;
 *   - >= 2 real connectivity-components -> re-split: label each ORIGINAL
 *     voxel-center point (clouds[i]) by which component it belongs to (nearest
 *     final vertex of its LOP image), then re-LOP + re-BPA each component's
 *     ORIGINAL points IN ISOLATION, so the merged-LOP artifacts are gone and the
 *     sheets enter Step 2 as clean single surfaces.
 *
 * A "real" component has >= RESPLIT_MIN_COMP_VERTS vertices; smaller fragments
 * are dropped. If a component's re-LOP/re-BPA fails, its original geometry is
 * emitted unchanged (never lose a sheet).
 *
 * clouds: index-aligned with `in` (cloud k describes mesh k); may be NULL, in
 * which case everything passes through (no re-split possible without originals).
 *
 * On success returns 0; *out is an arena-allocated ComponentMesh array of length
 * *n_out (>= n_in possible). comp_id is reassigned 1..*n_out.
 *
 * out_clouds (optional): per output mesh, the labelled original-point cloud that
 * produced it (index-aligned with *out). A later split (bridge-cut) re-LOPs its
 * pieces from these via MeshResplit_remesh_pieces. NULL to skip.
 */
int MeshResplit_run(Arena_T arena,
                    const ComponentMesh   *in, size_t n_in,
                    const MeshResplitCloud *clouds,
                    int                    n_threads,
                    const MeshResplitDump *dump,
                    ComponentMesh        **out, size_t *n_out,
                    MeshResplitCloud     **out_clouds);

/*
 * Re-LOP-from-original after ANY split. Given a parent component's original-point
 * cloud and the sub-meshes it was just split into (`pieces`, e.g. bridge-cut /
 * overlap output), label each original point to its nearest piece (by the LOP
 * image), then re-LOP + re-BPA each piece from ITS original points in isolation
 * (fall back to the piece unchanged if re-LOP/BPA fails). This is the same clean
 * re-mesh MeshResplit_run applies after the connectivity split, reused for step 4.
 *
 * cloud: the parent's MeshResplitCloud (may be NULL/empty -> pieces pass through).
 * out / out_clouds (optional): index-aligned re-meshed pieces and their labelled
 *   sub-clouds. dump_tag (or NULL) names the env-gated diagnostic point dumps.
 * Returns 0; *out arena-allocated, length *n_out (== n_pieces).
 */
int MeshResplit_remesh_pieces(Arena_T arena,
                              const MeshResplitCloud *cloud,
                              const ComponentMesh    *pieces, size_t n_pieces,
                              const char             *dump_tag,
                              const MeshResplitDump  *dump,
                              ComponentMesh         **out, size_t *n_out,
                              MeshResplitCloud      **out_clouds);

#endif
