#ifndef MESH_EXTRACT_INCLUDED
#define MESH_EXTRACT_INCLUDED

#include <stdint.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

/* ----------------------------------------------------------------
 * Exposed primitives (used by Step 0 itself).
 * ---------------------------------------------------------------- */

typedef struct {
    int32_t label;
    int32_t count;
} CompInfo;

/* Per-output-component snapshot of the raw + LOP'd point cloud, captured by
 * MeshExtract_run when out_clouds is non-NULL, for the post-extract connectivity
 * re-split (see src/extract/mesh_resplit.h). The arrays are index-aligned with
 * the out_meshes array (cloud k describes mesh k). n == 0 means "no snapshot". */
typedef struct {
    const float *orig_pts;       /* [n*3] pre-LOP voxel centers (cube-local z,y,x) */
    const float *lop_pts;        /* [n*3] their 1:1 LOP images (pre-weld)          */
    float        cell_origin[3]; /* cube_world_origin the LOP used (hash alignment)*/
    size_t       n;
} MeshResplitCloud;

/* 3D 6-connected connected components via BFS.
 * Writes labels[D*H*W] (0 = background, 1..n_comps = foreground labels).
 * Allocates and returns *out_comps as a sorted-by-count-desc array of length
 * *out_n_comps. All allocations are arena-backed. */
int cc_label_3d(Arena_T arena, const uint8_t *vol,
                int D, int H, int W,
                int32_t *labels,
                CompInfo **out_comps, int32_t *out_n_comps);

/*
 * MeshExtract_run -- Step 0 entry point.
 *
 * Loads the TIFF at tiff_path, extracts connected components, samples
 * voxel-center points, projects them to the sheet midline (MLS/LOP),
 * orients normals (Hoppe MST), triangulates with Ball-Pivoting, then
 * cleans up and trims.
 *
 * If halo_voxels > 0, loads a padded volume of size (cube_D + 2*halo)^3
 * via HaloLoader from <pred_dir>/<cube_id>.tif and its 26 grid neighbors;
 * vertex coordinates returned are in cube-local space (owned region at
 * [0, cube_D), halo regions at [-halo, 0) and [cube_D, cube_D + halo)).
 * Each ComponentMesh's pin_mask is allocated and set to 1 for verts in
 * halo region, 0 for verts in owned region. Components whose bounding
 * box lies entirely outside the owned region are dropped (they belong
 * to neighbor cubes).
 *
 * If halo_voxels == 0: tiff_path is loaded directly (legacy behavior);
 * pin_mask is NULL on output meshes; pred_dir is unused.
 *
 * dump_cube_dir / cube_id: if both non-NULL, writes pre- and post-QEM
 * OBJ meshes into {dump_cube_dir}/{cube_id}_step0_pre_qem/ and
 * {dump_cube_dir}/{cube_id}_step0/ respectively.
 *
 * Returns 0 on success, nonzero on failure.
 * On success, *out_meshes points to an arena-allocated array of
 * ComponentMesh structs, and *out_n_meshes is the count.
 *
 * out_clouds: optional (pass NULL to skip). If non-NULL, receives an
 * arena-allocated array of MeshResplitCloud (length *out_n_meshes, index-aligned
 * with *out_meshes) holding each component's pre-LOP voxel-center points and
 * their 1:1 LOP images, for the post-extract connectivity re-split. Adds a small
 * per-component arena snapshot; harmless when NULL.
 */
int MeshExtract_run(Arena_T          arena,
                    const char      *tiff_path,
                    const char      *pred_dir,        /* used when halo_voxels > 0 */
                    int              halo_voxels,
                    size_t           cube_D,
                    size_t           cube_H,
                    size_t           cube_W,
                    int              n_threads,
                    const char      *dump_cube_dir,
                    const char      *cube_id,
                    int              skip_qem,
                    ComponentMesh  **out_meshes,
                    size_t          *out_n_meshes,
                    MeshResplitCloud **out_clouds);

#endif
