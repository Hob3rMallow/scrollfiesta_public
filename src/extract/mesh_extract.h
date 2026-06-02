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
 */
/* In-memory input (optional): when vol_in != NULL the TIFF / halo loaders are
 * bypassed. vol_in is a padded (p_size_in)^3 uint8 buffer whose index (0,0,0)
 * is world voxel (cube_origin - halo_voxels), exactly like HaloLoader_load's
 * output; cube_origin_in is the owned cube origin (z,y,x). This lets a caller
 * stream a cube straight from a remote zarr instead of reading a TIFF. */
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
                    const uint8_t   *vol_in,          /* optional in-memory cube */
                    int              p_size_in,
                    const int64_t   *cube_origin_in,  /* owned origin (z,y,x) */
                    ComponentMesh  **out_meshes,
                    size_t          *out_n_meshes);

#endif
