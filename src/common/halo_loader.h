#ifndef HALO_LOADER_INCLUDED
#define HALO_LOADER_INCLUDED

#include <stdint.h>
#include "arena.h"

/*
 * HaloLoader -- load a padded volume from a cube's prediction TIFF plus its
 * up-to-26 neighbors (6 face + 12 edge + 8 corner), to support cross-cube
 * mesh extraction with overlap stitching.
 *
 * Reads <pred_dir>/<cube_id>.tif as the central cube and any of its 26
 * grid neighbors that exist as <pred_dir>/<neighbor_id>.tif, where each
 * neighbor_id is the cube_id with origin offset by (+/- cube_size) on
 * one or more axes. Missing neighbors (file absent or off-grid) leave
 * their region in the padded buffer zero-filled — equivalent to
 * "background outside the cube grid".
 *
 * Returns 0 on success, -1 on failure (cannot load own cube TIFF).
 *
 *   out_vol     : arena-allocated uint8[p_size^3]. Padded volume; the
 *                 central cube_size^3 region is the cube itself, and the
 *                 surrounding halo_voxels-thick shell holds neighbor data.
 *   out_p_size  : cube_size + 2*halo_voxels.
 *   out_origin  : (vz - halo, vy - halo, vx - halo) source-space coordinate
 *                 of padded_vol[0,0,0].
 */
int HaloLoader_load(Arena_T arena,
                    const char *pred_dir,
                    const char *cube_id,
                    int cube_size,
                    int halo_voxels,
                    uint8_t **out_vol,
                    int *out_p_size,
                    int64_t out_origin[3]);

#endif
