#ifndef TIFXYZ_WRITE_INCLUDED
#define TIFXYZ_WRITE_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * tifxyz_write.h -- rasterize a UV-parameterized mesh into the volume-
 * cartographer "tifxyz" surface format.
 *
 * A tifxyz is a directory:
 *   x.tif, y.tif, z.tif  -- 32-bit float, 1 sample/px, uncompressed; each grid
 *                           cell holds one world coordinate; invalid = -1.0.
 *   meta.json            -- { bbox, type:"seg", uuid, format:"tifxyz", scale }.
 *
 * IMPORTANT coordinate order: meshes here are (z,y,x); tifxyz/VC are (x,y,z).
 * We reorder on write: x.tif <- vert[..+2], y.tif <- vert[..+1], z.tif <- vert[..+0].
 * ==========================================================================*/

typedef struct {
    float px_per_vox;    /* grid cells per voxel, both axes (<=0 => default).
                          * Acts like vc_obj2tifxyz's stretch_factor. */
    int   max_dim;       /* cap on grid width/height (0 => default) */
    float contested_tol_vox; /* invalidate multi-covers farther apart than this
                              * (<=0 => 1.5 vox) */
} TifXYZOpts;

typedef struct {
    int    width, height;
    double scale_u, scale_v;   /* px per unwrapped-unit (meta.json "scale") */
    size_t cells_filled;
    size_t cells_contested;    /* invalidated multi-cover cells */
    double bbox_lo[3];         /* world x,y,z (tifxyz order) */
    double bbox_hi[3];
} TifXYZStats;

/* Rasterize and write the tifxyz directory.
 * verts: [nv*3] (z,y,x). uv: [nv*2] (u,v). faces: [nf*3] 0-based.
 * out_dir is created (with parents) if needed; uuid goes into meta.json.
 * out_stats may be NULL. Returns 0 on success, -1 on failure. */
int TifXYZ_write(Arena_T arena, const char *out_dir, const char *uuid,
                 const float *verts, size_t nv,
                 const int32_t *faces, size_t nf,
                 const float *uv, const TifXYZOpts *opts,
                 TifXYZStats *out_stats);

int TifXYZ_write_double(Arena_T arena, const char *out_dir, const char *uuid,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const double *uv, const TifXYZOpts *opts,
                        TifXYZStats *out_stats);

/* The writer also emits mask.tif (255 valid, 0 invalid). UV cells covered by
 * geometrically different surface points are invalidated rather than silently
 * keeping the first triangle; this is the Villa-safe contested-pixel contract.
 *
 * In-process unit tests (float-TIFF round-trip + tifxyz rasterization/reorder).
 * Returns 0 if all pass, else number of failures. Logs to stderr. */
int TifXYZ_selftest(void);

#endif
