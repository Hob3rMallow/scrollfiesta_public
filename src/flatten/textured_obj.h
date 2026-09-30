#ifndef TEXTURED_OBJ_INCLUDED
#define TEXTURED_OBJ_INCLUDED

#include "rawtex_bake.h"

/* Presentation-only export. VMESH/scientific UVs are never modified.
 * Reads prefix_rawtex.tif using bounded row access; writes 3D/flat OBJs,
 * a shared MTL and sparse overlapping texture tiles at the bake's EXACT
 * pixel density. Face topology/XYZ are unchanged, UV indices may differ
 * from vertex indices at material boundaries. No CT vertex-color sampling.
 * Faces outside a fixed bake window remain, with explicit no-data material.
 * Returns -1 rather than downsampling or silently clipping an oversized tile. */
int TexturedObj_write(const char *prefix, const float *verts, const float *uv,
                       const float *normals, size_t nv,
                       const int32_t *faces, size_t nf, const RawtexPlan *plan,
                       double du, double dv, int do_3d, int do_flat);

int TexturedObj_write_field(const char *prefix, const float *verts, RawtexUv uv,
                       const float *normals, size_t nv,
                       const int32_t *faces, size_t nf, const RawtexPlan *plan,
                       double du, double dv, int do_3d, int do_flat);
int TexturedObj_selftest(void);
/* Existing full-resolution TIFF to exact PNG tiles and a navigation pyramid. */
int TexturedObj_pyramid(const char *input,const char *output_prefix);

#endif
