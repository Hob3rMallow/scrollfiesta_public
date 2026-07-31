/* ============================================================================
 * rawtex_bake.h -- rasterize RAW CT intensity into an unrolled (u,v) TIFF.
 *
 * Extracted from obj_bake_raw.c so the unified unroll pipeline (scroll_unroll)
 * and the standalone obj_bake_raw tool share one rasterizer. Given a
 * UV-carrying mesh (verts in source-space z,y,x + per-vertex u,v) and a RAW
 * cube table, every output pixel center is located in its covering UV
 * triangle, the 3D point + normal barycentrically interpolated, and the volume
 * freshly sampled there (normal-max multi-tap). Never a per-vertex scatter.
 * ==========================================================================*/
#ifndef RAWTEX_BAKE_H
#define RAWTEX_BAKE_H

#include <stddef.h>
#include <stdint.h>

#include "../common/raw_sample.h"   /* CubeTable, sample_vertex */

/* Optional diagnostics rasterized in the SAME pass (diag != NULL):
 *   <prefix>_diagclass.tif   0 bg | 1 ok | 2 skip-uv | 3 skip-3d | 4 multi | 5 dark
 *   <prefix>_diagstretch.tif 128 + 42.5*log2(sigma_max)  (UV->3D compression)
 *   <prefix>_diagsquash.tif  128 + 42.5*log2(sigma_min)  (UV stretch)
 *   <prefix>_diagverr.tif    |vt.v - (z - zmin)| * 32     (axial-mapping error)
 * plus a worst-tile table on stderr; smear_obj (if set) dumps the skip-uv
 * faces in 3D + z/winding histograms. */
typedef struct DiagOpts {
    const char *prefix;        /* output path prefix */
    int         dark_thresh;   /* uint8 post-stretch dark cutoff (class 5) */
    const char *smear_obj;     /* if set: write skip-uv faces in 3D + histograms */
} DiagOpts;

/* Rasterize the (u,v) texture at (du,dv) vox/px and write `path` (+ a .png
 * sibling). Face gates skip UV-stretched (streak) and long-3D-edge (hole-fill
 * membrane) faces. `face_skip`: optional [nf] mask -- faces with face_skip[f]!=0
 * are omitted ENTIRELY (not painted, not diagnosed); NULL = consider all faces.
 * `lo`/`hi` are the pre-computed contrast window. Out-params (nullable where
 * noted): W,H image dims; fill fraction; multi = #pixels with >1 cover;
 * skip_uv/skip_3d = gated face counts. Returns 0 on success. */
int Rawtex_write_tif(const char *path, CubeTable *ct,
                     const float *verts, const float *uv, size_t nv,
                     const int32_t *faces, size_t nf,
                     const float *normals,
                     const uint8_t *face_skip,
                     double range, int nsteps,
                     double du, double dv, double lo, double hi,
                     double stretch_ratio, double stretch_floor,
                     double max_edge3d,
                     const DiagOpts *diag,
                     size_t *out_W, size_t *out_H,
                     double *out_fill, size_t *out_multi,
                     size_t *out_skip_uv, size_t *out_skip_3d);

/* Percentile contrast window (256-bin histogram) over the sampled population;
 * val[i]/has[i] as returned by sample_vertex, pct_lo/pct_hi e.g. 1/99. The
 * lo/hi it yields feed Rawtex_write_tif so papyrus is contrast-stretched. */
void Rawtex_stretch_window(const double *val, const uint8_t *has, size_t nv,
                           double pct_lo, double pct_hi,
                           double *out_lo, double *out_hi);

#endif /* RAWTEX_BAKE_H */
