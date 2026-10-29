#ifndef MARBLE_DIAG_INCLUDED
#define MARBLE_DIAG_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "rawtex_bake.h"

/* RAW/mesh alignment audit for a UV-carrying scroll surface.  This is not an
 * image-texture heuristic: sparse UV samples are mapped through the mesh into
 * source XYZ, a short profile is searched along the mesh normal, and the RAW
 * tangent tensor is compared with the mesh's intended U/V directions. */
typedef struct MarbleDiagOpts {
    int    stride_pixels;       /* coarse audit-cell edge in baked pixels (16) */
    double depth_range;         /* +/- normal search reach in voxels (4) */
    double depth_step;          /* normal-profile spacing in voxels (1) */
    double tensor_radius;       /* tangent-patch radius in voxels (2) */
    double score_threshold;     /* retained-region threshold [0,1] (.40) */
    int    minimum_region_cells;/* discard smaller 4-connected islands (12) */
} MarbleDiagOpts;

typedef struct MarbleDiagStats {
    size_t raster_width, raster_height;
    size_t audit_width, audit_height;
    size_t sampled_cells;
    size_t highlighted_cells;
    size_t highlighted_pixels;
    size_t highlighted_regions;
    double score_p50, score_p90, score_p95, score_p99;
} MarbleDiagStats;

void MarbleDiag_defaults(MarbleDiagOpts *opts);

/* `prefix` names all outputs; `rawtex_path` is the matching full-resolution
 * grayscale bake. Candidate normal-depth evidence is aggregated in UV before
 * a winning depth is selected, so isolated CT texture cannot masquerade as a
 * coherent sheet displacement. Outputs are:
 *   <prefix>_diagmarble_{offset,support_center,support_best,gain,axis,
 *                        axis_broad,disorder,surface_miss,score,mask}.png
 *   <prefix>_marble_mask.png
 *   <prefix>_marble_mask.pgm       (same binary mask, lossless solver input)
 *   <prefix>_marble_overlay.png       (RAW sheet, green tint on detections)
 *   <prefix>_marble_report.json
 * Returns 0 on success. */
int MarbleDiag_write(const char *prefix, const char *rawtex_path,
                     CubeTable *ct,
                     const float *verts, const float *uv, size_t nv,
                     const int32_t *faces, size_t nf,
                     const float *vertex_normals,
                     double raster_du, double raster_dv, size_t raster_max_px,
                     const RawtexWindow *window,
                     const MarbleDiagOpts *opts, MarbleDiagStats *stats);

int MarbleDiag_write_field(const char *prefix, const char *rawtex_path,
                     CubeTable *ct,
                     const float *verts, RawtexUv uv, size_t nv,
                     const int32_t *faces, size_t nf,
                     const float *vertex_normals,
                     double raster_du, double raster_dv, size_t raster_max_px,
                     const RawtexWindow *window,
                     const MarbleDiagOpts *opts, MarbleDiagStats *stats);

#endif
