/* ============================================================================
 * rawtex_layers.h -- a multi-layer CT surface volume on a sheet bake's grid.
 *
 * The ink detectors of the 9 um family (ink_9um, reader-v2) read a surface
 * volume: uint8 [depth, H, W] layers sampled along the surface normal, depth
 * zero at the reference surface, layer index increasing INWARD (toward the
 * scroll axis, the recto side -- docs/ink9um.md, the W043 control).  This
 * module renders that volume for a flattened, UV-carrying mesh on EXACTLY the
 * pixel grid obj_bake_raw uses for sheet_rawtex.tif (Rawtex_plan_field over
 * the mesh UV bbox, pixel centre = origin + (index + 0.5) * step, the same
 * UV-stretch face gate), so an ink map inferred from it overlays the big sheet
 * pixel for pixel.
 *
 * Every pixel centre is located in its covering UV triangle (first strict
 * interior owner wins, a boundary hit is only a fallback), its 3D point and
 * normal interpolated barycentrically, and layer k sampled trilinearly at
 * p + (k - (depth-1)/2) * step * n_inward.  Nothing is resampled in UV: one
 * pixel is one UV unit (native voxels for the assembler's layouts).
 *
 * Output: an uncompressed Zarr v2 array (dtype |u1, chunks [depth,C,C],
 * dimension separator '.', fill 0; chunks with no painted pixel are not
 * written) plus an H x W uint8 mask, 255 where the pixel is painted and all
 * depth layers are fully supported by present RAW (every trilinear corner).
 * Optionally (geom_dir) a float32 Zarr v2 [6,H,W] of the geometry every layer
 * column was sampled along: the depth-zero point and the unit inward normal.
 * ==========================================================================*/
#ifndef RAWTEX_LAYERS_H
#define RAWTEX_LAYERS_H

#include <stddef.h>
#include <stdint.h>

#include "rawtex_bake.h"            /* RawtexUv, RawtexPlan */
#include "../common/raw_sample.h"   /* CubeTable */

typedef struct RawtexLayersOpts {
    int    depth;          /* layers along the normal */
    double step;           /* vox between consecutive layers */
    double offset;         /* vox added to every layer offset (+ = inward) */
    int    tile_cols;      /* raster columns per work tile (multiple of chunk) */
    int    chunk;          /* Zarr chunk edge in y and x */
    double stretch_ratio;  /* obj_bake_raw's UV-stretch face gate */
    double stretch_floor;
    /* Seam fill (vmesh_tifxyz v2's rule): charts of neighbouring cubes abut
     * in UV without sharing vertices, leaving 1-2 px rows of unpainted
     * centres between them -- zero lines inside the model's 128 px context.
     * An unpainted pixel whose nearest painted pixels on opposite sides
     * (left/right or up/down) leave a gap of at most seam_reach pixels and lie
     * at most seam_tol vox apart in 3D takes their distance-weighted
     * interpolation (position and normal); fills never chain within a pass. */
    int    seam_reach;     /* 0 = off */
    double seam_tol;
    int    seam_passes;
    /* Optional geometry volume (NULL = off): an uncompressed float32 Zarr v2
     * [6, H, W] (dtype <f4, chunks [6,C,C], '.' separator, fill NaN; chunks
     * with no painted pixel are not written) holding, per rendered pixel,
     * the depth-zero point (z, y, x) and the unit inward normal (nz, ny, nx)
     * the layers were sampled along -- inside the mask, layer k sits at
     * point + offset_k * normal.  A degenerate interpolated normal is stored
     * as (0, 0, 0): its layers all sample the point.  It carries anything
     * drawn on another surface (e.g. ink labels) onto this grid through 3-D
     * correspondence.  The layer volume is byte-identical with or without it. */
    const char *geom_dir;
    /* Context apron (apron_px 0 = off).  The ink model reads a 128 px context
     * around every pixel; where our sheet stops -- holes, piece borders, box
     * faces -- that context is zeros, and reader-v2 on held-out PHerc0139 loses
     * ~0.11 AUC within 16 px of such an edge (vesuvius-c-f8, 2026-09-29).  With
     * the apron, every unpainted pixel within apron_px of the painted sheet gets
     * CONTEXT-ONLY geometry: a weighted quadric fit (Gaussian weights,
     * sigma = apron_fit_px / 2) of the painted 3-D position field within
     * apron_fit_px of its nearest painted pixel, evaluated at the pixel, so the
     * continuation keeps the sheet's curvature (a tangent step drifts 5-9 vox off
     * a radius-300 turn within 64 px).  Border pixels in one apron_cell_px cell
     * share a fit; fits sample every apron_step_px.  An extrapolation that fails
     * the isometry check (3-D distance to the border pixel between 0.5x and 1.6x
     * the pixel distance, tangent lengths 0.5..2 vox/px, normal within 60 deg of
     * the border pixel's) stays unrendered.  Apron pixels are written to the
     * layer volume only: never to the mask, never to the geometry volume. */
    int    apron_px;
    int    apron_fit_px;
    int    apron_cell_px;
    int    apron_step_px;
    /* Surface regularisation (smooth_sigma_px 0 = off).  Our depth-zero surface
     * follows the per-cube mesh, which wobbles at the few-pixel scale; reader-v2
     * on held-out PHerc0139 loses 0.07-0.17 AUC where our normal is > 8 deg off
     * the team's smooth surface or our depth jitters by >= 0.4 vox over 9 px
     * (away from edges and chart steps it otherwise reads our render better
     * than theirs; output/ink_sf_improve_20260929/geomdiag).  With this on, the
     * painted and seam-filled pixels' points and normals are low-passed with a
     * masked Gaussian (sigma smooth_sigma_px, separable, over painted pixels
     * only); every pixel then moves ONLY along its smoothed normal, by the
     * normal component of (smoothed point - own point) clamped to
     * +-smooth_cap_vox, and takes the smoothed normal as its depth axis.  The
     * in-plane placement of every pixel -- where its ink lands on the sheet --
     * is unchanged.  Runs before the apron, whose fits then see the smoothed
     * sheet. */
    double smooth_sigma_px;
    double smooth_cap_vox;
    /* Step blending (step_sigma_px 0 = off).  Where two charts meet with a 3-D
     * step (adjacent painted pixels > ~1 vox apart; on PHerc0139 mostly a
     * tangential misregistration, p50 2.6 vox) the texture tears; reader-v2
     * loses ~0.07 AUC within 32 px of such steps, whatever the depth.  Each
     * pixel's step strength is smoothstep((|dP| - 1) / 2) over its 4
     * neighbours; spread by a Gaussian of step_spread_px (scaled so a straight
     * seam line weighs 1) it becomes a blend weight w in [0, 1]; the pixel's
     * point then moves by w times (masked-Gaussian mean point, sigma
     * step_sigma_px, minus its own point) -- in 3-D, tangentially too -- and its
     * normal to the w-weighted mix with the mean normal.  Continuous in the
     * step size and the distance: the tear becomes a bend of the texture over
     * a few pixels; away from steps nothing moves.  Runs after the
     * regularisation, before the apron. */
    double step_sigma_px;
    double step_spread_px;
    /* Optional handedness census (NULL = off; filled by RawtexLayers_write,
     * release with RawtexHandedness_free).  See RawtexHandedness. */
    struct RawtexHandedness *hand;
    /* Absolute raster window (NULL = the mesh's UV bbox at 1 px per UV unit).
     * Pixel (x, y) samples (window->umin + x + 0.5, window->vmin + y + 0.5), so a
     * caller that sets the window to columns [c0, c1) of a sheet's full grid
     * (umin_full + c0 .. umin_full + c1, all its rows) renders exactly those
     * columns of the full sheet, whatever part of the mesh it passes -- the way
     * to render a sheet wider than the per-axis raster cap. */
    const RawtexWindow *window;
} RawtexLayersOpts;

/* Which way round each piece shows the recto.  Per painted pixel,
 * s = (dP/dcol x dP/drow) . n_inward, the cross product taken in the stored
 * (z, y, x) component order, from the depth-zero points of its painted
 * neighbours (within 2 vox; central where both sides are painted).  The team's
 * VC3D canvases -- the frame their text is read and labelled in -- always have
 * s > 0 (their layer axis is that cross product and reader-v2 needs it
 * inward); s < 0 is the recto mirror-imaged.  On PHerc0139 10^3 the sign
 * matches the mirror sign of our map onto the team's canvases on 99.998% of
 * 2.05 M pixels.  A piece is a run of painted columns (the ribbon gives each
 * piece its own column range); each is one handedness in practice (the layout
 * places some pieces with the scroll axis running up the rows, others down). */
typedef struct RawtexPieceHand {
    size_t col0, col1;          /* columns [col0, col1) */
    size_t row0, row1;          /* painted rows [row0, row1) */
    size_t px_team;             /* pixels with s > 0: the team's reading handedness */
    size_t px_mirror;           /* pixels with s < 0: mirror-imaged */
} RawtexPieceHand;

typedef struct RawtexHandedness {
    RawtexPieceHand *pieces;
    size_t npieces;
} RawtexHandedness;

void RawtexHandedness_free(RawtexHandedness *h);

typedef struct RawtexLayersStats {
    size_t W, H;               /* raster = the bake's grid */
    double umin, vmin;         /* grid origin (UV of pixel (0,0)'s lower edge) */
    size_t painted_px;         /* pixels with a covering face */
    size_t filled_px;          /* seam-fill pixels (not counted in painted_px) */
    size_t complete_px;        /* painted or filled, and every layer fully supported */
    size_t multi_px;           /* pixels with more than one strict-interior cover */
    size_t skip_uv_faces;      /* faces dropped by the UV-stretch gate */
    size_t apron_px;           /* context-only apron pixels rendered (not in the mask) */
    size_t apron_rejected;     /* apron pixels whose extrapolation failed the checks */
    size_t apron_fits;         /* quadric fits computed */
    size_t smooth_px;          /* pixels regularised (painted + seam-filled, core columns) */
    size_t smooth_capped_px;   /* of those, shifts clamped at smooth_cap_vox */
    double smooth_shift_rms;   /* RMS normal shift (vox) */
    double smooth_turn_mean;   /* mean angle between the interpolated and the smoothed normal (deg) */
    size_t step_px;            /* pixels moved by the step blending (weight > 0.01, core columns) */
    double step_move_rms;      /* RMS move of those pixels (vox) */
    size_t chunks_written;
    double seconds;
} RawtexLayersStats;

/* Orientation of per-vertex normals toward the scroll axis.
 *
 * inward[v*3..] is the unit direction from vertex v toward its axis foot
 * (perpendicular to the axis).  Components are the connected components of
 * the face graph (faces sharing a vertex; the assembler's charts share no
 * vertices, so these are charts or joined chart groups).  Every face votes
 * with dot(cross(e1,e2), mean inward) -- its area times the cosine to the
 * inward direction -- and a component whose total is negative has all of its
 * normals negated, so the normals of one sheet side stay one side.  A
 * per-pixel sign would split a chart wherever it runs tangent to the radius.
 * Normals are then smoothed `smooth_iters` times over the face one-ring.
 * `normals` must be the winding normals of (verts, faces) as MeshNormals_compute
 * returns them: the decision reads the face winding, not the array. */
typedef struct RawtexOrientStats {
    size_t components;         /* face-graph components with >= 1 face */
    size_t flipped;            /* components whose normals were negated */
    size_t weak;               /* components with |sum| < 0.25 * sum|vote| */
    double agree_area;         /* fraction of |vote| agreeing with its component */
    double flipped_area;       /* fraction of |vote| in flipped components */
} RawtexOrientStats;

int RawtexLayers_orient(float *normals, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf, const float *inward,
                        int smooth_iters, RawtexOrientStats *st);

/* Render.  normals must already be oriented inward (RawtexLayers_orient).
 * face_skip (nullable, [nf]) omits faces entirely, as in the bake.  `ct`
 * must be prewarmed (cubetable_prewarm_all) so reads are thread-safe.
 * Writes <zarr_dir>/.zgroup, <zarr_dir>/.zattrs (attrs_json verbatim, or {}),
 * <zarr_dir>/0/.zarray + chunks, and mask_path (uint8 TIFF); with
 * o->geom_dir also <geom_dir>/{.zgroup,.zattrs,0/...}.  Refuses an existing
 * zarr_dir or geom_dir.  Returns 0 on success. */
int RawtexLayers_write(const char *zarr_dir, const char *mask_path,
                       CubeTable *ct, const float *verts, RawtexUv uv,
                       size_t nv, const int32_t *faces, size_t nf,
                       const float *normals, const uint8_t *face_skip,
                       const RawtexLayersOpts *o, const char *attrs_json,
                       RawtexLayersStats *st);

/* The apron's weighted quadric fit of a W x H tile's painted position field
 * (state 1..3) around pixel (cx, cy), exposed for the selftest: coef[k*6 + j]
 * for coordinate k (z, y, x) and basis j (1, u, v, u^2, uv, v^2), with
 * u = (x - cx) / fit_px, v = (y - cy) / fit_px.  Returns 0 on success. */
int RawtexLayers_quadric_fit(const uint8_t *state, const float *pos, long W, long H, long cx, long cy,
                             int fit_px, int step_px, double coef[18]);

/* The surface regularisation on a W x H tile (state 1..3 = on the sheet), in
 * place, exposed for the selftest.  Returns 0 on success (-1: no memory);
 * out-params (nullable) accumulate the core-column statistics over columns
 * [c0, c1): pixels, clamped pixels, sum of squared shifts, sum of turn angles
 * (deg). */
int RawtexLayers_smooth(const uint8_t *state, float *pos, float *nrm, size_t W, size_t H,
                        double sigma_px, double cap_vox, size_t c0, size_t c1,
                        size_t *npx, size_t *ncapped, double *shift2, double *turn);

/* The step blending on a W x H tile (state 1..3 = on the sheet), in place,
 * exposed for the selftest.  Out-params (nullable) accumulate over core
 * columns [c0, c1): pixels with a weight > 0.01, and the sum of their squared
 * moves (vox^2).  Returns 0 on success (-1: no memory). */
int RawtexLayers_stepblend(const uint8_t *state, float *pos, float *nrm, size_t W, size_t H,
                           double sigma_px, double spread_px, size_t c0, size_t c1,
                           size_t *nmoved, double *move2);

int RawtexLayers_selftest(void);

#endif /* RAWTEX_LAYERS_H */
