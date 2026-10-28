/* ============================================================================
 * rawtex_tifxyz.h -- a flattened sheet as a VC3D tifxyz segment for villa's
 * vc_render_tifxyz (the team's official surface-volume renderer).
 *
 * Renderer contract (villa 6e53201ac; measured on a synthetic plane by
 * scripts/ink/vcrender_contract_check.py and read in QuadSurface::gen): with
 * --scale 1 --group-idx 0 a scale-1 tifxyz renders on a canvas of the grid's
 * own size, and canvas pixel (r, c)
 *   - samples the grid BILINEARLY at (r + 0.5, c + 0.5): nodes (r..r+1, c..c+1);
 *   - is valid when node (r + 1, c + 1) is (nearest, std::lround);
 *   - takes that node's normal N = normalize((P[r+1][c+2] - P[r+1][c]) x
 *     (P[r+2][c+1] - P[r][c+1])), (x, y, z) order, NaN unless all four
 *     neighbours are valid (and on the grid's border rows / columns);
 *   - renders slice k at point + (k - (n-1)/2) * step * N (--flip-normals: -N),
 *     trilinear, truncated to integer grey.
 * A pixel whose node (r+1, c+1) and its four neighbours are valid but whose
 * cell corner (r, c) is not is still rendered -- from a point bilinearly mixed
 * with the -1 sentinel.  RawtexTifxyz_canvas_mask leaves such pixels out.
 *
 * The sheet is therefore exported on its CORNER grid: node (r, c) holds the
 * sheet's point at the corner between bake pixels (r-1..r, c-1..c), i.e. the
 * review bake's grid shifted by half a pixel and one node wider and taller, so
 * canvas pixel (r, c) is bake pixel (r, c) (RawtexLayersOpts.window + pos_sink
 * produce it).
 *
 * Handedness: with --flip-normals the layer index grows along
 * (dP/dcol x dP/drow) in (z, y, x) order, which is inward only where the grid
 * has the team's reading handedness (RawtexHandedness s > 0).  Each
 * mirror-imaged piece (a run of painted node columns whose census is mostly
 * s < 0) is REFLECTED -- its node columns [col0, col1) reversed -- which makes
 * it team-handed and its text read the right way round; its canvas columns
 * [col0 - 1, col1) then show the sheet's columns in reverse order
 * (RawtexTifxyz_canvas_to_sheet undoes that for anything inferred on the
 * canvas).  Pixels whose own sign disagrees with their piece's majority would
 * still render outward; the canvas mask leaves them out and counts them.
 * ==========================================================================*/
#ifndef RAWTEX_TIFXYZ_H
#define RAWTEX_TIFXYZ_H

#include <stddef.h>
#include <stdint.h>

#include "../common/raw_sample.h"   /* CubeTable */

typedef struct RawtexTifxyzPiece {
    size_t col0, col1;          /* painted node columns [col0, col1) */
    size_t row0, row1;          /* painted node rows [row0, row1) */
    size_t px_team, px_mirror;  /* handedness census before any reflection */
    int    reflected;           /* 1 = node columns reversed: now team-handed */
} RawtexTifxyzPiece;

/* Reflect every piece with px_mirror > px_team: node columns [col0, col1) of the
 * three position planes (pos: plane-major (z, y, x), W * H each) and of sign
 * (nullable, W * H) are reversed, and the piece's signs negated.  Sets
 * .reflected (0 on the others).  Returns the number of reflected pieces. */
size_t RawtexTifxyz_reflect(float *pos, int8_t *sign, size_t W, size_t H,
                            RawtexTifxyzPiece *pieces, size_t npieces);

typedef struct RawtexTifxyzMaskStats {
    size_t rendered_px;     /* canvas pixels vc_render renders (valid node, finite normal) */
    size_t clean_px;        /* ...whose four cell corners are valid */
    size_t polluted_px;     /* ...rendered from a point mixed with an invalid corner */
    size_t reversed_px;     /* clean, but the node's sign says the layers run outward */
    size_t unsupported_px;  /* clean, not reversed, but a layer reaches missing RAW */
    size_t mask_px;         /* in the mask */
} RawtexTifxyzMaskStats;

/* The canvas mask of vc_render_tifxyz --scale 1 --num-slices depth --slice-step
 * step [--flip-normals when flip] on this grid (canvas W x H = the grid): 255
 * where the renderer samples a clean cell along a finite normal, the node's sign
 * (nullable) is not negative, and every layer's trilinear footprint is present
 * RAW (ct nullable = not checked; it must be prewarmed); 0 elsewhere.  pos is
 * plane-major (z, y, x) with NaN for invalid nodes.  Threads: OpenMP over rows.
 * Returns 0 on success, -1 on bad input. */
int RawtexTifxyz_canvas_mask(const float *pos, const int8_t *sign, size_t W, size_t H,
                             CubeTable *ct, int depth, double step, int flip,
                             uint8_t *mask, RawtexTifxyzMaskStats *st);

/* Write the segment directory <dir> (created with parents): x.tif, y.tif, z.tif
 * (float32, uncompressed, invalid = -1: NaN nodes are set to -1 in pos, which is
 * modified), mask.tif (uint8, 255 on valid nodes) and meta.json (type seg,
 * format tifxyz, uuid, scale [1, 1], bbox in (x, y, z), area_vx2 = clean cells)
 * with `extra_json` (nullable: members without braces, e.g. "\"k\": 1") appended.
 * Refuses an existing <dir>/meta.json.  Returns 0 on success. */
int RawtexTifxyz_write(const char *dir, const char *uuid, float *pos, size_t W, size_t H,
                       const char *extra_json);

/* RawtexTifxyz_write with deflate-compressed strips (TiffIO_save_float2d_deflate / TiffIO_save_deflate): the same
 * grids, values and meta.json, a few percent of the bytes on a sparse sheet canvas.  sheet_layers --tifxyz block
 * exports (2026-10-21); readers that memory-map a plane need tifffile's zarr view instead (scripts/tif_view.py). */
int RawtexTifxyz_write_packed(const char *dir, const char *uuid, float *pos, size_t W, size_t H,
                              const char *extra_json);

/* Map a raster on the canvas (cW x cH = the tifxyz grid) onto the bake grid
 * (sW x sH, at most (cW - 1) x (cH - 1)): sheet (r, t) = canvas (r, c) with
 * t = c, or t = c0 + c1 - 1 - c for canvas columns [c0, c1) = [col0 - 1, col1)
 * of a reflected piece; sheet pixels no canvas pixel maps to are set to 0.
 * Returns 0 on success, -1 on bad input. */
int RawtexTifxyz_canvas_to_sheet(const uint8_t *canvas, size_t cW, size_t cH,
                                 const RawtexTifxyzPiece *pieces, size_t npieces,
                                 uint8_t *sheet, size_t sW, size_t sH);

/* Map a layer volume rendered on the canvas -- an UNCOMPRESSED uint8 Zarr v2
 * array (vc_render_tifxyz --zarr-compressor none: shape [D, cH, cW], chunks
 * [D, C, C], either dimension separator) -- onto the bake grid, as
 * RawtexTifxyz_canvas_to_sheet does per layer: an uncompressed Zarr v2 group
 * <out_dir> ('.zgroup', '.zattrs' = attrs_json, array '0' [D, sH, sW], chunks
 * [D, 128, 128], '.' separator, fill 0) in sheet_layers' format.  Only chunks
 * holding a pixel of `sheet_mask` (sW x sH, nullable = every non-zero chunk) are
 * written.  Refuses an existing <out_dir>/.zgroup; publishes through
 * <out_dir>.partial.  Out-params (nullable): chunks written, canvas chunks
 * absent (read as fill 0).  Returns 0 on success. */
int RawtexTifxyz_canvas_to_sheet_zarr(const char *canvas_array_dir, const RawtexTifxyzPiece *pieces,
                                      size_t npieces, const uint8_t *sheet_mask, size_t sW, size_t sH,
                                      const char *out_dir, const char *attrs_json, size_t *chunks_written,
                                      size_t *canvas_chunks_absent);

/* The same mapping, with an exact k -> D-1-k permutation at nonzero pixels of
 * sheet_reverse (sheet coordinates, nullable). Both masks are on the OUTPUT
 * grid: a reflected piece must be mapped before constructing sheet_reverse.
 * Existing inward columns and world/UV coordinates are unchanged. */
int RawtexTifxyz_canvas_to_sheet_zarr_depth(const char *canvas_array_dir, const RawtexTifxyzPiece *pieces,
                                            size_t npieces, const uint8_t *sheet_mask,
                                            const uint8_t *sheet_reverse, size_t sW, size_t sH,
                                            const char *out_dir, const char *attrs_json, size_t *chunks_written,
                                            size_t *canvas_chunks_absent);
/* As above, with optional vertical Zarr v3 shards. shard_rows == 0 keeps v2;
 * 1..64 groups that many 128-row chunks, retaining only one row band and a
 * fixed index per output column. Vertical shards preserve 128-column block
 * alignment, so intake can link them without rewriting their contents. */
int RawtexTifxyz_canvas_to_sheet_zarr_storage(const char *canvas_array_dir, const RawtexTifxyzPiece *pieces,
                                              size_t npieces, const uint8_t *sheet_mask,
                                              const uint8_t *sheet_reverse, size_t sW, size_t sH,
                                              const char *out_dir, const char *attrs_json, size_t shard_rows,
                                              size_t *chunks_written, size_t *canvas_chunks_absent);

/* The review bake's contrast window (obj_bake_raw: the pct_lo / pct_hi
 * percentiles of a 256-bin histogram over the mask's pixels, grey = 255 *
 * clamp((v - lo) / (hi - lo)), rounded), applied in place to an 8-bit raster of
 * n pixels; pixels off the mask become 0.  *lo / *hi (nullable) receive the
 * window.  Returns 0 (-1 on bad input). */
int RawtexTifxyz_stretch(uint8_t *img, const uint8_t *mask, size_t n, double pct_lo, double pct_hi,
                         int *lo, int *hi);

/* The stretch in three parts, so one window can serve a sheet exported in
 * column blocks (sheet_layers --tifxyz --window): accumulate each block's
 * masked 256-bin histogram (RawtexTifxyz_histogram adds into hist), take the
 * window of the summed histogram (RawtexTifxyz_window: the percentile rule of
 * RawtexTifxyz_stretch, 0..255 for an empty histogram; -1 on bad input), and
 * apply it to every block (RawtexTifxyz_apply_window: the same rounding, 0 off
 * the mask; hi <= lo means 0..255).  RawtexTifxyz_stretch == the three on one
 * raster, byte for byte. */
void RawtexTifxyz_histogram(const uint8_t *img, const uint8_t *mask, size_t n, size_t hist[256]);
int RawtexTifxyz_window(const size_t hist[256], double pct_lo, double pct_hi, int *lo, int *hi);
void RawtexTifxyz_apply_window(uint8_t *img, const uint8_t *mask, size_t n, int lo, int hi);

/* Copy a block raster (bH x bW) into columns [col0, col0 + bW) of an H x W
 * sheet raster (rows must match): the columns of a sheet exported in blocks
 * put back on its full grid.  Returns -1 if the block does not fit. */
int RawtexTifxyz_compose(uint8_t *sheet, size_t W, size_t H, const uint8_t *block, size_t bW, size_t bH,
                         size_t col0);

/* A single-page 8-bit gray TIFF, stripped or tiled, any libtiff codec (the ink
 * runtime writes tiled LZW BigTIFFs): *img is malloc'd (caller frees).  0 = ok. */
int RawtexTifxyz_load_gray(const char *path, uint8_t **img, size_t *W, size_t *H);

/* Save an 8-bit gray raster as a 128 x 128 tiled, deflate-compressed TIFF
 * (BigTIFF above 2 GiB of pixels).  0 = ok. */
int RawtexTifxyz_save_gray(const char *path, const uint8_t *img, size_t W, size_t H);

/* In-process unit tests (reflection round trip, the canvas mask's vc_render
 * rules on a synthetic grid, the canvas-to-sheet map, the writer, the gray
 * TIFF round trip).  Returns 0 if all pass. */
int RawtexTifxyz_selftest(void);

#endif /* RAWTEX_TIFXYZ_H */
