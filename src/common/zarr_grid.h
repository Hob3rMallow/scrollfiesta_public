#ifndef ZARR_GRID_INCLUDED
#define ZARR_GRID_INCLUDED

/* Copy an in-bounds, 128-aligned region from a local Zarr v2 array directory
 * to the pipeline's world-indexed 128^3 TIFF cubes. Supports C-order uint8,
 * slash-separated chunks, compressor null or Blosc. Missing/corrupt chunks
 * fail; no absent data is silently replaced by zero. Binary maps nonzero
 * source values to 255. Every TIFF is decoded and compared before commit.
 * Existing output cubes are never overwritten. */
int ZarrGrid_export(const char *array_dir, const char *out_dir,
                    const long bbox[6], int binary);
int ZarrGrid_selftest(void);

#endif
