#ifndef ZARR_U8_INCLUDED
#define ZARR_U8_INCLUDED

#include <stddef.h>
#include <stdint.h>

/* ============================================================================
 * zarr_u8.h -- read-only access to an uncompressed uint8 zarr v2 volume
 * (the RAW CT: `<root>/0/.zarray` with dtype |u1, compressor null,
 * dimension_separator "/", C order; chunks at `<root>/0/<cz>/<cy>/<cx>`).
 * A small chunk cache keeps the working set of a lane stage in RAM; every
 * read outside the volume or of a missing chunk is the fill value (0).
 * ==========================================================================*/

#define ZARR_U8_CACHE_SLOTS 512   /* 128^3 chunks -> 1 GB at most */

typedef struct {
    char    root[1024];
    long    shape[3];        /* z, y, x */
    long    chunk[3];
    uint8_t *slot_data[ZARR_U8_CACHE_SLOTS];
    long    slot_key[ZARR_U8_CACHE_SLOTS][3];
    int     slot_used[ZARR_U8_CACHE_SLOTS];
    size_t  next_evict;
    size_t  misses, hits, missing;
} ZarrU8;

/* Parse <root>/0/.zarray.  Returns 0, or -1 when the array is not an
 * uncompressed uint8 zarr v2 with "/" separators. */
int  ZarrU8_open(ZarrU8 *z, const char *root);
void ZarrU8_close(ZarrU8 *z);

/* Nearest-voxel value at integer (iz, iy, ix); 0 outside / missing. */
uint8_t ZarrU8_get(ZarrU8 *z, long iz, long iy, long ix);

/* Nearest-voxel value at a float position (z, y, x). */
uint8_t ZarrU8_sample(ZarrU8 *z, double pz, double py, double px);

int ZarrU8_selftest(void);

#endif
