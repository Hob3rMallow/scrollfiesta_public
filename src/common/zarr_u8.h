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
    struct ZuSharedChunk *slot_entry[ZARR_U8_CACHE_SLOTS];
    struct ZuSharedCache *shared;
    long    slot_key[ZARR_U8_CACHE_SLOTS][3];
    int     slot_used[ZARR_U8_CACHE_SLOTS];
    size_t  next_evict;
    size_t  misses, hits, missing;
    size_t  loads, shared_hits; /* physical chunk load attempts / reused pinned buffers */
    size_t  slot_limit;      /* 0: all ZARR_U8_CACHE_SLOTS; else at most this many chunks cached (set after open) */
    /* the fast path (2026-10-22): power-of-two chunks by shifts and masks; the last chunk used remembered */
    int     pow2, csh[3];
    long    cmask[3];
    long    last_key[3];
    const uint8_t *last_data;
    int last_valid;          /* a cached missing chunk has valid NULL data */
} ZarrU8;

typedef struct ZarrU8CacheStats {
    size_t handles, chunks, peak_chunks, bytes, peak_bytes;
    size_t loads, shared_hits;
} ZarrU8CacheStats;

/* Parse <root>/0/.zarray.  Returns 0, or -1 when the array is not an
 * uncompressed uint8 zarr v2 with "/" separators. */
int  ZarrU8_open(ZarrU8 *z, const char *root);
/* Each handle is confined to one reading thread. Clones have independent
 * lookup/eviction tables but share immutable loaded chunks with their family.
 * The input volume must remain unchanged while this family is in use.
 * Clone may run concurrently against a live parent; clones outlive the parent.
 * Like open, dst must not already own resources. No global/root-path registry. */
int  ZarrU8_clone(ZarrU8 *dst, const ZarrU8 *src);
void ZarrU8_cache_stats(const ZarrU8 *z, ZarrU8CacheStats *out);
void ZarrU8_close(ZarrU8 *z);

/* Nearest-voxel value at integer (iz, iy, ix); 0 outside / missing. */
uint8_t ZarrU8_get(ZarrU8 *z, long iz, long iy, long ix);
/* The 2 x 2 x 2 block at (iz, iy, ix): v[4a + 2b + c] = ZarrU8_get(z, iz + a, iy + b, ix + c), one chunk lookup where the
 * block lies inside one chunk (2026-10-22: the trilinear reads of the ridge snaps). */
void ZarrU8_get8(ZarrU8 *z, long iz, long iy, long ix, uint8_t v[8]);

/* Nearest-voxel value at a float position (z, y, x). */
uint8_t ZarrU8_sample(ZarrU8 *z, double pz, double py, double px);

/* Present dark chunks are distinct from missing/short chunks. */
int ZarrU8_present(ZarrU8 *z, long iz, long iy, long ix);

int ZarrU8_selftest(void);

#endif
