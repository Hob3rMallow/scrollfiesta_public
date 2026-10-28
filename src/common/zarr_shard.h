#ifndef ZARR_SHARD_INCLUDED
#define ZARR_SHARD_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#else
#include "arena.h"
#endif

/* Uncompressed Zarr v3 sharding_indexed storage. Each shard contains up to
 * 64 adjacent, fixed-size uint8 chunks. The index is little-endian uint64
 * offset/length pairs, followed by CRC32C, at the end of the file.
 * https://zarr-specs.readthedocs.io/en/latest/v3/codecs/sharding-indexed/
 * No full-volume allocation: the writer borrows the caller's chunk band.
 * Readers own one FILE and a bounded index; one reader per calling thread. */
#define ZARR_SHARD_MAX_CHUNKS 64

typedef struct {
    FILE *file;
    size_t count, chunk_bytes;
    uint64_t offsets[ZARR_SHARD_MAX_CHUNKS];
    uint64_t lengths[ZARR_SHARD_MAX_CHUNKS];
} ZarrShard;

typedef struct {
    size_t depth, height, width, tile, rows, cols, count;
} ZarrShardArray;

/* Incremental writer for a shard whose chunks arrive in different row bands.
 * Only the index is retained; no open handle or volume-sized buffer. Calls
 * for one builder must be serial and always use the same fresh path. */
typedef struct {
    size_t count, chunk_bytes;
    uint64_t bytes, seen;
    uint8_t index[16 * ZARR_SHARD_MAX_CHUNKS];
    int closed, failed;
} ZarrShardBuilder;

int ZarrShard_begin(ZarrShardBuilder *builder, size_t chunk_bytes, size_t count);
/* Append one present chunk, including an explicitly present all-zero chunk.
 * Omitted indices are fill; writing an index twice is an error. */
int ZarrShard_append(ZarrShardBuilder *builder, const char *path, size_t index, const uint8_t *chunk);
int ZarrShard_finish(ZarrShardBuilder *builder, const char *path);

/* Write a fresh shard. supplied <= count; the omitted tail and all-zero
 * chunks are fill. An entirely empty shard has no file. Never replace a
 * preexisting file. Returns 0 or -1; written counts stored inner chunks. */
int ZarrShard_write(const char *path, const uint8_t *chunks, size_t chunk_bytes,
                    size_t supplied, size_t count, size_t *written);

/* 1: opened and validated, 0: missing shard, -1: I/O or invalid index.
 * A missing shard is fill; permission errors and corrupt files are errors.
 * Open requires an unowned handle. Close is valid after every open result. */
int ZarrShard_open(ZarrShard *shard, const char *path, size_t chunk_bytes, size_t count);
/* Return 0 or -1. A missing inner chunk writes exactly chunk_bytes zeroes. */
int ZarrShard_read(ZarrShard *shard, size_t index, uint8_t *out);
void ZarrShard_close(ZarrShard *shard);

/* Array metadata for [depth, height, width], inner chunks [depth,tile,tile],
 * outer shards [depth,tile,tile*count]. Directory already exists and is fresh.
 * Shard keys: c/0/<tile row>/<shard column>. */
int ZarrShard_array_meta(const char *directory, size_t depth, size_t height,
                         size_t width, size_t tile, size_t count);
/* Rectangular shards, rows * cols <= 64, with C-order inner chunk indices. */
int ZarrShard_array_meta_grid(const char *directory, size_t depth, size_t height,
                              size_t width, size_t tile, size_t rows, size_t cols);
/* Validate this uncompressed rectangular-shard profile. Unsupported
 * codecs, grids or fill values fail; never reinterpret them as raw bytes. */
#ifndef __cplusplus
int ZarrShard_array_open(Arena_T arena, const char *directory, ZarrShardArray *array);
#endif
/* Small temporary arena owned within C; usable by C++ readers too. */
int ZarrShard_array_load(const char *directory, ZarrShardArray *array);
/* Place a normalized block at any integer column in an existing vertical-
 * shard array. Reads at most two source chunks and one destination chunk at
 * once per worker. Replaces destination files atomically, never writes through
 * a hard link. Caller owns the unpublished destination and serializes blocks. */
int ZarrShard_merge_columns(const char *source, const char *target, size_t col0, size_t *files_written);
int ZarrShard_selftest(void);

#ifdef __cplusplus
}
#endif

#endif
