#include "zarr_shard.h"
#include "ves_platform.h"
#include "json_read.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

static void zs_u64(uint8_t *out, uint64_t value)
{
    for (unsigned int i = 0; i < 8; i++) out[i] = (uint8_t)(value >> (8 * i));
}

static uint64_t zs_read_u64(const uint8_t *in)
{
    uint64_t value = 0;
    for (unsigned int i = 0; i < 8; i++) value |= (uint64_t)in[i] << (8 * i);
    return value;
}

static uint32_t zs_crc(const uint8_t *data, size_t size)
{
    uint32_t crc = UINT32_MAX;
    for (size_t i = 0; i < size; i++) {
        crc ^= data[i];
        for (unsigned int b = 0; b < 8; b++)
            crc = (crc >> 1) ^ (0x82f63b78u & (0u - (crc & 1u)));
    }
    return ~crc;
}

static FILE *zs_fresh(const char *path)
{
    FILE *file = NULL;
#ifdef _WIN32
    int fd = _open(path, _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY, _S_IREAD | _S_IWRITE);
    if (fd >= 0) { file = _fdopen(fd, "wb"); if (file == NULL) _close(fd); }
#else
    int fd = open(path, O_CREAT | O_EXCL | O_WRONLY, 0666);
    if (fd >= 0) { file = fdopen(fd, "wb"); if (file == NULL) close(fd); }
#endif
    return file;
}

static int zs_seek(FILE *file, uint64_t offset)
{
    if (offset > INT64_MAX) return -1;
#ifdef _WIN32
    return _fseeki64(file, (int64_t)offset, SEEK_SET);
#else
    return fseeko(file, (off_t)offset, SEEK_SET);
#endif
}

static int zs_size(FILE *file, uint64_t *size)
{
#ifdef _WIN32
    if (_fseeki64(file, 0, SEEK_END)) return -1;
    int64_t end = _ftelli64(file);
#else
    if (fseeko(file, 0, SEEK_END)) return -1;
    off_t end = ftello(file);
#endif
    if (end < 0) return -1;
    *size = (uint64_t)end;
    return 0;
}

int ZarrShard_begin(ZarrShardBuilder *builder, size_t chunk_bytes, size_t count)
{
    if (builder == NULL) return -1;
    memset(builder, 0, sizeof *builder);
    builder->failed = 1;
    if (!chunk_bytes || !count || count > ZARR_SHARD_MAX_CHUNKS || chunk_bytes > (size_t)INT64_MAX / count)
        return -1;
    builder->count = count;
    builder->chunk_bytes = chunk_bytes;
    memset(builder->index, 255, sizeof builder->index);
    builder->failed = 0;
    return 0;
}

static FILE *zs_resume(const char *path, uint64_t expected)
{
    if (!expected) return zs_fresh(path);
    FILE *file = fopen(path, "r+b");
    uint64_t actual = 0;
    if (file && (zs_size(file, &actual) || actual != expected)) { fclose(file); return NULL; }
    return file;                         /* zs_size leaves the cursor at the end */
}

int ZarrShard_append(ZarrShardBuilder *builder, const char *path, size_t index, const uint8_t *chunk)
{
    if (!builder || !path || !chunk || builder->failed || builder->closed || !builder->count ||
        index >= builder->count || (builder->seen & ((uint64_t)1 << index))) return -1;
    FILE *file = zs_resume(path, builder->bytes);
    if (!file) { builder->failed = 1; return -1; }
    int failed = fwrite(chunk, 1, builder->chunk_bytes, file) != builder->chunk_bytes;
    if (fclose(file)) failed = 1;
    if (failed) { builder->failed = 1; return -1; }
    zs_u64(builder->index + 16 * index, builder->bytes);
    zs_u64(builder->index + 16 * index + 8, builder->chunk_bytes);
    builder->bytes += builder->chunk_bytes;
    builder->seen |= (uint64_t)1 << index;
    return 0;
}

int ZarrShard_finish(ZarrShardBuilder *builder, const char *path)
{
    if (!builder || !path || builder->failed || builder->closed || !builder->count) return -1;
    builder->closed = 1;
    if (!builder->bytes) {
        struct stat st;
        if (stat(path, &st) == 0 || errno != ENOENT) { builder->failed = 1; return -1; }
        return 0;
    }
    FILE *file = zs_resume(path, builder->bytes);
    if (!file) { builder->failed = 1; return -1; }
    uint8_t tail[4];
    uint32_t crc = zs_crc(builder->index, 16 * builder->count);
    for (unsigned int i = 0; i < 4; i++) tail[i] = (uint8_t)(crc >> (8 * i));
    int failed = fwrite(builder->index, 1, 16 * builder->count, file) != 16 * builder->count ||
                 fwrite(tail, 1, sizeof tail, file) != sizeof tail;
    if (fclose(file)) failed = 1;
    builder->failed = failed;
    return failed ? -1 : 0;
}

int ZarrShard_write(const char *path, const uint8_t *chunks, size_t chunk_bytes,
                    size_t supplied, size_t count, size_t *written)
{
    uint8_t index[16 * ZARR_SHARD_MAX_CHUNKS + 4];
    uint8_t present[ZARR_SHARD_MAX_CHUNKS] = {0};
    size_t stored = 0;
    uint64_t offset = 0;
    FILE *file = NULL;
    if (written != NULL) *written = 0;
    if (path == NULL || count == 0 || count > ZARR_SHARD_MAX_CHUNKS || supplied > count ||
        chunk_bytes == 0 || supplied > SIZE_MAX / chunk_bytes || (supplied && chunks == NULL)) return -1;
    memset(index, 255, 16 * count);
    for (size_t c = 0; c < supplied; c++) {
        const uint8_t *chunk = chunks + c * chunk_bytes;
        size_t i = 0;
        while (i < chunk_bytes && chunk[i] == 0) i++;
        present[c] = (uint8_t)(i != chunk_bytes);
        stored += present[c];
    }
    /* Even an all-fill write may not silently reuse an earlier output. */
    if (stored == 0) {
        struct stat st;
        if (stat(path, &st) == 0 || errno != ENOENT) return -1;
        return 0;
    }
    file = zs_fresh(path);
    if (file == NULL) return -1;
    int failed = 0;
    for (size_t c = 0; c < supplied; c++) {
        if (!present[c]) continue;
        if (offset > UINT64_MAX - chunk_bytes || fwrite(chunks + c * chunk_bytes, 1, chunk_bytes, file) != chunk_bytes) {
            failed = 1;
            break;
        }
        zs_u64(index + 16 * c, offset);
        zs_u64(index + 16 * c + 8, chunk_bytes);
        offset += chunk_bytes;
    }
    uint32_t crc = zs_crc(index, 16 * count);
    for (unsigned int i = 0; i < 4; i++) index[16 * count + i] = (uint8_t)(crc >> (8 * i));
    if (!failed && fwrite(index, 1, 16 * count + 4, file) != 16 * count + 4) failed = 1;
    if (fclose(file)) failed = 1;
    if (!failed && written != NULL) *written = stored;
    return failed ? -1 : 0;
}

int ZarrShard_open(ZarrShard *shard, const char *path, size_t chunk_bytes, size_t count)
{
    uint8_t index[16 * ZARR_SHARD_MAX_CHUNKS + 4];
    uint64_t size = 0;
    if (shard == NULL) return -1;
    memset(shard, 0, sizeof *shard);
    if (path == NULL || count == 0 || count > ZARR_SHARD_MAX_CHUNKS || chunk_bytes == 0) return -1;
    shard->count = count;
    shard->chunk_bytes = chunk_bytes;
    shard->file = fopen(path, "rb");
    if (shard->file == NULL) return errno == ENOENT ? 0 : -1;
    size_t index_size = 16 * count + 4;
    if (zs_size(shard->file, &size) || size < index_size || zs_seek(shard->file, size - index_size) ||
        fread(index, 1, index_size, shard->file) != index_size) goto fail;
    uint32_t actual = zs_crc(index, 16 * count), expected = 0;
    for (unsigned int i = 0; i < 4; i++) expected |= (uint32_t)index[16 * count + i] << (8 * i);
    if (actual != expected) goto fail;
    for (size_t c = 0; c < count; c++) {
        uint64_t offset = zs_read_u64(index + 16 * c), bytes = zs_read_u64(index + 16 * c + 8);
        if (offset == UINT64_MAX && bytes == UINT64_MAX) { shard->lengths[c] = 0; continue; }
        if (bytes != chunk_bytes || offset > size - index_size || bytes > size - index_size - offset) goto fail;
        shard->offsets[c] = offset;
        shard->lengths[c] = bytes;
    }
    return 1;
fail:
    ZarrShard_close(shard);
    return -1;
}

int ZarrShard_read(ZarrShard *shard, size_t index, uint8_t *out)
{
    if (shard == NULL || out == NULL || index >= shard->count || shard->chunk_bytes == 0) return -1;
    if (shard->lengths[index] == 0) { memset(out, 0, shard->chunk_bytes); return 0; }
    if (shard->file == NULL || zs_seek(shard->file, shard->offsets[index]) ||
        fread(out, 1, shard->chunk_bytes, shard->file) != shard->chunk_bytes) return -1;
    return 0;
}

void ZarrShard_close(ZarrShard *shard)
{
    if (shard == NULL) return;
    if (shard->file != NULL) fclose(shard->file);
    memset(shard, 0, sizeof *shard);
}

int ZarrShard_array_meta(const char *directory, size_t depth, size_t height,
                         size_t width, size_t tile, size_t count)
{
    return ZarrShard_array_meta_grid(directory, depth, height, width, tile, 1, count);
}

int ZarrShard_array_meta_grid(const char *directory, size_t depth, size_t height,
                              size_t width, size_t tile, size_t rows, size_t cols)
{
    char path[1600];
    if (directory == NULL || !depth || !height || !width || !tile || !rows || !cols ||
        rows > ZARR_SHARD_MAX_CHUNKS / cols || tile > SIZE_MAX / rows || tile > SIZE_MAX / cols) return -1;
    int n = snprintf(path, sizeof path, "%s/zarr.json", directory);
    if (n < 0 || (size_t)n >= sizeof path) return -1;
    FILE *file = zs_fresh(path);
    if (file == NULL) return -1;
    int ok = fprintf(file,
        "{\n \"zarr_format\": 3, \"node_type\": \"array\", \"data_type\": \"uint8\",\n"
        " \"shape\": [%zu, %zu, %zu],\n"
        " \"chunk_grid\": {\"name\": \"regular\", \"configuration\": {\"chunk_shape\": [%zu, %zu, %zu]}},\n"
        " \"chunk_key_encoding\": {\"name\": \"default\", \"configuration\": {\"separator\": \"/\"}},\n"
        " \"fill_value\": 0, \"attributes\": {},\n"
        " \"codecs\": [{\"name\": \"sharding_indexed\", \"configuration\": {\n"
        "  \"chunk_shape\": [%zu, %zu, %zu],\n"
        "  \"codecs\": [{\"name\": \"bytes\"}],\n"
        "  \"index_codecs\": [{\"name\": \"bytes\", \"configuration\": {\"endian\": \"little\"}}, "
        "{\"name\": \"crc32c\"}], \"index_location\": \"end\"}}]\n}\n",
        depth, height, width, depth, tile * rows, tile * cols, depth, tile, tile) >= 0;
    if (fclose(file)) ok = 0;
    return ok ? 0 : -1;
}

static int zs_name(const JsonValue *object, const char *key, const char *expected)
{
    const char *actual = Json_as_string(Json_object_get(object, key));
    return actual != NULL && strcmp(actual, expected) == 0;
}

static int zs_dimensions(const JsonValue *json, size_t shape[3])
{
    if (Json_array_len(json) != 3) return -1;
    for (size_t d = 0; d < 3; d++) {
        double n = Json_as_double(Json_array_get(json, d), 0);
        if (!isfinite(n) || n < 1 || n >= (double)SIZE_MAX || n != floor(n)) return -1;
        shape[d] = (size_t)n;
    }
    return 0;
}

int ZarrShard_array_open(Arena_T arena, const char *directory, ZarrShardArray *array)
{
    char path[1600];
    size_t shape[3], outer[3], inner[3];
    if (arena == NULL || directory == NULL || array == NULL) return -1;
    memset(array, 0, sizeof *array);
    int n = snprintf(path, sizeof path, "%s/zarr.json", directory);
    if (n < 0 || (size_t)n >= sizeof path) return -1;
    const JsonValue *meta = Json_parse_file(arena, path, NULL);
    const JsonValue *grid = Json_object_get(meta, "chunk_grid"), *key = Json_object_get(meta, "chunk_key_encoding");
    const JsonValue *codecs = Json_object_get(meta, "codecs");
    const JsonValue *sharding = Json_array_get(codecs, 0);
    const JsonValue *config = Json_object_get(sharding, "configuration");
    const JsonValue *body = Json_object_get(config, "codecs"), *index = Json_object_get(config, "index_codecs");
    const JsonValue *data_bytes = Json_array_get(body, 0), *index_bytes = Json_array_get(index, 0);
    const JsonValue *key_config = Json_object_get(key, "configuration");
    const JsonValue *location_value = Json_object_get(config, "index_location");
    const JsonValue *separator_value = Json_object_get(key_config, "separator");
    const char *location = Json_as_string(location_value);
    const char *separator = Json_as_string(separator_value);
    const JsonValue *storage = Json_object_get(meta, "storage_transformers");
    if (meta == NULL || Json_member_double(meta, "zarr_format", 0) != 3 || !zs_name(meta, "node_type", "array") ||
        !zs_name(meta, "data_type", "uint8") || Json_type(Json_object_get(meta, "fill_value")) != JSON_NUMBER ||
        Json_member_double(meta, "fill_value", -1) != 0 || !zs_name(grid, "name", "regular") ||
        !zs_name(key, "name", "default") || (separator_value && (!separator || strcmp(separator, "/"))) ||
        (storage && (Json_type(storage) != JSON_ARRAY || Json_array_len(storage))) ||
        Json_array_len(codecs) != 1 || !zs_name(sharding, "name", "sharding_indexed") ||
        (location_value && (!location || strcmp(location, "end"))) || Json_array_len(body) != 1 || !zs_name(data_bytes, "name", "bytes") ||
        Json_array_len(index) != 2 || !zs_name(index_bytes, "name", "bytes") ||
        !zs_name(Json_object_get(index_bytes, "configuration"), "endian", "little") ||
        !zs_name(Json_array_get(index, 1), "name", "crc32c") ||
        zs_dimensions(Json_object_get(meta, "shape"), shape) ||
        zs_dimensions(Json_object_get(Json_object_get(grid, "configuration"), "chunk_shape"), outer) ||
        zs_dimensions(Json_object_get(config, "chunk_shape"), inner)) return -1;
    if (inner[0] != shape[0] || inner[1] != inner[2] || outer[0] != inner[0] || outer[1] % inner[1] ||
        outer[2] % inner[2] || outer[1] < inner[1] || outer[2] < inner[2] ||
        outer[1] / inner[1] > ZARR_SHARD_MAX_CHUNKS / (outer[2] / inner[2]) ||
        inner[0] > SIZE_MAX / inner[1] || inner[0] * inner[1] > SIZE_MAX / inner[2]) return -1;
    array->depth = shape[0]; array->height = shape[1]; array->width = shape[2];
    array->tile = inner[1]; array->rows = outer[1] / inner[1]; array->cols = outer[2] / inner[2];
    array->count = array->rows * array->cols;
    return 0;
}

int ZarrShard_array_load(const char *directory, ZarrShardArray *array)
{
    Arena_T arena = Arena_new_sized(16384);
    int rc = ZarrShard_array_open(arena, directory, array);
    Arena_dispose(&arena);
    return rc;
}

int ZarrShard_merge_columns(const char *source, const char *target, size_t col0, size_t *files_written)
{
    ZarrShardArray src, dst;
    if (files_written) *files_written = 0;
    if (!source || !target || !strcmp(source, target) || ZarrShard_array_load(source, &src) ||
        ZarrShard_array_load(target, &dst) || src.cols != 1 || dst.cols != 1 || src.rows != dst.rows ||
        src.tile != dst.tile || src.depth != dst.depth || src.height > dst.height || col0 > dst.width ||
        src.width > dst.width - col0) return -1;
    const size_t C = src.tile, R = src.rows, bytes = src.depth * C * C;
    if (src.height > SIZE_MAX - C * R || dst.width > SIZE_MAX - C) return -1;
    const size_t first = col0 / C, last = (col0 + src.width - 1) / C, nx = last - first + 1;
    const size_t ny = (src.height + C * R - 1) / (C * R);
    if (ny > (size_t)INT64_MAX / nx) return -1;
    const int64_t groups = (int64_t)(nx * ny);
    size_t written = 0;
    int failed = 0;
    int64_t group = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(+:written) reduction(|:failed)
#endif
    for (group = 0; group < groups; group++) {
        const size_t sy = (size_t)group / nx, gx = first + (size_t)group % nx;
        const size_t g0 = gx * C > col0 ? gx * C : col0;
        const size_t g1 = (gx + 1) * C < col0 + src.width ? (gx + 1) * C : col0 + src.width;
        const size_t bx0 = (g0 - col0) / C, bx1 = (g1 - col0 - 1) / C;
        ZarrShard old = {0}, inputs[2] = {{0}, {0}};
        ZarrShardBuilder builder;
        uint8_t *out = (uint8_t *)malloc(bytes), *in = (uint8_t *)malloc(bytes);
        char path[2400], temporary[2440], input[2400];
        int bad = !out || !in;
        snprintf(path, sizeof path, "%s/c/0/%zu/%zu", target, sy, gx);
        snprintf(temporary, sizeof temporary, "%s.merge.partial", path);
        if (!bad && (ves_ensure_parent_dir(path) || ZarrShard_open(&old, path, bytes, R) < 0 ||
                     ZarrShard_begin(&builder, bytes, R))) bad = 1;
        for (size_t bx = bx0; !bad && bx <= bx1; bx++) {
            snprintf(input, sizeof input, "%s/c/0/%zu/%zu", source, sy, bx);
            if (ZarrShard_open(&inputs[bx - bx0], input, bytes, R) < 0) bad = 1;
        }
        for (size_t r = 0; !bad && r < R; r++) {
            if (ZarrShard_read(&old, r, out)) { bad = 1; break; }
            for (size_t bx = bx0; !bad && bx <= bx1; bx++) {
                if (ZarrShard_read(&inputs[bx - bx0], r, in)) { bad = 1; break; }
                const size_t a = g0 > col0 + bx * C ? g0 : col0 + bx * C;
                const size_t b = g1 < col0 + (bx + 1) * C ? g1 : col0 + (bx + 1) * C;
                for (size_t d = 0; d < src.depth; d++)
                    for (size_t y = 0; y < C; y++) {
                        uint8_t *dest = out + (d * C + y) * C + a - gx * C;
                        if ((sy * R + r) * C + y < src.height)
                            memcpy(dest, in + (d * C + y) * C + a - col0 - bx * C, b - a);
                        else memset(dest, 0, b - a);       /* source edge padding is outside its logical array */
                    }
            }
            size_t i = 0;
            while (i < bytes && out[i] == 0) i++;
            if (!bad && i < bytes && ZarrShard_append(&builder, temporary, r, out)) bad = 1;
        }
        ZarrShard_close(&old);
        for (size_t i = 0; i < 2; i++) ZarrShard_close(&inputs[i]);
        if (!bad && ZarrShard_finish(&builder, temporary)) bad = 1;
        if (!bad && builder.bytes) {
#ifdef _WIN32
            if (!MoveFileExA(temporary, path, MOVEFILE_REPLACE_EXISTING)) bad = 1;
#else
            if (rename(temporary, path)) bad = 1;
#endif
            if (!bad) written++;
        } else if (!bad && remove(path) && errno != ENOENT) bad = 1;
        free(out); free(in);
        if (bad) {
            fprintf(stderr, "ERROR: cannot merge shard %s\n", path);
            failed = 1;
        }
    }
    if (files_written) *files_written = written;
    return failed ? -1 : 0;
}

int ZarrShard_selftest(void)
{
    char directory[512], path[700], absent[700];
    uint8_t chunks[7 * 73], readback[73];
    size_t written = 0;
    ZarrShard shard;
    int fails = 0;
#define ZS_CHECK(c, message) do { if (!(c)) { fprintf(stderr, "zarr shard selftest: %s\n", message); fails++; } } while (0)
    if (ves_selftest_dir(directory, sizeof directory, "zarr_shard")) return 1;
    snprintf(path, sizeof path, "%s/data", directory);
    snprintf(absent, sizeof absent, "%s/empty", directory);
    for (size_t i = 0; i < sizeof chunks; i++) chunks[i] = (uint8_t)(i * 19 + 7);
    memset(chunks + 73, 0, 73);
    memset(chunks + 5 * 73, 0, 73);
    ZS_CHECK(zs_crc((const uint8_t *)"123456789", 9) == 0xe3069283u, "CRC32C known vector");
    ZS_CHECK(ZarrShard_write(path, chunks, 73, 7, 8, &written) == 0 && written == 5, "sparse partial shard");
    ZS_CHECK(ZarrShard_open(&shard, path, 73, 8) == 1, "open sparse shard");
    for (size_t c = 0; c < 8; c++) {
        memset(readback, 255, sizeof readback);
        ZS_CHECK(ZarrShard_read(&shard, c, readback) == 0, "read inner chunk");
        ZS_CHECK(memcmp(readback, c < 7 ? chunks + c * 73 : chunks + 73, 73) == 0, "exact bytes and missing tail");
    }
    ZS_CHECK(ZarrShard_read(&shard, 8, readback) < 0, "reject out of range");
    ZarrShard_close(&shard);
    ZS_CHECK(ZarrShard_write(path, chunks, 73, 7, 8, NULL) < 0, "refuse replacement");
    ZS_CHECK(ZarrShard_write(path, chunks + 73, 73, 1, 8, NULL) < 0, "empty write cannot reuse existing file");
    ZS_CHECK(ZarrShard_write(absent, chunks + 73, 73, 1, 8, &written) == 0 && written == 0, "empty shard not stored");
    ZS_CHECK(ZarrShard_open(&shard, absent, 73, 8) == 0 && ZarrShard_read(&shard, 7, readback) == 0 &&
             memcmp(readback, chunks + 73, 73) == 0, "missing shard reads fill");
    ZarrShard_close(&shard);
    ZS_CHECK(ZarrShard_open(&shard, path, 74, 8) < 0, "wrong encoded chunk length");
    ZarrShard_close(&shard);
    FILE *file = fopen(path, "r+b");
    ZS_CHECK(file != NULL, "open for corruption control");
    if (file != NULL) {
        zs_seek(file, 5 * 73 + 7);
        fputc(0x23, file);
        fclose(file);
        ZS_CHECK(ZarrShard_open(&shard, path, 73, 8) < 0, "reject corrupt index");
        ZarrShard_close(&shard);
    }
    file = zs_fresh(absent);
    if (file != NULL) { fwrite(chunks, 1, 3, file); fclose(file); }
    ZS_CHECK(ZarrShard_open(&shard, absent, 73, 8) < 0, "reject truncated index");
    ZarrShard_close(&shard);
    ZS_CHECK(ZarrShard_write(path, chunks, SIZE_MAX, 2, 8, NULL) < 0, "reject byte-count overflow");
    ZS_CHECK(ZarrShard_array_meta(directory, 21, 137, 269, 128, 8) == 0, "write standard array metadata");
    Arena_T arena = Arena_new();
    ZarrShardArray array;
    ZS_CHECK(ZarrShard_array_open(arena, directory, &array) == 0 && array.depth == 21 && array.height == 137 &&
             array.width == 269 && array.tile == 128 && array.rows == 1 && array.cols == 8 && array.count == 8,
             "read validated metadata profile");
    snprintf(path, sizeof path, "%s/zarr.json", directory);
    remove(path);
    ZS_CHECK(ZarrShard_array_meta_grid(directory, 21, 9000, 269, 128, 64, 1) == 0 &&
             ZarrShard_array_open(arena, directory, &array) == 0 && array.rows == 64 && array.cols == 1 &&
             array.count == 64, "vertical shard metadata");
    remove(path);
    ZS_CHECK(ZarrShard_array_meta_grid(directory, 3, 267, 281, 128, 2, 4) == 0 &&
             ZarrShard_array_open(arena, directory, &array) == 0 && array.rows == 2 && array.cols == 4 &&
             array.count == 8, "rectangular shard metadata");
    ZS_CHECK(ZarrShard_array_meta_grid(directory, 3, 1, 1, 128, 65, 1) < 0, "reject oversized shard grid");
    ZarrShardBuilder builder;
    snprintf(path, sizeof path, "%s/incremental", directory);
    ZS_CHECK(ZarrShard_begin(&builder, 73, 64) == 0, "begin bounded incremental shard");
    ZS_CHECK(ZarrShard_append(&builder, path, 63, chunks) == 0 &&
             ZarrShard_append(&builder, path, 2, chunks + 2 * 73) == 0 &&
             ZarrShard_append(&builder, path, 1, chunks + 73) == 0, "append sparse out-of-order chunks and explicit zero");
    ZS_CHECK(ZarrShard_append(&builder, path, 2, chunks) < 0, "reject repeated chunk index");
    ZS_CHECK(ZarrShard_finish(&builder, path) == 0 && ZarrShard_finish(&builder, path) < 0,
             "finish incremental shard only once");
    ZS_CHECK(ZarrShard_append(&builder, path, 3, chunks) < 0, "reject append after finish");
    ZS_CHECK(ZarrShard_open(&shard, path, 73, 64) == 1 && shard.lengths[1] == 73 && shard.lengths[62] == 0,
             "incremental index distinguishes present zero from absent");
    for (size_t c = 0; c < 64; c++) {
        const uint8_t *want = c == 63 ? chunks : c == 2 ? chunks + 2 * 73 : chunks + 73;
        ZS_CHECK(ZarrShard_read(&shard, c, readback) == 0 && !memcmp(readback, want, 73), "incremental exact bytes");
    }
    ZarrShard_close(&shard);
    ZS_CHECK(ZarrShard_begin(&builder, 73, 64) == 0 && ZarrShard_append(&builder, path, 0, chunks) < 0,
             "incremental writer refuses existing output");
    ZS_CHECK(ZarrShard_finish(&builder, path) < 0, "failed builder cannot publish an index");
    snprintf(path, sizeof path, "%s/incremental_empty", directory);
    ZS_CHECK(ZarrShard_begin(&builder, 73, 8) == 0 && ZarrShard_finish(&builder, path) == 0 &&
             ZarrShard_open(&shard, path, 73, 8) == 0, "empty incremental shard has no file");
    ZarrShard_close(&shard);
    ZS_CHECK(ZarrShard_begin(&builder, SIZE_MAX, 64) < 0, "reject incremental byte-count overflow");
    Arena_dispose(&arena);
    printf("zarr shard selftest: %s\n", fails ? "FAIL" : "PASS");
    if (!fails) ves_selftest_clean(directory);
#undef ZS_CHECK
    return fails;
}
