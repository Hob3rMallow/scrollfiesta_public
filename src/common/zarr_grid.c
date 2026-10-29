#include "zarr_grid.h"
#include "arena.h"
#include "json_read.h"
#include "tiff_io.h"
#include "ves_platform.h"

#include <blosc2.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { ZG_SLOTS = 8, ZG_PATH = 2048, ZG_CUBE = 128 };
typedef struct ZgSlot {
    long key[3];
    unsigned long long age;
    uint8_t *data;
} ZgSlot;
typedef struct ZgArray {
    Arena_T arena;
    char root[ZG_PATH];
    long shape[3], chunk[3];
    size_t bytes, encoded_capacity;
    int blosc;
    uint8_t *encoded;
    ZgSlot slot[ZG_SLOTS];
    unsigned long long clock, reads;
} ZgArray;

static int zg_string(const JsonValue *v, const char *expected)
{
    const char *s = Json_as_string(v);
    return s != NULL && strcmp(s, expected) == 0;
}

static int zg_open(Arena_T arena, const char *root, ZgArray *a)
{
    char path[ZG_PATH];
    const char *err = NULL;
    const JsonValue *j, *codec;
    memset(a, 0, sizeof *a);
    a->arena = arena;
    if (strlen(root) > 1024) return -1;
    snprintf(a->root, sizeof a->root, "%s", root);
    snprintf(path, sizeof path, "%s/.zarray", root);
    j = Json_parse_file(arena, path, &err);
    if (j == NULL || Json_member_long(j, "zarr_format", -1) != 2 ||
        !zg_string(Json_object_get(j, "dtype"), "|u1") ||
        !zg_string(Json_object_get(j, "order"), "C") ||
        !zg_string(Json_object_get(j, "dimension_separator"), "/") ||
        Json_type(Json_object_get(j, "filters")) != JSON_NULL ||
        Json_member_double(j, "fill_value", -1) != 0) return -1;
    a->bytes = 1;
    for (int k = 0; k < 3; k++) {
        const JsonValue *shape = Json_object_get(j, "shape");
        const JsonValue *chunk = Json_object_get(j, "chunks");
        const JsonValue *sv = Json_array_get(shape, (size_t)k);
        const JsonValue *cv = Json_array_get(chunk, (size_t)k);
        if (Json_array_len(shape) != 3 || Json_array_len(chunk) != 3 ||
            Json_type(sv) != JSON_NUMBER || Json_type(cv) != JSON_NUMBER)
            return -1;
        if (!(Json_as_double(sv,-1) >= 1 && Json_as_double(sv,-1) <= 99999) ||
            !(Json_as_double(cv,-1) >= 1 && Json_as_double(cv,-1) <= 256))
            return -1;
        a->shape[k] = Json_as_long(sv, -1);
        a->chunk[k] = Json_as_long(cv, -1);
        if (a->shape[k] < 1 || a->shape[k] > 99999 ||
            a->chunk[k] < 1 || a->chunk[k] > 256 ||
            Json_as_double(sv, -1) != (double)a->shape[k] ||
            Json_as_double(cv, -1) != (double)a->chunk[k]) return -1;
        a->bytes *= (size_t)a->chunk[k];
    }
    codec = Json_object_get(j, "compressor");
    if (codec == NULL) return -1;
    if (Json_type(codec) == JSON_NULL) a->blosc = 0;
    else if (Json_type(codec) == JSON_OBJECT &&
             zg_string(Json_object_get(codec, "id"), "blosc")) a->blosc = 1;
    else return -1;
    a->encoded_capacity = a->bytes + BLOSC2_MAX_OVERHEAD;
    a->encoded = ARENA_ALLOC(arena, (size_t)a->encoded_capacity);
    for (int i = 0; i < ZG_SLOTS; i++) {
        a->slot[i].age = 0;
        a->slot[i].data = ARENA_ALLOC(arena, (size_t)a->bytes);
    }
    return 0;
}

static const uint8_t *zg_chunk(ZgArray *a, const long key[3])
{
    ZgSlot *slot = &a->slot[0];
    char path[ZG_PATH];
    FILE *f;
    size_t bytes, decoded = 0;
    int extra, bad;
    for (int i = 0; i < ZG_SLOTS; i++) {
        ZgSlot *s = &a->slot[i];
        if (s->age && memcmp(s->key, key, sizeof s->key) == 0) {
            s->age = ++a->clock;
            return s->data;
        }
        if (s->age < slot->age) slot = s;
    }
    snprintf(path, sizeof path, "%s/%ld/%ld/%ld", a->root,
             key[0], key[1], key[2]);
    f = fopen(path, "rb");
    if (f == NULL) {
        fprintf(stderr, "[zarr grid] missing/unreadable source chunk: %s\n", path);
        return NULL;
    }
    bytes = fread(a->encoded, 1, a->encoded_capacity, f);
    extra = fgetc(f);
    bad = ferror(f);
    if (fclose(f)) bad = 1;
    if (bad || extra != EOF) return NULL;
    slot->age = 0;
    if (a->blosc) {
        if (blosc1_cbuffer_validate(a->encoded, bytes, &decoded) < 0 ||
            decoded != a->bytes ||
            blosc2_decompress(a->encoded, (int32_t)bytes, slot->data,
                              (int32_t)a->bytes) != (int)a->bytes) {
            fprintf(stderr, "[zarr grid] corrupt/incompatible Blosc chunk: %s\n", path);
            return NULL;
        }
    } else {
        if (bytes != a->bytes) {
            fprintf(stderr, "[zarr grid] wrong raw chunk length: %s\n", path);
            return NULL;
        }
        memcpy(slot->data, a->encoded, bytes);
    }
    memcpy(slot->key, key, sizeof slot->key);
    slot->age = ++a->clock;
    a->reads++;
    return slot->data;
}

static int zg_cube(ZgArray *a, const long origin[3], int side,
                    int binary, uint8_t *out)
{
    long first[3], last[3];
    if (side < 1 || side > ZG_CUBE || (binary != 0 && binary != 1)) return -1;
    for (int k = 0; k < 3; k++) {
        if (origin[k] < 0 || origin[k] > a->shape[k] - side) return -1;
        first[k] = origin[k] / a->chunk[k];
        last[k] = (origin[k] + side - 1) / a->chunk[k];
    }
    for (long iz = first[0]; iz <= last[0]; iz++)
    for (long iy = first[1]; iy <= last[1]; iy++)
    for (long ix = first[2]; ix <= last[2]; ix++) {
        const long key[3] = { iz, iy, ix };
        const uint8_t *src = zg_chunk(a, key);
        long lo[3], hi[3], base[3];
        if (src == NULL) return -1;
        for (int k = 0; k < 3; k++) {
            base[k] = key[k] * a->chunk[k];
            lo[k] = origin[k] > base[k] ? origin[k] : base[k];
            hi[k] = origin[k] + side < base[k] + a->chunk[k]
                    ? origin[k] + side : base[k] + a->chunk[k];
        }
        for (long z = lo[0]; z < hi[0]; z++)
        for (long y = lo[1]; y < hi[1]; y++) {
            size_t si = ((size_t)(z-base[0])*(size_t)a->chunk[1] +
                          (size_t)(y-base[1]))*(size_t)a->chunk[2] +
                          (size_t)(lo[2]-base[2]);
            size_t di = ((size_t)(z-origin[0])*(size_t)side +
                          (size_t)(y-origin[1]))*(size_t)side +
                          (size_t)(lo[2]-origin[2]);
            size_t n = (size_t)(hi[2]-lo[2]);
            if (binary) {
                for (size_t x = 0; x < n; x++) out[di+x] = src[si+x] ? 255 : 0;
            } else memcpy(out+di, src+si, n);
        }
    }
    return 0;
}

static int zg_exists(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL) return errno != ENOENT;
    fclose(f);
    return 1;
}

static int zg_export(ZgArray *a, const char *out_dir, const long bbox[6],
                      int side, int binary)
{
    Arena_T scratch = NULL;
    uint8_t *cube;
    const size_t count = (size_t)side*(size_t)side*(size_t)side;
    char path[ZG_PATH], temporary[ZG_PATH], report[ZG_PATH];
    size_t written = 0;
    int result = -1;
    if (side < 1 || side > ZG_CUBE || strlen(out_dir) > 1024 ||
        (binary != 0 && binary != 1)) return -1;
    for (int k = 0; k < 3; k++)
        if (bbox[2*k] < 0 || bbox[2*k] >= bbox[2*k+1] ||
            bbox[2*k+1] > a->shape[k] || bbox[2*k] % side ||
            bbox[2*k+1] % side) return -1;
    snprintf(report, sizeof report, "%s/carve_complete.json", out_dir);
    if (zg_exists(report)) return -1;
    /* Check every destination before creating any cube. */
    for (long z = bbox[0]; z < bbox[1]; z += side)
    for (long y = bbox[2]; y < bbox[3]; y += side)
    for (long x = bbox[4]; x < bbox[5]; x += side) {
        snprintf(path, sizeof path, "%s/z%05ld_y%05ld_x%05ld.tif", out_dir,z,y,x);
        snprintf(temporary, sizeof temporary, "%s.partial", path);
        if (zg_exists(path) || zg_exists(temporary)) {
            fprintf(stderr, "[zarr grid] refusing to overwrite: %s\n", path);
            return -1;
        }
    }
    if (ves_mkdir(out_dir) != 0 && errno != EEXIST) return -1;
    cube = ARENA_ALLOC(a->arena, (size_t)count);
    scratch = Arena_new();
    for (long z = bbox[0]; z < bbox[1]; z += side) {
        for (long y = bbox[2]; y < bbox[3]; y += side)
        for (long x = bbox[4]; x < bbox[5]; x += side) {
            const long origin[3] = { z,y,x };
            uint8_t *readback = NULL;
            int d = 0, h = 0, w = 0;
            snprintf(path,sizeof path,"%s/z%05ld_y%05ld_x%05ld.tif",out_dir,z,y,x);
            snprintf(temporary,sizeof temporary,"%s.partial",path);
            if (zg_cube(a,origin,side,binary,cube)) goto done;
            if (TiffIO_save(temporary,cube,side,side,side) ||
                TiffIO_load(scratch,temporary,&readback,&d,&h,&w) ||
                d != side || h != side || w != side ||
                memcmp(cube,readback,count) != 0 || rename(temporary,path)) {
                fprintf(stderr,"[zarr grid] write/readback failed: %s\n",path);
                goto done;
            }
            Arena_free(scratch);
            written++;
        }
        fprintf(stderr,"[zarr grid] z%05ld: %zu cubes copied and read back\n",z,written);
    }
    {
        FILE *f = fopen(report,"w");
        int bad;
        if (f == NULL) goto done;
        bad = fprintf(f,"{\n  \"schema\":\"native-zarr-grid-v1\",\n"
            "  \"complete\":true,\n  \"bbox_l0_zyx\":[%ld,%ld,%ld,%ld,%ld,%ld],\n"
            "  \"cube_size\":%d,\n  \"cubes_written\":%zu,\n"
            "  \"binary_nonzero_to_255\":%s,\n  \"missing_chunk_policy\":\"fail\",\n"
            "  \"all_output_voxels_readback_equal\":true,\n"
            "  \"source_chunk_reads\":%llu\n}\n",
            bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5],side,written,
            binary?"true":"false",a->reads) < 0;
        if (fclose(f)) bad = 1;
        if (bad) goto done;
    }
    result = 0;
done:
    Arena_dispose(&scratch);
    return result;
}

int ZarrGrid_export(const char *array_dir, const char *out_dir,
                    const long bbox[6], int binary)
{
    Arena_T arena = Arena_new();
    ZgArray a;
    int result = -1;
    if (zg_open(arena,array_dir,&a))
        fprintf(stderr,"[zarr grid] unsupported/invalid array metadata: %s\n",array_dir);
    else result = zg_export(&a,out_dir,bbox,ZG_CUBE,binary);
    Arena_dispose(&arena);
    return result;
}

#include "zarr_grid_test.inc"
