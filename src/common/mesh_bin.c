#include "mesh_bin.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

enum {
    MB_VERSION = 1,
    MB_HEADER_BYTES = 64
};

static const unsigned char mb_magic[8] = {
    'V', 'E', 'S', 'M', 'E', 'S', 'H', '1'
};

static int mb_mul_size(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a) return 0;
    *out = a * b;
    return 1;
}

static int mb_add_u64(uint64_t a, uint64_t b, uint64_t *out)
{
    if (b > UINT64_MAX - a) return 0;
    *out = a + b;
    return 1;
}

static int mb_write_bytes(FILE *f, const void *data, size_t nbytes)
{
    return nbytes == 0 || fwrite(data, 1, nbytes, f) == nbytes;
}

static int mb_read_bytes(FILE *f, void *data, size_t nbytes)
{
    return nbytes == 0 || fread(data, 1, nbytes, f) == nbytes;
}

static int mb_suffix_equal(const char *text, const char *suffix)
{
    size_t nt = strlen(text), ns = strlen(suffix);
    if (nt < ns) return 0;
    text += nt - ns;
    for (size_t i = 0; i < ns; i++) {
        unsigned char a = (unsigned char)text[i];
        unsigned char b = (unsigned char)suffix[i];
        if (a >= 'A' && a <= 'Z') a = (unsigned char)(a - 'A' + 'a');
        if (b >= 'A' && b <= 'Z') b = (unsigned char)(b - 'A' + 'a');
        if (a != b) return 0;
    }
    return 1;
}

int MeshBin_companion_path(const char *path, char *out, size_t out_cap)
{
    size_t n = 0, stem = 0;
    const char suffix[] = ".vmesh";
    if (path == NULL || out == NULL || out_cap == 0) return -1;
    n = strlen(path);
    stem = mb_suffix_equal(path, ".obj") ? n - 4
         : mb_suffix_equal(path, suffix) ? n - (sizeof suffix - 1)
         : n;
    if (stem > SIZE_MAX - (sizeof suffix - 1) ||
        stem + (sizeof suffix - 1) + 1 > out_cap)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof suffix);
    return 0;
}

int MeshBin_looks_complete(const char *path)
{
    unsigned char magic[8];
    uint32_t version = 0, endian = 0, flags = 0, header_bytes = 0;
    uint64_t nv = 0, nf = 0, vertices_offset = 0, uv_offset = 0;
    uint64_t faces_offset = 0, expected_faces = 0, expected_end = 0;
    uint64_t vertex_bytes = 0, uv_bytes = 0, face_bytes = 0;
    uint64_t actual_size = 0;
    FILE *f = NULL;
    int ok = 0;

    if (path == NULL) return 0;
    f = fopen(path, "rb");
    if (f == NULL) return 0;
    ok = mb_read_bytes(f, magic, sizeof magic) &&
         mb_read_bytes(f, &version, sizeof version) &&
         mb_read_bytes(f, &endian, sizeof endian) &&
         mb_read_bytes(f, &flags, sizeof flags) &&
         mb_read_bytes(f, &header_bytes, sizeof header_bytes) &&
         mb_read_bytes(f, &nv, sizeof nv) &&
         mb_read_bytes(f, &nf, sizeof nf) &&
         mb_read_bytes(f, &vertices_offset, sizeof vertices_offset) &&
         mb_read_bytes(f, &uv_offset, sizeof uv_offset) &&
         mb_read_bytes(f, &faces_offset, sizeof faces_offset);
    if (!ok || memcmp(magic, mb_magic, sizeof magic) != 0 ||
        version != MB_VERSION || endian != UINT32_C(0x01020304) ||
        (flags & ~MESH_BIN_HAS_UV) != 0 || header_bytes != MB_HEADER_BYTES ||
        nv == 0 || nf == 0 || nv > INT32_MAX ||
        nv > UINT64_MAX / (3u * sizeof(float)) ||
        nf > UINT64_MAX / (3u * sizeof(int32_t)))
        goto done;

    vertex_bytes = nv * 3u * sizeof(float);
    uv_bytes = nv * 2u * sizeof(float);
    face_bytes = nf * 3u * sizeof(int32_t);
    if (vertices_offset != MB_HEADER_BYTES ||
        !mb_add_u64(MB_HEADER_BYTES, vertex_bytes, &expected_faces))
        goto done;
    if (flags & MESH_BIN_HAS_UV) {
        if (uv_offset != expected_faces ||
            !mb_add_u64(expected_faces, uv_bytes, &expected_faces))
            goto done;
    } else if (uv_offset != 0) {
        goto done;
    }
    if (faces_offset != expected_faces ||
        !mb_add_u64(expected_faces, face_bytes, &expected_end))
        goto done;

#ifdef _WIN32
    if (_fseeki64(f, 0, SEEK_END) != 0) goto done;
    {
        __int64 n = _ftelli64(f);
        if (n < 0) goto done;
        actual_size = (uint64_t)n;
    }
#else
    if (fseek(f, 0, SEEK_END) != 0) goto done;
    {
        long n = ftell(f);
        if (n < 0) goto done;
        actual_size = (uint64_t)n;
    }
#endif
    ok = actual_size == expected_end;
done:
    fclose(f);
    return ok ? 1 : 0;
}

static int mb_replace(const char *temporary, const char *destination)
{
#ifdef _WIN32
    return MoveFileExA(temporary, destination,
                       MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)
         ? 0 : -1;
#else
    return rename(temporary, destination);
#endif
}

int MeshBin_write(const char *path, const float *verts, size_t nv,
                  const int32_t *faces, size_t nf, const float *uv)
{
    char temporary[2048];
    FILE *f = NULL;
    uint32_t version = MB_VERSION, endian = UINT32_C(0x01020304);
    uint32_t flags = uv != NULL ? MESH_BIN_HAS_UV : 0;
    uint32_t header_bytes = MB_HEADER_BYTES;
    uint64_t nv64 = (uint64_t)nv, nf64 = (uint64_t)nf;
    uint64_t vertices_offset = MB_HEADER_BYTES, uv_offset = 0, faces_offset = 0;
    uint64_t vertex_bytes = 0, uv_bytes = 0, face_bytes = 0, end_offset = 0;
    size_t vertex_values = 0, uv_values = 0, face_values = 0;
    size_t vertex_nbytes = 0, uv_nbytes = 0, face_nbytes = 0;
    int temporary_n;
    int ok = 0;
    if (path == NULL || verts == NULL || faces == NULL || nv == 0 || nf == 0 ||
        nv > (size_t)INT32_MAX ||
        !mb_mul_size(nv, 3, &vertex_values) ||
        !mb_mul_size(nv, 2, &uv_values) ||
        !mb_mul_size(nf, 3, &face_values) ||
        !mb_mul_size(vertex_values, sizeof(float), &vertex_nbytes) ||
        !mb_mul_size(uv_values, sizeof(float), &uv_nbytes) ||
        !mb_mul_size(face_values, sizeof(int32_t), &face_nbytes))
        return -1;
    for (size_t i = 0; i < face_values; i++)
        if (faces[i] < 0 || (size_t)faces[i] >= nv) return -1;
    if (vertex_values > UINT64_MAX / sizeof(float) ||
        uv_values > UINT64_MAX / sizeof(float) ||
        face_values > UINT64_MAX / sizeof(int32_t))
        return -1;
    vertex_bytes = (uint64_t)vertex_nbytes;
    uv_bytes = uv != NULL ? (uint64_t)uv_nbytes : 0;
    face_bytes = (uint64_t)face_nbytes;
    if (!mb_add_u64(vertices_offset, vertex_bytes, &end_offset)) return -1;
    if (uv != NULL) {
        uv_offset = end_offset;
        if (!mb_add_u64(end_offset, uv_bytes, &end_offset)) return -1;
    }
    faces_offset = end_offset;
    if (!mb_add_u64(end_offset, face_bytes, &end_offset)) return -1;
    (void)end_offset;

    temporary_n = snprintf(temporary, sizeof temporary, "%s.tmp", path);
    if (temporary_n < 0 || (size_t)temporary_n >= sizeof temporary)
        return -1;
    f = fopen(temporary, "wb");
    if (f == NULL) return -1;
    ok = mb_write_bytes(f, mb_magic, sizeof mb_magic) &&
         mb_write_bytes(f, &version, sizeof version) &&
         mb_write_bytes(f, &endian, sizeof endian) &&
         mb_write_bytes(f, &flags, sizeof flags) &&
         mb_write_bytes(f, &header_bytes, sizeof header_bytes) &&
         mb_write_bytes(f, &nv64, sizeof nv64) &&
         mb_write_bytes(f, &nf64, sizeof nf64) &&
         mb_write_bytes(f, &vertices_offset, sizeof vertices_offset) &&
         mb_write_bytes(f, &uv_offset, sizeof uv_offset) &&
         mb_write_bytes(f, &faces_offset, sizeof faces_offset) &&
         mb_write_bytes(f, verts, vertex_nbytes) &&
         mb_write_bytes(f, uv, uv != NULL ? uv_nbytes : 0) &&
         mb_write_bytes(f, faces, face_nbytes) &&
         fflush(f) == 0;
    if (fclose(f) != 0) ok = 0;
    if (!ok || mb_replace(temporary, path) != 0) {
        remove(temporary);
        return -1;
    }
    return 0;
}

static void *mb_alloc(Arena_T arena, int use_arena, size_t nbytes)
{
    if (use_arena) return ARENA_ALLOC(arena, nbytes ? nbytes : 1);
    return malloc(nbytes ? nbytes : 1);
}

static int mb_read(const char *path, Arena_T arena, int use_arena,
                   MeshBinData *out)
{
    unsigned char magic[8];
    uint32_t version = 0, endian = 0, flags = 0, header_bytes = 0;
    uint64_t nv64 = 0, nf64 = 0, vertices_offset = 0, uv_offset = 0;
    uint64_t faces_offset = 0, expected_faces = 0, expected_end = 0;
    size_t vertex_values = 0, uv_values = 0, face_values = 0;
    size_t vertex_nbytes = 0, uv_nbytes = 0, face_nbytes = 0;
    FILE *f = NULL;
    int ok = 0;
    if (path == NULL || out == NULL || (use_arena && arena == NULL)) return -1;
    memset(out, 0, sizeof *out);
    f = fopen(path, "rb");
    if (f == NULL) return -1;
    ok = mb_read_bytes(f, magic, sizeof magic) &&
         mb_read_bytes(f, &version, sizeof version) &&
         mb_read_bytes(f, &endian, sizeof endian) &&
         mb_read_bytes(f, &flags, sizeof flags) &&
         mb_read_bytes(f, &header_bytes, sizeof header_bytes) &&
         mb_read_bytes(f, &nv64, sizeof nv64) &&
         mb_read_bytes(f, &nf64, sizeof nf64) &&
         mb_read_bytes(f, &vertices_offset, sizeof vertices_offset) &&
         mb_read_bytes(f, &uv_offset, sizeof uv_offset) &&
         mb_read_bytes(f, &faces_offset, sizeof faces_offset);
    if (!ok || memcmp(magic, mb_magic, sizeof magic) != 0 ||
        version != MB_VERSION || endian != UINT32_C(0x01020304) ||
        (flags & ~MESH_BIN_HAS_UV) != 0 || header_bytes != MB_HEADER_BYTES ||
        nv64 == 0 || nf64 == 0 || nv64 > INT32_MAX ||
        nv64 > SIZE_MAX || nf64 > SIZE_MAX)
        goto done;
    out->nv = (size_t)nv64; out->nf = (size_t)nf64;
    if (!mb_mul_size(out->nv, 3, &vertex_values) ||
        !mb_mul_size(out->nv, 2, &uv_values) ||
        !mb_mul_size(out->nf, 3, &face_values) ||
        !mb_mul_size(vertex_values, sizeof(float), &vertex_nbytes) ||
        !mb_mul_size(uv_values, sizeof(float), &uv_nbytes) ||
        !mb_mul_size(face_values, sizeof(int32_t), &face_nbytes) ||
        vertex_values > UINT64_MAX / sizeof(float) ||
        uv_values > UINT64_MAX / sizeof(float) ||
        face_values > UINT64_MAX / sizeof(int32_t))
        goto done;
    if (!mb_add_u64(MB_HEADER_BYTES,
                    (uint64_t)vertex_values * sizeof(float), &expected_faces))
        goto done;
    if (flags & MESH_BIN_HAS_UV) {
        if (uv_offset != expected_faces ||
            !mb_add_u64(expected_faces,
                        (uint64_t)uv_values * sizeof(float), &expected_faces))
            goto done;
    } else if (uv_offset != 0) {
        goto done;
    }
    if (vertices_offset != MB_HEADER_BYTES || faces_offset != expected_faces ||
        !mb_add_u64(expected_faces,
                    (uint64_t)face_values * sizeof(int32_t), &expected_end))
        goto done;
    (void)expected_end;
    out->verts = (float *)mb_alloc(arena, use_arena, vertex_nbytes);
    if (flags & MESH_BIN_HAS_UV)
        out->uv = (float *)mb_alloc(arena, use_arena, uv_nbytes);
    out->faces = (int32_t *)mb_alloc(arena, use_arena, face_nbytes);
    if (out->verts == NULL || out->faces == NULL ||
        ((flags & MESH_BIN_HAS_UV) && out->uv == NULL))
        goto done;
    ok = mb_read_bytes(f, out->verts, vertex_nbytes) &&
         (!(flags & MESH_BIN_HAS_UV) ||
           mb_read_bytes(f, out->uv, uv_nbytes)) &&
         mb_read_bytes(f, out->faces, face_nbytes) &&
         fgetc(f) == EOF;
    if (!ok) goto done;
    for (size_t i = 0; i < face_values; i++)
        if (out->faces[i] < 0 || (size_t)out->faces[i] >= out->nv) {
            ok = 0; break;
        }
done:
    fclose(f);
    if (!ok) {
        if (!use_arena) MeshBin_dispose(out);
        else memset(out, 0, sizeof *out);
    }
    return ok ? 0 : -1;
}

int MeshBin_read_arena(Arena_T arena, const char *path, MeshBinData *out)
{
    return mb_read(path, arena, 1, out);
}

int MeshBin_read_malloc(const char *path, MeshBinData *out)
{
    return mb_read(path, NULL, 0, out);
}

void MeshBin_dispose(MeshBinData *mesh)
{
    if (mesh == NULL) return;
    free(mesh->verts); free(mesh->uv); free(mesh->faces);
    memset(mesh, 0, sizeof *mesh);
}

int MeshBin_selftest(void)
{
    const char *path = "mesh_bin_selftest.vmesh";
    float vertices[12] = {0,0,0, 1,0,0, 1,1,0, 0,1,0};
    float uv[8] = {0,0, 1,0, 1,1, 0,1};
    int32_t faces[6] = {0,1,2, 0,2,3};
    MeshBinData readback = {0};
    char companion[64];
    int failed = MeshBin_companion_path("foo.obj", companion,
                                        sizeof companion) != 0 ||
                 strcmp(companion, "foo.vmesh") != 0 ||
                 MeshBin_write(path, vertices, 4, faces, 2, uv) != 0 ||
                 !MeshBin_looks_complete(path) ||
                 MeshBin_read_malloc(path, &readback) != 0;
    if (!failed) {
        failed = readback.nv != 4 || readback.nf != 2 ||
                 readback.uv == NULL ||
                 memcmp(readback.verts, vertices, sizeof vertices) != 0 ||
                 memcmp(readback.uv, uv, sizeof uv) != 0 ||
                 memcmp(readback.faces, faces, sizeof faces) != 0;
        MeshBin_dispose(&readback);
    }
    remove(path);
    fprintf(stderr, "[mesh_bin selftest] %s\n", failed ? "FAILED" : "PASSED");
    return failed ? 1 : 0;
}
