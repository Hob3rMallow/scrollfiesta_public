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
    int header_read = ok;
    ok = 0;
    if (!header_read || memcmp(magic, mb_magic, sizeof magic) != 0 ||
        version != MB_VERSION || endian != UINT32_C(0x01020304) ||
        (flags & ~(MESH_BIN_HAS_UV|MESH_BIN_UV_F64)) != 0 ||
        ((flags & MESH_BIN_UV_F64) && !(flags & MESH_BIN_HAS_UV)) || header_bytes != MB_HEADER_BYTES ||
        nv == 0 || nf == 0 || nv > INT32_MAX ||
        nv > UINT64_MAX / (3u * sizeof(float)) ||
        nf > UINT64_MAX / (3u * sizeof(int32_t)))
        goto done;

    vertex_bytes = nv * 3u * sizeof(float);
    uv_bytes = nv * 2u * ((flags & MESH_BIN_UV_F64) ? sizeof(double) : sizeof(float));
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

static int mb_write(const char *path, const float *verts, size_t nv,
                     const int32_t *faces, size_t nf, const void *uv, size_t uv_size)
{
    char temporary[2048];
    FILE *f = NULL;
    uint32_t version = MB_VERSION, endian = UINT32_C(0x01020304);
    uint32_t flags = uv != NULL ? MESH_BIN_HAS_UV|(uv_size == sizeof(double) ? MESH_BIN_UV_F64 : 0) : 0;
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
        !mb_mul_size(uv_values, uv_size, &uv_nbytes) ||
        !mb_mul_size(face_values, sizeof(int32_t), &face_nbytes))
        return -1;
    for (size_t i = 0; i < face_values; i++)
        if (faces[i] < 0 || (size_t)faces[i] >= nv) return -1;
    if (vertex_values > UINT64_MAX / sizeof(float) ||
        uv_values > UINT64_MAX / uv_size ||
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

int MeshBin_write(const char *path, const float *verts, size_t nv,
                  const int32_t *faces, size_t nf, const float *uv)
{
    return mb_write(path,verts,nv,faces,nf,uv,sizeof(float));
}

int MeshBin_write_uv64(const char *path, const float *verts, size_t nv,
                        const int32_t *faces, size_t nf, const double *uv)
{
    return mb_write(path,verts,nv,faces,nf,uv,sizeof(double));
}

static void *mb_alloc(Arena_T arena, int use_arena, size_t nbytes)
{
    if (use_arena) return ARENA_ALLOC(arena, nbytes ? nbytes : 1);
    return malloc(nbytes ? nbytes : 1);
}

struct MeshBinStream_T {
    FILE *f;
    FILE *fu,*ff;          /* the uv and face sections on their own handles: every section writes sequentially */
    char path[2048],temporary[2048];
    uint64_t nv,nf,at_v,at_f,uvoff,foff;
    size_t uv_size;
    int failed;
};

static int mb_seek(FILE *f, uint64_t offset)
{
    if (offset>INT64_MAX) return -1;
#ifdef _WIN32
    return _fseeki64(f,(__int64)offset,SEEK_SET);
#else
    return fseeko(f,(off_t)offset,SEEK_SET);
#endif
}

static int mb_stream_open(Arena_T arena, const char *path, size_t nv, size_t nf,
                            int has_uv, size_t uv_size, MeshBinStream_T *out)
{
    MeshBinStream_T w=NULL;
    uint32_t fields[4]={MB_VERSION,UINT32_C(0x01020304),0,MB_HEADER_BYTES};
    uint64_t layout[5]={0};
    if (out==NULL) return -1;
    *out=NULL;
    if (arena==NULL || path==NULL || strlen(path)>=2040 || nv==0 || nf==0 ||
        nv>INT32_MAX || nf>(INT64_MAX-64-(uint64_t)nv*(12+2*uv_size))/12) return -1;
    w=(MeshBinStream_T)ARENA_CALLOC(arena,1,sizeof *w);
    memcpy(w->path,path,strlen(path)+1);
    snprintf(w->temporary,sizeof w->temporary,"%s.tmp",path);
    w->nv=nv; w->nf=nf; w->uv_size=uv_size;
    w->uvoff=has_uv ? MB_HEADER_BYTES+(uint64_t)nv*12 : 0;
    w->foff=MB_HEADER_BYTES+(uint64_t)nv*(12+(has_uv ? 2*uv_size : 0));
    fields[2]=has_uv ? MESH_BIN_HAS_UV|(uv_size == sizeof(double) ? MESH_BIN_UV_F64 : 0) : 0;
    layout[0]=nv; layout[1]=nf; layout[2]=MB_HEADER_BYTES;
    layout[3]=w->uvoff; layout[4]=w->foff;
    w->f=fopen(w->temporary,"wb");
    if (w->f==NULL) return -1;
    if (!mb_write_bytes(w->f,mb_magic,8) || !mb_write_bytes(w->f,fields,sizeof fields) ||
        !mb_write_bytes(w->f,layout,sizeof layout)) {
        MeshBin_stream_close(&w,0); return -1;
    }
    /* Appends used to seek between the three sections for every part, and
     * each seek flushed the buffer: a 1.2 GB export was ~70,000 small writes.
     * Two more handles keep each section sequential; if they cannot be
     * opened the single-handle seeks remain. The bytes are the same. */
    if (fflush(w->f)==0) {
        w->ff=fopen(w->temporary,"r+b");
        if (w->ff && mb_seek(w->ff,w->foff)!=0) { fclose(w->ff); w->ff=NULL; }
        if (w->ff && w->uvoff) {
            w->fu=fopen(w->temporary,"r+b");
            if (w->fu && mb_seek(w->fu,w->uvoff)!=0) { fclose(w->fu); w->fu=NULL; }
            if (!w->fu) { fclose(w->ff); w->ff=NULL; }
        }
        if (w->ff) {
            setvbuf(w->f,NULL,_IOFBF,(size_t)1<<20); setvbuf(w->ff,NULL,_IOFBF,(size_t)1<<20);
            if (w->fu) setvbuf(w->fu,NULL,_IOFBF,(size_t)1<<20);
        }
    }
    *out=w;
    return 0;
}

int MeshBin_stream_open(Arena_T arena, const char *path, size_t nv, size_t nf,
                        int has_uv, MeshBinStream_T *out)
{
    return mb_stream_open(arena,path,nv,nf,has_uv,sizeof(float),out);
}

int MeshBin_stream_open_uv64(Arena_T arena, const char *path, size_t nv, size_t nf,
                             MeshBinStream_T *out)
{
    return mb_stream_open(arena,path,nv,nf,1,sizeof(double),out);
}

int MeshBin_stream_append(MeshBinStream_T w, const MeshBinData *m)
{
    int32_t buffer[12288]={0};
    if (w==NULL || m==NULL || w->f==NULL || w->failed) return -1;
    const void *uv = w->uv_size == sizeof(double) ? (const void *)m->uv64 : (const void *)m->uv;
    if (m->nv>w->nv-w->at_v || m->nf>w->nf-w->at_f ||
        (m->nv && (m->verts==NULL || (w->uvoff && uv==NULL))) ||
        (m->nf && m->faces==NULL)) goto fail;
    for (size_t i=0; i<m->nf*3; i++)
        if (m->faces[i]<0 || (size_t)m->faces[i]>=m->nv) goto fail;
    FILE *faces=w->ff ? w->ff : w->f;
    if ((!w->ff && mb_seek(w->f,MB_HEADER_BYTES+w->at_v*12)!=0) ||
        !mb_write_bytes(w->f,m->verts,m->nv*3*sizeof(float))) goto fail;
    if (w->uvoff && ((!w->fu && mb_seek(w->f,w->uvoff+w->at_v*2*w->uv_size)!=0) ||
                     !mb_write_bytes(w->fu ? w->fu : w->f,uv,m->nv*2*w->uv_size))) goto fail;
    if (!w->ff && mb_seek(w->f,w->foff+w->at_f*12)!=0) goto fail;
    for (size_t i=0; i<m->nf*3; ) {
        size_t n=m->nf*3-i;
        if (n>sizeof buffer/sizeof *buffer) n=sizeof buffer/sizeof *buffer;
        for (size_t k=0; k<n; k++) buffer[k]=m->faces[i+k]+(int32_t)w->at_v;
        if (!mb_write_bytes(faces,buffer,n*sizeof(int32_t))) goto fail;
        i+=n;
    }
    w->at_v+=m->nv; w->at_f+=m->nf;
    return 0;
fail:
    w->failed=1;
    return -1;
}

int MeshBin_stream_close(MeshBinStream_T *writer, int publish)
{
    MeshBinStream_T w=writer ? *writer : NULL;
    int ok=0;
    if (w==NULL) return -1;
    *writer=NULL;
    ok=publish && !w->failed && w->at_v==w->nv && w->at_f==w->nf;
    if (w->fu && fclose(w->fu)!=0) ok=0;
    if (w->ff && fclose(w->ff)!=0) ok=0;
    if (w->f && fclose(w->f)!=0) ok=0;
    w->f=w->fu=w->ff=NULL;
    if (ok && mb_replace(w->temporary,w->path)!=0) ok=0;
    if (!ok) remove(w->temporary);
    return ok || !publish ? 0 : -1;
}

static int mb_read(const char *path, Arena_T arena, int use_arena, int accept_uv64,
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
    int header_read = ok;
    ok = 0;
    if (!header_read || memcmp(magic, mb_magic, sizeof magic) != 0 ||
        version != MB_VERSION || endian != UINT32_C(0x01020304) ||
        (flags & ~(MESH_BIN_HAS_UV|MESH_BIN_UV_F64)) != 0 ||
        ((flags & MESH_BIN_UV_F64) && (!(flags & MESH_BIN_HAS_UV) || !accept_uv64)) || header_bytes != MB_HEADER_BYTES ||
        nv64 == 0 || nf64 == 0 || nv64 > INT32_MAX ||
        nv64 > SIZE_MAX || nf64 > SIZE_MAX)
        goto done;
    out->nv = (size_t)nv64; out->nf = (size_t)nf64;
    size_t uv_size = (flags & MESH_BIN_UV_F64) ? sizeof(double) : sizeof(float);
    if (!mb_mul_size(out->nv, 3, &vertex_values) ||
        !mb_mul_size(out->nv, 2, &uv_values) ||
        !mb_mul_size(out->nf, 3, &face_values) ||
        !mb_mul_size(vertex_values, sizeof(float), &vertex_nbytes) ||
        !mb_mul_size(uv_values, uv_size, &uv_nbytes) ||
        !mb_mul_size(face_values, sizeof(int32_t), &face_nbytes) ||
        vertex_values > UINT64_MAX / sizeof(float) ||
        uv_values > UINT64_MAX / uv_size ||
        face_values > UINT64_MAX / sizeof(int32_t))
        goto done;
    if (!mb_add_u64(MB_HEADER_BYTES,
                    (uint64_t)vertex_values * sizeof(float), &expected_faces))
        goto done;
    if (flags & MESH_BIN_HAS_UV) {
        if (uv_offset != expected_faces ||
            !mb_add_u64(expected_faces,
                        (uint64_t)uv_values * uv_size, &expected_faces))
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
    if (flags & MESH_BIN_HAS_UV) {
        if (flags & MESH_BIN_UV_F64) out->uv64 = (double *)mb_alloc(arena,use_arena,uv_nbytes);
        else out->uv = (float *)mb_alloc(arena,use_arena,uv_nbytes);
    }
    out->faces = (int32_t *)mb_alloc(arena, use_arena, face_nbytes);
    if (out->verts == NULL || out->faces == NULL ||
        ((flags & MESH_BIN_HAS_UV) && out->uv == NULL && out->uv64 == NULL))
        goto done;
    ok = mb_read_bytes(f, out->verts, vertex_nbytes) &&
         (!(flags & MESH_BIN_HAS_UV) ||
           mb_read_bytes(f,out->uv64 ? (void *)out->uv64 : (void *)out->uv,uv_nbytes)) &&
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
    return mb_read(path, arena, 1, 0, out);
}

int MeshBin_read_malloc(const char *path, MeshBinData *out)
{
    return mb_read(path, NULL, 0, 0, out);
}

int MeshBin_read_precise_arena(Arena_T arena, const char *path, MeshBinData *out)
{
    return mb_read(path,arena,1,1,out);
}

int MeshBin_read_precise_malloc(const char *path, MeshBinData *out)
{
    return mb_read(path,NULL,0,1,out);
}

void MeshBin_dispose(MeshBinData *mesh)
{
    if (mesh == NULL) return;
    free(mesh->verts); free(mesh->uv); free(mesh->uv64); free(mesh->faces);
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
    /* A successful 64-byte read is not a valid header. Exercise every
     * structural field through both the cheap validator and full readers. */
    {
        const long offsets[] = {0,8,12,16,20,24,32,40,48,56};
        Arena_T arena = Arena_new();
        for (size_t k = 0; k < sizeof offsets/sizeof *offsets; k++) {
            if (MeshBin_write(path,vertices,4,faces,2,uv) != 0) { failed++; continue; }
            FILE *fp = fopen(path,"r+b");
            int byte = EOF, written = 0;
            if (fp && fseek(fp,offsets[k],SEEK_SET) == 0) byte = fgetc(fp);
            if (fp && byte != EOF && fseek(fp,offsets[k],SEEK_SET) == 0)
                written = fputc(byte^128,fp) != EOF;
            if (fp && fclose(fp) != 0) written = 0;
            if (!written) { failed++; continue; }
            if (MeshBin_looks_complete(path)) failed++;
            if (MeshBin_read_malloc(path,&readback) == 0) failed++;
            if (readback.nv || readback.nf || readback.verts || readback.uv || readback.uv64 || readback.faces) failed++;
            MeshBin_dispose(&readback);
            Arena_Mark mark = Arena_save(arena);
            if (MeshBin_read_arena(arena,path,&readback) == 0) failed++;
            if (readback.nv || readback.nf || readback.verts || readback.uv || readback.uv64 || readback.faces) failed++;
            Arena_restore(arena,mark);
        }
        Arena_dispose(&arena);
    }
    remove(path);
    {
        Arena_T arena=Arena_new();
        MeshBinStream_T writer=NULL;
        MeshBinData chunk={vertices,uv,faces,4,2};
        failed+=MeshBin_stream_open(arena,path,8,4,1,&writer)!=0;
        failed+=MeshBin_stream_append(writer,&chunk)!=0;
        failed+=MeshBin_looks_complete(path)!=0;
        failed+=MeshBin_stream_append(writer,&chunk)!=0;
        failed+=MeshBin_stream_close(&writer,1)!=0;
        failed+=MeshBin_read_malloc(path,&readback)!=0;
        if (readback.nv==8 && readback.nf==4 && readback.uv) {
            failed+=memcmp(readback.verts+12,vertices,sizeof vertices)!=0;
            failed+=memcmp(readback.uv+8,uv,sizeof uv)!=0;
            for (int i=0; i<6; i++) failed+=readback.faces[6+i]!=faces[i]+4;
        } else failed++;
        MeshBin_dispose(&readback);
        failed+=MeshBin_stream_open(arena,path,8,4,1,&writer)!=0;
        failed+=MeshBin_stream_append(writer,&chunk)!=0;
        failed+=MeshBin_stream_close(&writer,1)==0; /* partial cannot publish */
        failed+=!MeshBin_looks_complete(path); /* previous completed file survives */
        remove(path);
        Arena_dispose(&arena);
    }
    {
        /* Fine edges remain distinct at absolute coordinates which float32
         * cannot resolve. Both contiguous and bounded-memory exports retain
         * every bit; legacy readers must refuse an implicit precision loss. */
        double precise[] = {1e9+.1,2e9+.01,1e9+.2,2e9+.01,
                            1e9+.2,2e9+.02,1e9+.1,2e9+.02};
        Arena_T arena = Arena_new();
        MeshBinStream_T writer = NULL;
        MeshBinData chunk = {vertices,NULL,faces,4,2,precise};
        failed += MeshBin_write_uv64(path,vertices,4,faces,2,precise) != 0;
        failed += !MeshBin_looks_complete(path);
        failed += MeshBin_read_malloc(path,&readback) == 0;
        if (readback.verts || readback.uv || readback.uv64 || readback.faces || readback.nv || readback.nf) failed++;
        MeshBin_dispose(&readback);
        failed += MeshBin_read_precise_malloc(path,&readback) != 0;
        if (readback.nv == 4 && readback.nf == 2 && !readback.uv && readback.uv64) {
            failed += memcmp(readback.uv64,precise,sizeof precise) != 0;
            failed += memcmp(readback.verts,vertices,sizeof vertices) != 0;
            failed += memcmp(readback.faces,faces,sizeof faces) != 0;
            failed += !(MeshBin_uv(&readback,2)-MeshBin_uv(&readback,0) > .09);
        } else failed++;
        MeshBin_dispose(&readback);
        failed += MeshBin_stream_open_uv64(arena,path,8,4,&writer) != 0;
        failed += MeshBin_stream_append(writer,&chunk) != 0;
        failed += MeshBin_stream_append(writer,&chunk) != 0;
        failed += MeshBin_stream_close(&writer,1) != 0;
        failed += !MeshBin_looks_complete(path);
        failed += MeshBin_read_precise_arena(arena,path,&readback) != 0;
        if (readback.nv == 8 && readback.nf == 4 && !readback.uv && readback.uv64) {
            failed += memcmp(readback.uv64+8,precise,sizeof precise) != 0;
            for (int k = 0; k < 6; k++) failed += readback.faces[6+k] != faces[k]+4;
        } else failed++;
        /* The precise reader also retains existing float32 files exactly. */
        failed += MeshBin_write(path,vertices,4,faces,2,uv) != 0;
        failed += MeshBin_read_precise_arena(arena,path,&readback) != 0;
        if (readback.uv && !readback.uv64) failed += memcmp(readback.uv,uv,sizeof uv) != 0;
        else failed++;
        memset(&readback,0,sizeof readback);
        remove(path);
        Arena_dispose(&arena);
    }
    fprintf(stderr, "[mesh_bin selftest] %s\n", failed ? "FAILED" : "PASSED");
    return failed ? 1 : 0;
}
