#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#endif
#include "closeup_seed.h"
#include "../common/ves_platform.h"
#include "../common/mesh_bin.h"
#include "../common/json_read.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif

enum { CS_BLOCK = 256, CS_BATCH = 262144 };
typedef struct CsMap {
#ifdef _WIN32
    HANDLE file, mapping;
#else
    int fd;
#endif
    void *base;
    size_t length;
    uint64_t bytes;
} CsMap;
typedef struct CsBounds { double lo[2], hi[2]; } CsBounds;
typedef struct CsIndex { char magic[8]; uint64_t nv, nf, bytes, block; } CsIndex;

static void cs_unview(CsMap *m)
{
    if (!m->base) return;
#ifdef _WIN32
    UnmapViewOfFile(m->base);
#else
    munmap(m->base, m->length);
#endif
    m->base = NULL;
}

static int cs_open(CsMap *m, const char *path)
{
    memset(m, 0, sizeof *m);
#ifdef _WIN32
    m->file = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE,
                         NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (m->file == INVALID_HANDLE_VALUE) { m->file = NULL; return -1; }
    LARGE_INTEGER size;
    if (!GetFileSizeEx(m->file, &size) || size.QuadPart <= 0) { CloseHandle(m->file); m->file = NULL; return -1; }
    m->bytes = (uint64_t)size.QuadPart;
    m->mapping = CreateFileMappingA(m->file, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!m->mapping) { CloseHandle(m->file); m->file = NULL; return -1; }
#else
    m->fd = open(path, O_RDONLY);
    struct stat st;
    if (m->fd < 0) return -1;
    if (fstat(m->fd, &st) || st.st_size <= 0) { close(m->fd); m->fd = -1; return -1; }
    m->bytes = (uint64_t)st.st_size;
#endif
    return 0;
}

static const unsigned char *cs_view(CsMap *m, uint64_t offset, uint64_t length)
{
    uint64_t gran, aligned, skip;
    cs_unview(m);
    if (!length || offset > m->bytes || length > m->bytes-offset) return NULL;
#ifdef _WIN32
    SYSTEM_INFO sys;
    GetSystemInfo(&sys); gran = sys.dwAllocationGranularity;
#else
    gran = (uint64_t)sysconf(_SC_PAGE_SIZE);
#endif
    aligned = offset/gran*gran; skip = offset-aligned;
    if (length > SIZE_MAX-skip) return NULL;
    m->length = (size_t)(length+skip);
#ifdef _WIN32
    m->base = MapViewOfFile(m->mapping, FILE_MAP_READ, (DWORD)(aligned>>32), (DWORD)aligned, m->length);
#else
    m->base = mmap(NULL, m->length, PROT_READ, MAP_PRIVATE, m->fd, (off_t)aligned);
    if (m->base == MAP_FAILED) m->base = NULL;
#endif
    return m->base ? (const unsigned char *)m->base+(size_t)skip : NULL;
}

static void cs_close(CsMap *m)
{
    cs_unview(m);
#ifdef _WIN32
    if (m->mapping) CloseHandle(m->mapping);
    if (m->file) CloseHandle(m->file);
#else
    if (m->bytes) close(m->fd);
#endif
    memset(m, 0, sizeof *m);
}

static int cs_seek(FILE *f, uint64_t off)
{
#ifdef _WIN32
    return off > INT64_MAX ? -1 : _fseeki64(f, (__int64)off, SEEK_SET);
#else
    return off > INT64_MAX ? -1 : fseeko(f, (off_t)off, SEEK_SET);
#endif
}

/* a regular file at path: a directory where the index cache goes (a lost race's leftover) is no cache -- Linux's fopen
 * opens a directory and every read of it then fails, which refused the seed; Windows's fopen fails on it */
static int cs_is_file(const char *path)
{
#ifdef _WIN32
    struct _stat64 st;
    return _stat64(path, &st) == 0 && (st.st_mode & _S_IFREG) != 0;
#else
    struct stat st;
    return stat(path, &st) == 0 && S_ISREG(st.st_mode);
#endif
}

static int cs_replace(const char *from, const char *to)
{
#ifdef _WIN32
    return MoveFileExA(from, to, MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH) ? 0 : -1;
#else
    return rename(from, to);
#endif
}

static int cs_bytes(const char *path, const void *data, size_t n)
{
    char tmp[4096];
    if (snprintf(tmp, sizeof tmp, "%s.%d.tmp", path, ves_getpid()) >= (int)sizeof tmp) return -1;
    FILE *f = fopen(tmp, "wb");
    if (!f) return -1;
    int bad = n && fwrite(data, 1, n, f) != n;
    if (fclose(f)) bad = 1;
    return bad ? -1 : cs_replace(tmp, path);
}

static int cs_copy(const char *from, const char *to)
{
    CsMap m;
    if (cs_open(&m, from)) return -1;
    const void *data = cs_view(&m, 0, m.bytes);
    int rc = data ? cs_bytes(to, data, (size_t)m.bytes) : -1;
    cs_close(&m);
    return rc;
}

static int cs_cmp(const void *a, const void *b)
{
    int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
    return (x > y)-(x < y);
}

static size_t cs_find(const int32_t *ids, size_t n, int32_t id)
{
    size_t lo = 0, hi = n;
    while (lo < hi) { size_t mid = lo+(hi-lo)/2; if (ids[mid] < id) lo = mid+1; else hi = mid; }
    return lo;
}

static double cs_uv(const unsigned char *uv, int f64, size_t k)
{
    return f64 ? ((const double *)uv)[k] : (double)((const float *)uv)[k];
}

static int cs_meets(const CsBounds *b, const double box[4])
{
    return b->hi[0] >= box[0] && b->lo[0] < box[2] && b->hi[1] >= box[1] && b->lo[1] < box[3];
}

static int cs_count(const JsonValue *j, const char *key, uint64_t limit, uint64_t *out)
{
    double value = Json_member_double(j, key, -1);
    if (!isfinite(value) || value < 0 || value > (double)limit || value != floor(value)) return -1;
    *out = (uint64_t)value;
    return 0;
}

int CloseupSeed_write(const char *manifest)
{
    Arena_T arena = Arena_new();
    const char *err = NULL;
    const JsonValue *j = Json_parse_file(arena, manifest, &err);
    const char *mesh = Json_as_string(Json_object_get(j, "mesh"));
    const char *index = Json_as_string(Json_object_get(j, "index"));
    const char *out = Json_as_string(Json_object_get(j, "out"));
    const char *cache = Json_as_string(Json_object_get(j, "cache"));
    const char *ribbon = Json_as_string(Json_object_get(j, "ribbon"));
    const JsonValue *parts = Json_object_get(j, "parts"), *jb = Json_object_get(j, "box");
    MeshBinHeader h;
    CsMap mm = {0}, cm = {0};
    CsBounds *bounds = NULL;
    int32_t *batch = NULL, *faces = NULL, *ids = NULL, *charts = NULL, *allowed = NULL;
    float *xyz = NULL, *uv32 = NULL;
    double *uv64 = NULL;
    uint64_t *counts = NULL;
    FILE *fm = NULL, *fi = NULL;
    size_t nf = 0, cap = 0, nv = 0, na = 0, np = Json_array_len(parts), nb = 0;
    int rc = -1, cached = 0, f64 = 0;
    double box[4], margin = Json_member_double(j, "grow_px", -1), begun = ves_clock_sec();
    char path[4096], temp[4096];
    CsIndex ih;
    if (!mesh || !index || !out || !cache || !np || Json_array_len(jb) != 4 || !isfinite(margin) || margin < 0 ||
        MeshBin_header(mesh, &h) || !(h.flags&MESH_BIN_HAS_UV) || h.nv > INT32_MAX || !h.nf) goto done;
    for (size_t k = 0; k < 4; k++) {
        box[k] = Json_as_double(Json_array_get(jb, k), NAN);
        if (!isfinite(box[k])) goto done;
    }
    if (box[0] >= box[2] || box[1] >= box[3]) goto done;
    box[0] -= margin; box[1] -= margin; box[2] += margin; box[3] += margin;
    counts = calloc(np, sizeof *counts);
    uint64_t at = 0, source_faces = 0;
    for (size_t p = 0; p < np; p++) {
        const JsonValue *part = Json_array_get(parts, p);
        uint64_t a, b, nfaces;
        const char *chart = Json_as_string(Json_object_get(part, "chart"));
        if (cs_count(part, "start", h.nv, &a) || cs_count(part, "stop", h.nv, &b) ||
            cs_count(part, "faces", h.nf-source_faces, &nfaces) || a != at || b < a || !chart) goto done;
        if (b > a) {
            if (cs_open(&cm, chart) || cm.bytes != (b-a)*4) goto done;
            cs_close(&cm);
        }
        at = b; source_faces += nfaces;
    }
    if (!counts || at != h.nv || source_faces != h.nf) goto done;
    const JsonValue *list = Json_object_get(Json_parse_file(arena, index, &err), "charts");
    if (Json_type(list) != JSON_ARRAY) goto done;
    allowed = malloc((Json_array_len(list)+1)*sizeof *allowed);
    if (!allowed) goto done;
    for (size_t k = 0; k < Json_array_len(list); k++) {
        const JsonValue *ch = Json_array_get(list, k);
        long id = Json_member_long(ch, "id", -1);
        if (id >= 0 && Json_as_bool(Json_object_get(ch, "placed"), 0)) allowed[na++] = (int32_t)id;
    }
    qsort(allowed, na, sizeof *allowed, cs_cmp);
    nb = (size_t)((h.nf+CS_BLOCK-1)/CS_BLOCK);
    bounds = malloc(nb*sizeof *bounds); batch = malloc((size_t)CS_BATCH*3*sizeof *batch);
    if (!bounds || !batch || cs_open(&mm, mesh)) goto done;
    f64 = (h.flags&MESH_BIN_UV_F64) != 0;
    memcpy(ih.magic, "VCSIDX01", 8); ih.nv = h.nv; ih.nf = h.nf; ih.bytes = h.file_bytes; ih.block = CS_BLOCK;
    fi = cs_is_file(cache) ? fopen(cache, "rb") : NULL;
    if (fi) {
        CsIndex got;
        if (fread(&got, sizeof got, 1, fi) != 1 || memcmp(&got, &ih, sizeof ih) ||
            fread(bounds, sizeof *bounds, nb, fi) != nb || fgetc(fi) != EOF) goto done;
        fclose(fi); fi = NULL; cached = 1;
    }
    fm = fopen(mesh, "rb");
    if (!fm) goto done;
    size_t reads = 0;
    for (uint64_t start = 0; start < h.nf;) {
        if (cached && !cs_meets(&bounds[(size_t)(start/CS_BLOCK)], box)) { start += CS_BLOCK; continue; }
        uint64_t end = start+CS_BLOCK;
        while (end < h.nf && end-start < CS_BATCH && (!cached || cs_meets(&bounds[(size_t)(end/CS_BLOCK)], box)))
            end += CS_BLOCK;
        if (end > h.nf) end = h.nf;
        size_t n = (size_t)(end-start);
        if (cs_seek(fm, h.faces_offset+start*12) || fread(batch, 12, n, fm) != n) goto done;
        int32_t vmin = INT32_MAX, vmax = -1;
        for (size_t t = 0; t < n*3; t++) {
            if (batch[t] < 0 || (uint64_t)batch[t] >= h.nv) goto done;
            if (batch[t] < vmin) vmin = batch[t];
            if (batch[t] > vmax) vmax = batch[t];
        }
        uint64_t stride = f64 ? 16 : 8;
        const unsigned char *uv = cs_view(&mm, h.uv_offset+(uint64_t)vmin*stride, ((uint64_t)(vmax-vmin)+1)*stride);
        if (!uv) goto done;
        for (size_t t = 0; t < n; t++) {
            CsBounds *b = &bounds[(size_t)((start+t)/CS_BLOCK)];
            if (!cached && (start+t)%CS_BLOCK == 0) {
                b->lo[0] = b->lo[1] = INFINITY; b->hi[0] = b->hi[1] = -INFINITY;
            }
            double u = 0, v = 0;
            for (size_t k = 0; k < 3; k++) {
                size_t q = (size_t)(batch[t*3+k]-vmin)*2;
                u += cs_uv(uv, f64, q)/3; v += cs_uv(uv, f64, q+1)/3;
            }
            if (!isfinite(u) || !isfinite(v)) continue;
            if (!cached) {
                b->lo[0] = fmin(b->lo[0], u); b->lo[1] = fmin(b->lo[1], v);
                b->hi[0] = fmax(b->hi[0], u); b->hi[1] = fmax(b->hi[1], v);
            }
            if (u < box[0] || u >= box[2] || v < box[1] || v >= box[3]) continue;
            if (nf == cap) {
                size_t nc = cap ? cap*2 : 65536;
                if (nc < cap || nc > SIZE_MAX/12) goto done;
                int32_t *next = realloc(faces, nc*12);
                if (!next) goto done;
                faces = next; cap = nc;
            }
            memcpy(faces+nf*3, batch+t*3, 12); nf++;
        }
        cs_unview(&mm); reads++;
        if (reads%64 == 0) fprintf(stderr, "[closeup seed] scanned through %llu/%llu faces; selected %zu\n",
                                   (unsigned long long)end, (unsigned long long)h.nf, nf);
        start = end;
    }
    fclose(fm); fm = NULL;
    if (!cached) {
        if (snprintf(temp, sizeof temp, "%s.%d.tmp", cache, ves_getpid()) >= (int)sizeof temp) goto done;
        fi = fopen(temp, "wb");
        if (!fi || fwrite(&ih, sizeof ih, 1, fi) != 1 || fwrite(bounds, sizeof *bounds, nb, fi) != nb) goto done;
        if (fclose(fi)) { fi = NULL; goto done; } fi = NULL;
        /* (2026-10-24) The index is a shared CACHE, not a result: simultaneous close-ups of one sheet all build it (five
         * did at 20:02Z on the PHerc0358 21^3) and on Windows the replace of a file another process holds open fails.
         * A lost race drops the temp and keeps the seed; it refused close-up cu20261024T195830Z_00601c after a
         * complete selection. */
        if (cs_replace(temp, cache)) {
            remove(temp);
            fprintf(stderr, "[closeup seed] index cache not replaced (another close-up holds %s): the seed goes on without it\n", cache);
        }
    }
    if (!nf || nf > SIZE_MAX/12) goto done;
    ids = malloc(nf*12);
    if (!ids) goto done;
    memcpy(ids, faces, nf*12); qsort(ids, nf*3, sizeof *ids, cs_cmp);
    for (size_t k = 0; k < nf*3; k++) if (!nv || ids[k] != ids[nv-1]) ids[nv++] = ids[k];
    xyz = malloc(nv*12); charts = malloc(nv*4);
    if (f64) uv64 = malloc(nv*16); else uv32 = malloc(nv*8);
    if (!xyz || !charts || (f64 ? uv64 == NULL : uv32 == NULL)) goto done;
    const unsigned char *all = cs_view(&mm, 0, mm.bytes);
    if (!all) goto done;
    for (size_t k = 0; k < nv; k++) {
        memcpy(xyz+k*3, all+h.vertices_offset+(uint64_t)ids[k]*12, 12);
        if (f64) memcpy(uv64+k*2, all+h.uv_offset+(uint64_t)ids[k]*16, 16);
        else memcpy(uv32+k*2, all+h.uv_offset+(uint64_t)ids[k]*8, 8);
    }
    cs_close(&mm);
    size_t used = 0;
    for (size_t p = 0; p < np; p++) {
        const JsonValue *part = Json_array_get(parts, p);
        uint64_t a = (uint64_t)Json_member_double(part, "start", -1), b = (uint64_t)Json_member_double(part, "stop", -1);
        if (used == nv || (uint64_t)ids[used] >= b) continue;
        if (cs_open(&cm, Json_as_string(Json_object_get(part, "chart")))) goto done;
        const int32_t *owners = (const int32_t *)cs_view(&cm, 0, cm.bytes);
        if (!owners) goto done;
        while (used < nv && (uint64_t)ids[used] < b) {
            int32_t owner = owners[(size_t)((uint64_t)ids[used]-a)];
            size_t k = cs_find(allowed, na, owner);
            if (k == na || allowed[k] != owner) { fprintf(stderr, "ERROR: unregistered seed owner %d\n", owner); goto done; }
            charts[used++] = owner; counts[p]++;
        }
        cs_close(&cm);
    }
    if (used != nv) goto done;
    for (size_t k = 0; k < nf*3; k++) faces[k] = (int32_t)cs_find(ids, nv, faces[k]);
    if (snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", out) >= (int)sizeof path) goto done;
    if (f64 ? MeshBin_write_uv64(path, xyz, nv, faces, nf, uv64) : MeshBin_write(path, xyz, nv, faces, nf, uv32)) goto done;
    snprintf(path, sizeof path, "%s/stage6_repaired_chart.i32", out);
    if (cs_bytes(path, charts, nv*4)) goto done;
    snprintf(path, sizeof path, "%s/sheet_index.json", out);
    if (cs_copy(index, path)) goto done;
    if (ribbon) { snprintf(path, sizeof path, "%s/ribbon_layout.json", out); if (cs_copy(ribbon, path)) goto done; }
    snprintf(path, sizeof path, "%s/seed_native.json", out);
    fi = fopen(path, "wb");
    if (!fi) goto done;
    fprintf(fi, "{\"vertices\":%zu,\"faces\":%zu,\"indexed\":%s,\"face_batches_read\":%zu,\"vertices_by_part\":[",
            nv, nf, cached ? "true" : "false", reads);
    for (size_t p = 0; p < np; p++) fprintf(fi, "%s%llu", p ? "," : "", (unsigned long long)counts[p]);
    fprintf(fi, "],\"seconds\":%.3f}\n", ves_clock_sec()-begun);
    if (fclose(fi)) { fi = NULL; goto done; } fi = NULL;
    fprintf(stderr, "[closeup seed] %zu vertices, %zu faces, %zu batches (%s index), %.1f s\n", nv, nf, reads,
            cached ? "reused" : "built", ves_clock_sec()-begun);
    rc = 0;
done:
    if (rc) fprintf(stderr, "ERROR: rendered close-up seed refused (%s): %s\n", manifest, err ? err : "invalid input or I/O");
    if (fm) fclose(fm);
    if (fi) fclose(fi);
    cs_close(&mm); cs_close(&cm);
    free(bounds); free(batch); free(faces); free(ids); free(charts); free(allowed);
    free(xyz); free(uv32); free(uv64); free(counts); Arena_dispose(&arena);
    return rc;
}

int CloseupSeed_selftest(void)
{
    char dir[1024], mesh[1200], chart[1200], index[1200], out[1200], cfg[1200], cache[1200], path[1400];
    int fails = 0;
    if (ves_selftest_dir(dir, sizeof dir, "closeup_seed")) return 1;
    snprintf(mesh, sizeof mesh, "%s/input.vmesh", dir); snprintf(chart, sizeof chart, "%s/chart.i32", dir);
    snprintf(index, sizeof index, "%s/index.json", dir); snprintf(cfg, sizeof cfg, "%s/crop.json", dir);
    snprintf(out, sizeof out, "%s/out", dir); snprintf(cache, sizeof cache, "%s/seed.idx", dir); ves_mkdir(out);
    float xyz[27]; double uv[18]; int32_t faces[9], owners[9];
    for (int k = 0; k < 9; k++) {
        xyz[k*3] = (float)k; xyz[k*3+1] = 3; xyz[k*3+2] = 5;
        uv[k*2] = 1000000.001+(k/3 == 2 ? 20 : k/3*2)+(k%3 == 1 ? 1 : 0);
        uv[k*2+1] = k%3 == 2 ? 1 : 0; faces[k] = k; owners[k] = k/3 == 1 ? 9 : 2;
    }
    const char *js = "{\"charts\":[{\"id\":2,\"region\":7,\"placed\":true},{\"id\":9,\"region\":7,\"placed\":true}]}";
    if (MeshBin_write_uv64(mesh, xyz, 9, faces, 3, uv) || cs_bytes(chart, owners, sizeof owners) || cs_bytes(index, js, strlen(js))) fails++;
    FILE *f = fopen(cfg, "wb");
    if (!f) fails++;
    else {
        fprintf(f, "{\"mesh\":\"%s\",\"index\":\"%s\",\"out\":\"%s\",\"cache\":\"%s\","
                   "\"box\":[1000000,0,1000004,2],\"grow_px\":0,\"parts\":[{\"chart\":\"%s\",\"start\":0,\"stop\":9,\"faces\":3}]}",
                mesh,index,out,cache,chart); fclose(f);
    }
    for (int pass = 0; !fails && pass < 2; pass++) {
        MeshBinData got = {0};
        snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", out);
        if (CloseupSeed_write(cfg) || MeshBin_read_precise_malloc(path, &got)) fails++;
        else if (got.nv != 6 || got.nf != 2 || memcmp(got.verts, xyz, 6*12) || memcmp(got.uv64, uv, 6*16) ||
                 memcmp(got.faces, faces, 2*12)) fails++;
        MeshBin_dispose(&got);
        snprintf(path, sizeof path, "%s/stage6_repaired_chart.i32", out);
        int32_t read_owners[6]; f = fopen(path, "rb");
        if (!f || fread(read_owners, 4, 6, f) != 6 || memcmp(read_owners, owners, sizeof read_owners)) fails++;
        if (f) fclose(f);
    }
    /* (2026-10-24) a lost index-cache race must not refuse the seed: a directory sits where the cache goes, so the
     * replace fails; the seed is still written and the temp does not linger */
    if (!fails) {
        char cfg2[1200], blocked[1200], temp[1400], out2[1200];
        snprintf(blocked, sizeof blocked, "%s/blocked.idx", dir); ves_mkdir(blocked);
        snprintf(out2, sizeof out2, "%s/out_blocked", dir); ves_mkdir(out2);
        snprintf(cfg2, sizeof cfg2, "%s/crop_blocked.json", dir);
        f = fopen(cfg2, "wb");
        if (!f) fails++;
        else {
            fprintf(f, "{\"mesh\":\"%s\",\"index\":\"%s\",\"out\":\"%s\",\"cache\":\"%s\","
                       "\"box\":[1000000,0,1000004,2],\"grow_px\":0,\"parts\":[{\"chart\":\"%s\",\"start\":0,\"stop\":9,\"faces\":3}]}",
                    mesh,index,out2,blocked,chart); fclose(f);
            MeshBinData got = {0};
            snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", out2);
            if (CloseupSeed_write(cfg2) || MeshBin_read_precise_malloc(path, &got)) fails++;
            else if (got.nv != 6 || got.nf != 2) fails++;
            MeshBin_dispose(&got);
            snprintf(temp, sizeof temp, "%s.%d.tmp", blocked, ves_getpid());
            f = fopen(temp, "rb");
            if (f) { fclose(f); fails++; }
            fprintf(stderr, "[closeup seed] selftest: a lost index-cache replace keeps the seed: %s\n", fails ? "FAIL" : "ok");
        }
    }
    owners[4] = 77; cs_bytes(chart, owners, sizeof owners);
    if (!CloseupSeed_write(cfg)) fails++;
    fprintf(stderr, "[closeup seed] selftest %s\n", fails ? "FAIL" : "ok");
    if (!fails) ves_selftest_clean(dir);
    return fails ? 1 : 0;
}
