/* chart_pool.c -- an assembler chart store (ASMCHT02) as an ink_closeup adoption pool. */
#include "chart_pool.h"

#include "../assemble/asm_types.h"
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_chart_index.h"
#include "../common/ves_platform.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#endif

/* The chart store's records, as src/assemble/asm_store.c writes them (fwrite of its private structs, little-endian,
 * MSVC/GCC layout on x86_64): an 8-byte magic, the identity of the source mesh, its path, a chart count, then per
 * chart a fixed header followed by vid[nv] int32, xyz[nv*3] float, nrm[nv*3] float, uv[nv*2] float (has_uv),
 * faces[nf*3] int32, boundary[nv] uint8, repaired[nv] uint8 (has_repaired).  The sizes are checked below: a change of
 * the store's layout fails the build, never the read. */
#define CP_MAGIC "ASMCHT02"
enum { CP_MAX_NAME = 2048, CP_MAX_ELEMS = 1 << 28, CP_MAX_CHARTS = 1 << 26 };

typedef struct CpIdentity {
    uint64_t fingerprint;
    uint64_t src_size;
    int64_t  src_mtime;
    uint32_t path_len;
    uint32_t reserved;
} CpIdentity;

typedef struct CpHeader {
    int32_t  id, cube, comp, parent;
    uint32_t flags;
    int32_t  component, placed;
    uint32_t has_uv, has_repaired;
    uint64_t nv, nf, n_flipped;
    double   area3d, area_uv, sigma_lo, sigma_hi, stress_frac;
    double   centroid[3];
    double   pose_x, pose_y, pose_theta;
    float    bbox_lo[3], bbox_hi[3];
} CpHeader;

typedef char cp_identity_size_check[(sizeof(CpIdentity) == 32) ? 1 : -1];
typedef char cp_header_size_check[(sizeof(CpHeader) == 176) ? 1 : -1];

typedef struct CpCube {
    long oz, oy, ox;
    char name[64];           /* <cube id>.asc */
} CpCube;

/* "z#####_y#####_x#####.asc" -> origin; 0 on a match */
static int cp_parse_name(const char *name, long *oz, long *oy, long *ox)
{
    char tail[8];
    int got = 0;
    size_t n = strlen(name);
    if (n < 24 || strcmp(name + n - 4, ".asc") != 0) return -1;
    got = sscanf(name, "z%5ld_y%5ld_x%5ld%7s", oz, oy, ox, tail);
    return got == 4 && strcmp(tail, ".asc") == 0 ? 0 : -1;
}

static int cp_cmp_cube(const void *a, const void *b)
{
    const CpCube *p = (const CpCube *)a, *q = (const CpCube *)b;
    if (p->oz != q->oz) return p->oz < q->oz ? -1 : 1;
    if (p->oy != q->oy) return p->oy < q->oy ? -1 : 1;
    if (p->ox != q->ox) return p->ox < q->ox ? -1 : 1;
    return 0;
}

static int cp_push(CpCube **list, size_t *n, size_t *cap, const char *name, const long *box)
{
    long oz = 0, oy = 0, ox = 0;
    if (cp_parse_name(name, &oz, &oy, &ox) != 0) return 0;
    if (box != NULL && !(oz >= box[0] && oz < box[1] && oy >= box[2] && oy < box[3] && ox >= box[4] && ox < box[5]))
        return 0;
    if (strlen(name) >= sizeof (*list)[0].name) return 0;
    if (*n == *cap) {
        size_t nc = *cap ? *cap * 2 : 1024;
        CpCube *grown = (CpCube *)realloc(*list, nc * sizeof *grown);
        if (grown == NULL) return -1;
        *list = grown;
        *cap = nc;
    }
    (*list)[*n].oz = oz; (*list)[*n].oy = oy; (*list)[*n].ox = ox;
    memcpy((*list)[*n].name, name, strlen(name) + 1);
    (*n)++;
    return 0;
}

/* The store's cube files in the box, sorted by origin (malloc'd).  0, -1 when the directory cannot be listed. */
static int cp_list(const char *dir, const long *box, CpCube **out, size_t *n_out)
{
    CpCube *list = NULL;
    size_t n = 0, cap = 0;
    int rc = 0;
#ifdef _WIN32
    char glob[CP_MAX_NAME];
    WIN32_FIND_DATAA fd;
    HANDLE h;
    snprintf(glob, sizeof glob, "%s/*.asc", dir);
    h = FindFirstFileA(glob, &fd);
    if (h != INVALID_HANDLE_VALUE) {
        do {
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && cp_push(&list, &n, &cap, fd.cFileName, box) != 0)
                rc = -1;
        } while (rc == 0 && FindNextFileA(h, &fd));
        FindClose(h);
    } else if (GetLastError() != ERROR_FILE_NOT_FOUND) {
        rc = -1;
    }
#else
    DIR *d = opendir(dir);
    struct dirent *e = NULL;
    if (d == NULL) return -1;
    while (rc == 0 && (e = readdir(d)) != NULL) rc = cp_push(&list, &n, &cap, e->d_name, box);
    closedir(d);
#endif
    if (rc != 0) { free(list); return -1; }
    if (n > 1) qsort(list, n, sizeof *list, cp_cmp_cube);
    *out = list;
    *n_out = n;
    return 0;
}

static int cp_get(FILE *fp, void *p, size_t n)
{
    return n == 0 ? 0 : (fread(p, 1, n, fp) == n ? 0 : -1);
}

static int cp_skip(FILE *fp, uint64_t n)
{
#ifdef _WIN32
    return n == 0 ? 0 : (_fseeki64(fp, (__int64)n, SEEK_CUR) == 0 ? 0 : -1);
#else
    return n == 0 ? 0 : (fseeko(fp, (off_t)n, SEEK_CUR) == 0 ? 0 : -1);
#endif
}

/* Open a cube file at its first chart: the chart count, or -1 when it is not a chart store. */
static long cp_open(const char *path, FILE **fpo)
{
    FILE *fp = fopen(path, "rb");
    char magic[8];
    CpIdentity id;
    uint64_t count = 0;
    *fpo = NULL;
    if (fp == NULL) return -1;
    if (cp_get(fp, magic, 8) != 0 || memcmp(magic, CP_MAGIC, 8) != 0 || cp_get(fp, &id, sizeof id) != 0 ||
        id.path_len > 4096 || cp_skip(fp, id.path_len) != 0 || cp_get(fp, &count, sizeof count) != 0 ||
        count > (uint64_t)CP_MAX_CHARTS) {
        fclose(fp);
        return -1;
    }
    *fpo = fp;
    return (long)count;
}

static int cp_header_ok(const CpHeader *h)
{
    return h->nv <= (uint64_t)CP_MAX_ELEMS && h->nf <= (uint64_t)CP_MAX_ELEMS;
}

/* A chart the assembler lays out: none of the excluding flags, a flattening, at least one face. */
static int cp_keep(const CpHeader *h)
{
    return (h->flags & ASM_CHART_EXCLUDED_MASK) == 0u && h->has_uv != 0 && h->nv >= 3 && h->nf >= 1;
}

/* The bytes after the header of a chart record. */
static uint64_t cp_payload(const CpHeader *h)
{
    uint64_t nv = h->nv, nf = h->nf;
    return nv * 4u + nv * 12u + nv * 12u + (h->has_uv ? nv * 8u : 0u) + nf * 12u + nv + (h->has_repaired ? nv : 0u);
}

/* Pass 1 over one cube: the kept charts' vertex and face totals.  0, -1 when the file is malformed. */
static int cp_count(const char *path, size_t *charts, size_t *skipped, uint64_t *nv, uint64_t *nf)
{
    FILE *fp = NULL;
    long count = cp_open(path, &fp), k = 0;
    size_t c = 0, s = 0;
    uint64_t v = 0, f = 0;
    if (count < 0) return -1;
    for (k = 0; k < count; k++) {
        CpHeader h;
        if (cp_get(fp, &h, sizeof h) != 0 || !cp_header_ok(&h) || cp_skip(fp, cp_payload(&h)) != 0) {
            fclose(fp);
            return -1;
        }
        if (cp_keep(&h)) { c++; v += h.nv; f += h.nf; } else { s++; }
    }
    fclose(fp);
    *charts += c; *skipped += s; *nv += v; *nf += f;
    return 0;
}

/* Pass 2 over one cube: its kept charts appended to the pool (chart ids from *next_id).  0, -1 on error. */
static int cp_append(Arena_T arena, const char *path, MeshBinStream_T writer, FILE *chart_fp, int32_t *next_id,
                     size_t *nv_out, size_t *nf_out)
{
    FILE *fp = NULL;
    long count = cp_open(path, &fp), k = 0;
    int rc = 0;
    if (count < 0) return -1;
    for (k = 0; k < count && rc == 0; k++) {
        CpHeader h;
        Arena_Mark mark = Arena_save(arena);
        size_t nv = 0, nf = 0, i = 0;
        if (cp_get(fp, &h, sizeof h) != 0 || !cp_header_ok(&h)) { rc = -1; break; }
        nv = (size_t)h.nv; nf = (size_t)h.nf;
        if (!cp_keep(&h)) {
            rc = cp_skip(fp, cp_payload(&h));
            continue;
        }
        {
            float *xyz = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof *xyz);
            float *uvf = (float *)ARENA_ALLOC(arena, nv * 2 * sizeof *uvf);
            double *uv = (double *)ARENA_ALLOC(arena, nv * 2 * sizeof *uv);
            int32_t *faces = (int32_t *)ARENA_ALLOC(arena, nf * 3 * sizeof *faces);
            int32_t *ids = (int32_t *)ARENA_ALLOC(arena, nv * sizeof *ids);
            MeshBinData chunk;
            rc = cp_skip(fp, (uint64_t)nv * 4u);                                         /* vid */
            if (rc == 0) rc = cp_get(fp, xyz, nv * 3 * sizeof *xyz);
            if (rc == 0) rc = cp_skip(fp, (uint64_t)nv * 12u);                           /* nrm */
            if (rc == 0) rc = cp_get(fp, uvf, nv * 2 * sizeof *uvf);
            if (rc == 0) rc = cp_get(fp, faces, nf * 3 * sizeof *faces);
            if (rc == 0) rc = cp_skip(fp, (uint64_t)nv + (h.has_repaired ? (uint64_t)nv : 0u));   /* boundary, repaired */
            for (i = 0; rc == 0 && i < nf * 3; i++) if (faces[i] < 0 || (size_t)faces[i] >= nv) rc = -1;
            for (i = 0; rc == 0 && i < nv * 2; i++) {
                if (!isfinite(uvf[i])) rc = -1;
                uv[i] = (double)uvf[i];
            }
            for (i = 0; rc == 0 && i < nv * 3; i++) if (!isfinite(xyz[i])) rc = -1;
            if (rc == 0) {
                memset(&chunk, 0, sizeof chunk);
                chunk.verts = xyz; chunk.faces = faces; chunk.uv64 = uv; chunk.nv = nv; chunk.nf = nf;
                rc = MeshBin_stream_append(writer, &chunk);
            }
            for (i = 0; rc == 0 && i < nv; i++) ids[i] = *next_id;
            if (rc == 0 && fwrite(ids, sizeof *ids, nv, chart_fp) != nv) rc = -1;
            if (rc == 0) { (*next_id)++; *nv_out += nv; *nf_out += nf; }
        }
        Arena_restore(arena, mark);
    }
    fclose(fp);
    return rc;
}

static void cp_join(char *out, size_t cap, const char *dir, const char *name)
{
    snprintf(out, cap, "%s/%s", dir, name);
}

int ChartPool_from_store(const char *store_dir, const char *out_dir, const long *box, ChartPoolStats *st)
{
    char path[CP_MAX_NAME], vmesh[CP_MAX_NAME], chart[CP_MAX_NAME], ctmp[CP_MAX_NAME], idx[CP_MAX_NAME];
    CpCube *cubes = NULL;
    uint8_t *bad = NULL;
    size_t n = 0, k = 0, charts = 0, skipped = 0, nv_done = 0, nf_done = 0;
    uint64_t nv = 0, nf = 0;
    int32_t next_id = 0;
    Arena_T arena = NULL;
    MeshBinStream_T writer = NULL;
    FILE *cfp = NULL, *rep = NULL;
    double t0 = ves_clock_sec();
    int rc = -1;
    ChartPoolStats s0;
    if (st == NULL) st = &s0;
    memset(st, 0, sizeof *st);
    if (store_dir == NULL || out_dir == NULL) return -1;
    if (cp_list(store_dir, box, &cubes, &n) != 0) {
        fprintf(stderr, "ERROR: cannot list the chart store %s\n", store_dir);
        return -1;
    }
    st->cubes = n;
    bad = (uint8_t *)calloc(n ? n : 1, 1);
    if (bad == NULL) goto out;
    for (k = 0; k < n; k++) {                               /* pass 1: what the pool will hold */
        cp_join(path, sizeof path, store_dir, cubes[k].name);
        if (cp_count(path, &charts, &skipped, &nv, &nf) != 0) {
            bad[k] = 1;
            st->cubes_failed++;
            fprintf(stderr, "  [chart_pool] %s: not a readable chart store; skipped\n", cubes[k].name);
        }
    }
    st->charts_skipped = skipped;
    if (charts == 0 || nv == 0) { rc = 1; goto out; }
    if (nv > (uint64_t)INT32_MAX || charts > (size_t)INT32_MAX) {
        fprintf(stderr, "ERROR: the pool would hold %llu vertices: past int32 face indices\n", (unsigned long long)nv);
        goto out;
    }
    (void)ves_mkdir(out_dir);
    cp_join(vmesh, sizeof vmesh, out_dir, "stage6_repaired.vmesh");
    cp_join(chart, sizeof chart, out_dir, "stage6_repaired_chart.i32");
    snprintf(ctmp, sizeof ctmp, "%s.partial", chart);
    arena = Arena_new();
    cfp = fopen(ctmp, "wb");
    if (cfp == NULL || MeshBin_stream_open_uv64(arena, vmesh, (size_t)nv, (size_t)nf, &writer) != 0) {
        fprintf(stderr, "ERROR: cannot write the pool in %s\n", out_dir);
        goto out;
    }
    for (k = 0; k < n; k++) {                               /* pass 2: the charts, in cube order */
        if (bad[k]) continue;
        cp_join(path, sizeof path, store_dir, cubes[k].name);
        if (cp_append(arena, path, writer, cfp, &next_id, &nv_done, &nf_done) != 0) {
            fprintf(stderr, "ERROR: %s changed or broke between the passes\n", path);
            goto out;
        }
    }
    if ((uint64_t)nv_done != nv || (uint64_t)nf_done != nf) goto out;
    if (fclose(cfp) != 0) { cfp = NULL; goto out; }
    cfp = NULL;
    if (MeshBin_stream_close(&writer, 1) != 0) goto out;
    remove(chart);
    if (rename(ctmp, chart) != 0) goto out;
    if (MeshChartIndex_path(vmesh, idx, sizeof idx) != 0 || MeshChartIndex_build(vmesh, chart, idx) != 0) {
        fprintf(stderr, "ERROR: cannot index the pool %s\n", vmesh);
        goto out;
    }
    st->charts = (size_t)next_id;
    st->verts = nv_done;
    st->faces = nf_done;
    st->seconds = ves_clock_sec() - t0;
    cp_join(path, sizeof path, out_dir, "chart_pool_report.json");
    rep = fopen(path, "wb");
    if (rep != NULL) {
        fprintf(rep, "{\n  \"schema\": \"chart-pool-v1\",\n  \"cubes\": %zu,\n  \"cubes_failed\": %zu,\n  \"charts\": %zu,\n"
                "  \"charts_skipped\": %zu,\n  \"verts\": %zu,\n  \"faces\": %zu,\n  \"seconds\": %.1f,\n",
                st->cubes, st->cubes_failed, st->charts, st->charts_skipped, st->verts, st->faces, st->seconds);
        if (box != NULL)
            fprintf(rep, "  \"box_zyx\": [%ld, %ld, %ld, %ld, %ld, %ld],\n", box[0], box[1], box[2], box[3], box[4], box[5]);
        fprintf(rep, "  \"rule\": \"charts with no BLOB/EXTRAS/NONMANIFOLD/FLAT_FAILED/TINY flag and an intrinsic UV; "
                "ids 0..n-1 in cube order; UV = the chart's own flattening (float64)\"\n}\n");
        fclose(rep);
    }
    fprintf(stderr, "[chart_pool] %s: %zu charts (%zu skipped) of %zu cubes (%zu unreadable), %zu verts, %zu faces "
            "(%.1fs)\n", out_dir, st->charts, st->charts_skipped, st->cubes, st->cubes_failed, st->verts, st->faces,
            st->seconds);
    rc = 0;
out:
    if (writer != NULL) MeshBin_stream_close(&writer, 0);
    if (cfp != NULL) { fclose(cfp); remove(ctmp); }
    if (arena != NULL) Arena_dispose(&arena);
    free(cubes);
    free(bad);
    return rc;
}

/* ------------------------------------------------------------------ selftest */

/* One cube file in the store's own layout: count charts of 3 vertices / 1 face; chart k gets flags[k], UV unless
 * no_uv[k], vertices offset by k. */
static int cp_test_write(const char *path, int count, const uint32_t *flags, const int *no_uv, float base)
{
    FILE *fp = fopen(path, "wb");
    CpIdentity id;
    uint64_t c = (uint64_t)count;
    const char *src = "some/mesh.vmesh";
    int k = 0, ok = 1;
    if (fp == NULL) return -1;
    memset(&id, 0, sizeof id);
    id.path_len = (uint32_t)strlen(src);
    ok &= fwrite(CP_MAGIC, 1, 8, fp) == 8 && fwrite(&id, sizeof id, 1, fp) == 1 &&
          fwrite(src, 1, id.path_len, fp) == id.path_len && fwrite(&c, sizeof c, 1, fp) == 1;
    for (k = 0; k < count && ok; k++) {
        CpHeader h;
        int32_t vid[3] = { 0, 1, 2 }, faces[3] = { 0, 1, 2 };
        float xyz[9] = { 0 }, nrm[9] = { 0 }, uv[6] = { 0.f, 0.f, 1.f, 0.f, 0.f, 1.f };
        uint8_t b[3] = { 1, 1, 1 };
        int q = 0;
        for (q = 0; q < 3; q++) {
            xyz[q * 3] = base + (float)k; xyz[q * 3 + 1] = (float)q; xyz[q * 3 + 2] = (float)(q == 2);
            nrm[q * 3] = 1.f;
        }
        memset(&h, 0, sizeof h);
        h.id = k; h.flags = flags[k]; h.has_uv = no_uv[k] ? 0u : 1u; h.nv = 3; h.nf = 1;
        ok &= fwrite(&h, sizeof h, 1, fp) == 1 && fwrite(vid, sizeof vid, 1, fp) == 1 && fwrite(xyz, sizeof xyz, 1, fp) == 1 &&
              fwrite(nrm, sizeof nrm, 1, fp) == 1;
        if (!no_uv[k]) ok &= fwrite(uv, sizeof uv, 1, fp) == 1;
        ok &= fwrite(faces, sizeof faces, 1, fp) == 1 && fwrite(b, sizeof b, 1, fp) == 1;
    }
    if (fclose(fp) != 0) ok = 0;
    return ok ? 0 : -1;
}

int ChartPool_selftest(void)
{
    char dir[512], store[600], out[600], path[700];
    int fails = 0;
    ChartPoolStats st;
    MeshBinData m;
    MeshChartIndex ix;
    const uint32_t f1[3] = { 0u, ASM_CHART_TINY, 0u }, f2[2] = { ASM_CHART_SUSPECT, ASM_CHART_BLOB };
    const int u1[3] = { 0, 0, 1 }, u2[2] = { 0, 0 };
    long box[6] = { 0, 256, 0, 128, 0, 128 };
    memset(&m, 0, sizeof m);
    if (ves_selftest_dir(dir, sizeof dir, "chart_pool") != 0) return 1;
    snprintf(store, sizeof store, "%s/store", dir);
    snprintf(out, sizeof out, "%s/pool", dir);
    (void)ves_mkdir(store);
    /* cube A: kept, TINY (out), no UV (out); cube B: SUSPECT (kept: in layout), BLOB (out); cube C outside the box */
    snprintf(path, sizeof path, "%s/z00000_y00000_x00000.asc", store);
    if (cp_test_write(path, 3, f1, u1, 0.f) != 0) fails++;
    snprintf(path, sizeof path, "%s/z00128_y00000_x00000.asc", store);
    if (cp_test_write(path, 2, f2, u2, 100.f) != 0) fails++;
    snprintf(path, sizeof path, "%s/z00256_y00000_x00000.asc", store);
    if (cp_test_write(path, 2, f2, u2, 200.f) != 0) fails++;
    snprintf(path, sizeof path, "%s/notes.txt", store);
    { FILE *fp = fopen(path, "wb"); if (fp) { fputs("x", fp); fclose(fp); } }
    if (ChartPool_from_store(store, out, box, &st) != 0) { fprintf(stderr, "  chart_pool selftest FAIL: build\n"); fails++; }
    if (st.cubes != 2 || st.charts != 2 || st.charts_skipped != 3 || st.verts != 6 || st.faces != 2 || st.cubes_failed) {
        fprintf(stderr, "  chart_pool selftest FAIL: counts cubes %zu charts %zu skipped %zu verts %zu faces %zu\n",
                st.cubes, st.charts, st.charts_skipped, st.verts, st.faces);
        fails++;
    }
    snprintf(path, sizeof path, "%s/stage6_repaired.vmesh", out);
    if (MeshBin_read_precise_malloc(path, &m) != 0 || m.uv64 == NULL || m.nv != 6 || m.nf != 2) {
        fprintf(stderr, "  chart_pool selftest FAIL: read back\n");
        fails++;
    } else {
        int32_t ch[6];
        FILE *fp = NULL;
        char cpath[700];
        snprintf(cpath, sizeof cpath, "%s/stage6_repaired_chart.i32", out);
        fp = fopen(cpath, "rb");
        if (fp == NULL || fread(ch, sizeof ch[0], 6, fp) != 6) fails++;
        if (fp) fclose(fp);
        /* the faces rebased per chart; chart ids 0, 0, 0, 1, 1, 1; the first chart of cube B at z 100 */
        if (m.faces[3] != 3 || m.faces[5] != 5 || ch[0] != 0 || ch[2] != 0 || ch[3] != 1 || ch[5] != 1 ||
            fabs(m.verts[9] - 100.0) > 1e-6 || fabs(m.uv64[2] - 1.0) > 1e-12 || fabs(m.uv64[11] - 1.0) > 1e-12) {
            fprintf(stderr, "  chart_pool selftest FAIL: contents\n");
            fails++;
        }
        if (MeshChartIndex_path(path, cpath, sizeof cpath) != 0 || MeshChartIndex_load(cpath, path, &ix) != 0) {
            fprintf(stderr, "  chart_pool selftest FAIL: no chart index\n");
            fails++;
        } else {
            if (ix.n != 2 || ix.rec[1].chart != 1 || ix.rec[1].v0 != 3 || ix.rec[1].f0 != 1) {
                fprintf(stderr, "  chart_pool selftest FAIL: index records\n");
                fails++;
            }
            MeshChartIndex_free(&ix);
        }
    }
    MeshBin_dispose(&m);
    {   /* an empty box: nothing qualifies, nothing written */
        long none[6] = { 9000, 9128, 0, 128, 0, 128 };
        char out2[700];
        snprintf(out2, sizeof out2, "%s/empty", dir);
        if (ChartPool_from_store(store, out2, none, &st) != 1 || st.cubes != 0) {
            fprintf(stderr, "  chart_pool selftest FAIL: empty box\n");
            fails++;
        }
    }
    {   /* a malformed file is skipped whole, the rest still pools */
        char out3[700];
        long all[6] = { 0, 256, 0, 256, 0, 128 };
        snprintf(path, sizeof path, "%s/z00000_y00128_x00000.asc", store);
        { FILE *fp = fopen(path, "wb"); if (fp) { fputs("ASMCHT02 broken", fp); fclose(fp); } }
        snprintf(out3, sizeof out3, "%s/pool3", dir);
        if (ChartPool_from_store(store, out3, all, &st) != 0 || st.cubes_failed != 1 || st.charts != 2) {
            fprintf(stderr, "  chart_pool selftest FAIL: malformed file\n");
            fails++;
        }
    }
    fprintf(stderr, "  chart_pool selftest: %s\n", fails ? "FAIL" : "pass");
    if (fails == 0) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}
