/* asm_store.c -- per-cube chart store (fixed-width binary, atomic publish). */
#include "asm_store.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_platform.h"
#include "../common/sha256.h"

#include <sys/stat.h>
#ifdef _WIN32
#include <sys/types.h>
#include <windows.h>
#endif

#define AS_MAGIC "ASMCHT02"

/* identity of the source mesh + cleaning the record was made from */
typedef struct AsIdentity {
    uint64_t fingerprint;
    uint64_t src_size;
    int64_t  src_mtime;
    uint32_t path_len;
    uint32_t reserved;
} AsIdentity;

static void as_stat(const char *path, uint64_t *size, int64_t *mtime)
{
    *size = 0; *mtime = 0;
    if (path == NULL || path[0] == 0) return;
#ifdef _WIN32
    struct __stat64 st;
    if (_stat64(path, &st) == 0) { *size = (uint64_t)st.st_size; *mtime = (int64_t)st.st_mtime; }
#else
    struct stat st;
    if (stat(path, &st) == 0) { *size = (uint64_t)st.st_size; *mtime = (int64_t)st.st_mtime; }
#endif
}

uint64_t AsmStore_fingerprint(const char *text)
{
    uint64_t h = 1469598103934665603ull;
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) { h ^= (uint64_t)*p; h *= 1099511628211ull; }
    return h;
}

typedef struct AsHeader {
    int32_t  id, cube, comp, parent;
    uint32_t flags;
    int32_t  component, placed;
    uint32_t has_uv, has_repaired;
    uint64_t nv, nf, n_flipped;
    double   area3d, area_uv, sigma_lo, sigma_hi, stress_frac;
    double   centroid[3];
    double   pose_x, pose_y, pose_theta;
    float    bbox_lo[3], bbox_hi[3];
} AsHeader;

static void as_path(char *out, size_t cap, const char *dir, const char *cube_id, const char *suffix)
{
    snprintf(out, cap, "%s/%s.asc%s", dir, cube_id, suffix);
}

static int as_put(FILE *fp, const void *p, size_t n)
{
    return n == 0 ? 0 : (fwrite(p, 1, n, fp) == n ? 0 : -1);
}

static int as_get(FILE *fp, void *p, size_t n)
{
    return n == 0 ? 0 : (fread(p, 1, n, fp) == n ? 0 : -1);
}

int AsmStore_exists(const char *dir, const char *cube_id)
{
    char path[2048];
    as_path(path, sizeof path, dir, cube_id, "");
    FILE *fp = fopen(path, "rb");
    if (!fp) return 0;
    fclose(fp);
    return 1;
}

int AsmStore_write_cube(const char *dir, const char *cube_id,
                        const char *src_path, uint64_t fingerprint,
                        const AsmChart *charts, size_t n)
{
    char tmp[2048], path[2048];
    as_path(path, sizeof path, dir, cube_id, "");
    as_path(tmp, sizeof tmp, dir, cube_id, ".tmp");
    ves_ensure_parent_dir(path);
    FILE *fp = fopen(tmp, "wb");
    if (!fp) return -1;
    int rc = 0;
    uint64_t count = (uint64_t)n;
    AsIdentity ident;
    memset(&ident, 0, sizeof ident);
    ident.fingerprint = fingerprint;
    as_stat(src_path, &ident.src_size, &ident.src_mtime);
    ident.path_len = src_path ? (uint32_t)strlen(src_path) : 0u;
    rc |= as_put(fp, AS_MAGIC, 8);
    rc |= as_put(fp, &ident, sizeof ident);
    rc |= as_put(fp, src_path ? src_path : "", ident.path_len);
    rc |= as_put(fp, &count, sizeof count);
    for (size_t i = 0; i < n && rc == 0; i++) {
        const AsmChart *c = &charts[i];
        AsHeader h;
        memset(&h, 0, sizeof h);
        h.id = c->id; h.cube = c->cube; h.comp = c->comp; h.parent = c->parent;
        h.flags = c->flags; h.component = c->component; h.placed = c->placed;
        h.has_uv = c->uv != NULL; h.has_repaired = c->repaired != NULL;
        h.nv = (uint64_t)c->nv; h.nf = (uint64_t)c->nf; h.n_flipped = (uint64_t)c->n_flipped;
        h.area3d = c->area3d; h.area_uv = c->area_uv;
        h.sigma_lo = c->sigma_lo; h.sigma_hi = c->sigma_hi; h.stress_frac = c->stress_frac;
        memcpy(h.centroid, c->centroid, sizeof h.centroid);
        h.pose_x = c->pose_x; h.pose_y = c->pose_y; h.pose_theta = c->pose_theta;
        memcpy(h.bbox_lo, c->bbox_lo, sizeof h.bbox_lo);
        memcpy(h.bbox_hi, c->bbox_hi, sizeof h.bbox_hi);
        rc |= as_put(fp, &h, sizeof h);
        rc |= as_put(fp, c->vid, c->nv * sizeof(int32_t));
        rc |= as_put(fp, c->xyz, c->nv * 3 * sizeof(float));
        rc |= as_put(fp, c->nrm, c->nv * 3 * sizeof(float));
        if (c->uv) rc |= as_put(fp, c->uv, c->nv * 2 * sizeof(float));
        rc |= as_put(fp, c->faces, c->nf * 3 * sizeof(int32_t));
        rc |= as_put(fp, c->boundary, c->nv);
        if (c->repaired) rc |= as_put(fp, c->repaired, c->nv);
    }
    if (fclose(fp) != 0) rc = -1;
    if (rc != 0) { remove(tmp); return -1; }
    remove(path);
    if (rename(tmp, path) != 0) { remove(tmp); return -1; }
    return 0;
}

int AsmStore_read_cube(Arena_T arena, const char *dir, const char *cube_id,
                       const char *src_path, uint64_t fingerprint,
                       AsmChart **out_charts, size_t *out_n)
{
    *out_charts = NULL; *out_n = 0;
    char path[2048];
    as_path(path, sizeof path, dir, cube_id, "");
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    char magic[8];
    uint64_t count = 0;
    AsIdentity ident;
    if (as_get(fp, magic, 8) != 0 || memcmp(magic, AS_MAGIC, 8) != 0 ||
        as_get(fp, &ident, sizeof ident) != 0 || ident.path_len > 4096) { fclose(fp); return -1; }
    {
        char stored[4097];
        if (as_get(fp, stored, ident.path_len) != 0) { fclose(fp); return -1; }
        stored[ident.path_len] = 0;
        uint64_t size = 0; int64_t mtime = 0;
        as_stat(src_path, &size, &mtime);
        if (ident.fingerprint != fingerprint || ident.src_size != size || ident.src_mtime != mtime ||
            strcmp(stored, src_path ? src_path : "") != 0) { fclose(fp); return -2; }
    }
    if (as_get(fp, &count, sizeof count) != 0 || count > (1u << 26)) { fclose(fp); return -1; }
    AsmChart *charts = ARENA_CALLOC(arena, (size_t)count, sizeof(AsmChart));
    for (size_t i = 0; i < (size_t)count; i++) {
        AsHeader h;
        if (as_get(fp, &h, sizeof h) != 0) { fclose(fp); return -1; }
        if (h.nv > (1u << 28) || h.nf > (1u << 28)) { fclose(fp); return -1; }
        AsmChart *c = &charts[i];
        c->id = h.id; c->cube = h.cube; c->comp = h.comp; c->parent = h.parent;
        c->flags = h.flags; c->component = h.component; c->placed = h.placed;
        c->nv = (size_t)h.nv; c->nf = (size_t)h.nf; c->n_flipped = (size_t)h.n_flipped;
        c->area3d = h.area3d; c->area_uv = h.area_uv;
        c->sigma_lo = h.sigma_lo; c->sigma_hi = h.sigma_hi; c->stress_frac = h.stress_frac;
        memcpy(c->centroid, h.centroid, sizeof c->centroid);
        c->pose_x = h.pose_x; c->pose_y = h.pose_y; c->pose_theta = h.pose_theta;
        memcpy(c->bbox_lo, h.bbox_lo, sizeof c->bbox_lo);
        memcpy(c->bbox_hi, h.bbox_hi, sizeof c->bbox_hi);
        c->vid = ARENA_ALLOC(arena, c->nv * sizeof(int32_t));
        c->xyz = ARENA_ALLOC(arena, c->nv * 3 * sizeof(float));
        c->nrm = ARENA_ALLOC(arena, c->nv * 3 * sizeof(float));
        c->uv = h.has_uv ? ARENA_ALLOC(arena, c->nv * 2 * sizeof(float)) : NULL;
        c->faces = ARENA_ALLOC(arena, c->nf * 3 * sizeof(int32_t));
        c->boundary = ARENA_ALLOC(arena, c->nv);
        c->repaired = h.has_repaired ? ARENA_ALLOC(arena, c->nv) : NULL;
        int rc = 0;
        rc |= as_get(fp, c->vid, c->nv * sizeof(int32_t));
        rc |= as_get(fp, c->xyz, c->nv * 3 * sizeof(float));
        rc |= as_get(fp, c->nrm, c->nv * 3 * sizeof(float));
        if (c->uv) rc |= as_get(fp, c->uv, c->nv * 2 * sizeof(float));
        rc |= as_get(fp, c->faces, c->nf * 3 * sizeof(int32_t));
        rc |= as_get(fp, c->boundary, c->nv);
        if (c->repaired) rc |= as_get(fp, c->repaired, c->nv);
        if (rc != 0) { fclose(fp); return -1; }
    }
    fclose(fp);
    *out_charts = charts;
    *out_n = (size_t)count;
    return 0;
}

#include "asm_store_run.inc"

int AsmStore_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    char dir[1024];
    const char *tmpdir = getenv("TEMP");
    snprintf(dir, sizeof dir, "%s/asm_store_selftest_%d", tmpdir ? tmpdir : ".", (int)ves_getpid());
    AsmChart c[2];
    memset(c, 0, sizeof c);
    float xyz[9] = { 0,0,0, 0,1,0, 0,0,1 }, nrm[9] = { 1,0,0, 1,0,0, 1,0,0 }, uv[6] = { 0,0, 1,0, 0,1 };
    int32_t vid[3] = { 5, 6, 7 }, faces[3] = { 0, 1, 2 };
    uint8_t bnd[3] = { 1, 1, 1 };
    c[0].id = 3; c[0].cube = 1; c[0].flags = ASM_CHART_SUSPECT; c[0].nv = 3; c[0].nf = 1;
    c[0].vid = vid; c[0].xyz = xyz; c[0].nrm = nrm; c[0].uv = uv; c[0].faces = faces; c[0].boundary = bnd;
    c[0].area3d = 0.5; c[0].pose_theta = 0.25;
    c[1] = c[0]; c[1].id = 4; c[1].uv = NULL; c[1].flags = ASM_CHART_BLOB;
    uint64_t fp1 = AsmStore_fingerprint("opts-A"), fp2 = AsmStore_fingerprint("opts-B");
    if (fp1 == fp2 || fp1 != AsmStore_fingerprint("opts-A")) { fprintf(stderr, "  asm_store selftest FAIL: fingerprint\n"); fails++; }
    if (AsmStore_write_cube(dir, "z00000_y00000_x00000", "some/mesh.vmesh", fp1, c, 2) != 0) { fprintf(stderr, "  asm_store selftest FAIL: write\n"); fails++; }
    AsmChart *r = NULL; size_t n = 0;
    if (AsmStore_read_cube(arena, dir, "z00000_y00000_x00000", "some/mesh.vmesh", fp1, &r, &n) != 0 || n != 2) { fprintf(stderr, "  asm_store selftest FAIL: read\n"); fails++; }
    if (AsmStore_read_cube(arena, dir, "z00000_y00000_x00000", "some/mesh.vmesh", fp2, &r, &n) != -2) { fprintf(stderr, "  asm_store selftest FAIL: another fingerprint must be stale\n"); fails++; }
    if (AsmStore_read_cube(arena, dir, "z00000_y00000_x00000", "other/mesh.vmesh", fp1, &r, &n) != -2) { fprintf(stderr, "  asm_store selftest FAIL: another source must be stale\n"); fails++; }
    if (AsmStore_read_cube(arena, dir, "z00000_y00000_x00000", "some/mesh.vmesh", fp1, &r, &n) != 0 || n != 2) { fprintf(stderr, "  asm_store selftest FAIL: re-read\n"); fails++; }
    else {
        if (r[0].id != 3 || r[0].nv != 3 || r[0].uv == NULL || r[0].uv[2] != 1.0f || r[0].vid[2] != 7 ||
            r[0].pose_theta != 0.25 || r[0].flags != ASM_CHART_SUSPECT) { fprintf(stderr, "  asm_store selftest FAIL: chart 0\n"); fails++; }
        if (r[1].id != 4 || r[1].uv != NULL || r[1].flags != ASM_CHART_BLOB) { fprintf(stderr, "  asm_store selftest FAIL: chart 1\n"); fails++; }
    }
    if (!AsmStore_exists(dir, "z00000_y00000_x00000")) { fprintf(stderr, "  asm_store selftest FAIL: exists\n"); fails++; }
    {
        char path[1200];
        char checkpoint[2048];snprintf(checkpoint,sizeof checkpoint,"%s/complete.asr",dir);
        AsmChart charts[2]={c[0],c[1]};double precise[6]={0,0,1.000000000001,0,0,1};
        charts[0].id=0;charts[1].id=1;charts[0].placed_uv=precise;
        charts[0].placed=1;charts[0].placement_state=ASM_PLACE_SEAM;charts[0].placement_support=19;
        AsmRelation relation={0};relation.a=0;relation.b=1;relation.flags=ASM_REL_DROPPED|ASM_REL_WEAK;
        relation.continuity=16|32;relation.corr_count=1;AsmCorr corr={0};corr.valid=3;corr.gap_a[0]=1.25f;corr.trim_a.inverse[3]=.345678901234;
        AsmRun source={0};source.arena=arena;source.charts=charts;source.n_charts=2;source.rels=&relation;source.n_rels=1;
        source.corr=&corr;source.n_corr=1;source.continuity_ready=1;source.cube_size=256;
        AsmRun copy={0};copy.arena=arena;
        if(AsmStore_write_run(checkpoint,&source) || AsmStore_read_run(checkpoint,&copy) || copy.n_charts!=2 ||
            !copy.continuity_ready || copy.charts[1].uv || copy.charts[0].placement_support!=19 ||
            copy.charts[0].placement_state!=ASM_PLACE_SEAM || memcmp(copy.charts[0].placed_uv,precise,sizeof precise) ||
            memcmp(copy.charts[0].xyz,xyz,sizeof xyz) || memcmp(copy.charts[0].faces,faces,sizeof faces) ||
            memcmp(copy.rels,&relation,sizeof relation) || memcmp(copy.corr,&corr,sizeof corr))fails++;
        /* Private repair bytes round-trip without becoming chart UVs. They
         * share the complete checkpoint checksum and atomic replacement. */
        uint8_t resume[]={0,1,2,3,0xff,0};source.repair_resume=resume;source.repair_resume_size=sizeof resume;
        memset(&copy,0,sizeof copy);copy.arena=arena;
        if(AsmStore_write_run(checkpoint,&source) || AsmStore_read_run(checkpoint,&copy) ||
           copy.repair_resume_size!=sizeof resume || memcmp(copy.repair_resume,resume,sizeof resume) ||
           memcmp(copy.charts[0].placed_uv,precise,sizeof precise) || copy.charts[1].placed!=source.charts[1].placed)fails++;
        FILE *payload=fopen(checkpoint,"r+b");
        if(!payload || fseek(payload,-65,SEEK_END) || fputc(0x37,payload)==EOF)fails++;
        if(payload)fclose(payload);
        AsmRun corrupt_resume={0};corrupt_resume.arena=arena;
        if(!AsmStore_read_run(checkpoint,&corrupt_resume) || corrupt_resume.n_charts)fails++;
        if(AsmStore_write_run(checkpoint,&source))fails++;
        FILE *damaged=fopen(checkpoint,"r+b");
        if(!damaged || fseek(damaged,-1,SEEK_END) || fputc('!',damaged)==EOF)fails++;
        if(damaged)fclose(damaged);
        AsmRun rejected={0};rejected.arena=arena;
        if(!AsmStore_read_run(checkpoint,&rejected) || rejected.n_charts)fails++;
        /* A complete replacement also replaces a corrupt prior snapshot. */
        if(AsmStore_write_run(checkpoint,&source) || AsmStore_read_run(checkpoint,&rejected) || rejected.n_charts!=2)fails++;
        remove(checkpoint);
        snprintf(path, sizeof path, "%s/z00000_y00000_x00000.asc", dir);
        remove(path);
        ves_rmdir(dir);
    }
    Arena_dispose(&arena);
    if (fails == 0) fprintf(stderr, "  asm_store selftest: all passed\n");
    return fails;
}
