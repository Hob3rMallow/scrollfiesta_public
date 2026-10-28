#include "zarr_u8.h"
#include "ves_platform.h"
#include "pipeline_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <sys/mman.h>
#include <fcntl.h>
#endif

/* (2026-10-23) A chunk slot's buffer as whole pages from the system, so a closed handle's chunks go back to it at once: a
 * 2 MB malloc block went back to its thread's malloc arena, and on a 21^3 run 64 workers, each closing a piece's
 * handles after the piece, held 174 GiB of freed slots in their arenas.  NULL on failure. */
static uint8_t *zu_slot_alloc(size_t n)
{
#ifdef _WIN32
    return (uint8_t *)VirtualAlloc(NULL, n, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
#else
    void *p = mmap(NULL, n, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    return p == MAP_FAILED ? NULL : (uint8_t *)p;
#endif
}

static void zu_slot_free(uint8_t *p, size_t n)
{
    if (p == NULL) return;
#ifdef _WIN32
    (void)n;
    VirtualFree(p, 0, MEM_RELEASE);
#else
    munmap(p, n);
#endif
}

#if ZARR_U8_MAPPED_CHUNKS
/* Fibre samples touch a thin surface inside each immutable volume chunk. Mapping
 * faults only the pages sampled instead of reading all 2 MiB at every miss.
 * Keep the original buffered reader as a fallback when mapping is unavailable. */
static uint8_t *zu_map(const char *path, size_t bytes)
{
#ifdef _WIN32
    HANDLE file = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
                             NULL, OPEN_EXISTING, FILE_FLAG_RANDOM_ACCESS, NULL);
    HANDLE mapping = NULL;
    LARGE_INTEGER size;
    uint8_t *data = NULL;
    if (file == INVALID_HANDLE_VALUE) return NULL;
    if (GetFileSizeEx(file, &size) && size.QuadPart >= 0 && (uint64_t)size.QuadPart >= (uint64_t)bytes)
        mapping = CreateFileMappingA(file, NULL, PAGE_READONLY, 0, 0, NULL);
    if (mapping != NULL) {
        data = (uint8_t *)MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, bytes);
        CloseHandle(mapping);
    }
    CloseHandle(file);
    return data;
#else
    int fd = open(path, O_RDONLY);
    struct stat st;
    void *data = MAP_FAILED;
    if (fd < 0) return NULL;
    if (fstat(fd, &st) == 0 && st.st_size >= 0 && (uint64_t)st.st_size >= (uint64_t)bytes)
        data = mmap(NULL, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (data == MAP_FAILED) return NULL;
#ifdef MADV_RANDOM
    (void)madvise(data, bytes, MADV_RANDOM);
#endif
    return (uint8_t *)data;
#endif
}
#endif

static void zu_buffer_free(uint8_t *data, size_t bytes, int mapped)
{
    if (!data) return;
    if (!mapped) { zu_slot_free(data, bytes); return; }
#ifdef _WIN32
    (void)bytes;
    UnmapViewOfFile(data);
#else
    munmap(data, bytes);
#endif
}

#include "zarr_u8_shared.inc"

static const char *zu_find(const char *hay, const char *key)
{
    const char *p = strstr(hay, key);
    return p != NULL ? p + strlen(key) : NULL;
}

static int zu_read_triple(const char *p, long out[3])
{
    /* "[20974,6621,6621]" possibly with spaces */
    long a = 0, b = 0, c = 0;
    if (p == NULL) return -1;
    while (*p == ' ' || *p == ':') p++;
    if (sscanf(p, "[ %ld , %ld , %ld ]", &a, &b, &c) != 3) return -1;
    out[0] = a; out[1] = b; out[2] = c;
    return 0;
}

int ZarrU8_open(ZarrU8 *z, const char *root)
{
    char path[1400];
    char text[4096];
    FILE *f = NULL;
    size_t n = 0;
    if (z == NULL || root == NULL) return -1;
    memset(z, 0, sizeof *z);
    snprintf(z->root, sizeof z->root, "%s", root);
    snprintf(path, sizeof path, "%s/0/.zarray", root);
    f = fopen(path, "rb");
    if (f == NULL) return -1;
    n = fread(text, 1, sizeof text - 1, f);
    fclose(f);
    text[n] = '\0';
    if (zu_read_triple(zu_find(text, "\"shape\""), z->shape) != 0) return -1;
    if (zu_read_triple(zu_find(text, "\"chunks\""), z->chunk) != 0) return -1;
    if (strstr(text, "|u1") == NULL) return -1;
    if (strstr(text, "\"compressor\":null") == NULL &&
        strstr(text, "\"compressor\": null") == NULL) return -1;
    if (z->chunk[0] <= 0 || z->chunk[1] <= 0 || z->chunk[2] <= 0) return -1;
    z->pow2 = 1;                                  /* the fast path: shifts and masks where every chunk side is 2^k */
    for (int d = 0; d < 3; d++) {
        long cd = z->chunk[d]; int s = 0;
        if (cd & (cd - 1)) { z->pow2 = 0; break; }
        while ((1L << s) < cd) s++;
        z->csh[d] = s; z->cmask[d] = cd - 1;
    }
    size_t nbytes = 1;
    for (int d = 0; d < 3; d++) {
        if ((size_t)z->chunk[d] > SIZE_MAX/nbytes) return -1;
        nbytes *= (size_t)z->chunk[d];
    }
    z->shared = zu_shared_new(nbytes);
    return z->shared ? 0 : -1;
}

void ZarrU8_close(ZarrU8 *z)
{
    if (!z) return;
    for (size_t i = 0; i < ZARR_U8_CACHE_SLOTS; i++) zu_shared_release(z->shared,z->slot_entry[i]);
    zu_shared_close(z->shared);
    memset(z,0,sizeof *z);
}

static const uint8_t *zu_remember(ZarrU8 *z, size_t s)
{
    memcpy(z->last_key,z->slot_key[s],sizeof z->last_key);
    z->last_data = z->slot_data[s]; z->last_valid = 1;
    return z->last_data;
}

static const uint8_t *zu_chunk(ZarrU8 *z, long cz, long cy, long cx)
{
    size_t cap = (z->slot_limit && z->slot_limit < ZARR_U8_CACHE_SLOTS) ? z->slot_limit : ZARR_U8_CACHE_SLOTS;
    size_t h = zu_hash(cz,cy,cx)%cap;
    size_t s = h, probe = 0;
    for (probe = 0; probe < cap; probe++) {
        s = (h+probe)%cap;
        if (!z->slot_used[s]) break;
        if (z->slot_key[s][0] == cz && z->slot_key[s][1] == cy && z->slot_key[s][2] == cx) {
            z->hits++;
            return zu_remember(z,s);
        }
    }
    if (probe == cap) s = z->next_evict%cap;
    ZuSharedChunk *entry = zu_shared_acquire(z,cz,cy,cx);
    if (!entry) return NULL; /* The previous slot and last-data pin stay valid. */
    /* Clear the last pointer BEFORE an eviction can free its buffer. Missing
     * chunks use NULL, so validity cannot be inferred from the data pointer. */
    z->last_valid = 0; z->last_data = NULL;
    zu_shared_release(z->shared,z->slot_entry[s]);
    z->slot_entry[s] = entry; z->slot_data[s] = entry->data;
    z->slot_used[s] = 1;
    z->slot_key[s][0] = cz; z->slot_key[s][1] = cy; z->slot_key[s][2] = cx;
    if (probe == cap) z->next_evict = (s+1)%cap;
    z->misses++;
    if (entry->missing) z->missing++;
    return zu_remember(z,s);
}

uint8_t ZarrU8_get(ZarrU8 *z, long iz, long iy, long ix)
{
    const uint8_t *c = NULL;
    long cz = 0, cy = 0, cx = 0, oz = 0, oy = 0, ox = 0;
    if (z == NULL || iz < 0 || iy < 0 || ix < 0 ||
        iz >= z->shape[0] || iy >= z->shape[1] || ix >= z->shape[2])
        return 0;
    if (z->pow2) {
        cz = iz >> z->csh[0]; cy = iy >> z->csh[1]; cx = ix >> z->csh[2];
        oz = iz & z->cmask[0]; oy = iy & z->cmask[1]; ox = ix & z->cmask[2];
    } else {
        cz = iz / z->chunk[0]; cy = iy / z->chunk[1]; cx = ix / z->chunk[2];
        oz = iz % z->chunk[0]; oy = iy % z->chunk[1]; ox = ix % z->chunk[2];
    }
    if (z->last_valid && z->last_key[0] == cz && z->last_key[1] == cy && z->last_key[2] == cx) {
        c = z->last_data; z->hits++;              /* the chunk of the last read: no hash, no probe */
    } else {
        c = zu_chunk(z, cz, cy, cx);              /* may reload the last chunk's slot: the last chunk is this one now */
    }
    if (c == NULL) return 0;
    return c[((size_t)oz * (size_t)z->chunk[1] + (size_t)oy) * (size_t)z->chunk[2] + (size_t)ox];
}

void ZarrU8_get8(ZarrU8 *z, long iz, long iy, long ix, uint8_t v[8])
{
    if (z != NULL && z->pow2 && iz >= 0 && iy >= 0 && ix >= 0 &&
        iz + 1 < z->shape[0] && iy + 1 < z->shape[1] && ix + 1 < z->shape[2]) {
        long oz = iz & z->cmask[0], oy = iy & z->cmask[1], ox = ix & z->cmask[2];
        if (oz < z->cmask[0] && oy < z->cmask[1] && ox < z->cmask[2]) {   /* the block inside one chunk */
            long cz = iz >> z->csh[0], cy = iy >> z->csh[1], cx = ix >> z->csh[2];
            const uint8_t *ch = NULL;
            if (z->last_valid && z->last_key[0] == cz && z->last_key[1] == cy && z->last_key[2] == cx) {
                ch = z->last_data; z->hits++;
            } else {
                ch = zu_chunk(z, cz, cy, cx);
            }
            if (ch == NULL) { memset(v, 0, 8); return; }
            size_t sy = (size_t)z->chunk[2], sz = (size_t)z->chunk[1] * sy;
            size_t b0 = (size_t)oz * sz + (size_t)oy * sy + (size_t)ox;
            v[0] = ch[b0];           v[1] = ch[b0 + 1];
            v[2] = ch[b0 + sy];      v[3] = ch[b0 + sy + 1];
            v[4] = ch[b0 + sz];      v[5] = ch[b0 + sz + 1];
            v[6] = ch[b0 + sz + sy]; v[7] = ch[b0 + sz + sy + 1];
            return;
        }
    }
    for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int cc = 0; cc < 2; cc++)
        v[4 * a + 2 * b + cc] = ZarrU8_get(z, iz + a, iy + b, ix + cc);
}

uint8_t ZarrU8_sample(ZarrU8 *z, double pz, double py, double px)
{
    return ZarrU8_get(z, (long)floor(pz + 0.5), (long)floor(py + 0.5), (long)floor(px + 0.5));
}

int ZarrU8_present(ZarrU8 *z, long iz, long iy, long ix)
{
    long cz = 0, cy = 0, cx = 0;
    if (z == NULL || iz < 0 || iy < 0 || ix < 0 || iz >= z->shape[0] || iy >= z->shape[1] || ix >= z->shape[2]) return 0;
    if (z->pow2) { cz = iz >> z->csh[0]; cy = iy >> z->csh[1]; cx = ix >> z->csh[2]; }
    else { cz = iz / z->chunk[0]; cy = iy / z->chunk[1]; cx = ix / z->chunk[2]; }
    /* A missing or short shared entry has no data. A present dark chunk
     * has its own non-NULL buffer, even when every byte is zero. */
    return zu_chunk(z, cz, cy, cx) != NULL;
}

#include "zarr_u8_shared_test.inc"

int ZarrU8_selftest(void)
{
    /* Present, short and absent chunks, plus cache eviction in the selftest directory. */
    char root[512];
    char path[1400];
    FILE *f = NULL;
    int fails = zu_shared_selftest();
    ZarrU8 z;
    uint8_t *chunk = (uint8_t *)calloc(4 * 4 * 4, 1);
    if (chunk == NULL) return 1;
    if (ves_selftest_dir(root, sizeof root, "zarr_u8") != 0) { free(chunk); return 1; }
    snprintf(path, sizeof path, "%s/0", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/.zarray", root);
    f = fopen(path, "wb");
    if (f == NULL) { free(chunk); return 1; }
    fprintf(f, "{\"shape\":[%d,4,4],\"chunks\":[4,4,4],\"dtype\":\"|u1\",\"fill_value\":0,"
               "\"order\":\"C\",\"filters\":null,\"dimension_separator\":\"/\","
               "\"compressor\":null,\"zarr_format\":2}", (ZARR_U8_CACHE_SLOTS + 4) * 4);
    fclose(f);
    snprintf(path, sizeof path, "%s/0/1", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/1/0", root); ves_mkdir(path);
    for (int i = 0; i < 64; i++) chunk[i] = (uint8_t)(i * 3);
    snprintf(path, sizeof path, "%s/0/1/0/0", root);
    f = fopen(path, "wb");
    if (f != NULL) { fwrite(chunk, 1, 64, f); fclose(f); }
    snprintf(path, sizeof path, "%s/0/2", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/2/0", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/2/0/0", root);
    f = fopen(path, "wb");
    if (f != NULL) { fputc(99, f); fclose(f); }
    if (ZarrU8_open(&z, root) != 0) { fprintf(stderr, "  zarr_u8: open FAIL\n"); fails++; }
    else {
        if (z.shape[0] != (ZARR_U8_CACHE_SLOTS + 4) * 4 || z.chunk[2] != 4) { fprintf(stderr, "  zarr_u8: header FAIL\n"); fails++; }
        /* voxel (5,1,2) -> chunk (1,0,0), offset (1,1,2) -> index 1*16+1*4+2 = 22 -> 66 */
        if (ZarrU8_get(&z, 5, 1, 2) != 66) { fprintf(stderr, "  zarr_u8: value FAIL (%d)\n", ZarrU8_get(&z, 5, 1, 2)); fails++; }
        if (ZarrU8_get(&z, 1, 1, 2) != 0) { fprintf(stderr, "  zarr_u8: missing chunk should read 0\n"); fails++; }
        if (ZarrU8_get(&z, z.shape[0], 0, 0) != 0 || ZarrU8_get(&z, -1, 0, 0) != 0) {
            fprintf(stderr, "  zarr_u8: outside should read 0\n"); fails++;
        }
        if (ZarrU8_sample(&z, 4.6, 1.4, 2.2) != 66) { fprintf(stderr, "  zarr_u8: nearest sample FAIL\n"); fails++; }
        if (z.misses != 2) { fprintf(stderr, "  zarr_u8: cache misses %zu (want 2)\n", z.misses); fails++; }
        /* (2026-10-24) a missing chunk is told apart from a dark one, on the cache's first and later reads alike, without a reload */
        if (ZarrU8_present(&z, 5, 1, 2) != 1 || ZarrU8_present(&z, 1, 1, 2) != 0 || ZarrU8_present(&z, z.shape[0], 0, 0) != 0 ||
            ZarrU8_present(&z, 1, 1, 2) != 0 || ZarrU8_present(&z, 5, 1, 2) != 1) { fprintf(stderr, "  zarr_u8: present FAIL (missing chunk vs dark)\n"); fails++; }
        if (z.misses != 2) { fprintf(stderr, "  zarr_u8: present reloaded a chunk (misses %zu, want 2)\n", z.misses); fails++; }
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
                for (int c = 0; c < 4; c++)
                    if (ZarrU8_get(&z, a + 4, b, c) != chunk[(a * 4 + b) * 4 + c]) {
                        fprintf(stderr, "  zarr_u8: mapped bytes FAIL\n"); fails++;
                    }
        for (int k = 2; k < ZARR_U8_CACHE_SLOTS + 4; k++)
            if (ZarrU8_get(&z, k * 4, 0, 0) != 0) {
                fprintf(stderr, "  zarr_u8: short/missing chunk FAIL %d\n", k); fails++;
            }
        if (ZarrU8_get(&z, 5, 1, 2) != 66) {
            fprintf(stderr, "  zarr_u8: reload after eviction FAIL\n"); fails++;
        }
        ZarrU8_close(&z);
    }
    if (ZarrU8_open(&z, root) == 0) {     /* a one-chunk cache: every switch of chunk reloads, values unchanged */
        z.slot_limit = 1;
        for (int round = 0; round < 3; round++) {
            if (ZarrU8_get(&z, 5, 1, 2) != 66 || ZarrU8_get(&z, 1, 1, 2) != 0) { fprintf(stderr, "  zarr_u8: limited cache value FAIL\n"); fails++; }
        }
        if (z.misses != 6) { fprintf(stderr, "  zarr_u8: limited cache misses %zu (want 6)\n", z.misses); fails++; }
        for (size_t s = 1; s < ZARR_U8_CACHE_SLOTS; s++) if (z.slot_data[s]) { fprintf(stderr, "  zarr_u8: slot %zu used past the limit\n", s); fails++; break; }
        ZarrU8_close(&z);
    }
    if (ZarrU8_open(&z, root) == 0) {     /* the 2 x 2 x 2 reads: every block, inside a chunk, across the two, past the edges */
        int bad = 0;
        for (long iz = -1; iz < 9; iz++) for (long iy = -1; iy < 5; iy++) for (long ix = -1; ix < 5; ix++) {
            uint8_t v[8];
            ZarrU8_get8(&z, iz, iy, ix, v);
            for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int cc = 0; cc < 2; cc++)
                if (v[4 * a + 2 * b + cc] != ZarrU8_get(&z, iz + a, iy + b, ix + cc)) bad++;
        }
        if (bad) { fprintf(stderr, "  zarr_u8: 2x2x2 block reads FAIL (%d wrong)\n", bad); fails++; }
        ZarrU8_close(&z);
    }
    if (ZarrU8_open(&z, root) == 0) {     /* (fast path) reads alternating and random across the two chunks: the same values */
        unsigned long long r = 12345;
        int bad = 0;
        for (int i = 0; i < 2000; i++) {
            r = r * 6364136223846793005ULL + 1442695040888963407ULL;
            long iz = (long)((r >> 33) % 8), iy = (long)((r >> 40) % 4), ix = (long)((r >> 47) % 4);
            int want = iz >= 4 ? (int)((((iz - 4) * 4 + iy) * 4 + ix) * 3) : 0;
            if (ZarrU8_get(&z, iz, iy, ix) != want) bad++;
        }
        if (bad || !z.pow2) { fprintf(stderr, "  zarr_u8: fast path random reads FAIL (%d wrong, pow2 %d)\n", bad, z.pow2); fails++; }
        ZarrU8_close(&z);
    }
    {   /* a 6 x 3 x 3 array in 3 x 3 x 3 chunks (no power of two: the division path), chunk (1,0,0) written */
        char root3[600];
        snprintf(root3, sizeof root3, "%s/three.zarr", root);
        uint8_t c3[27];
        ves_mkdir(root3);
        snprintf(path, sizeof path, "%s/0", root3); ves_mkdir(path);
        snprintf(path, sizeof path, "%s/0/.zarray", root3);
        f = fopen(path, "wb");
        if (f != NULL) {
            fprintf(f, "{\"shape\":[6,3,3],\"chunks\":[3,3,3],\"dtype\":\"|u1\",\"fill_value\":0,\"order\":\"C\",\"filters\":null,"
                       "\"dimension_separator\":\"/\",\"compressor\":null,\"zarr_format\":2}");
            fclose(f);
        }
        snprintf(path, sizeof path, "%s/0/1", root3); ves_mkdir(path);
        snprintf(path, sizeof path, "%s/0/1/0", root3); ves_mkdir(path);
        for (int i = 0; i < 27; i++) c3[i] = (uint8_t)(i + 1);
        snprintf(path, sizeof path, "%s/0/1/0/0", root3);
        f = fopen(path, "wb");
        if (f != NULL) { fwrite(c3, 1, 27, f); fclose(f); }
        if (ZarrU8_open(&z, root3) == 0) {
            int bad = 0;
            for (long iz = 0; iz < 6; iz++) for (long iy = 0; iy < 3; iy++) for (long ix = 0; ix < 3; ix++) {
                int want = iz >= 3 ? (int)(((iz - 3) * 3 + iy) * 3 + ix + 1) : 0;
                if (ZarrU8_get(&z, iz, iy, ix) != want) bad++;
            }
            if (bad || z.pow2) { fprintf(stderr, "  zarr_u8: non-power-of-two chunks FAIL (%d wrong, pow2 %d)\n", bad, z.pow2); fails++; }
            ZarrU8_close(&z);
        } else { fprintf(stderr, "  zarr_u8: open (3^3 chunks) FAIL\n"); fails++; }
        snprintf(path, sizeof path, "%s/0/1/0/0", root3); remove(path);
        snprintf(path, sizeof path, "%s/0/.zarray", root3); remove(path);
        snprintf(path, sizeof path, "%s/0/1/0", root3); ves_rmdir(path);
        snprintf(path, sizeof path, "%s/0/1", root3); ves_rmdir(path);
        snprintf(path, sizeof path, "%s/0", root3); ves_rmdir(path);
        ves_rmdir(root3);
    }
    free(chunk);
    if (fails == 0) ves_selftest_clean(root);
    fprintf(stderr, "[selftest] zarr_u8 %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}
