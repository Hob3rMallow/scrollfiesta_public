#include "zarr_u8.h"
#include "ves_platform.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    return 0;
}

void ZarrU8_close(ZarrU8 *z)
{
    if (z == NULL) return;
    for (size_t i = 0; i < ZARR_U8_CACHE_SLOTS; i++) free(z->slot_data[i]);
    memset(z, 0, sizeof *z);
}

static const uint8_t *zu_chunk(ZarrU8 *z, long cz, long cy, long cx)
{
    size_t h = (size_t)(cz * 73856093L ^ cy * 19349663L ^ cx * 83492791L) % ZARR_U8_CACHE_SLOTS;
    size_t nbytes = (size_t)z->chunk[0] * (size_t)z->chunk[1] * (size_t)z->chunk[2];
    for (size_t probe = 0; probe < ZARR_U8_CACHE_SLOTS; probe++) {
        size_t s = (h + probe) % ZARR_U8_CACHE_SLOTS;
        if (!z->slot_used[s]) break;
        if (z->slot_key[s][0] == cz && z->slot_key[s][1] == cy && z->slot_key[s][2] == cx) {
            z->hits++;
            return z->slot_data[s];
        }
    }
    /* load into the first free slot from the hash, or evict round-robin */
    {
        size_t s = h, probe = 0;
        char path[1500];
        FILE *f = NULL;
        for (probe = 0; probe < ZARR_U8_CACHE_SLOTS; probe++) {
            s = (h + probe) % ZARR_U8_CACHE_SLOTS;
            if (!z->slot_used[s]) break;
        }
        if (probe == ZARR_U8_CACHE_SLOTS) {
            s = z->next_evict;
            z->next_evict = (z->next_evict + 1) % ZARR_U8_CACHE_SLOTS;
        }
        if (z->slot_data[s] == NULL) {
            z->slot_data[s] = (uint8_t *)malloc(nbytes);
            if (z->slot_data[s] == NULL) return NULL;
        }
        z->slot_used[s] = 1;
        z->slot_key[s][0] = cz; z->slot_key[s][1] = cy; z->slot_key[s][2] = cx;
        z->misses++;
        snprintf(path, sizeof path, "%s/0/%ld/%ld/%ld", z->root, cz, cy, cx);
        f = fopen(path, "rb");
        if (f == NULL || fread(z->slot_data[s], 1, nbytes, f) != nbytes) {
            memset(z->slot_data[s], 0, nbytes);
            z->missing++;
        }
        if (f != NULL) fclose(f);
        return z->slot_data[s];
    }
}

uint8_t ZarrU8_get(ZarrU8 *z, long iz, long iy, long ix)
{
    const uint8_t *c = NULL;
    long cz = 0, cy = 0, cx = 0, oz = 0, oy = 0, ox = 0;
    if (z == NULL || iz < 0 || iy < 0 || ix < 0 ||
        iz >= z->shape[0] || iy >= z->shape[1] || ix >= z->shape[2])
        return 0;
    cz = iz / z->chunk[0]; cy = iy / z->chunk[1]; cx = ix / z->chunk[2];
    oz = iz % z->chunk[0]; oy = iy % z->chunk[1]; ox = ix % z->chunk[2];
    c = zu_chunk(z, cz, cy, cx);
    if (c == NULL) return 0;
    return c[((size_t)oz * (size_t)z->chunk[1] + (size_t)oy) * (size_t)z->chunk[2] + (size_t)ox];
}

uint8_t ZarrU8_sample(ZarrU8 *z, double pz, double py, double px)
{
    return ZarrU8_get(z, (long)floor(pz + 0.5), (long)floor(py + 0.5), (long)floor(px + 0.5));
}

int ZarrU8_selftest(void)
{
    /* a two-chunk synthetic array written to the current directory */
    const char *root = "zarr_u8_selftest.zarr";
    char path[1400];
    FILE *f = NULL;
    int fails = 0;
    ZarrU8 z;
    uint8_t *chunk = (uint8_t *)calloc(4 * 4 * 4, 1);
    if (chunk == NULL) return 1;
    ves_mkdir(root);
    snprintf(path, sizeof path, "%s/0", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/.zarray", root);
    f = fopen(path, "wb");
    if (f == NULL) { free(chunk); return 1; }
    fprintf(f, "{\"shape\":[8,4,4],\"chunks\":[4,4,4],\"dtype\":\"|u1\",\"fill_value\":0,"
               "\"order\":\"C\",\"filters\":null,\"dimension_separator\":\"/\","
               "\"compressor\":null,\"zarr_format\":2}");
    fclose(f);
    snprintf(path, sizeof path, "%s/0/1", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/1/0", root); ves_mkdir(path);
    for (int i = 0; i < 64; i++) chunk[i] = (uint8_t)(i * 3);
    snprintf(path, sizeof path, "%s/0/1/0/0", root);
    f = fopen(path, "wb");
    if (f != NULL) { fwrite(chunk, 1, 64, f); fclose(f); }
    if (ZarrU8_open(&z, root) != 0) { fprintf(stderr, "  zarr_u8: open FAIL\n"); fails++; }
    else {
        if (z.shape[0] != 8 || z.chunk[2] != 4) { fprintf(stderr, "  zarr_u8: header FAIL\n"); fails++; }
        /* voxel (5,1,2) -> chunk (1,0,0), offset (1,1,2) -> index 1*16+1*4+2 = 22 -> 66 */
        if (ZarrU8_get(&z, 5, 1, 2) != 66) { fprintf(stderr, "  zarr_u8: value FAIL (%d)\n", ZarrU8_get(&z, 5, 1, 2)); fails++; }
        if (ZarrU8_get(&z, 1, 1, 2) != 0) { fprintf(stderr, "  zarr_u8: missing chunk should read 0\n"); fails++; }
        if (ZarrU8_get(&z, 9, 0, 0) != 0) { fprintf(stderr, "  zarr_u8: outside should read 0\n"); fails++; }
        if (ZarrU8_sample(&z, 4.6, 1.4, 2.2) != 66) { fprintf(stderr, "  zarr_u8: nearest sample FAIL\n"); fails++; }
        if (z.misses != 2) { fprintf(stderr, "  zarr_u8: cache misses %zu (want 2)\n", z.misses); fails++; }
        ZarrU8_close(&z);
    }
    free(chunk);
    snprintf(path, sizeof path, "%s/0/1/0/0", root); remove(path);
    snprintf(path, sizeof path, "%s/0/.zarray", root); remove(path);
    snprintf(path, sizeof path, "%s/0/1/0", root); ves_rmdir(path);
    snprintf(path, sizeof path, "%s/0/1", root); ves_rmdir(path);
    snprintf(path, sizeof path, "%s/0", root); ves_rmdir(path);
    ves_rmdir(root);
    fprintf(stderr, "[selftest] zarr_u8 %s (%d failures)\n", fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}
