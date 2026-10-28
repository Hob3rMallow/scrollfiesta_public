#include "ves_omp.h"

static void zu_shared_expect(int ok, const char *what, int *fails)
{
    if (!ok) { fprintf(stderr,"  zarr shared: %s FAIL\n",what); (*fails)++; }
}

static uint8_t zu_shared_expected(long iz, long iy, long ix)
{
    if (iz < 0 || iz >= 48 || iy < 0 || iy >= 4 || ix < 0 || ix >= 4 || iz/4 >= 10) return 0;
    return (uint8_t)(((iz/4)*19+((iz%4)*16+iy*4+ix)*7+13)%251);
}

static int zu_shared_selftest(void)
{
    const char *root = "zarr_u8_shared_selftest.zarr";
    char path[1500];
    ZarrU8 parent = {0}, keeper = {0}, separate = {0};
    ZarrU8 *reader = calloc(8,sizeof *reader);
    ZarrU8CacheStats cs = {0};
    uint8_t data[64] = {0};
    int fails = 0, bad[8] = {0}, t = 0;
    FILE *file = NULL;
    if (!reader) return 1;
    ves_mkdir(root);
    snprintf(path,sizeof path,"%s/0",root); ves_mkdir(path);
    snprintf(path,sizeof path,"%s/0/.zarray",root); file = fopen(path,"wb");
    if (!file) { fails++; goto done; }
    fprintf(file,"{\"shape\":[48,4,4],\"chunks\":[4,4,4],\"dtype\":\"|u1\",\"compressor\":null}");
    fclose(file); file = NULL;
    for (int c = 0; c < 12; c++) {
        snprintf(path,sizeof path,"%s/0/%d",root,c); ves_mkdir(path);
        snprintf(path,sizeof path,"%s/0/%d/0",root,c); ves_mkdir(path);
        if (c == 10) continue;
        for (int i = 0; i < 64; i++) data[i] = (uint8_t)((c*19+i*7+13)%251);
        snprintf(path,sizeof path,"%s/0/%d/0/0",root,c); file = fopen(path,"wb");
        if (!file) { fails++; goto done; }
        if (fwrite(data,1,c == 11 ? 19 : 64,file) != (size_t)(c == 11 ? 19 : 64)) fails++;
        fclose(file); file = NULL;
    }
    if (ZarrU8_open(&parent,root)) { fails++; goto done; }
    parent.slot_limit = 3;
    /* Concurrent clone and first-read publication, including several callers
     * reaching exactly the same previously absent shared entry. */
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
    for (t = 0; t < 8; t++) {
        if (ZarrU8_clone(reader+t,&parent)) bad[t]++;
        else if (ZarrU8_get(reader+t,1,2,3) != zu_shared_expected(1,2,3)) bad[t]++;
    }
    for (t = 0; t < 8; t++) fails += bad[t];
    if (fails) goto done;
    ZarrU8_cache_stats(&parent,&cs);
    zu_shared_expect(cs.loads == 1 && cs.shared_hits == 7 && cs.bytes == 64 && cs.handles == 9,
                     "concurrent single publication and exact buffer accounting",&fails);
    for (t = 1; t < 8; t++) zu_shared_expect(reader[t].last_data == reader[0].last_data,"one shared buffer",&fails);
    if (ZarrU8_clone(&keeper,reader)) { fails++; goto done; }
    ZarrU8_close(&parent);
    /* Clones survive the original handle. Each worker churns a bounded cache
     * across full, missing, and short chunks; all values are independently
     * generated from coordinates, including get8 across chunk/volume edges. */
    memset(bad,0,sizeof bad);
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
    for (t = 0; t < 8; t++) {
        uint64_t state = 1234u+(uint64_t)t;
        for (int q = 0; q < 4096; q++) {
            uint8_t block[8];
            state = state*6364136223846793005ULL+1442695040888963407ULL;
            long iz = (long)((state >> 17)%48), iy = (long)((state >> 31)%4), ix = (long)((state >> 42)%4);
            if (ZarrU8_get(reader+t,iz,iy,ix) != zu_shared_expected(iz,iy,ix)) bad[t]++;
            ZarrU8_get8(reader+t,iz,iy,ix,block);
            for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int c = 0; c < 2; c++)
                if (block[4*a+2*b+c] != zu_shared_expected(iz+a,iy+b,ix+c)) bad[t]++;
        }
        ZarrU8_close(reader+t);
    }
    for (t = 0; t < 8; t++) fails += bad[t];
    ZarrU8_cache_stats(&keeper,&cs);
    zu_shared_expect(cs.chunks == 0 && cs.bytes == 0 && cs.handles == 1 && cs.peak_chunks <= 8*4,
                     "bounded pins and release after concurrent churn",&fails);
    keeper.slot_limit = 1;
    zu_shared_expect(ZarrU8_get(&keeper,41,1,1) == 0,"missing chunk fill",&fails);
    ZarrU8_cache_stats(&keeper,&cs);
    zu_shared_expect(cs.chunks == 1 && cs.bytes == 0,"missing chunk has no zero buffer",&fails);
    {
        size_t misses = keeper.misses;
        uint8_t block[8];
        for (int q = 0; q < 10; q++) {
            ZarrU8_get8(&keeper,41,1,1,block);
            for (int k = 0; k < 8; k++) if (block[k]) fails++;
        }
        zu_shared_expect(keeper.misses == misses,"missing last-chunk cache",&fails);
    }
    for (int q = 0; q < 20; q++) {
        zu_shared_expect(ZarrU8_get(&keeper,0,0,0) == 13,"evicted real buffer reload",&fails);
        zu_shared_expect(ZarrU8_get(&keeper,44,0,0) == 0,"short chunk is entirely fill",&fails);
        zu_shared_expect(ZarrU8_get(&keeper,0,0,0) == 13,"last pointer after short-chunk eviction",&fails);
        zu_shared_expect(ZarrU8_get(&keeper,40,0,0) == 0,"last pointer after missing eviction",&fails);
    }
#if !ZARR_U8_MAPPED_CHUNKS
    /* Loaded buffers are snapshots, not writable/mapped aliases of the file.
     * Deliberately changing only this fixture proves family isolation. */
    zu_shared_expect(ZarrU8_get(&keeper,0,0,0) == 13,"snapshot initial value",&fails);
    memset(data,77,sizeof data);
    snprintf(path,sizeof path,"%s/0/0/0/0",root); file = fopen(path,"wb");
    if (!file) { fails++; goto done; }
    if (fwrite(data,1,sizeof data,file) != sizeof data) fails++;
    fclose(file); file = NULL;
    if (ZarrU8_clone(reader,&keeper) || ZarrU8_open(&separate,root)) { fails++; goto done; }
    zu_shared_expect(ZarrU8_get(&keeper,0,0,0) == 13 && ZarrU8_get(reader,0,0,0) == 13,
                     "pinned family snapshot after file rewrite",&fails);
    zu_shared_expect(ZarrU8_get(&separate,0,0,0) == 77 && separate.last_data != keeper.last_data,
                     "independent opens are isolated",&fails);
#else
    /* Mapping explicitly opts into immutable-file aliases, not snapshots. */
    if (ZarrU8_open(&separate,root)) { fails++; goto done; }
    zu_shared_expect(ZarrU8_get(&keeper,0,0,0) == 13 && ZarrU8_get(&separate,0,0,0) == 13 &&
                     keeper.last_data != separate.last_data,"independent mappings",&fails);
#endif
    zu_shared_expect(ZarrU8_clone(&keeper,&keeper) == -1,"self clone rejected",&fails);
done:
    if (file) fclose(file);
    ZarrU8_close(&parent); ZarrU8_close(&keeper); ZarrU8_close(&separate);
    for (t = 0; t < 8; t++) ZarrU8_close(reader+t);
    free(reader);
    if (!fails) {
        for (int c = 0; c < 12; c++) {
            snprintf(path,sizeof path,"%s/0/%d/0/0",root,c); if (c != 10) remove(path);
            snprintf(path,sizeof path,"%s/0/%d/0",root,c); ves_rmdir(path);
            snprintf(path,sizeof path,"%s/0/%d",root,c); ves_rmdir(path);
        }
        snprintf(path,sizeof path,"%s/0/.zarray",root); remove(path);
        snprintf(path,sizeof path,"%s/0",root); ves_rmdir(path);
        ves_rmdir(root);
    }
    fprintf(stderr,"[selftest] zarr shared: concurrent publication, clone lifetime, eviction, short/missing, snapshots; %d failures\n",fails);
    return fails;
}
