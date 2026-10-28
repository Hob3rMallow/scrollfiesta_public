/* Immutable buffers shared only by an explicitly cloned handle family.
 * Local hits need no synchronization. Each local slot pins its entry; a miss
 * acquires the replacement before releasing an evicted entry, so an allocation
 * failure cannot corrupt the old slot. At most one extra chunk per concurrent
 * miss is transiently live. Entries with no pins are immediately returned to
 * the OS; no unbounded secondary cache or malloc arena holds chunk pages.
 * Bucket locks serialize publication of a complete read, not voxel sampling.
 * Bucket -> stats is the only nested lock order. */
#define ZU_SHARED_BUCKETS 256

typedef struct ZuSharedChunk {
    struct ZuSharedChunk *next;
    long key[3];
    uint8_t *data;
    size_t refs;
    int missing, mapped;
} ZuSharedChunk;

typedef struct ZuSharedBucket {
    ves_mutex_t lock;
    ZuSharedChunk *head;
} ZuSharedBucket;

typedef struct ZuSharedCache {
    ZuSharedBucket bucket[ZU_SHARED_BUCKETS];
    ves_mutex_t stats_lock;
    ZarrU8CacheStats stats;
    size_t nbytes;
} ZuSharedCache;

static size_t zu_hash(long cz, long cy, long cx)
{
    /* Widen before multiplication: Windows long is 32 bits, and real scroll
     * chunk coordinates overflow the legacy signed-long products. */
    uint64_t h = (uint64_t)cz*73856093u ^ (uint64_t)cy*19349663u ^ (uint64_t)cx*83492791u;
    return (size_t)h;
}

static void zu_mutex_destroy(ves_mutex_t *lock)
{
#ifdef _WIN32
    (void)lock;
#else
    pthread_mutex_destroy(lock);
#endif
}

static ZuSharedCache *zu_shared_new(size_t nbytes)
{
    ZuSharedCache *cache = calloc(1, sizeof *cache);
    if (!cache) return NULL;
    cache->nbytes = nbytes;
    cache->stats.handles = 1;
    ves_mutex_init(&cache->stats_lock);
    for (size_t i = 0; i < ZU_SHARED_BUCKETS; i++) ves_mutex_init(&cache->bucket[i].lock);
    return cache;
}

static ZuSharedChunk *zu_shared_acquire(ZarrU8 *z, long cz, long cy, long cx)
{
    ZuSharedCache *cache = z->shared;
    ZuSharedBucket *bucket = cache->bucket + zu_hash(cz,cy,cx)%ZU_SHARED_BUCKETS;
    ZuSharedChunk *entry = NULL;
    ves_mutex_lock(&bucket->lock);
    for (entry = bucket->head; entry; entry = entry->next) {
        if (entry->key[0] == cz && entry->key[1] == cy && entry->key[2] == cx) {
            entry->refs++;
            z->shared_hits++;
            ves_mutex_lock(&cache->stats_lock);
            cache->stats.shared_hits++;
            ves_mutex_unlock(&cache->stats_lock);
            ves_mutex_unlock(&bucket->lock);
            return entry;
        }
    }
    entry = calloc(1, sizeof *entry);
    if (entry) {
        char path[1500];
        FILE *file = NULL;
        snprintf(path, sizeof path, "%s/0/%ld/%ld/%ld", z->root,cz,cy,cx);
#if ZARR_U8_MAPPED_CHUNKS
        entry->data = zu_map(path,cache->nbytes);
        entry->mapped = entry->data != NULL;
#endif
        if (!entry->data) {
        file = fopen(path,"rb");
        if (file) {
            entry->data = zu_slot_alloc(cache->nbytes);
            if (!entry->data) {
                fclose(file); free(entry); ves_mutex_unlock(&bucket->lock);
                return NULL;
            }
            if (fread(entry->data,1,cache->nbytes,file) != cache->nbytes) entry->missing = 1;
            fclose(file);
        } else entry->missing = 1;
        }
        if (entry->missing) {
            /* Legacy semantics: a short chunk, including its readable prefix,
             * is all fill. A missing entry needs no full zero-filled buffer. */
            zu_buffer_free(entry->data,cache->nbytes,entry->mapped); entry->data = NULL;
        }
        entry->key[0] = cz; entry->key[1] = cy; entry->key[2] = cx;
        entry->refs = 1;
        entry->next = bucket->head; bucket->head = entry;
        z->loads++;
        ves_mutex_lock(&cache->stats_lock);
        cache->stats.loads++;
        cache->stats.chunks++;
        if (entry->data) cache->stats.bytes += cache->nbytes;
        if (cache->stats.chunks > cache->stats.peak_chunks) cache->stats.peak_chunks = cache->stats.chunks;
        if (cache->stats.bytes > cache->stats.peak_bytes) cache->stats.peak_bytes = cache->stats.bytes;
        ves_mutex_unlock(&cache->stats_lock);
    }
    ves_mutex_unlock(&bucket->lock);
    return entry;
}

static void zu_shared_release(ZuSharedCache *cache, ZuSharedChunk *entry)
{
    ZuSharedBucket *bucket = NULL;
    if (!entry) return;
    bucket = cache->bucket + zu_hash(entry->key[0],entry->key[1],entry->key[2])%ZU_SHARED_BUCKETS;
    ves_mutex_lock(&bucket->lock);
    if (--entry->refs == 0) {
        ZuSharedChunk **link = &bucket->head;
        while (*link != entry) link = &(*link)->next;
        *link = entry->next;
        ves_mutex_lock(&cache->stats_lock);
        cache->stats.chunks--;
        if (entry->data) cache->stats.bytes -= cache->nbytes;
        ves_mutex_unlock(&cache->stats_lock);
        zu_buffer_free(entry->data,cache->nbytes,entry->mapped);
        free(entry);
    }
    ves_mutex_unlock(&bucket->lock);
}

int ZarrU8_clone(ZarrU8 *dst, const ZarrU8 *src)
{
    if (!dst || !src || dst == src || !src->shared) return -1;
    memset(dst,0,sizeof *dst);
    memcpy(dst->root,src->root,sizeof dst->root);
    memcpy(dst->shape,src->shape,sizeof dst->shape);
    memcpy(dst->chunk,src->chunk,sizeof dst->chunk);
    memcpy(dst->csh,src->csh,sizeof dst->csh);
    memcpy(dst->cmask,src->cmask,sizeof dst->cmask);
    dst->pow2 = src->pow2; dst->slot_limit = src->slot_limit;
    dst->shared = src->shared;
    ves_mutex_lock(&dst->shared->stats_lock);
    dst->shared->stats.handles++;
    ves_mutex_unlock(&dst->shared->stats_lock);
    return 0;
}

void ZarrU8_cache_stats(const ZarrU8 *z, ZarrU8CacheStats *out)
{
    if (!out) return;
    memset(out,0,sizeof *out);
    if (!z || !z->shared) return;
    ves_mutex_lock(&z->shared->stats_lock);
    *out = z->shared->stats;
    ves_mutex_unlock(&z->shared->stats_lock);
}

static void zu_shared_close(ZuSharedCache *cache)
{
    size_t handles = 0;
    if (!cache) return;
    ves_mutex_lock(&cache->stats_lock);
    handles = --cache->stats.handles;
    ves_mutex_unlock(&cache->stats_lock);
    if (handles) return;
    /* Each handle drops its slots before dropping its family reference. */
    for (size_t i = 0; i < ZU_SHARED_BUCKETS; i++) zu_mutex_destroy(&cache->bucket[i].lock);
    zu_mutex_destroy(&cache->stats_lock);
    free(cache);
}
