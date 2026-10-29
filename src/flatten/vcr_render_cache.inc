/* Bounded surface tiles and RAW chunk cache. Included after the resident
 * implementation, which remains an independent small-fixture test control. */
#define VCR_SURFACE_COLS 64
#define VCR_CACHE_MB 2048
#define VCR_CACHE_UPLOAD_CHUNKS 64

typedef struct {
    size_t uploaded, evicted, hits, splits, rectangles;
} VcrCacheStats;

typedef struct {
    CudaDrv cu;
    const Scene *scene;
    CudaDrv_ptr xyz[3], seen, count, keys, table, pool, offsets, output, changes;
    int *key, *table_host, *slot_key, *change;
    uint8_t *presence, *stage;
    uint64_t *used, *pinned, visit;
    unsigned stamp;
    size_t slots, grid, stage_chunks;
    int W, H, stride, flip, N, z0, z1;
    float omin, omax, step;
    VcrCacheStats stats;
    char *why;
    size_t why_len;
} VcrCache;

static int vcr_key_order(const void *a, const void *b)
{
    const int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

static int vcr_cache_error(VcrCache *c)
{
    snprintf(c->why, c->why_len, "%s", CudaDrv_error(c->cu));
    return -1;
}

/* Every needed resident slot is pinned before choosing any victims. Missing
 * RAW chunks are remembered separately; a read error is never fill data. */
static int vcr_cache_load(VcrCache *c, unsigned count)
{
    const Scene *s = c->scene;
    int errors = 0, i = 0, present = 0, changes = 0;
    qsort(c->key, count, sizeof *c->key, vcr_key_order);
    #pragma omp parallel for schedule(dynamic, 16) reduction(+ : errors)
    for (i = 0; i < (int)count; i++) {
        const int key = c->key[i];
        if (!c->presence[key]) {
            const long cz = key / (s->cg[1] * s->cg[2]);
            const long cy = (key / s->cg[2]) % s->cg[1], cx = key % s->cg[2];
            const int result = s->load(s->ctx, cz, cy, cx, NULL);
            if (result < 0) errors++;
            else c->presence[key] = (uint8_t)(result == 1 ? 2 : 1);
        }
    }
    if (errors) {
        snprintf(c->why, c->why_len, "volume chunk listing failed");
        return -1;
    }
    c->visit++;
    for (i = 0; i < (int)count; i++) {
        const int key = c->key[i], slot = c->table_host[key];
        if (c->presence[key] != 2) continue;
        present++;
        if (slot >= 0) {
            c->used[slot] = c->pinned[slot] = c->visit;
            c->stats.hits++;
        }
    }
    for (i = 0; i < (int)count; i++) {
        const int key = c->key[i];
        if (c->presence[key] != 2 || c->table_host[key] >= 0) continue;
        size_t victim = c->slots;
        for (size_t slot = 0; slot < c->slots; slot++) {
            if (c->pinned[slot] == c->visit) continue;
            if (victim == c->slots || c->used[slot] < c->used[victim]) victim = slot;
            if (c->slot_key[slot] < 0) { victim = slot; break; }
        }
        if (victim == c->slots) {
            snprintf(c->why, c->why_len, "RAW cache has no unpinned slot");
            return -1;
        }
        const int old_key = c->slot_key[victim];
        if (old_key >= 0) { c->table_host[old_key] = -1; c->stats.evicted++; }
        c->table_host[key] = (int)victim;
        c->slot_key[victim] = key;
        c->used[victim] = c->pinned[victim] = c->visit;
        c->change[3 * changes] = old_key;
        c->change[3 * changes + 1] = key;
        c->change[3 * changes + 2] = (int)victim;
        changes++;
    }
    for (size_t first = 0; first < (size_t)changes; first += c->stage_chunks) {
        const size_t remain = (size_t)changes - first;
        const int n = (int)(remain < c->stage_chunks ? remain : c->stage_chunks);
        errors = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+ : errors)
        for (i = 0; i < n; i++) {
            const int key = c->change[3 * (first + (size_t)i) + 1];
            const long cz = key / (s->cg[1] * s->cg[2]);
            const long cy = (key / s->cg[2]) % s->cg[1], cx = key % s->cg[2];
            if (s->load(s->ctx, cz, cy, cx, c->stage + (size_t)i * VCR_CHUNK3) != 1) errors++;
        }
        if (errors) {
            snprintf(c->why, c->why_len, "volume chunk read failed");
            return -1;
        }
        for (int start = 0; start < n;) {
            const int slot = c->change[3 * (first + (size_t)start) + 2];
            int end = start + 1;
            while (end < n && c->change[3 * (first + (size_t)end) + 2] == slot + end - start) end++;
            if (CudaDrv_upload(c->cu, c->pool + (CudaDrv_ptr)((size_t)slot * VCR_CHUNK3),
                               c->stage + (size_t)start * VCR_CHUNK3, (size_t)(end - start) * VCR_CHUNK3))
                return vcr_cache_error(c);
            start = end;
        }
    }
    if (changes) {
        void *args[] = {&c->table, &c->changes, &changes};
        if (CudaDrv_upload(c->cu, c->changes, c->change, (size_t)changes * 3 * sizeof(int)) ||
            CudaDrv_launch(c->cu, K_PATCH, ((unsigned)changes + 127) / 128, 1, 128, 1, args))
            return vcr_cache_error(c);
    }
    c->stats.uploaded += (size_t)changes;
    return present;
}

/* Mark only this rectangle, returning a bounded list of unique chunk keys.
 * The generation table avoids clearing or transferring a whole volume grid
 * per tile. Overflow causes subdivision; no data or sample is discarded. */
static int vcr_cache_rect(VcrCache *c, int r0, int c0, int rh, int cw)
{
    int sz = (int)c->scene->shape[0], sy = (int)c->scene->shape[1], sx = (int)c->scene->shape[2];
    int cgy = (int)c->scene->cg[1], cgx = (int)c->scene->cg[2];
    unsigned count = 0, capacity = (unsigned)c->slots;
    if (++c->stamp == 0) {
        if (CudaDrv_memset8(c->cu, c->seen, 0, c->grid * sizeof(unsigned))) return vcr_cache_error(c);
        c->stamp = 1;
    }
    void *mark[] = {&c->xyz[0], &c->xyz[1], &c->xyz[2], &c->W, &c->H, &c->flip, &c->omin, &c->omax,
                   &sz, &sy, &sx, &cgy, &cgx, &r0, &c0, &rh, &cw, &c->stamp, &c->seen, &c->count,
                   &c->keys, &capacity};
    if (CudaDrv_memset8(c->cu, c->count, 0, sizeof(unsigned)) ||
        CudaDrv_launch(c->cu, K_MARK_RECT, ((unsigned)cw + 127) / 128, ((unsigned)rh + 3) / 4, 128, 4, mark) ||
        CudaDrv_download(c->cu, &count, c->count, sizeof count)) return vcr_cache_error(c);
    if (count > capacity) {
        c->stats.splits++;
        if (cw == 1 && rh == 1) {
            snprintf(c->why, c->why_len, "one sample footprint exceeds the %zu-chunk RAW cache", c->slots);
            return -1;
        }
        if (cw >= rh && cw > 1) {
            const int half = cw / 2;
            if (vcr_cache_rect(c, r0, c0, rh, half)) return -1;
            return vcr_cache_rect(c, r0, c0 + half, rh, cw - half);
        }
        const int half = rh / 2;
        if (vcr_cache_rect(c, r0, c0, half, cw)) return -1;
        return vcr_cache_rect(c, r0 + half, c0, rh - half, cw);
    }
    c->stats.rectangles++;
    if (!count) return 0; /* output was zeroed once for the entire tile */
    if (CudaDrv_download(c->cu, c->key, c->keys, (size_t)count * sizeof(int))) return vcr_cache_error(c);
    const int present = vcr_cache_load(c, count);
    if (present <= 0) return present;
    if (c->N) {
        void *args[] = {&c->xyz[0], &c->xyz[1], &c->xyz[2], &c->W, &c->H, &c->flip, &r0, &c0, &rh, &cw,
                       &c->table, &c->pool, &sz, &sy, &sx, &cgy, &cgx, &c->offsets, &c->N, &c->output};
        if (CudaDrv_launch(c->cu, K_LAYERS_RECT, ((unsigned)cw + 127) / 128, ((unsigned)rh + 3) / 4,
                           128, 4, args)) return vcr_cache_error(c);
    } else {
        void *args[] = {&c->xyz[0], &c->xyz[1], &c->xyz[2], &c->W, &c->H, &c->flip, &r0, &c0, &rh, &cw,
                       &c->stride, &c->z0, &c->z1, &c->step, &c->table, &c->pool, &sz, &sy, &sx, &cgy, &cgx,
                       &c->output};
        if (CudaDrv_launch(c->cu, K_CTMAX_RECT, ((unsigned)cw + 127) / 128, ((unsigned)rh + 3) / 4,
                           128, 4, args)) return vcr_cache_error(c);
    }
    return 0;
}

static int gpu_render_cached(CudaDrv cu, const Scene *s, const VcrRender_job *job, BandFn band_fn, void *band_ctx,
                             uint8_t *comp_out, char *why, size_t why_len, int tile_cols, size_t slots,
                             VcrCacheStats *stats)
{
    Arena_T arena = Arena_new();
    VcrCache c;
    int result = -1;
    memset(&c, 0, sizeof c);
    c.cu = cu; c.scene = s; c.why = why; c.why_len = why_len;
    c.N = job->composite ? 0 : job->num_slices;
    c.flip = job->flip_normals ? 1 : 0;
    c.z0 = job->comp_start; c.z1 = job->comp_end; c.step = (float)job->slice_step;
    if (s->W < 1 || s->H < 1 || s->W > INT_MAX - 128 || s->H > INT_MAX - 128 ||
        tile_cols < 1 || tile_cols > VCR_SURFACE_COLS || !slots || c.N < 0 || c.N > 512 ||
        (!job->composite && !c.N)) {
        snprintf(why, why_len, "invalid cache render dimensions");
        goto done;
    }
    c.grid = 1;
    for (int k = 0; k < 3; k++) {
        if (s->cg[k] < 1 || (size_t)s->cg[k] > INT_MAX / c.grid) {
            snprintf(why, why_len, "RAW chunk grid exceeds the supported index range");
            goto done;
        }
        c.grid *= (size_t)s->cg[k];
    }
    c.slots = slots < c.grid ? slots : c.grid;
    c.stage_chunks = c.slots < VCR_CACHE_UPLOAD_CHUNKS ? c.slots : VCR_CACHE_UPLOAD_CHUNKS;
    const size_t max_nodes = (size_t)(tile_cols * VCR_CHUNK + 2) * (VCR_CHUNK + 2);
    const size_t output_bytes = (size_t)tile_cols * VCR_TILE * (size_t)(c.N ? c.N : 1);
    /* A few completed tiles restore the writer's independent-file parallelism
     * without returning to a whole-width band. At most four tiles / 128 MiB;
     * a single deeper tile may itself exceed that byte budget. */
    size_t write_tiles = ((size_t)128 << 20) / output_bytes;
    if (write_tiles < 1) write_tiles = 1;
    if (write_tiles > 4 || !c.N) write_tiles = c.N ? 4 : 1;
    const int write_cols = tile_cols * (int)write_tiles;
    float *xyz = (float *)ARENA_ALLOC(arena, max_nodes * 3 * sizeof(float));
    uint8_t *output = (uint8_t *)ARENA_ALLOC(arena, output_bytes * write_tiles);
    c.key = (int *)ARENA_ALLOC(arena, c.slots * sizeof(int));
    c.table_host = (int *)ARENA_ALLOC(arena, c.grid * sizeof(int));
    c.slot_key = (int *)ARENA_ALLOC(arena, c.slots * sizeof(int));
    c.change = (int *)ARENA_ALLOC(arena, c.slots * 3 * sizeof(int));
    c.presence = (uint8_t *)ARENA_CALLOC(arena, c.grid, 1);
    c.stage = (uint8_t *)ARENA_ALLOC(arena, c.stage_chunks * VCR_CHUNK3);
    c.used = (uint64_t *)ARENA_CALLOC(arena, c.slots, sizeof(uint64_t));
    c.pinned = (uint64_t *)ARENA_CALLOC(arena, c.slots, sizeof(uint64_t));
    for (size_t i = 0; i < c.grid; i++) c.table_host[i] = -1;
    for (size_t i = 0; i < c.slots; i++) c.slot_key[i] = -1;
    float offsets[512];
    if (c.N) {
        layer_offsets(c.N, job->slice_step, offsets);
        c.omin = c.omax = offsets[0];
        for (int i = 1; i < c.N; i++) {
            if (offsets[i] < c.omin) c.omin = offsets[i];
            if (offsets[i] > c.omax) c.omax = offsets[i];
        }
    } else {
        const float a = (float)c.z0 * c.step, b = (float)c.z1 * c.step;
        c.omin = a < b ? a : b;
        c.omax = a < b ? b : a;
    }
    for (int k = 0; k < 3; k++)
        if (CudaDrv_alloc(cu, &c.xyz[k], max_nodes * sizeof(float))) goto gpu_fail;
    if (CudaDrv_alloc(cu, &c.seen, c.grid * sizeof(unsigned)) ||
        CudaDrv_memset8(cu, c.seen, 0, c.grid * sizeof(unsigned)) ||
        CudaDrv_alloc(cu, &c.count, sizeof(unsigned)) || CudaDrv_alloc(cu, &c.keys, c.slots * sizeof(int)) ||
        CudaDrv_alloc(cu, &c.table, c.grid * sizeof(int)) ||
        CudaDrv_upload(cu, c.table, c.table_host, c.grid * sizeof(int)) ||
        CudaDrv_alloc(cu, &c.pool, c.slots * VCR_CHUNK3) ||
        CudaDrv_alloc(cu, &c.changes, c.slots * 3 * sizeof(int)) ||
        CudaDrv_alloc(cu, &c.output, output_bytes)) goto gpu_fail;
    if (c.N && (CudaDrv_alloc(cu, &c.offsets, (size_t)c.N * sizeof(float)) ||
                CudaDrv_upload(cu, c.offsets, offsets, (size_t)c.N * sizeof(float)))) goto gpu_fail;
    printf("[vcr] bounded RAW cache: %zu chunks (%zu MiB), surface tiles at most %d x %d; device buffers %zu bytes\n",
           c.slots, c.slots * VCR_CHUNK3 / ((size_t)1 << 20), tile_cols * VCR_CHUNK, VCR_CHUNK,
           c.slots * VCR_CHUNK3 + max_nodes * 12 + output_bytes + c.grid * 8 + c.slots * 16 + sizeof offsets + 4);
    fflush(stdout);
    const int ntx = (s->W + VCR_CHUNK - 1) / VCR_CHUNK, nty = (s->H + VCR_CHUNK - 1) / VCR_CHUNK;
    const float *source[3] = {s->X, s->Y, s->Z};
    double render_seconds = 0, write_seconds = 0;
    for (int ty = 0; ty < nty; ty++) {
        const int r0 = ty * VCR_CHUNK, rh = s->H - r0 < VCR_CHUNK ? s->H - r0 : VCR_CHUNK;
        int pending_cols = 0, pending_start = 0;
        for (int tx0 = 0; tx0 < ntx; tx0 += tile_cols) {
            const double start = wall_s();
            const int ncols = ntx - tx0 < tile_cols ? ntx - tx0 : tile_cols;
            const int c0 = tx0 * VCR_CHUNK, width = ncols * VCR_CHUNK;
            const int cw = s->W - c0 < width ? s->W - c0 : width;
            /* The normal at (r+1,c+1) needs only a two-node trailing halo.
             * At a true image edge, keep its actual size so border rejection
             * is identical to the resident renderer. */
            c.W = s->W - c0 < cw + 2 ? s->W - c0 : cw + 2;
            c.H = s->H - r0 < rh + 2 ? s->H - r0 : rh + 2;
            c.stride = width;
            if (s->nodes && s->nodes(s->node_ctx, r0, c0, c.H, c.W, xyz, max_nodes)) {
                snprintf(why, why_len, "coordinate TIFF read failed at (%d,%d)", r0, c0);
                goto done;
            }
            for (int k = 0; k < 3; k++) {
                float *plane = xyz + (size_t)k * max_nodes;
                if (!s->nodes)
                    for (int y = 0; y < c.H; y++)
                        memcpy(plane + (size_t)y * c.W, source[k] + (size_t)(r0 + y) * s->W + c0,
                               (size_t)c.W * sizeof(float));
                if (CudaDrv_upload(cu, c.xyz[k], plane, (size_t)c.W * c.H * sizeof(float))) goto gpu_fail;
            }
            const size_t bytes = (size_t)ncols * VCR_TILE * (size_t)(c.N ? c.N : 1);
            if (CudaDrv_memset8(cu, c.output, 0, bytes)) goto gpu_fail;
            if (vcr_cache_rect(&c, 0, 0, rh, cw)) goto done;
            uint8_t *tile_output = output + (size_t)pending_cols * VCR_TILE * (size_t)c.N;
            if (CudaDrv_download(cu, tile_output, c.output, bytes)) goto gpu_fail;
            render_seconds += wall_s() - start;
            const double write_start = wall_s();
            if (c.N) {
                pending_cols += ncols;
                if (pending_cols + tile_cols > write_cols || tx0 + ncols == ntx) {
                    if (band_fn(band_ctx, ty, pending_start, output, pending_cols, c.N)) {
                        snprintf(why, why_len, "writing tiles (%d,%d) failed", ty, pending_start);
                        goto done;
                    }
                    pending_start = tx0 + ncols;
                    pending_cols = 0;
                }
            } else {
                for (int y = 0; y < rh; y++)
                    memcpy(comp_out + (size_t)(r0 + y) * s->W + c0, output + (size_t)y * width, (size_t)cw);
            }
            write_seconds += wall_s() - write_start;
        }
        if (ty + 1 == nty || (ty + 1) % 8 == 0) {
            printf("[vcr] band %d/%d  render %.1fs  write %.1fs\n", ty + 1, nty, render_seconds, write_seconds);
            fflush(stdout);
        }
    }
    printf("[vcr] RAW cache: %zu uploads, %zu hits, %zu evictions, %zu splits, %zu rectangles\n",
           c.stats.uploaded, c.stats.hits, c.stats.evicted, c.stats.splits, c.stats.rectangles);
    result = 0;
    goto done;
gpu_fail:
    vcr_cache_error(&c);
done:
    for (int k = 0; k < 3; k++) CudaDrv_free(cu, c.xyz[k]);
    CudaDrv_free(cu, c.seen); CudaDrv_free(cu, c.count); CudaDrv_free(cu, c.keys);
    CudaDrv_free(cu, c.table); CudaDrv_free(cu, c.pool); CudaDrv_free(cu, c.offsets);
    CudaDrv_free(cu, c.output); CudaDrv_free(cu, c.changes);
    if (stats) *stats = c.stats;
    Arena_dispose(&arena);
    return result;
}

static int gpu_render(CudaDrv cu, const Scene *s, const VcrRender_job *job, BandFn band_fn, void *band_ctx,
                      uint8_t *comp_out, char *why, size_t why_len)
{
    int cols = VCR_SURFACE_COLS;
    if (job->zarr_shard_cols) cols -= cols % job->zarr_shard_cols;
    const size_t mb = job->gpu_cache_mb ? job->gpu_cache_mb : VCR_CACHE_MB;
    if (mb < 16 || mb > 65536) {
        snprintf(why, why_len, "GPU cache must be 16..65536 MiB");
        return -1;
    }
    return gpu_render_cached(cu, s, job, band_fn, band_ctx, comp_out, why, why_len, cols, mb / 2, NULL);
}
