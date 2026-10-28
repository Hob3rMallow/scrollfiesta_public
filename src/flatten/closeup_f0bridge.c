/* closeup_f0bridge.c -- see closeup_f0bridge.h. */
#include "closeup_f0bridge.h"
#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/pipeline_constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_omp.h"

#define CF_CUBE 128
#define CF_CUBE_VOX ((size_t)CF_CUBE * CF_CUBE * CF_CUBE)

/* ------------------------------------------------------------------ the cube cache */

typedef struct CfEntry {
    int64_t  key;        /* -1 empty */
    uint8_t *vol;        /* NULL: no cube in any directory (remembered), or evicted */
    uint64_t stamp;      /* last load or preload (eviction order) */
    int      present;    /* 1: a file exists */
} CfEntry;

struct CloseupF0 {
    char   **dirs;
    size_t   ndirs, max_cubes, loaded;
    CfEntry *tab;
    size_t   cap, n;
    uint64_t clock;
};

static size_t cf_hash(int64_t k, size_t cap)
{
    uint64_t x = (uint64_t)k * 0x9E3779B97F4A7C15ull;
    return (size_t)(x >> 19) & (cap - 1);
}

static int64_t cf_key(long cz, long cy, long cx)
{
    return ((int64_t)cz << 42) | ((int64_t)cy << 21) | (int64_t)cx;
}

static CfEntry *cf_find(const struct CloseupF0 *f, int64_t key)
{
    size_t i = 0;
    if (f->cap == 0) return NULL;
    for (i = cf_hash(key, f->cap); f->tab[i].key != -1; i = (i + 1) & (f->cap - 1))
        if (f->tab[i].key == key) return &f->tab[i];
    return NULL;
}

static int cf_grow(struct CloseupF0 *f)
{
    size_t cap = f->cap ? f->cap * 2 : 1024, i = 0;
    CfEntry *tab = (CfEntry *)malloc(cap * sizeof *tab);
    if (tab == NULL) return -1;
    for (i = 0; i < cap; i++) { tab[i].key = -1; tab[i].vol = NULL; tab[i].stamp = 0; tab[i].present = 0; }
    for (i = 0; i < f->cap; i++) {
        size_t t = 0;
        if (f->tab[i].key == -1) continue;
        for (t = cf_hash(f->tab[i].key, cap); tab[t].key != -1; t = (t + 1) & (cap - 1)) {}
        tab[t] = f->tab[i];
    }
    free(f->tab);
    f->tab = tab;
    f->cap = cap;
    return 0;
}

CloseupF0 CloseupF0_open(const char *const *dirs, size_t ndirs, size_t max_cubes)
{
    struct CloseupF0 *f = (struct CloseupF0 *)calloc(1, sizeof *f);
    size_t k = 0;
    if (f == NULL) return NULL;
    f->dirs = (char **)calloc(ndirs ? ndirs : 1, sizeof *f->dirs);
    if (f->dirs == NULL) { free(f); return NULL; }
    for (k = 0; k < ndirs; k++) {
        size_t len = strlen(dirs[k]);
        f->dirs[k] = (char *)malloc(len + 1);
        if (f->dirs[k] == NULL) { CloseupF0_close(&f); return NULL; }
        memcpy(f->dirs[k], dirs[k], len + 1);
        f->ndirs++;
    }
    f->max_cubes = max_cubes ? max_cubes : 1;
    if (cf_grow(f) != 0) { CloseupF0_close(&f); return NULL; }
    return f;
}

void CloseupF0_close(CloseupF0 *pf)
{
    struct CloseupF0 *f = pf ? *pf : NULL;
    size_t k = 0;
    if (f == NULL) return;
    for (k = 0; k < f->cap; k++) free(f->tab[k].vol);
    for (k = 0; k < f->ndirs; k++) free(f->dirs[k]);
    free(f->tab); free(f->dirs); free(f);
    *pf = NULL;
}

/* The decoded cube (cz, cy, cx) into the cache: 1 loaded (or already), 0 no file, -1 error (no memory, bad file). */
static int cf_load(struct CloseupF0 *f, long cz, long cy, long cx)
{
    const int64_t key = cf_key(cz, cy, cx);
    CfEntry *e = cf_find(f, key);
    char path[1200];
    size_t k = 0;
    if (cz < 0 || cy < 0 || cx < 0 || cz >= (1L << 20) || cy >= (1L << 20) || cx >= (1L << 20)) return 0;
    if (e != NULL && (e->vol != NULL || !e->present)) {
        e->stamp = ++f->clock;
        return e->vol != NULL;
    }
    if (e == NULL) {
        size_t i = 0;
        if (2 * (f->n + 1) > f->cap && cf_grow(f) != 0) return -1;
        for (i = cf_hash(key, f->cap); f->tab[i].key != -1; i = (i + 1) & (f->cap - 1)) {}
        e = &f->tab[i];
        e->key = key; e->vol = NULL; e->present = 0;
        f->n++;
    }
    e->stamp = ++f->clock;
    for (k = 0; k < f->ndirs; k++) {
        FILE *fp = NULL;
        snprintf(path, sizeof path, "%s/z%05ld_y%05ld_x%05ld.tif", f->dirs[k], cz * CF_CUBE, cy * CF_CUBE, cx * CF_CUBE);
        fp = fopen(path, "rb");
        if (fp == NULL) continue;
        fclose(fp);
        {
            Arena_T arena = Arena_new();
            uint8_t *vol = NULL, *keep = NULL;
            int D = 0, H = 0, W = 0;
            if (arena == NULL) return -1;
            if (TiffIO_load(arena, path, &vol, &D, &H, &W) != 0 || D != CF_CUBE || H != CF_CUBE || W != CF_CUBE) {
                Arena_dispose(&arena);
                fprintf(stderr, "  [f0bridge] %s: not a %d^3 uint8 cube; ignored\n", path, CF_CUBE);
                continue;
            }
            if (f->loaded >= f->max_cubes) {          /* evict the least recently used decoded cube */
                size_t i = 0, best = (size_t)-1;
                for (i = 0; i < f->cap; i++)
                    if (f->tab[i].key != -1 && f->tab[i].vol != NULL && (best == (size_t)-1 || f->tab[i].stamp < f->tab[best].stamp))
                        best = i;
                if (best != (size_t)-1) { free(f->tab[best].vol); f->tab[best].vol = NULL; f->loaded--; }
            }
            keep = (uint8_t *)malloc(CF_CUBE_VOX);
            if (keep == NULL) { Arena_dispose(&arena); return -1; }
            memcpy(keep, vol, CF_CUBE_VOX);
            Arena_dispose(&arena);
            e->vol = keep; e->present = 1;
            f->loaded++;
            return 1;
        }
    }
    e->present = 0;
    return 0;
}

long CloseupF0_preload(CloseupF0 f, const double lo[3], const double hi[3])
{
    long c0[3], c1[3], z = 0, y = 0, x = 0, have = 0;
    int q = 0;
    for (q = 0; q < 3; q++) {
        if (!(lo[q] <= hi[q])) return 0;
        c0[q] = (long)floor((lo[q] + 0.5) / CF_CUBE);
        c1[q] = (long)floor((hi[q] + 0.5) / CF_CUBE);
    }
    for (z = c0[0]; z <= c1[0]; z++)
        for (y = c0[1]; y <= c1[1]; y++)
            for (x = c0[2]; x <= c1[2]; x++) {
                int r = cf_load(f, z, y, x);
                if (r < 0) return -1;
                have += r;
            }
    return have;
}

int CloseupF0_at(const CloseupF0 f, double z, double y, double x)
{
    const long vz = (long)floor(z + 0.5), vy = (long)floor(y + 0.5), vx = (long)floor(x + 0.5);
    const CfEntry *e = NULL;
    long cz = 0, cy = 0, cx = 0;
    if (!isfinite(z) || !isfinite(y) || !isfinite(x) || vz < 0 || vy < 0 || vx < 0) return -1;
    cz = vz / CF_CUBE; cy = vy / CF_CUBE; cx = vx / CF_CUBE;
    e = cf_find(f, cf_key(cz, cy, cx));
    if (e == NULL || e->vol == NULL) return -1;
    return e->vol[((size_t)(vz - cz * CF_CUBE) * CF_CUBE + (size_t)(vy - cy * CF_CUBE)) * CF_CUBE +
                  (size_t)(vx - cx * CF_CUBE)] != 0;
}

/* Papyrus at p or one of its 6 face neighbours (a medial point read 1 vox off the band still counts): 1, 0, -1. */
static int cf_near(const CloseupF0 f, const double p[3])
{
    static const int d6[7][3] = { {0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1} };
    int k = 0, any = 0;
    for (k = 0; k < 7; k++) {
        int r = CloseupF0_at(f, p[0] + d6[k][0], p[1] + d6[k][1], p[2] + d6[k][2]);
        if (r > 0) return 1;
        if (r == 0) any = 1;
    }
    return any ? 0 : -1;
}

/* ------------------------------------------------------------------ the bridge */

void CloseupF0BridgeOpts_default(CloseupF0BridgeOpts *o)
{
    o->window_vox = CLOSEUP_F0_WINDOW_VOX;
    o->step_vox = 0.5;
    o->shift_vox = CLOSEUP_F0_SHIFT_VOX;
    o->max_run_vox = CLOSEUP_F0_MAX_RUN_VOX;
    o->depth_px = CLOSEUP_F0_DEPTH_PX;
    o->dup_vox = CLOSEUP_F0_DUP_VOX;
    o->dup_margin_px = CLOSEUP_F0_DUP_MARGIN_PX;
    o->ambig_vox = CLOSEUP_F0_AMBIG_VOX;
    o->layer_vox = CLOSEUP_F0_LAYER_VOX;
    o->rings = CLOSEUP_F0_RINGS;
    o->box[0] = o->box[1] = -INFINITY;
    o->box[2] = o->box[3] = INFINITY;
}

enum { CF_OK = 1, CF_NO_CUBE, CF_NO_RUN, CF_SHIFTED, CF_THICK, CF_BROKEN, CF_DUP, CF_AMBIG, CF_LAYER };

/* ---- the duplicate-layer guard: points of grid nodes in a hash of cells of `cell` vox (sorted keys) ---- */

typedef struct CfHash {
    int64_t *key;          /* sorted cell keys */
    int32_t *node;         /* the grid node of each point */
    size_t   n;
    double   cell;
} CfHash;

static int64_t cf_cell_key(const double p[3], double cell)
{
    int64_t c[3];
    int k = 0;
    for (k = 0; k < 3; k++) c[k] = (int64_t)floor(p[k] / cell) + ((int64_t)1 << 20);
    return (c[0] << 42) | (c[1] << 21) | c[2];
}

typedef struct CfKN { int64_t key; int32_t node; } CfKN;

static int cf_kn_cmp(const void *a, const void *b)
{
    const CfKN *x = (const CfKN *)a, *y = (const CfKN *)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    return x->node < y->node ? -1 : x->node > y->node ? 1 : 0;
}

/* The hash of every node k with use[k] (its point: pts[k*3]); 0, -1 on no memory. */
static int cf_hash_build(CfHash *hs, const float *pts, const uint8_t *use, size_t nn, double cell)
{
    size_t k = 0, n = 0;
    CfKN *kn = NULL;
    free(hs->key); free(hs->node);
    memset(hs, 0, sizeof *hs);
    hs->cell = cell;
    for (k = 0; k < nn; k++) n += use[k] != 0;
    kn = (CfKN *)malloc((n ? n : 1) * sizeof *kn);
    hs->key = (int64_t *)malloc((n ? n : 1) * sizeof *hs->key);
    hs->node = (int32_t *)malloc((n ? n : 1) * sizeof *hs->node);
    if (kn == NULL || hs->key == NULL || hs->node == NULL) { free(kn); return -1; }
    for (k = 0, n = 0; k < nn; k++) {
        double p[3];
        if (!use[k]) continue;
        p[0] = pts[k * 3]; p[1] = pts[k * 3 + 1]; p[2] = pts[k * 3 + 2];
        kn[n].key = cf_cell_key(p, cell);
        kn[n].node = (int32_t)k;
        n++;
    }
    qsort(kn, n, sizeof *kn, cf_kn_cmp);
    for (k = 0; k < n; k++) { hs->key[k] = kn[k].key; hs->node[k] = kn[k].node; }
    hs->n = n;
    free(kn);
    return 0;
}

static void cf_hash_free(CfHash *hs) { free(hs->key); free(hs->node); memset(hs, 0, sizeof *hs); }

/* The points committed in this ring (serial, in commit order): open-addressing cell -> chain of nodes. */
typedef struct CfRing {
    int64_t *cell;          /* cap slots: cell key, or INT64_MIN */
    int32_t *head;          /* cap slots: the first entry of that cell */
    int32_t *next, *node;   /* entries */
    size_t   cap, n, max;
    double   reach;
} CfRing;

static int cf_ring_init(CfRing *r, size_t max, double reach)
{
    size_t k = 0;
    memset(r, 0, sizeof *r);
    r->cap = 64;
    while (r->cap < 2 * max + 2) r->cap *= 2;
    r->max = max; r->reach = reach;
    r->cell = (int64_t *)malloc(r->cap * sizeof *r->cell);
    r->head = (int32_t *)malloc(r->cap * sizeof *r->head);
    r->next = (int32_t *)malloc((max ? max : 1) * sizeof *r->next);
    r->node = (int32_t *)malloc((max ? max : 1) * sizeof *r->node);
    if (!r->cell || !r->head || !r->next || !r->node) return -1;
    for (k = 0; k < r->cap; k++) { r->cell[k] = INT64_MIN; r->head[k] = -1; }
    return 0;
}

static void cf_ring_free(CfRing *r) { free(r->cell); free(r->head); free(r->next); free(r->node); memset(r, 0, sizeof *r); }

static size_t cf_ring_slot(const CfRing *r, int64_t key)
{
    size_t i = (size_t)(((uint64_t)key * 0x9E3779B97F4A7C15ull) >> 20) & (r->cap - 1);
    while (r->cell[i] != INT64_MIN && r->cell[i] != key) i = (i + 1) & (r->cap - 1);
    return i;
}

/* 1 when a point committed earlier in this ring (its node's point in pts) lies within reach of p and far from `node`
 * along the surface; else the point is added (0). */
static int cf_ring_check_add(CfRing *r, const float *pts, const CloseupGrid *g, size_t node, const double p[3],
                             double margin_px)
{
    const long i0 = (long)(node % g->ni), j0 = (long)(node / g->ni);
    int dz = 0, dy = 0, dx = 0;
    size_t s = 0;
    int64_t key = 0;
    for (dz = -1; dz <= 1; dz++)
        for (dy = -1; dy <= 1; dy++)
            for (dx = -1; dx <= 1; dx++) {
                double c[3];
                int32_t e = -1;
                c[0] = p[0] + dz * r->reach; c[1] = p[1] + dy * r->reach; c[2] = p[2] + dx * r->reach;
                s = cf_ring_slot(r, cf_cell_key(c, r->reach));
                for (e = r->cell[s] == INT64_MIN ? -1 : r->head[s]; e >= 0; e = r->next[e]) {
                    const size_t b = (size_t)r->node[e];
                    double d3 = 0.0;
                    int q = 0;
                    for (q = 0; q < 3; q++) d3 += ((double)pts[b * 3 + (size_t)q] - p[q]) * ((double)pts[b * 3 + (size_t)q] - p[q]);
                    d3 = sqrt(d3);
                    if (d3 <= r->reach && g->h * hypot((double)((long)(b % g->ni) - i0), (double)((long)(b / g->ni) - j0)) >
                                              3.0 * d3 + margin_px)
                        return 1;
                }
            }
    if (r->n >= r->max) return 0;
    key = cf_cell_key(p, r->reach);
    s = cf_ring_slot(r, key);
    if (r->cell[s] == INT64_MIN) r->cell[s] = key;
    r->node[r->n] = (int32_t)node;
    r->next[r->n] = r->head[s];
    r->head[s] = (int32_t)r->n;
    r->n++;
    return 0;
}

/* 1 when a hashed point within reach of p belongs to a node far from `node` along the surface (frame distance >
 * 3 x the 3-D distance + margin_px). */
static int cf_dup(const CfHash *hs, const float *pts, const CloseupGrid *g, size_t node, const double p[3], double reach,
                  double margin_px)
{
    const long ni = (long)g->ni, i0 = (long)(node % g->ni), j0 = (long)(node / g->ni);
    int dz = 0, dy = 0, dx = 0;
    if (hs->n == 0) return 0;
    for (dz = -1; dz <= 1; dz++)
        for (dy = -1; dy <= 1; dy++)
            for (dx = -1; dx <= 1; dx++) {
                double c[3];
                int64_t key = 0;
                size_t lo = 0, hi = hs->n;
                c[0] = p[0] + dz * hs->cell; c[1] = p[1] + dy * hs->cell; c[2] = p[2] + dx * hs->cell;
                key = cf_cell_key(c, hs->cell);
                while (lo < hi) { size_t mid = (lo + hi) / 2; if (hs->key[mid] < key) lo = mid + 1; else hi = mid; }
                for (; lo < hs->n && hs->key[lo] == key; lo++) {
                    const size_t b = (size_t)hs->node[lo];
                    double d3 = 0.0, duv = 0.0;
                    int q = 0;
                    for (q = 0; q < 3; q++) d3 += ((double)pts[b * 3 + (size_t)q] - p[q]) * ((double)pts[b * 3 + (size_t)q] - p[q]);
                    d3 = sqrt(d3);
                    if (d3 > reach) continue;
                    duv = g->h * hypot((double)((long)(b % g->ni) - i0), (double)((long)(b / g->ni) - j0));
                    if (duv > 3.0 * d3 + margin_px) return 1;
                }
            }
    (void)ni;
    return 0;
}

/* The layer guard: 1 when a hashed point lies within the hash's cell of p + t nrm for some t = +-t0, +-(t0 + 2), ...
 * while t + cell <= t1 (unit normal nrm: points up to t1 along it), is displaced from p more across than beside it
 * (|d . nrm| > d / 2) and belongs to a node far from `node` along the surface (frame distance > 3 x the 3-D distance +
 * margin_px): the snap would stack a second layer onto the surface. */
static int cf_layer(const CfHash *hs, const float *pts, const CloseupGrid *g, size_t node, const double p[3],
                    const double nrm[3], double t0, double t1, double margin_px)
{
    const long i0 = (long)(node % g->ni), j0 = (long)(node / g->ni);
    const double r2 = hs->cell * hs->cell;
    double t = 0.0;
    if (hs->n == 0 || !(t1 >= t0) || !(t0 > 0.0)) return 0;
    for (t = t0; t + hs->cell <= t1 + 1e-9; t += 2.0) {
        int sgn = 0;
        for (sgn = -1; sgn <= 1; sgn += 2) {
            double c[3];
            int dz = 0, dy = 0, dx = 0, q = 0;
            for (q = 0; q < 3; q++) c[q] = p[q] + (double)sgn * t * nrm[q];
            for (dz = -1; dz <= 1; dz++)
                for (dy = -1; dy <= 1; dy++)
                    for (dx = -1; dx <= 1; dx++) {
                        double cc[3];
                        int64_t key = 0;
                        size_t lo = 0, hi = hs->n;
                        cc[0] = c[0] + dz * hs->cell; cc[1] = c[1] + dy * hs->cell; cc[2] = c[2] + dx * hs->cell;
                        key = cf_cell_key(cc, hs->cell);
                        while (lo < hi) { size_t mid = (lo + hi) / 2; if (hs->key[mid] < key) lo = mid + 1; else hi = mid; }
                        for (; lo < hs->n && hs->key[lo] == key; lo++) {
                            const size_t b = (size_t)hs->node[lo];
                            double d3 = 0.0, dn = 0.0, ds = 0.0, duv = 0.0;
                            for (q = 0; q < 3; q++) {
                                const double v = (double)pts[b * 3 + (size_t)q];
                                d3 += (v - p[q]) * (v - p[q]);
                                dn += (v - p[q]) * nrm[q];
                                ds += (v - c[q]) * (v - c[q]);
                            }
                            if (ds > r2) continue;                              /* outside this sample's reach */
                            d3 = sqrt(d3);
                            if (!(fabs(dn) > 0.5 * d3)) continue;               /* beside the point, not across it */
                            duv = g->h * hypot((double)((long)(b % g->ni) - i0), (double)((long)(b / g->ni) - j0));
                            if (duv > 3.0 * d3 + margin_px) return 1;
                        }
                    }
        }
    }
    return 0;
}

/* The snap of node n from its anchor a (pa: the anchor's point) with the prediction pred and unit normal nrm:
 * CF_OK and *out, or the reason. */
static int cf_snap(const CloseupF0 f0, const CloseupF0BridgeOpts *o, const double pred[3], const double nrm[3],
                   const double pa[3], double out[3])
{
    const int K = (int)floor(o->window_vox / o->step_vox + 0.5), ns = 2 * K + 1;
    int vals[257], k = 0, best_s = -1, best_e = -1, s = -1;
    double best_c = 1e300, second_c = 1e300, c = 0.0, len = 0.0, seg = 0.0;
    int nseg = 0, bad = 0;
    if (ns > 257) return CF_NO_CUBE;
    for (k = 0; k < ns; k++) {
        const double t = o->step_vox * (double)(k - K);
        vals[k] = CloseupF0_at(f0, pred[0] + t * nrm[0], pred[1] + t * nrm[1], pred[2] + t * nrm[2]);
        if (vals[k] < 0) return CF_NO_CUBE;
    }
    for (k = 0; k <= ns; k++) {
        const int on = k < ns && vals[k] > 0;
        if (on && s < 0) s = k;
        if (!on && s >= 0) {
            c = o->step_vox * (0.5 * (double)(s + k - 1) - (double)K);
            if (fabs(c) < fabs(best_c)) { second_c = best_c; best_c = c; best_s = s; best_e = k - 1; }
            else if (fabs(c) < fabs(second_c)) second_c = c;
            s = -1;
        }
    }
    if (best_s < 0) return CF_NO_RUN;
    if (fabs(best_c) > o->shift_vox) return CF_SHIFTED;
    len = o->step_vox * (double)(best_e - best_s + 1);
    if (len > o->max_run_vox || best_s == 0 || best_e == ns - 1) return CF_THICK;   /* touching the window: unbounded */
    if (o->ambig_vox > 0.0 && second_c < 1e299 && fabs(second_c) - fabs(best_c) < o->ambig_vox) return CF_AMBIG;
    for (k = 0; k < 3; k++) out[k] = pred[k] + best_c * nrm[k];
    /* the same sheet: the segment from the anchor's point stays on the papyrus (one miss allowed: aliasing) */
    for (k = 0; k < 3; k++) seg += (out[k] - pa[k]) * (out[k] - pa[k]);
    seg = sqrt(seg);
    nseg = (int)ceil(seg / 0.5);
    for (k = 1; k < nseg; k++) {
        const double t = (double)k / (double)nseg;
        double p[3];
        int q = 0;
        for (q = 0; q < 3; q++) p[q] = pa[q] + t * (out[q] - pa[q]);
        if (cf_near(f0, p) != 1) bad++;
    }
    return bad > 1 ? CF_BROKEN : CF_OK;
}

typedef struct CfCand {
    size_t node;
    size_t anchor;       /* the 4-neighbour it continues (smallest depth) */
    double pred[3], nrm[3], pa[3], out[3];
    int    result;
} CfCand;

int CloseupF0Bridge_pass(const CloseupGrid *g, CloseupF0 f0, const CloseupF0BridgeOpts *o, int32_t f0_chart,
                         CloseupObs *obs, CloseupF0BridgeStats *st)
{
    const size_t ni = g->ni, nj = g->nj, nn = ni * nj;
    const double h = g->h;
    float *depth = NULL, *q = NULL, *corr = NULL;
    uint8_t *snap = NULL;
    size_t *queue = NULL, head = 0, tail = 0, k = 0, ncand = 0, cap_cand = 0;
    CfCand *cand = NULL;
    int32_t *vid = NULL, *faces = NULL, *fchart = NULL;
    float *verts = NULL;
    double *uv = NULL;
    int ring = 0, rc = -1;
    long ci = 0;
    /* the guards' hashes: the duplicate guard's reach as their cell (3 vox when it is off and the layer guard on) */
    const int guard = o->dup_vox > 0.0 || o->layer_vox > 0.0;
    const double hcell = o->dup_vox > 0.0 ? o->dup_vox : CLOSEUP_F0_DUP_VOX;
    CfHash hdata, hsnap;
    CfRing ringset;
    uint8_t *isdata = NULL;
    memset(&hdata, 0, sizeof hdata); memset(&hsnap, 0, sizeof hsnap); memset(&ringset, 0, sizeof ringset);
    memset(st, 0, sizeof *st);
    if (nn == 0) return 0;
    depth = (float *)malloc(nn * sizeof *depth);
    q = (float *)malloc(nn * 3 * sizeof *q);
    corr = (float *)calloc(nn * 3, sizeof *corr);
    snap = (uint8_t *)calloc(nn, 1);
    queue = (size_t *)malloc(nn * sizeof *queue);
    if (depth == NULL || q == NULL || corr == NULL || snap == NULL || queue == NULL) goto out;
    if (guard) {   /* the guards' hash of the data there now */
        isdata = (uint8_t *)malloc(nn);
        if (isdata == NULL) goto out;
        for (k = 0; k < nn; k++) isdata[k] = (uint8_t)(g->state[k] == CLOSEUP_NODE_DATA && isfinite(g->pos[k * 3]));
        if (cf_hash_build(&hdata, g->pos, isdata, nn, hcell) != 0) goto out;
    }
    /* depth beyond the chart data: BFS from chart-DATA nodes through bridge-DATA nodes */
    for (k = 0; k < nn; k++) {
        depth[k] = INFINITY;
        if (g->state[k] == CLOSEUP_NODE_DATA && g->chart[k] != f0_chart) { depth[k] = 0.0f; queue[tail++] = k; }
    }
    while (head < tail) {
        const size_t a = queue[head++], i = a % ni, j = a / ni;
        int d = 0;
        for (d = 0; d < 4; d++) {
            const long ii = (long)i + (d == 0) - (d == 1), jj = (long)j + (d == 2) - (d == 3);
            size_t b = 0;
            if (ii < 0 || jj < 0 || (size_t)ii >= ni || (size_t)jj >= nj) continue;
            b = (size_t)jj * ni + (size_t)ii;
            if (g->state[b] != CLOSEUP_NODE_DATA || g->chart[b] != f0_chart || !(depth[b] > depth[a] + (float)h)) continue;
            depth[b] = depth[a] + (float)h;
            queue[tail++] = b;
        }
    }
    /* the F0 cubes near the boundary: every FILL node within rings x h (+ 2 px) of the data, its profile reach */
    {
        const double reach = (double)o->rings * h + 2.0, pad = o->window_vox + (double)o->rings * h + 2.0;
        for (k = 0; k < nn; k++) {
            double lo[3], hi[3];
            int qq = 0;
            if (g->state[k] != CLOSEUP_NODE_FILL || !(g->dist[k] <= (float)reach)) continue;
            if (!isfinite(g->pos[k * 3])) continue;
            for (qq = 0; qq < 3; qq++) { lo[qq] = g->pos[k * 3 + (size_t)qq] - pad; hi[qq] = g->pos[k * 3 + (size_t)qq] + pad; }
            if (CloseupF0_preload(f0, lo, hi) < 0) goto out;
        }
    }
    for (ring = 0; ring < o->rings; ring++) {
        size_t c = 0, nsnap = 0;
        ncand = 0;
        for (k = 0; k < nn; k++) {
            const size_t i = k % ni, j = k / ni;
            const double u = g->u0 + h * (double)i, w = g->v0 + h * (double)j;
            size_t best = (size_t)-1;
            double ccount = 0.0, cm[3] = { 0.0, 0.0, 0.0 }, nl = 0.0;
            int d = 0, qq = 0;
            if (g->state[k] != CLOSEUP_NODE_FILL || snap[k]) continue;
            if (u < o->box[0] || u > o->box[2] || w < o->box[1] || w > o->box[3]) continue;
            for (d = 0; d < 4; d++) {
                const long ii = (long)i + (d == 0) - (d == 1), jj = (long)j + (d == 2) - (d == 3);
                size_t b = 0;
                if (ii < 0 || jj < 0 || (size_t)ii >= ni || (size_t)jj >= nj) continue;
                b = (size_t)jj * ni + (size_t)ii;
                if (!(g->state[b] == CLOSEUP_NODE_DATA || snap[b])) continue;
                if (snap[b]) { for (qq = 0; qq < 3; qq++) cm[qq] += corr[b * 3 + (size_t)qq]; }
                ccount += 1.0;
                if (best == (size_t)-1 || depth[b] < depth[best]) best = b;
            }
            if (best == (size_t)-1) continue;
            if (!(depth[best] + (float)h <= (float)o->depth_px)) { st->deep++; continue; }
            if (!isfinite(g->pos[k * 3]) || !isfinite(g->nrm[k * 3])) continue;
            if (ncand == cap_cand) {
                size_t nc = cap_cand ? cap_cand * 2 : 4096;
                CfCand *grown = (CfCand *)realloc(cand, nc * sizeof *cand);
                if (grown == NULL) goto out;
                cand = grown;
                cap_cand = nc;
            }
            cand[ncand].node = k;
            cand[ncand].anchor = best;
            for (qq = 0; qq < 3; qq++) {
                cand[ncand].pred[qq] = (double)g->pos[k * 3 + (size_t)qq] + cm[qq] / ccount;
                cand[ncand].nrm[qq] = (double)g->nrm[k * 3 + (size_t)qq];
                nl += cand[ncand].nrm[qq] * cand[ncand].nrm[qq];
                cand[ncand].pa[qq] = snap[best] ? (double)q[best * 3 + (size_t)qq] : (double)g->pos[best * 3 + (size_t)qq];
            }
            nl = sqrt(nl);
            if (!(nl > 1e-6)) continue;
            for (qq = 0; qq < 3; qq++) cand[ncand].nrm[qq] /= nl;
            cand[ncand].result = 0;
            ncand++;
        }
        st->candidates += ncand;
        if (ncand == 0) break;
        if (guard && ring > 0 && cf_hash_build(&hsnap, q, snap, nn, hcell) != 0) goto out;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (ci = 0; ci < (long)ncand; ci++) {
            CfCand *e = &cand[ci];
            e->result = cf_snap(f0, o, e->pred, e->nrm, e->pa, e->out);
            if (e->result == CF_OK && o->dup_vox > 0.0 &&
                (cf_dup(&hdata, g->pos, g, e->node, e->out, o->dup_vox, o->dup_margin_px) ||
                 (ring > 0 && cf_dup(&hsnap, q, g, e->node, e->out, o->dup_vox, o->dup_margin_px))))
                e->result = CF_DUP;
            if (e->result == CF_OK && o->layer_vox > 0.0 &&
                (cf_layer(&hdata, g->pos, g, e->node, e->out, e->nrm, hcell, o->layer_vox, o->dup_margin_px) ||
                 (ring > 0 && cf_layer(&hsnap, q, g, e->node, e->out, e->nrm, hcell, o->layer_vox, o->dup_margin_px))))
                e->result = CF_LAYER;
        }
        if (o->dup_vox > 0.0) {
            cf_ring_free(&ringset);
            if (cf_ring_init(&ringset, ncand, o->dup_vox) != 0) goto out;
        }
        for (c = 0; c < ncand; c++) {
            const CfCand *e = &cand[c];
            int qq = 0;
            if (e->result == CF_OK && o->dup_vox > 0.0) {   /* the fronts of this ring against each other */
                for (qq = 0; qq < 3; qq++) q[e->node * 3 + (size_t)qq] = (float)e->out[qq];
                if (cf_ring_check_add(&ringset, q, g, e->node, e->out, o->dup_margin_px)) { st->dup++; continue; }
            }
            switch (e->result) {
            case CF_OK: break;
            case CF_NO_CUBE: st->no_cube++; continue;
            case CF_NO_RUN: st->no_run++; continue;
            case CF_SHIFTED: st->shifted++; continue;
            case CF_THICK: st->thick++; continue;
            case CF_DUP: st->dup++; continue;
            case CF_AMBIG: st->ambiguous++; continue;
            case CF_LAYER: st->layered++; continue;
            default: st->broken++; continue;
            }
            snap[e->node] = 1;
            for (qq = 0; qq < 3; qq++) {
                q[e->node * 3 + (size_t)qq] = (float)e->out[qq];
                corr[e->node * 3 + (size_t)qq] = (float)(e->out[qq] - (double)g->pos[e->node * 3 + (size_t)qq]);
            }
            depth[e->node] = depth[e->anchor] + (float)h;
            nsnap++;
        }
        st->snapped += nsnap;
        if (nsnap == 0) break;
    }
    if (st->snapped > 0) {   /* observation faces: grid quads whose corners are all data or snapped, one snapped */
        size_t nv = 0, nf = 0, used = 0, gated = 0, i = 0, j = 0;
        CloseupObsGate gate;
        const size_t obs0 = obs->n;
        vid = (int32_t *)malloc(nn * sizeof *vid);
        if (vid == NULL) goto out;
        for (k = 0; k < nn; k++) vid[k] = -1;
        for (j = 0; j + 1 < nj; j++)
            for (i = 0; i + 1 < ni; i++) {
                const size_t c4[4] = { j * ni + i, j * ni + i + 1, (j + 1) * ni + i + 1, (j + 1) * ni + i };
                int a = 0, all = 1, one = 0;
                for (a = 0; a < 4; a++) {
                    all &= g->state[c4[a]] == CLOSEUP_NODE_DATA || snap[c4[a]];
                    one |= snap[c4[a]];
                }
                if (!all || !one) continue;
                for (a = 0; a < 4; a++) if (vid[c4[a]] < 0) vid[c4[a]] = (int32_t)nv++;
                nf += 2;
            }
        verts = (float *)malloc((nv ? nv : 1) * 3 * sizeof *verts);
        uv = (double *)malloc((nv ? nv : 1) * 2 * sizeof *uv);
        faces = (int32_t *)malloc((nf ? nf : 1) * 3 * sizeof *faces);
        fchart = (int32_t *)malloc((nf ? nf : 1) * sizeof *fchart);
        if (verts == NULL || uv == NULL || faces == NULL || fchart == NULL) goto out;
        for (k = 0; k < nn; k++) {
            int qq = 0;
            if (vid[k] < 0) continue;
            for (qq = 0; qq < 3; qq++)
                verts[(size_t)vid[k] * 3 + (size_t)qq] = snap[k] ? q[k * 3 + (size_t)qq] : g->pos[k * 3 + (size_t)qq];
            uv[(size_t)vid[k] * 2] = g->u0 + h * (double)(k % ni);
            uv[(size_t)vid[k] * 2 + 1] = g->v0 + h * (double)(k / ni);
        }
        nf = 0;
        for (j = 0; j + 1 < nj; j++)
            for (i = 0; i + 1 < ni; i++) {
                const size_t c4[4] = { j * ni + i, j * ni + i + 1, (j + 1) * ni + i + 1, (j + 1) * ni + i };
                int a = 0, all = 1, one = 0;
                for (a = 0; a < 4; a++) {
                    all &= g->state[c4[a]] == CLOSEUP_NODE_DATA || snap[c4[a]];
                    one |= snap[c4[a]];
                }
                if (!all || !one) continue;
                faces[nf * 3] = vid[c4[0]]; faces[nf * 3 + 1] = vid[c4[1]]; faces[nf * 3 + 2] = vid[c4[2]];
                faces[nf * 3 + 3] = vid[c4[0]]; faces[nf * 3 + 4] = vid[c4[2]]; faces[nf * 3 + 5] = vid[c4[3]];
                fchart[nf] = fchart[nf + 1] = f0_chart;
                nf += 2;
            }
        CloseupObsGate_default(&gate);
        if (CloseupObs_add_faces(obs, verts, uv, faces, nf, NULL, fchart, &gate, &used, &gated) != 0) goto out;
        st->faces = used;
        st->observations = obs->n - obs0;
    }
    rc = 0;
out:
    cf_hash_free(&hdata); cf_hash_free(&hsnap); cf_ring_free(&ringset); free(isdata);
    free(depth); free(q); free(corr); free(snap); free(queue); free(cand); free(vid); free(faces); free(fchart);
    free(verts); free(uv);
    return rc;
}

/* ------------------------------------------------------------------ selftest */

/* Two wraps of a cylinder sheet about the z axis line (y, x) = (128, 128), radii 80 (the sheet) and 72 (the next wrap
 * in), 3 vox thick, in four 128^3 cubes; over the angles [2.20, 2.45] rad the gap between them is papyrus too (the
 * wraps touch).  A seed observes the outer wrap over the angles [0.90, 1.30] for z 16..112; passes of the bridge, each
 * followed by a refit, must carry the data along the outer wrap both ways (angle down to <= 0.60 and up to >= 1.85),
 * every data node within 1.2 vox of radius 80, and never into the touching sector (no data beyond angle 2.25). */
static int cf_test_cubes(const char *dir)
{
    uint8_t *vol = (uint8_t *)malloc(CF_CUBE_VOX);
    char path[700];
    int cy = 0, cx = 0, fails = 0;
    if (vol == NULL) return 1;
    for (cy = 0; cy < 2; cy++)
        for (cx = 0; cx < 2; cx++) {
            size_t z = 0, y = 0, x = 0;
            for (z = 0; z < CF_CUBE; z++)
                for (y = 0; y < CF_CUBE; y++)
                    for (x = 0; x < CF_CUBE; x++) {
                        const double Y = (double)(cy * CF_CUBE) + (double)y - 128.0, X = (double)(cx * CF_CUBE) + (double)x - 128.0;
                        const double r = sqrt(Y * Y + X * X), th = atan2(Y, X);
                        int on = fabs(r - 80.0) <= 1.5 || fabs(r - 72.0) <= 1.5;
                        if (th >= 2.20 && th <= 2.45 && r >= 70.5 && r <= 81.5) on = 1;
                        vol[(z * CF_CUBE + y) * CF_CUBE + x] = on ? 255 : 0;
                    }
            snprintf(path, sizeof path, "%s/z00000_y%05d_x%05d.tif", dir, cy * CF_CUBE, cx * CF_CUBE);
            fails += TiffIO_save(path, vol, CF_CUBE, CF_CUBE, CF_CUBE) != 0;
        }
    free(vol);
    return fails;
}

/* The duplicate-layer guard: ONE wrap of a cylinder (radius 40, axis (y, x) = (64, 64), 3 vox thick) in one cube, a
 * seed over the angles [0, 0.5] for z 20..100 and a frame wide enough for both fronts to pass each other on the far
 * side (u -200..220 against a circumference of 251).  *dups = the data node pairs of one row within 2.5 vox in 3-D but
 * more than 3 x that + 24 px apart in the frame; *span = the row's data extent in px.  0, or 1 on a failure. */
static int cf_dup_test(const char *dir, int guard, size_t *dups, double *span)
{
    const char *dirs[1];
    const double R = 40.0;
    CloseupF0 f0 = NULL;
    CloseupObs obs;
    CloseupGrid g;
    CloseupFitOpts fit;
    CloseupFitStats fs;
    CloseupObsGate gate;
    CloseupF0BridgeOpts bo;
    CloseupF0BridgeStats bs;
    float *verts = NULL;
    double *uv = NULL;
    int32_t *faces = NULL;
    size_t nu = 21, nvv = 41, i = 0, j = 0, nf = 0, used = 0, gated = 0, a = 0, b = 0;
    double umin = 1e300, umax = -1e300;
    int fails = 0, pass = 0;
    memset(&obs, 0, sizeof obs); memset(&g, 0, sizeof g);
    *dups = 0; *span = 0.0;
    dirs[0] = dir;
    f0 = CloseupF0_open(dirs, 1, 8);
    verts = (float *)malloc(nu * nvv * 3 * sizeof *verts);
    uv = (double *)malloc(nu * nvv * 2 * sizeof *uv);
    faces = (int32_t *)malloc((nu - 1) * (nvv - 1) * 6 * sizeof *faces);
    if (f0 == NULL || verts == NULL || uv == NULL || faces == NULL) { fails++; goto done; }
    for (j = 0; j < nvv; j++)
        for (i = 0; i < nu; i++) {
            const double th = 0.5 * (double)i / (double)(nu - 1), z = 20.0 + 2.0 * (double)j;
            const size_t v = j * nu + i;
            verts[v * 3] = (float)z; verts[v * 3 + 1] = (float)(64.0 + R * sin(th)); verts[v * 3 + 2] = (float)(64.0 + R * cos(th));
            uv[v * 2] = R * th; uv[v * 2 + 1] = z;
        }
    for (j = 0; j + 1 < nvv; j++)
        for (i = 0; i + 1 < nu; i++) {
            const int32_t q = (int32_t)(j * nu + i);
            faces[nf * 3] = q; faces[nf * 3 + 1] = q + 1; faces[nf * 3 + 2] = q + (int32_t)nu + 1;
            faces[nf * 3 + 3] = q; faces[nf * 3 + 4] = q + (int32_t)nu + 1; faces[nf * 3 + 5] = q + (int32_t)nu;
            nf += 2;
        }
    CloseupObsGate_default(&gate);
    if (CloseupObs_add_faces(&obs, verts, uv, faces, nf, NULL, NULL, &gate, &used, &gated) != 0 || obs.n == 0) { fails++; goto done; }
    CloseupFitOpts_complete(&fit);
    if (CloseupGrid_alloc(&g, -200.0, 220.0, 14.0, 106.0, 2.0) != 0 || CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) { fails++; goto done; }
    CloseupF0BridgeOpts_default(&bo);
    bo.depth_px = 1000.0;
    if (!guard) bo.dup_vox = 0.0;
    for (pass = 0; pass < 16; pass++) {
        if (CloseupF0Bridge_pass(&g, f0, &bo, 7, &obs, &bs) != 0) { fails++; break; }
        if (bs.snapped == 0) break;
        if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) { fails++; break; }
    }
    j = (size_t)floor((60.0 - g.v0) / g.h + 0.5);          /* the row at z ~60 */
    for (a = 0; a < g.ni; a++) {
        const size_t ka = j * g.ni + a;
        if (g.state[ka] != CLOSEUP_NODE_DATA) continue;
        umin = fmin(umin, g.u0 + g.h * (double)a); umax = fmax(umax, g.u0 + g.h * (double)a);
        for (b = a + 1; b < g.ni; b++) {
            const size_t kb = j * g.ni + b;
            double d3 = 0.0;
            int q = 0;
            if (g.state[kb] != CLOSEUP_NODE_DATA) continue;
            for (q = 0; q < 3; q++) d3 += ((double)g.pos[ka * 3 + (size_t)q] - g.pos[kb * 3 + (size_t)q]) *
                                         ((double)g.pos[ka * 3 + (size_t)q] - g.pos[kb * 3 + (size_t)q]);
            d3 = sqrt(d3);
            if (d3 < 2.5 && g.h * (double)(b - a) > 3.0 * d3 + 24.0) (*dups)++;
        }
    }
    *span = umax > umin ? umax - umin : 0.0;
done:
    CloseupF0_close(&f0);
    CloseupObs_free(&obs); CloseupGrid_free(&g);
    free(verts); free(uv); free(faces);
    return fails;
}

/* cf_snap's run choice on two touching papyrus planes (voxels x 59-60 and 62-63: run centres 59.5 and 62.5, one voxel
 * apart), normal +x: a prediction at x 59.6 snaps to the first (the second is 2.9 vox away); one at x 61.2, nearly
 * between them (1.3 / 1.7 vox), is ambiguous; with the refusal off it takes the nearer (62.5). */
static int cf_ambig_test(const char *dir)
{
    const char *dirs[1];
    char path[700];
    uint8_t *vol = (uint8_t *)malloc(CF_CUBE_VOX);
    CloseupF0 f0 = NULL;
    CloseupF0BridgeOpts o;
    double nrm[3] = { 0.0, 0.0, 1.0 }, out[3], p1[3] = { 64.0, 64.0, 59.6 }, p2[3] = { 64.0, 64.0, 61.2 };
    const double lo[3] = { 0, 0, 0 }, hi[3] = { 127, 127, 127 };
    size_t z = 0, y = 0, x = 0;
    int fails = 0, r1 = 0, r2 = 0, r3 = 0;
    if (vol == NULL) return 1;
    for (z = 0; z < CF_CUBE; z++)
        for (y = 0; y < CF_CUBE; y++)
            for (x = 0; x < CF_CUBE; x++)
                vol[(z * CF_CUBE + y) * CF_CUBE + x] = (x >= 59 && x <= 60) || (x >= 62 && x <= 63) ? 255 : 0;
    snprintf(path, sizeof path, "%s/z00000_y00000_x00000.tif", dir);
    fails += TiffIO_save(path, vol, CF_CUBE, CF_CUBE, CF_CUBE) != 0;
    free(vol);
    dirs[0] = dir;
    f0 = CloseupF0_open(dirs, 1, 4);
    if (f0 == NULL || CloseupF0_preload(f0, lo, hi) != 1) { CloseupF0_close(&f0); return 1; }
    CloseupF0BridgeOpts_default(&o);
    r1 = cf_snap(f0, &o, p1, nrm, p1, out);
    fails += r1 != CF_OK || fabs(out[2] - 59.5) > 0.3;
    r2 = cf_snap(f0, &o, p2, nrm, p2, out);
    fails += r2 != CF_AMBIG;
    o.ambig_vox = 0.0;
    r3 = cf_snap(f0, &o, p2, nrm, p2, out);
    fails += r3 != CF_OK || fabs(out[2] - 62.5) > 0.3;
    if (fails) fprintf(stderr, "  f0bridge selftest FAIL: run choice %d %d %d (want ok at 59.5, ambiguous, ok at 62.5)\n", r1, r2, r3);
    CloseupF0_close(&f0);
    return fails;
}

/* cf_layer on a hand-made grid: data nodes (i < 60, j < 20) on the plane x = 50 at (z, y) = (j, i), h = 1.  A point
 * 7 vox across that plane (x 57) from node (30, 10) is a second layer for a node far along the surface (i 150) and not
 * for one beside the data (i 36: a tight bend); a point near the data's edge but beside it is not (it is the duplicate
 * guard's and the frame's business), nor one 14 vox across (beyond reach), nor anything with the guard's reach 0. */
static int cf_layer_test(void)
{
    CloseupGrid g;
    CfHash hs;
    uint8_t *use = NULL;
    size_t i = 0, j = 0, nn = 0;
    int fails = 0, r[6];
    const double n[3] = { 0.0, 0.0, 1.0 }, across[3] = { 10.0, 30.0, 57.0 }, beside[3] = { 10.0, 61.0, 51.0 },
                 far3[3] = { 10.0, 30.0, 64.0 };
    memset(&g, 0, sizeof g); memset(&hs, 0, sizeof hs);
    if (CloseupGrid_alloc(&g, 0.0, 199.0, 0.0, 19.0, 1.0) != 0) return 1;
    nn = g.ni * g.nj;
    use = (uint8_t *)calloc(nn, 1);
    if (use == NULL) { CloseupGrid_free(&g); return 1; }
    for (j = 0; j < g.nj; j++)
        for (i = 0; i < 60 && i < g.ni; i++) {
            const size_t k = j * g.ni + i;
            g.pos[k * 3] = (float)j; g.pos[k * 3 + 1] = (float)i; g.pos[k * 3 + 2] = 50.0f;
            use[k] = 1;
        }
    if (cf_hash_build(&hs, g.pos, use, nn, 3.0) != 0) { free(use); CloseupGrid_free(&g); return 1; }
    r[0] = cf_layer(&hs, g.pos, &g, 10 * g.ni + 150, across, n, 3.0, 10.0, 24.0);   /* far along it: a second layer */
    r[1] = cf_layer(&hs, g.pos, &g, 10 * g.ni + 36, across, n, 3.0, 10.0, 24.0);    /* near along it: a tight bend */
    r[2] = cf_layer(&hs, g.pos, &g, 10 * g.ni + 150, beside, n, 3.0, 10.0, 24.0);   /* beside the data, not across */
    r[3] = cf_layer(&hs, g.pos, &g, 10 * g.ni + 150, far3, n, 3.0, 10.0, 24.0);     /* 14 vox across: beyond reach */
    r[4] = cf_layer(&hs, g.pos, &g, 10 * g.ni + 150, across, n, 3.0, 0.0, 24.0);    /* reach 0: off */
    r[5] = cf_layer(&hs, g.pos, &g, 10 * g.ni + 150, across, n, 3.0, 6.0, 24.0);    /* reach 6 does not see 7 vox */
    r[4] |= cf_layer(&hs, g.pos, &g, 10 * g.ni + 150, across, n, 3.0, 5.0, 24.0);   /* reach below t0 + cell: none */
    fails += r[0] != 1 || r[1] != 0 || r[2] != 0 || r[3] != 0 || r[4] != 0 || r[5] != 0;
    if (fails) fprintf(stderr, "  f0bridge selftest FAIL: layer guard %d %d %d %d %d %d (want 1 0 0 0 0 0)\n",
                       r[0], r[1], r[2], r[3], r[4], r[5]);
    cf_hash_free(&hs);
    free(use);
    CloseupGrid_free(&g);
    return fails;
}

int CloseupF0Bridge_selftest(void)
{
    char dir[512];
    const char *dirs[2];
    const double R = 80.0;
    CloseupF0 f0 = NULL;
    CloseupObs obs;
    CloseupGrid g;
    CloseupFitOpts fit;
    CloseupFitStats fs;
    CloseupObsGate gate;
    CloseupF0BridgeOpts bo;
    CloseupF0BridgeStats bs;
    float *verts = NULL;
    double *uv = NULL;
    int32_t *faces = NULL;
    size_t nu = 41, nvv = 49, i = 0, j = 0, nf = 0, used = 0, gated = 0, k = 0, data = 0, off = 0, total_snapped = 0;
    double amin = 1e300, amax = -1e300, worst = 0.0;
    int fails = 0, pass = 0;
    memset(&obs, 0, sizeof obs); memset(&g, 0, sizeof g);
    if (ves_selftest_dir(dir, sizeof dir, "closeup_f0bridge") != 0) return 1;
    if (cf_test_cubes(dir) != 0) { fprintf(stderr, "  f0bridge selftest FAIL: cannot write the cubes\n"); fails++; goto done; }
    dirs[0] = "selftest/no_such_dir"; dirs[1] = dir;
    f0 = CloseupF0_open(dirs, 2, 3);           /* 3 slots for 4 cubes: eviction is exercised by the preload */
    if (f0 == NULL) { fails++; goto done; }
    {   /* the cache itself */
        const double lo[3] = { 0, 0, 0 }, hi[3] = { 127, 255, 255 };
        if (CloseupF0_preload(f0, lo, hi) != 4) fails++;
        CloseupF0_close(&f0);
        f0 = CloseupF0_open(dirs, 2, 64);
        if (f0 == NULL || CloseupF0_preload(f0, lo, hi) != 4) { fails++; goto done; }
        if (CloseupF0_at(f0, 50.0, 128.0, 208.0) != 1 || CloseupF0_at(f0, 50.0, 128.0, 204.0) != 0 ||
            CloseupF0_at(f0, 50.0, 128.0, 300.0) != -1 || CloseupF0_at(f0, -1.0, 128.0, 208.0) != -1) {
            fprintf(stderr, "  f0bridge selftest FAIL: cube reads\n");
            fails++;
        }
    }
    /* the seed: the outer wrap over angles [0.90, 1.30] (u = R theta), z 16..112 (v = z), vertices every 1 px */
    verts = (float *)malloc(nu * nvv * 3 * sizeof *verts);
    uv = (double *)malloc(nu * nvv * 2 * sizeof *uv);
    faces = (int32_t *)malloc((nu - 1) * (nvv - 1) * 6 * sizeof *faces);
    if (verts == NULL || uv == NULL || faces == NULL) { fails++; goto done; }
    for (j = 0; j < nvv; j++)
        for (i = 0; i < nu; i++) {
            const double th = 0.90 + 0.40 * (double)i / (double)(nu - 1), z = 16.0 + 2.0 * (double)j;
            const size_t v = j * nu + i;
            verts[v * 3] = (float)z; verts[v * 3 + 1] = (float)(128.0 + R * sin(th)); verts[v * 3 + 2] = (float)(128.0 + R * cos(th));
            uv[v * 2] = R * th; uv[v * 2 + 1] = z;
        }
    for (j = 0; j + 1 < nvv; j++)
        for (i = 0; i + 1 < nu; i++) {
            const int32_t a = (int32_t)(j * nu + i);
            faces[nf * 3] = a; faces[nf * 3 + 1] = a + 1; faces[nf * 3 + 2] = a + (int32_t)nu + 1;
            faces[nf * 3 + 3] = a; faces[nf * 3 + 4] = a + (int32_t)nu + 1; faces[nf * 3 + 5] = a + (int32_t)nu;
            nf += 2;
        }
    CloseupObsGate_default(&gate);
    if (CloseupObs_add_faces(&obs, verts, uv, faces, nf, NULL, NULL, &gate, &used, &gated) != 0 || obs.n == 0) { fails++; goto done; }
    CloseupFitOpts_complete(&fit);
    if (CloseupGrid_alloc(&g, R * 0.30, R * 2.60, 10.0, 118.0, 2.0) != 0) { fails++; goto done; }
    if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) { fails++; goto done; }
    CloseupF0BridgeOpts_default(&bo);
    bo.depth_px = 1000.0;
    for (pass = 0; pass < 12; pass++) {
        if (CloseupF0Bridge_pass(&g, f0, &bo, 7, &obs, &bs) != 0) { fails++; break; }
        total_snapped += bs.snapped;
        if (bs.snapped == 0) break;
        if (CloseupSurface_fit(&g, &obs, &fit, &fs) != 0) { fails++; break; }
    }
    for (k = 0; k < g.ni * g.nj; k++) {
        double Y = 0.0, X = 0.0, r = 0.0, th = 0.0, z = 0.0;
        if (g.state[k] != CLOSEUP_NODE_DATA) continue;
        z = g.pos[k * 3]; Y = g.pos[k * 3 + 1] - 128.0; X = g.pos[k * 3 + 2] - 128.0;
        if (z < 24.0 || z > 104.0) continue;            /* the rows the seed spans, away from its top and bottom */
        r = sqrt(Y * Y + X * X); th = atan2(Y, X);
        data++;
        worst = fmax(worst, fabs(r - R));
        off += fabs(r - R) > 1.2;
        amin = fmin(amin, th); amax = fmax(amax, th);
    }
    if (total_snapped == 0 || amin > 0.60 || amax < 1.85 || amax > 2.25 || off > data / 200) {
        fprintf(stderr, "  f0bridge selftest FAIL: %zu snapped; data angles [%.3f, %.3f] (want <= 0.60, >= 1.85, <= 2.25); "
                "%zu of %zu data nodes off the wrap (worst %.2f vox)\n", total_snapped, amin, amax, off, data, worst);
        fails++;
    }
    {   /* nothing to snap without cubes: every candidate is no_cube, the pass succeeds with 0 */
        const char *none[1] = { "selftest/no_such_dir" };
        CloseupF0 empty = CloseupF0_open(none, 1, 4);
        if (empty == NULL || CloseupF0Bridge_pass(&g, empty, &bo, 7, &obs, &bs) != 0 || bs.snapped != 0 || bs.no_cube == 0) fails++;
        CloseupF0_close(&empty);
    }
    fails += cf_layer_test();   /* the layer guard on a hand-made grid */
    {   /* the run choice and its ambiguity refusal */
        char adir[600];
        snprintf(adir, sizeof adir, "%s/two_planes", dir);
        (void)ves_mkdir(adir);
        fails += cf_ambig_test(adir);
    }
    {   /* the duplicate-layer guard: one wrap, two fronts meeting on the far side */
        char ddir[600], path[700];
        uint8_t *vol = (uint8_t *)malloc(CF_CUBE_VOX);
        size_t z = 0, y = 0, x = 0, dup_on = 0, dup_off = 0;
        double span_on = 0.0, span_off = 0.0;
        snprintf(ddir, sizeof ddir, "%s/one_wrap", dir);
        (void)ves_mkdir(ddir);
        if (vol == NULL) { fails++; goto done; }
        for (z = 0; z < CF_CUBE; z++)
            for (y = 0; y < CF_CUBE; y++)
                for (x = 0; x < CF_CUBE; x++)
                    vol[(z * CF_CUBE + y) * CF_CUBE + x] =
                        fabs(hypot((double)y - 64.0, (double)x - 64.0) - 40.0) <= 1.5 ? 255 : 0;
        snprintf(path, sizeof path, "%s/z00000_y00000_x00000.tif", ddir);
        fails += TiffIO_save(path, vol, CF_CUBE, CF_CUBE, CF_CUBE) != 0;
        free(vol);
        fails += cf_dup_test(ddir, 1, &dup_on, &span_on) + cf_dup_test(ddir, 0, &dup_off, &span_off);
        /* the guard: the fronts meet and stop (no duplicate pair, most of the 251 px circumference covered); off, they
         * overlap where they meet (duplicate pairs) -- the scene does test the guard */
        if (dup_on != 0 || span_on < 180.0 || span_on > 262.0 || dup_off == 0) {
            fprintf(stderr, "  f0bridge selftest FAIL: duplicate guard: on %zu pairs over %.0f px, off %zu pairs over %.0f px\n",
                    dup_on, span_on, dup_off, span_off);
            fails++;
        }
    }
done:
    CloseupF0_close(&f0);
    CloseupObs_free(&obs); CloseupGrid_free(&g);
    free(verts); free(uv); free(faces);
    fprintf(stderr, "  closeup_f0bridge selftest: %s (%zu snapped, data angles [%.3f, %.3f], worst %.2f vox)\n",
            fails ? "FAIL" : "pass", total_snapped, amin, amax, worst);
    if (fails == 0) (void)ves_selftest_clean(dir);
    return fails ? 1 : 0;
}
