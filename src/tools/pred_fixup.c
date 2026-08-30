/*
 * pred_fixup.c -- 2D prediction gap fixup over a whole cube grid.
 *
 *   pred_fixup <grid_dir> <out_dir> [options]
 *   pred_fixup --selftest
 *
 * Reads <grid_dir>/cubes_PRED/z#####_y#####_x#####.tif (128^3 uint8 binary
 * prediction cubes), assembles the WHOLE grid into one volume, closes visible
 * per-z-plane gaps by skeleton-endpoint matching (greedy stable matching on a
 * symmetric score + crossing bans; src/fixup/), enforces cross-plane
 * ConnectionTrack consistency (round-2 re-match + support filter), paints
 * cubic-Bezier joins (ADDITIVE only), and writes a parallel grid dir
 * <out_dir>/cubes_PRED/*.tif that grid_pipeline consumes unchanged
 * (manifest.json and present.json are copied when present).
 *
 * Full-grid z-planes, not per-cube: per-cube fixup would paint different
 * pixels in adjacent cubes' shared halo band and break the bit-identical
 * seam convergence the weld depends on. The whole 4x5x5 volume is 210 MB;
 * per-plane CC labels are the big allocation (int32, 840 MB). 21x21x21
 * (19.4 GB) needs the z-slab streaming variant — this tool refuses grids
 * over ~2.5 GB rather than thrash.
 *
 * Options:
 *   --umb-y F --umb-x F   world L0 umbilicus; arms the radial gate
 *   --reach-max F         hard join reach cap        (default 12)
 *   --reach-safe F        evidence-free reach tier   (default 6)
 *   --radial-dr F         max |r_a - r_b|            (default 4)
 *   --min-support N       ConnectionTrack planes     (default 3)
 *   --min-score F         score floor                (default 0.30)
 *   --paint-radius N      join stroke radius         (default 1)
 *   --no-tracks           round-1 only, no support filter (debug)
 *   --dry-run             no TIFF output (PNGs + manifests only)
 *   --png-all             overlay PNG for every plane
 *   --png-scale N         overlay upscale            (default 1)
 *   --png-max N           cap overlay PNG count      (default 256)
 *   --planes z0:z1        process local plane subrange [z0,z1)
 *   --threads N           OpenMP threads
 *
 * Exit: 0 ok, 1 IO error, 2 usage error, 3 selftest failure.
 */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#endif

#include "../common/arena.h"
#include "../common/pipeline_constants.h"
#include "../common/tiff_io.h"
#include "../common/ves_png.h"
#include "../fixup/bridge_scan.h"
#include "../fixup/fixup_viz.h"
#include "../fixup/join_paint.h"
#include "../fixup/join_tracks.h"
#include "../fixup/slab_split.h"
#include "../fixup/slice_match.h"
#include "../fixup/slice_trace.h"

/* ---------------------------------------------------------------- */
/* Small utilities                                                  */
/* ---------------------------------------------------------------- */

static double pf_now(void)
{
#ifdef _WIN32
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

typedef struct {
    char id[128];
    char path[1024];
    int  vz, vy, vx;   /* source-space voxel origin from the id */
} PfEntry;

static void pf_basename_noext(const char *path, char *out, size_t out_sz)
{
    const char *base = path;
    for (const char *p = path; *p; p++)
        if (*p == '/' || *p == '\\') base = p + 1;
    size_t len = strlen(base);
    if (len >= 4 && strcmp(base + len - 4, ".tif") == 0) len -= 4;
    if (len >= out_sz) len = out_sz - 1;
    memcpy(out, base, len);
    out[len] = '\0';
}

static int pf_entry_cmp(const void *a, const void *b)
{
    return strcmp(((const PfEntry *)a)->id, ((const PfEntry *)b)->id);
}

/* Collect every <dir>/*.tif with a parsable cube id (sorted, malloc'd). */
static int pf_collect_cubes(const char *dir, PfEntry **out, size_t *n_out)
{
    size_t cap = 256, n = 0;
    PfEntry *e = (PfEntry *)calloc(cap, sizeof(PfEntry));
    if (!e) return -1;

#define PF_PUSH(namestr)                                                     \
    do {                                                                     \
        if (n >= cap) {                                                      \
            cap *= 2;                                                        \
            PfEntry *g = (PfEntry *)realloc(e, cap * sizeof(PfEntry));       \
            if (!g) { free(e); return -1; }                                  \
            e = g;                                                           \
        }                                                                    \
        memset(&e[n], 0, sizeof(PfEntry));                                   \
        pf_basename_noext((namestr), e[n].id, sizeof(e[n].id));              \
        snprintf(e[n].path, sizeof(e[n].path), "%s/%s", dir, (namestr));     \
        if (sscanf(e[n].id, "z%d_y%d_x%d",                                   \
                   &e[n].vz, &e[n].vy, &e[n].vx) == 3)                       \
            n++;                                                             \
    } while (0)

#ifdef _WIN32
    char glob[1024];
    snprintf(glob, sizeof(glob), "%s/*.tif", dir);
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(glob, &fd);
    if (h == INVALID_HANDLE_VALUE) { free(e); return -1; }
    do {
        if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
        size_t nlen = strlen(fd.cFileName);
        if (nlen < 5 || strcmp(fd.cFileName + nlen - 4, ".tif") != 0) continue;
        PF_PUSH(fd.cFileName);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(dir);
    if (!d) { free(e); return -1; }
    struct dirent *de = NULL;
    while ((de = readdir(d)) != NULL) {
        size_t nlen = strlen(de->d_name);
        if (nlen < 5 || strcmp(de->d_name + nlen - 4, ".tif") != 0) continue;
        PF_PUSH(de->d_name);
    }
    closedir(d);
#endif
#undef PF_PUSH

    qsort(e, n, sizeof(PfEntry), pf_entry_cmp);
    *out = e;
    *n_out = n;
    return 0;
}

/* Byte-copy a file; returns 0 on success, -1 if src missing/unwritable. */
static int pf_copy_file(const char *src, const char *dst)
{
    FILE *fi = fopen(src, "rb");
    if (!fi) return -1;
    FILE *fo = fopen(dst, "wb");
    if (!fo) { fclose(fi); return -1; }
    char buf[65536];
    size_t r = 0;
    int rc = 0;
    while ((r = fread(buf, 1, sizeof(buf), fi)) > 0) {
        if (fwrite(buf, 1, r, fo) != r) { rc = -1; break; }
    }
    fclose(fi);
    if (fclose(fo) != 0) rc = -1;
    return rc;
}

/* ---------------------------------------------------------------- */
/* Options                                                          */
/* ---------------------------------------------------------------- */

typedef struct {
    const char *grid_dir;
    const char *out_dir;
    double umb_y, umb_x;     /* world L0; < 0 = not given */
    float  reach_max, reach_safe, radial_dr;
    float  min_score;
    int    min_support;
    int    paint_radius;
    int    no_tracks;
    int    cut_bridges;      /* erase certified persistent thin-neck welds */
    int    split_slabs;      /* geometric mid-plane carve of fused slabs
                              * (its OWN pass: joins/bridges never run) */
    float  slab_thick_r;     /* <= 0: SlabSplit default */
    float  slab_ext_r;       /* <= 0: SlabSplit default */
    int    slab_rounds;      /* <= 0: SlabSplit default */
    int    dry_run;
    int    png_all;
    int    png_scale;
    int    png_max;
    int    crop_max;         /* per-join magnified crop dumps */
    int    rej_crop_max;     /* sampled reject crop dumps */
    int    zlo, zhi;         /* local plane range; -1 = full */
    int    threads;
    int    tile;             /* y/x tile side in CUBES; 0 = auto (whole grid
                              * if it fits, else the largest side that does);
                              * tiles span FULL z so ConnectionTracks stay
                              * intact. Endpoints at tile edges fall in the
                              * existing border-exclusion band, so no join
                              * ever crosses a tile seam — the halo/weld
                              * machinery owns cube borders downstream. */
    int    quiet_viz;        /* selftest: skip PNG dumps */
} PfOpts;

/* One y/x tile of the cube grid (cube indices relative to the grid min). */
typedef struct {
    int  cy0, cy1, cx0, cx1; /* half-open cube ranges */
    int  append;             /* append to manifests instead of truncating */
    char tag[24];            /* viz filename prefix, "" for whole-grid */
} PfRegion;

/* Aggregated results for the grand summary in tiled runs. */
typedef struct {
    int64_t endpoints, candidates, joins, kept, far_kept, painted_px;
    int64_t bridges, certified, persistent, cut, cut_px;
    size_t  cubes_written;
} PfTotals;

static void pf_default_opts(PfOpts *o)
{
    memset(o, 0, sizeof(*o));
    o->umb_y = -1.0;
    o->umb_x = -1.0;
    o->reach_max = FIXUP_REACH_MAX_PX;
    o->reach_safe = FIXUP_REACH_SAFE_PX;
    o->radial_dr = FIXUP_RADIAL_DR_MAX;
    o->min_score = FIXUP_MIN_SCORE;
    o->min_support = -1;     /* auto: 0 when the radial gate is armed
                              * (support measured non-discriminating at 2.4um:
                              * dropped joins were 91.2% CONNECTED, and the
                              * filter-off mesh gate passed with components
                              * 44 -> 39), FIXUP_MIN_SUPPORT unarmed (then
                              * persistence IS a needed certificate). */
    o->paint_radius = FIXUP_PAINT_RADIUS;
    o->png_scale = 1;
    o->png_max = 256;
    o->crop_max = 200;
    o->rej_crop_max = 60;
    o->zlo = -1;
    o->zhi = -1;
    o->threads = 0;
}

/* ---------------------------------------------------------------- */
/* Per-plane results                                                */
/* ---------------------------------------------------------------- */

typedef struct {
    SliceTraceEndpoint *eps;
    int32_t             n_eps;
    int32_t             n_excluded;
    int32_t            *match;
    SliceMatchJoin     *joins;
    int32_t             n_joins;
    SliceMatchReject   *rejects;
    int32_t             n_rejects;
    SliceMatchStats     stats;
    /* final phase */
    int32_t            *support;      /* [n_joins] */
    int32_t            *conn_id;      /* [n_joins] */
    uint8_t            *kept;         /* [n_joins] */
    int32_t            *painted_px;   /* [n_joins] */
    /* thin-neck weld scan (detected on the ORIGINAL mask in round 1) */
    BridgeScanHit      *bridges;      /* [n_bridges] */
    int32_t             n_bridges;
    int32_t            *bridge_support; /* [n_bridges] planes with a nearby hit */
    int32_t            *bridge_cut_px;  /* [n_bridges] px erased (0 if not cut) */
} PfPlane;

static void *pf_persist(Arena_T pa, const void *src, size_t nbytes)
{
    if (nbytes == 0) return NULL;
    void *dst = ARENA_ALLOC(pa, nbytes);
    memcpy(dst, src, nbytes);
    return dst;
}

static const char *pf_reason_name(uint8_t reason)
{
    switch (reason) {
    case SLICEMATCH_REJ_SAME_CC:  return "same_cc";
    case SLICEMATCH_REJ_TANGENT:  return "tangent";
    case SLICEMATCH_REJ_RADIAL:   return "radial";
    case SLICEMATCH_REJ_EVIDENCE: return "evidence";
    case SLICEMATCH_REJ_SCORE:    return "score";
    case SLICEMATCH_REJ_ARC:      return "arc";
    case SLICEMATCH_REJ_MERGER:   return "merger";
    case SLICEMATCH_REJ_XSHEET:   return "cross_sheet";
    case SLICEMATCH_REJ_CROSSING: return "crossing";
    case SLICEMATCH_REJ_OCCUPIED: return "occupied";
    default:                      return "unknown";
    }
}

/* ---------------------------------------------------------------- */
/* The run                                                          */
/* ---------------------------------------------------------------- */

enum { PF_MAX_THREADS = 128 };
/* Whole-in-RAM ceiling. The dominant allocation is the int32 per-plane label
 * volume (4x vol) + vol + painted: ~6x vol total. 10x10x10 (1280^3 = 2.1 GB)
 * needs ~13 GB; 21x21x21 (2688^3 = 19.4 GB) needs ~120 GB — the dev box has
 * 256 GB, so whole-in-RAM covers the full ladder. The ceiling exists to fail
 * fast on anything larger (z-slab streaming would be the path there). */
#define PF_MAX_TOTAL_BYTES (200ull * 1024ull * 1024ull * 1024ull)

static int pf_run(const PfOpts *opts, const PfRegion *reg, PfTotals *tot)
{
    double t0 = pf_now();
    Arena_T main_arena = Arena_new();
    const char *tag = reg ? reg->tag : "";

    /* ---- scan the grid ---- */
    char pred_dir[1024];
    snprintf(pred_dir, sizeof(pred_dir), "%s/cubes_PRED", opts->grid_dir);
    PfEntry *entries = NULL;
    size_t n_entries = 0;
    if (pf_collect_cubes(pred_dir, &entries, &n_entries) != 0 ||
        n_entries == 0) {
        fprintf(stderr, "pred_fixup: no cubes in %s\n", pred_dir);
        Arena_dispose(&main_arena);
        free(entries);
        return 1;
    }


    /* chunk size from the first cube */
    int chunk = 0;
    {
        Arena_Mark mark = Arena_save(main_arena);
        uint8_t *cv = NULL;
        int cd = 0, ch = 0, cw = 0;
        if (TiffIO_load(main_arena, entries[0].path, &cv, &cd, &ch, &cw) != 0 ||
            cd != ch || ch != cw || cd <= 0) {
            fprintf(stderr, "pred_fixup: bad first cube %s (%dx%dx%d)\n",
                    entries[0].path, cd, ch, cw);
            Arena_dispose(&main_arena);
            free(entries);
            return 1;
        }
        chunk = cd;
        Arena_restore(main_arena, mark);
    }

    /* region filter: keep only this tile's cubes (ranges in cube units
     * relative to the grid minimum) */
    if (reg) {
        int gy0 = entries[0].vy, gx0 = entries[0].vx;
        for (size_t i = 0; i < n_entries; i++) {
            if (entries[i].vy < gy0) gy0 = entries[i].vy;
            if (entries[i].vx < gx0) gx0 = entries[i].vx;
        }
        size_t w = 0;
        for (size_t i = 0; i < n_entries; i++) {
            int cy = (entries[i].vy - gy0) / chunk;
            int cx = (entries[i].vx - gx0) / chunk;
            if (cy >= reg->cy0 && cy < reg->cy1 &&
                cx >= reg->cx0 && cx < reg->cx1)
                entries[w++] = entries[i];
        }
        n_entries = w;
        if (n_entries == 0) {
            Arena_dispose(&main_arena);
            free(entries);
            return 0;               /* empty tile (missing cubes) */
        }
    }

    int z0w = entries[0].vz, y0w = entries[0].vy, x0w = entries[0].vx;
    int z1w = z0w, y1w = y0w, x1w = x0w;
    for (size_t i = 0; i < n_entries; i++) {
        if (entries[i].vz < z0w) z0w = entries[i].vz;
        if (entries[i].vy < y0w) y0w = entries[i].vy;
        if (entries[i].vx < x0w) x0w = entries[i].vx;
        if (entries[i].vz > z1w) z1w = entries[i].vz;
        if (entries[i].vy > y1w) y1w = entries[i].vy;
        if (entries[i].vx > x1w) x1w = entries[i].vx;
    }
    int D = z1w - z0w + chunk;
    int H = y1w - y0w + chunk;
    int W = x1w - x0w + chunk;
    int nlz = D / chunk, nly = H / chunk, nlx = W / chunk;
    size_t plane_sz = (size_t)H * (size_t)W;
    size_t vol_sz = (size_t)D * plane_sz;

    if ((unsigned long long)vol_sz * 6ull > PF_MAX_TOTAL_BYTES) {
        fprintf(stderr,
                "pred_fixup: grid %dx%dx%d needs %.1f GB (vol+labels+paint); "
                "the z-slab streaming variant is required at this scale\n",
                D, H, W,
                (double)(vol_sz * 6) / (1024.0 * 1024.0 * 1024.0));
        Arena_dispose(&main_arena);
        free(entries);
        return 1;
    }

    fprintf(stderr,
            "pred_fixup: %zu cubes, chunk %d, grid %dx%dx%d "
            "(world z%d y%d x%d), vol %.0f MB\n",
            n_entries, chunk, D, H, W, z0w, y0w, x0w,
            (double)vol_sz / (1024.0 * 1024.0));

#ifdef _OPENMP
    if (opts->threads > 0) omp_set_num_threads(opts->threads);
#endif

    /* ---- output dirs ---- */
    char out_pred[1024], out_viz[1024];
    snprintf(out_pred, sizeof(out_pred), "%s/cubes_PRED", opts->out_dir);
    snprintf(out_viz, sizeof(out_viz), "%s/viz", opts->out_dir);
    ves_mkdir(opts->out_dir);
    if (!opts->dry_run) ves_mkdir(out_pred);
    if (!opts->quiet_viz) ves_mkdir(out_viz);

    /* ---- assemble the volume (thresholded to 0/1) ---- */
    uint8_t *vol = ARENA_CALLOC(main_arena, vol_sz, 1);
    uint8_t *present = ARENA_CALLOC(main_arena,
                                    (size_t)nlz * (size_t)nly * (size_t)nlx, 1);
    size_t fg_in = 0;
    for (size_t i = 0; i < n_entries; i++) {
        Arena_Mark mark = Arena_save(main_arena);
        uint8_t *cv = NULL;
        int cd = 0, chh = 0, cww = 0;
        if (TiffIO_load(main_arena, entries[i].path, &cv, &cd, &chh, &cww)
                != 0 || cd != chunk || chh != chunk || cww != chunk) {
            fprintf(stderr, "pred_fixup: skip bad cube %s\n", entries[i].id);
            Arena_restore(main_arena, mark);
            continue;
        }
        int lz = entries[i].vz - z0w, ly = entries[i].vy - y0w,
            lx = entries[i].vx - x0w;
        for (int z = 0; z < chunk; z++) {
            for (int y = 0; y < chunk; y++) {
                const uint8_t *src = cv + ((size_t)z * (size_t)chunk +
                                           (size_t)y) * (size_t)chunk;
                uint8_t *dst = vol + (size_t)(lz + z) * plane_sz +
                               (size_t)(ly + y) * (size_t)W + (size_t)lx;
                for (int x = 0; x < chunk; x++) {
                    uint8_t v = (uint8_t)(src[x] > 0 ? 1 : 0);
                    dst[x] = v;
                    fg_in += v;
                }
            }
        }
        present[((size_t)(lz / chunk) * (size_t)nly + (size_t)(ly / chunk)) *
                (size_t)nlx + (size_t)(lx / chunk)] = 1;
        Arena_restore(main_arena, mark);
    }

    /* absent-region layers (only if some grid cells are missing) */
    size_t n_cells = (size_t)nlz * (size_t)nly * (size_t)nlx;
    size_t n_present = 0;
    for (size_t i = 0; i < n_cells; i++) n_present += present[i];
    uint8_t **absent_layers = NULL;
    if (n_present < n_cells) {
        absent_layers = ARENA_CALLOC(main_arena, (size_t)nlz,
                                     sizeof(uint8_t *));
        for (int lz = 0; lz < nlz; lz++) {
            absent_layers[lz] = ARENA_CALLOC(main_arena, plane_sz, 1);
            for (int cy = 0; cy < nly; cy++)
                for (int cx = 0; cx < nlx; cx++) {
                    if (present[((size_t)lz * (size_t)nly + (size_t)cy) *
                                (size_t)nlx + (size_t)cx])
                        continue;
                    for (int y = cy * chunk; y < (cy + 1) * chunk; y++)
                        memset(absent_layers[lz] + (size_t)y * (size_t)W +
                               (size_t)(cx * chunk), 1, (size_t)chunk);
                }
        }
        fprintf(stderr, "pred_fixup: %zu/%zu cells present, absent mask armed\n",
                n_present, n_cells);
    }

    double t_load = pf_now();
    fprintf(stderr, "pred_fixup: loaded, fg %.2f%% (%.1fs)\n",
            100.0 * (double)fg_in / (double)vol_sz, t_load - t0);

    /* ---- per-plane CC labels for the whole volume ---- */
    int32_t *all_labels = ARENA_ALLOC(main_arena, vol_sz * sizeof(int32_t));
    Arena_T scratch_arenas[PF_MAX_THREADS];
    Arena_T persist_arenas[PF_MAX_THREADS];
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
    if (n_threads > PF_MAX_THREADS) n_threads = PF_MAX_THREADS;
#endif
    for (int t = 0; t < n_threads; t++) {
        scratch_arenas[t] = Arena_new();
        persist_arenas[t] = Arena_new();
    }

    {
        int z = 0;   /* MSVC OpenMP 2.0 (C mode): loop var declared outside */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for (z = 0; z < D; z++) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num() % PF_MAX_THREADS;
#endif
            SliceTrace_label_cc(scratch_arenas[tid],
                                vol + (size_t)z * plane_sz, H, W,
                                all_labels + (size_t)z * plane_sz);
        }
    }
    double t_label = pf_now();
    fprintf(stderr, "pred_fixup: labeled %d planes (%.1fs)\n", D,
            t_label - t_load);

    /* ---- plane range ---- */
    int zlo = (opts->zlo >= 0) ? opts->zlo : 0;
    int zhi = (opts->zhi >= 0 && opts->zhi <= D) ? opts->zhi : D;
    if (zlo < 0) zlo = 0;
    if (zlo > zhi) zlo = zhi;

    SliceMatchParams mp;
    memset(&mp, 0, sizeof(mp));
    mp.reach_safe = opts->reach_safe;
    mp.reach_max = opts->reach_max;
    mp.min_score = opts->min_score;
    mp.radial_dr_max = opts->radial_dr;
    if (opts->umb_y >= 0.0 && opts->umb_x >= 0.0) {
        mp.radial_armed = 1;
        mp.umb_py = (float)(opts->umb_y - (double)y0w);
        mp.umb_px = (float)(opts->umb_x - (double)x0w);
        fprintf(stderr,
                "pred_fixup: radial gate ARMED, umbilicus plane (%.1f, %.1f), "
                "dr_max %.1f\n",
                (double)mp.umb_py, (double)mp.umb_px, (double)mp.radial_dr_max);
    } else {
        fprintf(stderr,
                "pred_fixup: radial gate OFF (pass --umb-y/--umb-x to arm)\n");
    }
    int min_support = opts->min_support;
    if (min_support < 0)
        min_support = mp.radial_armed ? 0 : FIXUP_MIN_SUPPORT;
    fprintf(stderr, "pred_fixup: support filter min_support=%d%s\n",
            min_support, opts->min_support < 0 ? " (auto)" : "");

    PfPlane *r1 = ARENA_CALLOC(main_arena, (size_t)D, sizeof(PfPlane));
    PfPlane *rf = ARENA_CALLOC(main_arena, (size_t)D, sizeof(PfPlane));

    /* ---- round 1: trace + match ---- */
    {
    int r1_z = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (r1_z = zlo; r1_z < zhi; r1_z++) {
        int z = r1_z;
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num() % PF_MAX_THREADS;
#endif
        Arena_T sa = scratch_arenas[tid];
        Arena_T pa = persist_arenas[tid];
        Arena_Mark mark = Arena_save(sa);

        const uint8_t *absent =
            absent_layers ? absent_layers[z / chunk] : NULL;
        SliceTraceEndpoint *eps = NULL;
        int32_t n_eps = 0;
        uint8_t *skel = NULL;
        SliceTrace_run(sa, vol + (size_t)z * plane_sz, H, W, absent,
                       all_labels + (size_t)z * plane_sz, &skel,
                       &eps, &n_eps);

        /* thin-neck weld scan on the pristine mask (report; cut later) */
        BridgeScanHit *bhits = NULL;
        int32_t n_bhits = 0;
        BridgeScan_run(sa, vol + (size_t)z * plane_sz, skel, H, W,
                       mp.radial_armed, mp.umb_py, mp.umb_px,
                       &bhits, &n_bhits);
        r1[z].n_bridges = n_bhits;
        r1[z].bridges = pf_persist(pa, bhits,
                                   (size_t)n_bhits * sizeof(*bhits));

        int32_t *match = ARENA_ALLOC(sa, ((size_t)n_eps > 0 ? (size_t)n_eps
                                                            : 1) *
                                     sizeof(int32_t));
        SliceMatchJoin *joins = NULL;
        SliceMatchReject *rejects = NULL;
        int32_t n_joins = 0, n_rejects = 0;
        /* round 1 is provisional when a track round follows: its joins are
         * never painted, they only seed ConnectionTracks */
        SliceMatchParams mp1 = mp;
        mp1.provisional = opts->no_tracks ? 0 : 1;
        SliceMatch_run(sa, eps, n_eps, vol, all_labels, D, H, W, z,
                       all_labels + (size_t)z * plane_sz, NULL, &mp1, match,
                       &joins, &n_joins, &rejects, &n_rejects,
                       &r1[z].stats);

        r1[z].n_eps = n_eps;
        r1[z].eps = pf_persist(pa, eps, (size_t)n_eps * sizeof(*eps));
        r1[z].match = pf_persist(pa, match, (size_t)n_eps * sizeof(int32_t));
        r1[z].n_joins = n_joins;
        r1[z].joins = pf_persist(pa, joins, (size_t)n_joins * sizeof(*joins));
        r1[z].n_rejects = n_rejects;
        r1[z].rejects = pf_persist(pa, rejects,
                                   (size_t)n_rejects * sizeof(*rejects));
        for (int32_t i = 0; i < n_eps; i++)
            if (eps[i].excluded) r1[z].n_excluded++;

        Arena_restore(sa, mark);
    }
    }
    double t_r1 = pf_now();
    {
        int64_t sum_eps = 0, sum_joins = 0;
        for (int z = zlo; z < zhi; z++) {
            sum_eps += r1[z].n_eps;
            sum_joins += r1[z].n_joins;
        }
        fprintf(stderr,
                "pred_fixup: round 1 — %lld endpoints, %lld joins (%.1fs)\n",
                (long long)sum_eps, (long long)sum_joins, t_r1 - t_label);
    }

    /* ---- tracks + round 2 ---- */
    JoinTracks_T jt_final = NULL;
    if (!opts->no_tracks) {
        FixupPlaneView *views = ARENA_CALLOC(main_arena, (size_t)D,
                                             sizeof(FixupPlaneView));
        for (int z = 0; z < D; z++) {
            views[z].eps = r1[z].eps;
            views[z].match = r1[z].match;
            views[z].n_eps = r1[z].n_eps;
        }
        JoinTracks_T jt1 = JoinTracks_build(main_arena, views, D);
        fprintf(stderr,
                "pred_fixup: round-1 tracks — %d endpoint tracks, "
                "%d connection tracks\n",
                (int)JoinTracks_n_endpoint_tracks(jt1),
                (int)JoinTracks_n_connection_tracks(jt1));

        {
        int r2_z = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (r2_z = zlo; r2_z < zhi; r2_z++) {
            int z = r2_z;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num() % PF_MAX_THREADS;
#endif
            Arena_T sa = scratch_arenas[tid];
            Arena_T pa = persist_arenas[tid];
            Arena_Mark mark = Arena_save(sa);

            FixupPlaneView view;
            view.eps = r1[z].eps;
            view.match = r1[z].match;
            view.n_eps = r1[z].n_eps;
            float *s_track = JoinTracks_s_track_matrix(sa, jt1, z, &view);

            int32_t *match = ARENA_ALLOC(sa,
                                         ((size_t)view.n_eps > 0 ?
                                          (size_t)view.n_eps : 1) *
                                         sizeof(int32_t));
            SliceMatchJoin *joins = NULL;
            SliceMatchReject *rejects = NULL;
            int32_t n_joins = 0, n_rejects = 0;
            SliceMatch_run(sa, view.eps, view.n_eps, vol, all_labels, D, H, W,
                           z, all_labels + (size_t)z * plane_sz, s_track, &mp,
                           match, &joins, &n_joins, &rejects, &n_rejects,
                           &rf[z].stats);

            rf[z].eps = r1[z].eps;
            rf[z].n_eps = view.n_eps;
            rf[z].n_excluded = r1[z].n_excluded;
            rf[z].match = pf_persist(pa, match,
                                     (size_t)view.n_eps * sizeof(int32_t));
            rf[z].n_joins = n_joins;
            rf[z].joins = pf_persist(pa, joins,
                                     (size_t)n_joins * sizeof(*joins));
            rf[z].n_rejects = n_rejects;
            rf[z].rejects = pf_persist(pa, rejects,
                                       (size_t)n_rejects * sizeof(*rejects));
            rf[z].bridges = r1[z].bridges;
            rf[z].n_bridges = r1[z].n_bridges;

            Arena_restore(sa, mark);
        }
        }

        /* rebuild the registry on the round-2 matches: the support filter
         * must reflect the joins actually being considered */
        for (int z = 0; z < D; z++) {
            views[z].eps = rf[z].eps;
            views[z].match = rf[z].match;
            views[z].n_eps = rf[z].n_eps;
        }
        jt_final = JoinTracks_build(main_arena, views, D);
        double t_r2 = pf_now();
        int64_t sum_joins = 0;
        for (int z = zlo; z < zhi; z++) sum_joins += rf[z].n_joins;
        fprintf(stderr,
                "pred_fixup: round 2 — %lld joins; final tracks %d/%d (%.1fs)\n",
                (long long)sum_joins,
                (int)JoinTracks_n_endpoint_tracks(jt_final),
                (int)JoinTracks_n_connection_tracks(jt_final), t_r2 - t_r1);
    } else {
        for (int z = zlo; z < zhi; z++) rf[z] = r1[z];
        fprintf(stderr, "pred_fixup: --no-tracks: using round-1 joins\n");
    }

    /* ---- support filter + paint ---- */
    uint8_t *painted = ARENA_CALLOC(main_arena, vol_sz, 1);
    int64_t total_joins = 0, total_kept = 0, total_px = 0;
    int64_t support_dropped = 0;

    for (int z = zlo; z < zhi; z++) {
        PfPlane *p = &rf[z];
        if (p->n_joins == 0) continue;
        Arena_T pa = persist_arenas[0];
        p->support = ARENA_CALLOC(pa, (size_t)p->n_joins, sizeof(int32_t));
        p->conn_id = ARENA_CALLOC(pa, (size_t)p->n_joins, sizeof(int32_t));
        p->kept = ARENA_CALLOC(pa, (size_t)p->n_joins, 1);
        p->painted_px = ARENA_CALLOC(pa, (size_t)p->n_joins, sizeof(int32_t));
        for (int32_t k = 0; k < p->n_joins; k++) {
            if (jt_final) {
                p->support[k] = JoinTracks_pair_support(jt_final, z,
                                                        p->joins[k].a,
                                                        p->joins[k].b);
                p->conn_id[k] = JoinTracks_pair_conn_id(jt_final, z,
                                                        p->joins[k].a,
                                                        p->joins[k].b);
            } else {
                p->support[k] = min_support; /* --no-tracks: keep all */
                p->conn_id[k] = -1;
            }
            p->kept[k] = (uint8_t)(p->support[k] >= min_support);
            total_joins++;
            if (p->kept[k]) total_kept++;
            else support_dropped++;
        }
    }

    /* ---- thin-neck weld persistence + (opt-in) cutting ---- */
    uint8_t *cut_vol = NULL;
    int64_t total_bridges = 0, total_certified = 0, total_supported = 0;
    int64_t total_cut = 0, total_cut_px = 0;
    for (int z = zlo; z < zhi; z++) {
        PfPlane *p = &rf[z];
        if (p->n_bridges == 0) continue;
        Arena_T pa = persist_arenas[0];
        p->bridge_support = ARENA_CALLOC(pa, (size_t)p->n_bridges,
                                         sizeof(int32_t));
        p->bridge_cut_px = ARENA_CALLOC(pa, (size_t)p->n_bridges,
                                        sizeof(int32_t));
        for (int32_t k = 0; k < p->n_bridges; k++) {
            const BridgeScanHit *h = &p->bridges[k];
            int32_t sup = 0;
            for (int zz = z - 6; zz <= z + 6; zz++) {
                if (zz < zlo || zz >= zhi) continue;
                const PfPlane *q = &rf[zz];
                for (int32_t m = 0; m < q->n_bridges; m++) {
                    int ddy = q->bridges[m].my - h->my;
                    int ddx = q->bridges[m].mx - h->mx;
                    if (ddy * ddy + ddx * ddx <= 25) { sup++; break; }
                }
            }
            p->bridge_support[k] = sup;
            total_bridges++;
            if (h->certified) total_certified++;
            if (h->certified && sup >= FIXUP_BRIDGE_MIN_SUPPORT)
                total_supported++;
        }
    }
    if (opts->cut_bridges) {
        cut_vol = ARENA_CALLOC(main_arena, vol_sz, 1);
        for (int z = zlo; z < zhi; z++) {
            PfPlane *p = &rf[z];
            for (int32_t k = 0; k < p->n_bridges; k++) {
                const BridgeScanHit *h = &p->bridges[k];
                if (!h->certified ||
                    p->bridge_support[k] < FIXUP_BRIDGE_MIN_SUPPORT)
                    continue;
                p->bridge_cut_px[k] =
                    BridgeScan_cut(vol + (size_t)z * plane_sz,
                                   cut_vol + (size_t)z * plane_sz,
                                   H, W, NULL, h);
                total_cut++;
                total_cut_px += p->bridge_cut_px[k];
            }
        }
    }
    fprintf(stderr,
            "pred_fixup: bridges — %lld necks, %lld radial-certified, "
            "%lld persistent (support>=%d)%s\n",
            (long long)total_bridges, (long long)total_certified,
            (long long)total_supported, FIXUP_BRIDGE_MIN_SUPPORT,
            opts->cut_bridges ? " [CUTTING]" : " [report only]");
    if (opts->cut_bridges)
        fprintf(stderr, "pred_fixup: cut %lld welds, erased %lld px\n",
                (long long)total_cut, (long long)total_cut_px);

    {
        int pz = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4)
#endif
        for (pz = zlo; pz < zhi; pz++) {
            PfPlane *p = &rf[pz];
            for (int32_t k = 0; k < p->n_joins; k++) {
                if (!p->kept[k]) continue;
                p->painted_px[k] =
                    JoinPaint_draw(vol + (size_t)pz * plane_sz,
                                   painted + (size_t)pz * plane_sz, H, W,
                                   &p->eps[p->joins[k].a],
                                   &p->eps[p->joins[k].b],
                                   opts->paint_radius);
            }
        }
    }
    for (int z = zlo; z < zhi; z++)
        for (int32_t k = 0; k < rf[z].n_joins; k++)
            if (rf[z].kept && rf[z].kept[k]) total_px += rf[z].painted_px[k];
    double t_paint = pf_now();
    fprintf(stderr,
            "pred_fixup: kept %lld/%lld joins (%lld dropped by support<%d), "
            "painted %lld px\n",
            (long long)total_kept, (long long)total_joins,
            (long long)support_dropped, min_support,
            (long long)total_px);

    /* ---- manifests ---- */
    {
        const char *fmode = (reg && reg->append) ? "a" : "w";
        char path[1024];
        snprintf(path, sizeof(path), "%s/joins.jsonl", opts->out_dir);
        FILE *fj = fopen(path, fmode);
        snprintf(path, sizeof(path), "%s/rejected.jsonl", opts->out_dir);
        FILE *fr = fopen(path, fmode);
        snprintf(path, sizeof(path), "%s/bridges.jsonl", opts->out_dir);
        FILE *fb = fopen(path, fmode);
        snprintf(path, sizeof(path), "%s/planes.csv", opts->out_dir);
        FILE *fc = fopen(path, fmode);
        if (fc && !(reg && reg->append))
            fprintf(fc, "z_world,n_eps,n_excluded,in_reach,candidates,"
                        "joins,kept,painted_px,rej_same_cc,rej_tangent,"
                        "rej_radial,rej_evidence,rej_score,rej_arc,"
                        "rej_merger,rej_cross_sheet,rej_crossing,"
                        "rej_occupied\n");
        int64_t rej_rows = 0;
        enum { PF_MAX_REJ_ROWS = 200000 };
        for (int z = zlo; z < zhi; z++) {
            PfPlane *p = &rf[z];
            int zw = z0w + z;
            if (fj) {
                for (int32_t k = 0; k < p->n_joins; k++) {
                    const SliceMatchJoin *jn = &p->joins[k];
                    fprintf(fj,
                        "{\"z\":%d,\"a\":{\"y\":%d,\"x\":%d},"
                        "\"b\":{\"y\":%d,\"x\":%d},\"dist\":%.2f,"
                        "\"score\":%.3f,\"s_adj\":%.3f,\"s_track\":%.3f,"
                        "\"far_tier\":%d,\"support\":%d,\"conn_id\":%d,"
                        "\"kept\":%d,\"painted_px\":%d}\n",
                        zw,
                        y0w + (int)p->eps[jn->a].y, x0w + (int)p->eps[jn->a].x,
                        y0w + (int)p->eps[jn->b].y, x0w + (int)p->eps[jn->b].x,
                        (double)jn->dist, (double)jn->score,
                        (double)jn->s_adj, (double)jn->s_track,
                        (int)jn->far_tier,
                        p->support ? (int)p->support[k] : 0,
                        p->conn_id ? (int)p->conn_id[k] : -1,
                        p->kept ? (int)p->kept[k] : 0,
                        p->painted_px ? (int)p->painted_px[k] : 0);
                }
            }
            if (fr && rej_rows < PF_MAX_REJ_ROWS) {
                for (int32_t k = 0; k < p->n_rejects &&
                                    rej_rows < PF_MAX_REJ_ROWS; k++) {
                    const SliceMatchReject *rj = &p->rejects[k];
                    fprintf(fr,
                        "{\"z\":%d,\"a\":{\"y\":%d,\"x\":%d},"
                        "\"b\":{\"y\":%d,\"x\":%d},\"reason\":\"%s\"}\n",
                        zw,
                        y0w + (int)p->eps[rj->a].y, x0w + (int)p->eps[rj->a].x,
                        y0w + (int)p->eps[rj->b].y, x0w + (int)p->eps[rj->b].x,
                        pf_reason_name(rj->reason));
                    rej_rows++;
                }
            }
            if (fb) {
                for (int32_t k = 0; k < p->n_bridges; k++) {
                    const BridgeScanHit *h = &p->bridges[k];
                    fprintf(fb,
                        "{\"z\":%d,\"a\":{\"y\":%d,\"x\":%d},"
                        "\"b\":{\"y\":%d,\"x\":%d},\"len\":%d,"
                        "\"width\":%.2f,\"dr\":%.2f,\"radial_dot\":%.2f,"
                        "\"certified\":%d,\"support\":%d,\"cut_px\":%d}\n",
                        zw,
                        y0w + (int)h->ay, x0w + (int)h->ax,
                        y0w + (int)h->by, x0w + (int)h->bx,
                        (int)h->len, (double)h->width, (double)h->dr,
                        (double)h->radial_dot, (int)h->certified,
                        p->bridge_support ? (int)p->bridge_support[k] : 0,
                        p->bridge_cut_px ? (int)p->bridge_cut_px[k] : 0);
                }
            }
            if (fc) {
                int32_t kept = 0, px = 0;
                for (int32_t k = 0; k < p->n_joins; k++)
                    if (p->kept && p->kept[k]) {
                        kept++;
                        px += p->painted_px[k];
                    }
                fprintf(fc,
                        "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,"
                        "%d,%d\n",
                        zw, (int)p->n_eps, (int)p->n_excluded,
                        (int)p->stats.n_pairs_in_reach,
                        (int)p->stats.n_candidates, (int)p->n_joins,
                        (int)kept, (int)px,
                        (int)p->stats.rej_count[SLICEMATCH_REJ_SAME_CC],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_TANGENT],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_RADIAL],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_EVIDENCE],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_SCORE],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_ARC],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_MERGER],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_XSHEET],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_CROSSING],
                        (int)p->stats.rej_count[SLICEMATCH_REJ_OCCUPIED]);
            }
        }
        if (fj) fclose(fj);
        if (fr) fclose(fr);
        if (fb) fclose(fb);
        if (fc) fclose(fc);
    }

    /* ---- write the fixed grid ---- */
    size_t cubes_written = 0;
    if (!opts->dry_run) {
        int write_err = 0;
        int wi = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (wi = 0; wi < (int)n_entries; wi++) {
            int i = wi;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num() % PF_MAX_THREADS;
#endif
            Arena_T sa = scratch_arenas[tid];
            Arena_Mark mark = Arena_save(sa);
            uint8_t *cube = ARENA_ALLOC(sa, (size_t)chunk * (size_t)chunk *
                                        (size_t)chunk);
            int lz = entries[i].vz - z0w, ly = entries[i].vy - y0w,
                lx = entries[i].vx - x0w;
            for (int z = 0; z < chunk; z++) {
                for (int y = 0; y < chunk; y++) {
                    const uint8_t *src = vol + (size_t)(lz + z) * plane_sz +
                                         (size_t)(ly + y) * (size_t)W +
                                         (size_t)lx;
                    uint8_t *dst = cube + ((size_t)z * (size_t)chunk +
                                           (size_t)y) * (size_t)chunk;
                    for (int x = 0; x < chunk; x++)
                        dst[x] = (uint8_t)(src[x] ? 255 : 0);
                }
            }
            char path[1024];
            snprintf(path, sizeof(path), "%s/%s.tif", out_pred,
                     entries[i].id);
            if (TiffIO_save(path, cube, chunk, chunk, chunk) != 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    fprintf(stderr, "pred_fixup: FAILED to write %s\n", path);
                    write_err = 1;
                }
            }
            Arena_restore(sa, mark);
        }
        if (write_err) {
            fprintf(stderr, "pred_fixup: output grid INCOMPLETE\n");
        } else {
            cubes_written = n_entries;
        }

        /* manifest.json + present.json travel with the grid */
        char src[1024], dst[1024];
        snprintf(src, sizeof(src), "%s/manifest.json", opts->grid_dir);
        snprintf(dst, sizeof(dst), "%s/manifest.json", opts->out_dir);
        pf_copy_file(src, dst);
        snprintf(src, sizeof(src), "%s/cubes_PRED/present.json",
                 opts->grid_dir);
        snprintf(dst, sizeof(dst), "%s/cubes_PRED/present.json",
                 opts->out_dir);
        pf_copy_file(src, dst);
    }
    double t_write = pf_now();

    /* ---- viz ---- */
    int n_overlays = 0, n_crops = 0, n_rej_crops = 0, n_bridge_crops = 0;
    if (!opts->quiet_viz) {
        Arena_T va = scratch_arenas[0];
        char joins_dir[1024], rej_dir[1024], bridge_dir[1024];
        snprintf(joins_dir, sizeof(joins_dir), "%s/joins", out_viz);
        snprintf(rej_dir, sizeof(rej_dir), "%s/rejects", out_viz);
        snprintf(bridge_dir, sizeof(bridge_dir), "%s/bridges", out_viz);
        ves_mkdir(joins_dir);
        ves_mkdir(rej_dir);
        ves_mkdir(bridge_dir);
        for (int z = zlo; z < zhi; z++) {
            PfPlane *p = &rf[z];
            int has_kept = 0;
            for (int32_t k = 0; k < p->n_joins; k++)
                if (p->kept && p->kept[k]) { has_kept = 1; break; }
            int has_hot_bridge = 0;
            for (int32_t k = 0; k < p->n_bridges; k++)
                if (p->bridges[k].certified && p->bridge_support &&
                    p->bridge_support[k] >= FIXUP_BRIDGE_MIN_SUPPORT) {
                    has_hot_bridge = 1;
                    break;
                }
            int want_crops =
                (has_kept && (opts->crop_max <= 0 ||
                              n_crops < opts->crop_max)) ||
                (has_hot_bridge && n_bridge_crops < 60) ||
                (p->n_rejects > 0 && (opts->rej_crop_max <= 0 ||
                                      n_rej_crops < opts->rej_crop_max));
            int want = opts->png_all || p->stats.n_candidates > 0 ||
                       p->n_joins > 0;
            int want_overlay = want &&
                (opts->png_max <= 0 || n_overlays < opts->png_max);
            if (!want_overlay && !want_crops) continue;

            Arena_Mark mark = Arena_save(va);
            const uint8_t *absent =
                absent_layers ? absent_layers[z / chunk] : NULL;
            /* re-trace for the skeleton layer (deterministic, cheap);
             * the plane now contains paint, so trace the pre-paint mask */
            uint8_t *pre = ARENA_ALLOC(va, plane_sz);
            const uint8_t *pl = vol + (size_t)z * plane_sz;
            const uint8_t *pp = painted + (size_t)z * plane_sz;
            for (size_t i = 0; i < plane_sz; i++)
                pre[i] = (uint8_t)(pl[i] && !pp[i]);
            uint8_t *skel = NULL;
            SliceTraceEndpoint *eps2 = NULL;
            int32_t n2 = 0;
            SliceTrace_run(va, pre, H, W, absent,
                           all_labels + (size_t)z * plane_sz, &skel,
                           &eps2, &n2);

            char path[1024];
            int zw = z0w + z;
            if (want_overlay) {
                snprintf(path, sizeof(path), "%s/%sz%05d_overlay.png", out_viz,
                         tag, zw);
                FixupViz_plane_overlay(path, va, pl, pp, skel, H, W,
                                       p->eps, p->n_eps, p->joins, p->n_joins,
                                       p->support, min_support,
                                       p->rejects, p->n_rejects,
                                       opts->png_scale);
                n_overlays++;
                if (has_kept) {
                    char pb[1024], pa2[1024];
                    snprintf(pb, sizeof(pb), "%s/%sz%05d_before.png", out_viz,
                             tag, zw);
                    snprintf(pa2, sizeof(pa2), "%s/%sz%05d_after.png", out_viz,
                             tag, zw);
                    FixupViz_before_after(pb, pa2, va, pl, pp, H, W);
                }
            }

            /* magnified per-join crops — the viewer that shows a connection */
            for (int32_t k = 0; k < p->n_joins; k++) {
                if (!(p->kept && p->kept[k])) continue;
                if (opts->crop_max > 0 && n_crops >= opts->crop_max) break;
                snprintf(path, sizeof(path), "%s/joins/%sz%05d_j%02d.png",
                         out_viz, tag, zw, (int)k);
                FixupViz_pair_crop(path, va, pl, pp, skel, H, W,
                                   p->eps, p->n_eps, p->joins[k].a,
                                   p->joins[k].b, 1, 0, 255, 90,
                                   4, 32);
                n_crops++;
            }

            /* sampled reject crops (the gate-tuning views) */
            for (int32_t k = 0; k < p->n_rejects; k++) {
                if (opts->rej_crop_max > 0 &&
                    n_rej_crops >= opts->rej_crop_max) break;
                uint8_t rr = p->rejects[k].reason;
                if (rr == SLICEMATCH_REJ_SAME_CC ||
                    rr == SLICEMATCH_REJ_OCCUPIED)
                    continue;
                uint8_t lr2 = 255, lg2 = 150, lb2 = 40;
                if (rr == SLICEMATCH_REJ_MERGER ||
                    rr == SLICEMATCH_REJ_XSHEET ||
                    rr == SLICEMATCH_REJ_RADIAL) {
                    lr2 = 255; lg2 = 60; lb2 = 60;
                } else if (rr == SLICEMATCH_REJ_CROSSING) {
                    lr2 = 0; lg2 = 180; lb2 = 180;
                }
                snprintf(path, sizeof(path),
                         "%s/rejects/%sz%05d_r%02d_%s.png", out_viz, tag, zw,
                         (int)k, pf_reason_name(rr));
                FixupViz_pair_crop(path, va, pl, pp, skel, H, W,
                                   p->eps, p->n_eps, p->rejects[k].a,
                                   p->rejects[k].b, 0, lr2, lg2, lb2,
                                   4, 32);
                n_rej_crops++;
            }

            /* certified persistent thin-neck welds (red line) */
            for (int32_t k = 0; k < p->n_bridges &&
                                n_bridge_crops < 60; k++) {
                const BridgeScanHit *h = &p->bridges[k];
                if (!h->certified || !p->bridge_support ||
                    p->bridge_support[k] < FIXUP_BRIDGE_MIN_SUPPORT)
                    continue;
                SliceTraceEndpoint bp[2];
                memset(bp, 0, sizeof(bp));
                bp[0].y = h->ay;
                bp[0].x = h->ax;
                bp[1].y = h->by;
                bp[1].x = h->bx;
                snprintf(path, sizeof(path), "%s/bridges/%sz%05d_b%02d.png",
                         out_viz, tag, zw, (int)k);
                FixupViz_pair_crop(path, va, pl, pp, skel, H, W, bp, 2,
                                   0, 1, 0, 255, 60, 60, 4, 32);
                n_bridge_crops++;
            }
            Arena_restore(va, mark);
        }
        char path[1024];
        snprintf(path, sizeof(path), "%s/%ssummary_maxproj.png", out_viz, tag);
        FixupViz_grid_summary(path, scratch_arenas[0], vol, painted, D, H, W);
        fprintf(stderr,
                "pred_fixup: wrote %d overlays, %d join crops, %d reject "
                "crops, %d bridge crops + summary -> %s\n",
                n_overlays, n_crops, n_rej_crops, n_bridge_crops, out_viz);
    }

    /* ---- report (whole-grid runs only; tiled runs get the wrapper's
     * combined report) ---- */
    if (!reg) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/fixup_report.json", opts->out_dir);
        FILE *f = fopen(path, "w");
        if (f) {
            int64_t sum_eps = 0, sum_cand = 0;
            int64_t rejs[SLICEMATCH_REJ_N_REASONS];
            memset(rejs, 0, sizeof(rejs));
            for (int z = zlo; z < zhi; z++) {
                sum_eps += rf[z].n_eps;
                sum_cand += rf[z].stats.n_candidates;
                for (int rr = 0; rr < SLICEMATCH_REJ_N_REASONS; rr++)
                    rejs[rr] += rf[z].stats.rej_count[rr];
            }
            fprintf(f, "{\n");
            fprintf(f, "  \"grid\": \"%s\",\n", opts->grid_dir);
            fprintf(f, "  \"dims\": [%d, %d, %d],\n", D, H, W);
            fprintf(f, "  \"world_origin\": [%d, %d, %d],\n", z0w, y0w, x0w);
            fprintf(f, "  \"cubes\": %zu,\n", n_entries);
            fprintf(f, "  \"params\": {\"reach_safe\": %.1f, \"reach_max\": "
                       "%.1f, \"radial_armed\": %d, \"radial_dr\": %.1f, "
                       "\"min_support\": %d, \"min_score\": %.2f, "
                       "\"paint_radius\": %d, \"tracks\": %d},\n",
                    (double)opts->reach_safe, (double)opts->reach_max,
                    mp.radial_armed, (double)opts->radial_dr,
                    min_support, (double)opts->min_score,
                    opts->paint_radius, !opts->no_tracks);
            fprintf(f, "  \"endpoints\": %lld,\n", (long long)sum_eps);
            fprintf(f, "  \"candidates\": %lld,\n", (long long)sum_cand);
            fprintf(f, "  \"joins\": %lld,\n", (long long)total_joins);
            fprintf(f, "  \"kept\": %lld,\n", (long long)total_kept);
            fprintf(f, "  \"support_dropped\": %lld,\n",
                    (long long)support_dropped);
            fprintf(f, "  \"painted_px\": %lld,\n", (long long)total_px);
            fprintf(f, "  \"bridges\": {\"necks\": %lld, \"certified\": %lld, "
                       "\"persistent\": %lld, \"cut\": %lld, "
                       "\"cut_px\": %lld},\n",
                    (long long)total_bridges, (long long)total_certified,
                    (long long)total_supported, (long long)total_cut,
                    (long long)total_cut_px);
            fprintf(f, "  \"fg_in\": %zu,\n", fg_in);
            fprintf(f, "  \"rejects\": {");
            for (int rr = 1; rr < SLICEMATCH_REJ_N_REASONS; rr++)
                fprintf(f, "%s\"%s\": %lld", rr > 1 ? ", " : "",
                        pf_reason_name((uint8_t)rr), (long long)rejs[rr]);
            fprintf(f, "},\n");
            fprintf(f, "  \"cubes_written\": %zu,\n", cubes_written);
            fprintf(f, "  \"seconds\": %.1f\n", pf_now() - t0);
            fprintf(f, "}\n");
            fclose(f);
        }
    }

    /* ---- summary block ---- */
    {
        int64_t sum_eps = 0, sum_cand = 0, sum_far = 0;
        for (int z = zlo; z < zhi; z++) {
            sum_eps += rf[z].n_eps;
            sum_cand += rf[z].stats.n_candidates;
            for (int32_t k = 0; k < rf[z].n_joins; k++)
                if (rf[z].kept && rf[z].kept[k] && rf[z].joins[k].far_tier)
                    sum_far++;
        }
        if (tot) {
            tot->endpoints += sum_eps;
            tot->candidates += sum_cand;
            tot->joins += total_joins;
            tot->kept += total_kept;
            tot->far_kept += sum_far;
            tot->painted_px += total_px;
            tot->bridges += total_bridges;
            tot->certified += total_certified;
            tot->persistent += total_supported;
            tot->cut += total_cut;
            tot->cut_px += total_cut_px;
            tot->cubes_written += cubes_written;
        }
        fprintf(stderr,
                "\n=== pred_fixup summary %s===\n"
                "  grid          %dx%dx%d (%zu cubes, chunk %d)\n"
                "  planes        [%d, %d) of %d\n"
                "  endpoints     %lld\n"
                "  candidates    %lld\n"
                "  joins kept    %lld / %lld (%lld far-tier, %lld dropped "
                "by support<%d)\n"
                "  painted px    %lld (%.4f%% of fg)\n"
                "  weld necks    %lld found, %lld radial-certified, %lld "
                "persistent, %lld cut (%lld px)\n"
                "  output        %s%s\n"
                "  total         %.1fs (load %.1f, label %.1f, match %.1f, "
                "paint %.1f, write %.1f)\n",
                tag, D, H, W, n_entries, chunk, zlo, zhi, D,
                (long long)sum_eps, (long long)sum_cand,
                (long long)total_kept, (long long)total_joins,
                (long long)sum_far, (long long)support_dropped,
                min_support, (long long)total_px,
                fg_in > 0 ? 100.0 * (double)total_px / (double)fg_in : 0.0,
                (long long)total_bridges, (long long)total_certified,
                (long long)total_supported, (long long)total_cut,
                (long long)total_cut_px,
                opts->dry_run ? "(dry run) " : "", opts->out_dir,
                pf_now() - t0, t_load - t0, t_label - t_load, t_paint - t_label,
                t_paint - t_r1 > 0 ? t_paint - t_r1 : 0.0, t_write - t_paint);
    }

    for (int t = 0; t < n_threads; t++) {
        Arena_dispose(&scratch_arenas[t]);
        Arena_dispose(&persist_arenas[t]);
    }
    Arena_dispose(&main_arena);
    free(entries);
    return 0;
}

/* ---------------------------------------------------------------- */
/* Tiling wrapper                                                   */
/* ---------------------------------------------------------------- */

/* Whole grid when it fits; else y/x tiles (full z) sized to the budget.
 * Tile seams are no-join zones by construction (border exclusion). */
/* ---------------------------------------------------------------- */
/* --split-slabs: geometric mid-plane carve of pred-fused wrap slabs.
 * Deliberately its OWN pass (assemble -> carve -> viz -> write) so the
 * join/bridge machinery is never entered: joins could re-bridge a fresh
 * carve, and the carve invalidates the plane labels the joins rely on.
 * Chain two invocations to combine (split grid -> fixup on it).        */
/* ---------------------------------------------------------------- */

static int pf_run_split(const PfOpts *opts)
{
    double t0 = pf_now();
    Arena_T arena = Arena_new();

    char pred_dir[1024];
    snprintf(pred_dir, sizeof(pred_dir), "%s/cubes_PRED", opts->grid_dir);
    PfEntry *entries = NULL;
    size_t n_entries = 0;
    if (pf_collect_cubes(pred_dir, &entries, &n_entries) != 0 ||
        n_entries == 0) {
        fprintf(stderr, "pred_fixup: no cubes in %s\n", pred_dir);
        Arena_dispose(&arena);
        free(entries);
        return 1;
    }

    int chunk = 0;
    {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *cv = NULL;
        int cd = 0, ch = 0, cw = 0;
        if (TiffIO_load(arena, entries[0].path, &cv, &cd, &ch, &cw) != 0 ||
            cd != ch || ch != cw || cd <= 0) {
            fprintf(stderr, "pred_fixup: bad first cube %s\n",
                    entries[0].path);
            Arena_dispose(&arena);
            free(entries);
            return 1;
        }
        chunk = cd;
        Arena_restore(arena, mark);
    }

    int z0w = entries[0].vz, y0w = entries[0].vy, x0w = entries[0].vx;
    int z1w = z0w, y1w = y0w, x1w = x0w;
    for (size_t i = 0; i < n_entries; i++) {
        if (entries[i].vz < z0w) z0w = entries[i].vz;
        if (entries[i].vy < y0w) y0w = entries[i].vy;
        if (entries[i].vx < x0w) x0w = entries[i].vx;
        if (entries[i].vz > z1w) z1w = entries[i].vz;
        if (entries[i].vy > y1w) y1w = entries[i].vy;
        if (entries[i].vx > x1w) x1w = entries[i].vx;
    }
    int D = z1w - z0w + chunk;
    int H = y1w - y0w + chunk;
    int W = x1w - x0w + chunk;
    size_t plane_sz = (size_t)H * (size_t)W;
    size_t vol_sz = (size_t)D * plane_sz;
    /* vol + carve mask + SlabSplit's cand + dil = 4x */
    if ((unsigned long long)vol_sz * 4ull > PF_MAX_TOTAL_BYTES) {
        fprintf(stderr,
                "pred_fixup: grid %dx%dx%d needs %.1f GB for --split-slabs; "
                "a z-slab streaming variant is required at this scale\n",
                D, H, W,
                (double)(vol_sz * 4) / (1024.0 * 1024.0 * 1024.0));
        Arena_dispose(&arena);
        free(entries);
        return 1;
    }
    fprintf(stderr,
            "pred_fixup: --split-slabs, %zu cubes, grid %dx%dx%d "
            "(world z%d y%d x%d), vol %.0f MB\n",
            n_entries, D, H, W, z0w, y0w, x0w,
            (double)vol_sz / (1024.0 * 1024.0));

#ifdef _OPENMP
    if (opts->threads > 0) omp_set_num_threads(opts->threads);
#endif

    char out_pred[1024], out_viz[1024];
    snprintf(out_pred, sizeof(out_pred), "%s/cubes_PRED", opts->out_dir);
    snprintf(out_viz, sizeof(out_viz), "%s/viz", opts->out_dir);
    ves_mkdir(opts->out_dir);
    if (!opts->dry_run) ves_mkdir(out_pred);
    if (!opts->quiet_viz) ves_mkdir(out_viz);

    /* ---- assemble ---- */
    uint8_t *vol = ARENA_CALLOC(arena, vol_sz, 1);
    size_t fg_in = 0;
    for (size_t i = 0; i < n_entries; i++) {
        Arena_Mark mark = Arena_save(arena);
        uint8_t *cv = NULL;
        int cd = 0, chh = 0, cww = 0;
        if (TiffIO_load(arena, entries[i].path, &cv, &cd, &chh, &cww) != 0 ||
            cd != chunk || chh != chunk || cww != chunk) {
            fprintf(stderr, "pred_fixup: skip bad cube %s\n", entries[i].id);
            Arena_restore(arena, mark);
            continue;
        }
        int lz = entries[i].vz - z0w, ly = entries[i].vy - y0w,
            lx = entries[i].vx - x0w;
        for (int z = 0; z < chunk; z++) {
            for (int y = 0; y < chunk; y++) {
                const uint8_t *src = cv + ((size_t)z * (size_t)chunk +
                                           (size_t)y) * (size_t)chunk;
                uint8_t *dst = vol + (size_t)(lz + z) * plane_sz +
                               (size_t)(ly + y) * (size_t)W + (size_t)lx;
                for (int x = 0; x < chunk; x++) {
                    uint8_t v = (uint8_t)(src[x] > 0 ? 1 : 0);
                    dst[x] = v;
                    fg_in += v;
                }
            }
        }
        Arena_restore(arena, mark);
    }
    double t_load = pf_now();
    fprintf(stderr, "pred_fixup: loaded, fg %.2f%% (%.1fs)\n",
            100.0 * (double)fg_in / (double)vol_sz, t_load - t0);

    /* ---- carve ---- */
    SlabSplitParams sp;
    SlabSplit_params_default(&sp);
    if (opts->slab_thick_r > 0.0f) sp.thick_r = opts->slab_thick_r;
    if (opts->slab_ext_r > 0.0f) sp.ext_r = opts->slab_ext_r;
    if (opts->slab_rounds > 0) sp.rounds = opts->slab_rounds;
    fprintf(stderr,
            "pred_fixup: slab split thick_r=%.2f ext_r=%.2f grow_max=%d "
            "rounds=%d z_support=%d/%d cap=%.2f\n",
            (double)sp.thick_r, (double)sp.ext_r, sp.grow_max, sp.rounds,
            sp.z_support, sp.z_window, (double)sp.plane_cap_frac);
    uint8_t *carve_mask = ARENA_CALLOC(arena, vol_sz, 1);
    SlabSplitStats st;
    if (SlabSplit_run(arena, vol, carve_mask, D, H, W, &sp, &st) != 0) {
        fprintf(stderr, "pred_fixup: slab split FAILED (bad params)\n");
        Arena_dispose(&arena);
        free(entries);
        return 2;
    }
    double t_carve = pf_now();
    fprintf(stderr,
            "pred_fixup: slab split — primary %lld px, supported %lld px, "
            "carved %lld px (%.3f%% of fg) on %lld planes, %lld capped, "
            "%d round(s) (%.1fs)\n",
            (long long)st.primary_px, (long long)st.supported_px,
            (long long)st.carved_px,
            fg_in ? 100.0 * (double)st.carved_px / (double)fg_in : 0.0,
            (long long)st.planes_carved, (long long)st.planes_capped,
            st.rounds_run, t_carve - t_load);

    /* ---- per-plane manifest + viz ---- */
    {
        char path[1024];
        snprintf(path, sizeof(path), "%s/slabs.csv", opts->out_dir);
        FILE *fc = fopen(path, "w");
        if (fc) fprintf(fc, "z_world,carved_px\n");
        int n_png = 0;
        for (int z = 0; z < D; z++) {
            const uint8_t *cm = carve_mask + (size_t)z * plane_sz;
            size_t carved = 0;
            for (size_t i = 0; i < plane_sz; i++) carved += cm[i];
            if (carved == 0) continue;
            if (fc) fprintf(fc, "%d,%zu\n", z0w + z, carved);
            if (!opts->quiet_viz &&
                (opts->png_all || n_png < opts->png_max)) {
                Arena_Mark mark = Arena_save(arena);
                uint8_t *rgb = ARENA_ALLOC(arena, plane_sz * 3);
                const uint8_t *pl = vol + (size_t)z * plane_sz;
                for (size_t i = 0; i < plane_sz; i++) {
                    uint8_t g = pl[i] ? 160 : 0;
                    rgb[i * 3 + 0] = cm[i] ? 255 : g;
                    rgb[i * 3 + 1] = cm[i] ? 32 : g;
                    rgb[i * 3 + 2] = cm[i] ? 32 : g;
                }
                snprintf(path, sizeof(path), "%s/slab_z%05d.png", out_viz,
                         z0w + z);
                if (VesPng_write_rgb(path, rgb, W, H) == 0) n_png++;
                Arena_restore(arena, mark);
            }
        }
        if (fc) fclose(fc);
        fprintf(stderr, "pred_fixup: %d slab overlay PNG(s) in %s\n",
                n_png, out_viz);
    }

    /* ---- write the split grid ---- */
    size_t cubes_written = 0;
    if (!opts->dry_run) {
        int write_err = 0;
        Arena_T warenas[PF_MAX_THREADS];
        int n_threads = 1;
#ifdef _OPENMP
        n_threads = omp_get_max_threads();
        if (n_threads > PF_MAX_THREADS) n_threads = PF_MAX_THREADS;
#endif
        for (int t = 0; t < n_threads; t++) warenas[t] = Arena_new();
        int wi = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (wi = 0; wi < (int)n_entries; wi++) {
            int i = wi;
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num() % PF_MAX_THREADS;
#endif
            Arena_Mark mark = Arena_save(warenas[tid]);
            uint8_t *cube = ARENA_ALLOC(warenas[tid],
                                        (size_t)chunk * (size_t)chunk *
                                        (size_t)chunk);
            int lz = entries[i].vz - z0w, ly = entries[i].vy - y0w,
                lx = entries[i].vx - x0w;
            for (int z = 0; z < chunk; z++) {
                for (int y = 0; y < chunk; y++) {
                    const uint8_t *src = vol + (size_t)(lz + z) * plane_sz +
                                         (size_t)(ly + y) * (size_t)W +
                                         (size_t)lx;
                    uint8_t *dst = cube + ((size_t)z * (size_t)chunk +
                                           (size_t)y) * (size_t)chunk;
                    for (int x = 0; x < chunk; x++)
                        dst[x] = (uint8_t)(src[x] ? 255 : 0);
                }
            }
            char path[1024];
            snprintf(path, sizeof(path), "%s/%s.tif", out_pred,
                     entries[i].id);
            if (TiffIO_save(path, cube, chunk, chunk, chunk) != 0) {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    fprintf(stderr, "pred_fixup: FAILED to write %s\n", path);
                    write_err = 1;
                }
            }
            Arena_restore(warenas[tid], mark);
        }
        for (int t = 0; t < n_threads; t++) Arena_dispose(&warenas[t]);
        if (write_err) {
            fprintf(stderr, "pred_fixup: output grid INCOMPLETE\n");
        } else {
            cubes_written = n_entries;
        }
        char srcp[1024], dstp[1024];
        snprintf(srcp, sizeof(srcp), "%s/manifest.json", opts->grid_dir);
        snprintf(dstp, sizeof(dstp), "%s/manifest.json", opts->out_dir);
        pf_copy_file(srcp, dstp);
        snprintf(srcp, sizeof(srcp), "%s/cubes_PRED/present.json",
                 opts->grid_dir);
        snprintf(dstp, sizeof(dstp), "%s/cubes_PRED/present.json",
                 opts->out_dir);
        pf_copy_file(srcp, dstp);
    }
    double t_write = pf_now();
    fprintf(stderr,
            "pred_fixup: split summary — carved %lld px, %zu cubes written "
            "(load %.1fs, carve %.1fs, write %.1fs)\n",
            (long long)st.carved_px, cubes_written,
            t_load - t0, t_carve - t_load, t_write - t_carve);

    Arena_dispose(&arena);
    free(entries);
    return 0;
}

static int pf_run_auto(const PfOpts *opts)
{
    char pred_dir[1024];
    snprintf(pred_dir, sizeof(pred_dir), "%s/cubes_PRED", opts->grid_dir);
    PfEntry *entries = NULL;
    size_t n_entries = 0;
    if (pf_collect_cubes(pred_dir, &entries, &n_entries) != 0 ||
        n_entries == 0) {
        fprintf(stderr, "pred_fixup: no cubes in %s\n", pred_dir);
        free(entries);
        return 1;
    }
    int chunk = 0;
    {
        Arena_T probe = Arena_new();
        uint8_t *cv = NULL;
        int cd = 0, ch = 0, cw = 0;
        if (TiffIO_load(probe, entries[0].path, &cv, &cd, &ch, &cw) != 0 ||
            cd != ch || ch != cw || cd <= 0) {
            fprintf(stderr, "pred_fixup: bad first cube %s\n",
                    entries[0].path);
            Arena_dispose(&probe);
            free(entries);
            return 1;
        }
        chunk = cd;
        Arena_dispose(&probe);
    }
    int z0 = entries[0].vz, z1 = z0, y0 = entries[0].vy, y1 = y0;
    int x0 = entries[0].vx, x1 = x0;
    for (size_t i = 0; i < n_entries; i++) {
        if (entries[i].vz < z0) z0 = entries[i].vz;
        if (entries[i].vz > z1) z1 = entries[i].vz;
        if (entries[i].vy < y0) y0 = entries[i].vy;
        if (entries[i].vy > y1) y1 = entries[i].vy;
        if (entries[i].vx < x0) x0 = entries[i].vx;
        if (entries[i].vx > x1) x1 = entries[i].vx;
    }
    free(entries);
    int D = z1 - z0 + chunk;
    int ncy = (y1 - y0) / chunk + 1, ncx = (x1 - x0) / chunk + 1;
    unsigned long long whole =
        (unsigned long long)D * (unsigned long long)(ncy * chunk) *
        (unsigned long long)(ncx * chunk) * 6ull;
    const unsigned long long TILE_BUDGET = 16ull << 30;

    int T = opts->tile;
    if (T <= 0) {
        if (whole <= TILE_BUDGET) T = 0;
        else {
            for (T = (ncy > ncx ? ncy : ncx) - 1; T > 1; T--) {
                unsigned long long need =
                    (unsigned long long)D *
                    (unsigned long long)(T * chunk) *
                    (unsigned long long)(T * chunk) * 6ull;
                if (need <= TILE_BUDGET) break;
            }
            if (T < 1) T = 1;
        }
    }
    if (T >= ncy && T >= ncx) T = 0;
    if (T == 0)
        return pf_run(opts, NULL, NULL);

    fprintf(stderr,
            "pred_fixup: TILED — grid %dx%d cubes, tile side %d cubes, "
            "full z (%d planes); tile seams are no-join zones\n",
            ncy, ncx, T, D);
    PfTotals tot;
    memset(&tot, 0, sizeof(tot));
    int rc = 0, tile_idx = 0;
    for (int cy = 0; cy < ncy && rc == 0; cy += T) {
        for (int cx = 0; cx < ncx && rc == 0; cx += T) {
            PfRegion reg;
            memset(&reg, 0, sizeof(reg));
            reg.cy0 = cy;
            reg.cy1 = (cy + T < ncy) ? cy + T : ncy;
            reg.cx0 = cx;
            reg.cx1 = (cx + T < ncx) ? cx + T : ncx;
            reg.append = (tile_idx > 0);
            snprintf(reg.tag, sizeof(reg.tag), "t%02d%02d_", cy / T, cx / T);
            fprintf(stderr,
                    "\npred_fixup: ==== TILE %s y[%d,%d) x[%d,%d) cubes "
                    "====\n",
                    reg.tag, reg.cy0, reg.cy1, reg.cx0, reg.cx1);
            rc = pf_run(opts, &reg, &tot);
            tile_idx++;
        }
    }
    if (rc == 0) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/fixup_report.json", opts->out_dir);
        FILE *f = fopen(path, "w");
        if (f) {
            fprintf(f,
                "{\n  \"grid\": \"%s\",\n  \"tiled\": {\"side_cubes\": %d, "
                "\"tiles\": %d},\n  \"endpoints\": %lld,\n"
                "  \"candidates\": %lld,\n  \"joins\": %lld,\n"
                "  \"kept\": %lld,\n  \"far_kept\": %lld,\n"
                "  \"painted_px\": %lld,\n"
                "  \"bridges\": {\"necks\": %lld, \"certified\": %lld, "
                "\"persistent\": %lld, \"cut\": %lld, \"cut_px\": %lld},\n"
                "  \"cubes_written\": %zu\n}\n",
                opts->grid_dir, T, tile_idx,
                (long long)tot.endpoints, (long long)tot.candidates,
                (long long)tot.joins, (long long)tot.kept,
                (long long)tot.far_kept, (long long)tot.painted_px,
                (long long)tot.bridges, (long long)tot.certified,
                (long long)tot.persistent, (long long)tot.cut,
                (long long)tot.cut_px, tot.cubes_written);
            fclose(f);
        }
        fprintf(stderr,
                "\n=== pred_fixup GRAND TOTAL (%d tiles) ===\n"
                "  joins kept    %lld / %lld (%lld far-tier)\n"
                "  painted px    %lld\n"
                "  weld necks    %lld found, %lld certified, %lld persistent, "
                "%lld cut (%lld px)\n"
                "  cubes         %zu written\n",
                tile_idx, (long long)tot.kept, (long long)tot.joins,
                (long long)tot.far_kept, (long long)tot.painted_px,
                (long long)tot.bridges, (long long)tot.certified,
                (long long)tot.persistent, (long long)tot.cut,
                (long long)tot.cut_px, tot.cubes_written);
    }
    return rc;
}

/* ---------------------------------------------------------------- */
/* Selftest                                                         */
/* ---------------------------------------------------------------- */

/* End-to-end mini grid: two 16^3 cubes side by side in x, a 2px-thick line
 * crossing the cube boundary with a PERSISTENT ~7px far-tier gap AT the
 * boundary in every plane. Proves (a) full-plane processing sees across
 * cubes, and (b) the persistent-gap path works: round 1 provisionally
 * admits the far-tier pair (no corridor evidence exists — the gap persists
 * in every plane), tracks form, round 2 accepts on track support, the
 * support filter keeps it, painting bridges the seam. */
static int pf_e2e_selftest(void)
{
    const char *tmp = getenv("TEMP");
    if (!tmp) tmp = getenv("TMP");
    if (!tmp) tmp = "/tmp";
    char base[1024], grid[1024], pred[1024], out[1024];
    snprintf(base, sizeof(base), "%s/pred_fixup_selftest", tmp);
    snprintf(grid, sizeof(grid), "%s/grid", base);
    snprintf(pred, sizeof(pred), "%s/cubes_PRED", grid);
    snprintf(out, sizeof(out), "%s/out", base);
    ves_mkdir(base);
    ves_mkdir(grid);
    ves_mkdir(pred);

    enum { C = 16 };
    Arena_T arena = Arena_new();
    uint8_t *cube = ARENA_CALLOC(arena, (size_t)C * C * C, 1);
    size_t fg_in = 0;

    /* cube 0 (x 0..15): line rows y=8..9, x=1..12 (skeleton must stay >=
     * FIXUP_MIN_CURVE_PX after thinning retracts the tips) */
    for (int z = 0; z < C; z++)
        for (int y = 8; y <= 9; y++)
            for (int x = 1; x <= 12; x++) {
                cube[((size_t)z * C + (size_t)y) * C + (size_t)x] = 255;
                fg_in++;
            }
    char path[1024];
    snprintf(path, sizeof(path), "%s/z00000_y00000_x00000.tif", pred);
    if (TiffIO_save(path, cube, C, C, C) != 0) {
        fprintf(stderr, "[selftest] e2e: cannot write %s\n", path);
        Arena_dispose(&arena);
        return 3;
    }

    /* cube 1 (x 16..31): line rows y=8..9, x(world) = 19..30 -> local 3..14 */
    memset(cube, 0, (size_t)C * C * C);
    for (int z = 0; z < C; z++)
        for (int y = 8; y <= 9; y++)
            for (int x = 3; x <= 14; x++) {
                cube[((size_t)z * C + (size_t)y) * C + (size_t)x] = 255;
                fg_in++;
            }
    snprintf(path, sizeof(path), "%s/z00000_y00000_x00016.tif", pred);
    if (TiffIO_save(path, cube, C, C, C) != 0) {
        fprintf(stderr, "[selftest] e2e: cannot write %s\n", path);
        Arena_dispose(&arena);
        return 3;
    }

    PfOpts opts;
    pf_default_opts(&opts);
    opts.grid_dir = grid;
    opts.out_dir = out;
    opts.quiet_viz = 1;
    opts.png_max = 0;
    int rc = pf_run(&opts, NULL, NULL);
    if (rc != 0) {
        fprintf(stderr, "[selftest] e2e: pf_run rc=%d -> FAIL\n", rc);
        Arena_dispose(&arena);
        return 3;
    }

    /* reload the two output cubes; the gap must now be bridged */
    size_t fg_out = 0;
    int bridged = 1;
    const char *ids[2] = {"z00000_y00000_x00000", "z00000_y00000_x00016"};
    for (int i = 0; i < 2; i++) {
        snprintf(path, sizeof(path), "%s/cubes_PRED/%s.tif", out, ids[i]);
        uint8_t *cv = NULL;
        int cd = 0, chh = 0, cww = 0;
        if (TiffIO_load(arena, path, &cv, &cd, &chh, &cww) != 0 ||
            cd != C || chh != C || cww != C) {
            fprintf(stderr, "[selftest] e2e: cannot reload %s\n", path);
            Arena_dispose(&arena);
            return 3;
        }
        for (int k = 0; k < C * C * C; k++)
            if (cv[k]) fg_out++;
        /* gap world x = 14..17: cube0 local x 14..15, cube1 local x 0..1 */
        int xa = (i == 0) ? 14 : 0, xb = (i == 0) ? 15 : 1;
        for (int z = 2; z < C - 2; z++)
            for (int x = xa; x <= xb; x++)
                if (!cv[((size_t)z * C + 8) * C + (size_t)x]) bridged = 0;
    }
    int ok = (fg_out > fg_in && bridged);
    fprintf(stderr,
            "[selftest] pred_fixup e2e mini-grid (fg %zu -> %zu, bridged=%d) "
            "-> %s\n",
            fg_in, fg_out, bridged, ok ? "ok" : "FAIL");

    /* Tiled negative: --tile 1 puts the seam exactly at the gap; the
     * border-exclusion band must refuse the cross-tile join entirely. */
    char out2[1024];
    snprintf(out2, sizeof(out2), "%s/out_tiled", base);
    PfOpts topts;
    pf_default_opts(&topts);
    topts.grid_dir = grid;
    topts.out_dir = out2;
    topts.quiet_viz = 1;
    topts.png_max = 0;
    topts.tile = 1;
    int rc2 = pf_run_auto(&topts);
    size_t fg_tiled = 0;
    int tiled_ok = (rc2 == 0);
    for (int i = 0; i < 2 && tiled_ok; i++) {
        snprintf(path, sizeof(path), "%s/cubes_PRED/%s.tif", out2, ids[i]);
        uint8_t *cv = NULL;
        int cd = 0, chh = 0, cww = 0;
        if (TiffIO_load(arena, path, &cv, &cd, &chh, &cww) != 0) {
            tiled_ok = 0;
            break;
        }
        for (int k = 0; k < C * C * C; k++)
            if (cv[k]) fg_tiled++;
    }
    if (tiled_ok) tiled_ok = (fg_tiled == fg_in);
    fprintf(stderr,
            "[selftest] pred_fixup e2e tiled seam-refusal (fg %zu, "
            "unchanged=%d) -> %s\n",
            fg_tiled, fg_tiled == fg_in, tiled_ok ? "ok" : "FAIL");

    Arena_dispose(&arena);
    return (ok && tiled_ok) ? 0 : 3;
}

static int pf_selftest(void)
{
    int rc = 0;
    rc |= SliceTrace_selftest();
    rc |= JoinPaint_selftest();
    rc |= SliceMatch_selftest();
    rc |= JoinTracks_selftest();
    rc |= BridgeScan_selftest();
    rc |= FixupViz_selftest();
    rc |= SlabSplit_selftest();
    rc |= pf_e2e_selftest();
    fprintf(stderr, "=== pred_fixup selftest %s ===\n",
            rc ? "FAILED" : "passed");
    return rc ? 3 : 0;
}

/* ---------------------------------------------------------------- */
/* main                                                             */
/* ---------------------------------------------------------------- */

static void pf_usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s <grid_dir> <out_dir> [options]\n"
        "       %s --selftest\n"
        "Options: --umb-y F --umb-x F --reach-max F --reach-safe F\n"
        "         --radial-dr F --min-support N --min-score F\n"
        "         --paint-radius N --no-tracks --cut-bridges --dry-run --png-all\n"
        "         --png-scale N --png-max N --planes z0:z1 --threads N --tile N\n"
        "         --split-slabs [--slab-thick-r F --slab-ext-r F\n"
        "         --slab-rounds N]   (own pass: carve fused slabs, no joins)\n",
        prog, prog);
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return pf_selftest();
    if (argc < 3) {
        pf_usage(argv[0]);
        return 2;
    }

    PfOpts opts;
    pf_default_opts(&opts);
    opts.grid_dir = argv[1];
    opts.out_dir = argv[2];

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--umb-y") && i + 1 < argc)
            opts.umb_y = atof(argv[++i]);
        else if (!strcmp(argv[i], "--umb-x") && i + 1 < argc)
            opts.umb_x = atof(argv[++i]);
        else if (!strcmp(argv[i], "--reach-max") && i + 1 < argc)
            opts.reach_max = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--reach-safe") && i + 1 < argc)
            opts.reach_safe = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--radial-dr") && i + 1 < argc)
            opts.radial_dr = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--min-support") && i + 1 < argc)
            opts.min_support = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min-score") && i + 1 < argc)
            opts.min_score = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--paint-radius") && i + 1 < argc)
            opts.paint_radius = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--no-tracks"))
            opts.no_tracks = 1;
        else if (!strcmp(argv[i], "--cut-bridges"))
            opts.cut_bridges = 1;
        else if (!strcmp(argv[i], "--split-slabs"))
            opts.split_slabs = 1;
        else if (!strcmp(argv[i], "--slab-thick-r") && i + 1 < argc)
            opts.slab_thick_r = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--slab-ext-r") && i + 1 < argc)
            opts.slab_ext_r = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "--slab-rounds") && i + 1 < argc)
            opts.slab_rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--dry-run"))
            opts.dry_run = 1;
        else if (!strcmp(argv[i], "--png-all"))
            opts.png_all = 1;
        else if (!strcmp(argv[i], "--png-scale") && i + 1 < argc)
            opts.png_scale = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--png-max") && i + 1 < argc)
            opts.png_max = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--crop-max") && i + 1 < argc)
            opts.crop_max = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--rej-crop-max") && i + 1 < argc)
            opts.rej_crop_max = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--planes") && i + 1 < argc) {
            if (sscanf(argv[++i], "%d:%d", &opts.zlo, &opts.zhi) != 2) {
                fprintf(stderr, "bad --planes (want z0:z1)\n");
                return 2;
            }
        } else if (!strcmp(argv[i], "--threads") && i + 1 < argc)
            opts.threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--tile") && i + 1 < argc)
            opts.tile = atoi(argv[++i]);
        else {
            fprintf(stderr, "unknown arg: %s\n", argv[i]);
            pf_usage(argv[0]);
            return 2;
        }
    }

    if (opts.split_slabs) return pf_run_split(&opts);
    return pf_run_auto(&opts);
}
