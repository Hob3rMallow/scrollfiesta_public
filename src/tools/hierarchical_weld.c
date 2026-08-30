/*
 * hierarchical_weld.c -- mesh LOD-pyramid orchestrator.
 *
 * Builds a level-of-detail pyramid over a grid of already-meshed cubes. For
 * each level L = 1..N it partitions the level-(L-1) nodes into fanout^3 blocks,
 * WELDS each block (grid_weld.exe --subgrid), then SIMPLIFIES the welded block
 * a lot (qslim_obj.exe). Each decimated block becomes a node for level L+1.
 * Because every weld is bounded to <= fanout^3 inputs and each level is
 * decimated before the next weld, the pipeline never materializes a
 * full-resolution monolith -- sidestepping the whole-grid weld OOM -- while
 * producing a unified mesh at the top and a usable tile set at every tier.
 *
 * Usage:
 *   hierarchical_weld <leaf_dump_dir> <output_dir> [options]
 *   hierarchical_weld --selftest
 *
 * This is a PURE ORCHESTRATOR: directory enumeration + subprocess spawning +
 * text aggregation. It never holds a mesh in memory (no arena / obj_io / qem).
 *
 * ---- Global-anchoring rule (correctness-critical) ----
 * grid_weld's detect_planes finds seams only at GLOBAL multiples of --cube-size.
 * So every level-L node origin must be a global multiple of P_L = 128*fanout^L.
 * We anchor each block origin O = floor(child_origin / P_L) * P_L per axis, weld
 * with --subgrid over the present children and --cube-size = P_{L-1} (the child
 * pitch). This makes detect_planes find exactly the seams BETWEEN the grouped
 * children and skip the already-closed finer seams inside each child. A
 * min-relative / index grouping would put a seam at a non-multiple of the weld's
 * cube-size, so the two halves would concatenate UNWELDED with no error.
 *
 * ---- Merger safety ----
 * grid_weld caps bridge faces at 2*rho_max=6 vox < 7-vox inter-wrap clearance,
 * and QEM edge-collapse only simplifies existing topology (proximity anti-fusion
 * guard), so decimation cannot merge two wraps. The optional winding gate is
 * armed for every weld via SEAM_UMBILICUS_Y/X + SEAM_WRAP_PITCH env (set once in
 * main, inherited by all children).
 */
#include "../common/ves_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <omp.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
  #include <sys/wait.h>
  #include <unistd.h>
  #include <fcntl.h>
#endif

/* ================================================================
 * Types
 * ================================================================ */

typedef struct { long long z, y, x; } Origin;

typedef struct { Origin *o; size_t n, cap; } OriginList;

typedef struct {
    Origin      origin;     /* anchored block origin (global multiple of P_L) */
    Origin      lo, hi;     /* inclusive bbox over PRESENT children (--subgrid) */
    int         n_children;
} Block;

/* Per-block outcome, filled by run_one_block (one slot per block -> no races). */
#define RES_OK   0
#define RES_WARN 1   /* weld reported non-manifold/pinch but OBJ is usable */
#define RES_SKIP 2   /* final node already complete (resume) */
#define RES_FAIL 3   /* spawn fail / timeout / crash / qslim fail -> no output */

typedef struct {
    int           status;
    unsigned long weld_exit, qslim_exit;
    long long     weld_verts, weld_faces;             /* from weld_report.json */
    long long     unpaired, non_manifold, same_dir, pinch;
    long long     final_faces;                        /* counted post-qslim */
    int           n_children;
    int           no_holefill;                        /* 1 = weld retried w/o holefill */
    int           repair_passes;
    long long     shortcuts_before, shortcuts_after;
    long long     repair_cut_edges;
    long long     repair_cut_rank;
    long long     repair_components_delta;
    long long     repair_boundary_delta;
    int           topology_validated;
    double        seconds;
} BlockResult;

typedef struct {
    const char *leaf_dir;
    const char *out_dir;
    const char *leaf_stage;
    int         fanout;
    int         max_levels;
    int         stop_nodes;
    int         max_concurrent;
    double      keep_ratio;
    double      block_timeout_sec;
    int         skip_existing;
    int         simplify_top;
    int         no_decimate;   /* 1 = never qslim: weld-only every level (band-CVT
                                * welds already carry CVT quality; decimation only
                                * re-scars them with slivers). Trades bounded
                                * per-level memory for preserved quality. */
    int         cvt_simplify;  /* 1 = decimate each level with cvt_simplify (CVT
                                * per-component, boundary-preserving) instead of
                                * qslim: a proper mip pyramid that KEEPS CVT
                                * quality at every tier. Seams coarsen with the
                                * rest and are re-refined by the next weld. */
    int         dry_run;
    int         remesh;        /* 1 = pass --remesh to qslim (default on) */
    int         topology_first;/* preserve faces, disable generative closers,
                                * iteratively cut supported winding shortcuts,
                                * then exact/topology/winding validate each node */
    int         repair_passes;
    double      shortcut_tol;
    double      umb_y, umb_x, wrap_pitch;
    const char *axis_table;
    const char *grid_weld_exe;
    const char *qslim_exe;
    const char *cvt_simplify_exe;
    const char *manifold_exe;
    const char *obj_reorient_exe;
    const char *wind_audit_exe;
    const char *seam_audit_exe;
    const char *remesh_exe;
} Options;

/* Level-invariant context passed to per-block work. */
typedef struct {
    const Options *opt;
    int         level;
    const char *prev_dir;
    const char *child_stage;
    const char *this_dir;
    const char *this_stage;
    long long   child_pitch;   /* P_{level-1} = grid_weld --cube-size */
    Origin      domain_lo;     /* original leaf-domain lower bound */
    Origin      domain_hi;     /* original leaf-domain upper bound, exclusive */
    int         do_simplify;   /* 0 = weld only (terminal), 1 = weld + qslim */
} LevelCtx;

/* ================================================================
 * Pure helpers (unit-tested via --selftest)
 * ================================================================ */

/* P_L = 128 * fanout^L (integer). */
static long long level_pitch(int fanout, int level)
{
    long long p = 128;
    int i = 0;
    for (i = 0; i < level; i++) p *= (long long)fanout;
    return p;
}

/* Anchor an origin to the level-L pitch grid: floor(v/P)*P (origins are >= 0). */
static Origin block_origin_of(Origin c, int fanout, int level)
{
    long long P = level_pitch(fanout, level);
    Origin o;
    o.z = (c.z / P) * P;
    o.y = (c.y / P) * P;
    o.x = (c.x / P) * P;
    return o;
}

static int origin_eq(Origin a, Origin b)
{
    return a.z == b.z && a.y == b.y && a.x == b.x;
}

static void block_id_str(char *buf, size_t n, Origin o)
{
    snprintf(buf, n, "z%05lld_y%05lld_x%05lld",
             (long long)o.z, (long long)o.y, (long long)o.x);
}

static int parse_origin(const char *id, Origin *o)
{
    long long z = 0, y = 0, x = 0;
    if (sscanf(id, "z%lld_y%lld_x%lld", &z, &y, &x) != 3) return -1;
    o->z = z; o->y = y; o->x = x;
    return 0;
}

/* Flat per-level layout: <dir>/<id>/<id>_<stage>_all.obj -- ONE directory per
 * node, holding _weld.obj, the logs, and the final all.obj side by side (no
 * nested <id>_<stage>/ subdir). grid_weld's node resolver stats the old nested
 * <id>/<id>_<stage>/<id>_<stage>_all.obj form first, then falls back to this
 * flat form, so leaf dumps (nested) and LOD levels (flat) both read back. */
static void node_obj_path(char *buf, size_t n, const char *dir,
                          const char *id, const char *stage)
{
    snprintf(buf, n, "%s/%s/%s_%s_all.obj",
             dir, id, id, stage);
}

/* Group prev-level origins into anchored blocks; fill present-children bbox +
 * count. Returns malloc'd array (caller frees); deterministic (sorted). */
static int enumerate_level_blocks(const OriginList *prev, int fanout, int level,
                                  Block **out, size_t *out_n)
{
    Block *blks = NULL;
    size_t n = 0, cap = 0, i = 0, j = 0;

    for (i = 0; i < prev->n; i++) {
        Origin c = prev->o[i];
        Origin bo = block_origin_of(c, fanout, level);
        for (j = 0; j < n; j++)
            if (origin_eq(blks[j].origin, bo)) break;
        if (j == n) {
            if (n >= cap) {
                size_t ncap = cap ? cap * 2 : 64;
                Block *nb = (Block *)realloc(blks, ncap * sizeof(Block));
                if (!nb) { free(blks); return -1; }
                blks = nb; cap = ncap;
            }
            memset(&blks[n], 0, sizeof(Block));
            blks[n].origin = bo;
            blks[n].lo = c;
            blks[n].hi = c;
            blks[n].n_children = 0;
            j = n; n++;
        }
        if (c.z < blks[j].lo.z) blks[j].lo.z = c.z;
        if (c.y < blks[j].lo.y) blks[j].lo.y = c.y;
        if (c.x < blks[j].lo.x) blks[j].lo.x = c.x;
        if (c.z > blks[j].hi.z) blks[j].hi.z = c.z;
        if (c.y > blks[j].hi.y) blks[j].hi.y = c.y;
        if (c.x > blks[j].hi.x) blks[j].hi.x = c.x;
        blks[j].n_children++;
    }

    /* Sort by (z,y,x) for a stable summary and deterministic node ordering. */
    for (i = 0; i + 1 < n; i++) {
        for (j = i + 1; j < n; j++) {
            Origin a = blks[i].origin, b = blks[j].origin;
            int gt = (a.z > b.z) ||
                     (a.z == b.z && a.y > b.y) ||
                     (a.z == b.z && a.y == b.y && a.x > b.x);
            if (gt) { Block t = blks[i]; blks[i] = blks[j]; blks[j] = t; }
        }
    }

    *out = blks;
    *out_n = n;
    return 0;
}

/* A globally anchored hierarchy can plateau for one or more logical levels
 * when the occupied domain straddles the next coarser grid (PHerc1447 does so
 * at P=4096).  Re-running grid_weld on unary blocks is an exact no-op: there is
 * no sibling seam to weld, and the sole child already carries a topology
 * certificate.  Treat such a level as virtual and advance to the next
 * grouping pitch while retaining the last materialized node pitch. */
static int blocks_all_unary(const Block *blocks, size_t n)
{
    size_t i = 0;
    if (!blocks || n == 0) return 0;
    for (i = 0; i < n; i++)
        if (blocks[i].n_children != 1) return 0;
    return 1;
}

/* Clip the geometric support of the PRESENT children to the leaf domain.
 * The anchored parent origin is only a grouping key.  In particular, after a
 * virtual unary plateau it can be several child pitches away from the occupied
 * support and must not define the shortcut audit's boundary evidence. */
static void block_boundary_box(const Block *block, long long child_pitch,
                               Origin domain_lo, Origin domain_hi,
                               long long boundary[6])
{
    long long z1 = block->hi.z > LLONG_MAX - child_pitch ?
                   LLONG_MAX : block->hi.z + child_pitch;
    long long y1 = block->hi.y > LLONG_MAX - child_pitch ?
                   LLONG_MAX : block->hi.y + child_pitch;
    long long x1 = block->hi.x > LLONG_MAX - child_pitch ?
                   LLONG_MAX : block->hi.x + child_pitch;
    boundary[0] = block->lo.z > domain_lo.z ? block->lo.z : domain_lo.z;
    boundary[1] = z1 < domain_hi.z ? z1 : domain_hi.z;
    boundary[2] = block->lo.y > domain_lo.y ? block->lo.y : domain_lo.y;
    boundary[3] = y1 < domain_hi.y ? y1 : domain_hi.y;
    boundary[4] = block->lo.x > domain_lo.x ? block->lo.x : domain_lo.x;
    boundary[5] = x1 < domain_hi.x ? x1 : domain_hi.x;
}

/* ================================================================
 * Small I/O helpers
 * ================================================================ */

static void origins_init(OriginList *l) { l->o = NULL; l->n = 0; l->cap = 0; }
static void origins_free(OriginList *l) { free(l->o); l->o = NULL; l->n = l->cap = 0; }

static int origins_push(OriginList *l, Origin o)
{
    if (l->n >= l->cap) {
        size_t ncap = l->cap ? l->cap * 2 : 256;
        Origin *no = (Origin *)realloc(l->o, ncap * sizeof(Origin));
        if (!no) return -1;
        l->o = no; l->cap = ncap;
    }
    l->o[l->n++] = o;
    return 0;
}

static int ensure_dir(const char *path)
{
    char fake_child[1024];
    snprintf(fake_child, sizeof(fake_child), "%s/x", path);
    return ves_ensure_parent_dir(fake_child);
}

static int file_exists(const char *p)
{
#ifdef _WIN32
    return _access(p, 0) == 0;
#else
    return access(p, 0) == 0;
#endif
}

/* Resume an expensive derived artifact only when it is complete and at least
 * as new as its source.  A bare existence check can otherwise reuse a stale
 * topology precondition after a weld has been regenerated. */
static int file_not_older_than(const char *candidate, const char *source)
{
#ifdef _WIN32
    WIN32_FILE_ATTRIBUTE_DATA c, s;
    if (!GetFileAttributesExA(candidate, GetFileExInfoStandard, &c) ||
        !GetFileAttributesExA(source, GetFileExInfoStandard, &s))
        return 0;
    return CompareFileTime(&c.ftLastWriteTime, &s.ftLastWriteTime) >= 0;
#else
    struct stat c, s;
    if (stat(candidate, &c) != 0 || stat(source, &s) != 0) return 0;
    return c.st_mtime >= s.st_mtime;
#endif
}

/* A node OBJ counts as complete iff it exists, is non-trivial, and ends in a
 * clean newline (a process killed mid-dump leaves no trailing '\n'). Mirrors
 * grid_pipeline.c:obj_looks_complete. */
static int obj_looks_complete(const char *path)
{
    FILE *f = fopen(path, "rb");
    long sz = 0;
    int last = 0;
    if (!f) return 0;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return 0; }
    sz = ftell(f);
    if (sz < 128) { fclose(f); return 0; }
    if (fseek(f, -1, SEEK_END) != 0) { fclose(f); return 0; }
    last = fgetc(f);
    fclose(f);
    return last == '\n';
}

static int atomic_replace(const char *tmp, const char *dst)
{
#ifdef _WIN32
    return MoveFileExA(tmp, dst, MOVEFILE_REPLACE_EXISTING) ? 0 : -1;
#else
    return rename(tmp, dst);
#endif
}

static void set_env(const char *k, const char *v)
{
#ifdef _WIN32
    _putenv_s(k, v);
#else
    setenv(k, v, 1);
#endif
}

/* Count triangle faces in an OBJ (lines beginning "f "). Cheap summary stat. */
static long long count_obj_faces(const char *path)
{
    FILE *f = fopen(path, "rb");
    long long nf = 0;
    char line[512];
    if (!f) return 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) nf++;
    }
    fclose(f);
    return nf;
}

/* Extract an integer value for a quoted JSON key, e.g. "\"non_manifold\"". */
static long long json_get_ll(const char *buf, const char *quoted_key)
{
    const char *p = strstr(buf, quoted_key);
    if (!p) return 0;
    p += strlen(quoted_key);
    while (*p && *p != ':') p++;
    if (*p == ':') p++;
    while (*p == ' ' || *p == '\t') p++;
    return strtoll(p, NULL, 10);
}

/* Extract a floating-point JSON scalar. Missing keys return NaN so semantic
 * gates fail closed rather than silently treating an older report as zero. */
static double json_get_double(const char *buf, const char *quoted_key)
{
    const char *p = strstr(buf, quoted_key);
    if (!p) return NAN;
    p += strlen(quoted_key);
    while (*p && *p != ':') p++;
    if (*p == ':') p++;
    while (*p == ' ' || *p == '\t') p++;
    return strtod(p, NULL);
}

static void json_write_string(FILE *f, const char *s)
{
    const unsigned char *p = (const unsigned char *)(s ? s : "");
    fputc('"', f);
    while (*p) {
        switch (*p) {
        case '"': fputs("\\\"", f); break;
        case '\\': fputs("\\\\", f); break;
        case '\b': fputs("\\b", f); break;
        case '\f': fputs("\\f", f); break;
        case '\n': fputs("\\n", f); break;
        case '\r': fputs("\\r", f); break;
        case '\t': fputs("\\t", f); break;
        default:
            if (*p < 0x20) fprintf(f, "\\u%04x", (unsigned)*p);
            else fputc((int)*p, f);
            break;
        }
        p++;
    }
    fputc('"', f);
}

/* Evaluate the measured curved axis and its local tangent directly from the
 * small z,y,x CSV. The table is intentionally re-read once per hierarchy node:
 * this keeps worker state immutable, and the file has only O(10) rows. */
static int axis_table_eval_file(const char *path, double z,
                                double *y, double *x,
                                double *dy_dz, double *dx_dz)
{
    double zz[1024], yy[1024], xx[1024];
    size_t n = 0, i = 0;
    char line[1024];
    FILE *f = fopen(path, "r");
    if (f == NULL) return -1;
    while (n < 1024 && fgets(line, sizeof(line), f) != NULL) {
        double tz, ty, tx;
        if (sscanf(line, " %lf , %lf , %lf", &tz, &ty, &tx) != 3)
            continue;
        if (!isfinite(tz) || !isfinite(ty) || !isfinite(tx)) continue;
        if (n > 0 && tz <= zz[n-1]) { fclose(f); return -1; }
        zz[n] = tz; yy[n] = ty; xx[n] = tx; n++;
    }
    fclose(f);
    if (n < 2) return -1;
    if (z <= zz[0]) i = 0;
    else if (z >= zz[n-1]) i = n - 2;
    else {
        for (i = 0; i + 1 < n; i++)
            if (z <= zz[i+1]) break;
    }
    {
        double dz = zz[i+1] - zz[i];
        double t = (z - zz[i]) / dz;
        *dy_dz = (yy[i+1] - yy[i]) / dz;
        *dx_dz = (xx[i+1] - xx[i]) / dz;
        *y = yy[i] + t * (yy[i+1] - yy[i]);
        *x = xx[i] + t * (xx[i+1] - xx[i]);
    }
    return 0;
}

static int read_file_all(const char *path, char *buf, size_t bufsz)
{
    FILE *f = fopen(path, "rb");
    size_t n = 0;
    if (!f) return -1;
    n = fread(buf, 1, bufsz - 1, f);
    fclose(f);
    buf[n] = '\0';
    return 0;
}

typedef struct {
    long long faces, before_tracks, post_tracks, before_hits, post_hits;
    long long cuts, cut_rank, components_before, components_after;
    long long small_before, small_after, boundary_edges;
    double before_excess, post_excess;
} RepairPassMetrics;

/* Read only the bounded scalar prefix of wind_audit's JSON.  Large per-track
 * arrays follow these fields and are deliberately irrelevant to checkpoint
 * validation. */
static int read_repair_pass_metrics(const char *path, RepairPassMetrics *m)
{
    char buf[8192];
    static const char *required[] = {
        "\"faces\"", "\"shortcut_persistent_clusters\"",
        "\"shortcut_persistent_hits\"", "\"shortcut_persistent_excess\"",
        "\"repair_cut_edges\"", "\"repair_cut_rank\"",
        "\"repair_components_before\"", "\"repair_components_after\"",
        "\"repair_small_components_before\"",
        "\"repair_small_components_after\"",
        "\"repair_post_shortcut_clusters\"",
        "\"repair_post_shortcut_hits\"",
        "\"repair_post_shortcut_excess\"", "\"boundary_edges\""
    };
    if (m == NULL || read_file_all(path, buf, sizeof(buf)) != 0) return -1;
    for (size_t i = 0; i < sizeof(required) / sizeof(required[0]); i++)
        if (strstr(buf, required[i]) == NULL) return -1;
    memset(m, 0, sizeof(*m));
    m->faces = json_get_ll(buf, "\"faces\"");
    m->before_tracks = json_get_ll(buf, "\"shortcut_persistent_clusters\"");
    m->post_tracks = json_get_ll(buf, "\"repair_post_shortcut_clusters\"");
    m->before_hits = json_get_ll(buf, "\"shortcut_persistent_hits\"");
    m->post_hits = json_get_ll(buf, "\"repair_post_shortcut_hits\"");
    m->before_excess = json_get_double(buf, "\"shortcut_persistent_excess\"");
    m->post_excess = json_get_double(buf, "\"repair_post_shortcut_excess\"");
    m->cuts = json_get_ll(buf, "\"repair_cut_edges\"");
    m->cut_rank = json_get_ll(buf, "\"repair_cut_rank\"");
    m->components_before = json_get_ll(buf, "\"repair_components_before\"");
    m->components_after = json_get_ll(buf, "\"repair_components_after\"");
    m->small_before = json_get_ll(buf, "\"repair_small_components_before\"");
    m->small_after = json_get_ll(buf, "\"repair_small_components_after\"");
    m->boundary_edges = json_get_ll(buf, "\"boundary_edges\"");
    return 0;
}

/* OBJ round-tripping uses nine significant digits before float reload.  Hits
 * exactly on a tracking/join threshold can therefore merge or split one tiny
 * temporal track between an in-memory postcheck and the next pass's serialized
 * precheck.  Permit only that measured-scale discrepancy when validating a
 * checkpoint chain; timestamps, face counts, energy monotonicity, and all cut
 * budgets remain exact gates. */
static int serialized_track_counts_agree(long long a, long long b)
{
    long long hi, delta, tolerance;
    if (a < 0 || b < 0) return 0;
    hi = a > b ? a : b;
    delta = a > b ? a - b : b - a;
    tolerance = 2 + hi / 200; /* two tracks plus 0.5% */
    return delta <= tolerance;
}

/* ---- Enumerate node subdirs (names starting 'z', parseable origin) ---- */
static int enumerate_node_dirs(const char *dir, OriginList *out)
{
#ifdef _WIN32
    char glob[1024];
    WIN32_FIND_DATAA fd;
    HANDLE h;
    snprintf(glob, sizeof(glob), "%s/*", dir);
    h = FindFirstFileA(glob, &fd);
    if (h == INVALID_HANDLE_VALUE) return -1;
    do {
        Origin o;
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) continue;
        if (fd.cFileName[0] != 'z') continue;
        if (parse_origin(fd.cFileName, &o) != 0) continue;
        if (origins_push(out, o) != 0) { FindClose(h); return -1; }
    } while (FindNextFileA(h, &fd));
    FindClose(h);
    return 0;
#else
    DIR *d = opendir(dir);
    struct dirent *de = NULL;
    if (!d) return -1;
    while ((de = readdir(d)) != NULL) {
        Origin o;
        if (de->d_name[0] != 'z') continue;
        if (parse_origin(de->d_name, &o) != 0) continue;
        if (origins_push(out, o) != 0) { closedir(d); return -1; }
    }
    closedir(d);
    return 0;
#endif
}

/* ================================================================
 * Subprocess with per-child log capture + wall-clock timeout.
 * Returns 0 = ran to completion (*out_exit set), -1 = spawn failure,
 * -2 = timed out (child was killed).
 * ================================================================ */
static int spawn_logged(const char *exe, const char *const *argv,
                        const char *log_path, double timeout_sec,
                        unsigned long *out_exit)
{
#ifdef _WIN32
    char cmdline[8192];
    size_t pos = 0;
    int i = 0;
    SECURITY_ATTRIBUTES sa;
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    HANDLE hLog;
    BOOL ok;
    DWORD wait_ms, wr, code = 1;

    for (i = 0; argv[i] != NULL; i++) {
        size_t alen = strlen(argv[i]);
        if (i > 0 && pos < sizeof(cmdline) - 1) cmdline[pos++] = ' ';
        if (pos < sizeof(cmdline) - 1) cmdline[pos++] = '"';
        if (pos + alen < sizeof(cmdline) - 2) {
            memcpy(cmdline + pos, argv[i], alen);
            pos += alen;
        }
        if (pos < sizeof(cmdline) - 1) cmdline[pos++] = '"';
    }
    cmdline[pos] = '\0';

    sa.nLength = sizeof(sa);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;
    hLog = CreateFileA(log_path, GENERIC_WRITE,
                       FILE_SHARE_READ | FILE_SHARE_WRITE, &sa,
                       CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hLog == INVALID_HANDLE_VALUE) return -1;

    memset(&si, 0, sizeof(si));
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdOutput = hLog;
    si.hStdError  = hLog;
    si.hStdInput  = GetStdHandle(STD_INPUT_HANDLE);
    memset(&pi, 0, sizeof(pi));

    ok = CreateProcessA(exe, cmdline, NULL, NULL, TRUE, 0, NULL, NULL, &si, &pi);
    CloseHandle(hLog);
    if (!ok) return -1;

    wait_ms = (timeout_sec > 0) ? (DWORD)(timeout_sec * 1000.0) : INFINITE;
    wr = WaitForSingleObject(pi.hProcess, wait_ms);
    if (wr == WAIT_TIMEOUT) {
        TerminateProcess(pi.hProcess, 1);
        WaitForSingleObject(pi.hProcess, INFINITE);
        CloseHandle(pi.hProcess);
        CloseHandle(pi.hThread);
        return -2;
    }
    GetExitCodeProcess(pi.hProcess, &code);
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
    *out_exit = code;
    return 0;
#else
    int fd = open(log_path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    pid_t pid;
    double t0;
    if (fd < 0) return -1;
    pid = fork();
    if (pid < 0) { close(fd); return -1; }
    if (pid == 0) {
        dup2(fd, STDOUT_FILENO);
        dup2(fd, STDERR_FILENO);
        close(fd);
        execv(exe, (char *const *)argv);
        _exit(127);
    }
    close(fd);
    t0 = ves_clock_sec();
    for (;;) {
        int status = 0;
        pid_t r = waitpid(pid, &status, WNOHANG);
        if (r == pid) {
            if (WIFEXITED(status)) *out_exit = (unsigned long)WEXITSTATUS(status);
            else *out_exit = 1;
            return 0;
        }
        if (r < 0) return -1;
        if (timeout_sec > 0 && ves_clock_sec() - t0 > timeout_sec) {
            kill(pid, SIGKILL);
            waitpid(pid, &status, 0);
            return -2;
        }
        ves_sleep_ms(50);
    }
#endif
}

/* Normalize the output of a weld before semantic cutting. grid_weld already
 * performs both operations once, but exact conflict surgery can expose a few
 * same-direction knots, and knot cutting can expose a final geometric stab.
 * The convergent order orient -> exact repair -> orient removes those residuals
 * while each tool retains its own fail-closed budget. */
static int topology_precondition(const LevelCtx *lc, const Block *blk,
                                 const char *id, const char *input,
                                 char *output, size_t output_cap)
{
    const Options *opt = lc->opt;
    char orient1[1088], exact[1088], log[1088];
    char z_s[64], y_s[64], x_s[64], dy_s[64], dx_s[64];
    const char *argv[24];
    unsigned long ex = 1;
    int ac, st;
    double z, y, x, dy, dx;

    z = 0.5 * ((double)blk->lo.z + (double)blk->hi.z +
               (double)lc->child_pitch);
    if (axis_table_eval_file(opt->axis_table, z, &y, &x, &dy, &dx) != 0)
        return -1;
    snprintf(z_s, sizeof(z_s), "%.12g", z);
    snprintf(y_s, sizeof(y_s), "%.12g", y);
    snprintf(x_s, sizeof(x_s), "%.12g", x);
    snprintf(dy_s, sizeof(dy_s), "%.12g", dy);
    snprintf(dx_s, sizeof(dx_s), "%.12g", dx);
    snprintf(orient1, sizeof(orient1), "%s/%s/_pre_orient1.obj",
             lc->this_dir, id);
    snprintf(exact, sizeof(exact), "%s/%s/_pre_exact.obj",
             lc->this_dir, id);
    snprintf(output, output_cap, "%s/%s/_preconditioned.obj",
             lc->this_dir, id);
    (void)remove(orient1); (void)remove(exact); (void)remove(output);

    snprintf(log, sizeof(log), "%s/%s/_pre_orient1.log",
             lc->this_dir, id);
    ac = 0;
    argv[ac++] = opt->obj_reorient_exe;
    argv[ac++] = input; argv[ac++] = orient1;
    argv[ac++] = "--axis-point";
    argv[ac++] = z_s; argv[ac++] = y_s; argv[ac++] = x_s;
    argv[ac++] = "--axis-dir";
    argv[ac++] = "1"; argv[ac++] = dy_s; argv[ac++] = dx_s;
    argv[ac] = NULL;
    st = spawn_logged(opt->obj_reorient_exe, argv, log,
                      opt->block_timeout_sec, &ex);
    if (st != 0 || ex != 0 || !obj_looks_complete(orient1)) goto fail;

    snprintf(log, sizeof(log), "%s/%s/_pre_exact.log",
             lc->this_dir, id);
    ac = 0;
    argv[ac++] = opt->seam_audit_exe;
    argv[ac++] = orient1;
    argv[ac++] = "--repair"; argv[ac++] = exact;
    argv[ac++] = "--hinge";
    argv[ac++] = "--repair-max-fraction"; argv[ac++] = "0.01";
    argv[ac] = NULL;
    st = spawn_logged(opt->seam_audit_exe, argv, log,
                      opt->block_timeout_sec, &ex);
    if (st != 0 || ex != 0 || !obj_looks_complete(exact)) goto fail;

    snprintf(log, sizeof(log), "%s/%s/_pre_orient2.log",
             lc->this_dir, id);
    ac = 0;
    argv[ac++] = opt->obj_reorient_exe;
    argv[ac++] = exact; argv[ac++] = output;
    argv[ac++] = "--axis-point";
    argv[ac++] = z_s; argv[ac++] = y_s; argv[ac++] = x_s;
    argv[ac++] = "--axis-dir";
    argv[ac++] = "1"; argv[ac++] = dy_s; argv[ac++] = dx_s;
    argv[ac] = NULL;
    st = spawn_logged(opt->obj_reorient_exe, argv, log,
                      opt->block_timeout_sec, &ex);
    if (st != 0 || ex != 0 || !obj_looks_complete(output)) goto fail;

    (void)remove(orient1); (void)remove(exact);
    return 0;
fail:
    (void)remove(orient1); (void)remove(exact); (void)remove(output);
    return -1;
}

/* Iterative semantic repair for the topology-first hierarchy. Exact triangle
 * cleanup cannot see a smooth staircase connecting neighbouring plies, so each
 * welded node is sectioned in the measured winding coordinate. A supported
 * phase-isoline is cut and the mesh is re-audited. Integrated shortcut excess
 * may remain on a quantized plateau but must never increase; a finite pass cap
 * and cumulative cut budget make plateau traversal bounded. A regular-value
 * search avoids saddle/branch cuts and rising-energy midpoint cuts.
 * Mutual-nearest temporal tracking keeps nearby parallel bridges separate so
 * each receives its own seam. Faces and coordinates are never changed; only
 * vertex incidence along the cut is split.
 *
 * Acceptance is deliberately fail-closed:
 *   - <=1.5% cumulative cut edges;
 *   - component and small-fragment growth bounded by the embedded cut graph's
 *     conservative separation rank (chains + branches + endpoint pairs);
 *   - boundary growth exactly 2 per cut interior edge;
 *   - zero curved-axis shortcuts/full-turn mergers on an independent final run;
 *   - zero non-manifold/same-direction defects and zero exact conflicts.
 * A validation marker is written LAST and is required for resume. */
static int topology_first_repair(const LevelCtx *lc, const Block *blk,
                                 const char *id,
                                 const char *weld_path, const char *final_path,
                                 BlockResult *res)
{
    const Options *opt = lc->opt;
    char current[1088], next[1088], temp_a[1088], temp_b[1088];
    char json_path[1088], log_path[1088], final_json[1088];
    char manifold_log[1088], seam_log[1088], marker[1088];
    char pitch_s[32], tol_s[32], boundary_s[6][32];
    const char *argv[40];
    long long initial_faces = 0, initial_tracks = 0, initial_components = 0;
    long long initial_small = 0, initial_boundary = 0;
    long long last_components = 0, last_small = 0;
    long long total_cuts = 0, total_cut_rank = 0, final_boundary = 0;
    int pass = 0, clean = 0;
    long long boundary[6];

    block_boundary_box(blk, lc->child_pitch, lc->domain_lo, lc->domain_hi,
                       boundary);
    for (int k = 0; k < 6; k++)
        snprintf(boundary_s[k], sizeof(boundary_s[k]), "%lld", boundary[k]);

    snprintf(current, sizeof(current), "%s", weld_path);
    snprintf(temp_a, sizeof(temp_a), "%s/%s/_shortcut_a.obj",
             lc->this_dir, id);
    snprintf(temp_b, sizeof(temp_b), "%s/%s/_shortcut_b.obj",
             lc->this_dir, id);
    snprintf(pitch_s, sizeof(pitch_s), "%.9g", opt->wrap_pitch);
    snprintf(tol_s, sizeof(tol_s), "%.9g", opt->shortcut_tol);

    /* Resume from the newest pass whose candidate and complete scalar-report
     * chain can be tied to this exact preconditioned source.  The two candidate
     * paths alternate; requiring JSON >= candidate catches a newer partial
     * overwrite paired with an older report. */
    if (opt->skip_existing) {
        int checkpoint = 0;
        for (int p = opt->repair_passes; p >= 1; p--) {
            char candidate[1088], report[1088];
            RepairPassMetrics m;
            snprintf(candidate, sizeof(candidate), "%s", (p & 1) ? temp_a : temp_b);
            snprintf(report, sizeof(report), "%s/%s/_shortcut_pass%02d.json",
                     lc->this_dir, id, p);
            if (obj_looks_complete(candidate) &&
                file_not_older_than(candidate, weld_path) &&
                file_not_older_than(report, candidate) &&
                read_repair_pass_metrics(report, &m) == 0) {
                checkpoint = p;
                break;
            }
        }
        if (checkpoint > 0) {
            int valid = 1;
            long long previous_post = -1;
            RepairPassMetrics last;
            memset(&last, 0, sizeof(last));
            for (int p = 1; p <= checkpoint; p++) {
                char report[1088];
                RepairPassMetrics m;
                snprintf(report, sizeof(report),
                         "%s/%s/_shortcut_pass%02d.json",
                         lc->this_dir, id, p);
                if (!file_not_older_than(report, weld_path) ||
                    read_repair_pass_metrics(report, &m) != 0) {
                    valid = 0;
                    break;
                }
                if (p == 1) {
                    initial_faces = m.faces;
                    initial_tracks = m.before_tracks;
                    initial_components = m.components_before;
                    initial_small = m.small_before;
                    initial_boundary = m.boundary_edges;
                }
                if (m.faces != initial_faces || m.faces <= 0 ||
                    m.cuts < 0 || m.cut_rank < 0 ||
                    m.components_after <= 0 || m.small_after < 0 ||
                    m.before_tracks < 0 || m.post_tracks < 0 ||
                    m.before_hits < 0 || m.post_hits < 0 ||
                    !isfinite(m.before_excess) || !isfinite(m.post_excess) ||
                    (p > 1 && !serialized_track_counts_agree(
                        m.before_tracks, previous_post))) {
                    valid = 0;
                    break;
                }
                total_cuts += m.cuts;
                total_cut_rank += m.cut_rank;
                last_components = m.components_after;
                last_small = m.small_after;
                if (total_cuts > (initial_faces * 15 + 999) / 1000 ||
                    last_components > initial_components + total_cut_rank ||
                    last_small > initial_small + total_cut_rank ||
                    (m.post_tracks > 0 &&
                     (m.cuts <= 0 || m.post_excess > m.before_excess +
                      fmax(1e-4, 1e-4 * m.before_excess)))) {
                    valid = 0;
                    break;
                }
                previous_post = m.post_tracks;
                last = m;
            }
            if (valid) {
                snprintf(current, sizeof(current), "%s",
                         (checkpoint & 1) ? temp_a : temp_b);
                pass = checkpoint;
                clean = last.post_tracks == 0;
                fprintf(stderr,
                        "    semantic checkpoint: pass %d, tracks %lld->%lld, cuts %lld\n",
                        checkpoint, initial_tracks, last.post_tracks,
                        total_cuts);
            } else {
                total_cuts = total_cut_rank = 0;
                initial_faces = initial_tracks = initial_components = 0;
                initial_small = initial_boundary = 0;
                last_components = last_small = 0;
            }
        }
    }

    if (!clean) for (pass = pass + 1; pass <= opt->repair_passes; pass++) {
        unsigned long ex = 1;
        int ac = 0, st;
        char buf[8192];
        long long before_tracks, post_tracks, before_hits, post_hits;
        long long cuts, cut_rank, components_after, small_after, faces;
        double before_excess, post_excess;
        snprintf(next, sizeof(next), "%s", (pass & 1) ? temp_a : temp_b);
        (void)remove(next);
        snprintf(json_path, sizeof(json_path), "%s/%s/_shortcut_pass%02d.json",
                 lc->this_dir, id, pass);
        snprintf(log_path, sizeof(log_path), "%s/%s/_shortcut_pass%02d.log",
                 lc->this_dir, id, pass);

        argv[ac++] = opt->wind_audit_exe;
        argv[ac++] = current;
        argv[ac++] = "--axis-table"; argv[ac++] = opt->axis_table;
        argv[ac++] = "--pitch"; argv[ac++] = pitch_s;
        argv[ac++] = "--shortcut-tol"; argv[ac++] = tol_s;
        argv[ac++] = "--shortcut-length-pitches"; argv[ac++] = "2";
        argv[ac++] = "--shortcut-z-step"; argv[ac++] = "4";
        argv[ac++] = "--shortcut-min-planes"; argv[ac++] = "3";
        argv[ac++] = "--shortcut-track-mutual";
        argv[ac++] = "--no-split-audit";
        argv[ac++] = "--boundary-box";
        for (int k = 0; k < 6; k++) argv[ac++] = boundary_s[k];
        argv[ac++] = "--repair-auto-offset";
        argv[ac++] = "--repair-max-cut-fraction"; argv[ac++] = "0.015";
        argv[ac++] = "--repair-seams"; argv[ac++] = next;
        argv[ac++] = "--repair-emit-rejected";
        argv[ac++] = "--json"; argv[ac++] = json_path;
        argv[ac] = NULL;
        st = spawn_logged(opt->wind_audit_exe, argv, log_path,
                          opt->block_timeout_sec, &ex);
        if (st != 0 || (ex != 0 && ex != 5) ||
            !obj_looks_complete(next) ||
            read_file_all(json_path, buf, sizeof(buf)) != 0)
            goto fail;

        before_tracks = json_get_ll(buf, "\"shortcut_persistent_clusters\"");
        post_tracks = json_get_ll(buf, "\"repair_post_shortcut_clusters\"");
        before_hits = json_get_ll(buf, "\"shortcut_persistent_hits\"");
        post_hits = json_get_ll(buf, "\"repair_post_shortcut_hits\"");
        before_excess = json_get_double(buf, "\"shortcut_persistent_excess\"");
        post_excess = json_get_double(buf, "\"repair_post_shortcut_excess\"");
        cuts = json_get_ll(buf, "\"repair_cut_edges\"");
        cut_rank = json_get_ll(buf, "\"repair_cut_rank\"");
        faces = json_get_ll(buf, "\"faces\"");
        components_after = json_get_ll(buf, "\"repair_components_after\"");
        small_after = json_get_ll(buf, "\"repair_small_components_after\"");
        if (pass == 1) {
            initial_faces = faces;
            initial_tracks = before_tracks;
            initial_components = json_get_ll(buf, "\"repair_components_before\"");
            initial_small = json_get_ll(buf, "\"repair_small_components_before\"");
            initial_boundary = json_get_ll(buf, "\"boundary_edges\"");
        }
        if (faces != initial_faces || cuts < 0 || cut_rank < 0 ||
            components_after <= 0 ||
            small_after < 0 || post_tracks < 0 || before_tracks < 0 ||
            before_hits < 0 || post_hits < 0 || !isfinite(before_excess) ||
            !isfinite(post_excess))
            goto fail;
        total_cuts += cuts;
        total_cut_rank += cut_rank;
        last_components = components_after;
        last_small = small_after;
        if (total_cuts > (initial_faces * 15 + 999) / 1000 ||
            last_components > initial_components + total_cut_rank ||
            last_small > initial_small + total_cut_rank)
            goto fail;
        if (post_tracks > 0 &&
            (cuts <= 0 || post_excess > before_excess +
             fmax(1e-4, 1e-4 * before_excess)))
            goto fail;  /* no cut, or semantic energy rose */

        if (strcmp(current, weld_path) != 0) (void)remove(current);
        snprintf(current, sizeof(current), "%s", next);
        if (post_tracks == 0) { clean = 1; break; }
    }
    if (!clean) goto fail;

    /* Independent curved-axis audit of the emitted candidate (not the repair
     * routine's in-memory postcheck), including the full-turn edge gate. */
    {
        unsigned long ex = 1;
        int ac = 0, st;
        char buf[8192];
        snprintf(final_json, sizeof(final_json), "%s/%s/_shortcut_final.json",
                 lc->this_dir, id);
        snprintf(log_path, sizeof(log_path), "%s/%s/_shortcut_final.log",
                 lc->this_dir, id);
        argv[ac++] = opt->wind_audit_exe;
        argv[ac++] = current;
        argv[ac++] = "--axis-table"; argv[ac++] = opt->axis_table;
        argv[ac++] = "--pitch"; argv[ac++] = pitch_s;
        argv[ac++] = "--shortcut-tol"; argv[ac++] = tol_s;
        argv[ac++] = "--shortcut-length-pitches"; argv[ac++] = "2";
        argv[ac++] = "--shortcut-z-step"; argv[ac++] = "4";
        argv[ac++] = "--shortcut-min-planes"; argv[ac++] = "3";
        argv[ac++] = "--shortcut-track-mutual";
        argv[ac++] = "--boundary-box";
        for (int k = 0; k < 6; k++) argv[ac++] = boundary_s[k];
        argv[ac++] = "--fail-on-local";
        argv[ac++] = "--fail-on-shortcut";
        argv[ac++] = "--json"; argv[ac++] = final_json;
        argv[ac] = NULL;
        st = spawn_logged(opt->wind_audit_exe, argv, log_path,
                          opt->block_timeout_sec, &ex);
        if (st != 0 || ex != 0 ||
            read_file_all(final_json, buf, sizeof(buf)) != 0 ||
            json_get_ll(buf, "\"shortcut_persistent_clusters\"") != 0 ||
            json_get_ll(buf, "\"merge_ft\"") != 0 ||
            json_get_ll(buf, "\"faces\"") != initial_faces)
            goto fail;
        final_boundary = json_get_ll(buf, "\"boundary_edges\"");
        if (final_boundary - initial_boundary != 2 * total_cuts)
            goto fail;
    }

    /* Combinatorial manifold + orientation gate. manifold_check intentionally
     * exits zero for findings, so parse its brief record as well as exit code. */
    {
        unsigned long ex = 1;
        int ac = 0, st;
        char buf[4096];
        snprintf(manifold_log, sizeof(manifold_log),
                 "%s/%s/_shortcut_manifold.log", lc->this_dir, id);
        argv[ac++] = opt->manifold_exe;
        argv[ac++] = "--brief";
        argv[ac++] = current;
        argv[ac] = NULL;
        st = spawn_logged(opt->manifold_exe, argv, manifold_log,
                          opt->block_timeout_sec, &ex);
        if (st != 0 || ex != 0 ||
            read_file_all(manifold_log, buf, sizeof(buf)) != 0 ||
            strstr(buf, "nm_edge=0") == NULL ||
            strstr(buf, "nm_vert=0") == NULL ||
            strstr(buf, "same_dir=0") == NULL)
            goto fail;
    }

    /* Exact global intersection/fold-back gate. */
    {
        unsigned long ex = 1;
        int ac = 0, st;
        snprintf(seam_log, sizeof(seam_log),
                 "%s/%s/_shortcut_intersections.log", lc->this_dir, id);
        argv[ac++] = opt->seam_audit_exe;
        argv[ac++] = current;
        argv[ac++] = "--check";
        argv[ac++] = "--hinge";
        argv[ac] = NULL;
        st = spawn_logged(opt->seam_audit_exe, argv, seam_log,
                          opt->block_timeout_sec, &ex);
        if (st != 0 || ex != 0) goto fail;
    }

    if (atomic_replace(current, final_path) != 0) goto fail;
    (void)remove(temp_a); (void)remove(temp_b);
    snprintf(marker, sizeof(marker), "%s.topology.json", final_path);
    {
        FILE *f = fopen(marker, "w");
        if (f == NULL) goto fail;
        fprintf(f,
            "{\n"
            "  \"validated\": true,\n"
            "  \"axis_table\": ");
        json_write_string(f, opt->axis_table);
        fprintf(f,
            ",\n"
            "  \"pitch\": %.9g, \"shortcut_tolerance\": %.9g,\n"
            "  \"boundary_box_zyx\": [%lld, %lld, %lld, %lld, %lld, %lld],\n"
            "  \"repair_passes\": %d,\n"
            "  \"shortcuts_before\": %lld, \"shortcuts_after\": 0,\n"
            "  \"cut_edges\": %lld,\n"
            "  \"cut_graph_rank\": %lld,\n"
            "  \"faces\": %lld,\n"
            "  \"components_before\": %lld, \"components_after\": %lld,\n"
            "  \"small_components_before\": %lld, \"small_components_after\": %lld,\n"
            "  \"boundary_edges_before\": %lld, \"boundary_edges_after\": %lld,\n"
            "  \"manifold\": true, \"same_dir\": 0, \"exact_conflicts\": 0\n"
            "}\n",
            opt->wrap_pitch, opt->shortcut_tol,
            boundary[0], boundary[1], boundary[2], boundary[3],
            boundary[4], boundary[5], pass,
            initial_tracks, total_cuts, total_cut_rank, initial_faces,
            initial_components, last_components, initial_small, last_small,
            initial_boundary, final_boundary);
        if (fclose(f) != 0) goto fail;
    }
    res->repair_passes = pass;
    res->shortcuts_before = initial_tracks;
    res->shortcuts_after = 0;
    res->repair_cut_edges = total_cuts;
    res->repair_cut_rank = total_cut_rank;
    res->repair_components_delta = last_components - initial_components;
    res->repair_boundary_delta = final_boundary - initial_boundary;
    res->topology_validated = 1;
    return 0;

fail:
    if (strcmp(current, weld_path) != 0) (void)remove(current);
    (void)remove(temp_a); (void)remove(temp_b);
    return -1;
}

/* ================================================================
 * Per-block driver: weld -> (simplify) -> atomic finalize.
 * ================================================================ */
static void run_one_block(const LevelCtx *lc, const Block *blk, BlockResult *res)
{
    const Options *opt = lc->opt;
    char id[128], final[1024], weld_out[1024], report[1088], marker[1088];
    char preconditioned[1088];
    char weld_log[1024], qslim_log[1024], final_tmp[1088];
    char sz0[32], sz1[32], sy0[32], sy1[32], sx0[32], sx1[32];
    char pitch_s[32], ratio_s[32];
    const char *argv[40];
    int ac = 0;
    int ran_weld = 0;
    unsigned long ex = 1;
    int st = 0;
    double t0 = ves_clock_sec();

    memset(res, 0, sizeof(*res));
    res->n_children = blk->n_children;

    block_id_str(id, sizeof(id), blk->origin);
    node_obj_path(final, sizeof(final), lc->this_dir, id, lc->this_stage);
    snprintf(marker, sizeof(marker), "%s.topology.json", final);

    /* Resume: final node already complete. */
    if (opt->skip_existing && obj_looks_complete(final) &&
        (!opt->topology_first || file_exists(marker))) {
        res->status = RES_SKIP;
        res->topology_validated = opt->topology_first ? 1 : 0;
        res->final_faces = count_obj_faces(final);
        res->seconds = ves_clock_sec() - t0;
        return;
    }

    if (ves_ensure_parent_dir(final) != 0) { res->status = RES_FAIL; return; }

    /* Weld destination: for weld-only (terminal) levels, grid_weld writes the
     * final node directly; otherwise it writes _weld.obj which qslim reads. */
    if (lc->do_simplify || opt->topology_first)
        snprintf(weld_out, sizeof(weld_out), "%s/%s/_weld.obj", lc->this_dir, id);
    else
        snprintf(weld_out, sizeof(weld_out), "%s", final);
    if (ves_ensure_parent_dir(weld_out) != 0) { res->status = RES_FAIL; return; }
    snprintf(report, sizeof(report), "%s.weld_report.json", weld_out);

    snprintf(sz0, sizeof(sz0), "%lld", (long long)blk->lo.z);
    snprintf(sz1, sizeof(sz1), "%lld", (long long)blk->hi.z);
    snprintf(sy0, sizeof(sy0), "%lld", (long long)blk->lo.y);
    snprintf(sy1, sizeof(sy1), "%lld", (long long)blk->hi.y);
    snprintf(sx0, sizeof(sx0), "%lld", (long long)blk->lo.x);
    snprintf(sx1, sizeof(sx1), "%lld", (long long)blk->hi.x);
    snprintf(pitch_s, sizeof(pitch_s), "%lld", (long long)lc->child_pitch);

    /* Weld resume: grid_weld writes the report LAST, so report present => the
     * weld OBJ is complete. Skip re-welding in that case. */
    if (!(opt->skip_existing && file_exists(report) && obj_looks_complete(weld_out))) {
        ac = 0;
        argv[ac++] = opt->grid_weld_exe;
        argv[ac++] = lc->prev_dir;
        argv[ac++] = weld_out;
        argv[ac++] = "--subgrid";
        argv[ac++] = sz0; argv[ac++] = sz1;
        argv[ac++] = sy0; argv[ac++] = sy1;
        argv[ac++] = sx0; argv[ac++] = sx1;
        argv[ac++] = "--stage";     argv[ac++] = lc->child_stage;
        argv[ac++] = "--cube-size"; argv[ac++] = pitch_s;
        if (opt->axis_table != NULL) {
            argv[ac++] = "--axis-table";
            argv[ac++] = opt->axis_table;
        }
        /* Undecimated: the node-count cap (2M faces/node) is calibrated for
         * leaf cubes and far too small once a handful of nodes hold the whole
         * scroll. Give grid_weld a generous budget (under its 150M/80M
         * ceilings) so the terminal weld does not overflow. Costs upfront
         * allocation per weld -- pair with a lower --max-concurrent. */
        if (opt->no_decimate) {
            argv[ac++] = "--vert-cap"; argv[ac++] = "20000000";
            argv[ac++] = "--face-cap"; argv[ac++] = "30000000";
        }
        argv[ac] = NULL;
        ran_weld = 1;

        snprintf(weld_log, sizeof(weld_log), "%s/%s/_weld.log", lc->this_dir, id);
        st = spawn_logged(opt->grid_weld_exe, argv, weld_log,
                          opt->block_timeout_sec, &ex);
        res->weld_exit = ex;
        /* spawn failure (couldn't even launch) is unrecoverable. A TIMEOUT is
         * NOT fatal here -- grid_weld's interior hole-fill can hang the weld on
         * huge meshes -- so it falls through to the --no-holefill retry below. */
        if (st == -1) { res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return; }
    } else {
        res->weld_exit = 0;
    }

    /* Retry once with --no-holefill if the weld TIMED OUT (st==-2) or ran but
     * produced no complete output: grid_weld's interior hole-fill can hang /
     * crash / OOM on the huge, hole-riddled meshes at upper LOD levels (thousands
     * of loops, each with a 10 s safe_triangulate watchdog). Skipping it yields a
     * valid (if holey) coarse block instead of a gap in the pyramid. (A crash
     * AFTER a valid write leaves output present, so we never get here for that
     * case -- it's accepted as-is.) argv still holds the weld args, so append the
     * flag in the NULL slot. */
    if (ran_weld && (st == -2 || !(file_exists(report) && obj_looks_complete(weld_out)))) {
        argv[ac++] = "--no-holefill";
        argv[ac] = NULL;
        res->no_holefill = 1;
        snprintf(weld_log, sizeof(weld_log), "%s/%s/_weld_nohf.log", lc->this_dir, id);
        st = spawn_logged(opt->grid_weld_exe, argv, weld_log,
                          opt->block_timeout_sec, &ex);
        res->weld_exit = ex;
        if (st == -2 || st == -1) { res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return; }
    }

    /* Accept the weld iff it produced a COMPLETE OBJ + report, ignoring the exit
     * code. This (a) tolerates grid_weld crashing in arena teardown AFTER writing
     * valid output (seen at upper LOD levels on huge meshes -- exit > 0xFF but the
     * mesh is fine), and (b) rejects a weld that wrote nothing even after the
     * no-holefill retry. grid_weld writes the report LAST, so report-present =>
     * OBJ complete. */
    if (!(file_exists(report) && obj_looks_complete(weld_out))) {
        res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return;
    }

    /* Pull weld report stats. */
    {
        char rbuf[4096];
        if (read_file_all(report, rbuf, sizeof(rbuf)) == 0) {
            res->weld_verts   = json_get_ll(rbuf, "\"total_unique_verts\"");
            res->weld_faces   = json_get_ll(rbuf, "\"total_unique_faces\"");
            res->unpaired     = json_get_ll(rbuf, "\"unpaired\"");
            res->non_manifold = json_get_ll(rbuf, "\"non_manifold\"");
            res->same_dir     = json_get_ll(rbuf, "\"same_dir_pairs\"");
            res->pinch        = json_get_ll(rbuf, "\"pinch_verts\"");
        }
    }

    if (opt->topology_first) {
        snprintf(preconditioned, sizeof(preconditioned),
                 "%s/%s/_preconditioned.obj", lc->this_dir, id);
        if (!(opt->skip_existing &&
              obj_looks_complete(preconditioned) &&
              file_not_older_than(preconditioned, weld_out)) &&
            topology_precondition(lc, blk, id, weld_out,
                                  preconditioned,
                                  sizeof(preconditioned)) != 0) {
            res->status = RES_FAIL;
            res->seconds = ves_clock_sec() - t0;
            return;
        }
        if (topology_first_repair(lc, blk, id, preconditioned,
                                  final, res) != 0) {
            res->status = RES_FAIL;
            res->seconds = ves_clock_sec() - t0;
            return;
        }
        (void)remove(preconditioned);
        res->final_faces = count_obj_faces(final);
        res->status = RES_OK;
        res->seconds = ves_clock_sec() - t0;
        return;
    }

    if (lc->do_simplify && opt->cvt_simplify) {
        /* CVT-decimate (per-component, boundary-preserving) instead of qslim:
         * keeps CVT quality at every tier. Seams coarsen with the interior and
         * are re-refined by the next level's weld. */
        snprintf(final_tmp, sizeof(final_tmp), "%s.tmp", final);
        snprintf(ratio_s, sizeof(ratio_s), "%.6f", opt->keep_ratio);
        ac = 0;
        argv[ac++] = opt->cvt_simplify_exe;
        argv[ac++] = weld_out;
        argv[ac++] = final_tmp;
        argv[ac++] = "--keep-ratio";
        argv[ac++] = ratio_s;
        argv[ac] = NULL;
        snprintf(qslim_log, sizeof(qslim_log), "%s/%s/_cvtsimplify.log", lc->this_dir, id);
        st = spawn_logged(opt->cvt_simplify_exe, argv, qslim_log,
                          opt->block_timeout_sec, &ex);
        res->qslim_exit = ex;
        if (st != 0 || ex != 0) { res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return; }
        if (atomic_replace(final_tmp, final) != 0) { res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return; }
    } else if (lc->do_simplify) {
        snprintf(final_tmp, sizeof(final_tmp), "%s.tmp", final);
        snprintf(ratio_s, sizeof(ratio_s), "%.6f", opt->keep_ratio);
        ac = 0;
        argv[ac++] = opt->qslim_exe;
        argv[ac++] = weld_out;
        argv[ac++] = final_tmp;
        argv[ac++] = ratio_s;
        if (opt->remesh) argv[ac++] = "--remesh";
        argv[ac] = NULL;
        snprintf(qslim_log, sizeof(qslim_log), "%s/%s/_qslim.log", lc->this_dir, id);
        st = spawn_logged(opt->qslim_exe, argv, qslim_log,
                          opt->block_timeout_sec, &ex);
        res->qslim_exit = ex;
        if (st != 0 || ex != 0) { res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return; }
        if (atomic_replace(final_tmp, final) != 0) { res->status = RES_FAIL; res->seconds = ves_clock_sec() - t0; return; }
    } else if (opt->remesh && !opt->cvt_simplify) {
        /* Terminal weld-only apex: it never saw a remesh, and the BPA seam weld
         * leaves it sliver-heavy (measured 65% -> 47% faces <15deg after remesh).
         * Run the isotropic remesh as a cleanup (fail-closed, ~same face count).
         * On any spawn failure leave the weld-only final in place -- still valid. */
        snprintf(final_tmp, sizeof(final_tmp), "%s.tmp", final);
        ac = 0;
        argv[ac++] = opt->remesh_exe;
        argv[ac++] = final;
        argv[ac++] = final_tmp;
        argv[ac] = NULL;
        snprintf(qslim_log, sizeof(qslim_log), "%s/%s/_remesh.log", lc->this_dir, id);
        st = spawn_logged(opt->remesh_exe, argv, qslim_log,
                          opt->block_timeout_sec, &ex);
        if (st == 0 && ex == 0) (void)atomic_replace(final_tmp, final);
    }

    res->final_faces = count_obj_faces(final);
    res->status = (res->non_manifold > 0 || res->pinch > 0) ? RES_WARN : RES_OK;
    res->seconds = ves_clock_sec() - t0;
}

/* ================================================================
 * Terminal re-orient: qslim's global winding-repair can flip regions on a
 * wrapped surface and (unlike intermediate levels) no further weld re-orients
 * the top node. Re-run grid_weld on the single node with a huge cube-size so
 * detect_planes finds 0 seams (no bridge) but the orient/cleanup tail still
 * runs. Only used with --simplify-top.
 * ================================================================ */
static void reorient_node(const Options *opt, const char *dir,
                          const char *stage, Origin o)
{
    char id[128], node[1024], tmp[1088], log[1024], sz0[32], sz1[32];
    char sy0[32], sy1[32], sx0[32], sx1[32];
    const char *argv[24];
    int ac = 0;
    unsigned long ex = 1;

    block_id_str(id, sizeof(id), o);
    node_obj_path(node, sizeof(node), dir, id, stage);
    snprintf(tmp, sizeof(tmp), "%s.orient.tmp", node);
    snprintf(log, sizeof(log), "%s/%s/_orient.log", dir, id);
    snprintf(sz0, sizeof(sz0), "%lld", (long long)o.z);
    snprintf(sz1, sizeof(sz1), "%lld", (long long)o.z);
    snprintf(sy0, sizeof(sy0), "%lld", (long long)o.y);
    snprintf(sy1, sizeof(sy1), "%lld", (long long)o.y);
    snprintf(sx0, sizeof(sx0), "%lld", (long long)o.x);
    snprintf(sx1, sizeof(sx1), "%lld", (long long)o.x);

    ac = 0;
    argv[ac++] = opt->grid_weld_exe;
    argv[ac++] = dir;
    argv[ac++] = tmp;
    argv[ac++] = "--subgrid";
    argv[ac++] = sz0; argv[ac++] = sz1;
    argv[ac++] = sy0; argv[ac++] = sy1;
    argv[ac++] = sx0; argv[ac++] = sx1;
    argv[ac++] = "--stage";     argv[ac++] = stage;
    argv[ac++] = "--cube-size"; argv[ac++] = "100000000";
    argv[ac] = NULL;

    if (spawn_logged(opt->grid_weld_exe, argv, log, opt->block_timeout_sec, &ex) == 0
        && ex <= 0xFF && obj_looks_complete(tmp)) {
        atomic_replace(tmp, node);
        fprintf(stderr, "  [reorient] %s\n", id);
    } else {
        fprintf(stderr, "  [reorient] FAILED for %s (exit=%lu) -- keeping qslim output\n",
                id, ex);
    }
}

/* ================================================================
 * Level driver: fan out blocks, aggregate, verify.
 * ================================================================ */
static int run_level(const LevelCtx *lc, Block *blocks, size_t nblk,
                     OriginList *out_nodes)
{
    const Options *opt = lc->opt;
    BlockResult *res = (BlockResult *)calloc(nblk ? nblk : 1, sizeof(BlockResult));
    char csv_path[1024], list_path[1024], mlog[1024];
    FILE *cf = NULL, *lf = NULL;
    int i = 0, nthreads = 0;
    int done = 0;
    size_t k = 0;
    long long tot_final = 0;
    int n_ok = 0, n_warn = 0, n_skip = 0, n_fail = 0, worst_nm = 0;

    if (!res) return -1;
    if (ensure_dir(lc->this_dir) != 0) { free(res); return -1; }

    nthreads = opt->max_concurrent;
    if (nthreads < 1) nthreads = 1;
    if (nthreads > (int)nblk) nthreads = (int)nblk;
    if (nthreads < 1) nthreads = 1;

    fprintf(stderr, "\n=== Level %d: %zu blocks (cube-size %lld, keep %.3f, %s) ===\n",
            lc->level, nblk, (long long)lc->child_pitch, opt->keep_ratio,
            lc->do_simplify ? "weld+simplify" : "weld-only");

    #pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
    for (i = 0; i < (int)nblk; i++) {
        run_one_block(lc, &blocks[i], &res[i]);
        #pragma omp critical
        {
            char id[128];
            const char *tag = "OK";
            done++;
            block_id_str(id, sizeof(id), blocks[i].origin);
            if (res[i].status == RES_WARN) tag = "WARN";
            else if (res[i].status == RES_SKIP) tag = "SKIP";
            else if (res[i].status == RES_FAIL) tag = "FAIL";
            fprintf(stderr, "  [L%d %d/%zu] %-5s %s  ch=%d wf=%lld ff=%lld nm=%lld "
                    "short=%lld->%lld cut=%lld%s (%.1fs)\n",
                    lc->level, done, nblk, tag, id, res[i].n_children,
                    (long long)res[i].weld_faces, (long long)res[i].final_faces,
                    (long long)res[i].non_manifold,
                    (long long)res[i].shortcuts_before,
                    (long long)res[i].shortcuts_after,
                    (long long)res[i].repair_cut_edges,
                    res[i].no_holefill ? " nohf" : "", res[i].seconds);
        }
    }

    /* Per-level CSV. */
    snprintf(csv_path, sizeof(csv_path), "%s/_summary.csv", lc->this_dir);
    cf = fopen(csv_path, "w");
    if (cf) {
        fprintf(cf, "block_id,n_children,status,weld_exit,qslim_exit,"
                    "weld_verts,weld_faces,final_faces,unpaired,non_manifold,"
                    "same_dir,pinch,repair_passes,shortcuts_before,"
                    "shortcuts_after,repair_cut_edges,repair_components_delta,"
                    "repair_boundary_delta,topology_validated,seconds\n");
        for (k = 0; k < nblk; k++) {
            char id[128];
            const char *tag = "ok";
            block_id_str(id, sizeof(id), blocks[k].origin);
            if (res[k].status == RES_WARN) tag = "warn";
            else if (res[k].status == RES_SKIP) tag = "skip";
            else if (res[k].status == RES_FAIL) tag = "fail";
            fprintf(cf, "%s,%d,%s,%lu,%lu,%lld,%lld,%lld,%lld,%lld,%lld,%lld,"
                        "%d,%lld,%lld,%lld,%lld,%lld,%d,%.2f\n",
                    id, res[k].n_children, tag, res[k].weld_exit, res[k].qslim_exit,
                    (long long)res[k].weld_verts, (long long)res[k].weld_faces,
                    (long long)res[k].final_faces, (long long)res[k].unpaired,
                    (long long)res[k].non_manifold, (long long)res[k].same_dir,
                    (long long)res[k].pinch, res[k].repair_passes,
                    (long long)res[k].shortcuts_before,
                    (long long)res[k].shortcuts_after,
                    (long long)res[k].repair_cut_edges,
                    (long long)res[k].repair_components_delta,
                    (long long)res[k].repair_boundary_delta,
                    res[k].topology_validated, res[k].seconds);
        }
        fclose(cf);
    }

    /* _objs.txt for manifold_check + feed usable nodes forward. */
    snprintf(list_path, sizeof(list_path), "%s/_objs.txt", lc->this_dir);
    lf = fopen(list_path, "w");
    for (k = 0; k < nblk; k++) {
        char id[128], node[1024];
        if (res[k].status == RES_OK) n_ok++;
        else if (res[k].status == RES_WARN) n_warn++;
        else if (res[k].status == RES_SKIP) n_skip++;
        else { n_fail++; continue; }
        if (res[k].non_manifold > worst_nm) worst_nm = (int)res[k].non_manifold;
        tot_final += res[k].final_faces;
        block_id_str(id, sizeof(id), blocks[k].origin);
        node_obj_path(node, sizeof(node), lc->this_dir, id, lc->this_stage);
        if (lf) fprintf(lf, "%s\n", node);
        origins_push(out_nodes, blocks[k].origin);
    }
    if (lf) fclose(lf);

    /* Aggregate manifold audit over the level's node OBJs. */
    if (out_nodes->n > 0) {
        const char *margv[6];
        unsigned long mex = 1;
        int mac = 0;
        margv[mac++] = opt->manifold_exe;
        margv[mac++] = "--list";
        margv[mac++] = list_path;
        margv[mac] = NULL;
        snprintf(mlog, sizeof(mlog), "%s/_manifold.log", lc->this_dir);
        (void)spawn_logged(opt->manifold_exe, margv, mlog, opt->block_timeout_sec, &mex);
        fprintf(stderr, "  manifold_check: exit=%lu (see %s)\n", mex, mlog);
    }

    fprintf(stderr, "  Level %d done: ok=%d warn=%d skip=%d fail=%d  nodes=%zu"
            "  total_faces=%lld  worst_non_manifold=%d\n",
            lc->level, n_ok, n_warn, n_skip, n_fail, out_nodes->n,
            (long long)tot_final, worst_nm);

    free(res);
    (void)k;
    return n_fail;
}

/* ================================================================
 * Preview / dry-run: simulate the pyramid partitioning without spawning.
 * ================================================================ */
static void preview_pyramid(const Options *opt, const OriginList *level0)
{
    OriginList cur;
    int L = 0;
    origins_init(&cur);
    { size_t i; for (i = 0; i < level0->n; i++) origins_push(&cur, level0->o[i]); }

    fprintf(stderr, "Pyramid preview (fanout %d, stop_nodes %d):\n", opt->fanout, opt->stop_nodes);
    fprintf(stderr, "  L0: %zu nodes (leaf, pitch 128)\n", cur.n);
    for (L = 1; cur.n > (size_t)opt->stop_nodes && L <= opt->max_levels; L++) {
        Block *blocks = NULL;
        size_t nblk = 0;
        long long child_pitch = level_pitch(opt->fanout, L - 1);
        OriginList next;
        if (enumerate_level_blocks(&cur, opt->fanout, L, &blocks, &nblk) != 0) break;
        fprintf(stderr, "  L%d: %zu nodes (pitch %lld, cube-size %lld)\n",
                L, nblk, (long long)level_pitch(opt->fanout, L), (long long)child_pitch);
        if (opt->dry_run) {
            size_t b;
            for (b = 0; b < nblk && b < 4; b++) {
                char id[128];
                block_id_str(id, sizeof(id), blocks[b].origin);
                fprintf(stderr,
                    "      %s: grid_weld <dir> <out> --subgrid %lld %lld %lld %lld %lld %lld "
                    "--stage %s --cube-size %lld  (children=%d)\n",
                    id, (long long)blocks[b].lo.z, (long long)blocks[b].hi.z,
                    (long long)blocks[b].lo.y, (long long)blocks[b].hi.y,
                    (long long)blocks[b].lo.x, (long long)blocks[b].hi.x,
                    (L == 1) ? opt->leaf_stage : "L(prev)", (long long)child_pitch,
                    blocks[b].n_children);
            }
            if (nblk > 4) fprintf(stderr, "      ... (%zu more)\n", nblk - 4);
        }
        origins_init(&next);
        { size_t b; for (b = 0; b < nblk; b++) origins_push(&next, blocks[b].origin); }
        free(blocks);
        /* Equal node count can be a legitimate one-level plateau when the
         * crop straddles a global coarse-bin boundary. Advancing the pitch may
         * group those nodes at the next level; only expansion is invalid. */
        if (next.n > cur.n) {
            fprintf(stderr, "  (invalid expansion -- stopping)\n");
            origins_free(&next); break;
        }
        origins_free(&cur);
        cur = next;
    }
    origins_free(&cur);
}

/* ================================================================
 * Selftest (pure functions; no mesh I/O).
 * ================================================================ */
static int build_full_grid(OriginList *out)
{
    /* Synthetic full 4x21x21 grid (no rejects): z 4352 step128 x4, y 2048 x21,
     * x 1536 x21 -> 1764 origins. */
    int iz, iy, ix;
    origins_init(out);
    for (iz = 0; iz < 4; iz++)
        for (iy = 0; iy < 21; iy++)
            for (ix = 0; ix < 21; ix++) {
                Origin o;
                o.z = 4352 + 128 * iz;
                o.y = 2048 + 128 * iy;
                o.x = 1536 + 128 * ix;
                if (origins_push(out, o) != 0) return -1;
            }
    return 0;
}

static int selftest(void)
{
    int fails = 0;
    char buf[128];
    Origin o;
    OriginList cur;
    int L = 0;
    /* Expected node counts, global-anchored, fanout 2. */
    int expect[] = {242, 72, 12, 6, 4, 1};
    int ei = 0;

    #define CHECK(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); fails++; } } while (0)

    CHECK(level_pitch(2, 0) == 128, "level_pitch(2,0)");
    CHECK(level_pitch(2, 1) == 256, "level_pitch(2,1)");
    CHECK(level_pitch(2, 5) == 4096, "level_pitch(2,5)");
    CHECK(level_pitch(3, 2) == 128 * 9, "level_pitch(3,2)");
    CHECK(serialized_track_counts_agree(291, 290),
          "serialized one-track threshold tolerance");
    CHECK(!serialized_track_counts_agree(291, 280),
          "serialized track tolerance remains fail-closed");

    /* Semantic-pass checkpoints depend on a deliberately small scalar prefix,
     * not the potentially large per-track tail of wind_audit JSON. */
    {
        const char *path = "hierarchical_weld_selftest_repair.json";
        FILE *f = fopen(path, "wb");
        RepairPassMetrics m;
        CHECK(f != NULL, "repair checkpoint fixture open");
        if (f != NULL) {
            fprintf(f,
                "{\"faces\":42,\"boundary_edges\":9,"
                "\"shortcut_persistent_clusters\":7,"
                "\"shortcut_persistent_hits\":11,"
                "\"shortcut_persistent_excess\":3.5,"
                "\"repair_cut_edges\":4,\"repair_cut_rank\":2,"
                "\"repair_components_before\":3,"
                "\"repair_components_after\":4,"
                "\"repair_small_components_before\":1,"
                "\"repair_small_components_after\":1,"
                "\"repair_post_shortcut_clusters\":2,"
                "\"repair_post_shortcut_hits\":3,"
                "\"repair_post_shortcut_excess\":1.25,"
                "\"shortcut_tracks\":[{\"ignored\":true}]}\n");
            fclose(f);
            CHECK(read_repair_pass_metrics(path, &m) == 0,
                  "repair checkpoint scalar parse");
            CHECK(m.faces == 42 && m.before_tracks == 7 &&
                  m.post_tracks == 2 && m.cuts == 4 && m.cut_rank == 2 &&
                  fabs(m.before_excess - 3.5) < 1e-12 &&
                  fabs(m.post_excess - 1.25) < 1e-12,
                  "repair checkpoint scalar values");
            remove(path);
        }
    }

    { Origin c; c.z = 4480; c.y = 2176; c.x = 1664;
      o = block_origin_of(c, 2, 2);   /* P=512 */
      CHECK(o.z == 4096 && o.y == 2048 && o.x == 1536, "block_origin_of floor"); }

    /* Plateau regression: P=4096 keeps four boundary-straddling nodes in four
     * bins, but P=8192 groups them. A no-progress stop at the first level would
     * incorrectly strand the hierarchy without an apex. */
    {
        OriginList p;
        Block *b5 = NULL, *b6 = NULL;
        size_t n5 = 0, n6 = 0;
        Origin dlo, dhi;
        long long boundary[6];
        origins_init(&p);
        for (int iz = 0; iz < 2; iz++) for (int iy = 0; iy < 2; iy++) {
            Origin v;
            v.z = iz ? 12288 : 10240;
            v.y = iy ? 4096 : 2048;
            v.x = 2048;
            origins_push(&p, v);
        }
        CHECK(enumerate_level_blocks(&p, 2, 5, &b5, &n5) == 0 && n5 == 4,
              "global-anchor plateau at L5");
        CHECK(blocks_all_unary(b5, n5),
              "global-anchor plateau is unary and may be virtual");
        /* A virtual level retains the materialized L4 origins; L6 groups those
         * directly while their 2048-voxel support covers the whole crop. */
        CHECK(enumerate_level_blocks(&p, 2, 6, &b6, &n6) == 0 && n6 == 1,
              "plateau collapses at L6");
        dlo.z = 11904; dlo.y = 3712; dlo.x = 2432;
        dhi.z = 13184; dhi.y = 4992; dhi.x = 3712;
        if (n6 == 1) {
            block_boundary_box(&b6[0], 2048, dlo, dhi, boundary);
            CHECK(boundary[0] == 11904 && boundary[1] == 13184 &&
                  boundary[2] == 3712 && boundary[3] == 4992 &&
                  boundary[4] == 2432 && boundary[5] == 3712,
                  "virtual plateau preserves full occupied boundary box");
        }
        free(b5); free(b6); origins_free(&p);
    }

    o.z = 4352; o.y = 3200; o.x = 2816;
    block_id_str(buf, sizeof(buf), o);
    CHECK(strcmp(buf, "z04352_y03200_x02816") == 0, "block_id_str");
    { Origin p; CHECK(parse_origin(buf, &p) == 0 && origin_eq(p, o), "parse_origin round-trip"); }

    { char nb[512]; node_obj_path(nb, sizeof(nb), "d", "z00001_y00002_x00003", "L1");
      CHECK(strcmp(nb, "d/z00001_y00002_x00003/"
                       "z00001_y00002_x00003_L1_all.obj") == 0, "node_obj_path template"); }

    /* Headline: per-level node counts pin the global-anchoring math. */
    if (build_full_grid(&cur) != 0) { fprintf(stderr, "FAIL: build_full_grid\n"); return 1; }
    CHECK(cur.n == 1764, "full grid = 1764");
    for (L = 1; cur.n > 1 && ei < (int)(sizeof(expect) / sizeof(expect[0])); L++, ei++) {
        Block *blocks = NULL; size_t nblk = 0; OriginList next; size_t b;
        if (enumerate_level_blocks(&cur, 2, L, &blocks, &nblk) != 0) { fails++; break; }
        if ((int)nblk != expect[ei]) {
            fprintf(stderr, "FAIL: L%d nodes = %zu, expected %d\n", L, nblk, expect[ei]);
            fails++;
        }
        origins_init(&next);
        for (b = 0; b < nblk; b++) origins_push(&next, blocks[b].origin);
        free(blocks);
        origins_free(&cur);
        cur = next;
    }
    origins_free(&cur);

    #undef CHECK
    if (fails == 0) fprintf(stderr, "hierarchical_weld selftest: ALL PASS\n");
    else fprintf(stderr, "hierarchical_weld selftest: %d FAILURE(S)\n", fails);
    return fails ? 1 : 0;
}

/* ================================================================
 * Argument parsing + main
 * ================================================================ */
static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s <leaf_dump_dir> <output_dir> [options]\n"
        "  --leaf-stage NAME     leaf OBJ stage (default step12_final)\n"
        "  --fanout N            block edge in nodes (default 2 = 2x2x2)\n"
        "  --keep-ratio R        decimate to R of faces each level (default 0.25)\n"
        "  --no-decimate         weld-only every level, NO qslim (preserves band-CVT\n"
        "                        quality; face count stays ~constant up the pyramid,\n"
        "                        so pair with a low --max-concurrent for memory)\n"
        "  --cvt-simplify        decimate each level with cvt_simplify (CVT per-\n"
        "                        component, boundary-preserving) instead of qslim:\n"
        "                        a mip pyramid that KEEPS CVT quality at every tier\n"
        "  --topology-first      no decimation/holefill/recoarsen/band-CVT; after\n"
        "                        every weld iteratively cut supported winding\n"
        "                        shortcuts and require all exact/topology gates\n"
        "  --axis-table PATH     curved umbilicus z,y,x CSV (required above)\n"
        "  --repair-passes N     maximum supported-cut passes/node (default 24)\n"
        "  --shortcut-tol F      persistent local shortcut threshold (default .58)\n"
        "  --max-levels N        safety cap on levels (default 16)\n"
        "  --stop-nodes N        stop when a level has <= N nodes (default 1)\n"
        "  --max-concurrent K    concurrent blocks (default 16)\n"
        "  --umbilicus-y F       winding-gate umbilicus y (default 3405)\n"
        "  --umbilicus-x F       winding-gate umbilicus x (default 2878)\n"
        "  --wrap-pitch F        winding-gate wrap pitch (default 9.5)\n"
        "  --block-timeout S     per-block wall-clock timeout sec (default 900)\n"
        "  --grid-weld PATH      grid_weld.exe (default build/Release/grid_weld.exe)\n"
        "  --qslim PATH          qslim_obj.exe (default build/Release/qslim_obj.exe)\n"
        "  --manifold PATH       manifold_check.exe (default build/Release/manifold_check.exe)\n"
        "  --obj-reorient PATH    obj_reorient.exe (default build/Release/obj_reorient.exe)\n"
        "  --wind-audit PATH     wind_audit.exe (default build/Release/wind_audit.exe)\n"
        "  --seam-audit PATH     seam_audit.exe (default build/Release/seam_audit.exe)\n"
        "  --skip-existing       resume: skip blocks whose node OBJ is complete\n"
        "  --simplify-top        also decimate the terminal level (+ re-orient)\n"
        "  --no-remesh           skip the post-decimation isotropic remesh pass\n"
        "  --dry-run             print the pyramid + per-block argv, then exit\n"
        "  --selftest            run unit tests and exit\n",
        prog);
}

static int parse_args(int argc, char **argv, Options *o)
{
    int i = 0;
    memset(o, 0, sizeof(*o));
    o->leaf_stage = "step12_final";
    o->fanout = 2;
    o->keep_ratio = 0.25;
    o->max_levels = 16;
    o->stop_nodes = 1;
    o->max_concurrent = 16;
    o->block_timeout_sec = 900.0;
    o->repair_passes = 24;
    o->shortcut_tol = 0.58;
    o->umb_y = 3405.0; o->umb_x = 2878.0; o->wrap_pitch = 9.5;
    o->grid_weld_exe = "build/Release/grid_weld.exe";
    o->qslim_exe = "build/Release/qslim_obj.exe";
    o->cvt_simplify_exe = "build/Release/cvt_simplify.exe";
    o->manifold_exe = "build/Release/manifold_check.exe";
    o->obj_reorient_exe = "build/Release/obj_reorient.exe";
    o->wind_audit_exe = "build/Release/wind_audit.exe";
    o->seam_audit_exe = "build/Release/seam_audit.exe";
    o->remesh_exe = "build/Release/remesh_obj.exe";
    o->remesh = 1;   /* isotropic remesh quality pass on by default */

    if (argc < 3) { usage(argv[0]); return -1; }
    o->leaf_dir = argv[1];
    o->out_dir = argv[2];

    for (i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--leaf-stage") && i + 1 < argc) o->leaf_stage = argv[++i];
        else if (!strcmp(argv[i], "--fanout") && i + 1 < argc) o->fanout = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--keep-ratio") && i + 1 < argc) o->keep_ratio = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max-levels") && i + 1 < argc) o->max_levels = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--stop-nodes") && i + 1 < argc) o->stop_nodes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--max-concurrent") && i + 1 < argc) o->max_concurrent = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--umbilicus-y") && i + 1 < argc) o->umb_y = atof(argv[++i]);
        else if (!strcmp(argv[i], "--umbilicus-x") && i + 1 < argc) o->umb_x = atof(argv[++i]);
        else if (!strcmp(argv[i], "--wrap-pitch") && i + 1 < argc) o->wrap_pitch = atof(argv[++i]);
        else if (!strcmp(argv[i], "--block-timeout") && i + 1 < argc) o->block_timeout_sec = atof(argv[++i]);
        else if (!strcmp(argv[i], "--grid-weld") && i + 1 < argc) o->grid_weld_exe = argv[++i];
        else if (!strcmp(argv[i], "--qslim") && i + 1 < argc) o->qslim_exe = argv[++i];
        else if (!strcmp(argv[i], "--manifold") && i + 1 < argc) o->manifold_exe = argv[++i];
        else if (!strcmp(argv[i], "--obj-reorient") && i + 1 < argc) o->obj_reorient_exe = argv[++i];
        else if (!strcmp(argv[i], "--remesh-exe") && i + 1 < argc) o->remesh_exe = argv[++i];
        else if (!strcmp(argv[i], "--skip-existing") || !strcmp(argv[i], "--resume")) o->skip_existing = 1;
        else if (!strcmp(argv[i], "--simplify-top")) o->simplify_top = 1;
        else if (!strcmp(argv[i], "--no-decimate")) o->no_decimate = 1;
        else if (!strcmp(argv[i], "--cvt-simplify")) o->cvt_simplify = 1;
        else if (!strcmp(argv[i], "--cvt-simplify-exe") && i + 1 < argc) o->cvt_simplify_exe = argv[++i];
        else if (!strcmp(argv[i], "--topology-first")) o->topology_first = 1;
        else if (!strcmp(argv[i], "--axis-table") && i + 1 < argc) o->axis_table = argv[++i];
        else if (!strcmp(argv[i], "--repair-passes") && i + 1 < argc) o->repair_passes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-tol") && i + 1 < argc) o->shortcut_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--wind-audit") && i + 1 < argc) o->wind_audit_exe = argv[++i];
        else if (!strcmp(argv[i], "--seam-audit") && i + 1 < argc) o->seam_audit_exe = argv[++i];
        else if (!strcmp(argv[i], "--no-remesh")) o->remesh = 0;
        else if (!strcmp(argv[i], "--dry-run")) o->dry_run = 1;
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); usage(argv[0]); return -1; }
    }
    if (o->fanout < 2) { fprintf(stderr, "fanout must be >= 2\n"); return -1; }
    if (!(o->keep_ratio > 0.0 && o->keep_ratio < 1.0)) { fprintf(stderr, "keep-ratio must be in (0,1)\n"); return -1; }
    if (o->stop_nodes < 1) o->stop_nodes = 1;
    if (o->repair_passes < 1 || !(o->shortcut_tol > 0.0)) {
        fprintf(stderr, "repair-passes and shortcut-tol must be positive\n");
        return -1;
    }
    if (o->topology_first && o->axis_table == NULL) {
        fprintf(stderr, "--topology-first requires --axis-table z,y,x.csv\n");
        return -1;
    }
    if (o->topology_first) {
        o->no_decimate = 1;
        o->remesh = 0;
    }
    return 0;
}

int main(int argc, char **argv)
{
    Options opt;
    OriginList cur;
    const char *cur_dir = NULL, *cur_stage = NULL;
    long long cur_pitch = 128;
    char envbuf[64];
    int L = 0, total_fail = 0;
    Origin domain_lo, domain_hi;

    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return selftest();
    if (parse_args(argc, argv, &opt) != 0) return 1;

    /* Arm the winding gate for every child weld (inherited env). */
    snprintf(envbuf, sizeof(envbuf), "%lld", (long long)(opt.umb_y + 0.5)); set_env("SEAM_UMBILICUS_Y", envbuf);
    snprintf(envbuf, sizeof(envbuf), "%lld", (long long)(opt.umb_x + 0.5)); set_env("SEAM_UMBILICUS_X", envbuf);
    snprintf(envbuf, sizeof(envbuf), "%.4f", opt.wrap_pitch);               set_env("SEAM_WRAP_PITCH", envbuf);
    if (opt.topology_first) {
        set_env("SEAM_WIND_TOL", "0.25");
        set_env("SEAM_WIND_HARD_TOL", "0.40");
        set_env("SEAM_NO_HOLEFILL", "1");
        set_env("SEAM_NO_RECOARSEN", "1");
        set_env("SEAM_NO_BANDCVT", "1");
        fprintf(stderr,
            "topology-first: curved-axis bridge + iterative supported cuts; "
            "holefill/recoarsen/band-CVT/decimation disabled\n");
    }

    if (ensure_dir(opt.out_dir) != 0) { fprintf(stderr, "cannot create %s\n", opt.out_dir); return 1; }

    origins_init(&cur);
    if (enumerate_node_dirs(opt.leaf_dir, &cur) != 0 || cur.n == 0) {
        fprintf(stderr, "No leaf nodes found under %s\n", opt.leaf_dir);
        origins_free(&cur);
        return 1;
    }
    fprintf(stderr, "hierarchical_weld: %zu leaf nodes under %s\n", cur.n, opt.leaf_dir);

    domain_lo = cur.o[0];
    domain_hi = cur.o[0];
    for (size_t k = 1; k < cur.n; k++) {
        if (cur.o[k].z < domain_lo.z) domain_lo.z = cur.o[k].z;
        if (cur.o[k].y < domain_lo.y) domain_lo.y = cur.o[k].y;
        if (cur.o[k].x < domain_lo.x) domain_lo.x = cur.o[k].x;
        if (cur.o[k].z > domain_hi.z) domain_hi.z = cur.o[k].z;
        if (cur.o[k].y > domain_hi.y) domain_hi.y = cur.o[k].y;
        if (cur.o[k].x > domain_hi.x) domain_hi.x = cur.o[k].x;
    }
    domain_hi.z += level_pitch(opt.fanout, 0);
    domain_hi.y += level_pitch(opt.fanout, 0);
    domain_hi.x += level_pitch(opt.fanout, 0);

    preview_pyramid(&opt, &cur);
    if (opt.dry_run) { origins_free(&cur); return 0; }

    cur_dir = opt.leaf_dir;
    cur_stage = opt.leaf_stage;

    for (L = 1; cur.n > (size_t)opt.stop_nodes && L <= opt.max_levels; L++) {
        Block *blocks = NULL;
        size_t nblk = 0;
        size_t prev_n = cur.n;
        char this_dir[1024], this_stage[32];
        LevelCtx lc;
        OriginList next;
        int terminal = 0, rc = 0;

        if (enumerate_level_blocks(&cur, opt.fanout, L, &blocks, &nblk) != 0) {
            fprintf(stderr, "enumerate_level_blocks failed at L%d\n", L);
            break;
        }
        if (nblk == 0) { free(blocks); break; }
        if (nblk > prev_n) {
            fprintf(stderr, "L%d: invalid expansion (%zu -> %zu nodes) -- stopping\n",
                    L, prev_n, nblk);
            free(blocks);
            break;
        }

        /* Do not materialize a pure unary plateau unless it is the final level
         * allowed by the caller.  Keeping cur_dir/cur_stage/cur_pitch unchanged
         * makes the next real weld consume the certified children directly and
         * detect their seams on the correct (last materialized) lattice. */
        if (nblk == prev_n && blocks_all_unary(blocks, nblk) &&
            L < opt.max_levels) {
            fprintf(stderr,
                    "\n=== Level %d: %zu unary blocks -- virtual plateau "
                    "(retain cube-size %lld) ===\n",
                    L, nblk, (long long)cur_pitch);
            free(blocks);
            continue;
        }

        terminal = (nblk <= (size_t)opt.stop_nodes);
        snprintf(this_dir, sizeof(this_dir), "%s/level%d", opt.out_dir, L);
        snprintf(this_stage, sizeof(this_stage), "L%d", L);

        lc.opt = &opt;
        lc.level = L;
        lc.prev_dir = cur_dir;
        lc.child_stage = cur_stage;
        lc.this_dir = this_dir;
        lc.this_stage = this_stage;
        lc.child_pitch = cur_pitch;
        lc.domain_lo = domain_lo;
        lc.domain_hi = domain_hi;
        lc.do_simplify = ((!terminal) || opt.simplify_top) && !opt.no_decimate;

        origins_init(&next);
        rc = run_level(&lc, blocks, nblk, &next);
        total_fail += rc;

        if (terminal && opt.simplify_top) {
            size_t b;
            for (b = 0; b < next.n; b++) reorient_node(&opt, this_dir, this_stage, next.o[b]);
        }

        free(blocks);
        origins_free(&cur);
        cur = next;
        cur_pitch = level_pitch(opt.fanout, L);

        /* this_dir/this_stage live on the stack; copy into heap-stable storage
         * for the next iteration's prev_dir/prev_stage. */
        {
            static char kept_dir[1024];
            static char kept_stage[32];
            snprintf(kept_dir, sizeof(kept_dir), "%s", this_dir);
            snprintf(kept_stage, sizeof(kept_stage), "%s", this_stage);
            cur_dir = kept_dir;
            cur_stage = kept_stage;
        }
    }

    fprintf(stderr, "\n=== Pyramid complete: top level has %zu node(s), %d block failure(s) ===\n",
            cur.n, total_fail);
    if (cur.n >= 1) {
        char id[128], node[1024];
        block_id_str(id, sizeof(id), cur.o[0]);
        snprintf(node, sizeof(node), "%s/level%d/%s/%s_L%d_all.obj",
                 opt.out_dir, L - 1, id, id, L - 1);
        fprintf(stderr, "Top node: %s\n", node);
    }

    origins_free(&cur);
    return total_fail > 0 ? 2 : 0;
}
