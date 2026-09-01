/* ============================================================================
 * quadribbon_pipeline.c -- ONE command: pre-weld mesh pile in, big sheet out.
 *
 *   quadribbon <mesh_dir> <out_dir> [config.json] [--stop-after <stage>]
 *   quadribbon --selftest
 *
 * --stop-after concat|fit|reopt|untangle|reopt2|sheet (or 1..6) finishes the
 * named stage (including its debug bake) and exits cleanly, so an earlier
 * stage can be debugged without anything downstream running on top of it.
 * Re-running without the flag resumes from the completed artifacts.
 *
 * The canonical unwrap/flatten order of operations, in one tool:
 *   stage 1  concat  -- scan <mesh_dir> recursively for the per-cube pipeline's
 *            *_final_all.vmesh dumps (plus plain *.vmesh directly in the dir),
 *            deduplicate by cube id, and concatenate into one world-frame soup.
 *            NO WELDING: cross-cube continuity is the fit's job.
 *   stage 2  fit     -- scroll_ribbon claims-mode quadribbon fit (spawned).
 *   stage 3  reopt   -- scroll_ribbon --preserve-input-topology metric
 *            projection + TAUCS consistency solve (spawned).
 *   stage 4  untangle-- collision shell on the parameterized quadribbon
 *            (in-process; radial displacement only, topology/UV unchanged).
 *   stage 5  reopt   -- the same metric projection on the untangled geometry.
 *   stage 6  sheet   -- layer strips (spawned) -> RAW bakes (spawned) ->
 *            first-cover composite (in-process).  big_sheet_{pre,post,
 *            provenance,pre_vs_post}.png land at the out_dir root; "pre" is
 *            the pre-untangle parameterization, "post" the final one.
 *
 * All policy comes from the production config JSON + pipeline_constants.h;
 * there are no tuning flags.  Stages resume: a stage whose primary artifact
 * already exists is skipped, so a crashed run re-enters where it stopped.
 * ==========================================================================*/
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/mesh_bin.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include "../flatten/quadribbon_untangle.h"
#include "../flatten/sheet_composite.h"

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

enum { QP_MAX_PATH = 2048, QP_MAX_PILE = 65536 };

static const char QP_DEFAULT_CONFIG[] =
    "configs/flatten/pherc0139.production.json";

typedef struct {
    double axis_y, axis_x;
    double wrap_spacing;
    char raw_path[QP_MAX_PATH];
    long raw_chunk;
    long window_lo, window_hi, dark;
    int allow_missing_chunks;
} QpConfig;

/* ---- small path/file helpers --------------------------------------------- */

static void qp_join(char *out, size_t cap, const char *a, const char *b)
{
    snprintf(out, cap, "%s/%s", a, b);
}

/* ---- pipeline log: every milestone goes to stderr AND logs/pipeline.log,
 * prefixed with elapsed seconds, so a hand-run session leaves a timeline
 * even after the console scrolls away. ------------------------------------ */

static char qp_pipeline_log[QP_MAX_PATH] = "";
static double qp_log_t0 = 0.0;

static void qp_logf(const char *fmt, ...)
{
    va_list ap;
    char stamp[32];
    snprintf(stamp, sizeof stamp, "[%8.1fs] ",
             qp_log_t0 > 0.0 ? ves_clock_sec() - qp_log_t0 : 0.0);
    fputs(stamp, stderr);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fflush(stderr);
    if (qp_pipeline_log[0] != '\0') {
        FILE *f = fopen(qp_pipeline_log, "ab");
        if (f != NULL) {
            fputs(stamp, f);
            va_start(ap, fmt);
            vfprintf(f, fmt, ap);
            va_end(ap);
            fclose(f);
        }
    }
}

/* Stage ladder for --stop-after: name (or 1-based digit) -> stage index. */
static int qp_stage_index(const char *name)
{
    static const char *names[6] = {
        "concat", "fit", "reopt", "untangle", "reopt2", "sheet"
    };
    if (name == NULL) return 0;
    for (int i = 0; i < 6; i++)
        if (strcmp(name, names[i]) == 0) return i + 1;
    if (strlen(name) == 1 && name[0] >= '1' && name[0] <= '6')
        return name[0] - '0';
    return -1;
}

static int qp_exists(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL) return 0;
    fclose(f);
    return 1;
}

static int qp_mkdir(const char *path)
{
#ifdef _WIN32
    return _mkdir(path) == 0 || GetLastError() == ERROR_ALREADY_EXISTS ? 0 : -1;
#else
    return mkdir(path, 0775) == 0 || errno == EEXIST ? 0 : -1;
#endif
}

static int qp_copy_file(const char *src, const char *dst)
{
    FILE *in = fopen(src, "rb");
    FILE *out = NULL;
    char buf[1 << 16];
    size_t n = 0;
    if (in == NULL) return -1;
    out = fopen(dst, "wb");
    if (out == NULL) { fclose(in); return -1; }
    while ((n = fread(buf, 1, sizeof buf, in)) > 0) {
        if (fwrite(buf, 1, n, out) != n) {
            fclose(in); fclose(out); return -1;
        }
    }
    fclose(in);
    if (fclose(out) != 0) return -1;
    return 0;
}

/* ---- config -------------------------------------------------------------- */

static int qp_load_config(Arena_T arena, const char *path, QpConfig *cfg)
{
    const char *err = NULL;
    const JsonValue *root = Json_parse_file(arena, path, &err);
    const JsonValue *geom = NULL, *raw = NULL, *bake = NULL, *axis = NULL;
    const JsonValue *win = NULL;
    memset(cfg, 0, sizeof *cfg);
    if (root == NULL) {
        fprintf(stderr, "quadribbon: cannot parse config %s (%s)\n",
                path, err != NULL ? err : "io error");
        return -1;
    }
    geom = Json_object_get(root, "geometry");
    raw = Json_object_get(root, "raw_source");
    bake = Json_object_get(root, "bake");
    if (geom == NULL || raw == NULL || bake == NULL) {
        fprintf(stderr, "quadribbon: config missing geometry/raw_source/bake "
                "sections: %s\n", path);
        return -1;
    }
    axis = Json_object_get(geom, "cylindrical_axis_yx");
    if (axis == NULL || Json_array_len(axis) != 2) {
        fprintf(stderr, "quadribbon: config missing cylindrical_axis_yx\n");
        return -1;
    }
    cfg->axis_y = Json_as_double(Json_array_get(axis, 0), 0.0);
    cfg->axis_x = Json_as_double(Json_array_get(axis, 1), 0.0);
    cfg->wrap_spacing = Json_member_double(geom, "wrap_spacing_voxels", 0.0);
    {
        const char *rp = Json_as_string(Json_object_get(raw, "path"));
        if (rp == NULL) {
            fprintf(stderr, "quadribbon: config missing raw_source.path\n");
            return -1;
        }
        snprintf(cfg->raw_path, sizeof cfg->raw_path, "%s", rp);
    }
    cfg->raw_chunk = Json_member_long(raw, "chunk_size", 128);
    cfg->allow_missing_chunks =
        Json_as_bool(Json_object_get(raw, "allow_missing_chunks"), 0);
    win = Json_object_get(bake, "window_u8");
    if (win == NULL || Json_array_len(win) != 2) {
        fprintf(stderr, "quadribbon: config missing bake.window_u8\n");
        return -1;
    }
    cfg->window_lo = Json_as_long(Json_array_get(win, 0), 47);
    cfg->window_hi = Json_as_long(Json_array_get(win, 1), 193);
    cfg->dark = Json_member_long(bake, "dark_threshold_u8", 13);
    if (!(cfg->wrap_spacing > 0.0) || cfg->window_lo >= cfg->window_hi) {
        fprintf(stderr, "quadribbon: config geometry/window invalid\n");
        return -1;
    }
    return 0;
}

/* ---- stage 1: pile scan + concat ----------------------------------------- */

typedef struct {
    char path[QP_MAX_PATH];
    char cube_id[24];          /* "" for plain world-frame meshes */
    long oz, oy, ox;
} QpPileEntry;

static int qp_parse_cube_id(const char *name, char id[24],
                            long *oz, long *oy, long *ox)
{
    /* z#####_y#####_x##### prefix, zero-padded 5 digits */
    long z = 0, y = 0, x = 0;
    if (strlen(name) < 20) return -1;
    if (name[0] != 'z' || name[6] != '_' || name[7] != 'y' ||
        name[13] != '_' || name[14] != 'x')
        return -1;
    if (sscanf(name, "z%5ld_y%5ld_x%5ld", &z, &y, &x) != 3) return -1;
    memcpy(id, name, 20);
    id[20] = '\0';
    *oz = z; *oy = y; *ox = x;
    return 0;
}

static int qp_pile_push(QpPileEntry *pile, size_t *n, const char *dir,
                        const char *name)
{
    size_t len = strlen(name);
    QpPileEntry *e = NULL;
    char id[24] = "";
    long oz = 0, oy = 0, ox = 0;
    int has_id = qp_parse_cube_id(name, id, &oz, &oy, &ox) == 0;
    if (len < 7 || strcmp(name + len - 6, ".vmesh") != 0) return 0;
    /* per-cube final dumps carry the id AND the _final_all marker; a plain
     * world-frame mesh has neither restriction */
    if (has_id && strstr(name, "_final_all.vmesh") == NULL) return 0;
    if (has_id) {
        for (size_t i = 0; i < *n; i++)
            if (strcmp(pile[i].cube_id, id) == 0) return 0;  /* dedup by cube */
    }
    if (*n >= QP_MAX_PILE) return -1;
    e = &pile[(*n)++];
    qp_join(e->path, sizeof e->path, dir, name);
    memcpy(e->cube_id, id, sizeof id);
    e->oz = oz; e->oy = oy; e->ox = ox;
    return 0;
}

static int qp_scan_dir(const char *dir, int depth, QpPileEntry *pile,
                       size_t *n)
{
    if (depth > 6) return 0;
#ifdef _WIN32
    char glob[QP_MAX_PATH];
    WIN32_FIND_DATAA fd;
    HANDLE h;
    snprintf(glob, sizeof glob, "%s/*", dir);
    h = FindFirstFileA(glob, &fd);
    if (h == INVALID_HANDLE_VALUE) return -1;
    do {
        if (strcmp(fd.cFileName, ".") == 0 || strcmp(fd.cFileName, "..") == 0)
            continue;
        if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            char sub[QP_MAX_PATH];
            qp_join(sub, sizeof sub, dir, fd.cFileName);
            if (qp_scan_dir(sub, depth + 1, pile, n) != 0) {
                FindClose(h);
                return -1;
            }
        } else if (depth == 0 || strstr(fd.cFileName, "_final_all.vmesh")) {
            if (qp_pile_push(pile, n, dir, fd.cFileName) != 0) {
                FindClose(h);
                return -1;
            }
        }
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(dir);
    struct dirent *de = NULL;
    if (d == NULL) return -1;
    while ((de = readdir(d)) != NULL) {
        char sub[QP_MAX_PATH];
        struct stat st;
        if (strcmp(de->d_name, ".") == 0 || strcmp(de->d_name, "..") == 0)
            continue;
        qp_join(sub, sizeof sub, dir, de->d_name);
        if (stat(sub, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) {
            if (qp_scan_dir(sub, depth + 1, pile, n) != 0) {
                closedir(d);
                return -1;
            }
        } else if (depth == 0 || strstr(de->d_name, "_final_all.vmesh")) {
            if (qp_pile_push(pile, n, dir, de->d_name) != 0) {
                closedir(d);
                return -1;
            }
        }
    }
    closedir(d);
#endif
    return 0;
}

static int qp_stage_concat(Arena_T arena, const char *mesh_dir,
                           const char *soup_path, const char *table_path)
{
    QpPileEntry *pile = NULL;
    size_t np = 0, total_v = 0, total_f = 0, world_frame = 0, offset_used = 0;
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    FILE *table = NULL;
    if (MeshBin_looks_complete(soup_path)) {
        fprintf(stderr, "[concat] resume: %s already complete%s\n", soup_path,
                table_path != NULL && !qp_exists(table_path)
                    ? " (cube table missing; delete the soup to regenerate)"
                    : "");
        return 0;
    }
    pile = (QpPileEntry *)calloc(QP_MAX_PILE, sizeof *pile);
    if (pile == NULL) return -1;
    if (qp_scan_dir(mesh_dir, 0, pile, &np) != 0 || np == 0) {
        fprintf(stderr, "[concat] no input meshes under %s\n", mesh_dir);
        free(pile);
        return -1;
    }
    /* two passes so peak memory stays soup + one cube: pass 1 sizes only */
    for (size_t i = 0; i < np; i++) {
        MeshBinData m;
        Arena_Mark mark = Arena_save(arena);
        memset(&m, 0, sizeof m);
        if (MeshBin_read_arena(arena, pile[i].path, &m) != 0) {
            fprintf(stderr, "[concat] cannot read %s\n", pile[i].path);
            free(pile);
            return -1;
        }
        total_v += m.nv;
        total_f += m.nf;
        Arena_restore(arena, mark);
    }
    verts = (float *)ARENA_ALLOC(arena, total_v * 3 * sizeof(float));
    faces = (int32_t *)ARENA_ALLOC(arena, total_f * 3 * sizeof(int32_t));
    if (table_path != NULL) {
        table = fopen(table_path, "wb");
        if (table != NULL)
            fprintf(table, "cube_id\tnv\tnf\tframe\toz\toy\tox\tpath\n");
    }
    for (size_t i = 0; i < np; i++) {
        MeshBinData m;
        Arena_Mark mark = Arena_save(arena);
        double oz = 0.0, oy = 0.0, ox = 0.0;
        const char *frame = "plain";
        memset(&m, 0, sizeof m);
        if (MeshBin_read_arena(arena, pile[i].path, &m) != 0) {
            if (table != NULL) fclose(table);
            free(pile);
            return -1;
        }
        if (pile[i].cube_id[0] != '\0' && m.nv > 0) {
            /* per-cube dumps are usually already world-frame; apply the cube
             * offset only when the geometry is clearly cube-local */
            float zmax = 0.0f;
            for (size_t k = 0; k < m.nv; k++)
                if (m.verts[k * 3] > zmax) zmax = m.verts[k * 3];
            if (zmax < 512.0f && pile[i].oz >= 512) {
                oz = (double)pile[i].oz;
                oy = (double)pile[i].oy;
                ox = (double)pile[i].ox;
                offset_used++;
                frame = "offset";
            } else {
                world_frame++;
                frame = "world";
            }
        }
        if (table != NULL)
            fprintf(table, "%s\t%zu\t%zu\t%s\t%ld\t%ld\t%ld\t%s\n",
                    pile[i].cube_id[0] != '\0' ? pile[i].cube_id : "-",
                    m.nv, m.nf, frame, pile[i].oz, pile[i].oy, pile[i].ox,
                    pile[i].path);
        for (size_t k = 0; k < m.nv; k++) {
            verts[(nv + k) * 3 + 0] = (float)(m.verts[k * 3 + 0] + oz);
            verts[(nv + k) * 3 + 1] = (float)(m.verts[k * 3 + 1] + oy);
            verts[(nv + k) * 3 + 2] = (float)(m.verts[k * 3 + 2] + ox);
        }
        for (size_t k = 0; k < m.nf * 3; k++)
            faces[nf * 3 + k] = m.faces[k] + (int32_t)nv;
        nv += m.nv;
        nf += m.nf;
        Arena_restore(arena, mark);
        /* NOTE: the restore frees the cube read but ALSO everything after the
         * mark; verts/faces were allocated before the loop so they survive */
    }
    if (table != NULL) {
        fclose(table);
        fprintf(stderr, "[concat] per-cube table: %s\n", table_path);
    }
    fprintf(stderr,
            "[concat] %zu mesh(es) (%zu world-frame, %zu offset-applied) -> "
            "%zu verts, %zu faces\n",
            np, world_frame, offset_used, nv, nf);
    free(pile);
    if (nv == 0 || nf == 0) {
        fprintf(stderr, "[concat] empty soup\n");
        return -1;
    }
    return MeshBin_write(soup_path, verts, nv, faces, nf, NULL);
}

/* ---- per-stage debug bake: EVERY stage output renders the barycentric
 * sampled volume texture onto its UV coordinates (obj_bake_raw), so each
 * stage is judged by the actual sheet texture, not by proxies.  A root-level
 * copy of the rawtex PNG makes browsing one flat directory sufficient. ---- */

static int qp_debug_bake(const QpConfig *cfg, const char *vmesh,
                         const char *out_root, const char *tag);

/* untangle displacement: |radius change| painted over UV, 0..8 vox -> gray */
static int qp_dump_dr_png(const float *before, const float *after,
                          const float *uv, size_t nv, double axis_y,
                          double axis_x, const char *png)
{
    double ulo = 1e300, uhi = -1e300, vlo = 1e300, vhi = -1e300;
    for (size_t i = 0; i < nv; i++) {
        double u = uv[i * 2], v = uv[i * 2 + 1];
        if (u < ulo) ulo = u;
        if (u > uhi) uhi = u;
        if (v < vlo) vlo = v;
        if (v > vhi) vhi = v;
    }
    if (!(uhi > ulo) || !(vhi > vlo)) return -1;
    {
        double span_u = uhi - ulo, span_v = vhi - vlo;
        double scale = 4096.0 / (span_u > span_v ? span_u : span_v);
        size_t W = (size_t)(span_u * scale) + 2;
        size_t H = (size_t)(span_v * scale) + 2;
        uint8_t *img = (uint8_t *)calloc(W * H, 1);
        int rc = -1;
        if (img != NULL) {
            for (size_t i = 0; i < nv; i++) {
                double rb = hypot((double)before[i * 3 + 1] - axis_y,
                                  (double)before[i * 3 + 2] - axis_x);
                double ra = hypot((double)after[i * 3 + 1] - axis_y,
                                  (double)after[i * 3 + 2] - axis_x);
                double d = fabs(ra - rb);
                uint8_t g = d <= 0.0 ? 16
                          : (uint8_t)(32.0 + (d > 8.0 ? 8.0 : d) * 27.0);
                size_t x = (size_t)((uv[i * 2] - ulo) * scale);
                size_t y = (size_t)((uv[i * 2 + 1] - vlo) * scale);
                if (x < W && y < H && g > img[y * W + x])
                    img[y * W + x] = g;
            }
            rc = VesPng_write_gray(png, img, (int)W, (int)H);
            fprintf(stderr, "[debug] %s (untangle |dr| map)\n", png);
        }
        free(img);
        return rc;
    }
}

/* ---- spawn helper -------------------------------------------------------- */

static char qp_exe_dir[QP_MAX_PATH] = "";

static void qp_tool_path(char *out, size_t cap, const char *tool)
{
#ifdef _WIN32
    const char *suffix = ".exe";
#else
    const char *suffix = "";
#endif
    if (qp_exe_dir[0] != '\0')
        snprintf(out, cap, "%s/%s%s", qp_exe_dir, tool, suffix);
    else
        snprintf(out, cap, "%s%s", tool, suffix);
}

static char qp_log_dir[QP_MAX_PATH] = "";

/* A run's first spawn into a given log file truncates it: append-only logs
 * interleaved runs and made forensics unreadable.  Later spawns of the same
 * tool within THIS run still append (multi-spawn stages keep their order). */
static char qp_seen_logs[32][QP_MAX_PATH];
static size_t qp_n_seen_logs = 0;

static void qp_fresh_log(const char *log_path)
{
    for (size_t i = 0; i < qp_n_seen_logs; i++)
        if (strcmp(qp_seen_logs[i], log_path) == 0) return;
    remove(log_path);
    if (qp_n_seen_logs < sizeof qp_seen_logs / sizeof qp_seen_logs[0])
        snprintf(qp_seen_logs[qp_n_seen_logs++], QP_MAX_PATH, "%s", log_path);
}

static int qp_spawn(const char *tool, const char *const *argv)
{
    double t0 = ves_clock_sec();
    char log_path[QP_MAX_PATH];
    int rc = 0;
    snprintf(log_path, sizeof log_path, "%s/%s.log", qp_log_dir, tool);
    for (char *p = log_path + strlen(qp_log_dir) + 1; *p; p++)
        if (*p == '(' || *p == ')' || *p == ' ') *p = '_';
    if (qp_log_dir[0]) qp_fresh_log(log_path);
    {
        char line[QP_MAX_PATH];
        size_t off = 0;
        off += (size_t)snprintf(line + off, sizeof line - off, "%s", tool);
        for (size_t i = 1; argv[i] != NULL && off + 2 < sizeof line; i++)
            off += (size_t)snprintf(line + off, sizeof line - off, " %s",
                                    argv[i]);
        qp_logf("[spawn] %s\n  log: %s\n", line, log_path);
    }
    rc = ves_run_subprocess_logged(argv[0], argv, 0.0,
                                   qp_log_dir[0] ? log_path : NULL);
    qp_logf("[spawn] %s -> rc=%d (%.1fs)\n", tool, rc, ves_clock_sec() - t0);
    return rc;
}

static int qp_debug_bake(const QpConfig *cfg, const char *vmesh,
                         const char *out_root, const char *tag)
{
    char bake_dir[QP_MAX_PATH], tool[QP_MAX_PATH];
    char root_png[QP_MAX_PATH], stage_png[QP_MAX_PATH];
    char chunk_s[16], lo_s[16], hi_s[16], dark_s[16];
    snprintf(root_png, sizeof root_png, "%s/%s_bake.png", out_root, tag);
    if (qp_exists(root_png)) {
        fprintf(stderr, "[debug-bake] resume: %s\n", root_png);
        return 0;
    }
    snprintf(bake_dir, sizeof bake_dir, "%s/debug_bakes/%s", out_root, tag);
    {
        char parent[QP_MAX_PATH];
        snprintf(parent, sizeof parent, "%s/debug_bakes", out_root);
        qp_mkdir(parent);
    }
    qp_mkdir(bake_dir);
    snprintf(chunk_s, sizeof chunk_s, "%ld", cfg->raw_chunk);
    snprintf(lo_s, sizeof lo_s, "%ld", cfg->window_lo);
    snprintf(hi_s, sizeof hi_s, "%ld", cfg->window_hi);
    snprintf(dark_s, sizeof dark_s, "%ld", cfg->dark);
    qp_tool_path(tool, sizeof tool, "obj_bake_raw");
    {
        const char *argv[18];
        size_t a = 0;
        argv[a++] = tool;
        argv[a++] = vmesh;
        argv[a++] = cfg->raw_path;
        argv[a++] = bake_dir;
        argv[a++] = "--id"; argv[a++] = tag;
        argv[a++] = "--chunk"; argv[a++] = chunk_s;
        argv[a++] = cfg->allow_missing_chunks ? "--allow-incomplete-raw"
                                              : "--require-complete-raw";
        argv[a++] = "--window"; argv[a++] = lo_s; argv[a++] = hi_s;
        argv[a++] = "--diag-dark"; argv[a++] = dark_s;
        argv[a++] = "--raster-auto";
        argv[a] = NULL;
        if (qp_spawn("obj_bake_raw", argv) != 0) {
            fprintf(stderr, "[debug-bake] %s FAILED\n", tag);
            return -1;
        }
    }
    snprintf(stage_png, sizeof stage_png, "%s/%s_rawtex.png", bake_dir, tag);
    if (qp_copy_file(stage_png, root_png) == 0)
        fprintf(stderr, "[debug-bake] %s\n", root_png);
    return 0;
}

/* ---- stage 6: strips + bakes + composite --------------------------------- */

static int qp_strip_min_u(const char *strip_vmesh, double *out_min_u)
{
    MeshBinData m;
    Arena_T arena = Arena_new();
    double lo = 1e300;
    int rc = -1;
    memset(&m, 0, sizeof m);
    if (arena != NULL &&
        MeshBin_read_arena(arena, strip_vmesh, &m) == 0 && m.uv != NULL) {
        for (size_t i = 0; i < m.nv; i++)
            if (m.uv[i * 2] < lo) lo = m.uv[i * 2];
        rc = 0;
    }
    if (arena != NULL) Arena_dispose(&arena);
    *out_min_u = lo;
    return rc;
}

/* Pick the strips collision-raster --du for a ribbon: power of two, smallest
 * such that (u_span/du) * v_span fits the tool's 2^28-pixel cap with 2x
 * headroom.  1 at the 4x5x5 scale (unchanged behavior); larger rungs coarsen
 * only the collision granularity -- output UV is untouched by --du. */
static long qp_pick_strips_du(const char *ribbon_vmesh)
{
    MeshBinData m;
    Arena_T arena = Arena_new();
    double ulo = 1e300, uhi = -1e300, vlo = 1e300, vhi = -1e300;
    long du = 1;
    memset(&m, 0, sizeof m);
    if (arena != NULL &&
        MeshBin_read_arena(arena, ribbon_vmesh, &m) == 0 && m.uv != NULL) {
        for (size_t i = 0; i < m.nv; i++) {
            double u = m.uv[i * 2], v = m.uv[i * 2 + 1];
            if (u < ulo) ulo = u;
            if (u > uhi) uhi = u;
            if (v < vlo) vlo = v;
            if (v > vhi) vhi = v;
        }
        if (uhi > ulo && vhi > vlo) {
            double vspan = vhi - vlo + 1.0;
            while (du < 4096 &&
                   ((uhi - ulo) / (double)du + 1.0) * vspan >
                       134217728.0 /* 2^27 */)
                du *= 2;
        }
    }
    if (arena != NULL) Arena_dispose(&arena);
    return du;
}

static int qp_stage_sheet(Arena_T arena, const QpConfig *cfg,
                          const char *ribbon_vmesh, const char *labels_i32,
                          const char *stage_dir, const char *tag,
                          const char *out_root, uint8_t **out_tex,
                          size_t *out_w, size_t *out_h)
{
    char strips_dir[QP_MAX_PATH], strip_prefix[QP_MAX_PATH];
    char tool[QP_MAX_PATH];
    char tex_png[QP_MAX_PATH], prov_png[QP_MAX_PATH];
    char layer_png[SHEET_COMPOSITE_MAX_LAYERS][QP_MAX_PATH];
    const char *layer_paths[SHEET_COMPOSITE_MAX_LAYERS] = { NULL, NULL, NULL };
    SheetCompositeLayer layers[SHEET_COMPOSITE_MAX_LAYERS];
    SheetCompositeStats stats;
    double min_u[SHEET_COMPOSITE_MAX_LAYERS];
    uint8_t *stretch_buf[SHEET_COMPOSITE_MAX_LAYERS] = { NULL, NULL, NULL };
    uint8_t *squash_buf[SHEET_COMPOSITE_MAX_LAYERS] = { NULL, NULL, NULL };
    size_t n_layers = 0;
    memset(layers, 0, sizeof layers);

    qp_mkdir(stage_dir);
    qp_join(strips_dir, sizeof strips_dir, stage_dir, "strips");
    qp_mkdir(strips_dir);
    qp_join(strip_prefix, sizeof strip_prefix, strips_dir, "strip");

    {
        char manifest[QP_MAX_PATH];
        snprintf(manifest, sizeof manifest, "%s_manifest.json", strip_prefix);
        if (!qp_exists(manifest)) {
            char du[16] = "1", dv[8] = "1";
            snprintf(du, sizeof du, "%ld", qp_pick_strips_du(ribbon_vmesh));
            qp_tool_path(tool, sizeof tool, "quadribbon_strips");
            {
                const char *argv[12];
                size_t a = 0;
                argv[a++] = tool;
                argv[a++] = ribbon_vmesh;
                argv[a++] = strip_prefix;
                if (labels_i32 != NULL && qp_exists(labels_i32)) {
                    argv[a++] = "--labels";
                    argv[a++] = labels_i32;
                }
                argv[a++] = "--du"; argv[a++] = du;
                argv[a++] = "--dv"; argv[a++] = dv;
                argv[a] = NULL;
                if (qp_spawn("quadribbon_strips", argv) != 0) return -1;
            }
        } else {
            fprintf(stderr, "[sheet-%s] resume: strips already built\n", tag);
        }
    }

    for (size_t s = 0; s < SHEET_COMPOSITE_MAX_LAYERS; s++) {
        char strip_vmesh[QP_MAX_PATH], bake_dir[QP_MAX_PATH];
        char id[32], chunk_s[16], lo_s[16], hi_s[16], dark_s[16];
        char rawtex_tif[QP_MAX_PATH], cover_tif[QP_MAX_PATH];
        uint8_t *tex = NULL, *cov = NULL;
        int D = 0, H = 0, W = 0;
        snprintf(strip_vmesh, sizeof strip_vmesh, "%s_%02zu.vmesh",
                 strip_prefix, s);
        if (!qp_exists(strip_vmesh)) break;
        snprintf(id, sizeof id, "strip_%02zu", s);
        qp_join(bake_dir, sizeof bake_dir, stage_dir, id);
        qp_mkdir(bake_dir);
        snprintf(rawtex_tif, sizeof rawtex_tif, "%s/%s_rawtex.tif",
                 bake_dir, id);
        snprintf(cover_tif, sizeof cover_tif, "%s/%s_rawtex_coverage.tif",
                 bake_dir, id);
        if (!qp_exists(rawtex_tif) || !qp_exists(cover_tif)) {
            snprintf(chunk_s, sizeof chunk_s, "%ld", cfg->raw_chunk);
            snprintf(lo_s, sizeof lo_s, "%ld", cfg->window_lo);
            snprintf(hi_s, sizeof hi_s, "%ld", cfg->window_hi);
            snprintf(dark_s, sizeof dark_s, "%ld", cfg->dark);
            qp_tool_path(tool, sizeof tool, "obj_bake_raw");
            {
                const char *argv[18];
                size_t a = 0;
                argv[a++] = tool;
                argv[a++] = strip_vmesh;
                argv[a++] = cfg->raw_path;
                argv[a++] = bake_dir;
                argv[a++] = "--id"; argv[a++] = id;
                argv[a++] = "--chunk"; argv[a++] = chunk_s;
                argv[a++] = cfg->allow_missing_chunks
                          ? "--allow-incomplete-raw"
                          : "--require-complete-raw";
                argv[a++] = "--window"; argv[a++] = lo_s; argv[a++] = hi_s;
                argv[a++] = "--diag-dark"; argv[a++] = dark_s;
                argv[a++] = "--raster-auto";
                argv[a] = NULL;
                if (qp_spawn("obj_bake_raw", argv) != 0) return -1;
            }
        } else {
            fprintf(stderr, "[sheet-%s] resume: %s already baked\n", tag, id);
        }
        if (TiffIO_load(arena, rawtex_tif, &tex, &D, &H, &W) != 0 || D != 1) {
            fprintf(stderr, "[sheet-%s] cannot load %s\n", tag, rawtex_tif);
            return -1;
        }
        layers[n_layers].tex = tex;
        layers[n_layers].w = (size_t)W;
        layers[n_layers].h = (size_t)H;
        if (TiffIO_load(arena, cover_tif, &cov, &D, &H, &W) != 0 || D != 1 ||
            (size_t)W != layers[n_layers].w ||
            (size_t)H != layers[n_layers].h) {
            fprintf(stderr, "[sheet-%s] coverage mismatch %s\n", tag,
                    cover_tif);
            return -1;
        }
        layers[n_layers].cov = cov;
        {
            /* optional stretch/squash diags for the quality-masked variant */
            char diag_tif[QP_MAX_PATH];
            uint8_t *dbuf = NULL;
            snprintf(diag_tif, sizeof diag_tif, "%s/%s_diagstretch.tif",
                     bake_dir, id);
            if (qp_exists(diag_tif) &&
                TiffIO_load(arena, diag_tif, &dbuf, &D, &H, &W) == 0 &&
                D == 1 && (size_t)W == layers[n_layers].w &&
                (size_t)H == layers[n_layers].h)
                stretch_buf[n_layers] = dbuf;
            dbuf = NULL;
            snprintf(diag_tif, sizeof diag_tif, "%s/%s_diagsquash.tif",
                     bake_dir, id);
            if (qp_exists(diag_tif) &&
                TiffIO_load(arena, diag_tif, &dbuf, &D, &H, &W) == 0 &&
                D == 1 && (size_t)W == layers[n_layers].w &&
                (size_t)H == layers[n_layers].h)
                squash_buf[n_layers] = dbuf;
        }
        if (qp_strip_min_u(strip_vmesh, &min_u[n_layers]) != 0) {
            fprintf(stderr, "[sheet-%s] no UV in %s\n", tag, strip_vmesh);
            return -1;
        }
        n_layers++;
    }
    if (n_layers == 0) {
        fprintf(stderr, "[sheet-%s] no strips produced\n", tag);
        return -1;
    }
    for (size_t s = 0; s < n_layers; s++) {
        long off = (long)floor(min_u[s]) - (long)floor(min_u[0]);
        layers[s].off_u = off > 0 ? off : 0;
        layers[s].off_v = 0;
    }

    snprintf(tex_png, sizeof tex_png, "%s/big_sheet_%s.png", out_root, tag);
    snprintf(prov_png, sizeof prov_png, "%s/big_sheet_%s_provenance.png",
             out_root, tag);
    if (SheetComposite_run(layers, n_layers, tex_png,
                           prov_png, out_tex, &stats) != 0) {
        fprintf(stderr, "[sheet-%s] composite failed\n", tag);
        return -1;
    }
    for (size_t s = 0; s < n_layers; s++) {
        snprintf(layer_png[s], sizeof layer_png[s],
                 "%s/big_sheet_%s_provenance_%zu.png", out_root, tag, s);
        layer_paths[s] = layer_png[s];
    }
    if (SheetComposite_write_layer_views(
            layers, n_layers, layer_paths) != 0) {
        fprintf(stderr, "[sheet-%s] individual provenance views failed\n",
                tag);
        return -1;
    }
    if (strcmp(tag, "post") == 0) {
        char legacy_prov[QP_MAX_PATH];
        snprintf(legacy_prov, sizeof legacy_prov,
                 "%s/big_sheet_provenance.png", out_root);
        if (qp_copy_file(prov_png, legacy_prov) != 0) {
            fprintf(stderr, "[sheet-%s] legacy provenance alias failed\n",
                    tag);
            return -1;
        }
    }
    if (out_w != NULL) *out_w = stats.width;
    if (out_h != NULL) *out_h = stats.height;
    fprintf(stderr,
            "[sheet-%s] %zux%zu painted=%zu crack_fill=%zu empty=%zu "
            "layers=%zu (direct %zu/%zu/%zu)\n",
            tag, stats.width, stats.height, stats.painted_px,
            stats.composite_crack_fill_px, stats.empty_px, n_layers,
            stats.direct_px[0], stats.direct_px[1], stats.direct_px[2]);

    {
        /* Shipped-view doctrine (wrap-separation campaign): the deliverable
         * sheet is quality-masked -- texels stretched or squashed past 4x
         * (dsmax >= 213 / dsmin <= 43 in the 128+42.5*log2(sigma) encoding)
         * show as VOID, not smear.  The mask only removes coverage, so
         * first-cover falls through to lower layers; the unmasked sheet
         * above remains the A/B primary. */
        enum { QP_MASK_STRETCH_HI = 213, QP_MASK_SQUASH_LO = 43 };
        SheetCompositeLayer masked[SHEET_COMPOSITE_MAX_LAYERS];
        SheetCompositeStats mstats;
        char mask_png[QP_MAX_PATH];
        size_t masked_out = 0;
        memcpy(masked, layers, sizeof masked);
        for (size_t s = 0; s < n_layers; s++) {
            size_t npx = layers[s].w * layers[s].h;
            uint8_t *mc = (uint8_t *)ARENA_ALLOC(arena, npx);
            memcpy(mc, layers[s].cov, npx);
            if (stretch_buf[s] != NULL || squash_buf[s] != NULL) {
                for (size_t i = 0; i < npx; i++) {
                    int bad =
                        (stretch_buf[s] != NULL &&
                         stretch_buf[s][i] >= QP_MASK_STRETCH_HI) ||
                        (squash_buf[s] != NULL && squash_buf[s][i] != 0 &&
                         squash_buf[s][i] <= QP_MASK_SQUASH_LO);
                    if (bad && mc[i] != 0) {
                        mc[i] = 0;
                        masked_out++;
                    }
                }
            }
            masked[s].cov = mc;
        }
        snprintf(mask_png, sizeof mask_png, "%s/big_sheet_%s_masked.png",
                 out_root, tag);
        if (SheetComposite_run(masked, n_layers, mask_png, NULL, NULL,
                               &mstats) == 0)
            fprintf(stderr,
                    "[sheet-%s] quality mask: %zu px masked -> painted=%zu "
                    "empty=%zu (%s)\n",
                    tag, masked_out, mstats.painted_px, mstats.empty_px,
                    mask_png);
        else
            fprintf(stderr,
                    "[sheet-%s] masked composite failed (non-fatal)\n", tag);
    }
    return 0;
}

/* ---- main ---------------------------------------------------------------- */

static void qp_capture_exe_dir(const char *argv0)
{
    const char *a = strrchr(argv0, '/');
    const char *b = strrchr(argv0, '\\');
    const char *slash = a != NULL && (b == NULL || a > b) ? a : b;
    if (slash != NULL) {
        size_t n = (size_t)(slash - argv0);
        if (n >= sizeof qp_exe_dir) n = sizeof qp_exe_dir - 1;
        memcpy(qp_exe_dir, argv0, n);
        qp_exe_dir[n] = '\0';
    }
}

/* Stage 4: load the parameterized quadribbon, run the collision shell
 * (turn-order repair + elastic shell + settle; radial motion only), and
 * write the untangled geometry with the SAME topology and UV -- stage 5
 * re-solves honest arc length on the moved geometry. */
static int qp_stage_untangle(const QpConfig *cfg, const char *in_vmesh,
                             const char *out_vmesh)
{
    Arena_T arena = Arena_new();
    MeshBinData m;
    QuadribbonUntangleOpts o;
    QuadribbonUntangleStats st;
    int rc = -1;
    if (arena == NULL) return -1;
    memset(&m, 0, sizeof m);
    if (MeshBin_read_arena(arena, in_vmesh, &m) != 0 || m.uv == NULL) {
        fprintf(stderr, "[stage 4 untangle] cannot read %s (with UV)\n",
                in_vmesh);
        Arena_dispose(&arena);
        return -1;
    }
    QuadribbonUntangle_defaults(&o);
    o.wrap_pitch_hint = cfg->wrap_spacing;
    {
        float *before = (float *)ARENA_ALLOC(arena,
                                             m.nv * 3 * sizeof(float));
        memcpy(before, m.verts, m.nv * 3 * sizeof(float));
        rc = QuadribbonUntangle_run(m.verts, m.nv, m.faces, m.nf, m.uv, NULL,
                                    cfg->axis_y, cfg->axis_x, &o, &st);
        if (rc < 0) {
            fprintf(stderr, "[stage 4 untangle] FAILED\n");
            Arena_dispose(&arena);
            return -1;
        }
        fprintf(stderr,
                "[stage 4 untangle] long %zu->%zu complete=%d accepted=%d "
                "settle=%d/%d movement rms/max=%.3f/%.3f\n",
                st.input_long_conflicts, st.output_long_conflicts,
                st.complete, st.accepted_rounds, st.settle_accepted,
                st.settle_rounds_run, st.movement_rms, st.movement_max);
        {
            char png[QP_MAX_PATH];
            size_t len = strlen(out_vmesh);
            snprintf(png, sizeof png, "%.*s_dr.png",
                     (int)(len > 6 ? len - 6 : len), out_vmesh);
            qp_dump_dr_png(before, m.verts, m.uv, m.nv, cfg->axis_y,
                           cfg->axis_x, png);
        }
    }
    rc = MeshBin_write(out_vmesh, m.verts, m.nv, m.faces, m.nf, m.uv);
    Arena_dispose(&arena);
    return rc;
}

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    QpConfig cfg;
    char soup[QP_MAX_PATH], logs[QP_MAX_PATH];
    char fit_dir[QP_MAX_PATH], fit_vmesh[QP_MAX_PATH];
    char labels[QP_MAX_PATH];
    char reopt1_dir[QP_MAX_PATH], reopt1_vmesh[QP_MAX_PATH];
    char untangle_dir[QP_MAX_PATH], untangled_vmesh[QP_MAX_PATH];
    char reopt2_dir[QP_MAX_PATH], reopt2_vmesh[QP_MAX_PATH];
    char tool[QP_MAX_PATH];
    char axis_y_s[32], axis_x_s[32], pitch_s[32];
    const char *mesh_dir = NULL, *out_dir = NULL, *config_path = NULL;
    int stop_stage = 0;
    double t0 = 0.0;

    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        int fails = 0;
        if (SheetComposite_selftest() != 0) fails++;
        if (QuadribbonUntangle_selftest() != 0) fails++;
        fprintf(stderr, "quadribbon --selftest: %s\n",
                fails == 0 ? "PASS" : "FAIL");
        return fails == 0 ? 0 : 1;
    }
    {
        const char *stop_name = NULL;
        int bad = 0;
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--stop-after") == 0 && i + 1 < argc) {
                stop_name = argv[++i];
            } else if (mesh_dir == NULL) {
                mesh_dir = argv[i];
            } else if (out_dir == NULL) {
                out_dir = argv[i];
            } else if (config_path == NULL) {
                config_path = argv[i];
            } else {
                bad = 1;
            }
        }
        stop_stage = qp_stage_index(stop_name);
        if (bad || mesh_dir == NULL || out_dir == NULL || stop_stage < 0) {
            fprintf(stderr,
                    "usage: quadribbon <mesh_dir> <out_dir> [config.json] "
                    "[--stop-after <stage>]\n"
                    "       quadribbon --selftest\n"
                    "  <mesh_dir>: pre-weld per-cube pile (recursive "
                    "*_final_all.vmesh)\n"
                    "              and/or plain world-frame *.vmesh files\n"
                    "  <stage>: concat|fit|reopt|untangle|reopt2|sheet "
                    "(or 1..6) -- finish that\n"
                    "           stage (debug bake included), then exit; rerun "
                    "resumes from it\n"
                    "  policy: %s\n", QP_DEFAULT_CONFIG);
            return 2;
        }
        if (config_path == NULL) config_path = QP_DEFAULT_CONFIG;
    }
    qp_capture_exe_dir(argv[0]);
    t0 = ves_clock_sec();
    qp_log_t0 = t0;

    arena = Arena_new();
    if (arena == NULL) return 1;
    if (qp_load_config(arena, config_path, &cfg) != 0) return 1;
    snprintf(axis_y_s, sizeof axis_y_s, "%.3f", cfg.axis_y);
    snprintf(axis_x_s, sizeof axis_x_s, "%.3f", cfg.axis_x);
    snprintf(pitch_s, sizeof pitch_s, "%.3f", cfg.wrap_spacing);

    qp_mkdir(out_dir);
    qp_join(logs, sizeof logs, out_dir, "logs");
    qp_mkdir(logs);
    snprintf(qp_log_dir, sizeof qp_log_dir, "%s", logs);
    qp_join(qp_pipeline_log, sizeof qp_pipeline_log, logs, "pipeline.log");
    remove(qp_pipeline_log);   /* fresh timeline per run */
    {
        char stop_note[64] = "";
        if (stop_stage > 0)
            snprintf(stop_note, sizeof stop_note, " stop-after=stage-%d",
                     stop_stage);
        qp_logf("quadribbon: mesh_dir=%s out_dir=%s config=%s\n"
                "  axis=(%.1f,%.1f) pitch=%.3f raw=%s window=[%ld,%ld] "
                "dark=%ld%s\n",
                mesh_dir, out_dir, config_path, cfg.axis_y, cfg.axis_x,
                cfg.wrap_spacing, cfg.raw_path, cfg.window_lo, cfg.window_hi,
                cfg.dark, stop_note);
    }

    /* stage 1: concat */
    qp_join(soup, sizeof soup, out_dir, "stage1_soup.vmesh");
    {
        char table[QP_MAX_PATH];
        double ts = ves_clock_sec();
        Arena_T soup_arena = Arena_new();
        int rc = 0;
        qp_join(table, sizeof table, out_dir, "stage1_soup_cubes.tsv");
        rc = soup_arena == NULL ? -1
           : qp_stage_concat(soup_arena, mesh_dir, soup, table);
        if (soup_arena != NULL) Arena_dispose(&soup_arena);
        if (rc != 0) return 1;
        qp_logf("[stage 1 concat] done (%.1fs)\n", ves_clock_sec() - ts);
    }
    if (stop_stage == 1) {
        qp_logf("quadribbon: STOP after stage 1 concat as requested "
                "(%.1fs total)\n", ves_clock_sec() - t0);
        Arena_dispose(&arena);
        return 0;
    }

    /* stage 2: fit */
    qp_join(fit_dir, sizeof fit_dir, out_dir, "stage2_fit");
    qp_join(fit_vmesh, sizeof fit_vmesh, fit_dir, "fit_ribbon.vmesh");
    qp_join(labels, sizeof labels, fit_dir,
            "fit_ribbon_reconstruction_component.i32");
    if (!MeshBin_looks_complete(fit_vmesh)) {
        qp_mkdir(fit_dir);
        qp_tool_path(tool, sizeof tool, "scroll_ribbon");
        {
            const char *argv2[14];
            size_t a = 0;
            argv2[a++] = tool;
            argv2[a++] = soup;
            argv2[a++] = fit_dir;
            argv2[a++] = "--id"; argv2[a++] = "fit";
            argv2[a++] = "--component-global";
            argv2[a++] = "--axis-point";
            argv2[a++] = "0"; argv2[a++] = axis_y_s; argv2[a++] = axis_x_s;
            argv2[a++] = "--wrap-spacing"; argv2[a++] = pitch_s;
            argv2[a++] = "--scaffold-solve";
            argv2[a] = NULL;
            if (qp_spawn("scroll_ribbon(fit)", argv2) != 0) return 1;
        }
        /* the metric stages discover sidecars by FIXED names beside their
         * input; publish plain-name copies once, at the fit output */
        {
            static const char *pairs[][2] = {
                { "fit_ribbon_phase.f32", "ribbon_phase.f32" },
                { "fit_ribbon_material_identity.i32",
                  "ribbon_material_identity.i32" },
            };
            for (size_t i = 0; i < 2; i++) {
                char src[QP_MAX_PATH], dst[QP_MAX_PATH];
                qp_join(src, sizeof src, fit_dir, pairs[i][0]);
                qp_join(dst, sizeof dst, fit_dir, pairs[i][1]);
                if (qp_exists(src) && qp_copy_file(src, dst) != 0) {
                    fprintf(stderr, "quadribbon: cannot publish %s\n", dst);
                    return 1;
                }
            }
        }
    } else {
        qp_logf("[stage 2 fit] resume: already complete\n");
    }
    if (qp_debug_bake(&cfg, fit_vmesh, out_dir, "stage2_fit") != 0)
        return 1;
    qp_logf("[stage 2 fit] done\n");
    if (stop_stage == 2) {
        qp_logf("quadribbon: STOP after stage 2 fit as requested "
                "(%.1fs total)\n", ves_clock_sec() - t0);
        Arena_dispose(&arena);
        return 0;
    }

    /* stage 3: reoptimize UV (metric projection + consistency solve) */
    qp_join(reopt1_dir, sizeof reopt1_dir, out_dir, "stage3_reopt");
    qp_join(reopt1_vmesh, sizeof reopt1_vmesh, reopt1_dir,
            "reopt_ribbon.vmesh");
    if (!MeshBin_looks_complete(reopt1_vmesh)) {
        qp_mkdir(reopt1_dir);
        qp_tool_path(tool, sizeof tool, "scroll_ribbon");
        {
            const char *argv2[14];
            size_t a = 0;
            argv2[a++] = tool;
            argv2[a++] = fit_vmesh;
            argv2[a++] = reopt1_dir;
            argv2[a++] = "--id"; argv2[a++] = "reopt";
            argv2[a++] = "--umb-y"; argv2[a++] = axis_y_s;
            argv2[a++] = "--umb-x"; argv2[a++] = axis_x_s;
            argv2[a++] = "--wrap-spacing"; argv2[a++] = pitch_s;
            argv2[a++] = "--preserve-input-topology";
            argv2[a] = NULL;
            if (qp_spawn("scroll_ribbon(reopt)", argv2) != 0) return 1;
        }
    } else {
        qp_logf("[stage 3 reopt] resume: already complete\n");
    }
    if (qp_debug_bake(&cfg, reopt1_vmesh, out_dir, "stage3_reopt") != 0)
        return 1;
    qp_logf("[stage 3 reopt] done\n");
    if (stop_stage == 3) {
        qp_logf("quadribbon: STOP after stage 3 reopt as requested "
                "(%.1fs total)\n", ves_clock_sec() - t0);
        Arena_dispose(&arena);
        return 0;
    }

    /* stage 4: untangle */
    qp_join(untangle_dir, sizeof untangle_dir, out_dir, "stage4_untangle");
    qp_join(untangled_vmesh, sizeof untangled_vmesh, untangle_dir,
            "untangled_ribbon.vmesh");
    if (!MeshBin_looks_complete(untangled_vmesh)) {
        qp_mkdir(untangle_dir);
        if (qp_stage_untangle(&cfg, reopt1_vmesh, untangled_vmesh) != 0)
            return 1;
        /* sidecars travel beside the untangled mesh for the next reopt */
        {
            static const char *names[] = {
                "ribbon_phase.f32", "ribbon_material_identity.i32",
            };
            for (size_t i = 0; i < 2; i++) {
                char src[QP_MAX_PATH], dst[QP_MAX_PATH];
                qp_join(src, sizeof src, fit_dir, names[i]);
                qp_join(dst, sizeof dst, untangle_dir, names[i]);
                if (qp_exists(src) && qp_copy_file(src, dst) != 0) return 1;
            }
        }
    } else {
        qp_logf("[stage 4 untangle] resume: already complete\n");
    }
    if (qp_debug_bake(&cfg, untangled_vmesh, out_dir,
                      "stage4_untangle") != 0)
        return 1;
    qp_logf("[stage 4 untangle] done\n");
    if (stop_stage == 4) {
        qp_logf("quadribbon: STOP after stage 4 untangle as requested "
                "(%.1fs total)\n", ves_clock_sec() - t0);
        Arena_dispose(&arena);
        return 0;
    }

    /* stage 5: reoptimize again on the untangled geometry */
    qp_join(reopt2_dir, sizeof reopt2_dir, out_dir, "stage5_reopt");
    qp_join(reopt2_vmesh, sizeof reopt2_vmesh, reopt2_dir,
            "reopt_ribbon.vmesh");
    if (!MeshBin_looks_complete(reopt2_vmesh)) {
        qp_mkdir(reopt2_dir);
        qp_tool_path(tool, sizeof tool, "scroll_ribbon");
        {
            const char *argv2[14];
            size_t a = 0;
            argv2[a++] = tool;
            argv2[a++] = untangled_vmesh;
            argv2[a++] = reopt2_dir;
            argv2[a++] = "--id"; argv2[a++] = "reopt";
            argv2[a++] = "--umb-y"; argv2[a++] = axis_y_s;
            argv2[a++] = "--umb-x"; argv2[a++] = axis_x_s;
            argv2[a++] = "--wrap-spacing"; argv2[a++] = pitch_s;
            argv2[a++] = "--preserve-input-topology";
            argv2[a] = NULL;
            if (qp_spawn("scroll_ribbon(reopt2)", argv2) != 0) return 1;
        }
    } else {
        qp_logf("[stage 5 reopt] resume: already complete\n");
    }
    if (qp_debug_bake(&cfg, reopt2_vmesh, out_dir, "stage5_reopt") != 0)
        return 1;
    qp_logf("[stage 5 reopt] done\n");
    if (stop_stage == 5) {
        qp_logf("quadribbon: STOP after stage 5 reopt2 as requested "
                "(%.1fs total)\n", ves_clock_sec() - t0);
        Arena_dispose(&arena);
        return 0;
    }

    /* stage 6: sheets -- "pre" (before untangle) and "post" (final) */
    {
        char sheet_pre[QP_MAX_PATH], sheet_post[QP_MAX_PATH];
        uint8_t *tex_pre = NULL, *tex_post = NULL;
        size_t wpre = 0, hpre = 0, wpost = 0, hpost = 0;
        qp_join(sheet_pre, sizeof sheet_pre, out_dir, "stage6_sheet_pre");
        qp_join(sheet_post, sizeof sheet_post, out_dir, "stage6_sheet_post");
        if (qp_stage_sheet(arena, &cfg, reopt1_vmesh, labels, sheet_pre,
                           "pre", out_dir, &tex_pre, &wpre, &hpre) != 0) {
            qp_logf("[stage 6 sheet-pre] FAILED (see logs above)\n");
            return 1;
        }
        if (qp_stage_sheet(arena, &cfg, reopt2_vmesh, labels, sheet_post,
                           "post", out_dir, &tex_post, &wpost, &hpost) != 0) {
            qp_logf("[stage 6 sheet-post] FAILED (see logs above)\n");
            return 1;
        }
        if (tex_pre != NULL && tex_post != NULL) {
            size_t W = wpre > wpost ? wpre : wpost;
            size_t H = hpre + 8 + hpost;
            uint8_t *cmp = (uint8_t *)calloc(W * H, 3);
            if (cmp != NULL) {
                char cmp_png[QP_MAX_PATH];
                memset(cmp, 24, W * H * 3);
                for (size_t y = 0; y < hpre; y++)
                    memcpy(cmp + y * W * 3, tex_pre + y * wpre * 3,
                           wpre * 3);
                for (size_t y = 0; y < hpost; y++)
                    memcpy(cmp + (hpre + 8 + y) * W * 3,
                           tex_post + y * wpost * 3, wpost * 3);
                qp_join(cmp_png, sizeof cmp_png, out_dir,
                        "big_sheet_pre_vs_post.png");
                VesPng_write_rgb(cmp_png, cmp, (int)W, (int)H);
                free(cmp);
            }
        }
        free(tex_pre);
        free(tex_post);
    }

    qp_logf("quadribbon: OK (total %.1fs)\n", ves_clock_sec() - t0);
    Arena_dispose(&arena);
    return 0;
}
