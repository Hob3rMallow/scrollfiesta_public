/* ============================================================================
 * quadribbon_pipeline.c -- ONE command: meshed scroll region in, big sheet out.
 *
 *   quadribbon <mesh_dir|welded.vmesh> <out_dir> [config.json]
 *              [--trust-gauge] [--stop-after <stage>]
 *   quadribbon --selftest
 *
 * The canonical unwrap/flatten order of operations, in one tool (2026-09-01):
 *   stage 1  mesh     -- the meshed input: a welded world-frame VMESH (the
 *            measured 4x5x5 champion path) or a pre-weld per-cube pile
 *            (recursive *_final_all.vmesh, deduplicated by cube id, concatenated
 *            without welding; the certificate-limited 21x route). Before
 *            unwrapping, remove small disconnected components by physical area
 *            on each original mesh, then straighten/route the surviving mesh.
 *   stage 2  unwrap   -- scroll_ribbon --winding-only: initial winding field,
 *            winding registration + MRF, gauge sync and the arc-length U map,
 *            written as the winding certificate VMESH (topology preserved).
 *   stage 3  fit      -- scroll_ribbon --scaffold-solve on the certificate
 *            (plus --trust-gauge for the streaming carried-U mode):
 *            claims-mode lattice ribbon (peel bands, wrap-safe emission), the
 *            observations.  ribbon_verdict scores it; RAW bake.
 *   stage 4  solid    -- solidify EVERY collision chart emitted by stage 3 as
 *            an independent member of an authoritative atlas
 *            (src/flatten/quadribbon_solidify.c).  Charts retain the same
 *            absolute U/V coordinates and are never merged or arbitrated.
 *   stage 5  optimize -- rounds of untangle (radial collision shell, in-process)
 *            followed by a REFIT (the same solidify on the moved observations),
 *            stopping at a conflict-free state, the round cap, or no further
 *            gain; optional final metric re-solve (QS_REPARAM_FINAL).
 *   stage 6  sheet    -- layer strips (spawned) -> RAW bakes (spawned) ->
 *            first-cover composite (in-process).  big_sheet_{pre,post,
 *            provenance,pre_vs_post}.png land at the out_dir root; "pre" is
 *            the solid ribbon, "post" the optimized one.
 *
 * --stop-after mesh|unwrap|fit|solid|optimize|sheet (or 1..6) finishes the
 * named stage (verdict + debug bake included) and exits cleanly.  Re-running
 * without the flag resumes from the completed artifacts.
 *
 * All policy comes from the production config JSON + pipeline_constants.h;
 * there are no tuning flags.  Stages resume: a stage whose primary artifact
 * already exists is skipped, so a crashed run re-enters where it stopped.
 * ==========================================================================*/
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/union_find.h"
#include "../whole/axis_warp.h"
#include "../common/json_read.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_kibble.h"
#include "../common/pipeline_constants.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include "../flatten/quad_strip.h"
#include "../flatten/quadribbon_solidify.h"
#include "../flatten/ridge_track.h"
#include "../common/zarr_u8.h"
#include "../flatten/quadribbon_untangle.h"
#include "../flatten/sheet_composite.h"
#include "../flatten/sparse_solve.h"
#include "../flatten/scroll_coordinate.h"
#include "../flatten/scroll_source.h"
#include "../flatten/ribbon_lod.h"
#include "../flatten/quad_field.h"
#include "../flatten/quad_ribbon_fit.h"
#include "../flatten/scroll_model.h"
#include "../flatten/material_evidence.h"
#include "../flatten/material_front.h"
#include "../flatten/material_observation.h"
#include "../flatten/winding_register.h"

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

enum { QP_MAX_PATH = 2048, QP_MAX_PILE = 65536 };

static const char QP_DEFAULT_CONFIG[] =
    "configs/default.json";

typedef struct {
    double axis_y, axis_x;
    double wrap_spacing;
    char raw_path[QP_MAX_PATH];
    double bake_normal_reach;  /* bake.normal_reach_vox: max-sample +/- reach along the normal */
    long raw_chunk;
    long window_lo, window_hi, dark;
    int allow_missing_chunks;
    long threads;              /* compute.threads for the scroll_ribbon spawns */
    /* geometry.axis_table: sampled z,y,x umbilicus curve.  The scroll axis
     * DRIFTS (PHerc0139: 40 vox off the manifest umbilicus at z 4384, 600 by
     * z 6900); a fixed axis mis-unwraps everything above the core.  When set,
     * stage 1 STRAIGHTENS the mesh (each z-slice translated so the curve maps
     * onto the config umbilicus, mesh_world.vmesh kept), every later stage
     * runs on straightened coordinates, and un-straightened `*_world.vmesh`
     * copies feed the bakes and sheets (which sample the CT). */
    char axis_table_path[QP_MAX_PATH];
    AxisWarp axis_warp;
    int axis_warp_armed;
    long canonical_anchor[3]; /* geometry.canonical_block_anchor_zyx */
    int canonical_anchor_armed;
    int winding_sense;         /* geometry.winding_sense: +1/-1 pin, 0 auto */
    int overlap_family;        /* registration.overlap_family: 1 = arm the exact
                                * shared-cube OVERLAP relation family on piles
                                * (default), 0 = A/B switch: leave it unarmed */
    double core_wall_radius;   /* core.wall_radius_voxels: components lying
                                * entirely inside the innermost wall (r_max
                                * below this) are the crumpled roll end and
                                * are routed to extras before the certificate.
                                * 0 = off (weld/champion lanes byte-identical) */
    double kibble_min_area;    /* mesh_cleanup.min_component_area_vox2; 0 disables */
    char assembly_ribbon[QP_MAX_PATH]; /* assembly.ribbon: an EXTERNAL lattice ribbon (rows = z-slices,
                                * _support/_phase sidecars) replaces stages 1-3; the chart assembler's
                                * assembled_ribbon.vmesh is the measured case */
    char assembly_source[QP_MAX_PATH]; /* assembly.source: the mesh the verdict's coverage counts */
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
        "mesh", "unwrap", "fit", "solid", "optimize", "sheet"
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

static int qp_file_has_text(const char *path, const char *needle)
{
    FILE *f = fopen(path, "rb");
    char line[4096];
    if (f == NULL || needle == NULL) {
        if (f != NULL) fclose(f);
        return 0;
    }
    while (fgets(line, sizeof line, f) != NULL) {
        if (strstr(line, needle) != NULL) {
            fclose(f);
            return 1;
        }
    }
    fclose(f);
    return 0;
}

static uint64_t qp_file_stamp(const char *path)
{
#ifdef _WIN32
    WIN32_FILE_ATTRIBUTE_DATA data;
    ULARGE_INTEGER t;
    if (!GetFileAttributesExA(path, GetFileExInfoStandard, &data)) return 0;
    t.LowPart = data.ftLastWriteTime.dwLowDateTime;
    t.HighPart = data.ftLastWriteTime.dwHighDateTime;
    return (uint64_t)t.QuadPart;
#else
    struct stat st;
    if (stat(path, &st) != 0) return 0;
    return (uint64_t)st.st_mtime;
#endif
}

static int qp_source_newer(const char *source, const char *derived)
{
    uint64_t s = qp_file_stamp(source), d = qp_file_stamp(derived);
    return s != 0 && (d == 0 || s > d);
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

static int qp_replace_file(const char *source, const char *destination)
{
#ifdef _WIN32
    return MoveFileExA(source,destination,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH) ? 0 : -1;
#else
    return rename(source,destination)==0 ? 0 : -1;
#endif
}

static int qp_vector3(const JsonValue *value, double out[3])
{
    if (!value || Json_array_len(value)!=3) return -1;
    for (int a=0;a<3;a++) {
        out[a]=Json_as_double(Json_array_get(value,a),NAN);
        if (!isfinite(out[a])) return -1;
    }
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
        const JsonValue *cleanup = Json_object_get(root, "mesh_cleanup");
        cfg->kibble_min_area = Json_member_double(cleanup,
                                      "min_component_area_vox2", 1024.0);
        if (!isfinite(cfg->kibble_min_area) || cfg->kibble_min_area < 0) {
            fprintf(stderr, "quadribbon: mesh_cleanup.min_component_area_vox2 must be finite and >= 0\n");
            return -1;
        }
    }
    AxisWarp_init(&cfg->axis_warp);
    cfg->axis_warp_armed = 0;
    cfg->canonical_anchor_armed = 0;
    {
        const JsonValue *anchor =
            Json_object_get(geom, "canonical_block_anchor_zyx");
        if (anchor != NULL) {
            long az = Json_as_long(Json_array_get(anchor, 0), LONG_MIN);
            long ay = Json_as_long(Json_array_get(anchor, 1), LONG_MIN);
            long ax = Json_as_long(Json_array_get(anchor, 2), LONG_MIN);
            if (az == LONG_MIN || ay == LONG_MIN || ax == LONG_MIN ||
                Json_array_get(anchor, 3) != NULL) {
                fprintf(stderr, "quadribbon: geometry.canonical_block_anchor_zyx "
                        "must contain exactly three integer cube origins\n");
                return -1;
            }
            cfg->canonical_anchor[0] = az;
            cfg->canonical_anchor[1] = ay;
            cfg->canonical_anchor[2] = ax;
            cfg->canonical_anchor_armed = 1;
        }
    }
    cfg->winding_sense = (int)Json_member_long(geom, "winding_sense", 0);
    cfg->overlap_family = 0;   /* measured default OFF (2026-09-02 A/B) */
    {
        JsonValue *reg = Json_object_get(root, "registration");
        if (reg != NULL)
            cfg->overlap_family = Json_member_long(reg, "overlap_family", 0) != 0;
    }
    if (cfg->winding_sense > 0) cfg->winding_sense = 1;
    if (cfg->winding_sense < 0) cfg->winding_sense = -1;
    {
        const char *ap = Json_as_string(Json_object_get(geom, "axis_table"));
        const JsonValue *shaft=Json_object_get(geom,"shaft_coordinates");
        if (shaft && ap && ap[0]) {
            fprintf(stderr,"quadribbon: shaft_coordinates and axis_table are mutually exclusive\n");
            return -1;
        }
        if (shaft) {
            const JsonValue *rows=Json_object_get(shaft,"points_um_zyx");
            const JsonValue *origin=Json_object_get(shaft,"source_origin_um_zyx");
            const JsonValue *normal=Json_object_get(shaft,"initial_normal_zyx");
            AxisWarpPhysical physical={0};
            size_t count=rows ? Json_array_len(rows) : 0;
            if (count<3 || count>1048576 ||
                qp_vector3(Json_object_get(shaft,"source_voxel_um_zyx"),physical.source_voxel_um_zyx)!=0 ||
                (origin && qp_vector3(origin,physical.source_origin_um_zyx)!=0) ||
                (normal && qp_vector3(normal,physical.initial_normal_zyx)!=0)) {
                fprintf(stderr,"quadribbon: invalid physical shaft points or source coordinate units\n"); return -1;
            }
            physical.has_initial_normal=normal!=NULL;
            physical.metric_voxel_um=Json_member_double(shaft,"metric_voxel_um",NAN);
            physical.maximum_radius_um=Json_member_double(shaft,"maximum_radius_um",8500.);
            physical.minimum_jacobian=Json_member_double(shaft,"minimum_jacobian",.1);
            physical.ambiguity_um=Json_member_double(shaft,"ambiguity_um",1e-4);
            double *points=(double*)ARENA_ALLOC(arena,count*3*sizeof(double));
            for (size_t i=0;i<count;i++) if (qp_vector3(Json_array_get(rows,i),points+3*i)!=0) return -1;
            if (AxisWarp_create_physical(&cfg->axis_warp,points,count,&physical,cfg->axis_y,cfg->axis_x)!=0) {
                fprintf(stderr,"quadribbon: invalid or unresolved physical shaft geometry\n"); return -1;
            }
            cfg->axis_warp_armed=1;
            snprintf(cfg->axis_table_path,sizeof cfg->axis_table_path,"embedded physical shaft");
            if (!cfg->winding_sense) {
                fprintf(stderr,"quadribbon: shaft_coordinates requires geometry.winding_sense (+1/-1)\n"); return -1;
            }
        }
        if (ap != NULL && ap[0] != '\0') {
            snprintf(cfg->axis_table_path, sizeof cfg->axis_table_path, "%s", ap);
            if (AxisWarp_load_csv(&cfg->axis_warp, ap) != 0 ||
                !AxisWarp_valid(&cfg->axis_warp)) {
                fprintf(stderr, "quadribbon: cannot load geometry.axis_table %s\n", ap);
                return -1;
            }
            cfg->axis_warp.reference_y = cfg->axis_y;
            cfg->axis_warp.reference_x = cfg->axis_x;
            cfg->axis_warp_armed = 1;
            if (cfg->winding_sense == 0) {
                fprintf(stderr, "quadribbon: geometry.axis_table requires "
                        "geometry.winding_sense (+1/-1): the auto sense is a "
                        "coin flip once the axis moves\n");
                return -1;
            }
        }
    }
    {
        const char *rp = Json_as_string(Json_object_get(raw, "path"));
        if (rp == NULL) {
            fprintf(stderr, "quadribbon: config missing raw_source.path\n");
            return -1;
        }
        snprintf(cfg->raw_path, sizeof cfg->raw_path, "%s", rp);
    }
    cfg->raw_chunk = Json_member_long(raw, "chunk_size", 128);
    cfg->bake_normal_reach = Json_member_double(bake, "normal_reach_vox", 2.0);
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
    {
        const JsonValue *comp = Json_object_get(root, "compute");
        cfg->threads = comp != NULL ? Json_member_long(comp, "threads", 8) : 8;
        if (cfg->threads < 1) cfg->threads = 1;
    }
    {
        const JsonValue *core = Json_object_get(root, "core");
        cfg->core_wall_radius =
            core != NULL ? Json_member_double(core, "wall_radius_voxels", 0.0)
                         : 0.0;
        if (cfg->core_wall_radius < 0.0) cfg->core_wall_radius = 0.0;
    }
    {
        const JsonValue *asm_ = Json_object_get(root, "assembly");
        cfg->assembly_ribbon[0] = '\0';
        cfg->assembly_source[0] = '\0';
        if (asm_ != NULL) {
            const char *rb = Json_as_string(Json_object_get(asm_, "ribbon"));
            const char *sc = Json_as_string(Json_object_get(asm_, "source"));
            if (rb != NULL) snprintf(cfg->assembly_ribbon, sizeof cfg->assembly_ribbon, "%s", rb);
            if (sc != NULL) snprintf(cfg->assembly_source, sizeof cfg->assembly_source, "%s", sc);
            if (cfg->assembly_ribbon[0] != '\0' && !qp_exists(cfg->assembly_ribbon)) {
                fprintf(stderr, "quadribbon: assembly.ribbon not found: %s\n", cfg->assembly_ribbon);
                return -1;
            }
        }
    }
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

/* Checked source <-> metric conversion; topology and UV remain unchanged. */
static int qp_write_warped_copy(const QpConfig *cfg, const char *in_vmesh,
                                const char *out_vmesh, int direction)
{
    Arena_T arena = Arena_new();
    MeshBinData m;
    float *v2 = NULL;
    size_t i = 0;
    int rc = -1;
    char temporary[QP_MAX_PATH];
    if (arena == NULL) return -1;
    memset(&m, 0, sizeof m);
    if (MeshBin_read_arena(arena, in_vmesh, &m) != 0) { Arena_dispose(&arena); return -1; }
    v2 = (float *)ARENA_ALLOC(arena, (m.nv > 0 ? m.nv : 1) * 3 * sizeof(float));
    for (i = 0; i < m.nv; i++) {
        double in[3]={m.verts[3*i],m.verts[3*i+1],m.verts[3*i+2]},out[3];
        uint32_t flags=0;
        int status=direction>0 ? AxisWarp_to_metric(&cfg->axis_warp,in,out,&flags)
                               : AxisWarp_to_world(&cfg->axis_warp,in,out,&flags);
        if (status) {
            fprintf(stderr,"[axis] %s vertex %zu refused (status %d, flags %u); no converted mesh published\n",
                    in_vmesh,i,status,(unsigned)flags);
            Arena_dispose(&arena); return -1;
        }
        for (int a=0;a<3;a++) {
            if (!isfinite(out[a]) || fabs(out[a])>FLT_MAX) { Arena_dispose(&arena); return -1; }
            v2[3*i+a]=(float)out[a];
        }
    }
    if (snprintf(temporary,sizeof temporary,"%s.coordinate.tmp",out_vmesh)>=(int)sizeof temporary) {
        Arena_dispose(&arena); return -1;
    }
    rc = MeshBin_write(temporary, v2, m.nv, m.faces, m.nf, m.uv);
    if (rc==0) rc=qp_replace_file(temporary,out_vmesh);
    Arena_dispose(&arena);
    return rc;
}

/* Coordinate identity is checked before resume; paths and timestamps do not
 * identify a shaft, its transported frame, physical units or valid domain. */
#define QP_PHYSICAL_TILE_POLICY "swept_shaft_aabb_v1"
static int qp_coordinate_identity(const QpConfig *cfg, const char *path, const char *existing_mesh)
{
    char identity[24],temporary[QP_MAX_PATH];
    const AxisWarp *w=&cfg->axis_warp;
    snprintf(identity,sizeof identity,"%016llx",(unsigned long long)AxisWarp_fingerprint(w));
    if (qp_exists(path)) {
        Arena_T arena=Arena_new(); const char *error=NULL;
        const JsonValue *record=Json_parse_file(arena,path,&error);
        const char *schema=record ? Json_as_string(Json_object_get(record,"schema")) : NULL;
        const char *saved=record ? Json_as_string(Json_object_get(record,"fingerprint")) : NULL;
        int equal=schema && saved && !strcmp(schema,"vesuvius-coordinate-frame-v1") && !strcmp(saved,identity);
        if (w->physical) {
            const char *selection=record ? Json_as_string(Json_object_get(record,"physical_tile_selection")) : NULL;
            equal=equal && selection && !strcmp(selection,QP_PHYSICAL_TILE_POLICY);
        }
        Arena_dispose(&arena);
        if (!equal) {
            fprintf(stderr,"[axis] coordinate frame differs or is invalid; use a fresh output directory (%s)\n",path);
            return -1;
        }
        return 0;
    }
    if (cfg->axis_warp_armed && existing_mesh && qp_exists(existing_mesh)) {
        fprintf(stderr,"[axis] existing mesh has no coordinate identity; use a fresh output directory\n"); return -1;
    }
    if (snprintf(temporary,sizeof temporary,"%s.tmp",path)>=(int)sizeof temporary) return -1;
    FILE *f=fopen(temporary,"wb"); if (!f) return -1;
    fprintf(f,"{\n  \"schema\": \"vesuvius-coordinate-frame-v1\",\n  \"fingerprint\": \"%s\",\n"
              "  \"mode\": \"%s\",\n  \"reference_yx\": [%.17g, %.17g]",
              identity,w->physical ? "physical_arc" : cfg->axis_warp_armed ? "legacy_slice_translation" : "identity",
              w->reference_y,w->reference_x);
    if (w->physical) {
        const AxisWarpPhysical *c=&w->physical_config;
        fprintf(f,",\n  \"physical_tile_selection\": \"%s\"",QP_PHYSICAL_TILE_POLICY);
        fprintf(f,",\n  \"source_voxel_um_zyx\": [%.17g, %.17g, %.17g],\n"
                  "  \"source_origin_um_zyx\": [%.17g, %.17g, %.17g],\n"
                  "  \"metric_voxel_um\": %.17g,\n  \"maximum_radius_um\": %.17g,\n"
                  "  \"minimum_jacobian\": %.17g,\n  \"ambiguity_um\": %.17g,\n"
                  "  \"has_initial_normal\": %s,\n  \"initial_normal_zyx\": [%.17g, %.17g, %.17g]",
                  c->source_voxel_um_zyx[0],c->source_voxel_um_zyx[1],c->source_voxel_um_zyx[2],
                  c->source_origin_um_zyx[0],c->source_origin_um_zyx[1],c->source_origin_um_zyx[2],
                  c->metric_voxel_um,c->maximum_radius_um,c->minimum_jacobian,c->ambiguity_um,
                  c->has_initial_normal ? "true" : "false",c->initial_normal_zyx[0],c->initial_normal_zyx[1],c->initial_normal_zyx[2]);
    }
    fprintf(f,",\n  \"point_units\": \"%s\",\n  \"points_zyx\": [",w->physical ? "um" : "source_voxels");
    for (size_t i=0;i<w->n;i++) fprintf(f,"%s[%.17g, %.17g, %.17g]",i ? ",\n    " : "\n    ",w->z[i],w->y[i],w->x[i]);
    fprintf(f,"\n  ]\n}\n");
    int failed=ferror(f); if (fclose(f)!=0 || failed) return -1;
    return qp_replace_file(temporary,path);
}

/* World-coordinate twin of a straightened ribbon for CT-facing consumers
 * (bakes, sheets, cross-sections): <stem>_world.vmesh next to it.  Returns
 * the path to use; the input path when the warp is not armed. */
static const char *qp_world_twin(const QpConfig *cfg, const char *vmesh,
                                 char *buf, size_t cap)
{
    size_t n = strlen(vmesh);
    if (!cfg->axis_warp_armed) return vmesh;
    if (n < 6 || strcmp(vmesh + n - 6, ".vmesh") != 0) return vmesh;
    snprintf(buf, cap, "%.*s_world.vmesh", (int)(n - 6), vmesh);
    /* A complete twin can still be stale after its straightened parent was
     * rebuilt.  Regenerate exactly when its dependency is newer. */
    if ((!MeshBin_looks_complete(buf) || qp_source_newer(vmesh, buf)) &&
        qp_write_warped_copy(cfg, vmesh, buf, -1) != 0) {
        fprintf(stderr, "[axis] cannot write world twin %s; RAW consumer refused\n", buf);
        return NULL;
    }
    return buf;
}

static const AxisWarp *g_qp_axis_warp = NULL;   /* armed warp for --tile */

/* Optional block subgrid: inclusive source-voxel origin bounds (cube ids).
 * Set by --subgrid; when armed, the concat scan keeps only cubes whose origin
 * lies in the box.  This is the block-hierarchical primitive -- a block is a
 * subgrid of the full per-cube dump, certified and fit on its own, so the
 * whole-scroll monolith is never built.  Mirrors grid_weld --subgrid. */
static int  g_qp_subgrid_armed = 0;
static long g_qp_sg[6] = {0,0,0,0,0,0};  /* z0 z1 y0 y1 x0 x1 */
static int  g_qp_canonical_blocks = 0;

static int qp_in_subgrid(long oz, long oy, long ox)
{
    if (!g_qp_subgrid_armed) return 1;
    return oz >= g_qp_sg[0] && oz <= g_qp_sg[1] &&
           oy >= g_qp_sg[2] && oy <= g_qp_sg[3] &&
           ox >= g_qp_sg[4] && ox <= g_qp_sg[5];
}

/* Cylindrical TILE: the block-sizing law measured on 2026-08-30/09-01 is
 * radial (a tile must be <= ~10-12 wraps deep or it squashes, and a low
 * enough turn density that its lift stays one island), so the streaming
 * partition is an (r, theta, z) tile about the scroll axis, not an xyz box.
 * A tile keeps every WHOLE cube whose 128^3 footprint meets the annular
 * sector (sampled on a 5x5 grid of its xy footprint) with origin z in
 * [z0, z1]; adjacent tiles therefore SHARE their boundary cubes, and those
 * shared cubes are the exact registration evidence between tiles (same
 * dump file, same vertex order).  --tile r0 r1 theta0 theta1 z0 z1
 * (radii in vox, angles in degrees [0,360) about the config axis, z as cube
 * origins); combinable with --subgrid. Physical shafts instead use conservative
 * swept-shaft bounds over an arc interval, as described in qp_in_tile below. */
static int    g_qp_tile_armed = 0;
static double g_qp_tile[6] = {0,0,0,0,0,0};   /* r0 r1 th0 th1 z0 z1 */
static double g_qp_tile_axis[2] = {0.0, 0.0}; /* y, x (config) */
static int g_qp_tile_geometry_error=0;

static int qp_angle_in_sector(double deg, double th0, double th1)
{
    if (th0 <= th1) return deg >= th0 && deg <= th1;
    return deg >= th0 || deg <= th1;              /* wraps through 0 */
}

static int qp_in_tile(long oz, long oy, long ox)
{
    int i = 0, j = 0;
    if (!g_qp_tile_armed) return 1;
    double ay = g_qp_tile_axis[0], ax = g_qp_tile_axis[1];
    if (g_qp_axis_warp && g_qp_axis_warp->physical) {
        /* In physical-arc mode, the last interval is s in metric voxels.
         * Point probes can miss even a valid tube entirely inside the cube.
         * Bound the continuous shaft interval, expanded by the outer radius,
         * and retain intersecting WHOLE source cubes. This broad-phase bound
         * deliberately overincludes inner-radius and angular sectors; it is
         * not exact clipping, and vertex coordinate checks still apply. */
        const AxisWarp *w=g_qp_axis_warp;
        const AxisWarpPhysical *c=&w->physical_config;
        double begin=g_qp_tile[4]*c->metric_voxel_um,end=g_qp_tile[5]*c->metric_voxel_um;
        double radius=g_qp_tile[1]*c->metric_voxel_um,lower[3],upper[3];
        double origin[3]={(double)oz,(double)oy,(double)ox};
        for (int a=0;a<6;a++) if (!isfinite(g_qp_tile[a])) { g_qp_tile_geometry_error=1; return 0; }
        if (!isfinite(begin) || !isfinite(end) || !isfinite(radius)) { g_qp_tile_geometry_error=1; return 0; }
        if (radius<0. || g_qp_tile[0]>g_qp_tile[1] ||
            g_qp_tile[0]>c->maximum_radius_um/c->metric_voxel_um) return 0;
        begin=fmax(0.,begin); end=fmin(ShaftWarp_length_um(w->physical),end);
        if (begin>end) return 0;
        radius=fmin(radius,c->maximum_radius_um);
        if (ShaftWarp_bounds_arc(w->physical,begin,end,radius,lower,upper)) {
            g_qp_tile_geometry_error=1; return 0;
        }
        for (int a=0;a<3;a++) {
            double cube_lo=origin[a]*c->source_voxel_um_zyx[a]+c->source_origin_um_zyx[a];
            double cube_hi=(origin[a]+128.)*c->source_voxel_um_zyx[a]+c->source_origin_um_zyx[a];
            if (!isfinite(cube_lo) || !isfinite(cube_hi)) { g_qp_tile_geometry_error=1; return 0; }
            if (cube_lo>upper[a] || cube_hi<lower[a]) return 0;
        }
        return 1;
    }
    if ((double)oz < g_qp_tile[4] || (double)oz > g_qp_tile[5]) return 0;
    if (g_qp_axis_warp != NULL)      /* local umbilicus at the cube's centre z */
        AxisWarp_eval(g_qp_axis_warp, (double)oz + 64.0, &ay, &ax);
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 5; j++) {
            double y = (double)oy + 32.0 * (double)i;
            double x = (double)ox + 32.0 * (double)j;
            double dy = y - ay, dx = x - ax;
            double r = sqrt(dy * dy + dx * dx);
            double deg = atan2(dy, dx) * 180.0 / 3.14159265358979323846;
            if (deg < 0.0) deg += 360.0;
            if (r >= g_qp_tile[0] && r <= g_qp_tile[1] &&
                qp_angle_in_sector(deg, g_qp_tile[2], g_qp_tile[3]))
                return 1;
        }
    }
    return 0;
}

static int qp_tile_selftest(void)
{
    int fails = 0;
    g_qp_tile_armed = 1;
    g_qp_tile_axis[0] = 3405.0; g_qp_tile_axis[1] = 2878.0;
    /* band 400..600, sector 0..90 deg (dy >= 0, dx >= 0), z 4352..5120 */
    g_qp_tile[0] = 400.0; g_qp_tile[1] = 600.0;
    g_qp_tile[2] = 0.0;   g_qp_tile[3] = 90.0;
    g_qp_tile[4] = 4352.0; g_qp_tile[5] = 5120.0;
    /* cube at ~(r 500, 45 deg): origin y = 3405+354-64, x = 2878+354-64 */
    if (!qp_in_tile(4480, 3695, 3168)) fails++;
    if (qp_in_tile(5248, 3695, 3168)) fails++;          /* z out */
    if (qp_in_tile(4480, 3405 - 354 - 64, 3168)) fails++; /* 315 deg: out */
    if (qp_in_tile(4480, 3405 + 1000, 2878 + 1000)) fails++; /* r ~1414 out */
    /* the cube containing the axis is in every sector when the band starts at 0 */
    g_qp_tile[0] = 0.0; g_qp_tile[1] = 300.0;
    g_qp_tile[2] = 200.0; g_qp_tile[3] = 250.0;
    if (!qp_in_tile(4480, 3328, 2816)) fails++;
    /* wrapping sector 300..30 deg contains 10 deg and excludes 180 deg */
    g_qp_tile[0] = 400.0; g_qp_tile[1] = 600.0;
    g_qp_tile[2] = 300.0; g_qp_tile[3] = 30.0;
    if (!qp_in_tile(4480, 3405 + 87 - 64, 2878 + 492 - 64)) fails++;  /* 10 deg */
    if (qp_in_tile(4480, 3405 - 64, 2878 - 500 - 64)) fails++;        /* 180 deg */
    g_qp_tile_armed = 0;
    fprintf(stderr, "[selftest] quadribbon tile %s (%d failures)\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
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
    if (has_id && !qp_in_subgrid(oz, oy, ox)) return 0;  /* block filter */
    if (has_id && !qp_in_tile(oz, oy, ox)) return 0;     /* cylindrical tile */
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

static void qp_kibble_accumulate(MeshKibbleStats *total,
                                  const MeshKibbleStats *part)
{
    total->components += part->components;
    total->removed_components += part->removed_components;
    total->orphan_vertices += part->orphan_vertices;
    total->input_vertices += part->input_vertices;
    total->input_faces += part->input_faces;
    total->kept_vertices += part->kept_vertices;
    total->kept_faces += part->kept_faces;
    total->input_area += part->input_area;
    total->removed_area += part->removed_area;
}

static void qp_kibble_detail(FILE *f, const char *id, const MeshKibbleStats *s)
{
    if (f) fprintf(f, "%s\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%.17g\t%.17g\n",
        id, s->components, s->removed_components, s->orphan_vertices,
        s->input_vertices, s->kept_vertices, s->input_faces, s->kept_faces,
        s->input_area, s->removed_area);
}

static int qp_stage_concat(Arena_T arena, const char *mesh_dir,
                           const char *soup_path, const char *table_path,
                           double min_area, MeshKibbleStats *total, FILE *details)
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
        MeshBinData m, clean;
        MeshKibbleStats st;
        Arena_Mark mark = Arena_save(arena);
        memset(&m, 0, sizeof m);
        if (MeshBin_read_arena(arena, pile[i].path, &m) != 0) {
            fprintf(stderr, "[concat] cannot read %s\n", pile[i].path);
            free(pile);
            return -1;
        }
        if (MeshKibble_filter(arena, &m, min_area, &clean, &st) != 0) {
            free(pile); return -1;
        }
        total_v += clean.nv;
        total_f += clean.nf;
        Arena_restore(arena, mark);
    }
    if (total_v > INT32_MAX || total_v == 0 || total_f == 0) {
        fprintf(stderr, "[concat] cleanup leaves an empty or oversized mesh (min area %.9g)\n", min_area);
        free(pile); return -1;
    }
    verts = (float *)ARENA_ALLOC(arena, total_v * 3 * sizeof(float));
    faces = (int32_t *)ARENA_ALLOC(arena, total_f * 3 * sizeof(int32_t));
    if (table_path != NULL) {
        table = fopen(table_path, "wb");
        if (table != NULL)
            fprintf(table, "cube_id\tnv\tnf\tframe\toz\toy\tox\tpath\n");
    }
    for (size_t i = 0; i < np; i++) {
        MeshBinData m, clean;
        MeshKibbleStats st;
        Arena_Mark mark = Arena_save(arena);
        double oz = 0.0, oy = 0.0, ox = 0.0;
        const char *frame = "plain";
        memset(&m, 0, sizeof m);
        if (MeshBin_read_arena(arena, pile[i].path, &m) != 0) {
            if (table != NULL) fclose(table);
            free(pile);
            return -1;
        }
        if (MeshKibble_filter(arena, &m, min_area, &clean, &st) != 0) {
            if (table != NULL) fclose(table);
            free(pile); return -1;
        }
        qp_kibble_accumulate(total, &st);
        qp_kibble_detail(details, pile[i].cube_id[0] ? pile[i].cube_id : pile[i].path, &st);
        m = clean;
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
    char reach_s[32] = "", nsamp_s[32] = "";
    char world_buf[QP_MAX_PATH];
    vmesh = qp_world_twin(cfg, vmesh, world_buf, sizeof world_buf);
    if (!vmesh) return -1;
    snprintf(root_png, sizeof root_png, "%s/%s_bake.png", out_root, tag);
    if (qp_exists(root_png) && !qp_source_newer(vmesh, root_png)) {
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
    snprintf(reach_s, sizeof reach_s, "%.2f", cfg->bake_normal_reach);
    snprintf(nsamp_s, sizeof nsamp_s, "%d", (int)(2.0 * cfg->bake_normal_reach + 1.5));
    qp_tool_path(tool, sizeof tool, "obj_bake_raw");
    {
        const char *argv[32];
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
        argv[a++] = "--normal-range"; argv[a++] = reach_s;
        argv[a++] = "--normal-samples"; argv[a++] = nsamp_s;
        argv[a++] = "--raster-auto";
        argv[a] = NULL;
        if (qp_spawn("obj_bake_raw", argv) != 0) {
            fprintf(stderr, "[debug-bake] %s FAILED\n", tag);
            return -1;
        }
    }
    snprintf(stage_png, sizeof stage_png, "%s/%s_rawtex.png", bake_dir,
              tag);
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
    char world_buf[QP_MAX_PATH];
    ribbon_vmesh = qp_world_twin(cfg, ribbon_vmesh, world_buf, sizeof world_buf);
    if (!ribbon_vmesh) return -1;
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
                const char *argv[16];
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
                /* the sheet's v is arc length along the umbilicus, not world z:
                 * a v-edge spans sqrt(1+slope^2) voxels of surface per voxel of
                 * z, measured 1.01-1.07 across this scroll.  Only the strips
                 * (which is to say the bake) are remapped; the lattice and every
                 * gate keep v = z. */
                if (QS_ARCLEN_V && cfg->axis_table_path[0] != '\0') {
                    argv[a++] = "--axis-table";
                    argv[a++] = cfg->axis_table_path;
                }
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
        char reach_s[32] = "", nsamp_s[32] = "";
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
            snprintf(reach_s, sizeof reach_s, "%.2f", cfg->bake_normal_reach);
            snprintf(nsamp_s, sizeof nsamp_s, "%d", (int)(2.0 * cfg->bake_normal_reach + 1.5));
            qp_tool_path(tool, sizeof tool, "obj_bake_raw");
            {
                const char *argv[24];
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
                argv[a++] = "--normal-range"; argv[a++] = reach_s;
                argv[a++] = "--normal-samples"; argv[a++] = nsamp_s;
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

/* ---- untangle round: the collision shell (turn-order repair + elastic
 * shell + settle; radial motion only).  Topology and UV are unchanged, so the
 * refit that follows re-solves only the fills on the moved observations. */
static int qp_stage_untangle(const QpConfig *cfg, const char *in_vmesh,
                             const char *out_vmesh,
                             QuadribbonUntangleStats *out_stats)
{
    Arena_T arena = Arena_new();
    MeshBinData m;
    QuadribbonUntangleOpts o;
    QuadribbonUntangleStats st;
    int rc = -1;
    if (arena == NULL) return -1;
    memset(&m, 0, sizeof m);
    memset(&st, 0, sizeof st);
    if (MeshBin_read_arena(arena, in_vmesh, &m) != 0 || m.uv == NULL) {
        fprintf(stderr, "[untangle] cannot read %s (with UV)\n", in_vmesh);
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
            fprintf(stderr, "[untangle] FAILED\n");
            Arena_dispose(&arena);
            return -1;
        }
        qp_logf("[untangle] long %zu->%zu complete=%d accepted=%d "
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
    if (out_stats != NULL) *out_stats = st;
    return rc;
}

static int qp_write_untangle_json(const char *path,
                                  const QuadribbonUntangleStats *st)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f,
            "{ \"schema\": \"vesuvius-quadribbon-untangle-round-v1\",\n"
            "  \"input_long_conflicts\": %zu, \"output_long_conflicts\": %zu,\n"
            "  \"input_conflicts\": %zu, \"output_conflicts\": %zu,\n"
            "  \"complete\": %d, \"accepted_rounds\": %d, "
            "\"retained_partial\": %d,\n"
            "  \"contacts_built\": %zu, \"hard_contacts\": %zu,\n"
            "  \"pitch\": %.4f, \"clearance\": %.4f, \"movement_rms\": %.4f, "
            "\"movement_max\": %.4f,\n"
            "  \"settle_rounds_run\": %d, \"settle_accepted\": %d, "
            "\"turn_order_rounds\": %d }\n",
            st->input_long_conflicts, st->output_long_conflicts,
            st->input_conflicts, st->output_conflicts, st->complete,
            st->accepted_rounds, st->retained_partial, st->contacts_built,
            st->hard_contacts, st->pitch, st->clearance, st->movement_rms,
            st->movement_max, st->settle_rounds_run, st->settle_accepted,
            st->turn_order_rounds);
    return fclose(f) == 0 ? 0 : -1;
}

/* top-level integer member of a small JSON file, or def */
static long qp_read_json_long(const char *path, const char *key, long def)
{
    Arena_T arena = Arena_new();
    const char *err = NULL;
    const JsonValue *root = NULL;
    long v = def;
    if (arena == NULL) return def;
    root = Json_parse_file(arena, path, &err);
    if (root != NULL) v = Json_member_long(root, key, def);
    Arena_dispose(&arena);
    return v;
}

/* nested numeric member <obj>.<key> of a JSON file (obj NULL = top level) */
static double qp_read_json_double2(const char *path, const char *obj, const char *key, double def)
{
    Arena_T arena = Arena_new();
    const char *err = NULL;
    const JsonValue *root = NULL;
    double v = def;
    if (arena == NULL) return def;
    root = Json_parse_file(arena, path, &err);
    if (root != NULL) {
        const JsonValue *o = obj != NULL ? Json_object_get(root, obj) : root;
        if (o != NULL) v = Json_member_double(o, key, def);
    }
    Arena_dispose(&arena);
    return v;
}

static long long qp_file_size(const char *path)
{
    FILE *f = fopen(path, "rb");
    long long n = -1;
    if (f == NULL) return -1;
#ifdef _WIN32
    if (_fseeki64(f, 0, SEEK_END) == 0) n = _ftelli64(f);
#else
    if (fseek(f, 0, SEEK_END) == 0) n = (long long)ftell(f);
#endif
    fclose(f);
    return n;
}

static void qp_json_escape(const char *in, char *out, size_t cap)
{
    size_t o = 0;
    for (const char *p = in; *p != '\0' && o + 2 < cap; p++) {
        if (*p == '\\' || *p == '"') out[o++] = '\\';
        out[o++] = *p;
    }
    out[o] = '\0';
}

/* ---- sidecar plumbing ------------------------------------------------------
 * Every ribbon artifact travels with its provenance.  Consumers read the
 * sidecars BESIDE the mesh they consume (never a stale copy from an earlier
 * stage: the exact-length readers fail closed on any nv mismatch). */
static const char *const QP_SIDECARS[] = {
    "_support.u8", "_provenance.u8", "_phase.f32", "_material_identity.i32",
    "_reconstruction_component.i32", "_lane.i32", "_claimant_island.i32",
    "_claimant_chart.i32", "_atlas_column.i32", "_lattice_address.i32",
    "_dependency_kind.u8", "_dependency_endpoints.i32",
    "_stats.json", "_report.json", "_provenance.png"
};

static int qp_copy_sidecars(const char *src_stem, const char *dst_stem)
{
    for (size_t i = 0; i < sizeof QP_SIDECARS / sizeof QP_SIDECARS[0]; i++) {
        char src[QP_MAX_PATH], dst[QP_MAX_PATH];
        snprintf(src, sizeof src, "%s%s", src_stem, QP_SIDECARS[i]);
        snprintf(dst, sizeof dst, "%s%s", dst_stem, QP_SIDECARS[i]);
        if (!qp_exists(src)) continue;
        if (qp_copy_file(src, dst) != 0) {
            fprintf(stderr, "quadribbon: cannot copy %s -> %s\n", src, dst);
            return -1;
        }
    }
    return 0;
}

/* the metric stages discover sidecars by FIXED names beside their input */
static int qp_publish_plain_sidecars(const char *dir, const char *stem_base)
{
    static const char *const pairs[][2] = {
        { "_phase.f32", "ribbon_phase.f32" },
        { "_material_identity.i32", "ribbon_material_identity.i32" },
    };
    for (size_t i = 0; i < 2; i++) {
        char src[QP_MAX_PATH], dst[QP_MAX_PATH];
        snprintf(src, sizeof src, "%s/%s%s", dir, stem_base, pairs[i][0]);
        qp_join(dst, sizeof dst, dir, pairs[i][1]);
        if (qp_exists(src) && qp_copy_file(src, dst) != 0) {
            fprintf(stderr, "quadribbon: cannot publish %s\n", dst);
            return -1;
        }
    }
    return 0;
}

static int qp_publish_stem_plain_sidecars(const char *stem)
{
    const char *slash = NULL;
    char dir[QP_MAX_PATH];
    if (stem == NULL || stem[0] == '\0') return -1;
    for (const char *p = stem; *p != '\0'; p++)
        if (*p == '/' || *p == '\\') slash = p;
    if (slash == NULL) {
        snprintf(dir, sizeof dir, ".");
        return qp_publish_plain_sidecars(dir, stem);
    }
    if ((size_t)(slash - stem) >= sizeof dir) return -1;
    memcpy(dir, stem, (size_t)(slash - stem));
    dir[slash - stem] = '\0';
    return qp_publish_plain_sidecars(dir, slash + 1);
}

/* ---- verdict: the mechanical accept/reject gate on every ribbon ------------
 * Non-fatal by design: ribbon_verdict exits 4 on a failed gate and the
 * subprocess helper cannot tell that from a crash, so the JSON report is the
 * authority.  Every gate is logged; the run continues so the bakes exist. */
static int qp_verdict(const QpConfig *cfg, const char *tag, const char *vmesh,
                      const char *stem, const char *source_vmesh,
                      const char *report)
{
    char tool[QP_MAX_PATH], labels[QP_MAX_PATH], materials[QP_MAX_PATH];
    char support[QP_MAX_PATH], axis_y_s[32], axis_x_s[32], pitch_s[32];
    char core_s[32], spawn_tag[64];
    if (!qp_exists(report) || qp_source_newer(vmesh, report)) {
        const char *argv[24];
        size_t a = 0;
        qp_tool_path(tool, sizeof tool, "ribbon_verdict");
        snprintf(labels, sizeof labels, "%s_lane.i32", stem);
        if (!qp_exists(labels))
            snprintf(labels, sizeof labels, "%s_reconstruction_component.i32",
                     stem);
        snprintf(materials, sizeof materials, "%s_material_identity.i32", stem);
        snprintf(support, sizeof support, "%s_support.u8", stem);
        snprintf(axis_y_s, sizeof axis_y_s, "%.3f", cfg->axis_y);
        snprintf(axis_x_s, sizeof axis_x_s, "%.3f", cfg->axis_x);
        snprintf(pitch_s, sizeof pitch_s, "%.3f", cfg->wrap_spacing);
        snprintf(core_s, sizeof core_s, "%.1f", QP_VERDICT_CORE_RADIUS);
        argv[a++] = tool;
        argv[a++] = vmesh;
        if (qp_exists(labels)) { argv[a++] = "--labels"; argv[a++] = labels; }
        if (qp_exists(materials)) { argv[a++] = "--materials"; argv[a++] = materials; }
        if (qp_exists(support)) { argv[a++] = "--support"; argv[a++] = support; }
        argv[a++] = "--umb-y"; argv[a++] = axis_y_s;
        argv[a++] = "--umb-x"; argv[a++] = axis_x_s;
        argv[a++] = "--pitch"; argv[a++] = pitch_s;
        argv[a++] = "--core-radius"; argv[a++] = core_s;
        if (source_vmesh != NULL && qp_exists(source_vmesh)) {
            argv[a++] = "--source"; argv[a++] = source_vmesh;
        }
        argv[a++] = "--report"; argv[a++] = report;
        argv[a] = NULL;
        snprintf(spawn_tag, sizeof spawn_tag, "ribbon_verdict(%s)", tag);
        if (qp_spawn(spawn_tag, argv) != 0)
            qp_logf("[verdict %s] ribbon_verdict exited nonzero (a failed gate "
                    "or a crash; the report decides)\n", tag);
    } else {
        qp_logf("[verdict %s] resume: %s\n", tag, report);
    }
    {
        Arena_T arena = Arena_new();
        const char *err = NULL;
        const JsonValue *root = arena != NULL ? Json_parse_file(arena, report, &err)
                                              : NULL;
        if (root != NULL) {
            const char *verdict = Json_as_string(Json_object_get(root, "verdict"));
            long failed = Json_member_long(root, "gates_failed", -1);
            const JsonValue *gates = Json_object_get(root, "gates");
            size_t n = gates != NULL ? Json_array_len(gates) : 0;
            char names[768] = "";
            size_t off = 0;
            for (size_t i = 0; i < n; i++) {
                const JsonValue *g = Json_array_get(gates, i);
                const char *name = Json_as_string(Json_object_get(g, "name"));
                if (Json_as_bool(Json_object_get(g, "pass"), 1)) continue;
                if (name != NULL && off + strlen(name) + 3 < sizeof names)
                    off += (size_t)snprintf(names + off, sizeof names - off, "%s%s",
                                            off ? ", " : "", name);
            }
            qp_logf("[verdict %s] %s: %ld of %zu gates failed%s%s\n", tag,
                    verdict != NULL ? verdict : "?", failed, n,
                    names[0] ? " -- " : "", names);
        } else {
            qp_logf("[verdict %s] no report at %s (verdict unavailable)\n", tag,
                    report);
        }
        if (arena != NULL) Arena_dispose(&arena);
    }
    return 0;
}

/* ---- stage 1: the meshed input ----------------------------------------- */

static int qp_input_is_vmesh(const char *input)
{
    size_t n = strlen(input);
    return n > 6 && strcmp(input + n - 6, ".vmesh") == 0 &&
           MeshBin_looks_complete(input);
}

/* ---- core curl routing --------------------------------------------------
 * The umbilicus of a rolled scroll is a hollow tube whose wall is the
 * innermost wrap; the crumpled END of the roll sits inside that void as
 * separate mesh components (measured on the PHerc0139 5x3x3 core, 2026-09-02:
 * 45 components, 26% of the core mesh, r_max < 102 while the wall starts at
 * ~100).  Left in, the register chains them into the wall's turn and the fit
 * spends its coverage on them (r100-130 coverage 73.8% -> 99.1% once they are
 * out; distortion 0.36% -> 0.002%, coherence 50.6% -> 71.3%).  They are NOT
 * flattenable with the sheet: they belong in extras, labelled.  Rule: a
 * component whose maximum cylindrical radius lies below core.wall_radius_voxels
 * is entirely inside the wall -> routed.  The wall radius is a per-scroll
 * datum of the same kind as the umbilicus (config, never a flag).
 * --------------------------------------------------------------------------- */

typedef struct {
    size_t comps_total, comps_routed, verts_routed, faces_routed;
} QpCurlStats;

/* Classify every vertex of a soup: out_keep[v] = 1 when the vertex's component
 * reaches r_max >= wall_r (sheet material), 0 when it is routed. */
static int qp_core_curl_classify(Arena_T arena, const float *verts, size_t nv,
                                 const int32_t *faces, size_t nf,
                                 double axis_y, double axis_x, double wall_r,
                                 uint8_t *out_keep, QpCurlStats *st,
                                 FILE *ledger)
{
    UnionFind uf;
    int32_t *root = NULL;
    double *rmax = NULL;
    int32_t *cnv = NULL;
    size_t v = 0, f = 0;
    memset(st, 0, sizeof *st);
    if (nv == 0) return 0;
    if (nv > (size_t)INT32_MAX) return -1;
    uf = UF_new(arena, (int32_t)nv);
    for (f = 0; f < nf; f++) {
        uf_union(&uf, faces[f * 3 + 0], faces[f * 3 + 1]);
        uf_union(&uf, faces[f * 3 + 0], faces[f * 3 + 2]);
    }
    root = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    rmax = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    cnv = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (v = 0; v < nv; v++) { rmax[v] = -1.0; cnv[v] = 0; }
    for (v = 0; v < nv; v++) {
        double dy = (double)verts[v * 3 + 1] - axis_y;
        double dx = (double)verts[v * 3 + 2] - axis_x;
        double r = sqrt(dy * dy + dx * dx);
        int32_t rt = uf_find(&uf, (int32_t)v);
        root[v] = rt;
        if (r > rmax[rt]) rmax[rt] = r;
        cnv[rt]++;
    }
    for (v = 0; v < nv; v++) {
        if (root[v] == (int32_t)v) {
            st->comps_total++;
            if (rmax[v] < wall_r) {
                st->comps_routed++;
                st->verts_routed += (size_t)cnv[v];
                if (ledger != NULL)
                    fprintf(ledger, "%s    { \"root\": %zu, \"nv\": %d, "
                            "\"r_max\": %.2f }",
                            st->comps_routed > 1 ? ",\n" : "\n", v, cnv[v],
                            rmax[v]);
            }
        }
    }
    for (v = 0; v < nv; v++) out_keep[v] = rmax[root[v]] < wall_r ? 0 : 1;
    for (f = 0; f < nf; f++)
        if (!out_keep[faces[f * 3]]) st->faces_routed++;
    return 0;
}

/* Rewrite the concat's per-cube table after routing: kept vertices/faces stay
 * in their original order, so every cube's block is still contiguous and only
 * the nv/nf columns change.  Rows are "cube_id nv nf frame oz oy ox path". */
static int qp_core_curl_rewrite_table(const char *table_path,
                                      const uint8_t *keep, size_t nv,
                                      const int32_t *faces, size_t nf)
{
    char tmp[QP_MAX_PATH];
    char line[QP_MAX_PATH * 2];
    FILE *in = NULL, *out = NULL;
    size_t at_v = 0, at_f = 0;
    if (table_path == NULL || !qp_exists(table_path)) return 0;
    snprintf(tmp, sizeof tmp, "%s.tmp", table_path);
    in = fopen(table_path, "rb");
    if (in == NULL) return -1;
    out = fopen(tmp, "wb");
    if (out == NULL) { fclose(in); return -1; }
    if (fgets(line, sizeof line, in) != NULL) fputs(line, out);   /* header */
    while (fgets(line, sizeof line, in) != NULL) {
        char *t1 = strchr(line, '\t');
        char *t2 = t1 != NULL ? strchr(t1 + 1, '\t') : NULL;
        char *t3 = t2 != NULL ? strchr(t2 + 1, '\t') : NULL;
        unsigned long long cnv = 0, cnf = 0;
        size_t kv = 0, kf = 0, k = 0;
        if (t1 == NULL || t2 == NULL || t3 == NULL) { fputs(line, out); continue; }
        cnv = strtoull(t1 + 1, NULL, 10);
        cnf = strtoull(t2 + 1, NULL, 10);
        if (at_v + cnv > nv || at_f + cnf > nf) {
            fclose(in); fclose(out); remove(tmp);
            return -1;
        }
        for (k = 0; k < cnv; k++) kv += keep[at_v + k];
        for (k = 0; k < cnf; k++) kf += keep[faces[(at_f + k) * 3]];
        *t1 = '\0';
        fprintf(out, "%s\t%zu\t%zu%s", line, kv, kf, t3);
        at_v += cnv;
        at_f += cnf;
    }
    fclose(in);
    fclose(out);
    remove(table_path);
    return rename(tmp, table_path) == 0 ? 0 : -1;
}

/* Apply the routing to an on-disk soup: <stem>.vmesh keeps the sheet
 * material, <stage_dir>/mesh_core_curl_extras.vmesh receives the routed
 * components, mesh_core_curl.json is the ledger.  No-op when wall_r <= 0. */
static int qp_route_core_curl(Arena_T arena, const QpConfig *cfg,
                              const char *stage_dir, const char *mesh_vmesh,
                              const char *table_path)
{
    MeshBinData m;
    Arena_Mark mark;
    uint8_t *keep = NULL;
    QpCurlStats st;
    char extras[QP_MAX_PATH], ledger_path[QP_MAX_PATH];
    FILE *ledger = NULL;
    int rc = -1;
    if (cfg->core_wall_radius <= 0.0) return 0;
    mark = Arena_save(arena);
    memset(&m, 0, sizeof m);
    if (MeshBin_read_arena(arena, mesh_vmesh, &m) != 0) {
        qp_logf("[stage 1 mesh] core curl: cannot read %s\n", mesh_vmesh);
        Arena_restore(arena, mark);
        return -1;
    }
    qp_join(ledger_path, sizeof ledger_path, stage_dir, "mesh_core_curl.json");
    qp_join(extras, sizeof extras, stage_dir, "mesh_core_curl_extras.vmesh");
    ledger = fopen(ledger_path, "wb");
    if (ledger != NULL)
        fprintf(ledger, "{ \"schema\": \"vesuvius-quadribbon-core-curl-v1\",\n"
                        "  \"wall_radius_voxels\": %.2f,\n"
                        "  \"rule\": \"component r_max < wall_radius -> extras "
                        "(crumpled roll end inside the umbilicus void)\",\n"
                        "  \"routed_components\": [",
                cfg->core_wall_radius);
    keep = (uint8_t *)ARENA_ALLOC(arena, m.nv > 0 ? m.nv : 1);
    if (qp_core_curl_classify(arena, m.verts, m.nv, m.faces, m.nf,
                              cfg->axis_y, cfg->axis_x, cfg->core_wall_radius,
                              keep, &st, ledger) != 0) {
        qp_logf("[stage 1 mesh] core curl: classify failed (nv=%zu)\n", m.nv);
        if (ledger != NULL) fclose(ledger);
        Arena_restore(arena, mark);
        return -1;
    }
    if (ledger != NULL) {
        fprintf(ledger, "\n  ],\n  \"components_total\": %zu, "
                        "\"components_routed\": %zu,\n"
                        "  \"verts_total\": %zu, \"verts_routed\": %zu, "
                        "\"faces_total\": %zu, \"faces_routed\": %zu,\n"
                        "  \"extras\": \"mesh_core_curl_extras.vmesh\" }\n",
                st.comps_total, st.comps_routed, m.nv, st.verts_routed, m.nf,
                st.faces_routed);
        fclose(ledger);
    }
    qp_logf("[stage 1 mesh] core curl: wall_radius=%.1f routed %zu/%zu "
            "components, %zu verts (%.1f%%), %zu faces -> %s\n",
            cfg->core_wall_radius, st.comps_routed, st.comps_total,
            st.verts_routed,
            m.nv > 0 ? 100.0 * (double)st.verts_routed / (double)m.nv : 0.0,
            st.faces_routed, st.comps_routed > 0 ? extras : "(nothing)");
    if (st.comps_routed == 0) { Arena_restore(arena, mark); return 0; }
    {
        /* split: kept side rewrites the soup, routed side becomes extras */
        int32_t *newk = (int32_t *)ARENA_ALLOC(arena, m.nv * sizeof(int32_t));
        int32_t *newr = (int32_t *)ARENA_ALLOC(arena, m.nv * sizeof(int32_t));
        size_t nk = 0, nr = 0, fk = 0, fr = 0, v = 0, f = 0;
        float *vk = NULL, *vr = NULL, *uvk = NULL;
        int32_t *fkp = NULL, *frp = NULL;
        for (v = 0; v < m.nv; v++) {
            if (keep[v]) { newk[v] = (int32_t)nk++; newr[v] = -1; }
            else { newr[v] = (int32_t)nr++; newk[v] = -1; }
        }
        vk = (float *)ARENA_ALLOC(arena, (nk > 0 ? nk : 1) * 3 * sizeof(float));
        vr = (float *)ARENA_ALLOC(arena, (nr > 0 ? nr : 1) * 3 * sizeof(float));
        if (m.uv != NULL)
            uvk = (float *)ARENA_ALLOC(arena, (nk > 0 ? nk : 1) * 2 * sizeof(float));
        fkp = (int32_t *)ARENA_ALLOC(arena, (m.nf > 0 ? m.nf : 1) * 3 * sizeof(int32_t));
        frp = (int32_t *)ARENA_ALLOC(arena, (m.nf > 0 ? m.nf : 1) * 3 * sizeof(int32_t));
        for (v = 0; v < m.nv; v++) {
            if (keep[v]) {
                memcpy(vk + (size_t)newk[v] * 3, m.verts + v * 3, 3 * sizeof(float));
                if (uvk != NULL)
                    memcpy(uvk + (size_t)newk[v] * 2, m.uv + v * 2, 2 * sizeof(float));
            } else {
                memcpy(vr + (size_t)newr[v] * 3, m.verts + v * 3, 3 * sizeof(float));
            }
        }
        for (f = 0; f < m.nf; f++) {
            const int32_t *tri = m.faces + f * 3;
            if (keep[tri[0]]) {
                fkp[fk * 3 + 0] = newk[tri[0]];
                fkp[fk * 3 + 1] = newk[tri[1]];
                fkp[fk * 3 + 2] = newk[tri[2]];
                fk++;
            } else {
                frp[fr * 3 + 0] = newr[tri[0]];
                frp[fr * 3 + 1] = newr[tri[1]];
                frp[fr * 3 + 2] = newr[tri[2]];
                fr++;
            }
        }
        rc = 0;
        if (nr > 0 && MeshBin_write(extras, vr, nr, frp, fr, NULL) != 0) rc = -1;
        if (rc == 0 && nk > 0 &&
            MeshBin_write(mesh_vmesh, vk, nk, fkp, fk, uvk) != 0)
            rc = -1;
        if (rc == 0 && nk == 0) {
            qp_logf("[stage 1 mesh] core curl: EVERYTHING routed; wall radius "
                    "%.1f is wrong for this region\n", cfg->core_wall_radius);
            rc = -1;
        }
        if (rc == 0 && qp_core_curl_rewrite_table(table_path, keep, m.nv,
                                                  m.faces, m.nf) != 0) {
            qp_logf("[stage 1 mesh] core curl: cube table rewrite failed\n");
            rc = -1;
        }
    }
    Arena_restore(arena, mark);
    return rc;
}

/* selftest: a ring of triangles at r=150 (sheet) and a small fan at r=40
 * (curl) around axis (0,0); wall radius 100 must route exactly the fan. */
static int qp_core_curl_selftest(void)
{
    Arena_T arena = Arena_new();
    enum { RING_N = 32, FAN_N = 8 };
    float verts[(RING_N * 2 + FAN_N + 1) * 3];
    int32_t faces[(RING_N * 2 + FAN_N) * 3];
    uint8_t keep[RING_N * 2 + FAN_N + 1];
    QpCurlStats st;
    size_t nv = 0, nf = 0, i = 0;
    int fails = 0;
    if (arena == NULL) return 1;
    for (i = 0; i < RING_N; i++) {
        double a = 6.283185307179586 * (double)i / RING_N;
        verts[nv * 3 + 0] = 0.0f;
        verts[nv * 3 + 1] = (float)(150.0 * sin(a));
        verts[nv * 3 + 2] = (float)(150.0 * cos(a));
        nv++;
        verts[nv * 3 + 0] = 10.0f;
        verts[nv * 3 + 1] = (float)(150.0 * sin(a));
        verts[nv * 3 + 2] = (float)(150.0 * cos(a));
        nv++;
    }
    for (i = 0; i < RING_N; i++) {
        int32_t a0 = (int32_t)(2 * i), a1 = a0 + 1;
        int32_t b0 = (int32_t)(2 * ((i + 1) % RING_N)), b1 = b0 + 1;
        faces[nf * 3 + 0] = a0; faces[nf * 3 + 1] = a1; faces[nf * 3 + 2] = b0; nf++;
        faces[nf * 3 + 0] = a1; faces[nf * 3 + 1] = b1; faces[nf * 3 + 2] = b0; nf++;
    }
    {
        size_t c = nv;
        verts[nv * 3 + 0] = 5.0f; verts[nv * 3 + 1] = 0.0f; verts[nv * 3 + 2] = 0.0f;
        nv++;
        for (i = 0; i < FAN_N; i++) {
            double a = 6.283185307179586 * (double)i / FAN_N;
            verts[nv * 3 + 0] = 5.0f;
            verts[nv * 3 + 1] = (float)(40.0 * sin(a));
            verts[nv * 3 + 2] = (float)(40.0 * cos(a));
            nv++;
        }
        for (i = 0; i < FAN_N; i++) {
            faces[nf * 3 + 0] = (int32_t)c;
            faces[nf * 3 + 1] = (int32_t)(c + 1 + i);
            faces[nf * 3 + 2] = (int32_t)(c + 1 + (i + 1) % FAN_N);
            nf++;
        }
    }
    if (qp_core_curl_classify(arena, verts, nv, faces, nf, 0.0, 0.0, 100.0,
                              keep, &st, NULL) != 0) fails++;
    if (st.comps_total != 2 || st.comps_routed != 1 ||
        st.verts_routed != FAN_N + 1 || st.faces_routed != FAN_N) fails++;
    for (i = 0; i < nv; i++)
        if ((i < RING_N * 2) != (keep[i] != 0)) { fails++; break; }
    /* wall radius above everything routes both; zero routes nothing */
    if (qp_core_curl_classify(arena, verts, nv, faces, nf, 0.0, 0.0, 1000.0,
                              keep, &st, NULL) != 0 || st.comps_routed != 2)
        fails++;
    if (qp_core_curl_classify(arena, verts, nv, faces, nf, 0.0, 0.0, 0.0,
                              keep, &st, NULL) != 0 || st.comps_routed != 0)
        fails++;
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] quadribbon core_curl %s (%d failures)\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}

/* Stage 1b: straighten the soup about the axis curve.  mesh_world.vmesh keeps
 * the world coordinates; mesh.vmesh becomes the straightened drop-in that
 * every later stage (routing included) consumes.  No-op without an axis
 * table. */
static int qp_straighten_mesh(const QpConfig *cfg, const char *stage_dir,
                              const char *mesh_vmesh)
{
    char world[QP_MAX_PATH], tmp[QP_MAX_PATH];
    double y0 = 0.0, x0 = 0.0, y1 = 0.0, x1 = 0.0;
    if (!cfg->axis_warp_armed) return 0;
    qp_join(world, sizeof world, stage_dir, "mesh_world.vmesh");
    if (MeshBin_looks_complete(world)) {
        qp_logf("[stage 1 mesh] axis: already straightened (%s exists)\n", world);
        return 0;
    }
    snprintf(tmp, sizeof tmp, "%s.straight.tmp", mesh_vmesh);
    if (qp_write_warped_copy(cfg, mesh_vmesh, tmp, +1) != 0) return -1;
    if (qp_copy_file(mesh_vmesh, world) != 0) return -1;
    if (qp_replace_file(tmp, mesh_vmesh) != 0) return -1;
    if (cfg->axis_warp.physical) {
        qp_logf("[stage 1 mesh] physical arc: %zu ordered shaft samples, %.6f mm arc, %.9g um per metric voxel; world mesh kept\n",
                cfg->axis_warp.n,ShaftWarp_length_um(cfg->axis_warp.physical)/1000.,
                cfg->axis_warp.physical_config.metric_voxel_um);
        return 0;
    }
    AxisWarp_eval(&cfg->axis_warp, cfg->axis_warp.z[0], &y0, &x0);
    AxisWarp_eval(&cfg->axis_warp, cfg->axis_warp.z[cfg->axis_warp.n - 1], &y1, &x1);
    qp_logf("[stage 1 mesh] axis: straightened about %s (%zu samples, z %.0f..%.0f, "
            "curve (%.0f,%.0f)..(%.0f,%.0f) -> reference (%.1f,%.1f)); world mesh kept\n",
            cfg->axis_table_path, cfg->axis_warp.n, cfg->axis_warp.z[0],
            cfg->axis_warp.z[cfg->axis_warp.n - 1], y0, x0, y1, x1,
            cfg->axis_warp.reference_y, cfg->axis_warp.reference_x);
    return 0;
}

static int qp_stage_mesh(Arena_T arena, const QpConfig *cfg, const char *input,
                         int is_vmesh, const char *stage_dir,
                         const char *mesh_vmesh)
{
    char kibble_report[QP_MAX_PATH], kibble_details[QP_MAX_PATH];
    MeshKibbleStats total = {0};
    FILE *details = NULL;
    int rc = -1;
    qp_mkdir(stage_dir);
    qp_join(kibble_report, sizeof kibble_report, stage_dir, "mesh_kibble.json");
    qp_join(kibble_details, sizeof kibble_details, stage_dir, "mesh_kibble_inputs.tsv");
    if (MeshBin_looks_complete(mesh_vmesh)) {
        const char *err = NULL;
        Arena_Mark mark = Arena_save(arena);
        const JsonValue *report = Json_parse_file(arena, kibble_report, &err);
        const char *schema = report ? Json_as_string(Json_object_get(report, "schema")) : NULL;
        int matches = report && schema && !strcmp(schema, "vesuvius-pre-unwrap-kibble-v1") &&
            Json_member_double(report, "min_component_area_vox2", -1) == cfg->kibble_min_area;
        /* A legacy run is reproducible only with cleanup explicitly disabled.
         * Never silently reuse an unfiltered certificate or change its indices. */
        if (!matches && (report || cfg->kibble_min_area > 0 || cfg->axis_warp_armed || qp_exists(kibble_report))) {
            qp_logf("[stage 1 mesh] cleanup policy is missing or differs from this run; use a fresh output directory (mesh_cleanup.min_component_area_vox2=%.9g)\n", cfg->kibble_min_area);
            Arena_restore(arena, mark);
            return -1;
        }
        Arena_restore(arena, mark);
        qp_logf("[stage 1 mesh] resume: %s already complete\n", mesh_vmesh);
        return 0;
    }
    details = fopen(kibble_details, "wb");
    if (!details) return -1;
    fprintf(details, "input\tcomponents\tremoved_components\torphan_vertices\tinput_vertices\tkept_vertices\tinput_faces\tkept_faces\tinput_area_vox2\tremoved_area_vox2\n");
    if (is_vmesh) {
        char ledger[QP_MAX_PATH], esc[QP_MAX_PATH * 2];
        FILE *f = NULL;
        Arena_Mark mark = Arena_save(arena);
        MeshBinData m, clean;
        if (MeshBin_read_arena(arena, input, &m) != 0 ||
            MeshKibble_filter(arena, &m, cfg->kibble_min_area, &clean, &total) != 0 ||
            !clean.nv || !clean.nf ||
            MeshBin_write(mesh_vmesh, clean.verts, clean.nv, clean.faces, clean.nf, clean.uv) != 0) {
            fprintf(stderr, "[stage 1 mesh] cannot copy %s\n", input);
            Arena_restore(arena, mark);
            fclose(details);
            return -1;
        }
        Arena_restore(arena, mark);
        qp_kibble_detail(details, input, &total);
        qp_join(ledger, sizeof ledger, stage_dir, "mesh_source.json");
        qp_json_escape(input, esc, sizeof esc);
        f = fopen(ledger, "wb");
        if (f != NULL) {
            fprintf(f, "{ \"schema\": \"vesuvius-quadribbon-mesh-source-v1\", "
                       "\"kind\": \"welded_vmesh\", \"path\": \"%s\", "
                       "\"bytes\": %lld }\n", esc, qp_file_size(input));
            fclose(f);
        }
        qp_logf("[stage 1 mesh] welded mesh installed: %s\n", input);
        rc = qp_straighten_mesh(cfg, stage_dir, mesh_vmesh);
        if (rc == 0) rc = qp_route_core_curl(arena, cfg, stage_dir, mesh_vmesh, NULL);
    } else {
        char table[QP_MAX_PATH];
        Arena_T soup_arena = Arena_new();
        qp_join(table, sizeof table, stage_dir, "mesh_cubes.tsv");
        if (soup_arena != NULL)
            rc = qp_stage_concat(soup_arena, input, mesh_vmesh, table,
                                 cfg->kibble_min_area, &total, details);
        if (soup_arena != NULL) Arena_dispose(&soup_arena);
        if (g_qp_tile_geometry_error) { qp_logf("[tile] unresolved physical-coordinate projection\n"); rc=-1; }
        if (rc == 0) rc = qp_straighten_mesh(cfg, stage_dir, mesh_vmesh);
        if (rc == 0)
            rc = qp_route_core_curl(arena, cfg, stage_dir, mesh_vmesh, table);
    }
    if (fclose(details) != 0) return -1;
    if (rc != 0) return rc;
    /* This marker is published only after filtering, straightening and routing
     * have completed, so resume cannot accept partially prepared geometry. */
    char temp[QP_MAX_PATH];
    snprintf(temp, sizeof temp, "%s.tmp", kibble_report);
    FILE *report = fopen(temp, "wb");
    if (!report) return -1;
    fprintf(report, "{\n  \"schema\": \"vesuvius-pre-unwrap-kibble-v1\",\n"
        "  \"enabled\": %s,\n  \"min_component_area_vox2\": %.17g,\n"
        "  \"policy\": \"absolute area per original input mesh, before axis warp; surviving indices remain in source order\",\n"
        "  \"details\": \"mesh_kibble_inputs.tsv\",\n"
        "  \"components\": %zu,\n  \"removed_components\": %zu,\n  \"orphan_vertices\": %zu,\n"
        "  \"input_vertices\": %zu,\n  \"kept_vertices\": %zu,\n"
        "  \"input_faces\": %zu,\n  \"kept_faces\": %zu,\n"
        "  \"input_area_vox2\": %.17g,\n  \"removed_area_vox2\": %.17g\n}\n",
        cfg->kibble_min_area > 0 ? "true" : "false", cfg->kibble_min_area,
        total.components, total.removed_components, total.orphan_vertices,
        total.input_vertices, total.kept_vertices, total.input_faces, total.kept_faces,
        total.input_area, total.removed_area);
    if (fclose(report) != 0 || rename(temp, kibble_report) != 0) return -1;
    qp_logf("[stage 1 mesh] kibble: removed %zu/%zu components, %zu orphan vertices, %zu faces, %.6g voxel^2 (%.4f%% area); min area %.9g\n",
        total.removed_components, total.components, total.orphan_vertices,
        total.input_faces-total.kept_faces, total.removed_area,
        total.input_area > 0 ? 100*total.removed_area/total.input_area : 0,
        cfg->kibble_min_area);
    return 0;
}

/* ---- stages 2/3: scroll_ribbon spawns ----------------------------------- */

static int qp_spawn_scroll_ribbon(const QpConfig *cfg, const char *in_vmesh,
                                  const char *out_dir, const char *id, int fit,
                                  int trust_gauge,
                                  const char *boundary_winding,
                                  const char *boundary_u,
                                  const char *boundary_v,
                                  const char *boundary_material,
                                  const char *label)
{
    char tool[QP_MAX_PATH], axis_y_s[32], axis_x_s[32], pitch_s[32];
    char threads_s[16], grid_u_s[16], slice_h_s[16], sample_h_s[16];
    char iters_s[8], final_s[8], metric_s[8], peel_layers_s[8];
    char peel_share_s[32], peel_vertices_s[32];
    const char *argv[64];
    size_t a = 0;
    qp_tool_path(tool, sizeof tool, "scroll_ribbon");
    snprintf(axis_y_s, sizeof axis_y_s, "%.3f", cfg->axis_y);
    snprintf(axis_x_s, sizeof axis_x_s, "%.3f", cfg->axis_x);
    snprintf(pitch_s, sizeof pitch_s, "%.3f", cfg->wrap_spacing);
    snprintf(threads_s, sizeof threads_s, "%ld", cfg->threads);
    snprintf(grid_u_s, sizeof grid_u_s, "%.6g", QP_FIT_GRID_U);
    snprintf(slice_h_s, sizeof slice_h_s, "%.6g", QP_FIT_SLICE_H);
    snprintf(sample_h_s, sizeof sample_h_s, "%.6g", QP_FIT_SAMPLE_H);
    snprintf(iters_s, sizeof iters_s, "%d", QP_FIT_ITERS);
    snprintf(final_s, sizeof final_s, "%d", QP_FIT_FINAL_ITERS);
    snprintf(metric_s, sizeof metric_s, "%d", QP_FIT_METRIC_ITERS);
    snprintf(peel_layers_s, sizeof peel_layers_s, "%d", QP_FIT_PEEL_LAYERS);
    snprintf(peel_share_s, sizeof peel_share_s, "%.9g",
             QP_FIT_PROMOTE_PEEL_MIN_SHARE);
    snprintf(peel_vertices_s, sizeof peel_vertices_s, "%d",
             QP_FIT_PROMOTE_PEEL_MIN_VERTICES);
    argv[a++] = tool;
    argv[a++] = in_vmesh;
    argv[a++] = out_dir;
    argv[a++] = "--id"; argv[a++] = id;
    if (fit) {
        argv[a++] = "--scaffold-solve";
        if (trust_gauge) argv[a++] = "--trust-gauge";
        argv[a++] = "--component-global";
        argv[a++] = "--vmesh-only";
        argv[a++] = "--peel-layers";
        argv[a++] = peel_layers_s;
        argv[a++] = "--promote-peel-min-share";
        argv[a++] = peel_share_s;
        argv[a++] = "--promote-peel-min-vertices";
        argv[a++] = peel_vertices_s;
        argv[a++] = "--grid-u"; argv[a++] = grid_u_s;
        argv[a++] = "--slice-h"; argv[a++] = slice_h_s;
        argv[a++] = "--sample-h"; argv[a++] = sample_h_s;
        argv[a++] = "--iters"; argv[a++] = iters_s;
        argv[a++] = "--final-iters"; argv[a++] = final_s;
        argv[a++] = "--metric-iters"; argv[a++] = metric_s;
    } else {
        argv[a++] = "--winding-only";
        argv[a++] = "--vmesh-only";
    }
    argv[a++] = "--axis-point"; argv[a++] = "0"; argv[a++] = axis_y_s; argv[a++] = axis_x_s;
    argv[a++] = "--axis-dir"; argv[a++] = "1"; argv[a++] = "0"; argv[a++] = "0";
    argv[a++] = "--wrap-spacing"; argv[a++] = pitch_s;
    if (cfg->winding_sense != 0) {
        argv[a++] = "--winding-sense";
        argv[a++] = cfg->winding_sense > 0 ? "1" : "-1";
    }
    if (!cfg->overlap_family) argv[a++] = "--no-cube-table";
    if (boundary_winding != NULL) {
        argv[a++] = "--boundary-winding"; argv[a++] = boundary_winding;
        argv[a++] = "--boundary-u"; argv[a++] = boundary_u;
        argv[a++] = "--boundary-v"; argv[a++] = boundary_v;
        argv[a++] = "--boundary-material"; argv[a++] = boundary_material;
    }
    argv[a++] = "--threads"; argv[a++] = threads_s;
    argv[a] = NULL;
    return qp_spawn(label, argv);
}

/* ---- projective canonical blocks -----------------------------------------
 *
 * A crop-wide winding solve is not projective: adding one cube row changes
 * the global forest/MRF and can relabel vertices which received no new
 * evidence.  Canonical blocks make the finite solve domain part of the
 * algorithm instead.  The root block is the measured 4x5x5 domain about the
 * umbilicus.  Neighbours overlap by one cube, and every block has exactly one
 * canonical parent towards the root.  A parent is immutable; shared cubes
 * place only the child.  Consequently a requested crop and every extension of
 * it read the same source cube from the same owner block.
 *
 * The block lattice is expressed in source cube origins.  Z has four cube
 * rows with stride three; Y/X have five with stride four.  The root Y/X block
 * is centred on the configured straightened axis.  Its Z origin is the cube
 * containing the first axis-table sample, or an explicit per-scroll config
 * value.  It is never inferred from the requested crop. */
enum {
    QP_CB_Z_SIZE = 4, QP_CB_Y_SIZE = 5, QP_CB_X_SIZE = 5,
    QP_CB_Z_STRIDE = 3, QP_CB_Y_STRIDE = 4, QP_CB_X_STRIDE = 4,
    /* Reversible, collision-free 31-bit component key:
     * [ z:7 | y:6 | x:6 | local-label:12 ].  Bounds are checked; there is no
     * hash or crop-dependent collision resolution. */
    QP_CB_LABEL_LOCAL_BITS = 12,
    QP_CB_LABEL_LOCAL_MAX = (1 << QP_CB_LABEL_LOCAL_BITS) - 1,
    QP_CB_LABEL_Z_MIN = -64, QP_CB_LABEL_Z_MAX = 63,
    QP_CB_LABEL_Y_MIN = -32, QP_CB_LABEL_Y_MAX = 31,
    QP_CB_LABEL_X_MIN = -32, QP_CB_LABEL_X_MAX = 31
};

typedef struct {
    char id[24];
    size_t nv, nf, v0, f0;
    long oz, oy, ox;
} QpCubeRow;

typedef struct {
    QpCubeRow *row;
    size_t n, nv, nf;
} QpCubeTable;

typedef struct {
    int iz, iy, ix;
    int parent;                 /* plan index; -1 for the immutable root */
    long bounds[6];             /* inclusive z0 z1 y0 y1 x0 x1 */
    char id[64];
    char root[QP_MAX_PATH];
    int ready;
    int32_t *relation_global;
    int32_t *material_global;
    int32_t *mesh_global;
    int32_t *turn_offset;       /* indexed by local relation island */
    double *u_offset;           /* indexed by local relation island */
    size_t nrelation, nmaterial, nmesh;
    double v_offset;
    size_t shared_vertices;
    double turn_agreement;
} QpCanonicalBlock;

typedef struct {
    MeshBinData cert;
    QpCubeTable cubes;
    float *winding, *confidence, *field, *jump;
    int32_t *relation, *material, *mesh_component;
} QpBlockArtifact;

static void qp_cube_table_dispose(QpCubeTable *t)
{
    if (t == NULL) return;
    free(t->row);
    memset(t, 0, sizeof *t);
}

static int qp_read_cube_table(const char *path, QpCubeTable *out)
{
    FILE *f = NULL;
    char line[QP_MAX_PATH * 2];
    size_t cap = 0;
    QpCubeTable t;
    memset(&t, 0, sizeof t);
    f = fopen(path, "rb");
    if (f == NULL || fgets(line, sizeof line, f) == NULL) {
        if (f != NULL) fclose(f);
        return -1;
    }
    while (fgets(line, sizeof line, f) != NULL) {
        QpCubeRow r;
        char frame[32];
        memset(&r, 0, sizeof r);
        if (sscanf(line, "%23[^\t]\t%zu\t%zu\t%31[^\t]\t%ld\t%ld\t%ld",
                   r.id, &r.nv, &r.nf, frame, &r.oz, &r.oy, &r.ox) != 7)
            continue;
        if (t.n == cap) {
            size_t next = cap == 0 ? 128 : cap * 2;
            QpCubeRow *p = (QpCubeRow *)realloc(t.row, next * sizeof *p);
            if (p == NULL) {
                fclose(f); qp_cube_table_dispose(&t); return -1;
            }
            t.row = p; cap = next;
        }
        r.v0 = t.nv; r.f0 = t.nf;
        if (SIZE_MAX - t.nv < r.nv || SIZE_MAX - t.nf < r.nf) {
            fclose(f); qp_cube_table_dispose(&t); return -1;
        }
        t.nv += r.nv; t.nf += r.nf;
        t.row[t.n++] = r;
    }
    fclose(f);
    if (t.n == 0) { qp_cube_table_dispose(&t); return -1; }
    *out = t;
    return 0;
}

static const QpCubeRow *qp_find_cube(const QpCubeTable *t, const char *id)
{
    size_t i = 0;
    for (i = 0; i < t->n; i++)
        if (strcmp(t->row[i].id, id) == 0) return &t->row[i];
    return NULL;
}

static void *qp_read_raw_malloc(const char *path, size_t element_size,
                                 size_t count, int required)
{
    FILE *f = NULL;
    void *data = NULL;
    int extra = 0;
    if (element_size == 0 || count > SIZE_MAX / element_size) return NULL;
    f = fopen(path, "rb");
    if (f == NULL) return NULL;
    data = malloc((count > 0 ? count : 1) * element_size);
    if (data == NULL || fread(data, element_size, count, f) != count) {
        free(data); fclose(f); return NULL;
    }
    extra = fgetc(f);
    if (extra != EOF || ferror(f) || fclose(f) != 0) {
        free(data); return NULL;
    }
    (void)required;
    return data;
}

static int qp_write_raw_exact(const char *path, const void *data,
                              size_t element_size, size_t count)
{
    FILE *f = NULL;
    if (data == NULL || element_size == 0 || count > SIZE_MAX / element_size)
        return -1;
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    if (fwrite(data, element_size, count, f) != count || ferror(f)) {
        fclose(f); return -1;
    }
    return fclose(f) == 0 ? 0 : -1;
}

static void qp_block_artifact_dispose(QpBlockArtifact *a)
{
    if (a == NULL) return;
    MeshBin_dispose(&a->cert);
    qp_cube_table_dispose(&a->cubes);
    free(a->winding); free(a->confidence); free(a->field); free(a->jump);
    free(a->relation); free(a->material); free(a->mesh_component);
    memset(a, 0, sizeof *a);
}

static int qp_load_block_artifact(const QpCanonicalBlock *b,
                                  QpBlockArtifact *out)
{
    char path[QP_MAX_PATH];
    QpBlockArtifact a;
    size_t nv = 0;
    memset(&a, 0, sizeof a);
    snprintf(path, sizeof path, "%s/stage2_unwrap/unwrap_winding.vmesh", b->root);
    if (MeshBin_read_malloc(path, &a.cert) != 0 || a.cert.uv == NULL) goto fail;
    nv = a.cert.nv;
    snprintf(path, sizeof path, "%s/stage1_mesh/mesh_cubes.tsv", b->root);
    if (qp_read_cube_table(path, &a.cubes) != 0 || a.cubes.nv != nv ||
        a.cubes.nf != a.cert.nf)
        goto fail;
#define QP_CB_READ(member, suffix, type, required_value) do {                 \
    snprintf(path, sizeof path, "%s/stage2_unwrap/unwrap_winding" suffix,   \
             b->root);                                                        \
    a.member = (type *)qp_read_raw_malloc(path, sizeof(type), nv,             \
                                           required_value);                    \
    if ((required_value) && a.member == NULL) goto fail;                       \
} while (0)
    QP_CB_READ(winding, "_index.f32", float, 1);
    QP_CB_READ(confidence, "_confidence.f32", float, 1);
    QP_CB_READ(material, "_material_identity.i32", int32_t, 1);
    QP_CB_READ(relation, "_relation_island.i32", int32_t, 1);
    QP_CB_READ(mesh_component, "_mesh_component.i32", int32_t, 1);
    QP_CB_READ(field, "_field.f32", float, 0);
    QP_CB_READ(jump, "_jump.f32", float, 0);
#undef QP_CB_READ
    *out = a;
    return 0;
fail:
    qp_block_artifact_dispose(&a);
    return -1;
}

static long qp_floor_to(long value, long step)
{
    long q = value / step;
    long r = value % step;
    if (r < 0) q--;
    return q * step;
}

static int qp_cb_owner_axis(long origin, long anchor, long chunk,
                            int size, int stride)
{
    long delta = origin - anchor;
    long k = delta / chunk;
    if (delta % chunk != 0) return INT_MIN;
    if (k >= 0 && k < size) return 0;
    if (k >= size)
        return (int)((k - (size - 1) + stride - 1) / stride);
    return -(int)((-k + stride - 1) / stride);
}

static void qp_cb_parent_coords(int iz, int iy, int ix,
                                int *pz, int *py, int *px)
{
    *pz = iz; *py = iy; *px = ix;
    if (ix != 0) *px += ix > 0 ? -1 : 1;
    else if (iy != 0) *py += iy > 0 ? -1 : 1;
    else if (iz != 0) *pz += iz > 0 ? -1 : 1;
}

static int qp_cb_distance(const QpCanonicalBlock *b)
{
    return abs(b->iz) + abs(b->iy) + abs(b->ix);
}

static int qp_cb_compare(const void *pa, const void *pb)
{
    const QpCanonicalBlock *a = (const QpCanonicalBlock *)pa;
    const QpCanonicalBlock *b = (const QpCanonicalBlock *)pb;
    int da = qp_cb_distance(a), db = qp_cb_distance(b);
    if (da != db) return da < db ? -1 : 1;
    if (a->iz != b->iz) return a->iz < b->iz ? -1 : 1;
    if (a->iy != b->iy) return a->iy < b->iy ? -1 : 1;
    return a->ix == b->ix ? 0 : (a->ix < b->ix ? -1 : 1);
}

static int qp_cb_find(const QpCanonicalBlock *blocks, size_t n,
                      int iz, int iy, int ix)
{
    size_t i = 0;
    for (i = 0; i < n; i++)
        if (blocks[i].iz == iz && blocks[i].iy == iy && blocks[i].ix == ix)
            return (int)i;
    return -1;
}

static int qp_cb_stable_label(int iz, int iy, int ix, size_t local,
                              int32_t *out)
{
    uint32_t key = 0;
    if (out == NULL || iz < QP_CB_LABEL_Z_MIN || iz > QP_CB_LABEL_Z_MAX ||
        iy < QP_CB_LABEL_Y_MIN || iy > QP_CB_LABEL_Y_MAX ||
        ix < QP_CB_LABEL_X_MIN || ix > QP_CB_LABEL_X_MAX ||
        local > (size_t)QP_CB_LABEL_LOCAL_MAX)
        return -1;
    key = ((uint32_t)(iz - QP_CB_LABEL_Z_MIN) << 24) |
          ((uint32_t)(iy - QP_CB_LABEL_Y_MIN) << 18) |
          ((uint32_t)(ix - QP_CB_LABEL_X_MIN) << 12) |
          (uint32_t)local;
    if (key > (uint32_t)INT32_MAX) return -1;
    *out = (int32_t)key;
    return 0;
}

static int qp_cb_add_with_ancestors(QpCanonicalBlock **blocks, size_t *n,
                                    size_t *cap, int iz, int iy, int ix)
{
    for (;;) {
        QpCanonicalBlock *p = NULL;
        int pz = 0, py = 0, px = 0;
        if (qp_cb_find(*blocks, *n, iz, iy, ix) < 0) {
            if (*n == *cap) {
                size_t next = *cap == 0 ? 16 : *cap * 2;
                p = (QpCanonicalBlock *)realloc(*blocks, next * sizeof *p);
                if (p == NULL) return -1;
                *blocks = p; *cap = next;
            }
            p = &(*blocks)[(*n)++];
            memset(p, 0, sizeof *p);
            p->iz = iz; p->iy = iy; p->ix = ix; p->parent = -1;
        }
        if (iz == 0 && iy == 0 && ix == 0) return 0;
        qp_cb_parent_coords(iz, iy, ix, &pz, &py, &px);
        iz = pz; iy = py; ix = px;
    }
}

static int qp_canonical_block_selftest(void)
{
    int fails = 0, pz = 0, py = 0, px = 0;
    int32_t a0 = -1, a1 = -1, b0 = -1;
    long chunk = 128, a = 4352;
    if (qp_cb_owner_axis(a, a, chunk, 4, 3) != 0) fails++;
    if (qp_cb_owner_axis(a + 3 * chunk, a, chunk, 4, 3) != 0) fails++;
    if (qp_cb_owner_axis(a + 4 * chunk, a, chunk, 4, 3) != 1) fails++;
    if (qp_cb_owner_axis(a + 6 * chunk, a, chunk, 4, 3) != 1) fails++;
    if (qp_cb_owner_axis(a + 7 * chunk, a, chunk, 4, 3) != 2) fails++;
    if (qp_cb_owner_axis(a - chunk, a, chunk, 5, 4) != -1) fails++;
    if (qp_cb_owner_axis(a - 4 * chunk, a, chunk, 5, 4) != -1) fails++;
    if (qp_cb_owner_axis(a - 5 * chunk, a, chunk, 5, 4) != -2) fails++;
    qp_cb_parent_coords(2, -1, 1, &pz, &py, &px);
    if (pz != 2 || py != -1 || px != 0) fails++;
    qp_cb_parent_coords(2, -1, 0, &pz, &py, &px);
    if (pz != 2 || py != 0 || px != 0) fails++;
    qp_cb_parent_coords(2, 0, 0, &pz, &py, &px);
    if (pz != 1 || py != 0 || px != 0) fails++;
    if (qp_cb_stable_label(0, 0, 0, 0, &a0) != 0 ||
        qp_cb_stable_label(0, 0, 0, 1, &a1) != 0 ||
        qp_cb_stable_label(1, 0, 0, 0, &b0) != 0 ||
        a0 == a1 || a0 == b0 || a1 == b0)
        fails++;
    if (qp_cb_stable_label(QP_CB_LABEL_Z_MIN - 1, 0, 0, 0, &a0) == 0 ||
        qp_cb_stable_label(0, 0, 0,
                           (size_t)QP_CB_LABEL_LOCAL_MAX + 1, &a0) == 0)
        fails++;
    fprintf(stderr, "[selftest] quadribbon canonical_blocks %s (%d failures)\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}

typedef struct {
    size_t child_v, parent_v;
    int32_t child_relation;
    double dq, du, dv;
} QpOverlapSample;

typedef struct {
    int32_t child, parent;
} QpLabelPair;

static int qp_cmp_i32(const void *pa, const void *pb)
{
    int32_t a = *(const int32_t *)pa, b = *(const int32_t *)pb;
    return a == b ? 0 : (a < b ? -1 : 1);
}

static int qp_cmp_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a == b ? 0 : (a < b ? -1 : 1);
}

static int qp_cmp_label_pair(const void *pa, const void *pb)
{
    const QpLabelPair *a = (const QpLabelPair *)pa;
    const QpLabelPair *b = (const QpLabelPair *)pb;
    if (a->child != b->child) return a->child < b->child ? -1 : 1;
    return a->parent == b->parent ? 0 : (a->parent < b->parent ? -1 : 1);
}

static size_t qp_label_count(const int32_t *label, size_t n)
{
    int32_t hi = -1;
    size_t i = 0;
    for (i = 0; i < n; i++) if (label[i] > hi) hi = label[i];
    return hi < 0 ? 0 : (size_t)hi + 1;
}

static int32_t qp_global_label(const int32_t *map, size_t n, int32_t local)
{
    return local >= 0 && (size_t)local < n ? map[local] : -1;
}

/* Integer mode with the same contract as the Python registration gate:
 * candidates are represented rounded integers and an inlier is within a
 * quarter turn.  The return value is the winning inlier count. */
static size_t qp_turn_mode(const QpOverlapSample *sample, const size_t *order,
                           size_t n, int32_t *out_k, double *out_agreement)
{
    int32_t *v = NULL;
    size_t used = 0, finite = 0, best = 0, i = 0;
    int32_t best_k = 0;
    if (n == 0) { *out_k = 0; *out_agreement = 0.0; return 0; }
    v = (int32_t *)malloc(n * sizeof *v);
    if (v == NULL) { *out_k = 0; *out_agreement = 0.0; return 0; }
    for (i = 0; i < n; i++) {
        double d = sample[order != NULL ? order[i] : i].dq;
        double k = nearbyint(d);
        if (!isfinite(d) || k < (double)INT32_MIN || k > (double)INT32_MAX)
            continue;
        finite++;
        if (fabs(d - k) < 0.25) v[used++] = (int32_t)k;
    }
    qsort(v, used, sizeof *v, qp_cmp_i32);
    for (i = 0; i < used; ) {
        size_t j = i + 1;
        while (j < used && v[j] == v[i]) j++;
        if (j - i > best || (j - i == best && v[i] < best_k)) {
            best = j - i; best_k = v[i];
        }
        i = j;
    }
    free(v);
    *out_k = best_k;
    *out_agreement = finite > 0 ? (double)best / (double)finite : 0.0;
    return best;
}

static double qp_sample_median(const QpOverlapSample *sample,
                               const size_t *order, size_t n, int which,
                               int32_t turn_k, int require_turn_inlier,
                               size_t *out_n)
{
    double *v = NULL;
    size_t used = 0, i = 0;
    double answer = 0.0;
    v = (double *)malloc((n > 0 ? n : 1) * sizeof *v);
    if (v == NULL) { if (out_n != NULL) *out_n = 0; return 0.0; }
    for (i = 0; i < n; i++) {
        const QpOverlapSample *s = &sample[order != NULL ? order[i] : i];
        double value = which == 0 ? s->du : s->dv;
        if (!isfinite(value)) continue;
        if (require_turn_inlier &&
            (!isfinite(s->dq) || fabs(s->dq - (double)turn_k) >= 0.25))
            continue;
        v[used++] = value;
    }
    qsort(v, used, sizeof *v, qp_cmp_double);
    if (used > 0)
        answer = used & 1 ? v[used / 2]
                          : 0.5 * (v[used / 2 - 1] + v[used / 2]);
    free(v);
    if (out_n != NULL) *out_n = used;
    return answer;
}

static int qp_project_label_map(const int32_t *child_label,
                                const int32_t *parent_label,
                                const int32_t *parent_map,
                                size_t parent_map_n,
                                const QpOverlapSample *sample, size_t nsample,
                                size_t nlocal, int32_t **out_map,
                                const QpCanonicalBlock *block,
                                int require_unanimous)
{
    QpLabelPair *pair = NULL;
    int32_t *map = NULL;
    size_t npair = 0, i = 0;
    map = (int32_t *)malloc((nlocal > 0 ? nlocal : 1) * sizeof *map);
    pair = (QpLabelPair *)malloc((nsample > 0 ? nsample : 1) * sizeof *pair);
    if (map == NULL || pair == NULL) { free(map); free(pair); return -1; }
    for (i = 0; i < nlocal; i++) map[i] = -1;
    for (i = 0; i < nsample; i++) {
        int32_t c = child_label[sample[i].child_v];
        int32_t p = parent_label[sample[i].parent_v];
        int32_t g = qp_global_label(parent_map, parent_map_n, p);
        if (c < 0 || (size_t)c >= nlocal || g < 0) continue;
        pair[npair].child = c; pair[npair].parent = g; npair++;
    }
    qsort(pair, npair, sizeof *pair, qp_cmp_label_pair);
    for (i = 0; i < npair; ) {
        size_t child_end = i + 1, j = i, best = 0;
        int32_t best_parent = pair[i].parent;
        while (child_end < npair && pair[child_end].child == pair[i].child)
            child_end++;
        while (j < child_end) {
            size_t k = j + 1;
            while (k < child_end && pair[k].parent == pair[j].parent) k++;
            if (k - j > best ||
                (k - j == best && pair[j].parent < best_parent)) {
                best = k - j; best_parent = pair[j].parent;
            }
            j = k;
        }
        if (require_unanimous && best != child_end - i) {
            qp_logf("[canonical blocks] %s local label %d intersects "
                    "multiple immutable parent identities (%zu/%zu agree); "
                    "refusing a majority-map fusion\n", block->id,
                    pair[i].child, best, child_end - i);
            free(pair); free(map); return -1;
        }
        map[(size_t)pair[i].child] = best_parent;
        i = child_end;
    }
    for (i = 0; i < nlocal; i++) {
        if (map[i] >= 0) continue;
        if (qp_cb_stable_label(block->iz, block->iy, block->ix, i,
                               &map[i]) != 0) {
            qp_logf("[canonical blocks] %s has an unrepresentable stable "
                    "label (block [%d,%d,%d], local %zu); refusing a "
                    "crop-dependent fallback\n", block->id, block->iz,
                    block->iy, block->ix, i);
            free(pair); free(map); return -1;
        }
    }
    free(pair);
    *out_map = map;
    return 0;
}

static int qp_root_label_map(const QpCanonicalBlock *block, size_t n,
                             int32_t **out)
{
    int32_t *map = (int32_t *)malloc((n > 0 ? n : 1) * sizeof *map);
    size_t i = 0;
    if (map == NULL) return -1;
    for (i = 0; i < n; i++) {
        if (qp_cb_stable_label(block->iz, block->iy, block->ix, i,
                               &map[i]) != 0) {
            qp_logf("[canonical blocks] root has %zu labels, beyond the "
                    "stable-label contract (maximum %d)\n", n,
                    QP_CB_LABEL_LOCAL_MAX + 1);
            free(map); return -1;
        }
    }
    *out = map;
    return 0;
}

static int qp_collect_parent_overlap(const QpBlockArtifact *child,
                                     const QpBlockArtifact *parent,
                                     const QpCanonicalBlock *parent_block,
                                     QpOverlapSample **out_sample,
                                     size_t *out_n)
{
    QpOverlapSample *sample = NULL;
    size_t n = 0, at = 0, i = 0;
    (void)parent_block;
    for (i = 0; i < child->cubes.n; i++) {
        const QpCubeRow *cr = &child->cubes.row[i];
        const QpCubeRow *pr = qp_find_cube(&parent->cubes, cr->id);
        if (pr != NULL) {
            if (pr->nv != cr->nv || pr->nf != cr->nf) return -1;
            if (SIZE_MAX - n < cr->nv) return -1;
            n += cr->nv;
        }
    }
    sample = (QpOverlapSample *)malloc((n > 0 ? n : 1) * sizeof *sample);
    if (sample == NULL) return -1;
    for (i = 0; i < child->cubes.n; i++) {
        const QpCubeRow *cr = &child->cubes.row[i];
        const QpCubeRow *pr = qp_find_cube(&parent->cubes, cr->id);
        size_t k = 0;
        if (pr == NULL) continue;
        if (memcmp(child->cert.verts + cr->v0 * 3,
                   parent->cert.verts + pr->v0 * 3,
                   cr->nv * 3 * sizeof(float)) != 0) {
            qp_logf("[canonical blocks] shared cube %s changed geometry; "
                    "registration refused\n", cr->id);
            free(sample); return -1;
        }
        for (k = 0; k < cr->nv; k++) {
            size_t cv = cr->v0 + k, pv = pr->v0 + k;
            sample[at].child_v = cv;
            sample[at].parent_v = pv;
            sample[at].child_relation = child->relation[cv];
            sample[at].dq = (double)parent->winding[pv] -
                            (double)child->winding[cv];
            sample[at].du = (double)parent->cert.uv[pv * 2] -
                            (double)child->cert.uv[cv * 2];
            sample[at].dv = (double)parent->cert.uv[pv * 2 + 1] -
                            (double)child->cert.uv[cv * 2 + 1];
            at++;
        }
    }
    *out_sample = sample; *out_n = at;
    return 0;
}

static int qp_solve_block_alignment(QpCanonicalBlock *block,
                                    const QpBlockArtifact *child,
                                    const QpCanonicalBlock *parent_block,
                                    const QpBlockArtifact *parent)
{
    QpOverlapSample *sample = NULL;
    size_t nsample = 0, i = 0;
    block->nrelation = qp_label_count(child->relation, child->cert.nv);
    block->nmaterial = qp_label_count(child->material, child->cert.nv);
    block->nmesh = qp_label_count(child->mesh_component, child->cert.nv);
    block->turn_offset = (int32_t *)calloc(
        block->nrelation > 0 ? block->nrelation : 1, sizeof(int32_t));
    block->u_offset = (double *)calloc(
        block->nrelation > 0 ? block->nrelation : 1, sizeof(double));
    if (block->turn_offset == NULL || block->u_offset == NULL) return -1;

    if (parent_block == NULL || parent == NULL) {
        if (qp_root_label_map(block, block->nrelation,
                              &block->relation_global) != 0 ||
            qp_root_label_map(block, block->nmaterial,
                              &block->material_global) != 0 ||
            qp_root_label_map(block, block->nmesh,
                              &block->mesh_global) != 0)
            return -1;
        block->v_offset = 0.0;
        block->turn_agreement = 1.0;
        block->ready = 1;
        return 0;
    }
    if (qp_collect_parent_overlap(child, parent, parent_block,
                                  &sample, &nsample) != 0)
        return -1;
    block->shared_vertices = nsample;
    if (nsample == 0) goto fail;
    for (i = 0; i < nsample; i++) {
        if (sample[i].dq != 0.0 || sample[i].du != 0.0 ||
            sample[i].dv != 0.0) {
            qp_logf("[canonical blocks] %s <- %s overlap is not exact at "
                    "sample %zu (dq=%+.9g dU=%+.9g dV=%+.9g); refusing "
                    "post-hoc alignment\n", block->id, parent_block->id, i,
                    sample[i].dq, sample[i].du, sample[i].dv);
            goto fail;
        }
    }
    block->v_offset = 0.0;
    block->turn_agreement = 1.0;
    if (qp_project_label_map(child->relation, parent->relation,
                             parent_block->relation_global,
                             parent_block->nrelation, sample, nsample,
                             block->nrelation, &block->relation_global,
                             block, 0) != 0 ||
        qp_project_label_map(child->material, parent->material,
                             parent_block->material_global,
                             parent_block->nmaterial, sample, nsample,
                             block->nmaterial, &block->material_global,
                             block, 1) != 0 ||
        qp_project_label_map(child->mesh_component, parent->mesh_component,
                             parent_block->mesh_global,
                             parent_block->nmesh, sample, nsample,
                             block->nmesh, &block->mesh_global,
                             block, 0) != 0)
        goto fail;
    qp_logf("[canonical blocks] %s <- %s: %zu shared vertices are "
            "bit-identical in q/U/V; inherited material identities are "
            "unambiguous\n", block->id, parent_block->id, nsample);
    block->ready = 1;
    free(sample);
    return 0;
fail:
    free(sample);
    return -1;
}

typedef struct {
    long chunk, anchor_z, anchor_y, anchor_x;
} QpCanonicalGrid;

static void qp_cb_name(char *out, size_t cap, int iz, int iy, int ix)
{
    snprintf(out, cap, "z%c%03d_y%c%03d_x%c%03d",
             iz < 0 ? 'm' : 'p', abs(iz),
             iy < 0 ? 'm' : 'p', abs(iy),
             ix < 0 ? 'm' : 'p', abs(ix));
}

static int qp_build_canonical_plan(const QpConfig *cfg,
                                   const QpCubeTable *target,
                                   const char *blocks_root,
                                   QpCanonicalGrid *grid,
                                   QpCanonicalBlock **out_blocks,
                                   size_t *out_n)
{
    QpCanonicalBlock *blocks = NULL;
    size_t n = 0, cap = 0, i = 0;
    long chunk = cfg->raw_chunk > 0 ? cfg->raw_chunk : 128;
    if (target == NULL || target->n == 0 || chunk <= 0) return -1;
    grid->chunk = chunk;
    if (cfg->canonical_anchor_armed) {
        grid->anchor_z = cfg->canonical_anchor[0];
        grid->anchor_y = cfg->canonical_anchor[1];
        grid->anchor_x = cfg->canonical_anchor[2];
        if (grid->anchor_z % chunk != 0 || grid->anchor_y % chunk != 0 ||
            grid->anchor_x % chunk != 0) {
            qp_logf("[canonical blocks] geometry.canonical_block_anchor_zyx "
                    "must be aligned to raw.chunk_size=%ld\n", chunk);
            return -1;
        }
    } else if (cfg->axis_warp.physical) {
        const AxisWarp *w=&cfg->axis_warp;
        const AxisWarpPhysical *c=&w->physical_config;
        grid->anchor_z=qp_floor_to((long)floor((w->z[0]-c->source_origin_um_zyx[0])/c->source_voxel_um_zyx[0]),chunk);
        grid->anchor_y=qp_floor_to((long)floor((w->y[0]-c->source_origin_um_zyx[1])/c->source_voxel_um_zyx[1]),chunk)-2*chunk;
        grid->anchor_x=qp_floor_to((long)floor((w->x[0]-c->source_origin_um_zyx[2])/c->source_voxel_um_zyx[2]),chunk)-2*chunk;
    } else if (cfg->axis_warp_armed && cfg->axis_warp.n > 0) {
        grid->anchor_z =
            qp_floor_to((long)floor(cfg->axis_warp.z[0]), chunk);
        grid->anchor_y =
            qp_floor_to((long)floor(cfg->axis_y), chunk) - 2 * chunk;
        grid->anchor_x =
            qp_floor_to((long)floor(cfg->axis_x), chunk) - 2 * chunk;
    } else {
        qp_logf("[canonical blocks] no crop-independent origin: configure "
                "geometry.canonical_block_anchor_zyx or geometry.axis_table; "
                "refusing to derive alignment from the requested crop\n");
        return -1;
    }
    for (i = 0; i < target->n; i++) {
        const QpCubeRow *r = &target->row[i];
        int iz = qp_cb_owner_axis(r->oz, grid->anchor_z, chunk,
                                  QP_CB_Z_SIZE, QP_CB_Z_STRIDE);
        int iy = qp_cb_owner_axis(r->oy, grid->anchor_y, chunk,
                                  QP_CB_Y_SIZE, QP_CB_Y_STRIDE);
        int ix = qp_cb_owner_axis(r->ox, grid->anchor_x, chunk,
                                  QP_CB_X_SIZE, QP_CB_X_STRIDE);
        if (iz == INT_MIN || iy == INT_MIN || ix == INT_MIN ||
            qp_cb_add_with_ancestors(&blocks, &n, &cap, iz, iy, ix) != 0) {
            free(blocks); return -1;
        }
    }
    qsort(blocks, n, sizeof *blocks, qp_cb_compare);
    for (i = 0; i < n; i++) {
        QpCanonicalBlock *b = &blocks[i];
        int pz = 0, py = 0, px = 0;
        long z0 = grid->anchor_z +
                  (long)b->iz * QP_CB_Z_STRIDE * chunk;
        long y0 = grid->anchor_y +
                  (long)b->iy * QP_CB_Y_STRIDE * chunk;
        long x0 = grid->anchor_x +
                  (long)b->ix * QP_CB_X_STRIDE * chunk;
        b->bounds[0] = z0;
        b->bounds[1] = z0 + (QP_CB_Z_SIZE - 1) * chunk;
        b->bounds[2] = y0;
        b->bounds[3] = y0 + (QP_CB_Y_SIZE - 1) * chunk;
        b->bounds[4] = x0;
        b->bounds[5] = x0 + (QP_CB_X_SIZE - 1) * chunk;
        qp_cb_name(b->id, sizeof b->id, b->iz, b->iy, b->ix);
        snprintf(b->root, sizeof b->root, "%s/%s", blocks_root, b->id);
        if (b->iz == 0 && b->iy == 0 && b->ix == 0) {
            b->parent = -1;
        } else {
            qp_cb_parent_coords(b->iz, b->iy, b->ix, &pz, &py, &px);
            b->parent = qp_cb_find(blocks, n, pz, py, px);
            if (b->parent < 0 || (size_t)b->parent >= i) {
                free(blocks); return -1;
            }
        }
    }
    *out_blocks = blocks; *out_n = n;
    return 0;
}

/* Materialize exact parent evidence in the child's input-vertex order.  Cube
 * rows are immutable concat units, so matching cube id + byte-identical
 * geometry gives a one-to-one correspondence without a spatial search or a
 * fitted transform. */
static int qp_write_projective_boundary(
    const char *child_mesh_path, const char *child_table_path,
    const QpCanonicalBlock *parent_block,
    const QpBlockArtifact *parent,
    const char *winding_path, const char *u_path, const char *v_path,
    const char *material_path, size_t *out_shared)
{
    MeshBinData child_mesh;
    QpCubeTable child_cubes;
    float *winding = NULL, *u = NULL, *v = NULL;
    int32_t *material = NULL;
    size_t shared = 0;
    int rc = -1;
    memset(&child_mesh, 0, sizeof child_mesh);
    memset(&child_cubes, 0, sizeof child_cubes);
    if (parent_block == NULL || parent == NULL || !parent_block->ready ||
        MeshBin_read_malloc(child_mesh_path, &child_mesh) != 0 ||
        qp_read_cube_table(child_table_path, &child_cubes) != 0 ||
        child_cubes.nv != child_mesh.nv)
        goto done;
    winding = (float *)malloc((child_mesh.nv > 0 ? child_mesh.nv : 1) *
                              sizeof *winding);
    u = (float *)malloc((child_mesh.nv > 0 ? child_mesh.nv : 1) * sizeof *u);
    v = (float *)malloc((child_mesh.nv > 0 ? child_mesh.nv : 1) * sizeof *v);
    material = (int32_t *)malloc(
        (child_mesh.nv > 0 ? child_mesh.nv : 1) * sizeof *material);
    if (winding == NULL || u == NULL || v == NULL || material == NULL)
        goto done;
    for (size_t i = 0; i < child_mesh.nv; i++) {
        winding[i] = NAN; u[i] = NAN; v[i] = NAN; material[i] = -1;
    }
    for (size_t i = 0; i < child_cubes.n; i++) {
        const QpCubeRow *cr = &child_cubes.row[i];
        const QpCubeRow *pr = qp_find_cube(&parent->cubes, cr->id);
        if (pr == NULL) continue;
        if (pr->nv != cr->nv || pr->nf != cr->nf ||
            memcmp(child_mesh.verts + cr->v0 * 3,
                   parent->cert.verts + pr->v0 * 3,
                   cr->nv * 3 * sizeof(float)) != 0) {
            qp_logf("[canonical blocks] parent boundary cube %s changed "
                    "topology/geometry; refusing correspondence\n", cr->id);
            goto done;
        }
        for (size_t k = 0; k < cr->nv; k++) {
            size_t cv = cr->v0 + k, pv = pr->v0 + k;
            int32_t local_material = parent->material[pv];
            int32_t global_material = qp_global_label(
                parent_block->material_global, parent_block->nmaterial,
                local_material);
            if (!isfinite((double)parent->winding[pv]) ||
                !isfinite((double)parent->cert.uv[pv * 2]) ||
                !isfinite((double)parent->cert.uv[pv * 2 + 1]) ||
                global_material < 0)
                goto done;
            winding[cv] = parent->winding[pv];
            u[cv] = parent->cert.uv[pv * 2];
            v[cv] = parent->cert.uv[pv * 2 + 1];
            material[cv] = global_material;
            shared++;
        }
    }
    if (shared == 0) {
        qp_logf("[canonical blocks] child has no exact vertices shared with "
                "parent %s\n", parent_block->id);
        goto done;
    }
    if (qp_write_raw_exact(winding_path, winding, sizeof *winding,
                           child_mesh.nv) != 0 ||
        qp_write_raw_exact(u_path, u, sizeof *u, child_mesh.nv) != 0 ||
        qp_write_raw_exact(v_path, v, sizeof *v, child_mesh.nv) != 0 ||
        qp_write_raw_exact(material_path, material, sizeof *material,
                           child_mesh.nv) != 0)
        goto done;
    if (out_shared != NULL) *out_shared = shared;
    rc = 0;
done:
    MeshBin_dispose(&child_mesh);
    qp_cube_table_dispose(&child_cubes);
    free(winding); free(u); free(v); free(material);
    return rc;
}

static int qp_write_block_boundary_contract(
    const char *path, const QpCanonicalBlock *block,
    const QpCanonicalBlock *parent, size_t shared)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f,
            "{ \"schema\": \"vesuvius-canonical-boundary-v3\", "
            "\"block\": \"%s\", \"parent\": %s, "
            "\"shared_vertices\": %zu, "
            "\"rule\": \"component-locked q; inherited lineage; exact "
            "parent U(q),V overlap\" }\n",
            block->id, parent != NULL ? "\"set\"" : "null", shared);
    return fclose(f) == 0 ? 0 : -1;
}

static int qp_run_canonical_block(const QpConfig *cfg, const char *input,
                                  QpCanonicalBlock *block,
                                  const QpCanonicalBlock *parent_block)
{
    char stage1[QP_MAX_PATH], mesh[QP_MAX_PATH], table[QP_MAX_PATH];
    char stage2[QP_MAX_PATH], cert[QP_MAX_PATH], logs[QP_MAX_PATH];
    char boundary_winding[QP_MAX_PATH], boundary_u[QP_MAX_PATH];
    char boundary_v[QP_MAX_PATH], boundary_material[QP_MAX_PATH];
    char boundary_contract[QP_MAX_PATH];
    char saved_log[QP_MAX_PATH];
    long saved_sg[6];
    int saved_armed = g_qp_subgrid_armed;
    Arena_T arena = NULL;
    int rc = -1, k = 0;
    memcpy(saved_sg, g_qp_sg, sizeof saved_sg);
    snprintf(saved_log, sizeof saved_log, "%s", qp_log_dir);
    qp_mkdir(block->root);
    qp_join(stage1, sizeof stage1, block->root, "stage1_mesh");
    qp_join(mesh, sizeof mesh, stage1, "mesh.vmesh");
    qp_join(table, sizeof table, stage1, "mesh_cubes.tsv");
    qp_join(stage2, sizeof stage2, block->root, "stage2_unwrap");
    qp_join(cert, sizeof cert, stage2, "unwrap_winding.vmesh");
    qp_join(boundary_winding, sizeof boundary_winding, stage2,
            "parent_boundary_winding.f32");
    qp_join(boundary_u, sizeof boundary_u, stage2, "parent_boundary_u.f32");
    qp_join(boundary_v, sizeof boundary_v, stage2, "parent_boundary_v.f32");
    qp_join(boundary_material, sizeof boundary_material, stage2,
            "parent_boundary_material.i32");
    qp_join(boundary_contract, sizeof boundary_contract, stage2,
            "canonical_boundary.json");
    qp_join(logs, sizeof logs, block->root, "logs");
    qp_mkdir(logs);
    snprintf(qp_log_dir, sizeof qp_log_dir, "%s", logs);
    if (!MeshBin_looks_complete(cert) ||
        !qp_file_has_text(boundary_contract,
                          "vesuvius-canonical-boundary-v3")) {
        g_qp_subgrid_armed = 1;
        for (k = 0; k < 6; k++) g_qp_sg[k] = block->bounds[k];
        arena = Arena_new();
        qp_logf("[canonical blocks] build %s z=%ld..%ld y=%ld..%ld "
                "x=%ld..%ld\n", block->id, block->bounds[0], block->bounds[1],
                block->bounds[2], block->bounds[3], block->bounds[4],
                block->bounds[5]);
        if (arena == NULL ||
            qp_stage_mesh(arena, cfg, input, 0, stage1, mesh) != 0)
            goto done;
        g_qp_subgrid_armed = saved_armed;
        memcpy(g_qp_sg, saved_sg, sizeof saved_sg);
        if (arena != NULL) { Arena_dispose(&arena); arena = NULL; }
        qp_mkdir(stage2);
        size_t shared = 0;
        QpBlockArtifact parent;
        memset(&parent, 0, sizeof parent);
        if (parent_block != NULL &&
            (qp_load_block_artifact(parent_block, &parent) != 0 ||
             qp_write_projective_boundary(
                 mesh, table, parent_block, &parent,
                 boundary_winding, boundary_u, boundary_v,
                 boundary_material, &shared) != 0)) {
            qp_block_artifact_dispose(&parent);
            goto done;
        }
        if (qp_spawn_scroll_ribbon(cfg, mesh, stage2, "unwrap", 0, 0,
                                   parent_block != NULL ? boundary_winding : NULL,
                                   parent_block != NULL ? boundary_u : NULL,
                                   parent_block != NULL ? boundary_v : NULL,
                                   parent_block != NULL ? boundary_material : NULL,
                                   "scroll_ribbon(canonical_block)") != 0 ||
            !MeshBin_looks_complete(cert) ||
            qp_write_block_boundary_contract(
                boundary_contract, block, parent_block, shared) != 0) {
            qp_block_artifact_dispose(&parent);
            goto done;
        }
        qp_block_artifact_dispose(&parent);
    } else {
        qp_logf("[canonical blocks] resume %s\n", block->id);
    }
    rc = 0;
done:
    g_qp_subgrid_armed = saved_armed;
    memcpy(g_qp_sg, saved_sg, sizeof saved_sg);
    snprintf(qp_log_dir, sizeof qp_log_dir, "%s", saved_log);
    if (arena != NULL) Arena_dispose(&arena);
    return rc;
}

static int qp_run_and_align_canonical_blocks(
    const QpConfig *cfg, const char *input,
    QpCanonicalBlock *blocks, size_t n)
{
    size_t i = 0;
    for (i = 0; i < n; i++) {
        QpBlockArtifact child, parent;
        QpCanonicalBlock *pb = NULL;
        memset(&child, 0, sizeof child); memset(&parent, 0, sizeof parent);
        if (blocks[i].parent >= 0) {
            pb = &blocks[blocks[i].parent];
            if (!pb->ready) return -1;
        }
        if (qp_run_canonical_block(cfg, input, &blocks[i], pb) != 0)
            return -1;
        if (qp_load_block_artifact(&blocks[i], &child) != 0) {
            qp_logf("[canonical blocks] cannot load %s certificate\n",
                    blocks[i].id);
            return -1;
        }
        if (pb != NULL) {
            if (!pb->ready || qp_load_block_artifact(pb, &parent) != 0) {
                qp_block_artifact_dispose(&child);
                return -1;
            }
        }
        if (qp_solve_block_alignment(&blocks[i], &child, pb,
                                     pb != NULL ? &parent : NULL) != 0) {
            qp_block_artifact_dispose(&parent);
            qp_block_artifact_dispose(&child);
            return -1;
        }
        qp_block_artifact_dispose(&parent);
        qp_block_artifact_dispose(&child);
    }
    return 0;
}

static int qp_cb_owner_index(const QpCanonicalGrid *grid,
                             const QpCanonicalBlock *blocks, size_t n,
                             const QpCubeRow *r)
{
    int iz = qp_cb_owner_axis(r->oz, grid->anchor_z, grid->chunk,
                              QP_CB_Z_SIZE, QP_CB_Z_STRIDE);
    int iy = qp_cb_owner_axis(r->oy, grid->anchor_y, grid->chunk,
                              QP_CB_Y_SIZE, QP_CB_Y_STRIDE);
    int ix = qp_cb_owner_axis(r->ox, grid->anchor_x, grid->chunk,
                              QP_CB_X_SIZE, QP_CB_X_STRIDE);
    return qp_cb_find(blocks, n, iz, iy, ix);
}

static int qp_write_canonical_ledger(const char *unwrap_dir,
                                     const QpCanonicalGrid *grid,
                                     const QpCanonicalBlock *blocks, size_t n,
                                     size_t nv, size_t nf)
{
    char path[QP_MAX_PATH];
    FILE *f = NULL;
    size_t i = 0;
    qp_join(path, sizeof path, unwrap_dir, "canonical_blocks.json");
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f,
            "{\n  \"schema\": \"vesuvius-canonical-block-certificate-v3\",\n"
            "  \"rule\": \"immutable parent; component-locked winding and "
            "lineage; exact recursive U(q),V overlap; reversible block-local "
            "stable labels\",\n"
            "  \"label_key\": \"z7:y6:x6:local12 signed-bias packed i31\",\n"
            "  \"chunk\": %ld, \"anchor_zyx\": [%ld, %ld, %ld],\n"
            "  \"block_size_zyx\": [%d, %d, %d], "
            "\"stride_zyx\": [%d, %d, %d],\n"
            "  \"vertices\": %zu, \"faces\": %zu,\n  \"blocks\": [\n",
            grid->chunk, grid->anchor_z, grid->anchor_y, grid->anchor_x,
            QP_CB_Z_SIZE, QP_CB_Y_SIZE, QP_CB_X_SIZE,
            QP_CB_Z_STRIDE, QP_CB_Y_STRIDE, QP_CB_X_STRIDE, nv, nf);
    for (i = 0; i < n; i++) {
        const QpCanonicalBlock *b = &blocks[i];
        fprintf(f,
                "    %s{ \"id\": \"%s\", \"index_zyx\": [%d,%d,%d], "
                "\"parent\": %s, \"bounds\": [%ld,%ld,%ld,%ld,%ld,%ld], "
                "\"shared_vertices\": %zu, \"turn_agreement\": %.9g, "
                "\"v_offset\": %.9g }",
                i > 0 ? ",\n" : "", b->id, b->iz, b->iy, b->ix,
                b->parent >= 0 ? "\"set\"" : "null",
                b->bounds[0], b->bounds[1], b->bounds[2], b->bounds[3],
                b->bounds[4], b->bounds[5], b->shared_vertices,
                b->turn_agreement, b->v_offset);
    }
    fprintf(f, "\n  ]\n}\n");
    return fclose(f) == 0 ? 0 : -1;
}

static int qp_merge_canonical_blocks(const char *mesh_vmesh,
                                     const char *mesh_table,
                                     const char *unwrap_dir,
                                     const char *cert_vmesh,
                                     const QpCanonicalGrid *grid,
                                     const QpCanonicalBlock *blocks, size_t n,
                                     const QpCubeTable *target)
{
    MeshBinData mesh;
    float *uv = NULL, *winding = NULL, *confidence = NULL;
    float *field = NULL, *jump = NULL;
    int32_t *material = NULL, *relation = NULL, *component = NULL;
    int32_t *owner_block = NULL;
    uint8_t *assigned = NULL;
    int optional_complete = 1, rc = -1;
    size_t bi = 0, assigned_rows = 0;
    char path[QP_MAX_PATH];
    memset(&mesh, 0, sizeof mesh);
    if (MeshBin_read_malloc(mesh_vmesh, &mesh) != 0 || mesh.nv != target->nv ||
        mesh.nf != target->nf)
        goto done;
    uv = (float *)malloc((mesh.nv > 0 ? mesh.nv : 1) * 2 * sizeof *uv);
    winding = (float *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *winding);
    confidence = (float *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *confidence);
    field = (float *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *field);
    jump = (float *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *jump);
    material = (int32_t *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *material);
    relation = (int32_t *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *relation);
    component = (int32_t *)malloc((mesh.nv > 0 ? mesh.nv : 1) * sizeof *component);
    owner_block = (int32_t *)malloc(
        (mesh.nv > 0 ? mesh.nv : 1) * 3 * sizeof *owner_block);
    assigned = (uint8_t *)calloc(target->n > 0 ? target->n : 1, 1);
    if (uv == NULL || winding == NULL || confidence == NULL || field == NULL ||
        jump == NULL || material == NULL || relation == NULL ||
        component == NULL || owner_block == NULL || assigned == NULL)
        goto done;

    for (bi = 0; bi < n; bi++) {
        QpBlockArtifact a;
        size_t ti = 0;
        memset(&a, 0, sizeof a);
        if (qp_load_block_artifact(&blocks[bi], &a) != 0) goto done;
        if (a.field == NULL || a.jump == NULL) optional_complete = 0;
        for (ti = 0; ti < target->n; ti++) {
            const QpCubeRow *tr = &target->row[ti];
            const QpCubeRow *br = NULL;
            size_t k = 0;
            int owner = qp_cb_owner_index(grid, blocks, n, tr);
            if (owner != (int)bi) continue;
            br = qp_find_cube(&a.cubes, tr->id);
            if (br == NULL || br->nv != tr->nv || br->nf != tr->nf ||
                memcmp(mesh.verts + tr->v0 * 3, a.cert.verts + br->v0 * 3,
                       tr->nv * 3 * sizeof(float)) != 0) {
                qp_logf("[canonical blocks] owner %s cannot reproduce cube %s; "
                        "merge refused\n", blocks[bi].id, tr->id);
                qp_block_artifact_dispose(&a); goto done;
            }
            for (k = 0; k < tr->nv; k++) {
                size_t tv = tr->v0 + k, bv = br->v0 + k;
                int32_t rel = a.relation[bv];
                int32_t koff = rel >= 0 && (size_t)rel < blocks[bi].nrelation
                             ? blocks[bi].turn_offset[rel] : 0;
                double uoff = rel >= 0 && (size_t)rel < blocks[bi].nrelation
                            ? blocks[bi].u_offset[rel] : 0.0;
                winding[tv] = koff == 0 ? a.winding[bv]
                             : (float)((double)a.winding[bv] + (double)koff);
                uv[tv * 2] = uoff == 0.0 ? a.cert.uv[bv * 2]
                             : (float)((double)a.cert.uv[bv * 2] + uoff);
                uv[tv * 2 + 1] = blocks[bi].v_offset == 0.0
                                 ? a.cert.uv[bv * 2 + 1]
                                 : (float)((double)a.cert.uv[bv * 2 + 1] +
                                           blocks[bi].v_offset);
                confidence[tv] = a.confidence[bv];
                if (a.field != NULL) field[tv] = a.field[bv];
                if (a.jump != NULL) jump[tv] = a.jump[bv];
                relation[tv] = qp_global_label(blocks[bi].relation_global,
                                                blocks[bi].nrelation, rel);
                material[tv] = qp_global_label(blocks[bi].material_global,
                                                blocks[bi].nmaterial,
                                                a.material[bv]);
                component[tv] = qp_global_label(blocks[bi].mesh_global,
                                                 blocks[bi].nmesh,
                                                 a.mesh_component[bv]);
                owner_block[tv * 3] = (int32_t)blocks[bi].iz;
                owner_block[tv * 3 + 1] = (int32_t)blocks[bi].iy;
                owner_block[tv * 3 + 2] = (int32_t)blocks[bi].ix;
            }
            assigned[ti] = 1; assigned_rows++;
        }
        qp_block_artifact_dispose(&a);
    }
    if (assigned_rows != target->n) {
        qp_logf("[canonical blocks] only %zu/%zu requested cubes received an "
                "owner; merge refused\n", assigned_rows, target->n);
        goto done;
    }
    if (MeshBin_write(cert_vmesh, mesh.verts, mesh.nv, mesh.faces, mesh.nf,
                      uv) != 0)
        goto done;
#define QP_CB_WRITE(suffix, data, type) do {                                  \
    snprintf(path, sizeof path, "%s/unwrap_winding" suffix, unwrap_dir);     \
    if (qp_write_raw_exact(path, data, sizeof(type), mesh.nv) != 0) goto done; \
} while (0)
    QP_CB_WRITE("_index.f32", winding, float);
    QP_CB_WRITE("_confidence.f32", confidence, float);
    QP_CB_WRITE("_material_identity.i32", material, int32_t);
    QP_CB_WRITE("_relation_island.i32", relation, int32_t);
    QP_CB_WRITE("_mesh_component.i32", component, int32_t);
    snprintf(path, sizeof path, "%s/unwrap_winding_owner_block.i32",
             unwrap_dir);
    if (qp_write_raw_exact(path, owner_block, sizeof(int32_t),
                           mesh.nv * 3) != 0)
        goto done;
    if (optional_complete) {
        QP_CB_WRITE("_field.f32", field, float);
        QP_CB_WRITE("_jump.f32", jump, float);
    }
#undef QP_CB_WRITE
    qp_join(path, sizeof path, unwrap_dir, "mesh_cubes.tsv");
    if (qp_copy_file(mesh_table, path) != 0 ||
        qp_write_canonical_ledger(unwrap_dir, grid, blocks, n, mesh.nv,
                                  mesh.nf) != 0)
        goto done;
    qp_join(path, sizeof path, unwrap_dir,
            "unwrap_winding_visualization.json");
    {
        FILE *f = fopen(path, "wb");
        if (f == NULL) goto done;
        fprintf(f,
                "{ \"schema\": \"vesuvius-winding-certificate-v4-canonical-blocks\", "
                "\"vertices\": %zu, \"faces\": %zu, \"blocks\": %zu, "
                "\"ownership\": \"root-and-parent-prefix immutable\", "
                "\"owner_sidecar\": \"unwrap_winding_owner_block.i32\" }\n",
                mesh.nv, mesh.nf, n);
        if (fclose(f) != 0) goto done;
    }
    qp_logf("[canonical blocks] merged %zu requested cubes from %zu immutable "
            "blocks -> %s\n", target->n, n, cert_vmesh);
    rc = 0;
done:
    MeshBin_dispose(&mesh);
    free(uv); free(winding); free(confidence); free(field); free(jump);
    free(material); free(relation); free(component); free(owner_block);
    free(assigned);
    return rc;
}

static void qp_canonical_blocks_dispose(QpCanonicalBlock *blocks, size_t n)
{
    size_t i = 0;
    for (i = 0; i < n; i++) {
        free(blocks[i].relation_global); free(blocks[i].material_global);
        free(blocks[i].mesh_global); free(blocks[i].turn_offset);
        free(blocks[i].u_offset);
    }
    free(blocks);
}

static int qp_stage_unwrap_canonical(const QpConfig *cfg, const char *input,
                                     const char *mesh_vmesh,
                                     const char *mesh_table,
                                     const char *unwrap_dir,
                                     const char *cert_vmesh)
{
    QpCubeTable target;
    QpCanonicalGrid grid;
    QpCanonicalBlock *blocks = NULL;
    size_t n = 0, i = 0;
    char blocks_root[QP_MAX_PATH];
    int rc = -1;
    memset(&target, 0, sizeof target); memset(&grid, 0, sizeof grid);
    if (g_qp_tile_armed || qp_read_cube_table(mesh_table, &target) != 0)
        goto done;
    qp_join(blocks_root, sizeof blocks_root, unwrap_dir, "canonical_blocks");
    qp_mkdir(blocks_root);
    if (qp_build_canonical_plan(cfg, &target, blocks_root, &grid, &blocks,
                                &n) != 0)
        goto done;
    qp_logf("[canonical blocks] plan: %zu immutable 4x5x5 blocks, anchor "
            "(%ld,%ld,%ld), requested cubes=%zu\n", n, grid.anchor_z,
            grid.anchor_y, grid.anchor_x, target.n);
    if (qp_run_and_align_canonical_blocks(cfg, input, blocks, n) != 0)
        goto done;
    if (qp_merge_canonical_blocks(mesh_vmesh, mesh_table, unwrap_dir,
                                  cert_vmesh, &grid, blocks, n, &target) != 0)
        goto done;
    rc = 0;
done:
    qp_cube_table_dispose(&target);
    qp_canonical_blocks_dispose(blocks, n);
    return rc;
}

/* ---- stage 4: fit the solid quad ribbon ----------------------------------- */

static int qp_stage_solid(const QpConfig *cfg, const char *in_vmesh,
                          const char *out_stem, const char *tag)
{
    Arena_T arena = Arena_new();
    QsConfig qc;
    QsReport rep;
    int rc = -1;
    if (arena == NULL) return -1;
    {
        char frame[QP_MAX_PATH],existing[QP_MAX_PATH];
        if (snprintf(frame,sizeof frame,"%s.coordinate_frame.json",out_stem)>=(int)sizeof frame ||
            snprintf(existing,sizeof existing,"%s.vmesh",out_stem)>=(int)sizeof existing ||
            qp_coordinate_identity(cfg,frame,existing)!=0) { Arena_dispose(&arena); return -1; }
    }
    memset(&qc, 0, sizeof qc);
    qc.axis_y = cfg->axis_y;
    qc.axis_x = cfg->axis_x;
    qc.pitch = cfg->wrap_spacing;
    qc.projective_local = g_qp_canonical_blocks;
    qc.projective_copy_only = g_qp_canonical_blocks;
    qc.lattice_du = QP_FIT_GRID_U;
    qc.lattice_dv = QP_FIT_SLICE_H;
    qc.local_fill_radius = QS_BRIDGE_MAX_CELLS;
    qc.run_min_observed = QS_RUN_MIN_OBSERVED;
    qc.good_column_frac = QS_GOOD_COLUMN_FRAC;
    qc.fill_max_u_gap = QS_FILL_MAX_U_GAP;
    qc.skip_dr_pitch = QS_SKIP_DR_PITCH;
    qc.run_merge_gap = QS_RUN_MERGE_GAP_COLS;
    qc.raw_zarr = cfg->raw_path;      /* CT probes: MIDLINE / fill snap / ridge track */
    qc.coordinate_warp=cfg->axis_warp.physical ? &cfg->axis_warp : NULL;
    if (cfg->axis_warp_armed && AxisWarp_valid(&cfg->axis_warp)) {
        /* the stage works in straightened coordinates; give it the curve so a CT
         * probe can go back to world (measured 2026-09-03: without this every
         * probe read the CT up to the axis drift away, ~200 vox on the 10x3x3) */
        if (!cfg->axis_warp.physical) {
            qc.axis_curve_z = cfg->axis_warp.z;
            qc.axis_curve_y = cfg->axis_warp.y;
            qc.axis_curve_x = cfg->axis_warp.x;
            qc.axis_curve_n = cfg->axis_warp.n;
        }
        qc.axis_ref_y = cfg->axis_warp.reference_y;
        qc.axis_ref_x = cfg->axis_warp.reference_x;
    }
    rc = QuadribbonSolidify_run(arena, &qc, in_vmesh, out_stem, &rep);
    Arena_dispose(&arena);
    if (rc != 0) {
        qp_logf("[%s] solidify FAILED on %s\n", tag, in_vmesh);
        return -1;
    }
    qp_logf("[%s] %s.vmesh: %zu verts (%zu observed, %zu adopted band-1, "
            "%zu filled), %zu faces; runs %d (%d built, %d retried, %d "
            "observed-only); bands %d; duplicates %zu (|dr| p50 %.2f pitch); "
            "band-1 rejected %zu; not filled unbridged/skip/order/chord/material/bbox/core/"
            "unref %zu/%zu/%zu/%zu/%zu/%zu/%zu/%zu; confetti observed cells %zu; faces skipped "
            "gate/material %zu/%zu; phase mismatch max %.3f rad; %.1fs\n",
            tag, out_stem, rep.nv_out, rep.observed, rep.band1_adopted,
            rep.filled, rep.nf_out, rep.n_runs, rep.runs_built,
            rep.runs_retried, rep.runs_observed_only, rep.n_bands,
            rep.duplicates, rep.dup_dr_pitch_p50, rep.band1_rejected,
            rep.dropped_unbridged, rep.dropped_skip, rep.dropped_order, rep.dropped_chord,
            rep.dropped_material, rep.dropped_bbox, rep.dropped_core,
            rep.dropped_unreferenced, rep.dropped_confetti, rep.faces_skipped_gate,
            rep.faces_skipped_material, rep.phase_mismatch_max, rep.t_total);
    return 0;
}

/* A projective stage-3 result is a semantic atlas.  Collision colors are only
 * storage coordinates: layer 0 is <stem>.vmesh, layer 1 keeps the historical
 * <stem>_extras.vmesh name, and later layers are explicitly numbered.  Never
 * concatenate these meshes or choose a winner between coincident U/V samples.
 */
typedef struct {
    size_t layer;
    char input_vmesh[QP_MAX_PATH];
    char output_stem[QP_MAX_PATH];
    char output_vmesh[QP_MAX_PATH];
} QpProjectiveAtlasChart;

typedef struct {
    QpProjectiveAtlasChart *charts;
    unsigned char *has_surface;
    size_t n;
    size_t declared_layers;
    char source_stats[QP_MAX_PATH];
} QpProjectiveAtlas;

static void qp_projective_atlas_dispose(QpProjectiveAtlas *atlas)
{
    if (atlas == NULL) return;
    free(atlas->has_surface);
    free(atlas->charts);
    memset(atlas, 0, sizeof *atlas);
}

static int qp_projective_layer_stem(char *out, size_t cap,
                                    const char *base, size_t layer)
{
    int n = 0;
    if (layer == 0)
        n = snprintf(out, cap, "%s", base);
    else if (layer == 1)
        n = snprintf(out, cap, "%s_extras", base);
    else
        n = snprintf(out, cap, "%s_extras_layer%02zu", base, layer);
    return n >= 0 && (size_t)n < cap ? 0 : -1;
}

/* Read the declared color count from the stage-3 stats, then inspect only the
 * filenames prescribed by that declaration.  A missing secondary file is
 * legal: ribbon_grid_io deliberately does not emit a VMESH for a point-only
 * color.  A present but malformed file is an error, not an omitted chart. */
static int qp_projective_atlas_discover(const char *fit_stem,
                                        const char *solid_stem,
                                        QpProjectiveAtlas *atlas)
{
    Arena_T arena = NULL;
    const JsonValue *root = NULL, *ribbon = NULL;
    const char *err = NULL;
    long declared = -1;
    char primary_vmesh[QP_MAX_PATH];
    int rc = -1;
    if (fit_stem == NULL || solid_stem == NULL || atlas == NULL) return -1;
    memset(atlas, 0, sizeof *atlas);
    if (snprintf(atlas->source_stats, sizeof atlas->source_stats,
                 "%s_stats.json", fit_stem) < 0 ||
        strlen(atlas->source_stats) + 1 >= sizeof atlas->source_stats)
        return -1;
    arena = Arena_new();
    if (arena == NULL) return -1;
    root = Json_parse_file(arena, atlas->source_stats, &err);
    ribbon = root != NULL ? Json_object_get(root, "ribbon") : NULL;
    declared = ribbon != NULL ? Json_member_long(ribbon, "grid_layers", -1)
                              : -1;
    if (declared < 1 || declared > 100000) {
        qp_logf("[stage 4 atlas] cannot read a valid ribbon.grid_layers from "
                "%s%s%s\n", atlas->source_stats, err != NULL ? ": " : "",
                err != NULL ? err : "");
        goto done;
    }
    atlas->declared_layers = (size_t)declared;
    if (atlas->declared_layers > SIZE_MAX / sizeof *atlas->charts) goto done;
    atlas->charts = (QpProjectiveAtlasChart *)calloc(
        atlas->declared_layers, sizeof *atlas->charts);
    atlas->has_surface = (unsigned char *)calloc(atlas->declared_layers, 1);
    if (atlas->charts == NULL || atlas->has_surface == NULL) goto done;
    if (qp_projective_layer_stem(primary_vmesh, sizeof primary_vmesh,
                                 fit_stem, 0) != 0 ||
        strlen(primary_vmesh) + strlen(".vmesh") + 1 > sizeof primary_vmesh)
        goto done;
    strcat(primary_vmesh, ".vmesh");
    for (size_t layer = 0; layer < atlas->declared_layers; layer++) {
        QpProjectiveAtlasChart chart;
        char input_stem[QP_MAX_PATH];
        memset(&chart, 0, sizeof chart);
        if (qp_projective_layer_stem(input_stem, sizeof input_stem,
                                     fit_stem, layer) != 0 ||
            qp_projective_layer_stem(chart.output_stem,
                                     sizeof chart.output_stem,
                                     solid_stem, layer) != 0)
            goto done;
        if (snprintf(chart.input_vmesh, sizeof chart.input_vmesh,
                     "%s.vmesh", input_stem) < 0 ||
            strlen(chart.input_vmesh) + 1 >= sizeof chart.input_vmesh ||
            snprintf(chart.output_vmesh, sizeof chart.output_vmesh,
                     "%s.vmesh", chart.output_stem) < 0 ||
            strlen(chart.output_vmesh) + 1 >= sizeof chart.output_vmesh)
            goto done;
        if (!qp_exists(chart.input_vmesh)) {
            if (layer == 0) {
                qp_logf("[stage 4 atlas] missing primary chart %s\n",
                        chart.input_vmesh);
                goto done;
            }
            continue;
        }
        if (!MeshBin_looks_complete(chart.input_vmesh)) {
            qp_logf("[stage 4 atlas] malformed chart %s\n",
                    chart.input_vmesh);
            goto done;
        }
        chart.layer = layer;
        atlas->charts[atlas->n++] = chart;
        atlas->has_surface[layer] = 1;
    }
    if (atlas->n == 0 || strcmp(atlas->charts[0].input_vmesh,
                                primary_vmesh) != 0)
        goto done;
    rc = 0;
done:
    Arena_dispose(&arena);
    if (rc != 0) qp_projective_atlas_dispose(atlas);
    return rc;
}

static int qp_projective_atlas_write_manifest(
    const QpProjectiveAtlas *atlas, const char *solid_stem)
{
    char path[QP_MAX_PATH], stats_escaped[QP_MAX_PATH * 2];
    FILE *f = NULL;
    int close_rc = -1;
    if (snprintf(path, sizeof path, "%s_atlas.json", solid_stem) < 0 ||
        strlen(path) + 1 >= sizeof path)
        return -1;
    qp_json_escape(atlas->source_stats, stats_escaped, sizeof stats_escaped);
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f,
            "{\n  \"schema\": "
            "\"vesuvius-quadribbon-solid-projective-atlas-v1\",\n"
            "  \"source_stats\": \"%s\",\n"
            "  \"declared_grid_layers\": %zu,\n"
            "  \"surface_chart_count\": %zu,\n"
            "  \"authority\": \"charts are independent members at unchanged "
            "absolute U/V; storage layer is non-semantic\",\n"
            "  \"cross_chart_merge\": false,\n"
            "  \"charts\": [\n",
            stats_escaped, atlas->declared_layers, atlas->n);
    for (size_t i = 0; i < atlas->n; i++) {
        char input_escaped[QP_MAX_PATH * 2], output_escaped[QP_MAX_PATH * 2];
        qp_json_escape(atlas->charts[i].input_vmesh, input_escaped,
                       sizeof input_escaped);
        qp_json_escape(atlas->charts[i].output_vmesh, output_escaped,
                       sizeof output_escaped);
        fprintf(f,
                "    { \"layer\": %zu, \"input\": \"%s\", "
                "\"output\": \"%s\" }%s\n",
                atlas->charts[i].layer, input_escaped, output_escaped,
                i + 1 < atlas->n ? "," : "");
    }
    fprintf(f, "  ],\n  \"layers_without_surface_vmesh\": [");
    {
        int comma = 0;
        for (size_t layer = 0; layer < atlas->declared_layers; layer++) {
            if (atlas->has_surface[layer]) continue;
            fprintf(f, "%s%zu", comma ? ", " : "", layer);
            comma = 1;
        }
    }
    fprintf(f, "]\n}\n");
    close_rc = fclose(f);
    if (close_rc != 0) return -1;
    qp_logf("[stage 4 atlas] authority: %s (%zu surface chart(s), %zu "
            "declared layer(s))\n", path, atlas->n,
            atlas->declared_layers);
    return 0;
}

static int qp_projective_atlas_complete(const QpProjectiveAtlas *atlas,
                                        const char *solid_stem)
{
    Arena_T arena = Arena_new();
    const JsonValue *root = NULL, *charts = NULL;
    const char *err = NULL, *schema = NULL;
    char path[QP_MAX_PATH];
    int ok = 0;
    if (arena == NULL ||
        snprintf(path, sizeof path, "%s_atlas.json", solid_stem) < 0 ||
        strlen(path) + 1 >= sizeof path)
        goto done;
    root = Json_parse_file(arena, path, &err);
    schema = root != NULL
                 ? Json_as_string(Json_object_get(root, "schema"))
                 : NULL;
    charts = root != NULL ? Json_object_get(root, "charts") : NULL;
    if (schema == NULL ||
        strcmp(schema, "vesuvius-quadribbon-solid-projective-atlas-v1") != 0 ||
        Json_member_long(root, "declared_grid_layers", -1) !=
            (long)atlas->declared_layers ||
        Json_array_len(charts) != atlas->n)
        goto done;
    for (size_t i = 0; i < atlas->n; i++) {
        const JsonValue *entry = Json_array_get(charts, i);
        const char *input = Json_as_string(Json_object_get(entry, "input"));
        const char *output = Json_as_string(Json_object_get(entry, "output"));
        if (Json_member_long(entry, "layer", -1) !=
                (long)atlas->charts[i].layer ||
            input == NULL || strcmp(input, atlas->charts[i].input_vmesh) != 0 ||
            output == NULL || strcmp(output, atlas->charts[i].output_vmesh) != 0 ||
            !QuadribbonSolidify_complete(atlas->charts[i].output_stem))
            goto done;
    }
    ok = 1;
done:
    if (arena != NULL) Arena_dispose(&arena);
    return ok;
}

/* Regenerate every member unconditionally.  Copy-only projective solidify is
 * linear, while timestamp-based reuse cannot prove that sidecar evidence at a
 * stable path is unchanged.  Invalidate the old manifest first and publish a
 * new one only after every member is complete. */
static int qp_stage_solid_projective_atlas(const QpConfig *cfg,
                                           const char *fit_stem,
                                           const char *solid_stem,
                                           const char *tag)
{
    QpProjectiveAtlas atlas;
    char manifest[QP_MAX_PATH];
    int rc = -1;
    memset(&atlas, 0, sizeof atlas);
    if (qp_projective_atlas_discover(fit_stem, solid_stem, &atlas) != 0)
        return -1;
    if (snprintf(manifest, sizeof manifest, "%s_atlas.json", solid_stem) < 0 ||
        strlen(manifest) + 1 >= sizeof manifest)
        goto done;
    if (qp_exists(manifest) && remove(manifest) != 0) {
        qp_logf("[%s] cannot invalidate stale atlas manifest %s\n", tag,
                manifest);
        goto done;
    }
    for (size_t i = 0; i < atlas.n; i++) {
        char chart_tag[128];
        snprintf(chart_tag, sizeof chart_tag, "%s chart %zu", tag,
                 atlas.charts[i].layer);
        if (qp_stage_solid(cfg, atlas.charts[i].input_vmesh,
                           atlas.charts[i].output_stem, chart_tag) != 0 ||
            !QuadribbonSolidify_complete(atlas.charts[i].output_stem))
            goto done;
    }
    /* Each member writer publishes legacy directory-wide aliases.  The last
     * extras chart must not remain there: restore layer 0 after the full loop.
     * Atlas-aware consumers use the manifest/member sidecars directly. */
    if (qp_publish_stem_plain_sidecars(solid_stem) != 0) goto done;
    if (qp_projective_atlas_write_manifest(&atlas, solid_stem) != 0 ||
        !qp_projective_atlas_complete(&atlas, solid_stem)) {
        qp_logf("[%s] atlas postcondition FAILED\n", tag);
        goto done;
    }
    rc = 0;
done:
    qp_projective_atlas_dispose(&atlas);
    return rc;
}

/* ---- stage 5: optimize the solid ribbon, refitting as needed --------------
 * Each round untangles (radial collision shell on every vertex) and then
 * REFITS: the same solidify call on the moved ribbon, whose provenance
 * sidecar says which cells are observations (kept, now moved) and which were
 * fills (re-solved).  Rounds stop at a conflict-free state, at the cap, or
 * when long conflicts stop dropping. */
static int qp_stage_optimize(const QpConfig *cfg, const char *solid_stem,
                             const char *opt_dir, const char *cert_vmesh,
                             char *out_stem, size_t out_cap)
{
    char prev_stem[QP_MAX_PATH], report_path[QP_MAX_PATH], best_stem[QP_MAX_PATH];
    struct { long long long_in, long_out; int fresh, accepted; long gates; double cov, seams; } R[QS_OPTIMIZE_ROUNDS];
    int rounds_run = 0;
    long long long_prev = -1;
    const char *stop_reason = "round cap";
    /* the deliverable is the BEST ribbon by the verdict, never the last round:
     * measured 2026-09-08 on the F0 4x5x5 (output/quadribbon_f0_4x5x5_20260908)
     * the untangle raised its own long conflicts 49,623 -> 293,652 and the two
     * rounds took COVERAGE 98.9 -> 70.7 -> 63.6% and depth-seams 0.438 ->
     * 0.646 -> 0.887%, yet the last refit was published.  A round is kept only
     * when its conflicts drop AND the verdict does not regress (failed gates,
     * source coverage within 0.5 points, seam share within 0.05 points). */
    long best_gates = LONG_MAX;
    double best_cov = 0.0, best_seams = 1.0;
    snprintf(prev_stem, sizeof prev_stem, "%s", solid_stem);
    snprintf(best_stem, sizeof best_stem, "%s", solid_stem);
    memset(R, 0, sizeof R);
    {
        char sv[QP_MAX_PATH], sdir[QP_MAX_PATH];
        size_t n = strlen(solid_stem);
        const char *slash = NULL;
        for (size_t i = n; i > 0; i--)
            if (solid_stem[i - 1] == '/' || solid_stem[i - 1] == '\\') { slash = solid_stem + i - 1; break; }
        snprintf(sdir, sizeof sdir, "%.*s", slash != NULL ? (int)(slash - solid_stem) : 1,
                 slash != NULL ? solid_stem : ".");
        qp_join(sv, sizeof sv, sdir, "verdict.json");
        if (qp_exists(sv)) {
            best_gates = qp_read_json_long(sv, "gates_failed", LONG_MAX);
            best_cov = qp_read_json_double2(sv, "coverage", "source_frac", 0.0);
            best_seams = qp_read_json_double2(sv, "depth_seams", "seam_share", 1.0);
            qp_logf("[stage 5 optimize] baseline (solid): gates failed %ld, coverage %.3f, seam share %.4f\n",
                    best_gates, best_cov, best_seams);
        }
    }
    for (int r = 1; r <= QS_OPTIMIZE_ROUNDS; r++) {
        char round_dir[QP_MAX_PATH], unt_stem[QP_MAX_PATH], unt_vmesh[QP_MAX_PATH];
        char refit_stem[QP_MAX_PATH], refit_vmesh[QP_MAX_PATH];
        char prev_vmesh[QP_MAX_PATH], unt_json[QP_MAX_PATH], verdict[QP_MAX_PATH];
        char tag[32];
        QuadribbonUntangleStats st;
        int fresh = 0;
        memset(&st, 0, sizeof st);
        snprintf(round_dir, sizeof round_dir, "%s/round_%d", opt_dir, r);
        qp_mkdir(round_dir);
        qp_join(unt_stem, sizeof unt_stem, round_dir, "untangled_ribbon");
        snprintf(unt_vmesh, sizeof unt_vmesh, "%s.vmesh", unt_stem);
        qp_join(refit_stem, sizeof refit_stem, round_dir, "refit_ribbon");
        snprintf(refit_vmesh, sizeof refit_vmesh, "%s.vmesh", refit_stem);
        snprintf(prev_vmesh, sizeof prev_vmesh, "%s.vmesh", prev_stem);
        qp_join(unt_json, sizeof unt_json, round_dir, "untangle.json");
        qp_join(verdict, sizeof verdict, round_dir, "verdict.json");
        snprintf(tag, sizeof tag, "optimize-round-%d", r);
        if (!QuadribbonSolidify_complete(refit_stem)) {
            if (!MeshBin_looks_complete(unt_vmesh)) {
                double ts = ves_clock_sec();
                if (qp_stage_untangle(cfg, prev_vmesh, unt_vmesh, &st) != 0)
                    return -1;
                qp_write_untangle_json(unt_json, &st);
                fresh = 1;
                qp_logf("[stage 5 optimize] round %d untangle done (%.1fs)\n",
                        r, ves_clock_sec() - ts);
            } else {
                qp_logf("[stage 5 optimize] round %d resume: untangled\n", r);
            }
            if (qp_copy_sidecars(prev_stem, unt_stem) != 0) return -1;
            if (qp_stage_solid(cfg, unt_vmesh, refit_stem, "stage 5 refit") != 0)
                return -1;
        } else {
            qp_logf("[stage 5 optimize] round %d resume: refit complete\n", r);
        }
        R[r - 1].fresh = fresh;
        R[r - 1].long_in = fresh ? (long long)st.input_long_conflicts
                                 : qp_read_json_long(unt_json, "input_long_conflicts", -1);
        R[r - 1].long_out = fresh ? (long long)st.output_long_conflicts
                                  : qp_read_json_long(unt_json, "output_long_conflicts", -1);
        qp_verdict(cfg, tag, refit_vmesh, refit_stem, cert_vmesh, verdict);
        rounds_run = r;
        R[r - 1].gates = qp_read_json_long(verdict, "gates_failed", LONG_MAX);
        R[r - 1].cov = qp_read_json_double2(verdict, "coverage", "source_frac", 0.0);
        R[r - 1].seams = qp_read_json_double2(verdict, "depth_seams", "seam_share", 1.0);
        {
            int worse_conflicts = R[r - 1].long_in >= 0 && R[r - 1].long_out > R[r - 1].long_in;
            int regressed = R[r - 1].gates > best_gates ||
                            R[r - 1].cov < best_cov - 0.005 ||
                            R[r - 1].seams > best_seams + 0.0005;
            if (worse_conflicts || regressed) {
                R[r - 1].accepted = 0;
                stop_reason = worse_conflicts ? "round rejected: untangle raised its own long conflicts"
                                              : "round rejected: verdict regressed";
                qp_logf("[stage 5 optimize] round %d REJECTED (%s): conflicts %lld -> %lld, gates failed %ld (best %ld), "
                        "coverage %.3f (best %.3f), seam share %.4f (best %.4f); keeping %s\n",
                        r, stop_reason, R[r - 1].long_in, R[r - 1].long_out, R[r - 1].gates, best_gates,
                        R[r - 1].cov, best_cov, R[r - 1].seams, best_seams, best_stem);
                break;
            }
        }
        R[r - 1].accepted = 1;
        best_gates = R[r - 1].gates; best_cov = R[r - 1].cov; best_seams = R[r - 1].seams;
        snprintf(best_stem, sizeof best_stem, "%s", refit_stem);
        snprintf(prev_stem, sizeof prev_stem, "%s", refit_stem);
        qp_logf("[stage 5 optimize] round %d accepted: conflicts %lld -> %lld, gates failed %ld, coverage %.3f, seam share %.4f\n",
                r, R[r - 1].long_in, R[r - 1].long_out, R[r - 1].gates, R[r - 1].cov, R[r - 1].seams);
        if (R[r - 1].long_out == 0) { stop_reason = "conflict-free"; break; }
        if (r >= 2 && long_prev > 0 &&
            (double)R[r - 1].long_out >
                (1.0 - QS_OPTIMIZE_MIN_GAIN) * (double)long_prev) {
            stop_reason = "no further gain";
            break;
        }
        long_prev = R[r - 1].long_out;
    }
    /* publish the best, not the last */
    snprintf(prev_stem, sizeof prev_stem, "%s", best_stem);
#if QS_REPARAM_FINAL
    {
        /* metric re-solve of U on the optimized geometry (retained-gauge
         * lane keeps this OFF; here as the measured lever) */
        char rp_dir[QP_MAX_PATH], rp_stem[QP_MAX_PATH], rp_vmesh[QP_MAX_PATH];
        char prev_vmesh[QP_MAX_PATH], tool[QP_MAX_PATH];
        char axis_y_s[32], axis_x_s[32], pitch_s[32];
        qp_join(rp_dir, sizeof rp_dir, opt_dir, "reparam");
        qp_join(rp_stem, sizeof rp_stem, rp_dir, "reopt_ribbon");
        snprintf(rp_vmesh, sizeof rp_vmesh, "%s.vmesh", rp_stem);
        snprintf(prev_vmesh, sizeof prev_vmesh, "%s.vmesh", prev_stem);
        if (!MeshBin_looks_complete(rp_vmesh)) {
            const char *argv2[14];
            size_t a = 0;
            qp_mkdir(rp_dir);
            qp_tool_path(tool, sizeof tool, "scroll_ribbon");
            snprintf(axis_y_s, sizeof axis_y_s, "%.3f", cfg->axis_y);
            snprintf(axis_x_s, sizeof axis_x_s, "%.3f", cfg->axis_x);
            snprintf(pitch_s, sizeof pitch_s, "%.3f", cfg->wrap_spacing);
            argv2[a++] = tool;
            argv2[a++] = prev_vmesh;
            argv2[a++] = rp_dir;
            argv2[a++] = "--id"; argv2[a++] = "reopt";
            argv2[a++] = "--umb-y"; argv2[a++] = axis_y_s;
            argv2[a++] = "--umb-x"; argv2[a++] = axis_x_s;
            argv2[a++] = "--wrap-spacing"; argv2[a++] = pitch_s;
            argv2[a++] = "--preserve-input-topology";
            argv2[a] = NULL;
            if (qp_spawn("scroll_ribbon(reparam)", argv2) != 0) return -1;
        }
        if (qp_copy_sidecars(prev_stem, rp_stem) != 0) return -1;
        snprintf(prev_stem, sizeof prev_stem, "%s", rp_stem);
    }
#endif
    /* publish: sidecars first, the VMESH last, so completeness implies both */
    snprintf(out_stem, out_cap, "%s/optimized_ribbon", opt_dir);
    if (qp_copy_sidecars(prev_stem, out_stem) != 0) return -1;
    {
        char src[QP_MAX_PATH], dst[QP_MAX_PATH], src_dir[QP_MAX_PATH];
        size_t n = strlen(prev_stem);
        const char *slash = NULL;
        for (size_t i = n; i > 0; i--)
            if (prev_stem[i - 1] == '/' || prev_stem[i - 1] == '\\') { slash = prev_stem + i - 1; break; }
        snprintf(src_dir, sizeof src_dir, "%.*s",
                 slash != NULL ? (int)(slash - prev_stem) : 1,
                 slash != NULL ? prev_stem : ".");
        qp_join(src, sizeof src, src_dir, "verdict.json");
        snprintf(dst, sizeof dst, "%s_verdict.json", out_stem);
        if (qp_exists(src) && qp_copy_file(src, dst) != 0) return -1;
        snprintf(src, sizeof src, "%s.vmesh", prev_stem);
        snprintf(dst, sizeof dst, "%s.vmesh", out_stem);
        if (qp_copy_file(src, dst) != 0) return -1;
    }
    qp_join(report_path, sizeof report_path, opt_dir, "optimize_report.json");
    {
        FILE *f = fopen(report_path, "wb");
        if (f != NULL) {
            fprintf(f, "{ \"schema\": \"vesuvius-quadribbon-optimize-v1\",\n"
                       "  \"rounds_run\": %d, \"round_cap\": %d, \"min_gain\": %.3f, "
                       "\"stop\": \"%s\", \"reparam_final\": %d,\n  \"rounds\": [\n",
                    rounds_run, QS_OPTIMIZE_ROUNDS, QS_OPTIMIZE_MIN_GAIN,
                    stop_reason, QS_REPARAM_FINAL);
            for (int r = 0; r < rounds_run; r++)
                fprintf(f, "    { \"round\": %d, \"long_conflicts_in\": %lld, "
                           "\"long_conflicts_out\": %lld, \"fresh\": %d, \"accepted\": %d, "
                           "\"gates_failed\": %ld, \"coverage_source_frac\": %.6f, \"seam_share\": %.6f }%s\n",
                        r + 1, R[r].long_in, R[r].long_out, R[r].fresh, R[r].accepted,
                        R[r].gates, R[r].cov, R[r].seams,
                        r + 1 < rounds_run ? "," : "");
            fprintf(f, "  ],\n  \"published\": \"%s.vmesh\",\n  \"output\": \"%s.vmesh\" }\n", best_stem, out_stem);
            fclose(f);
        }
    }
    qp_logf("[stage 5 optimize] %d round(s), stop: %s -> %s.vmesh\n",
            rounds_run, stop_reason, out_stem);
    return 0;
}

/* ---- main ---------------------------------------------------------------- */

static void qp_stop(int stage, const char *name, double t0)
{
    qp_logf("quadribbon: STOP after stage %d %s as requested (%.1fs total)\n",
            stage, name, ves_clock_sec() - t0);
}

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    QpConfig cfg;
    char logs[QP_MAX_PATH];
    char mesh_dir[QP_MAX_PATH], mesh_vmesh[QP_MAX_PATH], mesh_table[QP_MAX_PATH];
    char unwrap_dir[QP_MAX_PATH], cert_vmesh[QP_MAX_PATH];
    char fit_dir[QP_MAX_PATH], fit_vmesh[QP_MAX_PATH], fit_stem[QP_MAX_PATH];
    char fit_verdict[QP_MAX_PATH];
    char solid_dir[QP_MAX_PATH], solid_stem[QP_MAX_PATH], solid_vmesh[QP_MAX_PATH];
    char solid_verdict[QP_MAX_PATH], solid_labels[QP_MAX_PATH];
    char opt_dir[QP_MAX_PATH], opt_stem[QP_MAX_PATH], opt_vmesh[QP_MAX_PATH];
    char opt_labels[QP_MAX_PATH];
    const char *input = NULL, *out_dir = NULL, *config_path = NULL, *scroll_model_dir = NULL;
    const char *quad_model_dir=NULL;
    const char *quad_field_path=NULL;
    int stop_stage = 0, is_vmesh = 0, trust_gauge = 0, build_scroll_model = 0;
    int observe_material=0;
    double t0 = 0.0;

    if (argc == 2 && strcmp(argv[1], "--selftest-material") == 0) {
        int failures=MaterialEvidence_selftest()!=0;
        failures+=MaterialFront_selftest()!=0;
        failures+=WindingMRF_selftest()!=0;
        return failures ? 1 : 0;
    }
    if (argc == 2 && strcmp(argv[1], "--selftest-coordinates") == 0) {
        int failures = ScrollCoordinate_selftest() != 0;
        failures += ScrollSource_selftest() != 0;
        failures += WindingRegister_selftest() != 0;
        return failures ? 1 : 0;
    }
    if (argc == 2 && strcmp(argv[1], "--selftest-lod") == 0)
        return RibbonLod_selftest() == 0 ? 0 : 1;
    if (argc == 2 && strcmp(argv[1], "--selftest-kibble") == 0)
        return MeshKibble_selftest() == 0 ? 0 : 1;
    if (argc == 2 && strcmp(argv[1], "--selftest-quad-field") == 0)
        return QuadField_selftest() == 0 ? 0 : 1;
    if (argc == 2 && strcmp(argv[1], "--selftest-global-ribbon-fit") == 0)
        return QuadRibbonFit_selftest() == 0 ? 0 : 1;
    if (argc == 4 && strcmp(argv[1], "--fit-global-ribbon") == 0)
        return QuadRibbonFit_write(argv[2],argv[3]) == 0 ? 0 : 1;
    if (argc == 5 && strcmp(argv[1], "--fit-global-atlas") == 0) {
        int rc=-1;
        arena=Arena_new();
        if (qp_load_config(arena,argv[4],&cfg)==0)
            rc=QuadRibbonFit_write_atlas(argv[2],argv[3],cfg.axis_warp_armed ? &cfg.axis_warp : NULL);
        Arena_dispose(&arena);
        return rc==0 ? 0 : 1;
    }
    if (argc == 8 && strcmp(argv[1], "--ribbon-lod") == 0)
        return RibbonLod_write(argv[2], argv[3], strtod(argv[4], NULL),
                               strtod(argv[5], NULL), atoi(argv[6]),
                               strtod(argv[7], NULL)) == 0 ? 0 : 1;
    if (argc == 2 && strcmp(argv[1], "--selftest-solidify") == 0) {
        int fails = QuadribbonSolidify_selftest();
        fprintf(stderr, "quadribbon --selftest-solidify: %s\n",
                fails == 0 ? "PASS" : "FAIL");
        return fails == 0 ? 0 : 1;
    }
    if (argc == 5 && (!strcmp(argv[1], "--axis-world-copy") || !strcmp(argv[1],"--axis-metric-copy"))) {
        int rc = -1;
        arena = Arena_new();
        if (arena == NULL) return 1;
        if (qp_load_config(arena, argv[4], &cfg) == 0) {
            if (!cfg.axis_warp_armed) {
                fprintf(stderr, "quadribbon: %s has no geometry.axis_table\n",
                        argv[4]);
            } else {
                rc = qp_write_warped_copy(&cfg, argv[2], argv[3], !strcmp(argv[1],"--axis-metric-copy") ? 1 : -1);
            }
        }
        Arena_dispose(&arena);
        return rc == 0 ? 0 : 1;
    }
    if (argc==5 && !strcmp(argv[1],"--probe-axis-ct")) {
        int rc=-1; MeshBinData mesh={0}; QsConfig qc={0};
        arena=Arena_new();
        if (!arena) return 1;
        if (qp_load_config(arena,argv[4],&cfg)==0 && MeshBin_read_arena(arena,argv[2],&mesh)==0) {
            double *points=(double*)ARENA_ALLOC(arena,(mesh.nv ? mesh.nv : 1)*3*sizeof(double));
            double *world=(double*)ARENA_ALLOC(arena,(mesh.nv ? mesh.nv : 1)*3*sizeof(double));
            uint8_t *values=(uint8_t*)ARENA_ALLOC(arena,mesh.nv ? mesh.nv : 1);
            for (size_t i=0;i<mesh.nv*3;i++) points[i]=mesh.verts[i];
            qc.coordinate_warp=cfg.axis_warp_armed ? &cfg.axis_warp : NULL; qc.raw_zarr=cfg.raw_path;
            if (QuadribbonSolidify_probe_ct(&qc,points,mesh.nv,world,values)==0) {
                char temporary[QP_MAX_PATH];
                if (snprintf(temporary,sizeof temporary,"%s.tmp",argv[3])<(int)sizeof temporary) {
                    FILE *f=fopen(temporary,"wb");
                    if (f) {
                        fprintf(f,"vertex\tworld_z\tworld_y\tworld_x\tintensity\n");
                        for (size_t i=0;i<mesh.nv;i++) fprintf(f,"%zu\t%.17g\t%.17g\t%.17g\t%u\n",i,world[3*i],world[3*i+1],world[3*i+2],(unsigned)values[i]);
                        int failed=ferror(f);
                        if (fclose(f)==0 && !failed) rc=qp_replace_file(temporary,argv[3]);
                    }
                }
            } else fprintf(stderr,"[axis CT probe] invalid inverse coordinate or unavailable RAW; no output published\n");
        }
        AxisWarp_dispose(&cfg.axis_warp); Arena_dispose(&arena); return rc==0 ? 0 : 1;
    }
    if (argc == 5 && strcmp(argv[1], "--solidify") == 0) {
        /* the production Stage-4 solid strip fitted to an EXTERNAL lattice ribbon
         * (rows = z-slices, sidecars _support.u8 / _phase.f32 [+ _lane / _material /
         * _provenance / _stats]): the assembly-to-solid-strip route for the chart
         * assembler's sheet (sheet_assemble.exe writes assembled_ribbon.*) */
        int rc = -1;
        qp_capture_exe_dir(argv[0]);
        arena = Arena_new();
        if (arena == NULL) return 1;
        if (qp_load_config(arena, argv[4], &cfg) == 0)
            rc = qp_stage_solid(&cfg, argv[2], argv[3], "solidify");
        Arena_dispose(&arena);
        return rc == 0 ? 0 : 1;
    }
    if (argc == 5 && strcmp(argv[1], "--projective-solidify") == 0) {
        int rc = -1;
        qp_capture_exe_dir(argv[0]);
        arena = Arena_new();
        if (arena == NULL) return 1;
        if (qp_load_config(arena, argv[4], &cfg) == 0) {
            /* This is the exact production Stage-4 authority path, exposed so
             * nested-domain tests can exercise it without also paying for CT
             * bakes and sheet presentation. */
            g_qp_canonical_blocks = 1;
            rc = qp_stage_solid(&cfg, argv[2], argv[3],
                                "projective solidify");
        }
        Arena_dispose(&arena);
        return rc == 0 ? 0 : 1;
    }
    if (argc == 5 && strcmp(argv[1], "--projective-solidify-atlas") == 0) {
        int rc = -1;
        qp_capture_exe_dir(argv[0]);
        arena = Arena_new();
        if (arena == NULL) return 1;
        if (qp_load_config(arena, argv[4], &cfg) == 0) {
            /* Stem-to-stem authority path.  The stage-3 stats declare the
             * collision-layer address space; every bakeable member is
             * solidified independently and recorded in one manifest. */
            g_qp_canonical_blocks = 1;
            rc = qp_stage_solid_projective_atlas(
                &cfg, argv[2], argv[3], "projective solidify atlas");
        }
        Arena_dispose(&arena);
        return rc == 0 ? 0 : 1;
    }
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        int fails = 0;
        if (SheetComposite_selftest() != 0) fails++;
        if (QuadribbonUntangle_selftest() != 0) fails++;
        if (Sparse_selftest() != 0) fails++;
        if (MeshBin_selftest() != 0) fails++;
        if (MeshKibble_selftest() != 0) fails++;
        if (ScrollCoordinate_selftest() != 0) fails++;
        if (ScrollSource_selftest() != 0) fails++;
        if (MaterialEvidence_selftest() != 0) fails++;
        if (MaterialFront_selftest() != 0) fails++;
        if (WindingMRF_selftest() != 0) fails++;
        if (RibbonLod_selftest() != 0) fails++;
        if (QuadField_selftest() != 0) fails++;
        if (QuadRibbonFit_selftest() != 0) fails++;
        if (WindingRegister_selftest() != 0) fails++;
        if (QuadStrip_selftest() != 0) fails++;
        if (QuadribbonSolidify_selftest() != 0) fails++;
        if (ZarrU8_selftest() != 0) fails++;
        if (RidgeTrack_selftest() != 0) fails++;
        if (qp_core_curl_selftest() != 0) fails++;
        if (qp_tile_selftest() != 0) fails++;
        if (qp_canonical_block_selftest() != 0) fails++;
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
            } else if (strcmp(argv[i], "--observe-material") == 0) {
                observe_material=1;
            } else if (strcmp(argv[i], "--build-scroll-model") == 0) {
                build_scroll_model = 1;
            } else if (strcmp(argv[i], "--scroll-model") == 0 && i+1<argc) {
                scroll_model_dir=argv[++i];
            } else if (strcmp(argv[i], "--build-scroll-quads") == 0 && i+1<argc) {
                quad_model_dir=argv[++i];
            } else if (strcmp(argv[i], "--scroll-quad-field") == 0 && i+1<argc) {
                quad_field_path=argv[++i];
            } else if (strcmp(argv[i], "--trust-gauge") == 0) {
                trust_gauge = 1;
            } else if (strcmp(argv[i], "--canonical-blocks") == 0) {
                g_qp_canonical_blocks = 1;
            } else if (strcmp(argv[i], "--subgrid") == 0 && i + 6 < argc) {
                for (int k = 0; k < 6; k++)
                    g_qp_sg[k] = strtol(argv[i + 1 + k], NULL, 10);
                i += 6;
                g_qp_subgrid_armed = 1;
            } else if (strcmp(argv[i], "--tile") == 0 && i + 6 < argc) {
                for (int k = 0; k < 6; k++)
                    g_qp_tile[k] = strtod(argv[i + 1 + k], NULL);
                i += 6;
                g_qp_tile_armed = 1;
            } else if (input == NULL) {
                input = argv[i];
            } else if (out_dir == NULL) {
                out_dir = argv[i];
            } else if (config_path == NULL) {
                config_path = argv[i];
            } else {
                bad = 1;
            }
        }
        stop_stage = qp_stage_index(stop_name);
        if (bad || input == NULL || out_dir == NULL || stop_stage < 0) {
            fprintf(stderr,
                    "usage: quadribbon <mesh_dir|welded.vmesh> <out_dir> "
                    "[config.json] [--trust-gauge] [--stop-after <stage>]\n"
                    "       --trust-gauge  retain the stage-2 carried U in "
                    "the fast streaming fit\n"
                    "       --canonical-blocks  immutable overlapping 4x5x5 "
                    "winding certificates (pile inputs)\n"
                    "       --subgrid z0 z1 y0 y1 x0 x1  requested cube origins\n"
                    "       --build-scroll-model  build a fixed-source model in out_dir (no area filter)\n"
                    "       --observe-material  independent crop source/soft-evidence diagnostics; no UV qualification\n"
                    "       --scroll-model DIR  query fixed coefficients; currently requires --stop-after unwrap\n"
                    "       --build-scroll-quads DIR  experimental full-source adaptive quad fit from fixed model DIR\n"
                    "       --scroll-quad-field FILE  query a saved global quad level with --scroll-model and --stop-after fit\n"
                    "       --tile r0 r1 theta0 theta1 z0 z1  cylindrical cube tile (block lane)\n"
                    "       quadribbon --selftest\n"
                    "       quadribbon --selftest-solidify\n"
                    "       quadribbon --fit-global-ribbon <observations_world.vmesh> <out_dir>\n"
                    "                   supported global geometry fit; preserves established observation UVs\n"
                    "       quadribbon --fit-global-atlas <fit_stem> <out_dir> <config.json>\n"
                    "                   all declared layers, with original inter-layer clearance envelopes\n"
                    "       quadribbon --projective-solidify <fit.vmesh> "
                    "<out_stem> <config.json>\n"
                    "       quadribbon --projective-solidify-atlas <fit_stem> "
                    "<out_stem> <config.json>\n"
                    "       quadribbon --axis-world-copy <straight.vmesh> "
                    "<world.vmesh> <config.json>\n"
                    "  <mesh_dir>:      pre-weld per-cube pile (recursive "
                    "*_final_all.vmesh)\n"
                    "                   and/or plain world-frame *.vmesh files\n"
                    "  <welded.vmesh>:  a welded world-frame mesh\n"
                    "  <stage>: mesh|unwrap|fit|solid|optimize|sheet (or 1..6) "
                    "-- finish that\n"
                    "           stage (verdict + debug bake included), then exit; "
                    "rerun resumes from it\n"
                    "  policy: %s\n", QP_DEFAULT_CONFIG);
            return 2;
        }
        if (config_path == NULL) config_path = QP_DEFAULT_CONFIG;
    }
    qp_capture_exe_dir(argv[0]);
    t0 = ves_clock_sec();
    qp_log_t0 = t0;
    is_vmesh = qp_input_is_vmesh(input);
    if (g_qp_canonical_blocks && is_vmesh) {
        fprintf(stderr, "quadribbon: --canonical-blocks requires a per-cube "
                        "mesh pile, not a welded VMESH\n");
        return 2;
    }
    if (build_scroll_model && (is_vmesh || g_qp_subgrid_armed ||
                              g_qp_tile_armed || g_qp_canonical_blocks)) {
        fprintf(stderr, "quadribbon: model construction needs the complete fixed source pile, without query filters or parent certificates\n");
        return 2;
    }
    if (g_qp_canonical_blocks && g_qp_tile_armed) {
        fprintf(stderr, "quadribbon: --canonical-blocks and --tile are "
                        "different partition contracts and cannot be combined\n");
        return 2;
    }

    arena = Arena_new();
    if (arena == NULL) return 1;
    if (qp_load_config(arena, config_path, &cfg) != 0) return 1;
    g_qp_tile_axis[0] = cfg.axis_y;
    g_qp_tile_axis[1] = cfg.axis_x;
    if (cfg.axis_warp_armed) {
        g_qp_axis_warp = &cfg.axis_warp;
        if (cfg.axis_warp.physical)
            fprintf(stderr,"quadribbon: physical shaft (%zu ordered samples, %.3f um arc)\n",
                    cfg.axis_warp.n,ShaftWarp_length_um(cfg.axis_warp.physical));
        else fprintf(stderr, "quadribbon: axis table %s (%zu samples) -> straightened lane\n",
                     cfg.axis_table_path, cfg.axis_warp.n);
    }
    if (g_qp_tile_armed && cfg.axis_warp.physical)
        fprintf(stderr,"quadribbon: conservative physical tile r<=%.3f s=[%.3f,%.3f] metric voxels; "
                "whole cubes may extend outside angular and inner-radius bounds\n",
                g_qp_tile[1],g_qp_tile[4],g_qp_tile[5]);
    else if (g_qp_tile_armed)
        fprintf(stderr, "quadribbon: tile r=[%.0f,%.0f] theta=[%.0f,%.0f] deg "
                "z=[%.0f,%.0f] about (%.1f,%.1f)\n",
                g_qp_tile[0], g_qp_tile[1], g_qp_tile[2], g_qp_tile[3],
                g_qp_tile[4], g_qp_tile[5], cfg.axis_y, cfg.axis_x);

    qp_mkdir(out_dir);
    {
        char frame[QP_MAX_PATH],old_mesh[QP_MAX_PATH];
        qp_join(frame,sizeof frame,out_dir,"coordinate_frame.json");
        qp_join(old_mesh,sizeof old_mesh,out_dir,"stage1_mesh/mesh.vmesh");
        if (qp_coordinate_identity(&cfg,frame,old_mesh)!=0) { Arena_dispose(&arena); return 1; }
    }
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
        qp_logf("quadribbon: input=%s (%s) out_dir=%s config=%s\n"
                "  axis=(%.1f,%.1f) pitch=%.3f threads=%ld raw=%s "
                "window=[%ld,%ld] dark=%ld core_wall=%.1f fit=%s%s\n",
                input, is_vmesh ? "welded vmesh" : "mesh pile", out_dir,
                config_path, cfg.axis_y, cfg.axis_x, cfg.wrap_spacing,
                cfg.threads, cfg.raw_path, cfg.window_lo, cfg.window_hi,
                 cfg.dark, cfg.core_wall_radius,
                 trust_gauge ? "trust-gauge" : "full-solve", stop_note);
        if (g_qp_canonical_blocks)
            qp_logf("  unwrap-domain=canonical-projective-blocks\n");
    }

    if (observe_material) {
        ScrollModelOptions options={0};
        int result=-1;
        if (build_scroll_model || scroll_model_dir || quad_model_dir || quad_field_path ||
            g_qp_tile_armed || g_qp_canonical_blocks || is_vmesh || (stop_stage && stop_stage!=1)) {
            qp_logf("[material observation] source/evidence stage only; incompatible modes refused\n");
            Arena_dispose(&arena); return 1;
        }
        options.axis=cfg.axis_warp_armed ? &cfg.axis_warp : NULL;
        options.axis_y=cfg.axis_y; options.axis_x=cfg.axis_x;
        options.pitch=cfg.wrap_spacing; options.core_radius=cfg.core_wall_radius;
        options.sense=cfg.winding_sense; options.chunk=cfg.raw_chunk;
        result=MaterialObservation_build(arena,input,out_dir,g_qp_subgrid_armed ? g_qp_sg : NULL,&options);
        Arena_dispose(&arena); return result==0 ? 0 : 1;
    }

    /* Fixed-source query route: no crop-local solve, certificate or packing.
     * Later sparse-fit stages will consume the same immutable model. */
    if (quad_model_dir) {
        ScrollModel_T model=NULL;
        int rc=-1;
        if (build_scroll_model || scroll_model_dir || quad_field_path || g_qp_subgrid_armed ||
            g_qp_tile_armed || g_qp_canonical_blocks || is_vmesh) {
            qp_logf("[quad field] requires the full fixed source catalog; query bounds and legacy modes refused\n");
            Arena_dispose(&arena); return 1;
        }
        rc=ScrollModel_load(arena,quad_model_dir,&model);
        if (rc==0 && ScrollModel_coordinate_fingerprint(model)!=AxisWarp_fingerprint(&cfg.axis_warp)) {
            qp_logf("[quad field] loaded model coordinate frame differs from config\n"); rc=-1;
        }
        if (rc==0) rc=ScrollModel_validate_catalog(arena,model,input);
        if (rc==0) rc=ScrollModel_fit_quads(arena,model,out_dir);
        qp_logf("[quad field] %s; physical/topology qualification remains required\n",rc==0 ? "sample-error convergence" : "INCOMPLETE; inspect retained levels");
        ScrollModel_release(model);
        Arena_dispose(&arena); return rc==0 ? 0 : 1;
    }
    if (scroll_model_dir) {
        ScrollModel_T model=NULL;
        int rc=-1;
        if (build_scroll_model || g_qp_tile_armed || g_qp_canonical_blocks || is_vmesh ||
            (quad_field_path ? stop_stage!=3 : stop_stage!=2)) {
            qp_logf("[scroll model] requires --stop-after unwrap, or a saved quad field and --stop-after fit; incompatible modes refused\n");
            Arena_dispose(&arena); return 1;
        }
        rc=ScrollModel_load(arena,scroll_model_dir,&model);
        if (rc==0 && ScrollModel_coordinate_fingerprint(model)!=AxisWarp_fingerprint(&cfg.axis_warp)) {
            qp_logf("[scroll model] loaded model coordinate frame differs from config\n"); rc=-1;
        }
        if (rc==0) rc=ScrollModel_validate_catalog(arena,model,input);
        if (rc==0) rc=ScrollModel_write_area(arena,model,g_qp_subgrid_armed ? g_qp_sg : NULL,out_dir);
        if (rc==0 && quad_field_path) rc=ScrollModel_write_quad_area(arena,model,quad_field_path,g_qp_subgrid_armed ? g_qp_sg : NULL,out_dir);
        qp_logf("[scroll model query] %s (physical qualification pending)\n",rc==0 ? "written" : "FAILED");
        ScrollModel_release(model);
        Arena_dispose(&arena); return rc==0 ? 0 : 1;
    }
    if (quad_field_path) {
        qp_logf("[quad field] --scroll-quad-field requires --scroll-model\n");
        Arena_dispose(&arena); return 1;
    }

    /* ASSEMBLY ROUTE: an external lattice ribbon (the chart assembler's
     * assembled_ribbon.vmesh with its sidecars) IS the fitted ribbon; stages
     * 1-3 (pile concat, certificate, claims fit) are not run.  The verdict's
     * coverage counts assembly.source. */
    if (cfg.assembly_ribbon[0] != '\0') {
        size_t n = strlen(cfg.assembly_ribbon);
        snprintf(fit_vmesh, sizeof fit_vmesh, "%s", cfg.assembly_ribbon);
        if (n > 6 && strcmp(cfg.assembly_ribbon + n - 6, ".vmesh") == 0)
            snprintf(fit_stem, sizeof fit_stem, "%.*s", (int)(n - 6), cfg.assembly_ribbon);
        else
            snprintf(fit_stem, sizeof fit_stem, "%s", cfg.assembly_ribbon);
        snprintf(cert_vmesh, sizeof cert_vmesh, "%s", cfg.assembly_source);
        qp_logf("[assembly] external ribbon %s (source %s): stages 1-3 skipped\n",
                fit_vmesh, cert_vmesh[0] != '\0' ? cert_vmesh : "(none)");
        goto stage4;
    }
    /* stage 1: mesh */
    if (build_scroll_model) {
        ScrollModelOptions model_options = {0};
        ScrollModel_T model = NULL;
        int rc = 0;
        model_options.axis = cfg.axis_warp_armed ? &cfg.axis_warp : NULL;
        model_options.axis_y = cfg.axis_y;
        model_options.axis_x = cfg.axis_x;
        model_options.pitch = cfg.wrap_spacing;
        model_options.core_radius = cfg.core_wall_radius;
        model_options.sense = cfg.winding_sense;
        model_options.chunk = cfg.raw_chunk;
        rc = ScrollModel_build(arena, input, out_dir, &model_options, &model);
        qp_logf("[scroll model] %s, %.2fs\n", rc == 0 ? "built (physical qualification pending)" : "FAILED", ves_clock_sec()-t0);
        Arena_dispose(&arena);
        return rc == 0 ? 0 : 1;
    }
    qp_join(mesh_dir, sizeof mesh_dir, out_dir, "stage1_mesh");
    qp_join(mesh_vmesh, sizeof mesh_vmesh, mesh_dir, "mesh.vmesh");
    qp_join(mesh_table, sizeof mesh_table, mesh_dir, "mesh_cubes.tsv");
    {
        double ts = ves_clock_sec();
        if (qp_stage_mesh(arena, &cfg, input, is_vmesh, mesh_dir, mesh_vmesh) != 0)
            return 1;
        qp_logf("[stage 1 mesh] done (%.1fs)\n", ves_clock_sec() - ts);
    }
    if (stop_stage == 1) { qp_stop(1, "mesh", t0); Arena_dispose(&arena); return 0; }

    /* stage 2: unwrap -- initial winding, registration/MRF, gauge sync and
     * the arc-length U map, all inside scroll_ribbon --winding-only */
    qp_join(unwrap_dir, sizeof unwrap_dir, out_dir, "stage2_unwrap");
    qp_join(cert_vmesh, sizeof cert_vmesh, unwrap_dir, "unwrap_winding.vmesh");
    {
    int stage2_ready = MeshBin_looks_complete(cert_vmesh);
    char canonical_ledger[QP_MAX_PATH];
    char canonical_owner[QP_MAX_PATH];
    qp_join(canonical_ledger, sizeof canonical_ledger, unwrap_dir,
            "canonical_blocks.json");
    qp_join(canonical_owner, sizeof canonical_owner, unwrap_dir,
            "unwrap_winding_owner_block.i32");
    if (stage2_ready && g_qp_canonical_blocks &&
        !qp_exists(canonical_ledger)) {
        qp_logf("[stage 2 unwrap] existing certificate is monolithic but "
                "--canonical-blocks was requested; use a fresh output "
                "directory\n");
        return 1;
    }
    if (stage2_ready && g_qp_canonical_blocks &&
        !qp_file_has_text(canonical_ledger,
                          "vesuvius-canonical-block-certificate-v3")) {
        qp_logf("[stage 2 unwrap] canonical certificate predates the "
                "projective boundary contract; rebuilding\n");
        stage2_ready = 0;
    }
    if (stage2_ready && g_qp_canonical_blocks && !qp_exists(canonical_owner)) {
        qp_logf("[stage 2 unwrap] upgrading canonical certificate with "
                "immutable owner-block provenance\n");
        stage2_ready = 0;
    }
    if (!stage2_ready) {
        qp_mkdir(unwrap_dir);
        if (g_qp_canonical_blocks &&
            qp_stage_unwrap_canonical(&cfg, input, mesh_vmesh, mesh_table,
                                      unwrap_dir, cert_vmesh) != 0)
            return 1;
        if (!g_qp_canonical_blocks &&
            qp_spawn_scroll_ribbon(&cfg, mesh_vmesh, unwrap_dir, "unwrap", 0, 0,
                                   NULL, NULL, NULL, NULL,
                                   "scroll_ribbon(unwrap)") != 0)
            return 1;
        if (!MeshBin_looks_complete(cert_vmesh)) {
            qp_logf("[stage 2 unwrap] no certificate written: %s\n", cert_vmesh);
            return 1;
        }
    } else {
        qp_logf("[stage 2 unwrap] resume: already complete\n");
    }
    }
    qp_logf("[stage 2 unwrap] done\n");
    if (stop_stage == 2) { qp_stop(2, "unwrap", t0); Arena_dispose(&arena); return 0; }

    /* stage 3: fit (parameterize the lattice ribbon on the certificate) */
    qp_join(fit_dir, sizeof fit_dir, out_dir, "stage3_fit");
    qp_join(fit_stem, sizeof fit_stem, fit_dir, "fit_ribbon");
    snprintf(fit_vmesh, sizeof fit_vmesh, "%s.vmesh", fit_stem);
    qp_join(fit_verdict, sizeof fit_verdict, fit_dir, "verdict.json");
    {
    int stage3_ready = MeshBin_looks_complete(fit_vmesh);
    char fit_owner[QP_MAX_PATH];
    char fit_contract[QP_MAX_PATH];
    /* scroll_ribbon names carried certificate sidecars from its --id, not
     * from the primary output stem (fit_winding_*, versus fit_ribbon_*). */
    qp_join(fit_owner, sizeof fit_owner, fit_dir,
            "fit_winding_owner_block.i32");
    qp_join(fit_contract, sizeof fit_contract, fit_dir,
            "projective_contract.json");
    if (stage3_ready && g_qp_canonical_blocks && trust_gauge &&
        !qp_exists(fit_owner)) {
        qp_logf("[stage 3 fit] existing ribbon predates canonical projective "
                "ownership; rebuilding\n");
        stage3_ready = 0;
    }
    if (stage3_ready && g_qp_canonical_blocks && trust_gauge &&
        !qp_file_has_text(fit_contract,
                          "vesuvius-projective-fit-contract-v6")) {
        qp_logf("[stage 3 fit] existing ribbon predates the lossless global "
                "lineage partition and literal projective coordinates; "
                "rebuilding\n");
        stage3_ready = 0;
    }
    if (!stage3_ready) {
        qp_mkdir(fit_dir);
        if (qp_spawn_scroll_ribbon(&cfg, cert_vmesh, fit_dir, "fit", 1,
                                   trust_gauge, NULL, NULL, NULL, NULL,
                                   "scroll_ribbon(fit)") != 0)
            return 1;
        if (!MeshBin_looks_complete(fit_vmesh)) {
            qp_logf("[stage 3 fit] no ribbon written: %s\n", fit_vmesh);
            return 1;
        }
        if (g_qp_canonical_blocks && trust_gauge) {
            FILE *f = fopen(fit_contract, "wb");
            if (f == NULL) return 1;
            fprintf(f,
                    "{ \"schema\": \"vesuvius-projective-fit-contract-v6\", "
                    "\"coordinates\": "
                    "\"literal carried U and absolute canonical V lattice\", "
                    "\"claimants\": "
                    "\"owner-DAG causal row selection; only equal or ancestor "
                    "owners may provide predecessor evidence\", "
                    "\"collision_layers\": "
                    "\"lossless global lineage coloring; one VMESH storage "
                    "chart per color at unchanged absolute U/V\", "
                    "\"semantic_sample_key\": "
                    "\"absolute U/V and lifted-phase bits plus stable material "
                    "identity and claimant chart; storage color is non-semantic\", "
                    "\"forbidden_dependencies\": "
                    "\"crop minimum,crop-wide coordinate consensus,sibling or "
                    "descendant coordinate votes,residual winner deletion,"
                    "crop-height coordinate offsets\", "
                    "\"semantic_ids\": \"reversible block-local i31\", "
                    "\"reconstruction_output\": "
                    "\"stable material identity; branch lanes and storage "
                    "colors internal\" }\n");
            if (fclose(f) != 0) return 1;
        }
    } else {
        qp_logf("[stage 3 fit] resume: already complete\n");
    }
    }
    if (qp_publish_plain_sidecars(fit_dir, "fit_ribbon") != 0) return 1;
    qp_verdict(&cfg, "fit", fit_vmesh, fit_stem, cert_vmesh, fit_verdict);
    (void)qp_debug_bake(&cfg, fit_vmesh, out_dir,
                        "stage3_fit");  /* diagnostic: never fatal */
    qp_logf("[stage 3 fit] done\n");
    if (stop_stage == 3) { qp_stop(3, "fit", t0); Arena_dispose(&arena); return 0; }

stage4:
    /* stage 4: solidify the observations.  Canonical/projective runs consume
     * the complete stage-3 collision atlas; the layer-0 VMESH remains the
     * compatibility presentation input, never the semantic authority. */
    qp_join(solid_dir, sizeof solid_dir, out_dir, "stage4_solid");
    qp_join(solid_stem, sizeof solid_stem, solid_dir, "solid_ribbon");
    snprintf(solid_vmesh, sizeof solid_vmesh, "%s.vmesh", solid_stem);
    snprintf(solid_labels, sizeof solid_labels,
             "%s_reconstruction_component.i32", solid_stem);
    qp_join(solid_verdict, sizeof solid_verdict, solid_dir, "verdict.json");
    qp_mkdir(solid_dir);
    if (g_qp_canonical_blocks) {
        QpProjectiveAtlas atlas;
        memset(&atlas, 0, sizeof atlas);
        if (qp_stage_solid_projective_atlas(
                &cfg, fit_stem, solid_stem, "stage 4 solid atlas") != 0)
            return 1;
        if (qp_projective_atlas_discover(fit_stem, solid_stem, &atlas) != 0)
            return 1;
        for (size_t i = 0; i < atlas.n; i++) {
            char chart_tag[128], verdict[QP_MAX_PATH];
            if (qp_projective_layer_stem(chart_tag, sizeof chart_tag,
                                         "stage4_solid",
                                         atlas.charts[i].layer) != 0) {
                qp_projective_atlas_dispose(&atlas);
                return 1;
            }
            if (atlas.charts[i].layer == 0)
                snprintf(verdict, sizeof verdict, "%s", solid_verdict);
            else
                snprintf(verdict, sizeof verdict, "%s_verdict.json",
                         atlas.charts[i].output_stem);
            qp_verdict(&cfg, chart_tag, atlas.charts[i].output_vmesh,
                       atlas.charts[i].output_stem, cert_vmesh, verdict);
            (void)qp_debug_bake(&cfg, atlas.charts[i].output_vmesh, out_dir,
                                chart_tag); /* diagnostic: never fatal */
        }
        qp_projective_atlas_dispose(&atlas);
        qp_logf("[stage 4 atlas] layer-0 remains the legacy sheet/optimizer "
                "presentation; solid_ribbon_atlas.json is authoritative\n");
    } else {
        int solid_ready = QuadribbonSolidify_complete(solid_stem);
        if (!solid_ready) {
            if (qp_stage_solid(&cfg, fit_vmesh, solid_stem,
                               "stage 4 solid") != 0)
                return 1;
        } else {
            qp_logf("[stage 4 solid] resume: already complete\n");
        }
        qp_verdict(&cfg, "solid", solid_vmesh, solid_stem, cert_vmesh,
                   solid_verdict);
        (void)qp_debug_bake(&cfg, solid_vmesh, out_dir,
                            "stage4_solid"); /* diagnostic: never fatal */
    }
    qp_logf("[stage 4 solid] done\n");

    /* stage 6a: the "pre" sheet (solid ribbon) is emitted BEFORE optimize:
     * the deliverable exists as soon as the solid ribbon does, and
     * --stop-after solid yields a sheet.  The untangler on a large pile
     * (21x3x3 tube: 2.2M verts, 280k collisions) runs for hours, and on
     * the 5x3x3 core it measured coverage -9 pts / seams x2 for +10 pts of
     * coherence -- the pre sheet must never wait for it. */
    char sheet_pre[QP_MAX_PATH];
    uint8_t *tex_pre = NULL;
    size_t wpre = 0, hpre = 0;
    qp_join(sheet_pre, sizeof sheet_pre, out_dir, "stage6_sheet_pre");
    if (qp_stage_sheet(arena, &cfg, solid_vmesh, solid_labels, sheet_pre,
                       "pre", out_dir, &tex_pre, &wpre, &hpre) != 0) {
        qp_logf("[stage 6 sheet-pre] FAILED (see logs above)\n");
        return 1;
    }
    qp_logf("[stage 6 sheet-pre] done (before optimize)\n");
    if (stop_stage == 4) { qp_stop(4, "solid", t0); Arena_dispose(&arena); return 0; }

    /* stage 5: optimize (untangle -> refit rounds) */
    qp_join(opt_dir, sizeof opt_dir, out_dir, "stage5_optimize");
    qp_join(opt_stem, sizeof opt_stem, opt_dir, "optimized_ribbon");
    snprintf(opt_vmesh, sizeof opt_vmesh, "%s.vmesh", opt_stem);
    snprintf(opt_labels, sizeof opt_labels, "%s_reconstruction_component.i32",
             opt_stem);
    if (!QuadribbonSolidify_complete(opt_stem)) {
        qp_mkdir(opt_dir);
        if (qp_stage_optimize(&cfg, solid_stem, opt_dir, cert_vmesh, opt_stem,
                              sizeof opt_stem) != 0)
            return 1;
        snprintf(opt_vmesh, sizeof opt_vmesh, "%s.vmesh", opt_stem);
        snprintf(opt_labels, sizeof opt_labels,
                 "%s_reconstruction_component.i32", opt_stem);
    } else {
        qp_logf("[stage 5 optimize] resume: already complete\n");
    }
    (void)qp_debug_bake(&cfg, opt_vmesh, out_dir,
                        "stage5_optimize");  /* diagnostic: never fatal */
    qp_logf("[stage 5 optimize] done\n");
    if (stop_stage == 5) { qp_stop(5, "optimize", t0); Arena_dispose(&arena); return 0; }

    /* stage 6b: the "post" sheet (optimized ribbon) + pre-vs-post composite;
     * the pre sheet was emitted after stage 4 */
    {
        char sheet_post[QP_MAX_PATH];
        uint8_t *tex_post = NULL;
        size_t wpost = 0, hpost = 0;
        qp_join(sheet_post, sizeof sheet_post, out_dir, "stage6_sheet_post");
        if (qp_stage_sheet(arena, &cfg, opt_vmesh, opt_labels, sheet_post,
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
