/* sheet_assemble.c -- entry point of the chart-assembly unwrapper.
 *
 *   sheet_assemble.exe <pile_dir> <out_dir> [config.json]
 *                      [--stop-after clean|relate|pose|clean2|place|discover|repair|audit|sheet]
 *   sheet_assemble.exe --selftest
 *
 * Stages (each resumable from its artifacts under <out_dir>):
 *   clean   per-cube charts, proactive cleaning, intrinsic flattening
 *           -> <out>/charts/<cube_id>.asc, stage1_charts.png
 *   axis    the scroll axis derived from the material (asm_axis_derive):
 *           stage1_axis.csv / .png / .json; frames pose, placement, verdict
 *   relate  seam correspondences + layer neighbours (asm_relate)
 *   pose    rigid pose graph (asm_pose)
 *   clean2  contradiction audit + cleaning rounds (asm_conflict, asm_cut)
 *   place   disconnected components by measured turn vectors (asm_place)
 *   repair  coupled original-metric / physical-seam / material-area repair (C)
 *   audit   complete source-face, metric, continuity and triangle-domain audit (C)
 *   sheet   precise uv mesh, encoded readback audit, bakes, verdict
 *
 * This file is ONLY the entry point: argument parsing, config, stage chaining
 * and timing.  Algorithms live in src/assemble/. */
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../common/ves_omp.h"

#include "../assemble/asm_clean.h"
#include "../assemble/asm_conflict.h"
#include "../assemble/asm_emit.h"
#include "../assemble/asm_place.h"
#include "../assemble/asm_flatten.h"
#include "../assemble/asm_layer_cut.h"
#include "../assemble/asm_pose.h"
#include "../assemble/asm_relate.h"
#include "../assemble/asm_report.h"
#include "../assemble/asm_axis.h"
#include "../assemble/asm_axis_derive.h"
#include "../assemble/asm_reading_order.h"
#include "../assemble/asm_store.h"
#include "../assemble/asm_types.h"
#include "../assemble/asm_verdict.h"
#include "../assemble/asm_repair.h"
#include "../assemble/asm_audit.h"
#include "../assemble/asm_contacts.h"
#include "../assemble/asm_discover.h"
#include "../flatten/active_set_qp.h"
#include "../flatten/sparse_solve.h"
#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_pile.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#define SA_MAX_CUBES 32768

typedef struct SaConfig {
    int          threads;
    int          resume;
    AsmCleanOpts clean;
    AsmRelateOpts relate;
    AsmPoseOpts  pose;
    AsmConflictOpts conflict;
    AsmPlaceOpts placeo;
    AsmEmitOpts  emit;
    AsmRepairOpts repair;
    AsmDiscoverOpts discover;
    char         raw_source[2048];
    char         bake_exe[2048];
    int          atlas_cell_px;
    size_t       atlas_max_cells;
    AsmVerdictOpts verdict;
    char         verdict_exe[2048];
    long         pile_bbox[6];   /* z0 z1 y0 y1 x0 x1 of cube ORIGINS to keep; all zero = no filter */
    char         axis_table[1024];   /* optional scroll axis table (z,y,x per z): the pose prior follows its direction over the bbox */
    char         verdict_axis_table[1024];   /* optional: the verdict measures in this table's local axis frame */
    int          axis_derive;        /* stage `axis`: 1 derive the axis from the material, 0 never, -1 auto (when no pose.axis_table) */
    char         axis_prior_table[1024];     /* optional prior for the derivation (default: the pose table / umbilicus line) */
    double       axis_slab_vox;
    AsmAxis     *axis_derived;       /* the derived polyline once the stage ran */
    AsmAxis     *axis_pose;      /* the loaded tables (or the umbilicus line): pose prior + placement, verdict */
    AsmAxis     *axis_verdict;
    int          have_bbox;
} SaConfig;

static void sa_usage(void)
{
    fprintf(stderr,
        "usage: sheet_assemble <pile_dir> <out_dir> [config.json]\n"
        "                      [--stop-after clean|axis|relate|pose|clean2|place|discover|repair|audit|sheet]\n"
        "                      [--resume-placement stage5_state.asr]\n"
        "       sheet_assemble --review checkpoint.asr output_dir [axis.csv]\n"
        "       sheet_assemble --winding-order checkpoint.asr axis.csv output_dir\n"
        "       sheet_assemble --selftest\n"
        "       sheet_assemble --help\n"
        "Default configuration: configs/default.json (paths relative to the working directory).\n");
}

static void sa_config_defaults(SaConfig *c)
{
    memset(c, 0, sizeof *c);
    c->threads = ves_cpu_count();
    if (c->threads < 1) c->threads = 1;
    c->resume = 1;
    AsmClean_default_opts(&c->clean);
    AsmRelate_default_opts(&c->relate);
    AsmPose_default_opts(&c->pose);
    AsmConflict_default_opts(&c->conflict);
    AsmPlace_default_opts(&c->placeo);
    AsmEmit_default_opts(&c->emit);
    AsmRepair_defaults(&c->repair);
    AsmDiscover_default_opts(&c->discover);
    AsmVerdict_default_opts(&c->verdict);
    snprintf(c->verdict_exe, sizeof c->verdict_exe, "%s", c->verdict.verdict_exe);
    c->verdict.verdict_exe = c->verdict_exe;
    c->raw_source[0] = 0;
    snprintf(c->bake_exe, sizeof c->bake_exe, "%s", c->emit.bake_exe);
    c->emit.bake_exe = c->bake_exe;
    c->atlas_cell_px = 96;
    c->atlas_max_cells = 4096;
    c->axis_derive = -1;
    c->axis_slab_vox = ASM_AXIS_DERIVE_SLAB_VOX;
}

static int sa_config_load(Arena_T arena, const char *path, SaConfig *c)
{
    const char *err = NULL;
    const JsonValue *root = Json_parse_file(arena, path, &err);
    if (!root) {
        fprintf(stderr, "sheet_assemble: cannot read config %s: %s\n", path, err ? err : "?");
        return -1;
    }
    const JsonValue *geometry = Json_object_get(root, "geometry");
    if (geometry) {
        const JsonValue *axis = Json_object_get(geometry, "cylindrical_axis_yx");
        if (axis && Json_array_len(axis) == 2) {
            c->verdict.umb_y = Json_as_double(Json_array_get(axis, 0), c->verdict.umb_y);
            c->verdict.umb_x = Json_as_double(Json_array_get(axis, 1), c->verdict.umb_x);
        }
        c->verdict.pitch = Json_member_double(geometry, "wrap_spacing_voxels", c->verdict.pitch);
        const char *table = Json_as_string(Json_object_get(geometry, "axis_table"));
        if (table) {
            snprintf(c->axis_table, sizeof c->axis_table, "%s", table);
            snprintf(c->verdict_axis_table, sizeof c->verdict_axis_table, "%s", table);
        }
    }
    const JsonValue *compute = Json_object_get(root, "compute");
    if (compute) c->threads = (int)Json_member_long(compute, "threads", c->threads);
    c->resume = Json_as_bool(Json_object_get(root, "resume"), c->resume);
    const JsonValue *repair = Json_object_get(root, "repair");
    if (repair) {
        c->repair.iterations = (int)Json_member_long(repair,"iterations",c->repair.iterations);
        c->repair.regional = (int)Json_member_long(repair,"regional",c->repair.regional);
        c->repair.admission = (int)Json_member_long(repair,"admission",c->repair.admission);
        c->repair.admission_rounds = (int)Json_member_long(repair,"admission_rounds",c->repair.admission_rounds);
        c->repair.closure_iterations = (int)Json_member_long(repair,"closure_iterations",c->repair.closure_iterations);
        c->repair.context_hops = (int)Json_member_long(repair,"context_hops",c->repair.context_hops);
        c->repair.compact = (int)Json_member_long(repair,"compact",c->repair.compact);
        c->repair.connectivity_first = (int)Json_member_long(repair,"connectivity_first",c->repair.connectivity_first);
        c->repair.patch_faces = (size_t)Json_member_long(repair,"patch_faces",(long)c->repair.patch_faces);
        c->repair.patch_field_faces = (size_t)Json_member_long(repair,"patch_field_faces",(long)c->repair.patch_field_faces);
        c->repair.patch_coordinates = (size_t)Json_member_long(repair,"patch_coordinates",(long)c->repair.patch_coordinates);
        c->repair.budget_seconds = Json_member_double(repair,"budget_seconds",c->repair.budget_seconds);
        c->repair.patch_seconds = Json_member_double(repair,"patch_seconds",c->repair.patch_seconds);
        c->repair.clearance_iterations = (int)Json_member_long(repair,"clearance_iterations",c->repair.clearance_iterations);
        c->repair.source_trials = (int)Json_member_long(repair,"source_trials",c->repair.source_trials);
        c->repair.metric_weight = Json_member_double(repair,"metric_weight",c->repair.metric_weight);
        c->repair.seam_weight = Json_member_double(repair,"seam_weight",c->repair.seam_weight);
        c->repair.contact_weight = Json_member_double(repair,"contact_weight",c->repair.contact_weight);
        c->repair.overlap_weight = Json_member_double(repair,"overlap_weight",c->repair.overlap_weight);
        c->repair.contact_length = Json_member_double(repair,"contact_length",c->repair.contact_length);
        c->repair.proximal = Json_member_double(repair,"proximal",c->repair.proximal);
        c->repair.proximal_length = Json_member_double(repair,"proximal_length",c->repair.proximal_length);
    }
    const JsonValue *dis = Json_object_get(root, "discover");
    if (dis) {
        c->discover.enabled = (int)Json_member_long(dis, "enabled", c->discover.enabled);
        c->discover.min_sides = (int)Json_member_long(dis, "min_sides", c->discover.min_sides);
        c->discover.soft_min_sides = (int)Json_member_long(dis, "soft_min_sides", c->discover.soft_min_sides);
        c->discover.rounds = (int)Json_member_long(dis, "rounds", c->discover.rounds);
        c->discover.selected_obligations = (int)Json_member_long(dis, "selected_obligations", c->discover.selected_obligations);
    }
    const JsonValue *axb = Json_object_get(root, "axis");
    if (axb) {
        c->axis_derive = Json_as_bool(Json_object_get(axb, "derive"), c->axis_derive);
        c->axis_slab_vox = Json_member_double(axb, "slab_vox", c->axis_slab_vox);
        const JsonValue *pt = Json_object_get(axb, "prior_table");
        if (pt && Json_as_string(pt)) { strncpy(c->axis_prior_table, Json_as_string(pt), sizeof(c->axis_prior_table) - 1); c->axis_prior_table[sizeof(c->axis_prior_table) - 1] = 0; }
    }
    const JsonValue *cl = Json_object_get(root, "clean");
    if (cl) {
        c->clean.blob_area_per_face = Json_member_double(cl, "blob_area_per_face", c->clean.blob_area_per_face);
        c->clean.blob_double_sided = Json_member_double(cl, "blob_double_sided", c->clean.blob_double_sided);
        c->clean.min_chart_area = Json_member_double(cl, "min_chart_area", c->clean.min_chart_area);
        c->clean.stress_band = Json_member_double(cl, "stress_band", c->clean.stress_band);
        c->clean.flatten_iters = (int)Json_member_long(cl, "flatten_iters", c->clean.flatten_iters);
    }
    const JsonValue *rl = Json_object_get(root, "relate");
    if (rl) {
        c->relate.cube_size = Json_member_double(rl, "cube_size", c->relate.cube_size);
        c->relate.seam_max_vox = Json_member_double(rl, "seam_max_vox", c->relate.seam_max_vox);
        c->relate.seam_band_vox = Json_member_double(rl, "seam_band_vox", c->relate.seam_band_vox);
        c->relate.seam_min_len = Json_member_double(rl, "seam_min_len", c->relate.seam_min_len);
        c->relate.normal_dot_min = Json_member_double(rl, "normal_dot_min", c->relate.normal_dot_min);
        c->relate.tangency_max_vox = Json_member_double(rl, "tangency_max_vox", c->relate.tangency_max_vox);
        c->relate.iso_tol = Json_member_double(rl, "iso_tol", c->relate.iso_tol);
        c->relate.seam_noise_vox = Json_member_double(rl, "seam_noise_vox", c->relate.seam_noise_vox);
        c->relate.layer_probe_vox = Json_member_double(rl, "layer_probe_vox", c->relate.layer_probe_vox);
        c->relate.layer_samples = (int)Json_member_long(rl, "layer_samples", c->relate.layer_samples);
    }
    const JsonValue *po = Json_object_get(root, "pose");
    if (po) {
        c->pose.gn_iters = (int)Json_member_long(po, "gn_iters", c->pose.gn_iters);
        c->pose.cauchy_c = Json_member_double(po, "cauchy_c", c->pose.cauchy_c);
        c->pose.switch_off = Json_member_double(po, "switch_off", c->pose.switch_off);
        c->pose.damping = Json_member_double(po, "damping", c->pose.damping);
        const JsonValue *ax = Json_object_get(po, "axis_table");
        if (ax && Json_as_string(ax)) { strncpy(c->axis_table, Json_as_string(ax), sizeof(c->axis_table) - 1); c->axis_table[sizeof(c->axis_table) - 1] = 0; }
    }
    const JsonValue *cf = Json_object_get(root, "conflict");
    if (cf) {
        c->conflict.cell = Json_member_double(cf, "cell", c->conflict.cell);
        c->conflict.contra_lo = Json_member_double(cf, "contra_lo", c->conflict.contra_lo);
        c->conflict.contra_hi = Json_member_double(cf, "contra_hi", c->conflict.contra_hi);
        c->conflict.rounds = (int)Json_member_long(cf, "rounds", c->conflict.rounds);
        c->conflict.max_path = (int)Json_member_long(cf, "max_path", c->conflict.max_path);
        c->conflict.drop_frac = Json_member_double(cf, "drop_frac", c->conflict.drop_frac);
        c->conflict.extras_ratio = Json_member_double(cf, "extras_ratio", c->conflict.extras_ratio);
    }
    const JsonValue *pl = Json_object_get(root, "place");
    if (pl) {
        c->placeo.min_anchors = (int)Json_member_long(pl, "min_anchors", c->placeo.min_anchors);
        c->placeo.max_rms = Json_member_double(pl, "max_rms", c->placeo.max_rms);
        c->placeo.rel_rms = Json_member_double(pl, "rel_rms", c->placeo.rel_rms);
        c->placeo.max_layers = (int)Json_member_long(pl, "max_layers", c->placeo.max_layers);
        c->placeo.layer_tol = Json_member_double(pl, "layer_tol", c->placeo.layer_tol);
        c->placeo.min_turn_hits = (int)Json_member_long(pl, "min_turn_hits", c->placeo.min_turn_hits);
        c->placeo.max_contra_ratio = Json_member_double(pl, "max_contra_ratio", c->placeo.max_contra_ratio);
        c->placeo.max_rot_refine = Json_member_double(pl, "max_rot_refine", c->placeo.max_rot_refine);
        c->placeo.max_rim_ratio = Json_member_double(pl, "max_rim_ratio", c->placeo.max_rim_ratio);
        c->placeo.crosswrap_anchor_veto = (int)Json_member_double(pl, "crosswrap_anchor_veto", (double)c->placeo.crosswrap_anchor_veto);
    }
    const JsonValue *rs = Json_object_get(root, "raw_source");
    if (rs) {
        const char *rp0 = Json_as_string(Json_object_get(rs, "path"));
        if (rp0) { snprintf(c->raw_source, sizeof c->raw_source, "%s", rp0); c->emit.raw_source = c->raw_source; }
    }
    const JsonValue *bk = Json_object_get(root, "bake");
    if (bk) {
        const char *exe = Json_as_string(Json_object_get(bk, "exe"));
        if (exe) { snprintf(c->bake_exe, sizeof c->bake_exe, "%s", exe); c->emit.bake_exe = c->bake_exe; }
        c->emit.normal_range = Json_member_double(bk, "normal_reach_vox", c->emit.normal_range);
        c->emit.raster_du = Json_member_double(bk, "raster_du", c->emit.raster_du);
        c->emit.raster_dv = Json_member_double(bk, "raster_dv", c->emit.raster_dv);
        c->emit.level_v_to_z = Json_as_bool(Json_object_get(bk, "level_v_to_z"), c->emit.level_v_to_z);
    }
    const JsonValue *pb = Json_object_get(root, "pile_bbox");
    if (pb && Json_type(pb) == JSON_ARRAY && Json_array_len(pb) == 6) {
        for (size_t k = 0; k < 6; k++) c->pile_bbox[k] = Json_as_long(Json_array_get(pb, k), 0);
        c->have_bbox = 1;
    }
    const JsonValue *vd = Json_object_get(root, "verdict");
    if (vd) {
        const char *exe = Json_as_string(Json_object_get(vd, "exe"));
        if (exe) { snprintf(c->verdict_exe, sizeof c->verdict_exe, "%s", exe); c->verdict.verdict_exe = c->verdict_exe; }
        c->verdict.du = Json_member_double(vd, "du", c->verdict.du);
        c->verdict.dv = Json_member_double(vd, "dv", c->verdict.dv);
        c->verdict.umb_y = Json_member_double(vd, "umb_y", c->verdict.umb_y);
        c->verdict.umb_x = Json_member_double(vd, "umb_x", c->verdict.umb_x);
        c->verdict.pitch = Json_member_double(vd, "pitch", c->verdict.pitch);
        c->verdict.core_radius = Json_member_double(vd, "core_radius", c->verdict.core_radius);
        c->verdict.include_blobs_in_source = Json_as_bool(Json_object_get(vd, "include_blobs_in_source"), 0);
        const char *vax = Json_as_string(Json_object_get(vd, "axis_table"));
        if (vax) { strncpy(c->verdict_axis_table, vax, sizeof(c->verdict_axis_table) - 1); c->verdict_axis_table[sizeof(c->verdict_axis_table) - 1] = 0; }
    }
    const JsonValue *rp = Json_object_get(root, "report");
    if (rp) {
        c->atlas_cell_px = (int)Json_member_long(rp, "atlas_cell_px", c->atlas_cell_px);
        c->atlas_max_cells = (size_t)Json_member_long(rp, "atlas_max_cells", (long)c->atlas_max_cells);
    }
    if (c->threads < 1) c->threads = 1;
    AsmField_threads(c->threads);
    return 0;
}

static int sa_bbox_filter(void *ctx, long oz, long oy, long ox)
{
    const SaConfig *c = ctx;
    return oz >= c->pile_bbox[0] && oz < c->pile_bbox[1] && oy >= c->pile_bbox[2] && oy < c->pile_bbox[3] &&
           ox >= c->pile_bbox[4] && ox < c->pile_bbox[5];
}

/* ---- stage 1: clean ------------------------------------------------------- */

typedef struct SaCubeResult {
    AsmChart *charts;
    size_t    n;
    int       from_store;
    int       failed;
} SaCubeResult;

static int sa_stage_clean(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    double t0 = ves_clock_sec();
    char store_dir[2048];
    snprintf(store_dir, sizeof store_dir, "%s/charts", out_dir);
    ves_mkdir(store_dir);

    int threads = cfg->threads;
    Arena_T *scratch = ARENA_ALLOC(run->arena, (size_t)threads * sizeof(Arena_T));
    Arena_T *persist = ARENA_ALLOC(run->arena, (size_t)threads * sizeof(Arena_T));
    AsmCleanStats *tstats = ARENA_CALLOC(run->arena, (size_t)threads, sizeof(AsmCleanStats));
    for (int t = 0; t < threads; t++) { scratch[t] = Arena_new(); persist[t] = Arena_new(); }
    SaCubeResult *res = ARENA_CALLOC(run->arena, run->n_cubes, sizeof(SaCubeResult));
    size_t n_cubes = run->n_cubes;
    int n_cubes_i = (int)n_cubes;
    size_t reused = 0, failed = 0, stale = 0;
    int ci = 0;
    /* the store identity: the cleaning options + the code version (the source mesh's
     * path / size / mtime are checked per cube) */
    uint64_t fingerprint = 0;
    {
        char text[512];
        snprintf(text, sizeof text, "%s|blob_area_per_face=%.9g|blob_double_sided=%.9g|min_chart_area=%.9g|stress_band=%.9g|flatten_iters=%d",
                 ASM_STORE_CODE_VERSION, cfg->clean.blob_area_per_face, cfg->clean.blob_double_sided, cfg->clean.min_chart_area,
                 cfg->clean.stress_band, cfg->clean.flatten_iters);
        fingerprint = AsmStore_fingerprint(text);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) reduction(+:reused, failed, stale)
#endif
    for (ci = 0; ci < n_cubes_i; ci++) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        const MeshPileEntry *e = &run->pile[ci];
        SaCubeResult *r = &res[ci];
        if (cfg->resume && AsmStore_exists(store_dir, e->cube_id)) {
            int src = AsmStore_read_cube(persist[tid], store_dir, e->cube_id, e->path, fingerprint, &r->charts, &r->n);
            if (src == 0) {
                r->from_store = 1;
                reused++;
                continue;
            }
            if (src == -2) stale++;   /* another source file or another cleaning: re-clean */
        }
        Arena_Mark mark = Arena_save(scratch[tid]);
        MeshBinData mesh;
        if (MeshBin_read_arena(scratch[tid], e->path, &mesh) != 0) {
            fprintf(stderr, "  [clean] cube %s: cannot read %s\n", e->cube_id, e->path);
            r->failed = 1;
            failed++;
            Arena_restore(scratch[tid], mark);
            continue;
        }
        AsmCleanStats st;
        AsmClean_cube(scratch[tid], persist[tid], &mesh, ci, &cfg->clean, &r->charts, &r->n, &st);
        Arena_restore(scratch[tid], mark);
        tstats[tid].components += st.components;
        tstats[tid].blobs += st.blobs; tstats[tid].tiny += st.tiny;
        tstats[tid].nonmanifold += st.nonmanifold; tstats[tid].handles += st.handles;
        tstats[tid].flat_failed += st.flat_failed; tstats[tid].suspect += st.suspect;
        tstats[tid].kept += st.kept;
        tstats[tid].verts_in += st.verts_in; tstats[tid].faces_in += st.faces_in;
        tstats[tid].verts_kept += st.verts_kept; tstats[tid].faces_kept += st.faces_kept;
        tstats[tid].verts_blob += st.verts_blob; tstats[tid].faces_blob += st.faces_blob;
        tstats[tid].area_in += st.area_in; tstats[tid].area_kept += st.area_kept;
        tstats[tid].area_blob += st.area_blob; tstats[tid].flatten_sec += st.flatten_sec;
        tstats[tid].degenerate_split += st.degenerate_split; tstats[tid].degenerate_faces += st.degenerate_faces;
        tstats[tid].crumple_split += st.crumple_split; tstats[tid].crumple_pieces += st.crumple_pieces;
        tstats[tid].crumple_faces += st.crumple_faces; tstats[tid].crumple_area_excluded += st.crumple_area_excluded;
        if (AsmStore_write_cube(store_dir, e->cube_id, e->path, fingerprint, r->charts, r->n) != 0)
            fprintf(stderr, "  [clean] cube %s: store write failed\n", e->cube_id);
    }

    /* merge in pile order with stable ids */
    AsmCleanStats tot;
    memset(&tot, 0, sizeof tot);
    for (int t = 0; t < threads; t++) {
        tot.components += tstats[t].components; tot.blobs += tstats[t].blobs; tot.tiny += tstats[t].tiny;
        tot.nonmanifold += tstats[t].nonmanifold; tot.handles += tstats[t].handles;
        tot.flat_failed += tstats[t].flat_failed; tot.suspect += tstats[t].suspect; tot.kept += tstats[t].kept;
        tot.verts_in += tstats[t].verts_in; tot.faces_in += tstats[t].faces_in;
        tot.verts_kept += tstats[t].verts_kept; tot.faces_kept += tstats[t].faces_kept;
        tot.verts_blob += tstats[t].verts_blob; tot.faces_blob += tstats[t].faces_blob;
        tot.area_in += tstats[t].area_in; tot.area_kept += tstats[t].area_kept;
        tot.area_blob += tstats[t].area_blob; tot.flatten_sec += tstats[t].flatten_sec;
        tot.degenerate_split += tstats[t].degenerate_split; tot.degenerate_faces += tstats[t].degenerate_faces;
        tot.crumple_split += tstats[t].crumple_split; tot.crumple_pieces += tstats[t].crumple_pieces;
        tot.crumple_faces += tstats[t].crumple_faces; tot.crumple_area_excluded += tstats[t].crumple_area_excluded;
    }
    size_t stored_kept = 0, stored_charts = 0;
    for (size_t cc = 0; cc < n_cubes; cc++) {
        for (size_t k = 0; k < res[cc].n; k++) {
            AsmChart c = res[cc].charts[k];
            c.id = (int32_t)run->n_charts;
            c.cube = (int32_t)cc;
            AsmRun_push_chart(run, &c);
            if (res[cc].from_store) {
                stored_charts++;
                if (AsmChart_in_layout(&c)) stored_kept++;
            }
        }
    }
    for (int t = 0; t < threads; t++) Arena_dispose(&scratch[t]);
    /* persist arenas hold the chart arrays for the run: keep them alive */

    size_t in_layout = 0, verts_layout = 0, faces_layout = 0;
    for (size_t i = 0; i < run->n_charts; i++) {
        if (!AsmChart_in_layout(&run->charts[i])) continue;
        in_layout++; verts_layout += run->charts[i].nv; faces_layout += run->charts[i].nf;
    }
    double wall = ves_clock_sec() - t0;
    fprintf(stderr,
        "[assemble clean] cubes %zu (reused %zu, stale re-cleaned %zu, failed %zu) charts %zu in_layout %zu | fresh: comps %zu kept %zu "
        "blobs %zu tiny %zu nonmanifold %zu handles %zu flat_failed %zu suspect %zu | "
        "verts in %zu kept %zu blob %zu | area in %.3e kept %.3e blob %.3e | flatten %.1f cpu-s | %.1f s\n",
        n_cubes, reused, stale, failed, run->n_charts, in_layout,
        tot.components, tot.kept, tot.blobs, tot.tiny, tot.nonmanifold, tot.handles, tot.flat_failed, tot.suspect,
        tot.verts_in, tot.verts_kept, tot.verts_blob, tot.area_in, tot.area_kept, tot.area_blob,
        tot.flatten_sec, wall);
    fprintf(stderr, "[assemble clean] zero-area faces: %zu charts flattened without them, %zu faces in excluded remainder charts"
                    " (fresh cubes)\n", tot.degenerate_split, tot.degenerate_faces);
    fprintf(stderr, "[assemble clean] crumpled charts: %zu rescued as %zu certified pieces; %zu faces (%.1f vox^2) in excluded remainder charts"
                    " (fresh cubes)\n", tot.crumple_split, tot.crumple_pieces, tot.crumple_faces, tot.crumple_area_excluded);
    fprintf(stderr, "[assemble clean] layout verts %zu faces %zu (stored charts %zu, stored in layout %zu)\n",
            verts_layout, faces_layout, stored_charts, stored_kept);

    char png[2048];
    snprintf(png, sizeof png, "%s/stage1_charts.png", out_dir);
    if (AsmReport_charts_atlas_png(run->arena, png, run->charts, run->n_charts,
                                   cfg->atlas_cell_px, cfg->atlas_max_cells) == 0)
        fprintf(stderr, "[assemble clean] wrote: %s\n", png);

    char rep[2048];
    snprintf(rep, sizeof rep, "%s/stage1_clean.json", out_dir);
    FILE *fp = fopen(rep, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"clean\",\n  \"cubes\": %zu,\n  \"cubes_reused\": %zu,\n  \"cubes_failed\": %zu,\n"
                    "  \"charts\": %zu,\n  \"charts_in_layout\": %zu,\n  \"layout_verts\": %zu,\n  \"layout_faces\": %zu,\n"
                    "  \"fresh\": {\"components\": %zu, \"kept\": %zu, \"blobs\": %zu, \"tiny\": %zu, \"nonmanifold\": %zu,"
                    " \"handles\": %zu, \"flat_failed\": %zu, \"suspect\": %zu,\n"
                    "            \"verts_in\": %zu, \"verts_kept\": %zu, \"verts_blob\": %zu,"
                    " \"area_in\": %.6e, \"area_kept\": %.6e, \"area_blob\": %.6e, \"flatten_cpu_sec\": %.3f},\n"
                    "  \"wall_sec\": %.3f\n}\n",
                n_cubes, reused, failed, run->n_charts, in_layout, verts_layout, faces_layout,
                tot.components, tot.kept, tot.blobs, tot.tiny, tot.nonmanifold, tot.handles, tot.flat_failed, tot.suspect,
                tot.verts_in, tot.verts_kept, tot.verts_blob, tot.area_in, tot.area_kept, tot.area_blob,
                tot.flatten_sec, wall);
        fclose(fp);
    }
    return 0;
}

/* ---- stage 2: relate --------------------------------------------------------- */

static int sa_stage_relate(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    AsmRelateStats st;
    char csv[2048];
    snprintf(csv, sizeof csv, "%s/stage2_pairs.csv", out_dir);
    AsmRelateOpts ro = cfg->relate;
    ro.diag_csv = csv;
    if (AsmRelate_run(run, &ro, cfg->threads, &st) != 0) return -1;
    /* Preserve the individual source observations used by placement, not
     * just the pair medians. Indices refer to the original clean-chart
     * vertices/faces; barycentric coordinates reconstruct the ray hit. */
    snprintf(csv,sizeof csv,"%s/stage2_layer_hits.csv",out_dir);
    FILE *hits=fopen(csv,"wb");
    if (!hits) return -1;
    fputs("layer_pair,hit,chart_a,vertex_a,chart_b,face_b,l0,l1,distance,side,order,count,d_median,pair_order,nsign,pair_nsign\n",hits);
    for (size_t p=0; p<run->n_layers; p++) {
        const AsmLayerPair *pair=&run->layers[p];
        for (int32_t k=0; k<pair->hit_count; k++) {
            size_t index=(size_t)pair->hit_first+(size_t)k;
            const AsmLayerHit *h=&run->layer_hits[index];
            fprintf(hits,"%zu,%zu,%d,%d,%d,%d,%.9g,%.9g,%.9g,%d,%d,%d,%.17g,%d,%d,%d\n",
                    p,index,h->a,h->va,h->b,h->fb,h->l0,h->l1,h->dist,(int)h->side,(int)h->order,
                    pair->count,pair->d_median,(int)pair->k,(int)h->nsign,(int)pair->nsign);
        }
    }
    int hit_error=ferror(hits);
    if (fclose(hits) || hit_error) return -1;
    fprintf(stderr,
        "[assemble relate] cube pairs %zu candidates %zu chart pairs %zu accepted %zu (parity-flipped %zu) corr %zu | "
        "rejected: len %zu tangency %zu iso %zu unimodal %zu fit %zu | kept for placement: weak %zu short %zu | rms p50 %.2f p95 %.2f vox | %.1f s\n",
        st.cube_pairs, st.candidates, st.chart_pairs, st.accepted, st.parity_flipped, st.corr_total,
        st.rej_len, st.rej_tangency, st.rej_iso, st.rej_unimodal, st.rej_fit, st.weak_kept, st.short_kept, st.rms_p50, st.rms_p95, st.seam_sec);
    fprintf(stderr, "[assemble relate] seam gap histogram of mutual-nearest boundary pairs (1-vox bins, normals agreeing): ");
    for (int h = 0; h < ASM_RELATE_HIST; h++) fprintf(stderr, "%zu ", st.gap_hist[h]);
    fprintf(stderr, "(pairs %zu, gate %.2f vox)\n", st.gap_pairs, cfg->relate.seam_max_vox);
    fprintf(stderr, "[assemble relate] mean per-pair median seam discrepancy %.3f over %zu pairs (gate %.3f)\n",
            st.iso_disc_n ? st.iso_disc_sum / (double)st.iso_disc_n : 0.0, st.iso_disc_n, cfg->relate.iso_tol);
    fprintf(stderr, "[assemble relate] layer pairs %zu charts with a layer distance %zu / %zu, global layer distance %.2f vox | %.1f s\n",
            st.layer_pairs, st.charts_with_layer, run->n_charts, st.layer_d_global, st.layer_sec);
    fprintf(stderr, "[assemble relate] layer hit histogram (2-vox bins from 0): ");
    for (int h = 0; h < ASM_RELATE_HIST; h++) fprintf(stderr, "%zu ", st.layer_hist[h]);
    fprintf(stderr, "(hits %zu)\n", st.layer_hits);
    char rep[2048];
    snprintf(rep, sizeof rep, "%s/stage2_relate.json", out_dir);
    FILE *fp = fopen(rep, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"relate\",\n  \"cube_pairs\": %zu,\n  \"candidates\": %zu,\n  \"chart_pairs\": %zu,\n"
                    "  \"accepted\": %zu,\n  \"parity_flipped\": %zu,\n  \"correspondences\": %zu,\n"
                    "  \"rejected\": {\"length\": %zu, \"tangency\": %zu, \"isometry\": %zu, \"unimodal\": %zu, \"fit\": %zu},\n"
                    "  \"kept_for_placement\": {\"weak\": %zu, \"short\": %zu, \"short_rejected_corr\": %zu, \"short_rejected_rms\": %zu},\n"
                    "  \"rms_p50\": %.4f,\n  \"rms_p95\": %.4f,\n  \"layer_pairs\": %zu,\n  \"charts_with_layer\": %zu,\n"
                    "  \"layer_d_global\": %.4f,\n  \"seam_sec\": %.3f,\n  \"layer_sec\": %.3f\n}\n",
                st.cube_pairs, st.candidates, st.chart_pairs, st.accepted, st.parity_flipped, st.corr_total,
                st.rej_len, st.rej_tangency, st.rej_iso, st.rej_unimodal, st.rej_fit, st.weak_kept, st.short_kept, st.short_rej_corr, st.short_rej_rms, st.rms_p50, st.rms_p95,
                st.layer_pairs, st.charts_with_layer, st.layer_d_global, st.seam_sec, st.layer_sec);
        fclose(fp);
    }
    return 0;
}

/* ---- stage 3: pose ----------------------------------------------------------- */


/* The scroll axis direction over the pile's z range from an axis table (rows z,y,x; '#' comments):
 * a straight-line fit of y(z) and x(z), normalized in (z,y,x).  World z when no table is given. */
/* The axis for the run: geometry.axis_table, optional legacy pose/verdict overrides,
 * or the configured umbilicus as a straight line along +z; each consumer projects its own points
 * onto the polyline (asm_axis.h).  Prints the axis' tilt over the pile bbox. */
static void sa_axis_setup(const AsmRun *run, SaConfig *c)
{
    Arena_T arena = run->arena;
    double z0 = c->have_bbox ? (double)c->pile_bbox[0] : 0.0, z1 = c->have_bbox ? (double)c->pile_bbox[1] + 128.0 : 0.0;
    /* A reusable config has no crop filter. Measure its axis over the actual
     * input cubes instead of the unrelated default interval around z=0. */
    if (!c->have_bbox) {
        int found = 0;
        for (size_t i = 0; i < run->n_cubes; i++) if (run->pile[i].has_id) {
            double lo = (double)run->pile[i].oz, hi = lo + c->relate.cube_size;
            if (!found || lo < z0) z0 = lo;
            if (!found || hi > z1) z1 = hi;
            found = 1;
        }
    }
    c->axis_pose = AsmAxis_load(arena, c->axis_table);
    if (c->axis_table[0] && c->axis_pose == NULL) fprintf(stderr, "[assemble] axis table %s: cannot read, the pose prior follows world z\n", c->axis_table);
    c->axis_verdict = c->verdict_axis_table[0] ? AsmAxis_load(arena, c->verdict_axis_table) : NULL;
    if (c->verdict_axis_table[0] && c->axis_verdict == NULL) fprintf(stderr, "[assemble] axis table %s: cannot read, the verdict measures about the umbilicus\n", c->verdict_axis_table);
    if (c->axis_verdict == NULL && c->verdict.umb_y > 0.0 && c->verdict.umb_x > 0.0) {
        double pt[3] = { 0.0, c->verdict.umb_y, c->verdict.umb_x }, dz[3] = { 1.0, 0.0, 0.0 };
        c->axis_verdict = AsmAxis_line(arena, pt, dz, z0, z1 > z0 ? z1 : 1.0e5);
    }
    if (c->axis_pose == NULL) c->axis_pose = c->axis_verdict;   /* the umbilicus line: s = z, r about the umbilicus */
    c->pose.axis = c->axis_pose;
    c->placeo.axis = c->axis_pose;
    c->verdict.axis = c->axis_verdict;
    c->pose.axis_dir[0] = 1.0; c->pose.axis_dir[1] = 0.0; c->pose.axis_dir[2] = 0.0;
    if (c->axis_pose) AsmAxis_mean_dir(c->axis_pose, z0 - 512.0, z1 + 512.0, c->pose.axis_dir);
    const AsmAxis *ap = c->axis_pose, *av = c->axis_verdict;
    fprintf(stderr, "[assemble] axis: pose + placement %s%s; verdict %s%s; mean direction over the box (%.4f, %.4f, %.4f) zyx, %.1f deg off vertical\n",
            ap ? (ap->is_line ? "umbilicus line" : "table ") : "none", ap && !ap->is_line ? c->axis_table : "",
            av ? (av->is_line ? "umbilicus line" : "table ") : "none (world frame)", av && !av->is_line ? c->verdict_axis_table : "",
            c->pose.axis_dir[0], c->pose.axis_dir[1], c->pose.axis_dir[2], acos(c->pose.axis_dir[0]) * 180.0 / 3.14159265358979323846);
    if (ap && !ap->is_line)
        fprintf(stderr, "[assemble] axis LOESS through %zu rows: rows off the curve %.2f vox rms, smallest curvature radius %.0f vox (the frame is rigid only well inside it)\n",
                ap->n_rows, ap->fit_rms, ap->r_curv_min);
}
/* ---- stage 1b: axis (derived from the material) ------------------------------- */

static int sa_stage_axis(AsmRun *run, SaConfig *cfg, const char *out_dir)
{
    int derive = cfg->axis_derive >= 0 ? cfg->axis_derive : (cfg->axis_table[0] == 0);
    if (!derive) {
        fprintf(stderr, "[assemble axis] derivation off (axis.derive %d, pose table %s): the configured frame stands\n", cfg->axis_derive, cfg->axis_table[0] ? cfg->axis_table : "none");
        return 0;
    }
    AsmAxisDeriveOpts o; AsmAxisDerive_default_opts(&o);
    o.slab_vox = cfg->axis_slab_vox;
    const AsmAxis *prior = NULL;
    if (cfg->axis_prior_table[0]) {
        prior = AsmAxis_load(run->arena, cfg->axis_prior_table);
        if (!prior) fprintf(stderr, "[assemble axis] prior table %s: cannot read\n", cfg->axis_prior_table);
    }
    if (!prior) prior = cfg->axis_pose;   /* the configured table, the umbilicus line, or NULL = seed */
    o.prior = prior;
    AsmAxisDeriveReport rep;
    AsmAxis *a = AsmAxis_derive_run(run->arena, run, &o, &rep);
    fprintf(stderr, "[assemble axis] direction from the normals alone (no centre, no circularity): (%.4f, %.4f, %.4f) zyx, %.1f deg off vertical, separation %.4f, %.0f%% of the area within 12 deg of perpendicular\n",
            rep.normal_dir[0], rep.normal_dir[1], rep.normal_dir[2], rep.normal_tilt_deg, rep.normal_cond, 100.0 * rep.normal_frac);
    char path[2048];
    snprintf(path, sizeof path, "%s/stage1_axis.csv", out_dir);
    if (AsmAxisDerive_write_csv(&rep, path) == 0) fprintf(stderr, "[assemble axis] wrote: %s\n", path);
    snprintf(path, sizeof path, "%s/stage1_axis.png", out_dir);
    if (a && AsmAxisDerive_write_png(run->arena, a, prior, &rep, path) == 0) fprintf(stderr, "[assemble axis] wrote: %s\n", path);
    fprintf(stderr, "[assemble axis] derived from %zu face samples: %zu of %zu slabs locked (%zu one-sided), prior %s, %d iterations (last move %.1f vox) | tilt %.1f deg off vertical, rows off the smoothed curve %.2f vox rms, curvature radius >= %.0f vox, cond min %.3f, aperture min %.0f deg, residual p50 max %.1f vox | vs the prior %.1f vox rms | %.1f s\n",
            rep.samples, rep.n_lock, rep.n_rows, rep.one_sided,
            rep.seeded == 0 ? (prior && prior->is_line ? "the umbilicus line" : "the configured table") : rep.seeded == 1 ? "a seed found on z-slabs" : "NONE (world-z fallback)",
            rep.iters, rep.last_move, rep.tilt_deg, rep.fit_rms, rep.r_curv_min < 1e299 ? rep.r_curv_min : 0.0, rep.cond_min < 1e299 ? rep.cond_min : 0.0,
            rep.aperture_min_deg < 1e299 ? rep.aperture_min_deg : 0.0, rep.resid_p50_max, rep.prior_rms, rep.sec);
    if (a) fprintf(stderr, "[assemble axis] the derived curve leans %.1f deg from the normals' own direction\n", rep.normal_disagree_deg);
    snprintf(path, sizeof path, "%s/stage1_axis.json", out_dir);
    FILE *fp = fopen(path, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"axis\",\n  \"derived\": %s,\n  \"accepted\": %s,\n  \"samples\": %zu,\n  \"slabs\": %zu,\n  \"locked\": %zu,\n  \"one_sided\": %zu,\n  \"seeded\": %d,\n  \"seed_slabs\": %zu,\n  \"iters\": %d,\n  \"last_move\": %.3f,\n"
                    "  \"tilt_deg\": %.3f,\n  \"fit_rms\": %.3f,\n  \"r_curv_min\": %.1f,\n  \"cond_min\": %.4f,\n  \"aperture_min_deg\": %.1f,\n  \"resid_p50_max\": %.2f,\n  \"prior_rms\": %.3f,\n"
                    "  \"normal_dir\": [%.6f, %.6f, %.6f],\n  \"normal_tilt_deg\": %.3f,\n  \"normal_cond\": %.6f,\n  \"normal_frac\": %.4f,\n  \"normal_disagree_deg\": %.3f,\n  \"sec\": %.3f\n}\n",
                a ? "true" : "false", a && rep.accepted ? "true" : "false", rep.samples, rep.n_rows, rep.n_lock, rep.one_sided, rep.seeded, rep.seed_slabs, rep.iters, rep.last_move,
                rep.tilt_deg, rep.fit_rms, rep.r_curv_min < 1e299 ? rep.r_curv_min : 0.0, rep.cond_min < 1e299 ? rep.cond_min : 0.0,
                rep.aperture_min_deg < 1e299 ? rep.aperture_min_deg : 0.0, rep.resid_p50_max, rep.prior_rms,
                rep.normal_dir[0], rep.normal_dir[1], rep.normal_dir[2], rep.normal_tilt_deg, rep.normal_cond, rep.normal_frac, rep.normal_disagree_deg, rep.sec);
        fclose(fp);
    }
    if (!a) { fprintf(stderr, "[assemble axis] no slab locked: the configured frame stands\n"); return 0; }
    if (!rep.accepted) {
        fprintf(stderr, "[assemble axis] NOT ACCEPTED (locked %zu of %zu slabs, rows off the curve %.1f vox rms; needs %.0f%% and <= %.0f vox): the configured frame stands\n",
                rep.n_lock, rep.n_rows, rep.fit_rms, 100.0 * ASM_AXIS_DERIVE_ACCEPT_LOCK, ASM_AXIS_DERIVE_ACCEPT_RMS);
        return 0;
    }
    /* the derived axis frames the pose prior, the placement and the verdict */
    double z0 = cfg->have_bbox ? (double)cfg->pile_bbox[0] : 0.0, z1 = cfg->have_bbox ? (double)cfg->pile_bbox[1] + 128.0 : 0.0;
    cfg->axis_derived = a;
    cfg->axis_pose = a; cfg->axis_verdict = a;
    cfg->pose.axis = a; cfg->placeo.axis = a; cfg->verdict.axis = a;
    AsmAxis_mean_dir(a, z0 - 512.0, z1 + 512.0, cfg->pose.axis_dir);
    fprintf(stderr, "[assemble axis] the derived axis frames the pose prior, the placement and the verdict: mean direction over the box (%.4f, %.4f, %.4f) zyx, %.1f deg off vertical\n",
            cfg->pose.axis_dir[0], cfg->pose.axis_dir[1], cfg->pose.axis_dir[2], acos(cfg->pose.axis_dir[0] > 1.0 ? 1.0 : cfg->pose.axis_dir[0]) * 180.0 / 3.14159265358979323846);
    return 0;
}

/* ---- stage 2b: layer cut (joins that close a zero-winding path between wraps) ---- */

static int sa_layer_cut(AsmRun *run, const AsmAxis *axis, const char *out_dir, const char *tag)
{
    AsmLayerCutOpts o;
    AsmLayerCut_defaults(&o);
    o.axis = axis;
    char path[2048];
    snprintf(path, sizeof path, "%s/%s_layer_cut.csv", out_dir, tag);
    FILE *ledger = fopen(path, "wb");
    if (!ledger) return -1;
    AsmLayerCutStats st;
    int rc = AsmLayerCut_run(run, &o, &st, ledger);
    if (fclose(ledger)) rc = -1;
    if (rc) return -1;
    fprintf(stderr, "[assemble layer cut] %zu layer pairs, %zu inside one join component of %zu joins: %zu zero-winding witness paths "
            "within %d joins -> %zu after dropping %zu joins in %d passes (%s) | %.1f s\n",
            st.pairs, st.pairs_joined, st.joins_first, st.witnesses_first, o.max_hops, st.witnesses_last, st.dropped, st.passes,
            axis ? "axis" : "NO AXIS: off", st.sec);
    snprintf(path, sizeof path, "%s/%s_layer_cut.json", out_dir, tag);
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    int okay = fprintf(fp, "{\n  \"stage\": \"layer_cut\",\n  \"axis\": %s,\n  \"layer_pairs\": %zu,\n  \"pairs_in_one_component\": %zu,\n"
                       "  \"joins\": %zu,\n  \"max_hops\": %d,\n  \"min_radius\": %.1f,\n  \"witnesses_first\": %zu,\n  \"witnesses_last\": %zu,\n"
                       "  \"dropped\": %zu,\n  \"passes\": %d,\n  \"sec\": %.3f\n}\n",
                       axis ? "true" : "false", st.pairs, st.pairs_joined, st.joins_first, o.max_hops, o.min_radius,
                       st.witnesses_first, st.witnesses_last, st.dropped, st.passes, st.sec) >= 0;
    if (fclose(fp)) okay = 0;
    return okay ? 0 : -1;
}

static int sa_stage_layer_cut(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    if (!ASM_LAYER_CUT) return 0;
    return sa_layer_cut(run, cfg->axis_pose, out_dir, "stage2b");
}

static int sa_stage_pose(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    AsmPoseStats st;
    ((SaConfig *)cfg)->pose.verbose_prior = 1;
    if (AsmPose_solve(run, &cfg->pose, &st) != 0) return -1;
    ((SaConfig *)cfg)->pose.verbose_prior = 0;
    fprintf(stderr,
        "[assemble pose] charts %zu rels %zu components %zu (largest %zu charts, area %.3e) parity conflicts %zu "
        "switched off %zu -> components after switching %zu | residual sigma p50 %.2f p95 %.2f max %.2f | xy p50 %.2f p95 %.2f vox | %d iters %.1f s\n",
        st.charts_active, st.rels_active, st.components, st.largest_component_charts, st.largest_component_area,
        st.parity_conflicts, st.switched_off, st.components_after_switch, st.resid_p50, st.resid_p95, st.resid_max, st.xy_p50, st.xy_p95,
        st.iters, st.sec);
    fprintf(stderr, "[assemble pose] axial plane: %zu charts, largest component v = %.1f + %.4f s + %.4f u, |v - plane| p50 %.1f p90 %.1f -> %.1f %.1f vox | frame directions: %zu components from their anchor, %zu from the charts' mean, %zu none; the largest component's from %s\n",
            st.vprior_charts, st.vprior_fit[0], st.vprior_fit[1], st.vprior_fit[2], st.vplane_p50_before, st.vplane_p90_before, st.vplane_p50_after, st.vplane_p90_after,
            st.frames_from_anchor, st.frames_from_mean, st.frames_none, st.largest_frame == 1 ? "its anchor" : st.largest_frame == 2 ? "the charts' mean (anchor outside the gradient band)" : "NOWHERE (no chart under the prior)");
    {
        char csv[2048];
        snprintf(csv, sizeof csv, "%s/stage3_poses.csv", out_dir);
        if (AsmReport_poses_csv(run, csv, "pose") == 0) fprintf(stderr, "[assemble pose] wrote: %s\n", csv);
    }
    fprintf(stderr, "[assemble pose] degree histogram (charts 0/1/2/3/4/5+): %zu %zu %zu %zu %zu %zu | area %.2e %.2e %.2e %.2e %.2e %.2e\n",
            st.degree_hist[0], st.degree_hist[1], st.degree_hist[2], st.degree_hist[3], st.degree_hist[4], st.degree_hist[5],
            st.degree_area[0], st.degree_area[1], st.degree_area[2], st.degree_area[3], st.degree_area[4], st.degree_area[5]);
    char png[2048];
    snprintf(png, sizeof png, "%s/stage3_layout.png", out_dir);
    if (AsmReport_layout_png(run->arena, png, run->charts, run->n_charts, NULL, 0, cfg->conflict.cell, 8192, 24) == 0)
        fprintf(stderr, "[assemble pose] wrote: %s\n", png);
    char rep[2048];
    snprintf(rep, sizeof rep, "%s/stage3_pose.json", out_dir);
    FILE *fp = fopen(rep, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"pose\",\n  \"charts_active\": %zu,\n  \"rels_active\": %zu,\n  \"components\": %zu,\n"
                    "  \"largest_component_charts\": %zu,\n  \"largest_component_area\": %.6e,\n  \"parity_conflicts\": %zu,\n"
                    "  \"switched_off\": %zu,\n  \"residual_sigma\": {\"p50\": %.4f, \"p95\": %.4f, \"max\": %.4f},\n"
                    "  \"xy_residual_vox\": {\"p50\": %.4f, \"p95\": %.4f},\n  \"iters\": %d,\n  \"sec\": %.3f\n}\n",
                st.charts_active, st.rels_active, st.components, st.largest_component_charts, st.largest_component_area,
                st.parity_conflicts, st.switched_off, st.resid_p50, st.resid_p95, st.resid_max, st.xy_p50, st.xy_p95,
                st.iters, st.sec);
        fclose(fp);
    }
    return 0;
}

/* ---- stage 4: clean2 (contradiction audit + cleaning) ------------------------ */

static int sa_stage_clean2(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    AsmConflictReport rep0;
    char pairs_csv[2048];
    snprintf(pairs_csv, sizeof pairs_csv, "%s/stage4_contra_pairs.csv", out_dir);
    AsmConflictOpts co = cfg->conflict;
    co.diag_pairs_csv = pairs_csv;   /* the FIRST audit's contradicting pairs, before any drop */
    if (AsmConflict_audit(run, &co, NULL, &rep0) != 0) return -1;
    fprintf(stderr, "[assemble clean2] audit before: cells %zu multi %zu contra %zu | covered %.4e contra %.4e score %.4e | pairs %zu (ledger %s)\n",
            rep0.cells, rep0.cells_multi, rep0.cells_contra, rep0.covered_area, rep0.contra_area, rep0.score, rep0.pairs, pairs_csv);
    char png[2048];
    snprintf(png, sizeof png, "%s/stage4_before.png", out_dir);
    if (AsmReport_layout_png(run->arena, png, run->charts, run->n_charts, rep0.contra_cells, rep0.n_contra_cells, cfg->conflict.cell, 8192, 24) == 0)
        fprintf(stderr, "[assemble clean2] wrote: %s\n", png);
    AsmConflictReport rep;
    AsmConflictStats st;
    if (AsmConflict_clean(run, &cfg->conflict, &cfg->pose, &rep, &st) != 0) return -1;
    fprintf(stderr, "[assemble clean2] rounds %d dropped %zu reverted %zu extras %zu splits refused %zu (%.3e vox^2 kept whole) | score %.4e -> %.4e covered %.4e -> %.4e contra %.4e -> %.4e pairs %zu -> %zu | %.1f s\n",
            st.rounds, st.rels_dropped, st.rels_reverted, st.charts_extras, st.splits_refused, st.split_area_refused, st.score_first, st.score_last,
            st.covered_first, st.covered_last, st.contra_first, st.contra_last, st.pairs_first, st.pairs_last, st.sec);
    snprintf(png, sizeof png, "%s/stage4_after.png", out_dir);
    if (AsmReport_layout_png(run->arena, png, run->charts, run->n_charts, rep.contra_cells, rep.n_contra_cells, cfg->conflict.cell, 8192, 24) == 0)
        fprintf(stderr, "[assemble clean2] wrote: %s\n", png);
    char path[2048];
    snprintf(path, sizeof path, "%s/stage4_clean2.json", out_dir);
    FILE *fp = fopen(path, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"clean2\",\n  \"rounds\": %d,\n  \"rels_dropped\": %zu,\n  \"rels_reverted\": %zu,\n  \"charts_extras\": %zu,\n"
                    "  \"splits_refused\": %zu,\n  \"split_area_refused\": %.6e,\n"
                    "  \"score\": {\"first\": %.6e, \"last\": %.6e},\n  \"covered\": {\"first\": %.6e, \"last\": %.6e},\n"
                    "  \"contra\": {\"first\": %.6e, \"last\": %.6e},\n  \"pairs\": {\"first\": %zu, \"last\": %zu},\n  \"sec\": %.3f\n}\n",
                st.rounds, st.rels_dropped, st.rels_reverted, st.charts_extras, st.splits_refused, st.split_area_refused, st.score_first, st.score_last,
                st.covered_first, st.covered_last, st.contra_first, st.contra_last, st.pairs_first, st.pairs_last, st.sec);
        fclose(fp);
    }
    return 0;
}

/* ---- stage 5: place ----------------------------------------------------------- */

static size_t js_joined_for_scores = 0;   /* the place stage's layout-join count, reported with the scores */
static AsmPlaceStats sa_place_for_scores;   /* the place stage's wobble numbers, reported with the scores */
static AsmRefineStats sa_refine_for_scores;  /* the seam re-solve's counts, reported with the scores */
static int sa_place_ran = 0;                 /* the place stage ran in this process (its axis verdict is known) */

static int sa_stage_place(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    AsmPlaceStats st;
    AsmPlaceOpts po = cfg->placeo;
    po.diag_dir = out_dir;
    po.stress_band = cfg->clean.stress_band;
    if (AsmPlace_run(run, &po, &cfg->conflict, &st) != 0) return -1;
    fprintf(stderr, "[assemble place] components %zu placed %zu unplaced %zu (mirrored %zu, reverted by audit %zu, hop retries %zu of which placed %zu) sweeps %zu anchors %zu (%zu from the radius) | never attempted %zu components %.3e (%.1f%% of the pile) | v from the axis: %zu pinned (max shift %.1f vox), %zu refused | area in %.3e primary %.3e placed %.3e (%.1f%%) | rms p50 %.2f max %.2f | %.1f s\n",
            st.components_in, st.placed, st.unplaced, st.oriented_mirrored, st.reverted_by_audit, st.hop_retries, st.hop_retry_placed, st.sweeps, st.anchors_used, st.anchors_radius,
            st.never_attempted, st.never_attempted_area, st.area_in > 0 ? 100.0 * st.never_attempted_area / st.area_in : 0.0,
            st.v_pinned, st.v_shift_max, st.v_refused, st.area_in, st.area_primary, st.area_placed,
            st.area_in > 0 ? 100.0 * st.area_placed / st.area_in : 0.0, st.rms_p50, st.rms_max, st.sec);
    fprintf(stderr, "[assemble place] confetti: single-chart components %zu (%.1f%% of the area): placed %zu (%.2e) reverted %zu (%.2e) refused %zu (%.2e) deferred %zu (%.2e) never %zu (%.2e) few %zu | placed by seam %zu by hop %zu (v kept from the seam %zu) | placed with rms > 30 vox %zu (%.2e) | reverted area %.2e (%.1f%%) | deferred at the end %zu components (%.2e, %.1f%%; %zu deferrals; last call placed %zu) | seam evidence on the unplaced: accepted-dropped %zu weak %zu short %zu on %zu components (%.2e), none on %zu (%.2e)\n",
            st.singles, st.area_in > 0 ? 100.0 * st.singles_area / st.area_in : 0.0,
            st.singles_by_status[1], st.singles_area_by_status[1], st.singles_by_status[3], st.singles_area_by_status[3], st.singles_by_status[2], st.singles_area_by_status[2],
            st.singles_by_status[6], st.singles_area_by_status[6], st.singles_by_status[0], st.singles_area_by_status[0], st.singles_by_status[5],
            st.seam_placed, st.hop_placed, st.v_seam_kept, st.placed_rms30, st.placed_rms30_area,
            st.reverted_area, st.area_in > 0 ? 100.0 * st.reverted_area / st.area_in : 0.0,
            st.deferred_comps, st.deferred_area, st.area_in > 0 ? 100.0 * st.deferred_area / st.area_in : 0.0, st.deferred, st.last_call_placed,
            st.evid_dropped, st.evid_weak, st.evid_short, st.evid_comps, st.evid_area, st.noevid_comps, st.noevid_area);
    AsmConflictReport rep;
    AsmConflict_audit(run, &cfg->conflict, NULL, &rep);
    fprintf(stderr, "[assemble place] audit after placement: covered %.4e contra %.4e (%.3f%%) pairs %zu\n",
            rep.covered_area, rep.contra_area, rep.covered_area > 0 ? 100.0 * rep.contra_area / rep.covered_area : 0.0, rep.pairs);
    AsmJoinStats js;
    {
        char csv[2048];
        snprintf(csv, sizeof csv, "%s/stage5_layout_joins.csv", out_dir);
        if (AsmConflict_confirm_joins(run, &cfg->conflict, csv, &js) != 0) return -1;
        js_joined_for_scores = js.pairs_joined;
        sa_place_for_scores = st;
        sa_place_ran = 1;
        fprintf(stderr, "[assemble place] wobble against v(s): primary |v - v(s,u)| p50 %.1f p90 %.1f vox | placed area beyond 30 vox %.1f%%, beyond 100 vox %.1f%% | spine steps per 2,000 vox of u: p50 %.1f p90 %.1f (%zu bins) | v span %.0f vox | outline: the box's top / bottom crop plane traces wave p2p %.0f / %.0f vox over %zu / %zu bins (the planes tilted %.1f deg against the axis at r %.0f predict %.0f)\n",
                st.wobble_primary_p50, st.wobble_primary_p90, 100.0 * st.placed_area_beyond30, 100.0 * st.placed_area_beyond100,
                st.spine_step_p50, st.spine_step_p90, st.spine_bins, st.v_span, st.outline_top_p2p, st.outline_bot_p2p, st.outline_bins_top, st.outline_bins_bot, st.axis_tilt_deg, st.outline_r_mean, st.outline_pred_p2p);
        fprintf(stderr, "[assemble place] layout-confirmed joins: pairs tested %zu joined %zu (refused: contradicting %zu, short %zu, contact-vetoed %zu) | agreeing cells %zu disagreeing %zu | lineages %zu -> %zu | aligned: %zu charts moved, max shift %.2f vox, max rotation %.4f rad | %.1f s\n",
                js.pairs_tested, js.pairs_joined, js.pairs_refused_disagree, js.pairs_refused_short, js.pairs_refused_veto, js.cells_agree, js.cells_disagree,
                js.lineages_before, js.lineages_after, js.charts_shifted, js.shift_max, js.rot_max, js.sec);
        if (ASM_REFINE_SEAMS) {
            AsmRefineStats rs;
            char rcsv[2048];
            snprintf(rcsv, sizeof rcsv, "%s/stage5_refine.csv", out_dir);
            AsmConflict_refine_seams(run, &cfg->conflict, &cfg->pose, rcsv, &rs);
            sa_refine_for_scores = rs;
            fprintf(stderr, "[assemble refine] seams: candidates %zu joins + %zu dropped %zu switched %zu weak %zu short -> readmitted %zu (dropped %zu switched %zu weak %zu short %zu; refused: gap %zu rotation %zu parity %zu no-corr %zu; %zu internal to the fixed lineage) | placed charts %zu, fixed %zu, free %zu, pinned %zu | rounds %d iters %d switched-off %zu re-refused %zu un-readmitted %zu | moved > 2 vox %zu charts (%.2e; shift p50 %.1f p90 %.1f max %.1f vox, rot max %.4f rad) | contra %.3e -> %.3e (%+.1f%%) %s | lineages %zu -> %zu | %.1f s\n",
                    rs.cand_join, rs.cand_by_kind[1], rs.cand_by_kind[2], rs.cand_by_kind[3], rs.cand_by_kind[4], rs.readmitted, rs.readmit_by_kind[1], rs.readmit_by_kind[2], rs.readmit_by_kind[3], rs.readmit_by_kind[4],
                    rs.refused_gap, rs.refused_rot, rs.refused_parity, rs.refused_nocorr, rs.internal_fixed, rs.charts_placed, rs.charts_fixed, rs.charts_free, rs.charts_pinned,
                    rs.rounds, rs.iters, rs.switched, rs.readmit_switched, rs.unreadmitted, rs.moved, rs.moved_area, rs.shift_p50, rs.shift_p90, rs.shift_max, rs.rot_max,
                    rs.contra_before, rs.contra_after, rs.contra_before > 0 ? 100.0 * (rs.contra_after - rs.contra_before) / rs.contra_before : 0.0,
                    rs.reverted ? "REVERTED" : (rs.kept ? "KEPT" : "nothing to readmit"), rs.lineages_before, rs.lineages_after, rs.sec);
            if (ASM_JOIN_AFTER_REFINE && rs.kept) {
                AsmJoinStats js2;
                snprintf(rcsv, sizeof rcsv, "%s/stage5_layout_joins2.csv", out_dir);
                AsmConflict_confirm_joins(run, &cfg->conflict, rcsv, &js2);
                js_joined_for_scores += js2.pairs_joined;
                fprintf(stderr, "[assemble place] layout-confirmed joins after the re-solve: pairs tested %zu joined %zu (refused: contradicting %zu, short %zu) | lineages %zu -> %zu | aligned: %zu charts moved, max shift %.2f vox | %.1f s\n",
                        js2.pairs_tested, js2.pairs_joined, js2.pairs_refused_disagree, js2.pairs_refused_short, js2.lineages_before, js2.lineages_after, js2.charts_shifted, js2.shift_max, js2.sec);
            }
        }
    }
    /* Report the geometry actually emitted, after every join decision. */
    AsmConflict_audit(run, &cfg->conflict, NULL, &rep);
    char png[2048];
    snprintf(png, sizeof png, "%s/stage5_placed.png", out_dir);
    if (AsmReport_layout_png(run->arena, png, run->charts, run->n_charts, rep.contra_cells, rep.n_contra_cells, cfg->conflict.cell, 8192, 8) == 0)
        fprintf(stderr, "[assemble place] wrote: %s\n", png);
    char path[2048];
    snprintf(path, sizeof path, "%s/stage5_place.json", out_dir);
    FILE *fp = fopen(path, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"place\",\n  \"components\": %zu,\n  \"placed\": %zu,\n  \"unplaced\": %zu,\n  \"mirrored\": %zu,\n"
                    "  \"anchors\": %zu,\n  \"area_in\": %.6e,\n  \"area_primary\": %.6e,\n  \"area_placed\": %.6e,\n"
                    "  \"rms_p50\": %.4f,\n  \"rms_max\": %.4f,\n  \"audit\": {\"covered\": %.6e, \"contra\": %.6e, \"pairs\": %zu},\n"
                    "  \"layout_joins\": {\"tested\": %zu, \"joined\": %zu, \"refused_contradicting\": %zu, \"refused_short\": %zu, \"lineages_before\": %zu, \"lineages_after\": %zu},\n  \"sec\": %.3f\n}\n",
                st.components_in, st.placed, st.unplaced, st.oriented_mirrored, st.anchors_used, st.area_in, st.area_primary, st.area_placed,
                st.rms_p50, st.rms_max, rep.covered_area, rep.contra_area, rep.pairs,
                js.pairs_tested, js.pairs_joined, js.pairs_refused_disagree, js.pairs_refused_short, js.lineages_before, js.lineages_after, st.sec);
        fclose(fp);
    }
    return 0;
}

/* ---- stage 5b: discover -------------------------------------------------------- */

static int sa_stage_discover(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    AsmDiscoverStats st;
    if (AsmDiscover_run(run, &cfg->discover, out_dir, &st) != 0) return -1;
    fprintf(stderr, "[assemble discover] control: %zu partner-implied poses of placed charts land p50 %.2f p90 %.2f max %.2f vox from the actual pose, %zu beyond %.0f vox | leftovers examined %zu, with a placed partner %zu | %s %zu charts (%.3e vox^2, %zu onto same-surface material) in %d rounds | refused: sides disagree %zu, too few sides %zu, rim %zu, another wrap %zu | %.1f s\n",
            st.control_poses, st.control_p50, st.control_p90, st.control_max, st.control_beyond, (double)ASM_DISCOVER_AGREE_VOX,
            st.examined, st.with_partner, cfg->discover.enabled ? "placed" : "would place (report only)", st.admitted, st.admitted_area, st.admitted_overlapping, st.rounds,
            st.refused_disagree, st.refused_sides, st.refused_rim, st.refused_wrap, st.sec);
    fprintf(stderr, "[assemble discover] wrote: %s/stage5_discovery.csv, %s/stage5_discovered.png\n", out_dir, out_dir);
    return 0;
}

/* ---- stage 6: sheet ----------------------------------------------------------- */

static int sa_stage_sheet(AsmRun *run, const SaConfig *cfg, const char *out_dir)
{
    AsmEmitStats st;int adapter_rc=0;
    if (AsmEmit_run(run, &cfg->emit, out_dir, &st) != 0) return -1;
    fprintf(stderr, "[assemble sheet] sheet charts %zu verts %zu faces %zu area %.3e (u span %.0f v span %.0f, gauge %.3f rad) | extras charts %zu verts %zu area %.3e | bake rc %d | %.1f s\n",
            st.sheet_charts, st.sheet_verts, st.sheet_faces, st.sheet_area, st.u_span, st.v_span, st.gauge_theta,
            st.extras_charts, st.extras_verts, st.extras_area, st.bake_rc, st.sec);
    fprintf(stderr, "[assemble sheet] wrote: %s/sheet.vmesh (+obj, sidecars), %s/sheet_extras.vmesh%s\n", out_dir, out_dir,
            st.bake_rc == 0 ? ", bake/" : "");
    {
        AsmVerdictStats vs;
        /* the arc-length gate reaches the verdict frame: a table that is not this sheet's axis
         * (placement disarmed it) must not frame the measurement either; the configured umbilicus
         * line takes over, as before round 3 */
        AsmVerdictOpts vo = cfg->verdict;
        if (sa_place_ran && !sa_place_for_scores.axis_armed && vo.axis != NULL && !vo.axis->is_line && cfg->verdict.umb_y > 0.0 && cfg->verdict.umb_x > 0.0) {
            double z0 = cfg->have_bbox ? (double)cfg->pile_bbox[0] : 0.0, z1 = cfg->have_bbox ? (double)cfg->pile_bbox[1] + 128.0 : 0.0;
            double pt[3] = { 0.0, cfg->verdict.umb_y, cfg->verdict.umb_x }, dz[3] = { 1.0, 0.0, 0.0 };
            vo.axis = AsmAxis_line(run->arena, pt, dz, z0, z1 > z0 ? z1 : 1.0e5);
            fprintf(stderr, "[assemble verdict] axis table DISARMED for the measurement frame (placement's b %.4f outside the arc-length band): the verdict measures about the umbilicus line (%.0f, %.0f)\n",
                    sa_place_for_scores.axis_b, cfg->verdict.umb_y, cfg->verdict.umb_x);
        }
        if (AsmVerdict_run(run, &vo, out_dir, &vs) == 0) {
            fprintf(stderr, "[assemble verdict] verdict lattice (v-rows, %s frame): verts %zu faces %zu rows %zu cols %zu | contested cells %zu -> stacked verts %zu | lineages %zu | cross-chart faces by relation %zu by geometry %zu, quads refused %zu (+%zu between lineages) | long-join edges %zu on %zu joins (%zu within 8 vox = the trim gap's tail) | bridged cells %zu | ribbon_verdict rc %d (report %s/verdict.json) | %.1f s\n",
                    vs.axis_frame ? "local-axis" : "world", vs.lattice_verts, vs.lattice_faces, vs.lattice_rows, vs.lattice_cols,
                    vs.cells_contested, vs.stacked_verts, vs.lineages, vs.faces_by_relation, vs.faces_by_geometry, vs.quads_refused, vs.quads_refused_lineage,
                    vs.long_join_edges, vs.long_join_pairs, vs.long_join_edges_tail, vs.cells_bridged_uv, vs.verdict_rc, out_dir, vs.sec);
            fprintf(stderr, "[assemble verdict] assembled ribbon for --solidify (z-rows, world): verts %zu faces %zu rows %zu cols %zu contested cells %zu seam-row cells %zu | source verts %zu faces %zu\n",
                    vs.solid_verts, vs.solid_faces, vs.solid_rows, vs.solid_cols, vs.solid_contested, vs.cells_bridged, vs.source_verts, vs.source_faces);
            char vlog[2048];
            snprintf(vlog, sizeof vlog, "%s/verdict.log", out_dir);
            FILE *lf = fopen(vlog, "rb");
            if (lf) {
                char line[512];
                while (fgets(line, sizeof line, lf)) if (strstr(line, "VERDICT") || strstr(line, "[gate") || strstr(line, "PASS") || strstr(line, "FAIL")) fputs(line, stderr);
                fclose(lf);
            }
            /* SCORES: the eight gates plus the discriminators the gates do not carry (2026-09-09:
             * the coherence gate measures the largest column run's band and rises when the sheet
             * fragments; stacked (u,z) cells expose a wrap placed one band off in v). */
            {
                char vpath[2048], spath[2048];
                snprintf(vpath, sizeof vpath, "%s/verdict.log", out_dir);
                double coverage = -1.0, coherence = -1.0, largest_overall = -1.0, seams = -1.0;
                size_t long_edges = 0, cross_wrap = 0, stacked = 0, run_cells = 0, components = 0, bands = 0;
                int gates_pass = 0, gates_n = 0, cross_wrap_measured = 1, depth_measured = 1;
                FILE *lf2 = fopen(vpath, "rb");
                if (lf2) {
                    char line[1024];
                    while (fgets(line, sizeof line, lf2)) {
                        const char *p = line;
                        while (*p == ' ') p++;
                        /* PASS / FAIL / ABSENT: an ABSENT gate counts, fails, and carries no number */
                        int gp = 0, gm = 1; const char *gv = NULL;
                        if (strncmp(p, "TANGLE/long-edges", 17) == 0) { gv = AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp; if (gv && gm) sscanf(gv, " %zu", &long_edges); }
                        else if (strncmp(p, "TANGLE/cross-wrap", 17) == 0) { gv = AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp; cross_wrap_measured = gm; if (gv && gm) sscanf(gv, " %zu", &cross_wrap); }
                        else if (strncmp(p, "WINDING/single-cover", 20) == 0) { gv = AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp; if (gv && gm) sscanf(gv, " %zu", &stacked); }
                        else if (strncmp(p, "COVERAGE/source", 15) == 0) { gv = AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp; if (gv && gm) sscanf(gv, " %lf%%", &coverage); }
                        else if (strncmp(p, "FRAGMENTATION/coherence", 23) == 0) {
                            AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp;
                            const char *h = strstr(p, "holds "); if (h) sscanf(h, "holds %zu cells, %lf%%", &run_cells, &coherence);
                            const char *c = strstr(p, "bands, "); if (c) sscanf(c, "bands, %zu components, largest %lf%%", &components, &largest_overall);
                            const char *b = strstr(p, " runs, "); if (b) sscanf(b, " runs, %zu bands", &bands);
                        }
                        else if (strncmp(p, "TEXTURE/depth-seams", 19) == 0) { AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp; depth_measured = gm; const char *s2 = strstr(p, "; "); if (s2 && gm) sscanf(s2, "; %lf%%", &seams); }
                        else if (strncmp(p, "WINDING/material-purity", 23) == 0 || strncmp(p, "DISTORTION/anisotropy", 21) == 0) { AsmVerdict_gate_word(p, &gp, &gm); gates_n++; gates_pass += gp; }
                    }
                    fclose(lf2);
                }
                fprintf(stderr, "[assemble scores] gates %d/%d | coverage %.1f%% | coherence %.1f%% of the largest run (%zu cells, %zu bands) | largest component %.1f%% overall of %zu | stacked (u,v) %zu (u,z) %zu | long edges %zu cross-wrap %zu | depth seams %.3f%% | lineages %zu layout joins %zu | seams readmitted %zu (%zu charts snapped)\n",
                        gates_pass, gates_n, coverage, coherence, run_cells, bands, largest_overall, components, stacked, vs.solid_contested, long_edges, cross_wrap, seams, vs.lineages, js_joined_for_scores,
                        sa_refine_for_scores.kept ? sa_refine_for_scores.readmitted - sa_refine_for_scores.readmit_switched : 0, sa_refine_for_scores.kept ? sa_refine_for_scores.moved : 0);
                /* The reason is the VERDICT's to give: an absent cross-wrap gate does mean no
                 * frame, but depth-seams is also skipped when the lattice exceeds its cell cap,
                 * and the 21x5x5 -- which has both a table and an umbilicus -- was being told
                 * it had neither (2026-09-17).  The gate line above carries "NOT MEASURED: ..." */
                if (!cross_wrap_measured || !depth_measured)
                    fprintf(stderr, "[assemble scores] %s%s%s ABSENT: counted as failed, and the number above is not a measurement; the gate line says why\n",
                            cross_wrap_measured ? "" : "cross-wrap", !cross_wrap_measured && !depth_measured ? " and " : "", depth_measured ? "" : "depth-seams");
                fprintf(stderr, "[assemble scores] wobble: primary |v - v(s)| p50 %.1f p90 %.1f vox, placed area beyond 30 / 100 vox %.1f%% / %.1f%%, spine step p50 %.1f p90 %.1f, v span %.0f, outline top/bottom wave %.0f / %.0f vox (crop tilt %.1f deg predicts %.0f) | gaps in u: within a wrap %zu (%.0f vox, p50 %.0f p90 %.0f), between wraps %zu (%.0f vox, p50 %.0f p90 %.0f)%s\n",
                        sa_place_for_scores.wobble_primary_p50, sa_place_for_scores.wobble_primary_p90, 100.0 * sa_place_for_scores.placed_area_beyond30, 100.0 * sa_place_for_scores.placed_area_beyond100,
                        sa_place_for_scores.spine_step_p50, sa_place_for_scores.spine_step_p90, sa_place_for_scores.v_span,
                        sa_place_for_scores.outline_top_p2p, sa_place_for_scores.outline_bot_p2p, sa_place_for_scores.axis_tilt_deg, sa_place_for_scores.outline_pred_p2p,
                        vs.gap_within_n, vs.gap_within_vox, vs.gap_within_p50, vs.gap_within_p90, vs.gap_between_n, vs.gap_between_vox, vs.gap_between_p50, vs.gap_between_p90,
                        vs.gap_measured ? "" : " (gaps UNMEASURED: no axis)");
                snprintf(spath, sizeof spath, "%s/scores.json", out_dir);
                FILE *sf = fopen(spath, "wb");
                if (sf) {
                    const AsmContinuityStats *cs = &vs.continuity;
                    fprintf(sf, "{\n  \"confetti_free\": %s,\n  \"continuity_audit_complete\": %s,\n  \"placement_budget_exhausted\": %s,\n"
                                "  \"unsupported_islands\": %zu,\n  \"unsupported_area\": %.6e,\n  \"unexplained_continuation_breaks\": %zu,\n"
                                "  \"metadata_only_joins\": %zu,\n  \"wrong_placements\": %zu,\n  \"unresolved_hypotheses\": %zu,\n  \"verified_joins\": %zu,\n",
                            cs->confetti_free && !sa_place_for_scores.budget_exhausted ? "true" : "false",
                            cs->audit_complete ? "true" : "false", sa_place_for_scores.budget_exhausted ? "true" : "false",
                            cs->unsupported_islands, cs->unsupported_area, cs->unexplained_breaks, cs->metadata_only_joins,
                            cs->wrong_placements, cs->unresolved_hypotheses, cs->verified_joins);
                    fprintf(sf, "  \"gates_pass\": %d,\n  \"gates_n\": %d,\n  \"coverage_pct\": %.2f,\n  \"coherence_pct\": %.2f,\n  \"largest_run_cells\": %zu,\n  \"bands\": %zu,\n"
                                "  \"largest_component_pct_overall\": %.2f,\n  \"components\": %zu,\n  \"stacked_uv_cells\": %zu,\n  \"stacked_uz_cells\": %zu,\n"
                                "  \"long_edges\": %zu,\n  \"long_join_edges\": %zu,\n  \"long_join_edges_tail\": %zu,\n  \"cross_wrap_edges\": %zu,\n  \"depth_seams_pct\": %.4f,\n"
                                "  \"cross_wrap_measured\": %d,\n  \"depth_seams_measured\": %d,\n"
                                "  \"lineages\": %zu,\n  \"layout_joins\": %zu,\n  \"readmitted\": %zu,\n  \"snapped\": %zu,\n  \"axis_frame\": %d,\n"
                                "  \"wobble\": {\"primary_p50\": %.2f, \"primary_p90\": %.2f, \"placed_area_beyond30\": %.4f, \"placed_area_beyond100\": %.4f, \"spine_step_p50\": %.2f, \"spine_step_p90\": %.2f, \"spine_bins\": %zu, \"v_span\": %.1f},\n"
                                "  \"outline\": {\"top_p2p\": %.1f, \"bottom_p2p\": %.1f, \"axis_tilt_deg\": %.2f, \"r_mean\": %.1f, \"predicted_p2p\": %.1f},\n"
                                "  \"gaps\": {\"measured\": %d, \"within_n\": %zu, \"within_vox\": %.0f, \"within_p50\": %.1f, \"within_p90\": %.1f, \"between_n\": %zu, \"between_vox\": %.0f, \"between_p50\": %.1f, \"between_p90\": %.1f}\n}\n",
                            gates_pass, gates_n, coverage, coherence, run_cells, bands, largest_overall, components, stacked, vs.solid_contested,
                            long_edges, vs.long_join_edges, vs.long_join_edges_tail, cross_wrap, seams, cross_wrap_measured, depth_measured, vs.lineages, js_joined_for_scores,
                            sa_refine_for_scores.kept ? sa_refine_for_scores.readmitted - sa_refine_for_scores.readmit_switched : (size_t)0, sa_refine_for_scores.kept ? sa_refine_for_scores.moved : (size_t)0, vs.axis_frame,
                            sa_place_for_scores.wobble_primary_p50, sa_place_for_scores.wobble_primary_p90, sa_place_for_scores.placed_area_beyond30, sa_place_for_scores.placed_area_beyond100,
                            sa_place_for_scores.spine_step_p50, sa_place_for_scores.spine_step_p90, sa_place_for_scores.spine_bins, sa_place_for_scores.v_span,
                            sa_place_for_scores.outline_top_p2p, sa_place_for_scores.outline_bot_p2p, sa_place_for_scores.axis_tilt_deg, sa_place_for_scores.outline_r_mean, sa_place_for_scores.outline_pred_p2p,
                            vs.gap_measured, vs.gap_within_n, vs.gap_within_vox, vs.gap_within_p50, vs.gap_within_p90, vs.gap_between_n, vs.gap_between_vox, vs.gap_between_p50, vs.gap_between_p90);
                    fclose(sf);
                }
            }
        } else {fprintf(stderr, "[assemble verdict] adapter failed\n");adapter_rc=-1;}
    }
    char path[2048];
    snprintf(path, sizeof path, "%s/stage6_sheet.json", out_dir);
    FILE *fp = fopen(path, "wb");
    if (fp) {
        fprintf(fp, "{\n  \"stage\": \"sheet\",\n  \"sheet\": {\"charts\": %zu, \"verts\": %zu, \"faces\": %zu, \"area\": %.6e, \"u_span\": %.2f, \"v_span\": %.2f},\n"
                    "  \"extras\": {\"charts\": %zu, \"verts\": %zu, \"faces\": %zu, \"area\": %.6e},\n  \"gauge_theta\": %.6f,\n  \"bake_rc\": %d,\n  \"sec\": %.3f\n}\n",
                st.sheet_charts, st.sheet_verts, st.sheet_faces, st.sheet_area, st.u_span, st.v_span,
                st.extras_charts, st.extras_verts, st.extras_faces, st.extras_area, st.gauge_theta, st.bake_rc, st.sec);
        fclose(fp);
    }
    return adapter_rc;
}

/* ---- selftest -------------------------------------------------------------- */

static int sa_selftest(void)
{
    int fails = 0;
    fprintf(stderr, "sheet_assemble --selftest\n");
    fails += AsmFlatten_selftest();
    fails += AsmClean_selftest();
    fails += AsmStore_selftest();
    fails += AsmAxis_selftest();
    fails += AsmAxisDerive_selftest();
    fails += AsmReadingOrder_selftest();
    fails += AsmVerdict_selftest();
    fails += AsmRelate_selftest();
    fails += AsmLayerCut_selftest();
    fails += AsmPose_selftest();
    fails += AsmConflict_selftest();
    fails += AsmPlace_selftest();
    fails += AsmContinuity_selftest();
    fails += AsmEmit_selftest();
    fails += MeshBin_selftest();
    fails += AsmContacts_selftest();
    fails += AsmField_selftest();
    fails += AsmAudit_selftest();
    fails += AsmRepair_selftest();
    fails += AsmDiscover_selftest();
    fprintf(stderr, "sheet_assemble selftest: %d failure(s)\n", fails);
    return fails;
}

/* ---- main ------------------------------------------------------------------- */

int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--help")) { sa_usage(); return 0; }
    if (argc == 2 && !strcmp(argv[1], "--selftest")) return sa_selftest();
    if ((argc == 4 || argc == 5) && !strcmp(argv[1], "--review"))
        return (argc == 5 ? AsmRepair_reading_review_case(argv[2], argv[3], argv[4]) :
                           AsmRepair_sheet_review_case(argv[2], argv[3])) != 0;
    if (argc == 5 && !strcmp(argv[1], "--winding-order"))
        return AsmRepair_winding_ribbon_case(argv[2], argv[3], argv[4]) != 0;
    const char *pile_dir = NULL, *out_dir = NULL, *config_path = NULL, *stop_after = NULL, *resume_placement = NULL;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--resume-placement") == 0 && i + 1 < argc) { resume_placement = argv[++i]; continue; }
        if (strcmp(argv[i], "--stop-after") == 0 && i + 1 < argc) { stop_after = argv[++i]; continue; }
        if (argv[i][0] == '-') { sa_usage(); return 2; }
        if (!pile_dir) pile_dir = argv[i];
        else if (!out_dir) out_dir = argv[i];
        else if (!config_path) config_path = argv[i];
        else { sa_usage(); return 2; }
    }
    if (!pile_dir || !out_dir) { sa_usage(); return 2; }
    if (stop_after && strcmp(stop_after,"clean") && strcmp(stop_after,"axis") && strcmp(stop_after,"relate") && strcmp(stop_after,"pose") &&
        strcmp(stop_after,"clean2") && strcmp(stop_after,"place") &&
        strcmp(stop_after,"discover") && strcmp(stop_after,"repair") && strcmp(stop_after,"audit") && strcmp(stop_after,"sheet")) { sa_usage(); return 2; }

    double t_start = ves_clock_sec();
    AsmRun run;
    memset(&run, 0, sizeof run);
    run.arena = Arena_new();
    SaConfig cfg;
    sa_config_defaults(&cfg);
    if (!config_path) config_path = "configs/default.json";
    if (sa_config_load(run.arena, config_path, &cfg) != 0) return 2;
    ves_mkdir(out_dir);

    run.pile = ARENA_ALLOC(run.arena, (size_t)SA_MAX_CUBES * sizeof(MeshPileEntry));
    if (MeshPile_scan(pile_dir, run.pile, SA_MAX_CUBES, &run.n_cubes, cfg.have_bbox ? sa_bbox_filter : NULL, &cfg) != 0 || run.n_cubes == 0) {
        fprintf(stderr, "sheet_assemble: no per-cube meshes under %s\n", pile_dir);
        return 3;
    }
    fprintf(stderr, "[assemble] pile %s: %zu cubes; out %s; threads %d; resume %d\n",
            pile_dir, run.n_cubes, out_dir, cfg.threads, cfg.resume);
    sa_axis_setup(&run, &cfg);

    double t0 = ves_clock_sec();
    if(resume_placement){
        if(AsmStore_read_run(resume_placement,&run))return 8;
        fprintf(stderr,"[assemble] resumed complete native placement: %s\n",resume_placement);
        if (sa_stage_axis(&run, &cfg, out_dir) != 0) return 8;   /* the verdict measures in the derived frame */
        goto discover;
    }
    if (sa_stage_clean(&run, &cfg, out_dir) != 0) return 4;
    fprintf(stderr, "[time] clean %.1f s\n", ves_clock_sec() - t0);
    if (stop_after && strcmp(stop_after, "clean") == 0) goto done;

    t0 = ves_clock_sec();
    if (sa_stage_axis(&run, &cfg, out_dir) != 0) return 4;
    fprintf(stderr, "[time] axis %.1f s\n", ves_clock_sec() - t0);
    if (stop_after && strcmp(stop_after, "axis") == 0) goto done;

    t0 = ves_clock_sec();
    if (sa_stage_relate(&run, &cfg, out_dir) != 0) return 5;
    fprintf(stderr, "[time] relate %.1f s\n", ves_clock_sec() - t0);
    if (stop_after && strcmp(stop_after, "relate") == 0) goto done;

    t0 = ves_clock_sec();
    if (sa_stage_layer_cut(&run, &cfg, out_dir) != 0) return 5;
    fprintf(stderr, "[time] layer cut %.1f s\n", ves_clock_sec() - t0);

    t0 = ves_clock_sec();
    if (sa_stage_pose(&run, &cfg, out_dir) != 0) return 6;
    fprintf(stderr, "[time] pose %.1f s\n", ves_clock_sec() - t0);
    if (stop_after && strcmp(stop_after, "pose") == 0) goto done;

    t0 = ves_clock_sec();
    if (sa_stage_clean2(&run, &cfg, out_dir) != 0) return 7;
    fprintf(stderr, "[time] clean2 %.1f s\n", ves_clock_sec() - t0);
    if (stop_after && strcmp(stop_after, "clean2") == 0) goto done;
    t0 = ves_clock_sec();
    if (sa_stage_place(&run, &cfg, out_dir) != 0) return 8;
    fprintf(stderr, "[time] place %.1f s\n", ves_clock_sec() - t0);
    {
        char checkpoint[2048];snprintf(checkpoint,sizeof checkpoint,"%s/stage5_state.asr",out_dir);
        if(AsmStore_write_run(checkpoint,&run))return 8;
    }
    if (stop_after && strcmp(stop_after, "place") == 0) goto done;

discover:
    /* Resume benchmarks include checkpoint loading and discovery. Reserve
     * 170 s for checkpoint/geometry output and the independent full audit. */
    if(cfg.repair.compact && cfg.repair.budget_seconds>0)
        cfg.repair.deadline=fmax(ves_clock_sec(),(resume_placement?t_start:ves_clock_sec())+cfg.repair.budget_seconds-170);
    t0 = ves_clock_sec();
    if (sa_stage_discover(&run, &cfg, out_dir) != 0) return 8;
    fprintf(stderr, "[time] discover %.1f s\n", ves_clock_sec() - t0);
    if (stop_after && strcmp(stop_after, "discover") == 0) goto done;

    t0 = ves_clock_sec();
    AsmRepairStats repair_stats;
    cfg.placeo.stress_band=cfg.clean.stress_band;
    cfg.repair.layer_place=&cfg.placeo;cfg.repair.layer_conflict=&cfg.conflict;
    cfg.repair.threads=cfg.threads;
    if (AsmRepair_run(&run,&cfg.repair,out_dir,&repair_stats)) {
        AsmAuditStats failed_audit; AsmAudit_run(&run,out_dir,&failed_audit);
        fprintf(stderr,"[assemble repair] failed; sheet export refused\n"); return 10;
    }
    fprintf(stderr,"[time] repair %.1f s\n",ves_clock_sec()-t0);
    {
        char png[2048]; snprintf(png,sizeof png,"%s/stage6_repaired.png",out_dir);
        Arena_Mark mark = Arena_save(run.arena);
        AsmChart *display = ARENA_ALLOC(run.arena,run.n_charts*sizeof(AsmChart));
        memcpy(display,run.charts,run.n_charts*sizeof(AsmChart));
        for (size_t c = 0; c < run.n_charts; c++) { display[c].placed = AsmChart_registered(run.charts+c); display[c].component = 0; }
        int png_rc = repair_stats.charts ? AsmReport_layout_png(run.arena,png,display,run.n_charts,NULL,0,cfg.conflict.cell,8192,1) : 0;
        Arena_restore(run.arena,mark);
        if (png_rc) return 10;
    }
    if (stop_after && strcmp(stop_after,"repair") == 0) goto done;

    t0 = ves_clock_sec();
    AsmAuditStats audit_stats;
    if (AsmAudit_run(&run,out_dir,&audit_stats)) return 11;
    fprintf(stderr,"[time] audit %.1f s\n",ves_clock_sec()-t0);
    if (stop_after && strcmp(stop_after,"audit") == 0) goto done;
    if (!audit_stats.source_preserved || audit_stats.invalid_faces) {
        fprintf(stderr,"[assemble audit] source/embedding invalid; sheet export refused\n"); return 11;
    }

    t0 = ves_clock_sec();
    if (sa_stage_sheet(&run, &cfg, out_dir) != 0) return 9;
    fprintf(stderr, "[time] sheet %.1f s\n", ves_clock_sec() - t0);

done:
    fprintf(stderr, "[time] total %.1f s\n", ves_clock_sec() - t_start);
    return 0;
}
