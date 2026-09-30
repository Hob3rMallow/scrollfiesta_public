/* sheet_layers -- the ink-model surface volume of a flattened sheet, on its bake's pixel grid.
 *
 *   sheet_layers <in.vmesh> <raw_dir|raw.zarr> <axis.csv> <out_dir> [--id NAME] [--geometry] [--apron]
 *                [--smooth] [--stepblend] [--v8in] [--window C0 C1]
 *   sheet_layers --selftest
 *
 * <in.vmesh> is a UV-carrying layout (sheet_assemble's stage6_repaired.vmesh), <raw> the RAW
 * source the review bake read (cubes_RAW TIFFs or an uncompressed uint8 Zarr v2, 128^3), and
 * <axis.csv> the z,y,x scroll-axis table that orients every chart's normals inward.  Writes
 *   <out_dir>/<id>_layers.zarr         uint8 [SHEET_LAYERS_DEPTH, H, W], layer index inward
 *   <out_dir>/<id>_layers_mask.tif     255 = painted and every layer fully supported by RAW
 *   <out_dir>/<id>_layers_report.json  grid (== the bake's sheet_rawtex_grid.json), orientation
 *   <out_dir>/<id>_handedness.json     per piece (run of painted columns): pixels showing the recto in
 *                                      the team's reading handedness vs mirror-imaged (RawtexHandedness)
 *   <out_dir>/<id>_geometry.zarr       (--geometry) float32 [6, H, W]: the depth-zero point
 *                                      (z,y,x) and unit inward normal behind every layer column,
 *                                      for carrying labels from other surfaces onto this grid
 * --apron renders a context-only apron (RawtexLayersOpts: holes and a SHEET_LAYERS_APRON_PX margin
 * beyond piece borders, on the sheet's own quadric continuation) into the layer volume; the mask
 * and the geometry volume are unchanged.
 * --smooth regularises the depth-zero surface before sampling (RawtexLayersOpts: a masked Gaussian of
 * SHEET_LAYERS_SMOOTH_SIGMA_PX over the sheet's points and normals; every pixel moves only along its
 * smoothed normal, by at most SHEET_LAYERS_SMOOTH_CAP_VOX): the mesh's few-pixel wobble tilts the
 * sampling axis and jitters depth zero.  The geometry volume records the regularised surface.
 * --stepblend bends the texture across 3-D steps between charts instead of tearing it (RawtexLayersOpts:
 * near steps, points move toward a SHEET_LAYERS_STEP_SIGMA_PX local linear fit, spread over
 * SHEET_LAYERS_STEP_SPREAD_PX).  The geometry volume records the blended surface.
 * --v8in renders the 24-layer contract of the ink-8um v8-in models (SHEET_LAYERS_V8IN_DEPTH: offsets
 * -11.5..+11.5 vox, layer index inward as ever) instead of the 21-layer reader-v2 volume.
 * --window C0 C1 renders only columns [C0, C1) of the sheet's full 1 px/UV grid (all its rows): the mesh is
 * cropped to the faces within SHEET_LAYERS_WINDOW_PAD_PX of those columns first (so only their RAW is
 * loaded), and pixel x of every output is column C0 + x of the full grid -- how a sheet wider than the
 * per-axis raster cap (a 21^3 layout at 1 px/vox) is rendered, window by window.
 * The grid is exactly obj_bake_raw's for the same mesh at 1 px per UV unit, so an ink map
 * inferred from the volume overlays sheet_rawtex.tif pixel for pixel.  Constants:
 * SHEET_LAYERS_* in pipeline_constants.h.  Threads: OMP_NUM_THREADS. */
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_normals.h"
#include "../common/pipeline_constants.h"
#include "../common/raw_sample.h"
#include "../common/ves_platform.h"
#include "../assemble/asm_axis.h"
#include "../flatten/rawtex_layers.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* JSON string literal (Windows paths carry backslashes) into buf. */
static void json_str(char *buf, size_t cap, const char *s)
{
    size_t o = 0;
    if (cap < 3) { if (cap) buf[0] = 0; return; }
    buf[o++] = '"';
    for (; *s && o + 8 < cap; s++) {
        if (*s == '"' || *s == '\\') { buf[o++] = '\\'; buf[o++] = *s; }
        else if ((unsigned char)*s < 0x20) o += (size_t)snprintf(buf + o, cap - o, "\\u%04x", (unsigned)(unsigned char)*s);
        else buf[o++] = *s;
    }
    buf[o++] = '"';
    buf[o] = 0;
}

/* Unit direction from each vertex to its axis foot (perpendicular to the axis); zero on the axis. */
static void sheet_inward(const AsmAxis *axis, const float *verts, size_t nv, float *inward)
{
    long i = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < (long)nv; i++) {
        double p[3], foot[3], d[3], l = 0.0;
        int k = 0;
        for (k = 0; k < 3; k++) p[k] = (double)verts[(size_t)i * 3 + (size_t)k];
        AsmAxis_project(axis, p, NULL, NULL, NULL, foot);
        for (k = 0; k < 3; k++) { d[k] = foot[k] - p[k]; l += d[k] * d[k]; }
        l = sqrt(l);
        for (k = 0; k < 3; k++) inward[(size_t)i * 3 + (size_t)k] = l > 1e-9 ? (float)(d[k] / l) : 0.0f;
    }
}

/* Keep the faces whose UV bounding box meets u in [ulo, uhi) and the vertices they use (compacted,
 * faces remapped); the replaced arrays are freed.  0 = ok, -1 = out of memory. */
static int sheet_crop_u(MeshBinData *m, double ulo, double uhi)
{
    int32_t *remap = (int32_t *)malloc((m->nv ? m->nv : 1) * sizeof *remap);
    uint8_t *keep = (uint8_t *)calloc(m->nf ? m->nf : 1, 1);
    float *nverts = NULL, *nuv = NULL;
    double *nuv64 = NULL;
    int32_t *nfaces = NULL;
    size_t f = 0, v = 0, kv = 0, kf = 0;
    if (remap == NULL || keep == NULL) { free(remap); free(keep); return -1; }
    for (v = 0; v < m->nv; v++) remap[v] = -1;
    for (f = 0; f < m->nf; f++) {
        double lo = 0.0, hi = 0.0;
        int k = 0;
        for (k = 0; k < 3; k++) {
            double u = MeshBin_uv(m, (size_t)m->faces[f * 3 + (size_t)k] * 2);
            if (k == 0 || u < lo) lo = u;
            if (k == 0 || u > hi) hi = u;
        }
        if (hi < ulo || lo >= uhi) continue;
        keep[f] = 1; kf++;
        for (k = 0; k < 3; k++) {
            size_t w = (size_t)m->faces[f * 3 + (size_t)k];
            if (remap[w] < 0) remap[w] = (int32_t)kv++;
        }
    }
    nverts = (float *)malloc((kv ? kv : 1) * 3 * sizeof *nverts);
    nfaces = (int32_t *)malloc((kf ? kf : 1) * 3 * sizeof *nfaces);
    if (m->uv64 != NULL) nuv64 = (double *)malloc((kv ? kv : 1) * 2 * sizeof *nuv64);
    else nuv = (float *)malloc((kv ? kv : 1) * 2 * sizeof *nuv);
    if (nverts == NULL || nfaces == NULL || (nuv64 == NULL && nuv == NULL)) {
        free(remap); free(keep); free(nverts); free(nfaces); free(nuv64); free(nuv);
        return -1;
    }
    for (v = 0; v < m->nv; v++) {
        size_t w = 0;
        if (remap[v] < 0) continue;
        w = (size_t)remap[v];
        memcpy(&nverts[w * 3], &m->verts[v * 3], 3 * sizeof *nverts);
        if (nuv64 != NULL) { nuv64[w * 2] = m->uv64[v * 2]; nuv64[w * 2 + 1] = m->uv64[v * 2 + 1]; }
        else { nuv[w * 2] = m->uv[v * 2]; nuv[w * 2 + 1] = m->uv[v * 2 + 1]; }
    }
    for (f = 0, kf = 0; f < m->nf; f++) {
        int k = 0;
        if (!keep[f]) continue;
        for (k = 0; k < 3; k++) nfaces[kf * 3 + (size_t)k] = remap[(size_t)m->faces[f * 3 + (size_t)k]];
        kf++;
    }
    free(m->verts); free(m->faces); free(m->uv); free(m->uv64);
    m->verts = nverts; m->faces = nfaces; m->uv = nuv; m->uv64 = nuv64; m->nv = kv; m->nf = kf;
    free(remap); free(keep);
    return 0;
}

/* Offset of layer k along the inward normal: layers centred on depth zero, SHEET_LAYERS_STEP_VOX apart. */
static double layer_offset(int k, int depth)
{
    return ((double)k - 0.5 * (double)(depth - 1)) * SHEET_LAYERS_STEP_VOX + SHEET_LAYERS_OFFSET_VOX;
}

static int selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    /* straight axis along +z through (y,x) = (0,0): a vertex at (5, 10, 0) looks toward -y */
    double point[3] = { 0.0, 0.0, 0.0 }, dir[3] = { 1.0, 0.0, 0.0 };
    AsmAxis *axis = AsmAxis_line(arena, point, dir, -100.0, 100.0);
    float v[6] = { 5.0f, 10.0f, 0.0f,   7.0f, -3.0f, 4.0f };
    float in[6];
    if (axis == NULL) { fprintf(stderr, "  FAIL: axis\n"); fails++; }
    else {
        sheet_inward(axis, v, 2, in);
        if (fabsf(in[0]) > 1e-6f || fabsf(in[1] + 1.0f) > 1e-6f || fabsf(in[2]) > 1e-6f) {
            fprintf(stderr, "  FAIL: inward of (5,10,0) is (%g,%g,%g), want (0,-1,0)\n", in[0], in[1], in[2]);
            fails++;
        }
        if (fabsf(in[3]) > 1e-6f || fabsf(in[4] - 0.6f) > 1e-6f || fabsf(in[5] + 0.8f) > 1e-6f) {
            fprintf(stderr, "  FAIL: inward of (7,-3,4) is (%g,%g,%g), want (0,0.6,-0.8)\n", in[3], in[4], in[5]);
            fails++;
        }
    }
    Arena_dispose(&arena);
    /* layer contracts: the default 21 layers sit at -10..+10, the v8-in 24 at -11.5..+11.5 */
    if (fabs(layer_offset(0, SHEET_LAYERS_DEPTH) + 10.0) > 1e-9 ||
        fabs(layer_offset(SHEET_LAYERS_DEPTH / 2, SHEET_LAYERS_DEPTH)) > 1e-9 ||
        fabs(layer_offset(0, SHEET_LAYERS_V8IN_DEPTH) + 11.5) > 1e-9 ||
        fabs(layer_offset(SHEET_LAYERS_V8IN_DEPTH - 1, SHEET_LAYERS_V8IN_DEPTH) - 11.5) > 1e-9 ||
        fabs(layer_offset(12, SHEET_LAYERS_V8IN_DEPTH) - 0.5) > 1e-9) {
        fprintf(stderr, "  FAIL: layer offsets (21: %g..%g, 24: %g..%g)\n", layer_offset(0, 21),
                layer_offset(20, 21), layer_offset(0, 24), layer_offset(23, 24));
        fails++;
    }
    {   /* window crop: of two triangles at u 0..1 and u 10..11, [5, 20) keeps the second, remapped */
        MeshBinData cm;
        memset(&cm, 0, sizeof cm);
        cm.nv = 6; cm.nf = 2;
        cm.verts = (float *)malloc(18 * sizeof(float));
        cm.uv64 = (double *)malloc(12 * sizeof(double));
        cm.faces = (int32_t *)malloc(6 * sizeof(int32_t));
        if (cm.verts == NULL || cm.uv64 == NULL || cm.faces == NULL) fails++;
        else {
            int k = 0;
            for (k = 0; k < 18; k++) cm.verts[k] = (float)k;
            for (k = 0; k < 6; k++) { cm.uv64[k * 2] = (k < 3 ? 0.0 : 10.0) + (k % 3 == 1); cm.uv64[k * 2 + 1] = k; }
            for (k = 0; k < 6; k++) cm.faces[k] = k;
            if (sheet_crop_u(&cm, 5.0, 20.0) != 0 || cm.nf != 1 || cm.nv != 3 || cm.faces[0] != 0 ||
                cm.faces[2] != 2 || cm.uv64[0] != 10.0 || cm.uv64[2] != 11.0 || cm.verts[0] != 9.0f) {
                fprintf(stderr, "  FAIL: window crop kept %zu faces / %zu verts\n", cm.nf, cm.nv);
                fails++;
            }
        }
        MeshBin_dispose(&cm);
    }
    fprintf(stderr, "[selftest] sheet_layers inward %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *raw_dir = NULL, *axis_path = NULL, *out_dir = NULL, *id = "sheet";
    char vmesh[2048], zarr[2400], mask[2400], report[2400], attrs[16384], geom[2400];
    char js_in[4200], js_raw[4200], js_axis[4200], js_zarr[4200], js_mask[4200], js_geom[4200];
    char smooth_js[512], hand_path[2400], hand_js[512], v8in_js[128], win_js[512];
    long win_c0 = 0, win_c1 = 0;
    int want_win = 0;
    RawtexWindow win;
    RawtexPlan full;
    size_t full_nv = 0, full_nf = 0;
    RawtexHandedness hand;
    size_t hand_team = 0, hand_mirror = 0, hand_mpieces = 0, hand_mpx = 0;
    int want_geom = 0, want_apron = 0, want_smooth = 0, want_step = 0, want_v8in = 0;
    int depth = SHEET_LAYERS_DEPTH;
    MeshBinData m;
    Arena_T arena = NULL;
    AsmAxis *axis = NULL;
    CubeTable ct;
    float *normals = NULL, *inward = NULL;
    uint8_t *face_skip = NULL, *reached = NULL;
    RawtexOrientStats os;
    RawtexLayersStats ls;
    RawtexLayersOpts opt;
    size_t faces_skipped = 0, reached_missing = 0, expected = 0;
    double reach = 0.0, t0 = ves_clock_sec(), t_load = 0.0, t_orient = 0.0, t_raw = 0.0;
    int i = 0, rc = 1;
    FILE *fp = NULL;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return RawtexLayers_selftest() == 0 && selftest() == 0 && AsmAxis_selftest() == 0 ? 0 : 1;
    if (argc < 5) {
        fprintf(stderr, "usage: sheet_layers <in.vmesh> <raw_dir|raw.zarr> <axis.csv> <out_dir> [--id NAME] "
                        "[--geometry] [--apron] [--smooth] [--stepblend] [--v8in] [--window C0 C1]\n"
                        "       sheet_layers --selftest\n");
        return 1;
    }
    in_path = argv[1]; raw_dir = argv[2]; axis_path = argv[3]; out_dir = argv[4];
    for (i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--id") == 0 && i + 1 < argc) id = argv[++i];
        else if (strcmp(argv[i], "--geometry") == 0) want_geom = 1;
        else if (strcmp(argv[i], "--apron") == 0) want_apron = 1;
        else if (strcmp(argv[i], "--smooth") == 0) want_smooth = 1;
        else if (strcmp(argv[i], "--stepblend") == 0) want_step = 1;
        else if (strcmp(argv[i], "--v8in") == 0) { want_v8in = 1; depth = SHEET_LAYERS_V8IN_DEPTH; }
        else if (strcmp(argv[i], "--window") == 0 && i + 2 < argc) {
            char *e0 = NULL, *e1 = NULL;
            win_c0 = strtol(argv[i + 1], &e0, 10); win_c1 = strtol(argv[i + 2], &e1, 10); i += 2;
            if (*e0 || *e1 || win_c0 < 0 || win_c1 <= win_c0) { fprintf(stderr, "bad --window\n"); return 1; }
            want_win = 1;
        }
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    memset(&m, 0, sizeof m);
    memset(&hand, 0, sizeof hand);
    if (MeshBin_companion_path(in_path, vmesh, sizeof vmesh) != 0 ||
        MeshBin_read_precise_malloc(vmesh, &m) != 0 || (m.uv == NULL && m.uv64 == NULL) ||
        m.nv == 0 || m.nf == 0) {
        fprintf(stderr, "ERROR: cannot read a UV-carrying mesh from %s\n", in_path);
        return 1;
    }
    for (size_t v = 0; v < m.nv * 3; v++)
        if (!isfinite(m.verts[v])) { fprintf(stderr, "ERROR: non-finite vertex coordinate\n"); return 1; }
    fprintf(stderr, "[sheet_layers] in=%s raw=%s axis=%s out=%s id=%s\n", vmesh, raw_dir, axis_path, out_dir, id);
    fprintf(stderr, "  loaded %zu verts, %zu faces\n", m.nv, m.nf);
    memset(&win, 0, sizeof win);
    memset(&full, 0, sizeof full);
    full_nv = m.nv; full_nf = m.nf;
    if (want_win) {       /* columns of the FULL grid; then only the faces near them */
        if (Rawtex_plan_field((RawtexUv){ m.uv64 ? NULL : m.uv, m.uv64 }, m.nv, 1.0, 1.0, (size_t)6000000000ULL,
                              NULL, &full) != 0 || (size_t)win_c1 > full.W) {
            fprintf(stderr, "ERROR: window [%ld, %ld) outside the sheet's %zu columns\n", win_c0, win_c1, full.W);
            MeshBin_dispose(&m);
            return 1;
        }
        win.umin = full.umin + (double)win_c0; win.umax = full.umin + (double)win_c1;
        win.vmin = full.vmin; win.vmax = full.vmin + (double)full.H;
        if (sheet_crop_u(&m, win.umin - SHEET_LAYERS_WINDOW_PAD_PX, win.umax + SHEET_LAYERS_WINDOW_PAD_PX) != 0 ||
            m.nf == 0) {
            fprintf(stderr, "ERROR: no faces in window [%ld, %ld) (or out of memory)\n", win_c0, win_c1);
            MeshBin_dispose(&m);
            return 1;
        }
        fprintf(stderr, "  window: columns [%ld, %ld) of the full %zu x %zu grid (origin u %.17g, v %.17g); "
                "kept %zu of %zu faces, %zu of %zu verts\n", win_c0, win_c1, full.W, full.H, full.umin, full.vmin,
                m.nf, full_nf, m.nv, full_nv);
    }
    arena = Arena_new();
    axis = AsmAxis_load(arena, axis_path);
    if (axis == NULL) { fprintf(stderr, "ERROR: cannot load axis table %s\n", axis_path); goto out; }
    t_load = ves_clock_sec() - t0;

    /* normals: winding normals, oriented inward per face-graph component, smoothed */
    t0 = ves_clock_sec();
    normals = MeshNormals_compute(m.verts, m.nv, m.faces, m.nf);
    inward = (float *)malloc(m.nv * 3 * sizeof *inward);
    if (normals == NULL || inward == NULL) { fprintf(stderr, "ERROR: out of memory\n"); goto out; }
    sheet_inward(axis, m.verts, m.nv, inward);
    if (RawtexLayers_orient(normals, m.verts, m.nv, m.faces, m.nf, inward, SHEET_LAYERS_NORMAL_SMOOTH, &os) != 0) {
        fprintf(stderr, "ERROR: normal orientation failed\n");
        goto out;
    }
    t_orient = ves_clock_sec() - t0;
    fprintf(stderr, "  orientation: %zu components, %zu flipped (%.1f%% of vote area), %zu weak; "
            "%.2f%% of vote area agrees with its component (%.1fs)\n", os.components, os.flipped,
            100.0 * os.flipped_area, os.weak, 100.0 * os.agree_area, t_orient);

    /* RAW: the table reaches every layer; faces the review bake left unpainted (a missing
     * chunk within its reach) are left unpainted here too */
    t0 = ves_clock_sec();
    reach = 0.5 * (double)(depth - 1) * SHEET_LAYERS_STEP_VOX + fabs(SHEET_LAYERS_OFFSET_VOX) + 2.0;
    /* an apron point may lie up to 1.6 x SHEET_LAYERS_APRON_PX from the mesh (the isometry gate) */
    if (want_apron) reach += 1.6 * SHEET_LAYERS_APRON_PX + 2.0;
    if (want_smooth) reach += SHEET_LAYERS_SMOOTH_CAP_VOX;     /* a regularised point moves at most the cap */
    if (want_step) reach += SHEET_LAYERS_STEP_REACH_VOX;       /* a blended point moves toward its neighbours' fit */
    if (cubetable_init(&ct, arena, raw_dir, SHEET_LAYERS_CHUNK, m.verts, m.nv, reach) != 0) {
        fprintf(stderr, "ERROR: cube table init failed for %s\n", raw_dir);
        goto out;
    }
    cubetable_prewarm_all(&ct);
    expected = cubetable_expected_chunks(&ct);
    if (!cubetable_is_complete(&ct) && expected > 0) {
        reached = (uint8_t *)ARENA_CALLOC(arena, expected, 1);
        reached_missing = cubetable_missing_reached(&ct, m.verts, m.faces, m.nf, SHEET_LAYERS_BAKE_REACH_VOX, reached);
        if (reached_missing > 0) {
            face_skip = (uint8_t *)calloc(m.nf, 1);
            if (face_skip == NULL) { fprintf(stderr, "ERROR: out of memory\n"); goto out; }
            faces_skipped = cubetable_faces_reaching(&ct, m.verts, m.faces, m.nf, SHEET_LAYERS_BAKE_REACH_VOX,
                                                     reached, face_skip);
        }
    }
    t_raw = ves_clock_sec() - t0;
    fprintf(stderr, "  RAW: %d chunks loaded, %d missing (%zu reached by a face, %zu faces unpainted), %.1fs\n",
            ct.n_loaded, ct.n_missing, reached_missing, faces_skipped, t_raw);

    memset(&opt, 0, sizeof opt);
    opt.depth = depth; opt.step = SHEET_LAYERS_STEP_VOX; opt.offset = SHEET_LAYERS_OFFSET_VOX;
    opt.tile_cols = SHEET_LAYERS_TILE_COLS; opt.chunk = SHEET_LAYERS_CHUNK;
    opt.stretch_ratio = SHEET_LAYERS_STRETCH_RATIO; opt.stretch_floor = SHEET_LAYERS_STRETCH_FLOOR;
    opt.seam_reach = SHEET_LAYERS_SEAM_REACH_PX; opt.seam_tol = SHEET_LAYERS_SEAM_TOL_VOX;
    opt.seam_passes = SHEET_LAYERS_SEAM_PASSES;
    snprintf(zarr, sizeof zarr, "%s/%s_layers.zarr", out_dir, id);
    snprintf(mask, sizeof mask, "%s/%s_layers_mask.tif", out_dir, id);
    snprintf(report, sizeof report, "%s/%s_layers_report.json", out_dir, id);
    snprintf(geom, sizeof geom, "%s/%s_geometry.zarr", out_dir, id);
    opt.geom_dir = want_geom ? geom : NULL;
    opt.hand = &hand;
    opt.window = want_win ? &win : NULL;
    snprintf(hand_path, sizeof hand_path, "%s/%s_handedness.json", out_dir, id);
    if (want_apron) {
        opt.apron_px = SHEET_LAYERS_APRON_PX; opt.apron_fit_px = SHEET_LAYERS_APRON_FIT_PX;
        opt.apron_cell_px = SHEET_LAYERS_APRON_CELL_PX; opt.apron_step_px = SHEET_LAYERS_APRON_STEP_PX;
    }
    if (want_smooth) {
        opt.smooth_sigma_px = SHEET_LAYERS_SMOOTH_SIGMA_PX;
        opt.smooth_cap_vox = SHEET_LAYERS_SMOOTH_CAP_VOX;
    }
    if (want_step) {
        opt.step_sigma_px = SHEET_LAYERS_STEP_SIGMA_PX;
        opt.step_spread_px = SHEET_LAYERS_STEP_SPREAD_PX;
    }
    json_str(js_in, sizeof js_in, vmesh);
    json_str(js_raw, sizeof js_raw, raw_dir);
    json_str(js_axis, sizeof js_axis, axis_path);
    if (want_geom) json_str(js_geom, sizeof js_geom, geom);
    else snprintf(js_geom, sizeof js_geom, "null");
    {
        size_t o = 0;
        int k = 0;
        o += (size_t)snprintf(attrs + o, sizeof attrs - o,
            "{\n  \"format\": \"scrollfiesta-sheet-layers-v1\",\n  \"mesh\": %s,\n  \"raw\": %s,\n"
            "  \"axis\": %s,\n  \"axis_order\": \"depth,y,x\",\n"
            "  \"layer_direction\": \"index increases inward, toward the scroll axis (recto side)\",\n"
            "  \"depth_zero\": \"the mesh surface (F0 recto-surface midline)\",\n"
            "  \"pixel_grid\": \"obj_bake_raw sheet_rawtex grid, 1 px per UV unit, centre = origin + index + 0.5\",\n"
            "  \"layer_offsets_voxels\": [", js_in, js_raw, js_axis);
        for (k = 0; k < depth; k++)
            o += (size_t)snprintf(attrs + o, sizeof attrs - o, "%s%.3f", k ? ", " : "", layer_offset(k, depth));
        o += (size_t)snprintf(attrs + o, sizeof attrs - o, "]");
        if (want_v8in)        /* default volumes keep their exact attributes */
            o += (size_t)snprintf(attrs + o, sizeof attrs - o,
                ",\n  \"layer_contract\": \"ink-8um v8-in: %d layers, offsets %.1f..%.1f vox along the inward "
                "normal\"", depth, layer_offset(0, depth), layer_offset(depth - 1, depth));
        if (want_apron)       /* default volumes keep their exact attributes */
            o += (size_t)snprintf(attrs + o, sizeof attrs - o,
                ",\n  \"context_apron\": {\"px\": %d, \"fit_px\": %d, \"cell_px\": %d, \"step_px\": %d, "
                "\"meaning\": \"pixels outside the mask carry context-only layers on a quadric continuation\"}",
                SHEET_LAYERS_APRON_PX, SHEET_LAYERS_APRON_FIT_PX, SHEET_LAYERS_APRON_CELL_PX,
                SHEET_LAYERS_APRON_STEP_PX);
        if (want_smooth)
            o += (size_t)snprintf(attrs + o, sizeof attrs - o,
                ",\n  \"surface_regularisation\": {\"sigma_px\": %.2f, \"cap_vox\": %.2f, "
                "\"meaning\": \"depth zero and normals low-passed over the sheet; pixels moved along the normal only\"}",
                SHEET_LAYERS_SMOOTH_SIGMA_PX, SHEET_LAYERS_SMOOTH_CAP_VOX);
        if (want_step)
            o += (size_t)snprintf(attrs + o, sizeof attrs - o,
                ",\n  \"step_blending\": {\"sigma_px\": %.2f, \"spread_px\": %.2f, "
                "\"meaning\": \"near 3-D steps between charts, points moved toward a local linear fit of their neighbours\"}",
                SHEET_LAYERS_STEP_SIGMA_PX, SHEET_LAYERS_STEP_SPREAD_PX);
        if (want_win)
            o += (size_t)snprintf(attrs + o, sizeof attrs - o,
                ",\n  \"window\": {\"cols\": [%ld, %ld], \"full_origin_uv\": [%.17g, %.17g], \"full_shape_hw\": [%zu, %zu], "
                "\"meaning\": \"pixel x is column cols[0] + x of the sheet's full 1 px/UV grid\"}",
                win_c0, win_c1, full.umin, full.vmin, full.H, full.W);
        o += (size_t)snprintf(attrs + o, sizeof attrs - o, "\n}\n");
    }
    if (RawtexLayers_write(zarr, mask, &ct, m.verts,
                           (RawtexUv){ m.uv64 ? NULL : m.uv, m.uv64 }, m.nv, m.faces, m.nf,
                           normals, face_skip, &opt, attrs, &ls) != 0) {
        fprintf(stderr, "ERROR: layer volume not written\n");
        goto out;
    }
    fprintf(stderr, "  wrote %s: %d x %zu x %zu, %zu painted px + %zu seam-filled, %zu complete (%.2f%%), "
            "%zu multi-cover, %zu uv-stretched faces skipped, %zu chunks (%.1fs)\n", zarr, depth, ls.H,
            ls.W, ls.painted_px, ls.filled_px, ls.complete_px,
            ls.painted_px + ls.filled_px
                ? 100.0 * (double)ls.complete_px / (double)(ls.painted_px + ls.filled_px) : 0.0,
            ls.multi_px, ls.skip_uv_faces, ls.chunks_written, ls.seconds);
    if (want_apron)
        fprintf(stderr, "  context apron: %zu px rendered outside the mask, %zu rejected, %zu quadric fits\n",
                ls.apron_px, ls.apron_rejected, ls.apron_fits);
    if (want_smooth)
        fprintf(stderr, "  surface regularisation: %zu px, normal shift RMS %.3f vox, %zu capped at %.1f vox, "
                "mean normal turn %.2f deg\n", ls.smooth_px, ls.smooth_shift_rms, ls.smooth_capped_px,
                SHEET_LAYERS_SMOOTH_CAP_VOX, ls.smooth_turn_mean);
    if (want_step)
        fprintf(stderr, "  step blending: %zu px moved, RMS move %.3f vox\n", ls.step_px, ls.step_move_rms);
    {   /* handedness: which pieces show the recto mirror-imaged */
        size_t k = 0;
        for (k = 0; k < hand.npieces; k++) {
            const RawtexPieceHand *hp = &hand.pieces[k];
            hand_team += hp->px_team;
            hand_mirror += hp->px_mirror;
            if (hp->px_mirror > hp->px_team) { hand_mpieces++; hand_mpx += hp->px_team + hp->px_mirror; }
        }
        fprintf(stderr, "  handedness: %zu pieces, %zu mirror-imaged (%.1f%% of judged px) vs the team's reading "
                "canvases\n", hand.npieces, hand_mpieces,
                hand_team + hand_mirror ? 100.0 * (double)hand_mpx / (double)(hand_team + hand_mirror) : 0.0);
        if (ves_ensure_parent_dir(hand_path) != 0 || (fp = fopen(hand_path, "wb")) == NULL) {
            fprintf(stderr, "ERROR: cannot write %s\n", hand_path);
            goto out;
        }
        fprintf(fp, "{\n  \"schema\": \"sheet-handedness-v1\",\n"
                "  \"convention\": \"s = (dP/dcol x dP/drow) . n_inward, (z,y,x) component order; s > 0 = the "
                "handedness of the team's VC3D reading canvases, s < 0 = the recto mirror-imaged\",\n"
                "  \"grid_hw\": [%zu, %zu],\n  \"px_team\": %zu,\n  \"px_mirror\": %zu,\n"
                "  \"mirrored_pieces\": %zu,\n  \"piece_fields\": [\"col0\", \"col1\", \"row0\", \"row1\", "
                "\"px_team\", \"px_mirror\"],\n  \"pieces\": [", ls.H, ls.W, hand_team, hand_mirror,
                hand_mpieces);
        for (k = 0; k < hand.npieces; k++) {
            const RawtexPieceHand *hp = &hand.pieces[k];
            fprintf(fp, "%s\n    [%zu, %zu, %zu, %zu, %zu, %zu]", k ? "," : "", hp->col0, hp->col1, hp->row0,
                    hp->row1, hp->px_team, hp->px_mirror);
        }
        fprintf(fp, "\n  ]\n}\n");
        if (fclose(fp) != 0) { fp = NULL; fprintf(stderr, "ERROR: cannot write %s\n", hand_path); goto out; }
        fp = NULL;
        snprintf(hand_js, sizeof hand_js,
                 "  \"handedness\": {\"pieces\": %zu, \"mirrored_pieces\": %zu, \"px_team\": %zu, "
                 "\"px_mirror\": %zu},\n", hand.npieces, hand_mpieces, hand_team, hand_mirror);
    }
    json_str(js_zarr, sizeof js_zarr, zarr);
    json_str(js_mask, sizeof js_mask, mask);
    smooth_js[0] = 0;
    v8in_js[0] = 0;
    win_js[0] = 0;
    if (want_win)          /* default reports keep their bytes */
        snprintf(win_js, sizeof win_js,
                 "  \"window\": {\"cols\": [%ld, %ld], \"full_origin_uv\": [%.17g, %.17g], \"full_shape_hw\": "
                 "[%zu, %zu], \"faces_kept\": %zu, \"faces_full\": %zu},\n", win_c0, win_c1, full.umin, full.vmin,
                 full.H, full.W, m.nf, full_nf);
    if (want_v8in)         /* default reports keep their bytes */
        snprintf(v8in_js, sizeof v8in_js, "  \"layer_contract\": \"ink-8um-v8in-%d\",\n", depth);
    if (want_smooth)       /* default reports keep their bytes */
        snprintf(smooth_js, sizeof smooth_js,
                 "  \"smooth\": {\"sigma_px\": %.2f, \"cap_vox\": %.2f, \"px\": %zu, \"shift_rms_vox\": %.4f, "
                 "\"capped_px\": %zu, \"turn_mean_deg\": %.3f},\n", SHEET_LAYERS_SMOOTH_SIGMA_PX,
                 SHEET_LAYERS_SMOOTH_CAP_VOX, ls.smooth_px, ls.smooth_shift_rms, ls.smooth_capped_px,
                 ls.smooth_turn_mean);
    if (ves_ensure_parent_dir(report) != 0 || (fp = fopen(report, "wb")) == NULL) {
        fprintf(stderr, "ERROR: cannot write %s\n", report);
        goto out;
    }
    fprintf(fp,
        "{\n  \"schema\": \"sheet-layers-report-v1\",\n  \"mesh\": %s,\n  \"raw\": %s,\n  \"axis\": %s,\n"
        "  \"grid\": {\"schema\": \"vesuvius-rawtex-grid-v1\", \"coordinate_frame\": \"absolute_uv\", "
        "\"origin_uv\": [%.17g, %.17g], \"step_uv\": [1, 1], \"shape_hw\": [%zu, %zu], "
        "\"pixel_center_rule\": \"origin + (index + 0.5) * step\"},\n"
        "  \"depth\": %d,\n  \"step_vox\": %.3f,\n  \"offset_vox\": %.3f,\n  \"normal_smooth_iters\": %d,\n"
        "  \"orientation\": {\"components\": %zu, \"flipped\": %zu, \"weak\": %zu, \"agree_area\": %.6f, "
        "\"flipped_area\": %.6f},\n"
        "  \"raw_chunks_loaded\": %d,\n  \"raw_chunks_missing\": %d,\n  \"missing_chunks_reached\": %zu,\n"
        "  \"faces\": %zu,\n  \"faces_unpainted_missing_raw\": %zu,\n  \"faces_skipped_uv_stretch\": %zu,\n"
        "  \"seam\": {\"reach_px\": %d, \"tol_vox\": %.2f, \"passes\": %d},\n"
        "  \"painted_px\": %zu,\n  \"seam_filled_px\": %zu,\n  \"complete_px\": %zu,\n  \"multi_cover_px\": %zu,\n"
        "  \"chunks_written\": %zu,\n"
        "  \"apron\": {\"on\": %s, \"px\": %d, \"fit_px\": %d, \"rendered_px\": %zu, \"rejected_px\": %zu, "
        "\"fits\": %zu},\n%s%s%s%s"
        "  \"zarr\": %s,\n  \"mask\": %s,\n  \"geometry\": %s,\n"
        "  \"seconds\": {\"load\": %.2f, \"orient\": %.2f, \"raw\": %.2f, \"render\": %.2f}\n}\n",
        js_in, js_raw, js_axis, ls.umin, ls.vmin, ls.H, ls.W,
        depth, SHEET_LAYERS_STEP_VOX, SHEET_LAYERS_OFFSET_VOX, SHEET_LAYERS_NORMAL_SMOOTH,
        os.components, os.flipped, os.weak, os.agree_area, os.flipped_area,
        ct.n_loaded, ct.n_missing, reached_missing, m.nf, faces_skipped, ls.skip_uv_faces,
        SHEET_LAYERS_SEAM_REACH_PX, SHEET_LAYERS_SEAM_TOL_VOX, SHEET_LAYERS_SEAM_PASSES,
        ls.painted_px, ls.filled_px, ls.complete_px, ls.multi_px, ls.chunks_written,
        want_apron ? "true" : "false", want_apron ? SHEET_LAYERS_APRON_PX : 0,
        want_apron ? SHEET_LAYERS_APRON_FIT_PX : 0, ls.apron_px, ls.apron_rejected, ls.apron_fits, smooth_js, v8in_js,
        hand_js, win_js,
        js_zarr, js_mask, js_geom,
        t_load, t_orient, t_raw, ls.seconds);
    if (fclose(fp) != 0) { fprintf(stderr, "ERROR: cannot write %s\n", report); goto out; }
    fprintf(stderr, "[sheet_layers] SUMMARY id=%s grid=%zux%zu painted=%zu filled=%zu complete=%zu flipped=%zu/%zu "
            "agree=%.4f\n", id, ls.W, ls.H, ls.painted_px, ls.filled_px, ls.complete_px, os.flipped, os.components,
            os.agree_area);
    rc = 0;
out:
    free(normals); free(inward); free(face_skip);
    RawtexHandedness_free(&hand);
    MeshBin_dispose(&m);
    if (arena != NULL) Arena_dispose(&arena);
    return rc;
}
