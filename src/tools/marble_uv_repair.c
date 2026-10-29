/* marble_uv_repair.c -- transactional StrokeStrip UV repair driver.
 *
 * The matching RAW bake and detector mask select a low-trust part of an
 * already structured quad ribbon.  XYZ, topology, and axial V remain frozen;
 * MarbleStripUv solves only the common, monotone material coordinate U.
 */

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../flatten/marble_strip_uv.h"

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct GridInfo {
    int H, W;
    double umin, umax, vmin, vmax;
    double du, dv;
} GridInfo;

typedef struct MaskImage {
    int W, H;
    uint8_t *pixels;
} MaskImage;

typedef struct Crop {
    int r0, r1, c0, c1, H, W;
    float *verts;
    double *reference_uv;
    double *trust;
    uint8_t *pin;
    size_t core_vertices, movable_vertices;
} Crop;

static void usage(void)
{
    fprintf(stderr,
        "usage: marble_uv_repair <input.vmesh> <rawtex.tif> <mask.pgm> "
        "<output.vmesh> [options]\n"
        "The rawtex must have a strict, complete _raw_coverage.json companion,\n"
        "and the mask must be its matching _marble_mask.pgm.\n"
        "options:\n"
        "  --rect U0 V0 U1 V1      inclusive UV rectangle (required)\n"
        "  --feather F              movable halo outside detections (default 16)\n"
        "  --halo F                 additional pinned context (default 128)\n"
        "  --stride-u N             coarse material stride (default 4)\n"
        "  --stride-v N             coarse axial-row stride (default 4)\n"
        "  --coarsest N             compatibility alias: set both strides\n"
        "  --l1-rounds N            robust L1 initialization rounds (default 8)\n"
        "  --iterations N           total relaxed global-solve cap (default 1000)\n"
        "  --sweeps N               compatibility alias for --iterations\n"
        "  --final-iterations N     optional unknown-cross-section stage (default 0)\n"
        "  --final-damping F        local/global damping (default .1)\n"
        "  --tolerance F            max U movement stop (default 1e-5)\n"
        "  --sigma F                adjacency sigma; <=0 uses span/30\n"
        "  --lambda-length F        shared arc-length weight (default 1)\n"
        "  --lambda-align F         cross-section weight (default 1)\n"
        "  --lambda F               compatibility alias for --lambda-align\n"
        "  --lambda-local F         weak per-row metric weight (default .05)\n"
        "  --anchor-weight F        trusted-chart prior weight (default 100)\n"
        "  --monotone-fraction F    minimum U/arc speed (default .05)\n"
        "  --membership-floor F     detector-core prior (default .05)\n"
        "  --match-radius F         lateral XYZ adjacency reach (default 12)\n"
        "  --match-angle F          tangent gate in degrees (default 35)\n"
        "  --match-columns N        target-row search radius (default 24)\n"
        "  --topology-fallback F    same-column alternative prior (default .02)\n"
        "  --length-winsor F        bad/latent edge cap (default 3)\n"
        "  --max-disp F             transactional U trust radius (default 64)\n"
        "  --report PATH            JSON report (default <output>.json)\n"
        "  --verbose                print global-solve progress\n"
        "  --selftest               run variational and VMESH tests\n");
}

static int pgm_token(FILE *f, char *token, size_t capacity)
{
    int ch;
    size_t n = 0;
    if (f == NULL || token == NULL || capacity < 2) return -1;
    do {
        ch = fgetc(f);
        if (ch == '#') do ch = fgetc(f); while (ch != '\n' && ch != EOF);
    } while (ch != EOF && isspace((unsigned char)ch));
    if (ch == EOF) return -1;
    do {
        if (n + 1 >= capacity) return -1;
        token[n++] = (char)ch;
        ch = fgetc(f);
    } while (ch != EOF && !isspace((unsigned char)ch));
    token[n] = '\0';
    return 0;
}

static int mask_read_pgm(const char *path, MaskImage *mask)
{
    FILE *f = NULL;
    char token[64];
    long W, H, maximum;
    size_t n;
    memset(mask, 0, sizeof *mask);
    f = fopen(path, "rb");
    if (f == NULL) return -1;
    if (pgm_token(f, token, sizeof token) != 0 || strcmp(token, "P5") != 0 ||
        pgm_token(f, token, sizeof token) != 0 || (W = strtol(token, NULL, 10)) < 1 ||
        pgm_token(f, token, sizeof token) != 0 || (H = strtol(token, NULL, 10)) < 1 ||
        pgm_token(f, token, sizeof token) != 0 ||
        (maximum = strtol(token, NULL, 10)) != 255 || W > INT32_MAX ||
        H > INT32_MAX || (size_t)W > SIZE_MAX / (size_t)H) {
        fclose(f); return -1;
    }
    n = (size_t)W * (size_t)H;
    mask->pixels = (uint8_t *)malloc(n);
    if (mask->pixels == NULL || fread(mask->pixels, 1, n, f) != n ||
        fclose(f) != 0) {
        free(mask->pixels); memset(mask, 0, sizeof *mask); return -1;
    }
    mask->W = (int)W; mask->H = (int)H;
    return 0;
}

static void mask_dispose(MaskImage *mask)
{
    if (mask == NULL) return;
    free(mask->pixels); memset(mask, 0, sizeof *mask);
}

static void crop_dispose(Crop *crop)
{
    if (crop == NULL) return;
    free(crop->verts); free(crop->reference_uv);
    free(crop->trust); free(crop->pin);
    memset(crop, 0, sizeof *crop);
}

static int ends_with(const char *text, const char *suffix)
{
    size_t n = strlen(text), m = strlen(suffix);
    return n >= m && strcmp(text + n - m, suffix) == 0;
}

static char *read_small_file(const char *path)
{
    FILE *f = fopen(path, "rb");
    long length;
    char *text;
    if (f == NULL || fseek(f, 0, SEEK_END) != 0 ||
        (length = ftell(f)) < 0 || length > 1024 * 1024 ||
        fseek(f, 0, SEEK_SET) != 0) {
        if (f != NULL) fclose(f);
        return NULL;
    }
    text = (char *)malloc((size_t)length + 1);
    if (text == NULL || fread(text, 1, (size_t)length, f) != (size_t)length ||
        fclose(f) != 0) {
        free(text); return NULL;
    }
    text[length] = '\0';
    return text;
}

static const char *json_value(const char *text, const char *key)
{
    const char *p = strstr(text, key);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(key), ':');
    if (p == NULL) return NULL;
    do p++; while (*p != '\0' && isspace((unsigned char)*p));
    return p;
}

/* Refuse a visually plausible but unrelated bake.  Filename pairing catches
 * rawtex/mask swaps; the coverage companion catches partial/wrong RAW stores. */
static int validate_bake_contract(const char *rawtex, const char *mask,
                                  char *coverage, size_t capacity)
{
    static const char suffix[] = "_rawtex.tif";
    char expected_mask[2600];
    char *json = NULL;
    const char *value;
    size_t prefix;
    int ok = 0;
    if (!ends_with(rawtex, suffix)) return -1;
    prefix = strlen(rawtex) - (sizeof suffix - 1);
    if (prefix > INT32_MAX ||
        snprintf(expected_mask, sizeof expected_mask, "%.*s_marble_mask.pgm",
                 (int)prefix, rawtex) >= (int)sizeof expected_mask ||
        snprintf(coverage, capacity, "%.*s_raw_coverage.json",
                 (int)prefix, rawtex) >= (int)capacity ||
        strcmp(expected_mask, mask) != 0)
        return -1;
    json = read_small_file(coverage);
    if (json == NULL) return -1;
    value = json_value(json, "\"source_kind\"");
    if (value == NULL || strncmp(value, "\"zarr_v2\"", 9) != 0) goto done;
    value = json_value(json, "\"complete\"");
    if (value == NULL || strncmp(value, "true", 4) != 0) goto done;
    value = json_value(json, "\"strict\"");
    if (value == NULL || strncmp(value, "true", 4) != 0) goto done;
    value = json_value(json, "\"missing_chunks\"");
    if (value == NULL || strtol(value, NULL, 10) != 0) goto done;
    ok = 1;
done:
    free(json);
    return ok ? 0 : -1;
}

/* Infer dimensions from canonical row-major topology, not from an assumed
 * regular UV lattice.  This remains valid after the first repair/refit. */
static int structured_grid(const MeshBinData *mesh, GridInfo *grid)
{
    size_t W, H;
    double previous_v = -DBL_MAX;
    if (mesh == NULL || mesh->verts == NULL || mesh->uv == NULL ||
        mesh->faces == NULL || mesh->nv < 4 || mesh->nf < 2 ||
        mesh->faces[0] != 0 || mesh->faces[1] != 1 || mesh->faces[2] < 3)
        return -1;
    W = (size_t)mesh->faces[2] - 1u;
    if (W < 2 || W > INT32_MAX || mesh->nv % W != 0) return -1;
    H = mesh->nv / W;
    if (H < 2 || H > INT32_MAX || mesh->nf != 2u * (H - 1u) * (W - 1u))
        return -1;
    for (size_t r = 0; r + 1 < H; r++) for (size_t c = 0; c + 1 < W; c++) {
        size_t cell = r * (W - 1u) + c, f = cell * 6u;
        int32_t a = (int32_t)(r * W + c), b = a + 1;
        int32_t q = a + (int32_t)W, d = q + 1;
        if (mesh->faces[f] != a || mesh->faces[f + 1] != b ||
            mesh->faces[f + 2] != d || mesh->faces[f + 3] != a ||
            mesh->faces[f + 4] != d || mesh->faces[f + 5] != q)
            return -1;
    }
    grid->umin = grid->umax = mesh->uv[0];
    grid->vmin = grid->vmax = mesh->uv[1];
    for (size_t r = 0; r < H; r++) {
        double row_v = mesh->uv[(r * W) * 2 + 1];
        if (!isfinite(row_v) || (r > 0 && row_v <= previous_v + 1e-7)) return -1;
        previous_v = row_v;
        for (size_t c = 0; c < W; c++) {
            size_t i = r * W + c;
            double u = mesh->uv[i * 2], v = mesh->uv[i * 2 + 1];
            if (!isfinite(u) || !isfinite(v) || fabs(v - row_v) > 1e-3 ||
                (c > 0 && u <= (double)mesh->uv[(i - 1) * 2] + 1e-7))
                return -1;
            if (u < grid->umin) grid->umin = u;
            if (u > grid->umax) grid->umax = u;
            if (v < grid->vmin) grid->vmin = v;
            if (v > grid->vmax) grid->vmax = v;
        }
    }
    grid->H = (int)H; grid->W = (int)W;
    grid->du = (grid->umax - grid->umin) / (double)(W - 1u);
    grid->dv = (grid->vmax - grid->vmin) / (double)(H - 1u);
    return grid->du > 0.0 && grid->dv > 0.0 ? 0 : -1;
}

static int mask_at_uv(const MaskImage *mask, const GridInfo *grid,
                      double raster_du, double raster_dv, double u, double v)
{
    int x = (int)floor((u - grid->umin) / raster_du);
    int y = (int)floor((v - grid->vmin) / raster_dv);
    if (x < 0) x = 0; else if (x >= mask->W) x = mask->W - 1;
    if (y < 0) y = 0; else if (y >= mask->H) y = mask->H - 1;
    return mask->pixels[(size_t)y * (size_t)mask->W + (size_t)x] != 0;
}

static int crop_build(const MeshBinData *mesh, const GridInfo *grid,
                      const MaskImage *mask, const double rect[4],
                      double raster_du, double raster_dv,
                      double feather, double halo, Crop *crop)
{
    int core_r0 = grid->H, core_r1 = -1, core_c0 = grid->W, core_c1 = -1;
    int pad_r, pad_c;
    double *distance = NULL;
    const double inf = DBL_MAX / 16.0;
    const double diagonal = hypot(grid->du, grid->dv);
    size_t n;
    memset(crop, 0, sizeof *crop);

    for (int r = 0; r < grid->H; r++) for (int c = 0; c < grid->W; c++) {
        size_t i = (size_t)r * (size_t)grid->W + (size_t)c;
        double u = mesh->uv[i * 2], v = mesh->uv[i * 2 + 1];
        if (u >= rect[0] && u <= rect[2] && v >= rect[1] && v <= rect[3] &&
            mask_at_uv(mask, grid, raster_du, raster_dv, u, v)) {
            if (r < core_r0) core_r0 = r; if (r > core_r1) core_r1 = r;
            if (c < core_c0) core_c0 = c; if (c > core_c1) core_c1 = c;
        }
    }
    if (core_r1 < core_r0 || core_c1 < core_c0) return -1;
    pad_r = (int)ceil((feather + halo) / grid->dv);
    pad_c = (int)ceil((feather + halo) / grid->du);
    crop->r0 = core_r0 - pad_r; if (crop->r0 < 0) crop->r0 = 0;
    crop->r1 = core_r1 + pad_r; if (crop->r1 >= grid->H) crop->r1 = grid->H - 1;
    crop->c0 = core_c0 - pad_c; if (crop->c0 < 0) crop->c0 = 0;
    crop->c1 = core_c1 + pad_c; if (crop->c1 >= grid->W) crop->c1 = grid->W - 1;
    crop->H = crop->r1 - crop->r0 + 1;
    crop->W = crop->c1 - crop->c0 + 1;
    n = (size_t)crop->H * (size_t)crop->W;
    if (crop->H < 2 || crop->W < 2 || n > SIZE_MAX / (3 * sizeof(float)) ||
        n > SIZE_MAX / (2 * sizeof(double))) return -1;
    crop->verts = (float *)malloc(n * 3 * sizeof *crop->verts);
    crop->reference_uv = (double *)malloc(n * 2 * sizeof *crop->reference_uv);
    crop->trust = (double *)malloc(n * sizeof *crop->trust);
    crop->pin = (uint8_t *)malloc(n);
    distance = (double *)malloc(n * sizeof *distance);
    if (crop->verts == NULL || crop->reference_uv == NULL ||
        crop->trust == NULL || crop->pin == NULL || distance == NULL)
        goto fail;
    for (int r = 0; r < crop->H; r++) for (int c = 0; c < crop->W; c++) {
        int gr = crop->r0 + r, gc = crop->c0 + c;
        size_t gi = (size_t)gr * (size_t)grid->W + (size_t)gc;
        size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
        double u = mesh->uv[gi * 2], v = mesh->uv[gi * 2 + 1];
        int core = u >= rect[0] && u <= rect[2] && v >= rect[1] && v <= rect[3] &&
                   mask_at_uv(mask, grid, raster_du, raster_dv, u, v);
        memcpy(crop->verts + i * 3, mesh->verts + gi * 3, 3 * sizeof(float));
        crop->reference_uv[i * 2] = u;
        crop->reference_uv[i * 2 + 1] = v;
        distance[i] = core ? 0.0 : inf;
        if (core) crop->core_vertices++;
    }
#define RELAX_DIST(I,CANDIDATE) do { double _q=(CANDIDATE); \
    if (_q < distance[(I)]) distance[(I)] = _q; } while (0)
    for (int r = 0; r < crop->H; r++) for (int c = 0; c < crop->W; c++) {
        size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
        if (c > 0) RELAX_DIST(i, distance[i - 1] + grid->du);
        if (r > 0) {
            RELAX_DIST(i, distance[i - (size_t)crop->W] + grid->dv);
            if (c > 0) RELAX_DIST(i, distance[i - (size_t)crop->W - 1] + diagonal);
            if (c + 1 < crop->W)
                RELAX_DIST(i, distance[i - (size_t)crop->W + 1] + diagonal);
        }
    }
    for (int r = crop->H - 1; r >= 0; r--) for (int c = crop->W - 1; c >= 0; c--) {
        size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
        if (c + 1 < crop->W) RELAX_DIST(i, distance[i + 1] + grid->du);
        if (r + 1 < crop->H) {
            RELAX_DIST(i, distance[i + (size_t)crop->W] + grid->dv);
            if (c > 0) RELAX_DIST(i, distance[i + (size_t)crop->W - 1] + diagonal);
            if (c + 1 < crop->W)
                RELAX_DIST(i, distance[i + (size_t)crop->W + 1] + diagonal);
        }
    }
#undef RELAX_DIST
    for (int r = 0; r < crop->H; r++) for (int c = 0; c < crop->W; c++) {
        size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
        int boundary = r == 0 || c == 0 || r + 1 == crop->H || c + 1 == crop->W;
        int active = distance[i] <= feather;
        double t = feather > 0.0 ? distance[i] / feather : (distance[i] == 0.0 ? 0.0 : 1.0);
        if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
        crop->trust[i] = t * t * (3.0 - 2.0 * t);
        crop->pin[i] = (uint8_t)(boundary || !active);
        if (!crop->pin[i]) crop->movable_vertices++;
    }
    free(distance);
    if (crop->core_vertices == 0 || crop->movable_vertices == 0) goto fail_after_distance;
    return 0;
fail:
    free(distance);
fail_after_distance:
    crop_dispose(crop); return -1;
}

static int commit_crop(MeshBinData *mesh, const GridInfo *grid,
                       const Crop *crop, const double *uv)
{
    for (int r = 0; r < crop->H; r++) for (int c = 0; c < crop->W; c++) {
        size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
        size_t gi = (size_t)(crop->r0 + r) * (size_t)grid->W +
                    (size_t)(crop->c0 + c);
        mesh->uv[gi * 2] = (float)uv[i * 2];
        mesh->uv[gi * 2 + 1] = (float)uv[i * 2 + 1];
    }
    GridInfo check;
    return structured_grid(mesh, &check);
}

static void json_string(FILE *f, const char *text)
{
    fputc('"', f);
    for (; *text; text++) {
        unsigned char ch = (unsigned char)*text;
        if (ch == '"' || ch == '\\') { fputc('\\', f); fputc(ch, f); }
        else if (ch == '\n') fputs("\\n", f);
        else if (ch == '\r') fputs("\\r", f);
        else if (ch == '\t') fputs("\\t", f);
        else if (ch < 32) fprintf(f, "\\u%04x", (unsigned)ch);
        else fputc(ch, f);
    }
    fputc('"', f);
}

int main(int argc, char **argv)
{
    const char *input, *rawtex_path, *mask_path, *output, *report_path = NULL;
    char resolved[2600], default_report[2600], coverage_path[2600];
    double rect[4] = {0,0,0,0}, feather = 16.0, halo = 128.0;
    int have_rect = 0;
    MarbleStripUvOptions options;
    MarbleStripUvStats stats;
    MeshBinData mesh;
    GridInfo grid;
    MaskImage mask;
    Crop crop;
    Arena_T image_arena = NULL, solve_arena = NULL;
    uint8_t *rawtex = NULL;
    double *solved_uv = NULL;
    int depth = 0, texture_h = 0, texture_w = 0;
    double raster_du, raster_dv, started, seconds;

    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return MarbleStripUv_selftest() + MeshBin_selftest() == 0 ? 0 : 1;
    if (argc < 5) { usage(); return 2; }
    input = argv[1]; rawtex_path = argv[2]; mask_path = argv[3]; output = argv[4];
    MarbleStripUvOptions_default(&options);
    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--rect") == 0 && i + 4 < argc) {
            for (int k = 0; k < 4; k++) rect[k] = strtod(argv[++i], NULL);
            have_rect = 1;
        } else if (strcmp(argv[i], "--feather") == 0 && i + 1 < argc)
            feather = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--halo") == 0 && i + 1 < argc)
            halo = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--stride-u") == 0 && i + 1 < argc)
            options.sample_stride_u = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--stride-v") == 0 && i + 1 < argc)
            options.sample_stride_v = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--coarsest") == 0 && i + 1 < argc) {
            int n = (int)strtol(argv[++i], NULL, 10);
            options.sample_stride_u = options.sample_stride_v = n;
        } else if (strcmp(argv[i], "--l1-rounds") == 0 && i + 1 < argc)
            options.l1_iterations = (int)strtol(argv[++i], NULL, 10);
        else if ((strcmp(argv[i], "--iterations") == 0 ||
                  strcmp(argv[i], "--sweeps") == 0) && i + 1 < argc)
            options.max_irls_iterations = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--final-iterations") == 0 && i + 1 < argc)
            options.final_iterations = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--final-damping") == 0 && i + 1 < argc)
            options.final_damping = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--tolerance") == 0 && i + 1 < argc)
            options.movement_tolerance = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--sigma") == 0 && i + 1 < argc)
            options.likelihood_sigma = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--lambda-length") == 0 && i + 1 < argc)
            options.lambda_length = strtod(argv[++i], NULL);
        else if ((strcmp(argv[i], "--lambda-align") == 0 ||
                  strcmp(argv[i], "--lambda") == 0) && i + 1 < argc)
            options.lambda_align = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--lambda-local") == 0 && i + 1 < argc)
            options.lambda_local = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--anchor-weight") == 0 && i + 1 < argc)
            options.lambda_prior = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--monotone-fraction") == 0 && i + 1 < argc)
            options.monotone_fraction = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--membership-floor") == 0 && i + 1 < argc)
            options.membership_floor = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--match-radius") == 0 && i + 1 < argc)
            options.match_radius = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--match-angle") == 0 && i + 1 < argc)
            options.match_angle_deg = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--match-columns") == 0 && i + 1 < argc)
            options.match_search_columns = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--topology-fallback") == 0 && i + 1 < argc)
            options.topology_fallback = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--length-winsor") == 0 && i + 1 < argc)
            options.length_winsor = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--max-disp") == 0 && i + 1 < argc)
            options.max_displacement = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc)
            report_path = argv[++i];
        else if (strcmp(argv[i], "--verbose") == 0)
            options.verbose = 1;
        else {
            fprintf(stderr, "unknown or incomplete option: %s\n", argv[i]);
            usage(); return 2;
        }
    }
    if (!have_rect || !isfinite(rect[0]) || !isfinite(rect[1]) ||
        !isfinite(rect[2]) || !isfinite(rect[3]) || rect[0] > rect[2] ||
        rect[1] > rect[3] || !(feather > 0.0) || halo < 0.0 ||
        strcmp(input, output) == 0) {
        usage(); return 2;
    }
    if (report_path == NULL) {
        if (snprintf(default_report, sizeof default_report, "%s.json", output) >=
            (int)sizeof default_report) return 2;
        report_path = default_report;
    }
    setvbuf(stderr, NULL, _IONBF, 0);
    memset(&mesh, 0, sizeof mesh); memset(&mask, 0, sizeof mask);
    memset(&crop, 0, sizeof crop); memset(&stats, 0, sizeof stats);

    if (validate_bake_contract(rawtex_path, mask_path, coverage_path,
                               sizeof coverage_path) != 0) {
        fprintf(stderr, "ERROR: rawtex/mask/strict-coverage contract failed; "
                        "refusing an unpaired or incomplete RAW bake\n");
        goto fail;
    }
    if (MeshBin_companion_path(input, resolved, sizeof resolved) != 0 ||
        MeshBin_read_malloc(resolved, &mesh) != 0 ||
        structured_grid(&mesh, &grid) != 0) {
        fprintf(stderr, "ERROR: input is not a canonical monotone row-major ribbon\n");
        goto fail;
    }
    if (mask_read_pgm(mask_path, &mask) != 0) {
        fprintf(stderr, "ERROR: cannot read binary P5/255 mask %s\n", mask_path);
        goto fail;
    }
    image_arena = Arena_new();
    if (image_arena == NULL ||
        TiffIO_load(image_arena, rawtex_path, &rawtex, &depth,
                    &texture_h, &texture_w) != 0 || depth != 1 ||
        mask.W != texture_w || mask.H != texture_h) {
        fprintf(stderr, "ERROR: matching mask/rawtex raster contract failed\n");
        goto fail;
    }
    raster_du = (grid.umax - grid.umin) / (double)texture_w;
    raster_dv = (grid.vmax - grid.vmin) / (double)texture_h;
    if (!(raster_du > 0.0) || !(raster_dv > 0.0) ||
        crop_build(&mesh, &grid, &mask, rect, raster_du, raster_dv,
                   feather, halo, &crop) != 0) {
        fprintf(stderr, "ERROR: detector rectangle has no movable core and clean halo\n");
        goto fail;
    }
    solved_uv = (double *)malloc((size_t)crop.H * (size_t)crop.W * 2u *
                                 sizeof *solved_uv);
    solve_arena = Arena_new();
    if (solved_uv == NULL || solve_arena == NULL) goto fail;

    fprintf(stderr,
        "[marble_uv_repair] StrokeStrip variational pass\n"
        "  input=%s grid=%dx%d bbox=[%.6g %.6g]x[%.6g %.6g]\n"
        "  raw contract=%s; raster=%dx%d du/dv=%.6g/%.6g\n"
        "  rect=[%.3f %.3f %.3f %.3f] crop=r%d..%d,c%d..%d (%dx%d) "
        "core=%zu movable=%zu\n"
        "  coarse stride=%dx%d relaxed cap=%d tolerance=%.3g; V/XYZ/topology frozen\n",
        resolved, grid.H, grid.W, grid.umin, grid.umax, grid.vmin, grid.vmax,
        coverage_path, texture_w, texture_h, raster_du, raster_dv,
        rect[0], rect[1], rect[2], rect[3], crop.r0, crop.r1, crop.c0, crop.c1,
        crop.H, crop.W, crop.core_vertices, crop.movable_vertices,
        options.sample_stride_u, options.sample_stride_v,
        options.max_irls_iterations, options.movement_tolerance);
    started = ves_clock_sec();
    if (MarbleStripUv_solve(solve_arena, crop.verts, crop.reference_uv,
            crop.trust, crop.pin, crop.H, crop.W, &options,
            solved_uv, &stats) != 0) {
        fprintf(stderr, "ERROR: StrokeStrip variational solve failed\n");
        goto fail;
    }
    seconds = ves_clock_sec() - started;
    if (commit_crop(&mesh, &grid, &crop, solved_uv) != 0) {
        fprintf(stderr, "ERROR: solver violated its monotone-U/fixed-V invariant\n");
        goto fail;
    }
    if (ves_ensure_parent_dir(output) != 0 ||
        ves_ensure_parent_dir(report_path) != 0 ||
        MeshBin_write(output, mesh.verts, mesh.nv, mesh.faces, mesh.nf,
                      mesh.uv) != 0) {
        fprintf(stderr, "ERROR: cannot checkpoint %s\n", output);
        goto fail;
    }
    {
        FILE *report = fopen(report_path, "wb");
        if (report == NULL) goto fail;
        fputs("{\n \"schema\":\"marble-uv-strokestrip-v2\",\n \"input\":", report);
        json_string(report, resolved);
        fputs(",\n \"rawtex\":", report); json_string(report, rawtex_path);
        fputs(",\n \"repair_mask\":", report); json_string(report, mask_path);
        fputs(",\n \"raw_coverage_contract\":", report); json_string(report, coverage_path);
        fputs(",\n \"output\":", report); json_string(report, output);
        fprintf(report,
            ",\n \"contract\":{\"xyz_frozen\":true,\"topology_frozen\":true,"
            "\"v_frozen\":true,\"uv_only\":true,\"post_geometry_steps\":false,"
            "\"strict_complete_raw\":true},\n"
            " \"grid\":{\"height\":%d,\"width\":%d,\"bbox\":[%.17g,%.17g,%.17g,%.17g]},\n"
            " \"raster\":{\"width\":%d,\"height\":%d,\"du\":%.17g,\"dv\":%.17g},\n"
            " \"repair_rectangle_uv\":[%.17g,%.17g,%.17g,%.17g],\n"
            " \"crop\":{\"r0\":%d,\"r1\":%d,\"c0\":%d,\"c1\":%d,"
            "\"height\":%d,\"width\":%d,\"core_vertices\":%zu,"
            "\"movable_vertices\":%zu,\"feather_uv\":%.17g,\"clean_halo_uv\":%.17g},\n"
            " \"objective\":{\"kind\":\"StrokeStrip relaxed/final variational U\","
            "\"stride_u\":%d,\"stride_v\":%d,\"l1_rounds\":%d,"
            "\"relaxed_solve_cap\":%d,\"final_solve_cap\":%d,"
            "\"final_damping\":%.17g,\"movement_tolerance\":%.17g,"
            "\"likelihood_sigma_requested\":%.17g,"
            "\"lambda_length\":%.17g,\"lambda_align\":%.17g,"
            "\"lambda_local\":%.17g,\"lambda_prior\":%.17g,"
            "\"monotone_fraction\":%.17g,\"membership_floor\":%.17g,"
            "\"match_radius\":%.17g,\"match_angle_deg\":%.17g,"
            "\"match_search_columns\":%d,\"topology_fallback\":%.17g,"
            "\"length_winsor\":%.17g,\"max_displacement\":%.17g},\n"
            " \"solve\":{\"coarse_height\":%d,\"coarse_width\":%d,"
            "\"coarse_samples\":%zu,\"coarse_cross_sections\":%zu,"
            "\"coarse_members\":%zu,\"exact_anchors\":%zu,"
            "\"relaxed_qp_solves\":%d,\"final_qp_solves\":%d,"
            "\"converged\":%s,\"last_max_u_change\":%.17g,"
            "\"likelihood_sigma\":%.17g,\"downweighted_members\":%zu,"
            "\"membership_min_mean_max\":[%.17g,%.17g,%.17g],"
            "\"length_rms_before_after\":[%.17g,%.17g],"
            "\"align_rms_before_after\":[%.17g,%.17g],"
            "\"local_speed_rms_before_after\":[%.17g,%.17g],"
            "\"min_monotone_ratio\":%.17g,\"safe_commit_alpha\":%.17g,"
            "\"moved_vertices\":%zu,\"rms_u_displacement\":%.17g,"
            "\"max_u_displacement\":%.17g,\"seconds\":%.17g},\n"
            " \"orientation_invariant\":\"strict rowwise U and fixed ordered V\"\n}\n",
            grid.H, grid.W, grid.umin, grid.umax, grid.vmin, grid.vmax,
            texture_w, texture_h, raster_du, raster_dv,
            rect[0], rect[1], rect[2], rect[3],
            crop.r0, crop.r1, crop.c0, crop.c1, crop.H, crop.W,
            crop.core_vertices, crop.movable_vertices, feather, halo,
            options.sample_stride_u, options.sample_stride_v,
            options.l1_iterations, options.max_irls_iterations,
            options.final_iterations, options.final_damping,
            options.movement_tolerance,
            options.likelihood_sigma, options.lambda_length,
            options.lambda_align, options.lambda_local, options.lambda_prior,
            options.monotone_fraction, options.membership_floor,
            options.match_radius, options.match_angle_deg,
            options.match_search_columns, options.topology_fallback,
            options.length_winsor, options.max_displacement,
            stats.coarse_h, stats.coarse_w, stats.coarse_samples,
            stats.coarse_cross_sections, stats.coarse_members,
            stats.exact_anchors, stats.robust_solves, stats.final_solves,
            stats.converged ? "true" : "false", stats.last_max_u_change,
            stats.likelihood_sigma, stats.downweighted_members,
            stats.membership_min, stats.membership_mean, stats.membership_max,
            stats.length_rms_before, stats.length_rms_after,
            stats.align_rms_before, stats.align_rms_after,
            stats.local_speed_rms_before, stats.local_speed_rms_after,
            stats.min_monotone_ratio, stats.safe_commit_alpha,
            stats.moved_vertices, stats.rms_u_displacement,
            stats.max_u_displacement, seconds);
        if (fclose(report) != 0) goto fail;
    }
    fprintf(stderr,
        "[marble_uv_repair] checkpointed %s\n"
        "  relaxed/final solves=%d/%d converged=%s last-dU=%.6g; "
        "moved=%zu rms/max=%.6g/%.6g alpha=%.6g (%.1fs)\n",
        output, stats.robust_solves, stats.final_solves,
        stats.converged ? "yes" : "no", stats.last_max_u_change,
        stats.moved_vertices, stats.rms_u_displacement,
        stats.max_u_displacement, stats.safe_commit_alpha, seconds);
    free(solved_uv); crop_dispose(&crop); mask_dispose(&mask);
    MeshBin_dispose(&mesh); Arena_dispose(&solve_arena); Arena_dispose(&image_arena);
    return 0;
fail:
    free(solved_uv); crop_dispose(&crop); mask_dispose(&mask);
    MeshBin_dispose(&mesh);
    if (solve_arena != NULL) Arena_dispose(&solve_arena);
    if (image_arena != NULL) Arena_dispose(&image_arena);
    return 1;
}
