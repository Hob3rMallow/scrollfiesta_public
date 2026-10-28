/*
 * closeup_strip_turn -- a close-up's canvas turned rigidly onto its strip piece's frame (2026-10-21).
 *
 *   closeup_strip_turn <job.json>
 *   closeup_strip_turn --selftest
 *
 * The user, 2026-10-21, for the fibre-framed sheets (PHerc0358 bk2ff: each piece already turned so its fibre cross lies
 * along the axes): "Do it" -- turn each close-up to match its strip piece instead of aligning it to its own fibre
 * estimate.  That estimate asked for non-rigid warps of p50 50-430 px, which the reframe gate refused, and one close-up
 * (cu20261021T221752Z_839632) came out a quarter turn off its strip, which the stage cannot paste.
 *
 * job.json:
 *   "seg"      the close-up's tifxyz export (mask.tif: the canvas and its valid nodes; meta.json, for --reuv)
 *   "matches"  a float64 file of (col, row, x, y) per match: a close-up canvas node and the strip crop px of the same
 *              papyrus, matched in 3-D within MATCH_VOX (closeup_splice.strip_matches: the stage's own matching)
 *   "count"    the number of matches in that file
 *   "out"      the output directory
 * Tuning: pipeline_constants.h STRIP_TURN_*.
 *
 * The turn is the 2-D Kabsch rotation carrying the canvas onto the strip, refitted STRIP_TURN_TRIM_ROUNDS times over
 * the best STRIP_TURN_TRIM_KEEP of the residuals.  When a reflection relates the two better than any rotation (the
 * strip shows the piece reflected for reading: a mirror-handed piece reads flip_lr / flip_ud), the fit is that
 * reflection, R(t) diag(1, -1), and the written turn is the rotation that leaves only a reflection about one canvas
 * axis -- R(t) diag(1, -1) = diag(-1, 1) R(180 - t) = diag(1, -1) R(-t), the smaller turn -- which the stage pastes
 * (closeup_splice.canvas_turn; the user, 2026-10-23: "it's unclear to me that mirroring *should* cause rejection").
 * Refused (exit 3; strip_turn.json says why and no map is left behind): fewer than STRIP_TURN_MIN_MATCHES matches, or
 * an rms residual over STRIP_TURN_MAX_RESIDUAL x the matches' spread.
 *
 * Writes the fiber_flatten map contract that ink_closeup --reuv reads: <out>/fiber_flatten.json (seg, grid, new_page,
 * map.affine) and <out>/fiber_map_u.tif / fiber_map_v.tif, the rotation at every node of a grid reaching
 * STRIP_TURN_PAD_CELLS past the canvas.  Every vertex is then interpolated bilinearly, and the bilinear interpolation of
 * a rotation is the rotation: each triangle keeps its 2-D metric and orientation.  <out>/strip_turn.json: the fit.
 * Exit 0 turned, 1 error, 2 usage, 3 refused.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/pipeline_constants.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"

#define ST_PI 3.14159265358979323846
#define ST_REFUSED 3

typedef struct StripTurnFit {
    double theta;              /* the fit's turn (radians) carrying the canvas (col, row) -- reflected first when
                                * reflect -- onto the strip (x, y) */
    double cos_t, sin_t;
    int reflect;               /* 1: the fit is R(theta) diag(1, -1) (a mirrored piece) */
    double turn;               /* the written turn (radians): theta, or for a reflection the smaller rotation that
                                * leaves only the reflection about one canvas axis named by mirror */
    const char *mirror;        /* NULL, "flip_lr" (columns reversed) or "flip_ud" (rows reversed) after the turn */
    double centre_a[2];        /* the kept matches' centre on the canvas ... */
    double centre_b[2];        /* ... and on the strip */
    double rms_px;             /* the rigid fit's rms residual over the kept matches */
    double spread_px;          /* the kept matches' rms distance from their canvas centre */
    double score_rotation;     /* |sum b . R a| of the best rotation and of the best reflection, over all matches */
    double score_reflection;
    size_t matches, used;
    const char *refused;       /* NULL, or why no turn is written */
} StripTurnFit;

/* a JSON string literal: quotes, backslashes (Windows paths) and control characters escaped */
static void st_json_str(FILE *f, const char *s)
{
    fputc('"', f);
    for (; s != NULL && *s != '\0'; s++) {
        const unsigned char ch = (unsigned char)*s;
        if (ch == '"' || ch == '\\') { fputc('\\', f); fputc((int)ch, f); }
        else if (ch < 0x20) fprintf(f, "\\u%04x", (unsigned)ch);
        else fputc((int)ch, f);
    }
    fputc('"', f);
}

static int st_cmp_double(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* The best rotation of the kept canvas points -- reflected first (ay -> -ay) when fit->reflect -- onto their strip
 * points (2-D Kabsch about both centres), and the scores of the best rotation and the best reflection of the
 * unreflected points over the same matches. */
static void st_kabsch(const double *m, size_t n, const uint8_t *keep, StripTurnFit *fit)
{
    double sa[2] = { 0.0, 0.0 }, sb[2] = { 0.0, 0.0 };
    double h00 = 0.0, h01 = 0.0, h10 = 0.0, h11 = 0.0;
    size_t i = 0, k = 0;
    for (i = 0; i < n; i++) {
        if (!keep[i]) continue;
        sa[0] += m[4 * i]; sa[1] += m[4 * i + 1];
        sb[0] += m[4 * i + 2]; sb[1] += m[4 * i + 3];
        k++;
    }
    if (k == 0) return;
    fit->centre_a[0] = sa[0] / (double)k; fit->centre_a[1] = sa[1] / (double)k;
    fit->centre_b[0] = sb[0] / (double)k; fit->centre_b[1] = sb[1] / (double)k;
    for (i = 0; i < n; i++) {
        double ax = 0.0, ay = 0.0, bx = 0.0, by = 0.0;
        if (!keep[i]) continue;
        ax = m[4 * i] - fit->centre_a[0]; ay = m[4 * i + 1] - fit->centre_a[1];
        bx = m[4 * i + 2] - fit->centre_b[0]; by = m[4 * i + 3] - fit->centre_b[1];
        h00 += ax * bx; h01 += ax * by; h10 += ay * bx; h11 += ay * by;
    }
    /* sum b . R(t) a = cos t (h00 + h11) + sin t (h01 - h10); a reflection first (ay -> -ay) gives
     * cos t (h00 - h11) + sin t (h01 + h10) */
    fit->theta = fit->reflect ? atan2(h01 + h10, h00 - h11) : atan2(h01 - h10, h00 + h11);
    fit->cos_t = cos(fit->theta);
    fit->sin_t = sin(fit->theta);
    fit->score_rotation = hypot(h00 + h11, h01 - h10);
    fit->score_reflection = hypot(h00 - h11, h01 + h10);
    fit->used = k;
}

/* res[i] = |b_i - R a_i| about the current centres (a_i reflected first when fit->reflect), every match */
static void st_residuals(const double *m, size_t n, const StripTurnFit *fit, double *res)
{
    const double sy = fit->reflect ? -1.0 : 1.0;
    size_t i = 0;
    for (i = 0; i < n; i++) {
        const double ax = m[4 * i] - fit->centre_a[0], ay = sy * (m[4 * i + 1] - fit->centre_a[1]);
        const double bx = m[4 * i + 2] - fit->centre_b[0], by = m[4 * i + 3] - fit->centre_b[1];
        const double ex = bx - (fit->cos_t * ax - fit->sin_t * ay), ey = by - (fit->sin_t * ax + fit->cos_t * ay);
        res[i] = hypot(ex, ey);
    }
}

static double st_wrap(double a)                       /* an angle in (-pi, pi] */
{
    while (a > ST_PI) a -= 2.0 * ST_PI;
    while (a <= -ST_PI) a += 2.0 * ST_PI;
    return a;
}

/* The written turn: the fit's rotation, or for a reflection R(t) diag(1, -1) the smaller of the two rotations that
 * leave a reflection about one canvas axis: diag(1, -1) R(-t) (rows reversed, flip_ud) or diag(-1, 1) R(180 - t)
 * (columns reversed, flip_lr). */
static void st_turn(StripTurnFit *fit)
{
    double ud = 0.0, lr = 0.0;
    fit->mirror = NULL;
    fit->turn = fit->theta;
    if (!fit->reflect) return;
    ud = st_wrap(-fit->theta);
    lr = st_wrap(ST_PI - fit->theta);
    if (fabs(lr) < fabs(ud)) { fit->turn = lr; fit->mirror = "flip_lr"; }
    else { fit->turn = ud; fit->mirror = "flip_ud"; }
}

/* The trimmed rigid fit (see the file comment).  0 with fit->refused NULL (a turn) or set (refused). */
static int StripTurn_fit(Arena_T arena, const double *m, size_t n, StripTurnFit *fit)
{
    uint8_t *keep = NULL;
    double *res = NULL, *sorted = NULL, sq = 0.0, sp = 0.0;
    size_t i = 0;
    int round = 0;
    memset(fit, 0, sizeof *fit);
    fit->matches = n;
    fit->cos_t = 1.0;
    if (n < (size_t)STRIP_TURN_MIN_MATCHES) {
        fit->refused = "too few matched samples between the strip box and the close-up";
        return 0;
    }
    keep = (uint8_t *)ARENA_ALLOC(arena, n);
    res = (double *)ARENA_ALLOC(arena, n * sizeof *res);
    sorted = (double *)ARENA_ALLOC(arena, n * sizeof *sorted);
    memset(keep, 1, n);
    st_kabsch(m, n, keep, fit);
    if (fit->score_reflection > fit->score_rotation) {   /* a mirrored piece: fit the reflection */
        const double rot = fit->score_rotation, ref = fit->score_reflection;
        fit->reflect = 1;
        st_kabsch(m, n, keep, fit);
        fit->score_rotation = rot;
        fit->score_reflection = ref;
    }
    for (round = 0; round < STRIP_TURN_TRIM_ROUNDS; round++) {
        const double rot = fit->score_rotation, ref = fit->score_reflection;
        size_t cut = 0;
        st_residuals(m, n, fit, res);
        memcpy(sorted, res, n * sizeof *res);
        qsort(sorted, n, sizeof *sorted, st_cmp_double);
        cut = (size_t)floor(STRIP_TURN_TRIM_KEEP * (double)(n - 1));
        for (i = 0; i < n; i++) keep[i] = (uint8_t)(res[i] <= sorted[cut]);
        st_kabsch(m, n, keep, fit);
        fit->score_rotation = rot;          /* the mirror test's scores stay those of every match */
        fit->score_reflection = ref;
    }
    st_residuals(m, n, fit, res);
    for (i = 0; i < n; i++) {
        double ax = 0.0, ay = 0.0;
        if (!keep[i]) continue;
        ax = m[4 * i] - fit->centre_a[0]; ay = m[4 * i + 1] - fit->centre_a[1];
        sq += res[i] * res[i];
        sp += ax * ax + ay * ay;
    }
    fit->rms_px = fit->used ? sqrt(sq / (double)fit->used) : 0.0;
    fit->spread_px = fit->used ? sqrt(sp / (double)fit->used) : 0.0;
    st_turn(fit);
    if (!(fit->spread_px > 0.0) || fit->rms_px > STRIP_TURN_MAX_RESIDUAL * fit->spread_px)
        fit->refused = "no single turn relates the close-up to its strip (rigid fit residual too large)";
    return 0;
}

static int st_read_matches(Arena_T arena, const char *path, size_t n, double **out)
{
    FILE *fp = fopen(path, "rb");
    double *m = NULL;
    int extra = 0;
    if (fp == NULL) { fprintf(stderr, "ERROR: cannot open %s\n", path); return -1; }
    m = (double *)ARENA_ALLOC(arena, (n ? n : 1) * 4 * sizeof *m);
    if (n && fread(m, sizeof *m, n * 4, fp) != n * 4) {
        fprintf(stderr, "ERROR: %s holds fewer than %zu matches\n", path, n);
        fclose(fp);
        return -1;
    }
    extra = fgetc(fp) != EOF;
    fclose(fp);
    if (extra) { fprintf(stderr, "ERROR: %s holds more than %zu matches\n", path, n); return -1; }
    *out = m;
    return 0;
}

static void st_unlink_map(const char *out)
{
    static const char *const names[3] = { "fiber_flatten.json", "fiber_map_u.tif", "fiber_map_v.tif" };
    char path[2048];
    int k = 0;
    for (k = 0; k < 3; k++) {
        snprintf(path, sizeof path, "%s/%s", out, names[k]);
        (void)ves_unlink(path);
    }
}

static int st_write_report(const char *out, const char *seg, const char *matches, const StripTurnFit *fit, int W,
                           int H, const double page[4])
{
    char path[2048];
    FILE *fp = NULL;
    snprintf(path, sizeof path, "%s/strip_turn.json", out);
    fp = fopen(path, "wb");
    if (fp == NULL) { fprintf(stderr, "ERROR: cannot write %s\n", path); return -1; }
    fprintf(fp, "{\n  \"schema\": \"closeup-strip-turn-v1\",\n  \"seg\": ");
    st_json_str(fp, seg);
    fprintf(fp, ",\n  \"matches\": ");
    st_json_str(fp, matches);
    fprintf(fp, ",\n  \"count\": %zu,\n  \"used\": %zu,\n  \"canvas\": [%d, %d],\n", fit->matches, fit->used, W, H);
    fprintf(fp, "  \"turn_deg\": %.6f,\n  \"fit_deg\": %.6f,\n  \"reflection\": %s,\n  \"mirror\": ",
            fit->turn * 180.0 / ST_PI, fit->theta * 180.0 / ST_PI, fit->reflect ? "true" : "false");
    if (fit->mirror != NULL) st_json_str(fp, fit->mirror);
    else fputs("null", fp);
    fprintf(fp, ",\n  \"rms_px\": %.4f,\n  \"spread_px\": %.4f,\n", fit->rms_px, fit->spread_px);
    fprintf(fp, "  \"score_rotation\": %.6g,\n  \"score_reflection\": %.6g,\n", fit->score_rotation,
            fit->score_reflection);
    if (page != NULL)
        fprintf(fp, "  \"new_page\": {\"u_min\": %.6f, \"v_min\": %.6f, \"width\": %.0f, \"height\": %.0f},\n", page[0],
                page[1], page[2], page[3]);
    fprintf(fp, "  \"constants\": {\"min_matches\": %d, \"trim_rounds\": %d, \"trim_keep\": %.3f, \"max_residual\": %.3f, "
                "\"cell_px\": %d, \"pad_cells\": %d},\n", STRIP_TURN_MIN_MATCHES, STRIP_TURN_TRIM_ROUNDS,
            STRIP_TURN_TRIM_KEEP, STRIP_TURN_MAX_RESIDUAL, STRIP_TURN_CELL_PX, STRIP_TURN_PAD_CELLS);
    fprintf(fp, "  \"refused\": ");
    if (fit->refused != NULL) st_json_str(fp, fit->refused);
    else fputs("null", fp);
    fprintf(fp, "\n}\n");
    return fclose(fp) == 0 ? 0 : -1;
}

/* The rotation at every node of a grid STRIP_TURN_PAD_CELLS past the W x H canvas, and the new page's origin over the
 * valid nodes: page = {u_min, v_min, width, height}. */
static int st_write_map(Arena_T arena, const char *out, const char *seg, const StripTurnFit *fit, const uint8_t *mask,
                        int W, int H, double page[4])
{
    const double cell = (double)STRIP_TURN_CELL_PX, o = -(double)STRIP_TURN_PAD_CELLS * cell;
    const int rows = (int)ceil((double)H / cell) + 1 + 2 * STRIP_TURN_PAD_CELLS;
    const int cols = (int)ceil((double)W / cell) + 1 + 2 * STRIP_TURN_PAD_CELLS;
    const double c = cos(fit->turn), s = sin(fit->turn);    /* the written turn, a rotation even for a mirror */
    double umin = 1e300, vmin = 1e300, umax = -1e300, vmax = -1e300;
    float *mu = NULL, *mv = NULL;
    char path[2048];
    FILE *fp = NULL;
    size_t q = 0, valid = 0;
    int i = 0, j = 0;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            double u = 0.0, v = 0.0;
            if (!mask[(size_t)i * (size_t)W + (size_t)j]) continue;
            u = c * (double)j - s * (double)i;
            v = s * (double)j + c * (double)i;
            if (u < umin) umin = u;
            if (u > umax) umax = u;
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
            valid++;
        }
    if (valid == 0) { fprintf(stderr, "ERROR: %s/mask.tif has no valid node\n", seg); return -1; }
    page[0] = floor(umin) - 1.0;
    page[1] = floor(vmin) - 1.0;
    page[2] = ceil(umax - page[0]) + 2.0;
    page[3] = ceil(vmax - page[1]) + 2.0;
    mu = (float *)ARENA_ALLOC(arena, (size_t)rows * (size_t)cols * sizeof *mu);
    mv = (float *)ARENA_ALLOC(arena, (size_t)rows * (size_t)cols * sizeof *mv);
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++) {
            const double x = o + cell * (double)j, y = o + cell * (double)i;
            q = (size_t)i * (size_t)cols + (size_t)j;
            mu[q] = (float)(c * x - s * y);
            mv[q] = (float)(s * x + c * y);
        }
    snprintf(path, sizeof path, "%s/fiber_map_u.tif", out);
    if (TiffIO_save_float2d(path, mu, cols, rows) != 0) { fprintf(stderr, "ERROR: cannot write %s\n", path); return -1; }
    snprintf(path, sizeof path, "%s/fiber_map_v.tif", out);
    if (TiffIO_save_float2d(path, mv, cols, rows) != 0) { fprintf(stderr, "ERROR: cannot write %s\n", path); return -1; }
    snprintf(path, sizeof path, "%s/fiber_flatten.json", out);
    fp = fopen(path, "wb");
    if (fp == NULL) { fprintf(stderr, "ERROR: cannot write %s\n", path); return -1; }
    fprintf(fp, "{\n  \"schema\": \"closeup-strip-turn-v1\",\n  \"seg\": ");
    st_json_str(fp, seg);
    fprintf(fp, ",\n  \"grid\": {\"rows\": %d, \"cols\": %d, \"cell\": %.1f, \"row0\": %.1f, \"col0\": %.1f},\n", rows,
            cols, cell, o, o);
    fprintf(fp, "  \"new_page\": {\"u_min\": %.6f, \"v_min\": %.6f, \"width\": %.0f, \"height\": %.0f},\n", page[0],
            page[1], page[2], page[3]);
    fprintf(fp, "  \"map\": {\"affine\": [[%.17g, %.17g, 0.0], [%.17g, %.17g, 0.0]]},\n", c, -s, s, c);
    fprintf(fp, "  \"turn_deg\": %.6f,\n  \"mirror\": ", fit->turn * 180.0 / ST_PI);
    if (fit->mirror != NULL) st_json_str(fp, fit->mirror);
    else fputs("null", fp);
    fprintf(fp, ",\n  \"strip_turn\": \"strip_turn.json\"\n}\n");
    return fclose(fp) == 0 ? 0 : -1;
}

/* the job (see the file comment): 0 turned, 1 error, ST_REFUSED refused */
static int st_run(const char *job)
{
    Arena_T a = Arena_new();
    const char *err = NULL, *seg = NULL, *matches = NULL, *out = NULL;
    const JsonValue *j = NULL;
    double *m = NULL, page[4] = { 0.0, 0.0, 0.0, 0.0 };
    uint8_t *mask = NULL;
    int D = 0, H = 0, W = 0, rc = 1;
    long count = 0;
    StripTurnFit fit;
    char path[2048];
    memset(&fit, 0, sizeof fit);
    if (a == NULL) return 1;
    j = Json_parse_file(a, job, &err);
    if (j == NULL) { fprintf(stderr, "ERROR: cannot parse %s (%s)\n", job, err ? err : "?"); goto out; }
    seg = Json_as_string(Json_object_get(j, "seg"));
    matches = Json_as_string(Json_object_get(j, "matches"));
    out = Json_as_string(Json_object_get(j, "out"));
    count = Json_member_long(j, "count", -1);
    if (seg == NULL || matches == NULL || out == NULL || count < 0) {
        fprintf(stderr, "ERROR: %s lacks seg/matches/count/out\n", job);
        goto out;
    }
    printf("[orch] stage load\n");
    fflush(stdout);
    snprintf(path, sizeof path, "%s/mask.tif", seg);
    if (TiffIO_load(a, path, &mask, &D, &H, &W) != 0 || D != 1 || W <= 0 || H <= 0) {
        fprintf(stderr, "ERROR: cannot read the canvas mask %s\n", path);
        goto out;
    }
    if (st_read_matches(a, matches, (size_t)count, &m) != 0) goto out;
    printf("[orch] stage fit\n");
    fflush(stdout);
    if (StripTurn_fit(a, m, (size_t)count, &fit) != 0) goto out;
    snprintf(path, sizeof path, "%s/x", out);
    if (ves_ensure_parent_dir(path) != 0) { fprintf(stderr, "ERROR: cannot create %s\n", out); goto out; }
    st_unlink_map(out);                              /* a refusal never leaves an older turn behind */
    if (fit.refused != NULL) {
        if (st_write_report(out, seg, matches, &fit, W, H, NULL) != 0) goto out;
        fprintf(stderr, "REFUSED: %s (%zu matches; %s %.2f deg, rms %.1f px over a %.1f px spread)\n", fit.refused,
                fit.matches, fit.reflect ? "reflection" : "turn", fit.theta * 180.0 / ST_PI, fit.rms_px, fit.spread_px);
        rc = ST_REFUSED;
        goto out;
    }
    printf("[orch] stage write\n");
    fflush(stdout);
    if (st_write_map(a, out, seg, &fit, mask, W, H, page) != 0) { st_unlink_map(out); goto out; }
    if (st_write_report(out, seg, matches, &fit, W, H, page) != 0) { st_unlink_map(out); goto out; }
    printf("[closeup_strip_turn] turn %.3f deg%s%s from %zu of %zu matches: rms %.2f px over a %.1f px spread; new "
           "page %.0f x %.0f\n", fit.turn * 180.0 / ST_PI, fit.mirror ? ", then mirrored for the stage: " : "",
           fit.mirror ? fit.mirror : "", fit.used, fit.matches, fit.rms_px, fit.spread_px, page[2], page[3]);
    printf("[orch] stage done\n");
    fflush(stdout);
    rc = 0;
out:
    Arena_dispose(&a);
    return rc;
}

/* ---------------------------------------------------------------------------------------------------------- selftest */

static uint64_t st_rand_state = 0x9e3779b97f4a7c15ULL;

static double st_rand(void)                          /* uniform [0, 1), a fixed sequence */
{
    st_rand_state = st_rand_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(st_rand_state >> 11) / 9007199254740992.0;
}

/* A selftest case: n matches of valid canvas nodes carried by a turn (mirror: after a reflection) plus outliers;
 * returns the job's exit code. */
static int st_case(const char *dir, const char *seg, const uint8_t *mask, int W, int H, size_t n, double turn_deg,
                   int mirror, int unrelated, const char *name)
{
    char job[1024], mpath[1024], out[1024];
    const double t = turn_deg * ST_PI / 180.0, c = cos(t), s = sin(t);
    FILE *fp = NULL;
    size_t k = 0;
    snprintf(mpath, sizeof mpath, "%s/%s_matches.f64", dir, name);
    snprintf(out, sizeof out, "%s/%s_out", dir, name);
    snprintf(job, sizeof job, "%s/%s_job.json", dir, name);
    fp = fopen(mpath, "wb");
    if (fp == NULL) return -1;
    while (k < n) {
        const int col = (int)(st_rand() * (double)W), row = (int)(st_rand() * (double)H);
        double rec[4] = { 0.0, 0.0, 0.0, 0.0 }, y = 0.0;
        if (!mask[(size_t)row * (size_t)W + (size_t)col]) continue;
        y = mirror ? -(double)row : (double)row;
        rec[0] = (double)col; rec[1] = (double)row;
        rec[2] = c * (double)col - s * y + 500.0 + (st_rand() - 0.5);
        rec[3] = s * (double)col + c * y - 40.0 + (st_rand() - 0.5);
        if (unrelated || st_rand() < 0.15) {        /* 15% wrong-wrap matches anywhere on the strip */
            rec[2] = st_rand() * 2000.0;
            rec[3] = st_rand() * 2000.0;
        }
        if (fwrite(rec, sizeof rec[0], 4, fp) != 4) { fclose(fp); return -1; }
        k++;
    }
    fclose(fp);
    fp = fopen(job, "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "{\"seg\": ");
    st_json_str(fp, seg);
    fprintf(fp, ", \"matches\": ");
    st_json_str(fp, mpath);
    fprintf(fp, ", \"count\": %zu, \"out\": ", n);
    st_json_str(fp, out);
    fprintf(fp, "}\n");
    fclose(fp);
    return st_run(job);
}

/* A written map (<dir>/<name>_out): the fiber_flatten.json contract, its turn want_deg and mirror, and ink_closeup
 * --reuv's bilinear lookup at every valid node landing on the turned position, on the new page.  Failures. */
static int st_check_map(Arena_T a, const char *dir, const char *name, const char *seg, const uint8_t *mask, int W,
                        int H, double want_deg, const char *want_mirror)
{
    char path[1024];
    const char *err = NULL, *mirror = NULL;
    const JsonValue *fj = NULL, *gr = NULL, *np = NULL, *af = NULL;
    float *mu = NULL, *mv = NULL;
    int mw = 0, mh = 0, rows = 0, cols = 0, bad = 0, fails = 0, i = 0, j = 0;
    double cell = 0.0, o = 0.0, umin = 0.0, vmin = 0.0, got = 0.0;
    snprintf(path, sizeof path, "%s/%s_out/fiber_flatten.json", dir, name);
    fj = Json_parse_file(a, path, &err);
    gr = fj ? Json_object_get(fj, "grid") : NULL;
    np = fj ? Json_object_get(fj, "new_page") : NULL;
    af = fj ? Json_object_get(Json_object_get(fj, "map"), "affine") : NULL;
    if (gr == NULL || np == NULL || af == NULL || Json_as_string(Json_object_get(fj, "seg")) == NULL ||
        strcmp(Json_as_string(Json_object_get(fj, "seg")), seg) != 0) {
        fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s fiber_flatten.json contract\n", name);
        return 1;
    }
    rows = (int)Json_member_long(gr, "rows", 0); cols = (int)Json_member_long(gr, "cols", 0);
    cell = Json_member_double(gr, "cell", 0.0); o = Json_member_double(gr, "row0", 0.0);
    umin = Json_member_double(np, "u_min", 0.0); vmin = Json_member_double(np, "v_min", 0.0);
    got = Json_member_double(fj, "turn_deg", 0.0);
    mirror = Json_as_string(Json_object_get(fj, "mirror"));
    if (fabs(got - want_deg) > 0.05) {
        fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s turn %.4f, want %.1f\n", name, got, want_deg);
        fails++;
    }
    if ((mirror == NULL) != (want_mirror == NULL) || (mirror != NULL && strcmp(mirror, want_mirror) != 0)) {
        fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s mirror %s, want %s\n", name, mirror ? mirror : "null",
                want_mirror ? want_mirror : "null");
        fails++;
    }
    if (fabs(Json_as_double(Json_array_get(Json_array_get(af, 1), 0), 0.0) - sin(want_deg * ST_PI / 180.0)) > 1e-3) {
        fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s map.affine\n", name);
        fails++;
    }
    snprintf(path, sizeof path, "%s/%s_out/fiber_map_u.tif", dir, name);
    if (TiffIO_load_float2d(a, path, &mu, &mw, &mh) != 0 || mw != cols || mh != rows) bad++;
    snprintf(path, sizeof path, "%s/%s_out/fiber_map_v.tif", dir, name);
    if (TiffIO_load_float2d(a, path, &mv, &mw, &mh) != 0 || mw != cols || mh != rows) bad++;
    if (bad) {
        fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s map files\n", name);
        return fails + 1;
    }
    {
        const double c = cos(got * ST_PI / 180.0), s = sin(got * ST_PI / 180.0);
        double worst = 0.0, low = 1e300;
        for (i = 0; i < H; i++)
            for (j = 0; j < W; j++) {
                const double gi = ((double)i - o) / cell, gj = ((double)j - o) / cell;
                const int i0 = (int)floor(gi), j0 = (int)floor(gj);
                const double fy = gi - (double)i0, fx = gj - (double)j0;
                double u = 0.0, v = 0.0;
                size_t q00 = 0, q01 = 0, q10 = 0, q11 = 0;
                if (!mask[(size_t)i * (size_t)W + (size_t)j]) continue;
                if (i0 < 0 || j0 < 0 || i0 + 1 >= rows || j0 + 1 >= cols) { worst = 1e9; continue; }
                q00 = (size_t)i0 * (size_t)cols + (size_t)j0; q01 = q00 + 1;
                q10 = q00 + (size_t)cols; q11 = q10 + 1;
                u = (1 - fy) * ((1 - fx) * mu[q00] + fx * mu[q01]) + fy * ((1 - fx) * mu[q10] + fx * mu[q11]);
                v = (1 - fy) * ((1 - fx) * mv[q00] + fx * mv[q01]) + fy * ((1 - fx) * mv[q10] + fx * mv[q11]);
                worst = fmax(worst, hypot(u - (c * j - s * i), v - (s * j + c * i)));
                low = fmin(low, fmin(u - umin, v - vmin));
            }
        if (worst > 2e-3 || low < 0.0) {
            fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s reuv lookup off by %.4g px, new page low %.3f\n",
                    name, worst, low);
            fails++;
        }
    }
    return fails;
}

/* The written turn followed by the named reflection carries the canvas onto the strip by a translation: the matches
 * of <dir>/<name>_matches.f64 within 2 px of their median offset, at least 80% (15% are outliers).  Failures. */
static int st_check_flip(const char *dir, const char *name, size_t n, double turn_deg, const char *mirror)
{
    char path[1024];
    const double c = cos(turn_deg * ST_PI / 180.0), s = sin(turn_deg * ST_PI / 180.0);
    double *m = NULL, *dx = NULL, *dy = NULL, mx = 0.0, my = 0.0;
    size_t i = 0, kept = 0;
    FILE *fp = NULL;
    int fails = 0;
    snprintf(path, sizeof path, "%s/%s_matches.f64", dir, name);
    m = (double *)malloc(n * 4 * sizeof *m);
    dx = (double *)malloc(n * sizeof *dx);
    dy = (double *)malloc(n * sizeof *dy);
    fp = fopen(path, "rb");
    if (m == NULL || dx == NULL || dy == NULL || fp == NULL || fread(m, sizeof *m, n * 4, fp) != n * 4) {
        fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s matches unreadable\n", name);
        fails++;
    } else {
        for (i = 0; i < n; i++) {
            double px = c * m[4 * i] - s * m[4 * i + 1], py = s * m[4 * i] + c * m[4 * i + 1];
            if (strcmp(mirror, "flip_lr") == 0) px = -px;
            else py = -py;
            dx[i] = m[4 * i + 2] - px;
            dy[i] = m[4 * i + 3] - py;
        }
        {
            double *sx = (double *)malloc(n * sizeof *sx), *sy = (double *)malloc(n * sizeof *sy);
            if (sx == NULL || sy == NULL) fails++;
            else {
                memcpy(sx, dx, n * sizeof *sx);
                memcpy(sy, dy, n * sizeof *sy);
                qsort(sx, n, sizeof *sx, st_cmp_double);
                qsort(sy, n, sizeof *sy, st_cmp_double);
                mx = sx[n / 2];
                my = sy[n / 2];
                for (i = 0; i < n; i++) kept += hypot(dx[i] - mx, dy[i] - my) <= 2.0;
                if ((double)kept < 0.8 * (double)n) {
                    fprintf(stderr, "[selftest] closeup_strip_turn FAIL: %s turned and %s, only %zu of %zu matches are a "
                            "translation of the strip\n", name, mirror, kept, n);
                    fails++;
                }
            }
            free(sx);
            free(sy);
        }
    }
    if (fp != NULL) fclose(fp);
    free(m);
    free(dx);
    free(dy);
    return fails;
}

static int st_selftest(void)
{
    char dir[512], seg[700], path[1024];
    const int W = 320, H = 240;
    const double turn = 92.8;
    Arena_T a = Arena_new();
    uint8_t *mask = NULL;
    int fails = 0, rc = 0, i = 0, j = 0;
    if (a == NULL || ves_selftest_dir(dir, sizeof dir, "closeup_strip_turn") != 0) return 1;
    fails += Json_read_selftest() != 0;
    snprintf(seg, sizeof seg, "%s/seg", dir);
    (void)ves_mkdir(seg);
    mask = (uint8_t *)ARENA_ALLOC(a, (size_t)W * (size_t)H);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            const double ex = ((double)j - 160.0) / 140.0, ey = ((double)i - 120.0) / 100.0;
            mask[(size_t)i * (size_t)W + (size_t)j] = (uint8_t)(ex * ex + ey * ey <= 1.0 ? 255 : 0);
        }
    snprintf(path, sizeof path, "%s/mask.tif", seg);
    if (TiffIO_save(path, mask, 1, H, W) != 0) { fprintf(stderr, "[selftest] closeup_strip_turn FAIL: mask\n"); fails++; }

    /* a quarter turn (cu...839632's 92.8 deg) with 15% outliers: the turn, an exact map, a valid new page */
    rc = st_case(dir, seg, mask, W, H, 1000, turn, 0, 0, "turn");
    if (rc != 0) { fprintf(stderr, "[selftest] closeup_strip_turn FAIL: turn case exit %d\n", rc); fails++; }
    else fails += st_check_map(a, dir, "turn", seg, mask, W, H, turn, NULL);
    /* mirrored pieces (2026-10-23): the reflection R(t) diag(1, -1) is fitted and written as the smaller turn that
     * leaves one axis reversed -- t = 30: R(-30) then rows reversed; t = 170: R(10) then columns reversed -- and that
     * turn followed by that reflection carries the canvas onto the strip by a translation */
    rc = st_case(dir, seg, mask, W, H, 1000, 30.0, 1, 0, "mirror");
    if (rc != 0) { fprintf(stderr, "[selftest] closeup_strip_turn FAIL: mirror case exit %d\n", rc); fails++; }
    else fails += st_check_map(a, dir, "mirror", seg, mask, W, H, -30.0, "flip_ud") +
                  st_check_flip(dir, "mirror", 1000, -30.0, "flip_ud");
    rc = st_case(dir, seg, mask, W, H, 1000, 170.0, 1, 0, "mirror_lr");
    if (rc != 0) { fprintf(stderr, "[selftest] closeup_strip_turn FAIL: mirror_lr case exit %d\n", rc); fails++; }
    else fails += st_check_map(a, dir, "mirror_lr", seg, mask, W, H, 10.0, "flip_lr") +
                  st_check_flip(dir, "mirror_lr", 1000, 10.0, "flip_lr");
    /* too few matches, and matches no turn relates: refused, and no map left behind */
    rc = st_case(dir, seg, mask, W, H, (size_t)STRIP_TURN_MIN_MATCHES - 1, 10.0, 0, 0, "few");
    if (rc != ST_REFUSED) { fprintf(stderr, "[selftest] closeup_strip_turn FAIL: few exit %d\n", rc); fails++; }
    rc = st_case(dir, seg, mask, W, H, 1000, 0.0, 0, 1, "unrelated");
    if (rc != ST_REFUSED) { fprintf(stderr, "[selftest] closeup_strip_turn FAIL: unrelated exit %d\n", rc); fails++; }
    /* a refusal after a turn removes the earlier turn's map */
    rc = st_case(dir, seg, mask, W, H, 1000, 0.0, 0, 1, "turn");
    snprintf(path, sizeof path, "%s/turn_out/fiber_flatten.json", dir);
    {
        FILE *probe = fopen(path, "rb");
        if (rc != ST_REFUSED || probe != NULL) {
            fprintf(stderr, "[selftest] closeup_strip_turn FAIL: a refusal left the earlier map\n");
            fails++;
        }
        if (probe != NULL) fclose(probe);
    }
    Arena_dispose(&a);
    if (fails == 0) ves_selftest_clean(dir);
    fprintf(stderr, "[selftest] closeup_strip_turn %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) return st_selftest();
    if (argc != 2) {
        fprintf(stderr, "usage: closeup_strip_turn <job.json>\n       closeup_strip_turn --selftest\n");
        return 2;
    }
    return st_run(argv[1]);
}
