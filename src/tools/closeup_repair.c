/* ============================================================================
 * closeup_repair -- a close-up surface without its duplicate and wrong-side
 * geometry (src/flatten/closeup_repair.c; 2026-10-12, built on PHerc0211
 * close-up surfaces).
 *
 *   closeup_repair <in_dir> <stem> <axis.csv> <out_dir>
 *   closeup_repair --selftest
 *
 * in_dir, stem  an ink_closeup surface pass: <stem>_surface.vmesh (UV in the
 *               close-up frame, float64), <stem>_surface_attr.u8 (1 observed,
 *               2 fill, 3 observed on an adopted chart), <stem>_surface_chart.i32,
 *               <stem>_surface_nrm.f32 -- read only
 * axis.csv      the z,y,x scroll-axis table (AsmAxis_load)
 * out_dir       created; refused when it already holds <stem>_surface.vmesh.
 *               The kept surface under the input's names; the observed flaps
 *               (never deleted) as <stem>_extras_surface.* with
 *               <stem>_extras_source_vertex.i32 (input vertex of each);
 *               <stem>_repair_reason.u8 (every input vertex's CLOSEUP_REPAIR_*);
 *               <stem>_repair_report.json (inputs with sha256, constants, counts);
 *               <stem>_repair_map.png (the frame at half scale: grey observed,
 *               blue fill kept; green duplicate, yellow wrong-side fill,
 *               magenta flap, cyan orphan).
 * ==========================================================================*/
#include "../flatten/closeup_repair.h"
#include "../assemble/asm_axis.h"
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/pipeline_constants.h"
#include "../common/sha256.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cr_read(const char *path, void *buf, size_t bytes)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    size_t got = fread(buf, 1, bytes, f);
    int extra = fgetc(f) != EOF;
    fclose(f);
    return got == bytes && !extra ? 0 : -1;
}

static int cr_write(const char *path, const void *buf, size_t bytes)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    size_t put = bytes ? fwrite(buf, 1, bytes, f) : 0;
    int rc = fclose(f);
    return put == bytes && rc == 0 ? 0 : -1;
}

static int cr_exists(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    fclose(f);
    return 1;
}

static int cr_sha256(const char *path, char hex[65], long long *bytes)
{
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    unsigned char *buf = malloc((size_t)1 << 20);
    if (!buf) { fclose(f); return -1; }
    VesSha256 s;
    VesSha256_init(&s);
    size_t n = 0;
    *bytes = 0;
    while ((n = fread(buf, 1, (size_t)1 << 20, f)) > 0) { VesSha256_update(&s, buf, n); *bytes += (long long)n; }
    fclose(f);
    free(buf);
    VesSha256_hex(&s, hex);
    return 0;
}

/* The faces whose three corners carry `want`, with their vertices compacted in input order, written as
 * <out_dir>/<name>_surface.vmesh + _surface_attr.u8 + _surface_chart.i32 + _surface_nrm.f32 (and
 * <name>_source_vertex.i32 with with_source). No file when no face qualifies. */
static int cr_write_part(const char *out_dir, const char *name, const MeshBinData *m, const uint8_t *attr,
                         const int32_t *chart, const float *nrm, const uint8_t *reason, uint8_t want, int with_source,
                         size_t *nv_out, size_t *nf_out)
{
    size_t nv = m->nv, nf = m->nf, kv = 0, kf = 0;
    int rc = -1;
    *nv_out = *nf_out = 0;
    int32_t *map = malloc(nv * sizeof *map);
    if (!map) return -1;
    for (size_t v = 0; v < nv; v++) map[v] = -1;
    for (size_t f = 0; f < nf; f++) {
        const int32_t *t = m->faces + 3*f;
        if (reason[t[0]] != want || reason[t[1]] != want || reason[t[2]] != want) continue;
        kf++;
        map[t[0]] = map[t[1]] = map[t[2]] = -2;   /* used by a selected face */
    }
    if (!kf) { free(map); return 0; }
    for (size_t v = 0; v < nv; v++) if (map[v] == -2) map[v] = (int32_t)kv++;
    float *verts = malloc(kv * 3 * sizeof *verts), *n3 = malloc(kv * 3 * sizeof *n3);
    double *uv = malloc(kv * 2 * sizeof *uv);
    int32_t *faces = malloc(kf * 3 * sizeof *faces), *ch = malloc(kv * sizeof *ch), *src = malloc(kv * sizeof *src);
    uint8_t *at = malloc(kv);
    if (!verts || !n3 || !uv || !faces || !ch || !src || !at) goto done;
    for (size_t v = 0; v < nv; v++) {
        if (map[v] < 0) continue;
        size_t k = (size_t)map[v];
        memcpy(verts + 3*k, m->verts + 3*v, 3 * sizeof *verts);
        memcpy(n3 + 3*k, nrm + 3*v, 3 * sizeof *n3);
        uv[2*k] = m->uv64[2*v]; uv[2*k+1] = m->uv64[2*v+1];
        ch[k] = chart[v]; at[k] = attr[v]; src[k] = (int32_t)v;
    }
    kf = 0;
    for (size_t f = 0; f < nf; f++) {
        const int32_t *t = m->faces + 3*f;
        if (reason[t[0]] != want || reason[t[1]] != want || reason[t[2]] != want) continue;
        for (int k = 0; k < 3; k++) faces[3*kf + (size_t)k] = map[t[k]];
        kf++;
    }
    char path[4096];
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", out_dir, name);
    if (MeshBin_write_uv64(path, verts, kv, faces, kf, uv) != 0) { fprintf(stderr, "ERROR: cannot write %s\n", path); goto done; }
    snprintf(path, sizeof path, "%s/%s_surface_attr.u8", out_dir, name);
    if (cr_write(path, at, kv)) goto done;
    snprintf(path, sizeof path, "%s/%s_surface_chart.i32", out_dir, name);
    if (cr_write(path, ch, kv * sizeof *ch)) goto done;
    snprintf(path, sizeof path, "%s/%s_surface_nrm.f32", out_dir, name);
    if (cr_write(path, n3, kv * 3 * sizeof *n3)) goto done;
    if (with_source) {
        snprintf(path, sizeof path, "%s/%s_source_vertex.i32", out_dir, name);
        if (cr_write(path, src, kv * sizeof *src)) goto done;
    }
    *nv_out = kv; *nf_out = kf; rc = 0;
done:
    free(map); free(verts); free(n3); free(uv); free(faces); free(ch); free(src); free(at);
    return rc;
}

/* the frame at half scale, every vertex one pixel: kept first, removals over them */
static int cr_map_png(const char *path, const MeshBinData *m, const uint8_t *attr, const uint8_t *reason)
{
    double u0 = INFINITY, v0 = INFINITY, u1 = -INFINITY, v1 = -INFINITY;
    for (size_t v = 0; v < m->nv; v++) {
        u0 = fmin(u0, m->uv64[2*v]); u1 = fmax(u1, m->uv64[2*v]);
        v0 = fmin(v0, m->uv64[2*v+1]); v1 = fmax(v1, m->uv64[2*v+1]);
    }
    if (!(u1 >= u0) || !(v1 >= v0)) return -1;
    int W = (int)((u1 - u0) / 2.0) + 1, H = (int)((v1 - v0) / 2.0) + 1;
    uint8_t *img = calloc((size_t)W * (size_t)H * 3, 1);
    if (!img) return -1;
    static const uint8_t color[5][3] = { { 0, 0, 0 }, { 0, 255, 120 }, { 255, 220, 0 }, { 255, 0, 255 }, { 0, 200, 255 } };
    for (int pass = 0; pass < 2; pass++) for (size_t v = 0; v < m->nv; v++) {
        if ((reason[v] != CLOSEUP_REPAIR_KEEP) != pass) continue;
        int x = (int)((m->uv64[2*v] - u0) / 2.0), y = (int)((m->uv64[2*v+1] - v0) / 2.0);
        uint8_t *px = img + 3 * ((size_t)y * (size_t)W + (size_t)x);
        if (!pass) {
            uint8_t g = attr[v] == 2 ? 0 : 110;
            px[0] = attr[v] == 2 ? 40 : g; px[1] = attr[v] == 2 ? 60 : g; px[2] = attr[v] == 2 ? 140 : g;
        } else memcpy(px, color[reason[v] <= 4 ? reason[v] : 0], 3);
    }
    int rc = VesPng_write_rgb(path, img, W, H);
    free(img);
    return rc;
}

static int cr_run(const char *in_dir, const char *stem, const char *axis_path, const char *out_dir)
{
    char vmesh[4096], attr_p[4096], chart_p[4096], nrm_p[4096], path[4096];
    snprintf(vmesh, sizeof vmesh, "%s/%s_surface.vmesh", in_dir, stem);
    snprintf(attr_p, sizeof attr_p, "%s/%s_surface_attr.u8", in_dir, stem);
    snprintf(chart_p, sizeof chart_p, "%s/%s_surface_chart.i32", in_dir, stem);
    snprintf(nrm_p, sizeof nrm_p, "%s/%s_surface_nrm.f32", in_dir, stem);
    snprintf(path, sizeof path, "%s/%s_surface.vmesh", out_dir, stem);
    if (cr_exists(path)) { fprintf(stderr, "ERROR: %s exists; a repair never overwrites a surface\n", path); return 3; }
    (void)ves_mkdir(out_dir);
    double t0 = ves_clock_sec();
    MeshBinData m;
    memset(&m, 0, sizeof m);
    if (MeshBin_read_precise_malloc(vmesh, &m) != 0 || !m.uv64) { fprintf(stderr, "ERROR: cannot read %s (float64 UV)\n", vmesh); return 1; }
    size_t nv = m.nv, nf = m.nf;
    uint8_t *attr = malloc(nv), *reason = malloc(nv);
    int32_t *chart = malloc(nv * sizeof *chart);
    float *nrm = malloc(nv * 3 * sizeof *nrm), *inward = malloc(nv * 3 * sizeof *inward);
    Arena_T arena = Arena_new();
    int rc = 1;
    if (!attr || !reason || !chart || !nrm || !inward) goto done;
    if (cr_read(attr_p, attr, nv) || cr_read(chart_p, chart, nv * sizeof *chart) || cr_read(nrm_p, nrm, nv * 3 * sizeof *nrm)) {
        fprintf(stderr, "ERROR: the per-vertex files of %s/%s do not match its %zu vertices\n", in_dir, stem, nv);
        goto done;
    }
    AsmAxis *axis = AsmAxis_load(arena, axis_path);
    if (!axis) { fprintf(stderr, "ERROR: cannot read the axis table %s\n", axis_path); goto done; }
    long n_long = (long)nv, i = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (i = 0; i < n_long; i++) {
        double q[3], s = 0.0, foot[3], d[3], l = 0.0;
        for (int k = 0; k < 3; k++) q[k] = (double)m.verts[(size_t)i * 3 + (size_t)k];
        float *o = inward + (size_t)i * 3;
        o[0] = o[1] = o[2] = 0.0f;
        if (!isfinite(q[0]) || !isfinite(q[1]) || !isfinite(q[2])) continue;
        AsmAxis_project(axis, q, &s, NULL, NULL, foot);
        for (int k = 0; k < 3; k++) { d[k] = foot[k] - q[k]; l += d[k] * d[k]; }
        l = sqrt(l);
        if (l > 1e-9) for (int k = 0; k < 3; k++) o[k] = (float)(d[k] / l);
    }
    CloseupRepairStats st;
    if (CloseupRepair_mark(m.verts, m.uv64, nv, m.faces, nf, attr, nrm, inward, reason, &st) != 0) {
        fprintf(stderr, "ERROR: the repair refused %s (bad input)\n", vmesh);
        goto done;
    }
    size_t kv = 0, kf = 0, xv = 0, xf = 0;
    char extras[1024];
    snprintf(extras, sizeof extras, "%s_extras", stem);
    if (cr_write_part(out_dir, stem, &m, attr, chart, nrm, reason, CLOSEUP_REPAIR_KEEP, 0, &kv, &kf) || !kf) {
        fprintf(stderr, "ERROR: no kept surface written for %s\n", vmesh);
        goto done;
    }
    if (cr_write_part(out_dir, extras, &m, attr, chart, nrm, reason, CLOSEUP_REPAIR_FLAP, 1, &xv, &xf)) goto done;
    snprintf(path, sizeof path, "%s/%s_repair_reason.u8", out_dir, stem);
    if (cr_write(path, reason, nv)) goto done;
    snprintf(path, sizeof path, "%s/%s_repair_map.png", out_dir, stem);
    if (cr_map_png(path, &m, attr, reason)) fprintf(stderr, "WARNING: cannot write %s\n", path);
    /* the report: inputs (sha256), constants, counts */
    const char *inputs[5] = { vmesh, attr_p, chart_p, nrm_p, axis_path };
    const char *names[5] = { "vmesh", "attr", "chart", "nrm", "axis" };
    snprintf(path, sizeof path, "%s/%s_repair_report.json", out_dir, stem);
    FILE *fp = fopen(path, "wb");
    if (!fp) goto done;
    fprintf(fp, "{\"schema\":\"closeup-repair-v1\",\"stem\":\"%s\",\"inputs\":{", stem);
    for (int k = 0; k < 5; k++) {
        char hex[65] = { 0 };
        long long bytes = 0;
        if (cr_sha256(inputs[k], hex, &bytes)) { fclose(fp); goto done; }
        fprintf(fp, "%s\"%s\":{\"path\":\"", k ? "," : "", names[k]);
        for (const char *c = inputs[k]; *c; c++) fputc(*c == '\\' ? '/' : *c, fp);
        fprintf(fp, "\",\"bytes\":%lld,\"sha256\":\"%s\"}", bytes, hex);
    }
    fprintf(fp, "},\"constants\":{\"dup_vox\":%.17g,\"dup_px\":%.17g,\"flap_min_verts\":%d,\"flap_facing\":%.17g,\"fill_facing\":%.17g},"
        "\"vertices\":%zu,\"faces\":%zu,\"observed\":%zu,\"fill\":%zu,\"undefined_direction\":%zu,"
        "\"inward_flipped\":%s,\"observed_facing_axis\":%.6f,"
        "\"removed\":{\"duplicate\":%zu,\"wrong_side_fill\":%zu,\"flap\":%zu,\"orphan\":%zu},"
        "\"duplicate_pairs\":{\"fill_over_observed\":%zu,\"fill_over_fill\":%zu,\"observed_over_observed_kept\":%zu},"
        "\"flap_components\":%zu,\"wrong_side_observed_specks_kept\":%zu,\"orphan_components\":%zu,"
        "\"kept\":{\"vertices\":%zu,\"faces\":%zu},\"extras\":{\"vertices\":%zu,\"faces\":%zu},"
        "\"mark_seconds\":%.3f,\"seconds\":%.3f}\n",
        CLOSEUP_REPAIR_DUP_VOX, CLOSEUP_REPAIR_DUP_PX, CLOSEUP_REPAIR_FLAP_MIN_VERTS, CLOSEUP_REPAIR_FLAP_FACING,
        CLOSEUP_REPAIR_FILL_FACING, nv, nf, st.observed, st.fill, st.undefined, st.inward_flipped ? "true" : "false",
        st.observed_facing_axis, st.duplicate, st.wrong_fill, st.flap, st.orphan, st.pairs_fill_observed,
        st.pairs_fill_fill, st.pairs_observed_observed, st.flap_components, st.flap_specks, st.orphan_components,
        kv, kf, xv, xf, st.seconds, ves_clock_sec() - t0);
    if (fclose(fp)) goto done;
    fprintf(stderr, "[closeup_repair] %s: %zu vertices (%zu observed, %zu fill); removed %zu duplicate, %zu wrong-side fill, "
        "%zu orphan; %zu flap vertices in %zu flaps to extras (%zu specks kept); observed facing the axis %.4f%s; kept %zu faces "
        "of %zu (%.1f s)\n", stem, nv, st.observed, st.fill, st.duplicate, st.wrong_fill, st.orphan, st.flap,
        st.flap_components, st.flap_specks, st.observed_facing_axis, st.inward_flipped ? " (normals read flipped)" : "",
        kf, nf, ves_clock_sec() - t0);
    rc = 0;
done:
    MeshBin_dispose(&m);
    free(attr); free(reason); free(chart); free(nrm); free(inward);
    Arena_dispose(&arena);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        int fails = CloseupRepair_selftest();
        fprintf(stderr, "closeup_repair selftest: %d failure(s)\n", fails);
        return fails ? 1 : 0;
    }
    if (argc != 5) {
        fprintf(stderr, "usage: closeup_repair <in_dir> <stem> <axis.csv> <out_dir>\n       closeup_repair --selftest\n");
        return 2;
    }
    return cr_run(argv[1], argv[2], argv[3], argv[4]);
}
