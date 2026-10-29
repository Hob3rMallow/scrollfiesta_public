/* Fitted ribbon OBJ emission shared by fast parameterizers. */
#define _USE_MATH_DEFINES
#include "ribbon.h"
#include "../common/mesh_bin.h"
#include "../common/ves_platform.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif


/* Ribbon grid -> OBJ/VMESH.  The parameterizer owns topology: every defined
 * four-corner grid cell is emitted as exactly two triangles.  This writer does
 * not reinterpret geometry, identity, provenance, or intersections. */

static int ribbon_replace_file(const char *temporary, const char *destination)
{
#ifdef _WIN32
    return MoveFileExA(temporary, destination,
                       MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)
         ? 0 : -1;
#else
    return rename(temporary, destination);
#endif
}

static int ribbon_claimant_path(const char *path, char *out, size_t capacity)
{
    const char suffix[] = "_claimant_chart.i32";
    size_t n, stem;
    if (path == NULL || out == NULL || capacity == 0) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof(suffix) || stem + sizeof(suffix) > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof(suffix));
    return 0;
}

static int ribbon_island_path(const char *path, char *out, size_t capacity)
{
    /* Compatibility alias retained for older consumers.  New code should use
     * the explicit reconstruction-component companion below. */
    const char suffix[] = "_claimant_island.i32";
    size_t n, stem;
    if (path == NULL || out == NULL || capacity == 0) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof(suffix) || stem + sizeof(suffix) > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof(suffix));
    return 0;
}

static int ribbon_reconstruction_component_path(const char *path, char *out,
                                                size_t capacity)
{
    const char suffix[] = "_reconstruction_component.i32";
    size_t n, stem;
    if (path == NULL || out == NULL || capacity == 0) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof(suffix) || stem + sizeof(suffix) > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof(suffix));
    return 0;
}

static int ribbon_material_identity_path(const char *path, char *out,
                                         size_t capacity)
{
    const char suffix[] = "_material_identity.i32";
    size_t n, stem;
    if (path == NULL || out == NULL || capacity == 0) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof(suffix) || stem + sizeof(suffix) > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof(suffix));
    return 0;
}

static int ribbon_support_path(const char *path, char *out, size_t capacity)
{
    const char suffix[] = "_support.u8";
    size_t n, stem;
    if (path == NULL || out == NULL || capacity == 0) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof(suffix) || stem + sizeof(suffix) > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof(suffix));
    return 0;
}

static int ribbon_phase_path(const char *path, char *out, size_t capacity)
{
    const char suffix[] = "_phase.f32";
    size_t n, stem;
    if (path == NULL || out == NULL || capacity == 0) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof(suffix) || stem + sizeof(suffix) > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof(suffix));
    return 0;
}

static int ribbon_projective_layer_path(const char *path, size_t layer,
                                        char *out, size_t capacity)
{
    size_t len = 0, stem = 0;
    const char *dot = NULL;
    int written = 0;
    if (path == NULL || out == NULL || capacity == 0 || layer == 0) return -1;
    len = strlen(path);
    dot = strrchr(path, '.');
    stem = dot != NULL ? (size_t)(dot - path) : len;
    if (stem >= capacity) return -1;
    memcpy(out, path, stem);
    out[stem] = '\0';
    if (layer == 1)
        written = snprintf(out + stem, capacity - stem, "_extras%s",
                           path + stem);
    else
        written = snprintf(out + stem, capacity - stem,
                           "_extras_layer%02zu%s", layer, path + stem);
    return written >= 0 && (size_t)written < capacity - stem ? 0 : -1;
}

typedef struct { size_t vertices, faces, point_only_vertices; } RibbonLayerCount;

/* A point claimant is evidence even when it cannot own a triangle. Keep a
 * separate evidence record; never manufacture a degenerate VMESH face. Float
 * bits use the same final f32 conversion as the surface writer. */
typedef struct {
    uint32_t u, v, z, y, x, phase;
    int32_t chart, material, reconstruction;
    unsigned support;
} RibbonPointRecord;

static uint32_t ribbon_f32_bits(float value)
{
    uint32_t bits;
    memcpy(&bits, &value, sizeof bits);
    return bits;
}

static int ribbon_point_record(const RibbonResult *R, size_t layer, size_t np,
                               size_t k, size_t j, RibbonPointRecord *record)
{
    size_t row, slot;
    const float *p;
    int32_t chart = -1;
    if (!R || !record || !np || !R->nu || !R->grid_projective ||
        R->nk % np || layer >= R->nk / np || k >= np || j >= R->nu ||
        R->nk > SIZE_MAX / R->nu || !R->grid_pos || !R->grid_present ||
        !R->grid_phi || !R->grid_material || !R->grid_island) return -1;
    row = layer * np + k;
    slot = row * R->nu + j;
    if (!R->grid_present[slot]) return -1;
    p = R->grid_pos + 3 * slot;
    if (!isfinite(p[0]) || !isfinite(p[1]) || !isfinite(p[2]) ||
        !isfinite(R->grid_phi[slot])) return -1;
    record->support = !R->grid_valid || R->grid_valid[slot] ? 255 : 0;
    if (record->support) {
        size_t first, last;
        if (!R->grid_chart_row_offsets || !R->grid_chart_runs) return -1;
        first = R->grid_chart_row_offsets[row];
        last = R->grid_chart_row_offsets[row + 1];
        if (first > last || last > R->grid_chart_run_count) return -1;
        for (size_t q = first; q < last; q++) {
            const RibbonChartRun *run = R->grid_chart_runs + q;
            if (run->first_col > j || run->last_col < j) continue;
            if (run->source_chart < 0 || (chart >= 0 && chart != run->source_chart)) return -1;
            chart = run->source_chart;
        }
        if (chart < 0) return -1;
    }
    record->u = ribbon_f32_bits((float)(R->grid_u_origin + (double)j * R->grid_du));
    record->v = ribbon_f32_bits((float)(R->grid_v_origin + (double)k * R->grid_dv));
    record->z = ribbon_f32_bits(p[0]); record->y = ribbon_f32_bits(p[1]);
    record->x = ribbon_f32_bits(p[2]); record->phase = ribbon_f32_bits(R->grid_phi[slot]);
    record->chart = chart; record->material = R->grid_material[slot];
    record->reconstruction = R->grid_island[slot];
    return 0;
}

static int ribbon_point_evidence_path(const char *path, char *out, size_t capacity)
{
    size_t n = strlen(path), stem = n > 4 && !strcmp(path + n - 4, ".obj") ? n - 4 : n;
    const char suffix[] = "_point_evidence.json";
    if (stem >= capacity || sizeof suffix > capacity - stem) return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof suffix);
    return 0;
}

static int ribbon_write_point_evidence(const char *path, const RibbonResult *R,
                                       size_t layer, size_t np, size_t expected)
{
    char final_path[2048], temporary[2060];
    FILE *file;
    size_t count = 0;
    int rc = -1, n;
    if (ribbon_point_evidence_path(path, final_path, sizeof final_path)) return -1;
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary) return -1;
    file = fopen(temporary, "wb");
    if (!file) return -1;
    fprintf(file, "{\"schema\":\"ribbon-point-evidence-v1\",\"source_layer\":%zu,"
            "\"coordinate_frame\":\"carried_unwrap\",\"count\":%zu,\"points\":[",
            R->grid_layer_base + layer, expected);
    for (size_t k = 0; k < np; k++) for (size_t j = 0; j < R->nu; j++) {
        RibbonPointRecord point;
        if (!R->grid_present[((layer * np + k) * R->nu) + j]) continue;
        if (ribbon_point_record(R, layer, np, k, j, &point)) goto done;
        fprintf(file, "%s{\"u\":%u,\"v\":%u,\"z\":%u,\"y\":%u,\"x\":%u,"
                "\"phase\":%u,\"chart\":%d,\"material\":%d,\"reconstruction\":%d,\"support\":%u}",
                count ? "," : "", point.u, point.v, point.z, point.y, point.x,
                point.phase, point.chart, point.material, point.reconstruction, point.support);
        count++;
    }
    if (count != expected) goto done;
    fprintf(file, "],\"float_encoding\":\"IEEE754-f32-bits-as-u32\",\"complete\":true}\n");
    rc = ferror(file) ? -1 : 0;
done:
    if (fclose(file)) rc = -1;
    if (!rc) rc = ribbon_replace_file(temporary, final_path);
    if (rc) remove(temporary);
    return rc;
}

static int ribbon_inventory_path(const char *path, char *out, size_t capacity)
{
    const char suffix[] = "_layers.json";
    size_t n, stem;
    if (path == NULL || out == NULL) return -1;
    n = strlen(path);
    stem = n > 4 && strcmp(path + n - 4, ".obj") == 0 ? n - 4 : n;
    if (stem > SIZE_MAX - sizeof suffix || stem + sizeof suffix > capacity)
        return -1;
    memcpy(out, path, stem);
    memcpy(out + stem, suffix, sizeof suffix);
    return 0;
}

/* Publish only after every layer has been considered and every surface write
 * succeeded. A claimant multiplicity histogram is NOT an emission inventory. */
static int ribbon_write_inventory(const char *path, const RibbonLayerCount *counts,
                                  size_t layers, size_t first, size_t source_layers)
{
    char temporary[2056];
    FILE *file;
    int n = snprintf(temporary, sizeof temporary, "%s.tmp", path), rc = -1;
    if (n < 0 || (size_t)n >= sizeof temporary) return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    fprintf(file, "{\"schema\":\"ribbon-emitted-layers-v1\",\"grid_layers\":%zu,", layers);
    if (source_layers > layers)
        fprintf(file, "\"partial_atlas\":true,\"source_peel_layer\":%zu,\"source_peel_layers\":%zu,",
                first, source_layers);
    fprintf(file, "\"layers\":[");
    for (size_t i = 0; i < layers; i++)
        fprintf(file, "%s{\"layer\":%zu,\"vertices\":%zu,\"faces\":%zu,\"point_only_vertices\":%zu,\"point_evidence_preserved\":true}",
                i ? "," : "", i, counts[i].vertices, counts[i].faces,
                counts[i].point_only_vertices);
    fprintf(file, "],\"complete\":true}\n");
    rc = ferror(file) ? -1 : 0;
    if (fclose(file) != 0) rc = -1;
    if (rc == 0) rc = ribbon_replace_file(temporary, path);
    if (rc != 0) remove(temporary);
    return rc;
}

/* A point-only color has no VMESH surface, but an older run at the same path
 * may have emitted one.  Remove only this exact output family so the Stage-4
 * atlas cannot mistake stale geometry for the current Stage-3 result. */
static int ribbon_remove_output_family(const char *path)
{
    char target[8][2048];
    size_t n = 0;
    if (path == NULL) return -1;
    if (MeshBin_companion_path(path, target[n], sizeof target[n]) != 0)
        return -1;
    n++;
    if (ribbon_support_path(path, target[n], sizeof target[n]) != 0) return -1;
    n++;
    if (ribbon_phase_path(path, target[n], sizeof target[n]) != 0) return -1;
    n++;
    if (ribbon_claimant_path(path, target[n], sizeof target[n]) != 0) return -1;
    n++;
    if (ribbon_reconstruction_component_path(
            path, target[n], sizeof target[n]) != 0)
        return -1;
    n++;
    if (ribbon_material_identity_path(path, target[n], sizeof target[n]) != 0)
        return -1;
    n++;
    if (ribbon_island_path(path, target[n], sizeof target[n]) != 0) return -1;
    n++;
    if (ribbon_point_evidence_path(path, target[n], sizeof target[n]) != 0) return -1;
    n++;
    for (size_t i = 0; i < n; i++)
        if (remove(target[i]) != 0 && errno != ENOENT) return -1;
    if (remove(path) != 0 && errno != ENOENT) return -1;
    return 0;
}

/* One little-endian int32 per final VMESH vertex, in the exact same order.
 * -1 means derived/unsupported (including center-fan vertices).  Raw payload
 * plus strict VMESH vertex-count validation keeps the sidecar compact enough
 * for whole-scroll runs while making stale/mismatched files fail closed. */
static int ribbon_write_claimant_charts(const char *path,
                                        const int32_t *chart, size_t count)
{
    char final_path[2048], temporary[2060];
    FILE *file = NULL;
    int n;
    if (ribbon_claimant_path(path, final_path, sizeof final_path) != 0)
        return -1;
    if (chart == NULL) {
        if (remove(final_path) != 0 && errno != ENOENT) return -1;
        return 0;
    }
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary ||
        count > SIZE_MAX / sizeof(*chart))
        return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    if (fwrite(chart, sizeof(*chart), count, file) != count ||
        fflush(file) != 0) {
        fclose(file);
        file = NULL;
        remove(temporary);
        return -1;
    }
    if (fclose(file) != 0) {
        file = NULL;
        remove(temporary);
        return -1;
    }
    file = NULL;
    if (ribbon_replace_file(temporary, final_path) != 0) {
        remove(temporary);
        return -1;
    }
    fprintf(stderr,
            "  ribbon output claimant provenance: %zu vertex label(s) -> %s\n",
            count, final_path);
    return 0;
}

static int ribbon_write_claimant_islands(const char *path,
                                          const int32_t *island, size_t count)
{
    char final_path[2048], temporary[2060];
    FILE *file = NULL;
    int n;
    if (ribbon_island_path(path, final_path, sizeof final_path) != 0)
        return -1;
    if (island == NULL) {
        if (remove(final_path) != 0 && errno != ENOENT) return -1;
        return 0;
    }
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary ||
        count > SIZE_MAX / sizeof(*island))
        return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    if (fwrite(island, sizeof(*island), count, file) != count ||
        fflush(file) != 0) {
        fclose(file);
        file = NULL;
        remove(temporary);
        return -1;
    }
    if (fclose(file) != 0) {
        file = NULL;
        remove(temporary);
        return -1;
    }
    file = NULL;
    if (ribbon_replace_file(temporary, final_path) != 0) {
        remove(temporary);
        return -1;
    }
    fprintf(stderr,
            "  ribbon output legacy identity alias: %zu vertex label(s) -> %s\n",
            count, final_path);
    return 0;
}

static int ribbon_write_reconstruction_components(const char *path,
                                                   const int32_t *component,
                                                   size_t count)
{
    char final_path[2048], temporary[2060];
    FILE *file = NULL;
    int n;
    if (ribbon_reconstruction_component_path(
            path, final_path, sizeof final_path) != 0)
        return -1;
    if (component == NULL) {
        if (remove(final_path) != 0 && errno != ENOENT) return -1;
        return 0;
    }
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary ||
        count > SIZE_MAX / sizeof(*component))
        return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    if (fwrite(component, sizeof(*component), count, file) != count ||
        fflush(file) != 0) {
        fclose(file);
        remove(temporary);
        return -1;
    }
    if (fclose(file) != 0) {
        remove(temporary);
        return -1;
    }
    if (ribbon_replace_file(temporary, final_path) != 0) {
        remove(temporary);
        return -1;
    }
    fprintf(stderr,
            "  ribbon output reconstruction identity: %zu vertex label(s) -> %s\n",
            count, final_path);
    return 0;
}

static int ribbon_write_material_identities(const char *path,
                                            const int32_t *material,
                                            size_t count)
{
    char final_path[2048], temporary[2060];
    FILE *file = NULL;
    int n;
    if (ribbon_material_identity_path(path, final_path, sizeof final_path) != 0)
        return -1;
    if (material == NULL) {
        if (remove(final_path) != 0 && errno != ENOENT) return -1;
        return 0;
    }
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary ||
        count > SIZE_MAX / sizeof(*material))
        return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    if (fwrite(material, sizeof(*material), count, file) != count ||
        fflush(file) != 0) {
        fclose(file);
        remove(temporary);
        return -1;
    }
    if (fclose(file) != 0) {
        remove(temporary);
        return -1;
    }
    if (ribbon_replace_file(temporary, final_path) != 0) {
        remove(temporary);
        return -1;
    }
    fprintf(stderr,
            "  ribbon output material lineage: %zu vertex label(s) -> %s\n",
            count, final_path);
    return 0;
}

/* One byte per final VMESH vertex in exactly the VMESH order: 255 is an
 * unambiguous measured Dirichlet sample, 0 is generated continuation.  The
 * consumer validates the payload length against the authoritative VMESH, so a
 * stale confidence file fails closed instead of silently pinning filled data. */
static int ribbon_write_support(const char *path, const uint8_t *support,
                                size_t count)
{
    char final_path[2048], temporary[2060];
    FILE *file = NULL;
    int n;
    size_t supported = 0;
    if (ribbon_support_path(path, final_path, sizeof final_path) != 0 ||
        support == NULL)
        return -1;
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary) return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    if (fwrite(support, sizeof(*support), count, file) != count ||
        fflush(file) != 0) {
        fclose(file);
        remove(temporary);
        return -1;
    }
    if (fclose(file) != 0) {
        remove(temporary);
        return -1;
    }
    if (ribbon_replace_file(temporary, final_path) != 0) {
        remove(temporary);
        return -1;
    }
    for (size_t v = 0; v < count; v++) supported += support[v] == 255;
    fprintf(stderr,
            "  ribbon output confidence: %zu fixed, %zu generated -> %s\n",
            supported, count - supported, final_path);
    return 0;
}

/* One finite float32 radian value per final VMESH vertex.  The confidence
 * companion, not a floating-point sentinel, distinguishes measured
 * StrokeStrip claims from generated continuation. */
static int ribbon_write_phase(const char *path, const float *phase, size_t count)
{
    char final_path[2048], temporary[2060];
    FILE *file = NULL;
    int n;
    if (ribbon_phase_path(path, final_path, sizeof final_path) != 0)
        return -1;
    if (phase == NULL) {
        if (remove(final_path) != 0 && errno != ENOENT) return -1;
        return 0;
    }
    for (size_t v = 0; v < count; v++)
        if (!isfinite((double)phase[v])) return -1;
    n = snprintf(temporary, sizeof temporary, "%s.tmp", final_path);
    if (n < 0 || (size_t)n >= sizeof temporary ||
        count > SIZE_MAX / sizeof(*phase)) return -1;
    file = fopen(temporary, "wb");
    if (file == NULL) return -1;
    if (fwrite(phase, sizeof(*phase), count, file) != count ||
        fflush(file) != 0) {
        fclose(file); remove(temporary); return -1;
    }
    if (fclose(file) != 0) { remove(temporary); return -1; }
    if (ribbon_replace_file(temporary, final_path) != 0) {
        remove(temporary); return -1;
    }
    fprintf(stderr,
            "  ribbon output lifted phase: %zu finite value(s) -> %s\n",
            count, final_path);
    return 0;
}

static int write_ribbon_obj(const char *path, const RibbonResult *R,
                            int write_obj, int v_rebase,
                            size_t *out_nv, size_t *out_nf,
                            size_t *out_supported, size_t *out_generated,
                            size_t *out_atlas_cols,
                            size_t *out_atlas_runs,
                            size_t *out_empty_cols_removed)
{
    size_t nu, nk, nslot, nvout = 0, nfout = 0;
    size_t atlas_cols = 0, atlas_runs = 0, empty_cols_removed = 0;
    int32_t *vid = NULL, *faces = NULL;
    int32_t *vertex_chart = NULL, *vertex_island = NULL;
    int32_t *vertex_material = NULL;
    int32_t *column_map = NULL, *run_of_col = NULL;
    uint8_t *column_used = NULL, *vertex_support = NULL;
    float *verts = NULL, *uv = NULL, *vertex_phase = NULL;
    FILE *ff = NULL;
    int rc = -1;
    double total_t0 = ves_clock_sec(), stage_t0 = total_t0;

    *out_nv = 0; *out_nf = 0;
    *out_supported = 0; *out_generated = 0;
    *out_atlas_cols = 0; *out_atlas_runs = 0;
    *out_empty_cols_removed = 0;
    if (R->grid_pos == NULL || R->grid_present == NULL ||
        R->nu < 2 || R->nk < 2) return -1;
    if (path == NULL) return -1;
    {
        char points[2048];
        if (ribbon_point_evidence_path(path, points, sizeof points) != 0 ||
            (remove(points) != 0 && errno != ENOENT)) return -1;
    }
    nu = R->nu; nk = R->nk;
    if (nu > SIZE_MAX / nk) return -1;
    nslot = nu * nk;
    if (nslot > SIZE_MAX / sizeof(*vid)) return -1;
    vid = (int32_t *)malloc(nslot * sizeof(*vid));
    if (!vid) goto done;
    for (size_t q = 0; q < nslot; q++) vid[q] = -1;

    /* Compact the parameterizer's occupied grid vertices.  Occupancy is
     * explicit: geometry values are never interpreted as sentinels. */
    for (size_t k = 0; k < nk; k++) {
        for (size_t j = 0; j < nu; j++) {
            const float *p = &R->grid_pos[(k*nu + j)*3];
            int present = R->grid_present[k*nu + j] != 0;
            if (!present) continue;
            if (!isfinite((double)p[0]) || !isfinite((double)p[1]) ||
                !isfinite((double)p[2])) {
                fprintf(stderr,
                        "ribbon: non-finite geometry at occupied grid cell "
                        "(%zu,%zu)\n", k, j);
                goto done;
            }
            if (nvout >= (size_t)INT32_MAX) goto done;
            vid[k*nu + j] = (int32_t)nvout++;
        }
    }
    if (nvout == 0 || nvout > SIZE_MAX / (3 * sizeof(*verts)) ||
        nvout > SIZE_MAX / (2 * sizeof(*uv))) goto done;
    verts = (float *)malloc(nvout * 3 * sizeof(*verts));
    uv = (float *)malloc(nvout * 2 * sizeof(*uv));
    vertex_support = (uint8_t *)malloc(nvout * sizeof(*vertex_support));
    if (R->grid_phi != NULL)
        vertex_phase = (float *)malloc(nvout * sizeof(*vertex_phase));
    if (!verts || !uv || !vertex_support ||
        (R->grid_phi != NULL && vertex_phase == NULL)) goto done;
    if (R->grid_island != NULL) {
        vertex_island = (int32_t *)malloc(nvout * sizeof(*vertex_island));
        if (vertex_island == NULL) goto done;
    }
    if (R->grid_material != NULL) {
        vertex_material =
            (int32_t *)malloc(nvout * sizeof(*vertex_material));
        if (vertex_material == NULL) goto done;
    }
    if ((R->grid_chart_row_offsets == NULL) !=
        (R->grid_chart_runs == NULL)) goto done;
    if (R->grid_chart_row_offsets != NULL) {
        vertex_chart = (int32_t *)malloc(nvout * sizeof(*vertex_chart));
        if (vertex_chart == NULL) goto done;
        for (size_t v = 0; v < nvout; v++) vertex_chart[v] = -1;
    }
    for (size_t k = 0; k < nk; k++) {
        for (size_t j = 0; j < nu; j++) {
            int32_t v = vid[k*nu+j];
            if (v < 0) continue;
            memcpy(&verts[(size_t)v*3],&R->grid_pos[(k*nu+j)*3],3*sizeof(float));
            uv[(size_t)v*2] = (float)(R->grid_projective
                                  ? R->grid_u_origin +
                                    (double)j * (double)R->grid_du
                                  : (double)j * (double)R->grid_du);
            uv[(size_t)v*2+1] = (float)(R->grid_projective
                                    ? R->grid_v_origin +
                                      (double)k * (double)R->grid_dv
                                    : (double)k * (double)R->grid_dv);
            vertex_support[(size_t)v] =
                R->grid_valid == NULL || R->grid_valid[k*nu+j]
                    ? (uint8_t)255 : (uint8_t)0;
            if (vertex_phase != NULL)
                vertex_phase[(size_t)v] = R->grid_phi[k*nu+j];
            if (vertex_island != NULL)
                vertex_island[(size_t)v] = R->grid_island[k*nu+j];
            if (vertex_material != NULL)
                vertex_material[(size_t)v] = R->grid_material[k*nu+j];
        }
    }
    if (vertex_chart != NULL) {
        if (R->grid_chart_row_offsets[nk] != R->grid_chart_run_count)
            goto done;
        for (size_t k = 0; k < nk; k++) {
            size_t r0 = R->grid_chart_row_offsets[k];
            size_t r1 = R->grid_chart_row_offsets[k + 1];
            if (r0 > r1 || r1 > R->grid_chart_run_count) goto done;
            for (size_t r = r0; r < r1; r++) {
                const RibbonChartRun *run = &R->grid_chart_runs[r];
                if (run->source_chart < 0 || run->first_col > run->last_col ||
                    run->last_col >= nu) goto done;
                for (size_t j = run->first_col; j <= run->last_col; j++) {
                    int32_t v = vid[k*nu + j];
                    if (v < 0 || (R->grid_valid != NULL &&
                                  !R->grid_valid[k*nu + j]) ||
                        (vertex_chart[(size_t)v] >= 0 &&
                         vertex_chart[(size_t)v] != run->source_chart))
                        goto done;
                    vertex_chart[(size_t)v] = run->source_chart;
                }
            }
        }
    }
    /* Face emission is a direct serialization of fitted-grid topology.  The
     * parameterizer decides which triangles a cell may emit (R->grid_quad,
     * RIB_CELL_*); this writer maps corners through `vid`, validates them, and
     * writes.  It applies no geometry policy of its own -- presence alone used
     * to be the entire test, and that is what drew 35,447 edges across the gap
     * between physical wraps on the 4x5x5. */
    if (R->grid_quad == NULL) {
        fprintf(stderr,
                "BUG: fitted grid carries no emitted-topology decision; "
                "refusing to guess one in the writer\n");
        goto done;
    }
    {
        /* corner order a,b,c,d; one row per RIB_CELL_* bit */
        static const int tri[4][3] = { {0,1,2}, {1,3,2}, {0,1,3}, {0,3,2} };
        int32_t corner[4];
        for (size_t k = 0; k + 1 < nk; k++)
            for (size_t j = 0; j + 1 < nu; j++) {
                uint8_t code = R->grid_quad[k*(nu-1) + j];
                int t = 0;
                for (t = 0; t < 4; t++) if (code & (1 << t)) nfout++;
            }
        if (nfout == 0 || nfout > SIZE_MAX / (3 * sizeof(*faces))) goto done;
        faces = (int32_t *)malloc(nfout * 3 * sizeof(*faces));
        if (!faces) goto done;
        nfout = 0;
        for (size_t k = 0; k + 1 < nk; k++) {
            for (size_t j = 0; j + 1 < nu; j++) {
                uint8_t code = R->grid_quad[k*(nu-1) + j];
                int t = 0;
                if (code == 0) continue;
                corner[0] = vid[k*nu + j];
                corner[1] = vid[k*nu + j + 1];
                corner[2] = vid[(k+1)*nu + j];
                corner[3] = vid[(k+1)*nu + j + 1];
                for (t = 0; t < 4; t++) {
                    int e = 0;
                    if (!(code & (1 << t))) continue;
                    for (e = 0; e < 3; e++) {
                        int32_t v = corner[tri[t][e]];
                        if (v < 0) {
                            fprintf(stderr,
                                    "BUG: emitted cell (%zu,%zu) references an "
                                    "absent corner\n", k, j);
                            goto done;
                        }
                        faces[nfout*3 + (size_t)e] = v;
                    }
                    nfout++;
                }
            }
        }
    }
    free(vid); vid = NULL;
    fprintf(stderr,
            "  ribbon output grid materialize: %zu vertices, %zu faces (%.2fs)\n",
            nvout, nfout, ves_clock_sec()-stage_t0);

    /* Pack the parameterizer's occupied columns.  A completely unused U column
     * cannot be crossed by a fitted
     * face (grid faces only join adjacent columns), so removing such columns is
     * an additive translation per disconnected run: every triangle metric and
     * the winding order inside every surviving ribbon piece remain exact.
     * Two blank grid columns separate runs in the final atlas, enough to keep
     * raster footprints disjoint without paying the old ten-column fit gutter
     * hundreds of times for micro-fragments that emitted nothing. */
    column_used = (uint8_t *)calloc(nu ? nu : 1, 1);
    column_map = (int32_t *)malloc((nu ? nu : 1) * sizeof(*column_map));
    if (!column_used || !column_map) goto done;
    for (size_t j = 0; j < nu; j++) column_map[j] = -1;
    for (size_t v = 0; v < nvout; v++) {
        double relative_u = (double)uv[v*2] -
            (R->grid_projective ? R->grid_u_origin : 0.0);
        long j = lround(relative_u / (double)R->grid_du);
        if (j < 0 || (size_t)j >= nu) goto done;
        column_used[(size_t)j] = 1;
    }
    /* atlas cursor can exceed nu by 2 gutter columns per run */
    run_of_col = (int32_t *)malloc((nu * 3 + 2) * sizeof(*run_of_col));
    if (!run_of_col) goto done;
    for (size_t j = 0; j < nu * 3 + 2; j++) run_of_col[j] = -1;
    {
        const size_t gutter_columns = 2;
        size_t cursor = 0;
        int inside = 0;
        for (size_t j = 0; j < nu; j++) {
            if (!column_used[j]) { inside = 0; continue; }
            if (!inside) {
                if (!R->grid_projective && atlas_runs > 0)
                    cursor += gutter_columns;
                atlas_runs++;
                inside = 1;
            }
            if (R->grid_projective) {
                if (j > (size_t)INT32_MAX) goto done;
                run_of_col[j] = (int32_t)(atlas_runs - 1);
                column_map[j] = (int32_t)j;
            } else {
                if (cursor > (size_t)INT32_MAX) goto done;
                run_of_col[cursor] = (int32_t)(atlas_runs - 1);
                column_map[j] = (int32_t)cursor++;
            }
        }
        atlas_cols = R->grid_projective ? nu : cursor;
        empty_cols_removed = R->grid_projective ? 0
                                                : (nu > atlas_cols
                                                   ? nu - atlas_cols : 0);
    }
    for (size_t v = 0; v < nvout; v++) {
        double relative_u = (double)uv[v*2] -
            (R->grid_projective ? R->grid_u_origin : 0.0);
        long old_col = lround(relative_u / (double)R->grid_du);
        int32_t new_col = column_map[(size_t)old_col];
        if (new_col < 0) goto done;
        if (!R->grid_projective)
            uv[v*2] = (float)new_col * R->grid_du;
    }
    char runs_path[2048];
    {
        const char suffix[] = "_runs.json";
        size_t plen = strlen(path);
        size_t stem = plen > 4 && strcmp(path + plen - 4, ".obj") == 0
                    ? plen - 4 : plen;
        if (stem + sizeof suffix > sizeof runs_path) goto done;
        memcpy(runs_path, path, stem);
        memcpy(runs_path + stem, suffix, sizeof suffix);
    }
    if (!v_rebase && remove(runs_path) != 0 && errno != ENOENT) {
        fprintf(stderr, "  ribbon output: cannot remove stale v-rebase "
                "sidecar %s\n", runs_path);
        goto done;
    }
    if (v_rebase && R->grid_projective) {
        fprintf(stderr,
                "ribbon: v_rebase is incompatible with a projective grid\n");
        goto done;
    }
    if (v_rebase && atlas_runs > 0) {
        /* Per-run v rebase: faces never span runs (the gutter columns are
         * empty), so dropping each run to its own v minimum is a pure per-run
         * translation.  The sidecar records the offsets for z recovery. */
        double *run_vmin = (double *)malloc(atlas_runs * sizeof(*run_vmin));
        int32_t *run_c0 = (int32_t *)malloc(atlas_runs * sizeof(*run_c0));
        int32_t *run_c1 = (int32_t *)malloc(atlas_runs * sizeof(*run_c1));
        if (!run_vmin || !run_c0 || !run_c1) {
            free(run_vmin); free(run_c0); free(run_c1);
            goto done;
        }
        for (size_t r = 0; r < atlas_runs; r++) {
            run_vmin[r] = 1e300;
            run_c0[r] = INT32_MAX;
            run_c1[r] = -1;
        }
        for (size_t v = 0; v < nvout; v++) {
            long col = lround((double)uv[v*2] / (double)R->grid_du);
            int32_t r = run_of_col[(size_t)col];
            if (r < 0 || (size_t)r >= atlas_runs) continue;
            if ((double)uv[v*2+1] < run_vmin[r]) run_vmin[r] = (double)uv[v*2+1];
            if ((int32_t)col < run_c0[r]) run_c0[r] = (int32_t)col;
            if ((int32_t)col > run_c1[r]) run_c1[r] = (int32_t)col;
        }
        for (size_t v = 0; v < nvout; v++) {
            long col = lround((double)uv[v*2] / (double)R->grid_du);
            int32_t r = run_of_col[(size_t)col];
            if (r < 0 || (size_t)r >= atlas_runs || run_vmin[r] > 1e299)
                continue;
            uv[v*2+1] = (float)((double)uv[v*2+1] - run_vmin[r]);
        }
        {
            FILE *jf = NULL;
            jf = fopen(runs_path, "w");
            if (jf != NULL) {
                fprintf(jf, "{ \"grid_du\": %.9g, \"grid_dv\": %.9g, "
                        "\"runs\": [\n",
                        (double)R->grid_du, (double)R->grid_dv);
                for (size_t r = 0; r < atlas_runs; r++)
                    fprintf(jf, "  { \"run\": %zu, \"col0\": %d, "
                            "\"col1\": %d, \"v_offset\": %.9g }%s\n",
                            r, run_c0[r], run_c1[r],
                            run_vmin[r] > 1e299 ? 0.0 : run_vmin[r],
                            r + 1 < atlas_runs ? "," : "");
                fprintf(jf, "] }\n");
                fclose(jf);
            }
            fprintf(stderr,
                "  ribbon output v-rebase: %zu run(s) dropped to their own "
                "v minimum (sidecar %s)\n", atlas_runs, runs_path);
        }
        free(run_vmin); free(run_c0); free(run_c1);
    }

    if (vertex_island == NULL && vertex_chart != NULL &&
        R->source_chart_island != NULL) {
        if (R->source_chart_island_count == 0 ||
            nvout > SIZE_MAX / sizeof(*vertex_island))
            goto done;
        vertex_island = (int32_t *)malloc(nvout * sizeof(*vertex_island));
        if (vertex_island == NULL) goto done;
        for (size_t v = 0; v < nvout; v++) {
            int32_t source_chart = vertex_chart[v];
            if (source_chart < 0) vertex_island[v] = -1;
            else if ((size_t)source_chart >= R->source_chart_island_count ||
                     R->source_chart_island[source_chart] < 0)
                goto done;
            else
                vertex_island[v] = R->source_chart_island[source_chart];
        }
    }
    if (vertex_material == NULL && vertex_chart != NULL &&
        R->source_chart_island != NULL) {
        if (R->source_chart_island_count == 0 ||
            nvout > SIZE_MAX / sizeof(*vertex_material))
            goto done;
        vertex_material =
            (int32_t *)malloc(nvout * sizeof(*vertex_material));
        if (vertex_material == NULL) goto done;
        for (size_t v = 0; v < nvout; v++) {
            int32_t source_chart = vertex_chart[v];
            if (source_chart < 0) vertex_material[v] = -1;
            else if ((size_t)source_chart >= R->source_chart_island_count ||
                     R->source_chart_island[source_chart] < 0)
                goto done;
            else
                vertex_material[v] = R->source_chart_island[source_chart];
        }
    }

    /* Commit all required companions before publishing the authoritative
     * VMESH.  A failed confidence/provenance write therefore cannot leave a
     * fresh VMESH paired with stale authority metadata. */
    if (vertex_phase != NULL) {
        for (size_t v = 0; v < nvout; v++) {
            if (!isfinite((double)vertex_phase[v])) {
                fprintf(stderr,
                        "BUG: non-finite lifted phase at emitted vertex %zu\n",
                        v);
                goto done;
            }
        }
    }
    if (ribbon_write_support(path, vertex_support, nvout) != 0) {
        fprintf(stderr,
                "ribbon: cannot write confidence companion for %s\n", path);
        goto done;
    }
    if (ribbon_write_phase(path, vertex_phase, nvout) != 0) {
        fprintf(stderr,
                "ribbon: cannot write lifted-phase companion for %s\n", path);
        goto done;
    }
    if (ribbon_write_claimant_charts(path, vertex_chart, nvout) != 0) {
        fprintf(stderr,
                "ribbon: cannot write claimant-provenance companion for %s\n",
                path);
        goto done;
    }
    if (ribbon_write_reconstruction_components(
            path, vertex_island, nvout) != 0) {
        fprintf(stderr,
                "ribbon: cannot write reconstruction-component companion for %s\n",
                path);
        goto done;
    }
    if (ribbon_write_material_identities(
            path, vertex_material, nvout) != 0) {
        fprintf(stderr,
                "ribbon: cannot write material-lineage companion for %s\n",
                path);
        goto done;
    }
    if (ribbon_write_claimant_islands(path, vertex_island, nvout) != 0) {
        fprintf(stderr,
                "ribbon: cannot write legacy identity alias for %s\n",
                path);
        goto done;
    }
    {
        char binary_path[2048];
        stage_t0 = ves_clock_sec();
        if (MeshBin_companion_path(path, binary_path, sizeof binary_path) != 0 ||
            strcmp(binary_path, path) == 0 ||
            MeshBin_write(binary_path, verts, nvout, faces, nfout, uv) != 0) {
            fprintf(stderr,
                    "ribbon: cannot write authoritative binary companion for %s\n",
                    path);
            goto done;
        }
        fprintf(stderr,
                "  ribbon output authoritative VMESH: %.2fs\n",
                ves_clock_sec()-stage_t0);
    }
    if (write_obj) {
        stage_t0 = ves_clock_sec();
        ff = fopen(path,"wb");
        if (!ff) goto done;
        fprintf(ff,"# scroll_ribbon fitted grid: %zu source cols -> %zu emitted-atlas cols "
                   "in %zu occupied run(s), %zu rows (v)\n",
                nu,atlas_cols,atlas_runs,nk);
        fprintf(ff,"# faces are the direct fixed-diagonal triangulation of defined grid quads\n");
        for (size_t v = 0; v < nvout; v++) {
            fprintf(ff,"v %.6f %.6f %.6f\n",(double)verts[v*3],
                    (double)verts[v*3+1],(double)verts[v*3+2]);
            fprintf(ff,"vt %.6f %.6f\n",(double)uv[v*2],(double)uv[v*2+1]);
        }
        for (size_t f = 0; f < nfout; f++) {
            int a=faces[f*3]+1,b=faces[f*3+1]+1,c=faces[f*3+2]+1;
            fprintf(ff,"f %d/%d %d/%d %d/%d\n",a,a,b,b,c,c);
        }
        if (fclose(ff) != 0) { ff = NULL; goto done; }
        ff = NULL;
        fprintf(stderr,
                "  ribbon output companion OBJ: %.2fs (total output %.2fs)\n",
                ves_clock_sec()-stage_t0,ves_clock_sec()-total_t0);
    } else {
        fprintf(stderr,
                "  ribbon output: vmesh companion only (total output %.2fs)\n",
                ves_clock_sec()-total_t0);
    }
    *out_nv=nvout; *out_nf=nfout;
    for (size_t v = 0; v < nvout; v++)
        *out_supported += vertex_support[v] == 255;
    *out_generated = nvout - *out_supported;
    *out_atlas_cols=atlas_cols;
    *out_atlas_runs=atlas_runs;
    *out_empty_cols_removed=empty_cols_removed;
    rc=0;
done:
    if (ff) fclose(ff);
    free(run_of_col);
    free(column_map); free(column_used);
    free(vertex_phase); free(vertex_support); free(vertex_island);
    free(vertex_material); free(vertex_chart);
    free(faces); free(uv); free(verts); free(vid);
    return rc;
}

typedef struct {
    size_t layer;
    size_t first_col;
    size_t last_col;
    size_t vertices;
    int promoted;
} RibbonPeelRun;

typedef struct {
    size_t layer;
    size_t first_col;
    size_t last_col;
    size_t dest_col;
} RibbonAtlasSegment;

typedef struct {
    RibbonResult view;
    float *pos;
    uint8_t *present;
    uint8_t *valid;
    float *phi;
    int32_t *island;
    int32_t *material;
    uint8_t *quad;
    size_t *chart_offsets;
    RibbonChartRun *chart_runs;
} RibbonOwnedGrid;

static void ribbon_owned_grid_dispose(RibbonOwnedGrid *owned)
{
    if (owned == NULL) return;
    free(owned->chart_runs);
    free(owned->chart_offsets);
    free(owned->quad);
    free(owned->material);
    free(owned->island);
    free(owned->phi);
    free(owned->valid);
    free(owned->present);
    free(owned->pos);
    memset(owned, 0, sizeof(*owned));
}

static int ribbon_collect_peel_runs(const RibbonResult *source, size_t rows,
                                    RibbonPeelRun **out_runs,
                                    size_t *out_count,
                                    size_t *out_primary_vertices)
{
    RibbonPeelRun *runs = NULL;
    size_t capacity, count = 0, primary_vertices = 0;
    if (out_runs == NULL || out_count == NULL ||
        out_primary_vertices == NULL || source == NULL || rows == 0 ||
        source->grid_layers < 2 || source->nu == 0 ||
        source->grid_present == NULL)
        return -1;
    *out_runs = NULL;
    *out_count = 0;
    *out_primary_vertices = 0;
    if (source->grid_layers - 1 > SIZE_MAX / source->nu) return -1;
    capacity = (source->grid_layers - 1) * source->nu;
    if (capacity > SIZE_MAX / sizeof(*runs)) return -1;
    runs = (RibbonPeelRun *)calloc(capacity ? capacity : 1, sizeof(*runs));
    if (runs == NULL) return -1;

    for (size_t k = 0; k < rows; k++)
        for (size_t j = 0; j < source->nu; j++)
            primary_vertices += source->grid_present[k * source->nu + j] != 0;

    for (size_t layer = 1; layer < source->grid_layers; layer++) {
        int inside = 0;
        for (size_t j = 0; j < source->nu; j++) {
            size_t column_vertices = 0;
            for (size_t k = 0; k < rows; k++) {
                size_t row = layer * rows + k;
                column_vertices +=
                    source->grid_present[row * source->nu + j] != 0;
            }
            if (column_vertices == 0) {
                inside = 0;
                continue;
            }
            if (!inside) {
                if (count >= capacity) {
                    free(runs);
                    return -1;
                }
                runs[count].layer = layer;
                runs[count].first_col = j;
                runs[count].last_col = j;
                count++;
                inside = 1;
            }
            runs[count - 1].last_col = j;
            runs[count - 1].vertices += column_vertices;
        }
    }
    *out_runs = runs;
    *out_count = count;
    *out_primary_vertices = primary_vertices;
    return 0;
}

static int ribbon_select_peel_runs(RibbonPeelRun *runs, size_t run_count,
                                   size_t primary_vertices,
                                   double min_primary_share,
                                   size_t min_vertices,
                                   size_t *out_threshold,
                                   size_t *out_promoted_runs,
                                   size_t *out_promoted_vertices,
                                   size_t *out_remaining_runs,
                                   size_t *out_remaining_vertices)
{
    size_t threshold = 0, promoted_runs = 0, promoted_vertices = 0;
    size_t remaining_runs = 0, remaining_vertices = 0;
    if ((runs == NULL && run_count > 0) || out_threshold == NULL ||
        out_promoted_runs == NULL || out_promoted_vertices == NULL ||
        out_remaining_runs == NULL || out_remaining_vertices == NULL ||
        !isfinite(min_primary_share) || min_primary_share < 0.0)
        return -1;
    if (min_primary_share > 0.0 && primary_vertices > 0) {
        long double wanted = ceill((long double)min_primary_share *
                                   (long double)primary_vertices);
        if (wanted > (long double)SIZE_MAX) return -1;
        threshold = (size_t)wanted;
        if (threshold == 0) threshold = 1;
    }
    if (min_vertices > 0 &&
        (threshold == 0 || min_vertices < threshold))
        threshold = min_vertices;
    for (size_t r = 0; r < run_count; r++) {
        runs[r].promoted = threshold > 0 && runs[r].vertices >= threshold;
        if (runs[r].promoted) {
            promoted_runs++;
            promoted_vertices += runs[r].vertices;
        } else {
            remaining_runs++;
            remaining_vertices += runs[r].vertices;
        }
    }
    *out_threshold = threshold;
    *out_promoted_runs = promoted_runs;
    *out_promoted_vertices = promoted_vertices;
    *out_remaining_runs = remaining_runs;
    *out_remaining_vertices = remaining_vertices;
    return 0;
}

static int ribbon_alloc_grid_arrays(const RibbonResult *source,
                                    size_t slots, size_t cells,
                                    RibbonOwnedGrid *owned)
{
    if (slots == 0 || slots > SIZE_MAX / (3 * sizeof(*owned->pos)))
        return -1;
    owned->pos = (float *)calloc(slots * 3, sizeof(*owned->pos));
    owned->present = (uint8_t *)calloc(slots, sizeof(*owned->present));
    if (source->grid_valid != NULL)
        owned->valid = (uint8_t *)calloc(slots, sizeof(*owned->valid));
    if (source->grid_phi != NULL) {
        if (slots > SIZE_MAX / sizeof(*owned->phi)) return -1;
        owned->phi = (float *)calloc(slots, sizeof(*owned->phi));
    }
    if (source->grid_island != NULL) {
        if (slots > SIZE_MAX / sizeof(*owned->island)) return -1;
        owned->island = (int32_t *)calloc(slots, sizeof(*owned->island));
    }
    if (source->grid_material != NULL) {
        if (slots > SIZE_MAX / sizeof(*owned->material)) return -1;
        owned->material =
            (int32_t *)calloc(slots, sizeof(*owned->material));
    }
    owned->quad = (uint8_t *)calloc(cells ? cells : 1, sizeof(*owned->quad));
    if (owned->pos == NULL || owned->present == NULL || owned->quad == NULL ||
        (source->grid_valid != NULL && owned->valid == NULL) ||
        (source->grid_phi != NULL && owned->phi == NULL) ||
        (source->grid_island != NULL && owned->island == NULL) ||
        (source->grid_material != NULL && owned->material == NULL))
        return -1;
    return 0;
}

/* Re-layout selected peel bands beside one another in U.  A blank source-grid
 * column separates every segment, so the ordinary writer's compactor will
 * retain a two-column atlas gutter and can never synthesize a cross-segment
 * face.  Source topology and every per-vertex authority channel are copied
 * verbatim. */
static int ribbon_build_atlas_view(const RibbonResult *source, size_t rows,
                                   RibbonAtlasSegment *segments,
                                   size_t segment_count,
                                   RibbonOwnedGrid *owned)
{
    size_t dest_nu = 0, slots, cells;
    RibbonResult *dest;
    if (source == NULL || segments == NULL || segment_count == 0 ||
        owned == NULL || rows < 2 || source->nu < 2)
        return -1;
    memset(owned, 0, sizeof(*owned));
    dest = &owned->view;
    *dest = *source;
    for (size_t s = 0; s < segment_count; s++) {
        size_t width;
        if (segments[s].layer >= source->grid_layers ||
            segments[s].first_col > segments[s].last_col ||
            segments[s].last_col >= source->nu)
            goto fail;
        width = segments[s].last_col - segments[s].first_col + 1;
        if (s > 0) {
            if (dest_nu == SIZE_MAX) goto fail;
            dest_nu++;
        }
        segments[s].dest_col = dest_nu;
        if (width > SIZE_MAX - dest_nu) goto fail;
        dest_nu += width;
    }
    if (dest_nu < 2 || dest_nu > SIZE_MAX / rows) goto fail;
    slots = dest_nu * rows;
    if (dest_nu - 1 > SIZE_MAX / (rows - 1)) goto fail;
    cells = (dest_nu - 1) * (rows - 1);
    if (ribbon_alloc_grid_arrays(source, slots, cells, owned) != 0)
        goto fail;

    for (size_t s = 0; s < segment_count; s++) {
        const RibbonAtlasSegment *segment = &segments[s];
        size_t width = segment->last_col - segment->first_col + 1;
        for (size_t k = 0; k < rows; k++) {
            size_t source_row = segment->layer * rows + k;
            for (size_t d = 0; d < width; d++) {
                size_t sq = source_row * source->nu +
                            segment->first_col + d;
                size_t dq = k * dest_nu + segment->dest_col + d;
                memcpy(&owned->pos[dq * 3], &source->grid_pos[sq * 3],
                       3 * sizeof(*owned->pos));
                owned->present[dq] = source->grid_present[sq];
                if (owned->valid != NULL)
                    owned->valid[dq] = source->grid_valid[sq];
                if (owned->phi != NULL)
                    owned->phi[dq] = source->grid_phi[sq];
                if (owned->island != NULL)
                    owned->island[dq] = source->grid_island[sq];
                if (owned->material != NULL)
                    owned->material[dq] = source->grid_material[sq];
            }
        }
        if (source->grid_quad == NULL) goto fail;
        for (size_t k = 0; k + 1 < rows; k++) {
            size_t source_row = segment->layer * rows + k;
            for (size_t d = 0; d + 1 < width; d++) {
                size_t sj = segment->first_col + d;
                size_t dj = segment->dest_col + d;
                owned->quad[k * (dest_nu - 1) + dj] =
                    source->grid_quad[source_row * (source->nu - 1) + sj];
            }
        }
    }

    if ((source->grid_chart_row_offsets == NULL) !=
        (source->grid_chart_runs == NULL))
        goto fail;
    if (source->grid_chart_row_offsets != NULL) {
        size_t emitted_runs = 0;
        if (rows + 1 > SIZE_MAX / sizeof(*owned->chart_offsets) ||
            source->grid_chart_run_count >
                SIZE_MAX / sizeof(*owned->chart_runs))
            goto fail;
        owned->chart_offsets =
            (size_t *)calloc(rows + 1, sizeof(*owned->chart_offsets));
        owned->chart_runs = (RibbonChartRun *)calloc(
            source->grid_chart_run_count ? source->grid_chart_run_count : 1,
            sizeof(*owned->chart_runs));
        if (owned->chart_offsets == NULL || owned->chart_runs == NULL)
            goto fail;
        for (size_t k = 0; k < rows; k++) {
            owned->chart_offsets[k] = emitted_runs;
            for (size_t s = 0; s < segment_count; s++) {
                const RibbonAtlasSegment *segment = &segments[s];
                size_t source_row = segment->layer * rows + k;
                size_t r0 = source->grid_chart_row_offsets[source_row];
                size_t r1 = source->grid_chart_row_offsets[source_row + 1];
                if (r0 > r1 || r1 > source->grid_chart_run_count) goto fail;
                for (size_t r = r0; r < r1; r++) {
                    const RibbonChartRun *in = &source->grid_chart_runs[r];
                    RibbonChartRun *out;
                    size_t first, last;
                    if (in->first_col > in->last_col ||
                        in->last_col >= source->nu)
                        goto fail;
                    first = in->first_col > segment->first_col
                          ? in->first_col : segment->first_col;
                    last = in->last_col < segment->last_col
                         ? in->last_col : segment->last_col;
                    if (first > last) continue;
                    if (emitted_runs >= source->grid_chart_run_count)
                        goto fail;
                    out = &owned->chart_runs[emitted_runs++];
                    out->first_col = segment->dest_col +
                                     first - segment->first_col;
                    out->last_col = segment->dest_col +
                                    last - segment->first_col;
                    out->source_chart = in->source_chart;
                }
            }
        }
        owned->chart_offsets[rows] = emitted_runs;
        owned->view.grid_chart_row_offsets = owned->chart_offsets;
        owned->view.grid_chart_runs = owned->chart_runs;
        owned->view.grid_chart_run_count = emitted_runs;
    }

    dest->grid_pos = owned->pos;
    dest->grid_present = owned->present;
    dest->grid_valid = owned->valid;
    dest->grid_phi = owned->phi;
    dest->grid_island = owned->island;
    dest->grid_material = owned->material;
    dest->grid_quad = owned->quad;
    dest->grid_cut_h = NULL;
    dest->grid_cut_v = NULL;
    dest->nu = dest_nu;
    dest->nk = rows;
    dest->grid_layers = 1;
    if (source->grid_chart_row_offsets == NULL) {
        dest->grid_chart_row_offsets = NULL;
        dest->grid_chart_runs = NULL;
        dest->grid_chart_run_count = 0;
    } else {
        dest->grid_chart_row_offsets = owned->chart_offsets;
        dest->grid_chart_runs = owned->chart_runs;
    }
    return 0;

fail:
    ribbon_owned_grid_dispose(owned);
    return -1;
}

int Ribbon_write_obj_promote_peels(const char *path,
                                   const RibbonResult *result,
                                   int write_obj, int v_rebase,
                                   double min_primary_share,
                                   size_t min_vertices,
                                   RibbonWriteStats *stats)
{
    RibbonWriteStats local;
    char single_inventory_path[2048];
    memset(&local, 0, sizeof local);
    int rc = 0;

    if (!isfinite(min_primary_share) || min_primary_share < 0.0)
        return -1;
    if (result != NULL && result->grid_projective &&
        (min_primary_share > 0.0 || min_vertices > 0)) {
        fprintf(stderr,
                "  ribbon projective atlas: peel promotion disabled; "
                "secondary claims remain in extras\n");
        min_primary_share = 0.0;
        min_vertices = 0;
    }

    /* Peel bands beyond layer 0 are RECOVERED CLAIMS, not the sheet: their
     * geometry is real but their atlas position is a parking spot (the same
     * u as the winning wrap, one band down in v).  Shipping them inside the
     * deliverable's v-space buries the sheet under a second strip of
     * interleaved fragments and corrupts every whole-atlas statistic
     * (2026-09-01: 11.5%% of the 4x5x5 champion's vertices sat in band 1,
     * 100%% of them directly over occupied band-0 cells).  The normal policy
     * therefore writes band 0 to the primary output and the remaining bands
     * to "<stem>_extras.*".  When explicitly armed, a substantial contiguous
     * peel run is instead packed BESIDE band 0 in U.  This recovers genuine
     * coverage without overlaying two claims in one UV cell; the unpromoted
     * partition remains in extras with every authority sidecar intact. */
    if (result != NULL && result->grid_layers > 1 &&
        result->nk >= result->grid_layers &&
        result->nk % result->grid_layers == 0) {
        size_t np = result->nk / result->grid_layers;
        size_t nu = result->nu;
        RibbonResult primary = *result;
        RibbonResult extras = *result;
        const RibbonResult *published = &primary;
        const RibbonResult *remaining = &extras;
        RibbonOwnedGrid published_owned, remaining_owned;
        RibbonPeelRun *runs = NULL;
        RibbonAtlasSegment *published_segments = NULL;
        RibbonAtlasSegment *remaining_segments = NULL;
        RibbonLayerCount *layer_counts = NULL;
        char inventory_path[2048];
        size_t run_count = 0, primary_vertices = 0, threshold = 0;
        size_t promoted_runs = 0, promoted_vertices = 0;
        size_t remaining_runs = 0, remaining_vertices = 0;
        int have_extras = 0;
        rc = -1;
        memset(&published_owned, 0, sizeof published_owned);
        memset(&remaining_owned, 0, sizeof remaining_owned);
        primary.nk = np;
        primary.grid_layers = 1;
        if (primary.grid_chart_row_offsets != NULL)
            primary.grid_chart_run_count =
                primary.grid_chart_row_offsets[np];
        extras.nk = result->nk - np;
        extras.grid_layers = 1;
        extras.grid_pos = result->grid_pos + np * nu * 3;
        extras.grid_present = result->grid_present + np * nu;
        if (result->grid_phi != NULL)
            extras.grid_phi = result->grid_phi + np * nu;
        if (result->grid_valid != NULL)
            extras.grid_valid = result->grid_valid + np * nu;
        if (result->grid_island != NULL)
            extras.grid_island = result->grid_island + np * nu;
        if (result->grid_material != NULL)
            extras.grid_material = result->grid_material + np * nu;
        if (result->grid_quad != NULL)
            extras.grid_quad = result->grid_quad + np * (nu - 1);
            /* grid_quad is a CELL array: row stride nu-1, not nu */
        if (result->grid_chart_row_offsets != NULL)
            extras.grid_chart_row_offsets =
                result->grid_chart_row_offsets + np;

        if (ribbon_collect_peel_runs(result, np, &runs, &run_count,
                                     &primary_vertices) != 0)
            goto layered_done;
        if (ribbon_select_peel_runs(
                runs, run_count, primary_vertices, min_primary_share,
                min_vertices,
                &threshold, &promoted_runs, &promoted_vertices,
                &remaining_runs, &remaining_vertices) != 0)
            goto layered_done;
        have_extras = remaining_vertices > 0;
        local.promoted_peel_runs = promoted_runs;
        local.promoted_peel_vertices = promoted_vertices;
        local.extra_peel_runs = remaining_runs;
        local.extra_peel_vertices = remaining_vertices;
        local.primary_peel_vertices = primary_vertices;
        local.peel_promotion_threshold = threshold;

        /* A projective peel is an additional coordinate chart, not a strip
         * parked below the current crop.  Concatenating layers along V encoded
         * layer 2 at `V + rows*dv`; enlarging A to B then moved every sample in
         * that layer even though its evidence was unchanged.  Publish one
         * VMESH per layer instead.  Every file uses the same absolute U,V
         * lattice; the filename is only a storage color, not a coordinate.
         * Layer 1 retains the historical `_extras` name for consumers which
         * only know about one secondary chart; later layers are numbered. */
        if (result->grid_projective && promoted_runs == 0) {
            if (result->grid_layers > SIZE_MAX / sizeof *layer_counts ||
                ribbon_inventory_path(path, inventory_path, sizeof inventory_path) != 0)
                goto layered_done;
            /* Invalidate an old completion record before replacing any layer. */
            if (remove(inventory_path) != 0 && errno != ENOENT) goto layered_done;
            layer_counts = (RibbonLayerCount *)calloc(result->grid_layers, sizeof *layer_counts);
            if (layer_counts == NULL) goto layered_done;
            rc = write_ribbon_obj(
                path, &primary, write_obj, v_rebase,
                &local.vertices, &local.faces,
                &local.supported_vertices, &local.generated_vertices,
                &local.atlas_columns,
                &local.atlas_runs, &local.empty_columns_removed);
            if (rc != 0) goto layered_done;
            layer_counts[0].vertices = local.vertices;
            layer_counts[0].faces = local.faces;
            for (size_t layer = 1; layer < result->grid_layers; layer++) {
                RibbonResult band = *result;
                size_t occupied = 0, triangles = 0;
                char xpath[2048];
                size_t xnv = 0, xnf = 0, xs = 0, xg = 0, xc = 0, xr = 0,
                       xe = 0;
                if (ribbon_projective_layer_path(
                        path, layer, xpath, sizeof xpath) != 0) {
                    rc = -1;
                    goto layered_done;
                }
                for (size_t k = 0; k < np; k++)
                    for (size_t j = 0; j < nu; j++)
                        occupied += result->grid_present[
                            (layer * np + k) * nu + j] != 0;
                if (result->grid_quad != NULL)
                    for (size_t k = 0; k + 1 < np; k++)
                        for (size_t j = 0; j + 1 < nu; j++) {
                            uint8_t code = result->grid_quad[
                                (layer * np + k) * (nu - 1) + j];
                            for (int bit = 0; bit < 4; bit++)
                                triangles += (code & (uint8_t)(1u << bit)) != 0;
                        }
                /* VMESH is a surface format and deliberately rejects a
                 * vertex-only artifact.  A color can be required by a
                 * sub-cell claimant which never owns a complete triangle;
                 * that diagnostic band has no bakeable output.  Invalidate
                 * an older family at the same exact layer path either way. */
                if (occupied == 0 || triangles == 0) {
                    if (occupied == 0 && triangles != 0) { rc = -1; goto layered_done; }
                    layer_counts[layer].point_only_vertices = occupied;
                    if (ribbon_remove_output_family(xpath) != 0) {
                        fprintf(stderr,
                                "ribbon: cannot invalidate stale point-only "
                                "projective layer %zu at %s\n", layer, xpath);
                        rc = -1;
                        goto layered_done;
                    }
                    if (ribbon_write_point_evidence(xpath, result, layer, np, occupied) != 0) {
                        rc = -1;
                        goto layered_done;
                    }
                    fprintf(stderr,
                            "  ribbon projective peel layer %zu: preserved "
                            "%zu point-only record(s) separately (no emitted surface)\n",
                            layer, occupied);
                    continue;
                }
                band.nk = np;
                band.grid_layers = 1;
                band.grid_pos = result->grid_pos + layer * np * nu * 3;
                band.grid_present = result->grid_present + layer * np * nu;
                if (result->grid_phi != NULL)
                    band.grid_phi = result->grid_phi + layer * np * nu;
                if (result->grid_valid != NULL)
                    band.grid_valid = result->grid_valid + layer * np * nu;
                if (result->grid_island != NULL)
                    band.grid_island = result->grid_island + layer * np * nu;
                if (result->grid_material != NULL)
                    band.grid_material = result->grid_material + layer * np * nu;
                if (result->grid_quad != NULL)
                    band.grid_quad = result->grid_quad +
                                     layer * np * (nu - 1);
                if (result->grid_chart_row_offsets != NULL) {
                    band.grid_chart_row_offsets =
                        result->grid_chart_row_offsets + layer * np;
                    band.grid_chart_run_count =
                        band.grid_chart_row_offsets[np];
                }
                if (write_ribbon_obj(
                        xpath, &band, 0, v_rebase,
                        &xnv, &xnf, &xs, &xg, &xc, &xr, &xe) != 0) {
                    fprintf(stderr,
                            "ribbon: projective peel layer %zu write failed "
                            "for %s\n", layer, xpath);
                    rc = -1;
                    goto layered_done;
                }
                layer_counts[layer].vertices = xnv;
                layer_counts[layer].faces = xnf;
                fprintf(stderr,
                        "  ribbon projective peel layer %zu: %zu verts, "
                        "%zu faces -> %s\n", layer, xnv, xnf, xpath);
            }
            rc = ribbon_write_inventory(inventory_path, layer_counts, result->grid_layers,
                                         result->grid_layer_base, result->grid_source_layers);
            goto layered_done;
        }

        if (promoted_runs > 0) {
            size_t ps = 1, xs = 0;
            if (promoted_runs == SIZE_MAX ||
                promoted_runs + 1 >
                    SIZE_MAX / sizeof(*published_segments) ||
                remaining_runs >
                    SIZE_MAX / sizeof(*remaining_segments))
                goto layered_done;
            published_segments = (RibbonAtlasSegment *)calloc(
                promoted_runs + 1, sizeof(*published_segments));
            remaining_segments = (RibbonAtlasSegment *)calloc(
                remaining_runs ? remaining_runs : 1,
                sizeof(*remaining_segments));
            if (published_segments == NULL || remaining_segments == NULL)
                goto layered_done;
            published_segments[0].layer = 0;
            published_segments[0].first_col = 0;
            published_segments[0].last_col = nu - 1;
            for (size_t r = 0; r < run_count; r++) {
                RibbonAtlasSegment segment;
                memset(&segment, 0, sizeof segment);
                segment.layer = runs[r].layer;
                segment.first_col = runs[r].first_col;
                segment.last_col = runs[r].last_col;
                if (runs[r].promoted)
                    published_segments[ps++] = segment;
                else
                    remaining_segments[xs++] = segment;
            }
            if (ps != promoted_runs + 1 || xs != remaining_runs ||
                ribbon_build_atlas_view(result, np, published_segments, ps,
                                        &published_owned) != 0)
                goto layered_done;
            published = &published_owned.view;
            if (remaining_runs > 0) {
                if (ribbon_build_atlas_view(
                        result, np, remaining_segments, remaining_runs,
                        &remaining_owned) != 0)
                    goto layered_done;
                remaining = &remaining_owned.view;
            }
            fprintf(stderr,
                    "  ribbon peel promotion: %zu/%zu run(s), %zu vertices "
                    "published beside the primary (minimum %zu; %.4g%% of "
                    "%zu primary vertices OR %zu absolute)\n",
                    promoted_runs, run_count, promoted_vertices, threshold,
                    min_primary_share * 100.0, primary_vertices, min_vertices);
            for (size_t r = 0; r < run_count; r++)
                fprintf(stderr,
                        "    peel layer %zu columns [%zu,%zu]: %zu vertices "
                        "-> %s\n",
                        runs[r].layer, runs[r].first_col,
                        runs[r].last_col, runs[r].vertices,
                        runs[r].promoted ? "PRIMARY" : "extras");
        }

        rc = write_ribbon_obj(
            path, published, write_obj, v_rebase,
            &local.vertices, &local.faces,
            &local.supported_vertices, &local.generated_vertices,
            &local.atlas_columns,
            &local.atlas_runs, &local.empty_columns_removed);
        if (rc == 0 && have_extras) {
            char xpath[2048];
            size_t len = strlen(path);
            const char *dot = strrchr(path, '.');
            size_t stem = dot != NULL ? (size_t)(dot - path) : len;
            if (stem + 12 < sizeof xpath) {
                size_t xnv = 0, xnf = 0, xs = 0, xg = 0, xc = 0, xr = 0,
                       xe = 0;
                memcpy(xpath, path, stem);
                memcpy(xpath + stem, "_extras", 7);
                memcpy(xpath + stem + 7, path + stem, len - stem + 1);
                if (write_ribbon_obj(
                        xpath, remaining, 0, v_rebase,
                        &xnv, &xnf, &xs, &xg, &xc, &xr, &xe) != 0) {
                    fprintf(stderr,
                            "ribbon: extras band write failed for %s\n",
                            xpath);
                    rc = -1;
                } else {
                    fprintf(stderr,
                            "  ribbon extras band: %zu verts, %zu faces -> "
                            "%s (unpromoted recovered claims)\n",
                            xnv, xnf, xpath);
                }
            } else rc = -1;
        }
layered_done:
        free(layer_counts);
        free(remaining_segments);
        free(published_segments);
        free(runs);
        ribbon_owned_grid_dispose(&remaining_owned);
        ribbon_owned_grid_dispose(&published_owned);
        if (stats != NULL) *stats = local;
        return rc;
    }

    if (result != NULL && result->grid_projective) {
        if (ribbon_inventory_path(path, single_inventory_path, sizeof single_inventory_path) != 0 ||
            (remove(single_inventory_path) != 0 && errno != ENOENT)) return -1;
        size_t occupied = 0, triangles = 0;
        if (result->grid_present == NULL || result->nu == 0 || result->nk == 0 ||
            result->nk > SIZE_MAX / result->nu) return -1;
        for (size_t q = 0; q < result->nk * result->nu; q++)
            occupied += result->grid_present[q] != 0;
        if (result->grid_quad != NULL)
            for (size_t q = 0; q < (result->nk - 1) * (result->nu - 1); q++)
                for (int bit = 0; bit < 4; bit++)
                    triangles += (result->grid_quad[q] & (uint8_t)(1u << bit)) != 0;
        if (occupied == 0 || (result->grid_quad != NULL && triangles == 0)) {
            const RibbonLayerCount count = {0, 0, occupied};
            if (occupied == 0 && triangles != 0) return -1;
            if (ribbon_remove_output_family(path) != 0) return -1;
            if (ribbon_write_point_evidence(path, result, 0, result->nk, occupied) != 0) return -1;
            rc = ribbon_write_inventory(single_inventory_path, &count, 1,
                                         result->grid_layer_base, result->grid_source_layers);
            if (stats != NULL) *stats = local;
            fprintf(stderr, "  ribbon selected peel %zu: %zu point-only records preserved, no emitted surface\n",
                    result->grid_layer_base, occupied);
            return rc;
        }
    }
    rc = write_ribbon_obj(
        path, result, write_obj, v_rebase,
        &local.vertices, &local.faces,
        &local.supported_vertices, &local.generated_vertices,
        &local.atlas_columns,
        &local.atlas_runs, &local.empty_columns_removed);
    if (rc == 0 && result != NULL && result->grid_projective) {
        const RibbonLayerCount count = {local.vertices, local.faces, 0};
        rc = ribbon_write_inventory(single_inventory_path, &count, 1,
                                     result->grid_layer_base, result->grid_source_layers);
    }
    if (stats != NULL) *stats = local;
    return rc;
}

int Ribbon_write_obj(const char *path, const RibbonResult *result,
                     int write_obj, int v_rebase, RibbonWriteStats *stats)
{
    return Ribbon_write_obj_promote_peels(
        path, result, write_obj, v_rebase, 0.0, 0, stats);
}

int RibbonGridIO_selftest(void)
{
    enum { ROWS = 3, COLS = 8, LAYERS = 2, SLOTS = 48, CELLS = 35 };
    float pos[SLOTS * 3], phi[SLOTS];
    uint8_t present[SLOTS], valid[SLOTS], quad[CELLS];
    int32_t island[SLOTS], material[SLOTS], chart_island[31];
    size_t chart_offsets[ROWS * LAYERS + 1] = { 0, 1, 2, 3, 5, 7, 9 };
    RibbonChartRun chart_runs[9];
    RibbonResult source;
    RibbonPeelRun *runs = NULL;
    RibbonOwnedGrid published, remaining;
    RibbonAtlasSegment published_segments[2] = {
        { 0, 0, COLS - 1, 0 }, { 1, 4, 7, 0 }
    };
    RibbonAtlasSegment remaining_segments[1] = {
        { 1, 0, 1, 0 }
    };
    size_t run_count = 0, primary_vertices = 0, threshold = 0;
    size_t promoted_runs = 0, promoted_vertices = 0;
    size_t remaining_runs = 0, remaining_vertices = 0;
    size_t published_present = 0, residual_present = 0;
    size_t published_cells = 0, residual_cells = 0;
    int fails = 0;

    memset(&source, 0, sizeof source);
    memset(&published, 0, sizeof published);
    memset(&remaining, 0, sizeof remaining);
    memset(present, 0, sizeof present);
    memset(valid, 0, sizeof valid);
    memset(quad, 0, sizeof quad);
    for (size_t q = 0; q < SLOTS; q++) {
        pos[q * 3] = (float)(q / COLS);
        pos[q * 3 + 1] = (float)(q % COLS);
        pos[q * 3 + 2] = (float)q + 0.25f;
        phi[q] = (float)q * 0.01f;
        island[q] = (int32_t)(1000 + q);
        material[q] = (int32_t)(2000 + q);
    }
    for (size_t q = 0; q < 31; q++) chart_island[q] = (int32_t)q;
    for (size_t k = 0; k < ROWS; k++) {
        for (size_t j = 0; j <= 2; j++)
            present[k * COLS + j] = valid[k * COLS + j] = 1;
        for (size_t j = 0; j <= 1; j++)
            present[(ROWS + k) * COLS + j] =
                valid[(ROWS + k) * COLS + j] = 1;
        for (size_t j = 4; j <= 7; j++)
            present[(ROWS + k) * COLS + j] =
                valid[(ROWS + k) * COLS + j] = 1;
        chart_runs[k].first_col = 0;
        chart_runs[k].last_col = 2;
        chart_runs[k].source_chart = 10;
        chart_runs[3 + k * 2].first_col = 0;
        chart_runs[3 + k * 2].last_col = 1;
        chart_runs[3 + k * 2].source_chart = 20;
        chart_runs[4 + k * 2].first_col = 4;
        chart_runs[4 + k * 2].last_col = 7;
        chart_runs[4 + k * 2].source_chart = 30;
    }
    for (size_t k = 0; k + 1 < ROWS; k++) {
        for (size_t j = 0; j < 2; j++)
            quad[k * (COLS - 1) + j] = RIB_CELL_MAIN;
        quad[(ROWS + k) * (COLS - 1)] = RIB_CELL_MAIN;
        for (size_t j = 4; j < 7; j++)
            quad[(ROWS + k) * (COLS - 1) + j] = RIB_CELL_MAIN;
    }
    source.grid_pos = pos;
    source.grid_present = present;
    source.grid_valid = valid;
    source.grid_phi = phi;
    source.grid_island = island;
    source.grid_material = material;
    source.grid_quad = quad;
    source.grid_chart_row_offsets = chart_offsets;
    source.grid_chart_runs = chart_runs;
    source.grid_chart_run_count = 9;
    source.source_chart_island = chart_island;
    source.source_chart_island_count = 31;
    source.nu = COLS;
    source.nk = ROWS * LAYERS;
    source.grid_layers = LAYERS;
    source.grid_du = 2.0f;
    source.grid_dv = 2.0f;

    {
        RibbonResult point_source = source;
        RibbonPointRecord point;
        point_source.grid_projective = 1;
        point_source.grid_u_origin = -12.0;
        point_source.grid_v_origin = 99.25;
        if (ribbon_point_record(&point_source, 1, ROWS, 1, 0, &point) ||
            point.u != ribbon_f32_bits(-12.0f) ||
            point.v != ribbon_f32_bits(101.25f) ||
            point.z != ribbon_f32_bits(pos[(ROWS + 1) * COLS * 3]) ||
            point.phase != ribbon_f32_bits(phi[(ROWS + 1) * COLS]) ||
            point.chart != 20 || point.material != material[(ROWS + 1) * COLS] ||
            point.reconstruction != island[(ROWS + 1) * COLS] || point.support != 255 ||
            ribbon_point_record(&point_source, 2, ROWS, 0, 0, &point) == 0 ||
            ribbon_point_record(&point_source, 1, ROWS, 0, 3, &point) == 0) {
            fprintf(stderr, "  FAIL gridio point-only evidence bits / physical V / bounds\n");
            fails++;
        }
        valid[(ROWS + 1) * COLS] = 0;
        if (ribbon_point_record(&point_source, 1, ROWS, 1, 0, &point) ||
            point.support != 0 || point.chart != -1) {
            fprintf(stderr, "  FAIL gridio generated point evidence support\n");
            fails++;
        }
        valid[(ROWS + 1) * COLS] = 1;
        point_source.grid_chart_runs = NULL;
        if (ribbon_point_record(&point_source, 1, ROWS, 1, 0, &point) == 0) {
            fprintf(stderr, "  FAIL gridio point evidence missing chart accepted\n");
            fails++;
        }
        fprintf(stderr, "  %s gridio isolated point evidence\n", fails ? "FAIL" : "PASS");
    }

    if (ribbon_collect_peel_runs(&source, ROWS, &runs, &run_count,
                                 &primary_vertices) != 0 ||
        run_count != 2 || primary_vertices != 9 || runs == NULL ||
        runs[0].first_col != 0 || runs[0].last_col != 1 ||
        runs[0].vertices != 6 || runs[1].first_col != 4 ||
        runs[1].last_col != 7 || runs[1].vertices != 12) {
        fprintf(stderr, "  FAIL gridio peel run census\n");
        fails++;
        goto done;
    }
    if (ribbon_select_peel_runs(
            runs, run_count, primary_vertices, 0.75, 0, &threshold,
            &promoted_runs, &promoted_vertices, &remaining_runs,
            &remaining_vertices) != 0 ||
        threshold != 7 || promoted_runs != 1 || promoted_vertices != 12 ||
        remaining_runs != 1 || remaining_vertices != 6 ||
        runs[0].promoted || !runs[1].promoted) {
        fprintf(stderr, "  FAIL gridio peel threshold partition\n");
        fails++;
        goto done;
    }
    if (ribbon_select_peel_runs(
            runs, run_count, primary_vertices, 0.75, 5, &threshold,
            &promoted_runs, &promoted_vertices, &remaining_runs,
            &remaining_vertices) != 0 ||
        threshold != 5 || promoted_runs != 2 || promoted_vertices != 18 ||
        remaining_runs != 0 || remaining_vertices != 0 ||
        !runs[0].promoted || !runs[1].promoted) {
        fprintf(stderr, "  FAIL gridio absolute peel threshold criterion\n");
        fails++;
        goto done;
    }
    if (ribbon_build_atlas_view(&source, ROWS, published_segments, 2,
                                &published) != 0 ||
        ribbon_build_atlas_view(&source, ROWS, remaining_segments, 1,
                                &remaining) != 0) {
        fprintf(stderr, "  FAIL gridio peel horizontal repack\n");
        fails++;
        goto done;
    }
    for (size_t q = 0; q < published.view.nu * published.view.nk; q++)
        published_present += published.view.grid_present[q] != 0;
    for (size_t q = 0; q < remaining.view.nu * remaining.view.nk; q++)
        residual_present += remaining.view.grid_present[q] != 0;
    for (size_t q = 0;
         q < (published.view.nu - 1) * (published.view.nk - 1); q++)
        published_cells += published.view.grid_quad[q] != 0;
    for (size_t q = 0;
         q < (remaining.view.nu - 1) * (remaining.view.nk - 1); q++)
        residual_cells += remaining.view.grid_quad[q] != 0;
    if (published.view.nu != 13 || remaining.view.nu != 2 ||
        published_present != 21 || residual_present != 6 ||
        published_cells != 10 || residual_cells != 2 ||
        published.view.grid_present[8] != 0 ||
        published.view.grid_present[9] == 0 ||
        published.view.grid_quad[7] != 0 ||
        published.view.grid_quad[8] != 0 ||
        published.view.grid_chart_run_count != 6 ||
        remaining.view.grid_chart_run_count != 3 ||
        published.view.grid_chart_row_offsets[3] != 6 ||
        published.view.grid_chart_runs[1].first_col != 9 ||
        published.view.grid_chart_runs[1].last_col != 12 ||
        published.view.grid_chart_runs[1].source_chart != 30 ||
        remaining.view.grid_chart_runs[0].source_chart != 20 ||
        published.view.grid_material[9] != material[ROWS * COLS + 4] ||
        published.view.grid_island[9] != island[ROWS * COLS + 4] ||
        published.view.grid_phi[9] != phi[ROWS * COLS + 4]) {
        fprintf(stderr, "  FAIL gridio peel topology/provenance preservation\n");
        fails++;
    }

done:
    free(runs);
    ribbon_owned_grid_dispose(&remaining);
    ribbon_owned_grid_dispose(&published);
    fprintf(stderr, "  %s gridio peel atlas partition\n",
            fails == 0 ? "PASS" : "FAIL");
    return fails;
}
