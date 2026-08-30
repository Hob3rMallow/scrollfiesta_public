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
            uv[(size_t)v*2] = (float)j * R->grid_du;
            uv[(size_t)v*2+1] = (float)k * R->grid_dv;
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
    /* Face emission is a direct serialization of fitted-grid topology.  A
     * defined quad is two fixed-diagonal triangles; an undefined corner means
     * the parameterizer did not create that quad. */
    for (size_t k = 0; k + 1 < nk; k++) {
        for (size_t j = 0; j + 1 < nu; j++) {
            int32_t a = vid[k*nu + j],     b = vid[k*nu + j + 1];
            int32_t c = vid[(k+1)*nu + j], d = vid[(k+1)*nu + j + 1];
            if (a >= 0 && b >= 0 && c >= 0 && d >= 0) nfout += 2;
        }
    }
    if (nfout == 0 || nfout > SIZE_MAX / (3 * sizeof(*faces))) goto done;
    faces = (int32_t *)malloc(nfout * 3 * sizeof(*faces));
    if (!faces) goto done;
    nfout = 0;
    for (size_t k = 0; k + 1 < nk; k++) {
        for (size_t j = 0; j + 1 < nu; j++) {
            int32_t a = vid[k*nu + j],     b = vid[k*nu + j + 1];
            int32_t c = vid[(k+1)*nu + j], d = vid[(k+1)*nu + j + 1];
            if (a < 0 || b < 0 || c < 0 || d < 0) continue;
            faces[nfout*3] = a;
            faces[nfout*3+1] = b;
            faces[nfout*3+2] = c;
            nfout++;
            faces[nfout*3] = b;
            faces[nfout*3+1] = d;
            faces[nfout*3+2] = c;
            nfout++;
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
        long j = lround((double)uv[v*2] / (double)R->grid_du);
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
                if (atlas_runs > 0) cursor += gutter_columns;
                atlas_runs++;
                inside = 1;
            }
            if (cursor > (size_t)INT32_MAX) goto done;
            run_of_col[cursor] = (int32_t)(atlas_runs - 1);
            column_map[j] = (int32_t)cursor++;
        }
        atlas_cols = cursor;
        empty_cols_removed = nu > atlas_cols ? nu - atlas_cols : 0;
    }
    for (size_t v = 0; v < nvout; v++) {
        long old_col = lround((double)uv[v*2] / (double)R->grid_du);
        int32_t new_col = column_map[(size_t)old_col];
        if (new_col < 0) goto done;
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

int Ribbon_write_obj(const char *path, const RibbonResult *result,
                     int write_obj, int v_rebase, RibbonWriteStats *stats)
{
    RibbonWriteStats local;
    memset(&local, 0, sizeof local);
    int rc = write_ribbon_obj(
        path, result, write_obj, v_rebase,
        &local.vertices, &local.faces,
        &local.supported_vertices, &local.generated_vertices,
        &local.atlas_columns,
        &local.atlas_runs, &local.empty_columns_removed);
    if (stats != NULL) *stats = local;
    return rc;
}
