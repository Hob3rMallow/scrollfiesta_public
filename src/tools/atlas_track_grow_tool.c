/*
 * atlas_track_grow_tool.c
 *
 * Experimental evidence-first atlas outputter.  This is intentionally a
 * separate executable from atlas_ribbon_fit: it consumes the same frozen
 * section evidence, grows only four-edge cycle-supported cells, and writes
 * a new parameterized patch without altering any existing checkpoint.
 */

#include "../common/arena.h"
#include "../common/ves_platform.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"
#include "../whole/atlas_ribbon_fit.h"
#include "../whole/atlas_solution.h"
#include "../whole/atlas_track_export.h"
#include "../whole/atlas_track_grow.h"

#include <ctype.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#include <windows.h>
#else
#include <dirent.h>
#endif

#define ATGT_PATH_CAP 2048

/* Checkpoint 58 predates AtlasSolution v2 by one field: v1 chart records do
 * not carry the final continuous wind.  The general AtlasSolution reader is
 * intentionally strict and must keep rejecting v1.  This experimental tool
 * has a narrower adapter: it validates the same PieceSet fingerprint, reads
 * the v1 payload, and reconstructs chart wind from that checkpoint's own
 * charts.csv (w_phase plus the per-group rounded gauge used by the producer).
 * No state is borrowed from a newer checkpoint. */
static int atgt_read_exact(FILE *fp, void *data, size_t bytes)
{
    return bytes == 0 || fread(data, 1, bytes, fp) == bytes ? 0 : -1;
}

static int atgt_compare_double(const void *left, const void *right)
{
    double a = *(const double *)left;
    double b = *(const double *)right;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static size_t atgt_split_csv(char *line, char **field, size_t capacity)
{
    if (line == NULL || field == NULL || capacity == 0) return 0;
    size_t count = 1;
    field[0] = line;
    for (char *p = line; *p != '\0'; p++) {
        if (*p == '\r' || *p == '\n') {
            *p = '\0';
            break;
        }
        if (*p == ',') {
            *p = '\0';
            if (count >= capacity) return 0;
            field[count++] = p + 1;
        }
    }
    return count;
}

static int atgt_parse_long(const char *text, long *out)
{
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (end == text || *end != '\0') return -1;
    *out = value;
    return 0;
}

static int atgt_parse_double(const char *text, double *out)
{
    char *end = NULL;
    double value = strtod(text, &end);
    if (end == text || *end != '\0' || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

typedef struct {
    int32_t a;
    int32_t b;
} AtgtChartPair;

static int atgt_chart_pair_compare(const void *left, const void *right)
{
    const AtgtChartPair *a = (const AtgtChartPair *)left;
    const AtgtChartPair *b = (const AtgtChartPair *)right;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : (a->b > b->b ? 1 : 0);
}

static int atgt_find_chart_shard(
    const char *dir,
    const char *cube_id,
    char path[ATGT_PATH_CAP])
{
    if (dir == NULL || cube_id == NULL || path == NULL) return -1;
#ifdef _MSC_VER
    char pattern[ATGT_PATH_CAP];
    int n = snprintf(pattern, sizeof pattern, "%s/%s.*.gwshard",
                     dir, cube_id);
    if (n < 0 || (size_t)n >= sizeof pattern) return -1;
    WIN32_FIND_DATAA data;
    HANDLE handle = FindFirstFileA(pattern, &data);
    if (handle == INVALID_HANDLE_VALUE) return -1;
    int found = 0;
    do {
        if (data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
        if (found) {
            FindClose(handle);
            return -1;
        }
        n = snprintf(path, ATGT_PATH_CAP, "%s/%s", dir, data.cFileName);
        if (n < 0 || n >= ATGT_PATH_CAP) {
            FindClose(handle);
            return -1;
        }
        found = 1;
    } while (FindNextFileA(handle, &data));
    FindClose(handle);
    return found ? 0 : -1;
#else
    DIR *directory = opendir(dir);
    if (directory == NULL) return -1;
    size_t prefix = strlen(cube_id);
    const char *suffix = ".gwshard";
    size_t suffix_length = strlen(suffix);
    int found = 0;
    struct dirent *entry;
    while ((entry = readdir(directory)) != NULL) {
        size_t length = strlen(entry->d_name);
        if (length <= prefix + suffix_length + 1 ||
            strncmp(entry->d_name, cube_id, prefix) != 0 ||
            entry->d_name[prefix] != '.' ||
            strcmp(entry->d_name + length - suffix_length, suffix) != 0)
            continue;
        if (found) {
            closedir(directory);
            return -1;
        }
        int n = snprintf(path, ATGT_PATH_CAP, "%s/%s",
                         dir, entry->d_name);
        if (n < 0 || n >= ATGT_PATH_CAP) {
            closedir(directory);
            return -1;
        }
        found = 1;
    }
    closedir(directory);
    return found ? 0 : -1;
#endif
}

static int atgt_map_cube_shard(
    const char *shard_dir,
    size_t cube,
    const PieceSet *ps,
    const AtlasSolution *solution,
    const size_t *node_cube,
    const int32_t *node_local_chart,
    size_t nnode,
    int32_t *chart_node,
    size_t *node_member_count)
{
    static const unsigned char expected_magic[8] = {
        'G', 'W', 'S', 'H', 'R', 'D', '1', '\n'
    };
    char path[ATGT_PATH_CAP], shard_id[48];
    unsigned char magic[8];
    uint32_t version = 0, flags = 0;
    uint64_t nv64 = 0, nf64 = 0, ncharts64 = 0;
    FILE *fp = NULL;
    int32_t *chart = NULL, *parent0 = NULL, *parent1 = NULL;
    int32_t *local_node = NULL;
    int result = -1;
    if (shard_dir == NULL || cube >= ps->n_cubes ||
        node_cube == NULL || node_local_chart == NULL ||
        chart_node == NULL || node_member_count == NULL ||
        atgt_find_chart_shard(shard_dir, ps->ids[cube], path) != 0)
        return -1;
    fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    if (atgt_read_exact(fp, magic, sizeof magic) != 0 ||
        atgt_read_exact(fp, &version, sizeof version) != 0 ||
        atgt_read_exact(fp, &flags, sizeof flags) != 0 ||
        atgt_read_exact(fp, shard_id, sizeof shard_id) != 0 ||
        atgt_read_exact(fp, &nv64, sizeof nv64) != 0 ||
        atgt_read_exact(fp, &nf64, sizeof nf64) != 0 ||
        atgt_read_exact(fp, &ncharts64, sizeof ncharts64) != 0 ||
        memcmp(magic, expected_magic, sizeof magic) != 0 ||
        version != 1 || !(flags & 1) ||
        memchr(shard_id, '\0', sizeof shard_id) == NULL ||
        strcmp(shard_id, ps->ids[cube]) != 0 || nv64 == 0 || nf64 == 0 ||
        nv64 > SIZE_MAX / sizeof(int32_t) ||
        ncharts64 == 0 || ncharts64 > (uint64_t)INT32_MAX ||
        ncharts64 > SIZE_MAX / sizeof(int32_t) ||
        nv64 + nf64 > (uint64_t)LONG_MAX / 12)
        goto cleanup;
    size_t nv = (size_t)nv64;
    size_t ncharts = (size_t)ncharts64;
    local_node = (int32_t *)malloc(ncharts * sizeof(*local_node));
    if (local_node == NULL) goto cleanup;
    for (size_t local = 0; local < ncharts; local++) local_node[local] = -1;
    for (size_t node = 0; node < nnode; node++) {
        if (node_cube[node] != cube) continue;
        int32_t local = node_local_chart[node];
        if (local < 0 || (size_t)local >= ncharts ||
            local_node[local] >= 0)
            goto cleanup;
        local_node[local] = (int32_t)node;
    }

    long skip = (long)((nv64 + nf64) * 12);
    if (fseek(fp, skip, SEEK_CUR) != 0) goto cleanup;
    chart = (int32_t *)malloc(nv * sizeof(*chart));
    parent0 = (int32_t *)malloc(nv * sizeof(*parent0));
    parent1 = (int32_t *)malloc(nv * sizeof(*parent1));
    if (chart == NULL || parent0 == NULL || parent1 == NULL ||
        atgt_read_exact(fp, chart, nv * sizeof(*chart)) != 0 ||
        atgt_read_exact(fp, parent0, nv * sizeof(*parent0)) != 0 ||
        atgt_read_exact(fp, parent1, nv * sizeof(*parent1)) != 0 ||
        fgetc(fp) != EOF)
        goto cleanup;

    size_t cube_vertices = ps->cube_voff[cube + 1] - ps->cube_voff[cube];
    for (size_t v = 0; v < nv; v++) {
        int32_t local = chart[v];
        if (local < -1 || (local >= 0 && (size_t)local >= ncharts))
            goto cleanup;
        if (local < 0 || local_node[local] < 0 ||
            parent0[v] > -2 || parent1[v] < 0)
            continue;
        if ((size_t)parent1[v] >= cube_vertices) goto cleanup;
        int32_t candidate = solution->vertex_chart[
            ps->cube_voff[cube] + (size_t)parent1[v]];
        if (candidate < 0) continue;
        if ((size_t)candidate >= solution->ncharts ||
            solution->chart[candidate].chart != candidate)
            goto cleanup;
        int32_t node = local_node[local];
        if (chart_node[candidate] >= 0 && chart_node[candidate] != node)
            goto cleanup;
        if (chart_node[candidate] < 0) {
            chart_node[candidate] = node;
            node_member_count[node]++;
        }
    }
    result = 0;

cleanup:
    if (fp != NULL) fclose(fp);
    free(chart);
    free(parent0);
    free(parent1);
    free(local_node);
    return result;
}

static int atgt_load_direct_chart_relations(
    Arena_T arena,
    const char *path,
    const char *shard_dir,
    const PieceSet *ps,
    const AtlasSolution *solution,
    int32_t **out_a,
    int32_t **out_b,
    size_t *out_count,
    size_t *out_nodes,
    size_t *out_charts,
    double **out_layout_u)
{
    static const unsigned char expected_magic[8] = {
        'G', 'W', 'L', 'A', 'Y', 'T', '1', '\n'
    };
    unsigned char magic[8];
    uint32_t version = 0, flags = 0;
    uint64_t nnode64 = 0, nedge64 = 0, ncycle64 = 0;
    FILE *fp = NULL;
    size_t *node_cube = NULL, *node_member_count = NULL;
    int32_t *node_local_chart = NULL;
    int32_t *node_representative = NULL;
    int32_t *node_component = NULL, *node_turn = NULL;
    double *node_v_reference = NULL, *node_v_step = NULL;
    double *node_shift_v = NULL;
    int32_t *chart_node = NULL;
    uint8_t *cube_seen = NULL;
    uint32_t *edge_a = NULL, *edge_b = NULL;
    AtgtChartPair *pair = NULL;
    int result = -1;

    if (arena == NULL || path == NULL || ps == NULL || solution == NULL ||
        out_a == NULL || out_b == NULL || out_count == NULL ||
        out_nodes == NULL || out_charts == NULL || out_layout_u == NULL)
        return -1;
    *out_a = NULL;
    *out_b = NULL;
    *out_count = 0;
    *out_nodes = 0;
    *out_charts = 0;
    *out_layout_u = NULL;
    fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    if (atgt_read_exact(fp, magic, sizeof magic) != 0 ||
        atgt_read_exact(fp, &version, sizeof version) != 0 ||
        atgt_read_exact(fp, &flags, sizeof flags) != 0 ||
        atgt_read_exact(fp, &nnode64, sizeof nnode64) != 0 ||
        atgt_read_exact(fp, &nedge64, sizeof nedge64) != 0 ||
        atgt_read_exact(fp, &ncycle64, sizeof ncycle64) != 0 ||
        memcmp(magic, expected_magic, sizeof magic) != 0 ||
        version != 4 || flags != 0 || nnode64 == 0 ||
        nnode64 > ps->nv || nnode64 > SIZE_MAX ||
        nedge64 == 0 || nedge64 > SIZE_MAX || ncycle64 > nedge64)
        goto cleanup;

    size_t nnode = (size_t)nnode64;
    size_t nedge = (size_t)nedge64;
    node_cube = (size_t *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_cube));
    node_local_chart = (int32_t *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_local_chart));
    node_representative = (int32_t *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_representative));
    node_component = (int32_t *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_component));
    node_turn = (int32_t *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_turn));
    node_v_reference = (double *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_v_reference));
    node_v_step = (double *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_v_step));
    node_shift_v = (double *)ARENA_ALLOC(
        arena, nnode * 9 * sizeof(*node_shift_v));
    edge_a = (uint32_t *)ARENA_ALLOC(arena, nedge * sizeof(*edge_a));
    edge_b = (uint32_t *)ARENA_ALLOC(arena, nedge * sizeof(*edge_b));
    cube_seen = (uint8_t *)ARENA_CALLOC(arena, ps->n_cubes, 1);

    size_t last_cube = SIZE_MAX;
    for (size_t i = 0; i < nnode; i++) {
        char cube_id[48];
        int32_t local_chart = -1, representative = -1;
        double v_reference = 0.0, v_step = 0.0, shift_v[9];
        int32_t winding_turn = 0, relation_component = -1;
        if (atgt_read_exact(fp, cube_id, sizeof cube_id) != 0 ||
            atgt_read_exact(fp, &local_chart, sizeof local_chart) != 0 ||
            atgt_read_exact(fp, &representative, sizeof representative) != 0 ||
            atgt_read_exact(fp, &v_reference, sizeof v_reference) != 0 ||
            atgt_read_exact(fp, &v_step, sizeof v_step) != 0 ||
            atgt_read_exact(fp, shift_v, sizeof shift_v) != 0 ||
            atgt_read_exact(fp, &winding_turn, sizeof winding_turn) != 0 ||
            atgt_read_exact(fp, &relation_component,
                            sizeof relation_component) != 0 ||
            memchr(cube_id, '\0', sizeof cube_id) == NULL ||
            local_chart < 0 || representative < 0 ||
            relation_component < 0 ||
            (uint64_t)relation_component >= nnode64 ||
            !isfinite(v_reference) || !isfinite(v_step) || v_step <= 0.0) {
            fprintf(stderr,
                    "[atlas_track_grow] invalid layout node record %zu\n", i);
            goto cleanup;
        }
        for (size_t knot = 0; knot < 9; knot++)
            if (!isfinite(shift_v[knot])) {
                fprintf(stderr,
                        "[atlas_track_grow] invalid layout knot node=%zu "
                        "knot=%zu\n", i, knot);
                goto cleanup;
            }

        size_t cube = ps->n_cubes;
        if (last_cube < ps->n_cubes &&
            strcmp(ps->ids[last_cube], cube_id) == 0) {
            cube = last_cube;
        } else {
            for (size_t c = 0; c < ps->n_cubes; c++)
                if (strcmp(ps->ids[c], cube_id) == 0) {
                    cube = c;
                    break;
                }
            last_cube = cube;
        }
        if (cube >= ps->n_cubes ||
            (size_t)representative >=
                ps->cube_voff[cube + 1] - ps->cube_voff[cube]) {
            fprintf(stderr,
                    "[atlas_track_grow] layout node provenance miss node=%zu "
                    "cube=%s local_chart=%d representative=%d\n",
                    i, cube_id, local_chart, representative);
            goto cleanup;
        }
        node_cube[i] = cube;
        node_local_chart[i] = local_chart;
        node_representative[i] = representative;
        node_component[i] = relation_component;
        node_turn[i] = winding_turn;
        node_v_reference[i] = v_reference;
        node_v_step[i] = v_step;
        memcpy(node_shift_v + i * 9, shift_v, sizeof shift_v);
        cube_seen[cube] = 1;
    }

    for (size_t e = 0; e < nedge; e++) {
        uint32_t na = 0, nb = 0;
        double support = 0.0;
        if (atgt_read_exact(fp, &na, sizeof na) != 0 ||
            atgt_read_exact(fp, &nb, sizeof nb) != 0 ||
            atgt_read_exact(fp, &support, sizeof support) != 0 ||
            na >= nnode || nb >= nnode || na == nb ||
            !isfinite(support) || support <= 0.0 ||
            node_component[na] != node_component[nb] ||
            node_turn[na] != node_turn[nb]) {
            fprintf(stderr,
                    "[atlas_track_grow] invalid layout relation %zu "
                    "nodes=%u/%u\n", e, na, nb);
            goto cleanup;
        }
        edge_a[e] = na;
        edge_b[e] = nb;
    }
    if (fgetc(fp) != EOF) goto cleanup;

    chart_node = (int32_t *)ARENA_ALLOC(
        arena, solution->ncharts * sizeof(*chart_node));
    node_member_count = (size_t *)ARENA_CALLOC(
        arena, nnode, sizeof(*node_member_count));
    for (size_t chart = 0; chart < solution->ncharts; chart++)
        chart_node[chart] = -1;
    for (size_t cube = 0; cube < ps->n_cubes; cube++)
        if (cube_seen[cube] && atgt_map_cube_shard(
                shard_dir, cube, ps, solution,
                node_cube, node_local_chart, nnode,
                chart_node, node_member_count) != 0) {
            fprintf(stderr,
                    "[atlas_track_grow] cube-sheet membership mismatch "
                    "cube=%s\n", ps->ids[cube]);
            goto cleanup;
        }
    for (size_t node = 0; node < nnode; node++) {
        size_t vertex = ps->cube_voff[node_cube[node]] +
                        (size_t)node_representative[node];
        int32_t chart = solution->vertex_chart[vertex];
        if (chart >= 0 && ((size_t)chart >= solution->ncharts ||
            chart_node[chart] != (int32_t)node)) {
            fprintf(stderr,
                    "[atlas_track_grow] representative membership mismatch "
                    "node=%zu chart=%d\n", node, chart);
            goto cleanup;
        }
    }

    size_t mapped_nodes = 0, mapped_charts = 0;
    size_t *node_offset = (size_t *)ARENA_ALLOC(
        arena, (nnode + 1) * sizeof(*node_offset));
    node_offset[0] = 0;
    for (size_t node = 0; node < nnode; node++) {
        if (node_member_count[node] > 0) mapped_nodes++;
        if (node_offset[node] > SIZE_MAX - node_member_count[node])
            goto cleanup;
        node_offset[node + 1] = node_offset[node] + node_member_count[node];
    }
    mapped_charts = node_offset[nnode];
    if (mapped_charts == 0) goto cleanup;
    int32_t *node_member = (int32_t *)ARENA_ALLOC(
        arena, mapped_charts * sizeof(*node_member));
    size_t *node_cursor = (size_t *)ARENA_ALLOC(
        arena, nnode * sizeof(*node_cursor));
    memcpy(node_cursor, node_offset, nnode * sizeof(*node_cursor));
    for (size_t chart = 0; chart < solution->ncharts; chart++) {
        int32_t node = chart_node[chart];
        if (node >= 0)
            node_member[node_cursor[node]++] = (int32_t)chart;
    }

    size_t pair_count = 0;
    for (size_t node = 0; node < nnode; node++) {
        size_t count = node_member_count[node];
        if (count > 1) {
            if (count > SIZE_MAX / (count - 1)) goto cleanup;
            size_t add = count * (count - 1) / 2;
            if (pair_count > SIZE_MAX - add) goto cleanup;
            pair_count += add;
        }
    }
    for (size_t e = 0; e < nedge; e++) {
        size_t ca = node_member_count[edge_a[e]];
        size_t cb = node_member_count[edge_b[e]];
        if (ca != 0 && cb > SIZE_MAX / ca) goto cleanup;
        size_t add = ca * cb;
        if (pair_count > SIZE_MAX - add) goto cleanup;
        pair_count += add;
    }
    if (pair_count == 0 || pair_count > SIZE_MAX / sizeof(*pair))
        goto cleanup;
    pair = (AtgtChartPair *)ARENA_ALLOC(
        arena, pair_count * sizeof(*pair));
    size_t at = 0;
    for (size_t node = 0; node < nnode; node++)
        for (size_t i = node_offset[node]; i < node_offset[node + 1]; i++)
            for (size_t j = i + 1; j < node_offset[node + 1]; j++) {
                pair[at].a = node_member[i];
                pair[at].b = node_member[j];
                at++;
            }
    for (size_t e = 0; e < nedge; e++)
        for (size_t i = node_offset[edge_a[e]];
             i < node_offset[edge_a[e] + 1]; i++)
            for (size_t j = node_offset[edge_b[e]];
                 j < node_offset[edge_b[e] + 1]; j++) {
                int32_t a = node_member[i], b = node_member[j];
                if (b < a) {
                    int32_t swap = a;
                    a = b;
                    b = swap;
                }
                pair[at].a = a;
                pair[at].b = b;
                at++;
            }
    if (at != pair_count) goto cleanup;
    qsort(pair, pair_count, sizeof(*pair), atgt_chart_pair_compare);
    for (size_t e = 1; e < pair_count; e++)
        if (pair[e - 1].a == pair[e].a && pair[e - 1].b == pair[e].b) {
            fprintf(stderr,
                    "[atlas_track_grow] duplicate mapped relation %zu "
                    "charts=%d/%d\n", e, pair[e].a, pair[e].b);
            goto cleanup;
        }

    int32_t *a = (int32_t *)ARENA_ALLOC(
        arena, pair_count * sizeof(*a));
    int32_t *b = (int32_t *)ARENA_ALLOC(
        arena, pair_count * sizeof(*b));
    for (size_t e = 0; e < pair_count; e++) {
        a[e] = pair[e].a;
        b[e] = pair[e].b;
    }
    double *layout_u = (double *)ARENA_ALLOC(
        arena, solution->nvertices * sizeof(*layout_u));
    for (size_t vertex = 0; vertex < solution->nvertices; vertex++) {
        int32_t chart = solution->vertex_chart[vertex];
        int32_t node = chart >= 0 && (size_t)chart < solution->ncharts
                     ? chart_node[chart] : -1;
        if (node < 0) {
            layout_u[vertex] = solution->u[vertex];
            continue;
        }
        double base = ps->uv[vertex * 2];
        double z = ps->verts[vertex * 3];
        double t = (z - node_v_reference[node]) / node_v_step[node];
        double shift;
        if (!isfinite(base) || !isfinite(z) || !isfinite(t)) goto cleanup;
        if (t <= 0.0) {
            shift = node_shift_v[(size_t)node * 9];
        } else if (t >= 8.0) {
            shift = node_shift_v[(size_t)node * 9 + 8];
        } else {
            size_t knot = (size_t)floor(t);
            double fraction = t - (double)knot;
            shift = (1.0 - fraction) *
                        node_shift_v[(size_t)node * 9 + knot] +
                    fraction *
                        node_shift_v[(size_t)node * 9 + knot + 1];
        }
        layout_u[vertex] = base + shift;
        if (!isfinite(layout_u[vertex])) goto cleanup;
    }
    *out_a = a;
    *out_b = b;
    *out_count = pair_count;
    *out_nodes = mapped_nodes;
    *out_charts = mapped_charts;
    *out_layout_u = layout_u;
    result = 0;

cleanup:
    if (fp != NULL) fclose(fp);
    return result;
}

static int atgt_load_v1_chart_wind(
    Arena_T arena,
    const char *solution_path,
    AtlasSolution *solution)
{
    char csv_path[ATGT_PATH_CAP];
    int n = snprintf(csv_path, sizeof csv_path, "%s", solution_path);
    if (n < 0 || (size_t)n >= sizeof csv_path) return -1;
    char *slash = strrchr(csv_path, '/');
    char *backslash = strrchr(csv_path, '\\');
    char *separator = slash;
    if (backslash != NULL &&
        (separator == NULL || backslash > separator))
        separator = backslash;
    if (separator == NULL) return -1;
    int tail = snprintf(separator + 1,
                        sizeof csv_path - (size_t)(separator + 1 - csv_path),
                        "charts.csv");
    if (tail < 0 || (size_t)tail >=
            sizeof csv_path - (size_t)(separator + 1 - csv_path))
        return -1;

    FILE *fp = fopen(csv_path, "rb");
    if (fp == NULL) return -1;
    char line[4096];
    if (fgets(line, sizeof line, fp) == NULL ||
        strstr(line, "w_phi0") == NULL || strstr(line, "w_phase") == NULL) {
        fclose(fp);
        return -1;
    }
    double *w_phi0 = (double *)ARENA_ALLOC(
        arena, solution->ncharts * sizeof(*w_phi0));
    double *w_phase = (double *)ARENA_ALLOC(
        arena, solution->ncharts * sizeof(*w_phase));
    uint8_t *seen = (uint8_t *)ARENA_CALLOC(arena, solution->ncharts, 1);
    size_t rows = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
        char *field[32];
        size_t fields = atgt_split_csv(line, field, 32);
        long chart = -1, group = -1, k = 0;
        double phi0 = 0.0, phase = 0.0;
        if (fields < 18 || atgt_parse_long(field[0], &chart) != 0 ||
            atgt_parse_long(field[2], &group) != 0 ||
            atgt_parse_double(field[12], &phi0) != 0 ||
            atgt_parse_double(field[14], &phase) != 0 ||
            atgt_parse_long(field[17], &k) != 0 ||
            chart < 0 || (size_t)chart >= solution->ncharts ||
            seen[chart] || solution->chart[chart].chart != chart ||
            solution->chart[chart].group != group ||
            solution->chart[chart].winding != k) {
            fclose(fp);
            return -1;
        }
        w_phi0[chart] = phi0;
        w_phase[chart] = phase;
        seen[chart] = 1;
        rows++;
    }
    fclose(fp);
    if (rows != solution->ncharts) return -1;
    int32_t max_group = -1;
    for (size_t c = 0; c < solution->ncharts; c++) {
        if (!seen[c] || solution->chart[c].group < 0 ||
            (size_t)solution->chart[c].group >= solution->ncharts)
            return -1;
        if (solution->chart[c].group > max_group)
            max_group = solution->chart[c].group;
    }
    double *gather = (double *)ARENA_ALLOC(
        arena, solution->ncharts * sizeof(*gather));
    for (int32_t group = 0; group <= max_group; group++) {
        size_t count = 0;
        for (size_t c = 0; c < solution->ncharts; c++)
            if (solution->chart[c].group == group)
                gather[count++] = w_phi0[c] +
                    (double)solution->chart[c].winding - w_phase[c];
        if (count == 0) continue;
        qsort(gather, count, sizeof(*gather), atgt_compare_double);
        double gauge = (double)lround(gather[count / 2]);
        for (size_t c = 0; c < solution->ncharts; c++)
            if (solution->chart[c].group == group)
                solution->chart[c].wind = w_phase[c] + gauge;
    }
    fprintf(stderr,
            "[atlas_track_grow] adapted legacy v1 wind from %s "
            "(%zu charts, %d gauge groups)\n",
            csv_path, solution->ncharts, max_group + 1);
    return 0;
}

static int atgt_read_solution_v1(
    Arena_T arena,
    const char *path,
    const PieceSet *ps,
    AtlasSolution *out)
{
    static const unsigned char magic_expected[8] = {
        'A', 'T', 'L', 'S', 'O', 'L', '1', 0
    };
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    unsigned char magic[8];
    uint32_t version = 0, endian = 0;
    uint64_t nv = 0, nf = 0, nc = 0, nr = 0, fingerprint = 0;
    int io = 0;
    io |= atgt_read_exact(fp, magic, sizeof magic);
    io |= atgt_read_exact(fp, &version, sizeof version);
    io |= atgt_read_exact(fp, &endian, sizeof endian);
    io |= atgt_read_exact(fp, &nv, sizeof nv);
    io |= atgt_read_exact(fp, &nf, sizeof nf);
    io |= atgt_read_exact(fp, &nc, sizeof nc);
    io |= atgt_read_exact(fp, &nr, sizeof nr);
    io |= atgt_read_exact(fp, &fingerprint, sizeof fingerprint);
    if (io != 0 || memcmp(magic, magic_expected, sizeof magic) != 0 ||
        version != 1 || endian != UINT32_C(0x01020304) ||
        nv != (uint64_t)ps->nv || nf != (uint64_t)ps->nf ||
        nc == 0 || nc > (uint64_t)ps->nv || nr > (uint64_t)ps->nf * 8 ||
        fingerprint != AtlasSolution_piece_fingerprint(ps)) {
        fclose(fp);
        return -1;
    }
    memset(out, 0, sizeof *out);
    out->nvertices = (size_t)nv;
    out->nfaces = (size_t)nf;
    out->ncharts = (size_t)nc;
    out->nresiduals = (size_t)nr;
    out->piece_fingerprint = fingerprint;
    out->u = (double *)ARENA_ALLOC(arena, out->nvertices * sizeof(*out->u));
    out->v = (double *)ARENA_ALLOC(arena, out->nvertices * sizeof(*out->v));
    out->face_keep = (uint8_t *)ARENA_ALLOC(
        arena, out->nfaces * sizeof(*out->face_keep));
    out->vertex_chart = (int32_t *)ARENA_ALLOC(
        arena, out->nvertices * sizeof(*out->vertex_chart));
    out->chart = (AtlasSolutionChart *)ARENA_ALLOC(
        arena, out->ncharts * sizeof(*out->chart));
    out->residual = out->nresiduals > 0
        ? (AtlasSolutionResidual *)ARENA_ALLOC(
            arena, out->nresiduals * sizeof(*out->residual))
        : NULL;
    io = 0;
    io |= atgt_read_exact(fp, out->u,
                          out->nvertices * sizeof(*out->u));
    io |= atgt_read_exact(fp, out->v,
                          out->nvertices * sizeof(*out->v));
    io |= atgt_read_exact(fp, out->face_keep,
                          out->nfaces * sizeof(*out->face_keep));
    io |= atgt_read_exact(fp, out->vertex_chart,
                          out->nvertices * sizeof(*out->vertex_chart));
    for (size_t c = 0; c < out->ncharts; c++) {
        int32_t fields[3];
        uint32_t reserved = 0;
        uint64_t rank = 0;
        io |= atgt_read_exact(fp, fields, sizeof fields);
        io |= atgt_read_exact(fp, &reserved, sizeof reserved);
        io |= atgt_read_exact(fp, &rank, sizeof rank);
        io |= atgt_read_exact(fp, &out->chart[c].qc_max,
                              sizeof out->chart[c].qc_max);
        out->chart[c].chart = fields[0];
        out->chart[c].group = fields[1];
        out->chart[c].winding = fields[2];
        out->chart[c].rank = rank;
        out->chart[c].wind = NAN;
        if (reserved != 0 || rank >= nc || fields[0] < 0 ||
            (uint64_t)fields[0] >= nc)
            io = -1;
    }
    for (size_t i = 0; i < out->nresiduals; i++) {
        int32_t fields[6];
        io |= atgt_read_exact(fp, fields, sizeof fields);
        io |= atgt_read_exact(fp, &out->residual[i].min_xyz,
                              sizeof out->residual[i].min_xyz);
        io |= atgt_read_exact(fp, &out->residual[i].radius_delta,
                              sizeof out->residual[i].radius_delta);
        out->residual[i].chart0 = fields[0];
        out->residual[i].chart1 = fields[1];
        out->residual[i].face0 = fields[2];
        out->residual[i].face1 = fields[3];
        out->residual[i].allowed = fields[4];
        out->residual[i].reason = fields[5];
    }
    if (io == 0 && fgetc(fp) != EOF) io = -1;
    fclose(fp);
    if (io != 0 || atgt_load_v1_chart_wind(arena, path, out) != 0)
        return -1;
    return 0;
}

static int atgt_solution_version(const char *path)
{
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    unsigned char magic[8];
    uint32_t version = 0;
    int ok = atgt_read_exact(fp, magic, sizeof magic) == 0 &&
             atgt_read_exact(fp, &version, sizeof version) == 0 &&
             memcmp(magic, "ATLSOL1", 7) == 0;
    fclose(fp);
    return ok && version <= INT32_MAX ? (int)version : -1;
}

static void usage(const char *program)
{
    fprintf(stderr,
        "usage:\n"
        "  %s --selftest\n"
        "  %s <placed_dir> <atlas_solution.bin> <out_dir> <prefix> [options]\n\n"
        "Evidence extraction:\n"
        "  --slice-spacing F          constant-v row spacing (default 4)\n"
        "  --observation-u F           source-U sample spacing (default 2)\n"
        "  --local-xyz F               duplicate cluster radius (default 8)\n"
        "  --tangent-dot F             extraction and growth tangent gate (.35)\n\n"
        "Cycle growth:\n"
        "  --horizontal-gap F          source-U neighbor window (default 24)\n"
        "  --horizontal-stretch F      max chord/source-U scale (default 1.6)\n"
        "  --horizontal-slack F        horizontal metric slack (default 2.5)\n"
        "  --horizontal-min-frac F     min chord/source-U scale (default .35)\n"
        "  --vertical-search-u F       adjacent-row source-U window (default 12)\n"
        "  --vertical-radius F         adjacent-row physical radius (default 8)\n"
        "  --vertical-winding-delta N  max winding change across a physical V link\n"
        "  --quad-normal-dot F         minimum cell triangle cosine (default 0)\n"
        "  --all-local-layers          retain every locally clustered physical\n"
        "                              layer at a collided UV sample\n"
        "  --physical-neighbours       discover sheet H/V edges from local XYZ\n"
        "                              and tangent evidence\n"
        "  --chart-layout PATH         require cross-chart joins to be direct\n"
        "                              cube-sheet relations from this v4 layout\n"
        "  --chart-shards DIR          exact cube-sheet membership used by layout\n"
        "  --vertical-chain-residual F cut cross-row links inconsistent with\n"
        "                              the physical row-chain layout (0=off)\n"
        "  --min-component-quads N     minimum published sheet support (64)\n"
        "  --gap-fill-rounds N         UV->XYZ frontier completion rounds (2)\n"
        "  --gap-fit-rms F             local affine fit RMS gate (default 1.5)\n"
        "  --gap-fit-max F             local affine fit max gate (default 4)\n"
        "  --uv-bridge-u F             max eight-observation ribbon gap (0=off)\n"
        "  --uv-bridge-fit-rms F       bank-fit RMS gate (default .75)\n"
        "  --uv-bridge-fit-max F       bank-fit max gate (default 2)\n"
        "  --stats-only                write metrics + component table, no meshes\n"
        "  --atlas-only                write patch/atlas OBJs + metrics for bake QA\n"
        "  --vmesh-only                write binary UV mesh + metrics for fitting\n"
        "  --no-node-csv               omit the redundant per-vertex CSV dump\n"
        "  --no-winding-export         skip per-(winding x slab) tifxyz atlas\n"
        "  --winding-wraps N           windings per machine piece (default 1)\n"
        "  --winding-slab F            axial slab height (default 4096)\n"
        "  --winding-du F              machine atlas U vox/px (default 1)\n"
        "  --winding-dv F              machine atlas V vox/px (default 1)\n"
        "  --parameter-stay F          checkpoint-U anchor weight (default .01)\n"
        "  --parameter-bridge F        cross-piece evidence weight (default .25)\n"
        "  --parameter-iters N         evidence-coordinate sweeps (default 160)\n",
        program, program);
}

static int valid_prefix(const char *prefix)
{
    if (prefix == NULL || *prefix == '\0') return 0;
    for (const unsigned char *p = (const unsigned char *)prefix; *p; p++)
        if (!isalnum(*p) && *p != '_' && *p != '-' && *p != '.')
            return 0;
    return 1;
}

static FILE *open_output(
    const char *dir,
    const char *prefix,
    const char *suffix)
{
    char path[ATGT_PATH_CAP];
    int n = snprintf(path, sizeof path, "%s/%s%s", dir, prefix, suffix);
    if (n < 0 || (size_t)n >= sizeof path ||
        ves_ensure_parent_dir(path) != 0)
        return NULL;
    return fopen(path, "wb");
}

static int close_output(FILE *fp)
{
    int failed = ferror(fp);
    if (fclose(fp) != 0) failed = 1;
    return failed ? -1 : 0;
}

static int write_exact(FILE *fp, const void *data, size_t bytes)
{
    return bytes == 0 || fwrite(data, 1, bytes, fp) == bytes ? 0 : -1;
}

static int write_patch_vmesh(
    const char *dir,
    const char *prefix,
    const AtlasRibbonObservationSet *set,
    const AtlasTrackGrowResult *result)
{
    FILE *fp = open_output(dir, prefix, "_track_patch.vmesh");
    if (fp == NULL) return -1;
    const unsigned char magic[8] = {'V','E','S','M','E','S','H','1'};
    const uint32_t version = 1;
    const uint32_t endian = UINT32_C(0x01020304);
    const uint32_t flags = 1;
    const uint32_t header_bytes = 64;
    const uint64_t nv = (uint64_t)result->nv;
    const uint64_t nf = (uint64_t)result->nf;
    const uint64_t vertex_offset = header_bytes;
    if (nv > (UINT64_MAX - vertex_offset) / (3 * sizeof(float))) {
        fclose(fp);
        return -1;
    }
    const uint64_t uv_offset = vertex_offset + nv * 3 * sizeof(float);
    if (nv > (UINT64_MAX - uv_offset) / (2 * sizeof(float))) {
        fclose(fp);
        return -1;
    }
    const uint64_t face_offset = uv_offset + nv * 2 * sizeof(float);
    int failed =
        write_exact(fp, magic, sizeof magic) ||
        write_exact(fp, &version, sizeof version) ||
        write_exact(fp, &endian, sizeof endian) ||
        write_exact(fp, &flags, sizeof flags) ||
        write_exact(fp, &header_bytes, sizeof header_bytes) ||
        write_exact(fp, &nv, sizeof nv) ||
        write_exact(fp, &nf, sizeof nf) ||
        write_exact(fp, &vertex_offset, sizeof vertex_offset) ||
        write_exact(fp, &uv_offset, sizeof uv_offset) ||
        write_exact(fp, &face_offset, sizeof face_offset);
    enum { CHUNK = 4096 };
    float xyz[CHUNK * 3];
    float uv[CHUNK * 2];
    for (size_t first = 0; !failed && first < result->nv; first += CHUNK) {
        size_t count = result->nv - first;
        if (count > CHUNK) count = CHUNK;
        for (size_t i = 0; i < count * 3; i++)
            xyz[i] = (float)result->xyz[first * 3 + i];
        failed = write_exact(fp, xyz, count * 3 * sizeof(*xyz));
    }
    for (size_t first = 0; !failed && first < result->nv; first += CHUNK) {
        size_t count = result->nv - first;
        if (count > CHUNK) count = CHUNK;
        for (size_t i = 0; i < count * 2; i++)
            uv[i] = (float)result->uv[first * 2 + i];
        failed = write_exact(fp, uv, count * 2 * sizeof(*uv));
    }
    if (!failed)
        failed = write_exact(
            fp, result->faces, result->nf * 3 * sizeof(*result->faces));
    if (failed) {
        fclose(fp);
        return -1;
    }
    if (close_output(fp) != 0) return -1;
    fp = open_output(dir, prefix, "_track_patch_claimant_island.i32");
    if (fp == NULL) return -1;
    if (write_exact(
            fp, result->vertex_component,
            result->nv * sizeof(*result->vertex_component)) != 0) {
        fclose(fp);
        return -1;
    }
    if (close_output(fp) != 0) return -1;
    fp = open_output(dir, prefix, "_track_patch_source_target.i32");
    if (fp == NULL) return -1;
    if (write_exact(
            fp, result->source_target,
            result->nv * sizeof(*result->source_target)) != 0) {
        fclose(fp);
        return -1;
    }
    if (close_output(fp) != 0) return -1;
    fp = open_output(dir, prefix, "_track_patch_winding.i32");
    if (fp == NULL) return -1;
    if (write_exact(
            fp, result->vertex_winding,
            result->nv * sizeof(*result->vertex_winding)) != 0) {
        fclose(fp);
        return -1;
    }
    if (close_output(fp) != 0) return -1;
    fp = open_output(dir, prefix, "_track_patch_rank.u64");
    if (fp == NULL) return -1;
    if (write_exact(
            fp, result->vertex_rank,
            result->nv * sizeof(*result->vertex_rank)) != 0) {
        fclose(fp);
        return -1;
    }
    if (close_output(fp) != 0) return -1;
    fp = open_output(dir, prefix, "_track_patch_chart.i32");
    if (fp == NULL) return -1;
    int32_t chart[CHUNK];
    for (size_t first = 0; first < result->nv; first += CHUNK) {
        size_t count = result->nv - first;
        if (count > CHUNK) count = CHUNK;
        for (size_t i = 0; i < count; i++) {
            int32_t source = result->source_target[first + i];
            if (source < 0 || (size_t)source >= set->ntarget) {
                fclose(fp);
                return -1;
            }
            chart[i] = set->target[(size_t)source].chart0;
        }
        if (write_exact(fp, chart, count * sizeof(*chart)) != 0) {
            fclose(fp);
            return -1;
        }
    }
    return close_output(fp);
}

static int write_patch_obj(
    const char *dir,
    const char *prefix,
    const AtlasTrackGrowResult *result)
{
    FILE *fp = open_output(dir, prefix, "_track_patch.obj");
    if (fp == NULL) return -1;
    fprintf(fp,
        "# Cycle-supported evidence-grown sheet.\n"
        "# Source-space coordinates retain repository order (z,y,x).\n"
        "# vt coordinates retain the checkpoint-packed global U gauge.\n");
    for (size_t i = 0; i < result->nv; i++) {
        const double *p = &result->xyz[i * 3];
        fprintf(fp, "v %.17g %.17g %.17g 0.18 0.82 0.42\n",
                p[0], p[1], p[2]);
    }
    for (size_t i = 0; i < result->nv; i++)
        fprintf(fp, "vt %.17g %.17g\n",
                result->uv[i * 2], result->uv[i * 2 + 1]);
    fprintf(fp, "s 1\n");
    for (size_t ci = 0; ci < result->ncomponents; ci++) {
        const AtlasTrackGrowComponent *component = &result->component[ci];
        if (!component->selected) continue;
        size_t end = component->first_face + component->faces;
        int previous_kind = -1;
        for (size_t f = component->first_face; f < end; f++) {
            int kind = result->face_inferred[f] ? 1 : 0;
            if (kind != previous_kind) {
                fprintf(fp, "g sheet_%d_w%d_%d_b%d_%s\n", component->id,
                        component->winding_min, component->winding_max,
                        component->atlas_band, kind ? "uvxyz_fill" : "core");
                previous_kind = kind;
            }
            int32_t a = result->faces[f * 3] + 1;
            int32_t b = result->faces[f * 3 + 1] + 1;
            int32_t c = result->faces[f * 3 + 2] + 1;
            fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
    }
    return close_output(fp);
}

static int write_winding_obj(
    const char *dir,
    const char *prefix,
    const AtlasTrackGrowResult *result)
{
    FILE *fp = open_output(dir, prefix, "_track_winding.obj");
    if (fp == NULL) return -1;
    fprintf(fp,
        "# Evidence-grown sheets in the shared source winding/depth frame.\n"
        "# Overlap here is semantic evidence, not atlas packing.\n");
    for (size_t i = 0; i < result->nv; i++) {
        const double *p = &result->xyz[i * 3];
        fprintf(fp, "v %.17g %.17g %.17g\n", p[0], p[1], p[2]);
    }
    for (size_t i = 0; i < result->nv; i++)
        fprintf(fp, "vt %.17g %.17g\n",
                result->global_uv[i * 2], result->global_uv[i * 2 + 1]);
    fprintf(fp, "s 1\n");
    for (size_t ci = 0; ci < result->ncomponents; ci++) {
        const AtlasTrackGrowComponent *component = &result->component[ci];
        if (!component->selected) continue;
        fprintf(fp, "g sheet_%d_w%d_%d\n", component->id,
                component->winding_min, component->winding_max);
        size_t end = component->first_face + component->faces;
        for (size_t f = component->first_face; f < end; f++) {
            int32_t a = result->faces[f * 3] + 1;
            int32_t b = result->faces[f * 3 + 1] + 1;
            int32_t c = result->faces[f * 3 + 2] + 1;
            fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        }
    }
    return close_output(fp);
}

static int write_atlas_obj(
    const char *dir,
    const char *prefix,
    const AtlasTrackGrowResult *result)
{
    FILE *fp = open_output(dir, prefix, "_track_atlas.obj");
    if (fp == NULL) return -1;
    fprintf(fp,
        "# Winding/depth-organized planar view of evidence-grown sheets.\n");
    for (size_t i = 0; i < result->nv; i++)
        fprintf(fp, "v %.17g %.17g 0 0.20 0.65 0.95\n",
                result->uv[i * 2], result->uv[i * 2 + 1]);
    for (size_t ci = 0; ci < result->ncomponents; ci++) {
        const AtlasTrackGrowComponent *component = &result->component[ci];
        if (!component->selected) continue;
        fprintf(fp, "g sheet_%d_w%d_%d_b%d\n", component->id,
                component->winding_min, component->winding_max,
                component->atlas_band);
        size_t end = component->first_face + component->faces;
        for (size_t f = component->first_face; f < end; f++)
            fprintf(fp, "f %d %d %d\n",
                    result->faces[f * 3] + 1,
                    result->faces[f * 3 + 1] + 1,
                    result->faces[f * 3 + 2] + 1);
    }
    return close_output(fp);
}

static int write_nodes_csv(
    const char *dir,
    const char *prefix,
    const AtlasRibbonObservationSet *set,
    const AtlasTrackGrowResult *result)
{
    FILE *fp = open_output(dir, prefix, "_track_nodes.csv");
    if (fp == NULL) return -1;
    fprintf(fp,
        "output_vertex,component,target,row,column,chart0,chart1,source_u,"
        "source_v,local_u,local_v,global_u,global_v,atlas_u,atlas_v,"
        "winding,rank,atlas_band,p0,p1,tangent0,tangent1,z,y,x\n");
    for (size_t i = 0; i < result->nv; i++) {
        int32_t source = result->source_target[i];
        const AtlasRibbonTarget *t = &set->target[source];
        const double *xyz = &result->xyz[i * 3];
        fprintf(fp,
            "%zu,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%.17g,%.17g,%d,%" PRIu64 ",%d,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            i, result->vertex_component[i], source,
            t->row, t->column, t->chart0, t->chart1,
            t->u, t->v,
            result->local_uv[i * 2], result->local_uv[i * 2 + 1],
            result->global_uv[i * 2], result->global_uv[i * 2 + 1],
            result->uv[i * 2], result->uv[i * 2 + 1],
            result->vertex_winding[i], result->vertex_rank[i],
            result->component[result->vertex_component[i]].atlas_band,
            t->p[0], t->p[1], t->tangent[0], t->tangent[1],
            xyz[0], xyz[1], xyz[2]);
    }
    return close_output(fp);
}

static int write_components_csv(
    const char *dir,
    const char *prefix,
    const AtlasTrackGrowResult *result)
{
    FILE *fp = open_output(dir, prefix, "_track_components.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "component,quads,core_quads,inferred_quads,selected,first_vertex,vertices,first_face,"
                "faces,local_u_span,local_v_span,global_u0,global_v0,"
                "atlas_u0,atlas_v0,atlas_band,winding_min,winding_max,"
                "rank_min,rank_max\n");
    for (size_t i = 0; i < result->ncomponents; i++)
        fprintf(fp, "%d,%zu,%zu,%zu,%d,%zu,%zu,%zu,%zu,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%d,%d,%d,%" PRIu64 ",%" PRIu64 "\n",
                result->component[i].id,
                result->component[i].quads,
                result->component[i].core_quads,
                result->component[i].inferred_quads,
                result->component[i].selected,
                result->component[i].first_vertex,
                result->component[i].vertices,
                result->component[i].first_face,
                result->component[i].faces,
                result->component[i].local_u_span,
                result->component[i].local_v_span,
                result->component[i].global_u0,
                result->component[i].global_v0,
                result->component[i].atlas_u0,
                result->component[i].atlas_v0,
                result->component[i].atlas_band,
                result->component[i].winding_min,
                result->component[i].winding_max,
                result->component[i].rank_min,
                result->component[i].rank_max);
    return close_output(fp);
}

static void write_json_string(FILE *fp, const char *text)
{
    fputc('"', fp);
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
        if (*p == '"' || *p == '\\') {
            fputc('\\', fp);
            fputc(*p, fp);
        } else if (*p == '\n') {
            fputs("\\n", fp);
        } else if (*p == '\r') {
            fputs("\\r", fp);
        } else if (*p == '\t') {
            fputs("\\t", fp);
        } else if (*p < 32) {
            fprintf(fp, "\\u%04x", (unsigned)*p);
        } else {
            fputc(*p, fp);
        }
    }
    fputc('"', fp);
}

static int write_stats_json(
    const char *dir,
    const char *prefix,
    const char *placed_dir,
    const char *solution_path,
    const PieceSet *ps,
    const AtlasSolution *solution,
    size_t kept_faces,
    const AtlasRibbonFitOptions *evidence_options,
    const AtlasTrackGrowOptions *grow_options,
    const AtlasRibbonObservationSet *set,
    const AtlasTrackGrowResult *result,
    double load_seconds,
    double evidence_seconds,
    double grow_seconds,
    double artifact_seconds)
{
    FILE *fp = open_output(dir, prefix, "_track_grow_stats.json");
    if (fp == NULL) return -1;
    double accepted_coverage = result->accepted_targets > 0
        ? (double)result->selected_cycle_nodes /
          (double)result->accepted_targets : 0.0;
    double cycle_selection = result->cycle_quads > 0
        ? (double)result->selected_core_quads /
          (double)result->cycle_quads : 0.0;
    fprintf(fp, "{\n  \"schema\": \"atlas_track_grow_v8\",\n");
    fprintf(fp, "  \"placed_dir\": ");
    write_json_string(fp, placed_dir);
    fprintf(fp, ",\n  \"solution_path\": ");
    write_json_string(fp, solution_path);
    fprintf(fp, ",\n  \"prefix\": ");
    write_json_string(fp, prefix);
    fprintf(fp,
        ",\n  \"input\": {\n"
        "    \"cubes\": %zu,\n"
        "    \"vertices\": %zu,\n"
        "    \"faces\": %zu,\n"
        "    \"kept_faces\": %zu,\n"
        "    \"charts\": %zu,\n"
        "    \"piece_fingerprint\": \"%016" PRIx64 "\"\n"
        "  },\n",
        ps->n_cubes, ps->nv, ps->nf, kept_faces, solution->ncharts,
        solution->piece_fingerprint);
    fprintf(fp,
        "  \"evidence_options\": {\n"
        "    \"slice_spacing\": %.17g,\n"
        "    \"observation_u_spacing\": %.17g,\n"
        "    \"local_xyz_tolerance\": %.17g,\n"
        "    \"tangent_dot_min\": %.17g\n"
        "  },\n",
        evidence_options->slice_spacing,
        evidence_options->observation_u_spacing,
        evidence_options->local_xyz_tolerance,
        evidence_options->tangent_dot_min);
    fprintf(fp,
        "  \"grow_options\": {\n"
        "    \"horizontal_max_gap\": %.17g,\n"
        "    \"horizontal_stretch\": %.17g,\n"
        "    \"horizontal_slack\": %.17g,\n"
        "    \"horizontal_min_fraction\": %.17g,\n"
        "    \"vertical_search_u\": %.17g,\n"
        "    \"vertical_radius\": %.17g,\n"
        "    \"vertical_winding_delta_max\": %d,\n"
        "    \"tangent_dot_min\": %.17g,\n"
        "    \"physical_neighbours\": %s,\n"
        "    \"direct_relation_count\": %zu,\n"
        "    \"preserve_layout_parameterization\": %s,\n"
        "    \"vertical_chain_residual_max\": %.17g,\n"
        "    \"quad_normal_dot_min\": %.17g,\n"
        "    \"min_component_quads\": %d,\n"
        "    \"gap_fill_rounds\": %d,\n"
        "    \"gap_fit_rms_max\": %.17g,\n"
        "    \"gap_fit_max_max\": %.17g,\n"
        "    \"uv_bridge_max_u\": %.17g,\n"
        "    \"uv_bridge_fit_rms_max\": %.17g,\n"
        "    \"uv_bridge_fit_max_max\": %.17g,\n"
        "    \"parameter_stay_weight\": %.17g,\n"
        "    \"parameter_bridge_weight\": %.17g,\n"
        "    \"parameter_iterations\": %d\n"
        "  },\n",
        grow_options->horizontal_max_gap,
        grow_options->horizontal_stretch,
        grow_options->horizontal_slack,
        grow_options->horizontal_min_fraction,
        grow_options->vertical_search_u,
        grow_options->vertical_radius,
        grow_options->vertical_winding_delta_max,
        grow_options->tangent_dot_min,
        grow_options->physical_neighbours ? "true" : "false",
        grow_options->direct_relation_count,
        grow_options->preserve_layout_parameterization ? "true" : "false",
        grow_options->vertical_chain_residual_max,
        grow_options->quad_normal_dot_min,
        grow_options->min_component_quads,
        grow_options->gap_fill_rounds,
        grow_options->gap_fit_rms_max,
        grow_options->gap_fit_max_max,
        grow_options->uv_bridge_max_u,
        grow_options->uv_bridge_fit_rms_max,
        grow_options->uv_bridge_fit_max_max,
        grow_options->parameter_stay_weight,
        grow_options->parameter_bridge_weight,
        grow_options->parameter_iterations);
    fprintf(fp,
        "  \"evidence\": {\n"
        "    \"input_faces\": %zu,\n"
        "    \"sliced_faces\": %zu,\n"
        "    \"slice_segments\": %zu,\n"
        "    \"observations\": %zu,\n"
        "    \"targets\": %zu,\n"
        "    \"accepted_targets\": %zu,\n"
        "    \"duplicate_targets\": %zu,\n"
        "    \"conflict_targets\": %zu,\n"
        "    \"conflict_clusters\": %zu,\n"
        "    \"known_topology_conflicts\": %zu,\n"
        "    \"layer_samples\": %zu,\n"
        "    \"remote_layer_samples\": %zu,\n"
        "    \"rows\": %zu\n"
        "  },\n",
        set->input_faces, set->sliced_faces, set->slice_segments,
        set->nobservation, set->ntarget, set->accepted_targets,
        set->duplicate_targets, set->conflict_targets,
        set->conflict_clusters, set->known_topology_conflicts,
        set->nlayer_samples, set->remote_layer_samples, set->nrows);
    fprintf(fp,
        "  \"growth\": {\n"
        "    \"accepted_targets\": %zu,\n"
        "    \"conflict_targets\": %zu,\n"
        "    \"horizontal_candidates\": %zu,\n"
        "    \"horizontal_relation_rejects\": %zu,\n"
        "    \"horizontal_edges\": %zu,\n"
        "    \"vertical_candidates\": %zu,\n"
        "    \"vertical_relation_rejects\": %zu,\n"
        "    \"vertical_winding_rejects\": %zu,\n"
        "    \"vertical_edges\": %zu,\n"
        "    \"vertical_chain_checks\": %zu,\n"
        "    \"vertical_chain_pairs\": %zu,\n"
        "    \"vertical_chain_cuts\": %zu,\n"
        "    \"vertical_chain_residual_rms\": %.17g,\n"
        "    \"vertical_chain_residual_max\": %.17g,\n"
        "    \"vertical_global_checks\": %zu,\n"
        "    \"vertical_global_cuts\": %zu,\n"
        "    \"vertical_global_rounds\": %d,\n"
        "    \"vertical_global_residual_rms\": %.17g,\n"
        "    \"vertical_global_residual_max\": %.17g,\n"
        "    \"quad_candidates\": %zu,\n"
        "    \"quad_twist_rejects\": %zu,\n"
        "    \"cycle_quads\": %zu,\n"
        "    \"gap_fill_candidates\": %zu,\n"
        "    \"gap_fill_evidence_rejects\": %zu,\n"
        "    \"gap_fill_fit_rejects\": %zu,\n"
        "    \"gap_fill_topology_rejects\": %zu,\n"
        "    \"gap_fill_quads\": %zu,\n"
        "    \"gap_fill_pair_quads\": %zu,\n"
        "    \"gap_fill_horizontal_edges\": %zu,\n"
        "    \"gap_fill_vertical_edges\": %zu,\n"
        "    \"uv_bridge_candidates\": %zu,\n"
        "    \"uv_bridge_evidence_rejects\": %zu,\n"
        "    \"uv_bridge_fit_rejects\": %zu,\n"
        "    \"uv_bridge_topology_rejects\": %zu,\n"
        "    \"uv_bridge_quads\": %zu,\n"
        "    \"uv_bridge_horizontal_edges\": %zu,\n"
        "    \"uv_bridge_fit_rms\": %.17g,\n"
        "    \"uv_bridge_fit_max\": %.17g,\n"
        "    \"uv_bridge_u_span_max\": %.17g,\n"
        "    \"selected_core_quads\": %zu,\n"
        "    \"selected_inferred_quads\": %zu,\n"
        "    \"core_uvxyz_fit_rms\": %.17g,\n"
        "    \"core_uvxyz_fit_max\": %.17g,\n"
        "    \"gap_fill_uvxyz_fit_rms\": %.17g,\n"
        "    \"gap_fill_uvxyz_fit_max\": %.17g,\n"
        "    \"cycle_nodes\": %zu,\n"
        "    \"components\": %zu,\n"
        "    \"selected_components\": %zu,\n"
        "    \"selected_quads\": %zu,\n"
        "    \"selected_cycle_nodes\": %zu,\n"
        "    \"atlas_vertices\": %zu,\n"
        "    \"duplicated_point_or_fan_vertices\": %zu,\n"
        "    \"selected_faces\": %zu,\n"
        "    \"selected_rows\": %zu,\n"
        "    \"selected_charts\": %zu,\n"
        "    \"selected_row_min\": %d,\n"
        "    \"selected_row_max\": %d,\n"
        "    \"selected_column_min\": %d,\n"
        "    \"selected_column_max\": %d,\n"
        "    \"accepted_target_coverage\": %.17g,\n"
        "    \"cycle_quad_selection_fraction\": %.17g\n"
        "  },\n",
        result->accepted_targets, result->conflict_targets,
        result->horizontal_candidates, result->horizontal_relation_rejects,
        result->horizontal_edges,
        result->vertical_candidates, result->vertical_relation_rejects,
        result->vertical_winding_rejects,
        result->vertical_edges,
        result->vertical_chain_checks,
        result->vertical_chain_pairs,
        result->vertical_chain_cuts,
        result->vertical_chain_residual_rms,
        result->vertical_chain_residual_max,
        result->vertical_global_checks,
        result->vertical_global_cuts,
        result->vertical_global_rounds,
        result->vertical_global_residual_rms,
        result->vertical_global_residual_max,
        result->quad_candidates, result->quad_twist_rejects,
        result->cycle_quads,
        result->gap_fill_candidates,
        result->gap_fill_evidence_rejects,
        result->gap_fill_fit_rejects,
        result->gap_fill_topology_rejects,
        result->gap_fill_quads,
        result->gap_fill_pair_quads,
        result->gap_fill_horizontal_edges,
        result->gap_fill_vertical_edges,
        result->uv_bridge_candidates,
        result->uv_bridge_evidence_rejects,
        result->uv_bridge_fit_rejects,
        result->uv_bridge_topology_rejects,
        result->uv_bridge_quads,
        result->uv_bridge_horizontal_edges,
        result->uv_bridge_fit_rms,
        result->uv_bridge_fit_max,
        result->uv_bridge_u_span_max,
        result->selected_core_quads,
        result->selected_inferred_quads,
        result->core_fit_rms, result->core_fit_max,
        result->gap_fill_fit_rms, result->gap_fill_fit_max,
        result->cycle_nodes, result->ncomponents,
        result->selected_components, result->selected_quads,
        result->selected_cycle_nodes, result->nv,
        result->nv - result->selected_cycle_nodes, result->nf,
        result->selected_rows, result->selected_charts,
        result->selected_row_min, result->selected_row_max,
        result->selected_column_min, result->selected_column_max,
        accepted_coverage, cycle_selection);
    fprintf(fp,
        "  \"parameterization\": {\n"
        "    \"source_u_span\": %.17g,\n"
        "    \"largest_local_u_span\": %.17g,\n"
        "    \"atlas_u_span\": %.17g,\n"
        "    \"atlas_v_span\": %.17g,\n"
        "    \"atlas_bands\": %d,\n"
        "    \"placed_box_overlaps\": %zu,\n"
        "    \"semantic_winding_depth_overlaps\": %zu,\n"
        "    \"coupling_edges\": %zu,\n"
        "    \"horizontal_bridges\": %zu,\n"
        "    \"vertical_bridges\": %zu,\n"
        "    \"point_bridges\": %zu,\n"
        "    \"parameter_islands\": %zu,\n"
        "    \"winding_direction\": %d,\n"
        "    \"winding_min\": %d,\n"
        "    \"winding_max\": %d,\n"
        "    \"surface_area\": %.17g,\n"
        "    \"horizontal_residual_rms\": %.17g,\n"
        "    \"horizontal_residual_max\": %.17g,\n"
        "    \"vertical_residual_rms\": %.17g,\n"
        "    \"vertical_residual_max\": %.17g,\n"
        "    \"vertical_residual_hist\": [%zu, %zu, %zu, %zu, %zu, "
            "%zu, %zu, %zu, %zu, %zu],\n"
        "    \"bridge_residual_rms\": %.17g,\n"
        "    \"bridge_residual_max\": %.17g,\n"
        "    \"anchor_residual_rms\": %.17g,\n"
        "    \"anchor_residual_max\": %.17g,\n"
        "    \"iterations\": %d,\n"
        "    \"max_change\": %.17g\n"
        "  },\n",
        result->source_u_span, result->local_u_span_max,
        result->atlas_u_span, result->atlas_v_span,
        result->atlas_bands, result->atlas_placed_overlaps,
        result->atlas_winding_overlaps,
        result->parameter_bridge_edges,
        result->parameter_horizontal_bridges,
        result->parameter_vertical_bridges,
        result->parameter_point_bridges,
        result->parameter_islands,
        result->winding_direction,
        result->winding_min, result->winding_max,
        result->surface_area,
        result->horizontal_residual_rms,
        result->horizontal_residual_max,
        result->vertical_residual_rms,
        result->vertical_residual_max,
        result->vertical_residual_hist[0],
        result->vertical_residual_hist[1],
        result->vertical_residual_hist[2],
        result->vertical_residual_hist[3],
        result->vertical_residual_hist[4],
        result->vertical_residual_hist[5],
        result->vertical_residual_hist[6],
        result->vertical_residual_hist[7],
        result->vertical_residual_hist[8],
        result->vertical_residual_hist[9],
        result->parameter_bridge_residual_rms,
        result->parameter_bridge_residual_max,
        result->parameter_anchor_rms,
        result->parameter_anchor_max,
        result->parameter_iterations, result->parameter_max_change);
    fprintf(fp,
        "  \"timing_seconds\": {\n"
        "    \"load\": %.9g,\n"
        "    \"evidence\": %.9g,\n"
        "    \"grow\": %.9g,\n"
        "    \"artifacts_before_stats\": %.9g,\n"
        "    \"total_before_stats\": %.9g\n"
        "  }\n"
        "}\n",
        load_seconds, evidence_seconds, grow_seconds, artifact_seconds,
        load_seconds + evidence_seconds + grow_seconds + artifact_seconds);
    return close_output(fp);
}

static int parse_integer(const char *text, int *out)
{
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 0 || value > INT_MAX)
        return -1;
    *out = (int)value;
    return 0;
}

static int parse_real(const char *text, double *out)
{
    char *end = NULL;
    double value = strtod(text, &end);
    if (end == text || *end != '\0' || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

/* The ordinary ribbon target stream intentionally rejects a UV key whenever
 * remote material lands on that key.  That is correct for a single-valued
 * ribbon field, but wrong for constructing the cube-sheet topology: every
 * locally coherent physical cluster is a real observation and must survive so
 * the sheets can be separated and laid out later.  Layer samples are already
 * sorted by (row, source_column, cluster rank), so this shallow view only
 * expands the target stream; it does not merge or invent geometry. */
static int build_all_local_layer_view(
    Arena_T arena,
    const AtlasRibbonObservationSet *source,
    AtlasRibbonObservationSet *view)
{
    if (arena == NULL || source == NULL || view == NULL ||
        (source->nlayer_samples > 0 && source->layer_sample == NULL))
        return -1;
    *view = *source;
    view->target = (AtlasRibbonTarget *)ARENA_ALLOC(
        arena, (source->nlayer_samples ? source->nlayer_samples : 1) *
               sizeof(*view->target));
    view->ntarget = source->nlayer_samples;
    view->accepted_targets = source->nlayer_samples;
    view->duplicate_targets = 0;
    view->conflict_targets = 0;
    view->conflict_clusters = 0;
    view->known_topology_conflicts = 0;
    for (size_t i = 0; i < source->nlayer_samples; i++) {
        const AtlasRibbonLayerSample *sample = &source->layer_sample[i];
        AtlasRibbonTarget *target = &view->target[i];
        memset(target, 0, sizeof(*target));
        target->row = sample->row;
        target->column = sample->source_column;
        target->chart0 = sample->chart;
        target->chart1 = -1;
        target->face0 = sample->face;
        target->face1 = -1;
        target->known_reason = -1;
        target->observations = sample->charts;
        target->charts = sample->charts;
        target->clusters = 1;
        target->accepted = 1;
        target->u = sample->source_u;
        target->v = sample->v;
        target->p[0] = sample->p[0];
        target->p[1] = sample->p[1];
        target->tangent[0] = sample->tangent[0];
        target->tangent[1] = sample->tangent[1];
    }
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return AtlasTrackGrow_selftest() == 0 ? 0 : 1;
    if (argc < 5) {
        usage(argv[0]);
        return 1;
    }
    const char *placed_dir = argv[1];
    const char *solution_path = argv[2];
    const char *out_dir = argv[3];
    const char *prefix = argv[4];
    if (!valid_prefix(prefix)) {
        fprintf(stderr, "%s: prefix may contain only [A-Za-z0-9_.-]\n",
                argv[0]);
        return 1;
    }

    AtlasRibbonFitOptions evidence_options;
    AtlasRibbonFitOptions_default(&evidence_options);
    AtlasTrackGrowOptions grow_options;
    AtlasTrackGrowOptions_default(&grow_options);
    int winding_export = 1;
    int stats_only = 0;
    int atlas_only = 0;
    int vmesh_only = 0;
    int node_csv = 1;
    int all_local_layers = 0;
    const char *chart_layout_path = NULL;
    const char *chart_shard_dir = NULL;
    int winding_wraps = 1;
    double winding_slab = 4096.0;
    double winding_du = 1.0, winding_dv = 1.0;
    evidence_options.tangent_dot_min = grow_options.tangent_dot_min;
    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--stats-only") == 0) {
            stats_only = 1;
            continue;
        }
        if (strcmp(argv[i], "--atlas-only") == 0) {
            atlas_only = 1;
            continue;
        }
        if (strcmp(argv[i], "--vmesh-only") == 0) {
            vmesh_only = 1;
            continue;
        }
        if (strcmp(argv[i], "--no-node-csv") == 0) {
            node_csv = 0;
            continue;
        }
        if (strcmp(argv[i], "--all-local-layers") == 0) {
            all_local_layers = 1;
            continue;
        }
        if (strcmp(argv[i], "--physical-neighbours") == 0) {
            grow_options.physical_neighbours = 1;
            continue;
        }
        if (strcmp(argv[i], "--chart-layout") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s: missing --chart-layout path\n", argv[0]);
                return 1;
            }
            chart_layout_path = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "--chart-shards") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s: missing --chart-shards directory\n",
                        argv[0]);
                return 1;
            }
            chart_shard_dir = argv[++i];
            continue;
        }
        if (strcmp(argv[i], "--no-winding-export") == 0) {
            winding_export = 0;
            continue;
        }
        if (strcmp(argv[i], "--winding-wraps") == 0) {
            if (i + 1 >= argc || parse_integer(argv[++i], &winding_wraps) ||
                winding_wraps < 1) {
                fprintf(stderr, "%s: invalid --winding-wraps\n", argv[0]);
                return 1;
            }
            continue;
        }
        if (strcmp(argv[i], "--parameter-iters") == 0 ||
            strcmp(argv[i], "--min-component-quads") == 0 ||
            strcmp(argv[i], "--gap-fill-rounds") == 0 ||
            strcmp(argv[i], "--vertical-winding-delta") == 0) {
            const char *option = argv[i];
            int *destination = NULL;
            if (strcmp(option, "--parameter-iters") == 0)
                destination = &grow_options.parameter_iterations;
            else if (strcmp(option, "--gap-fill-rounds") == 0)
                destination = &grow_options.gap_fill_rounds;
            else if (strcmp(option, "--vertical-winding-delta") == 0)
                destination = &grow_options.vertical_winding_delta_max;
            else
                destination = &grow_options.min_component_quads;
            if (i + 1 >= argc ||
                parse_integer(argv[++i], destination) ||
                (destination == &grow_options.min_component_quads &&
                 *destination < 1) ||
                (destination == &grow_options.parameter_iterations &&
                 *destination > 10000) ||
                (destination == &grow_options.gap_fill_rounds &&
                 (*destination < 0 || *destination > 16))) {
                fprintf(stderr, "%s: invalid %s\n", argv[0], option);
                return 1;
            }
            continue;
        }
        double *destination = NULL;
        int shared_tangent = 0;
        if (strcmp(argv[i], "--slice-spacing") == 0)
            destination = &evidence_options.slice_spacing;
        else if (strcmp(argv[i], "--observation-u") == 0)
            destination = &evidence_options.observation_u_spacing;
        else if (strcmp(argv[i], "--local-xyz") == 0)
            destination = &evidence_options.local_xyz_tolerance;
        else if (strcmp(argv[i], "--tangent-dot") == 0) {
            destination = &grow_options.tangent_dot_min;
            shared_tangent = 1;
        } else if (strcmp(argv[i], "--horizontal-gap") == 0)
            destination = &grow_options.horizontal_max_gap;
        else if (strcmp(argv[i], "--horizontal-stretch") == 0)
            destination = &grow_options.horizontal_stretch;
        else if (strcmp(argv[i], "--horizontal-slack") == 0)
            destination = &grow_options.horizontal_slack;
        else if (strcmp(argv[i], "--horizontal-min-frac") == 0)
            destination = &grow_options.horizontal_min_fraction;
        else if (strcmp(argv[i], "--vertical-search-u") == 0)
            destination = &grow_options.vertical_search_u;
        else if (strcmp(argv[i], "--vertical-radius") == 0)
            destination = &grow_options.vertical_radius;
        else if (strcmp(argv[i], "--vertical-chain-residual") == 0)
            destination = &grow_options.vertical_chain_residual_max;
        else if (strcmp(argv[i], "--quad-normal-dot") == 0)
            destination = &grow_options.quad_normal_dot_min;
        else if (strcmp(argv[i], "--gap-fit-rms") == 0)
            destination = &grow_options.gap_fit_rms_max;
        else if (strcmp(argv[i], "--gap-fit-max") == 0)
            destination = &grow_options.gap_fit_max_max;
        else if (strcmp(argv[i], "--uv-bridge-u") == 0)
            destination = &grow_options.uv_bridge_max_u;
        else if (strcmp(argv[i], "--uv-bridge-fit-rms") == 0)
            destination = &grow_options.uv_bridge_fit_rms_max;
        else if (strcmp(argv[i], "--uv-bridge-fit-max") == 0)
            destination = &grow_options.uv_bridge_fit_max_max;
        else if (strcmp(argv[i], "--winding-slab") == 0)
            destination = &winding_slab;
        else if (strcmp(argv[i], "--winding-du") == 0)
            destination = &winding_du;
        else if (strcmp(argv[i], "--winding-dv") == 0)
            destination = &winding_dv;
        else if (strcmp(argv[i], "--parameter-stay") == 0)
            destination = &grow_options.parameter_stay_weight;
        else if (strcmp(argv[i], "--parameter-bridge") == 0)
            destination = &grow_options.parameter_bridge_weight;
        const char *option = argv[i];
        if (destination == NULL || i + 1 >= argc ||
            parse_real(argv[++i], destination) != 0) {
            fprintf(stderr, "%s: unknown or invalid option %s\n",
                    argv[0], option);
            return 1;
        }
        if (shared_tangent)
            evidence_options.tangent_dot_min = grow_options.tangent_dot_min;
    }
    if (!(winding_slab > 0.0) || !(winding_du > 0.0) ||
        !(winding_dv > 0.0)) {
        fprintf(stderr, "%s: winding slab/du/dv must be positive\n", argv[0]);
        return 1;
    }
    if ((chart_layout_path == NULL) != (chart_shard_dir == NULL)) {
        fprintf(stderr,
                "%s: --chart-layout and --chart-shards must be supplied "
                "together\n", argv[0]);
        return 1;
    }
    if (stats_only && atlas_only) {
        fprintf(stderr, "%s: --stats-only and --atlas-only are exclusive\n",
                argv[0]);
        return 1;
    }

    char probe[ATGT_PATH_CAP];
    int probe_n = snprintf(probe, sizeof probe, "%s/%s.probe", out_dir, prefix);
    if (probe_n < 0 || (size_t)probe_n >= sizeof probe ||
        ves_ensure_parent_dir(probe) != 0) {
        fprintf(stderr, "%s: cannot create output directory %s\n",
                argv[0], out_dir);
        return 1;
    }

    double start = ves_clock_sec();
    Arena_T arena = Arena_new();
    if (arena == NULL) return 1;
    PieceSet ps;
    if (PieceSet_build(arena, placed_dir, &ps) != 0) {
        fprintf(stderr, "%s: cannot load placed directory %s\n",
                argv[0], placed_dir);
        Arena_dispose(&arena);
        return 1;
    }
    ScaffoldCalib calibration;
    if (Scaffold_read_calib(placed_dir, &calibration) != 0) {
        fprintf(stderr, "%s: missing placed calibration in %s\n",
                argv[0], placed_dir);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr,
        "[atlas_track_grow] loaded cubes=%zu vertices=%zu faces=%zu "
        "fingerprint=%016" PRIx64 "\n",
        ps.n_cubes, ps.nv, ps.nf, AtlasSolution_piece_fingerprint(&ps));
    AtlasSolution solution;
    int solution_version = atgt_solution_version(solution_path);
    int solution_read = solution_version == 1
        ? atgt_read_solution_v1(arena, solution_path, &ps, &solution)
        : AtlasSolution_read(arena, solution_path, &ps, &solution);
    if (solution_read != 0) {
        fprintf(stderr, "%s: checkpoint does not match placed input: %s\n",
                argv[0], solution_path);
        Arena_dispose(&arena);
        return 1;
    }
    int wind_ranked = AtlasSolution_rank_by_wind(arena, &solution);
    if (wind_ranked < 0) {
        fprintf(stderr, "%s: checkpoint wind field invalid: %s\n",
                argv[0], solution_path);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr, "[atlas_track_grow] chart ranks: %s\n",
            wind_ranked ? "wind order" : "producer relayout order");
    if (chart_layout_path != NULL) {
        int32_t *relation_a = NULL, *relation_b = NULL;
        double *layout_u = NULL;
        size_t relation_count = 0, layout_nodes = 0, layout_charts = 0;
        if (atgt_load_direct_chart_relations(
                arena, chart_layout_path, chart_shard_dir, &ps, &solution,
                &relation_a, &relation_b,
                &relation_count, &layout_nodes, &layout_charts,
                &layout_u) != 0) {
            fprintf(stderr,
                    "%s: chart layout does not match placed input: %s\n",
                    argv[0], chart_layout_path);
            Arena_dispose(&arena);
            return 1;
        }
        grow_options.direct_relation_a = relation_a;
        grow_options.direct_relation_b = relation_b;
        grow_options.direct_relation_count = relation_count;
        grow_options.preserve_layout_parameterization = 1;
        solution.u = layout_u;
        fprintf(stderr,
                "[atlas_track_grow] cube-sheet graph: live_nodes=%zu "
                "mapped_charts=%zu direct_chart_pairs=%zu orphan_charts=%zu\n",
                layout_nodes, layout_charts, relation_count,
                solution.ncharts - layout_charts);
    }
    size_t kept_faces = 0;
    for (size_t f = 0; f < solution.nfaces; f++)
        if (solution.face_keep[f]) kept_faces++;
    double loaded = ves_clock_sec();
    fprintf(stderr,
        "[atlas_track_grow] cubes=%zu vertices=%zu faces=%zu kept=%zu "
        "charts=%zu\n",
        ps.n_cubes, ps.nv, ps.nf, kept_faces, solution.ncharts);

    AtlasRibbonObservationSet evidence;
    if (AtlasRibbonFit_build_observations(
            arena, &ps, &solution, &calibration,
            &evidence_options, &evidence) != 0) {
        fprintf(stderr, "%s: section evidence extraction failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    double extracted = ves_clock_sec();
    fprintf(stderr,
        "[atlas_track_grow] observations=%zu targets=%zu accepted=%zu "
        "conflicts=%zu layers=%zu rows=%zu\n",
        evidence.nobservation, evidence.ntarget, evidence.accepted_targets,
        evidence.conflict_targets, evidence.nlayer_samples, evidence.nrows);

    AtlasRibbonObservationSet growth_evidence = evidence;
    if (all_local_layers &&
        build_all_local_layer_view(arena, &evidence, &growth_evidence) != 0) {
        fprintf(stderr, "%s: local-layer target expansion failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    if (all_local_layers)
        fprintf(stderr,
            "[atlas_track_grow] growth view: all %zu locally clustered "
            "physical layers retained\n",
            growth_evidence.ntarget);

    AtlasTrackGrowResult result;
    if (AtlasTrackGrow_build(
            arena, &growth_evidence, &calibration,
            &grow_options, &result) != 0) {
        fprintf(stderr, "%s: evidence track growth failed (check options)\n",
                argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    double grown = ves_clock_sec();
    fprintf(stderr,
        "[atlas_track_grow] H=%zu/%zu Hrel=%zu V=%zu/%zu Vrel=%zu "
        "Vcut=%zu+%zu/%zu core_quads=%zu "
        "fill=%zu (pair=%zu H=%zu V=%zu fit_rms=%.4g) twist_reject=%zu "
        "uv_bank=%zu/%zu (edges=%zu fit_rms=%.4g max_u=%.4g) "
        "components=%zu selected_components=%zu selected_quads=%zu\n",
        result.horizontal_edges, result.horizontal_candidates,
        result.horizontal_relation_rejects,
        result.vertical_edges, result.vertical_candidates,
        result.vertical_relation_rejects,
        result.vertical_chain_cuts, result.vertical_global_cuts,
        result.vertical_chain_checks,
        result.cycle_quads, result.gap_fill_quads,
        result.gap_fill_pair_quads,
        result.gap_fill_horizontal_edges,
        result.gap_fill_vertical_edges,
        result.gap_fill_fit_rms, result.quad_twist_rejects,
        result.uv_bridge_quads, result.uv_bridge_candidates,
        result.uv_bridge_horizontal_edges,
        result.uv_bridge_fit_rms, result.uv_bridge_u_span_max,
        result.ncomponents, result.selected_components,
        result.selected_quads);
    fprintf(stderr,
        "[atlas_track_grow] output vertices=%zu faces=%zu rows=%zu charts=%zu "
        "windings=%d..%d source_U_span=%.6g max_local_U=%.6g "
        "global_atlas=%.6gx%.6g islands=%zu bridges=%zu "
        "H_rms=%.6g V_rms=%.6g bridge_rms=%.6g\n",
        result.nv, result.nf, result.selected_rows, result.selected_charts,
        result.winding_min, result.winding_max,
        result.source_u_span, result.local_u_span_max,
        result.atlas_u_span, result.atlas_v_span,
        result.parameter_islands, result.parameter_bridge_edges,
        result.horizontal_residual_rms, result.vertical_residual_rms,
        result.parameter_bridge_residual_rms);

    int io = 0;
    AtlasTrackExportStats winding_stats;
    memset(&winding_stats, 0, sizeof winding_stats);
    if (stats_only) {
        io |= write_components_csv(out_dir, prefix, &result);
    } else if (vmesh_only) {
        io |= write_patch_vmesh(
            out_dir, prefix, &growth_evidence, &result);
        io |= write_components_csv(out_dir, prefix, &result);
    } else if (atlas_only) {
        io |= write_patch_obj(out_dir, prefix, &result);
        io |= write_atlas_obj(out_dir, prefix, &result);
        io |= write_components_csv(out_dir, prefix, &result);
    } else {
        io |= write_patch_obj(out_dir, prefix, &result);
        io |= write_winding_obj(out_dir, prefix, &result);
        io |= write_atlas_obj(out_dir, prefix, &result);
        if (node_csv)
            io |= write_nodes_csv(
                out_dir, prefix, &growth_evidence, &result);
        io |= write_components_csv(out_dir, prefix, &result);
    }
    if (!stats_only && !atlas_only && !vmesh_only && winding_export) {
        char winding_root[ATGT_PATH_CAP];
        int n = snprintf(winding_root, sizeof winding_root,
                         "%s/%s_track_winding_atlas", out_dir, prefix);
        if (n < 0 || (size_t)n >= sizeof winding_root ||
            AtlasTrackExportWinding_run(
                arena, &growth_evidence, &calibration, &result,
                winding_root, prefix, winding_wraps, winding_slab,
                winding_du, winding_dv, &winding_stats) != 0) {
            fprintf(stderr, "%s: winding atlas export failed\n", argv[0]);
            io = -1;
        } else {
            fprintf(stderr,
                "[atlas_track_grow] winding atlas: %zu/%zu pieces written "
                "valid=%zu conflict=%zu (%.4f%%) quarantine_faces=%zu "
                "in %.3fs\n",
                winding_stats.written, winding_stats.pieces,
                winding_stats.valid_pixels, winding_stats.conflict_pixels,
                100.0 * winding_stats.conflict_fraction,
                winding_stats.quarantine_faces, winding_stats.seconds);
        }
    }
    double artifacts_done = ves_clock_sec();
    io |= write_stats_json(
        out_dir, prefix, placed_dir, solution_path, &ps, &solution,
        kept_faces, &evidence_options, &grow_options, &evidence, &result,
        loaded - start, extracted - loaded, grown - extracted,
        artifacts_done - grown);
    if (io != 0) {
        fprintf(stderr, "%s: one or more output writes failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr, "[atlas_track_grow] wrote %s/%s_track_%s (%.3f s)\n",
            out_dir, prefix,
            stats_only ? "metrics" : (atlas_only ? "atlas+metrics" : "*"),
            ves_clock_sec() - start);
    Arena_dispose(&arena);
    return 0;
}
