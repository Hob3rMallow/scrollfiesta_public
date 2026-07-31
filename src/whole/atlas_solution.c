#include "atlas_solution.h"

#include "../common/ves_platform.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    /* v2: per-chart record grew `wind` (absolute continuous wind, turns).
     * v1 checkpoints carried only the round-local tabu k as "winding" --
     * all zero at --rounds convergence, which under-counted every
     * multi-turn spiral fill downstream.  v1 is rejected; regenerate by
     * rerunning atlas_overlap_fix. */
    ATLAS_SOLUTION_VERSION = 2
};

static const unsigned char atlas_solution_magic[8] = {
    'A', 'T', 'L', 'S', 'O', 'L', '1', 0
};

static uint64_t as_fnv_bytes(uint64_t hash, const void *data, size_t size)
{
    const unsigned char *p = (const unsigned char *)data;
    for (size_t i = 0; i < size; i++) {
        hash ^= (uint64_t)p[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

uint64_t AtlasSolution_piece_fingerprint(const PieceSet *ps)
{
    if (ps == NULL) return 0;
    uint64_t hash = UINT64_C(1469598103934665603);
    uint64_t count[3] = {
        (uint64_t)ps->nv, (uint64_t)ps->nf, (uint64_t)ps->n_cubes
    };
    hash = as_fnv_bytes(hash, count, sizeof count);
    if (ps->verts != NULL)
        hash = as_fnv_bytes(hash, ps->verts, ps->nv * 3 * sizeof(float));
    if (ps->faces != NULL)
        hash = as_fnv_bytes(hash, ps->faces, ps->nf * 3 * sizeof(int32_t));
    if (ps->face_cube != NULL)
        hash = as_fnv_bytes(hash, ps->face_cube,
                            ps->nf * sizeof(int32_t));
    if (ps->cube_voff != NULL)
        hash = as_fnv_bytes(hash, ps->cube_voff,
                            (ps->n_cubes + 1) * sizeof(size_t));
    if (ps->ids != NULL)
        hash = as_fnv_bytes(hash, ps->ids, ps->n_cubes * sizeof ps->ids[0]);
    return hash;
}

static int as_write(FILE *fp, const void *data, size_t size)
{
    return size == 0 || fwrite(data, 1, size, fp) == size ? 0 : -1;
}

static int as_read(FILE *fp, void *data, size_t size)
{
    return size == 0 || fread(data, 1, size, fp) == size ? 0 : -1;
}

static int as_write_u32(FILE *fp, uint32_t value)
{
    return as_write(fp, &value, sizeof value);
}

static int as_write_u64(FILE *fp, uint64_t value)
{
    return as_write(fp, &value, sizeof value);
}

static int as_read_u32(FILE *fp, uint32_t *value)
{
    return as_read(fp, value, sizeof *value);
}

static int as_read_u64(FILE *fp, uint64_t *value)
{
    return as_read(fp, value, sizeof *value);
}

int AtlasSolution_write(const char *path,
                        const PieceSet *ps,
                        const double *u,
                        const double *v,
                        const uint8_t *face_keep,
                        const AtlasOverlapFixTabuReport *report)
{
    if (path == NULL || ps == NULL || ps->nv == 0 || ps->nf == 0 ||
        u == NULL || v == NULL || face_keep == NULL || report == NULL ||
        report->nchart == 0 || report->chart == NULL ||
        report->vertex_chart == NULL)
        return -1;
    if (ves_ensure_parent_dir(path) != 0) return -1;

    FILE *fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    uint32_t endian = UINT32_C(0x01020304);
    uint64_t fingerprint = AtlasSolution_piece_fingerprint(ps);
    uint64_t nresidual = report->relayout_residual != NULL
                       ? (uint64_t)report->nrelayout_residual : 0;

    int io = 0;
    io |= as_write(fp, atlas_solution_magic, sizeof atlas_solution_magic);
    io |= as_write_u32(fp, ATLAS_SOLUTION_VERSION);
    io |= as_write_u32(fp, endian);
    io |= as_write_u64(fp, (uint64_t)ps->nv);
    io |= as_write_u64(fp, (uint64_t)ps->nf);
    io |= as_write_u64(fp, (uint64_t)report->nchart);
    io |= as_write_u64(fp, nresidual);
    io |= as_write_u64(fp, fingerprint);
    io |= as_write(fp, u, ps->nv * sizeof(double));
    io |= as_write(fp, v, ps->nv * sizeof(double));
    io |= as_write(fp, face_keep, ps->nf * sizeof(uint8_t));
    io |= as_write(fp, report->vertex_chart, ps->nv * sizeof(int32_t));

    for (size_t c = 0; c < report->nchart; c++) {
        const AtlasOverlapFixChart *src = &report->chart[c];
        int32_t fields[3] = {src->chart, src->group, src->k};
        uint32_t reserved = 0;
        io |= as_write(fp, fields, sizeof fields);
        io |= as_write_u32(fp, reserved);
        io |= as_write_u64(fp, (uint64_t)src->relayout_rank);
        io |= as_write(fp, &src->qc_max, sizeof src->qc_max);
        io |= as_write(fp, &src->wind, sizeof src->wind);
    }
    for (size_t i = 0; i < (size_t)nresidual; i++) {
        const AtlasOverlapFixRelayoutResidual *src =
            &report->relayout_residual[i];
        int32_t fields[6] = {
            src->chart0, src->chart1, src->face0, src->face1,
            src->allowed, src->reason
        };
        io |= as_write(fp, fields, sizeof fields);
        io |= as_write(fp, &src->min_xyz, sizeof src->min_xyz);
        io |= as_write(fp, &src->radius_delta, sizeof src->radius_delta);
    }

    if (fclose(fp) != 0) io = -1;
    return io == 0 ? 0 : -1;
}

int AtlasSolution_read(Arena_T arena,
                       const char *path,
                       const PieceSet *ps,
                       AtlasSolution *out)
{
    if (arena == NULL || path == NULL || ps == NULL || out == NULL)
        return -1;
    memset(out, 0, sizeof *out);
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;

    unsigned char magic[8];
    uint32_t version = 0, endian = 0;
    uint64_t nv = 0, nf = 0, nc = 0, nr = 0, fingerprint = 0;
    int io = 0;
    io |= as_read(fp, magic, sizeof magic);
    io |= as_read_u32(fp, &version);
    io |= as_read_u32(fp, &endian);
    io |= as_read_u64(fp, &nv);
    io |= as_read_u64(fp, &nf);
    io |= as_read_u64(fp, &nc);
    io |= as_read_u64(fp, &nr);
    io |= as_read_u64(fp, &fingerprint);
    if (io != 0 ||
        memcmp(magic, atlas_solution_magic, sizeof magic) != 0 ||
        version != ATLAS_SOLUTION_VERSION ||
        endian != UINT32_C(0x01020304) ||
        nv != (uint64_t)ps->nv || nf != (uint64_t)ps->nf ||
        nc == 0 || nc > (uint64_t)SIZE_MAX ||
        nr > (uint64_t)SIZE_MAX ||
        fingerprint != AtlasSolution_piece_fingerprint(ps)) {
        fclose(fp);
        return -1;
    }

    out->nvertices = (size_t)nv;
    out->nfaces = (size_t)nf;
    out->ncharts = (size_t)nc;
    out->nresiduals = (size_t)nr;
    out->piece_fingerprint = fingerprint;
    out->u = (double *)ARENA_ALLOC(arena, out->nvertices * sizeof(double));
    out->v = (double *)ARENA_ALLOC(arena, out->nvertices * sizeof(double));
    out->face_keep = (uint8_t *)ARENA_ALLOC(
        arena, out->nfaces * sizeof(uint8_t));
    out->vertex_chart = (int32_t *)ARENA_ALLOC(
        arena, out->nvertices * sizeof(int32_t));
    out->chart = (AtlasSolutionChart *)ARENA_ALLOC(
        arena, out->ncharts * sizeof(AtlasSolutionChart));
    out->residual = out->nresiduals > 0
        ? (AtlasSolutionResidual *)ARENA_ALLOC(
            arena, out->nresiduals * sizeof(AtlasSolutionResidual))
        : NULL;

    io = 0;
    io |= as_read(fp, out->u, out->nvertices * sizeof(double));
    io |= as_read(fp, out->v, out->nvertices * sizeof(double));
    io |= as_read(fp, out->face_keep,
                  out->nfaces * sizeof(uint8_t));
    io |= as_read(fp, out->vertex_chart,
                  out->nvertices * sizeof(int32_t));
    for (size_t c = 0; c < out->ncharts; c++) {
        int32_t fields[3];
        uint32_t reserved = 0;
        uint64_t rank = 0;
        io |= as_read(fp, fields, sizeof fields);
        io |= as_read_u32(fp, &reserved);
        io |= as_read_u64(fp, &rank);
        io |= as_read(fp, &out->chart[c].qc_max,
                      sizeof out->chart[c].qc_max);
        io |= as_read(fp, &out->chart[c].wind,
                      sizeof out->chart[c].wind);
        out->chart[c].chart = fields[0];
        out->chart[c].group = fields[1];
        out->chart[c].winding = fields[2];
        out->chart[c].rank = rank;
        if (reserved != 0 || rank >= nc ||
            fields[0] < 0 || (uint64_t)fields[0] >= nc ||
            !isfinite(out->chart[c].wind))
            io = -1;
    }
    for (size_t i = 0; i < out->nresiduals; i++) {
        int32_t fields[6];
        io |= as_read(fp, fields, sizeof fields);
        io |= as_read(fp, &out->residual[i].min_xyz,
                      sizeof out->residual[i].min_xyz);
        io |= as_read(fp, &out->residual[i].radius_delta,
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
    return io == 0 ? 0 : -1;
}

typedef struct {
    double wind;
    uint64_t old_rank;
    size_t chart;
} AsWindRank;

static int as_compare_wind_rank(const void *pa, const void *pb)
{
    const AsWindRank *a = (const AsWindRank *)pa;
    const AsWindRank *b = (const AsWindRank *)pb;
    if (a->wind != b->wind) return a->wind < b->wind ? -1 : 1;
    if (a->old_rank != b->old_rank) return a->old_rank < b->old_rank ? -1 : 1;
    return 0;
}

int AtlasSolution_rank_by_wind(Arena_T arena, AtlasSolution *solution)
{
    if (arena == NULL || solution == NULL || solution->ncharts == 0 ||
        solution->chart == NULL)
        return -1;
    double wind_min = DBL_MAX, wind_max = -DBL_MAX;
    for (size_t c = 0; c < solution->ncharts; c++) {
        double w = solution->chart[c].wind;
        if (!isfinite(w)) return -1;
        if (w < wind_min) wind_min = w;
        if (w > wind_max) wind_max = w;
    }
    /* A degenerate wind field (single wrap, or a fixture) carries no order
     * information; keep the producer's relayout ranks untouched. */
    if (!(wind_max - wind_min > 0.5)) return 0;
    AsWindRank *order = (AsWindRank *)ARENA_ALLOC(
        arena, solution->ncharts * sizeof(*order));
    for (size_t c = 0; c < solution->ncharts; c++) {
        order[c].wind = solution->chart[c].wind;
        order[c].old_rank = solution->chart[c].rank;
        order[c].chart = c;
    }
    qsort(order, solution->ncharts, sizeof(*order), as_compare_wind_rank);
    for (size_t r = 0; r < solution->ncharts; r++)
        solution->chart[order[r].chart].rank = (uint64_t)r;
    return 1;
}

int AtlasSolution_selftest(Arena_T arena, const char *path)
{
    if (arena == NULL || path == NULL) return 1;
    float vertices[9] = {
        0, 0, 0,
        0, 1, 0,
        0, 0, 1
    };
    int32_t faces[3] = {0, 1, 2};
    int32_t face_cube[1] = {0};
    char ids[1][48];
    memset(ids, 0, sizeof ids);
    memcpy(ids[0], "z00000_y00000_x00000", 21);
    size_t cube_voff[2] = {0, 3};
    PieceSet ps;
    memset(&ps, 0, sizeof ps);
    ps.verts = vertices;
    ps.faces = faces;
    ps.face_cube = face_cube;
    ps.ids = ids;
    ps.cube_voff = cube_voff;
    ps.nv = 3;
    ps.nf = 1;
    ps.n_cubes = 1;

    double u[3] = {10, 11, 10};
    double v[3] = {20, 20, 21};
    uint8_t keep[1] = {1};
    int32_t vertex_chart[3] = {0, 0, 0};
    AtlasOverlapFixChart chart;
    memset(&chart, 0, sizeof chart);
    chart.chart = 0;
    chart.group = 7;
    chart.k = -2;
    chart.wind = -2.25; /* non-integer on purpose: continuous turns */
    chart.qc_max = 3.5;
    chart.relayout_rank = 0;
    AtlasOverlapFixRelayoutResidual residual;
    memset(&residual, 0, sizeof residual);
    residual.chart0 = 0;
    residual.chart1 = 0;
    residual.face0 = 0;
    residual.face1 = 0;
    residual.reason = 3;
    residual.min_xyz = 0.0;
    AtlasOverlapFixTabuReport report;
    memset(&report, 0, sizeof report);
    report.chart = &chart;
    report.nchart = 1;
    report.vertex_chart = vertex_chart;
    report.relayout_residual = &residual;
    report.nrelayout_residual = 1;

    if (AtlasSolution_write(path, &ps, u, v, keep, &report) != 0)
        return 1;
    AtlasSolution loaded;
    if (AtlasSolution_read(arena, path, &ps, &loaded) != 0)
        return 1;
    if (loaded.nvertices != 3 || loaded.nfaces != 1 ||
        loaded.ncharts != 1 || loaded.nresiduals != 1 ||
        memcmp(loaded.u, u, sizeof u) != 0 ||
        memcmp(loaded.v, v, sizeof v) != 0 ||
        loaded.face_keep[0] != 1 ||
        loaded.vertex_chart[2] != 0 ||
        loaded.chart[0].group != 7 ||
        loaded.chart[0].winding != -2 ||
        loaded.chart[0].wind != -2.25 ||
        loaded.chart[0].qc_max != 3.5 ||
        loaded.residual[0].reason != 3)
        return 1;
    fprintf(stderr, "[atlas_solution selftest] PASS\n");
    return 0;
}
