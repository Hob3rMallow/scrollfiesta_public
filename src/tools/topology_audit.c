/*
 * topology_audit.c -- native CLI for the component topology invariant suite.
 *
 * Usage:
 *   topology_audit mesh.obj [--json report.json]
 *       [--generators generators.obj]
 *       [--max-generators-per-component N] [--top N]
 *       [--fail-nondisk] [--require-consistent-winding]
 *   topology_audit --selftest
 *
 * Exit codes: 0 success, 1 IO/analysis failure, 2 usage, 3 selftest failure,
 *             4 requested topology gate failed.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "../flatten/topology_invariants.h"
#include "../flatten/disk_topology_repair.h"
#include "../flatten/seam_cut.h"
#include "../common/arena.h"
#include "../common/obj_io.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *input;
    const char *json;
    const char *generators;
    const char *repair_disk;
    size_t max_generators_per_component;
    size_t top;
    int fail_nondisk;
    int require_consistent_winding;
} CliOptions;

typedef struct {
    uint32_t bit;
    const char *name;
} DefectName;

static const DefectName DEFECT_NAME[] = {
    { TOPOLOGY_DEFECT_NONMANIFOLD_EDGE, "nonmanifold_edge" },
    { TOPOLOGY_DEFECT_NONMANIFOLD_VERTEX, "nonmanifold_vertex" },
    { TOPOLOGY_DEFECT_IRREGULAR_BOUNDARY, "irregular_boundary" },
    { TOPOLOGY_DEFECT_NONORIENTABLE, "nonorientable" },
    { TOPOLOGY_DEFECT_BOUNDARY_COUNT, "boundary_count" },
    { TOPOLOGY_DEFECT_NONTRIVIAL_H1, "nontrivial_h1" },
    { TOPOLOGY_DEFECT_NONTRIVIAL_H2, "nontrivial_h2" },
    { TOPOLOGY_DEFECT_INCONSISTENT_WINDING, "inconsistent_winding" },
    { TOPOLOGY_DEFECT_GENERATOR_FAILURE, "generator_failure" },
    { TOPOLOGY_DEFECT_ALGEBRA_MISMATCH, "algebra_mismatch" }
};

static void usage(FILE *stream)
{
    fprintf(stream,
        "usage:\n"
        "  topology_audit mesh.obj [--json report.json]\n"
        "      [--generators generators.obj]\n"
        "      [--repair-disk repaired.obj]\n"
        "      [--max-generators-per-component N] [--top N]\n"
        "      [--fail-nondisk] [--require-consistent-winding]\n"
        "  topology_audit --selftest\n");
}

static int parse_size(const char *text, size_t *out)
{
    char *end = NULL;
    unsigned long long value;
    if (!text || !text[0] || text[0] == '-') return 0;
    errno = 0;
    value = strtoull(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' ||
        value > (unsigned long long)SIZE_MAX) return 0;
    *out = (size_t)value;
    return 1;
}

static int parse_cli(int argc, char **argv, CliOptions *options)
{
    memset(options, 0, sizeof(*options));
    options->top = 25;
    if (argc < 2) return -1;
    options->input = argv[1];
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--json") && i + 1 < argc) {
            options->json = argv[++i];
        } else if (!strcmp(argv[i], "--generators") && i + 1 < argc) {
            options->generators = argv[++i];
        } else if (!strcmp(argv[i], "--repair-disk") && i + 1 < argc) {
            options->repair_disk = argv[++i];
        } else if (!strcmp(argv[i], "--max-generators-per-component") &&
                   i + 1 < argc) {
            if (!parse_size(argv[++i],
                            &options->max_generators_per_component))
                return -1;
        } else if (!strcmp(argv[i], "--top") && i + 1 < argc) {
            if (!parse_size(argv[++i], &options->top)) return -1;
        } else if (!strcmp(argv[i], "--fail-nondisk")) {
            options->fail_nondisk = 1;
        } else if (!strcmp(argv[i], "--require-consistent-winding") ||
                   !strcmp(argv[i], "--require-oriented")) {
            options->require_consistent_winding = 1;
        } else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
            return 1;
        } else {
            return -1;
        }
    }
    return 0;
}

static void print_betti(FILE *stream, int64_t value)
{
    if (value < 0) fputs("?", stream);
    else fprintf(stream, "%lld", (long long)value);
}

static void print_defects(FILE *stream, uint32_t mask)
{
    int first = 1;
    if (!mask) {
        fputs("-", stream);
        return;
    }
    for (size_t i = 0;
         i < sizeof(DEFECT_NAME) / sizeof(DEFECT_NAME[0]); i++) {
        if (!(mask & DEFECT_NAME[i].bit)) continue;
        fprintf(stream, "%s%s", first ? "" : ",", DEFECT_NAME[i].name);
        first = 0;
    }
}

static void json_string(FILE *stream, const char *text)
{
    fputc('"', stream);
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
        switch (*p) {
        case '"': fputs("\\\"", stream); break;
        case '\\': fputs("\\\\", stream); break;
        case '\b': fputs("\\b", stream); break;
        case '\f': fputs("\\f", stream); break;
        case '\n': fputs("\\n", stream); break;
        case '\r': fputs("\\r", stream); break;
        case '\t': fputs("\\t", stream); break;
        default:
            if (*p < 0x20) fprintf(stream, "\\u%04x", (unsigned)*p);
            else fputc(*p, stream);
        }
    }
    fputc('"', stream);
}

static void json_betti(FILE *stream, int64_t value)
{
    if (value < 0) fputs("null", stream);
    else fprintf(stream, "%lld", (long long)value);
}

static void json_defects(FILE *stream, uint32_t mask)
{
    int first = 1;
    fputc('[', stream);
    for (size_t i = 0;
         i < sizeof(DEFECT_NAME) / sizeof(DEFECT_NAME[0]); i++) {
        if (!(mask & DEFECT_NAME[i].bit)) continue;
        if (!first) fputs(", ", stream);
        json_string(stream, DEFECT_NAME[i].name);
        first = 0;
    }
    fputc(']', stream);
}

static int write_json(const char *path, const char *input,
                      const TopologyAuditReport *report)
{
    FILE *stream = fopen(path, "wb");
    if (!stream) return -1;
    fputs("{\n  \"schema\": \"scrollfiesta.topology_audit.v2\",\n"
          "  \"coefficients\": \"F2\",\n"
          "  \"generator_method\": "
          "\"shortest-path tree / maximum-weight extended dual cotree\",\n"
          "  \"generator_optimality\": "
          "\"minimum cardinality; minimum rooted weight for the chosen primal tree\",\n"
          "  \"input\": ", stream);
    json_string(stream, input);
    fprintf(stream,
        ",\n  \"vertices\": %zu,\n"
        "  \"edges\": %zu,\n"
        "  \"faces\": %zu,\n"
        "  \"face_components\": %zu,\n"
        "  \"isolated_vertices\": %zu,\n"
        "  \"topological_components\": %zu,\n"
        "  \"euler_characteristic\": %lld,\n"
        "  \"betti\": [",
        report->vertices, report->edges, report->faces,
        report->face_components, report->isolated_vertices,
        report->topological_components,
        (long long)report->euler_characteristic);
    json_betti(stream, report->beta_0); fputs(", ", stream);
    json_betti(stream, report->beta_1); fputs(", ", stream);
    json_betti(stream, report->beta_2);
    fprintf(stream,
        "],\n  \"betti_complete\": %s,\n"
        "  \"disk_components\": %zu,\n"
        "  \"nondisk_components\": %zu,\n"
        "  \"invalid_surface_components\": %zu,\n"
        "  \"nonorientable_components\": %zu,\n"
        "  \"inconsistent_winding_components\": %zu,\n"
        "  \"all_components_are_disks\": %s,\n"
        "  \"minimal_generator_rank\": %zu,\n"
        "  \"emitted_generators\": %zu,\n"
        "  \"generators_truncated\": %s,\n"
        "  \"generator_basis_complete\": %s,\n"
        "  \"components\": [\n",
        report->betti_complete ? "true" : "false",
        report->disk_components, report->nondisk_components,
        report->invalid_surface_components,
        report->nonorientable_components,
        report->inconsistent_winding_components,
        report->all_components_are_disks ? "true" : "false",
        report->minimal_generator_rank,
        report->emitted_generators,
        report->generators_truncated ? "true" : "false",
        report->generator_basis_complete ? "true" : "false");

    for (size_t i = 0; i < report->face_components; i++) {
        const TopologyComponentInvariant *c = &report->component[i];
        fprintf(stream,
            "    {\"rank\": %zu, \"root_vertex\": %d, "
            "\"vertices\": %zu, \"edges\": %zu, \"faces\": %zu, "
            "\"bbox_min_zyx\": [%.17g, %.17g, %.17g], "
            "\"bbox_max_zyx\": [%.17g, %.17g, %.17g], "
            "\"centroid_zyx\": [%.17g, %.17g, %.17g], "
            "\"boundary_edges\": %zu, \"boundary_loops\": %zu, "
            "\"boundary_irregular_vertices\": %zu, "
            "\"nonmanifold_edges\": %zu, \"nonmanifold_vertices\": %zu, "
            "\"same_direction_edges\": %zu, "
            "\"euler_characteristic\": %lld, \"orientable\": %s, "
            "\"input_winding_consistent\": %s, \"surface_valid\": %s, "
            "\"betti\": [",
            c->rank, c->root_vertex, c->vertices, c->edges, c->faces,
            c->bbox_min[0], c->bbox_min[1], c->bbox_min[2],
            c->bbox_max[0], c->bbox_max[1], c->bbox_max[2],
            c->centroid[0], c->centroid[1], c->centroid[2],
            c->boundary_edges, c->boundary_loops,
            c->boundary_irregular_vertices, c->nonmanifold_edges,
            c->nonmanifold_vertices, c->same_direction_edges,
            (long long)c->euler_characteristic,
            c->orientable ? "true" : "false",
            c->input_winding_consistent ? "true" : "false",
            c->surface_valid ? "true" : "false");
        json_betti(stream, c->beta_0); fputs(", ", stream);
        json_betti(stream, c->beta_1); fputs(", ", stream);
        json_betti(stream, c->beta_2);
        fputs("], \"orientable_genus\": ", stream);
        json_betti(stream, c->orientable_genus);
        fputs(", \"crosscap_number\": ", stream);
        json_betti(stream, c->crosscap_number);
        fprintf(stream,
            ", \"homeomorphic_to_disk\": %s, \"defect_mask\": %u, "
            "\"defects\": ",
            c->homeomorphic_to_disk ? "true" : "false",
            (unsigned)c->defect_mask);
        json_defects(stream, c->defect_mask);
        fprintf(stream,
            ", \"minimal_generator_rank\": %zu, "
            "\"emitted_generators\": %zu, "
            "\"generator_basis_complete\": %s}%s\n",
            c->minimal_generator_rank, c->emitted_generators,
            c->generator_basis_complete ? "true" : "false",
            i + 1 == report->face_components ? "" : ",");
    }
    fputs("  ],\n  \"boundary_components\": [\n", stream);
    for (size_t i = 0; i < report->boundary_components; i++) {
        const TopologyBoundaryInvariant *b = &report->boundary[i];
        fprintf(stream,
            "    {\"rank\": %zu, \"component\": %zu, \"ordinal\": %zu, "
            "\"root_vertex\": %d, \"vertices\": %zu, \"edges\": %zu, "
            "\"irregular_vertices\": %zu, \"simple_cycle\": %s, "
            "\"component_perimeter\": %s, "
            "\"bbox_min_zyx\": [%.17g, %.17g, %.17g], "
            "\"bbox_max_zyx\": [%.17g, %.17g, %.17g], "
            "\"length\": %.17g, \"max_edge_length\": %.17g, "
            "\"diameter\": %.17g, \"diameter_exact\": %s}%s\n",
            b->rank, b->component, b->ordinal, b->root_vertex,
            b->vertices, b->edges, b->irregular_vertices,
            b->simple_cycle ? "true" : "false",
            b->component_perimeter ? "true" : "false",
            b->bbox_min[0], b->bbox_min[1], b->bbox_min[2],
            b->bbox_max[0], b->bbox_max[1], b->bbox_max[2],
            b->length, b->max_edge_length, b->diameter,
            b->diameter_exact ? "true" : "false",
            i + 1 == report->boundary_components ? "" : ",");
    }
    fputs("  ],\n  \"generators\": [\n", stream);
    for (size_t i = 0; i < report->emitted_generators; i++) {
        const TopologyGenerator *g = &report->generator[i];
        fprintf(stream,
            "    {\"component\": %zu, \"ordinal\": %zu, "
            "\"closing_edge\": [%d, %d], \"vertices\": %zu, "
            "\"length\": %.17g, \"rooted_length\": %.17g}%s\n",
            g->component, g->ordinal, g->closing_edge_a,
            g->closing_edge_b, g->nvertices, g->length,
            g->rooted_length,
            i + 1 == report->emitted_generators ? "" : ",");
    }
    fputs("  ]\n}\n", stream);
    return fclose(stream) == 0 ? 0 : -1;
}

static int write_generators_obj(const char *path, const float *verts,
                                const TopologyAuditReport *report)
{
    static const float palette[][3] = {
        {0.95f,0.20f,0.20f}, {0.20f,0.75f,1.00f},
        {0.25f,0.95f,0.35f}, {1.00f,0.70f,0.15f},
        {0.75f,0.30f,1.00f}, {0.05f,0.95f,0.85f},
        {1.00f,0.30f,0.70f}, {0.75f,0.95f,0.15f},
        {0.35f,0.45f,1.00f}, {1.00f,0.50f,0.25f},
        {0.30f,1.00f,0.65f}, {0.85f,0.35f,0.55f}
    };
    FILE *stream = fopen(path, "wb");
    size_t vertex_base = 1;
    if (!stream) return -1;
    fprintf(stream,
        "# topology_audit minimum-cardinality H1 generator representatives\n"
        "# %zu of %zu generators emitted; colored OBJ vertices plus line edges\n",
        report->emitted_generators, report->minimal_generator_rank);
    for (size_t i = 0; i < report->emitted_generators; i++) {
        const TopologyGenerator *g = &report->generator[i];
        const float *color =
            palette[(g->component * 5 + g->ordinal) %
                    (sizeof(palette) / sizeof(palette[0]))];
        fprintf(stream,
            "o component_%zu_generator_%zu\n"
            "# simple_length %.17g rooted_length %.17g\n",
            g->component, g->ordinal, g->length, g->rooted_length);
        for (size_t k = 0; k < g->nvertices; k++) {
            int32_t v = g->vertices[k];
            fprintf(stream, "v %.9g %.9g %.9g %.3f %.3f %.3f\n",
                    verts[(size_t)v*3+0], verts[(size_t)v*3+1],
                    verts[(size_t)v*3+2],
                    color[0], color[1], color[2]);
        }
        for (size_t k = 1; k < g->nvertices; k++)
            fprintf(stream, "l %zu %zu\n",
                    vertex_base + k - 1, vertex_base + k);
        if (g->nvertices > 1)
            fprintf(stream, "l %zu %zu\n",
                    vertex_base + g->nvertices - 1, vertex_base);
        vertex_base += g->nvertices;
    }
    return fclose(stream) == 0 ? 0 : -1;
}

static void print_report(const TopologyAuditReport *report, size_t top)
{
    printf("topology_audit: V=%zu E=%zu F=%zu chi=%lld\n",
           report->vertices, report->edges, report->faces,
           (long long)report->euler_characteristic);
    printf("  components=%zu face + %zu isolated = %zu, Betti_F2=(",
           report->face_components, report->isolated_vertices,
           report->topological_components);
    print_betti(stdout, report->beta_0); fputs(", ", stdout);
    print_betti(stdout, report->beta_1); fputs(", ", stdout);
    print_betti(stdout, report->beta_2); fputs(")\n", stdout);
    printf("  disks=%zu nondisks=%zu invalid_surfaces=%zu "
           "nonorientable=%zu inconsistent_winding=%zu\n",
           report->disk_components, report->nondisk_components,
           report->invalid_surface_components,
           report->nonorientable_components,
           report->inconsistent_winding_components);
    printf("  H1 minimum generators=%zu, basis=%s, representatives=%zu%s\n",
           report->minimal_generator_rank,
           report->generator_basis_complete ? "verified" : "INCOMPLETE",
           report->emitted_generators,
           report->generators_truncated ? " (truncated)" : "");
    printf("  disk certificate: %s\n",
           report->all_components_are_disks ? "PASS" : "FAIL");
    if (top > report->face_components) top = report->face_components;
    if (top == 0) return;
    puts("  rank       faces       verts       edges   bnd  chi   "
         "g/k   Betti_F2    disk  defects");
    for (size_t i = 0; i < top; i++) {
        const TopologyComponentInvariant *c = &report->component[i];
        char genus[32];
        if (c->orientable_genus >= 0)
            snprintf(genus, sizeof(genus), "g%lld",
                     (long long)c->orientable_genus);
        else if (c->crosscap_number >= 0)
            snprintf(genus, sizeof(genus), "k%lld",
                     (long long)c->crosscap_number);
        else strcpy(genus, "?");
        printf("  %4zu %11zu %11zu %11zu %5zu %4lld %5s   (",
               c->rank, c->faces, c->vertices, c->edges,
               c->boundary_loops, (long long)c->euler_characteristic,
               genus);
        print_betti(stdout, c->beta_0); fputc(',', stdout);
        print_betti(stdout, c->beta_1); fputc(',', stdout);
        print_betti(stdout, c->beta_2);
        printf(")   %s   ", c->homeomorphic_to_disk ? "yes" : "NO ");
        print_defects(stdout, c->defect_mask);
        fputc('\n', stdout);
    }
}

int main(int argc, char **argv)
{
    CliOptions cli;
    TopologyAuditOptions audit_options;
    TopologyAuditReport report;
    Arena_T arena;
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    int parsed, gate_failed = 0;

    if (argc == 2 && !strcmp(argv[1], "--selftest")) {
        int failures = TopologyAudit_selftest();
        failures += SeamCut_selftest();
        failures += DiskTopologyRepair_selftest();
        return failures ? 3 : 0;
    }
    if (argc == 2 &&
        (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h"))) {
        usage(stdout);
        return 0;
    }
    parsed = parse_cli(argc, argv, &cli);
    if (parsed != 0) {
        usage(parsed > 0 ? stdout : stderr);
        return parsed > 0 ? 0 : 2;
    }

    arena = Arena_new();
    if (ObjIO_read(arena, cli.input, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "topology_audit: cannot read %s\n", cli.input);
        Arena_dispose(&arena);
        return 1;
    }
    TopologyAudit_options_default(&audit_options);
    audit_options.emit_generators = cli.generators != NULL;
    audit_options.max_generators_per_component =
        cli.max_generators_per_component;
    if (TopologyAudit_analyze(verts, nv, faces, nf,
                              &audit_options, &report) != 0) {
        fprintf(stderr, "topology_audit: %s\n", report.error);
        Arena_dispose(&arena);
        return 1;
    }

    print_report(&report, cli.top);
    if (cli.json) {
        if (write_json(cli.json, cli.input, &report) != 0) {
            fprintf(stderr, "topology_audit: cannot write %s\n", cli.json);
            TopologyAudit_dispose(&report);
            Arena_dispose(&arena);
            return 1;
        }
        printf("  wrote %s\n", cli.json);
    }
    if (cli.generators) {
        if (write_generators_obj(cli.generators, verts, &report) != 0) {
            fprintf(stderr, "topology_audit: cannot write %s\n",
                    cli.generators);
            TopologyAudit_dispose(&report);
            Arena_dispose(&arena);
            return 1;
        }
        printf("  wrote %s\n", cli.generators);
    }
    if (cli.repair_disk) {
        float *repair_verts = NULL;
        int32_t *repair_faces = NULL, *repair_source = NULL;
        size_t repair_nv = 0, repair_nf = 0;
        DiskTopologyRepairStats stats;
        TopologyAuditReport repaired;
        if (DiskTopologyRepair_process(arena, verts, nv, faces, nf,
                                       &repair_verts, &repair_nv,
                                       &repair_faces, &repair_nf,
                                       &repair_source, &stats) != 0) {
            fprintf(stderr,
                    "topology_audit: disk repair failed after %zu/%zu "
                    "component(s)\n",
                    stats.repaired_components + stats.already_disks,
                    stats.components);
            TopologyAudit_dispose(&report);
            Arena_dispose(&arena);
            return 1;
        }
        (void)repair_source;
        if (TopologyAudit_analyze(repair_verts, repair_nv,
                                  repair_faces, repair_nf,
                                  NULL, &repaired) != 0 ||
            !repaired.all_components_are_disks) {
            fprintf(stderr,
                    "topology_audit: repaired mesh failed exact disk gate%s%s\n",
                    repaired.error[0] ? ": " : "",
                    repaired.error[0] ? repaired.error : "");
            TopologyAudit_dispose(&repaired);
            TopologyAudit_dispose(&report);
            Arena_dispose(&arena);
            return 1;
        }
        if (ObjIO_write(cli.repair_disk, repair_verts, repair_nv,
                        repair_faces, repair_nf) != 0) {
            fprintf(stderr, "topology_audit: cannot write %s\n",
                    cli.repair_disk);
            TopologyAudit_dispose(&repaired);
            TopologyAudit_dispose(&report);
            Arena_dispose(&arena);
            return 1;
        }
        printf("  disk repair: %zu already + %zu repaired; "
               "loops %zu -> %zu, handles=%zu, seam_edges=%zu, "
               "V %zu -> %zu, F %zu -> %zu\n",
               stats.already_disks, stats.repaired_components,
               stats.boundary_loops_before, stats.boundary_loops_after,
               stats.handles_opened, stats.seam_edges,
               stats.vertices_in, stats.vertices_out,
               stats.faces_in, stats.faces_out);
        printf("  exact repaired disk certificate: PASS (%zu/%zu)\n",
               repaired.disk_components, repaired.face_components);
        printf("  wrote %s\n", cli.repair_disk);
        TopologyAudit_dispose(&repaired);
    }

    if (cli.fail_nondisk && !report.all_components_are_disks)
        gate_failed = 1;
    if (cli.require_consistent_winding &&
        report.inconsistent_winding_components)
        gate_failed = 1;
    TopologyAudit_dispose(&report);
    Arena_dispose(&arena);
    return gate_failed ? 4 : 0;
}
