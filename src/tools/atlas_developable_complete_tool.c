/*
 * Standalone checkpoint-to-completion experiment.
 *
 * The executable consumes placed geometry plus an authoritative
 * atlas_solution.bin.  It never reads a ribbon OBJ.  The trusted evidence
 * core and all completion patches therefore retain checkpoint chart, face,
 * winding, rank, and collision provenance in one process.
 */

#include "../common/arena.h"
#include "../common/ves_platform.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"
#include "../whole/atlas_developable_complete.h"
#include "../whole/atlas_ribbon_fit.h"
#include "../whole/atlas_solution.h"
#include "../whole/atlas_track_export.h"
#include "../whole/atlas_track_grow.h"

#include <ctype.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ADCT_PATH_CAP 2048

static void adct_usage(const char *program)
{
    fprintf(stderr,
        "usage:\n"
        "  %s --selftest\n"
        "  %s --fingerprint <placed_dir>\n"
        "  %s <placed_dir> <atlas_solution.bin> <out_dir> <prefix> "
        "[options]\n\n"
        "Evidence/core:\n"
        "  --slice-spacing F          constant-v evidence rows (default 2)\n"
        "  --observation-u F          checkpoint-U sampling (default 2)\n"
        "  --local-xyz F              duplicate cluster radius (default 8)\n"
        "  --min-component-quads N    trusted core threshold (default 96)\n"
        "  --gap-fill-rounds N        evidence-only frontier rounds (default 8)\n"
        "  --parameter-iters N        shared-atlas coupling sweeps (default 64)\n\n"
        "Developable completion:\n"
        "  --min-gap-u F              smallest bridge span (default 4)\n"
        "  --max-gap-u F              largest bridge span (default 96)\n"
        "  --min-gap-v F              smallest axial bridge span (default 4)\n"
        "  --max-gap-v F              largest axial bridge span (default 96)\n"
        "  --min-patch-rows N         coherent rows required (default 4)\n"
        "  --bank-shift-columns N     boundary drift per row (default 6)\n"
        "  --sample-columns N         fixed patch columns; 0=atlas spacing\n"
        "  --fill-internal-gaps       heal bounded holes inside a component\n"
        "  --fill-vertical-gaps       heal gaps across evidence rows\n"
        "  --allow-component-drift    track one bank across fragment IDs\n"
        "  --shared-lattice-fill      union accepted U/V fills on one grid\n"
        "  --fill-aligned-gaps        close aligned U/V cells and short strips\n"
        "  --short-strip-max-span N   set both bracketed fallback limits\n"
        "  --short-strip-max-u N      missing-U fallback limit (default 8)\n"
        "  --short-strip-max-v N      missing-V fallback limit (default 8)\n"
        "  --normal-cone-dot F        local height-chart gate (default .25)\n"
        "  --bank-fit-rms F           tangential bank RMS gate (default 2)\n"
        "  --bank-fit-max F           tangential bank max gate (default 6)\n"
        "  --hermite-weight F         evidence-tangent prior (default .2)\n"
        "  --nuclear-lambda F         convex warm-start weight (default .05)\n"
        "  --rank-tail-lambda F       smaller-curvature weight (default .2)\n"
        "  --admm-rho F               ADMM penalty (default 1)\n"
        "  --nuclear-iters N          convex iterations (default 8)\n"
        "  --rank-tail-iters N        rank-one iterations (default 24)\n"
        "  --pcg-iters N              inner solve cap (default 160)\n"
        "  --max-patches N            deterministic patch cap (default 4096)\n\n"
        "Export:\n"
        "  --no-winding-export        skip machine winding organizer\n"
        "  --winding-wraps N          windings per piece (default 1)\n"
        "  --winding-slab F           axial slab height (default 4096)\n"
        "  --winding-du F             organizer U vox/px (default 1)\n"
        "  --winding-dv F             organizer V vox/px (default 1)\n",
        program, program, program);
}

static int adct_valid_prefix(const char *prefix)
{
    if (prefix == NULL || *prefix == '\0') return 0;
    for (const unsigned char *p = (const unsigned char *)prefix; *p; p++)
        if (!isalnum(*p) && *p != '_' && *p != '-' && *p != '.')
            return 0;
    return 1;
}

static int adct_parse_integer(const char *text, int *out)
{
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (end == text || *end != '\0' ||
        value < INT32_MIN || value > INT32_MAX)
        return -1;
    *out = (int)value;
    return 0;
}

static int adct_parse_real(const char *text, double *out)
{
    char *end = NULL;
    double value = strtod(text, &end);
    if (end == text || *end != '\0' || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

static FILE *adct_open_output(
    const char *directory,
    const char *prefix,
    const char *suffix)
{
    char path[ADCT_PATH_CAP];
    int count = snprintf(
        path, sizeof path, "%s/%s%s", directory, prefix, suffix);
    if (count < 0 || (size_t)count >= sizeof path ||
        ves_ensure_parent_dir(path) != 0)
        return NULL;
    return fopen(path, "wb");
}

static int adct_close_output(FILE *file)
{
    int failed = ferror(file);
    if (fclose(file) != 0) failed = 1;
    return failed ? -1 : 0;
}

static void adct_json_string(FILE *file, const char *text)
{
    fputc('"', file);
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
        if (*p == '"' || *p == '\\') {
            fputc('\\', file);
            fputc((int)*p, file);
        } else if (*p == '\n') {
            fputs("\\n", file);
        } else if (*p == '\r') {
            fputs("\\r", file);
        } else if (*p == '\t') {
            fputs("\\t", file);
        } else if (*p >= 0x20) {
            fputc((int)*p, file);
        }
    }
    fputc('"', file);
}

static int adct_write_surface_obj(
    const char *directory,
    const char *prefix,
    const char *suffix,
    const AtlasTrackGrowResult *mesh,
    size_t first_face)
{
    FILE *file = adct_open_output(directory, prefix, suffix);
    if (file == NULL) return -1;
    fprintf(file,
        "# Direct atlas-checkpoint developable completion.\n"
        "# Core vertices are green; inferred bridge vertices are magenta.\n"
        "# Coordinates remain source-space (z,y,x); vt is the shared atlas.\n");
    for (size_t i = 0; i < mesh->nv; i++) {
        const double *point = &mesh->xyz[i * 3];
        if (mesh->source_target[i] < 0)
            fprintf(file, "v %.17g %.17g %.17g 0.95 0.24 0.72\n",
                    point[0], point[1], point[2]);
        else
            fprintf(file, "v %.17g %.17g %.17g 0.16 0.82 0.44\n",
                    point[0], point[1], point[2]);
    }
    for (size_t i = 0; i < mesh->nv; i++)
        fprintf(file, "vt %.17g %.17g\n",
                mesh->uv[i * 2], mesh->uv[i * 2 + 1]);
    fprintf(file, "s 1\n");
    int previous = -1;
    for (size_t face = first_face; face < mesh->nf; face++) {
        int inferred = mesh->face_inferred[face] ? 1 : 0;
        if (inferred != previous) {
            fprintf(file, "g %s\n",
                    inferred ? "developable_completion" : "evidence_core");
            previous = inferred;
        }
        int32_t a = mesh->faces[face * 3] + 1;
        int32_t b = mesh->faces[face * 3 + 1] + 1;
        int32_t c = mesh->faces[face * 3 + 2] + 1;
        fprintf(file, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    return adct_close_output(file);
}

static int adct_write_atlas_obj(
    const char *directory,
    const char *prefix,
    const AtlasTrackGrowResult *mesh)
{
    FILE *file = adct_open_output(
        directory, prefix, "_developable_atlas.obj");
    if (file == NULL) return -1;
    fprintf(file,
        "# Planar view in the checkpoint's shared packed atlas gauge.\n"
        "# Core is blue; developable completion is magenta.\n");
    for (size_t i = 0; i < mesh->nv; i++) {
        if (mesh->source_target[i] < 0)
            fprintf(file, "v %.17g %.17g 0 0.95 0.24 0.72\n",
                    mesh->uv[i * 2], mesh->uv[i * 2 + 1]);
        else
            fprintf(file, "v %.17g %.17g 0 0.20 0.64 0.96\n",
                    mesh->uv[i * 2], mesh->uv[i * 2 + 1]);
    }
    int previous = -1;
    for (size_t face = 0; face < mesh->nf; face++) {
        int inferred = mesh->face_inferred[face] ? 1 : 0;
        if (inferred != previous) {
            fprintf(file, "g %s\n",
                    inferred ? "developable_completion" : "evidence_core");
            previous = inferred;
        }
        fprintf(file, "f %d %d %d\n",
                mesh->faces[face * 3] + 1,
                mesh->faces[face * 3 + 1] + 1,
                mesh->faces[face * 3 + 2] + 1);
    }
    return adct_close_output(file);
}

static int adct_write_patches_csv(
    const char *directory,
    const char *prefix,
    const AtlasDevelopableCompleteResult *result)
{
    FILE *file = adct_open_output(
        directory, prefix, "_developable_patches.csv");
    if (file == NULL) return -1;
    fprintf(file,
        "patch,direction,accepted,reject_reason,component_left,component_right,"
        "row_first,row_last,rows,columns,column_left_min,column_left_max,"
        "column_right_min,column_right_max,added_vertices,added_faces,"
        "component_pair_changes,"
        "gap_span_mean,gap_span_max,chord_ratio_mean,bank_fit_rms,bank_fit_max,"
        "tail_before,tail_after,primal,dual,pcg_iterations\n");
    for (size_t i = 0; i < result->npatches; i++) {
        const AtlasDevelopablePatch *patch = &result->patch[i];
        fprintf(file,
            "%d,%d,%d,%d,%d,%d,%d,%d,%zu,%zu,%d,%d,%d,%d,%zu,%zu,%zu,"
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",
            patch->id, patch->direction,
            patch->accepted, patch->reject_reason,
            patch->component_left, patch->component_right,
            patch->row_first, patch->row_last, patch->rows, patch->columns,
            patch->column_left_min, patch->column_left_max,
            patch->column_right_min, patch->column_right_max,
            patch->added_vertices, patch->added_faces,
            patch->component_pair_changes,
            patch->gap_u_mean, patch->gap_u_max, patch->chord_ratio_mean,
            patch->bank_fit_rms, patch->bank_fit_max,
            patch->tail_before, patch->tail_after,
            patch->primal_residual, patch->dual_residual,
            patch->pcg_iterations);
    }
    return adct_close_output(file);
}

static int adct_write_stats(
    const char *directory,
    const char *prefix,
    const char *placed_directory,
    const char *solution_path,
    const PieceSet *pieces,
    size_t kept_faces,
    const AtlasRibbonFitOptions *evidence_options,
    const AtlasTrackGrowOptions *grow_options,
    const AtlasDevelopableCompleteOptions *complete_options,
    const AtlasRibbonObservationSet *evidence,
    const AtlasTrackGrowResult *seed,
    const AtlasDevelopableCompleteResult *result,
    const AtlasTrackExportStats *winding,
    const double timing[5])
{
    FILE *file = adct_open_output(
        directory, prefix, "_developable_stats.json");
    if (file == NULL) return -1;
    fputs("{\n  \"schema\": \"atlas_developable_complete_v1\",\n"
          "  \"placed_dir\": ", file);
    adct_json_string(file, placed_directory);
    fputs(",\n  \"solution_path\": ", file);
    adct_json_string(file, solution_path);
    fprintf(file,
        ",\n  \"input\": {\"cubes\": %zu, \"vertices\": %zu, "
        "\"faces\": %zu, \"kept_faces\": %zu},\n"
        "  \"evidence_options\": {\"slice_spacing\": %.17g, "
        "\"observation_u_spacing\": %.17g, \"local_xyz_tolerance\": %.17g},\n"
        "  \"grow_options\": {\"min_component_quads\": %d, "
        "\"gap_fill_rounds\": %d, \"uv_bridge_max_u\": %.17g, "
        "\"parameter_iterations\": %d},\n"
        "  \"completion_options\": {\"min_gap_u\": %.17g, "
        "\"max_gap_u\": %.17g, \"min_gap_v\": %.17g, "
        "\"max_gap_v\": %.17g, \"min_patch_rows\": %d, "
        "\"max_bank_shift_columns\": %d, \"sample_columns\": %d, "
        "\"fill_internal_gaps\": %d, "
        "\"fill_vertical_gaps\": %d, "
        "\"allow_component_drift\": %d, "
        "\"shared_lattice_fill\": %d, "
        "\"fill_aligned_gaps\": %d, "
        "\"short_strip_max_span_u\": %d, "
        "\"short_strip_max_span_v\": %d, "
        "\"normal_cone_dot_min\": %.17g, \"bank_fit_rms_max\": %.17g, "
        "\"bank_fit_max_max\": %.17g, \"hermite_prior_weight\": %.17g, "
        "\"nuclear_lambda\": %.17g, \"rank_tail_lambda\": %.17g, "
        "\"admm_rho\": %.17g, \"nuclear_iterations\": %d, "
        "\"rank_tail_iterations\": %d, \"pcg_iterations\": %d, "
        "\"max_patches\": %d},\n",
        pieces->n_cubes, pieces->nv, pieces->nf, kept_faces,
        evidence_options->slice_spacing,
        evidence_options->observation_u_spacing,
        evidence_options->local_xyz_tolerance,
        grow_options->min_component_quads,
        grow_options->gap_fill_rounds,
        grow_options->uv_bridge_max_u,
        grow_options->parameter_iterations,
        complete_options->min_gap_u, complete_options->max_gap_u,
        complete_options->min_gap_v, complete_options->max_gap_v,
        complete_options->min_patch_rows,
        complete_options->max_bank_shift_columns,
        complete_options->sample_columns,
        complete_options->fill_internal_gaps,
        complete_options->fill_vertical_gaps,
        complete_options->allow_component_drift,
        complete_options->shared_lattice_fill,
        complete_options->fill_aligned_gaps,
        complete_options->short_strip_max_span_u,
        complete_options->short_strip_max_span_v,
        complete_options->normal_cone_dot_min,
        complete_options->bank_fit_rms_max,
        complete_options->bank_fit_max_max,
        complete_options->hermite_prior_weight,
        complete_options->nuclear_lambda,
        complete_options->rank_tail_lambda,
        complete_options->admm_rho,
        complete_options->nuclear_iterations,
        complete_options->rank_tail_iterations,
        complete_options->pcg_iterations,
        complete_options->max_patches);
    fprintf(file,
        "  \"evidence\": {\"observations\": %zu, \"targets\": %zu, "
        "\"accepted_targets\": %zu, \"conflict_targets\": %zu, "
        "\"rows\": %zu},\n"
        "  \"seed\": {\"components\": %zu, \"vertices\": %zu, "
        "\"faces\": %zu, \"core_quads\": %zu, \"inferred_quads\": %zu},\n"
        "  \"completion\": {\"candidate_rows\": %zu, "
        "\"candidate_pairs\": %zu, \"patches\": %zu, "
        "\"accepted_patches\": %zu, \"rejected_patches\": %zu, "
        "\"components_before\": %zu, \"components_after\": %zu, "
        "\"added_vertices\": %zu, \"added_faces\": %zu, "
        "\"lattice_nodes\": %zu, \"lattice_cells\": %zu, "
        "\"lattice_conflicts\": %zu, "
        "\"lattice_rejected_cells\": %zu, "
        "\"aligned_gap_lines\": %zu, \"aligned_gap_cells\": %zu, "
        "\"adjacent_gap_cells\": %zu, "
        "\"adjacent_gap_u_cells\": %zu, "
        "\"adjacent_gap_v_cells\": %zu, "
        "\"short_strip_patches\": %zu, "
        "\"short_strip_cells\": %zu, "
        "\"short_strip_u_cells\": %zu, "
        "\"short_strip_v_cells\": %zu, "
        "\"vertices\": %zu, \"faces\": %zu, "
        "\"tail_before\": %.17g, \"tail_after\": %.17g, "
        "\"primal_residual_max\": %.17g, \"dual_residual_max\": %.17g, "
        "\"pcg_iterations\": %zu},\n",
        evidence->nobservation, evidence->ntarget,
        evidence->accepted_targets, evidence->conflict_targets,
        evidence->nrows,
        seed->selected_components, seed->nv, seed->nf,
        seed->selected_core_quads, seed->selected_inferred_quads,
        result->candidate_rows, result->candidate_pairs, result->npatches,
        result->accepted_patches, result->rejected_patches,
        result->seed_components, result->completed_components,
        result->added_vertices, result->added_faces,
        result->lattice_nodes, result->lattice_cells,
        result->lattice_conflicts, result->lattice_rejected_cells,
        result->aligned_gap_lines, result->aligned_gap_cells,
        result->adjacent_gap_cells,
        result->adjacent_gap_u_cells, result->adjacent_gap_v_cells,
        result->short_strip_patches, result->short_strip_cells,
        result->short_strip_u_cells, result->short_strip_v_cells,
        result->mesh.nv, result->mesh.nf,
        result->tail_before, result->tail_after,
        result->primal_residual_max, result->dual_residual_max,
        result->pcg_iterations);
    fprintf(file,
        "  \"winding_export\": {\"pieces\": %zu, \"written\": %zu, "
        "\"valid_pixels\": %zu, \"conflict_pixels\": %zu, "
        "\"conflict_fraction\": %.17g},\n"
        "  \"timing_seconds\": {\"load\": %.6f, \"evidence\": %.6f, "
        "\"core\": %.6f, \"completion\": %.6f, \"artifacts\": %.6f, "
        "\"total\": %.6f}\n}\n",
        winding->pieces, winding->written, winding->valid_pixels,
        winding->conflict_pixels, winding->conflict_fraction,
        timing[0], timing[1], timing[2], timing[3], timing[4],
        timing[0] + timing[1] + timing[2] + timing[3] + timing[4]);
    return adct_close_output(file);
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return AtlasDevelopableComplete_selftest() == 0 ? 0 : 1;
    if (argc == 3 && strcmp(argv[1], "--fingerprint") == 0) {
        Arena_T probe_arena = Arena_new();
        PieceSet probe_pieces;
        if (probe_arena == NULL ||
            PieceSet_build(probe_arena, argv[2], &probe_pieces) != 0) {
            fprintf(stderr, "%s: cannot load placed directory %s\n",
                    argv[0], argv[2]);
            if (probe_arena != NULL) Arena_dispose(&probe_arena);
            return 1;
        }
        printf("cubes=%zu vertices=%zu faces=%zu fingerprint=%016" PRIx64
               " placed=%s\n",
               probe_pieces.n_cubes, probe_pieces.nv, probe_pieces.nf,
               AtlasSolution_piece_fingerprint(&probe_pieces), argv[2]);
        Arena_dispose(&probe_arena);
        return 0;
    }
    if (argc < 5) {
        adct_usage(argv[0]);
        return 1;
    }
    const char *placed_directory = argv[1];
    const char *solution_path = argv[2];
    const char *output_directory = argv[3];
    const char *prefix = argv[4];
    if (!adct_valid_prefix(prefix)) {
        fprintf(stderr, "%s: invalid prefix\n", argv[0]);
        return 1;
    }

    AtlasRibbonFitOptions evidence_options;
    AtlasRibbonFitOptions_default(&evidence_options);
    evidence_options.slice_spacing = 2.0;
    evidence_options.observation_u_spacing = 2.0;
    AtlasTrackGrowOptions grow_options;
    AtlasTrackGrowOptions_default(&grow_options);
    grow_options.min_component_quads = 96;
    grow_options.gap_fill_rounds = 8;
    grow_options.uv_bridge_max_u = 0.0;
    grow_options.parameter_iterations = 64;
    evidence_options.tangent_dot_min = grow_options.tangent_dot_min;
    AtlasDevelopableCompleteOptions complete_options;
    AtlasDevelopableCompleteOptions_default(&complete_options);
    int winding_export = 1;
    int winding_wraps = 1;
    double winding_slab = 4096.0;
    double winding_du = 1.0, winding_dv = 1.0;

    for (int i = 5; i < argc; i++) {
        const char *option = argv[i];
        if (strcmp(option, "--no-winding-export") == 0) {
            winding_export = 0;
            continue;
        }
        if (strcmp(option, "--fill-internal-gaps") == 0) {
            complete_options.fill_internal_gaps = 1;
            continue;
        }
        if (strcmp(option, "--fill-vertical-gaps") == 0) {
            complete_options.fill_vertical_gaps = 1;
            continue;
        }
        if (strcmp(option, "--allow-component-drift") == 0) {
            complete_options.allow_component_drift = 1;
            continue;
        }
        if (strcmp(option, "--shared-lattice-fill") == 0) {
            complete_options.shared_lattice_fill = 1;
            continue;
        }
        if (strcmp(option, "--fill-aligned-gaps") == 0) {
            complete_options.shared_lattice_fill = 1;
            complete_options.fill_aligned_gaps = 1;
            continue;
        }
        if (strcmp(option, "--short-strip-max-span") == 0) {
            int limit = 0;
            if (i + 1 >= argc ||
                adct_parse_integer(argv[++i], &limit) != 0) {
                fprintf(stderr, "%s: invalid %s\n", argv[0], option);
                return 1;
            }
            complete_options.short_strip_max_span_u = limit;
            complete_options.short_strip_max_span_v = limit;
            continue;
        }
        int *integer = NULL;
        if (strcmp(option, "--min-component-quads") == 0)
            integer = &grow_options.min_component_quads;
        else if (strcmp(option, "--gap-fill-rounds") == 0)
            integer = &grow_options.gap_fill_rounds;
        else if (strcmp(option, "--parameter-iters") == 0)
            integer = &grow_options.parameter_iterations;
        else if (strcmp(option, "--min-patch-rows") == 0)
            integer = &complete_options.min_patch_rows;
        else if (strcmp(option, "--bank-shift-columns") == 0)
            integer = &complete_options.max_bank_shift_columns;
        else if (strcmp(option, "--sample-columns") == 0)
            integer = &complete_options.sample_columns;
        else if (strcmp(option, "--short-strip-max-u") == 0)
            integer = &complete_options.short_strip_max_span_u;
        else if (strcmp(option, "--short-strip-max-v") == 0)
            integer = &complete_options.short_strip_max_span_v;
        else if (strcmp(option, "--nuclear-iters") == 0)
            integer = &complete_options.nuclear_iterations;
        else if (strcmp(option, "--rank-tail-iters") == 0)
            integer = &complete_options.rank_tail_iterations;
        else if (strcmp(option, "--pcg-iters") == 0)
            integer = &complete_options.pcg_iterations;
        else if (strcmp(option, "--max-patches") == 0)
            integer = &complete_options.max_patches;
        else if (strcmp(option, "--winding-wraps") == 0)
            integer = &winding_wraps;
        if (integer != NULL) {
            if (i + 1 >= argc || adct_parse_integer(argv[++i], integer) != 0) {
                fprintf(stderr, "%s: invalid %s\n", argv[0], option);
                return 1;
            }
            continue;
        }

        double *real = NULL;
        if (strcmp(option, "--slice-spacing") == 0)
            real = &evidence_options.slice_spacing;
        else if (strcmp(option, "--observation-u") == 0)
            real = &evidence_options.observation_u_spacing;
        else if (strcmp(option, "--local-xyz") == 0)
            real = &evidence_options.local_xyz_tolerance;
        else if (strcmp(option, "--min-gap-u") == 0)
            real = &complete_options.min_gap_u;
        else if (strcmp(option, "--max-gap-u") == 0)
            real = &complete_options.max_gap_u;
        else if (strcmp(option, "--min-gap-v") == 0)
            real = &complete_options.min_gap_v;
        else if (strcmp(option, "--max-gap-v") == 0)
            real = &complete_options.max_gap_v;
        else if (strcmp(option, "--normal-cone-dot") == 0)
            real = &complete_options.normal_cone_dot_min;
        else if (strcmp(option, "--bank-fit-rms") == 0)
            real = &complete_options.bank_fit_rms_max;
        else if (strcmp(option, "--bank-fit-max") == 0)
            real = &complete_options.bank_fit_max_max;
        else if (strcmp(option, "--hermite-weight") == 0)
            real = &complete_options.hermite_prior_weight;
        else if (strcmp(option, "--nuclear-lambda") == 0)
            real = &complete_options.nuclear_lambda;
        else if (strcmp(option, "--rank-tail-lambda") == 0)
            real = &complete_options.rank_tail_lambda;
        else if (strcmp(option, "--admm-rho") == 0)
            real = &complete_options.admm_rho;
        else if (strcmp(option, "--winding-slab") == 0)
            real = &winding_slab;
        else if (strcmp(option, "--winding-du") == 0)
            real = &winding_du;
        else if (strcmp(option, "--winding-dv") == 0)
            real = &winding_dv;
        if (real == NULL || i + 1 >= argc ||
            adct_parse_real(argv[++i], real) != 0) {
            fprintf(stderr, "%s: unknown or invalid option %s\n",
                    argv[0], option);
            return 1;
        }
    }
    if (winding_wraps < 1 || !(winding_slab > 0.0) ||
        !(winding_du > 0.0) || !(winding_dv > 0.0)) {
        fprintf(stderr, "%s: invalid winding export dimensions\n", argv[0]);
        return 1;
    }

    char probe[ADCT_PATH_CAP];
    int probe_count = snprintf(
        probe, sizeof probe, "%s/%s.probe", output_directory, prefix);
    if (probe_count < 0 || (size_t)probe_count >= sizeof probe ||
        ves_ensure_parent_dir(probe) != 0) {
        fprintf(stderr, "%s: cannot create output directory %s\n",
                argv[0], output_directory);
        return 1;
    }

    double start = ves_clock_sec();
    Arena_T arena = Arena_new();
    if (arena == NULL) return 1;
    PieceSet pieces;
    if (PieceSet_build(arena, placed_directory, &pieces) != 0) {
        fprintf(stderr, "%s: cannot load placed directory %s\n",
                argv[0], placed_directory);
        Arena_dispose(&arena);
        return 1;
    }
    ScaffoldCalib calibration;
    if (Scaffold_read_calib(placed_directory, &calibration) != 0) {
        fprintf(stderr, "%s: missing placed calibration in %s\n",
                argv[0], placed_directory);
        Arena_dispose(&arena);
        return 1;
    }
    AtlasSolution solution;
    if (AtlasSolution_read(
            arena, solution_path, &pieces, &solution) != 0) {
        fprintf(stderr,
                "%s: atlas checkpoint does not match placed input: %s\n",
                argv[0], solution_path);
        Arena_dispose(&arena);
        return 1;
    }
    int wind_ranked = AtlasSolution_rank_by_wind(arena, &solution);
    if (wind_ranked < 0) {
        fprintf(stderr, "%s: invalid checkpoint wind field\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    size_t kept_faces = 0;
    for (size_t face = 0; face < solution.nfaces; face++)
        if (solution.face_keep[face]) kept_faces++;
    double loaded = ves_clock_sec();
    fprintf(stderr,
        "[atlas_developable_complete] checkpoint cubes=%zu vertices=%zu "
        "faces=%zu kept=%zu charts=%zu ranks=%s fingerprint=%016" PRIx64 "\n",
        pieces.n_cubes, pieces.nv, pieces.nf, kept_faces, solution.ncharts,
        wind_ranked ? "wind" : "producer",
        AtlasSolution_piece_fingerprint(&pieces));

    AtlasRibbonObservationSet evidence;
    if (AtlasRibbonFit_build_observations(
            arena, &pieces, &solution, &calibration,
            &evidence_options, &evidence) != 0) {
        fprintf(stderr, "%s: evidence extraction failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    double extracted = ves_clock_sec();
    fprintf(stderr,
        "[atlas_developable_complete] evidence observations=%zu targets=%zu "
        "accepted=%zu conflicts=%zu rows=%zu\n",
        evidence.nobservation, evidence.ntarget, evidence.accepted_targets,
        evidence.conflict_targets, evidence.nrows);

    AtlasTrackGrowResult seed;
    if (AtlasTrackGrow_build(
            arena, &evidence, &calibration, &grow_options, &seed) != 0) {
        fprintf(stderr, "%s: trusted evidence core failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    double core_done = ves_clock_sec();
    fprintf(stderr,
        "[atlas_developable_complete] core components=%zu vertices=%zu "
        "faces=%zu core_quads=%zu inferred_quads=%zu\n",
        seed.selected_components, seed.nv, seed.nf,
        seed.selected_core_quads, seed.selected_inferred_quads);

    AtlasDevelopableCompleteResult result;
    if (AtlasDevelopableComplete_build(
            arena, &evidence, &seed, &complete_options, &result) != 0) {
        fprintf(stderr, "%s: developable completion failed\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    double completion_done = ves_clock_sec();
    fprintf(stderr,
        "[atlas_developable_complete] patches=%zu/%zu rows=%zu pairs=%zu "
        "added=%zu vertices/%zu faces adjacent=%zu (u=%zu v=%zu) "
        "short=%zu/%zu (u=%zu v=%zu) "
        "components=%zu->%zu "
        "rank_tail=%.6g->%.6g pcg=%zu\n",
        result.accepted_patches, result.npatches,
        result.candidate_rows, result.candidate_pairs,
        result.added_vertices, result.added_faces,
        result.adjacent_gap_cells,
        result.adjacent_gap_u_cells, result.adjacent_gap_v_cells,
        result.short_strip_patches, result.short_strip_cells,
        result.short_strip_u_cells, result.short_strip_v_cells,
        result.seed_components, result.completed_components,
        result.tail_before, result.tail_after, result.pcg_iterations);

    int io = 0;
    io |= adct_write_surface_obj(
        output_directory, prefix, "_developable_seed.obj", &seed, 0);
    io |= adct_write_surface_obj(
        output_directory, prefix, "_developable.obj", &result.mesh, 0);
    io |= adct_write_surface_obj(
        output_directory, prefix, "_developable_completion_only.obj",
        &result.mesh, seed.nf);
    io |= adct_write_atlas_obj(output_directory, prefix, &result.mesh);
    io |= adct_write_patches_csv(output_directory, prefix, &result);

    AtlasTrackExportStats winding_stats;
    memset(&winding_stats, 0, sizeof winding_stats);
    if (winding_export && result.mesh.nv > 0 && result.mesh.nf > 0) {
        char winding_root[ADCT_PATH_CAP];
        int count = snprintf(
            winding_root, sizeof winding_root,
            "%s/%s_developable_winding_atlas", output_directory, prefix);
        if (count < 0 || (size_t)count >= sizeof winding_root ||
            AtlasTrackExportWinding_run(
                arena, &evidence, &calibration, &result.mesh,
                winding_root, prefix, winding_wraps, winding_slab,
                winding_du, winding_dv, &winding_stats) != 0) {
            fprintf(stderr, "%s: winding organizer export failed\n", argv[0]);
            io = -1;
        } else {
            fprintf(stderr,
                "[atlas_developable_complete] winding atlas=%zu/%zu "
                "valid=%zu conflict=%zu (%.4f%%)\n",
                winding_stats.written, winding_stats.pieces,
                winding_stats.valid_pixels, winding_stats.conflict_pixels,
                100.0 * winding_stats.conflict_fraction);
        }
    }
    double artifacts_done = ves_clock_sec();
    double timing[5] = {
        loaded - start,
        extracted - loaded,
        core_done - extracted,
        completion_done - core_done,
        artifacts_done - completion_done
    };
    io |= adct_write_stats(
        output_directory, prefix, placed_directory, solution_path,
        &pieces, kept_faces, &evidence_options, &grow_options,
        &complete_options, &evidence, &seed, &result, &winding_stats, timing);
    fprintf(stderr,
        "[atlas_developable_complete] timing load=%.3fs evidence=%.3fs "
        "core=%.3fs completion=%.3fs artifacts=%.3fs total=%.3fs\n",
        timing[0], timing[1], timing[2], timing[3], timing[4],
        timing[0] + timing[1] + timing[2] + timing[3] + timing[4]);
    Arena_dispose(&arena);
    return io == 0 ? 0 : 1;
}
