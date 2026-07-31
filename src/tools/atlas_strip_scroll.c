/*
 * atlas_strip_scroll.c
 *
 * Geometry-first adapter from a scroll_whole placed directory to the
 * StrokeStrip-style coarse atlas solve.  This tool deliberately writes a
 * diagnostic checkpoint before solving: candidate admission is based on
 * retained 3-D geometry, never on the already-registered uv/phi field.
 */

#include "../common/arena.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"
#include "../whole/atlas_component_lift.h"
#include "../whole/atlas_boxcut.h"
#include "../whole/atlas_boxcut_layout.h"
#include "../whole/atlas_overlap_audit.h"
#include "../whole/atlas_radial_order.h"
#include "../whole/atlas_register.h"
#include "../whole/atlas_seam_audit.h"
#include "../whole/atlas_place_search.h"
#include "../whole/atlas_sheet_split.h"
#include "../whole/atlas_turn_advance.h"
#include "../whole/atlas_warp.h"
#include "../whole/atlas_xyz_weld_audit.h"
#include "../whole/atlas_winding_sync.h"
#include "../whole/atlas_candidates.h"
#include "../whole/atlas_field_apply.h"
#include "../whole/atlas_field_refine.h"
#include "../whole/atlas_strip.h"
#include "../whole/monotone_qp.h"
#include "../split/multicut_wrap.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AS_PATH_CAP 2048

typedef struct {
    const char *placed_dir;
    const char *out_dir;
    AtlasCandidateOptions candidate;
    AtlasStripOptions strip;
    AtlasStripRobustOptions robust;
    MonotoneQpOptions qp;
    const char *field_from_csv;
    const char *rank_pack_labels;
    const char *boxcut_shifts_from_csv;
    const char *boxcut_export_labels;
    AtlasCandidateRefineOptions refine;
    AtlasRadialOrderOptions order;
    AtlasWarpOptions warp;
    AtlasSheetSplitOptions sheet;
    int isolate_support;   /* -1 = off; else the one support component to work on */
    int sheet_no_shift;    /* 1 = detect and dump, but leave u alone */
    int sheet_no_weld;     /* 1 = fragments stay separate; skip the BPA merge */
    int sheet_uniform_turns; /* non-zero: move the WHOLE block by this many
                              * windings instead of one per sheet */
    int sheet_all;         /* 1 = every support component at once */
    int sheet_dump_piece;  /* >=0: dump this welded piece for inspection */
    double sheet_wind_tol; /* |dw| above which a weld bridge edge is a
                            * cross-wrap short-circuit; 0 disables the cut */
    int sheet_search;      /* 1 = search the wind per piece instead of taking
                            * it from the (unregistered) phi turn graph */
    int sheet_search_span; /* candidate winds are -span..+span (4) */
    double sheet_search_cell;
    double sheet_search_lambda;
    double sheet_search_pair;   /* <0 = module default */
    int sheet_rigid;       /* 1 = one constant shift, not the phi-dependent
                            * shear (the two differ by b*k*dphi across a block) */
    int strip_wind_gate;   /* 1 = reject cross-wrap relations at construction */
    double strip_wind_tol; /* |dwind| admission tolerance, turns */
    int strip_geo_gauge;   /* 1 = orient/gauge strokes from the radius-anchored
                            * winding instead of the raw per-cube u */
    int strip_wind_cut;    /* 1 = segment slice chains at wrap crossings */
    int user_continuation_radius; /* CLI override seen */
    int user_match_radius;        /* CLI override seen */
    int build_only;
    int boxcut_layout_only;
    int boxcut_component_blocks;
    int diag_boxcut_rigid;
    int diag_component_lift;
    int diag_classified_resolve;
    int order_only;
    int prune_fused_links;
    int warp_per_sample;
    int collapse_dump_count;

    /* Run-local diagnostic plumbing.  Copies of ScrollConfig intentionally
     * share the counter so every stage lands in one chronological directory. */
    const char *pipeline_trace_dir;
    int *pipeline_trace_next_frame;
    const float *pipeline_trace_raw_u;
} ScrollConfig;

static int build_boxcut_quality_filtered_piece_set(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *source,
    const ScaffoldCalib *cal, const double *first_u, PieceSet *filtered);
static int sheet_in_target(int32_t target, int32_t support);

static int as_path(char path[AS_PATH_CAP], const char *dir, const char *name)
{
    int n = snprintf(path, AS_PATH_CAP, "%s/%s", dir, name);
    return n < 0 || n >= AS_PATH_CAP ? -1 : 0;
}

static FILE *as_open(const char *dir, const char *name, const char *mode)
{
    char path[AS_PATH_CAP];
    if (as_path(path, dir, name) != 0 || ves_ensure_parent_dir(path) != 0)
        return NULL;
    return fopen(path, mode);
}

static void as_json_number(char out[64], double value)
{
    if (isfinite(value)) snprintf(out, 64, "%.17g", value);
    else snprintf(out, 64, "null");
}

static void as_write_json_string(FILE *fp, const char *value)
{
    fputc('"', fp);
    for (const unsigned char *p = (const unsigned char *)value; *p; p++) {
        switch (*p) {
        case '"': fputs("\\\"", fp); break;
        case '\\': fputs("\\\\", fp); break;
        case '\b': fputs("\\b", fp); break;
        case '\f': fputs("\\f", fp); break;
        case '\n': fputs("\\n", fp); break;
        case '\r': fputs("\\r", fp); break;
        case '\t': fputs("\\t", fp); break;
        default:
            if (*p < 0x20) fprintf(fp, "\\u%04x", (unsigned)*p);
            else fputc(*p, fp);
        }
    }
    fputc('"', fp);
}

static int load_raw_u(Arena_T arena, const char *dir,
                      const PieceSet *ps, float **out_u)
{
    float *u = (float *)ARENA_ALLOC(arena, ps->nv * sizeof(float));
    for (size_t c = 0; c < ps->n_cubes; c++) {
        char path[AS_PATH_CAP];
        size_t nv = ps->cube_voff[c + 1] - ps->cube_voff[c];
        if (snprintf(path, sizeof path, "%s/%s_uvphi_raw.f32",
                     dir, ps->ids[c]) < 0)
            return -1;
        FILE *fp = fopen(path, "rb");
        if (fp == NULL) return -1;
        float triple[3];
        for (size_t i = 0; i < nv; i++) {
            if (fread(triple, sizeof(float), 3, fp) != 3) {
                fclose(fp);
                return -1;
            }
            u[ps->cube_voff[c] + i] = triple[0];
        }
        if (fgetc(fp) != EOF) {
            fclose(fp);
            return -1;
        }
        fclose(fp);
    }
    *out_u = u;
    return 0;
}

static int load_chart_shifts_csv(Arena_T arena, const char *path,
                                 size_t ncharts, double **out_shift)
{
    if (arena == NULL || path == NULL || ncharts == 0 ||
        out_shift == NULL)
        return -1;
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    char line[4096];
    if (fgets(line, sizeof line, fp) == NULL) {
        fclose(fp);
        return -1;
    }
    double *shift = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(double));
    uint8_t *seen = (uint8_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(uint8_t));
    size_t count = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
        char *end_chart = NULL;
        unsigned long long chart = strtoull(line, &end_chart, 10);
        if (end_chart == line || *end_chart != ',' || chart >= ncharts ||
            seen[chart]) {
            fclose(fp);
            return -1;
        }
        char *value_start = end_chart + 1;
        char *end_shift = NULL;
        double value = strtod(value_start, &end_shift);
        if (end_shift == value_start) {
            fclose(fp);
            return -1;
        }
        while (*end_shift == ' ' || *end_shift == '\t' ||
               *end_shift == '\r' || *end_shift == '\n')
            end_shift++;
        if (*end_shift != '\0' || !isfinite(value)) {
            fclose(fp);
            return -1;
        }
        shift[chart] = value;
        seen[chart] = 1;
        count++;
    }
    fclose(fp);
    if (count != ncharts) return -1;
    *out_shift = shift;
    return 0;
}

static int load_boxcut_face_labels(Arena_T arena, const char *path,
                                   size_t nfaces, int32_t **out_label,
                                   size_t *out_charts)
{
    if (arena == NULL || path == NULL || nfaces == 0 || out_label == NULL ||
        out_charts == NULL || nfaces > SIZE_MAX / sizeof(int32_t))
        return -1;
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    int32_t *label = (int32_t *)ARENA_ALLOC(
        arena, nfaces * sizeof(int32_t));
    if (fread(label, sizeof(int32_t), nfaces, fp) != nfaces ||
        fgetc(fp) != EOF) {
        fclose(fp);
        return -1;
    }
    if (fclose(fp) != 0) return -1;
    size_t charts = 0;
    for (size_t f = 0; f < nfaces; f++) {
        if (label[f] < 0 || (size_t)label[f] >= (size_t)INT32_MAX)
            return -1;
        size_t count = (size_t)label[f] + 1;
        if (count > charts) charts = count;
    }
    if (charts == 0) return -1;
    *out_label = label;
    *out_charts = charts;
    return 0;
}

static void component_color(int32_t component, double rgb[3])
{
    unsigned x = (unsigned)(component >= 0 ? component : 0) * 2654435761u;
    rgb[0] = 0.25 + 0.70 * (double)((x >>  0) & 255u) / 255.0;
    rgb[1] = 0.25 + 0.70 * (double)((x >>  8) & 255u) / 255.0;
    rgb[2] = 0.25 + 0.70 * (double)((x >> 16) & 255u) / 255.0;
}

static int write_strokes_obj(const char *dir, const char *name,
                             const AtlasCandidateSet *set,
                             const double *u, int parameter_space)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    for (size_t i = 0; i < set->problem.nsamples; i++) {
        int32_t sid = set->samples[i].stroke;
        double rgb[3];
        component_color(set->strokes[sid].component, rgb);
        if (parameter_space) {
            fprintf(fp, "v %.9g %.9g 0 %.4f %.4f %.4f\n",
                    u[i], set->sample_ref[i].axial,
                    rgb[0], rgb[1], rgb[2]);
        } else {
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                    set->samples[i].p[0], set->samples[i].p[1],
                    set->samples[i].p[2], rgb[0], rgb[1], rgb[2]);
        }
    }
    for (size_t s = 0; s < set->problem.nstrokes; s++) {
        const AtlasStripStroke *st = &set->strokes[s];
        fprintf(fp, "g stroke_%zu_component_%d\n", s, st->component);
        for (int32_t j = 1; j < st->count; j++)
            fprintf(fp, "l %d %d\n", st->first + j,
                    st->first + j + 1);
    }
    fclose(fp);
    return 0;
}

static int write_continuations_obj(const char *dir,
                                   const AtlasCandidateSet *set)
{
    FILE *fp = as_open(dir, "continuations_world.obj", "wb");
    if (fp == NULL) return -1;
    size_t vi = 1;
    for (size_t i = 0; i < set->problem.ncontinuations; i++) {
        const AtlasStripContinuation *c = &set->continuations[i];
        const AtlasStripSample *a = &set->samples[c->a];
        const AtlasStripSample *b = &set->samples[c->b];
        int active = set->continuation_state == NULL ||
                     set->continuation_state[i] != ATLAS_CANDIDATE_INACTIVE;
        fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f 0.1\n",
                a->p[0], a->p[1], a->p[2],
                active ? 0.1 : 1.0, active ? 0.9 : 0.1);
        fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f 0.1\n",
                b->p[0], b->p[1], b->p[2],
                active ? 0.1 : 1.0, active ? 0.9 : 0.1);
        fprintf(fp, "l %zu %zu\n", vi, vi + 1);
        vi += 2;
    }
    fclose(fp);
    return 0;
}

static int write_candidates_obj(const char *dir, const char *name,
                                const AtlasCandidateSet *set,
                                const double *membership)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    size_t vi = 1;
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        if (cs->count < 2) continue;
        const AtlasStripMember *source = &set->members[cs->first];
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            const AtlasStripMember *m = &set->members[mi];
            double q = membership != NULL && m->membership > 0.0
                     ? membership[mi] / m->membership
                     : (m->membership > 0.0 ? 1.0 : 0.0);
            if (q < 0.0) q = 0.0;
            if (q > 1.0) q = 1.0;
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f 0.1\n",
                    source->p[0], source->p[1], source->p[2],
                    1.0 - q, q);
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f 0.1\n",
                    m->p[0], m->p[1], m->p[2], 1.0 - q, q);
            fprintf(fp, "l %zu %zu\n", vi, vi + 1);
            vi += 2;
        }
    }
    fclose(fp);
    return 0;
}

static int write_samples_csv(const char *dir, const PieceSet *ps,
                             const AtlasCandidateSet *set,
                             const double *relaxed, const double *final)
{
    FILE *fp = as_open(dir, "samples.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "sample,stroke,ordinal,component,cube,cube_id,plane,axial,"
                "z,y,x,s,raw_u,initial_u,relaxed_u,final_u,edge0,edge1,edge_t\n");
    for (size_t i = 0; i < set->problem.nsamples; i++) {
        const AtlasStripSample *s = &set->samples[i];
        const AtlasCandidateSampleRef *r = &set->sample_ref[i];
        const char *id = r->cube >= 0 && (size_t)r->cube < ps->n_cubes
                       ? ps->ids[r->cube] : "?";
        fprintf(fp, "%zu,%d,%d,%d,%d,%s,%d,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g\n",
                i, s->stroke, s->ordinal, set->strokes[s->stroke].component,
                r->cube, id, r->plane, r->axial,
                s->p[0], s->p[1], s->p[2], s->s, r->raw_u,
                set->initial_u[i], relaxed != NULL ? relaxed[i] : NAN,
                final != NULL ? final[i] : NAN,
                r->mesh_vertex[0], r->mesh_vertex[1], r->mesh_t);
    }
    fclose(fp);
    return 0;
}

static int write_constraints_csv(const char *dir,
                                 const AtlasCandidateSet *set,
                                 const double *membership,
                                 const double *relaxed, const double *final)
{
    FILE *fp = as_open(dir, "cross_sections.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "cross,member,observation,value0,value1,value_t,prior,"
                "membership,raw_value,relaxed_value,final_value\n");
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        for (int32_t j = 0; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            const AtlasStripMember *m = &set->members[mi];
            fprintf(fp, "%zu,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    c, j, m->observation, m->value0, m->value1, m->value_t,
                    m->membership, membership != NULL ? membership[mi] : NAN,
                    AtlasStrip_member_value(m, set->initial_u),
                    relaxed != NULL ? AtlasStrip_member_value(m, relaxed) : NAN,
                    final != NULL ? AtlasStrip_member_value(m, final) : NAN);
        }
    }
    fclose(fp);
    return 0;
}

static const char *candidate_state_name(uint8_t state)
{
    switch ((AtlasCandidateState)state) {
    case ATLAS_CANDIDATE_SOURCE: return "source";
    case ATLAS_CANDIDATE_BASELINE: return "baseline";
    case ATLAS_CANDIDATE_TOPOLOGY: return "topology";
    case ATLAS_CANDIDATE_REFINED_DELAMINATION: return "delamination";
    case ATLAS_CANDIDATE_REFINED_SHEET_SEAM: return "same_sheet_seam";
    default: return "inactive";
    }
}

static const char *bundle_decision_name(AtlasCandidateBundleDecision decision)
{
    switch (decision) {
    case ATLAS_BUNDLE_SEPARATE_SHEET: return "separate_sheet";
    case ATLAS_BUNDLE_DELAMINATION: return "delamination";
    case ATLAS_BUNDLE_SAME_SHEET_SEAM: return "same_sheet_seam";
    default: return "unclassified";
    }
}

static int write_selection_csv(const char *dir,
                               const AtlasCandidateSet *set,
                               const double *frozen_u)
{
    FILE *fp = as_open(dir, "candidate_selection.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "cross,member,source_component,target_component,prior,"
                "active_membership,state,bundle,bundle_decision,"
                "frozen_residual\n");
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        const AtlasStripMember *source = &set->members[cs->first];
        int32_t source_component =
            set->sample_ref[source->value0].mesh_component;
        double source_u = frozen_u != NULL
                        ? AtlasStrip_member_value(source, frozen_u) : NAN;
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            const AtlasStripMember *m = &set->members[mi];
            int32_t target_component =
                set->sample_ref[m->value0].mesh_component;
            int32_t bi = set->member_bundle != NULL
                       ? set->member_bundle[mi] : -1;
            const char *decision = bi >= 0 && (size_t)bi < set->nbundles
                ? bundle_decision_name(set->bundles[bi].decision) : "none";
            double residual = frozen_u != NULL
                ? AtlasStrip_member_value(m, frozen_u) - source_u : NAN;
            fprintf(fp, "%zu,%d,%d,%d,%.17g,%.17g,%s,%d,%s,%.17g\n",
                    c, j, source_component, target_component,
                    set->member_prior[mi], m->membership,
                    candidate_state_name(set->member_state[mi]), bi,
                    decision, residual);
        }
    }
    fclose(fp);

    fp = as_open(dir, "continuation_selection.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "continuation,component_a,component_b,plane,prior,"
                "active_membership,state,bundle,bundle_decision,"
                "frozen_residual\n");
    for (size_t i = 0; i < set->problem.ncontinuations; i++) {
        const AtlasStripContinuation *cn = &set->continuations[i];
        int32_t ca = set->sample_ref[cn->a].mesh_component;
        int32_t cb = set->sample_ref[cn->b].mesh_component;
        int32_t bi = set->continuation_bundle != NULL
                   ? set->continuation_bundle[i] : -1;
        const char *decision = bi >= 0 && (size_t)bi < set->nbundles
            ? bundle_decision_name(set->bundles[bi].decision) : "none";
        double residual = frozen_u != NULL
            ? frozen_u[cn->a] - frozen_u[cn->b] - cn->target : NAN;
        fprintf(fp, "%zu,%d,%d,%d,%.17g,%.17g,%s,%d,%s,%.17g\n",
                i, ca, cb, set->sample_ref[cn->a].plane,
                set->continuation_prior[i], cn->membership,
                candidate_state_name(set->continuation_state[i]), bi,
                decision, residual);
    }
    fclose(fp);
    return 0;
}

static int write_bundles_csv(const char *dir,
                             const AtlasCandidateSet *set)
{
    FILE *fp = as_open(dir, "component_bundles.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "bundle,component0,component1,links,cross_links,"
                "continuation_links,planes,plane_min,plane_max,"
                "component0_lifespan,component1_lifespan,"
                "component0_extension,component1_extension,"
                "contact_overlap_fraction,initial_residual_median,"
                "initial_residual_mad,frozen_residual_median,"
                "frozen_residual_mad,prior_median,evidence,"
                "top_prior_fraction,topology_alternative_fraction,"
                "arc_correlation,first_plane_gap,"
                "last_plane_gap,component0_runner_up_ratio,"
                "component1_runner_up_ratio,mutual_best,decision\n");
    for (size_t i = 0; i < set->nbundles; i++) {
        const AtlasCandidateBundle *b = &set->bundles[i];
        fprintf(fp, "%zu,%d,%d,%zu,%zu,%zu,%zu,%d,%d,%d,%d,%d,%d,"
                    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%.17g,%d,%s\n",
                i, b->component0, b->component1, b->links,
                b->cross_links, b->continuation_links, b->planes,
                b->plane_min, b->plane_max, b->component0_lifespan,
                b->component1_lifespan, b->component0_extension,
                b->component1_extension, b->contact_overlap_fraction,
                b->initial_residual_median, b->initial_residual_mad,
                b->frozen_residual_median, b->frozen_residual_mad,
                b->prior_median, b->evidence, b->top_prior_fraction,
                b->topology_alternative_fraction, b->arc_correlation,
                b->first_plane_gap,
                b->last_plane_gap, b->component0_runner_up_ratio,
                b->component1_runner_up_ratio, b->mutual_best,
                bundle_decision_name(b->decision));
    }
    fclose(fp);
    return 0;
}

static int write_bundle_decisions_obj(const char *dir,
                                      const AtlasCandidateSet *set)
{
    FILE *fp = as_open(dir, "overlap_classification_world.obj", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "# vertex colors: red=missed sheet (removed), "
        "cyan=delamination (kept), green=same-sheet seam (kept), "
        "gray=inconclusive (kept)\n");
    size_t vi = 1;
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        const AtlasStripMember *source = &set->members[cs->first];
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            int32_t bi = set->member_bundle[mi];
            if (bi < 0 || (size_t)bi >= set->nbundles) continue;
            AtlasCandidateBundleDecision decision = set->bundles[bi].decision;
            double r = 0.45, g = 0.45, b = 0.45;
            if (decision == ATLAS_BUNDLE_SEPARATE_SHEET) {
                r = 1.0; g = 0.05; b = 0.05;
            } else if (decision == ATLAS_BUNDLE_DELAMINATION) {
                r = 0.05; g = 0.9; b = 1.0;
            } else if (decision == ATLAS_BUNDLE_SAME_SHEET_SEAM) {
                r = 0.1; g = 1.0; b = 0.1;
            }
            const AtlasStripMember *target = &set->members[mi];
            fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f %.3f\n",
                    source->p[0], source->p[1], source->p[2], r, g, b);
            fprintf(fp, "v %.9g %.9g %.9g %.3f %.3f %.3f\n",
                    target->p[0], target->p[1], target->p[2], r, g, b);
            fprintf(fp, "l %zu %zu\n", vi, vi + 1);
            vi += 2;
        }
    }
    fclose(fp);
    return 0;
}

static int write_robust_trace(const char *dir,
                              const AtlasStripRobustTrace *trace)
{
    FILE *fp = as_open(dir, "robust_solver_trace.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "iteration,phase,qp_rc,objective,max_u_change,"
                "max_membership_change,residual_rms,max_residual,"
                "membership_min,membership_mean,membership_max,downweighted\n");
    for (size_t i = 0; i < trace->count; i++) {
        const AtlasStripRobustTraceEntry *e = &trace->entry[i];
        fprintf(fp, "%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%zu\n",
                e->iteration, e->phase, e->qp_rc, e->objective,
                e->max_u_change, e->max_membership_change,
                e->residual_rms, e->max_residual, e->membership_min,
                e->membership_mean, e->membership_max,
                e->downweighted_members);
    }
    fclose(fp);
    return 0;
}

static int write_summary(const ScrollConfig *cfg, const PieceSet *ps,
                         const ScaffoldCalib *cal,
                         const AtlasCandidateSet *set, int robust_rc,
                         int final_rc, const AtlasStripRobustStats *rs,
                         const MonotoneQpStats *fs,
                         const AtlasStripMetrics *rm,
                         const AtlasStripMetrics *fm,
                         int fixed_head_membership)
{
    FILE *fp = as_open(cfg->out_dir, "summary.json", "wb");
    if (fp == NULL) return -1;
    const AtlasCandidateStats *s = &set->stats;
    char rsigma[64], fobjective[64], rmetric[6][64], fmetric[6][64];
    as_json_number(rsigma, rs ? rs->likelihood_sigma : NAN);
    as_json_number(fobjective, fs ? fs->objective_final : NAN);
    as_json_number(rmetric[0], rm ? rm->rms_length_row : NAN);
    as_json_number(rmetric[1], rm ? rm->rms_align_row : NAN);
    as_json_number(rmetric[2], rm ? rm->rms_continuation : NAN);
    as_json_number(rmetric[3], rm ? rm->max_align_residual : NAN);
    as_json_number(rmetric[4], rm ? rm->max_continuation_residual : NAN);
    as_json_number(rmetric[5], rm ? rm->min_monotone_ratio : NAN);
    as_json_number(fmetric[0], fm ? fm->rms_length_row : NAN);
    as_json_number(fmetric[1], fm ? fm->rms_align_row : NAN);
    as_json_number(fmetric[2], fm ? fm->rms_continuation : NAN);
    as_json_number(fmetric[3], fm ? fm->max_align_residual : NAN);
    as_json_number(fmetric[4], fm ? fm->max_continuation_residual : NAN);
    as_json_number(fmetric[5], fm ? fm->min_monotone_ratio : NAN);
    fputs("{\n  \"placed_dir\": ", fp);
    as_write_json_string(fp, cfg->placed_dir);
    fprintf(fp,
        ",\n  \"build_only\": %s,\n"
        "  \"mesh\": {\"cubes\": %zu, \"vertices\": %zu, \"faces\": %zu},\n"
        "  \"calibration\": {\"pitch\": %.17g, \"spiral_a\": %.17g, "
        "\"spiral_b\": %.17g, \"sense\": %d},\n"
        "  \"candidate_options\": {\"slice_spacing\": %.17g, "
        "\"min_stroke_length\": %.17g, \"sample_spacing\": %.17g, "
        "\"match_radius\": %.17g, "
        "\"match_angle_deg\": %.17g, \"max_slice_stride\": %d, "
        "\"max_candidates\": %d, \"continuation_radius\": %.17g, "
        "\"continuation_angle_deg\": %.17g},\n"
        "  \"candidates\": {\"intersection_nodes\": %zu, "
        "\"intersection_segments\": %zu, \"branch_nodes\": %zu, "
        "\"short_strokes_dropped\": %zu, \"strokes\": %zu, "
        "\"samples\": %zu, \"continuations\": %zu, "
        "\"cross_sections\": %zu, \"members\": %zu, "
        "\"matched_source_samples\": %zu, \"ambiguous_cross_sections\": %zu, "
        "\"truncated_ball_queries\": %zu, \"support_components\": %zu, "
        "\"largest_component_samples\": %zu, \"source_match_fraction\": %.17g},\n"
        "  \"selection\": {\"mesh_components\": %zu, "
        "\"candidate_links\": %zu, \"topology_links\": %zu, "
        "\"topology_cross_sections\": %zu, \"topology_ambiguous\": %zu, "
        "\"topology_continuations\": %zu, \"baseline_unmatched\": %zu, "
        "\"bundles\": %zu, \"mutual_bundles\": %zu, "
        "\"separate_sheet_bundles\": %zu, \"delamination_bundles\": %zu, "
        "\"same_sheet_seam_bundles\": %zu, \"inconclusive_bundles\": %zu, "
        "\"pruned_candidate_links\": %zu, "
        "\"pruned_continuations\": %zu, \"remaining_active_cross_sections\": %zu},\n"
        "  \"solve\": {\"attempted\": %s, \"membership_mode\": \"%s\", "
        "\"robust_rc\": %d, \"final_rc\": %d, "
        "\"robust_qp_solves\": %d, \"robust_sigma\": %s, "
        "\"final_iterations\": %d, \"final_objective\": %s},\n"
        "  \"relaxed_metrics\": {\"rms_length\": %s, "
        "\"rms_align\": %s, \"rms_continuation\": %s, "
        "\"max_align\": %s, \"max_continuation\": %s, "
        "\"min_monotone_ratio\": %s},\n"
        "  \"final_metrics\": {\"rms_length\": %s, "
        "\"rms_align\": %s, \"rms_continuation\": %s, "
        "\"max_align\": %s, \"max_continuation\": %s, "
        "\"min_monotone_ratio\": %s}\n}\n",
        cfg->build_only ? "true" : "false",
        ps->n_cubes, ps->nv, ps->nf,
        cal->pitch, cal->spiral_a, cal->spiral_b, cal->sense,
        cfg->candidate.slice_spacing, cfg->candidate.min_stroke_length,
        cfg->candidate.sample_spacing, cfg->candidate.match_radius,
        cfg->candidate.match_angle_deg,
        cfg->candidate.max_slice_stride, cfg->candidate.max_candidates,
        cfg->candidate.continuation_radius,
        cfg->candidate.continuation_angle_deg,
        s->intersection_nodes, s->intersection_segments, s->branch_nodes,
        s->short_strokes_dropped, s->strokes, s->samples,
        s->continuations, s->cross_sections, s->members,
        s->matched_source_samples, s->ambiguous_cross_sections,
        s->truncated_ball_queries, s->support_components,
        s->largest_component_samples, s->source_match_fraction,
        set->selection.mesh_components, set->selection.candidate_links,
        set->selection.topology_candidate_links,
        set->selection.topology_selected_cross_sections,
        set->selection.topology_ambiguous_cross_sections,
        set->selection.topology_selected_continuations,
        set->selection.unmatched_cross_sections,
        set->selection.refinement_bundles,
        set->selection.refinement_mutual_bundles,
        set->selection.refinement_separate_sheet_bundles,
        set->selection.refinement_delamination_bundles,
        set->selection.refinement_sheet_seam_bundles,
        set->selection.refinement_inconclusive_bundles,
        set->selection.refinement_pruned_candidate_links,
        set->selection.refinement_pruned_continuations,
        set->selection.refinement_remaining_active_cross_sections,
        cfg->build_only ? "false" : "true",
        fixed_head_membership ? "frozen_head" : "robust_irls",
        robust_rc, final_rc,
        rs ? rs->total_qp_solves : 0, rsigma,
        fs ? fs->iterations : 0, fobjective,
        rmetric[0], rmetric[1], rmetric[2], rmetric[3], rmetric[4], rmetric[5],
        fmetric[0], fmetric[1], fmetric[2], fmetric[3], fmetric[4], fmetric[5]);
    fclose(fp);
    return 0;
}

static int load_sample_solution_csv(const char *path, size_t nsample,
                                    double *target)
{
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    char line[4096];
    if (fgets(line, sizeof line, fp) == NULL) {
        fclose(fp);
        return -1;
    }
    size_t n = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
        char *end_index = NULL;
        unsigned long long index = strtoull(line, &end_index, 10);
        if (end_index == line || *end_index != ',' || index != n ||
            n >= nsample) {
            fclose(fp);
            return -1;
        }
        char *p = line;
        for (int column = 0; column < 15; column++) {
            p = strchr(p, ',');
            if (p == NULL) {
                fclose(fp);
                return -1;
            }
            p++;
        }
        char *end_value = NULL;
        double value = strtod(p, &end_value);
        if (end_value == p || !isfinite(value)) {
            fclose(fp);
            return -1;
        }
        target[n++] = value;
    }
    fclose(fp);
    return n == nsample ? 0 : -1;
}

static int companion_csv_path(char out[AS_PATH_CAP], const char *path,
                              const char *name)
{
    if (out == NULL || path == NULL || name == NULL) return -1;
    size_t length = strlen(path);
    if (length >= AS_PATH_CAP) return -1;
    memcpy(out, path, length + 1);
    char *slash0 = strrchr(out, '/');
    char *slash1 = strrchr(out, '\\');
    char *slash = slash0 != NULL && (slash1 == NULL || slash0 > slash1)
                ? slash0 : slash1;
    size_t prefix = slash != NULL ? (size_t)(slash - out + 1) : 0;
    size_t name_length = strlen(name);
    if (prefix + name_length >= AS_PATH_CAP) return -1;
    memcpy(out + prefix, name, name_length + 1);
    return 0;
}

static int load_member_membership_csv(const char *path,
                                      const AtlasCandidateSet *set,
                                      double *target)
{
    if (path == NULL || set == NULL || target == NULL) return -1;
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    char line[4096];
    if (fgets(line, sizeof line, fp) == NULL) {
        fclose(fp);
        return -1;
    }
    size_t n = 0;
    while (fgets(line, sizeof line, fp) != NULL) {
        char *end_cross = NULL;
        unsigned long long cross = strtoull(line, &end_cross, 10);
        if (end_cross == line || *end_cross != ',' ||
            cross >= set->problem.ncross_sections) {
            fclose(fp);
            return -1;
        }
        char *end_member = NULL;
        unsigned long member = strtoul(end_cross + 1, &end_member, 10);
        const AtlasStripCrossSection *cs = &set->cross_sections[cross];
        if (end_member == end_cross + 1 || *end_member != ',' ||
            member >= (unsigned long)cs->count ||
            cs->first + (size_t)member != n || n >= set->problem.nmembers) {
            fclose(fp);
            return -1;
        }
        const AtlasStripMember *expected = &set->members[n];
        char *field = end_member + 1;
        char *end_field = NULL;
        long observation = strtol(field, &end_field, 10);
        if (end_field == field || *end_field != ',' ||
            observation != expected->observation) {
            fclose(fp);
            return -1;
        }
        field = end_field + 1;
        long value0 = strtol(field, &end_field, 10);
        if (end_field == field || *end_field != ',' ||
            value0 != expected->value0) {
            fclose(fp);
            return -1;
        }
        field = end_field + 1;
        long value1 = strtol(field, &end_field, 10);
        if (end_field == field || *end_field != ',' ||
            value1 != expected->value1) {
            fclose(fp);
            return -1;
        }
        field = end_field + 1;
        double value_t = strtod(field, &end_field);
        if (end_field == field || *end_field != ',' ||
            !isfinite(value_t) ||
            fabs(value_t - expected->value_t) >
                1e-12 * (1.0 + fabs(expected->value_t))) {
            fclose(fp);
            return -1;
        }
        field = end_field + 1;
        double prior = strtod(field, &end_field);
        if (end_field == field || *end_field != ',' || !isfinite(prior) ||
            fabs(prior - set->member_prior[n]) >
                1e-12 * (1.0 + fabs(set->member_prior[n]))) {
            fclose(fp);
            return -1;
        }
        char *p = line;
        for (int column = 0; column < 7; column++) {
            p = strchr(p, ',');
            if (p == NULL) {
                fclose(fp);
                return -1;
            }
            p++;
        }
        char *end_value = NULL;
        double value = strtod(p, &end_value);
        if (end_value == p || !isfinite(value) || value < 0.0) {
            fclose(fp);
            return -1;
        }
        target[n++] = value;
    }
    fclose(fp);
    return n == set->problem.nmembers ? 0 : -1;
}

static int copy_file(const char *source, const char *destination)
{
    FILE *in = fopen(source, "rb");
    if (in == NULL) return -1;
    if (ves_ensure_parent_dir(destination) != 0) {
        fclose(in);
        return -1;
    }
    FILE *out = fopen(destination, "wb");
    if (out == NULL) {
        fclose(in);
        return -1;
    }
    unsigned char buffer[65536];
    int rc = 0;
    for (;;) {
        size_t n = fread(buffer, 1, sizeof buffer, in);
        if (n > 0 && fwrite(buffer, 1, n, out) != n) {
            rc = -1;
            break;
        }
        if (n < sizeof buffer) {
            if (ferror(in)) rc = -1;
            break;
        }
    }
    if (fclose(out) != 0) rc = -1;
    fclose(in);
    return rc;
}

static double axis_coordinate(const ScaffoldCalib *cal, const float *p)
{
    double axis[3] = {(double)cal->axis_dir[0],
                      (double)cal->axis_dir[1],
                      (double)cal->axis_dir[2]};
    double norm = sqrt(axis[0] * axis[0] + axis[1] * axis[1] +
                       axis[2] * axis[2]);
    if (norm <= 1e-15) return 0.0;
    return ((double)p[0] - (double)cal->axis_point[0]) * axis[0] / norm +
           ((double)p[1] - (double)cal->axis_point[1]) * axis[1] / norm +
           ((double)p[2] - (double)cal->axis_point[2]) * axis[2] / norm;
}

/* The pinned map itself: u(phi) = a*phi + b*phi^2/(4pi).  Its derivative is
 * du/dphi = a + b*phi/(2pi) = r(phi), so u IS arclength along the spiral and
 * u(phi_b) - u(phi_a) is exactly the length of papyrus between two angles. */
static double spiral_u(const ScaffoldCalib *cal, double phi)
{
    const double pi = 3.14159265358979323846;
    return cal->spiral_a * phi + cal->spiral_b * phi * phi / (4.0 * pi);
}

static int export_field_placed(const ScrollConfig *cfg, const PieceSet *ps,
                               const ScaffoldCalib *cal,
                               const double *field_u);
static int write_global_mesh_xyz_obj(const char *dir, const char *name,
                                     const PieceSet *ps,
                                     const ScaffoldCalib *cal,
                                     const double *u);

static int write_global_mesh_obj(const char *dir, const char *name,
                                 const PieceSet *ps,
                                 const ScaffoldCalib *cal,
                                 const double *u, int flat)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    for (size_t i = 0; i < ps->nv; i++) {
        double v = axis_coordinate(cal, &ps->verts[i * 3]);
        if (flat)
            fprintf(fp, "v %.9g %.9g 0\n", u[i], v);
        else
            fprintf(fp, "v %.9g %.9g %.9g\n",
                    (double)ps->verts[i * 3 + 0],
                    (double)ps->verts[i * 3 + 1],
                    (double)ps->verts[i * 3 + 2]);
        fprintf(fp, "vt %.9g %.9g\n", u[i], v);
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int a = ps->faces[f * 3 + 0] + 1;
        int b = ps->faces[f * 3 + 1] + 1;
        int c = ps->faces[f * 3 + 2] + 1;
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(fp);
    return 0;
}

static int write_pipeline_trace_mesh(
    const ScrollConfig *cfg, const PieceSet *ps, const ScaffoldCalib *cal,
    const char *stage, const double *u,
    const AtlasRegisterIterationStats *registration_iteration,
    const AtlasFieldApplyStats *field_stats);
static int write_pipeline_trace_strokes(
    const ScrollConfig *cfg, const AtlasCandidateSet *set,
    const char *stage, const double *u,
    const AtlasStripRobustTraceEntry *robust_iteration);

static int smooth_seam_bundle_eligible(const AtlasSeamBundle *b)
{
    return b != NULL && b->component0 != b->component1 && b->pairs >= 2 &&
           b->target_shift_mad <= 2.0 &&
           fabs(b->phase_residual_median) <= 0.10 &&
           b->phase_residual_mad <= 0.05;
}

static int run_smooth_seam_stage(Arena_T arena, const ScrollConfig *cfg,
                                 const PieceSet *ps,
                                 const ScaffoldCalib *cal,
                                 const AtlasCandidateSet *set,
                                 const double *head_field,
                                 double **out_smooth_field)
{
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "smooth_seams") != 0) return -1;
    if (write_pipeline_trace_mesh(cfg, ps, cal, "smooth_seam_input",
                                  head_field, NULL, NULL) != 0)
        return -1;
    AtlasSeamAudit seam;
    if (AtlasSeamAudit_build(
            arena, ps, head_field, set->vertex_mesh_component,
            set->mesh_components, &seam) != 0)
        return -1;

    uint8_t *eligible_bundle = (uint8_t *)ARENA_CALLOC(
        arena, (seam.nbundles ? seam.nbundles : 1), sizeof(uint8_t));
    size_t neligible_bundle = 0;
    for (size_t i = 0; i < seam.nbundles; i++) {
        eligible_bundle[i] =
            (uint8_t)smooth_seam_bundle_eligible(&seam.bundles[i]);
        if (eligible_bundle[i]) neligible_bundle++;
    }
    size_t neligible_pair = 0;
    for (size_t i = 0; i < seam.npairs; i++) {
        int32_t b = seam.pairs[i].bundle;
        if (b < 0 || (size_t)b >= seam.nbundles) return -1;
        if (eligible_bundle[b]) neligible_pair++;
    }
    if (neligible_pair == 0) return -1;

    AtlasFieldConstraint *constraint =
        (AtlasFieldConstraint *)ARENA_ALLOC(
            arena, neligible_pair * sizeof(AtlasFieldConstraint));
    AtlasFieldConstraintCoeff *coeff =
        (AtlasFieldConstraintCoeff *)ARENA_ALLOC(
            arena, 2 * neligible_pair * sizeof(AtlasFieldConstraintCoeff));
    size_t nr = 0;
    const double seam_weight = 1000.0;
    for (size_t i = 0; i < seam.npairs; i++) {
        const AtlasSeamPair *p = &seam.pairs[i];
        if (!eligible_bundle[p->bundle]) continue;
        AtlasFieldConstraint *row = &constraint[nr];
        row->first = 2 * nr;
        row->count = 2;
        row->target = 0.0;
        row->weight = seam_weight;
        row->kind = ATLAS_FIELD_ROW_SEAM;
        row->source = i <= (size_t)INT32_MAX ? (int32_t)i : INT32_MAX;
        coeff[2 * nr].variable = p->vertex0;
        coeff[2 * nr].coefficient = -1.0;
        coeff[2 * nr + 1].variable = p->vertex1;
        coeff[2 * nr + 1].coefficient = 1.0;
        nr++;
    }
    if (nr != neligible_pair) return -1;

    double *smooth = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    AtlasFieldRefineStats stats;
    double t0 = ves_clock_sec();
    int solve_rc = AtlasFieldRefine_solve(
        arena, ps->verts, ps->nv, ps->faces, ps->nf,
        set->vertex_mesh_component, set->mesh_components, head_field,
        constraint, neligible_pair, coeff, 2 * neligible_pair,
        0, NULL, &cfg->qp, smooth, NULL, &stats);
    fprintf(stderr,
        "[atlas_strip_scroll] smooth seams %.2fs: rc=%d "
        "bundles=%zu pairs=%zu quotient=%zu seam_rms=%.4g->%.4g "
        "correction_rms=%.4g grad_rms=%.4g\n",
        ves_clock_sec() - t0, solve_rc, neligible_bundle, neligible_pair,
        stats.quotient_components, stats.seam_rms_before,
        stats.seam_rms_after, stats.correction_rms,
        stats.edge_gradient_delta_rms);
    if (solve_rc != 0) return -1;
    if (write_pipeline_trace_mesh(cfg, ps, cal, "smooth_seam_output",
                                  smooth, NULL, NULL) != 0)
        return -1;

    double *axial = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    float *registered = (float *)ARENA_ALLOC(
        arena, ps->nv * sizeof(float));
    for (size_t i = 0; i < ps->nv; i++) {
        axial[i] = axis_coordinate(cal, &ps->verts[i * 3]);
        registered[i] = ps->uv[i * 2];
    }
    AtlasOverlapAudit overlap_before, overlap_after;
    if (AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, head_field, axial,
            registered, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &overlap_before) != 0 ||
        AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, smooth, axial,
            registered, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &overlap_after) != 0)
        return -1;

    size_t comparable = 0, relative_flips = 0, degenerate = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = ps->faces[f * 3];
        int32_t b = ps->faces[f * 3 + 1];
        int32_t c = ps->faces[f * 3 + 2];
        double head_det = (head_field[b] - head_field[a]) *
                          (axial[c] - axial[a]) -
                          (head_field[c] - head_field[a]) *
                          (axial[b] - axial[a]);
        double new_det = (smooth[b] - smooth[a]) *
                         (axial[c] - axial[a]) -
                         (smooth[c] - smooth[a]) *
                         (axial[b] - axial[a]);
        if (fabs(new_det) < 1e-10) degenerate++;
        if (fabs(head_det) >= 1e-10) {
            comparable++;
            if (head_det * new_det < 0.0) relative_flips++;
        }
    }

    FILE *fp = as_open(dir, "seam_pairs.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "pair,bundle,eligible,vertex0,vertex1,component0,component1,"
                "distance,phase_residual,head_du,final_du\n");
    for (size_t i = 0; i < seam.npairs; i++) {
        const AtlasSeamPair *p = &seam.pairs[i];
        fprintf(fp, "%zu,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g\n",
                i, p->bundle, eligible_bundle[p->bundle],
                p->vertex0, p->vertex1, p->component0, p->component1,
                p->distance, p->phase_residual,
                head_field[p->vertex1] - head_field[p->vertex0],
                smooth[p->vertex1] - smooth[p->vertex0]);
    }
    fclose(fp);

    fp = as_open(dir, "seam_bundles.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "bundle,eligible,component0,component1,pairs,"
                "target_shift_median,target_shift_mad,distance_median,"
                "phase_residual_median,phase_residual_mad,"
                "registered_du_median,registered_du_mad\n");
    for (size_t i = 0; i < seam.nbundles; i++) {
        const AtlasSeamBundle *b = &seam.bundles[i];
        fprintf(fp, "%zu,%d,%d,%d,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g\n",
                i, eligible_bundle[i], b->component0, b->component1,
                b->pairs, b->target_shift_median, b->target_shift_mad,
                b->distance_median, b->phase_residual_median,
                b->phase_residual_mad, b->registered_du_median,
                b->registered_du_mad);
    }
    fclose(fp);

    ScrollConfig smooth_cfg = *cfg;
    smooth_cfg.out_dir = dir;
    if (write_global_mesh_obj(dir, "atlas_head_flat.obj", ps, cal,
                              head_field, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_smooth_seam_flat.obj", ps, cal,
                              smooth, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_smooth_seam_world_uv.obj", ps, cal,
                              smooth, 0) != 0 ||
        write_global_mesh_xyz_obj(dir, "atlas_head_flat_xyzcolor.obj", ps,
                                  cal, head_field) != 0 ||
        write_global_mesh_xyz_obj(dir, "atlas_smooth_seam_flat_xyzcolor.obj",
                                  ps, cal, smooth) != 0 ||
        export_field_placed(&smooth_cfg, ps, cal, smooth) != 0)
        return -1;

    fp = as_open(dir, "smooth_seam_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n  \"operation\": \"smooth_per_vertex_seam_quotient\",\n"
        "  \"base\": \"immutable_HEAD\",\n"
        "  \"faces_before\": %zu,\n  \"faces_after\": %zu,\n"
        "  \"seam_weight\": %.17g,\n"
        "  \"boundary_vertices\": %zu,\n"
        "  \"mutual_pairs\": %zu,\n"
        "  \"eligible_bundles\": %zu,\n"
        "  \"eligible_pairs\": %zu,\n"
        "  \"quotient_components\": %zu,\n"
        "  \"anchors\": %zu,\n"
        "  \"seam_rms_before\": %.17g,\n"
        "  \"seam_rms_after\": %.17g,\n"
        "  \"seam_max_after\": %.17g,\n"
        "  \"correction_rms\": %.17g,\n"
        "  \"correction_max\": %.17g,\n"
        "  \"edge_gradient_delta_rms\": %.17g,\n"
        "  \"edge_gradient_delta_max\": %.17g,\n"
        "  \"relative_flips\": %zu,\n"
        "  \"comparable_faces\": %zu,\n"
        "  \"degenerate_faces\": %zu,\n"
        "  \"exact_overlap_pairs_before\": %zu,\n"
        "  \"exact_overlap_pairs_after\": %zu,\n"
        "  \"cross_component_overlap_pairs_before\": %zu,\n"
        "  \"cross_component_overlap_pairs_after\": %zu,\n"
        "  \"qp_objective_initial\": %.17g,\n"
        "  \"qp_objective_final\": %.17g,\n"
        "  \"qp_linear_residual\": %.17g\n}\n",
        ps->nf, ps->nf, seam_weight, seam.boundary_vertices,
        seam.mutual_pairs, neligible_bundle, neligible_pair,
        stats.quotient_components, stats.anchors,
        stats.seam_rms_before, stats.seam_rms_after,
        stats.seam_max_after, stats.correction_rms, stats.correction_max,
        stats.edge_gradient_delta_rms, stats.edge_gradient_delta_max,
        relative_flips, comparable, degenerate,
        overlap_before.exact_face_pairs, overlap_after.exact_face_pairs,
        overlap_before.cross_component_pairs,
        overlap_after.cross_component_pairs,
        stats.qp.objective_initial, stats.qp.objective_final,
        stats.qp.max_reduced_linear_residual);
    fclose(fp);
    if (out_smooth_field != NULL) *out_smooth_field = smooth;
    return 0;
}

static int write_global_mesh_xyz_obj(const char *dir, const char *name,
                                     const PieceSet *ps,
                                     const ScaffoldCalib *cal,
                                     const double *u)
{
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    for (size_t i = 0; i < ps->nv; i++) {
        for (int d = 0; d < 3; d++) {
            double x = (double)ps->verts[i * 3 + (size_t)d];
            if (x < lo[d]) lo[d] = x;
            if (x > hi[d]) hi[d] = x;
        }
    }
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    for (size_t i = 0; i < ps->nv; i++) {
        double rgb[3];
        for (int d = 0; d < 3; d++) {
            double span = hi[d] - lo[d];
            rgb[d] = span > 1e-12
                   ? ((double)ps->verts[i * 3 + (size_t)d] - lo[d]) / span
                   : 0.5;
        }
        fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                u[i], axis_coordinate(cal, &ps->verts[i * 3]),
                rgb[0], rgb[1], rgb[2]);
    }
    for (size_t f = 0; f < ps->nf; f++)
        fprintf(fp, "f %d %d %d\n",
                ps->faces[f * 3 + 0] + 1,
                ps->faces[f * 3 + 1] + 1,
                ps->faces[f * 3 + 2] + 1);
    return fclose(fp) == 0 ? 0 : -1;
}

typedef struct {
    int32_t bundle;
    int32_t turn;
    double phase_residual;
    double registered_du;
    double head_du;
} AsLiftObservation;

typedef struct {
    size_t observations;
    int32_t turn_mode;
    double turn_agreement;
    double phase_residual_median;
    double phase_residual_mad;
    double registered_du_median;
    double registered_du_mad;
    double head_du_median;
    double head_du_mad;
    int eligible;
} AsLiftBundleAudit;

static int compare_lift_observation(const void *pa, const void *pb)
{
    const AsLiftObservation *a = (const AsLiftObservation *)pa;
    const AsLiftObservation *b = (const AsLiftObservation *)pb;
    if (a->bundle != b->bundle) return a->bundle < b->bundle ? -1 : 1;
    if (a->turn != b->turn) return a->turn < b->turn ? -1 : 1;
    return 0;
}

static int compare_lift_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double lift_median(double *value, size_t count)
{
    if (count == 0) return NAN;
    qsort(value, count, sizeof(double), compare_lift_double);
    if (count & 1u) return value[count / 2];
    return 0.5 * (value[count / 2 - 1] + value[count / 2]);
}

typedef struct {
    int32_t inner_label;
    int32_t outer_label;
    double correction_target;
    double desired_du;
    double base_du;
    int local_order_consistent;
} AsBoxcutOrderObservation;

typedef struct {
    int32_t inner_label;
    int32_t outer_label;
    size_t pairs;
    double correction_target_median;
    double correction_target_mad;
    double desired_du_median;
    double base_du_median;
    double local_radial_agreement;
} AsBoxcutOrderBundle;

static int compare_boxcut_order_observation(const void *pa, const void *pb)
{
    const AsBoxcutOrderObservation *a =
        (const AsBoxcutOrderObservation *)pa;
    const AsBoxcutOrderObservation *b =
        (const AsBoxcutOrderObservation *)pb;
    if (a->inner_label != b->inner_label)
        return a->inner_label < b->inner_label ? -1 : 1;
    if (a->outer_label != b->outer_label)
        return a->outer_label < b->outer_label ? -1 : 1;
    return a->correction_target < b->correction_target ? -1 :
           (a->correction_target > b->correction_target ? 1 : 0);
}

static double boxcut_face_mean(const PieceSet *ps, size_t face,
                               const double *value)
{
    return (value[ps->faces[face * 3]] +
            value[ps->faces[face * 3 + 1]] +
            value[ps->faces[face * 3 + 2]]) / 3.0;
}

static double boxcut_face_registered_u(const PieceSet *ps, size_t face)
{
    return ((double)ps->uv[(size_t)ps->faces[face * 3] * 2] +
            (double)ps->uv[(size_t)ps->faces[face * 3 + 1] * 2] +
            (double)ps->uv[(size_t)ps->faces[face * 3 + 2] * 2]) / 3.0;
}

static double boxcut_face_radius(const PieceSet *ps,
                                 const ScaffoldCalib *cal, size_t face)
{
    double axis[3] = {(double)cal->axis_dir[0],
                      (double)cal->axis_dir[1],
                      (double)cal->axis_dir[2]};
    double norm2 = axis[0] * axis[0] + axis[1] * axis[1] +
                   axis[2] * axis[2];
    double radius = 0.0;
    for (int k = 0; k < 3; k++) {
        int32_t vertex = ps->faces[face * 3 + (size_t)k];
        const float *p = &ps->verts[(size_t)vertex * 3];
        double d[3] = {(double)p[0] - (double)cal->axis_point[0],
                       (double)p[1] - (double)cal->axis_point[1],
                       (double)p[2] - (double)cal->axis_point[2]};
        double t = (d[0] * axis[0] + d[1] * axis[1] +
                    d[2] * axis[2]) / norm2;
        double q0 = d[0] - t * axis[0];
        double q1 = d[1] - t * axis[1];
        double q2 = d[2] - t * axis[2];
        radius += sqrt(q0 * q0 + q1 * q1 + q2 * q2);
    }
    return radius / 3.0;
}

static uint64_t boxcut_face_pair_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static int boxcut_compare_u64(const void *pa, const void *pb)
{
    uint64_t a = *(const uint64_t *)pa;
    uint64_t b = *(const uint64_t *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int boxcut_has_u64(const uint64_t *value, size_t count, uint64_t key)
{
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (value[mid] < key) lo = mid + 1;
        else hi = mid;
    }
    return lo < count && value[lo] == key;
}

static int boxcut_overlap_pair_is_delamination(
    const ScrollConfig *cfg, const PieceSet *ps,
    const AtlasCandidateSet *set, const double *axial,
    const AtlasOverlapPair *pair)
{
    if (cfg == NULL || ps == NULL || set == NULL || axial == NULL ||
        pair == NULL || pair->face0 < 0 || pair->face1 < 0 ||
        (size_t)pair->face0 >= ps->nf || (size_t)pair->face1 >= ps->nf ||
        !isfinite(cfg->candidate.slice_spacing) ||
        !(cfg->candidate.slice_spacing > 0.0))
        return 0;
    int32_t component0 = pair->component0 < pair->component1
                       ? pair->component0 : pair->component1;
    int32_t component1 = pair->component0 < pair->component1
                       ? pair->component1 : pair->component0;
    double minimum0 = DBL_MAX, maximum0 = -DBL_MAX;
    double minimum1 = DBL_MAX, maximum1 = -DBL_MAX;
    for (int k = 0; k < 3; k++) {
        double a = axial[ps->faces[(size_t)pair->face0 * 3 + (size_t)k]];
        double b = axial[ps->faces[(size_t)pair->face1 * 3 + (size_t)k]];
        if (a < minimum0) minimum0 = a;
        if (a > maximum0) maximum0 = a;
        if (b < minimum1) minimum1 = b;
        if (b > maximum1) maximum1 = b;
    }
    double contact_min = fmax(minimum0, minimum1);
    double contact_max = fmin(maximum0, maximum1);
    if (contact_min > contact_max) return 0;
    double contact_midpoint = 0.5 * (contact_min + contact_max);
    double spacing = cfg->candidate.slice_spacing;
    for (size_t i = 0; i < set->nbundles; i++) {
        const AtlasCandidateBundle *bundle = &set->bundles[i];
        if (bundle->decision != ATLAS_BUNDLE_DELAMINATION) continue;
        int32_t b0 = bundle->component0 < bundle->component1
                   ? bundle->component0 : bundle->component1;
        int32_t b1 = bundle->component0 < bundle->component1
                   ? bundle->component1 : bundle->component0;
        if (b0 != component0 || b1 != component1) continue;
        /* Candidate planes are centred at (k+1/2)h. The enclosing cell
         * boundaries give a conservative axial support for the affirmative
         * bounded-front classification. */
        double axial_min = (double)bundle->plane_min * spacing;
        double axial_max = ((double)bundle->plane_max + 1.0) * spacing;
        if (contact_midpoint >= axial_min && contact_midpoint <= axial_max)
            return 1;
    }
    return 0;
}

static int build_boxcut_component_face_labels(
    Arena_T arena, const PieceSet *ps, const AtlasCandidateSet *set,
    int32_t **out_face_label, size_t *out_components)
{
    if (arena == NULL || ps == NULL || set == NULL || out_face_label == NULL ||
        out_components == NULL || set->vertex_mesh_component == NULL ||
        set->mesh_components == 0 || set->mesh_components > (size_t)INT32_MAX)
        return -1;
    uint8_t *active = (uint8_t *)ARENA_CALLOC(
        arena, set->mesh_components, sizeof(uint8_t));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t component = -1;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            int32_t current = set->vertex_mesh_component[vertex];
            if (current < 0 || (size_t)current >= set->mesh_components ||
                (k > 0 && current != component))
                return -1;
            component = current;
        }
        active[component] = 1;
    }
    int32_t *dense = (int32_t *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(int32_t));
    size_t count = 0;
    for (size_t component = 0; component < set->mesh_components; component++) {
        dense[component] = active[component] ? (int32_t)count++ : -1;
    }
    if (count == 0 || count > (size_t)INT32_MAX) return -1;
    int32_t *label = (int32_t *)ARENA_ALLOC(
        arena, ps->nf * sizeof(int32_t));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t vertex = ps->faces[f * 3];
        int32_t component = set->vertex_mesh_component[vertex];
        label[f] = dense[component];
        if (label[f] < 0) return -1;
    }
    *out_face_label = label;
    *out_components = count;
    return 0;
}

static int write_boxcut_chart_obj(const char *dir, const char *name,
                                  const PieceSet *ps,
                                  const ScaffoldCalib *cal,
                                  const double *u,
                                  const int32_t *face_label)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    size_t vertex = 1;
    for (size_t f = 0; f < ps->nf; f++) {
        double rgb[3];
        component_color(face_label[f], rgb);
        for (int k = 0; k < 3; k++) {
            int32_t vi = ps->faces[f * 3 + (size_t)k];
            fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                    u[vi], axis_coordinate(cal, &ps->verts[(size_t)vi * 3]),
                    rgb[0], rgb[1], rgb[2]);
        }
        fprintf(fp, "f %zu %zu %zu\n", vertex, vertex + 1, vertex + 2);
        vertex += 3;
    }
    return fclose(fp) == 0 ? 0 : -1;
}

typedef struct {
    uint64_t key;
    size_t corner;
} AsBoxcutSplitCorner;

typedef struct {
    float *verts;
    int32_t *faces;
    double *u;
    double *v;
    float *registered_u;
    float *phi;
    int32_t *component;
    int32_t *source_vertex;
    size_t nv;
    size_t nf;
} AsBoxcutSplitMesh;

static int compare_boxcut_split_corner(const void *pa, const void *pb)
{
    const AsBoxcutSplitCorner *a = (const AsBoxcutSplitCorner *)pa;
    const AsBoxcutSplitCorner *b = (const AsBoxcutSplitCorner *)pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->corner < b->corner ? -1 : (a->corner > b->corner ? 1 : 0);
}

static int build_boxcut_split_mesh(
    Arena_T arena, const PieceSet *ps, const double *axial,
    const double *base_u, const int32_t *face_label, size_t ncharts,
    const double *chart_shift, AsBoxcutSplitMesh *out)
{
    if (arena == NULL || ps == NULL || axial == NULL || base_u == NULL ||
        face_label == NULL || ncharts == 0 || chart_shift == NULL ||
        out == NULL || ps->nf > SIZE_MAX / 3)
        return -1;
    memset(out, 0, sizeof(*out));
    size_t ncorner = 3 * ps->nf;
    AsBoxcutSplitCorner *corner = (AsBoxcutSplitCorner *)ARENA_ALLOC(
        arena, ncorner * sizeof(AsBoxcutSplitCorner));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_label[f];
        if (chart < 0 || (size_t)chart >= ncharts) return -1;
        for (int k = 0; k < 3; k++) {
            size_t ci = f * 3 + (size_t)k;
            int32_t vertex = ps->faces[ci];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            corner[ci].key = ((uint64_t)(uint32_t)chart << 32) |
                             (uint64_t)(uint32_t)vertex;
            corner[ci].corner = ci;
        }
    }
    qsort(corner, ncorner, sizeof(AsBoxcutSplitCorner),
          compare_boxcut_split_corner);

    size_t nsplit = 0;
    for (size_t i = 0; i < ncorner; i++)
        if (i == 0 || corner[i].key != corner[i - 1].key) nsplit++;
    if (nsplit == 0 || nsplit > (size_t)INT32_MAX) return -1;
    out->verts = (float *)ARENA_ALLOC(arena, nsplit * 3 * sizeof(float));
    out->faces = (int32_t *)ARENA_ALLOC(
        arena, ps->nf * 3 * sizeof(int32_t));
    out->u = (double *)ARENA_ALLOC(arena, nsplit * sizeof(double));
    out->v = (double *)ARENA_ALLOC(arena, nsplit * sizeof(double));
    out->registered_u = (float *)ARENA_ALLOC(
        arena, nsplit * sizeof(float));
    out->phi = (float *)ARENA_ALLOC(arena, nsplit * sizeof(float));
    out->component = (int32_t *)ARENA_ALLOC(
        arena, nsplit * sizeof(int32_t));
    out->source_vertex = (int32_t *)ARENA_ALLOC(
        arena, nsplit * sizeof(int32_t));

    size_t split = 0;
    for (size_t first = 0; first < ncorner;) {
        size_t last = first + 1;
        while (last < ncorner && corner[last].key == corner[first].key)
            last++;
        int32_t chart = (int32_t)(corner[first].key >> 32);
        int32_t source = (int32_t)(corner[first].key & UINT32_MAX);
        if (chart < 0 || (size_t)chart >= ncharts || source < 0 ||
            (size_t)source >= ps->nv || split > (size_t)INT32_MAX)
            return -1;
        memcpy(&out->verts[split * 3], &ps->verts[(size_t)source * 3],
               3 * sizeof(float));
        out->u[split] = base_u[source] + chart_shift[chart];
        out->v[split] = axial[source];
        out->registered_u[split] = ps->uv[(size_t)source * 2];
        out->phi[split] = ps->phi[source];
        out->component[split] = chart;
        out->source_vertex[split] = source;
        for (size_t i = first; i < last; i++)
            out->faces[corner[i].corner] = (int32_t)split;
        split++;
        first = last;
    }
    if (split != nsplit) return -1;
    out->nv = nsplit;
    out->nf = ps->nf;
    return 0;
}

enum {
    AS_BOXCUT_COLOR_NONE = 0,
    AS_BOXCUT_COLOR_XYZ = 1,
    AS_BOXCUT_COLOR_CHART = 2
};

static int write_boxcut_split_obj(const char *dir, const char *name,
                                  const AsBoxcutSplitMesh *mesh,
                                  int world, int color_mode)
{
    double lo[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double hi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
    if (color_mode == AS_BOXCUT_COLOR_XYZ) {
        for (size_t i = 0; i < mesh->nv; i++) {
            for (int d = 0; d < 3; d++) {
                double value = (double)mesh->verts[i * 3 + (size_t)d];
                if (value < lo[d]) lo[d] = value;
                if (value > hi[d]) hi[d] = value;
            }
        }
    }
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    for (size_t i = 0; i < mesh->nv; i++) {
        double x = world ? (double)mesh->verts[i * 3] : mesh->u[i];
        double y = world ? (double)mesh->verts[i * 3 + 1] : mesh->v[i];
        double z = world ? (double)mesh->verts[i * 3 + 2] : 0.0;
        if (color_mode == AS_BOXCUT_COLOR_NONE) {
            fprintf(fp, "v %.9g %.9g %.9g\n", x, y, z);
        } else {
            double rgb[3];
            if (color_mode == AS_BOXCUT_COLOR_CHART) {
                component_color(mesh->component[i], rgb);
            } else {
                for (int d = 0; d < 3; d++) {
                    double span = hi[d] - lo[d];
                    rgb[d] = span > 1e-12
                        ? ((double)mesh->verts[i * 3 + (size_t)d] - lo[d]) /
                              span
                        : 0.5;
                }
            }
            fprintf(fp, "v %.9g %.9g %.9g %.6f %.6f %.6f\n",
                    x, y, z, rgb[0], rgb[1], rgb[2]);
        }
        fprintf(fp, "vt %.9g %.9g\n", mesh->u[i], mesh->v[i]);
    }
    for (size_t f = 0; f < mesh->nf; f++) {
        int32_t a = mesh->faces[f * 3] + 1;
        int32_t b = mesh->faces[f * 3 + 1] + 1;
        int32_t c = mesh->faces[f * 3 + 2] + 1;
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

enum {
    AS_XYZ_WELD_OBJ_ATLAS = 0,
    AS_XYZ_WELD_OBJ_BRIDGE = 1,
    AS_XYZ_WELD_OBJ_COMBINED = 2
};

static int write_xyz_weld_audit_mesh_obj(
    const char *dir, const char *name, const AtlasXyzWeldAuditMesh *mesh,
    int selection)
{
    if (dir == NULL || name == NULL || mesh == NULL || mesh->faces == NULL ||
        mesh->u == NULL || mesh->v == NULL || mesh->component == NULL ||
        mesh->nfaces != mesh->base_faces + mesh->bridge_faces ||
        mesh->nvertices != 3 * mesh->nfaces ||
        selection < AS_XYZ_WELD_OBJ_ATLAS ||
        selection > AS_XYZ_WELD_OBJ_COMBINED)
        return -1;
    size_t first = selection == AS_XYZ_WELD_OBJ_BRIDGE
        ? mesh->base_faces : 0;
    size_t last = selection == AS_XYZ_WELD_OBJ_ATLAS
        ? mesh->base_faces : mesh->nfaces;
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
            "# Flat BoxCutter chart placement with XYZ-only BPA weld faces\n"
            "# Coordinates are exact audit UV positions; no display offset.\n"
            "# Atlas faces are chart-colored; weld faces are magenta.\n");
    size_t output_vertex = 1;
    int in_bridge = -1;
    for (size_t f = first; f < last; f++) {
        int bridge = f >= mesh->base_faces;
        if (bridge != in_bridge) {
            fprintf(fp, "o %s\ng %s\n",
                    bridge ? "xyz_bpa_weld" : "atlas_boxcut_charts",
                    bridge ? "xyz_bpa_weld" : "atlas_boxcut_charts");
            in_bridge = bridge;
        }
        for (int k = 0; k < 3; k++) {
            int32_t source = mesh->faces[f * 3 + (size_t)k];
            if (source < 0 || (size_t)source >= mesh->nvertices) {
                fclose(fp);
                return -1;
            }
            double rgb[3];
            if (bridge) {
                rgb[0] = 1.0; rgb[1] = 0.05; rgb[2] = 0.75;
            } else {
                component_color(mesh->component[source], rgb);
            }
            fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                    mesh->u[source], mesh->v[source],
                    rgb[0], rgb[1], rgb[2]);
            fprintf(fp, "vt %.9g %.9g\n",
                    mesh->u[source], mesh->v[source]);
        }
        fprintf(fp, "f %zu/%zu %zu/%zu %zu/%zu\n",
                output_vertex, output_vertex,
                output_vertex + 1, output_vertex + 1,
                output_vertex + 2, output_vertex + 2);
        output_vertex += 3;
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int run_saved_boxcut_weld_export(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *ps,
    const ScaffoldCalib *cal)
{
    double *base_u = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    double *axial = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    for (size_t i = 0; i < ps->nv; i++) {
        base_u[i] = ps->uv[i * 2];
        axial[i] = axis_coordinate(cal, &ps->verts[i * 3]);
    }
    PieceSet filtered_storage;
    if (build_boxcut_quality_filtered_piece_set(
            arena, cfg, ps, cal, base_u, &filtered_storage) != 0) {
        fprintf(stderr,
                "atlas_strip_scroll: BoxCutter quality filter failed\n");
        return -1;
    }
    const PieceSet *filtered = &filtered_storage;
    int32_t *face_label = NULL;
    size_t charts = 0;
    if (load_boxcut_face_labels(arena, cfg->boxcut_export_labels, filtered->nf,
                                &face_label, &charts) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot load BoxCutter labels %s\n",
                cfg->boxcut_export_labels);
        return -1;
    }
    double *shift = NULL;
    if (cfg->boxcut_shifts_from_csv != NULL) {
        if (load_chart_shifts_csv(arena, cfg->boxcut_shifts_from_csv, charts,
                                  &shift) != 0) {
            fprintf(stderr, "atlas_strip_scroll: cannot load BoxCutter shifts %s\n",
                    cfg->boxcut_shifts_from_csv);
            return -1;
        }
    } else {
        shift = (double *)ARENA_CALLOC(arena, charts, sizeof(double));
    }
    BpaBridgeGate gate;
    memset(&gate, 0, sizeof(gate));
    gate.umb_y = cal->axis_point[1];
    gate.umb_x = cal->axis_point[2];
    gate.pitch = cal->pitch;
    gate.tol = SEAM_WIND_TOL_DEFAULT_TURNS;
    gate.hard = SEAM_WIND_HARD_TOL_DEFAULT_TURNS;
    AtlasXyzWeldTopology topology;
    memset(&topology, 0, sizeof(topology));
    if (AtlasXyzWeldTopology_build(
            arena, filtered, 128.0f, 1.5f, 0.0f, 6.0f, &gate,
            &topology) != 0 ||
        topology.nfaces == 0) {
        fprintf(stderr, "atlas_strip_scroll: XYZ/BPA weld topology failed\n");
        return -1;
    }
    AtlasXyzWeldAuditStats stats;
    AtlasXyzWeldAuditMesh mesh;
    if (AtlasXyzWeldAudit_evaluate_with_mesh(
            arena, filtered, &topology, base_u, axial, face_label, charts,
            shift, &stats, &mesh) != 0) {
        fprintf(stderr, "atlas_strip_scroll: XYZ/BPA weld embedding failed\n");
        return -1;
    }
    AtlasXyzWeldConnection *connection = NULL;
    size_t nconnection = 0;
    AtlasXyzWeldConnectionStats connection_stats;
    memset(&connection_stats, 0, sizeof(connection_stats));
    if (AtlasXyzWeldTopology_collect_connections(
            arena, filtered, &topology, face_label, charts, &connection,
            &nconnection, &connection_stats) != 0) {
        fprintf(stderr,
                "atlas_strip_scroll: XYZ/BPA weld connection export failed\n");
        return -1;
    }
    if (write_xyz_weld_audit_mesh_obj(
            cfg->out_dir, "atlas_boxcut_charts_flat.obj", &mesh,
            AS_XYZ_WELD_OBJ_ATLAS) != 0 ||
        write_xyz_weld_audit_mesh_obj(
            cfg->out_dir, "xyz_bpa_weld_flat.obj", &mesh,
            AS_XYZ_WELD_OBJ_BRIDGE) != 0 ||
        write_xyz_weld_audit_mesh_obj(
            cfg->out_dir, "atlas_boxcut_plus_xyz_bpa_weld_flat.obj", &mesh,
            AS_XYZ_WELD_OBJ_COMBINED) != 0)
        return -1;
    FILE *connection_fp = as_open(
        cfg->out_dir, "xyz_bpa_weld_connections.csv", "wb");
    if (connection_fp == NULL) return -1;
    fprintf(connection_fp,
            "relation,chart0,chart1,cross_cube_edges,"
            "total_xyz_edge_length,solved_delta\n");
    for (size_t i = 0; i < nconnection; i++) {
        const AtlasXyzWeldConnection *item = &connection[i];
        double delta = shift[item->chart1] - shift[item->chart0];
        fprintf(connection_fp, "%zu,%d,%d,%zu,%.17g,%.17g\n",
                i, item->chart0, item->chart1, item->cross_cube_edges,
                item->total_xyz_edge_length, delta);
    }
    if (fclose(connection_fp) != 0) return -1;
    FILE *fp = as_open(cfg->out_dir, "export_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "{\n  \"operation\": \"saved_boxcut_xyz_bpa_weld_obj_export\",\n");
    fprintf(fp, "  \"placed_dir\": ");
    as_write_json_string(fp, cfg->placed_dir);
    fprintf(fp, ",\n  \"face_labels\": ");
    as_write_json_string(fp, cfg->boxcut_export_labels);
    fprintf(fp, ",\n  \"chart_shifts\": ");
    if (cfg->boxcut_shifts_from_csv != NULL)
        as_write_json_string(fp, cfg->boxcut_shifts_from_csv);
    else
        fputs("null", fp);
    fprintf(fp,
            ",\n  \"charts\": %zu, \"atlas_faces\": %zu, "
            "\"xyz_bpa_weld_faces\": %zu,\n"
            "  \"exact_overlap_pairs\": %zu, \"base_base_overlap_pairs\": %zu,\n"
            "  \"bridge_base_overlap_pairs\": %zu, \"bridge_bridge_overlap_pairs\": %zu,\n"
            "  \"xyz_weld_connections\": %zu, "
            "\"xyz_weld_cross_chart_edges\": %zu,\n"
            "  \"xyz_weld_edge_stretch_p95\": %.17g, "
            "\"xyz_weld_face_stretch_p95\": %.17g\n}\n",
            charts, mesh.base_faces, mesh.bridge_faces,
            stats.exact_overlap_pairs, stats.base_base_overlap_pairs,
            stats.bridge_base_overlap_pairs, stats.bridge_bridge_overlap_pairs,
            nconnection, connection_stats.cross_chart_edges,
            stats.best_cross_cube_edge_stretch_p95,
            stats.best_bridge_face_symmetric_stretch_p95);
    if (fclose(fp) != 0) return -1;
    fprintf(stderr,
            "[atlas_strip_scroll] exported BoxCutter + XYZ/BPA weld: "
            "charts=%zu atlas_faces=%zu bridge_faces=%zu exact=%zu\n",
            charts, mesh.base_faces, mesh.bridge_faces,
            stats.exact_overlap_pairs);
    return 0;
}

static int export_boxcut_split_placed(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *ps,
    const AsBoxcutSplitMesh *mesh, const char *dir)
{
    char placed[AS_PATH_CAP], source[AS_PATH_CAP], destination[AS_PATH_CAP];
    if (as_path(placed, dir, "placed") != 0 ||
        as_path(source, cfg->placed_dir, "placed_index.json") != 0 ||
        as_path(destination, placed, "placed_index.json") != 0 ||
        copy_file(source, destination) != 0)
        return -1;
    float *raw_u = NULL;
    if (load_raw_u(arena, cfg->placed_dir, ps, &raw_u) != 0) return -1;

    size_t face_cursor = 0;
    for (size_t cube = 0; cube < ps->n_cubes; cube++) {
        size_t first_face = face_cursor;
        while (face_cursor < ps->nf &&
               ps->face_cube[face_cursor] == (int32_t)cube)
            face_cursor++;
        size_t nface = face_cursor - first_face;
        if (nface == 0) return -1;
        Arena_Mark mark = Arena_save(arena);
        int32_t *local_index = (int32_t *)ARENA_ALLOC(
            arena, mesh->nv * sizeof(int32_t));
        for (size_t i = 0; i < mesh->nv; i++) local_index[i] = -1;
        int32_t *local_vertex = (int32_t *)ARENA_ALLOC(
            arena, 3 * nface * sizeof(int32_t));
        size_t nlocal = 0;
        for (size_t f = first_face; f < face_cursor; f++) {
            for (int k = 0; k < 3; k++) {
                int32_t split = mesh->faces[f * 3 + (size_t)k];
                if (split < 0 || (size_t)split >= mesh->nv) return -1;
                if (local_index[split] < 0) {
                    if (nlocal > (size_t)INT32_MAX) return -1;
                    local_index[split] = (int32_t)nlocal;
                    local_vertex[nlocal++] = split;
                }
            }
        }

        char name[96];
        snprintf(name, sizeof(name), "%s_mesh.obj", ps->ids[cube]);
        FILE *fp = as_open(placed, name, "wb");
        if (fp == NULL) return -1;
        for (size_t i = 0; i < nlocal; i++) {
            int32_t split = local_vertex[i];
            fprintf(fp, "v %.9g %.9g %.9g\n",
                    (double)mesh->verts[(size_t)split * 3],
                    (double)mesh->verts[(size_t)split * 3 + 1],
                    (double)mesh->verts[(size_t)split * 3 + 2]);
        }
        for (size_t f = first_face; f < face_cursor; f++) {
            int32_t a = local_index[mesh->faces[f * 3]] + 1;
            int32_t b = local_index[mesh->faces[f * 3 + 1]] + 1;
            int32_t c = local_index[mesh->faces[f * 3 + 2]] + 1;
            fprintf(fp, "f %d %d %d\n", a, b, c);
        }
        if (fclose(fp) != 0) return -1;

        snprintf(name, sizeof(name), "%s_facekeep.u8", ps->ids[cube]);
        fp = as_open(placed, name, "wb");
        if (fp == NULL) return -1;
        for (size_t f = 0; f < nface; f++)
            if (fputc(1, fp) == EOF) { fclose(fp); return -1; }
        if (fclose(fp) != 0) return -1;

        snprintf(name, sizeof(name), "%s_group.i32", ps->ids[cube]);
        fp = as_open(placed, name, "wb");
        if (fp == NULL) return -1;
        for (size_t i = 0; i < nlocal; i++) {
            int32_t source_vertex = mesh->source_vertex[local_vertex[i]];
            int32_t gid = ps->gid[source_vertex];
            if (fwrite(&gid, sizeof(gid), 1, fp) != 1) {
                fclose(fp); return -1;
            }
        }
        if (fclose(fp) != 0) return -1;

        snprintf(name, sizeof(name), "%s_uvphi_raw.f32", ps->ids[cube]);
        fp = as_open(placed, name, "wb");
        if (fp == NULL) return -1;
        for (size_t i = 0; i < nlocal; i++) {
            int32_t split = local_vertex[i];
            int32_t source_vertex = mesh->source_vertex[split];
            float triple[3] = {raw_u[source_vertex], (float)mesh->v[split],
                               mesh->phi[split]};
            if (fwrite(triple, sizeof(float), 3, fp) != 3) {
                fclose(fp); return -1;
            }
        }
        if (fclose(fp) != 0) return -1;

        snprintf(name, sizeof(name), "%s_uvphi.f32", ps->ids[cube]);
        fp = as_open(placed, name, "wb");
        if (fp == NULL) return -1;
        for (size_t i = 0; i < nlocal; i++) {
            int32_t split = local_vertex[i];
            float triple[3] = {(float)mesh->u[split], (float)mesh->v[split],
                               mesh->phi[split]};
            if (fwrite(triple, sizeof(float), 3, fp) != 3) {
                fclose(fp); return -1;
            }
        }
        if (fclose(fp) != 0) return -1;
        Arena_restore(arena, mark);
    }
    return face_cursor == ps->nf ? 0 : -1;
}

typedef struct {
    int32_t chart;
    size_t faces;
    double radius_median;
    double radius_min;
    double radius_max;
    double local_u_min;
    double local_u_max;
    double packed_u_min;
    double packed_u_max;
} AsBoxcutRankEntry;

static int compare_boxcut_rank_entry(const void *pa, const void *pb)
{
    const AsBoxcutRankEntry *a = (const AsBoxcutRankEntry *)pa;
    const AsBoxcutRankEntry *b = (const AsBoxcutRankEntry *)pb;
    if (a->radius_median != b->radius_median)
        return a->radius_median < b->radius_median ? -1 : 1;
    if (a->radius_min != b->radius_min)
        return a->radius_min < b->radius_min ? -1 : 1;
    return a->chart < b->chart ? -1 : (a->chart > b->chart ? 1 : 0);
}

static int run_boxcut_rank_pack(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *ps,
    const double *axial, const double *base_field,
    const double *face_radius, const int32_t *face_label, size_t ncharts,
    double average_uv_edge_length)
{
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "boxcut_rankpack") != 0) return -1;
    size_t *offset = (size_t *)ARENA_CALLOC(
        arena, ncharts + 1, sizeof(size_t));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_label[f];
        if (chart < 0 || (size_t)chart >= ncharts) return -1;
        offset[(size_t)chart + 1]++;
    }
    for (size_t i = 0; i < ncharts; i++) offset[i + 1] += offset[i];
    size_t *cursor = (size_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(size_t));
    memcpy(cursor, offset, ncharts * sizeof(size_t));
    double *radius_value = (double *)ARENA_ALLOC(
        arena, ps->nf * sizeof(double));
    for (size_t f = 0; f < ps->nf; f++)
        radius_value[cursor[face_label[f]]++] = face_radius[f];

    AsBoxcutRankEntry *rank = (AsBoxcutRankEntry *)ARENA_ALLOC(
        arena, ncharts * sizeof(AsBoxcutRankEntry));
    for (size_t i = 0; i < ncharts; i++) {
        rank[i].chart = (int32_t)i;
        rank[i].faces = offset[i + 1] - offset[i];
        if (rank[i].faces == 0) return -1;
        rank[i].radius_median = lift_median(
            &radius_value[offset[i]], rank[i].faces);
        rank[i].radius_min = DBL_MAX;
        rank[i].radius_max = -DBL_MAX;
        rank[i].local_u_min = DBL_MAX;
        rank[i].local_u_max = -DBL_MAX;
        rank[i].packed_u_min = rank[i].packed_u_max = NAN;
    }
    for (size_t f = 0; f < ps->nf; f++) {
        AsBoxcutRankEntry *r = &rank[face_label[f]];
        if (face_radius[f] < r->radius_min) r->radius_min = face_radius[f];
        if (face_radius[f] > r->radius_max) r->radius_max = face_radius[f];
        for (int k = 0; k < 3; k++) {
            double u = base_field[ps->faces[f * 3 + (size_t)k]];
            if (u < r->local_u_min) r->local_u_min = u;
            if (u > r->local_u_max) r->local_u_max = u;
        }
    }
    qsort(rank, ncharts, sizeof(AsBoxcutRankEntry),
          compare_boxcut_rank_entry);
    double gap = fmax(32.0, 4.0 * average_uv_edge_length);
    double packed_cursor = 0.0;
    double *chart_shift = (double *)ARENA_ALLOC(
        arena, ncharts * sizeof(double));
    for (size_t i = 0; i < ncharts; i++) {
        AsBoxcutRankEntry *r = &rank[i];
        double width = r->local_u_max - r->local_u_min;
        if (!isfinite(width) || width < 0.0) return -1;
        chart_shift[r->chart] = packed_cursor - r->local_u_min;
        r->packed_u_min = packed_cursor;
        r->packed_u_max = packed_cursor + width;
        packed_cursor = r->packed_u_max + gap;
    }
    double packed_span = packed_cursor - gap;

    AsBoxcutSplitMesh split;
    if (build_boxcut_split_mesh(arena, ps, axial, base_field, face_label,
                                ncharts, chart_shift, &split) != 0)
        return -1;
    AtlasOverlapAudit overlap;
    if (AtlasOverlapAudit_build(
            arena, split.faces, split.nf, split.nv, split.u, split.v,
            split.registered_u, split.phi, split.component, ncharts,
            &overlap) != 0)
        return -1;

    FILE *fp = as_open(dir, "boxcut_rank_order.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "rank,chart,faces,radius_median,radius_min,radius_max,"
                "local_u_min,local_u_max,local_u_span,packed_u_min,"
                "packed_u_max,shift\n");
    for (size_t i = 0; i < ncharts; i++) {
        const AsBoxcutRankEntry *r = &rank[i];
        fprintf(fp, "%zu,%d,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g\n",
                i, r->chart, r->faces, r->radius_median, r->radius_min,
                r->radius_max, r->local_u_min, r->local_u_max,
                r->local_u_max - r->local_u_min, r->packed_u_min,
                r->packed_u_max, chart_shift[r->chart]);
    }
    fclose(fp);
    if (write_boxcut_split_obj(dir, "atlas_boxcut_rankpack_flat.obj", &split,
                               0, AS_BOXCUT_COLOR_NONE) != 0 ||
        write_boxcut_split_obj(dir,
                               "atlas_boxcut_rankpack_flat_xyzcolor.obj",
                               &split, 0, AS_BOXCUT_COLOR_XYZ) != 0 ||
        write_boxcut_split_obj(dir,
                               "atlas_boxcut_rankpack_flat_chartcolor.obj",
                               &split, 0, AS_BOXCUT_COLOR_CHART) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_rankpack_world_uv.obj",
                               &split, 1, AS_BOXCUT_COLOR_NONE) != 0 ||
        export_boxcut_split_placed(arena, cfg, ps, &split, dir) != 0)
        return -1;
    fp = as_open(dir, "boxcut_rankpack_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n  \"operation\": \"boxcutter_radius_rank_pack_no_global_u\",\n"
        "  \"ordering\": \"ascending_median_3d_radius\",\n"
        "  \"charts\": %zu,\n  \"faces\": %zu,\n"
        "  \"vertices_after_chart_split\": %zu,\n"
        "  \"fixed_gap\": %.17g,\n  \"packed_u_span\": %.17g,\n"
        "  \"exact_overlap_pairs_after\": %zu,\n"
        "  \"within_chart_overlap_pairs_after\": %zu,\n"
        "  \"cross_chart_overlap_pairs_after\": %zu\n}\n",
        ncharts, split.nf, split.nv, gap, packed_span,
        overlap.exact_face_pairs, overlap.same_component_pairs,
        overlap.cross_component_pairs);
    fclose(fp);
    fprintf(stderr,
        "[atlas_strip_scroll] BoxCutter rank-pack: charts=%zu span=%.0f "
        "overlap=%zu within=%zu cross=%zu\n",
        ncharts, packed_span, overlap.exact_face_pairs,
        overlap.same_component_pairs, overlap.cross_component_pairs);
    return overlap.cross_component_pairs == 0 ? 0 : -1;
}

static int run_saved_boxcut_rank_pack(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *ps,
    const ScaffoldCalib *cal)
{
    FILE *fp = fopen(cfg->rank_pack_labels, "rb");
    if (fp == NULL) {
        fprintf(stderr, "atlas_strip_scroll: cannot open saved BoxCutter "
                        "labels %s\n", cfg->rank_pack_labels);
        return -1;
    }
    int32_t *face_label = (int32_t *)ARENA_ALLOC(
        arena, ps->nf * sizeof(int32_t));
    if (fread(face_label, sizeof(int32_t), ps->nf, fp) != ps->nf ||
        fgetc(fp) != EOF) {
        fclose(fp);
        fprintf(stderr, "atlas_strip_scroll: saved BoxCutter labels do not "
                        "match the source mesh (%zu faces)\n", ps->nf);
        return -1;
    }
    fclose(fp);

    int32_t maximum_label = -1;
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_label[f] < 0) return -1;
        if (face_label[f] > maximum_label) maximum_label = face_label[f];
    }
    if (maximum_label < 0) return -1;
    size_t ncharts = (size_t)maximum_label + 1;
    uint8_t *present = (uint8_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(uint8_t));
    for (size_t f = 0; f < ps->nf; f++)
        present[face_label[f]] = 1;
    for (size_t chart = 0; chart < ncharts; chart++) {
        if (!present[chart]) {
            fprintf(stderr, "atlas_strip_scroll: saved BoxCutter labels are "
                            "not dense (missing chart %zu)\n", chart);
            return -1;
        }
    }

    double *axial = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    double *base_field = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    for (size_t i = 0; i < ps->nv; i++) {
        axial[i] = axis_coordinate(cal, &ps->verts[i * 3]);
        base_field[i] = (double)ps->uv[i * 2];
    }
    double *face_radius = (double *)ARENA_ALLOC(
        arena, ps->nf * sizeof(double));
    for (size_t f = 0; f < ps->nf; f++)
        face_radius[f] = boxcut_face_radius(ps, cal, f);

    /* Audit the saved parameterization before applying the deliberately
     * destructive radius-rank diagnostic.  This is a cheap way to run the
     * production overlap predicate on a previously exported atlas, and keeps
     * validator failures separate from layout experiments. */
    double *zero_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(double));
    AsBoxcutSplitMesh input_split;
    AtlasOverlapAudit input_overlap;
    if (build_boxcut_split_mesh(arena, ps, axial, base_field, face_label,
                                ncharts, zero_shift, &input_split) != 0 ||
        AtlasOverlapAudit_build(
            arena, input_split.faces, input_split.nf, input_split.nv,
            input_split.u, input_split.v, input_split.registered_u,
            input_split.phi, input_split.component, ncharts,
            &input_overlap) != 0 ||
        !input_overlap.broad_phase_complete ||
        input_overlap.indexed_faces != input_split.nf)
        return -1;
    FILE *input_summary = as_open(
        cfg->out_dir, "saved_parameter_overlap_summary.json", "wb");
    if (input_summary == NULL) return -1;
    fprintf(input_summary,
        "{\n  \"operation\": \"saved_parameter_overlap_audit\",\n"
        "  \"faces\": %zu,\n  \"charts\": %zu,\n"
        "  \"indexed_faces\": %zu,\n"
        "  \"broad_phase_complete\": true,\n"
        "  \"broad_phase_records\": %zu,\n"
        "  \"broad_phase_cells\": %zu,\n"
        "  \"broad_phase_candidate_pairs\": %zu,\n"
        "  \"broad_phase_cell_size\": %.17g,\n"
        "  \"exact_overlap_pairs\": %zu,\n"
        "  \"within_chart_overlap_pairs\": %zu,\n"
        "  \"cross_chart_overlap_pairs\": %zu\n}\n",
        input_split.nf, ncharts, input_overlap.indexed_faces,
        input_overlap.broad_phase_records, input_overlap.broad_phase_cells,
        input_overlap.broad_phase_candidate_pairs,
        input_overlap.broad_phase_cell_size,
        input_overlap.exact_face_pairs,
        input_overlap.same_component_pairs,
        input_overlap.cross_component_pairs);
    if (fclose(input_summary) != 0) return -1;
    fprintf(stderr,
        "[atlas_strip_scroll] saved parameter audit: faces=%zu/%zu "
        "records=%zu cells=%zu candidates=%zu overlap=%zu "
        "(within=%zu cross=%zu)\n",
        input_overlap.indexed_faces, input_split.nf,
        input_overlap.broad_phase_records, input_overlap.broad_phase_cells,
        input_overlap.broad_phase_candidate_pairs,
        input_overlap.exact_face_pairs,
        input_overlap.same_component_pairs,
        input_overlap.cross_component_pairs);

    fprintf(stderr,
        "[atlas_strip_scroll] saved BoxCutter diagnostic: labels=%s "
        "charts=%zu; ordering uses only median 3-D radius\n",
        cfg->rank_pack_labels, ncharts);
    /* Passing zero deliberately makes the inter-chart gap a fixed 32 units.
     * Neither chart ordering nor spacing uses the previous global U layout. */
    return run_boxcut_rank_pack(arena, cfg, ps, axial, base_field,
                                face_radius, face_label, ncharts, 0.0);
}

typedef struct {
    double maximum_uv_edge;
    double maximum_xyz_edge;
    double maximum_edge_stretch;
    double uv_aspect;
    int reject_stretch;
    int reject_aspect;
} AsBoxcutFaceQuality;

static AsBoxcutFaceQuality boxcut_face_quality(
    const PieceSet *ps, const ScaffoldCalib *cal, const double *first_u,
    size_t face)
{
    AsBoxcutFaceQuality q;
    memset(&q, 0, sizeof(q));
    double edge2[3] = {0.0, 0.0, 0.0};
    double uv[3][2];
    for (int k = 0; k < 3; k++) {
        int32_t vertex = ps->faces[face * 3 + (size_t)k];
        uv[k][0] = first_u[vertex];
        uv[k][1] = axis_coordinate(cal, &ps->verts[(size_t)vertex * 3]);
    }
    for (int e = 0; e < 3; e++) {
        int32_t a = ps->faces[face * 3 + (size_t)e];
        int32_t b = ps->faces[face * 3 + (size_t)((e + 1) % 3)];
        double du = uv[(e + 1) % 3][0] - uv[e][0];
        double dv = uv[(e + 1) % 3][1] - uv[e][1];
        double euv = sqrt(du * du + dv * dv);
        double d0 = (double)ps->verts[(size_t)b * 3] -
                    (double)ps->verts[(size_t)a * 3];
        double d1 = (double)ps->verts[(size_t)b * 3 + 1] -
                    (double)ps->verts[(size_t)a * 3 + 1];
        double d2 = (double)ps->verts[(size_t)b * 3 + 2] -
                    (double)ps->verts[(size_t)a * 3 + 2];
        double exyz = sqrt(d0 * d0 + d1 * d1 + d2 * d2);
        edge2[e] = euv * euv;
        if (euv > q.maximum_uv_edge) q.maximum_uv_edge = euv;
        if (exyz > q.maximum_xyz_edge) q.maximum_xyz_edge = exyz;
        if (exyz > 1e-12 && euv / exyz > q.maximum_edge_stretch)
            q.maximum_edge_stretch = euv / exyz;
        if (!isfinite(euv) || !isfinite(exyz) || exyz <= 1e-12 ||
            euv > fmax(4.0 * exyz, 25.0))
            q.reject_stretch = 1;
    }
    double area = 0.5 * fabs(
        (uv[1][0] - uv[0][0]) * (uv[2][1] - uv[0][1]) -
        (uv[2][0] - uv[0][0]) * (uv[1][1] - uv[0][1]));
    q.uv_aspect = area > 1e-12
        ? (edge2[0] + edge2[1] + edge2[2]) /
              (4.0 * sqrt(3.0) * area)
        : DBL_MAX;
    if (!isfinite(q.uv_aspect) || q.uv_aspect > 100.0)
        q.reject_aspect = 1;
    return q;
}

static int build_boxcut_quality_filtered_piece_set(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *source,
    const ScaffoldCalib *cal, const double *first_u, PieceSet *filtered)
{
    if (arena == NULL || cfg == NULL || source == NULL || cal == NULL ||
        first_u == NULL || filtered == NULL)
        return -1;
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "boxcut_input_filter") != 0) return -1;
    uint8_t *keep = (uint8_t *)ARENA_ALLOC(arena, source->nf);
    double *aspect = (double *)ARENA_ALLOC(
        arena, source->nf * sizeof(double));
    size_t rejected = 0, rejected_stretch = 0, rejected_aspect = 0,
           rejected_both = 0;
    FILE *csv = as_open(dir, "rejected_faces.csv", "wb");
    FILE *world = as_open(dir, "rejected_faces_world.obj", "wb");
    FILE *flat = as_open(dir, "rejected_faces_first_uv.obj", "wb");
    if (csv == NULL || world == NULL || flat == NULL) return -1;
    fprintf(csv, "face,cube,cube_id,max_uv_edge,max_xyz_edge,"
                 "max_edge_stretch,uv_aspect,reject_stretch,reject_aspect\n");
    size_t output_vertex = 1;
    for (size_t f = 0; f < source->nf; f++) {
        AsBoxcutFaceQuality q = boxcut_face_quality(
            source, cal, first_u, f);
        aspect[f] = q.uv_aspect;
        int reject = q.reject_stretch || q.reject_aspect;
        keep[f] = (uint8_t)!reject;
        if (!reject) continue;
        rejected++;
        if (q.reject_stretch) rejected_stretch++;
        if (q.reject_aspect) rejected_aspect++;
        if (q.reject_stretch && q.reject_aspect) rejected_both++;
        int32_t cube = source->face_cube[f];
        const char *cube_id = cube >= 0 && (size_t)cube < source->n_cubes
            ? source->ids[cube] : "?";
        fprintf(csv, "%zu,%d,%s,%.17g,%.17g,%.17g,%.17g,%d,%d\n",
                f, cube, cube_id, q.maximum_uv_edge, q.maximum_xyz_edge,
                q.maximum_edge_stretch, q.uv_aspect,
                q.reject_stretch, q.reject_aspect);
        for (int k = 0; k < 3; k++) {
            int32_t vertex = source->faces[f * 3 + (size_t)k];
            const float *p = &source->verts[(size_t)vertex * 3];
            fprintf(world, "v %.9g %.9g %.9g 1 0 0\n",
                    (double)p[0], (double)p[1], (double)p[2]);
            fprintf(flat, "v %.9g %.9g 0 1 0 0\n", first_u[vertex],
                    axis_coordinate(cal, p));
        }
        fprintf(world, "f %zu %zu %zu\n", output_vertex,
                output_vertex + 1, output_vertex + 2);
        fprintf(flat, "f %zu %zu %zu\n", output_vertex,
                output_vertex + 1, output_vertex + 2);
        output_vertex += 3;
    }
    if (fclose(csv) != 0 || fclose(world) != 0 || fclose(flat) != 0)
        return -1;

    *filtered = *source;
    filtered->faces = (int32_t *)ARENA_ALLOC(
        arena, source->nf * 3 * sizeof(int32_t));
    filtered->face_cube = (int32_t *)ARENA_ALLOC(
        arena, source->nf * sizeof(int32_t));
    memcpy(filtered->faces, source->faces,
           source->nf * 3 * sizeof(int32_t));
    memcpy(filtered->face_cube, source->face_cube,
           source->nf * sizeof(int32_t));
    PieceSet_apply_facekeep(filtered, keep);
    if (filtered->nf == 0 || filtered->nf + rejected != source->nf)
        return -1;
    size_t *cube_faces = (size_t *)ARENA_CALLOC(
        arena, source->n_cubes, sizeof(size_t));
    for (size_t f = 0; f < filtered->nf; f++) {
        int32_t cube = filtered->face_cube[f];
        if (cube < 0 || (size_t)cube >= source->n_cubes) return -1;
        cube_faces[cube]++;
    }
    for (size_t cube = 0; cube < source->n_cubes; cube++)
        if (cube_faces[cube] == 0) return -1;

    qsort(aspect, source->nf, sizeof(double), compare_lift_double);
    double p99 = aspect[(size_t)(0.99 * (double)(source->nf - 1))];
    double p999 = aspect[(size_t)(0.999 * (double)(source->nf - 1))];
    FILE *summary = as_open(dir, "filter_summary.json", "wb");
    if (summary == NULL) return -1;
    fprintf(summary,
        "{\n  \"operation\": \"first_parameterization_uv_quality_filter\",\n"
        "  \"faces_before\": %zu,\n  \"faces_after\": %zu,\n"
        "  \"faces_rejected\": %zu,\n"
        "  \"stretch_rejected\": %zu,\n"
        "  \"aspect_rejected\": %zu,\n"
        "  \"rejected_by_both\": %zu,\n"
        "  \"stretch_rule\": \"uv_edge > max(4 * xyz_edge, 25)\",\n"
        "  \"aspect_rule\": \"sum(edge_squared) / (4 sqrt(3) area) > 100\",\n"
        "  \"uv_aspect_p99\": %.17g,\n"
        "  \"uv_aspect_p999\": %.17g\n}\n",
        source->nf, filtered->nf, rejected, rejected_stretch,
        rejected_aspect, rejected_both, p99, p999);
    fclose(summary);
    fprintf(stderr,
        "[atlas_strip_scroll] first-UV filter: %zu -> %zu faces "
        "(rejected=%zu stretch=%zu aspect=%zu both=%zu)\n",
        source->nf, filtered->nf, rejected, rejected_stretch,
        rejected_aspect, rejected_both);
    return 0;
}

typedef struct {
    MonotoneQpStats qp;
    size_t skipped_degenerate_triangles;
    double model_rms_before;
    double model_rms_after;
    double model_max_after;
    double correction_rms;
    double correction_max;
} AsWindingFieldStats;

typedef struct {
    int32_t island;
    int32_t correction_before;
    int32_t correction_after;
    size_t cross_overlap_pairs_before;
    size_t cross_overlap_pairs_after;
    double local_energy_change;
} AsWindingRefineMove;

typedef struct {
    size_t initial_exact_overlap_pairs;
    size_t initial_cross_overlap_pairs;
    size_t final_exact_overlap_pairs;
    size_t final_cross_overlap_pairs;
    size_t initial_overlap_charts;
    size_t head_anchor_charts;
    size_t movable_overlap_charts;
    size_t moved_charts;
    size_t moved_faces;
    size_t candidate_evaluations;
    size_t iterations;
    size_t moves;
} AsWindingRefineStats;

static double winding_model_value(const ScaffoldCalib *cal, double phi)
{
    return cal->spiral_a * phi +
           cal->spiral_b * phi * phi / (2.0 * SCAFFOLD_2PI);
}

static void fill_winding_model_island(
    const PieceSet *ps, const ScaffoldCalib *cal,
    const AsBoxcutSplitMesh *split, int32_t island, int32_t correction,
    double *model)
{
    for (size_t i = 0; i < split->nv; i++) {
        if (island >= 0 && split->component[i] != island) continue;
        int32_t source = split->source_vertex[i];
        double corrected_phi = (double)ps->phi[source] +
            (double)cal->sense * SCAFFOLD_2PI * (double)correction;
        model[i] = winding_model_value(cal, corrected_phi);
    }
}

static int winding_model_overlap_audit(
    Arena_T arena, const AsBoxcutSplitMesh *split, const double *model,
    size_t nislands, size_t *out_exact, size_t *out_cross,
    size_t *incident)
{
    Arena_Mark mark = Arena_save(arena);
    AtlasOverlapAudit audit;
    int rc = AtlasOverlapAudit_build(
        arena, split->faces, split->nf, split->nv, model, split->v,
        split->registered_u, split->phi, split->component, nislands, &audit);
    if (rc != 0 || audit.pair_buffer_truncated) {
        Arena_restore(arena, mark);
        return -1;
    }
    if (out_exact != NULL) *out_exact = audit.exact_face_pairs;
    if (out_cross != NULL) *out_cross = audit.cross_component_pairs;
    if (incident != NULL) {
        memset(incident, 0, nislands * sizeof(size_t));
        for (size_t i = 0; i < audit.npairs; i++) {
            int32_t a = audit.pairs[i].component0;
            int32_t b = audit.pairs[i].component1;
            if (a == b) continue;
            if (a < 0 || b < 0 ||
                (size_t)a >= nislands || (size_t)b >= nislands) {
                Arena_restore(arena, mark);
                return -1;
            }
            incident[a]++;
            incident[b]++;
        }
    }
    Arena_restore(arena, mark);
    return 0;
}

static double winding_local_label_energy(
    int32_t island, int32_t candidate, const int32_t *correction,
    const double *prior, const double *prior_mad,
    const size_t *island_faces,
    const AtlasWindingSyncRelation *relation, size_t nrelation)
{
    double support = sqrt((double)island_faces[island]);
    double prior_weight =
        (1.0 + fmin(8.0, support / 8.0)) /
        (1.0 + 2.0 * fmin(2.0, prior_mad[island]));
    double delta = (double)candidate - prior[island];
    double energy = prior_weight * delta * delta;
    for (size_t i = 0; i < nrelation; i++) {
        const AtlasWindingSyncRelation *r = &relation[i];
        if (!r->eligible ||
            (r->island0 != island && r->island1 != island))
            continue;
        double solved_delta =
            r->island0 == island
                ? (double)correction[r->island1] - (double)candidate
                : (double)candidate - (double)correction[r->island0];
        double residual = solved_delta -
                          (double)r->target_turn_correction;
        double square = residual * residual;
        energy += r->weight * fmin(square, 4.0);
    }
    return energy;
}

/*
 * The synchronization graph supplies a good first integer lift but is
 * intentionally local and can settle on an inconsistent cycle.  Refine only
 * islands which the exact overlap detector proves are still colliding.
 * Candidate labels stay close to the absolute radial prior; the exact
 * post-cut overlap count is the primary term and the original seam/cross-
 * section graph is a soft continuity guard.
 */
static int refine_winding_from_exact_overlaps(
    Arena_T arena, const PieceSet *ps, const ScaffoldCalib *cal,
    const AsBoxcutSplitMesh *split, size_t nislands,
    const double *prior, const double *prior_mad,
    const AtlasWindingSyncRelation *relation, size_t nrelation,
    int32_t *correction, AsWindingRefineMove *move,
    size_t move_capacity, AsWindingRefineStats *stats)
{
    if (arena == NULL || ps == NULL || cal == NULL || split == NULL ||
        nislands == 0 || prior == NULL || prior_mad == NULL ||
        relation == NULL || correction == NULL || move == NULL ||
        move_capacity == 0 || stats == NULL)
        return -1;
    memset(stats, 0, sizeof(*stats));
    double *model = (double *)ARENA_ALLOC(
        arena, split->nv * sizeof(double));
    for (size_t island = 0; island < nislands; island++)
        fill_winding_model_island(
            ps, cal, split, (int32_t)island, correction[island], model);
    size_t *incident = (size_t *)ARENA_ALLOC(
        arena, nislands * sizeof(size_t));
    size_t *island_faces = (size_t *)ARENA_CALLOC(
        arena, nislands, sizeof(size_t));
    uint8_t *selected = (uint8_t *)ARENA_ALLOC(
        arena, nislands * sizeof(uint8_t));
    for (size_t f = 0; f < split->nf; f++) {
        int32_t island = split->component[split->faces[f * 3]];
        if (island < 0 || (size_t)island >= nislands ||
            split->component[split->faces[f * 3 + 1]] != island ||
            split->component[split->faces[f * 3 + 2]] != island)
            return -1;
        island_faces[island]++;
    }
    size_t current_exact = 0, current_cross = 0;
    if (winding_model_overlap_audit(
            arena, split, model, nislands, &current_exact, &current_cross,
            incident) != 0)
        return -1;
    stats->initial_exact_overlap_pairs = current_exact;
    stats->initial_cross_overlap_pairs = current_cross;

    for (size_t iteration = 0;
         iteration < move_capacity && current_cross > 0; iteration++) {
        if (winding_model_overlap_audit(
                arena, split, model, nislands, &current_exact, &current_cross,
                incident) != 0)
            return -1;
        memset(selected, 0, nislands * sizeof(uint8_t));
        int have_best = 0;
        int32_t best_island = -1, best_correction = 0;
        size_t best_exact = current_exact, best_cross = current_cross;
        double best_gain = 0.0, best_energy_change = 0.0;
        for (int rank = 0; rank < 4; rank++) {
            int32_t island = -1;
            size_t maximum_incident = 0;
            for (size_t i = 0; i < nislands; i++) {
                if (!selected[i] && incident[i] > maximum_incident) {
                    maximum_incident = incident[i];
                    island = (int32_t)i;
                }
            }
            if (island < 0 || maximum_incident == 0) break;
            selected[island] = 1;
            int32_t center = (int32_t)lround(prior[island]);
            double current_energy = winding_local_label_energy(
                island, correction[island], correction, prior, prior_mad,
                island_faces, relation, nrelation);
            for (int32_t candidate = center - 3;
                 candidate <= center + 3; candidate++) {
                if (candidate == correction[island]) continue;
                fill_winding_model_island(
                    ps, cal, split, island, candidate, model);
                size_t candidate_exact = 0, candidate_cross = 0;
                if (winding_model_overlap_audit(
                        arena, split, model, nislands, &candidate_exact,
                        &candidate_cross, NULL) != 0)
                    return -1;
                stats->candidate_evaluations++;
                double candidate_energy = winding_local_label_energy(
                    island, candidate, correction, prior, prior_mad,
                    island_faces, relation, nrelation);
                double energy_change = candidate_energy - current_energy;
                double gain = (double)current_cross -
                              (double)candidate_cross -
                              0.25 * fmax(0.0, energy_change);
                if (gain > best_gain + 1e-9 ||
                    (fabs(gain - best_gain) <= 1e-9 &&
                     candidate_cross < best_cross) ||
                    (fabs(gain - best_gain) <= 1e-9 &&
                     candidate_cross == best_cross &&
                     energy_change < best_energy_change)) {
                    have_best = 1;
                    best_island = island;
                    best_correction = candidate;
                    best_exact = candidate_exact;
                    best_cross = candidate_cross;
                    best_gain = gain;
                    best_energy_change = energy_change;
                }
            }
            fill_winding_model_island(
                ps, cal, split, island, correction[island], model);
        }
        if (!have_best || best_gain <= 0.5 || best_cross >= current_cross)
            break;
        AsWindingRefineMove *m = &move[stats->moves++];
        m->island = best_island;
        m->correction_before = correction[best_island];
        m->correction_after = best_correction;
        m->cross_overlap_pairs_before = current_cross;
        m->cross_overlap_pairs_after = best_cross;
        m->local_energy_change = best_energy_change;
        correction[best_island] = best_correction;
        fill_winding_model_island(
            ps, cal, split, best_island, best_correction, model);
        current_exact = best_exact;
        current_cross = best_cross;
        stats->iterations = iteration + 1;
    }
    if (winding_model_overlap_audit(
            arena, split, model, nislands, &current_exact, &current_cross,
            NULL) != 0)
        return -1;
    stats->final_exact_overlap_pairs = current_exact;
    stats->final_cross_overlap_pairs = current_cross;
    return 0;
}

static int solve_winding_continuous_field(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *ps,
    const ScaffoldCalib *cal, const int32_t *island_correction,
    AsBoxcutSplitMesh *split, double **out_model,
    AsWindingFieldStats *stats)
{
    if (arena == NULL || cfg == NULL || ps == NULL || cal == NULL ||
        island_correction == NULL || split == NULL || out_model == NULL ||
        stats == NULL)
        return -1;
    memset(stats, 0, sizeof(*stats));
    double *position = (double *)ARENA_ALLOC(
        arena, split->nv * 3 * sizeof(double));
    double *model = (double *)ARENA_ALLOC(
        arena, split->nv * sizeof(double));
    AtlasFieldObservation *observation =
        (AtlasFieldObservation *)ARENA_ALLOC(
            arena, split->nv * sizeof(AtlasFieldObservation));
    for (size_t i = 0; i < split->nv; i++) {
        int32_t island = split->component[i];
        int32_t source = split->source_vertex[i];
        if (island < 0 || source < 0 || (size_t)source >= ps->nv) return -1;
        for (int d = 0; d < 3; d++)
            position[i * 3 + (size_t)d] =
                (double)split->verts[i * 3 + (size_t)d];
        double corrected_phi = (double)ps->phi[source] +
            (double)cal->sense * SCAFFOLD_2PI *
                (double)island_correction[island];
        model[i] = winding_model_value(cal, corrected_phi);
        if (!isfinite(model[i])) return -1;
        observation[i].vertex[0] = (int32_t)i;
        observation[i].vertex[1] = (int32_t)i;
        observation[i].vertex[2] = (int32_t)i;
        observation[i].bary[0] = 1.0;
        observation[i].bary[1] = 0.0;
        observation[i].bary[2] = 0.0;
        observation[i].target = model[i];
        observation[i].weight = 80.0;
        observation[i].source = (int32_t)i;
    }
    AtlasFieldProblem field;
    memset(&field, 0, sizeof(field));
    field.position = position;
    field.u0 = split->u;
    field.nvertex = split->nv;
    field.triangle = split->faces;
    field.ntriangle = split->nf;
    field.observation = observation;
    field.nobservation = split->nv;
    AtlasFieldSystem system;
    if (AtlasField_build(arena, &field, &system) != 0) return -1;
    stats->skipped_degenerate_triangles =
        system.skipped_degenerate_triangles;
    double *solution = (double *)ARENA_ALLOC(
        arena, split->nv * sizeof(double));
    memcpy(solution, split->u, split->nv * sizeof(double));
    if (MonotoneQp_solve(arena, &system.qp, &cfg->qp, solution, NULL,
                         &stats->qp) != 0)
        return -1;
    double before2 = 0.0, after2 = 0.0, correction2 = 0.0;
    for (size_t i = 0; i < split->nv; i++) {
        double before = split->u[i] - model[i];
        double after = solution[i] - model[i];
        double correction = solution[i] - split->u[i];
        before2 += before * before;
        after2 += after * after;
        correction2 += correction * correction;
        if (fabs(after) > stats->model_max_after)
            stats->model_max_after = fabs(after);
        if (fabs(correction) > stats->correction_max)
            stats->correction_max = fabs(correction);
        split->u[i] = solution[i];
    }
    stats->model_rms_before = sqrt(before2 / (double)split->nv);
    stats->model_rms_after = sqrt(after2 / (double)split->nv);
    stats->correction_rms = sqrt(correction2 / (double)split->nv);
    *out_model = model;
    return 0;
}

typedef struct {
    double *constant;
    double *linear;
    double quadratic;
    size_t *faces;
    size_t *samples;
    int32_t *sync_turn;
    double *sync_face_agreement;
} AsWindingChartGauge;

static int build_winding_chart_gauge(
    Arena_T arena, const PieceSet *ps, const ScaffoldCalib *cal,
    const double *base_field, const int32_t *face_chart, size_t ncharts,
    const int32_t *face_island, size_t nislands,
    const int32_t *island_correction, AsWindingChartGauge *out)
{
    if (arena == NULL || ps == NULL || cal == NULL || base_field == NULL ||
        face_chart == NULL || ncharts == 0 || face_island == NULL ||
        nislands == 0 || island_correction == NULL || out == NULL)
        return -1;
    memset(out, 0, sizeof(*out));
    out->constant = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(double));
    out->linear = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(double));
    out->faces = (size_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(size_t));
    out->samples = (size_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(size_t));
    out->sync_turn = (int32_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int32_t));
    out->sync_face_agreement = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(double));
    int32_t minimum_turn = INT32_MAX, maximum_turn = INT32_MIN;
    for (size_t island = 0; island < nislands; island++) {
        if (island_correction[island] < minimum_turn)
            minimum_turn = island_correction[island];
        if (island_correction[island] > maximum_turn)
            maximum_turn = island_correction[island];
    }
    int64_t turn_bins64 = (int64_t)maximum_turn -
                          (int64_t)minimum_turn + 1;
    if (turn_bins64 <= 0 || turn_bins64 > 1024 ||
        ncharts > SIZE_MAX / (size_t)turn_bins64)
        return -1;
    size_t turn_bins = (size_t)turn_bins64;
    size_t *turn_histogram = (size_t *)ARENA_CALLOC(
        arena, ncharts * turn_bins, sizeof(size_t));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_chart[f], island = face_island[f];
        if (chart < 0 || (size_t)chart >= ncharts || island < 0 ||
            (size_t)island >= nislands)
            return -1;
        int64_t bin64 = (int64_t)island_correction[island] - minimum_turn;
        if (bin64 < 0 || bin64 >= turn_bins64) return -1;
        out->faces[chart]++;
        turn_histogram[(size_t)chart * turn_bins + (size_t)bin64]++;
    }
    for (size_t chart = 0; chart < ncharts; chart++) {
        if (out->faces[chart] == 0) return -1;
        size_t cumulative = 0;
        for (size_t bin = 0; bin < turn_bins; bin++) {
            size_t count = turn_histogram[chart * turn_bins + bin];
            cumulative += count;
            if (2 * cumulative >= out->faces[chart]) {
                out->sync_turn[chart] =
                    minimum_turn + (int32_t)bin;
                out->sync_face_agreement[chart] =
                    (double)count / (double)out->faces[chart];
                break;
            }
        }
    }
    /* The accepted HEAD atlas is the absolute reference.  The synchronized
     * turn is evidence for choosing a label, not a displacement that is
     * automatically applied.  Therefore shift(0) is exactly zero. */
    double whole_turn = (double)cal->sense * SCAFFOLD_2PI;
    out->quadratic = cal->spiral_b * whole_turn * whole_turn /
                     (2.0 * SCAFFOLD_2PI);
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t chart = face_chart[f], island = face_island[f];
        if (chart < 0 || (size_t)chart >= ncharts || island < 0 ||
            (size_t)island >= nislands)
            return -1;
        for (int k = 0; k < 3; k++) {
            int32_t vertex = ps->faces[f * 3 + (size_t)k];
            if (vertex < 0 || (size_t)vertex >= ps->nv) return -1;
            double phi = (double)ps->phi[vertex];
            out->linear[chart] += cal->spiral_a * whole_turn +
                cal->spiral_b * phi * whole_turn / SCAFFOLD_2PI;
            out->samples[chart]++;
        }
    }
    for (size_t chart = 0; chart < ncharts; chart++) {
        if (out->faces[chart] == 0 || out->samples[chart] == 0) return -1;
        out->linear[chart] /= (double)out->samples[chart];
    }
    return 0;
}

static double winding_chart_shift_at(const AsWindingChartGauge *gauge,
                                     size_t chart, int32_t chart_turn)
{
    double turn = (double)chart_turn;
    return gauge->constant[chart] + gauge->linear[chart] * turn +
           gauge->quadratic * turn * turn;
}

static void fill_winding_chart_gauge(
    const double *base_field, AsBoxcutSplitMesh *split,
    const AsWindingChartGauge *gauge, int32_t chart, int32_t chart_turn)
{
    double shift = winding_chart_shift_at(gauge, (size_t)chart, chart_turn);
    for (size_t i = 0; i < split->nv; i++) {
        if (split->component[i] != chart) continue;
        int32_t source = split->source_vertex[i];
        split->u[i] = base_field[source] + shift;
    }
}

static int refine_winding_chart_gauges(
    Arena_T arena, const double *base_field, AsBoxcutSplitMesh *split,
    const AsWindingChartGauge *gauge, size_t ncharts,
    int32_t *chart_turn, AsWindingRefineMove *move, size_t move_capacity,
    uint8_t *initial_overlap, uint8_t *head_anchor,
    AsWindingRefineStats *stats)
{
    if (arena == NULL || base_field == NULL || split == NULL ||
        gauge == NULL || ncharts == 0 || chart_turn == NULL || move == NULL ||
        move_capacity == 0 || initial_overlap == NULL ||
        head_anchor == NULL || stats == NULL || ncharts > (size_t)INT32_MAX)
        return -1;
    memset(stats, 0, sizeof(*stats));
    memset(initial_overlap, 0, ncharts * sizeof(uint8_t));
    memset(head_anchor, 0, ncharts * sizeof(uint8_t));
    for (size_t chart = 0; chart < ncharts; chart++)
        fill_winding_chart_gauge(base_field, split, gauge, (int32_t)chart,
                                 chart_turn[chart]);
    size_t *incident = (size_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(size_t));
    uint8_t *selected = (uint8_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(uint8_t));
    AtlasOverlapAudit initial_audit;
    if (AtlasOverlapAudit_build(
            arena, split->faces, split->nf, split->nv, split->u, split->v,
            split->registered_u, split->phi, split->component, ncharts,
            &initial_audit) != 0 || initial_audit.pair_buffer_truncated)
        return -1;
    UnionFind overlap_group = UF_new(arena, (int32_t)ncharts);
    for (size_t i = 0; i < initial_audit.nbundles; i++) {
        const AtlasOverlapBundle *bundle = &initial_audit.bundles[i];
        int32_t a = bundle->component0, b = bundle->component1;
        if (a == b) continue;
        if (a < 0 || b < 0 || (size_t)a >= ncharts ||
            (size_t)b >= ncharts)
            return -1;
        initial_overlap[a] = initial_overlap[b] = 1;
        uf_union(&overlap_group, a, b);
    }
    int32_t *root_anchor = (int32_t *)ARENA_ALLOC(
        arena, ncharts * sizeof(int32_t));
    for (size_t i = 0; i < ncharts; i++) root_anchor[i] = -1;
    for (size_t chart = 0; chart < ncharts; chart++) {
        if (!initial_overlap[chart]) continue;
        stats->initial_overlap_charts++;
        int32_t root = uf_find(&overlap_group, (int32_t)chart);
        int32_t previous = root_anchor[root];
        if (previous < 0 || gauge->faces[chart] > gauge->faces[previous] ||
            (gauge->faces[chart] == gauge->faces[previous] &&
             chart < (size_t)previous))
            root_anchor[root] = (int32_t)chart;
    }
    for (size_t root = 0; root < ncharts; root++) {
        int32_t anchor = root_anchor[root];
        if (anchor < 0) continue;
        head_anchor[anchor] = 1;
        stats->head_anchor_charts++;
    }
    stats->movable_overlap_charts = stats->initial_overlap_charts -
                                    stats->head_anchor_charts;
    size_t current_exact = initial_audit.exact_face_pairs;
    size_t current_cross = initial_audit.cross_component_pairs;
    stats->initial_exact_overlap_pairs = current_exact;
    stats->initial_cross_overlap_pairs = current_cross;

    for (size_t iteration = 0;
         iteration < move_capacity && current_cross > 0; iteration++) {
        if (winding_model_overlap_audit(
                arena, split, split->u, ncharts, &current_exact,
                &current_cross, incident) != 0)
            return -1;
        memset(selected, 0, ncharts * sizeof(uint8_t));
        int have_best = 0;
        int32_t best_chart = -1, best_turn = 0;
        size_t best_exact = current_exact, best_cross = current_cross;
        double best_gain = 0.0, best_penalty_change = 0.0;
        for (int rank = 0; rank < 4; rank++) {
            int32_t chart = -1;
            size_t maximum_incident = 0;
            for (size_t i = 0; i < ncharts; i++) {
                if (initial_overlap[i] && !head_anchor[i] && !selected[i] &&
                    incident[i] > maximum_incident) {
                    maximum_incident = incident[i];
                    chart = (int32_t)i;
                }
            }
            if (chart < 0 || maximum_incident == 0) break;
            selected[chart] = 1;
            int32_t current_turn = chart_turn[chart];
            for (int32_t candidate = -8; candidate <= 8; candidate++) {
                if (candidate == current_turn) continue;
                fill_winding_chart_gauge(
                    base_field, split, gauge, chart, candidate);
                size_t candidate_exact = 0, candidate_cross = 0;
                if (winding_model_overlap_audit(
                        arena, split, split->u, ncharts, &candidate_exact,
                        &candidate_cross, NULL) != 0)
                    return -1;
                stats->candidate_evaluations++;
                double movement_change =
                    (double)candidate * (double)candidate -
                    (double)current_turn * (double)current_turn;
                double sync_before = (double)current_turn -
                                     (double)gauge->sync_turn[chart];
                double sync_after = (double)candidate -
                                    (double)gauge->sync_turn[chart];
                double sync_change = sync_after * sync_after -
                                     sync_before * sync_before;
                double penalty_change = movement_change + 0.2 * sync_change;
                double gain = (double)current_cross -
                              (double)candidate_cross -
                              0.05 * fmax(0.0, movement_change) -
                              0.01 * sync_change;
                if (gain > best_gain + 1e-9 ||
                    (fabs(gain - best_gain) <= 1e-9 &&
                     candidate_cross < best_cross) ||
                    (fabs(gain - best_gain) <= 1e-9 &&
                     candidate_cross == best_cross &&
                     penalty_change < best_penalty_change)) {
                    have_best = 1;
                    best_chart = chart;
                    best_turn = candidate;
                    best_exact = candidate_exact;
                    best_cross = candidate_cross;
                    best_gain = gain;
                    best_penalty_change = penalty_change;
                }
            }
            fill_winding_chart_gauge(
                base_field, split, gauge, chart, current_turn);
        }
        if (!have_best || best_gain <= 0.5 || best_cross >= current_cross)
            break;
        AsWindingRefineMove *m = &move[stats->moves++];
        m->island = best_chart;
        m->correction_before = chart_turn[best_chart];
        m->correction_after = best_turn;
        m->cross_overlap_pairs_before = current_cross;
        m->cross_overlap_pairs_after = best_cross;
        m->local_energy_change = best_penalty_change;
        chart_turn[best_chart] = best_turn;
        fill_winding_chart_gauge(
            base_field, split, gauge, best_chart, best_turn);
        current_exact = best_exact;
        current_cross = best_cross;
        stats->iterations = iteration + 1;
    }
    if (winding_model_overlap_audit(
            arena, split, split->u, ncharts, &current_exact, &current_cross,
            NULL) != 0)
        return -1;
    stats->final_exact_overlap_pairs = current_exact;
    stats->final_cross_overlap_pairs = current_cross;
    for (size_t chart = 0; chart < ncharts; chart++) {
        if (chart_turn[chart] == 0) continue;
        stats->moved_charts++;
        stats->moved_faces += gauge->faces[chart];
    }
    return 0;
}

static int build_winding_chart_model_diagnostic(
    Arena_T arena, const PieceSet *ps, const ScaffoldCalib *cal,
    const double *base_field, const AsWindingChartGauge *gauge,
    size_t ncharts, const int32_t *chart_turn, AsBoxcutSplitMesh *split,
    double **out_model, AsWindingFieldStats *stats)
{
    if (arena == NULL || ps == NULL || cal == NULL || base_field == NULL ||
        gauge == NULL || ncharts == 0 || chart_turn == NULL ||
        split == NULL || out_model == NULL || stats == NULL)
        return -1;
    memset(stats, 0, sizeof(*stats));
    double *model = (double *)ARENA_ALLOC(
        arena, split->nv * sizeof(double));
    double before2 = 0.0, after2 = 0.0, correction2 = 0.0;
    for (size_t i = 0; i < split->nv; i++) {
        int32_t chart = split->component[i];
        int32_t source = split->source_vertex[i];
        if (chart < 0 || (size_t)chart >= ncharts || source < 0) return -1;
        int32_t total_turn = chart_turn[chart];
        double original_phi = (double)ps->phi[source];
        double phi = original_phi + (double)cal->sense * SCAFFOLD_2PI *
                         (double)total_turn;
        double target_shift = winding_model_value(cal, phi) -
                              winding_model_value(cal, original_phi);
        model[i] = base_field[source] + target_shift;
        double before = base_field[source] - model[i];
        double after = split->u[i] - model[i];
        double correction = split->u[i] - base_field[source];
        before2 += before * before;
        after2 += after * after;
        correction2 += correction * correction;
        if (fabs(after) > stats->model_max_after)
            stats->model_max_after = fabs(after);
        if (fabs(correction) > stats->correction_max)
            stats->correction_max = fabs(correction);
    }
    stats->model_rms_before = sqrt(before2 / (double)split->nv);
    stats->model_rms_after = sqrt(after2 / (double)split->nv);
    stats->correction_rms = sqrt(correction2 / (double)split->nv);
    *out_model = model;
    return 0;
}

static int run_boxcut_winding_stage(
    Arena_T arena, const ScrollConfig *cfg, const PieceSet *ps,
    const ScaffoldCalib *cal, const double *axial,
    const double *base_field, const double *face_radius,
    const int32_t *face_chart, size_t ncharts,
    const AtlasBoxcutStats *box_stats)
{
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "boxcut_winding") != 0) return -1;
    double *face_axial = (double *)ARENA_ALLOC(
        arena, ps->nf * sizeof(double));
    for (size_t f = 0; f < ps->nf; f++)
        face_axial[f] = boxcut_face_mean(ps, f, axial);

    AtlasWindingSyncProblem problem;
    memset(&problem, 0, sizeof(problem));
    problem.vertices = ps->verts;
    problem.nvertices = ps->nv;
    problem.faces = ps->faces;
    problem.nfaces = ps->nf;
    problem.phi = ps->phi;
    problem.face_radius = face_radius;
    problem.face_axial = face_axial;
    problem.face_chart = face_chart;
    problem.ncharts = ncharts;
    problem.adjacency_face0 = box_stats->adjacency_face0;
    problem.adjacency_face1 = box_stats->adjacency_face1;
    problem.nintrinsic_adjacency = box_stats->intrinsic_adjacency_edges;
    problem.nadjacency = box_stats->original_graph_edges;
    problem.spiral_a = cal->spiral_a;
    problem.spiral_b = cal->spiral_b;
    problem.pitch = cal->pitch;
    problem.sense = cal->sense;
    problem.axial_bin_spacing = 8.0;
    problem.phase_bins = 64;
    int32_t *face_island = NULL, *correction = NULL;
    double *prior = NULL, *prior_mad = NULL;
    size_t nislands = 0, nrelation = 0;
    AtlasWindingSyncRelation *relation = NULL;
    AtlasWindingSyncStats sync;
    double t0 = ves_clock_sec();
    if (AtlasWindingSync_solve(
            arena, &problem, &face_island, &nislands, &correction,
            &prior, &prior_mad, &relation, &nrelation, &sync) != 0)
        return -1;
    fprintf(stderr,
        "[atlas_strip_scroll] winding sync %.2fs: charts=%zu islands=%zu "
        "relations=%zu/%zu satisfied=%.1f%% k=[%d,%d]\n",
        ves_clock_sec() - t0, ncharts, nislands,
        sync.eligible_relations, sync.relations,
        100.0 * sync.relation_satisfaction_fraction,
        sync.minimum_correction, sync.maximum_correction);

    size_t *island_faces = (size_t *)ARENA_CALLOC(
        arena, nislands, sizeof(size_t));
    int32_t *island_chart = (int32_t *)ARENA_ALLOC(
        arena, nislands * sizeof(int32_t));
    for (size_t i = 0; i < nislands; i++) island_chart[i] = -1;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t island = face_island[f];
        island_faces[island]++;
        if (island_chart[island] < 0) island_chart[island] = face_chart[f];
        else if (island_chart[island] != face_chart[f]) return -1;
    }

    double *zero_shift = (double *)ARENA_CALLOC(
        arena, ncharts, sizeof(double));
    AsBoxcutSplitMesh split;
    if (build_boxcut_split_mesh(arena, ps, axial, base_field, face_chart,
                                ncharts, zero_shift, &split) != 0)
        return -1;
    AsWindingChartGauge chart_gauge;
    if (build_winding_chart_gauge(
            arena, ps, cal, base_field, face_chart, ncharts, face_island,
            nislands, correction, &chart_gauge) != 0)
        return -1;
    int32_t *chart_turn = (int32_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(int32_t));
    const size_t refine_move_capacity = 64;
    AsWindingRefineMove *refine_move =
        (AsWindingRefineMove *)ARENA_ALLOC(
            arena, refine_move_capacity * sizeof(AsWindingRefineMove));
    uint8_t *initial_overlap = (uint8_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(uint8_t));
    uint8_t *head_anchor = (uint8_t *)ARENA_CALLOC(
        arena, ncharts, sizeof(uint8_t));
    AsWindingRefineStats refine_stats;
    t0 = ves_clock_sec();
    if (refine_winding_chart_gauges(
            arena, base_field, &split, &chart_gauge, ncharts, chart_turn,
            refine_move, refine_move_capacity, initial_overlap, head_anchor,
            &refine_stats) != 0)
        return -1;
    fprintf(stderr,
        "[atlas_strip_scroll] anchored chart refinement %.2fs: "
        "cross=%zu->%zu overlap_charts=%zu anchors=%zu moved=%zu/%zu "
        "faces=%zu evaluations=%zu\n",
        ves_clock_sec() - t0, refine_stats.initial_cross_overlap_pairs,
        refine_stats.final_cross_overlap_pairs,
        refine_stats.initial_overlap_charts, refine_stats.head_anchor_charts,
        refine_stats.moved_charts, refine_stats.moves,
        refine_stats.moved_faces,
        refine_stats.candidate_evaluations);

    int32_t *final_correction = (int32_t *)ARENA_ALLOC(
        arena, nislands * sizeof(int32_t));
    size_t relation_satisfied_after = 0, nonzero_after = 0;
    double relation_residual2_after = 0.0, prior_residual2_after = 0.0;
    int32_t minimum_correction_after = INT32_MAX;
    int32_t maximum_correction_after = INT32_MIN;
    for (size_t i = 0; i < nislands; i++) {
        int32_t chart = island_chart[i];
        if (chart < 0 || (size_t)chart >= ncharts) return -1;
        final_correction[i] = chart_turn[chart];
        double residual = (double)final_correction[i] - prior[i];
        prior_residual2_after += residual * residual;
        if (final_correction[i] != 0) nonzero_after++;
        if (final_correction[i] < minimum_correction_after)
            minimum_correction_after = final_correction[i];
        if (final_correction[i] > maximum_correction_after)
            maximum_correction_after = final_correction[i];
    }
    for (size_t i = 0; i < nrelation; i++) {
        AtlasWindingSyncRelation *r = &relation[i];
        r->solved_turn_correction =
            final_correction[r->island1] - final_correction[r->island0];
        r->final_residual = r->solved_turn_correction -
                            r->target_turn_correction;
        if (!r->eligible) continue;
        relation_residual2_after +=
            (double)r->final_residual * (double)r->final_residual;
        if (r->final_residual == 0) relation_satisfied_after++;
    }
    double relation_satisfaction_after = sync.eligible_relations
        ? (double)relation_satisfied_after /
              (double)sync.eligible_relations
        : 0.0;
    double relation_residual_rms_after = sync.eligible_relations
        ? sqrt(relation_residual2_after /
               (double)sync.eligible_relations)
        : 0.0;
    double prior_residual_rms_after =
        sqrt(prior_residual2_after / (double)nislands);
    double *model = NULL;
    AsWindingFieldStats field_stats;
    if (build_winding_chart_model_diagnostic(
            arena, ps, cal, base_field, &chart_gauge, ncharts, chart_turn,
            &split, &model, &field_stats) != 0)
        return -1;
    fprintf(stderr,
        "[atlas_strip_scroll] chart gauge lift: model_rms %.3g->%.3g "
        "gauge_rms=%.3g (HEAD differentials preserved exactly)\n",
        field_stats.model_rms_before, field_stats.model_rms_after,
        field_stats.correction_rms);

    AtlasOverlapAudit overlap;
    if (AtlasOverlapAudit_build(
            arena, split.faces, split.nf, split.nv, split.u, split.v,
            split.registered_u, split.phi, split.component, ncharts,
            &overlap) != 0)
        return -1;
    size_t relative_flips = 0, comparable = 0, degenerate = 0;
    double maximum_local_u_edge_error = 0.0;
    for (size_t f = 0; f < split.nf; f++) {
        int32_t a = split.faces[f * 3], b = split.faces[f * 3 + 1],
                c = split.faces[f * 3 + 2];
        int32_t sa = split.source_vertex[a], sb = split.source_vertex[b],
                sc = split.source_vertex[c];
        double old_det = (base_field[sb] - base_field[sa]) *
                             (split.v[c] - split.v[a]) -
                         (base_field[sc] - base_field[sa]) *
                             (split.v[b] - split.v[a]);
        double new_det = (split.u[b] - split.u[a]) *
                             (split.v[c] - split.v[a]) -
                         (split.u[c] - split.u[a]) *
                             (split.v[b] - split.v[a]);
        double edge_error_ab = fabs((split.u[b] - split.u[a]) -
                                    (base_field[sb] - base_field[sa]));
        double edge_error_ac = fabs((split.u[c] - split.u[a]) -
                                    (base_field[sc] - base_field[sa]));
        if (edge_error_ab > maximum_local_u_edge_error)
            maximum_local_u_edge_error = edge_error_ab;
        if (edge_error_ac > maximum_local_u_edge_error)
            maximum_local_u_edge_error = edge_error_ac;
        if (fabs(new_det) < 1e-10) degenerate++;
        if (fabs(old_det) >= 1e-10) {
            comparable++;
            if (old_det * new_det < 0.0) relative_flips++;
        }
    }

    FILE *fp = as_open(dir, "winding_islands.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "island,boxcut_chart,faces,prior_turn_correction,"
                "prior_mad,sync_island_turn_correction,"
                "chart_sync_turn_vote,final_chart_turn\n");
    for (size_t i = 0; i < nislands; i++)
        fprintf(fp, "%zu,%d,%zu,%.17g,%.17g,%d,%d,%d\n", i,
                island_chart[i], island_faces[i], prior[i], prior_mad[i],
                correction[i], chart_gauge.sync_turn[island_chart[i]],
                final_correction[i]);
    fclose(fp);
    fp = as_open(dir, "winding_chart_gauges.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "boxcut_chart,faces,samples,chart_sync_turn_vote,"
                "sync_face_agreement,initial_overlap,head_anchor,"
                "final_chart_turn,gauge_shift\n");
    for (size_t chart = 0; chart < ncharts; chart++)
        fprintf(fp, "%zu,%zu,%zu,%d,%.17g,%d,%d,%d,%.17g\n", chart,
                chart_gauge.faces[chart], chart_gauge.samples[chart],
                chart_gauge.sync_turn[chart],
                chart_gauge.sync_face_agreement[chart],
                initial_overlap[chart], head_anchor[chart], chart_turn[chart],
                winding_chart_shift_at(&chart_gauge, chart,
                                       chart_turn[chart]));
    fclose(fp);
    fp = as_open(dir, "winding_overlap_refinement.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "move,boxcut_chart,chart_turn_before,chart_turn_after,"
                "cross_overlap_pairs_before,cross_overlap_pairs_after,"
                "turn_penalty_change\n");
    for (size_t i = 0; i < refine_stats.moves; i++) {
        const AsWindingRefineMove *m = &refine_move[i];
        fprintf(fp, "%zu,%d,%d,%d,%zu,%zu,%.17g\n", i, m->island,
                m->correction_before, m->correction_after,
                m->cross_overlap_pairs_before,
                m->cross_overlap_pairs_after, m->local_energy_change);
    }
    fclose(fp);
    fp = as_open(dir, "winding_relations.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "relation,island0,island1,observations,cross_section,"
                "seam,target,mode_agreement,residual_median,weight,eligible,"
                "solved_delta,final_residual\n");
    for (size_t i = 0; i < nrelation; i++) {
        const AtlasWindingSyncRelation *r = &relation[i];
        fprintf(fp, "%zu,%d,%d,%zu,%zu,%zu,%d,%.17g,%.17g,%.17g,%d,%d,%d\n",
                i, r->island0, r->island1, r->observations,
                r->cross_section_observations, r->seam_observations,
                r->target_turn_correction, r->mode_agreement,
                r->residual_median, r->weight, r->eligible,
                r->solved_turn_correction, r->final_residual);
    }
    fclose(fp);
    fp = as_open(dir, "boxcut_face_labels.i32", "wb");
    if (fp == NULL ||
        fwrite(face_chart, sizeof(int32_t), ps->nf, fp) != ps->nf ||
        fclose(fp) != 0)
        return -1;
    fp = as_open(dir, "winding_face_islands.i32", "wb");
    if (fp == NULL ||
        fwrite(face_island, sizeof(int32_t), ps->nf, fp) != ps->nf ||
        fclose(fp) != 0)
        return -1;

    AsBoxcutSplitMesh model_mesh = split;
    model_mesh.u = model;
    if (write_boxcut_chart_obj(
            dir, "atlas_head_reference_filtered_chartcolor.obj", ps, cal,
            base_field, face_chart) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_winding_model_flat.obj",
                               &model_mesh, 0, AS_BOXCUT_COLOR_NONE) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_winding_flat.obj", &split,
                               0, AS_BOXCUT_COLOR_NONE) != 0 ||
        write_boxcut_split_obj(dir,
                               "atlas_boxcut_winding_flat_xyzcolor.obj",
                               &split, 0, AS_BOXCUT_COLOR_XYZ) != 0 ||
        write_boxcut_split_obj(dir,
                               "atlas_boxcut_winding_flat_chartcolor.obj",
                               &split, 0, AS_BOXCUT_COLOR_CHART) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_winding_world_uv.obj",
                               &split, 1, AS_BOXCUT_COLOR_NONE) != 0 ||
        export_boxcut_split_placed(arena, cfg, ps, &split, dir) != 0)
        return -1;

    fp = as_open(dir, "boxcut_winding_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n  \"operation\": \"head_anchored_post_multicut_integer_winding\",\n"
        "  \"boxcut_charts\": %zu,\n  \"intrinsic_islands\": %zu,\n"
        "  \"faces\": %zu,\n  \"vertices_after_chart_split\": %zu,\n"
        "  \"pitch\": %.17g,\n  \"phase_bins\": %d,\n"
        "  \"axial_bin_spacing\": %.17g,\n"
        "  \"cylindrical_bins\": %zu,\n"
        "  \"strand_samples\": %zu,\n"
        "  \"winding_observations\": %zu,\n"
        "  \"cross_section_observations\": %zu,\n"
        "  \"seam_observations\": %zu,\n"
        "  \"relations\": %zu,\n  \"eligible_relations\": %zu,\n"
        "  \"islands_with_relations\": %zu,\n"
        "  \"relation_graph_components\": %zu,\n"
        "  \"relation_satisfaction_fraction_before_refine\": %.17g,\n"
        "  \"relation_satisfaction_fraction_after_refine\": %.17g,\n"
        "  \"chart_exact_overlap_pairs_before_refine\": %zu,\n"
        "  \"chart_cross_overlap_pairs_before_refine\": %zu,\n"
        "  \"chart_exact_overlap_pairs_after_refine\": %zu,\n"
        "  \"chart_cross_overlap_pairs_after_refine\": %zu,\n"
        "  \"initial_overlap_charts\": %zu,\n"
        "  \"head_anchor_charts\": %zu,\n"
        "  \"movable_overlap_charts\": %zu,\n"
        "  \"moved_charts\": %zu,\n"
        "  \"moved_faces\": %zu,\n"
        "  \"whole_chart_refinement_moves\": %zu,\n"
        "  \"whole_chart_refinement_candidate_evaluations\": %zu,\n"
        "  \"turn_correction_min\": %d,\n"
        "  \"turn_correction_max\": %d,\n"
        "  \"nonzero_turn_corrections\": %zu,\n"
        "  \"prior_residual_rms_before_refine\": %.17g,\n"
        "  \"prior_residual_rms_after_refine\": %.17g,\n"
        "  \"relation_residual_rms_before_refine\": %.17g,\n"
        "  \"relation_residual_rms_after_refine\": %.17g,\n"
        "  \"relative_winding_model_rms_before_gauge\": %.17g,\n"
        "  \"relative_winding_model_rms_after_gauge\": %.17g,\n"
        "  \"chart_gauge_shift_rms\": %.17g,\n"
        "  \"chart_gauge_shift_max\": %.17g,\n"
        "  \"local_head_differentials_preserved\": true,\n"
        "  \"maximum_local_u_edge_error\": %.17g,\n"
        "  \"relative_flips\": %zu,\n"
        "  \"comparable_faces\": %zu,\n"
        "  \"degenerate_faces\": %zu,\n"
        "  \"exact_overlap_pairs_after\": %zu,\n"
        "  \"within_chart_overlap_pairs_after\": %zu,\n"
        "  \"cross_chart_overlap_pairs_after\": %zu\n}\n",
        ncharts, nislands, split.nf, split.nv, cal->pitch,
        problem.phase_bins, problem.axial_bin_spacing,
        sync.cylindrical_bins, sync.strand_samples, sync.observations,
        sync.cross_section_observations, sync.seam_observations,
        sync.relations, sync.eligible_relations,
        sync.islands_with_relations, sync.relation_graph_components,
        sync.relation_satisfaction_fraction, relation_satisfaction_after,
        refine_stats.initial_exact_overlap_pairs,
        refine_stats.initial_cross_overlap_pairs,
        refine_stats.final_exact_overlap_pairs,
        refine_stats.final_cross_overlap_pairs,
        refine_stats.initial_overlap_charts, refine_stats.head_anchor_charts,
        refine_stats.movable_overlap_charts, refine_stats.moved_charts,
        refine_stats.moved_faces, refine_stats.moves,
        refine_stats.candidate_evaluations, minimum_correction_after,
        maximum_correction_after, nonzero_after,
        sync.prior_residual_rms, prior_residual_rms_after,
        sync.final_relation_residual_rms, relation_residual_rms_after,
        field_stats.model_rms_before, field_stats.model_rms_after,
        field_stats.correction_rms, field_stats.correction_max,
        maximum_local_u_edge_error, relative_flips, comparable, degenerate,
        overlap.exact_face_pairs,
        overlap.same_component_pairs, overlap.cross_component_pairs);
    fclose(fp);
    fprintf(stderr,
        "[atlas_strip_scroll] winding atlas: overlap=%zu "
        "(within=%zu cross=%zu) flips=%zu local_du_error=%.3g output=%s\n",
        overlap.exact_face_pairs, overlap.same_component_pairs,
        overlap.cross_component_pairs, relative_flips,
        maximum_local_u_edge_error, dir);
    return 0;
}

static int run_boxcut_rigid_stage(Arena_T arena, const ScrollConfig *cfg,
                                  const PieceSet *ps,
                                  const ScaffoldCalib *cal,
                                  const AtlasCandidateSet *set,
                                  const double *first_field,
                                  const double *base_field)
{
    if (first_field == NULL || base_field == NULL) return -1;
    PieceSet filtered;
    if (build_boxcut_quality_filtered_piece_set(
            arena, cfg, ps, cal, first_field, &filtered) != 0)
        return -1;
    ps = &filtered;
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "boxcut_rigid") != 0) return -1;
    double *axial = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    float *registered = (float *)ARENA_ALLOC(arena, ps->nv * sizeof(float));
    for (size_t i = 0; i < ps->nv; i++) {
        axial[i] = axis_coordinate(cal, &ps->verts[i * 3]);
        registered[i] = ps->uv[i * 2];
    }

    AtlasOverlapAudit overlap_before;
    if (AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, base_field, axial,
            registered, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &overlap_before) != 0)
        return -1;
    uint8_t *allowed_delamination_pair = (uint8_t *)ARENA_CALLOC(
        arena, (overlap_before.npairs ? overlap_before.npairs : 1),
        sizeof(uint8_t));
    uint64_t *allowed_delamination_key = (uint64_t *)ARENA_ALLOC(
        arena, (overlap_before.npairs ? overlap_before.npairs : 1) *
                   sizeof(uint64_t));
    size_t nallowed_delamination_pair = 0;
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        if (!boxcut_overlap_pair_is_delamination(
                cfg, ps, set, axial, &overlap_before.pairs[i]))
            continue;
        allowed_delamination_pair[i] = 1;
        allowed_delamination_key[nallowed_delamination_pair++] =
            boxcut_face_pair_key(overlap_before.pairs[i].face0,
                                 overlap_before.pairs[i].face1);
    }
    qsort(allowed_delamination_key, nallowed_delamination_pair,
          sizeof(uint64_t), boxcut_compare_u64);
    AtlasSeamAudit seam;
    if (AtlasSeamAudit_build(arena, ps, base_field,
                             set->vertex_mesh_component,
                             set->mesh_components, &seam) != 0)
        return -1;
    uint8_t *eligible_seam = (uint8_t *)ARENA_CALLOC(
        arena, (seam.nbundles ? seam.nbundles : 1), sizeof(uint8_t));
    size_t neligible_seam_bundle = 0;
    for (size_t i = 0; i < seam.nbundles; i++) {
        eligible_seam[i] =
            (uint8_t)smooth_seam_bundle_eligible(&seam.bundles[i]);
        if (eligible_seam[i]) neligible_seam_bundle++;
    }

    /*
     * Partition fail-closed: an overlap is not presumed to be a delamination
     * merely because the two face centroids have similar radius.  Real
     * delaminations require the explicit bounded-front whitelist built above.
     * Keep the positive gate only for the softer radial-order observations
     * used by the continuous placement QP.
     */
    const double partition_radius_gate = 0.0;
    const double layout_radius_gate =
        cal->pitch > 0.0 ? 0.5 * cal->pitch : 4.75;
    AtlasBoxcutProblem box_problem;
    memset(&box_problem, 0, sizeof(box_problem));
    box_problem.vertices = ps->verts;
    box_problem.nvertices = ps->nv;
    box_problem.faces = ps->faces;
    box_problem.nfaces = ps->nf;
    box_problem.u = base_field;
    box_problem.v = axial;
    box_problem.axis_point = cal->axis_point;
    box_problem.axis_dir = cal->axis_dir;
    box_problem.overlap = &overlap_before;
    box_problem.allowed_overlap_pair = allowed_delamination_pair;
    box_problem.seam = &seam;
    box_problem.eligible_seam_bundle = eligible_seam;
    box_problem.radius_gate = partition_radius_gate;
    int32_t *face_label = NULL;
    AtlasBoxcutStats box_stats;
    double t0 = ves_clock_sec();
    int box_rc = AtlasBoxcut_partition(
        arena, &box_problem, &face_label, &box_stats);
    if (box_rc != 0 || face_label == NULL || box_stats.charts == 0 ||
        box_stats.repulsive_edges_joined != 0) {
        fprintf(stderr,
            "[atlas_strip_scroll] BoxCutter partition rejected: rc=%d "
            "charts=%zu repulsive=%zu cut=%zu joined=%zu ignored=%zu\n",
            box_rc, box_stats.charts, box_stats.radial_repulsive_edges,
            box_stats.repulsive_edges_cut,
            box_stats.repulsive_edges_joined,
            box_stats.near_radius_pairs_ignored);
        return -1;
    }
    size_t multicut_charts = box_stats.charts;
    size_t layout_blocks = multicut_charts;
    if (cfg->boxcut_component_blocks &&
        build_boxcut_component_face_labels(
            arena, ps, set, &face_label, &layout_blocks) != 0) {
        fprintf(stderr,
                "[atlas_strip_scroll] cannot build retained-component "
                "BoxCutter block labels\n");
        return -1;
    }
    size_t within_layout_block_overlap_pairs_before = 0;
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        const AtlasOverlapPair *pair = &overlap_before.pairs[i];
        if (face_label[pair->face0] == face_label[pair->face1])
            within_layout_block_overlap_pairs_before++;
    }
    size_t cross_layout_block_overlap_pairs_before =
        overlap_before.npairs - within_layout_block_overlap_pairs_before;

    double *face_radius = (double *)ARENA_ALLOC(
        arena, ps->nf * sizeof(double));
    for (size_t f = 0; f < ps->nf; f++)
        face_radius[f] = boxcut_face_radius(ps, cal, f);
    AtlasBoxcutLayoutProblem layout_problem;
    memset(&layout_problem, 0, sizeof(layout_problem));
    layout_problem.faces = ps->faces;
    layout_problem.nfaces = ps->nf;
    layout_problem.nvertices = ps->nv;
    layout_problem.base_u = base_field;
    layout_problem.v = axial;
    layout_problem.registered_u = registered;
    layout_problem.phi = ps->phi;
    layout_problem.face_radius = face_radius;
    layout_problem.face_label = face_label;
    layout_problem.ncharts = layout_blocks;
    layout_problem.adjacency_face0 = box_stats.adjacency_face0;
    layout_problem.adjacency_face1 = box_stats.adjacency_face1;
    layout_problem.adjacency_weight = box_stats.adjacency_weight;
    layout_problem.intrinsic_adjacency_edges =
        box_stats.intrinsic_adjacency_edges;
    layout_problem.adjacency_edges = box_stats.original_graph_edges;
    layout_problem.overlap = &overlap_before;
    layout_problem.allowed_overlap_key = allowed_delamination_key;
    layout_problem.allowed_overlap_keys = nallowed_delamination_pair;
    layout_problem.radius_gate = layout_radius_gate;
    layout_problem.average_uv_edge_length = box_stats.average_uv_edge_length;
    layout_problem.coherence_weight_scale = 64.0;
    double *chart_shift = NULL;
    AtlasBoxcutLayoutRelation *relation = NULL;
    size_t nrelation = 0;
    AtlasBoxcutLayoutCoherence *coherence = NULL;
    size_t ncoherence = 0;
    AtlasBoxcutLayoutDisjunction *disjunction = NULL;
    size_t ndisjunction = 0;
    AtlasBoxcutLayoutStats layout_stats;
    memset(&layout_stats, 0, sizeof(layout_stats));
    int external_layout = cfg->boxcut_shifts_from_csv != NULL;
    if (external_layout) {
        if (load_chart_shifts_csv(
                arena, cfg->boxcut_shifts_from_csv, layout_blocks,
                &chart_shift) != 0) {
            fprintf(stderr,
                    "[atlas_strip_scroll] cannot load external BoxCutter "
                    "shifts %s\n", cfg->boxcut_shifts_from_csv);
            return -1;
        }
        uint8_t *chart_moved = (uint8_t *)ARENA_CALLOC(
            arena, layout_blocks, sizeof(uint8_t));
        for (size_t chart = 0; chart < layout_blocks; chart++) {
            if (fabs(chart_shift[chart]) > 1.0e-10) {
                chart_moved[chart] = 1;
                layout_stats.moved_charts++;
            }
            if (fabs(chart_shift[chart]) > layout_stats.maximum_abs_shift)
                layout_stats.maximum_abs_shift = fabs(chart_shift[chart]);
        }
        for (size_t f = 0; f < ps->nf; f++)
            if (chart_moved[face_label[f]]) layout_stats.moved_faces++;
    } else if (AtlasBoxcutLayout_solve(
                   arena, &layout_problem, &chart_shift, &relation,
                   &nrelation, &coherence, &ncoherence, &disjunction,
                   &ndisjunction, &layout_stats) != 0) {
        return -1;
    }

    AsBoxcutSplitMesh split;
    if (build_boxcut_split_mesh(arena, ps, axial, base_field, face_label,
                                layout_blocks, chart_shift, &split) != 0)
        return -1;
    AtlasOverlapAudit overlap_after;
    if (AtlasOverlapAudit_build(
            arena, split.faces, split.nf, split.nv, split.u, split.v,
            split.registered_u, split.phi, split.component,
            layout_blocks, &overlap_after) != 0)
        return -1;
    size_t final_allowed_delamination_pairs = 0;
    for (size_t i = 0; i < overlap_after.npairs; i++) {
        uint64_t key = boxcut_face_pair_key(overlap_after.pairs[i].face0,
                                            overlap_after.pairs[i].face1);
        if (boxcut_has_u64(allowed_delamination_key,
                           nallowed_delamination_pair, key))
            final_allowed_delamination_pairs++;
    }
    AtlasBoxcutLayoutDisjunction *remaining_disjunction = NULL;
    size_t nremaining_disjunction = 0;
    if (AtlasBoxcutLayout_collect_disjunctions(
            arena, &layout_problem, &overlap_after,
            &remaining_disjunction, &nremaining_disjunction) != 0)
        return -1;

    /* Held-out physical oracle.  The bridge topology is generated once from
     * XYZ and the original BPA seam-weld only; neither the UV field nor the
     * chart labels can influence which faces are proposed.  We then score the
     * same frozen topology against the pre-placement and final block layouts. */
    AtlasXyzWeldAuditStats xyz_weld_initial, xyz_weld_final;
    memset(&xyz_weld_initial, 0, sizeof(xyz_weld_initial));
    memset(&xyz_weld_final, 0, sizeof(xyz_weld_final));
    AtlasXyzWeldTopology xyz_weld_topology;
    memset(&xyz_weld_topology, 0, sizeof(xyz_weld_topology));
    AtlasXyzWeldConnection *xyz_weld_connection = NULL;
    size_t nxyz_weld_connection = 0;
    AtlasXyzWeldConnectionStats xyz_weld_connection_stats;
    memset(&xyz_weld_connection_stats, 0, sizeof(xyz_weld_connection_stats));
    int xyz_weld_connection_rc = -1;
    int xyz_weld_oracle_rc = -1;
    Arena_T xyz_topology_arena = Arena_new();
    BpaBridgeGate xyz_gate;
    memset(&xyz_gate, 0, sizeof(xyz_gate));
    xyz_gate.umb_y = cal->axis_point[1];
    xyz_gate.umb_x = cal->axis_point[2];
    xyz_gate.pitch = cal->pitch;
    xyz_gate.tol = SEAM_WIND_TOL_DEFAULT_TURNS;
    xyz_gate.hard = SEAM_WIND_HARD_TOL_DEFAULT_TURNS;
    if (AtlasXyzWeldTopology_build(
            xyz_topology_arena, ps, 128.0f, 1.5f, 0.0f, 6.0f,
            &xyz_gate, &xyz_weld_topology) == 0 &&
        xyz_weld_topology.nfaces > 0) {
        xyz_weld_connection_rc = AtlasXyzWeldTopology_collect_connections(
            xyz_topology_arena, ps, &xyz_weld_topology, face_label,
            layout_blocks, &xyz_weld_connection, &nxyz_weld_connection,
            &xyz_weld_connection_stats);
        Arena_T xyz_initial_arena = Arena_new();
        double *zero_shift = (double *)ARENA_CALLOC(
            xyz_initial_arena, layout_blocks, sizeof(double));
        int initial_rc = AtlasXyzWeldAudit_evaluate(
            xyz_initial_arena, ps, &xyz_weld_topology, base_field, axial,
            face_label, layout_blocks, zero_shift, &xyz_weld_initial);
        Arena_dispose(&xyz_initial_arena);
        Arena_T xyz_final_arena = Arena_new();
        int final_rc = AtlasXyzWeldAudit_evaluate(
            xyz_final_arena, ps, &xyz_weld_topology, base_field, axial,
            face_label, layout_blocks, chart_shift, &xyz_weld_final);
        Arena_dispose(&xyz_final_arena);
        xyz_weld_oracle_rc = initial_rc == 0 && final_rc == 0 ? 0 : -1;
    }

    FILE *xyz_fp = as_open(dir, "xyz_bpa_weld_oracle_summary.json", "wb");
    if (xyz_fp == NULL) {
        Arena_dispose(&xyz_topology_arena);
        return -1;
    }
#define AS_WRITE_XYZ_WELD_STATS(name, s, comma) do {                       \
        fprintf(xyz_fp, "  \"%s\": {\n", (name));                        \
        fprintf(xyz_fp, "    \"bridge_faces\": %zu, \"bridge_vertices\": %zu, \"bridge_edges\": %zu, \"cross_cube_bridge_edges\": %zu,\n", \
                (s).bridge_faces, (s).bridge_vertices, (s).bridge_edges,    \
                (s).cross_cube_bridge_edges);                              \
        fprintf(xyz_fp, "    \"incident_chart_positions\": %zu, \"split_bridge_vertices\": %zu,\n", \
                (s).incident_chart_positions, (s).split_bridge_vertices);  \
        fprintf(xyz_fp, "    \"incident_u_spread\": {\"median\": %.17g, \"p95\": %.17g, \"max\": %.17g},\n", \
                (s).incident_u_spread_median, (s).incident_u_spread_p95,   \
                (s).incident_u_spread_max);                                \
        fprintf(xyz_fp, "    \"cross_cube_edge_stretch\": {\"median\": %.17g, \"p95\": %.17g, \"max\": %.17g, \"exploded\": %zu, \"compressed\": %zu},\n", \
                (s).best_cross_cube_edge_stretch_median,                   \
                (s).best_cross_cube_edge_stretch_p95,                      \
                (s).best_cross_cube_edge_stretch_max,                      \
                (s).exploded_cross_cube_edges,                             \
                (s).compressed_cross_cube_edges);                          \
        fprintf(xyz_fp, "    \"cross_cube_edge_abs_error\": {\"median\": %.17g, \"p95\": %.17g, \"max\": %.17g},\n", \
                (s).best_cross_cube_edge_abs_error_median,                 \
                (s).best_cross_cube_edge_abs_error_p95,                    \
                (s).best_cross_cube_edge_abs_error_max);                   \
        fprintf(xyz_fp, "    \"bridge_face_symmetric_stretch\": {\"median\": %.17g, \"p95\": %.17g, \"max\": %.17g, \"exploded\": %zu, \"degenerate\": %zu},\n", \
                (s).best_bridge_face_symmetric_stretch_median,             \
                (s).best_bridge_face_symmetric_stretch_p95,                \
                (s).best_bridge_face_symmetric_stretch_max,                \
                (s).exploded_bridge_faces, (s).degenerate_bridge_faces);   \
        fprintf(xyz_fp, "    \"selected_embedding\": {\"vertices\": %zu, \"discontinuous_vertices\": %zu, \"u_spread_median\": %.17g, \"u_spread_p95\": %.17g, \"u_spread_max\": %.17g},\n", \
                (s).embedded_bridge_vertices,                              \
                (s).discontinuous_embedded_bridge_vertices,               \
                (s).embedded_u_spread_median,                              \
                (s).embedded_u_spread_p95, (s).embedded_u_spread_max);     \
        fprintf(xyz_fp, "    \"overlap\": {\"exact\": %zu, \"base_base\": %zu, \"bridge_base\": %zu, \"bridge_bridge\": %zu},\n", \
                (s).exact_overlap_pairs, (s).base_base_overlap_pairs,      \
                (s).bridge_base_overlap_pairs,                             \
                (s).bridge_bridge_overlap_pairs);                          \
        fprintf(xyz_fp, "    \"broad_phase\": {\"complete\": %s, \"records\": %zu, \"cells\": %zu, \"candidate_pairs\": %zu, \"cell_size\": %.17g}\n", \
                (s).broad_phase_complete ? "true" : "false",             \
                (s).broad_phase_records, (s).broad_phase_cells,            \
                (s).broad_phase_candidate_pairs,                           \
                (s).broad_phase_cell_size);                                \
        fprintf(xyz_fp, "  }%s\n", (comma));                              \
    } while (0)
    fprintf(xyz_fp,
        "{\n  \"operation\": \"held_out_xyz_bpa_weld_oracle\",\n"
        "  \"topology_uses_uv\": false,\n"
        "  \"oracle_rc\": %d,\n"
        "  \"connection_rc\": %d,\n"
        "  \"xyz_bridge_faces\": %zu,\n",
        xyz_weld_oracle_rc, xyz_weld_connection_rc,
        xyz_weld_topology.nfaces);
    fprintf(xyz_fp,
        "  \"connection_graph\": {\"cross_cube_edges\": %zu, \"same_chart_edges\": %zu, \"cross_chart_edges\": %zu, \"ambiguous_chart_edges\": %zu, \"relations\": %zu},\n",
        xyz_weld_connection_stats.cross_cube_edges,
        xyz_weld_connection_stats.same_chart_edges,
        xyz_weld_connection_stats.cross_chart_edges,
        xyz_weld_connection_stats.ambiguous_chart_edges,
        xyz_weld_connection_stats.relations);
    AS_WRITE_XYZ_WELD_STATS("initial_layout", xyz_weld_initial, ",");
    AS_WRITE_XYZ_WELD_STATS("final_layout", xyz_weld_final, "");
    fprintf(xyz_fp, "}\n");
#undef AS_WRITE_XYZ_WELD_STATS
    if (fclose(xyz_fp) != 0) {
        Arena_dispose(&xyz_topology_arena);
        return -1;
    }
    FILE *xyz_connection_fp = as_open(
        dir, "xyz_bpa_weld_connections.csv", "wb");
    if (xyz_connection_fp == NULL) {
        Arena_dispose(&xyz_topology_arena);
        return -1;
    }
    fprintf(xyz_connection_fp,
            "relation,chart0,chart1,cross_cube_edges,total_xyz_edge_length,solved_delta\n");
    for (size_t i = 0; i < nxyz_weld_connection; i++) {
        const AtlasXyzWeldConnection *connection = &xyz_weld_connection[i];
        double delta = chart_shift[connection->chart1] -
                       chart_shift[connection->chart0];
        fprintf(xyz_connection_fp, "%zu,%d,%d,%zu,%.17g,%.17g\n",
                i, connection->chart0, connection->chart1,
                connection->cross_cube_edges,
                connection->total_xyz_edge_length, delta);
    }
    if (fclose(xyz_connection_fp) != 0) {
        Arena_dispose(&xyz_topology_arena);
        return -1;
    }
    Arena_dispose(&xyz_topology_arena);
    fprintf(stderr,
        "[atlas_strip_scroll] XYZ-BPA weld oracle: rc=%d bridge_faces=%zu "
        "edge_p95=%.3g face_p95=%.3g split=%zu discontinuous=%zu "
        "new_overlap=%zu\n",
        xyz_weld_oracle_rc, xyz_weld_topology.nfaces,
        xyz_weld_final.best_cross_cube_edge_stretch_p95,
        xyz_weld_final.best_bridge_face_symmetric_stretch_p95,
        xyz_weld_final.split_bridge_vertices,
        xyz_weld_final.discontinuous_embedded_bridge_vertices,
        xyz_weld_final.bridge_base_overlap_pairs +
            xyz_weld_final.bridge_bridge_overlap_pairs);

    size_t relative_flips = 0, comparable = 0, degenerate = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a0 = ps->faces[f * 3], b0 = ps->faces[f * 3 + 1],
                c0 = ps->faces[f * 3 + 2];
        int32_t a1 = split.faces[f * 3], b1 = split.faces[f * 3 + 1],
                c1 = split.faces[f * 3 + 2];
        double d0 = (base_field[b0] - base_field[a0]) *
                        (axial[c0] - axial[a0]) -
                    (base_field[c0] - base_field[a0]) *
                        (axial[b0] - axial[a0]);
        double d1 = (split.u[b1] - split.u[a1]) *
                        (split.v[c1] - split.v[a1]) -
                    (split.u[c1] - split.u[a1]) *
                        (split.v[b1] - split.v[a1]);
        if (fabs(d1) < 1e-10) degenerate++;
        if (fabs(d0) >= 1e-10) {
            comparable++;
            if (d0 * d1 < 0.0) relative_flips++;
        }
    }
    size_t seam_edges_cut = 0;
    for (size_t i = box_stats.intrinsic_adjacency_edges;
         i < box_stats.original_graph_edges; i++)
        if (face_label[box_stats.adjacency_face0[i]] !=
            face_label[box_stats.adjacency_face1[i]])
            seam_edges_cut++;

    FILE *fp = as_open(dir, "boxcut_face_labels.i32", "wb");
    if (fp == NULL ||
        fwrite(face_label, sizeof(int32_t), ps->nf, fp) != ps->nf ||
        fclose(fp) != 0)
        return -1;
    fp = as_open(dir, "boxcut_chart_shifts.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "chart,shift\n");
    for (size_t i = 0; i < layout_blocks; i++)
        fprintf(fp, "%zu,%.17g\n", i, chart_shift[i]);
    fclose(fp);
    fp = as_open(dir, "boxcut_layout_remaining_disjunctions.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "disjunction,chart0,chart1,overlap_pairs,"
                "chart1_left_upper,chart1_right_lower\n");
    for (size_t i = 0; i < nremaining_disjunction; i++) {
        const AtlasBoxcutLayoutDisjunction *d = &remaining_disjunction[i];
        fprintf(fp, "%zu,%d,%d,%zu,%.17g,%.17g\n", i, d->chart0,
                d->chart1, d->overlap_pairs, d->chart1_left_upper,
                d->chart1_right_lower);
    }
    fclose(fp);
    fp = as_open(dir, "boxcut_layout_disjunctions.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "disjunction,chart0,chart1,overlap_pairs,"
                "chart1_left_upper,chart1_right_lower\n");
    for (size_t i = 0; i < ndisjunction; i++) {
        const AtlasBoxcutLayoutDisjunction *d = &disjunction[i];
        fprintf(fp, "%zu,%d,%d,%zu,%.17g,%.17g\n", i, d->chart0,
                d->chart1, d->overlap_pairs, d->chart1_left_upper,
                d->chart1_right_lower);
    }
    fclose(fp);
    fp = as_open(dir, "boxcut_layout_coherence.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "relation,chart0,chart1,intrinsic_edges,seam_edges,"
                "source_weight,qp_weight,solved_delta,residual\n");
    for (size_t i = 0; i < ncoherence; i++) {
        const AtlasBoxcutLayoutCoherence *r = &coherence[i];
        fprintf(fp, "%zu,%d,%d,%zu,%zu,%.17g,%.17g,%.17g,%.17g\n",
                i, r->chart0, r->chart1, r->intrinsic_edges,
                r->seam_edges, r->source_weight, r->qp_weight,
                r->solved_delta, r->residual);
    }
    fclose(fp);
    fp = as_open(dir, "boxcut_layout_relations.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "relation,chart0,chart1,pairs,positive_pairs,negative_pairs,"
                "target_median,target_mad,sign_agreement,separation_lower,"
                "qp_weight,eligible,selected,"
                "solved_delta,residual\n");
    for (size_t i = 0; i < nrelation; i++) {
        const AtlasBoxcutLayoutRelation *r = &relation[i];
        fprintf(fp, "%zu,%d,%d,%zu,%zu,%zu,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%d,%d,%.17g,%.17g\n",
                 i, r->chart0, r->chart1, r->pairs, r->positive_pairs,
                 r->negative_pairs, r->target_median, r->target_mad,
                 r->sign_agreement, r->separation_lower, r->qp_weight,
                 r->eligible, r->selected,
                 r->solved_delta, r->residual);
    }
    fclose(fp);

    if (write_boxcut_chart_obj(dir, "atlas_boxcut_charts_flat.obj", ps, cal,
                               base_field, face_label) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_rigid_flat.obj", &split,
                               0, AS_BOXCUT_COLOR_NONE) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_rigid_flat_xyzcolor.obj",
                               &split, 0, AS_BOXCUT_COLOR_XYZ) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_rigid_flat_chartcolor.obj",
                               &split, 0, AS_BOXCUT_COLOR_CHART) != 0 ||
        write_boxcut_split_obj(dir, "atlas_boxcut_rigid_world_uv.obj", &split,
                               1, AS_BOXCUT_COLOR_NONE) != 0 ||
        export_boxcut_split_placed(arena, cfg, ps, &split, dir) != 0)
        return -1;

    fp = as_open(dir, "boxcut_rigid_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "{\n");
    fprintf(fp, "  \"operation\": \"boxcutter_multicut_global_block_qp\",\n");
    fprintf(fp, "  \"base\": \"lifted_components_plus_smooth_seams\",\n");
    fprintf(fp, "  \"block_basis\": \"%s\",\n",
            cfg->boxcut_component_blocks ? "retained_mesh_components" :
                                           "multicut_charts");
    fprintf(fp, "  \"layout_source\": \"%s\",\n",
            external_layout ? "external_chart_shifts" : "internal_qp");
    fprintf(fp, "  \"collision_order\": \"%s\",\n",
            external_layout ? "external_disjunctive_layout" :
                              "initial_parameter_chart_mean");
    fprintf(fp, "  \"slab_envelope_constraints_enabled\": false,\n");
    fprintf(fp, "  \"faces_before\": %zu, \"faces_after\": %zu,\n",
            ps->nf, split.nf);
    fprintf(fp, "  \"vertices_before\": %zu, \"vertices_after_chart_split\": %zu,\n",
            ps->nv, split.nv);
    fprintf(fp, "  \"pitch\": %.17g,\n", cal->pitch);
    fprintf(fp, "  \"partition_radius_gate\": %.17g, \"layout_radius_gate\": %.17g,\n",
            partition_radius_gate, layout_radius_gate);
    fprintf(fp, "  \"exact_overlap_pairs_before\": %zu,\n",
            overlap_before.exact_face_pairs);
    fprintf(fp, "  \"within_layout_block_overlap_pairs_before\": %zu, \"cross_layout_block_overlap_pairs_before\": %zu,\n",
            within_layout_block_overlap_pairs_before,
            cross_layout_block_overlap_pairs_before);
    fprintf(fp, "  \"radial_repulsive_edges\": %zu, \"repulsive_edges_cut\": %zu, \"repulsive_edges_joined\": %zu,\n",
            box_stats.radial_repulsive_edges, box_stats.repulsive_edges_cut,
            box_stats.repulsive_edges_joined);
    fprintf(fp, "  \"near_radius_pairs_ignored_by_partition\": %zu,\n",
            box_stats.near_radius_pairs_ignored);
    fprintf(fp, "  \"allowed_delamination_pairs_before\": %zu, \"partition_allowed_delamination_pairs\": %zu,\n",
            nallowed_delamination_pair,
            box_stats.allowed_delamination_pairs);
    fprintf(fp, "  \"charts\": %zu, \"multicut_charts\": %zu, \"layout_blocks\": %zu,\n",
            layout_blocks, multicut_charts, layout_blocks);
    fprintf(fp, "  \"eligible_seam_bundles\": %zu, \"verified_seam_edges_across_blocks\": %zu,\n",
            neligible_seam_bundle, seam_edges_cut);
    fprintf(fp, "  \"verified_seam_edges_cut\": %zu,\n", seam_edges_cut);
    fprintf(fp, "  \"layout_observations\": %zu, \"layout_relations\": %zu, \"eligible_relations\": %zu, \"qp_relation_rows\": %zu,\n",
             layout_stats.observations, layout_stats.relations,
             layout_stats.eligible_relations, layout_stats.selected_relations);
    fprintf(fp, "  \"coherence_weight_scale\": %.17g, \"coherence_relations\": %zu,\n",
            layout_problem.coherence_weight_scale,
            layout_stats.coherence_relations);
    fprintf(fp, "  \"coherence_edges\": %zu, \"coherence_intrinsic_edges\": %zu, \"coherence_seam_edges\": %zu, \"coherence_direct_overlap_edges_excluded\": %zu,\n",
            layout_stats.coherence_adjacency_edges,
            layout_stats.coherence_intrinsic_edges,
            layout_stats.coherence_seam_edges,
            layout_stats.coherence_direct_overlap_edges_excluded);
    fprintf(fp, "  \"mixed_direction_relations\": %zu, \"radial_direction_relations_rejected\": %zu,\n",
             layout_stats.mixed_direction_relations,
             layout_stats.collision_relations_rejected);
    fprintf(fp, "  \"collision_bounds\": %zu, \"qp_iterations\": %d, \"qp_active_bounds\": %d,\n",
            layout_stats.collision_bounds, layout_stats.qp_iterations,
            layout_stats.qp_active_bounds);
    fprintf(fp, "  \"collision_disjunctions\": %zu,\n",
            layout_stats.collision_disjunctions);
    fprintf(fp, "  \"remaining_collision_disjunctions\": %zu,\n",
            nremaining_disjunction);
    fprintf(fp, "  \"collision_outer_iterations\": %zu, \"collision_bounds_added\": %zu,\n",
            layout_stats.collision_outer_iterations,
            layout_stats.collision_bounds_added);
    fprintf(fp, "  \"slab_height\": %.17g, \"slab_records\": %zu, \"slabs\": %zu, \"slab_bounds\": %zu,\n",
            layout_stats.slab_height, layout_stats.slab_records,
            layout_stats.slab_count, layout_stats.slab_bounds);
    fprintf(fp, "  \"layout_initial_exact_overlap_pairs\": %zu, \"layout_initial_cross_chart_overlap_pairs\": %zu,\n",
            layout_stats.initial_exact_overlap_pairs,
            layout_stats.initial_cross_overlap_pairs);
    fprintf(fp, "  \"layout_final_exact_overlap_pairs\": %zu, \"layout_final_cross_chart_overlap_pairs\": %zu,\n",
            layout_stats.final_exact_overlap_pairs,
            layout_stats.final_cross_overlap_pairs);
    fprintf(fp, "  \"qp_objective\": %.17g, \"qp_stationarity\": %.17g, \"qp_max_violation\": %.17g,\n",
            layout_stats.qp_objective, layout_stats.qp_stationarity,
            layout_stats.qp_max_violation);
    fprintf(fp, "  \"radial_order_pair_fraction\": %.17g, \"target_pair_fraction\": %.17g,\n",
            layout_stats.radial_order_pair_fraction,
            layout_stats.target_pair_fraction);
    fprintf(fp, "  \"moved_charts\": %zu, \"moved_faces\": %zu, \"maximum_abs_shift\": %.17g,\n",
            layout_stats.moved_charts, layout_stats.moved_faces,
            layout_stats.maximum_abs_shift);
    fprintf(fp, "  \"relation_residual_rms\": %.17g, \"relation_residual_max\": %.17g,\n",
            layout_stats.relation_residual_rms,
            layout_stats.relation_residual_max);
    fprintf(fp, "  \"coherence_relations_preserved\": %zu, \"coherence_edges_preserved\": %zu,\n",
            layout_stats.coherence_relations_preserved,
            layout_stats.coherence_edges_preserved);
    fprintf(fp, "  \"coherence_residual_rms\": %.17g, \"coherence_weighted_residual_rms\": %.17g,\n",
            layout_stats.coherence_residual_rms,
            layout_stats.coherence_weighted_residual_rms);
    fprintf(fp, "  \"coherence_residual_median\": %.17g, \"coherence_residual_p95\": %.17g, \"coherence_residual_max\": %.17g,\n",
            layout_stats.coherence_residual_median,
            layout_stats.coherence_residual_p95,
            layout_stats.coherence_residual_max);
    fprintf(fp, "  \"relative_flips\": %zu, \"comparable_faces\": %zu, \"degenerate_faces\": %zu,\n",
            relative_flips, comparable, degenerate);
    fprintf(fp, "  \"exact_overlap_pairs_after\": %zu, \"within_layout_block_overlap_pairs_after\": %zu, \"cross_layout_block_overlap_pairs_after\": %zu,\n",
            overlap_after.exact_face_pairs, overlap_after.same_component_pairs,
            overlap_after.cross_component_pairs);
    fprintf(fp, "  \"cross_chart_overlap_pairs_after\": %zu,\n",
            overlap_after.cross_component_pairs);
    fprintf(fp, "  \"allowed_delamination_pairs_after\": %zu, \"unwhitelisted_overlap_pairs_after\": %zu,\n",
            final_allowed_delamination_pairs,
            overlap_after.exact_face_pairs - final_allowed_delamination_pairs);
    fprintf(fp, "  \"xyz_bpa_weld_oracle_rc\": %d, \"xyz_bridge_faces\": %zu,\n",
            xyz_weld_oracle_rc, xyz_weld_topology.nfaces);
    fprintf(fp, "  \"xyz_weld_final_edge_stretch_p95\": %.17g, \"xyz_weld_final_face_stretch_p95\": %.17g,\n",
            xyz_weld_final.best_cross_cube_edge_stretch_p95,
            xyz_weld_final.best_bridge_face_symmetric_stretch_p95);
    fprintf(fp, "  \"xyz_weld_final_split_vertices\": %zu, \"xyz_weld_final_discontinuous_vertices\": %zu,\n",
            xyz_weld_final.split_bridge_vertices,
            xyz_weld_final.discontinuous_embedded_bridge_vertices);
    fprintf(fp, "  \"xyz_weld_final_introduced_overlap_pairs\": %zu,\n",
            xyz_weld_final.bridge_base_overlap_pairs +
                xyz_weld_final.bridge_bridge_overlap_pairs);
    fprintf(fp, "  \"elapsed_seconds\": %.17g\n}\n",
            ves_clock_sec() - t0);
    fclose(fp);
    fprintf(stderr,
        "[atlas_strip_scroll] BoxCutter rigid %.2fs: charts=%zu split_v=%zu "
        "repulsive=%zu/%zu layout=%zu/%zu radial=%.1f%% "
        "coherence=%zu/%zu p95=%.3g overlap=%zu->%zu "
        "(allowed=%zu) flips=%zu "
        "seam_cuts=%zu\n",
        ves_clock_sec() - t0, layout_blocks, split.nv,
        box_stats.repulsive_edges_cut, box_stats.radial_repulsive_edges,
        layout_stats.selected_relations, layout_stats.relations,
        100.0 * layout_stats.radial_order_pair_fraction,
        layout_stats.coherence_relations_preserved,
        layout_stats.coherence_relations,
        layout_stats.coherence_residual_p95,
        overlap_before.exact_face_pairs, overlap_after.exact_face_pairs,
        final_allowed_delamination_pairs,
        relative_flips, seam_edges_cut);
    if (cfg->boxcut_layout_only)
        return relative_flips == 0 ? 0 : -1;
    int rankpack_rc = run_boxcut_rank_pack(
        arena, cfg, ps, axial, base_field, face_radius, face_label,
        layout_blocks, box_stats.average_uv_edge_length);
    int winding_rc = run_boxcut_winding_stage(
        arena, cfg, ps, cal, axial, first_field, face_radius, face_label,
        layout_blocks, &box_stats);
    return relative_flips == 0 && rankpack_rc == 0 && winding_rc == 0
        ? 0 : -1;
}

static int run_boxcut_radial_stage(Arena_T arena, const ScrollConfig *cfg,
                                   const PieceSet *ps,
                                   const ScaffoldCalib *cal,
                                   const AtlasCandidateSet *set,
                                   const double *base_field)
{
    if (base_field == NULL) return -1;
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "boxcut_radial") != 0) return -1;

    double *axial = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    float *registered = (float *)ARENA_ALLOC(
        arena, ps->nv * sizeof(float));
    for (size_t i = 0; i < ps->nv; i++) {
        axial[i] = axis_coordinate(cal, &ps->verts[i * 3]);
        registered[i] = ps->uv[i * 2];
    }

    AtlasOverlapAudit overlap_before;
    if (AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, base_field, axial,
            registered, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &overlap_before) != 0)
        return -1;
    AtlasSeamAudit seam;
    if (AtlasSeamAudit_build(
            arena, ps, base_field, set->vertex_mesh_component,
            set->mesh_components, &seam) != 0)
        return -1;
    uint8_t *eligible_seam = (uint8_t *)ARENA_CALLOC(
        arena, (seam.nbundles ? seam.nbundles : 1), sizeof(uint8_t));
    size_t neligible_seam_bundle = 0, neligible_seam_pair = 0;
    for (size_t i = 0; i < seam.nbundles; i++) {
        eligible_seam[i] =
            (uint8_t)smooth_seam_bundle_eligible(&seam.bundles[i]);
        if (eligible_seam[i]) neligible_seam_bundle++;
    }
    for (size_t i = 0; i < seam.npairs; i++)
        if (eligible_seam[seam.pairs[i].bundle]) neligible_seam_pair++;

    double radius_gate = cal->pitch > 0.0 ? 0.5 * cal->pitch : 4.75;
    AtlasBoxcutProblem box_problem;
    memset(&box_problem, 0, sizeof box_problem);
    box_problem.vertices = ps->verts;
    box_problem.nvertices = ps->nv;
    box_problem.faces = ps->faces;
    box_problem.nfaces = ps->nf;
    box_problem.u = base_field;
    box_problem.v = axial;
    box_problem.axis_point = cal->axis_point;
    box_problem.axis_dir = cal->axis_dir;
    box_problem.overlap = &overlap_before;
    box_problem.seam = &seam;
    box_problem.eligible_seam_bundle = eligible_seam;
    box_problem.radius_gate = radius_gate;
    int32_t *face_label = NULL;
    AtlasBoxcutStats box_stats;
    double t0 = ves_clock_sec();
    int box_rc = AtlasBoxcut_partition(
        arena, &box_problem, &face_label, &box_stats);
    fprintf(stderr,
        "[atlas_strip_scroll] BoxCutter %.2fs: rc=%d charts=%zu "
        "adj=%zu+%zu overlap=%zu radial=%zu cut=%zu joined=%zu\n",
        ves_clock_sec() - t0, box_rc, box_stats.charts,
        box_stats.intrinsic_adjacency_edges,
        box_stats.seam_adjacency_edges, box_stats.exact_overlap_pairs,
        box_stats.radial_repulsive_edges, box_stats.repulsive_edges_cut,
        box_stats.repulsive_edges_joined);
    if (box_rc != 0 || box_stats.charts == 0 || face_label == NULL ||
        box_stats.radial_repulsive_edges == 0)
        return -1;

    double *face_radius = (double *)ARENA_ALLOC(
        arena, ps->nf * sizeof(double));
    size_t *label_faces = (size_t *)ARENA_CALLOC(
        arena, box_stats.charts, sizeof(size_t));
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_label[f] < 0 ||
            (size_t)face_label[f] >= box_stats.charts)
            return -1;
        face_radius[f] = boxcut_face_radius(ps, cal, f);
        label_faces[face_label[f]]++;
    }

    /* BoxCutter charts can be long spiral segments, so a single chart-wide
     * translation has no meaningful radial order.  Seed the actual collision
     * faces, grow two same-chart topology rings, and use the connected grown
     * regions as local StrokeStrip support patches. */
    const int support_rings = 2;
    uint8_t *support = (uint8_t *)ARENA_CALLOC(
        arena, ps->nf, sizeof(uint8_t));
    size_t *label_collision_endpoints = (size_t *)ARENA_CALLOC(
        arena, box_stats.charts, sizeof(size_t));
    size_t ncollision_endpoint = 0;
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        const AtlasOverlapPair *p = &overlap_before.pairs[i];
        if (fabs(face_radius[p->face1] - face_radius[p->face0]) <=
            radius_gate)
            continue;
        int32_t l0 = face_label[p->face0], l1 = face_label[p->face1];
        if (l0 == l1) continue;
        support[p->face0] = support[p->face1] = 1;
        label_collision_endpoints[l0]++;
        label_collision_endpoints[l1]++;
        ncollision_endpoint += 2;
    }
    for (int ring = 0; ring < support_rings; ring++) {
        uint8_t *grown = (uint8_t *)ARENA_ALLOC(
            arena, ps->nf * sizeof(uint8_t));
        memcpy(grown, support, ps->nf * sizeof(uint8_t));
        for (size_t i = 0; i < box_stats.original_graph_edges; i++) {
            int32_t a = box_stats.adjacency_face0[i];
            int32_t b = box_stats.adjacency_face1[i];
            if (face_label[a] != face_label[b]) continue;
            if (support[a] || support[b]) grown[a] = grown[b] = 1;
        }
        support = grown;
    }
    UnionFind support_uf = UF_new(arena, (int32_t)ps->nf);
    for (size_t i = 0; i < box_stats.original_graph_edges; i++) {
        int32_t a = box_stats.adjacency_face0[i];
        int32_t b = box_stats.adjacency_face1[i];
        if (support[a] && support[b] && face_label[a] == face_label[b])
            uf_union(&support_uf, a, b);
    }
    int32_t *root_patch = (int32_t *)ARENA_ALLOC(
        arena, ps->nf * sizeof(int32_t));
    int32_t *face_patch = (int32_t *)ARENA_ALLOC(
        arena, ps->nf * sizeof(int32_t));
    for (size_t f = 0; f < ps->nf; f++) {
        root_patch[f] = -1;
        face_patch[f] = -1;
    }
    size_t npatch = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        if (!support[f]) continue;
        int32_t root = uf_find(&support_uf, (int32_t)f);
        if (root_patch[root] < 0) root_patch[root] = (int32_t)npatch++;
        face_patch[f] = root_patch[root];
    }
    if (npatch == 0) return -1;
    size_t *patch_faces = (size_t *)ARENA_CALLOC(
        arena, npatch, sizeof(size_t));
    size_t *patch_collision_endpoints = (size_t *)ARENA_CALLOC(
        arena, npatch, sizeof(size_t));
    int32_t *patch_chart = (int32_t *)ARENA_ALLOC(
        arena, npatch * sizeof(int32_t));
    for (size_t i = 0; i < npatch; i++) patch_chart[i] = -1;
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_patch[f] < 0) continue;
        int32_t patch = face_patch[f];
        patch_faces[patch]++;
        if (patch_chart[patch] < 0) patch_chart[patch] = face_label[f];
        else if (patch_chart[patch] != face_label[f]) return -1;
    }

    size_t *radius_offset = (size_t *)ARENA_CALLOC(
        arena, npatch + 1, sizeof(size_t));
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        const AtlasOverlapPair *p = &overlap_before.pairs[i];
        if (fabs(face_radius[p->face1] - face_radius[p->face0]) <=
            radius_gate || face_label[p->face0] == face_label[p->face1])
            continue;
        int32_t p0 = face_patch[p->face0], p1 = face_patch[p->face1];
        if (p0 < 0 || p1 < 0 || p0 == p1) return -1;
        radius_offset[(size_t)p0 + 1]++;
        radius_offset[(size_t)p1 + 1]++;
        patch_collision_endpoints[p0]++;
        patch_collision_endpoints[p1]++;
    }
    for (size_t i = 0; i < npatch; i++)
        radius_offset[i + 1] += radius_offset[i];
    size_t *radius_cursor = (size_t *)ARENA_ALLOC(
        arena, npatch * sizeof(size_t));
    memcpy(radius_cursor, radius_offset, npatch * sizeof(size_t));
    double *radius_value = (double *)ARENA_ALLOC(
        arena, (ncollision_endpoint ? ncollision_endpoint : 1) *
                   sizeof(double));
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        const AtlasOverlapPair *p = &overlap_before.pairs[i];
        if (fabs(face_radius[p->face1] - face_radius[p->face0]) <=
            radius_gate || face_label[p->face0] == face_label[p->face1])
            continue;
        int32_t p0 = face_patch[p->face0], p1 = face_patch[p->face1];
        radius_value[radius_cursor[p0]++] = face_radius[p->face0];
        radius_value[radius_cursor[p1]++] = face_radius[p->face1];
    }
    double *patch_radius = (double *)ARENA_ALLOC(
        arena, npatch * sizeof(double));
    for (size_t patch = 0; patch < npatch; patch++) {
        size_t count = radius_offset[patch + 1] - radius_offset[patch];
        patch_radius[patch] = count
            ? lift_median(&radius_value[radius_offset[patch]], count)
            : NAN;
    }

    AsBoxcutOrderObservation *observation =
        (AsBoxcutOrderObservation *)ARENA_ALLOC(
            arena, box_stats.radial_repulsive_edges *
                       sizeof(AsBoxcutOrderObservation));
    size_t nobservation = 0;
    const double two_pi = 6.283185307179586476925286766559;
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        const AtlasOverlapPair *p = &overlap_before.pairs[i];
        double r0 = face_radius[p->face0], r1 = face_radius[p->face1];
        if (fabs(r1 - r0) <= radius_gate) continue;
        int32_t l0 = face_label[p->face0], l1 = face_label[p->face1];
        if (l0 == l1) continue;
        int32_t p0 = face_patch[p->face0], p1 = face_patch[p->face1];
        if (p0 < 0 || p1 < 0 || p0 == p1) return -1;
        int p0_inner = patch_radius[p0] < patch_radius[p1];
        if (patch_radius[p0] == patch_radius[p1]) p0_inner = r0 < r1;
        int32_t inner = p0_inner ? p0 : p1;
        int32_t outer = p0_inner ? p1 : p0;
        size_t inner_face = (size_t)(p0_inner ? p->face0 : p->face1);
        size_t outer_face = (size_t)(p0_inner ? p->face1 : p->face0);
        double base_du = boxcut_face_mean(ps, outer_face, base_field) -
                         boxcut_face_mean(ps, inner_face, base_field);
        double desired_du = fabs(
            boxcut_face_registered_u(ps, outer_face) -
            boxcut_face_registered_u(ps, inner_face));
        double predicted_du = fabs((double)p->turn) * two_pi *
                              0.5 * (r0 + r1);
        if (desired_du < 4.0 * box_stats.average_uv_edge_length &&
            predicted_du > desired_du)
            desired_du = predicted_du;
        double correction = desired_du - base_du;
        double clearance = 2.0 * box_stats.average_uv_edge_length - base_du;
        if (correction < clearance) correction = clearance;
        if (!isfinite(correction) || correction <= 0.0) continue;
        AsBoxcutOrderObservation *o = &observation[nobservation++];
        o->inner_label = inner;
        o->outer_label = outer;
        o->correction_target = correction;
        o->desired_du = desired_du;
        o->base_du = base_du;
        o->local_order_consistent = p0_inner ? r0 < r1 : r1 < r0;
    }
    if (nobservation == 0) return -1;
    qsort(observation, nobservation, sizeof(AsBoxcutOrderObservation),
          compare_boxcut_order_observation);
    AsBoxcutOrderBundle *order = (AsBoxcutOrderBundle *)ARENA_ALLOC(
        arena, nobservation * sizeof(AsBoxcutOrderBundle));
    double *scratch = (double *)ARENA_ALLOC(
        arena, nobservation * sizeof(double));
    size_t norder = 0;
    for (size_t first = 0; first < nobservation;) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].inner_label ==
                   observation[first].inner_label &&
               observation[last].outer_label ==
                   observation[first].outer_label)
            last++;
        AsBoxcutOrderBundle *b = &order[norder++];
        b->inner_label = observation[first].inner_label;
        b->outer_label = observation[first].outer_label;
        b->pairs = last - first;
        size_t count = 0, consistent = 0;
        for (size_t i = first; i < last; i++) {
            scratch[count++] = observation[i].correction_target;
            consistent += (size_t)observation[i].local_order_consistent;
        }
        b->correction_target_median = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = fabs(observation[i].correction_target -
                                      b->correction_target_median);
        b->correction_target_mad = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = observation[i].desired_du;
        b->desired_du_median = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = observation[i].base_du;
        b->base_du_median = lift_median(scratch, count);
        b->local_radial_agreement = (double)consistent / (double)b->pairs;
        first = last;
    }

    uint8_t *patch_active = (uint8_t *)ARENA_CALLOC(
        arena, npatch, sizeof(uint8_t));
    for (size_t i = 0; i < norder; i++) {
        patch_active[order[i].inner_label] = 1;
        patch_active[order[i].outer_label] = 1;
    }
    int32_t *patch_aux = (int32_t *)ARENA_ALLOC(
        arena, npatch * sizeof(int32_t));
    size_t nauxiliary = 0;
    for (size_t i = 0; i < npatch; i++)
        patch_aux[i] = patch_active[i] ? (int32_t)nauxiliary++ : -1;
    if (ps->nv > (size_t)INT32_MAX - nauxiliary) return -1;
    size_t nsupport = 0;
    for (size_t f = 0; f < ps->nf; f++)
        if (face_patch[f] >= 0 && patch_active[face_patch[f]]) nsupport++;
    size_t nconstraint = neligible_seam_pair + nsupport + norder;
    size_t ncoeff = 2 * neligible_seam_pair + 4 * nsupport + 2 * norder;
    AtlasFieldConstraint *constraint =
        (AtlasFieldConstraint *)ARENA_ALLOC(
            arena, nconstraint * sizeof(AtlasFieldConstraint));
    AtlasFieldConstraintCoeff *coeff =
        (AtlasFieldConstraintCoeff *)ARENA_ALLOC(
            arena, ncoeff * sizeof(AtlasFieldConstraintCoeff));
    size_t nr = 0, nc = 0;
    const double seam_weight = 1000.0;
    const double support_weight = 1.0;
    const double order_weight = 1000.0;
    for (size_t i = 0; i < seam.npairs; i++) {
        const AtlasSeamPair *p = &seam.pairs[i];
        if (!eligible_seam[p->bundle]) continue;
        AtlasFieldConstraint *row = &constraint[nr++];
        row->first = nc; row->count = 2; row->target = 0.0;
        row->weight = seam_weight; row->kind = ATLAS_FIELD_ROW_SEAM;
        row->source = i <= (size_t)INT32_MAX ? (int32_t)i : INT32_MAX;
        coeff[nc].variable = p->vertex0; coeff[nc++].coefficient = -1.0;
        coeff[nc].variable = p->vertex1; coeff[nc++].coefficient = 1.0;
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t patch = face_patch[f];
        if (patch < 0 || !patch_active[patch]) continue;
        AtlasFieldConstraint *row = &constraint[nr++];
        row->first = nc; row->count = 4;
        row->target = boxcut_face_mean(ps, f, base_field);
        row->weight = support_weight; row->kind = ATLAS_FIELD_ROW_SUPPORT;
        row->source = f <= (size_t)INT32_MAX ? (int32_t)f : INT32_MAX;
        for (int k = 0; k < 3; k++) {
            coeff[nc].variable = ps->faces[f * 3 + (size_t)k];
            coeff[nc++].coefficient = 1.0 / 3.0;
        }
        coeff[nc].variable = (int32_t)ps->nv + patch_aux[patch];
        coeff[nc++].coefficient = -1.0;
    }
    for (size_t i = 0; i < norder; i++) {
        AtlasFieldConstraint *row = &constraint[nr++];
        row->first = nc; row->count = 2;
        row->target = order[i].correction_target_median;
        row->weight = order_weight; row->kind = ATLAS_FIELD_ROW_ORDER;
        row->source = i <= (size_t)INT32_MAX ? (int32_t)i : INT32_MAX;
        coeff[nc].variable = (int32_t)ps->nv +
            patch_aux[order[i].inner_label];
        coeff[nc++].coefficient = -1.0;
        coeff[nc].variable = (int32_t)ps->nv +
            patch_aux[order[i].outer_label];
        coeff[nc++].coefficient = 1.0;
    }
    if (nr != nconstraint || nc != ncoeff) return -1;

    double *aux_initial = (double *)ARENA_CALLOC(
        arena, nauxiliary, sizeof(double));
    double *aux_final = (double *)ARENA_ALLOC(
        arena, nauxiliary * sizeof(double));
    double *final_field = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    AtlasFieldRefineStats refine;
    t0 = ves_clock_sec();
    int solve_rc = AtlasFieldRefine_solve(
        arena, ps->verts, ps->nv, ps->faces, ps->nf,
        set->vertex_mesh_component, set->mesh_components, base_field,
        constraint, nconstraint, coeff, ncoeff,
        nauxiliary, aux_initial, &cfg->qp,
        final_field, aux_final, &refine);
    fprintf(stderr,
        "[atlas_strip_scroll] BoxCutter radial solve %.2fs: rc=%d "
        "support_patches=%zu active=%zu faces=%zu order=%zu "
        "rms=%.4g->%.4g "
        "seam=%.4g->%.4g correction=%.4g grad=%.4g\n",
        ves_clock_sec() - t0, solve_rc, npatch, nauxiliary, nsupport, norder,
        refine.order_rms_before, refine.order_rms_after,
        refine.seam_rms_before, refine.seam_rms_after,
        refine.correction_rms, refine.edge_gradient_delta_rms);
    if (solve_rc != 0) return -1;

    AtlasOverlapAudit overlap_after;
    if (AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, final_field, axial,
            registered, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &overlap_after) != 0)
        return -1;
    size_t comparable = 0, relative_flips = 0, degenerate = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = ps->faces[f * 3];
        int32_t b = ps->faces[f * 3 + 1];
        int32_t c = ps->faces[f * 3 + 2];
        double d0 = (base_field[b] - base_field[a]) *
                    (axial[c] - axial[a]) -
                    (base_field[c] - base_field[a]) *
                    (axial[b] - axial[a]);
        double d1 = (final_field[b] - final_field[a]) *
                    (axial[c] - axial[a]) -
                    (final_field[c] - final_field[a]) *
                    (axial[b] - axial[a]);
        if (fabs(d1) < 1e-10) degenerate++;
        if (fabs(d0) >= 1e-10) {
            comparable++;
            if (d0 * d1 < 0.0) relative_flips++;
        }
    }

    FILE *fp = as_open(dir, "boxcut_face_labels.i32", "wb");
    if (fp == NULL ||
        fwrite(face_label, sizeof(int32_t), ps->nf, fp) != ps->nf ||
        fclose(fp) != 0)
        return -1;
    fp = as_open(dir, "boxcut_support_patch.i32", "wb");
    if (fp == NULL ||
        fwrite(face_patch, sizeof(int32_t), ps->nf, fp) != ps->nf ||
        fclose(fp) != 0)
        return -1;
    fp = as_open(dir, "boxcut_charts.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "chart,faces,collision_endpoints\n");
    for (size_t i = 0; i < box_stats.charts; i++)
        fprintf(fp, "%zu,%zu,%zu\n", i, label_faces[i],
                label_collision_endpoints[i]);
    fclose(fp);
    fp = as_open(dir, "boxcut_support_patches.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "patch,chart,faces,collision_endpoints,median_radius,"
                "active,auxiliary,solved_shift\n");
    for (size_t i = 0; i < npatch; i++)
        fprintf(fp, "%zu,%d,%zu,%zu,%.17g,%d,%d,%.17g\n",
                i, patch_chart[i], patch_faces[i],
                patch_collision_endpoints[i], patch_radius[i],
                patch_active[i], patch_aux[i],
                patch_active[i] ? aux_final[patch_aux[i]] : 0.0);
    fclose(fp);
    fp = as_open(dir, "boxcut_radial_order.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "relation,inner_patch,outer_patch,inner_chart,outer_chart,"
                "pairs,inner_radius,outer_radius,local_radial_agreement,"
                "base_du,desired_du,correction_target,correction_mad,"
                "solved_correction,residual\n");
    for (size_t i = 0; i < norder; i++) {
        double solved = aux_final[patch_aux[order[i].outer_label]] -
                        aux_final[patch_aux[order[i].inner_label]];
        fprintf(fp, "%zu,%d,%d,%d,%d,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%.17g\n",
                i, order[i].inner_label, order[i].outer_label,
                patch_chart[order[i].inner_label],
                patch_chart[order[i].outer_label], order[i].pairs,
                patch_radius[order[i].inner_label],
                patch_radius[order[i].outer_label],
                order[i].local_radial_agreement,
                order[i].base_du_median, order[i].desired_du_median,
                order[i].correction_target_median,
                order[i].correction_target_mad, solved,
                solved - order[i].correction_target_median);
    }
    fclose(fp);

    ScrollConfig final_cfg = *cfg;
    final_cfg.out_dir = dir;
    if (write_global_mesh_obj(dir, "atlas_smooth_seam_flat.obj", ps, cal,
                              base_field, 1) != 0 ||
        write_boxcut_chart_obj(dir, "atlas_boxcut_charts_flat.obj", ps, cal,
                               base_field, face_label) != 0 ||
        write_global_mesh_obj(dir, "atlas_boxcut_radial_flat.obj", ps, cal,
                              final_field, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_boxcut_radial_world_uv.obj", ps,
                              cal, final_field, 0) != 0 ||
        write_global_mesh_xyz_obj(dir,
                                  "atlas_boxcut_radial_flat_xyzcolor.obj",
                                  ps, cal, final_field) != 0 ||
        export_field_placed(&final_cfg, ps, cal, final_field) != 0)
        return -1;

    fp = as_open(dir, "boxcut_radial_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n  \"operation\": \"boxcutter_multicut_plus_radial_fem\",\n"
        "  \"base\": \"smooth_seams\",\n"
        "  \"faces_before\": %zu,\n  \"faces_after\": %zu,\n"
        "  \"pitch\": %.17g,\n  \"radius_gate\": %.17g,\n"
        "  \"intrinsic_adjacency_edges\": %zu,\n"
        "  \"seam_adjacency_edges\": %zu,\n"
        "  \"eligible_seam_bundles\": %zu,\n"
        "  \"eligible_seam_pairs\": %zu,\n"
        "  \"exact_overlap_pairs_before\": %zu,\n"
        "  \"radial_repulsive_edges\": %zu,\n"
        "  \"near_radius_pairs_ignored\": %zu,\n"
        "  \"repulsive_edges_cut\": %zu,\n"
        "  \"repulsive_edges_joined\": %zu,\n"
        "  \"boxcut_charts\": %zu,\n"
        "  \"support_growth_rings\": %d,\n"
        "  \"support_patches\": %zu,\n"
        "  \"active_support_patches\": %zu,\n"
        "  \"support_rows\": %zu,\n"
        "  \"radial_order_rows\": %zu,\n"
        "  \"order_rms_before\": %.17g,\n"
        "  \"order_rms_after\": %.17g,\n"
        "  \"order_max_after\": %.17g,\n"
        "  \"seam_rms_before\": %.17g,\n"
        "  \"seam_rms_after\": %.17g,\n"
        "  \"correction_rms\": %.17g,\n"
        "  \"correction_max\": %.17g,\n"
        "  \"edge_gradient_delta_rms\": %.17g,\n"
        "  \"edge_gradient_delta_max\": %.17g,\n"
        "  \"relative_flips\": %zu,\n"
        "  \"comparable_faces\": %zu,\n"
        "  \"degenerate_faces\": %zu,\n"
        "  \"exact_overlap_pairs_after\": %zu,\n"
        "  \"cross_component_overlap_pairs_after\": %zu,\n"
        "  \"qp_objective_initial\": %.17g,\n"
        "  \"qp_objective_final\": %.17g,\n"
        "  \"qp_linear_residual\": %.17g\n}\n",
        ps->nf, ps->nf, cal->pitch, radius_gate,
        box_stats.intrinsic_adjacency_edges,
        box_stats.seam_adjacency_edges, neligible_seam_bundle,
        neligible_seam_pair, overlap_before.exact_face_pairs,
        box_stats.radial_repulsive_edges,
        box_stats.near_radius_pairs_ignored,
        box_stats.repulsive_edges_cut, box_stats.repulsive_edges_joined,
        box_stats.charts, support_rings, npatch, nauxiliary, nsupport, norder,
        refine.order_rms_before, refine.order_rms_after,
        refine.order_max_after, refine.seam_rms_before,
        refine.seam_rms_after, refine.correction_rms,
        refine.correction_max, refine.edge_gradient_delta_rms,
        refine.edge_gradient_delta_max, relative_flips, comparable,
        degenerate, overlap_after.exact_face_pairs,
        overlap_after.cross_component_pairs,
        refine.qp.objective_initial, refine.qp.objective_final,
        refine.qp.max_reduced_linear_residual);
    fclose(fp);
    return 0;
}

typedef struct {
    double value;
    double weight;
} AsLiftWeightedValue;

static int compare_lift_weighted_value(const void *pa, const void *pb)
{
    const AsLiftWeightedValue *a = (const AsLiftWeightedValue *)pa;
    const AsLiftWeightedValue *b = (const AsLiftWeightedValue *)pb;
    return a->value < b->value ? -1 : (a->value > b->value ? 1 : 0);
}

/* Build the same protected quotient used by AtlasComponentLift_solve, then
 * assign one robust registered-minus-HEAD gauge prior to every whole group.
 * A protected chart must have one translation even when its constituent mesh
 * components have noisy or contradictory raw absolute priors. */
static int build_lift_protected_groups(
    Arena_T arena, size_t ncomponents,
    const AtlasComponentLiftEdge *protected_edge, size_t nprotected,
    const double *component_weight, const double *component_prior,
    int32_t **out_component_group, size_t *out_ngroup,
    double **out_group_prior, double **out_component_group_prior)
{
    if (ncomponents == 0 || ncomponents > (size_t)INT32_MAX ||
        component_weight == NULL || component_prior == NULL ||
        out_component_group == NULL || out_ngroup == NULL ||
        out_group_prior == NULL || out_component_group_prior == NULL)
        return -1;
    UnionFind uf = UF_new(arena, (int32_t)ncomponents);
    for (size_t i = 0; i < nprotected; i++) {
        int32_t a = protected_edge[i].component0;
        int32_t b = protected_edge[i].component1;
        if (a < 0 || b < 0 || (size_t)a >= ncomponents ||
            (size_t)b >= ncomponents)
            return -1;
        uf_union(&uf, a, b);
    }
    int32_t *root_group = (int32_t *)ARENA_ALLOC(
        arena, ncomponents * sizeof(int32_t));
    int32_t *component_group = (int32_t *)ARENA_ALLOC(
        arena, ncomponents * sizeof(int32_t));
    for (size_t i = 0; i < ncomponents; i++) root_group[i] = -1;
    size_t ngroup = 0;
    for (size_t i = 0; i < ncomponents; i++) {
        int32_t root = uf_find(&uf, (int32_t)i);
        if (root_group[root] < 0) root_group[root] = (int32_t)ngroup++;
        component_group[i] = root_group[root];
    }
    double *group_prior = (double *)ARENA_ALLOC(
        arena, ngroup * sizeof(double));
    double *component_group_prior = (double *)ARENA_ALLOC(
        arena, ncomponents * sizeof(double));
    AsLiftWeightedValue *scratch = (AsLiftWeightedValue *)ARENA_ALLOC(
        arena, ncomponents * sizeof(AsLiftWeightedValue));
    for (size_t g = 0; g < ngroup; g++) {
        size_t count = 0;
        double total_weight = 0.0;
        for (size_t c = 0; c < ncomponents; c++) {
            if (component_group[c] != (int32_t)g) continue;
            if (!isfinite(component_prior[c]) ||
                !isfinite(component_weight[c]) || component_weight[c] < 0.0)
                return -1;
            scratch[count].value = component_prior[c];
            scratch[count].weight = component_weight[c];
            total_weight += component_weight[c];
            count++;
        }
        if (count == 0) return -1;
        if (total_weight <= 0.0) {
            total_weight = (double)count;
            for (size_t i = 0; i < count; i++) scratch[i].weight = 1.0;
        }
        qsort(scratch, count, sizeof(AsLiftWeightedValue),
              compare_lift_weighted_value);
        double cumulative = 0.0;
        group_prior[g] = scratch[count - 1].value;
        for (size_t i = 0; i < count; i++) {
            cumulative += scratch[i].weight;
            if (2.0 * cumulative >= total_weight) {
                group_prior[g] = scratch[i].value;
                break;
            }
        }
    }
    for (size_t c = 0; c < ncomponents; c++)
        component_group_prior[c] = group_prior[component_group[c]];
    *out_component_group = component_group;
    *out_ngroup = ngroup;
    *out_group_prior = group_prior;
    *out_component_group_prior = component_group_prior;
    return 0;
}

typedef enum {
    AS_OVERLAP_SELECTED = 0,
    AS_OVERLAP_SAME_COMPONENT,
    AS_OVERLAP_SAME_PROTECTED_GROUP,
    AS_OVERLAP_WEAK_PAIR_SUPPORT,
    AS_OVERLAP_WEAK_FACE_SUPPORT,
    AS_OVERLAP_ZERO_TURN,
    AS_OVERLAP_TURN_DISAGREEMENT,
    AS_OVERLAP_PHASE_RESIDUAL,
    AS_OVERLAP_PHASE_MAD,
    AS_OVERLAP_REGISTERED_MAD,
    AS_OVERLAP_SHIFT_MAD,
    AS_OVERLAP_SEPARATION_CYCLE,
    AS_OVERLAP_TRUNCATED_AUDIT
} AsOverlapLiftDecision;

static const char *overlap_lift_decision_name(AsOverlapLiftDecision d)
{
    switch (d) {
    case AS_OVERLAP_SELECTED: return "selected";
    case AS_OVERLAP_SAME_COMPONENT: return "same_component";
    case AS_OVERLAP_SAME_PROTECTED_GROUP: return "same_protected_group";
    case AS_OVERLAP_WEAK_PAIR_SUPPORT: return "weak_pair_support";
    case AS_OVERLAP_WEAK_FACE_SUPPORT: return "weak_face_support";
    case AS_OVERLAP_ZERO_TURN: return "zero_turn";
    case AS_OVERLAP_TURN_DISAGREEMENT: return "turn_disagreement";
    case AS_OVERLAP_PHASE_RESIDUAL: return "phase_residual";
    case AS_OVERLAP_PHASE_MAD: return "phase_mad";
    case AS_OVERLAP_REGISTERED_MAD: return "registered_mad";
    case AS_OVERLAP_SHIFT_MAD: return "shift_mad";
    case AS_OVERLAP_SEPARATION_CYCLE: return "separation_cycle";
    case AS_OVERLAP_TRUNCATED_AUDIT: return "truncated_audit";
    default: return "unknown";
    }
}

typedef struct {
    int32_t component0, component1;
    double target;
    double confidence;
    size_t source_index;
    int source_is_overlap;
    int priority;
} AsLiftProposal;

static int compare_lift_proposal(const void *pa, const void *pb)
{
    const AsLiftProposal *a = (const AsLiftProposal *)pa;
    const AsLiftProposal *b = (const AsLiftProposal *)pb;
    if (a->priority != b->priority)
        return a->priority > b->priority ? -1 : 1;
    if (a->confidence != b->confidence)
        return a->confidence > b->confidence ? -1 : 1;
    if (a->source_is_overlap != b->source_is_overlap)
        return a->source_is_overlap < b->source_is_overlap ? -1 : 1;
    return a->source_index < b->source_index ? -1 :
           (a->source_index > b->source_index ? 1 : 0);
}

static AsOverlapLiftDecision classify_overlap_lift_evidence(
    const AtlasOverlapAudit *audit, const AtlasOverlapBundle *b,
    int32_t protected_group0, int32_t protected_group1,
    size_t minimum_pairs, size_t minimum_faces,
    double minimum_turn_agreement, double maximum_phase_residual,
    double maximum_phase_mad, double maximum_registered_mad,
    double maximum_shift_mad)
{
    if (audit->pair_buffer_truncated) return AS_OVERLAP_TRUNCATED_AUDIT;
    if (b->component0 == b->component1) return AS_OVERLAP_SAME_COMPONENT;
    if (protected_group0 == protected_group1)
        return AS_OVERLAP_SAME_PROTECTED_GROUP;
    if (b->face_pairs < minimum_pairs) return AS_OVERLAP_WEAK_PAIR_SUPPORT;
    if (b->unique_faces0 < minimum_faces ||
        b->unique_faces1 < minimum_faces)
        return AS_OVERLAP_WEAK_FACE_SUPPORT;
    if (b->turn_mode == 0) return AS_OVERLAP_ZERO_TURN;
    if (b->turn_agreement < minimum_turn_agreement)
        return AS_OVERLAP_TURN_DISAGREEMENT;
    if (fabs(b->phase_residual_median) > maximum_phase_residual)
        return AS_OVERLAP_PHASE_RESIDUAL;
    if (b->phase_residual_mad > maximum_phase_mad)
        return AS_OVERLAP_PHASE_MAD;
    if (b->registered_du_mad > maximum_registered_mad)
        return AS_OVERLAP_REGISTERED_MAD;
    if (!isfinite(b->registered_parameter_shift_median) ||
        b->registered_parameter_shift_mad > maximum_shift_mad)
        return AS_OVERLAP_SHIFT_MAD;
    return AS_OVERLAP_SELECTED;
}

#define AS_LIFT_REFINEMENT_ROUNDS 2

typedef struct {
    AtlasOverlapAudit audit;
    AsOverlapLiftDecision *decision;
    double *absolute_target;
    size_t evidence;
    size_t selected;
} AsLiftRefinementRound;

static double lift_member_value(const AtlasStripMember *member,
                                const double *sample_value)
{
    double value = sample_value[member->value0];
    if (member->value1 >= 0)
        value = (1.0 - member->value_t) * value +
                member->value_t * sample_value[member->value1];
    return value;
}

static int build_lift_bundle_audit(Arena_T arena, const ScrollConfig *cfg,
                                   const PieceSet *ps,
                                   const AtlasCandidateSet *set,
                                   const double *head_sample,
                                   AsLiftBundleAudit **out_audit)
{
    if (set->nbundles == 0 || head_sample == NULL || out_audit == NULL)
        return -1;
    size_t ns = set->problem.nsamples;
    double *sample_phi = (double *)ARENA_ALLOC(arena, ns * sizeof(double));
    double *sample_registered = (double *)ARENA_ALLOC(
        arena, ns * sizeof(double));
    for (size_t i = 0; i < ns; i++) {
        const AtlasCandidateSampleRef *r = &set->sample_ref[i];
        int32_t a = r->mesh_vertex[0], b = r->mesh_vertex[1];
        if (a < 0 || b < 0 || (size_t)a >= ps->nv || (size_t)b >= ps->nv ||
            !isfinite(r->mesh_t) || r->mesh_t < 0.0 || r->mesh_t > 1.0)
            return -1;
        sample_phi[i] = (1.0 - r->mesh_t) * (double)ps->phi[a] +
                        r->mesh_t * (double)ps->phi[b];
        sample_registered[i] =
            (1.0 - r->mesh_t) * (double)ps->uv[(size_t)a * 2] +
            r->mesh_t * (double)ps->uv[(size_t)b * 2];
    }

    AsLiftObservation *observation = (AsLiftObservation *)ARENA_ALLOC(
        arena, set->problem.nmembers * sizeof(AsLiftObservation));
    size_t no = 0;
    const double two_pi = 6.28318530717958647692;
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->cross_sections[c];
        const AtlasStripMember *source = &set->members[cs->first];
        int32_t source_component =
            set->sample_ref[source->value0].mesh_component;
        double source_phi = lift_member_value(source, sample_phi);
        double source_registered = lift_member_value(source, sample_registered);
        double source_head = AtlasStrip_member_value(source, head_sample);
        for (int32_t j = 1; j < cs->count; j++) {
            size_t mi = cs->first + (size_t)j;
            int32_t bi = set->member_bundle[mi];
            if (bi < 0 || (size_t)bi >= set->nbundles) continue;
            const AtlasCandidateBundle *bundle = &set->bundles[bi];
            const AtlasStripMember *target = &set->members[mi];
            int32_t target_component =
                set->sample_ref[target->value0].mesh_component;
            double target_phi = lift_member_value(target, sample_phi);
            double target_registered =
                lift_member_value(target, sample_registered);
            double target_head = AtlasStrip_member_value(target, head_sample);
            double phi0, phi1, registered0, registered1, head0, head1;
            if (source_component == bundle->component0 &&
                target_component == bundle->component1) {
                phi0 = source_phi; phi1 = target_phi;
                registered0 = source_registered;
                registered1 = target_registered;
                head0 = source_head; head1 = target_head;
            } else if (source_component == bundle->component1 &&
                       target_component == bundle->component0) {
                phi0 = target_phi; phi1 = source_phi;
                registered0 = target_registered;
                registered1 = source_registered;
                head0 = target_head; head1 = source_head;
            } else return -1;
            double turn_value = (phi1 - phi0) / two_pi;
            long turn = lround(turn_value);
            if (turn < INT32_MIN || turn > INT32_MAX) return -1;
            observation[no].bundle = bi;
            observation[no].turn = (int32_t)turn;
            observation[no].phase_residual =
                phi1 - phi0 - two_pi * (double)turn;
            observation[no].registered_du = registered1 - registered0;
            observation[no].head_du = head1 - head0;
            no++;
        }
    }
    qsort(observation, no, sizeof(AsLiftObservation),
          compare_lift_observation);
    AsLiftBundleAudit *audit = (AsLiftBundleAudit *)ARENA_CALLOC(
        arena, set->nbundles, sizeof(AsLiftBundleAudit));
    for (size_t i = 0; i < set->nbundles; i++) {
        audit[i].phase_residual_median = NAN;
        audit[i].phase_residual_mad = NAN;
        audit[i].registered_du_median = NAN;
        audit[i].registered_du_mad = NAN;
        audit[i].head_du_median = NAN;
        audit[i].head_du_mad = NAN;
    }
    double *scratch = (double *)ARENA_ALLOC(
        arena, (no ? no : 1) * sizeof(double));
    size_t first = 0;
    while (first < no) {
        size_t last = first + 1;
        while (last < no && observation[last].bundle == observation[first].bundle)
            last++;
        int32_t bi = observation[first].bundle;
        AsLiftBundleAudit *a = &audit[bi];
        a->observations = last - first;
        size_t mode_count = 0;
        for (size_t r0 = first; r0 < last;) {
            size_t r1 = r0 + 1;
            while (r1 < last && observation[r1].turn == observation[r0].turn)
                r1++;
            if (r1 - r0 > mode_count) {
                mode_count = r1 - r0;
                a->turn_mode = observation[r0].turn;
            }
            r0 = r1;
        }
        a->turn_agreement = (double)mode_count / (double)a->observations;
        size_t count = 0;
        for (size_t i = first; i < last; i++)
            if (observation[i].turn == a->turn_mode)
                scratch[count++] = observation[i].phase_residual;
        a->phase_residual_median = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            if (observation[i].turn == a->turn_mode)
                scratch[count++] = fabs(observation[i].phase_residual -
                                        a->phase_residual_median);
        a->phase_residual_mad = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = observation[i].registered_du;
        a->registered_du_median = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = fabs(observation[i].registered_du -
                                    a->registered_du_median);
        a->registered_du_mad = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = observation[i].head_du;
        a->head_du_median = lift_median(scratch, count);
        count = 0;
        for (size_t i = first; i < last; i++)
            scratch[count++] = fabs(observation[i].head_du -
                                    a->head_du_median);
        a->head_du_mad = lift_median(scratch, count);
        const AtlasCandidateBundle *b = &set->bundles[bi];
        double registered_tolerance = fmax(
            32.0, 2.0 * cfg->refine.initial_coherence_mad_limit);
        a->eligible = b->decision == ATLAS_BUNDLE_SEPARATE_SHEET &&
            a->turn_mode != 0 && a->turn_agreement >= 0.95 &&
            fabs(a->phase_residual_median) <= 0.10 &&
            a->phase_residual_mad <= 0.05 &&
            isfinite(b->initial_residual_median) &&
            fabs(a->registered_du_median - b->initial_residual_median) <=
                registered_tolerance;
        first = last;
    }
    *out_audit = audit;
    return 0;
}

static int run_component_lift_stage(Arena_T arena, const ScrollConfig *cfg,
                                    const PieceSet *ps,
                                    const ScaffoldCalib *cal,
                                    const AtlasCandidateSet *set,
                                    const double *head_sample,
                                    const double *head_field,
                                    double **out_lifted_field)
{
    if (out_lifted_field == NULL) return -1;
    *out_lifted_field = NULL;
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "lifted_components") != 0) return -1;
    AsLiftBundleAudit *audit = NULL;
    if (build_lift_bundle_audit(arena, cfg, ps, set, head_sample,
                                &audit) != 0)
        return -1;

    double *component_weight = (double *)ARENA_CALLOC(
        arena, set->mesh_components, sizeof(double));
    size_t *component_faces = (size_t *)ARENA_CALLOC(
        arena, set->mesh_components, sizeof(size_t));
    size_t *component_vertices = (size_t *)ARENA_CALLOC(
        arena, set->mesh_components, sizeof(size_t));
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t c = set->vertex_mesh_component[i];
        if (c >= 0) {
            if ((size_t)c >= set->mesh_components) return -1;
            component_vertices[c]++;
        }
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t ia = ps->faces[f * 3 + 0];
        int32_t ib = ps->faces[f * 3 + 1];
        int32_t ic = ps->faces[f * 3 + 2];
        int32_t component = set->vertex_mesh_component[ia];
        if (component < 0 || set->vertex_mesh_component[ib] != component ||
            set->vertex_mesh_component[ic] != component)
            return -1;
        const float *a = &ps->verts[(size_t)ia * 3];
        const float *b = &ps->verts[(size_t)ib * 3];
        const float *c = &ps->verts[(size_t)ic * 3];
        double x0 = (double)b[0] - (double)a[0];
        double x1 = (double)b[1] - (double)a[1];
        double x2 = (double)b[2] - (double)a[2];
        double y0 = (double)c[0] - (double)a[0];
        double y1 = (double)c[1] - (double)a[1];
        double y2 = (double)c[2] - (double)a[2];
        double cx = x1 * y2 - x2 * y1;
        double cy = x2 * y0 - x0 * y2;
        double cz = x0 * y1 - x1 * y0;
        component_weight[component] +=
            0.5 * sqrt(cx * cx + cy * cy + cz * cz);
        component_faces[component]++;
    }

    /*
     * HEAD has the useful local differential; the registered field has a much
     * better absolute winding placement.  Record only the robust constant
     * registered-minus-HEAD gauge of each connected mesh component.  The lift
     * solver may use this to choose the one otherwise-free translation of a
     * separation graph, never to deform a component or change a relative
     * separation constraint.
     */
    size_t *prior_offset = (size_t *)ARENA_ALLOC(
        arena, (set->mesh_components + 1) * sizeof(size_t));
    size_t *prior_cursor = (size_t *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(size_t));
    prior_offset[0] = 0;
    size_t max_component_vertices = 0;
    for (size_t c = 0; c < set->mesh_components; c++) {
        prior_offset[c + 1] = prior_offset[c] + component_vertices[c];
        prior_cursor[c] = prior_offset[c];
        if (component_vertices[c] > max_component_vertices)
            max_component_vertices = component_vertices[c];
    }
    size_t used_vertex_count = prior_offset[set->mesh_components];
    double *prior_observation = (double *)ARENA_ALLOC(
        arena, (used_vertex_count ? used_vertex_count : 1) * sizeof(double));
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t c = set->vertex_mesh_component[i];
        if (c < 0) continue;
        double value = (double)ps->uv[i * 2] - head_field[i];
        if ((size_t)c >= set->mesh_components || !isfinite(value))
            return -1;
        prior_observation[prior_cursor[c]++] = value;
    }
    double *component_prior = (double *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(double));
    double *component_prior_mad = (double *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(double));
    double *prior_scratch = (double *)ARENA_ALLOC(
        arena, (max_component_vertices ? max_component_vertices : 1) *
                   sizeof(double));
    for (size_t c = 0; c < set->mesh_components; c++) {
        size_t count = component_vertices[c];
        if (count == 0) return -1;
        double *value = &prior_observation[prior_offset[c]];
        component_prior[c] = lift_median(value, count);
        for (size_t i = 0; i < count; i++)
            prior_scratch[i] = fabs(value[i] - component_prior[c]);
        component_prior_mad[c] = lift_median(prior_scratch, count);
    }

    AtlasComponentLiftEdge *protected_edge =
        (AtlasComponentLiftEdge *)ARENA_ALLOC(
            arena, (set->nbundles ? set->nbundles : 1) *
                       sizeof(AtlasComponentLiftEdge));
    size_t nprotected = 0;
    for (size_t i = 0; i < set->nbundles; i++) {
        const AtlasCandidateBundle *b = &set->bundles[i];
        int raw_coherent = b->decision == ATLAS_BUNDLE_UNCLASSIFIED &&
            isfinite(b->initial_residual_median) &&
            isfinite(b->initial_residual_mad) &&
            fabs(b->initial_residual_median) <=
                cfg->refine.overlap_residual_limit &&
            b->initial_residual_mad <=
                cfg->refine.overlap_coherence_mad_limit;
        if (b->decision == ATLAS_BUNDLE_SAME_SHEET_SEAM ||
            b->decision == ATLAS_BUNDLE_DELAMINATION || raw_coherent) {
            AtlasComponentLiftEdge *e = &protected_edge[nprotected++];
            e->component0 = b->component0;
            e->component1 = b->component1;
            e->target = 0.0;
            e->confidence = fmax(0.0, b->evidence);
            e->source = (int32_t)i;
        }
    }

    int32_t *prior_protected_group = NULL;
    size_t prior_protected_groups = 0;
    double *protected_group_prior = NULL;
    double *component_group_prior = NULL;
    if (build_lift_protected_groups(
            arena, set->mesh_components, protected_edge, nprotected,
            component_weight, component_prior, &prior_protected_group,
            &prior_protected_groups, &protected_group_prior,
            &component_group_prior) != 0)
        return -1;

    AtlasSeamAudit seam;
    if (AtlasSeamAudit_build(
            arena, ps, head_field, set->vertex_mesh_component,
            set->mesh_components, &seam) != 0)
        return -1;
    size_t nseam_eligible = 0, nseam_eligible_pairs = 0;
    size_t nseam_selected = 0;
    uint8_t *seam_selected = (uint8_t *)ARENA_CALLOC(
        arena, (seam.nbundles ? seam.nbundles : 1), sizeof(uint8_t));
    uint8_t *seam_cycle = (uint8_t *)ARENA_CALLOC(
        arena, (seam.nbundles ? seam.nbundles : 1), sizeof(uint8_t));
    for (size_t i = 0; i < seam.nbundles; i++) {
        const AtlasSeamBundle *b = &seam.bundles[i];
        int32_t g0 = prior_protected_group[b->component0];
        int32_t g1 = prior_protected_group[b->component1];
        int eligible = b->component0 != b->component1 && g0 != g1 &&
            b->pairs >= 2 && b->target_shift_mad <= 2.0 &&
            fabs(b->phase_residual_median) <= 0.10 &&
            b->phase_residual_mad <= 0.05;
        if (eligible) {
            nseam_eligible++;
            nseam_eligible_pairs += b->pairs;
        }
    }

    /* The raw HEAD field itself tells us where chart images collide.  The
     * registered field and phase are observations used to distinguish two
     * wraps from an ordinary same-sheet contact; neither is copied vertex by
     * vertex into the result. */
    double *axial_v = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    float *registered_u = (float *)ARENA_ALLOC(arena, ps->nv * sizeof(float));
    for (size_t i = 0; i < ps->nv; i++) {
        axial_v[i] = axis_coordinate(cal, &ps->verts[i * 3]);
        registered_u[i] = ps->uv[i * 2];
    }
    AtlasOverlapAudit overlap;
    if (AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, head_field, axial_v,
            registered_u, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &overlap) != 0)
        return -1;
    if (set->nbundles > SIZE_MAX - overlap.nbundles) return -1;
    size_t separation_capacity = set->nbundles + overlap.nbundles;
    AtlasComponentLiftEdge *separation_edge =
        (AtlasComponentLiftEdge *)ARENA_ALLOC(
            arena, (separation_capacity ? separation_capacity : 1) *
                       sizeof(AtlasComponentLiftEdge));
    uint8_t *candidate_selected = (uint8_t *)ARENA_CALLOC(
        arena, (set->nbundles ? set->nbundles : 1), sizeof(uint8_t));
    double *candidate_target = (double *)ARENA_ALLOC(
        arena, (set->nbundles ? set->nbundles : 1) * sizeof(double));
    AsOverlapLiftDecision *overlap_decision =
        (AsOverlapLiftDecision *)ARENA_ALLOC(
            arena, (overlap.nbundles ? overlap.nbundles : 1) *
                       sizeof(AsOverlapLiftDecision));
    double *overlap_target = (double *)ARENA_ALLOC(
        arena, (overlap.nbundles ? overlap.nbundles : 1) * sizeof(double));
    AsLiftProposal *proposal = (AsLiftProposal *)ARENA_ALLOC(
        arena, (separation_capacity ? separation_capacity : 1) *
                   sizeof(AsLiftProposal));
    size_t nproposal = 0, nseparation = 0;
    size_t ncandidate_evidence = 0, ncandidate_selected = 0;
    size_t noverlap_evidence = 0, noverlap_selected = 0;
    for (size_t i = 0; i < set->nbundles; i++) {
        candidate_target[i] = NAN;
        if (audit[i].eligible) {
            ncandidate_evidence++;
            const AtlasCandidateBundle *b = &set->bundles[i];
            int32_t g0 = prior_protected_group[b->component0];
            int32_t g1 = prior_protected_group[b->component1];
            double target = b->initial_residual_median;
            candidate_target[i] = target;
            if (g0 == g1 || !isfinite(target)) continue;
            AsLiftProposal *p = &proposal[nproposal++];
            p->component0 = b->component0;
            p->component1 = b->component1;
            p->target = target;
            p->confidence = fmax(0.0, b->evidence);
            p->source_index = i;
            p->source_is_overlap = 0;
            p->priority = 2;
        }
    }
    const size_t minimum_overlap_pairs = 100;
    const size_t minimum_overlap_faces = 20;
    const double minimum_turn_agreement = 0.95;
    const double maximum_phase_residual = 0.10;
    const double maximum_phase_mad = 0.05;
    const double maximum_registered_mad =
        cfg->refine.initial_coherence_mad_limit;
    const double maximum_shift_mad =
        cfg->refine.initial_coherence_mad_limit;
    for (size_t i = 0; i < overlap.nbundles; i++) {
        const AtlasOverlapBundle *b = &overlap.bundles[i];
        int32_t g0 = prior_protected_group[b->component0];
        int32_t g1 = prior_protected_group[b->component1];
        double target = b->registered_parameter_shift_median;
        overlap_target[i] = target;
        AsOverlapLiftDecision decision = classify_overlap_lift_evidence(
            &overlap, b, g0, g1, minimum_overlap_pairs,
            minimum_overlap_faces, minimum_turn_agreement,
            maximum_phase_residual, maximum_phase_mad,
            maximum_registered_mad, maximum_shift_mad);
        if (decision == AS_OVERLAP_SELECTED) {
            noverlap_evidence++;
            AsLiftProposal *p = &proposal[nproposal++];
            p->component0 = b->component0;
            p->component1 = b->component1;
            p->target = target;
            p->confidence = (double)b->face_pairs;
            p->source_index = i;
            p->source_is_overlap = 1;
            p->priority = 1;
        }
        overlap_decision[i] = decision;
    }

    /* Exact physical seam copies are the strongest relations.  Close those
     * first, then admit classified wrap and overlap relations on the remaining
     * quotient.  Every tier is a confidence-ordered forest: redundant cycles
     * remain measurable diagnostics and cannot destabilize a coherent chart. */
    AsLiftProposal *seam_proposal = (AsLiftProposal *)ARENA_ALLOC(
        arena, (nseam_eligible ? nseam_eligible : 1) *
                   sizeof(AsLiftProposal));
    size_t nseam_proposal = 0;
    for (size_t i = 0; i < seam.nbundles; i++) {
        const AtlasSeamBundle *b = &seam.bundles[i];
        int32_t g0 = prior_protected_group[b->component0];
        int32_t g1 = prior_protected_group[b->component1];
        int eligible = b->component0 != b->component1 && g0 != g1 &&
            b->pairs >= 2 && b->target_shift_mad <= 2.0 &&
            fabs(b->phase_residual_median) <= 0.10 &&
            b->phase_residual_mad <= 0.05;
        if (!eligible) continue;
        AsLiftProposal *p = &seam_proposal[nseam_proposal++];
        p->component0 = b->component0;
        p->component1 = b->component1;
        p->target = b->target_shift_median;
        p->confidence = (double)b->pairs;
        p->source_index = i;
        p->source_is_overlap = 0;
        p->priority = 3;
    }
    qsort(seam_proposal, nseam_proposal, sizeof(AsLiftProposal),
          compare_lift_proposal);
    UnionFind separation_forest = UF_new(
        arena, (int32_t)prior_protected_groups);
    for (size_t i = 0; i < nseam_proposal; i++) {
        const AsLiftProposal *p = &seam_proposal[i];
        int32_t g0 = prior_protected_group[p->component0];
        int32_t g1 = prior_protected_group[p->component1];
        if (uf_find(&separation_forest, g0) ==
            uf_find(&separation_forest, g1)) {
            seam_cycle[p->source_index] = 1;
            continue;
        }
        if (nseparation >= separation_capacity) return -1;
        uf_union(&separation_forest, g0, g1);
        AtlasComponentLiftEdge *e = &separation_edge[nseparation++];
        e->component0 = p->component0;
        e->component1 = p->component1;
        e->target = p->target;
        e->confidence = p->confidence;
        e->source = INT32_MAX - 1;
        seam_selected[p->source_index] = 1;
        nseam_selected++;
    }
    qsort(proposal, nproposal, sizeof(AsLiftProposal),
          compare_lift_proposal);
    for (size_t i = 0; i < nproposal; i++) {
        const AsLiftProposal *p = &proposal[i];
        int32_t g0 = prior_protected_group[p->component0];
        int32_t g1 = prior_protected_group[p->component1];
        if (uf_find(&separation_forest, g0) ==
            uf_find(&separation_forest, g1)) {
            if (p->source_is_overlap)
                overlap_decision[p->source_index] =
                    AS_OVERLAP_SEPARATION_CYCLE;
            continue;
        }
        uf_union(&separation_forest, g0, g1);
        AtlasComponentLiftEdge *e = &separation_edge[nseparation++];
        e->component0 = p->component0;
        e->component1 = p->component1;
        e->target = p->target;
        e->confidence = p->confidence;
        if (p->source_is_overlap) {
            size_t source = set->nbundles + p->source_index;
            e->source = source <= (size_t)INT32_MAX ?
                (int32_t)source : INT32_MAX;
            noverlap_selected++;
        } else {
            e->source = (int32_t)p->source_index;
            candidate_selected[p->source_index] = 1;
            ncandidate_selected++;
        }
    }
    AtlasComponentLiftProblem problem;
    memset(&problem, 0, sizeof problem);
    problem.ncomponents = set->mesh_components;
    problem.component_weight = component_weight;
    /* Registered U supplies relative wrap evidence only.  Giving every
     * disconnected separation graph an absolute registered gauge scatters
     * otherwise-good HEAD charts into the old patchy layout.  With no absolute
     * prior, AtlasComponentLift keeps the largest protected chart in each
     * graph exactly at HEAD and moves only its aliased neighbours. */
    problem.component_prior_shift = NULL;
    problem.protected_edges = protected_edge;
    problem.nprotected_edges = nprotected;
    problem.separation_edges = separation_edge;
    problem.nseparation_edges = nseparation;
    problem.cycle_tolerance = fmax(
        16.0, 2.0 * cfg->refine.initial_coherence_mad_limit);
    double *component_shift = (double *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(double));
    int32_t *protected_group = (int32_t *)ARENA_ALLOC(
        arena, set->mesh_components * sizeof(int32_t));
    AtlasComponentLiftStats stats;
    if (AtlasComponentLift_solve(arena, &problem, component_shift,
                                 protected_group, &stats) != 0)
        return -1;
    if (stats.protected_groups != prior_protected_groups) return -1;
    for (size_t i = 0; i < set->mesh_components; i++)
        if (protected_group[i] != prior_protected_group[i]) return -1;

    AsLiftRefinementRound refinement[AS_LIFT_REFINEMENT_ROUNDS];
    memset(refinement, 0, sizeof refinement);
    size_t completed_refinement_rounds = 0;
    for (size_t round = 0; round < AS_LIFT_REFINEMENT_ROUNDS; round++) {
        double *current_field = (double *)ARENA_ALLOC(
            arena, ps->nv * sizeof(double));
        for (size_t i = 0; i < ps->nv; i++) {
            int32_t c = set->vertex_mesh_component[i];
            current_field[i] = head_field[i] +
                (c >= 0 ? component_shift[c] : 0.0);
        }
        AsLiftRefinementRound *rr = &refinement[round];
        if (AtlasOverlapAudit_build(
                arena, ps->faces, ps->nf, ps->nv, current_field, axial_v,
                registered_u, ps->phi, set->vertex_mesh_component,
                set->mesh_components, &rr->audit) != 0)
            return -1;
        rr->decision = (AsOverlapLiftDecision *)ARENA_ALLOC(
            arena, (rr->audit.nbundles ? rr->audit.nbundles : 1) *
                       sizeof(AsOverlapLiftDecision));
        rr->absolute_target = (double *)ARENA_ALLOC(
            arena, (rr->audit.nbundles ? rr->audit.nbundles : 1) *
                       sizeof(double));
        AsLiftProposal *round_proposal = (AsLiftProposal *)ARENA_ALLOC(
            arena, (rr->audit.nbundles ? rr->audit.nbundles : 1) *
                       sizeof(AsLiftProposal));
        size_t nround_proposal = 0;
        for (size_t i = 0; i < rr->audit.nbundles; i++) {
            const AtlasOverlapBundle *b = &rr->audit.bundles[i];
            int32_t g0 = prior_protected_group[b->component0];
            int32_t g1 = prior_protected_group[b->component1];
            double current_difference =
                component_shift[b->component1] -
                component_shift[b->component0];
            rr->absolute_target[i] = current_difference +
                b->registered_parameter_shift_median;
            AsOverlapLiftDecision decision = classify_overlap_lift_evidence(
                &rr->audit, b, g0, g1, minimum_overlap_pairs,
                minimum_overlap_faces, minimum_turn_agreement,
                maximum_phase_residual, maximum_phase_mad,
                maximum_registered_mad, maximum_shift_mad);
            rr->decision[i] = decision;
            if (decision != AS_OVERLAP_SELECTED) continue;
            rr->evidence++;
            AsLiftProposal *p = &round_proposal[nround_proposal++];
            p->component0 = b->component0;
            p->component1 = b->component1;
            p->target = rr->absolute_target[i];
            p->confidence = (double)b->face_pairs;
            p->source_index = i;
            p->source_is_overlap = 1;
            p->priority = 1;
        }
        qsort(round_proposal, nround_proposal, sizeof(AsLiftProposal),
              compare_lift_proposal);
        for (size_t i = 0; i < nround_proposal; i++) {
            const AsLiftProposal *p = &round_proposal[i];
            int32_t g0 = prior_protected_group[p->component0];
            int32_t g1 = prior_protected_group[p->component1];
            if (uf_find(&separation_forest, g0) ==
                uf_find(&separation_forest, g1)) {
                rr->decision[p->source_index] =
                    AS_OVERLAP_SEPARATION_CYCLE;
                continue;
            }
            if (nseparation >= separation_capacity) return -1;
            uf_union(&separation_forest, g0, g1);
            AtlasComponentLiftEdge *e = &separation_edge[nseparation++];
            e->component0 = p->component0;
            e->component1 = p->component1;
            e->target = p->target;
            e->confidence = p->confidence;
            e->source = INT32_MAX;
            rr->selected++;
        }
        completed_refinement_rounds++;
        if (rr->selected == 0) break;
        problem.nseparation_edges = nseparation;
        if (AtlasComponentLift_solve(arena, &problem, component_shift,
                                     protected_group, &stats) != 0)
            return -1;
        if (stats.protected_groups != prior_protected_groups ||
            stats.conflicting_graphs != 0)
            return -1;
        for (size_t i = 0; i < set->mesh_components; i++)
            if (protected_group[i] != prior_protected_group[i]) return -1;
    }

    double *lifted_field = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    for (size_t i = 0; i < ps->nv; i++) {
        int32_t component = set->vertex_mesh_component[i];
        lifted_field[i] = head_field[i] +
            (component >= 0 ? component_shift[component] : 0.0);
    }
    double *lifted_sample = (double *)ARENA_ALLOC(
        arena, set->problem.nsamples * sizeof(double));
    for (size_t i = 0; i < set->problem.nsamples; i++)
        lifted_sample[i] = head_sample[i] +
            component_shift[set->sample_ref[i].mesh_component];

    AtlasOverlapAudit final_overlap;
    if (AtlasOverlapAudit_build(
            arena, ps->faces, ps->nf, ps->nv, lifted_field, axial_v,
            registered_u, ps->phi, set->vertex_mesh_component,
            set->mesh_components, &final_overlap) != 0)
        return -1;

    double max_local_delta = 0.0, max_protected_residual = 0.0;
    double max_separation_residual = 0.0;
    size_t moved_faces = 0, affected_unclassified = 0;
    size_t affected_unclassified_links = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t v[3] = {ps->faces[f * 3 + 0], ps->faces[f * 3 + 1],
                        ps->faces[f * 3 + 2]};
        for (int e = 0; e < 3; e++) {
            int32_t a = v[e], b = v[(e + 1) % 3];
            double d = (lifted_field[b] - lifted_field[a]) -
                       (head_field[b] - head_field[a]);
            if (fabs(d) > max_local_delta) max_local_delta = fabs(d);
        }
    }
    for (size_t i = 0; i < set->mesh_components; i++)
        if (fabs(component_shift[i]) > 1e-10)
            moved_faces += component_faces[i];
    for (size_t i = 0; i < nprotected; i++) {
        const AtlasComponentLiftEdge *e = &protected_edge[i];
        double r = component_shift[e->component1] -
                   component_shift[e->component0];
        if (fabs(r) > max_protected_residual)
            max_protected_residual = fabs(r);
    }
    for (size_t i = 0; i < nseparation; i++) {
        const AtlasComponentLiftEdge *e = &separation_edge[i];
        double r = component_shift[e->component1] -
                   component_shift[e->component0] - e->target;
        if (fabs(r) > max_separation_residual)
            max_separation_residual = fabs(r);
    }
    for (size_t i = 0; i < set->nbundles; i++) {
        const AtlasCandidateBundle *b = &set->bundles[i];
        double d = component_shift[b->component1] -
                   component_shift[b->component0];
        if (b->decision == ATLAS_BUNDLE_UNCLASSIFIED && fabs(d) > 1e-10) {
            affected_unclassified++;
            affected_unclassified_links += b->cross_links;
        }
    }

    FILE *fp = as_open(dir, "component_lift_constraints.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "bundle,component0,component1,protected_group0,"
                "protected_group1,decision,observations,"
                "turn_mode,turn_agreement,phase_residual_median,"
                "phase_residual_mad,registered_du_median,registered_du_mad,"
                "head_du_median,head_du_mad,initial_du_median,audit_eligible,"
                "selected,group_target,shift_difference,predicted_head_du,"
                "target_residual\n");
    for (size_t i = 0; i < set->nbundles; i++) {
        const AtlasCandidateBundle *b = &set->bundles[i];
        const AsLiftBundleAudit *a = &audit[i];
        double d = component_shift[b->component1] -
                   component_shift[b->component0];
        fprintf(fp,
            "%zu,%d,%d,%d,%d,%s,%zu,%d,",
            i, b->component0, b->component1,
            protected_group[b->component0], protected_group[b->component1],
            bundle_decision_name(b->decision), a->observations,
            a->turn_mode);
        fprintf(fp,
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,"
            "%.17g,%.17g,%.17g,%.17g\n",
            a->turn_agreement, a->phase_residual_median,
            a->phase_residual_mad, a->registered_du_median,
            a->registered_du_mad, a->head_du_median, a->head_du_mad,
            b->initial_residual_median, a->eligible, candidate_selected[i],
            candidate_target[i], d,
            a->head_du_median + d,
            candidate_selected[i] ? d - candidate_target[i] : NAN);
    }
    fclose(fp);

    fp = as_open(dir, "component_overlap_bundles.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "overlap_bundle,component0,component1,protected_group0,"
                "protected_group1,face_pairs,unique_faces0,unique_faces1,"
                "component_faces0,component_faces1,support_fraction0,"
                "support_fraction1,turn_mode,turn_agreement,"
                "phase_residual_median,phase_residual_mad,"
                "parameter_du_median,parameter_du_mad,registered_du_median,"
                "registered_du_mad,required_shift,required_shift_mad,"
                "protected_group_prior_target,group_prior_minus_required,"
                "evidence_pass,selected,reason,"
                "shift_difference,target_residual\n");
    for (size_t i = 0; i < overlap.nbundles; i++) {
        const AtlasOverlapBundle *b = &overlap.bundles[i];
        double support0 = b->component_faces0 > 0 ?
            (double)b->unique_faces0 / (double)b->component_faces0 : 0.0;
        double support1 = b->component_faces1 > 0 ?
            (double)b->unique_faces1 / (double)b->component_faces1 : 0.0;
        double d = component_shift[b->component1] -
                   component_shift[b->component0];
        double group_prior_target =
            protected_group_prior[protected_group[b->component1]] -
            protected_group_prior[protected_group[b->component0]];
        int evidence_pass = overlap_decision[i] == AS_OVERLAP_SELECTED ||
            overlap_decision[i] == AS_OVERLAP_SEPARATION_CYCLE;
        int selected = overlap_decision[i] == AS_OVERLAP_SELECTED;
        fprintf(fp, "%zu,%d,%d,%d,%d,%zu,%zu,%zu,%zu,%zu,",
                i, b->component0, b->component1,
                protected_group[b->component0], protected_group[b->component1],
                b->face_pairs, b->unique_faces0, b->unique_faces1,
                b->component_faces0, b->component_faces1);
        fprintf(fp,
                "%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%s,%.17g,%.17g\n",
                support0, support1, b->turn_mode, b->turn_agreement,
                b->phase_residual_median, b->phase_residual_mad,
                b->parameter_du_median, b->parameter_du_mad,
                b->registered_du_median, b->registered_du_mad,
                overlap_target[i], b->registered_parameter_shift_mad,
                group_prior_target, group_prior_target - overlap_target[i],
                evidence_pass, selected,
                overlap_lift_decision_name(overlap_decision[i]), d,
                evidence_pass ? d - overlap_target[i] : NAN);
    }
    fclose(fp);

    fp = as_open(dir, "component_seam_bundles.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "seam_bundle,component0,component1,protected_group0,"
                "protected_group1,pairs,target_shift_median,target_shift_mad,"
                "distance_median,phase_residual_median,phase_residual_mad,"
                "registered_du_median,registered_du_mad,eligible,"
                "selected,cycle,reason,final_shift_difference,"
                "final_seam_residual\n");
    for (size_t i = 0; i < seam.nbundles; i++) {
        const AtlasSeamBundle *b = &seam.bundles[i];
        int32_t g0 = protected_group[b->component0];
        int32_t g1 = protected_group[b->component1];
        int eligible = b->component0 != b->component1 && g0 != g1 &&
            b->pairs >= 2 && b->target_shift_mad <= 2.0 &&
            fabs(b->phase_residual_median) <= 0.10 &&
            b->phase_residual_mad <= 0.05;
        double d = component_shift[b->component1] -
                   component_shift[b->component0];
        const char *reason = seam_selected[i] ? "selected" :
            (seam_cycle[i] ? "separation_cycle" :
             (eligible ? "not_selected" : "rejected"));
        fprintf(fp, "%zu,%d,%d,%d,%d,%zu,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%d,%d,%d,%s,%.17g,%.17g\n",
                i, b->component0, b->component1, g0, g1, b->pairs,
                b->target_shift_median, b->target_shift_mad,
                b->distance_median, b->phase_residual_median,
                b->phase_residual_mad, b->registered_du_median,
                b->registered_du_mad, eligible, seam_selected[i],
                seam_cycle[i], reason, d,
                d - b->target_shift_median);
    }
    fclose(fp);

    for (size_t round = 0; round < completed_refinement_rounds; round++) {
        char name[96];
        if (snprintf(name, sizeof name,
                     "component_overlap_refine_round_%zu.csv", round + 1) < 0)
            return -1;
        fp = as_open(dir, name, "wb");
        if (fp == NULL) return -1;
        fprintf(fp, "overlap_bundle,component0,component1,protected_group0,"
                    "protected_group1,face_pairs,unique_faces0,unique_faces1,"
                    "turn_mode,turn_agreement,phase_residual_median,"
                    "phase_residual_mad,parameter_du_median,"
                    "registered_du_median,additional_required_shift,"
                    "additional_required_shift_mad,absolute_target,"
                    "evidence_pass,selected,reason,final_shift_difference,"
                    "final_target_residual\n");
        AsLiftRefinementRound *rr = &refinement[round];
        for (size_t i = 0; i < rr->audit.nbundles; i++) {
            const AtlasOverlapBundle *b = &rr->audit.bundles[i];
            double d = component_shift[b->component1] -
                       component_shift[b->component0];
            int evidence_pass = rr->decision[i] == AS_OVERLAP_SELECTED ||
                rr->decision[i] == AS_OVERLAP_SEPARATION_CYCLE;
            int selected = rr->decision[i] == AS_OVERLAP_SELECTED;
            fprintf(fp, "%zu,%d,%d,%d,%d,%zu,%zu,%zu,%d,",
                    i, b->component0, b->component1,
                    protected_group[b->component0],
                    protected_group[b->component1], b->face_pairs,
                    b->unique_faces0, b->unique_faces1, b->turn_mode);
            fprintf(fp,
                    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%d,%d,%s,%.17g,%.17g\n",
                    b->turn_agreement, b->phase_residual_median,
                    b->phase_residual_mad, b->parameter_du_median,
                    b->registered_du_median,
                    b->registered_parameter_shift_median,
                    b->registered_parameter_shift_mad,
                    rr->absolute_target[i], evidence_pass, selected,
                    overlap_lift_decision_name(rr->decision[i]), d,
                    evidence_pass ? d - rr->absolute_target[i] : NAN);
        }
        fclose(fp);
    }

    fp = as_open(dir, "component_overlap_final.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "overlap_bundle,component0,component1,protected_group0,"
                "protected_group1,face_pairs,unique_faces0,unique_faces1,"
                "turn_mode,turn_agreement,phase_residual_median,"
                "phase_residual_mad,parameter_du_median,registered_du_median,"
                "additional_required_shift,additional_required_shift_mad,"
                "evidence_pass,separation_graph_cycle,reason\n");
    for (size_t i = 0; i < final_overlap.nbundles; i++) {
        const AtlasOverlapBundle *b = &final_overlap.bundles[i];
        int32_t g0 = protected_group[b->component0];
        int32_t g1 = protected_group[b->component1];
        AsOverlapLiftDecision decision = classify_overlap_lift_evidence(
            &final_overlap, b, g0, g1, minimum_overlap_pairs,
            minimum_overlap_faces, minimum_turn_agreement,
            maximum_phase_residual, maximum_phase_mad,
            maximum_registered_mad, maximum_shift_mad);
        int evidence_pass = decision == AS_OVERLAP_SELECTED;
        int graph_cycle = evidence_pass &&
            uf_find(&separation_forest, g0) ==
            uf_find(&separation_forest, g1);
        if (graph_cycle) decision = AS_OVERLAP_SEPARATION_CYCLE;
        fprintf(fp, "%zu,%d,%d,%d,%d,%zu,%zu,%zu,%d,",
                i, b->component0, b->component1, g0, g1, b->face_pairs,
                b->unique_faces0, b->unique_faces1, b->turn_mode);
        fprintf(fp,
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%s\n",
                b->turn_agreement, b->phase_residual_median,
                b->phase_residual_mad, b->parameter_du_median,
                b->registered_du_median,
                b->registered_parameter_shift_median,
                b->registered_parameter_shift_mad, evidence_pass,
                graph_cycle, overlap_lift_decision_name(decision));
    }
    fclose(fp);

    fp = as_open(dir, "component_lift_components.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "component,protected_group,vertices,faces,area,"
                "registered_head_prior,registered_head_prior_mad,"
                "protected_group_prior,shift,shift_minus_component_prior,"
                "shift_minus_group_prior\n");
    for (size_t i = 0; i < set->mesh_components; i++)
        fprintf(fp, "%zu,%d,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g\n",
                i, protected_group[i], component_vertices[i],
                component_faces[i], component_weight[i], component_prior[i],
                component_prior_mad[i], component_group_prior[i],
                component_shift[i], component_shift[i] - component_prior[i],
                component_shift[i] - component_group_prior[i]);
    fclose(fp);

    ScrollConfig lifted_cfg = *cfg;
    lifted_cfg.out_dir = dir;
    if (write_strokes_obj(dir, "lifted_parameter.obj", set,
                          lifted_sample, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_head_flat.obj", ps, cal,
                              head_field, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_lifted_flat.obj", ps, cal,
                              lifted_field, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_lifted_world_uv.obj", ps, cal,
                              lifted_field, 0) != 0 ||
        write_global_mesh_xyz_obj(dir, "atlas_head_flat_xyzcolor.obj", ps,
                                  cal, head_field) != 0 ||
        write_global_mesh_xyz_obj(dir, "atlas_lifted_flat_xyzcolor.obj", ps,
                                  cal, lifted_field) != 0 ||
        export_field_placed(&lifted_cfg, ps, cal, lifted_field) != 0)
        return -1;

    fp = as_open(dir, "component_lift_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n  \"baseline\": \"immutable_HEAD\",\n"
        "  \"operation\": \"constant_component_gauge_only\",\n"
        "  \"absolute_gauge\": "
        "\"largest_protected_chart_kept_at_HEAD\",\n"
        "  \"faces_before\": %zu,\n  \"faces_after\": %zu,\n"
        "  \"candidate_evidence_bundles\": %zu,\n"
        "  \"candidate_selected_bundles\": %zu,\n"
        "  \"seam_boundary_vertices\": %zu,\n"
        "  \"seam_mutual_pairs\": %zu,\n"
        "  \"seam_component_bundles\": %zu,\n"
        "  \"seam_eligible_bundles\": %zu,\n"
        "  \"seam_eligible_pairs\": %zu,\n"
        "  \"seam_selected_bundles\": %zu,\n"
        "  \"seam_constraints_applied\": true,\n"
        "  \"initial_exact_overlap_face_pairs\": %zu,\n"
        "  \"initial_same_component_overlap_pairs\": %zu,\n"
        "  \"initial_cross_component_overlap_pairs\": %zu,\n"
        "  \"initial_overlap_component_bundles\": %zu,\n"
        "  \"initial_overlap_evidence_bundles\": %zu,\n"
        "  \"initial_overlap_selected_bundles\": %zu,\n"
        "  \"overlap_pair_buffer_truncated\": %s,\n"
        "  \"completed_refinement_rounds\": %zu,\n"
        "  \"refinement_round_1_exact_pairs\": %zu,\n"
        "  \"refinement_round_1_evidence_bundles\": %zu,\n"
        "  \"refinement_round_1_selected_bundles\": %zu,\n"
        "  \"refinement_round_2_exact_pairs\": %zu,\n"
        "  \"refinement_round_2_evidence_bundles\": %zu,\n"
        "  \"refinement_round_2_selected_bundles\": %zu,\n"
        "  \"final_exact_overlap_face_pairs\": %zu,\n"
        "  \"final_same_component_overlap_pairs\": %zu,\n"
        "  \"final_cross_component_overlap_pairs\": %zu,\n"
        "  \"final_overlap_component_bundles\": %zu,\n"
        "  \"total_separation_edges\": %zu,\n"
        "  \"protected_edges\": %zu,\n"
        "  \"protected_groups\": %zu,\n"
        "  \"separation_graphs\": %zu,\n"
        "  \"solved_graphs\": %zu,\n"
        "  \"conflicting_graphs\": %zu,\n"
        "  \"moved_components\": %zu,\n"
        "  \"moved_faces\": %zu,\n"
        "  \"moved_area\": %.17g,\n"
        "  \"max_abs_shift\": %.17g,\n"
        "  \"max_cycle_residual\": %.17g,\n"
        "  \"max_local_edge_delta\": %.17g,\n"
        "  \"max_protected_residual\": %.17g,\n"
        "  \"max_separation_residual\": %.17g,\n"
        "  \"affected_unclassified_bundles\": %zu,\n"
        "  \"affected_unclassified_links\": %zu\n}\n",
        ps->nf, ps->nf, ncandidate_evidence, ncandidate_selected,
        seam.boundary_vertices, seam.mutual_pairs, seam.nbundles,
        nseam_eligible, nseam_eligible_pairs, nseam_selected,
        overlap.exact_face_pairs, overlap.same_component_pairs,
        overlap.cross_component_pairs, overlap.nbundles,
        noverlap_evidence, noverlap_selected,
        overlap.pair_buffer_truncated ? "true" : "false",
        completed_refinement_rounds,
        refinement[0].audit.exact_face_pairs, refinement[0].evidence,
        refinement[0].selected,
        refinement[1].audit.exact_face_pairs, refinement[1].evidence,
        refinement[1].selected,
        final_overlap.exact_face_pairs, final_overlap.same_component_pairs,
        final_overlap.cross_component_pairs, final_overlap.nbundles,
        nseparation, nprotected, stats.protected_groups,
        stats.separation_graphs, stats.solved_graphs,
        stats.conflicting_graphs, stats.moved_components, moved_faces,
        stats.moved_weight, stats.max_abs_shift,
        stats.max_cycle_residual, max_local_delta,
        max_protected_residual, max_separation_residual,
        affected_unclassified, affected_unclassified_links);
    fclose(fp);
    fprintf(stderr,
        "[atlas_strip_scroll] component lift: candidate=%zu/%zu "
        "overlap=%zu/%zu exact_pairs=%zu protected=%zu "
        "graphs=%zu/%zu conflicts=%zu moved=%zu components/%zu faces "
        "local_delta=%.3e protected_residual=%.3e\n",
        ncandidate_selected, ncandidate_evidence,
        noverlap_selected, noverlap_evidence, overlap.exact_face_pairs,
        nprotected, stats.solved_graphs,
        stats.separation_graphs, stats.conflicting_graphs,
        stats.moved_components, moved_faces, max_local_delta,
        max_protected_residual);
    int valid = stats.conflicting_graphs == 0 &&
        max_local_delta <= 1e-8 && max_protected_residual <= 1e-8 &&
        max_separation_residual <= problem.cycle_tolerance;
    if (valid) *out_lifted_field = lifted_field;
    return valid ? 0 : -1;
}

static int export_field_placed(const ScrollConfig *cfg, const PieceSet *ps,
                               const ScaffoldCalib *cal,
                               const double *field_u)
{
    char placed[AS_PATH_CAP], source[AS_PATH_CAP], destination[AS_PATH_CAP];
    if (snprintf(placed, sizeof placed, "%s/placed", cfg->out_dir) < 0)
        return -1;
    if (snprintf(source, sizeof source, "%s/placed_index.json",
                 cfg->placed_dir) < 0 ||
        snprintf(destination, sizeof destination, "%s/placed_index.json",
                 placed) < 0 || copy_file(source, destination) != 0)
        return -1;

    static const char *required_suffix[] = {
        "mesh.obj", "facekeep.u8", "group.i32", "uvphi_raw.f32"
    };
    for (size_t c = 0; c < ps->n_cubes; c++) {
        for (size_t s = 0; s < sizeof required_suffix /
                                sizeof required_suffix[0]; s++) {
            if (snprintf(source, sizeof source, "%s/%s_%s",
                         cfg->placed_dir, ps->ids[c], required_suffix[s]) < 0 ||
                snprintf(destination, sizeof destination, "%s/%s_%s",
                         placed, ps->ids[c], required_suffix[s]) < 0 ||
                copy_file(source, destination) != 0)
                return -1;
        }
        if (snprintf(destination, sizeof destination, "%s/%s_uvphi.f32",
                     placed, ps->ids[c]) < 0 ||
            ves_ensure_parent_dir(destination) != 0)
            return -1;
        FILE *fp = fopen(destination, "wb");
        if (fp == NULL) return -1;
        for (size_t i = ps->cube_voff[c]; i < ps->cube_voff[c + 1]; i++) {
            float uvphi[3] = {
                (float)field_u[i],
                (float)axis_coordinate(cal, &ps->verts[i * 3]),
                ps->phi[i]
            };
            if (fwrite(uvphi, sizeof(float), 3, fp) != 3) {
                fclose(fp);
                return -1;
            }
        }
        if (fclose(fp) != 0) return -1;
    }
    return 0;
}

static int compare_double_value(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static void write_csv_quoted(FILE *fp, const char *text)
{
    fputc('"', fp);
    if (text != NULL) {
        for (const char *p = text; *p != '\0'; p++) {
            if (*p == '"') fputc('"', fp);
            fputc(*p, fp);
        }
    }
    fputc('"', fp);
}

static int initialize_pipeline_trace(ScrollConfig *cfg,
                                     char trace_dir[AS_PATH_CAP],
                                     int *next_frame)
{
    if (cfg == NULL || trace_dir == NULL || next_frame == NULL ||
        as_path(trace_dir, cfg->out_dir, "pipeline_trace") != 0)
        return -1;
    FILE *fp = as_open(trace_dir, "index.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "frame,obj,metadata,domain,stage,source_output_dir\n");
    if (fclose(fp) != 0) return -1;
    fp = as_open(trace_dir, "README.txt", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "Files in this directory are one chronological execution trace.\n"
        "Sort by the numeric prefix, or use index.csv.\n"
        "mesh OBJ: triangles in (u, axial, 0).\n"
        "strokes OBJ: ordered coarse samples joined by line segments.\n"
        "Every OBJ has a same-named JSON sidecar with validity and scale "
        "diagnostics.\n"
        "The index is authoritative if this output directory is reused.\n");
    if (fclose(fp) != 0) return -1;
    *next_frame = 0;
    cfg->pipeline_trace_dir = trace_dir;
    cfg->pipeline_trace_next_frame = next_frame;
    return 0;
}

static int append_pipeline_trace_index(const ScrollConfig *cfg, int frame,
                                       const char *obj_name,
                                       const char *json_name,
                                       const char *domain,
                                       const char *stage)
{
    FILE *fp = as_open(cfg->pipeline_trace_dir, "index.csv", "ab");
    if (fp == NULL) return -1;
    fprintf(fp, "%d,", frame);
    write_csv_quoted(fp, obj_name);
    fputc(',', fp);
    write_csv_quoted(fp, json_name);
    fputc(',', fp);
    write_csv_quoted(fp, domain);
    fputc(',', fp);
    write_csv_quoted(fp, stage);
    fputc(',', fp);
    write_csv_quoted(fp, cfg->out_dir);
    fputc('\n', fp);
    return fclose(fp) == 0 ? 0 : -1;
}

typedef struct {
    size_t finite_vertices;
    size_t nonfinite_vertices;
    size_t comparable_faces;
    size_t relative_flips;
    size_t degenerate_faces;
    size_t aspect_count;
    size_t aspect_over_50;
    size_t stretch_count;
    size_t stretch_over_4;
    double u_min;
    double u_max;
    double delta_raw_rms;
    double delta_raw_max;
    double stretch_p50;
    double stretch_p95;
    double stretch_p99;
    double stretch_max;
    double aspect_p99;
    double aspect_p999;
    double aspect_max;
} AsPipelineMeshQuality;

static int measure_pipeline_mesh_quality(
    const ScrollConfig *cfg, const PieceSet *ps, const ScaffoldCalib *cal,
    const double *u, AsPipelineMeshQuality *quality)
{
    if (cfg == NULL || ps == NULL || cal == NULL || u == NULL ||
        quality == NULL || (ps->nf > 0 && ps->nf > SIZE_MAX / 3))
        return -1;
    memset(quality, 0, sizeof(*quality));
    quality->u_min = DBL_MAX;
    quality->u_max = -DBL_MAX;
    size_t capacity = 3 * ps->nf;
    double *stretch = (double *)malloc(
        (capacity ? capacity : 1) * sizeof(*stretch));
    double *aspect = (double *)malloc(
        (ps->nf ? ps->nf : 1) * sizeof(*aspect));
    if (stretch == NULL || aspect == NULL) {
        free(stretch);
        free(aspect);
        return -1;
    }
    double delta2 = 0.0;
    size_t ndelta = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        if (!isfinite(u[i])) {
            quality->nonfinite_vertices++;
            continue;
        }
        quality->finite_vertices++;
        if (u[i] < quality->u_min) quality->u_min = u[i];
        if (u[i] > quality->u_max) quality->u_max = u[i];
        if (cfg->pipeline_trace_raw_u != NULL &&
            isfinite(cfg->pipeline_trace_raw_u[i])) {
            double delta = u[i] - (double)cfg->pipeline_trace_raw_u[i];
            delta2 += delta * delta;
            if (fabs(delta) > quality->delta_raw_max)
                quality->delta_raw_max = fabs(delta);
            ndelta++;
        }
    }
    if (ndelta > 0)
        quality->delta_raw_rms = sqrt(delta2 / (double)ndelta);
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t vertex[3] = {
            ps->faces[3 * f], ps->faces[3 * f + 1], ps->faces[3 * f + 2]
        };
        double axial[3];
        int valid = 1;
        for (int k = 0; k < 3; k++) {
            if (vertex[k] < 0 || (size_t)vertex[k] >= ps->nv ||
                !isfinite(u[vertex[k]])) {
                valid = 0;
                break;
            }
            axial[k] = axis_coordinate(cal,
                                       &ps->verts[(size_t)vertex[k] * 3]);
        }
        if (!valid) continue;
        double new_det = (u[vertex[1]] - u[vertex[0]]) *
                         (axial[2] - axial[0]) -
                         (u[vertex[2]] - u[vertex[0]]) *
                         (axial[1] - axial[0]);
        if (fabs(new_det) < 1e-10) quality->degenerate_faces++;
        if (cfg->pipeline_trace_raw_u != NULL &&
            isfinite(cfg->pipeline_trace_raw_u[vertex[0]]) &&
            isfinite(cfg->pipeline_trace_raw_u[vertex[1]]) &&
            isfinite(cfg->pipeline_trace_raw_u[vertex[2]])) {
            double raw_det =
                ((double)cfg->pipeline_trace_raw_u[vertex[1]] -
                 (double)cfg->pipeline_trace_raw_u[vertex[0]]) *
                    (axial[2] - axial[0]) -
                ((double)cfg->pipeline_trace_raw_u[vertex[2]] -
                 (double)cfg->pipeline_trace_raw_u[vertex[0]]) *
                    (axial[1] - axial[0]);
            if (fabs(raw_det) >= 1e-10) {
                quality->comparable_faces++;
                if (raw_det * new_det < 0.0) quality->relative_flips++;
            }
        }
        double uv_edge_max = 0.0;
        for (int e = 0; e < 3; e++) {
            int32_t a = vertex[e], b = vertex[(e + 1) % 3];
            const float *pa = &ps->verts[(size_t)a * 3];
            const float *pb = &ps->verts[(size_t)b * 3];
            double d0 = (double)pa[0] - (double)pb[0];
            double d1 = (double)pa[1] - (double)pb[1];
            double d2 = (double)pa[2] - (double)pb[2];
            double length = sqrt(d0 * d0 + d1 * d1 + d2 * d2);
            if (length <= 1e-12) continue;
            double du = u[a] - u[b];
            double dv = axial[e] - axial[(e + 1) % 3];
            double uv_length = sqrt(du * du + dv * dv);
            if (uv_length > uv_edge_max) uv_edge_max = uv_length;
            double ratio = uv_length / length;
            if (!isfinite(ratio)) continue;
            stretch[quality->stretch_count++] = ratio;
            if (ratio > 4.0) quality->stretch_over_4++;
        }
        if (fabs(new_det) >= 1e-10 && isfinite(uv_edge_max)) {
            /* Lmax / minimum altitude == Lmax^2 / (2*area), while
             * fabs(new_det) is exactly twice the UV triangle area. */
            double ratio = uv_edge_max * uv_edge_max / fabs(new_det);
            if (isfinite(ratio)) {
                aspect[quality->aspect_count++] = ratio;
                if (ratio > 50.0) quality->aspect_over_50++;
            }
        }
    }
    qsort(stretch, quality->stretch_count, sizeof(*stretch),
          compare_double_value);
    size_t n = quality->stretch_count;
    if (n > 0) {
        quality->stretch_p50 = stretch[(size_t)(0.50 * (double)(n - 1))];
        quality->stretch_p95 = stretch[(size_t)(0.95 * (double)(n - 1))];
        quality->stretch_p99 = stretch[(size_t)(0.99 * (double)(n - 1))];
        quality->stretch_max = stretch[n - 1];
    }
    qsort(aspect, quality->aspect_count, sizeof(*aspect),
          compare_double_value);
    n = quality->aspect_count;
    if (n > 0) {
        quality->aspect_p99 = aspect[(size_t)(0.99 * (double)(n - 1))];
        quality->aspect_p999 = aspect[(size_t)(0.999 * (double)(n - 1))];
        quality->aspect_max = aspect[n - 1];
    }
    free(stretch);
    free(aspect);
    return 0;
}

static int write_pipeline_trace_mesh(
    const ScrollConfig *cfg, const PieceSet *ps, const ScaffoldCalib *cal,
    const char *stage, const double *u,
    const AtlasRegisterIterationStats *registration_iteration,
    const AtlasFieldApplyStats *field_stats)
{
    if (cfg == NULL || cfg->pipeline_trace_dir == NULL ||
        cfg->pipeline_trace_next_frame == NULL)
        return 0;
    if (ps == NULL || cal == NULL || stage == NULL || u == NULL) return -1;
    int frame = *cfg->pipeline_trace_next_frame;
    char obj_name[256], json_name[256];
    int no = snprintf(obj_name, sizeof(obj_name), "%04d_%s.obj", frame, stage);
    int nj = snprintf(json_name, sizeof(json_name), "%04d_%s.json", frame,
                      stage);
    if (no < 0 || (size_t)no >= sizeof(obj_name) ||
        nj < 0 || (size_t)nj >= sizeof(json_name) ||
        write_global_mesh_obj(cfg->pipeline_trace_dir, obj_name, ps, cal,
                              u, 1) != 0)
        return -1;
    AsPipelineMeshQuality q;
    if (measure_pipeline_mesh_quality(cfg, ps, cal, u, &q) != 0) return -1;
    FILE *fp = as_open(cfg->pipeline_trace_dir, json_name, "wb");
    if (fp == NULL) return -1;
    char umin[64], umax[64], uspan[64];
    as_json_number(umin, q.finite_vertices ? q.u_min : (double)NAN);
    as_json_number(umax, q.finite_vertices ? q.u_max : (double)NAN);
    as_json_number(uspan, q.finite_vertices ? q.u_max - q.u_min :
                                              (double)NAN);
    fprintf(fp, "{\n  \"frame\": %d,\n  \"domain\": \"mesh\",\n"
                "  \"stage\": ", frame);
    as_write_json_string(fp, stage);
    fprintf(fp, ",\n  \"source_output_dir\": ");
    as_write_json_string(fp, cfg->out_dir);
    fprintf(fp,
        ",\n  \"representation\": \"triangle mesh in (u, axial, 0)\",\n"
        "  \"vertices\": %zu,\n  \"faces\": %zu,\n"
        "  \"finite_vertices\": %zu,\n"
        "  \"nonfinite_vertices\": %zu,\n"
        "  \"u_min\": %s,\n  \"u_max\": %s,\n  \"u_span\": %s,\n"
        "  \"delta_from_raw_rms\": %.17g,\n"
        "  \"delta_from_raw_max\": %.17g,\n"
        "  \"comparable_faces\": %zu,\n"
        "  \"relative_flips_from_raw\": %zu,\n"
        "  \"degenerate_faces\": %zu,\n"
        "  \"aspect\": {\"count\": %zu, \"p99\": %.17g, "
        "\"p999\": %.17g, \"max\": %.17g, \"over_50\": %zu},\n"
        "  \"stretch\": {\"count\": %zu, \"p50\": %.17g, "
        "\"p95\": %.17g, \"p99\": %.17g, \"max\": %.17g, "
        "\"over_4\": %zu}",
        ps->nv, ps->nf, q.finite_vertices, q.nonfinite_vertices,
        umin, umax, uspan, q.delta_raw_rms, q.delta_raw_max,
        q.comparable_faces, q.relative_flips, q.degenerate_faces,
        q.aspect_count, q.aspect_p99, q.aspect_p999, q.aspect_max,
        q.aspect_over_50,
        q.stretch_count, q.stretch_p50, q.stretch_p95, q.stretch_p99,
        q.stretch_max, q.stretch_over_4);
    if (registration_iteration != NULL) {
        fprintf(fp,
            ",\n  \"registration_irls\": {\"iteration\": %d, "
            "\"frozen_energy_before\": %.17g, "
            "\"frozen_energy_after\": %.17g, "
            "\"robust_energy_before\": %.17g, "
            "\"robust_energy_after\": %.17g, "
            "\"lambda_min_equilibrated\": %.17g, "
            "\"lambda_max_equilibrated\": %.17g, "
            "\"condition_equilibrated\": %.17g, "
            "\"downweighted_edges\": %zu, "
            "\"weld_residual_rms\": %.17g, "
            "\"radial_residual_rms\": %.17g}",
            registration_iteration->iteration,
            registration_iteration->frozen_energy_before,
            registration_iteration->frozen_energy_after,
            registration_iteration->robust_energy_before,
            registration_iteration->robust_energy_after,
            registration_iteration->lambda_min,
            registration_iteration->lambda_max,
            registration_iteration->condition_estimate,
            registration_iteration->downweighted_edges,
            registration_iteration->weld_residual_rms,
            registration_iteration->radial_residual_rms);
    }
    if (field_stats != NULL) {
        fprintf(fp,
            ",\n  \"field_solve\": {\"observations\": %zu, "
            "\"observation_rms\": %.17g, "
            "\"observation_max\": %.17g, "
            "\"gauge_consistency_rms\": %.17g, "
            "\"edge_gradient_delta_rms\": %.17g, "
            "\"irls_rounds\": %d, "
            "\"robust_energy_initial\": %.17g, "
            "\"robust_energy_final\": %.17g, "
            "\"condition_min_equilibrated\": %.17g, "
            "\"condition_max_equilibrated\": %.17g, "
            "\"local_residual_rms_before\": %.17g, "
            "\"local_residual_rms_after\": %.17g}",
            field_stats->observations, field_stats->observation_rms,
            field_stats->observation_max,
            field_stats->gauge_consistency_rms,
            field_stats->edge_gradient_delta_rms,
            field_stats->gauge_solver.irls_rounds_run,
            field_stats->gauge_solver.robust_energy_initial,
            field_stats->gauge_solver.robust_energy_final,
            field_stats->gauge_solver.condition_min,
            field_stats->gauge_solver.condition_max,
            field_stats->gauge_solver.local_residual_rms_before,
            field_stats->gauge_solver.local_residual_rms_after);
    }
    fprintf(fp, "\n}\n");
    if (fclose(fp) != 0 ||
        append_pipeline_trace_index(cfg, frame, obj_name, json_name, "mesh",
                                    stage) != 0)
        return -1;
    (*cfg->pipeline_trace_next_frame)++;
    fprintf(stderr,
        "[atlas_strip_scroll] TRACE %04d %-32s mesh span=%.6g "
        "flips=%zu deg=%zu aspect(p99=%.4g max=%.4g over50=%zu) "
        "stretch(p95=%.4g max=%.4g over4=%zu)\n",
        frame, stage, q.finite_vertices ? q.u_max - q.u_min : (double)NAN,
        q.relative_flips, q.degenerate_faces, q.aspect_p99, q.aspect_max,
        q.aspect_over_50, q.stretch_p95, q.stretch_max,
        q.stretch_over_4);
    return 0;
}

static int write_pipeline_trace_strokes(
    const ScrollConfig *cfg, const AtlasCandidateSet *set,
    const char *stage, const double *u,
    const AtlasStripRobustTraceEntry *robust_iteration)
{
    if (cfg == NULL || cfg->pipeline_trace_dir == NULL ||
        cfg->pipeline_trace_next_frame == NULL)
        return 0;
    if (set == NULL || stage == NULL || u == NULL) return -1;
    int frame = *cfg->pipeline_trace_next_frame;
    char obj_name[256], json_name[256];
    int no = snprintf(obj_name, sizeof(obj_name), "%04d_%s.obj", frame, stage);
    int nj = snprintf(json_name, sizeof(json_name), "%04d_%s.json", frame,
                      stage);
    if (no < 0 || (size_t)no >= sizeof(obj_name) ||
        nj < 0 || (size_t)nj >= sizeof(json_name) ||
        write_strokes_obj(cfg->pipeline_trace_dir, obj_name, set, u, 1) != 0)
        return -1;
    size_t ns = set->problem.nsamples;
    double *ratio = (double *)malloc((ns ? ns : 1) * sizeof(*ratio));
    if (ratio == NULL) return -1;
    size_t finite = 0, nonfinite = 0, nratio = 0, reversed = 0;
    double umin = DBL_MAX, umax = -DBL_MAX;
    double delta2 = 0.0, delta_max = 0.0;
    size_t ndelta = 0;
    for (size_t i = 0; i < ns; i++) {
        if (!isfinite(u[i])) {
            nonfinite++;
            continue;
        }
        finite++;
        if (u[i] < umin) umin = u[i];
        if (u[i] > umax) umax = u[i];
        if (set->initial_u != NULL && isfinite(set->initial_u[i])) {
            double delta = u[i] - set->initial_u[i];
            delta2 += delta * delta;
            if (fabs(delta) > delta_max) delta_max = fabs(delta);
            ndelta++;
        }
    }
    for (size_t s = 0; s < set->problem.nstrokes; s++) {
        const AtlasStripStroke *stroke = &set->strokes[s];
        for (int32_t j = 1; j < stroke->count; j++) {
            size_t a = (size_t)stroke->first + (size_t)j - 1;
            size_t b = a + 1;
            double ds = set->samples[b].s - set->samples[a].s;
            double du = u[b] - u[a];
            if (!isfinite(ds) || !isfinite(du) || fabs(ds) <= 1e-12)
                continue;
            ratio[nratio++] = fabs(du / ds);
            if (du * ds < 0.0) reversed++;
        }
    }
    qsort(ratio, nratio, sizeof(*ratio), compare_double_value);
    double p50 = nratio ? ratio[(size_t)(0.50 * (double)(nratio - 1))] : 0.0;
    double p95 = nratio ? ratio[(size_t)(0.95 * (double)(nratio - 1))] : 0.0;
    double p99 = nratio ? ratio[(size_t)(0.99 * (double)(nratio - 1))] : 0.0;
    double rmax = nratio ? ratio[nratio - 1] : 0.0;
    free(ratio);
    FILE *fp = as_open(cfg->pipeline_trace_dir, json_name, "wb");
    if (fp == NULL) return -1;
    char umin_text[64], umax_text[64], uspan_text[64];
    as_json_number(umin_text, finite ? umin : (double)NAN);
    as_json_number(umax_text, finite ? umax : (double)NAN);
    as_json_number(uspan_text, finite ? umax - umin : (double)NAN);
    fprintf(fp, "{\n  \"frame\": %d,\n  \"domain\": \"strokes\",\n"
                "  \"stage\": ", frame);
    as_write_json_string(fp, stage);
    fprintf(fp, ",\n  \"source_output_dir\": ");
    as_write_json_string(fp, cfg->out_dir);
    fprintf(fp,
        ",\n  \"representation\": \"ordered coarse samples in (u, axial, 0)\",\n"
        "  \"samples\": %zu,\n  \"strokes\": %zu,\n"
        "  \"support_components\": %zu,\n"
        "  \"finite_samples\": %zu,\n  \"nonfinite_samples\": %zu,\n"
        "  \"u_min\": %s,\n  \"u_max\": %s,\n  \"u_span\": %s,\n"
        "  \"delta_from_candidate_initial_rms\": %.17g,\n"
        "  \"delta_from_candidate_initial_max\": %.17g,\n"
        "  \"stroke_segments\": {\"count\": %zu, \"reversed\": %zu, "
        "\"abs_du_over_ds_p50\": %.17g, \"abs_du_over_ds_p95\": %.17g, "
        "\"abs_du_over_ds_p99\": %.17g, \"abs_du_over_ds_max\": %.17g}",
        ns, set->problem.nstrokes, set->stats.support_components,
        finite, nonfinite, umin_text, umax_text, uspan_text,
        ndelta ? sqrt(delta2 / (double)ndelta) : 0.0, delta_max,
        nratio, reversed, p50, p95, p99, rmax);
    if (robust_iteration != NULL) {
        fprintf(fp,
            ",\n  \"coarse_robust_solve\": {\"iteration\": %d, "
            "\"phase\": %d, \"qp_rc\": %d, \"objective\": %.17g, "
            "\"max_u_change\": %.17g, "
            "\"max_membership_change\": %.17g, "
            "\"residual_rms\": %.17g, \"max_residual\": %.17g, "
            "\"membership_min\": %.17g, "
            "\"membership_mean\": %.17g, "
            "\"membership_max\": %.17g, "
            "\"downweighted_members\": %zu}",
            robust_iteration->iteration, robust_iteration->phase,
            robust_iteration->qp_rc, robust_iteration->objective,
            robust_iteration->max_u_change,
            robust_iteration->max_membership_change,
            robust_iteration->residual_rms, robust_iteration->max_residual,
            robust_iteration->membership_min,
            robust_iteration->membership_mean,
            robust_iteration->membership_max,
            robust_iteration->downweighted_members);
    }
    fprintf(fp, "\n}\n");
    if (fclose(fp) != 0 ||
        append_pipeline_trace_index(cfg, frame, obj_name, json_name,
                                    "strokes", stage) != 0)
        return -1;
    (*cfg->pipeline_trace_next_frame)++;
    fprintf(stderr,
        "[atlas_strip_scroll] TRACE %04d %-32s strokes span=%.6g "
        "reversed=%zu speed(p95=%.4g max=%.4g)\n",
        frame, stage, finite ? umax - umin : (double)NAN,
        reversed, p95, rmax);
    return 0;
}

static int write_field_summary(Arena_T arena, const ScrollConfig *cfg,
                               const PieceSet *ps,
                               const ScaffoldCalib *cal,
                               const float *raw_u, const double *field_u,
                               const AtlasFieldApplyStats *stats)
{
    double umin = DBL_MAX, umax = -DBL_MAX;
    size_t nedge = 3 * ps->nf;
    double *stretch = (double *)ARENA_ALLOC(
        arena, (nedge ? nedge : 1) * sizeof(double));
    size_t ns = 0, over4 = 0, flip = 0, comparable = 0, degenerate = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        if (field_u[i] < umin) umin = field_u[i];
        if (field_u[i] > umax) umax = field_u[i];
    }
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t v[3] = {ps->faces[f * 3 + 0], ps->faces[f * 3 + 1],
                        ps->faces[f * 3 + 2]};
        double axial[3];
        for (int k = 0; k < 3; k++)
            axial[k] = axis_coordinate(cal, &ps->verts[(size_t)v[k] * 3]);
        double raw_det = ((double)raw_u[v[1]] - (double)raw_u[v[0]]) *
                         (axial[2] - axial[0]) -
                         ((double)raw_u[v[2]] - (double)raw_u[v[0]]) *
                         (axial[1] - axial[0]);
        double new_det = (field_u[v[1]] - field_u[v[0]]) *
                         (axial[2] - axial[0]) -
                         (field_u[v[2]] - field_u[v[0]]) *
                         (axial[1] - axial[0]);
        if (fabs(new_det) < 1e-10) degenerate++;
        if (fabs(raw_det) >= 1e-10) {
            comparable++;
            if (raw_det * new_det < 0.0) flip++;
        }
        for (int e = 0; e < 3; e++) {
            int32_t a = v[e], b = v[(e + 1) % 3];
            const float *pa = &ps->verts[(size_t)a * 3];
            const float *pb = &ps->verts[(size_t)b * 3];
            double dz = (double)pa[0] - (double)pb[0];
            double dy = (double)pa[1] - (double)pb[1];
            double dx = (double)pa[2] - (double)pb[2];
            double length = sqrt(dz * dz + dy * dy + dx * dx);
            if (length <= 1e-12) continue;
            double du = field_u[a] - field_u[b];
            double dv = axial[e] - axial[(e + 1) % 3];
            double ratio = sqrt(du * du + dv * dv) / length;
            stretch[ns++] = ratio;
            if (ratio > 4.0) over4++;
        }
    }
    qsort(stretch, ns, sizeof(double), compare_double_value);
    double p50 = ns ? stretch[(size_t)(0.50 * (double)(ns - 1))] : 0.0;
    double p95 = ns ? stretch[(size_t)(0.95 * (double)(ns - 1))] : 0.0;
    double p99 = ns ? stretch[(size_t)(0.99 * (double)(ns - 1))] : 0.0;
    double smax = ns ? stretch[ns - 1] : 0.0;

    FILE *fp = as_open(cfg->out_dir, "field_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n  \"operation\": \"rigid_support_chart_registration\",\n"
        "  \"mesh_components\": %zu, \"observed_components\": %zu, "
        "\"anchored_unobserved_components\": %zu,\n"
        "  \"observations\": %zu, \"observation_rms\": %.17g, "
        "\"observation_max\": %.17g,\n"
        "  \"gauge_reconciliation\": {\"graph_components\": %zu, "
        "\"correction_rms\": %.17g, \"correction_max\": %.17g, "
        "\"consistency_rms\": %.17g, \"reference_shift_rms\": %.17g, "
        "\"reference_shift_max\": %.17g, \"reference_residual_rms\": %.17g},\n"
        "  \"edge_gradient_delta_rms\": %.17g, "
        "\"skipped_degenerate_triangles\": %zu,\n"
        "  \"gauge_solver\": {\"solver\": "
        "\"grounded_graph_laplacian_jacobi_sparse_spd_cholesky\", "
        "\"variables\": %zu, \"local_edges\": %zu, "
        "\"irls_rounds\": %d, \"downweighted_edges\": %zu, "
        "\"robust_energy_initial\": %.17g, "
        "\"robust_energy_final\": %.17g, "
        "\"condition_min_equilibrated\": %.17g, "
        "\"condition_max_equilibrated\": %.17g, "
        "\"raw_diagonal_min\": %.17g, "
        "\"raw_diagonal_max\": %.17g, "
        "\"local_residual_rms_before\": %.17g, "
        "\"local_residual_rms_after\": %.17g},\n"
        "  \"atlas\": {\"u_min\": %.17g, \"u_max\": %.17g, "
        "\"u_span\": %.17g, \"relative_flips\": %zu, "
        "\"comparable_faces\": %zu, \"degenerate_faces\": %zu, "
        "\"stretch_p50\": %.17g, \"stretch_p95\": %.17g, "
        "\"stretch_p99\": %.17g, \"stretch_max\": %.17g, "
        "\"stretch_over_4\": %zu}\n}\n",
        stats->mesh_components, stats->observed_components,
        stats->anchored_unobserved_components,
        stats->observations, stats->observation_rms,
        stats->observation_max, stats->gauge_graph_components,
        stats->gauge_correction_rms, stats->gauge_correction_max,
        stats->gauge_consistency_rms, stats->reference_gauge_shift_rms,
        stats->reference_gauge_shift_max, stats->reference_residual_rms,
        stats->edge_gradient_delta_rms,
        stats->skipped_degenerate_triangles,
        stats->gauge_solver.variables, stats->gauge_solver.local_edges,
        stats->gauge_solver.irls_rounds_run,
        stats->gauge_solver.downweighted_edges,
        stats->gauge_solver.robust_energy_initial,
        stats->gauge_solver.robust_energy_final,
        stats->gauge_solver.condition_min,
        stats->gauge_solver.condition_max,
        stats->gauge_solver.raw_diagonal_min,
        stats->gauge_solver.raw_diagonal_max,
        stats->gauge_solver.local_residual_rms_before,
        stats->gauge_solver.local_residual_rms_after,
        umin, umax, umax - umin, flip, comparable, degenerate,
        p50, p95, p99, smax, over4);
    fclose(fp);
    return 0;
}

typedef struct {
    const ScrollConfig *cfg;
    const PieceSet *ps;
    const ScaffoldCalib *cal;
    const AtlasCandidateSet *set;
} AsFieldPipelineTrace;

static int dump_field_registration_iteration(
    void *opaque, const double *corrected_sample, const double *field_u,
    const AtlasRegisterIterationStats *iteration)
{
    AsFieldPipelineTrace *trace = (AsFieldPipelineTrace *)opaque;
    if (trace == NULL || trace->cfg == NULL || trace->ps == NULL ||
        trace->cal == NULL || trace->set == NULL || corrected_sample == NULL ||
        field_u == NULL || iteration == NULL)
        return -1;
    char target_stage[128], mesh_stage[128];
    int nt = snprintf(target_stage, sizeof(target_stage),
                      "field_gauge_irls_%02d_target", iteration->iteration);
    int nm = snprintf(mesh_stage, sizeof(mesh_stage),
                      "field_gauge_irls_%02d_mesh", iteration->iteration);
    if (nt < 0 || (size_t)nt >= sizeof(target_stage) ||
        nm < 0 || (size_t)nm >= sizeof(mesh_stage) ||
        write_pipeline_trace_strokes(trace->cfg, trace->set, target_stage,
                                     corrected_sample, NULL) != 0 ||
        write_pipeline_trace_mesh(trace->cfg, trace->ps, trace->cal,
                                  mesh_stage, field_u, iteration, NULL) != 0)
        return -1;
    return 0;
}

static int run_field_stage(Arena_T arena, const ScrollConfig *cfg,
                           const PieceSet *ps, const ScaffoldCalib *cal,
                           const float *raw_u, const float *reference_u,
                           const AtlasCandidateSet *set,
                           const double *sample_target,
                           double **out_corrected_target,
                           double **out_field_u)
{
    double *field_u = (double *)ARENA_ALLOC(
        arena, ps->nv * sizeof(double));
    double *corrected_target = (double *)ARENA_ALLOC(
        arena, set->problem.nsamples * sizeof(double));
    int32_t *sample_component = (int32_t *)ARENA_ALLOC(
        arena, set->problem.nsamples * sizeof(int32_t));
    for (size_t i = 0; i < set->problem.nsamples; i++) {
        int32_t stroke = set->samples[i].stroke;
        sample_component[i] = set->strokes[stroke].component;
    }
    AtlasFieldApplyStats stats;
    memset(&stats, 0, sizeof(stats));
    if (write_pipeline_trace_strokes(cfg, set, "field_apply_input_target",
                                     sample_target, NULL) != 0)
        return -1;
    AsFieldPipelineTrace pipeline_trace = {cfg, ps, cal, set};
    double t0 = ves_clock_sec();
    int rc = AtlasFieldApply_solve(
        arena, ps->verts, ps->nv, ps->faces, ps->nf, raw_u, reference_u,
        set->sample_ref, sample_target, sample_component,
        set->problem.nsamples, set->stats.support_components,
        dump_field_registration_iteration, &pipeline_trace,
        corrected_target, field_u, &stats);
    fprintf(stderr,
        "[atlas_strip_scroll] rigid field registration %.2fs: rc=%d "
        "components=%zu observed=%zu anchors=%zu gauge_rms=%.4g "
        "consistency=%.4g local_rms=%.4g->%.4g condition=%.4g..%.4g "
        "energy=%.6g->%.6g grad_delta=%.3e\n",
        ves_clock_sec() - t0, rc, stats.mesh_components,
        stats.observed_components, stats.anchored_unobserved_components,
        stats.gauge_correction_rms, stats.gauge_consistency_rms,
        stats.gauge_solver.local_residual_rms_before,
        stats.gauge_solver.local_residual_rms_after,
        stats.gauge_solver.condition_min, stats.gauge_solver.condition_max,
        stats.gauge_solver.robust_energy_initial,
        stats.gauge_solver.robust_energy_final,
        stats.edge_gradient_delta_rms);
    if (rc != 0) return -1;
    if (write_pipeline_trace_strokes(cfg, set, "field_gauge_corrected_target",
                                     corrected_target, NULL) != 0 ||
        write_pipeline_trace_mesh(cfg, ps, cal, "field_rigid_output",
                                  field_u, NULL, &stats) != 0)
        return -1;
    if (write_global_mesh_obj(cfg->out_dir, "atlas_flat.obj", ps, cal,
                              field_u, 1) != 0 ||
        write_global_mesh_obj(cfg->out_dir, "atlas_world_uv.obj", ps, cal,
                              field_u, 0) != 0 ||
        export_field_placed(cfg, ps, cal, field_u) != 0 ||
        write_field_summary(arena, cfg, ps, cal, raw_u, field_u,
                            &stats) != 0)
        return -1;
    fprintf(stderr,
        "[atlas_strip_scroll] flattened mesh and placed sidecars written to %s\n",
        cfg->out_dir);
    if (out_corrected_target != NULL)
        *out_corrected_target = corrected_target;
    if (out_field_u != NULL) *out_field_u = field_u;
    return 0;
}

typedef struct {
    double *sample_final;
    double *member_membership;
    double *corrected_sample;
    double *field_u;
    size_t fixed_membership_rows_preserved;
    size_t fixed_membership_rows_above_floor;
    size_t fixed_membership_rows_zeroed;
    int used_fixed_membership;
    AtlasStripMetrics final_metrics;
} AtlasStageResult;

typedef struct {
    const ScrollConfig *cfg;
    const AtlasCandidateSet *set;
} AsStripPipelineTrace;

static int dump_strip_robust_iteration(
    void *context, const AtlasStripProblem *problem, const double *x,
    const AtlasStripMember *members,
    const AtlasStripRobustTraceEntry *iteration)
{
    (void)members;
    AsStripPipelineTrace *trace = (AsStripPipelineTrace *)context;
    if (trace == NULL || trace->cfg == NULL || trace->set == NULL ||
        problem == NULL || x == NULL || iteration == NULL ||
        problem->nsamples != trace->set->problem.nsamples)
        return -1;
    const char *phase = iteration->phase == ATLAS_STRIP_ROBUST_INITIAL
                      ? "initial"
                      : (iteration->phase == ATLAS_STRIP_ROBUST_L1
                         ? "l1" : "likelihood");
    char stage[128];
    int n = snprintf(stage, sizeof(stage), "coarse_robust_%02d_%s",
                     iteration->iteration, phase);
    if (n < 0 || (size_t)n >= sizeof(stage)) return -1;
    return write_pipeline_trace_strokes(trace->cfg, trace->set, stage, x,
                                        iteration);
}

static int run_atlas_stage(Arena_T arena, const ScrollConfig *cfg,
                           const PieceSet *ps, const ScaffoldCalib *cal,
                           const float *raw_u, const float *reference_u,
                           AtlasCandidateSet *set,
                           const double *initial_sample,
                           const double *fixed_membership,
                           const double *frozen_diagnostic,
                           AtlasStageResult *result)
{
    if (result == NULL) return -1;
    memset(result, 0, sizeof *result);
    if (write_strokes_obj(cfg->out_dir, "strokes_world.obj", set,
                          set->initial_u, 0) != 0 ||
        write_strokes_obj(cfg->out_dir, "initial_parameter.obj", set,
                          set->initial_u, 1) != 0 ||
        (initial_sample != NULL &&
         write_strokes_obj(cfg->out_dir, "solver_start_parameter.obj", set,
                           initial_sample, 1) != 0) ||
        write_continuations_obj(cfg->out_dir, set) != 0 ||
        write_candidates_obj(cfg->out_dir, "candidate_links_world.obj",
                             set, NULL) != 0 ||
        write_samples_csv(cfg->out_dir, ps, set, NULL, NULL) != 0 ||
        write_constraints_csv(cfg->out_dir, set, NULL, NULL, NULL) != 0 ||
        write_selection_csv(cfg->out_dir, set, frozen_diagnostic) != 0 ||
        write_bundles_csv(cfg->out_dir, set) != 0 ||
        write_bundle_decisions_obj(cfg->out_dir, set) != 0 ||
        write_summary(cfg, ps, cal, set, 0, 0,
                      NULL, NULL, NULL, NULL,
                      fixed_membership != NULL) != 0)
        return -1;
    if (cfg->build_only) return 0;

    size_t nvar = set->problem.nsamples + set->problem.ncross_sections;
    double *relaxed = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    memcpy(relaxed,
           initial_sample != NULL ? initial_sample : set->initial_u,
           set->problem.nsamples * sizeof(double));
    AtlasStrip_initialize_intercepts(&set->problem, &cfg->strip, relaxed);
    if (write_pipeline_trace_strokes(cfg, set, "coarse_solver_start",
                                     relaxed, NULL) != 0)
        return -1;
    double *membership = (double *)ARENA_ALLOC(
        arena, set->problem.nmembers * sizeof(double));
    size_t trace_cap = (size_t)(1 + cfg->robust.l1_iterations +
                                cfg->robust.likelihood_iterations);
    AtlasStripRobustTrace trace;
    trace.entry = (AtlasStripRobustTraceEntry *)ARENA_ALLOC(
        arena, (trace_cap ? trace_cap : 1) * sizeof(*trace.entry));
    trace.capacity = trace_cap;
    trace.count = 0;
    AsStripPipelineTrace pipeline_trace;
    pipeline_trace.cfg = cfg;
    pipeline_trace.set = set;
    trace.iteration_fn = dump_strip_robust_iteration;
    trace.iteration_context = &pipeline_trace;
    AtlasStripRobustStats robust_stats;
    memset(&robust_stats, 0, sizeof robust_stats);
    double t0 = ves_clock_sec();
    int robust_rc = 0;
    if (fixed_membership != NULL) {
        size_t preserved = 0, above_floor = 0, zeroed = 0;
        for (size_t i = 0; i < set->problem.nmembers; i++) {
            if (!isfinite(fixed_membership[i]) || fixed_membership[i] < 0.0)
                return -1;
            if (set->member_state[i] == ATLAS_CANDIDATE_INACTIVE) {
                membership[i] = 0.0;
                if (fixed_membership[i] > 0.0) zeroed++;
            } else {
                membership[i] = fixed_membership[i];
                preserved++;
                if (membership[i] > cfg->strip.membership_floor)
                    above_floor++;
            }
        }
        fprintf(stderr,
            "[atlas_strip_scroll] %s fixed HEAD memberships: "
            "preserved=%zu (%zu above floor) zeroed=%zu "
            "(no IRLS restart)\n",
            cfg->out_dir, preserved, above_floor, zeroed);
        result->fixed_membership_rows_preserved = preserved;
        result->fixed_membership_rows_above_floor = above_floor;
        result->fixed_membership_rows_zeroed = zeroed;
        result->used_fixed_membership = 1;
    } else {
        robust_rc = AtlasStrip_solve_robust(
            arena, &set->problem, &cfg->strip, &cfg->robust, &cfg->qp,
            relaxed, membership, &trace, &robust_stats);
        fprintf(stderr,
            "[atlas_strip_scroll] %s robust %.2fs: rc=%d solves=%d "
            "objective=%.9g sigma=%.6g downweighted=%zu\n",
            cfg->out_dir, ves_clock_sec() - t0, robust_rc,
            robust_stats.total_qp_solves,
            robust_stats.final_qp.objective_final,
            robust_stats.likelihood_sigma, robust_stats.downweighted_members);
    }
    if (robust_rc != 0) {
        write_robust_trace(cfg->out_dir, &trace);
        write_summary(cfg, ps, cal, set, robust_rc, -99,
                      &robust_stats, NULL, NULL, NULL,
                      fixed_membership != NULL);
        return -1;
    }

    AtlasStripMember *final_members = (AtlasStripMember *)ARENA_ALLOC(
        arena, set->problem.nmembers * sizeof(AtlasStripMember));
    memcpy(final_members, set->problem.members,
           set->problem.nmembers * sizeof(AtlasStripMember));
    for (size_t i = 0; i < set->problem.nmembers; i++)
        final_members[i].membership = membership[i];
    AtlasStripProblem final_problem = set->problem;
    final_problem.members = final_members;
    AtlasStripOptions final_opts = cfg->strip;
    final_opts.mode = ATLAS_STRIP_FINAL;
    double *final = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    memcpy(final, relaxed, nvar * sizeof(double));
    AtlasStrip_initialize_intercepts(&final_problem, &final_opts, final);
    AtlasStripSystem final_system;
    int final_build_rc = AtlasStrip_build(
        arena, &final_problem, &final_opts, &final_system);
    MonotoneQpStats final_stats;
    memset(&final_stats, 0, sizeof final_stats);
    t0 = ves_clock_sec();
    int final_rc = -1;
    if (final_build_rc == 0 && cfg->robust.use_active_set) {
        final_rc = MonotoneQp_solve(arena, &final_system.qp, &cfg->qp,
                                    final, NULL, &final_stats);
    } else if (final_build_rc == 0) {
        double lr = 0.0;
        final_rc = MonotoneQp_solve_ls(arena, &final_system.qp, final, &lr);
        if (final_rc == 0) {
            AtlasStrip_project_monotone(arena, &final_problem, &final_opts,
                                        final);
            AtlasStrip_initialize_intercepts(&final_problem, &final_opts,
                                             final);
        }
        final_stats.objective_final =
            MonotoneQp_objective(&final_system.qp, final);
        final_stats.max_reduced_linear_residual = lr;
    }
    fprintf(stderr,
        "[atlas_strip_scroll] %s final %.2fs: rc=%d iterations=%d "
        "objective=%.9g min_slack=%.3e stationarity=%.3e\n",
        cfg->out_dir, ves_clock_sec() - t0, final_rc,
        final_stats.iterations, final_stats.objective_final,
        final_stats.min_slack, final_stats.stationarity_residual);

    AtlasStripMetrics relaxed_metrics, final_metrics;
    memset(&relaxed_metrics, 0, sizeof relaxed_metrics);
    memset(&final_metrics, 0, sizeof final_metrics);
    AtlasStrip_measure(&final_problem, &cfg->strip, relaxed,
                       &relaxed_metrics);
    if (final_rc == 0)
        AtlasStrip_measure(&final_problem, &final_opts, final,
                           &final_metrics);
    write_robust_trace(cfg->out_dir, &trace);
    write_strokes_obj(cfg->out_dir, "relaxed_parameter.obj", set,
                      relaxed, 1);
    if (final_rc == 0)
        write_strokes_obj(cfg->out_dir, "final_parameter.obj", set,
                          final, 1);
    if (write_pipeline_trace_strokes(cfg, set, "coarse_relaxed_output",
                                     relaxed, NULL) != 0 ||
        (final_rc == 0 &&
         write_pipeline_trace_strokes(cfg, set, "coarse_final_output",
                                      final, NULL) != 0))
        return -1;
    write_candidates_obj(cfg->out_dir, "candidate_likelihoods_world.obj",
                         set, membership);
    write_samples_csv(cfg->out_dir, ps, set, relaxed,
                      final_rc == 0 ? final : NULL);
    write_constraints_csv(cfg->out_dir, set, membership, relaxed,
                          final_rc == 0 ? final : NULL);
    write_summary(cfg, ps, cal, set, robust_rc, final_rc,
                  &robust_stats, &final_stats, &relaxed_metrics,
                  final_rc == 0 ? &final_metrics : NULL,
                  fixed_membership != NULL);
    if (final_rc != 0) return -1;

    double *corrected = NULL, *field_u = NULL;
    int field_rc = run_field_stage(
        arena, cfg, ps, cal, raw_u, reference_u, set, final,
        &corrected, &field_u);
    fprintf(stderr,
        "[atlas_strip_scroll] %s %s; final rms length=%.4g align=%.4g "
        "continuation=%.4g min-ratio=%.6g\n",
        cfg->out_dir, field_rc == 0 ? "PASSED" : "FAILED",
        final_metrics.rms_length_row, final_metrics.rms_align_row,
        final_metrics.rms_continuation, final_metrics.min_monotone_ratio);
    if (field_rc != 0) return -1;
    result->sample_final = final;
    result->member_membership = membership;
    result->corrected_sample = corrected;
    result->field_u = field_u;
    result->final_metrics = final_metrics;
    return 0;
}

/* ==========================================================================
 * Radially-ordered warp.
 *
 * Everything upstream of here is unchanged HEAD machinery.  This stage takes
 * the composed vertex field, re-solves the coarse atlas with radial ordering
 * added, extends the result back to every vertex, and reports whether the
 * exchange was worth it -- overlap removed versus physical weld disturbed.
 * The weld oracle is scored but never fed back: it is held out precisely so
 * it can reject this stage the way it rejected rigid chart layout.
 * ======================================================================== */

/* Sample the vertex field at a coarse sample's exact edge position, so the
 * warp starts from the composed field rather than re-deriving it. */
static double warp_sample_from_field(const AtlasCandidateSet *set,
                                     const double *field, size_t i)
{
    const AtlasCandidateSampleRef *ref = &set->sample_ref[i];
    double a = field[ref->mesh_vertex[0]];
    double b = field[ref->mesh_vertex[1]];
    return a + (b - a) * ref->mesh_t;
}

typedef struct {
    size_t total;
    size_t whitelisted;
    size_t unwhitelisted;
    size_t cross_component;
    size_t same_component;
    size_t turn_nonzero;
} WarpOverlapCount;

/* Count overlaps the way the gate table reads them: delamination is a physical
 * fact about the scroll, not a defect, so it is excluded from the number the
 * stage is judged on -- but it is always reported alongside. */
static void warp_count_overlaps(const AtlasOverlapAudit *audit,
                                const uint64_t *delam_keys, size_t ndelam,
                                WarpOverlapCount *out)
{
    memset(out, 0, sizeof *out);
    out->total = audit->npairs;
    for (size_t i = 0; i < audit->npairs; i++) {
        const AtlasOverlapPair *p = &audit->pairs[i];
        uint64_t key = boxcut_face_pair_key(p->face0, p->face1);
        int allowed = ndelam > 0 &&
                      bsearch(&key, delam_keys, ndelam, sizeof(uint64_t),
                              boxcut_compare_u64) != NULL;
        if (allowed) out->whitelisted++;
        else out->unwhitelisted++;
        if (p->component0 == p->component1) out->same_component++;
        else out->cross_component++;
        if (p->turn != 0) out->turn_nonzero++;
    }
}

typedef struct {
    WarpOverlapCount before;
    WarpOverlapCount after;
    AtlasXyzWeldAuditStats weld_before;
    AtlasXyzWeldAuditStats weld_after;
    int weld_rc;
    AtlasWarpStats warp;
    AtlasRadialOrderStats order;
    double seam_rms_before;
    double seam_rms_after;
    size_t relative_flips;
} WarpStageReport;

/* Field-level orientation check: a triangle whose UV winding flips relative to
 * the incoming field means the warp folded the atlas over itself. */
static size_t warp_count_relative_flips(const PieceSet *ps,
                                        const double *before,
                                        const double *after,
                                        const double *axial)
{
    size_t flips = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = ps->faces[f * 3 + 0];
        int32_t b = ps->faces[f * 3 + 1];
        int32_t c = ps->faces[f * 3 + 2];
        double b0 = (before[b] - before[a]) * (axial[c] - axial[a]) -
                    (before[c] - before[a]) * (axial[b] - axial[a]);
        double a0 = (after[b] - after[a]) * (axial[c] - axial[a]) -
                    (after[c] - after[a]) * (axial[b] - axial[a]);
        if (b0 != 0.0 && a0 != 0.0 && ((b0 > 0.0) != (a0 > 0.0))) flips++;
    }
    return flips;
}

static double warp_seam_rms(Arena_T arena, const PieceSet *ps,
                            const AtlasCandidateSet *set, const double *field)
{
    AtlasSeamAudit seam;
    if (AtlasSeamAudit_build(arena, ps, field, set->vertex_mesh_component,
                             set->mesh_components, &seam) != 0)
        return -1.0;
    double sum = 0.0;
    size_t n = 0;
    for (size_t i = 0; i < seam.npairs; i++) {
        int32_t b = seam.pairs[i].bundle;
        if (b < 0 || (size_t)b >= seam.nbundles) continue;
        if (!smooth_seam_bundle_eligible(&seam.bundles[b])) continue;
        double d = field[seam.pairs[i].vertex1] - field[seam.pairs[i].vertex0];
        sum += d * d;
        n++;
    }
    return n > 0 ? sqrt(sum / (double)n) : 0.0;
}

/* One frozen XYZ-derived weld topology scored against two UV fields.  The
 * topology is built from geometry alone, so neither field can influence which
 * bridge faces exist -- only how far they have to stretch to fit. */
static int warp_weld_oracle(const PieceSet *ps, const ScaffoldCalib *cal,
                            const double *before, const double *after,
                            const double *axial,
                            AtlasXyzWeldAuditStats *out_before,
                            AtlasXyzWeldAuditStats *out_after)
{
    memset(out_before, 0, sizeof *out_before);
    memset(out_after, 0, sizeof *out_after);
    Arena_T topo = Arena_new();
    BpaBridgeGate gate;
    memset(&gate, 0, sizeof gate);
    gate.umb_y = cal->axis_point[1];
    gate.umb_x = cal->axis_point[2];
    gate.pitch = cal->pitch;
    gate.tol = SEAM_WIND_TOL_DEFAULT_TURNS;
    gate.hard = SEAM_WIND_HARD_TOL_DEFAULT_TURNS;
    AtlasXyzWeldTopology topology;
    memset(&topology, 0, sizeof topology);
    int rc = -1;
    if (AtlasXyzWeldTopology_build(topo, ps, 128.0f, 1.5f, 0.0f, 6.0f,
                                   &gate, &topology) == 0 &&
        topology.nfaces > 0) {
        /* One block, zero shift: the warp moves vertices, not charts, so the
         * oracle scores the field itself. */
        int32_t *one_label = (int32_t *)ARENA_CALLOC(topo, ps->nf ? ps->nf : 1,
                                                     sizeof(int32_t));
        double zero = 0.0;
        Arena_T a0 = Arena_new();
        int rc0 = AtlasXyzWeldAudit_evaluate(a0, ps, &topology, before, axial,
                                             one_label, 1, &zero, out_before);
        Arena_dispose(&a0);
        Arena_T a1 = Arena_new();
        int rc1 = AtlasXyzWeldAudit_evaluate(a1, ps, &topology, after, axial,
                                             one_label, 1, &zero, out_after);
        Arena_dispose(&a1);
        rc = rc0 == 0 && rc1 == 0 ? 0 : -1;
    }
    Arena_dispose(&topo);
    return rc;
}

/*
 * Dump the support components that carry the collapse, in both spaces.
 *
 * A support component is one connected piece of the candidate graph, and
 * inside it the strip rows determine u up to a single constant -- so a
 * component whose own wrap pairs sit on top of each other in u cannot be
 * repaired by moving anything.  Seeing it is the point: the XYZ dump shows
 * concentric wraps that are unmistakably separate in space, and the UV dump of
 * the same triangles shows them landing on each other.
 */
/*
 * Colour the dump by CONNECTED COMPONENT IN XYZ -- pieces that are physically
 * separate surfaces, welded to nothing.  The same colour array is used for
 * both the world and the atlas dump, so a piece can be followed from one to
 * the other: in XYZ the colours are nested concentric wraps, and in UV the
 * ones that landed on each other are immediately readable by their colour.
 */
static void collapse_hue_rgb(size_t index, double rgb[3])
{
    /* Golden-angle hue walk keeps neighbouring indices far apart in colour. */
    double h = fmod((double)index * 0.6180339887498949, 1.0) * 6.0;
    double s = 0.78, v = 0.96;
    int i = (int)h;
    double f = h - (double)i;
    double p = v * (1.0 - s);
    double q = v * (1.0 - s * f);
    double t = v * (1.0 - s * (1.0 - f));
    switch (i % 6) {
        case 0: rgb[0] = v; rgb[1] = t; rgb[2] = p; break;
        case 1: rgb[0] = q; rgb[1] = v; rgb[2] = p; break;
        case 2: rgb[0] = p; rgb[1] = v; rgb[2] = t; break;
        case 3: rgb[0] = p; rgb[1] = q; rgb[2] = v; break;
        case 4: rgb[0] = t; rgb[1] = p; rgb[2] = v; break;
        default: rgb[0] = v; rgb[1] = p; rgb[2] = q; break;
    }
}

typedef struct { int32_t root; size_t faces; } CollapseCcEntry;

static int collapse_cc_cmp(const void *va, const void *vb)
{
    const CollapseCcEntry *a = (const CollapseCcEntry *)va;
    const CollapseCcEntry *b = (const CollapseCcEntry *)vb;
    if (a->faces > b->faces) return -1;
    if (a->faces < b->faces) return 1;
    return a->root < b->root ? -1 : (a->root > b->root ? 1 : 0);
}

/* Union-find over the dumped faces; colours are assigned largest piece first
 * so the dominant surface keeps a stable colour between runs. */
static int collapse_component_colors(Arena_T arena, const PieceSet *ps,
                                     const int32_t *vertex_support,
                                     int32_t target, double *rgb,
                                     size_t *out_ncomp)
{
    Arena_Mark mark = Arena_save(arena);
    UnionFind uf = UF_new(arena, (int32_t)ps->nv);
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = ps->faces[f * 3 + 0];
        int32_t b = ps->faces[f * 3 + 1];
        int32_t c = ps->faces[f * 3 + 2];
        if (vertex_support[a] != target || vertex_support[b] != target ||
            vertex_support[c] != target)
            continue;
        uf_union(&uf, a, b);
        uf_union(&uf, b, c);
    }
    size_t *face_count = (size_t *)ARENA_CALLOC(arena, ps->nv, sizeof(size_t));
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = ps->faces[f * 3 + 0];
        if (vertex_support[a] != target ||
            vertex_support[ps->faces[f * 3 + 1]] != target ||
            vertex_support[ps->faces[f * 3 + 2]] != target)
            continue;
        face_count[uf_find(&uf, a)]++;
    }
    CollapseCcEntry *entry =
        (CollapseCcEntry *)ARENA_ALLOC(arena, ps->nv * sizeof *entry);
    size_t ncomp = 0;
    for (size_t v = 0; v < ps->nv; v++)
        if (face_count[v] > 0) {
            entry[ncomp].root = (int32_t)v;
            entry[ncomp].faces = face_count[v];
            ncomp++;
        }
    qsort(entry, ncomp, sizeof *entry, collapse_cc_cmp);
    int32_t *rank_of = (int32_t *)ARENA_ALLOC(arena, ps->nv * sizeof(int32_t));
    for (size_t v = 0; v < ps->nv; v++) rank_of[v] = -1;
    for (size_t i = 0; i < ncomp; i++) rank_of[entry[i].root] = (int32_t)i;

    for (size_t v = 0; v < ps->nv; v++) {
        rgb[v * 3 + 0] = rgb[v * 3 + 1] = rgb[v * 3 + 2] = 0.5;
        if (vertex_support[v] != target) continue;
        int32_t r = rank_of[uf_find(&uf, (int32_t)v)];
        if (r < 0) continue;
        collapse_hue_rgb((size_t)r, &rgb[v * 3]);
    }
    if (out_ncomp != NULL) *out_ncomp = ncomp;
    Arena_restore(arena, mark);
    return 0;
}

static int write_component_obj(const char *dir, const char *name,
                               const PieceSet *ps, const double *u,
                               const double *axial,
                               const int32_t *vertex_support,
                               const double *rgb,
                               int32_t target, int flat, size_t *out_faces)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    int32_t *remap = (int32_t *)malloc(ps->nv * sizeof(int32_t));
    if (remap == NULL) { fclose(fp); return -1; }
    for (size_t i = 0; i < ps->nv; i++) remap[i] = -1;

    int32_t next = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        if (vertex_support[i] != target) continue;
        remap[i] = ++next;
        if (flat)
            fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                    u[i], axial[i],
                    rgb[i * 3 + 0], rgb[i * 3 + 1], rgb[i * 3 + 2]);
        else
            fprintf(fp, "v %.9g %.9g %.9g %.6f %.6f %.6f\n",
                    (double)ps->verts[i * 3 + 0],
                    (double)ps->verts[i * 3 + 1],
                    (double)ps->verts[i * 3 + 2],
                    rgb[i * 3 + 0], rgb[i * 3 + 1], rgb[i * 3 + 2]);
        fprintf(fp, "vt %.9g %.9g\n", u[i], axial[i]);
    }
    size_t nf = 0;
    for (size_t f = 0; f < ps->nf; f++) {
        int32_t a = remap[ps->faces[f * 3 + 0]];
        int32_t b = remap[ps->faces[f * 3 + 1]];
        int32_t c = remap[ps->faces[f * 3 + 2]];
        if (a < 0 || b < 0 || c < 0) continue;
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
        nf++;
    }
    free(remap);
    fclose(fp);
    if (out_faces != NULL) *out_faces = nf;
    return 0;
}

/* The collapsed pairs themselves, as line segments joining the two samples
 * that should be one circumference apart and are not. */
static int write_collapse_pairs_obj(const char *dir, const char *name,
                                    const PieceSet *ps,
                                    const AtlasCandidateSet *set,
                                    const AtlasRadialOrderSet *order,
                                    const double *sample_u,
                                    const int32_t *vertex_support,
                                    int32_t target, int flat)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    size_t n = 0;
    for (size_t i = 0; i < order->npairs; i++) {
        const AtlasRadialOrderPair *p = &order->pairs[i];
        const AtlasCandidateSampleRef *ri = &set->sample_ref[p->inner];
        const AtlasCandidateSampleRef *ro = &set->sample_ref[p->outer];
        if (vertex_support[ri->mesh_vertex[0]] != target ||
            vertex_support[ro->mesh_vertex[0]] != target)
            continue;
        double t = fabs(p->spacing_target);
        if (!(t > 0.0)) continue;
        double du = sample_u[p->outer] - sample_u[p->inner];
        if (fabs(du) >= 0.5 * t) continue;          /* only the collapsed ones */
        for (int e = 0; e < 2; e++) {
            const AtlasCandidateSampleRef *r = e == 0 ? ri : ro;
            int32_t s = e == 0 ? p->inner : p->outer;
            if (flat) {
                double a0 = (double)ps->verts[r->mesh_vertex[0] * 3];
                double a1 = (double)ps->verts[r->mesh_vertex[1] * 3];
                (void)a0; (void)a1;
                fprintf(fp, "v %.9g %.9g 0\n", sample_u[s], r->axial);
            } else {
                const float *v0 = &ps->verts[r->mesh_vertex[0] * 3];
                const float *v1 = &ps->verts[r->mesh_vertex[1] * 3];
                fprintf(fp, "v %.9g %.9g %.9g\n",
                        (double)v0[0] + ((double)v1[0] - v0[0]) * r->mesh_t,
                        (double)v0[1] + ((double)v1[1] - v0[1]) * r->mesh_t,
                        (double)v0[2] + ((double)v1[2] - v0[2]) * r->mesh_t);
            }
        }
        n++;
    }
    for (size_t i = 0; i < n; i++)
        fprintf(fp, "l %zu %zu\n", 2 * i + 1, 2 * i + 2);
    fclose(fp);
    return 0;
}

static int prune_fused_candidate_links(AtlasCandidateSet *set,
                                       const double *sample_phi,
                                       size_t *out_members,
                                       size_t *out_cross_sections,
                                       size_t *out_continuations);
static int write_collapse_component_objs(Arena_T arena,
                                         const char *dir, const PieceSet *ps,
                                         const ScaffoldCalib *cal,
                                         const AtlasCandidateSet *set,
                                         const AtlasRadialOrderSet *order,
                                         const double *field_u,
                                         const double *sample_u,
                                         const double *axial,
                                         int count)
{
    (void)cal;
    if (count <= 0) return 0;

    /* Support components are a property of strokes; carry them to vertices
     * through the mesh component each sample sits on. */
    size_t nmesh = set->mesh_components;
    int32_t *m2s = (int32_t *)malloc((nmesh ? nmesh : 1) * sizeof(int32_t));
    int32_t *vertex_support =
        (int32_t *)malloc(ps->nv * sizeof(int32_t));
    double *rgb = (double *)malloc(ps->nv * 3 * sizeof(double));
    if (m2s == NULL || vertex_support == NULL || rgb == NULL) {
        free(m2s); free(vertex_support); free(rgb); return -1;
    }
    /*
     * A mesh component can carry strokes from several support components --
     * 24 of 623 do on the 4x5x5 fixture, one of them spanning nine -- so the
     * label is assigned by majority vote over that component's samples rather
     * than by whichever sample happens to be written last.  The ambiguous
     * count is reported: those vertices are attributed to their dominant
     * component and the dump is approximate exactly there.
     */
    size_t nsupport_slots = set->stats.support_components + 1;
    size_t *votes = (size_t *)calloc(nmesh ? nmesh * nsupport_slots : 1,
                                     sizeof(size_t));
    if (votes == NULL) {
        free(m2s); free(vertex_support); free(rgb); return -1;
    }
    for (size_t i = 0; i < set->problem.nsamples; i++) {
        int32_t mc = set->sample_ref[i].mesh_component;
        int32_t sc = set->strokes[set->samples[i].stroke].component;
        if (mc < 0 || (size_t)mc >= nmesh) continue;
        if (sc < 0 || (size_t)sc >= nsupport_slots) continue;
        votes[(size_t)mc * nsupport_slots + (size_t)sc]++;
    }
    size_t ambiguous = 0;
    for (size_t c = 0; c < nmesh; c++) {
        size_t best_n = 0, distinct = 0;
        int32_t best_s = -1;
        for (size_t s = 0; s < nsupport_slots; s++) {
            size_t n = votes[c * nsupport_slots + s];
            if (n == 0) continue;
            distinct++;
            if (n > best_n) { best_n = n; best_s = (int32_t)s; }
        }
        m2s[c] = best_s;
        if (distinct > 1) ambiguous++;
    }
    free(votes);
    for (size_t v = 0; v < ps->nv; v++) {
        int32_t mc = set->vertex_mesh_component[v];
        vertex_support[v] = mc >= 0 && (size_t)mc < nmesh ? m2s[mc] : -1;
    }
    if (ambiguous > 0)
        fprintf(stderr,
            "[atlas_strip_scroll] collapse dump: %zu/%zu mesh components carry "
            "strokes from more than one support component; labelled by "
            "majority\n", ambiguous, nmesh);

    /* Rank support components by how many of their own wrap pairs collapsed. */
    size_t nsupport = set->stats.support_components + 1;
    size_t *coincident = (size_t *)calloc(nsupport, sizeof(size_t));
    size_t *pairs = (size_t *)calloc(nsupport, sizeof(size_t));
    if (coincident == NULL || pairs == NULL) {
        free(m2s); free(vertex_support); free(rgb);
        free(coincident); free(pairs);
        return -1;
    }
    for (size_t i = 0; i < order->npairs; i++) {
        const AtlasRadialOrderPair *p = &order->pairs[i];
        int32_t si = set->strokes[set->samples[p->inner].stroke].component;
        int32_t so = set->strokes[set->samples[p->outer].stroke].component;
        if (si != so || si < 0 || (size_t)si >= nsupport) continue;
        double t = fabs(p->spacing_target);
        if (!(t > 0.0)) continue;
        pairs[si]++;
        if (fabs(sample_u[p->outer] - sample_u[p->inner]) < 0.5 * t)
            coincident[si]++;
    }

    for (int rank = 0; rank < count; rank++) {
        int32_t best = -1;
        size_t best_n = 0;
        for (size_t s = 0; s < nsupport; s++)
            if (coincident[s] > best_n) { best_n = coincident[s]; best = (int32_t)s; }
        if (best < 0) break;
        char name[64];
        size_t nf_xyz = 0, nf_uv = 0, ncc = 0;
        /* One colour array shared by both dumps: the point is to follow a
         * physically separate piece from world space into the atlas. */
        collapse_component_colors(arena, ps, vertex_support, best, rgb, &ncc);
        snprintf(name, sizeof name, "collapse_support%02d_xyz.obj", best);
        if (write_component_obj(dir, name, ps, field_u, axial, vertex_support,
                                rgb, best, 0, &nf_xyz) != 0) break;
        snprintf(name, sizeof name, "collapse_support%02d_uv.obj", best);
        if (write_component_obj(dir, name, ps, field_u, axial, vertex_support,
                                rgb, best, 1, &nf_uv) != 0) break;
        snprintf(name, sizeof name, "collapse_support%02d_pairs_xyz.obj", best);
        write_collapse_pairs_obj(dir, name, ps, set, order, sample_u,
                                 vertex_support, best, 0);
        snprintf(name, sizeof name, "collapse_support%02d_pairs_uv.obj", best);
        write_collapse_pairs_obj(dir, name, ps, set, order, sample_u,
                                 vertex_support, best, 1);
        fprintf(stderr,
            "[atlas_strip_scroll] collapse dump: support %d -- %zu/%zu of its "
            "own wrap pairs coincident (%.1f%%), %zu faces in %zu XYZ "
            "connected components -> collapse_support%02d_{xyz,uv}.obj\n",
            best, coincident[best], pairs[best],
            pairs[best] ? 100.0 * (double)coincident[best] /
                          (double)pairs[best] : 0.0,
            nf_xyz, ncc, best);
        coincident[best] = 0;
    }

    free(m2s);
    free(vertex_support);
    free(rgb);
    free(coincident);
    free(pairs);
    return 0;
}

/* ==========================================================================
 * One support component at a time: find its sheets, then move each one out by
 * the whole number of turns it is short.
 *
 * This is the smallest honest test of the whole idea.  A support component is
 * where the collapse lives -- 80% of the coincident radial pairs are inside
 * one -- and inside one the strip solve pins u up to a single constant, so
 * nothing that moves whole components can help.  What CAN help is moving the
 * physically separate fragments the component is made of, and for that each
 * fragment needs an integer: which sheet of papyrus is it.
 *
 * The move is the exact spiral shear, evaluated at every vertex's own phi.
 * Collapsing it to a per-fragment constant is what turned earlier attempts
 * into teleportation (chart shifts of 24,513 u-units, weld stretch p95 1663).
 * v is not touched at all: the axial coordinate is the pinned calibration's
 * projection and was never a variable.
 * ========================================================================== */
static int write_sheet_csv(const char *dir, const AtlasSheetSplitResult *res)
{
    FILE *fp = as_open(dir, "sheets.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "weld_component,island,local_level,sheet,samples,pairs,"
                "pairs_internal,mean_radius,radius_min,radius_max,mean_phi,"
                "turn_span\n");
    for (size_t i = 0; i < res->ncomponents; i++) {
        const AtlasSheetComponent *c = &res->components[i];
        fprintf(fp, "%d,%d,%d,%d,%zu,%zu,%zu,%.6f,%.6f,%.6f,%.6f,%.4f\n",
                c->mesh_component, c->island, c->local_level, c->sheet,
                c->samples, c->pairs, c->pairs_internal, c->mean_radius,
                c->radius_min, c->radius_max, c->mean_phi, c->turn_span);
    }
    fclose(fp);
    fp = as_open(dir, "sheet_islands.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "island,components,samples,level_min,level_max,intercept,"
                "offset,residual_pitches,mean_radius\n");
    for (size_t i = 0; i < res->nislands; i++) {
        const AtlasSheetIsland *s = &res->islands[i];
        fprintf(fp, "%d,%zu,%zu,%d,%d,%.6f,%d,%.6f,%.6f\n",
                s->island, s->components, s->samples, s->level_min,
                s->level_max, s->intercept, s->offset, s->residual,
                s->mean_radius);
    }
    fclose(fp);
    return 0;
}

/*
 * Before and after in ONE file, stacked in v at the same scale.
 *
 * Rendered separately they auto-fit to their own bounding boxes, so a 549-vox
 * clump and a 30,477-vox spread come out the same size on screen and the
 * picture says nothing.  Here the two copies share a coordinate system and the
 * only difference visible is the one that matters.
 */
static int write_sheet_compare_obj(const char *dir, const char *name,
                                   const PieceSet *ps, const double *before,
                                   const double *after, const double *axial,
                                   const int32_t *vertex_in, const double *rgb)
{
    FILE *fp = as_open(dir, name, "wb");
    if (fp == NULL) return -1;
    double vlo = 0.0, vhi = 0.0;
    size_t seen = 0;
    for (size_t i = 0; i < ps->nv; i++) {
        if (!vertex_in[i]) continue;
        if (seen == 0) { vlo = vhi = axial[i]; }
        else {
            if (axial[i] < vlo) vlo = axial[i];
            if (axial[i] > vhi) vhi = axial[i];
        }
        seen++;
    }
    double lift = (vhi - vlo) * 1.6 + 1.0;
    int32_t *remap = (int32_t *)malloc(ps->nv * sizeof(int32_t));
    if (remap == NULL) { fclose(fp); return -1; }
    for (int copy = 0; copy < 2; copy++) {
        const double *u = copy == 0 ? before : after;
        double dv = copy == 0 ? lift : 0.0;
        for (size_t i = 0; i < ps->nv; i++) {
            if (!vertex_in[i]) continue;
            fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                    u[i], axial[i] + dv,
                    rgb[i * 3 + 0], rgb[i * 3 + 1], rgb[i * 3 + 2]);
        }
    }
    for (int copy = 0; copy < 2; copy++) {
        int32_t next = 0;
        for (size_t i = 0; i < ps->nv; i++)
            remap[i] = vertex_in[i] ? ++next : -1;
        int32_t base = copy * next;
        for (size_t f = 0; f < ps->nf; f++) {
            int32_t a = remap[ps->faces[f * 3 + 0]];
            int32_t b = remap[ps->faces[f * 3 + 1]];
            int32_t c = remap[ps->faces[f * 3 + 2]];
            if (a < 0 || b < 0 || c < 0) continue;
            fprintf(fp, "f %d %d %d\n", base + a, base + b, base + c);
        }
    }
    free(remap);
    fclose(fp);
    return 0;
}

typedef struct {
    size_t verts;
    size_t faces;
    double u_span_before;
    double u_span_after;
    double u_owed_arclength;   /* sum over pieces of |u(phi_max)-u(phi_min)| */
    double u_held_arclength;   /* sum over pieces of their actual u extent */
    double shift_min;
    double shift_max;
    size_t pairs_in_support;
    size_t coincident_before;
    size_t coincident_after;
    size_t satisfied_before;
    size_t satisfied_after;

    AtlasTurnAdvanceEdgeStats advance_stroke;
    AtlasTurnAdvanceEdgeStats advance_cross;
    AtlasTurnAdvanceEdgeStats advance_cont;
    AtlasTurnAdvanceEdgeStats advance_cross_xcube;
    AtlasTurnAdvanceEdgeStats advance_cont_xcube;
    size_t flat_pairs;
    size_t flat_torn;
    double flat_jump_median;
    double flat_jump_max;
    size_t flat_torn_after;
    double flat_jump_median_after;
    double flat_jump_max_after;
    size_t search_collisions_before;
    size_t search_collisions_after;
    size_t search_moved;
    int    search_converged;
    size_t search_pairs;
    size_t search_pairs_intra;
    size_t search_pairs_ok_before;
    size_t search_pairs_ok_after;
    size_t exact_before;
    size_t exact_after;
    size_t exact_touch_before;
    size_t exact_touch_after;
    double advance_fraction;
    double advance_step_median;
    double advance_circumference;
    size_t advance_steps;
    size_t advance_steps_backward;
} SheetStageReport;

static int write_sheet_summary(const char *dir,
                               const AtlasSheetSplitResult *res,
                               const SheetStageReport *rep,
                               int32_t target, int applied)
{
    FILE *fp = as_open(dir, "sheet_split_summary.json", "wb");
    if (fp == NULL) return -1;
    const AtlasSheetSplitStats *s = &res->stats;
    fprintf(fp,
        "{\n  \"support_component\": %d,\n"
        "  \"shift_applied\": %d,\n"
        "  \"weld_components\": %zu,\n"
        "  \"weld_components_with_radial_pairs\": %zu,\n"
        "  \"islands\": %zu,\n"
        "  \"singleton_islands\": %zu,\n"
        "  \"radial_pairs_used\": %zu,\n"
        "  \"radial_pairs_intra_component\": %zu,\n"
        "  \"component_edges\": %zu,\n"
        "  \"component_edges_dropped_weak\": %zu,\n"
        "  \"component_edges_in_forest\": %zu,\n"
        "  \"component_edges_consistent\": %zu,\n"
        "  \"component_edges_contradictory\": %zu,\n"
        "  \"sheets\": %zu,\n"
        "  \"sheet_min\": %d,\n"
        "  \"sheet_max\": %d,\n"
        "  \"island_align_residual_rms_pitches\": %.17g,\n"
        "  \"island_align_residual_max_pitches\": %.17g,\n"
        "  \"island_align_unreliable\": %zu,\n"
        "  \"phi_cross_check_agree\": %zu,\n"
        "  \"phi_cross_check_disagree\": %zu,\n"
        "  \"spiral_components\": %zu,\n"
        "  \"turn_span_max\": %.17g,\n"
        "  \"turn_span_total\": %.17g,\n"
        "  \"dumped_vertices\": %zu,\n"
        "  \"dumped_faces\": %zu,\n"
        "  \"u_span_before\": %.17g,\n"
        "  \"u_span_after\": %.17g,\n"
        "  \"u_owed_arclength\": %.17g,\n"
        "  \"u_held_arclength\": %.17g,\n"
        "  \"shift_min\": %.17g,\n"
        "  \"shift_max\": %.17g,\n"
        "  \"radial_pairs_in_support\": %zu,\n"
        "  \"coincident_before\": %zu,\n"
        "  \"coincident_after\": %zu,\n"
        "  \"satisfied_before\": %zu,\n"
        "  \"satisfied_after\": %zu,\n"
        "  \"advance_fraction_per_turn\": %.17g,\n"
        "  \"advance_step_median\": %.17g,\n"
        "  \"advance_circumference_median\": %.17g,\n"
        "  \"advance_steps\": %zu,\n"
        "  \"advance_steps_backward\": %zu,\n"
        "  \"advance_stroke\": {\"edges\": %zu, \"owed\": %.17g, "
        "\"delivered\": %.17g, \"ratio_median\": %.17g, \"collapsed\": %zu},\n"
        "  \"advance_cross_section\": {\"edges\": %zu, \"owed\": %.17g, "
        "\"delivered\": %.17g, \"ratio_median\": %.17g, \"collapsed\": %zu},\n"
        "  \"advance_continuation\": {\"edges\": %zu, \"owed\": %.17g, "
        "\"delivered\": %.17g, \"ratio_median\": %.17g, \"collapsed\": %zu},\n"
        "  \"advance_cross_section_cross_cube\": {\"edges\": %zu, "
        "\"owed\": %.17g, \"delivered\": %.17g, \"collapsed\": %zu},\n"
        "  \"advance_continuation_cross_cube\": {\"edges\": %zu, "
        "\"owed\": %.17g, \"delivered\": %.17g, \"collapsed\": %zu},\n"
        "  \"flat_test\": {\"adjacent_cube_pairs\": %zu, \"torn\": %zu, "
        "\"u_jump_median\": %.17g, \"u_jump_max\": %.17g, "
        "\"torn_after\": %zu, \"u_jump_median_after\": %.17g, "
        "\"u_jump_max_after\": %.17g},\n"
        "  \"place_search\": {\"collisions_before\": %zu, "
        "\"collisions_after\": %zu, \"pieces_moved\": %zu, "
        "\"converged\": %d, \"radial_pairs\": %zu, "
        "\"radial_pairs_unreachable\": %zu, "
        "\"radial_pairs_satisfied_before\": %zu, "
        "\"radial_pairs_satisfied_after\": %zu},\n"
        "  \"exact_overlap_audit\": {\"atlas_before\": %zu, "
        "\"atlas_after\": %zu, \"involving_support_before\": %zu, "
        "\"involving_support_after\": %zu}\n"
        "}\n",
        target, applied, s->components, s->components_with_pairs, s->islands,
        s->singleton_islands, s->pairs_used, s->pairs_intra_component,
        s->edges, s->edges_dropped_weak, s->edges_in_forest,
        s->edges_consistent, s->edges_contradictory, s->sheets, s->sheet_min,
        s->sheet_max, s->align_residual_rms, s->align_residual_max,
        s->align_unreliable, s->phi_agree, s->phi_disagree,
        s->spiral_components, s->turn_span_max, s->turn_span_total,
        rep->verts, rep->faces, rep->u_span_before, rep->u_span_after,
        rep->u_owed_arclength, rep->u_held_arclength,
        rep->shift_min, rep->shift_max, rep->pairs_in_support,
        rep->coincident_before, rep->coincident_after,
        rep->satisfied_before, rep->satisfied_after,
        rep->advance_fraction, rep->advance_step_median,
        rep->advance_circumference, rep->advance_steps,
        rep->advance_steps_backward,
        rep->advance_stroke.edges, rep->advance_stroke.owed_total,
        rep->advance_stroke.delivered_total, rep->advance_stroke.ratio_median,
        rep->advance_stroke.collapsed,
        rep->advance_cross.edges, rep->advance_cross.owed_total,
        rep->advance_cross.delivered_total, rep->advance_cross.ratio_median,
        rep->advance_cross.collapsed,
        rep->advance_cont.edges, rep->advance_cont.owed_total,
        rep->advance_cont.delivered_total, rep->advance_cont.ratio_median,
        rep->advance_cont.collapsed,
        rep->advance_cross_xcube.edges, rep->advance_cross_xcube.owed_total,
        rep->advance_cross_xcube.delivered_total,
        rep->advance_cross_xcube.collapsed,
        rep->advance_cont_xcube.edges, rep->advance_cont_xcube.owed_total,
        rep->advance_cont_xcube.delivered_total,
        rep->advance_cont_xcube.collapsed,
        rep->flat_pairs, rep->flat_torn, rep->flat_jump_median,
        rep->flat_jump_max, rep->flat_torn_after,
        rep->flat_jump_median_after, rep->flat_jump_max_after,
        rep->search_collisions_before, rep->search_collisions_after,
        rep->search_moved, rep->search_converged, rep->search_pairs,
        rep->search_pairs_intra, rep->search_pairs_ok_before,
        rep->search_pairs_ok_after, rep->exact_before, rep->exact_after,
        rep->exact_touch_before, rep->exact_touch_after);
    fclose(fp);
    return 0;
}

/*
 * Which fragments are physically the same sheet?  Ask the welder.
 *
 * The atlas fragments the scroll twice over: once because the surface really
 * is broken, and once because every cube is trimmed to its owned box so
 * adjacent cubes' charts never touch.  The second kind is an artifact, and the
 * XYZ/BPA seam weld is exactly the machine that repairs it -- it rolls a ball
 * across the ~1-vox trim gap at a radius capped so that 2*rho_max stays under
 * the 7-vox inter-wrap clearance.  That cap is the whole reason this is safe
 * to believe: the weld physically cannot join two different wraps, so every
 * bridge it builds is evidence of one continuous sheet.
 *
 * Using it here inverts the previous approach.  Before, fragments were the
 * atoms and their sheet indices had to be inferred, which left the turn graph
 * in six disconnected islands whose relative offsets could only be guessed
 * from radius -- and radius cannot carry that, because the axis wanders about
 * 85 vox (~7 pitches), so the guessed offsets were reading wander and inflated
 * seven sheets into eleven.  Welding first removes the guess: fragments the
 * weld joins are one sheet by construction, and only the survivors need a
 * turn count.
 *
 * The weld is built from XYZ alone.  No atlas quantity is admitted, so it
 * cannot inherit the error it is being used to diagnose.
 */
/* A refined vertex added by the seam refiner is a midpoint; walk back to a
 * source vertex so it can be looked up in the mesh-component labelling. */
static int32_t sheet_source_vertex(const AtlasXyzWeldTopology *t,
                                   size_t nv_src, int32_t v)
{
    int guard = 0;
    while (v >= 0 && (size_t)v >= nv_src && guard++ < 64) {
        size_t k = (size_t)v - nv_src;
        if (t->parent0 == NULL) return -1;
        v = t->parent0[k];
    }
    return v;
}

/*
 * The winding coordinate, branch-cut-free when compared LOCALLY.
 *
 *     w = sense * r / pitch - theta / 2pi
 *
 * is constant along one wrap and steps by ~1 across the inter-wrap gap, so for
 * two points a voxel or two apart |dw| ~ 1 means the edge between them is a
 * radial short-circuit rather than a piece of surface.  Comparing it only
 * across a short edge is what keeps theta's branch cut out of it -- the same
 * reason `wind_cut` works down at the r=12 core where no length threshold can.
 */
static double sheet_winding_at(const ScaffoldCalib *cal, const double e0[3],
                               const double e1[3], const double org[3],
                               const float *p)
{
    double d[3];
    for (int k = 0; k < 3; k++) d[k] = (double)p[k] - org[k];
    double a = d[0] * e0[0] + d[1] * e0[1] + d[2] * e0[2];
    double b = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2];
    double r = sqrt(a * a + b * b);
    double th = atan2(b, a);
    return (double)cal->sense * r / cal->pitch -
           th / (2.0 * 3.14159265358979323846);
}

static int sheet_weld_components(Arena_T arena, const PieceSet *ps,
                                 const ScaffoldCalib *cal,
                                 const AtlasCandidateSet *set,
                                 double wind_tol,
                                 int32_t *weld_of, size_t *out_nweld,
                                 size_t *out_bridges, size_t *out_relations,
                                 size_t *out_cut)
{
    size_t nmesh = set->mesh_components;
    for (size_t i = 0; i < nmesh; i++) weld_of[i] = (int32_t)i;
    *out_nweld = nmesh;
    *out_bridges = 0;
    *out_relations = 0;
    if (nmesh == 0 || ps->n_cubes < 2) return 0;

    Arena_T topo = Arena_new();
    BpaBridgeGate gate;
    memset(&gate, 0, sizeof gate);
    gate.umb_y = cal->axis_point[1];
    gate.umb_x = cal->axis_point[2];
    gate.pitch = cal->pitch;
    gate.tol = SEAM_WIND_TOL_DEFAULT_TURNS;
    gate.hard = SEAM_WIND_HARD_TOL_DEFAULT_TURNS;
    AtlasXyzWeldTopology topology;
    memset(&topology, 0, sizeof topology);
    int rc = -1;
    if (AtlasXyzWeldTopology_build(topo, ps, 128.0f, 1.5f, 0.0f, 6.0f,
                                   &gate, &topology) == 0 &&
        topology.nfaces > 0) {
        int32_t *face_chart = (int32_t *)ARENA_ALLOC(
            topo, (ps->nf ? ps->nf : 1) * sizeof(int32_t));
        for (size_t f = 0; f < ps->nf; f++)
            face_chart[f] = set->vertex_mesh_component[ps->faces[f * 3]];
        /*
         * Union the charts one BRIDGE EDGE at a time, not one aggregated
         * relation at a time, so the winding test can be applied where the
         * fusion is actually created.  The edge's two endpoints are a voxel or
         * two apart, which is exactly the regime in which the local dw test is
         * branch-cut free; an edge with |dw| past the tolerance is spanning
         * the inter-wrap gap rather than lying on the surface, and letting it
         * union its two charts is what builds a 25-turn body out of 328
         * fragments.
         */
        double e0[3], e1[3], org[3];
        {
            double ax[3] = {(double)cal->axis_dir[0], (double)cal->axis_dir[1],
                            (double)cal->axis_dir[2]};
            double n = sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
            if (!(n > 1e-15)) { Arena_dispose(&topo); return -1; }
            ax[0] /= n; ax[1] /= n; ax[2] /= n;
            double t[3] = {1.0, 0.0, 0.0};
            if (fabs(ax[0]) > 0.9) { t[0] = 0.0; t[1] = 1.0; }
            double d = t[0] * ax[0] + t[1] * ax[1] + t[2] * ax[2];
            for (int k = 0; k < 3; k++) e0[k] = t[k] - d * ax[k];
            double m = sqrt(e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2]);
            for (int k = 0; k < 3; k++) e0[k] /= m;
            e1[0] = ax[1] * e0[2] - ax[2] * e0[1];
            e1[1] = ax[2] * e0[0] - ax[0] * e0[2];
            e1[2] = ax[0] * e0[1] - ax[1] * e0[0];
            for (int k = 0; k < 3; k++) org[k] = (double)cal->axis_point[k];
        }
        UnionFind uf = UF_new(topo, (int32_t)nmesh);
        size_t joined = 0, cut = 0;
        for (size_t f = 0; f < topology.nfaces; f++) {
            for (int e = 0; e < 3; e++) {
                int32_t va = topology.faces[f * 3 + (size_t)e];
                int32_t vb = topology.faces[f * 3 + (size_t)((e + 1) % 3)];
                int32_t sa = sheet_source_vertex(&topology, ps->nv, va);
                int32_t sb = sheet_source_vertex(&topology, ps->nv, vb);
                if (sa < 0 || sb < 0) continue;
                if ((size_t)sa >= ps->nv || (size_t)sb >= ps->nv) continue;
                int32_t ca = set->vertex_mesh_component[sa];
                int32_t cb = set->vertex_mesh_component[sb];
                if (ca < 0 || cb < 0 || ca == cb) continue;
                if ((size_t)ca >= nmesh || (size_t)cb >= nmesh) continue;
                if (wind_tol > 0.0) {
                    double wa = sheet_winding_at(
                        cal, e0, e1, org, &topology.vertices[(size_t)va * 3]);
                    double wb = sheet_winding_at(
                        cal, e0, e1, org, &topology.vertices[(size_t)vb * 3]);
                    if (fabs(wa - wb) > wind_tol) { cut++; continue; }
                }
                uf_union(&uf, ca, cb);
                joined++;
            }
        }
        {
            int32_t *label = (int32_t *)ARENA_ALLOC(
                topo, nmesh * sizeof(int32_t));
            for (size_t i = 0; i < nmesh; i++) label[i] = -1;
            size_t next = 0;
            for (size_t i = 0; i < nmesh; i++) {
                int32_t r = uf_find(&uf, (int32_t)i);
                if (label[r] < 0) label[r] = (int32_t)next++;
                weld_of[i] = label[r];
            }
            *out_nweld = next;
            *out_bridges = topology.nfaces;
            *out_relations = joined;
            if (out_cut != NULL) *out_cut = cut;
            rc = 0;
        }
    }
    Arena_dispose(&topo);
    (void)arena;
    return rc;
}

/*
 * The flat test: cubes that touch in XYZ should touch in u.
 *
 * Inside one welded piece the surface is physically continuous, so two cubes
 * sharing a face hold papyrus a couple of voxels apart and their u medians
 * should differ by roughly the width of a cube, not by a circumference.  This
 * needs no phi and no spiral -- just cube origins from the ids and the field
 * itself -- so unlike the advance audit it cannot be blamed on registration.
 * A jump near a full circumference is a tear, and a tear inside a physically
 * continuous piece is a bug wherever it came from.
 *
 * Run it on the field going in and the field coming out: that A/B is what says
 * whether a move introduced tearing or merely failed to remove it.
 */
typedef struct {
    size_t pairs;
    size_t torn;
    double jump_median;
    double jump_p95;
    double jump_max;
} FlatTestStats;

static int sheet_flat_test(Arena_T arena, const PieceSet *ps,
                           const AtlasCandidateSet *set,
                           const int32_t *sample_support,
                           const int32_t *sample_unit, size_t ns,
                           size_t nweld, int32_t target, const double *u,
                           double tear_threshold,
                           size_t *per_piece_pairs, size_t *per_piece_torn,
                           double *per_piece_max, FlatTestStats *out)
{
    memset(out, 0, sizeof *out);
    size_t ncube = ps->n_cubes;
    if (ncube == 0 || nweld == 0) return 0;
    Arena_Mark mark = Arena_save(arena);
    double *cu = (double *)ARENA_ALLOC(arena, ncube * nweld * sizeof(double));
    size_t *cn = (size_t *)ARENA_CALLOC(arena, ncube * nweld, sizeof(size_t));
    for (size_t i = 0; i < ns; i++) {
        if (!sheet_in_target(target, sample_support[i])) continue;
        int32_t w = sample_unit[i];
        int32_t c = set->sample_ref[i].cube;
        if (w < 0 || (size_t)w >= nweld) continue;
        if (c < 0 || (size_t)c >= ncube) continue;
        size_t slot = (size_t)w * ncube + (size_t)c;
        if (cn[slot]++ == 0) cu[slot] = u[i];
        else cu[slot] += (u[i] - cu[slot]) / (double)cn[slot];
    }
    double *jump = (double *)ARENA_ALLOC(
        arena, (ncube * ncube + 1) * sizeof(double));
    size_t njump = 0;
    for (size_t w = 0; w < nweld; w++) {
        if (per_piece_pairs != NULL) per_piece_pairs[w] = 0;
        if (per_piece_torn != NULL) per_piece_torn[w] = 0;
        if (per_piece_max != NULL) per_piece_max[w] = 0.0;
        for (size_t a = 0; a < ncube; a++) {
            if (cn[w * ncube + a] == 0) continue;
            for (size_t b = a + 1; b < ncube; b++) {
                if (cn[w * ncube + b] == 0) continue;
                long d0 = labs(ps->cube_org[a][0] - ps->cube_org[b][0]);
                long d1 = labs(ps->cube_org[a][1] - ps->cube_org[b][1]);
                long d2 = labs(ps->cube_org[a][2] - ps->cube_org[b][2]);
                if (!((d0 == 128 && d1 == 0 && d2 == 0) ||
                      (d1 == 128 && d0 == 0 && d2 == 0) ||
                      (d2 == 128 && d0 == 0 && d1 == 0)))
                    continue;
                double dj = fabs(cu[w * ncube + a] - cu[w * ncube + b]);
                out->pairs++;
                jump[njump++] = dj;
                if (dj > out->jump_max) out->jump_max = dj;
                if (per_piece_pairs != NULL) per_piece_pairs[w]++;
                if (per_piece_max != NULL && dj > per_piece_max[w])
                    per_piece_max[w] = dj;
                if (dj > tear_threshold) {
                    out->torn++;
                    if (per_piece_torn != NULL) per_piece_torn[w]++;
                }
            }
        }
    }
    qsort(jump, njump, sizeof(double), compare_double_value);
    out->jump_median = njump ? jump[njump / 2] : 0.0;
    out->jump_p95 = njump ? jump[(size_t)((double)njump * 0.95)] : 0.0;
    Arena_restore(arena, mark);
    return 0;
}

/* A negative target means every support component at once: the welded pieces
 * are a global partition anyway, and the placement search has always scored
 * against the whole atlas, so running them all together is the case it was
 * built for rather than a generalisation of it. */
/*
 * Dump one welded piece for inspection, coloured by what it was BEFORE the
 * weld.
 *
 * The interesting question about a piece that overlaps itself is which join
 * made it one piece, and the pre-weld mesh component is exactly that record:
 * every colour boundary inside the dump is a place the BPA seam bridge decided
 * two fragments were the same surface.  Alongside it, `_tears.obj` draws a
 * segment across every pair of cubes that touch in XYZ inside this piece yet
 * sit more than half a circumference apart in u.  Those segments are the
 * suspect welds -- a physically continuous surface cannot jump a wrap between
 * adjacent cubes, so either the weld joined two wraps or the atlas is torn
 * there, and both are visible at the segment's ends.
 */
static int write_piece_dump(Arena_T arena, const char *dir, int32_t slot,
                            const PieceSet *ps, const AtlasCandidateSet *set,
                            const int32_t *vertex_slot, const double *u,
                            const double *axial, double tear,
                            size_t *out_verts, size_t *out_faces,
                            size_t *out_fragments, size_t *out_tears)
{
    Arena_Mark mark = Arena_save(arena);
    size_t nmesh = set->mesh_components;

    /* Rank the pre-weld fragments by size so colours are stable run to run. */
    size_t *fcount = (size_t *)ARENA_CALLOC(arena, nmesh ? nmesh : 1,
                                            sizeof(size_t));
    for (size_t v = 0; v < ps->nv; v++) {
        if (vertex_slot[v] != slot) continue;
        int32_t mc = set->vertex_mesh_component[v];
        if (mc >= 0 && (size_t)mc < nmesh) fcount[mc]++;
    }
    int32_t *frank = (int32_t *)ARENA_ALLOC(arena, (nmesh ? nmesh : 1) *
                                                   sizeof(int32_t));
    size_t nfrag = 0;
    for (size_t i = 0; i < nmesh; i++) frank[i] = fcount[i] > 0 ? 0 : -1;
    for (;;) {
        size_t best = nmesh, bn = 0;
        for (size_t i = 0; i < nmesh; i++)
            if (frank[i] == 0 && fcount[i] > bn) { bn = fcount[i]; best = i; }
        if (best == nmesh) break;
        frank[best] = (int32_t)++nfrag;
    }

    char name[64];
    size_t nv_out = 0, nf_out = 0;
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena, ps->nv * sizeof(int32_t));
    for (int flat = 0; flat < 2; flat++) {
        snprintf(name, sizeof name, "piece%02d_%s.obj", slot,
                 flat ? "uv" : "xyz");
        FILE *fp = as_open(dir, name, "wb");
        if (fp == NULL) { Arena_restore(arena, mark); return -1; }
        int32_t next = 0;
        for (size_t v = 0; v < ps->nv; v++) {
            remap[v] = -1;
            if (vertex_slot[v] != slot) continue;
            remap[v] = ++next;
            double rgb[3] = {0.5, 0.5, 0.5};
            int32_t mc = set->vertex_mesh_component[v];
            if (mc >= 0 && (size_t)mc < nmesh && frank[mc] > 0)
                collapse_hue_rgb((size_t)(frank[mc] - 1), rgb);
            if (flat)
                fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                        u[v], axial[v], rgb[0], rgb[1], rgb[2]);
            else
                fprintf(fp, "v %.9g %.9g %.9g %.6f %.6f %.6f\n",
                        (double)ps->verts[v * 3 + 0],
                        (double)ps->verts[v * 3 + 1],
                        (double)ps->verts[v * 3 + 2],
                        rgb[0], rgb[1], rgb[2]);
        }
        nv_out = (size_t)next;
        nf_out = 0;
        for (size_t f = 0; f < ps->nf; f++) {
            int32_t a = remap[ps->faces[f * 3 + 0]];
            int32_t b = remap[ps->faces[f * 3 + 1]];
            int32_t c = remap[ps->faces[f * 3 + 2]];
            if (a < 0 || b < 0 || c < 0) continue;
            fprintf(fp, "f %d %d %d\n", a, b, c);
            nf_out++;
        }
        fclose(fp);
    }

    /* --- the suspect welds ---------------------------------------------- */
    size_t ncube = ps->n_cubes;
    double *cu = (double *)ARENA_CALLOC(arena, ncube ? ncube : 1,
                                        sizeof(double));
    double *cv = (double *)ARENA_CALLOC(arena, ncube ? ncube : 1,
                                        sizeof(double));
    double *cx = (double *)ARENA_CALLOC(arena, ncube * 3 + 1, sizeof(double));
    size_t *cn = (size_t *)ARENA_CALLOC(arena, ncube ? ncube : 1,
                                        sizeof(size_t));
    for (size_t v = 0; v < ps->nv; v++) {
        if (vertex_slot[v] != slot) continue;
        int32_t c = -1;
        for (size_t k = 0; k < ncube; k++)
            if (v >= ps->cube_voff[k] && v < ps->cube_voff[k + 1]) {
                c = (int32_t)k;
                break;
            }
        if (c < 0) continue;
        size_t n = ++cn[c];
        cu[c] += (u[v] - cu[c]) / (double)n;
        cv[c] += (axial[v] - cv[c]) / (double)n;
        for (int k = 0; k < 3; k++)
            cx[(size_t)c * 3 + (size_t)k] +=
                ((double)ps->verts[v * 3 + (size_t)k] -
                 cx[(size_t)c * 3 + (size_t)k]) / (double)n;
    }
    size_t ntear = 0;
    for (int flat = 0; flat < 2; flat++) {
        snprintf(name, sizeof name, "piece%02d_tears_%s.obj", slot,
                 flat ? "uv" : "xyz");
        FILE *fp = as_open(dir, name, "wb");
        if (fp == NULL) { Arena_restore(arena, mark); return -1; }
        size_t n = 0;
        for (size_t a = 0; a < ncube; a++) {
            if (cn[a] == 0) continue;
            for (size_t b = a + 1; b < ncube; b++) {
                if (cn[b] == 0) continue;
                long d0 = labs(ps->cube_org[a][0] - ps->cube_org[b][0]);
                long d1 = labs(ps->cube_org[a][1] - ps->cube_org[b][1]);
                long d2 = labs(ps->cube_org[a][2] - ps->cube_org[b][2]);
                if (!((d0 == 128 && d1 == 0 && d2 == 0) ||
                      (d1 == 128 && d0 == 0 && d2 == 0) ||
                      (d2 == 128 && d0 == 0 && d1 == 0)))
                    continue;
                if (fabs(cu[a] - cu[b]) <= tear) continue;
                for (int e = 0; e < 2; e++) {
                    size_t c = e == 0 ? a : b;
                    if (flat)
                        fprintf(fp, "v %.9g %.9g 0\n", cu[c], cv[c]);
                    else
                        fprintf(fp, "v %.9g %.9g %.9g\n",
                                cx[c * 3 + 0], cx[c * 3 + 1], cx[c * 3 + 2]);
                }
                n++;
            }
        }
        for (size_t i = 0; i < n; i++)
            fprintf(fp, "l %zu %zu\n", 2 * i + 1, 2 * i + 2);
        fclose(fp);
        ntear = n;
    }

    if (out_verts != NULL) *out_verts = nv_out;
    if (out_faces != NULL) *out_faces = nf_out;
    if (out_fragments != NULL) *out_fragments = nfrag;
    if (out_tears != NULL) *out_tears = ntear;
    Arena_restore(arena, mark);
    return 0;
}

/* Red (t=0) through yellow/green/cyan to blue (t=1). */
static void sheet_ramp_rgb(double t, double rgb[3])
{
    if (!(t > 0.0)) t = 0.0;
    if (t > 1.0) t = 1.0;
    double h = t * 240.0;   /* hue 0 = red at t=0, hue 240 = blue at t=1 */
    double s = 1.0, v = 1.0;
    double c = v * s;
    double hp = h / 60.0;
    double x = c * (1.0 - fabs(fmod(hp, 2.0) - 1.0));
    double r = 0.0, g = 0.0, b = 0.0;
    if (hp < 1.0)      { r = c; g = x; }
    else if (hp < 2.0) { r = x; g = c; }
    else if (hp < 3.0) { g = c; b = x; }
    else if (hp < 4.0) { g = x; b = c; }
    else if (hp < 5.0) { r = x; b = c; }
    else               { r = c; b = x; }
    rgb[0] = r; rgb[1] = g; rgb[2] = b;
}

/*
 * The arclength picture: what u IS against what u OWES.
 *
 * Under the pinned map du/dphi = r(phi), so u is arclength and
 * spiral_u(phi_v) is exactly where vertex v belongs once the piece is panned
 * out flat.  Colouring the piece by the ACTUAL field on the scale of the OWED
 * one therefore shows the compression directly: if the solve advances u by 3%
 * of the arclength it owes, the whole piece comes out red and the blue end of
 * the ramp is simply never reached.  The companion `_expected` dump paints the
 * same geometry with spiral_u itself, which is what a correct solve would look
 * like, and `_actual` re-normalises to the field's own range so whatever
 * structure u does have is still visible instead of collapsing to one colour.
 *
 * phi is per-cube registered, which matters here only as a caveat on the
 * absolute span: across twenty-odd turns the real winding dominates the
 * registration steps, so the ramp is honest about the shape even where it is
 * approximate about the ends.
 */
static int write_piece_arclength_dump(Arena_T arena, const char *dir,
                                      int32_t slot, const PieceSet *ps,
                                      const ScaffoldCalib *cal,
                                      const int32_t *vertex_slot,
                                      const double *u, const double *axial,
                                      double *out_owed, double *out_held,
                                      double *out_fraction)
{
    Arena_Mark mark = Arena_save(arena);
    double axis[3] = {(double)cal->axis_dir[0], (double)cal->axis_dir[1],
                      (double)cal->axis_dir[2]};
    double an = sqrt(axis[0] * axis[0] + axis[1] * axis[1] +
                     axis[2] * axis[2]);
    if (!(an > 1e-15)) { Arena_restore(arena, mark); return -1; }
    axis[0] /= an; axis[1] /= an; axis[2] /= an;

    /*
     * Both maps share ONE origin and ONE scale, so the two pictures can be
     * laid side by side and read against each other.  The origin is the
     * core-most vertex and the far end is the rim-most, which orients t from 0
     * to 1 along the spiral no matter which way the pinned sense runs u.
     */
    /* In-plane frame, so theta is available alongside radius. */
    double f0[3], f1[3];
    {
        double t[3] = {1.0, 0.0, 0.0};
        if (fabs(axis[0]) > 0.9) { t[0] = 0.0; t[1] = 1.0; }
        double d = t[0] * axis[0] + t[1] * axis[1] + t[2] * axis[2];
        for (int k = 0; k < 3; k++) f0[k] = t[k] - d * axis[k];
        double m = sqrt(f0[0] * f0[0] + f0[1] * f0[1] + f0[2] * f0[2]);
        for (int k = 0; k < 3; k++) f0[k] /= m;
        f1[0] = axis[1] * f0[2] - axis[2] * f0[1];
        f1[1] = axis[2] * f0[0] - axis[0] * f0[2];
        f1[2] = axis[0] * f0[1] - axis[1] * f0[0];
    }

    double *want = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    /*
     * A SECOND reference that never touches phi.
     *
     * phi's integer part is a per-cube rounding the cross-cube registration
     * never applied, and it disagrees with radius on 43% of radial pairs -- so
     * an "expected" built from phi inherits exactly the defect the picture is
     * supposed to expose.  This one re-derives the winding from geometry with
     * ONE globally uniform rule, `ribbon.c`'s own anchor
     *     phi_geo = theta + 2pi * round(sense*r/pitch - theta/2pi)
     * applied per vertex instead of per chain group per cube.  It is still a
     * rounding, so axis wander can still move it, but it cannot be
     * INCONSISTENT between neighbouring cubes, which is the failure mode at
     * issue.  Disagreement between the two is reported in turns.
     */
    double *want_geo = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    double *radius = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    double *dturn = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    size_t vref = 0, vend = 0;
    double rmin = 0.0, rmax = 0.0, alo = 0.0, ahi = 0.0;
    size_t n = 0;
    const double two_pi = 2.0 * 3.14159265358979323846;
    for (size_t v = 0; v < ps->nv; v++) {
        if (vertex_slot[v] != slot) continue;
        want[v] = spiral_u(cal, (double)ps->phi[v]);
        double d[3], along = 0.0;
        for (int k = 0; k < 3; k++) {
            d[k] = (double)ps->verts[v * 3 + (size_t)k] -
                   (double)cal->axis_point[k];
            along += d[k] * axis[k];
        }
        double pa = 0.0, pb = 0.0;
        for (int k = 0; k < 3; k++) {
            double q = d[k] - along * axis[k];
            pa += q * f0[k];
            pb += q * f1[k];
        }
        double r = sqrt(pa * pa + pb * pb);
        radius[v] = r;
        double th = atan2(pb, pa);
        double wgeo = (double)cal->sense * r / cal->pitch - th / two_pi;
        double phi_geo = th + two_pi * floor(wgeo + 0.5);
        want_geo[v] = spiral_u(cal, phi_geo);
        dturn[n] = ((double)ps->phi[v] - phi_geo) / two_pi;
        if (n == 0) {
            rmin = rmax = r; vref = vend = v; alo = ahi = u[v];
        } else {
            if (r < rmin) { rmin = r; vref = v; }
            if (r > rmax) { rmax = r; vend = v; }
            if (u[v] < alo) alo = u[v];
            if (u[v] > ahi) ahi = u[v];
        }
        n++;
    }
    if (n == 0) { Arena_restore(arena, mark); return 0; }
    double denom = want[vend] - want[vref];   /* signed: core -> rim */
    if (fabs(denom) < 1e-9) denom = 1.0;
    double wlo = want[vref], uref = u[vref];
    double owed = fabs(denom), held = ahi - alo;
    double dgeo = want_geo[vend] - want_geo[vref];
    if (fabs(dgeo) < 1e-9) dgeo = 1.0;

    /* How far apart are the two references, in whole turns? */
    {
        qsort(dturn, n, sizeof(double), compare_double_value);
        double med = dturn[n / 2];
        size_t within_half = 0;
        double lo95 = dturn[(size_t)((double)n * 0.025)];
        double hi95 = dturn[(size_t)((double)n * 0.975)];
        for (size_t i = 0; i < n; i++)
            if (fabs(dturn[i] - med) <= 0.5) within_half++;
        fprintf(stderr,
            "[atlas_strip_scroll] REFERENCE CHECK slot %d: phi vs the "
            "radius-anchored winding differ by median %+.2f turns, 95%% band "
            "[%+.2f %+.2f] = %.1f turns wide; only %zu/%zu vertices (%.1f%%) "
            "sit within half a turn of the median. %s\n",
            slot, med, lo95, hi95, hi95 - lo95, within_half, n,
            100.0 * (double)within_half / (double)n,
            hi95 - lo95 > 1.0
                ? "phi is NOT a usable reference here -- the expected map "
                  "built from it is unreliable"
                : "the two references agree; phi is usable here");
        fprintf(stderr,
            "[atlas_strip_scroll]   owed span: from phi %.0f vox, from the "
            "radius anchor %.0f vox\n", fabs(denom), fabs(dgeo));
    }

    static const char *mode_name[4] = {"owed", "expected", "actual",
                                       "expected_geo"};
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena, ps->nv * sizeof(int32_t));
    for (int mode = 0; mode < 4; mode++) {
        for (int flat = 0; flat < 2; flat++) {
            char name[80];
            snprintf(name, sizeof name, "piece%02d_arc_%s_%s.obj", slot,
                     mode_name[mode], flat ? "uv" : "xyz");
            FILE *fp = as_open(dir, name, "wb");
            if (fp == NULL) { Arena_restore(arena, mark); return -1; }
            int32_t next = 0;
            for (size_t v = 0; v < ps->nv; v++) {
                remap[v] = -1;
                if (vertex_slot[v] != slot) continue;
                remap[v] = ++next;
                double t;
                if (mode == 0)      t = (u[v] - uref) / denom;    /* actual, on the owed scale */
                else if (mode == 1) t = (want[v] - wlo) / denom;  /* target, from phi */
                else if (mode == 3) t = (want_geo[v] - want_geo[vref]) / dgeo;
                else t = held > 1e-9 ? (u[v] - alo) / held : 0.0; /* actual, own range */
                double rgb[3];
                sheet_ramp_rgb(t, rgb);
                if (flat)
                    fprintf(fp, "v %.9g %.9g 0 %.6f %.6f %.6f\n",
                            u[v], axial[v], rgb[0], rgb[1], rgb[2]);
                else
                    fprintf(fp, "v %.9g %.9g %.9g %.6f %.6f %.6f\n",
                            (double)ps->verts[v * 3 + 0],
                            (double)ps->verts[v * 3 + 1],
                            (double)ps->verts[v * 3 + 2],
                            rgb[0], rgb[1], rgb[2]);
            }
            for (size_t f = 0; f < ps->nf; f++) {
                int32_t a = remap[ps->faces[f * 3 + 0]];
                int32_t b = remap[ps->faces[f * 3 + 1]];
                int32_t c = remap[ps->faces[f * 3 + 2]];
                if (a < 0 || b < 0 || c < 0) continue;
                fprintf(fp, "f %d %d %d\n", a, b, c);
            }
            fclose(fp);
        }
    }
    /*
     * Is it a SHAPE error or a SCALE error?
     *
     * Fit actual = alpha * expected + beta over the piece's vertices.  If the
     * field is the right function of the winding and merely under-driven, the
     * fit is near-perfect and alpha is the gain that is missing -- a very
     * different bug from a field that has the wrong shape, and a far easier
     * one.  The residual rms after the fit is what separates the two, so both
     * are reported rather than inferred from the picture.
     */
    /*
     * Which way does u run?
     *
     * A global sign on u is only a gauge, but the PINNED map fixes a
     * convention and `CubeReg_deltaU` follows it, so a field running the other
     * way means every winding shift computed from the pinned map points
     * backwards.  Regress both the field and the reference against radius; if
     * the two slopes disagree in sign, that is the whole explanation for a
     * negative fit slope, and it is a fact about orientation rather than about
     * shape.
     */
    for (int ref = 0; ref < 3; ref++) {
        const double *x = ref == 0 ? want : ref == 1 ? want_geo : radius;
        double sx = 0.0, sy = 0.0;
        for (size_t v = 0; v < ps->nv; v++) {
            if (vertex_slot[v] != slot) continue;
            sx += x[v];
            sy += u[v];
        }
        double mx = sx / (double)n, my = sy / (double)n;
        double sxx = 0.0, sxy = 0.0, syy = 0.0;
        for (size_t v = 0; v < ps->nv; v++) {
            if (vertex_slot[v] != slot) continue;
            double dx = x[v] - mx, dy = u[v] - my;
            sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
        }
        double alpha = sxx > 1e-12 ? sxy / sxx : 0.0;
        double beta = my - alpha * mx;
        double ss = 0.0, worst = 0.0;
        for (size_t v = 0; v < ps->nv; v++) {
            if (vertex_slot[v] != slot) continue;
            double e = u[v] - (alpha * x[v] + beta);
            ss += e * e;
            if (fabs(e) > worst) worst = fabs(e);
        }
        double rms = sqrt(ss / (double)n);
        double r2 = (sxx > 1e-12 && syy > 1e-12)
                  ? (sxy * sxy) / (sxx * syy) : 0.0;
        fprintf(stderr,
            "[atlas_strip_scroll] ARCLENGTH FIT slot %d vs %-13s: actual = "
            "%+.4f * expected %+.1f, r2 = %.4f, residual rms %.1f vox "
            "(worst %.0f). %s\n",
            slot, ref == 0 ? "phi" : ref == 1 ? "radius anchor" : "RADIUS",
            alpha, beta, r2, rms,
            worst,
            r2 > 0.95
                ? "SHAPE IS RIGHT -- a GAIN error, not a topology one"
                : "shape does not fit a pure gain");
    }

    if (out_owed != NULL) *out_owed = owed;
    if (out_held != NULL) *out_held = held;
    if (out_fraction != NULL) *out_fraction = held / owed;
    Arena_restore(arena, mark);
    return 0;
}

static int sheet_in_target(int32_t target, int32_t support)
{
    return target < 0 || support == target;
}

static int run_sheet_split_stage(Arena_T arena, const ScrollConfig *cfg,
                                 const PieceSet *ps, const ScaffoldCalib *cal,
                                 const AtlasCandidateSet *set,
                                 const AtlasRadialOrderSet *order,
                                 const int32_t *sample_support,
                                 const double *sample_phi,
                                 const double *base_sample,
                                 const double *base_field,
                                 const double *axial,
                                 double **out_field)
{
    char dir[AS_PATH_CAP], probe[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "sheet_split") != 0 ||
        as_path(probe, dir, ".probe") != 0 ||
        ves_ensure_parent_dir(probe) != 0)
        return -1;

    size_t ns = set->problem.nsamples;
    int32_t target = cfg->sheet_all ? -1 : (int32_t)cfg->isolate_support;

    /* Sample radius: pure geometry off the pinned axis, the same quantity the
     * ray builder ranks by. */
    double axis[3] = {(double)cal->axis_dir[0], (double)cal->axis_dir[1],
                      (double)cal->axis_dir[2]};
    double anorm = sqrt(axis[0] * axis[0] + axis[1] * axis[1] +
                        axis[2] * axis[2]);
    if (!(anorm > 1e-15)) return -1;
    axis[0] /= anorm; axis[1] /= anorm; axis[2] /= anorm;

    int32_t *sample_component = (int32_t *)ARENA_ALLOC(
        arena, ns * sizeof(int32_t));
    double *sample_radius = (double *)ARENA_ALLOC(arena, ns * sizeof(double));
    for (size_t i = 0; i < ns; i++) {
        const AtlasCandidateSampleRef *ref = &set->sample_ref[i];
        sample_component[i] = ref->mesh_component;
        const float *v0 = &ps->verts[ref->mesh_vertex[0] * 3];
        const float *v1 = &ps->verts[ref->mesh_vertex[1] * 3];
        double d[3], along = 0.0;
        for (int k = 0; k < 3; k++) {
            double p = (double)v0[k] + ((double)v1[k] - (double)v0[k]) *
                                       ref->mesh_t;
            d[k] = p - (double)cal->axis_point[k];
            along += d[k] * axis[k];
        }
        double r2 = 0.0;
        for (int k = 0; k < 3; k++) {
            double t = d[k] - along * axis[k];
            r2 += t * t;
        }
        sample_radius[i] = sqrt(r2);
    }

    /* --- physical connectivity first, sheet indices second --------------- */
    size_t nmesh = set->mesh_components;
    int32_t *weld_of = (int32_t *)ARENA_ALLOC(
        arena, (nmesh ? nmesh : 1) * sizeof(int32_t));
    size_t nweld = nmesh, bridge_faces = 0, weld_relations = 0, weld_cut = 0;
    double t0 = ves_clock_sec();
    if (cfg->sheet_no_weld) {
        for (size_t i = 0; i < nmesh; i++) weld_of[i] = (int32_t)i;
    } else if (sheet_weld_components(arena, ps, cal, set, cfg->sheet_wind_tol,
                                     weld_of, &nweld, &bridge_faces,
                                     &weld_relations, &weld_cut) != 0) {
        fprintf(stderr, "atlas_strip_scroll: XYZ weld topology failed\n");
        return -1;
    }
    /* How much did the weld actually merge inside the target support? */
    size_t mesh_in_target = 0, weld_in_target = 0;
    {
        uint8_t *mseen = (uint8_t *)ARENA_CALLOC(arena, nmesh ? nmesh : 1,
                                                 sizeof(uint8_t));
        uint8_t *wseen = (uint8_t *)ARENA_CALLOC(arena, nweld ? nweld : 1,
                                                 sizeof(uint8_t));
        for (size_t i = 0; i < ns; i++) {
            if (!sheet_in_target(target, sample_support[i])) continue;
            int32_t mc = sample_component[i];
            if (mc < 0 || (size_t)mc >= nmesh) continue;
            if (!mseen[mc]) { mseen[mc] = 1; mesh_in_target++; }
            int32_t wc = weld_of[mc];
            if (wc >= 0 && (size_t)wc < nweld && !wseen[wc]) {
                wseen[wc] = 1;
                weld_in_target++;
            }
        }
    }
    fprintf(stderr,
        "[atlas_strip_scroll] XYZ weld %.2fs: %zu bridge faces, %zu chart "
        "joins accepted, %zu CUT as cross-wrap (|dw| > %.2f); mesh components "
        "%zu -> %zu welded overall, %zu -> %zu inside support %d\n",
        ves_clock_sec() - t0, bridge_faces, weld_relations, weld_cut,
        cfg->sheet_wind_tol, nmesh, nweld,
        mesh_in_target, weld_in_target, target);

    /* The sheet solve now runs over WELDED pieces, not raw fragments. */
    int32_t *sample_unit = (int32_t *)ARENA_ALLOC(arena, ns * sizeof(int32_t));
    for (size_t i = 0; i < ns; i++) {
        int32_t mc = sample_component[i];
        sample_unit[i] = mc >= 0 && (size_t)mc < nmesh ? weld_of[mc] : -1;
    }

    AtlasSheetSplitResult res;
    t0 = ves_clock_sec();
    if (AtlasSheetSplit_solve(arena, order->pairs, order->npairs,
                              sample_unit, sample_support, sample_radius,
                              sample_phi, ns, nweld, target,
                              cal->pitch, (double)cal->sense, &cfg->sheet,
                              &res) != 0) {
        fprintf(stderr, "atlas_strip_scroll: sheet split failed\n");
        return -1;
    }
    const AtlasSheetSplitStats *st = &res.stats;
    fprintf(stderr,
        "[atlas_strip_scroll] sheet split %.2fs on support %d: %zu welded "
        "pieces (%zu carry radial evidence), %zu islands (%zu singleton), "
        "%zu radial pairs -> %zu piece edges (%zu dropped as too weak, "
        "%zu in the forest, %zu consistent, %zu CONTRADICTORY)\n",
        ves_clock_sec() - t0, target, st->components,
        st->components_with_pairs, st->islands, st->singleton_islands,
        st->pairs_used, st->edges, st->edges_dropped_weak,
        st->edges_in_forest, st->edges_consistent, st->edges_contradictory);
    fprintf(stderr,
        "[atlas_strip_scroll] SHEETS DETECTED: %zu distinct, indices %d..%d; "
        "island alignment residual rms %.3f max %.3f pitches (%zu islands "
        "unreliable); phi cross-check agrees on %zu/%zu components\n",
        st->sheets, st->sheet_min, st->sheet_max, st->align_residual_rms,
        st->align_residual_max, st->align_unreliable, st->phi_agree,
        st->phi_agree + st->phi_disagree);
    if (st->pairs_intra_component > 0 || st->spiral_components > 0)
        fprintf(stderr,
            "[atlas_strip_scroll] SPIRAL CHECK: %zu/%zu welded pieces span a "
            "whole turn or more (worst %.2f turns, %.2f turns of papyrus in "
            "total); %zu radial pairs have BOTH ends inside one piece. A "
            "piece that spirals past itself has no single sheet index -- its "
            "wrap number changes along its own length\n",
            st->spiral_components, st->components, st->turn_span_max,
            st->turn_span_total, st->pairs_intra_component);
    /*
     * Is the atlas even the right SIZE?
     *
     * The pinned map u(phi) = a*phi + b*phi^2/(4pi) has du/dphi = r(phi), so u
     * is arclength along the spiral and the u a piece is owed is just
     * u(phi_max) - u(phi_min).  Comparing that against the u it actually
     * occupies separates two very different failures: pieces sitting on top of
     * each other (a placement error, which whole-turn shifts fix) and a piece
     * compressed along its own length (a scale error, which they cannot).
     */
    double *piece_u_lo = (double *)ARENA_ALLOC(
        arena, (res.ncomponents ? res.ncomponents : 1) * sizeof(double));
    double *piece_u_hi = (double *)ARENA_ALLOC(
        arena, (res.ncomponents ? res.ncomponents : 1) * sizeof(double));
    size_t *piece_seen = (size_t *)ARENA_CALLOC(
        arena, res.ncomponents ? res.ncomponents : 1, sizeof(size_t));
    int32_t *unit_slot = (int32_t *)ARENA_ALLOC(
        arena, (nweld ? nweld : 1) * sizeof(int32_t));
    for (size_t i = 0; i < nweld; i++) unit_slot[i] = -1;
    for (size_t i = 0; i < res.ncomponents; i++)
        unit_slot[res.components[i].mesh_component] = (int32_t)i;
    for (size_t i = 0; i < ns; i++) {
        if (!sheet_in_target(target, sample_support[i])) continue;
        int32_t u0 = sample_unit[i];
        if (u0 < 0 || (size_t)u0 >= nweld) continue;
        int32_t s = unit_slot[u0];
        if (s < 0) continue;
        if (piece_seen[s]++ == 0) {
            piece_u_lo[s] = piece_u_hi[s] = base_sample[i];
        } else {
            if (base_sample[i] < piece_u_lo[s]) piece_u_lo[s] = base_sample[i];
            if (base_sample[i] > piece_u_hi[s]) piece_u_hi[s] = base_sample[i];
        }
    }
    double owed_total = 0.0, held_total = 0.0;
    for (size_t i = 0; i < res.ncomponents; i++) {
        const AtlasSheetComponent *c = &res.components[i];
        double owed = fabs(spiral_u(cal, c->phi_max) -
                           spiral_u(cal, c->phi_min));
        double held = piece_seen[i] ? piece_u_hi[i] - piece_u_lo[i] : 0.0;
        owed_total += owed;
        held_total += held;
        fprintf(stderr,
            "[atlas_strip_scroll]   piece %d: sheet %d, %zu samples, r "
            "%.1f..%.1f, %.2f turns of phi, %zu radial pairs (%zu internal); "
            "u owed %.0f, holds %.0f -- compressed %.0fx\n",
            c->mesh_component, c->sheet, c->samples, c->radius_min,
            c->radius_max, c->turn_span, c->pairs, c->pairs_internal,
            owed, held, held > 1e-9 ? owed / held : 0.0);
    }
    for (size_t i = 0; i < res.nislands; i++) {
        const AtlasSheetIsland *is = &res.islands[i];
        fprintf(stderr,
            "[atlas_strip_scroll]   island %d: %zu comps, %zu samples, levels "
            "%+d..%+d, mean r %.1f, intercept %.1f, offset %+d (residual "
            "%+.3f pitches)\n",
            is->island, is->components, is->samples, is->level_min,
            is->level_max, is->mean_radius, is->intercept, is->offset,
            is->residual);
    }

    /* --- select the component's geometry --------------------------------- */
    int32_t *vertex_in = (int32_t *)ARENA_ALLOC(arena,
                                                ps->nv * sizeof(int32_t));
    int32_t *vertex_sheet = (int32_t *)ARENA_ALLOC(arena,
                                                   ps->nv * sizeof(int32_t));
    int32_t *vertex_slot = (int32_t *)ARENA_ALLOC(arena,
                                                  ps->nv * sizeof(int32_t));
    for (size_t v = 0; v < ps->nv; v++) {
        int32_t mc = set->vertex_mesh_component[v];
        int32_t wc = mc >= 0 && (size_t)mc < nmesh ? weld_of[mc] : -1;
        int32_t sh = wc >= 0 && (size_t)wc < res.nmesh ? res.sheet_of[wc] : -1;
        vertex_sheet[v] = sh;
        vertex_slot[v] = wc >= 0 && (size_t)wc < nweld ? unit_slot[wc] : -1;
        vertex_in[v] = sh >= 0 ? 1 : 0;
    }

    double *rgb = (double *)ARENA_ALLOC(arena, ps->nv * 3 * sizeof(double));
    for (size_t v = 0; v < ps->nv; v++) {
        rgb[v * 3 + 0] = rgb[v * 3 + 1] = rgb[v * 3 + 2] = 0.5;
        if (vertex_sheet[v] < 0) continue;
        collapse_hue_rgb((size_t)vertex_sheet[v], &rgb[v * 3]);
    }

    /* --- the move -------------------------------------------------------- */
    double *shifted = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    SheetStageReport rep;
    memset(&rep, 0, sizeof rep);
    rep.shift_min = 0.0;
    rep.shift_max = 0.0;
    int applied = cfg->sheet_no_shift ? 0 : 1;
    /*
     * A one-winding move is NOT a translation.  Under the pinned spiral,
     * u(phi + 2pi k) - u(phi) = 2pi a k + b k phi + pi b k^2, and the middle
     * term varies with phi -- so the same "one turn out" is a SHEAR of
     * b*k*dphi across a block, ~479 vox over the six turns one of these
     * pieces spans.  `--sheet-rigid` evaluates it once at the block's median
     * phi instead, which is the honest way to ask what a rigid move does.
     */
    size_t nshift = 0;
    /*
     * Per-piece reference phi.  The move is RIGID, so each piece needs exactly
     * one phi at which to evaluate the winding step -- its own median, not the
     * per-vertex value.  Using the per-vertex value is what injected phi's
     * unrepaired per-cube rounding into the shift and degraded cube adjacency
     * 5x; evaluating once per piece removes phi from the interior of the move
     * entirely, which is why the rigid variant preserves adjacency exactly.
     */
    double *piece_phi = (double *)ARENA_ALLOC(
        arena, (res.ncomponents ? res.ncomponents : 1) * sizeof(double));
    {
        double *pv = (double *)ARENA_ALLOC(arena, (ps->nv ? ps->nv : 1) *
                                                  sizeof(double));
        for (size_t i = 0; i < res.ncomponents; i++) {
            size_t np = 0;
            int32_t want = res.components[i].mesh_component;
            for (size_t vtx = 0; vtx < ps->nv; vtx++) {
                int32_t mc = set->vertex_mesh_component[vtx];
                if (mc < 0 || (size_t)mc >= nmesh) continue;
                if (weld_of[mc] != want) continue;
                pv[np++] = (double)ps->phi[vtx];
            }
            if (np == 0) { piece_phi[i] = res.components[i].mean_phi; continue; }
            qsort(pv, np, sizeof(double), compare_double_value);
            piece_phi[i] = np % 2 ? pv[np / 2]
                                  : 0.5 * (pv[np / 2 - 1] + pv[np / 2]);
        }
    }

    /*
     * Search the wind rather than inherit it.
     *
     * The turn graph's integers come from phi, and phi's integer part is a
     * per-cube rounding that the cross-cube registration never repaired --
     * 622 group nodes, 258 edges, 373 components, 235 frustrated, zero
     * corrections applied.  So the wind is chosen here by asking the only
     * question that matters: at which whole wind does this piece land on
     * nothing?  Scored against the WHOLE atlas, because a piece cleared of its
     * own neighbours can still land on another support component.
     */
    int32_t *searched = NULL;
    if (cfg->sheet_search && res.ncomponents > 0) {
        size_t np = res.ncomponents;
        int span = cfg->sheet_search_span > 0 ? cfg->sheet_search_span : 4;
        size_t ncand = (size_t)(2 * span + 1);
        double *cand_shift = (double *)ARENA_ALLOC(
            arena, np * ncand * sizeof(double));
        double *cand_cost = (double *)ARENA_ALLOC(
            arena, np * ncand * sizeof(double));
        for (size_t p = 0; p < np; p++)
            for (size_t j = 0; j < ncand; j++) {
                int32_t k = (int32_t)j - span;
                cand_shift[p * ncand + j] = k == 0 ? 0.0
                    : CubeReg_deltaU(cal->spiral_a, cal->spiral_b,
                                     piece_phi[p],
                                     (int32_t)(k * cal->sense));
                cand_cost[p * ncand + j] = fabs((double)k);
            }
        int32_t *vertex_piece = (int32_t *)ARENA_ALLOC(
            arena, ps->nv * sizeof(int32_t));
        for (size_t vtx = 0; vtx < ps->nv; vtx++) {
            int32_t mc = set->vertex_mesh_component[vtx];
            int32_t wc = mc >= 0 && (size_t)mc < nmesh ? weld_of[mc] : -1;
            int32_t slot = wc >= 0 && (size_t)wc < nweld ? unit_slot[wc] : -1;
            vertex_piece[vtx] = slot;   /* -1 = held fixed as background */
        }
        /*
         * The radial evidence rides along as the second objective term.
         * Collision-freedom alone is degenerate -- any permutation that tiles
         * scores perfectly -- and the pairs are the thing that says which
         * permutation is physically right.  86% of them are intra-cube, so
         * they do not inherit the unrepaired cross-cube winding registration
         * that makes phi's own integer untrustworthy.
         */
        AtlasPlacePair *ppair = (AtlasPlacePair *)ARENA_ALLOC(
            arena, (order->npairs ? order->npairs : 1) *
                   sizeof(AtlasPlacePair));
        size_t nppair = 0;
        for (size_t i = 0; i < order->npairs; i++) {
            const AtlasRadialOrderPair *rp = &order->pairs[i];
            if (!(fabs(rp->spacing_target) > 0.0)) continue;
            int32_t wi = sample_unit[rp->inner];
            int32_t wo = sample_unit[rp->outer];
            int32_t si2 = wi >= 0 && (size_t)wi < nweld ? unit_slot[wi] : -1;
            int32_t so2 = wo >= 0 && (size_t)wo < nweld ? unit_slot[wo] : -1;
            if (si2 < 0 && so2 < 0) continue;   /* neither end can move */
            ppair[nppair].piece_a = si2;
            ppair[nppair].piece_b = so2;
            ppair[nppair].du_base =
                base_sample[rp->outer] - base_sample[rp->inner];
            ppair[nppair].target = rp->spacing_target;
            nppair++;
        }

        AtlasPlaceSearchOptions so;
        AtlasPlaceSearchOptions_default(&so);
        if (cfg->sheet_search_cell > 0.0) so.cell = cfg->sheet_search_cell;
        if (cfg->sheet_search_lambda > 0.0)
            so.lambda_move = cfg->sheet_search_lambda;
        if (cfg->sheet_search_pair >= 0.0)
            so.lambda_pair = cfg->sheet_search_pair;
        int32_t *choice = (int32_t *)ARENA_ALLOC(arena, np * sizeof(int32_t));
        AtlasPlaceSearchStats ss;
        double tsearch = ves_clock_sec();
        if (AtlasPlaceSearch_solve(arena, base_field, axial, ps->nv,
                                   vertex_piece, np, ncand, cand_shift,
                                   cand_cost, ppair, nppair, &so, choice,
                                   &ss) != 0) {
            fprintf(stderr, "atlas_strip_scroll: placement search failed\n");
            return -1;
        }
        searched = (int32_t *)ARENA_ALLOC(arena, np * sizeof(int32_t));
        for (size_t p = 0; p < np; p++)
            searched[p] = (int32_t)choice[p] - span;
        fprintf(stderr,
            "[atlas_strip_scroll] PLACE SEARCH %.2fs: %zu pieces x %zu winds, "
            "%zu candidates scored; raster cell %.1f, %zu background cells, "
            "%zu moving cells; collisions %zu -> %zu (%zu against the rest of "
            "the atlas); %zu pieces moved; %d rounds%s\n",
            ves_clock_sec() - tsearch, ss.pieces, ncand,
            ss.candidates_evaluated, so.cell, ss.static_cells,
            ss.moving_cells, ss.collisions_before, ss.collisions_after,
            ss.collisions_static_after, ss.moved_pieces, ss.rounds_run,
            ss.converged ? ", converged" : ", HIT THE ROUND CAP");
        fprintf(stderr,
            "[atlas_strip_scroll] PLACE SEARCH radial term (lambda=%.1f): "
            "%zu pairs, %zu unreachable (both ends in one piece); satisfied "
            "%zu -> %zu\n",
            so.lambda_pair, ss.pairs, ss.pairs_intra,
            ss.pairs_satisfied_before, ss.pairs_satisfied_after);
        for (size_t p = 0; p < np; p++)
            fprintf(stderr,
                "[atlas_strip_scroll]   piece %d: phi turn graph said sheet "
                "%+d, search chose wind %+d (shift %.0f)\n",
                res.components[p].mesh_component, res.components[p].sheet,
                searched[p], cand_shift[p * ncand + (size_t)choice[p]]);
        rep.search_collisions_before = ss.collisions_before;
        rep.search_collisions_after = ss.collisions_after;
        rep.search_moved = ss.moved_pieces;
        rep.search_converged = ss.converged;
        rep.search_pairs = ss.pairs;
        rep.search_pairs_intra = ss.pairs_intra;
        rep.search_pairs_ok_before = ss.pairs_satisfied_before;
        rep.search_pairs_ok_after = ss.pairs_satisfied_after;
    }

    double phi_ref = 0.0;
    if (cfg->sheet_rigid) {
        double *pv = (double *)ARENA_ALLOC(arena, (ps->nv ? ps->nv : 1) *
                                                  sizeof(double));
        size_t np = 0;
        for (size_t v = 0; v < ps->nv; v++)
            if (vertex_sheet[v] >= 0) pv[np++] = (double)ps->phi[v];
        if (np > 0) {
            qsort(pv, np, sizeof(double), compare_double_value);
            phi_ref = np % 2 ? pv[np / 2]
                             : 0.5 * (pv[np / 2 - 1] + pv[np / 2]);
        }
    }
    for (size_t v = 0; v < ps->nv; v++) {
        shifted[v] = base_field[v];
        if (!applied || vertex_sheet[v] < 0) continue;
        int32_t slot = vertex_slot[v];
        int32_t turns = cfg->sheet_uniform_turns != 0
                      ? cfg->sheet_uniform_turns
                      : (searched != NULL && slot >= 0 ? searched[slot]
                                                       : vertex_sheet[v]);
        if (turns == 0) continue;
        int32_t k = (int32_t)(turns * cal->sense);
        /* Rigid: one phi per piece, so the move carries no phi structure into
         * the piece's interior and cube adjacency survives untouched. */
        double at = cfg->sheet_rigid ? phi_ref
                  : (searched != NULL && slot >= 0 ? piece_phi[slot]
                                                   : (double)ps->phi[v]);
        double dv = CubeReg_deltaU(cal->spiral_a, cal->spiral_b, at, k);
        shifted[v] += dv;
        /* Seed from the first moved vertex: the shifts are all one sign, so
         * initialising the extremes to zero silently reports the wrong
         * spread. */
        if (nshift++ == 0) { rep.shift_min = rep.shift_max = dv; }
        else {
            if (dv < rep.shift_min) rep.shift_min = dv;
            if (dv > rep.shift_max) rep.shift_max = dv;
        }
    }
    if (cfg->sheet_uniform_turns != 0)
        fprintf(stderr,
            "[atlas_strip_scroll] UNIFORM MOVE: the whole block by %d "
            "winding(s), %s; shift range %.1f .. %.1f (spread %.1f vox -- "
            "that spread IS the shear, b*k*dphi)\n",
            cfg->sheet_uniform_turns,
            cfg->sheet_rigid ? "rigid, at the block's median phi"
                             : "the exact phi-dependent shear",
            rep.shift_min, rep.shift_max, rep.shift_max - rep.shift_min);

    /* --- did it help?  measured on this component's own radial pairs ----- */
    for (size_t i = 0; i < order->npairs; i++) {
        const AtlasRadialOrderPair *p = &order->pairs[i];
        if (!sheet_in_target(target, sample_support[p->inner]) ||
            !sheet_in_target(target, sample_support[p->outer]))
            continue;
        double want = fabs(p->spacing_target);
        if (!(want > 0.0)) continue;
        rep.pairs_in_support++;
        double before = base_sample[p->outer] - base_sample[p->inner];
        int32_t si = sample_unit[p->inner] >= 0
                   ? res.sheet_of[sample_unit[p->inner]] : -1;
        int32_t so = sample_unit[p->outer] >= 0
                   ? res.sheet_of[sample_unit[p->outer]] : -1;
        double after = before;
        if (applied) {
            /* Score whatever move was actually applied, not the default one:
             * a uniform move gives every sample the same turn count and so
             * changes no relative position at all. */
            int32_t sli = sample_unit[p->inner] >= 0
                        ? unit_slot[sample_unit[p->inner]] : -1;
            int32_t slo = sample_unit[p->outer] >= 0
                        ? unit_slot[sample_unit[p->outer]] : -1;
            int32_t ti = cfg->sheet_uniform_turns != 0
                       ? cfg->sheet_uniform_turns
                       : (searched != NULL && sli >= 0 ? searched[sli] : si);
            int32_t to = cfg->sheet_uniform_turns != 0
                       ? cfg->sheet_uniform_turns
                       : (searched != NULL && slo >= 0 ? searched[slo] : so);
            double ai = cfg->sheet_rigid ? phi_ref
                      : (searched != NULL && sli >= 0 ? piece_phi[sli]
                                                      : sample_phi[p->inner]);
            double ao = cfg->sheet_rigid ? phi_ref
                      : (searched != NULL && slo >= 0 ? piece_phi[slo]
                                                      : sample_phi[p->outer]);
            double di = si >= 0 && ti != 0
                ? CubeReg_deltaU(cal->spiral_a, cal->spiral_b, ai,
                                 (int32_t)(ti * cal->sense)) : 0.0;
            double dobj = so >= 0 && to != 0
                ? CubeReg_deltaU(cal->spiral_a, cal->spiral_b, ao,
                                 (int32_t)(to * cal->sense)) : 0.0;
            after = before + dobj - di;
        }
        if (fabs(before) < 0.5 * want) rep.coincident_before++;
        else if (fabs(before - p->spacing_target) < 0.5 * want)
            rep.satisfied_before++;
        if (fabs(after) < 0.5 * want) rep.coincident_after++;
        else if (fabs(after - p->spacing_target) < 0.5 * want)
            rep.satisfied_after++;
    }

    double lo_b = 0.0, hi_b = 0.0, lo_a = 0.0, hi_a = 0.0;
    for (size_t v = 0; v < ps->nv; v++) {
        if (!vertex_in[v]) continue;
        if (rep.verts == 0) {
            lo_b = hi_b = base_field[v];
            lo_a = hi_a = shifted[v];
        } else {
            if (base_field[v] < lo_b) lo_b = base_field[v];
            if (base_field[v] > hi_b) hi_b = base_field[v];
            if (shifted[v] < lo_a) lo_a = shifted[v];
            if (shifted[v] > hi_a) hi_a = shifted[v];
        }
        rep.verts++;
    }
    rep.u_span_before = hi_b - lo_b;
    rep.u_span_after = hi_a - lo_a;
    rep.u_owed_arclength = owed_total;
    rep.u_held_arclength = held_total;

    /*
     * Where does the missing advance go?
     *
     * Three kinds of relation carry u between samples and they promise
     * different things, so each is scored against its own target: strokes owe
     * the spiral arc they span, continuations owe their own measured physical
     * gap, and cross-sections owe NOTHING by design -- they are alignment
     * rows.  That last one is why they are scored inverted: a cross-section
     * whose members span real angle is not failing to advance u, it is
     * actively instructing the solve to destroy that advance, and the amount
     * it destroys is the arc its members span.
     */
    {
        int32_t *ea = (int32_t *)ARENA_ALLOC(
            arena, (set->problem.nmembers + set->problem.nsamples +
                    set->problem.ncontinuations + 1) * sizeof(int32_t));
        int32_t *eb = (int32_t *)ARENA_ALLOC(
            arena, (set->problem.nmembers + set->problem.nsamples +
                    set->problem.ncontinuations + 1) * sizeof(int32_t));
        double *etarget = (double *)ARENA_ALLOC(
            arena, (set->problem.ncontinuations + 1) * sizeof(double));
        AtlasTurnAdvanceEdgeStats st_stroke, st_cross, st_cont;
        memset(&st_stroke, 0, sizeof st_stroke);
        memset(&st_cross, 0, sizeof st_cross);
        memset(&st_cont, 0, sizeof st_cont);

        /* strokes: consecutive samples, the unit-speed rows */
        size_t ne = 0;
        for (size_t s = 0; s < set->problem.nstrokes; s++) {
            const AtlasStripStroke *sk = &set->problem.strokes[s];
            for (int32_t j = 0; j + 1 < sk->count; j++) {
                int32_t a = sk->first + j;
                if (!sheet_in_target(target, sample_support[a])) continue;
                ea[ne] = a;
                eb[ne] = a + 1;
                ne++;
            }
        }
        AtlasTurnAdvance_edges(arena, base_sample, sample_phi, ea, eb, NULL,
                               ne, ns, cal->spiral_a, cal->spiral_b, 0,
                               &st_stroke);

        /*
         * Cross-sections and continuations are scored TWICE, split by whether
         * the relation stays inside one cube.  phi is registered per cube
         * upstream, so a relation crossing a seam can look a whole turn wide
         * purely from registration -- that false positive is what overturned
         * the previous diagnosis, and cutting on it tore the physical weld
         * (oracle p95 3.0 -> 2089).  Only the same-cube column is evidence.
         */
        int32_t *xa = (int32_t *)ARENA_ALLOC(
            arena, (set->problem.nmembers + 1) * sizeof(int32_t));
        int32_t *xb = (int32_t *)ARENA_ALLOC(
            arena, (set->problem.nmembers + 1) * sizeof(int32_t));
        AtlasTurnAdvanceEdgeStats st_cross_x, st_cont_x;
        memset(&st_cross_x, 0, sizeof st_cross_x);
        memset(&st_cont_x, 0, sizeof st_cont_x);

        ne = 0;
        size_t nx = 0;
        for (size_t c = 0; c < set->problem.ncross_sections; c++) {
            const AtlasStripCrossSection *cs = &set->problem.cross_sections[c];
            int32_t prev = -1;
            for (int32_t j = 0; j < cs->count; j++) {
                size_t mi = cs->first + (size_t)j;
                if (set->member_state[mi] == ATLAS_CANDIDATE_INACTIVE) continue;
                const AtlasStripMember *m = &set->problem.members[mi];
                int32_t v = m->value1 >= 0 && m->value_t > 0.5
                          ? m->value1 : m->value0;
                if (v < 0 || (size_t)v >= ns) continue;
                if (!sheet_in_target(target, sample_support[v])) { prev = -1; continue; }
                if (prev >= 0) {
                    if (set->sample_ref[prev].cube == set->sample_ref[v].cube) {
                        ea[ne] = prev; eb[ne] = v; ne++;
                    } else {
                        xa[nx] = prev; xb[nx] = v; nx++;
                    }
                }
                prev = v;
            }
        }
        AtlasTurnAdvance_edges(arena, base_sample, sample_phi, ea, eb, NULL,
                               ne, ns, cal->spiral_a, cal->spiral_b, 1,
                               &st_cross);
        AtlasTurnAdvance_edges(arena, base_sample, sample_phi, xa, xb, NULL,
                               nx, ns, cal->spiral_a, cal->spiral_b, 1,
                               &st_cross_x);

        /* continuations: each carries its own measured gap */
        ne = 0;
        nx = 0;
        for (size_t i = 0; i < set->problem.ncontinuations; i++) {
            if (set->continuation_state[i] == ATLAS_CANDIDATE_INACTIVE)
                continue;
            const AtlasStripContinuation *cn = &set->problem.continuations[i];
            if (cn->a < 0 || (size_t)cn->a >= ns) continue;
            if (cn->b < 0 || (size_t)cn->b >= ns) continue;
            if (!sheet_in_target(target, sample_support[cn->a])) continue;
            if (set->sample_ref[cn->a].cube == set->sample_ref[cn->b].cube) {
                ea[ne] = cn->a;
                eb[ne] = cn->b;
                etarget[ne] = -cn->target;  /* the row is u[a]-u[b] = target */
                ne++;
            } else {
                xa[nx] = cn->a; xb[nx] = cn->b; nx++;
            }
        }
        AtlasTurnAdvance_edges(arena, base_sample, sample_phi, ea, eb,
                               etarget, ne, ns, cal->spiral_a, cal->spiral_b,
                               0, &st_cont);
        AtlasTurnAdvance_edges(arena, base_sample, sample_phi, xa, xb, NULL,
                               nx, ns, cal->spiral_a, cal->spiral_b, 0,
                               &st_cont_x);

        const char *label[5] = {"stroke            ",
                                "cross-sec same-cube",
                                "cross-sec CROSS-cube",
                                "continue  same-cube",
                                "continue  CROSS-cube"};
        const AtlasTurnAdvanceEdgeStats *all[5] = {
            &st_stroke, &st_cross, &st_cross_x, &st_cont, &st_cont_x};
        for (int k = 0; k < 5; k++) {
            const AtlasTurnAdvanceEdgeStats *e = all[k];
            fprintf(stderr,
                "[atlas_strip_scroll] ADVANCE %s: %zu edges, median span "
                "%.3f turns, owes %.0f delivers %.0f (%.1f%%), median ratio "
                "%.3f, %zu collapsed, %zu reversed\n",
                label[k], e->edges, e->dphi_turns_median, e->owed_total,
                e->delivered_total,
                e->owed_total > 1e-9
                    ? 100.0 * e->delivered_total / e->owed_total : 0.0,
                e->ratio_median, e->collapsed, e->reversed);
        }
        rep.advance_stroke = st_stroke;
        rep.advance_cross = st_cross;
        rep.advance_cont = st_cont;
        rep.advance_cross_xcube = st_cross_x;
        rep.advance_cont_xcube = st_cont_x;
    }

    /*
     * The flat test: cubes that touch in XYZ should touch in u.
     *
     * Inside one welded piece the surface is physically continuous, so two
     * cubes sharing a face hold papyrus a couple of voxels apart and their u
     * medians should differ by roughly the width of a cube, not by a
     * circumference.  This needs no phi and no spiral -- just cube origins
     * from the ids and the field itself -- so unlike everything above it
     * cannot be blamed on registration.  A jump near a full circumference is
     * a tear, and a tear inside a physically continuous piece is a bug.
     */
    {
        /* Same-piece samples carried through the same move as the vertices. */
        double *shifted_sample = (double *)ARENA_ALLOC(
            arena, ns * sizeof(double));
        for (size_t i = 0; i < ns; i++) {
            shifted_sample[i] = base_sample[i];
            if (!applied || !sheet_in_target(target, sample_support[i])) continue;
            int32_t w = sample_unit[i];
            int32_t slot = w >= 0 && (size_t)w < nweld ? unit_slot[w] : -1;
            if (slot < 0) continue;
            int32_t turns = cfg->sheet_uniform_turns != 0
                          ? cfg->sheet_uniform_turns
                          : (searched != NULL ? searched[slot]
                                              : res.components[slot].sheet);
            if (turns == 0) continue;
            double at = cfg->sheet_rigid ? phi_ref
                      : (searched != NULL ? piece_phi[slot] : sample_phi[i]);
            shifted_sample[i] += CubeReg_deltaU(
                cal->spiral_a, cal->spiral_b, at,
                (int32_t)(turns * cal->sense));
        }
        /* Half a local circumference: unambiguously a wrap, not a cube's
         * worth of surface. */
        double tear = 0.5 * 2.0 * 3.14159265358979323846 * 350.0;
        size_t *pp = (size_t *)ARENA_ALLOC(arena, (nweld ? nweld : 1) *
                                                  sizeof(size_t));
        size_t *pt = (size_t *)ARENA_ALLOC(arena, (nweld ? nweld : 1) *
                                                  sizeof(size_t));
        double *pm = (double *)ARENA_ALLOC(arena, (nweld ? nweld : 1) *
                                                  sizeof(double));
        size_t *pt_a = (size_t *)ARENA_ALLOC(arena, (nweld ? nweld : 1) *
                                                    sizeof(size_t));
        double *pm_a = (double *)ARENA_ALLOC(arena, (nweld ? nweld : 1) *
                                                    sizeof(double));
        FlatTestStats before, after;
        sheet_flat_test(arena, ps, set, sample_support, sample_unit, ns,
                        nweld, target, base_sample, tear, pp, pt, pm,
                        &before);
        sheet_flat_test(arena, ps, set, sample_support, sample_unit, ns,
                        nweld, target, shifted_sample, tear, NULL, pt_a,
                        pm_a, &after);
        for (size_t w = 0; w < nweld; w++) {
            if (pp[w] == 0) continue;
            int32_t slot = unit_slot[w];
            fprintf(stderr,
                "[atlas_strip_scroll]   FLAT piece %zu (sheet %d): %zu "
                "face-adjacent cube pairs; torn %zu -> %zu, worst jump "
                "%.0f -> %.0f vox\n",
                w, slot >= 0 ? res.components[slot].sheet : -1, pp[w],
                pt[w], pt_a[w], pm[w], pm_a[w]);
        }
        fprintf(stderr,
            "[atlas_strip_scroll] FLAT TEST on support %d: %zu face-adjacent "
            "cube pairs inside a single welded piece. BEFORE median %.1f p95 "
            "%.1f worst %.0f, %zu torn. AFTER median %.1f p95 %.1f worst "
            "%.0f, %zu torn\n",
            target, before.pairs, before.jump_median, before.jump_p95,
            before.jump_max, before.torn, after.jump_median, after.jump_p95,
            after.jump_max, after.torn);
        rep.flat_pairs = before.pairs;
        rep.flat_torn = before.torn;
        rep.flat_jump_median = before.jump_median;
        rep.flat_jump_max = before.jump_max;
        rep.flat_torn_after = after.torn;
        rep.flat_jump_median_after = after.jump_median;
        rep.flat_jump_max_after = after.jump_max;
    }

    /* Per-turn advance: how far does u move between adjacent turns? */
    {
        int32_t *sample_piece = (int32_t *)ARENA_ALLOC(
            arena, ns * sizeof(int32_t));
        for (size_t i = 0; i < ns; i++)
            sample_piece[i] = sheet_in_target(target, sample_support[i])
                            ? sample_unit[i] : -1;
        AtlasTurnAdvanceResult adv;
        if (AtlasTurnAdvance_bins(arena, base_sample, sample_phi,
                                  sample_radius, sample_piece, ns, nweld,
                                  (double)cal->sense, &adv) == 0) {
            fprintf(stderr,
                "[atlas_strip_scroll] ADVANCE per turn: %zu turn bins, %zu "
                "adjacent-turn steps (%zu backward); median step %.1f vox "
                "against a %.1f vox circumference = %.1f%% of one turn\n",
                adv.nbins, adv.steps, adv.steps_backward, adv.step_median,
                adv.circumference_median, 100.0 * adv.advance_fraction);
            rep.advance_fraction = adv.advance_fraction;
            rep.advance_steps = adv.steps;
            rep.advance_steps_backward = adv.steps_backward;
            rep.advance_step_median = adv.step_median;
            rep.advance_circumference = adv.circumference_median;
        }
    }

    fprintf(stderr,
        "[atlas_strip_scroll] SCALE on support %d: its pieces are owed %.0f "
        "vox of u along their own arclength and hold %.0f -- compressed "
        "%.0fx. The whole-turn shifts add %.0f BETWEEN pieces, which is the "
        "correct size; it is the blocks that are ~%.0fx too narrow, so the "
        "gaps only LOOK explosive\n",
        target, owed_total, held_total,
        held_total > 1e-9 ? owed_total / held_total : 0.0,
        rep.u_span_after - rep.u_span_before,
        held_total > 1e-9 ? owed_total / held_total : 0.0);

    /*
     * Verify with the exact tri-tri audit.
     *
     * The search scores a 4-vox occupancy raster, which is a proxy: it counts
     * cells two pieces share, not triangles that actually intersect.  The
     * audit is the project's ground truth for overlap, so the placement the
     * proxy chose is checked against it before any of these numbers are
     * quoted.  It runs on the WHOLE atlas, so it also catches a support-3
     * piece that cleared its own neighbours and landed on someone else's.
     */
    {
        float *registered = (float *)ARENA_ALLOC(arena,
                                                 ps->nv * sizeof(float));
        for (size_t v = 0; v < ps->nv; v++) registered[v] = ps->uv[v * 2];
        AtlasOverlapAudit oa_before, oa_after;
        double taudit = ves_clock_sec();
        if (AtlasOverlapAudit_build(arena, ps->faces, ps->nf, ps->nv,
                                    base_field, axial, registered, ps->phi,
                                    set->vertex_mesh_component,
                                    set->mesh_components, &oa_before) == 0 &&
            AtlasOverlapAudit_build(arena, ps->faces, ps->nf, ps->nv,
                                    shifted, axial, registered, ps->phi,
                                    set->vertex_mesh_component,
                                    set->mesh_components, &oa_after) == 0) {
            /* Split by whether the pair involves the pieces we moved: a rigid
             * move cannot change overlap among faces that all stayed put. */
            size_t touch_before = 0, touch_after = 0;
            for (size_t i = 0; i < oa_before.npairs; i++) {
                int32_t f0 = oa_before.pairs[i].face0;
                int32_t f1 = oa_before.pairs[i].face1;
                if (vertex_in[ps->faces[(size_t)f0 * 3]] ||
                    vertex_in[ps->faces[(size_t)f1 * 3]]) touch_before++;
            }
            for (size_t i = 0; i < oa_after.npairs; i++) {
                int32_t f0 = oa_after.pairs[i].face0;
                int32_t f1 = oa_after.pairs[i].face1;
                if (vertex_in[ps->faces[(size_t)f0 * 3]] ||
                    vertex_in[ps->faces[(size_t)f1 * 3]]) touch_after++;
            }
            fprintf(stderr,
                "[atlas_strip_scroll] EXACT OVERLAP AUDIT %.2fs: whole atlas "
                "%zu -> %zu face pairs; pairs involving support %d: %zu -> "
                "%zu\n",
                ves_clock_sec() - taudit, oa_before.npairs, oa_after.npairs,
                target, touch_before, touch_after);
            rep.exact_before = oa_before.npairs;
            rep.exact_after = oa_after.npairs;
            rep.exact_touch_before = touch_before;
            rep.exact_touch_after = touch_after;

            /*
             * Whose overlap is it?
             *
             * A total is not actionable.  Attribute every surviving pair to
             * the welded piece(s) its two faces belong to and split SELF from
             * CROSS: self-overlap is a piece folded onto itself, which no
             * rigid move of that piece can ever fix and which points at the
             * weld having joined two wraps; cross-overlap is a placement
             * failure between two pieces, which the search can in principle
             * still resolve.  The two want completely different repairs, so
             * the split is the first thing to look at.
             */
            size_t np2 = res.ncomponents ? res.ncomponents : 1;
            size_t *ov_self = (size_t *)ARENA_CALLOC(arena, np2,
                                                     sizeof(size_t));
            size_t *ov_cross = (size_t *)ARENA_CALLOC(arena, np2,
                                                      sizeof(size_t));
            /* Self-overlap at turn 0 is two plies of ONE wrap lying on each
             * other -- delamination, which is legitimate and must not be
             * separated.  Only turn != 0 is a piece folded across wraps. */
            size_t *ov_delam = (size_t *)ARENA_CALLOC(arena, np2,
                                                      sizeof(size_t));
            size_t *pfaces = (size_t *)ARENA_CALLOC(arena, np2,
                                                    sizeof(size_t));
            double *pu_lo = (double *)ARENA_ALLOC(arena, np2 * sizeof(double));
            double *pu_hi = (double *)ARENA_ALLOC(arena, np2 * sizeof(double));
            for (size_t i = 0; i < np2; i++) { pu_lo[i] = 0.0; pu_hi[i] = 0.0; }
            for (size_t f = 0; f < ps->nf; f++) {
                int32_t s = vertex_slot[ps->faces[f * 3]];
                if (s < 0 || (size_t)s >= np2) continue;
                if (pfaces[s]++ == 0) {
                    pu_lo[s] = pu_hi[s] = shifted[ps->faces[f * 3]];
                }
                for (int k = 0; k < 3; k++) {
                    double x = shifted[ps->faces[f * 3 + (size_t)k]];
                    if (x < pu_lo[s]) pu_lo[s] = x;
                    if (x > pu_hi[s]) pu_hi[s] = x;
                }
            }
            for (size_t i = 0; i < oa_after.npairs; i++) {
                int32_t a = vertex_slot[ps->faces[
                    (size_t)oa_after.pairs[i].face0 * 3]];
                int32_t b = vertex_slot[ps->faces[
                    (size_t)oa_after.pairs[i].face1 * 3]];
                int oka = a >= 0 && (size_t)a < np2;
                int okb = b >= 0 && (size_t)b < np2;
                if (oka && okb && a == b) {
                    if (oa_after.pairs[i].turn == 0) ov_delam[a]++;
                    else ov_self[a]++;
                } else {
                    if (oka) ov_cross[a]++;
                    if (okb) ov_cross[b]++;
                }
            }
            /* rank by total owned, worst first */
            size_t *rank = (size_t *)ARENA_ALLOC(arena, np2 * sizeof(size_t));
            for (size_t i = 0; i < np2; i++) rank[i] = i;
            for (size_t i = 1; i < np2; i++) {
                size_t key = rank[i], j = i;
                while (j > 0 &&
                       ov_self[rank[j - 1]] + ov_cross[rank[j - 1]] +
                           ov_delam[rank[j - 1]] <
                       ov_self[key] + ov_cross[key] + ov_delam[key]) {
                    rank[j] = rank[j - 1];
                    j--;
                }
                rank[j] = key;
            }
            fprintf(stderr,
                "[atlas_strip_scroll] OVERLAP BY PIECE (worst first). SELF = "
                "folded across wraps, no rigid move can fix it. DELAM = self "
                "at turn 0, two plies of ONE wrap: legitimate, leave it.\n");
            size_t show = np2 < 12 ? np2 : 12;
            size_t tot_self = 0, tot_cross = 0, tot_delam = 0;
            for (size_t s = 0; s < np2; s++) {
                tot_self += ov_self[s];
                tot_cross += ov_cross[s];
                tot_delam += ov_delam[s];
            }
            for (size_t i = 0; i < show; i++) {
                size_t s = rank[i];
                const AtlasSheetComponent *c = &res.components[s];
                fprintf(stderr,
                    "[atlas_strip_scroll]   piece %-3d slot %-2zu %7zu faces, "
                    "u [%9.0f %9.0f], %5.2f turns, self=%-7zu delam=%-7zu "
                    "cross=%-7zu, %zu pairs internal\n",
                    c->mesh_component, s, pfaces[s], pu_lo[s], pu_hi[s],
                    c->turn_span, ov_self[s], ov_delam[s], ov_cross[s],
                    c->pairs_internal);
            }
            fprintf(stderr,
                "[atlas_strip_scroll] OVERLAP TOTALS: self=%zu delam=%zu "
                "cross=%zu (cross counts each pair twice)\n",
                tot_self, tot_delam, tot_cross);
        } else {
            fprintf(stderr,
                "[atlas_strip_scroll] EXACT OVERLAP AUDIT unavailable\n");
        }
    }

    if (cfg->sheet_dump_piece >= 0 &&
        (size_t)cfg->sheet_dump_piece < res.ncomponents) {
        size_t dv = 0, df = 0, dfrag = 0, dtear = 0;
        double tear = 0.5 * 2.0 * 3.14159265358979323846 * 350.0;
        if (write_piece_dump(arena, dir, (int32_t)cfg->sheet_dump_piece, ps,
                             set, vertex_slot, shifted, axial, tear,
                             &dv, &df, &dfrag, &dtear) != 0)
            return -1;
        fprintf(stderr,
            "[atlas_strip_scroll] PIECE DUMP slot %d: %zu verts, %zu faces, "
            "%zu pre-weld fragments, %zu suspect welds -> "
            "piece%02d_{xyz,uv,tears_xyz,tears_uv}.obj\n",
            cfg->sheet_dump_piece, dv, df, dfrag, dtear,
            cfg->sheet_dump_piece);
        double owed = 0.0, held = 0.0, frac = 0.0;
        if (write_piece_arclength_dump(arena, dir,
                                       (int32_t)cfg->sheet_dump_piece, ps,
                                       cal, vertex_slot, base_field, axial,
                                       &owed, &held, &frac) != 0)
            return -1;
        fprintf(stderr,
            "[atlas_strip_scroll] ARCLENGTH slot %d: owes %.0f vox of u, "
            "holds %.0f = %.2f%% of the ramp. Red is where u starts, blue is "
            "where it belongs fully panned out -> "
            "piece%02d_arc_{owed,expected,actual}_{xyz,uv}.obj\n",
            cfg->sheet_dump_piece, owed, held, 100.0 * frac,
            cfg->sheet_dump_piece);
    }

    /* --- dumps ----------------------------------------------------------- */
    size_t nf = 0;
    if (write_component_obj(dir, "sheet_xyz.obj", ps, base_field, axial,
                            vertex_in, rgb, 1, 0, &nf) != 0 ||
        write_component_obj(dir, "sheet_uv_before.obj", ps, base_field, axial,
                            vertex_in, rgb, 1, 1, NULL) != 0 ||
        write_component_obj(dir, "sheet_uv_after.obj", ps, shifted, axial,
                            vertex_in, rgb, 1, 1, NULL) != 0 ||
        write_sheet_compare_obj(dir, "sheet_uv_compare.obj", ps, base_field,
                                shifted, axial, vertex_in, rgb) != 0)
        return -1;
    rep.faces = nf;

    fprintf(stderr,
        "[atlas_strip_scroll] sheet dump: %zu vertices, %zu faces; u span "
        "%.6g -> %.6g; shift range %.6g .. %.6g\n",
        rep.verts, rep.faces, rep.u_span_before, rep.u_span_after,
        rep.shift_min, rep.shift_max);
    fprintf(stderr,
        "[atlas_strip_scroll] SHEET GATES on support %d: %zu own radial "
        "pairs, coincident %zu -> %zu, correctly separated %zu -> %zu\n",
        target, rep.pairs_in_support, rep.coincident_before,
        rep.coincident_after, rep.satisfied_before, rep.satisfied_after);

    if (write_sheet_csv(dir, &res) != 0 ||
        write_sheet_summary(dir, &res, &rep, target, applied) != 0)
        return -1;

    /* Placed sidecars for the moved field, so the RAW texture can be baked
     * through this atlas and looked at as papyrus rather than as colour. */
    {
        ScrollConfig sheet_cfg = *cfg;
        sheet_cfg.out_dir = dir;
        if (export_field_placed(&sheet_cfg, ps, cal, shifted) != 0) return -1;
        if (write_global_mesh_obj(dir, "atlas_flat.obj", ps, cal,
                                  shifted, 1) != 0)
            return -1;
        fprintf(stderr,
            "[atlas_strip_scroll] placed sidecars written to %s/placed\n", dir);
    }
    if (out_field != NULL) *out_field = shifted;
    return 0;
}

static int write_warp_summary(const char *dir, const WarpStageReport *rep)
{
    FILE *fp = as_open(dir, "radial_warp_summary.json", "wb");
    if (fp == NULL) return -1;
    const AtlasWarpStats *w = &rep->warp;
    const AtlasRadialOrderStats *o = &rep->order;
    fprintf(fp,
        "{\n"
        "  \"radial_order\": {\n"
        "    \"pairs\": %zu, \"planes\": %zu, \"median_radius_gap\": %.17g,\n"
        "    \"turn_source_agreement\": {\"radius_agrees\": %zu, "
        "\"radius_disagrees\": %zu},\n"
        "    \"incoming\": {\"coincident\": %zu, \"satisfied\": %zu, "
        "\"reversed\": %zu},\n"
        "    \"rejected\": {\"core\": %zu, \"no_neighbour\": %zu, "
        "\"same_wrap\": %zu, \"ambiguous\": %zu, \"wide\": %zu, "
        "\"delamination\": %zu}\n  },\n",
        o->admitted, o->planes_used, o->median_radius_gap,
        o->turn_agree_radius, o->turn_disagree_radius,
        o->coincident_pairs, o->satisfied_pairs, o->reversed_pairs,
        o->rejected_core, o->rejected_no_neighbour, o->rejected_same_wrap,
        o->rejected_ambiguous, o->rejected_wide, o->rejected_delamination);
    fprintf(fp,
        "  \"warp\": {\n"
        "    \"variables\": %zu, \"rows\": %zu, \"bounds\": %zu, "
        "\"ray_bounds\": %zu,\n"
        "    \"ray_bounds_violated_before\": %zu, "
        "\"ray_bounds_violated_after\": %zu,\n"
        "    \"ray_bounds_disabled\": %zu, \"spacing_rows\": %zu,\n"
        "    \"bootstrap\": {\"moved\": %zu, \"max_shift\": %.17g, "
        "\"anchor_blocked\": %zu, \"cycle_broken\": %zu},\n"
        "    \"qp\": {\"rc\": %d, \"iterations\": %d, \"active\": %d, "
        "\"objective_initial\": %.17g, \"objective_final\": %.17g, "
        "\"min_slack\": %.17g, \"stationarity\": %.17g},\n"
        "    \"gauge_shift\": %.17g,\n"
        "    \"shift\": {\"rms\": %.17g, \"p95\": %.17g, \"max\": %.17g},\n"
        "    \"spacing_residual\": {\"rms_before\": %.17g, \"rms_after\": %.17g,"
        " \"p95_before\": %.17g, \"p95_after\": %.17g},\n"
        "    \"irls\": {\"rounds\": %d, \"sigma\": %.17g, "
        "\"downweighted\": %zu, \"weight_retained\": %.17g},\n"
        "    \"span_before\": %.17g, \"span_after\": %.17g\n  },\n",
        w->nvar, w->nrows, w->nbounds, w->ray_bounds,
        w->ray_bounds_violated_before, w->ray_bounds_violated_after,
        w->ray_bounds_disabled, w->spacing_rows,
        w->bootstrap.moved_variables, w->bootstrap.max_shift,
        w->bootstrap.anchor_blocked, w->bootstrap.cycle_broken,
        w->qp_rc, w->qp_iterations, w->qp_active_final,
        w->objective_initial, w->objective_final, w->min_slack,
        w->stationarity, w->gauge_shift,
        w->shift_rms, w->shift_p95, w->shift_max,
        w->spacing_rms_before, w->spacing_rms_after,
        w->spacing_p95_before, w->spacing_p95_after,
        w->irls_rounds_run, w->irls_sigma, w->spacing_rows_downweighted,
        w->spacing_weight_retained,
        w->span_before, w->span_after);
#define AS_WARP_OVERLAP(name, c, comma)                                       \
    fprintf(fp,                                                               \
        "  \"%s\": {\"total\": %zu, \"whitelisted_delamination\": %zu, "       \
        "\"unwhitelisted\": %zu, \"cross_component\": %zu, "                   \
        "\"same_component\": %zu, \"turn_nonzero\": %zu}%s\n",                 \
        (name), (c).total, (c).whitelisted, (c).unwhitelisted,                \
        (c).cross_component, (c).same_component, (c).turn_nonzero, (comma))
    AS_WARP_OVERLAP("overlap_before", rep->before, ",");
    AS_WARP_OVERLAP("overlap_after", rep->after, ",");
#undef AS_WARP_OVERLAP
    fprintf(fp,
        "  \"relative_flips\": %zu,\n"
        "  \"seam_rms_before\": %.17g,\n"
        "  \"seam_rms_after\": %.17g,\n"
        "  \"weld_oracle_rc\": %d,\n",
        rep->relative_flips, rep->seam_rms_before, rep->seam_rms_after,
        rep->weld_rc);
#define AS_WARP_WELD(name, s, comma)                                          \
    fprintf(fp,                                                               \
        "  \"%s\": {\"cross_cube_bridge_edges\": %zu, "                        \
        "\"edge_stretch\": {\"median\": %.17g, \"p95\": %.17g, "               \
        "\"max\": %.17g}, \"exploded\": %zu, \"compressed\": %zu, "            \
        "\"face_stretch_p95\": %.17g, \"overlap_bridge_base\": %zu, "          \
        "\"overlap_bridge_bridge\": %zu}%s\n",                                 \
        (name), (s).cross_cube_bridge_edges,                                  \
        (s).best_cross_cube_edge_stretch_median,                              \
        (s).best_cross_cube_edge_stretch_p95,                                 \
        (s).best_cross_cube_edge_stretch_max,                                 \
        (s).exploded_cross_cube_edges, (s).compressed_cross_cube_edges,       \
        (s).best_bridge_face_symmetric_stretch_p95,                           \
        (s).bridge_base_overlap_pairs,                                        \
        (s).bridge_bridge_overlap_pairs, (comma))
    AS_WARP_WELD("weld_before", rep->weld_before, ",");
    AS_WARP_WELD("weld_after", rep->weld_after, "");
#undef AS_WARP_WELD
    fprintf(fp, "}\n");
    fclose(fp);
    return 0;
}

/* Per-pair record so a disputed constraint can be traced back to geometry
 * rather than argued about from aggregates. */
static int write_warp_pairs_csv(const char *dir,
                                const AtlasRadialOrderSet *order,
                                const AtlasCandidateSet *set,
                                const double *before, const double *after)
{
    FILE *fp = as_open(dir, "radial_order_pairs.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "inner,outer,plane,component_inner,component_outer,"
                "same_component,support_inner,support_outer,same_support,"
                "turns,radius_inner,radius_gap,spacing_target,"
                "du_before,du_after,weight\n");
    for (size_t i = 0; i < order->npairs; i++) {
        const AtlasRadialOrderPair *p = &order->pairs[i];
        int32_t lo = order->bounds[i].lo, hi = order->bounds[i].hi;
        int32_t ci = set->sample_ref[p->inner].mesh_component;
        int32_t co = set->sample_ref[p->outer].mesh_component;
        /* The support component is the connected piece of the candidate graph:
         * its gauge is the strip solve's null space, so a pair INSIDE one is
         * evidence no per-component shift could ever act on. */
        int32_t si = set->strokes[set->samples[p->inner].stroke].component;
        int32_t so = set->strokes[set->samples[p->outer].stroke].component;
        fprintf(fp,
                "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                p->inner, p->outer, p->plane, ci, co, ci == co ? 1 : 0,
                si, so, si == so ? 1 : 0,
                p->turns, p->radius_inner, p->radius_gap, p->spacing_target,
                before[hi] - before[lo], after[hi] - after[lo], p->weight);
    }
    fclose(fp);
    return 0;
}

static int write_collapse_component_objs(Arena_T arena,
                                         const char *dir, const PieceSet *ps,
                                         const ScaffoldCalib *cal,
                                         const AtlasCandidateSet *set,
                                         const AtlasRadialOrderSet *order,
                                         const double *field_u,
                                         const double *sample_u,
                                         const double *axial,
                                         int count);

typedef struct {
    int32_t component;
    double value;
} AsRegistrationPrior;

static int compare_registration_prior(const void *left, const void *right)
{
    const AsRegistrationPrior *a = (const AsRegistrationPrior *)left;
    const AsRegistrationPrior *b = (const AsRegistrationPrior *)right;
    if (a->component < b->component) return -1;
    if (a->component > b->component) return 1;
    if (a->value < b->value) return -1;
    if (a->value > b->value) return 1;
    return 0;
}

typedef struct {
    const char *dir;
    const ScrollConfig *cfg;
    const PieceSet *ps;
    const ScaffoldCalib *cal;
    const AtlasCandidateSet *set;
    const double *base_field;
    double *iteration_field;
} AsRegistrationTrace;

static int dump_registration_iteration(
    void *context, const AtlasRegisterProblem *problem,
    const double *absolute_shift,
    const AtlasRegisterIterationStats *iteration)
{
    AsRegistrationTrace *trace = (AsRegistrationTrace *)context;
    if (trace == NULL || problem == NULL || absolute_shift == NULL ||
        iteration == NULL || problem->ncomponents != trace->set->mesh_components)
        return -1;
    for (size_t vertex = 0; vertex < trace->ps->nv; vertex++) {
        int32_t component = trace->set->vertex_mesh_component[vertex];
        trace->iteration_field[vertex] = trace->base_field[vertex];
        if (component >= 0 && (size_t)component < problem->ncomponents)
            trace->iteration_field[vertex] += absolute_shift[component];
    }
    char name[128];
    int n = snprintf(name, sizeof(name), "iteration_%02d_flat.obj",
                     iteration->iteration);
    if (n < 0 || (size_t)n >= sizeof(name) ||
        write_global_mesh_obj(trace->dir, name, trace->ps, trace->cal,
                              trace->iteration_field, 1) != 0)
        return -1;
    n = snprintf(name, sizeof(name), "iteration_%02d.json",
                 iteration->iteration);
    if (n < 0 || (size_t)n >= sizeof(name)) return -1;
    FILE *fp = as_open(trace->dir, name, "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"iteration\": %d,\n"
        "  \"frozen_energy_before\": %.17g,\n"
        "  \"frozen_energy_after\": %.17g,\n"
        "  \"robust_energy_before\": %.17g,\n"
        "  \"robust_energy_after\": %.17g,\n"
        "  \"lambda_min_equilibrated\": %.17g,\n"
        "  \"lambda_max_equilibrated\": %.17g,\n"
        "  \"condition_estimate_equilibrated\": %.17g,\n"
        "  \"raw_diagonal_min\": %.17g,\n"
        "  \"raw_diagonal_max\": %.17g,\n"
        "  \"max_relative_step\": %.17g,\n"
        "  \"downweighted_edges\": %zu,\n"
        "  \"weld_residual_rms\": %.17g,\n"
        "  \"radial_residual_rms\": %.17g\n"
        "}\n",
        iteration->iteration, iteration->frozen_energy_before,
        iteration->frozen_energy_after, iteration->robust_energy_before,
        iteration->robust_energy_after, iteration->lambda_min,
        iteration->lambda_max, iteration->condition_estimate,
        iteration->raw_diagonal_min, iteration->raw_diagonal_max,
        iteration->max_relative_step, iteration->downweighted_edges,
        iteration->weld_residual_rms, iteration->radial_residual_rms);
    if (fclose(fp) != 0) return -1;
    n = snprintf(name, sizeof(name), "registration_irls_%02d",
                 iteration->iteration);
    if (n < 0 || (size_t)n >= sizeof(name)) return -1;
    return write_pipeline_trace_mesh(trace->cfg, trace->ps, trace->cal,
                                     name, trace->iteration_field,
                                     iteration, NULL);
}

static int write_registration_outputs(
    const char *dir, const AtlasRegisterProblem *problem,
    const AtlasXyzWeldShiftConstraint *weld, size_t nweld,
    const double *shift, const AtlasRegisterStats *stats,
    int audit_available, const AtlasXyzWeldAuditStats *audit_before,
    const AtlasXyzWeldAuditStats *audit_after)
{
    FILE *fp = as_open(dir, "component_shifts.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "component,shift,geo_prior,prior_weight\n");
    for (size_t component = 0; component < problem->ncomponents; component++) {
        double prior = problem->prior_shift != NULL
                           ? problem->prior_shift[component] : (double)NAN;
        double weight = problem->prior_weight != NULL
                            ? problem->prior_weight[component] : 0.0;
        fprintf(fp, "%zu,%.17g,", component, shift[component]);
        if (isfinite(prior)) fprintf(fp, "%.17g", prior);
        fprintf(fp, ",%.17g\n", weight);
    }
    if (fclose(fp) != 0) return -1;

    fp = as_open(dir, "weld_constraints.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "constraint,component0,component1,cross_cube_edges,"
                "total_xyz_edge_length,target_shift_median,target_shift_mad\n");
    for (size_t i = 0; i < nweld; i++)
        fprintf(fp, "%zu,%d,%d,%zu,%.17g,%.17g,%.17g\n", i,
                weld[i].chart0, weld[i].chart1, weld[i].cross_cube_edges,
                weld[i].total_xyz_edge_length,
                weld[i].target_shift_median, weld[i].target_shift_mad);
    if (fclose(fp) != 0) return -1;

    fp = as_open(dir, "registration_iterations.csv", "wb");
    if (fp == NULL) return -1;
    fprintf(fp, "iteration,frozen_energy_before,frozen_energy_after,"
                "robust_energy_before,robust_energy_after,lambda_min,lambda_max,"
                "condition,raw_diagonal_min,raw_diagonal_max,max_step,"
                "downweighted,weld_rms,radial_rms\n");
    for (int i = 0; i < stats->irls_rounds_run; i++) {
        const AtlasRegisterIterationStats *it = &stats->iteration[i];
        fprintf(fp,
            "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%zu,%.17g,%.17g\n",
            it->iteration, it->frozen_energy_before, it->frozen_energy_after,
            it->robust_energy_before, it->robust_energy_after,
            it->lambda_min, it->lambda_max, it->condition_estimate,
            it->raw_diagonal_min, it->raw_diagonal_max,
            it->max_relative_step, it->downweighted_edges,
            it->weld_residual_rms, it->radial_residual_rms);
    }
    if (fclose(fp) != 0) return -1;

    fp = as_open(dir, "component_registration_summary.json", "wb");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"operation\": \"rigid_mesh_component_registration\",\n"
        "  \"unknown\": \"one additive u shift per connected input mesh component\",\n"
        "  \"solver\": \"grounded graph Laplacian / symmetric Jacobi / sparse SPD Cholesky\",\n"
        "  \"gauge_policy\": \"one null mode eliminated per evidence island; u_geo applied only after solve\",\n"
        "  \"components\": %zu,\n"
        "  \"evidence_islands\": %zu,\n"
        "  \"variables_after_nullspace_elimination\": %zu,\n"
        "  \"usable_edges\": %zu,\n"
        "  \"rejected_edges\": %zu,\n"
        "  \"weld_edges\": %zu,\n"
        "  \"radial_edges\": %zu,\n"
        "  \"trusted_geo_priors\": %zu,\n"
        "  \"irls_rounds\": %d,\n"
        "  \"downweighted_edges_final\": %zu,\n"
        "  \"robust_energy_initial\": %.17g,\n"
        "  \"robust_energy_final\": %.17g,\n"
        "  \"condition_min_equilibrated\": %.17g,\n"
        "  \"condition_max_equilibrated\": %.17g,\n"
        "  \"raw_diagonal_min\": %.17g,\n"
        "  \"raw_diagonal_max\": %.17g,\n"
        "  \"weld_residual_rms_before\": %.17g,\n"
        "  \"weld_residual_rms_after\": %.17g,\n"
        "  \"radial_residual_rms_before\": %.17g,\n"
        "  \"radial_residual_rms_after\": %.17g,\n"
        "  \"shift_rms\": %.17g,\n"
        "  \"shift_max\": %.17g,\n"
        "  \"audit_available\": %s",
        stats->components, stats->evidence_islands, stats->variables,
        stats->usable_edges, stats->rejected_edges, stats->weld_edges,
        stats->radial_edges, stats->trusted_priors,
        stats->irls_rounds_run, stats->downweighted_edges,
        stats->robust_energy_initial, stats->robust_energy_final,
        stats->condition_min, stats->condition_max,
        stats->raw_diagonal_min, stats->raw_diagonal_max,
        stats->weld_residual_rms_before, stats->weld_residual_rms_after,
        stats->radial_residual_rms_before, stats->radial_residual_rms_after,
        stats->shift_rms, stats->shift_max,
        audit_available ? "true" : "false");
    if (audit_available) {
        fprintf(fp,
            ",\n  \"held_out_xyz_weld\": {\n"
            "    \"cross_cube_edges\": %zu,\n"
            "    \"exploded_before\": %zu,\n"
            "    \"exploded_after\": %zu,\n"
            "    \"stretch_p95_before\": %.17g,\n"
            "    \"stretch_p95_after\": %.17g,\n"
            "    \"stretch_max_before\": %.17g,\n"
            "    \"stretch_max_after\": %.17g\n"
            "  }",
            audit_after->cross_cube_bridge_edges,
            audit_before->exploded_cross_cube_edges,
            audit_after->exploded_cross_cube_edges,
            audit_before->best_cross_cube_edge_stretch_p95,
            audit_after->best_cross_cube_edge_stretch_p95,
            audit_before->best_cross_cube_edge_stretch_max,
            audit_after->best_cross_cube_edge_stretch_max);
    }
    fprintf(fp, "\n}\n");
    return fclose(fp) == 0 ? 0 : -1;
}

static int run_component_registration_stage(
    Arena_T result_arena, const ScrollConfig *cfg, const char *dir,
    const PieceSet *ps,
    const ScaffoldCalib *cal, const AtlasCandidateSet *set,
    const AtlasRadialOrderSet *order, const double *base_sample,
    const double *base_field, const double *axial,
    double **out_sample, double **out_field)
{
    if (result_arena == NULL || cfg == NULL || dir == NULL || ps == NULL ||
        cal == NULL ||
        set == NULL || order == NULL || base_sample == NULL ||
        base_field == NULL || axial == NULL || out_sample == NULL ||
        out_field == NULL || set->mesh_components == 0)
        return -1;
    Arena_T scratch = Arena_new();
    if (scratch == NULL) return -1;
    int rc = -1;
    size_t nmesh = set->mesh_components;
    size_t ns = set->problem.nsamples;

    int32_t *face_chart = (int32_t *)ARENA_ALLOC(
        scratch, (ps->nf ? ps->nf : 1) * sizeof(*face_chart));
    for (size_t face = 0; face < ps->nf; face++) {
        int32_t vertex = ps->faces[face * 3];
        if (vertex < 0 || (size_t)vertex >= ps->nv) goto cleanup;
        face_chart[face] = set->vertex_mesh_component[vertex];
        if (face_chart[face] < 0 || (size_t)face_chart[face] >= nmesh)
            goto cleanup;
    }

    AtlasXyzWeldTopology topology;
    memset(&topology, 0, sizeof(topology));
    AtlasXyzWeldShiftConstraint *weld = NULL;
    size_t nweld = 0;
    AtlasXyzWeldConnectionStats weld_stats;
    memset(&weld_stats, 0, sizeof(weld_stats));
    if (ps->n_cubes >= 2) {
        BpaBridgeGate gate;
        memset(&gate, 0, sizeof(gate));
        gate.umb_y = cal->axis_point[1];
        gate.umb_x = cal->axis_point[2];
        gate.pitch = cal->pitch;
        gate.tol = SEAM_WIND_TOL_DEFAULT_TURNS;
        gate.hard = SEAM_WIND_HARD_TOL_DEFAULT_TURNS;
        if (AtlasXyzWeldTopology_build(
                scratch, ps, 128.0f, 1.5f, 0.0f, 6.0f, &gate,
                &topology) != 0 || topology.nfaces == 0 ||
            AtlasXyzWeldTopology_collect_shift_constraints(
                scratch, ps, &topology, base_field, face_chart, nmesh,
                &weld, &nweld, &weld_stats) != 0)
            goto cleanup;
    }

    if (nweld > SIZE_MAX - order->npairs) goto cleanup;
    size_t edge_capacity = nweld + order->npairs;
    AtlasRegisterEdge *edge = (AtlasRegisterEdge *)ARENA_ALLOC(
        scratch, (edge_capacity ? edge_capacity : 1) * sizeof(*edge));
    size_t nedge = 0;
    for (size_t i = 0; i < nweld; i++) {
        double sigma = 1.4826 * weld[i].target_shift_mad;
        edge[nedge].component0 = weld[i].chart0;
        edge[nedge].component1 = weld[i].chart1;
        edge[nedge].target = weld[i].target_shift_median;
        edge[nedge].weight = (double)weld[i].cross_cube_edges;
        edge[nedge].robust_scale = fmax(8.0, 4.0 * sigma);
        edge[nedge].kind = ATLAS_REGISTER_WELD;
        nedge++;
    }
    size_t radial_intra = 0;
    for (size_t i = 0; i < order->npairs; i++) {
        const AtlasRadialOrderPair *pair = &order->pairs[i];
        if (pair->inner < 0 || pair->outer < 0 ||
            (size_t)pair->inner >= ns || (size_t)pair->outer >= ns)
            goto cleanup;
        int32_t inner = set->sample_ref[pair->inner].mesh_component;
        int32_t outer = set->sample_ref[pair->outer].mesh_component;
        if (inner == outer) {
            radial_intra++;
            continue;
        }
        edge[nedge].component0 = inner;
        edge[nedge].component1 = outer;
        edge[nedge].target = pair->spacing_target -
            (base_sample[pair->outer] - base_sample[pair->inner]);
        edge[nedge].weight = pair->weight;
        edge[nedge].robust_scale = 0.35 * fabs(pair->spacing_target);
        edge[nedge].kind = ATLAS_REGISTER_RADIAL;
        nedge++;
    }

    double *prior_shift = (double *)ARENA_ALLOC(
        scratch, nmesh * sizeof(*prior_shift));
    double *prior_weight = (double *)ARENA_CALLOC(
        scratch, nmesh, sizeof(*prior_weight));
    for (size_t component = 0; component < nmesh; component++)
        prior_shift[component] = (double)NAN;
    AsRegistrationPrior *observation = (AsRegistrationPrior *)ARENA_ALLOC(
        scratch, (ns ? ns : 1) * sizeof(*observation));
    size_t nobservation = 0;
    if (set->geo_u != NULL) {
        for (size_t sample = 0; sample < ns; sample++) {
            int32_t component = set->sample_ref[sample].mesh_component;
            if (component < 0 || (size_t)component >= nmesh ||
                !isfinite(set->geo_u[sample]))
                continue;
            observation[nobservation].component = component;
            observation[nobservation].value =
                set->geo_u[sample] - base_sample[sample];
            nobservation++;
        }
    }
    qsort(observation, nobservation, sizeof(*observation),
          compare_registration_prior);
    for (size_t first = 0; first < nobservation; ) {
        size_t last = first + 1;
        while (last < nobservation &&
               observation[last].component == observation[first].component)
            last++;
        size_t count = last - first;
        double median = observation[first + count / 2].value;
        if ((count & 1u) == 0)
            median = 0.5 * (median +
                            observation[first + count / 2 - 1].value);
        prior_shift[observation[first].component] = median;
        prior_weight[observation[first].component] = 1.0;
        first = last;
    }

    AtlasRegisterProblem problem;
    problem.ncomponents = nmesh;
    problem.edges = edge;
    problem.nedges = nedge;
    problem.prior_shift = prior_shift;
    problem.prior_weight = prior_weight;
    double *shift = (double *)ARENA_ALLOC(scratch, nmesh * sizeof(*shift));
    double *iteration_field = (double *)ARENA_ALLOC(
        scratch, ps->nv * sizeof(*iteration_field));
    char trace_dir[AS_PATH_CAP];
    if (as_path(trace_dir, dir, "irls") != 0 ||
        write_global_mesh_obj(trace_dir, "iteration_00_input_flat.obj", ps,
                              cal, base_field, 1) != 0 ||
        write_pipeline_trace_mesh(cfg, ps, cal,
                                  "registration_irls_00_input", base_field,
                                  NULL, NULL) != 0)
        goto cleanup;
    AsRegistrationTrace trace;
    trace.dir = trace_dir;
    trace.cfg = cfg;
    trace.ps = ps;
    trace.cal = cal;
    trace.set = set;
    trace.base_field = base_field;
    trace.iteration_field = iteration_field;
    AtlasRegisterOptions options;
    AtlasRegisterOptions_default(&options);
    options.iteration_fn = dump_registration_iteration;
    options.iteration_context = &trace;
    AtlasRegisterStats stats;
    memset(&stats, 0, sizeof(stats));
    double start = ves_clock_sec();
    if (AtlasRegister_solve(scratch, &problem, &options, shift, &stats) != 0) {
        fprintf(stderr,
            "atlas_strip_scroll: component registration failed "
            "(solve=%d condition=%d energy=%d trace=%d)\n",
            stats.solve_failed, stats.condition_failed, stats.energy_failed,
            stats.trace_failed);
        goto cleanup;
    }

    AtlasXyzWeldAuditStats audit_before, audit_after;
    memset(&audit_before, 0, sizeof(audit_before));
    memset(&audit_after, 0, sizeof(audit_after));
    int audit_available = 0;
    if (topology.nfaces > 0) {
        double *zero = (double *)ARENA_CALLOC(
            scratch, nmesh, sizeof(*zero));
        Arena_Mark audit_mark = Arena_save(scratch);
        if (AtlasXyzWeldAudit_evaluate(
                scratch, ps, &topology, base_field, axial, face_chart, nmesh,
                zero, &audit_before) == 0) {
            Arena_restore(scratch, audit_mark);
            audit_mark = Arena_save(scratch);
            if (AtlasXyzWeldAudit_evaluate(
                    scratch, ps, &topology, base_field, axial, face_chart,
                    nmesh, shift, &audit_after) == 0)
                audit_available = 1;
            Arena_restore(scratch, audit_mark);
        } else {
            Arena_restore(scratch, audit_mark);
        }
    }

    fprintf(stderr,
        "[atlas_strip_scroll] COMPONENT REGISTRATION %.2fs: %zu mesh "
        "components, %zu islands, %zu variables; weld=%zu radial=%zu "
        "(radial intra=%zu); IRLS=%d downweighted=%zu; robust energy "
        "%.6g -> %.6g; equilibrated condition %.4g .. %.4g\n",
        ves_clock_sec() - start, stats.components, stats.evidence_islands,
        stats.variables, stats.weld_edges, stats.radial_edges, radial_intra,
        stats.irls_rounds_run, stats.downweighted_edges,
        stats.robust_energy_initial, stats.robust_energy_final,
        stats.condition_min, stats.condition_max);
    fprintf(stderr,
        "[atlas_strip_scroll] registration residual rms: weld %.4g -> %.4g, "
        "radial %.4g -> %.4g; shift rms=%.4g max=%.4g; %d IRLS meshes "
        "written under %s\n",
        stats.weld_residual_rms_before, stats.weld_residual_rms_after,
        stats.radial_residual_rms_before, stats.radial_residual_rms_after,
        stats.shift_rms, stats.shift_max, stats.irls_rounds_run, trace_dir);
    if (audit_available)
        fprintf(stderr,
            "[atlas_strip_scroll] HELD-OUT XYZ WELD: exploded edges %zu -> "
            "%zu/%zu, stretch p95 %.4g -> %.4g, max %.4g -> %.4g\n",
            audit_before.exploded_cross_cube_edges,
            audit_after.exploded_cross_cube_edges,
            audit_after.cross_cube_bridge_edges,
            audit_before.best_cross_cube_edge_stretch_p95,
            audit_after.best_cross_cube_edge_stretch_p95,
            audit_before.best_cross_cube_edge_stretch_max,
            audit_after.best_cross_cube_edge_stretch_max);

    if (write_registration_outputs(
            dir, &problem, weld, nweld, shift, &stats, audit_available,
            &audit_before, &audit_after) != 0)
        goto cleanup;

    double *registered_sample = (double *)ARENA_ALLOC(
        result_arena, ns * sizeof(*registered_sample));
    for (size_t sample = 0; sample < ns; sample++) {
        int32_t component = set->sample_ref[sample].mesh_component;
        if (component < 0 || (size_t)component >= nmesh) goto cleanup;
        registered_sample[sample] = base_sample[sample] + shift[component];
    }
    double *registered_field = (double *)ARENA_ALLOC(
        result_arena, ps->nv * sizeof(*registered_field));
    for (size_t vertex = 0; vertex < ps->nv; vertex++) {
        int32_t component = set->vertex_mesh_component[vertex];
        registered_field[vertex] = base_field[vertex];
        if (component >= 0 && (size_t)component < nmesh)
            registered_field[vertex] += shift[component];
    }
    if (write_pipeline_trace_mesh(cfg, ps, cal, "registration_output",
                                  registered_field, NULL, NULL) != 0)
        goto cleanup;
    *out_sample = registered_sample;
    *out_field = registered_field;
    rc = 0;

cleanup:
    Arena_dispose(&scratch);
    return rc;
}

static int run_radial_warp_stage(Arena_T arena, const ScrollConfig *cfg,
                                 const PieceSet *ps, const ScaffoldCalib *cal,
                                 const float *raw_u,
                                 AtlasCandidateSet *set,
                                 const double *base_field,
                                 const double *frozen_membership,
                                 double **out_field)
{
    if (base_field == NULL) return -1;
    char dir[AS_PATH_CAP];
    if (as_path(dir, cfg->out_dir, "radial_warp") != 0) return -1;
    if (write_pipeline_trace_mesh(cfg, ps, cal, "radial_stage_input",
                                  base_field, NULL, NULL) != 0)
        return -1;

    size_t ns = set->problem.nsamples;
    double *axial = (double *)ARENA_ALLOC(arena, ps->nv * sizeof(double));
    float *registered = (float *)ARENA_ALLOC(arena, ps->nv * sizeof(float));
    for (size_t i = 0; i < ps->nv; i++) {
        axial[i] = axis_coordinate(cal, &ps->verts[i * 3]);
        registered[i] = ps->uv[i * 2];
    }

    /* --- overlap state of the incoming field, and the delamination whitelist */
    AtlasOverlapAudit overlap_before;
    if (AtlasOverlapAudit_build(arena, ps->faces, ps->nf, ps->nv, base_field,
                                axial, registered, ps->phi,
                                set->vertex_mesh_component,
                                set->mesh_components, &overlap_before) != 0)
        return -1;
    uint64_t *delam_key = (uint64_t *)ARENA_ALLOC(
        arena, (overlap_before.npairs ? overlap_before.npairs : 1) *
                   sizeof(uint64_t));
    size_t ndelam = 0;
    for (size_t i = 0; i < overlap_before.npairs; i++) {
        if (!boxcut_overlap_pair_is_delamination(cfg, ps, set, axial,
                                                 &overlap_before.pairs[i]))
            continue;
        delam_key[ndelam++] =
            boxcut_face_pair_key(overlap_before.pairs[i].face0,
                                 overlap_before.pairs[i].face1);
    }
    qsort(delam_key, ndelam, sizeof(uint64_t), boxcut_compare_u64);

    /* Component pairs a delamination bundle vouches for, so the ray builder
     * can skip them even where eccentricity spreads the plies radially. */
    uint64_t *delam_component = (uint64_t *)ARENA_ALLOC(
        arena, (set->nbundles ? set->nbundles : 1) * sizeof(uint64_t));
    size_t ndelam_component = 0;
    for (size_t i = 0; i < set->nbundles; i++) {
        if (set->bundles[i].decision != ATLAS_BUNDLE_DELAMINATION) continue;
        int32_t c0 = set->bundles[i].component0;
        int32_t c1 = set->bundles[i].component1;
        uint32_t lo = (uint32_t)(c0 < c1 ? c0 : c1);
        uint32_t hi = (uint32_t)(c0 < c1 ? c1 : c0);
        delam_component[ndelam_component++] =
            ((uint64_t)lo << 32) | (uint64_t)hi;
    }
    qsort(delam_component, ndelam_component, sizeof(uint64_t),
          boxcut_compare_u64);

    /* --- ray construction on the composed field ------------------------- */
    double *base_sample = (double *)ARENA_ALLOC(arena, ns * sizeof(double));
    int32_t *sample_plane = (int32_t *)ARENA_ALLOC(arena,
                                                   ns * sizeof(int32_t));
    int32_t *sample_component = (int32_t *)ARENA_ALLOC(arena,
                                                       ns * sizeof(int32_t));
    int32_t *sample_cube = (int32_t *)ARENA_ALLOC(arena,
                                                  ns * sizeof(int32_t));
    double *sample_phi = (double *)ARENA_ALLOC(arena, ns * sizeof(double));
    for (size_t i = 0; i < ns; i++) {
        const AtlasCandidateSampleRef *ref = &set->sample_ref[i];
        base_sample[i] = warp_sample_from_field(set, base_field, i);
        sample_plane[i] = ref->plane;
        sample_component[i] = ref->mesh_component;
        sample_cube[i] = ref->cube;
        double p0 = (double)ps->phi[ref->mesh_vertex[0]];
        double p1 = (double)ps->phi[ref->mesh_vertex[1]];
        sample_phi[i] = p0 + (p1 - p0) * ref->mesh_t;
    }
    if (write_pipeline_trace_strokes(cfg, set, "radial_sample_input",
                                     base_sample, NULL) != 0)
        return -1;

    double t0 = ves_clock_sec();
    AtlasRadialOrderSet order;
    if (AtlasRadialOrder_build(arena, set->samples, sample_plane,
                               sample_component, sample_phi, sample_cube, ns,
                               base_sample, cal, delam_component,
                               ndelam_component,
                               &cfg->order, &order) != 0) {
        fprintf(stderr, "atlas_strip_scroll: radial order build failed\n");
        return -1;
    }
    size_t pair_cross_cube = 0;
    for (size_t i = 0; i < order.npairs; i++)
        if (set->sample_ref[order.pairs[i].inner].cube !=
            set->sample_ref[order.pairs[i].outer].cube)
            pair_cross_cube++;
    const AtlasRadialOrderStats *os = &order.stats;
    size_t turn_votes = os->turn_agree_radius + os->turn_disagree_radius;
    fprintf(stderr,
        "[atlas_strip_scroll] radial order %.2fs: pairs=%zu planes=%zu "
        "median_gap=%.3f local/geometric turn agreement=%.1f%% | "
        "rejected(core=%zu none=%zu same_wrap=%zu ambiguous=%zu wide=%zu "
        "delam=%zu)\n",
        ves_clock_sec() - t0, order.npairs, os->planes_used,
        os->median_radius_gap,
        turn_votes ? 100.0 * (double)os->turn_agree_radius /
                     (double)turn_votes : 0.0,
        os->rejected_core, os->rejected_no_neighbour, os->rejected_same_wrap,
        os->rejected_ambiguous, os->rejected_wide, os->rejected_delamination);
    fprintf(stderr,
        "[atlas_strip_scroll] radial pairs crossing a cube boundary: %zu/%zu "
        "(%.1f%%) -- these use branch-safe radius/angle differentials, not "
        "the independently gauged cube phi values\n",
        pair_cross_cube, order.npairs,
        order.npairs ? 100.0 * (double)pair_cross_cube / (double)order.npairs
                     : 0.0);
    fprintf(stderr,
        "[atlas_strip_scroll] incoming field vs radial evidence: "
        "coincident=%zu (%.1f%%) satisfied=%zu (%.1f%%) reversed=%zu (%.1f%%)"
        " -- coincident pairs ARE the overlap\n",
        os->coincident_pairs,
        order.npairs ? 100.0 * (double)os->coincident_pairs /
                       (double)order.npairs : 0.0,
        os->satisfied_pairs,
        order.npairs ? 100.0 * (double)os->satisfied_pairs /
                       (double)order.npairs : 0.0,
        os->reversed_pairs,
        order.npairs ? 100.0 * (double)os->reversed_pairs /
                       (double)order.npairs : 0.0);

    /*
     * The local strip/FEM field is already the desired shape.  Its remaining
     * ambiguity is one additive u constant per physically connected input
     * mesh component.  Resolve exactly those constants here; never give the
     * post-process enough freedom to turn valid triangles into independent
     * samples or to shear a physically continuous patch.
     */
    double *registered_sample = NULL;
    double *registered_field = NULL;
    if (run_component_registration_stage(
            arena, cfg, dir, ps, cal, set, &order, base_sample, base_field,
            axial,
            &registered_sample, &registered_field) != 0)
        return -1;

    AtlasOverlapAudit registration_overlap_after;
    if (AtlasOverlapAudit_build(arena, ps->faces, ps->nf, ps->nv,
                                registered_field, axial, registered, ps->phi,
                                set->vertex_mesh_component,
                                set->mesh_components,
                                &registration_overlap_after) != 0)
        return -1;
    WarpOverlapCount before_count, after_count;
    warp_count_overlaps(&overlap_before, delam_key, ndelam, &before_count);
    warp_count_overlaps(&registration_overlap_after, delam_key, ndelam,
                        &after_count);
    size_t relative_flips = warp_count_relative_flips(
        ps, base_field, registered_field, axial);
    double seam_before = warp_seam_rms(arena, ps, set, base_field);
    double seam_after = warp_seam_rms(arena, ps, set, registered_field);
    fprintf(stderr,
        "[atlas_strip_scroll] REGISTRATION GATES: overlap %zu -> %zu total, "
        "%zu -> %zu excluding %zu whitelisted delaminations; relative "
        "triangle flips=%zu; seam rms %.6g -> %.6g\n",
        before_count.total, after_count.total,
        before_count.unwhitelisted, after_count.unwhitelisted, ndelam,
        relative_flips, seam_before, seam_after);
    if (relative_flips != 0) {
        fprintf(stderr,
            "atlas_strip_scroll: rigid registration changed triangle "
            "orientation; mesh-component labeling is inconsistent\n");
        return -1;
    }

    FILE *gate_fp = as_open(dir, "registration_gates.json", "wb");
    if (gate_fp == NULL) return -1;
    fprintf(gate_fp,
        "{\n"
        "  \"overlap_total_before\": %zu,\n"
        "  \"overlap_total_after\": %zu,\n"
        "  \"overlap_unwhitelisted_before\": %zu,\n"
        "  \"overlap_unwhitelisted_after\": %zu,\n"
        "  \"delamination_whitelist_pairs\": %zu,\n"
        "  \"relative_triangle_flips\": %zu,\n"
        "  \"seam_rms_before\": %.17g,\n"
        "  \"seam_rms_after\": %.17g\n"
        "}\n",
        before_count.total, after_count.total,
        before_count.unwhitelisted, after_count.unwhitelisted, ndelam,
        relative_flips, seam_before, seam_after);
    if (fclose(gate_fp) != 0 ||
        write_warp_pairs_csv(dir, &order, set, base_sample,
                             registered_sample) != 0 ||
        write_global_mesh_obj(dir, "atlas_registered_flat.obj", ps, cal,
                              registered_field, 1) != 0 ||
        write_global_mesh_obj(dir, "atlas_registered_world_uv.obj", ps, cal,
                              registered_field, 0) != 0 ||
        write_global_mesh_xyz_obj(dir, "atlas_registered_flat_xyzcolor.obj",
                                  ps, cal, registered_field) != 0)
        return -1;

    ScrollConfig registered_cfg = *cfg;
    registered_cfg.out_dir = dir;
    if (export_field_placed(&registered_cfg, ps, cal, registered_field) != 0)
        return -1;

    /* Keep the legacy sheet reports for fixture comparability, but retire
     * their second placement.  That path independently shifted welded pieces
     * by phi-derived turns after registration and is the operation that made
     * atlas_strip_geowind_v6 look like triangle soup. */
    if (cfg->isolate_support >= 0 || cfg->sheet_all) {
        int32_t *sample_support = (int32_t *)ARENA_ALLOC(
            arena, ns * sizeof(*sample_support));
        for (size_t i = 0; i < ns; i++)
            sample_support[i] =
                set->strokes[set->samples[i].stroke].component;
        ScrollConfig diagnostic_cfg = *cfg;
        diagnostic_cfg.sheet_no_shift = 1;
        diagnostic_cfg.sheet_search = 0;
        double *sheet_field = NULL;
        fprintf(stderr,
            "[atlas_strip_scroll] legacy sheet placement is diagnostic-only; "
            "the rigid registered field will be exported unchanged\n");
        if (run_sheet_split_stage(
                arena, &diagnostic_cfg, ps, cal, set, &order,
                sample_support, sample_phi, registered_sample,
                registered_field, axial, &sheet_field) != 0)
            return -1;
        if (write_pipeline_trace_mesh(cfg, ps, cal,
                                      "legacy_sheet_diagnostic_output",
                                      sheet_field, NULL, NULL) != 0)
            return -1;
        if (out_field != NULL) *out_field = sheet_field;
        return 0;
    }

    if (out_field != NULL) *out_field = registered_field;
    return 0;

}

/*
 * Cut the candidate links that fuse wraps.
 *
 * A cross-section models ONE material ruling seen from several slices, and a
 * continuation joins two fragments of ONE stroke; both are contracts about
 * staying on a single wrap.  Admission checks only 3-D proximity and tangent
 * agreement, so nothing enforced them, and on the 4x5x5 fixture 9.2% of active
 * cross-sections and 88% of continuations join samples more than half a turn
 * apart.  Those rows do not merely fail to help -- their unit-speed and
 * alignment terms actively drag two wraps onto one u, and once that is inside
 * a support component no downstream stage can undo it, because there the strip
 * solve fixes u up to a single constant.  Measured: 80% of the collapsed
 * radial pairs lie inside one support component.
 *
 * The minority side is what gets cut, so a cross-section always keeps at least
 * one live member.
 *
 * WRONG, AND KEPT ONLY AS A RECORD.  Almost every link this flags is a false
 * positive: 11,390 of the 11,667 cross-sections and 4,171 of 4,171
 * continuations it cuts join samples in DIFFERENT CUBES, and phi is registered
 * per cube upstream, so a jump there is a registration turn-off rather than a
 * change of wrap.  Cutting them severs genuine physical continuity across the
 * cube seam, and the held-out weld oracle says so unambiguously: cross-cube
 * edge stretch median 0.93 -> 632, p95 3.0 -> 2089, exploded edges 591 ->
 * 56,544 of 65,143.  Restricted to same-cube links the flag has almost nothing
 * left to cut, which is the real finding: the candidate graph is close to
 * clean and the collapse is not caused by wrap-fusing links.
 *
 * INCOMPLETE -- the flag is experimental and off by default.  Cutting links
 * splits the candidate graph into new support components, but the exact gauge
 * anchors were chosen by AtlasCandidates_select_baseline BEFORE the cut, so a
 * newly separated piece has nothing pinning it and the strip system goes
 * singular (measured: robust solve rc=-3 after 11,995 members and 4,171
 * continuations were cut).  Finishing this means recomputing support
 * components and anchors after the cut, which belongs inside atlas_candidates
 * next to the existing select_* passes rather than here.
 */
static int prune_fused_candidate_links(AtlasCandidateSet *set,
                                       const double *sample_phi,
                                       size_t *out_members,
                                       size_t *out_cross_sections,
                                       size_t *out_continuations)
{
    const double pi = 3.14159265358979323846;
    size_t cut_members = 0, cut_cs = 0, cut_cont = 0;
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->problem.cross_sections[c];
        /* Reference wrap = the member closest to the section's own median phi,
         * so the majority side survives however the members are ordered. */
        double best = 0.0;
        double best_score = -1.0;
        int32_t best_cube = -1;
        for (int32_t j = 0; j < cs->count; j++) {
            size_t mj = cs->first + (size_t)j;
            if (set->member_state[mj] == ATLAS_CANDIDATE_INACTIVE) continue;
            const AtlasStripMember *mm = &set->problem.members[mj];
            int32_t vj = mm->value1 >= 0 && mm->value_t > 0.5
                       ? mm->value1 : mm->value0;
            double agree = 0.0;
            for (int32_t k = 0; k < cs->count; k++) {
                size_t mk = cs->first + (size_t)k;
                if (set->member_state[mk] == ATLAS_CANDIDATE_INACTIVE) continue;
                const AtlasStripMember *mn = &set->problem.members[mk];
                int32_t vk = mn->value1 >= 0 && mn->value_t > 0.5
                           ? mn->value1 : mn->value0;
                if (fabs(sample_phi[vk] - sample_phi[vj]) <= pi) agree += 1.0;
            }
            if (agree > best_score) {
                best_score = agree;
                best = sample_phi[vj];
                best_cube = set->sample_ref[vj].cube;
            }
        }
        if (best_score < 0.0) continue;
        int32_t cube_ref = best_cube;
        size_t cut_here = 0;
        for (int32_t j = 0; j < cs->count; j++) {
            size_t mj = cs->first + (size_t)j;
            if (set->member_state[mj] == ATLAS_CANDIDATE_INACTIVE) continue;
            const AtlasStripMember *mm = &set->problem.members[mj];
            int32_t vj = mm->value1 >= 0 && mm->value_t > 0.5
                       ? mm->value1 : mm->value0;
            if (fabs(sample_phi[vj] - best) <= pi) continue;
            /* Only ever cut within one cube: across a seam the phi jump is
             * registration, not geometry. */
            if (set->sample_ref[vj].cube != cube_ref) continue;
            set->member_state[mj] = ATLAS_CANDIDATE_INACTIVE;
            /* member_state alone is inert on the fresh-solve path: the robust
             * solve reads membership directly and only consults the state
             * array when a frozen membership vector is supplied. */
            set->members[mj].membership = 0.0;
            cut_here++;
        }
        if (cut_here > 0) { cut_cs++; cut_members += cut_here; }
    }
    /*
     * Every cross-section owns a latent intercept q whether or not any of its
     * members survive, and a q that appears in no row leaves a zero column in
     * the Hessian -- a singular system reported only as a reduced-SPD failure.
     *
     * The test has to be the BUILD's own criterion, not member_state: the
     * robust solve multiplies every membership by initial_likelihood (1e-5)
     * before the first pass, so a member is only assembled when
     * membership*base_weight*1e-5 clears the 1e-8 floor -- an effective
     * threshold a thousand times stricter than "state says active".
     */
    const double first_pass_floor = 1e-8 / 1e-5;
    size_t empty_cs = 0;
    for (size_t c = 0; c < set->problem.ncross_sections; c++) {
        const AtlasStripCrossSection *cs = &set->problem.cross_sections[c];
        size_t live = 0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m =
                &set->problem.members[cs->first + (size_t)j];
            if (m->membership * m->base_weight > first_pass_floor) live++;
        }
        if (live == 0) empty_cs++;
    }
    if (empty_cs > 0)
        fprintf(stderr,
            "[atlas_strip_scroll] WARNING: %zu cross-sections have no live "
            "member after pruning; their latent intercepts are unconstrained\n",
            empty_cs);
    for (size_t i = 0; i < set->problem.ncontinuations; i++) {
        if (set->continuation_state[i] == ATLAS_CANDIDATE_INACTIVE) continue;
        const AtlasStripContinuation *cn = &set->problem.continuations[i];
        if (fabs(sample_phi[cn->a] - sample_phi[cn->b]) <= pi) continue;
        if (set->sample_ref[cn->a].cube != set->sample_ref[cn->b].cube)
            continue;                      /* registration, not a wrap change */
        set->continuation_state[i] = ATLAS_CANDIDATE_INACTIVE;
        set->continuations[i].membership = 0.0;
        cut_cont++;
    }
    if (out_members != NULL) *out_members = cut_members;
    if (out_cross_sections != NULL) *out_cross_sections = cut_cs;
    if (out_continuations != NULL) *out_continuations = cut_cont;
    return 0;
}

static void usage(const char *exe)
{
    fprintf(stderr,
        "Usage: %s PLACED_DIR OUT_DIR [options]\n"
        "  --build-only                  stop after geometry candidates\n"
        "  --order-only                  stop after the radial-order audit (fast probe)\n"
        "  --isolate-support N           work on support component N alone: detect its\n"
        "                                sheets, move each out by its own whole number of\n"
        "                                turns, dump sheet_split/sheet_{xyz,uv_before,uv_after}\n"
        "  --sheet-no-shift              with --isolate-support, detect and dump only\n"
        "  --sheet-no-weld               skip the XYZ/BPA merge; treat every fragment as\n"
        "                                its own piece (A/B control)\n"
        "  --sheet-uniform-turns K       move the WHOLE block by K windings instead of\n"
        "                                one per sheet (K may be negative)\n"
        "  --sheet-rigid                 evaluate the move once at the block's median phi\n"
        "                                instead of per vertex (no shear)\n"
        "  --sheet-all                   run the sheet stage over EVERY support component\n"
        "                                at once (pieces from different components then\n"
        "                                compete for the same winds)\n"
        "  --sheet-wind-tol X            cut weld bridge edges whose local winding delta\n"
        "                                |dw| exceeds X (wind_cut's test, applied to the\n"
        "                                weld); 0 = off, 0.7 is wind_cut's default\n"
        "  --sheet-dump-piece N          dump welded piece N (slot index from the overlap\n"
        "                                table) coloured by PRE-WELD fragment, plus the\n"
        "                                suspect welds as line segments\n"
        "  --sheet-search                choose each piece's wind by SEARCH against atlas\n"
        "                                occupancy instead of taking it from phi\n"
        "  --sheet-search-span N         candidate winds are -N..+N (4)\n"
        "  --sheet-search-cell X         occupancy raster cell, vox (4)\n"
        "  --sheet-search-lambda X       collision cells one wind of travel is worth (8)\n"
        "  --sheet-search-pair X         collision cells one unsatisfied radial pair is\n"
        "                                worth (1); 0 disables the radial term\n"
        "  --sheet-min-edge-pairs N      drop component links carrying fewer pairs (2)\n"
        "  --sheet-align-residual X      island alignment called unreliable past X pitches (0.35)\n"
        "  --prune-fused-links           EXPERIMENTAL: cut candidate links joining different\n"
        "                                wraps; needs anchor recompute, currently goes singular\n"
        "  --order-arc-window X          ray angular half-window as arc length (6)\n"
        "  --order-separation X          |dr| >= X*pitch is cross-wrap (0.6)\n"
        "  --order-gap-min X             hard ordering bound, vox (8)\n"
        "  --order-turn-tolerance X      |dr/pitch - m| gate for spacing (0.3)\n"
        "  --order-core-pitches X        skip pairs inside X pitches (2)\n"
        "  --order-max-turns N           drop pairs more than N wraps out (4)\n"
        "  --warp-prior X                proximity-to-incoming weight (0.05)\n"
        "  --warp-spacing X              radial spacing row weight (0.25)\n"
        "  --warp-keep-anchors           keep AtlasStrip's exact per-component gauges\n"
        "  --warp-hard-bounds            enable ordering inequalities (slow; off by default)\n"
        "  --warp-irls N                 robust rounds over the radial evidence (4)\n"
        "  --warp-irls-scale X           sigma = X * residual MAD (3)\n"
        "  --diag-boxcut-rigid           also run the retired rigid chart layout (A/B only)\n"
        "  --diag-component-lift         also run the retired constant per-component lift\n"
        "  --diag-classified-resolve     also re-solve with classified memberships\n"
        "  --boxcut-layout-only          stop after the post-multicut chart layout\n"
        "  --boxcut-component-blocks     diagnostic-only: move retained XYZ mesh components as rigid blocks\n"
        "  --slice-spacing X             axial plane spacing (default 8)\n"
        "  --min-stroke-length X         discard shorter fragments\n"
        "  --sample-spacing X             coarse arclength sample spacing\n"
        "  --match-radius X              cross-slice geometry radius\n"
        "  --match-angle X               cross-slice tangent gate, degrees\n"
        "  --max-slice-stride N          bridge up to N empty planes\n"
        "  --max-candidates N            alternatives per sample\n"
        "  --continuation-radius X       same-plane endpoint gap radius\n"
        "  --continuation-angle X        endpoint tangent gate, degrees\n"
        "  --strip-no-wind-gate          admit cross-wrap relations (legacy behavior;\n"
        "                                the gate drops candidate pairs whose pairwise\n"
        "                                winding delta exceeds the tolerance)\n"
        "  --strip-wind-tol X            |dwind| admission tolerance in turns (0.35);\n"
        "                                also the per-stroke prior trust fraction\n"
        "  --strip-no-geo-gauge          orient/gauge strokes from raw per-cube u\n"
        "                                instead of the radius-anchored winding\n"
        "  --strip-no-wind-cut           do not segment slice chains at wrap\n"
        "                                crossings (snake strokes carry the collapse\n"
        "                                through in-mesh fusions)\n"
        "  --strip-active-set            use the historical MonotoneQp active-set\n"
        "                                solver instead of the anchors-only exact LS\n"
        "                                solve + per-stroke PAVA projection\n"
        "  --strip-lambda-prior X        weak spiral prior weight (1e-3); 0 disables.\n"
        "                                Decides which wind an undetermined fragment\n"
        "                                sits on, never its shape\n"
        "  --lambda-align X              latent ruling weight\n"
        "  --lambda-continuation X       fragment continuation weight\n"
        "  --l1-iterations N             robust L1 rounds\n"
        "  --likelihood-iterations N     likelihood rounds\n"
        "  --likelihood-sigma X          0 = 10*pitch when calibrated (else stroke-span\n"
        "                                derivation); <0 forces the legacy derivation\n"
        "  --field-from SAMPLES.csv      reuse a saved coarse final solution\n"
        "  --rank-pack-labels LABELS.i32 diagnostic-only: pack saved charts by 3-D radius\n"
        "  --boxcut-shifts-from FILE.csv diagnostic-only: evaluate external rigid chart shifts\n"
        "  --boxcut-export-labels LABELS.i32 export saved BoxCutter charts and XYZ/BPA weld OBJs\n"
        "  --gross-initial-residual X    pre-fit sheet-offset gate (default 512)\n"
        "  --initial-coherence-mad X     pre-fit bundle MAD gate (default 16)\n"
        "  --overlap-residual-limit X    post-lift overlap gate (default 4)\n"
        "  --overlap-coherence-mad X     post-lift bundle MAD gate (default 4)\n"
        "  --topology-alternative-fraction X  missed-sheet alternative gate (default .5)\n",
        exe);
}

static int parse_double_arg(int argc, char **argv, int *i, double *out)
{
    if (*i + 1 >= argc) return -1;
    char *end = NULL;
    double value = strtod(argv[++*i], &end);
    if (end == argv[*i] || *end != '\0' || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

static int parse_int_arg(int argc, char **argv, int *i, int *out)
{
    if (*i + 1 >= argc) return -1;
    char *end = NULL;
    long value = strtol(argv[++*i], &end, 10);
    if (end == argv[*i] || *end != '\0' || value < -2147483647L ||
        value > 2147483647L) return -1;
    *out = (int)value;
    return 0;
}

static int parse_config(int argc, char **argv, ScrollConfig *cfg)
{
    memset(cfg, 0, sizeof *cfg);
    AtlasCandidateOptions_default(&cfg->candidate);
    AtlasStripOptions_default(&cfg->strip);
    AtlasStripRobustOptions_default(&cfg->robust);
    MonotoneQpOptions_default(&cfg->qp);
    AtlasCandidateRefineOptions_default(&cfg->refine);
    AtlasRadialOrderOptions_default(&cfg->order);
    AtlasWarpOptions_default(&cfg->warp);
    AtlasSheetSplitOptions_default(&cfg->sheet);
    cfg->collapse_dump_count = 3;
    cfg->isolate_support = -1;
    cfg->sheet_dump_piece = -1;
    cfg->sheet_search_pair = -1.0;
    cfg->strip_wind_gate = 1;
    cfg->strip_wind_tol = 0.35;
    cfg->strip_geo_gauge = 1;
    cfg->strip_wind_cut = 1;
    cfg->strip.lambda_prior = 1e-3;
    if (argc < 3) return -1;
    cfg->placed_dir = argv[1];
    cfg->out_dir = argv[2];
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--build-only") == 0) cfg->build_only = 1;
        else if (strcmp(argv[i], "--boxcut-layout-only") == 0)
            cfg->boxcut_layout_only = 1;
        else if (strcmp(argv[i], "--boxcut-component-blocks") == 0)
            cfg->boxcut_component_blocks = 1;
        else if (strcmp(argv[i], "--order-only") == 0)
            cfg->order_only = 1;
        else if (strcmp(argv[i], "--prune-fused-links") == 0)
            cfg->prune_fused_links = 1;
        else if (strcmp(argv[i], "--warp-per-sample") == 0)
            cfg->warp_per_sample = 1;
        else if (strcmp(argv[i], "--collapse-dump") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->collapse_dump_count) != 0) return -1;
        }
        else if (strcmp(argv[i], "--isolate-support") == 0) {
            if (parse_int_arg(argc, argv, &i, &cfg->isolate_support) != 0)
                return -1;
        }
        else if (strcmp(argv[i], "--sheet-no-shift") == 0)
            cfg->sheet_no_shift = 1;
        else if (strcmp(argv[i], "--sheet-no-weld") == 0)
            cfg->sheet_no_weld = 1;
        else if (strcmp(argv[i], "--sheet-rigid") == 0)
            cfg->sheet_rigid = 1;
        else if (strcmp(argv[i], "--sheet-all") == 0)
            cfg->sheet_all = 1;
        else if (strcmp(argv[i], "--sheet-wind-tol") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->sheet_wind_tol) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-dump-piece") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->sheet_dump_piece) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-search") == 0)
            cfg->sheet_search = 1;
        else if (strcmp(argv[i], "--sheet-search-span") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->sheet_search_span) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-search-cell") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->sheet_search_cell) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-search-lambda") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->sheet_search_lambda) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-search-pair") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->sheet_search_pair) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-uniform-turns") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->sheet_uniform_turns) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-min-edge-pairs") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->sheet.min_edge_pairs) != 0) return -1;
        }
        else if (strcmp(argv[i], "--sheet-align-residual") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->sheet.align_max_residual) != 0)
                return -1;
        }
        else if (strcmp(argv[i], "--diag-boxcut-rigid") == 0)
            cfg->diag_boxcut_rigid = 1;
        else if (strcmp(argv[i], "--diag-component-lift") == 0)
            cfg->diag_component_lift = 1;
        else if (strcmp(argv[i], "--diag-classified-resolve") == 0)
            cfg->diag_classified_resolve = 1;
        else if (strcmp(argv[i], "--order-arc-window") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->order.arc_window) != 0) return -1;
        } else if (strcmp(argv[i], "--order-separation") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->order.separation_min_ratio) != 0)
                return -1;
        } else if (strcmp(argv[i], "--order-gap-min") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->order.gap_min) != 0) return -1;
        } else if (strcmp(argv[i], "--order-turn-tolerance") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->order.turn_tolerance) != 0) return -1;
        } else if (strcmp(argv[i], "--order-core-pitches") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->order.core_radius_pitches) != 0)
                return -1;
        } else if (strcmp(argv[i], "--order-max-turns") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->order.max_turns) != 0) return -1;
        } else if (strcmp(argv[i], "--warp-prior") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->warp.lambda_prior) != 0) return -1;
        } else if (strcmp(argv[i], "--warp-spacing") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->warp.lambda_spacing) != 0) return -1;
        } else if (strcmp(argv[i], "--warp-keep-anchors") == 0) {
            cfg->warp.keep_anchors = 1;
        } else if (strcmp(argv[i], "--warp-hard-bounds") == 0) {
            cfg->warp.hard_bounds = 1;
        } else if (strcmp(argv[i], "--warp-irls") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->warp.irls_rounds) != 0) return -1;
        } else if (strcmp(argv[i], "--warp-irls-scale") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->warp.irls_scale) != 0) return -1;
        }
        else if (strcmp(argv[i], "--slice-spacing") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.slice_spacing) != 0) return -1;
        } else if (strcmp(argv[i], "--min-stroke-length") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.min_stroke_length) != 0) return -1;
        } else if (strcmp(argv[i], "--sample-spacing") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.sample_spacing) != 0) return -1;
        } else if (strcmp(argv[i], "--match-radius") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.match_radius) != 0) return -1;
            cfg->user_match_radius = 1;
        } else if (strcmp(argv[i], "--match-angle") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.match_angle_deg) != 0) return -1;
        } else if (strcmp(argv[i], "--max-slice-stride") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->candidate.max_slice_stride) != 0) return -1;
        } else if (strcmp(argv[i], "--max-candidates") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->candidate.max_candidates) != 0) return -1;
        } else if (strcmp(argv[i], "--continuation-radius") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.continuation_radius) != 0) return -1;
            cfg->user_continuation_radius = 1;
        } else if (strcmp(argv[i], "--continuation-angle") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->candidate.continuation_angle_deg) != 0) return -1;
        } else if (strcmp(argv[i], "--strip-no-wind-gate") == 0) {
            cfg->strip_wind_gate = 0;
        } else if (strcmp(argv[i], "--strip-wind-tol") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->strip_wind_tol) != 0) return -1;
        } else if (strcmp(argv[i], "--strip-no-geo-gauge") == 0) {
            cfg->strip_geo_gauge = 0;
        } else if (strcmp(argv[i], "--strip-no-wind-cut") == 0) {
            cfg->strip_wind_cut = 0;
        } else if (strcmp(argv[i], "--strip-active-set") == 0) {
            cfg->robust.use_active_set = 1;
        } else if (strcmp(argv[i], "--strip-lambda-prior") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->strip.lambda_prior) != 0) return -1;
        } else if (strcmp(argv[i], "--lambda-align") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->strip.lambda_align) != 0) return -1;
        } else if (strcmp(argv[i], "--lambda-continuation") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->strip.lambda_continuation) != 0) return -1;
        } else if (strcmp(argv[i], "--l1-iterations") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->robust.l1_iterations) != 0) return -1;
        } else if (strcmp(argv[i], "--likelihood-iterations") == 0) {
            if (parse_int_arg(argc, argv, &i,
                              &cfg->robust.likelihood_iterations) != 0) return -1;
        } else if (strcmp(argv[i], "--likelihood-sigma") == 0) {
            if (parse_double_arg(argc, argv, &i,
                                 &cfg->robust.likelihood_sigma) != 0) return -1;
        } else if (strcmp(argv[i], "--field-from") == 0 && i + 1 < argc) {
            cfg->field_from_csv = argv[++i];
        } else if (strcmp(argv[i], "--rank-pack-labels") == 0 &&
                   i + 1 < argc) {
            cfg->rank_pack_labels = argv[++i];
        } else if (strcmp(argv[i], "--boxcut-shifts-from") == 0 &&
                   i + 1 < argc) {
            cfg->boxcut_shifts_from_csv = argv[++i];
        } else if (strcmp(argv[i], "--boxcut-export-labels") == 0 &&
                   i + 1 < argc) {
            cfg->boxcut_export_labels = argv[++i];
        } else if (strcmp(argv[i], "--gross-initial-residual") == 0) {
            if (parse_double_arg(argc, argv, &i,
                    &cfg->refine.gross_initial_residual) != 0) return -1;
        } else if (strcmp(argv[i], "--initial-coherence-mad") == 0) {
            if (parse_double_arg(argc, argv, &i,
                    &cfg->refine.initial_coherence_mad_limit) != 0) return -1;
        } else if (strcmp(argv[i], "--overlap-residual-limit") == 0) {
            if (parse_double_arg(argc, argv, &i,
                    &cfg->refine.overlap_residual_limit) != 0) return -1;
        } else if (strcmp(argv[i], "--overlap-coherence-mad") == 0) {
            if (parse_double_arg(argc, argv, &i,
                    &cfg->refine.overlap_coherence_mad_limit) != 0) return -1;
        } else if (strcmp(argv[i], "--topology-alternative-fraction") == 0) {
            if (parse_double_arg(
                    argc, argv, &i,
                    &cfg->refine.minimum_topology_alternative_fraction) != 0)
                return -1;
        } else if (strcmp(argv[i], "--help") == 0 ||
                   strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            exit(0);
        } else return -1;
    }
    if (!isfinite(cfg->refine.gross_initial_residual) ||
        cfg->refine.gross_initial_residual <= 0.0 ||
        !isfinite(cfg->refine.initial_coherence_mad_limit) ||
        cfg->refine.initial_coherence_mad_limit <= 0.0 ||
        !isfinite(cfg->refine.overlap_residual_limit) ||
        cfg->refine.overlap_residual_limit <= 0.0 ||
        !isfinite(cfg->refine.overlap_coherence_mad_limit) ||
        cfg->refine.overlap_coherence_mad_limit <= 0.0 ||
        !isfinite(cfg->refine.minimum_topology_alternative_fraction) ||
        cfg->refine.minimum_topology_alternative_fraction < 0.0 ||
        cfg->refine.minimum_topology_alternative_fraction > 1.0)
        return -1;
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        int failures = 0;
        failures += AtlasOverlapAudit_selftest();
        failures += AtlasRadialOrder_selftest();
        failures += AtlasRegister_selftest();
        failures += AtlasFieldApply_selftest();
        failures += AtlasPlaceSearch_selftest();
        failures += AtlasSheetSplit_selftest();
        failures += AtlasTurnAdvance_selftest();
        failures += AtlasWarp_selftest();
        failures += AtlasBoxcutLayout_selftest();
        failures += LiftedMulticut_selftest();
        fprintf(stderr, "[atlas_strip_scroll] selftest %s (%d failures)\n",
                failures == 0 ? "PASSED" : "FAILED", failures);
        return failures == 0 ? 0 : 1;
    }
    ScrollConfig cfg;
    if (parse_config(argc, argv, &cfg) != 0) {
        usage(argv[0]);
        return 2;
    }
    char probe[AS_PATH_CAP];
    if (as_path(probe, cfg.out_dir, ".atlas_strip_scroll_probe") != 0 ||
        ves_ensure_parent_dir(probe) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot create %s\n", cfg.out_dir);
        return 1;
    }
    char pipeline_trace_dir[AS_PATH_CAP];
    int pipeline_trace_next_frame = 0;
    if (initialize_pipeline_trace(&cfg, pipeline_trace_dir,
                                  &pipeline_trace_next_frame) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot initialize ordered trace\n");
        return 1;
    }

    Arena_T arena = Arena_new();
    PieceSet ps;
    if (PieceSet_build(arena, cfg.placed_dir, &ps) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot load placed dir %s\n",
                cfg.placed_dir);
        Arena_dispose(&arena);
        return 1;
    }
    ScaffoldCalib cal;
    if (Scaffold_read_calib(cfg.placed_dir, &cal) != 0) {
        fprintf(stderr, "atlas_strip_scroll: missing placed calibration\n");
        Arena_dispose(&arena);
        return 1;
    }
    /* Winding-aware candidate admission runs off the pinned calibration. */
    cfg.candidate.pitch = cal.pitch;
    cfg.candidate.spiral_a = cal.spiral_a;
    cfg.candidate.spiral_b = cal.spiral_b;
    cfg.candidate.sense = cal.sense;
    cfg.candidate.wind_gate = cfg.strip_wind_gate;
    cfg.candidate.wind_tol = cfg.strip_wind_tol;
    cfg.candidate.geo_gauge = cfg.strip_geo_gauge;
    cfg.candidate.wind_cut = cfg.strip_wind_cut;
    if (cal.pitch > 0.0 && cfg.strip_wind_gate) {
        /*
         * The dwind gate, not "radius < pitch", is now the wrap protection.
         * The historical 4-vox radii existed because proximity alone could
         * not tell a same-wrap neighbour from the next wrap in; with the
         * pairwise winding test at admission they can widen past the CVT
         * edge length, restoring the connectivity the wrap-crossing cuts
         * and dropped slivers cost.
         */
        if (!cfg.user_continuation_radius)
            cfg.candidate.continuation_radius = 0.5 * cal.pitch;
        if (!cfg.user_match_radius)
            cfg.candidate.match_radius = 6.0;
    }
    if (cfg.robust.likelihood_sigma == 0.0 && cal.pitch > 0.0)
        cfg.robust.likelihood_sigma = 10.0 * cal.pitch;
    else if (cfg.robust.likelihood_sigma < 0.0)
        cfg.robust.likelihood_sigma = 0.0; /* legacy stroke-span derivation */
    if (cfg.boxcut_export_labels != NULL) {
        int export_rc = run_saved_boxcut_weld_export(
            arena, &cfg, &ps, &cal);
        Arena_dispose(&arena);
        return export_rc == 0 ? 0 : 1;
    }
    float *raw_u = NULL;
    if (load_raw_u(arena, cfg.placed_dir, &ps, &raw_u) != 0) {
        fprintf(stderr, "atlas_strip_scroll: raw uv sidecar mismatch\n");
        Arena_dispose(&arena);
        return 1;
    }
    float *reference_u = (float *)ARENA_ALLOC(
        arena, ps.nv * sizeof(float));
    for (size_t i = 0; i < ps.nv; i++) {
        reference_u[i] = ps.uv[i * 2];
        if (!isfinite(reference_u[i])) {
            fprintf(stderr, "atlas_strip_scroll: non-finite registered u\n");
            Arena_dispose(&arena);
            return 1;
        }
    }
    cfg.pipeline_trace_raw_u = raw_u;
    double *raw_trace_u = (double *)ARENA_ALLOC(
        arena, ps.nv * sizeof(*raw_trace_u));
    double *registered_trace_u = (double *)ARENA_ALLOC(
        arena, ps.nv * sizeof(*registered_trace_u));
    for (size_t i = 0; i < ps.nv; i++) {
        raw_trace_u[i] = raw_u[i];
        registered_trace_u[i] = reference_u[i];
    }
    if (write_pipeline_trace_mesh(&cfg, &ps, &cal, "input_raw_uv",
                                  raw_trace_u, NULL, NULL) != 0 ||
        write_pipeline_trace_mesh(&cfg, &ps, &cal, "input_registered_uv",
                                  registered_trace_u, NULL, NULL) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot write input trace\n");
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr, "[atlas_strip_scroll] mesh: %zu cubes, %zu vertices, "
                    "%zu retained faces; pitch=%.6g\n",
            ps.n_cubes, ps.nv, ps.nf, cal.pitch);

    if (cfg.rank_pack_labels != NULL) {
        int rank_rc = run_saved_boxcut_rank_pack(arena, &cfg, &ps, &cal);
        Arena_dispose(&arena);
        return rank_rc == 0 ? 0 : 1;
    }

    AtlasCandidateSet set;
    double t0 = ves_clock_sec();
    int build_rc = AtlasCandidates_build(
        arena, ps.verts, ps.nv, ps.faces, ps.face_cube, ps.nf, raw_u,
        cal.axis_point, cal.axis_dir, &cfg.candidate, &set);
    if (build_rc != 0) {
        fprintf(stderr, "atlas_strip_scroll: candidate construction failed\n");
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr,
        "[atlas_strip_scroll] candidates %.2fs: mesh_components=%zu "
        "strokes=%zu samples=%zu cross=%zu continuations=%zu "
        "matched=%.1f%% branches=%zu truncated=%zu\n",
        ves_clock_sec() - t0, set.mesh_components, set.stats.strokes,
        set.stats.samples, set.stats.cross_sections,
        set.stats.continuations, 100.0 * set.stats.source_match_fraction,
        set.stats.branch_nodes, set.stats.truncated_ball_queries);
    if (cal.pitch > 0.0) {
        const AtlasCandidateStats *acs = &set.stats;
        fprintf(stderr,
            "[atlas_strip_scroll] CANDIDATE WINDING (GEO) members |dw|: "
            "<0.1:%zu 0.1-0.35:%zu 0.35-0.65:%zu 0.65-1.5:%zu >1.5:%zu; "
            "fusing=%zu rejected=%zu\n",
            acs->member_wind_hist[0], acs->member_wind_hist[1],
            acs->member_wind_hist[2], acs->member_wind_hist[3],
            acs->member_wind_hist[4], acs->member_wind_fusing,
            acs->member_wind_rejected);
        fprintf(stderr,
            "[atlas_strip_scroll] CANDIDATE WINDING (GEO) continuations "
            "|dw|: <0.1:%zu 0.1-0.35:%zu 0.35-0.65:%zu 0.65-1.5:%zu "
            ">1.5:%zu; fusing=%zu rejected=%zu\n",
            acs->continuation_wind_hist[0], acs->continuation_wind_hist[1],
            acs->continuation_wind_hist[2], acs->continuation_wind_hist[3],
            acs->continuation_wind_hist[4], acs->continuation_wind_fusing,
            acs->continuation_wind_rejected);
        fprintf(stderr,
            "[atlas_strip_scroll] STROKE WINDING (GEO): %zu slice chains "
            "cut at %zu wrap crossings; %zu/%zu emitted strokes still span "
            ">=0.5 turn internally (max span %.2f turns); %zu/%zu "
            "prior-untrusted (u_geo drifts off arclength; no prior rows "
            "there)\n",
            acs->chains_wind_split, acs->chain_wind_cuts,
            acs->strokes_wind_bridge, acs->strokes,
            acs->stroke_wind_span_max,
            acs->strokes_prior_untrusted, acs->strokes);
    }

    if (AtlasCandidates_select_baseline(
            arena, cfg.refine.support_floor, &set) != 0) {
        fprintf(stderr, "atlas_strip_scroll: HEAD baseline selection failed\n");
        Arena_dispose(&arena);
        return 1;
    }
    /* Pruning happens inside the warp stage, after the baseline solve -- see
     * the comment there for why doing it before is singular. */

    size_t baseline_support = set.stats.support_components;
    size_t baseline_candidate_links = set.selection.candidate_links;
    fprintf(stderr,
        "[atlas_strip_scroll] HEAD baseline: active=%zu links, "
        "cross_sections=%zu continuations=%zu support=%zu\n",
        set.selection.candidate_links,
        set.problem.ncross_sections, set.problem.ncontinuations,
        baseline_support);
    if (write_pipeline_trace_strokes(&cfg, &set, "candidate_initial",
                                     set.initial_u, NULL) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot write candidate trace\n");
        Arena_dispose(&arena);
        return 1;
    }

    char baseline_dir[AS_PATH_CAP], classified_dir[AS_PATH_CAP];
    if (as_path(baseline_dir, cfg.out_dir, "baseline_head") != 0 ||
        as_path(classified_dir, cfg.out_dir, "classified") != 0) {
        Arena_dispose(&arena);
        return 1;
    }
    ScrollConfig baseline_cfg = cfg;
    baseline_cfg.out_dir = baseline_dir;

    AtlasStageResult baseline_result;
    memset(&baseline_result, 0, sizeof baseline_result);
    if (cfg.field_from_csv != NULL) {
        double *saved = (double *)ARENA_ALLOC(
            arena, set.problem.nsamples * sizeof(double));
        double *saved_membership = (double *)ARENA_ALLOC(
            arena, set.problem.nmembers * sizeof(double));
        char membership_path[AS_PATH_CAP];
        if (load_sample_solution_csv(cfg.field_from_csv,
                                     set.problem.nsamples, saved) != 0 ||
            companion_csv_path(membership_path, cfg.field_from_csv,
                               "cross_sections.csv") != 0 ||
            load_member_membership_csv(membership_path, &set,
                                       saved_membership) != 0) {
            fprintf(stderr, "atlas_strip_scroll: cannot load %s\n",
                    cfg.field_from_csv);
            Arena_dispose(&arena);
            return 1;
        }
        baseline_result.sample_final = saved;
        baseline_result.member_membership = saved_membership;
        if (write_pipeline_trace_strokes(&baseline_cfg, &set,
                                         "coarse_reused_final", saved,
                                         NULL) != 0) {
            fprintf(stderr,
                    "atlas_strip_scroll: cannot write reused solution trace\n");
            Arena_dispose(&arena);
            return 1;
        }
        if (write_strokes_obj(baseline_dir, "final_parameter.obj", &set,
                              saved, 1) != 0 ||
            write_samples_csv(baseline_dir, &ps, &set, saved, saved) != 0 ||
            write_constraints_csv(baseline_dir, &set, saved_membership,
                                  saved, saved) != 0 ||
            write_selection_csv(baseline_dir, &set, saved) != 0 ||
            write_bundles_csv(baseline_dir, &set) != 0 ||
            run_field_stage(arena, &baseline_cfg, &ps, &cal,
                            raw_u, reference_u, &set, saved,
                            &baseline_result.corrected_sample,
                            &baseline_result.field_u) != 0) {
            fprintf(stderr, "atlas_strip_scroll: reused HEAD baseline failed\n");
            Arena_dispose(&arena);
            return 1;
        }
        fprintf(stderr,
                "[atlas_strip_scroll] reused solved HEAD baseline from %s\n",
                cfg.field_from_csv);
    } else {
        if (run_atlas_stage(arena, &baseline_cfg, &ps, &cal,
                            raw_u, reference_u, &set, NULL, NULL, NULL,
                            &baseline_result) != 0) {
            Arena_dispose(&arena);
            return 1;
        }
    }
    if (cfg.build_only) {
        fprintf(stderr, "[atlas_strip_scroll] HEAD baseline build-only=%s\n",
                baseline_dir);
        Arena_dispose(&arena);
        return 0;
    }

    if (AtlasCandidates_select_refined(
            arena, baseline_result.corrected_sample,
            &cfg.refine, &set) != 0) {
        fprintf(stderr, "atlas_strip_scroll: overlap classification failed\n");
        Arena_dispose(&arena);
        return 1;
    }
    size_t classified_support = set.stats.support_components;
    fprintf(stderr,
        "[atlas_strip_scroll] overlap classification: bundles=%zu "
        "separate_sheet=%zu delamination=%zu same_sheet_seam=%zu "
        "inconclusive=%zu "
        "pruned_links=%zu pruned_continuations=%zu "
        "active_cross_sections=%zu support=%zu\n",
        set.selection.refinement_bundles,
        set.selection.refinement_separate_sheet_bundles,
        set.selection.refinement_delamination_bundles,
        set.selection.refinement_sheet_seam_bundles,
        set.selection.refinement_inconclusive_bundles,
        set.selection.refinement_pruned_candidate_links,
        set.selection.refinement_pruned_continuations,
        set.selection.refinement_remaining_active_cross_sections,
        classified_support);

    /* Narrow missed-sheet cuts should remove only redundant false candidates,
     * so the HEAD quotient normally remains connected.  If a cut nevertheless
     * creates a new gauge, keep it at its solved HEAD position rather than
     * silently tiling it by a raw per-cube origin. */
    size_t head_reanchored_components = 0;
    double maximum_anchor_change = 0.0;
    for (size_t i = 0; i < set.problem.nanchors; i++) {
        int32_t var = set.anchors[i].var;
        double change = baseline_result.sample_final[var] -
                        set.anchors[i].value;
        if (fabs(change) > 1e-8) {
            set.anchors[i].value = baseline_result.sample_final[var];
            head_reanchored_components++;
        }
        if (fabs(change) > maximum_anchor_change)
            maximum_anchor_change = fabs(change);
    }
    fprintf(stderr,
            "[atlas_strip_scroll] HEAD gauge preservation: reanchored=%zu/%zu "
            "max_change=%.6g\n",
            head_reanchored_components, set.problem.nanchors,
            maximum_anchor_change);

    float *baseline_reference_u = (float *)ARENA_ALLOC(
        arena, ps.nv * sizeof(float));
    for (size_t i = 0; i < ps.nv; i++)
        baseline_reference_u[i] = (float)baseline_result.field_u[i];
    /*
     * The chain is a linear fold: each stage consumes the field its
     * predecessor produced, and none of them reaches back into
     * baseline_result.  That is not cosmetic -- the previous winding stage was
     * handed the raw HEAD field while the layout it post-processed had been
     * solved on the composed one, so the two overlap-resolution passes never
     * actually composed.
     */
    ScrollConfig classified_cfg = cfg;
    classified_cfg.out_dir = classified_dir;
    AtlasStageResult classified_result;
    memset(&classified_result, 0, sizeof classified_result);
    int classified_rc = 0;
    if (cfg.diag_classified_resolve)
        classified_rc = run_atlas_stage(
            arena, &classified_cfg, &ps, &cal, raw_u, baseline_reference_u,
            &set, baseline_result.sample_final,
            baseline_result.member_membership,
            baseline_result.corrected_sample, &classified_result);

    const double *chain_field = baseline_result.field_u;
    double *lifted_field = NULL;
    int lift_rc = 0;
    if (cfg.diag_component_lift) {
        lift_rc = run_component_lift_stage(
            arena, &cfg, &ps, &cal, &set, baseline_result.corrected_sample,
            chain_field, &lifted_field);
        if (lift_rc == 0) chain_field = lifted_field;
    }

    double *smooth_seam_field = NULL;
    int smooth_seam_rc = run_smooth_seam_stage(
        arena, &cfg, &ps, &cal, &set, chain_field, &smooth_seam_field);
    if (smooth_seam_rc == 0) chain_field = smooth_seam_field;

    double *warp_field = NULL;
    int warp_rc = run_radial_warp_stage(
        arena, &cfg, &ps, &cal, raw_u, &set, chain_field,
        baseline_result.member_membership, &warp_field);
    if (warp_rc == 0) chain_field = warp_field;

    int boxcut_rigid_rc = 0;
    if (cfg.diag_boxcut_rigid)
        boxcut_rigid_rc = run_boxcut_rigid_stage(
            arena, &cfg, &ps, &cal, &set, baseline_result.field_u,
            chain_field);
    if (write_pipeline_trace_mesh(&cfg, &ps, &cal, "pipeline_final",
                                  chain_field, NULL, NULL) != 0) {
        fprintf(stderr, "atlas_strip_scroll: cannot write final trace frame\n");
        Arena_dispose(&arena);
        return 1;
    }

    FILE *summary = as_open(cfg.out_dir, "two_pass_summary.json", "wb");
    if (summary != NULL) {
        fprintf(summary,
            "{\n  \"baseline_dir\": \"baseline_head\",\n"
            "  \"smooth_seam_dir\": \"smooth_seams\",\n"
            "  \"radial_warp_dir\": \"radial_warp\",\n"
            "  \"diagnostic_dirs\": {\"classified\": \"classified\", "
            "\"lifted\": \"lifted_components\", "
            "\"boxcut_rigid\": \"boxcut_rigid\"},\n"
            "  \"baseline_candidate_links\": %zu,\n"
            "  \"baseline_support_components\": %zu,\n"
            "  \"overlap_bundles\": %zu,\n"
            "  \"separate_sheet_bundles\": %zu,\n"
            "  \"delamination_bundles\": %zu,\n"
            "  \"same_sheet_seam_bundles\": %zu,\n"
            "  \"inconclusive_bundles\": %zu,\n"
            "  \"pruned_candidate_links\": %zu,\n"
            "  \"pruned_continuations\": %zu,\n"
            "  \"remaining_active_cross_sections\": %zu,\n"
            "  \"classified_support_components\": %zu,\n"
            "  \"membership_mode\": \"frozen_head\",\n"
            "  \"frozen_membership_rows_preserved\": %zu,\n"
            "  \"frozen_membership_rows_above_floor\": %zu,\n"
            "  \"frozen_membership_rows_zeroed\": %zu,\n"
            "  \"head_reanchored_components\": %zu,\n"
            "  \"maximum_anchor_change\": %.17g,\n"
            "  \"classified_rc\": %d,\n"
            "  \"component_lift_rc\": %d,\n"
            "  \"smooth_seam_rc\": %d,\n"
            "  \"radial_warp_rc\": %d,\n"
            "  \"boxcut_rigid_rc\": %d\n}\n",
            baseline_candidate_links, baseline_support,
            set.selection.refinement_bundles,
            set.selection.refinement_separate_sheet_bundles,
            set.selection.refinement_delamination_bundles,
            set.selection.refinement_sheet_seam_bundles,
            set.selection.refinement_inconclusive_bundles,
            set.selection.refinement_pruned_candidate_links,
            set.selection.refinement_pruned_continuations,
            set.selection.refinement_remaining_active_cross_sections,
            classified_support,
            classified_result.fixed_membership_rows_preserved,
            classified_result.fixed_membership_rows_above_floor,
            classified_result.fixed_membership_rows_zeroed,
            head_reanchored_components,
            maximum_anchor_change, classified_rc, lift_rc, smooth_seam_rc,
            warp_rc, boxcut_rigid_rc);
        fclose(summary);
    }
    fprintf(stderr,
        "[atlas_strip_scroll] baseline/smooth-seam/radial-warp %s; "
        "baseline=%s classified=%s lifted=%s/lifted_components "
        "smooth=%s/smooth_seams warp=%s/radial_warp\n",
        classified_rc == 0 && lift_rc == 0 && smooth_seam_rc == 0 &&
                warp_rc == 0 && boxcut_rigid_rc == 0
            ? "PASSED" : "FAILED",
        baseline_dir, classified_dir, cfg.out_dir, cfg.out_dir, cfg.out_dir);
    Arena_dispose(&arena);
    return classified_rc == 0 && lift_rc == 0 && smooth_seam_rc == 0 &&
           warp_rc == 0 && boxcut_rigid_rc == 0
        ? 0 : 1;
}
