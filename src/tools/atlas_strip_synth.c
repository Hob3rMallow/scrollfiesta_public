/*
 * atlas_strip_synth.c
 *
 * Inspectable reference tests for the two-level whole-scroll atlas:
 *
 *   1. StrokeStrip-style latent cross-section solve on ordered slice curves.
 *   2. Local mesh-differential extension from the resolved sample atlas.
 *
 * Every stage writes human-inspectable OBJ plus CSV/MatrixMarket/JSON data.
 */

#include "../whole/atlas_strip.h"
#include "../whole/atlas_field.h"
#include "../whole/atlas_field_refine.h"
#include "../whole/atlas_candidates.h"
#include "../whole/atlas_component_lift.h"
#include "../whole/atlas_overlap_audit.h"
#include "../whole/atlas_seam_audit.h"
#include "../whole/monotone_qp.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SYNTH_PI
#define SYNTH_PI 3.14159265358979323846
#endif

#define SYNTH_PATH_CAP 1024
#define SYNTH_NS 33
#define SYNTH_NV 5
#define DENSITY_NDENSE 13
#define DENSITY_NCOARSE 5
#define DENSITY_NFALSE 13
#define DENSITY_NSAMPLE (DENSITY_NDENSE + DENSITY_NCOARSE + DENSITY_NFALSE)
#define DENSITY_NCROSS DENSITY_NDENSE
#define DENSITY_NMEMBER (3 * DENSITY_NCROSS)

enum {
    ROLE_PRIMARY = 0,
    ROLE_PEEL = 1,
    ROLE_OTHER = 2
};

typedef struct {
    int32_t path_index;
    int32_t slice;
    int32_t role;
    int32_t resolved;
    double u_true;
    double u_raw;
} SampleMeta;

typedef struct {
    AtlasStripProblem problem;
    AtlasStripSample *samples;
    AtlasStripStroke *strokes;
    AtlasStripMember *members;
    AtlasStripCrossSection *cross;
    MonotoneQpAnchor *anchors;
    SampleMeta *meta;
    int32_t *primary_map;
    int32_t *peel_map;
    int32_t *other_map;
    double curve_p[SYNTH_NS][3];
    double curve_t[SYNTH_NS][3];
    double curve_n[SYNTH_NS][3];
    double curve_s[SYNTH_NS];
    size_t nsamples;
    size_t nstrokes;
    size_t nmembers;
    size_t ncross;
} StripFixture;

typedef struct {
    int32_t a, b;
    double weight;
    double target;
    int32_t patch;
} MeshEdge;

typedef struct {
    int32_t vertex;
    double target;
    double weight;
    int32_t source;
} MeshObservation;

typedef struct {
    double *p;
    double *vcoord;
    double *u0;
    double *utrue;
    int32_t *faces;
    int32_t *patch;
    int32_t *iu;
    int32_t *iv;
    MeshEdge *edges;
    MeshObservation *obs;
    size_t nv, nf, ne, no;
    size_t vertex_cap, face_cap, edge_cap, obs_cap;
} MeshFixture;

typedef struct {
    AtlasStripProblem problem;
    AtlasStripSample samples[DENSITY_NSAMPLE];
    AtlasStripStroke strokes[3];
    AtlasStripMember members[DENSITY_NMEMBER];
    AtlasStripCrossSection cross[DENSITY_NCROSS];
    MonotoneQpAnchor anchor;
    int32_t dense[DENSITY_NDENSE];
    int32_t coarse[DENSITY_NCOARSE];
    int32_t false_neighbor[DENSITY_NFALSE];
    int32_t role[DENSITY_NSAMPLE];
    double truth[DENSITY_NSAMPLE];
    double raw[DENSITY_NSAMPLE];
} DensityFixture;

static const double PALETTE[8][3] = {
    {0.12, 0.55, 0.95},
    {0.95, 0.38, 0.18},
    {0.20, 0.76, 0.42},
    {0.72, 0.35, 0.88},
    {0.94, 0.72, 0.16},
    {0.18, 0.78, 0.78},
    {0.92, 0.32, 0.66},
    {0.55, 0.58, 0.64}
};

static double d3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static double dist3(const double a[3], const double b[3])
{
    double x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
    return sqrt(x * x + y * y + z * z);
}

static void normalize3(double v[3])
{
    double n = sqrt(d3(v, v));
    if (n < 1e-30) {
        v[0] = 1.0; v[1] = 0.0; v[2] = 0.0;
        return;
    }
    v[0] /= n; v[1] /= n; v[2] /= n;
}

static int compare_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static int make_dir(const char *dir)
{
    char probe[SYNTH_PATH_CAP];
    if (snprintf(probe, sizeof probe, "%s/.atlas_strip_probe", dir) < 0)
        return -1;
    return ves_ensure_parent_dir(probe);
}

static int make_path(char path[SYNTH_PATH_CAP], const char *dir,
                     const char *name)
{
    int n = snprintf(path, SYNTH_PATH_CAP, "%s/%s", dir, name);
    return (n < 0 || n >= SYNTH_PATH_CAP) ? -1 : 0;
}

static FILE *open_output(const char *dir, const char *name,
                         char path[SYNTH_PATH_CAP])
{
    if (make_path(path, dir, name) != 0) return NULL;
    if (ves_ensure_parent_dir(path) != 0) return NULL;
    return fopen(path, "wb");
}

static void build_curve(StripFixture *f)
{
    for (int i = 0; i < SYNTH_NS; i++) {
        double t = 0.72 * (double)i;
        f->curve_p[i][0] = t;
        f->curve_p[i][1] = 3.2 * sin(0.34 * t) +
                           0.35 * sin(0.91 * t);
        f->curve_p[i][2] = 0.0;
    }
    f->curve_s[0] = 0.0;
    for (int i = 1; i < SYNTH_NS; i++)
        f->curve_s[i] = f->curve_s[i - 1] +
                        dist3(f->curve_p[i], f->curve_p[i - 1]);
    for (int i = 0; i < SYNTH_NS; i++) {
        int lo = i > 0 ? i - 1 : i;
        int hi = i + 1 < SYNTH_NS ? i + 1 : i;
        f->curve_t[i][0] = f->curve_p[hi][0] - f->curve_p[lo][0];
        f->curve_t[i][1] = f->curve_p[hi][1] - f->curve_p[lo][1];
        f->curve_t[i][2] = 0.0;
        normalize3(f->curve_t[i]);
        f->curve_n[i][0] = -f->curve_t[i][1];
        f->curve_n[i][1] = f->curve_t[i][0];
        f->curve_n[i][2] = 0.0;
    }
}

static int32_t add_stroke(StripFixture *f, int slice, int role,
                          int first_path, int last_path,
                          double normal_offset, double raw_shift,
                          double true_offset, int resolved)
{
    int32_t stroke_id = (int32_t)f->nstrokes;
    AtlasStripStroke *st = &f->strokes[f->nstrokes++];
    st->first = (int32_t)f->nsamples;
    st->count = last_path - first_path + 1;
    st->component = resolved ? 0 : 1;
    st->resolved = resolved;

    for (int j = first_path; j <= last_path; j++) {
        int32_t id = (int32_t)f->nsamples++;
        AtlasStripSample *sample = &f->samples[id];
        SampleMeta *meta = &f->meta[id];
        memset(sample, 0, sizeof *sample);
        sample->p[0] = f->curve_p[j][0] +
                       normal_offset * f->curve_n[j][0];
        sample->p[1] = f->curve_p[j][1] +
                       normal_offset * f->curve_n[j][1];
        sample->p[2] = (double)slice;
        memcpy(sample->tangent, f->curve_t[j], 3 * sizeof(double));
        sample->s = f->curve_s[j] - f->curve_s[first_path];
        sample->stroke = stroke_id;
        sample->ordinal = j - first_path;
        sample->confidence = role == ROLE_PEEL ? 0.45 : 1.0;

        meta->path_index = j;
        meta->slice = slice;
        meta->role = role;
        meta->resolved = resolved;
        meta->u_true = f->curve_s[j] + true_offset;
        meta->u_raw = meta->u_true + raw_shift;

        if (role == ROLE_PRIMARY)
            f->primary_map[slice * SYNTH_NS + j] = id;
        else if (role == ROLE_PEEL)
            f->peel_map[j] = id;
        else
            f->other_map[j] = id;
    }
    return stroke_id;
}

static void member_derivative(const StripFixture *f, int32_t sample_id,
                              int32_t *lo, int32_t *hi, double *length)
{
    const AtlasStripSample *sample = &f->samples[sample_id];
    const AtlasStripStroke *stroke = &f->strokes[sample->stroke];
    if (sample->ordinal > 0) {
        *lo = sample_id - 1;
        *hi = sample_id;
    } else {
        *lo = sample_id;
        *hi = sample_id + 1;
    }
    *length = f->samples[*hi].s - f->samples[*lo].s;
    if (*length <= 0.0 || *hi >= stroke->first + stroke->count)
        *length = 1.0;
}

static void add_member(StripFixture *f, int32_t sample_id,
                       double base_weight, int observation)
{
    AtlasStripMember *m = &f->members[f->nmembers++];
    memset(m, 0, sizeof *m);
    m->value0 = sample_id;
    m->value1 = -1;
    m->value_t = 0.0;
    member_derivative(f, sample_id, &m->deriv_lo, &m->deriv_hi,
                      &m->deriv_length);
    memcpy(m->p, f->samples[sample_id].p, sizeof m->p);
    memcpy(m->tangent, f->samples[sample_id].tangent, sizeof m->tangent);
    m->dual_width = 1.0;
    m->membership = 1.0;
    m->base_weight = base_weight;
    m->observation = observation;
}

static int build_strip_fixture(Arena_T arena, StripFixture *f)
{
    memset(f, 0, sizeof *f);
    f->samples = (AtlasStripSample *)ARENA_ALLOC(
        arena, 256 * sizeof(AtlasStripSample));
    f->strokes = (AtlasStripStroke *)ARENA_ALLOC(
        arena, 16 * sizeof(AtlasStripStroke));
    f->members = (AtlasStripMember *)ARENA_ALLOC(
        arena, 320 * sizeof(AtlasStripMember));
    f->cross = (AtlasStripCrossSection *)ARENA_ALLOC(
        arena, SYNTH_NS * sizeof(AtlasStripCrossSection));
    f->anchors = (MonotoneQpAnchor *)ARENA_ALLOC(
        arena, 4 * sizeof(MonotoneQpAnchor));
    f->meta = (SampleMeta *)ARENA_ALLOC(arena, 256 * sizeof(SampleMeta));
    f->primary_map = (int32_t *)ARENA_ALLOC(
        arena, SYNTH_NV * SYNTH_NS * sizeof(int32_t));
    f->peel_map = (int32_t *)ARENA_ALLOC(
        arena, SYNTH_NS * sizeof(int32_t));
    f->other_map = (int32_t *)ARENA_ALLOC(
        arena, SYNTH_NS * sizeof(int32_t));
    for (int i = 0; i < SYNTH_NV * SYNTH_NS; i++)
        f->primary_map[i] = -1;
    for (int i = 0; i < SYNTH_NS; i++) {
        f->peel_map[i] = -1;
        f->other_map[i] = -1;
    }

    build_curve(f);

    add_stroke(f, 0, ROLE_PRIMARY, 0, SYNTH_NS - 1,
               0.0, 0.0, 0.0, 1);
    add_stroke(f, 1, ROLE_PRIMARY, 0, SYNTH_NS - 1,
               0.0, 18.0, 0.0, 1);
    add_stroke(f, 2, ROLE_PRIMARY, 0, SYNTH_NS - 1,
               0.0, -11.0, 0.0, 1);
    add_stroke(f, 3, ROLE_PRIMARY, 0, 12,
               0.0, 27.0, 0.0, 1);
    add_stroke(f, 3, ROLE_PRIMARY, 19, SYNTH_NS - 1,
               0.0, 27.0, 0.0, 1);
    add_stroke(f, 4, ROLE_PRIMARY, 0, SYNTH_NS - 1,
               0.0, -23.0, 0.0, 1);

    /* A short normally displaced duplicate: same latent material coordinate. */
    add_stroke(f, 2, ROLE_PEEL, 9, 21,
               0.8, 41.0, 0.0, 1);

    /* Nearby but unrelated layer.  It is deliberately absent from C_j. */
    add_stroke(f, 2, ROLE_OTHER, 5, 27,
               7.0, 0.0, 120.0, 0);

    int observation = 0;
    for (int j = 0; j < SYNTH_NS; j++) {
        AtlasStripCrossSection *cs = &f->cross[f->ncross++];
        cs->first = f->nmembers;
        cs->count = 0;
        memcpy(cs->tangent, f->curve_t[j], sizeof cs->tangent);
        cs->weight = 1.0;
        cs->id = j;

        for (int v = 0; v < SYNTH_NV; v++) {
            int32_t id = f->primary_map[v * SYNTH_NS + j];
            if (id < 0) continue;
            double base = (v == 2 && f->peel_map[j] >= 0) ? 0.5 : 1.0;
            add_member(f, id, base, observation++);
            cs->count++;
        }
        if (f->peel_map[j] >= 0) {
            add_member(f, f->peel_map[j], 0.5, observation++);
            cs->count++;
        }
    }

    f->anchors[0].var = f->primary_map[0];
    f->anchors[0].value = 0.0;
    f->anchors[0].component = 0;
    f->anchors[1].var = f->other_map[5];
    f->anchors[1].value = f->meta[f->other_map[5]].u_true;
    f->anchors[1].component = 1;

    f->problem.samples = f->samples;
    f->problem.nsamples = f->nsamples;
    f->problem.strokes = f->strokes;
    f->problem.nstrokes = f->nstrokes;
    f->problem.members = f->members;
    f->problem.nmembers = f->nmembers;
    f->problem.cross_sections = f->cross;
    f->problem.ncross_sections = f->ncross;
    f->problem.anchors = f->anchors;
    f->problem.nanchors = 2;
    return 0;
}

static void role_color(int role, int stroke, double color[3])
{
    if (role == ROLE_PEEL) {
        color[0] = 0.96; color[1] = 0.28; color[2] = 0.72;
    } else if (role == ROLE_OTHER) {
        color[0] = 0.88; color[1] = 0.12; color[2] = 0.10;
    } else {
        const double *p = PALETTE[stroke & 7];
        color[0] = p[0]; color[1] = p[1]; color[2] = p[2];
    }
}

static int write_strokes_world(const char *dir, const char *name,
                               const StripFixture *f)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "# Ordered slice observations. Pink=peel, red=unrelated.\n");
    for (size_t i = 0; i < f->nsamples; i++) {
        double color[3];
        role_color(f->meta[i].role, f->samples[i].stroke, color);
        fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                f->samples[i].p[0], f->samples[i].p[1],
                f->samples[i].p[2], color[0], color[1], color[2]);
    }
    for (size_t s = 0; s < f->nstrokes; s++) {
        const AtlasStripStroke *st = &f->strokes[s];
        fprintf(fp, "g stroke_%zu\nl", s);
        for (int32_t i = 0; i < st->count; i++)
            fprintf(fp, " %d", st->first + i + 1);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}

static int write_parameter_strokes(const char *dir, const char *name,
                                   const StripFixture *f, const double *x)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "# x=u, y=axial slice, z separates duplicate roles visually.\n");
    for (size_t i = 0; i < f->nsamples; i++) {
        double color[3], lift = 0.0;
        role_color(f->meta[i].role, f->samples[i].stroke, color);
        if (f->meta[i].role == ROLE_PEEL) lift = 0.18;
        if (f->meta[i].role == ROLE_OTHER) lift = 0.36;
        fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                x[i], f->samples[i].p[2], lift,
                color[0], color[1], color[2]);
    }
    for (size_t s = 0; s < f->nstrokes; s++) {
        const AtlasStripStroke *st = &f->strokes[s];
        fprintf(fp, "g stroke_%zu\nl", s);
        for (int32_t i = 0; i < st->count; i++)
            fprintf(fp, " %d", st->first + i + 1);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}

static int write_cross_sections_world(const char *dir, const char *name,
                                      const StripFixture *f)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "# Latent cross-sections: each observation connects to a centroid.\n");
    size_t vertex = 0;
    for (size_t c = 0; c < f->ncross; c++) {
        const AtlasStripCrossSection *cs = &f->cross[c];
        double center[3] = {0.0, 0.0, 0.0}, wsum = 0.0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m = &f->members[cs->first + (size_t)j];
            double w = m->membership * m->base_weight;
            center[0] += w * m->p[0];
            center[1] += w * m->p[1];
            center[2] += w * m->p[2];
            wsum += w;
        }
        if (wsum <= 0.0) continue;
        center[0] /= wsum; center[1] /= wsum; center[2] /= wsum;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m = &f->members[cs->first + (size_t)j];
            const SampleMeta *meta = &f->meta[m->value0];
            double color[3];
            role_color(meta->role, f->samples[m->value0].stroke, color);
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                    m->p[0], m->p[1], m->p[2],
                    color[0], color[1], color[2]);
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                    center[0], center[1], center[2],
                    0.18, 0.90, 0.30);
            fprintf(fp, "l %zu %zu\n", vertex + 1, vertex + 2);
            vertex += 2;
        }
    }
    fclose(fp);
    return 0;
}

static int write_cross_sections_parameter(const char *dir, const char *name,
                                          const StripFixture *f,
                                          const AtlasStripOptions *opts,
                                          const double *x)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    size_t vertex = 0;
    for (size_t c = 0; c < f->ncross; c++) {
        const AtlasStripCrossSection *cs = &f->cross[c];
        double center_u = 0.0, center_v = 0.0, wsum = 0.0;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m = &f->members[cs->first + (size_t)j];
            double w = m->membership * m->base_weight;
            center_u += w * AtlasStrip_member_value(m, x);
            center_v += w * m->p[2];
            wsum += w;
        }
        if (wsum <= 0.0) continue;
        center_u /= wsum;
        center_v /= wsum;
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m = &f->members[cs->first + (size_t)j];
            double residual;
            if (opts->mode == ATLAS_STRIP_FINAL)
                residual = AtlasStrip_member_value(m, x) -
                           x[f->nsamples + c] -
                           d3(cs->tangent, m->p);
            else
                residual = AtlasStrip_member_value(m, x) -
                           x[f->nsamples + c];
            double scale = fmin(1.0, fabs(residual) / 2.0);
            double red = scale, green = 1.0 - 0.65 * scale;
            fprintf(fp, "v %.9g %.9g 0.05 %.4f %.4f %.4f\n",
                    AtlasStrip_member_value(m, x), m->p[2],
                    red, green, 0.12);
            fprintf(fp, "v %.9g %.9g 0.05 %.4f %.4f %.4f\n",
                    center_u, center_v, red, green, 0.12);
            fprintf(fp, "l %zu %zu\n", vertex + 1, vertex + 2);
            vertex += 2;
        }
    }
    fclose(fp);
    return 0;
}

static int write_samples_csv(const char *dir, const StripFixture *f,
                             const double *raw, const double *relaxed,
                             const double *final)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "samples.csv", path);
    if (fp == NULL) return -1;
    fprintf(fp, "sample,stroke,ordinal,path_index,slice,role,resolved,"
                "x,y,z,s,u_true,u_raw,u_relaxed,u_final,"
                "error_relaxed,error_final\n");
    for (size_t i = 0; i < f->nsamples; i++) {
        fprintf(fp,
                "%zu,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                i, f->samples[i].stroke, f->samples[i].ordinal,
                f->meta[i].path_index, f->meta[i].slice,
                f->meta[i].role, f->meta[i].resolved,
                f->samples[i].p[0], f->samples[i].p[1],
                f->samples[i].p[2], f->samples[i].s,
                f->meta[i].u_true, raw[i], relaxed[i], final[i],
                relaxed[i] - f->meta[i].u_true,
                final[i] - f->meta[i].u_true);
    }
    fclose(fp);
    return 0;
}

static int write_cross_sections_csv(const char *dir,
                                    const StripFixture *f,
                                    const double *raw,
                                    const double *relaxed,
                                    const double *final)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "cross_sections.csv", path);
    if (fp == NULL) return -1;
    fprintf(fp, "cross,member,value0,value1,value_t,role,membership,"
                "base_weight,u_raw,u_relaxed,u_final,residual_relaxed,"
                "residual_final\n");
    for (size_t c = 0; c < f->ncross; c++) {
        const AtlasStripCrossSection *cs = &f->cross[c];
        for (int32_t j = 0; j < cs->count; j++) {
            const AtlasStripMember *m = &f->members[cs->first + (size_t)j];
            double rr = AtlasStrip_member_value(m, relaxed) -
                        relaxed[f->nsamples + c];
            double rf = AtlasStrip_member_value(m, final) -
                        final[f->nsamples + c] -
                        d3(cs->tangent, m->p);
            fprintf(fp,
                    "%zu,%d,%d,%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g\n",
                    c, j, m->value0, m->value1, m->value_t,
                    f->meta[m->value0].role,
                    m->membership, m->base_weight,
                    AtlasStrip_member_value(m, raw),
                    AtlasStrip_member_value(m, relaxed),
                    AtlasStrip_member_value(m, final),
                    rr, rf);
        }
    }
    fclose(fp);
    return 0;
}

static int write_rows_csv(const char *dir, const char *name,
                          const MonotoneQpProblem *qp)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "row,kind,owner,target,weight,coefficient,var,value\n");
    for (size_t r = 0; r < qp->nrows; r++) {
        const MonotoneQpRow *row = &qp->rows[r];
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &qp->coeff[row->first + (size_t)j];
            fprintf(fp, "%zu,%d,%d,%.17g,%.17g,%d,%d,%.17g\n",
                    r, row->kind, row->owner, row->target, row->weight,
                    j, c->var, c->value);
        }
    }
    fclose(fp);
    return 0;
}

static int write_bounds_csv(const char *dir, const char *name,
                            const MonotoneQpProblem *qp,
                            const double *x)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "bound,stroke,edge,lo,hi,lower,difference,slack\n");
    for (size_t e = 0; e < qp->nbounds; e++) {
        const MonotoneQpBound *b = &qp->bounds[e];
        double difference = x[b->hi] - x[b->lo];
        fprintf(fp, "%zu,%d,%d,%d,%d,%.17g,%.17g,%.17g\n",
                e, b->stroke, b->edge, b->lo, b->hi, b->lower,
                difference, difference - b->lower);
    }
    fclose(fp);
    return 0;
}

static int write_trace_csv(const char *dir, const char *name,
                           const MonotoneQpTrace *trace)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "iteration,active,added,dropped,objective,max_violation,"
                "direction_norm,step,reduced_linear_residual\n");
    for (size_t i = 0; i < trace->count; i++) {
        const MonotoneQpTraceEntry *e = &trace->entry[i];
        fprintf(fp, "%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                e->iteration, e->n_active, e->added_bound,
                e->dropped_bound, e->objective, e->max_violation,
                e->direction_norm, e->step,
                e->reduced_linear_residual);
    }
    fclose(fp);
    return 0;
}

static int write_robust_trace_csv(const char *dir, const char *name,
                                  const AtlasStripRobustTrace *trace)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "iteration,phase,qp_rc,objective,max_u_change,"
                "max_membership_change,residual_rms,max_residual,"
                "membership_min,membership_mean,membership_max,"
                "downweighted_members\n");
    for (size_t i = 0; i < trace->count; i++) {
        const AtlasStripRobustTraceEntry *e = &trace->entry[i];
        fprintf(fp, "%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%zu\n",
                e->iteration, e->phase, e->qp_rc, e->objective,
                e->max_u_change, e->max_membership_change,
                e->residual_rms, e->max_residual,
                e->membership_min, e->membership_mean,
                e->membership_max, e->downweighted_members);
    }
    fclose(fp);
    return 0;
}

static int write_matrix_market(Arena_T arena, const char *dir,
                               const char *name,
                               const MonotoneQpProblem *qp)
{
    if (qp->nvar > 4096) return -1;
    Arena_Mark mark = Arena_save(arena);
    size_t n = qp->nvar;
    double *dense = (double *)ARENA_CALLOC(arena, n * n, sizeof(double));
    for (size_t r = 0; r < qp->nrows; r++) {
        const MonotoneQpRow *row = &qp->rows[r];
        for (int32_t i = 0; i < row->count; i++) {
            const MonotoneQpCoeff *a = &qp->coeff[row->first + (size_t)i];
            for (int32_t j = 0; j < row->count; j++) {
                const MonotoneQpCoeff *b =
                    &qp->coeff[row->first + (size_t)j];
                dense[(size_t)a->var * n + (size_t)b->var] +=
                    row->weight * a->value * b->value;
            }
        }
    }
    size_t nnz = 0;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j <= i; j++)
            if (fabs(dense[i * n + j]) > 1e-15) nnz++;

    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) { Arena_restore(arena, mark); return -1; }
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real symmetric\n");
    fprintf(fp, "%% Full ungauged Hessian; exact anchors are eliminated by solver.\n");
    fprintf(fp, "%zu %zu %zu\n", n, n, nnz);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j <= i; j++)
            if (fabs(dense[i * n + j]) > 1e-15)
                fprintf(fp, "%zu %zu %.17g\n", i + 1, j + 1,
                        dense[i * n + j]);
    fclose(fp);
    Arena_restore(arena, mark);
    return 0;
}

static int write_rhs_csv(const char *dir, const char *name,
                         const MonotoneQpProblem *qp)
{
    double *rhs = (double *)calloc(qp->nvar, sizeof(double));
    if (rhs == NULL) return -1;
    for (size_t r = 0; r < qp->nrows; r++) {
        const MonotoneQpRow *row = &qp->rows[r];
        for (int32_t j = 0; j < row->count; j++) {
            const MonotoneQpCoeff *c = &qp->coeff[row->first + (size_t)j];
            rhs[c->var] += row->weight * c->value * row->target;
        }
    }
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) { free(rhs); return -1; }
    fprintf(fp, "var,rhs\n");
    for (size_t i = 0; i < qp->nvar; i++)
        fprintf(fp, "%zu,%.17g\n", i, rhs[i]);
    fclose(fp);
    free(rhs);
    return 0;
}

static void mesh_init(Arena_T arena, MeshFixture *m)
{
    memset(m, 0, sizeof *m);
    m->vertex_cap = 512;
    m->face_cap = 1024;
    m->edge_cap = 2048;
    m->obs_cap = 1024;
    m->p = (double *)ARENA_ALLOC(arena, m->vertex_cap * 3 * sizeof(double));
    m->vcoord = (double *)ARENA_ALLOC(arena, m->vertex_cap * sizeof(double));
    m->u0 = (double *)ARENA_ALLOC(arena, m->vertex_cap * sizeof(double));
    m->utrue = (double *)ARENA_ALLOC(arena, m->vertex_cap * sizeof(double));
    m->faces = (int32_t *)ARENA_ALLOC(arena, m->face_cap * 3 * sizeof(int32_t));
    m->patch = (int32_t *)ARENA_ALLOC(arena, m->vertex_cap * sizeof(int32_t));
    m->iu = (int32_t *)ARENA_ALLOC(arena, m->vertex_cap * sizeof(int32_t));
    m->iv = (int32_t *)ARENA_ALLOC(arena, m->vertex_cap * sizeof(int32_t));
    m->edges = (MeshEdge *)ARENA_ALLOC(arena, m->edge_cap * sizeof(MeshEdge));
    m->obs = (MeshObservation *)ARENA_ALLOC(
        arena, m->obs_cap * sizeof(MeshObservation));
}

static int mesh_add_edge(MeshFixture *m, int32_t a, int32_t b,
                         int patch)
{
    if (m->ne >= m->edge_cap) return -1;
    MeshEdge *e = &m->edges[m->ne++];
    e->a = a; e->b = b; e->patch = patch;
    e->target = m->u0[b] - m->u0[a];
    double length = dist3(&m->p[(size_t)a * 3], &m->p[(size_t)b * 3]);
    e->weight = 1.0 / fmax(length, 1e-6);
    return 0;
}

static int mesh_add_observation(MeshFixture *m, int32_t vertex,
                                double target, double weight, int source)
{
    if (m->no >= m->obs_cap) return -1;
    MeshObservation *o = &m->obs[m->no++];
    o->vertex = vertex;
    o->target = target;
    o->weight = weight;
    o->source = source;
    return 0;
}

static int mesh_add_patch(MeshFixture *m, const StripFixture *f,
                          const double *coarse_u,
                          int patch_id, int role,
                          int first_path, int last_path,
                          const double *axial, int nrow,
                          double normal_offset, double raw_shift,
                          double true_offset)
{
    int nu = last_path - first_path + 1;
    size_t base = m->nv;
    size_t addv = (size_t)nu * (size_t)nrow;
    size_t addf = (size_t)(nu - 1) * (size_t)(nrow - 1) * 2;
    if (m->nv + addv > m->vertex_cap || m->nf + addf > m->face_cap)
        return -1;

    for (int r = 0; r < nrow; r++) {
        for (int i = 0; i < nu; i++) {
            int j = first_path + i;
            size_t v = base + (size_t)r * (size_t)nu + (size_t)i;
            m->p[v * 3 + 0] = f->curve_p[j][0] +
                              normal_offset * f->curve_n[j][0];
            m->p[v * 3 + 1] = f->curve_p[j][1] +
                              normal_offset * f->curve_n[j][1];
            m->p[v * 3 + 2] = axial[r];
            m->vcoord[v] = axial[r];
            m->u0[v] = f->curve_s[j] + true_offset + raw_shift;
            m->utrue[v] = f->curve_s[j] + true_offset;
            m->patch[v] = patch_id;
            m->iu[v] = j;
            m->iv[v] = r;

            int32_t sample = -1;
            if (role == ROLE_PRIMARY) {
                int slice = (int)lround(axial[r]);
                if (fabs(axial[r] - (double)slice) < 1e-9 &&
                    slice >= 0 && slice < SYNTH_NV)
                    sample = f->primary_map[slice * SYNTH_NS + j];
            } else if (role == ROLE_PEEL) {
                sample = f->peel_map[j];
            } else {
                sample = f->other_map[j];
            }
            if (sample >= 0 &&
                (i == 0 || i == nu - 1 || (j - first_path) % 4 == 0))
                if (mesh_add_observation(m, (int32_t)v, coarse_u[sample],
                                         80.0, sample) != 0)
                    return -1;
        }
    }

    for (int r = 0; r + 1 < nrow; r++) {
        for (int i = 0; i + 1 < nu; i++) {
            int32_t a = (int32_t)(base + (size_t)r * (size_t)nu + (size_t)i);
            int32_t b = a + 1;
            int32_t c = a + nu;
            int32_t d = c + 1;
            m->faces[m->nf * 3 + 0] = a;
            m->faces[m->nf * 3 + 1] = b;
            m->faces[m->nf * 3 + 2] = d;
            m->nf++;
            m->faces[m->nf * 3 + 0] = a;
            m->faces[m->nf * 3 + 1] = d;
            m->faces[m->nf * 3 + 2] = c;
            m->nf++;
        }
    }
    for (int r = 0; r < nrow; r++)
        for (int i = 0; i + 1 < nu; i++) {
            int32_t a = (int32_t)(base + (size_t)r * (size_t)nu + (size_t)i);
            if (mesh_add_edge(m, a, a + 1, patch_id) != 0) return -1;
        }
    for (int r = 0; r + 1 < nrow; r++)
        for (int i = 0; i < nu; i++) {
            int32_t a = (int32_t)(base + (size_t)r * (size_t)nu + (size_t)i);
            if (mesh_add_edge(m, a, a + nu, patch_id) != 0) return -1;
        }
    for (int r = 0; r + 1 < nrow; r++)
        for (int i = 0; i + 1 < nu; i++) {
            int32_t a = (int32_t)(base + (size_t)r * (size_t)nu + (size_t)i);
            if (mesh_add_edge(m, a, a + nu + 1, patch_id) != 0) return -1;
        }
    m->nv += addv;
    return 0;
}

static int build_mesh_fixture(Arena_T arena, const StripFixture *f,
                              const double *coarse_u, MeshFixture *m)
{
    mesh_init(arena, m);
    const double rows0[3] = {0.0, 1.0, 2.0};
    const double rows1[3] = {2.0, 3.0, 4.0};
    const double rows_peel[2] = {1.8, 2.2};
    const double rows_other[2] = {1.5, 2.5};
    if (mesh_add_patch(m, f, coarse_u, 0, ROLE_PRIMARY,
                       0, SYNTH_NS - 1, rows0, 3, 0.0, 0.0, 0.0) != 0 ||
        mesh_add_patch(m, f, coarse_u, 1, ROLE_PRIMARY,
                       0, SYNTH_NS - 1, rows1, 3, 0.0, 18.0, 0.0) != 0 ||
        mesh_add_patch(m, f, coarse_u, 2, ROLE_PEEL,
                       9, 21, rows_peel, 2, 0.8, 41.0, 0.0) != 0 ||
        mesh_add_patch(m, f, coarse_u, 3, ROLE_OTHER,
                       5, 27, rows_other, 2, 7.0, 0.0, 120.0) != 0)
        return -1;
    return 0;
}

static int build_mesh_qp(Arena_T arena, const MeshFixture *m,
                          MonotoneQpProblem *qp)
{
    AtlasFieldObservation *observation = (AtlasFieldObservation *)ARENA_ALLOC(
        arena, m->no * sizeof(AtlasFieldObservation));
    for (size_t o = 0; o < m->no; o++) {
        observation[o].vertex[0] = m->obs[o].vertex;
        observation[o].vertex[1] = m->obs[o].vertex;
        observation[o].vertex[2] = m->obs[o].vertex;
        observation[o].bary[0] = 1.0;
        observation[o].bary[1] = 0.0;
        observation[o].bary[2] = 0.0;
        observation[o].target = m->obs[o].target;
        observation[o].weight = m->obs[o].weight;
        observation[o].source = m->obs[o].source;
    }
    AtlasFieldProblem problem;
    memset(&problem, 0, sizeof problem);
    problem.position = m->p;
    problem.u0 = m->u0;
    problem.nvertex = m->nv;
    problem.triangle = m->faces;
    problem.ntriangle = m->nf;
    problem.observation = observation;
    problem.nobservation = m->no;
    AtlasFieldSystem system;
    if (AtlasField_build(arena, &problem, &system) != 0) return -1;
    *qp = system.qp;
    return 0;
}

static int write_mesh_obj(const char *dir, const char *name,
                          const MeshFixture *m, const double *u,
                          int flat)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    for (size_t i = 0; i < m->nv; i++) {
        const double *color = PALETTE[m->patch[i] & 7];
        double x = flat ? u[i] : m->p[i * 3 + 0];
        double y = flat ? m->vcoord[i] : m->p[i * 3 + 1];
        double z = flat ? 0.03 * (double)m->patch[i] : m->p[i * 3 + 2];
        fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                x, y, z, color[0], color[1], color[2]);
    }
    for (size_t t = 0; t < m->nf; t++)
        fprintf(fp, "f %d %d %d\n",
                m->faces[t * 3 + 0] + 1,
                m->faces[t * 3 + 1] + 1,
                m->faces[t * 3 + 2] + 1);
    fclose(fp);
    return 0;
}

static int write_mesh_obj_patch_limit(const char *dir, const char *name,
                                      const MeshFixture *m, const double *u,
                                      int flat, int max_patch)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    int32_t *remap = (int32_t *)malloc(m->nv * sizeof(int32_t));
    if (remap == NULL) { fclose(fp); return -1; }
    int32_t next = 1;
    for (size_t i = 0; i < m->nv; i++) {
        if (m->patch[i] > max_patch) {
            remap[i] = -1;
            continue;
        }
        remap[i] = next++;
        const double *color = PALETTE[m->patch[i] & 7];
        double x = flat ? u[i] : m->p[i * 3 + 0];
        double y = flat ? m->vcoord[i] : m->p[i * 3 + 1];
        double z = flat ? 0.03 * (double)m->patch[i] : m->p[i * 3 + 2];
        fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                x, y, z, color[0], color[1], color[2]);
    }
    for (size_t t = 0; t < m->nf; t++) {
        int32_t a = m->faces[t * 3 + 0];
        int32_t b = m->faces[t * 3 + 1];
        int32_t c = m->faces[t * 3 + 2];
        if (remap[a] < 0 || remap[b] < 0 || remap[c] < 0) continue;
        fprintf(fp, "f %d %d %d\n", remap[a], remap[b], remap[c]);
    }
    free(remap);
    fclose(fp);
    return 0;
}

static int write_mesh_csv(const char *dir, const MeshFixture *m,
                          const double *u)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "mesh_vertices.csv", path);
    if (fp == NULL) return -1;
    fprintf(fp, "vertex,patch,path_index,row,x,y,z,v,u_true,u_raw,u_final,error\n");
    for (size_t i = 0; i < m->nv; i++)
        fprintf(fp, "%zu,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g\n",
                i, m->patch[i], m->iu[i], m->iv[i],
                m->p[i * 3 + 0], m->p[i * 3 + 1], m->p[i * 3 + 2],
                m->vcoord[i], m->utrue[i], m->u0[i], u[i],
                u[i] - m->utrue[i]);
    fclose(fp);
    return 0;
}

typedef struct {
    double max_resolved_error;
    double max_peel_error;
    double max_unresolved_error;
    double mesh_observation_rms;
    double mesh_local_delta_rms;
    double overlap_raw_rms;
    double overlap_final_rms;
    double min_jacobian;
    size_t flipped_triangles;
} SynthMetrics;

static void measure_sample_errors(const StripFixture *f, const double *u,
                                  SynthMetrics *m)
{
    for (size_t i = 0; i < f->nsamples; i++) {
        double error = fabs(u[i] - f->meta[i].u_true);
        if (f->meta[i].resolved && error > m->max_resolved_error)
            m->max_resolved_error = error;
        if (f->meta[i].role == ROLE_PEEL && error > m->max_peel_error)
            m->max_peel_error = error;
        if (!f->meta[i].resolved && error > m->max_unresolved_error)
            m->max_unresolved_error = error;
    }
}

static double overlap_rms(const MeshFixture *m, const double *u)
{
    double sum = 0.0;
    size_t count = 0;
    for (size_t a = 0; a < m->nv; a++) {
        if (m->patch[a] != 0 || fabs(m->vcoord[a] - 2.0) > 1e-9) continue;
        for (size_t b = 0; b < m->nv; b++) {
            if (m->patch[b] != 1 || fabs(m->vcoord[b] - 2.0) > 1e-9 ||
                m->iu[a] != m->iu[b]) continue;
            double d = u[a] - u[b];
            sum += d * d;
            count++;
            break;
        }
    }
    return count ? sqrt(sum / (double)count) : DBL_MAX;
}

static void measure_mesh(const MeshFixture *m, const double *u,
                          SynthMetrics *out)
{
    double obs2 = 0.0, local2 = 0.0;
    for (size_t i = 0; i < m->no; i++) {
        double r = u[m->obs[i].vertex] - m->obs[i].target;
        obs2 += r * r;
    }
    for (size_t e = 0; e < m->ne; e++) {
        double r = (u[m->edges[e].b] - u[m->edges[e].a]) -
                   m->edges[e].target;
        local2 += r * r;
    }
    out->mesh_observation_rms = m->no ? sqrt(obs2 / (double)m->no) : 0.0;
    out->mesh_local_delta_rms = m->ne ? sqrt(local2 / (double)m->ne) : 0.0;
    out->overlap_raw_rms = overlap_rms(m, m->u0);
    out->overlap_final_rms = overlap_rms(m, u);
    out->min_jacobian = DBL_MAX;
    for (size_t t = 0; t < m->nf; t++) {
        int32_t a = m->faces[t * 3 + 0];
        int32_t b = m->faces[t * 3 + 1];
        int32_t c = m->faces[t * 3 + 2];
        double jac = (u[b] - u[a]) * (m->vcoord[c] - m->vcoord[a]) -
                     (u[c] - u[a]) * (m->vcoord[b] - m->vcoord[a]);
        if (jac < out->min_jacobian) out->min_jacobian = jac;
        if (jac <= 0.0) out->flipped_triangles++;
    }
    if (out->min_jacobian == DBL_MAX) out->min_jacobian = 0.0;
}

static void density_set_sample(DensityFixture *f, int32_t id,
                               int32_t stroke, int32_t ordinal,
                               double s, double z, int role,
                               double truth, double raw)
{
    AtlasStripSample *sample = &f->samples[id];
    memset(sample, 0, sizeof *sample);
    sample->p[0] = s;
    sample->p[2] = z;
    sample->tangent[0] = 1.0;
    sample->s = s;
    sample->stroke = stroke;
    sample->ordinal = ordinal;
    sample->confidence = 1.0;
    f->role[id] = role;
    f->truth[id] = truth;
    f->raw[id] = raw;
}

static void density_set_member(DensityFixture *f, AtlasStripMember *m,
                               int32_t value0, int32_t value1, double value_t,
                               int32_t deriv_lo, int32_t deriv_hi,
                               double px, double pz, int observation)
{
    memset(m, 0, sizeof *m);
    m->value0 = value0;
    m->value1 = value1;
    m->value_t = value_t;
    m->deriv_lo = deriv_lo;
    m->deriv_hi = deriv_hi;
    m->deriv_length = f->samples[deriv_hi].s - f->samples[deriv_lo].s;
    m->p[0] = px;
    m->p[2] = pz;
    m->tangent[0] = 1.0;
    m->dual_width = 1.0;
    m->membership = 1.0;
    m->base_weight = 1.0;
    m->observation = observation;
}

static void build_density_fixture(DensityFixture *f)
{
    memset(f, 0, sizeof *f);
    int32_t next = 0;

    f->strokes[0].first = next;
    f->strokes[0].count = DENSITY_NDENSE;
    f->strokes[0].component = 0;
    f->strokes[0].resolved = 1;
    for (int j = 0; j < DENSITY_NDENSE; j++) {
        int32_t id = next++;
        f->dense[j] = id;
        density_set_sample(f, id, 0, j, (double)j, 0.0, 0,
                           (double)j, (double)j);
    }

    f->strokes[1].first = next;
    f->strokes[1].count = DENSITY_NCOARSE;
    f->strokes[1].component = 0;
    f->strokes[1].resolved = 1;
    for (int j = 0; j < DENSITY_NCOARSE; j++) {
        int32_t id = next++;
        double s = 3.0 * (double)j;
        f->coarse[j] = id;
        density_set_sample(f, id, 1, j, s, 1.0, 1,
                           s, s + 18.0);
    }

    f->strokes[2].first = next;
    f->strokes[2].count = DENSITY_NFALSE;
    f->strokes[2].component = 1;
    f->strokes[2].resolved = 0;
    for (int j = 0; j < DENSITY_NFALSE; j++) {
        int32_t id = next++;
        f->false_neighbor[j] = id;
        density_set_sample(f, id, 2, j, (double)j, 2.0, 2,
                           60.0 + (double)j, 5.0 + (double)j);
    }

    for (int j = 0; j < DENSITY_NCROSS; j++) {
        AtlasStripCrossSection *cs = &f->cross[j];
        cs->first = (size_t)(3 * j);
        cs->count = 3;
        cs->tangent[0] = 1.0;
        cs->weight = 1.0;
        cs->id = j;

        int32_t da = f->dense[j];
        int32_t da_lo = j > 0 ? f->dense[j - 1] : f->dense[0];
        int32_t da_hi = j > 0 ? f->dense[j] : f->dense[1];
        density_set_member(f, &f->members[3 * j], da, -1, 0.0,
                           da_lo, da_hi, (double)j, 0.0, 3 * j);

        int coarse_edge = j / 3;
        if (coarse_edge >= DENSITY_NCOARSE - 1)
            coarse_edge = DENSITY_NCOARSE - 2;
        double t = ((double)j - 3.0 * (double)coarse_edge) / 3.0;
        density_set_member(f, &f->members[3 * j + 1],
                           f->coarse[coarse_edge],
                           f->coarse[coarse_edge + 1], t,
                           f->coarse[coarse_edge],
                           f->coarse[coarse_edge + 1],
                           (double)j, 1.0, 3 * j + 1);

        int k = DENSITY_NFALSE - 1 - j;
        int32_t df = f->false_neighbor[k];
        int32_t df_lo = k > 0 ? f->false_neighbor[k - 1]
                              : f->false_neighbor[0];
        int32_t df_hi = k > 0 ? f->false_neighbor[k]
                              : f->false_neighbor[1];
        density_set_member(f, &f->members[3 * j + 2], df, -1, 0.0,
                           df_lo, df_hi, (double)k, 2.0, 3 * j + 2);
    }

    f->anchor.var = f->dense[0];
    f->anchor.value = 0.0;
    f->anchor.component = 0;
    f->problem.samples = f->samples;
    f->problem.nsamples = DENSITY_NSAMPLE;
    f->problem.strokes = f->strokes;
    f->problem.nstrokes = 3;
    f->problem.members = f->members;
    f->problem.nmembers = DENSITY_NMEMBER;
    f->problem.cross_sections = f->cross;
    f->problem.ncross_sections = DENSITY_NCROSS;
    f->problem.anchors = &f->anchor;
    f->problem.nanchors = 1;
}

static const double *density_role_color(int role)
{
    static const double colors[3][3] = {
        {0.10, 0.56, 0.96},
        {0.18, 0.78, 0.38},
        {0.92, 0.16, 0.12}
    };
    return colors[role < 0 ? 0 : (role > 2 ? 2 : role)];
}

static int write_density_strokes(const char *dir, const char *name,
                                 const DensityFixture *f,
                                 const double *u, int world)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    for (int i = 0; i < DENSITY_NSAMPLE; i++) {
        const double *color = density_role_color(f->role[i]);
        double x = world ? f->samples[i].p[0] : u[i];
        double y = world ? f->samples[i].p[1] : f->samples[i].p[2];
        double z = world ? f->samples[i].p[2] : 0.0;
        fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                x, y, z, color[0], color[1], color[2]);
    }
    for (int s = 0; s < 3; s++) {
        fprintf(fp, "g stroke_%d\nl", s);
        for (int32_t i = 0; i < f->strokes[s].count; i++)
            fprintf(fp, " %d", f->strokes[s].first + i + 1);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}

static int write_density_candidates(const char *dir, const char *name,
                                    const DensityFixture *f,
                                    const double *u,
                                    const double *membership,
                                    int world)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    size_t vertex = 0;
    for (int c = 0; c < DENSITY_NCROSS; c++) {
        const AtlasStripMember *a = &f->members[3 * c];
        for (int j = 1; j < 3; j++) {
            const AtlasStripMember *b = &f->members[3 * c + j];
            double likelihood = membership != NULL ? membership[3 * c + j]
                                                    : 1.0;
            if (likelihood < 0.0) likelihood = 0.0;
            if (likelihood > 1.0) likelihood = 1.0;
            double color[3] = {1.0 - likelihood, likelihood, 0.08};
            double ax = world ? a->p[0] : AtlasStrip_member_value(a, u);
            double ay = world ? a->p[1] : a->p[2];
            double az = world ? a->p[2] : 0.0;
            double bx = world ? b->p[0] : AtlasStrip_member_value(b, u);
            double by = world ? b->p[1] : b->p[2];
            double bz = world ? b->p[2] : 0.0;
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                    ax, ay, az, color[0], color[1], color[2]);
            fprintf(fp, "v %.9g %.9g %.9g %.4f %.4f %.4f\n",
                    bx, by, bz, color[0], color[1], color[2]);
            fprintf(fp, "l %zu %zu\n", vertex + 1, vertex + 2);
            vertex += 2;
        }
    }
    fclose(fp);
    return 0;
}

static int write_density_csv(const char *dir, const DensityFixture *f,
                             const double *relaxed, const double *final,
                             const double *membership)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "samples.csv", path);
    if (fp == NULL) return -1;
    fprintf(fp, "sample,stroke,ordinal,role,s,u_truth,u_raw,u_relaxed,"
                "u_final,error_relaxed,error_final\n");
    for (int i = 0; i < DENSITY_NSAMPLE; i++)
        fprintf(fp, "%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g\n",
                i, f->samples[i].stroke, f->samples[i].ordinal, f->role[i],
                f->samples[i].s, f->truth[i], f->raw[i], relaxed[i], final[i],
                relaxed[i] - f->truth[i], final[i] - f->truth[i]);
    fclose(fp);

    fp = open_output(dir, "members.csv", path);
    if (fp == NULL) return -1;
    fprintf(fp, "cross,member,kind,value0,value1,value_t,deriv_lo,deriv_hi,"
                "u_relaxed,u_final,q_relaxed,q_final,membership\n");
    for (int c = 0; c < DENSITY_NCROSS; c++)
        for (int j = 0; j < 3; j++) {
            int mi = 3 * c + j;
            const AtlasStripMember *m = &f->members[mi];
            fprintf(fp, "%d,%d,%d,%d,%d,%.17g,%d,%d,%.17g,%.17g,"
                        "%.17g,%.17g,%.17g\n",
                    c, j, j, m->value0, m->value1, m->value_t,
                    m->deriv_lo, m->deriv_hi,
                    AtlasStrip_member_value(m, relaxed),
                    AtlasStrip_member_value(m, final),
                    relaxed[DENSITY_NSAMPLE + c],
                    final[DENSITY_NSAMPLE + c], membership[mi]);
        }
    fclose(fp);
    return 0;
}

static int write_density_readme(const char *dir)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "README.txt", path);
    if (fp == NULL) return -1;
    fprintf(fp,
        "Unequal-density robust cross-section test\n"
        "=========================================\n\n"
        "Blue is a unit-spaced reference stroke. Green samples the same strip\n"
        "three times more coarsely and starts translated by +18. Red is an\n"
        "unrelated monotone stroke deliberately inserted into every candidate\n"
        "cross-section in reverse order. Candidate lines are green when their\n"
        "final relaxed adjacency likelihood is high and red when it is low.\n\n"
        "00_observed_world.obj                 Input curves.\n"
        "01_candidate_cross_sections_world.obj Unfiltered candidates.\n"
        "02_raw_parameter.obj                  Independent raw gauges.\n"
        "03_robust_relaxed_parameter.obj       Local/global relaxed solution.\n"
        "04_relaxed_likelihoods.obj            Classified candidate lines.\n"
        "05_final_parameter.obj                Final Eq. 7 solve.\n"
        "06_final_likelihoods.obj              Frozen relaxed memberships.\n");
    fclose(fp);
    return 0;
}

static int run_density_robust_case(Arena_T arena, const char *root)
{
    char dir[SYNTH_PATH_CAP];
    if (snprintf(dir, sizeof dir, "%s/02_barycentric_robust", root) < 0 ||
        make_dir(dir) != 0) return 1;

    DensityFixture fixture;
    build_density_fixture(&fixture);
    const size_t nvar = DENSITY_NSAMPLE + DENSITY_NCROSS;
    double relaxed[DENSITY_NSAMPLE + DENSITY_NCROSS];
    double final[DENSITY_NSAMPLE + DENSITY_NCROSS];
    double membership[DENSITY_NMEMBER];
    for (int i = 0; i < DENSITY_NSAMPLE; i++) relaxed[i] = fixture.raw[i];
    for (int i = DENSITY_NSAMPLE; i < (int)nvar; i++) relaxed[i] = 0.0;

    AtlasStripOptions relaxed_opts;
    AtlasStripOptions_default(&relaxed_opts);
    relaxed_opts.mode = ATLAS_STRIP_RELAXED;
    relaxed_opts.lambda_length = 1.0;
    relaxed_opts.lambda_align = 4.0;
    relaxed_opts.lambda_local = 0.05;
    relaxed_opts.monotone_fraction = 0.5;
    relaxed_opts.membership_floor = 1e-10;

    AtlasStripRobustOptions robust_opts;
    AtlasStripRobustOptions_default(&robust_opts);
    robust_opts.l1_iterations = 4;
    robust_opts.likelihood_iterations = 20;
    robust_opts.likelihood_sigma = 0.4;
    robust_opts.likelihood_floor = 1e-7;
    robust_opts.convergence_tolerance = 1e-11;
    MonotoneQpOptions qp_opts;
    MonotoneQpOptions_default(&qp_opts);
    qp_opts.max_active_iterations = 256;

    AtlasStripRobustTraceEntry robust_entries[32];
    AtlasStripRobustTrace robust_trace = {robust_entries, 32, 0, NULL, NULL};
    AtlasStripRobustStats robust_stats;
    int robust_rc = AtlasStrip_solve_robust(
        arena, &fixture.problem, &relaxed_opts, &robust_opts, &qp_opts,
        relaxed, membership, &robust_trace, &robust_stats);

    AtlasStripMember weighted_members[DENSITY_NMEMBER];
    memcpy(weighted_members, fixture.members, sizeof weighted_members);
    for (int i = 0; i < DENSITY_NMEMBER; i++)
        weighted_members[i].membership = membership[i];
    AtlasStripProblem weighted_problem = fixture.problem;
    weighted_problem.members = weighted_members;

    AtlasStripSystem relaxed_system;
    int relaxed_build_rc = AtlasStrip_build(
        arena, &weighted_problem, &relaxed_opts, &relaxed_system);

    memcpy(final, relaxed, sizeof final);
    AtlasStripOptions final_opts = relaxed_opts;
    final_opts.mode = ATLAS_STRIP_FINAL;
    AtlasStrip_initialize_intercepts(&weighted_problem, &final_opts, final);
    AtlasStripSystem final_system;
    int final_build_rc = AtlasStrip_build(
        arena, &weighted_problem, &final_opts, &final_system);
    MonotoneQpTraceEntry final_entries[64];
    MonotoneQpTrace final_trace = {final_entries, 64, 0};
    MonotoneQpStats final_stats;
    memset(&final_stats, 0, sizeof final_stats);
    int final_rc = final_build_rc == 0
                 ? MonotoneQp_solve(arena, &final_system.qp, &qp_opts,
                                    final, &final_trace, &final_stats)
                 : -1;

    double max_dense_relaxed = 0.0, max_dense_final = 0.0;
    double max_coarse_relaxed = 0.0, max_coarse_final = 0.0;
    double max_bary_relaxed = 0.0, max_bary_final = 0.0;
    double true_min = DBL_MAX, false_values[DENSITY_NFALSE];
    for (int j = 0; j < DENSITY_NDENSE; j++) {
        double er = fabs(relaxed[fixture.dense[j]] - (double)j);
        double ef = fabs(final[fixture.dense[j]] - (double)j);
        if (er > max_dense_relaxed) max_dense_relaxed = er;
        if (ef > max_dense_final) max_dense_final = ef;
        const AtlasStripMember *coarse_member = &fixture.members[3 * j + 1];
        double br = fabs(AtlasStrip_member_value(coarse_member, relaxed) -
                         relaxed[fixture.dense[j]]);
        double bf = fabs(AtlasStrip_member_value(coarse_member, final) -
                         final[fixture.dense[j]]);
        if (br > max_bary_relaxed) max_bary_relaxed = br;
        if (bf > max_bary_final) max_bary_final = bf;
        double ma = membership[3 * j];
        double mb = membership[3 * j + 1];
        if (ma < true_min) true_min = ma;
        if (mb < true_min) true_min = mb;
        false_values[j] = membership[3 * j + 2];
    }
    for (int j = 0; j < DENSITY_NCOARSE; j++) {
        double truth = 3.0 * (double)j;
        double er = fabs(relaxed[fixture.coarse[j]] - truth);
        double ef = fabs(final[fixture.coarse[j]] - truth);
        if (er > max_coarse_relaxed) max_coarse_relaxed = er;
        if (ef > max_coarse_final) max_coarse_final = ef;
    }
    qsort(false_values, DENSITY_NFALSE, sizeof(double), compare_double);
    double false_median = false_values[DENSITY_NFALSE / 2];
    double false_max = false_values[DENSITY_NFALSE - 1];
    int false_high = 0;
    for (int j = 0; j < DENSITY_NFALSE; j++)
        if (false_values[j] >= 0.5) false_high++;

    AtlasStripMetrics final_metrics;
    AtlasStrip_measure(&weighted_problem, &final_opts, final, &final_metrics);
    int passed = robust_rc == 0 && relaxed_build_rc == 0 && final_rc == 0 &&
                 max_dense_relaxed < 1e-3 && max_coarse_relaxed < 1e-3 &&
                 max_bary_relaxed < 1e-3 && max_dense_final < 1e-3 &&
                 max_coarse_final < 1e-3 && max_bary_final < 1e-3 &&
                 true_min > 0.9 && false_median < 0.1 && false_high <= 3 &&
                 final_metrics.min_monotone_ratio >= 0.5 - 1e-10;

    write_density_readme(dir);
    write_density_strokes(dir, "00_observed_world.obj", &fixture, relaxed, 1);
    write_density_candidates(dir, "01_candidate_cross_sections_world.obj",
                             &fixture, relaxed, NULL, 1);
    write_density_strokes(dir, "02_raw_parameter.obj", &fixture,
                          fixture.raw, 0);
    write_density_strokes(dir, "03_robust_relaxed_parameter.obj", &fixture,
                          relaxed, 0);
    write_density_candidates(dir, "04_relaxed_likelihoods.obj", &fixture,
                             relaxed, membership, 0);
    write_density_strokes(dir, "05_final_parameter.obj", &fixture, final, 0);
    write_density_candidates(dir, "06_final_likelihoods.obj", &fixture,
                             final, membership, 0);
    write_density_csv(dir, &fixture, relaxed, final, membership);
    write_robust_trace_csv(dir, "robust_solver_trace.csv", &robust_trace);
    if (relaxed_build_rc == 0) {
        write_rows_csv(dir, "relaxed_rows.csv", &relaxed_system.qp);
        write_bounds_csv(dir, "relaxed_monotonicity.csv",
                         &relaxed_system.qp, relaxed);
        write_matrix_market(arena, dir, "relaxed_H.mtx", &relaxed_system.qp);
        write_rhs_csv(dir, "relaxed_rhs.csv", &relaxed_system.qp);
    }
    if (final_build_rc == 0) {
        write_rows_csv(dir, "final_rows.csv", &final_system.qp);
        write_bounds_csv(dir, "final_monotonicity.csv",
                         &final_system.qp, final);
        write_trace_csv(dir, "final_solver_trace.csv", &final_trace);
        write_matrix_market(arena, dir, "final_H.mtx", &final_system.qp);
        write_rhs_csv(dir, "final_rhs.csv", &final_system.qp);
    }

    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "metrics.json", path);
    if (fp != NULL) {
        fprintf(fp,
            "{\n  \"case\": \"barycentric_robust\",\n"
            "  \"passed\": %s,\n"
            "  \"domain\": {\"dense_samples\": %d, "
            "\"coarse_samples\": %d, \"density_ratio\": 3, "
            "\"false_candidates\": %d},\n"
            "  \"robust\": {\"rc\": %d, \"qp_solves\": %d, "
            "\"initial_solves\": %d, \"l1_solves\": %d, "
            "\"likelihood_solves\": %d, \"sigma\": %.17g, "
            "\"max_u_change_last\": %.17g, "
            "\"max_membership_change_last\": %.17g},\n"
            "  \"relaxed_error\": {\"dense_max\": %.17g, "
            "\"coarse_max\": %.17g, \"barycentric_match_max\": %.17g},\n"
            "  \"memberships\": {\"true_min\": %.17g, "
            "\"false_median\": %.17g, \"false_max\": %.17g, "
            "\"false_ge_0_5\": %d},\n"
            "  \"final\": {\"rc\": %d, \"dense_max\": %.17g, "
            "\"coarse_max\": %.17g, \"barycentric_match_max\": %.17g, "
            "\"min_monotone_ratio\": %.17g}\n}\n",
            passed ? "true" : "false",
            DENSITY_NDENSE, DENSITY_NCOARSE, DENSITY_NFALSE,
            robust_rc, robust_stats.total_qp_solves,
            robust_stats.initial_solves, robust_stats.l1_solves,
            robust_stats.likelihood_solves, robust_stats.likelihood_sigma,
            robust_stats.max_u_change, robust_stats.max_membership_change,
            max_dense_relaxed, max_coarse_relaxed, max_bary_relaxed,
            true_min, false_median, false_max, false_high,
            final_rc, max_dense_final, max_coarse_final, max_bary_final,
            final_metrics.min_monotone_ratio);
        fclose(fp);
    }

    fprintf(stderr,
            "[atlas_strip_synth] barycentric robust: %s "
            "(dense %.3e, coarse %.3e, bary %.3e, "
            "membership true-min %.3f false-med %.3g high %d)\n",
            passed ? "ok" : "FAIL", max_dense_final, max_coarse_final,
            max_bary_final, true_min, false_median, false_high);
    return passed ? 0 : 1;
}

static int run_active_set_reference(Arena_T arena, const char *root)
{
    char dir[SYNTH_PATH_CAP];
    if (snprintf(dir, sizeof dir, "%s/00_active_set_reference", root) < 0 ||
        make_dir(dir) != 0) return 1;
    int fails = 0;
    MonotoneQpOptions opts;
    MonotoneQpOptions_default(&opts);
    opts.max_active_iterations = 64;

    MonotoneQpCoeff coeff[1] = {{1, 1.0}};
    MonotoneQpRow row = {0, 1, -1.0, 1.0, 100, 0};
    MonotoneQpBound bound = {0, 1, 0.5, 0, 0};
    MonotoneQpAnchor anchor = {0, 0.0, 0};
    MonotoneQpProblem qp;
    memset(&qp, 0, sizeof qp);
    qp.nvar = 2; qp.rows = &row; qp.nrows = 1;
    qp.coeff = coeff; qp.ncoeff = 1;
    qp.bounds = &bound; qp.nbounds = 1;
    qp.anchors = &anchor; qp.nanchors = 1;

    MonotoneQpTraceEntry entries_a[32], entries_b[32];
    MonotoneQpTrace trace_a = {entries_a, 32, 0};
    MonotoneQpTrace trace_b = {entries_b, 32, 0};
    MonotoneQpStats stats_a, stats_b;
    double xa[2] = {0.0, 1.0};
    int rc_a = MonotoneQp_solve(arena, &qp, &opts, xa, &trace_a, &stats_a);
    if (rc_a != 0 || fabs(xa[1] - 0.5) > 1e-9 || stats_a.active_final != 1)
        fails++;

    row.target = 2.0;
    double xb[2] = {0.0, 0.5};
    int rc_b = MonotoneQp_solve(arena, &qp, &opts, xb, &trace_b, &stats_b);
    if (rc_b != 0 || fabs(xb[1] - 2.0) > 1e-9 ||
        stats_b.dropped_total < 1 || stats_b.active_final != 0)
        fails++;

    /*
     * A five-variable oracle exercises several simultaneously active chain
     * bounds.  With y_i = x_i - i/2, this is ordinary isotonic regression
     * with adjusted targets (2.5, -1, 2.5, -1), whose pooled solution is
     * y_1 = ... = y_4 = 0.75.  Thus x = (0, 1.25, 1.75, 2.25, 2.75).
     * The boundary start must release the first constraint.  The interior
     * start must discover the two strongly active pair constraints.  The
     * middle tight constraint has zero multiplier, so retaining it in the
     * working set is deliberately not required.
     */
    MonotoneQpCoeff chain_coeff[4];
    MonotoneQpRow chain_rows[4];
    MonotoneQpBound chain_bounds[4];
    const double chain_targets[4] = {3.0, 0.0, 4.0, 1.0};
    const double chain_oracle[5] = {0.0, 1.25, 1.75, 2.25, 2.75};
    for (int i = 0; i < 4; i++) {
        chain_coeff[i].var = i + 1;
        chain_coeff[i].value = 1.0;
        chain_rows[i].first = (size_t)i;
        chain_rows[i].count = 1;
        chain_rows[i].target = chain_targets[i];
        chain_rows[i].weight = 1.0;
        chain_rows[i].kind = 101;
        chain_rows[i].owner = i;
        chain_bounds[i].lo = i;
        chain_bounds[i].hi = i + 1;
        chain_bounds[i].lower = 0.5;
        chain_bounds[i].stroke = 0;
        chain_bounds[i].edge = i;
    }
    MonotoneQpProblem chain_qp;
    memset(&chain_qp, 0, sizeof chain_qp);
    chain_qp.nvar = 5;
    chain_qp.rows = chain_rows;
    chain_qp.nrows = 4;
    chain_qp.coeff = chain_coeff;
    chain_qp.ncoeff = 4;
    chain_qp.bounds = chain_bounds;
    chain_qp.nbounds = 4;
    chain_qp.anchors = &anchor;
    chain_qp.nanchors = 1;

    MonotoneQpTraceEntry entries_c[64], entries_d[64];
    MonotoneQpTrace trace_c = {entries_c, 64, 0};
    MonotoneQpTrace trace_d = {entries_d, 64, 0};
    MonotoneQpStats stats_c, stats_d;
    double xc[5] = {0.0, 0.5, 1.0, 1.5, 2.0};
    double xd[5] = {0.0, 2.0, 3.0, 4.0, 5.0};
    int rc_c = MonotoneQp_solve(arena, &chain_qp, &opts,
                                xc, &trace_c, &stats_c);
    int rc_d = MonotoneQp_solve(arena, &chain_qp, &opts,
                                xd, &trace_d, &stats_d);
    double chain_boundary_error = 0.0;
    double chain_interior_error = 0.0;
    for (int i = 0; i < 5; i++) {
        double ec = fabs(xc[i] - chain_oracle[i]);
        double ed = fabs(xd[i] - chain_oracle[i]);
        if (ec > chain_boundary_error) chain_boundary_error = ec;
        if (ed > chain_interior_error) chain_interior_error = ed;
    }
    if (rc_c != 0 || chain_boundary_error > 1e-9 ||
        stats_c.dropped_total < 1 || stats_c.active_final != 3)
        fails++;
    if (rc_d != 0 || chain_interior_error > 1e-9 ||
        stats_d.added_total < 2 || stats_d.active_final < 2)
        fails++;

    write_trace_csv(dir, "01_bound_activation_trace.csv", &trace_a);
    write_trace_csv(dir, "02_bound_release_trace.csv", &trace_b);
    write_trace_csv(dir, "03_chain_from_boundary_trace.csv", &trace_c);
    write_trace_csv(dir, "04_chain_from_interior_trace.csv", &trace_d);
    char path[SYNTH_PATH_CAP];
    FILE *solution_fp = open_output(dir, "chain_solution.csv", path);
    if (solution_fp != NULL) {
        fprintf(solution_fp, "variable,oracle,boundary_start,interior_start\n");
        for (int i = 0; i < 5; i++)
            fprintf(solution_fp, "%d,%.17g,%.17g,%.17g\n",
                    i, chain_oracle[i], xc[i], xd[i]);
        fclose(solution_fp);
    }
    FILE *fp = open_output(dir, "metrics.json", path);
    if (fp != NULL) {
        fprintf(fp,
                "{\n  \"passed\": %s,\n"
                "  \"activation\": {\"rc\": %d, \"x\": [%.17g, %.17g], "
                "\"active_final\": %d},\n"
                "  \"release\": {\"rc\": %d, \"x\": [%.17g, %.17g], "
                "\"dropped\": %d, \"active_final\": %d},\n"
                "  \"chain_from_boundary\": {\"rc\": %d, "
                "\"max_oracle_error\": %.17g, \"added\": %d, "
                "\"dropped\": %d, \"active_final\": %d},\n"
                "  \"chain_from_interior\": {\"rc\": %d, "
                "\"max_oracle_error\": %.17g, \"added\": %d, "
                "\"dropped\": %d, \"active_final\": %d}\n}\n",
                fails ? "false" : "true", rc_a, xa[0], xa[1],
                stats_a.active_final, rc_b, xb[0], xb[1],
                stats_b.dropped_total, stats_b.active_final,
                rc_c, chain_boundary_error, stats_c.added_total,
                stats_c.dropped_total, stats_c.active_final,
                rc_d, chain_interior_error, stats_d.added_total,
                stats_d.dropped_total, stats_d.active_final);
        fclose(fp);
    }
    fprintf(stderr, "[atlas_strip_synth] active-set reference: %s\n",
            fails ? "FAIL" : "ok");
    return fails;
}

static int write_case_readme(const char *dir)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "README.txt", path);
    if (fp == NULL) return -1;
    fprintf(fp,
        "Two-level latent-strip synthetic\n"
        "================================\n\n"
        "00_observed_strokes_world.obj     Ordered 3D observations.\n"
        "01_cross_sections_world.obj       Known latent material rulings.\n"
        "02_raw_parameter.obj              Independently translated charts.\n"
        "03_relaxed_parameter.obj          StrokeStrip relaxed solve.\n"
        "04_relaxed_constraints.obj        Relaxed residual lines.\n"
        "05_final_parameter.obj            Final geometric alignment solve.\n"
        "06_final_constraints.obj          Final residual lines.\n"
        "07_mesh_world.obj                 Disjoint observed triangle patches.\n"
        "08_mesh_raw_parameter.obj         Mesh-only/raw chart gauges.\n"
        "08a_mesh_raw_resolved.obj         Raw resolved patches only.\n"
        "09_mesh_hybrid_parameter.obj      FEM extension from coarse atlas.\n"
        "10_mesh_hybrid_resolved.obj       Resolved patches only, for viewing.\n\n"
        "Pink is a normally displaced peeled duplicate. Red is a nearby but\n"
        "unrelated layer, deliberately excluded from the latent cross-sections.\n"
        "The primary slice at v=3 has a six-sample hole and is split into two\n"
        "independent strokes; other cross-sections must carry its global gauge.\n");
    fclose(fp);
    return 0;
}

static int write_case_metrics(const char *dir, int passed,
                              const StripFixture *fixture,
                              const AtlasStripMetrics *relaxed_metrics,
                              const AtlasStripMetrics *final_metrics,
                              const MonotoneQpStats *relaxed_stats,
                              const MonotoneQpStats *final_stats,
                              const MonotoneQpStats *mesh_stats,
                              const SynthMetrics *synth)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "metrics.json", path);
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"case\": \"translated_peel_hole\",\n"
        "  \"passed\": %s,\n"
        "  \"domain\": {\"samples\": %zu, \"strokes\": %zu, "
        "\"cross_sections\": %zu, \"members\": %zu, "
        "\"quotient_components\": 2, \"resolved_components\": 1},\n"
        "  \"relaxed\": {\n"
        "    \"iterations\": %d, \"spd_solves\": %d, "
        "\"objective_initial\": %.17g, \"objective_final\": %.17g,\n"
        "    \"energy_length\": %.17g, \"energy_align\": %.17g, "
        "\"energy_local\": %.17g,\n"
        "    \"rms_length_row\": %.17g, \"rms_align_row\": %.17g, "
        "\"min_monotone_ratio\": %.17g,\n"
        "    \"stationarity_residual\": %.17g, "
        "\"max_linear_residual\": %.17g\n"
        "  },\n"
        "  \"final\": {\n"
        "    \"iterations\": %d, \"spd_solves\": %d, "
        "\"objective_initial\": %.17g, \"objective_final\": %.17g,\n"
        "    \"energy_length\": %.17g, \"energy_align\": %.17g, "
        "\"energy_local\": %.17g,\n"
        "    \"rms_length_row\": %.17g, \"rms_align_row\": %.17g, "
        "\"max_align_residual\": %.17g, "
        "\"min_monotone_ratio\": %.17g,\n"
        "    \"stationarity_residual\": %.17g, "
        "\"max_linear_residual\": %.17g\n"
        "  },\n"
        "  \"ground_truth\": {\"max_resolved_error\": %.17g, "
        "\"max_peel_error\": %.17g, \"max_unresolved_internal_error\": %.17g},\n"
        "  \"mesh_extension\": {\n"
        "    \"iterations\": %d, \"objective_final\": %.17g, "
        "\"stationarity_residual\": %.17g,\n"
        "    \"observation_rms\": %.17g, \"local_delta_rms\": %.17g,\n"
        "    \"overlap_raw_rms\": %.17g, \"overlap_final_rms\": %.17g,\n"
        "    \"min_jacobian\": %.17g, \"flipped_triangles\": %zu\n"
        "  }\n"
        "}\n",
        passed ? "true" : "false",
        fixture->nsamples, fixture->nstrokes, fixture->ncross,
        fixture->nmembers,
        relaxed_stats->iterations, relaxed_stats->spd_solves,
        relaxed_stats->objective_initial, relaxed_stats->objective_final,
        relaxed_metrics->energy_length, relaxed_metrics->energy_align,
        relaxed_metrics->energy_local, relaxed_metrics->rms_length_row,
        relaxed_metrics->rms_align_row, relaxed_metrics->min_monotone_ratio,
        relaxed_stats->stationarity_residual,
        relaxed_stats->max_reduced_linear_residual,
        final_stats->iterations, final_stats->spd_solves,
        final_stats->objective_initial, final_stats->objective_final,
        final_metrics->energy_length, final_metrics->energy_align,
        final_metrics->energy_local, final_metrics->rms_length_row,
        final_metrics->rms_align_row, final_metrics->max_align_residual,
        final_metrics->min_monotone_ratio,
        final_stats->stationarity_residual,
        final_stats->max_reduced_linear_residual,
        synth->max_resolved_error, synth->max_peel_error,
        synth->max_unresolved_error,
        mesh_stats->iterations, mesh_stats->objective_final,
        mesh_stats->stationarity_residual,
        synth->mesh_observation_rms, synth->mesh_local_delta_rms,
        synth->overlap_raw_rms, synth->overlap_final_rms,
        synth->min_jacobian, synth->flipped_triangles);
    fclose(fp);
    return 0;
}

static int run_latent_strip_case(Arena_T arena, const char *root)
{
    char dir[SYNTH_PATH_CAP];
    if (snprintf(dir, sizeof dir, "%s/01_translated_peel_hole", root) < 0 ||
        make_dir(dir) != 0) return 1;

    StripFixture fixture;
    if (build_strip_fixture(arena, &fixture) != 0) return 1;
    size_t nvar = fixture.nsamples + fixture.ncross;
    double *raw = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    double *relaxed = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    double *final = (double *)ARENA_ALLOC(arena, nvar * sizeof(double));
    for (size_t i = 0; i < fixture.nsamples; i++) raw[i] = fixture.meta[i].u_raw;

    AtlasStripOptions relaxed_opts;
    AtlasStripOptions_default(&relaxed_opts);
    relaxed_opts.mode = ATLAS_STRIP_RELAXED;
    relaxed_opts.lambda_length = 1.0;
    relaxed_opts.lambda_align = 4.0;
    relaxed_opts.lambda_local = 0.05;
    relaxed_opts.monotone_fraction = 0.5;
    AtlasStrip_initialize_intercepts(&fixture.problem, &relaxed_opts, raw);
    memcpy(relaxed, raw, nvar * sizeof(double));

    AtlasStripSystem relaxed_system;
    if (AtlasStrip_build(arena, &fixture.problem, &relaxed_opts,
                         &relaxed_system) != 0) return 1;
    MonotoneQpOptions qp_opts;
    MonotoneQpOptions_default(&qp_opts);
    qp_opts.max_active_iterations = 512;
    MonotoneQpTraceEntry *relaxed_entries =
        (MonotoneQpTraceEntry *)ARENA_ALLOC(
            arena, 512 * sizeof(MonotoneQpTraceEntry));
    MonotoneQpTrace relaxed_trace = {relaxed_entries, 512, 0};
    MonotoneQpStats relaxed_stats;
    int relaxed_rc = MonotoneQp_solve(arena, &relaxed_system.qp, &qp_opts,
                                      relaxed, &relaxed_trace, &relaxed_stats);
    AtlasStripMetrics relaxed_metrics;
    AtlasStrip_measure(&fixture.problem, &relaxed_opts, relaxed,
                       &relaxed_metrics);

    AtlasStripOptions final_opts = relaxed_opts;
    final_opts.mode = ATLAS_STRIP_FINAL;
    memcpy(final, relaxed, nvar * sizeof(double));
    AtlasStrip_initialize_intercepts(&fixture.problem, &final_opts, final);
    AtlasStripSystem final_system;
    if (AtlasStrip_build(arena, &fixture.problem, &final_opts,
                         &final_system) != 0) return 1;
    MonotoneQpTraceEntry *final_entries =
        (MonotoneQpTraceEntry *)ARENA_ALLOC(
            arena, 512 * sizeof(MonotoneQpTraceEntry));
    MonotoneQpTrace final_trace = {final_entries, 512, 0};
    MonotoneQpStats final_stats;
    int final_rc = MonotoneQp_solve(arena, &final_system.qp, &qp_opts,
                                    final, &final_trace, &final_stats);
    AtlasStripMetrics final_metrics;
    AtlasStrip_measure(&fixture.problem, &final_opts, final, &final_metrics);

    MeshFixture mesh;
    int mesh_build_rc = build_mesh_fixture(arena, &fixture, final, &mesh);
    MonotoneQpProblem mesh_qp;
    int mesh_qp_rc = mesh_build_rc == 0 ? build_mesh_qp(arena, &mesh, &mesh_qp) : -1;
    double *mesh_u = mesh_build_rc == 0
                   ? (double *)ARENA_ALLOC(arena, mesh.nv * sizeof(double)) : NULL;
    if (mesh_u != NULL) memcpy(mesh_u, mesh.u0, mesh.nv * sizeof(double));
    MonotoneQpTraceEntry *mesh_entries =
        (MonotoneQpTraceEntry *)ARENA_ALLOC(
            arena, 512 * sizeof(MonotoneQpTraceEntry));
    MonotoneQpTrace mesh_trace = {mesh_entries, 512, 0};
    MonotoneQpStats mesh_stats;
    memset(&mesh_stats, 0, sizeof mesh_stats);
    int mesh_solve_rc = mesh_qp_rc == 0
                      ? MonotoneQp_solve(arena, &mesh_qp, &qp_opts,
                                         mesh_u, &mesh_trace, &mesh_stats)
                      : -1;

    SynthMetrics synth;
    memset(&synth, 0, sizeof synth);
    measure_sample_errors(&fixture, final, &synth);
    if (mesh_solve_rc == 0) measure_mesh(&mesh, mesh_u, &synth);

    int passed = relaxed_rc == 0 && final_rc == 0 && mesh_solve_rc == 0 &&
                 synth.max_resolved_error < 1e-7 &&
                 synth.max_peel_error < 1e-7 &&
                 final_metrics.min_monotone_ratio >= 0.5 - 1e-10 &&
                 synth.overlap_raw_rms > 10.0 &&
                 synth.overlap_final_rms < 1e-6 &&
                 synth.flipped_triangles == 0 &&
                 synth.min_jacobian > 0.0;

    write_case_readme(dir);
    write_strokes_world(dir, "00_observed_strokes_world.obj", &fixture);
    write_cross_sections_world(dir, "01_cross_sections_world.obj", &fixture);
    write_parameter_strokes(dir, "02_raw_parameter.obj", &fixture, raw);
    write_parameter_strokes(dir, "03_relaxed_parameter.obj", &fixture, relaxed);
    write_cross_sections_parameter(dir, "04_relaxed_constraints.obj",
                                   &fixture, &relaxed_opts, relaxed);
    write_parameter_strokes(dir, "05_final_parameter.obj", &fixture, final);
    write_cross_sections_parameter(dir, "06_final_constraints.obj",
                                   &fixture, &final_opts, final);
    if (mesh_build_rc == 0) {
        write_mesh_obj(dir, "07_mesh_world.obj", &mesh, mesh.u0, 0);
        write_mesh_obj(dir, "08_mesh_raw_parameter.obj", &mesh, mesh.u0, 1);
        write_mesh_obj_patch_limit(dir, "08a_mesh_raw_resolved.obj",
                                   &mesh, mesh.u0, 1, 2);
        if (mesh_solve_rc == 0)
            write_mesh_obj(dir, "09_mesh_hybrid_parameter.obj", &mesh,
                           mesh_u, 1);
        if (mesh_solve_rc == 0)
            write_mesh_obj_patch_limit(dir, "10_mesh_hybrid_resolved.obj",
                                       &mesh, mesh_u, 1, 2);
        write_mesh_csv(dir, &mesh, mesh_solve_rc == 0 ? mesh_u : mesh.u0);
    }
    write_samples_csv(dir, &fixture, raw, relaxed, final);
    write_cross_sections_csv(dir, &fixture, raw, relaxed, final);
    write_rows_csv(dir, "relaxed_rows.csv", &relaxed_system.qp);
    write_rows_csv(dir, "final_rows.csv", &final_system.qp);
    write_bounds_csv(dir, "final_monotonicity.csv", &final_system.qp, final);
    write_trace_csv(dir, "relaxed_solver_trace.csv", &relaxed_trace);
    write_trace_csv(dir, "final_solver_trace.csv", &final_trace);
    write_matrix_market(arena, dir, "relaxed_H.mtx", &relaxed_system.qp);
    write_matrix_market(arena, dir, "final_H.mtx", &final_system.qp);
    write_rhs_csv(dir, "relaxed_rhs.csv", &relaxed_system.qp);
    write_rhs_csv(dir, "final_rhs.csv", &final_system.qp);
    if (mesh_qp_rc == 0) {
        write_rows_csv(dir, "mesh_rows.csv", &mesh_qp);
        write_trace_csv(dir, "mesh_solver_trace.csv", &mesh_trace);
        write_matrix_market(arena, dir, "mesh_H.mtx", &mesh_qp);
        write_rhs_csv(dir, "mesh_rhs.csv", &mesh_qp);
    }
    write_case_metrics(dir, passed, &fixture, &relaxed_metrics, &final_metrics,
                       &relaxed_stats, &final_stats, &mesh_stats, &synth);

    fprintf(stderr,
            "[atlas_strip_synth] latent strip: %s "
            "(sample max %.3e, peel %.3e, overlap %.3f -> %.3e, "
            "Jmin %.3e, flips %zu)\n",
            passed ? "ok" : "FAIL", synth.max_resolved_error,
            synth.max_peel_error, synth.overlap_raw_rms,
            synth.overlap_final_rms, synth.min_jacobian,
            synth.flipped_triangles);
    return passed ? 0 : 1;
}

static int write_lift_texture(const char *dir)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "checker.ppm", path);
    if (fp == NULL) return -1;
    fprintf(fp, "P6\n64 64\n255\n");
    for (int y = 0; y < 64; y++) {
        for (int x = 0; x < 64; x++) {
            int checker = ((x / 8) ^ (y / 8)) & 1;
            unsigned char rgb[3] = {
                (unsigned char)(checker ? 235 : 45),
                (unsigned char)(checker ? 185 : 75),
                (unsigned char)(checker ? 55 : 205)
            };
            if (x < 3) { rgb[0] = 255; rgb[1] = 40; rgb[2] = 40; }
            if (y < 3) { rgb[0] = 40; rgb[1] = 255; rgb[2] = 40; }
            if (fwrite(rgb, 1, 3, fp) != 3) {
                fclose(fp);
                return -1;
            }
        }
    }
    if (fclose(fp) != 0) return -1;
    fp = open_output(dir, "component_lift.mtl", path);
    if (fp == NULL) return -1;
    fprintf(fp, "newmtl checker\nKd 1 1 1\nKa 0 0 0\nmap_Kd checker.ppm\n");
    return fclose(fp) == 0 ? 0 : -1;
}

static int write_lift_flat_obj(const char *dir, const char *name,
                               const double shift[3])
{
    static const double u0[3][2] = {{0.0, 8.0}, {0.0, 8.0}, {8.0, 12.0}};
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "mtllib component_lift.mtl\nusemtl checker\n");
    for (int c = 0; c < 3; c++) {
        double lo = u0[c][0] + shift[c], hi = u0[c][1] + shift[c];
        fprintf(fp, "g component_%d\n", c);
        fprintf(fp, "v %.9g 0 0\n", lo);
        fprintf(fp, "v %.9g 0 0\n", hi);
        fprintf(fp, "v %.9g 4 0\n", hi);
        fprintf(fp, "v %.9g 4 0\n", lo);
    }
    for (int c = 0; c < 3; c++) {
        fprintf(fp, "vt 0 0\nvt 1 0\nvt 1 1\nvt 0 1\n");
    }
    for (int c = 0; c < 3; c++) {
        int b = 4 * c + 1;
        fprintf(fp, "g component_%d\n", c);
        fprintf(fp, "f %d/%d %d/%d %d/%d\n",
                b, b, b + 1, b + 1, b + 2, b + 2);
        fprintf(fp, "f %d/%d %d/%d %d/%d\n",
                b, b, b + 2, b + 2, b + 3, b + 3);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

static int run_component_lift_case(Arena_T arena, const char *root)
{
    char dir[SYNTH_PATH_CAP];
    if (snprintf(dir, sizeof dir, "%s/03_component_gauge_lift", root) < 0 ||
        make_dir(dir) != 0) return 1;

    double weight[3] = {100.0, 10.0, 5.0};
    AtlasComponentLiftEdge protected_edge = {1, 2, 0.0, 1.0, 0};
    AtlasComponentLiftEdge separation_edge = {0, 1, 20.0, 1.0, 1};
    AtlasComponentLiftProblem problem = {
        3, weight, NULL, &protected_edge, 1, &separation_edge, 1, 1e-8
    };
    double shift[3] = {0.0, 0.0, 0.0};
    int32_t group[3] = {-1, -1, -1};
    AtlasComponentLiftStats stats;
    int rc = AtlasComponentLift_solve(
        arena, &problem, shift, group, &stats);

    double max_edge_change = 0.0;
    static const double span[3] = {8.0, 8.0, 4.0};
    for (int c = 0; c < 3; c++) {
        double after = (span[c] + shift[c]) - shift[c];
        double change = fabs(after - span[c]);
        if (change > max_edge_change) max_edge_change = change;
    }
    int overlap_before = 1;
    int overlap_after = !(8.0 <= shift[1] || shift[1] + 8.0 <= 0.0);
    int passed = rc == 0 && fabs(shift[0]) < 1e-12 &&
                 fabs(shift[1] - 20.0) < 1e-12 &&
                 fabs(shift[2] - 20.0) < 1e-12 &&
                 max_edge_change < 1e-12 && overlap_before && !overlap_after &&
                 stats.moved_components == 2;

    double zero[3] = {0.0, 0.0, 0.0};
    write_lift_texture(dir);
    write_lift_flat_obj(dir, "00_head_collapsed_flat_textured.obj", zero);
    write_lift_flat_obj(dir, "01_lifted_flat_textured.obj", shift);

    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "components.csv", path);
    if (fp != NULL) {
        fprintf(fp, "component,weight,protected_group,shift,faces_before,faces_after\n");
        for (int c = 0; c < 3; c++)
            fprintf(fp, "%d,%.17g,%d,%.17g,2,2\n",
                    c, weight[c], group[c], shift[c]);
        fclose(fp);
    }
    fp = open_output(dir, "constraints.csv", path);
    if (fp != NULL) {
        fprintf(fp, "kind,component0,component1,target,residual_after\n");
        fprintf(fp, "protected,1,2,0,%.17g\n", shift[2] - shift[1]);
        fprintf(fp, "separate,0,1,20,%.17g\n", shift[1] - shift[0] - 20.0);
        fclose(fp);
    }
    fp = open_output(dir, "metrics.json", path);
    if (fp != NULL) {
        fprintf(fp,
            "{\n  \"case\": \"component_gauge_lift\",\n"
            "  \"passed\": %s,\n  \"solve_rc\": %d,\n"
            "  \"faces_before\": 6,\n  \"faces_after\": 6,\n"
            "  \"overlap_before\": %s,\n  \"overlap_after\": %s,\n"
            "  \"max_local_edge_change\": %.17g,\n"
            "  \"shift\": [%.17g, %.17g, %.17g],\n"
            "  \"protected_groups\": %zu,\n"
            "  \"moved_components\": %zu\n}\n",
            passed ? "true" : "false", rc,
            overlap_before ? "true" : "false",
            overlap_after ? "true" : "false", max_edge_change,
            shift[0], shift[1], shift[2], stats.protected_groups,
            stats.moved_components);
        fclose(fp);
    }
    fp = open_output(dir, "README.txt", path);
    if (fp != NULL) {
        fprintf(fp,
            "Component gauge-lift regression\n"
            "===============================\n\n"
            "00 is the deliberately collapsed HEAD atlas: components 0 and 1\n"
            "occupy the same UV rectangle. Component 2 is a protected, attached\n"
            "continuation of component 1. 01 applies one constant +20 gauge to\n"
            "the entire protected {1,2} chart. All six faces and all local edge\n"
            "differences are unchanged; only the false global overlap is gone.\n"
            "Both OBJs use checker.ppm so continuity can be inspected directly.\n");
        fclose(fp);
    }
    fprintf(stderr,
            "[atlas_strip_synth] component gauge lift: %s "
            "(faces 6 -> 6, overlap %d -> %d, local delta %.3e)\n",
            passed ? "ok" : "FAIL", overlap_before, overlap_after,
            max_edge_change);
    return passed ? 0 : 1;
}

enum {
    SMOOTH_NX = 5,
    SMOOTH_NY = 5,
    SMOOTH_NPATCH = 4,
    SMOOTH_NV = SMOOTH_NX * SMOOTH_NY * SMOOTH_NPATCH,
    SMOOTH_NF = (SMOOTH_NX - 1) * (SMOOTH_NY - 1) * 2 * SMOOTH_NPATCH
};

static int32_t smooth_vertex(int patch, int y, int x)
{
    return (int32_t)(patch * SMOOTH_NX * SMOOTH_NY +
                     y * SMOOTH_NX + x);
}

static int write_smooth_flat_obj(const char *dir, const char *name,
                                 const double *position,
                                 const double *parameter,
                                 const double *vcoord,
                                 const int32_t *faces)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, name, path);
    if (fp == NULL) return -1;
    fprintf(fp, "mtllib component_lift.mtl\nusemtl checker\n");
    for (int i = 0; i < SMOOTH_NV; i++) {
        int patch = i / (SMOOTH_NX * SMOOTH_NY);
        fprintf(fp, "v %.17g %.17g %.17g\n",
                parameter[i], vcoord[i], 0.03 * (double)patch);
    }
    for (int i = 0; i < SMOOTH_NV; i++)
        fprintf(fp, "vt %.17g %.17g\n",
                position[(size_t)i * 3] / 8.0,
                position[(size_t)i * 3 + 1] / 4.0);
    for (int f = 0; f < SMOOTH_NF; f++) {
        int a = faces[(size_t)f * 3] + 1;
        int b = faces[(size_t)f * 3 + 1] + 1;
        int c = faces[(size_t)f * 3 + 2] + 1;
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    return fclose(fp) == 0 ? 0 : -1;
}

/* A deliberately over-constrained chart fixture which a constant component
 * lift cannot solve.  Each physical ribbon is split at x=4.  The right chart
 * has both a bad gauge and a y-dependent drift, so all five exact seam copies
 * cannot be closed by one rigid shift.  A second ribbon starts on top of the
 * first in UV; sparse latent support offsets ask it to occupy the next wrap.
 * The FEM correction must close both seams, vary smoothly within the right
 * charts, and separate the two complete ribbons without dropping a face. */
static int run_smooth_quotient_case(Arena_T arena, const char *root)
{
    char dir[SYNTH_PATH_CAP];
    if (snprintf(dir, sizeof dir, "%s/04_smooth_seam_quotient", root) < 0 ||
        make_dir(dir) != 0) return 1;

    double *position = (double *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NV * 3 * sizeof(double));
    double *vcoord = (double *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NV * sizeof(double));
    double *u0 = (double *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NV * sizeof(double));
    int32_t *vertex_component = (int32_t *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NV * sizeof(int32_t));
    int32_t *faces = (int32_t *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NF * 3 * sizeof(int32_t));
    size_t nf = 0;
    for (int patch = 0; patch < SMOOTH_NPATCH; patch++) {
        int right = patch & 1;
        int layer = patch >> 1;
        for (int y = 0; y < SMOOTH_NY; y++) {
            for (int x = 0; x < SMOOTH_NX; x++) {
                int32_t v = smooth_vertex(patch, y, x);
                double px = (double)(right ? 4 + x : x);
                double py = (double)y;
                position[(size_t)v * 3] = px;
                position[(size_t)v * 3 + 1] = py;
                position[(size_t)v * 3 + 2] = 2.0 * (double)layer;
                vcoord[v] = py;
                vertex_component[v] = patch;
                if (patch == 0 || patch == 2)
                    u0[v] = px;
                else
                    u0[v] = px + 12.0 + 0.5 * py;
            }
        }
        for (int y = 0; y + 1 < SMOOTH_NY; y++) {
            for (int x = 0; x + 1 < SMOOTH_NX; x++) {
                int32_t a = smooth_vertex(patch, y, x);
                int32_t b = smooth_vertex(patch, y, x + 1);
                int32_t c = smooth_vertex(patch, y + 1, x);
                int32_t d = smooth_vertex(patch, y + 1, x + 1);
                faces[nf * 3] = a;
                faces[nf * 3 + 1] = b;
                faces[nf * 3 + 2] = d;
                nf++;
                faces[nf * 3] = a;
                faces[nf * 3 + 1] = d;
                faces[nf * 3 + 2] = c;
                nf++;
            }
        }
    }
    if (nf != SMOOTH_NF) return 1;

    const size_t seam_rows = 2 * SMOOTH_NY;
    const size_t support_rows = 4 * SMOOTH_NY;
    const size_t nconstraint = seam_rows + support_rows + 1;
    const size_t ncoeff = 2 * nconstraint;
    AtlasFieldConstraint *constraint = (AtlasFieldConstraint *)ARENA_ALLOC(
        arena, nconstraint * sizeof(AtlasFieldConstraint));
    AtlasFieldConstraintCoeff *coeff =
        (AtlasFieldConstraintCoeff *)ARENA_ALLOC(
            arena, ncoeff * sizeof(AtlasFieldConstraintCoeff));
    size_t nr = 0, nk = 0;
#define SMOOTH_PAIR_ROW(kind_value, source_value, target_value, weight_value, \
                        var0, coeff0, var1, coeff1) do {                    \
        constraint[nr].first = nk;                                         \
        constraint[nr].count = 2;                                          \
        constraint[nr].target = (target_value);                             \
        constraint[nr].weight = (weight_value);                             \
        constraint[nr].kind = (kind_value);                                 \
        constraint[nr].source = (source_value);                             \
        coeff[nk].variable = (var0); coeff[nk].coefficient = (coeff0); nk++; \
        coeff[nk].variable = (var1); coeff[nk].coefficient = (coeff1); nk++; \
        nr++;                                                               \
    } while (0)
    for (int y = 0; y < SMOOTH_NY; y++) {
        int32_t a = smooth_vertex(0, y, SMOOTH_NX - 1);
        int32_t b = smooth_vertex(1, y, 0);
        SMOOTH_PAIR_ROW(ATLAS_FIELD_ROW_SEAM, y, 0.0, 1000.0,
                        a, -1.0, b, 1.0);
        a = smooth_vertex(2, y, SMOOTH_NX - 1);
        b = smooth_vertex(3, y, 0);
        SMOOTH_PAIR_ROW(ATLAS_FIELD_ROW_SEAM, SMOOTH_NY + y,
                        0.0, 1000.0, a, -1.0, b, 1.0);
    }
    int32_t q_main = SMOOTH_NV;
    int32_t q_alias = SMOOTH_NV + 1;
    for (int y = 0; y < SMOOTH_NY; y++) {
        for (int x = SMOOTH_NX - 2; x < SMOOTH_NX; x++) {
            int32_t main_v = smooth_vertex(1, y, x);
            int32_t alias_v = smooth_vertex(3, y, x);
            SMOOTH_PAIR_ROW(ATLAS_FIELD_ROW_OBSERVATION, main_v,
                            u0[main_v], 10.0,
                            main_v, 1.0, q_main, -1.0);
            SMOOTH_PAIR_ROW(ATLAS_FIELD_ROW_OBSERVATION, alias_v,
                            u0[alias_v], 10.0,
                            alias_v, 1.0, q_alias, -1.0);
        }
    }
    SMOOTH_PAIR_ROW(ATLAS_FIELD_ROW_ORDER, 0, 16.0, 100.0,
                    q_main, -1.0, q_alias, 1.0);
#undef SMOOTH_PAIR_ROW
    if (nr != nconstraint || nk != ncoeff) return 1;

    MonotoneQpAnchor anchor = {0, u0[0], 0};
    AtlasFieldProblem problem;
    memset(&problem, 0, sizeof problem);
    problem.position = position;
    problem.u0 = u0;
    problem.nvertex = SMOOTH_NV;
    problem.nauxiliary = 2;
    problem.triangle = faces;
    problem.ntriangle = SMOOTH_NF;
    problem.constraint = constraint;
    problem.nconstraint = nconstraint;
    problem.constraint_coeff = coeff;
    problem.nconstraint_coeff = ncoeff;
    problem.anchor = &anchor;
    problem.nanchor = 1;
    AtlasFieldSystem system;
    int build_rc = AtlasField_build(arena, &problem, &system);
    double *solution = (double *)ARENA_ALLOC(
        arena, ((size_t)SMOOTH_NV + 2) * sizeof(double));
    memcpy(solution, u0, (size_t)SMOOTH_NV * sizeof(double));
    solution[q_main] = 0.0;
    solution[q_alias] = 0.0;
    MonotoneQpOptions options;
    MonotoneQpOptions_default(&options);
    MonotoneQpTraceEntry trace_entry[32];
    MonotoneQpTrace trace = {trace_entry, 32, 0};
    MonotoneQpStats stats;
    memset(&stats, 0, sizeof stats);
    int solve_rc = build_rc == 0 ? MonotoneQp_solve(
        arena, &system.qp, &options, solution, &trace, &stats) : -1;
    double *helper_u = (double *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NV * sizeof(double));
    double helper_aux[2] = {0.0, 0.0};
    double auxiliary_initial[2] = {0.0, 0.0};
    AtlasFieldRefineStats helper_stats;
    int helper_rc;
    /* The helper consumes float production geometry. Keep a literal float
     * copy here so this regression exercises the exact adapter path. */
    float *float_position = (float *)ARENA_ALLOC(
        arena, (size_t)SMOOTH_NV * 3 * sizeof(float));
    for (int i = 0; i < SMOOTH_NV * 3; i++)
        float_position[i] = (float)position[i];
    helper_rc = AtlasFieldRefine_solve(
        arena, float_position, SMOOTH_NV, faces, SMOOTH_NF,
        vertex_component, SMOOTH_NPATCH, u0,
        constraint, nconstraint, coeff, ncoeff,
        2, auxiliary_initial, &options, helper_u, helper_aux,
        &helper_stats);
    double helper_max_difference = 0.0;
    if (helper_rc == 0)
        for (int i = 0; i < SMOOTH_NV; i++) {
            double d = fabs(helper_u[i] - solution[i]);
            if (d > helper_max_difference) helper_max_difference = d;
        }

    double seam_raw2 = 0.0, seam_final2 = 0.0;
    size_t nseam = 0;
    for (int layer = 0; layer < 2; layer++) {
        for (int y = 0; y < SMOOTH_NY; y++) {
            int32_t a = smooth_vertex(2 * layer, y, SMOOTH_NX - 1);
            int32_t b = smooth_vertex(2 * layer + 1, y, 0);
            double r0 = u0[b] - u0[a];
            double r1 = solution[b] - solution[a];
            seam_raw2 += r0 * r0;
            seam_final2 += r1 * r1;
            nseam++;
        }
    }
    double seam_raw_rms = sqrt(seam_raw2 / (double)nseam);
    double seam_final_rms = sqrt(seam_final2 / (double)nseam);
    double order_value = solution[q_alias] - solution[q_main];
    double support_raw_difference = 0.0;
    double support_final_difference = 0.0;
    size_t nsupport = 0;
    for (int y = 0; y < SMOOTH_NY; y++) {
        for (int x = SMOOTH_NX - 2; x < SMOOTH_NX; x++) {
            int32_t main_v = smooth_vertex(1, y, x);
            int32_t alias_v = smooth_vertex(3, y, x);
            support_raw_difference += u0[alias_v] - u0[main_v];
            support_final_difference +=
                solution[alias_v] - solution[main_v];
            nsupport++;
        }
    }
    support_raw_difference /= (double)nsupport;
    support_final_difference /= (double)nsupport;
    double cmin = DBL_MAX, cmax = -DBL_MAX;
    for (int y = 0; y < SMOOTH_NY; y++) {
        for (int x = 0; x < SMOOTH_NX; x++) {
            int32_t v = smooth_vertex(1, y, x);
            double d = solution[v] - u0[v];
            if (d < cmin) cmin = d;
            if (d > cmax) cmax = d;
        }
    }
    double correction_span = cmax - cmin;
    double min_jacobian = DBL_MAX;
    size_t flips = 0;
    for (int f = 0; f < SMOOTH_NF; f++) {
        int32_t a = faces[(size_t)f * 3];
        int32_t b = faces[(size_t)f * 3 + 1];
        int32_t c = faces[(size_t)f * 3 + 2];
        double jac = (solution[b] - solution[a]) *
                     (vcoord[c] - vcoord[a]) -
                     (solution[c] - solution[a]) *
                     (vcoord[b] - vcoord[a]);
        if (jac < min_jacobian) min_jacobian = jac;
        if (jac <= 0.0) flips++;
    }
    int passed = build_rc == 0 && solve_rc == 0 && helper_rc == 0 &&
                 helper_max_difference < 1e-9 &&
                 seam_final_rms < 0.05 &&
                 seam_final_rms < 0.01 * seam_raw_rms &&
                 fabs(order_value - 16.0) < 0.1 &&
                 fabs(support_raw_difference) < 1e-12 &&
                 fabs(support_final_difference - 16.0) < 0.1 &&
                 correction_span > 0.25 && flips == 0 &&
                 min_jacobian > 0.0;

    write_lift_texture(dir);
    write_smooth_flat_obj(dir, "00_head_collapsed_flat_textured.obj",
                          position, u0, vcoord, faces);
    if (solve_rc == 0)
        write_smooth_flat_obj(dir, "01_smooth_quotient_flat_textured.obj",
                              position, solution, vcoord, faces);
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(dir, "vertices.csv", path);
    if (fp != NULL) {
        fprintf(fp, "vertex,patch,x,y,z,u_head,u_final,correction\n");
        for (int i = 0; i < SMOOTH_NV; i++)
            fprintf(fp, "%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    i, i / (SMOOTH_NX * SMOOTH_NY),
                    position[(size_t)i * 3], position[(size_t)i * 3 + 1],
                    position[(size_t)i * 3 + 2], u0[i], solution[i],
                    solution[i] - u0[i]);
        fclose(fp);
    }
    fp = open_output(dir, "metrics.json", path);
    if (fp != NULL) {
        fprintf(fp,
            "{\n  \"case\": \"smooth_seam_quotient\",\n"
            "  \"passed\": %s,\n  \"build_rc\": %d,\n"
            "  \"solve_rc\": %d,\n"
            "  \"helper_rc\": %d,\n"
            "  \"helper_max_difference\": %.17g,\n"
            "  \"faces_before\": %d,\n"
            "  \"faces_after\": %d,\n"
            "  \"seam_rms_before\": %.17g,\n"
            "  \"seam_rms_after\": %.17g,\n"
            "  \"order_target\": 16,\n"
            "  \"order_value\": %.17g,\n"
            "  \"support_separation_before\": %.17g,\n"
            "  \"support_separation_after\": %.17g,\n"
            "  \"right_chart_correction_span\": %.17g,\n"
            "  \"minimum_jacobian\": %.17g,\n"
            "  \"flipped_triangles\": %zu,\n"
            "  \"objective_initial\": %.17g,\n"
            "  \"objective_final\": %.17g,\n"
            "  \"linear_residual\": %.17g\n}\n",
            passed ? "true" : "false", build_rc, solve_rc,
            helper_rc, helper_max_difference, SMOOTH_NF, SMOOTH_NF,
            seam_raw_rms, seam_final_rms,
            order_value, support_raw_difference, support_final_difference,
            correction_span, min_jacobian, flips,
            stats.objective_initial, stats.objective_final,
            stats.max_reduced_linear_residual);
        fclose(fp);
    }
    fp = open_output(dir, "README.txt", path);
    if (fp != NULL) {
        fprintf(fp,
            "Smooth seam-quotient regression\n"
            "===============================\n\n"
            "00 shows two complete ribbons collapsed in UV, with a varying\n"
            "gauge error across each x=4 cube seam. 01 is the exact same 128\n"
            "faces after the sparse FEM solve. Five pointwise seam springs per\n"
            "ribbon close the physical boundary; two latent support offsets\n"
            "move the second coherent ribbon to the next wrap. The correction\n"
            "within the right chart is intentionally non-constant, so a rigid\n"
            "component lift cannot pass this test. Both OBJs use checker.ppm.\n");
        fclose(fp);
    }
    if (build_rc == 0) {
        write_rows_csv(dir, "field_rows.csv", &system.qp);
        write_matrix_market(arena, dir, "field_H.mtx", &system.qp);
        write_rhs_csv(dir, "field_rhs.csv", &system.qp);
        write_trace_csv(dir, "field_solver_trace.csv", &trace);
    }
    fprintf(stderr,
            "[atlas_strip_synth] smooth seam quotient: %s "
            "(seam %.3f -> %.3e, support %.3f -> %.3f, "
            "correction span %.3f, "
            "faces %d -> %d, flips %zu)\n",
            passed ? "ok" : "FAIL", seam_raw_rms, seam_final_rms,
            support_raw_difference, support_final_difference,
            correction_span, SMOOTH_NF, SMOOTH_NF, flips);
    return passed ? 0 : 1;
}

static int write_manifest(const char *root, int active_fails,
                           int strip_fails, int density_fails,
                           int candidate_fails, int lift_fails,
                           int smooth_fails)
{
    char path[SYNTH_PATH_CAP];
    FILE *fp = open_output(root, "manifest.json", path);
    if (fp == NULL) return -1;
    int passed = active_fails == 0 && strip_fails == 0 &&
                 density_fails == 0 && candidate_fails == 0 &&
                 lift_fails == 0 && smooth_fails == 0;
    fprintf(fp,
            "{\n  \"suite\": \"atlas_strip_synth\",\n"
            "  \"passed\": %s,\n"
            "  \"cases\": [\n"
            "    {\"name\": \"active_set_reference\", \"passed\": %s, "
            "\"directory\": \"00_active_set_reference\"},\n"
            "    {\"name\": \"translated_peel_hole\", \"passed\": %s, "
            "\"directory\": \"01_translated_peel_hole\"},\n"
            "    {\"name\": \"barycentric_robust\", \"passed\": %s, "
            "\"directory\": \"02_barycentric_robust\"},\n"
            "    {\"name\": \"geometry_candidate_adapter\", \"passed\": %s},\n"
            "    {\"name\": \"component_gauge_lift\", \"passed\": %s, "
            "\"directory\": \"03_component_gauge_lift\"},\n"
            "    {\"name\": \"smooth_seam_quotient\", \"passed\": %s, "
            "\"directory\": \"04_smooth_seam_quotient\"}\n"
            "  ]\n}\n",
            passed ? "true" : "false",
            active_fails ? "false" : "true",
            strip_fails ? "false" : "true",
            density_fails ? "false" : "true",
            candidate_fails ? "false" : "true",
            lift_fails ? "false" : "true",
            smooth_fails ? "false" : "true");
    fclose(fp);
    return 0;
}

static void usage(const char *exe)
{
    fprintf(stderr, "Usage: %s [--out DIR]\n", exe);
}

int main(int argc, char **argv)
{
    const char *root = "output/atlas_strip_synth";
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            root = argv[++i];
        } else if (strcmp(argv[i], "--help") == 0 ||
                   strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else {
            usage(argv[0]);
            return 2;
        }
    }
    if (make_dir(root) != 0) {
        fprintf(stderr, "atlas_strip_synth: cannot create %s\n", root);
        return 1;
    }

    Arena_T arena = Arena_new();
    int active_fails = run_active_set_reference(arena, root);
    int strip_fails = run_latent_strip_case(arena, root);
    int density_fails = run_density_robust_case(arena, root);
    int candidate_fails = AtlasCandidates_selftest();
    int lift_fails = AtlasComponentLift_selftest();
    lift_fails += AtlasOverlapAudit_selftest();
    lift_fails += AtlasSeamAudit_selftest();
    lift_fails += run_component_lift_case(arena, root);
    int smooth_fails = run_smooth_quotient_case(arena, root);
    write_manifest(root, active_fails, strip_fails, density_fails,
                   candidate_fails, lift_fails, smooth_fails);
    Arena_dispose(&arena);

    int fails = active_fails + strip_fails + density_fails + candidate_fails +
                lift_fails + smooth_fails;
    fprintf(stderr, "[atlas_strip_synth] %s (%d failure%s), output=%s\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s",
            root);
    return fails ? 1 : 0;
}
