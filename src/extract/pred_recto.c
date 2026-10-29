#include "pred_recto.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../common/pipeline_constants.h"

/* Binary occupancy interpolated at voxel centres; an unavailable sample
 * invalidates a ray instead of inventing an edge at the volume boundary. */
static int sample(const PredRectoField *f, const float p[3], double *out)
{
    int k[3] = {0, 0, 0};
    int dims[3] = {f->d, f->h, f->w};
    double frac[3] = {0.0, 0.0, 0.0};
    double result = 0.0;
    for (int a = 0; a < 3; a++) {
        double q = (double)p[a] + (double)f->halo - 0.5;
        if (!isfinite(q) || q < 0.0 || q >= (double)(dims[a] - 1)) return 0;
        k[a] = (int)floor(q);
        frac[a] = q - (double)k[a];
    }
    for (int dz = 0; dz < 2; dz++)
        for (int dy = 0; dy < 2; dy++)
            for (int dx = 0; dx < 2; dx++) {
                size_t i = ((size_t)(k[0] + dz) * (size_t)f->h +
                            (size_t)(k[1] + dy)) * (size_t)f->w + (size_t)(k[2] + dx);
                double wt = (dz ? frac[0] : 1.0 - frac[0]) *
                            (dy ? frac[1] : 1.0 - frac[1]) *
                            (dx ? frac[2] : 1.0 - frac[2]);
                if (f->vol[i]) result += wt;
            }
    *out = result;
    return 1;
}

int PredRecto_target(const PredRectoField *f, const float pos[3],
                     const float normal[3], float inward[3], float *distance,
                     PredRectoStats *s)
{
    double values[PR_RECTO_SAMPLES] = {0.0};
    double ay = 0.0, ax = 0.0, len = 0.0, radial = 0.0, dot = 0.0;
    double best = 1e30, edge = 0.0, outer = 0.0;
    int seed = -1, lo = 0, hi = 0;
    s->tested++;
    *distance = 0.0f;
    for (int a = 0; a < 3; a++) len += (double)normal[a] * normal[a];
    BpaReconGate_axis_at(&f->axis, (double)pos[0] + f->axis.origin_z, &ay, &ax);
    ay -= (double)pos[1] + f->axis.origin_y;
    ax -= (double)pos[2] + f->axis.origin_x;
    radial = hypot(ay, ax);
    if (!isfinite(len) || len < 1e-12 || !isfinite(radial) || radial < 1e-6) {
        s->ambiguous++; return 0;
    }
    len = sqrt(len);
    dot = ((double)normal[1] * ay + (double)normal[2] * ax) / (len * radial);
    if (!isfinite(dot) || fabs(dot) < PR_RECTO_MIN_RADIAL_COS) {
        s->ambiguous++; return 0;
    }
    for (int a = 0; a < 3; a++) inward[a] = (float)((dot < 0.0 ? -1.0 : 1.0) * normal[a] / len);
    for (int i = 0; i < PR_RECTO_SAMPLES; i++) {
        double t = ((double)i - PR_RECTO_HALF_SAMPLES) * PR_RECTO_STEP;
        float q[3] = {0.0f, 0.0f, 0.0f};
        for (int a = 0; a < 3; a++) q[a] = (float)((double)pos[a] + t * inward[a]);
        if (!sample(f, q, &values[i])) { s->outside++; return 0; }
        if (values[i] >= 0.5 && fabs(t) < best) { best = fabs(t); seed = i; }
    }
    if (seed < 0 || best > PR_RECTO_SEED_REACH) { s->missing++; return 0; }
    lo = hi = seed;
    while (lo > 0 && values[lo - 1] >= 0.5) lo--;
    while (hi + 1 < PR_RECTO_SAMPLES && values[hi + 1] >= 0.5) hi++;
    if (lo == 0 || hi == PR_RECTO_SAMPLES - 1) { s->broad++; return 0; }
    /* Reject equally near competing runs rather than choosing by loop order. */
    for (int i = 0; i < PR_RECTO_SAMPLES; i++) {
        double t = ((double)i - PR_RECTO_HALF_SAMPLES) * PR_RECTO_STEP;
        if ((i < lo || i > hi) && values[i] >= 0.5 && fabs(t) <= best + PR_RECTO_STEP) {
            s->ambiguous++; return 0;
        }
    }
    edge = ((double)hi - PR_RECTO_HALF_SAMPLES +
            (values[hi] - 0.5) / (values[hi] - values[hi + 1])) * PR_RECTO_STEP;
    outer = ((double)lo - 1.0 - PR_RECTO_HALF_SAMPLES +
             (0.5 - values[lo - 1]) / (values[lo] - values[lo - 1])) * PR_RECTO_STEP;
    if (edge - outer > PR_RECTO_MAX_THICKNESS || fabs(edge) > PR_RECTO_MAX_MOVE) {
        s->broad++; return 0;
    }
    *distance = (float)edge;
    s->supported++;
    s->sum_abs += fabs(edge);
    s->sum_signed += edge;
    if (fabs(edge) > s->max_abs) s->max_abs = fabs(edge);
    return 1;
}

void PredRecto_cloud(const PredRectoField *field, float *verts,
                     const float *normals, size_t nv, PredRectoStats *stats)
{
    for (size_t i = 0; i < nv; i++) {
        float inward[3] = {0.0f, 0.0f, 0.0f}, distance = 0.0f;
        if (PredRecto_target(field, verts + 3*i, normals + 3*i, inward, &distance, stats))
            for (int a = 0; a < 3; a++) verts[3*i + (size_t)a] += distance * inward[a];
    }
}

void PredRecto_report(const char *stage, const PredRectoStats *s)
{
    fprintf(stderr, "  F0-inner %s: tested=%zu supported=%zu outside=%zu ambiguous=%zu missing=%zu broad=%zu "
                    "reverted=%zu mean_abs=%.4f mean_signed=%.4f max_abs=%.4f\n",
            stage, s->tested, s->supported, s->outside, s->ambiguous, s->missing, s->broad, s->reverted,
            s->supported ? s->sum_abs / (double)s->supported : 0.0,
            s->supported ? s->sum_signed / (double)s->supported : 0.0, s->max_abs);
}

static void face_normal(const float *v, const int32_t *f, double n[3])
{
    double u[3] = {0.0}, w[3] = {0.0};
    for (int a = 0; a < 3; a++) {
        u[a] = (double)v[3*(size_t)f[1] + (size_t)a] - v[3*(size_t)f[0] + (size_t)a];
        w[a] = (double)v[3*(size_t)f[2] + (size_t)a] - v[3*(size_t)f[0] + (size_t)a];
    }
    for (int a = 0; a < 3; a++) n[a] = u[(a+1)%3]*w[(a+2)%3] - u[(a+2)%3]*w[(a+1)%3];
}

void PredRecto_meshes(Arena_T arena, const PredRectoField *field,
                      ComponentMesh *meshes, size_t count, const char *stage)
{
    PredRectoStats stats = {0};
    for (size_t c = 0; c < count; c++) {
        ComponentMesh *m = meshes + c;
        Arena_Mark mark = Arena_save(arena);
        float *normals = (float *)ARENA_CALLOC(arena, m->nv*3, sizeof(float));
        float *next = (float *)ARENA_ALLOC(arena, m->nv*3*sizeof(float));
        size_t supported_before = stats.supported;
        int safe = 1;
        memcpy(next, m->verts, m->nv*3*sizeof(float));
        for (size_t j = 0; j < m->nf; j++) {
            const int32_t *f = m->faces + j*3;
            double n[3] = {0.0};
            face_normal(m->verts, f, n);
            for (int k = 0; k < 3; k++)
                for (int a = 0; a < 3; a++) normals[3*(size_t)f[k] + (size_t)a] += (float)n[a];
        }
        PredRecto_cloud(field, next, normals, m->nv, &stats);
        for (size_t j = 0; j < m->nf && safe; j++) {
            double a[3] = {0.0}, b[3] = {0.0}, dot = 0.0, aa = 0.0;
            face_normal(m->verts, m->faces + 3*j, a);
            face_normal(next, m->faces + 3*j, b);
            for (int k = 0; k < 3; k++) { dot += a[k]*b[k]; aa += a[k]*a[k]; }
            if (!isfinite(dot) || dot < PR_RECTO_MIN_FACE_RATIO*aa) safe = 0;
        }
        if (safe) {
            memcpy(m->verts, next, m->nv*3*sizeof(float));
            m->vert_normals = NULL;
        } else stats.reverted += stats.supported - supported_before;
        Arena_restore(arena, mark);
    }
    PredRecto_report(stage, &stats);
}

int PredRecto_selftest(void)
{
    Arena_T arena = Arena_new();
    PredRectoField field = {0};
    uint8_t *v = (uint8_t *)ARENA_CALLOC(arena, 32*32*32, sizeof(uint8_t));
    float p[3] = {16.5f, 16.5f, 16.5f}, n[3] = {0.0f, 0.0f, -1.0f};
    float inward[3] = {0.0f}, d = 0.0f;
    PredRectoStats s = {0};
    int failures = 0;
    field.vol = v; field.d = field.h = field.w = 32;
    field.axis.umb_y = 16.5; field.axis.umb_x = 0.0;
    for (int z = 0; z < 32; z++) for (int y = 0; y < 32; y++)
        for (int x = 15; x <= 17; x++) v[(z*32+y)*32+x] = 1;
    /* A three-voxel layer: centres 15.5..17.5, inward boundary 15.0. */
    if (!PredRecto_target(&field, p, n, inward, &d, &s) || fabs(d - 1.5f) > 1e-5) failures++;
    n[2] = 1.0f;
    if (!PredRecto_target(&field, p, n, inward, &d, &s) || fabs(d - 1.5f) > 1e-5 || inward[2] > 0.0f) failures++;
    p[2] = 15.0f;
    if (!PredRecto_target(&field, p, n, inward, &d, &s) || fabs(d) > 1e-5) failures++;
    /* Inward and outward edges differ; repeated placement is idempotent. */
    field.axis.umb_x = 32.0;
    p[2] = 16.5f;
    if (!PredRecto_target(&field, p, n, inward, &d, &s) || fabs(d - 1.5f) > 1e-5 || inward[2] < 0.0f) failures++;
    p[2] = 3.0f;
    if (PredRecto_target(&field, p, n, inward, &d, &s)) failures++;
    p[2] = 8.0f;
    if (PredRecto_target(&field, p, n, inward, &d, &s)) failures++;
    p[2] = 16.5f; n[0] = 1.0f; n[2] = 0.0f;
    if (PredRecto_target(&field, p, n, inward, &d, &s)) failures++;
    memset(v, 1, 32*32*32); n[0] = 0.0f; n[2] = 1.0f;
    if (PredRecto_target(&field, p, n, inward, &d, &s)) failures++;
    PredRecto_cloud(&field, NULL, NULL, 0, &s);
    PredRecto_meshes(arena, &field, NULL, 0, "empty-selftest");
    Arena_dispose(&arena);
    fprintf(stderr, "[pred_recto selftest] %s (%d failures)\n", failures ? "FAIL" : "PASS", failures);
    return failures;
}
