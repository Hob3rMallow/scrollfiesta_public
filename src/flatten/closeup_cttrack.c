#include "closeup_cttrack.h"

#include "../common/arena.h"
#include "../common/eig3.h"
#include "../common/kdtree.h"
#include "../common/mesh_bin.h"
#include "../common/pipeline_constants.h"
#include "../common/raw_sample.h"
#include "../common/ves_platform.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double (*CctSample)(void *ctx, const double p[3]);

typedef struct CctPoint {
    double p[3], n[3];
    double profile[CLOSEUP_CTTRACK_PROFILE_N];
    double corr, rms, shift, cycle;
} CctPoint;

typedef struct CctSeed {
    CctPoint point;
    double uv[2];
    size_t face;
} CctSeed;

enum { CCT_OK = 0, CCT_RAW = 1, CCT_TANGENT = 2, CCT_TENSOR = 3,
       CCT_PROFILE = 4, CCT_MATCH = 5, CCT_CYCLE = 6, CCT_BOX = 7 };

static double cct_dot(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static int cct_unit(double a[3])
{
    double l = sqrt(cct_dot(a, a));
    if (!(l > 1e-12) || !isfinite(l)) return -1;
    for (int k = 0; k < 3; k++) a[k] /= l;
    return 0;
}

static void cct_cross(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/* Project the requested axial direction into the local tangent plane. */
static int cct_basis(const double n[3], int direction, double t[3], double u[3])
{
    for (int k = 0; k < 3; k++)
        t[k] = (double)direction * ((k == 0 ? 1.0 : 0.0) - n[0] * n[k]);
    if (sqrt(cct_dot(t, t)) < CLOSEUP_CTTRACK_AXIAL_MIN || cct_unit(t)) return -1;
    cct_cross(n, t, u);
    return cct_unit(u);
}

/* A normal profile averaged over a 3x3 tangent neighborhood. */
static int cct_profile(CctSample sample, void *ctx, CctPoint *point)
{
    double t[3] = {0}, u[3] = {0}, mean = 0.0, ss = 0.0;
    if (cct_basis(point->n, 1, t, u)) return CCT_TANGENT;
    for (int j = 0; j < CLOSEUP_CTTRACK_PROFILE_N; j++) {
        double d = (double)(j - CLOSEUP_CTTRACK_PROFILE_N / 2) * CLOSEUP_CTTRACK_PROFILE_STEP;
        double sum = 0.0;
        for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) {
            double p[3] = {0};
            for (int k = 0; k < 3; k++)
                p[k] = point->p[k] + d * point->n[k] + CLOSEUP_CTTRACK_PATCH_STEP * ((double)a * t[k] + (double)b * u[k]);
            double value = sample(ctx, p);
            if (!isfinite(value) || value < 0.0) return CCT_RAW;
            sum += value;
        }
        point->profile[j] = sum / 9.0;
        mean += point->profile[j];
    }
    mean /= (double)CLOSEUP_CTTRACK_PROFILE_N;
    for (int j = 0; j < CLOSEUP_CTTRACK_PROFILE_N; j++) {
        double d = point->profile[j] - mean;
        ss += d * d;
    }
    point->rms = sqrt(ss / (double)CLOSEUP_CTTRACK_PROFILE_N);
    return point->rms >= CLOSEUP_CTTRACK_PROFILE_RMS ? CCT_OK : CCT_PROFILE;
}

static double cct_corr(const CctPoint *a, const CctPoint *b)
{
    double ma = 0.0, mb = 0.0, aa = 0.0, bb = 0.0, ab = 0.0;
    for (int j = 0; j < CLOSEUP_CTTRACK_PROFILE_N; j++) { ma += a->profile[j]; mb += b->profile[j]; }
    ma /= (double)CLOSEUP_CTTRACK_PROFILE_N; mb /= (double)CLOSEUP_CTTRACK_PROFILE_N;
    for (int j = 0; j < CLOSEUP_CTTRACK_PROFILE_N; j++) {
        double x = a->profile[j] - ma, y = b->profile[j] - mb;
        aa += x * x; bb += y * y; ab += x * y;
    }
    return aa > 1e-12 && bb > 1e-12 ? ab / sqrt(aa * bb) : -1.0;
}

/* The dominant CT-gradient direction supplies a local normal. The small
 * neighborhood and per-step angular gate prevent wholesale ridge switching;
 * they are not a proof of layer identity, which the withheld audit must test. */
static int cct_normal(CctSample sample, void *ctx, const double p[3],
                       const double old[3], double out[3])
{
    double t[3] = {0}, u[3] = {0}, tensor[3][3] = {{0}}, eval[3] = {0}, evec[3][3] = {{0}};
    if (cct_basis(old, 1, t, u)) return CCT_TANGENT;
    for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) for (int c = -1; c <= 1; c++) {
        double q[3] = {0}, grad[3] = {0};
        for (int k = 0; k < 3; k++) q[k] = p[k] + CLOSEUP_CTTRACK_TENSOR_STEP *
            ((double)a * t[k] + (double)b * u[k] + (double)c * old[k]);
        for (int k = 0; k < 3; k++) {
            double lo[3] = {q[0], q[1], q[2]}, hi[3] = {q[0], q[1], q[2]};
            lo[k] -= 0.5; hi[k] += 0.5;
            double x = sample(ctx, lo), y = sample(ctx, hi);
            if (!isfinite(x) || !isfinite(y) || x < 0.0 || y < 0.0) return CCT_RAW;
            grad[k] = y - x;
        }
        for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) tensor[k][l] += grad[k] * grad[l];
    }
    Eig3_sym(tensor, eval, evec);
    if (!(eval[2] > 1e-8) || (eval[2] - eval[1]) / eval[2] < CLOSEUP_CTTRACK_COHERENCE) return CCT_TENSOR;
    for (int k = 0; k < 3; k++) out[k] = evec[k][2];
    if (cct_dot(out, old) < 0.0) for (int k = 0; k < 3; k++) out[k] = -out[k];
    if (cct_dot(out, old) < CLOSEUP_CTTRACK_NORMAL_COS) return CCT_TENSOR;
    return CCT_OK;
}

static int cct_advance(CctSample sample, void *ctx, const CctPoint *old,
                        const CctPoint *anchor, int direction, CctPoint *out)
{
    double t[3] = {0}, u[3] = {0}, prediction[3] = {0}, n[3] = {0}, best = -1e100;
    CctPoint candidate = {0};
    int status = CCT_MATCH, found = 0;
    if (cct_basis(old->n, direction, t, u)) return CCT_TANGENT;
    for (int k = 0; k < 3; k++) prediction[k] = old->p[k] + CLOSEUP_CTTRACK_STEP * t[k];
    status = cct_normal(sample, ctx, prediction, old->n, n);
    if (status != CCT_OK) return status;
    for (int k = 0; k < 3; k++) candidate.n[k] = n[k];
    for (int j = -CLOSEUP_CTTRACK_SEARCH_N; j <= CLOSEUP_CTTRACK_SEARCH_N; j++) {
        double d = (double)j * CLOSEUP_CTTRACK_SEARCH_STEP;
        for (int k = 0; k < 3; k++) candidate.p[k] = prediction[k] + d * n[k];
        if ((double)direction * (candidate.p[0] - old->p[0]) <= 0.0) continue;
        if (cct_profile(sample, ctx, &candidate) != CCT_OK) continue;
        double current = cct_corr(old, &candidate), initial = cct_corr(anchor, &candidate);
        double score = 0.5 * (current + initial) - CLOSEUP_CTTRACK_OFFSET_PRIOR * d * d;
        if (current < CLOSEUP_CTTRACK_CORR || initial < CLOSEUP_CTTRACK_CORR ||
            candidate.rms < anchor->rms * CLOSEUP_CTTRACK_GAIN_MIN ||
            candidate.rms > anchor->rms * CLOSEUP_CTTRACK_GAIN_MAX) continue;
        if (score > best) {
            best = score; candidate.corr = fmin(current, initial); candidate.shift = d;
            *out = candidate; found = 1;
        }
    }
    return found ? CCT_OK : CCT_MATCH;
}

static int cct_step(CctSample sample, void *ctx, const CctPoint *old,
                     const CctPoint *anchor, CctPoint *out)
{
    CctPoint reverse = {0};
    int rc = cct_advance(sample, ctx, old, anchor, -1, out);
    if (rc != CCT_OK) return rc;
    rc = cct_advance(sample, ctx, out, anchor, 1, &reverse);
    if (rc != CCT_OK) return CCT_CYCLE;
    double d[3] = {0};
    for (int k = 0; k < 3; k++) d[k] = reverse.p[k] - old->p[k];
    out->cycle = sqrt(cct_dot(d, d));
    return out->cycle <= CLOSEUP_CTTRACK_CYCLE ? CCT_OK : CCT_CYCLE;
}

static int cct_inside(const double p[3], const double box[6])
{
    for (int k = 0; k < 3; k++) if (!(p[k] >= box[2 * k] && p[k] <= box[2 * k + 1])) return 0;
    return 1;
}

static double cct_raw(void *ctx, const double p[3])
{
    return sample_trilinear((CubeTable *)ctx, p[0], p[1], p[2]);
}

static int cct_read(const char *path, void *value, size_t size)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    int rc = fread(value, 1, size, f) == size && fgetc(f) == EOF ? 0 : -1;
    fclose(f);
    return rc;
}

static int cct_seed_cmp(const void *va, const void *vb)
{
    const CctSeed *a = (const CctSeed *)va, *b = (const CctSeed *)vb;
    if (a->uv[0] != b->uv[0]) return a->uv[0] < b->uv[0] ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

/* Triangle/plane intersection uses only observed corners, before anything
 * withheld is read by the tracker. Sampling is uniform in the existing U. */
static int cct_seed_face(const MeshBinData *m, const float *normal, const uint8_t *attr,
                          size_t f, double z, const double box[6], CctSeed *seed)
{
    int count = 0;
    memset(seed, 0, sizeof *seed);
    for (int k = 0; k < 3; k++) {
        size_t v = (size_t)m->faces[3 * f + (size_t)k];
        if (attr[v] != 1 && attr[v] != 3) return 0;
    }
    for (int k = 0; k < 3; k++) {
        size_t a = (size_t)m->faces[3 * f + (size_t)k], b = (size_t)m->faces[3 * f + (size_t)((k + 1) % 3)];
        double za = m->verts[a * 3], zb = m->verts[b * 3];
        if ((za <= z && zb > z) || (zb <= z && za > z)) {
            double t = (z - za) / (zb - za);
            for (int d = 0; d < 3; d++) {
                seed->point.p[d] += (1.0 - t) * m->verts[a * 3 + (size_t)d] + t * m->verts[b * 3 + (size_t)d];
                seed->point.n[d] += (1.0 - t) * normal[a * 3 + (size_t)d] + t * normal[b * 3 + (size_t)d];
            }
            for (int d = 0; d < 2; d++) seed->uv[d] += (1.0 - t) * MeshBin_uv(m, a * 2 + (size_t)d) + t * MeshBin_uv(m, b * 2 + (size_t)d);
            count++;
        }
    }
    if (count != 2) return 0;
    for (int d = 0; d < 3; d++) seed->point.p[d] *= 0.5;
    for (int d = 0; d < 2; d++) seed->uv[d] *= 0.5;
    if (!cct_inside(seed->point.p, box) || cct_unit(seed->point.n)) return 0;
    seed->face = f;
    return 1;
}

int CloseupCtTrack_audit(const char *surface_stem, const char *raw,
                         const char *out_stem, const double box[6], double seed_z)
{
    Arena_T arena = Arena_new();
    MeshBinData m = {0};
    CubeTable ct = {0};
    KDTree_T truth = NULL;
    char path[2600];
    CctSeed *seeds = NULL, candidate = {0};
    float *normal = NULL, *observed = NULL, bounds[6] = {0};
    uint8_t *attr = NULL;
    FILE *csv = NULL, *report = NULL;
    size_t nseeds = 0, ntruth = 0, selected = 0, accepted = 0, stops[8] = {0};
    double start = ves_clock_sec(), last_u = -1e100;
    int rc = -1;
    if (surface_stem == NULL || raw == NULL || out_stem == NULL || box == NULL || !isfinite(seed_z)) goto out;
    for (int k = 0; k < 3; k++) if (!isfinite(box[2 * k]) || !isfinite(box[2 * k + 1]) || box[2 * k + 1] <= box[2 * k]) goto out;
    if (seed_z < box[0] || seed_z > box[1]) goto out;
    snprintf(path, sizeof path, "%s_surface.vmesh", surface_stem);
    if (MeshBin_read_precise_arena(arena, path, &m) != 0 || m.nv == 0 || (m.uv == NULL && m.uv64 == NULL)) goto out;
    normal = (float *)ARENA_ALLOC(arena, m.nv * 3 * sizeof *normal);
    attr = (uint8_t *)ARENA_ALLOC(arena, m.nv);
    snprintf(path, sizeof path, "%s_surface_nrm.f32", surface_stem);
    if (cct_read(path, normal, m.nv * 3 * sizeof *normal)) goto out;
    snprintf(path, sizeof path, "%s_surface_attr.u8", surface_stem);
    if (cct_read(path, attr, m.nv)) goto out;
    for (size_t f = 0; f < m.nf; f++) nseeds += (size_t)cct_seed_face(&m, normal, attr, f, seed_z, box, &candidate);
    if (nseeds == 0) { fprintf(stderr, "[cttrack] no observed seed intersections\n"); goto out; }
    seeds = (CctSeed *)ARENA_ALLOC(arena, nseeds * sizeof *seeds);
    nseeds = 0;
    for (size_t f = 0; f < m.nf; f++) if (cct_seed_face(&m, normal, attr, f, seed_z, box, &candidate)) seeds[nseeds++] = candidate;
    qsort(seeds, nseeds, sizeof *seeds, cct_seed_cmp);
    for (size_t v = 0; v < m.nv; v++) {
        double p[3] = {m.verts[v * 3], m.verts[v * 3 + 1], m.verts[v * 3 + 2]};
        if ((attr[v] == 1 || attr[v] == 3) && cct_inside(p, box)) ntruth++;
    }
    if (ntruth == 0) goto out;
    observed = (float *)ARENA_ALLOC(arena, ntruth * 3 * sizeof *observed);
    ntruth = 0;
    for (size_t v = 0; v < m.nv; v++) {
        double p[3] = {m.verts[v * 3], m.verts[v * 3 + 1], m.verts[v * 3 + 2]};
        if ((attr[v] == 1 || attr[v] == 3) && cct_inside(p, box)) {
            for (int k = 0; k < 3; k++) observed[ntruth * 3 + (size_t)k] = m.verts[v * 3 + (size_t)k];
            ntruth++;
        }
    }
    truth = KDTree_new(arena, observed, ntruth);
    for (int k = 0; k < 3; k++) { bounds[k] = (float)box[2 * k]; bounds[k + 3] = (float)box[2 * k + 1]; }
    if (cubetable_init(&ct, arena, raw, 128, bounds, 2, CLOSEUP_CTTRACK_RAW_PAD)) goto out;
    cubetable_prewarm_all(&ct);
    if (!cubetable_is_complete(&ct)) { fprintf(stderr, "[cttrack] incomplete raw CT: %d missing\n", ct.n_missing); goto out; }
    snprintf(path, sizeof path, "%s_tracks.csv", out_stem);
    ves_ensure_parent_dir(path);
    csv = fopen(path, "wb");
    if (csv == NULL) goto out;
    fprintf(csv, "seed,step,seed_u,seed_v,z,y,x,nz,ny,nx,corr,profile_rms,offset,cycle_error,nearest_observed_vertex,status\n");
    for (size_t s = 0; s < nseeds; s++) {
        CctPoint current = seeds[s].point, anchor = {0}, next = {0};
        int stop = CCT_OK, steps = 0;
        if (seeds[s].uv[0] < last_u + CLOSEUP_CTTRACK_SEED_PX) continue;
        last_u = seeds[s].uv[0];
        stop = cct_profile(cct_raw, &ct, &current);
        anchor = current;
        current.corr = 1.0;
        for (int step = 0; step <= CLOSEUP_CTTRACK_STEPS; step++) {
            float query[3] = {(float)current.p[0], (float)current.p[1], (float)current.p[2]}, d2 = 0.0f;
            (void)KDTree_nearest(truth, query, &d2);
            fprintf(csv, "%zu,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%d\n",
                selected, step, seeds[s].uv[0], seeds[s].uv[1], current.p[0], current.p[1], current.p[2],
                current.n[0], current.n[1], current.n[2], current.corr, current.rms, current.shift, current.cycle, sqrt((double)d2), stop);
            if (stop != CCT_OK || step == CLOSEUP_CTTRACK_STEPS) break;
            stop = cct_step(cct_raw, &ct, &current, &anchor, &next);
            if (stop == CCT_OK && !cct_inside(next.p, box)) stop = CCT_BOX;
            if (stop == CCT_OK) { current = next; steps++; accepted++; }
        }
        stops[stop]++;
        fprintf(stderr, "[cttrack] seed %zu at u %.1f: %d steps, stop %d, z %.3f\n", selected, last_u, steps, stop, current.p[0]);
        selected++;
    }
    fclose(csv); csv = NULL;
    snprintf(path, sizeof path, "%s_report.json", out_stem);
    report = fopen(path, "wb");
    if (report == NULL) goto out;
    fprintf(report, "{\n  \"schema\": \"native-closeup-cttrack-audit-v1\",\n  \"seed_z\": %.8f,\n"
        "  \"box_zyx\": [%.6f,%.6f,%.6f,%.6f,%.6f,%.6f],\n  \"seed_intersections\": %zu,\n"
        "  \"selected_seeds\": %zu,\n  \"accepted_steps\": %zu,\n  \"raw_loaded\": %d,\n  \"raw_missing\": %d,\n"
        "  \"observed_reference_vertices\": %zu,\n  \"reference_is_used_only_after_tracking\": true,\n"
        "  \"surface_modified\": false,\n  \"stops\": {", seed_z, box[0], box[1], box[2], box[3], box[4], box[5],
        nseeds, selected, accepted, ct.n_loaded, ct.n_missing, ntruth);
    for (int k = 0; k < 8; k++) fprintf(report, "%s\"%d\":%zu", k ? "," : "", k, stops[k]);
    fprintf(report, "},\n  \"seconds\": %.3f\n}\n", ves_clock_sec() - start);
    rc = 0;
out:
    if (csv != NULL) fclose(csv);
    if (report != NULL) fclose(report);
    Arena_dispose(&arena);
    return rc;
}

typedef struct CctSynthetic { int bend, missing; } CctSynthetic;

static double cct_synthetic(void *ctx, const double p[3])
{
    const CctSynthetic *c = (const CctSynthetic *)ctx;
    double mid = c->bend ? 0.0025 * p[0] * p[0] : 0.0;
    double d = p[2] - mid, adjacent = p[2] - mid - 7.0;
    double target = c->missing && p[0] < -16.0 ? 0.0 : 160.0 * exp(-0.5 * d * d);
    return 20.0 + target + 150.0 * exp(-0.5 * adjacent * adjacent);
}

int CloseupCtTrack_selftest(void)
{
    int failures = 0;
    for (int test = 0; test < 3; test++) {
        CctSynthetic sample = {test == 1, test == 2};
        CctPoint current = {0}, anchor = {0}, next = {0};
        double max_error = 0.0;
        int stop = 0, steps = 0;
        current.n[2] = 1.0;
        stop = cct_profile(cct_synthetic, &sample, &current);
        anchor = current;
        if (stop != CCT_OK) { failures++; continue; }
        for (int s = 0; s < CLOSEUP_CTTRACK_STEPS; s++) {
            stop = cct_step(cct_synthetic, &sample, &current, &anchor, &next);
            if (stop != CCT_OK) break;
            current = next; steps++;
            double expected = sample.bend ? 0.0025 * current.p[0] * current.p[0] : 0.0;
            max_error = fmax(max_error, fabs(current.p[2] - expected));
        }
        int bad = test < 2 ? (steps != CLOSEUP_CTTRACK_STEPS || max_error > 0.5) :
            (steps >= 20 || steps < 10 || fabs(current.p[2]) > 0.5 || stop == CCT_OK);
        fprintf(stderr, "  cttrack synthetic %d: steps %d stop %d max error %.4f %s\n", test, steps, stop, max_error, bad ? "FAIL" : "PASS");
        failures += bad;
    }
    return failures ? 1 : 0;
}
