#include "recto_refine.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "../common/csr.h"
#include "../common/kdtree.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"

#define RR_MAX_HITS 16

static double rr_clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

void RectoRefineOpts_default(RectoRefineOpts *o)
{
    assert(o);
    memset(o, 0, sizeof *o);
    o->axis_point[1] = 3405.0f;
    o->axis_point[2] = 2878.0f;
    o->axis_dir[0] = 1.0f;
    o->window_lo = 0.0;
    o->window_hi = 255.0;
    o->band_frac = 0.35;
    o->range = 3.0;
    o->sample_step = 0.25;
    o->edge_half = 0.5;
    o->min_edge_frac = 0.08;
    o->distance_bias = 4.0;
    o->normal_blend = 0.85;
    o->data_weight = 12.0;
    o->anchor_weight = 24.0;
    o->smooth_mu = 1.0;
    o->max_slope = 0.55;
    o->damping = 0.70;
    o->max_step = 0.85;
    o->max_total = 3.0;
    o->occ_thresh = 0.9;
    o->local_r = 2.0;
    o->outer_iters = 4;
    o->inner_iters = 48;
    o->tol = 0.02;
    o->threads = 32;
}

/* Stable outward direction: radial defines the recto/verso sense; the mesh
 * normal, flipped outward, improves the local across-sheet line without being
 * allowed to reverse that sense. */
static int rr_direction(const float *p, const float *normal,
                        const RectoRefineOpts *o, double out[3])
{
    double axis[3] = { o->axis_dir[0], o->axis_dir[1], o->axis_dir[2] };
    double al = sqrt(axis[0] * axis[0] + axis[1] * axis[1]
                     + axis[2] * axis[2]);
    double d[3], radial[3], axial = 0.0, rl = 0.0;
    double nn = 0.0, blend = rr_clamp(o->normal_blend, 0.0, 1.0);
    int k = 0;

    if (al < 1e-12) {
        axis[0] = 1.0;
        axis[1] = axis[2] = 0.0;
    } else {
        for (k = 0; k < 3; k++) axis[k] /= al;
    }
    for (k = 0; k < 3; k++) {
        d[k] = (double)p[k] - (double)o->axis_point[k];
        axial += d[k] * axis[k];
    }
    for (k = 0; k < 3; k++) {
        radial[k] = d[k] - axial * axis[k];
        rl += radial[k] * radial[k];
    }
    rl = sqrt(rl);
    if (rl > 1e-9)
        for (k = 0; k < 3; k++) radial[k] /= rl;

    if (normal != NULL)
        nn = sqrt((double)normal[0] * normal[0]
                  + (double)normal[1] * normal[1]
                  + (double)normal[2] * normal[2]);
    if (rl <= 1e-9 && nn <= 0.5) {
        out[0] = out[1] = out[2] = 0.0;
        return -1;
    }
    if (rl <= 1e-9) {
        for (k = 0; k < 3; k++) out[k] = (double)normal[k] / nn;
        return 0;
    }
    if (nn > 0.5 && blend > 0.0) {
        double nu[3], dot = 0.0, len = 0.0;
        for (k = 0; k < 3; k++) {
            nu[k] = (double)normal[k] / nn;
            dot += nu[k] * radial[k];
        }
        if (dot < 0.0)
            for (k = 0; k < 3; k++) nu[k] = -nu[k];
        for (k = 0; k < 3; k++) {
            out[k] = (1.0 - blend) * radial[k] + blend * nu[k];
            len += out[k] * out[k];
        }
        len = sqrt(len);
        if (len > 1e-9) {
            for (k = 0; k < 3; k++) out[k] /= len;
            return 0;
        }
    }
    for (k = 0; k < 3; k++) out[k] = radial[k];
    return 0;
}

/* Positive derivative along OUTWARD is the inner-facing (recto) edge.  The
 * opposite boundary of the same bright band has a negative derivative and
 * therefore cannot win, regardless of its peak intensity. */
static int rr_find_target(CubeTable *ct, const double p[3],
                          const double dir[3], const RectoRefineOpts *o,
                          double band, double min_edge,
                          double *out_t, double *out_grad)
{
    double range = o->range > 0.0 ? o->range : 3.0;
    double step = o->sample_step > 0.0 ? o->sample_step : 0.25;
    double half = o->edge_half > 0.0 ? o->edge_half : 0.5;
    int ns = (int)floor((2.0 * range) / step + 0.5) + 1;
    double best_score = -1e30, best_t = 0.0, best_grad = 0.0;
    int found = 0, si = 0;

    if (ns < 3) ns = 3;
    if (ns > 1025) ns = 1025;
    for (si = 0; si < ns; si++) {
        double t = -range + (2.0 * range * (double)si) / (double)(ns - 1);
        double qm[3], qp[3], sm, sp, grad, score, tc = t;
        int crossing = 0, k = 0;

        for (k = 0; k < 3; k++) {
            qm[k] = p[k] + (t - half) * dir[k];
            qp[k] = p[k] + (t + half) * dir[k];
        }
        sm = sample_trilinear(ct, qm[0], qm[1], qm[2]);
        sp = sample_trilinear(ct, qp[0], qp[1], qp[2]);
        if (sm < 0.0 || sp < 0.0) continue;
        grad = sp - sm;
        if (grad < min_edge) continue;
        crossing = sm <= band && sp >= band && sp > sm + 1e-9;
        if (crossing) {
            double f = rr_clamp((band - sm) / (sp - sm), 0.0, 1.0);
            tc = t - half + 2.0 * half * f;
        }
        score = grad - o->distance_bias * fabs(tc);
        if (!crossing) {
            double mismatch = fabs(0.5 * (sm + sp) - band);
            score -= 0.5 * mismatch;
        }
        if (!found || score > best_score) {
            best_score = score;
            best_t = tc;
            best_grad = grad;
            found = 1;
        }
    }
    *out_t = best_t;
    *out_grad = best_grad;
    return found;
}

static int rr_is_one_ring(const int32_t *off, const int32_t *tgt,
                          size_t i, size_t j)
{
    if (i == j) return 1;
    for (int32_t e = off[i]; e < off[i + 1]; e++)
        if ((size_t)tgt[e] == j) return 1;
    return 0;
}

static int rr_hits_other(KDTree_T occ, const float *start,
                         const int32_t *off, const int32_t *tgt, size_t i,
                         const float q[3], double occ2, double local2)
{
    int32_t hit[RR_MAX_HITS];
    size_t nh = KDTree_ball_query(occ, q, (float)occ2, hit, RR_MAX_HITS);
    size_t h = 0;

    if (nh > RR_MAX_HITS) nh = RR_MAX_HITS;
    for (h = 0; h < nh; h++) {
        size_t j = (size_t)hit[h];
        double dz = (double)start[j * 3 + 0] - start[i * 3 + 0];
        double dy = (double)start[j * 3 + 1] - start[i * 3 + 1];
        double dx = (double)start[j * 3 + 2] - start[i * 3 + 2];
        double dsq = dz * dz + dy * dy + dx * dx;
        if (dsq > local2 || !rr_is_one_ring(off, tgt, i, j)) return 1;
    }
    return 0;
}

int RectoRefine_run(Arena_T arena,
                    CubeTable *ct,
                    float *verts,
                    size_t nv,
                    const int32_t *faces,
                    size_t nf,
                    const RectoRefineOpts *opts,
                    RectoRefineStats *out)
{
    Arena_Mark mark;
    CSR_T adj;
    const int32_t *off = NULL, *tgt = NULL;
    KDTree_T occ;
    float *start = NULL, *dirs = NULL, *target = NULL, *weight = NULL;
    float *base_disp = NULL;
    float *disp_a = NULL, *disp_b = NULL;
    uint8_t *ever_supported = NULL, *slope_hit = NULL;
    RectoRefineOpts o;
    double contrast = 0.0, band = 0.0, min_edge = 0.0;
    double occ2 = 0.0, local2 = 0.0, t0 = 0.0;
    int outer = 0, nthreads = 0;

    assert(arena && ct && verts && faces && opts && out);
    memset(out, 0, sizeof *out);
    if (nv == 0 || nf == 0) return 0;
    if (nv > (size_t)INT_MAX
        || nv > (size_t)LONG_MAX / (3 * sizeof(float))) return -1;

    o = *opts;
    contrast = o.window_hi - o.window_lo;
    if (contrast <= 1e-6) {
        o.window_lo = 0.0;
        o.window_hi = 255.0;
        contrast = 255.0;
    }
    band = o.window_lo + rr_clamp(o.band_frac, 0.05, 0.95) * contrast;
    min_edge = (o.min_edge_frac > 0.0 ? o.min_edge_frac : 0.08) * contrast;
    occ2 = o.occ_thresh > 0.0 ? o.occ_thresh : 0.9;
    occ2 *= occ2;
    local2 = o.local_r > 0.0 ? o.local_r : 2.0;
    local2 *= local2;
    nthreads = o.threads > 0 ? o.threads : 32;
    t0 = ves_clock_sec();

    /* Pre-warm before the scratch mark: cube buffers belong to CubeTable and
     * must outlive this call; restoring the scratch mark must not invalidate
     * slots loaded by concurrent profile sampling. */
    (void)cubetable_prewarm_all(ct);
    mark = Arena_save(arena);
    adj = CSR_from_faces(arena, faces, nf, nv);
    off = CSR_offset(adj);
    tgt = CSR_target(adj);
    start = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    memcpy(start, verts, nv * 3 * sizeof(float));
    occ = KDTree_new(arena, start, nv);
    dirs = (float *)ARENA_CALLOC(arena, (long)nv, 3 * sizeof(float));
    target = (float *)ARENA_CALLOC(arena, (long)nv, sizeof(float));
    weight = (float *)ARENA_CALLOC(arena, (long)nv, sizeof(float));
    base_disp = (float *)ARENA_CALLOC(arena, (long)nv, sizeof(float));
    disp_a = (float *)ARENA_CALLOC(arena, (long)nv, sizeof(float));
    disp_b = (float *)ARENA_CALLOC(arena, (long)nv, sizeof(float));
    ever_supported = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);
    slope_hit = (uint8_t *)ARENA_CALLOC(arena, (long)nv, 1);

#ifdef _OPENMP
    ves_omp_set_threads(nthreads);
#else
    (void)nthreads;
#endif

    for (outer = 0; outer < (o.outer_iters > 0 ? o.outer_iters : 4); outer++) {
        float *normals = vertex_normals(verts, nv, faces, nf);
        long long supported = 0;
        double grad_sum = 0.0;
        int n = (int)nv, i = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:supported,grad_sum)
#endif
        for (i = 0; i < n; i++) {
            double p[3] = { verts[(size_t)i * 3 + 0],
                            verts[(size_t)i * 3 + 1],
                            verts[(size_t)i * 3 + 2] };
            double d[3] = { 0.0, 0.0, 0.0 }, toff = 0.0, grad = 0.0;
            int ok = 0;

            target[i] = 0.0f;
            base_disp[i] = 0.0f;
            weight[i] = (float)(o.anchor_weight > 0.0 ? o.anchor_weight : 24.0);
            if (off[i] == off[i + 1]
                || rr_direction(&verts[(size_t)i * 3],
                                &normals[(size_t)i * 3], &o, d) != 0) {
                dirs[(size_t)i * 3 + 0] = 0.0f;
                dirs[(size_t)i * 3 + 1] = 0.0f;
                dirs[(size_t)i * 3 + 2] = 0.0f;
                continue;
            }
            dirs[(size_t)i * 3 + 0] = (float)d[0];
            dirs[(size_t)i * 3 + 1] = (float)d[1];
            dirs[(size_t)i * 3 + 2] = (float)d[2];
            base_disp[i] = (float)(
                ((double)verts[(size_t)i * 3 + 0]
                 - start[(size_t)i * 3 + 0]) * d[0]
                + ((double)verts[(size_t)i * 3 + 1]
                   - start[(size_t)i * 3 + 1]) * d[1]
                + ((double)verts[(size_t)i * 3 + 2]
                   - start[(size_t)i * 3 + 2]) * d[2]);
            ok = rr_find_target(ct, p, d, &o, band, min_edge, &toff, &grad);
            if (ok) {
                double confidence = rr_clamp(
                    (grad - min_edge) / (0.25 * contrast + 1e-9), 0.0, 1.0);
                target[i] = (float)((double)base_disp[i] + toff);
                weight[i] = (float)((o.data_weight > 0.0 ? o.data_weight : 12.0)
                                    * (0.35 + 0.65 * confidence));
                ever_supported[i] = 1;
                supported++;
                grad_sum += grad;
            }
        }
        free(normals);
        out->n_supported = (size_t)supported;
        out->mean_recto_gradient =
            supported > 0 ? grad_sum / (double)supported : 0.0;

        memcpy(disp_a, base_disp, nv * sizeof(float));
        memcpy(disp_b, base_disp, nv * sizeof(float));
        for (int inner = 0;
             inner < (o.inner_iters > 0 ? o.inner_iters : 48);
             inner++) {
            float *src = (inner & 1) ? disp_b : disp_a;
            float *dst = (inner & 1) ? disp_a : disp_b;
            double mu = o.smooth_mu > 0.0 ? o.smooth_mu : 1.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (i = 0; i < n; i++) {
                int32_t e = 0;
                double sum = 0.0, val = 0.0;
                double lo = -1e30, hi = 1e30;
                int32_t degree = off[i + 1] - off[i];

                if (degree <= 0) {
                    dst[i] = 0.0f;
                    continue;
                }
                for (e = off[i]; e < off[i + 1]; e++) {
                    int32_t j = tgt[e];
                    double dz = (double)start[(size_t)j * 3 + 0]
                                - start[(size_t)i * 3 + 0];
                    double dy = (double)start[(size_t)j * 3 + 1]
                                - start[(size_t)i * 3 + 1];
                    double dx = (double)start[(size_t)j * 3 + 2]
                                - start[(size_t)i * 3 + 2];
                    double edge = sqrt(dz * dz + dy * dy + dx * dx);
                    double cap = (o.max_slope > 0.0 ? o.max_slope : 0.55)
                                 * (edge > 1e-3 ? edge : 1e-3);
                    double dj = src[j];
                    sum += dj;
                    if (dj - cap > lo) lo = dj - cap;
                    if (dj + cap < hi) hi = dj + cap;
                }
                val = ((double)weight[i] * target[i] + mu * sum)
                      / ((double)weight[i] + mu * (double)degree);
                if (lo <= hi) {
                    double clipped = rr_clamp(val, lo, hi);
                    if (fabs(clipped - val) > 1e-7) slope_hit[i] = 1;
                    val = clipped;
                } else {
                    val = sum / (double)degree;
                    slope_hit[i] = 1;
                }
                dst[i] = (float)val;
            }
        }

        {
            int inner_count = o.inner_iters > 0 ? o.inner_iters : 48;
            float *disp = (inner_count & 1) ? disp_b : disp_a;
            long long moved_step = 0, reverted = 0;
            double abs_step_sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    reduction(+:moved_step,reverted,abs_step_sum)
#endif
            for (i = 0; i < n; i++) {
                double ds = (o.damping > 0.0 ? o.damping : 0.70)
                          * ((double)disp[i] - (double)base_disp[i]);
                double max_step = o.max_step > 0.0 ? o.max_step : 0.85;
                double max_total = o.max_total > 0.0 ? o.max_total : 3.0;
                float q[3];
                double total[3], tl = 0.0;
                int k = 0;

                ds = rr_clamp(ds, -max_step, max_step);
                if (fabs(ds) < 1e-6) continue;
                for (k = 0; k < 3; k++)
                    q[k] = (float)((double)verts[(size_t)i * 3 + k]
                                   + ds * dirs[(size_t)i * 3 + k]);
                for (k = 0; k < 3; k++) {
                    total[k] = (double)q[k] - start[(size_t)i * 3 + k];
                    tl += total[k] * total[k];
                }
                tl = sqrt(tl);
                if (tl > max_total) {
                    double s = max_total / tl;
                    for (k = 0; k < 3; k++)
                        q[k] = (float)((double)start[(size_t)i * 3 + k]
                                      + s * total[k]);
                }
                if (!isfinite(q[0]) || !isfinite(q[1]) || !isfinite(q[2])
                    || rr_hits_other(occ, start, off, tgt, (size_t)i,
                                     q, occ2, local2)) {
                    reverted++;
                    continue;
                }
                {
                    double dz = (double)q[0] - verts[(size_t)i * 3 + 0];
                    double dy = (double)q[1] - verts[(size_t)i * 3 + 1];
                    double dx = (double)q[2] - verts[(size_t)i * 3 + 2];
                    double actual = sqrt(dz * dz + dy * dy + dx * dx);
                    verts[(size_t)i * 3 + 0] = q[0];
                    verts[(size_t)i * 3 + 1] = q[1];
                    verts[(size_t)i * 3 + 2] = q[2];
                    moved_step++;
                    abs_step_sum += actual;
                }
            }
            out->n_reverted += (size_t)reverted;
            out->mean_abs_step =
                moved_step > 0 ? abs_step_sum / (double)moved_step : 0.0;
            out->iterations = outer + 1;
            if (moved_step == 0 || out->mean_abs_step < o.tol) break;
        }
    }

    {
        double sum = 0.0, maxd = 0.0;
        size_t i = 0;
        out->n_supported = 0;
        for (i = 0; i < nv; i++) {
            double dz = (double)verts[i * 3 + 0] - start[i * 3 + 0];
            double dy = (double)verts[i * 3 + 1] - start[i * 3 + 1];
            double dx = (double)verts[i * 3 + 2] - start[i * 3 + 2];
            double d = sqrt(dz * dz + dy * dy + dx * dx);
            if (ever_supported[i]) out->n_supported++;
            if (slope_hit[i]) out->n_slope_limited++;
            if (d > 1e-5) {
                out->n_moved++;
                sum += d;
                if (d > maxd) maxd = d;
            }
        }
        out->mean_disp = out->n_moved > 0 ? sum / (double)out->n_moved : 0.0;
        out->max_disp = maxd;
    }
    out->sec_total = ves_clock_sec() - t0;
    if (o.verbose) {
        fprintf(stderr,
                "[recto] iter=%d supported=%zu moved=%zu reverted=%zu "
                "slope_limited=%zu mean=%.3f max=%.3f step=%.3f "
                "recto_grad=%.2f (%.2fs)\n",
                out->iterations, out->n_supported, out->n_moved,
                out->n_reverted, out->n_slope_limited, out->mean_disp,
                out->max_disp, out->mean_abs_step,
                out->mean_recto_gradient, out->sec_total);
    }
    Arena_restore(arena, mark);
    return 0;
}

/* ---------------------------------------------------------------- self-test */

#define RR_CHECK(cond, msg) \
    do { \
        if (!(cond)) { \
            fprintf(stderr, "  FAIL: %s\n", (msg)); \
            fails++; \
        } \
    } while (0)

int RectoRefine_selftest(void)
{
    enum { CHUNK = 48, SIDE = 9 };
    const char *dir = "output/_selftest_recto_refine/raw";
    Arena_T arena = Arena_new();
    uint8_t *vol = (uint8_t *)malloc(
        (size_t)CHUNK * CHUNK * CHUNK * sizeof(uint8_t));
    float *verts = (float *)malloc((size_t)SIDE * SIDE * 3 * sizeof(float));
    int32_t *faces = (int32_t *)malloc(
        (size_t)(SIDE - 1) * (SIDE - 1) * 2 * 3 * sizeof(int32_t));
    size_t nv = (size_t)SIDE * SIDE;
    size_t nf = 0, i = 0;
    int fails = 0;
    char path[512];

    fprintf(stderr, "[selftest] recto_refine\n");
    if (vol == NULL || verts == NULL || faces == NULL) {
        free(vol);
        free(verts);
        free(faces);
        Arena_dispose(&arena);
        return 1;
    }
    for (int z = 0; z < CHUNK; z++)
        for (int y = 0; y < CHUNK; y++)
            for (int x = 0; x < CHUNK; x++)
                vol[((size_t)z * CHUNK + (size_t)y) * CHUNK + (size_t)x] =
                    (uint8_t)((y >= 20 && y <= 24) ? 200 : 20);
    snprintf(path, sizeof path, "%s/z00000_y00000_x00000.tif", dir);
    ves_ensure_parent_dir(path);
    RR_CHECK(TiffIO_save(path, vol, CHUNK, CHUNK, CHUNK) == 0,
             "fixture write");
    free(vol);

    for (int z = 0; z < SIDE; z++) {
        for (int x = 0; x < SIDE; x++) {
            size_t vi = (size_t)z * SIDE + (size_t)x;
            verts[vi * 3 + 0] = (float)(12 + z);
            verts[vi * 3 + 1] = 22.0f; /* papyrus core */
            verts[vi * 3 + 2] = (float)(20 + x);
        }
    }
    for (int z = 0; z < SIDE - 1; z++) {
        for (int x = 0; x < SIDE - 1; x++) {
            int32_t a = z * SIDE + x, b = a + 1, c = a + SIDE, d = c + 1;
            faces[nf * 3 + 0] = a;
            faces[nf * 3 + 1] = b;
            faces[nf * 3 + 2] = c;
            nf++;
            faces[nf * 3 + 0] = b;
            faces[nf * 3 + 1] = d;
            faces[nf * 3 + 2] = c;
            nf++;
        }
    }

    {
        CubeTable ct;
        RectoRefineOpts o;
        RectoRefineStats st1, st2;
        double mean_y = 0.0, max_jump = 0.0;
        RectoRefineOpts_default(&o);
        o.axis_point[0] = 0.0f;
        o.axis_point[1] = 0.0f;
        o.axis_point[2] = 24.0f;
        o.axis_dir[0] = 1.0f;
        o.axis_dir[1] = o.axis_dir[2] = 0.0f;
        o.window_lo = 20.0;
        o.window_hi = 200.0;
        o.normal_blend = 0.0; /* exact +y outward in the fixture */
        o.threads = 2;
        o.verbose = 0;
        RR_CHECK(cubetable_init(&ct, arena, dir, CHUNK, verts, nv, 8.0) == 0,
                 "cubetable init");
        RR_CHECK(RectoRefine_run(arena, &ct, verts, nv, faces, nf,
                                 &o, &st1) == 0, "first run");
        for (i = 0; i < nv; i++) mean_y += verts[i * 3 + 1];
        mean_y /= (double)nv;
        for (i = 0; i < nv; i++) {
            int32_t vi = (int32_t)i;
            if ((vi % SIDE) + 1 < SIDE) {
                double d = fabs((double)verts[i * 3 + 1]
                                - verts[(i + 1) * 3 + 1]);
                if (d > max_jump) max_jump = d;
            }
            if ((vi / SIDE) + 1 < SIDE) {
                double d = fabs((double)verts[i * 3 + 1]
                                - verts[(i + SIDE) * 3 + 1]);
                if (d > max_jump) max_jump = d;
            }
        }
        RR_CHECK(st1.n_supported == nv && st1.n_moved > 0,
                 "all core vertices find the oriented recto edge");
        RR_CHECK(mean_y < 20.7 && mean_y > 19.2,
                 "core moves inward to the positive-gradient inner edge");
        RR_CHECK(max_jump < 0.1, "global displacement stays smooth");

        RR_CHECK(RectoRefine_run(arena, &ct, verts, nv, faces, nf,
                                 &o, &st2) == 0, "second run");
        RR_CHECK(st2.mean_disp < 0.20,
                 "second pass is near-idempotent at the recto edge");
    }

    free(verts);
    free(faces);
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] recto_refine %s (%d failure%s)\n",
            fails == 0 ? "PASS" : "FAIL", fails, fails == 1 ? "" : "s");
    return fails;
}
