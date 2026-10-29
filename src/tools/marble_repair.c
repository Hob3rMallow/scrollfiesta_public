/* marble_repair.c -- local, checkpointed repair of a marbled UV->XYZ patch.
 *
 * A marble is an extrinsic bending failure: the structured UV lattice still
 * exists, but its tangent frame crumples through the source volume.  Metric
 * ARAP alone cannot prefer the smooth scroll because isometric bending is free.
 * This tool crops one positively identified UV region, supplies a screened
 * low-frequency tangent-frame target to metric ARAP, strongly anchors a clean
 * halo, and atomically writes a new authoritative VMESH.  No CT snap, quilt,
 * topology change, or post-ARAP geometry step is performed. */

#include "../flatten/quad_strip.h"
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/ves_platform.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct GridInfo {
    int H, W;
    double u0, v0, du, dv;
} GridInfo;

typedef struct MaskImage {
    int W, H;
    uint8_t *pixels;
} MaskImage;

typedef struct Crop {
    int r0, r1, c0, c1;       /* inclusive full-grid vertex indices */
    int H, W;
    double commit_halo;
    float *initial;
    float *work;
    float *uv;
    float *weight;
    uint8_t *filled;
    uint8_t *fixed;
} Crop;

static void usage(void)
{
    fprintf(stderr,
        "usage: marble_repair <input.obj|input.vmesh> <output.vmesh> [options]\n"
        "options:\n"
        "  --rect U0 V0 U1 V1       inclusive detected marble rectangle (required)\n"
        "  --mask-pgm PATH           detector mask; rect restricts selected components\n"
        "  --cyl-loft Y X            replace detected data by winding-aware loft about axis\n"
        "  --loft-mask-threshold F   mask weight removed before loft (default .999)\n"
        "  --feather F               soft repair falloff from mask/rectangle (default 128 UV)\n"
        "  --halo F                  additional clean anchor halo (default 128 UV)\n"
        "  --report PATH             JSON report (default <output>.json)\n"
        "  --max-iterations N         ARAP cap (default 1000)\n"
        "  --tolerance F              max-movement stop (default 1e-5 vox)\n"
        "  --source-weight F          clean observation anchor (default 100)\n"
        "  --repair-source-scale F    anchor fraction at repair=1 (default .0001)\n"
        "  --frame-blend F            fair-frame local-step blend (default .85)\n"
        "  --frame-screen F           screened-frame data term (default .002)\n"
        "  --frame-sweeps N           fair-frame diffusion sweeps (default 128)\n"
        "  --cg-iterations N          PCG cap per coordinate (default 256)\n"
        "  --cg-tolerance F           final relative residual (default 1e-8)\n"
        "  --mg-levels N              Galerkin levels including fine (default 16)\n"
        "  --mg-cycles N              V-cycles per PCG application (default 1)\n"
        "  --max-step F               per-iteration trust radius (default 2 vox)\n"
        "  --minimum-area-ratio F     orientation barrier ratio (default .05)\n"
        "  --quiet                     suppress per-iteration progress\n"
        "  --selftest                  run quad-strip and VMESH tests\n");
}

static int pgm_token(FILE *f, char *token, size_t capacity)
{
    int ch; size_t n = 0;
    if (!f || !token || capacity < 2) return -1;
    do {
        ch = fgetc(f);
        if (ch == '#') do ch = fgetc(f); while (ch != '\n' && ch != EOF);
    } while (ch != EOF && (ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n'));
    if (ch == EOF) return -1;
    do {
        if (n + 1 >= capacity) return -1;
        token[n++] = (char)ch;
        ch = fgetc(f);
    } while (ch != EOF && ch != ' ' && ch != '\t' && ch != '\r' && ch != '\n');
    token[n] = '\0';
    return 0;
}

static int mask_read_pgm(const char *path, MaskImage *mask)
{
    FILE *f; char token[64]; long W, H, maximum; size_t n;
    memset(mask, 0, sizeof *mask);
    f = fopen(path, "rb"); if (!f) return -1;
    if (pgm_token(f, token, sizeof token) != 0 || strcmp(token, "P5") != 0 ||
        pgm_token(f, token, sizeof token) != 0 || (W = strtol(token, NULL, 10)) < 1 ||
        pgm_token(f, token, sizeof token) != 0 || (H = strtol(token, NULL, 10)) < 1 ||
        pgm_token(f, token, sizeof token) != 0 || (maximum = strtol(token, NULL, 10)) != 255 ||
        W > INT32_MAX || H > INT32_MAX || (size_t)W > SIZE_MAX / (size_t)H) {
        fclose(f); return -1;
    }
    n = (size_t)W * (size_t)H;
    mask->pixels = (uint8_t *)malloc(n);
    if (!mask->pixels || fread(mask->pixels, 1, n, f) != n || fclose(f) != 0) {
        free(mask->pixels); memset(mask, 0, sizeof *mask); return -1;
    }
    mask->W = (int)W; mask->H = (int)H;
    return 0;
}

static void mask_dispose(MaskImage *mask)
{
    if (!mask) return;
    free(mask->pixels); memset(mask, 0, sizeof *mask);
}

static double mask_sample(const MaskImage *mask, const GridInfo *g,
                          double u, double v)
{
    if (!mask || !mask->pixels) return 0.0;
    double eu = g->du * (double)(g->W - 1), ev = g->dv * (double)(g->H - 1);
    double tx = eu > 0.0 ? (u - g->u0) / eu : 0.0;
    double ty = ev > 0.0 ? (v - g->v0) / ev : 0.0;
    int x = (int)floor(tx * (double)mask->W);
    int y = (int)floor(ty * (double)mask->H);
    if (x < 0) x = 0; if (x >= mask->W) x = mask->W - 1;
    if (y < 0) y = 0; if (y >= mask->H) y = mask->H - 1;
    return (double)mask->pixels[(size_t)y * (size_t)mask->W + (size_t)x] / 255.0;
}

static int structured_grid(const MeshBinData *m, GridInfo *g)
{
    if (!m || !m->verts || !m->uv || !m->faces || m->nv < 4 || m->nf < 2)
        return -1;
    double v0 = m->uv[1];
    size_t W = 1;
    while (W < m->nv && fabs((double)m->uv[W * 2 + 1] - v0) <= 1e-5)
        W++;
    if (W < 2 || m->nv % W != 0 || W > (size_t)INT32_MAX)
        return -1;
    size_t H = m->nv / W;
    if (H < 2 || H > (size_t)INT32_MAX ||
        m->nf != 2u * (H - 1u) * (W - 1u))
        return -1;
    double u0 = m->uv[0];
    double du = (double)m->uv[2] - u0;
    double dv = (double)m->uv[W * 2 + 1] - v0;
    if (!(du > 0.0) || !(dv > 0.0) || !isfinite(du) || !isfinite(dv))
        return -1;
    for (size_t r = 0; r < H; r++) for (size_t c = 0; c < W; c++) {
        size_t i = r * W + c;
        double eu = u0 + (double)c * du, ev = v0 + (double)r * dv;
        if (fabs((double)m->uv[i * 2] - eu) > 1e-3 ||
            fabs((double)m->uv[i * 2 + 1] - ev) > 1e-3)
            return -1;
    }
    for (size_t r = 0; r + 1 < H; r++) for (size_t c = 0; c + 1 < W; c++) {
        size_t cell = r * (W - 1u) + c;
        size_t f = cell * 6u;
        int32_t a = (int32_t)(r * W + c), b = a + 1;
        int32_t q = a + (int32_t)W, d = q + 1;
        if (m->faces[f] != a || m->faces[f + 1] != b || m->faces[f + 2] != d ||
            m->faces[f + 3] != a || m->faces[f + 4] != d || m->faces[f + 5] != q)
            return -1;
    }
    g->H = (int)H; g->W = (int)W;
    g->u0 = u0; g->v0 = v0; g->du = du; g->dv = dv;
    return 0;
}

static int lower_grid_index(double value, double origin, double spacing)
{
    return (int)ceil((value - origin) / spacing - 1e-8);
}

static int upper_grid_index(double value, double origin, double spacing)
{
    return (int)floor((value - origin) / spacing + 1e-8);
}

static double repair_weight(double u, double v, const double rect[4],
                            double feather)
{
    double ou = u < rect[0] ? rect[0] - u : u > rect[2] ? u - rect[2] : 0.0;
    double ov = v < rect[1] ? rect[1] - v : v > rect[3] ? v - rect[3] : 0.0;
    double distance = hypot(ou, ov);
    if (distance <= 0.0) return 1.0;
    if (!(feather > 0.0) || distance >= feather) return 0.0;
    return 0.5 * (1.0 + cos(3.14159265358979323846 * distance / feather));
}

static void crop_dispose(Crop *c)
{
    if (!c) return;
    free(c->initial); free(c->work); free(c->uv); free(c->weight);
    free(c->filled); free(c->fixed);
    memset(c, 0, sizeof *c);
}

static int crop_build(const MeshBinData *m, const GridInfo *g,
                      const double rect[4], const MaskImage *mask,
                      double feather, double halo,
                      Crop *crop, size_t *out_active, double *out_weight_sum)
{
    memset(crop, 0, sizeof *crop);
    int box_c0 = lower_grid_index(rect[0], g->u0, g->du);
    int box_c1 = upper_grid_index(rect[2], g->u0, g->du);
    int box_r0 = lower_grid_index(rect[1], g->v0, g->dv);
    int box_r1 = upper_grid_index(rect[3], g->v0, g->dv);
    if (box_c0 < 0) box_c0 = 0; if (box_c1 >= g->W) box_c1 = g->W - 1;
    if (box_r0 < 0) box_r0 = 0; if (box_r1 >= g->H) box_r1 = g->H - 1;
    if (box_c0 > box_c1 || box_r0 > box_r1) return -1;
    int pad_c = (int)ceil((feather + halo) / g->du);
    int pad_r = (int)ceil((feather + halo) / g->dv);
    crop->c0 = box_c0 - pad_c; if (crop->c0 < 0) crop->c0 = 0;
    crop->c1 = box_c1 + pad_c; if (crop->c1 >= g->W) crop->c1 = g->W - 1;
    crop->r0 = box_r0 - pad_r; if (crop->r0 < 0) crop->r0 = 0;
    crop->r1 = box_r1 + pad_r; if (crop->r1 >= g->H) crop->r1 = g->H - 1;
    crop->W = crop->c1 - crop->c0 + 1;
    crop->H = crop->r1 - crop->r0 + 1;
    crop->commit_halo = 0.0; /* the solve itself now has an exact Dirichlet seam */
    size_t n = (size_t)crop->H * (size_t)crop->W;
    if (n > SIZE_MAX / (3 * sizeof(float)) || n > SIZE_MAX / (2 * sizeof(float)))
        return -1;
    crop->initial = (float *)malloc(n * 3 * sizeof *crop->initial);
    crop->work = (float *)malloc(n * 3 * sizeof *crop->work);
    crop->uv = (float *)malloc(n * 2 * sizeof *crop->uv);
    crop->weight = (float *)malloc(n * sizeof *crop->weight);
    crop->filled = (uint8_t *)calloc(n, 1);
    crop->fixed = (uint8_t *)calloc(n, 1);
    if (!crop->initial || !crop->work || !crop->uv || !crop->weight ||
        !crop->filled || !crop->fixed) {
        crop_dispose(crop); return -1;
    }
    size_t active = 0, clean = 0; double weight_sum = 0.0;
    for (int r = 0; r < crop->H; r++) for (int c = 0; c < crop->W; c++) {
        int gr = crop->r0 + r, gc = crop->c0 + c;
        size_t gi = (size_t)gr * (size_t)g->W + (size_t)gc;
        size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
        memcpy(&crop->initial[i * 3], &m->verts[gi * 3], 3 * sizeof(float));
        memcpy(&crop->work[i * 3], &m->verts[gi * 3], 3 * sizeof(float));
        memcpy(&crop->uv[i * 2], &m->uv[gi * 2], 2 * sizeof(float));
        double u = m->uv[gi * 2], v = m->uv[gi * 2 + 1];
        double w;
        if (mask && mask->pixels) {
            int in_rect = u >= rect[0] - 1e-8 && u <= rect[2] + 1e-8 &&
                          v >= rect[1] - 1e-8 && v <= rect[3] + 1e-8;
            w = in_rect ? mask_sample(mask, g, u, v) : 0.0;
        } else {
            w = repair_weight(u, v, rect, feather);
        }
        crop->weight[i] = (float)w;
        /* Only artificial crop boundaries are fixed. A repair touching a
         * physical sheet boundary is allowed to move that physical edge. */
        crop->fixed[i] = (uint8_t)(
            (crop->c0 > 0 && c == 0) ||
            (crop->c1 < g->W - 1 && c == crop->W - 1) ||
            (crop->r0 > 0 && r == 0) ||
            (crop->r1 < g->H - 1 && r == crop->H - 1));
    }

    /* A detector mask is a classification, not a hard cut in the variational
     * problem. Extend it by an approximate Euclidean chamfer distance and a
     * cosine falloff. Core detections stay weight one; clean tracks elsewhere
     * in the component bbox remain ordinary ARAP observations. */
    if (mask && mask->pixels && feather > 0.0) {
        double *distance = (double *)malloc(n * sizeof *distance);
        if (!distance) { crop_dispose(crop); return -1; }
        const double inf = DBL_MAX / 16.0, dd = hypot(g->du, g->dv);
        for (size_t i = 0; i < n; i++) distance[i] = crop->weight[i] > 0.0f ? 0.0 : inf;
#define RELAX_DIST(INDEX,CANDIDATE) do { double _d=(CANDIDATE); \
        if(_d<distance[(INDEX)])distance[(INDEX)]=_d; } while(0)
        for (int r = 0; r < crop->H; r++) for (int c = 0; c < crop->W; c++) {
            size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
            if (c > 0) RELAX_DIST(i, distance[i - 1] + g->du);
            if (r > 0) {
                RELAX_DIST(i, distance[i - (size_t)crop->W] + g->dv);
                if (c > 0) RELAX_DIST(i, distance[i - (size_t)crop->W - 1] + dd);
                if (c + 1 < crop->W) RELAX_DIST(i, distance[i - (size_t)crop->W + 1] + dd);
            }
        }
        for (int r = crop->H - 1; r >= 0; r--) for (int c = crop->W - 1; c >= 0; c--) {
            size_t i = (size_t)r * (size_t)crop->W + (size_t)c;
            if (c + 1 < crop->W) RELAX_DIST(i, distance[i + 1] + g->du);
            if (r + 1 < crop->H) {
                RELAX_DIST(i, distance[i + (size_t)crop->W] + g->dv);
                if (c > 0) RELAX_DIST(i, distance[i + (size_t)crop->W - 1] + dd);
                if (c + 1 < crop->W) RELAX_DIST(i, distance[i + (size_t)crop->W + 1] + dd);
            }
        }
#undef RELAX_DIST
        for (size_t i = 0; i < n; i++) if (crop->weight[i] <= 0.0f && distance[i] < feather)
            crop->weight[i] = (float)(0.5 * (1.0 + cos(3.14159265358979323846 *
                                                      distance[i] / feather)));
        free(distance);
    }
    for (size_t i = 0; i < n; i++) {
        double w = crop->weight[i];
        if (w > 0.0) { active++; weight_sum += w; } else clean++;
    }
    if (active == 0 || clean == 0) { crop_dispose(crop); return -1; }
    if (out_active) *out_active = active;
    if (out_weight_sum) *out_weight_sum = weight_sum;
    return 0;
}

/* Treat the positively detected core as absent and reconstruct it using the
 * same winding-aware initializer as the main ribbon pipeline. Cartesian
 * Laplacian fill cuts chords through a roll; (z,r,unwrapped theta) continues a
 * one-sided sheet tail and preserves its winding gauge. The soft detector
 * feather blends that loft into untouched observations. This is an ARAP
 * initializer/anchor, never a post-ARAP geometry operation. */
static int crop_cylindrical_loft(Crop *crop, const GridInfo *g,
                                 double axis_y, double axis_x,
                                 double mask_threshold,
                                 size_t *out_removed, double *out_rms,
                                 double *out_max)
{
    if (!crop || !g || !isfinite(axis_y) || !isfinite(axis_x) ||
        !(mask_threshold > 0.0 && mask_threshold <= 1.0)) return -1;
    size_t n = (size_t)crop->H * (size_t)crop->W, removed = 0, clean = 0;
    uint8_t *domain = (uint8_t *)malloc(n);
    if (!domain) return -1;
    for (size_t i = 0; i < n; i++) {
        domain[i] = (uint8_t)((double)crop->weight[i] < mask_threshold);
        if (domain[i]) clean++; else removed++;
    }
    if (removed == 0 || clean == 0) { free(domain); return -1; }
    Arena_T arena = Arena_new();
    if (!arena) { free(domain); return -1; }
    QuadStripOpts q; QuadStrip_defaults(&q);
    q.mode = QUAD_STRIP_RECT; q.cylindrical_fill = 1;
    q.axis_y = axis_y; q.axis_x = axis_x;
    q.pde_max_iter = 4096; q.pde_tol = 1e-7; q.exact_fitted = 1;
    float *loft = NULL; size_t nv = 0;
    int rc = QuadStrip_build(arena, domain, crop->initial, crop->H, crop->W,
                             crop->r0, crop->c0, g->du, &q,
                             &loft, &nv, NULL, NULL, NULL, NULL);
    free(domain);
    if (rc != 0 || !loft || nv != n) { Arena_dispose(&arena); return -1; }
    double ss = 0.0, mx = 0.0;
    for (size_t i = 0; i < n; i++) {
        double w = fmin(fmax((double)crop->weight[i], 0.0), 1.0), m2 = 0.0;
        for (int ch = 0; ch < 3; ch++) {
            double before = crop->initial[i*3+(size_t)ch];
            double after = (1.0-w)*before + w*(double)loft[i*3+(size_t)ch];
            crop->work[i*3+(size_t)ch] = (float)after;
            double d = after-before; m2 += d*d;
        }
        ss += m2; if (sqrt(m2) > mx) mx = sqrt(m2);
    }
    Arena_dispose(&arena);
    if (out_removed) *out_removed = removed;
    if (out_rms) *out_rms = sqrt(ss/(double)n);
    if (out_max) *out_max = mx;
    return 0;
}

static double commit_factor(const GridInfo *g, const Crop *crop, int gr, int gc)
{
    if (!(crop->commit_halo > 0.0)) return 1.0;
    double f = 1.0;
#define MARBLE_TAPER(DISTANCE) do { \
        double _t=(DISTANCE)/crop->commit_halo; \
        if(_t<0.0)_t=0.0;else if(_t>1.0)_t=1.0; \
        _t=_t*_t*(3.0-2.0*_t);if(_t<f)f=_t; \
    } while(0)
    if (crop->c0 > 0) MARBLE_TAPER((double)(gc - crop->c0) * g->du);
    if (crop->c1 < g->W - 1) MARBLE_TAPER((double)(crop->c1 - gc) * g->du);
    if (crop->r0 > 0) MARBLE_TAPER((double)(gr - crop->r0) * g->dv);
    if (crop->r1 < g->H - 1) MARBLE_TAPER((double)(crop->r1 - gr) * g->dv);
#undef MARBLE_TAPER
    return f;
}

static void trial_position(const MeshBinData *m, const GridInfo *g,
                           const Crop *crop, int gr, int gc, double alpha,
                           double out[3])
{
    size_t gi = (size_t)gr * (size_t)g->W + (size_t)gc;
    for (int ch = 0; ch < 3; ch++) out[ch] = m->verts[gi * 3 + (size_t)ch];
    if (gr < crop->r0 || gr > crop->r1 || gc < crop->c0 || gc > crop->c1)
        return;
    size_t i = (size_t)(gr - crop->r0) * (size_t)crop->W +
               (size_t)(gc - crop->c0);
    double factor = commit_factor(g, crop, gr, gc) * alpha;
    for (int ch = 0; ch < 3; ch++)
        out[ch] += factor * ((double)crop->work[i * 3 + (size_t)ch] -
                            (double)crop->initial[i * 3 + (size_t)ch]);
}

static int commit_is_valid(const MeshBinData *m, const GridInfo *g,
                           const Crop *crop, double alpha, double minimum_ratio)
{
    int r0 = crop->r0 > 0 ? crop->r0 - 1 : 0;
    int c0 = crop->c0 > 0 ? crop->c0 - 1 : 0;
    int r1 = crop->r1 < g->H - 1 ? crop->r1 : g->H - 2;
    int c1 = crop->c1 < g->W - 1 ? crop->c1 : g->W - 2;
    for (int r = r0; r <= r1; r++) for (int c = c0; c <= c1; c++) {
        int vr[4] = {r, r, r + 1, r + 1};
        int vc[4] = {c, c + 1, c, c + 1};
        int tri[2][3] = {{0, 1, 3}, {0, 3, 2}};
        for (int t = 0; t < 2; t++) {
            double oldp[3][3], newp[3][3];
            for (int q = 0; q < 3; q++) {
                int v = tri[t][q];
                size_t gi = (size_t)vr[v] * (size_t)g->W + (size_t)vc[v];
                for (int ch = 0; ch < 3; ch++)
                    oldp[q][ch] = m->verts[gi * 3 + (size_t)ch];
                trial_position(m, g, crop, vr[v], vc[v], alpha, newp[q]);
            }
            double e10[3], e20[3], e11[3], e21[3], n0[3], n1[3];
            for (int ch = 0; ch < 3; ch++) {
                e10[ch] = oldp[1][ch] - oldp[0][ch];
                e20[ch] = oldp[2][ch] - oldp[0][ch];
                e11[ch] = newp[1][ch] - newp[0][ch];
                e21[ch] = newp[2][ch] - newp[0][ch];
            }
            n0[0] = e10[1]*e20[2] - e10[2]*e20[1];
            n0[1] = e10[2]*e20[0] - e10[0]*e20[2];
            n0[2] = e10[0]*e20[1] - e10[1]*e20[0];
            n1[0] = e11[1]*e21[2] - e11[2]*e21[1];
            n1[1] = e11[2]*e21[0] - e11[0]*e21[2];
            n1[2] = e11[0]*e21[1] - e11[1]*e21[0];
            double a0 = n0[0]*n0[0] + n0[1]*n0[1] + n0[2]*n0[2];
            double a1 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
            double dot = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
            int tapered = 0;
            for (int q = 0; q < 3; q++) {
                int v = tri[t][q];
                if (vr[v] < crop->r0 || vr[v] > crop->r1 ||
                    vc[v] < crop->c0 || vc[v] > crop->c1 ||
                    commit_factor(g, crop, vr[v], vc[v]) < 1.0 - 1e-9) {
                    tapered = 1; break;
                }
            }
            /* The ARAP barrier certified every interior update along its
             * nonlinear path.  Its endpoint need not lie in the same convex
             * normal hemisphere as the original after a large smooth bend.
             * Only newly coupled crop/outside faces require the direct
             * original-normal transaction check. */
            if (!isfinite(a1) || a1 < 1e-24 ||
                (tapered && a0 >= 1e-24 && dot < minimum_ratio*a0))
                return 0;
        }
    }
    return 1;
}

static void json_string(FILE *f, const char *s)
{
    fputc('"', f);
    for (; *s; s++) {
        unsigned char c = (unsigned char)*s;
        if (c == '"' || c == '\\') { fputc('\\', f); fputc(c, f); }
        else if (c == '\n') fputs("\\n", f);
        else if (c == '\r') fputs("\\r", f);
        else if (c == '\t') fputs("\\t", f);
        else if (c < 32) fprintf(f, "\\u%04x", (unsigned)c);
        else fputc(c, f);
    }
    fputc('"', f);
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return QuadStrip_selftest() == 0 && MeshBin_selftest() == 0 ? 0 : 1;
    if (argc < 4) { usage(); return 1; }
    const char *input = argv[1], *output = argv[2], *report_path = NULL;
    const char *mask_path = NULL;
    char default_report[2600], resolved[2600];
    double rect[4] = {0,0,0,0}, feather = 128.0, halo = 128.0;
    double loft_axis_y = 0.0, loft_axis_x = 0.0, loft_mask_threshold = 0.999;
    int have_rect = 0, quiet = 0, have_cyl_loft = 0;
    QuadStripArapOpts opts; QuadStripArap_defaults(&opts);
    opts.source_weight = 100.0;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--rect") == 0 && i + 4 < argc) {
            for (int q = 0; q < 4; q++) rect[q] = strtod(argv[++i], NULL);
            have_rect = 1;
        } else if (strcmp(argv[i], "--mask-pgm") == 0 && i + 1 < argc)
            mask_path = argv[++i];
        else if (strcmp(argv[i], "--cyl-loft") == 0 && i + 2 < argc) {
            loft_axis_y = strtod(argv[++i], NULL);
            loft_axis_x = strtod(argv[++i], NULL);
            have_cyl_loft = 1;
        } else if (strcmp(argv[i], "--loft-mask-threshold") == 0 && i + 1 < argc)
            loft_mask_threshold = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--feather") == 0 && i + 1 < argc)
            feather = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--halo") == 0 && i + 1 < argc)
            halo = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc)
            report_path = argv[++i];
        else if (strcmp(argv[i], "--max-iterations") == 0 && i + 1 < argc)
            opts.max_iterations = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--tolerance") == 0 && i + 1 < argc)
            opts.movement_tolerance = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--source-weight") == 0 && i + 1 < argc)
            opts.source_weight = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--repair-source-scale") == 0 && i + 1 < argc)
            opts.repair_source_scale = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--frame-blend") == 0 && i + 1 < argc)
            opts.repair_frame_blend = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--frame-screen") == 0 && i + 1 < argc)
            opts.repair_frame_screen = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--frame-sweeps") == 0 && i + 1 < argc)
            opts.repair_frame_sweeps = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--cg-iterations") == 0 && i + 1 < argc)
            opts.cg_max_iterations = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--cg-tolerance") == 0 && i + 1 < argc)
            opts.cg_relative_tolerance = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--mg-levels") == 0 && i + 1 < argc)
            opts.multigrid_levels = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--mg-cycles") == 0 && i + 1 < argc)
            opts.multigrid_cycles = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--max-step") == 0 && i + 1 < argc)
            opts.maximum_vertex_step = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--minimum-area-ratio") == 0 && i + 1 < argc)
            opts.minimum_area_ratio = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--quiet") == 0) quiet = 1;
        else { fprintf(stderr, "unknown or incomplete option: %s\n", argv[i]); usage(); return 1; }
    }
    if (!report_path) {
        if (snprintf(default_report, sizeof default_report, "%s.json", output) >= (int)sizeof default_report)
            return 1;
        report_path = default_report;
    }
    if (!have_rect || !isfinite(rect[0]) || !isfinite(rect[1]) ||
        !isfinite(rect[2]) || !isfinite(rect[3]) || rect[0] > rect[2] ||
        rect[1] > rect[3] || feather < 0.0 || halo < 0.0 ||
        (have_cyl_loft && (!mask_path || !isfinite(loft_axis_y) ||
                           !isfinite(loft_axis_x) || !(loft_mask_threshold > 0.0) ||
                           loft_mask_threshold > 1.0)) ||
        strcmp(input, output) == 0) {
        usage(); return 1;
    }
    opts.verbose = quiet ? 0 : 1;

    MeshBinData mesh; memset(&mesh, 0, sizeof mesh);
    if (MeshBin_companion_path(input, resolved, sizeof resolved) != 0 ||
        MeshBin_read_malloc(resolved, &mesh) != 0) {
        fprintf(stderr, "ERROR: cannot read authoritative input for %s\n", input);
        return 1;
    }
    GridInfo grid;
    if (structured_grid(&mesh, &grid) != 0) {
        fprintf(stderr, "ERROR: input is not one exact row-major RECT ribbon\n");
        MeshBin_dispose(&mesh); return 1;
    }
    MaskImage mask; memset(&mask, 0, sizeof mask);
    if (mask_path && mask_read_pgm(mask_path, &mask) != 0) {
        fprintf(stderr, "ERROR: cannot read P5/255 detector mask %s\n", mask_path);
        MeshBin_dispose(&mesh); return 1;
    }
    Crop crop; size_t selected = 0; double weight_sum = 0.0;
    if (crop_build(&mesh, &grid, rect, mask_path ? &mask : NULL,
                   feather, halo, &crop,
                   &selected, &weight_sum) != 0) {
        fprintf(stderr, "ERROR: repair rectangle has no selected cells or no clean anchor halo\n");
        mask_dispose(&mask); MeshBin_dispose(&mesh); return 1;
    }
    mask_dispose(&mask);
    size_t loft_removed = 0; double loft_rms = 0.0, loft_max = 0.0;
    if (have_cyl_loft && crop_cylindrical_loft(&crop, &grid,
            loft_axis_y, loft_axis_x, loft_mask_threshold,
            &loft_removed, &loft_rms, &loft_max) != 0) {
        fprintf(stderr, "ERROR: winding-aware marble loft failed\n");
        crop_dispose(&crop); MeshBin_dispose(&mesh); return 1;
    }
    fprintf(stderr,
        "[marble_repair] input=%s grid=%dx%d du/dv=%.6g/%.6g\n"
        "  rect=[%.3f %.3f %.3f %.3f], crop rows %d..%d cols %d..%d (%dx%d), "
        "selected=%zu weight-sum=%.3f mask=%s\n",
        resolved, grid.H, grid.W, grid.du, grid.dv,
        rect[0], rect[1], rect[2], rect[3], crop.r0, crop.r1, crop.c0, crop.c1,
        crop.H, crop.W, selected, weight_sum, mask_path ? mask_path : "<rectangle>");
    if (have_cyl_loft) fprintf(stderr,
        "  cylindrical loft initializer: axis=(%.3f,%.3f), removed=%zu, "
        "blend rms/max=%.6g/%.6g\n",
        loft_axis_y, loft_axis_x, loft_removed, loft_rms, loft_max);

    opts.repair_weight = crop.weight;
    opts.fixed_vertices = crop.fixed;
    opts.nonlinear_multigrid_levels = 1; /* the crop is already the spatial multigrid unit */
    QuadStripArapStats stats; double t0 = ves_clock_sec();
    int rc = QuadStrip_metric_arap(crop.work, crop.uv, crop.filled,
                                   crop.H, crop.W, &opts, &stats);
    double solve_seconds = ves_clock_sec() - t0;
    if (rc != 0) {
        fprintf(stderr, "ERROR: marble-aware metric ARAP failed, code %d; no output written\n", rc);
        crop_dispose(&crop); MeshBin_dispose(&mesh); return 1;
    }

    double commit_alpha = 1.0; int commit_halvings = 0;
    while (commit_halvings < 32 &&
           !commit_is_valid(&mesh, &grid, &crop, commit_alpha, opts.minimum_area_ratio)) {
        commit_alpha *= 0.5; commit_halvings++;
    }
    if (!commit_is_valid(&mesh, &grid, &crop, commit_alpha, opts.minimum_area_ratio)) {
        fprintf(stderr, "ERROR: no orientation-safe feathered commit; no output written\n");
        crop_dispose(&crop); MeshBin_dispose(&mesh); return 1;
    }
    double sumsq = 0.0, maxmove = 0.0; size_t moved = 0;
    for (int r = 0; r < crop.H; r++) for (int c = 0; c < crop.W; c++) {
        size_t i = (size_t)r * (size_t)crop.W + (size_t)c;
        size_t gi = (size_t)(crop.r0 + r) * (size_t)grid.W + (size_t)(crop.c0 + c);
        double m2 = 0.0;
        for (int ch = 0; ch < 3; ch++) {
            double d = commit_alpha * commit_factor(&grid, &crop,
                       crop.r0 + r, crop.c0 + c) *
                       ((double)crop.work[i * 3 + (size_t)ch] -
                        (double)crop.initial[i * 3 + (size_t)ch]);
            mesh.verts[gi * 3 + (size_t)ch] = (float)((double)mesh.verts[gi * 3 + (size_t)ch] + d);
            m2 += d*d;
        }
        double movement = sqrt(m2); sumsq += m2; if (movement > maxmove) maxmove = movement;
        if (movement > 0.0) moved++;
    }
    if (ves_ensure_parent_dir(output) != 0 || ves_ensure_parent_dir(report_path) != 0 ||
        MeshBin_write(output, mesh.verts, mesh.nv, mesh.faces, mesh.nf, mesh.uv) != 0) {
        fprintf(stderr, "ERROR: cannot atomically checkpoint %s\n", output);
        crop_dispose(&crop); MeshBin_dispose(&mesh); return 1;
    }
    FILE *report = fopen(report_path, "wb");
    if (!report) {
        fprintf(stderr, "ERROR: geometry checkpoint exists but report cannot be written: %s\n", report_path);
        crop_dispose(&crop); MeshBin_dispose(&mesh); return 1;
    }
    fputs("{\n \"schema\":\"marble-frame-arap-repair-v1\",\n \"input\":", report);
    json_string(report, resolved); fputs(",\n \"output\":", report); json_string(report, output);
    fputs(",\n \"repair_mask\":", report);
    if (mask_path) json_string(report, mask_path); else fputs("null", report);
    fprintf(report,
        ",\n \"geometry_terminal\":\"marble-aware metric ARAP; no post-ARAP geometry\",\n"
        " \"initializer\":{\"type\":\"%s\",\"axis_y\":%.17g,\"axis_x\":%.17g,"
        "\"mask_threshold\":%.17g,\"removed_vertices\":%zu,"
        "\"blend_rms_displacement\":%.17g,\"blend_max_displacement\":%.17g},\n"
        " \"grid\":{\"height\":%d,\"width\":%d,\"du\":%.17g,\"dv\":%.17g},\n"
        " \"repair_rectangle_uv\":[%.17g,%.17g,%.17g,%.17g],\n"
        " \"crop\":{\"r0\":%d,\"r1\":%d,\"c0\":%d,\"c1\":%d,"
        "\"height\":%d,\"width\":%d},\n"
        " \"selection\":{\"vertices\":%zu,\"weight_sum\":%.17g,"
        "\"feather_uv\":%.17g,\"clean_halo_uv\":%.17g},\n"
        " \"solver\":{\"iterations\":%d,\"converged\":%s,"
        "\"movement_tolerance\":%.17g,\"last_max_movement\":%.17g,"
        "\"last_rms_movement\":%.17g,\"barrier_rejections\":%d,"
        "\"barrier_zero_fallbacks\":%d,\"barrier_stalled\":%s,"
        "\"minimum_step_scale\":%.17g,\"edge_log_rms_before\":%.17g,"
        "\"edge_log_rms_after\":%.17g,\"edge_log_max_before\":%.17g,"
        "\"edge_log_max_after\":%.17g,\"repair_frame_sweeps\":%d,"
        "\"repair_frame_blend\":%.17g,\"repair_frame_screen\":%.17g,"
        "\"source_weight\":%.17g,\"repair_source_scale\":%.17g,"
        "\"seconds\":%.17g},\n"
        " \"commit\":{\"orientation_safe\":true,\"alpha\":%.17g,"
        "\"halvings\":%d,\"moved_vertices\":%zu,\"rms_displacement\":%.17g,"
        "\"max_displacement\":%.17g},\n"
        " \"topology_changed\":false,\n \"uv_changed\":false\n}\n",
        have_cyl_loft ? "cylindrical_laplacian_loft" : "input_geometry",
        loft_axis_y, loft_axis_x, loft_mask_threshold, loft_removed, loft_rms, loft_max,
        grid.H, grid.W, grid.du, grid.dv, rect[0], rect[1], rect[2], rect[3],
        crop.r0, crop.r1, crop.c0, crop.c1, crop.H, crop.W,
        selected, weight_sum, feather, halo,
        stats.iterations, stats.converged ? "true" : "false",
        opts.movement_tolerance, stats.maximum_vertex_movement,
        stats.rms_vertex_movement, stats.barrier_rejections,
        stats.barrier_zero_fallbacks, stats.barrier_stalled ? "true" : "false",
        stats.minimum_step_scale, stats.initial_edge_log_rms, stats.final_edge_log_rms,
        stats.initial_edge_log_max, stats.final_edge_log_max,
        opts.repair_frame_sweeps, opts.repair_frame_blend, opts.repair_frame_screen,
        opts.source_weight, opts.repair_source_scale, solve_seconds,
        commit_alpha, commit_halvings, moved,
        moved ? sqrt(sumsq/(double)moved) : 0.0, maxmove);
    fclose(report);
    fprintf(stderr,
        "[marble_repair] checkpointed %s; iterations=%d converged=%s "
        "edge-log-rms %.6g->%.6g commit alpha=%.6g moved=%zu rms/max=%.6g/%.6g (%.1fs)\n",
        output, stats.iterations, stats.converged ? "yes" : "no",
        stats.initial_edge_log_rms, stats.final_edge_log_rms,
        commit_alpha, moved, moved ? sqrt(sumsq/(double)moved) : 0.0,
        maxmove, solve_seconds);
    crop_dispose(&crop); MeshBin_dispose(&mesh);
    return 0;
}
