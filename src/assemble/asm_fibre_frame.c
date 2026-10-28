/* asm_fibre_frame.c -- see asm_fibre_frame.h. */
#include "asm_fibre_frame.h"
#include "../common/zarr_u8.h"
#include "../common/ves_platform.h"
#include "../common/pipeline_constants.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AFF_PI 3.14159265358979323846

typedef struct { uint64_t key; size_t chart; } AffOrder;

static int aff_order(const void *a, const void *b)
{
    const AffOrder *x = a, *y = b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    return x->chart < y->chart ? -1 : x->chart > y->chart;
}

/* The CT along the normal at one surface point: nearest voxel at -1, 0, +1 vox, averaged; NAN when the centre reads
 * no datum (0: outside the volume, a missing chunk or the mask of a masked scan). */
static float aff_ct(ZarrU8 *z, const double p[3], const double n[3])
{
    double s = 0;
    int k = 0;
    if (ZarrU8_sample(z, p[0], p[1], p[2]) == 0) return NAN;
    for (int t = -1; t <= 1; t++) {
        uint8_t v = ZarrU8_sample(z, p[0] + t*n[0], p[1] + t*n[1], p[2] + t*n[2]);
        if (v) { s += v; k++; }
    }
    return k ? (float)(s / k) : NAN;
}

/* Sampled Gaussian derivative kernels (correlation form) at sigma AFF_SIGMA, truncated at 4 sigma: g sums to 1, g1 and
 * g2 are the analytic first and second derivatives of that same g (g2 made zero-mean).  What matters for the angle is
 * that g1 x g1 and g2 x g agree at the fibres' frequencies: second differences on a small stencil read a 20 deg cross
 * as 15 deg, and kernels truncated at 3 sigma and renormalized to polynomial moments still read 18.3; these read every
 * angle within 0.01 deg (period 12 vox, 0-45 deg).  Sigma 4 peaks near an 18 vox period, the scale of the papyrus
 * fibre bundles in the 9 um scans. */
#define AFF_SIGMA 4.0
#define AFF_RADIUS 16
static void aff_kernels(double g[2*AFF_RADIUS+1], double g1[2*AFF_RADIUS+1], double g2[2*AFF_RADIUS+1])
{
    double s = 0, z2 = 0, s2 = AFF_SIGMA*AFF_SIGMA;
    for (int k = -AFF_RADIUS; k <= AFF_RADIUS; k++) { g[k+AFF_RADIUS] = exp(-(double)(k*k)/(2*s2)); s += g[k+AFF_RADIUS]; }
    for (int k = -AFF_RADIUS; k <= AFF_RADIUS; k++) {
        double v = g[k+AFF_RADIUS] /= s;
        g1[k+AFF_RADIUS] = (double)k/s2*v;
        g2[k+AFF_RADIUS] = ((double)(k*k)/(s2*s2) - 1/s2)*v; z2 += g2[k+AFF_RADIUS];
    }
    for (int k = -AFF_RADIUS; k <= AFF_RADIUS; k++) g2[k+AFF_RADIUS] -= z2*g[k+AFF_RADIUS];
}

/* One separable pass: out_q[n] = sum_k ker_q[k] in[n + k*stride] (q < nq <= 3) along a line of len nodes, where the
 * whole window is known (NAN elsewhere).  bad: a prefix count of unknown nodes along the line (len + 1 entries). */
static void aff_pass(const float *in, size_t stride, size_t len, int nq, const double *const *ker, float *const *out, size_t *bad)
{
    bad[0] = 0;
    for (size_t n = 0; n < len; n++) bad[n+1] = bad[n] + !isfinite(in[n*stride]);
    for (size_t n = 0; n < len; n++) {
        for (int q = 0; q < nq; q++) out[q][n*stride] = NAN;
        if (n < AFF_RADIUS || n + AFF_RADIUS >= len || bad[n+AFF_RADIUS+1] != bad[n-AFF_RADIUS]) continue;
        double acc[3] = { 0, 0, 0 };
        for (int k = -AFF_RADIUS; k <= AFF_RADIUS; k++) {
            double v = in[(size_t)((long)n + k)*stride];
            for (int q = 0; q < nq; q++) acc[q] += ker[q][k+AFF_RADIUS]*v;
        }
        for (int q = 0; q < nq; q++) out[q][n*stride] = (float)acc[q];
    }
}

/* One chart in its current UV frame: raster at 1 vox, CT, two smoothing passes, gradient moments. */
static void aff_chart(ZarrU8 *z, Arena_T arena, const AsmChart *c, AsmFibreMoments *m,
                       size_t *face_samples, int physical_patch)
{
    Arena_Mark mark = Arena_save(arena);
    double u0 = INFINITY, v0 = INFINITY, u1 = -INFINITY, v1 = -INFINITY;
    for (size_t v = 0; v < c->nv; v++) {
        double p[2];
        AsmChart_point(c, v, p);
        u0 = fmin(u0, p[0]); u1 = fmax(u1, p[0]); v0 = fmin(v0, p[1]); v1 = fmax(v1, p[1]);
    }
    if (!(u1 - u0 > 2*AFF_RADIUS) || !(v1 - v0 > 2*AFF_RADIUS) || !isfinite(u1 - u0) || !isfinite(v1 - v0)) { Arena_restore(arena, mark); return; }
    size_t W = (size_t)ceil(u1 - u0) + 1, H = (size_t)ceil(v1 - v0) + 1, N = W*H;
    float *img = ARENA_ALLOC(arena, N*sizeof *img);
    uint8_t *mapped = physical_patch ? ARENA_CALLOC(arena,N,sizeof *mapped) : NULL;
    double *world = physical_patch ? ARENA_ALLOC(arena,3*N*sizeof *world) : NULL;
    for (size_t i = 0; i < N; i++) img[i] = NAN;
    for (size_t t = 0; t < c->nf; t++) {
        const int32_t *f = c->faces + 3*t;
        double q[3][2];
        for (int a = 0; a < 3; a++) AsmChart_point(c, (size_t)f[a], q[a]);
        double det = (q[1][0]-q[0][0])*(q[2][1]-q[0][1]) - (q[2][0]-q[0][0])*(q[1][1]-q[0][1]);
        if (!(fabs(det) > 1e-12)) continue;
        double fn[3] = { 0, 0, 0 };                                 /* the face normal, when the chart has none */
        if (!c->nrm) {
            double e[3], g[3];
            for (int d = 0; d < 3; d++) {
                e[d] = (double)c->xyz[3*(size_t)f[1]+(size_t)d] - c->xyz[3*(size_t)f[0]+(size_t)d];
                g[d] = (double)c->xyz[3*(size_t)f[2]+(size_t)d] - c->xyz[3*(size_t)f[0]+(size_t)d];
            }
            fn[0] = e[1]*g[2]-e[2]*g[1]; fn[1] = e[2]*g[0]-e[0]*g[2]; fn[2] = e[0]*g[1]-e[1]*g[0];
        }
        long i0 = (long)ceil(fmin(q[0][0], fmin(q[1][0], q[2][0])) - u0), i1 = (long)floor(fmax(q[0][0], fmax(q[1][0], q[2][0])) - u0);
        long j0 = (long)ceil(fmin(q[0][1], fmin(q[1][1], q[2][1])) - v0), j1 = (long)floor(fmax(q[0][1], fmax(q[1][1], q[2][1])) - v0);
        for (long j = j0 < 0 ? 0 : j0; j <= j1 && j < (long)H; j++) for (long i = i0 < 0 ? 0 : i0; i <= i1 && i < (long)W; i++) {
            size_t node = (size_t)j*W + (size_t)i;
            if ((!physical_patch && isfinite(img[node])) || (mapped && mapped[node] == 2)) continue;
            double x = u0 + (double)i, y = v0 + (double)j;
            double l1 = ((x-q[0][0])*(q[2][1]-q[0][1]) - (q[2][0]-q[0][0])*(y-q[0][1]))/det;
            double l2 = ((q[1][0]-q[0][0])*(y-q[0][1]) - (x-q[0][0])*(q[1][1]-q[0][1]))/det;
            double l0 = 1 - l1 - l2;
            if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
            double p[3], n[3], len = 0;
            for (int d = 0; d < 3; d++) {
                p[d] = l0*c->xyz[3*(size_t)f[0]+(size_t)d] + l1*c->xyz[3*(size_t)f[1]+(size_t)d] + l2*c->xyz[3*(size_t)f[2]+(size_t)d];
                n[d] = c->nrm ? l0*c->nrm[3*(size_t)f[0]+(size_t)d] + l1*c->nrm[3*(size_t)f[1]+(size_t)d] + l2*c->nrm[3*(size_t)f[2]+(size_t)d] : fn[d];
                len += n[d]*n[d];
            }
            if (mapped) {
                if (mapped[node]) {
                    double error = 0, scale = 1;
                    for (int d = 0; d < 3; d++) {
                        error = fmax(error,fabs(p[d]-world[3*node+(size_t)d]));
                        scale += fabs(p[d])+fabs(world[3*node+(size_t)d]);
                    }
                    if (error > 64*DBL_EPSILON*scale) { mapped[node] = 2; img[node] = NAN; }
                    continue;
                }
                mapped[node] = 1;
                for (int d = 0; d < 3; d++) world[3*node+(size_t)d] = p[d];
            }
            len = sqrt(len);
            for (int d = 0; d < 3; d++) n[d] = len > 1e-9 ? n[d]/len : 0;
            img[node] = aff_ct(z, p, n);
        }
    }
    /* Gaussian derivatives (separable: along u with g, g1, g2, then along v), then the oriented second-derivative
     * energy: D(theta) = P + Q cos 2theta + R sin 2theta with P = (Suu + Svv)/2, Q = (Suu - Svv)/2, R = Suv, and
     * D^2 = w + Re(m2 e^{-i2theta}) + Re(m4 e^{-i4theta}) for w = P^2 + (Q^2 + R^2)/2, m2 = 2P(Q + iR), m4 = (Q + iR)^2/2.
     * Each family adds to D linearly with cos^2 / sin^2 weights, so crossing families reinforce the 4-fold term. */
    double g[2*AFF_RADIUS+1], g1[2*AFF_RADIUS+1], g2[2*AFF_RADIUS+1];
    aff_kernels(g, g1, g2);
    float *a0 = ARENA_ALLOC(arena, N*sizeof *a0), *a1 = ARENA_ALLOC(arena, N*sizeof *a1), *a2 = ARENA_ALLOC(arena, N*sizeof *a2);
    float *suu = ARENA_ALLOC(arena, N*sizeof *suu), *svv = ARENA_ALLOC(arena, N*sizeof *svv), *suv = ARENA_ALLOC(arena, N*sizeof *suv);
    size_t *bad = ARENA_ALLOC(arena, ((W > H ? W : H) + 1)*sizeof *bad);
    const double *rows[3] = { g, g1, g2 }, *smooth[1] = { g }, *first[1] = { g1 }, *second[1] = { g2 };
    for (size_t j = 0; j < H; j++) {                  /* along u: smoothed, d/du, d2/du2 */
        float *out[3] = { a0 + j*W, a1 + j*W, a2 + j*W };
        aff_pass(img + j*W, 1, W, 3, rows, out, bad);
    }
    for (size_t i = 0; i < W; i++) {                  /* along v: Suu = g(v) a2, Svv = g2(v) a0, Suv = g1(v) a1 */
        float *ouu[1] = { suu + i }, *ovv[1] = { svv + i }, *ouv[1] = { suv + i };
        aff_pass(a2 + i, W, H, 1, smooth, ouu, bad);
        aff_pass(a0 + i, W, H, 1, second, ovv, bad);
        aff_pass(a1 + i, W, H, 1, first, ouv, bad);
    }
    for (size_t n = 0; n < N; n++) {
        if (!isfinite(suu[n]) || !isfinite(svv[n]) || !isfinite(suv[n])) continue;
        double P = .5*((double)suu[n] + svv[n]), Q = .5*((double)suu[n] - svv[n]), R = suv[n];
        m->w += P*P + .5*(Q*Q + R*R);
        m->m2[0] += 2*P*Q; m->m2[1] += 2*P*R;
        m->m4[0] += .5*(Q*Q - R*R); m->m4[1] += Q*R;
        m->px++;
    }
    if (face_samples) {
        /* Include all faces sharing a valid edge/vertex sample. Assigning it
         * to only the first triangle loses half the support on a regular
         * two-voxel triangular grid. Ambiguous/missing samples have already
         * invalidated their full Gaussian derivative stencil. */
        for (size_t t = 0; t < c->nf; t++) {
            const int32_t *f = c->faces+3*t; double q[3][2];
            for (int a = 0; a < 3; a++) AsmChart_point(c,(size_t)f[a],q[a]);
            double det = (q[1][0]-q[0][0])*(q[2][1]-q[0][1])-(q[2][0]-q[0][0])*(q[1][1]-q[0][1]);
            if (!(fabs(det)>1e-12)) continue;
            long i0=(long)ceil(fmin(q[0][0],fmin(q[1][0],q[2][0]))-u0),i1=(long)floor(fmax(q[0][0],fmax(q[1][0],q[2][0]))-u0);
            long j0=(long)ceil(fmin(q[0][1],fmin(q[1][1],q[2][1]))-v0),j1=(long)floor(fmax(q[0][1],fmax(q[1][1],q[2][1]))-v0);
            for (long j=j0<0?0:j0;j<=j1 && j<(long)H;j++) for (long i=i0<0?0:i0;i<=i1 && i<(long)W;i++) {
                size_t node=(size_t)j*W+(size_t)i;
                if (!isfinite(suu[node]) || !isfinite(svv[node]) || !isfinite(suv[node])) continue;
                double x=u0+(double)i,y=v0+(double)j;
                double l1=((x-q[0][0])*(q[2][1]-q[0][1])-(q[2][0]-q[0][0])*(y-q[0][1]))/det;
                double l2=((q[1][0]-q[0][0])*(y-q[0][1])-(x-q[0][0])*(q[1][1]-q[0][1]))/det;
                if (l1>=-1e-9 && l2>=-1e-9 && 1-l1-l2>=-1e-9) face_samples[t]++;
            }
        }
    }
    Arena_restore(arena, mark);
}

void AsmFibre_measure_patch(ZarrU8 *volume, Arena_T arena, const AsmChart *chart,
                            AsmFibreMoments *out, size_t *face_samples)
{
    if (!out) return;
    memset(out, 0, sizeof *out);
    if (face_samples && chart) memset(face_samples,0,chart->nf*sizeof *face_samples);
    if (!volume || !arena || !chart || !chart->nv || !chart->nf) return;
    aff_chart(volume, arena, chart, out, face_samples, 1);
}

int AsmFibre_measure(const AsmRun *run, const char *raw, AsmFibreMoments *out, AsmFibreStats *st)
{
    double started = ves_clock_sec();
    AsmFibreStats s;
    memset(&s, 0, sizeof s);
    if (!run || !raw || !out) return -1;
    memset(out, 0, (run->n_charts ? run->n_charts : 1)*sizeof *out);
    ZarrU8 *z = calloc(1, sizeof *z);
    AffOrder *order = malloc((run->n_charts ? run->n_charts : 1)*sizeof *order);
    if (!z || !order || ZarrU8_open(z, raw)) {
        fprintf(stderr, "[fibre frame] cannot read the CT at %s (an uncompressed uint8 zarr v2 is required)\n", raw);
        free(z); free(order); return -1;
    }
    size_t n = 0;
    for (size_t c = 0; c < run->n_charts; c++) {             /* chunk-major order: the cache keeps its working set */
        const AsmChart *ch = run->charts + c;
        if (!AsmChart_registered(ch) || !ch->nf) continue;
        uint64_t key = 0;
        for (int d = 0; d < 3; d++) {
            double q = floor(ch->centroid[d] / 128.0);
            key = (key << 21) | (uint64_t)(q < 0 ? 0 : q > 2097151 ? 2097151 : q);
        }
        order[n].key = key; order[n].chart = c; n++;
    }
    qsort(order, n, sizeof *order, aff_order);
    Arena_T arena = Arena_new();
    for (size_t k = 0; k < n; k++) {
        size_t c = order[k].chart;
        aff_chart(z, arena, run->charts + c, out + c, NULL, 0);
        s.charts++; s.px += out[c].px;
        if (!out[c].px) s.skipped++;
    }
    s.chunk_misses = z->misses; s.chunk_missing = z->missing;
    Arena_dispose(&arena);
    ZarrU8_close(z); free(z); free(order);
    s.seconds = ves_clock_sec() - started;
    if (st) *st = s;
    return 0;
}

void AsmFibre_transform(AsmFibreMoments *m, int mirror, double cosine, double sine)
{
    if (!m) return;
    if (mirror) { m->m2[1] = -m->m2[1]; m->m4[1] = -m->m4[1]; }     /* u -> -u: phi -> pi - phi */
    double c2 = cosine*cosine - sine*sine, s2 = 2*cosine*sine;     /* e^{i2 rho} */
    double c4 = c2*c2 - s2*s2, s4 = 2*c2*s2;                        /* e^{i4 rho} */
    double a = m->m2[0], b = m->m2[1];
    m->m2[0] = a*c2 - b*s2; m->m2[1] = a*s2 + b*c2;
    a = m->m4[0]; b = m->m4[1];
    m->m4[0] = a*c4 - b*s4; m->m4[1] = a*s4 + b*c4;
}

void AsmFibre_turn(const AsmFibreMoments *m, AsmFibreTurn *o)
{
    if (!o) return;
    memset(o, 0, sizeof *o);
    if (!m || !(m->w > 0) || m->px < ASM_FIBRE_MIN_NODES) return;
    double cross = atan2(m->m4[1], m->m4[0]) / 4;                 /* (-45, 45] deg: the families' normals, so their fibres too */
    double c2 = m->m2[0]*cos(2*cross) + m->m2[1]*sin(2*cross);    /* Re(m2 e^{-i2 cross}) */
    double c4 = m->m4[0]*cos(4*cross) + m->m4[1]*sin(4*cross);    /* |m4| */
    o->coherence = 3*hypot(m->m4[0], m->m4[1]) / m->w;            /* 1 for ideal stripes or an ideal cross, 0 isotropic */
    o->cross_deg = cross * 180 / AFF_PI;
    o->energy_near = m->w - c2 + c4;                               /* D^2 across the fibres at cross (normal cross + 90) */
    o->energy_far = m->w + c2 + c4;                                /* D^2 across the fibres at cross + 90 (normal cross) */
    if (!(o->coherence >= ASM_FIBRE_MIN_COHERENCE)) return;
    double chosen = o->cross_deg;
    if (fabs(o->cross_deg) <= 45 - ASM_FIBRE_TIE_DEG) o->decided_by = 1;
    else {
        o->decided_by = 2;
        if (o->energy_far > o->energy_near) chosen = o->cross_deg > 0 ? o->cross_deg - 90 : o->cross_deg + 90;
    }
    o->turn_deg = -chosen;
}

/* ------------------------------------------------------------------ selftest */

/* A crosshatch in the (x, y) plane, constant in z: the strong family's fibres at alpha + 90 deg, the weak one's at
 * alpha, written as a one-chunk uint8 zarr. */
static int aff_write_volume(const char *root, double alpha_deg)
{
    enum { D = 16, S = 160 };
    char path[1400];
    uint8_t *vol = malloc((size_t)D*S*S);
    if (!vol) return -1;
    double a = alpha_deg*AFF_PI/180, ca = cos(a), sa = sin(a);
    for (int y = 0; y < S; y++) for (int x = 0; x < S; x++) {
        double v = 128 + 50*cos(2*AFF_PI*(x*ca + y*sa)/12.0) + 25*cos(2*AFF_PI*(-x*sa + y*ca)/12.0);
        uint8_t b = (uint8_t)(v < 1 ? 1 : v > 255 ? 255 : v);
        for (int zz = 0; zz < D; zz++) vol[((size_t)zz*S + (size_t)y)*S + (size_t)x] = b;
    }
    ves_mkdir(root);
    snprintf(path, sizeof path, "%s/0", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/.zarray", root);
    FILE *f = fopen(path, "wb");
    if (!f) { free(vol); return -1; }
    fprintf(f, "{\"shape\":[%d,%d,%d],\"chunks\":[%d,%d,%d],\"dtype\":\"|u1\",\"fill_value\":0,\"order\":\"C\",\"filters\":null,"
               "\"dimension_separator\":\"/\",\"compressor\":null,\"zarr_format\":2}", D, S, S, D, S, S);
    fclose(f);
    snprintf(path, sizeof path, "%s/0/0", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/0/0", root); ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/0/0/0", root);
    f = fopen(path, "wb");
    int ok = f && fwrite(vol, 1, (size_t)D*S*S, f) == (size_t)D*S*S;
    if (f && fclose(f)) ok = 0;
    free(vol);
    return ok ? 0 : -1;
}

/* A flat chart in the plane z = 8 over x, y in [20, 140], uv = (x, y). */
static void aff_flat_chart(Arena_T arena, AsmChart *ch)
{
    enum { K = 31 };
    memset(ch, 0, sizeof *ch);
    float *xyz = ARENA_ALLOC(arena, 3*K*K*sizeof *xyz), *nrm = ARENA_ALLOC(arena, 3*K*K*sizeof *nrm), *uv = ARENA_ALLOC(arena, 2*K*K*sizeof *uv);
    int32_t *faces = ARENA_ALLOC(arena, 6*(K-1)*(K-1)*sizeof *faces);
    size_t nf = 0;
    for (int j = 0; j < K; j++) for (int i = 0; i < K; i++) {
        size_t v = (size_t)(j*K + i);
        float x = 20.f + 4.f*(float)i, y = 20.f + 4.f*(float)j;
        xyz[3*v] = 8; xyz[3*v+1] = y; xyz[3*v+2] = x;
        nrm[3*v] = 1; nrm[3*v+1] = 0; nrm[3*v+2] = 0;
        uv[2*v] = x; uv[2*v+1] = y;
    }
    for (int j = 0; j+1 < K; j++) for (int i = 0; i+1 < K; i++) {
        int32_t a = j*K + i, b = a + 1, c = a + K + 1, d = a + K;
        int32_t t[6] = { a, b, c, a, c, d };
        memcpy(faces + 3*nf, t, sizeof t); nf += 2;
    }
    ch->nv = K*K; ch->nf = nf; ch->xyz = xyz; ch->nrm = nrm; ch->uv = uv; ch->faces = faces;
    ch->placed = 1; ch->placement_state = ASM_PLACE_ROOT;
    ch->centroid[0] = 8; ch->centroid[1] = 80; ch->centroid[2] = 80;
}

int AsmFibre_selftest(void)
{
    int fails = 0;
    char dir[512], root[640];
    if (ves_selftest_dir(dir, sizeof dir, "asm_fibre_frame")) return 1;
    Arena_T arena = Arena_new();
    struct { double alpha, cross, turn; int decided; } cases[2] = {
        { 20, 20, -20, 1 },     /* the weak family nearest +u, clearly: lengthwise wins */
        { 40, 40, 50, 2 },      /* a near-diagonal cross: the strong family (fibres at 130 = -50 deg) takes +u */
    };
    for (int k = 0; k < 2; k++) {
        snprintf(root, sizeof root, "%s/cross%d.zarr", dir, k);
        if (aff_write_volume(root, cases[k].alpha)) { fails++; continue; }
        AsmChart ch; aff_flat_chart(arena, &ch);
        AsmRun run; memset(&run, 0, sizeof run); run.arena = arena; run.charts = &ch; run.n_charts = 1;
        AsmFibreMoments m; AsmFibreStats st; AsmFibreTurn t;
        if (AsmFibre_measure(&run, root, &m, &st) || !m.px) { fails++; continue; }
        AsmFibre_turn(&m, &t);
        int bad = fabs(t.cross_deg - cases[k].cross) > 1.0 || t.decided_by != cases[k].decided ||
                  fabs(t.turn_deg - cases[k].turn) > 1.0 || !(t.energy_far > t.energy_near) || !(t.coherence > .2);
        /* the same moments reflected (u -> -u) and turned by 30 deg: cross 20 -> -20 -> 10; 40 -> -40 -> -10 */
        AsmFibreMoments r = m; AsmFibreTurn tr;
        AsmFibre_transform(&r, 1, cos(30*AFF_PI/180), sin(30*AFF_PI/180));
        AsmFibre_turn(&r, &tr);
        double want = k == 0 ? 10 : -10;
        bad |= fabs(tr.cross_deg - want) > 1.0 || fabs(tr.coherence - t.coherence) > 1e-9 ||
               fabs(tr.energy_near - t.energy_near) > 1e-6*m.w || fabs(tr.energy_far - t.energy_far) > 1e-6*m.w;
        fprintf(stderr, "  fibre frame: cross at %.0f deg -> measured %.2f, coherence %.3f, energies near %.3g far %.3g, "
                        "decided %d, turn %.2f; reflected and turned 30 deg: %.2f (%zu px): %s\n",
                cases[k].alpha, t.cross_deg, t.coherence, t.energy_near, t.energy_far, t.decided_by, t.turn_deg,
                tr.cross_deg, m.px, bad ? "FAIL" : "ok");
        fails += bad;
    }
    {   /* no texture: no turn */
        AsmFibreMoments m; memset(&m, 0, sizeof m); m.w = 1; m.px = 100000; AsmFibreTurn t;
        AsmFibre_turn(&m, &t);
        if (t.decided_by || t.turn_deg != 0) { fprintf(stderr, "  fibre frame: a flat image must not turn\n"); fails++; }
    }
    Arena_dispose(&arena);
    if (!fails) ves_selftest_clean(dir);
    return fails;
}
