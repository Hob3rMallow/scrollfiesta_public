/* scroll_raster.c -- u-band-tiled whole-grid strip raster. See scroll_raster.h.
 * The per-pixel semantics mirror Rawtex_write_tif (rawtex_bake.c); that proven
 * single-mesh path is left untouched -- this is the absolute-frame, unbounded-
 * width engine, gated by banded-vs-single-shot byte identity. */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/tiff_io.h"
#include "../common/tif_strip.h"
#include "../common/ves_png.h"
#include "../flatten/rawtex_bake.h"   /* Rawtex_stretch_window */
#include "scroll_raster.h"

void RasterOpts_default(RasterOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    o->du = 1.0;
    o->dv = 1.0;
    o->band_cols = 32768;
    o->range = 2.0;
    o->nsteps = 5;
    o->stretch_ratio = 4.0;
    o->stretch_floor = 25.0;
    /* Source faces are already accepted mesh topology and their valid edge
     * scale follows the upstream remesher (the CVT-coarse 4x5x5 is ~10 vox
     * median, with legitimate edges beyond 20). The historical 6-vox cap was
     * designed for fitted/synthetic hole membranes, not real faces. */
    o->max_edge3d = 0.0;
    o->synth_max_edge3d = 6.0;
    o->pct_lo = 1.0;
    o->pct_hi = 99.0;
    o->write_diag = 1;
    o->write_xyzmap = 1;
    o->strip_w = 0;
    o->synth_f0 = (size_t)-1;   /* memset 0 would mean ALL faces synthetic */
}

static double sr_gray(double v, double lo, double hi)
{
    double g = 0.0;
    if (hi <= lo) return 0.0;
    g = (v - lo) / (hi - lo);
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    return g;
}

static uint8_t sr_norm8(double v, double lo, double hi)
{
    return (uint8_t)(sr_gray(v, lo, hi) * 255.0 + 0.5);
}

/* RGB analog of tif_strip's stacker (PNG only -- TiffIO_save is grayscale):
 * slice W x H `rgb` (3 B/px) into ceil(W/sw) vertical strips stacked with `gap`
 * navy rows between, write <prefix>_xyzmap_strip.png. Debug convenience: any
 * failure is silently ignored (never load-bearing). */
static void sr_strip_rgb(const char *path, const uint8_t *rgb,
                         size_t W, size_t H, int target_strip_w)
{
    int gap = 6;
    const uint8_t navy[3] = { 12, 14, 40 };
    if (path == NULL || rgb == NULL || W == 0 || H == 0) return;
    if (target_strip_w <= 0) target_strip_w = 3000;
    size_t nstrips = (W + (size_t)target_strip_w - 1) / (size_t)target_strip_w;
    if (nstrips < 1) nstrips = 1;
    size_t sw = (W + nstrips - 1) / nstrips;
    size_t oh = nstrips * (H + (size_t)gap) - (size_t)gap;
    if (oh < H) oh = H;
    if (oh > 100000000 || sw * oh > ((size_t)1 << 30)) return;   /* runaway */
    uint8_t *out = (uint8_t *)malloc(sw * oh * 3);
    if (out == NULL) return;
    for (size_t i = 0; i < sw * oh; i++) {
        out[i * 3 + 0] = navy[0];
        out[i * 3 + 1] = navy[1];
        out[i * 3 + 2] = navy[2];
    }
    for (size_t s = 0; s < nstrips; s++) {
        size_t x0 = s * sw, oy0 = s * (H + (size_t)gap);
        for (size_t y = 0; y < H; y++) {
            for (size_t x = 0; x < sw; x++) {
                size_t sx = x0 + x;
                if (sx >= W) continue;              /* ragged last strip */
                size_t o = ((oy0 + y) * sw + x) * 3;
                size_t p = (y * W + sx) * 3;
                out[o + 0] = rgb[p + 0];
                out[o + 1] = rgb[p + 1];
                out[o + 2] = rgb[p + 2];
            }
        }
    }
    VesPng_write_rgb(path, out, (int)sw, (int)oh);
    free(out);
}

/* diagclass codes (rawtex_bake convention; 5 = its "dark", reserved) */
enum { SR_BG = 0, SR_OK = 1, SR_SKIPUV = 2, SR_SKIP3D = 3, SR_MULTI = 4,
       SR_SKIPOWN = 6, SR_FILL = 7 };

/* Splat faces bface[from..to) into the band accumulators. Real faces
 * (synth=0) mark their gate class into cuv/c3d/cown so uncovered px keep a
 * diagnostic; synthetic faces (synth=1) are pure paint -- a gated synth face
 * is skipped outright so it can neither paint nor overwrite a diagnostic. */
typedef struct {
    const PieceSet *ps;
    const uint8_t *gate;
    CubeTable *ct;
    const RasterOpts *opts;
    long X0, Y0;
    size_t H;
    double *sum;
    uint32_t *cnt;
    uint8_t *cuv, *c3d, *cown;
    uint8_t *touch;   /* synth pass only: geometry-present marker (px inside
                         a synth face even when the RAW sample fails) */
    float *xsum;      /* [bw*H*3] world-pos accum for the xyzmap (NULL = off).
                         Accumulated for EVERY covered px, sampled or not, so
                         mid-air ribbon still gets a position color. */
    uint32_t *xcnt;   /* [bw*H] cover count for xsum */
} SrBand;

static void sr_splat(const SrBand *sb, const int32_t *bface,
                     int32_t from, int32_t to,
                     size_t bx0, size_t bw, int synth)
{
    const PieceSet *ps = sb->ps;
    const RasterOpts *opts = sb->opts;
    const double du = opts->du, dv = opts->dv;
    const long X0 = sb->X0, Y0 = sb->Y0;
    const size_t H = sb->H;
    for (int32_t bf = from; bf < to; bf++) {
        size_t f = (size_t)bface[bf];
        int g = sb->gate[f];
        if (synth && g != 0) continue;
        size_t a = (size_t)ps->faces[f * 3 + 0];
        size_t b = (size_t)ps->faces[f * 3 + 1];
        size_t c = (size_t)ps->faces[f * 3 + 2];
        double ua = (double)ps->uv[a * 2 + 0] / du - (double)X0;
        double va = (double)ps->uv[a * 2 + 1] / dv - (double)Y0;
        double ub = (double)ps->uv[b * 2 + 0] / du - (double)X0;
        double vb = (double)ps->uv[b * 2 + 1] / dv - (double)Y0;
        double uc = (double)ps->uv[c * 2 + 0] / du - (double)X0;
        double vc = (double)ps->uv[c * 2 + 1] / dv - (double)Y0;
        double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        if (fabs(A2) < 1e-12) continue;
        double lox = ua < ub ? (ua < uc ? ua : uc) : (ub < uc ? ub : uc);
        double hix = ua > ub ? (ua > uc ? ua : uc) : (ub > uc ? ub : uc);
        double loy = va < vb ? (va < vc ? va : vc) : (vb < vc ? vb : vc);
        double hiy = va > vb ? (va > vc ? va : vc) : (vb > vc ? vb : vc);
        long x0 = (long)ceil(lox - 0.001), x1 = (long)floor(hix + 0.001);
        long y0 = (long)ceil(loy - 0.001), y1 = (long)floor(hiy + 0.001);
        if (x0 < (long)bx0) x0 = (long)bx0;
        if (x1 >= (long)(bx0 + bw)) x1 = (long)(bx0 + bw) - 1;
        if (y0 < 0) y0 = 0;
        if (y1 >= (long)H) y1 = (long)H - 1;
        for (long yy = y0; yy <= y1; yy++) {
            for (long xx = x0; xx <= x1; xx++) {
                double pu = (double)xx, pv = (double)yy;
                double l1 = ((pu - ua) * (vc - va) - (pv - va) * (uc - ua)) / A2;
                double l2 = ((ub - ua) * (pv - va) - (vb - va) * (pu - ua)) / A2;
                double l0 = 1.0 - l1 - l2;
                if (l0 < -1e-6 || l1 < -1e-6 || l2 < -1e-6) continue;
                size_t bi = (size_t)yy * bw + ((size_t)xx - bx0);
                if (g == SR_SKIPOWN) { sb->cown[bi] = 1; continue; }
                if (g == SR_SKIPUV) { sb->cuv[bi] = 1; continue; }
                if (g == SR_SKIP3D) { sb->c3d[bi] = 1; continue; }
                if (sb->touch != NULL) sb->touch[bi] = 1;
                float p3[3], n3[3];
                double nn = 0.0;
                for (int k = 0; k < 3; k++) {
                    p3[k] = (float)(l0 * (double)ps->verts[a * 3 + k]
                                    + l1 * (double)ps->verts[b * 3 + k]
                                    + l2 * (double)ps->verts[c * 3 + k]);
                    n3[k] = (float)(l0 * (double)ps->normals[a * 3 + k]
                                    + l1 * (double)ps->normals[b * 3 + k]
                                    + l2 * (double)ps->normals[c * 3 + k]);
                }
                nn = sqrt((double)n3[0] * n3[0] + (double)n3[1] * n3[1]
                          + (double)n3[2] * n3[2]);
                if (nn > 1e-6) {
                    n3[0] = (float)((double)n3[0] / nn);
                    n3[1] = (float)((double)n3[1] / nn);
                    n3[2] = (float)((double)n3[2] / nn);
                }
                if (sb->xsum != NULL) {   /* xyz before the sample gate: color
                                             even where RAW is missing */
                    sb->xsum[bi * 3 + 0] += p3[0];
                    sb->xsum[bi * 3 + 1] += p3[1];
                    sb->xsum[bi * 3 + 2] += p3[2];
                    sb->xcnt[bi]++;
                }
                double s = sample_vertex(sb->ct, p3, nn > 0.5 ? n3 : NULL,
                                         opts->range, opts->nsteps);
                if (s < 0.0) continue;
                sb->sum[bi] += s;
                sb->cnt[bi]++;
            }
        }
    }
}

int ScrollRaster_run(Arena_T arena, const PieceSet *ps,
                     const char *raw_dir, long raw_chunk,
                     const char *prefix, const RasterOpts *opts,
                     RasterStats *out)
{
    assert(arena && ps && prefix && opts && out);
    memset(out, 0, sizeof(*out));
    if (ps->nf == 0 || ps->nv == 0 || opts->du <= 0.0 || opts->dv <= 0.0)
        return -1;

    const double du = opts->du, dv = opts->dv;
    /* canvas origin/pixel mapping: px x = round(u/du) - X0 */
    long X0 = (long)floor(ps->u_min / du + 0.5);
    long Y0 = (long)floor(ps->v_min / dv + 0.5);
    long X1 = (long)floor(ps->u_max / du + 0.5);
    long Y1 = (long)floor(ps->v_max / dv + 0.5);
    size_t W = (size_t)(X1 - X0) + 1;
    size_t H = (size_t)(Y1 - Y0) + 1;
    if (W < 1 || H < 1 || W > (size_t)1 << 31 || H > (size_t)1 << 24 ||
        W * H > (size_t)6 << 30) {
        fprintf(stderr, "scroll_raster: ERROR canvas %zux%zu unreasonable "
                "(raise --du/--dv?)\n", W, H);
        return -1;
    }
    out->W = W;
    out->H = H;
    out->u0 = (double)X0 * du;
    out->v0 = (double)Y0 * dv;

    size_t B = (opts->band_cols == 0 || opts->band_cols >= W) ? W
                                                              : opts->band_cols;
    size_t n_bands = (W + B - 1) / B;
    out->n_bands = n_bands;

    /* RAW table: the caller's shared one, else a lazy local for this run */
    CubeTable ct_local;
    CubeTable *ct = opts->ct;
    if (ct == NULL) {
        if (cubetable_init(&ct_local, arena, raw_dir, raw_chunk, ps->verts,
                           ps->nv, opts->range + 2.0) != 0) {
            fprintf(stderr, "scroll_raster: ERROR cubetable_init failed (%s)\n",
                    raw_dir);
            return -1;
        }
        ct = &ct_local;
    }

    /* contrast window over a vertex subsample (same recipe as obj_bake_raw).
     * NO Arena_save/restore around this: sample_vertex lazily loads RAW cube
     * volumes into ct->arena == arena, and a restore here would free them
     * while ct->slot[] still points at them -- the band splat then reads
     * recycled memory (silent garbage paint at small scale, AV once the
     * freed chunks exceed the arena freelist). The few MB of val/has scratch
     * just stay in the arena. A caller re-baking a stage passes lo/hi_fixed
     * so every stage TIF shares one window and stays differenceable. */
    if (opts->hi_fixed > opts->lo_fixed) {
        out->lo = opts->lo_fixed;
        out->hi = opts->hi_fixed;
    } else {
        size_t target = 400000;
        size_t stride = ps->nv > target ? ps->nv / target : 1;
        size_t nprobe = ps->nv / stride + 1;
        double *val = (double *)ARENA_ALLOC(arena, nprobe * sizeof(double));
        uint8_t *has = (uint8_t *)ARENA_CALLOC(arena, nprobe, 1);
        size_t q = 0;
        for (size_t i = 0; i < ps->nv; i += stride) {
            double s = sample_vertex(ct, &ps->verts[i * 3],
                                     &ps->normals[i * 3], opts->range,
                                     opts->nsteps);
            if (q >= nprobe) break;
            val[q] = s < 0.0 ? 0.0 : s;
            has[q] = s >= 0.0;
            q++;
        }
        Rawtex_stretch_window(val, has, q, opts->pct_lo, opts->pct_hi,
                              &out->lo, &out->hi);
    }

    /* full-size canvases (uint8), band-size accumulators */
    int have_synth = opts->synth_f0 != (size_t)-1;
    uint8_t *img = (uint8_t *)calloc(W * H, 1);
    uint8_t *cls = opts->write_diag ? (uint8_t *)calloc(W * H, 1) : NULL;
    double *sum = (double *)malloc(B * H * sizeof(double));
    uint32_t *cnt = (uint32_t *)malloc(B * H * sizeof(uint32_t));
    uint8_t *cuv = (uint8_t *)malloc(B * H);
    uint8_t *c3d = (uint8_t *)malloc(B * H);
    uint8_t *cown = (uint8_t *)malloc(B * H);
    uint8_t *occ = have_synth ? (uint8_t *)malloc(B * H) : NULL;
    uint8_t *touch = have_synth ? (uint8_t *)malloc(B * H) : NULL;
    if (img == NULL || sum == NULL || cnt == NULL || cuv == NULL ||
        c3d == NULL || cown == NULL || (opts->write_diag && cls == NULL) ||
        (have_synth && (occ == NULL || touch == NULL))) {
        fprintf(stderr, "scroll_raster: ERROR out of memory (canvas %zux%zu)\n",
                W, H);
        free(img); free(cls); free(sum); free(cnt); free(cuv); free(c3d);
        free(cown); free(occ); free(touch);
        return -1;
    }

    /* xyzmap: full-size RGB canvas + band-size world-pos accumulators. Non-
     * fatal -- a grid too big for these just skips the map (rawtex still bakes).
     * xocc marks real xyz-claim so the synth pass never overpaints real
     * geometry (only needed when a synth pass exists). */
    int have_xyz = opts->write_xyzmap;
    uint8_t *xyz = NULL;
    float *xsum = NULL;
    uint32_t *xcnt = NULL;
    uint8_t *xocc = NULL;
    if (have_xyz) {
        xyz = (uint8_t *)calloc(W * H * 3, 1);
        xsum = (float *)malloc(B * H * 3 * sizeof(float));
        xcnt = (uint32_t *)malloc(B * H * sizeof(uint32_t));
        xocc = have_synth ? (uint8_t *)malloc(B * H) : NULL;
        if (xyz == NULL || xsum == NULL || xcnt == NULL ||
            (have_synth && xocc == NULL)) {
            fprintf(stderr, "scroll_raster: WARN xyzmap disabled (out of "
                    "memory, canvas %zux%zu)\n", W, H);
            free(xyz); free(xsum); free(xcnt); free(xocc);
            xyz = NULL; xsum = NULL; xcnt = NULL; xocc = NULL;
            have_xyz = 0;
        }
    }
    if (have_xyz) {   /* normalization bbox over REAL-face verts (z,y,x) */
        double lo[3] = { 1e300, 1e300, 1e300 };
        double hi[3] = { -1e300, -1e300, -1e300 };
        size_t nreal = ps->nf < opts->synth_f0 ? ps->nf : opts->synth_f0;
        for (size_t f = 0; f < nreal; f++)
            for (int e = 0; e < 3; e++) {
                size_t v = (size_t)ps->faces[f * 3 + e];
                for (int k = 0; k < 3; k++) {
                    double cc = (double)ps->verts[v * 3 + k];
                    if (cc < lo[k]) lo[k] = cc;
                    if (cc > hi[k]) hi[k] = cc;
                }
            }
        for (int k = 0; k < 3; k++) {
            if (hi[k] <= lo[k]) hi[k] = lo[k] + 1.0;
            out->xyz_lo[k] = lo[k];
            out->xyz_hi[k] = hi[k];
        }
    }

    /* bin faces to bands by their u-bbox (a face may straddle two bands) */
    Arena_Mark mark = Arena_save(arena);
    int32_t *boff = (int32_t *)ARENA_CALLOC(arena, n_bands + 2, sizeof(int32_t));
    long *fb0 = (long *)ARENA_ALLOC(arena, (ps->nf + 1) * sizeof(long));
    long *fb1 = (long *)ARENA_ALLOC(arena, (ps->nf + 1) * sizeof(long));
    for (size_t f = 0; f < ps->nf; f++) {
        double ulo = 1e300, uhi = -1e300;
        for (int e = 0; e < 3; e++) {
            double u = (double)ps->uv[(size_t)ps->faces[f * 3 + e] * 2 + 0];
            if (u < ulo) ulo = u;
            if (u > uhi) uhi = u;
        }
        long px0 = (long)ceil(ulo / du - 0.001) - X0;
        long px1 = (long)floor(uhi / du + 0.001) - X0;
        if (px0 < 0) px0 = 0;
        if (px1 >= (long)W) px1 = (long)W - 1;
        if (px1 < px0) { fb0[f] = 1; fb1[f] = 0; continue; }   /* no pixels */
        fb0[f] = px0 / (long)B;
        fb1[f] = px1 / (long)B;
        for (long b = fb0[f]; b <= fb1[f]; b++) boff[b + 1]++;
    }
    for (size_t b = 0; b < n_bands; b++) boff[b + 1] += boff[b];
    int32_t *bface = (int32_t *)ARENA_ALLOC(arena,
                                            ((size_t)boff[n_bands] + 1) *
                                            sizeof(int32_t));
    {
        int32_t *cur = (int32_t *)ARENA_ALLOC(arena,
                                              (n_bands + 1) * sizeof(int32_t));
        memcpy(cur, boff, n_bands * sizeof(int32_t));
        for (size_t f = 0; f < ps->nf; f++)
            for (long b = fb0[f]; b <= fb1[f]; b++)
                bface[cur[b]++] = (int32_t)f;
    }

    /* face gates are per-face (band-independent): precompute once so a face
     * straddling a band boundary is judged identically in both bands */
    uint8_t *gate = (uint8_t *)ARENA_CALLOC(arena, ps->nf + 1, 1);
    for (size_t f = 0; f < ps->nf; f++) {
        int real = f < opts->synth_f0;   /* synth gates apply but don't count */
        double edge_cap = real ? opts->max_edge3d : opts->synth_max_edge3d;
        if (opts->face_keep != NULL && !opts->face_keep[f]) {
            gate[f] = SR_SKIPOWN;
            if (real) out->skip_own++;
            continue;
        }
        size_t a = (size_t)ps->faces[f * 3 + 0];
        size_t b = (size_t)ps->faces[f * 3 + 1];
        size_t c = (size_t)ps->faces[f * 3 + 2];
        size_t vi[3];
        vi[0] = a; vi[1] = b; vi[2] = c;
        int bad_uv = 0, bad_3d = 0;
        for (int e = 0; e < 3; e++) {
            size_t p = vi[e], q = vi[(e + 1) % 3];
            double d3 = 0.0;
            for (int k = 0; k < 3; k++) {
                double dd = (double)ps->verts[q * 3 + k]
                            - (double)ps->verts[p * 3 + k];
                d3 += dd * dd;
            }
            d3 = sqrt(d3);
            double due = (double)ps->uv[q * 2 + 0] - (double)ps->uv[p * 2 + 0];
            double dve = (double)ps->uv[q * 2 + 1] - (double)ps->uv[p * 2 + 1];
            double e2d = sqrt(due * due + dve * dve);
            if (edge_cap > 0.0 && d3 > edge_cap) bad_3d = 1;
            if (opts->stretch_ratio > 0.0
                && e2d > (opts->stretch_ratio * d3 > opts->stretch_floor
                          ? opts->stretch_ratio * d3 : opts->stretch_floor))
                bad_uv = 1;
        }
        if (bad_3d) { gate[f] = SR_SKIP3D; if (real) out->skip_3d++; }
        else if (bad_uv) { gate[f] = SR_SKIPUV; if (real) out->skip_uv++; }
    }

    SrBand sb;
    sb.ps = ps;
    sb.gate = gate;
    sb.ct = ct;
    sb.opts = opts;
    sb.X0 = X0;
    sb.Y0 = Y0;
    sb.H = H;
    sb.sum = sum;
    sb.cnt = cnt;
    sb.cuv = cuv;
    sb.c3d = c3d;
    sb.cown = cown;
    sb.touch = NULL;   /* set per synth pass only */
    sb.xsum = have_xyz ? xsum : NULL;
    sb.xcnt = xcnt;

    double t0 = ves_clock_sec();
    for (size_t band = 0; band < n_bands; band++) {
        size_t bx0 = band * B;
        size_t bw = (bx0 + B <= W) ? B : W - bx0;
        memset(sum, 0, bw * H * sizeof(double));
        memset(cnt, 0, bw * H * sizeof(uint32_t));
        memset(cuv, 0, bw * H);
        memset(c3d, 0, bw * H);
        memset(cown, 0, bw * H);
        if (have_xyz) {
            memset(xsum, 0, bw * H * 3 * sizeof(float));
            memset(xcnt, 0, bw * H * sizeof(uint32_t));
        }

        /* real/synth split point: bface is filled by ascending f per band */
        int32_t bf_lo = boff[band], bf_hi = boff[band + 1];
        int32_t bf_split = bf_hi;
        if (have_synth) {
            int32_t a = bf_lo, b = bf_hi;
            while (a < b) {
                int32_t m = a + (b - a) / 2;
                if ((size_t)bface[m] >= opts->synth_f0) b = m;
                else a = m + 1;
            }
            bf_split = a;
        }

        sr_splat(&sb, bface, bf_lo, bf_split, bx0, bw, 0);

        /* collapse the real band into the canvases */
        for (size_t yy = 0; yy < H; yy++) {
            for (size_t x = 0; x < bw; x++) {
                size_t bi = yy * bw + x;
                size_t pi = yy * W + (bx0 + x);
                if (cnt[bi] > 0) {
                    double gy = sr_gray(sum[bi] / (double)cnt[bi],
                                        out->lo, out->hi);
                    img[pi] = (uint8_t)(gy * 255.0 + 0.5);
                    out->filled++;
                    if (cnt[bi] > 1) out->multi++;
                    if (cls != NULL)
                        cls[pi] = cnt[bi] > 1 ? SR_MULTI : SR_OK;
                } else if (cls != NULL) {
                    if (cown[bi]) cls[pi] = SR_SKIPOWN;
                    else if (cuv[bi]) cls[pi] = SR_SKIPUV;
                    else if (c3d[bi]) cls[pi] = SR_SKIP3D;
                }
                if (have_synth)
                    occ[bi] = (uint8_t)(cnt[bi] > 0 || cuv[bi] || c3d[bi]
                                        || cown[bi]);
                if (have_xyz && xcnt[bi] > 0) {
                    double inv = 1.0 / (double)xcnt[bi];
                    xyz[pi * 3 + 0] = sr_norm8((double)xsum[bi * 3 + 2] * inv,
                                               out->xyz_lo[2], out->xyz_hi[2]);
                    xyz[pi * 3 + 1] = sr_norm8((double)xsum[bi * 3 + 1] * inv,
                                               out->xyz_lo[1], out->xyz_hi[1]);
                    xyz[pi * 3 + 2] = sr_norm8((double)xsum[bi * 3 + 0] * inv,
                                               out->xyz_lo[0], out->xyz_hi[0]);
                }
                if (have_xyz && xocc != NULL)
                    xocc[bi] = (uint8_t)(occ[bi] || xcnt[bi] > 0);
            }
        }

        /* synth pass: paint ONLY px with no real claim (paint or gate).
         * Geometry-touched px whose RAW sample failed (outside the loaded
         * cubes) still get class 7 with gray 0 -- the green tint then reads
         * as "ribbon exists here", sampled or not. */
        if (bf_split < bf_hi) {
            memset(sum, 0, bw * H * sizeof(double));
            memset(cnt, 0, bw * H * sizeof(uint32_t));
            memset(touch, 0, bw * H);
            if (have_xyz) {
                memset(xsum, 0, bw * H * 3 * sizeof(float));
                memset(xcnt, 0, bw * H * sizeof(uint32_t));
            }
            sb.touch = touch;
            sr_splat(&sb, bface, bf_split, bf_hi, bx0, bw, 1);
            sb.touch = NULL;
            for (size_t yy = 0; yy < H; yy++) {
                for (size_t x = 0; x < bw; x++) {
                    size_t bi = yy * bw + x;
                    /* color synth-only px (real xyz-claim wins): before the
                     * occ-continue, since a sample-failed real px has occ=0
                     * but must still block synth xyz overpaint via xocc */
                    if (have_xyz && xocc != NULL && !xocc[bi] && xcnt[bi] > 0) {
                        size_t pi = yy * W + (bx0 + x);
                        double inv = 1.0 / (double)xcnt[bi];
                        xyz[pi * 3 + 0] = sr_norm8((double)xsum[bi * 3 + 2] * inv,
                                                   out->xyz_lo[2], out->xyz_hi[2]);
                        xyz[pi * 3 + 1] = sr_norm8((double)xsum[bi * 3 + 1] * inv,
                                                   out->xyz_lo[1], out->xyz_hi[1]);
                        xyz[pi * 3 + 2] = sr_norm8((double)xsum[bi * 3 + 0] * inv,
                                                   out->xyz_lo[0], out->xyz_hi[0]);
                    }
                    if (occ[bi]) continue;
                    if (cnt[bi] > 0) {
                        size_t pi = yy * W + (bx0 + x);
                        double gy = sr_gray(sum[bi] / (double)cnt[bi],
                                            out->lo, out->hi);
                        img[pi] = (uint8_t)(gy * 255.0 + 0.5);
                        if (cls != NULL) cls[pi] = SR_FILL;
                        out->filled_synth++;
                    } else if (touch[bi]) {
                        if (cls != NULL)
                            cls[yy * W + (bx0 + x)] = SR_FILL;
                        out->synth_nosample++;
                    }
                }
            }
        }
        if (n_bands > 1)
            fprintf(stderr, "  [raster] band %zu/%zu (cols %zu..%zu) "
                    "raw cubes %d (%.1fs)\n", band + 1, n_bands, bx0,
                    bx0 + bw - 1, ct->n_loaded, ves_clock_sec() - t0);
    }
    Arena_restore(arena, mark);
    out->fill = (double)out->filled / (double)(W * H);
    out->raw_loaded = ct->n_loaded;
    out->raw_missing = ct->n_missing;

    /* outputs */
    int rc = 0;
    char path[2048];
    snprintf(path, sizeof(path), "%s_rawtex.tif", prefix);
    if (TiffIO_save(path, img, 1, (int)H, (int)W) != 0) {
        fprintf(stderr, "scroll_raster: ERROR writing %s\n", path);
        rc = -1;
    }
    if (rc == 0 && W <= 100000) {
        snprintf(path, sizeof(path), "%s_rawtex.png", prefix);
        VesPng_write_gray(path, img, (int)W, (int)H);
    }
    if (rc == 0) {
        snprintf(path, sizeof(path), "%s_rawtex_strip.tif", prefix);
        TifStrip_write_gray(path, img, (int)W, (int)H, opts->strip_w, -1);
    }
    if (rc == 0 && cls != NULL) {
        snprintf(path, sizeof(path), "%s_diagclass.tif", prefix);
        if (TiffIO_save(path, cls, 1, (int)H, (int)W) != 0)
            fprintf(stderr, "scroll_raster: WARN diagclass write failed\n");
    }
    if (rc == 0 && have_xyz) {
        if (W <= 100000) {
            snprintf(path, sizeof(path), "%s_xyzmap.png", prefix);
            VesPng_write_rgb(path, xyz, (int)W, (int)H);
        }
        snprintf(path, sizeof(path), "%s_xyzmap_strip.png", prefix);
        sr_strip_rgb(path, xyz, W, H, opts->strip_w);
    }
    free(img);
    free(cls);
    free(sum);
    free(cnt);
    free(cuv);
    free(c3d);
    free(cown);
    free(occ);
    free(touch);
    free(xyz);
    free(xsum);
    free(xcnt);
    free(xocc);
    return rc;
}

/* ============================================================================
 * Self-test: injected CubeTable, banding byte-identity, gates.
 * ==========================================================================*/

static int sr_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[scroll_raster selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* slurp a whole file into a malloc'd buffer (caller frees); NULL on failure.
 * stb's PNG encoder is deterministic, so identical pixels -> identical bytes,
 * which lets the xyzmap banding-identity check compare PNGs without a decoder. */
static uint8_t *sr_slurp(const char *path, size_t *out_len)
{
    FILE *f = fopen(path, "rb");
    *out_len = 0;
    if (f == NULL) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long n = ftell(f);
    if (n < 0) { fclose(f); return NULL; }
    rewind(f);
    uint8_t *buf = (uint8_t *)malloc((size_t)n + 1);
    if (buf == NULL) { fclose(f); return NULL; }
    size_t got = fread(buf, 1, (size_t)n, f);
    fclose(f);
    if (got != (size_t)n) { free(buf); return NULL; }
    *out_len = got;
    return buf;
}

/* run the raster over a synthetic 2-triangle piece set into a scratch dir;
 * returns the canvas via out (arena) for comparisons */
static int sr_run_case(Arena_T arena, size_t band_cols, const char *tag,
                       uint8_t **out_img, size_t *out_W, size_t *out_H,
                       RasterStats *st)
{
    /* a 20x20-vox flat quad at z=10..30, y=20, x=10..30 mapped to
     * u=[-1000,-980], v=[10,30] (absolute negative frame on purpose) */
    PieceSet ps;
    memset(&ps, 0, sizeof(ps));
    float verts[4 * 3] = {
        10, 20, 10,   30, 20, 10,   10, 20, 30,   30, 20, 30
    };
    float uv[4 * 2] = {
        -1000, 10,   -1000, 30,   -980, 10,   -980, 30
    };
    float normals[4 * 3] = {
        0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0
    };
    int32_t faces[2 * 3] = { 0, 1, 2, 1, 3, 2 };
    int32_t face_cube[2] = { 0, 0 };
    ps.verts = verts;
    ps.uv = uv;
    ps.normals = normals;
    ps.faces = faces;
    ps.face_cube = face_cube;
    ps.nv = 4;
    ps.nf = 2;
    ps.n_cubes = 1;
    ps.u_min = -1000; ps.u_max = -980;
    ps.v_min = 10; ps.v_max = 30;

    RasterOpts o;
    RasterOpts_default(&o);
    o.band_cols = band_cols;
    o.write_diag = 1;
    o.pct_lo = 0.0;      /* keep the tiny sample's window deterministic */
    o.pct_hi = 100.0;
    o.max_edge3d = 100.0;   /* the synthetic quad has 20-vox edges; the 6-vox
                             * membrane gate (tuned for real ~2-4 vox meshes)
                             * would rightly skip it */

    char prefix[1024];
    snprintf(prefix, sizeof(prefix), "output/_selftest_scroll_unroll/%s", tag);
    ves_ensure_parent_dir(prefix);

    /* write the synthetic gradient volume to disk once as a real cube TIFF
     * (ScrollRaster_run builds its own lazy CubeTable from the dir; this also
     * exercises the TiffIO round trip) */
    static int vol_written = 0;
    if (!vol_written) {
        Arena_Mark mark = Arena_save(arena);
        long chunk = 64;
        uint8_t *vol = (uint8_t *)ARENA_ALLOC(arena, (size_t)(chunk * chunk * chunk));
        for (long z = 0; z < chunk; z++)
            for (long y = 0; y < chunk; y++)
                for (long x = 0; x < chunk; x++)
                    vol[z * chunk * chunk + y * chunk + x] =
                        (uint8_t)(y * 255 / (chunk - 1));
        char volpath[1024];
        snprintf(volpath, sizeof(volpath),
                 "output/_selftest_scroll_unroll/raw/z00000_y00000_x00000.tif");
        ves_ensure_parent_dir(volpath);
        if (TiffIO_save(volpath, vol, (int)chunk, (int)chunk, (int)chunk) != 0)
            return -1;
        Arena_restore(arena, mark);
        vol_written = 1;
    }

    int rc = ScrollRaster_run(arena, &ps, "output/_selftest_scroll_unroll/raw",
                              64, prefix, &o, st);
    if (rc != 0) return rc;

    /* read the canvas back for byte comparisons */
    char path[1200];
    snprintf(path, sizeof(path), "%s_rawtex.tif", prefix);
    uint8_t *img = NULL;
    int d = 0, h = 0, w = 0;
    if (TiffIO_load(arena, path, &img, &d, &h, &w) != 0) return -1;
    *out_img = img;
    *out_W = (size_t)w;
    *out_H = (size_t)h;
    return 0;
}

int ScrollRaster_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    uint8_t *img1 = NULL, *img2 = NULL;
    size_t W1 = 0, H1 = 0, W2 = 0, H2 = 0;
    RasterStats s1, s2;
    int rc1 = sr_run_case(arena, 0, "single", &img1, &W1, &H1, &s1);
    int rc2 = sr_run_case(arena, 7, "banded", &img2, &W2, &H2, &s2);
    fprintf(stderr, "[scroll_raster selftest] single rc=%d %zux%zu fill=%.2f "
            "(lo=%.0f hi=%.0f) | banded rc=%d %zux%zu fill=%.2f bands=%zu\n",
            rc1, W1, H1, s1.fill, s1.lo, s1.hi, rc2, W2, H2, s2.fill,
            s2.n_bands);
    sr_check(rc1 == 0 && rc2 == 0, "raster rc", &fails);
    sr_check(W1 == 21 && H1 == 21, "canvas dims from absolute uv", &fails);
    sr_check(s1.fill > 0.9, "quad covers the canvas", &fails);
    sr_check(s2.n_bands == 3, "banded run banded", &fails);
    sr_check(W1 == W2 && H1 == H2 &&
             memcmp(img1, img2, W1 * H1) == 0,
             "banding byte-identity", &fails);
    /* the volume value equals y=20 everywhere on the quad (gradient = y):
     * expect a flat mid-gray strip -- verify the center pixel is nonzero and
     * uniform along u */
    if (rc1 == 0 && W1 == 21 && H1 == 21) {
        uint8_t center = img1[10 * W1 + 10];
        int uniform = 1;
        for (size_t x = 2; x < W1 - 2; x++)
            if (img1[10 * W1 + x] != center) uniform = 0;
        sr_check(center > 0 && uniform, "flat quad samples uniformly", &fails);
    }

    /* xyzmap: single-shot vs banded PNG byte-identical (stb deterministic),
     * and non-black (the quad spans x=10..30, z=10..30 so R and B vary) */
    {
        size_t l1 = 0, l2 = 0;
        uint8_t *p1 = sr_slurp("output/_selftest_scroll_unroll/single_xyzmap.png",
                               &l1);
        uint8_t *p2 = sr_slurp("output/_selftest_scroll_unroll/banded_xyzmap.png",
                               &l2);
        sr_check(p1 != NULL && p2 != NULL, "xyzmap png written", &fails);
        if (p1 != NULL && p2 != NULL)
            sr_check(l1 == l2 && memcmp(p1, p2, l1) == 0,
                     "xyzmap banding byte-identity", &fails);
        int nonblack = 0;
        for (size_t i = 0; p1 != NULL && i < l1; i++)
            if (p1[i] != 0) { nonblack = 1; break; }
        sr_check(nonblack, "xyzmap non-black", &fails);
        free(p1);
        free(p2);
    }

    /* write_xyzmap = 0 emits no xyzmap file */
    {
        RasterStats st;
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        float verts[4 * 3] = { 10, 20, 10, 30, 20, 10, 10, 20, 30, 30, 20, 30 };
        float uv[4 * 2] = { -1000, 10, -1000, 30, -980, 10, -980, 30 };
        float nrm[4 * 3] = { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 };
        int32_t faces[2 * 3] = { 0, 1, 2, 1, 3, 2 };
        int32_t fc[2] = { 0, 0 };
        ps.verts = verts; ps.uv = uv; ps.normals = nrm;
        ps.faces = faces; ps.face_cube = fc;
        ps.nv = 4; ps.nf = 2; ps.n_cubes = 1;
        ps.u_min = -1000; ps.u_max = -980; ps.v_min = 10; ps.v_max = 30;
        RasterOpts o;
        RasterOpts_default(&o);
        o.write_xyzmap = 0;
        o.pct_lo = 0.0; o.pct_hi = 100.0; o.max_edge3d = 100.0;
        int rc = ScrollRaster_run(arena, &ps,
                                  "output/_selftest_scroll_unroll/raw", 64,
                                  "output/_selftest_scroll_unroll/noxyz", &o,
                                  &st);
        size_t nl = 0;
        uint8_t *nb = sr_slurp("output/_selftest_scroll_unroll/noxyz_xyzmap.png",
                               &nl);
        sr_check(rc == 0 && nb == NULL, "write_xyzmap=0 emits nothing", &fails);
        free(nb);
    }

    /* face_keep: drop one of the quad's two triangles -> half fill, class 6
     * where only the dropped face covered */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        float verts[4 * 3] = {
            10, 20, 10,   30, 20, 10,   10, 20, 30,   30, 20, 30
        };
        float uv[4 * 2] = { -1000, 10, -1000, 30, -980, 10, -980, 30 };
        float normals[4 * 3] = { 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 };
        int32_t faces[2 * 3] = { 0, 1, 2, 1, 3, 2 };
        int32_t fc[2] = { 0, 0 };
        ps.verts = verts; ps.uv = uv; ps.normals = normals;
        ps.faces = faces; ps.face_cube = fc;
        ps.nv = 4; ps.nf = 2; ps.n_cubes = 1;
        ps.u_min = -1000; ps.u_max = -980; ps.v_min = 10; ps.v_max = 30;
        RasterOpts o;
        RasterOpts_default(&o);
        o.band_cols = 0;
        o.pct_lo = 0.0; o.pct_hi = 100.0;
        o.max_edge3d = 100.0;
        uint8_t km[2] = { 1, 0 };
        o.face_keep = km;
        RasterStats st;
        int rc = ScrollRaster_run(arena, &ps,
                                  "output/_selftest_scroll_unroll/raw", 64,
                                  "output/_selftest_scroll_unroll/own", &o,
                                  &st);
        sr_check(rc == 0 && st.skip_own == 1, "face_keep counts the drop",
                 &fails);
        sr_check(st.filled > 0 && st.filled < W1 * H1,
                 "face_keep halves the fill", &fails);
        uint8_t *dcls = NULL;
        int dd = 0, dh = 0, dw = 0;
        if (TiffIO_load(arena,
                        "output/_selftest_scroll_unroll/own_diagclass.tif",
                        &dcls, &dd, &dh, &dw) == 0) {
            size_t n6 = 0;
            for (long i = 0; i < (long)dw * dh; i++)
                if (dcls[i] == 6) n6++;
            sr_check(n6 > 0, "dropped-face sole cover -> class 6", &fails);
        } else {
            sr_check(0, "own_diagclass load", &fails);
        }
    }

    /* face gates: stretch one triangle's uv far beyond ratio*3d -> skip_uv */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        float verts[3 * 3] = { 10, 20, 10, 12, 20, 10, 10, 20, 12 };
        float uv[3 * 2] = { 0, 0, 300, 0, 0, 2 };   /* 2-vox edge -> 300 px */
        float normals[3 * 3] = { 0, 1, 0, 0, 1, 0, 0, 1, 0 };
        int32_t faces[3] = { 0, 1, 2 };
        int32_t fc[1] = { 0 };
        ps.verts = verts; ps.uv = uv; ps.normals = normals;
        ps.faces = faces; ps.face_cube = fc;
        ps.nv = 3; ps.nf = 1; ps.n_cubes = 1;
        ps.u_min = 0; ps.u_max = 300; ps.v_min = 0; ps.v_max = 2;
        RasterOpts o;
        RasterOpts_default(&o);
        o.band_cols = 0;
        RasterStats st;
        int rc = ScrollRaster_run(arena, &ps,
                                  "output/_selftest_scroll_unroll/raw", 64,
                                  "output/_selftest_scroll_unroll/gate", &o,
                                  &st);
        sr_check(rc == 0 && st.skip_uv == 1 && st.filled == 0,
                 "stretch gate skips the streak face", &fails);
    }

    /* Edge-scale regression: a coarse real face is valid at the default,
     * while an equally large synthetic face is still treated as a membrane.
     * An explicit real cap restores the old opt-in behaviour. */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        float verts[6 * 3] = {
            10, 20, 10,  30, 20, 10,  10, 20, 30,
            10, 40, 35,  30, 40, 35,  10, 40, 55
        };
        float uv[6 * 2] = {
             0, 10,   0, 30,  20, 10,
            25, 10,  25, 30,  45, 10
        };
        float normals[6 * 3] = {
            0, 1, 0,  0, 1, 0,  0, 1, 0,
            0, 1, 0,  0, 1, 0,  0, 1, 0
        };
        int32_t faces[2 * 3] = { 0, 1, 2, 3, 4, 5 };
        int32_t fc[2] = { 0, 0 };
        ps.verts = verts; ps.uv = uv; ps.normals = normals;
        ps.faces = faces; ps.face_cube = fc;
        ps.nv = 6; ps.nf = 2; ps.n_cubes = 1;
        ps.u_min = 0; ps.u_max = 45; ps.v_min = 10; ps.v_max = 30;

        RasterOpts o;
        RasterOpts_default(&o);
        o.band_cols = 0;
        o.pct_lo = 0.0; o.pct_hi = 100.0;
        o.synth_f0 = 1;
        RasterStats st;
        int rc = ScrollRaster_run(arena, &ps,
                                  "output/_selftest_scroll_unroll/raw", 64,
                                  "output/_selftest_scroll_unroll/edgeclass",
                                  &o, &st);
        sr_check(rc == 0 && st.skip_3d == 0 && st.filled > 0
                 && st.filled_synth == 0,
                 "edge gate: coarse real kept, coarse synth rejected", &fails);

        o.max_edge3d = 6.0;
        rc = ScrollRaster_run(arena, &ps,
                              "output/_selftest_scroll_unroll/raw", 64,
                              "output/_selftest_scroll_unroll/edgeclass_cap",
                              &o, &st);
        sr_check(rc == 0 && st.skip_3d == 1 && st.filled == 0,
                 "edge gate: explicit real cap remains available", &fails);
    }

    /* synth_f0: 2 real faces + 2 synthetic faces whose quad overlaps the real
     * one on px 15..20 -- synth paints class 7 only where the real quad
     * didn't, real px stay untouched, banded == single-shot byte-identical */
    {
        uint8_t *simg[2] = { NULL, NULL };
        uint8_t *scls[2] = { NULL, NULL };
        size_t sW[2] = { 0, 0 }, sH[2] = { 0, 0 };
        RasterStats st[2];
        size_t bands[2] = { 0, 7 };
        int rc_all = 0;
        for (int r = 0; r < 2 && rc_all == 0; r++) {
            PieceSet ps;
            memset(&ps, 0, sizeof(ps));
            /* real quad at y=20 (u px 0..20), synth quad at y=40 (brighter
             * in the gradient volume) spanning u px 15..35 */
            float verts[8 * 3] = {
                10, 20, 10,   30, 20, 10,   10, 20, 30,   30, 20, 30,
                10, 40, 35,   30, 40, 35,   10, 40, 55,   30, 40, 55
            };
            float uv[8 * 2] = {
                -1000, 10,  -1000, 30,  -980, 10,  -980, 30,
                -985, 10,   -985, 30,   -965, 10,  -965, 30
            };
            float normals[8 * 3] = {
                0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,
                0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0
            };
            int32_t faces[4 * 3] = { 0, 1, 2, 1, 3, 2, 4, 5, 6, 5, 7, 6 };
            int32_t fc[4] = { 0, 0, 0, 0 };
            ps.verts = verts; ps.uv = uv; ps.normals = normals;
            ps.faces = faces; ps.face_cube = fc;
            ps.nv = 8; ps.nf = 4; ps.n_cubes = 1;
            ps.u_min = -1000; ps.u_max = -965;
            ps.v_min = 10; ps.v_max = 30;
            RasterOpts o;
            RasterOpts_default(&o);
            o.band_cols = bands[r];
            o.pct_lo = 0.0; o.pct_hi = 100.0;
            o.max_edge3d = 100.0;
            o.synth_max_edge3d = 100.0;
            o.synth_f0 = 2;
            char pre[1024];
            snprintf(pre, sizeof(pre),
                     "output/_selftest_scroll_unroll/synth%d", r);
            if (ScrollRaster_run(arena, &ps,
                                 "output/_selftest_scroll_unroll/raw", 64,
                                 pre, &o, &st[r]) != 0) { rc_all = -1; break; }
            char path[1200];
            int d = 0, h = 0, w = 0;
            snprintf(path, sizeof(path), "%s_rawtex.tif", pre);
            if (TiffIO_load(arena, path, &simg[r], &d, &h, &w) != 0) {
                rc_all = -1;
                break;
            }
            sW[r] = (size_t)w;
            sH[r] = (size_t)h;
            snprintf(path, sizeof(path), "%s_diagclass.tif", pre);
            if (TiffIO_load(arena, path, &scls[r], &d, &h, &w) != 0)
                rc_all = -1;
        }
        sr_check(rc_all == 0, "synth raster rc", &fails);
        if (rc_all == 0) {
            sr_check(sW[0] == 36 && sH[0] == 21, "synth canvas dims", &fails);
            sr_check(st[0].filled == 441, "synth: real fill unchanged",
                     &fails);
            sr_check(st[0].filled_synth == 315, "synth: 15x21 fill px",
                     &fails);
            sr_check(sW[0] == sW[1] && sH[0] == sH[1]
                     && memcmp(simg[0], simg[1], sW[0] * sH[0]) == 0
                     && memcmp(scls[0], scls[1], sW[0] * sH[0]) == 0
                     && st[0].filled_synth == st[1].filled_synth,
                     "synth banding byte-identity", &fails);
            size_t Wn = sW[0];
            sr_check(scls[0][2 * Wn + 2] == 1, "real px class 1", &fails);
            sr_check(scls[0][2 * Wn + 16] == 1, "overlap px stays real",
                     &fails);
            sr_check(scls[0][10 * Wn + 30] == 7, "synth px class 7", &fails);
            sr_check(simg[0][10 * Wn + 30] == 255, "synth px sampled bright",
                     &fails);
        }
    }

    /* synth geometry OUTSIDE the raw volume: sample fails but class 7 must
     * still mark the ribbon (gray stays 0) */
    {
        PieceSet ps;
        memset(&ps, 0, sizeof(ps));
        /* real quad in-volume (u px 0..20), synth quad at y=500 (no cube) */
        float verts[8 * 3] = {
            10, 20, 10,   30, 20, 10,   10, 20, 30,   30, 20, 30,
            10, 500, 10,  30, 500, 10,  10, 500, 30,  30, 500, 30
        };
        float uv[8 * 2] = {
            0, 10,   0, 30,   20, 10,   20, 30,
            25, 10,  25, 30,  35, 10,   35, 30
        };
        float nrm[8 * 3];
        for (int i = 0; i < 8; i++) {
            nrm[i * 3 + 0] = 0; nrm[i * 3 + 1] = 1; nrm[i * 3 + 2] = 0;
        }
        int32_t faces[4 * 3] = { 0, 1, 3,  0, 3, 2,  4, 5, 7,  4, 7, 6 };
        int32_t fc[4] = { 0, 0, 0, 0 };
        ps.verts = verts; ps.uv = uv; ps.normals = nrm;
        ps.faces = faces; ps.face_cube = fc;
        ps.nv = 8; ps.nf = 4;
        ps.u_min = 0; ps.u_max = 35; ps.v_min = 10; ps.v_max = 30;
        RasterStats st;
        RasterOpts o;
        RasterOpts_default(&o);
        o.pct_lo = 0.0; o.pct_hi = 100.0;
        o.max_edge3d = 100.0;
        o.synth_max_edge3d = 100.0;
        o.synth_f0 = 2;
        int rc = ScrollRaster_run(arena, &ps,
                                  "output/_selftest_scroll_unroll/raw", 64,
                                  "output/_selftest_scroll_unroll/synthns",
                                  &o, &st);
        sr_check(rc == 0, "nosample raster rc", &fails);
        sr_check(st.filled_synth == 0 && st.synth_nosample > 0,
                 "nosample: marked but not painted", &fails);
        uint8_t *nsimg = NULL, *nscls = NULL;
        int d = 0, h = 0, w = 0;
        if (TiffIO_load(arena, "output/_selftest_scroll_unroll/"
                        "synthns_rawtex.tif", &nsimg, &d, &h, &w) == 0
            && TiffIO_load(arena, "output/_selftest_scroll_unroll/"
                           "synthns_diagclass.tif", &nscls, &d, &h, &w) == 0) {
            sr_check(nscls[10 * (size_t)w + 30] == 7
                     && nsimg[10 * (size_t)w + 30] == 0,
                     "nosample px: class 7, gray 0", &fails);
        } else {
            sr_check(0, "nosample tif readback", &fails);
        }
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[scroll_raster selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
