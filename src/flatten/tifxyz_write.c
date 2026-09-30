#define _USE_MATH_DEFINES
#include "tifxyz_write.h"

#include "../common/ves_platform.h"
#include "../common/tiff_io.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TIFXYZ_INVALID  (-1.0f)

static double dmin3(double a, double b, double c)
{
    double m = a; if (b < m) m = b; if (c < m) m = c; return m;
}
static double dmax3(double a, double b, double c)
{
    double m = a; if (b > m) m = b; if (c > m) m = c; return m;
}

int TifXYZ_write_double(Arena_T arena, const char *out_dir, const char *uuid,
                        const float *verts, size_t nv,
                        const int32_t *faces, size_t nf,
                        const double *uv, const TifXYZOpts *opts,
                        TifXYZStats *out_stats)
{
    assert(arena);
    assert(out_dir && uuid);

    if (nv < 3 || nf < 1 || verts == NULL || faces == NULL || uv == NULL) {
        return -1;
    }

    float px_per_vox = (opts && opts->px_per_vox > 0.0f)
                           ? opts->px_per_vox : FLATTEN_PX_PER_VOX;
    int   max_dim    = (opts && opts->max_dim > 0)
                           ? opts->max_dim : FLATTEN_MAX_GRID_DIM;
    double contested_tol = (opts && opts->contested_tol_vox > 0.0f)
                               ? (double)opts->contested_tol_vox : 1.5;
    double contested_tol2 = contested_tol * contested_tol;

    /* --- UV bounds + grid sizing. --- */
    double umin = INFINITY, umax = -INFINITY, vmin = INFINITY, vmax = -INFINITY;
    for (size_t i = 0; i < nv; i++) {
        double u = (double)uv[i * 2 + 0];
        double v = (double)uv[i * 2 + 1];
        if (u < umin) umin = u;
        if (u > umax) umax = u;
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    double uspan = umax - umin; if (uspan < 1e-9) uspan = 1e-9;
    double vspan = vmax - vmin; if (vspan < 1e-9) vspan = 1e-9;

    /* Grid sized ~1 cell per voxel (vc_obj2tifxyz metric mode: ceil(range)+1). */
    long Wl = (long)ceil(uspan * (double)px_per_vox) + 1;
    long Hl = (long)ceil(vspan * (double)px_per_vox) + 1;
    int W = (Wl < 2) ? 2 : (Wl > max_dim ? max_dim : (int)Wl);
    int H = (Hl < 2) ? 2 : (Hl > max_dim ? max_dim : (int)Hl);

    /* meta.json "scale" = world units per pixel (matches vc_obj2tifxyz). */
    double scale_u = uspan / (double)(W - 1);
    double scale_v = vspan / (double)(H - 1);
    /* UV -> pixel factor (distinct from the metric scale above). */
    double pxf_u = (double)(W - 1) / uspan;
    double pxf_v = (double)(H - 1) / vspan;

    size_t ncell = (size_t)W * (size_t)H;
    float *xg = (float *)ARENA_ALLOC(arena, (size_t)(ncell * sizeof(float)));
    float *yg = (float *)ARENA_ALLOC(arena, (size_t)(ncell * sizeof(float)));
    float *zg = (float *)ARENA_ALLOC(arena, (size_t)(ncell * sizeof(float)));
    uint8_t *state = (uint8_t *)ARENA_CALLOC(arena, (size_t)ncell, 1L);
    for (size_t i = 0; i < ncell; i++) {
        xg[i] = TIFXYZ_INVALID; yg[i] = TIFXYZ_INVALID; zg[i] = TIFXYZ_INVALID;
    }

    /* --- Rasterize triangles in UV space (barycentric, contested-safe). --- */
    double lo[3] = { INFINITY, INFINITY, INFINITY };
    double hi[3] = { -INFINITY, -INFINITY, -INFINITY };
    size_t filled = 0, contested = 0;
    const double eps = 1e-6;

    for (size_t fcur = 0; fcur < nf; fcur++) {
        int32_t i0 = faces[fcur * 3 + 0];
        int32_t i1 = faces[fcur * 3 + 1];
        int32_t i2 = faces[fcur * 3 + 2];
        if (i0 < 0 || i1 < 0 || i2 < 0 ||
            (size_t)i0 >= nv || (size_t)i1 >= nv || (size_t)i2 >= nv) {
            continue;
        }

        /* Pixel-space triangle (col=fu, row=fv). */
        double fu0 = ((double)uv[i0 * 2 + 0] - umin) * pxf_u;
        double fv0 = ((double)uv[i0 * 2 + 1] - vmin) * pxf_v;
        double fu1 = ((double)uv[i1 * 2 + 0] - umin) * pxf_u;
        double fv1 = ((double)uv[i1 * 2 + 1] - vmin) * pxf_v;
        double fu2 = ((double)uv[i2 * 2 + 0] - umin) * pxf_u;
        double fv2 = ((double)uv[i2 * 2 + 1] - vmin) * pxf_v;

        double d = (fv1 - fv2) * (fu0 - fu2) + (fu2 - fu1) * (fv0 - fv2);
        if (fabs(d) < 1e-9) continue;  /* degenerate in UV */

        /* World coords, reordered (z,y,x) -> (x,y,z). */
        double wx0 = (double)verts[i0 * 3 + 2], wy0 = (double)verts[i0 * 3 + 1],
               wz0 = (double)verts[i0 * 3 + 0];
        double wx1 = (double)verts[i1 * 3 + 2], wy1 = (double)verts[i1 * 3 + 1],
               wz1 = (double)verts[i1 * 3 + 0];
        double wx2 = (double)verts[i2 * 3 + 2], wy2 = (double)verts[i2 * 3 + 1],
               wz2 = (double)verts[i2 * 3 + 0];

        /* Pad the pixel bbox by 1 (matches vc_obj2tifxyz: catches edge slivers). */
        int minx = (int)floor(dmin3(fu0, fu1, fu2)) - 1;
        int maxx = (int)ceil(dmax3(fu0, fu1, fu2)) + 1;
        int miny = (int)floor(dmin3(fv0, fv1, fv2)) - 1;
        int maxy = (int)ceil(dmax3(fv0, fv1, fv2)) + 1;
        if (minx < 0) minx = 0;
        if (miny < 0) miny = 0;
        if (maxx > W - 1) maxx = W - 1;
        if (maxy > H - 1) maxy = H - 1;

        for (int y = miny; y <= maxy; y++) {
            for (int x = minx; x <= maxx; x++) {
                double px = (double)x, py = (double)y;
                double l0 = ((fv1 - fv2) * (px - fu2) + (fu2 - fu1) * (py - fv2)) / d;
                double l1 = ((fv2 - fv0) * (px - fu2) + (fu0 - fu2) * (py - fv2)) / d;
                double l2 = 1.0 - l0 - l1;
                if (l0 < -eps || l1 < -eps || l2 < -eps) continue;

                size_t cell = (size_t)y * (size_t)W + (size_t)x;
                double wx = l0 * wx0 + l1 * wx1 + l2 * wx2;
                double wy = l0 * wy0 + l1 * wy1 + l2 * wy2;
                double wz = l0 * wz0 + l1 * wz1 + l2 * wz2;
                if (state[cell] == 2) continue; /* already known contested */
                if (state[cell] == 1) {
                    double dx = wx - (double)xg[cell];
                    double dy = wy - (double)yg[cell];
                    double dz = wz - (double)zg[cell];
                    if (dx*dx + dy*dy + dz*dz > contested_tol2) {
                        /* Two distinct sheets landed on one UV pixel. Keeping
                         * either one would manufacture connectivity in Villa. */
                        xg[cell] = yg[cell] = zg[cell] = TIFXYZ_INVALID;
                        state[cell] = 2;
                        contested++;
                        if (filled > 0) filled--;
                    }
                    continue;
                }
                xg[cell] = (float)wx;
                yg[cell] = (float)wy;
                zg[cell] = (float)wz;
                state[cell] = 1;
                filled++;
                if (wx < lo[0]) lo[0] = wx; if (wx > hi[0]) hi[0] = wx;
                if (wy < lo[1]) lo[1] = wy; if (wy > hi[1]) hi[1] = wy;
                if (wz < lo[2]) lo[2] = wz; if (wz > hi[2]) hi[2] = wz;
            }
        }
    }

    /* Recompute the bbox from surviving cells. A point first written and later
     * invalidated as contested must not expand the advertised surface bounds. */
    lo[0] = lo[1] = lo[2] = INFINITY;
    hi[0] = hi[1] = hi[2] = -INFINITY;
    for (size_t i = 0; i < ncell; i++) {
        if (state[i] != 1) continue;
        if (xg[i] < lo[0]) lo[0] = xg[i]; if (xg[i] > hi[0]) hi[0] = xg[i];
        if (yg[i] < lo[1]) lo[1] = yg[i]; if (yg[i] > hi[1]) hi[1] = yg[i];
        if (zg[i] < lo[2]) lo[2] = zg[i]; if (zg[i] > hi[2]) hi[2] = zg[i];
    }
    if (filled == 0) {
        lo[0] = lo[1] = lo[2] = 0.0;
        hi[0] = hi[1] = hi[2] = 0.0;
    }

    /* --- Write x/y/z/mask.tif + meta.json. --- */
    char path[4096];
    snprintf(path, sizeof path, "%s/x.tif", out_dir);
    ves_ensure_parent_dir(path);   /* creates out_dir + parents */

    int rc = 0;
    rc |= TiffIO_save_float2d(path, xg, W, H);
    snprintf(path, sizeof path, "%s/y.tif", out_dir);
    rc |= TiffIO_save_float2d(path, yg, W, H);
    snprintf(path, sizeof path, "%s/z.tif", out_dir);
    rc |= TiffIO_save_float2d(path, zg, W, H);
    for (size_t i = 0; i < ncell; i++)
        state[i] = (state[i] == 1) ? (uint8_t)255 : (uint8_t)0;
    snprintf(path, sizeof path, "%s/mask.tif", out_dir);
    rc |= TiffIO_save(path, state, 1, H, W);
    if (rc != 0) return -1;

    snprintf(path, sizeof path, "%s/meta.json", out_dir);
    FILE *m = fopen(path, "w");
    if (m == NULL) return -1;
    fprintf(m,
        "{\n"
        "  \"bbox\": [[%.6f, %.6f, %.6f], [%.6f, %.6f, %.6f]],\n"
        "  \"type\": \"seg\",\n"
        "  \"uuid\": \"%s\",\n"
        "  \"format\": \"tifxyz\",\n"
        "  \"scale\": [%.8f, %.8f],\n"
        "  \"contested_tolerance_vox\": %.4f,\n"
        "  \"invalid_contested_cells\": %zu\n"
        "}\n",
        lo[0], lo[1], lo[2], hi[0], hi[1], hi[2],
        uuid, scale_u, scale_v, contested_tol, contested);
    fclose(m);

    if (out_stats) {
        out_stats->width = W;
        out_stats->height = H;
        out_stats->scale_u = scale_u;
        out_stats->scale_v = scale_v;
        out_stats->cells_filled = filled;
        out_stats->cells_contested = contested;
        for (int i = 0; i < 3; i++) {
            out_stats->bbox_lo[i] = lo[i];
            out_stats->bbox_hi[i] = hi[i];
        }
    }
    return 0;
}

int TifXYZ_write(Arena_T arena, const char *out_dir, const char *uuid,
                 const float *verts, size_t nv,
                 const int32_t *faces, size_t nf,
                 const float *uv, const TifXYZOpts *opts,
                 TifXYZStats *out_stats)
{
    Arena_Mark mark;
    double *uv_double = NULL;
    int result = -1;
    size_t i = 0;
    if (arena == NULL || (uv == NULL && nv != 0)) return -1;
    mark = Arena_save(arena);
    uv_double = (double *)ARENA_ALLOC(arena, nv * 2 * sizeof(*uv_double));
    for (i = 0; i < nv * 2; i++) uv_double[i] = (double)uv[i];
    result = TifXYZ_write_double(arena, out_dir, uuid, verts, nv, faces, nf,
                                 uv_double, opts, out_stats);
    Arena_restore(arena, mark);
    return result;
}

/* ============================================================================
 * Self-test.
 * ==========================================================================*/

int TifXYZ_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* (1) Float-TIFF round-trip (incl. negative sentinel values). */
    {
        const int W = 7, H = 5;
        float img[35];
        for (int i = 0; i < W * H; i++) {
            img[i] = (i % 3 == 0) ? -1.0f : (float)i * 0.5f - 2.25f;
        }
        const char *p = "output/flatten_selftest/float2d.tif";
        ves_ensure_parent_dir(p);
        int wrc = TiffIO_save_float2d(p, img, W, H);
        float *back = NULL; int bw = 0, bh = 0;
        int rrc = TiffIO_load_float2d(arena, p, &back, &bw, &bh);
        int ok = (wrc == 0 && rrc == 0 && bw == W && bh == H && back != NULL);
        if (ok) {
            for (int i = 0; i < W * H; i++) {
                if (back[i] != img[i]) { ok = 0; break; }
            }
        }
        if (!ok) {
            fprintf(stderr, "[tifxyz selftest] FLOAT-TIFF FAIL "
                    "wrc=%d rrc=%d bw=%d bh=%d\n", wrc, rrc, bw, bh);
            fails++;
        } else {
            fprintf(stderr, "[tifxyz selftest] float-tiff round-trip OK\n");
        }
    }

    /* (2) Rasterize one triangle; verify reorder + barycentric + sentinel. */
    {
        /* verts in (z,y,x); uv chosen so the triangle is the lower-left half. */
        float v[9] = {
            10.0f, 20.0f, 30.0f,   /* uv (0,0) */
            10.0f, 20.0f, 40.0f,   /* uv (1,0) -> world x grows with u */
            60.0f, 20.0f, 30.0f    /* uv (0,1) -> world z grows with v */
        };
        int32_t f[3] = { 0, 1, 2 };
        float uv[6] = { 0.0f, 0.0f,  1.0f, 0.0f,  0.0f, 1.0f };

        TifXYZOpts opt = { 32.0f, 0, 1.5f };  /* 32 px/unit -> ~33x33 grid for span 1 */
        TifXYZStats st = { 0 };
        const char *dir = "output/flatten_selftest/tifxyz";
        int wrc = TifXYZ_write(arena, dir, "selftest", v, 3, f, 1, uv, &opt, &st);

        float *xg = NULL, *yg = NULL, *zg = NULL;
        uint8_t *mask = NULL;
        int gw = 0, gh = 0, gw2 = 0, gh2 = 0, gw3 = 0, gh3 = 0;
        char p[4096];
        snprintf(p, sizeof p, "%s/x.tif", dir);
        int r1 = TiffIO_load_float2d(arena, p, &xg, &gw, &gh);
        snprintf(p, sizeof p, "%s/y.tif", dir);
        int r2 = TiffIO_load_float2d(arena, p, &yg, &gw2, &gh2);
        snprintf(p, sizeof p, "%s/z.tif", dir);
        int r3 = TiffIO_load_float2d(arena, p, &zg, &gw3, &gh3);
        int md = 0, mh = 0, mw = 0;
        snprintf(p, sizeof p, "%s/mask.tif", dir);
        int rm = TiffIO_load(arena, p, &mask, &md, &mh, &mw);

        int ok = (wrc == 0 && r1 == 0 && r2 == 0 && r3 == 0 && rm == 0 &&
                  gw == st.width && gh == st.height && st.cells_filled > 0);

        /* Inside pixel near uv (0.25, 0.25): u+v=0.5 < 1 -> inside. */
        int ix = (int)(0.25 * (gw - 1));
        int iy = (int)(0.25 * (gh - 1));
        size_t ci = (size_t)iy * (size_t)gw + (size_t)ix;
        if (ok) {
            ok = ok && (xg[ci] > 29.0f && xg[ci] < 41.0f);   /* world x 30..40 */
            ok = ok && (yg[ci] > 19.0f && yg[ci] < 21.0f);   /* world y == 20 */
            ok = ok && (zg[ci] > 9.0f && zg[ci] < 61.0f);    /* world z 10..60 */
            ok = ok && mask[ci] == 255;
        }
        /* Outside pixel near uv (0.9, 0.9): u+v=1.8 > 1 -> unmapped (-1). */
        int ox = (int)(0.9 * (gw - 1));
        int oy = (int)(0.9 * (gh - 1));
        size_t co = (size_t)oy * (size_t)gw + (size_t)ox;
        if (ok) {
            ok = ok && (xg[co] == TIFXYZ_INVALID);
            ok = ok && mask[co] == 0;
        }
        if (!ok) {
            double in_x = (xg != NULL) ? (double)xg[ci] : 0.0;
            double in_y = (yg != NULL) ? (double)yg[ci] : 0.0;
            double in_z = (zg != NULL) ? (double)zg[ci] : 0.0;
            double out_x = (xg != NULL) ? (double)xg[co] : 0.0;
            fprintf(stderr, "[tifxyz selftest] RASTER FAIL wrc=%d filled=%zu "
                    "grid=%dx%d in(x=%.2f y=%.2f z=%.2f) out_x=%.2f\n",
                    wrc, st.cells_filled, gw, gh, in_x, in_y, in_z, out_x);
            fails++;
        } else {
            fprintf(stderr, "[tifxyz selftest] raster + reorder OK "
                    "(grid %dx%d, filled=%zu)\n", gw, gh, st.cells_filled);
        }
    }

    /* (3) Two distinct sheets with identical UV must invalidate every shared
     * cell instead of exposing an arbitrary first-triangle winner. */
    {
        float v[18] = {
            10.0f, 20.0f, 30.0f,  10.0f, 20.0f, 40.0f,
            60.0f, 20.0f, 30.0f,
            20.0f, 20.0f, 30.0f,  20.0f, 20.0f, 40.0f,
            70.0f, 20.0f, 30.0f
        };
        int32_t f[6] = { 0, 1, 2,  3, 4, 5 };
        float uv[12] = {
            0.0f, 0.0f,  1.0f, 0.0f,  0.0f, 1.0f,
            0.0f, 0.0f,  1.0f, 0.0f,  0.0f, 1.0f
        };
        TifXYZOpts opt = { 16.0f, 0, 1.5f };
        TifXYZStats st = { 0 };
        const char *dir = "output/flatten_selftest/tifxyz_contested";
        int wrc = TifXYZ_write(arena, dir, "contested", v, 6, f, 2,
                               uv, &opt, &st);
        float *xg = NULL; uint8_t *mask = NULL;
        int gw = 0, gh = 0, md = 0, mh = 0, mw = 0;
        char p[4096];
        snprintf(p, sizeof p, "%s/x.tif", dir);
        int rx = TiffIO_load_float2d(arena, p, &xg, &gw, &gh);
        snprintf(p, sizeof p, "%s/mask.tif", dir);
        int rm = TiffIO_load(arena, p, &mask, &md, &mh, &mw);
        int ok = (wrc == 0 && rx == 0 && rm == 0 && st.cells_contested > 0
                  && st.cells_filled == 0 && gw == mw && gh == mh);
        if (ok) {
            for (int i = 0; i < gw * gh; i++) {
                if (mask[i] != 0 || xg[i] != TIFXYZ_INVALID) {
                    ok = 0; break;
                }
            }
        }
        if (!ok) {
            fprintf(stderr, "[tifxyz selftest] CONTESTED FAIL wrc=%d "
                    "filled=%zu contested=%zu grid=%dx%d\n",
                    wrc, st.cells_filled, st.cells_contested, gw, gh);
            fails++;
        } else {
            fprintf(stderr, "[tifxyz selftest] contested overlap invalidated "
                    "(%zu cells)\n", st.cells_contested);
        }
    }

    Arena_dispose(&arena);
    return fails;
}
