/*
 * tiff_io_stub.c — TiffIO_* stubs for TIFF-less library builds.
 *
 * Compiled instead of tiff_io.c when SCROLLFIESTA_WITH_TIFF=OFF (the
 * embedded/library configuration, where volumes arrive in memory through
 * the C API rather than from prediction TIFFs on disk). Every entry point
 * fails cleanly with -1, which the callers (halo_loader, mesh_extract's
 * file-path mode, tifxyz writers) already treat as "input unavailable".
 */
#include "tiff_io.h"

#include <stdio.h>

static int tiff_stub_fail(const char *fn, const char *path)
{
    fprintf(stderr, "%s: ScrollFiesta was built without TIFF support "
                    "(SCROLLFIESTA_WITH_TIFF=OFF); cannot access '%s'\n",
            fn, path ? path : "(null)");
    return -1;
}

int TiffIO_load(Arena_T arena, const char *path,
                uint8_t **out_vol, int *out_D, int *out_H, int *out_W)
{
    (void)arena;
    if (out_vol) *out_vol = NULL;
    if (out_D) *out_D = 0;
    if (out_H) *out_H = 0;
    if (out_W) *out_W = 0;
    return tiff_stub_fail("TiffIO_load", path);
}

int TiffIO_save(const char *path,
                const uint8_t *vol, int D, int H, int W)
{
    (void)vol; (void)D; (void)H; (void)W;
    return tiff_stub_fail("TiffIO_save", path);
}

int TiffIO_save_float2d(const char *path, const float *img, int W, int H)
{
    (void)img; (void)W; (void)H;
    return tiff_stub_fail("TiffIO_save_float2d", path);
}

int TiffIO_load_float2d(Arena_T arena, const char *path,
                        float **out_img, int *out_W, int *out_H)
{
    (void)arena;
    if (out_img) *out_img = NULL;
    if (out_W) *out_W = 0;
    if (out_H) *out_H = 0;
    return tiff_stub_fail("TiffIO_load_float2d", path);
}
