/* ves_png.c -- thin wrappers over the vendored stb_image_write single-header.
 * This is the ONE translation unit that instantiates the implementation. */

#include "ves_png.h"
#include <limits.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/* stb uses signed int for filtered bytes, row offsets, entropy estimates,
 * and its doubling output buffer. Oversized previews must not reach it;
 * full-resolution TIFFs and the native tile pyramid remain available. */
static int ves_png_size_ok(int w, int h, int channels)
{
    uint64_t row;
    if (w <= 0 || h <= 0) return 0;
    row = (uint64_t)w * channels;
    return row <= INT_MAX / 128 && (row + 1) * (uint64_t)h <= INT_MAX / 2;
}

int VesPng_write_gray(const char *path, const uint8_t *data, int w, int h)
{
    if (path == NULL || data == NULL || !ves_png_size_ok(w, h, 1)) return -1;
    return stbi_write_png(path, w, h, 1, data, w) ? 0 : -1;
}

int VesPng_write_rgb(const char *path, const uint8_t *data, int w, int h)
{
    if (path == NULL || data == NULL || !ves_png_size_ok(w, h, 3)) return -1;
    return stbi_write_png(path, w, h, 3, data, w * 3) ? 0 : -1;
}
