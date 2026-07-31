/* ves_png.c -- thin wrappers over the vendored stb_image_write single-header.
 * This is the ONE translation unit that instantiates the implementation. */

#include "ves_png.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int VesPng_write_gray(const char *path, const uint8_t *data, int w, int h)
{
    if (path == NULL || data == NULL || w <= 0 || h <= 0) return -1;
    return stbi_write_png(path, w, h, 1, data, w) ? 0 : -1;
}

int VesPng_write_rgb(const char *path, const uint8_t *data, int w, int h)
{
    if (path == NULL || data == NULL || w <= 0 || h <= 0) return -1;
    return stbi_write_png(path, w, h, 3, data, w * 3) ? 0 : -1;
}
