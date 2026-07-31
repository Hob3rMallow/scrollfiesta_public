#ifndef VES_PNG_INCLUDED
#define VES_PNG_INCLUDED

#include <stdint.h>

/* Direct PNG writers (vendored stb_image_write, deps/include/stb_image_write.h)
 * so the pipeline never round-trips images through GDI+/PowerShell -- GDI+
 * silently truncates the very wide (>~3k-col) diagnostic/texture rasters on
 * read. All are row-major, top-left origin. Return 0 on success, -1 on failure. */

/* 8-bit grayscale, 1 byte/pixel, w*h bytes. */
int VesPng_write_gray(const char *path, const uint8_t *data, int w, int h);

/* 24-bit RGB, 3 bytes/pixel in R,G,B order, w*h*3 bytes. */
int VesPng_write_rgb(const char *path, const uint8_t *data, int w, int h);

#endif
