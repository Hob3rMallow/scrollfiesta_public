#ifndef PNG_WRITE_INCLUDED
#define PNG_WRITE_INCLUDED

#include <stdint.h>

/* ============================================================================
 * png_write.h -- minimal, dependency-free 8-bit RGB PNG writer.
 *
 * Uses uncompressed ("stored") DEFLATE blocks so it needs no zlib -- only a
 * hand-rolled CRC-32 and Adler-32. Intended for small generated textures (e.g.
 * the UV checkerboard), not photographic data. rgb is [w*h*3], row-major,
 * top-to-bottom, R,G,B per pixel. Returns 0 on success, -1 on failure.
 * ==========================================================================*/
int PngWrite_rgb(const char *path, int w, int h, const uint8_t *rgb);

#endif
