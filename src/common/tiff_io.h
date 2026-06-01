#ifndef TIFF_IO_INCLUDED
#define TIFF_IO_INCLUDED

#include <stdint.h>
#include "arena.h"

/* Load multi-page TIFF as flat uint8 volume.
 * Arena-allocates *out_vol of size D*H*W.
 * Returns 0 on success, -1 on failure.
 * Page 0 = z=0, page 1 = z=1, etc. */
int TiffIO_load(Arena_T arena, const char *path,
                uint8_t **out_vol, int *out_D, int *out_H, int *out_W);

/* Save flat uint8 volume as multi-page TIFF.
 * Returns 0 on success, -1 on failure. */
int TiffIO_save(const char *path,
                const uint8_t *vol, int D, int H, int W);

#endif
