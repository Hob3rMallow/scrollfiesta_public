#include "ves_platform.h"
#include "tiff_io.h"
#include <tiff.h>
#include <tiffio.h>
#include <assert.h>
#include <string.h>

/* Mutex-serialize all libtiff calls [common.md §9, c-style-guide §11.2] */
static ves_mutex_t tiff_lock = VES_MUTEX_INITIALIZER;

int TiffIO_load(Arena_T arena, const char *path,
                uint8_t **out_vol, int *out_D, int *out_H, int *out_W)
{
    assert(arena);
    assert(path);
    assert(out_vol && out_D && out_H && out_W);

    ves_mutex_lock(&tiff_lock);

    TIFF *tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    /* Read dimensions from first page */
    uint32_t w = 0, h = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);

    /* Count pages */
    int n_pages = 0;
    do {
        n_pages++;
    } while (TIFFReadDirectory(tif));

    /* Allocate volume */
    size_t vol_size = (size_t)n_pages * (size_t)h * (size_t)w;
    uint8_t *vol = (uint8_t *)ARENA_CALLOC(arena, (long)vol_size, 1L);

    /* Re-read from beginning */
    TIFFSetDirectory(tif, 0);

    for (int z = 0; z < n_pages; z++) {
        /* Verify page matches first page dimensions */
        uint32_t pw = 0, ph = 0;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &pw);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ph);

        if (pw != w || ph != h) {
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }

        /* Read using strip-based API */
        tsize_t strip_size = TIFFStripSize(tif);
        tstrip_t n_strips = TIFFNumberOfStrips(tif);
        uint8_t *page_start = vol + (size_t)z * (size_t)h * (size_t)w;
        tsize_t offset = 0;

        for (tstrip_t s = 0; s < n_strips; s++) {
            tsize_t bytes_read = TIFFReadEncodedStrip(tif, s,
                                                       page_start + offset,
                                                       strip_size);
            if (bytes_read < 0) {
                TIFFClose(tif);
                ves_mutex_unlock(&tiff_lock);
                return -1;
            }
            offset += bytes_read;
        }

        /* Verify page ordering */
        assert((int)TIFFCurrentDirectory(tif) == z);

        if (z < n_pages - 1) {
            if (!TIFFReadDirectory(tif)) {
                TIFFClose(tif);
                ves_mutex_unlock(&tiff_lock);
                return -1;
            }
        }
    }

    /* Advise kernel to drop cached pages for this file — prevents page cache
     * pollution when processing many volumes sequentially. */
    ves_fadvise(TIFFFileno(tif), 0, 0, VES_FADVISE_DONTNEED);

    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);

    *out_vol = vol;
    *out_D = n_pages;
    *out_H = (int)h;
    *out_W = (int)w;
    return 0;
}

int TiffIO_save(const char *path,
                const uint8_t *vol, int D, int H, int W)
{
    assert(path);
    assert(D >= 0 && H >= 0 && W >= 0);

    if (D == 0 || H == 0 || W == 0) {
        return -1;
    }

    assert(vol);

    ves_mutex_lock(&tiff_lock);

    TIFF *tif = TIFFOpen(path, "w");
    if (tif == NULL) {
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    for (int z = 0; z < D; z++) {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32_t)W);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32_t)H);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, (uint32_t)H);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

        /* Multi-page TIFF subfile type */
        TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        TIFFSetField(tif, TIFFTAG_PAGENUMBER, (uint16_t)z, (uint16_t)D);

        const uint8_t *page = vol + (size_t)z * (size_t)H * (size_t)W;

        for (int y = 0; y < H; y++) {
            if (TIFFWriteScanline(tif, (void *)(page + y * W), (uint32_t)y, 0) < 0) {
                TIFFClose(tif);
                ves_mutex_unlock(&tiff_lock);
                return -1;
            }
        }

        TIFFWriteDirectory(tif);
    }

    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);
    return 0;
}
