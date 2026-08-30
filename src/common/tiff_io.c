#include "ves_platform.h"
#include "tiff_io.h"
#include <tiff.h>
#include <tiffio.h>
#include <assert.h>
#include <stdio.h>
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

int TiffIO_save_float2d(const char *path, const float *img, int W, int H)
{
    assert(path);
    assert(W >= 0 && H >= 0);

    if (W == 0 || H == 0) {
        return -1;
    }

    assert(img);

    ves_mutex_lock(&tiff_lock);

    TIFF *tif = TIFFOpen(path, "w");
    if (tif == NULL) {
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32_t)W);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32_t)H);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, (uint32_t)H);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

    for (int y = 0; y < H; y++) {
        /* TIFFWriteScanline takes a non-const buffer; the data is not modified. */
        if (TIFFWriteScanline(tif, (void *)(img + (size_t)y * (size_t)W),
                              (uint32_t)y, 0) < 0) {
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }
    }

    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);
    return 0;
}

/* 64-bit seek: MSVC has no fseeko. */
static int rows_seek(FILE *f, uint64_t off)
{
#if defined(_WIN32)
    return _fseeki64(f, (long long)off, SEEK_SET);
#else
    return fseeko(f, (off_t)off, SEEK_SET);
#endif
}

int TiffIO_rows_open(Arena_T arena, const char *path, TiffRowReader *rd)
{
    assert(arena);
    assert(path);
    assert(rd);

    memset(rd, 0, sizeof(*rd));

    ves_mutex_lock(&tiff_lock);
    TIFF *tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    uint32_t w = 0, h = 0, rps = 0;
    uint16_t bps = 0, spp = 0, comp = 0, fmt = 0;
    uint64_t *offs = NULL;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &comp);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &fmt);
    TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rps);

    if (w == 0 || h == 0 || bps != 8 || spp != 1
        || comp != COMPRESSION_NONE
        || (fmt != SAMPLEFORMAT_UINT && fmt != 0)
        || TIFFIsTiled(tif)
        || !TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offs) || offs == NULL) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    if (rps == 0 || rps > h) rps = h;   /* single-strip default */

    uint32_t n_strips = (uint32_t)TIFFNumberOfStrips(tif);
    if (n_strips == 0 || (uint64_t)n_strips * rps < h) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    uint64_t *off_copy = (uint64_t *)ARENA_ALLOC(arena,
                             (size_t)n_strips * sizeof(uint64_t));
    for (uint32_t s = 0; s < n_strips; s++)
        off_copy[s] = offs[s];

    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);

    FILE *f = fopen(path, "rb");
    if (f == NULL)
        return -1;

    rd->f = f;
    rd->W = (int)w;
    rd->H = (int)h;
    rd->rows_per_strip = rps;
    rd->strip_off = off_copy;
    rd->n_strips = n_strips;
    return 0;
}

int TiffIO_rows_read(TiffRowReader *rd, int y, int x0, int w, uint8_t *dst)
{
    assert(rd);
    assert(dst);

    if (rd->f == NULL || y < 0 || y >= rd->H || x0 < 0 || w < 0
        || x0 + w > rd->W)
        return -1;
    if (w == 0)
        return 0;

    uint32_t strip = (uint32_t)y / rd->rows_per_strip;
    uint32_t row_in = (uint32_t)y % rd->rows_per_strip;
    if (strip >= rd->n_strips)
        return -1;

    uint64_t off = rd->strip_off[strip]
                 + (uint64_t)row_in * (uint64_t)rd->W + (uint64_t)x0;
    if (rows_seek(rd->f, off) != 0)
        return -1;
    if (fread(dst, 1, (size_t)w, rd->f) != (size_t)w)
        return -1;
    return 0;
}

void TiffIO_rows_close(TiffRowReader *rd)
{
    if (rd == NULL)
        return;
    if (rd->f != NULL) {
        fclose(rd->f);
        rd->f = NULL;
    }
    rd->strip_off = NULL;   /* arena-owned */
    rd->n_strips = 0;
}

int TiffIO_load_float2d(Arena_T arena, const char *path,
                        float **out_img, int *out_W, int *out_H)
{
    assert(arena);
    assert(path);
    assert(out_img && out_W && out_H);

    ves_mutex_lock(&tiff_lock);

    TIFF *tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    uint32_t w = 0, h = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    if (w == 0 || h == 0) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }

    float *img = (float *)ARENA_ALLOC(arena,
                    (long)((size_t)w * (size_t)h * sizeof(float)));

    for (uint32_t y = 0; y < h; y++) {
        if (TIFFReadScanline(tif, img + (size_t)y * (size_t)w, y, 0) < 0) {
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }
    }

    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);

    *out_img = img;
    *out_W = (int)w;
    *out_H = (int)h;
    return 0;
}

int TiffIO_load_int32_2d(Arena_T arena, const char *path,
                         int32_t **out_img, int *out_W, int *out_H)
{
    assert(arena);
    assert(path);
    assert(out_img && out_W && out_H);

    ves_mutex_lock(&tiff_lock);
    TIFF *tif = TIFFOpen(path, "r");
    if (tif == NULL) {
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    uint32_t w = 0, h = 0;
    uint16_t bps = 0, spp = 0, format = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &format);
    if (w == 0 || h == 0 || bps != 32 || spp != 1 ||
        format != SAMPLEFORMAT_INT) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    int32_t *img = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)w * (size_t)h * sizeof *img));
    for (uint32_t y = 0; y < h; y++) {
        if (TIFFReadScanline(tif, img + (size_t)y * (size_t)w, y, 0) < 0) {
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }
    }
    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);
    *out_img = img;
    *out_W = (int)w;
    *out_H = (int)h;
    return 0;
}
