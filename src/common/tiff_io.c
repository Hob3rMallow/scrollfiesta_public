#include "ves_platform.h"
#include "tiff_io.h"
#include "pipeline_constants.h"
#include <tiff.h>
#include <tiffio.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifdef VES_TIFF_PARALLEL_DEFLATE
#include <zlib.h>
#include "ves_omp.h"
#endif

/* Mutex-serialize all libtiff calls [common.md §9, c-style-guide §11.2] */
static ves_mutex_t tiff_lock = VES_MUTEX_INITIALIZER;

#include "tiff_float_window.inc"

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
    uint8_t *vol = (uint8_t *)ARENA_CALLOC(arena, vol_size, 1);

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

/* Leave space for directories and strip tables before classic TIFF's 4 GiB
 * offset limit. A long full-resolution ribbon can exceed it in one page. */
static const char *tiff_write_mode(size_t bytes, int pages)
{
    const uint64_t limit = UINT64_C(0xffffffff) - UINT64_C(1048576);
    return bytes >= limit || (uint64_t)pages * 4096 >= limit - bytes
         ? "w8" : "w";
}

static uint32_t tiff_strip_rows(int H, int W, int bytes_per_pixel)
{
    size_t rows = TIFF_OUTPUT_STRIP_BYTES / ((size_t)W * bytes_per_pixel);
    if (rows == 0) rows = 1;
    if (rows > (size_t)H) rows = (size_t)H;
    return (uint32_t)rows;
}

/* Compression of the writers: none (the historical files), or deflate strips (TiffIO_save_deflate,
 * TiffIO_save_float2d_deflate) at zlib level 1, no predictor -- the same strips and rows, so every strip-based reader
 * decodes them unchanged. */
static void tiff_set_compression(TIFF *tif, int deflate)
{
    TIFFSetField(tif, TIFFTAG_COMPRESSION, deflate ? COMPRESSION_ADOBE_DEFLATE : COMPRESSION_NONE);
    if (deflate)
        TIFFSetField(tif, TIFFTAG_ZIPQUALITY, TIFF_DEFLATE_LEVEL);
}

#ifdef VES_TIFF_PARALLEL_DEFLATE
typedef struct TiffDeflateSlot {
    z_stream stream;
    uint8_t *data;
    uint8_t *predicted;
    size_t bytes;
    int status;
} TiffDeflateSlot;

/* TIFF predictor 3, one sample/pixel, 32-bit IEEE input. It rearranges each
 * row into most-significant-byte planes, then differences the entire byte
 * row (including plane boundaries), exactly as libtiff's fpDiff specifies.
 * No floating-point arithmetic; the caller's pixels are never modified. */
static void tiff_predict_float_rows(uint8_t *out, const uint8_t *in, size_t rows, size_t width)
{
    const uint16_t one = 1;
    const int little = *(const uint8_t *)&one;
    const size_t bytes = width * 4;
    for (size_t y = 0; y < rows; y++) {
        uint8_t *dst = out + y * bytes;
        const uint8_t *src = in + y * bytes;
        for (size_t b = 0; b < 4; b++) {
            uint8_t *plane = dst + (little ? 3 - b : b) * width;
            for (size_t x = 0; x < width; x++) plane[x] = src[x * 4 + b];
        }
        for (size_t i = bytes; i-- > 1;) dst[i] = (uint8_t)(dst[i] - dst[i - 1]);
    }
}

/* libtiff remains serialized by tiff_lock. Only independent zlib streams run
 * in parallel; the calling thread writes their strips in original order.
 * 1 = use the scanline path (one worker, one strip, or an oversized row),
 * 0 = written, -1 = failure. Never allocate in proportion to image height. */
static int tiff_deflate_strips(TIFF *tif, const uint8_t *pixels, int H, int W, int px_bytes, int predictor)
{
    uint32_t rows = 0;
    size_t row_bytes = (size_t)W * (size_t)px_bytes, strip_bytes = 0, nstrips = 0;
    size_t bound = 0, slot_bytes = 0, workers = 1, start = 0, i = 0, initialized = 0;
    TiffDeflateSlot *slots = NULL;
    uint8_t *buffer = NULL;
    int rc = -1;
    Arena_T arena = NULL;
#ifdef _OPENMP
    if (!omp_in_parallel()) workers = (size_t)omp_get_max_threads();
#endif
    if (workers < 2 || TIFFIsByteSwapped(tif) ||
        !TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rows) || rows == 0) return 1;
    nstrips = ((size_t)H + rows - 1) / rows;
    if (nstrips < 2 || row_bytes > UINT_MAX / (size_t)rows) return 1;
    strip_bytes = row_bytes * rows;
    bound = (size_t)compressBound((uLong)strip_bytes);
    if (bound > TIFF_DEFLATE_BUFFER_BYTES || bound > UINT_MAX) return 1;
    slot_bytes = bound;
    if (predictor == PREDICTOR_FLOATINGPOINT) {
        if (px_bytes != 4 || strip_bytes > TIFF_DEFLATE_BUFFER_BYTES - bound) return 1;
        slot_bytes += strip_bytes;
    } else if (predictor != PREDICTOR_NONE) return 1;
    if (workers > TIFF_DEFLATE_BUFFER_BYTES / slot_bytes) workers = TIFF_DEFLATE_BUFFER_BYTES / slot_bytes;
    if (workers > nstrips) workers = nstrips;
    if (workers < 2) return 1;
    arena = Arena_new();
    slots = (TiffDeflateSlot *)ARENA_CALLOC(arena, workers, sizeof *slots);
    buffer = (uint8_t *)ARENA_ALLOC(arena, workers * slot_bytes);
    for (i = 0; i < workers; i++) {
        slots[i].data = buffer + i * slot_bytes;
        slots[i].predicted = predictor == PREDICTOR_FLOATINGPOINT ? slots[i].data + bound : NULL;
        if (deflateInit(&slots[i].stream, TIFF_DEFLATE_LEVEL) != Z_OK) goto done;
        initialized++;
    }
    for (start = 0; start < nstrips; start += workers) {
        const size_t count = nstrips - start < workers ? nstrips - start : workers;
        long j = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads((int)workers)
#endif
        for (j = 0; j < (long)count; j++) {
            TiffDeflateSlot *slot = &slots[j];
            const size_t y = (start + (size_t)j) * rows;
            const size_t nr = (size_t)H - y < rows ? (size_t)H - y : rows;
            slot->status = deflateReset(&slot->stream);
            if (slot->status != Z_OK) continue;
            const uint8_t *input = pixels + y * row_bytes;
            if (slot->predicted) {
                tiff_predict_float_rows(slot->predicted, input, nr, (size_t)W);
                input = slot->predicted;
            }
            slot->stream.next_in = (Bytef *)input;
            slot->stream.avail_in = (uInt)(nr * row_bytes);
            slot->stream.next_out = slot->data;
            slot->stream.avail_out = (uInt)bound;
            slot->status = deflate(&slot->stream, Z_FINISH);
            slot->bytes = (size_t)slot->stream.total_out;
        }
        for (i = 0; i < count; i++) {
            const TiffDeflateSlot *slot = &slots[i];
            if (slot->status != Z_STREAM_END ||
                TIFFWriteRawStrip(tif, (tstrip_t)(start + i), slot->data, (tmsize_t)slot->bytes)
                    != (tmsize_t)slot->bytes) goto done;
        }
    }
    rc = 0;
done:
    for (i = 0; i < initialized; i++) deflateEnd(&slots[i].stream);
    Arena_dispose(&arena);
    return rc;
}
#endif

static int tiff_save_u8(const char *path, const uint8_t *vol,
                        int D, int H, int W, const char *mode, int deflate)
{
    assert(path);
    assert(D >= 0 && H >= 0 && W >= 0);

    if (D == 0 || H == 0 || W == 0) {
        return -1;
    }

    assert(vol);

    ves_mutex_lock(&tiff_lock);

    TIFF *tif = TIFFOpen(path, mode);
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
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, tiff_strip_rows(H, W, 1));
        tiff_set_compression(tif, deflate);

        /* Multi-page TIFF subfile type, for volumes only.  A single 2-D raster (a baked sheet) must stay a plain
         * image: viewers (IrfanView, Paint) take the multipage path for a page-tagged file and refuse a page above
         * ~1 GiB of pixels, while the same pixels without the page tags open (measured 2026-09-27 on a
         * 1048576 x 1798 sheet). */
        if (D > 1) {
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER, (uint16_t)z, (uint16_t)D);
        }

        const uint8_t *page = vol + (size_t)z * (size_t)H * (size_t)W;
        int scanlines = 1;
#ifdef VES_TIFF_PARALLEL_DEFLATE
        if (deflate) scanlines = tiff_deflate_strips(tif, page, H, W, 1, PREDICTOR_NONE);
#endif
        if (scanlines < 0) {
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }
        for (int y = 0; scanlines && y < H; y++) {
            /* size_t row offset: a 1-vox 21x21x21 sheet is ~51k x 54k, and
             * y * W in int overflows past row 2^31 / W. */
            if (TIFFWriteScanline(tif, (void *)(page + (size_t)y * (size_t)W), (uint32_t)y, 0) < 0) {
                TIFFClose(tif);
                ves_mutex_unlock(&tiff_lock);
                return -1;
            }
        }

        if (!TIFFWriteDirectory(tif)) {
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }
    }

    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);
    return 0;
}

int TiffIO_save(const char *path, const uint8_t *vol, int D, int H, int W)
{
    size_t bytes;
    if (D <= 0 || H <= 0 || W <= 0 || (size_t)H > SIZE_MAX / (size_t)W)
        return -1;
    bytes = (size_t)H * (size_t)W;
    if ((size_t)D > SIZE_MAX / bytes) return -1;
    bytes *= (size_t)D;
    return tiff_save_u8(path, vol, D, H, W, tiff_write_mode(bytes, D), 0);
}

int TiffIO_save_deflate(const char *path, const uint8_t *vol, int D, int H, int W)
{
    size_t bytes;
    if (D <= 0 || H <= 0 || W <= 0 || (size_t)H > SIZE_MAX / (size_t)W)
        return -1;
    bytes = (size_t)H * (size_t)W;
    if ((size_t)D > SIZE_MAX / bytes) return -1;
    bytes *= (size_t)D;
    return tiff_save_u8(path, vol, D, H, W, tiff_write_mode(bytes, D), 1);
}

static int tiff_save_f32(const char *path, const float *img, int W, int H, int deflate, int predictor)
{
    assert(path);
    assert(W >= 0 && H >= 0);

    if (W <= 0 || H <= 0 || (size_t)W > SIZE_MAX / sizeof(float) / (size_t)H) {
        return -1;
    }

    assert(img);

    ves_mutex_lock(&tiff_lock);

    TIFF *tif = TIFFOpen(path,
        tiff_write_mode((size_t)W * (size_t)H * sizeof(float), 1));
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
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, tiff_strip_rows(H, W, 4));
    tiff_set_compression(tif, deflate);
    if (predictor != PREDICTOR_NONE && (!deflate || !TIFFSetField(tif, TIFFTAG_PREDICTOR, predictor))) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    int scanlines = 1;
#ifdef VES_TIFF_PARALLEL_DEFLATE
    if (deflate) scanlines = tiff_deflate_strips(tif, (const uint8_t *)img, H, W, 4, predictor);
#endif
    if (scanlines < 0) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    /* Older libtiff predictor codecs transform TIFFWriteScanline's input in
     * place. Keep the public const-input contract independently of version. */
    float *row_copy = scanlines && predictor != PREDICTOR_NONE ? (float *)malloc((size_t)W * sizeof(float)) : NULL;
    if (scanlines && predictor != PREDICTOR_NONE && !row_copy) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    for (int y = 0; scanlines && y < H; y++) {
        const float *row = img + (size_t)y * (size_t)W;
        if (row_copy) { memcpy(row_copy, row, (size_t)W * sizeof(float)); row = row_copy; }
        if (TIFFWriteScanline(tif, (void *)row, (uint32_t)y, 0) < 0) {
            free(row_copy);
            TIFFClose(tif);
            ves_mutex_unlock(&tiff_lock);
            return -1;
        }
    }
    free(row_copy);

    if (!TIFFWriteDirectory(tif)) {
        TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        return -1;
    }
    TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);
    return 0;
}

int TiffIO_save_float2d(const char *path, const float *img, int W, int H)
{
    return tiff_save_f32(path, img, W, H, 0, PREDICTOR_NONE);
}

int TiffIO_save_float2d_deflate(const char *path, const float *img, int W, int H)
{
    return tiff_save_f32(path, img, W, H, 1, PREDICTOR_NONE);
}

int TiffIO_save_float2d_predicted(const char *path, const float *img, int W, int H)
{
    return tiff_save_f32(path, img, W, H, 1, PREDICTOR_FLOATINGPOINT);
}

/* Compare decoded strips to the source without allocating another plane. */
static int tiff_verify_plane(const char *path, const void *pixels, int W, int H, int pixel_bytes)
{
    uint8_t *buffer = NULL;
    int failed = 0;
    uint32_t w = 0, h = 0, rows = 0;
    uint16_t bits = 0, samples = 0;
    ves_mutex_lock(&tiff_lock);
    TIFF *tif = TIFFOpen(path, "r");
    if (!tif) { failed = 1; goto done; }
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &samples);
    TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rows);
    if (w != (uint32_t)W || h != (uint32_t)H || bits != pixel_bytes * 8 || samples != 1 ||
        !rows || TIFFIsTiled(tif) || TIFFNumberOfDirectories(tif) != 1) { failed = 1; goto done; }
    if (rows > h) rows = h;
    const size_t row_bytes = (size_t)W * (size_t)pixel_bytes;
    if (rows > SIZE_MAX / row_bytes) { failed = 1; goto done; }
    buffer = (uint8_t *)malloc((size_t)rows * row_bytes);
    if (!buffer) { failed = 1; goto done; }
    for (uint32_t y = 0; y < h; y += rows) {
        const uint32_t nr = h - y < rows ? h - y : rows;
        const size_t bytes = (size_t)nr * row_bytes;
        if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, y, 0), buffer, (tmsize_t)bytes) != (tmsize_t)bytes ||
            memcmp(buffer, (const uint8_t *)pixels + (size_t)y * row_bytes, bytes)) { failed = 1; break; }
    }
done:
    free(buffer);
    if (tif) TIFFClose(tif);
    ves_mutex_unlock(&tiff_lock);
    return failed ? -1 : 0;
}

int TiffIO_storage_probe(const char *input, const char *fresh_directory)
{
    uint16_t bits = 0, samples = 0, format = 0;
    int single = 0, W = 0, H = 0, D = 0, rc = -1;
    uint8_t *mask = NULL;
    float *plane = NULL;
    double start = ves_clock_sec();
    if (!input || !fresh_directory) return -1;
    ves_mutex_lock(&tiff_lock);
    TIFF *tif = TIFFOpen(input, "r");
    if (tif) {
        TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits);
        TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &samples);
        TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &format);
        single = TIFFNumberOfDirectories(tif) == 1;
        TIFFClose(tif);
    }
    ves_mutex_unlock(&tiff_lock);
    if (!single || samples != 1 || !((bits == 32 && format == SAMPLEFORMAT_IEEEFP) ||
                                    (bits == 8 && format == SAMPLEFORMAT_UINT))) return -1;
    Arena_T arena = Arena_new();
    if (bits == 32 ? TiffIO_load_float2d(arena, input, &plane, &W, &H) :
                    (TiffIO_load(arena, input, &mask, &D, &H, &W) || D != 1)) goto done;
    if (ves_ensure_parent_dir(fresh_directory) || ves_mkdir(fresh_directory)) goto done;
    printf("{\"read_seconds\": %.6f, \"shape\": [%d, %d], \"pixel_bytes\": %u}\n",
           ves_clock_sec() - start, H, W, bits / 8);
    fflush(stdout);
    for (int arm = 0; arm < 4; arm++) {
        const int candidate = arm == 1 || arm == 2;
        char path[2400];
        snprintf(path, sizeof path, "%s/arm%d_%s.tif", fresh_directory, arm, candidate ? "candidate" : "legacy");
        start = ves_clock_sec();
        int result = bits == 32 ? (candidate ? TiffIO_save_float2d_predicted(path, plane, W, H) :
                                              TiffIO_save_float2d_deflate(path, plane, W, H)) :
                                 (candidate ? TiffIO_save_deflate(path, mask, 1, H, W) :
                                              TiffIO_save(path, mask, 1, H, W));
        double write_seconds = ves_clock_sec() - start;
        if (result) goto done;
        start = ves_clock_sec();
        if (tiff_verify_plane(path, plane ? (const void *)plane : (const void *)mask, W, H, bits / 8)) goto done;
        double verify_seconds = ves_clock_sec() - start;
        FILE *file = fopen(path, "rb");
        if (!file) goto done;
#ifdef _WIN32
        int failed = _fseeki64(file, 0, SEEK_END);
        int64_t bytes = _ftelli64(file);
#else
        int failed = fseeko(file, 0, SEEK_END);
        int64_t bytes = (int64_t)ftello(file);
#endif
        fclose(file);
        if (failed || bytes < 0) goto done;
        printf("{\"arm\": %d, \"candidate\": %s, \"bytes\": %llu, \"write_seconds\": %.6f, "
               "\"verify_seconds\": %.6f, \"all_pixel_bits_identical\": true}\n", arm, candidate ? "true" : "false",
               (unsigned long long)bytes, write_seconds, verify_seconds);
        fflush(stdout);
    }
    rc = 0;
done:
    Arena_dispose(&arena);
    return rc;
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

    int px_bytes = 0;
    if (bps == 8 && (fmt == SAMPLEFORMAT_UINT || fmt == 0)) px_bytes = 1;
    else if (bps == 32 && fmt == SAMPLEFORMAT_IEEEFP) px_bytes = 4;
    if (w == 0 || h == 0 || w > INT32_MAX || h > INT32_MAX || px_bytes == 0 || spp != 1
        || comp != COMPRESSION_NONE
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
    rd->px_bytes = px_bytes;
    rd->rows_per_strip = rps;
    rd->strip_off = off_copy;
    rd->n_strips = n_strips;
    return 0;
}

static int rows_read_bytes(TiffRowReader *rd, int y, int x0, int w, void *dst, int px_bytes)
{
    if (rd->f == NULL || rd->px_bytes != px_bytes || y < 0 || y >= rd->H || x0 < 0 || w < 0
        || x0 > rd->W || w > rd->W - x0)
        return -1;
    if (w == 0)
        return 0;

    uint32_t strip = (uint32_t)y / rd->rows_per_strip;
    uint32_t row_in = (uint32_t)y % rd->rows_per_strip;
    if (strip >= rd->n_strips)
        return -1;

    uint64_t off = rd->strip_off[strip]
                 + ((uint64_t)row_in * (uint64_t)rd->W + (uint64_t)x0) * (uint64_t)px_bytes;
    if (rows_seek(rd->f, off) != 0)
        return -1;
    if (fread(dst, (size_t)px_bytes, (size_t)w, rd->f) != (size_t)w)
        return -1;
    return 0;
}

int TiffIO_rows_read(TiffRowReader *rd, int y, int x0, int w, uint8_t *dst)
{
    assert(rd);
    assert(dst);
    return rows_read_bytes(rd, y, x0, w, dst, 1);
}

int TiffIO_rows_read_f32(TiffRowReader *rd, int y, int x0, int w, float *dst)
{
    assert(rd);
    assert(dst);
    return rows_read_bytes(rd, y, x0, w, dst, 4);
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
                    (size_t)w * (size_t)h * sizeof(float));

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
        arena, (size_t)w * (size_t)h * sizeof *img);
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

static int tiff_window_matches(const char *path, const float *pixels, int W, int H, size_t budget)
{
    int w = 0, h = 0, failed = 0;
    TiffFloatWindow *r = TiffIO_float_window_open(path, budget, &w, &h);
    if (!r) return 1;
    float got[71 * 35];
    if (w != W || h != H) failed++;
    /* Alternate edges, then cross every strip/tile row boundary. A two-unit
     * budget forces eviction on the tiled fixture; the final partial unit and
     * arbitrary IEEE bit patterns must also copy unchanged. */
    for (int turn = 0; turn < 2; turn++) {
        const int bands = (H - 1) / r->unit_H + 1;
        for (int band = -2; band < bands; band++) {
            int y = band == -2 ? H - 1 : band == -1 ? 0 : band * r->unit_H;
            if (y > 0) y--;
            const int rh = H - y < 35 ? H - y : 35;
            for (int col = 0; col < 9; col++) {
                const int x = col == 0 ? W - 1 : (col - 1) * ((W - 1) / 7);
                const int cw = W - x < 71 ? W - x : 71;
                if (TiffIO_float_window_read(r, y, x, rh, cw, got)) { failed++; continue; }
                for (int row = 0; row < rh; row++)
                    if (memcmp(got + (size_t)row * cw, pixels + (size_t)(y + row) * W + x,
                               (size_t)cw * sizeof(float))) failed++;
            }
        }
    }
    if (!TiffIO_float_window_read(r, 0, W, 1, 1, got) ||
        !TiffIO_float_window_read(r, H - 1, 0, 2, 1, got) ||
        !TiffIO_float_window_read(r, -1, 0, 1, 1, got)) failed++;
    TiffIO_float_window_close(r);
    if (failed) fprintf(stderr, "[tiff window selftest] %s: %d failures\n", path, failed);
    return failed;
}

int TiffIO_selftest(const char *dir)
{
    const int W = 4099, H = 2049, D = 2;
    const size_t bytes = (size_t)W * H * D;
    uint8_t *pixels = (uint8_t *)malloc(bytes);
    int fails = 0;
    if (!pixels) return 1;
    for (size_t i = 0; i < bytes; i++)
        pixels[i] = (uint8_t)((i * 31 + i / W * 17) % 251);
    if (strcmp(tiff_write_mode((size_t)971188 * 5219, 1), "w8") != 0 ||
        strcmp(tiff_write_mode((size_t)117407 * 2181, 1), "w") != 0 ||
        strcmp(tiff_write_mode((size_t)30000 * 40000 * 4, 1), "w8") != 0)
        fails++;
    /* Exercise both file formats with multiple strips and pages, without a
     * multi-gigabyte allocation in the routine selftest. Production 21x
     * bakes additionally exercise real offsets beyond classic TIFF's limit. */
    for (int big = 0; big < 2; big++) {
        char path[2048];
        uint8_t *loaded = NULL, row[29];
        int d = 0, h = 0, w = 0;
        Arena_T arena = Arena_new_sized(bytes + 4096);
        TiffRowReader rd;
        snprintf(path, sizeof path, "%s/tiff_%s.tif", dir, big ? "big" : "classic");
        int rc = big ? tiff_save_u8(path, pixels, D, H, W, "w8", 0)
                     : TiffIO_save(path, pixels, D, H, W);
        if (rc != 0) { fails++; Arena_dispose(&arena); continue; }
        ves_mutex_lock(&tiff_lock);
        TIFF *tif = TIFFOpen(path, "r");
        if (!tif || !!TIFFIsBigTIFF(tif) != big) fails++;
        if (tif) TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
        if (TiffIO_rows_open(arena, path, &rd) != 0) fails++;
        else {
            if (rd.W != W || rd.H != H || rd.n_strips < 2) fails++;
            int rows[4] = {0, (int)rd.rows_per_strip - 1, (int)rd.rows_per_strip, H - 1};
            for (int j = 0; j < 4; j++) {
                if (TiffIO_rows_read(&rd, rows[j], W - 29, 29, row) != 0 ||
                    memcmp(row, pixels + (size_t)rows[j] * W + W - 29, 29)) fails++;
            }
            if (TiffIO_rows_read(&rd, 0, INT32_MAX, 29, row) == 0) fails++;
            TiffIO_rows_close(&rd);
        }
        if (TiffIO_load(arena, path, &loaded, &d, &h, &w) != 0 ||
            d != D || h != H || w != W || memcmp(loaded, pixels, bytes)) fails++;
        Arena_dispose(&arena);
    }
    /* The row reader on a 32-bit float plane (a tifxyz coordinate plane): a window of a row equals the saved
     * image; the 8-bit read on a float reader is refused. */
    {
        Arena_T arena = Arena_new();
        enum { FW = 301, FH = 37 };
        float *img = (float *)ARENA_ALLOC(arena, (size_t)FW * FH * sizeof(float));
        float got[29];
        uint8_t bad[29];
        TiffRowReader rd;
        char path[512];
        int i = 0, rows[3] = { 0, 17, FH - 1 }, j = 0;
        for (i = 0; i < FW * FH; i++) img[i] = (float)i * 0.25f - 1.0f;
        snprintf(path, sizeof path, "%s/tiff_rows_f32.tif", dir);
        if (TiffIO_save_float2d(path, img, FW, FH) != 0 || TiffIO_rows_open(arena, path, &rd) != 0) fails++;
        else {
            if (rd.W != FW || rd.H != FH || rd.px_bytes != 4) fails++;
            for (j = 0; j < 3; j++)
                if (TiffIO_rows_read_f32(&rd, rows[j], FW - 29, 29, got) != 0
                    || memcmp(got, img + (size_t)rows[j] * FW + FW - 29, 29 * sizeof(float))) fails++;
            if (TiffIO_rows_read(&rd, 0, 0, 29, bad) == 0) fails++;
            TiffIO_rows_close(&rd);
        }
        fails += tiff_window_matches(path, img, FW, FH, (size_t)64 << 20);
        /* A compressed tiled float image, including right/bottom padding. */
        snprintf(path, sizeof path, "%s/tiff_window_tiled_f32.tif", dir);
        ves_mutex_lock(&tiff_lock);
        TIFF *tiled = TIFFOpen(path, "w");
        int tiled_bad = !tiled;
        if (tiled) {
            if (!TIFFSetField(tiled, TIFFTAG_IMAGEWIDTH, (uint32_t)FW) ||
                !TIFFSetField(tiled, TIFFTAG_IMAGELENGTH, (uint32_t)FH) ||
                !TIFFSetField(tiled, TIFFTAG_BITSPERSAMPLE, 32) ||
                !TIFFSetField(tiled, TIFFTAG_SAMPLESPERPIXEL, 1) ||
                !TIFFSetField(tiled, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP) ||
                !TIFFSetField(tiled, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG) ||
                !TIFFSetField(tiled, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK) ||
                !TIFFSetField(tiled, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE) ||
                !TIFFSetField(tiled, TIFFTAG_PREDICTOR, PREDICTOR_FLOATINGPOINT) ||
                !TIFFSetField(tiled, TIFFTAG_TILEWIDTH, 64) || !TIFFSetField(tiled, TIFFTAG_TILELENGTH, 32)) tiled_bad++;
            float tile[64 * 32];
            for (int y = 0; !tiled_bad && y < FH; y += 32)
                for (int x = 0; !tiled_bad && x < FW; x += 64) {
                    memset(tile, 0, sizeof tile);
                    const int cw = FW - x < 64 ? FW - x : 64;
                    for (int row = 0; row < 32 && y + row < FH; row++)
                        memcpy(tile + row * 64, img + (size_t)(y + row) * FW + x, (size_t)cw * sizeof(float));
                    if (TIFFWriteEncodedTile(tiled, TIFFComputeTile(tiled, (uint32_t)x, (uint32_t)y, 0, 0),
                                             tile, sizeof tile) != sizeof tile) tiled_bad++;
                }
            TIFFClose(tiled);
        }
        ves_mutex_unlock(&tiff_lock);
        fails += tiled_bad;
        if (!tiled_bad) fails += tiff_window_matches(path, img, FW, FH, 2 * 64 * 32 * sizeof(float));
        int ignored_w = 0, ignored_h = 0;
        TiffFloatWindow *too_small = TiffIO_float_window_open(path, 128, &ignored_w, &ignored_h);
        if (too_small) { fails++; TiffIO_float_window_close(too_small); }
        Arena_dispose(&arena);
    }
    /* A single 2-D raster is a plain image; only a volume is page-tagged (viewers refuse large page-tagged
     * sheets). */
    for (int pages = 1; pages <= 2; pages++) {
        char path[2048];
        snprintf(path, sizeof path, "%s/tiff_pages%d.tif", dir, pages);
        if (TiffIO_save(path, pixels, pages, 64, 64) != 0) { fails++; continue; }
        ves_mutex_lock(&tiff_lock);
        TIFF *tif = TIFFOpen(path, "r");
        uint32_t subfile = 0;
        uint16_t first = 0, count = 0;
        int typed = tif && TIFFGetField(tif, TIFFTAG_SUBFILETYPE, &subfile) && subfile == FILETYPE_PAGE;
        int numbered = tif && TIFFGetField(tif, TIFFTAG_PAGENUMBER, &first, &count);
        if (!tif || typed != (pages > 1) || numbered != (pages > 1)) fails++;
        if (tif) TIFFClose(tif);
        ves_mutex_unlock(&tiff_lock);
    }
    /* Deflate strips (2026-10-21): a sparse canvas -- filler (0 / -1) with one painted patch -- reads back exactly
     * through TiffIO_load / TiffIO_load_float2d, is tagged deflate, takes a few percent of the plain file's bytes, and
     * the uncompressed-only row reader declines it (its callers fall back to TiffIO_load). */
    {
        enum { SW = 4099, SH = 2049 };
        Arena_T arena = Arena_new_sized((size_t)SW * SH * 6 + 4096);
        uint8_t *mask = (uint8_t *)ARENA_CALLOC(arena, (size_t)SW * SH, 1), *got8 = NULL;
        float *xyz = (float *)ARENA_ALLOC(arena, (size_t)SW * SH * sizeof(float)), *gotf = NULL;
        char path[2048], plain[2048];
        int d = 0, h = 0, w = 0;
        TiffRowReader rd;
        for (size_t i = 0; i < (size_t)SW * SH; i++) xyz[i] = -1.0f;
        for (int y = 300; y < 700; y++)
            for (int x = 1000; x < 1900; x++) {
                const size_t i = (size_t)y * SW + (size_t)x;
                mask[i] = 255;
                xyz[i] = 1234.5f + 0.37f * (float)x - 0.11f * (float)y;
            }
        snprintf(path, sizeof path, "%s/tiff_deflate_u8.tif", dir);
        snprintf(plain, sizeof plain, "%s/tiff_deflate_u8_plain.tif", dir);
        if (TiffIO_save_deflate(path, mask, 1, SH, SW) != 0 || TiffIO_save(plain, mask, 1, SH, SW) != 0 ||
            TiffIO_load(arena, path, &got8, &d, &h, &w) != 0 || d != 1 || h != SH || w != SW ||
            memcmp(got8, mask, (size_t)SW * SH) != 0) fails++;
        if (TiffIO_rows_open(arena, path, &rd) == 0) { fails++; TiffIO_rows_close(&rd); }
        ves_mutex_lock(&tiff_lock);
        {
            TIFF *tif = TIFFOpen(path, "r");
            uint16_t comp = 0;
            if (!tif || !TIFFGetField(tif, TIFFTAG_COMPRESSION, &comp) || comp != COMPRESSION_ADOBE_DEFLATE ||
                TIFFIsTiled(tif)) fails++;
            if (tif) TIFFClose(tif);
        }
        ves_mutex_unlock(&tiff_lock);
        {
            FILE *a = fopen(path, "rb"), *b = fopen(plain, "rb");
            long sa = -1, sb = -1;
            if (a) { fseek(a, 0, SEEK_END); sa = ftell(a); fclose(a); }
            if (b) { fseek(b, 0, SEEK_END); sb = ftell(b); fclose(b); }
            if (sa <= 0 || sb <= 0 || sa * 20 > sb) fails++;
        }
        snprintf(path, sizeof path, "%s/tiff_deflate_f32.tif", dir);
        snprintf(plain, sizeof plain, "%s/tiff_deflate_f32_plain.tif", dir);
        if (TiffIO_save_float2d_deflate(path, xyz, SW, SH) != 0 || TiffIO_save_float2d(plain, xyz, SW, SH) != 0 ||
            TiffIO_load_float2d(arena, path, &gotf, &w, &h) != 0 || w != SW || h != SH ||
            memcmp(gotf, xyz, (size_t)SW * SH * sizeof(float)) != 0) fails++;
        if (TiffIO_rows_open(arena, path, &rd) == 0) { fails++; TiffIO_rows_close(&rd); }
        fails += tiff_window_matches(path, xyz, SW, SH, (size_t)18 << 20);
        fails += tiff_window_matches(plain, xyz, SW, SH, (size_t)18 << 20);
        {
            FILE *a = fopen(path, "rb"), *b = fopen(plain, "rb");
            long sa = -1, sb = -1;
            if (a) { fseek(a, 0, SEEK_END); sa = ftell(a); fclose(a); }
            if (b) { fseek(b, 0, SEEK_END); sb = ftell(b); fclose(b); }
            if (sa <= 0 || sb <= 0 || sa * 10 > sb) fails++;
        }
        /* Multi-page, full and partial strips, and incompressible bytes must
         * retain page/strip order. With two workers the float fixture above
         * also crosses batches. Preserve float bits, including NaN payloads. */
        {
            const size_t n = (size_t)SW * SH;
            uint8_t *pages = (uint8_t *)ARENA_ALLOC(arena, 2 * n), *back = NULL;
            uint32_t seed = 17u;
            memcpy(pages, mask, n);
            for (size_t i = 0; i < n; i++) {
                seed = seed * 1664525u + 1013904223u;
                pages[n + i] = (uint8_t)(seed >> 24);
            }
            snprintf(path, sizeof path, "%s/tiff_deflate_pages.tif", dir);
            if (TiffIO_save_deflate(path, pages, 2, SH, SW) != 0 ||
                TiffIO_load(arena, path, &back, &d, &h, &w) != 0 || d != 2 || h != SH || w != SW ||
                memcmp(pages, back, 2 * n) != 0) fails++;
            {
                const uint32_t bits[4] = { UINT32_C(0x7fc00001), UINT32_C(0xff800000),
                                          UINT32_C(0x80000000), UINT32_C(0x00000001) };
                memcpy(xyz, bits, sizeof bits);
                memcpy(xyz + n - 4, bits, sizeof bits);
            }
            snprintf(path, sizeof path, "%s/tiff_deflate_float_bits.tif", dir);
            if (TiffIO_save_float2d_deflate(path, xyz, SW, SH) != 0 ||
                TiffIO_load_float2d(arena, path, &gotf, &w, &h) != 0 || w != SW || h != SH ||
                memcmp(gotf, xyz, n * sizeof(float)) != 0) fails++;
            snprintf(path, sizeof path, "%s/tiff_predicted_float_bits.tif", dir);
            if (TiffIO_save_float2d_predicted(path, xyz, SW, SH) != 0 ||
                memcmp(gotf, xyz, n * sizeof(float)) != 0 ||
                tiff_verify_plane(path, xyz, SW, SH, 4) != 0) fails++;
            fails += tiff_window_matches(path, xyz, SW, SH, (size_t)18 << 20);
            ves_mutex_lock(&tiff_lock);
            TIFF *predicted = TIFFOpen(path, "r");
            uint16_t predictor = 0;
            if (!predicted || !TIFFGetField(predicted, TIFFTAG_PREDICTOR, &predictor) ||
                predictor != PREDICTOR_FLOATINGPOINT) fails++;
            if (predicted) TIFFClose(predicted);
            ves_mutex_unlock(&tiff_lock);
            /* An awkward one-row case takes libtiff's own scanline encoder,
             * including byte-plane boundaries and arbitrary float bit patterns. */
            uint32_t tiny[67], seed_bits = 287u;
            for (size_t i = 0; i < 67; i++) { seed_bits = seed_bits * 1664525u + 1013904223u; tiny[i] = seed_bits; }
            snprintf(path, sizeof path, "%s/tiff_predicted_one_row.tif", dir);
            if (TiffIO_save_float2d_predicted(path, (const float *)tiny, 67, 1) ||
                tiff_verify_plane(path, tiny, 67, 1, 4)) fails++;
        }
        Arena_dispose(&arena);
    }
    free(pixels);
    fprintf(stderr, "[tiff I/O selftest] %s (%d failures)\n", fails ? "FAIL" : "ok", fails);
    return fails;
}
