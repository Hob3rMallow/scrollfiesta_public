/* libtiff-owned stream and bounded decoded-unit cache; included by tiff_io.c. */
struct TiffFloatWindow {
    TIFF *tif;
    uint8_t *pixels;
    size_t unit_bytes, count;
    uint64_t clock, used[8];
    uint32_t key[8];
    int W, H, unit_W, unit_H, tiled;
};

TiffFloatWindow *TiffIO_float_window_open(const char *path, size_t cache_bytes, int *W, int *H)
{
    if (!path || !W || !H || !cache_bytes) return NULL;
    TiffFloatWindow *r = (TiffFloatWindow *)calloc(1, sizeof *r);
    if (!r) return NULL;
    ves_mutex_lock(&tiff_lock);
    r->tif = TIFFOpen(path, "r");
    uint32_t width = 0, height = 0, unit_w = 0, unit_h = 0;
    uint16_t bits = 0, samples = 0, format = 0, planar = 0, orientation = 0;
    if (!r->tif || !TIFFGetField(r->tif, TIFFTAG_IMAGEWIDTH, &width) ||
        !TIFFGetField(r->tif, TIFFTAG_IMAGELENGTH, &height) ||
        !TIFFGetFieldDefaulted(r->tif, TIFFTAG_BITSPERSAMPLE, &bits) ||
        !TIFFGetFieldDefaulted(r->tif, TIFFTAG_SAMPLESPERPIXEL, &samples) ||
        !TIFFGetFieldDefaulted(r->tif, TIFFTAG_SAMPLEFORMAT, &format) ||
        !TIFFGetFieldDefaulted(r->tif, TIFFTAG_PLANARCONFIG, &planar) ||
        !TIFFGetFieldDefaulted(r->tif, TIFFTAG_ORIENTATION, &orientation) ||
        bits != 32 || samples != 1 || format != SAMPLEFORMAT_IEEEFP || planar != PLANARCONFIG_CONTIG ||
        orientation != ORIENTATION_TOPLEFT || !width || !height || width > INT_MAX || height > INT_MAX ||
        !TIFFLastDirectory(r->tif)) goto fail;
    r->tiled = TIFFIsTiled(r->tif);
    if (r->tiled) {
        if (!TIFFGetField(r->tif, TIFFTAG_TILEWIDTH, &unit_w) ||
            !TIFFGetField(r->tif, TIFFTAG_TILELENGTH, &unit_h)) goto fail;
    } else {
        unit_w = width;
        if (!TIFFGetFieldDefaulted(r->tif, TIFFTAG_ROWSPERSTRIP, &unit_h)) goto fail;
        if (unit_h > height) unit_h = height;
    }
    if (!unit_w || !unit_h || unit_w > INT_MAX || unit_h > INT_MAX ||
        (size_t)unit_w > SIZE_MAX / sizeof(float) / unit_h) goto fail;
    r->unit_bytes = (size_t)unit_w * unit_h * sizeof(float);
    const tmsize_t decoded = r->tiled ? TIFFTileSize(r->tif) : TIFFStripSize(r->tif);
    if (decoded <= 0 || (size_t)decoded != r->unit_bytes || r->unit_bytes > cache_bytes) goto fail;
    /* Refuse oversized encoded units too, before libtiff allocates its raw
     * decode buffer. Metadata tables are still proportional to TIFF units. */
    uint64_t *byte_counts = NULL;
    const uint32_t units = r->tiled ? TIFFNumberOfTiles(r->tif) : TIFFNumberOfStrips(r->tif);
    if (!units || !TIFFGetField(r->tif, r->tiled ? TIFFTAG_TILEBYTECOUNTS : TIFFTAG_STRIPBYTECOUNTS, &byte_counts))
        goto fail;
    for (uint32_t i = 0; i < units; i++)
        if (!byte_counts[i] || byte_counts[i] > cache_bytes) goto fail;
    r->count = cache_bytes / r->unit_bytes;
    if (r->count > 8) r->count = 8;
    if (r->count > units) r->count = units;
    r->pixels = (uint8_t *)malloc(r->count * r->unit_bytes);
    if (!r->pixels) goto fail;
    r->W = (int)width; r->H = (int)height;
    r->unit_W = (int)unit_w; r->unit_H = (int)unit_h;
    for (size_t i = 0; i < r->count; i++) r->key[i] = UINT32_MAX;
    *W = r->W; *H = r->H;
    ves_mutex_unlock(&tiff_lock);
    return r;
fail:
    if (r->tif) TIFFClose(r->tif);
    free(r->pixels);
    free(r);
    ves_mutex_unlock(&tiff_lock);
    return NULL;
}

/* Called under tiff_lock. Returned storage belongs to this reader. */
static const uint8_t *tiff_float_unit(TiffFloatWindow *r, int y, int x)
{
    const uint32_t key = r->tiled ? TIFFComputeTile(r->tif, (uint32_t)x, (uint32_t)y, 0, 0) :
                                  TIFFComputeStrip(r->tif, (uint32_t)y, 0);
    size_t slot = 0;
    for (size_t i = 0; i < r->count; i++) {
        if (r->key[i] == key) { r->used[i] = ++r->clock; return r->pixels + i * r->unit_bytes; }
        if (r->used[i] < r->used[slot]) slot = i;
    }
    uint8_t *data = r->pixels + slot * r->unit_bytes;
    r->key[slot] = UINT32_MAX;
    const tmsize_t got = r->tiled ? TIFFReadEncodedTile(r->tif, key, data, (tmsize_t)r->unit_bytes) :
                                  TIFFReadEncodedStrip(r->tif, key, data, (tmsize_t)r->unit_bytes);
    const int rows = r->H - y < r->unit_H ? r->H - y : r->unit_H;
    const size_t want = r->tiled ? r->unit_bytes : (size_t)rows * r->unit_W * sizeof(float);
    if (got < 0 || (size_t)got != want) return NULL;
    r->key[slot] = key;
    r->used[slot] = ++r->clock;
    return data;
}

int TiffIO_float_window_read(TiffFloatWindow *r, int y, int x, int h, int w, float *dst)
{
    if (!r || !dst || x < 0 || y < 0 || w <= 0 || h <= 0 || x >= r->W || y >= r->H ||
        w > r->W - x || h > r->H - y) return -1;
    int failed = 0;
    ves_mutex_lock(&tiff_lock);
    for (int uy = y / r->unit_H; !failed && uy <= (y + h - 1) / r->unit_H; uy++) {
        const int top = uy * r->unit_H;
        const int first_y = y > top ? y : top;
        const int rows = r->H - top < r->unit_H ? r->H - top : r->unit_H;
        const int last_y = y + h < top + rows ? y + h : top + rows;
        for (int ux = x / r->unit_W; ux <= (x + w - 1) / r->unit_W; ux++) {
            const int left = ux * r->unit_W;
            const int first_x = x > left ? x : left;
            const int cols = r->W - left < r->unit_W ? r->W - left : r->unit_W;
            const int last_x = x + w < left + cols ? x + w : left + cols;
            const uint8_t *unit = tiff_float_unit(r, top, left);
            if (!unit) { failed = 1; break; }
            for (int row = first_y; row < last_y; row++)
                memcpy(dst + (size_t)(row - y) * w + first_x - x,
                       unit + ((size_t)(row - top) * r->unit_W + first_x - left) * sizeof(float),
                       (size_t)(last_x - first_x) * sizeof(float));
        }
    }
    ves_mutex_unlock(&tiff_lock);
    return failed ? -1 : 0;
}

void TiffIO_float_window_close(TiffFloatWindow *r)
{
    if (!r) return;
    ves_mutex_lock(&tiff_lock);
    TIFFClose(r->tif);
    ves_mutex_unlock(&tiff_lock);
    free(r->pixels);
    free(r);
}
