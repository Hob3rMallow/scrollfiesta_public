/* rawtex_tifxyz.c -- see rawtex_tifxyz.h. */
#include "rawtex_tifxyz.h"

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/zarr_shard.h"

#include "tiffio.h"

#include "../common/ves_omp.h"

/* ------------------------------------------------------------- reflection */

size_t RawtexTifxyz_reflect(float *pos, int8_t *sign, size_t W, size_t H,
                            RawtexTifxyzPiece *pieces, size_t npieces)
{
    const size_t plane = W * H;
    size_t k = 0, done = 0;
    if (pos == NULL || pieces == NULL) return 0;
    for (k = 0; k < npieces; k++) {
        RawtexTifxyzPiece *p = &pieces[k];
        size_t r = 0, r1 = p->row1 < H ? p->row1 : H;
        p->reflected = 0;
        if (p->px_mirror <= p->px_team || p->col0 >= p->col1 || p->col1 > W) continue;
        for (r = p->row0; r < r1; r++) {
            const size_t row = r * W;
            size_t a = p->col0, b = p->col1 - 1;
            for (; a < b; a++, b--) {
                int ch = 0;
                for (ch = 0; ch < 3; ch++) {
                    float *q = pos + (size_t)ch * plane + row;
                    float t = q[a];
                    q[a] = q[b];
                    q[b] = t;
                }
                if (sign != NULL) {
                    int8_t t = sign[row + a];
                    sign[row + a] = sign[row + b];
                    sign[row + b] = t;
                }
            }
            if (sign != NULL)
                for (a = p->col0; a < p->col1; a++) sign[row + a] = (int8_t)-sign[row + a];
        }
        p->reflected = 1;
        done++;
    }
    return done;
}

/* ------------------------------------------------------------ canvas mask */

static int rt_valid(const float *pos, size_t i)
{
    return pos[i] == pos[i];            /* NaN = invalid node */
}

/* All eight trilinear corners of (z, y, x) are present RAW. */
static int rt_sample_supported(CubeTable *ct, double z, double y, double x)
{
    const long iz = (long)floor(z), iy = (long)floor(y), ix = (long)floor(x);
    int k = 0;
    for (k = 0; k < 8; k++)
        if (!cube_supported(ct, iz + (k & 1), iy + ((k >> 1) & 1), ix + ((k >> 2) & 1))) return 0;
    return 1;
}

/* Every layer of the column p + off_k * n (off_k = (k - (depth-1)/2) * step)
 * has its trilinear footprint in present RAW.  Fast path: the chunks meeting
 * the footprint's bounding box are all present; otherwise sample by sample. */
static int rt_column_supported(CubeTable *ct, const double p[3], const double n[3], int depth, double step)
{
    const double h = 0.5 * (double)(depth - 1) * step;
    const long C = ct->chunk;
    long lo[3], hi[3], cz = 0, cy = 0, cx = 0;
    int a = 0, all = 1, k = 0;
    for (a = 0; a < 3; a++) {
        double e0 = p[a] - h * n[a], e1 = p[a] + h * n[a];
        lo[a] = (long)floor(e0 < e1 ? e0 : e1);
        hi[a] = (long)floor(e0 < e1 ? e1 : e0) + 1;
    }
    for (cz = lo[0] - ((lo[0] % C) + C) % C; cz <= hi[0] && all; cz += C)
        for (cy = lo[1] - ((lo[1] % C) + C) % C; cy <= hi[1] && all; cy += C)
            for (cx = lo[2] - ((lo[2] % C) + C) % C; cx <= hi[2] && all; cx += C)
                if (!cube_supported(ct, cz, cy, cx)) all = 0;
    if (all) return 1;
    for (k = 0; k < depth; k++) {
        const double off = ((double)k - 0.5 * (double)(depth - 1)) * step;
        if (!rt_sample_supported(ct, p[0] + off * n[0], p[1] + off * n[1], p[2] + off * n[2])) return 0;
    }
    return 1;
}

int RawtexTifxyz_canvas_mask(const float *pos, const int8_t *sign, size_t W, size_t H,
                             CubeTable *ct, int depth, double step, int flip,
                             uint8_t *mask, RawtexTifxyzMaskStats *st)
{
    const size_t plane = W * H;
    size_t rendered = 0, clean = 0, polluted = 0, reversed = 0, unsupported = 0, inmask = 0;
    long r = 0;
    if (st != NULL) memset(st, 0, sizeof *st);
    if (pos == NULL || mask == NULL || depth < 1 || !(step > 0.0) || W == 0 || H == 0) return -1;
    memset(mask, 0, plane);
    if (W < 3 || H < 3) return 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) reduction(+:rendered,clean,polluted,reversed,unsupported,inmask)
#endif
    for (r = 0; r < (long)H - 2; r++) {
        size_t c = 0;
        for (c = 0; c + 2 < W; c++) {
            const size_t rr = (size_t)r;
            const size_t n11 = (rr + 1) * W + c + 1, nl = (rr + 1) * W + c, nr = (rr + 1) * W + c + 2;
            const size_t nu = rr * W + c + 1, nd = (rr + 2) * W + c + 1, n00 = rr * W + c;
            double xv[3], yv[3], N[3], len = 0.0, p[3], n[3];
            int a = 0;
            if (!rt_valid(pos, n11) || !rt_valid(pos, nl) || !rt_valid(pos, nr) || !rt_valid(pos, nu) ||
                !rt_valid(pos, nd))
                continue;
            /* vc_render's grid normal, (x, y, z) order: planes are (z, y, x) */
            for (a = 0; a < 3; a++) {
                const size_t pl = (size_t)(2 - a) * plane;
                xv[a] = (double)pos[pl + nr] - (double)pos[pl + nl];
                yv[a] = (double)pos[pl + nd] - (double)pos[pl + nu];
            }
            N[0] = xv[1] * yv[2] - xv[2] * yv[1];
            N[1] = xv[2] * yv[0] - xv[0] * yv[2];
            N[2] = xv[0] * yv[1] - xv[1] * yv[0];
            len = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
            if (!(len > 0.0)) continue;
            rendered++;
            if (!rt_valid(pos, n00)) { polluted++; continue; }
            clean++;
            if (sign != NULL && sign[n11] < 0) { reversed++; continue; }
            if (ct != NULL) {
                /* the bilinear cell centre and the (flipped) unit normal, in (z, y, x) */
                for (a = 0; a < 3; a++) {
                    const size_t pl = (size_t)a * plane;
                    p[a] = 0.25 * ((double)pos[pl + n00] + (double)pos[pl + nu] + (double)pos[pl + nl] +
                                   (double)pos[pl + n11]);
                    n[a] = (flip ? -1.0 : 1.0) * N[2 - a] / len;
                }
                if (!rt_column_supported(ct, p, n, depth, step)) { unsupported++; continue; }
            }
            mask[rr * W + c] = 255;
            inmask++;
        }
    }
    if (st != NULL) {
        st->rendered_px = rendered; st->clean_px = clean; st->polluted_px = polluted;
        st->reversed_px = reversed; st->unsupported_px = unsupported; st->mask_px = inmask;
    }
    return 0;
}

/* ----------------------------------------------------------------- writer */

static int rt_write_segment(const char *dir, const char *uuid, float *pos, size_t W, size_t H,
                            const char *extra_json, int packed)
{
    int (*save_f32)(const char *, const float *, int, int) = packed ? TiffIO_save_float2d_deflate
                                                                    : TiffIO_save_float2d;
    int (*save_u8)(const char *, const uint8_t *, int, int, int) = packed ? TiffIO_save_deflate : TiffIO_save;
    const size_t plane = W * H;
    char path[2400];
    uint8_t *m = NULL;
    double lo[3] = { INFINITY, INFINITY, INFINITY }, hi[3] = { -INFINITY, -INFINITY, -INFINITY };
    size_t i = 0, valid = 0, cells = 0, r = 0, c = 0;
    FILE *fp = NULL;
    int a = 0, rc = -1;
    if (dir == NULL || uuid == NULL || pos == NULL || W == 0 || H == 0 || W > 0x7fffffff || H > 0x7fffffff)
        return -1;
    snprintf(path, sizeof path, "%s/meta.json", dir);
    fp = fopen(path, "rb");
    if (fp != NULL) {
        fclose(fp);
        fprintf(stderr, "ERROR: refusing to replace the tifxyz segment %s\n", dir);
        return -1;
    }
    if (ves_ensure_parent_dir(path) != 0) return -1;
    m = (uint8_t *)malloc(plane);
    if (m == NULL) return -1;
    for (r = 0; r + 1 < H; r++)          /* clean cells: the isometric area at 1 node per voxel */
        for (c = 0; c + 1 < W; c++) {
            const size_t j = r * W + c;
            cells += rt_valid(pos, j) && rt_valid(pos, j + 1) && rt_valid(pos, j + W) && rt_valid(pos, j + W + 1);
        }
    for (i = 0; i < plane; i++) {
        if (!rt_valid(pos, i) || !(pos[plane + i] == pos[plane + i]) || !(pos[2 * plane + i] == pos[2 * plane + i])) {
            pos[i] = pos[plane + i] = pos[2 * plane + i] = -1.0f;
            m[i] = 0;
            continue;
        }
        m[i] = 255;
        valid++;
        for (a = 0; a < 3; a++) {         /* bbox in (x, y, z): plane 2 - a */
            const double v = (double)pos[(size_t)(2 - a) * plane + i];
            if (v < lo[a]) lo[a] = v;
            if (v > hi[a]) hi[a] = v;
        }
    }
    if (valid == 0) { free(m); return -1; }
    snprintf(path, sizeof path, "%s/x.tif", dir);
    if (save_f32(path, pos + 2 * plane, (int)W, (int)H) != 0) goto out;
    snprintf(path, sizeof path, "%s/y.tif", dir);
    if (save_f32(path, pos + plane, (int)W, (int)H) != 0) goto out;
    snprintf(path, sizeof path, "%s/z.tif", dir);
    if (save_f32(path, pos, (int)W, (int)H) != 0) goto out;
    snprintf(path, sizeof path, "%s/mask.tif", dir);
    if (save_u8(path, m, 1, (int)H, (int)W) != 0) goto out;
    snprintf(path, sizeof path, "%s/meta.json", dir);
    fp = fopen(path, "wb");
    if (fp == NULL) goto out;
    fprintf(fp, "{\n  \"type\": \"seg\",\n  \"uuid\": \"%s\",\n  \"format\": \"tifxyz\",\n  \"scale\": [1.0, 1.0],\n"
            "  \"bbox\": [[%.3f, %.3f, %.3f], [%.3f, %.3f, %.3f]],\n  \"area_vx2\": %zu,\n  \"valid_nodes\": %zu%s%s\n}\n",
            uuid, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], cells, valid,
            extra_json != NULL && extra_json[0] ? ",\n  " : "", extra_json != NULL ? extra_json : "");
    if (fclose(fp) != 0) goto out;
    rc = 0;
out:
    free(m);
    if (rc != 0) fprintf(stderr, "ERROR: cannot write the tifxyz segment %s\n", dir);
    return rc;
}

int RawtexTifxyz_write(const char *dir, const char *uuid, float *pos, size_t W, size_t H,
                       const char *extra_json)
{
    return rt_write_segment(dir, uuid, pos, W, H, extra_json, 0);
}

int RawtexTifxyz_write_packed(const char *dir, const char *uuid, float *pos, size_t W, size_t H,
                              const char *extra_json)
{
    return rt_write_segment(dir, uuid, pos, W, H, extra_json, 1);
}

/* -------------------------------------------------------- canvas -> sheet */

int RawtexTifxyz_canvas_to_sheet(const uint8_t *canvas, size_t cW, size_t cH,
                                 const RawtexTifxyzPiece *pieces, size_t npieces,
                                 uint8_t *sheet, size_t sW, size_t sH)
{
    long *to = NULL;
    size_t c = 0, k = 0;
    long r = 0;
    if (canvas == NULL || sheet == NULL || (npieces > 0 && pieces == NULL) || sW > cW || sH > cH || cW == 0)
        return -1;
    to = (long *)malloc(cW * sizeof *to);
    if (to == NULL) return -1;
    for (c = 0; c < cW; c++) to[c] = (long)c;
    for (k = 0; k < npieces; k++) {
        const RawtexTifxyzPiece *p = &pieces[k];
        const long c0 = (long)p->col0 - 1, c1 = (long)p->col1;
        long cc = 0;
        if (!p->reflected) continue;
        if (p->col0 >= p->col1 || p->col1 > cW) { free(to); return -1; }
        for (cc = c0 < 0 ? 0 : c0; cc < c1; cc++) to[cc] = c0 + c1 - 1 - cc;
    }
    memset(sheet, 0, sW * sH);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (r = 0; r < (long)sH; r++) {
        const uint8_t *src = canvas + (size_t)r * cW;
        uint8_t *dst = sheet + (size_t)r * sW;
        size_t cc = 0;
        for (cc = 0; cc < cW; cc++) {
            const long t = to[cc];
            if (t >= 0 && (size_t)t < sW) dst[t] = src[cc];
        }
    }
    free(to);
    return 0;
}

/* ------------------------------------------------- canvas volume -> sheet */

#define RT_OUT_CHUNK 128

/* canvas column -> sheet column (identity, reversed inside a reflected piece's [col0 - 1, col1)); -1 = none */
static long *rt_column_map(const RawtexTifxyzPiece *pieces, size_t npieces, size_t cW)
{
    long *to = (long *)malloc((cW ? cW : 1) * sizeof *to);
    size_t c = 0, k = 0;
    if (to == NULL) return NULL;
    for (c = 0; c < cW; c++) to[c] = (long)c;
    for (k = 0; k < npieces; k++) {
        const RawtexTifxyzPiece *p = &pieces[k];
        const long c0 = (long)p->col0 - 1, c1 = (long)p->col1;
        long cc = 0;
        if (!p->reflected) continue;
        if (p->col0 >= p->col1 || p->col1 > cW) { free(to); return NULL; }
        for (cc = c0 < 0 ? 0 : c0; cc < c1; cc++) to[cc] = c0 + c1 - 1 - cc;
    }
    return to;
}

static int rt_write_file(const char *path, const void *buf, size_t n)
{
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    if (n && fwrite(buf, 1, n, fp) != n) { fclose(fp); return -1; }
    return fclose(fp) == 0 ? 0 : -1;
}

static int rt_publish_directory(const char *partial, const char *out_dir)
{
#ifdef _WIN32
    DWORD error = ERROR_SUCCESS;
    int attempt = 0;
    /* A scanner can briefly open a just-written chunk without FILE_SHARE_DELETE.
     * Retry only publication; never regenerate the finished volume or replace an
     * existing destination. These are I/O retries, not numerical work limits. */
    for (attempt = 0; attempt < 40; attempt++) {
        if (MoveFileExA(partial, out_dir, MOVEFILE_WRITE_THROUGH)) return 0;
        error = GetLastError();
        if (error != ERROR_SHARING_VIOLATION && error != ERROR_ACCESS_DENIED && error != ERROR_LOCK_VIOLATION)
            break;
        if (GetFileAttributesA(out_dir) != INVALID_FILE_ATTRIBUTES) break;
        if (attempt < 39) ves_sleep_ms(50);
    }
    fprintf(stderr, "ERROR: directory publication Windows error %lu after %d attempts\n",
            (unsigned long)error, attempt < 40 ? attempt + 1 : 40);
    return -1;
#else
    return rename(partial, out_dir);
#endif
}

int RawtexTifxyz_canvas_to_sheet_zarr_storage(const char *canvas_array_dir, const RawtexTifxyzPiece *pieces,
                                            size_t npieces, const uint8_t *sheet_mask,
                                            const uint8_t *sheet_reverse, size_t sW, size_t sH,
                                            const char *out_dir, const char *attrs_json, size_t shard_rows, size_t *chunks_written,
                                            size_t *canvas_chunks_absent)
{
    Arena_T a = NULL;
    const JsonValue *j = NULL, *shp = NULL, *chk = NULL, *comp = NULL;
    const char *dt = NULL, *sep = ".";
    char sepc = '.';                    /* the separator outlives the arena */
    char path[2400], partial[2400], text[1024];
    size_t D = 0, cH = 0, cW = 0, C = 0, band = 0, nbands = 0, ncx = 0, written = 0, absent = 0;
    ZarrShardArray shards = {0};
    ZarrShardBuilder *builders = NULL;
    int sharded = 0;
    long *to = NULL, *src = NULL;
    uint8_t *buf = NULL;
    FILE *probe = NULL;
    int failed = 0;
    size_t t = 0;
    if (chunks_written != NULL) *chunks_written = 0;
    if (canvas_chunks_absent != NULL) *canvas_chunks_absent = 0;
    if (canvas_array_dir == NULL || out_dir == NULL || (npieces > 0 && pieces == NULL) || sW == 0 || sH == 0 ||
        shard_rows > ZARR_SHARD_MAX_CHUNKS)
        return -1;
    a = Arena_new();
    snprintf(path, sizeof path, "%s/.zarray", canvas_array_dir);
    j = Json_parse_file(a, path, NULL);
    if (j == NULL) {
        /* Existing but malformed v2 metadata must not fall through to v3. */
        probe = fopen(path, "rb");
        if (probe != NULL) { fclose(probe); Arena_dispose(&a); return -1; }
        if (ZarrShard_array_open(a, canvas_array_dir, &shards)) { Arena_dispose(&a); return -1; }
        D = shards.depth; cH = shards.height; cW = shards.width; C = shards.tile;
        sharded = 1;
    }
    if (!sharded) {
        shp = Json_object_get(j, "shape");
        chk = Json_object_get(j, "chunks");
        comp = Json_object_get(j, "compressor");
        dt = Json_as_string(Json_object_get(j, "dtype"));
        if (Json_object_get(j, "dimension_separator") != NULL) sep = Json_as_string(Json_object_get(j, "dimension_separator"));
        if (j == NULL || Json_array_len(shp) != 3 || Json_array_len(chk) != 3 || Json_type(comp) != JSON_NULL ||
            dt == NULL || (strcmp(dt, "|u1") != 0 && strcmp(dt, "uint8") != 0) || sep == NULL ||
            (strcmp(sep, ".") != 0 && strcmp(sep, "/") != 0)) {
            fprintf(stderr, "ERROR: %s is not an uncompressed uint8 Zarr v2 array\n", canvas_array_dir);
            Arena_dispose(&a);
            return -1;
        }
        D = (size_t)Json_as_double(Json_array_get(shp, 0), 0); cH = (size_t)Json_as_double(Json_array_get(shp, 1), 0);
        cW = (size_t)Json_as_double(Json_array_get(shp, 2), 0);
        C = (size_t)Json_as_double(Json_array_get(chk, 1), 0);
        sepc = sep[0];
    }
    if (D == 0 || D > 255 || (!sharded && (size_t)Json_as_double(Json_array_get(chk, 0), 0) != D) || C == 0 ||
        (!sharded && (size_t)Json_as_double(Json_array_get(chk, 2), 0) != C) || sW > cW || sH > cH ||
        cW > LONG_MAX || C > SIZE_MAX / C || D > SIZE_MAX / (C * C) ||
        cW > SIZE_MAX / ((size_t)RT_OUT_CHUNK * D)) {
        fprintf(stderr, "ERROR: canvas array %zu x %zu x %zu (chunks %zu deep, %zu wide) cannot hold a %zu x %zu sheet\n",
                D, cH, cW, (size_t)Json_as_double(Json_array_get(chk, 0), 0), C, sH, sW);
        Arena_dispose(&a);
        return -1;
    }
    Arena_dispose(&a);
    snprintf(path, sizeof path, "%s/.zgroup", out_dir);
    if ((probe = fopen(path, "rb")) != NULL) {
        fclose(probe);
        fprintf(stderr, "ERROR: refusing to replace existing %s\n", out_dir);
        return -1;
    }
    snprintf(path, sizeof path, "%s/zarr.json", out_dir);
    if ((probe = fopen(path, "rb")) != NULL) {
        fclose(probe);
        fprintf(stderr, "ERROR: refusing to replace existing %s\n", out_dir);
        return -1;
    }
    snprintf(partial, sizeof partial, "%s.partial", out_dir);
    snprintf(path, sizeof path, "%s/0", partial);
    if (ves_ensure_parent_dir(partial) != 0 || ves_mkdir(partial) != 0 || ves_mkdir(path) != 0) {
        fprintf(stderr, "ERROR: cannot create a fresh %s (remove a stale partial first)\n", path);
        return -1;
    }
    to = rt_column_map(pieces, npieces, cW);
    src = (long *)malloc(sW * sizeof *src);
    buf = (uint8_t *)malloc(D * (size_t)RT_OUT_CHUNK * cW);
    if (to == NULL || src == NULL || buf == NULL) { failed = 1; goto done; }
    for (t = 0; t < sW; t++) src[t] = -1;
    for (t = 0; t < cW; t++)
        if (to[t] >= 0 && (size_t)to[t] < sW) src[to[t]] = (long)t;
    ncx = (cW + C - 1) / C;
    nbands = (sH + RT_OUT_CHUNK - 1) / RT_OUT_CHUNK;
    const long nout = (long)((sW + RT_OUT_CHUNK - 1) / RT_OUT_CHUNK);
    if (shard_rows) {
        builders = (ZarrShardBuilder *)calloc((size_t)nout, sizeof *builders);
        if (!builders) { failed = 1; goto done; }
    }
    for (band = 0; band < nbands && !failed; band++) {
        const size_t y0 = band * RT_OUT_CHUNK, yh = sH - y0 < RT_OUT_CHUNK ? sH - y0 : RT_OUT_CHUNK;
        long group = 0, oc = 0;
        const size_t group_cols = sharded ? shards.cols : 1;
        const long ngroups = (long)((ncx + group_cols - 1) / group_cols);
        if (shard_rows && band % shard_rows == 0) {
            snprintf(path, sizeof path, "%s/0/c/0/%zu/placeholder", partial, band / shard_rows);
            if (ves_ensure_parent_dir(path)) { failed = 1; break; }
            for (oc = 0; oc < nout; oc++)
                if (ZarrShard_begin(&builders[oc], D * RT_OUT_CHUNK * RT_OUT_CHUNK, shard_rows)) failed = 1;
            if (failed) break;
        }
        memset(buf, 0, D * (size_t)RT_OUT_CHUNK * cW);
        /* the canvas rows [y0, y0 + yh) of every layer, all columns: buf[(k * 128 + row) * cW + col] */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(+:absent)
#endif
        for (group = 0; group < ngroups; group++) {
            size_t cy = y0 / C;
            const size_t first = (size_t)group * group_cols;
            const size_t last = first + group_cols < ncx ? first + group_cols : ncx;
            uint8_t *ch = (uint8_t *)malloc(D * C * C);
            if (ch == NULL) {
#ifdef _OPENMP
#pragma omp critical(rt_fail)
#endif
                failed = 1;
                continue;
            }
            for (; cy * C < y0 + yh; cy++) {
                char cp[2400];
                FILE *fp = NULL;
                size_t k = 0, r = 0, xw = 0;
                ZarrShard shard;
                memset(&shard, 0, sizeof shard);
                int opened = 0;
                if (sharded) {
                    snprintf(cp, sizeof cp, "%s/c/0/%zu/%ld", canvas_array_dir, cy / shards.rows, group);
                    opened = ZarrShard_open(&shard, cp, D * C * C, shards.count);
                    if (opened < 0) {
#ifdef _OPENMP
#pragma omp critical(rt_fail)
#endif
                        failed = 1;
                        break;
                    }
                }
                for (size_t cx = first; cx < last; cx++) {
                    const size_t cx0 = cx * C;
                    int bad_read = 0;
                    if (sharded) {
                        const size_t index = (cy % shards.rows) * group_cols + cx - first;
                        if (!opened || !shard.lengths[index]) { absent++; continue; }
                        bad_read = ZarrShard_read(&shard, index, ch) != 0;
                    } else {
                        if (sepc == '/') snprintf(cp, sizeof cp, "%s/0/%zu/%zu", canvas_array_dir, cy, cx);
                        else snprintf(cp, sizeof cp, "%s/0.%zu.%zu", canvas_array_dir, cy, cx);
                        fp = fopen(cp, "rb");
                        if (fp == NULL) { absent++; continue; }                 /* fill value 0 */
                        bad_read = fread(ch, 1, D * C * C, fp) != D * C * C;
                        fclose(fp);
                    }
                    if (bad_read) {
#ifdef _OPENMP
#pragma omp critical(rt_fail)
#endif
                        failed = 1;
                        break;
                    }
                    xw = cW - cx0 < C ? cW - cx0 : C;
                    for (k = 0; k < D; k++)
                        for (r = 0; r < C; r++) {
                            const size_t row = cy * C + r;
                            if (row < y0 || row >= y0 + yh) continue;
                            memcpy(buf + ((size_t)k * RT_OUT_CHUNK + (row - y0)) * cW + cx0, ch + (k * C + r) * C, xw);
                        }
                }
                ZarrShard_close(&shard);
            }
            free(ch);
        }
        if (failed) break;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) reduction(+:written)
#endif
        for (oc = 0; oc < nout; oc++) {
            const size_t x0 = (size_t)oc * RT_OUT_CHUNK, xw = sW - x0 < RT_OUT_CHUNK ? sW - x0 : RT_OUT_CHUNK;
            uint8_t *out = (uint8_t *)calloc(D * RT_OUT_CHUNK * RT_OUT_CHUNK, 1);
            size_t k = 0, r = 0, x = 0;
            int any = 0;
            char op[2400];
            if (out == NULL) {
#ifdef _OPENMP
#pragma omp critical(rt_fail)
#endif
                failed = 1;
                continue;
            }
            for (r = 0; r < yh; r++)
                for (x = 0; x < xw; x++) {
                    const size_t sx = x0 + x;
                    const long c = src[sx];
                    if (c < 0) continue;
                    if (sheet_mask != NULL && !sheet_mask[(y0 + r) * sW + sx]) continue;
                    for (k = 0; k < D; k++) {
                        const size_t depth = sheet_reverse != NULL && sheet_reverse[(y0 + r) * sW + sx] ? D-1-k : k;
                        const uint8_t v = buf[(depth * RT_OUT_CHUNK + r) * cW + (size_t)c];
                        out[(k * RT_OUT_CHUNK + r) * RT_OUT_CHUNK + x] = v;
                        any |= v != 0 || sheet_mask != NULL;
                    }
                }
            int bad_write = 0;
            if (shard_rows) {
                snprintf(op, sizeof op, "%s/0/c/0/%zu/%ld", partial, band / shard_rows, oc);
                if (any) bad_write = ZarrShard_append(&builders[oc], op, band % shard_rows, out);
                if (!bad_write && ((band + 1) % shard_rows == 0 || band + 1 == nbands))
                    bad_write = ZarrShard_finish(&builders[oc], op);
            } else if (any) {
                snprintf(op, sizeof op, "%s/0/0.%zu.%ld", partial, band, oc);
                bad_write = rt_write_file(op, out, D * RT_OUT_CHUNK * RT_OUT_CHUNK);
            }
            if (bad_write) {
#ifdef _OPENMP
#pragma omp critical(rt_fail)
#endif
                failed = 1;
            } else if (any) written++;
            free(out);
        }
    }
    if (failed) { fprintf(stderr, "ERROR: canvas-to-sheet volume failed\n"); goto done; }
    if (shard_rows) {
        snprintf(path, sizeof path, "%s/0", partial);
        if (ZarrShard_array_meta_grid(path, D, sH, sW, RT_OUT_CHUNK, shard_rows, 1)) { failed = 1; goto done; }
        snprintf(path, sizeof path, "%s/zarr.json", partial);
        FILE *meta = fopen(path, "wb");
        if (!meta) { failed = 1; goto done; }
        failed = fprintf(meta, "{\n\"zarr_format\":3,\"node_type\":\"group\",\"attributes\":%s\n}\n",
                         attrs_json ? attrs_json : "{}") < 0;
        if (fclose(meta)) failed = 1;
        if (failed) goto done;
    } else {
    snprintf(path, sizeof path, "%s/0/.zarray", partial);
    snprintf(text, sizeof text,
             "{\n  \"zarr_format\": 2,\n  \"shape\": [%zu, %zu, %zu],\n  \"chunks\": [%zu, %d, %d],\n"
             "  \"dtype\": \"|u1\",\n  \"compressor\": null,\n  \"fill_value\": 0,\n"
             "  \"order\": \"C\",\n  \"filters\": null,\n  \"dimension_separator\": \".\"\n}\n",
             D, sH, sW, D, RT_OUT_CHUNK, RT_OUT_CHUNK);
    if (rt_write_file(path, text, strlen(text)) != 0) { failed = 1; goto done; }
    snprintf(path, sizeof path, "%s/.zgroup", partial);
    if (rt_write_file(path, "{\n  \"zarr_format\": 2\n}\n", strlen("{\n  \"zarr_format\": 2\n}\n")) != 0) {
        failed = 1;
        goto done;
    }
    snprintf(path, sizeof path, "%s/.zattrs", partial);
    if (rt_write_file(path, attrs_json != NULL ? attrs_json : "{}\n", strlen(attrs_json != NULL ? attrs_json : "{}\n")) != 0) {
        failed = 1;
        goto done;
    }
    }
    if (rt_publish_directory(partial, out_dir) != 0) {
        fprintf(stderr, "ERROR: cannot publish %s -> %s\n", partial, out_dir);
        failed = 1;
    }
done:
    free(to); free(src); free(buf); free(builders);
    if (chunks_written != NULL) *chunks_written = written;
    if (canvas_chunks_absent != NULL) *canvas_chunks_absent = absent;
    return failed ? -1 : 0;
}

int RawtexTifxyz_canvas_to_sheet_zarr_depth(const char *canvas_array_dir, const RawtexTifxyzPiece *pieces,
                                            size_t npieces, const uint8_t *sheet_mask,
                                            const uint8_t *sheet_reverse, size_t sW, size_t sH,
                                            const char *out_dir, const char *attrs_json, size_t *chunks_written,
                                            size_t *canvas_chunks_absent)
{
    return RawtexTifxyz_canvas_to_sheet_zarr_storage(canvas_array_dir, pieces, npieces, sheet_mask, sheet_reverse,
                                                     sW, sH, out_dir, attrs_json, 0, chunks_written, canvas_chunks_absent);
}

int RawtexTifxyz_canvas_to_sheet_zarr(const char *canvas_array_dir, const RawtexTifxyzPiece *pieces,
                                      size_t npieces, const uint8_t *sheet_mask, size_t sW, size_t sH,
                                      const char *out_dir, const char *attrs_json, size_t *chunks_written,
                                      size_t *canvas_chunks_absent)
{
    return RawtexTifxyz_canvas_to_sheet_zarr_depth(canvas_array_dir,pieces,npieces,sheet_mask,NULL,sW,sH,
                                                   out_dir,attrs_json,chunks_written,canvas_chunks_absent);
}

/* ---------------------------------------------------------------- contrast */

void RawtexTifxyz_histogram(const uint8_t *img, const uint8_t *mask, size_t n, size_t hist[256])
{
    size_t i = 0;
    if (img == NULL || mask == NULL || hist == NULL) return;
    for (i = 0; i < n; i++)
        if (mask[i]) hist[img[i]]++;
}

int RawtexTifxyz_window(const size_t hist[256], double pct_lo, double pct_hi, int *lo, int *hi)
{
    size_t cnt = 0, cum = 0, tlo = 0, thi = 0;
    int b = 0, l = 0, h = 255;
    if (hist == NULL || !(pct_lo >= 0.0) || !(pct_hi <= 100.0) || !(pct_lo < pct_hi)) return -1;
    for (b = 0; b < 256; b++) cnt += hist[b];
    if (cnt > 0) {                        /* obj_bake_raw's stretch_window, bin for bin */
        tlo = (size_t)(pct_lo / 100.0 * (double)cnt);
        thi = (size_t)(pct_hi / 100.0 * (double)cnt);
        for (b = 0, cum = 0; b < 256; b++) { cum += hist[b]; if (cum > tlo) { l = b; break; } }
        for (b = 0, cum = 0; b < 256; b++) { cum += hist[b]; if (cum >= thi) { h = b; break; } }
        if (h <= l) { l = 0; h = 255; }
    }
    if (lo != NULL) *lo = l;
    if (hi != NULL) *hi = h;
    return 0;
}

void RawtexTifxyz_apply_window(uint8_t *img, const uint8_t *mask, size_t n, int lo, int hi)
{
    uint8_t lut[256];
    size_t i = 0;
    int b = 0;
    if (img == NULL || mask == NULL) return;
    if (hi <= lo) { lo = 0; hi = 255; }
    for (b = 0; b < 256; b++) {
        double g = ((double)b - (double)lo) / (double)(hi - lo);
        g = g < 0.0 ? 0.0 : (g > 1.0 ? 1.0 : g);
        lut[b] = (uint8_t)floor(255.0 * g + 0.5);
    }
    for (i = 0; i < n; i++) img[i] = mask[i] ? lut[img[i]] : 0;
}

int RawtexTifxyz_stretch(uint8_t *img, const uint8_t *mask, size_t n, double pct_lo, double pct_hi,
                         int *lo, int *hi)
{
    size_t hist[256];
    int l = 0, h = 255;
    if (img == NULL || mask == NULL || !(pct_lo >= 0.0) || !(pct_hi <= 100.0) || !(pct_lo < pct_hi)) return -1;
    memset(hist, 0, sizeof hist);
    RawtexTifxyz_histogram(img, mask, n, hist);
    if (RawtexTifxyz_window(hist, pct_lo, pct_hi, &l, &h) != 0) return -1;
    RawtexTifxyz_apply_window(img, mask, n, l, h);
    if (lo != NULL) *lo = l;
    if (hi != NULL) *hi = h;
    return 0;
}

int RawtexTifxyz_compose(uint8_t *sheet, size_t W, size_t H, const uint8_t *block, size_t bW, size_t bH,
                         size_t col0)
{
    size_t r = 0;
    if (sheet == NULL || block == NULL || bH != H || bW == 0 || col0 > W || bW > W - col0) return -1;
    for (r = 0; r < H; r++) memcpy(sheet + r * W + col0, block + r * bW, bW);
    return 0;
}

/* ------------------------------------------------------------- gray TIFFs */

int RawtexTifxyz_load_gray(const char *path, uint8_t **img, size_t *W, size_t *H)
{
    TIFF *tif = NULL;
    uint32_t w = 0, h = 0, tw = 0, th = 0;
    uint16_t bps = 0, spp = 0;
    uint8_t *out = NULL, *buf = NULL;
    int ok = 0;
    if (path == NULL || img == NULL || W == NULL || H == NULL) return -1;
    *img = NULL;
    tif = TIFFOpen(path, "r");
    if (tif == NULL) { fprintf(stderr, "ERROR: cannot open %s\n", path); return -1; }
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    if (w == 0 || h == 0 || bps != 8 || spp != 1 || TIFFNumberOfDirectories(tif) != 1) {
        fprintf(stderr, "ERROR: %s is not a single-page 8-bit gray TIFF\n", path);
        TIFFClose(tif);
        return -1;
    }
    out = (uint8_t *)malloc((size_t)w * (size_t)h);
    if (out == NULL) { TIFFClose(tif); return -1; }
    if (TIFFIsTiled(tif)) {
        uint32_t x = 0, y = 0;
        TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
        TIFFGetField(tif, TIFFTAG_TILELENGTH, &th);
        buf = tw && th ? (uint8_t *)malloc((size_t)TIFFTileSize(tif)) : NULL;
        ok = buf != NULL;
        for (y = 0; ok && y < h; y += th)
            for (x = 0; ok && x < w; x += tw) {
                uint32_t cw = w - x < tw ? w - x : tw, ch = h - y < th ? h - y : th, r = 0;
                if (TIFFReadTile(tif, buf, x, y, 0, 0) < 0) { ok = 0; break; }
                for (r = 0; r < ch; r++) memcpy(out + (size_t)(y + r) * w + x, buf + (size_t)r * tw, cw);
            }
    } else {
        uint32_t rps = 0, y = 0;
        TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rps);
        if (rps == 0 || rps > h) rps = h;
        ok = 1;
        for (y = 0; ok && y < h; y += rps) {
            const uint32_t rows = h - y < rps ? h - y : rps;
            const tmsize_t want = (tmsize_t)rows * (tmsize_t)w;
            if (TIFFReadEncodedStrip(tif, TIFFComputeStrip(tif, y, 0), out + (size_t)y * w, want) != want) ok = 0;
        }
    }
    free(buf);
    TIFFClose(tif);
    if (!ok) { fprintf(stderr, "ERROR: cannot read %s\n", path); free(out); return -1; }
    *img = out; *W = w; *H = h;
    return 0;
}

int RawtexTifxyz_save_gray(const char *path, const uint8_t *img, size_t W, size_t H)
{
    const uint32_t T = 128;
    TIFF *tif = NULL;
    uint8_t *tile = NULL;
    size_t x = 0, y = 0;
    int ok = 1;
    if (path == NULL || img == NULL || W == 0 || H == 0 || W > 0xffffffffu || H > 0xffffffffu) return -1;
    if (ves_ensure_parent_dir(path) != 0) return -1;
    tif = TIFFOpen(path, W * H > ((size_t)1 << 31) ? "w8" : "w");
    tile = (uint8_t *)malloc((size_t)T * T);
    if (tif == NULL || tile == NULL) { if (tif) TIFFClose(tif); free(tile); return -1; }
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32_t)W);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32_t)H);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_TILEWIDTH, T);
    TIFFSetField(tif, TIFFTAG_TILELENGTH, T);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
    for (y = 0; ok && y < H; y += T)
        for (x = 0; ok && x < W; x += T) {
            const size_t cw = W - x < T ? W - x : T, ch = H - y < T ? H - y : T;
            size_t r = 0;
            memset(tile, 0, (size_t)T * T);
            for (r = 0; r < ch; r++) memcpy(tile + r * T, img + (y + r) * W + x, cw);
            if (TIFFWriteTile(tif, tile, (uint32_t)x, (uint32_t)y, 0, 0) < 0) ok = 0;
        }
    TIFFClose(tif);
    free(tile);
    if (!ok) fprintf(stderr, "ERROR: cannot write %s\n", path);
    return ok ? 0 : -1;
}

/* --------------------------------------------------------------- selftest */

#define RT_CHECK(cond, msg) do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", msg); fails++; } } while (0)

/* synthetic RAW: a 32^3 nested Zarr (chunk 16), value = x; chunk (cz, cy, cx) = (1, 0, 0) absent */
static int rt_make_zarr(const char *dir)
{
    char path[2400];
    int cz = 0, cy = 0, cx = 0;
    uint8_t *cube = (uint8_t *)malloc(16 * 16 * 16);
    if (cube == NULL) return -1;
    ves_mkdir(dir);
    snprintf(path, sizeof path, "%s/0", dir);
    ves_mkdir(path);
    snprintf(path, sizeof path, "%s/0/.zarray", dir);
    {
        FILE *fp = fopen(path, "wb");
        if (fp == NULL) { free(cube); return -1; }
        fputs("{\"zarr_format\": 2, \"shape\": [32, 32, 32], \"chunks\": [16, 16, 16], \"dtype\": \"|u1\", "
              "\"compressor\": null, \"fill_value\": 0, \"order\": \"C\", \"filters\": null, "
              "\"dimension_separator\": \"/\"}", fp);
        fclose(fp);
    }
    for (cz = 0; cz < 2; cz++)
        for (cy = 0; cy < 2; cy++)
            for (cx = 0; cx < 2; cx++) {
                int i = 0;
                FILE *fp = NULL;
                if (cz == 1 && cy == 0 && cx == 0) continue;
                snprintf(path, sizeof path, "%s/0/%d", dir, cz); ves_mkdir(path);
                snprintf(path, sizeof path, "%s/0/%d/%d", dir, cz, cy); ves_mkdir(path);
                for (i = 0; i < 16 * 16 * 16; i++) cube[i] = (uint8_t)(cx * 16 + i % 16);
                snprintf(path, sizeof path, "%s/0/%d/%d/%d", dir, cz, cy, cx);
                fp = fopen(path, "wb");
                if (fp == NULL) { free(cube); return -1; }
                fwrite(cube, 1, 16 * 16 * 16, fp);
                fclose(fp);
            }
    free(cube);
    return 0;
}

static void rt_rm_zarr(const char *dir)
{
    char path[2400];
    int cz = 0, cy = 0, cx = 0;
    for (cz = 0; cz < 2; cz++) {
        for (cy = 0; cy < 2; cy++) {
            for (cx = 0; cx < 2; cx++) {
                snprintf(path, sizeof path, "%s/0/%d/%d/%d", dir, cz, cy, cx);
                remove(path);
            }
            snprintf(path, sizeof path, "%s/0/%d/%d", dir, cz, cy); ves_rmdir(path);
        }
        snprintf(path, sizeof path, "%s/0/%d", dir, cz); ves_rmdir(path);
    }
    snprintf(path, sizeof path, "%s/0/.zarray", dir); remove(path);
    snprintf(path, sizeof path, "%s/0", dir); ves_rmdir(path);
    ves_rmdir(dir);
}

static void rt_rm_seg(const char *dir)
{
    const char *names[5] = { "x.tif", "y.tif", "z.tif", "mask.tif", "meta.json" };
    char path[2400];
    int i = 0;
    for (i = 0; i < 5; i++) {
        snprintf(path, sizeof path, "%s/%s", dir, names[i]);
        remove(path);
    }
    ves_rmdir(dir);
}

/* A W x H plane grid, plane-major (z, y, x): node (r, c) = (z0 + r, y0 + c, x0). */
static void rt_plane(float *pos, size_t W, size_t H, float z0, float y0, float x0)
{
    size_t r = 0, c = 0, plane = W * H;
    for (r = 0; r < H; r++)
        for (c = 0; c < W; c++) {
            pos[r * W + c] = z0 + (float)r;
            pos[plane + r * W + c] = y0 + (float)c;
            pos[2 * plane + r * W + c] = x0;
        }
}

#ifdef _WIN32
static DWORD WINAPI rt_test_release_directory(LPVOID opaque)
{
    HANDLE directory = (HANDLE)opaque;
    ves_sleep_ms(150);
    CloseHandle(directory);
    return 0;
}
#endif

int RawtexTifxyz_selftest(void)
{
    int fails = ZarrShard_selftest();
    char root[128];
    if (ves_selftest_dir(root, sizeof root, "rawtex_tifxyz")) return 1;
    char raw[256], seg[256];
    snprintf(raw, sizeof raw, "%s/raw.zarr", root);
    snprintf(seg, sizeof seg, "%s/seg", root);
    ves_mkdir(root);
    rt_rm_zarr(raw);
    rt_rm_seg(seg);

    /* A completed directory survives a short sharing violation, and a completed
     * destination is never replaced by a competing publication. */
    {
        char pending[300], published[300], original[350], final[350];
        char bytes[16] = {0};
        FILE *fp = NULL;
        snprintf(pending, sizeof pending, "%s/publish.partial", root);
        snprintf(published, sizeof published, "%s/publish", root);
        snprintf(original, sizeof original, "%s/marker", pending);
        snprintf(final, sizeof final, "%s/marker", published);
        remove(original); remove(final); ves_rmdir(pending); ves_rmdir(published);
        ves_mkdir(pending);
        RT_CHECK(rt_write_file(original, "complete", 8) == 0, "publication marker written");
#ifdef _WIN32
        {
            HANDLE locked = CreateFileA(pending, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE,
                                        NULL, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, NULL);
            HANDLE thread = NULL;
            RT_CHECK(locked != INVALID_HANDLE_VALUE, "directory lock opened");
            if (locked != INVALID_HANDLE_VALUE) {
                RT_CHECK(!MoveFileExA(pending, published, MOVEFILE_WRITE_THROUGH), "test lock blocks publication");
                thread = CreateThread(NULL, 0, rt_test_release_directory, locked, 0, NULL);
                RT_CHECK(thread != NULL, "directory release thread started");
                if (thread == NULL) CloseHandle(locked);
            }
            RT_CHECK(rt_publish_directory(pending, published) == 0, "publication retries a sharing violation");
            if (thread != NULL) { WaitForSingleObject(thread, INFINITE); CloseHandle(thread); }
        }
#else
        RT_CHECK(rt_publish_directory(pending, published) == 0, "completed directory published");
#endif
        fp = fopen(final, "rb");
        if (fp != NULL) { (void)fread(bytes, 1, 8, fp); fclose(fp); }
        RT_CHECK(strcmp(bytes, "complete") == 0, "published bytes preserved");
        ves_mkdir(pending);
        RT_CHECK(rt_write_file(original, "other", 5) == 0, "competing marker written");
        RT_CHECK(rt_publish_directory(pending, published) != 0, "completed destination retained");
        memset(bytes, 0, sizeof bytes);
        fp = fopen(final, "rb");
        if (fp != NULL) { (void)fread(bytes, 1, 8, fp); fclose(fp); }
        RT_CHECK(strcmp(bytes, "complete") == 0, "publication conflict leaves original bytes");
        remove(original); remove(final); ves_rmdir(pending); ves_rmdir(published);
    }

    /* 1. reflection: a mirror-majority piece's columns are reversed (and its
     * signs negated), a team piece's untouched; reflecting twice restores both */
    {
        enum { W = 12, H = 3 };
        float pos[3 * W * H], orig[3 * W * H];
        int8_t sign[W * H];
        RawtexTifxyzPiece pc[2];
        size_t i = 0, n = 0;
        int bad = 0;
        for (i = 0; i < 3 * W * H; i++) pos[i] = NAN;
        for (i = 0; i < W * H; i++) {
            size_t col = i % W;
            sign[i] = 0;
            if ((col >= 1 && col < 5) || (col >= 7 && col < 11)) {
                pos[i] = (float)i; pos[W * H + i] = (float)(100 + i); pos[2 * W * H + i] = (float)(200 + i);
                sign[i] = (int8_t)(col < 5 ? -1 : 1);
            }
        }
        memcpy(orig, pos, sizeof pos);
        memset(pc, 0, sizeof pc);
        pc[0].col0 = 1; pc[0].col1 = 5; pc[0].row0 = 0; pc[0].row1 = H; pc[0].px_mirror = 12;
        pc[1].col0 = 7; pc[1].col1 = 11; pc[1].row0 = 0; pc[1].row1 = H; pc[1].px_team = 12;
        n = RawtexTifxyz_reflect(pos, sign, W, H, pc, 2);
        RT_CHECK(n == 1 && pc[0].reflected == 1 && pc[1].reflected == 0, "one piece reflected");
        for (i = 0; i < W * H; i++) {
            size_t col = i % W, row = i / W, src = i;
            if (col >= 1 && col < 5) src = row * W + (1 + 5 - 1 - col);
            if (!(pos[i] == orig[src] || (isnan(pos[i]) && isnan(orig[src])))) bad++;
            if (col >= 1 && col < 11 && !(col >= 5 && col < 7) && sign[i] != 1) bad++;
        }
        RT_CHECK(bad == 0, "reflected columns reversed, team piece untouched, all signs now +1");
        RawtexTifxyz_reflect(pos, sign, W, H, pc, 1);   /* the census is still the input's: reflect again */
        for (i = 0, bad = 0; i < 3 * W * H; i++)
            if (!(pos[i] == orig[i] || (isnan(pos[i]) && isnan(orig[i])))) bad++;
        RT_CHECK(bad == 0, "reflecting twice restores the grid");
    }

    /* 2. the canvas mask follows vc_render's rules: the plane z = 10 + r,
     * y = 4 + c, x = 12.3 on a 10 x 10 grid; (dP/dcol x dP/drow) in (x, y, z) is
     * +x, so --flip-normals renders along -x.  Canvas pixels (r, c) with
     * r, c <= 7 render (their node (r+1, c+1) has four neighbours) = 64.
     * Removing node (3, 3) removes the six pixels that use it; (3, 3) itself is
     * still RENDERED by vc_render (its node (4, 4) and neighbours are valid) --
     * polluted.  A negative sign at node (2, 6) drops pixel (1, 5).  With RAW:
     * layers at x = 12.3 +- 2 stay in chunk cx 0, but rows r >= 5 have cell
     * centres z = 10.5 + r >= 15.5 whose trilinear corners reach the absent
     * chunk (1, 0, 0). */
    {
        enum { W = 10, H = 10 };
        float pos[3 * W * H];
        int8_t sign[W * H];
        uint8_t mask[W * H];
        RawtexTifxyzMaskStats ms;
        Arena_T arena = NULL;
        CubeTable ct;
        size_t i = 0;
        int on = 0, bad = 0;
        rt_plane(pos, W, H, 10.0f, 4.0f, 12.3f);
        memset(sign, 1, sizeof sign);
        RT_CHECK(RawtexTifxyz_canvas_mask(pos, NULL, W, H, NULL, 5, 1.0, 1, mask, &ms) == 0 && ms.rendered_px == 64 &&
                 ms.clean_px == 64 && ms.mask_px == 64, "a full plane: (W-2) x (H-2) canvas pixels render");
        for (i = 0; i < W * H; i++) {
            size_t r = i / W, c = i % W;
            on += mask[i] == 255;
            if ((mask[i] == 255) != (r <= 7 && c <= 7)) bad++;
        }
        RT_CHECK(on == 64 && bad == 0, "the rendered block is rows/cols 0..7");
        pos[3 * W + 3] = NAN; pos[W * H + 3 * W + 3] = NAN; pos[2 * W * H + 3 * W + 3] = NAN;
        sign[2 * W + 6] = -1;
        RT_CHECK(RawtexTifxyz_canvas_mask(pos, sign, W, H, NULL, 5, 1.0, 1, mask, &ms) == 0 && ms.rendered_px == 59 &&
                 ms.polluted_px == 1 && ms.clean_px == 58 && ms.reversed_px == 1 && ms.mask_px == 57,
                 "an invalid node: 5 pixels lose their normal, 1 renders polluted; one reversed pixel");
        RT_CHECK(mask[3 * W + 3] == 0 && mask[3 * W + 2] == 0 && mask[2 * W + 3] == 0 && mask[2 * W + 2] == 0 &&
                 mask[2 * W + 1] == 0 && mask[1 * W + 2] == 0 && mask[1 * W + 5] == 0 && mask[4 * W + 4] == 255,
                 "the pixels that use node (3, 3), and the reversed one, are out");
        rt_plane(pos, W, H, 10.0f, 4.0f, 12.3f);
        RT_CHECK(rt_make_zarr(raw) == 0, "synthetic RAW");
        arena = Arena_new();
        {   /* the cube table over the grid's points, padded by the layer reach */
            float v[3 * W * H];
            for (i = 0; i < W * H; i++) { v[i * 3] = pos[i]; v[i * 3 + 1] = pos[W * H + i]; v[i * 3 + 2] = pos[2 * W * H + i]; }
            RT_CHECK(cubetable_init(&ct, arena, raw, 16, v, W * H, 4.0) == 0, "cube table");
        }
        cubetable_prewarm_all(&ct);
        RT_CHECK(RawtexTifxyz_canvas_mask(pos, NULL, W, H, &ct, 5, 1.0, 1, mask, &ms) == 0 && ms.clean_px == 64 &&
                 ms.unsupported_px == 24 && ms.mask_px == 40, "rows reaching the absent chunk are unsupported");
        for (i = 0, bad = 0; i < W * H; i++) {
            size_t r = i / W, c = i % W;
            if ((mask[i] == 255) != (r <= 4 && c <= 7)) bad++;
        }
        RT_CHECK(bad == 0, "supported rows 0..4");
        rt_plane(pos, W, H, 2.0f, 4.0f, 14.5f);         /* layers x 12.5..16.5 cross into chunk cx 1: present */
        RT_CHECK(RawtexTifxyz_canvas_mask(pos, NULL, W, H, &ct, 5, 1.0, 1, mask, &ms) == 0 && ms.mask_px == 64,
                 "a column crossing into a present chunk is supported");
        {   /* a diagonal column whose bounding box meets the absent chunk but whose samples do not:
             * (z, x) = (12.5, 15.5) .. (16.5, 19.5) -- the fast path refuses, the per-sample pass accepts;
             * a column along x at z = 15 needs the corner z = 16 with x < 16: absent */
            const double s2 = sqrt(2.0);
            const double p[3] = { 14.5, 8.0, 17.5 }, n[3] = { 1.0 / s2, 0.0, 1.0 / s2 };
            const double q[3] = { 15.0, 8.0, 12.0 }, e[3] = { 0.0, 0.0, 1.0 };
            RT_CHECK(rt_column_supported(&ct, p, n, 5, s2) == 1, "a diagonal column beside an absent chunk");
            RT_CHECK(rt_column_supported(&ct, q, e, 3, 1.0) == 0, "a column touching an absent chunk");
        }
        Arena_dispose(&arena);
        rt_rm_zarr(raw);
    }

    /* 3. canvas -> sheet: a reflected piece's canvas columns [col0 - 1, col1)
     * are reversed, everything else is copied; the last canvas row/column drop */
    {
        enum { CW = 12, CH = 3, SW = 11, SH = 2 };
        uint8_t canvas[CW * CH], sheet[SW * SH];
        RawtexTifxyzPiece pc[2];
        size_t i = 0;
        int bad = 0;
        for (i = 0; i < CW * CH; i++) canvas[i] = (uint8_t)(10 * (i / CW) + i % CW);
        memset(pc, 0, sizeof pc);
        pc[0].col0 = 2; pc[0].col1 = 6; pc[0].reflected = 1;
        pc[1].col0 = 8; pc[1].col1 = 11; pc[1].reflected = 0;
        RT_CHECK(RawtexTifxyz_canvas_to_sheet(canvas, CW, CH, pc, 2, sheet, SW, SH) == 0, "canvas to sheet runs");
        for (i = 0; i < SW * SH; i++) {
            size_t r = i / SW, t = i % SW, c = t;
            if (t >= 1 && t < 6) c = 1 + 6 - 1 - t;
            if (sheet[i] != canvas[r * CW + c]) bad++;
        }
        RT_CHECK(bad == 0, "reflected canvas columns land reversed on the sheet");
        pc[0].col0 = 0; pc[0].col1 = 4;             /* a piece at column 0: canvas [-1, 4) */
        RT_CHECK(RawtexTifxyz_canvas_to_sheet(canvas, CW, CH, pc, 1, sheet, SW, SH) == 0 && sheet[3] == 0 &&
                 sheet[2] == canvas[0] && sheet[0] == canvas[2] && sheet[4] == canvas[4],
                 "a piece at column 0: the sheet column with no canvas source is 0");
    }

    /* 4. the writer: -1 sentinels, a readable grid, a meta.json with scale 1; no overwrite */
    {
        enum { W = 6, H = 5 };
        float pos[3 * W * H];
        Arena_T a = Arena_new();
        float *x = NULL;
        int w = 0, h = 0;
        FILE *fp = NULL;
        char path[300], text[1024];
        size_t got = 0;
        rt_plane(pos, W, H, 100.0f, 200.0f, 300.0f);
        pos[W + 2] = NAN; pos[W * H + W + 2] = NAN; pos[2 * W * H + W + 2] = NAN;
        RT_CHECK(RawtexTifxyz_write(seg, "unit_seg", pos, W, H, "\"note\": \"selftest\"") == 0, "segment written");
        snprintf(path, sizeof path, "%s/x.tif", seg);
        RT_CHECK(TiffIO_load_float2d(a, path, &x, &w, &h) == 0 && w == W && h == H && x[W + 2] == -1.0f &&
                 x[0] == 300.0f, "x.tif: x plane, -1 at the invalid node");
        snprintf(path, sizeof path, "%s/meta.json", seg);
        fp = fopen(path, "rb");
        if (fp != NULL) { got = fread(text, 1, sizeof text - 1, fp); fclose(fp); }
        text[got] = 0;
        RT_CHECK(strstr(text, "\"scale\": [1.0, 1.0]") != NULL && strstr(text, "\"uuid\": \"unit_seg\"") != NULL &&
                 strstr(text, "\"format\": \"tifxyz\"") != NULL && strstr(text, "\"note\": \"selftest\"") != NULL &&
                 strstr(text, "\"area_vx2\": 16") != NULL, "meta.json: scale 1, uuid, format, extras, area");
        rt_plane(pos, W, H, 100.0f, 200.0f, 300.0f);
        RT_CHECK(RawtexTifxyz_write(seg, "unit_seg", pos, W, H, NULL) != 0, "refuses to replace a segment");
        /* the packed writer (2026-10-21): the same planes, mask and meta.json in deflate strips */
        {
            static const char *planes[4] = { "x.tif", "y.tif", "z.tif", "mask.tif" };
            char seg2[300], p2[400];
            int k = 0;
            snprintf(seg2, sizeof seg2, "%s_packed", seg);
            rt_plane(pos, W, H, 100.0f, 200.0f, 300.0f);
            pos[W + 2] = NAN; pos[W * H + W + 2] = NAN; pos[2 * W * H + W + 2] = NAN;
            RT_CHECK(RawtexTifxyz_write_packed(seg2, "unit_seg", pos, W, H, "\"note\": \"selftest\"") == 0,
                     "packed segment written");
            for (k = 0; k < 4; k++) {
                uint8_t *u0 = NULL, *u1 = NULL;
                float *f0 = NULL, *f1 = NULL;
                int w0 = 0, h0 = 0, w1 = 0, h1 = 0, d0 = 0, d1 = 0, same = 0;
                uint16_t comp = 0;
                TIFF *tif = NULL;
                snprintf(path, sizeof path, "%s/%s", seg, planes[k]);
                snprintf(p2, sizeof p2, "%s/%s", seg2, planes[k]);
                if (k < 3)
                    same = TiffIO_load_float2d(a, path, &f0, &w0, &h0) == 0 &&
                           TiffIO_load_float2d(a, p2, &f1, &w1, &h1) == 0 && w0 == w1 && h0 == h1 &&
                           memcmp(f0, f1, (size_t)w0 * (size_t)h0 * sizeof(float)) == 0;
                else
                    same = TiffIO_load(a, path, &u0, &d0, &h0, &w0) == 0 && TiffIO_load(a, p2, &u1, &d1, &h1, &w1) == 0 &&
                           d0 == d1 && w0 == w1 && h0 == h1 && memcmp(u0, u1, (size_t)w0 * (size_t)h0) == 0;
                tif = TIFFOpen(p2, "r");
                RT_CHECK(same && tif != NULL && TIFFGetField(tif, TIFFTAG_COMPRESSION, &comp) &&
                         comp == COMPRESSION_ADOBE_DEFLATE, "packed plane: deflate strips, equal to the plain plane");
                if (tif != NULL) TIFFClose(tif);
            }
            snprintf(p2, sizeof p2, "%s/meta.json", seg2);
            fp = fopen(p2, "rb");
            got = 0;
            if (fp != NULL) { got = fread(text, 1, sizeof text - 1, fp); fclose(fp); }
            text[got] = 0;
            RT_CHECK(strstr(text, "\"area_vx2\": 16") != NULL && strstr(text, "\"note\": \"selftest\"") != NULL,
                     "packed meta.json: the plain writer's");
            rt_rm_seg(seg2);
        }
        Arena_dispose(&a);
        rt_rm_seg(seg);
    }

    /* 6. canvas volume -> sheet volume: a 3-layer 10 x 6 canvas in 4 x 4 chunks ('/' separator, one chunk
     * absent) holding value 10k + canvas column (+100 on row 5); pieces [2, 6) reflected and [7, 9) not; the
     * sheet is 9 x 5 and masked at (1, 1).  Sheet (r, t) must hold canvas (r, map(t)) of every layer, 0 at the
     * masked pixel and where the absent chunk (rows 4-5, columns 8-9) was. */
    {
        enum { D = 3, CW = 10, CH = 6, C = 4, SW = 9, SH = 5 };
        char cz[300], oz[300], p[400];
        RawtexTifxyzPiece pc[2];
        uint8_t mask[SW * SH], chunk[D * C * C];
        uint8_t *out = (uint8_t *)malloc(D * 128 * 128);
        size_t written = 0, absent = 0, t = 0, r = 0, cy = 0, cx = 0, k = 0;
        int bad = 0;
        FILE *fp = NULL;
        snprintf(cz, sizeof cz, "%s/canvas", root);
        snprintf(oz, sizeof oz, "%s/sheet.zarr", root);
        {   /* leftovers of an interrupted earlier selftest: the output and its partial */
            const char *leaf[5] = { "0/0.0.0", "0/.zarray", "0", ".zgroup", ".zattrs" };
            int s = 0, l = 0;
            for (s = 0; s < 2; s++) {
                char base[320];
                snprintf(base, sizeof base, "%s%s", oz, s ? ".partial" : "");
                for (l = 0; l < 5; l++) {
                    snprintf(p, sizeof p, "%s/%s", base, leaf[l]);
                    if (l == 2) ves_rmdir(p); else remove(p);
                }
                ves_rmdir(base);
            }
        }
        ves_mkdir(cz);
        snprintf(p, sizeof p, "%s/.zarray", cz);
        fp = fopen(p, "wb");
        if (fp != NULL) {
            fputs("{\"zarr_format\": 2, \"shape\": [3, 6, 10], \"chunks\": [3, 4, 4], \"dtype\": \"|u1\", "
                  "\"compressor\": null, \"fill_value\": 0, \"order\": \"C\", \"filters\": null, "
                  "\"dimension_separator\": \"/\"}", fp);
            fclose(fp);
        }
        snprintf(p, sizeof p, "%s/0", cz); ves_mkdir(p);
        for (cy = 0; cy < 2; cy++) {
            snprintf(p, sizeof p, "%s/0/%zu", cz, cy); ves_mkdir(p);
            for (cx = 0; cx < 3; cx++) {
                if (cy == 1 && cx == 2) continue;               /* absent: fill 0 */
                for (k = 0; k < D; k++)
                    for (r = 0; r < C; r++)
                        for (t = 0; t < C; t++)
                            chunk[(k * C + r) * C + t] = (uint8_t)(10 * k + cx * C + t + (cy * C + r == 5 ? 100 : 0));
                snprintf(p, sizeof p, "%s/0/%zu/%zu", cz, cy, cx);
                fp = fopen(p, "wb");
                if (fp != NULL) { fwrite(chunk, 1, sizeof chunk, fp); fclose(fp); }
            }
        }
        memset(pc, 0, sizeof pc);
        pc[0].col0 = 2; pc[0].col1 = 6; pc[0].reflected = 1;
        pc[1].col0 = 7; pc[1].col1 = 9; pc[1].reflected = 0;
        memset(mask, 255, sizeof mask);
        mask[1 * SW + 1] = 0;
        RT_CHECK(RawtexTifxyz_canvas_to_sheet_zarr(cz, pc, 2, mask, SW, SH, oz, NULL, &written, &absent) == 0 &&
                 written == 1 && absent == 1, "canvas volume mapped");
        snprintf(p, sizeof p, "%s/0/0.0.0", oz);
        fp = out != NULL ? fopen(p, "rb") : NULL;
        RT_CHECK(fp != NULL && fread(out, 1, D * 128 * 128, fp) == D * 128 * 128, "sheet chunk readable");
        if (fp != NULL) fclose(fp);
        for (k = 0; out != NULL && k < D; k++)
            for (r = 0; r < SH; r++)
                for (t = 0; t < SW; t++) {
                    size_t c = (t >= 1 && t < 6) ? 1 + 6 - 1 - t : t;
                    int want = (int)(10 * k + c);
                    if (r == 4 && c >= 8) want = 0;            /* the absent chunk covers rows 4-7, columns 8-11 */
                    if (r == 1 && t == 1) want = 0;            /* masked */
                    if (out[(k * 128 + r) * 128 + t] != want) bad++;
                }
        RT_CHECK(bad == 0, "sheet volume = canvas columns un-reflected, masked, absent chunks zero");
        RT_CHECK(RawtexTifxyz_canvas_to_sheet_zarr(cz, pc, 2, mask, SW, SH, oz, NULL, NULL, NULL) != 0,
                 "refuses to replace a sheet volume");
        /* Standard shards retain the same missing-chunk, reflection, mask and
         * depth-reversal behavior even when the source tile differs from 128. */
        for (size_t sc = 1; sc <= 4; sc *= 2) {
            char zs[400], so[400], file[600];
            uint8_t payload[4 * D * C * C], actual[D * 128 * 128], reverse[SW * SH];
            snprintf(zs, sizeof zs, "%s/shards%zu", root, sc);
            snprintf(so, sizeof so, "%s/shards_sheet%zu", root, sc);
            ves_mkdir(zs);
            RT_CHECK(ZarrShard_array_meta(zs, D, CH, CW, C, sc) == 0, "sharded canvas metadata");
            snprintf(file, sizeof file, "%s/c", zs); ves_mkdir(file);
            snprintf(file, sizeof file, "%s/c/0", zs); ves_mkdir(file);
            for (cy = 0; cy < 2; cy++) {
                snprintf(file, sizeof file, "%s/c/0/%zu", zs, cy); ves_mkdir(file);
                for (size_t sg = 0; sg < (3 + sc - 1) / sc; sg++) {
                    memset(payload, 0, sizeof payload);
                    size_t supplied = 3 - sg * sc < sc ? 3 - sg * sc : sc;
                    for (size_t i = 0; i < supplied; i++) {
                        cx = sg * sc + i;
                        if (cy == 1 && cx == 2) continue;
                        for (k = 0; k < D; k++)
                            for (r = 0; r < C; r++)
                                for (t = 0; t < C; t++)
                                    payload[i * D * C * C + (k * C + r) * C + t] =
                                        (uint8_t)(10 * k + cx * C + t + (cy * C + r == 5 ? 100 : 0));
                    }
                    snprintf(file, sizeof file, "%s/c/0/%zu/%zu", zs, cy, sg);
                    RT_CHECK(ZarrShard_write(file, payload, D * C * C, supplied, sc, NULL) == 0,
                             "write sharded canvas");
                }
            }
            for (t = 0; t < SW * SH; t++) reverse[t] = (uint8_t)(sc == 4 && t % 3 == 0);
            RT_CHECK(RawtexTifxyz_canvas_to_sheet_zarr_depth(zs, pc, 2, mask, reverse, SW, SH, so, NULL,
                                                           &written, &absent) == 0 && written == 1 && absent == 1,
                     "sharded volume mapped with exact missing-chunk count");
            snprintf(file, sizeof file, "%s/0/0.0.0", so);
            fp = fopen(file, "rb");
            RT_CHECK(fp != NULL && fread(actual, 1, sizeof actual, fp) == sizeof actual, "read sharded conversion");
            if (fp != NULL) fclose(fp);
            bad = 0;
            for (k = 0; out != NULL && k < D; k++)
                for (r = 0; r < 128; r++)
                    for (t = 0; t < 128; t++) {
                        size_t depth = r < SH && t < SW && reverse[r * SW + t] ? D - 1 - k : k;
                        bad += actual[(k * 128 + r) * 128 + t] != out[(depth * 128 + r) * 128 + t];
                    }
            RT_CHECK(bad == 0, "sharded conversion equals original including padding and depth reversal");
        }
        free(out);
        for (cy = 0; cy < 2; cy++) {
            for (cx = 0; cx < 3; cx++) { snprintf(p, sizeof p, "%s/0/%zu/%zu", cz, cy, cx); remove(p); }
            snprintf(p, sizeof p, "%s/0/%zu", cz, cy); ves_rmdir(p);
        }
        snprintf(p, sizeof p, "%s/0", cz); ves_rmdir(p);
        snprintf(p, sizeof p, "%s/.zarray", cz); remove(p);
        ves_rmdir(cz);
        snprintf(p, sizeof p, "%s/0/0.0.0", oz); remove(p);
        snprintf(p, sizeof p, "%s/0/.zarray", oz); remove(p);
        snprintf(p, sizeof p, "%s/0", oz); ves_rmdir(p);
        snprintf(p, sizeof p, "%s/.zgroup", oz); remove(p);
        snprintf(p, sizeof p, "%s/.zattrs", oz); remove(p);
        ves_rmdir(oz);
    }

    /* 7. the bake's contrast window: values 0..99 on the mask (one each) and 250 off it; 10/90 percentiles
     * give the window [10, 89]: 10 -> 0, 50 -> 129, 89 and above -> 255, off the mask 0 */
    {
        uint8_t img[101], mask[101];
        int lo = -1, hi = -1, i = 0;
        for (i = 0; i < 100; i++) { img[i] = (uint8_t)i; mask[i] = 1; }
        img[100] = 250; mask[100] = 0;
        RT_CHECK(RawtexTifxyz_stretch(img, mask, 101, 10.0, 90.0, &lo, &hi) == 0 && lo == 10 && hi == 89,
                 "contrast window = the bake's percentile rule");
        RT_CHECK(img[10] == 0 && img[5] == 0 && img[89] == 255 && img[99] == 255 && img[100] == 0 &&
                 img[50] == (uint8_t)floor(255.0 * 40.0 / 79.0 + 0.5), "stretch: clamp, linear, off-mask zero");
    }

    /* 8. the stretch in parts (a sheet in column blocks): a 1000-pixel pseudo-random raster split into blocks
     * of 333/333/334; summed block histograms give the whole raster's window, and the window applied per
     * block gives the whole raster's stretch byte for byte; compose puts the blocks back; bad fits refused */
    {
        uint8_t img[1000], mask[1000], whole[1000], parts[1000], sheet[4 * 250];
        size_t hist[256], k = 0, i = 0;
        int lo = -1, hi = -1, l2 = -2, h2 = -2;
        uint32_t s = 12345u;
        for (i = 0; i < 1000; i++) {
            s = s * 1664525u + 1013904223u;
            img[i] = (uint8_t)(s >> 24);
            mask[i] = (uint8_t)((s >> 8) & 3 ? 1 : 0);
        }
        memcpy(whole, img, sizeof whole);
        memcpy(parts, img, sizeof parts);
        RT_CHECK(RawtexTifxyz_stretch(whole, mask, 1000, 1.0, 99.0, &lo, &hi) == 0, "whole stretch");
        memset(hist, 0, sizeof hist);
        for (k = 0; k < 3; k++) {
            const size_t a = k * 333, n = k == 2 ? 334 : 333;
            RawtexTifxyz_histogram(parts + a, mask + a, n, hist);
        }
        RT_CHECK(RawtexTifxyz_window(hist, 1.0, 99.0, &l2, &h2) == 0 && l2 == lo && h2 == hi,
                 "summed block histograms give the whole raster's window");
        for (k = 0; k < 3; k++) {
            const size_t a = k * 333, n = k == 2 ? 334 : 333;
            RawtexTifxyz_apply_window(parts + a, mask + a, n, l2, h2);
        }
        RT_CHECK(memcmp(parts, whole, sizeof whole) == 0, "the window applied per block == the whole stretch");
        RT_CHECK(RawtexTifxyz_window(hist, 50.0, 40.0, NULL, NULL) != 0, "window refuses lo >= hi percentiles");
        memset(hist, 0, sizeof hist);
        RT_CHECK(RawtexTifxyz_window(hist, 1.0, 99.0, &l2, &h2) == 0 && l2 == 0 && h2 == 255,
                 "an empty histogram gives 0..255");
        /* a 4 x 250 sheet from blocks of 100 + 150 columns, each 4 rows */
        memset(sheet, 0, sizeof sheet);
        RT_CHECK(RawtexTifxyz_compose(sheet, 250, 4, img, 100, 4, 0) == 0 &&
                 RawtexTifxyz_compose(sheet, 250, 4, img + 400, 150, 4, 100) == 0, "compose two blocks");
        for (i = 0, k = 0; i < 4; i++)
            k += memcmp(sheet + i * 250, img + i * 100, 100) != 0 ||
                 memcmp(sheet + i * 250 + 100, img + 400 + i * 150, 150) != 0;
        RT_CHECK(k == 0, "composed rows are the blocks' rows at their offsets");
        RT_CHECK(RawtexTifxyz_compose(sheet, 250, 4, img, 100, 4, 151) != 0 &&
                 RawtexTifxyz_compose(sheet, 250, 4, img, 100, 3, 0) != 0, "compose refuses an overhang or row mismatch");
    }

    /* 5. gray TIFFs: a 300 x 200 ramp through the tiled deflate writer (partial
     * edge tiles) and back; a stripped TiffIO_save raster reads the same */
    {
        const size_t W = 300, H = 200;
        uint8_t *img = (uint8_t *)malloc(W * H), *back = NULL;
        size_t w = 0, h = 0, i = 0;
        char tp[300], sp[300];
        snprintf(tp, sizeof tp, "%s/tiled.tif", root);
        snprintf(sp, sizeof sp, "%s/strips.tif", root);
        RT_CHECK(img != NULL, "ramp allocation");
        for (i = 0; img != NULL && i < W * H; i++) img[i] = (uint8_t)((i % W) * 7 + (i / W) * 3);
        RT_CHECK(img != NULL && RawtexTifxyz_save_gray(tp, img, W, H) == 0, "tiled gray write");
        RT_CHECK(RawtexTifxyz_load_gray(tp, &back, &w, &h) == 0 && w == W && h == H && back != NULL &&
                 memcmp(back, img, W * H) == 0, "tiled gray round trip");
        free(back); back = NULL;
        RT_CHECK(img != NULL && TiffIO_save(sp, img, 1, (int)H, (int)W) == 0 &&
                 RawtexTifxyz_load_gray(sp, &back, &w, &h) == 0 && w == W && h == H && back != NULL &&
                 memcmp(back, img, W * H) == 0, "stripped gray read");
        free(back); free(img);
        remove(tp); remove(sp);
    }

    if (!fails) ves_selftest_clean(root);
    fprintf(stderr, "[selftest] rawtex_tifxyz %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}
