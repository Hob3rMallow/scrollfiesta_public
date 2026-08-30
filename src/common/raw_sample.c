/* raw_sample.c -- see raw_sample.h. Lifted verbatim from obj_bake_raw.c so the
 * RAW read-back tool and the overlap-repair module share one sampler. */

#include "raw_sample.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tiff_io.h"

/* slot[] sentinel: tried and missing (grid edge). */
#define CUBE_ABSENT ((uint8_t *)(uintptr_t)1)

static long floor_div(long a, long b)
{
    long q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) q--;
    return q;
}

static int parse_long3(const char *json, const char *key, long out[3])
{
    const char *p = strstr(json, key);
    if (!p) return -1;
    p = strchr(p, '[');
    if (!p) return -1;
    p++;
    for (int i = 0; i < 3; i++) {
        char *end = NULL;
        while (*p && !((*p >= '0' && *p <= '9') || *p == '-')) p++;
        if (!*p) return -1;
        out[i] = strtol(p, &end, 10);
        if (end == p) return -1;
        p = end;
    }
    return 0;
}

/* Return 1 for a supported raw Zarr, 0 when this is an ordinary TIFF-cube
 * directory, and -1 when a Zarr exists but would require a codec/dtype we do
 * not silently reinterpret. */
static int detect_raw_zarr(CubeTable *ct)
{
    char path[2300];
    char text[4096];
    long chunks[3] = {0,0,0};
    FILE *f;
    size_t n;
    snprintf(path, sizeof path, "%s/0/.zarray", ct->dir);
    f = fopen(path, "rb");
    if (!f) return 0;
    n = fread(text, 1, sizeof text - 1, f);
    fclose(f);
    text[n] = 0;
    const char *comp = strstr(text, "\"compressor\"");
    if (comp) comp = strchr(comp, ':');
    if (comp) { comp++; while (*comp == ' ' || *comp == '\t' || *comp == '\r' || *comp == '\n') comp++; }
    if (!strstr(text, "\"dtype\"") ||
        !(strstr(text, "\"|u1\"") || strstr(text, "\"<u1\"")) ||
        !comp || strncmp(comp, "null", 4) != 0 ||
        parse_long3(text, "\"shape\"", ct->shape) != 0 ||
        parse_long3(text, "\"chunks\"", chunks) != 0 ||
        chunks[0] != ct->chunk || chunks[1] != ct->chunk || chunks[2] != ct->chunk) {
        fprintf(stderr, "ERROR: RAW Zarr must be uncompressed uint8 with %ld^3 chunks: %s\n",
                ct->chunk, path);
        return -1;
    }
    ct->is_zarr = 1;
    return 1;
}

static void cube_path(const CubeTable *ct, long cz, long cy, long cx,
                      char *path, size_t cap)
{
    if (ct->is_zarr)
        snprintf(path, cap, "%s/0/%ld/%ld/%ld", ct->dir, cz, cy, cx);
    else
        snprintf(path, cap, "%s/z%05ld_y%05ld_x%05ld.tif",
                 ct->dir, cz * ct->chunk, cy * ct->chunk, cx * ct->chunk);
}

int cubetable_init(CubeTable *ct, Arena_T arena, const char *dir,
                   long chunk, const float *verts, size_t nv, double pad)
{
    double mn[3] = { 1e30, 1e30, 1e30 }, mx[3] = { -1e30, -1e30, -1e30 };
    size_t i = 0, total = 0, requested = 0;
    long c1[3] = { 0, 0, 0 };
    long requested_dims[3] = { 0, 0, 0 };
    int a = 0;

    memset(ct, 0, sizeof *ct);
    if (nv == 0 || chunk <= 0) return -1;
    ct->arena = arena;
    ct->chunk = chunk;
    snprintf(ct->dir, sizeof ct->dir, "%s", dir);
    if (detect_raw_zarr(ct) < 0) return -1;
    for (i = 0; i < nv; i++) {
        for (a = 0; a < 3; a++) {
            double v = (double)verts[i * 3 + a];
            if (v < mn[a]) mn[a] = v;
            if (v > mx[a]) mx[a] = v;
        }
    }
    ct->cz0 = floor_div((long)floor(mn[0] - pad) - 1, chunk);
    ct->cy0 = floor_div((long)floor(mn[1] - pad) - 1, chunk);
    ct->cx0 = floor_div((long)floor(mn[2] - pad) - 1, chunk);
    c1[0] = floor_div((long)ceil(mx[0] + pad) + 1, chunk);
    c1[1] = floor_div((long)ceil(mx[1] + pad) + 1, chunk);
    c1[2] = floor_div((long)ceil(mx[2] + pad) + 1, chunk);
    requested_dims[0] = c1[0] - ct->cz0 + 1;
    requested_dims[1] = c1[1] - ct->cy0 + 1;
    requested_dims[2] = c1[2] - ct->cx0 + 1;
    if (requested_dims[0] <= 0 || requested_dims[1] <= 0 ||
        requested_dims[2] <= 0 ||
        (size_t)requested_dims[0] > ((size_t)1 << 24) ||
        (size_t)requested_dims[1] > ((size_t)1 << 24) ||
        (size_t)requested_dims[2] > ((size_t)1 << 24) ||
        (size_t)requested_dims[0] > ((size_t)1 << 24) /
                                      (size_t)requested_dims[1] ||
        (size_t)requested_dims[0] * (size_t)requested_dims[1] >
            ((size_t)1 << 24) / (size_t)requested_dims[2]) {
        fprintf(stderr, "ERROR: requested cube table %ldx%ldx%ld unreasonable\n",
                requested_dims[0], requested_dims[1], requested_dims[2]);
        return -1;
    }
    requested = (size_t)requested_dims[0] * (size_t)requested_dims[1] *
                (size_t)requested_dims[2];
    if (ct->is_zarr) {
        long *lo[3] = { &ct->cz0, &ct->cy0, &ct->cx0 };
        for (a = 0; a < 3; a++) {
            long last = (ct->shape[a] - 1) / chunk;
            if (*lo[a] < 0) *lo[a] = 0;
            if (c1[a] > last) c1[a] = last;
        }
        if (c1[0] < ct->cz0 || c1[1] < ct->cy0 || c1[2] < ct->cx0) {
            ct->nz = ct->ny = ct->nx = 0;
            ct->n_outside = requested;
            return 0;
        }
    }
    ct->nz = c1[0] - ct->cz0 + 1;
    ct->ny = c1[1] - ct->cy0 + 1;
    ct->nx = c1[2] - ct->cx0 + 1;
    total = (size_t)ct->nz * (size_t)ct->ny * (size_t)ct->nx;
    if (total == 0 || total > (size_t)1 << 24) {
        fprintf(stderr, "ERROR: cube table %ldx%ldx%ld unreasonable\n",
                ct->nz, ct->ny, ct->nx);
        return -1;
    }
    ct->n_outside = requested - total;
    ct->slot = (uint8_t **)ARENA_CALLOC(arena, (long)total,
                                        (long)sizeof(uint8_t *));
    return 0;
}

int cube_fetch(CubeTable *ct, long iz, long iy, long ix)
{
    long cz = floor_div(iz, ct->chunk);
    long cy = floor_div(iy, ct->chunk);
    long cx = floor_div(ix, ct->chunk);
    long tz = cz - ct->cz0, ty = cy - ct->cy0, tx = cx - ct->cx0;
    size_t si = 0;
    uint8_t *buf = NULL;

    if (ct->is_zarr && (iz < 0 || iy < 0 || ix < 0 ||
                        iz >= ct->shape[0] || iy >= ct->shape[1] ||
                        ix >= ct->shape[2])) return -1;
    if (tz < 0 || tz >= ct->nz || ty < 0 || ty >= ct->ny
        || tx < 0 || tx >= ct->nx) return -1;
    si = ((size_t)tz * (size_t)ct->ny + (size_t)ty) * (size_t)ct->nx
         + (size_t)tx;
    buf = ct->slot[si];
    if (buf == NULL) {   /* lazy load */
        char path[2600];
        uint8_t *vol = NULL;
        int D = 0, H = 0, W = 0;
        size_t chunk_bytes = (size_t)ct->chunk * (size_t)ct->chunk * (size_t)ct->chunk;
        cube_path(ct, cz, cy, cx, path, sizeof path);
        /* Missing halo cubes are expected at a cropped volume boundary.
         * Probe quietly before libtiff: TIFFOpen otherwise emits one scary
         * diagnostic per absent neighbor (hundreds in a normal full run). */
        FILE *probe = NULL;
        if (ct->is_zarr && (iz < 0 || iy < 0 || ix < 0 ||
                            iz >= ct->shape[0] || iy >= ct->shape[1] || ix >= ct->shape[2])) {
            probe = NULL;
        } else {
            probe = fopen(path, "rb");
        }
        if (probe == NULL) {
            buf = CUBE_ABSENT;
            ct->n_missing++;
        } else {
            if (ct->is_zarr) {
                vol = (uint8_t *)ARENA_ALLOC(ct->arena, (long)chunk_bytes);
                size_t got = fread(vol, 1, chunk_bytes, probe);
                fclose(probe);
                if (got == chunk_bytes) {
                    buf = vol;
                    ct->n_loaded++;
                } else {
                    buf = CUBE_ABSENT;
                    ct->n_missing++;
                }
            } else {
                fclose(probe);
                if (TiffIO_load(ct->arena, path, &vol, &D, &H, &W) == 0
                    && D == (int)ct->chunk && H == (int)ct->chunk
                    && W == (int)ct->chunk) {
                    buf = vol;
                    ct->n_loaded++;
                } else {
                    buf = CUBE_ABSENT;
                    ct->n_missing++;
                }
            }
        }
        ct->slot[si] = buf;
    }
    if (buf == CUBE_ABSENT) return -1;
    return (int)buf[((size_t)(iz - cz * ct->chunk) * (size_t)ct->chunk
                     + (size_t)(iy - cy * ct->chunk)) * (size_t)ct->chunk
                    + (size_t)(ix - cx * ct->chunk)];
}

int cubetable_prewarm_all(CubeTable *ct)
{
    long before = ct->n_loaded;
    long total = ct->nz * ct->ny * ct->nx;
#ifdef _OPENMP
    /* Phase 1: parallel page-cache warm.  The serial TIFF parse below is
     * memcpy-speed once the bytes are in the OS cache; reading ~2k cube
     * files one at a time was 80+ s of pure disk latency (the whole bake
     * cost on the 4x21x21).  Neither libtiff nor the arena is touched here,
     * so this phase is trivially thread-safe. */
    {
        long t = 0;
#pragma omp parallel for schedule(dynamic, 4)
        for (t = 0; t < total; t++) {
            long tz = t / (ct->ny * ct->nx);
            long ty = (t / ct->nx) % ct->ny;
            long tx = t % ct->nx;
            char path[2600];
            FILE *f = NULL;
            char buf[32768];
            cube_path(ct, ct->cz0 + tz, ct->cy0 + ty, ct->cx0 + tx,
                      path, sizeof path);
            f = fopen(path, "rb");
            if (f == NULL) continue;
            while (fread(buf, 1, sizeof buf, f) == sizeof buf) { /* warm */ }
            fclose(f);
        }
    }
#endif
    for (long tz = 0; tz < ct->nz; tz++)
        for (long ty = 0; ty < ct->ny; ty++)
            for (long tx = 0; tx < ct->nx; tx++)
                (void)cube_fetch(ct, (ct->cz0 + tz) * ct->chunk,
                                 (ct->cy0 + ty) * ct->chunk,
                                 (ct->cx0 + tx) * ct->chunk);
    (void)total;
    return (int)(ct->n_loaded - before);
}

size_t cubetable_expected_chunks(const CubeTable *ct)
{
    if (ct == NULL || ct->nz <= 0 || ct->ny <= 0 || ct->nx <= 0) return 0;
    return (size_t)ct->nz * (size_t)ct->ny * (size_t)ct->nx;
}

int cubetable_is_complete(const CubeTable *ct)
{
    size_t expected = cubetable_expected_chunks(ct);
    return ct != NULL && ct->n_missing == 0 &&
           (size_t)ct->n_loaded == expected;
}

double sample_trilinear(CubeTable *ct, double z, double y, double x)
{
    double fz = floor(z), fy = floor(y), fx = floor(x);
    long iz = (long)fz, iy = (long)fy, ix = (long)fx;
    double dz = z - fz, dy = y - fy, dx = x - fx;
    double acc = 0.0, wsum = 0.0;
    int k = 0;

    for (k = 0; k < 8; k++) {
        long oz = k & 1, oy = (k >> 1) & 1, ox = (k >> 2) & 1;
        int v = cube_fetch(ct, iz + oz, iy + oy, ix + ox);
        double w = 0.0;
        if (v < 0) continue;
        w = (oz ? dz : 1.0 - dz) * (oy ? dy : 1.0 - dy)
            * (ox ? dx : 1.0 - dx);
        acc += w * (double)v;
        wsum += w;
    }
    if (wsum < 1e-9) return -1.0;
    return acc / wsum;
}

int sample_tangent_tensor(CubeTable *ct, const double p[3],
                          const double tu_in[3], const double tv_in[3],
                          double radius, RawTangentTensor *out)
{
    double tu[3], tv[3], ul = 0.0, vl = 0.0, dot = 0.0;
    double img[25];
    uint8_t has[25];
    double Sxx = 0.0, Sxy = 0.0, Syy = 0.0, wsum = 0.0;
    double step;
    int i, j, ngrad = 0;

    if (ct == NULL || p == NULL || tu_in == NULL || tv_in == NULL || out == NULL)
        return -1;
    memset(out, 0, sizeof *out);
    if (radius <= 0.0) radius = 2.0;
    step = radius * 0.5;

    for (i = 0; i < 3; i++) ul += tu_in[i] * tu_in[i];
    ul = sqrt(ul);
    if (ul < 1e-9) return -1;
    for (i = 0; i < 3; i++) tu[i] = tu_in[i] / ul;
    for (i = 0; i < 3; i++) dot += tv_in[i] * tu[i];
    for (i = 0; i < 3; i++) {
        tv[i] = tv_in[i] - dot * tu[i];
        vl += tv[i] * tv[i];
    }
    vl = sqrt(vl);
    if (vl < 1e-9) return -1;
    for (i = 0; i < 3; i++) tv[i] /= vl;

    for (j = -2; j <= 2; j++) {
        for (i = -2; i <= 2; i++) {
            double q[3];
            double s;
            int k, pi = (j + 2) * 5 + (i + 2);
            for (k = 0; k < 3; k++)
                q[k] = p[k] + step * ((double)i * tu[k] + (double)j * tv[k]);
            s = sample_trilinear(ct, q[0], q[1], q[2]);
            img[pi] = s;
            has[pi] = (s >= 0.0) ? (uint8_t)1 : (uint8_t)0;
        }
    }

    /* Tensor integration over the inner 3x3. A light binomial weight makes the
     * score stable without paying for another Gaussian sampling halo. */
    for (j = 1; j <= 3; j++) {
        for (i = 1; i <= 3; i++) {
            int p0 = j * 5 + i;
            int pl = p0 - 1, pr = p0 + 1, pu = p0 - 5, pd = p0 + 5;
            double gx, gy, w;
            if (!has[pl] || !has[pr] || !has[pu] || !has[pd]) continue;
            gx = (img[pr] - img[pl]) / (2.0 * step);
            gy = (img[pd] - img[pu]) / (2.0 * step);
            w = (i == 2 ? 2.0 : 1.0) * (j == 2 ? 2.0 : 1.0);
            Sxx += w * gx * gx;
            Sxy += w * gx * gy;
            Syy += w * gy * gy;
            wsum += w;
            ngrad++;
        }
    }
    if (ngrad < 5 || wsum <= 0.0) return -1;
    Sxx /= wsum; Sxy /= wsum; Syy /= wsum;
    {
        double E = Sxx + Syy;
        double A = Sxx - Syy, B = 2.0 * Sxy;
        double mag = sqrt(A * A + B * B);
        double coh = E > 1e-12 ? mag / E : 0.0;
        double axis = 0.5;
        double rms = sqrt(E > 0.0 ? E : 0.0);
        double strength = rms / (rms + 8.0); /* smooth saturation in CT units */
        double tn[3] = { tu[1] * tv[2] - tu[2] * tv[1],
                         tu[2] * tv[0] - tu[0] * tv[2],
                         tu[0] * tv[1] - tu[1] * tv[0] };
        double qm[3], qp[3], sm, sp, asym = 255.0, ridge = 0.0;
        for (i = 0; i < 3; i++) {
            qm[i] = p[i] - 0.75 * tn[i];
            qp[i] = p[i] + 0.75 * tn[i];
        }
        sm = sample_trilinear(ct, qm[0], qm[1], qm[2]);
        sp = sample_trilinear(ct, qp[0], qp[1], qp[2]);
        if (sm >= 0.0 && sp >= 0.0) {
            double a;
            asym = fabs(sp - sm);
            a = asym / 16.0;
            ridge = 1.0 / (1.0 + a * a);
        }
        if (mag > 1e-12) {
            double cos4 = (A * A - B * B) / (mag * mag);
            axis = 0.5 * (1.0 + cos4);
        }
        if (coh < 0.0) coh = 0.0; else if (coh > 1.0) coh = 1.0;
        if (axis < 0.0) axis = 0.0; else if (axis > 1.0) axis = 1.0;
        out->energy = E;
        out->rms_gradient = rms;
        out->coherence = coh;
        out->axis_alignment = axis;
        out->normal_asymmetry = asym;
        out->ridge_center = ridge;
        /* Energy keeps balanced cross-hatch alive; coherence/alignment promote
         * a clean axial/circumferential fiber family over arbitrary texture.
         * The ridge factor is essential: without it, a tangent frame tilted a
         * few degrees at a sheet shoulder aliases the much stronger normal
         * bright/dark transition into a convincing but fake tangent tensor. */
        out->quality = ridge * strength * (0.25 + 0.75 * coh * axis);
        if (out->quality < 0.0) out->quality = 0.0;
        if (out->quality > 1.0) out->quality = 1.0;
    }
    return 0;
}

double sample_vertex(CubeTable *ct, const float *p, const float *n,
                     double range, int nsteps)
{
    double nn = 0.0, best = -1.0;
    int k = 0;

    if (n != NULL)
        nn = sqrt((double)n[0] * (double)n[0] + (double)n[1] * (double)n[1]
                  + (double)n[2] * (double)n[2]);
    if (range <= 0.0 || nsteps < 2 || nn < 0.5)
        return sample_trilinear(ct, (double)p[0], (double)p[1], (double)p[2]);
    for (k = 0; k < nsteps; k++) {
        double t = -range + (2.0 * range * (double)k) / (double)(nsteps - 1);
        double s = sample_trilinear(ct,
                                    (double)p[0] + t * (double)n[0],
                                    (double)p[1] + t * (double)n[1],
                                    (double)p[2] + t * (double)n[2]);
        if (s > best) best = s;
    }
    return best;
}
