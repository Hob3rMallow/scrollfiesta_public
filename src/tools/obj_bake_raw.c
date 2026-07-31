/* ============================================================================
 * obj_bake_raw -- bake RAW CT intensity onto a UV-carrying OBJ.
 *
 * Standalone validation tool for a UV-carrying scroll mesh: "read the original volumetric
 * information back from the scroll and write it out to an OBJ -- if this
 * works, we should see nice, clean papyrus."
 *
 * Input: an OBJ whose vertices are SOURCE-SPACE voxel coordinates (z,y,x)
 * with per-vertex "vt u v" texture coordinates (as written by the unroll tools:
 * <id>_uv.obj or <id>_ribbon.obj), plus the grid's cubes_RAW/ directory of
 * chunk^3 uint8 multipage TIFFs named z#####_y#####_x#####.tif (the name is
 * the cube's source-space voxel origin, so it doubles as the placement).
 *
 * Per vertex: trilinear-sample the RAW volume, taking the MAX over a few
 * steps along +/- the vertex normal (default +/-2 vox, 5 steps -- below half
 * the 7-vox inter-wrap clearance, so a sample can never read the NEIGHBORING
 * wrap; the max picks up the bright papyrus core even when the recto boundary
 * vertex sits half in air). Intensities are contrast-stretched (percentile
 * window over the sampled population) and emitted as grayscale per-vertex
 * colors:
 *
 *   <out>/<id>_raw3d.obj    original 3D positions + colors
 *   <out>/<id>_rawflat.obj  vertices at (u, v, 0) + colors -- the unrolled
 *                           papyrus; render with mesh_render --vcolor
 *   <out>/<id>_rawtex.tif   RASTERIZED unrolled texture (default 1 px/vox):
 *                           every pixel center is located in its covering UV
 *                           triangle, the 3D point + normal barycentrically
 *                           interpolated, and the volume freshly sampled
 *                           there. Never a per-vertex scatter.
 *
 * Missing cubes (grid edge) degrade gracefully: trilinear weights are
 * renormalized over the corners that exist; vertices with no support at all
 * are counted and colored black.
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include "../common/raw_sample.h"
#include "../flatten/rawtex_bake.h"   /* Rawtex_write_tif + DiagOpts (extracted) */

/* ------------------------------------------------------------------ util */

static void *xmalloc(size_t nbytes)
{
    void *p = malloc(nbytes > 0 ? nbytes : 1);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu bytes)\n", nbytes);
        exit(1);
    }
    return p;
}

static void *xcalloc(size_t count, size_t size)
{
    void *p = calloc(count > 0 ? count : 1, size);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu x %zu)\n", count, size);
        exit(1);
    }
    return p;
}

static void *xgrow(void *p, size_t need_elems, size_t *cap_elems, size_t elem)
{
    size_t ncap = *cap_elems;
    if (need_elems <= ncap) return p;
    ncap = (ncap == 0) ? 4096 : ncap;
    while (ncap < need_elems) ncap *= 2;
    p = realloc(p, ncap * elem);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu bytes)\n", ncap * elem);
        exit(1);
    }
    *cap_elems = ncap;
    return p;
}

/* -------------------------------------------------- OBJ with vt parsing */

typedef struct UvMesh {
    float   *verts;   /* [nv*3] (z,y,x) as stored in the file */
    float   *uv;      /* [nvt*2] */
    int32_t *faces;   /* [nf*3] 0-based */
    size_t   nv, nvt, nf;
} UvMesh;

/* Parse v / vt / f lines. Face corners "a", "a/b" and "a/b/c" all accepted
 * (only the vertex index is kept; the unroll tools write a/a so vt index ==
 * v index by construction). Polygons are fan-triangulated. Vertex color
 * floats after v x y z are ignored. Returns 0 on success. */
static int parse_uv_obj(const char *path, UvMesh *m)
{
    FILE *f = NULL;
    char line[1024];
    size_t cap_v = 0, cap_t = 0, cap_f = 0;

    memset(m, 0, sizeof *m);
    f = fopen(path, "rb");
    if (f == NULL) {
        fprintf(stderr, "ERROR: cannot open %s\n", path);
        return -1;
    }
    while (fgets(line, sizeof line, f) != NULL) {
        if (line[0] == 'v' && line[1] == ' ') {
            char *s = line + 2, *end = NULL;
            double a = strtod(s, &end); s = end;
            double b = strtod(s, &end); s = end;
            double c = strtod(s, &end);
            m->verts = (float *)xgrow(m->verts, (m->nv + 1) * 3, &cap_v,
                                      sizeof(float));
            m->verts[m->nv * 3 + 0] = (float)a;
            m->verts[m->nv * 3 + 1] = (float)b;
            m->verts[m->nv * 3 + 2] = (float)c;
            m->nv++;
        } else if (line[0] == 'v' && line[1] == 't' && line[2] == ' ') {
            char *s = line + 3, *end = NULL;
            double u = strtod(s, &end); s = end;
            double v = strtod(s, &end);
            m->uv = (float *)xgrow(m->uv, (m->nvt + 1) * 2, &cap_t,
                                   sizeof(float));
            m->uv[m->nvt * 2 + 0] = (float)u;
            m->uv[m->nvt * 2 + 1] = (float)v;
            m->nvt++;
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char *s = line + 2;
            long idx[16];
            int nidx = 0, k = 0;
            while (nidx < 16) {
                char *end = NULL;
                long v = 0;
                while (*s == ' ' || *s == '\t') s++;
                if (*s == '\0' || *s == '\n' || *s == '\r') break;
                v = strtol(s, &end, 10);
                if (end == s) break;
                idx[nidx++] = v;
                s = end;
                while (*s != '\0' && *s != ' ' && *s != '\t'
                       && *s != '\n' && *s != '\r') s++;  /* skip /vt/vn */
            }
            for (k = 2; k < nidx; k++) {
                m->faces = (int32_t *)xgrow(m->faces, (m->nf + 1) * 3, &cap_f,
                                            sizeof(int32_t));
                m->faces[m->nf * 3 + 0] = (int32_t)(idx[0] - 1);
                m->faces[m->nf * 3 + 1] = (int32_t)(idx[k - 1] - 1);
                m->faces[m->nf * 3 + 2] = (int32_t)(idx[k] - 1);
                m->nf++;
            }
        }
    }
    fclose(f);
    if (m->nv == 0) {
        fprintf(stderr, "ERROR: %s has no vertices\n", path);
        return -1;
    }
    {   /* validate face indices (negative/relative indices unsupported) */
        size_t t = 0, bad = 0;
        for (t = 0; t < m->nf * 3; t++)
            if (m->faces[t] < 0 || (size_t)m->faces[t] >= m->nv) bad++;
        if (bad > 0) {
            fprintf(stderr, "ERROR: %zu face indices out of range in %s\n",
                    bad, path);
            return -1;
        }
    }
    return 0;
}

static void uvmesh_free(UvMesh *m)
{
    free(m->verts); free(m->uv); free(m->faces);
    memset(m, 0, sizeof *m);
}

/* cube table + RAW sampling now live in src/common/raw_sample.{c,h}
 * (CubeTable, cubetable_init, cube_fetch, sample_trilinear, sample_vertex,
 * vertex_normals) -- shared with the overlap-repair module. */

/* --------------------------------------------------- stretch + writers */

/* Percentile window over the sampled population (256-bin histogram). */
static void stretch_window(const double *val, const uint8_t *has, size_t nv,
                           double pct_lo, double pct_hi,
                           double *out_lo, double *out_hi)
{
    size_t hist[256];
    size_t i = 0, n = 0, cum = 0, tlo = 0, thi = 0;
    int b = 0, lo = 0, hi = 255;

    memset(hist, 0, sizeof hist);
    for (i = 0; i < nv; i++) {
        int bin = 0;
        if (!has[i]) continue;
        bin = (int)(val[i] + 0.5);
        if (bin < 0) bin = 0;
        if (bin > 255) bin = 255;
        hist[bin]++;
        n++;
    }
    *out_lo = 0.0; *out_hi = 255.0;
    if (n == 0) return;
    tlo = (size_t)(pct_lo / 100.0 * (double)n);
    thi = (size_t)(pct_hi / 100.0 * (double)n);
    cum = 0;
    for (b = 0; b < 256; b++) {
        cum += hist[b];
        if (cum > tlo) { lo = b; break; }
    }
    cum = 0;
    for (b = 0; b < 256; b++) {
        cum += hist[b];
        if (cum >= thi) { hi = b; break; }
    }
    if (hi <= lo) { lo = 0; hi = 255; }
    *out_lo = (double)lo;
    *out_hi = (double)hi;
}

static double gray_of(double v, double lo, double hi)
{
    double g = (v - lo) / (hi - lo);
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    return g;
}

/* Colored OBJ: "v x y z r g b" (grayscale). flat != 0 places vertices at
 * (u, v, 0) from uv instead of their 3D position. Unsampled verts get 0. */
static int write_colored_obj(const char *path, const float *verts, size_t nv,
                             const int32_t *faces, size_t nf, const float *uv,
                             const double *val, const uint8_t *has,
                             double lo, double hi, int flat)
{
    FILE *f = fopen(path, "wb");
    size_t i = 0, t = 0;

    if (f == NULL) return -1;
    for (i = 0; i < nv; i++) {
        double g = has[i] ? gray_of(val[i], lo, hi) : 0.0;
        if (flat)
            fprintf(f, "v %.4f %.4f 0 %.4f %.4f %.4f\n",
                    (double)uv[i * 2 + 0], (double)uv[i * 2 + 1], g, g, g);
        else
            fprintf(f, "v %.4f %.4f %.4f %.4f %.4f %.4f\n",
                    (double)verts[i * 3 + 0], (double)verts[i * 3 + 1],
                    (double)verts[i * 3 + 2], g, g, g);
    }
    for (t = 0; t < nf; t++)
        fprintf(f, "f %d %d %d\n", faces[t * 3 + 0] + 1,
                faces[t * 3 + 1] + 1, faces[t * 3 + 2] + 1);
    fclose(f);
    return 0;
}


/* -------------------------------------------------------------- selftest */

#define CHECK(cond, msg) \
    do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", (msg)); nfail++; } \
         else { fprintf(stderr, "  ok: %s\n", (msg)); } } while (0)

static int selftest(void)
{
    int nfail = 0;
    const char *dir = "output/obj_bake_raw_selftest";
    Arena_T arena = Arena_new();

    fprintf(stderr, "[selftest] obj_bake_raw\n");

    /* 1. OBJ parse: v/vt/f with a/a corners + quad fan */
    {
        char path[512];
        FILE *f = NULL;
        UvMesh m;
        snprintf(path, sizeof path, "%s/t1.obj", dir);
        ves_ensure_parent_dir(path);
        f = fopen(path, "wb");
        if (f != NULL) {
            fprintf(f, "# hdr\nv 1 2 3\nv 4 5 6\nv 7 8 9\nv 10 11 12\n"
                       "vt 0.5 1.5\nvt 2.5 3.5\nvt 4.5 5.5\nvt 6.5 7.5\n"
                       "f 1/1 2/2 3/3\nf 1/1 2/2 3/3 4/4\n");
            fclose(f);
        }
        CHECK(parse_uv_obj(path, &m) == 0, "t1 parse ok");
        CHECK(m.nv == 4 && m.nvt == 4 && m.nf == 3,
              "t1 counts (4v 4vt 1+2f fan)");
        CHECK(fabs((double)m.verts[3 * 3 + 2] - 12.0) < 1e-6
              && fabs((double)m.uv[2 * 2 + 0] - 4.5) < 1e-6,
              "t1 values exact");
        CHECK(m.faces[2 * 3 + 0] == 0 && m.faces[2 * 3 + 1] == 2
              && m.faces[2 * 3 + 2] == 3, "t1 fan triangulation");
        uvmesh_free(&m);
    }

    /* 2. trilinear exact on a linear ramp (injected cube, no file) */
    {
        CubeTable ct;
        long chunk = 8;
        uint8_t *buf = (uint8_t *)ARENA_ALLOC(arena, chunk * chunk * chunk);
        long z = 0, y = 0, x = 0;
        float fake[6] = { 1.0f, 1.0f, 1.0f, 6.0f, 6.0f, 6.0f };
        double s = 0.0;
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++)
                    buf[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)(z + 2 * y + 3 * x);
        CHECK(cubetable_init(&ct, arena, dir, chunk, fake, 2, 0.0) == 0,
              "t2 table init");
        ct.slot[0] = buf;   /* inject: table spans exactly cube (0,0,0) */
        CHECK(ct.nz == 1 && ct.ny == 1 && ct.nx == 1, "t2 table dims 1x1x1");
        s = sample_trilinear(&ct, 1.5, 2.25, 3.75);
        CHECK(fabs(s - (1.5 + 4.5 + 11.25)) < 1e-9, "t2 trilinear exact");
        s = sample_trilinear(&ct, 3.0, 4.0, 5.0);
        CHECK(fabs(s - 26.0) < 1e-9, "t2 integer point exact");
    }

    /* 3. cross-cube continuity + absent-corner renormalization */
    {
        CubeTable ct;
        long chunk = 4;
        uint8_t *b0 = (uint8_t *)ARENA_ALLOC(arena, chunk * chunk * chunk);
        uint8_t *b1 = (uint8_t *)ARENA_ALLOC(arena, chunk * chunk * chunk);
        long z = 0, y = 0, x = 0;
        float fake[6] = { 1.0f, 1.0f, 1.0f, 2.0f, 2.0f, 6.0f };
        double s = 0.0;
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++) {
                    b0[(size_t)((z * chunk + y) * chunk + x)] = (uint8_t)x;
                    b1[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)(x + chunk);
                }
        CHECK(cubetable_init(&ct, arena, dir, chunk, fake, 2, 0.0) == 0
              && ct.nx == 2 && ct.nz == 1 && ct.ny == 1,
              "t3 table 1x1x2");
        ct.slot[0] = b0;
        ct.slot[1] = b1;
        s = sample_trilinear(&ct, 2.0, 2.0, 3.5);   /* straddles the seam */
        CHECK(fabs(s - 3.5) < 1e-9, "t3 cross-cube seam exact");
        s = sample_trilinear(&ct, 2.0, 2.0, 6.25);
        CHECK(fabs(s - 6.25) < 1e-9, "t3 second cube exact");
        /* corner beyond x=7 leaves the table: renormalize to x=7 plane */
        s = sample_trilinear(&ct, 2.0, 2.0, 7.4);
        CHECK(fabs(s - 7.0) < 1e-9, "t3 renormalized at volume edge");
        s = sample_trilinear(&ct, 2.0, 2.0, -3.0);
        CHECK(s < 0.0, "t3 fully outside -> -1");
    }

    /* 4. normal-max finds the bright slab */
    {
        CubeTable ct;
        long chunk = 8;
        uint8_t *buf = (uint8_t *)ARENA_ALLOC(arena, chunk * chunk * chunk);
        long z = 0, y = 0, x = 0;
        float fake[6] = { 1.0f, 1.0f, 1.0f, 6.0f, 6.0f, 6.0f };
        float p[3] = { 2.5f, 4.0f, 4.0f };
        float n[3] = { 1.0f, 0.0f, 0.0f };  /* +z in (z,y,x) */
        double s = 0.0;
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++)
                    buf[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)((z == 4 || z == 5) ? 200 : 20);
        cubetable_init(&ct, arena, dir, chunk, fake, 2, 0.0);
        ct.slot[0] = buf;
        s = sample_vertex(&ct, p, n, 2.0, 5);
        CHECK(fabs(s - 200.0) < 1e-9, "t4 normal-max hits slab");
        s = sample_vertex(&ct, p, n, 0.0, 5);   /* range 0 -> point sample */
        CHECK(fabs(s - 20.0) < 1e-9, "t4 point sample misses slab");
        {
            float zn[3] = { 0.0f, 0.0f, 0.0f };  /* degenerate normal */
            s = sample_vertex(&ct, p, zn, 2.0, 5);
            CHECK(fabs(s - 20.0) < 1e-9, "t4 zero normal -> point fallback");
        }
    }

    /* 5. lazy TIFF loading through the real path */
    {
        CubeTable ct;
        long chunk = 8;
        char cdir[512], path[600];
        uint8_t *buf = (uint8_t *)xmalloc((size_t)(chunk * chunk * chunk));
        long z = 0, y = 0, x = 0;
        float fake[6] = { 1.0f, 1.0f, 1.0f, 6.0f, 6.0f, 14.0f };
        int v = 0;
        snprintf(cdir, sizeof cdir, "%s/cubes", dir);
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++)
                    buf[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)(10 + x);
        snprintf(path, sizeof path, "%s/z00000_y00000_x00000.tif", cdir);
        ves_ensure_parent_dir(path);
        CHECK(TiffIO_save(path, buf, (int)chunk, (int)chunk, (int)chunk) == 0,
              "t5 save cube A");
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++)
                    buf[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)(100 + x);
        snprintf(path, sizeof path, "%s/z00000_y00000_x00008.tif", cdir);
        CHECK(TiffIO_save(path, buf, (int)chunk, (int)chunk, (int)chunk) == 0,
              "t5 save cube B");
        free(buf);
        CHECK(cubetable_init(&ct, arena, cdir, chunk, fake, 2, 0.0) == 0,
              "t5 table init");
        v = cube_fetch(&ct, 2, 2, 3);
        CHECK(v == 13, "t5 fetch cube A voxel");
        v = cube_fetch(&ct, 2, 2, 11);
        CHECK(v == 103, "t5 fetch cube B voxel (lazy load)");
        CHECK(ct.n_loaded == 2, "t5 two cubes loaded");
        v = cube_fetch(&ct, 2, 9, 3);   /* y=9 -> cube y0=8: no file */
        (void)v;
        CHECK(ct.n_missing >= 0 && cube_fetch(&ct, 2, 9, 3) == -1,
              "t5 missing neighbor -> -1");
    }

    /* 6. stretch window percentiles */
    {
        double val[100];
        uint8_t has[100];
        double lo = 0.0, hi = 0.0;
        int i = 0;
        for (i = 0; i < 100; i++) { val[i] = (double)(i + 78); has[i] = 1; }
        /* values 78..177; 2%/98% -> ~80 / ~175 */
        stretch_window(val, has, 100, 2.0, 98.0, &lo, &hi);
        CHECK(lo >= 79.0 && lo <= 81.0 && hi >= 174.0 && hi <= 176.0,
              "t6 percentile window");
        memset(has, 0, sizeof has);
        stretch_window(val, has, 100, 2.0, 98.0, &lo, &hi);
        CHECK(lo == 0.0 && hi == 255.0, "t6 empty -> full window");
    }

    /* 7. colored OBJ writer round-trip */
    {
        char path[512];
        float verts[12] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0 };
        float uv[8] = { 0, 0, 2, 0, 0, 2, 2, 2 };
        int32_t faces[3] = { 0, 1, 2 };
        double val[4] = { 0.0, 100.0, 200.0, 255.0 };
        uint8_t has[4] = { 1, 1, 1, 1 };
        snprintf(path, sizeof path, "%s/t7_flat.obj", dir);
        CHECK(write_colored_obj(path, verts, 4, faces, 1, uv, val, has,
                                0.0, 255.0, 1) == 0, "t7 write flat obj");
        {   /* re-parse: 6-float v lines must still read as 3 coords */
            UvMesh m;
            CHECK(parse_uv_obj(path, &m) == 0 && m.nv == 4 && m.nf == 1,
                  "t7 colored obj re-parses");
            CHECK(fabs((double)m.verts[1 * 3 + 0] - 2.0) < 1e-6
                  && fabs((double)m.verts[1 * 3 + 1] - 0.0) < 1e-6,
                  "t7 flat coords are (u,v,0)");
            uvmesh_free(&m);
        }
    }

    /* 8. rasterized texture: exact on a linear field, degenerate faces skip */
    {
        CubeTable ct;
        long chunk = 8;
        uint8_t *buf = (uint8_t *)ARENA_ALLOC(arena, chunk * chunk * chunk);
        long z = 0, y = 0, x = 0;
        float fake[6] = { 1.0f, 1.0f, 1.0f, 6.0f, 6.0f, 6.0f };
        /* square sheet on plane z=2: 3D (2, 1+v, 1+u), uv spans [0,4]^2 */
        float verts[12] = { 2, 1, 1,  2, 1, 5,  2, 5, 1,  2, 5, 5 };
        float uvq[8] = { 0, 0,  4, 0,  0, 4,  4, 4 };
        int32_t fq[6] = { 0, 1, 2,  1, 3, 2 };
        char path[512];
        size_t W = 0, H = 0, multi = 0;
        double fill = 0.0;
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++)
                    buf[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)(z + 2 * y + 3 * x);
        cubetable_init(&ct, arena, dir, chunk, fake, 2, 0.0);
        ct.slot[0] = buf;
        snprintf(path, sizeof path, "%s/t8_raster.tif", dir);
        ves_ensure_parent_dir(path);
        CHECK(Rawtex_write_tif(path, &ct, verts, uvq, 4, fq, 2, NULL, NULL,
                               0.0, 1, 1.0, 1.0, 0.0, 255.0,
                               4.0, 25.0, 6.0, NULL,
                               &W, &H, &fill, &multi, NULL, NULL) == 0,
              "t8 raster runs");
        CHECK(W == 5 && H == 5 && fabs(fill - 1.0) < 1e-9,
              "t8 raster 5x5, fully covered");
        {   /* field(2, 1+pv, 1+pu) = 2 + 2(1+pv) + 3(1+pu) = 7 + 2pv + 3pu */
            uint8_t *img = NULL;
            int D = 0, IH = 0, IW = 0;
            int okpix = 1, py = 0, pxx = 0;
            CHECK(TiffIO_load(arena, path, &img, &D, &IH, &IW) == 0
                  && D == 1 && IH == 5 && IW == 5, "t8 raster loads back");
            for (py = 0; py < 5 && img != NULL; py++)
                for (pxx = 0; pxx < 5; pxx++)
                    if (img[py * 5 + pxx]
                        != (uint8_t)(7 + 2 * py + 3 * pxx)) okpix = 0;
            CHECK(okpix, "t8 every pixel exact on linear field");
            CHECK(multi <= 5, "t8 multi-cover only on the shared diagonal");
        }
        {   /* degenerate UV face (zero area) is skipped cleanly */
            int32_t fdeg[3] = { 0, 0, 1 };
            size_t W2 = 0, H2 = 0, m2 = 0;
            double f2 = 0.0;
            snprintf(path, sizeof path, "%s/t8_degen.tif", dir);
            CHECK(Rawtex_write_tif(path, &ct, verts, uvq, 4, fdeg, 1, NULL, NULL,
                                   0.0, 1, 1.0, 1.0, 0.0, 255.0,
                                   4.0, 25.0, 6.0, NULL,
                                   &W2, &H2, &f2, &m2, NULL, NULL) == 0
                  && f2 == 0.0,
                  "t8 degenerate face skipped, empty raster");
        }
        {   /* UV-STRETCH gate: same 3D square, uv stretched 20x in u */
            float uvs[8] = { 0, 0,  80, 0,  0, 4,  80, 4 };
            size_t W2 = 0, H2 = 0, m2 = 0, suv = 0, s3d = 0;
            double f2 = 0.0;
            snprintf(path, sizeof path, "%s/t8_stretch.tif", dir);
            CHECK(Rawtex_write_tif(path, &ct, verts, uvs, 4, fq, 2, NULL, NULL,
                                   0.0, 1, 1.0, 1.0, 0.0, 255.0,
                                   4.0, 25.0, 6.0, NULL,
                                   &W2, &H2, &f2, &m2, &suv, &s3d) == 0
                  && f2 == 0.0 && suv == 2 && s3d == 0,
                  "t8 uv-stretched faces skipped (paint nothing)");
        }
        {   /* MAX-EDGE-3D gate: 10-vox 3D edges (a fill membrane) */
            float vbig[12] = { 2, 1, 1,  2, 1, 11,  2, 5, 1,  2, 5, 11 };
            size_t W2 = 0, H2 = 0, m2 = 0, suv = 0, s3d = 0;
            double f2 = 0.0;
            snprintf(path, sizeof path, "%s/t8_bigedge.tif", dir);
            CHECK(Rawtex_write_tif(path, &ct, vbig, uvq, 4, fq, 2, NULL, NULL,
                                   0.0, 1, 1.0, 1.0, 0.0, 255.0,
                                   4.0, 25.0, 6.0, NULL,
                                   &W2, &H2, &f2, &m2, &suv, &s3d) == 0
                  && f2 == 0.0 && s3d == 2,
                  "t8 big-3D-edge faces skipped (fill membranes)");
        }
    }

    /* 8b. diagnostics: class/dark split, verr + stretch flag a broken vertex */
    {
        CubeTable ct;
        long chunk = 16;
        uint8_t *buf = (uint8_t *)ARENA_ALLOC(arena, chunk * chunk * chunk);
        long z = 0, y = 0, x = 0;
        float fake[6] = { 1.0f, 1.0f, 1.0f, 14.0f, 6.0f, 14.0f };
        /* 9x9 strip: verts (z=j, y=4, x=2+i), uv (i+0.25, j+0.25) -> the
         * quarter offset keeps every pixel center strictly INSIDE one
         * triangle (centers on vertices/edges would classify as multi).
         * v == z - const exactly, u == arclen exactly => isometric, verr 0.
         * Field: dark for x<6, bright for x>=6. */
        float verts[9 * 9 * 3], uvq[9 * 9 * 2];
        int32_t fq[8 * 8 * 2 * 3];
        size_t nvq = 0, nfq = 0;
        int ii = 0, jj = 0;
        char path[512], dpath[600];
        size_t W = 0, H = 0, multi = 0;
        double fill = 0.0;
        DiagOpts dopts;
        for (z = 0; z < chunk; z++)
            for (y = 0; y < chunk; y++)
                for (x = 0; x < chunk; x++)
                    buf[(size_t)((z * chunk + y) * chunk + x)]
                        = (uint8_t)(x < 6 ? 0 : 200);
        cubetable_init(&ct, arena, dir, chunk, fake, 2, 0.0);
        ct.slot[0] = buf;
        for (jj = 0; jj < 9; jj++) {
            for (ii = 0; ii < 9; ii++) {
                verts[nvq * 3 + 0] = (float)jj;
                verts[nvq * 3 + 1] = 4.0f;
                verts[nvq * 3 + 2] = (float)(2 + ii);
                uvq[nvq * 2 + 0] = (float)ii + 0.25f;
                uvq[nvq * 2 + 1] = (float)jj + 0.25f;
                nvq++;
            }
        }
        /* break ONE interior vertex's v (mapping error, not geometry) */
        uvq[(4 * 9 + 6) * 2 + 1] += 2.0f;
        for (jj = 0; jj < 8; jj++) {
            for (ii = 0; ii < 8; ii++) {
                int32_t a = (int32_t)(jj * 9 + ii), b = a + 1;
                int32_t c = a + 9, d = c + 1;
                fq[nfq * 3 + 0] = a; fq[nfq * 3 + 1] = b; fq[nfq * 3 + 2] = c;
                nfq++;
                fq[nfq * 3 + 0] = b; fq[nfq * 3 + 1] = d; fq[nfq * 3 + 2] = c;
                nfq++;
            }
        }
        snprintf(path, sizeof path, "%s/t8b_tex.tif", dir);
        snprintf(dpath, sizeof dpath, "%s/t8b", dir);
        memset(&dopts, 0, sizeof dopts);
        dopts.prefix = dpath;
        dopts.dark_thresh = 13;
        CHECK(Rawtex_write_tif(path, &ct, verts, uvq, nvq, fq, nfq, NULL, NULL,
                               0.0, 1, 1.0, 1.0, 0.0, 255.0,
                               4.0, 25.0, 6.0, &dopts,
                               &W, &H, &fill, &multi, NULL, NULL) == 0
              && W == 9 && H == 9, "t8b diag raster runs");
        {
            uint8_t *cls = NULL, *verr = NULL, *str = NULL;
            int D = 0, IH = 0, IW = 0;
            snprintf(dpath, sizeof dpath, "%s/t8b_diagclass.tif", dir);
            CHECK(TiffIO_load(arena, dpath, &cls, &D, &IH, &IW) == 0
                  && IW == 9 && IH == 9, "t8b class map loads");
            /* x = 1.75+u_px: px u=1 -> x=2.75 dark(5); u=7 -> x=8.75 ok(1)
             * (row 2 is away from the broken vertex at (6.25, 6.25)) */
            CHECK(cls != NULL && cls[2 * 9 + 1] == 5 && cls[2 * 9 + 7] == 1,
                  "t8b dark/ok split at the field boundary");
            snprintf(dpath, sizeof dpath, "%s/t8b_diagverr.tif", dir);
            CHECK(TiffIO_load(arena, dpath, &verr, &D, &IH, &IW) == 0
                  && verr != NULL && verr[1 * 9 + 1] == 0
                  && verr[4 * 9 + 6] == 64,
                  "t8b verr: 0 on clean grid, 2 vox at broken vertex");
            snprintf(dpath, sizeof dpath, "%s/t8b_diagstretch.tif", dir);
            CHECK(TiffIO_load(arena, dpath, &str, &D, &IH, &IW) == 0
                  && str != NULL && str[1 * 9 + 1] == 128
                  && str[4 * 9 + 6] > 132,
                  "t8b stretch: isometric 128 clean, >128 at broken vertex");
        }
    }

    /* 9. degenerate inputs (rule 18) */
    {
        UvMesh m;
        size_t W = 0, H = 0;
        double fill = 0.0;
        float *n = NULL;
        size_t mu = 0;
        CHECK(parse_uv_obj("output/obj_bake_raw_selftest/nonexistent.obj",
                           &m) != 0, "t9 missing file -> error");
        CHECK(Rawtex_write_tif("x", NULL, NULL, NULL, 0, NULL, 0, NULL, NULL,
                               0.0, 1, 1.0, 1.0, 0.0, 255.0, 4.0, 25.0, 6.0,
                               NULL, &W, &H, &fill, &mu, NULL, NULL) != 0,
              "t9 empty raster -> error");
        n = vertex_normals((const float *)&fill, 0, NULL, 0);
        CHECK(n != NULL, "t9 zero-vert normals ok");
        free(n);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n",
            nfail == 0 ? "ALL PASS" : "FAILURES", nfail);
    return nfail;
}

/* ------------------------------------------------------------------ main */

static void usage(void)
{
    fprintf(stderr,
        "usage: obj_bake_raw <uv.obj> <raw_cube_dir> <out_dir> [options]\n"
        "       obj_bake_raw --selftest\n"
        "  Bake RAW CT intensity as grayscale per-vertex colors onto a mesh\n"
        "  whose OBJ carries source-space (z,y,x) verts + per-vertex vt (u,v).\n"
        "options:\n"
        "  --id S              output prefix (default: bake)\n"
        "  --chunk N           cube edge, vox (default 128)\n"
        "  --normal-range F    max-sample +/-F vox along the vertex normal\n"
        "                      (default 2.0; 0 = single trilinear sample)\n"
        "  --normal-samples N  steps along the normal (default 5)\n"
        "  --pct-lo F          contrast stretch low percentile (default 1)\n"
        "  --pct-hi F          contrast stretch high percentile (default 99)\n"
        "  --no-stretch        disable contrast stretch (window 0..255)\n"
        "  --require-contrast  fail if the sampled percentile window is constant\n"
        "  --no-raster         skip the rasterized texture TIFF (on by\n"
        "                      default: every pixel gets a fresh volume\n"
        "                      sample at its barycentric surface point;\n"
        "                      per-vertex scattering is never used)\n"
        "  --raster-du F       raster pixel size in u, vox (default 1.0)\n"
        "  --raster-dv F       raster pixel size in v, vox (default 1.0)\n"
        "  --raster-stretch-ratio F / --raster-stretch-floor F\n"
        "                      skip faces with a UV edge > max(ratio*3d_len,\n"
        "                      floor) -- collapsed-core streak paint (4.0/25)\n"
        "  --raster-max-edge F skip faces with a 3D edge > F vox -- hole-fill\n"
        "                      membranes off the true sheet (default 0=off)\n"
        "  --no-diag           skip the diagnostic rasters (written by default\n"
        "                      with the texture: <id>_diagclass.tif 0=bg 1=ok\n"
        "                      2=skip-uv 3=skip-3d 4=multi 5=dark;\n"
        "                      _diagstretch/_diagsquash = 128+42.5*log2(sigma);\n"
        "                      _diagverr = |vt.v-(z-zmin)|*32) + hot-tile table\n"
        "  --diag-dark N       dark-pixel threshold, post-stretch u8 (def 13)\n"
        "  --no-3d / --no-flat skip that output\n");
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *raw_dir = NULL, *out_dir = NULL;
    const char *id = "bake";
    long chunk = 128;
    double normal_range = 2.0;
    int normal_samples = 5;
    double pct_lo = 1.0, pct_hi = 99.0;
    int do_stretch = 1, do_raster = 1, do_3d = 1, do_flat = 1;
    int require_contrast = 0, output_failures = 0;
    double raster_du = 1.0, raster_dv = 1.0;
    double raster_stretch_ratio = 4.0, raster_stretch_floor = 25.0;
    /* Real input-mesh edge scale belongs to the upstream remesher. A fixed
     * 6-voxel cap hid most valid faces in coarse welded scroll meshes; callers
     * rasterizing known synthetic membrane fill can opt back in explicitly. */
    double raster_max_edge = 0.0;
    int do_diag = 1, diag_dark = 13;
    const char *dump_smear = NULL;
    int i = 0;
    UvMesh m;
    CubeTable ct;
    Arena_T arena = NULL;
    float *normals = NULL;
    double *val = NULL;
    uint8_t *has = NULL;
    size_t nsamp = 0, k = 0;
    double lo = 0.0, hi = 255.0;
    double sample_min = 1.0e300, sample_max = -1.0e300;
    double t0 = 0.0, t_parse = 0.0, t_sample = 0.0, t_write = 0.0;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return selftest() == 0 ? 0 : 1;
    if (argc < 4) { usage(); return 1; }
    in_path = argv[1];
    raw_dir = argv[2];
    out_dir = argv[3];
    for (i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--id") == 0 && i + 1 < argc) id = argv[++i];
        else if (strcmp(argv[i], "--chunk") == 0 && i + 1 < argc)
            chunk = strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--normal-range") == 0 && i + 1 < argc)
            normal_range = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--normal-samples") == 0 && i + 1 < argc)
            normal_samples = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--pct-lo") == 0 && i + 1 < argc)
            pct_lo = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--pct-hi") == 0 && i + 1 < argc)
            pct_hi = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--no-stretch") == 0) do_stretch = 0;
        else if (strcmp(argv[i], "--require-contrast") == 0) require_contrast = 1;
        else if (strcmp(argv[i], "--raster") == 0) do_raster = 1;  /* default */
        else if (strcmp(argv[i], "--no-raster") == 0) do_raster = 0;
        else if (strcmp(argv[i], "--raster-du") == 0 && i + 1 < argc)
            raster_du = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-dv") == 0 && i + 1 < argc)
            raster_dv = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-stretch-ratio") == 0 && i + 1 < argc)
            raster_stretch_ratio = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-stretch-floor") == 0 && i + 1 < argc)
            raster_stretch_floor = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-max-edge") == 0 && i + 1 < argc)
            raster_max_edge = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--no-diag") == 0) do_diag = 0;
        else if (strcmp(argv[i], "--diag-dark") == 0 && i + 1 < argc)
            diag_dark = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--dump-smear") == 0 && i + 1 < argc)
            dump_smear = argv[++i];
        else if (strcmp(argv[i], "--no-3d") == 0) do_3d = 0;
        else if (strcmp(argv[i], "--no-flat") == 0) do_flat = 0;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); usage(); return 1; }
    }
    if (chunk <= 0 || normal_samples < 1) { usage(); return 1; }

    fprintf(stderr, "[obj_bake_raw] in=%s raw=%s out=%s id=%s\n",
            in_path, raw_dir, out_dir, id);
    t0 = ves_clock_sec();
    if (parse_uv_obj(in_path, &m) != 0) return 1;
    t_parse = ves_clock_sec() - t0;
    fprintf(stderr, "  parsed: %zu verts, %zu vt, %zu faces (%.1fs)\n",
            m.nv, m.nvt, m.nf, t_parse);
    if (m.nvt != m.nv) {
        fprintf(stderr, "  WARN: nvt (%zu) != nv (%zu) -- flat/raster outputs "
                "disabled\n", m.nvt, m.nv);
        do_flat = 0;
        do_raster = 0;
    }

    arena = Arena_new();
    if (cubetable_init(&ct, arena, raw_dir, chunk, m.verts, m.nv,
                       normal_range + 2.0) != 0) {
        fprintf(stderr, "ERROR: cube table init failed\n");
        return 1;
    }
    fprintf(stderr, "  cube table: %ldx%ldx%ld cubes (chunk %ld) over bbox\n",
            ct.nz, ct.ny, ct.nx, chunk);

    if (m.nf > 0) normals = vertex_normals(m.verts, m.nv, m.faces, m.nf);
    val = (double *)xmalloc(m.nv * sizeof(double));
    has = (uint8_t *)xcalloc(m.nv, sizeof(uint8_t));

    t0 = ves_clock_sec();
    for (k = 0; k < m.nv; k++) {
        double s = sample_vertex(&ct, &m.verts[k * 3],
                                 normals != NULL ? &normals[k * 3] : NULL,
                                 normal_range, normal_samples);
        val[k] = (s >= 0.0) ? s : 0.0;
        if (s >= 0.0) {
            has[k] = 1;
            nsamp++;
            if (s < sample_min) sample_min = s;
            if (s > sample_max) sample_max = s;
        }
    }
    t_sample = ves_clock_sec() - t0;
    fprintf(stderr, "  sampled: %zu/%zu verts (%.2f%%), %d cubes loaded, "
            "%d missing (%.1fs)\n", nsamp, m.nv,
            m.nv > 0 ? 100.0 * (double)nsamp / (double)m.nv : 0.0,
            ct.n_loaded, ct.n_missing, t_sample);

    if (do_stretch) stretch_window(val, has, m.nv, pct_lo, pct_hi, &lo, &hi);
    fprintf(stderr, "  contrast window: [%.0f, %.0f] (stretch %s)\n",
            lo, hi, do_stretch ? "on" : "off");
    if (require_contrast && (nsamp == 0 || !(sample_max > sample_min))) {
        fprintf(stderr, "ERROR: sampled RAW data has no intensity contrast "
                "(range [%.0f, %.0f])\n", sample_min, sample_max);
        free(normals); free(val); free(has);
        uvmesh_free(&m);
        Arena_dispose(&arena);
        return 1;
    }
    {   /* value distribution over sampled verts */
        double q[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
        double pcts[5] = { 1.0, 25.0, 50.0, 75.0, 99.0 };
        int qi = 0;
        for (qi = 0; qi < 5; qi++) {
            double wlo = 0.0, whi = 0.0;
            stretch_window(val, has, m.nv, pcts[qi], 100.0, &wlo, &whi);
            q[qi] = wlo;
        }
        fprintf(stderr, "  intensity percentiles p1/p25/p50/p75/p99: "
                "%.0f / %.0f / %.0f / %.0f / %.0f\n",
                q[0], q[1], q[2], q[3], q[4]);
    }

    t0 = ves_clock_sec();
    if (do_3d) {
        char path[2600];
        snprintf(path, sizeof path, "%s/%s_raw3d.obj", out_dir, id);
        ves_ensure_parent_dir(path);
        if (write_colored_obj(path, m.verts, m.nv, m.faces, m.nf, m.uv,
                              val, has, lo, hi, 0) != 0) {
            fprintf(stderr, "  WARN: cannot write %s\n", path);
            output_failures++;
        } else
            fprintf(stderr, "  wrote %s\n", path);
    }
    if (do_flat) {
        char path[2600];
        snprintf(path, sizeof path, "%s/%s_rawflat.obj", out_dir, id);
        ves_ensure_parent_dir(path);
        if (write_colored_obj(path, m.verts, m.nv, m.faces, m.nf, m.uv,
                              val, has, lo, hi, 1) != 0) {
            fprintf(stderr, "  WARN: cannot write %s\n", path);
            output_failures++;
        } else
            fprintf(stderr, "  wrote %s\n", path);
    }
    if (do_raster) {
        char path[2600], dprefix[2600];
        size_t W = 0, H = 0, multi = 0, skuv = 0, sk3d = 0;
        double fill = 0.0, tr = ves_clock_sec();
        DiagOpts dopts;
        memset(&dopts, 0, sizeof dopts);
        snprintf(path, sizeof path, "%s/%s_rawtex.tif", out_dir, id);
        snprintf(dprefix, sizeof dprefix, "%s/%s", out_dir, id);
        dopts.prefix = dprefix;
        dopts.dark_thresh = diag_dark;
        dopts.smear_obj = dump_smear;
        ves_ensure_parent_dir(path);
        if (Rawtex_write_tif(path, &ct, m.verts, m.uv, m.nv, m.faces, m.nf,
                             normals, NULL /*face_skip*/, normal_range, normal_samples,
                             raster_du, raster_dv, lo, hi,
                             raster_stretch_ratio, raster_stretch_floor,
                             raster_max_edge,
                             do_diag ? &dopts : NULL,
                             &W, &H, &fill, &multi, &skuv, &sk3d) != 0) {
            fprintf(stderr, "  WARN: cannot write %s\n", path);
            output_failures++;
        } else
            fprintf(stderr, "  wrote %s (%zux%zu px, %.1f%% filled, "
                    "%zu multi-cover px, skipped %zu uv-stretched + %zu "
                    "big-3d faces, %.1fs)\n",
                    path, W, H, 100.0 * fill, multi, skuv, sk3d,
                    ves_clock_sec() - tr);
    }
    t_write = ves_clock_sec() - t0;

    fprintf(stderr, "[obj_bake_raw] SUMMARY id=%s verts=%zu sampled=%.2f%% "
            "cubes=%d/%d window=[%.0f,%.0f] parse=%.1fs sample=%.1fs "
            "write=%.1fs\n", id, m.nv,
            m.nv > 0 ? 100.0 * (double)nsamp / (double)m.nv : 0.0,
            ct.n_loaded, ct.n_loaded + ct.n_missing, lo, hi,
            t_parse, t_sample, t_write);

    free(normals); free(val); free(has);
    uvmesh_free(&m);
    Arena_dispose(&arena);
    return output_failures == 0 ? 0 : 1;
}
