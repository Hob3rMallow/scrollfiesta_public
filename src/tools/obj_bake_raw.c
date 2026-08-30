/* ============================================================================
 * obj_bake_raw -- bake RAW CT intensity onto a UV-carrying OBJ.
 *
 * Validation step for the scroll_ribbon unwrap: "read the original volumetric
 * information back from the scroll and write it out to an OBJ -- if this
 * works, we should see nice, clean papyrus."
 *
 * Input: an OBJ whose vertices are SOURCE-SPACE voxel coordinates (z,y,x)
 * with per-vertex "vt u v" texture coordinates (as written by scroll_ribbon:
 * <id>_uv.obj or <id>_ribbon.obj locator whose authoritative .vmesh companion
 * is read, plus the grid's cubes_RAW/ directory of
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
#include "../common/mesh_normals.h"
#include "../common/mesh_bin.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"
#include "../common/raw_sample.h"
#include "../flatten/rawtex_bake.h"   /* Rawtex_write_tif + DiagOpts (extracted) */
#include "../flatten/marble_diag.h"

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
 * (only the vertex index is kept; scroll_ribbon writes a/a so vt index ==
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

/* The OBJ argument remains a convenient name for the human-readable sibling,
 * but production baking consumes only its VESMESH1 companion. */
static int load_uv_mesh_binary(const char *locator, UvMesh *m,
                               char *resolved, size_t resolved_cap)
{
    MeshBinData mesh;
    memset(m, 0, sizeof *m);
    if (MeshBin_companion_path(locator, resolved, resolved_cap) != 0 ||
        MeshBin_read_malloc(resolved, &mesh) != 0) {
        fprintf(stderr,
                "ERROR: cannot load authoritative mesh container for %s "
                "(OBJ fallback is disabled)\n", locator);
        return -1;
    }
    if (mesh.uv == NULL) {
        fprintf(stderr, "ERROR: %s has no per-vertex UV field\n", resolved);
        MeshBin_dispose(&mesh);
        return -1;
    }
    m->verts = mesh.verts;
    m->uv = mesh.uv;
    m->faces = mesh.faces;
    m->nv = m->nvt = mesh.nv;
    m->nf = mesh.nf;
    return 0;
}

/* cube table + RAW sampling now live in src/common/raw_sample.{c,h}
 * (CubeTable, cubetable_init, cube_fetch, sample_trilinear, sample_vertex,
 * sampling primitives) -- shared with the overlap-repair module. */

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

/* Colored, parameterized OBJ: "v x y z r g b" plus one vt per vertex.
 * flat != 0 places vertices at (u, v, 0) from uv instead of their 3D
 * position. Unsampled verts get 0. */
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
    for (i = 0; i < nv; i++)
        fprintf(f, "vt %.9g %.9g\n", (double)uv[i * 2 + 0],
                (double)uv[i * 2 + 1]);
    for (t = 0; t < nf; t++)
        fprintf(f, "f %d/%d %d/%d %d/%d\n",
                faces[t * 3 + 0] + 1, faces[t * 3 + 0] + 1,
                faces[t * 3 + 1] + 1, faces[t * 3 + 1] + 1,
                faces[t * 3 + 2] + 1, faces[t * 3 + 2] + 1);
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
        CHECK(cubetable_expected_chunks(&ct) == 2 && cubetable_is_complete(&ct),
              "t5 complete two-chunk bbox is production-safe");
        v = cube_fetch(&ct, 2, 9, 3);   /* y=9 -> cube y0=8: no file */
        (void)v;
        CHECK(ct.n_missing >= 0 && cube_fetch(&ct, 2, 9, 3) == -1,
              "t5 missing neighbor -> -1");
        {
            CubeTable partial;
            float wider[6] = { 1.0f, 1.0f, 1.0f, 6.0f, 14.0f, 14.0f };
            CHECK(cubetable_init(&partial, arena, cdir, chunk, wider, 2, 0.0) == 0,
                  "t5 strict-coverage table init");
            cubetable_prewarm_all(&partial);
            CHECK(cubetable_expected_chunks(&partial) == 4 &&
                  partial.n_loaded == 2 && partial.n_missing == 2 &&
                  !cubetable_is_complete(&partial),
                  "t5 partial bbox is rejected by production coverage");
        }
    }

    /* 6. Zarr coverage clips the padded geometry bbox to the physical array.
     * Out-of-volume halo is absence of evidence, not a missing in-volume
     * chunk; sampling it must still return no datum. */
    {
        CubeTable ct, outside;
        const long chunk = 8;
        const size_t bytes = (size_t)chunk * (size_t)chunk * (size_t)chunk;
        char zdir[512], path[700];
        uint8_t *buf = (uint8_t *)xmalloc(bytes);
        float crossing[6] = { -4.0f, -4.0f, -4.0f, 12.0f, 12.0f, 14.0f };
        float beyond[6] = { 20.0f, 20.0f, 20.0f, 24.0f, 24.0f, 24.0f };
        FILE *f = NULL;
        snprintf(zdir, sizeof zdir, "%s/zarr_boundary", dir);
        snprintf(path, sizeof path, "%s/0/.zarray", zdir);
        ves_ensure_parent_dir(path);
        f = fopen(path, "wb");
        if (f != NULL) {
            fputs("{\"zarr_format\":2,\"shape\":[8,8,10],"
                  "\"chunks\":[8,8,8],\"dtype\":\"|u1\","
                  "\"compressor\":null,\"fill_value\":0,"
                  "\"order\":\"C\",\"filters\":null}\n", f);
            fclose(f);
        }
        memset(buf, 31, bytes);
        snprintf(path, sizeof path, "%s/0/0/0/0", zdir);
        ves_ensure_parent_dir(path);
        f = fopen(path, "wb");
        if (f != NULL) { fwrite(buf, 1, bytes, f); fclose(f); }
        memset(buf, 197, bytes);
        snprintf(path, sizeof path, "%s/0/0/0/1", zdir);
        ves_ensure_parent_dir(path);
        f = fopen(path, "wb");
        if (f != NULL) { fwrite(buf, 1, bytes, f); fclose(f); }
        free(buf);

        CHECK(cubetable_init(&ct, arena, zdir, chunk, crossing, 2, 0.0) == 0,
              "t6 boundary-crossing Zarr table init");
        cubetable_prewarm_all(&ct);
        CHECK(cubetable_expected_chunks(&ct) == 2 && ct.n_loaded == 2 &&
              ct.n_missing == 0 && ct.n_outside == 25 &&
              cubetable_is_complete(&ct),
              "t6 only in-volume Zarr chunks participate in strict coverage");
        CHECK(cube_fetch(&ct, 2, 2, 9) == 197,
              "t6 valid voxel in partial boundary chunk loads");
        CHECK(cube_fetch(&ct, 2, 2, 10) == -1 &&
              cube_fetch(&ct, -1, 2, 2) == -1,
              "t6 physical Zarr boundary returns no datum");

        CHECK(cubetable_init(&outside, arena, zdir, chunk, beyond, 2, 0.0) == 0,
              "t6 fully out-of-volume Zarr table init");
        cubetable_prewarm_all(&outside);
        CHECK(cubetable_expected_chunks(&outside) == 0 &&
              outside.n_loaded == 0 && outside.n_missing == 0 &&
              outside.n_outside == 8 && cubetable_is_complete(&outside),
              "t6 fully out-of-volume bbox is complete but has no evidence");
        CHECK(cube_fetch(&outside, 22, 22, 22) == -1,
              "t6 fully out-of-volume fetch returns no datum");
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
                               4.0, 25.0, 6.0, 0, NULL,
                               &W, &H, &fill, &multi, NULL, NULL) == 0,
              "t8 raster runs");
        CHECK(W == 4 && H == 4 && fabs(fill - 1.0) < 1e-9,
              "t8 raster 4x4 (half-texel cells), fully covered");
        {   /* centers at +0.5: field(2, 1.5+pv, 1.5+pu) = 9.5 + 2pv + 3pu;
             * the writer rounds g*255+0.5 -> 10 + 2pv + 3pu exactly */
            uint8_t *img = NULL;
            int D = 0, IH = 0, IW = 0;
            int okpix = 1, py = 0, pxx = 0;
            CHECK(TiffIO_load(arena, path, &img, &D, &IH, &IW) == 0
                  && D == 1 && IH == 4 && IW == 4, "t8 raster loads back");
            for (py = 0; py < 4 && img != NULL; py++)
                for (pxx = 0; pxx < 4; pxx++)
                    if (img[py * 4 + pxx]
                        != (uint8_t)(10 + 2 * py + 3 * pxx)) okpix = 0;
            CHECK(okpix, "t8 every pixel exact on linear field");
            CHECK(multi == 0,
                  "t8 shared diagonals/vertices have one raster owner");
        }
        {   /* The raster must be sized from the UV BOUNDING BOX, not [0,max].
             * A winding frame lifted from registration starts wherever the
             * spiral starts -- on the 4x21x21 u is negative almost everywhere
             * and v is world z (~4352..4864).  Anchoring at (0,0) dropped 99.98%
             * of that mesh and asked for a 341208x4865 raster. */
            float uvoff[8];
            size_t W3 = 0, H3 = 0, multi3 = 0;
            double fill3 = 0.0;
            char path3[512];
            uint8_t *img3 = NULL;
            int D3 = 0, IH3 = 0, IW3 = 0, okpix3 = 1, py3 = 0, px3 = 0, m3 = 0;
            for (m3 = 0; m3 < 4; m3++) {
                uvoff[m3 * 2 + 0] = uvq[m3 * 2 + 0] - 1000.0f;
                uvoff[m3 * 2 + 1] = uvq[m3 * 2 + 1] + 4352.0f;
            }
            snprintf(path3, sizeof path3, "%s/t8c_offset.tif", dir);
            CHECK(Rawtex_write_tif(path3, &ct, verts, uvoff, 4, fq, 2, NULL,
                                   NULL, 0.0, 1, 1.0, 1.0, 0.0, 255.0,
                                   4.0, 25.0, 6.0, 0, NULL,
                                   &W3, &H3, &fill3, &multi3, NULL, NULL) == 0,
                  "t8c offset-origin raster runs");
            CHECK(W3 == 4 && H3 == 4 && fabs(fill3 - 1.0) < 1e-9,
                  "t8c offset UV still gives a 4x4 fully covered raster");
            CHECK(TiffIO_load(arena, path3, &img3, &D3, &IH3, &IW3) == 0
                  && D3 == 1 && IH3 == 4 && IW3 == 4,
                  "t8c offset raster loads back");
            for (py3 = 0; py3 < 4 && img3 != NULL; py3++)
                for (px3 = 0; px3 < 4; px3++)
                    if (img3[py3 * 4 + px3]
                        != (uint8_t)(10 + 2 * py3 + 3 * px3)) okpix3 = 0;
            CHECK(okpix3,
                  "t8c offset raster is pixel-identical to the origin raster");
        }
        {   /* degenerate UV face (zero area) is skipped cleanly */
            int32_t fdeg[3] = { 0, 0, 1 };
            size_t W2 = 0, H2 = 0, m2 = 0;
            double f2 = 0.0;
            snprintf(path, sizeof path, "%s/t8_degen.tif", dir);
            CHECK(Rawtex_write_tif(path, &ct, verts, uvq, 4, fdeg, 1, NULL, NULL,
                                   0.0, 1, 1.0, 1.0, 0.0, 255.0,
                                   4.0, 25.0, 6.0, 0, NULL,
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
                                   4.0, 25.0, 6.0, 0, NULL,
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
                                   4.0, 25.0, 6.0, 0, NULL,
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
                               4.0, 25.0, 6.0, 0, &dopts,
                               &W, &H, &fill, &multi, NULL, NULL) == 0
              && W == 8 && H == 8, "t8b diag raster runs");
        {
            uint8_t *cls = NULL, *verr = NULL, *str = NULL;
            int D = 0, IH = 0, IW = 0;
            int qx = 0, qy = 0;
            uint8_t verr_max = 0, str_max = 0;
            snprintf(dpath, sizeof dpath, "%s/t8b_diagclass.tif", dir);
            CHECK(TiffIO_load(arena, dpath, &cls, &D, &IH, &IW) == 0
                  && IW == 8 && IH == 8, "t8b class map loads");
            /* half-texel center px -> x = 2.5+px: px=1 -> x=3.5 dark(5);
             * px=6 -> x=8.5 ok(1) (row 2 is away from the broken vertex) */
            CHECK(cls != NULL && cls[2 * 8 + 1] == 5 && cls[2 * 8 + 6] == 1,
                  "t8b dark/ok split at the field boundary");
            snprintf(dpath, sizeof dpath, "%s/t8b_diagverr.tif", dir);
            CHECK(TiffIO_load(arena, dpath, &verr, &D, &IH, &IW) == 0
                  && verr != NULL && verr[1 * 8 + 1] == 0,
                  "t8b verr: 0 on the clean grid");
            snprintf(dpath, sizeof dpath, "%s/t8b_diagstretch.tif", dir);
            CHECK(TiffIO_load(arena, dpath, &str, &D, &IH, &IW) == 0
                  && str != NULL && str[1 * 8 + 1] == 128,
                  "t8b stretch: isometric 128 on the clean grid");
            for (qy = 0; qy < 8; qy++)
                for (qx = 0; qx < 8; qx++) {
                    if (verr != NULL && verr[qy * 8 + qx] > verr_max)
                        verr_max = verr[qy * 8 + qx];
                    if (str != NULL && str[qy * 8 + qx] > str_max)
                        str_max = str[qy * 8 + qx];
                }
            /* interpolated centers no longer sit ON the broken vertex, so the
             * 2-vox verr peak (64) reads back attenuated but must still be
             * loud, and the stretch map must flag the distortion */
            CHECK(verr_max >= 32 && str_max > 132,
                  "t8b broken vertex flagged by verr + stretch maxima");
        }
        {   /* RAW/mesh alignment audit writes evidence maps + green overlay. */
            MarbleDiagOpts mopts;
            MarbleDiagStats mstats;
            FILE *probe_file = NULL;
            MarbleDiag_defaults(&mopts);
            mopts.stride_pixels = 2;
            mopts.depth_range = 1.0;
            mopts.depth_step = 1.0;
            mopts.tensor_radius = 1.0;
            mopts.score_threshold = 1.0;
            mopts.minimum_region_cells = 1;
            snprintf(dpath, sizeof dpath, "%s/t8b_marble", dir);
            CHECK(MarbleDiag_write(dpath, path, &ct, verts, uvq, nvq,
                                   fq, nfq, NULL, 1.0, 1.0, 0,
                                   &mopts, &mstats) == 0 &&
                  mstats.audit_width == 4 && mstats.audit_height == 4 &&
                  mstats.sampled_cells > 0,
                  "t8b RAW/mesh marble audit runs on sparse UV grid");
            snprintf(dpath, sizeof dpath,
                     "%s/t8b_marble_marble_overlay.png", dir);
            probe_file = fopen(dpath, "rb");
            CHECK(probe_file != NULL,
                  "t8b marble audit writes the green-overlay PNG");
            if (probe_file != NULL) fclose(probe_file);
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
                               0, NULL, &W, &H, &fill, &mu, NULL, NULL) != 0,
              "t9 empty raster -> error");
        n = MeshNormals_compute((const float *)&fill, 0, NULL, 0);
        CHECK(n != NULL, "t9 zero-vert normals ok");
        free(n);
    }

    /* 10. raster preflight: cheap sizing BEFORE sampling, per-axis advice.
     * The 21x21x21 bake burned 446s of vertex sampling and THEN rejected the
     * raster; the old retry advice also scaled dv by sqrt(area), silently
     * downsampling a full-height v axis because u had exploded. */
    {
        RawtexPlan rp, rp2;
        float uvbig[64 * 2];
        float uvq2[8] = { 0, 0,  4, 0,  0, 4,  4, 4 };
        int t = 0;
        for (t = 0; t < 64; t++) {          /* u spans 3e6, v spans 500 */
            uvbig[t * 2 + 0] = (float)t / 63.0f * 3.0e6f;
            uvbig[t * 2 + 1] = (float)t / 63.0f * 500.0f;
        }
        CHECK(Rawtex_plan(uvbig, 64, 1.0, 1.0, 0, &rp) == 0 && rp.ok == 0,
              "t10 huge-u raster rejected by preflight");
        CHECK(rp.need_du > 1.0 && fabs(rp.need_dv - 1.0) < 1e-12,
              "t10 advice shrinks u only; v keeps full resolution");
        CHECK(Rawtex_plan(uvbig, 64, rp.need_du, rp.need_dv, 0, &rp2) == 0
              && rp2.ok == 1,
              "t10 suggested steps fit the caps");
        CHECK(Rawtex_plan(uvq2, 4, 1.0, 1.0, 0, &rp) == 0 && rp.ok == 1
              && rp.W == 4 && rp.H == 4,
              "t10 plan dims match the rasterizer (4x4 quad)");
        CHECK(Rawtex_plan(uvq2, 4, 1.0, 1.0, 8, &rp) == 0 && rp.ok == 0
              && Rawtex_plan(uvq2, 4, rp.need_du, rp.need_dv, 8, &rp2) == 0
              && rp2.ok == 1 && rp2.W * rp2.H <= 8,
              "t10 custom --raster-max-px cap honored with fitting advice");
        CHECK(Rawtex_plan(NULL, 0, 1.0, 1.0, 0, &rp) != 0,
              "t10 empty input -> error");
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
        "usage: obj_bake_raw <uv.obj> <raw_cube_dir_or_zarr> <out_dir> [options]\n"
        "       (production reads only the authoritative companion .vmesh)\n"
        "       obj_bake_raw --selftest\n"
        "       obj_bake_raw --migrate-obj <legacy.obj> [output.vmesh]\n"
        "  Bake RAW CT intensity as grayscale per-vertex colors onto a mesh\n"
        "  whose OBJ carries source-space (z,y,x) verts + per-vertex vt (u,v).\n"
        "options:\n"
        "  --id S              output prefix (default: bake)\n"
        "  --chunk N           cube edge, vox (default 128)\n"
        "  --require-complete-raw  fail if any required RAW chunk is absent (default)\n"
        "  --allow-incomplete-raw  diagnostics only; never publish this output\n"
        "  --normal-range F    max-sample +/-F vox along the vertex normal\n"
        "                      (default 2.0; 0 = single trilinear sample)\n"
        "  --normal-samples N  steps along the normal (default 5)\n"
        "  --pct-lo F          contrast stretch low percentile (default 1)\n"
        "  --pct-hi F          contrast stretch high percentile (default 99)\n"
        "  --window LO HI      use an explicit shared contrast window\n"
        "  --no-stretch        disable contrast stretch (window 0..255)\n"
        "  --no-raster         skip the rasterized texture TIFF (on by\n"
        "                      default: every pixel gets a fresh volume\n"
        "                      sample at its barycentric surface point;\n"
        "                      per-vertex scattering is never used)\n"
        "  --raster-du F       raster pixel size in u, vox (default 1.0)\n"
        "  --raster-dv F       raster pixel size in v, vox (default 1.0)\n"
        "  --raster-auto       if (du,dv) would exceed the raster caps, adopt\n"
        "                      the smallest per-axis steps that fit instead of\n"
        "                      erroring (the preflight runs BEFORE sampling)\n"
        "  --raster-max-px N   total-pixel cap (default 268435456 = 1<<28);\n"
        "                      raise on big-RAM boxes to avoid downsampling\n"
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
        "  --marble-diag       audit RAW/mesh UV alignment and emit green overlay\n"
        "  --marble-stride N   audit-cell edge in baked pixels (default 16)\n"
        "  --marble-depth F    normal search reach +/-F vox (default 4)\n"
        "  --marble-depth-step F  normal search spacing (default 1)\n"
        "  --marble-radius F   RAW tangent-patch radius (default 2)\n"
        "  --marble-threshold F  green-mask score threshold (default .40)\n"
        "  --marble-min-cells N  remove smaller audit islands (default 12)\n"
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
    double fixed_lo = 0.0, fixed_hi = 255.0;
    int fixed_window = 0;
    int do_stretch = 1, do_raster = 1, do_3d = 1, do_flat = 1;
    double raster_du = 1.0, raster_dv = 1.0;
    int raster_auto = 0;
    size_t raster_max_px = 0;   /* 0 = Rawtex default (1<<28) */
    double raster_stretch_ratio = 4.0, raster_stretch_floor = 25.0;
    /* Real input-mesh edge scale belongs to the upstream remesher. A fixed
     * 6-voxel cap hid most valid faces in coarse welded scroll meshes; callers
     * rasterizing known synthetic membrane fill can opt back in explicitly. */
    double raster_max_edge = 0.0;
    int do_diag = 1, diag_dark = 13;
    int do_marble = 0;
    MarbleDiagOpts marble_opts;
    int require_complete_raw = 1;
    int need_vertex_samples = 1;
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
    double t0 = 0.0, t_parse = 0.0, t_sample = 0.0, t_write = 0.0;

    MarbleDiag_defaults(&marble_opts);

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return selftest() == 0 && MeshBin_selftest() == 0 ? 0 : 1;
    /* Explicit one-time migration valve for artifacts predating VESMESH1.
     * No normal pipeline path calls the text parser. */
    if (argc >= 3 && strcmp(argv[1], "--migrate-obj") == 0) {
        char companion[2048];
        const char *destination = argc >= 4 ? argv[3] : companion;
        if (argc > 4 ||
            (argc < 4 && MeshBin_companion_path(
                argv[2], companion, sizeof companion) != 0)) {
            usage(); return 1;
        }
        if (parse_uv_obj(argv[2], &m) != 0) return 1;
        if (MeshBin_write(destination, m.verts, m.nv, m.faces, m.nf,
                          m.nvt == m.nv ? m.uv : NULL) != 0) {
            fprintf(stderr, "ERROR: cannot migrate %s to %s\n",
                    argv[2], destination);
            uvmesh_free(&m);
            return 1;
        }
        fprintf(stderr, "migrated %s -> %s (%zu vertices, %zu faces%s)\n",
                argv[2], destination, m.nv, m.nf,
                m.nvt == m.nv ? ", UV" : "");
        uvmesh_free(&m);
        return 0;
    }
    if (argc < 4) { usage(); return 1; }
    in_path = argv[1];
    raw_dir = argv[2];
    out_dir = argv[3];
    for (i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--id") == 0 && i + 1 < argc) id = argv[++i];
        else if (strcmp(argv[i], "--chunk") == 0 && i + 1 < argc)
            chunk = strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--require-complete-raw") == 0)
            require_complete_raw = 1;
        else if (strcmp(argv[i], "--allow-incomplete-raw") == 0)
            require_complete_raw = 0;
        else if (strcmp(argv[i], "--normal-range") == 0 && i + 1 < argc)
            normal_range = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--normal-samples") == 0 && i + 1 < argc)
            normal_samples = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--pct-lo") == 0 && i + 1 < argc)
            pct_lo = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--pct-hi") == 0 && i + 1 < argc)
            pct_hi = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--window") == 0 && i + 2 < argc) {
            fixed_lo = strtod(argv[++i], NULL);
            fixed_hi = strtod(argv[++i], NULL);
            fixed_window = 1;
        }
        else if (strcmp(argv[i], "--no-stretch") == 0) do_stretch = 0;
        else if (strcmp(argv[i], "--raster") == 0) do_raster = 1;  /* default */
        else if (strcmp(argv[i], "--no-raster") == 0) do_raster = 0;
        else if (strcmp(argv[i], "--raster-du") == 0 && i + 1 < argc)
            raster_du = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-dv") == 0 && i + 1 < argc)
            raster_dv = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-auto") == 0) raster_auto = 1;
        else if (strcmp(argv[i], "--raster-max-px") == 0 && i + 1 < argc)
            raster_max_px = (size_t)strtoull(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--raster-stretch-ratio") == 0 && i + 1 < argc)
            raster_stretch_ratio = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-stretch-floor") == 0 && i + 1 < argc)
            raster_stretch_floor = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--raster-max-edge") == 0 && i + 1 < argc)
            raster_max_edge = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--no-diag") == 0) do_diag = 0;
        else if (strcmp(argv[i], "--diag-dark") == 0 && i + 1 < argc)
            diag_dark = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--marble-diag") == 0) do_marble = 1;
        else if (strcmp(argv[i], "--marble-stride") == 0 && i + 1 < argc)
            marble_opts.stride_pixels = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--marble-depth") == 0 && i + 1 < argc)
            marble_opts.depth_range = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--marble-depth-step") == 0 && i + 1 < argc)
            marble_opts.depth_step = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--marble-radius") == 0 && i + 1 < argc)
            marble_opts.tensor_radius = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--marble-threshold") == 0 && i + 1 < argc)
            marble_opts.score_threshold = strtod(argv[++i], NULL);
        else if (strcmp(argv[i], "--marble-min-cells") == 0 && i + 1 < argc)
            marble_opts.minimum_region_cells = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--dump-smear") == 0 && i + 1 < argc)
            dump_smear = argv[++i];
        else if (strcmp(argv[i], "--no-3d") == 0) do_3d = 0;
        else if (strcmp(argv[i], "--no-flat") == 0) do_flat = 0;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); usage(); return 1; }
    }
    if (chunk <= 0 || normal_samples < 1 ||
        (do_marble && (marble_opts.stride_pixels < 2 ||
                       marble_opts.depth_range <= 0.0 ||
                       marble_opts.depth_step <= 0.0 ||
                       marble_opts.tensor_radius <= 0.0 ||
                       marble_opts.score_threshold < 0.0 ||
                       marble_opts.score_threshold > 1.0 ||
                       marble_opts.minimum_region_cells < 1)) ||
        (fixed_window && (!isfinite(fixed_lo) || !isfinite(fixed_hi) ||
                          fixed_hi <= fixed_lo))) {
        usage(); return 1;
    }

    fprintf(stderr, "[obj_bake_raw] in=%s raw=%s out=%s id=%s\n",
            in_path, raw_dir, out_dir, id);
    if (!require_complete_raw)
        fprintf(stderr, "WARNING: incomplete RAW is allowed: this bake is diagnostic "
                        "and must not be published as a production sheet\n");
    t0 = ves_clock_sec();
    {
        char binary_path[2048];
        if (load_uv_mesh_binary(in_path, &m, binary_path,
                                sizeof binary_path) != 0)
            return 1;
        fprintf(stderr, "  authoritative input: %s\n", binary_path);
    }
    t_parse = ves_clock_sec() - t0;
    fprintf(stderr, "  loaded: %zu verts, %zu vt, %zu faces (%.1fs)\n",
            m.nv, m.nvt, m.nf, t_parse);
    {
        double lo3[3] = { HUGE_VAL, HUGE_VAL, HUGE_VAL };
        double hi3[3] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL };
        for (size_t v = 0; v < m.nv; v++) {
            for (int d = 0; d < 3; d++) {
                double x = (double)m.verts[v*3 + (size_t)d];
                if (!isfinite(x)) {
                    fprintf(stderr,
                            "BUG: non-finite geometry at input vertex %zu "
                            "axis %d\n", v, d);
                    uvmesh_free(&m);
                    return 1;
                }
                if (x < lo3[d]) lo3[d] = x;
                if (x > hi3[d]) hi3[d] = x;
            }
        }
        fprintf(stderr,
                "  geometry bbox z=[%.3f,%.3f] y=[%.3f,%.3f] "
                "x=[%.3f,%.3f]\n",
                lo3[0], hi3[0], lo3[1], hi3[1], lo3[2], hi3[2]);
    }
    if (m.nvt != m.nv) {
        fprintf(stderr, "  WARN: nvt (%zu) != nv (%zu) -- flat/raster outputs "
                "disabled\n", m.nvt, m.nv);
        do_flat = 0;
        do_raster = 0;
    }

    /* Raster preflight: size + cap check from the UV bbox alone, BEFORE the
     * cube table and the (minutes-long) vertex sampling pass.  The 21x21x21
     * chain sampled for 446s and only then learned the raster was rejected. */
    if (do_raster) {
        RawtexPlan rp;
        if (Rawtex_plan(m.uv, m.nv, raster_du, raster_dv,
                        raster_max_px, &rp) != 0) {
            fprintf(stderr, "ERROR: raster preflight failed (no UV?)\n");
            uvmesh_free(&m);
            return 1;
        }
        if (!rp.ok && raster_auto) {
            fprintf(stderr,
                "  raster preflight: %zux%zu at du=%.3g dv=%.3g exceeds caps; "
                "auto-adopting du=%.4g dv=%.4g\n",
                rp.W, rp.H, raster_du, raster_dv, rp.need_du, rp.need_dv);
            raster_du = rp.need_du;
            raster_dv = rp.need_dv;
            if (Rawtex_plan(m.uv, m.nv, raster_du, raster_dv,
                            raster_max_px, &rp) != 0 || !rp.ok) {
                fprintf(stderr, "ERROR: raster preflight advice does not fit "
                        "(%zux%zu, cap %zu px)\n", rp.W, rp.H, rp.max_px);
                uvmesh_free(&m);
                return 1;
            }
        }
        if (!rp.ok) {
            fprintf(stderr,
                "ERROR: raster %zux%zu unreasonable for uv u=[%.1f,%.1f] "
                "v=[%.1f,%.1f] at du=%.3g dv=%.3g; retry with --raster-du "
                "%.4g --raster-dv %.4g, or pass --raster-auto / "
                "--raster-max-px N (nothing sampled yet -- no time wasted)\n",
                rp.W, rp.H, rp.umin, rp.umax, rp.vmin, rp.vmax,
                raster_du, raster_dv, rp.need_du, rp.need_dv);
            uvmesh_free(&m);
            return 1;
        }
        fprintf(stderr, "  raster preflight: %zux%zu px at du=%.3g dv=%.3g "
                "(u=[%.1f,%.1f] v=[%.1f,%.1f], cap %zu px)\n",
                rp.W, rp.H, raster_du, raster_dv,
                rp.umin, rp.umax, rp.vmin, rp.vmax, rp.max_px);
    }

    arena = Arena_new();
    {
        double raw_pad = normal_range + 2.0;
        if (do_marble) {
            double marble_pad = marble_opts.depth_range
                              + 2.0 * marble_opts.tensor_radius + 2.0;
            if (marble_pad > raw_pad) raw_pad = marble_pad;
        }
        if (cubetable_init(&ct, arena, raw_dir, chunk, m.verts, m.nv,
                           raw_pad) != 0) {
            fprintf(stderr, "ERROR: cube table init failed\n");
            return 1;
        }
    }
    fprintf(stderr, "  cube table: %ldx%ldx%ld cubes (chunk %ld) over bbox\n",
            ct.nz, ct.ny, ct.nx, chunk);

    if (m.nf > 0)
        normals = MeshNormals_compute(m.verts, m.nv, m.faces, m.nf);
    need_vertex_samples = do_3d || do_flat || !fixed_window;
    if (need_vertex_samples) {
        val = (double *)xmalloc(m.nv * sizeof(double));
        has = (uint8_t *)xcalloc(m.nv, sizeof(uint8_t));
    }

    t0 = ves_clock_sec();
    /* prewarm freezes the cube table (raw_sample.h: reads are then OMP-safe),
     * so the 30M-vertex sample loop parallelizes cleanly */
    cubetable_prewarm_all(&ct);
    {
        size_t expected = cubetable_expected_chunks(&ct);
        char coverage_path[2300];
        FILE *coverage = NULL;
        snprintf(coverage_path, sizeof coverage_path,
                 "%s/%s_raw_coverage.json", out_dir, id);
        ves_ensure_parent_dir(coverage_path);
        coverage = fopen(coverage_path, "wb");
        if (coverage == NULL) {
            fprintf(stderr, "ERROR: cannot write RAW coverage report %s\n",
                    coverage_path);
            return 1;
        }
        fprintf(coverage,
                "{\n  \"schema\": \"obj-bake-raw-coverage-v1\",\n"
                "  \"source_kind\": \"%s\",\n"
                "  \"chunk_size\": %ld,\n"
                "  \"table_shape_zyx\": [%ld, %ld, %ld],\n"
                "  \"expected_chunks\": %zu,\n"
                "  \"loaded_chunks\": %d,\n"
                "  \"missing_chunks\": %d,\n"
                "  \"outside_volume_chunks\": %zu,\n"
                "  \"complete\": %s,\n"
                "  \"strict\": %s\n}\n",
                ct.is_zarr ? "zarr_v2" : "cube_tiffs", chunk,
                ct.nz, ct.ny, ct.nx, expected, ct.n_loaded, ct.n_missing,
                ct.n_outside, cubetable_is_complete(&ct) ? "true" : "false",
                require_complete_raw ? "true" : "false");
        fclose(coverage);
        fprintf(stderr,
                "  RAW coverage: %d/%zu in-volume chunks loaded, %d missing, "
                "%zu out-of-volume bbox chunks excluded (%.1fs); report=%s\n",
                ct.n_loaded, expected, ct.n_missing, ct.n_outside,
                ves_clock_sec() - t0, coverage_path);
        if (require_complete_raw && !cubetable_is_complete(&ct)) {
            fprintf(stderr,
                    "ERROR: incomplete RAW source for bake bbox: %d/%zu chunks "
                    "loaded, %d missing, %zu out-of-volume bbox chunks excluded "
                    "(source=%s). Refusing to sample or emit "
                    "a production texture. Use the authoritative full RAW source; "
                    "--allow-incomplete-raw is diagnostics-only.\n",
                    ct.n_loaded, expected, ct.n_missing, ct.n_outside, raw_dir);
            return 1;
        }
    }
    if (need_vertex_samples) {
        int kk = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (kk = 0; kk < (int)m.nv; kk++) {
            double s = sample_vertex(&ct, &m.verts[kk * 3],
                                     normals != NULL ? &normals[kk * 3] : NULL,
                                     normal_range, normal_samples);
            val[kk] = (s >= 0.0) ? s : 0.0;
            if (s >= 0.0) has[kk] = 1;
        }
    }
    if (need_vertex_samples)
        for (k = 0; k < m.nv; k++) nsamp += has[k];
    t_sample = ves_clock_sec() - t0;
    if (need_vertex_samples)
        fprintf(stderr, "  sampled: %zu/%zu verts (%.2f%%), %d cubes loaded, "
                "%d missing (%.1fs)\n", nsamp, m.nv,
                m.nv > 0 ? 100.0 * (double)nsamp / (double)m.nv : 0.0,
                ct.n_loaded, ct.n_missing, t_sample);
    else
        fprintf(stderr,
                "  skipped redundant per-vertex sampling; raster samples RAW "
                "directly with the fixed window (%.1fs)\n", t_sample);

    if (fixed_window) { lo = fixed_lo; hi = fixed_hi; }
    else if (do_stretch) stretch_window(val, has, m.nv, pct_lo, pct_hi, &lo, &hi);
    fprintf(stderr, "  contrast window: [%.0f, %.0f] (%s)\n",
            lo, hi, fixed_window ? "fixed" : (do_stretch ? "stretch on" : "stretch off"));
    if (need_vertex_samples) {   /* value distribution over sampled verts */
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
                              val, has, lo, hi, 0) != 0)
            fprintf(stderr, "  WARN: cannot write %s\n", path);
        else
            fprintf(stderr, "  wrote %s\n", path);
    }
    if (do_flat) {
        char path[2600];
        snprintf(path, sizeof path, "%s/%s_rawflat.obj", out_dir, id);
        ves_ensure_parent_dir(path);
        if (write_colored_obj(path, m.verts, m.nv, m.faces, m.nf, m.uv,
                              val, has, lo, hi, 1) != 0)
            fprintf(stderr, "  WARN: cannot write %s\n", path);
        else
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
                             raster_max_edge, raster_max_px,
                             do_diag ? &dopts : NULL,
                             &W, &H, &fill, &multi, &skuv, &sk3d) != 0)
            fprintf(stderr, "  WARN: cannot write %s\n", path);
        else {
            fprintf(stderr, "  wrote %s (%zux%zu px, %.1f%% filled, "
                    "%zu multi-cover px, skipped %zu uv-stretched + %zu "
                    "big-3d faces, %.1fs)\n",
                    path, W, H, 100.0 * fill, multi, skuv, sk3d,
                    ves_clock_sec() - tr);
            if (do_marble) {
                MarbleDiagStats mstats;
                if (MarbleDiag_write(dprefix, path, &ct,
                                     m.verts, m.uv, m.nv, m.faces, m.nf,
                                     normals, raster_du, raster_dv,
                                     raster_max_px, &marble_opts, &mstats) != 0)
                    fprintf(stderr, "  WARN: RAW/mesh marble audit failed\n");
                else
                    fprintf(stderr, "  wrote %s_marble_overlay.png\n", dprefix);
            }
        }
    }
    t_write = ves_clock_sec() - t0;

    fprintf(stderr, "[obj_bake_raw] SUMMARY id=%s verts=%zu sampled=%.2f%% "
            "cubes=%d/%d window=[%.0f,%.0f] parse=%.1fs sample=%.1fs "
            "write=%.1fs\n", id, m.nv,
            need_vertex_samples && m.nv > 0
                ? 100.0 * (double)nsamp / (double)m.nv : 0.0,
            ct.n_loaded, ct.n_loaded + ct.n_missing, lo, hi,
            t_parse, t_sample, t_write);

    free(normals); free(val); free(has);
    uvmesh_free(&m);
    Arena_dispose(&arena);
    return 0;
}
