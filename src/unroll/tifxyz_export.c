/* tifxyz_export.c -- unrolled strip -> VC3D tifxyz segment. See
 * tifxyz_export.h. The canvas mapping and face gates mirror scroll_raster.c
 * so the segment stays pixel-aligned with the stage review TIFs; coverage is
 * first-cover-wins (never a positional average -- overlap covers can
 * straddle wraps and an average would land in mid-air). */
#include "../common/ves_platform.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/tiff_io.h"
#include "tifxyz_export.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* classic (non-Big) TIFF caps a file at 4 GB; one float32 band is W*H*4 B.
 * Stay well clear so libtiff never truncates silently. */
#define TXZ_MAX_PX ((size_t)950000000)

void TifxyzOpts_default(TifxyzOpts *o)
{
    assert(o);
    memset(o, 0, sizeof(*o));
    o->du = 1.0;
    o->dv = 1.0;
    o->stretch_ratio = 4.0;
    o->stretch_floor = 25.0;
    o->max_edge3d = 0.0;
    o->conflict_dist = 2.0;
    o->write_provenance = 1;
}

/* px cover state (sticky max) */
enum { TXZ_EMPTY = 0, TXZ_ONE = 1, TXZ_MULTI = 2, TXZ_CONFLICT = 3 };

/* Per-face gate, identical semantics to scroll_raster: face_keep first, then
 * membrane (any 3D edge > max_edge3d), then uv stretch. NaN anywhere in a
 * face's verts/uv is gated as 3d (it is broken geometry). Returns 0 = keep,
 * else the stat bucket it was counted into. */
static int txz_gate_face(const PieceSet *ps, const TifxyzOpts *opts, size_t f,
                         TifxyzStats *st)
{
    if (opts->face_keep != NULL && !opts->face_keep[f]) {
        st->skip_own++;
        return 1;
    }
    size_t vi[3];
    vi[0] = (size_t)ps->faces[f * 3 + 0];
    vi[1] = (size_t)ps->faces[f * 3 + 1];
    vi[2] = (size_t)ps->faces[f * 3 + 2];
    int bad_uv = 0, bad_3d = 0;
    for (int e = 0; e < 3; e++) {
        size_t p = vi[e], q = vi[(e + 1) % 3];
        double d3 = 0.0;
        for (int k = 0; k < 3; k++) {
            double dd = (double)ps->verts[q * 3 + k]
                        - (double)ps->verts[p * 3 + k];
            d3 += dd * dd;
        }
        d3 = sqrt(d3);
        double due = (double)ps->uv[q * 2 + 0] - (double)ps->uv[p * 2 + 0];
        double dve = (double)ps->uv[q * 2 + 1] - (double)ps->uv[p * 2 + 1];
        double e2d = sqrt(due * due + dve * dve);
        if (!(d3 == d3) || !(e2d == e2d)) bad_3d = 1;   /* NaN */
        if (opts->max_edge3d > 0.0 && d3 > opts->max_edge3d) bad_3d = 1;
        if (opts->stretch_ratio > 0.0
            && e2d > (opts->stretch_ratio * d3 > opts->stretch_floor
                      ? opts->stretch_ratio * d3 : opts->stretch_floor))
            bad_uv = 1;
    }
    if (bad_3d) { st->skip_3d++; return 1; }
    if (bad_uv) { st->skip_uv++; return 1; }
    return 0;
}

typedef struct {
    const PieceSet *ps;
    const TifxyzOpts *opts;
    long X0, Y0;
    size_t W, H;
    float *px, *py, *pz;      /* [W*H] world x/y/z planes, -1 = invalid */
    float *pw;                /* [W*H] winding (phi/2pi), NULL = off */
    uint8_t *state;           /* [W*H] TXZ_* cover state */
} TxzCanvas;

/* uv (vox) -> canvas px coords, flips applied. */
static void txz_px(const TxzCanvas *cv, size_t v, double *out_x, double *out_y)
{
    const TifxyzOpts *o = cv->opts;
    double fx = (double)cv->ps->uv[v * 2 + 0] / o->du - (double)cv->X0;
    double fy = (double)cv->ps->uv[v * 2 + 1] / o->dv - (double)cv->Y0;
    if (o->flip_u) fx = (double)(cv->W - 1) - fx;
    if (o->flip_v) fy = (double)(cv->H - 1) - fy;
    *out_x = fx;
    *out_y = fy;
}

/* Splat one kept face; first cover writes the px, later covers only bump
 * multi/conflict state. Same bbox/epsilon/inside rules as sr_splat (so px
 * on a shared mesh edge count as multi, exactly like the raster metric). */
static void txz_splat_face(TxzCanvas *cv, size_t f)
{
    const PieceSet *ps = cv->ps;
    size_t a = (size_t)ps->faces[f * 3 + 0];
    size_t b = (size_t)ps->faces[f * 3 + 1];
    size_t c = (size_t)ps->faces[f * 3 + 2];
    double ua = 0.0, va = 0.0, ub = 0.0, vb = 0.0, uc = 0.0, vc = 0.0;
    txz_px(cv, a, &ua, &va);
    txz_px(cv, b, &ub, &vb);
    txz_px(cv, c, &uc, &vc);
    double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
    if (fabs(A2) < 1e-12) return;
    double lox = ua < ub ? (ua < uc ? ua : uc) : (ub < uc ? ub : uc);
    double hix = ua > ub ? (ua > uc ? ua : uc) : (ub > uc ? ub : uc);
    double loy = va < vb ? (va < vc ? va : vc) : (vb < vc ? vb : vc);
    double hiy = va > vb ? (va > vc ? va : vc) : (vb > vc ? vb : vc);
    long x0 = (long)ceil(lox - 0.001), x1 = (long)floor(hix + 0.001);
    long y0 = (long)ceil(loy - 0.001), y1 = (long)floor(hiy + 0.001);
    if (x0 < 0) x0 = 0;
    if (x1 >= (long)cv->W) x1 = (long)cv->W - 1;
    if (y0 < 0) y0 = 0;
    if (y1 >= (long)cv->H) y1 = (long)cv->H - 1;
    for (long yy = y0; yy <= y1; yy++) {
        for (long xx = x0; xx <= x1; xx++) {
            double pu = (double)xx, pv = (double)yy;
            double l1 = ((pu - ua) * (vc - va) - (pv - va) * (uc - ua)) / A2;
            double l2 = ((ub - ua) * (pv - va) - (vb - va) * (pu - ua)) / A2;
            double l0 = 1.0 - l1 - l2;
            if (l0 < -1e-6 || l1 < -1e-6 || l2 < -1e-6) continue;
            size_t pi = (size_t)yy * cv->W + (size_t)xx;
            /* native (z,y,x) -> world (x,y,z) */
            double wx = l0 * (double)ps->verts[a * 3 + 2]
                        + l1 * (double)ps->verts[b * 3 + 2]
                        + l2 * (double)ps->verts[c * 3 + 2];
            double wy = l0 * (double)ps->verts[a * 3 + 1]
                        + l1 * (double)ps->verts[b * 3 + 1]
                        + l2 * (double)ps->verts[c * 3 + 1];
            double wz = l0 * (double)ps->verts[a * 3 + 0]
                        + l1 * (double)ps->verts[b * 3 + 0]
                        + l2 * (double)ps->verts[c * 3 + 0];
            if (cv->state[pi] == TXZ_EMPTY) {
                cv->px[pi] = (float)wx;
                cv->py[pi] = (float)wy;
                cv->pz[pi] = (float)wz;
                if (cv->pw != NULL && ps->phi != NULL) {
                    double phi = l0 * (double)ps->phi[a]
                                 + l1 * (double)ps->phi[b]
                                 + l2 * (double)ps->phi[c];
                    cv->pw[pi] = (float)(phi / (2.0 * M_PI));
                }
                cv->state[pi] = TXZ_ONE;
            } else {
                double dx = wx - (double)cv->px[pi];
                double dy = wy - (double)cv->py[pi];
                double dz = wz - (double)cv->pz[pi];
                double d = sqrt(dx * dx + dy * dy + dz * dz);
                if (d > cv->opts->conflict_dist)
                    cv->state[pi] = TXZ_CONFLICT;
                else if (cv->state[pi] == TXZ_ONE)
                    cv->state[pi] = TXZ_MULTI;
            }
        }
    }
}

static int txz_write_meta(const char *seg_dir, const char *uuid,
                          const PieceSet *ps, const TifxyzOpts *opts,
                          const TifxyzStats *st)
{
    char path[2048];
    snprintf(path, sizeof(path), "%s/meta.json", seg_dir);
    FILE *f = fopen(path, "w");
    if (f == NULL) return -1;
    fprintf(f,
        "{\n"
        "    \"scale\": [%.8f, %.8f],\n"
        "    \"bbox\": [[%.3f, %.3f, %.3f], [%.3f, %.3f, %.3f]],\n"
        "    \"type\": \"seg\",\n"
        "    \"format\": \"tifxyz\",\n"
        "    \"uuid\": \"%s\",\n"
        "    \"source\": \"scrollfiesta\",\n"
        "    \"area_vx2\": %.1f,\n"
        "    \"scrollfiesta\": {\n"
        "        \"du\": %.6f, \"dv\": %.6f, \"flip_u\": %d, \"flip_v\": %d,\n"
        "        \"origin_uv\": [%.3f, %.3f],\n"
        "        \"n_cubes\": %zu,\n"
        "        \"valid_px\": %zu, \"multi_px\": %zu, \"conflict_px\": %zu,\n"
        "        \"skip_uv_faces\": %zu, \"skip_3d_faces\": %zu, "
        "\"skip_own_faces\": %zu%s\n"
        "    }\n"
        "}\n",
        1.0 / opts->du, 1.0 / opts->dv,
        st->bbox_lo[0], st->bbox_lo[1], st->bbox_lo[2],
        st->bbox_hi[0], st->bbox_hi[1], st->bbox_hi[2],
        uuid, st->area_vx2,
        opts->du, opts->dv, opts->flip_u ? 1 : 0, opts->flip_v ? 1 : 0,
        st->u0, st->v0, ps->n_cubes,
        st->valid, st->multi, st->conflicts,
        st->skip_uv, st->skip_3d, st->skip_own,
        opts->write_provenance
            ? ",\n        \"provenance\": \"provenance.tif: 0=empty 1=single "
              "2=overlap 3=contested\""
            : "");
    int bad = ferror(f) != 0;
    if (fclose(f) != 0) bad = 1;
    return bad ? -1 : 0;
}

int TifxyzExport_run(Arena_T arena, const PieceSet *ps, const char *seg_dir,
                     const char *uuid, const TifxyzOpts *opts,
                     TifxyzStats *out)
{
    assert(arena && ps && seg_dir && uuid && opts && out);
    (void)arena;   /* planes are transient multi-GB: malloc/free like
                      scroll_raster's canvases, never arena-resident */
    memset(out, 0, sizeof(*out));
    if (ps->nf == 0 || ps->nv == 0 || opts->du <= 0.0 || opts->dv <= 0.0) {
        fprintf(stderr, "tifxyz_export: ERROR empty piece set or bad step\n");
        return -1;
    }
    if (uuid[0] == '\0') {
        fprintf(stderr, "tifxyz_export: ERROR empty uuid\n");
        return -1;
    }

    /* canvas mapping identical to scroll_raster: px x = round(u/du) - X0.
     * A forced window crops to one atlas piece; without it the canvas spans
     * the whole piece set. Either way X0/Y0 come from floor(u/du+0.5) so all
     * pieces share one integer lattice. */
    const double du = opts->du, dv = opts->dv;
    double u_lo = opts->have_window ? opts->win_u0 : ps->u_min;
    double u_hi = opts->have_window ? opts->win_u1 : ps->u_max;
    double v_lo = opts->have_window ? opts->win_v0 : ps->v_min;
    double v_hi = opts->have_window ? opts->win_v1 : ps->v_max;
    long X0 = (long)floor(u_lo / du + 0.5);
    long Y0 = (long)floor(v_lo / dv + 0.5);
    long X1 = (long)floor(u_hi / du + 0.5);
    long Y1 = (long)floor(v_hi / dv + 0.5);
    size_t W = (size_t)(X1 - X0) + 1;
    size_t H = (size_t)(Y1 - Y0) + 1;
    if (W < 1 || H < 1 || W > (size_t)1 << 31 || H > (size_t)1 << 24) {
        fprintf(stderr, "tifxyz_export: ERROR canvas %zux%zu unreasonable "
                "(raise du/dv?)\n", W, H);
        return -1;
    }
    if (W * H > TXZ_MAX_PX) {
        fprintf(stderr, "tifxyz_export: ERROR canvas %zux%zu = %zu px "
                "exceeds the classic-TIFF 4 GB band limit (%zu px) -- "
                "raise du/dv\n", W, H, W * H, TXZ_MAX_PX);
        return -1;
    }
    out->W = W;
    out->H = H;
    out->u0 = (double)X0 * du;
    out->v0 = (double)Y0 * dv;

    const size_t npx = W * H;
    int want_winding = opts->write_winding && ps->phi != NULL;
    float *px = (float *)malloc(npx * sizeof(float));
    float *py = (float *)malloc(npx * sizeof(float));
    float *pz = (float *)malloc(npx * sizeof(float));
    float *pw = want_winding ? (float *)malloc(npx * sizeof(float)) : NULL;
    uint8_t *state = (uint8_t *)calloc(npx, 1);
    if (px == NULL || py == NULL || pz == NULL || state == NULL ||
        (want_winding && pw == NULL)) {
        fprintf(stderr, "tifxyz_export: ERROR out of memory "
                "(canvas %zux%zu)\n", W, H);
        free(px); free(py); free(pz); free(pw); free(state);
        return -1;
    }
    for (size_t i = 0; i < npx; i++) {
        px[i] = -1.0f;
        py[i] = -1.0f;
        pz[i] = -1.0f;
    }
    if (pw != NULL)
        memset(pw, 0, npx * sizeof(float));

    TxzCanvas cv;
    cv.ps = ps;
    cv.opts = opts;
    cv.X0 = X0;
    cv.Y0 = Y0;
    cv.W = W;
    cv.H = H;
    cv.px = px;
    cv.py = py;
    cv.pz = pz;
    cv.pw = pw;
    cv.state = state;

    /* ascending face order => first-cover-wins is deterministic */
    for (size_t f = 0; f < ps->nf; f++) {
        if (txz_gate_face(ps, opts, f, out) != 0) continue;
        txz_splat_face(&cv, f);
    }

    /* stats + world bbox (x,y,z order) over covered px */
    double lo[3] = { 1e300, 1e300, 1e300 };
    double hi[3] = { -1e300, -1e300, -1e300 };
    for (size_t i = 0; i < npx; i++) {
        if (state[i] == TXZ_EMPTY) continue;
        out->valid++;
        if (state[i] >= TXZ_MULTI) out->multi++;
        if (state[i] == TXZ_CONFLICT) out->conflicts++;
        double w3[3];
        w3[0] = (double)px[i];
        w3[1] = (double)py[i];
        w3[2] = (double)pz[i];
        for (int k = 0; k < 3; k++) {
            if (w3[k] < lo[k]) lo[k] = w3[k];
            if (w3[k] > hi[k]) hi[k] = w3[k];
        }
    }
    if (out->valid == 0) {
        fprintf(stderr, "tifxyz_export: ERROR zero covered pixels "
                "(all faces gated?)\n");
        free(px); free(py); free(pz); free(pw); free(state);
        return -1;
    }
    for (int k = 0; k < 3; k++) {
        out->bbox_lo[k] = lo[k];
        out->bbox_hi[k] = hi[k];
    }

    /* area = valid 2x2 quads * du*dv (the python reference's area_vx2) */
    {
        size_t quads = 0;
        for (size_t y = 0; y + 1 < H; y++)
            for (size_t x = 0; x + 1 < W; x++)
                if (state[y * W + x] != TXZ_EMPTY &&
                    state[y * W + x + 1] != TXZ_EMPTY &&
                    state[(y + 1) * W + x] != TXZ_EMPTY &&
                    state[(y + 1) * W + x + 1] != TXZ_EMPTY)
                    quads++;
        out->area_vx2 = (double)quads * du * dv;
    }

    /* write the segment: dir, bands, mask, meta */
    int rc = 0;
    char path[2048];
    snprintf(path, sizeof(path), "%s/x.tif", seg_dir);
    if (ves_ensure_parent_dir(path) != 0) {
        fprintf(stderr, "tifxyz_export: ERROR cannot create %s\n", seg_dir);
        rc = -1;
    }
    if (rc == 0 && TiffIO_save_float2d(path, px, (int)W, (int)H) != 0) {
        fprintf(stderr, "tifxyz_export: ERROR writing %s\n", path);
        rc = -1;
    }
    if (rc == 0) {
        snprintf(path, sizeof(path), "%s/y.tif", seg_dir);
        if (TiffIO_save_float2d(path, py, (int)W, (int)H) != 0) {
            fprintf(stderr, "tifxyz_export: ERROR writing %s\n", path);
            rc = -1;
        }
    }
    if (rc == 0) {
        snprintf(path, sizeof(path), "%s/z.tif", seg_dir);
        if (TiffIO_save_float2d(path, pz, (int)W, (int)H) != 0) {
            fprintf(stderr, "tifxyz_export: ERROR writing %s\n", path);
            rc = -1;
        }
    }
    if (rc == 0 && pw != NULL) {
        snprintf(path, sizeof(path), "%s/winding.tif", seg_dir);
        if (TiffIO_save_float2d(path, pw, (int)W, (int)H) != 0) {
            fprintf(stderr, "tifxyz_export: ERROR writing %s\n", path);
            rc = -1;
        }
    }
    if (rc == 0 && opts->write_provenance) {
        /* provenance.tif: per-px cover state (0 empty / 1 single / 2 overlap /
         * 3 contested) -- emitted BEFORE the mask block below reuses `state`
         * in place. Marks the winding-collapse core (state 3 = two covers
         * > conflict_dist apart in 3D) WITHOUT dropping it: the mask stays 255
         * on those px, so nothing is hidden -- the contest is recorded, not
         * silently resolved. */
        snprintf(path, sizeof(path), "%s/provenance.tif", seg_dir);
        if (TiffIO_save(path, state, 1, (int)H, (int)W) != 0) {
            fprintf(stderr, "tifxyz_export: ERROR writing %s\n", path);
            rc = -1;
        }
    }
    if (rc == 0) {
        /* mask: 255 = keep (VC3D invalidates any px whose mask < 255) */
        uint8_t *mask = state;   /* reuse in place */
        for (size_t i = 0; i < npx; i++)
            mask[i] = mask[i] != TXZ_EMPTY ? (uint8_t)255 : (uint8_t)0;
        snprintf(path, sizeof(path), "%s/mask.tif", seg_dir);
        if (TiffIO_save(path, mask, 1, (int)H, (int)W) != 0) {
            fprintf(stderr, "tifxyz_export: ERROR writing %s\n", path);
            rc = -1;
        }
    }
    if (rc == 0 && txz_write_meta(seg_dir, uuid, ps, opts, out) != 0) {
        fprintf(stderr, "tifxyz_export: ERROR writing %s/meta.json\n",
                seg_dir);
        rc = -1;
    }

    free(px);
    free(py);
    free(pz);
    free(pw);
    free(state);
    return rc;
}

/* ============================================================================
 * Self-test
 * ==========================================================================*/

static int txz_check(int cond, const char *what, int *fails)
{
    if (!cond) {
        fprintf(stderr, "[tifxyz selftest]   FAIL: %s\n", what);
        (*fails)++;
    }
    return cond;
}

/* the scroll_raster selftest quad: 20x20-vox flat sheet at y=20 with
 * u = x - 1010 (so world x = u + 1010), v = z, absolute-negative uv frame */
static void txz_fixture(PieceSet *ps, float *verts, float *uv, float *phi,
                        float *normals, int32_t *faces, int32_t *face_cube)
{
    static const float V[4 * 3] = {
        10, 20, 10,   30, 20, 10,   10, 20, 30,   30, 20, 30
    };
    static const float U[4 * 2] = {
        -1000, 10,   -1000, 30,   -980, 10,   -980, 30
    };
    static const int32_t F[2 * 3] = { 0, 1, 2, 1, 3, 2 };
    memset(ps, 0, sizeof(*ps));
    memcpy(verts, V, sizeof(V));
    memcpy(uv, U, sizeof(U));
    /* phi = 2pi * (x-10)/20 -> winding = (u px)/20 */
    for (int i = 0; i < 4; i++) {
        phi[i] = (float)(2.0 * M_PI * ((double)verts[i * 3 + 2] - 10.0) / 20.0);
        normals[i * 3 + 0] = 0;
        normals[i * 3 + 1] = 1;
        normals[i * 3 + 2] = 0;
    }
    memcpy(faces, F, sizeof(F));
    face_cube[0] = 0;
    face_cube[1] = 0;
    ps->verts = verts;
    ps->uv = uv;
    ps->phi = phi;
    ps->normals = normals;
    ps->faces = faces;
    ps->face_cube = face_cube;
    ps->nv = 4;
    ps->nf = 2;
    ps->n_cubes = 1;
    ps->u_min = -1000; ps->u_max = -980;
    ps->v_min = 10; ps->v_max = 30;
}

/* slurp a small text file; returns arena string or NULL */
static char *txz_slurp(Arena_T arena, const char *path)
{
    FILE *f = fopen(path, "rb");
    if (f == NULL) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return NULL; }
    long n = ftell(f);
    if (n < 0 || n > 1000000) { fclose(f); return NULL; }
    rewind(f);
    char *buf = (char *)ARENA_ALLOC(arena, (size_t)n + 1);
    size_t got = fread(buf, 1, (size_t)n, f);
    fclose(f);
    buf[got] = '\0';
    return buf;
}

int TifxyzExport_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();
    float verts[4 * 3], uv[4 * 2], phi[4], normals[4 * 3];
    int32_t faces[2 * 3], face_cube[2];
    const char *base = "output/_selftest_scroll_unroll";

    /* ---- A: flat quad roundtrip (axis swap, mask, winding, meta) --------- */
    {
        PieceSet ps;
        txz_fixture(&ps, verts, uv, phi, normals, faces, face_cube);
        TifxyzOpts o;
        TifxyzOpts_default(&o);
        o.max_edge3d = 100.0;   /* fixture edges are 20 vox */
        o.write_winding = 1;
        TifxyzStats st;
        char dir[512];
        snprintf(dir, sizeof(dir), "%s/tifxyz_a", base);
        int rc = TifxyzExport_run(arena, &ps, dir, "tifxyz_a", &o, &st);
        txz_check(rc == 0, "A: rc", &fails);
        txz_check(st.W == 21 && st.H == 21, "A: canvas dims", &fails);
        /* the quad's two triangles share the 21-px diagonal (inclusive
         * barycentric edge, same as the raster) -> multi = 21, 0 conflicts */
        txz_check(st.valid == 441 && st.multi == 21 && st.conflicts == 0,
                  "A: full cover, diagonal-only multi", &fails);
        txz_check(fabs(st.area_vx2 - 400.0) < 1e-9, "A: area 20x20", &fails);
        txz_check(fabs(st.bbox_lo[0] - 10.0) < 1e-3
                  && fabs(st.bbox_hi[0] - 30.0) < 1e-3
                  && fabs(st.bbox_lo[1] - 20.0) < 1e-3
                  && fabs(st.bbox_hi[1] - 20.0) < 1e-3
                  && fabs(st.bbox_lo[2] - 10.0) < 1e-3
                  && fabs(st.bbox_hi[2] - 30.0) < 1e-3,
                  "A: bbox in (x,y,z) order", &fails);

        float *bx = NULL, *by = NULL, *bz = NULL, *bw = NULL;
        int w = 0, h = 0;
        char p[600];
        snprintf(p, sizeof(p), "%s/x.tif", dir);
        int ok = TiffIO_load_float2d(arena, p, &bx, &w, &h) == 0;
        snprintf(p, sizeof(p), "%s/y.tif", dir);
        ok = ok && TiffIO_load_float2d(arena, p, &by, &w, &h) == 0;
        snprintf(p, sizeof(p), "%s/z.tif", dir);
        ok = ok && TiffIO_load_float2d(arena, p, &bz, &w, &h) == 0;
        snprintf(p, sizeof(p), "%s/winding.tif", dir);
        ok = ok && TiffIO_load_float2d(arena, p, &bw, &w, &h) == 0;
        txz_check(ok && w == 21 && h == 21, "A: band readback", &fails);
        if (ok && w == 21 && h == 21) {
            /* px (row cy, col cx): u = -1000+cx, v = 10+cy
             * => world x = u+1010 = 10+cx, y = 20, z = v = 10+cy */
            int exact = 1;
            static const int probe[4][2] = {
                { 0, 0 }, { 10, 10 }, { 20, 5 }, { 3, 17 }
            };
            for (int t = 0; t < 4; t++) {
                int cy = probe[t][0], cx = probe[t][1];
                size_t i = (size_t)cy * 21 + (size_t)cx;
                if (fabs((double)bx[i] - (10.0 + cx)) > 1e-3) exact = 0;
                if (fabs((double)by[i] - 20.0) > 1e-3) exact = 0;
                if (fabs((double)bz[i] - (10.0 + cy)) > 1e-3) exact = 0;
            }
            txz_check(exact, "A: x/y/z axis swap + interpolation exact",
                      &fails);
            size_t ci = (size_t)10 * 21 + 10;
            txz_check(fabs((double)bw[ci] - 0.5) < 1e-3,
                      "A: winding = phi/2pi", &fails);
        }
        uint8_t *mask = NULL;
        int md = 0, mh = 0, mw = 0;
        snprintf(p, sizeof(p), "%s/mask.tif", dir);
        if (txz_check(TiffIO_load(arena, p, &mask, &md, &mh, &mw) == 0,
                      "A: mask readback", &fails)) {
            int all255 = (mh == 21 && mw == 21);
            for (long i = 0; all255 && i < 21 * 21; i++)
                if (mask[i] != 255) all255 = 0;
            txz_check(all255, "A: mask all 255", &fails);
        }
        /* provenance.tif: single cover = 1, the 21 shared-diagonal px = 2
         * (overlap within conflict_dist), no contested px. Ties the band to
         * st.multi/st.conflicts exactly. */
        uint8_t *prov = NULL;
        int pd = 0, ph = 0, pw = 0;
        snprintf(p, sizeof(p), "%s/provenance.tif", dir);
        if (txz_check(TiffIO_load(arena, p, &prov, &pd, &ph, &pw) == 0,
                      "A: provenance readback", &fails)) {
            size_t n0 = 0, n1 = 0, n2 = 0, n3 = 0;
            for (long i = 0; i < 21 * 21; i++)
                switch (prov[i]) {
                    case 0: n0++; break; case 1: n1++; break;
                    case 2: n2++; break; default: n3++; break;
                }
            txz_check(pw == 21 && ph == 21 && n0 == 0 && n1 == 420
                      && n2 == 21 && n3 == 0,
                      "A: provenance 1=single(420) 2=overlap(21) 3=none",
                      &fails);
        }
        snprintf(p, sizeof(p), "%s/meta.json", dir);
        char *meta = txz_slurp(arena, p);
        txz_check(meta != NULL
                  && strstr(meta, "\"scale\": [1.00000000, 1.00000000]") != NULL
                  && strstr(meta, "\"format\": \"tifxyz\"") != NULL
                  && strstr(meta, "\"type\": \"seg\"") != NULL
                  && strstr(meta, "\"provenance\": \"provenance.tif") != NULL
                  && strstr(meta, "\"uuid\": \"tifxyz_a\"") != NULL,
                  "A: meta.json required fields", &fails);
    }

    /* ---- B: face_keep drop -> partial cover, -1 sentinel, mask 0 --------- */
    {
        PieceSet ps;
        txz_fixture(&ps, verts, uv, phi, normals, faces, face_cube);
        TifxyzOpts o;
        TifxyzOpts_default(&o);
        o.max_edge3d = 100.0;
        uint8_t keep[2] = { 1, 0 };
        o.face_keep = keep;
        TifxyzStats st;
        char dir[512];
        snprintf(dir, sizeof(dir), "%s/tifxyz_b", base);
        int rc = TifxyzExport_run(arena, &ps, dir, "tifxyz_b", &o, &st);
        txz_check(rc == 0 && st.skip_own == 1, "B: face_keep counted", &fails);
        txz_check(st.valid > 0 && st.valid < 441, "B: partial cover", &fails);
        float *bx = NULL;
        int w = 0, h = 0;
        char p[600];
        snprintf(p, sizeof(p), "%s/x.tif", dir);
        if (txz_check(TiffIO_load_float2d(arena, p, &bx, &w, &h) == 0,
                      "B: x readback", &fails)) {
            /* kept face 0 = px corners (0,0),(0,20),(20,0): covers
             * cx+cy <= 20. Strictly-outside probe = (row 20, col 20);
             * strictly-inside probe = (0, 0). */
            size_t i_out = (size_t)20 * (size_t)w + 20;
            size_t i_in = 0;
            txz_check(bx[i_out] == -1.0f, "B: uncovered px = -1", &fails);
            txz_check(bx[i_in] != -1.0f, "B: covered px real", &fails);
        }
        uint8_t *mask = NULL;
        int md = 0, mh = 0, mw = 0;
        snprintf(p, sizeof(p), "%s/mask.tif", dir);
        if (TiffIO_load(arena, p, &mask, &md, &mh, &mw) == 0) {
            txz_check(mask[(size_t)20 * (size_t)mw + 20] == 0
                      && mask[0] == 255,
                      "B: mask matches cover", &fails);
        } else {
            txz_check(0, "B: mask readback", &fails);
        }
    }

    /* ---- C: flip_u mirrors columns --------------------------------------- */
    {
        PieceSet ps;
        txz_fixture(&ps, verts, uv, phi, normals, faces, face_cube);
        TifxyzOpts o;
        TifxyzOpts_default(&o);
        o.max_edge3d = 100.0;
        o.flip_u = 1;
        TifxyzStats st;
        char dir[512];
        snprintf(dir, sizeof(dir), "%s/tifxyz_c", base);
        int rc = TifxyzExport_run(arena, &ps, dir, "tifxyz_c", &o, &st);
        txz_check(rc == 0 && st.valid == 441, "C: rc + full cover", &fails);
        float *bx = NULL;
        int w = 0, h = 0;
        char p[600];
        snprintf(p, sizeof(p), "%s/x.tif", dir);
        if (txz_check(TiffIO_load_float2d(arena, p, &bx, &w, &h) == 0,
                      "C: x readback", &fails)) {
            /* unflipped x(cy,cx) = 10+cx; flipped = 10+(W-1-cx) */
            int exact = 1;
            for (int cx = 0; cx < 21 && exact; cx += 5) {
                size_t i = (size_t)10 * 21 + (size_t)cx;
                if (fabs((double)bx[i] - (10.0 + (20 - cx))) > 1e-3) exact = 0;
            }
            txz_check(exact, "C: flip_u mirrors world x", &fails);
        }
    }

    /* ---- D: coincident-uv second face far in 3D -> conflict, first wins -- */
    {
        float v2[8 * 3], u2[8 * 2], ph2[8], n2[8 * 3];
        int32_t f2[4 * 3], fc2[4];
        PieceSet ps;
        txz_fixture(&ps, v2, u2, ph2, n2, f2, fc2);
        /* duplicate the quad 10 vox away in y, same uv */
        for (int i = 0; i < 4; i++) {
            v2[(4 + i) * 3 + 0] = v2[i * 3 + 0];
            v2[(4 + i) * 3 + 1] = v2[i * 3 + 1] + 10.0f;
            v2[(4 + i) * 3 + 2] = v2[i * 3 + 2];
            u2[(4 + i) * 2 + 0] = u2[i * 2 + 0];
            u2[(4 + i) * 2 + 1] = u2[i * 2 + 1];
            ph2[4 + i] = ph2[i];
            n2[(4 + i) * 3 + 0] = 0;
            n2[(4 + i) * 3 + 1] = 1;
            n2[(4 + i) * 3 + 2] = 0;
        }
        for (int i = 0; i < 6; i++) f2[6 + i] = f2[i] + 4;
        fc2[2] = 0; fc2[3] = 0;
        ps.nv = 8;
        ps.nf = 4;
        TifxyzOpts o;
        TifxyzOpts_default(&o);
        o.max_edge3d = 100.0;
        TifxyzStats st;
        char dir[512];
        snprintf(dir, sizeof(dir), "%s/tifxyz_d", base);
        int rc = TifxyzExport_run(arena, &ps, dir, "tifxyz_d", &o, &st);
        txz_check(rc == 0, "D: rc", &fails);
        txz_check(st.valid == 441 && st.multi == 441 && st.conflicts == 441,
                  "D: all px conflicted", &fails);
        float *by = NULL;
        int w = 0, h = 0;
        char p[600];
        snprintf(p, sizeof(p), "%s/y.tif", dir);
        if (txz_check(TiffIO_load_float2d(arena, p, &by, &w, &h) == 0,
                      "D: y readback", &fails)) {
            size_t ci = (size_t)10 * 21 + 10;
            txz_check(fabs((double)by[ci] - 20.0) < 1e-3,
                      "D: first cover wins (y=20 not 30)", &fails);
        }
        uint8_t *prov = NULL;
        int pd = 0, ph = 0, pw = 0;
        snprintf(p, sizeof(p), "%s/provenance.tif", dir);
        if (txz_check(TiffIO_load(arena, p, &prov, &pd, &ph, &pw) == 0,
                      "D: provenance readback", &fails)) {
            int all3 = (pw == 21 && ph == 21);
            for (size_t i = 0; i < (size_t)pw * (size_t)ph; i++)
                if (prov[i] != 3) all3 = 0;
            txz_check(all3, "D: provenance all contested (3)", &fails);
        }
    }

    /* ---- E: degenerate inputs error cleanly ------------------------------ */
    {
        PieceSet ps;
        txz_fixture(&ps, verts, uv, phi, normals, faces, face_cube);
        TifxyzOpts o;
        TifxyzOpts_default(&o);
        o.max_edge3d = 100.0;
        TifxyzStats st;
        char dir[512];
        snprintf(dir, sizeof(dir), "%s/tifxyz_e", base);

        PieceSet empty = ps;
        empty.nf = 0;
        txz_check(TifxyzExport_run(arena, &empty, dir, "tifxyz_e", &o,
                                   &st) != 0,
                  "E: zero faces -> error", &fails);

        uint8_t keep_none[2] = { 0, 0 };
        o.face_keep = keep_none;
        txz_check(TifxyzExport_run(arena, &ps, dir, "tifxyz_e", &o, &st) != 0
                  && st.skip_own == 2,
                  "E: all faces dropped -> error", &fails);
        o.face_keep = NULL;

        /* Real/source geometry is not length-gated by default: its edge scale
         * follows the upstream remesher. The cap remains explicitly opt-in. */
        TifxyzOpts oreal;
        TifxyzOpts_default(&oreal);
        txz_check(TifxyzExport_run(arena, &ps, dir, "tifxyz_e", &oreal,
                                   &st) == 0
                  && st.skip_3d == 0 && st.valid > 0,
                  "E: default keeps coarse real/source faces", &fails);

        TifxyzOpts o6;
        TifxyzOpts_default(&o6);
        o6.max_edge3d = 6.0;
        txz_check(TifxyzExport_run(arena, &ps, dir, "tifxyz_e", &o6,
                                   &st) != 0
                  && st.skip_3d == 2,
                  "E: explicit real edge cap -> zero cover error", &fails);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[tifxyz selftest] %s (%d failure%s)\n",
            fails == 0 ? "PASSED" : "FAILED", fails, fails == 1 ? "" : "s");
    return fails;
}
