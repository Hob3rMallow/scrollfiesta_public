/* ============================================================================
 * solid_quad_ribbon.c -- assemble the SOLID quad-ribbon deliverable.
 *
 * One solid quad strip per scroll component (a track, or a registration
 * super-component = several fragments of one wrap): a full num_v x num_u lattice
 * whose positions use the dense fitted/PDE field as a continuous initializer.
 * Production trust is finite, unambiguous direct support with nonzero fitted
 * confidence and a universal-cover phase; the legacy certified-mask policy is
 * retained as an explicit alternate mode.  Everything else remains movable,
 * and coherent high-distortion patches are excised and rebuilt by the
 * winding-aware cylindrical continuation.  Components are never bridged; the
 * empty inter-component domain is never filled -- no inter-wrap merger.
 *
 * Reads the same inputs as the Python assembler (assemble_solid_ribbon.py):
 *   <aggregate>/report.json                        bbox_zyx[0] = z_lo
 *   <aggregate>/material_tracks/manifest.json       tracks[].global_track,
 *                                                   tracks[].members[0].global_u_offset
 *   <tracks_root>/track_XXXXX/tifxyz/{field_z,field_y,field_x,mask,support,
 *                           confidence,interpolation_domain,reference_phase}.tif
 *   <tracks_root>/track_XXXXX/geometry_report.json  origin_zu[0]
 * Writes <out.obj> AND its <out.vmesh> companion (VESMESH1, authoritative for
 * the bake), plus an optional per-component topology/coverage PNG.
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/tiff_io.h"
#include "../common/obj_io.h"
#include "../common/mesh_bin.h"
#include "../common/ves_png.h"
#include "../flatten/quad_strip.h"
#include "../flatten/quad_strip_snap.h"
#include "../remesh/intersection_cleanup.h"

/* ---- one loaded component fragment (a manifest track) -------------------- */
typedef struct {
    long     tid;
    uint8_t *topology;/* [H*W] finite support/interpolation surface footprint */
    uint8_t *domain;  /* [H*W] finite initializer samples; may be quarantined */
    uint8_t *trusted; /* [H*W] high-confidence direct observations only */
    float   *field;   /* [H*W*3] (z,y,x) */
    float   *phase;   /* [H*W] upstream lifted phase; finite on domain */
    int      H, W;
    int      gr0, c0; /* lattice origin: u=(c0+i)*grid_du, v=gr0+j */
} Loaded;

/* Restrict one fitted track to a half-open range of global material-U lattice
 * columns.  This is deliberately done before initializer quarantine and the
 * cylindrical solve: a debug slice is its own boundary-value problem, not a
 * crop of an already solved/self-intersecting mesh. */
static int crop_loaded_columns(Arena_T arena, Loaded *ld,
                               int first_col, int end_col) {
    if (!arena || !ld || !ld->topology || !ld->domain || !ld->trusted || !ld->field ||
        ld->H < 1 || ld->W < 1 || end_col <= first_col) return -1;
    int lo = ld->c0 > first_col ? ld->c0 : first_col;
    int hi0 = ld->c0 + ld->W;
    int hi = hi0 < end_col ? hi0 : end_col;
    if (hi <= lo) return 1;
    if (lo == ld->c0 && hi == hi0) return 0;

    int offset = lo - ld->c0, new_w = hi - lo;
    size_t n = (size_t)ld->H * (size_t)new_w;
    uint8_t *topology = ARENA_ALLOC(arena, n);
    uint8_t *domain = ARENA_ALLOC(arena, n);
    uint8_t *trusted = ARENA_ALLOC(arena, n);
    float *field = ARENA_ALLOC(arena, n * 3u * sizeof *field);
    float *phase = ld->phase ? ARENA_ALLOC(arena, n * sizeof *phase) : NULL;
    for (int row = 0; row < ld->H; row++) {
        size_t src = (size_t)row * (size_t)ld->W + (size_t)offset;
        size_t dst = (size_t)row * (size_t)new_w;
        memcpy(topology + dst, ld->topology + src, (size_t)new_w);
        memcpy(domain + dst, ld->domain + src, (size_t)new_w);
        memcpy(trusted + dst, ld->trusted + src, (size_t)new_w);
        memcpy(field + dst * 3u, ld->field + src * 3u,
               (size_t)new_w * 3u * sizeof *field);
        if (phase)
            memcpy(phase + dst, ld->phase + src,
                   (size_t)new_w * sizeof *phase);
    }
    ld->topology = topology;
    ld->domain = domain;
    ld->trusted = trusted;
    ld->field = field;
    ld->phase = phase;
    ld->W = new_w;
    ld->c0 = lo;
    return 0;
}

static int loaded_crop_selftest(void) {
    Arena_T arena = Arena_new();
    Loaded ld = {0};
    ld.H = 2; ld.W = 5; ld.gr0 = 7; ld.c0 = 10;
    ld.topology = ARENA_ALLOC(arena, 10);
    ld.domain = ARENA_ALLOC(arena, 10);
    ld.trusted = ARENA_ALLOC(arena, 10);
    ld.field = ARENA_ALLOC(arena, 30 * sizeof *ld.field);
    ld.phase = ARENA_ALLOC(arena, 10 * sizeof *ld.phase);
    for (size_t k = 0; k < 10; k++) {
        ld.topology[k] = (uint8_t)(200u-k);
        ld.domain[k] = (uint8_t)(k + 1u);
        ld.trusted[k] = (uint8_t)(100u + k);
        ld.phase[k] = (float)(200u + k);
        for (int ch = 0; ch < 3; ch++)
            ld.field[k * 3u + (size_t)ch] =
                (float)(1000u + 10u * k + (size_t)ch);
    }
    int fail = crop_loaded_columns(arena, &ld, 11, 14) != 0 ||
               ld.c0 != 11 || ld.W != 3 || ld.H != 2;
    const size_t expected[6] = {1, 2, 3, 6, 7, 8};
    for (size_t k = 0; k < 6 && !fail; k++) {
        size_t src = expected[k];
        fail |= ld.topology[k] != (uint8_t)(200u-src);
        fail |= ld.domain[k] != (uint8_t)(src + 1u);
        fail |= ld.trusted[k] != (uint8_t)(100u + src);
        fail |= ld.phase[k] != (float)(200u + src);
        for (int ch = 0; ch < 3; ch++)
            fail |= ld.field[k * 3u + (size_t)ch] !=
                    (float)(1000u + 10u * src + (size_t)ch);
    }
    fail |= crop_loaded_columns(arena, &ld, 20, 22) != 1;
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] loaded UV-column crop: %s\n",
            fail ? "FAIL" : "PASS");
    return fail;
}

/* ---- confidence-preserving regular-grid ribbon coarsening ----------------
 *
 * The fitted StrokeStrip field is intentionally built at its native lattice
 * first: all upstream observations therefore participate in the cylindrical
 * continuation.  The expensive turn-order, collision, snap, and ARAP stages
 * can then operate on this substantially smaller regular grid.  A coarse
 * vertex is not a point sample.  It is the value at the coarse UV node of a
 * confidence-weighted affine fit over the fine vertices in that node's nearest
 * (Voronoi) block.  Direct trusted observations receive four times the weight
 * of continuation samples.  The coarse node remains hard-trusted only when at
 * least half of its contributing fine samples were hard-trusted; otherwise it
 * is movable, but its initializer still incorporates every fine sample.
 */
typedef struct RibbonCoarsenStats {
    int fine_h, fine_w, coarse_h, coarse_w;
    size_t fine_vertices, coarse_vertices, fine_faces, coarse_faces;
    size_t direct_samples, coarse_with_direct, coarse_trusted;
} RibbonCoarsenStats;

static int coarse_axis_count(int n, int stride) {
    if (n < 1 || stride < 1) return 0;
    if (n == 1) return 1;
    return (n - 1 + stride - 1) / stride + 1;
}

static void coarse_axis_indices(int n, int stride, int count, int *index) {
    for (int k = 0; k < count; k++) {
        long long q = (long long)k * (long long)stride;
        index[k] = q < n - 1 ? (int)q : n - 1;
    }
}

static int coarse_block_lo(const int *index, int k) {
    return k == 0 ? 0 : (index[k - 1] + index[k]) / 2 + 1;
}

static int coarse_block_hi(const int *index, int count, int k, int fine_n) {
    return k + 1 == count ? fine_n
                          : (index[k] + index[k + 1]) / 2 + 1;
}

static double coarse_affine_value(
        double target_u, double target_v, double mean_u, double mean_v,
        double mean_value, double suu, double suv, double svv,
        double su_value, double sv_value) {
    double det = suu * svv - suv * suv;
    double scale = (suu + svv) * (suu + svv);
    if (!(det > 1e-12 * fmax(1.0, scale))) return mean_value;
    double grad_u = (svv * su_value - suv * sv_value) / det;
    double grad_v = (suu * sv_value - suv * su_value) / det;
    return mean_value + grad_u * (target_u - mean_u)
                      + grad_v * (target_v - mean_v);
}

static int coarsen_rect_ribbon(
        Arena_T arena, int u_stride, int v_stride,
        float **verts_io, size_t *nv_io, int32_t **faces_io, size_t *nf_io,
        float **uv_io, uint8_t **filled_io, float **phase_io,
        int *H_io, int *W_io, RibbonCoarsenStats *out) {
    RibbonCoarsenStats s; memset(&s, 0, sizeof s);
    if (out) *out = s;
    if (!arena || !verts_io || !*verts_io || !nv_io || !faces_io ||
        !nf_io || !uv_io || !*uv_io || !filled_io || !*filled_io ||
        !H_io || !W_io || *H_io < 2 || *W_io < 2 ||
        u_stride < 1 || v_stride < 1 ||
        *nv_io != (size_t)*H_io * (size_t)*W_io) return -1;
    if (u_stride == 1 && v_stride == 1) return 0;

    const int FH = *H_io, FW = *W_io;
    const int CH = coarse_axis_count(FH, v_stride);
    const int CW = coarse_axis_count(FW, u_stride);
    if (CH < 2 || CW < 2) return -1;
    const size_t CN = (size_t)CH * (size_t)CW;
    const size_t CF = 2u * (size_t)(CH - 1) * (size_t)(CW - 1);
    int *rows = ARENA_ALLOC(arena, (size_t)CH * sizeof *rows);
    int *cols = ARENA_ALLOC(arena, (size_t)CW * sizeof *cols);
    float *cv = ARENA_ALLOC(arena, CN * 3u * sizeof *cv);
    float *cuv = ARENA_ALLOC(arena, CN * 2u * sizeof *cuv);
    uint8_t *cfilled = ARENA_ALLOC(arena, CN);
    float *cphase = phase_io && *phase_io
                  ? ARENA_ALLOC(arena, CN * sizeof *cphase) : NULL;
    int32_t *cf = ARENA_ALLOC(arena, CF * 3u * sizeof *cf);
    coarse_axis_indices(FH, v_stride, CH, rows);
    coarse_axis_indices(FW, u_stride, CW, cols);

    const float *fv = *verts_io, *fuv = *uv_io, *fphase = cphase ? *phase_io : NULL;
    const uint8_t *ffilled = *filled_io;
    for (int cr = 0; cr < CH; cr++) for (int cc = 0; cc < CW; cc++) {
        const size_t dst = (size_t)cr * (size_t)CW + (size_t)cc;
        const size_t target = (size_t)rows[cr] * (size_t)FW + (size_t)cols[cc];
        const double tu = fuv[target * 2u], tv = fuv[target * 2u + 1u];
        const int r0 = coarse_block_lo(rows, cr);
        const int r1 = coarse_block_hi(rows, CH, cr, FH);
        const int c0 = coarse_block_lo(cols, cc);
        const int c1 = coarse_block_hi(cols, CW, cc, FW);
        double sw = 0.0, su = 0.0, sv = 0.0, sp[3] = {0.0, 0.0, 0.0};
        double psw = 0.0, psu = 0.0, psv = 0.0, sph = 0.0;
        size_t samples = 0, direct = 0, phase_samples = 0;
        for (int r = r0; r < r1; r++) for (int c = c0; c < c1; c++) {
            size_t src = (size_t)r * (size_t)FW + (size_t)c;
            int is_direct = ffilled[src] == 0;
            double w = is_direct ? 4.0 : 1.0;
            double u = fuv[src * 2u], v = fuv[src * 2u + 1u];
            sw += w; su += w * u; sv += w * v;
            for (int ch = 0; ch < 3; ch++) sp[ch] += w * fv[src * 3u + (size_t)ch];
            samples++; direct += (size_t)is_direct;
            if (is_direct && fphase && isfinite((double)fphase[src])) {
                psw += 1.0; psu += u; psv += v; sph += fphase[src];
                phase_samples++;
            }
        }
        if (!(sw > 0.0) || samples == 0) return -1;
        double mu = su / sw, mv = sv / sw, mp[3];
        for (int ch = 0; ch < 3; ch++) mp[ch] = sp[ch] / sw;
        double suu = 0.0, suv = 0.0, svv = 0.0;
        double sup[3] = {0.0, 0.0, 0.0}, svp[3] = {0.0, 0.0, 0.0};
        for (int r = r0; r < r1; r++) for (int c = c0; c < c1; c++) {
            size_t src = (size_t)r * (size_t)FW + (size_t)c;
            double w = ffilled[src] == 0 ? 4.0 : 1.0;
            double du = (double)fuv[src * 2u] - mu;
            double dv = (double)fuv[src * 2u + 1u] - mv;
            suu += w * du * du; suv += w * du * dv; svv += w * dv * dv;
            for (int ch = 0; ch < 3; ch++) {
                double dp = (double)fv[src * 3u + (size_t)ch] - mp[ch];
                sup[ch] += w * du * dp; svp[ch] += w * dv * dp;
            }
        }
        cuv[dst * 2u] = (float)tu; cuv[dst * 2u + 1u] = (float)tv;
        for (int ch = 0; ch < 3; ch++)
            cv[dst * 3u + (size_t)ch] = (float)coarse_affine_value(
                tu, tv, mu, mv, mp[ch], suu, suv, svv, sup[ch], svp[ch]);
        cfilled[dst] = (uint8_t)(direct > 0 && 2u * direct >= samples ? 0 : 1);
        s.direct_samples += direct;
        s.coarse_with_direct += (size_t)(direct > 0);
        s.coarse_trusted += (size_t)(cfilled[dst] == 0);

        if (cphase) {
            cphase[dst] = NAN;
            if (phase_samples > 0) {
                double pmu = psu / psw, pmv = psv / psw, pm = sph / psw;
                double psuu = 0.0, psuv = 0.0, psvv = 0.0;
                double psuph = 0.0, psvph = 0.0;
                for (int r = r0; r < r1; r++) for (int c = c0; c < c1; c++) {
                    size_t src = (size_t)r * (size_t)FW + (size_t)c;
                    if (ffilled[src] || !isfinite((double)fphase[src])) continue;
                    double du = (double)fuv[src * 2u] - pmu;
                    double dv = (double)fuv[src * 2u + 1u] - pmv;
                    double dp = (double)fphase[src] - pm;
                    psuu += du * du; psuv += du * dv; psvv += dv * dv;
                    psuph += du * dp; psvph += dv * dp;
                }
                cphase[dst] = (float)coarse_affine_value(
                    tu, tv, pmu, pmv, pm, psuu, psuv, psvv, psuph, psvph);
            }
        }
    }

    size_t fn = 0;
    for (int r = 0; r + 1 < CH; r++) for (int c = 0; c + 1 < CW; c++) {
        int32_t a = (int32_t)((size_t)r * (size_t)CW + (size_t)c);
        int32_t b = a + 1, d = a + (int32_t)CW + 1, q = a + (int32_t)CW;
        cf[fn * 3u] = a; cf[fn * 3u + 1u] = b; cf[fn * 3u + 2u] = d; fn++;
        cf[fn * 3u] = a; cf[fn * 3u + 1u] = d; cf[fn * 3u + 2u] = q; fn++;
    }
    s.fine_h = FH; s.fine_w = FW; s.coarse_h = CH; s.coarse_w = CW;
    s.fine_vertices = *nv_io; s.coarse_vertices = CN;
    s.fine_faces = *nf_io; s.coarse_faces = CF;
    *verts_io = cv; *nv_io = CN; *faces_io = cf; *nf_io = CF;
    *uv_io = cuv; *filled_io = cfilled;
    if (phase_io) *phase_io = cphase;
    *H_io = CH; *W_io = CW;
    if (out) *out = s;
    return 0;
}

/* Compact/ragged coarsening is vertex clustering in the exported UV lattice,
 * not reconstruction of a coarse rectangle.  Every emitted coarse triangle
 * must be the non-degenerate image of at least one real fine triangle.  Thus a
 * missing alpha region, a material split, or a narrow exterior notch cannot be
 * bridged merely because its opposite rims occupy neighbouring coarse bins. */
typedef struct MaskCoarseKey {
    long long row, col;
    size_t fine;
} MaskCoarseKey;

typedef struct MaskCoarseFace {
    int32_t a, b, c;
    int32_t lo, mid, hi;
    size_t source_face;
} MaskCoarseFace;

static long long floor_div_ll(long long value, int divisor) {
    if (value >= 0) return value / divisor;
    return -(((-value) + divisor - 1) / divisor);
}

static int mask_coarse_key_cmp(const void *aa, const void *bb) {
    const MaskCoarseKey *a = (const MaskCoarseKey *)aa;
    const MaskCoarseKey *b = (const MaskCoarseKey *)bb;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->col != b->col) return a->col < b->col ? -1 : 1;
    return a->fine < b->fine ? -1 : a->fine > b->fine ? 1 : 0;
}

static int mask_coarse_face_cmp(const void *aa, const void *bb) {
    const MaskCoarseFace *a = (const MaskCoarseFace *)aa;
    const MaskCoarseFace *b = (const MaskCoarseFace *)bb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->mid != b->mid) return a->mid < b->mid ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return a->source_face < b->source_face ? -1 :
           a->source_face > b->source_face ? 1 : 0;
}

static void sort_face_key3(int32_t a, int32_t b, int32_t c,
                           int32_t *lo, int32_t *mid, int32_t *hi) {
    if (a > b) { int32_t t = a; a = b; b = t; }
    if (b > c) { int32_t t = b; b = c; c = t; }
    if (a > b) { int32_t t = a; a = b; b = t; }
    *lo = a; *mid = b; *hi = c;
}

static double uv_area2(const float *uv, int32_t a, int32_t b, int32_t c) {
    double ax = uv[(size_t)a * 2u], ay = uv[(size_t)a * 2u + 1u];
    double bx = uv[(size_t)b * 2u], by = uv[(size_t)b * 2u + 1u];
    double cx = uv[(size_t)c * 2u], cy = uv[(size_t)c * 2u + 1u];
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

static int coarsen_masked_ribbon(
        Arena_T arena, int u_stride, int v_stride, double fine_grid_du,
        float **verts_io, size_t *nv_io, int32_t **faces_io, size_t *nf_io,
        float **uv_io, uint8_t **filled_io, float **phase_io,
        int *H_io, int *W_io, RibbonCoarsenStats *out) {
    RibbonCoarsenStats s; memset(&s, 0, sizeof s);
    if (out) *out = s;
    if (!arena || !verts_io || !*verts_io || !nv_io || !*nv_io ||
        !faces_io || !*faces_io || !nf_io || !*nf_io || !uv_io || !*uv_io ||
        !filled_io || !*filled_io || !H_io || !W_io || *H_io < 1 ||
        *W_io < 1 || u_stride < 1 || v_stride < 1 ||
        !(fine_grid_du > 0.0) || !isfinite(fine_grid_du)) return -1;

    const size_t FN = *nv_io, FF = *nf_io;
    const float *fv = *verts_io, *fuv = *uv_io;
    const uint8_t *ffilled = *filled_io;
    const float *fphase = phase_io && *phase_io ? *phase_io : NULL;
    MaskCoarseKey *key = ARENA_ALLOC(arena, FN * sizeof *key);
    size_t *fine_to_coarse = ARENA_ALLOC(arena, FN * sizeof *fine_to_coarse);
    float *cv = ARENA_ALLOC(arena, FN * 3u * sizeof *cv);
    float *cuv = ARENA_ALLOC(arena, FN * 2u * sizeof *cuv);
    uint8_t *cfilled = ARENA_ALLOC(arena, FN);
    float *cphase = fphase ? ARENA_ALLOC(arena, FN * sizeof *cphase) : NULL;
    MaskCoarseFace *candidate = ARENA_ALLOC(arena, FF * sizeof *candidate);
    int32_t *cf = ARENA_ALLOC(arena, FF * 3u * sizeof *cf);
    uint8_t *used = ARENA_ALLOC(arena, FN);
    size_t *compact = ARENA_ALLOC(arena, FN * sizeof *compact);
    if (!key || !fine_to_coarse || !cv || !cuv || !cfilled ||
        (fphase && !cphase) || !candidate || !cf || !used || !compact)
        return -1;

    for (size_t v = 0; v < FN; v++) {
        double u = fuv[v * 2u], vv = fuv[v * 2u + 1u];
        if (!isfinite(u) || !isfinite(vv)) return -1;
        long long fine_col = llround(u / fine_grid_du);
        long long fine_row = llround(vv);
        if (fabs(u - (double)fine_col * fine_grid_du) > 1e-3 ||
            fabs(vv - (double)fine_row) > 1e-3) return -1;
        key[v].row = floor_div_ll(fine_row, v_stride);
        key[v].col = floor_div_ll(fine_col, u_stride);
        key[v].fine = v;
    }
    qsort(key, FN, sizeof *key, mask_coarse_key_cmp);

    size_t CN = 0;
    long long min_row = LLONG_MAX, max_row = LLONG_MIN;
    long long min_col = LLONG_MAX, max_col = LLONG_MIN;
    for (size_t begin = 0; begin < FN; ) {
        size_t end = begin + 1;
        while (end < FN && key[end].row == key[begin].row &&
               key[end].col == key[begin].col) end++;
        const size_t dst = CN++;
        if (key[begin].row < min_row) min_row = key[begin].row;
        if (key[begin].row > max_row) max_row = key[begin].row;
        if (key[begin].col < min_col) min_col = key[begin].col;
        if (key[begin].col > max_col) max_col = key[begin].col;

        /* Upstream observations are high-weight, but not immutable.  This is
         * the coarsening analogue of the 64x source attachment in CT fitting. */
        double sw = 0.0, su = 0.0, sv = 0.0, sp[3] = {0.0, 0.0, 0.0};
        double psw = 0.0, psu = 0.0, psv = 0.0, sph = 0.0;
        size_t direct = 0, phase_samples = 0;
        for (size_t q = begin; q < end; q++) {
            size_t src = key[q].fine;
            int is_direct = ffilled[src] == 0;
            double w = is_direct ? 64.0 : 1.0;
            double u = fuv[src * 2u], vv = fuv[src * 2u + 1u];
            sw += w; su += w * u; sv += w * vv;
            for (int ch = 0; ch < 3; ch++)
                sp[ch] += w * fv[src * 3u + (size_t)ch];
            direct += (size_t)is_direct;
            if (is_direct && fphase && isfinite((double)fphase[src])) {
                psw += 1.0; psu += u; psv += vv; sph += fphase[src];
                phase_samples++;
            }
            fine_to_coarse[src] = dst;
        }
        if (!(sw > 0.0)) return -1;
        double mu = su / sw, mv = sv / sw, mp[3];
        for (int ch = 0; ch < 3; ch++) mp[ch] = sp[ch] / sw;

        /* Evaluate the affine fit at a real fine-lattice sample nearest the
         * weighted centroid, so UV remains an exact lattice coordinate. */
        size_t target = key[begin].fine;
        double best = DBL_MAX;
        for (size_t q = begin; q < end; q++) {
            size_t src = key[q].fine;
            double du = ((double)fuv[src * 2u] - mu) / fine_grid_du;
            double dv = (double)fuv[src * 2u + 1u] - mv;
            double d2 = du * du + dv * dv;
            if (d2 < best || (d2 == best && ffilled[src] < ffilled[target])) {
                best = d2; target = src;
            }
        }
        double tu = fuv[target * 2u], tv = fuv[target * 2u + 1u];
        double suu = 0.0, suv = 0.0, svv = 0.0;
        double sup[3] = {0.0, 0.0, 0.0}, svp[3] = {0.0, 0.0, 0.0};
        for (size_t q = begin; q < end; q++) {
            size_t src = key[q].fine;
            double w = ffilled[src] == 0 ? 64.0 : 1.0;
            double du = (double)fuv[src * 2u] - mu;
            double dv = (double)fuv[src * 2u + 1u] - mv;
            suu += w * du * du; suv += w * du * dv; svv += w * dv * dv;
            for (int ch = 0; ch < 3; ch++) {
                double dp = (double)fv[src * 3u + (size_t)ch] - mp[ch];
                sup[ch] += w * du * dp; svp[ch] += w * dv * dp;
            }
        }
        cuv[dst * 2u] = (float)tu; cuv[dst * 2u + 1u] = (float)tv;
        for (int ch = 0; ch < 3; ch++)
            cv[dst * 3u + (size_t)ch] = (float)coarse_affine_value(
                tu, tv, mu, mv, mp[ch], suu, suv, svv, sup[ch], svp[ch]);
        size_t samples = end - begin;
        cfilled[dst] = (uint8_t)(direct > 0 && 2u * direct >= samples ? 0 : 1);
        s.direct_samples += direct;
        s.coarse_with_direct += (size_t)(direct > 0);
        s.coarse_trusted += (size_t)(cfilled[dst] == 0);

        if (cphase) {
            cphase[dst] = NAN;
            if (phase_samples > 0) {
                double pmu = psu / psw, pmv = psv / psw, pm = sph / psw;
                double psuu = 0.0, psuv = 0.0, psvv = 0.0;
                double psuph = 0.0, psvph = 0.0;
                for (size_t q = begin; q < end; q++) {
                    size_t src = key[q].fine;
                    if (ffilled[src] || !isfinite((double)fphase[src])) continue;
                    double du = (double)fuv[src * 2u] - pmu;
                    double dv = (double)fuv[src * 2u + 1u] - pmv;
                    double dp = (double)fphase[src] - pm;
                    psuu += du * du; psuv += du * dv; psvv += dv * dv;
                    psuph += du * dp; psvph += dv * dp;
                }
                cphase[dst] = (float)coarse_affine_value(
                    tu, tv, pmu, pmv, pm, psuu, psuv, psvv, psuph, psvph);
            }
        }
        begin = end;
    }

    size_t NC = 0;
    for (size_t f = 0; f < FF; f++) {
        int32_t fa = (*faces_io)[f * 3u], fb = (*faces_io)[f * 3u + 1u];
        int32_t fc0 = (*faces_io)[f * 3u + 2u];
        if (fa < 0 || fb < 0 || fc0 < 0 || (size_t)fa >= FN ||
            (size_t)fb >= FN || (size_t)fc0 >= FN) return -1;
        int32_t a = (int32_t)fine_to_coarse[(size_t)fa];
        int32_t b = (int32_t)fine_to_coarse[(size_t)fb];
        int32_t c = (int32_t)fine_to_coarse[(size_t)fc0];
        if (a == b || b == c || a == c) continue;
        double fine_area = uv_area2(fuv, fa, fb, fc0);
        double coarse_area = uv_area2(cuv, a, b, c);
        if (fabs(fine_area) <= 1e-9 || fabs(coarse_area) <= 1e-9 ||
            fine_area * coarse_area <= 0.0) continue;
        MaskCoarseFace *q = &candidate[NC++];
        q->a = a; q->b = b; q->c = c; q->source_face = f;
        sort_face_key3(a, b, c, &q->lo, &q->mid, &q->hi);
    }
    if (!NC) return -1;
    qsort(candidate, NC, sizeof *candidate, mask_coarse_face_cmp);
    size_t CF = 0;
    for (size_t q = 0; q < NC; q++) {
        if (q && candidate[q].lo == candidate[q - 1].lo &&
                 candidate[q].mid == candidate[q - 1].mid &&
                 candidate[q].hi == candidate[q - 1].hi) continue;
        cf[CF * 3u] = candidate[q].a;
        cf[CF * 3u + 1u] = candidate[q].b;
        cf[CF * 3u + 2u] = candidate[q].c;
        CF++;
    }

    memset(used, 0, CN);
    for (size_t k = 0; k < CF * 3u; k++) used[(size_t)cf[k]] = 1;
    size_t compact_n = 0;
    for (size_t v = 0; v < CN; v++) {
        if (!used[v]) { compact[v] = SIZE_MAX; continue; }
        compact[v] = compact_n;
        if (compact_n != v) {
            memcpy(&cv[compact_n * 3u], &cv[v * 3u], 3u * sizeof *cv);
            memcpy(&cuv[compact_n * 2u], &cuv[v * 2u], 2u * sizeof *cuv);
            cfilled[compact_n] = cfilled[v];
            if (cphase) cphase[compact_n] = cphase[v];
        }
        compact_n++;
    }
    for (size_t k = 0; k < CF * 3u; k++) {
        size_t v = compact[(size_t)cf[k]];
        if (v == SIZE_MAX || v > INT32_MAX) return -1;
        cf[k] = (int32_t)v;
    }

    long long ch = max_row - min_row + 1, cw = max_col - min_col + 1;
    if (ch < 1 || cw < 1 || ch > INT_MAX || cw > INT_MAX) return -1;
    s.fine_h = *H_io; s.fine_w = *W_io;
    s.coarse_h = (int)ch; s.coarse_w = (int)cw;
    s.fine_vertices = FN; s.coarse_vertices = compact_n;
    s.fine_faces = FF; s.coarse_faces = CF;
    *verts_io = cv; *nv_io = compact_n; *faces_io = cf; *nf_io = CF;
    *uv_io = cuv; *filled_io = cfilled;
    if (phase_io) *phase_io = cphase;
    *H_io = (int)ch; *W_io = (int)cw;
    if (out) *out = s;
    return 0;
}

static int coarsen_ribbon(
        Arena_T arena, int u_stride, int v_stride, double fine_grid_du,
        float **verts_io, size_t *nv_io, int32_t **faces_io, size_t *nf_io,
        float **uv_io, uint8_t **filled_io, float **phase_io,
        int *H_io, int *W_io, RibbonCoarsenStats *out) {
    if (*nv_io == (size_t)*H_io * (size_t)*W_io)
        return coarsen_rect_ribbon(arena, u_stride, v_stride, verts_io, nv_io,
            faces_io, nf_io, uv_io, filled_io, phase_io, H_io, W_io, out);
    return coarsen_masked_ribbon(arena, u_stride, v_stride, fine_grid_du,
        verts_io, nv_io, faces_io, nf_io, uv_io, filled_io, phase_io,
        H_io, W_io, out);
}

static int ribbon_coarsen_selftest(void) {
    enum { FH = 9, FW = 10 };
    float verts[FH * FW * 3], uv[FH * FW * 2], phase[FH * FW];
    uint8_t filled[FH * FW];
    int32_t dummy_faces[1] = {0};
    for (int r = 0; r < FH; r++) for (int c = 0; c < FW; c++) {
        size_t v = (size_t)r * FW + (size_t)c;
        double u = 2.0 * c, vv = (double)r;
        uv[v * 2] = (float)u; uv[v * 2 + 1] = (float)vv;
        verts[v * 3] = (float)(1.0 + 2.0 * u + 3.0 * vv);
        verts[v * 3 + 1] = (float)(-4.0 + 0.5 * u - vv);
        verts[v * 3 + 2] = (float)(7.0 - u + 0.25 * vv);
        filled[v] = (uint8_t)(r >= 7 && c >= 8);
        phase[v] = filled[v] ? NAN : (float)(0.1 * u + 0.2 * vv);
    }
    Arena_T arena = Arena_new();
    float *v = verts, *u = uv, *p = phase;
    uint8_t *f = filled; int32_t *face = dummy_faces;
    size_t nv = FH * FW, nf = 0; int H = FH, W = FW;
    RibbonCoarsenStats s;
    int fail = coarsen_rect_ribbon(arena, 3, 4, &v, &nv, &face, &nf,
                                   &u, &f, &p, &H, &W, &s) != 0;
    fail |= H != 3 || W != 4 || nv != 12 || nf != 12;
    fail |= fabs((double)u[0]) > 1e-6 || fabs((double)u[1]) > 1e-6;
    fail |= fabs((double)u[(nv - 1) * 2] - 18.0) > 1e-6;
    fail |= fabs((double)u[(nv - 1) * 2 + 1] - 8.0) > 1e-6;
    for (size_t k = 0; k < nv && !fail; k++) {
        double uu = u[k * 2], vv = u[k * 2 + 1];
        fail |= fabs((double)v[k * 3] - (1.0 + 2.0 * uu + 3.0 * vv)) > 1e-4;
        fail |= fabs((double)v[k * 3 + 1] - (-4.0 + 0.5 * uu - vv)) > 1e-4;
        fail |= fabs((double)v[k * 3 + 2] - (7.0 - uu + 0.25 * vv)) > 1e-4;
    }
    fail |= f[nv - 1] != 1 || isfinite((double)p[nv - 1]);
    fail |= s.coarse_vertices != 12 || s.coarse_faces != 12 ||
            s.coarse_with_direct != 11 || s.coarse_trusted != 11;
    for (size_t k = 0; k < nf * 3 && !fail; k++)
        fail |= face[k] < 0 || (size_t)face[k] >= nv;
    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] confidence-preserving ribbon coarsening: %s\n",
            fail ? "FAIL" : "PASS");
    return fail;
}

static int ribbon_masked_coarsen_selftest(void) {
    enum { RH = 7, RW = 10, RN = RH * RW, RFMAX = 2 * (RH - 1) * (RW - 1) };
    float verts[RN * 3], uv[RN * 2], phase[RN];
    uint8_t filled[RN];
    int32_t faces[RFMAX * 3], map[RN];
    for (int i = 0; i < RN; i++) map[i] = -1;
    size_t nv = 0, nf = 0;
    for (int r = 0; r < RH; r++) for (int c = 0; c < RW; c++) {
        /* Two sheets separated by a two-column alpha slit. */
        if (c == 4 || c == 5) continue;
        size_t v = nv++;
        map[r * RW + c] = (int32_t)v;
        double u = 2.0 * c, vv = (double)r;
        uv[v * 2u] = (float)u; uv[v * 2u + 1u] = (float)vv;
        verts[v * 3u] = (float)(3.0 + vv);
        verts[v * 3u + 1u] = (float)(10.0 + 0.25 * u);
        verts[v * 3u + 2u] = (float)(-2.0 + u - 0.5 * vv);
        filled[v] = (uint8_t)(((r + c) % 3) != 0);
        phase[v] = filled[v] ? NAN : (float)(0.2 * u + 0.1 * vv);
    }
    for (int r = 0; r + 1 < RH; r++) for (int c = 0; c + 1 < RW; c++) {
        int32_t a = map[r * RW + c], b = map[r * RW + c + 1];
        int32_t q = map[(r + 1) * RW + c], d = map[(r + 1) * RW + c + 1];
        if (a < 0 || b < 0 || q < 0 || d < 0) continue;
        faces[nf * 3u] = a; faces[nf * 3u + 1u] = b;
        faces[nf * 3u + 2u] = d; nf++;
        faces[nf * 3u] = a; faces[nf * 3u + 1u] = d;
        faces[nf * 3u + 2u] = q; nf++;
    }
    Arena_T arena = Arena_new();
    float *v = verts, *u = uv, *p = phase;
    uint8_t *f = filled; int32_t *face = faces;
    int H = RH, W = RW; RibbonCoarsenStats s;
    int fail = coarsen_ribbon(arena, 2, 2, 2.0, &v, &nv, &face, &nf,
                              &u, &f, &p, &H, &W, &s) != 0;
    fail |= H != 4 || W != 5 || nv != 16 || nf != 12;
    int parent[RN];
    for (size_t i = 0; i < nv; i++) parent[i] = (int)i;
    for (size_t k = 0; k < nf && !fail; k++) {
        int32_t a = face[k * 3u], b = face[k * 3u + 1u], c = face[k * 3u + 2u];
        fail |= a < 0 || b < 0 || c < 0 || (size_t)a >= nv ||
                (size_t)b >= nv || (size_t)c >= nv;
        fail |= uv_area2(u, a, b, c) <= 0.0;
        int crosses = 0;
        int32_t tri[3] = {a, b, c};
        for (int i = 0; i < 3; i++) for (int j = i + 1; j < 3; j++)
            crosses |= (u[(size_t)tri[i] * 2u] <= 6.0f &&
                        u[(size_t)tri[j] * 2u] >= 12.0f) ||
                       (u[(size_t)tri[j] * 2u] <= 6.0f &&
                        u[(size_t)tri[i] * 2u] >= 12.0f);
        fail |= crosses;
        for (int e = 1; e < 3; e++) {
            int x = tri[0], y = tri[e];
            while (parent[x] != x) x = parent[x];
            while (parent[y] != y) y = parent[y];
            if (x != y) parent[y] = x;
        }
    }
    int components = 0;
    for (size_t i = 0; i < nv; i++) {
        int x = (int)i;
        while (parent[x] != x) x = parent[x];
        components += x == (int)i;
    }
    fail |= components != 2 || s.coarse_vertices != 16 ||
            s.coarse_faces != 12;
    Arena_dispose(&arena);
    fprintf(stderr,
        "[selftest] mask-aware ribbon coarsening preserves two alpha-separated sheets: %s\n",
        fail ? "FAIL" : "PASS");
    return fail;
}

/* ---- growable output buffers (malloc; one-time assembly, not a hot path) -- */
typedef struct { float   *v; size_t n, cap; } FVec;
typedef struct { int32_t *v; size_t n, cap; } IVec;
typedef struct { uint8_t *v; size_t n, cap; } BVec;

typedef struct PassthroughStats {
    size_t components;
    size_t vertices;
    size_t supported_vertices;
    size_t generated_vertices;
    size_t faces;
} PassthroughStats;

typedef struct PassthroughAtlas {
    size_t component_index;
    long label;
    int height, width;
    double u_min, u_max, v_min, v_max;
} PassthroughAtlas;

typedef struct FinalFoldRepairStats {
    size_t components_audited;
    size_t components_repaired;
    size_t input_folds;
    size_t faces_removed;
    size_t output_folds;
} FinalFoldRepairStats;

/* The proximal solver's exact gate is deliberately monotone: it refuses to
 * create a fold but cannot erase a fold already present in the initializer.
 * Once geometry is final, remove only a minimum-area cover of shared-boundary
 * fold pairs.  Ordinary inter-layer overlap and stab pairs are excluded by the
 * kind mask, and a second exact audit must certify zero residual folds. */
static int repair_final_fold_flaps(
        float *verts, size_t nv, int32_t *faces, size_t *pnf,
        int known_count_valid, size_t known_count,
        const char *domain, long identity, FinalFoldRepairStats *totals) {
    IntersectionCleanupParams p;
    IntersectionCleanupStats before = {0}, after = {0};
    size_t nf_before;
    if (!verts || !faces || !pnf || !totals) return -1;
    totals->components_audited++;
    if (*pnf < 2 || (known_count_valid && known_count == 0)) return 0;
    nf_before = *pnf;
    IntersectionCleanup_default_params(&p);
    p.gap_max = 1.0;
    p.parallel_angle_deg = 20.0;
    p.max_conflicts = *pnf > 1000000u ? *pnf : 1000000u;
    p.include_hinges = 0;
    p.hit_kind_mask = INTERSECTION_HIT_MASK_FOLD;
    if (IntersectionCleanup_process(
            verts, nv, faces, pnf, NULL, &p, &before) != 0) {
        fprintf(stderr,
                "ERROR: final fold-only repair exceeded its bounded deletion "
                "budget for %s %ld (%zu input fold pairs)\n",
                domain, identity, before.fold_pairs);
        return -1;
    }
    if (IntersectionCleanup_audit(
            verts, nv, faces, *pnf, NULL, &p, NULL, &after) != 0 ||
        after.fold_pairs != 0) {
        fprintf(stderr,
                "ERROR: final fold-only repair did not certify %s %ld "
                "(%zu residual fold pairs)\n",
                domain, identity, after.fold_pairs);
        return -1;
    }
    totals->input_folds += before.fold_pairs;
    totals->faces_removed += before.faces_deleted;
    totals->output_folds += after.fold_pairs;
    totals->components_repaired += before.faces_deleted != 0;
    if (before.fold_pairs != 0 || nf_before != *pnf) {
        fprintf(stderr,
                "[fold-repair] %s=%ld exact folds %zu->%zu; removed %zu "
                "small flap face(s), faces %zu->%zu\n",
                domain, identity, before.fold_pairs, after.fold_pairs,
                before.faces_deleted, nf_before, *pnf);
    }
    return 0;
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (!q) { fprintf(stderr, "solid_quad_ribbon: OOM (%zu bytes)\n", n); exit(1); }
    return q;
}
static void fv_reserve(FVec *a, size_t extra) {
    if (a->n + extra > a->cap) { a->cap = (a->n + extra) * 2; a->v = xrealloc(a->v, a->cap * sizeof *a->v); }
}
static void iv_reserve(IVec *a, size_t extra) {
    if (a->n + extra > a->cap) { a->cap = (a->n + extra) * 2; a->v = xrealloc(a->v, a->cap * sizeof *a->v); }
}
static void bv_reserve(BVec *a, size_t extra) {
    if (a->n + extra > a->cap) { a->cap = (a->n + extra) * 2; a->v = xrealloc(a->v, a->cap * sizeof *a->v); }
}

/* QRCORR1 is deliberately smaller and stricter than a JSON vertex manifest.
 * One 24-byte record accompanies every VESMESH1 vertex in the same order.  UV
 * coordinates are captured before shelf packing, so independently packed
 * geometry and bake-proxy meshes can meet again by exact common-lattice key.
 * Namespace 0 means material lineage; namespace 1 means an unfitted aggregate
 * topology label. */
static int write_correspondence(
        const char *path, double grid_du, size_t nv,
        const FVec *canonical_uv, const IVec *name_space,
        const IVec *lineage) {
    static const unsigned char magic[8] = {
        'Q','R','C','O','R','R','1','\0'
    };
    const uint32_t version = 1u, endian = UINT32_C(0x01020304);
    const uint32_t header_bytes = 64u, record_bytes = 24u;
    const uint64_t records = (uint64_t)nv, zero = 0;
    char temporary[2048];
    FILE *f = NULL;
    int n = 0, ok = 1;
    if (!path || !canonical_uv || !name_space || !lineage ||
        !isfinite(grid_du) || grid_du <= 0.0 ||
        canonical_uv->n != nv * 2u || name_space->n != nv ||
        lineage->n != nv) {
        fprintf(stderr, "invalid correspondence payload\n");
        return -1;
    }
    n = snprintf(temporary, sizeof temporary, "%s.tmp", path);
    if (n < 0 || (size_t)n >= sizeof temporary) return -1;
    f = fopen(temporary, "wb");
    if (!f) {
        fprintf(stderr, "cannot write correspondence %s\n", temporary);
        return -1;
    }
#define CORR_WRITE(value) \
    do { if (fwrite(&(value), sizeof(value), 1, f) != 1) ok = 0; } while (0)
    if (fwrite(magic, 1, sizeof magic, f) != sizeof magic) ok = 0;
    CORR_WRITE(version); CORR_WRITE(endian); CORR_WRITE(header_bytes);
    CORR_WRITE(record_bytes); CORR_WRITE(records); CORR_WRITE(grid_du);
    CORR_WRITE(zero); CORR_WRITE(zero); CORR_WRITE(zero);
    for (size_t i = 0; i < nv && ok; i++) {
        int32_t ns = name_space->v[i], id = lineage->v[i];
        int64_t uc = (int64_t)llround(
            (double)canonical_uv->v[i * 2u] / grid_du);
        int64_t vr = (int64_t)llround(
            (double)canonical_uv->v[i * 2u + 1u]);
        double ue = fabs((double)canonical_uv->v[i * 2u] -
                         (double)uc * grid_du);
        double ve = fabs((double)canonical_uv->v[i * 2u + 1u] -
                         (double)vr);
        if (ue > 1e-4 || ve > 1e-4) {
            fprintf(stderr,
                "correspondence vertex %zu is off the common lattice "
                "(uerr=%.9g verr=%.9g)\n", i, ue, ve);
            ok = 0;
            break;
        }
        CORR_WRITE(ns); CORR_WRITE(id); CORR_WRITE(uc); CORR_WRITE(vr);
    }
#undef CORR_WRITE
    if (fclose(f) != 0) ok = 0;
    if (!ok) {
        remove(temporary);
        return -1;
    }
    remove(path);
    if (rename(temporary, path) != 0) {
        fprintf(stderr, "cannot publish correspondence %s\n", path);
        remove(temporary);
        return -1;
    }
    return 0;
}

static int make_output_dir(const char *path) {
#ifdef _WIN32
    int rc = _mkdir(path);
#else
    int rc = mkdir(path, 0777);
#endif
    return rc == 0 || errno == EEXIST ? 0 : -1;
}

/* Apply one image-space bake/refit mask before QuadStrip_build.  The mask is a
 * local HxW vertex raster emitted by audit_quad_strip_bake.py.  Removing those
 * Dirichlet samples makes the ordinary cylindrical/harmonic initializer refit
 * the broad failed-prediction patch from its closed material boundary; setting
 * `filled` only after the build would retain the wrong 3-D initializer and is
 * therefore deliberately not supported. */
static int apply_bake_refit_mask(Arena_T arena, const char *dir, long group,
                                 uint8_t *domain, uint8_t *trusted, int H, int W,
                                 size_t *out_requested, size_t *out_removed) {
    char path[1600];
    FILE *probe = NULL;
    uint8_t *mask = NULL;
    int D = 0, MH = 0, MW = 0;
    size_t requested = 0, removed = 0, boundary = 0, remaining = 0;
    *out_requested = 0;
    *out_removed = 0;
    if (!dir || !domain || !trusted || H < 2 || W < 2) return 0;
    snprintf(path, sizeof path, "%s/group_%ld.tif", dir, group);
    probe = fopen(path, "rb");
    if (!probe) return 0; /* No broad-dark candidate for this group. */
    fclose(probe);
    if (TiffIO_load(arena, path, &mask, &D, &MH, &MW) != 0 ||
        D != 1 || MH != H || MW != W) {
        fprintf(stderr, "[bake-refit] invalid mask %s: got %dx%dx%d, expected 1x%dx%d\n",
                path, D, MH, MW, H, W);
        return -1;
    }
    for (int r = 0; r < H; r++) for (int c = 0; c < W; c++) {
        size_t v = (size_t)r*(size_t)W+(size_t)c;
        if (!mask[v]) continue;
        requested++;
        if (r == 0 || c == 0 || r + 1 == H || c + 1 == W) boundary++;
    }
    if (boundary > 0) {
        fprintf(stderr, "[bake-refit] refusing %s: %zu masked vertices touch the "
                        "outer lattice boundary\n", path, boundary);
        return -1;
    }
    if (requested > (size_t)H*(size_t)W/3u) {
        fprintf(stderr, "[bake-refit] refusing %s: %zu/%zu vertices exceeds the "
                        "one-third safety cap\n", path, requested,
                (size_t)H*(size_t)W);
        return -1;
    }
    for (size_t v = 0; v < (size_t)H*(size_t)W; v++) {
        if (mask[v] && domain[v]) { domain[v] = 0; trusted[v] = 0; removed++; }
        if (domain[v]) remaining++;
    }
    if (remaining < 4) {
        fprintf(stderr, "[bake-refit] refusing %s: only %zu fitted boundary samples remain\n",
                path, remaining);
        return -1;
    }
    *out_requested = requested;
    *out_removed = removed;
    fprintf(stderr, "[bake-refit] group=%ld mask=%s requested=%zu removed-evidence=%zu "
                    "remaining=%zu\n", group, path, requested, removed, remaining);
    return 1;
}

typedef struct InitQuarantineStats {
    size_t measured_triangles, severe_triangles;
    size_t marked_tiles, dilated_tiles;
    size_t removed_initializer, removed_trusted;
} InitQuarantineStats;

typedef struct PhaseGateStats {
    size_t observations,edges,ambiguous_observations,discontinuous_edges;
    size_t removed;
    double mismatch_rms,mismatch_max;
} PhaseGateStats;

static int init_field_vertex_finite(const float *field,size_t v) {
    return isfinite((double)field[v*3]) &&
           isfinite((double)field[v*3+1]) &&
           isfinite((double)field[v*3+2]);
}

static double init_triangle_sander(const float *field,size_t a,size_t b,size_t c,
                                   double sa,double ta,double sb,double tb,
                                   double sc,double tc) {
    double A2=(sb-sa)*(tc-ta)-(sc-sa)*(tb-ta);
    if(fabs(A2)<1e-12)return 1e9;
    double aa=0.0,bb=0.0,cc=0.0;
    for(int k=0;k<3;k++) {
        double q1=(double)field[a*3+(size_t)k];
        double q2=(double)field[b*3+(size_t)k];
        double q3=(double)field[c*3+(size_t)k];
        double ds=(q1*(tb-tc)+q2*(tc-ta)+q3*(ta-tb))/A2;
        double dt=(q1*(sc-sb)+q2*(sa-sc)+q3*(sb-sa))/A2;
        aa+=ds*ds;bb+=ds*dt;cc+=dt*dt;
    }
    double disc=sqrt(fmax(0.0,(aa-cc)*(aa-cc)+4.0*bb*bb));
    double lmax=0.5*(aa+cc+disc),lmin=0.5*(aa+cc-disc);
    if(lmin<=1e-12)return 1e9;
    double fwd=sqrt(0.5*(lmax+lmin));
    double inv=sqrt(0.5*(1.0/lmax+1.0/lmin));
    return fwd>inv?fwd:inv;
}

/* Excise broad initializer patches whose emitted-triangle metric is already
 * badly folded or stretched.  Measure the dense initializer domain, including
 * accepted interpolation rather than only direct source support: StrokeStrip
 * support commonly occupies alternating z rows, so a support-only test contains
 * literally no triangles and silently certifies every false-clean marbled strip.
 * Never measure merely-finite values outside the explicit topology.  Tile density
 * distinguishes a coherent failed parameterization from a single naturally
 * compressed cell; dilation supplies a low-distortion boundary for the
 * cylindrical refill. */
static int quarantine_initializer_geometry(
        Arena_T arena,uint8_t *domain,uint8_t *trusted,const float *field,
        int H,int W,double grid_du,double threshold,int tile,double fraction,
        int dilation,InitQuarantineStats *stats) {
    InitQuarantineStats s={0};if(stats)*stats=s;
    if(!arena||!domain||!trusted||!field||H<2||W<2||grid_du<=0.0||
       threshold<=0.0||tile<1||fraction<=0.0||dilation<0)return 0;
    int TH=((H-1)+tile-1)/tile,TW=((W-1)+tile-1)/tile;
    size_t TN=(size_t)TH*(size_t)TW,HW=(size_t)H*(size_t)W;
    size_t *count=ARENA_CALLOC(arena,TN,sizeof *count);
    size_t *severe=ARENA_CALLOC(arena,TN,sizeof *severe);
    uint8_t *mark=ARENA_CALLOC(arena,TN,1);
    uint8_t *tmp=ARENA_CALLOC(arena,TN,1);
    for(int r=0;r+1<H;r++)for(int c=0;c+1<W;c++) {
        size_t a=(size_t)r*(size_t)W+(size_t)c,b=a+1;
        size_t cc=a+(size_t)W,d=cc+1;
        if(!domain[a]||!domain[b]||!domain[cc]||!domain[d]||
           !init_field_vertex_finite(field,a)||!init_field_vertex_finite(field,b)||
           !init_field_vertex_finite(field,cc)||!init_field_vertex_finite(field,d))continue;
        size_t t=(size_t)(r/tile)*(size_t)TW+(size_t)(c/tile);
        size_t tri[2][3]={{a,b,d},{a,d,cc}};
        double uv[4][2]={{(double)c*grid_du,(double)r},
                         {(double)(c+1)*grid_du,(double)r},
                         {(double)c*grid_du,(double)(r+1)},
                         {(double)(c+1)*grid_du,(double)(r+1)}};
        int corner[2][3]={{0,1,3},{0,3,2}};
        for(int q=0;q<2;q++) {
            double D=init_triangle_sander(field,tri[q][0],tri[q][1],tri[q][2],
                uv[corner[q][0]][0],uv[corner[q][0]][1],
                uv[corner[q][1]][0],uv[corner[q][1]][1],
                uv[corner[q][2]][0],uv[corner[q][2]][1]);
            count[t]++;s.measured_triangles++;
            if(D>=threshold){severe[t]++;s.severe_triangles++;}
        }
    }
    size_t min_samples=(size_t)tile*(size_t)tile/2u;if(min_samples<16)min_samples=16;
    for(size_t t=0;t<TN;t++)if(count[t]>=min_samples&&
        (double)severe[t]>=(double)count[t]*fraction){mark[t]=1;s.marked_tiles++;}
    for(int it=0;it<dilation;it++) {
        memcpy(tmp,mark,TN);
        for(int r=0;r<TH;r++)for(int c=0;c<TW;c++)if(mark[(size_t)r*TW+c])
            for(int dr=-1;dr<=1;dr++)for(int dc=-1;dc<=1;dc++) {
                int rr=r+dr,cc=c+dc;if(rr>=0&&rr<TH&&cc>=0&&cc<TW)
                    tmp[(size_t)rr*TW+cc]=1;
            }
        memcpy(mark,tmp,TN);
    }
    for(size_t t=0;t<TN;t++)s.dilated_tiles+=mark[t]!=0;
    size_t would_remove=0,domain_before=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        int tr=(r<H-1?r:H-2)/tile,tc=(c<W-1?c:W-2)/tile;
        size_t t=(size_t)tr*(size_t)TW+(size_t)tc;
        size_t v=(size_t)r*(size_t)W+(size_t)c;
        domain_before+=domain[v]!=0;
        if(mark[t]&&domain[v])would_remove++;
    }
    if(would_remove>HW/3u) {
        fprintf(stderr,"[init-quarantine] refusing %zu/%zu initializer removals "
                       "above one-third safety cap\n",would_remove,HW);
        return -1;
    }
    if(domain_before-would_remove<4u) {
        fprintf(stderr,"[init-quarantine] refusing removal: only %zu fitted "
                       "observations would remain\n",domain_before-would_remove);
        return -1;
    }
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        int tr=(r<H-1?r:H-2)/tile,tc=(c<W-1?c:W-2)/tile;
        size_t t=(size_t)tr*(size_t)TW+(size_t)tc;
        size_t v=(size_t)r*(size_t)W+(size_t)c;
        if(!mark[t]||!domain[v])continue;
        domain[v]=0;s.removed_initializer++;
        if(trusted[v]){trusted[v]=0;s.removed_trusted++;}
    }
    if(stats)*stats=s;
    fprintf(stderr,"[init-quarantine] dense-domain triangles=%zu severe=%zu (>=%.2fx); "
        "tiles=%zu->%zu; removed initializer/trusted=%zu/%zu\n",
        s.measured_triangles,s.severe_triangles,threshold,s.marked_tiles,
        s.dilated_tiles,s.removed_initializer,s.removed_trusted);
    return 0;
}

/* Verify that the upstream universal-cover phase really selects an unambiguous
 * atan2 branch for the CURRENT fitted XYZ.  A normal atan2 seam changes the
 * integer branch but leaves lifted theta continuous, so it is not rejected.
 * Ambiguous observations and endpoints of impossible >pi observed edges are
 * removed from the Dirichlet set; they become cylindrical PDE unknowns. */
static int quarantine_phase_gauge(
        Arena_T arena,uint8_t *domain,uint8_t *trusted,const float *field,
        const float *phase,int H,int W,double axis_y,double axis_x,
        PhaseGateStats *stats) {
    const double pi=3.1415926535897932384626433832795;
    const double two_pi=6.283185307179586476925286766559;
    size_t HW=(size_t)H*(size_t)W,domain_before=0;
    uint8_t *remove=ARENA_CALLOC(arena,HW,1);
    double *theta=ARENA_ALLOC(arena,HW*sizeof *theta);
    PhaseGateStats s={0};double mismatch_sum2=0.0;
    if(!phase)return -1;
    for(size_t k=0;k<HW;k++) {
        theta[k]=NAN;
        if(!domain[k])continue;
        domain_before++;
        double ref=(double)phase[k];
        double dy=(double)field[k*3+1]-axis_y;
        double dx=(double)field[k*3+2]-axis_x;
        if(!isfinite(ref)||!isfinite(dy)||!isfinite(dx))return -1;
        double raw=atan2(dx,dy);
        theta[k]=raw+nearbyint((ref-raw)/two_pi)*two_pi;
        double mismatch=fabs(theta[k]-ref);
        mismatch_sum2+=mismatch*mismatch;
        if(mismatch>s.mismatch_max)s.mismatch_max=mismatch;
        if(mismatch>0.5*pi){remove[k]=1;s.ambiguous_observations++;}
        s.observations++;
    }
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t a=(size_t)r*(size_t)W+(size_t)c;
        if(!domain[a])continue;
        size_t neighbour[2];int count=0;
        if(c+1<W)neighbour[count++]=a+1;
        if(r+1<H)neighbour[count++]=a+(size_t)W;
        for(int q=0;q<count;q++) {
            size_t b=neighbour[q];
            if(!domain[b])continue;
            s.edges++;
            if(fabs(theta[b]-theta[a])>pi) {
                remove[a]=remove[b]=1;
                s.discontinuous_edges++;
            }
        }
    }
    for(size_t k=0;k<HW;k++)s.removed+=domain[k]&&remove[k];
    s.mismatch_rms=s.observations
                  ?sqrt(mismatch_sum2/(double)s.observations):0.0;
    if(s.removed>domain_before/3u||domain_before-s.removed<4u) {
        fprintf(stderr,
                "[phase-gate] refusing unsafe removal %zu/%zu observations\n",
                s.removed,domain_before);
        return -1;
    }
    for(size_t k=0;k<HW;k++)if(domain[k]&&remove[k]) {
        domain[k]=0;trusted[k]=0;
    }
    fprintf(stderr,
            "[phase-gate] observations=%zu mismatch_rms/max=%.4f/%.4f rad "
            "ambiguous=%zu; local_edges=%zu discontinuous=%zu; removed=%zu\n",
            s.observations,s.mismatch_rms,s.mismatch_max,
            s.ambiguous_observations,s.edges,s.discontinuous_edges,s.removed);
    if(stats)*stats=s;
    return 0;
}

/* Emit one RECT component as variable-length, independently loadable ribbon
 * snips.  A snip owns a half-open run of quad cells and includes both endpoint
 * vertex columns.  Consecutive snips therefore share exactly one identical
 * vertex column, but never duplicate faces.  The manifest's global lattice
 * coordinates make that seam deterministic for a future streaming joiner;
 * MeshLab can append the OBJ layers and merge coincident seam vertices. */
/* ---- CT-occupancy continuation trim --------------------------------------
 * User rule: harmonic fill over BLACK CT is deleted; fill with a glimmer of
 * material is kept and recovered by the snap.  A quad cell is "black" when
 * all four corners are PDE fill AND no bake-style sample in it was ever
 * non-dark AND no dark sample ever earned a supported bright alternative
 * (both accumulated across every snap pass, including rejected ones --
 * deletion is the unfixable direction, so the evidence bias is toward KEEP).
 * Cells within the rim guard of any fitted-corner cell are kept for seam
 * integrity, and pepper regions below the minimum size are kept. */
typedef struct TrimStats {
    size_t black_cells, regions, trimmed_cells, faces_removed;
    size_t kept_small_region, kept_rim_guard;
} TrimStats;

static void trim_black_continuation(
        const uint8_t *filled, const uint8_t *cell_lit,
        const uint8_t *cell_cand, int H, int W,
        int32_t *faces, size_t *pnf,
        size_t min_region, int rim_guard, TrimStats *out) {
    memset(out, 0, sizeof *out);
    if (H < 2 || W < 2 || filled == NULL || cell_lit == NULL ||
        cell_cand == NULL || faces == NULL || pnf == NULL) return;
    size_t cw = (size_t)(W - 1), ch = (size_t)(H - 1), nc = cw * ch;
    uint8_t *black = (uint8_t *)calloc(nc, 1);
    uint8_t *fitted_cell = (uint8_t *)calloc(nc, 1);
    int32_t *region = (int32_t *)malloc(nc * sizeof *region);
    int32_t *queue = (int32_t *)malloc(nc * sizeof *queue);
    if (black == NULL || fitted_cell == NULL || region == NULL ||
        queue == NULL) {
        free(black); free(fitted_cell); free(region); free(queue);
        return;
    }
    for (size_t r = 0; r < ch; r++) for (size_t c = 0; c < cw; c++) {
        size_t cell = r * cw + c;
        size_t a = r * (size_t)W + c, b = a + 1;
        size_t d = a + (size_t)W, e = d + 1;
        int nfill = (filled[a] != 0) + (filled[b] != 0) +
                    (filled[d] != 0) + (filled[e] != 0);
        fitted_cell[cell] = nfill < 4;
        if (nfill == 4 && !cell_lit[cell] && !cell_cand[cell]) {
            black[cell] = 1;
            out->black_cells++;
        }
    }
    /* rim guard: keep black cells near any fitted-corner cell */
    if (rim_guard > 0) {
        for (size_t r = 0; r < ch; r++) for (size_t c = 0; c < cw; c++) {
            size_t cell = r * cw + c;
            if (!black[cell]) continue;
            int keep = 0;
            for (int dr = -rim_guard; dr <= rim_guard && !keep; dr++) {
                long rr = (long)r + dr;
                if (rr < 0 || rr >= (long)ch) continue;
                for (int dc = -rim_guard; dc <= rim_guard; dc++) {
                    long cc = (long)c + dc;
                    if (cc < 0 || cc >= (long)cw) continue;
                    if (fitted_cell[(size_t)rr * cw + (size_t)cc]) {
                        keep = 1; break;
                    }
                }
            }
            if (keep) { black[cell] = 0; out->kept_rim_guard++; }
        }
    }
    /* 4-connected regions; keep pepper regions below the size floor */
    for (size_t i = 0; i < nc; i++) region[i] = -1;
    int32_t nregion = 0;
    for (size_t seed = 0; seed < nc; seed++) {
        if (!black[seed] || region[seed] >= 0) continue;
        size_t head = 0, tail = 0, count = 0;
        queue[tail++] = (int32_t)seed; region[seed] = nregion;
        while (head < tail) {
            size_t cell = (size_t)queue[head++]; count++;
            size_t r = cell / cw, c = cell % cw;
            size_t nb[4]; int nn = 0;
            if (c > 0) nb[nn++] = cell - 1;
            if (c + 1 < cw) nb[nn++] = cell + 1;
            if (r > 0) nb[nn++] = cell - cw;
            if (r + 1 < ch) nb[nn++] = cell + cw;
            for (int k = 0; k < nn; k++)
                if (black[nb[k]] && region[nb[k]] < 0) {
                    region[nb[k]] = nregion; queue[tail++] = (int32_t)nb[k];
                }
        }
        out->regions++;
        if (count < min_region) {
            for (size_t i = 0; i < nc; i++)
                if (region[i] == nregion) black[i] = 0;
            out->kept_small_region += count;
        } else {
            out->trimmed_cells += count;
        }
        nregion++;
    }
    /* remove both faces of each trimmed cell (face -> cell via f/2) */
    size_t nf = *pnf, w = 0;
    for (size_t f = 0; f < nf; f++) {
        size_t cell = f / 2u;
        if (cell < nc && black[cell]) { out->faces_removed++; continue; }
        faces[w*3+0] = faces[f*3+0];
        faces[w*3+1] = faces[f*3+1];
        faces[w*3+2] = faces[f*3+2];
        w++;
    }
    *pnf = w;
    free(black); free(fitted_cell); free(region); free(queue);
}

static int continuation_trim_selftest(void) {
    int fail = 0;
    /* 7x9 vertex rect (6x8 cells).  Vertex columns 0-1 are fitted, so cells
     * c<=1 have fitted corners; everything else is filled.  One lit cell and
     * one candidate cell sit inside the black area and must survive; the rim
     * guard (2) protects black cells at c<=3. */
    enum { TH = 7, TW = 9, TCH = TH - 1, TCW = TW - 1 };
    uint8_t filled[TH * TW];
    for (int r = 0; r < TH; r++) for (int c = 0; c < TW; c++)
        filled[r * TW + c] = c >= 2;   /* columns 0-1 fitted */
    uint8_t lit[TCH * TCW]; memset(lit, 0, sizeof lit);
    uint8_t cand[TCH * TCW]; memset(cand, 0, sizeof cand);
    lit[1 * TCW + 6] = 1;              /* glimmer keeps this cell */
    cand[2 * TCW + 6] = 1;             /* supported alternative keeps this */
    /* faces: full rect triangulation, 2 per cell, face f -> cell f/2 */
    int32_t faces[TCH * TCW * 2 * 3];
    size_t nf = 0;
    for (int r = 0; r < TCH; r++) for (int c = 0; c < TCW; c++) {
        int32_t a = r * TW + c, b = a + 1, d = a + TW, e = d + 1;
        faces[nf*3+0]=a; faces[nf*3+1]=b; faces[nf*3+2]=e; nf++;
        faces[nf*3+0]=a; faces[nf*3+1]=e; faces[nf*3+2]=d; nf++;
    }
    TrimStats ts = {0};
    size_t nf_out = nf;
    trim_black_continuation(filled, lit, cand, TH, TW, faces, &nf_out,
                            4, 2, &ts);
    size_t expected_black = 0, expected_trim = 0;
    for (int r = 0; r < TCH; r++) for (int c = 0; c < TCW; c++) {
        int all4 = c >= 2;             /* cell filled when its c>=2 (corners c..c+1) */
        int is_lit = lit[r * TCW + c] || cand[r * TCW + c];
        if (all4 && !is_lit) expected_black++;
        /* fitted-corner cells are c<=1; rim guard 2 keeps black c<=3 */
        if (all4 && !is_lit && c >= 4) expected_trim++;
    }
    int ok1 = ts.black_cells == expected_black;
    int ok2 = ts.trimmed_cells == expected_trim;
    int ok3 = nf_out == nf - 2 * ts.trimmed_cells;
    int ok4 = ts.kept_rim_guard > 0;
    fprintf(stderr,
        "  %s: continuation trim classifies black fill (black=%zu/%zu)\n",
        ok1 ? "ok" : "FAIL", ts.black_cells, expected_black);
    fprintf(stderr,
        "  %s: rim guard + lit/candidate cells survive; %zu cells trimmed "
        "(expect %zu), faces %zu -> %zu\n",
        (ok2 && ok3 && ok4) ? "ok" : "FAIL",
        ts.trimmed_cells, expected_trim, nf, nf_out);
    fail += !ok1; fail += !(ok2 && ok3 && ok4);
    return fail;
}

static int emit_rect_snips(const char *dir, FILE *manifest, int *manifest_first,
                           size_t *serial_io, int snip_cell_cols,
                           size_t component_index, long group, long representative_track,
                           int H, int W, int gr0, int c0, double grid_du,
                           int row_stride, int col_stride,
                           const float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const float *uv, const uint8_t *filled,
                           size_t *out_nv, size_t *out_nf) {
    if (!dir || !manifest || H < 1 || W < 1 || snip_cell_cols < 1 ||
        row_stride < 1 || col_stride < 1) return -1;
    int cells = W > 1 ? W - 1 : 0;
    int nparts = cells > 0 ? (cells + snip_cell_cols - 1) / snip_cell_cols : 1;
    int structured = nv == (size_t)H * (size_t)W;
    size_t serial0 = *serial_io;
    size_t *part_serial = (size_t *)malloc((size_t)nparts*sizeof *part_serial);
    if (!part_serial) return -1;
    for (int part=0;part<nparts;part++) part_serial[part]=SIZE_MAX;

    /* Sparse bounded topology can leave an entire longitudinal partition empty.
     * Assign serials only to partitions containing emitted triangles. */
    size_t active_parts=0;
    for (int part=0;part<nparts;part++) {
        int cs=cells>0?part*snip_cell_cols:0;
        int ce=cells>0?(part+1)*snip_cell_cols:0;
        if(ce>cells)ce=cells;
        int gc0=c0+cs*col_stride,gc1=c0+ce*col_stride,active=0;
        for(size_t f=0;f<nf&&!active;f++) {
            active=1;
            for(int q=0;q<3;q++) {
                size_t v=(size_t)faces[f*3+(size_t)q];
                if (structured) {
                    int lc=(int)(v%(size_t)W);
                    if(lc<cs||lc>ce){active=0;break;}
                } else {
                    int gc=(int)lround((double)uv[v*2]/grid_du);
                    if(gc<gc0||gc>gc1){active=0;break;}
                }
            }
        }
        if(active)part_serial[part]=serial0+active_parts++;
    }
    for (int part = 0; part < nparts; part++) {
        if(part_serial[part]==SIZE_MAX)continue;
        int cs = cells > 0 ? part * snip_cell_cols : 0;
        int ce = cells > 0 ? (part + 1) * snip_cell_cols : 0;
        if (ce > cells) ce = cells;
        int SW = ce - cs + 1;
        int gc0=c0+cs*col_stride, gc1=c0+ce*col_stride;
        int32_t *map=(int32_t *)malloc(nv*sizeof *map);
        if(!map){free(part_serial);return -1;}
        size_t snv=0;
        for(size_t v=0;v<nv;v++) {
            int inside;
            if (structured) {
                int lc=(int)(v%(size_t)W);
                inside=lc>=cs&&lc<=ce;
            } else {
                int gc=(int)lround((double)uv[v*2]/grid_du);
                inside=gc>=gc0&&gc<=gc1;
            }
            map[v]=inside?(int32_t)snv++:-1;
        }
        size_t snf=0;
        for(size_t f=0;f<nf;f++)
            if(map[(size_t)faces[f*3]]>=0&&map[(size_t)faces[f*3+1]]>=0&&
               map[(size_t)faces[f*3+2]]>=0)snf++;
        float *sv = (float *)malloc(snv * 3 * sizeof *sv);
        float *su = (float *)malloc(snv * 2 * sizeof *su);
        int32_t *sf = (int32_t *)malloc((snf ? snf : 1) * 3 * sizeof *sf);
        if (!sv || !su || !sf) {
            free(map);free(sv);free(su);free(sf);free(part_serial);return -1;
        }
        size_t nfill = 0;
        double u_min=DBL_MAX,u_max=-DBL_MAX,v_min=DBL_MAX,v_max=-DBL_MAX;
        for (size_t src=0;src<nv;src++) if(map[src]>=0) {
            size_t dst=(size_t)map[src];
            memcpy(&sv[dst*3],&verts[src*3],3*sizeof *sv);
            memcpy(&su[dst*2],&uv[src*2],2*sizeof *su);
            if(su[dst*2]<u_min)u_min=su[dst*2];
            if(su[dst*2]>u_max)u_max=su[dst*2];
            if(su[dst*2+1]<v_min)v_min=su[dst*2+1];
            if(su[dst*2+1]>v_max)v_max=su[dst*2+1];
            nfill+=filled[src]?1u:0u;
        }
        size_t fn=0;
        for(size_t f=0;f<nf;f++) {
            int32_t a=map[(size_t)faces[f*3]],b=map[(size_t)faces[f*3+1]];
            int32_t cc=map[(size_t)faces[f*3+2]];
            if(a<0||b<0||cc<0)continue;
            sf[fn*3]=a;sf[fn*3+1]=b;sf[fn*3+2]=cc;fn++;
        }

        size_t serial = part_serial[part];
        char stem[256], obj[1600], vmesh[1600];
        snprintf(stem, sizeof stem, "snip_%06zu_g%ld_p%04d", serial, group, part);
        snprintf(obj, sizeof obj, "%s/%s.obj", dir, stem);
        snprintf(vmesh, sizeof vmesh, "%s/%s.vmesh", dir, stem);
        if (ObjIO_write_uv(obj, sv, snv, sf, snf, su) != 0 ||
            MeshBin_write(vmesh, sv, snv, sf, snf, su) != 0) {
            free(map);free(sv);free(su);free(sf);free(part_serial);return -1;
        }

        if (!*manifest_first) fprintf(manifest, ",\n");
        *manifest_first = 0;
        fprintf(manifest,
                "  {\"id\":%zu,\"component_index\":%zu,\"group\":%ld,"
                "\"representative_track\":%ld,\"part\":%d,\"parts\":%d,"
                "\"obj\":\"%s.obj\",\"vmesh\":\"%s.vmesh\","
                "\"global_row0\":%d,\"global_col0\":%d,\"height\":%d,\"width\":%d,"
                "\"row_stride\":%d,\"col_stride\":%d,"
                "\"cell_col_begin\":%d,\"cell_col_end\":%d,"
                "\"left_snip\":", serial, component_index, group, representative_track,
                part, nparts, stem, stem, gr0, gc0, H, SW,
                row_stride, col_stride, gc0, gc1);
        if (part > 0 && part_serial[part-1] != SIZE_MAX)
            fprintf(manifest, "%zu", part_serial[part-1]);
        else fputs("null", manifest);
        fputs(",\"right_snip\":", manifest);
        if (part+1 < nparts && part_serial[part+1] != SIZE_MAX)
            fprintf(manifest, "%zu", part_serial[part+1]);
        else fputs("null", manifest);
        fputs(",\"shared_left_col\":", manifest);
        if (part > 0 && part_serial[part-1] != SIZE_MAX)
            fprintf(manifest, "%d", gc0);
        else fputs("null", manifest);
        fputs(",\"shared_right_col\":", manifest);
        if (part+1 < nparts && part_serial[part+1] != SIZE_MAX)
            fprintf(manifest, "%d", gc1);
        else fputs("null", manifest);
        fprintf(manifest,
                ",\"u_min\":%.9g,\"u_max\":%.9g,\"v_min\":%.9g,\"v_max\":%.9g,"
                "\"vertices\":%zu,\"faces\":%zu,\"fitted_vertices\":%zu,"
                "\"filled_vertices\":%zu}",
                u_min,u_max,v_min,v_max,snv,snf,snv-nfill,nfill);
        fflush(manifest);
        *out_nv += snv; *out_nf += snf;
        free(map);free(sv);free(su);free(sf);
    }
    *serial_io += active_parts;
    free(part_serial);
    return 0;
}

/* One independently reconstructable component is also one independently
 * bakeable artifact.  Snips remain the streaming representation, while this
 * full component mesh is the direct input for a component Sheet PNG.  Keep its
 * original/global UV here; the combined review mesh may be shelf-packed later. */
static int emit_component_artifact(const char *dir, size_t component_index,
                                   long group, const float *verts, size_t nv,
                                   const int32_t *faces, size_t nf,
                                   const float *uv)
{
    char stem[256], obj[1600], vmesh[1600];
    if (dir == NULL || verts == NULL || faces == NULL || uv == NULL ||
        nv == 0 || nf == 0)
        return -1;
    snprintf(stem, sizeof stem, "component_%06zu_g%ld",
             component_index, group);
    snprintf(obj, sizeof obj, "%s/%s.obj", dir, stem);
    snprintf(vmesh, sizeof vmesh, "%s/%s.vmesh", dir, stem);
    return ObjIO_write_uv(obj, verts, nv, faces, nf, uv) == 0 &&
           MeshBin_write(vmesh, verts, nv, faces, nf, uv) == 0 ? 0 : -1;
}

/* Reproduce the aggregate's deliberately un-fitted topology domains exactly.
 * These small components are pass-through geometry, not new continuation
 * solves: retain their upstream XYZ (while preserving support/generated
 * identity) and connect only four same-label corners.  As in the Python
 * reference assembler, an adjacent-row quad has priority; gap-two is merely a
 * fallback for alternating axial samples, and any one top-left cell connects
 * down at most once. */
static size_t passthrough_faces_for_label(
        const size_t *cells, size_t ncells, int H, int W,
        int32_t *index, uint8_t *connected, int32_t *faces) {
    size_t nf = 0;
    for (int gap = 1; gap <= 2; gap++) {
        for (size_t k = 0; k < ncells; k++) {
            size_t cell = cells[k];
            int r = (int)(cell / (size_t)W);
            int c = (int)(cell % (size_t)W);
            if (connected[cell] || r + gap >= H || c + 1 >= W) continue;
            size_t pa = cell, pb = cell + 1u;
            size_t pc = cell + (size_t)gap * (size_t)W;
            size_t pd = pc + 1u;
            int32_t a = index[pa], b = index[pb];
            int32_t cc = index[pc], d = index[pd];
            if (a < 0 || b < 0 || cc < 0 || d < 0) continue;
            faces[nf * 3u] = a;
            faces[nf * 3u + 1u] = b;
            faces[nf * 3u + 2u] = d;
            nf++;
            faces[nf * 3u] = a;
            faces[nf * 3u + 1u] = d;
            faces[nf * 3u + 2u] = cc;
            nf++;
            connected[cell] = 1;
        }
    }
    return nf;
}

static int passthrough_faces_selftest(void) {
    enum { H = 4, W = 3, N = H * W };
    /* Same label at rows 0,2,3 and columns 0,1.  The adjacent 2->3 quad is
     * selected first, then the gap-two 0->2 fallback: four triangles total. */
    size_t cells[] = {0, 1, 6, 7, 9, 10};
    int32_t index[N];
    uint8_t connected[N] = {0};
    int32_t faces[12];
    for (int i = 0; i < N; i++) index[i] = -1;
    for (size_t i = 0; i < sizeof cells / sizeof cells[0]; i++)
        index[cells[i]] = (int32_t)i;
    size_t nf = passthrough_faces_for_label(
        cells, sizeof cells / sizeof cells[0], H, W,
        index, connected, faces);
    int fail = nf != 4 || !connected[0] || !connected[6] || connected[9];
    fprintf(stderr,
        "[selftest] passthrough topology keeps adjacent-row priority and "
        "gap-two fallback: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}

/* Append aggregate domains below the fitting threshold.  The report producer
 * currently writes one aggregate tifxyz shard (and the Python reference makes
 * the same contract explicit), so fail closed rather than silently omitting
 * data if that format changes. */
static int append_passthrough_components(
        const JsonValue *report, const char *aggregate, double grid_du,
        int u_column_range, int u_column_first, int u_column_end,
        int u_column_halo, int omit_mesh, int repair_folds,
        const char *snips_dir, FILE *snip_manifest,
        int *snip_manifest_first, size_t *snip_serial, int snip_cols,
        size_t component_base, double pack_width,
        double *pack_u, double *pack_v, double *pack_row_h,
        FVec *verts, FVec *uv, IVec *faces, BVec *filled, IVec *compv,
        FVec *corr_uv, IVec *corr_namespace, IVec *corr_lineage,
        size_t *streamed_nv, size_t *streamed_nf,
        PassthroughAtlas **out_atlas, size_t *out_atlas_n,
        PassthroughStats *out_stats, FinalFoldRepairStats *fold_stats) {
    PassthroughStats stats = {0};
    const JsonValue *labels_json = Json_object_get(report, "passthrough_components");
    size_t nlabels = Json_array_len(labels_json);
    *out_atlas = NULL;
    *out_atlas_n = 0;
    *out_stats = stats;
    if (nlabels == 0) return 0;

    long max_label = -1;
    for (size_t i = 0; i < nlabels; i++) {
        long label = Json_as_long(Json_array_get(labels_json, i), -1);
        if (label < 0 || label > INT32_MAX) {
            fprintf(stderr, "[passthrough] invalid topology label %ld\n", label);
            return -1;
        }
        if (label > max_label) max_label = label;
    }
    /* Topology labels are dense component ordinals.  This also prevents a bad
     * manifest from turning one label into an unbounded allocation request. */
    if (max_label > 100000000L) {
        fprintf(stderr, "[passthrough] unreasonable maximum label %ld\n", max_label);
        return -1;
    }
    size_t label_slots = (size_t)max_label + 1u;
    uint8_t *wanted = (uint8_t *)calloc(label_slots, 1);
    size_t *counts = (size_t *)calloc(label_slots, sizeof *counts);
    size_t *offsets = (size_t *)calloc(label_slots + 1u, sizeof *offsets);
    size_t *cursor = (size_t *)calloc(label_slots, sizeof *cursor);
    PassthroughAtlas *atlas = (PassthroughAtlas *)calloc(
        nlabels ? nlabels : 1u, sizeof *atlas);
    if (!wanted || !counts || !offsets || !cursor || !atlas) {
        fprintf(stderr, "[passthrough] out of memory for label index\n");
        free(wanted); free(counts); free(offsets); free(cursor); free(atlas);
        return -1;
    }
    for (size_t i = 0; i < nlabels; i++) {
        size_t label = (size_t)Json_as_long(Json_array_get(labels_json, i), -1);
        if (wanted[label]) {
            fprintf(stderr, "[passthrough] duplicate topology label %zu\n", label);
            free(wanted); free(counts); free(offsets); free(cursor); free(atlas);
            return -1;
        }
        wanted[label] = 1;
    }

    const JsonValue *shards = Json_object_get(report, "shards");
    if (Json_array_len(shards) != 1) {
        fprintf(stderr,
            "[passthrough] expected exactly one aggregate tifxyz shard, got %zu\n",
            Json_array_len(shards));
        free(wanted); free(counts); free(offsets); free(cursor); free(atlas);
        return -1;
    }
    const char *reported_shard = Json_as_string(
        Json_object_get(Json_array_get(shards, 0), "path"));
    if (!reported_shard || !reported_shard[0]) {
        fprintf(stderr, "[passthrough] aggregate shard has no path\n");
        free(wanted); free(counts); free(offsets); free(cursor); free(atlas);
        return -1;
    }
    char shard[1600];
    int absolute = reported_shard[0] == '/' || reported_shard[0] == '\\' ||
                   (reported_shard[0] && reported_shard[1] == ':');
    if (absolute) snprintf(shard, sizeof shard, "%s", reported_shard);
    else snprintf(shard, sizeof shard, "%s/%s", aggregate, reported_shard);

    Arena_T image_arena = Arena_new();
    int32_t *turn = NULL;
    uint8_t *mask = NULL;
    float *gx = NULL, *gy = NULL, *gz = NULL;
    int W = 0, H = 0, iw = 0, ih = 0, md = 0, mh = 0, mw = 0;
    char path[1800];
    snprintf(path, sizeof path, "%s/tifxyz/turn.tif", shard);
    if (TiffIO_load_int32_2d(image_arena, path, &turn, &W, &H) != 0) {
        fprintf(stderr, "[passthrough] cannot load signed labels %s\n", path);
        goto fail;
    }
    snprintf(path, sizeof path, "%s/tifxyz/mask.tif", shard);
    if (TiffIO_load(image_arena, path, &mask, &md, &mh, &mw) != 0 ||
        md != 1 || mw != W || mh != H) {
        fprintf(stderr, "[passthrough] invalid support mask %s\n", path);
        goto fail;
    }
    snprintf(path, sizeof path, "%s/tifxyz/x.tif", shard);
    if (TiffIO_load_float2d(image_arena, path, &gx, &iw, &ih) != 0 ||
        iw != W || ih != H) {
        fprintf(stderr, "[passthrough] invalid x coordinate map %s\n", path);
        goto fail;
    }
    snprintf(path, sizeof path, "%s/tifxyz/y.tif", shard);
    if (TiffIO_load_float2d(image_arena, path, &gy, &iw, &ih) != 0 ||
        iw != W || ih != H) {
        fprintf(stderr, "[passthrough] invalid y coordinate map %s\n", path);
        goto fail;
    }
    snprintf(path, sizeof path, "%s/tifxyz/z.tif", shard);
    if (TiffIO_load_float2d(image_arena, path, &gz, &iw, &ih) != 0 ||
        iw != W || ih != H) {
        fprintf(stderr, "[passthrough] invalid z coordinate map %s\n", path);
        goto fail;
    }

    long long range_first = INT_MIN, range_end = INT_MAX;
    if (u_column_range) {
        range_first = (long long)u_column_first - (long long)u_column_halo;
        range_end = (long long)u_column_end + (long long)u_column_halo;
    }
    size_t N = (size_t)H * (size_t)W, total = 0, supported_total = 0;
    for (size_t pos = 0; pos < N; pos++) {
        int32_t label = turn[pos];
        int c = (int)(pos % (size_t)W);
        if (label < 0 || (size_t)label >= label_slots || !wanted[label] ||
            (long long)c < range_first || (long long)c >= range_end) continue;
        if (!isfinite((double)gx[pos]) || !isfinite((double)gy[pos]) ||
            !isfinite((double)gz[pos])) {
            fprintf(stderr,
                "[passthrough] label %d has non-finite XYZ at row=%zu col=%d\n",
                label, pos / (size_t)W, c);
            goto fail;
        }
        counts[label]++;
        total++;
        supported_total += mask[pos] != 0;
    }
    if (!u_column_range) {
        long expected = Json_member_long(report, "passthrough_cells", -1);
        if (expected < 0 || (size_t)expected != supported_total) {
            fprintf(stderr,
                "[passthrough] support/report count mismatch: loaded=%zu "
                "reported=%ld (geometry cells=%zu)\n",
                supported_total, expected, total);
            goto fail;
        }
    }
    for (size_t label = 0; label < label_slots; label++)
        offsets[label + 1u] = offsets[label] + counts[label];
    memcpy(cursor, offsets, label_slots * sizeof *cursor);
    size_t *cells = (size_t *)malloc((total ? total : 1u) * sizeof *cells);
    int32_t *index = (int32_t *)malloc((N ? N : 1u) * sizeof *index);
    uint8_t *connected = (uint8_t *)calloc(N ? N : 1u, 1);
    if (!cells || !index || !connected) {
        fprintf(stderr, "[passthrough] out of memory for %zu source cells\n", total);
        free(cells); free(index); free(connected);
        goto fail;
    }
    memset(index, 0xff, N * sizeof *index);
    for (size_t pos = 0; pos < N; pos++) {
        int32_t label = turn[pos];
        int c = (int)(pos % (size_t)W);
        if (label >= 0 && (size_t)label < label_slots && wanted[label] &&
            (long long)c >= range_first && (long long)c < range_end)
            cells[cursor[label]++] = pos;
    }

    for (size_t ordinal = 0; ordinal < nlabels; ordinal++) {
        size_t label = (size_t)Json_as_long(Json_array_get(labels_json, ordinal), -1);
        size_t ncells = counts[label];
        if (ncells == 0) continue;
        if (ncells > (size_t)INT32_MAX) {
            fprintf(stderr, "[passthrough] label %zu exceeds 32-bit mesh indexing\n", label);
            free(cells); free(index); free(connected);
            goto fail;
        }
        float *lv = (float *)malloc(ncells * 3u * sizeof *lv);
        float *lu = (float *)malloc(ncells * 2u * sizeof *lu);
        uint8_t *lf = (uint8_t *)calloc(ncells, 1);
        int32_t *lfaces = (int32_t *)malloc(ncells * 6u * sizeof *lfaces);
        if (!lv || !lu || !lf || !lfaces) {
            fprintf(stderr, "[passthrough] out of memory for label %zu\n", label);
            free(lv); free(lu); free(lf); free(lfaces);
            free(cells); free(index); free(connected);
            goto fail;
        }
        int min_r = H, max_r = -1, min_c = W, max_c = -1;
        for (size_t k = 0; k < ncells; k++) {
            size_t pos = cells[offsets[label] + k];
            int r = (int)(pos / (size_t)W), c = (int)(pos % (size_t)W);
            index[pos] = (int32_t)k;
            lv[k * 3u] = gz[pos];
            lv[k * 3u + 1u] = gy[pos];
            lv[k * 3u + 2u] = gx[pos];
            lu[k * 2u] = (float)((double)c * grid_du);
            lu[k * 2u + 1u] = (float)r;
            lf[k] = mask[pos] ? 0 : 1;
            if (r < min_r) min_r = r; if (r > max_r) max_r = r;
            if (c < min_c) min_c = c; if (c > max_c) max_c = c;
        }
        size_t nf = passthrough_faces_for_label(
            cells + offsets[label], ncells, H, W,
            index, connected, lfaces);
        size_t component_index = component_base + ordinal;
        if (repair_folds && repair_final_fold_flaps(
                lv, ncells, lfaces, &nf, 0, 0,
                "passthrough-label", (long)label, fold_stats) != 0) {
            free(lv); free(lu); free(lf); free(lfaces);
            free(cells); free(index); free(connected);
            goto fail;
        }

        if (snip_manifest && emit_component_artifact(
                snips_dir, component_index, -(long)label - 1L,
                lv, ncells, lfaces, nf, lu) != 0) {
            fprintf(stderr,
                    "[passthrough] failed to emit label %zu component artifact\n",
                    label);
            free(lv); free(lu); free(lf); free(lfaces);
            free(cells); free(index); free(connected);
            goto fail;
        }
        if (snip_manifest && emit_rect_snips(
                snips_dir, snip_manifest, snip_manifest_first, snip_serial,
                snip_cols, component_index, -(long)label - 1L, -1,
                max_r - min_r + 1, max_c - min_c + 1,
                min_r, min_c, grid_du, 1, 1,
                lv, ncells, lfaces, nf, lu, lf,
                streamed_nv, streamed_nf) != 0) {
            fprintf(stderr, "[passthrough] failed to emit label %zu snips\n", label);
            free(lv); free(lu); free(lf); free(lfaces);
            free(cells); free(index); free(connected);
            goto fail;
        }

        double u0 = lu[0], u1 = lu[0], v0 = lu[1], v1 = lu[1];
        for (size_t k = 1; k < ncells; k++) {
            double u = lu[k * 2u], v = lu[k * 2u + 1u];
            if (u < u0) u0 = u; if (u > u1) u1 = u;
            if (v < v0) v0 = v; if (v > v1) v1 = v;
        }
        if (pack_width > 0.0 && !omit_mesh) {
            double pw = u1 - u0, ph = v1 - v0, gap = 4.0;
            if (*pack_u > 0.0 && *pack_u + pw > pack_width) {
                *pack_u = 0.0;
                *pack_v += *pack_row_h + gap;
                *pack_row_h = 0.0;
            }
            double su = *pack_u - u0, sv = *pack_v - v0;
            for (size_t k = 0; k < ncells; k++) {
                lu[k * 2u] += (float)su;
                lu[k * 2u + 1u] += (float)sv;
            }
            u0 = *pack_u; u1 = *pack_u + pw;
            v0 = *pack_v; v1 = *pack_v + ph;
            *pack_u += pw + gap;
            if (ph > *pack_row_h) *pack_row_h = ph;
        }
        PassthroughAtlas *pa = &atlas[*out_atlas_n];
        pa->component_index = component_index;
        pa->label = (long)label;
        pa->height = max_r - min_r + 1;
        pa->width = max_c - min_c + 1;
        pa->u_min = u0; pa->u_max = u1;
        pa->v_min = v0; pa->v_max = v1;
        (*out_atlas_n)++;

        if (!omit_mesh) {
            size_t base = verts->n / 3u;
            if (base > (size_t)INT32_MAX || ncells > (size_t)INT32_MAX - base) {
                fprintf(stderr, "[passthrough] combined mesh exceeds 32-bit indexing\n");
                free(lv); free(lu); free(lf); free(lfaces);
                free(cells); free(index); free(connected);
                goto fail;
            }
            fv_reserve(verts, ncells * 3u);
            memcpy(verts->v + verts->n, lv, ncells * 3u * sizeof *lv);
            verts->n += ncells * 3u;
            fv_reserve(uv, ncells * 2u);
            memcpy(uv->v + uv->n, lu, ncells * 2u * sizeof *lu);
            uv->n += ncells * 2u;
            bv_reserve(filled, ncells);
            memcpy(filled->v + filled->n, lf, ncells * sizeof *lf);
            filled->n += ncells;
            iv_reserve(compv, ncells);
            for (size_t k = 0; k < ncells; k++)
                compv->v[compv->n + k] = (int32_t)component_index;
            compv->n += ncells;
            if (corr_uv && corr_namespace && corr_lineage) {
                fv_reserve(corr_uv, ncells * 2u);
                iv_reserve(corr_namespace, ncells);
                iv_reserve(corr_lineage, ncells);
                for (size_t k = 0; k < ncells; k++) {
                    size_t pos = cells[offsets[label] + k];
                    int r = (int)(pos / (size_t)W);
                    int c = (int)(pos % (size_t)W);
                    corr_uv->v[corr_uv->n++] = (float)((double)c * grid_du);
                    corr_uv->v[corr_uv->n++] = (float)r;
                    corr_namespace->v[corr_namespace->n++] = 1;
                    corr_lineage->v[corr_lineage->n++] = (int32_t)label;
                }
            }
            iv_reserve(faces, nf * 3u);
            for (size_t k = 0; k < nf * 3u; k++)
                faces->v[faces->n + k] = lfaces[k] + (int32_t)base;
            faces->n += nf * 3u;
        }

        for (size_t k = 0; k < ncells; k++) {
            size_t pos = cells[offsets[label] + k];
            index[pos] = -1;
            connected[pos] = 0;
        }
        stats.components++;
        stats.vertices += ncells;
        for (size_t k = 0; k < ncells; k++) {
            stats.supported_vertices += lf[k] == 0;
            stats.generated_vertices += lf[k] != 0;
        }
        stats.faces += nf;
        free(lv); free(lu); free(lf); free(lfaces);
    }

    free(cells); free(index); free(connected);
    Arena_dispose(&image_arena);
    free(wanted); free(counts); free(offsets); free(cursor);
    *out_atlas = atlas;
    *out_stats = stats;
    fprintf(stderr,
        "[passthrough] retained %zu topology domains: %zu vertices "
        "(%zu supported, %zu generated), %zu faces\n",
        stats.components, stats.vertices, stats.supported_vertices,
        stats.generated_vertices, stats.faces);
    return 0;

fail:
    Arena_dispose(&image_arena);
    free(wanted); free(counts); free(offsets); free(cursor); free(atlas);
    return -1;
}

/* Rasterize coarse mesh triangles back onto the native material vertex
 * lattice used by the topology/bake audit.  This keeps a coarse solver mesh
 * from masquerading as holes merely because vertices no longer exist at every
 * native UV sample. */
static void topology_raster_triangle(
        uint8_t *rgb,long image_w,long image_h,long min_gc,long min_gr,
        double grid_du,const float *uv,const uint8_t *filled,
        const int32_t *tri,int32_t component) {
    double x[3],y[3];
    for(int k=0;k<3;k++) {
        size_t v=(size_t)tri[k];
        x[k]=(double)uv[v*2]/grid_du-(double)min_gc;
        y[k]=(double)uv[v*2+1]-(double)min_gr;
    }
    double den=(y[1]-y[2])*(x[0]-x[2])+(x[2]-x[1])*(y[0]-y[2]);
    if(fabs(den)<1e-12)return;
    long x0=(long)floor(fmin(x[0],fmin(x[1],x[2]))-1e-8);
    long x1=(long)ceil (fmax(x[0],fmax(x[1],x[2]))+1e-8);
    long y0=(long)floor(fmin(y[0],fmin(y[1],y[2]))-1e-8);
    long y1=(long)ceil (fmax(y[0],fmax(y[1],y[2]))+1e-8);
    if(x0<0)x0=0;if(y0<0)y0=0;
    if(x1>=image_w)x1=image_w-1;if(y1>=image_h)y1=image_h-1;
    unsigned h=(unsigned)component*2654435761u;
    uint8_t cr=(uint8_t)(60+(h&127u));
    uint8_t cg=(uint8_t)(60+((h>>8)&127u));
    uint8_t cb=(uint8_t)(60+((h>>16)&127u));
    for(long py=y0;py<=y1;py++)for(long px=x0;px<=x1;px++) {
        double a=((y[1]-y[2])*((double)px-x[2])+
                  (x[2]-x[1])*((double)py-y[2]))/den;
        double b=((y[2]-y[0])*((double)px-x[2])+
                  (x[0]-x[2])*((double)py-y[2]))/den;
        double c=1.0-a-b;
        if(a < -1e-7 || b < -1e-7 || c < -1e-7)continue;
        double fill=a*(double)(filled[(size_t)tri[0]]!=0)+
                    b*(double)(filled[(size_t)tri[1]]!=0)+
                    c*(double)(filled[(size_t)tri[2]]!=0);
        size_t o=((size_t)py*(size_t)image_w+(size_t)px)*3u;
        if(fill>=0.5){rgb[o]=40;rgb[o+1]=210;rgb[o+2]=70;}
        else {rgb[o]=cr;rgb[o+1]=cg;rgb[o+2]=cb;}
    }
}

static void usage(const char *p) {
    fprintf(stderr,
        "usage: %s <aggregate> <tracks_root> <out.obj> [options]\n"
        "  --fallback-tracks-root <dir>   second root to resolve tracks from\n"
        "  --registration <reg.json>      merge same-sheet charts into super-components\n"
        "  --registration-max-gap <n>     max developed-grid gap allowed in a merge (default 8)\n"
        "  --legacy-v92-registration      allow only the audited v92 chart-merge-v1 registration\n"
        "  --correspondence <out.qrcorr>  write exact pre-pack lattice key per VMesh vertex\n"
        "  --bake-refit-masks <dir>       broad-dark HxW masks from bake audit\n"
        "  --only-group <n>               assemble one registration group (A/B/debug)\n"
        "  --u-column-range <c0> <c1>     solve only half-open global UV-column range [c0,c1)\n"
        "  --u-column-halo <n>             include n hidden solve columns on each side of the requested range\n"
        "  --giant-ribbon                 experimental single global u/v lattice\n"
        "  --snips-dir <dir>              write variable-length OBJ/VMESH ribbon snips\n"
        "  --snip-cols <n>                max quad-cell columns per snip (default 512)\n"
        "  --snips-only                   do not retain/write the monolithic mesh\n"
        "  --report-only                  run assembly/CT diagnostics without writing meshes\n"
        "  --pack-components <u-width>    shelf-pack disconnected output charts (0=global UV)\n"
        "  --grid-du <f>                  voxels per u-column (default 2.0)\n"
        "  --ribbon-u-stride <n>          confidence-fit coarse ribbon every n U columns (default 1)\n"
        "  --ribbon-v-stride <n>          confidence-fit coarse ribbon every n V rows (default 1)\n"
        "  --solid-mode rect|span|bounded full rectangle (default), span hull, or fitted footprint+closed holes\n"
        "  --max-hole-distance <n>        max closed-hole depth in grid cells (default 16; 0 disables hole fill)\n"
        "  --max-row-gap <n>              nearest later live row usable for quads (default 1; ragged fits use 2)\n"
        "  --fit-stiffness <f>            fitted-cell anchor weight (default 1e4)\n"
        "  --cyl-axis <y> <x>             winding-aware harmonic fill around scroll axis (requires lifted phase)\n"
        "  --snap-axis <y> <x>            orient CT/ray correction around this axis without cylindrical refill\n"
        "  --init-quarantine-sander <f>   excise broad initializer distortion >= f (default 0=off)\n"
        "  --init-quarantine-tile <n>     quarantine tile edge in grid cells (default 16)\n"
        "  --init-quarantine-fraction <f> severe-triangle fraction per tile (default .20)\n"
        "  --init-quarantine-dilate <n>   tile-ring dilation around defects (default 3)\n"
        "  --require-collision-free       abort unless elastic-shell repair removes every long/stab contact\n"
        "                                 (default retains its best exact-audited partial state)\n"
        "  --repair-folds                 delete only smaller edge-fold flaps, then require zero residual folds\n"
        "  --collision-patience <n>        accepted non-best shell rounds before partial stop (default 12; 0=all rounds)\n"
        "  --collision-rounds <n>          maximum elastic-shell active-set rounds (default 64; max 256)\n"
        "  --collision-collar <f>          base speculative contact collar, vox (default 4.0; stall escalates to base+8)\n"
        "  --collision-displacement-bound <f>  radial displacement bound, vox (default 128)\n"
        "  --collision-settle-rounds <n>   settle-toward-rest rounds after a conflict-free shell (default 16; 0=off)\n"
        "  --collision-settle-beta <f>     per-round pull of radial displacement toward the rest coil (default 0.70)\n"
        "  --continuation-trim             delete continuation quads whose CT stayed black through every\n"
        "                                 snap pass (no non-dark sample, no supported alternative); glimmer kept\n"
        "  --trim-min-region <n>          smallest black region to trim, cells (default 24)\n"
        "  --trim-rim-guard <n>           keep black cells within this Chebyshev ring of fitted cells (default 2)\n"
        "  --relax-rounds <n>             ARAP developable fill relaxation (default 0)\n"
        "  --arap-only                    trust unwrap support, metric-ARAP full ribbon, then stop geometry\n"
        "  --trust-unwrap-domain          use finite nonzero-confidence direct support; gates remain diagnostics only\n"
        "  --arap-report <out.json>       metric/error/convergence report by component\n"
        "  --arap-max-iterations <n>      local/global iterations (default 1000)\n"
        "  --arap-tolerance <f>           stop at maximum XYZ movement (default 1e-5 vox)\n"
        "  --arap-tolerance-rel <f>       scale-aware stop: max(abs, rel*first-iter move) (default 1e-3; 0=legacy)\n"
        "  --arap-stall-window <n>        stop as stalled after n non-improving iterations (default 50; 0=off)\n"
        "  --arap-stall-fraction <f>      stall improvement floor (default 0.01)\n"
        "  --arap-cg-iterations <n>       PCG steps per coordinate/global step (default 256)\n"
        "  --arap-cg-tolerance <f>        PCG relative residual target (default 1e-8)\n"
        "  --arap-mg-levels <n>           Galerkin hierarchy depth incl. fine (default 16; 1=off)\n"
        "  --arap-mg-cycles <n>           V-cycles per PCG preconditioner application (default 1)\n"
        "  --arap-mg-pre-sweeps <n>       forward four-color smoothing sweeps (default 2)\n"
        "  --arap-mg-post-sweeps <n>      reverse four-color smoothing sweeps (default 2)\n"
        "  --arap-mg-coarse-sweeps <n>    coarsest symmetric sweeps (default 32)\n"
        "  --arap-nlmg-levels <n>          nonlinear coarse-to-fine levels incl. fine (default 4)\n"
        "  --arap-nlmg-iterations <n>      ARAP rounds on each coarse level (default 128)\n"
        "  --arap-nlmg-min-vertices <n>    activation threshold on fine grid (default 4000000)\n"
        "  --arap-source-weight <f>       robust direct-support attachment (default 1)\n"
        "  --arap-fill-weight <f>         weak initializer attachment (default .01)\n"
        "  --arap-huber <f>               support attachment Huber displacement (default 2 vox)\n"
        "  --arap-min-source-scale <f>    robust source-weight floor fraction (default .05)\n"
        "  --arap-min-area-ratio <f>      per-step oriented triangle barrier (default .05)\n"
        "  --arap-max-step <f>            local XYZ trust radius per outer round (default 2 vox)\n"
        "  --arap-axis-normal-weight <f>  soft radial-normal prior (0..1) forbidding axial 'bucket-lid'\n"
        "                                 collapse; requires --cyl-axis; 0 = ordinary ARAP (default 0)\n"
        "  --arap-axis-normal-radius <f>  fade the radial-normal prior to 0 past this radius (default 0=const)\n"
        "  --arap-axis-fold-limit <f>     no-fold guard: reject a step raising |n.z|/|n| past this (default 0=off)\n"
        "  --post-snap-arap               re-fit the accepted CT snap with metric ARAP; every\n"
        "                                 snapped vertex is a soft anchor, z is exact with --cyl-axis,\n"
        "                                 and an exact no-new-intersections transaction selects the step\n"
        "  --raw-cubes <dir>              cubes_RAW or uncompressed RAW Zarr; CT snap\n"
        "  --require-complete-raw         fail if any required RAW chunk is absent (default)\n"
        "  --allow-incomplete-raw         diagnostics only; never publish this output\n"
        "  --snap-reach <f>               snap march reach, vox (default 6)\n"
        "  --snap-bright-min <f>          min ridge intensity 0..255 (default auto)\n"
        "  --snap-band-frac <f>           auto floor fraction in fitted p1..p99 (default .30)\n"
        "  --snap-gain <f>                min intensity gain to snap (default 20)\n"
        "  --snap-ridge-lock <f>          fitted local-ridge lock radius (default .75 vox)\n"
        "  --snap-dark-probe-frac <f>     fitted p1..p99 bake-dark diagnostic (default .20; 0=off)\n"
        "  --snap-guided-reach <f>        winding-gated residual search reach (default 12)\n"
        "  --snap-winding-tube <f>        minimum adaptive half-wrap corridor (default 3)\n"
        "  --snap-patch-recover           recover connected raster-dark fill patches as one winding decision\n"
        "  --snap-no-patch-recover        disable connected dark-patch recovery (library default)\n"
        "  --snap-patch-buffer <n>        fixed/movable collar rings around each dark patch (default 2)\n"
        "  --snap-patch-min-cells <n>     minimum connected dark-cell area (default 1)\n"
        "  --snap-patch-ray-reach <f>     bidirectional patch ray reach in voxels (default 48)\n"
        "  --snap-patch-ray-margin <f>    +/- first-hit decisiveness margin (default 1 voxel)\n"
        "  --snap-patch-relax-sweeps <n> screened-Laplacian seam-blend sweeps (default 64)\n"
        "  --snap-patch-target-weight <f> direct ray-target screen weight (default 16)\n"
        "  --snap-patch-fitted-dark-slack <n> fitted seam threshold toggles allowed (default 512)\n"
        "  --snap-no-local-contrast       disable profile-local transition normalization\n"
        "  --snap-bake-window <lo> <hi>   fixed RAW window for in-pass raster gate (default 47 193)\n"
        "  --snap-bake-dark <n>           post-window dark threshold for pass gate (default 13)\n"
        "  --snap-side midline|recto|nearest  CT target policy: midline = center of the bright papyrus band\n"
        "                                 (default midline with --cyl-axis); recto = legacy inward bright->dark face\n"
        "  --snap-quilt-smooth <f>        physical 3-D offset neighbour weight (default 4)\n"
        "  --snap-quilt-seam-smooth <f>   fitted/filled edge multiplier (default 1)\n"
        "  --snap-quilt-seam-refine <f>   diagnostic legacy seam feedback (default 1=off)\n"
        "  --snap-quilt-anchor <f>        no-ridge zero-depth anchor (default .01)\n"
        "  --snap-quilt-huber <f>         raw-offset coherence Huber scale (default 1.5 vox)\n"
        "  --snap-quilt-sweeps <n>        red/black offset-field sweeps (default 4)\n"
        "  --snap-rounds <n>              gentle iterated snap rounds (default 8)\n"
        "  --snap-passes <n>              proximal refit/snap passes; reset init anchor each pass (default 2)\n"
        "  --snap-policy modern|v92       select modern defaults or the frozen v92 bake contract\n"
        "  --snap-report <out.json>       machine-readable pass accept/rollback diagnostics\n"
        "  --snap-refine-anchor <w>       later-pass soft on-target anchor multiplier (default 4; 0=off)\n"
        "  --snap-fitted-source-anchor <w> soft original-position multiplier for fitted support\n"
        "                                 (default 0=legacy hard freeze; production uses 64)\n"
        "  --snap-sweeps <n>              red/black position sweeps per round (default 2)\n"
        "  --snap-omega <f>               red/black under-relaxation (default .5)\n"
        "  --snap-mg-levels <n>           hierarchy levels incl. fine grid (default 1=off; max 16)\n"
        "  --snap-mg-cycles <n>           coarse corrections per scalar/position solve (default 1)\n"
        "  --snap-mg-patch <n>            coarse-cell correction patch width (default 16; 0=global)\n"
        "  --remesh-motion <f>            advise remesh at tangential p95 cells (default .5)\n"
        "  --remesh-motion-max <f>        advise on localized tangential motion (default 4 cells)\n"
        "  --remesh-sander <f>            advise remesh at Sander p95 (default 2)\n"
        "  --remesh-over4 <percent>       advise when this percent is >=4x (default .01)\n"
        "  --remesh-sweeps <n>            objective-selected correspondence sweeps at a fired gate (default 1; 0=off)\n"
        "  --snap-w-target/-orig/-smooth  gentle snap energy weights\n"
        "  --topology-map <out.png>       per-component coverage PNG\n"
        "  --selftest\n", p);
}

/* Resolve a track directory across the given roots (mirrors the Python
 * resolve_track_dir variants); writes the found dir into out[cap].
 * report.json is quad_ribbon's atomic commit marker, written only after frozen
 * geometry and strict RAW sampling finish.  Field TIFFs alone may be debris
 * from a power loss and must never shadow a complete fallback fit. */
static int resolve_track_dir(const char **roots, int nroots, long tid,
                             char *out, size_t cap) {
    const char *suff[] = { "track_%03ld", "track_%05ld/quad_ribbon", "track_%05ld" };
    for (int r = 0; r < nroots; r++) {
        for (int s = 0; s < 3; s++) {
            char rel[256];
            snprintf(rel, sizeof rel, suff[s], tid);
            snprintf(out, cap, "%s/%s", roots[r], rel);
            char probe[1200];
            snprintf(probe, sizeof probe, "%s/report.json", out);
            FILE *commit = fopen(probe, "rb");
            if (!commit) continue;
            fclose(commit);
            snprintf(probe, sizeof probe, "%s/tifxyz/field_x.tif", out);
            FILE *f = fopen(probe, "rb");
            if (f) { fclose(f); return 0; }
        }
    }
    return -1;
}

/* Interpolation certification is deliberately separable from the core fitted
 * geometry.  A track with clean frozen source/PDE geometry may still have no
 * holdout evidence for its optional interpolation bank; production can consume
 * its certified mask while ignoring that bank.  Any other failed gate means the
 * upstream surface/parameterization itself is not publishable. */
/* Return 0=reject, 1=certified domain, 2=source-support-only fallback. */
static int geometry_domain_policy(const JsonValue *geom, long tid) {
    if (!geom || !Json_as_bool(Json_object_get(geom, "geometry_frozen"), 0)) {
        fprintf(stderr, "[geometry-gate] reject track %ld: geometry is not frozen\n", tid);
        return 0;
    }
    const JsonValue *gates = Json_object_get(geom, "gates");
    size_t n = Json_array_len(gates);
    if (n == 0) {
        fprintf(stderr, "[geometry-gate] reject track %ld: no geometry gates\n", tid);
        return 0;
    }
    int source_only = 0;
    for (size_t i = 0; i < n; i++) {
        const JsonValue *gate = Json_array_get(gates, i);
        if (Json_as_bool(Json_object_get(gate, "pass"), 0)) continue;
        const char *name = Json_as_string(Json_object_get(gate, "name"));
        if (name && !strcmp(name, "bounded_interpolation_is_certified"))
            continue; /* those optional cells are absent from mask.tif */
        if (name && (!strcmp(name, "structured_holdout_points") ||
                     !strcmp(name, "certified_safe_extrapolation_voxels") ||
                     !strcmp(name, "holdout_p95_at_safe_distance_voxels") ||
                     !strcmp(name, "gap_matched_shadow_holdout_p95_voxels") ||
                     !strcmp(name, "certified_pde_continuation_cells"))) {
            /* No continuation certificate: retain direct source support, but
             * do not consume any PDE/interpolation continuation. */
            source_only = 1;
            continue;
        }
        fprintf(stderr, "[geometry-gate] reject track %ld: failed %s\n",
                tid, name ? name : "unnamed gate");
        return 0;
    }
    if (source_only)
        fprintf(stderr, "[geometry-gate] track %ld: continuation uncertified; "
                        "using source support only\n", tid);
    return source_only ? 2 : 1;
}

/* Load one track fragment into `scratch`.  Returns 0 on success, -2 for a
 * deliberate upstream-geometry quarantine, and -1 for malformed/missing data. */
static int unwrap_support_authoritative(uint8_t support,uint8_t confidence) {
    /* Pass-2 distortion masks and the robust fit both express rejection by
     * driving output confidence to exactly zero while deliberately retaining
     * support for raster/domain bookkeeping.  Treating support alone as a pin
     * silently resurrects those rejected observations in the solid-strip
     * initializer.  Positive confidence remains a Boolean authority gate here;
     * its magnitude belongs to the upstream variational fit, not to this fill. */
    return support==255&&confidence>0;
}

static int unwrap_support_authority_selftest(void) {
    int fail=0;
    fail+=unwrap_support_authoritative(255,255)!=1;
    fail+=unwrap_support_authoritative(255,1)!=1;
    fail+=unwrap_support_authoritative(255,0)!=0;
    fail+=unwrap_support_authoritative(0,255)!=0;
    fprintf(stderr,"[selftest] confidence-weighted unwrap authority: %s\n",
            fail?"FAIL":"PASS");
    return fail;
}

static int load_track(Arena_T scratch, const char *dir,
                       const JsonValue *mtrack, double grid_du, long z_lo,
                       int trust_unwrap_domain, int require_lifted_phase,
                       Loaded *out) {
    char p[1300];
    snprintf(p, sizeof p, "%s/geometry_report.json", dir);
    const JsonValue *geom = Json_parse_file(scratch, p, NULL);
    if (!geom && trust_unwrap_domain) {
        /* Certification/reporting is diagnostic in this mode.  Some legacy
         * reports contain non-JSON NaN spellings, while the compact field
         * metadata remains valid and carries the same lattice origin. */
        snprintf(p, sizeof p, "%s/tifxyz/meta.json", dir);
        geom = Json_parse_file(scratch, p, NULL);
        if (geom)
            fprintf(stderr, "[trust] track %ld: diagnostic geometry report is unparsable; using field metadata for origin only\n",
                    out->tid);
    }
    if (!geom) {
        fprintf(stderr, "[geometry-gate] reject track %ld: cannot read geometry origin from %s\n",
                out->tid, p);
        return -1;
    }
    int geometry_policy = trust_unwrap_domain ? 1 : geometry_domain_policy(geom, out->tid);
    if (!geometry_policy) return -2;

    float *Z = NULL, *Y = NULL, *X = NULL;
    int W = 0, H = 0, W2 = 0, H2 = 0;
    snprintf(p, sizeof p, "%s/tifxyz/field_z.tif", dir);
    if (TiffIO_load_float2d(scratch, p, &Z, &W, &H) != 0) return -1;
    snprintf(p, sizeof p, "%s/tifxyz/field_y.tif", dir);
    if (TiffIO_load_float2d(scratch, p, &Y, &W2, &H2) != 0 || W2 != W || H2 != H) return -1;
    snprintf(p, sizeof p, "%s/tifxyz/field_x.tif", dir);
    if (TiffIO_load_float2d(scratch, p, &X, &W2, &H2) != 0 || W2 != W || H2 != H) return -1;
    float *phase = NULL;
    if (require_lifted_phase) {
        snprintf(p, sizeof p, "%s/tifxyz/reference_phase.tif", dir);
        if (TiffIO_load_float2d(scratch, p, &phase, &W2, &H2) != 0 ||
            W2 != W || H2 != H) {
            fprintf(stderr,
                    "[phase-gate] reject track %ld: missing/invalid upstream "
                    "lifted phase %s\n",out->tid,p);
            return -3;
        }
    }

    uint8_t *cert = NULL, *sup = NULL, *itp = NULL, *confidence = NULL;
    if (!trust_unwrap_domain) {
        int Dc = 0, Hc = 0, Wc = 0;
        snprintf(p, sizeof p, "%s/tifxyz/mask.tif", dir);
        if (TiffIO_load(scratch, p, &cert, &Dc, &Hc, &Wc) != 0 ||
            Dc != 1 || Hc != H || Wc != W) {
            fprintf(stderr, "[trust] refusing track %ld: missing/invalid certified mask %s\n",
                    out->tid, p);
            return -1;
        }
    }
    int Ds = 0, Hs = 0, Ws = 0;
    snprintf(p, sizeof p, "%s/tifxyz/support.tif", dir);
    if (TiffIO_load(scratch, p, &sup, &Ds, &Hs, &Ws) != 0 ||
        Ds != 1 || Hs != H || Ws != W) return -1;
    if (trust_unwrap_domain) {
        int Dq = 0, Hq = 0, Wq = 0;
        snprintf(p, sizeof p, "%s/tifxyz/confidence.tif", dir);
        if (TiffIO_load(scratch, p, &confidence, &Dq, &Hq, &Wq) != 0 ||
            Dq != 1 || Hq != H || Wq != W) {
            fprintf(stderr,
                    "[trust] refusing track %ld: missing/invalid fitted confidence %s\n",
                    out->tid,p);
            return -1;
        }
    }
    int Di = 0, Hi = 0, Wi = 0;
    snprintf(p, sizeof p, "%s/tifxyz/interpolation_domain.tif", dir);
    if (TiffIO_load(scratch, p, &itp, &Di, &Hi, &Wi) != 0 ||
        Di != 1 || Hi != H || Wi != W) itp = NULL;

    size_t HW = (size_t)H * (size_t)W;
    uint8_t *topology = ARENA_ALLOC(scratch, HW);
    uint8_t *domain = ARENA_ALLOC(scratch, HW);
    uint8_t *trusted = ARENA_ALLOC(scratch, HW);
    float   *field  = ARENA_ALLOC(scratch, HW * 3 * sizeof *field);
    size_t ntopology = 0, ndom = 0, ntrusted = 0, ncandidate = 0;
    size_t nignored_policy = 0, nnonfinite_candidate = 0;
    size_t nzero_confidence_support = 0;
    for (size_t c = 0; c < HW; c++) {
        int candidate = (sup[c] == 255) || (itp && itp[c] > 0);
        int finite_field = isfinite((double)Z[c]) && isfinite((double)Y[c]) &&
                           isfinite((double)X[c]);
        int finite_phase = !require_lifted_phase || isfinite((double)phase[c]);
        int policy_footprint = trust_unwrap_domain
                             ? candidate
                             : (geometry_policy == 2
                                  ? sup[c] == 255
                                  : candidate && cert[c] > 0);
        int topo = policy_footprint && finite_field && finite_phase;
        /* The dense fitted/interpolated field is a good initializer and defines
         * where this track actually exists.  It is not a hard observation.
         * Only direct source support (and, when available, positive confidence)
         * receives the high-weight/frozen treatment in the downstream solve. */
        int trust = topo && sup[c] == 255 &&
                    (!trust_unwrap_domain || confidence[c] > 0);
        int d = topo;
        topology[c] = (uint8_t)(topo ? 1 : 0);
        domain[c] = (uint8_t)(d ? 1 : 0);
        trusted[c] = (uint8_t)(trust ? 1 : 0);
        ntopology += (size_t)topo;
        ndom += (size_t)d;
        ntrusted += (size_t)trust;
        ncandidate += (size_t)candidate;
        nignored_policy += (size_t)(candidate && !policy_footprint);
        nnonfinite_candidate +=
            (size_t)(policy_footprint && (!finite_field || !finite_phase));
        nzero_confidence_support +=
            (size_t)(trust_unwrap_domain&&sup[c]==255&&confidence[c]==0);
        field[c * 3 + 0] = Z[c];
        field[c * 3 + 1] = Y[c];
        field[c * 3 + 2] = X[c];
    }
    if (ntopology == 0 || ndom == 0) return -1;
    fprintf(stderr,
            "[trust] track=%ld topology/initializer/trusted=%zu/%zu/%zu of %zu; "
            "support-or-interpolation candidates=%zu, policy-excluded=%zu, "
            "nonfinite-excluded=%zu, zero-confidence-support=%zu; policy=%s\n",
            out->tid,ntopology,ndom,ntrusted,HW,ncandidate,nignored_policy,
            nnonfinite_candidate,nzero_confidence_support,
            trust_unwrap_domain?"full-finite-fitted-field/direct-positive-confidence":
              (geometry_policy==2?"source-support-only":"certified-fitted-field"));

    /* origin_zu[0] from the track's geometry report -> gr0 relative to z_lo. */
    const JsonValue *ozu = Json_object_get(geom, "origin_zu");
    double origin_z = Json_as_double(Json_array_get(ozu, 0), (double)z_lo);

    const JsonValue *members = Json_object_get(mtrack, "members");
    const JsonValue *m0 = Json_array_get(members, 0);
    double u_off = Json_member_double(m0, "global_u_offset", 0.0);

    out->topology = topology;
    out->domain = domain;
    out->trusted = trusted;
    out->field = field;
    out->phase = phase;
    out->H = H;
    out->W = W;
    out->gr0 = (int)lround(origin_z - (double)z_lo);
    out->c0 = (int)lround(-u_off / grid_du);
    return 0;
}

/* Registration proximity alone can join neighboring WRAPS rather than two
 * fragments of one material sheet.  A decisive test is their existing atlas
 * overlap: if coincident grid cells are many voxels apart in 3-D, merging the
 * fields creates a hard last-writer seam. */
static void loaded_overlap_audit(const Loaded *a,const Loaded *b,double compatible,
                                 size_t *out_count,size_t *out_conflict,
                                 double *out_mean,double *out_max) {
    int r0=a->gr0>b->gr0?a->gr0:b->gr0;
    int c0=a->c0>b->c0?a->c0:b->c0;
    int r1=a->gr0+a->H<b->gr0+b->H?a->gr0+a->H:b->gr0+b->H;
    int c1=a->c0+a->W<b->c0+b->W?a->c0+a->W:b->c0+b->W;
    size_t n=0,bad=0;double sum=0.0,mx=0.0;
    for(int r=r0;r<r1;r++)for(int c=c0;c<c1;c++){
        size_t ia=(size_t)(r-a->gr0)*a->W+(size_t)(c-a->c0);
        size_t ib=(size_t)(r-b->gr0)*b->W+(size_t)(c-b->c0);
        if(!a->domain[ia]||!b->domain[ib])continue;
        double dz=(double)a->field[ia*3]-b->field[ib*3];
        double dy=(double)a->field[ia*3+1]-b->field[ib*3+1];
        double dx=(double)a->field[ia*3+2]-b->field[ib*3+2];
        double gap=sqrt(dz*dz+dy*dy+dx*dx);
        n++;sum+=gap;if(gap>mx)mx=gap;if(gap>compatible)bad++;
    }
    if(out_count)*out_count=n;if(out_conflict)*out_conflict=bad;
    if(out_mean)*out_mean=n?sum/(double)n:0.0;if(out_max)*out_max=mx;
}

/* Minimum separation of two half-open developed-grid rectangles.  Zero means
 * overlap or direct edge contact.  A large value is decisive evidence that a
 * 3-D near-contact is a fold/wrap encounter, not a seam the assembler may loft. */
static int loaded_developed_gap(const Loaded *a,const Loaded *b,
                                int *out_row_gap,int *out_col_gap) {
    int ar1=a->gr0+a->H,br1=b->gr0+b->H;
    int ac1=a->c0+a->W,bc1=b->c0+b->W;
    int rg=a->gr0>=br1?a->gr0-br1:(b->gr0>=ar1?b->gr0-ar1:0);
    int cg=a->c0>=bc1?a->c0-bc1:(b->c0>=ac1?b->c0-ac1:0);
    if(out_row_gap)*out_row_gap=rg;if(out_col_gap)*out_col_gap=cg;
    return rg>cg?rg:cg;
}

static int solid_overlap_selftest(void) {
    uint8_t ad[6]={1,1,1,1,1,1},bd[6]={1,1,1,1,1,1};
    float af[18]={0},bf[18]={0};
    Loaded a={0},b={0};
    a.domain=ad;a.field=af;a.H=2;a.W=3;a.gr0=0;a.c0=0;
    b.domain=bd;b.field=bf;b.H=2;b.W=3;b.gr0=0;b.c0=1;
    /* The two overlap columns have gaps 0 and 10 in both rows. */
    bf[(size_t)1*3+2]=10.0f;
    bf[((size_t)3+1)*3+2]=10.0f;
    size_t n=0,bad=0;double mean=0.0,mx=0.0;
    loaded_overlap_audit(&a,&b,5.0,&n,&bad,&mean,&mx);
    int ok=n==4 && bad==2 && fabs(mean-5.0)<1e-9 && fabs(mx-10.0)<1e-9;
    fprintf(stderr,"  %s: solid false-merge overlap audit\n",ok?"ok":"FAIL");
    return ok?0:1;
}

static const double SNAP_TXN_P95_STEP = 0.030;
static const double SNAP_TXN_P99_STEP = 0.100;
static const double SNAP_TXN_OVER4_FACTOR = 1.5;
static const double SNAP_TXN_OVER4_STEP = 0.0001;
static const double SNAP_TXN_P95_ENVELOPE = 0.050;
static const double SNAP_TXN_P99_ENVELOPE = 0.150;
static const double SNAP_TXN_OVER4_ENVELOPE = 0.0002;
static const double SNAP_TXN_REFIT_DARK_MIN_FRACTION = 0.05;
static const double SNAP_TXN_REFIT_DARK_P95_STEP = 0.080;
static const double SNAP_TXN_REFIT_DARK_P99_STEP = 0.250;
static const double SNAP_TXN_REFIT_DARK_P95_ENVELOPE = 0.100;
static const double SNAP_TXN_REFIT_DARK_P99_ENVELOPE = 0.300;
static const double SNAP_TXN_FIRST_DARK_P95_STEP = 0.200;
static const double SNAP_TXN_FIRST_DARK_P99_STEP = 0.350;
static const double SNAP_TXN_FIRST_DARK_OVER4_STEP = 0.00025;
static const double SNAP_TXN_MATURE_P95_EPS = 0.005;
static const double SNAP_TXN_MATURE_P99_EPS = 0.005;
static const double SNAP_TXN_MATURE_QUILT_GAIN = 0.025;
static const int SNAP_TXN_MATURE_CT_GAIN = 2;
static const double SNAP_QUILT_SEAM_TRIGGER_RATIO = 1.25;
static const double SNAP_QUILT_SEAM_TRIGGER_STEP = 0.030;
static const size_t SNAP_QUILT_SEAM_MIN_EDGES = 128;

static int snap_intersection_audit(const float *verts,size_t nv,
                                   const int32_t *faces,size_t nf,
                                   IntersectionCleanupStats *out) {
    IntersectionCleanupParams p;
    IntersectionCleanup_default_params(&p);
    p.gap_max=1.0;
    p.parallel_angle_deg=20.0;
    p.max_conflicts=nf>1000000u?nf:1000000u;
    p.include_hinges=0;
    return IntersectionCleanup_audit_visit_parallel(
        verts,nv,faces,nf,NULL,&p,NULL,NULL,NULL,out);
}

static int snap_intersection_audit_degrees(const float *verts,size_t nv,
                                           const int32_t *faces,size_t nf,
                                           size_t *face_conflict_degree,
                                           IntersectionCleanupStats *out) {
    IntersectionCleanupParams p;
    IntersectionCleanup_default_params(&p);
    p.gap_max=1.0;p.parallel_angle_deg=20.0;
    p.max_conflicts=nf>1000000u?nf:1000000u;p.include_hinges=0;
    return IntersectionCleanup_audit_visit_parallel(verts,nv,faces,nf,NULL,&p,
                                     face_conflict_degree,NULL,NULL,out);
}

typedef struct PostSnapArapTxn {
    int attempted,accepted,solver_code,trial_audits,orientation_rejects;
    size_t locally_rolled_back_vertices;
    double selected_scale,displacement_rms,displacement_max,axial_drift_max;
    double input_edge_log_rms,selected_edge_log_rms;
    size_t input_intersections,input_overlap,input_stab,input_fold;
    size_t output_intersections,output_overlap,output_stab,output_fold;
    QuadStripArapStats solver;
} PostSnapArapTxn;

static double post_snap_edge_log_rms(const float *p,const float *uv,
                                     int H,int W) {
    double ss=0.0;size_t count=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++){
        size_t i=(size_t)r*(size_t)W+(size_t)c,js[3];int ne=0;
        if(c+1<W)js[ne++]=i+1;
        if(r+1<H)js[ne++]=i+(size_t)W;
        if(r+1<H&&c+1<W)js[ne++]=i+(size_t)W+1;
        for(int e=0;e<ne;e++){
            size_t j=js[e];double dz=(double)p[j*3]-p[i*3];
            double dy=(double)p[j*3+1]-p[i*3+1];
            double dx=(double)p[j*3+2]-p[i*3+2];
            double du=(double)uv[j*2]-uv[i*2];
            double dv=(double)uv[j*2+1]-uv[i*2+1];
            double got=sqrt(dz*dz+dy*dy+dx*dx),rest=hypot(du,dv);
            if(got>0.0&&rest>0.0&&isfinite(got)){
                double q=log(got/rest);ss+=q*q;count++;
            }
        }
    }
    return count?sqrt(ss/(double)count):INFINITY;
}

static int post_snap_face_orientation_bad(const float *p,const float *anchor,
                                          size_t nv,const int32_t *face) {
        int32_t ii[3]={face[0],face[1],face[2]};
        if(ii[0]<0||ii[1]<0||ii[2]<0||(size_t)ii[0]>=nv||
           (size_t)ii[1]>=nv||(size_t)ii[2]>=nv)return 1;
        size_t a=(size_t)ii[0],b=(size_t)ii[1],c=(size_t)ii[2];
        double oe0[3],oe1[3],ce0[3],ce1[3],on[3],cn[3];
        for(int k=0;k<3;k++){
            oe0[k]=(double)anchor[b*3+(size_t)k]-anchor[a*3+(size_t)k];
            oe1[k]=(double)anchor[c*3+(size_t)k]-anchor[a*3+(size_t)k];
            ce0[k]=(double)p[b*3+(size_t)k]-p[a*3+(size_t)k];
            ce1[k]=(double)p[c*3+(size_t)k]-p[a*3+(size_t)k];
        }
        on[0]=oe0[1]*oe1[2]-oe0[2]*oe1[1];
        on[1]=oe0[2]*oe1[0]-oe0[0]*oe1[2];
        on[2]=oe0[0]*oe1[1]-oe0[1]*oe1[0];
        cn[0]=ce0[1]*ce1[2]-ce0[2]*ce1[1];
        cn[1]=ce0[2]*ce1[0]-ce0[0]*ce1[2];
        cn[2]=ce0[0]*ce1[1]-ce0[1]*ce1[0];
        double ref2=on[0]*on[0]+on[1]*on[1]+on[2]*on[2];
        double dot=on[0]*cn[0]+on[1]*cn[1]+on[2]*cn[2];
        return !(ref2>1e-18)||!(dot>0.05*ref2);
}

/* Faces already degenerate or inverted in the transaction/round input can
 * never be repaired by zeroing candidate motion: their identity trial fails
 * the same preflight, every alpha dies on them, and the whole transaction
 * aborts with "no exact-audit-safe response remained" (full-sheet v120).
 * Flag them once per round against the round input; candidate gating then
 * reads "no NEW orientation violations".  Returns the flagged count. */
static size_t post_snap_orientation_baseline(const float *anchor,size_t nv,
                                             const int32_t *faces,size_t nf,
                                             uint8_t *input_bad) {
    size_t bad=0;
    ptrdiff_t f;
#pragma omp parallel for reduction(+:bad) schedule(static)
    for(f=0;f<(ptrdiff_t)nf;f++) {
        uint8_t b=(uint8_t)post_snap_face_orientation_bad(
            anchor,anchor,nv,&faces[(size_t)f*3]);
        input_bad[(size_t)f]=b;bad+=b;
    }
    return bad;
}

static size_t post_snap_orientation_violations(const float *p,const float *anchor,
                                               size_t nv,const int32_t *faces,
                                               size_t nf,const uint8_t *input_bad,
                                               float *motion_scale,
                                               size_t *newly_zeroed) {
    size_t bad=0,nnew=0;
    if(!motion_scale) {
        /* Candidate preflights scan tens of millions of independent faces and
         * often reject several motion scales before one exact BVH audit. */
        ptrdiff_t f;
#pragma omp parallel for reduction(+:bad) schedule(static)
        for(f=0;f<(ptrdiff_t)nf;f++) {
            if(input_bad&&input_bad[(size_t)f])continue;
            bad+=(size_t)post_snap_face_orientation_bad(
                p,anchor,nv,&faces[(size_t)f*3]);
        }
        if(newly_zeroed)*newly_zeroed=0;
        return bad;
    }
    for(size_t f=0;f<nf;f++){
        if(input_bad&&input_bad[f])continue;
        if(post_snap_face_orientation_bad(p,anchor,nv,&faces[f*3])){
            const int32_t *ii=&faces[f*3];bad++;
            for(int k=0;k<3;k++){
                size_t v=(size_t)ii[k];
                if(motion_scale[v]>0.0f){motion_scale[v]=0.0f;nnew++;}
            }
        }
    }
    if(newly_zeroed)*newly_zeroed=nnew;
    return bad;
}

/* Face-list preflight for band-confined candidates: faces outside the list
 * are bitwise identical to the round input, so with the baseline mask their
 * verdict cannot change and the restricted scan equals the full scan. */
static size_t post_snap_orientation_violations_list(
        const float *p,const float *anchor,size_t nv,const int32_t *faces,
        const int32_t *face_list,size_t nfl,const uint8_t *input_bad,
        float *motion_scale,size_t *newly_zeroed) {
    size_t bad=0,nnew=0;
    if(!motion_scale) {
        ptrdiff_t i;
#pragma omp parallel for reduction(+:bad) schedule(static)
        for(i=0;i<(ptrdiff_t)nfl;i++) {
            size_t f=(size_t)face_list[i];
            if(input_bad&&input_bad[f])continue;
            bad+=(size_t)post_snap_face_orientation_bad(
                p,anchor,nv,&faces[f*3]);
        }
        if(newly_zeroed)*newly_zeroed=0;
        return bad;
    }
    for(size_t i=0;i<nfl;i++){
        size_t f=(size_t)face_list[i];
        if(input_bad&&input_bad[f])continue;
        if(post_snap_face_orientation_bad(p,anchor,nv,&faces[f*3])){
            const int32_t *ii=&faces[f*3];bad++;
            for(int k=0;k<3;k++){
                size_t v=(size_t)ii[k];
                if(motion_scale[v]>0.0f){motion_scale[v]=0.0f;nnew++;}
            }
        }
    }
    if(newly_zeroed)*newly_zeroed=nnew;
    return bad;
}

static void feather_motion_scale(float *scale,float *temporary,size_t nv,
                                 int H,int W,int rings) {
    if(!scale||!temporary||H<1||W<1||rings<1||
       nv!=(size_t)H*(size_t)W)return;
    for(int ring=0;ring<rings;ring++) {
        memcpy(temporary,scale,nv*sizeof *temporary);
        int r;
        /* Each output cell reads only the immutable scale from the preceding
         * ring, so rows are independent and deterministic. */
#pragma omp parallel for schedule(static)
        for(r=0;r<H;r++)for(int c=0;c<W;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            float a=scale[v],b;
            if(c>0){b=scale[v-1]+0.125f;if(b<a)a=b;}
            if(c+1<W){b=scale[v+1]+0.125f;if(b<a)a=b;}
            if(r>0){b=scale[v-(size_t)W]+0.125f;if(b<a)a=b;}
            if(r+1<H){b=scale[v+(size_t)W]+0.125f;if(b<a)a=b;}
            temporary[v]=a<1.0f?a:1.0f;
        }
        memcpy(scale,temporary,nv*sizeof *scale);
    }
}

/* A remaining scroll pass-through has a very specific signature: two exact
 * triangle hits whose material-U centroids are about one winding apart.  Do
 * not turn that into a global cylinder constraint.  The measured portions of
 * the sheet are already good; collect only those nonlocal pairs and repair the
 * small generated patches which participate in them. */
typedef struct TurnOrderConflict {
    size_t face_lo,face_hi;
    double u_lo,u_hi;
    int hit_kind;
} TurnOrderConflict;

typedef struct TurnOrderConflictScan {
    const float *uv;
    const int32_t *faces;
    double minimum_u_separation;
    size_t *local_face_degree;
    TurnOrderConflict *pair;
    size_t count,capacity;
    int allocation_failed;
} TurnOrderConflictScan;

typedef struct TurnOrderRepairStats {
    int attempted,accepted_rounds;
    size_t input_conflicts,input_long_conflicts;
    size_t output_conflicts,output_long_conflicts;
    size_t seed_vertices,active_vertices;
    double pitch_median,desired_gap,movement_rms,movement_max;
} TurnOrderRepairStats;

static double turn_face_mean_u(const float *uv,const int32_t *face) {
    return ((double)uv[(size_t)face[0]*2]+(double)uv[(size_t)face[1]*2]+
            (double)uv[(size_t)face[2]*2])/3.0;
}

static double turn_face_mean_radius(const float *verts,const int32_t *face,
                                    double axis_y,double axis_x) {
    double sum=0.0;
    for(int k=0;k<3;k++) {
        size_t v=(size_t)face[k];
        sum+=hypot((double)verts[v*3+1]-axis_y,
                   (double)verts[v*3+2]-axis_x);
    }
    return sum/3.0;
}

static int turn_order_conflict_visitor(size_t face_a,size_t face_b,
                                       int hit_kind,void *context) {
    TurnOrderConflictScan *s=(TurnOrderConflictScan *)context;
    const int32_t *fa=&s->faces[face_a*3],*fb=&s->faces[face_b*3];
    double ua=turn_face_mean_u(s->uv,fa),ub=turn_face_mean_u(s->uv,fb);
    (void)hit_kind;
    if(fabs(ub-ua)<s->minimum_u_separation) {
        if(s->local_face_degree) {
            s->local_face_degree[face_a]++;
            s->local_face_degree[face_b]++;
        }
        return 0;
    }
    if(s->count==s->capacity) {
        size_t cap=s->capacity?s->capacity*2u:1024u;
        if(cap<s->capacity||cap>SIZE_MAX/sizeof *s->pair) {
            s->allocation_failed=1;return -1;
        }
        TurnOrderConflict *p=(TurnOrderConflict *)realloc(
            s->pair,cap*sizeof *s->pair);
        if(!p){s->allocation_failed=1;return -1;}
        s->pair=p;s->capacity=cap;
    }
    TurnOrderConflict *p=&s->pair[s->count++];
    if(ua<=ub){p->face_lo=face_a;p->face_hi=face_b;p->u_lo=ua;p->u_hi=ub;}
    else{p->face_lo=face_b;p->face_hi=face_a;p->u_lo=ub;p->u_hi=ua;}
    p->hit_kind=hit_kind;
    return 0;
}

/* Wall-clock accounting for the exact audits: they dominate every transaction
 * on large sheets, so each call is logged and per-transaction totals are
 * reported by the drivers (file-private instrumentation only). */
static size_t audit_call_count=0;
static double audit_total_ms=0.0;

static void audit_totals_reset(void){audit_call_count=0;audit_total_ms=0.0;}

static void audit_totals_report(const char *tag) {
    fprintf(stderr,"[audit] %s totals: calls=%zu wall=%.0f ms (mean %.0f ms)\n",
            tag,audit_call_count,audit_total_ms,
            audit_call_count?audit_total_ms/(double)audit_call_count:0.0);
}

static int turn_order_audit(const float *verts,size_t nv,
                            const int32_t *faces,size_t nf,const float *uv,
                            double minimum_u_separation,
                            size_t *face_degree,size_t *local_face_degree,
                            TurnOrderConflictScan *scan,
                            IntersectionCleanupStats *stats) {
    IntersectionCleanupParams p;
    memset(scan,0,sizeof *scan);
    scan->uv=uv;scan->faces=faces;
    scan->minimum_u_separation=minimum_u_separation;
    scan->local_face_degree=local_face_degree;
    if(local_face_degree)memset(local_face_degree,0,nf*sizeof *local_face_degree);
    IntersectionCleanup_default_params(&p);
    p.gap_max=1.0;p.parallel_angle_deg=20.0;
    /* A fixed million-pair ceiling rejects the 10x10x10 initializer before it
     * can report or repair anything.  The audit visitor streams pairs, and the
     * retained nonlocal ledger is already allocation-checked, so use a linear
     * input-sized safety bound.  This remains fail-closed while scaling with a
     * legitimately larger sheet instead of with an arbitrary old fixture. */
    p.max_conflicts=nf>1000000u?nf:1000000u;p.include_hinges=0;
    clock_t audit_t0=clock();
    int audit_rc=IntersectionCleanup_audit_visit_parallel(
        verts,nv,faces,nf,NULL,&p,face_degree,turn_order_conflict_visitor,
        scan,stats);
    {
        double ms=1000.0*(double)(clock()-audit_t0)/(double)CLOCKS_PER_SEC;
        audit_call_count++;audit_total_ms+=ms;
        fprintf(stderr,"[audit] faces=%zu conflicts=%zu candidates=%zu ms=%.0f\n",
                nf,stats?stats->conflicts:0,
                stats?stats->candidate_pairs:0,ms);
    }
    if(audit_rc!=0||scan->allocation_failed) {
        fprintf(stderr,
                "[turn-order] exact audit incomplete: rc=%d conflicts=%zu "
                "cap=%zu candidates=%zu retained-long=%zu allocation-failed=%d\n",
                audit_rc,stats?stats->conflicts:0,p.max_conflicts,
                stats?stats->candidate_pairs:0,scan->count,
                scan->allocation_failed);
        free(scan->pair);memset(scan,0,sizeof *scan);return -1;
    }
    return 0;
}

static int turn_double_cmp(const void *aa,const void *bb) {
    double a=*(const double *)aa,b=*(const double *)bb;
    return a<b?-1:a>b?1:0;
}

typedef struct TurnPhasePoint {
    double phase;
    int column;
    size_t vertex;
} TurnPhasePoint;

static int turn_phase_point_cmp(const void *aa,const void *bb) {
    const TurnPhasePoint *a=(const TurnPhasePoint *)aa;
    const TurnPhasePoint *b=(const TurnPhasePoint *)bb;
    return a->phase<b->phase?-1:a->phase>b->phase?1:0;
}

static size_t turn_phase_lower_bound(const TurnPhasePoint *p,size_t n,double x) {
    size_t lo=0,hi=n;
    while(lo<hi){size_t mid=lo+(hi-lo)/2;
        if(p[mid].phase<x)lo=mid+1;else hi=mid;}
    return lo;
}

/* Infer both winding direction and pitch from pairs of vertices on the same
 * row and at the UV period actually observed in the conflict graph.  The first
 * pass uses hard-certified samples only.  A fallback to the whole initializer
 * is permitted solely when there are too few certified pairs to establish a
 * sign; it is logged as such. */
static int turn_order_infer_pitch(const float *verts,const float *phase,
                                  const uint8_t *movable,int H,int W,
                                  double grid_du,double axis_y,double axis_x,
                                  const TurnOrderConflictScan *scan,
                                  int *out_dcol,double *out_pitch) {
    if(!scan||scan->count==0||H<1||W<3||!(grid_du>0.0))return -1;
    const double two_pi=6.2831853071795864769;
    size_t capacity=(size_t)H*(size_t)W;
    if(scan->count>capacity)capacity=scan->count;
    double *span=(double *)malloc(capacity*sizeof *span);
    double *sample=(double *)malloc(capacity*sizeof *sample);
    TurnPhasePoint *row=(TurnPhasePoint *)malloc((size_t)W*sizeof *row);
    if(!span||!sample||!row){free(span);free(sample);free(row);return -1;}
    size_t ns=0;int fallback=0;
    if(phase)for(int r=0;r<H;r++) {
        size_t nr=0;
        for(int c=0;c<W;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            if((movable&&movable[v])||!isfinite((double)phase[v]))continue;
            row[nr].phase=(double)phase[v];row[nr].column=c;
            row[nr].vertex=v;nr++;
        }
        qsort(row,nr,sizeof *row,turn_phase_point_cmp);
        for(size_t i=0;i<nr;i++) {
            double target=row[i].phase+two_pi,error=1e300;size_t best=nr;
            size_t j=turn_phase_lower_bound(row,nr,target);
            if(j<nr){best=j;error=fabs(row[j].phase-target);}
            if(j>0&&fabs(row[j-1].phase-target)<error) {
                best=j-1;error=fabs(row[j-1].phase-target);
            }
            if(best>=nr||best==i||error>0.05)continue;
            size_t a=row[i].vertex,b=row[best].vertex;
            double ra=hypot((double)verts[a*3+1]-axis_y,
                            (double)verts[a*3+2]-axis_x);
            double rb=hypot((double)verts[b*3+1]-axis_y,
                            (double)verts[b*3+2]-axis_x);
            if(!isfinite(ra)||!isfinite(rb))continue;
            int dc=row[best].column-row[i].column;
            sample[ns]=dc>=0?rb-ra:ra-rb;span[ns]=(double)(dc>=0?dc:-dc);
            ns++;
        }
    }
    /* Malformed/absent phase is a diagnostic fallback, never the preferred
     * winding model: use the median conflict period on same-row geometry. */
    if(ns<16) {
        fallback=1;ns=0;
        for(size_t q=0;q<scan->count;q++)
            span[q]=(scan->pair[q].u_hi-scan->pair[q].u_lo)/grid_du;
        qsort(span,scan->count,sizeof *span,turn_double_cmp);
        int dc=(int)lround(span[scan->count/2]);
        if(dc<2)dc=2;if(dc>=W)dc=W-1;
        for(int r=0;r<H;r++)for(int c=0;c+dc<W;c++) {
            size_t a=(size_t)r*(size_t)W+(size_t)c,b=a+(size_t)dc;
            if(movable&&(movable[a]||movable[b]))continue;
            double ra=hypot((double)verts[a*3+1]-axis_y,
                            (double)verts[a*3+2]-axis_x);
            double rb=hypot((double)verts[b*3+1]-axis_y,
                            (double)verts[b*3+2]-axis_x);
            if(isfinite(ra)&&isfinite(rb)){sample[ns]=rb-ra;span[ns]=(double)dc;ns++;}
        }
    }
    if(ns<16){free(span);free(sample);free(row);return -1;}
    qsort(span,ns,sizeof *span,turn_double_cmp);
    int dcol=(int)lround(span[ns/2]);
    if(dcol<2)dcol=2;if(dcol>=W)dcol=W-1;
    qsort(sample,ns,sizeof *sample,turn_double_cmp);
    double pitch=sample[ns/2];
    size_t agrees=0;
    for(size_t k=0;k<ns;k++)agrees+=(sample[k]*pitch)>0.0;
    double agreement=(double)agrees/(double)ns;
    fprintf(stderr,
        "[turn-order] inferred period=%d columns (%.1f UV); pitch "
        "p10/p50/p90=%.3f/%.3f/%.3f vox, sign agreement %.1f%%, "
        "samples=%zu%s\n",
        dcol,(double)dcol*grid_du,sample[ns/10],pitch,
        sample[(ns*9)/10],100.0*agreement,ns,
        fallback?" (initializer fallback)":" (hard-certified)");
    free(span);free(sample);free(row);
    if(!isfinite(pitch)||fabs(pitch)<2.0||agreement<0.60)return -1;
    *out_dcol=dcol;*out_pitch=pitch;return 0;
}

static size_t turn_mark_face_seed(const int32_t *face,const uint8_t *movable,
                                  double value,double weight,
                                  double *seed_sum,double *seed_weight) {
    size_t n=0;
    for(int k=0;k<3;k++)n+=(size_t)(movable[(size_t)face[k]]!=0);
    if(!n)return 0;
    for(int k=0;k<3;k++) {
        size_t v=(size_t)face[k];if(!movable[v])continue;
        seed_sum[v]+=weight*value;seed_weight[v]+=weight;
    }
    return n;
}

/* Construct a compact screened-harmonic radial correction.  The active set is
 * only a 32-grid-ring halo around exact long-range collision faces, and hard
 * observations are zero-Dirichlet holes in that set. */
static int turn_order_build_correction(
        const float *verts,size_t nv,const int32_t *faces,
        const uint8_t *movable,int H,int W,double axis_y,double axis_x,
        const TurnOrderConflictScan *scan,double pitch,double desired_gap,
        float *correction,size_t *out_seed,size_t *out_active,
        int *out_iterations,double *out_max) {
    double *seed_sum=(double *)calloc(nv,sizeof *seed_sum);
    double *seed_weight=(double *)calloc(nv,sizeof *seed_weight);
    uint8_t *active=(uint8_t *)calloc(nv,1),*tmp=(uint8_t *)calloc(nv,1);
    if(!seed_sum||!seed_weight||!active||!tmp) {
        free(seed_sum);free(seed_weight);free(active);free(tmp);return -1;
    }
    const double order=pitch>=0.0?1.0:-1.0;
    const double deficit_cap=fmax(desired_gap,1.5*fabs(pitch));
    for(size_t q=0;q<scan->count;q++) {
        const int32_t *lo=&faces[scan->pair[q].face_lo*3];
        const int32_t *hi=&faces[scan->pair[q].face_hi*3];
        double rlo=turn_face_mean_radius(verts,lo,axis_y,axis_x);
        double rhi=turn_face_mean_radius(verts,hi,axis_y,axis_x);
        double deficit=desired_gap-order*(rhi-rlo);
        if(!(deficit>0.0))continue;
        if(deficit>deficit_cap)deficit=deficit_cap;
        size_t nlo=0,nhi=0;
        for(int k=0;k<3;k++){
            nlo+=(size_t)(movable[(size_t)lo[k]]!=0);
            nhi+=(size_t)(movable[(size_t)hi[k]]!=0);
        }
        if(!nlo&&!nhi)continue;
        double lo_share=nlo?(nhi?0.5:1.0):0.0;
        double hi_share=nhi?(nlo?0.5:1.0):0.0;
        double severity=1.0+deficit/fmax(desired_gap,1.0);
        if(lo_share)turn_mark_face_seed(lo,movable,-order*deficit*lo_share,
                                        severity,seed_sum,seed_weight);
        if(hi_share)turn_mark_face_seed(hi,movable, order*deficit*hi_share,
                                        severity,seed_sum,seed_weight);
    }
    size_t seeds=0;
    for(size_t v=0;v<nv;v++)if(seed_weight[v]>0.0){active[v]=1;seeds++;}
    if(!seeds){free(seed_sum);free(seed_weight);free(active);free(tmp);return 1;}
    for(int ring=0;ring<32;ring++) {
        memcpy(tmp,active,nv);
        for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
            size_t v=(size_t)r*(size_t)W+(size_t)c;
            if(!active[v])continue;
            size_t n[4];int nn=0;
            if(c>0)n[nn++]=v-1;if(c+1<W)n[nn++]=v+1;
            if(r>0)n[nn++]=v-(size_t)W;
            if(r+1<H)n[nn++]=v+(size_t)W;
            for(int k=0;k<nn;k++)if(movable[n[k]])tmp[n[k]]=1;
        }
        memcpy(active,tmp,nv);
    }
    size_t nactive=0;for(size_t v=0;v<nv;v++)nactive+=active[v]!=0;
    memset(correction,0,nv*sizeof *correction);
    const double screen=0.002;
    int iterations=0;
    for(iterations=0;iterations<384;iterations++) {
        double maximum_change=0.0;
        for(int parity=0;parity<2;parity++) {
#pragma omp parallel
            {
                double thread_maximum_change=0.0;
                int r;
#pragma omp for schedule(static)
                for(r=0;r<H;r++)for(int c=0;c<W;c++) {
                if(((r+c)&1)!=parity)continue;
                size_t v=(size_t)r*(size_t)W+(size_t)c;
                if(!active[v]||!movable[v])continue;
                double sum=0.0,degree=0.0;
                if(c>0){sum+=(double)correction[v-1];degree+=1.0;}
                if(c+1<W){sum+=(double)correction[v+1];degree+=1.0;}
                if(r>0){sum+=(double)correction[v-(size_t)W];degree+=1.0;}
                if(r+1<H){sum+=(double)correction[v+(size_t)W];degree+=1.0;}
                double sw=seed_weight[v]>0.0
                         ?20.0*fmin(seed_weight[v],5.0):0.0;
                double target=seed_weight[v]>0.0
                             ?seed_sum[v]/seed_weight[v]:0.0;
                double value=(sum+sw*target)/(degree+screen+sw);
                double change=fabs(value-(double)correction[v]);
                if(change>thread_maximum_change)thread_maximum_change=change;
                correction[v]=(float)value;
                }
#pragma omp critical(turn_order_maximum_change)
                if(thread_maximum_change>maximum_change)
                    maximum_change=thread_maximum_change;
            }
        }
        if(maximum_change<=1e-5){iterations++;break;}
    }
    double mx=0.0;
    for(size_t v=0;v<nv;v++)if(fabs((double)correction[v])>mx)
        mx=fabs((double)correction[v]);
    *out_seed=seeds;*out_active=nactive;*out_iterations=iterations;*out_max=mx;
    free(seed_sum);free(seed_weight);free(active);free(tmp);return 0;
}

static void turn_order_apply_radial(float *out,const float *base,size_t nv,
                                    const uint8_t *movable,
                                    const float *correction,const float *keep,
                                    double scale,
                                    double axis_y,double axis_x,
                                    double *out_rms,double *out_max) {
    double ss=0.0,mx=0.0;size_t moved=0;
    memcpy(out,base,nv*3*sizeof *out);
    for(size_t v=0;v<nv;v++) {
        if(!movable[v]||correction[v]==0.0f)continue;
        double dy=(double)base[v*3+1]-axis_y;
        double dx=(double)base[v*3+2]-axis_x;
        double retained=keep?(double)keep[v]:1.0;
        double radius=hypot(dy,dx),d=scale*retained*(double)correction[v];
        if(!(radius>1e-9)||radius+d<=1e-6)continue;
        double factor=(radius+d)/radius;
        out[v*3+1]=(float)(axis_y+factor*dy);
        out[v*3+2]=(float)(axis_x+factor*dx);
        ss+=d*d;if(fabs(d)>mx)mx=fabs(d);moved++;
    }
    *out_rms=moved?sqrt(ss/(double)moved):0.0;*out_max=mx;
}

static int repair_cylindrical_turn_order(
        float *verts,size_t nv,const int32_t *faces,size_t nf,const float *uv,
        const float *phase,const uint8_t *movable,int H,int W,double grid_du,
        double axis_y,double axis_x,double minimum_u_separation,
        TurnOrderRepairStats *out) {
    TurnOrderRepairStats s;memset(&s,0,sizeof s);s.attempted=1;
    if(out)*out=s;
    if(!verts||!faces||!uv||!movable||H<2||W<3||
       nv!=(size_t)H*(size_t)W)return -1;
    float *base=(float *)malloc(nv*3*sizeof *base);
    float *candidate=(float *)malloc(nv*3*sizeof *candidate);
    float *correction=(float *)malloc(nv*sizeof *correction);
    float *keep=(float *)malloc(nv*sizeof *keep);
    float *keep_tmp=(float *)malloc(nv*sizeof *keep_tmp);
    size_t *degree_before=(size_t *)calloc(nf,sizeof *degree_before);
    size_t *degree_after=(size_t *)calloc(nf,sizeof *degree_after);
    size_t *local_before=(size_t *)calloc(nf,sizeof *local_before);
    size_t *local_after=(size_t *)calloc(nf,sizeof *local_after);
    uint8_t *orientation_input_bad=(uint8_t *)malloc(nf);
    if(!base||!candidate||!correction||!keep||!keep_tmp||!degree_before||
       !degree_after||!local_before||!local_after||!orientation_input_bad){
        free(base);free(candidate);free(correction);free(keep);free(keep_tmp);
        free(degree_before);free(degree_after);free(local_before);free(local_after);
        free(orientation_input_bad);
        return -1;
    }
    audit_totals_reset();
    TurnOrderConflictScan cached_scan={0};
    IntersectionCleanupStats cached_stats={0};
    int cached_audit=0;
    for(int round=0;round<10;round++) {
        TurnOrderConflictScan before_scan={0};
        IntersectionCleanupStats before={0};
        if(cached_audit) {
            before_scan=cached_scan;before=cached_stats;
            memset(&cached_scan,0,sizeof cached_scan);cached_audit=0;
            fprintf(stderr,
                "[turn-order] reusing accepted exact audit for next solve\n");
        } else if(turn_order_audit(
                    verts,nv,faces,nf,uv,minimum_u_separation,
                    degree_before,local_before,&before_scan,&before)!=0) {
            free(base);free(candidate);free(correction);free(keep);free(keep_tmp);
            free(degree_before);free(degree_after);free(local_before);free(local_after);
            free(orientation_input_bad);
            return -1;
        }
        if(round==0){s.input_conflicts=before.conflicts;
                    s.input_long_conflicts=before_scan.count;}
        fprintf(stderr,
            "[turn-order] round %d exact conflicts=%zu (long-range=%zu), "
            "folds=%zu\n",round+1,before.conflicts,before_scan.count,
            before.fold_pairs);
        if(before_scan.count==0){
            cached_scan=before_scan;cached_stats=before;cached_audit=1;
            break;
        }
        int dcol=0;double pitch=0.0;
        if(turn_order_infer_pitch(verts,phase,movable,H,W,grid_du,axis_y,axis_x,
                                  &before_scan,&dcol,&pitch)!=0) {
            fprintf(stderr,"[turn-order] refusing ambiguous winding direction/pitch\n");
            cached_scan=before_scan;cached_stats=before;cached_audit=1;break;
        }
        (void)dcol;
        double desired_gap=fmax(3.0,fmin(16.0,0.75*fabs(pitch)));
        size_t seeds=0,active=0;int solve_iterations=0;double correction_max=0.0;
        int crc=turn_order_build_correction(
            verts,nv,faces,movable,H,W,axis_y,axis_x,&before_scan,pitch,
            desired_gap,correction,&seeds,&active,&solve_iterations,&correction_max);
        if(crc!=0) {
            if(crc>0)fprintf(stderr,"[turn-order] no movable vertices in long-range pairs\n");
            if(crc<0){
                free(before_scan.pair);
                free(base);free(candidate);free(correction);free(keep);free(keep_tmp);
                free(degree_before);free(degree_after);free(local_before);free(local_after);
                free(orientation_input_bad);
                return -1;
            }
            cached_scan=before_scan;cached_stats=before;cached_audit=1;
            break;
        }
        fprintf(stderr,
            "[turn-order] sparse radial solve: desired gap %.3f, seeds=%zu, "
            "active=%zu/%zu, iterations=%d, max correction %.3f vox\n",
            desired_gap,seeds,active,nv,solve_iterations,correction_max);
        memcpy(base,verts,nv*3*sizeof *base);
        {
            size_t input_bad_faces=post_snap_orientation_baseline(
                base,nv,faces,nf,orientation_input_bad);
            if(input_bad_faces)fprintf(stderr,
                "[turn-order] input orientation-bad faces=%zu excluded from "
                "candidate preflight\n",input_bad_faces);
        }
        /* Dense intersection fronts often need a sub-quarter transaction once
         * the large easy crossings have left the sheet.  Try those scales
         * before escalating; a rejected trial is read-only and exact-audited. */
        const double trial_scale_dense[7]={0.25,0.125,0.0625,0.375,0.5,0.75,1.0};
        const double trial_scale_sparse[7]={0.375,0.5,0.25,0.125,0.0625,0.75,1.0};
        const double *trial_scale=before_scan.count>3000
                                  ?trial_scale_dense:trial_scale_sparse;
        /* Large sheets can leave a useful globally improving correction with a
         * thin orientation/front collar after five feathered retries.  Keep
         * shrinking only that collar; every retry still has to improve exact
         * total/long/local counts, preserve folds, and pass the orientation
         * audit before it can be committed. */
        const int maximum_local_rollbacks=12;
        int accepted=0;
        for(int trial=0;trial<7;trial++) {
            for(size_t v=0;v<nv;v++)keep[v]=1.0f;
            for(int rollback=0;rollback<=maximum_local_rollbacks;rollback++) {
                double move_rms=0.0,move_max=0.0;
                turn_order_apply_radial(candidate,base,nv,movable,correction,keep,
                                        trial_scale[trial],axis_y,axis_x,
                                        &move_rms,&move_max);
                size_t orientation_bad=post_snap_orientation_violations(
                    candidate,base,nv,faces,nf,orientation_input_bad,NULL,NULL);
                /* A candidate with an orientation reversal cannot possibly
                 * commit.  On the 36-million-face sheet the exact BVH audit is
                 * orders of magnitude dearer than this linear preflight, so
                 * first freeze and feather the known-bad faces.  A candidate
                 * still receives the complete intersection audit before every
                 * acceptance; this only removes audits whose answer cannot
                 * affect the transaction. */
                if(orientation_bad) {
                    fprintf(stderr,
                        "[turn-order] reject trial alpha %.3f rollback %d "
                        "before exact audit: orientation violations=%zu\n",
                        trial_scale[trial],rollback,orientation_bad);
                    if(rollback==maximum_local_rollbacks)break;
                    size_t orientation_marked=0;
                    (void)post_snap_orientation_violations(
                        candidate,base,nv,faces,nf,orientation_input_bad,
                        keep,&orientation_marked);
                    if(orientation_marked==0)break;
                    feather_motion_scale(keep,keep_tmp,nv,H,W,8);
                    fprintf(stderr,
                        "[turn-order] preflight rolled back %zu orientation "
                        "vertices and feathered 8 grid rings\n",
                        orientation_marked);
                    continue;
                }
                TurnOrderConflictScan after_scan={0};
                IntersectionCleanupStats after={0};
                if(turn_order_audit(candidate,nv,faces,nf,uv,minimum_u_separation,
                                    degree_after,local_after,&after_scan,&after)!=0) {
                    free(before_scan.pair);free(base);free(candidate);free(correction);
                    free(keep);free(keep_tmp);free(degree_before);free(degree_after);
                    free(local_before);free(local_after);
                    free(orientation_input_bad);return -1;
                }
                size_t before_local=before.conflicts-before_scan.count;
                size_t after_local=after.conflicts-after_scan.count;
                int safe=after.conflicts<before.conflicts&&
                         after_scan.count<before_scan.count&&
                         after_local<=before_local&&
                         after.fold_pairs<=before.fold_pairs&&orientation_bad==0;
                fprintf(stderr,
                    "[turn-order] %s trial alpha %.3f rollback %d: move rms/max "
                    "%.3f/%.3f; conflicts %zu->%zu, long %zu->%zu, local "
                    "%zu->%zu, folds %zu->%zu, orientation violations=%zu\n",
                    safe?"accept":"reject",trial_scale[trial],rollback,
                    move_rms,move_max,before.conflicts,after.conflicts,
                    before_scan.count,after_scan.count,before_local,after_local,
                    before.fold_pairs,after.fold_pairs,orientation_bad);
                if(safe) {
                    memcpy(verts,candidate,nv*3*sizeof *verts);accepted=1;
                    s.accepted_rounds++;s.pitch_median=pitch;s.desired_gap=desired_gap;
                    s.seed_vertices=seeds;s.active_vertices=active;
                    s.movement_rms=move_rms;s.movement_max=move_max;
                    cached_scan=after_scan;cached_stats=after;cached_audit=1;
                    {
                        size_t *tmp=degree_before;
                        degree_before=degree_after;degree_after=tmp;
                        tmp=local_before;
                        local_before=local_after;local_after=tmp;
                    }
                    break;
                }
                if(rollback==maximum_local_rollbacks||
                   after.conflicts>=before.conflicts||
                   after_scan.count>=before_scan.count) {
                    free(after_scan.pair);break;
                }
                size_t marked=0,orientation_marked=0;
                for(size_t f=0;f<nf;f++)if(local_after[f]>local_before[f])
                    for(int k=0;k<3;k++) {
                        size_t v=(size_t)faces[f*3+(size_t)k];
                        if(keep[v]>0.0f){keep[v]=0.0f;marked++;}
                    }
                (void)post_snap_orientation_violations(
                    candidate,base,nv,faces,nf,orientation_input_bad,
                    keep,&orientation_marked);
                marked+=orientation_marked;
                free(after_scan.pair);
                if(marked==0)break;
                feather_motion_scale(keep,keep_tmp,nv,H,W,8);
                fprintf(stderr,
                    "[turn-order] locally rolled back %zu seed/boundary vertices "
                    "and feathered 8 grid rings\n",marked);
            }
            if(accepted)break;
        }
        if(accepted)free(before_scan.pair);
        if(!accepted){
            cached_scan=before_scan;cached_stats=before;cached_audit=1;
            break;
        }
    }
    {
        TurnOrderConflictScan final_scan={0};IntersectionCleanupStats final={0};
        if(cached_audit) {
            final_scan=cached_scan;final=cached_stats;
            memset(&cached_scan,0,sizeof cached_scan);cached_audit=0;
        } else if(turn_order_audit(
                    verts,nv,faces,nf,uv,minimum_u_separation,
                    NULL,NULL,&final_scan,&final)!=0) {
            free(base);free(candidate);free(correction);free(keep);free(keep_tmp);
            free(degree_before);free(degree_after);free(local_before);free(local_after);
            free(orientation_input_bad);
            return -1;
        }
        s.output_conflicts=final.conflicts;s.output_long_conflicts=final_scan.count;
        fprintf(stderr,
            "[turn-order] transaction result: rounds=%d, exact conflicts %zu->%zu, "
            "long-range %zu->%zu, folds=%zu\n",
            s.accepted_rounds,s.input_conflicts,s.output_conflicts,
            s.input_long_conflicts,s.output_long_conflicts,final.fold_pairs);
        audit_totals_report("turn-order");
        free(final_scan.pair);
    }
    free(base);free(candidate);free(correction);free(keep);free(keep_tmp);
    free(degree_before);free(degree_after);free(local_before);free(local_after);
    free(orientation_input_bad);
    if(out)*out=s;return 0;
}

/* Exact self-contact response for the structured ribbon.  This is a static
 * elastic-shell solve, not a dynamics simulation.  The unknown is one scalar
 * radial displacement per vertex measured from the pre-escape rest ribbon
 * (the measured coil); Z is fixed by lattice V, so all motion is radial Y/X.
 *
 * Each exact BVH audit's long-range triangle conflicts become speculative
 * quad-cell contact pairs (a +/-3-cell UV inflation kept within a 3-D AABB
 * collar), merged into a persistent ledger so a constraint survives after the
 * intersection curve migrates off it.  Inside a 24-ring active band the solve
 * alternates screened-Laplacian red/black smoothing of the displacement field
 * with active-support projection of the single per-pair inequality
 *   min(order*r_hi) >= max(order*r_lo) + clearance,
 * then rebuilds the candidate FROM the rest ribbon -- open the penetration,
 * settle back toward the measured coil.  A trust-region alpha ladder with
 * feathered local rollbacks gates every candidate through a baseline-relative
 * orientation preflight and a complete exact BVH audit; acceptance requires
 * phase-violation progress under a bounded pair-count churn (the escape is
 * deliberately non-monotone in pair counts) and never admits a new fold or
 * orientation reversal.  An accepted audit is carried directly into the next
 * active-set round, which avoids rebuilding the unchanged accepted state
 * while retaining the exact gate the one-shot harmonic correction lacked.
 * The lexicographically best audited state is what a partial solve publishes;
 * --require-collision-free instead rolls the whole transaction back. */
typedef struct ClothTurnContact {
    size_t cell_lo,cell_hi;
} ClothTurnContact;

typedef struct ClothCollisionStats {
    int attempted,accepted_rounds,elastic_sweeps,contact_sweeps;
    int complete,retained_partial;
    size_t input_conflicts,input_long_conflicts;
    size_t output_conflicts,output_long_conflicts;
    size_t contacts_built,active_vertices,hard_contacts;
    double pitch,clearance,movement_rms,movement_max;
    int settle_rounds_run,settle_accepted;
    double settle_rms_before,settle_rms_after;
    double settle_max_before,settle_max_after;
} ClothCollisionStats;

/* Front-stall escalation: when accepted rounds stop shrinking the long-range
 * set, the intersection curve is travelling slower than the speculative
 * constraints extend.  Widen the collar/inflation/band so constraints stay
 * ahead of the front; drop back to base after a decisive improvement.
 * Returns 1 when the escalation state changed. */
static int cloth_stall_escalation_step(size_t prev_long,size_t new_long,
                                       double collar_base,double collar_max,
                                       int *stall_rounds,double *collar,
                                       int *inflate,int *rings) {
    int changed=0;
    double improvement=prev_long
        ?1.0-(double)new_long/(double)prev_long:1.0;
    if(improvement>=0.20) {
        if(*collar!=collar_base||*inflate!=3||*rings!=24)changed=1;
        *collar=collar_base;*inflate=3;*rings=24;*stall_rounds=0;
    } else if(improvement<0.02) {
        (*stall_rounds)++;
        if(*stall_rounds>=3) {
            *stall_rounds=0;
            if(*collar<collar_max) {
                *collar=fmin(*collar+4.0,collar_max);changed=1;
            }
            if(*inflate<5){*inflate=5;changed=1;}
            if(*rings<48){*rings=48;changed=1;}
        }
    } else {
        *stall_rounds=0;
    }
    return changed;
}

typedef struct ClothFacePair {
    size_t cell_lo,cell_hi;
} ClothFacePair;

static int cloth_face_pair_cmp(const void *aa,const void *bb) {
    const ClothFacePair *a=(const ClothFacePair *)aa;
    const ClothFacePair *b=(const ClothFacePair *)bb;
    if(a->cell_lo<b->cell_lo)return -1;if(a->cell_lo>b->cell_lo)return 1;
    if(a->cell_hi<b->cell_hi)return -1;if(a->cell_hi>b->cell_hi)return 1;
    return 0;
}

static void cloth_cell_vertices(size_t cell,int W,size_t vertex[4]) {
    size_t row=cell/(size_t)(W-1),col=cell%(size_t)(W-1);
    vertex[0]=row*(size_t)W+col;
    vertex[1]=vertex[0]+1u;
    vertex[2]=vertex[0]+(size_t)W;
    vertex[3]=vertex[2]+1u;
}

static int cloth_cells_within_collar(const float *verts,int W,
                                     size_t ca,size_t cb,double collar) {
    size_t av[4],bv[4];
    cloth_cell_vertices(ca,W,av);cloth_cell_vertices(cb,W,bv);
    double distance2=0.0;
    for(int k=0;k<3;k++) {
        double alo=DBL_MAX,ahi=-DBL_MAX,blo=DBL_MAX,bhi=-DBL_MAX;
        for(int j=0;j<4;j++) {
            double a=(double)verts[av[j]*3+(size_t)k];
            double b=(double)verts[bv[j]*3+(size_t)k];
            if(a<alo)alo=a;if(a>ahi)ahi=a;
            if(b<blo)blo=b;if(b>bhi)bhi=b;
        }
        double gap=ahi<blo?blo-ahi:bhi<alo?alo-bhi:0.0;
        distance2+=gap*gap;
    }
    return distance2<=collar*collar;
}

/* Inflate every exact contact by a small UV neighbourhood, then retain only
 * face pairs whose 3-D AABBs lie inside the speculative collar.  Collision
 * engines do this before response: without the collar an intersection can
 * simply migrate from the constrained face pair to its unconstrained neighbour.
 * The two sides use the same row offset because structured V is exact world Z. */
static ClothTurnContact *cloth_build_contacts(
        const float *verts,const TurnOrderConflictScan *scan,int H,int W,
        double collar,int inflate,
    size_t *out_count) {
    if(out_count)*out_count=0;
    if(!scan||scan->count==0||H<2||W<2||inflate<1)return NULL;
    size_t reserve_multiplier=scan->count<1024?64u:8u;
    if(scan->count>SIZE_MAX/reserve_multiplier)return NULL;
    size_t cap=scan->count*reserve_multiplier;
    if(cap<1024)cap=1024;
    if(cap>SIZE_MAX/sizeof(ClothFacePair))return NULL;
    ClothFacePair *pair=(ClothFacePair *)malloc(cap*sizeof *pair);
    if(!pair)return NULL;
    size_t count=0;const int cells_w=W-1;
    for(size_t q=0;q<scan->count;q++) {
        size_t qlo=scan->pair[q].face_lo/2u,qhi=scan->pair[q].face_hi/2u;
        int rlo=(int)(qlo/(size_t)cells_w),clo=(int)(qlo%(size_t)cells_w);
        int rhi=(int)(qhi/(size_t)cells_w),chi=(int)(qhi%(size_t)cells_w);
        for(int dr=-inflate;dr<=inflate;dr++) {
            int rl=rlo+dr,rh=rhi+dr;
            if(rl<0||rl>=H-1||rh<0||rh>=H-1)continue;
            for(int dl=-inflate;dl<=inflate;dl++)
            for(int dh=-inflate;dh<=inflate;dh++) {
                int cl=clo+dl,ch=chi+dh;
                if(cl<0||cl>=cells_w||ch<0||ch>=cells_w)continue;
                size_t cell_lo=(size_t)rl*(size_t)cells_w+(size_t)cl;
                size_t cell_hi=(size_t)rh*(size_t)cells_w+(size_t)ch;
                if(!cloth_cells_within_collar(
                        verts,W,cell_lo,cell_hi,collar))continue;
                if(count==cap) {
                    size_t next=cap*2u;
                    if(next<cap||next>SIZE_MAX/sizeof *pair){free(pair);return NULL;}
                    ClothFacePair *grown=(ClothFacePair *)realloc(pair,next*sizeof *pair);
                    if(!grown){free(pair);return NULL;}pair=grown;cap=next;
                }
                pair[count].cell_lo=cell_lo;pair[count].cell_hi=cell_hi;
                count++;
            }
        }
    }
    if(!count){free(pair);return NULL;}
    qsort(pair,count,sizeof *pair,cloth_face_pair_cmp);
    size_t unique=0;
    for(size_t q=0;q<count;q++) {
        if(unique&&pair[q].cell_lo==pair[unique-1].cell_lo&&
           pair[q].cell_hi==pair[unique-1].cell_hi)continue;
        pair[unique++]=pair[q];
    }
    if(unique>SIZE_MAX/sizeof(ClothTurnContact)){free(pair);return NULL;}
    ClothTurnContact *contact=(ClothTurnContact *)malloc(unique*sizeof *contact);
    if(!contact){free(pair);return NULL;}
    for(size_t q=0;q<unique;q++) {
        contact[q].cell_lo=pair[q].cell_lo;contact[q].cell_hi=pair[q].cell_hi;
    }
    free(pair);if(out_count)*out_count=unique;return contact;
}

static int cloth_merge_contact_ledger(
        const ClothTurnContact *old_contact,size_t nold,
        const ClothTurnContact *new_contact,size_t nnew,
        ClothTurnContact **out_contact,size_t *out_count) {
    if(!out_contact||!out_count||
       (nold&& !old_contact)||(nnew&& !new_contact)||
       nold>SIZE_MAX-nnew||(nold+nnew)>SIZE_MAX/sizeof(ClothTurnContact))return -1;
    ClothTurnContact *merged=(ClothTurnContact *)malloc(
        (nold+nnew?nold+nnew:1u)*sizeof *merged);
    if(!merged)return -1;
    size_t i=0,j=0,n=0;
    while(i<nold||j<nnew) {
        const ClothTurnContact *take=NULL;
        if(j>=nnew)take=&old_contact[i++];
        else if(i>=nold)take=&new_contact[j++];
        else {
            const ClothTurnContact *a=&old_contact[i],*b=&new_contact[j];
            if(a->cell_lo<b->cell_lo||
               (a->cell_lo==b->cell_lo&&a->cell_hi<b->cell_hi))take=&old_contact[i++];
            else if(b->cell_lo<a->cell_lo||
                    (b->cell_lo==a->cell_lo&&b->cell_hi<a->cell_hi))take=&new_contact[j++];
            else{take=&old_contact[i++];j++;}
        }
        if(n&&merged[n-1].cell_lo==take->cell_lo&&
           merged[n-1].cell_hi==take->cell_hi)continue;
        merged[n++]=*take;
    }
    *out_contact=merged;*out_count=n;return 0;
}

/* Frontier BFS produces the identical activation set as the old full-grid
 * ring dilation (activation is idempotent and each vertex expands exactly
 * once, at its first-activation ring) without rings x H x W scans and their
 * two full-grid memcpys per ring. */
static uint8_t *cloth_active_contact_band(
        size_t nv,const ClothTurnContact *contact,
        size_t ncontact,const uint8_t *movable,int H,int W,int rings,
        size_t *out_active) {
    uint8_t *active=(uint8_t *)calloc(nv,1);
    int32_t *frontier=(int32_t *)malloc((nv?nv:1u)*sizeof *frontier);
    int32_t *next=(int32_t *)malloc((nv?nv:1u)*sizeof *next);
    if(!active||!frontier||!next){
        free(active);free(frontier);free(next);return NULL;
    }
    size_t nfrontier=0;
    for(size_t q=0;q<ncontact;q++)for(int side=0;side<2;side++) {
        size_t cell=side?contact[q].cell_hi:contact[q].cell_lo,vertex[4];
        cloth_cell_vertices(cell,W,vertex);
        for(int k=0;k<4;k++)if(vertex[k]<nv&&!active[vertex[k]]) {
            active[vertex[k]]=1;frontier[nfrontier++]=(int32_t)vertex[k];
        }
    }
    size_t count=nfrontier;
    for(int ring=0;ring<rings&&nfrontier;ring++) {
        size_t nnext=0;
        for(size_t i=0;i<nfrontier;i++) {
            size_t v=(size_t)frontier[i];
            int r=(int)(v/(size_t)W),col=(int)(v%(size_t)W);
            size_t n[4];int nn=0;
            if(col>0)n[nn++]=v-1;if(col+1<W)n[nn++]=v+1;
            if(r>0)n[nn++]=v-(size_t)W;
            if(r+1<H)n[nn++]=v+(size_t)W;
            for(int k=0;k<nn;k++)if(movable[n[k]]&&!active[n[k]]) {
                active[n[k]]=1;next[nnext++]=(int32_t)n[k];
            }
        }
        int32_t *swap=frontier;frontier=next;next=swap;
        nfrontier=nnext;count+=nnext;
    }
    free(frontier);free(next);
    if(out_active)*out_active=count;return active;
}

/* Row-major red/black index lists over movable-and-active vertices.  Red/black
 * Gauss-Seidel updates read only the opposite parity, so sweeping these lists
 * is bitwise identical to the full-grid sweep while visiting only the band. */
static int cloth_build_parity_lists(const uint8_t *movable,
                                    const uint8_t *active,int H,int W,
                                    int32_t **out_red,size_t *out_nred,
                                    int32_t **out_black,size_t *out_nblack) {
    size_t nv=(size_t)H*(size_t)W,count=0;
    for(size_t v=0;v<nv;v++)count+=(size_t)(movable[v]&&active[v]);
    int32_t *red=(int32_t *)malloc((count?count:1u)*sizeof *red);
    int32_t *black=(int32_t *)malloc((count?count:1u)*sizeof *black);
    if(!red||!black){free(red);free(black);return -1;}
    size_t nred=0,nblack=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t v=(size_t)r*(size_t)W+(size_t)c;
        if(!movable[v]||!active[v])continue;
        if(((r+c)&1)==0)red[nred++]=(int32_t)v;
        else black[nblack++]=(int32_t)v;
    }
    *out_red=red;*out_nred=nred;*out_black=black;*out_nblack=nblack;
    return 0;
}

/* Faces with at least one movable-and-active vertex: the only faces a
 * band-confined trial can alter, hence the only faces whose orientation the
 * candidate preflight must re-examine. */
static int32_t *cloth_band_face_list(const int32_t *faces,size_t nf,size_t nv,
                                     const uint8_t *movable,
                                     const uint8_t *active,size_t *out_count) {
    size_t count=0;
    for(size_t f=0;f<nf;f++)for(int k=0;k<3;k++) {
        size_t v=(size_t)faces[f*3+(size_t)k];
        if(v<nv&&movable[v]&&active[v]){count++;break;}
    }
    int32_t *list=(int32_t *)malloc((count?count:1u)*sizeof *list);
    if(!list)return NULL;
    size_t n=0;
    for(size_t f=0;f<nf;f++)for(int k=0;k<3;k++) {
        size_t v=(size_t)faces[f*3+(size_t)k];
        if(v<nv&&movable[v]&&active[v]){list[n++]=(int32_t)f;break;}
    }
    if(out_count)*out_count=n;
    return list;
}

static double cloth_vertex_radius(const float *verts,size_t v,
                                  double axis_y,double axis_x) {
    return hypot((double)verts[v*3+1]-axis_y,
                 (double)verts[v*3+2]-axis_x);
}

static void cloth_move_vertex_radial(float *verts,size_t v,double amount,
                                     double axis_y,double axis_x) {
    double dy=(double)verts[v*3+1]-axis_y;
    double dx=(double)verts[v*3+2]-axis_x;
    double radius=hypot(dy,dx);
    if(!(radius>1e-9)||radius+amount<=1e-6)return;
    double scale=(radius+amount)/radius;
    verts[v*3+1]=(float)(axis_y+scale*dy);
    verts[v*3+2]=(float)(axis_x+scale*dx);
}

static size_t cloth_project_radial_constraints(
        float *displacement,const double *rest_radius,size_t nv,
        const ClothTurnContact *contact,size_t ncontact,int W,
        const uint8_t *movable,const uint8_t *active,double pitch,
        double clearance,double relaxation,double maximum_impulse,
        double maximum_displacement,int reverse,
        size_t *hard_contacts,double *maximum_correction) {
    const double order=pitch>=0.0?1.0:-1.0;
    size_t projected=0,hard=0;double mx=0.0;
    for(size_t at=0;at<ncontact;at++) {
        size_t q=reverse?ncontact-1u-at:at;
        size_t lo[4],hi[4];
        cloth_cell_vertices(contact[q].cell_lo,W,lo);
        cloth_cell_vertices(contact[q].cell_hi,W,hi);
        /* All sixteen quad-vertex inequalities are equivalent to the single
         * support inequality
         *   min(order*r_hi) >= max(order*r_lo) + clearance.
         * Project its currently active extrema.  Re-selecting them every
         * sweep is the usual active-support contact solve and avoids testing
         * duplicate triangle pairs billions of times on a whole scroll. */
        size_t vl=SIZE_MAX,vh=SIZE_MAX;
        double support_lo=-DBL_MAX,support_hi=DBL_MAX;
        for(int kl=0;kl<4;kl++) {
            size_t v=lo[kl];if(v>=nv)continue;
            double s=order*(rest_radius[v]+(double)displacement[v]);
            if(s>support_lo){support_lo=s;vl=v;}
        }
        for(int kh=0;kh<4;kh++) {
            size_t v=hi[kh];if(v>=nv)continue;
            double s=order*(rest_radius[v]+(double)displacement[v]);
            if(s<support_hi){support_hi=s;vh=v;}
        }
        if(vl==SIZE_MAX||vh==SIZE_MAX)continue;
        double C=support_hi-support_lo-clearance;
        if(C>=0.0)continue;
        double wl=(movable[vl]&&active[vl])?1.0:0.0;
        double wh=(movable[vh]&&active[vh])?1.0:0.0;
        if(wl+wh<=0.0){hard++;continue;}
        double impulse=relaxation*(-C)/(wl+wh+1e-6);
        if(impulse>maximum_impulse)impulse=maximum_impulse;
        double dl=-order*wl*impulse,dh=order*wh*impulse;
        if(wl){
            double value=(double)displacement[vl]+dl;
            if(value>maximum_displacement)value=maximum_displacement;
            if(value< -maximum_displacement)value=-maximum_displacement;
            displacement[vl]=(float)value;if(fabs(dl)>mx)mx=fabs(dl);
        }
        if(wh){
            double value=(double)displacement[vh]+dh;
            if(value>maximum_displacement)value=maximum_displacement;
            if(value< -maximum_displacement)value=-maximum_displacement;
            displacement[vh]=(float)value;if(fabs(dh)>mx)mx=fabs(dh);
        }
        projected++;
    }
    if(hard_contacts)*hard_contacts+=hard;
    if(maximum_correction&&mx>*maximum_correction)*maximum_correction=mx;
    return projected;
}

static void cloth_smooth_radial_displacement(
        float *displacement,const int32_t *red,size_t nred,
        const int32_t *black,size_t nblack,int H,int W,
        double screen,double relaxation) {
    for(int parity=0;parity<2;parity++) {
        const int32_t *list=parity?black:red;
        size_t n=parity?nblack:nred;
        ptrdiff_t i;
        /* Same-parity grid points never share an edge, so every update in one
         * phase reads only the frozen opposite parity: identical arithmetic
         * to the full-grid red/black sweep, restricted to the active band. */
#pragma omp parallel for schedule(static)
        for(i=0;i<(ptrdiff_t)n;i++) {
            size_t v=(size_t)list[i];
            int r=(int)(v/(size_t)W),c=(int)(v%(size_t)W);
            double sum=0.0,degree=0.0;
            if(c>0){sum+=(double)displacement[v-1];degree+=1.0;}
            if(c+1<W){sum+=(double)displacement[v+1];degree+=1.0;}
            if(r>0){sum+=(double)displacement[v-(size_t)W];degree+=1.0;}
            if(r+1<H){sum+=(double)displacement[v+(size_t)W];degree+=1.0;}
            double target=sum/(degree+screen);
            displacement[v]=(float)((1.0-relaxation)*(double)displacement[v]+
                                    relaxation*target);
        }
    }
}

static void cloth_apply_radial_displacement(float *out,const float *rest,
                                            size_t nv,const float *displacement,
                                            double axis_y,double axis_x) {
    memcpy(out,rest,nv*3*sizeof *out);
    for(size_t v=0;v<nv;v++)if(displacement[v]!=0.0f)
        cloth_move_vertex_radial(out,v,(double)displacement[v],axis_y,axis_x);
}

static void cloth_constraint_violation(
        const float *verts,size_t nv,const ClothTurnContact *contact,
        size_t ncontact,int W,double pitch,
        double clearance,double axis_y,double axis_x,
        size_t *out_count,double *out_rms,double *out_maximum) {
    const double order=pitch>=0.0?1.0:-1.0;
    size_t count=0;double ss=0.0,mx=0.0;
    for(size_t q=0;q<ncontact;q++) {
        size_t lo[4],hi[4];double support_lo=-DBL_MAX,support_hi=DBL_MAX;
        cloth_cell_vertices(contact[q].cell_lo,W,lo);
        cloth_cell_vertices(contact[q].cell_hi,W,hi);
        for(int k=0;k<4;k++)if(lo[k]<nv) {
            double s=order*cloth_vertex_radius(verts,lo[k],axis_y,axis_x);
            if(s>support_lo)support_lo=s;
        }
        for(int k=0;k<4;k++)if(hi[k]<nv) {
            double s=order*cloth_vertex_radius(verts,hi[k],axis_y,axis_x);
            if(s<support_hi)support_hi=s;
        }
        double deficit=clearance-(support_hi-support_lo);
        if(!(deficit>0.0))continue;
        ss+=deficit*deficit;if(deficit>mx)mx=deficit;count++;
    }
    if(out_count)*out_count=count;
    if(out_rms)*out_rms=count?sqrt(ss/(double)count):0.0;
    if(out_maximum)*out_maximum=mx;
}

static void cloth_displacement_bounds(const float *displacement,size_t nv,
                                      double bound,double *out_min,
                                      double *out_max,size_t *out_at_bound) {
    double lo=DBL_MAX,hi=-DBL_MAX;size_t at=0;
    for(size_t v=0;v<nv;v++) {
        double d=(double)displacement[v];
        if(d<lo)lo=d;if(d>hi)hi=d;
        if(fabs(d)>=bound-1e-4)at++;
    }
    if(out_min)*out_min=lo;if(out_max)*out_max=hi;
    if(out_at_bound)*out_at_bound=at;
}

static void cloth_log_conflict_grid_bbox(const TurnOrderConflictScan *scan,
                                         int W) {
    if(!scan||!scan->count||W<2)return;
    size_t rlo0=SIZE_MAX,rlo1=0,clo0=SIZE_MAX,clo1=0;
    size_t rhi0=SIZE_MAX,rhi1=0,chi0=SIZE_MAX,chi1=0;
    size_t cells_w=(size_t)(W-1);
    for(size_t q=0;q<scan->count;q++) {
        size_t a=scan->pair[q].face_lo/2u,b=scan->pair[q].face_hi/2u;
        size_t ra=a/cells_w,ca=a%cells_w,rb=b/cells_w,cb=b%cells_w;
        if(ra<rlo0)rlo0=ra;if(ra>rlo1)rlo1=ra;
        if(ca<clo0)clo0=ca;if(ca>clo1)clo1=ca;
        if(rb<rhi0)rhi0=rb;if(rb>rhi1)rhi1=rb;
        if(cb<chi0)chi0=cb;if(cb>chi1)chi1=cb;
    }
    fprintf(stderr,
        "[cloth-collision] exact contact UV cells: low rows %zu..%zu cols "
        "%zu..%zu; high rows %zu..%zu cols %zu..%zu\n",
        rlo0,rlo1,clo0,clo1,rhi0,rhi1,chi0,chi1);
}

static void cloth_displacement_stats(const float *p,const float *anchor,
                                     size_t nv,double *rms,double *maximum) {
    double ss=0.0,mx=0.0;for(size_t v=0;v<nv;v++) {
        double d2=0.0;for(int k=0;k<3;k++){
            double d=(double)p[v*3+(size_t)k]-anchor[v*3+(size_t)k];d2+=d*d;
        }
        double d=sqrt(d2);ss+=d2;if(d>mx)mx=d;
    }
    *rms=nv?sqrt(ss/(double)nv):0.0;*maximum=mx;
}

static void cloth_radial_field_stats(const float *displacement,size_t nv,
                                     double *rms,double *maximum) {
    double ss=0.0,mx=0.0;
    for(size_t v=0;v<nv;v++) {
        double d=fabs((double)displacement[v]);
        ss+=d*d;if(d>mx)mx=d;
    }
    if(rms)*rms=nv?sqrt(ss/(double)nv):0.0;
    if(maximum)*maximum=mx;
}

/* Settle toward the measured coil.  Untangling deliberately over-opens the
 * sheet (crop evidence: ~10 vox radial RMS from the pre-escape ribbon while
 * the ledger clearance only requires ~1.25), and the downstream CT snap then
 * starts far off the data -- the 6.49%->17.10% dark regression.  Once the
 * shell is conflict-free, repeatedly pull the radial displacement toward
 * zero (beta per round), re-project the persistent contact ledger to keep
 * the inter-turn clearance, and accept only candidates that a full exact
 * audit certifies as still conflict-free with strictly decreasing radial
 * RMS.  A rejected pull weakens beta ((1+beta)/2) up to three times before
 * the settle stops; the shell's published state is always the last accepted
 * conflict-free geometry. */
static void cloth_settle_toward_rest(
        float *verts,size_t nv,const int32_t *faces,size_t nf,const float *uv,
        const uint8_t *movable,const float *rest_target,
        const double *rest_radius,
        ClothTurnContact **contact_ledger,size_t *ncontact_ledger,
        int H,int W,double pitch,double clearance,
        double axis_y,double axis_x,double minimum_u_separation,
        size_t transaction_input_fold,
        int settle_rounds,double settle_beta,double displacement_bound,
        float *base,float *response,float *trial,float *displacement,
        uint8_t *orientation_input_bad,
        TurnOrderConflictScan *cached_scan,
        IntersectionCleanupStats *cached_stats,
        ClothCollisionStats *s) {
    if(settle_rounds<1||!contact_ledger||!ncontact_ledger||
       *ncontact_ledger==0||!(settle_beta>0.0&&settle_beta<1.0))
        return;
    double beta=settle_beta;
    int betas_weakened=0,accepted_rounds=0;
    size_t previous_reject_conflicts=SIZE_MAX;
    int stale_merges=0;
    double previous_rms=-1.0;
    for(int round=0;round<settle_rounds;round++) {
        const ClothTurnContact *contact=*contact_ledger;
        size_t ncontact=*ncontact_ledger;
        memcpy(base,verts,nv*3*sizeof *base);
        size_t nactive=0;
        uint8_t *active=cloth_active_contact_band(
            nv,contact,ncontact,movable,H,W,24,&nactive);
        if(!active)return;
        for(size_t v=0;v<nv;v++) {
            double current_r=cloth_vertex_radius(base,v,axis_y,axis_x);
            displacement[v]=(float)(current_r-rest_radius[v]);
            if(movable[v]&&fabs((double)displacement[v])>1e-5)active[v]=1;
        }
        double rms_before=0.0,max_before=0.0;
        cloth_radial_field_stats(displacement,nv,&rms_before,&max_before);
        if(round==0) {
            s->settle_rms_before=rms_before;s->settle_max_before=max_before;
            s->settle_rms_after=rms_before;s->settle_max_after=max_before;
            previous_rms=rms_before;
        }
        if(!(rms_before>1e-6)){free(active);break;}
        int32_t *red_list=NULL,*black_list=NULL,*band_face=NULL;
        size_t nred=0,nblack=0,nband_face=0;
        if(cloth_build_parity_lists(movable,active,H,W,
                                    &red_list,&nred,&black_list,&nblack)!=0) {
            free(active);return;
        }
        band_face=cloth_band_face_list(faces,nf,nv,movable,active,&nband_face);
        if(!band_face){free(red_list);free(black_list);free(active);return;}
        size_t input_bad_faces=post_snap_orientation_baseline(
            base,nv,faces,nf,orientation_input_bad);
        (void)input_bad_faces;
        for(int half=0;half<2;half++) {
            const int32_t *wl=half?black_list:red_list;
            ptrdiff_t wn=(ptrdiff_t)(half?nblack:nred),bi;
#pragma omp parallel for schedule(static)
            for(bi=0;bi<wn;bi++) {
                size_t v=(size_t)wl[bi];
                displacement[v]=(float)((double)displacement[v]*beta);
            }
        }
        size_t hard=0,projected=0;double max_impulse=0.0;
        for(int sweep=0;sweep<96;sweep++) {
            cloth_smooth_radial_displacement(
                displacement,red_list,nred,black_list,nblack,H,W,0.05,0.70);
            for(int cs=0;cs<2;cs++)projected+=
                cloth_project_radial_constraints(
                    displacement,rest_radius,nv,contact,ncontact,W,movable,
                    active,pitch,clearance,0.65,0.35,displacement_bound,
                    (sweep+cs)&1,&hard,&max_impulse);
        }
        for(int cs=0;cs<32;cs++)projected+=
            cloth_project_radial_constraints(
                displacement,rest_radius,nv,contact,ncontact,W,movable,active,
                pitch,clearance,0.90,0.5,displacement_bound,cs&1,
                &hard,&max_impulse);
        cloth_apply_radial_displacement(
            response,rest_target,nv,displacement,axis_y,axis_x);
        memcpy(trial,base,nv*3*sizeof *trial);
        for(int half=0;half<2;half++) {
            const int32_t *wl=half?black_list:red_list;
            ptrdiff_t wn=(ptrdiff_t)(half?nblack:nred),bi;
#pragma omp parallel for schedule(static)
            for(bi=0;bi<wn;bi++) {
                size_t v=(size_t)wl[bi];
                for(int k=0;k<3;k++)
                    trial[v*3+(size_t)k]=response[v*3+(size_t)k];
            }
        }
        double rms_after=0.0,max_after=0.0;
        cloth_radial_field_stats(displacement,nv,&rms_after,&max_after);
        size_t orientation_bad=post_snap_orientation_violations_list(
            trial,base,nv,faces,band_face,nband_face,
            orientation_input_bad,NULL,NULL);
        int accept=0;
        TurnOrderConflictScan after_scan={0};
        IntersectionCleanupStats after={0};
        int audited=0;
        if(orientation_bad==0&&rms_after+1e-9<0.995*previous_rms) {
            if(turn_order_audit(trial,nv,faces,nf,uv,minimum_u_separation,
                                NULL,NULL,&after_scan,&after)!=0) {
                free(red_list);free(black_list);free(band_face);free(active);
                return;
            }
            audited=1;
            accept=after_scan.count==0&&after.stab_pairs==0&&
                   after.fold_pairs<=transaction_input_fold;
        }
        fprintf(stderr,
            "[elastic-shell] settle round %d: beta %.4f, displacement rms/max "
            "%.4f/%.4f -> %.4f/%.4f, exact long/stab/fold=%zu/%zu/%zu (%s)\n",
            round+1,beta,rms_before,max_before,rms_after,max_after,
            audited?after_scan.count:(size_t)0,
            audited?after.stab_pairs:(size_t)0,
            audited?after.fold_pairs:(size_t)0,
            accept?"accept":(orientation_bad?"reject-orientation":
                             (audited?"reject-audit":"reject-no-progress")));
        s->settle_rounds_run++;
        if(accept) {
            memcpy(verts,trial,nv*3*sizeof *verts);
            free(cached_scan->pair);
            *cached_scan=after_scan;*cached_stats=after;
            s->settle_rms_after=rms_after;s->settle_max_after=max_after;
            accepted_rounds++;
            previous_rms=rms_after;
            previous_reject_conflicts=SIZE_MAX;stale_merges=0;
            free(red_list);free(black_list);free(band_face);free(active);
            continue;
        }
        if(audited&&(after_scan.count>0||after.stab_pairs>0)) {
            /* The return path crossed pairs the escape ledger never saw --
             * the pull re-collides where the coil originally interpenetrated.
             * Adopt the failed trial's exact conflicts as contacts and retry:
             * the projection then holds clearance exactly there, so the next
             * pull relaxes only through genuinely free space.  Only when
             * merging stops shrinking the rejection do we weaken the pull. */
            int made_progress=after.conflicts<previous_reject_conflicts;
            previous_reject_conflicts=after.conflicts;
            size_t nnew=0;
            ClothTurnContact *new_contact=cloth_build_contacts(
                trial,&after_scan,H,W,4.0,3,&nnew);
            free(after_scan.pair);
            if(new_contact) {
                ClothTurnContact *merged=NULL;size_t nmerged=0;
                if(cloth_merge_contact_ledger(
                        *contact_ledger,*ncontact_ledger,new_contact,nnew,
                        &merged,&nmerged)==0) {
                    free(*contact_ledger);
                    *contact_ledger=merged;*ncontact_ledger=nmerged;
                    fprintf(stderr,
                        "[elastic-shell] settle adopted %zu return-path "
                        "contacts (ledger %zu)\n",nnew,nmerged);
                }
                free(new_contact);
            }
            free(red_list);free(black_list);free(band_face);free(active);
            if(made_progress){stale_merges=0;continue;}
            if(++stale_merges<2)continue;
            stale_merges=0;
            betas_weakened++;
            if(betas_weakened>3)break;
            beta=(1.0+beta)*0.5;
            continue;
        }
        if(audited)free(after_scan.pair);
        free(red_list);free(black_list);free(band_face);free(active);
        /* A rejected pull is either blocked by the ledger clearance floor or
         * pulled too hard through it; weaken toward the identity and retry,
         * then stop -- the accepted state already satisfies every gate. */
        betas_weakened++;
        if(betas_weakened>3)break;
        beta=(1.0+beta)*0.5;
    }
    s->settle_accepted=accepted_rounds;
    fprintf(stderr,
        "[elastic-shell] settle result: rounds=%d accepted=%d, displacement "
        "rms %.4f->%.4f, max %.4f->%.4f\n",
        s->settle_rounds_run,accepted_rounds,
        s->settle_rms_before,s->settle_rms_after,
        s->settle_max_before,s->settle_max_after);
}

static int resolve_cylindrical_self_collisions(
        float *verts,size_t nv,const int32_t *faces,size_t nf,const float *uv,
        const float *phase,const uint8_t *movable,const float *rest_target,
        int H,int W,double grid_du,
        double axis_y,double axis_x,double minimum_u_separation,
        int require_collision_free,int collision_patience,int collision_rounds,
        double collision_collar,double collision_displacement_bound,
        int collision_settle_rounds,double collision_settle_beta,
        ClothCollisionStats *out) {
    ClothCollisionStats s;memset(&s,0,sizeof s);s.attempted=1;
    if(out)*out=s;
    if(!verts||!faces||!uv||!movable||!rest_target||H<2||W<3||
       nv!=(size_t)H*(size_t)W)return -1;
    float *base=(float *)malloc(nv*3*sizeof *base);
    float *response=(float *)malloc(nv*3*sizeof *response);
    float *trial=(float *)malloc(nv*3*sizeof *trial);
    float *transaction_start=(float *)malloc(nv*3*sizeof *transaction_start);
    float *displacement=(float *)calloc(nv,sizeof *displacement);
    double *rest_radius=(double *)malloc(nv*sizeof *rest_radius);
    float *keep=(float *)malloc(nv*sizeof *keep);
    float *keep_tmp=(float *)malloc(nv*sizeof *keep_tmp);
    size_t *degree_before=(size_t *)calloc(nf,sizeof *degree_before);
    size_t *degree_after=(size_t *)calloc(nf,sizeof *degree_after);
    size_t *local_before=(size_t *)calloc(nf,sizeof *local_before);
    size_t *local_after=(size_t *)calloc(nf,sizeof *local_after);
    uint8_t *orientation_input_bad=(uint8_t *)malloc(nf);
    if(!base||!response||!trial||!transaction_start||!displacement||
       !rest_radius||
       !keep||!keep_tmp||
       !degree_before||!degree_after||
       !local_before||!local_after||!orientation_input_bad) {
        free(base);free(response);free(trial);free(transaction_start);free(displacement);
        free(rest_radius);
        free(keep);free(keep_tmp);
        free(degree_before);free(degree_after);
        free(local_before);free(local_after);
        free(orientation_input_bad);return -1;
    }
    audit_totals_reset();
    /* Strict mode retains the entry state here for atomic rollback.  Default
     * inspectable-output mode updates the same buffer only when a fully audited
     * candidate improves the exact collision objective. */
    memcpy(transaction_start,verts,nv*3*sizeof *transaction_start);
    for(size_t v=0;v<nv;v++)
        rest_radius[v]=cloth_vertex_radius(rest_target,v,axis_y,axis_x);
    double pitch=0.0,clearance=1.25;int dcol=0;
    size_t transaction_input_fold=0;
    size_t best_long=SIZE_MAX;
    size_t best_violation_count=SIZE_MAX;
    IntersectionCleanupStats best_stats={0};
    int best_initialized=0,rounds_since_best=0;
    const double collar_base=collision_collar>0.0?collision_collar:4.0;
    const double collar_max=collar_base+8.0;
    double collar=collar_base;
    int contact_inflate=3,band_rings=24,stall_rounds=0;
    size_t prev_accepted_long=0;int have_prev_accepted_long=0;
    ClothTurnContact *contact_ledger=NULL;size_t ncontact_ledger=0;
    TurnOrderConflictScan cached_scan={0};
    IntersectionCleanupStats cached_stats={0};
    int cached_audit=0;
    for(int round=0;round<collision_rounds;round++) {
        TurnOrderConflictScan before_scan={0};IntersectionCleanupStats before={0};
        if(cached_audit) {
            before_scan=cached_scan;
            before=cached_stats;
            memset(&cached_scan,0,sizeof cached_scan);
            memset(&cached_stats,0,sizeof cached_stats);
            cached_audit=0;
            fprintf(stderr,
                "[cloth-collision] reusing accepted exact audit for next active set\n");
        } else if(turn_order_audit(
                    verts,nv,faces,nf,uv,minimum_u_separation,
                    degree_before,local_before,&before_scan,&before)!=0) {
            goto fail;
        }
        if(round==0){
            s.input_conflicts=before.conflicts;
            s.input_long_conflicts=before_scan.count;
            transaction_input_fold=before.fold_pairs;
            best_long=before_scan.count;best_stats=before;best_initialized=1;
            if(before_scan.count&&turn_order_infer_pitch(
                    verts,phase,movable,H,W,grid_du,axis_y,axis_x,
                    &before_scan,&dcol,&pitch)!=0) {
                fprintf(stderr,"[cloth-collision] refusing ambiguous phase ordering\n");
                free(before_scan.pair);break;
            }
            s.pitch=pitch;
            /* Numerical mid-surface clearance, deliberately much smaller than
             * the measured inter-turn pitch.  The rest energy, not clearance,
             * determines where the metal settles after it opens. */
            clearance=1.25;
            s.clearance=clearance;
        }
        fprintf(stderr,
            "[cloth-collision] active-set round %d: exact=%zu, long=%zu "
            "[overlap=%zu stab=%zu fold=%zu], clearance=%.3f\n",
            round+1,before.conflicts,before_scan.count,before.overlap_pairs,
            before.stab_pairs,before.fold_pairs,clearance);
        cloth_log_conflict_grid_bbox(&before_scan,W);
        if(before_scan.count==0){
            cached_scan=before_scan;cached_stats=before;cached_audit=1;
            break;
        }
        size_t nnew_contact=0;
        ClothTurnContact *new_contact=cloth_build_contacts(
            verts,&before_scan,H,W,collar,contact_inflate,&nnew_contact);
        if(!new_contact){free(before_scan.pair);goto fail;}
        ClothTurnContact *merged_contact=NULL;size_t nmerged_contact=0;
        if(cloth_merge_contact_ledger(
                contact_ledger,ncontact_ledger,new_contact,nnew_contact,
                &merged_contact,&nmerged_contact)!=0) {
            free(new_contact);free(before_scan.pair);goto fail;
        }
        free(new_contact);free(contact_ledger);contact_ledger=merged_contact;
        ncontact_ledger=nmerged_contact;
        ClothTurnContact *contact=contact_ledger;
        size_t ncontact=ncontact_ledger;
        size_t nactive=0;
        uint8_t *active=cloth_active_contact_band(
            nv,contact,ncontact,movable,H,W,band_rings,&nactive);
        if(!active){free(before_scan.pair);goto fail;}
        memcpy(base,verts,nv*3*sizeof *base);
        {
            size_t input_bad_faces=post_snap_orientation_baseline(
                base,nv,faces,nf,orientation_input_bad);
            if(input_bad_faces)fprintf(stderr,
                "[cloth-collision] input orientation-bad faces=%zu excluded "
                "from candidate preflight\n",input_bad_faces);
        }
        nactive=0;
        for(size_t v=0;v<nv;v++) {
            double current_r=cloth_vertex_radius(base,v,axis_y,axis_x);
            double target_r=rest_radius[v];
            displacement[v]=(float)(current_r-target_r);
            /* The old coarse escape is only an initializer.  Every vertex it
             * moved belongs to the elastic return solve even if it is no longer
             * on today's exact contact curve. */
            if(movable[v]&&fabs(current_r-target_r)>1e-5)active[v]=1;
            nactive+=active[v]!=0;
        }
        int32_t *red_list=NULL,*black_list=NULL,*band_face=NULL;
        size_t nred=0,nblack=0,nband_face=0;
        if(cloth_build_parity_lists(movable,active,H,W,
                                    &red_list,&nred,&black_list,&nblack)!=0) {
            free(active);free(before_scan.pair);goto fail;
        }
        band_face=cloth_band_face_list(faces,nf,nv,movable,active,&nband_face);
        if(!band_face) {
            free(red_list);free(black_list);
            free(active);free(before_scan.pair);goto fail;
        }
        size_t hard=0,projected=0;double max_impulse=0.0;
        const double displacement_bound=
            collision_displacement_bound>0.0?collision_displacement_bound:128.0;
        for(int sweep=0;sweep<512;sweep++) {
            cloth_smooth_radial_displacement(
                displacement,red_list,nred,black_list,nblack,H,W,0.006,0.70);
            for(int cs=0;cs<2;cs++)projected+=
                cloth_project_radial_constraints(
                    displacement,rest_radius,nv,contact,ncontact,W,movable,active,
                    pitch,clearance,0.65,0.35,displacement_bound,
                    (sweep+cs)&1,&hard,&max_impulse);
        }
        /* Finish on the inequalities, not on a smoothing sweep. */
        for(int cs=0;cs<128;cs++)projected+=
            cloth_project_radial_constraints(
                displacement,rest_radius,nv,contact,ncontact,W,movable,active,
                pitch,clearance,0.90,0.5,displacement_bound,cs&1,
                &hard,&max_impulse);
        cloth_apply_radial_displacement(
            response,rest_target,nv,displacement,axis_y,axis_x);
        s.contacts_built+=ncontact;s.active_vertices=nactive;
        s.hard_contacts+=hard;s.elastic_sweeps+=512;s.contact_sweeps+=1152;
        double full_rms=0.0,full_max=0.0;
        cloth_displacement_stats(response,base,nv,&full_rms,&full_max);
        double disp_min=0.0,disp_max=0.0;size_t disp_at_bound=0;
        cloth_displacement_bounds(displacement,nv,displacement_bound,
                                  &disp_min,&disp_max,&disp_at_bound);
        size_t violation_count_before=0;double violation_rms_before=0.0;
        double violation_max_before=0.0;
        cloth_constraint_violation(
            base,nv,contact,ncontact,W,pitch,clearance,axis_y,axis_x,
            &violation_count_before,&violation_rms_before,&violation_max_before);
        fprintf(stderr,
            "[cloth-collision] contacts=%zu, active=%zu/%zu, projections=%zu, "
            "hard=%zu, response rms/max=%.6f/%.6f, max impulse %.6f; "
            "constraint violations=%zu rms/max %.6f/%.6f\n",
            ncontact,nactive,nv,projected,hard,full_rms,full_max,
            max_impulse,violation_count_before,violation_rms_before,
            violation_max_before);
        fprintf(stderr,
            "[elastic-shell] radial displacement from original min/max "
            "%.6f/%.6f vox; at +/-%.1f bound=%zu\n",
            disp_min,disp_max,displacement_bound,disp_at_bound);
        const double alpha[7]={1.0,0.75,0.5,0.375,0.25,0.125,0.0625};
        /* Whole-sheet evidence shows that a large, useful response can be down
         * to only a handful of unsafe front vertices at rollback six.  Cutting
         * it off there discards the response and replaces one local repair with
         * several more global contact rounds.  Continue the monotonically
         * feathered local transaction; every retry is still independently
         * orientation- and exact-intersection-audited before acceptance. */
        const int maximum_local_rollbacks=12;
        int accepted=0;
        /* Band-confined trials: static vertices remain bitwise equal to the
         * round input, so one base copy serves every attempt and rejected
         * trials need no restore -- the next trial rewrites the same band. */
        memcpy(trial,base,nv*3*sizeof *trial);
        for(int attempt=0;attempt<7;attempt++) {
            for(size_t v=0;v<nv;v++)keep[v]=1.0f;
            for(int rollback=0;rollback<=maximum_local_rollbacks;rollback++) {
                for(int half=0;half<2;half++) {
                    const int32_t *wl=half?black_list:red_list;
                    ptrdiff_t wn=(ptrdiff_t)(half?nblack:nred),bi;
#pragma omp parallel for schedule(static)
                    for(bi=0;bi<wn;bi++) {
                        size_t v=(size_t)wl[bi];
                        for(int k=0;k<3;k++)
                            trial[v*3+(size_t)k]=
                                (float)((double)base[v*3+(size_t)k]+
                                    alpha[attempt]*(double)keep[v]*
                                    ((double)response[v*3+(size_t)k]-
                                     base[v*3+(size_t)k]));
                    }
                }
                size_t orientation_bad=post_snap_orientation_violations_list(
                    trial,base,nv,faces,band_face,nband_face,
                    orientation_input_bad,NULL,NULL);
                if(orientation_bad) {
                    fprintf(stderr,
                        "[cloth-collision] reject alpha %.4f rollback %d "
                        "before exact audit: orientation violations=%zu\n",
                        alpha[attempt],rollback,orientation_bad);
                    if(rollback==maximum_local_rollbacks)break;
                    size_t orientation_marked=0;
                    (void)post_snap_orientation_violations_list(
                        trial,base,nv,faces,band_face,nband_face,
                        orientation_input_bad,keep,&orientation_marked);
                    if(orientation_marked==0)break;
                    feather_motion_scale(keep,keep_tmp,nv,H,W,8);
                    fprintf(stderr,
                        "[cloth-collision] preflight rolled back %zu orientation "
                        "vertices and feathered 8 grid rings\n",
                        orientation_marked);
                    continue;
                }
                TurnOrderConflictScan after_scan={0};
                IntersectionCleanupStats after={0};
                if(turn_order_audit(trial,nv,faces,nf,uv,minimum_u_separation,
                                    degree_after,local_after,&after_scan,&after)!=0) {
                    free(red_list);free(black_list);free(band_face);
                    free(active);free(before_scan.pair);goto fail;
                }
                size_t before_local=before.conflicts-before_scan.count;
                size_t after_local=after.conflicts-after_scan.count;
                size_t violation_count_after=0;double violation_rms_after=0.0;
                double violation_max_after=0.0;
                cloth_constraint_violation(
                    trial,nv,contact,ncontact,W,pitch,clearance,axis_y,axis_x,
                    &violation_count_after,&violation_rms_after,&violation_max_after);
                /* Existing penetration has no collision-free continuous endpoint
                 * step: its intersection curve must sweep across neighbouring
                 * tessellation faces before leaving the sheet.  Permit that pair
                 * identity/count churn inside this atomic transaction; phase
                 * violation and orientation, rather than pair monotonicity, are
                 * the progress/safety measures. */
                size_t churn_cap=before.conflicts<=SIZE_MAX/8u
                                ?before.conflicts*8u:SIZE_MAX-before.conflicts;
                if(churn_cap<4096u)churn_cap=4096u;
                if(churn_cap>SIZE_MAX-before.conflicts)
                    churn_cap=SIZE_MAX-before.conflicts;
                int phase_progress=violation_count_after<violation_count_before||
                    violation_rms_after+1e-6<0.995*violation_rms_before||
                    after_scan.count<before_scan.count;
                int bounded_churn=after.conflicts<=before.conflicts+churn_cap;
                int safe=phase_progress&&bounded_churn&&
                         after.fold_pairs<=before.fold_pairs&&orientation_bad==0;
                fprintf(stderr,
                    "[cloth-collision] %s alpha %.4f rollback %d: exact "
                    "%zu->%zu, long %zu->%zu, local %zu->%zu, folds "
                    "%zu->%zu, orientation=%zu\n",
                    safe?"accept":"reject",alpha[attempt],rollback,
                    before.conflicts,after.conflicts,before_scan.count,
                    after_scan.count,before_local,after_local,before.fold_pairs,
                    after.fold_pairs,orientation_bad);
                fprintf(stderr,
                    "[cloth-collision]   phase violations %zu->%zu, rms "
                    "%.6f->%.6f, max %.6f->%.6f (temporary churn cap %zu)\n",
                    violation_count_before,violation_count_after,
                    violation_rms_before,violation_rms_after,
                    violation_max_before,violation_max_after,
                    before.conflicts+churn_cap);
                if(safe){
                    memcpy(verts,trial,nv*3*sizeof *verts);accepted=1;
                    s.accepted_rounds++;
                    double rms=0.0,mx=0.0;
                    cloth_displacement_stats(verts,base,nv,&rms,&mx);
                    s.movement_rms=rms;s.movement_max=mx;
                    if(!require_collision_free) {
                        int lex_better=after_scan.count<best_long||
                            (after_scan.count==best_long&&
                             (after.stab_pairs<best_stats.stab_pairs||
                              (after.stab_pairs==best_stats.stab_pairs&&
                               after.conflicts<best_stats.conflicts)));
                        /* The escape is deliberately non-monotone in pair
                         * counts: while the intersection front travels, phase
                         * violations are the trustworthy progress signal
                         * (full-sheet v139 collapsed 661k->89k violations in
                         * one round yet was stopped by lexicographic-only
                         * patience).  Either kind of new best resets patience;
                         * only the lexicographic best is ever published. */
                        int phase_better=
                            violation_count_after<best_violation_count;
                        if(lex_better) {
                            /* In inspectable-output mode transaction_start is
                             * the best fully audited state, not merely the
                             * entry state.  A partial solve can therefore
                             * never publish a later active-set oscillation
                             * that is worse than its input. */
                            memcpy(transaction_start,trial,nv*3*sizeof *trial);
                            best_long=after_scan.count;best_stats=after;
                        }
                        if(lex_better||phase_better)rounds_since_best=0;
                        else rounds_since_best++;
                        if(phase_better)
                            best_violation_count=violation_count_after;
                        fprintf(stderr,
                            "[cloth-collision] patience %d/%d; best "
                            "long/stab/total=%zu/%zu/%zu, best "
                            "phase-violations=%zu\n",
                            rounds_since_best,collision_patience,
                            best_long,best_stats.stab_pairs,
                            best_stats.conflicts,best_violation_count);
                    }
                    if(have_prev_accepted_long&&
                       cloth_stall_escalation_step(
                           prev_accepted_long,after_scan.count,
                           collar_base,collar_max,&stall_rounds,
                           &collar,&contact_inflate,&band_rings))
                        fprintf(stderr,
                            "[cloth-collision] front stall response: "
                            "collar=%.1f inflate=%d rings=%d\n",
                            collar,contact_inflate,band_rings);
                    prev_accepted_long=after_scan.count;
                    have_prev_accepted_long=1;
                    /* This exact audit describes the accepted geometry.  Carry
                     * both its scan and per-face degrees into the next active
                     * set instead of immediately rebuilding the same BVH and
                     * rediscovering the same contacts. */
                    cached_scan=after_scan;cached_stats=after;cached_audit=1;
                    {
                        size_t *tmp=degree_before;
                        degree_before=degree_after;degree_after=tmp;
                        tmp=local_before;
                        local_before=local_after;local_after=tmp;
                    }
                    break;
                }
                if(rollback==maximum_local_rollbacks||
                   !phase_progress||!bounded_churn) {
                    free(after_scan.pair);break;
                }
                /* A full-sheet contact move can safely untangle one winding
                 * while folding a few cells at the moving intersection front.
                 * Freeze only vertices on newly-created local conflicts or
                 * orientation failures, feather that freeze into the response,
                 * and re-audit.  Long-range pair churn is deliberately not
                 * marked: it is the intersection curve travelling to the edge. */
                size_t marked=0,orientation_marked=0;
                for(size_t f=0;f<nf;f++)if(local_after[f]>local_before[f])
                    for(int k=0;k<3;k++) {
                        size_t v=(size_t)faces[f*3+(size_t)k];
                        if(keep[v]>0.0f){keep[v]=0.0f;marked++;}
                    }
                (void)post_snap_orientation_violations_list(
                    trial,base,nv,faces,band_face,nband_face,
                    orientation_input_bad,keep,&orientation_marked);
                marked+=orientation_marked;
                free(after_scan.pair);
                if(marked==0)break;
                feather_motion_scale(keep,keep_tmp,nv,H,W,8);
                fprintf(stderr,
                    "[cloth-collision] locally rolled back %zu fold/orientation "
                    "vertices and feathered 8 grid rings\n",marked);
            }
            if(accepted)break;
        }
        free(red_list);free(black_list);free(band_face);
        free(active);
        if(accepted)free(before_scan.pair);
        if(accepted&&!require_collision_free&&collision_patience>0&&
           rounds_since_best>=collision_patience) {
            fprintf(stderr,
                "[cloth-collision] stopping after %d accepted rounds without "
                "a better exact-audited topology (patience=%d); retaining "
                "the best state\n",
                rounds_since_best,collision_patience);
            break;
        }
        if(!accepted){
            fprintf(stderr,"[cloth-collision] no exact-audit-safe response remained\n");
            /* Geometry did not change, so the round-input audit is also the
             * transaction-final audit. */
            cached_scan=before_scan;cached_stats=before;cached_audit=1;
            break;
        }
    }
    if(cached_audit&&cached_scan.count==0&&cached_stats.stab_pairs==0&&
       ncontact_ledger>0&&pitch!=0.0)
        cloth_settle_toward_rest(
            verts,nv,faces,nf,uv,movable,rest_target,rest_radius,
            &contact_ledger,&ncontact_ledger,H,W,pitch,clearance,
            axis_y,axis_x,minimum_u_separation,transaction_input_fold,
            collision_settle_rounds,collision_settle_beta,
            collision_displacement_bound>0.0?collision_displacement_bound:128.0,
            base,response,trial,displacement,orientation_input_bad,
            &cached_scan,&cached_stats,&s);
    {
        TurnOrderConflictScan final_scan={0};IntersectionCleanupStats final={0};
        if(cached_audit) {
            final_scan=cached_scan;final=cached_stats;
            memset(&cached_scan,0,sizeof cached_scan);cached_audit=0;
        } else if(turn_order_audit(
                    verts,nv,faces,nf,uv,minimum_u_separation,
                    NULL,NULL,&final_scan,&final)!=0) {
            goto fail;
        }
        int complete=final_scan.count==0&&final.stab_pairs==0&&
                     final.fold_pairs<=transaction_input_fold;
        if(!complete&&!require_collision_free&&best_initialized) {
            free(final_scan.pair);memset(&final_scan,0,sizeof final_scan);
            memcpy(verts,transaction_start,nv*3*sizeof *verts);
            final=best_stats;final_scan.count=best_long;
            s.retained_partial=1;
            fprintf(stderr,
                "[cloth-collision] PARTIAL: retained best exact-audited state; "
                "remaining long=%zu, stab=%zu, folds=%zu (input folds=%zu)\n",
                final_scan.count,final.stab_pairs,final.fold_pairs,
                transaction_input_fold);
        }
        s.complete=complete;
        s.output_conflicts=final.conflicts;s.output_long_conflicts=final_scan.count;
        fprintf(stderr,
            "[cloth-collision] result: accepted rounds=%d, exact %zu->%zu, "
            "long %zu->%zu [overlap=%zu stab=%zu fold=%zu]\n",
            s.accepted_rounds,s.input_conflicts,s.output_conflicts,
            s.input_long_conflicts,s.output_long_conflicts,final.overlap_pairs,
            final.stab_pairs,final.fold_pairs);
        audit_totals_report("elastic-shell");
        {
            double rest_rms=0.0,rest_max=0.0;
            cloth_displacement_stats(verts,rest_target,nv,&rest_rms,&rest_max);
            fprintf(stderr,
                "[elastic-shell] constrained return distance to original "
                "rms/max %.6f/%.6f vox\n",rest_rms,rest_max);
        }
        if(!complete&&require_collision_free) {
            fprintf(stderr,
                "[cloth-collision] ROLLBACK: the elastic-shell transaction did "
                "not reach a non-interpenetrating state; restoring its input\n");
            memcpy(verts,transaction_start,nv*3*sizeof *verts);
            free(final_scan.pair);
            free(contact_ledger);
            free(base);free(response);free(trial);free(transaction_start);
            free(displacement);free(rest_radius);free(keep);free(keep_tmp);
            free(degree_before);free(degree_after);
            free(local_before);free(local_after);
            free(orientation_input_bad);if(out)*out=s;return 1;
        }
        free(final_scan.pair);
    }
    free(contact_ledger);
    free(base);free(response);free(trial);free(transaction_start);free(displacement);
    free(rest_radius);
    free(keep);free(keep_tmp);
    free(degree_before);free(degree_after);
    free(local_before);free(local_after);
    free(orientation_input_bad);if(out)*out=s;return 0;
fail:
    free(cached_scan.pair);
    free(contact_ledger);
    free(base);free(response);free(trial);free(transaction_start);free(displacement);
    free(rest_radius);
    free(keep);free(keep_tmp);
    free(degree_before);free(degree_after);
    free(local_before);free(local_after);
    free(orientation_input_bad);return -1;
}

static int cloth_support_selftest(void) {
    const int W=4;
    const double rest_radius[8]={10.0,12.0,9.0,15.0,11.0,10.0,14.0,13.0};
    const ClothTurnContact contact={0,2};
    const uint8_t movable[8]={1,1,1,1,1,1,1,1};
    const uint8_t active[8]={1,1,1,1,1,1,1,1};
    float displacement[8]={0,0,0,0,0,0,0,0};
    size_t projected=0,hard=0;double maximum=0.0;
    for(int sweep=0;sweep<64;sweep++)
        projected+=cloth_project_radial_constraints(
            displacement,rest_radius,8,&contact,1,W,movable,active,
            8.0,1.25,1.0,100.0,128.0,sweep&1,&hard,&maximum);
    double max_lo=-DBL_MAX,min_hi=DBL_MAX;
    const size_t lo_vertex[4]={0,1,4,5},hi_vertex[4]={2,3,6,7};
    for(int k=0;k<4;k++) {
        double lo=rest_radius[lo_vertex[k]]+(double)displacement[lo_vertex[k]];
        double hi=rest_radius[hi_vertex[k]]+(double)displacement[hi_vertex[k]];
        if(lo>max_lo)max_lo=lo;if(hi<min_hi)min_hi=hi;
    }
    int fail=projected==0||hard!=0||maximum<=0.0||
             min_hi-max_lo<1.25-1e-5;
    fprintf(stderr,"  %s: quad-compressed cloth projection satisfies all sixteen inequalities\n",
            fail?"FAIL":"ok");
    return fail;
}

static int cloth_settle_selftest(void) {
    const int H=2,W=4;const size_t nv=8;
    const double rest_radius[8]={10.0,12.0,9.0,15.0,11.0,10.0,14.0,13.0};
    const ClothTurnContact contact={0,2};
    const uint8_t movable[8]={1,1,1,1,1,1,1,1};
    const uint8_t active[8]={1,1,1,1,1,1,1,1};
    const double pitch=8.0,clearance=1.25,beta=0.70;
    /* Over-opened exactly the way the escape leaves the sheet: the hi cell
     * is 10 vox beyond its rest coil while the constraint only needs the
     * clearance. */
    float displacement[8]={-2,-2,10,10,-2,-2,10,10};
    double rms_before=0.0,mx_before=0.0;
    cloth_radial_field_stats(displacement,nv,&rms_before,&mx_before);
    int32_t *red=NULL,*black=NULL;size_t nred=0,nblack=0;
    if(cloth_build_parity_lists(movable,active,H,W,
                                &red,&nred,&black,&nblack)!=0)return 1;
    for(size_t v=0;v<nv;v++)
        displacement[v]=(float)((double)displacement[v]*beta);
    size_t projected=0,hard=0;double mximp=0.0;
    for(int sweep=0;sweep<96;sweep++) {
        cloth_smooth_radial_displacement(displacement,red,nred,black,nblack,
                                         H,W,0.05,0.70);
        for(int cs=0;cs<2;cs++)projected+=
            cloth_project_radial_constraints(
                displacement,rest_radius,nv,&contact,1,W,movable,active,
                pitch,clearance,0.65,0.35,128.0,(sweep+cs)&1,&hard,&mximp);
    }
    for(int cs=0;cs<32;cs++)projected+=
        cloth_project_radial_constraints(
            displacement,rest_radius,nv,&contact,1,W,movable,active,
            pitch,clearance,0.90,0.5,128.0,cs&1,&hard,&mximp);
    free(red);free(black);
    double rms_after=0.0,mx_after=0.0;
    cloth_radial_field_stats(displacement,nv,&rms_after,&mx_after);
    double max_lo=-DBL_MAX,min_hi=DBL_MAX;
    const size_t lo_vertex[4]={0,1,4,5},hi_vertex[4]={2,3,6,7};
    for(int k=0;k<4;k++) {
        double lo=rest_radius[lo_vertex[k]]+(double)displacement[lo_vertex[k]];
        double hi=rest_radius[hi_vertex[k]]+(double)displacement[hi_vertex[k]];
        if(lo>max_lo)max_lo=lo;if(hi<min_hi)min_hi=hi;
    }
    int fail=!(rms_after<rms_before)||min_hi-max_lo<clearance-1e-5;
    fprintf(stderr,
        "  %s: settle pull reduces radial rms (%.3f->%.3f) and keeps the "
        "clearance (%.3f)\n",fail?"FAIL":"ok",rms_before,rms_after,
        min_hi-max_lo);
    return fail;
}

static int cloth_stall_escalation_selftest(void) {
    double collar=4.0;int inflate=3,rings=24,stall=0,fail=0;
    /* Three consecutive <2% rounds escalate once. */
    int changed=0;
    changed|=cloth_stall_escalation_step(1000,995,4.0,12.0,&stall,
                                         &collar,&inflate,&rings);
    changed|=cloth_stall_escalation_step(995,990,4.0,12.0,&stall,
                                         &collar,&inflate,&rings);
    if(changed||collar!=4.0){fail++;}
    changed=cloth_stall_escalation_step(990,985,4.0,12.0,&stall,
                                        &collar,&inflate,&rings);
    if(!changed||collar!=8.0||inflate!=5||rings!=48)fail++;
    /* A middling round resets the stall counter without de-escalating. */
    changed=cloth_stall_escalation_step(985,900,4.0,12.0,&stall,
                                        &collar,&inflate,&rings);
    if(changed||stall!=0||collar!=8.0)fail++;
    /* A decisive round returns to base. */
    changed=cloth_stall_escalation_step(900,600,4.0,12.0,&stall,
                                        &collar,&inflate,&rings);
    if(!changed||collar!=4.0||inflate!=3||rings!=24)fail++;
    fprintf(stderr,
        "  %s: stall escalation ladder escalates and de-escalates\n",
        fail?"FAIL":"ok");
    return fail;
}

static int cloth_band_frontier_selftest(void) {
    const int H=6,W=7;const size_t nv=(size_t)H*(size_t)W;
    uint8_t movable[42];memset(movable,1,sizeof movable);
    for(int r=0;r<H;r++)movable[(size_t)r*(size_t)W+3]=0;
    ClothTurnContact contact={2,(size_t)4*(size_t)(W-1)+4};
    int fail=0;
    for(int rings=1;rings<=4;rings++) {
        size_t got_count=0;
        uint8_t *got=cloth_active_contact_band(
            nv,&contact,1,movable,H,W,rings,&got_count);
        uint8_t want[42];memset(want,0,sizeof want);
        size_t vertex[4];
        cloth_cell_vertices(contact.cell_lo,W,vertex);
        for(int k=0;k<4;k++)if(vertex[k]<nv)want[vertex[k]]=1;
        cloth_cell_vertices(contact.cell_hi,W,vertex);
        for(int k=0;k<4;k++)if(vertex[k]<nv)want[vertex[k]]=1;
        for(int ring=0;ring<rings;ring++) {
            uint8_t next[42];memcpy(next,want,sizeof next);
            for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
                size_t v=(size_t)r*(size_t)W+(size_t)c;
                if(!want[v])continue;
                if(c>0&&movable[v-1])next[v-1]=1;
                if(c+1<W&&movable[v+1])next[v+1]=1;
                if(r>0&&movable[v-(size_t)W])next[v-(size_t)W]=1;
                if(r+1<H&&movable[v+(size_t)W])next[v+(size_t)W]=1;
            }
            memcpy(want,next,sizeof want);
        }
        size_t want_count=0;for(size_t v=0;v<nv;v++)want_count+=want[v]!=0;
        int bad=!got||memcmp(got,want,nv)!=0||got_count!=want_count;
        free(got);
        if(bad){fail++;fprintf(stderr,
            "  FAIL: frontier band != reference ring dilation at rings=%d\n",
            rings);}
    }
    if(!fail)fprintf(stderr,
        "  ok: frontier band equals reference ring dilation (rings 1..4)\n");
    return fail;
}

static int cloth_smooth_list_selftest(void) {
    const int H=8,W=8;const size_t nv=(size_t)H*(size_t)W;
    uint8_t movable[64],active[64];
    float a[64],b[64];
    for(size_t v=0;v<nv;v++) {
        movable[v]=(v%5u)!=0;
        active[v]=(v%3u)!=0;
        a[v]=b[v]=(float)((double)((v*2654435761u)&1023u)/64.0-8.0);
    }
    for(int sweep=0;sweep<4;sweep++)
        for(int parity=0;parity<2;parity++)
            for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
                if(((r+c)&1)!=parity)continue;
                size_t v=(size_t)r*(size_t)W+(size_t)c;
                if(!movable[v]||!active[v])continue;
                double sum=0.0,degree=0.0;
                if(c>0){sum+=(double)a[v-1];degree+=1.0;}
                if(c+1<W){sum+=(double)a[v+1];degree+=1.0;}
                if(r>0){sum+=(double)a[v-(size_t)W];degree+=1.0;}
                if(r+1<H){sum+=(double)a[v+(size_t)W];degree+=1.0;}
                double target=sum/(degree+0.006);
                a[v]=(float)((1.0-0.70)*(double)a[v]+0.70*target);
            }
    int32_t *red=NULL,*black=NULL;size_t nred=0,nblack=0;
    if(cloth_build_parity_lists(movable,active,H,W,
                                &red,&nred,&black,&nblack)!=0)return 1;
    for(int sweep=0;sweep<4;sweep++)
        cloth_smooth_radial_displacement(b,red,nred,black,nblack,H,W,
                                         0.006,0.70);
    free(red);free(black);
    int fail=memcmp(a,b,sizeof a)!=0;
    fprintf(stderr,
        "  %s: list smoother is bitwise equal to the full-grid sweep\n",
        fail?"FAIL":"ok");
    return fail;
}

static int orientation_baseline_selftest(void) {
    int fail=0;
    enum{OB_NV=6,OB_NF=4};
    /* 2x3 lattice; v2 duplicates v1 so face (1,2,5) is degenerate at input. */
    float anchor[OB_NV*3]={
        0,0,0, 0,0,1, 0,0,1,
        0,1,0, 0,1,1, 0,1,2};
    const int32_t faces[OB_NF*3]={0,1,4, 0,4,3, 1,2,5, 1,5,4};
    uint8_t input_bad[OB_NF];
    size_t base_bad=post_snap_orientation_baseline(
        anchor,OB_NV,faces,OB_NF,input_bad);
    int a=base_bad==1&&input_bad[2]==1&&
          input_bad[0]==0&&input_bad[1]==0&&input_bad[3]==0;
    fprintf(stderr,
        "  %s: orientation baseline flags exactly the degenerate input face\n",
        a?"ok":"FAIL");fail+=!a;
    size_t legacy=post_snap_orientation_violations(
        anchor,anchor,OB_NV,faces,OB_NF,NULL,NULL,NULL);
    size_t masked=post_snap_orientation_violations(
        anchor,anchor,OB_NV,faces,OB_NF,input_bad,NULL,NULL);
    int b=legacy==1&&masked==0;
    fprintf(stderr,
        "  %s: identity candidate passes the baseline-relative gate "
        "(legacy=%zu masked=%zu)\n",b?"ok":"FAIL",legacy,masked);fail+=!b;
    float flipped[OB_NV*3];memcpy(flipped,anchor,sizeof flipped);
    flipped[4*3+1]=-1.0f;
    size_t flip_masked=post_snap_orientation_violations(
        flipped,anchor,OB_NV,faces,OB_NF,input_bad,NULL,NULL);
    int c=flip_masked>=1;
    fprintf(stderr,
        "  %s: genuine candidate flip still rejected under the baseline "
        "mask (%zu)\n",c?"ok":"FAIL",flip_masked);fail+=!c;
    return fail;
}

/* ARAP is deliberately a separate transaction after CT snapping.  Every current
 * vertex is a soft data anchor, irrespective of its pre-snap fitted/fill label.
 * We solve one full candidate, then repeatedly roll back only vertices on faces
 * whose exact conflict degree grew, feathering that mask through four grid
 * rings.  A global half-step is only the fail-safe after local repair stalls.
 * This lets fairing remove old collisions but never manufacture a new fold or
 * exchange neighbouring windings. */
static int post_snap_arap_relax(float *verts,size_t nv,const int32_t *faces,
                                size_t nf,const float *uv,int H,int W,
                                const QuadStripArapOpts *opts,
                                PostSnapArapTxn *out) {
    PostSnapArapTxn s;memset(&s,0,sizeof s);s.attempted=1;
    if(out)*out=s;
    if(!verts||!faces||!uv||!opts||H<2||W<2||nv!=(size_t)H*(size_t)W)
        return -1;
    float *anchor=(float *)malloc(nv*3*sizeof *anchor);
    float *full=(float *)malloc(nv*3*sizeof *full);
    float *motion_scale=(float *)malloc(nv*sizeof *motion_scale);
    float *motion_tmp=(float *)malloc(nv*sizeof *motion_tmp);
    size_t *degree_before=(size_t *)calloc(nf,sizeof *degree_before);
    size_t *degree_candidate=(size_t *)calloc(nf,sizeof *degree_candidate);
    uint8_t *all_source=(uint8_t *)calloc(nv,1);
    uint8_t *orientation_input_bad=(uint8_t *)malloc(nf);
    if(!anchor||!full||!motion_scale||!motion_tmp||!degree_before||
       !degree_candidate||!all_source||!orientation_input_bad){
        free(anchor);free(full);free(motion_scale);free(motion_tmp);
        free(degree_before);free(degree_candidate);free(all_source);
        free(orientation_input_bad);return -1;
    }
    memcpy(anchor,verts,nv*3*sizeof *anchor);
    {
        size_t input_bad_faces=post_snap_orientation_baseline(
            anchor,nv,faces,nf,orientation_input_bad);
        if(input_bad_faces)fprintf(stderr,
            "[post-arap] input orientation-bad faces=%zu excluded from "
            "candidate preflight\n",input_bad_faces);
    }
    IntersectionCleanupStats before={0};
    if(snap_intersection_audit_degrees(anchor,nv,faces,nf,degree_before,&before)!=0){
        free(anchor);free(full);free(motion_scale);free(motion_tmp);
        free(degree_before);free(degree_candidate);free(all_source);
        free(orientation_input_bad);return -1;
    }
    s.input_intersections=before.conflicts;s.input_overlap=before.overlap_pairs;
    s.input_stab=before.stab_pairs;s.input_fold=before.fold_pairs;
    s.input_edge_log_rms=post_snap_edge_log_rms(anchor,uv,H,W);

    QuadStripArapOpts o=*opts;
    o.fixed_vertices=NULL;o.repair_weight=NULL;
    s.solver_code=QuadStrip_metric_arap(verts,uv,all_source,H,W,&o,&s.solver);
    if(s.solver_code!=0){
        memcpy(verts,anchor,nv*3*sizeof *verts);
        if(out)*out=s;
        free(anchor);free(full);free(motion_scale);free(motion_tmp);
        free(degree_before);free(degree_candidate);free(all_source);
        free(orientation_input_bad);return 0;
    }
    memcpy(full,verts,nv*3*sizeof *full);

    for(size_t i=0;i<nv;i++)motion_scale[i]=1.0f;
    int accepted=0;IntersectionCleanupStats accepted_ix=before;
    double accepted_edge=s.input_edge_log_rms,accepted_scale=0.0;
    for(int attempt=0;attempt<=12;attempt++){
        size_t zeroed=0;double scale_mean=1.0;
        if(attempt>0){
            size_t marked=0;
            for(size_t f=0;f<nf;f++)if(degree_candidate[f]>degree_before[f])
                for(int k=0;k<3;k++){
                    int32_t vi=faces[f*3+(size_t)k];
                    if(vi>=0&&(size_t)vi<nv&&motion_scale[(size_t)vi]>0.0f){
                        motion_scale[(size_t)vi]=0.0f;marked++;
                    }
                }
            size_t orientation_marked=0;
            (void)post_snap_orientation_violations(
                verts,anchor,nv,faces,nf,orientation_input_bad,
                motion_scale,&orientation_marked);
            marked+=orientation_marked;
            if(marked==0||attempt>6)
                for(size_t i=0;i<nv;i++)motion_scale[i]*=0.5f;
            for(int ring=0;ring<4;ring++){
                for(int r=0;r<H;r++)for(int c=0;c<W;c++){
                    size_t i=(size_t)r*(size_t)W+(size_t)c;
                    float a=motion_scale[i],b;
                    if(c>0){b=motion_scale[i-1]+0.25f;if(b<a)a=b;}
                    if(c+1<W){b=motion_scale[i+1]+0.25f;if(b<a)a=b;}
                    if(r>0){b=motion_scale[i-(size_t)W]+0.25f;if(b<a)a=b;}
                    if(r+1<H){b=motion_scale[i+(size_t)W]+0.25f;if(b<a)a=b;}
                    motion_tmp[i]=a<1.0f?a:1.0f;
                }
                memcpy(motion_scale,motion_tmp,nv*sizeof *motion_scale);
            }
            double sum=0.0;
            for(size_t i=0;i<nv;i++){
                double a=(double)motion_scale[i];sum+=a;
                if(!(a>0.0))zeroed++;
                for(int k=0;k<3;k++)
                    verts[i*3+(size_t)k]=(float)((double)anchor[i*3+(size_t)k]+
                        a*((double)full[i*3+(size_t)k]-anchor[i*3+(size_t)k]));
            }
            scale_mean=nv?sum/(double)nv:0.0;
            if(zeroed>s.locally_rolled_back_vertices)
                s.locally_rolled_back_vertices=zeroed;
            if(o.preserve_axial)
                for(size_t i=0;i<nv;i++)verts[i*3]=anchor[i*3];
        }
        size_t orientation_bad=post_snap_orientation_violations(
            verts,anchor,nv,faces,nf,orientation_input_bad,NULL,NULL);
        IntersectionCleanupStats ix={0};s.trial_audits++;
        memset(degree_candidate,0,nf*sizeof *degree_candidate);
        if(snap_intersection_audit_degrees(verts,nv,faces,nf,
                                           degree_candidate,&ix)!=0){
            memcpy(verts,anchor,nv*3*sizeof *verts);
            free(anchor);free(full);free(motion_scale);free(motion_tmp);
            free(degree_before);free(degree_candidate);free(all_source);
            free(orientation_input_bad);return -1;
        }
        if(orientation_bad)s.orientation_rejects++;
        double edge=post_snap_edge_log_rms(verts,uv,H,W);
        int geometry_safe=ix.conflicts<=before.conflicts&&
                          ix.fold_pairs<=before.fold_pairs&&orientation_bad==0;
        int progress=ix.conflicts<before.conflicts||edge+1e-9<s.input_edge_log_rms;
        if(o.verbose)fprintf(stderr,
            "[post-arap] %s attempt %d: retained %.6f, local zero %zu, "
            "edge-log RMS %.6g->%.6g, conflicts %zu->%zu, folds %zu->%zu, "
            "orientation violations %zu\n",
            geometry_safe&&progress?"accept":"reject",attempt,scale_mean,zeroed,
            s.input_edge_log_rms,edge,before.conflicts,ix.conflicts,
            before.fold_pairs,ix.fold_pairs,orientation_bad);
        if(geometry_safe&&progress){
            accepted=1;accepted_ix=ix;accepted_edge=edge;
            accepted_scale=scale_mean;break;
        }
    }
    if(!accepted){
        memcpy(verts,anchor,nv*3*sizeof *verts);
        s.selected_edge_log_rms=s.input_edge_log_rms;
        s.output_intersections=s.input_intersections;s.output_overlap=s.input_overlap;
        s.output_stab=s.input_stab;s.output_fold=s.input_fold;
    }else{
        s.accepted=1;s.selected_scale=accepted_scale;
        s.selected_edge_log_rms=accepted_edge;
        s.output_intersections=accepted_ix.conflicts;
        s.output_overlap=accepted_ix.overlap_pairs;
        s.output_stab=accepted_ix.stab_pairs;s.output_fold=accepted_ix.fold_pairs;
        double ss=0.0,mx=0.0,zmx=0.0;
        for(size_t i=0;i<nv;i++){
            double d2=0.0;
            for(int k=0;k<3;k++){
                double d=(double)verts[i*3+(size_t)k]-anchor[i*3+(size_t)k];
                d2+=d*d;
            }
            double d=sqrt(d2);ss+=d2;if(d>mx)mx=d;
            double zd=fabs((double)verts[i*3]-anchor[i*3]);if(zd>zmx)zmx=zd;
        }
        s.displacement_rms=sqrt(ss/(double)nv);s.displacement_max=mx;
        s.axial_drift_max=zmx;
    }
    if(out)*out=s;
    free(anchor);free(full);free(motion_scale);free(motion_tmp);
    free(degree_before);free(degree_candidate);free(all_source);
    free(orientation_input_bad);return 0;
}

static size_t snap_bake_dark_total(const QuadStripSnapStats *s) {
    return s->output_bake_fitted_dark+s->output_bake_filled_dark;
}

static size_t snap_bake_samples_total(const QuadStripSnapStats *s) {
    return s->output_bake_fitted_samples+s->output_bake_filled_samples;
}

static size_t snap_dark_slack(size_t dark) {
    size_t proportional=(dark+999u)/1000u; /* 0.1%, plus an 8-pixel floor */
    return proportional>8u?proportional:8u;
}

/* Pass one used to be a mandatory checkpoint: only folds could restore the
 * initializer.  Gate it against statistics measured inside the same run before
 * any motion, including fixed-window triangle-interior RAW samples.  Absolute
 * Sander thresholds are diagnostics, not authority: a difficult whole-scroll
 * initializer must still produce inspectable geometry when the transaction is
 * finite, orientation-safe, and non-regressing relative to its own input. */
static int snap_first_pass_accept(const QuadStripSnapStats *s,
                                  char *reason,size_t reason_cap) {
    const double eps=1e-12;
    int dark_progress=0;
    if(!s||!isfinite(s->sander_p95)||!isfinite(s->sander_p99)||
       s->flipped||s->degenerate){
        snprintf(reason,reason_cap,"invalid first-pass geometry");return 0;
    }
    if(!s->intersection_audit_complete) {
        snprintf(reason,reason_cap,"missing exact intersection audit");return 0;
    }
    if(s->output_intersections>s->input_intersections) {
        snprintf(reason,reason_cap,"exact intersections grew %zu->%zu",
                 s->input_intersections,s->output_intersections);return 0;
    }
    if(s->output_fold>s->input_fold) {
        snprintf(reason,reason_cap,"exact fold conflicts grew %zu->%zu",
                 s->input_fold,s->output_fold);return 0;
    }
    size_t in_samples=s->input_bake_fitted_samples+s->input_bake_filled_samples;
    size_t out_samples=snap_bake_samples_total(s);
    if(in_samples&&out_samples){
        size_t before=s->input_bake_fitted_dark+s->input_bake_filled_dark;
        size_t after=snap_bake_dark_total(s);
        if(s->output_bake_missing>s->input_bake_missing){
            snprintf(reason,reason_cap,"first-pass bake missing grew %zu->%zu",
                s->input_bake_missing,s->output_bake_missing);return 0;
        }
        if(after>before+snap_dark_slack(before)){
            snprintf(reason,reason_cap,"first-pass fixed-window dark grew %zu->%zu",
                before,after);return 0;
        }
        size_t fit_slack=snap_dark_slack(s->input_bake_fitted_dark)*2u;
        if(s->output_bake_fitted_dark>s->input_bake_fitted_dark+fit_slack &&
           after>=before){
            snprintf(reason,reason_cap,"first-pass fitted dark grew %zu->%zu",
                s->input_bake_fitted_dark,s->output_bake_fitted_dark);return 0;
        }
        dark_progress=after+snap_dark_slack(before)<before;
    }
    /* A directly measured raster recovery may earn a bounded first-recovery
     * envelope.  Without that exception, the old
     * mandatory first pass could roll back a huge dark-material reduction over
     * a few hundredths of central Sander.  The recovery envelope was bracketed
     * against the complete 20-component 4x candidate fleet (all good V88-like
     * recoveries stay inside +.20/+.35 p95/p99 and +.025% severe tail).  No
     * darkness evidence can excuse a fold, missing RAW samples, or leaving it. */
    {
        double limit=dark_progress?SNAP_TXN_FIRST_DARK_P95_STEP:SNAP_TXN_P95_STEP;
        if(s->sander_p95>s->input_sander_p95+limit+eps){
            snprintf(reason,reason_cap,"first-pass Sander p95 %+0.3f > +%.3f",
                s->sander_p95-s->input_sander_p95,limit);return 0;
        }
    }
    {
        double limit=dark_progress?SNAP_TXN_FIRST_DARK_P99_STEP:SNAP_TXN_P99_STEP;
        if(s->sander_p99>s->input_sander_p99+limit+eps){
            snprintf(reason,reason_cap,"first-pass Sander p99 %+0.3f > +%.3f",
                s->sander_p99-s->input_sander_p99,limit);return 0;
        }
    }
    {
        double factor=dark_progress?3.0:SNAP_TXN_OVER4_FACTOR;
        double step=dark_progress?SNAP_TXN_FIRST_DARK_OVER4_STEP:SNAP_TXN_OVER4_STEP;
        double limit=fmax(s->input_over4_fraction*factor,
                          s->input_over4_fraction+step);
        if(s->over4_fraction>limit+eps){
            snprintf(reason,reason_cap,"first-pass >=4x %.4f%% exceeds %.4f%%",
                100.0*s->over4_fraction,100.0*limit);return 0;
        }
    }
    snprintf(reason,reason_cap,dark_progress?
        "fixed-window recovery within initializer safety envelope":
        "within initializer geometry/bake envelope");
    return 1;
}

/* Pass one is an unbiased measurement.  Only a statistically supported seam
 * whose applied signed-depth jump is both 25% and .03 vox above the larger
 * fitted/filled interior RMS earns stronger mixed-edge coupling on later
 * proximal passes.  V52 brackets this cleanly: groups 5/15 are 1.53x/1.70x;
 * the next component is 1.17x. */
static int snap_quilt_seam_needs_refine(const QuadStripSnapStats *s,
                                        double *out_ratio,
                                        double *out_excess) {
    if (!s) return 0;
    const QuadStripSnapRegionEdgeStats *r = &s->quilt_applied_regions;
    double interior = fmax(r->fitted_fitted_rms, r->filled_filled_rms);
    double ratio = interior > 0.0 ? r->fitted_filled_rms / interior : 0.0;
    double excess = r->fitted_filled_rms - interior;
    if (out_ratio) *out_ratio = ratio;
    if (out_excess) *out_excess = excess;
    return r->fitted_filled_edges >= SNAP_QUILT_SEAM_MIN_EDGES &&
           interior > 0.0 && ratio > SNAP_QUILT_SEAM_TRIGGER_RATIO &&
           excess > SNAP_QUILT_SEAM_TRIGGER_STEP;
}

/* A later proximal pass is a speculative refit, not an instruction to keep
 * moving forever.  The pass may include a structured remesh before it starts,
 * so compare UV->scroll conditioning against both the last accepted state and
 * the first-pass envelope.  Small geometry costs are allowed for new
 * recto evidence; folds, a growing severe tail, and cumulative drift are not.
 * The thresholds below were bracketed by isolated A/Bs: track 1's pass two must
 * roll back, while track 14's useful active-set refinement must survive. */
static int snap_transaction_accept(const QuadStripSnapStats *prev,
                                   const QuadStripSnapStats *first,
                                   const QuadStripSnapStats *cand,
                                   int pass_number,
                                   char *reason, size_t reason_cap) {
    const double eps = 1e-12;
    if (!prev || !first || !cand) {
        snprintf(reason, reason_cap, "missing transaction statistics"); return 0;
    }
    if(!cand->intersection_audit_complete||!prev->intersection_audit_complete) {
        snprintf(reason,reason_cap,"missing exact intersection audit");return 0;
    }
    if(cand->output_intersections>prev->output_intersections) {
        snprintf(reason,reason_cap,"exact intersections grew %zu->%zu",
                 prev->output_intersections,cand->output_intersections);return 0;
    }
    if(cand->output_fold>prev->output_fold) {
        snprintf(reason,reason_cap,"exact fold conflicts grew %zu->%zu",
                 prev->output_fold,cand->output_fold);return 0;
    }
    if (!isfinite(cand->sander_p95) || !isfinite(cand->sander_p99) ||
        !isfinite(cand->over4_fraction)) {
        snprintf(reason, reason_cap, "non-finite audit statistic"); return 0;
    }
    if (cand->flipped > 0) {
        snprintf(reason, reason_cap, "%zu flipped triangles", cand->flipped); return 0;
    }
    if (cand->degenerate > 0) {
        snprintf(reason, reason_cap, "%zu degenerate triangles", cand->degenerate); return 0;
    }
    int measured_dark_progress=0;
    if(pass_number==2&&snap_bake_samples_total(prev)>0&&
       snap_bake_samples_total(cand)>0&&
       cand->output_bake_missing<=prev->output_bake_missing) {
        size_t pd=snap_bake_dark_total(prev),cd=snap_bake_dark_total(cand);
        double gain=pd?(double)(pd>cd?pd-cd:0)/(double)pd:0.0;
        measured_dark_progress=cd<pd&&gain>=SNAP_TXN_REFIT_DARK_MIN_FRACTION;
    }
    double p95_step=measured_dark_progress?
        SNAP_TXN_REFIT_DARK_P95_STEP:SNAP_TXN_P95_STEP;
    double p99_step=measured_dark_progress?
        SNAP_TXN_REFIT_DARK_P99_STEP:SNAP_TXN_P99_STEP;
    double p95_envelope=measured_dark_progress?
        SNAP_TXN_REFIT_DARK_P95_ENVELOPE:SNAP_TXN_P95_ENVELOPE;
    double p99_envelope=measured_dark_progress?
        SNAP_TXN_REFIT_DARK_P99_ENVELOPE:SNAP_TXN_P99_ENVELOPE;
    if (cand->sander_p95 > prev->sander_p95 + p95_step + eps) {
        snprintf(reason, reason_cap, "Sander p95 step %.3f > +%.3f",
                 cand->sander_p95-prev->sander_p95, p95_step); return 0;
    }
    if (cand->sander_p95 > first->sander_p95 + p95_envelope + eps) {
        snprintf(reason, reason_cap, "Sander p95 left first-pass envelope by %.3f",
                 cand->sander_p95-first->sander_p95); return 0;
    }
    {
        int central_improved = cand->sander_p95 + SNAP_TXN_P95_STEP < prev->sander_p95;
        if (!central_improved && cand->sander_p99 > prev->sander_p99 + p99_step + eps) {
            snprintf(reason, reason_cap, "Sander p99 step %.3f > +%.3f",
                     cand->sander_p99-prev->sander_p99, p99_step); return 0;
        }
        double tail_limit = fmax(prev->over4_fraction*SNAP_TXN_OVER4_FACTOR,
                                 prev->over4_fraction+SNAP_TXN_OVER4_STEP);
        if (!central_improved && cand->over4_fraction > tail_limit + eps) {
            snprintf(reason, reason_cap, ">=4x tail %.4f%% exceeds %.4f%%",
                     100.0*cand->over4_fraction, 100.0*tail_limit); return 0;
        }
    }
    if (cand->sander_p99 > first->sander_p99 + p99_envelope + eps) {
        snprintf(reason, reason_cap, "Sander p99 left first-pass envelope by %.3f",
                 cand->sander_p99-first->sander_p99); return 0;
    }
    {
        double tail_limit = fmax(first->over4_fraction*2.0,
                                 first->over4_fraction+SNAP_TXN_OVER4_ENVELOPE);
        if (cand->over4_fraction > tail_limit + eps) {
            snprintf(reason, reason_cap, ">=4x tail left first-pass envelope (%.4f%% > %.4f%%)",
                     100.0*cand->over4_fraction, 100.0*tail_limit); return 0;
        }
    }
    if(snap_bake_samples_total(prev)>0&&snap_bake_samples_total(cand)>0){
        size_t pd=snap_bake_dark_total(prev),cd=snap_bake_dark_total(cand);
        if(cand->output_bake_missing>prev->output_bake_missing){
            snprintf(reason,reason_cap,"fixed-window bake missing grew %zu->%zu",
                prev->output_bake_missing,cand->output_bake_missing);return 0;
        }
        if(cd>pd+snap_dark_slack(pd)){
            snprintf(reason,reason_cap,"fixed-window dark grew %zu->%zu",pd,cd);
            return 0;
        }
        if(snap_bake_samples_total(first)>0){
            size_t fd=snap_bake_dark_total(first);
            if(cd>fd+2u*snap_dark_slack(fd)){
                snprintf(reason,reason_cap,
                    "fixed-window dark left first-pass envelope (%zu>%zu)",cd,fd);
                return 0;
            }
        }
    }

    /* Pass two is the useful active-set refit: the isolated track-14 A/B
     * demonstrates that it can earn a small, bounded geometry cost.  Once two
     * states have already been accepted, however, merely remaining inside the
     * broad safety envelope is not convergence.  Mature passes must stop
     * walking the sheet unless they are Pareto-nonworse in the primary Sander
     * statistics and make measurable geometric or coupled CT/quilt progress.
     * Roughly one severe-tail triangle is treated as histogram noise. */
    if (pass_number >= 3) {
        double tail_eps = cand->sander_samples > 0
                        ? 1.5/(double)cand->sander_samples : eps;
        int mature_dark_progress=0;
        if(snap_bake_samples_total(prev)>0&&snap_bake_samples_total(cand)>0) {
            size_t pd=snap_bake_dark_total(prev),cd=snap_bake_dark_total(cand);
            double gain=pd?(double)(pd>cd?pd-cd:0)/(double)pd:0.0;
            mature_dark_progress=cd<pd&&
                gain>=SNAP_TXN_REFIT_DARK_MIN_FRACTION;
        }
        if (cand->sander_p95 > prev->sander_p95 + SNAP_TXN_MATURE_P95_EPS + eps) {
            snprintf(reason, reason_cap, "mature pass worsened Sander p95 by %.3f",
                     cand->sander_p95-prev->sander_p95); return 0;
        }
        if (cand->sander_p99 > prev->sander_p99 + SNAP_TXN_MATURE_P99_EPS + eps) {
            snprintf(reason, reason_cap, "mature pass worsened Sander p99 by %.3f",
                     cand->sander_p99-prev->sander_p99); return 0;
        }
        if (!mature_dark_progress&&
            cand->over4_fraction > prev->over4_fraction + tail_eps + eps) {
            snprintf(reason, reason_cap, "mature pass grew >=4x tail by %.0f triangles",
                     (cand->over4_fraction-prev->over4_fraction)
                         *(double)cand->sander_samples);
            return 0;
        }

        int ct_prev = prev->filled_ct_p25 + prev->filled_ct_p50 + prev->filled_ct_p75;
        int ct_cand = cand->filled_ct_p25 + cand->filled_ct_p50 + cand->filled_ct_p75;
        int ct_nondecreasing = cand->filled_ct_p25 >= prev->filled_ct_p25 &&
                               cand->filled_ct_p50 >= prev->filled_ct_p50 &&
                               cand->filled_ct_p75 >= prev->filled_ct_p75;
        int ct_progress = ct_nondecreasing &&
                          ct_cand >= ct_prev + SNAP_TXN_MATURE_CT_GAIN;
        int geometry_progress =
            cand->sander_p95 + SNAP_TXN_MATURE_P95_EPS <= prev->sander_p95 ||
            cand->sander_p99 + 0.010 <= prev->sander_p99 ||
            cand->over4_fraction + tail_eps < prev->over4_fraction;
        int evidence_nondecreasing = ct_cand >= ct_prev-1;
        int dark_progress=mature_dark_progress||(
            snap_bake_samples_total(prev)>0&&snap_bake_samples_total(cand)>0&&
            snap_bake_dark_total(cand)+snap_dark_slack(snap_bake_dark_total(prev))<
                snap_bake_dark_total(prev));
        if (!((geometry_progress && evidence_nondecreasing) ||
              ct_progress || dark_progress)) {
            snprintf(reason, reason_cap,
                     "mature pass made no geometry/CT progress (CT %+d)",
                     ct_cand-ct_prev);
            return 0;
        }
        snprintf(reason, reason_cap, mature_dark_progress
            ?"mature measured-dark progress within topology/central-geometry envelope"
            :"mature Pareto progress");
        return 1;
    }
    snprintf(reason, reason_cap, measured_dark_progress?
        "measured dark recovery within refit safety envelope":
        "within geometry/CT envelope");
    return 1;
}

static int snap_transaction_selftest(void) {
    QuadStripSnapStats a={0},bad={0},good={0}; char why[160]; int failures=0;
    a.intersection_audit_complete=1;
    a.sander_p95=1.249;a.sander_p99=1.585;a.over4_fraction=0.00009;
    a.sander_samples=10000;
    a.quilt_edge_rms=0.177;a.quilt_edges=100;
    bad=a;bad.sander_p95=1.290;bad.sander_p99=1.551;
    bad.over4_fraction=0.00023;bad.quilt_edge_rms=0.172;
    if(snap_transaction_accept(&a,&a,&bad,2,why,sizeof why))failures++;
    fprintf(stderr,"  %s: transaction rejects track-1 geometry regression (%s)\n",
            failures?"FAIL":"ok",why);
    a.sander_p95=1.120;a.sander_p99=1.235;a.over4_fraction=0.0;
    a.quilt_edge_rms=0.174;a.quilt_edges=100;
    good=a;good.sander_p95=1.133;good.sander_p99=1.276;
    good.over4_fraction=0.00002;good.quilt_edge_rms=0.190;
    {
        int ok=snap_transaction_accept(&a,&a,&good,2,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: transaction accepts track-14 recto refinement (%s)\n",
                ok?"ok":"FAIL",why);
    }
    good.flipped=1;
    {
        int ok=!snap_transaction_accept(&a,&a,&good,2,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: transaction never accepts a fold (%s)\n",
                ok?"ok":"FAIL",why);
    }
    good.flipped=0;
    {
        QuadStripSnapStats prev={0},refit={0};
        prev.intersection_audit_complete=1;
        prev.output_intersections=2452;
        prev.sander_p95=1.262;prev.sander_p99=1.453;
        prev.over4_fraction=0.0000119;prev.sander_samples=1261872;
        prev.output_bake_filled_samples=1261872;
        prev.output_bake_filled_dark=104901;
        refit=prev;refit.output_intersections=2354;
        refit.sander_p95=1.304;refit.sander_p99=1.637;
        refit.over4_fraction=0.000052;
        refit.output_bake_filled_dark=85322;
        int ok=snap_transaction_accept(&prev,&prev,&refit,2,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: measured pass-2 dark recovery earns bounded geometry cost (%s)\n",
                ok?"ok":"FAIL",why);
        refit.output_bake_filled_dark=100000; /* <5%: not enough evidence */
        ok=!snap_transaction_accept(&prev,&prev,&refit,2,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: weak pass-2 dark change cannot excuse geometry drift (%s)\n",
                ok?"ok":"FAIL",why);
    }
    a.sander_p95=1.290;a.sander_p99=1.437;a.over4_fraction=0.0;
    a.quilt_edge_rms=0.149;a.quilt_edges=100;a.sander_samples=12480;a.quads=6240;
    a.filled_ct_p25=48;a.filled_ct_p50=77;a.filled_ct_p75=125;
    good=a;good.sander_p99=1.422;good.quilt_edge_rms=0.134;
    good.filled_ct_p25=49;
    {
        int ok=snap_transaction_accept(&a,&a,&good,3,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: mature transaction accepts Pareto progress (%s)\n",
                ok?"ok":"FAIL",why);
    }
    good=a;good.sander_p99=1.451;good.quilt_edge_rms=0.140;
    good.filled_ct_p25=50;good.filled_ct_p50=80;
    {
        int ok=!snap_transaction_accept(&a,&a,&good,3,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: mature transaction rejects renewed drift (%s)\n",
                 ok?"ok":"FAIL",why);
    }
    {
        QuadStripSnapStats prev={0},cand={0};
        prev.intersection_audit_complete=1;prev.output_intersections=120;
        prev.output_fold=4;prev.sander_p95=4.433;prev.sander_p99=5.688;
        prev.over4_fraction=0.059017;prev.sander_samples=708914;
        prev.output_bake_fitted_samples=272045;
        prev.output_bake_filled_samples=85331;
        prev.output_bake_fitted_dark=250643;
        prev.output_bake_filled_dark=521424;
        prev.filled_ct_p25=49;prev.filled_ct_p50=75;prev.filled_ct_p75=124;
        cand=prev;cand.output_intersections=110;cand.output_fold=3;
        cand.over4_fraction=prev.over4_fraction+30.0/(double)cand.sander_samples;
        cand.output_bake_fitted_dark=214312;cand.output_bake_filled_dark=507020;
        int ok=snap_transaction_accept(&prev,&prev,&cand,3,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: mature measured dark recovery can absorb a tiny severe-tail histogram cost (%s)\n",
                ok?"ok":"FAIL",why);
    }
    {
        QuadStripSnapStats seam={0}; double ratio=0.0,excess=0.0;
        seam.quilt_applied_regions.fitted_fitted_rms=0.21;
        seam.quilt_applied_regions.filled_filled_rms=0.24;
        seam.quilt_applied_regions.fitted_filled_rms=0.37;
        seam.quilt_applied_regions.fitted_filled_edges=1097;
        int ok=snap_quilt_seam_needs_refine(&seam,&ratio,&excess);
        if(!ok)failures++;
        fprintf(stderr,"  %s: measured quilt seam outlier activates refinement "
                       "(ratio %.2f, excess %.3f)\n",ok?"ok":"FAIL",ratio,excess);
        seam.quilt_applied_regions.fitted_filled_rms=0.27;
        ok=!snap_quilt_seam_needs_refine(&seam,&ratio,&excess);
        if(!ok)failures++;
        fprintf(stderr,"  %s: ordinary quilt variation stays uniformly coupled "
                       "(ratio %.2f, excess %.3f)\n",ok?"ok":"FAIL",ratio,excess);
    }
    {
        QuadStripSnapStats first={0};
        first.intersection_audit_complete=1;
        first.input_sander_p95=first.sander_p95=1.10;
        first.input_sander_p99=first.sander_p99=1.20;
        first.sander_samples=1000;
        first.input_bake_fitted_samples=first.output_bake_fitted_samples=1000;
        first.input_bake_fitted_dark=100;
        first.output_bake_fitted_dark=110;
        int ok=!snap_first_pass_accept(&first,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: first transaction rejects fixed-window dark growth (%s)\n",
                ok?"ok":"FAIL",why);
        first.output_bake_fitted_dark=90;
        first.sander_p95=1.14;
        ok=snap_first_pass_accept(&first,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: first transaction accepts bake improvement (%s)\n",
                ok?"ok":"FAIL",why);
        /* Absolute quality is reported downstream, but must never erase a
         * relative, exact-audited improvement on a difficult full sheet. */
        first.input_sander_p95=3.10;first.sander_p95=3.00;
        first.input_sander_p99=5.30;first.sander_p99=5.20;
        first.input_over4_fraction=0.032;first.over4_fraction=0.031;
        first.input_intersections=4;first.output_intersections=3;
        first.input_fold=2;first.output_fold=2;
        ok=snap_first_pass_accept(&first,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: first transaction emits an exact-audited "
                       "improvement above diagnostic Sander limits (%s)\n",
                ok?"ok":"FAIL",why);
        first.input_sander_p95=1.10;first.sander_p95=1.14;
        first.input_sander_p99=first.sander_p99=1.20;
        first.input_over4_fraction=first.over4_fraction=0.0;
        first.input_intersections=first.output_intersections=0;
        first.input_fold=first.output_fold=0;
        first.output_bake_fitted_dark=100;
        ok=!snap_first_pass_accept(&first,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: first transaction requires evidence for geometry cost (%s)\n",
                ok?"ok":"FAIL",why);
        first.output_bake_fitted_dark=90;
        first.sander_p95=1.31;
        ok=!snap_first_pass_accept(&first,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: first recovery cannot leave its hard geometry envelope (%s)\n",
                ok?"ok":"FAIL",why);
        first.sander_p95=1.10;
        first.output_bake_fitted_dark=90;
        first.input_intersections=4;first.output_intersections=5;
        ok=!snap_first_pass_accept(&first,why,sizeof why);
        if(!ok)failures++;
        fprintf(stderr,"  %s: first transaction rejects exact intersection growth (%s)\n",
                ok?"ok":"FAIL",why);
    }
    return failures;
}

static int initializer_quarantine_selftest(void) {
    const int H=9,W=17;
    size_t HW=(size_t)H*(size_t)W;
    Arena_T arena=Arena_new();
    uint8_t *domain=ARENA_ALLOC(arena,HW);
    uint8_t *trusted=ARENA_ALLOC(arena,HW);
    float *field=ARENA_ALLOC(arena,HW*3*sizeof *field);
    int failures=0;
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t k=(size_t)r*(size_t)W+(size_t)c;
        domain[k]=1;                         /* dense fitted initializer */
        trusted[k]=(uint8_t)((r&1)==0);      /* no all-supported quad */
        double x=c<4?2.0*(double)c:
                 c<=8?8.0+6.0*(double)(c-4):
                      32.0+2.0*(double)(c-8);
        field[k*3]=(float)r;
        field[k*3+1]=0.0f;
        field[k*3+2]=(float)x;
    }
    InitQuarantineStats s={0};
    int rc=quarantine_initializer_geometry(
            arena,domain,trusted,field,H,W,2.0,1.5,4,0.20,0,&s);
    int boundary_removed=!domain[5]&&!domain[(size_t)(H-1)*W+5];
    int unaffected=domain[0]&&domain[(size_t)(H-1)*W+(W-1)];
    int ok=rc==0&&s.measured_triangles==2u*(size_t)(H-1)*(size_t)(W-1)&&
           s.severe_triangles>0&&s.removed_initializer>0&&boundary_removed&&
           unaffected;
    if(!ok)failures++;
    fprintf(stderr,
            "  %s: dense initializer quarantine sees alternating-row "
            "distortion and removes boundary support\n",ok?"ok":"FAIL");
    Arena_dispose(&arena);
    return failures;
}

static int phase_gauge_selftest(void) {
    const int H=3,W=10;const double two_pi=6.2831853071795864769;
    size_t HW=(size_t)H*(size_t)W;
    Arena_T arena=Arena_new();int failures=0;
    uint8_t *domain=ARENA_ALLOC(arena,HW),*trusted=ARENA_ALLOC(arena,HW);
    float *field=ARENA_ALLOC(arena,HW*3*sizeof *field);
    float *phase=ARENA_ALLOC(arena,HW*sizeof *phase);
    for(int r=0;r<H;r++)for(int c=0;c<W;c++) {
        size_t k=(size_t)r*(size_t)W+(size_t)c;
        double a=3.0+0.05*(double)c;
        domain[k]=trusted[k]=1;
        field[k*3]=(float)r;field[k*3+1]=(float)(10.0*cos(a));
        field[k*3+2]=(float)(10.0*sin(a));phase[k]=(float)a;
    }
    PhaseGateStats s={0};
    int rc=quarantine_phase_gauge(
            arena,domain,trusted,field,phase,H,W,0.0,0.0,&s);
    int seam_ok=rc==0&&s.removed==0&&s.discontinuous_edges==0;
    if(!seam_ok)failures++;
    fprintf(stderr,"  %s: phase gauge accepts a continuous atan2 branch crossing\n",
            seam_ok?"ok":"FAIL");

    for(size_t k=0;k<HW;k++)domain[k]=trusted[k]=1;
    phase[(size_t)1*W+5]+=(float)two_pi;
    memset(&s,0,sizeof s);
    rc=quarantine_phase_gauge(
            arena,domain,trusted,field,phase,H,W,0.0,0.0,&s);
    int jump_ok=rc==0&&s.discontinuous_edges==4&&s.removed==5&&
                !domain[(size_t)1*W+5];
    if(!jump_ok)failures++;
    fprintf(stderr,"  %s: phase gauge removes a false whole-turn island and rim\n",
            jump_ok?"ok":"FAIL");
    Arena_dispose(&arena);
    return failures;
}

static void json_write_string(FILE *f,const char *s) {
    fputc('"',f);
    for(;s&&*s;s++){
        unsigned char ch=(unsigned char)*s;
        if(ch=='"'||ch=='\\'){fputc('\\',f);fputc(ch,f);}
        else if(ch=='\n')fputs("\\n",f);
        else if(ch=='\r')fputs("\\r",f);
        else if(ch=='\t')fputs("\\t",f);
        else if(ch<0x20)fprintf(f,"\\u%04x",(unsigned)ch);
        else fputc(ch,f);
    }
    fputc('"',f);
}

static void snap_report_region_edge_stats(
        FILE *f, const QuadStripSnapRegionEdgeStats *s) {
    fprintf(f,
            "{\"all\":{\"rms\":%.9g,\"edges\":%zu},"
            "\"fitted_fitted\":{\"rms\":%.9g,\"edges\":%zu},"
            "\"filled_filled\":{\"rms\":%.9g,\"edges\":%zu},"
            "\"fitted_filled\":{\"rms\":%.9g,\"edges\":%zu}}",
            s->all_rms,s->all_edges,
            s->fitted_fitted_rms,s->fitted_fitted_edges,
            s->filled_filled_rms,s->filled_filled_edges,
            s->fitted_filled_rms,s->fitted_filled_edges);
}

static void snap_report_entry(FILE *f,int *first,size_t component,long group,
                              long track,int pass,int accepted,const char *reason,
                              const QuadStripSnapStats *s) {
    if(!f)return;
    if(!*first)fputs(",\n",f);*first=0;
    fprintf(f,"  {\"component\":%zu,\"group\":%ld,\"representative_track\":%ld,"
              "\"pass\":%d,\"accepted\":%s,\"reason\":",
            component,group,track,pass,accepted?"true":"false");
    json_write_string(f,reason?reason:"");
    if(s){
        fprintf(f,
          ",\"geometry\":{\"sander_basis\":\"bake_triangles\","
          "\"input_sander_p95\":%.9g,\"input_sander_p99\":%.9g,"
          "\"input_over4_fraction\":%.9g,"
          "\"displacement_rms\":%.9g,\"tangential_p95\":%.9g,"
          "\"tangential_max\":%.9g,\"sander_p95\":%.9g,\"sander_p99\":%.9g,"
          "\"sander_max\":%.9g,\"over2_fraction\":%.9g,\"over4_fraction\":%.9g,"
          "\"samples\":%zu,\"quads\":%zu,\"flipped\":%zu,\"degenerate\":%zu,"
          "\"axial_drift_max\":%.9g,"
          "\"quad_average\":{\"sander_p95\":%.9g,\"sander_p99\":%.9g,"
          "\"sander_max\":%.9g,\"over2_fraction\":%.9g,"
          "\"over4_fraction\":%.9g},"
          "\"intersection_audit\":{\"complete\":%s,"
          "\"input\":{\"conflicts\":%zu,\"overlap\":%zu,\"stab\":%zu,\"fold\":%zu},"
          "\"output\":{\"conflicts\":%zu,\"overlap\":%zu,\"stab\":%zu,\"fold\":%zu},"
          "\"round_guard\":{\"audits\":%d,\"backtracks\":%d,"
          "\"round_rollbacks\":%d,\"orientation_rejects\":%d,"
          "\"locally_rolled_back_vertices\":%zu}}},"
          "\"quilt\":{\"offset_basis\":\"world_zyx_vector\","
          "\"raw_ridge\":%d,\"recto_ridge\":%d,"
          "\"band_fallback\":%d,\"band_abstain\":%d,"
          "\"near_zero_ridge\":%d,"
          "\"snapped\":%d,\"nopap\":%d,\"raw_no_candidate\":%d,"
          "\"raw_candidate_applied\":%d,\"quilt_only_applied\":%d,"
          "\"guided_ridge\":%d,\"wide_ridge\":%d,"
          "\"local_contrast_ridge\":%d,\"winding_rejected\":%d,"
          "\"dark_seed_forced_movable\":%d,\"dark_seed_quilt_only\":%d,"
          "\"dark_face_samples\":%zu,\"dark_face_candidates\":%zu,"
          "\"dark_face_outward_candidates\":%zu,\"dark_face_no_candidate\":%zu,"
          "\"adaptive_winding_vertices\":%zu,"
          "\"adaptive_winding_half_width_p05_p50_p95\":[%.9g,%.9g,%.9g],"
          "\"downweighted_candidates\":%d,\"candidate_weight_mean\":%.9g,"
          "\"raw_edge_rms\":%.9g,\"edge_rms\":%.9g,"
          "\"raw_edges\":%zu,\"edges\":%zu,\"seam_smooth\":%.9g},"
          "\"raw_coverage\":{\"expected_chunks\":%zu,\"loaded_chunks\":%zu,"
          "\"missing_chunks\":%zu,\"outside_volume_chunks\":%zu,"
          "\"complete\":%s},"
          "\"ct\":{\"fitted_samples\":%zu,\"filled_samples\":%zu,"
          "\"fitted_p1_p25_p50_p75_p99\":[%d,%d,%d,%d,%d],"
          "\"filled_p1_p25_p50_p75_p99\":[%d,%d,%d,%d,%d]},"
          "\"dark_probe\":{\"threshold\":%.9g,"
          "\"fitted_bake_dark\":%zu,\"fitted_recto_target\":%zu,"
          "\"filled_bake_dark\":%zu,\"filled_recto_target\":%zu}",
          s->input_sander_p95,s->input_sander_p99,s->input_over4_fraction,
          s->displacement_rms,s->tangential_p95,s->tangential_max,
          s->sander_p95,s->sander_p99,s->sander_max,
          s->over2_fraction,s->over4_fraction,s->sander_samples,s->quads,
          s->flipped,s->degenerate,s->axial_drift_max,
          s->quad_sander_p95,s->quad_sander_p99,s->quad_sander_max,
          s->quad_over2_fraction,s->quad_over4_fraction,
          s->intersection_audit_complete?"true":"false",
          s->input_intersections,s->input_overlap,s->input_stab,s->input_fold,
          s->output_intersections,s->output_overlap,s->output_stab,s->output_fold,
          s->intersection_guard_audits,s->intersection_guard_backtracks,
          s->intersection_guard_round_rollbacks,
          s->intersection_guard_orientation_rejects,
          s->intersection_guard_locally_rolled_back_vertices,
          s->raw_ridge,s->recto_ridge,s->band_fallback,s->band_abstain,
          s->near_zero_ridge,s->snapped,s->nopap,
          s->raw_no_candidate,s->raw_candidate_applied,s->quilt_only_applied,
          s->guided_ridge,s->wide_ridge,s->local_contrast_ridge,s->winding_rejected,
          s->dark_seed_forced_movable,s->dark_seed_quilt_only,
          s->dark_face_samples,s->dark_face_candidates,
          s->dark_face_outward_candidates,s->dark_face_no_candidate,
          s->adaptive_winding_vertices,s->adaptive_winding_p05,
          s->adaptive_winding_p50,s->adaptive_winding_p95,
          s->quilt_downweighted,s->quilt_candidate_weight_mean,
          s->quilt_raw_edge_rms,s->quilt_edge_rms,s->quilt_raw_edges,s->quilt_edges,
          s->quilt_seam_smooth,
          s->raw_chunks_expected,s->raw_chunks_loaded,s->raw_chunks_missing,
          s->raw_chunks_outside_volume,
          s->raw_chunks_missing == 0 ? "true" : "false",
          s->fitted_ct_samples,s->filled_ct_samples,
          s->fitted_ct_p1,s->fitted_ct_p25,s->fitted_ct_p50,s->fitted_ct_p75,s->fitted_ct_p99,
          s->filled_ct_p1,s->filled_ct_p25,s->filled_ct_p50,s->filled_ct_p75,s->filled_ct_p99,
          s->dark_probe_threshold,s->fitted_bake_dark,s->fitted_dark_recto_target,
          s->filled_bake_dark,s->filled_dark_recto_target);
        fprintf(f,
          ",\"patch_recovery\":{\"regions\":%zu,\"regions_small\":%zu,"
          "\"border_consensus\":%zu,\"ray_consensus\":%zu,"
          "\"border_selected\":%zu,\"ambiguous_selected\":%zu,"
          "\"wobbly\":%zu,\"no_consensus\":%zu,\"overlap_skipped\":%zu,"
          "\"proposed\":%zu,\"applied\":%zu,\"attenuated\":%zu,\"geometry_rejected\":%zu,"
          "\"vertices_proposed\":%zu,\"vertices_retained\":%zu,"
          "\"filled_dark_before\":%zu,\"filled_dark_after\":%zu,"
          "\"transaction_accepted\":%s,\"retained_scale\":%.9g}",
          s->patch_regions,s->patch_regions_small,s->patch_border_consensus,
          s->patch_ray_consensus,s->patch_border_selected,s->patch_ambiguous_selected,
          s->patch_wobbly,s->patch_no_consensus,
          s->patch_overlap_skipped,s->patch_proposed,s->patch_applied,
          s->patch_attenuated,s->patch_geometry_rejected,
          s->patch_vertices,s->patch_vertices_retained,
          s->patch_dark_before,s->patch_dark_after,
          s->patch_transaction_accepted?"true":"false",s->patch_retained_scale);
        fprintf(f,
          ",\"bake_audit\":{\"basis\":\"fixed_window_triangle_interiors\","
          "\"input\":{\"fitted_samples\":%zu,\"filled_samples\":%zu,"
          "\"fitted_dark\":%zu,\"filled_dark\":%zu,\"missing\":%zu},"
          "\"output\":{\"fitted_samples\":%zu,\"filled_samples\":%zu,"
          "\"fitted_dark\":%zu,\"filled_dark\":%zu,\"missing\":%zu}}",
          s->input_bake_fitted_samples,s->input_bake_filled_samples,
          s->input_bake_fitted_dark,s->input_bake_filled_dark,s->input_bake_missing,
          s->output_bake_fitted_samples,s->output_bake_filled_samples,
          s->output_bake_fitted_dark,s->output_bake_filled_dark,s->output_bake_missing);
        fprintf(f,
          ",\"multigrid\":{\"quilt\":{\"levels\":%d,\"corrections\":%d,"
          "\"residual_before\":%.9g,\"residual_after\":%.9g},"
          "\"position\":{\"levels\":%d,\"corrections\":%d,"
          "\"residual_before\":%.9g,\"residual_after\":%.9g,"
          "\"backtracks\":%d,\"rejects\":%d,"
          "\"patches_accepted\":%d,\"patches_rejected\":%d}}",
          s->quilt_mg_levels,s->quilt_mg_corrections,
          s->quilt_mg_residual_before,s->quilt_mg_residual_after,
          s->position_mg_levels,s->position_mg_corrections,
          s->position_mg_residual_before,s->position_mg_residual_after,
          s->position_mg_backtracks,s->position_mg_rejects,
          s->position_mg_patches_accepted,s->position_mg_patches_rejected);
        fputs(",\"quilt_regions\":{\"raw\":",f);
        snap_report_region_edge_stats(f,&s->quilt_raw_regions);
        fputs(",\"solved\":",f);
        snap_report_region_edge_stats(f,&s->quilt_solved_regions);
        fputs(",\"applied\":",f);
        snap_report_region_edge_stats(f,&s->quilt_applied_regions);
        fputc('}',f);
    }
    fputc('}',f);fflush(f);
}

static void snap_post_arap_report_entry(FILE *f,int *first,size_t component,
                                        long group,long track,
                                        const PostSnapArapTxn *s) {
    if(!f||!s)return;
    if(!*first)fputs(",\n",f);*first=0;
    fprintf(f,
      "  {\"component\":%zu,\"group\":%ld,\"representative_track\":%ld,"
      "\"stage\":\"post_snap_metric_arap\",\"accepted\":%s,"
      "\"solver_code\":%d,\"iterations\":%d,\"converged\":%s,"
      "\"stop_reason\":%d,\"effective_movement_tolerance\":%.17g,"
      "\"barrier_stalled\":%s,\"barrier_rejections\":%d,"
      "\"barrier_zero_fallbacks\":%d,\"minimum_step_scale\":%.17g,"
      "\"selected_transaction_scale\":%.17g,\"trial_exact_audits\":%d,"
      "\"orientation_rejects\":%d,\"locally_rolled_back_vertices\":%zu,"
      "\"movement_rms\":%.17g,\"movement_max\":%.17g,"
      "\"axial_drift_max\":%.17g,"
      "\"edge_log_rms\":{\"input\":%.17g,\"full_solve\":%.17g,"
      "\"selected\":%.17g},"
      "\"intersection_audit\":{"
      "\"input\":{\"conflicts\":%zu,\"overlap\":%zu,\"stab\":%zu,\"fold\":%zu},"
      "\"output\":{\"conflicts\":%zu,\"overlap\":%zu,\"stab\":%zu,\"fold\":%zu}}}",
      component,group,track,s->accepted?"true":"false",s->solver_code,
      s->solver.iterations,s->solver.converged?"true":"false",
      s->solver.stop_reason,s->solver.effective_movement_tolerance,
      s->solver.barrier_stalled?"true":"false",s->solver.barrier_rejections,
      s->solver.barrier_zero_fallbacks,s->solver.minimum_step_scale,
      s->selected_scale,s->trial_audits,s->orientation_rejects,
      s->locally_rolled_back_vertices,
      s->displacement_rms,s->displacement_max,s->axial_drift_max,
      s->input_edge_log_rms,s->solver.final_edge_log_rms,
      s->selected_edge_log_rms,
      s->input_intersections,s->input_overlap,s->input_stab,s->input_fold,
      s->output_intersections,s->output_overlap,s->output_stab,s->output_fold);
    fflush(f);
}

int main(int argc, char **argv) {
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        int f = 0;
        f += Json_read_selftest();
        f += QuadStrip_selftest();
        f += QuadStripSnap_selftest();
        f += MeshBin_selftest();
        f += solid_overlap_selftest();
        f += snap_transaction_selftest();
        f += initializer_quarantine_selftest();
        f += phase_gauge_selftest();
        f += loaded_crop_selftest();
        f += ribbon_coarsen_selftest();
        f += ribbon_masked_coarsen_selftest();
        f += unwrap_support_authority_selftest();
        f += cloth_support_selftest();
        f += cloth_band_frontier_selftest();
        f += cloth_smooth_list_selftest();
        f += cloth_settle_selftest();
        f += cloth_stall_escalation_selftest();
        f += continuation_trim_selftest();
        f += passthrough_faces_selftest();
        f += orientation_baseline_selftest();
        f += IntersectionCleanup_selftest();
        fprintf(stderr, "[selftest] solid_quad_ribbon: %s (%d failures)\n",
                f == 0 ? "ALL PASS" : "FAILURES", f);
        return f == 0 ? 0 : 1;
    }
    if (argc < 4) { usage(argv[0]); return 2; }

    const char *aggregate = argv[1];
    const char *tracks_root = argv[2];
    const char *out_obj = argv[3];
    const char *fallback = NULL, *reg_path = NULL, *topo_png = NULL, *raw_cubes = NULL;
    const char *correspondence_path = NULL;
    const char *bake_refit_masks = NULL;
    const char *snap_report_path = NULL;
    const char *arap_report_path = NULL;
    long only_group = LONG_MIN;
    int u_column_range = 0, u_column_first = 0, u_column_end = 0;
    int u_column_halo = 0;
    int ribbon_u_stride = 1, ribbon_v_stride = 1;
    int giant_ribbon = 0, metric_arap = 0, post_snap_arap = 0;
    int v92_snap_policy = 0, legacy_v92_registration = 0;
    int repair_folds = 0;
    int trust_unwrap_domain = 0;
    int require_collision_free = 0, collision_patience = 12, collision_rounds = 64;
    double collision_collar = 4.0, collision_displacement_bound = 128.0;
    int collision_settle_rounds = 16;
    double collision_settle_beta = 0.70;
    int continuation_trim = 0, trim_min_region = 24, trim_rim_guard = 2;
    int registration_max_gap = 8;
    const char *snips_dir = NULL;
    int snip_cols = 512, snips_only = 0, report_only = 0;
    int snap_side_set = 0, snap_axis_set = 0, snap_passes = 2;
    int remesh_sweeps = 1;
    double grid_du = 2.0, refine_anchor = 4.0, pack_width = 0.0;
    double seam_refine_multiplier = 1.0;
    double init_quarantine_sander=0.0,init_quarantine_fraction=0.20;
    int init_quarantine_tile=16,init_quarantine_dilate=3;
    QuadStripOpts qopts; QuadStrip_defaults(&qopts);
    QuadStripArapOpts arap_opts; QuadStripArap_defaults(&arap_opts);
    QuadStripSnapOpts sopts; QuadStripSnap_defaults(&sopts);
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--fallback-tracks-root") && i + 1 < argc) fallback = argv[++i];
        else if (!strcmp(argv[i], "--registration") && i + 1 < argc) reg_path = argv[++i];
        else if (!strcmp(argv[i], "--legacy-v92-registration")) legacy_v92_registration = 1;
        else if (!strcmp(argv[i], "--correspondence") && i + 1 < argc)
            correspondence_path = argv[++i];
        else if (!strcmp(argv[i], "--registration-max-gap") && i + 1 < argc)
            registration_max_gap=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bake-refit-masks") && i + 1 < argc) bake_refit_masks = argv[++i];
        else if (!strcmp(argv[i], "--only-group") && i + 1 < argc) only_group = atol(argv[++i]);
        else if (!strcmp(argv[i], "--u-column-range") && i + 2 < argc) {
            long first = strtol(argv[++i], NULL, 10);
            long end = strtol(argv[++i], NULL, 10);
            if (first < INT_MIN || first > INT_MAX || end < INT_MIN || end > INT_MAX) {
                fprintf(stderr, "--u-column-range values exceed integer range\n");
                return 2;
            }
            u_column_first = (int)first;
            u_column_end = (int)end;
            u_column_range = 1;
        }
        else if (!strcmp(argv[i], "--u-column-halo") && i + 1 < argc)
            u_column_halo = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--giant-ribbon")) giant_ribbon = 1;
        else if (!strcmp(argv[i], "--snips-dir") && i + 1 < argc) snips_dir = argv[++i];
        else if (!strcmp(argv[i], "--snip-cols") && i + 1 < argc) snip_cols = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snips-only")) snips_only = 1;
        else if (!strcmp(argv[i], "--report-only")) report_only = 1;
        else if (!strcmp(argv[i], "--pack-components") && i + 1 < argc) pack_width = atof(argv[++i]);
        else if (!strcmp(argv[i], "--grid-du") && i + 1 < argc) grid_du = atof(argv[++i]);
        else if (!strcmp(argv[i], "--ribbon-u-stride") && i + 1 < argc)
            ribbon_u_stride = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--ribbon-v-stride") && i + 1 < argc)
            ribbon_v_stride = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--topology-map") && i + 1 < argc) topo_png = argv[++i];
        else if (!strcmp(argv[i], "--fit-stiffness") && i + 1 < argc) qopts.fit_stiffness = atof(argv[++i]);
        else if (!strcmp(argv[i], "--cyl-axis") && i + 2 < argc) {
            qopts.cylindrical_fill = 1;
            qopts.axis_y = atof(argv[++i]);
            qopts.axis_x = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--snap-axis") && i + 2 < argc) {
            sopts.axis_y = atof(argv[++i]);
            sopts.axis_x = atof(argv[++i]);
            snap_axis_set = 1;
        }
        else if (!strcmp(argv[i], "--relax-rounds") && i + 1 < argc) qopts.relax_rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-only")) {
            metric_arap = 1; trust_unwrap_domain = 1;
        }
        else if (!strcmp(argv[i], "--trust-unwrap-domain")) trust_unwrap_domain = 1;
        else if (!strcmp(argv[i], "--arap-report") && i + 1 < argc) arap_report_path = argv[++i];
        else if (!strcmp(argv[i], "--arap-max-iterations") && i + 1 < argc) arap_opts.max_iterations=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-tolerance") && i + 1 < argc) arap_opts.movement_tolerance=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-tolerance-rel") && i + 1 < argc) arap_opts.movement_tolerance_rel=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-stall-window") && i + 1 < argc) arap_opts.movement_stall_window=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-stall-fraction") && i + 1 < argc) arap_opts.movement_stall_fraction=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-cg-iterations") && i + 1 < argc) arap_opts.cg_max_iterations=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-cg-tolerance") && i + 1 < argc) arap_opts.cg_relative_tolerance=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-mg-levels") && i + 1 < argc) arap_opts.multigrid_levels=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-mg-cycles") && i + 1 < argc) arap_opts.multigrid_cycles=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-mg-pre-sweeps") && i + 1 < argc) arap_opts.multigrid_pre_sweeps=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-mg-post-sweeps") && i + 1 < argc) arap_opts.multigrid_post_sweeps=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-mg-coarse-sweeps") && i + 1 < argc) arap_opts.multigrid_coarse_sweeps=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-nlmg-levels") && i + 1 < argc) arap_opts.nonlinear_multigrid_levels=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-nlmg-iterations") && i + 1 < argc) arap_opts.nonlinear_coarse_iterations=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-nlmg-min-vertices") && i + 1 < argc) arap_opts.nonlinear_min_vertices=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arap-source-weight") && i + 1 < argc) arap_opts.source_weight=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-fill-weight") && i + 1 < argc) arap_opts.fill_weight=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-huber") && i + 1 < argc) arap_opts.huber_delta=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-min-source-scale") && i + 1 < argc) arap_opts.minimum_source_scale=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-min-area-ratio") && i + 1 < argc) arap_opts.minimum_area_ratio=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-max-step") && i + 1 < argc) arap_opts.maximum_vertex_step=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-axis-normal-weight") && i + 1 < argc) arap_opts.axis_normal_weight=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-axis-normal-radius") && i + 1 < argc) arap_opts.axis_normal_radius=atof(argv[++i]);
        else if (!strcmp(argv[i], "--arap-axis-fold-limit") && i + 1 < argc) arap_opts.axis_fold_limit=atof(argv[++i]);
        else if (!strcmp(argv[i], "--post-snap-arap")) post_snap_arap = 1;
        else if (!strcmp(argv[i], "--init-quarantine-sander") && i + 1 < argc) init_quarantine_sander=atof(argv[++i]);
        else if (!strcmp(argv[i], "--init-quarantine-tile") && i + 1 < argc) init_quarantine_tile=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--init-quarantine-fraction") && i + 1 < argc) init_quarantine_fraction=atof(argv[++i]);
        else if (!strcmp(argv[i], "--init-quarantine-dilate") && i + 1 < argc) init_quarantine_dilate=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--require-collision-free")) require_collision_free = 1;
        else if (!strcmp(argv[i], "--repair-folds")) repair_folds = 1;
        else if (!strcmp(argv[i], "--collision-patience") && i + 1 < argc)
            collision_patience = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--collision-rounds") && i + 1 < argc)
            collision_rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--collision-collar") && i + 1 < argc)
            collision_collar = atof(argv[++i]);
        else if (!strcmp(argv[i], "--collision-displacement-bound") && i + 1 < argc)
            collision_displacement_bound = atof(argv[++i]);
        else if (!strcmp(argv[i], "--collision-settle-rounds") && i + 1 < argc)
            collision_settle_rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--collision-settle-beta") && i + 1 < argc)
            collision_settle_beta = atof(argv[++i]);
        else if (!strcmp(argv[i], "--continuation-trim"))
            continuation_trim = 1;
        else if (!strcmp(argv[i], "--trim-min-region") && i + 1 < argc)
            trim_min_region = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--trim-rim-guard") && i + 1 < argc)
            trim_rim_guard = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--raw-cubes") && i + 1 < argc) raw_cubes = argv[++i];
        else if (!strcmp(argv[i], "--require-complete-raw")) sopts.require_complete_raw = 1;
        else if (!strcmp(argv[i], "--allow-incomplete-raw")) sopts.require_complete_raw = 0;
        else if (!strcmp(argv[i], "--snap-reach") && i + 1 < argc) sopts.reach = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-bright-min") && i + 1 < argc) sopts.bright_min = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-band-frac") && i + 1 < argc) sopts.band_frac = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-gain") && i + 1 < argc) sopts.min_gain = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-ridge-lock") && i + 1 < argc) sopts.ridge_lock = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-dark-probe-frac") && i + 1 < argc) sopts.dark_probe_frac = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-guided-reach") && i + 1 < argc) sopts.guided_reach = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-winding-tube") && i + 1 < argc) sopts.winding_tube = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-recover")) sopts.patch_recover = 1;
        else if (!strcmp(argv[i], "--snap-no-patch-recover")) sopts.patch_recover = 0;
        else if (!strcmp(argv[i], "--snap-patch-buffer") && i + 1 < argc) sopts.patch_buffer = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-min-cells") && i + 1 < argc) sopts.patch_min_cells = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-ray-reach") && i + 1 < argc) sopts.patch_ray_reach = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-ray-margin") && i + 1 < argc) sopts.patch_ray_margin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-relax-sweeps") && i + 1 < argc) sopts.patch_relax_sweeps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-target-weight") && i + 1 < argc) sopts.patch_target_weight = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-patch-fitted-dark-slack") && i + 1 < argc) sopts.patch_fitted_dark_slack = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-no-local-contrast")) sopts.local_contrast = 0;
        else if (!strcmp(argv[i], "--snap-bake-window") && i + 2 < argc) {
            sopts.bake_window_low=atof(argv[++i]);sopts.bake_window_high=atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--snap-bake-dark") && i + 1 < argc) sopts.bake_dark_u8=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-side") && i + 1 < argc) {
            const char *side = argv[++i];
            if (!strcmp(side, "midline")) sopts.snap_side = QUAD_SNAP_SIDE_MIDLINE;
            else if (!strcmp(side, "recto")) sopts.snap_side = QUAD_SNAP_SIDE_RECTO;
            else if (!strcmp(side, "nearest")) sopts.snap_side = QUAD_SNAP_SIDE_NEAREST;
            else { fprintf(stderr, "--snap-side must be midline, recto, or nearest\n"); return 2; }
            snap_side_set = 1;
        }
        else if (!strcmp(argv[i], "--snap-quilt-smooth") && i + 1 < argc) sopts.quilt_smooth = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-quilt-seam-smooth") && i + 1 < argc) sopts.quilt_seam_smooth = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-quilt-seam-refine") && i + 1 < argc) seam_refine_multiplier = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-quilt-anchor") && i + 1 < argc) sopts.quilt_anchor = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-quilt-huber") && i + 1 < argc) sopts.quilt_huber = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-quilt-sweeps") && i + 1 < argc) sopts.quilt_sweeps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-rounds") && i + 1 < argc) sopts.rounds = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-passes") && i + 1 < argc) snap_passes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-policy") && i + 1 < argc) {
            const char *policy = argv[++i];
            if (!strcmp(policy, "modern")) v92_snap_policy = 0;
            else if (!strcmp(policy, "v92")) v92_snap_policy = 1;
            else {
                fprintf(stderr, "--snap-policy must be modern or v92\n");
                return 2;
            }
        }
        else if (!strcmp(argv[i], "--snap-report") && i + 1 < argc) snap_report_path = argv[++i];
        else if (!strcmp(argv[i], "--snap-refine-anchor") && i + 1 < argc) refine_anchor = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-fitted-source-anchor") && i + 1 < argc)
            sopts.fitted_source_anchor = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-sweeps") && i + 1 < argc) sopts.solve_sweeps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-omega") && i + 1 < argc) sopts.solve_omega = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-mg-levels") && i + 1 < argc) sopts.multigrid_levels = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-mg-cycles") && i + 1 < argc) sopts.multigrid_cycles = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-mg-patch") && i + 1 < argc) sopts.multigrid_patch = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--remesh-motion") && i + 1 < argc) sopts.remesh_motion = atof(argv[++i]);
        else if (!strcmp(argv[i], "--remesh-motion-max") && i + 1 < argc) sopts.remesh_motion_max = atof(argv[++i]);
        else if (!strcmp(argv[i], "--remesh-sander") && i + 1 < argc) sopts.remesh_sander = atof(argv[++i]);
        else if (!strcmp(argv[i], "--remesh-over4") && i + 1 < argc) sopts.remesh_over4 = 0.01 * atof(argv[++i]);
        else if (!strcmp(argv[i], "--remesh-sweeps") && i + 1 < argc) remesh_sweeps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--snap-w-target") && i + 1 < argc) sopts.w_target = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-w-orig") && i + 1 < argc) sopts.w_orig = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-w-smooth") && i + 1 < argc) sopts.w_smooth = atof(argv[++i]);
        else if (!strcmp(argv[i], "--snap-w-dev") && i + 1 < argc) sopts.w_dev = atof(argv[++i]);
        else if (!strcmp(argv[i], "--max-hole-distance") && i + 1 < argc)
            qopts.max_hole_distance=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--max-row-gap") && i + 1 < argc)
            qopts.max_row_gap=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--solid-mode") && i + 1 < argc) {
            const char *mode=argv[++i];
            if(!strcmp(mode,"bounded"))qopts.mode=QUAD_STRIP_BOUNDED;
            else if(!strcmp(mode,"rect"))qopts.mode=QUAD_STRIP_RECT;
            else if(!strcmp(mode,"span"))qopts.mode=QUAD_STRIP_SPAN;
            else {fprintf(stderr,"--solid-mode must be bounded, rect, or span\n");return 2;}
        } else { fprintf(stderr, "unknown arg: %s\n", argv[i]); usage(argv[0]); return 2; }
    }
    /* This is a named, immutable compatibility contract rather than a loose
     * bag of caller-selected switches.  Apply it after option parsing so a
     * reordered command line cannot silently produce a v92-labelled hybrid. */
    if (v92_snap_policy) {
        qopts.mode = QUAD_STRIP_RECT;
        qopts.relax_rounds = 0;
        trust_unwrap_domain = 0;
        metric_arap = 0;
        post_snap_arap = 0;
        repair_folds = 0;
        continuation_trim = 0;
        collision_settle_rounds = 0;
        sopts.require_complete_raw = 1;
        sopts.snap_side = QUAD_SNAP_SIDE_RECTO;
        snap_side_set = 1;
        sopts.patch_recover = 0;
        sopts.local_contrast = 1;
        sopts.guided_reach = 12.0;
        sopts.winding_tube = 3.0;
        sopts.bake_window_low = 47.0;
        sopts.bake_window_high = 193.0;
        sopts.bake_dark_u8 = 13;
        sopts.quilt_smooth = 4.0;
        sopts.quilt_seam_smooth = 1.0;
        sopts.quilt_huber = 1.5;
        sopts.quilt_sweeps = 4;
        sopts.fitted_source_anchor = 0.0;
        sopts.multigrid_levels = 8;
        sopts.multigrid_cycles = 1;
        sopts.multigrid_patch = 16;
        seam_refine_multiplier = 2.0;
        snap_passes = 4;
        sopts.rounds = 8;
        remesh_sweeps = 1;
    }
    if (giant_ribbon && only_group != LONG_MIN) {
        fprintf(stderr, "--giant-ribbon and --only-group are mutually exclusive\n");
        return 2;
    }
    if (u_column_range &&
        (long long)u_column_end - (long long)u_column_first < 2LL) {
        fprintf(stderr, "--u-column-range must contain at least two columns\n");
        return 2;
    }
    if (u_column_halo < 0 || (u_column_halo > 0 && !u_column_range)) {
        fprintf(stderr, "--u-column-halo must be non-negative and requires --u-column-range\n");
        return 2;
    }
    if (collision_patience < 0 || collision_patience > 256) {
        fprintf(stderr, "--collision-patience must be in [0,256]\n");
        return 2;
    }
    if (collision_collar < 1.0 || collision_collar > 32.0) {
        fprintf(stderr, "--collision-collar must be in [1,32]\n");
        return 2;
    }
    if (collision_displacement_bound < 8.0 ||
        collision_displacement_bound > 1024.0) {
        fprintf(stderr, "--collision-displacement-bound must be in [8,1024]\n");
        return 2;
    }
    if (collision_settle_rounds < 0 || collision_settle_rounds > 64) {
        fprintf(stderr, "--collision-settle-rounds must be in [0,64]\n");
        return 2;
    }
    if (trim_min_region < 1 || trim_min_region > 100000 ||
        trim_rim_guard < 0 || trim_rim_guard > 32) {
        fprintf(stderr, "--trim-min-region must be in [1,100000] and "
                        "--trim-rim-guard in [0,32]\n");
        return 2;
    }
    if (!(collision_settle_beta > 0.0) || !(collision_settle_beta < 1.0)) {
        fprintf(stderr, "--collision-settle-beta must be in (0,1)\n");
        return 2;
    }
    if (collision_rounds < 1 || collision_rounds > 256) {
        fprintf(stderr, "--collision-rounds must be in [1,256]\n");
        return 2;
    }
    if (ribbon_u_stride < 1 || ribbon_u_stride > 64 ||
        ribbon_v_stride < 1 || ribbon_v_stride > 64) {
        fprintf(stderr, "--ribbon-u-stride and --ribbon-v-stride must be in [1,64]\n");
        return 2;
    }
    if (arap_report_path && !metric_arap) {
        fprintf(stderr, "--arap-report requires --arap-only\n");
        return 2;
    }
    if (metric_arap && raw_cubes) {
        fprintf(stderr, "--arap-only stops after metric ARAP and cannot be combined with --raw-cubes\n");
        return 2;
    }
    if (post_snap_arap && (!raw_cubes || metric_arap)) {
        fprintf(stderr, "--post-snap-arap requires --raw-cubes and cannot be combined with --arap-only\n");
        return 2;
    }
    if (v92_snap_policy && (!raw_cubes || !qopts.cylindrical_fill)) {
        fprintf(stderr, "--snap-policy v92 requires --raw-cubes and --cyl-axis\n");
        return 2;
    }
    if (v92_snap_policy && bake_refit_masks) {
        fprintf(stderr, "--snap-policy v92 is the frozen no-bake-refit contract\n");
        return 2;
    }
    if (legacy_v92_registration && !v92_snap_policy) {
        fprintf(stderr, "--legacy-v92-registration requires --snap-policy v92\n");
        return 2;
    }
    if (post_snap_arap && qopts.mode != QUAD_STRIP_RECT) {
        fprintf(stderr, "--post-snap-arap requires complete rect topology\n");
        return 2;
    }
    if (metric_arap && (qopts.mode != QUAD_STRIP_RECT || qopts.relax_rounds != 0 ||
                        init_quarantine_sander > 0.0 || bake_refit_masks)) {
        fprintf(stderr, "--arap-only requires rect topology, --relax-rounds 0, no initializer quarantine, and no bake-refit masks\n");
        return 2;
    }
    if ((metric_arap || post_snap_arap) &&
                        (arap_opts.max_iterations < 1 || arap_opts.max_iterations > 1000 ||
                        !isfinite(arap_opts.movement_tolerance) || arap_opts.movement_tolerance <= 0.0 ||
                        !isfinite(arap_opts.movement_tolerance_rel) || arap_opts.movement_tolerance_rel < 0.0 ||
                        arap_opts.movement_tolerance_rel > 1.0 ||
                        arap_opts.movement_stall_window < 0 || arap_opts.movement_stall_window > 1000 ||
                        !isfinite(arap_opts.movement_stall_fraction) ||
                        arap_opts.movement_stall_fraction <= 0.0 || arap_opts.movement_stall_fraction >= 1.0 ||
                        arap_opts.cg_max_iterations < 1 || arap_opts.cg_max_iterations > 12000 ||
                        !isfinite(arap_opts.cg_relative_tolerance) || arap_opts.cg_relative_tolerance <= 0.0 ||
                        arap_opts.multigrid_levels < 1 || arap_opts.multigrid_levels > 16 ||
                        arap_opts.multigrid_cycles < 0 || arap_opts.multigrid_cycles > 8 ||
                        arap_opts.multigrid_pre_sweeps < 0 || arap_opts.multigrid_pre_sweeps > 16 ||
                        arap_opts.multigrid_post_sweeps < 0 || arap_opts.multigrid_post_sweeps > 16 ||
                        arap_opts.multigrid_coarse_sweeps < 1 || arap_opts.multigrid_coarse_sweeps > 1024 ||
                        arap_opts.nonlinear_multigrid_levels < 1 || arap_opts.nonlinear_multigrid_levels > 8 ||
                        arap_opts.nonlinear_coarse_iterations < 1 || arap_opts.nonlinear_coarse_iterations > 1000 ||
                        arap_opts.nonlinear_min_vertices < 0 ||
                        !isfinite(arap_opts.source_weight) || arap_opts.source_weight < 0.0 ||
                        !isfinite(arap_opts.fill_weight) || arap_opts.fill_weight < 0.0 ||
                        !(arap_opts.source_weight > 0.0 || arap_opts.fill_weight > 0.0) ||
                        !isfinite(arap_opts.huber_delta) || arap_opts.huber_delta <= 0.0 ||
                        !isfinite(arap_opts.minimum_source_scale) || arap_opts.minimum_source_scale < 0.0 ||
                        arap_opts.minimum_source_scale > 1.0 ||
                        !isfinite(arap_opts.minimum_area_ratio) || arap_opts.minimum_area_ratio <= 0.0 ||
                        arap_opts.minimum_area_ratio > 1.0 ||
                        !isfinite(arap_opts.maximum_vertex_step) || arap_opts.maximum_vertex_step <= 0.0 ||
                        !isfinite(arap_opts.axis_normal_weight) || arap_opts.axis_normal_weight < 0.0 || arap_opts.axis_normal_weight > 1.0 ||
                        !isfinite(arap_opts.axis_normal_radius) || arap_opts.axis_normal_radius < 0.0 ||
                         !isfinite(arap_opts.axis_fold_limit) || arap_opts.axis_fold_limit < 0.0 || arap_opts.axis_fold_limit >= 1.0)) {
        fprintf(stderr, "metric ARAP options are outside their valid range (outer cap is 1000)\n");
        return 2;
    }
    if ((metric_arap || post_snap_arap) && arap_opts.axis_normal_weight > 0.0 &&
        !qopts.cylindrical_fill && !snap_axis_set) {
        fprintf(stderr, "--arap-axis-normal-weight requires --cyl-axis or --snap-axis <y> <x>\n");
        return 2;
    }
    if (registration_max_gap < 0 || registration_max_gap > 1024) {
        fprintf(stderr, "--registration-max-gap must be in 0..1024 grid cells\n");
        return 2;
    }
    if (giant_ribbon && bake_refit_masks) {
        fprintf(stderr, "--bake-refit-masks is not yet defined for --giant-ribbon\n");
        return 2;
    }
    if (snips_only && !snips_dir) {
        fprintf(stderr, "--snips-only requires --snips-dir\n");
        return 2;
    }
    if (report_only && (!raw_cubes || !snap_report_path)) {
        fprintf(stderr, "--report-only requires --raw-cubes and --snap-report\n");
        return 2;
    }
    if (report_only && (snips_dir || topo_png)) {
        fprintf(stderr, "--report-only cannot be combined with --snips-dir or --topology-map\n");
        return 2;
    }
    int omit_mesh = snips_only || report_only;
    if (correspondence_path && omit_mesh) {
        fprintf(stderr, "--correspondence requires the monolithic VMesh output\n");
        return 2;
    }
    if (legacy_v92_registration && correspondence_path) {
        fprintf(stderr, "legacy v92 manifests have no material lineage and cannot emit correspondence\n");
        return 2;
    }
    if (snap_report_path && !raw_cubes) {
        fprintf(stderr, "--snap-report requires --raw-cubes\n");
        return 2;
    }
    if (raw_cubes && !sopts.require_complete_raw)
        fprintf(stderr, "WARNING: incomplete RAW is allowed: this run is diagnostic "
                        "and must not be published as a production sheet\n");
    if (snap_passes < 1 || snap_passes > 64) {
        fprintf(stderr, "--snap-passes must be in 1..64\n");
        return 2;
    }
    if (remesh_sweeps < 0 || remesh_sweeps > 64) {
        fprintf(stderr, "--remesh-sweeps must be in 0..64\n");
        return 2;
    }
    if (sopts.multigrid_levels < 1 || sopts.multigrid_levels > 16 ||
        sopts.multigrid_cycles < 0 || sopts.multigrid_cycles > 16 ||
        sopts.multigrid_patch < 0 || sopts.multigrid_patch > 4096) {
        fprintf(stderr, "--snap-mg-levels must be in 1..16, --snap-mg-cycles in 0..16, and --snap-mg-patch in 0..4096\n");
        return 2;
    }
    if(!isfinite(sopts.guided_reach)||sopts.guided_reach<sopts.reach||
       sopts.guided_reach>64.0||!isfinite(sopts.winding_tube)||
       sopts.winding_tube<0.0||sopts.winding_tube>4.5){
        fprintf(stderr,"--snap-guided-reach must be in [snap-reach,64] and --snap-winding-tube in 0..4.5\n");
        return 2;
    }
    if(sopts.patch_buffer<0||sopts.patch_buffer>8||
       sopts.patch_min_cells<1||
       !isfinite(sopts.patch_ray_reach)||sopts.patch_ray_reach<2.0||
       sopts.patch_ray_reach>64.0||
       !isfinite(sopts.patch_ray_margin)||sopts.patch_ray_margin<0.0||
       sopts.patch_ray_margin>sopts.patch_ray_reach||
       sopts.patch_relax_sweeps<1||sopts.patch_relax_sweeps>512||
       !isfinite(sopts.patch_target_weight)||sopts.patch_target_weight<=0.0||
       sopts.patch_target_weight>1000.0||sopts.patch_fitted_dark_slack<0||
       sopts.patch_fitted_dark_slack>1000000) {
        fprintf(stderr,"invalid --snap-patch-* option (buffer 0..8, min-cells >=1, "
                       "reach 2..64, margin 0..reach, sweeps 1..512, weight 0..1000, "
                       "fitted-dark-slack 0..1000000)\n");
        return 2;
    }
    if(!isfinite(sopts.bake_window_low)||!isfinite(sopts.bake_window_high)||
       sopts.bake_window_low>=sopts.bake_window_high||sopts.bake_dark_u8<0||
       sopts.bake_dark_u8>255){
        fprintf(stderr,"--snap-bake-window requires lo < hi and --snap-bake-dark in 0..255\n");
        return 2;
    }
    if (!isfinite(refine_anchor) ||
        refine_anchor < 0.0 || refine_anchor > 1000.0) {
        fprintf(stderr, "--snap-refine-anchor must be in 0..1000\n");
        return 2;
    }
    if (!isfinite(sopts.fitted_source_anchor) ||
        sopts.fitted_source_anchor < 0.0 ||
        sopts.fitted_source_anchor > 1000.0) {
        fprintf(stderr, "--snap-fitted-source-anchor must be in 0..1000\n");
        return 2;
    }
    if (!isfinite(sopts.dark_probe_frac) ||
        sopts.dark_probe_frac < 0.0 || sopts.dark_probe_frac > 1.0) {
        fprintf(stderr, "--snap-dark-probe-frac must be in 0..1\n");
        return 2;
    }
    if (!isfinite(sopts.quilt_seam_smooth) ||
        sopts.quilt_seam_smooth < 0.0 || sopts.quilt_seam_smooth > 100.0) {
        fprintf(stderr, "--snap-quilt-seam-smooth must be in 0..100\n");
        return 2;
    }
    if (!isfinite(sopts.quilt_huber) ||
        sopts.quilt_huber < 0.0 || sopts.quilt_huber > 64.0) {
        fprintf(stderr, "--snap-quilt-huber must be in 0..64 voxels\n");
        return 2;
    }
    if (!isfinite(seam_refine_multiplier) ||
        seam_refine_multiplier < 1.0 || seam_refine_multiplier > 100.0) {
        fprintf(stderr, "--snap-quilt-seam-refine must be in 1..100\n");
        return 2;
    }
    if (!isfinite(pack_width) || pack_width < 0.0) {
        fprintf(stderr, "--pack-components must be >= 0\n");
        return 2;
    }
    if(!isfinite(init_quarantine_sander)||init_quarantine_sander<0.0||
       init_quarantine_sander>64.0||init_quarantine_tile<1||
       init_quarantine_tile>1024||!isfinite(init_quarantine_fraction)||
       init_quarantine_fraction<=0.0||init_quarantine_fraction>1.0||
       init_quarantine_dilate<0||init_quarantine_dilate>32) {
        fprintf(stderr,"initializer quarantine options are outside their valid range\n");
        return 2;
    }
    if (!isfinite(sopts.remesh_motion) || sopts.remesh_motion < 0.0 ||
        !isfinite(sopts.remesh_motion_max) || sopts.remesh_motion_max < 0.0 ||
        !isfinite(sopts.remesh_sander) || sopts.remesh_sander < 0.0 ||
        !isfinite(sopts.remesh_over4) ||
        sopts.remesh_over4 < 0.0 || sopts.remesh_over4 > 1.0) {
        fprintf(stderr, "remesh thresholds must be finite/nonnegative and --remesh-over4 must be 0..100 percent\n");
        return 2;
    }
    if (qopts.max_hole_distance < 0 || qopts.max_hole_distance > 1024) {
        fprintf(stderr, "--max-hole-distance must be in 0..1024\n");
        return 2;
    }
    if (qopts.max_row_gap < 1 || qopts.max_row_gap > 64) {
        fprintf(stderr, "--max-row-gap must be in 1..64\n");
        return 2;
    }
    if ((qopts.cylindrical_fill&&(!isfinite(qopts.axis_y)||!isfinite(qopts.axis_x)))||
        (snap_axis_set&&(!isfinite(sopts.axis_y)||!isfinite(sopts.axis_x)))) {
        fprintf(stderr, "scroll axis coordinates must be finite\n");
        return 2;
    }
    if (snips_dir && snip_cols < 1) {
        fprintf(stderr, "snip output requires --snip-cols >= 1\n");
        return 2;
    }
    if (qopts.cylindrical_fill) {
        if (!snap_side_set) sopts.snap_side = QUAD_SNAP_SIDE_MIDLINE;
        sopts.axis_y = qopts.axis_y; sopts.axis_x = qopts.axis_x;
        sopts.preserve_axial = 1;
        arap_opts.axis_y = qopts.axis_y; arap_opts.axis_x = qopts.axis_x;
        arap_opts.preserve_axial = 1;
    } else if (snap_axis_set) {
        if (!snap_side_set) sopts.snap_side = QUAD_SNAP_SIDE_MIDLINE;
        sopts.preserve_axial = 1;
        arap_opts.axis_y = sopts.axis_y; arap_opts.axis_x = sopts.axis_x;
        arap_opts.preserve_axial = 1;
    } else if (sopts.snap_side != QUAD_SNAP_SIDE_NEAREST||sopts.patch_recover) {
        fprintf(stderr, "oriented snap/patch recovery requires --cyl-axis or --snap-axis <y> <x>\n");
        return 2;
    }
    sopts.verbose = 1;
    arap_opts.verbose = 1;

    const char *roots[2]; int nroots = 0;
    roots[nroots++] = tracks_root;
    if (fallback) roots[nroots++] = fallback;

    Arena_T cfg = Arena_new();

    /* z_lo from report.json */
    char path[1400];
    snprintf(path, sizeof path, "%s/report.json", aggregate);
    const char *err = NULL;
    const JsonValue *report = Json_parse_file(cfg, path, &err);
    if (!report) { fprintf(stderr, "cannot read %s: %s\n", path, err ? err : "?"); return 1; }
    long z_lo = Json_as_long(Json_array_get(Json_object_get(report, "bbox_zyx"), 0), 0);

    /* manifest tracks */
    snprintf(path, sizeof path, "%s/material_tracks/manifest.json", aggregate);
    const JsonValue *manifest = Json_parse_file(cfg, path, &err);
    if (!manifest) { fprintf(stderr, "cannot read %s: %s\n", path, err ? err : "?"); return 1; }
    const JsonValue *tracks = Json_object_get(manifest, "tracks");
    size_t ntr = Json_array_len(tracks);
    long chunk = Json_member_long(manifest, "chunk_size", 128);

    /* optional registration: tid -> super_component */
    const JsonValue *charts = NULL;
    const JsonValue *registration_root = NULL;
    if (reg_path) {
        const JsonValue *reg = Json_parse_file(cfg, reg_path, &err);
        if (!reg) { fprintf(stderr, "cannot read %s: %s\n", reg_path, err ? err : "?"); return 1; }
        registration_root = reg;
        charts = Json_object_get(reg, "charts");
    }
    if (legacy_v92_registration) {
        const char *schema = registration_root
            ? Json_as_string(Json_object_get(registration_root, "schema")) : NULL;
        const JsonValue *edges = registration_root
            ? Json_object_get(registration_root, "edges") : NULL;
        int saw_2_9 = 0, saw_20_24 = 0;
        if (!schema || strcmp(schema, "chart-merge-v1") != 0 ||
            Json_array_len(edges) != 2) {
            fprintf(stderr,
                    "--legacy-v92-registration accepts only the audited "
                    "two-edge chart-merge-v1 file\n");
            return 1;
        }
        for (size_t e = 0; e < Json_array_len(edges); e++) {
            const JsonValue *edge = Json_array_get(edges, e);
            long a = Json_member_long(edge, "a", LONG_MIN);
            long b = Json_member_long(edge, "b", LONG_MIN);
            if (a > b) { long t = a; a = b; b = t; }
            if (a == 2 && b == 9 && !saw_2_9) saw_2_9 = 1;
            else if (a == 20 && b == 24 && !saw_20_24) saw_20_24 = 1;
            else {
                fprintf(stderr,
                        "--legacy-v92-registration found an unaudited edge %ld/%ld\n",
                        a, b);
                return 1;
            }
        }
        if (!saw_2_9 || !saw_20_24) {
            fprintf(stderr, "--legacy-v92-registration edge set is incomplete\n");
            return 1;
        }
    }

    /* Solve/emit rectangular fit charts independently.  A semantic
     * reconstruction component may contain many charts around local branches;
     * grouping those charts here would span-fill the gaps and recreate the
     * false horizontal weld that chart splitting was meant to remove. */
    long *tid_of = ARENA_ALLOC(cfg, ntr * sizeof *tid_of);
    long *grp_of = ARENA_ALLOC(cfg, ntr * sizeof *grp_of);
    long *material_of = ARENA_ALLOC(cfg, ntr * sizeof *material_of);
    unsigned char *material_explicit = ARENA_ALLOC(
        cfg, ntr * sizeof *material_explicit);
    long next_unique = 1;
    /* find max super to keep unique keys disjoint from super ids */
    long max_super = -1;
    for (size_t i = 0; i < ntr; i++) {
        const JsonValue *tr = Json_array_get(tracks, i);
        long tid = Json_member_long(tr, "global_track", -1);
        tid_of[i] = tid;
        const JsonValue *fit_component_value =
            Json_object_get(tr, "fit_component");
        long fit_component = Json_member_long(tr, "fit_component", tid);
        long reconstruction = Json_member_long(
            tr, "reconstruction_component", tid);
        const JsonValue *material_value = Json_object_get(
            tr, "material_identity");
        material_explicit[i] = (unsigned char)(material_value != NULL);
        material_of[i] = Json_member_long(
            tr, "material_identity",
            Json_member_long(tr, "relation_identity", reconstruction));
        if (legacy_v92_registration && material_explicit[i]) {
            fprintf(stderr,
                    "--legacy-v92-registration is only for pre-lineage manifests; "
                    "track %ld already has material_identity\n", tid);
            return 1;
        }
        if (correspondence_path &&
            (!material_explicit[i] || material_of[i] < INT32_MIN ||
             material_of[i] > INT32_MAX)) {
            fprintf(stderr,
                    "--correspondence requires an explicit int32 material_identity "
                    "for every track (track %ld)\n", tid);
            return 1;
        }
        long super = -1;
        if (charts) {
            char key[32];
            snprintf(key, sizeof key, "%ld", tid);
            const JsonValue *ce = Json_object_get(charts, key);
            if (ce) super = Json_member_long(ce, "super_component", -1);
        }
        grp_of[i] = super >= 0 ? super
                  : fit_component_value != NULL ? fit_component
                  : reconstruction;
        if (super > max_super) max_super = super;
    }
    for (size_t i = 0; i < ntr; i++)
        if (grp_of[i] > max_super) max_super = grp_of[i];
    for (size_t i = 0; i < ntr; i++) {
        if (grp_of[i] < 0) grp_of[i] = max_super + (next_unique++);
    }
    if (!giant_ribbon) {
        for (size_t i = 0; i < ntr; i++) for (size_t j = i + 1; j < ntr; j++) {
            if (grp_of[i] == grp_of[j] && !legacy_v92_registration &&
                (material_explicit[i] != material_explicit[j] ||
                 material_of[i] != material_of[j])) {
                fprintf(stderr,
                        "[registration] refusing group %ld: tracks %ld/%ld cross "
                        "material identities %ld/%ld\n",
                        grp_of[i], tid_of[i], tid_of[j],
                        material_of[i], material_of[j]);
                return 1;
            }
        }
    }
    if (!giant_ribbon && charts) {
        Arena_T split_arena = Arena_new();
        for (size_t i = 0; i < ntr; i++) for (size_t j = i + 1; j < ntr; j++) {
            if (grp_of[i] != grp_of[j]) continue;
            char da[1200], db[1200];
            if (resolve_track_dir(roots, nroots, tid_of[i], da, sizeof da) != 0
             || resolve_track_dir(roots, nroots, tid_of[j], db, sizeof db) != 0) continue;
            Arena_free(split_arena);
            Loaded a, b; a.tid = tid_of[i]; b.tid = tid_of[j];
            int lra=load_track(split_arena, da, Json_array_get(tracks, i), grid_du, z_lo,
                               trust_unwrap_domain,qopts.cylindrical_fill,&a);
            int lrb=load_track(split_arena, db, Json_array_get(tracks, j), grid_du, z_lo,
                               trust_unwrap_domain,qopts.cylindrical_fill,&b);
            if(lra==-3||lrb==-3)return 1;
            if(lra!=0||lrb!=0)continue;
            size_t n = 0, bad = 0; double mean = 0.0, mx = 0.0;
            int row_gap=0,col_gap=0;
            loaded_overlap_audit(&a, &b, 5.0, &n, &bad, &mean, &mx);
            int developed_gap=loaded_developed_gap(&a,&b,&row_gap,&col_gap);
            int bad_overlap=n>=25&&bad*10>=n*9;
            if (developed_gap>registration_max_gap || bad_overlap) {
                long old = grp_of[j];
                grp_of[j] = max_super + (next_unique++);
                fprintf(stderr, "[registration] split tracks %ld/%ld from group %ld: "
                                "developed row/col gap %d/%d cells (limit %d); "
                                "%zu overlapping certified cells, %.1f%% conflict, "
                                "mean/max 3-D gap %.2f/%.2f vox\n",
                        tid_of[i], tid_of[j], old,
                        row_gap,col_gap,registration_max_gap,
                        n,n?100.0*(double)bad/(double)n:0.0,mean,mx);
            }
        }
        Arena_dispose(&split_arena);
    }
    if (giant_ribbon) {
        for (size_t i = 0; i < ntr; i++) grp_of[i] = 0;
        fprintf(stderr, "[giant] assembling all fitted evidence on one global u/v lattice\n");
    }

    /* distinct group keys, in first-appearance order */
    long *groups = ARENA_ALLOC(cfg, ntr * sizeof *groups);
    size_t ngroups = 0;
    for (size_t i = 0; i < ntr; i++) {
        int seen = 0;
        for (size_t g = 0; g < ngroups; g++) if (groups[g] == grp_of[i]) { seen = 1; break; }
        if (!seen) groups[ngroups++] = grp_of[i];
    }

    /* output accumulators */
    FVec verts = {0}, uv = {0}, corr_uv = {0};
    IVec faces = {0}, corr_namespace = {0}, corr_lineage = {0};
    BVec filled = {0};
    IVec compv = {0};  /* per-vertex component id (for the PNG) */
    size_t used = 0, skipped = 0;
    size_t fitted_total = 0, filled_total = 0;
    size_t snap_target_total = 0, snap_nopap_total = 0;
    size_t trim_cells_total = 0, trim_faces_total = 0;
    size_t snap_pass_accepted_total = 0, snap_pass_rejected_total = 0;
    size_t arap_component_total = 0, arap_converged_total = 0;
    size_t arap_iteration_total = 0, arap_barrier_total = 0;
    size_t arap_barrier_zero_total = 0, arap_barrier_stall_total = 0;
    size_t arap_cap_total = 0;
    InitQuarantineStats init_quarantine_total = {0};
    size_t streamed_nv = 0, streamed_nf = 0, snip_serial = 0;
    PassthroughStats passthrough_stats = {0};
    FinalFoldRepairStats fold_repair_total = {0};
    PassthroughAtlas *passthrough_atlas = NULL;
    size_t passthrough_atlas_n = 0;
    double pack_u = 0.0, pack_v = 0.0, pack_row_h = 0.0;
    double *atlas_box_u0 = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_box_u0);
    double *atlas_box_u1 = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_box_u1);
    double *atlas_box_v0 = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_box_v0);
    double *atlas_box_v1 = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_box_v1);
    long *atlas_track = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_track);
    int *atlas_height = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_height);
    int *atlas_width = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_width);
    unsigned char *atlas_written = ARENA_ALLOC(cfg, ngroups * sizeof *atlas_written);
    memset(atlas_written, 0, ngroups * sizeof *atlas_written);
    FILE *snip_manifest = NULL;
    int snip_manifest_first = 1;
    FILE *snap_report = NULL;
    int snap_report_first = 1;
    FILE *arap_report = NULL;
    int arap_report_first = 1;
    if (snips_dir) {
        if (make_output_dir(snips_dir) != 0) {
            fprintf(stderr, "cannot create snips directory %s\n", snips_dir); return 1;
        }
        char mp[1600]; snprintf(mp, sizeof mp, "%s/manifest.json", snips_dir);
        snip_manifest = fopen(mp, "wb");
        if (!snip_manifest) { fprintf(stderr, "cannot write %s\n", mp); return 1; }
        fprintf(snip_manifest,
                "{\n \"schema\":\"quad-ribbon-snips-v2\",\n \"grid_du\":%.17g,\n"
                " \"ribbon_stride_uv\":[%d,%d],\n"
                " \"shared_boundary\":\"adjacent snips duplicate one vertex column; faces are unique; join by (group,global_row,global_col)\",\n"
                " \"snips\":[\n", grid_du,ribbon_u_stride,ribbon_v_stride);
    }
    if (arap_report_path) {
        arap_report=fopen(arap_report_path,"wb");
        if(!arap_report){fprintf(stderr,"cannot write %s\n",arap_report_path);return 1;}
        fprintf(arap_report,
          "{\n \"schema\":\"quad-ribbon-metric-arap-v1\",\n"
          " \"policy\":{\"geometry_terminal\":\"metric_arap\","
          "\"unwrap_domain\":\"finite-direct-support\","
          "\"rest_metric\":\"regular-uv-triangle-edge-lengths\","
          "\"ribbon_stride_uv\":[%d,%d],"
          "\"all_vertices_movable\":true,\"statistical_certification_masks_topology\":false,"
          "\"max_iterations\":%d,\"movement_tolerance\":%.17g,"
          "\"movement_tolerance_rel\":%.17g,"
          "\"movement_stall_window\":%d,\"movement_stall_fraction\":%.17g,"
          "\"cg_max_iterations\":%d,\"cg_relative_tolerance\":%.17g,"
          "\"multigrid_levels\":%d,\"multigrid_cycles\":%d,"
          "\"multigrid_pre_sweeps\":%d,\"multigrid_post_sweeps\":%d,"
          "\"multigrid_coarse_sweeps\":%d,"
          "\"nonlinear_multigrid_levels\":%d,"
          "\"nonlinear_coarse_iterations\":%d,"
          "\"nonlinear_min_vertices\":%d,"
          "\"inexact_global_forcing\":\"motion-scaled-then-final-cg-tolerance\","
          "\"source_weight\":%.17g,\"fill_weight\":%.17g,"
          "\"huber_delta\":%.17g,\"minimum_source_scale\":%.17g,"
          "\"minimum_area_ratio\":%.17g,\"maximum_vertex_step\":%.17g},\n \"components\":[\n",
          ribbon_u_stride,ribbon_v_stride,
          arap_opts.max_iterations,arap_opts.movement_tolerance,
          arap_opts.movement_tolerance_rel,
          arap_opts.movement_stall_window,arap_opts.movement_stall_fraction,
          arap_opts.cg_max_iterations,arap_opts.cg_relative_tolerance,
          arap_opts.multigrid_levels,arap_opts.multigrid_cycles,
          arap_opts.multigrid_pre_sweeps,arap_opts.multigrid_post_sweeps,
          arap_opts.multigrid_coarse_sweeps,
          arap_opts.nonlinear_multigrid_levels,
          arap_opts.nonlinear_coarse_iterations,
          arap_opts.nonlinear_min_vertices,
          arap_opts.source_weight,arap_opts.fill_weight,arap_opts.huber_delta,
          arap_opts.minimum_source_scale,arap_opts.minimum_area_ratio,
          arap_opts.maximum_vertex_step);
    }
    if (snap_report_path) {
        snap_report=fopen(snap_report_path,"wb");
        if(!snap_report){fprintf(stderr,"cannot write %s\n",snap_report_path);return 1;}
        fprintf(snap_report,
          "{\n \"schema\":\"quad-ribbon-snap-transactions-v13\",\n"
          " \"policy\":{\"bake_contract\":\"%s\","
          "\"initializer_domain\":\"%s\","
          "\"topology_domain\":\"%s\","
          "\"hard_trust_domain\":\"%s\","
          "\"snap_side\":\"%s\",\"refine_anchor\":%.9g,"
          "\"fitted_source_anchor\":%.9g,"
          "\"fitted_source_anchor_reference\":\"persistent-pre-ct-upstream-position\","
          "\"require_complete_raw\":%s,"
          "\"upstream_geometry_frozen\":%s,"
          "\"initializer_quarantine\":{\"sander_threshold\":%.9g,"
          "\"tile_cells\":%d,\"severe_fraction\":%.9g,\"dilate_tiles\":%d},"
          "\"dark_probe_frac\":%.9g,"
          "\"guided_reach\":%.9g,\"position_displacement_cap\":%.9g,"
          "\"winding_guard\":\"adaptive-local-45pct-adjacent-wrap\","
          "\"winding_min_half_width\":%.9g,"
          "\"preserve_axial\":%s,"
          "\"patch_recovery\":{\"enabled\":%s,\"buffer_rings\":%d,"
          "\"minimum_dark_cells\":%d,\"ray_reach\":%.9g,"
          "\"ray_margin\":%.9g,\"relax_sweeps\":%d,"
          "\"target_weight\":%.9g,"
          "\"agreement\":\"coherent-rim-majority-else-unanimous-or-wobbly\","
          "\"geometry_repair\":\"whole-patch-attenuation-to-one-eighth\","
          "\"fitted_dark_slack_samples\":%d,"
          "\"boundary\":\"trusted-xyz-region-zero-displacement-dirichlet\"},"
          "\"collision_output\":\"%s\","
          "\"collision_patience\":%d,\"collision_rounds\":%d,"
          "\"collision_collar\":%.9g,\"collision_displacement_bound\":%.9g,"
          "\"collision_settle_rounds\":%d,\"collision_settle_beta\":%.9g,"
          "\"continuation_trim\":{\"enabled\":%s,\"min_region_cells\":%d,"
          "\"rim_guard_cells\":%d,"
          "\"rule\":\"all-filled-corners-and-never-lit-and-no-supported-alternative\"},"
          "\"u_column_halo\":%d,"
          "\"ribbon_stride_uv\":[%d,%d],"
          "\"final_fold_repair\":{\"enabled\":%s,"
          "\"scope\":\"edge-adjacent-overlap-only\","
          "\"ordinary_layer_overlaps_preserved\":true,"
          "\"acceptance\":\"zero-residual-folds\"},"
          "\"exact_intersection_gate\":\"whole-component-no-growth\","
          "\"exact_audit_reuse\":\"snap-transaction-proof-reused-by-outer-pass-gate\","
          "\"post_snap_metric_arap\":{\"enabled\":%s,"
          "\"anchors\":\"all-current-snap-vertices\","
          "\"source_weight\":%.9g,\"huber_delta\":%.9g,"
          "\"minimum_source_scale\":%.9g,\"axis_normal_weight\":%.9g,"
          "\"preserve_axial\":%s,"
          "\"acceptance\":\"face-local-rollback-then-exact-no-growth-plus-metric-or-conflict-progress\"},"
          "\"local_contrast\":%s,\"bake_window\":[%.9g,%.9g],"
          "\"bake_dark_u8\":%d,"
          "\"quilt_basis\":\"world_zyx_vector\","
          "\"quilt_smooth\":%.9g,\"quilt_seam_smooth\":%.9g,"
          "\"quilt_huber\":%.9g,"
          "\"quilt_sweeps\":%d,\"quilt_conf_min\":%.9g,"
          "\"multigrid_levels\":%d,\"multigrid_cycles\":%d,"
          "\"multigrid_patch\":%d,"
          "\"quilt_seam_refine\":%.9g,\"quilt_seam_trigger_ratio\":%.9g,"
          "\"quilt_seam_trigger_step\":%.9g,\"quilt_seam_min_edges\":%zu,"
          "\"p95_step\":%.9g,\"p99_step\":%.9g,"
          "\"over4_factor\":%.9g,\"over4_step\":%.9g,"
          "\"p95_first_envelope\":%.9g,\"p99_first_envelope\":%.9g,"
          "\"over4_first_envelope\":%.9g,\"mature_from_pass\":3,"
          "\"mature_p95_epsilon\":%.9g,\"mature_p99_epsilon\":%.9g,"
          "\"mature_quilt_gain\":%.9g,\"mature_ct_quartile_sum_gain\":%d},"
          "\n \"passes\":[\n",
          v92_snap_policy?"v92":"modern",
          trust_unwrap_domain?"finite-nonzero-confidence-source-support":
                              "certified-or-source-only",
          qopts.mode==QUAD_STRIP_RECT?"component-rectangle":
            (qopts.mode==QUAD_STRIP_SPAN?"row-column-span-hull":
             "closed-single-island-bounded-holes"),
          sopts.fitted_source_anchor>0.0
            ? (trust_unwrap_domain
                ? "nonzero-confidence-upstream-support-soft-attached"
                : "authoritative-upstream-geometry-soft-attached")
            : (trust_unwrap_domain
                ? "nonzero-confidence-upstream-support-frozen"
                : "all-authoritative-upstream-geometry-frozen"),
          sopts.snap_side==QUAD_SNAP_SIDE_MIDLINE?"midline":
            (sopts.snap_side==QUAD_SNAP_SIDE_RECTO?"recto":"nearest"),
          refine_anchor,
          sopts.fitted_source_anchor,
          sopts.require_complete_raw?"true":"false",
          sopts.fitted_source_anchor>0.0?"false":"true",
          init_quarantine_sander,init_quarantine_tile,
          init_quarantine_fraction,init_quarantine_dilate,
          sopts.dark_probe_frac,sopts.guided_reach,sopts.guided_reach,
          sopts.winding_tube,sopts.preserve_axial?"true":"false",
          sopts.patch_recover?"true":"false",sopts.patch_buffer,
          sopts.patch_min_cells,sopts.patch_ray_reach,sopts.patch_ray_margin,
          sopts.patch_relax_sweeps,sopts.patch_target_weight,
          sopts.patch_fitted_dark_slack,
          v92_snap_policy?"disabled-v92-contract":
            (require_collision_free?"require-zero-long-and-stab":"best-exact-audited-partial"),
          collision_patience,collision_rounds,
          collision_collar,collision_displacement_bound,
          collision_settle_rounds,collision_settle_beta,
          continuation_trim?"true":"false",trim_min_region,trim_rim_guard,
          u_column_halo,
          ribbon_u_stride,ribbon_v_stride,repair_folds?"true":"false",
          post_snap_arap?"true":"false",arap_opts.source_weight,
          arap_opts.huber_delta,arap_opts.minimum_source_scale,
          arap_opts.axis_normal_weight,arap_opts.preserve_axial?"true":"false",
          sopts.local_contrast?"true":"false",sopts.bake_window_low,
          sopts.bake_window_high,sopts.bake_dark_u8,
          sopts.quilt_smooth,sopts.quilt_seam_smooth,sopts.quilt_huber,
          sopts.quilt_sweeps,sopts.quilt_conf_min,
          sopts.multigrid_levels,sopts.multigrid_cycles,sopts.multigrid_patch,
          seam_refine_multiplier,SNAP_QUILT_SEAM_TRIGGER_RATIO,
          SNAP_QUILT_SEAM_TRIGGER_STEP,SNAP_QUILT_SEAM_MIN_EDGES,
          SNAP_TXN_P95_STEP,SNAP_TXN_P99_STEP,SNAP_TXN_OVER4_FACTOR,
          SNAP_TXN_OVER4_STEP,SNAP_TXN_P95_ENVELOPE,SNAP_TXN_P99_ENVELOPE,
          SNAP_TXN_OVER4_ENVELOPE,SNAP_TXN_MATURE_P95_EPS,
          SNAP_TXN_MATURE_P99_EPS,SNAP_TXN_MATURE_QUILT_GAIN,
          SNAP_TXN_MATURE_CT_GAIN);
    }

    Arena_T scratch = Arena_new();
    Loaded *mem = ARENA_ALLOC(cfg, ntr * sizeof *mem);  /* per-group member buffer */

    for (size_t g = 0; g < ngroups; g++) {
        if (only_group != LONG_MIN && groups[g] != only_group) continue;
        Arena_free(scratch);
        size_t nm = 0;
        for (size_t i = 0; i < ntr; i++) {
            if (grp_of[i] != groups[g]) continue;
            char dir[1200];
            if (resolve_track_dir(roots, nroots, tid_of[i], dir, sizeof dir) != 0) { skipped++; continue; }
            const JsonValue *tr = Json_array_get(tracks, i);
            Loaded ld; ld.tid = tid_of[i];
            int load_rc=load_track(scratch, dir, tr, grid_du, z_lo,
                                   trust_unwrap_domain,qopts.cylindrical_fill,&ld);
            if(load_rc==-3)return 1;
            if(load_rc!=0){skipped++;continue;}
            if (u_column_range) {
                int original_c0 = ld.c0, original_w = ld.W;
                long long solve_first_ll=(long long)u_column_first-u_column_halo;
                long long solve_end_ll=(long long)u_column_end+u_column_halo;
                int solve_first=solve_first_ll<INT_MIN?INT_MIN:
                                solve_first_ll>INT_MAX?INT_MAX:(int)solve_first_ll;
                int solve_end=solve_end_ll<INT_MIN?INT_MIN:
                              solve_end_ll>INT_MAX?INT_MAX:(int)solve_end_ll;
                int crop_rc = crop_loaded_columns(
                    scratch, &ld, solve_first, solve_end);
                if (crop_rc > 0) continue;
                if (crop_rc < 0) {
                    fprintf(stderr, "UV-column crop failed for track %ld\n", ld.tid);
                    return 1;
                }
                fprintf(stderr,
                        "[uv-crop] track=%ld requested=[%d,%d) halo=%d "
                        "solve=[%d,%d) source=[%d,%d) retained=[%d,%d) "
                        "grid=%dx%d\n",
                        ld.tid, u_column_first, u_column_end,u_column_halo,
                        solve_first,solve_end,
                        original_c0, original_c0 + original_w,
                        ld.c0, ld.c0 + ld.W, ld.H, ld.W);
            }
            /* Diagnose each fitted initializer while its dense field is still
             * intact.  The merged group stores only authoritative observations;
             * evaluating that sparse buffer would either miss alternating-row
             * failures or inspect uninitialised off-domain memory. */
            if(init_quarantine_sander>0.0) {
                InitQuarantineStats iq={0};
                if(quarantine_initializer_geometry(
                        scratch,ld.domain,ld.trusted,ld.field,ld.H,ld.W,grid_du,
                        init_quarantine_sander,init_quarantine_tile,
                        init_quarantine_fraction,init_quarantine_dilate,&iq)!=0) {
                    fprintf(stderr,
                            "initializer quarantine failed for track %ld "
                            "(group %ld)\n",ld.tid,groups[g]);
                    return 1;
                }
                init_quarantine_total.measured_triangles+=iq.measured_triangles;
                init_quarantine_total.severe_triangles+=iq.severe_triangles;
                init_quarantine_total.marked_tiles+=iq.marked_tiles;
                init_quarantine_total.dilated_tiles+=iq.dilated_tiles;
                init_quarantine_total.removed_initializer+=iq.removed_initializer;
                init_quarantine_total.removed_trusted+=iq.removed_trusted;
            }
            mem[nm++] = ld;
        }
        if (nm == 0) continue;

        /* union bbox over the group's fragments */
        int gr0 = mem[0].gr0, c0 = mem[0].c0, gr1 = mem[0].gr0 + mem[0].H, c1 = mem[0].c0 + mem[0].W;
        for (size_t k = 1; k < nm; k++) {
            if (mem[k].gr0 < gr0) gr0 = mem[k].gr0;
            if (mem[k].c0 < c0) c0 = mem[k].c0;
            if (mem[k].gr0 + mem[k].H > gr1) gr1 = mem[k].gr0 + mem[k].H;
            if (mem[k].c0 + mem[k].W > c1) c1 = mem[k].c0 + mem[k].W;
        }
        int H = gr1 - gr0, W = c1 - c0;
        size_t HW = (size_t)H * (size_t)W;
        fprintf(stderr, "[component] index=%zu group=%ld representative_track=%ld "
                        "members=%zu grid=%dx%d origin=(v=%d,u_col=%d)\n",
                g, groups[g], mem[0].tid, nm, H, W, gr0, c0);
        uint8_t *gtopology = ARENA_ALLOC(scratch, HW);
        uint8_t *gdomain = ARENA_ALLOC(scratch, HW);
        uint8_t *gtrusted = ARENA_ALLOC(scratch, HW);
        uint8_t *gconflict = ARENA_CALLOC(scratch, HW, 1);
        float   *gfield  = ARENA_ALLOC(scratch, HW * 3 * sizeof *gfield);
        float   *gphase = qopts.cylindrical_fill
                        ? ARENA_ALLOC(scratch,HW*sizeof *gphase) : NULL;
        uint16_t *gcount = giant_ribbon ? ARENA_CALLOC(scratch, HW, sizeof *gcount) : NULL;
        size_t overlap_samples = 0, overlap_compatible = 0, overlap_conflict = 0;
        size_t overlap_conflict_cells = 0;
        double overlap_gap_sum = 0.0, overlap_gap_max = 0.0;
        memset(gtopology, 0, HW);
        memset(gdomain, 0, HW);
        memset(gtrusted, 0, HW);
        for (size_t k = 0; k < nm; k++) {
            const Loaded *m = &mem[k];
            int dr = m->gr0 - gr0, dc = m->c0 - c0;
            for (int j = 0; j < m->H; j++) for (int ii = 0; ii < m->W; ii++) {
                size_t sc = (size_t)j * (size_t)m->W + (size_t)ii;
                size_t gc = (size_t)(j + dr) * (size_t)W + (size_t)(ii + dc);
                /* Topology is a union of same-material footprints and survives
                 * initializer conflicts/quarantine.  Those operations turn a
                 * cell into a movable PDE unknown; they must not cut the sheet. */
                if (m->topology[sc]) gtopology[gc] = 1;
                if (!m->domain[sc]) continue;
                int incoming_trusted=m->trusted[sc]!=0;
                double incoming_phase=gphase?(double)m->phase[sc]:0.0;
                /* Once two physical geometries disagree at this UV cell, no
                 * later member may repopulate it.  It is a PDE unknown, not an
                 * ownership contest. */
                if (gconflict[gc]) continue;
                if (gdomain[gc]) {
                    double dz = (double)gfield[gc*3+0] - (double)m->field[sc*3+0];
                    double dy = (double)gfield[gc*3+1] - (double)m->field[sc*3+1];
                    double dx = (double)gfield[gc*3+2] - (double)m->field[sc*3+2];
                    double gap = sqrt(dz*dz + dy*dy + dx*dx);
                    overlap_samples++;
                    overlap_gap_sum += gap;
                    if (gap > overlap_gap_max) overlap_gap_max = gap;
                    if (gap <= 5.0) {
                        overlap_compatible++;
                        if(gphase) {
                            const double two_pi=6.2831853071795864769;
                            incoming_phase+=nearbyint(
                                ((double)gphase[gc]-incoming_phase)/two_pi)*two_pi;
                        }
                    } else {
                        /* Never choose, average, or overwrite distinct wraps.
                         * Erase the cell from the Dirichlet domain and let the
                         * cylindrical r/phi continuation reconstruct it solely
                         * from unambiguous neighbours. */
                        overlap_conflict++;
                        overlap_conflict_cells++;
                        gconflict[gc] = 1;
                        gdomain[gc] = 0;
                        gtrusted[gc] = 0;
                        if (gcount) gcount[gc] = 0;
                        continue;
                    }
                    if (giant_ribbon) {
                        if (gap <= 5.0) {
                            unsigned cnt = gcount[gc] ? gcount[gc] : 1u;
                            double den = (double)(cnt + 1u);
                            for (int ch = 0; ch < 3; ch++)
                                gfield[gc*3+(size_t)ch] = (float)(((double)gfield[gc*3+(size_t)ch]*(double)cnt
                                                                  + (double)m->field[sc*3+(size_t)ch]) / den);
                            if(gphase&&m->phase){
                                gphase[gc]=(float)(((double)gphase[gc]*(double)cnt+
                                                   incoming_phase)/den);
                            }
                            if (cnt < 65535u) gcount[gc] = (uint16_t)(cnt + 1u);
                        }
                        if(incoming_trusted)gtrusted[gc]=1;
                        continue;
                    }
                    if(gtrusted[gc]&&!incoming_trusted)continue;
                }
                gdomain[gc] = 1;
                if(incoming_trusted)gtrusted[gc]=1;
                if (gcount) gcount[gc] = 1;
                gfield[gc * 3 + 0] = m->field[sc * 3 + 0];
                gfield[gc * 3 + 1] = m->field[sc * 3 + 1];
                gfield[gc * 3 + 2] = m->field[sc * 3 + 2];
                if(gphase)gphase[gc]=(float)incoming_phase;
            }
        }
        if (overlap_samples)
            fprintf(stderr, "[%s] fitted overlap samples=%zu compatible=%zu conflict=%zu discarded_cells=%zu mean_gap=%.2f max_gap=%.2f vox\n",
                    giant_ribbon ? "giant" : "component",
                    overlap_samples, overlap_compatible, overlap_conflict,
                    overlap_conflict_cells,
                    overlap_samples ? overlap_gap_sum/(double)overlap_samples : 0.0,
                    overlap_gap_max);

        if(qopts.cylindrical_fill) {
            PhaseGateStats pg={0};
            if(quarantine_phase_gauge(
                    scratch,gdomain,gtrusted,gfield,gphase,H,W,
                    qopts.axis_y,qopts.axis_x,&pg)!=0) {
                fprintf(stderr,
                        "lifted-phase gauge preflight failed for component %zu "
                        "(group %ld)\n",g,groups[g]);
                return 1;
            }
        }

        if (bake_refit_masks) {
            size_t mask_requested = 0, mask_removed = 0;
            int mask_rc = apply_bake_refit_mask(scratch, bake_refit_masks, groups[g],
                                                gdomain, gtrusted, H, W,
                                                &mask_requested, &mask_removed);
            if (mask_rc < 0) return 1;
            (void)mask_requested; (void)mask_removed;
        }

        size_t topology_cells=0,initializer_cells=0,trusted_cells=0;
        int mask_invariant=1;
        for(size_t k=0;k<HW;k++) {
            topology_cells+=gtopology[k]!=0;
            initializer_cells+=gdomain[k]!=0;
            trusted_cells+=gtrusted[k]!=0;
            if((gdomain[k]&&!gtopology[k])||(gtrusted[k]&&!gdomain[k]))
                mask_invariant=0;
        }
        if(!mask_invariant) {
            fprintf(stderr,"[component] internal mask invariant failed for group %ld\n",
                    groups[g]);
            return 1;
        }
        fprintf(stderr,
                "[topology] component=%zu group=%ld footprint/initializer/trusted="
                "%zu/%zu/%zu grid-cells; movable-or-refit=%zu\n",
                g,groups[g],topology_cells,initializer_cells,trusted_cells,
                topology_cells>trusted_cells?topology_cells-trusted_cells:0u);

        float *gv = NULL, *guv = NULL; int32_t *gf = NULL; uint8_t *gfl = NULL;
        size_t gnv = 0, gnf = 0;
        if (QuadStrip_build_topology_with_phase(
                scratch,gtopology,gdomain,gfield,gphase,H,W,gr0,c0,
                grid_du,&qopts,&gv,&gnv,&gf,&gnf,&guv,&gfl) != 0) { continue; }

        /* Reclassify from direct high-confidence observations.  Interpolated
         * initializer samples remain part of the topology but are movable. */
        size_t hard_trusted = 0;
        for (size_t v = 0; v < gnv; v++) {
            int cc = (int)lround((double)guv[v*2] / grid_du) - c0;
            int rr = (int)lround((double)guv[v*2+1]) - gr0;
            int trust = rr >= 0 && rr < H && cc >= 0 && cc < W &&
                        gtrusted[(size_t)rr*(size_t)W+(size_t)cc] != 0;
            gfl[v] = (uint8_t)(trust ? 0 : 1);
            hard_trusted += (size_t)trust;
        }
        fprintf(stderr,"[trust] component=%zu group=%ld hard-certified=%zu/%zu; "
                       "movable-initializer=%zu\n",
                g,groups[g],hard_trusted,gnv,gnv-hard_trusted);

        double solve_grid_du = grid_du;
        if (ribbon_u_stride > 1 || ribbon_v_stride > 1) {
            RibbonCoarsenStats cs={0};
            if (coarsen_ribbon(
                    scratch,ribbon_u_stride,ribbon_v_stride,grid_du,
                    &gv,&gnv,&gf,&gnf,&guv,&gfl,&gphase,&H,&W,&cs)!=0) {
                fprintf(stderr,
                    "ERROR: confidence-preserving ribbon coarsening failed "
                    "for component %zu (group %ld)\n",g,groups[g]);
                return 1;
            }
            HW=(size_t)H*(size_t)W;
            solve_grid_du=grid_du*(double)ribbon_u_stride;
            fprintf(stderr,
                "[coarsen] component=%zu group=%ld stride UxV=%dx%d; "
                "grid %dx%d -> %dx%d, vertices %zu -> %zu (%.2f%%), "
                "faces %zu -> %zu; fine direct samples=%zu, coarse "
                "with-direct/trusted=%zu/%zu\n",
                g,groups[g],ribbon_u_stride,ribbon_v_stride,
                cs.fine_h,cs.fine_w,cs.coarse_h,cs.coarse_w,
                cs.fine_vertices,cs.coarse_vertices,
                100.0*(double)cs.coarse_vertices/(double)cs.fine_vertices,
                cs.fine_faces,cs.coarse_faces,cs.direct_samples,
                cs.coarse_with_direct,cs.coarse_trusted);
        }

        /* Keep the data term tied to the same upstream UV->XYZ observation on
         * every proximal pass.  Using each pass's already-moved input as its
         * new anchor would turn a high soft weight into cumulative drift. */
        float *fitted_source_positions=ARENA_ALLOC(
            scratch,gnv*3*sizeof *fitted_source_positions);
        memcpy(fitted_source_positions,gv,gnv*3*sizeof *fitted_source_positions);

        size_t expected_structured_faces = H > 1 && W > 1
            ? 2u * (size_t)(H - 1) * (size_t)(W - 1) : 0u;
        int structured_mesh = gnv == HW && gnf == expected_structured_faces;

        if(qopts.cylindrical_fill&&structured_mesh&&!v92_snap_policy) {
            TurnOrderRepairStats tr={0};
            ClothCollisionStats cloth={0};
            float *turn_rest=(float *)malloc(gnv*3*sizeof *turn_rest);
            if(!turn_rest){fprintf(stderr,"ERROR: turn-rest allocation failed\n");return 1;}
            memcpy(turn_rest,gv,gnv*3*sizeof *turn_rest);
            double long_u=fmax(100.0,64.0*grid_du);
            if(repair_cylindrical_turn_order(
                    gv,gnv,gf,gnf,guv,gphase,gfl,H,W,solve_grid_du,
                    qopts.axis_y,qopts.axis_x,long_u,&tr)!=0) {
                free(turn_rest);
                fprintf(stderr,
                    "ERROR: exact cylindrical turn-order transaction failed "
                    "for component %zu (group %ld)\n",g,groups[g]);
                return 1;
            }
            if(resolve_cylindrical_self_collisions(
                    gv,gnv,gf,gnf,guv,gphase,gfl,turn_rest,H,W,solve_grid_du,
                    qopts.axis_y,qopts.axis_x,long_u,
                    require_collision_free,collision_patience,collision_rounds,
                    collision_collar,collision_displacement_bound,
                    collision_settle_rounds,collision_settle_beta,
                    &cloth)!=0) {
                free(turn_rest);
                fprintf(stderr,
                    "ERROR: phase-ordered cloth collision transaction failed "
                    "for component %zu (group %ld)\n",g,groups[g]);
                return 1;
            }
            free(turn_rest);
        }

        if (metric_arap) {
            if (!structured_mesh) {
                fprintf(stderr, "ERROR: metric ARAP requires the complete row-major %dx%d rectangle; got %zu vertices\n",
                        H, W, gnv);
                return 1;
            }
            QuadStripArapStats astats;
            int arc=QuadStrip_metric_arap(gv,guv,gfl,H,W,&arap_opts,&astats);
            if (arc != 0) {
                fprintf(stderr, "ERROR: metric ARAP failed for component %zu (group %ld), code %d; no downstream geometry is permitted\n",
                        g,groups[g],arc);
                return 1;
            }
            fprintf(stderr,
                    "[metric-arap] component=%zu group=%ld iterations=%d converged=%s "
                    "max/rms-move=%.9g/%.9g edge-log-rms %.6g->%.6g "
                    "edge-log-max %.6g->%.6g local-barrier-rounds=%d min-step=%.6g "
                    "zero-fallbacks=%d barrier-stalled=%s "
                    "MG-L%d cycles=%d last=%.3g->%.3g max-cg=%d relres=%.3g "
                    "NLMG-coarse-levels=%d coarse-iters=%d caps=%d min-prolong=%.6g\n",
                    g,groups[g],astats.iterations,astats.converged?"yes":"no",
                    astats.maximum_vertex_movement,astats.rms_vertex_movement,
                    astats.initial_edge_log_rms,astats.final_edge_log_rms,
                    astats.initial_edge_log_max,astats.final_edge_log_max,
                    astats.barrier_rejections,astats.minimum_step_scale,
                    astats.barrier_zero_fallbacks,astats.barrier_stalled?"yes":"no",
                    astats.multigrid_levels_used,
                    astats.multigrid_cycles,astats.final_multigrid_residual_before,
                    astats.final_multigrid_residual_after,astats.maximum_cg_iterations,
                    astats.final_cg_relative_residual,
                    astats.nonlinear_coarse_levels,astats.nonlinear_coarse_iterations,
                    astats.nonlinear_coarse_caps,astats.nonlinear_min_prolongation_scale);
            if (arap_report) {
                if(!arap_report_first)fputs(",\n",arap_report);
                arap_report_first=0;
                fprintf(arap_report,
                        "  {\"component_index\":%zu,\"group\":%ld,"
                        "\"representative_track\":%ld,\"height\":%d,\"width\":%d,"
                        "\"vertices\":%zu,\"iterations\":%d,\"converged\":%s,"
                        "\"stop_reason\":%d,"
                        "\"effective_movement_tolerance\":%.17g,"
                        "\"barrier_rejections\":%d,\"barrier_zero_fallbacks\":%d,"
                        "\"barrier_stalled\":%s,\"minimum_step_scale\":%.17g,"
                        "\"maximum_cg_iterations\":%d,"
                        "\"multigrid_levels_used\":%d,\"multigrid_cycles\":%d,"
                        "\"nonlinear_coarse_levels\":%d,"
                        "\"nonlinear_coarse_iterations\":%d,"
                        "\"nonlinear_coarse_caps\":%d,"
                        "\"nonlinear_min_prolongation_scale\":%.17g,"
                        "\"final_multigrid_residual_before\":%.17g,"
                        "\"final_multigrid_residual_after\":%.17g,"
                        "\"final_cg_relative_residual\":%.17g,"
                        "\"maximum_vertex_movement\":%.17g,"
                        "\"rms_vertex_movement\":%.17g,"
                        "\"initial_edge_log_rms\":%.17g,\"final_edge_log_rms\":%.17g,"
                        "\"initial_edge_log_max\":%.17g,\"final_edge_log_max\":%.17g}",
                        g,groups[g],mem[0].tid,H,W,gnv,astats.iterations,
                        astats.converged?"true":"false",
                        astats.stop_reason,astats.effective_movement_tolerance,
                        astats.barrier_rejections,
                        astats.barrier_zero_fallbacks,
                        astats.barrier_stalled?"true":"false",
                        astats.minimum_step_scale,astats.maximum_cg_iterations,
                        astats.multigrid_levels_used,
                        astats.multigrid_cycles,astats.nonlinear_coarse_levels,
                        astats.nonlinear_coarse_iterations,astats.nonlinear_coarse_caps,
                        astats.nonlinear_min_prolongation_scale,
                        astats.final_multigrid_residual_before,
                        astats.final_multigrid_residual_after,
                        astats.final_cg_relative_residual,
                        astats.maximum_vertex_movement,astats.rms_vertex_movement,
                        astats.initial_edge_log_rms,astats.final_edge_log_rms,
                        astats.initial_edge_log_max,astats.final_edge_log_max);
            }
            arap_component_total++;
            arap_converged_total+=(size_t)(astats.converged!=0);
            arap_iteration_total+=(size_t)astats.iterations;
            arap_barrier_total+=(size_t)astats.barrier_rejections;
            arap_barrier_zero_total+=(size_t)astats.barrier_zero_fallbacks;
            arap_barrier_stall_total+=(size_t)(astats.barrier_stalled!=0);
            arap_cap_total+=(size_t)(!astats.converged&&
                                     astats.iterations>=arap_opts.max_iterations);
        }

        /* Components are disconnected by definition, so relax each one while it
         * is still a compact scratch-resident "snip".  The old path appended the
         * whole scroll first, then allocated another whole-scroll solver state;
         * it produced the same disconnected mathematics at much higher peak RAM.
         * RECT output is row-major HxW and takes the structured four-neighbour
         * red/black path.  SPAN remains supported through the face-CSR fallback. */
        int publish_component=1;
        int final_fold_hint_valid=0;
        size_t final_fold_hint=0;
        if (raw_cubes) {
            int gsnap = 0, gnopap = 0, snap_ok = 1;
            int grid_h = structured_mesh ? H : 0;
            int grid_w = grid_h ? W : 0;
            float *accepted_v = (float *)malloc(gnv*3*sizeof *accepted_v);
            uint8_t *accepted_fl = (uint8_t *)malloc(gnv*sizeof *accepted_fl);
            QuadStripSnapStats accepted_stats={0}, first_stats={0};
            int have_accepted = 0;
            int seam_refine_active = 0;
            if (!accepted_v || !accepted_fl) {
                fprintf(stderr, "warning: cannot allocate CT-pass rollback checkpoint for "
                                "component %zu (%zu vertices)\n", g, gnv);
                free(accepted_v); free(accepted_fl); snap_ok = 0;
            } else {
                memcpy(accepted_v, gv, gnv*3*sizeof *accepted_v);
                memcpy(accepted_fl, gfl, gnv*sizeof *accepted_fl);
            }
            /* CT-occupancy accumulators for the continuation trim, OR-ed
             * across every pass (including rejected ones). */
            uint8_t *trim_lit = NULL, *trim_cand = NULL;
            if (continuation_trim && structured_mesh && H > 1 && W > 1) {
                size_t ncell = (size_t)(H - 1) * (size_t)(W - 1);
                trim_lit = (uint8_t *)calloc(ncell, 1);
                trim_cand = (uint8_t *)calloc(ncell, 1);
            }
            for (int pass = 0; snap_ok && pass < snap_passes; pass++) {
                IntersectionCleanupStats ix_before={0},ix_after={0};
                Arena_T snap_arena = Arena_new();
                QuadStripSnapOpts pass_opts = sopts;
                /* Pass one discovers recto support.  Subsequent passes softly
                 * anchor those recovered filled sites and concentrate motion on
                 * the unresolved active set after the bounded remesh. */
                pass_opts.filled_target_anchor = pass > 0 ? refine_anchor : 0.0;
                pass_opts.fitted_source_positions = fitted_source_positions;
                pass_opts.cell_lit_accum = trim_lit;
                pass_opts.cell_candidate_accum = trim_cand;
                if (seam_refine_active &&
                    seam_refine_multiplier > pass_opts.quilt_seam_smooth)
                    pass_opts.quilt_seam_smooth = seam_refine_multiplier;
                fprintf(stderr, "[component] CT proximal pass %d/%d\n",
                        pass + 1, snap_passes);
                int gremesh = 0;
                QuadStripSnapStats gstats;
                int rc = QuadStripSnap_run(snap_arena, gv, gnv, gf, gnf, gfl, guv,
                                           grid_h, grid_w, raw_cubes, chunk, &pass_opts,
                                           &gsnap, &gnopap, &gremesh, &gstats);
                Arena_dispose(&snap_arena);
                if(rc==0) {
                    /* QuadStripSnap_run transactionally exact-audits its input
                     * and every committed state.  Reuse that proof instead of
                     * rebuilding the same whole-component BVHs before and after
                     * every proximal pass.  Keep a fallback for callers that
                     * explicitly disable the internal intersection guard. */
                    if(!gstats.intersection_audit_complete) {
                        if(snap_intersection_audit(accepted_v,gnv,gf,gnf,&ix_before)!=0||
                           snap_intersection_audit(gv,gnv,gf,gnf,&ix_after)!=0) {
                            fprintf(stderr,"ERROR: exact fallback pass audit failed "
                                           "for component %zu pass %d\n",g,pass+1);
                            rc=QUAD_STRIP_SNAP_ERROR;
                        } else {
                            gstats.intersection_audit_complete=1;
                            gstats.input_intersections=ix_before.conflicts;
                            gstats.output_intersections=ix_after.conflicts;
                            gstats.input_overlap=ix_before.overlap_pairs;
                            gstats.input_stab=ix_before.stab_pairs;
                            gstats.input_fold=ix_before.fold_pairs;
                            gstats.output_overlap=ix_after.overlap_pairs;
                            gstats.output_stab=ix_after.stab_pairs;
                            gstats.output_fold=ix_after.fold_pairs;
                        }
                    }
                    if(rc==0) {
                        if(pass==0) {
                            final_fold_hint_valid=1;
                            final_fold_hint=gstats.input_fold;
                        }
                        fprintf(stderr,
                            "[component] exact intersections pass %d: %zu -> %zu "
                            "[overlap %zu->%zu, stab %zu->%zu, fold %zu->%zu]\n",
                            pass+1,gstats.input_intersections,
                            gstats.output_intersections,
                            gstats.input_overlap,gstats.output_overlap,
                            gstats.input_stab,gstats.output_stab,
                            gstats.input_fold,gstats.output_fold);
                    }
                }
                if (rc != 0) {
                    snap_report_entry(snap_report,&snap_report_first,g,groups[g],
                                      mem[0].tid,pass+1,0,
                                      rc == QUAD_STRIP_SNAP_RAW_INCOMPLETE
                                          ? "RAW coverage incomplete"
                                          : "proximal pass failed",
                                      &gstats);
                    memcpy(gv, accepted_v, gnv*3*sizeof *gv);
                    memcpy(gfl, accepted_fl, gnv*sizeof *gfl);
                    if (sopts.require_complete_raw) {
                        fprintf(stderr,
                                "ERROR: aborting solid ribbon: strict RAW preflight "
                                "failed for component %zu (group %ld)\n",
                                g, groups[g]);
                        free(accepted_v); free(accepted_fl);
                        return 1;
                    }
                    if (have_accepted) {
                        gsnap=accepted_stats.snapped;gnopap=accepted_stats.nopap;
                        fprintf(stderr, "[component] pass %d failed; restored accepted pass %d\n",
                                pass+1,pass);
                    } else {
                        snap_ok=0;
                    }
                    break;
                }
                char decision_reason[256]="";
                if (!have_accepted) {
                    char why[256];
                    int accept=snap_first_pass_accept(&gstats,why,sizeof why);
                    fprintf(stderr,
                        "[component] %s pass 1 vs initializer: %s; Sander p95/p99 "
                        "%.3f/%.3f -> %.3f/%.3f, >=4x %.4f%% -> %.4f%%; "
                        "fixed-window dark fit/fill %zu/%zu -> %zu/%zu\n",
                        accept?"ACCEPT":"REJECT",why,
                        gstats.input_sander_p95,gstats.input_sander_p99,
                        gstats.sander_p95,gstats.sander_p99,
                        100.0*gstats.input_over4_fraction,100.0*gstats.over4_fraction,
                        gstats.input_bake_fitted_dark,gstats.input_bake_filled_dark,
                        gstats.output_bake_fitted_dark,gstats.output_bake_filled_dark);
                    if(!accept){
                        snap_report_entry(snap_report,&snap_report_first,g,groups[g],
                            mem[0].tid,1,0,why,&gstats);
                        memcpy(gv,accepted_v,gnv*3*sizeof *gv);
                        memcpy(gfl,accepted_fl,gnv*sizeof *gfl);
                        gsnap=gnopap=0;snap_pass_rejected_total++;break;
                    }
                    snprintf(decision_reason,sizeof decision_reason,"%s",why);
                }
                if (have_accepted) {
                    char why[256];
                    int accept=snap_transaction_accept(&accepted_stats,&first_stats,
                                                       &gstats,pass+1,why,sizeof why);
                    fprintf(stderr,
                            "[component] %s pass %d transaction: %s; Sander p95/p99 "
                            "%.3f/%.3f -> %.3f/%.3f, >=4x %.4f%% -> %.4f%%; "
                            "fill CT p25/p50/p75 %d/%d/%d -> %d/%d/%d; quilt %.3f -> %.3f\n",
                            accept?"ACCEPT":"REJECT",pass+1,why,
                            accepted_stats.sander_p95,accepted_stats.sander_p99,
                            gstats.sander_p95,gstats.sander_p99,
                            100.0*accepted_stats.over4_fraction,100.0*gstats.over4_fraction,
                            accepted_stats.filled_ct_p25,accepted_stats.filled_ct_p50,
                            accepted_stats.filled_ct_p75,gstats.filled_ct_p25,
                            gstats.filled_ct_p50,gstats.filled_ct_p75,
                            accepted_stats.quilt_edge_rms,gstats.quilt_edge_rms);
                    fprintf(stderr,
                            "[component] fixed-window bake dark fit/fill %zu/%zu -> %zu/%zu; "
                            "raw candidates/applied/quilt-only/no-candidate %d/%d/%d/%d\n",
                            accepted_stats.output_bake_fitted_dark,
                            accepted_stats.output_bake_filled_dark,
                            gstats.output_bake_fitted_dark,
                            gstats.output_bake_filled_dark,
                            gstats.raw_ridge,gstats.raw_candidate_applied,
                            gstats.quilt_only_applied,gstats.raw_no_candidate);
                    if (!accept) {
                        snap_report_entry(snap_report,&snap_report_first,g,groups[g],
                                          mem[0].tid,pass+1,0,why,&gstats);
                        memcpy(gv, accepted_v, gnv*3*sizeof *gv);
                        memcpy(gfl, accepted_fl, gnv*sizeof *gfl);
                        gsnap=accepted_stats.snapped;gnopap=accepted_stats.nopap;
                        snap_pass_rejected_total++;
                        break;
                    }
                    snprintf(decision_reason,sizeof decision_reason,"%s",why);
                }
                if (!have_accepted) first_stats=gstats;
                accepted_stats=gstats;have_accepted=1;
                final_fold_hint_valid=gstats.intersection_audit_complete;
                final_fold_hint=gstats.output_fold;
                memcpy(accepted_v, gv, gnv*3*sizeof *accepted_v);
                memcpy(accepted_fl, gfl, gnv*sizeof *accepted_fl);
                snap_pass_accepted_total++;
                snap_report_entry(snap_report,&snap_report_first,g,groups[g],
                                  mem[0].tid,pass+1,1,decision_reason,&gstats);
                if (!seam_refine_active && pass + 1 < snap_passes &&
                    seam_refine_multiplier > pass_opts.quilt_seam_smooth) {
                    double seam_ratio=0.0,seam_excess=0.0;
                    if (snap_quilt_seam_needs_refine(&gstats,&seam_ratio,&seam_excess)) {
                        seam_refine_active=1;
                        fprintf(stderr,
                                "[component] measured fitted/filled quilt seam %.2fx "
                                "(+%.3f vox, %zu edges); use %.2fx mixed-edge coupling "
                                "from pass %d\n",
                                seam_ratio,seam_excess,
                                gstats.quilt_applied_regions.fitted_filled_edges,
                                seam_refine_multiplier,pass+2);
                    }
                }
                printf("component %zu (group %ld) pass %d/%d: CT-relax %d ridge, %d no-ridge\n",
                       g, groups[g], pass + 1, snap_passes, gsnap, gnopap);
                if (gstats.displacement_rms < 0.01 && pass + 1 < snap_passes) {
                    fprintf(stderr, "[component] converged after pass %d: %.4f-voxel RMS "
                                    "movement; stopping iteration\n",
                            pass+1,gstats.displacement_rms);
                    break;
                }
                if (gremesh && pass + 1 < snap_passes && remesh_sweeps > 0) {
                    QuadStripRemeshStats rstats = {0};
                    int rrc = grid_h
                        ? QuadStripSnap_remesh_structured(gv, gfl, guv,
                                                         grid_h, grid_w,
                                                         remesh_sweeps, &rstats)
                        : -1;
                    if (rrc == 0 && rstats.applied_sweeps > 0) {
                        if(sopts.preserve_axial)
                            for(size_t v=0;v<gnv;v++)gv[v*3]=accepted_v[v*3];
                        fprintf(stderr,
                                "[component] correspondence remesh applied %d/%d sweep(s) "
                                "before pass %d: %s metric, %s-first, alpha %.5g; "
                                "Sander p95 %.6f->%.6f p99 %.6f->%.6f, "
                                ">=2 %zu->%zu >=4 %zu->%zu, edge-log RMS %.6f->%.6f\n",
                                rstats.applied_sweeps, remesh_sweeps, pass + 2,
                                rstats.used_robust_metric ? "robust" : "mean",
                                rstats.used_v_first ? "v" : "u", rstats.alpha,
                                rstats.before_p95, rstats.after_p95,
                                rstats.before_p99, rstats.after_p99,
                                rstats.before_over2, rstats.after_over2,
                                rstats.before_over4, rstats.after_over4,
                                rstats.before_edge_log_rms,
                                rstats.after_edge_log_rms);
                    } else if (rrc == 0) {
                        fprintf(stderr,
                                "[component] correspondence gate fired, but every shared "
                                "map regressed an emitted-triangle objective; retaining the "
                                "current lattice for pass %d\n",
                                pass + 2);
                    } else {
                        memcpy(gv, accepted_v, gnv*3*sizeof *gv);
                        memcpy(gfl, accepted_fl, gnv*sizeof *gfl);
                        fprintf(stderr, "[component] remesh gate fired but no supported "
                                        "remesh was available; stopping after pass %d\n",
                                pass + 1);
                        break;
                    }
                } else if (gremesh && pass + 1 < snap_passes) {
                    fprintf(stderr, "[component] remesh gate fired but remeshing is "
                                    "disabled; continuing to pass %d on the current lattice\n",
                            pass + 2);
                }
            }
            if(post_snap_arap&&snap_ok&&have_accepted){
                PostSnapArapTxn pa;
                fprintf(stderr,"[component] post-snap metric ARAP: all %zu current "
                               "vertices are soft CT anchors; preserve-z=%s\n",
                        gnv,arap_opts.preserve_axial?"yes":"no");
                int prc=post_snap_arap_relax(gv,gnv,gf,gnf,guv,grid_h,grid_w,
                                              &arap_opts,&pa);
                if(prc!=0){
                    memcpy(gv,accepted_v,gnv*3*sizeof *gv);
                    memset(&pa,0,sizeof pa);pa.attempted=1;pa.solver_code=-1;
                    fprintf(stderr,"[component] post-snap metric ARAP transaction "
                                   "could not be audited; restored accepted CT snap\n");
                }else{
                    fprintf(stderr,
                        "[component] %s post-snap metric ARAP: solver rc=%d iters=%d "
                        "converged=%s, selected alpha %.6g, move rms/max %.6g/%.6g, "
                        "axial drift %.9g\n",
                        pa.accepted?"ACCEPT":"REJECT",pa.solver_code,
                        pa.solver.iterations,pa.solver.converged?"yes":"no",
                        pa.selected_scale,pa.displacement_rms,pa.displacement_max,
                        pa.axial_drift_max);
                    fprintf(stderr,
                        "[component] post-snap edge-log RMS %.6g -> %.6g "
                        "(full solve %.6g); exact conflicts %zu -> %zu "
                        "[overlap %zu->%zu, stab %zu->%zu, fold %zu->%zu]; "
                        "%d exact trial audits, %d orientation rejects, "
                        "%zu locally rolled-back vertices\n",
                        pa.input_edge_log_rms,pa.selected_edge_log_rms,
                        pa.solver.final_edge_log_rms,
                        pa.input_intersections,pa.output_intersections,
                        pa.input_overlap,pa.output_overlap,pa.input_stab,
                        pa.output_stab,pa.input_fold,pa.output_fold,
                        pa.trial_audits,pa.orientation_rejects,
                        pa.locally_rolled_back_vertices);
                    if(pa.accepted) {
                        memcpy(accepted_v,gv,gnv*3*sizeof *accepted_v);
                        final_fold_hint_valid=1;
                        final_fold_hint=pa.output_fold;
                    }
                    else memcpy(gv,accepted_v,gnv*3*sizeof *gv);
                }
                snap_post_arap_report_entry(snap_report,&snap_report_first,g,
                                            groups[g],mem[0].tid,&pa);
            }
            if (continuation_trim && snap_ok && have_accepted &&
                trim_lit != NULL && trim_cand != NULL) {
                TrimStats ts;
                trim_black_continuation(gfl, trim_lit, trim_cand, H, W,
                                        gf, &gnf,
                                        (size_t)trim_min_region,
                                        trim_rim_guard, &ts);
                fprintf(stderr,
                    "[trim] group %ld: black cells=%zu regions=%zu "
                    "trimmed=%zu (kept: small-region %zu, rim-guard %zu); "
                    "%zu faces removed, %zu remain\n",
                    groups[g], ts.black_cells, ts.regions, ts.trimmed_cells,
                    ts.kept_small_region, ts.kept_rim_guard,
                    ts.faces_removed, gnf);
                trim_cells_total += ts.trimmed_cells;
                trim_faces_total += ts.faces_removed;
            }
            free(trim_lit); free(trim_cand);
            trim_lit = trim_cand = NULL;
            free(accepted_v); free(accepted_fl);
            if (snap_ok && have_accepted) {
                snap_target_total += (size_t)gsnap;
                snap_nopap_total += (size_t)gnopap;
            } else if(!snap_ok) {
                fprintf(stderr, "warning: CT-relax failed for component %zu "
                                "(group %ld, raw=%s)\n", g, groups[g], raw_cubes);
                publish_component=0;
            } else {
                fprintf(stderr,"[component] retained exact-audited frozen initializer "
                               "after rejecting pass 1 for group %ld\n",groups[g]);
            }
        }

        if(!publish_component) {
            fprintf(stderr,"[component] quarantine group %ld after a solver or "
                           "input-integrity failure; output remains alpha\n",groups[g]);
            skipped++;
            continue;
        }

        if (repair_folds && repair_final_fold_flaps(
                gv, gnv, gf, &gnf, final_fold_hint_valid, final_fold_hint,
                "group", groups[g], &fold_repair_total) != 0) {
            fprintf(stderr,
                    "ERROR: refusing to publish group %ld after fold repair failure\n",
                    groups[g]);
            return 1;
        }

        if (snip_manifest && emit_component_artifact(
                snips_dir, g, groups[g], gv, gnv, gf, gnf, guv) != 0) {
            fprintf(stderr,
                    "failed to emit component artifact %zu (group %ld)\n",
                    g, groups[g]);
            fclose(snip_manifest); return 1;
        }
        if (snip_manifest && emit_rect_snips(snips_dir, snip_manifest,
                                             &snip_manifest_first, &snip_serial,
                                             snip_cols, g, groups[g], mem[0].tid,
                                             H, W, gr0, c0, grid_du,
                                             ribbon_v_stride,ribbon_u_stride,
                                             gv, gnv, gf, gnf, guv, gfl,
                                             &streamed_nv, &streamed_nf) != 0) {
            fprintf(stderr, "failed to emit snips for component %zu (group %ld)\n",
                    g, groups[g]);
            fclose(snip_manifest); return 1;
        }

        /* Capture the common material lattice before shelf packing.  This is
         * the only coordinate system shared by the branch-safe geometry lane
         * and the deliberately coarser v92 bake-proxy lane. */
        if (correspondence_path && !omit_mesh) {
            long group_material = LONG_MIN;
            for (size_t i = 0; i < ntr; i++) if (grp_of[i] == groups[g]) {
                if (!material_explicit[i] || material_of[i] < INT32_MIN ||
                    material_of[i] > INT32_MAX ||
                    (group_material != LONG_MIN &&
                     group_material != material_of[i])) {
                    fprintf(stderr,
                            "cannot key group %ld by one explicit material lineage\n",
                            groups[g]);
                    return 1;
                }
                group_material = material_of[i];
            }
            if (group_material == LONG_MIN) {
                fprintf(stderr, "cannot determine material lineage for group %ld\n",
                        groups[g]);
                return 1;
            }
            fv_reserve(&corr_uv, gnv * 2u);
            memcpy(corr_uv.v + corr_uv.n, guv, gnv * 2u * sizeof *guv);
            corr_uv.n += gnv * 2u;
            iv_reserve(&corr_namespace, gnv);
            iv_reserve(&corr_lineage, gnv);
            for (size_t t = 0; t < gnv; t++) {
                corr_namespace.v[corr_namespace.n++] = 0;
                corr_lineage.v[corr_lineage.n++] = (int32_t)group_material;
            }
        }

        /* The snip manifest above retains material/global UV.  The monolithic
         * review atlas may instead shelf-pack disconnected charts so physically
         * distinct wraps never multi-cover or mean-blend in obj_bake_raw. */
        double atlas_u0=guv[0],atlas_u1=guv[0],atlas_v0=guv[1],atlas_v1=guv[1];
        for(size_t t=1;t<gnv;t++){
            double u=guv[t*2],v=guv[t*2+1];
            if(u<atlas_u0)atlas_u0=u;if(u>atlas_u1)atlas_u1=u;
            if(v<atlas_v0)atlas_v0=v;if(v>atlas_v1)atlas_v1=v;
        }
        if (pack_width > 0.0 && !omit_mesh) {
            double u0=atlas_u0,u1=atlas_u1,v0=atlas_v0,v1=atlas_v1;
            double pw=u1-u0,ph=v1-v0,gap=4.0;
            if(pack_u>0.0&&pack_u+pw>pack_width){
                pack_u=0.0;pack_v+=pack_row_h+gap;pack_row_h=0.0;
            }
            double su=pack_u-u0,sv=pack_v-v0;
            for(size_t t=0;t<gnv;t++){guv[t*2]+=(float)su;guv[t*2+1]+=(float)sv;}
            fprintf(stderr,"[pack] component %zu group=%ld box %.0fx%.0f at (u=%.0f,v=%.0f)\n",
                    g,groups[g],pw,ph,pack_u,pack_v);
            atlas_u0=pack_u;atlas_u1=pack_u+pw;
            atlas_v0=pack_v;atlas_v1=pack_v+ph;
            pack_u+=pw+gap;if(ph>pack_row_h)pack_row_h=ph;
        }
        if (snip_manifest) {
            atlas_written[g]=1;atlas_box_u0[g]=atlas_u0;atlas_box_u1[g]=atlas_u1;
            atlas_box_v0[g]=atlas_v0;atlas_box_v1[g]=atlas_v1;
            atlas_track[g]=mem[0].tid;atlas_height[g]=H;atlas_width[g]=W;
        }

        /* append with a running vertex-base offset on the faces */
        if (!omit_mesh) {
            size_t base = verts.n / 3;
            fv_reserve(&verts, gnv * 3);  memcpy(verts.v + verts.n, gv, gnv * 3 * sizeof *gv);  verts.n += gnv * 3;
            fv_reserve(&uv, gnv * 2);     memcpy(uv.v + uv.n, guv, gnv * 2 * sizeof *guv);       uv.n += gnv * 2;
            bv_reserve(&filled, gnv);     memcpy(filled.v + filled.n, gfl, gnv);                 filled.n += gnv;
            iv_reserve(&compv, gnv);      for (size_t t = 0; t < gnv; t++) compv.v[compv.n + t] = (int32_t)mem[0].tid;  compv.n += gnv;
            iv_reserve(&faces, gnf * 3);
            for (size_t t = 0; t < gnf * 3; t++) faces.v[faces.n + t] = gf[t] + (int32_t)base;
            faces.n += gnf * 3;
        }

        size_t gfit = 0; for (size_t t = 0; t < gnv; t++) gfit += (gfl[t] == 0);
        fitted_total += gfit;
        filled_total += gnv - gfit;
        used++;
    }
    Arena_dispose(&scratch);

    /* Small source domains never enter the fitted-track manifest.  Preserve
     * them as exact, disconnected observations in a complete whole-volume
     * assembly; a one-group diagnostic remains intentionally one group. */
    if (only_group == LONG_MIN && append_passthrough_components(
            report, aggregate, grid_du,
            u_column_range, u_column_first, u_column_end, u_column_halo,
            omit_mesh, repair_folds, snips_dir, snip_manifest,
            &snip_manifest_first, &snip_serial, snip_cols, ngroups,
            pack_width, &pack_u, &pack_v, &pack_row_h,
            &verts, &uv, &faces, &filled, &compv,
            correspondence_path ? &corr_uv : NULL,
            correspondence_path ? &corr_namespace : NULL,
            correspondence_path ? &corr_lineage : NULL,
            &streamed_nv, &streamed_nf,
            &passthrough_atlas, &passthrough_atlas_n,
            &passthrough_stats, &fold_repair_total) != 0) {
        if (snip_manifest) fclose(snip_manifest);
        if (snap_report) fclose(snap_report);
        if (arap_report) fclose(arap_report);
        return 1;
    }
    used += passthrough_stats.components;
    fitted_total += passthrough_stats.supported_vertices;
    filled_total += passthrough_stats.generated_vertices;

    if (arap_report) {
        fprintf(arap_report,
                "\n ],\n \"totals\":{\"components\":%zu,\"converged\":%zu,"
                "\"hit_iteration_cap\":%zu,\"iterations\":%zu,"
                "\"barrier_rejections\":%zu,\"barrier_zero_fallbacks\":%zu,"
                "\"components_with_barrier_stall\":%zu}\n}\n",
                arap_component_total,arap_converged_total,
                arap_cap_total,arap_iteration_total,arap_barrier_total,
                arap_barrier_zero_total,arap_barrier_stall_total);
        fclose(arap_report);
        printf("wrote metric ARAP report %s\n",arap_report_path);
    }

    if (snap_report) {
        fprintf(snap_report,
                "\n ],\n \"totals\":{\"accepted_passes\":%zu,\"rejected_passes\":%zu,"
                "\"initializer_quarantine\":{\"measured_triangles\":%zu,"
                "\"severe_triangles\":%zu,\"marked_tiles\":%zu,"
                "\"dilated_tiles\":%zu,\"removed_initializer\":%zu,"
                "\"removed_trusted\":%zu},"
                "\"passthrough\":{\"components\":%zu,\"vertices\":%zu,"
                "\"supported_vertices\":%zu,\"generated_vertices\":%zu,"
                "\"faces\":%zu},"
                "\"fold_repair\":{\"enabled\":%s,"
                "\"components_audited\":%zu,\"components_repaired\":%zu,"
                "\"input_folds\":%zu,\"faces_removed\":%zu,"
                "\"output_folds\":%zu}}\n}\n",
                snap_pass_accepted_total,snap_pass_rejected_total,
                init_quarantine_total.measured_triangles,
                init_quarantine_total.severe_triangles,
                init_quarantine_total.marked_tiles,
                init_quarantine_total.dilated_tiles,
                init_quarantine_total.removed_initializer,
                init_quarantine_total.removed_trusted,
                passthrough_stats.components,passthrough_stats.vertices,
                passthrough_stats.supported_vertices,
                passthrough_stats.generated_vertices,passthrough_stats.faces,
                repair_folds?"true":"false",
                fold_repair_total.components_audited,
                fold_repair_total.components_repaired,
                fold_repair_total.input_folds,
                fold_repair_total.faces_removed,
                fold_repair_total.output_folds);
        fclose(snap_report);
        printf("wrote snap transaction report %s\n",snap_report_path);
    }

    if (snip_manifest) {
        fputs("\n ],\n \"atlas_components\":[\n", snip_manifest);
        int atlas_first = 1;
        for (size_t g = 0; g < ngroups; g++) if (atlas_written[g]) {
            long atlas_material = -1;
            if (!atlas_first) fputs(",\n", snip_manifest);
            atlas_first = 0;
            for (size_t i = 0; i < ntr; i++)
                if (grp_of[i] == groups[g]) {
                    atlas_material = material_of[i];
                    break;
                }
            fprintf(snip_manifest,
                    "  {\"component_index\":%zu,\"group\":%ld,"
                    "\"representative_track\":%ld,\"material_identity\":",
                    g, groups[g], atlas_track[g]);
            fprintf(snip_manifest, "%ld,\"member_tracks\":[", atlas_material);
            {
                int member_first = 1;
                for (size_t i = 0; i < ntr; i++) if (grp_of[i] == groups[g]) {
                    if (!member_first) fputc(',', snip_manifest);
                    member_first = 0;
                    fprintf(snip_manifest, "%ld", tid_of[i]);
                }
            }
            fprintf(snip_manifest,
                    "],\"height\":%d,\"width\":%d,"
                    "\"artifact_obj\":\"component_%06zu_g%ld.obj\","
                    "\"artifact_vmesh\":\"component_%06zu_g%ld.vmesh\","
                    "\"u_min\":%.17g,\"u_max\":%.17g,"
                    "\"v_min\":%.17g,\"v_max\":%.17g}",
                    atlas_height[g], atlas_width[g],
                    g, groups[g], g, groups[g],
                    atlas_box_u0[g], atlas_box_u1[g], atlas_box_v0[g], atlas_box_v1[g]);
        }
        for (size_t i = 0; i < passthrough_atlas_n; i++) {
            const PassthroughAtlas *pa = &passthrough_atlas[i];
            if (!atlas_first) fputs(",\n", snip_manifest);
            atlas_first = 0;
            fprintf(snip_manifest,
                    "  {\"component_index\":%zu,\"group\":%ld,"
                    "\"representative_track\":null,\"topology_label\":%ld,"
                    "\"material_identity\":null,\"member_tracks\":[],"
                    "\"height\":%d,\"width\":%d,"
                    "\"artifact_obj\":\"component_%06zu_g%ld.obj\","
                    "\"artifact_vmesh\":\"component_%06zu_g%ld.vmesh\","
                    "\"u_min\":%.17g,\"u_max\":%.17g,"
                    "\"v_min\":%.17g,\"v_max\":%.17g}",
                    pa->component_index, -pa->label - 1L, pa->label,
                    pa->height, pa->width,
                    pa->component_index, -pa->label - 1L,
                    pa->component_index, -pa->label - 1L,
                    pa->u_min, pa->u_max, pa->v_min, pa->v_max);
        }
        fprintf(snip_manifest,
                "\n ],\n \"packing\":{\"packed\":%s,\"pack_width\":%.17g,\"gap\":4},\n"
                " \"totals\":{\"snips\":%zu,\"stored_vertices_with_shared_seams\":%zu,"
                "\"stored_faces\":%zu,\"fitted_vertices\":%zu,\"filled_vertices\":%zu}\n}\n",
                (pack_width > 0.0 && !omit_mesh) ? "true" : "false", pack_width,
                snip_serial, streamed_nv, streamed_nf, fitted_total, filled_total);
        fclose(snip_manifest);
        printf("wrote %zu ribbon snips to %s (manifest.json)\n", snip_serial, snips_dir);
    }
    free(passthrough_atlas);

    size_t total_nv = verts.n / 3, total_nf = faces.n / 3;
    if ((report_only ? used : (snips_only ? streamed_nv : total_nv)) == 0) {
        fprintf(stderr, "no components assembled\n"); return 1;
    }

    if (raw_cubes)
        printf("CT-relax total: %zu movable verts with quilt target, %zu movable "
               "without a raw ridge candidate; accepted/rejected proximal passes %zu/%zu\n",
               snap_target_total, snap_nopap_total,
               snap_pass_accepted_total,snap_pass_rejected_total);
    if (continuation_trim)
        printf("continuation trim: %zu black-fill cells deleted "
               "(%zu faces); glimmer and fitted geometry untouched\n",
               trim_cells_total, trim_faces_total);
    if (repair_folds)
        printf("fold-only repair: audited %zu components, repaired %zu; "
               "%zu fold pairs -> %zu, %zu flap faces removed\n",
               fold_repair_total.components_audited,
               fold_repair_total.components_repaired,
               fold_repair_total.input_folds,
               fold_repair_total.output_folds,
               fold_repair_total.faces_removed);

    if (!omit_mesh) {
        if (ObjIO_write_uv(out_obj, verts.v, total_nv, faces.v, total_nf, uv.v) != 0) {
            fprintf(stderr, "failed to write %s\n", out_obj); return 1;
        }
        char vmesh[1500];
        if (MeshBin_companion_path(out_obj, vmesh, sizeof vmesh) != 0) {
            fprintf(stderr, "failed to derive VMesh path from %s\n", out_obj);
            return 1;
        }
        if (MeshBin_write(vmesh, verts.v, total_nv, faces.v, total_nf, uv.v) != 0) {
            fprintf(stderr, "failed to write %s\n", vmesh);
            return 1;
        }
        if (correspondence_path && write_correspondence(
                correspondence_path, grid_du, total_nv,
                &corr_uv, &corr_namespace, &corr_lineage) != 0) {
            fprintf(stderr, "failed to write %s\n", correspondence_path);
            return 1;
        }
    }

    printf("used %zu components, skipped %zu; fitted cells %zu, hole-filled %zu (+%.2f%%)\n",
           used, skipped, fitted_total, filled_total,
           100.0 * (double)filled_total / (double)(fitted_total ? fitted_total : 1));
    if (!omit_mesh)
        printf("wrote %s: %zu verts, %zu faces\n", out_obj, total_nv, total_nf);

    /* topology / coverage PNG: per-component tint, filled cells green. */
    if (topo_png && omit_mesh)
        fprintf(stderr, "warning: --topology-map skipped when mesh output is disabled\n");
    if (topo_png && !omit_mesh && total_nv > 0) {
        long min_gc = 0, max_gc = 0, min_gr = 0, max_gr = 0;
        int first = 1;
        for (size_t t = 0; t < total_nv; t++) {
            long gc = lround((double)uv.v[t * 2 + 0] / grid_du);
            long gr = lround((double)uv.v[t * 2 + 1]);
            if (first) { min_gc = max_gc = gc; min_gr = max_gr = gr; first = 0; }
            if (gc < min_gc) min_gc = gc; if (gc > max_gc) max_gc = gc;
            if (gr < min_gr) min_gr = gr; if (gr > max_gr) max_gr = gr;
        }
        long Wc = max_gc - min_gc + 1, Hc = max_gr - min_gr + 1;
        if (Wc > 0 && Hc > 0 && (double)Wc * (double)Hc < 6.0e8) {
            uint8_t *rgb = (uint8_t *)calloc((size_t)Wc * (size_t)Hc * 3, 1);
            for(size_t f=0;f<total_nf;f++) {
                const int32_t *tri=&faces.v[f*3];
                topology_raster_triangle(rgb,Wc,Hc,min_gc,min_gr,grid_du,
                    uv.v,filled.v,tri,compv.v[(size_t)tri[0]]);
            }
            /* Preserve the exact classification at every actual vertex after
             * triangle interpolation paints the native lattice between them. */
            for (size_t t = 0; t < total_nv; t++) {
                long gc = lround((double)uv.v[t * 2 + 0] / grid_du) - min_gc;
                long gr = lround((double)uv.v[t * 2 + 1]) - min_gr;
                size_t o = ((size_t)gr * (size_t)Wc + (size_t)gc) * 3;
                if (filled.v[t]) { rgb[o] = 40; rgb[o + 1] = 210; rgb[o + 2] = 70; }
                else {
                    unsigned h = (unsigned)compv.v[t] * 2654435761u;
                    rgb[o] = (uint8_t)(60 + (h & 127u));
                    rgb[o + 1] = (uint8_t)(60 + ((h >> 8) & 127u));
                    rgb[o + 2] = (uint8_t)(60 + ((h >> 16) & 127u));
                }
            }
            if (VesPng_write_rgb(topo_png, rgb, (int)Wc, (int)Hc) == 0)
                printf("wrote %s (%ldx%ld)\n", topo_png, Wc, Hc);
            else
                fprintf(stderr, "warning: failed to write %s\n", topo_png);
            free(rgb);
        } else {
            fprintf(stderr, "topology map too large (%ldx%ld); skipped\n", Wc, Hc);
        }
    }

    free(verts.v); free(uv.v); free(faces.v); free(filled.v); free(compv.v);
    free(corr_uv.v); free(corr_namespace.v); free(corr_lineage.v);
    Arena_dispose(&cfg);
    return 0;
}
