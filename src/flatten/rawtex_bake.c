/* rawtex_bake.c -- see rawtex_bake.h. Extracted verbatim from obj_bake_raw.c
 * (the rasterizer + its private helpers), plus an optional face_skip mask. */
#include "rawtex_bake.h"
#include "../common/pipeline_constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint64_t key; size_t hits; } RawtexPair;
typedef struct {
    Arena_T arena;
    RawtexPair *table;
    size_t capacity,used,byte_limit,peak_payload;
} RawtexPairs;

static size_t rawtex_pair_slot(uint64_t key,size_t capacity)
{
    key^=key>>30; key*=UINT64_C(0xbf58476d1ce4e5b9);
    key^=key>>27; key*=UINT64_C(0x94d049bb133111eb);
    key^=key>>31;
    return (size_t)key&(capacity-1);
}

static int rawtex_pair_grow(RawtexPairs *pairs)
{
    size_t capacity=0;
    size_t limit=pairs->byte_limit/sizeof(RawtexPair),payload=0;
    Arena_T arena=NULL;
    RawtexPair *table=NULL;
    if (pairs->capacity>SIZE_MAX/2) return -1;
    capacity=pairs->capacity ? 2*pairs->capacity : 64;
    if (capacity>limit ||
        pairs->capacity>limit-capacity) return -1;
    /* Bound both live tables during rehash, not just the eventual table.
     * Arena chunk overhead and other raster buffers are separate. */
    payload=(capacity+pairs->capacity)*sizeof(RawtexPair);
    arena=Arena_new();
    table=(RawtexPair *)ARENA_CALLOC(arena,capacity,sizeof *table);
    for (size_t i=0;i<pairs->capacity;i++) if (pairs->table[i].key) {
        size_t slot=rawtex_pair_slot(pairs->table[i].key,capacity);
        while (table[slot].key) slot=(slot+1)&(capacity-1);
        table[slot]=pairs->table[i];
    }
    if (pairs->arena) Arena_dispose(&pairs->arena);
    pairs->arena=arena; pairs->table=table; pairs->capacity=capacity;
    if (payload>pairs->peak_payload) pairs->peak_payload=payload;
    return 0;
}

static int rawtex_pair_add(RawtexPairs *pairs,int32_t a,int32_t b)
{
    uint64_t key=0;
    size_t slot=0;
    if (!pairs || a<0 || b<0 || a==b) return -1;
    if (a>b) { int32_t swap=a; a=b; b=swap; }
    key=((uint64_t)((uint32_t)a+1)<<32)|((uint32_t)b+1);
    if (pairs->capacity) {
        slot=rawtex_pair_slot(key,pairs->capacity);
        while (pairs->table[slot].key && pairs->table[slot].key!=key)
            slot=(slot+1)&(pairs->capacity-1);
        if (pairs->table[slot].key) {
            if (pairs->table[slot].hits==SIZE_MAX) return -1;
            pairs->table[slot].hits++; return 0;
        }
    }
    if (pairs->used>=pairs->capacity/2) {
        if (rawtex_pair_grow(pairs)!=0) return -1;
        slot=rawtex_pair_slot(key,pairs->capacity);
        while (pairs->table[slot].key) slot=(slot+1)&(pairs->capacity-1);
    }
    pairs->table[slot]=(RawtexPair){key,1}; pairs->used++;
    return 0;
}

static int rawtex_pair_compare(const void *pa,const void *pb)
{
    const RawtexPair *a=(const RawtexPair *)pa,*b=(const RawtexPair *)pb;
    return a->key<b->key ? -1 : a->key>b->key;
}

/* The legacy diagnostic prints unordered pairs in source-label order. Compact
 * and sort the exact sparse counts in place only after rasterization ends. */
static void rawtex_pairs_report(RawtexPairs *pairs,size_t components)
{
    size_t count=0;
    for (size_t i=0;i<pairs->capacity;i++) if (pairs->table[i].key)
        pairs->table[count++]=pairs->table[i];
    if (count) qsort(pairs->table,count,sizeof *pairs->table,rawtex_pair_compare);
    for (size_t i=0;i<count;i++)
        fprintf(stderr,"    component conflict %zu-%zu support=%zu\n",
            (size_t)(pairs->table[i].key>>32)-1,
            (size_t)(uint32_t)pairs->table[i].key-1,pairs->table[i].hits);
    fprintf(stderr,"  component conflict graph: vertices=%zu edges=%zu\n",components,count);
    fprintf(stderr,"  component conflict storage: sparse pairs=%zu table=%zu peak-rehash-payload=%zu cap=%zu bytes\n",
        count,pairs->capacity*sizeof *pairs->table,pairs->peak_payload,pairs->byte_limit);
}

static int32_t rawtex_component_find(int32_t *parent, int32_t vertex)
{
    int32_t root = vertex;
    while (parent[root] >= 0) root = parent[root];
    while (vertex != root) {
        int32_t next = parent[vertex];
        parent[vertex] = root;
        vertex = next;
    }
    return root;
}

static void rawtex_component_union(int32_t *parent, int32_t a, int32_t b)
{
    a = rawtex_component_find(parent, a);
    b = rawtex_component_find(parent, b);
    if (a == b) return;
    if (parent[a] > parent[b]) { int32_t swap = a; a = b; b = swap; }
    parent[a] += parent[b];
    parent[b] = a;
}
#include <math.h>
#include <stdint.h>

#include "../common/tiff_io.h"
#include "../common/ves_png.h"
#include "../common/tif_strip.h"
#include "../common/raw_sample.h"

/* ------------------------------------------------------------------ util
 * Trivial malloc wrappers + the rasterizer's private math helpers. These are
 * intentionally file-local duplicates of the same statics in obj_bake_raw.c
 * (5-line wrappers not worth a shared alloc header). */
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

static double gray_of(double v, double lo, double hi)
{
    double g = (v - lo) / (hi - lo);
    if (g < 0.0) g = 0.0;
    if (g > 1.0) g = 1.0;
    return g;
}

static double face_sigma(const float *verts, RawtexUv uv,
                         size_t a, size_t b, size_t c, double *out_smin)
{
    /* Sander et al. 2001 singular values of the UV->3D map */
    double s1 = (double)Rawtex_uv(uv,a * 2 + 0), t1 = (double)Rawtex_uv(uv,a * 2 + 1);
    double s2 = (double)Rawtex_uv(uv,b * 2 + 0), t2 = (double)Rawtex_uv(uv,b * 2 + 1);
    double s3 = (double)Rawtex_uv(uv,c * 2 + 0), t3 = (double)Rawtex_uv(uv,c * 2 + 1);
    double A2 = (s2 - s1) * (t3 - t1) - (s3 - s1) * (t2 - t1);
    double Ss[3], St[3];
    double aa = 0.0, bb = 0.0, cc = 0.0, disc = 0.0;
    int k = 0;
    *out_smin = 0.0;
    if (fabs(A2) < 1e-12) return 1e6;   /* UV-degenerate: infinite compression */
    for (k = 0; k < 3; k++) {
        double q1 = (double)verts[a * 3 + k], q2 = (double)verts[b * 3 + k];
        double q3 = (double)verts[c * 3 + k];
        Ss[k] = (q1 * (t2 - t3) + q2 * (t3 - t1) + q3 * (t1 - t2)) / A2;
        St[k] = (q1 * (s3 - s2) + q2 * (s1 - s3) + q3 * (s2 - s1)) / A2;
    }
    aa = Ss[0] * Ss[0] + Ss[1] * Ss[1] + Ss[2] * Ss[2];
    bb = Ss[0] * St[0] + Ss[1] * St[1] + Ss[2] * St[2];
    cc = St[0] * St[0] + St[1] * St[1] + St[2] * St[2];
    disc = sqrt((aa - cc) * (aa - cc) + 4.0 * bb * bb);
    *out_smin = sqrt(((aa + cc - disc) > 0.0 ? (aa + cc - disc) : 0.0) / 2.0);
    return sqrt((aa + cc + disc) / 2.0);
}

static uint8_t log2_to_u8(double sigma)
{
    double t = 0.0;
    if (sigma < 1e-9) return 0;
    t = 128.0 + 42.5 * (log(sigma) / 0.6931471805599453);
    if (t < 0.0) t = 0.0;
    if (t > 255.0) t = 255.0;
    return (uint8_t)(t + 0.5);
}

static int cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

void Rawtex_stretch_window(const double *val, const uint8_t *has, size_t nv,
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

/* Per-axis raster cap (image dimension) and default total-pixel budget. */
#define RAWTEX_AXIS_CAP       ((size_t)1 << 20)
#define RAWTEX_MAX_PX_DEFAULT ((size_t)1 << 28)

/* Half-texel convention: pixel k covers [k*step, (k+1)*step] and samples at
 * its CENTER (k+0.5)*step.  The fitted ribbon emits uv on an exact grid_du
 * lattice; the old corner-sample convention put every pixel center ON a
 * lattice line, so one dropped grid cell painted a 1-px black crack tracing
 * every section boundary, and shared edges beat against the sampler. */
static size_t plan_dim(double span, double step)
{
    double cells = ceil(span / step);
    if (cells < 1.0) return 1;
    if (!isfinite(cells) || cells > (double)SIZE_MAX) return SIZE_MAX;
    return (size_t)cells;
}

int Rawtex_plan_field(RawtexUv uv, size_t nv, double du, double dv,
                       size_t max_px, const RawtexWindow *window,
                       RawtexPlan *out)
{
    size_t i = 0;
    double span_u = 0.0, span_v = 0.0;
    int pass = 0;

    if (out == NULL) return -1;
    memset(out, 0, sizeof *out);
    if (nv == 0 || !Rawtex_has_uv(uv) || !isfinite(du) || !isfinite(dv) ||
        du <= 0.0 || dv <= 0.0)
        return -1;
    for (i = 0; i < nv; i++) {
        double uu = (double)Rawtex_uv(uv,i * 2 + 0), vv = (double)Rawtex_uv(uv,i * 2 + 1);
        if (!isfinite(uu) || !isfinite(vv)) return -1;
        if (window != NULL) continue;
        if (i == 0) {
            out->umin = out->umax = uu;
            out->vmin = out->vmax = vv;
            continue;
        }
        if (uu < out->umin) out->umin = uu;
        if (uu > out->umax) out->umax = uu;
        if (vv < out->vmin) out->vmin = vv;
        if (vv > out->vmax) out->vmax = vv;
    }
    if (window != NULL) {
        if (!isfinite(window->umin) || !isfinite(window->umax) ||
            !isfinite(window->vmin) || !isfinite(window->vmax) ||
            window->umax <= window->umin || window->vmax <= window->vmin)
            return -1;
        out->umin = window->umin;
        out->umax = window->umax;
        out->vmin = window->vmin;
        out->vmax = window->vmax;
    }
    out->max_px = max_px > 0 ? max_px : RAWTEX_MAX_PX_DEFAULT;
    span_u = out->umax - out->umin;
    span_v = out->vmax - out->vmin;
    out->W = plan_dim(span_u, du);
    out->H = plan_dim(span_v, dv);
    out->ok = out->W <= RAWTEX_AXIS_CAP && out->H <= RAWTEX_AXIS_CAP &&
              out->W * out->H <= out->max_px;

    /* Smallest steps that fit: satisfy each axis cap independently, then
     * shrink only the LARGER axis to meet the total budget -- the old advice
     * scaled both axes by sqrt(area) and silently downsampled a full-height
     * 512-row v axis because u had exploded. */
    out->need_du = du;
    out->need_dv = dv;
    if (out->W > RAWTEX_AXIS_CAP)
        out->need_du = span_u / (double)RAWTEX_AXIS_CAP * 1.0000001;
    if (out->H > RAWTEX_AXIS_CAP)
        out->need_dv = span_v / (double)RAWTEX_AXIS_CAP * 1.0000001;
    for (pass = 0; pass < 8; pass++) {
        size_t w = plan_dim(span_u, out->need_du);
        size_t h = plan_dim(span_v, out->need_dv);
        double target = 0.0;
        if (w <= RAWTEX_AXIS_CAP && h <= RAWTEX_AXIS_CAP &&
            w * h <= out->max_px)
            break;
        if (w >= h) {
            target = (double)(out->max_px / (h > 0 ? h : 1));
            if (target < 1.0) target = 1.0;
            out->need_du = span_u / target * 1.001;
        } else {
            target = (double)(out->max_px / (w > 0 ? w : 1));
            if (target < 1.0) target = 1.0;
            out->need_dv = span_v / target * 1.001;
        }
    }
    return 0;
}

int Rawtex_plan(const float *uv, size_t nv, double du, double dv,
                size_t max_px, RawtexPlan *out)
{
    return Rawtex_plan_window(uv, nv, du, dv, max_px, NULL, out);
}

int Rawtex_write_tif_field(const char *path, CubeTable *ct,
                            const float *verts, RawtexUv uv, size_t nv,
                            const int32_t *faces, size_t nf,
                            const float *normals,
                            const uint8_t *face_skip,
                            double range, int nsteps,
                            double du, double dv, double lo, double hi,
                            double stretch_ratio, double stretch_floor,
                            double max_edge3d, size_t max_px,
                            const RawtexWindow *window,
                            const DiagOpts *diag,
                            size_t *out_W, size_t *out_H,
                            double *out_fill, size_t *out_multi,
                            size_t *out_skip_uv, size_t *out_skip_3d)
{
    /* A boundary-only hit is ownership, not a second surface layer.  The
     * fitted ribbon deliberately puts vertices on its raster lattice, so a
     * pixel can lie on a shared diagonal (or even on a grid vertex touched by
     * several triangles).  Counting every closed-triangle incidence made an
     * exactly injective ribbon look almost entirely multi-covered.  Keep one
     * boundary owner until a strict-interior hit appears; only multiple strict
     * interiors are evidence of an overlap.  The high bit avoids another
     * W*H owner raster.  Exact UV intersection auditing remains the authority
     * for measure-zero edge coincidences. */
    const uint32_t boundary_hit = UINT32_C(0x80000000);
    size_t skip_uv = 0, skip_3d = 0;
    double umin = 0.0, umax = 0.0, vmin = 0.0, vmax = 0.0;
    size_t W = 0, H = 0, i = 0, f = 0, px = 0, filled = 0, multi = 0;
    size_t centroid_stamps = 0;
    double *sum = NULL;
    uint32_t *cnt = NULL;
    uint8_t *img = NULL;
    /* diag buffers */
    uint8_t *dsmax = NULL, *dsmin = NULL, *dverr = NULL, *dcls = NULL;
    uint8_t *cover_uv = NULL, *cover_3d = NULL;
    uint8_t *is_smear = NULL;   /* per-face skip-uv flag for the 3D dump */
    int32_t *component_parent = NULL, *component_label = NULL;
    int32_t *owner_component = NULL;
    float *owner_position = NULL;
    RawtexPairs component_pairs={0};
    size_t component_count = 0;
    uint8_t *component_active = NULL;
    size_t *component_vertices = NULL;
    double *component_u_lo = NULL, *component_u_hi = NULL;
    uint8_t *multi_topology = NULL; /* bit 0: same CC, bit 1: cross CC */
    uint8_t *multi_physical = NULL; /* bit 0: <=6 vox, bit 1: >6 vox */
    double zoff = 0.0;
    int rc = 0, coverage_rc = 0;

    *out_W = 0; *out_H = 0; *out_fill = 0.0; *out_multi = 0;
    if (out_skip_uv != NULL) *out_skip_uv = 0;
    if (out_skip_3d != NULL) *out_skip_3d = 0;
    if (nv == 0 || nf == 0 || !Rawtex_has_uv(uv) || du <= 0.0 || dv <= 0.0) return -1;
    /* Size from the UV BOUNDING BOX.  A winding frame lifted from registration
     * has its natural origin wherever the spiral starts -- on the 4x21x21 that
     * is u in [-1174233, +341207] and v in [4352, 4864], v being world z -- so
     * anchoring the raster at (0,0) silently drops every negative-u vertex,
     * which there is 99.98% of the mesh, and wastes 4352 empty rows in v.
     * Sizing + caps live in Rawtex_plan so callers can preflight BEFORE the
     * expensive sampling pass; a plan that reports ok never rejects here. */
    {
        RawtexPlan plan;
        if (Rawtex_plan_field(uv, nv, du, dv, max_px, window, &plan) != 0)
            return -1;
        if (!plan.ok) {
            fprintf(stderr,
                "ERROR: raster %zux%zu unreasonable for uv u=[%.1f,%.1f] "
                "v=[%.1f,%.1f] at du=%.3g dv=%.3g; retry with --raster-du "
                "%.3g --raster-dv %.3g, or --raster-auto / --raster-max-px\n",
                plan.W, plan.H, plan.umin, plan.umax, plan.vmin, plan.vmax,
                du, dv, plan.need_du, plan.need_dv);
            return -1;
        }
        umin = plan.umin; umax = plan.umax;
        vmin = plan.vmin; vmax = plan.vmax;
        W = plan.W; H = plan.H;
    }
    sum = (double *)xcalloc(W * H, sizeof(double));
    cnt = (uint32_t *)xcalloc(W * H, sizeof(uint32_t));
    img = (uint8_t *)xcalloc(W * H, sizeof(uint8_t));
    if (diag != NULL) {
        size_t n_probe = 0, stride = 0;
        double *zv = NULL;
        dsmax = (uint8_t *)xcalloc(W * H, 1);
        dsmin = (uint8_t *)xcalloc(W * H, 1);
        dverr = (uint8_t *)xcalloc(W * H, 1);
        dcls = (uint8_t *)xcalloc(W * H, 1);
        cover_uv = (uint8_t *)xcalloc(W * H, 1);
        cover_3d = (uint8_t *)xcalloc(W * H, 1);
        component_parent = (int32_t *)xmalloc(nv * sizeof(*component_parent));
        component_label = (int32_t *)xmalloc(nv * sizeof(*component_label));
        owner_component = (int32_t *)xmalloc(W * H * sizeof(*owner_component));
        owner_position = (float *)xmalloc(W * H * 3 * sizeof(*owner_position));
        multi_topology = (uint8_t *)xcalloc(W * H, 1);
        multi_physical = (uint8_t *)xcalloc(W * H, 1);
        /* INT32_MIN marks a vertex that is not referenced by any face.  It is
         * not a singleton mesh component: compact diagnostic/probe meshes can
         * deliberately retain unused input vertices, and counting those as
         * components makes the dense conflict matrix grow quadratically. */
        for (i = 0; i < nv; i++) component_parent[i] = INT32_MIN;
        for (f = 0; f < nf; f++) {
            int32_t a = faces[f * 3];
            int32_t b = faces[f * 3 + 1];
            int32_t c = faces[f * 3 + 2];
            if (component_parent[a] == INT32_MIN) component_parent[a] = -1;
            if (component_parent[b] == INT32_MIN) component_parent[b] = -1;
            if (component_parent[c] == INT32_MIN) component_parent[c] = -1;
            rawtex_component_union(component_parent, faces[f * 3],
                                    faces[f * 3 + 1]);
            rawtex_component_union(component_parent, faces[f * 3],
                                    faces[f * 3 + 2]);
        }
        for (i = 0; i < nv; i++) component_label[i] = -1;
        for (i = 0; i < nv; i++)
            if (component_parent[i] < 0 && component_parent[i] != INT32_MIN)
                component_label[i] = (int32_t)component_count++;
        component_pairs.byte_limit=RAWTEX_CONFLICT_HASH_BYTES;
        component_active=(uint8_t *)xcalloc(component_count,1);
        component_vertices = (size_t *)xcalloc(
            component_count, sizeof(*component_vertices));
        component_u_lo = (double *)xmalloc(
            component_count * sizeof(*component_u_lo));
        component_u_hi = (double *)xmalloc(
            component_count * sizeof(*component_u_hi));
        for (i = 0; i < component_count; i++) {
            component_u_lo[i] = 1e300;
            component_u_hi[i] = -1e300;
        }
        for (i = 0; i < nv; i++) {
            int32_t label;
            double u = (double)Rawtex_uv(uv,i * 2);
            if (component_parent[i] == INT32_MIN) continue;
            label = component_label[
                rawtex_component_find(component_parent, (int32_t)i)];
            component_vertices[label]++;
            if (u < component_u_lo[label]) component_u_lo[label] = u;
            if (u > component_u_hi[label]) component_u_hi[label] = u;
        }
        for (px = 0; px < W * H; px++) owner_component[px] = -1;
        /* zmin estimate: median of (z - vt.v) over a vertex sample -- exact
         * on correctly-mapped verts, robust to the broken ones */
        stride = (nv > 200000) ? nv / 200000 : 1;
        zv = (double *)xmalloc((nv / stride + 1) * sizeof(double));
        for (i = 0; i < nv; i += stride)
            zv[n_probe++] = (double)verts[i * 3 + 0] - (double)Rawtex_uv(uv,i * 2 + 1);
        qsort(zv, n_probe, sizeof(double), cmp_double);
        zoff = zv[n_probe / 2];
        free(zv);
        if (diag->smear_obj != NULL)
            is_smear = (uint8_t *)xcalloc(nf, 1);
    }

    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3 + 0];
        size_t b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        int32_t face_component = component_parent != NULL
                               ? component_label[rawtex_component_find(
                                     component_parent, (int32_t)a)]
                               : -1;
        double ua = ((double)Rawtex_uv(uv,a * 2 + 0) - umin) / du;
        double va = ((double)Rawtex_uv(uv,a * 2 + 1) - vmin) / dv;
        double ub = ((double)Rawtex_uv(uv,b * 2 + 0) - umin) / du;
        double vb = ((double)Rawtex_uv(uv,b * 2 + 1) - vmin) / dv;
        double uc = ((double)Rawtex_uv(uv,c * 2 + 0) - umin) / du;
        double vc = ((double)Rawtex_uv(uv,c * 2 + 1) - vmin) / dv;
        double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        double lox = ua < ub ? (ua < uc ? ua : uc) : (ub < uc ? ub : uc);
        double hix = ua > ub ? (ua > uc ? ua : uc) : (ub > uc ? ub : uc);
        double loy = va < vb ? (va < vc ? va : vc) : (vb < vc ? vb : vc);
        double hiy = va > vb ? (va > vc ? va : vc) : (vb > vc ? vb : vc);
        long x0 = 0, x1 = 0, y0 = 0, y1 = 0, xx = 0, yy = 0;
        int bad_uv = 0, bad_3d = 0;
        uint8_t f_smax = 0, f_smin = 0, f_verr = 0;
        if (face_skip != NULL && face_skip[f]) continue;   /* omitted entirely */
        if (fabs(A2) < 1e-12) continue;   /* degenerate in UV */
        {   /* face gates: UV stretch + max 3D edge */
            size_t vi[3];
            int e = 0;
            vi[0] = a; vi[1] = b; vi[2] = c;
            for (e = 0; e < 3; e++) {
                size_t p = vi[e], q = vi[(e + 1) % 3];
                double d3 = 0.0, dueg = 0.0, dveg = 0.0, e2d = 0.0;
                int k = 0;
                for (k = 0; k < 3; k++) {
                    double dd = (double)verts[q * 3 + k]
                                - (double)verts[p * 3 + k];
                    d3 += dd * dd;
                }
                d3 = sqrt(d3);
                dueg = (double)Rawtex_uv(uv,q * 2 + 0) - (double)Rawtex_uv(uv,p * 2 + 0);
                dveg = (double)Rawtex_uv(uv,q * 2 + 1) - (double)Rawtex_uv(uv,p * 2 + 1);
                e2d = sqrt(dueg * dueg + dveg * dveg);
                if (max_edge3d > 0.0 && d3 > max_edge3d) bad_3d = 1;
                if (stretch_ratio > 0.0
                    && e2d > (stretch_ratio * d3 > stretch_floor
                              ? stretch_ratio * d3 : stretch_floor)) bad_uv = 1;
            }
            if (bad_3d) skip_3d++;
            else if (bad_uv) skip_uv++;
            if (is_smear != NULL && bad_uv && !bad_3d) is_smear[f] = 1;
            if ((bad_3d || bad_uv) && diag == NULL) continue;
        }
        if (diag != NULL) {
            double smin = 0.0, smax = face_sigma(verts, uv, a, b, c, &smin);
            double e0 = fabs(((double)verts[a * 3 + 0] - zoff)
                             - (double)Rawtex_uv(uv,a * 2 + 1));
            double e1 = fabs(((double)verts[b * 3 + 0] - zoff)
                             - (double)Rawtex_uv(uv,b * 2 + 1));
            double e2 = fabs(((double)verts[c * 3 + 0] - zoff)
                             - (double)Rawtex_uv(uv,c * 2 + 1));
            double emax = e0 > e1 ? (e0 > e2 ? e0 : e2) : (e1 > e2 ? e1 : e2);
            f_smax = log2_to_u8(smax);
            f_smin = log2_to_u8(smin);
            emax *= 32.0;
            f_verr = (uint8_t)(emax > 255.0 ? 255.0 : emax);
        }
        /* Centers are sampled at xx+0.5, so the first covered pixel is
         * ceil(lo - 0.5).  The former ceil(lo - 0.001) bound implicitly
         * assumed integer-lattice UVs (every historical quadribbon) and
         * silently dropped each face's first row/column at generic UV
         * phases -- the metric-projection bake lost ~20% of covered pixels
         * to that in smoothly swirling moire bands.  The upper bound only
         * overscans; the barycentric test rejects those centers. */
        x0 = (long)ceil(lox - 0.501); x1 = (long)floor(hix + 0.001);
        y0 = (long)ceil(loy - 0.501); y1 = (long)floor(hiy + 0.001);
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= (long)W) x1 = (long)W - 1;
        if (y1 >= (long)H) y1 = (long)H - 1;
        int face_hit = 0;
        for (yy = y0; yy <= y1; yy++) {
            for (xx = x0; xx <= x1; xx++) {
                /* half-texel: sample the pixel's CENTER, never a lattice
                 * line (see plan_dim) */
                double pu = (double)xx + 0.5, pv = (double)yy + 0.5;
                double l1 = ((pu - ua) * (vc - va) - (pv - va) * (uc - ua)) / A2;
                double l2 = ((ub - ua) * (pv - va) - (vb - va) * (pu - ua)) / A2;
                double l0 = 1.0 - l1 - l2;
                size_t pi = 0;
                float p3[3], n3[3];
                double nn = 0.0, s = 0.0;
                int k = 0;
                int on_boundary = 0;
                if (l0 < -1e-9 || l1 < -1e-9 || l2 < -1e-9) continue;
                face_hit = 1;
                on_boundary = l0 <= 1e-9 || l1 <= 1e-9 || l2 <= 1e-9;
                pi = (size_t)yy * W + (size_t)xx;
                if (diag != NULL) {
                    if (f_smax > dsmax[pi]) dsmax[pi] = f_smax;
                    if (dsmin[pi] == 0 || (f_smin != 0 && f_smin < dsmin[pi]))
                        dsmin[pi] = f_smin;
                    if (f_verr > dverr[pi]) dverr[pi] = f_verr;
                    if (bad_uv) cover_uv[pi] = 1;
                    else if (bad_3d) cover_3d[pi] = 1;
                }
                if (bad_uv || bad_3d) continue;   /* diagnosed, not painted */
                for (k = 0; k < 3; k++) {
                    p3[k] = (float)(l0 * (double)verts[a * 3 + k]
                                    + l1 * (double)verts[b * 3 + k]
                                    + l2 * (double)verts[c * 3 + k]);
                    n3[k] = (normals != NULL)
                            ? (float)(l0 * (double)normals[a * 3 + k]
                                      + l1 * (double)normals[b * 3 + k]
                                      + l2 * (double)normals[c * 3 + k])
                            : 0.0f;
                }
                nn = sqrt((double)n3[0] * n3[0] + (double)n3[1] * n3[1]
                          + (double)n3[2] * n3[2]);
                if (nn > 1e-6) {
                    n3[0] = (float)((double)n3[0] / nn);
                    n3[1] = (float)((double)n3[1] / nn);
                    n3[2] = (float)((double)n3[2] / nn);
                }
                s = sample_vertex(ct, p3, nn > 0.5 ? n3 : NULL, range, nsteps);
                if (s < 0.0) continue;
                if (on_boundary) {
                    /* Shared edges/vertices have one deterministic first
                     * owner.  A strict interior already stored at this pixel
                     * is stronger evidence and is left untouched. */
                    if (cnt[pi] == 0) {
                        sum[pi] = s;
                        cnt[pi] = boundary_hit | UINT32_C(1);
                    }
                } else if ((cnt[pi] & boundary_hit) != 0) {
                    /* Replace a boundary fallback by the unique strict
                     * interior sample. */
                    sum[pi] = s;
                    cnt[pi] = 1;
                    if (owner_component != NULL)
                        owner_component[pi] = face_component;
                    if (owner_position != NULL)
                        memcpy(&owner_position[pi * 3], p3, 3 * sizeof(*p3));
                } else {
                    if (owner_component != NULL && cnt[pi] != 0) {
                        if (owner_component[pi] == face_component)
                            multi_topology[pi] |= 1;
                        else {
                            double dz = (double)p3[0] -
                                        (double)owner_position[pi * 3];
                            double dy = (double)p3[1] -
                                        (double)owner_position[pi * 3 + 1];
                            double dx = (double)p3[2] -
                                        (double)owner_position[pi * 3 + 2];
                            multi_topology[pi] |= 2;
                            if (dz * dz + dy * dy + dx * dx <= 36.0)
                                multi_physical[pi] |= 1;
                            else
                                multi_physical[pi] |= 2;
                            if (rawtex_pair_add(&component_pairs,owner_component[pi],face_component)!=0) {
                                fprintf(stderr,"ERROR: exact component conflict counts exceed the sparse storage/counter budget; bake not published\n");
                                rc=-1; goto done;
                            }
                            component_active[owner_component[pi]]=1;
                            component_active[face_component]=1;
                        }
                    }
                    sum[pi] += s;
                    cnt[pi]++;
                    if (owner_component != NULL && cnt[pi] == 1)
                        owner_component[pi] = face_component;
                    if (owner_position != NULL && cnt[pi] == 1)
                        memcpy(&owner_position[pi * 3], p3, 3 * sizeof(*p3));
                }
            }
        }
        if (!face_hit && !bad_uv && !bad_3d) {
            /* Sub-pixel face: its bbox caught no pixel CENTER (coarse du/dv
             * can make a roughly 2-vox face a half-pixel footprint).
             * Stamp the face centroid as a weak boundary-class sample so
             * coverage survives any raster resolution; a strict interior
             * owner still wins the pixel. */
            double cu = (ua + ub + uc) / 3.0, cv = (va + vb + vc) / 3.0;
            long sx2 = (long)cu, sy2 = (long)cv;
            if (sx2 >= 0 && sy2 >= 0 && sx2 < (long)W && sy2 < (long)H) {
                size_t pi = (size_t)sy2 * W + (size_t)sx2;
                if (cnt[pi] == 0) {
                    float p3[3], n3[3];
                    double nn = 0.0, s = 0.0;
                    int k = 0;
                    for (k = 0; k < 3; k++) {
                        p3[k] = (float)(((double)verts[a * 3 + k]
                                         + (double)verts[b * 3 + k]
                                         + (double)verts[c * 3 + k]) / 3.0);
                        n3[k] = (normals != NULL)
                                ? (float)(((double)normals[a * 3 + k]
                                           + (double)normals[b * 3 + k]
                                           + (double)normals[c * 3 + k]) / 3.0)
                                : 0.0f;
                    }
                    nn = sqrt((double)n3[0] * n3[0] + (double)n3[1] * n3[1]
                              + (double)n3[2] * n3[2]);
                    if (nn > 1e-6) {
                        n3[0] = (float)((double)n3[0] / nn);
                        n3[1] = (float)((double)n3[1] / nn);
                        n3[2] = (float)((double)n3[2] / nn);
                    }
                    s = sample_vertex(ct, p3, nn > 0.5 ? n3 : NULL,
                                      range, nsteps);
                    if (s >= 0.0) {
                        sum[pi] = s;
                        cnt[pi] = boundary_hit | UINT32_C(1);
                        centroid_stamps++;
                    }
                }
            }
        }
    }
    for (px = 0; px < W * H; px++) {
        uint32_t ncover = cnt[px] & ~boundary_hit;
        if (ncover > 0) {
            double g = gray_of(sum[px] / (double)ncover, lo, hi);
            img[px] = (uint8_t)(g * 255.0 + 0.5);
            filled++;
            if (ncover > 1) multi++;
        }
    }
    if (centroid_stamps > 0)
        fprintf(stderr, "  sub-pixel faces stamped at centroid: %zu\n",
                centroid_stamps);
    if (multi_topology != NULL) {
        size_t same = 0, cross = 0, mixed = 0;
        for (px = 0; px < W * H; px++) {
            if (multi_topology[px] == 1) same++;
            else if (multi_topology[px] == 2) cross++;
            else if (multi_topology[px] == 3) mixed++;
        }
        fprintf(stderr,
                "  multi-cover topology: same-component=%zu "
                "cross-component=%zu mixed=%zu\n",
                same, cross, mixed);
        if (multi_physical != NULL) {
            size_t near = 0, far = 0, both = 0;
            for (px = 0; px < W * H; px++) {
                if (multi_physical[px] == 1) near++;
                else if (multi_physical[px] == 2) far++;
                else if (multi_physical[px] == 3) both++;
            }
            fprintf(stderr,
                    "  cross-component physical separation: near=%zu "
                    "far=%zu mixed=%zu (gate %.1f vox)\n",
                    near, far, both, 6.0);
        }
        if (component_active != NULL) {
            rawtex_pairs_report(&component_pairs,component_count);
            for (size_t a = 0; a < component_count; a++) {
                if (component_active[a])
                    fprintf(stderr,
                            "    component %zu vertices=%zu U=[%.1f,%.1f] "
                            "width=%.1f\n",
                            a, component_vertices[a], component_u_lo[a],
                            component_u_hi[a],
                            component_u_hi[a] - component_u_lo[a]);
            }
        }
    }
    /* Bounded void-fill: an uncovered pixel with >=5 of 8 supported
     * neighbours takes their mean.  Two passes close the 1-2 px cracks and
     * pinholes that a dropped grid cell leaves, without inventing texture at
     * real-hole scale.  The coverage sibling keeps display honest: painted
     * (255) vs void-filled (128) vs background (0); out_fill counts only
     * genuinely painted pixels. */
    {
        uint8_t *cover = (uint8_t *)xcalloc(W * H, 1);
        size_t vf_total = 0;
        int pass2 = 0;
        for (px = 0; px < W * H; px++)
            if ((cnt[px] & ~boundary_hit) > 0) cover[px] = 1;
        for (pass2 = 0; pass2 < 2; pass2++) {
            size_t cap_fill = 4096, nfill = 0;
            size_t *fill_px = (size_t *)xmalloc(cap_fill * sizeof(size_t));
            uint8_t *fill_val = (uint8_t *)xmalloc(cap_fill);
            size_t yy2 = 0, xx2 = 0;
            for (yy2 = 0; yy2 < H; yy2++) {
                for (xx2 = 0; xx2 < W; xx2++) {
                    size_t p2 = yy2 * W + xx2;
                    int nsupp = 0;
                    unsigned int acc = 0;
                    int dy2 = 0, dx2 = 0;
                    int up = 0, down = 0, left = 0, right = 0;
                    if (cover[p2] != 0) continue;
                    for (dy2 = -1; dy2 <= 1; dy2++) {
                        for (dx2 = -1; dx2 <= 1; dx2++) {
                            long ny2 = (long)yy2 + dy2;
                            long nx2 = (long)xx2 + dx2;
                            size_t np2 = 0;
                            if ((dy2 == 0 && dx2 == 0) || ny2 < 0 ||
                                nx2 < 0 || ny2 >= (long)H || nx2 >= (long)W)
                                continue;
                            np2 = (size_t)ny2 * W + (size_t)nx2;
                            if (cover[np2] != 0) {
                                nsupp++;
                                acc += img[np2];
                                if (dy2 < 0 && dx2 == 0) up = 1;
                                if (dy2 > 0 && dx2 == 0) down = 1;
                                if (dy2 == 0 && dx2 < 0) left = 1;
                                if (dy2 == 0 && dx2 > 0) right = 1;
                            }
                        }
                    }
                    /* A 1-2 px slit with paint on BOTH sides is always
                     * artificial: adjacent wraps sit a full pitch (~9.5 vox =
                     * 9+ px at dv=1) apart in v and real u-holes split into
                     * gutter runs.  One dropped GRID cell is TWO raster rows
                     * (grid_dv 2 vox, raster dv 1), so the reach is 2. */
                    if (!up && yy2 >= 2)
                        up = cover[p2 - 2 * W] != 0;
                    if (!down && yy2 + 2 < H)
                        down = cover[p2 + 2 * W] != 0;
                    if (!left && xx2 >= 2)
                        left = cover[p2 - 2] != 0;
                    if (!right && xx2 + 2 < W)
                        right = cover[p2 + 2] != 0;
                    if (nsupp < 5 && !(up && down) && !(left && right))
                        continue;
                    if (nsupp == 0) continue;
                    if (nfill == cap_fill) {
                        size_t ncap = cap_fill * 2;
                        size_t *npx = (size_t *)xmalloc(ncap * sizeof(size_t));
                        uint8_t *nvl = (uint8_t *)xmalloc(ncap);
                        memcpy(npx, fill_px, nfill * sizeof(size_t));
                        memcpy(nvl, fill_val, nfill);
                        free(fill_px); free(fill_val);
                        fill_px = npx; fill_val = nvl;
                        cap_fill = ncap;
                    }
                    fill_px[nfill] = p2;
                    fill_val[nfill] = (uint8_t)(acc / (unsigned int)nsupp);
                    nfill++;
                }
            }
            for (px = 0; px < nfill; px++) {
                img[fill_px[px]] = fill_val[px];
                cover[fill_px[px]] = 2;
            }
            vf_total += nfill;
            free(fill_px);
            free(fill_val);
            if (nfill == 0) break;
        }
        if (vf_total > 0)
            fprintf(stderr, "  void-filled %zu crack/pinhole px "
                    "(coverage mask tags them 128)\n", vf_total);
        {   /* coverage sibling: <base>_coverage.png */
            char cpath[2600];
            size_t plen = strlen(path);
            uint8_t *cimg = (uint8_t *)xmalloc(W * H);
            for (px = 0; px < W * H; px++)
                cimg[px] = cover[px] == 1 ? 255 : (cover[px] == 2 ? 128 : 0);
            snprintf(cpath, sizeof cpath, "%s", path);
            if (plen > 4 && strcmp(cpath + plen - 4, ".tif") == 0)
                snprintf(cpath + plen - 4, sizeof cpath - (plen - 4),
                         "_coverage.png");
            else
                snprintf(cpath + plen, sizeof cpath - plen, "_coverage.png");
            if (VesPng_write_gray(cpath, cimg, (int)W, (int)H) != 0)
                fprintf(stderr, "  coverage PNG unavailable; retaining TIFF and tiles\n");
            /* TIF sibling: the C sheet compositor reads coverage via TiffIO
             * (ves_png is write-only) */
            snprintf(cpath, sizeof cpath, "%s", path);
            if (plen > 4 && strcmp(cpath + plen - 4, ".tif") == 0)
                snprintf(cpath + plen - 4, sizeof cpath - (plen - 4),
                         "_coverage.tif");
            else
                snprintf(cpath + plen, sizeof cpath - plen, "_coverage.tif");
            coverage_rc = TiffIO_save(cpath, cimg, 1, (int)H, (int)W);
            free(cimg);
        }
        free(cover);
    }
    rc = TiffIO_save(path, img, 1, (int)H, (int)W);
    if (coverage_rc != 0) rc = -1;
    {   /* also a full-res grayscale PNG sibling (<path>.png), written from C
         * via stb -- no GDI+ round-trip (which truncates >~3k-col images). */
        char pngpath[2600];
        size_t plen = strlen(path);
        snprintf(pngpath, sizeof pngpath, "%s", path);
        if (plen > 4 && strcmp(pngpath + plen - 4, ".tif") == 0)
            snprintf(pngpath + plen - 4, sizeof pngpath - (plen - 4), ".png");
        else
            snprintf(pngpath + plen, sizeof pngpath - plen, ".png");
        if (VesPng_write_gray(pngpath, img, (int)W, (int)H) != 0)
            fprintf(stderr, "  texture PNG unavailable; retaining TIFF and tiles\n");
    }
    {   /* multi-line strip sibling (<base>_strip.tif + .png): a readable
         * stacked view of the very wide unroll (see tif_strip.h). Every bake
         * stage thus delivers both <stage>_rawtex.tif and _rawtex_strip.tif. */
        char sp[2600];
        size_t plen = strlen(path);
        snprintf(sp, sizeof sp, "%s", path);
        if (plen > 4 && strcmp(sp + plen - 4, ".tif") == 0)
            snprintf(sp + plen - 4, sizeof sp - (plen - 4), "_strip.tif");
        else
            snprintf(sp + plen, sizeof sp - plen, "_strip.tif");
        TifStrip_write_gray(sp, img, (int)W, (int)H, 0, -1);
    }
    *out_W = W; *out_H = H;
    *out_fill = (double)filled / (double)(W * H);
    *out_multi = multi;
    if (out_skip_uv != NULL) *out_skip_uv = skip_uv;
    if (out_skip_3d != NULL) *out_skip_3d = skip_3d;

    if (diag != NULL && rc == 0) {
        char dpath[2600];
        size_t n_bg = 0, n_ok = 0, n_suv = 0, n_s3d = 0, n_multi = 0, n_dark = 0;
        for (px = 0; px < W * H; px++) {
            uint8_t cl = 0;
            if (cover_uv[px]) cl = 2;
            else if (cover_3d[px]) cl = 3;
            else if ((cnt[px] & ~boundary_hit) > 1) cl = 4;
            else if ((cnt[px] & ~boundary_hit) == 1)
                cl = (img[px] < diag->dark_thresh) ? 5 : 1;
            dcls[px] = cl;
            switch (cl) {
            case 0: n_bg++; break;
            case 1: n_ok++; break;
            case 2: n_suv++; break;
            case 3: n_s3d++; break;
            case 4: n_multi++; break;
            default: n_dark++; break;
            }
        }
        fprintf(stderr, "  [diag] class px: ok=%zu (%.2f%%) skip-uv=%zu "
                "(%.2f%%) skip-3d=%zu (%.2f%%) dark=%zu (%.2f%%) multi=%zu "
                "bg=%zu  (%% of non-bg non-multi)\n",
                n_ok, 100.0 * (double)n_ok
                      / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_suv, 100.0 * (double)n_suv
                       / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_s3d, 100.0 * (double)n_s3d
                       / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_dark, 100.0 * (double)n_dark
                        / (double)(n_ok + n_suv + n_s3d + n_dark + 1),
                n_multi, n_bg);
        {   /* worst-tile table: defect fraction per 256x64-px tile */
            size_t tw = 256, th = 64;
            size_t ntx = (W + tw - 1) / tw, nty = (H + th - 1) / th;
            typedef struct { double bad; size_t tx, ty, nbad, ncov; } Tile;
            Tile *tiles = (Tile *)xcalloc(ntx * nty, sizeof(Tile));
            size_t ti = 0, shown = 0;
            for (px = 0; px < W * H; px++) {
                size_t x = px % W, y = px / W;
                uint8_t cl = dcls[px];
                Tile *T = &tiles[(y / th) * ntx + (x / tw)];
                T->tx = x / tw; T->ty = y / th;
                if (cl == 1) T->ncov++;
                else if (cl == 2 || cl == 3 || cl == 5) { T->ncov++; T->nbad++; }
            }
            for (ti = 0; ti < ntx * nty; ti++)
                tiles[ti].bad = (tiles[ti].ncov > 300)
                    ? (double)tiles[ti].nbad / (double)tiles[ti].ncov : 0.0;
            /* selection sort of top 12 (tiny n) */
            for (shown = 0; shown < 12 && shown < ntx * nty; shown++) {
                size_t best = shown;
                Tile tmp;
                for (ti = shown; ti < ntx * nty; ti++)
                    if (tiles[ti].bad > tiles[best].bad) best = ti;
                tmp = tiles[shown]; tiles[shown] = tiles[best];
                tiles[best] = tmp;
                if (tiles[shown].bad <= 0.0) break;
                fprintf(stderr, "  [diag] hot tile %2zu: px(%zu..%zu, %zu..%zu)"
                        " u=[%.0f,%.0f] v=[%.0f,%.0f]  bad=%.1f%% (%zu/%zu)\n",
                        shown,
                        tiles[shown].tx * tw, (tiles[shown].tx + 1) * tw,
                        tiles[shown].ty * th, (tiles[shown].ty + 1) * th,
                         umin + (double)(tiles[shown].tx * tw) * du,
                         umin + (double)((tiles[shown].tx + 1) * tw) * du,
                         vmin + (double)(tiles[shown].ty * th) * dv,
                         vmin + (double)((tiles[shown].ty + 1) * th) * dv,
                        100.0 * tiles[shown].bad,
                        tiles[shown].nbad, tiles[shown].ncov);
            }
            free(tiles);
        }
        snprintf(dpath, sizeof dpath, "%s_diagclass.tif", diag->prefix);
        TiffIO_save(dpath, dcls, 1, (int)H, (int)W);
        snprintf(dpath, sizeof dpath, "%s_diagstretch.tif", diag->prefix);
        TiffIO_save(dpath, dsmax, 1, (int)H, (int)W);
        snprintf(dpath, sizeof dpath, "%s_diagsquash.tif", diag->prefix);
        TiffIO_save(dpath, dsmin, 1, (int)H, (int)W);
        snprintf(dpath, sizeof dpath, "%s_diagverr.tif", diag->prefix);
        TiffIO_save(dpath, dverr, 1, (int)H, (int)W);
        fprintf(stderr, "  [diag] wrote %s_diag{class,stretch,squash,verr}.tif"
                " (zmin est %.2f)\n", diag->prefix, zoff);
    }
    if (is_smear != NULL) {
        /* 3D dump of the skip-uv (smear) faces + a z-slice histogram, to
         * localize where the u-stretch defects sit on the scroll. */
        FILE *fp = fopen(diag->smear_obj, "wb");
        size_t nsm = 0;
        long zhist[16];
        double zmin = 1e300, zmax = -1e300;
        int hb = 0;
        memset(zhist, 0, sizeof zhist);
        for (i = 0; i < nv; i++) {
            double z = (double)verts[i * 3 + 0];
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        /* winding-jump test: for each smear face, is max|du| across its edges
         * close to an integer multiple of 2*pi*r (a winding-index disagreement)
         * or arbitrary (a QP/registration offset)? r = mean radial distance of
         * the face from the axis (axis_point z,y,x; scroll axis = Z). */
        long ratio_hist[8];   /* du/(2*pi*r) rounded: 0,1,2,3,4,5,6,>=7 */
        memset(ratio_hist, 0, sizeof ratio_hist);
        if (fp != NULL) {
            for (i = 0; i < nv; i++)
                fprintf(fp, "v %.4f %.4f %.4f\n", (double)verts[i * 3 + 0],
                        (double)verts[i * 3 + 1], (double)verts[i * 3 + 2]);
            for (f = 0; f < nf; f++) {
                double zc = 0.0, umn = 1e300, umx = -1e300, rsum = 0.0, twopir = 0.0;
                int bidx = 0, m = 0, rr = 0;
                if (!is_smear[f]) continue;
                fprintf(fp, "f %d %d %d\n", faces[f * 3 + 0] + 1,
                        faces[f * 3 + 1] + 1, faces[f * 3 + 2] + 1);
                nsm++;
                zc = ((double)verts[(size_t)faces[f*3+0]*3]
                      + (double)verts[(size_t)faces[f*3+1]*3]
                      + (double)verts[(size_t)faces[f*3+2]*3]) / 3.0;
                bidx = (zmax > zmin)
                     ? (int)(15.0 * (zc - zmin) / (zmax - zmin)) : 0;
                if (bidx < 0) bidx = 0;
                if (bidx > 15) bidx = 15;
                zhist[bidx]++;
                for (m = 0; m < 3; m++) {
                    size_t vi = (size_t)faces[f*3+m];
                    double uu = (double)Rawtex_uv(uv,vi*2+0);
                    double dy = (double)verts[vi*3+1] - 3405.0;
                    double dx = (double)verts[vi*3+2] - 2878.0;
                    if (uu < umn) umn = uu;
                    if (uu > umx) umx = uu;
                    rsum += sqrt(dy*dy + dx*dx);
                }
                twopir = 6.283185307 * (rsum / 3.0);
                rr = (twopir > 1.0) ? (int)((umx - umn) / twopir + 0.5) : 7;
                if (rr < 0) rr = 0;
                if (rr > 7) rr = 7;
                ratio_hist[rr]++;
            }
            fclose(fp);
            fprintf(stderr, "  [diag] wrote %s (%zu smear faces, 3D)\n",
                    diag->smear_obj, nsm);
            fprintf(stderr, "  [diag] smear z-histogram over [%.0f,%.0f] "
                    "(16 bands):\n    ", zmin, zmax);
            for (hb = 0; hb < 16; hb++)
                fprintf(stderr, "%ld ", zhist[hb]);
            fprintf(stderr, "\n  [diag] smear du/(2pi r) rounded 0..6,>=7: ");
            for (hb = 0; hb < 8; hb++)
                fprintf(stderr, "%ld ", ratio_hist[hb]);
            fprintf(stderr, "  (>=1 = winding-index jump; 0 = sub-wrap "
                    "registration/QP)\n");
        }
    }
done:
    free(is_smear);
    free(sum); free(cnt); free(img);
    free(dsmax); free(dsmin); free(dverr); free(dcls);
    free(cover_uv); free(cover_3d);
    free(component_parent); free(component_label); free(owner_component);
    free(owner_position);
    if (component_pairs.arena) Arena_dispose(&component_pairs.arena);
    free(component_active); free(multi_topology);
    free(multi_physical);
    free(component_vertices); free(component_u_lo); free(component_u_hi);
    return rc;
}

int Rawtex_write_tif(const char *path, CubeTable *ct,
                     const float *verts, const float *uv, size_t nv,
                     const int32_t *faces, size_t nf,
                     const float *normals,
                     const uint8_t *face_skip,
                     double range, int nsteps,
                     double du, double dv, double lo, double hi,
                     double stretch_ratio, double stretch_floor,
                     double max_edge3d, size_t max_px,
                     const DiagOpts *diag,
                     size_t *out_W, size_t *out_H,
                     double *out_fill, size_t *out_multi,
                     size_t *out_skip_uv, size_t *out_skip_3d)
{
    return Rawtex_write_tif_window(
        path, ct, verts, uv, nv, faces, nf, normals, face_skip,
        range, nsteps, du, dv, lo, hi, stretch_ratio, stretch_floor,
        max_edge3d, max_px, NULL, diag, out_W, out_H, out_fill, out_multi,
        out_skip_uv, out_skip_3d);
}

int Rawtex_plan_window(const float *uv, size_t nv, double du, double dv,
                       size_t max_px, const RawtexWindow *window,
                       RawtexPlan *out)
{
    return Rawtex_plan_field((RawtexUv){uv,NULL},nv,du,dv,max_px,window,out);
}

int Rawtex_write_tif_window(const char *path, CubeTable *ct,
                            const float *verts, const float *uv, size_t nv,
                            const int32_t *faces, size_t nf,
                            const float *normals,
                            const uint8_t *face_skip,
                            double range, int nsteps,
                            double du, double dv, double lo, double hi,
                            double stretch_ratio, double stretch_floor,
                            double max_edge3d, size_t max_px,
                            const RawtexWindow *window,
                            const DiagOpts *diag,
                            size_t *out_W, size_t *out_H,
                            double *out_fill, size_t *out_multi,
                            size_t *out_skip_uv, size_t *out_skip_3d)
{
    return Rawtex_write_tif_field(path,ct,verts,(RawtexUv){uv,NULL},nv,faces,nf,normals,face_skip,
        range,nsteps,du,dv,lo,hi,stretch_ratio,stretch_floor,max_edge3d,max_px,window,diag,
        out_W,out_H,out_fill,out_multi,out_skip_uv,out_skip_3d);
}

int Rawtex_selftest(void)
{
    enum { N=16 };
    size_t dense[N*N]={0},nonzero=0,total=0;
    RawtexPairs pairs={0},limited={0};
    uint32_t state=813735u;
    int fail=0;
    pairs.byte_limit=1024u*1024u;
    for (size_t i=0;i<10000;i++) {
        int32_t a=0,b=0;
        state=1664525u*state+1013904223u; a=(int32_t)((state>>16)%N);
        state=1664525u*state+1013904223u; b=(int32_t)((state>>16)%N);
        if (a==b) continue;
        dense[(size_t)a*N+b]++; total++;
        fail+=rawtex_pair_add(&pairs,a,b)!=0;
    }
    for (size_t a=0;a<N;a++) for (size_t b=a+1;b<N;b++)
        nonzero+=dense[a*N+b]+dense[b*N+a]>0;
    fail+=pairs.used!=nonzero || pairs.peak_payload>pairs.byte_limit;
    for (size_t i=0;i<pairs.capacity;i++) if (pairs.table[i].key) {
        size_t a=(size_t)(pairs.table[i].key>>32)-1;
        size_t b=(size_t)(uint32_t)pairs.table[i].key-1;
        fail+=a>=N || b>=N;
        if (a<N && b<N) fail+=pairs.table[i].hits!=dense[a*N+b]+dense[b*N+a];
        total-=pairs.table[i].hits;
    }
    fail+=total!=0;
    fail+=rawtex_pair_add(&pairs,3984599,INT32_MAX)!=0;
    fail+=rawtex_pair_add(&pairs,INT32_MAX,3984599)!=0;
    fail+=rawtex_pair_add(&pairs,-1,0)==0 || rawtex_pair_add(&pairs,0,0)==0;
    limited.byte_limit=64*sizeof(RawtexPair);
    for (int32_t i=1;i<=32;i++) fail+=rawtex_pair_add(&limited,0,i)!=0;
    fail+=rawtex_pair_add(&limited,0,33)==0 || limited.used!=32;
    /* Existing counters stay usable at capacity; no accepted prefix is
     * published by the caller if a NEW pair cannot be counted. */
    fail+=rawtex_pair_add(&limited,1,0)!=0 || limited.used!=32;
    if (limited.capacity) {
        size_t slot=rawtex_pair_slot((UINT64_C(1)<<32)|2,limited.capacity);
        while (limited.table[slot].key!=((UINT64_C(1)<<32)|2))
            slot=(slot+1)&(limited.capacity-1);
        fail+=limited.table[slot].hits!=2;
        limited.table[slot].hits=SIZE_MAX;
        fail+=rawtex_pair_add(&limited,0,1)==0 || limited.table[slot].hits!=SIZE_MAX;
    }
    if (pairs.arena) Arena_dispose(&pairs.arena);
    if (limited.arena) Arena_dispose(&limited.arena);
    fprintf(stderr,"[selftest] rawtex sparse component counts %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    return fail ? -1 : 0;
}
