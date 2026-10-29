/* wind_audit.c -- winding-coordinate topology auditor for welded scroll meshes.
 *
 * Two topology defects matter and mesh_quality / manifold_check cannot see
 * either (both are per-triangle or per-edge-manifold tests, blind to how the
 * sheet is wired):
 *
 *   MERGER (inter-wrap fusion): an edge that jumps ACROSS the ~7-vox gap between
 *   two turns of the spiral -- a radial short-circuit. Fatal: it welds turn N to
 *   turn N+1 directly instead of via the 360-deg traverse, so the spiral can no
 *   longer be unwrapped. Signature: winding coordinate w = r/pitch - theta/2pi
 *   jumps by ~1 across a single edge (large radial delta at nearly constant
 *   angle).
 *
 *   SPLIT CONTACT (possible intra-sheet break): two DIFFERENT connected
 *   components occupy the same turn (same w) and sit within a small physical
 *   gap.  This is deliberately diagnostic rather than a hard defect gate:
 *   disk normalization and pinch repair duplicate boundary vertices on purpose,
 *   and clipped/source fragments can also approach without a certified zipper.
 *   Only a topology- and geometry-certified boundary transaction is repairable.
 *   Signature: a close vertex pair in two components with |dw| ~ 0.
 *
 *   BROAD FUSION (local-cube multi-wrap wall): a connected component spans more
 *   than one local winding even though every individual edge has a small |dw|.
 *   This catches the common "staircase" failure where BPA joins adjacent turns
 *   through many short triangles: no edge is a full-turn jump, so MERGER alone
 *   misses it. The test is intentionally opt-in (--local-span-tol) because a
 *   whole-scroll mesh is supposed to be one sheet spanning many windings. It is
 *   a semantic gate for leaf-cube BPA stages and larger exported regions away
 *   from the umbilicus.
 *
 *   LOCAL SHORTCUT (compact gradual bridge): an axial mesh section changes
 *   winding by a substantial fraction of a turn within only a short GEODESIC
 *   walk along its contour. This is the failure that defeats both tests above:
 *   every edge can advance by only 0.02 turn, and the whole 3-D component can
 *   stay below the broad-fusion threshold, while a 5-20 voxel staircase still
 *   joins neighbouring plies. Axial sections also remove the z-dependent axis
 *   wander that makes a global radius span a poor discriminator.
 *
 * The winding coordinate is LOCAL and branch-cut-free when computed pairwise:
 *   dr  = hypot(dy_v,dx_v) - hypot(dy_u,dx_u)      (dy,dx = coord - umbilicus)
 *   dth = wrap_to_pmpi( atan2(dy_v,dx_v) - atan2(dy_u,dx_u) )
 *   dw  = dr/pitch - dth/(2*pi)
 * exactly the seam gate's own test (seam_weld.c:536-540), so a merger edge here
 * is an edge the gate SHOULD have rejected (or that a non-bridge stage created).
 *
 * Vertex order is (z,y,x): coord[0]=z (axis), coord[1]=y, coord[2]=x. Radius is
 * in the (y,x) plane about the umbilicus. Turn count spanned by a component is
 * (r_max - r_min)/pitch -- the single stat that says whether a giant component
 * is a legit continuous spiral (few merger edges) or a fused pile (many).
 *
 *   wind_audit <in.obj> [--axis-table z_y_x.csv | --umb-y Y --umb-x X]
 *              [--pitch P=9.5] [--top N=25]
 *              [--merge-tol T=0.40] [--split-wtol T=0.25] [--split-gap G=14]
 *              [--no-split-audit]
 *              [--core-pitches C=3] [--min-faces M=64]
 *              [--local-span-tol T] [--fail-on-local]
 *              [--shortcut-tol T] [--shortcut-length-pitches L=2]
 *              [--shortcut-z-step Z=4] [--shortcut-min-planes N=3]
 *              [--shortcut-track-mutual]
 *              [--fail-on-shortcut]
 *              [--dump-mergers f.obj] [--dump-splits f.obj]
 *              [--dump-shortcuts f.obj]
 *              [--repair-seams f.obj] [--repair-radius-pitches R=1.5]
 *              [--repair-dry-run] [--repair-max-cut-fraction F=.05]
 *              [--repair-cover-all-supported] [--repair-min-chain-support N=1]
 *              [--repair-monotone-step]
 *              [--repair-phase-offset W=0] [--repair-auto-offset]
 *              [--repair-max-auto-tracks N=4096; 0=unlimited]
 *   wind_audit --selftest
 *
 * Standalone C99, own OBJ parser (shared shape with obj_cc_color), no deps.
 * Exit: 0 ok / selftest pass, 1 IO error, 2 usage error, 3 selftest fail,
 *       4 audit gate failed, 5 requested repair rejected.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <errno.h>

#include "../common/pitch_table.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)

static int parse_size_arg(const char *text, size_t *out)
{
    char *end = NULL;
    unsigned long long value;
    if (text == NULL || text[0] == '\0' || text[0] == '-') return 0;
    errno = 0;
    value = strtoull(text, &end, 10);
    if (errno == ERANGE || end == text || *end != '\0' ||
        value > (unsigned long long)SIZE_MAX) return 0;
    *out = (size_t)value;
    return 1;
}

/* ---------- growing arrays ---------- */
typedef struct { float *v; size_t n, cap; } FVec;   /* verts: 3 floats each */
typedef struct { int32_t *f; size_t n, cap; } IVec; /* faces: 3 ints each   */

static int fv_push(FVec *a, float x, float y, float z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        float *np = (float *)realloc(a->v, nc * 3 * sizeof(float)); if (!np) return -1; a->v = np; a->cap = nc; }
    a->v[a->n*3+0] = x; a->v[a->n*3+1] = y; a->v[a->n*3+2] = z; a->n++; return 0;
}
static int iv_push(IVec *a, int32_t x, int32_t y, int32_t z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        int32_t *np = (int32_t *)realloc(a->f, nc * 3 * sizeof(int32_t)); if (!np) return -1; a->f = np; a->cap = nc; }
    a->f[a->n*3+0] = x; a->f[a->n*3+1] = y; a->f[a->n*3+2] = z; a->n++; return 0;
}

static int32_t parse_findex(const char *tok) { return (int32_t)atol(tok); }

static int read_obj(const char *path, FVec *V, IVec *F)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "wind_audit: cannot open %s\n", path); return -1; }
    size_t lcap = 1u << 16; char *line = (char *)malloc(lcap);
    int rc = 0;
    if (!line) { fclose(fp); return -1; }
    while (fgets(line, (int)lcap, fp)) {
        while (!strchr(line, '\n') && !feof(fp)) {
            size_t len = strlen(line);
            char *nl = (char *)realloc(line, lcap * 2); if (!nl) { rc = -1; goto done; }
            line = nl; lcap *= 2;
            if (!fgets(line + len, (int)(lcap - len), fp)) break;
        }
        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            double x = 0, y = 0, z = 0;
            if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3)
                if (fv_push(V, (float)x, (float)y, (float)z) != 0) { rc = -1; goto done; }
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char *tok = strtok(line + 2, " \t\r\n");
            int32_t idx[3]; int n = 0;
            while (tok && n < 3) { idx[n++] = parse_findex(tok); tok = strtok(NULL, " \t\r\n"); }
            if (n == 3) {
                int32_t a = idx[0], b = idx[1], c = idx[2];
                if (a < 0) a = (int32_t)V->n + a + 1;
                if (b < 0) b = (int32_t)V->n + b + 1;
                if (c < 0) c = (int32_t)V->n + c + 1;
                if (iv_push(F, a - 1, b - 1, c - 1) != 0) { rc = -1; goto done; }
            }
        }
    }
done:
    free(line); fclose(fp);
    return rc;
}

static int write_obj(const char *path, const FVec *V, const IVec *F)
{
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) {
        fprintf(stderr, "wind_audit: cannot write %s\n", path);
        return -1;
    }
    fprintf(fp, "# wind_audit supported local-shortcut seam repair\n");
    for (size_t i = 0; i < V->n; i++)
        fprintf(fp, "v %.9g %.9g %.9g\n",
                V->v[i*3+0], V->v[i*3+1], V->v[i*3+2]);
    for (size_t i = 0; i < F->n; i++)
        fprintf(fp, "f %d %d %d\n", F->f[i*3+0]+1,
                F->f[i*3+1]+1, F->f[i*3+2]+1);
    if (fclose(fp) != 0) return -1;
    return 0;
}

static void write_octahedron_marker(FILE *fp, const double p[3],
                                    double radius,
                                    double red, double green, double blue,
                                    size_t *vertex_count)
{
    static const int offset[6][3] = {
        { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 },
        { 0,-1, 0 }, {  0, 0, 1 }, { 0, 0,-1 }
    };
    static const int face[8][3] = {
        {0,2,4}, {0,4,3}, {0,3,5}, {0,5,2},
        {1,4,2}, {1,3,4}, {1,5,3}, {1,2,5}
    };
    size_t base = *vertex_count;
    for (int i = 0; i < 6; i++)
        fprintf(fp, "v %.7g %.7g %.7g %.2f %.2f %.2f\n",
                p[0] + radius * offset[i][0],
                p[1] + radius * offset[i][1],
                p[2] + radius * offset[i][2], red, green, blue);
    for (int i = 0; i < 8; i++)
        fprintf(fp, "f %zu %zu %zu\n", base+(size_t)face[i][0]+1,
                base+(size_t)face[i][1]+1,
                base+(size_t)face[i][2]+1);
    *vertex_count += 6;
}

/* ---------- union-find (path halving + union by size) ---------- */
static int32_t uf_find(int32_t *p, int32_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
static void uf_union(int32_t *p, int32_t *sz, int32_t a, int32_t b)
{
    a = uf_find(p, a); b = uf_find(p, b);
    if (a == b) return;
    if (sz[a] < sz[b]) { int32_t t = a; a = b; b = t; }
    p[b] = a; sz[a] += sz[b];
}

/* ---------- winding helpers (coord order z,y,x; umbilicus in y,x) ---------- */
typedef struct {
    double *z, *y, *x;
    size_t n;
} AxisTable;

static void axis_table_free(AxisTable *table)
{
    free(table->z); free(table->y); free(table->x);
    memset(table, 0, sizeof(*table));
}

static int axis_table_parse_line(const char *line,
                                 double *z, double *y, double *x)
{
    char *end = NULL;
    const char *p = line;
    *z = strtod(p, &end);
    if (end == p) return 0; /* header/comment */
    p = end;
    while (*p == ' ' || *p == '\t') p++;
    if (*p++ != ',') return 0;
    *y = strtod(p, &end);
    if (end == p) return 0;
    p = end;
    while (*p == ' ' || *p == '\t') p++;
    if (*p++ != ',') return 0;
    *x = strtod(p, &end);
    if (end == p) return 0;
    return isfinite(*z) && isfinite(*y) && isfinite(*x);
}

/* Load z,y,x rows from CSV. Rows are sorted so either raw samples or a
 * pre-smoothed table can be supplied. Duplicate z samples are rejected. */
static int axis_table_load(const char *path, AxisTable *table)
{
    FILE *f;
    size_t cap = 32;
    char line[1024];
    memset(table, 0, sizeof(*table));
    f = fopen(path, "r");
    if (f == NULL) return -1;
    table->z = (double *)malloc(cap * sizeof(double));
    table->y = (double *)malloc(cap * sizeof(double));
    table->x = (double *)malloc(cap * sizeof(double));
    if (table->z == NULL || table->y == NULL || table->x == NULL) {
        fclose(f); axis_table_free(table); return -1;
    }
    while (fgets(line, sizeof(line), f) != NULL) {
        double z = 0.0, y = 0.0, x = 0.0;
        if (!axis_table_parse_line(line, &z, &y, &x)) continue;
        if (table->n == cap) {
            double *next;
            cap *= 2;
            next = (double *)realloc(table->z, cap * sizeof(double));
            if (next == NULL) {
                fclose(f); axis_table_free(table); return -1;
            }
            table->z = next;
            next = (double *)realloc(table->y, cap * sizeof(double));
            if (next == NULL) {
                fclose(f); axis_table_free(table); return -1;
            }
            table->y = next;
            next = (double *)realloc(table->x, cap * sizeof(double));
            if (next == NULL) {
                fclose(f); axis_table_free(table); return -1;
            }
            table->x = next;
        }
        table->z[table->n] = z;
        table->y[table->n] = y;
        table->x[table->n] = x;
        table->n++;
    }
    fclose(f);
    if (table->n < 2) {
        axis_table_free(table);
        return -1;
    }
    for (size_t i = 1; i < table->n; i++) {
        double z = table->z[i], y = table->y[i], x = table->x[i];
        size_t j = i;
        while (j > 0 && table->z[j-1] > z) {
            table->z[j] = table->z[j-1];
            table->y[j] = table->y[j-1];
            table->x[j] = table->x[j-1];
            j--;
        }
        table->z[j] = z; table->y[j] = y; table->x[j] = x;
    }
    for (size_t i = 1; i < table->n; i++) {
        if (!(table->z[i] > table->z[i-1])) {
            axis_table_free(table);
            return -1;
        }
    }
    return 0;
}

/* Piecewise-linear axis with endpoint tangent extrapolation. The PHerc crop
 * extends half a cube beyond its first/last centre samples; clamping there
 * creates an artificial kink in radius and winding. */
static void axis_table_eval(const AxisTable *table, double z,
                            double fallback_y, double fallback_x,
                            double *out_y, double *out_x)
{
    size_t lo, hi;
    double t;
    if (table == NULL || table->n < 2) {
        *out_y = fallback_y; *out_x = fallback_x; return;
    }
    if (z <= table->z[0]) {
        lo = 0; hi = 1;
    } else if (z >= table->z[table->n-1]) {
        lo = table->n - 2; hi = table->n - 1;
    } else {
        lo = 0; hi = table->n - 1;
        while (hi - lo > 1) {
            size_t mid = lo + (hi - lo) / 2;
            if (table->z[mid] <= z) lo = mid;
            else hi = mid;
        }
    }
    t = (z - table->z[lo]) / (table->z[hi] - table->z[lo]);
    *out_y = table->y[lo] + t * (table->y[hi] - table->y[lo]);
    *out_x = table->x[lo] + t * (table->x[hi] - table->x[lo]);
}

static double wrap_pmpi(double a) { while (a > M_PI) a -= TWO_PI; while (a <= -M_PI) a += TWO_PI; return a; }

/* per-vertex radius about the axis (y,x plane) */
static double vradius(const float *V, size_t i, double uy, double ux,
                      const AxisTable *axis)
{
    axis_table_eval(axis, (double)V[i*3+0], uy, ux, &uy, &ux);
    double dy = (double)V[i*3+1] - uy, dx = (double)V[i*3+2] - ux;
    return hypot(dy, dx);
}
/* pairwise winding delta w_v - w_u (branch-cut-free); the radial term is the
 * measured-table integral when a table is armed, else dr/pitch exactly. */
static double dwind(const float *V, size_t u, size_t v, double uy, double ux,
                    double pitch, const AxisTable *axis, PitchTable ptab)
{
    double uyu, uxu, uyv, uxv;
    axis_table_eval(axis, (double)V[u*3+0], uy, ux, &uyu, &uxu);
    axis_table_eval(axis, (double)V[v*3+0], uy, ux, &uyv, &uxv);
    double dyu = (double)V[u*3+1] - uyu, dxu = (double)V[u*3+2] - uxu;
    double dyv = (double)V[v*3+1] - uyv, dxv = (double)V[v*3+2] - uxv;
    double ru  = hypot(dyu, dxu), rv = hypot(dyv, dxv);
    double dth = wrap_pmpi(atan2(dyv, dxv) - atan2(dyu, dxu));
    return PitchTable_dturns(ptab, ru, rv, pitch) - dth / TWO_PI;
}

/* ---------- component table ---------- */
typedef struct {
    int32_t root; size_t faces; int32_t rank;
    size_t vertices, edges, boundary_edges, boundary_loops;
    size_t nonmanifold_edges, boundary_irregular_vertices;
    int64_t euler_chi;
    double genus;
    double rmin, rmax; float zmin, zmax;
    size_t merger_edges;   /* reliable (outer) merger edges wholly inside comp */
    double phase_min, phase_max; /* adjacency-unwrapped r/pitch - theta/(2pi) */
    size_t outer_vertices;
} Comp;
static int comp_cmp(const void *pa, const void *pb)
{
    const Comp *a = (const Comp *)pa, *b = (const Comp *)pb;
    if (a->faces != b->faces) return a->faces < b->faces ? 1 : -1;
    return a->root < b->root ? -1 : 1;
}

/* ---------- spatial hash over vertices (cell = split_gap) ---------- */
typedef struct { int64_t key; int32_t head; } HCell;
typedef struct { HCell *cell; size_t mask; int32_t *next; double inv; int64_t off; } SHash;

static int64_t sh_key(const SHash *h, double a, double b, double c)
{
    int64_t ia = (int64_t)floor(a * h->inv) + h->off;
    int64_t ib = (int64_t)floor(b * h->inv) + h->off;
    int64_t ic = (int64_t)floor(c * h->inv) + h->off;
    /* pack into 63 bits: 21 bits each (range 0..2M cells per axis) */
    return ((ia & 0x1FFFFF) << 42) | ((ib & 0x1FFFFF) << 21) | (ic & 0x1FFFFF);
}
static size_t sh_hash(int64_t k, size_t mask)
{
    uint64_t x = (uint64_t)k * 0x9E3779B97F4A7C15ull;
    x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ull; x ^= x >> 32;
    return (size_t)x & mask;
}
static int32_t *sh_slot(SHash *h, int64_t key, int create)
{
    size_t i = sh_hash(key, h->mask);
    for (;;) {
        if (h->cell[i].head == -2) { /* empty */
            if (!create) return NULL;
            h->cell[i].key = key; h->cell[i].head = -1; return &h->cell[i].head;
        }
        if (h->cell[i].key == key) return &h->cell[i].head;
        i = (i + 1) & h->mask;
    }
}

/* ---------- reporting core ---------- */
typedef struct {
    double umb_y, umb_x, pitch;
    const AxisTable *axis;
    const char *axis_table_path;
    PitchTable ptab;               /* measured radial pitch (NULL = scalar) */
    const char *pitch_table_path;
    double merge_tol, split_wtol, split_gap, core_pitches;
    double local_span_tol;       /* 0 disables component broad-fusion test */
    double shortcut_tol;         /* 0 disables axial local-shortcut test */
    double shortcut_length_pitches;
    double shortcut_z_step;
    double shortcut_cluster_pitches;
    double repair_radius_pitches;
    double repair_support_radius_pitches;
    double repair_band_margin;
    double repair_phase_offset;
    double repair_max_cut_fraction;
    double seam_pitch, seam_tol;   /* leaf-cube seam planes (multiples of seam_pitch) */
    double boundary_box[6];        /* z0,z1,y0,y1,x0,x1 for clipped-track support */
    size_t min_faces, shortcut_min_points, repair_min_chain_support;
    size_t repair_max_auto_tracks;
    int top, shortcut_min_planes;
    int fail_on_local, fail_on_shortcut, repair_dry_run;
    int repair_emit_rejected, repair_cover_all_supported, repair_auto_offset;
    int repair_monotone_step;
    int shortcut_track_mutual, boundary_box_set, boundary_box_from_mesh, split_audit;
    const char *dump_mergers, *dump_splits, *dump_shortcuts, *json_path;
    const char *repair_chain_csv;
    const char *repair_path;
} Params;

/* absolute winding potential of a radius under the active gauge.  The zero
 * point is arbitrary; only differences and spans are consumed, so the scalar
 * and table gauges are interchangeable consumers of this helper. */
static double rturns(const Params *P, double r)
{
    return PitchTable_dturns(P->ptab, 0.0, r, P->pitch);
}

/* ---------- axial-section local-shortcut audit ---------- */
typedef struct {
    double p[3];                 /* midpoint of the worst short geodesic */
    double a[3], b[3];           /* its two contour endpoints */
    double delta_w, signed_delta_w, distance, plane;
    size_t curve_points;
    int persistent, boundary_clipped;
    int32_t cluster_id;
} ShortcutHit;

typedef struct {
    size_t hits, planes, boundary_hits;
    double z_min, z_max, max_delta_w, min_distance, excess;
    double centroid[3], bbox_min[3], bbox_max[3];
    int persistent, boundary_supported;
    int32_t source_id;           /* pre-sort id; remaps hit.cluster_id */
} ShortcutCluster;

typedef struct {
    ShortcutHit *hit;
    size_t nhit, hit_cap;
    ShortcutCluster *cluster;
    size_t ncluster;
    size_t planes_tested, curves_tested, section_face_incidents;
    size_t simple_curves, fallback_curves;
    size_t persistent_clusters, persistent_hits;
    double max_delta_w, persistent_max_delta_w, persistent_excess;
} ShortcutAudit;

typedef struct { int32_t a, b; double length; } SectionEdge;
typedef struct { double d; int32_t v; } SectionHeapNode;

/* A shortcut audit samples many parallel z sections.  Scanning every face for
 * every section made the terminal hierarchy O(nplanes * nfaces), even though a
 * triangle normally crosses only one or two sampled planes.  This exact
 * incidence index stores each face in precisely the planes satisfying the
 * original strict predicate (zmin < plane && zmax > plane). */
typedef struct {
    double *plane;
    size_t *row;
    size_t *face;
    size_t nplane;
} SectionFaceIndex;

static void section_face_index_free(SectionFaceIndex *index)
{
    if (index == NULL) return;
    free(index->plane);
    free(index->row);
    free(index->face);
    memset(index, 0, sizeof(*index));
}

static void section_plane_span(const SectionFaceIndex *index,
                               double first, double step,
                               double lo, double hi,
                               size_t *begin, size_t *end)
{
    double gb = floor((lo - first) / step) + 1.0;
    double ge = ceil((hi - first) / step);
    size_t b = gb <= 0.0 ? 0 :
               (gb >= (double)index->nplane ? index->nplane : (size_t)gb);
    size_t e = ge <= 0.0 ? 0 :
               (ge >= (double)index->nplane ? index->nplane : (size_t)ge);

    /* Correct the O(1) estimates against the stored plane values.  The final
     * comparisons are exactly those used by the old full scan, including its
     * strict treatment of coplanar vertices. */
    while (b > 0 && index->plane[b - 1] > lo) b--;
    while (b < index->nplane && !(index->plane[b] > lo)) b++;
    while (e > 0 && index->plane[e - 1] >= hi) e--;
    while (e < index->nplane && index->plane[e] < hi) e++;
    if (e < b) e = b;
    *begin = b;
    *end = e;
}

static int section_face_index_build(const FVec *V, const IVec *F,
                                    double first, double step, double z_max,
                                    SectionFaceIndex *index)
{
    size_t *count = NULL;
    size_t total = 0;
    double plane = first;
    memset(index, 0, sizeof(*index));

    while (plane < z_max) {
        double next = plane + step;
        if (!(next > plane) || index->nplane == SIZE_MAX) return -1;
        index->nplane++;
        plane = next;
    }
    if (index->nplane == 0) return 0;
    if (index->nplane > SIZE_MAX / sizeof(*index->plane) ||
        index->nplane + 1 < index->nplane ||
        index->nplane + 1 > SIZE_MAX / sizeof(*index->row))
        return -1;
    index->plane = (double *)malloc(index->nplane * sizeof(*index->plane));
    index->row = (size_t *)malloc((index->nplane + 1) * sizeof(*index->row));
    count = (size_t *)calloc(index->nplane, sizeof(*count));
    if (index->plane == NULL || index->row == NULL || count == NULL)
        goto oom;
    plane = first;
    for (size_t p = 0; p < index->nplane; p++) {
        index->plane[p] = plane;
        plane += step;
    }

    for (size_t f = 0; f < F->n; f++) {
        int32_t a = F->f[f*3+0], b = F->f[f*3+1], c = F->f[f*3+2];
        double lo, hi, zb, zc;
        size_t begin, end;
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= V->n || (size_t)b >= V->n ||
            (size_t)c >= V->n) continue;
        lo = hi = (double)V->v[(size_t)a*3+0];
        zb = (double)V->v[(size_t)b*3+0];
        zc = (double)V->v[(size_t)c*3+0];
        if (zb < lo) lo = zb;
        if (zb > hi) hi = zb;
        if (zc < lo) lo = zc;
        if (zc > hi) hi = zc;
        if (!(lo < hi) || !isfinite(lo) || !isfinite(hi)) continue;
        section_plane_span(index, first, step, lo, hi, &begin, &end);
        for (size_t p = begin; p < end; p++) {
            if (!(lo < index->plane[p] && hi > index->plane[p])) continue;
            if (count[p] == SIZE_MAX) goto oom;
            count[p]++;
        }
    }
    index->row[0] = 0;
    for (size_t p = 0; p < index->nplane; p++) {
        if (count[p] > SIZE_MAX - total) goto oom;
        total += count[p];
        index->row[p + 1] = total;
        count[p] = 0;
    }
    if (total > SIZE_MAX / sizeof(*index->face)) goto oom;
    index->face = (size_t *)malloc((total ? total : 1) *
                                   sizeof(*index->face));
    if (index->face == NULL) goto oom;

    for (size_t f = 0; f < F->n; f++) {
        int32_t a = F->f[f*3+0], b = F->f[f*3+1], c = F->f[f*3+2];
        double lo, hi, zb, zc;
        size_t begin, end;
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= V->n || (size_t)b >= V->n ||
            (size_t)c >= V->n) continue;
        lo = hi = (double)V->v[(size_t)a*3+0];
        zb = (double)V->v[(size_t)b*3+0];
        zc = (double)V->v[(size_t)c*3+0];
        if (zb < lo) lo = zb;
        if (zb > hi) hi = zb;
        if (zc < lo) lo = zc;
        if (zc > hi) hi = zc;
        if (!(lo < hi) || !isfinite(lo) || !isfinite(hi)) continue;
        section_plane_span(index, first, step, lo, hi, &begin, &end);
        for (size_t p = begin; p < end; p++) {
            if (!(lo < index->plane[p] && hi > index->plane[p])) continue;
            index->face[index->row[p] + count[p]++] = f;
        }
    }
    free(count);
    return 0;

oom:
    free(count);
    section_face_index_free(index);
    return -1;
}

/* A manifold triangle mesh intersects a regular plane in disjoint paths and
 * cycles.  On those degree-<=2 graphs, truncated shortest paths are obtained by
 * walking in the two contour directions; running a binary-heap Dijkstra from
 * every sample is equivalent but needlessly superlinear.  Return 1 when the
 * exact simple-contour path was used and 0 when a branch, duplicate adjacency,
 * or short-cycle ambiguity requires the general Dijkstra fallback. */
static int section_simple_shortcut(const int32_t *member,
                                   size_t member_begin, size_t member_end,
                                   const size_t *row, const int32_t *col,
                                   const double *weight, const uint8_t *outer,
                                   const double *radius, const double *theta,
                                   double max_path, double pitch,
                                   double *best_delta,
                                   double *best_signed_delta,
                                   double *best_distance,
                                   int32_t *best_a, int32_t *best_b)
{
    size_t endpoints = 0;
    double total_length = 0.0;

    for (size_t at = member_begin; at < member_end; at++) {
        int32_t u = member[at], neighbor[2] = {-1, -1};
        size_t ndegree = 0;
        for (size_t e = row[(size_t)u]; e < row[(size_t)u + 1]; e++) {
            int32_t v = col[e];
            if (!outer[(size_t)v]) continue;
            if (ndegree >= 2 || (ndegree == 1 && neighbor[0] == v)) return 0;
            neighbor[ndegree++] = v;
            total_length += weight[e];
        }
        if (ndegree == 1) endpoints++;
        else if (ndegree != 2) return 0;
    }
    if (endpoints != 0 && endpoints != 2) return 0;
    if (endpoints == 0) {
        /* Each cycle edge was counted in both directions.  If the admissible
         * radius reaches a half-cycle, equal-length routes make the old heap's
         * tie choice observable; retain the fallback in that rare case. */
        total_length *= 0.5;
        if (!(total_length > 0.0) ||
            max_path >= 0.5 * total_length - 1e-9)
            return 0;
    }

    for (size_t sat = member_begin; sat < member_end; sat++) {
        int32_t source = member[sat];
        for (size_t first_edge = row[(size_t)source];
             first_edge < row[(size_t)source + 1]; first_edge++) {
            int32_t previous = source;
            int32_t u = col[first_edge];
            double distance, path_phase;
            if (!outer[(size_t)u]) continue;
            distance = weight[first_edge];
            path_phase =
                (radius[(size_t)u] - radius[(size_t)source]) / pitch -
                wrap_pmpi(theta[(size_t)u] - theta[(size_t)source]) / TWO_PI;
            while (distance <= max_path) {
                double delta = fabs(path_phase);
                int32_t next = -1;
                double next_weight = 0.0;
                if (delta > *best_delta) {
                    *best_delta = delta;
                    *best_signed_delta = path_phase;
                    *best_distance = distance;
                    *best_a = source;
                    *best_b = u;
                }
                for (size_t e = row[(size_t)u]; e < row[(size_t)u + 1]; e++) {
                    int32_t v = col[e];
                    if (!outer[(size_t)v] || v == previous) continue;
                    next = v;
                    next_weight = weight[e];
                    break;
                }
                if (next < 0 || next == source ||
                    distance + next_weight > max_path)
                    break;
                path_phase +=
                    (radius[(size_t)next] - radius[(size_t)u]) / pitch -
                    wrap_pmpi(theta[(size_t)next] - theta[(size_t)u]) / TWO_PI;
                distance += next_weight;
                previous = u;
                u = next;
            }
        }
    }
    return 1;
}

static int near_seam(double z, double y, double x, double pitch, double tol);
static int near_boundary_box(double z, double y, double x, const Params *P);

static size_t hash_u64(uint64_t key, size_t mask)
{
    key ^= key >> 30;
    key *= UINT64_C(0xBF58476D1CE4E5B9);
    key ^= key >> 27;
    key *= UINT64_C(0x94D049BB133111EB);
    key ^= key >> 31;
    return (size_t)key & mask;
}

static int shortcut_hit_push(ShortcutAudit *a, const ShortcutHit *h)
{
    if (a->nhit == a->hit_cap) {
        size_t cap = a->hit_cap ? a->hit_cap * 2 : 256;
        ShortcutHit *next = (ShortcutHit *)realloc(
            a->hit, cap * sizeof(*next));
        if (next == NULL) return -1;
        a->hit = next;
        a->hit_cap = cap;
    }
    a->hit[a->nhit++] = *h;
    return 0;
}

static int section_heap_push(SectionHeapNode **heap, size_t *n, size_t *cap,
                             double d, int32_t v)
{
    if (*n == *cap) {
        size_t next_cap = *cap ? *cap * 2 : 256;
        SectionHeapNode *next = (SectionHeapNode *)realloc(
            *heap, next_cap * sizeof(*next));
        if (next == NULL) return -1;
        *heap = next;
        *cap = next_cap;
    }
    size_t i = (*n)++;
    (*heap)[i].d = d;
    (*heap)[i].v = v;
    while (i > 0) {
        size_t parent = (i - 1) / 2;
        if ((*heap)[parent].d <= (*heap)[i].d) break;
        SectionHeapNode tmp = (*heap)[parent];
        (*heap)[parent] = (*heap)[i];
        (*heap)[i] = tmp;
        i = parent;
    }
    return 0;
}

static int section_heap_pop(SectionHeapNode *heap, size_t *n,
                            SectionHeapNode *out)
{
    if (*n == 0) return 0;
    *out = heap[0];
    (*n)--;
    if (*n == 0) return 1;
    heap[0] = heap[*n];
    size_t i = 0;
    for (;;) {
        size_t left = i * 2 + 1, right = left + 1, best = i;
        if (left < *n && heap[left].d < heap[best].d) best = left;
        if (right < *n && heap[right].d < heap[best].d) best = right;
        if (best == i) break;
        SectionHeapNode tmp = heap[i];
        heap[i] = heap[best];
        heap[best] = tmp;
        i = best;
    }
    return 1;
}

static int shortcut_cluster_cmp(const void *pa, const void *pb)
{
    const ShortcutCluster *a = (const ShortcutCluster *)pa;
    const ShortcutCluster *b = (const ShortcutCluster *)pb;
    if (a->persistent != b->persistent)
        return a->persistent ? -1 : 1;
    if (a->planes != b->planes) return a->planes < b->planes ? 1 : -1;
    if (a->max_delta_w != b->max_delta_w)
        return a->max_delta_w < b->max_delta_w ? 1 : -1;
    return a->hits < b->hits ? 1 : (a->hits > b->hits ? -1 : 0);
}

static void shortcut_audit_free(ShortcutAudit *a)
{
    free(a->hit);
    free(a->cluster);
    memset(a, 0, sizeof(*a));
}

/* Intersect the mesh with axial planes, join section segments through shared
 * mesh edges, then search each contour for the largest branch-local winding
 * change reachable within a bounded geodesic distance. A false bridge is a
 * winding SHORTCUT; a real spiral continuation must travel around the scroll
 * and cannot realize the same delta inside this local distance budget. */
static int shortcut_audit(const FVec *V, const IVec *F, const Params *P,
                          ShortcutAudit *out)
{
    memset(out, 0, sizeof(*out));
    if (P->shortcut_tol <= 0.0 || V->n == 0 || F->n == 0) return 0;

    double z_min = DBL_MAX, z_max = -DBL_MAX;
    for (size_t i = 0; i < V->n; i++) {
        double z = (double)V->v[i * 3 + 0];
        if (z < z_min) z_min = z;
        if (z > z_max) z_max = z;
    }
    double step = P->shortcut_z_step;
    double max_path = P->shortcut_length_pitches * P->pitch;
    double core_r = P->core_pitches * P->pitch;
    if (step <= 0.0 || max_path <= 0.0 || z_max <= z_min) return 0;

    /* Keep sections away from integer and half-integer mesh layers, avoiding
     * coplanar-edge special cases without changing the sampling interval. */
    double first = floor(z_min) + 0.371;
    while (first <= z_min) first += step;

    SectionFaceIndex section_index;
    if (section_face_index_build(V, F, first, step, z_max,
                                 &section_index) != 0) {
        shortcut_audit_free(out);
        return -1;
    }
    if (section_index.nplane > 0)
        out->section_face_incidents = section_index.row[section_index.nplane];

    for (size_t plane_i = 0; plane_i < section_index.nplane; plane_i++) {
        double plane = section_index.plane[plane_i];
        size_t face_begin = section_index.row[plane_i];
        size_t face_end = section_index.row[plane_i + 1];
        size_t ncross = face_end - face_begin;
        if (ncross == 0) continue;
        out->planes_tested++;

        size_t max_points = ncross * 2;
        double *point = (double *)malloc(max_points * 3 * sizeof(*point));
        SectionEdge *edge = (SectionEdge *)malloc(ncross * sizeof(*edge));
        size_t slots = 1;
        while (slots < ncross * 8) slots <<= 1;
        uint64_t *edge_key = (uint64_t *)malloc(slots * sizeof(*edge_key));
        int32_t *edge_point = (int32_t *)malloc(slots * sizeof(*edge_point));
        if (point == NULL || edge == NULL || edge_key == NULL ||
            edge_point == NULL) {
            free(point); free(edge); free(edge_key); free(edge_point);
            section_face_index_free(&section_index);
            shortcut_audit_free(out);
            return -1;
        }
        for (size_t i = 0; i < slots; i++) edge_key[i] = UINT64_MAX;
        size_t npoint = 0, nedge = 0, mask = slots - 1;

        for (size_t active = face_begin; active < face_end; active++) {
            size_t f = section_index.face[active];
            int32_t tri[3] = {
                F->f[f*3+0], F->f[f*3+1], F->f[f*3+2]
            };
            if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0 ||
                (size_t)tri[0] >= V->n || (size_t)tri[1] >= V->n ||
                (size_t)tri[2] >= V->n) continue;
            int32_t hit[3];
            int nhit = 0;
            for (int k = 0; k < 3; k++) {
                int32_t a = tri[k], b = tri[(k + 1) % 3];
                double za = (double)V->v[(size_t)a*3+0];
                double zb = (double)V->v[(size_t)b*3+0];
                if (!((za < plane && zb > plane) ||
                      (zb < plane && za > plane))) continue;
                uint32_t lo = (uint32_t)(a < b ? a : b);
                uint32_t hi = (uint32_t)(a < b ? b : a);
                uint64_t key = ((uint64_t)lo << 32) | (uint64_t)hi;
                size_t slot = hash_u64(key, mask);
                while (edge_key[slot] != UINT64_MAX &&
                       edge_key[slot] != key)
                    slot = (slot + 1) & mask;
                int32_t pi = -1;
                if (edge_key[slot] == key) {
                    pi = edge_point[slot];
                } else {
                    if (npoint >= max_points) continue;
                    edge_key[slot] = key;
                    pi = (int32_t)npoint;
                    edge_point[slot] = pi;
                    double t = (plane - za) / (zb - za);
                    for (int axis = 0; axis < 3; axis++) {
                        double va = (double)V->v[(size_t)a*3+(size_t)axis];
                        double vb = (double)V->v[(size_t)b*3+(size_t)axis];
                        point[npoint*3+(size_t)axis] = va + t * (vb - va);
                    }
                    npoint++;
                }
                if (nhit < 3) hit[nhit++] = pi;
            }
            if (nhit == 2 && hit[0] != hit[1] && nedge < ncross) {
                double dz = point[(size_t)hit[0]*3+0] -
                            point[(size_t)hit[1]*3+0];
                double dy = point[(size_t)hit[0]*3+1] -
                            point[(size_t)hit[1]*3+1];
                double dx = point[(size_t)hit[0]*3+2] -
                            point[(size_t)hit[1]*3+2];
                edge[nedge].a = hit[0];
                edge[nedge].b = hit[1];
                edge[nedge].length = sqrt(dz*dz + dy*dy + dx*dx);
                nedge++;
            }
        }
        free(edge_key);
        free(edge_point);
        if (npoint == 0 || nedge == 0) {
            free(point); free(edge);
            continue;
        }

        int32_t *par = (int32_t *)malloc(npoint * sizeof(*par));
        int32_t *usz = (int32_t *)malloc(npoint * sizeof(*usz));
        size_t *degree = (size_t *)calloc(npoint, sizeof(*degree));
        size_t *row = (size_t *)malloc((npoint + 1) * sizeof(*row));
        int32_t *col = (int32_t *)malloc(nedge * 2 * sizeof(*col));
        double *weight = (double *)malloc(nedge * 2 * sizeof(*weight));
        size_t *member_count = (size_t *)calloc(npoint, sizeof(*member_count));
        double *radius = (double *)malloc(npoint * sizeof(*radius));
        double *theta = (double *)malloc(npoint * sizeof(*theta));
        uint8_t *outer = (uint8_t *)malloc(npoint);
        if (par == NULL || usz == NULL || degree == NULL || row == NULL ||
            col == NULL || weight == NULL || member_count == NULL ||
            radius == NULL || theta == NULL || outer == NULL) {
            free(par); free(usz); free(degree); free(row); free(col);
            free(weight); free(member_count); free(radius); free(theta);
            free(outer); free(point); free(edge);
            section_face_index_free(&section_index);
            shortcut_audit_free(out);
            return -1;
        }
        for (size_t i = 0; i < npoint; i++) {
            double uy, ux;
            par[i] = (int32_t)i;
            usz[i] = 1;
            axis_table_eval(P->axis, point[i*3+0], P->umb_y, P->umb_x,
                            &uy, &ux);
            double dy = point[i*3+1] - uy;
            double dx = point[i*3+2] - ux;
            radius[i] = hypot(dy, dx);
            theta[i] = atan2(dy, dx);
            outer[i] = radius[i] >= core_r;
        }
        for (size_t i = 0; i < nedge; i++) {
            if (outer[(size_t)edge[i].a] && outer[(size_t)edge[i].b])
                uf_union(par, usz, edge[i].a, edge[i].b);
            degree[(size_t)edge[i].a]++;
            degree[(size_t)edge[i].b]++;
        }
        row[0] = 0;
        for (size_t i = 0; i < npoint; i++) row[i+1] = row[i] + degree[i];
        size_t *fill = (size_t *)malloc(npoint * sizeof(*fill));
        if (fill == NULL) {
            free(par); free(usz); free(degree); free(row); free(col);
            free(weight); free(member_count); free(radius); free(theta);
            free(outer); free(point); free(edge);
            section_face_index_free(&section_index);
            shortcut_audit_free(out);
            return -1;
        }
        memcpy(fill, row, npoint * sizeof(*fill));
        for (size_t i = 0; i < nedge; i++) {
            int32_t a = edge[i].a, b = edge[i].b;
            size_t ia = fill[(size_t)a]++, ib = fill[(size_t)b]++;
            col[ia] = b; weight[ia] = edge[i].length;
            col[ib] = a; weight[ib] = edge[i].length;
        }
        free(fill);
        free(edge);

        for (size_t i = 0; i < npoint; i++) {
            int32_t root = uf_find(par, (int32_t)i);
            member_count[(size_t)root]++;
        }
        size_t *member_row = (size_t *)calloc(
            npoint + 1, sizeof(*member_row));
        int32_t *member = (int32_t *)malloc(npoint * sizeof(*member));
        double *dist = (double *)malloc(npoint * sizeof(*dist));
        double *path_w = (double *)malloc(npoint * sizeof(*path_w));
        if (member_row == NULL || member == NULL || dist == NULL ||
            path_w == NULL) {
            free(member_row); free(member); free(dist); free(path_w);
            free(par); free(usz); free(degree); free(row); free(col);
            free(weight); free(member_count); free(radius); free(theta);
            free(outer); free(point);
            section_face_index_free(&section_index);
            shortcut_audit_free(out);
            return -1;
        }
        for (size_t i = 0; i < npoint; i++)
            member_row[i+1] = member_row[i] + member_count[i];
        size_t *member_fill = (size_t *)malloc(
            npoint * sizeof(*member_fill));
        if (member_fill == NULL) {
            free(member_row); free(member); free(dist); free(path_w);
            free(par); free(usz); free(degree); free(row); free(col);
            free(weight); free(member_count); free(radius); free(theta);
            free(outer); free(point);
            section_face_index_free(&section_index);
            shortcut_audit_free(out);
            return -1;
        }
        memcpy(member_fill, member_row, npoint * sizeof(*member_fill));
        for (size_t i = 0; i < npoint; i++) {
            int32_t root = uf_find(par, (int32_t)i);
            member[member_fill[(size_t)root]++] = (int32_t)i;
        }
        free(member_fill);

        SectionHeapNode *heap = NULL;
        size_t heap_cap = 0;
        int heap_failed = 0;
        for (size_t root = 0; root < npoint && !heap_failed; root++) {
            size_t count = member_count[root];
            if (count < P->shortcut_min_points || !outer[root])
                continue;
            out->curves_tested++;
            int boundary_clipped = 0;
            for (size_t at = member_row[root];
                 at < member_row[root+1]; at++) {
                size_t i = (size_t)member[at];
                if (degree[i] <= 1 &&
                    (P->boundary_box_set ?
                     near_boundary_box(point[i*3+0], point[i*3+1],
                                       point[i*3+2], P) :
                     near_seam(point[i*3+0], point[i*3+1], point[i*3+2],
                               P->seam_pitch, P->seam_tol))) {
                    boundary_clipped = 1;
                    break;
                }
            }
            double best_delta = 0.0, best_signed_delta = 0.0;
            double best_distance = 0.0;
            int32_t best_a = -1, best_b = -1;
            int simple = section_simple_shortcut(
                member, member_row[root], member_row[root+1], row, col,
                weight, outer, radius, theta, max_path, P->pitch,
                &best_delta, &best_signed_delta, &best_distance,
                &best_a, &best_b);
            if (simple) {
                out->simple_curves++;
            } else {
                out->fallback_curves++;
                for (size_t sat = member_row[root];
                     sat < member_row[root+1] && !heap_failed; sat++) {
                    int32_t source = member[sat];
                    for (size_t at = member_row[root];
                         at < member_row[root+1]; at++)
                        dist[(size_t)member[at]] = DBL_MAX;
                    dist[(size_t)source] = 0.0;
                    path_w[(size_t)source] = 0.0;
                    size_t heap_n = 0;
                    if (section_heap_push(&heap, &heap_n, &heap_cap,
                                          0.0, source) != 0) {
                        heap_failed = 1;
                        break;
                    }
                    SectionHeapNode hn;
                    while (section_heap_pop(heap, &heap_n, &hn)) {
                        int32_t u = hn.v;
                        if (hn.d > dist[(size_t)u]) continue;
                        double delta = fabs(path_w[(size_t)u]);
                        if (delta > best_delta) {
                            best_delta = delta;
                            best_signed_delta = path_w[(size_t)u];
                            best_distance = hn.d;
                            best_a = source;
                            best_b = u;
                        }
                        for (size_t e = row[(size_t)u];
                             e < row[(size_t)u+1]; e++) {
                            int32_t v = col[e];
                            if (!outer[(size_t)v]) continue;
                            double next = hn.d + weight[e];
                            if (next <= max_path && next < dist[(size_t)v]) {
                                dist[(size_t)v] = next;
                                path_w[(size_t)v] = path_w[(size_t)u] +
                                    PitchTable_dturns(
                                        P->ptab, radius[(size_t)u],
                                        radius[(size_t)v], P->pitch) -
                                    wrap_pmpi(theta[(size_t)v] -
                                              theta[(size_t)u]) / TWO_PI;
                                if (section_heap_push(
                                        &heap, &heap_n, &heap_cap,
                                        next, v) != 0) {
                                    heap_failed = 1;
                                    break;
                                }
                            }
                        }
                        if (heap_failed) break;
                    }
                }
            }
            if (heap_failed) break;
            if (best_delta >= P->shortcut_tol &&
                best_a >= 0 && best_b >= 0) {
                ShortcutHit h;
                memset(&h, 0, sizeof(h));
                h.delta_w = best_delta;
                h.signed_delta_w = best_signed_delta;
                h.distance = best_distance;
                h.plane = plane;
                h.curve_points = count;
                h.boundary_clipped = boundary_clipped;
                for (int axis = 0; axis < 3; axis++) {
                    h.a[axis] = point[(size_t)best_a*3+(size_t)axis];
                    h.b[axis] = point[(size_t)best_b*3+(size_t)axis];
                    h.p[axis] = 0.5 * (h.a[axis] + h.b[axis]);
                }
                if (shortcut_hit_push(out, &h) != 0) {
                    heap_failed = 1;
                    break;
                }
                if (best_delta > out->max_delta_w)
                    out->max_delta_w = best_delta;
            }
        }
        free(heap);
        free(member_row); free(member); free(dist); free(path_w);
        free(par); free(usz); free(degree); free(row); free(col);
        free(weight); free(member_count); free(radius); free(theta);
        free(outer); free(point);
        if (heap_failed) {
            section_face_index_free(&section_index);
            shortcut_audit_free(out);
            return -1;
        }
    }

    section_face_index_free(&section_index);

    if (out->nhit == 0) return 0;

    /* Join nearby section hits into 3-D tracks. Persistent tracks suppress an
     * isolated sliver while retaining a bridge wall that survives through z. */
    size_t n = out->nhit;
    int32_t *par = (int32_t *)malloc(n * sizeof(*par));
    int32_t *usz = (int32_t *)malloc(n * sizeof(*usz));
    int32_t *root_cluster = (int32_t *)malloc(n * sizeof(*root_cluster));
    if (par == NULL || usz == NULL || root_cluster == NULL) {
        free(par); free(usz); free(root_cluster);
        shortcut_audit_free(out);
        return -1;
    }
    for (size_t i = 0; i < n; i++) {
        par[i] = (int32_t)i;
        usz[i] = 1;
        root_cluster[i] = -1;
    }
    double join = P->shortcut_cluster_pitches * P->pitch;
    double join2 = join * join;
    if (P->shortcut_track_mutual) {
        /* A proximity union merges parallel bridges through same-plane and
         * transitive links.  Mutual nearest-neighbour temporal matching instead
         * gives every section hit at most one predecessor and one successor.
         * This preserves distinct bridge walls as distinct tracks, which is
         * essential when each wall needs its own separating seam. */
        int32_t *next = (int32_t *)malloc(n * sizeof(*next));
        int32_t *prev = (int32_t *)malloc(n * sizeof(*prev));
        double *next_d2 = (double *)malloc(n * sizeof(*next_d2));
        double *prev_d2 = (double *)malloc(n * sizeof(*prev_d2));
        if (next == NULL || prev == NULL || next_d2 == NULL ||
            prev_d2 == NULL) {
            free(next); free(prev); free(next_d2); free(prev_d2);
            free(par); free(usz); free(root_cluster);
            shortcut_audit_free(out);
            return -1;
        }
        for (size_t i = 0; i < n; i++) {
            next[i] = prev[i] = -1;
            next_d2[i] = prev_d2[i] = DBL_MAX;
        }
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                double dz = out->hit[j].p[0] - out->hit[i].p[0];
                if (dz > join) break; /* hits are appended in plane order */
                if (dz <= 1e-6) continue; /* never merge distinct same-plane hits */
                double dy = out->hit[j].p[1] - out->hit[i].p[1];
                double dx = out->hit[j].p[2] - out->hit[i].p[2];
                double d2 = dz*dz + dy*dy + dx*dx;
                if (d2 > join2) continue;
                if (d2 < next_d2[i]) {
                    next_d2[i] = d2;
                    next[i] = (int32_t)j;
                }
                if (d2 < prev_d2[j]) {
                    prev_d2[j] = d2;
                    prev[j] = (int32_t)i;
                }
            }
        }
        for (size_t i = 0; i < n; i++) {
            int32_t j = next[i];
            if (j >= 0 && prev[(size_t)j] == (int32_t)i)
                uf_union(par, usz, (int32_t)i, j);
        }
        free(next); free(prev); free(next_d2); free(prev_d2);
    } else {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                double dz = out->hit[j].p[0] - out->hit[i].p[0];
                if (dz > join) break; /* hits are appended in plane order */
                double dy = out->hit[j].p[1] - out->hit[i].p[1];
                double dx = out->hit[j].p[2] - out->hit[i].p[2];
                if (dz*dz + dy*dy + dx*dx <= join2)
                    uf_union(par, usz, (int32_t)i, (int32_t)j);
            }
        }
    }
    size_t ncluster = 0;
    for (size_t i = 0; i < n; i++) {
        int32_t root = uf_find(par, (int32_t)i);
        if (root_cluster[(size_t)root] < 0)
            root_cluster[(size_t)root] = (int32_t)ncluster++;
    }
    out->cluster = (ShortcutCluster *)calloc(
        ncluster, sizeof(*out->cluster));
    double *last_plane = (double *)malloc(ncluster * sizeof(*last_plane));
    if (out->cluster == NULL || last_plane == NULL) {
        free(last_plane); free(par); free(usz); free(root_cluster);
        shortcut_audit_free(out);
        return -1;
    }
    out->ncluster = ncluster;
    for (size_t c = 0; c < ncluster; c++) {
        ShortcutCluster *g = &out->cluster[c];
        g->source_id = (int32_t)c;
        g->z_min = DBL_MAX; g->z_max = -DBL_MAX;
        g->min_distance = DBL_MAX;
        last_plane[c] = -DBL_MAX;
        for (int axis = 0; axis < 3; axis++) {
            g->bbox_min[axis] = DBL_MAX;
            g->bbox_max[axis] = -DBL_MAX;
        }
    }
    for (size_t i = 0; i < n; i++) {
        int32_t root = uf_find(par, (int32_t)i);
        int32_t ci = root_cluster[(size_t)root];
        ShortcutCluster *g = &out->cluster[(size_t)ci];
        ShortcutHit *h = &out->hit[i];
        h->cluster_id = ci;
        g->hits++;
        if (h->boundary_clipped) g->boundary_hits++;
        if (fabs(h->plane - last_plane[(size_t)ci]) > 1e-6) {
            g->planes++;
            last_plane[(size_t)ci] = h->plane;
        }
        if (h->plane < g->z_min) g->z_min = h->plane;
        if (h->plane > g->z_max) g->z_max = h->plane;
        if (h->delta_w > g->max_delta_w) g->max_delta_w = h->delta_w;
        if (h->distance < g->min_distance) g->min_distance = h->distance;
        g->excess += h->delta_w - P->shortcut_tol;
        for (int axis = 0; axis < 3; axis++) {
            g->centroid[axis] += h->p[axis];
            double low = h->a[axis] < h->b[axis] ?
                         h->a[axis] : h->b[axis];
            double high = h->a[axis] > h->b[axis] ?
                          h->a[axis] : h->b[axis];
            if (low < g->bbox_min[axis]) g->bbox_min[axis] = low;
            if (high > g->bbox_max[axis]) g->bbox_max[axis] = high;
        }
    }
    for (size_t c = 0; c < ncluster; c++) {
        ShortcutCluster *g = &out->cluster[c];
        for (int axis = 0; axis < 3; axis++)
            g->centroid[axis] /= (double)g->hits;
        g->boundary_supported =
            (int)g->planes < P->shortcut_min_planes &&
            g->planes >= 2 && g->boundary_hits >= 2;
        g->persistent = (int)g->planes >= P->shortcut_min_planes ||
                        g->boundary_supported;
        if (g->persistent) {
            out->persistent_clusters++;
            out->persistent_hits += g->hits;
            out->persistent_excess += g->excess;
            if (g->max_delta_w > out->persistent_max_delta_w)
                out->persistent_max_delta_w = g->max_delta_w;
        }
    }
    for (size_t i = 0; i < n; i++) {
        int32_t root = uf_find(par, (int32_t)i);
        int32_t ci = root_cluster[(size_t)root];
        if (out->cluster[(size_t)ci].persistent)
            out->hit[i].persistent = 1;
    }
    qsort(out->cluster, out->ncluster, sizeof(*out->cluster),
          shortcut_cluster_cmp);
    /* Clusters are severity-sorted for reporting. Hits retain the pre-sort
     * identity until this point, so remap them before repair consumes the
     * track id. root_cluster has n slots and ncluster <= n. */
    for (size_t c = 0; c < ncluster; c++)
        root_cluster[(size_t)out->cluster[c].source_id] = (int32_t)c;
    for (size_t i = 0; i < n; i++)
        out->hit[i].cluster_id =
            root_cluster[(size_t)out->hit[i].cluster_id];
    free(last_plane);
    free(par); free(usz); free(root_cluster);
    return 0;
}

/* ---------- supported shortcut seam repair ---------- */
typedef struct {
    size_t components, small_components, biggest_faces;
} MeshComponentSummary;

typedef struct {
    uint64_t key;
    int32_t head, count, track;
    int selected;
} RepairEdgeSlot;

typedef struct {
    int32_t face, corner_lo, corner_hi, next;
} RepairEdgeOccurrence;

typedef struct {
    size_t hits;
    double z_mean, mid_mean, intercept, slope, half_mean;
    int active;
} RepairTrack;

typedef struct {
    size_t labelled_low_faces, labelled_high_faces;
    size_t candidate_edges, selected_chains, cut_edges;
    size_t selected_cut_components;
    size_t selected_chain_endpoints, selected_chain_boundary_endpoints;
    size_t selected_chain_branch_vertices, selected_chain_branch_excess;
    size_t cut_rank;
    size_t vertices_before, vertices_after, faces;
    MeshComponentSummary before, after;
    ShortcutAudit post;
    double phase_offset;
    int auto_trials;
    int auto_decisive;
    int auto_suppressed;
    int accepted, step_accepted, output_written;
    char reason[160];
} ShortcutRepair;

static int repair_trial_monotone(const ShortcutRepair *r,
                                 const ShortcutAudit *before,
                                 const Params *P, size_t nfaces);
static size_t positive_size_delta(size_t after, size_t before);

static double point_phase(const double p[3], const Params *P)
{
    double uy, ux;
    axis_table_eval(P->axis, p[0], P->umb_y, P->umb_x, &uy, &ux);
    double dy = p[1] - uy;
    double dx = p[2] - ux;
    return rturns(P, hypot(dy, dx)) - atan2(dy, dx) / TWO_PI;
}

static double vertex_phase(const FVec *V, size_t i, const Params *P)
{
    double p[3] = {
        (double)V->v[i*3+0], (double)V->v[i*3+1],
        (double)V->v[i*3+2]
    };
    return point_phase(p, P);
}

static int mesh_component_summary(const FVec *V, const IVec *F,
                                  size_t small_limit,
                                  MeshComponentSummary *out)
{
    memset(out, 0, sizeof(*out));
    int32_t *par = (int32_t *)malloc(V->n * sizeof(*par));
    int32_t *usz = (int32_t *)malloc(V->n * sizeof(*usz));
    size_t *faces = (size_t *)calloc(V->n, sizeof(*faces));
    if (par == NULL || usz == NULL || faces == NULL) {
        free(par); free(usz); free(faces);
        return -1;
    }
    for (size_t i = 0; i < V->n; i++) {
        par[i] = (int32_t)i;
        usz[i] = 1;
    }
    for (size_t f = 0; f < F->n; f++) {
        int32_t a = F->f[f*3+0], b = F->f[f*3+1], c = F->f[f*3+2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= V->n || (size_t)b >= V->n ||
            (size_t)c >= V->n) {
            free(par); free(usz); free(faces);
            return -1;
        }
        uf_union(par, usz, a, b);
        uf_union(par, usz, b, c);
    }
    for (size_t f = 0; f < F->n; f++)
        faces[(size_t)uf_find(par, F->f[f*3+0])]++;
    for (size_t i = 0; i < V->n; i++) {
        if (faces[i] == 0) continue;
        out->components++;
        if (faces[i] < small_limit) out->small_components++;
        if (faces[i] > out->biggest_faces) out->biggest_faces = faces[i];
    }
    free(par); free(usz); free(faces);
    return 0;
}

static int shortcut_seam_repair(const FVec *V, const IVec *F,
                                const Params *P,
                                const ShortcutAudit *audit,
                                ShortcutRepair *out)
{
    RepairTrack *model = NULL;
    int8_t *face_label = NULL;
    int32_t *face_track = NULL;
    RepairEdgeSlot *edge = NULL;
    RepairEdgeOccurrence *occ = NULL;
    size_t *candidate_slot = NULL;
    int32_t *cpar = NULL, *csz = NULL;
    uint64_t *tvkey = NULL;
    int32_t *tvvalue = NULL;
    uint16_t *tvdegree = NULL;
    size_t *chain_edges = NULL, *chain_support = NULL;
    size_t *chain_endpoints = NULL, *chain_boundary_endpoints = NULL;
    size_t *chain_branches = NULL;
    uint8_t *chain_selected = NULL, *boundary_vertex = NULL;
    double *chain_zmin = NULL, *chain_zmax = NULL;
    int32_t *best_root = NULL;
    int32_t *corner_par = NULL, *corner_size = NULL;
    int32_t *first_corner = NULL, *root_new = NULL;
    int32_t *cut_par = NULL, *cut_size = NULL;
    uint32_t *cut_degree = NULL;
    uint8_t *cut_vertex = NULL;
    FVec repaired_v = {0};
    IVec repaired_f = {0};
    int rc = -1;

    memset(out, 0, sizeof(*out));
    out->phase_offset = P->repair_phase_offset;
    out->vertices_before = V->n;
    out->faces = F->n;
    if (mesh_component_summary(V, F, P->min_faces, &out->before) != 0) {
        snprintf(out->reason, sizeof(out->reason), "invalid input mesh");
        goto cleanup;
    }
    if (audit->persistent_clusters == 0) {
        out->vertices_after = V->n;
        out->after = out->before;
        out->accepted = 1;
        snprintf(out->reason, sizeof(out->reason),
                 "already clean; copied without repair");
        if (P->repair_path != NULL) {
            if (write_obj(P->repair_path, V, F) != 0) goto cleanup;
            out->output_written = 1;
        }
        return 0;
    }

    size_t ntrack = audit->ncluster;
    model = (RepairTrack *)calloc(ntrack ? ntrack : 1, sizeof(*model));
    face_label = (int8_t *)calloc(F->n, sizeof(*face_label));
    face_track = (int32_t *)malloc(F->n * sizeof(*face_track));
    if (model == NULL || face_label == NULL || face_track == NULL) goto oom;
    for (size_t f = 0; f < F->n; f++) face_track[f] = -1;

    /* Fit one smooth phase-midpoint line w_mid(z) per persistent track. This
     * suppresses the many broken isolines produced by fitting every section
     * independently while allowing slow axial drift of the real sheet. */
    for (size_t i = 0; i < audit->nhit; i++) {
        const ShortcutHit *h = &audit->hit[i];
        if (!h->persistent || h->cluster_id < 0 ||
            (size_t)h->cluster_id >= ntrack) continue;
        RepairTrack *m = &model[(size_t)h->cluster_id];
        double mid = point_phase(h->a, P) + 0.5 * h->signed_delta_w;
        if (m->hits > 0) {
            double ref = m->mid_mean / (double)m->hits;
            mid += floor(ref - mid + 0.5);
        }
        m->hits++;
        m->z_mean += h->p[0];
        m->mid_mean += mid;
        m->half_mean += 0.5 * h->delta_w;
        m->active = 1;
    }
    for (size_t t = 0; t < ntrack; t++) {
        if (!model[t].active) continue;
        model[t].z_mean /= (double)model[t].hits;
        model[t].mid_mean /= (double)model[t].hits;
        model[t].half_mean /= (double)model[t].hits;
    }
    for (size_t i = 0; i < audit->nhit; i++) {
        const ShortcutHit *h = &audit->hit[i];
        if (!h->persistent || h->cluster_id < 0 ||
            (size_t)h->cluster_id >= ntrack) continue;
        RepairTrack *m = &model[(size_t)h->cluster_id];
        double mid = point_phase(h->a, P) + 0.5 * h->signed_delta_w;
        mid += floor(m->mid_mean - mid + 0.5);
        double dz = h->p[0] - m->z_mean;
        m->slope += dz * (mid - m->mid_mean);
        m->intercept += dz * dz;
    }
    for (size_t t = 0; t < ntrack; t++) {
        RepairTrack *m = &model[t];
        if (!m->active) continue;
        m->slope = m->intercept > 0.0 ? m->slope / m->intercept : 0.0;
        m->intercept = m->mid_mean - m->slope * m->z_mean;
    }

    double radius = P->repair_radius_pitches * P->pitch;
    double radius2 = radius * radius;
    for (size_t f = 0; f < F->n; f++) {
        int32_t tri[3] = {
            F->f[f*3+0], F->f[f*3+1], F->f[f*3+2]
        };
        double centroid[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; k++)
            for (int axis = 0; axis < 3; axis++)
                centroid[axis] +=
                    (double)V->v[(size_t)tri[k]*3+(size_t)axis] / 3.0;
        size_t best = SIZE_MAX;
        double best_d2 = DBL_MAX;
        for (size_t i = 0; i < audit->nhit; i++) {
            const ShortcutHit *h = &audit->hit[i];
            if (!h->persistent) continue;
            double dz = centroid[0] - h->p[0];
            double dy = centroid[1] - h->p[1];
            double dx = centroid[2] - h->p[2];
            double d2 = dz*dz + dy*dy + dx*dx;
            if (d2 < best_d2) { best_d2 = d2; best = i; }
        }
        if (best == SIZE_MAX || best_d2 > radius2) continue;
        int32_t track = audit->hit[best].cluster_id;
        if (track < 0 || (size_t)track >= ntrack ||
            !model[(size_t)track].active) continue;
        RepairTrack *m = &model[(size_t)track];
        double mid = m->intercept + m->slope * centroid[0];
        double face_phase = 0.0;
        for (int k = 0; k < 3; k++) {
            double w = vertex_phase(V, (size_t)tri[k], P);
            w += floor(mid - w + 0.5);
            face_phase += w / 3.0;
        }
        if (fabs(face_phase - mid) >
            m->half_mean + P->repair_band_margin) continue;
        face_label[f] = face_phase >= mid + P->repair_phase_offset ? 1 : -1;
        face_track[f] = track;
        if (face_label[f] > 0) out->labelled_high_faces++;
        else out->labelled_low_faces++;
    }

    size_t edge_slots = 1;
    while (edge_slots < F->n * 8) edge_slots <<= 1;
    edge = (RepairEdgeSlot *)calloc(edge_slots, sizeof(*edge));
    occ = (RepairEdgeOccurrence *)malloc(F->n * 3 * sizeof(*occ));
    if (edge == NULL || occ == NULL) goto oom;
    for (size_t s = 0; s < edge_slots; s++) {
        edge[s].key = UINT64_MAX;
        edge[s].head = -1;
        edge[s].track = -1;
    }
    size_t edge_mask = edge_slots - 1, nocc = 0;
    for (size_t f = 0; f < F->n; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = F->f[f*3+(size_t)k];
            int32_t b = F->f[f*3+(size_t)((k+1)%3)];
            int32_t ca = (int32_t)(f*3+(size_t)k);
            int32_t cb = (int32_t)(f*3+(size_t)((k+1)%3));
            if (a < 0 || b < 0 || (size_t)a >= V->n || (size_t)b >= V->n)
                goto invalid;
            if (a > b) {
                int32_t tmp = a; a = b; b = tmp;
                tmp = ca; ca = cb; cb = tmp;
            }
            uint64_t key = ((uint64_t)(uint32_t)a << 32) |
                           (uint64_t)(uint32_t)b;
            size_t s = hash_u64(key, edge_mask);
            while (edge[s].key != UINT64_MAX && edge[s].key != key)
                s = (s + 1) & edge_mask;
            if (edge[s].key == UINT64_MAX) edge[s].key = key;
            occ[nocc].face = (int32_t)f;
            occ[nocc].corner_lo = ca;
            occ[nocc].corner_hi = cb;
            occ[nocc].next = edge[s].head;
            edge[s].head = (int32_t)nocc;
            edge[s].count++;
            nocc++;
        }
    }

    candidate_slot = (size_t *)malloc(F->n * 3 * sizeof(*candidate_slot));
    if (candidate_slot == NULL) goto oom;
    size_t ncandidate = 0;
    for (size_t s = 0; s < edge_slots; s++) {
        if (edge[s].key == UINT64_MAX || edge[s].count != 2) continue;
        int32_t o0 = edge[s].head, o1 = occ[(size_t)o0].next;
        int32_t f0 = occ[(size_t)o0].face, f1 = occ[(size_t)o1].face;
        if (face_label[(size_t)f0] * face_label[(size_t)f1] < 0 &&
            face_track[(size_t)f0] >= 0 &&
            face_track[(size_t)f0] == face_track[(size_t)f1]) {
            edge[s].track = face_track[(size_t)f0];
            candidate_slot[ncandidate++] = s;
        }
    }
    out->candidate_edges = ncandidate;
    if (ncandidate == 0) {
        snprintf(out->reason, sizeof(out->reason),
                 "no supported phase-isoline edge candidates");
        rc = 0;
        goto cleanup;
    }

    cpar = (int32_t *)malloc(ncandidate * sizeof(*cpar));
    csz = (int32_t *)malloc(ncandidate * sizeof(*csz));
    size_t tvslots = 1;
    while (tvslots < ncandidate * 8) tvslots <<= 1;
    tvkey = (uint64_t *)malloc(tvslots * sizeof(*tvkey));
    tvvalue = (int32_t *)malloc(tvslots * sizeof(*tvvalue));
    tvdegree = (uint16_t *)calloc(tvslots, sizeof(*tvdegree));
    chain_edges = (size_t *)calloc(ncandidate, sizeof(*chain_edges));
    chain_support = (size_t *)calloc(ncandidate, sizeof(*chain_support));
    chain_endpoints = (size_t *)calloc(ncandidate, sizeof(*chain_endpoints));
    chain_boundary_endpoints =
        (size_t *)calloc(ncandidate, sizeof(*chain_boundary_endpoints));
    chain_branches = (size_t *)calloc(ncandidate, sizeof(*chain_branches));
    chain_selected = (uint8_t *)calloc(ncandidate, 1);
    boundary_vertex = (uint8_t *)calloc(V->n, 1);
    chain_zmin = (double *)malloc(ncandidate * sizeof(*chain_zmin));
    chain_zmax = (double *)malloc(ncandidate * sizeof(*chain_zmax));
    best_root = (int32_t *)malloc((ntrack ? ntrack : 1) * sizeof(*best_root));
    if (cpar == NULL || csz == NULL || tvkey == NULL || tvvalue == NULL ||
        tvdegree == NULL || chain_edges == NULL || chain_support == NULL ||
        chain_endpoints == NULL || chain_boundary_endpoints == NULL ||
        chain_branches == NULL || chain_selected == NULL ||
        boundary_vertex == NULL || chain_zmin == NULL || chain_zmax == NULL ||
        best_root == NULL) goto oom;
    for (size_t s = 0; s < edge_slots; s++) {
        if (edge[s].key == UINT64_MAX || edge[s].count != 1) continue;
        boundary_vertex[(size_t)(uint32_t)(edge[s].key >> 32)] = 1;
        boundary_vertex[(size_t)(uint32_t)edge[s].key] = 1;
    }
    for (size_t i = 0; i < ncandidate; i++) {
        cpar[i] = (int32_t)i; csz[i] = 1;
        chain_zmin[i] = DBL_MAX; chain_zmax[i] = -DBL_MAX;
    }
    for (size_t i = 0; i < tvslots; i++) tvkey[i] = UINT64_MAX;
    size_t tvmask = tvslots - 1;
    for (size_t i = 0; i < ncandidate; i++) {
        RepairEdgeSlot *e = &edge[candidate_slot[i]];
        uint32_t verts[2] = {
            (uint32_t)(e->key >> 32), (uint32_t)e->key
        };
        for (int k = 0; k < 2; k++) {
            uint64_t key = ((uint64_t)(uint32_t)e->track << 32) | verts[k];
            size_t s = hash_u64(key, tvmask);
            while (tvkey[s] != UINT64_MAX && tvkey[s] != key)
                s = (s + 1) & tvmask;
            if (tvkey[s] == UINT64_MAX) {
                tvkey[s] = key;
                tvvalue[s] = (int32_t)i;
                tvdegree[s] = 1;
            } else {
                if (tvdegree[s] < UINT16_MAX) tvdegree[s]++;
                uf_union(cpar, csz, (int32_t)i, tvvalue[s]);
            }
        }
    }
    for (size_t i = 0; i < ncandidate; i++) {
        int32_t root = uf_find(cpar, (int32_t)i);
        RepairEdgeSlot *e = &edge[candidate_slot[i]];
        size_t a = (size_t)(uint32_t)(e->key >> 32);
        size_t b = (size_t)(uint32_t)e->key;
        double za = (double)V->v[a*3+0], zb = (double)V->v[b*3+0];
        chain_edges[(size_t)root]++;
        if (za < chain_zmin[(size_t)root]) chain_zmin[(size_t)root] = za;
        if (zb < chain_zmin[(size_t)root]) chain_zmin[(size_t)root] = zb;
        if (za > chain_zmax[(size_t)root]) chain_zmax[(size_t)root] = za;
        if (zb > chain_zmax[(size_t)root]) chain_zmax[(size_t)root] = zb;
    }
    for (size_t s = 0; s < tvslots; s++) {
        if (tvkey[s] == UINT64_MAX) continue;
        int32_t root = uf_find(cpar, tvvalue[s]);
        size_t ri = (size_t)root;
        if (tvdegree[s] == 1) {
            size_t vertex = (size_t)(uint32_t)tvkey[s];
            chain_endpoints[ri]++;
            if (boundary_vertex[vertex]) chain_boundary_endpoints[ri]++;
        } else if (tvdegree[s] > 2) {
            chain_branches[ri]++;
        }
    }
    double support_r = P->repair_support_radius_pitches * P->pitch;
    double support_r2 = support_r * support_r;
    for (size_t hidx = 0; hidx < audit->nhit; hidx++) {
        const ShortcutHit *h = &audit->hit[hidx];
        if (!h->persistent) continue;
        size_t best = SIZE_MAX;
        double best_d2 = DBL_MAX;
        for (size_t i = 0; i < ncandidate; i++) {
            RepairEdgeSlot *e = &edge[candidate_slot[i]];
            if (e->track != h->cluster_id) continue;
            size_t a = (size_t)(uint32_t)(e->key >> 32);
            size_t b = (size_t)(uint32_t)e->key;
            double dz = 0.5*((double)V->v[a*3+0] + V->v[b*3+0]) - h->p[0];
            double dy = 0.5*((double)V->v[a*3+1] + V->v[b*3+1]) - h->p[1];
            double dx = 0.5*((double)V->v[a*3+2] + V->v[b*3+2]) - h->p[2];
            double d2 = dz*dz + dy*dy + dx*dx;
            if (d2 < best_d2) { best_d2 = d2; best = i; }
        }
        if (best != SIZE_MAX && best_d2 <= support_r2) {
            int32_t root = uf_find(cpar, (int32_t)best);
            chain_support[(size_t)root]++;
        }
    }
    for (size_t t = 0; t < ntrack; t++) best_root[t] = -1;
    for (size_t i = 0; i < ncandidate; i++) {
        int32_t root = uf_find(cpar, (int32_t)i);
        if ((size_t)root != i || chain_edges[i] == 0) continue;
        int32_t track = edge[candidate_slot[i]].track;
        int32_t old = best_root[(size_t)track];
        double span = chain_zmax[i] - chain_zmin[i];
        double old_span = old >= 0 ?
            chain_zmax[(size_t)old] - chain_zmin[(size_t)old] : -1.0;
        if (old < 0 || chain_support[i] > chain_support[(size_t)old] ||
            (chain_support[i] == chain_support[(size_t)old] &&
              (span > old_span + 1e-9 ||
               (fabs(span - old_span) <= 1e-9 &&
                chain_edges[i] > chain_edges[(size_t)old]))))
            best_root[(size_t)track] = root;
    }
    if (P->repair_cover_all_supported) {
        for (size_t i = 0; i < ncandidate; i++) {
            int32_t root = uf_find(cpar, (int32_t)i);
            if ((size_t)root == i &&
                chain_support[i] >= P->repair_min_chain_support)
                out->selected_chains++;
        }
    } else {
        for (size_t t = 0; t < ntrack; t++)
            if (best_root[t] >= 0 && chain_support[(size_t)best_root[t]] > 0)
                out->selected_chains++;
    }
    for (size_t i = 0; i < ncandidate; i++) {
        RepairEdgeSlot *e = &edge[candidate_slot[i]];
        int32_t root = uf_find(cpar, (int32_t)i);
        if (chain_support[(size_t)root] >=
                (P->repair_cover_all_supported ?
                 P->repair_min_chain_support : 1) &&
            (P->repair_cover_all_supported ||
             best_root[(size_t)e->track] == root)) {
            e->selected = 1;
            chain_selected[(size_t)root] = 1;
            out->cut_edges++;
        }
    }
    /* The selection graph above is keyed by (track, vertex), intentionally
     * keeping nearby temporal tracks independent.  The TOPOLOGICAL cut graph
     * must instead be assembled globally: isolines from different tracks can
     * meet at a mesh vertex.  Treating them as disjoint hid those cross-track
     * saddles and badly underestimated separation rank on dense blocks. */
    cut_par = (int32_t *)malloc(V->n * sizeof(*cut_par));
    cut_size = (int32_t *)malloc(V->n * sizeof(*cut_size));
    cut_degree = (uint32_t *)calloc(V->n, sizeof(*cut_degree));
    if (cut_par == NULL || cut_size == NULL || cut_degree == NULL) goto oom;
    for (size_t v = 0; v < V->n; v++) {
        cut_par[v] = (int32_t)v;
        cut_size[v] = 1;
    }
    for (size_t i = 0; i < ncandidate; i++) {
        RepairEdgeSlot *e = &edge[candidate_slot[i]];
        if (!e->selected) continue;
        int32_t a = (int32_t)(uint32_t)(e->key >> 32);
        int32_t b = (int32_t)(uint32_t)e->key;
        if (cut_degree[(size_t)a] != UINT32_MAX) cut_degree[(size_t)a]++;
        if (cut_degree[(size_t)b] != UINT32_MAX) cut_degree[(size_t)b]++;
        uf_union(cut_par, cut_size, a, b);
    }
    for (size_t v = 0; v < V->n; v++) {
        uint32_t degree = cut_degree[v];
        if (degree == 0) continue;
        if (uf_find(cut_par, (int32_t)v) == (int32_t)v)
            out->selected_cut_components++;
        if (degree == 1) {
            out->selected_chain_endpoints++;
            if (boundary_vertex[v]) out->selected_chain_boundary_endpoints++;
        } else if (degree > 2) {
            out->selected_chain_branch_vertices++;
            out->selected_chain_branch_excess += (size_t)degree - 2;
        }
    }
    /* Conservative separation capacity of an embedded graph. A cycle carries
     * one unit; branches and open endpoints add local complement sectors. A
     * boundary endpoint gets one additional unit: on these cropped, already
     * disconnected sheets its two incident boundary fans can separate
     * independently. Thus a boundary-to-boundary arc has capacity five. This
     * bound is empirical but topologically motivated, and remains far below the
     * raw edge count used by an indiscriminate fragmentation allowance. */
    out->cut_rank = out->selected_cut_components +
                    out->selected_chain_branch_excess +
                    out->selected_chain_endpoints +
                    out->selected_chain_boundary_endpoints;
    if (P->repair_chain_csv != NULL) {
        FILE *csv = fopen(P->repair_chain_csv, "wb");
        if (csv == NULL) {
            snprintf(out->reason, sizeof(out->reason),
                     "cannot write repair chain csv");
            goto cleanup;
        }
        fprintf(csv,
                "track,root,edges,support,z_min,z_max,z_span,endpoints,"
                "boundary_endpoints,branch_vertices,selected\n");
        for (size_t i = 0; i < ncandidate; i++) {
            int32_t root = uf_find(cpar, (int32_t)i);
            if ((size_t)root != i || chain_edges[i] == 0) continue;
            int32_t track = edge[candidate_slot[i]].track;
            fprintf(csv,
                    "%d,%zu,%zu,%zu,%.9g,%.9g,%.9g,%zu,%zu,%zu,%d\n",
                    track, i, chain_edges[i], chain_support[i],
                    chain_zmin[i], chain_zmax[i],
                    chain_zmax[i] - chain_zmin[i], chain_endpoints[i],
                    chain_boundary_endpoints[i], chain_branches[i],
                    chain_selected[i] ? 1 : 0);
        }
        if (fclose(csv) != 0) goto cleanup;
    }
    if (out->cut_edges == 0) {
        snprintf(out->reason, sizeof(out->reason),
                 "no shortcut-supported continuous isoline chain");
        rc = 0;
        goto cleanup;
    }

    size_t ncorner = F->n * 3;
    corner_par = (int32_t *)malloc(ncorner * sizeof(*corner_par));
    corner_size = (int32_t *)malloc(ncorner * sizeof(*corner_size));
    first_corner = (int32_t *)malloc(V->n * sizeof(*first_corner));
    root_new = (int32_t *)malloc(ncorner * sizeof(*root_new));
    cut_vertex = (uint8_t *)calloc(V->n, 1);
    if (corner_par == NULL || corner_size == NULL || first_corner == NULL ||
        root_new == NULL || cut_vertex == NULL) goto oom;
    for (size_t i = 0; i < ncorner; i++) {
        corner_par[i] = (int32_t)i;
        corner_size[i] = 1;
        root_new[i] = -1;
    }
    for (size_t i = 0; i < V->n; i++) first_corner[i] = -1;
    for (size_t s = 0; s < edge_slots; s++) {
        if (edge[s].key == UINT64_MAX) continue;
        size_t a = (size_t)(uint32_t)(edge[s].key >> 32);
        size_t b = (size_t)(uint32_t)edge[s].key;
        if (edge[s].selected) {
            cut_vertex[a] = 1; cut_vertex[b] = 1;
            continue;
        }
        int32_t first = edge[s].head;
        if (first < 0) continue;
        for (int32_t oi = occ[(size_t)first].next; oi >= 0;
             oi = occ[(size_t)oi].next) {
            uf_union(corner_par, corner_size,
                     occ[(size_t)first].corner_lo,
                     occ[(size_t)oi].corner_lo);
            uf_union(corner_par, corner_size,
                     occ[(size_t)first].corner_hi,
                     occ[(size_t)oi].corner_hi);
        }
    }
    for (size_t f = 0; f < F->n; f++) {
        for (int k = 0; k < 3; k++) {
            size_t corner = f*3+(size_t)k;
            int32_t v = F->f[corner];
            if (cut_vertex[(size_t)v]) continue;
            if (first_corner[(size_t)v] < 0)
                first_corner[(size_t)v] = (int32_t)corner;
            else
                uf_union(corner_par, corner_size,
                         first_corner[(size_t)v], (int32_t)corner);
        }
    }
    for (size_t f = 0; f < F->n; f++) {
        int32_t tri[3];
        for (int k = 0; k < 3; k++) {
            size_t corner = f*3+(size_t)k;
            int32_t root = uf_find(corner_par, (int32_t)corner);
            if (root_new[(size_t)root] < 0) {
                int32_t old = F->f[corner];
                if (fv_push(&repaired_v, V->v[(size_t)old*3+0],
                            V->v[(size_t)old*3+1],
                            V->v[(size_t)old*3+2]) != 0) goto oom;
                root_new[(size_t)root] = (int32_t)(repaired_v.n - 1);
            }
            tri[k] = root_new[(size_t)root];
        }
        if (iv_push(&repaired_f, tri[0], tri[1], tri[2]) != 0) goto oom;
    }
    out->vertices_after = repaired_v.n;
    if (mesh_component_summary(&repaired_v, &repaired_f, P->min_faces,
                               &out->after) != 0) goto invalid;
    if (shortcut_audit(&repaired_v, &repaired_f, P, &out->post) != 0) goto oom;

    size_t cut_limit = (size_t)ceil(P->repair_max_cut_fraction *
                                    (double)F->n);
    if (out->cut_edges > cut_limit) {
        snprintf(out->reason, sizeof(out->reason),
                 "cut budget exceeded: %zu > %zu edges", out->cut_edges,
                 cut_limit);
    } else if (out->selected_chain_branch_vertices != 0) {
        snprintf(out->reason, sizeof(out->reason),
                 "branched cut graph: %zu branch vertex/vertices; "
                 "choose a regular phase offset",
                 out->selected_chain_branch_vertices);
    } else if (out->after.small_components > out->before.small_components) {
        snprintf(out->reason, sizeof(out->reason),
                 "new small components forbidden: %zu -> %zu below %zu faces",
                 out->before.small_components, out->after.small_components,
                 P->min_faces);
    } else if (positive_size_delta(out->after.components,
                                   out->before.components) >
               out->selected_chains) {
        snprintf(out->reason, sizeof(out->reason),
                 "component growth exceeds chain count: %zu -> %zu "
                 "from %zu chain(s)", out->before.components,
                 out->after.components, out->selected_chains);
    } else if (out->after.components >
               out->before.components + out->cut_rank) {
        snprintf(out->reason, sizeof(out->reason),
                 "cut-graph rank exceeded: components %zu -> %zu (rank %zu)",
                 out->before.components, out->after.components,
                 out->cut_rank);
    } else if (out->post.persistent_clusters == 0) {
        out->accepted = 1;
        snprintf(out->reason, sizeof(out->reason),
                 "all shortcut tracks cleared; faces and coordinates preserved");
    } else if (P->repair_monotone_step &&
               repair_trial_monotone(out, audit, P, F->n)) {
        out->step_accepted = 1;
        snprintf(out->reason, sizeof(out->reason),
                 "bounded monotone step: tracks %zu->%zu, hits %zu->%zu, "
                 "excess %.3f->%.3f",
                 audit->persistent_clusters, out->post.persistent_clusters,
                 audit->persistent_hits, out->post.persistent_hits,
                 audit->persistent_excess, out->post.persistent_excess);
    } else {
        snprintf(out->reason, sizeof(out->reason),
                 "post-repair audit still has %zu persistent track(s)",
                 out->post.persistent_clusters);
    }
    if ((out->accepted || out->step_accepted) && P->repair_path != NULL) {
        if (write_obj(P->repair_path, &repaired_v, &repaired_f) != 0)
            goto cleanup;
        out->output_written = 1;
    }
    /* Research-only continuation lane: emit a bounded, topology-preserving
     * candidate even when one-pass acceptance fails, so repeated supported
     * cuts can be evaluated empirically. The default remains fail-closed. */
    if (!out->accepted && P->repair_emit_rejected && P->repair_path != NULL) {
        if (write_obj(P->repair_path, &repaired_v, &repaired_f) != 0)
            goto cleanup;
        out->output_written = 1;
    }
    rc = 0;
    goto cleanup;

invalid:
    snprintf(out->reason, sizeof(out->reason), "invalid mesh topology");
    goto cleanup;
oom:
    snprintf(out->reason, sizeof(out->reason), "out of memory");
cleanup:
    free(model); free(face_label); free(face_track); free(edge); free(occ);
    free(candidate_slot); free(cpar); free(csz); free(tvkey); free(tvvalue);
    free(tvdegree); free(chain_edges); free(chain_support);
    free(chain_endpoints); free(chain_boundary_endpoints);
    free(chain_branches); free(chain_selected); free(boundary_vertex);
    free(chain_zmin); free(chain_zmax);
    free(best_root); free(cut_par); free(cut_size); free(cut_degree);
    free(corner_par); free(corner_size); free(first_corner);
    free(root_new); free(cut_vertex); free(repaired_v.v); free(repaired_f.f);
    return rc;
}

static size_t positive_size_delta(size_t after, size_t before)
{
    return after > before ? after - before : 0;
}

/* A trial is admissible when it is a bounded topological cut and either
 * clears the shortcut field or strictly lowers its integrated excess.  A
 * zero-cut early exit has no post-audit, so it must not masquerade as clean. */
static int repair_trial_admissible(const ShortcutRepair *r,
                                   const ShortcutAudit *before,
                                   const Params *P, size_t nfaces)
{
    size_t cut_limit = (size_t)ceil(P->repair_max_cut_fraction *
                                    (double)nfaces);
    double eps = 1e-9 * fmax(1.0, before->persistent_excess);
    if (before->persistent_clusters != 0 && r->cut_edges == 0) return 0;
    if (r->cut_edges > cut_limit) return 0;
    if (r->selected_chain_branch_vertices != 0) return 0;
    if (r->after.small_components > r->before.small_components) return 0;
    if (positive_size_delta(r->after.components, r->before.components) >
        r->selected_chains) return 0;
    if (positive_size_delta(r->after.components, r->before.components) >
        r->cut_rank) return 0;
    if (r->post.persistent_clusters != 0 &&
        !(r->post.persistent_excess < before->persistent_excess - eps))
        return 0;
    return 1;
}

/* A partial cut may be materialized only when every headline shortcut metric
 * is monotone, the worst local event does not grow, and integrated excess
 * strictly falls. Structural cut/fragmentation bounds remain mandatory. */
static int repair_trial_monotone(const ShortcutRepair *r,
                                 const ShortcutAudit *before,
                                 const Params *P, size_t nfaces)
{
    double eps = 1e-9 * fmax(1.0, before->persistent_excess);
    double max_eps = 1e-9 * fmax(1.0, before->persistent_max_delta_w);
    if (!repair_trial_admissible(r, before, P, nfaces)) return 0;
    if (r->post.persistent_clusters > before->persistent_clusters) return 0;
    if (r->post.persistent_hits > before->persistent_hits) return 0;
    if (r->post.persistent_max_delta_w >
        before->persistent_max_delta_w + max_eps) return 0;
    if (r->after.small_components > r->before.small_components) return 0;
    return r->post.persistent_excess < before->persistent_excess - eps;
}

static int repair_trial_selectable(const ShortcutRepair *r,
                                   const ShortcutAudit *before,
                                   const Params *P, size_t nfaces)
{
    if (!repair_trial_admissible(r, before, P, nfaces)) return 0;
    return !P->repair_monotone_step ||
           repair_trial_monotone(r, before, P, nfaces);
}

static int repair_auto_search_allowed(const Params *P,
                                      const ShortcutAudit *before)
{
    return P->repair_max_auto_tracks == 0 ||
           before->persistent_clusters <= P->repair_max_auto_tracks;
}

/* Branch vertices are never publishable: they create several complement
 * sectors at one mesh vertex and were the source of visually ugly component
 * confetti in large-region trials. The caller samples nearby regular values
 * until it finds a simple embedded cut graph. A branch-free base may skip that
 * sweep when it removes at least 60% of tracks, hits, and integrated excess
 * with half of every structural budget still in reserve. */
static int repair_trial_decisive(const ShortcutRepair *r,
                                 const ShortcutAudit *before,
                                 const Params *P, size_t nfaces)
{
    size_t cut_limit;
    if (!repair_trial_selectable(r, before, P, nfaces)) return 0;
    if (before->persistent_clusters == 0 ||
        r->post.persistent_clusters == 0) return 1;
    if (before->persistent_hits == 0 ||
        !(before->persistent_excess > 0.0)) return 0;
    if ((double)r->post.persistent_clusters >
            0.40 * (double)before->persistent_clusters ||
        (double)r->post.persistent_hits >
            0.40 * (double)before->persistent_hits ||
        r->post.persistent_excess > 0.40 * before->persistent_excess)
        return 0;
    cut_limit = (size_t)ceil(P->repair_max_cut_fraction * (double)nfaces);
    if ((double)r->cut_edges > 0.5 * (double)cut_limit ||
        (double)positive_size_delta(r->after.components,
                                    r->before.components) >
            0.5 * (double)r->cut_rank ||
        (double)positive_size_delta(r->after.small_components,
                                    r->before.small_components) >
            0.5 * (double)r->cut_rank)
        return 0;
    return 1;
}

static int repair_needs_offset_search(const ShortcutRepair *r,
                                      const ShortcutAudit *before,
                                      const Params *P, size_t nfaces)
{
    if (!repair_trial_selectable(r, before, P, nfaces)) return 1;
    /* A safe but weak base still gets a local regular-value sweep. Only a
     * decisive base can skip it; non-selectable trials, including every
     * branch, always request the sweep. */
    return !repair_trial_decisive(r, before, P, nfaces);
}

static size_t repair_effective_tracks(const ShortcutRepair *r,
                                      const ShortcutAudit *before)
{
    return (before->persistent_clusters != 0 && r->cut_edges == 0) ?
           before->persistent_clusters : r->post.persistent_clusters;
}

static size_t repair_effective_hits(const ShortcutRepair *r,
                                    const ShortcutAudit *before)
{
    return (before->persistent_clusters != 0 && r->cut_edges == 0) ?
           before->persistent_hits : r->post.persistent_hits;
}

static double repair_effective_excess(const ShortcutRepair *r,
                                      const ShortcutAudit *before)
{
    return (before->persistent_clusters != 0 && r->cut_edges == 0) ?
           before->persistent_excess : r->post.persistent_excess;
}

/* Lexicographic geometry-processing objective: honor hard topology/size
 * bounds first; then remove semantic shortcut energy; then minimize the new
 * fragmentation and finally the topological/metric size of the cut. */
static int repair_trial_better(const ShortcutRepair *a,
                               const ShortcutRepair *b,
                               const ShortcutAudit *before,
                               const Params *P, size_t nfaces)
{
    int av = repair_trial_selectable(a, before, P, nfaces);
    int bv = repair_trial_selectable(b, before, P, nfaces);
    size_t at = repair_effective_tracks(a, before);
    size_t bt = repair_effective_tracks(b, before);
    size_t ah = repair_effective_hits(a, before);
    size_t bh = repair_effective_hits(b, before);
    double ae = repair_effective_excess(a, before);
    double be = repair_effective_excess(b, before);
    size_t af = positive_size_delta(a->after.components,
                                    a->before.components) +
                positive_size_delta(a->after.small_components,
                                    a->before.small_components);
    size_t bf = positive_size_delta(b->after.components,
                                    b->before.components) +
                positive_size_delta(b->after.small_components,
                                    b->before.small_components);
    if (av != bv) return av > bv;
    if ((at == 0) != (bt == 0)) return at == 0;
    if (fabs(ae - be) > 1e-9) return ae < be;
    if (at != bt) return at < bt;
    if (ah != bh) return ah < bh;
    if (af != bf) return af < bf;
    if (a->cut_rank != b->cut_rank) return a->cut_rank < b->cut_rank;
    if (a->cut_edges != b->cut_edges) return a->cut_edges < b->cut_edges;
    return fabs(a->phase_offset) < fabs(b->phase_offset);
}

static char *repair_trial_path(const char *base, const char *tag)
{
    size_t n;
    char *path;
    if (base == NULL) return NULL;
    n = strlen(base) + strlen(tag) + 1;
    path = (char *)malloc(n);
    if (path != NULL) snprintf(path, n, "%s%s", base, tag);
    return path;
}

/* is coord within seam_tol of a multiple of seam_pitch on ANY of the 3 axes? */
static int near_seam(double z, double y, double x, double pitch, double tol)
{
    if (pitch <= 0) return 0;
    double az = fabs(z - pitch * floor(z / pitch + 0.5));
    double ay = fabs(y - pitch * floor(y / pitch + 0.5));
    double ax = fabs(x - pitch * floor(x / pitch + 0.5));
    return (az <= tol || ay <= tol || ax <= tol);
}

/* A short two-plane event is promoted to a persistent bridge only when its
 * open section contour is clipped by the actual node domain. Hierarchy-internal
 * child seams are not domain boundaries and therefore must not supply this
 * evidence. Standalone callers without a box retain periodic seam behavior. */
static int near_boundary_box(double z, double y, double x, const Params *P)
{
    const double q[3] = { z, y, x };
    for (int axis = 0; axis < 3; axis++) {
        if (fabs(q[axis] - P->boundary_box[2*axis]) <= P->seam_tol ||
            fabs(q[axis] - P->boundary_box[2*axis+1]) <= P->seam_tol)
            return 1;
    }
    return 0;
}

/* Unwrap theta over mesh adjacency before measuring phase span.  The previous
 * component-mean reference was intentionally leaf-cube-local and produced a
 * false one-turn jump when a valid large region crossed atan2's branch cut (or
 * covered more than pi radians).  Adjacency unwrapping removes only angular
 * 2*pi jumps; a radial inter-wrap fusion keeps its ~integer r/pitch offset and
 * therefore remains visible to the broad-fusion gate. */
static int broad_phase_spans(const FVec *V, const IVec *F, const Params *P,
                             const int32_t *vrank, Comp *comp, size_t ncomp,
                             double core_r)
{
    size_t nv = V->n;
    size_t *degree = (size_t *)calloc(nv, sizeof(*degree));
    size_t *row = (size_t *)malloc((nv + 1) * sizeof(*row));
    size_t *cursor = (size_t *)malloc(nv * sizeof(*cursor));
    double *raw = (double *)malloc(nv * sizeof(*raw));
    double *theta = (double *)malloc(nv * sizeof(*theta));
    double *radius = (double *)malloc(nv * sizeof(*radius));
    size_t *queue = (size_t *)malloc(nv * sizeof(*queue));
    uint8_t *seen = (uint8_t *)calloc(nv, 1);
    int32_t *adj = NULL;
    if (!degree || !row || !cursor || !raw || !theta || !radius ||
        !queue || !seen) goto oom;
    for (size_t f = 0; f < F->n; f++) {
        int32_t tri[3] = {F->f[f*3+0], F->f[f*3+1], F->f[f*3+2]};
        if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0 ||
            (size_t)tri[0] >= nv || (size_t)tri[1] >= nv ||
            (size_t)tri[2] >= nv) continue;
        for (int e = 0; e < 3; e++) {
            degree[(size_t)tri[e]]++;
            degree[(size_t)tri[(e+1)%3]]++;
        }
    }
    row[0] = 0;
    for (size_t i = 0; i < nv; i++) row[i+1] = row[i] + degree[i];
    adj = (int32_t *)malloc((row[nv] ? row[nv] : 1) * sizeof(*adj));
    if (!adj) goto oom;
    memcpy(cursor, row, nv * sizeof(*cursor));
    for (size_t f = 0; f < F->n; f++) {
        int32_t tri[3] = {F->f[f*3+0], F->f[f*3+1], F->f[f*3+2]};
        if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0 ||
            (size_t)tri[0] >= nv || (size_t)tri[1] >= nv ||
            (size_t)tri[2] >= nv) continue;
        for (int e = 0; e < 3; e++) {
            size_t a = (size_t)tri[e], b = (size_t)tri[(e+1)%3];
            adj[cursor[a]++] = (int32_t)b;
            adj[cursor[b]++] = (int32_t)a;
        }
    }
    for (size_t i = 0; i < nv; i++) {
        double uy, ux;
        axis_table_eval(P->axis, (double)V->v[i*3+0],
                        P->umb_y, P->umb_x, &uy, &ux);
        double dy = (double)V->v[i*3+1] - uy;
        double dx = (double)V->v[i*3+2] - ux;
        raw[i] = atan2(dy, dx);
        radius[i] = hypot(dy, dx);
    }
    for (size_t seed = 0; seed < nv; seed++) {
        if (seen[seed] || vrank[seed] < 0) continue;
        size_t head = 0, tail = 0;
        seen[seed] = 1; theta[seed] = raw[seed]; queue[tail++] = seed;
        while (head < tail) {
            size_t u = queue[head++];
            for (size_t at = row[u]; at < row[u+1]; at++) {
                size_t v = (size_t)adj[at];
                if (vrank[v] != vrank[u] || seen[v]) continue;
                double turns = floor((theta[u] - raw[v]) / TWO_PI + 0.5);
                theta[v] = raw[v] + turns * TWO_PI;
                seen[v] = 1; queue[tail++] = v;
            }
        }
    }
    for (size_t r = 0; r < ncomp; r++) {
        comp[r].phase_min = 1e30; comp[r].phase_max = -1e30;
        comp[r].outer_vertices = 0;
    }
    for (size_t i = 0; i < nv; i++) {
        int32_t rk = vrank[i];
        if (rk < 0 || radius[i] < core_r) continue;
        double phase = rturns(P, radius[i]) - theta[i] / TWO_PI;
        if (phase < comp[rk].phase_min) comp[rk].phase_min = phase;
        if (phase > comp[rk].phase_max) comp[rk].phase_max = phase;
        comp[rk].outer_vertices++;
    }
    free(adj); free(seen); free(queue); free(radius); free(theta); free(raw);
    free(cursor); free(row); free(degree);
    return 0;
oom:
    free(adj); free(seen); free(queue); free(radius); free(theta); free(raw);
    free(cursor); free(row); free(degree);
    return -1;
}

/* a detected split candidate pair, aggregated per (rankA<rankB) */
typedef struct { int32_t ra, rb; size_t n; double min_gap, sum_gap; size_t rep_u, rep_v; double rep_dw; } SplitAgg;

static int run(const char *in, const Params *P)
{
    FVec V = {0}; IVec F = {0};
    if (read_obj(in, &V, &F) != 0) { free(V.v); free(F.f); return 1; }
    size_t nv = V.n, nf = F.n;
    if (nv == 0 || nf == 0) { fprintf(stderr, "wind_audit: empty mesh\n"); free(V.v); free(F.f); return 1; }
    Params effective;
    if (P->boundary_box_from_mesh) {
        effective = *P;
        for (int axis = 0; axis < 3; axis++) {
            double lo = DBL_MAX, hi = -DBL_MAX;
            for (size_t i = 0; i < nv; i++) {
                double value = V.v[i*3+(size_t)axis];
                if (!isfinite(value)) {
                    fprintf(stderr, "wind_audit: non-finite OBJ vertex\n");
                    free(V.v); free(F.f); return 1;
                }
                if (value < lo) lo = value;
                if (value > hi) hi = value;
            }
            if (hi - lo < 1e-6) { lo -= 0.5; hi += 0.5; }
            effective.boundary_box[axis*2] = lo;
            effective.boundary_box[axis*2+1] = hi;
        }
        effective.boundary_box_set = 1;
        P = &effective;
    }
    double core_r = P->core_pitches * P->pitch;

    /* --- connected components --- */
    int32_t *par = (int32_t *)malloc(nv * sizeof(int32_t));
    int32_t *usz = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!par || !usz) { fprintf(stderr, "wind_audit: OOM\n"); free(par); free(usz); free(V.v); free(F.f); return 1; }
    for (size_t i = 0; i < nv; i++) { par[i] = (int32_t)i; usz[i] = 1; }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F.f[f*3+0], b = F.f[f*3+1], c = F.f[f*3+2];
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv) continue;
        uf_union(par, usz, a, b); uf_union(par, usz, b, c);
    }
    size_t *root_faces = (size_t *)calloc(nv, sizeof(size_t));
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F.f[f*3+0];
        if (a < 0 || (size_t)a >= nv) continue;
        root_faces[uf_find(par, a)]++;
    }
    size_t ncomp = 0;
    for (size_t i = 0; i < nv; i++) if (root_faces[i] > 0) ncomp++;
    Comp *comp = (Comp *)malloc((ncomp ? ncomp : 1) * sizeof(Comp));
    size_t ci = 0;
    for (size_t i = 0; i < nv; i++) if (root_faces[i] > 0) {
        comp[ci].root = (int32_t)i; comp[ci].faces = root_faces[i]; comp[ci].rank = 0;
        comp[ci].vertices = (size_t)usz[i];
        comp[ci].edges = comp[ci].boundary_edges = comp[ci].boundary_loops = 0;
        comp[ci].nonmanifold_edges = comp[ci].boundary_irregular_vertices = 0;
        comp[ci].euler_chi = 0; comp[ci].genus = 0.0;
        comp[ci].rmin = 1e30; comp[ci].rmax = -1e30; comp[ci].zmin = 1e30f; comp[ci].zmax = -1e30f;
        comp[ci].merger_edges = 0;
        comp[ci].phase_min = 1e30; comp[ci].phase_max = -1e30;
        comp[ci].outer_vertices = 0;
        ci++;
    }
    qsort(comp, ncomp, sizeof(Comp), comp_cmp);
    int32_t *root2rank = (int32_t *)malloc(nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) root2rank[i] = -1;
    for (size_t r = 0; r < ncomp; r++) { comp[r].rank = (int32_t)r; root2rank[comp[r].root] = (int32_t)r; }

    /* per-vertex rank + radius; accumulate comp r/z bounds */
    int32_t *vrank = (int32_t *)malloc(nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) {
        int32_t rk = root2rank[uf_find(par, (int32_t)i)];
        vrank[i] = rk;
        if (rk < 0) continue;
        double r = vradius(V.v, i, P->umb_y, P->umb_x, P->axis);
        float z = V.v[i*3+0];
        if (r < comp[rk].rmin) comp[rk].rmin = r;
        if (r > comp[rk].rmax) comp[rk].rmax = r;
        if (z < comp[rk].zmin) comp[rk].zmin = z;
        if (z > comp[rk].zmax) comp[rk].zmax = z;
    }

    if (broad_phase_spans(&V, &F, P, vrank, comp, ncomp, core_r) != 0) {
        fprintf(stderr, "wind_audit: OOM unwrapping large-region phase\n");
        free(vrank); free(par); free(usz); free(root_faces); free(comp);
        free(root2rank); free(V.v); free(F.f);
        return 1;
    }
    size_t broad_fusion_components = 0, broad_fusion_faces = 0;
    double broad_fusion_max_span = 0.0;
    if (P->local_span_tol > 0.0) {
        for (size_t r = 0; r < ncomp; r++) {
            double span = comp[r].phase_max - comp[r].phase_min;
            /* Components touching the core may legitimately turn tightly, so
             * the broad phase-span gate only judges wholly outer components. */
            if (comp[r].faces >= P->min_faces && comp[r].rmin >= core_r
                && comp[r].outer_vertices > 0 && span > P->local_span_tol) {
                broad_fusion_components++;
                broad_fusion_faces += comp[r].faces;
                if (span > broad_fusion_max_span) broad_fusion_max_span = span;
            }
        }
    }

    /* --- MERGER edges: dedup unique edges, test |dw| --- */
    /* dedup edge set (open addressing over packed u<<32|v) */
    size_t eslots = 1; while (eslots < nf * 4) eslots <<= 1;
    int64_t *eset = (int64_t *)malloc(eslots * sizeof(int64_t));
    uint8_t *ecount = (uint8_t *)calloc(eslots, 1);
    if (!eset || !ecount) {
        fprintf(stderr, "wind_audit: OOM allocating edge table\n");
        free(eset); free(ecount);
        free(vrank); free(par); free(usz); free(root_faces); free(comp); free(root2rank);
        free(V.v); free(F.f);
        return 1;
    }
    for (size_t i = 0; i < eslots; i++) eset[i] = -1;
    size_t emask = eslots - 1;
    size_t merge_outer = 0, merge_core = 0, edges_total = 0;
    /* merger-edge geometry: are they ~one-pitch radial jumps? */
    double mlen_sum = 0, mdr_sum = 0;           /* outer merger edges */
    double olen_sum = 0; size_t outer_edges = 0; /* all outer edges (contrast) */
    size_t dwhist[5] = {0,0,0,0,0};             /* |dw| bins over outer edges */
    size_t merge_ft = 0, merge_ft_seam = 0;     /* full-turn (|dw| in [0.7,1.3]) fusions, seam-classified */
    int artifact_failed = 0;
    uint8_t *face_flag = NULL;
    if (P->dump_mergers) {
        face_flag = (uint8_t *)calloc(nf, 1);
        if (face_flag == NULL) {
            fprintf(stderr,
                    "wind_audit: cannot allocate mergers overlay state\n");
            artifact_failed = 1;
        }
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t vv[3] = { F.f[f*3+0], F.f[f*3+1], F.f[f*3+2] };
        if (vv[0]<0||vv[1]<0||vv[2]<0||(size_t)vv[0]>=nv||(size_t)vv[1]>=nv||(size_t)vv[2]>=nv) continue;
        for (int e = 0; e < 3; e++) {
            int32_t u = vv[e], w = vv[(e+1)%3];
            int32_t a = u < w ? u : w, b = u < w ? w : u;
            int64_t key = ((int64_t)a << 32) | (uint32_t)b;
            size_t i = sh_hash(key, emask);
            int seen = 0;
            for (;;) {
                if (eset[i] == -1) { eset[i] = key; ecount[i] = 1; break; }
                if (eset[i] == key) {
                    if (ecount[i] < 255) ecount[i]++;
                    seen = 1; break;
                }
                i = (i+1)&emask;
            }
            if (seen) continue;
            edges_total++;
            double dw = dwind(V.v, (size_t)a, (size_t)b, P->umb_y,
                              P->umb_x, P->pitch, P->axis, P->ptab);
            double ra = vradius(V.v, (size_t)a, P->umb_y, P->umb_x,
                                P->axis);
            double rb = vradius(V.v, (size_t)b, P->umb_y, P->umb_x,
                                P->axis);
            double rmin = ra < rb ? ra : rb;
            double adw = fabs(dw);
            if (rmin >= core_r) {              /* outer edge: tally geometry + hist */
                double dz = (double)V.v[a*3+0]-V.v[b*3+0], dy = (double)V.v[a*3+1]-V.v[b*3+1], dx = (double)V.v[a*3+2]-V.v[b*3+2];
                double len = sqrt(dz*dz + dy*dy + dx*dx);
                outer_edges++; olen_sum += len;
                if      (adw < 0.40) dwhist[0]++;
                else if (adw < 0.70) dwhist[1]++;
                else if (adw < 1.30) dwhist[2]++;   /* ~one-turn radial jump */
                else if (adw < 2.30) dwhist[3]++;   /* ~two-turn */
                else                 dwhist[4]++;
                if (adw > P->merge_tol) { mlen_sum += len; mdr_sum += fabs(rb - ra); }
                if (adw >= 0.70 && adw <= 1.30) {   /* unambiguous full-turn fusion */
                    merge_ft++;
                    double mz = 0.5*((double)V.v[a*3+0]+V.v[b*3+0]);
                    double my = 0.5*((double)V.v[a*3+1]+V.v[b*3+1]);
                    double mx = 0.5*((double)V.v[a*3+2]+V.v[b*3+2]);
                    if (merge_ft <= 8)
                        printf("    FULL-TURN edge #%zu at (z %.0f, y %.0f, "
                               "x %.0f)  r %.1f->%.1f  |dw|=%.2f len=%.1f\n",
                               merge_ft, mz, my, mx, ra, rb, adw, len);
                    if (near_seam(mz, my, mx, P->seam_pitch, P->seam_tol)) merge_ft_seam++;
                    int32_t rk = vrank[a];       /* per-comp column = reliable full-turn count */
                    if (rk >= 0 && rk == vrank[b]) comp[rk].merger_edges++;
                }
            }
            if (adw > P->merge_tol) {
                if (rmin < core_r) merge_core++;
                else merge_outer++;
                if (face_flag) face_flag[f] = 1;
            }
        }
    }
    /* Per-component Euler topology.  For a connected orientable manifold
     * with boundary, chi = V-E+F = 2-2g-b.  The boundary-edge graph gives b;
     * grid_weld's manifold gate separately guarantees that its connected
     * components are cycles rather than branched boundary graphs. */
    int32_t *bpar = (int32_t *)malloc(nv * sizeof(int32_t));
    int32_t *bsz = (int32_t *)malloc(nv * sizeof(int32_t));
    uint16_t *bdegree = (uint16_t *)calloc(nv, sizeof(*bdegree));
    if (!bpar || !bsz || !bdegree) {
        fprintf(stderr, "wind_audit: OOM allocating topology audit\n");
        free(bpar); free(bsz); free(bdegree); free(ecount); free(eset);
        free(face_flag); free(vrank); free(par); free(usz); free(root_faces);
        free(comp); free(root2rank); free(V.v); free(F.f);
        return 1;
    }
    for (size_t i = 0; i < nv; i++) { bpar[i] = (int32_t)i; bsz[i] = 1; }
    size_t boundary_edges = 0;
    for (size_t i = 0; i < eslots; i++) if (eset[i] != -1) {
        uint64_t key = (uint64_t)eset[i];
        int32_t a = (int32_t)(uint32_t)(key >> 32);
        int32_t b = (int32_t)(uint32_t)key;
        int32_t rk = vrank[a];
        if (rk < 0 || rk != vrank[b]) continue;
        comp[rk].edges++;
        if (ecount[i] > 2) comp[rk].nonmanifold_edges++;
        if (ecount[i] == 1) {
            boundary_edges++;
            comp[rk].boundary_edges++;
            if (bdegree[a] < UINT16_MAX) bdegree[a]++;
            if (bdegree[b] < UINT16_MAX) bdegree[b]++;
            uf_union(bpar, bsz, a, b);
        }
    }
    for (size_t i = 0; i < nv; i++) if (bdegree[i] > 0) {
        int32_t rk = vrank[i];
        if (rk < 0) continue;
        if (bdegree[i] != 2) comp[rk].boundary_irregular_vertices++;
        if (uf_find(bpar, (int32_t)i) == (int32_t)i)
            comp[rk].boundary_loops++;
    }
    for (size_t r = 0; r < ncomp; r++) {
        comp[r].euler_chi = (int64_t)comp[r].vertices -
                            (int64_t)comp[r].edges +
                            (int64_t)comp[r].faces;
        comp[r].genus = 0.5 * (2.0 - (double)comp[r].boundary_loops -
                               (double)comp[r].euler_chi);
    }
    free(bdegree); free(bsz); free(bpar);
    free(ecount);
    free(eset);

    /* --- SPLIT CONTACTS: spatial hash, cross-component close low-dw pairs ---
     * The hierarchy's iterative bridge-cut lane does not consume this separate
     * diagnostic (and repeats it unchanged at every pass), so it may opt out;
     * the independent final audit retains it.  Counts are vertex-pair density,
     * not a number of missing seams, and are not a repair authorization. */
    SHash H; memset(&H, 0, sizeof(H));
    size_t hslots = 0, aslots = 0;
    SplitAgg *agg = NULL;
    size_t nagg = 0;
    double gap2 = P->split_gap * P->split_gap;
    size_t split_pairs_total = 0, split_pairs_seam = 0;
    if (P->split_audit) {
        hslots = 1; while (hslots < nv * 2) hslots <<= 1;
        H.cell = (HCell *)malloc(hslots * sizeof(HCell));
        H.next = (int32_t *)malloc(nv * sizeof(int32_t));
        aslots = 4096;
        agg = (SplitAgg *)calloc(aslots, sizeof(SplitAgg));
        if (H.cell == NULL || H.next == NULL || agg == NULL) {
            fprintf(stderr, "wind_audit: OOM allocating split audit\n");
            free(agg); free(H.cell); free(H.next); free(face_flag);
            free(vrank); free(par); free(usz); free(root_faces); free(comp);
            free(root2rank); free(V.v); free(F.f);
            return 1;
        }
        H.mask = hslots - 1; H.inv = 1.0 / P->split_gap; H.off = 1 << 20;
        for (size_t i = 0; i < hslots; i++) {
            H.cell[i].head = -2; H.cell[i].key = 0;
        }
        for (size_t i = 0; i < nv; i++) {
            if (vrank[i] < 0) { H.next[i] = -1; continue; }
            int64_t key = sh_key(&H, V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
            int32_t *slot = sh_slot(&H, key, 1);
            H.next[i] = *slot; *slot = (int32_t)i;
        }
        for (size_t i = 0; i < aslots; i++) {
            agg[i].ra = -1; agg[i].rb = -1;
        }
        for (size_t i = 0; i < nv; i++) {
            int32_t rki = vrank[i]; if (rki < 0) continue;
            double az = V.v[i*3+0], ay = V.v[i*3+1], ax = V.v[i*3+2];
            for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
                int64_t key = sh_key(&H, az + dz*P->split_gap, ay + dy*P->split_gap, ax + dx*P->split_gap);
                int32_t *slot = sh_slot(&H, key, 0);
                if (!slot) continue;
                int scan = 0;
                for (int32_t j = *slot; j >= 0 && scan < 96; j = H.next[j], scan++) {
                    if ((size_t)j <= i) continue;          /* unordered, once */
                    int32_t rkj = vrank[j]; if (rkj < 0 || rkj == rki) continue;
                    double d2 = (az-V.v[j*3+0])*(az-V.v[j*3+0]) + (ay-V.v[j*3+1])*(ay-V.v[j*3+1]) + (ax-V.v[j*3+2])*(ax-V.v[j*3+2]);
                    if (d2 > gap2) continue;
                    double dw = dwind(V.v, i, (size_t)j, P->umb_y, P->umb_x,
                                      P->pitch, P->axis, P->ptab);
                    if (fabs(dw) > P->split_wtol) continue;
                    split_pairs_total++;
                    if (near_seam(0.5*(az+V.v[j*3+0]), 0.5*(ay+V.v[j*3+1]), 0.5*(ax+V.v[j*3+2]), P->seam_pitch, P->seam_tol)) split_pairs_seam++;
                    int32_t ra = rki < rkj ? rki : rkj, rb = rki < rkj ? rkj : rki;
                    uint64_t hk = ((uint64_t)ra * 2654435761u) ^ ((uint64_t)rb * 40503u);
                    size_t s = (size_t)hk & (aslots - 1), guard = 0;
                    for (;;) {
                        if (agg[s].ra == -1) { agg[s].ra = ra; agg[s].rb = rb; agg[s].n = 0; agg[s].min_gap = 1e30; agg[s].sum_gap = 0; nagg++; }
                        if (agg[s].ra == ra && agg[s].rb == rb) break;
                        s = (s + 1) & (aslots - 1);
                        if (++guard >= aslots) { s = SIZE_MAX; break; }
                    }
                    if (s == SIZE_MAX) continue;
                    double gap = sqrt(d2);
                    agg[s].n++; agg[s].sum_gap += gap;
                    if (gap < agg[s].min_gap) { agg[s].min_gap = gap; agg[s].rep_u = i; agg[s].rep_v = (size_t)j; agg[s].rep_dw = dw; }
                }
            }
        }
    }

    ShortcutAudit shortcut;
    int shortcut_rc = shortcut_audit(&V, &F, P, &shortcut);
    if (shortcut_rc != 0)
        fprintf(stderr, "wind_audit: local-shortcut audit failed (OOM)\n");
    ShortcutRepair repair;
    memset(&repair, 0, sizeof(repair));
    int repair_rc = 0;
    int repair_requested = P->repair_path != NULL || P->repair_dry_run;
    if (repair_requested && shortcut_rc == 0) {
        if (!P->repair_auto_offset) {
            repair_rc = shortcut_seam_repair(&V, &F, P, &shortcut, &repair);
            repair.auto_trials = 1;
        } else {
            ShortcutRepair trial[17];
            int trial_ok[17] = {0};
            const double offset_delta[17] = {
                0.0, -0.05, 0.05, -0.10, 0.10, -0.15, 0.15,
                -0.20, 0.20, -0.25, 0.25, -0.30, 0.30,
                -0.35, 0.35, -0.40, 0.40
            };
            const char *trial_tag[17] = {
                "", ".auto_m05", ".auto_p05", ".auto_m10",
                ".auto_p10", ".auto_m15", ".auto_p15",
                ".auto_m20", ".auto_p20", ".auto_m25", ".auto_p25",
                ".auto_m30", ".auto_p30", ".auto_m35", ".auto_p35",
                ".auto_m40", ".auto_p40"
            };
            double offsets[17];
            char *mesh_path[17] = {NULL};
            char *csv_path[17] = {NULL};
            int ntrial = 1, chosen = 0, base_decisive = 0;
            int base_needs_search = 0, auto_suppressed = 0;
            Params TP = *P;
            memset(trial, 0, sizeof(trial));
            for (int t = 0; t < 17; t++)
                offsets[t] = P->repair_phase_offset + offset_delta[t];
            mesh_path[0] = (char *)P->repair_path;
            csv_path[0] = (char *)P->repair_chain_csv;
            repair_rc = shortcut_seam_repair(&V, &F, &TP, &shortcut,
                                              &trial[0]);
            trial_ok[0] = repair_rc == 0;
            base_decisive = trial_ok[0] &&
                repair_trial_decisive(&trial[0], &shortcut, P, F.n);
            base_needs_search = trial_ok[0] &&
                repair_needs_offset_search(&trial[0], &shortcut, P, F.n);
            auto_suppressed = base_needs_search &&
                !repair_auto_search_allowed(P, &shortcut);
            if (base_needs_search && !auto_suppressed) {
                int paths_ok = 1;
                for (int t = 1; t < 7; t++) {
                    mesh_path[t] = repair_trial_path(P->repair_path,
                                                     trial_tag[t]);
                    csv_path[t] = repair_trial_path(P->repair_chain_csv,
                                                    trial_tag[t]);
                    if ((P->repair_path != NULL && mesh_path[t] == NULL) ||
                        (P->repair_chain_csv != NULL && csv_path[t] == NULL))
                        paths_ok = 0;
                }
                if (!paths_ok) {
                    repair_rc = -1;
                    snprintf(trial[0].reason, sizeof(trial[0].reason),
                             "out of memory preparing phase-offset trials");
                } else {
                    for (int t = 1; t < 7; t++) {
                        if (fabs(offsets[t]) >= 0.5) continue;
                        TP = *P;
                        TP.repair_phase_offset = offsets[t];
                        TP.repair_path = mesh_path[t];
                        TP.repair_chain_csv = csv_path[t];
                        if (mesh_path[t] != NULL) remove(mesh_path[t]);
                        if (csv_path[t] != NULL) remove(csv_path[t]);
                        trial_ok[t] =
                            shortcut_seam_repair(&V, &F, &TP, &shortcut,
                                                 &trial[t]) == 0;
                        ntrial++;
                    }
                    for (int t = 1; t < 7; t++) {
                        if (trial_ok[t] &&
                            repair_trial_better(&trial[t], &trial[chosen],
                                                &shortcut, P, F.n))
                            chosen = t;
                    }
                    /* Backtrack to more distant regular values only when the
                     * local stencil cannot produce an admissible monotone cut.
                     * The normal path therefore remains seven trials. */
                    if (!repair_trial_selectable(&trial[chosen], &shortcut,
                                                 P, F.n)) {
                        paths_ok = 1;
                        for (int t = 7; t < 17; t++) {
                            mesh_path[t] = repair_trial_path(P->repair_path,
                                                             trial_tag[t]);
                            csv_path[t] = repair_trial_path(
                                P->repair_chain_csv, trial_tag[t]);
                            if ((P->repair_path != NULL &&
                                 mesh_path[t] == NULL) ||
                                (P->repair_chain_csv != NULL &&
                                 csv_path[t] == NULL)) paths_ok = 0;
                        }
                        if (!paths_ok) {
                            repair_rc = -1;
                            trial_ok[chosen] = 0;
                            snprintf(trial[chosen].reason,
                                     sizeof(trial[chosen].reason),
                                     "out of memory preparing wide offset trials");
                        } else {
                            for (int t = 7; t < 17; t++) {
                                if (fabs(offsets[t]) >= 0.5) continue;
                                TP = *P;
                                TP.repair_phase_offset = offsets[t];
                                TP.repair_path = mesh_path[t];
                                TP.repair_chain_csv = csv_path[t];
                                if (mesh_path[t] != NULL) remove(mesh_path[t]);
                                if (csv_path[t] != NULL) remove(csv_path[t]);
                                trial_ok[t] = shortcut_seam_repair(
                                    &V, &F, &TP, &shortcut, &trial[t]) == 0;
                                ntrial++;
                            }
                            for (int t = 7; t < 17; t++) {
                                if (trial_ok[t] && repair_trial_better(
                                        &trial[t], &trial[chosen], &shortcut,
                                        P, F.n)) chosen = t;
                            }
                        }
                    }
                    repair_rc = trial_ok[chosen] ? 0 : -1;
                }
            }
            repair = trial[chosen];
            repair.auto_trials = ntrial;
            repair.auto_suppressed = auto_suppressed;
            repair.auto_decisive = chosen == 0 && ntrial == 1 &&
                base_decisive;
            memset(&trial[chosen].post, 0, sizeof(trial[chosen].post));

            if (chosen != 0 && P->repair_path != NULL) {
                remove(P->repair_path);
                if (repair.output_written) {
                    if (rename(mesh_path[chosen], P->repair_path) != 0) {
                        repair_rc = -1;
                        repair.output_written = 0;
                        snprintf(repair.reason, sizeof(repair.reason),
                                 "cannot promote selected phase-offset mesh");
                    }
                }
            }
            if (chosen != 0 && P->repair_chain_csv != NULL) {
                FILE *probe = fopen(csv_path[chosen], "rb");
                remove(P->repair_chain_csv);
                if (probe != NULL) {
                    fclose(probe);
                    if (rename(csv_path[chosen], P->repair_chain_csv) != 0) {
                        repair_rc = -1;
                        snprintf(repair.reason, sizeof(repair.reason),
                                 "cannot promote selected phase-offset CSV");
                    }
                }
            }
            for (int t = 1; t < 17; t++) {
                if (mesh_path[t] != NULL) remove(mesh_path[t]);
                if (csv_path[t] != NULL) remove(csv_path[t]);
                free(mesh_path[t]);
                free(csv_path[t]);
            }
            for (int t = 0; t < 17; t++)
                shortcut_audit_free(&trial[t].post);
        }
        if (repair_rc != 0)
            fprintf(stderr, "wind_audit: shortcut seam repair failed: %s\n",
                    repair.reason);
    }

    /* ---- report ---- */
    printf("wind_audit: %s\n", in);
    if (P->axis != NULL) {
        printf("  umbilicus=curve[%zu samples, z %.1f..%.1f]  pitch=%.3f  merge_tol=%.2f  split_wtol=%.2f  split_gap=%.1f  core<%.1f vox\n",
               P->axis->n, P->axis->z[0], P->axis->z[P->axis->n-1],
               P->pitch, P->merge_tol, P->split_wtol, P->split_gap, core_r);
        printf("  axis_table=%s (linear interpolation + endpoint tangent extrapolation)\n",
               P->axis_table_path != NULL ? P->axis_table_path : "<memory>");
    } else {
        printf("  umbilicus=(y %.1f, x %.1f)  pitch=%.3f  merge_tol=%.2f  split_wtol=%.2f  split_gap=%.1f  core<%.1f vox\n",
               P->umb_y, P->umb_x, P->pitch, P->merge_tol,
               P->split_wtol, P->split_gap, core_r);
    }
    printf("  verts=%zu  faces=%zu  components=%zu  unique_edges=%zu\n", nv, nf, ncomp, edges_total);
    printf("\n  == MERGERS (radial short-circuits, |dw|>%.2f across one edge) ==\n", P->merge_tol);
    printf("    outer (r>=%.0f, reliable): %zu edges (%.2f%% of %zu outer)    core (r<%.0f): %zu edges\n",
           core_r, merge_outer, outer_edges ? 100.0*(double)merge_outer/(double)outer_edges : 0.0,
           outer_edges, core_r, merge_core);
    printf("    outer |dw| histogram:  <0.40 (ok)=%zu  0.40-0.70=%zu  0.70-1.30 (~1 turn)=%zu  1.30-2.30 (~2)=%zu  >2.30=%zu\n",
           dwhist[0], dwhist[1], dwhist[2], dwhist[3], dwhist[4]);
    printf("    merger-edge geometry:  mean_len=%.2f vox  mean|dr|=%.2f vox   (all outer edges mean_len=%.2f vox; pitch=%.2f)\n",
           merge_outer ? mlen_sum/(double)merge_outer : 0.0, merge_outer ? mdr_sum/(double)merge_outer : 0.0,
           outer_edges ? olen_sum/(double)outer_edges : 0.0, P->pitch);
    printf("    FULL-TURN fusions (|dw| in [0.7,1.3], the reliable count): %zu   of which near a %.0f-vox seam: %zu (%.1f%%)\n",
           merge_ft, P->seam_pitch, merge_ft_seam, merge_ft ? 100.0*(double)merge_ft_seam/(double)merge_ft : 0.0);
    printf("    mesh boundary edges: %zu\n", boundary_edges);
    if (P->local_span_tol > 0.0) {
        printf("\n  == BROAD FUSIONS (adjacency-unwrapped phase span > %.2f turns) ==\n",
               P->local_span_tol);
        printf("    wholly-outer components: %zu  faces involved: %zu  max phase span: %.3f turns\n",
               broad_fusion_components, broad_fusion_faces, broad_fusion_max_span);
    }
    if (P->shortcut_tol > 0.0) {
        printf("\n  == LOCAL SHORTCUTS (axial contour delta-w >= %.2f within %.1f pitches) ==\n",
               P->shortcut_tol, P->shortcut_length_pitches);
        if (P->boundary_box_set)
            printf("    clipped-track domain: z[%.1f,%.1f] y[%.1f,%.1f] x[%.1f,%.1f] (tol %.1f)\n",
                   P->boundary_box[0], P->boundary_box[1],
                   P->boundary_box[2], P->boundary_box[3],
                   P->boundary_box[4], P->boundary_box[5], P->seam_tol);
        printf("    tracker: %s\n", P->shortcut_track_mutual ?
               "mutual-nearest temporal chains" : "proximity-union clusters");
        printf("    sampled planes: %zu  face-plane incidences: %zu  outer contours: %zu (simple %zu, fallback %zu)  section hits: %zu  tracks: %zu\n",
               shortcut.planes_tested, shortcut.section_face_incidents,
               shortcut.curves_tested,
               shortcut.simple_curves, shortcut.fallback_curves,
               shortcut.nhit, shortcut.ncluster);
        printf("    persistent tracks (>= %d planes, or >=2 boundary-clipped): %zu  hits: %zu  max delta-w: %.3f  integrated excess: %.3f turn-sections\n",
               P->shortcut_min_planes, shortcut.persistent_clusters,
               shortcut.persistent_hits, shortcut.persistent_max_delta_w,
               shortcut.persistent_excess);
        printf("    %-5s %7s %7s %8s %8s %8s %5s   %s\n",
               "track", "planes", "hits", "max_dw", "min_len", "z_span", "clip",
               "centroid (z,y,x)  bbox (z,y,x)..(z,y,x)");
        int shown = 0;
        for (size_t r = 0; r < shortcut.ncluster && shown < P->top; r++) {
            const ShortcutCluster *g = &shortcut.cluster[r];
            if (!g->persistent) continue;
            printf("    #%-4zu %7zu %7zu %8.3f %8.2f %8.1f %5s   (%.1f,%.1f,%.1f)  (%.1f,%.1f,%.1f)..(%.1f,%.1f,%.1f)\n",
                   r, g->planes, g->hits, g->max_delta_w, g->min_distance,
                   g->z_max - g->z_min,
                   g->boundary_supported ? "yes" : "no",
                   g->centroid[0], g->centroid[1], g->centroid[2],
                   g->bbox_min[0], g->bbox_min[1], g->bbox_min[2],
                   g->bbox_max[0], g->bbox_max[1], g->bbox_max[2]);
            shown++;
        }
    }
    if (repair_requested) {
        printf("\n  == SUPPORTED SHORTCUT SEAM REPAIR ==\n");
        printf("    chain policy: %s\n",
               P->repair_cover_all_supported ?
               "all hit-supported isoline components" :
               "single best-supported component per track");
        printf("    phase offset from fitted bridge midpoint: %+.3f turns",
               repair.phase_offset);
        if (P->repair_auto_offset) {
            if (repair.auto_suppressed)
                printf(" (phase sweep suppressed: %zu tracks > %zu limit)",
                       shortcut.persistent_clusters,
                       P->repair_max_auto_tracks);
            else
                printf(repair.auto_decisive ?
                       " (decisive base; phase sweep skipped)" :
                       " (selected from %d trial%s)", repair.auto_trials,
                       repair.auto_trials == 1 ? "" : "s");
        }
        printf("\n");
        if (P->repair_cover_all_supported)
            printf("    minimum hit support per selected chain: %zu\n",
                   P->repair_min_chain_support);
        printf("    labelled faces: low=%zu high=%zu  isoline candidates=%zu\n",
               repair.labelled_low_faces, repair.labelled_high_faces,
               repair.candidate_edges);
        printf("    selected chains=%zu  seam edges=%zu (%.3f%% of faces)  vertices=%zu -> %zu\n",
               repair.selected_chains, repair.cut_edges,
               nf ? 100.0*(double)repair.cut_edges/(double)nf : 0.0,
               repair.vertices_before, repair.vertices_after);
        printf("    cut graph: components=%zu  endpoints=%zu (%zu on boundary)  branch vertices=%zu (excess=%zu)  separation rank=%zu\n",
               repair.selected_cut_components,
               repair.selected_chain_endpoints,
               repair.selected_chain_boundary_endpoints,
               repair.selected_chain_branch_vertices,
               repair.selected_chain_branch_excess, repair.cut_rank);
        printf("    components=%zu -> %zu  small(<%zu faces)=%zu -> %zu  post-tracks=%zu\n",
               repair.before.components, repair.after.components, P->min_faces,
               repair.before.small_components, repair.after.small_components,
               repair.post.persistent_clusters);
        printf("    verdict: %s -- %s\n",
               repair.accepted ? "ACCEPT" :
               (repair.step_accepted ? "STEP" : "REJECT"), repair.reason);
        if (repair.output_written)
            printf("    wrote repaired mesh -> %s\n", P->repair_path);
    }

    printf("\n  == COMPONENTS (top %d by faces) ==\n", P->top);
    printf("    %-4s %10s %6s   %8s %8s %8s  %9s %9s  %8s %8s   %-23s %s\n",
           "rank", "faces", "%", "r_min", "r_max", "r_span", "turns", "phase", "z_min", "z_max", "topology V/E/b/chi/g", "ft_fusions");
    int rows = P->top < (int)ncomp ? P->top : (int)ncomp;
    for (int r = 0; r < rows; r++) {
        double span = comp[r].rmax - comp[r].rmin;
        double pspan = comp[r].outer_vertices ?
            comp[r].phase_max - comp[r].phase_min : 0.0;
        printf("    #%-3d %10zu %5.1f%%   %8.1f %8.1f %8.1f  %8.1f %8.2f  %8.0f %8.0f   %zu/%zu/%zu/%lld/%.1f %zu\n",
               r, comp[r].faces, 100.0*(double)comp[r].faces/(double)nf,
               comp[r].rmin, comp[r].rmax, span,
               PitchTable_dturns(P->ptab, comp[r].rmin, comp[r].rmax,
                                 P->pitch), pspan,
               (double)comp[r].zmin, (double)comp[r].zmax,
               comp[r].vertices, comp[r].edges, comp[r].boundary_loops,
               (long long)comp[r].euler_chi, comp[r].genus,
               comp[r].merger_edges);
    }

    /* rank split aggregates by pair count */
    printf("\n  == SPLIT CONTACTS (same-turn |dw|<%.2f, gap<%.1f vox, DIFFERENT components) ==\n", P->split_wtol, P->split_gap);
    if (!P->split_audit) printf("    disabled for iterative shortcut repair\n");
    printf("    diagnostic vertex pairs: %zu (near a %.0f-vox seam: %zu = %.1f%%)   distinct component-pairs: %zu\n",
           split_pairs_total, P->seam_pitch, split_pairs_seam,
           split_pairs_total ? 100.0*(double)split_pairs_seam/(double)split_pairs_total : 0.0, nagg);
    printf("    NOTE: proximity density is not a missing-seam count; disk cuts, pinch splits,\n"
           "          and unsupported clipped contacts remain separate by design.\n");
    /* collect + sort agg */
    SplitAgg *list = (SplitAgg *)malloc((nagg ? nagg : 1) * sizeof(SplitAgg));
    size_t nl = 0;
    for (size_t s = 0; s < aslots; s++) if (agg[s].ra != -1) list[nl++] = agg[s];
    /* simple selection of top by n */
    int srows = P->top < (int)nl ? P->top : (int)nl;
    printf("    %-9s %8s %8s %8s   %s\n", "A<->B", "pairs", "min_gap", "mean_gap", "representative (z,y,x) r  dw");
    for (int r = 0; r < srows; r++) {
        size_t best = 0; size_t bestn = 0; int have = 0;
        for (size_t s = 0; s < nl; s++) if (list[s].n > bestn) { bestn = list[s].n; best = s; have = 1; }
        if (!have) break;
        SplitAgg *g = &list[best];
        double rr = vradius(V.v, g->rep_u, P->umb_y, P->umb_x, P->axis);
        printf("    #%-3d<->#%-3d %8zu %8.2f %8.2f   (%.0f,%.0f,%.0f) r=%.0f dw=%+.3f\n",
               g->ra, g->rb, g->n, g->min_gap, g->sum_gap / (double)g->n,
               V.v[g->rep_u*3+0], V.v[g->rep_u*3+1], V.v[g->rep_u*3+2], rr, g->rep_dw);
        list[best].n = 0; /* consume */
    }

    /* ---- optional machine-readable summary (regression gate) ---- */
    if (P->json_path) {
        FILE *j = fopen(P->json_path, "wb");
        if (j == NULL) {
            fprintf(stderr, "wind_audit: cannot write json: %s\n",
                    P->json_path);
            artifact_failed = 1;
        } else {
            fprintf(j, "{\n");
            fprintf(j, "  \"verts\": %zu, \"faces\": %zu, \"components\": %zu, \"outer_edges\": %zu, \"boundary_edges\": %zu,\n",
                    nv, nf, ncomp, outer_edges, boundary_edges);
            fprintf(j, "  \"merge_ft\": %zu, \"merge_ft_seam\": %zu, \"merge_ft_seam_frac\": %.4f,\n",
                    merge_ft, merge_ft_seam, merge_ft ? (double)merge_ft_seam/(double)merge_ft : 0.0);
            fprintf(j, "  \"merge_ft_per_outer_edge\": %.6f,\n", outer_edges ? (double)merge_ft/(double)outer_edges : 0.0);
            fprintf(j, "  \"split_audit\": %s, \"split_pairs\": %zu, \"split_pairs_seam\": %zu, \"split_pairs_seam_frac\": %.4f, \"split_comp_pairs\": %zu,\n",
                    P->split_audit ? "true" : "false",
                    split_pairs_total, split_pairs_seam, split_pairs_total ? (double)split_pairs_seam/(double)split_pairs_total : 0.0, nagg);
            fprintf(j, "  \"split_semantics\": \"diagnostic_proximity_contacts_not_certified_missing_seams\",\n");
            fprintf(j, "  \"local_span_tol\": %.4f, \"broad_fusion_components\": %zu, \"broad_fusion_faces\": %zu, \"broad_fusion_max_span\": %.4f,\n",
                    P->local_span_tol, broad_fusion_components, broad_fusion_faces, broad_fusion_max_span);
            fprintf(j, "  \"shortcut_tol\": %.4f, \"shortcut_length_pitches\": %.3f, \"shortcut_z_step\": %.3f,\n",
                    P->shortcut_tol, P->shortcut_length_pitches,
                    P->shortcut_z_step);
            fprintf(j, "  \"shortcut_boundary_mode\": \"%s\",\n",
                    P->boundary_box_from_mesh ? "mesh_box" :
                    (P->boundary_box_set ? "node_box" : "periodic_seams"));
            if (P->boundary_box_set)
                fprintf(j, "  \"shortcut_boundary_box_zyx\": [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f],\n",
                        P->boundary_box[0], P->boundary_box[1],
                        P->boundary_box[2], P->boundary_box[3],
                        P->boundary_box[4], P->boundary_box[5]);
            fprintf(j, "  \"shortcut_track_policy\": \"%s\",\n",
                    P->shortcut_track_mutual ? "mutual_temporal" : "proximity_union");
            fprintf(j, "  \"shortcut_planes\": %zu, \"shortcut_section_face_incidents\": %zu, \"shortcut_curves\": %zu, \"shortcut_simple_curves\": %zu, \"shortcut_fallback_curves\": %zu, \"shortcut_hits\": %zu, \"shortcut_clusters\": %zu,\n",
                    shortcut.planes_tested, shortcut.section_face_incidents,
                    shortcut.curves_tested,
                    shortcut.simple_curves, shortcut.fallback_curves,
                    shortcut.nhit, shortcut.ncluster);
            fprintf(j, "  \"shortcut_persistent_clusters\": %zu, \"shortcut_persistent_hits\": %zu, \"shortcut_max_delta_w\": %.4f, \"shortcut_persistent_max_delta_w\": %.4f, \"shortcut_persistent_excess\": %.4f,\n",
                    shortcut.persistent_clusters, shortcut.persistent_hits,
                    shortcut.max_delta_w, shortcut.persistent_max_delta_w,
                    shortcut.persistent_excess);
            fprintf(j, "  \"repair_requested\": %s, \"repair_accepted\": %s, \"repair_step_accepted\": %s, \"repair_output_written\": %s,\n",
                    repair_requested ? "true" : "false",
                    repair.accepted ? "true" : "false",
                    repair.step_accepted ? "true" : "false",
                    repair.output_written ? "true" : "false");
            fprintf(j, "  \"repair_chain_policy\": \"%s\",\n",
                    P->repair_cover_all_supported ? "cover_all_supported" : "best_per_track");
            fprintf(j, "  \"repair_min_chain_support\": %zu,\n",
                    P->repair_min_chain_support);
            fprintf(j, "  \"repair_auto_offset\": %s, \"repair_auto_trials\": %d, \"repair_auto_decisive\": %s, \"repair_auto_suppressed\": %s, \"repair_max_auto_tracks\": %zu,\n",
                    P->repair_auto_offset ? "true" : "false",
                    repair.auto_trials,
                    repair.auto_decisive ? "true" : "false",
                    repair.auto_suppressed ? "true" : "false",
                    P->repair_max_auto_tracks);
            fprintf(j, "  \"repair_phase_offset\": %.6f,\n",
                    repair.phase_offset);
            fprintf(j, "  \"repair_candidate_edges\": %zu, \"repair_selected_chains\": %zu, \"repair_cut_edges\": %zu, \"repair_vertices_before\": %zu, \"repair_vertices_after\": %zu,\n",
                    repair.candidate_edges, repair.selected_chains,
                    repair.cut_edges, repair.vertices_before,
                    repair.vertices_after);
            fprintf(j, "  \"repair_selected_cut_components\": %zu, \"repair_selected_chain_endpoints\": %zu, \"repair_selected_chain_boundary_endpoints\": %zu, \"repair_selected_chain_branch_vertices\": %zu, \"repair_selected_chain_branch_excess\": %zu, \"repair_cut_rank\": %zu,\n",
                    repair.selected_cut_components,
                    repair.selected_chain_endpoints,
                    repair.selected_chain_boundary_endpoints,
                    repair.selected_chain_branch_vertices,
                    repair.selected_chain_branch_excess,
                    repair.cut_rank);
            fprintf(j, "  \"repair_components_before\": %zu, \"repair_components_after\": %zu, \"repair_small_components_before\": %zu, \"repair_small_components_after\": %zu, \"repair_post_shortcut_clusters\": %zu, \"repair_post_shortcut_hits\": %zu, \"repair_post_shortcut_excess\": %.4f, \"repair_post_shortcut_max_delta_w\": %.4f,\n",
                    repair.before.components, repair.after.components,
                    repair.before.small_components,
                    repair.after.small_components,
                    repair.post.persistent_clusters,
                    repair.post.persistent_hits,
                    repair.post.persistent_excess,
                    repair.post.persistent_max_delta_w);
            fprintf(j, "  \"repair_reason\": \"%s\",\n", repair.reason);
            fprintf(j, "  \"top_components\": [\n");
            int jr = (rows < 8) ? rows : 8;
            for (int r = 0; r < jr; r++)
                fprintf(j, "    {\"rank\": %d, \"faces\": %zu, \"vertices\": %zu, \"edges\": %zu, \"boundary_edges\": %zu, \"boundary_loops\": %zu, \"nonmanifold_edges\": %zu, \"boundary_irregular_vertices\": %zu, \"topology_valid\": %s, \"euler_chi\": %lld, \"genus\": %.3f, \"turns\": %.1f, \"phase_span\": %.3f, \"ft_fusions\": %zu}%s\n",
                        r, comp[r].faces, comp[r].vertices, comp[r].edges,
                        comp[r].boundary_edges, comp[r].boundary_loops,
                        comp[r].nonmanifold_edges,
                        comp[r].boundary_irregular_vertices,
                        (comp[r].nonmanifold_edges == 0 &&
                         comp[r].boundary_irregular_vertices == 0) ? "true" : "false",
                        (long long)comp[r].euler_chi, comp[r].genus,
                        PitchTable_dturns(P->ptab, comp[r].rmin,
                                          comp[r].rmax, P->pitch),
                        comp[r].outer_vertices ? comp[r].phase_max - comp[r].phase_min : 0.0,
                        comp[r].merger_edges, (r+1<jr)?",":"");
            fprintf(j, "  ],\n  \"top_shortcut_clusters\": [\n");
            size_t jq = shortcut.ncluster < 8 ? shortcut.ncluster : 8;
            for (size_t r = 0; r < jq; r++) {
                const ShortcutCluster *g = &shortcut.cluster[r];
                fprintf(j, "    {\"rank\": %zu, \"persistent\": %s, \"boundary_supported\": %s, \"boundary_hits\": %zu, \"planes\": %zu, \"hits\": %zu, \"max_delta_w\": %.4f, \"min_distance\": %.3f, \"z_span\": %.3f, \"excess\": %.4f, \"centroid_zyx\": [%.3f, %.3f, %.3f], \"bbox_min_zyx\": [%.3f, %.3f, %.3f], \"bbox_max_zyx\": [%.3f, %.3f, %.3f]}%s\n",
                        r, g->persistent ? "true" : "false",
                        g->boundary_supported ? "true" : "false",
                        g->boundary_hits, g->planes, g->hits, g->max_delta_w,
                        g->min_distance, g->z_max - g->z_min, g->excess,
                        g->centroid[0], g->centroid[1], g->centroid[2],
                        g->bbox_min[0], g->bbox_min[1], g->bbox_min[2],
                        g->bbox_max[0], g->bbox_max[1], g->bbox_max[2],
                        (r + 1 < jq) ? "," : "");
            }
            fprintf(j, "  ],\n  \"shortcut_tracks\": [\n");
            size_t jt = 0;
            for (size_t r = 0; r < shortcut.ncluster; r++) {
                const ShortcutCluster *g = &shortcut.cluster[r];
                if (!g->persistent) continue;
                fprintf(j, "    {\"rank\": %zu, \"boundary_supported\": %s, \"boundary_hits\": %zu, \"planes\": %zu, \"hits\": %zu, \"max_delta_w\": %.4f, \"min_distance\": %.3f, \"z_span\": %.3f, \"excess\": %.4f, \"centroid_zyx\": [%.3f, %.3f, %.3f], \"bbox_min_zyx\": [%.3f, %.3f, %.3f], \"bbox_max_zyx\": [%.3f, %.3f, %.3f]}%s\n",
                        r, g->boundary_supported ? "true" : "false",
                        g->boundary_hits, g->planes, g->hits, g->max_delta_w,
                        g->min_distance, g->z_max - g->z_min, g->excess,
                        g->centroid[0], g->centroid[1], g->centroid[2],
                        g->bbox_min[0], g->bbox_min[1], g->bbox_min[2],
                        g->bbox_max[0], g->bbox_max[1], g->bbox_max[2],
                        (jt + 1 < shortcut.persistent_clusters) ? "," : "");
                jt++;
            }
            fprintf(j, "  ]\n}\n");
            int write_failed = ferror(j);
            if (fclose(j) != 0) write_failed = 1;
            if (write_failed) {
                fprintf(stderr, "wind_audit: failed writing json: %s\n",
                        P->json_path);
                artifact_failed = 1;
            } else {
                printf("  wrote json -> %s\n", P->json_path);
            }
        }
    }

    /* ---- optional dumps ---- */
    if (P->dump_mergers && face_flag) {
        FILE *o = fopen(P->dump_mergers, "wb");
        if (o == NULL) {
            fprintf(stderr, "wind_audit: cannot write mergers overlay: %s\n",
                    P->dump_mergers);
            artifact_failed = 1;
        } else {
            fprintf(o, "# wind_audit mergers: flagged faces red, rest grey\n");
            for (size_t i = 0; i < nv; i++)
                fprintf(o, "v %.5g %.5g %.5g 0.30 0.30 0.30\n", V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
            /* re-emit: flagged faces get bright verts via a second pass would need
             * per-vertex; simplest: emit flagged faces with a marker colour by
             * writing them last as separate verts */
            for (size_t f = 0; f < nf; f++)
                fprintf(o, "f %d %d %d\n", F.f[f*3+0]+1, F.f[f*3+1]+1, F.f[f*3+2]+1);
            /* overlay flagged faces as bright-red duplicated triangles */
            size_t base = nv;
            for (size_t f = 0; f < nf; f++) if (face_flag[f]) {
                for (int k = 0; k < 3; k++) { int32_t vi = F.f[f*3+k];
                    fprintf(o, "v %.5g %.5g %.5g 1.0 0.05 0.05\n", V.v[vi*3+0], V.v[vi*3+1], V.v[vi*3+2]); }
            }
            size_t idx = base;
            for (size_t f = 0; f < nf; f++) if (face_flag[f]) { fprintf(o, "f %zu %zu %zu\n", idx+1, idx+2, idx+3); idx += 3; }
            int write_failed = ferror(o);
            if (fclose(o) != 0) write_failed = 1;
            if (write_failed) {
                fprintf(stderr,
                        "wind_audit: failed writing mergers overlay: %s\n",
                        P->dump_mergers);
                artifact_failed = 1;
            } else {
                printf("\n  wrote mergers overlay -> %s\n",
                       P->dump_mergers);
            }
        }
    }
    if (P->dump_splits) {
        if (!P->split_audit) {
            fprintf(stderr,
                    "wind_audit: --dump-splits requires split auditing\n");
            artifact_failed = 1;
        } else {
            FILE *o = fopen(P->dump_splits, "wb");
            if (o == NULL) {
                fprintf(stderr,
                        "wind_audit: cannot write splits overlay: %s\n",
                        P->dump_splits);
                artifact_failed = 1;
            } else {
                fprintf(o, "# wind_audit splits: full mesh grey, split-pair verts bright, l-lines across gaps\n");
                for (size_t i = 0; i < nv; i++)
                    fprintf(o, "v %.5g %.5g %.5g 0.28 0.28 0.28\n", V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
                for (size_t f = 0; f < nf; f++)
                    fprintf(o, "f %d %d %d\n", F.f[f*3+0]+1, F.f[f*3+1]+1, F.f[f*3+2]+1);
                /* bright endpoints + connecting lines (scan agg: list[].n was
                 * consumed by the top-K report above) */
                size_t idx = nv;
                for (size_t s = 0; s < aslots; s++) if (agg[s].ra != -1) {
                    size_t u = agg[s].rep_u, v = agg[s].rep_v;
                    fprintf(o, "v %.5g %.5g %.5g 0.1 1.0 0.1\n", V.v[u*3+0], V.v[u*3+1], V.v[u*3+2]);
                    fprintf(o, "v %.5g %.5g %.5g 1.0 0.9 0.1\n", V.v[v*3+0], V.v[v*3+1], V.v[v*3+2]);
                    fprintf(o, "l %zu %zu\n", idx+1, idx+2); idx += 2;
                }
                int write_failed = ferror(o);
                if (fclose(o) != 0) write_failed = 1;
                if (write_failed) {
                    fprintf(stderr,
                            "wind_audit: failed writing splits overlay: %s\n",
                            P->dump_splits);
                    artifact_failed = 1;
                } else {
                    printf("  wrote splits overlay -> %s\n",
                           P->dump_splits);
                }
            }
        }
    }
    if (P->dump_shortcuts) {
        if (P->shortcut_tol <= 0.0) {
            fprintf(stderr,
                    "wind_audit: --dump-shortcuts requires shortcut auditing\n");
            artifact_failed = 1;
        } else {
            FILE *o = fopen(P->dump_shortcuts, "wb");
            if (o == NULL) {
                fprintf(stderr,
                        "wind_audit: cannot write shortcut overlay: %s\n",
                        P->dump_shortcuts);
                artifact_failed = 1;
            } else {
                fprintf(o, "# wind_audit persistent local shortcuts: mesh grey; endpoint pairs red/yellow\n");
                for (size_t i = 0; i < nv; i++)
                    fprintf(o, "v %.5g %.5g %.5g 0.24 0.24 0.24\n",
                            V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
                for (size_t f = 0; f < nf; f++)
                    fprintf(o, "f %d %d %d\n", F.f[f*3+0]+1,
                            F.f[f*3+1]+1, F.f[f*3+2]+1);
                size_t idx = nv;
                for (size_t i = 0; i < shortcut.nhit; i++) {
                    const ShortcutHit *h = &shortcut.hit[i];
                    if (!h->persistent) continue;
                    write_octahedron_marker(o, h->a, 0.9,
                                            1.0, 0.05, 0.05, &idx);
                    write_octahedron_marker(o, h->b, 0.9,
                                            1.0, 0.85, 0.05, &idx);
                }
                int write_failed = ferror(o);
                if (fclose(o) != 0) write_failed = 1;
                if (write_failed) {
                    fprintf(stderr,
                            "wind_audit: failed writing shortcut overlay: %s\n",
                            P->dump_shortcuts);
                    artifact_failed = 1;
                } else {
                    printf("  wrote shortcut overlay -> %s\n",
                           P->dump_shortcuts);
                }
            }
        }
    }

    size_t gate_shortcuts = shortcut.persistent_clusters;
    if (repair_requested && (repair.accepted || repair.step_accepted))
        gate_shortcuts = repair.post.persistent_clusters;
    int repair_rejected = repair_requested &&
                          !repair.accepted && !repair.step_accepted;
    free(list); free(agg); free(H.cell); free(H.next);
    free(face_flag);
    free(vrank); free(par); free(usz); free(root_faces); free(comp); free(root2rank);
    free(V.v); free(F.f);
    int local_failed = merge_ft > 0 || broad_fusion_components > 0 ||
                       (P->shortcut_tol > 0.0 && gate_shortcuts > 0);
    int shortcut_failed = P->fail_on_shortcut &&
                          gate_shortcuts > 0;
    shortcut_audit_free(&repair.post);
    shortcut_audit_free(&shortcut);
    if (shortcut_rc != 0 || repair_rc != 0 || artifact_failed) return 1;
    if (repair_rejected) return 5;
    return ((P->fail_on_local && local_failed) || shortcut_failed) ? 4 : 0;
}

/* ---------- selftest ---------- */
static int approx(double a, double b, double e) { return fabs(a - b) <= e; }

static int selftest(void)
{
    int fails = 0;

    /* A branchy base must search for a regular value; a simple embedded
     * admissible cut does not need a phase sweep. */
    {
        Params P; ShortcutAudit before; ShortcutRepair r;
        memset(&P, 0, sizeof(P)); memset(&before, 0, sizeof(before));
        memset(&r, 0, sizeof(r));
        P.repair_max_cut_fraction = 0.015;
        before.persistent_clusters = 100;
        before.persistent_hits = 1000;
        before.persistent_excess = 100.0;
        r.cut_edges = 500; r.cut_rank = 100;
        r.before.components = 10; r.after.components = 20;
        r.before.small_components = 2; r.after.small_components = 7;
        r.post.persistent_clusters = 40;
        r.post.persistent_hits = 400;
        r.post.persistent_excess = 40.0;
        r.post.persistent_max_delta_w = 0.8;
        before.persistent_max_delta_w = 1.0;
        r.selected_chains = 20;
        r.selected_chain_branch_vertices = 3;
        if (!repair_needs_offset_search(&r, &before, &P, 100000)) {
            fprintf(stderr, "selftest: branchy base skipped phase sweep\n");
            fails++;
        }
        r.selected_chain_branch_vertices = 0;
        if (repair_trial_admissible(&r, &before, &P, 100000)) {
            fprintf(stderr, "selftest: small-fragment-producing cut accepted\n");
            fails++;
        }
        r.after.small_components = r.before.small_components;
        if (repair_needs_offset_search(&r, &before, &P, 100000)) {
            fprintf(stderr, "selftest: branch-free admissible base searched\n");
            fails++;
        }

        P.repair_monotone_step = 1;
        r.after.small_components = r.before.small_components;
        r.post.persistent_clusters = 40;
        r.post.persistent_hits = 400;
        r.post.persistent_excess = 40.0;
        r.post.persistent_max_delta_w = 0.8;
        if (!repair_trial_monotone(&r, &before, &P, 100000)) {
            fprintf(stderr, "selftest: safe monotone repair step rejected\n");
            fails++;
        }
        r.post.persistent_hits = 1001;
        if (repair_trial_monotone(&r, &before, &P, 100000)) {
            fprintf(stderr, "selftest: rising-hit repair step accepted\n");
            fails++;
        }
        r.post.persistent_hits = 400;
        P.repair_max_auto_tracks = 4096;
        before.persistent_clusters = 4097;
        if (repair_auto_search_allowed(&P, &before)) {
            fprintf(stderr, "selftest: oversized phase sweep allowed\n");
            fails++;
        }
        before.persistent_clusters = 4096;
        if (!repair_auto_search_allowed(&P, &before)) {
            fprintf(stderr, "selftest: bounded phase sweep suppressed\n");
            fails++;
        }
        {
            size_t parsed = 1;
            if (!parse_size_arg("4096", &parsed) || parsed != 4096 ||
                !parse_size_arg("0", &parsed) || parsed != 0 ||
                parse_size_arg("-1", &parsed) ||
                parse_size_arg("unlimited", &parsed)) {
                fprintf(stderr, "selftest: unsafe auto-track limit parsing\n");
                fails++;
            }
        }
    }

    /* The section incidence accelerator must be bit-for-bit equivalent to the
     * strict brute-force crossing predicate, including faces whose extrema lie
     * exactly on a sampled plane. */
    {
        float vv[] = {
            -1,0,0,  0,0,0,  1,0,0,
             .5f,0,0,  1,0,0,  2.5f,0,0,
            -2,0,0,  5,0,0,  6,0,0,
             2.5f,0,0,  3,0,0,  4.5f,0,0,
             .49f,0,0, .51f,0,0, .5f,0,0
        };
        int32_t ff[] = {
             0, 1, 2,  3, 4, 5,  6, 7, 8,
             9,10,11, 12,13,14
        };
        FVec Vt = { vv, 15, 15 };
        IVec Ft = { ff, 5, 5 };
        SectionFaceIndex index;
        if (section_face_index_build(&Vt, &Ft, .5, 2.0, 6.0,
                                     &index) != 0) {
            fprintf(stderr, "selftest: section face index build failed\n");
            fails++;
        } else {
            for (size_t p = 0; p < index.nplane; p++) {
                uint8_t seen[5] = {0};
                for (size_t at = index.row[p]; at < index.row[p + 1]; at++) {
                    size_t f = index.face[at];
                    if (f >= Ft.n || seen[f]) {
                        fprintf(stderr,
                                "selftest: section face index duplicate/range\n");
                        fails++;
                        break;
                    }
                    seen[f] = 1;
                }
                for (size_t f = 0; f < Ft.n; f++) {
                    double lo = DBL_MAX, hi = -DBL_MAX;
                    for (int k = 0; k < 3; k++) {
                        int32_t v = Ft.f[f*3+(size_t)k];
                        double z = Vt.v[(size_t)v*3];
                        if (z < lo) lo = z;
                        if (z > hi) hi = z;
                    }
                    int expected = lo < index.plane[p] && hi > index.plane[p];
                    if (seen[f] != (uint8_t)expected) {
                        fprintf(stderr,
                                "selftest: section index p=%zu f=%zu got=%u want=%d\n",
                                p, f, (unsigned)seen[f], expected);
                        fails++;
                    }
                }
            }
            section_face_index_free(&index);
        }
    }

    /* winding math: a full-pitch radial jump at constant angle -> dw ~ 1 */
    {
        float Vt[6] = { 0, 50.0f, 0, 0, 59.5f, 0 };  /* (z,y,x): r 50 -> 59.5, same angle */
        double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5, NULL, NULL);
        if (!approx(dw, 1.0, 0.02)) { fprintf(stderr, "selftest: radial-jump dw=%.4f want ~1.0\n", dw); fails++; }
    }
    /* Large-region phase must remain continuous across atan2 and beyond a
     * half-turn, while a true radial full-pitch excursion stays visible. */
    {
        const int columns = 81;
        FVec Vt = {0}; IVec Ft = {0};
        Params P; Comp c; memset(&P, 0, sizeof(P)); memset(&c, 0, sizeof(c));
        P.pitch = 9.5;
        for (int row = 0; row < 2; row++) {
            for (int col = 0; col < columns; col++) {
                double theta = -0.75*TWO_PI + 1.5*TWO_PI*(double)col/(double)(columns-1);
                double radius = P.pitch*(6.0 + theta/TWO_PI);
                if (fv_push(&Vt, (float)row, (float)(radius*sin(theta)),
                            (float)(radius*cos(theta))) != 0) fails++;
            }
        }
        for (int col = 0; col + 1 < columns; col++) {
            int32_t a=col, b=col+1, d=columns+col, e=columns+col+1;
            if (iv_push(&Ft,a,b,d) != 0 || iv_push(&Ft,b,e,d) != 0) fails++;
        }
        int32_t *rank = (int32_t *)malloc(Vt.n*sizeof(*rank));
        if (rank == NULL) fails++;
        else {
            for (size_t i = 0; i < Vt.n; i++) rank[i] = 0;
            if (broad_phase_spans(&Vt,&Ft,&P,rank,&c,1,20.0) != 0 ||
                c.phase_max-c.phase_min > 2e-3) {
                fprintf(stderr,
                        "selftest: large valid spiral phase span %.6f want ~0\n",
                        c.phase_max-c.phase_min);
                fails++;
            }
            int mid = columns/2;
            for (int row = 0; row < 2; row++) {
                size_t v = (size_t)(row*columns+mid);
                double theta = 0.0, radius = P.pitch*7.0;
                Vt.v[v*3+1]=(float)(radius*sin(theta));
                Vt.v[v*3+2]=(float)(radius*cos(theta));
            }
            if (broad_phase_spans(&Vt,&Ft,&P,rank,&c,1,20.0) != 0 ||
                c.phase_max-c.phase_min < 0.9) {
                fprintf(stderr,
                        "selftest: radial excursion phase span %.6f want >=0.9\n",
                        c.phase_max-c.phase_min);
                fails++;
            }
            free(rank);
        }
        free(Vt.v); free(Ft.f);
    }
    /* winding math: same turn, small angular step -> dw ~ 0 */
    {
        /* two points on r=50 circle 0.1 rad apart: r const, |dth|=0.1 -> |dw|=0.1/2pi ~ 0.0159 */
        float Vt[6] = { 0, 50.0f, 0, 0, (float)(50*cos(0.1)), (float)(50*sin(0.1)) };
        double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5, NULL, NULL);
        if (!approx(fabs(dw), 0.1/TWO_PI, 0.01)) { fprintf(stderr, "selftest: tangential |dw|=%.4f want ~%.4f\n", fabs(dw), 0.1/TWO_PI); fails++; }
    }
    /* branch-cut: crossing theta=0 ray on the same circle stays ~0 */
    {
        float Vt[6] = { 0, (float)(50*sin(-0.05)), (float)(50*cos(-0.05)),
                        0, (float)(50*sin(+0.05)), (float)(50*cos(+0.05)) };
        double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5, NULL, NULL);
        if (fabs(dw) > 0.02) { fprintf(stderr, "selftest: branch-cut dw=%.4f want ~0\n", dw); fails++; }
    }

    /* Curved-axis regression: input rows may be unordered; interpolation and
     * endpoint extrapolation must preserve a fixed local polar coordinate as
     * the umbilicus moves with z. A static axis would report a false jump. */
    {
        const char *tin = "wind_audit_selftest_axis.csv";
        FILE *fp = fopen(tin, "wb");
        AxisTable axis;
        memset(&axis, 0, sizeof(axis));
        if (fp == NULL) {
            fprintf(stderr, "selftest: axis tmp write fail\n");
            return 3;
        }
        fprintf(fp, "z_l0,y_l0,x_l0\n20,40,-20\n0,0,0\n10,20,-10\n");
        fclose(fp);
        if (axis_table_load(tin, &axis) != 0) {
            fprintf(stderr, "selftest: axis table load fail\n");
            fails++;
        } else {
            double y, x;
            axis_table_eval(&axis, 5.0, 999.0, 999.0, &y, &x);
            if (!approx(y, 10.0, 1e-12) || !approx(x, -5.0, 1e-12)) {
                fprintf(stderr, "selftest: axis interpolation=(%.4f,%.4f)\n", y, x);
                fails++;
            }
            axis_table_eval(&axis, -5.0, 999.0, 999.0, &y, &x);
            if (!approx(y, -10.0, 1e-12) || !approx(x, 5.0, 1e-12)) {
                fprintf(stderr, "selftest: axis low extrapolation=(%.4f,%.4f)\n", y, x);
                fails++;
            }
            axis_table_eval(&axis, 25.0, 999.0, 999.0, &y, &x);
            if (!approx(y, 50.0, 1e-12) || !approx(x, -25.0, 1e-12)) {
                fprintf(stderr, "selftest: axis high extrapolation=(%.4f,%.4f)\n", y, x);
                fails++;
            }
            {
                float Vt[6] = { 0, 50, 0, 10, 70, -10 };
                double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5, &axis, NULL);
                if (!approx(dw, 0.0, 1e-6)) {
                    fprintf(stderr, "selftest: moving-axis dw=%.6f want 0\n", dw);
                    fails++;
                }
            }
        }
        axis_table_free(&axis);
        remove(tin);
    }

    /* end-to-end: build a tiny mesh with a known merger + a known split,
     * write it, run(), and check the printed-path invariants via return code +
     * internal recomputation. We verify detection logic directly here. */
    {
        /* Case A: two concentric ring-arcs (r=50, r=57.6) = ADJACENT turns.
         *   dw ~ 7.6/9.5 = 0.80 -> NOT a split, NOT... they are separate comps
         *   with a ~7.6 gap. With split_wtol 0.35 they must NOT be flagged. */
        double dw_adj = 7.6 / 9.5;
        if (dw_adj <= 0.35) { fprintf(stderr, "selftest: adjacent-turn dw=%.3f should exceed split_wtol\n", dw_adj); fails++; }

        /* Case B: same ring severed -> two arcs at r=50, angular gap; dw ~ 0,
         *   physical gap ~ small -> MUST be flagged. Represented by two verts
         *   at r=50, 3 vox apart tangentially. */
        float p0y = 50.0f, p0x = 0.0f;
        float p1y = (float)(50*cos(3.0/50)), p1x = (float)(50*sin(3.0/50)); /* ~3 vox arc */
        float Vt[6] = { 0, p0y, p0x, 0, p1y, p1x };
        double dw_sev = dwind(Vt, 0, 1, 0.0, 0.0, 9.5, NULL, NULL);
        double gap = hypot((double)p1y-p0y, (double)p1x-p0x);
        if (fabs(dw_sev) > 0.35) { fprintf(stderr, "selftest: severed dw=%.3f should be < split_wtol\n", dw_sev); fails++; }
        if (gap > 14.0) { fprintf(stderr, "selftest: severed gap=%.2f should be < split_gap\n", gap); fails++; }
    }
    /* A staircase fusion can span a full turn while every constituent edge is
     * below the abrupt-merger threshold. This is exactly why the component
     * phase-span gate exists. */
    {
        float Vt[5 * 3] = {
            0, 50.000f, 0,  0, 52.375f, 0,  0, 54.750f, 0,
            0, 57.125f, 0,  0, 59.500f, 0
        };
        double pmin = 1e30, pmax = -1e30;
        for (size_t i = 0; i < 5; i++) {
            double ph = vradius(Vt, i, 0.0, 0.0, NULL) / 9.5;
            if (ph < pmin) pmin = ph;
            if (ph > pmax) pmax = ph;
            if (i > 0 && fabs(dwind(Vt, i-1, i, 0.0, 0.0, 9.5, NULL, NULL)) >= 0.40) {
                fprintf(stderr, "selftest: staircase edge %zu looked abrupt\n", i);
                fails++;
            }
        }
        if (!approx(pmax - pmin, 1.0, 0.02)) {
            fprintf(stderr, "selftest: staircase phase span %.3f want 1.0\n", pmax-pmin);
            fails++;
        }
    }

    /* The new local metric must catch a gradual one-turn wall even though its
     * five constituent edges each advance only 0.2 turn. The same topology
     * laid along an Archimedean spiral has constant winding and must pass. */
    {
        float verts[42 * 3];
        int32_t faces[40 * 3];
        for (int i = 0; i < 6; i++) {
            double r = 50.0 + 9.5 * (double)i / 5.0;
            verts[(size_t)(2*i)*3+0] = 0.0f;
            verts[(size_t)(2*i)*3+1] = (float)r;
            verts[(size_t)(2*i)*3+2] = 0.0f;
            verts[(size_t)(2*i+1)*3+0] = 4.0f;
            verts[(size_t)(2*i+1)*3+1] = (float)r;
            verts[(size_t)(2*i+1)*3+2] = 0.0f;
        }
        for (int i = 0; i < 5; i++) {
            int32_t a = 2*i, b = a+1, c = a+2, d = a+3;
            faces[(size_t)(2*i)*3+0] = a;
            faces[(size_t)(2*i)*3+1] = c;
            faces[(size_t)(2*i)*3+2] = b;
            faces[(size_t)(2*i+1)*3+0] = b;
            faces[(size_t)(2*i+1)*3+1] = c;
            faces[(size_t)(2*i+1)*3+2] = d;
        }
        FVec Vt = { verts, 12, 42 };
        IVec Ft = { faces, 10, 40 };
        Params P; memset(&P, 0, sizeof(P));
        P.pitch = 9.5; P.core_pitches = 3.0;
        P.shortcut_tol = 0.55; P.shortcut_length_pitches = 2.0;
        P.shortcut_z_step = 1.0; P.shortcut_cluster_pitches = 1.5;
        P.shortcut_min_points = 4; P.shortcut_min_planes = 3;
        ShortcutAudit A;
        if (shortcut_audit(&Vt, &Ft, &P, &A) != 0 ||
            A.persistent_clusters != 1 || A.persistent_max_delta_w < 0.95) {
            fprintf(stderr,
                    "selftest: local shortcut wall tracks=%zu max=%.3f want 1/~1\n",
                    A.persistent_clusters, A.persistent_max_delta_w);
            fails++;
        }
        P.repair_path = "wind_audit_selftest_repair.obj";
        P.repair_radius_pitches = 1.5;
        P.repair_support_radius_pitches = 1.0;
        P.repair_band_margin = 0.12;
        P.repair_max_cut_fraction = 0.5;
        P.min_faces = 2;
        ShortcutRepair R;
        if (shortcut_seam_repair(&Vt, &Ft, &P, &A, &R) != 0 ||
            !R.accepted || !R.output_written || R.cut_edges == 0 ||
            R.post.persistent_clusters != 0 ||
            R.before.components != 1 || R.after.components != 2) {
            fprintf(stderr,
                    "selftest: seam repair accepted=%d cut=%zu post=%zu comps=%zu->%zu\n",
                    R.accepted, R.cut_edges, R.post.persistent_clusters,
                    R.before.components, R.after.components);
            fails++;
        }
        shortcut_audit_free(&R.post);
        remove(P.repair_path);
        P.repair_path = NULL;
        shortcut_audit_free(&A);

        /* A clipped bridge may occupy only two sampled planes. It remains a
         * hard failure when its open section contour terminates on a known
         * cube boundary; without that evidence, the same two-plane event is
         * only a warning. */
        for (int i = 0; i < 6; i++)
            verts[(size_t)(2*i+1)*3+0] = 2.2f;
        if (shortcut_audit(&Vt, &Ft, &P, &A) != 0 ||
            A.persistent_clusters != 0 || A.nhit == 0) {
            fprintf(stderr,
                    "selftest: isolated two-plane event tracks=%zu hits=%zu want 0/>0\n",
                    A.persistent_clusters, A.nhit);
            fails++;
        }
        shortcut_audit_free(&A);
        P.seam_pitch = 128.0;
        P.seam_tol = 3.0;
        if (shortcut_audit(&Vt, &Ft, &P, &A) != 0 ||
            A.persistent_clusters != 1 || A.ncluster == 0 ||
            !A.cluster[0].boundary_supported) {
            fprintf(stderr,
                    "selftest: clipped two-plane event tracks=%zu boundary=%d want 1/1\n",
                    A.persistent_clusters,
                    A.ncluster ? A.cluster[0].boundary_supported : 0);
            fails++;
        }
        shortcut_audit_free(&A);
        P.boundary_box_set = 1;
        P.boundary_box[0] = 10.0; P.boundary_box[1] = 20.0;
        P.boundary_box[2] = -100.0; P.boundary_box[3] = 100.0;
        P.boundary_box[4] = -100.0; P.boundary_box[5] = 100.0;
        if (shortcut_audit(&Vt, &Ft, &P, &A) != 0 ||
            A.persistent_clusters != 0 || A.nhit == 0) {
            fprintf(stderr,
                    "selftest: internal seam event tracks=%zu hits=%zu want 0/>0\n",
                    A.persistent_clusters, A.nhit);
            fails++;
        }
        shortcut_audit_free(&A);
        P.boundary_box_set = 0;

        /* A valid clipped section need not make a complete turn. This open
         * 0.38-turn Archimedean segment crosses atan2's branch cut and lies
         * on a cube boundary, but its local winding coordinate is constant. */
        P.shortcut_length_pitches = 6.0;
        for (int i = 0; i < 21; i++) {
            double theta = 2.50 + 0.12 * (double)i;
            double r = 35.0 + 9.5 * theta / TWO_PI;
            float y = (float)(r * sin(theta));
            float x = (float)(r * cos(theta));
            verts[(size_t)(2*i)*3+0] = 0.0f;
            verts[(size_t)(2*i+1)*3+0] = 4.0f;
            verts[(size_t)(2*i)*3+1] = y;
            verts[(size_t)(2*i)*3+2] = x;
            verts[(size_t)(2*i+1)*3+1] = y;
            verts[(size_t)(2*i+1)*3+2] = x;
        }
        for (int i = 0; i < 20; i++) {
            int32_t a = 2*i, b = a+1, c = a+2, d = a+3;
            faces[(size_t)(2*i)*3+0] = a;
            faces[(size_t)(2*i)*3+1] = c;
            faces[(size_t)(2*i)*3+2] = b;
            faces[(size_t)(2*i+1)*3+0] = b;
            faces[(size_t)(2*i+1)*3+1] = c;
            faces[(size_t)(2*i+1)*3+2] = d;
        }
        Vt.n = 42;
        Ft.n = 40;
        if (shortcut_audit(&Vt, &Ft, &P, &A) != 0 ||
            A.persistent_clusters != 0) {
            fprintf(stderr,
                    "selftest: valid spiral shortcut tracks=%zu want 0\n",
                    A.persistent_clusters);
            fails++;
        }
        shortcut_audit_free(&A);
    }

    /* Severity sorting must preserve the hit-to-track association used by
     * repair. Put a 0.75-turn wall first and a 1.0-turn wall second so qsort
     * reverses their construction order, then verify every hit still lies in
     * the bbox of the cluster named by hit.cluster_id. */
    {
        float verts[24 * 3];
        int32_t faces[20 * 3];
        for (int wall = 0; wall < 2; wall++) {
            double base = wall == 0 ? 50.0 : 100.0;
            double span = wall == 0 ? 0.75 * 9.5 : 9.5;
            int vo = wall * 12, fo = wall * 10;
            for (int i = 0; i < 6; i++) {
                double r = base + span * (double)i / 5.0;
                verts[(size_t)(vo+2*i)*3+0] = 0.0f;
                verts[(size_t)(vo+2*i+1)*3+0] = 4.0f;
                verts[(size_t)(vo+2*i)*3+1] = wall == 0 ? (float)r : 0.0f;
                verts[(size_t)(vo+2*i+1)*3+1] = wall == 0 ? (float)r : 0.0f;
                verts[(size_t)(vo+2*i)*3+2] = wall == 0 ? 0.0f : (float)r;
                verts[(size_t)(vo+2*i+1)*3+2] = wall == 0 ? 0.0f : (float)r;
            }
            for (int i = 0; i < 5; i++) {
                int32_t a = vo+2*i, b = a+1, c = a+2, d = a+3;
                faces[(size_t)(fo+2*i)*3+0] = a;
                faces[(size_t)(fo+2*i)*3+1] = c;
                faces[(size_t)(fo+2*i)*3+2] = b;
                faces[(size_t)(fo+2*i+1)*3+0] = b;
                faces[(size_t)(fo+2*i+1)*3+1] = c;
                faces[(size_t)(fo+2*i+1)*3+2] = d;
            }
        }
        FVec Vt = { verts, 24, 24 };
        IVec Ft = { faces, 20, 20 };
        Params P; memset(&P, 0, sizeof(P));
        P.pitch = 9.5; P.core_pitches = 3.0;
        P.shortcut_tol = 0.55; P.shortcut_length_pitches = 2.0;
        P.shortcut_z_step = 1.0; P.shortcut_cluster_pitches = 1.5;
        P.shortcut_min_points = 4; P.shortcut_min_planes = 3;
        ShortcutAudit A;
        int bad_map = 0;
        if (shortcut_audit(&Vt, &Ft, &P, &A) != 0 ||
            A.persistent_clusters != 2) {
            fprintf(stderr,
                    "selftest: sorted track map tracks=%zu want 2\n",
                    A.persistent_clusters);
            fails++;
        } else {
            for (size_t i = 0; i < A.nhit; i++) {
                int32_t ci = A.hit[i].cluster_id;
                if (ci < 0 || (size_t)ci >= A.ncluster) {
                    bad_map = 1;
                    break;
                }
                const ShortcutCluster *g = &A.cluster[(size_t)ci];
                for (int axis = 0; axis < 3; axis++)
                    if (A.hit[i].p[axis] < g->bbox_min[axis] - 1e-6 ||
                        A.hit[i].p[axis] > g->bbox_max[axis] + 1e-6)
                        bad_map = 1;
            }
            if (bad_map) {
                fprintf(stderr,
                        "selftest: severity sort broke hit-to-track mapping\n");
                fails++;
            }
        }
        shortcut_audit_free(&A);
    }

    /* full run smoke on a written OBJ: 2 disjoint triangles at different radii */
    {
        const char *tin = "wind_audit_selftest_in.obj";
        FILE *fp = fopen(tin, "wb");
        if (!fp) { fprintf(stderr, "selftest: tmp write fail\n"); return 3; }
        /* comp0 near r~50, comp1 near r~200 (far apart, no split, no merger) */
        fprintf(fp, "v 0 50 0\nv 0 51 1\nv 0 50 2\n");
        fprintf(fp, "v 0 200 0\nv 0 201 1\nv 0 200 2\n");
        fprintf(fp, "f 1 2 3\nf 4 5 6\n");
        fclose(fp);
        Params P; memset(&P, 0, sizeof(P));
        P.umb_y = 0; P.umb_x = 0; P.pitch = 9.5; P.merge_tol = 0.40; P.split_wtol = 0.35;
        P.split_gap = 14.0; P.core_pitches = 3.0; P.min_faces = 0; P.top = 5;
        P.split_audit = 1;
        int rc = run(tin, &P);
        remove(tin);
        if (rc != 0) { fprintf(stderr, "selftest: run rc=%d\n", rc); fails++; }
    }

    /* Integrated broad-fusion gate: a connected radial staircase spans one
     * phase turn, but each edge advances only one quarter turn. */
    {
        const char *tin = "wind_audit_selftest_local.obj";
        FILE *fp = fopen(tin, "wb");
        if (!fp) { fprintf(stderr, "selftest: local tmp write fail\n"); return 3; }
        for (int i = 0; i < 5; i++) {
            double r = 50.0 + 2.375 * (double)i;
            fprintf(fp, "v 0 %.6f 0\nv 1 %.6f 0\n", r, r);
        }
        for (int i = 0; i < 4; i++) {
            int a = 2*i + 1, b = a + 1, c = a + 2, d = a + 3;
            fprintf(fp, "f %d %d %d\nf %d %d %d\n", a, c, b, b, c, d);
        }
        fclose(fp);
        Params P; memset(&P, 0, sizeof(P));
        P.umb_y = 0; P.umb_x = 0; P.pitch = 9.5;
        P.merge_tol = 0.40; P.split_wtol = 0.35; P.split_gap = 14.0;
        P.core_pitches = 3.0; P.min_faces = 1; P.top = 5;
        P.local_span_tol = 0.9; P.fail_on_local = 1;
        int rc = run(tin, &P);
        remove(tin);
        if (rc != 4) {
            fprintf(stderr, "selftest: broad-fusion gate rc=%d want 4\n", rc);
            fails++;
        }
    }

    fprintf(stderr, "wind_audit selftest: %s\n", fails ? "FAIL" : "PASS");
    return fails ? 3 : 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--selftest")) {
        int pt = PitchTable_selftest();
        int st = selftest();
        return st != 0 ? st : (pt != 0 ? 3 : 0);
    }
    if (argc < 2) {
        fprintf(stderr,
            "usage: %s <in.obj> [--axis-table z_y_x.csv | --umb-y Y --umb-x X]\n"
            "         [--pitch P=9.5] [--top N=25]\n"
            "         [--merge-tol T=0.40] [--split-wtol T=0.25] [--split-gap G=14] [--no-split-audit]\n"
            "         [--core-pitches C=3] [--min-faces M=64]\n"
            "         [--local-span-tol T] [--fail-on-local]\n"
            "         [--shortcut-tol T] [--shortcut-length-pitches L=2]\n"
            "         [--shortcut-z-step Z=4] [--shortcut-min-planes N=3]\n"
            "         [--shortcut-track-mutual]\n"
            "         [--boundary-box Z0 Z1 Y0 Y1 X0 X1 | --boundary-box-mesh]\n"
            "         [--fail-on-shortcut] [--dump-shortcuts f.obj]\n"
            "         [--repair-seams f.obj] [--repair-radius-pitches R=1.5]\n"
            "         [--repair-dry-run] [--repair-emit-rejected]\n"
            "         [--repair-cover-all-supported]\n"
            "         [--repair-monotone-step]\n"
            "         [--repair-min-chain-support N=1]\n"
            "         [--repair-chain-csv chains.csv]\n"
            "         [--repair-max-cut-fraction F=.05]\n"
            "         [--repair-support-radius-pitches R=1] [--repair-band-margin M=.12]\n"
            "         [--repair-phase-offset W=0] [--repair-auto-offset]\n"
            "         [--repair-max-auto-tracks N=4096; 0=unlimited]\n"
            "         [--dump-mergers f.obj] [--dump-splits f.obj]\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    Params P; AxisTable axis;
    memset(&P, 0, sizeof(P));
    memset(&axis, 0, sizeof(axis));
    /* split_wtol default = the seam gate's SOFT tol (0.25): a split below it is
     * one the gate WOULD have allowed -> a genuine weld failure, not a gate-
     * correct rejection of a ~1/3-turn-apart chart. */
    P.umb_y = 3405; P.umb_x = 2878; P.pitch = 9.5; P.merge_tol = 0.40; P.split_wtol = 0.25;
    P.split_gap = 14.0; P.core_pitches = 3.0; P.min_faces = 64; P.top = 25;
    P.split_audit = 1;
    P.seam_pitch = 128.0; P.seam_tol = 3.0;
    P.shortcut_length_pitches = 2.0; P.shortcut_z_step = 4.0;
    P.shortcut_cluster_pitches = 1.5; P.shortcut_min_points = 4;
    P.shortcut_min_planes = 3;
    P.repair_radius_pitches = 1.5;
    P.repair_support_radius_pitches = 1.0;
    P.repair_band_margin = 0.12;
    P.repair_max_cut_fraction = 0.05;
    P.repair_min_chain_support = 1;
    P.repair_max_auto_tracks = 4096;
    const char *in = argv[1];
    for (int i = 2; i < argc; i++) {
        if      (!strcmp(argv[i], "--axis-table") && i+1 < argc) P.axis_table_path = argv[++i];
        else if (!strcmp(argv[i], "--umb-y") && i+1 < argc) P.umb_y = atof(argv[++i]);
        else if (!strcmp(argv[i], "--umb-x") && i+1 < argc) P.umb_x = atof(argv[++i]);
        else if (!strcmp(argv[i], "--pitch") && i+1 < argc) P.pitch = atof(argv[++i]);
        else if (!strcmp(argv[i], "--pitch-table") && i+1 < argc) P.pitch_table_path = argv[++i];
        else if (!strcmp(argv[i], "--top") && i+1 < argc) P.top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--merge-tol") && i+1 < argc) P.merge_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--split-wtol") && i+1 < argc) P.split_wtol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--split-gap") && i+1 < argc) P.split_gap = atof(argv[++i]);
        else if (!strcmp(argv[i], "--no-split-audit")) P.split_audit = 0;
        else if (!strcmp(argv[i], "--core-pitches") && i+1 < argc) P.core_pitches = atof(argv[++i]);
        else if (!strcmp(argv[i], "--local-span-tol") && i+1 < argc) P.local_span_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--fail-on-local")) P.fail_on_local = 1;
        else if (!strcmp(argv[i], "--shortcut-tol") && i+1 < argc) P.shortcut_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-length-pitches") && i+1 < argc) P.shortcut_length_pitches = atof(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-z-step") && i+1 < argc) P.shortcut_z_step = atof(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-cluster-pitches") && i+1 < argc) P.shortcut_cluster_pitches = atof(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-min-points") && i+1 < argc) P.shortcut_min_points = (size_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-min-planes") && i+1 < argc) P.shortcut_min_planes = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--shortcut-track-mutual")) P.shortcut_track_mutual = 1;
        else if (!strcmp(argv[i], "--boundary-box") && i+6 < argc) {
            for (int k = 0; k < 6; k++) P.boundary_box[k] = atof(argv[++i]);
            P.boundary_box_set = 1;
        }
        else if (!strcmp(argv[i], "--boundary-box-mesh")) P.boundary_box_from_mesh = 1;
        else if (!strcmp(argv[i], "--fail-on-shortcut")) P.fail_on_shortcut = 1;
        else if (!strcmp(argv[i], "--repair-seams") && i+1 < argc) P.repair_path = argv[++i];
        else if (!strcmp(argv[i], "--repair-dry-run")) P.repair_dry_run = 1;
        else if (!strcmp(argv[i], "--repair-emit-rejected")) P.repair_emit_rejected = 1;
        else if (!strcmp(argv[i], "--repair-cover-all-supported")) P.repair_cover_all_supported = 1;
        else if (!strcmp(argv[i], "--repair-monotone-step")) P.repair_monotone_step = 1;
        else if (!strcmp(argv[i], "--repair-min-chain-support") && i+1 < argc) P.repair_min_chain_support = (size_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--repair-chain-csv") && i+1 < argc) P.repair_chain_csv = argv[++i];
        else if (!strcmp(argv[i], "--repair-radius-pitches") && i+1 < argc) P.repair_radius_pitches = atof(argv[++i]);
        else if (!strcmp(argv[i], "--repair-support-radius-pitches") && i+1 < argc) P.repair_support_radius_pitches = atof(argv[++i]);
        else if (!strcmp(argv[i], "--repair-band-margin") && i+1 < argc) P.repair_band_margin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--repair-phase-offset") && i+1 < argc) P.repair_phase_offset = atof(argv[++i]);
        else if (!strcmp(argv[i], "--repair-auto-offset")) P.repair_auto_offset = 1;
        else if (!strcmp(argv[i], "--repair-max-auto-tracks") && i+1 < argc) {
            if (!parse_size_arg(argv[++i], &P.repair_max_auto_tracks)) {
                fprintf(stderr,
                        "wind_audit: invalid --repair-max-auto-tracks value\n");
                return 2;
            }
        }
        else if (!strcmp(argv[i], "--repair-max-cut-fraction") && i+1 < argc) P.repair_max_cut_fraction = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seam-pitch") && i+1 < argc) P.seam_pitch = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seam-tol") && i+1 < argc) P.seam_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min-faces") && i+1 < argc) P.min_faces = (size_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--dump-mergers") && i+1 < argc) P.dump_mergers = argv[++i];
        else if (!strcmp(argv[i], "--dump-splits") && i+1 < argc) P.dump_splits = argv[++i];
        else if (!strcmp(argv[i], "--dump-shortcuts") && i+1 < argc) P.dump_shortcuts = argv[++i];
        else if (!strcmp(argv[i], "--json") && i+1 < argc) P.json_path = argv[++i];
        else { fprintf(stderr, "wind_audit: unknown arg %s\n", argv[i]); return 2; }
    }
    if (P.pitch <= 0 || P.split_gap <= 0 || P.shortcut_z_step <= 0 ||
        P.shortcut_length_pitches <= 0 || P.shortcut_cluster_pitches <= 0 ||
        P.shortcut_min_points < 2 || P.shortcut_min_planes < 1 ||
        P.repair_radius_pitches <= 0 ||
        P.repair_support_radius_pitches <= 0 ||
        P.repair_band_margin < 0 || fabs(P.repair_phase_offset) >= 0.5 ||
        P.repair_max_cut_fraction <= 0 ||
        P.repair_max_cut_fraction > 1 || P.repair_min_chain_support < 1) {
        fprintf(stderr, "wind_audit: invalid positive distance/count option\n");
        return 2;
    }
    if (P.boundary_box_set &&
        (!isfinite(P.boundary_box[0]) || !isfinite(P.boundary_box[1]) ||
         !isfinite(P.boundary_box[2]) || !isfinite(P.boundary_box[3]) ||
         !isfinite(P.boundary_box[4]) || !isfinite(P.boundary_box[5]) ||
         P.boundary_box[0] >= P.boundary_box[1] ||
         P.boundary_box[2] >= P.boundary_box[3] ||
         P.boundary_box[4] >= P.boundary_box[5])) {
        fprintf(stderr, "wind_audit: invalid --boundary-box extents\n");
        return 2;
    }
    if (P.boundary_box_set && P.boundary_box_from_mesh) {
        fprintf(stderr, "wind_audit: choose --boundary-box or --boundary-box-mesh\n");
        return 2;
    }
    if ((P.repair_path != NULL || P.repair_dry_run) &&
        P.shortcut_tol <= 0.0) {
        fprintf(stderr, "wind_audit: shortcut repair requires --shortcut-tol\n");
        return 2;
    }
    if (P.repair_emit_rejected && P.repair_path == NULL) {
        fprintf(stderr,
                "wind_audit: --repair-emit-rejected requires --repair-seams FILE\n");
        return 2;
    }
    if (P.repair_monotone_step &&
        P.repair_path == NULL && !P.repair_dry_run) {
        fprintf(stderr,
                "wind_audit: --repair-monotone-step requires repair output or dry-run\n");
        return 2;
    }
    if (P.axis_table_path != NULL) {
        if (axis_table_load(P.axis_table_path, &axis) != 0) {
            fprintf(stderr,
                    "wind_audit: cannot load axis table %s (need >=2 unique finite z,y,x rows)\n",
                    P.axis_table_path);
            return 2;
        }
        P.axis = &axis;
    }
    if (P.pitch_table_path != NULL) {
        if (PitchTable_load(P.pitch_table_path, &P.ptab) != 0) {
            fprintf(stderr, "wind_audit: cannot load pitch table %s\n",
                    P.pitch_table_path);
            return 2;
        }
        printf("wind_audit: radial pitch table %s (%d knots); scalar %.2f is "
               "the out-of-range fallback only\n",
               P.pitch_table_path, PitchTable_knots(P.ptab), P.pitch);
    }
    {
        int rc = run(in, &P);
        axis_table_free(&axis);
        PitchTable_free(P.ptab);
        return rc;
    }
}
