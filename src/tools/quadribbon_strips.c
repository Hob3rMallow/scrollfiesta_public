/* quadribbon_strips.c -- expose non-injective quadribbon branches as layers.
 *
 * StrokeStrip supplies one continuous metric coordinate (u,v).  A reconstructed
 * quadribbon can nevertheless contain several physically distinct charts at
 * the same (u,v).  Moving those charts in u would violate the StrokeStrip
 * contract; averaging them in a bake produces the familiar crossed texture.
 * This tool therefore assigns WHOLE charts to at most three discrete sheets.
 * It never changes a coordinate and it assigns every input face exactly once.
 *
 * Strict triangle interiors are rasterized with the same half-texel convention
 * as rawtex_bake.  Shared edges are not conflicts.  The primary sheet is a
 * deterministic, support-ordered independent set in the chart-conflict graph;
 * the remaining charts are split between two auxiliary sheets by a weighted
 * max-cut heuristic.
 *
 *   quadribbon_strips <in.vmesh|obj> <out_prefix>
 *       [--labels reconstruction_component.i32]
 *       [--du F=1] [--dv F=1] [--max-strips N=3] [--write-obj]
 *   quadribbon_strips --selftest
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/mesh_bin.h"
#include "../common/obj_io.h"
#include "../common/ves_platform.h"

enum { QR_MAX_SHEETS = 3, QR_PIXEL_SLOTS = 4 };

/* A conflict edge separates two charts into different sheets only when the
 * collision covers a material fraction of the smaller chart.  A genuinely
 * stacked ply or displaced chart collides over most of its area (ratio near
 * 1); an exact-metric chart expanded against a compressed carried frame only
 * kisses its neighbours at interval boundaries (ratio well under 1e-2, e.g.
 * 19 px on an 85k-vertex chart).  Exiling a whole chart over a boundary kiss
 * manufactures a missing strip in the primary sheet.  Soft edges are still
 * measured, reported per edge in the manifest, and counted per sheet -- they
 * are recorded ambiguity, never hidden. */
static const double QR_HARD_CONFLICT_RATIO = 0.05;

typedef struct QrChartOrder {
    size_t chart;
    uint64_t pixels;
    uint64_t faces;
    uint64_t degree;
} QrChartOrder;

typedef struct QrGraph {
    size_t charts;
    uint64_t *edge;         /* upper triangular entries in a dense matrix */
    uint64_t *chart_pixels; /* strict-interior pixels touched by each chart */
    size_t W, H;
    double umin, vmin, du, dv;
    uint64_t occupied_pixels;
    uint64_t overflow_hits;
    int32_t *slots;         /* W*H*QR_PIXEL_SLOTS unique chart ids */
} QrGraph;

static void *qr_malloc(size_t n)
{
    void *p = malloc(n != 0 ? n : 1);
    if (p == NULL) fprintf(stderr, "quadribbon_strips: out of memory (%zu bytes)\n", n);
    return p;
}

static void *qr_calloc(size_t n, size_t s)
{
    void *p;
    if (s != 0 && n > SIZE_MAX / s) return NULL;
    p = calloc(n != 0 ? n : 1, s);
    if (p == NULL)
        fprintf(stderr, "quadribbon_strips: out of memory (%zu x %zu bytes)\n", n, s);
    return p;
}

static int qr_mul_size(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a) return -1;
    *out = a * b;
    return 0;
}

static int32_t qr_find(int32_t *parent, int32_t v)
{
    int32_t root = v;
    while (parent[root] >= 0) root = parent[root];
    while (v != root) {
        int32_t next = parent[v];
        parent[v] = root;
        v = next;
    }
    return root;
}

static void qr_union(int32_t *parent, int32_t a, int32_t b)
{
    a = qr_find(parent, a);
    b = qr_find(parent, b);
    if (a == b) return;
    if (parent[a] > parent[b]) {
        int32_t t = a;
        a = b;
        b = t;
    }
    parent[a] += parent[b];
    parent[b] = a;
}

static int qr_geometry_components(const MeshBinData *mesh, int32_t **out_parent,
                                  int32_t **out_root_component,
                                  size_t *out_components)
{
    int32_t *parent = NULL, *root_component = NULL;
    size_t i, f, components = 0;

    if (mesh->nv > (size_t)INT32_MAX) {
        fprintf(stderr, "quadribbon_strips: too many vertices for VMESH indices\n");
        return -1;
    }
    parent = (int32_t *)qr_malloc(mesh->nv * sizeof(*parent));
    root_component = (int32_t *)qr_malloc(mesh->nv * sizeof(*root_component));
    if (parent == NULL || root_component == NULL) goto fail;
    for (i = 0; i < mesh->nv; i++) {
        parent[i] = INT32_MIN;          /* not referenced by a face */
        root_component[i] = -1;
    }
    for (f = 0; f < mesh->nf; f++) {
        int32_t a = mesh->faces[f * 3], b = mesh->faces[f * 3 + 1];
        int32_t c = mesh->faces[f * 3 + 2];
        if (parent[a] == INT32_MIN) parent[a] = -1;
        if (parent[b] == INT32_MIN) parent[b] = -1;
        if (parent[c] == INT32_MIN) parent[c] = -1;
        qr_union(parent, a, b);
        qr_union(parent, a, c);
    }
    for (i = 0; i < mesh->nv; i++) {
        if (parent[i] < 0 && parent[i] != INT32_MIN) {
            if (components >= (size_t)INT32_MAX) goto fail;
            root_component[i] = (int32_t)components++;
        }
    }
    *out_parent = parent;
    *out_root_component = root_component;
    *out_components = components;
    return 0;

fail:
    free(parent);
    free(root_component);
    return -1;
}

static int qr_read_labels(const char *path, size_t count, int32_t **out)
{
    FILE *file = fopen(path, "rb");
    int32_t *labels;
    int tail;
    if (file == NULL) {
        fprintf(stderr, "quadribbon_strips: cannot open labels %s: %s\n",
                path, strerror(errno));
        return -1;
    }
    labels = (int32_t *)qr_malloc(count * sizeof(*labels));
    if (labels == NULL) {
        fclose(file);
        return -1;
    }
    if (fread(labels, sizeof(*labels), count, file) != count) {
        fprintf(stderr, "quadribbon_strips: labels %s do not contain %zu int32 values\n",
                path, count);
        free(labels);
        fclose(file);
        return -1;
    }
    tail = fgetc(file);
    if (tail != EOF) {
        fprintf(stderr, "quadribbon_strips: labels %s have trailing data\n", path);
        free(labels);
        fclose(file);
        return -1;
    }
    fclose(file);
    *out = labels;
    return 0;
}

/* A chart is always one exact connected face component.  A reconstruction
 * label is provenance attached to that component, never its identity: the
 * same label can occur on several disconnected pieces, and merging those
 * pieces would hide precisely the UV collisions this tool must separate.
 * Mixed-label faces (or labels changing within a connected component) are a
 * provenance bug; silently choosing one endpoint would be a face-emission
 * rule. */
static int qr_build_face_charts(const MeshBinData *mesh, const char *label_path,
                                int32_t **out_face_chart,
                                int32_t **out_source_label,
                                uint64_t **out_face_count,
                                size_t *out_charts, size_t *out_components)
{
    int32_t *parent = NULL, *root_component = NULL, *labels = NULL;
    int32_t *face_chart = NULL, *source_label = NULL;
    uint64_t *face_count = NULL;
    size_t components = 0, charts = 0, f;

    if (qr_geometry_components(mesh, &parent, &root_component, &components) != 0)
        return -1;
    face_chart = (int32_t *)qr_malloc(mesh->nf * sizeof(*face_chart));
    if (face_chart == NULL) goto fail;
    charts = components;
    source_label = (int32_t *)qr_malloc(
        (charts != 0 ? charts : 1) * sizeof(*source_label));
    if (source_label == NULL) goto fail;
    for (size_t i = 0; i < charts; i++) source_label[i] = INT32_MIN;

    if (label_path != NULL) {
        if (qr_read_labels(label_path, mesh->nv, &labels) != 0) goto fail;
        for (f = 0; f < mesh->nf; f++) {
            int32_t a = mesh->faces[f * 3], b = mesh->faces[f * 3 + 1];
            int32_t c = mesh->faces[f * 3 + 2];
            int32_t label = labels[a];
            int32_t geom = root_component[qr_find(parent, a)];
            if (labels[b] != label || labels[c] != label) {
                fprintf(stderr,
                        "quadribbon_strips: face %zu crosses reconstruction labels %d/%d/%d\n",
                        f, label, labels[b], labels[c]);
                goto fail;
            }
            if (source_label[geom] == INT32_MIN) {
                source_label[geom] = label;
            } else if (source_label[geom] != label) {
                fprintf(stderr,
                        "quadribbon_strips: geometry component %d changes reconstruction label %d -> %d\n",
                        geom, source_label[geom], label);
                goto fail;
            }
            face_chart[f] = geom;
        }
    } else {
        for (size_t i = 0; i < charts; i++) source_label[i] = (int32_t)i;
        for (f = 0; f < mesh->nf; f++) {
            int32_t a = mesh->faces[f * 3];
            face_chart[f] = root_component[qr_find(parent, a)];
        }
    }

    face_count = (uint64_t *)qr_calloc(charts, sizeof(*face_count));
    if (face_count == NULL) goto fail;
    for (f = 0; f < mesh->nf; f++) face_count[face_chart[f]]++;
    free(parent);
    free(root_component);
    free(labels);
    *out_face_chart = face_chart;
    *out_source_label = source_label;
    *out_face_count = face_count;
    *out_charts = charts;
    *out_components = components;
    return 0;

fail:
    free(parent);
    free(root_component);
    free(labels);
    free(face_chart);
    free(source_label);
    free(face_count);
    return -1;
}

static int qr_graph_build(const MeshBinData *mesh, const int32_t *face_chart,
                          size_t charts, double du, double dv, QrGraph *graph)
{
    size_t i, f, pixels, slot_count, matrix_count;
    double umax, vmax;

    memset(graph, 0, sizeof(*graph));
    if (mesh->nv == 0 || mesh->nf == 0 || mesh->uv == NULL ||
        !(du >= 1e-6) || !(dv >= 1e-6) || charts == 0)
        return -1;
    graph->umin = umax = mesh->uv[0];
    graph->vmin = vmax = mesh->uv[1];
    for (i = 0; i < mesh->nv; i++) {
        double u = mesh->uv[i * 2], v = mesh->uv[i * 2 + 1];
        double z = mesh->verts[i * 3], y = mesh->verts[i * 3 + 1];
        double x = mesh->verts[i * 3 + 2];
        if (!isfinite(u) || !isfinite(v) || !isfinite(z) ||
            !isfinite(y) || !isfinite(x)) {
            fprintf(stderr, "quadribbon_strips: non-finite geometry at vertex %zu\n", i);
            return -1;
        }
        if (u < graph->umin) graph->umin = u;
        if (u > umax) umax = u;
        if (v < graph->vmin) graph->vmin = v;
        if (v > vmax) vmax = v;
    }
    graph->W = (size_t)ceil((umax - graph->umin) / du);
    graph->H = (size_t)ceil((vmax - graph->vmin) / dv);
    if (graph->W == 0) graph->W = 1;
    if (graph->H == 0) graph->H = 1;
    graph->du = du;
    graph->dv = dv;
    if (qr_mul_size(graph->W, graph->H, &pixels) != 0 ||
        pixels > ((size_t)1 << 28)) {
        fprintf(stderr,
                "quadribbon_strips: collision raster %zux%zu exceeds 2^28 pixels\n",
                graph->W, graph->H);
        return -1;
    }
    if (qr_mul_size(pixels, QR_PIXEL_SLOTS, &slot_count) != 0 ||
        qr_mul_size(charts, charts, &matrix_count) != 0)
        return -1;
    graph->slots = (int32_t *)qr_malloc(slot_count * sizeof(*graph->slots));
    graph->edge = (uint64_t *)qr_calloc(matrix_count, sizeof(*graph->edge));
    graph->chart_pixels = (uint64_t *)qr_calloc(charts, sizeof(*graph->chart_pixels));
    if (graph->slots == NULL || graph->edge == NULL || graph->chart_pixels == NULL)
        return -1;
    for (i = 0; i < slot_count; i++) graph->slots[i] = -1;
    graph->charts = charts;

    for (f = 0; f < mesh->nf; f++) {
        size_t a = (size_t)mesh->faces[f * 3], b = (size_t)mesh->faces[f * 3 + 1];
        size_t c = (size_t)mesh->faces[f * 3 + 2];
        int32_t chart = face_chart[f];
        double ua = ((double)mesh->uv[a * 2] - graph->umin) / du;
        double va = ((double)mesh->uv[a * 2 + 1] - graph->vmin) / dv;
        double ub = ((double)mesh->uv[b * 2] - graph->umin) / du;
        double vb = ((double)mesh->uv[b * 2 + 1] - graph->vmin) / dv;
        double uc = ((double)mesh->uv[c * 2] - graph->umin) / du;
        double vc = ((double)mesh->uv[c * 2 + 1] - graph->vmin) / dv;
        double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
        double lox, hix, loy, hiy;
        long x0, x1, y0, y1, xx, yy;
        if (fabs(A2) < 1e-12) continue;
        lox = fmin(ua, fmin(ub, uc)); hix = fmax(ua, fmax(ub, uc));
        loy = fmin(va, fmin(vb, vc)); hiy = fmax(va, fmax(vb, vc));
        /* Centers are sampled at xx+0.5: the first covered center is
         * ceil(lo - 0.5).  ceil(lo - 0.001) assumed integer-lattice UVs and
         * dropped each face's first row/column at generic UV phases. */
        x0 = (long)ceil(lox - 0.501); x1 = (long)floor(hix + 0.001);
        y0 = (long)ceil(loy - 0.501); y1 = (long)floor(hiy + 0.001);
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 >= (long)graph->W) x1 = (long)graph->W - 1;
        if (y1 >= (long)graph->H) y1 = (long)graph->H - 1;
        for (yy = y0; yy <= y1; yy++) {
            for (xx = x0; xx <= x1; xx++) {
                double pu = (double)xx + 0.5, pv = (double)yy + 0.5;
                double l1 = ((pu - ua) * (vc - va) -
                             (pv - va) * (uc - ua)) / A2;
                double l2 = ((ub - ua) * (pv - va) -
                             (vb - va) * (pu - ua)) / A2;
                double l0 = 1.0 - l1 - l2;
                size_t pi, base, s, empty = QR_PIXEL_SLOTS;
                int duplicate = 0;
                if (l0 <= 1e-9 || l1 <= 1e-9 || l2 <= 1e-9) continue;
                pi = (size_t)yy * graph->W + (size_t)xx;
                base = pi * QR_PIXEL_SLOTS;
                for (s = 0; s < QR_PIXEL_SLOTS; s++) {
                    int32_t other = graph->slots[base + s];
                    if (other == chart) { duplicate = 1; break; }
                    if (other < 0 && empty == QR_PIXEL_SLOTS) empty = s;
                }
                if (duplicate) continue;
                for (s = 0; s < QR_PIXEL_SLOTS; s++) {
                    int32_t other = graph->slots[base + s];
                    size_t lo, hi;
                    if (other < 0) continue;
                    lo = (size_t)(other < chart ? other : chart);
                    hi = (size_t)(other < chart ? chart : other);
                    graph->edge[lo * charts + hi]++;
                }
                if (empty < QR_PIXEL_SLOTS) {
                    if (graph->slots[base] < 0) graph->occupied_pixels++;
                    graph->slots[base + empty] = chart;
                    graph->chart_pixels[chart]++;
                } else {
                    graph->overflow_hits++;
                }
            }
        }
    }
    return 0;
}

static void qr_graph_dispose(QrGraph *graph)
{
    free(graph->edge);
    free(graph->chart_pixels);
    free(graph->slots);
    memset(graph, 0, sizeof(*graph));
}

static const QrChartOrder *qr_sort_order = NULL;

static int qr_order_cmp(const void *aa, const void *bb)
{
    const QrChartOrder *a = (const QrChartOrder *)aa;
    const QrChartOrder *b = (const QrChartOrder *)bb;
    (void)qr_sort_order;
    if (a->pixels != b->pixels) return a->pixels > b->pixels ? -1 : 1;
    if (a->faces != b->faces) return a->faces > b->faces ? -1 : 1;
    return a->chart < b->chart ? -1 : a->chart > b->chart;
}

static int qr_aux_cmp(const void *aa, const void *bb)
{
    const QrChartOrder *a = (const QrChartOrder *)aa;
    const QrChartOrder *b = (const QrChartOrder *)bb;
    if (a->degree != b->degree) return a->degree > b->degree ? -1 : 1;
    return qr_order_cmp(aa, bb);
}

static uint64_t qr_edge(const QrGraph *g, size_t a, size_t b)
{
    size_t lo, hi;
    if (a == b) return 0;
    lo = a < b ? a : b;
    hi = a < b ? b : a;
    return g->edge[lo * g->charts + hi];
}

/* Hard = sheet-separating.  Soft = boundary kiss, kept in one sheet. */
static int qr_edge_hard(const QrGraph *g, size_t a, size_t b)
{
    uint64_t px = qr_edge(g, a, b);
    uint64_t m;
    if (px == 0) return 0;
    m = g->chart_pixels[a] < g->chart_pixels[b]
      ? g->chart_pixels[a] : g->chart_pixels[b];
    if (m == 0) return 1;
    return (double)px > QR_HARD_CONFLICT_RATIO * (double)m;
}

static int qr_assign_sheets(const QrGraph *g, const uint64_t *face_count,
                            int max_sheets, int8_t *sheet,
                            uint64_t sheet_pixels[QR_MAX_SHEETS],
                            uint64_t intra_hard[QR_MAX_SHEETS],
                            uint64_t intra_soft[QR_MAX_SHEETS])
{
    QrChartOrder *order = NULL, *aux = NULL;
    size_t n = g->charts, naux = 0, i, j;
    int used_sheets = 1;

    order = (QrChartOrder *)qr_malloc(n * sizeof(*order));
    aux = (QrChartOrder *)qr_malloc(n * sizeof(*aux));
    if (order == NULL || aux == NULL) goto fail;
    memset(sheet, -1, n * sizeof(*sheet));
    memset(sheet_pixels, 0, QR_MAX_SHEETS * sizeof(*sheet_pixels));
    memset(intra_hard, 0, QR_MAX_SHEETS * sizeof(*intra_hard));
    memset(intra_soft, 0, QR_MAX_SHEETS * sizeof(*intra_soft));
    for (i = 0; i < n; i++) {
        order[i].chart = i;
        order[i].pixels = g->chart_pixels[i];
        order[i].faces = face_count[i];
        order[i].degree = 0;
    }
    qsort(order, n, sizeof(*order), qr_order_cmp);

    /* Hard-collision-free primary independent set, biased toward material
     * pixels.  Soft boundary kisses never displace a chart. */
    for (i = 0; i < n; i++) {
        size_t c = order[i].chart;
        int compatible = 1;
        for (j = 0; j < n; j++) {
            if (sheet[j] == 0 && qr_edge_hard(g, c, j)) {
                compatible = 0;
                break;
            }
        }
        if (compatible) sheet[c] = 0;
    }
    if (max_sheets == 1) {
        for (i = 0; i < n; i++) sheet[i] = 0;
        goto count;
    }

    for (i = 0; i < n; i++) {
        size_t c = order[i].chart;
        uint64_t degree = 0;
        if (sheet[c] == 0) continue;
        for (j = 0; j < n; j++)
            if (sheet[j] != 0 && qr_edge_hard(g, c, j))
                degree += qr_edge(g, c, j);
        aux[naux] = order[i];
        aux[naux].degree = degree;
        naux++;
    }
    qsort(aux, naux, sizeof(*aux), qr_aux_cmp);
    if (naux != 0) used_sheets = 2;
    for (i = 0; i < naux; i++) {
        size_t c = aux[i].chart;
        uint64_t cost1 = 0, cost2 = 0, load1 = 0, load2 = 0;
        if (max_sheets == 2) {
            sheet[c] = 1;
            continue;
        }
        for (j = 0; j < n; j++) {
            if (!qr_edge_hard(g, c, j)) continue;
            if (sheet[j] == 1) { cost1 += qr_edge(g, c, j); load1 += g->chart_pixels[j]; }
            if (sheet[j] == 2) { cost2 += qr_edge(g, c, j); load2 += g->chart_pixels[j]; }
        }
        if (cost1 == 0 && cost2 == 0) {
            for (j = 0; j < n; j++) {
                if (sheet[j] == 1) load1 += g->chart_pixels[j];
                if (sheet[j] == 2) load2 += g->chart_pixels[j];
            }
        }
        sheet[c] = cost1 < cost2 ? 1 : cost2 < cost1 ? 2 : load1 <= load2 ? 1 : 2;
        if (sheet[c] == 2) used_sheets = 3;
    }

    /* Strictly improving flips minimize residual hard auxiliary collisions. */
    if (max_sheets >= 3) {
        int changed = 1;
        for (int pass = 0; pass < 32 && changed; pass++) {
            changed = 0;
            for (i = 0; i < naux; i++) {
                size_t c = aux[i].chart;
                int8_t current = sheet[c], other = current == 1 ? 2 : 1;
                uint64_t current_cost = 0, other_cost = 0;
                for (j = 0; j < n; j++) {
                    if (!qr_edge_hard(g, c, j)) continue;
                    if (sheet[j] == current) current_cost += qr_edge(g, c, j);
                    if (sheet[j] == other) other_cost += qr_edge(g, c, j);
                }
                if (other_cost < current_cost) {
                    sheet[c] = other;
                    changed = 1;
                    if (other == 2) used_sheets = 3;
                }
            }
        }
    }

count:
    for (i = 0; i < n; i++) sheet_pixels[sheet[i]] += g->chart_pixels[i];
    for (i = 0; i < n; i++)
        for (j = i + 1; j < n; j++) {
            if (sheet[i] != sheet[j]) continue;
            if (qr_edge_hard(g, i, j))
                intra_hard[sheet[i]] += qr_edge(g, i, j);
            else
                intra_soft[sheet[i]] += qr_edge(g, i, j);
        }
    free(order);
    free(aux);
    return used_sheets;

fail:
    free(order);
    free(aux);
    return -1;
}

static void qr_json_string(FILE *file, const char *text)
{
    const unsigned char *p = (const unsigned char *)text;
    fputc('"', file);
    for (; *p; p++) {
        if (*p == '"' || *p == '\\') fprintf(file, "\\%c", *p);
        else if (*p == '\n') fputs("\\n", file);
        else if (*p == '\r') fputs("\\r", file);
        else if (*p == '\t') fputs("\\t", file);
        else if (*p < 32) fprintf(file, "\\u%04x", (unsigned)*p);
        else fputc(*p, file);
    }
    fputc('"', file);
}

static int qr_write_sheet(const MeshBinData *mesh, const int32_t *face_chart,
                          const int8_t *sheet, int which, const char *prefix,
                          int write_obj, size_t *out_nv, size_t *out_nf)
{
    int32_t *map = NULL, *faces = NULL;
    float *verts = NULL, *uv = NULL;
    size_t i, f, nv = 0, nf = 0, fi = 0;
    char vmesh_path[2048], obj_path[2048];

    if (snprintf(vmesh_path, sizeof(vmesh_path), "%s_%02d.vmesh", prefix, which) < 0 ||
        snprintf(obj_path, sizeof(obj_path), "%s_%02d.obj", prefix, which) < 0)
        return -1;
    map = (int32_t *)qr_malloc(mesh->nv * sizeof(*map));
    if (map == NULL) return -1;
    for (i = 0; i < mesh->nv; i++) map[i] = -1;
    for (f = 0; f < mesh->nf; f++) {
        if (sheet[face_chart[f]] != which) continue;
        nf++;
        for (int k = 0; k < 3; k++) map[mesh->faces[f * 3 + (size_t)k]] = 0;
    }
    for (i = 0; i < mesh->nv; i++)
        if (map[i] == 0) map[i] = (int32_t)nv++;
    verts = (float *)qr_malloc(nv * 3 * sizeof(*verts));
    uv = (float *)qr_malloc(nv * 2 * sizeof(*uv));
    faces = (int32_t *)qr_malloc(nf * 3 * sizeof(*faces));
    if (verts == NULL || uv == NULL || faces == NULL) goto fail;
    for (i = 0; i < mesh->nv; i++) {
        int32_t j = map[i];
        if (j < 0) continue;
        memcpy(&verts[(size_t)j * 3], &mesh->verts[i * 3], 3 * sizeof(*verts));
        memcpy(&uv[(size_t)j * 2], &mesh->uv[i * 2], 2 * sizeof(*uv));
    }
    for (f = 0; f < mesh->nf; f++) {
        if (sheet[face_chart[f]] != which) continue;
        faces[fi * 3] = map[mesh->faces[f * 3]];
        faces[fi * 3 + 1] = map[mesh->faces[f * 3 + 1]];
        faces[fi * 3 + 2] = map[mesh->faces[f * 3 + 2]];
        fi++;
    }
    if (fi != nf || ves_ensure_parent_dir(vmesh_path) != 0 ||
        MeshBin_write(vmesh_path, verts, nv, faces, nf, uv) != 0) {
        fprintf(stderr, "quadribbon_strips: cannot write %s\n", vmesh_path);
        goto fail;
    }
    if (write_obj && (ves_ensure_parent_dir(obj_path) != 0 ||
                      ObjIO_write_uv(obj_path, verts, nv, faces, nf, uv) != 0)) {
        fprintf(stderr, "quadribbon_strips: cannot write %s\n", obj_path);
        goto fail;
    }
    fprintf(stderr, "  sheet %d: %zu vertices, %zu faces -> %s%s\n",
            which, nv, nf, vmesh_path, write_obj ? " + OBJ" : "");
    *out_nv = nv;
    *out_nf = nf;
    free(map); free(verts); free(uv); free(faces);
    return 0;

fail:
    free(map); free(verts); free(uv); free(faces);
    return -1;
}

static int qr_write_manifest(const char *prefix, const char *input,
                             const char *labels, const QrGraph *graph,
                             const int32_t *source_label,
                             const uint64_t *face_count, const int8_t *sheet,
                             int sheets, const size_t *sheet_nv,
                             const size_t *sheet_nf,
                             const uint64_t sheet_pixels[QR_MAX_SHEETS],
                             const uint64_t intra_hard[QR_MAX_SHEETS],
                             const uint64_t intra_soft[QR_MAX_SHEETS],
                             size_t input_nv, size_t input_nf,
                             size_t geometry_components)
{
    char path[2048], mesh_path[2048];
    FILE *file;
    size_t conflict_edges = 0, conflict_edges_hard = 0;
    if (snprintf(path, sizeof(path), "%s_manifest.json", prefix) < 0 ||
        ves_ensure_parent_dir(path) != 0)
        return -1;
    file = fopen(path, "wb");
    if (file == NULL) return -1;
    for (size_t a = 0; a < graph->charts; a++)
        for (size_t b = a + 1; b < graph->charts; b++)
            if (qr_edge(graph, a, b) != 0) {
                conflict_edges++;
                if (qr_edge_hard(graph, a, b)) conflict_edges_hard++;
            }
    fputs("{\n  \"contract\": \"strokestrip-discrete-branch-sheets-v1\",\n", file);
    fputs("  \"input\": ", file); qr_json_string(file, input); fputs(",\n", file);
    fputs("  \"labels\": ", file);
    if (labels != NULL) qr_json_string(file, labels); else fputs("null", file);
    fprintf(file,
            ",\n  \"coordinate_rule\": \"XYZ and UV copied bit-for-bit; no cuts, packing, or deformation\",\n"
            "  \"input_vertices\": %zu,\n  \"input_faces\": %zu,\n"
            "  \"geometry_components\": %zu,\n  \"charts\": %zu,\n"
            "  \"raster\": {\"width\": %zu, \"height\": %zu, "
            "\"du\": %.9g, \"dv\": %.9g, \"occupied_strict_pixels\": %" PRIu64
            ", \"slot_overflow_hits\": %" PRIu64 "},\n"
            "  \"conflict_edges\": %zu,\n"
            "  \"conflict_edges_hard\": %zu,\n"
            "  \"hard_conflict_ratio\": %.9g,\n"
            "  \"conflict_edge_list\": [",
            input_nv, input_nf, geometry_components, graph->charts,
            graph->W, graph->H, graph->du, graph->dv,
            graph->occupied_pixels, graph->overflow_hits, conflict_edges,
            conflict_edges_hard, QR_HARD_CONFLICT_RATIO);
    {
        int first_edge = 1;
        for (size_t a = 0; a < graph->charts; a++)
            for (size_t b = a + 1; b < graph->charts; b++) {
                uint64_t px = qr_edge(graph, a, b);
                uint64_t m;
                if (px == 0) continue;
                m = graph->chart_pixels[a] < graph->chart_pixels[b]
                  ? graph->chart_pixels[a] : graph->chart_pixels[b];
                fprintf(file,
                        "%s\n    {\"chart_a\": %zu, \"chart_b\": %zu, "
                        "\"label_a\": %d, \"label_b\": %d, "
                        "\"pixels\": %" PRIu64 ", \"smaller_chart_pixels\": %"
                        PRIu64 ", \"hard\": %s, "
                        "\"sheet_a\": %d, \"sheet_b\": %d}",
                        first_edge ? "" : ",", a, b,
                        source_label[a], source_label[b], px, m,
                        qr_edge_hard(graph, a, b) ? "true" : "false",
                        (int)sheet[a], (int)sheet[b]);
                first_edge = 0;
            }
        fputs(first_edge ? "],\n" : "\n  ],\n", file);
    }
    fputs("  \"sheets\": [\n", file);
    for (int s = 0; s < sheets; s++) {
        snprintf(mesh_path, sizeof(mesh_path), "%s_%02d.vmesh", prefix, s);
        fprintf(file, "    {\"index\": %d, \"mesh\": ", s);
        qr_json_string(file, mesh_path);
        fprintf(file,
                ", \"vertices\": %zu, \"faces\": %zu, "
                "\"strict_pixel_support_sum\": %" PRIu64
                ", \"intra_sheet_hard_conflict_support\": %" PRIu64
                ", \"intra_sheet_soft_conflict_support\": %" PRIu64
                ", \"source_labels\": [",
                sheet_nv[s], sheet_nf[s], sheet_pixels[s],
                intra_hard[s], intra_soft[s]);
        {
            int first = 1;
            for (size_t c = 0; c < graph->charts; c++) {
                if (sheet[c] != s) continue;
                if (!first) fputs(", ", file);
                fprintf(file, "%d", source_label[c]);
                first = 0;
            }
        }
        fputs("]}", file);
        fputs(s + 1 == sheets ? "\n" : ",\n", file);
    }
    {
        size_t total_faces = 0;
        for (int s = 0; s < sheets; s++) total_faces += sheet_nf[s];
        fprintf(file,
                "  ],\n  \"face_conservation\": {\"output_faces\": %zu, "
                "\"exact\": %s}\n}\n",
                total_faces, total_faces == input_nf ? "true" : "false");
    }
    if (fclose(file) != 0) return -1;
    fprintf(stderr, "  manifest -> %s\n", path);
    (void)face_count;
    return 0;
}

static int qr_run(const char *input, const char *prefix, const char *labels,
                  double du, double dv, int max_sheets, int write_obj)
{
    MeshBinData mesh;
    char resolved[2048];
    int32_t *face_chart = NULL, *source_label = NULL;
    uint64_t *face_count = NULL;
    int8_t *sheet = NULL;
    size_t charts = 0, components = 0, total_faces = 0;
    size_t sheet_nv[QR_MAX_SHEETS] = {0}, sheet_nf[QR_MAX_SHEETS] = {0};
    uint64_t sheet_pixels[QR_MAX_SHEETS];
    uint64_t intra_hard[QR_MAX_SHEETS], intra_soft[QR_MAX_SHEETS];
    QrGraph graph;
    int sheets, rc = 1;
    double t0 = ves_clock_sec();

    memset(&mesh, 0, sizeof(mesh));
    memset(&graph, 0, sizeof(graph));
    if (MeshBin_companion_path(input, resolved, sizeof(resolved)) != 0 ||
        MeshBin_read_malloc(resolved, &mesh) != 0 || mesh.uv == NULL) {
        fprintf(stderr, "quadribbon_strips: cannot read UV VMESH %s\n", input);
        return 1;
    }
    fprintf(stderr, "quadribbon_strips: %s: %zu vertices, %zu faces\n",
            resolved, mesh.nv, mesh.nf);
    if (qr_build_face_charts(&mesh, labels, &face_chart, &source_label,
                             &face_count, &charts, &components) != 0)
        goto done;
    fprintf(stderr, "  identity: %zu geometry components -> %zu whole charts (%s)\n",
            components, charts, labels != NULL ? "reconstruction sidecar" : "geometry");
    if (qr_graph_build(&mesh, face_chart, charts, du, dv, &graph) != 0)
        goto done;
    fprintf(stderr,
            "  strict collision raster: %zux%zu, %" PRIu64
            " occupied pixels, %" PRIu64 " >%d-layer hits\n",
            graph.W, graph.H, graph.occupied_pixels,
            graph.overflow_hits, QR_PIXEL_SLOTS);
    sheet = (int8_t *)qr_malloc(charts * sizeof(*sheet));
    if (sheet == NULL) goto done;
    sheets = qr_assign_sheets(&graph, face_count, max_sheets, sheet,
                              sheet_pixels, intra_hard, intra_soft);
    if (sheets < 1) goto done;
    for (int s = 0; s < sheets; s++) {
        if (qr_write_sheet(&mesh, face_chart, sheet, s, prefix, write_obj,
                           &sheet_nv[s], &sheet_nf[s]) != 0)
            goto done;
        total_faces += sheet_nf[s];
        fprintf(stderr,
                "    sheet %d graph audit: strict support=%" PRIu64
                " intra-conflict hard=%" PRIu64 " soft-kiss=%" PRIu64 "\n",
                s, sheet_pixels[s], intra_hard[s], intra_soft[s]);
    }
    if (total_faces != mesh.nf) {
        fprintf(stderr,
                "quadribbon_strips: face conservation failed: %zu != %zu\n",
                total_faces, mesh.nf);
        goto done;
    }
    if (qr_write_manifest(prefix, resolved, labels, &graph, source_label,
                          face_count, sheet, sheets, sheet_nv, sheet_nf,
                          sheet_pixels, intra_hard, intra_soft,
                          mesh.nv, mesh.nf, components) != 0)
        goto done;
    fprintf(stderr,
            "quadribbon_strips: OK: %d sheet(s), all %zu faces conserved, "
            "primary hard conflicts=%" PRIu64 " (soft kisses=%" PRIu64
            ") (%.2fs)\n",
            sheets, mesh.nf, intra_hard[0], intra_soft[0],
            ves_clock_sec() - t0);
    rc = intra_hard[0] == 0 ? 0 : 1;

done:
    qr_graph_dispose(&graph);
    free(face_chart);
    free(source_label);
    free(face_count);
    free(sheet);
    MeshBin_dispose(&mesh);
    return rc;
}

static int qr_selftest(void)
{
    QrGraph g;
    uint64_t face_count[5] = {100, 50, 40, 20, 10};
    int8_t sheet[5];
    uint64_t sp[QR_MAX_SHEETS];
    uint64_t ih[QR_MAX_SHEETS], is_[QR_MAX_SHEETS];
    int sheets, fail = 0;
    memset(&g, 0, sizeof(g));
    g.charts = 5;
    g.edge = (uint64_t *)calloc(25, sizeof(*g.edge));
    g.chart_pixels = (uint64_t *)calloc(5, sizeof(*g.chart_pixels));
    if (g.edge == NULL || g.chart_pixels == NULL) return 1;
    for (size_t i = 0; i < 5; i++) g.chart_pixels[i] = face_count[i];
    g.edge[0 * 5 + 1] = 12;
    g.edge[0 * 5 + 2] = 9;
    g.edge[1 * 5 + 2] = 4;
    g.edge[2 * 5 + 3] = 2;
    sheets = qr_assign_sheets(&g, face_count, 3, sheet, sp, ih, is_);
    fail |= sheets < 2 || sheets > 3;
    fail |= sheet[0] != 0 || sheet[1] == 0 || sheet[2] == 0;
    fail |= ih[0] != 0;
    for (size_t i = 0; i < 5; i++) fail |= sheet[i] < 0 || sheet[i] >= 3;
    qr_graph_dispose(&g);

    /* Two large charts joined only by a boundary kiss far below the hard
     * ratio must SHARE the primary sheet; the kiss is counted as soft. */
    {
        QrGraph k;
        uint64_t kiss_faces[2] = {9000, 7000};
        int8_t kiss_sheet[2];
        uint64_t ksp[QR_MAX_SHEETS], kih[QR_MAX_SHEETS], kis[QR_MAX_SHEETS];
        memset(&k, 0, sizeof(k));
        k.charts = 2;
        k.edge = (uint64_t *)calloc(4, sizeof(*k.edge));
        k.chart_pixels = (uint64_t *)calloc(2, sizeof(*k.chart_pixels));
        if (k.edge == NULL || k.chart_pixels == NULL) {
            qr_graph_dispose(&k);
            return 1;
        }
        k.chart_pixels[0] = 10000;
        k.chart_pixels[1] = 8000;
        k.edge[0 * 2 + 1] = 19;   /* the measured 4x5x5 boundary kiss scale */
        sheets = qr_assign_sheets(&k, kiss_faces, 3, kiss_sheet,
                                  ksp, kih, kis);
        fail |= sheets != 1;
        fail |= kiss_sheet[0] != 0 || kiss_sheet[1] != 0;
        fail |= kih[0] != 0 || kis[0] != 19;
        qr_graph_dispose(&k);
    }

    /* Raster half-texel regression: a quad whose UV origin sits at a
     * FRACTIONAL raster phase (u=10.3, v=0.2) must still cover every pixel
     * whose center lies inside it.  The old ceil(lo-0.001) bound skipped the
     * first column of every fractionally-phased face (the metric-projection
     * moire).  Quad A pins the raster origin; its integer-lattice diagonal
     * passes exactly through two pixel centers, which the strict-interior
     * rule excludes (2 of 4 px), the documented lattice pathology.  Quad B's
     * four covered centers avoid its diagonal and must ALL survive. */
    {
        static const float quv[16] = {
            0.0f, 0.0f,  2.0f, 0.0f,  0.0f, 2.0f,  2.0f, 2.0f,
            10.3f, 0.2f, 12.3f, 0.2f, 10.3f, 2.2f, 12.3f, 2.2f
        };
        static const int32_t qf[12] = { 0,1,2, 1,3,2, 4,5,6, 5,7,6 };
        float qv[24];
        int32_t fchart[4] = {0, 0, 0, 0};
        MeshBinData m;
        QrGraph rg;
        memset(qv, 0, sizeof(qv));
        memset(&m, 0, sizeof(m));
        m.verts = qv;
        m.uv = (float *)quv;
        m.faces = (int32_t *)qf;
        m.nv = 8;
        m.nf = 4;
        if (qr_graph_build(&m, fchart, 1, 1.0, 1.0, &rg) != 0) {
            fail |= 1;
        } else {
            fail |= rg.occupied_pixels != 6;   /* A: 2 of 4; B: all 4 */
            fail |= rg.chart_pixels[0] != 6;
            qr_graph_dispose(&rg);
        }
    }
    fprintf(stderr, "quadribbon_strips selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}

int main(int argc, char **argv)
{
    const char *labels = NULL;
    /* Audit at one source voxel.  A two-voxel grid can miss a genuine overlap
     * solely because the raster origin moves by half a texel. */
    double du = 1.0, dv = 1.0;
    int max_sheets = 3, write_obj = 0;
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return qr_selftest() ? 3 : 0;
    if (argc < 3) {
        fprintf(stderr,
                "usage: %s <in.vmesh|obj> <out_prefix> "
                "[--labels file.i32] [--du F] [--dv F] "
                "[--max-strips 1..3] [--write-obj]\n"
                "       %s --selftest\n",
                argv[0], argv[0]);
        return 2;
    }
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--labels") == 0 && i + 1 < argc)
            labels = argv[++i];
        else if (strcmp(argv[i], "--du") == 0 && i + 1 < argc)
            du = atof(argv[++i]);
        else if (strcmp(argv[i], "--dv") == 0 && i + 1 < argc)
            dv = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-strips") == 0 && i + 1 < argc)
            max_sheets = atoi(argv[++i]);
        else if (strcmp(argv[i], "--write-obj") == 0)
            write_obj = 1;
        else {
            fprintf(stderr, "quadribbon_strips: unknown/incomplete option %s\n", argv[i]);
            return 2;
        }
    }
    if (!(du >= 1e-6) || !(dv >= 1e-6) || max_sheets < 1 || max_sheets > 3) {
        fprintf(stderr,
                "quadribbon_strips: du/dv must be >=1e-6 and max-strips in [1,3]\n");
        return 2;
    }
    return qr_run(argv[1], argv[2], labels, du, dv, max_sheets, write_obj);
}
