/* asm_report.c -- PNG previews for the assembly stages. */
#include "asm_report.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/ves_png.h"

AsmCanvas AsmCanvas_new(Arena_T arena, int w, int h, uint8_t r, uint8_t g, uint8_t b)
{
    AsmCanvas c;
    c.w = w; c.h = h;
    c.rgb = ARENA_ALLOC(arena, (size_t)w * (size_t)h * 3);
    for (size_t i = 0; i < (size_t)w * (size_t)h; i++) {
        c.rgb[i*3] = r; c.rgb[i*3+1] = g; c.rgb[i*3+2] = b;
    }
    return c;
}

static double ar_edge(double ax, double ay, double bx, double by, double px, double py)
{
    return (bx - ax) * (py - ay) - (by - ay) * (px - ax);
}

static void ar_put(AsmCanvas *c, int x, int y, const uint8_t rgb[3])
{
    if (x < 0 || y < 0 || x >= c->w || y >= c->h) return;
    size_t k = ((size_t)y * (size_t)c->w + (size_t)x) * 3;
    c->rgb[k] = rgb[0]; c->rgb[k+1] = rgb[1]; c->rgb[k+2] = rgb[2];
}

void AsmCanvas_tri(AsmCanvas *c, double x0, double y0, double x1, double y1,
                   double x2, double y2, const uint8_t rgb[3])
{
    double minx = fmin(x0, fmin(x1, x2)), maxx = fmax(x0, fmax(x1, x2));
    double miny = fmin(y0, fmin(y1, y2)), maxy = fmax(y0, fmax(y1, y2));
    int ix0 = (int)floor(minx), ix1 = (int)ceil(maxx);
    int iy0 = (int)floor(miny), iy1 = (int)ceil(maxy);
    if (ix0 < 0) ix0 = 0;
    if (iy0 < 0) iy0 = 0;
    if (ix1 >= c->w) ix1 = c->w - 1;
    if (iy1 >= c->h) iy1 = c->h - 1;
    double area = ar_edge(x0, y0, x1, y1, x2, y2);
    if (fabs(area) < 1e-12) { ar_put(c, (int)x0, (int)y0, rgb); return; }
    double sgn = area > 0 ? 1.0 : -1.0;
    for (int y = iy0; y <= iy1; y++) {
        double py = y + 0.5;
        for (int x = ix0; x <= ix1; x++) {
            double px = x + 0.5;
            double w0 = sgn * ar_edge(x1, y1, x2, y2, px, py);
            double w1 = sgn * ar_edge(x2, y2, x0, y0, px, py);
            double w2 = sgn * ar_edge(x0, y0, x1, y1, px, py);
            if (w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9) continue;
            ar_put(c, x, y, rgb);
        }
    }
    ar_put(c, (int)((x0 + x1 + x2) / 3.0), (int)((y0 + y1 + y2) / 3.0), rgb);
}

void AsmCanvas_line(AsmCanvas *c, double x0, double y0, double x1, double y1,
                    const uint8_t rgb[3])
{
    double dx = x1 - x0, dy = y1 - y0;
    double len = fmax(fabs(dx), fabs(dy));
    int n = (int)ceil(len) + 1;
    for (int i = 0; i <= n; i++) {
        double t = n > 0 ? (double)i / (double)n : 0.0;
        ar_put(c, (int)(x0 + t * dx), (int)(y0 + t * dy), rgb);
    }
}

void AsmCanvas_rect(AsmCanvas *c, double x0, double y0, double x1, double y1,
                    const uint8_t rgb[3])
{
    int ix0 = (int)floor(fmin(x0, x1)), ix1 = (int)ceil(fmax(x0, x1));
    int iy0 = (int)floor(fmin(y0, y1)), iy1 = (int)ceil(fmax(y0, y1));
    for (int y = iy0; y < iy1; y++) for (int x = ix0; x < ix1; x++) ar_put(c, x, y, rgb);
    if (ix1 <= ix0 || iy1 <= iy0) ar_put(c, ix0, iy0, rgb);
}

int AsmCanvas_write(const AsmCanvas *c, const char *path)
{
    return VesPng_write_rgb(path, c->rgb, c->w, c->h);
}

void AsmReport_colour(uint32_t key, uint8_t rgb[3])
{
    double h = fmod((double)key * 0.61803398875, 1.0) * 6.0;
    double s = 0.55, v = 0.85;
    int i = (int)floor(h);
    double f = h - i;
    double p = v * (1 - s), q = v * (1 - s * f), t = v * (1 - s * (1 - f));
    double r, g, b;
    switch (i % 6) {
    case 0: r = v; g = t; b = p; break;
    case 1: r = q; g = v; b = p; break;
    case 2: r = p; g = v; b = t; break;
    case 3: r = p; g = q; b = v; break;
    case 4: r = t; g = p; b = v; break;
    default: r = v; g = p; b = q; break;
    }
    rgb[0] = (uint8_t)(r * 255.0); rgb[1] = (uint8_t)(g * 255.0); rgb[2] = (uint8_t)(b * 255.0);
}

static void ar_chart_colour(const AsmChart *c, uint8_t rgb[3])
{
    AsmReport_colour((uint32_t)c->cube * 7919u + 13u, rgb);
    double s = c->stress_frac * 5.0;
    if (s > 1.0) s = 1.0;
    if (s > 0.0) {
        rgb[0] = (uint8_t)(rgb[0] + (255 - rgb[0]) * s);
        rgb[1] = (uint8_t)(rgb[1] * (1.0 - 0.7 * s));
        rgb[2] = (uint8_t)(rgb[2] * (1.0 - 0.7 * s));
    }
}

typedef struct ArOrder { double area; int32_t idx; } ArOrder;

static int ar_cmp_area_desc(const void *a, const void *b)
{
    const ArOrder *x = a, *y = b;
    if (x->area > y->area) return -1;
    if (x->area < y->area) return 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

int AsmReport_charts_atlas_png(Arena_T arena, const char *path,
                               const AsmChart *charts, size_t n,
                               int cell_px, size_t max_cells)
{
    Arena_Mark mark = Arena_save(arena);
    ArOrder *ord = ARENA_ALLOC(arena, (n ? n : 1) * sizeof(ArOrder));
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        if (!AsmChart_in_layout(&charts[i])) continue;
        ord[m].area = charts[i].area3d; ord[m].idx = (int32_t)i; m++;
    }
    if (m == 0) { Arena_restore(arena, mark); return -1; }
    qsort(ord, m, sizeof(ArOrder), ar_cmp_area_desc);
    if (m > max_cells) m = max_cells;
    int cols = (int)ceil(sqrt((double)m));
    int rows = (int)((m + (size_t)cols - 1) / (size_t)cols);
    AsmCanvas cv = AsmCanvas_new(arena, cols * cell_px, rows * cell_px, 24, 24, 28);
    for (size_t k = 0; k < m; k++) {
        const AsmChart *c = &charts[(size_t)ord[k].idx];
        double lo0 = 1e300, lo1 = 1e300, hi0 = -1e300, hi1 = -1e300;
        for (size_t i = 0; i < c->nv; i++) {
            double u = c->uv[i*2], v = c->uv[i*2+1];
            if (u < lo0) lo0 = u; if (u > hi0) hi0 = u;
            if (v < lo1) lo1 = v; if (v > hi1) hi1 = v;
        }
        double span = fmax(hi0 - lo0, hi1 - lo1);
        if (span <= 0.0) span = 1.0;
        double scale = (cell_px - 4) / span;
        double ox = (double)((int)k % cols) * cell_px + 2.0 - lo0 * scale;
        double oy = (double)((int)k / cols) * cell_px + 2.0 - lo1 * scale;
        uint8_t rgb[3];
        ar_chart_colour(c, rgb);
        const uint8_t magenta[3] = { 255, 0, 255 };
        for (size_t f = 0; f < c->nf; f++) {
            const int32_t *fv = &c->faces[f*3];
            double x0 = c->uv[(size_t)fv[0]*2] * scale + ox, y0 = c->uv[(size_t)fv[0]*2+1] * scale + oy;
            double x1 = c->uv[(size_t)fv[1]*2] * scale + ox, y1 = c->uv[(size_t)fv[1]*2+1] * scale + oy;
            double x2 = c->uv[(size_t)fv[2]*2] * scale + ox, y2 = c->uv[(size_t)fv[2]*2+1] * scale + oy;
            double det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
            AsmCanvas_tri(&cv, x0, y0, x1, y1, x2, y2, det < 0 ? magenta : rgb);
        }
    }
    int rc = AsmCanvas_write(&cv, path);
    Arena_restore(arena, mark);
    return rc;
}

static void ar_pose_apply(const AsmChart *c, double u, double v, double *gx, double *gy)
{
    if (c->flags & ASM_CHART_MIRROR) u = -u;
    double ct = cos(c->pose_theta), st = sin(c->pose_theta);
    *gx = ct * u - st * v + c->pose_x;
    *gy = st * u + ct * v + c->pose_y;
}

typedef struct ArComp { int32_t id; double area; double lo0, lo1, hi0, hi1; double off_u, off_v; int drawn; } ArComp;

static int ar_cmp_comp_desc(const void *a, const void *b)
{
    const ArComp *x = a, *y = b;
    if (x->area > y->area) return -1;
    if (x->area < y->area) return 1;
    return (x->id > y->id) - (x->id < y->id);
}

int AsmReport_layout_png(Arena_T arena, const char *path,
                         const AsmChart *charts, size_t n,
                         const AsmConflictCell *cells, size_t n_cells, double cell_size,
                         int max_px, size_t max_components)
{
    Arena_Mark mark = Arena_save(arena);
    int32_t ncomp = 0;
    for (size_t i = 0; i < n; i++) if (charts[i].placed && charts[i].component + 1 > ncomp) ncomp = charts[i].component + 1;
    if (ncomp == 0) { Arena_restore(arena, mark); return -1; }
    ArComp *comp = ARENA_CALLOC(arena, (size_t)ncomp, sizeof(ArComp));
    for (int32_t c = 0; c < ncomp; c++) { comp[c].id = c; comp[c].lo0 = comp[c].lo1 = 1e300; comp[c].hi0 = comp[c].hi1 = -1e300; }
    for (size_t i = 0; i < n; i++) {
        const AsmChart *ch = &charts[i];
        if (!AsmChart_in_layout(ch) || !ch->placed || ch->component < 0) continue;
        ArComp *cp = &comp[ch->component];
        cp->area += ch->area3d;
        for (size_t k = 0; k < ch->nv; k++) {
            double gx, gy;
            double point[2]; AsmChart_point(ch,k,point); gx = point[0]; gy = point[1];
            if (gx < cp->lo0) cp->lo0 = gx; if (gx > cp->hi0) cp->hi0 = gx;
            if (gy < cp->lo1) cp->lo1 = gy; if (gy > cp->hi1) cp->hi1 = gy;
        }
    }
    ArComp *order = ARENA_ALLOC(arena, (size_t)ncomp * sizeof(ArComp));
    memcpy(order, comp, (size_t)ncomp * sizeof(ArComp));
    qsort(order, (size_t)ncomp, sizeof(ArComp), ar_cmp_comp_desc);
    /* stack components vertically, largest first */
    double gap = 24.0, width = 0.0, height = 0.0;
    size_t drawn = 0;
    for (int32_t k = 0; k < ncomp && drawn < max_components; k++) {
        ArComp *cp = &order[k];
        if (cp->area <= 0.0) continue;
        double w = cp->hi0 - cp->lo0, h = cp->hi1 - cp->lo1;
        comp[cp->id].off_u = -cp->lo0;
        comp[cp->id].off_v = height - cp->lo1;
        comp[cp->id].drawn = 1;
        height += h + gap;
        if (w > width) width = w;
        drawn++;
    }
    if (drawn == 0) { Arena_restore(arena, mark); return -1; }
    double span = fmax(width, height);
    double scale = (double)(max_px - 8) / span;
    if (scale > 1.0) scale = 1.0;
    int W = (int)ceil(width * scale) + 8, H = (int)ceil(height * scale) + 8;
    if (W < 8) W = 8;
    if (H < 8) H = 8;
    AsmCanvas cv = AsmCanvas_new(arena, W, H, 24, 24, 28);
    for (size_t i = 0; i < n; i++) {
        const AsmChart *ch = &charts[i];
        if (!AsmChart_in_layout(ch) || !ch->placed || ch->component < 0 || !comp[ch->component].drawn) continue;
        const ArComp *cp = &comp[ch->component];
        uint8_t rgb[3];
        ar_chart_colour(ch, rgb);
        for (size_t f = 0; f < ch->nf; f++) {
            const int32_t *fv = &ch->faces[f*3];
            double p[6];
            for (int k = 0; k < 3; k++) {
                double gx, gy;
                double point[2]; AsmChart_point(ch,(size_t)fv[k],point); gx = point[0]; gy = point[1];
                p[k*2] = (gx + cp->off_u) * scale + 4.0;
                p[k*2+1] = (gy + cp->off_v) * scale + 4.0;
            }
            AsmCanvas_tri(&cv, p[0], p[1], p[2], p[3], p[4], p[5], rgb);
        }
    }
    for (size_t k = 0; k < n_cells; k++) {
        const AsmConflictCell *cl = &cells[k];
        if (cl->component < 0 || cl->component >= ncomp || !comp[cl->component].drawn) continue;
        const ArComp *cp = &comp[cl->component];
        double x0 = (cl->cx * cell_size + cp->off_u) * scale + 4.0, y0 = (cl->cy * cell_size + cp->off_v) * scale + 4.0;
        double x1 = ((cl->cx + 1) * cell_size + cp->off_u) * scale + 4.0, y1 = ((cl->cy + 1) * cell_size + cp->off_v) * scale + 4.0;
        uint8_t red[3] = { (uint8_t)(120 + 135 * cl->w), 0, 0 };
        AsmCanvas_rect(&cv, x0, y0, x1, y1, red);
    }
    int rc = AsmCanvas_write(&cv, path);
    Arena_restore(arena, mark);
    return rc;
}

/* ---- pose ledger ------------------------------------------------------------------ */

int AsmReport_poses_csv(const AsmRun *run, const char *path, const char *stage)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    fputs("stage,chart,component,u,v,theta,mirror,cu,cv,area,cube,cz,cy,cx\n", fp);
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c)) continue;
        double su = 0.0, sv = 0.0, p[2] = { 0.0, 0.0 };
        if (c->nv) {
            for (size_t v = 0; v < c->nv; v++) { su += c->uv[2*v]; sv += c->uv[2*v+1]; }
            /* the layout point of the mean uv: the pose is affine, so this is the layout centroid */
            double u = su / (double)c->nv, w = sv / (double)c->nv;
            if (c->flags & ASM_CHART_MIRROR) u = -u;
            double ct = cos(c->pose_theta), sn = sin(c->pose_theta);
            p[0] = ct*u - sn*w + c->pose_x; p[1] = sn*u + ct*w + c->pose_y;
        }
        fprintf(fp, "%s,%zu,%d,%.9g,%.9g,%.12g,%d,%.3f,%.3f,%.6e,%d,%.1f,%.1f,%.1f\n", stage, i, c->component,
                c->pose_x, c->pose_y, c->pose_theta, (c->flags & ASM_CHART_MIRROR) != 0, p[0], p[1], c->area3d,
                c->cube, c->centroid[0], c->centroid[1], c->centroid[2]);
    }
    return fclose(fp) ? -1 : 0;
}
