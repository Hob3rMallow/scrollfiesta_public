/* scroll_axis_track -- derive the scroll's umbilicus CURVE from a per-cube
 * mesh pile (or one world-frame vmesh), streaming one cube at a time.
 *
 *   scroll_axis_track <dump_dir|mesh.vmesh> <out_dir> --tag NAME
 *                     [--slab 64] [--manifest Y X]
 *                     [--subgrid z0 z1 y0 y1 x0 x1]
 *                     [--prior-table CSV]   pass 2: keep only samples within
 *                                           AXIS_TRACE_RESERVOIR_R of the prior curve
 *   scroll_axis_track --selftest
 *
 * Two passes on a big grid: pass 1 with the uniform per-slab reservoir finds
 * the curve coarsely (a 21x21 slab spreads the reservoir thin, ~7k samples
 * in the local window); pass 2 streams again keeping a DENSE local
 * reservoir around the pass-1 curve.  Memory stays one cube + the slabs.
 *
 * Outputs under <out_dir>:
 *   <tag>.axis.csv       the table (z,y,x,... with a provenance header); the
 *                        first three columns load through AxisWarp_load_csv
 *   <tag>.axis.json      every row with both estimators' diagnostics
 *   <tag>_contact.png    per-slab occupancy + polar panels, z order
 *   <tag>_slab_z#####.png the same panels one file per slab
 *
 * Doctrine: the manifest point is REPORTED, never used as the seed; the
 * seed is derived on the bottom slab from a grid of starts.  A slab where
 * the primary position estimator abstains is interpolated (short gaps) or
 * dropped, never held at the previous value; the secondary is diagnostic. */

#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_pile.h"
#include "../common/png_write.h"
#include "../common/pipeline_constants.h"
#include "../whole/axis_track.h"
#include "../whole/axis_warp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SAT_MAX_PILE 65536
#define SAT_MAX_SLABS 512
#define SAT_MAX_CLUSTERS 32

typedef struct { long sg[6]; int armed; } SatSubgrid;

static int sat_filter(void *ctx, long oz, long oy, long ox)
{
    const SatSubgrid *g = (const SatSubgrid *)ctx;
    if (!g->armed) return 1;
    return oz >= g->sg[0] && oz <= g->sg[1] && oy >= g->sg[2] &&
           oy <= g->sg[3] && ox >= g->sg[4] && ox <= g->sg[5];
}

static int sat_ends_with(const char *s, const char *suf)
{
    size_t ls = strlen(s), lf = strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

/* ---- drawing --------------------------------------------------------------*/

static void sat_px(uint8_t *rgb, int w, int h, int x, int y,
                   uint8_t r, uint8_t g, uint8_t b)
{
    if (x < 0 || y < 0 || x >= w || y >= h) return;
    rgb[((size_t)y * (size_t)w + (size_t)x) * 3 + 0] = r;
    rgb[((size_t)y * (size_t)w + (size_t)x) * 3 + 1] = g;
    rgb[((size_t)y * (size_t)w + (size_t)x) * 3 + 2] = b;
}

static void sat_cross(uint8_t *rgb, int w, int h, int x, int y, int arm,
                      uint8_t r, uint8_t g, uint8_t b)
{
    for (int d = -arm; d <= arm; d++) {
        sat_px(rgb, w, h, x + d, y, r, g, b);
        sat_px(rgb, w, h, x, y + d, r, g, b);
    }
}

static void sat_circle(uint8_t *rgb, int w, int h, int x, int y, int rad,
                       uint8_t r, uint8_t g, uint8_t b)
{
    for (int a = 0; a < 64; a++) {
        double t = 2.0 * M_PI * a / 64.0;
        sat_px(rgb, w, h, x + (int)lround(rad * cos(t)), y + (int)lround(rad * sin(t)), r, g, b);
    }
}

static void sat_square(uint8_t *rgb, int w, int h, int x, int y, int half,
                       uint8_t r, uint8_t g, uint8_t b)
{
    for (int d = -half; d <= half; d++) {
        sat_px(rgb, w, h, x + d, y - half, r, g, b);
        sat_px(rgb, w, h, x + d, y + half, r, g, b);
        sat_px(rgb, w, h, x - half, y + d, r, g, b);
        sat_px(rgb, w, h, x + half, y + d, r, g, b);
    }
}

/* Occupancy panel (S x S over the global bbox) stacked over a polar panel
 * (180 theta columns x 300 r rows, 2 vox per row) about the fused centre. */
#define SAT_POLAR_W 180
#define SAT_POLAR_H 300

static void sat_render_slab(const AxisSlab *s, const AxisRow *row,
                            const AxisRow *prev_row, double ylo, double xlo,
                            double cell, int S, int have_manifest,
                            double man_y, double man_x, uint8_t *rgb, int W,
                            int H, int ox0, int oy0)
{
    uint16_t *occ = (uint16_t *)calloc((size_t)S * (size_t)S, sizeof *occ);
    uint16_t *pol = (uint16_t *)calloc((size_t)SAT_POLAR_W * SAT_POLAR_H, sizeof *pol);
    double cy = row->y, cx = row->x;
    int use_centre = row->status != AXIS_ROW_ABSTAIN;
    if (occ == NULL || pol == NULL) { free(occ); free(pol); return; }
    if (!use_centre) { cy = row->y2; cx = row->x2; }
    for (size_t i = 0; i < s->nv; i++) {
        int px = (int)((s->vx[i] - xlo) / cell), py = (int)((s->vy[i] - ylo) / cell);
        double dy = s->vy[i] - cy, dx = s->vx[i] - cx;
        double r = sqrt(dy * dy + dx * dx);
        int pr = (int)(r / 2.0);
        int pt = (int)((atan2(dx, dy) + M_PI) / (2.0 * M_PI) * SAT_POLAR_W);
        if (px >= 0 && py >= 0 && px < S && py < S && occ[(size_t)py * S + px] < 65535)
            occ[(size_t)py * S + px]++;
        if (pt < 0) pt = 0;
        if (pt >= SAT_POLAR_W) pt = SAT_POLAR_W - 1;
        if (pr >= 0 && pr < SAT_POLAR_H && pol[(size_t)pr * SAT_POLAR_W + pt] < 65535)
            pol[(size_t)pr * SAT_POLAR_W + pt]++;
    }
    for (int y = 0; y < S; y++)
        for (int x = 0; x < S; x++) {
            double v = occ[(size_t)y * S + x];
            uint8_t g = (uint8_t)(v > 0 ? fmin(255.0, 60.0 + 40.0 * log(1.0 + v)) : 0);
            sat_px(rgb, W, H, ox0 + x, oy0 + y, g, g, g);
        }
    for (int y = 0; y < SAT_POLAR_H; y++)
        for (int x = 0; x < SAT_POLAR_W; x++) {
            double v = pol[(size_t)y * SAT_POLAR_W + x];
            uint8_t g = (uint8_t)(v > 0 ? fmin(255.0, 60.0 + 40.0 * log(1.0 + v)) : 0);
            sat_px(rgb, W, H, ox0 + x, oy0 + S + 4 + y, g, g, g);
        }
    /* markers: manifest grey square, primary red cross, secondary cyan
     * circle, fused yellow cross, previous slab's centre magenta dot */
    if (have_manifest)
        sat_square(rgb, W, H, ox0 + (int)((man_x - xlo) / cell),
                   oy0 + (int)((man_y - ylo) / cell), 5, 160, 160, 160);
    if (row->n_in > 0)
        sat_cross(rgb, W, H, ox0 + (int)((row->x1 - xlo) / cell),
                  oy0 + (int)((row->y1 - ylo) / cell), 6, 255, 60, 60);
    if (row->ridges > 0)
        sat_circle(rgb, W, H, ox0 + (int)((row->x2 - xlo) / cell),
                   oy0 + (int)((row->y2 - ylo) / cell), 7, 60, 220, 255);
    if (use_centre)
        sat_cross(rgb, W, H, ox0 + (int)((row->x - xlo) / cell),
                  oy0 + (int)((row->y - ylo) / cell), 9, 255, 230, 40);
    if (prev_row != NULL && prev_row->status != AXIS_ROW_ABSTAIN)
        sat_circle(rgb, W, H, ox0 + (int)((prev_row->x - xlo) / cell),
                   oy0 + (int)((prev_row->y - ylo) / cell), 2, 255, 60, 255);
    /* status bar under the polar panel: yellow lock, cyan interp, orange
     * extrap, red abstain */
    {
        uint8_t r = 255, g = 60, b = 60;
        if (row->status == AXIS_ROW_LOCK) { r = 255; g = 230; b = 40; }
        else if (row->status == AXIS_ROW_INTERP) { r = 60; g = 220; b = 255; }
        else if (row->status == AXIS_ROW_EXTRAP) { r = 255; g = 140; b = 40; }
        for (int x = 0; x < S; x++)
            for (int y = 0; y < 3; y++)
                sat_px(rgb, W, H, ox0 + x, oy0 + S + 4 + SAT_POLAR_H + 2 + y, r, g, b);
    }
    free(occ); free(pol);
}

static const char *sat_status(int st)
{
    switch (st) {
    case AXIS_ROW_LOCK: return "lock";
    case AXIS_ROW_INTERP: return "interp";
    case AXIS_ROW_EXTRAP: return "extrap";
    default: return "abstain";
    }
}

/* ---- main ------------------------------------------------------------------*/

static void sat_usage(const char *argv0)
{
    fprintf(stderr,
            "usage: %s <dump_dir|mesh.vmesh> <out_dir> --tag NAME [--slab 64]\n"
            "          [--manifest Y X] [--subgrid z0 z1 y0 y1 x0 x1]\n"
            "       %s --selftest\n", argv0, argv0);
}

static int sat_selftest(void)
{
    int fails = 0;
    fails += MeshPile_selftest();
    fails += AxisWarp_selftest();
    fails += AxisTrack_selftest();
    fprintf(stderr, "scroll_axis_track --selftest: %s (%d failures)\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails == 0 ? 0 : 1;
}

int main(int argc, char **argv)
{
    const char *input = NULL, *out_dir = NULL, *tag = "scroll";
    double slab_vox = AXIS_TRACE_SLAB_VOX;
    int have_manifest = 0, dump_polar = 0, have_prior = 0;
    double man_y = 0.0, man_x = 0.0, dump_y = 0.0, dump_x = 0.0;
    const char *prior_path = NULL;
    AxisWarp prior;
    size_t kept_v = 0, kept_f = 0;
    AxisWarp_init(&prior);
    SatSubgrid grid;
    MeshPileEntry *pile = NULL;
    size_t npile = 0;
    double zmin = 1e30, zmax = -1e30;
    size_t nslab = 0;
    AxisSlab *slabs = NULL;
    AxisRow *rows = NULL, *rows_b = NULL;
    size_t total_v = 0, total_f = 0;
    char path[MESH_PILE_MAX_PATH];
    clock_t t0 = clock();
    memset(&grid, 0, sizeof grid);
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return sat_selftest();
    if (argc < 3) { sat_usage(argv[0]); return 1; }
    input = argv[1]; out_dir = argv[2];
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--tag") == 0 && i + 1 < argc) tag = argv[++i];
        else if (strcmp(argv[i], "--slab") == 0 && i + 1 < argc) slab_vox = atof(argv[++i]);
        else if (strcmp(argv[i], "--manifest") == 0 && i + 2 < argc) {
            man_y = atof(argv[++i]); man_x = atof(argv[++i]); have_manifest = 1;
        } else if (strcmp(argv[i], "--prior-table") == 0 && i + 1 < argc) {
            prior_path = argv[++i];
        } else if (strcmp(argv[i], "--dump-polar") == 0 && i + 2 < argc) {
            dump_y = atof(argv[++i]); dump_x = atof(argv[++i]); dump_polar = 1;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            AxisTrack_verbose = 1;
        } else if (strcmp(argv[i], "--trace") == 0) {
            AxisTrack_verbose = 2;
        } else if (strcmp(argv[i], "--subgrid") == 0 && i + 6 < argc) {
            for (int k = 0; k < 6; k++) grid.sg[k] = strtol(argv[i + 1 + k], NULL, 10);
            i += 6; grid.armed = 1;
        } else { sat_usage(argv[0]); return 1; }
    }
    if (!(slab_vox >= 16.0)) { fprintf(stderr, "slab must be >= 16 vox\n"); return 1; }
    if (prior_path != NULL) {
        if (AxisWarp_load_csv(&prior, prior_path) != 0 || !AxisWarp_valid(&prior)) {
            fprintf(stderr, "scroll_axis_track: cannot load the prior table %s\n", prior_path);
            return 1;
        }
        have_prior = 1;
        fprintf(stderr, "scroll_axis_track: pass 2 -- reservoir local to the prior curve %s (%zu rows, within %.0f vox)\n",
                prior_path, prior.n, AXIS_TRACE_RESERVOIR_R);
    }
    ves_mkdir(out_dir);

    /* 1. the pile */
    pile = (MeshPileEntry *)malloc(SAT_MAX_PILE * sizeof *pile);
    if (pile == NULL) return 1;
    if (sat_ends_with(input, ".vmesh")) {
        memset(&pile[0], 0, sizeof pile[0]);
        snprintf(pile[0].path, sizeof pile[0].path, "%s", input);
        npile = 1;
    } else if (MeshPile_scan(input, pile, SAT_MAX_PILE, &npile, sat_filter, &grid) != 0) {
        fprintf(stderr, "scroll_axis_track: cannot scan %s\n", input);
        return 1;
    }
    if (npile == 0) { fprintf(stderr, "scroll_axis_track: empty pile\n"); return 1; }
    for (size_t i = 0; i < npile; i++) {
        if (!pile[i].has_id) continue;
        if ((double)pile[i].oz < zmin) zmin = (double)pile[i].oz;
        if ((double)pile[i].oz + 128.0 > zmax) zmax = (double)pile[i].oz + 128.0;
    }
    if (!(zmax > zmin)) {
        /* a plain mesh: one pass over it for the z-range */
        Arena_T ar = Arena_new();
        MeshBinData md;
        const float *v = NULL; size_t nv = 0;
        memset(&md, 0, sizeof md);
        if (MeshBin_read_arena(ar, pile[0].path, &md) != 0 || md.nv == 0) {
            fprintf(stderr, "scroll_axis_track: cannot read %s\n", pile[0].path);
            return 1;
        }
        v = md.verts; nv = md.nv;
        for (size_t i = 0; i < nv; i++) {
            if (v[i * 3] < zmin) zmin = v[i * 3];
            if (v[i * 3] > zmax) zmax = v[i * 3];
        }
        zmax += 1.0;
        Arena_dispose(&ar);
    }
    nslab = (size_t)ceil((zmax - zmin) / slab_vox);
    if (nslab == 0 || nslab > SAT_MAX_SLABS) {
        fprintf(stderr, "scroll_axis_track: %zu slabs over z [%.0f, %.0f) is out of range\n",
                nslab, zmin, zmax);
        return 1;
    }
    fprintf(stderr, "scroll_axis_track: %zu pile entries, z [%.0f, %.0f) -> %zu slabs of %.0f vox "
            "(reservoir %d faces + %d verts per slab, ~%.0f MB)\n",
            npile, zmin, zmax, nslab, slab_vox, AXIS_TRACE_RESERVOIR_FACES,
            AXIS_TRACE_RESERVOIR_VERTS,
            (double)nslab * (AXIS_TRACE_RESERVOIR_FACES * 32.0 + AXIS_TRACE_RESERVOIR_VERTS * 16.0) / 1048576.0);
    slabs = (AxisSlab *)calloc(nslab, sizeof *slabs);
    rows = (AxisRow *)calloc(nslab, sizeof *rows);
    rows_b = (AxisRow *)calloc(nslab, sizeof *rows_b);
    if (slabs == NULL || rows == NULL || rows_b == NULL) return 1;
    for (size_t k = 0; k < nslab; k++) {
        if (AxisSlab_init(&slabs[k], zmin + slab_vox * (double)k,
                          zmin + slab_vox * (double)(k + 1),
                          AXIS_TRACE_RESERVOIR_FACES, AXIS_TRACE_RESERVOIR_VERTS,
                          (uint32_t)(0x1234567u + 7919u * (uint32_t)k)) != 0) {
            fprintf(stderr, "scroll_axis_track: out of memory for slab %zu\n", k);
            return 1;
        }
        rows[k].z = zmin + slab_vox * ((double)k + 0.5);
    }

    /* 2. stream the pile into the slab reservoirs */
    for (size_t i = 0; i < npile; i++) {
        Arena_T ar = Arena_new();
        MeshBinData md;
        const float *v = NULL; const int32_t *f = NULL; size_t nv = 0, nf = 0;
        memset(&md, 0, sizeof md);
        if (MeshBin_read_arena(ar, pile[i].path, &md) != 0) {
            fprintf(stderr, "  skip (unreadable) %s\n", pile[i].path);
            Arena_dispose(&ar);
            continue;
        }
        v = md.verts; f = md.faces; nv = md.nv; nf = md.nf;
        for (size_t j = 0; j < nv; j++) {
            long k = (long)floor(((double)v[j * 3] - zmin) / slab_vox);
            if (k < 0 || (size_t)k >= nslab) continue;
            if (have_prior) {
                double py = 0.0, px = 0.0;
                AxisWarp_eval(&prior, (double)v[j * 3], &py, &px);
                if (hypot((double)v[j * 3 + 1] - py, (double)v[j * 3 + 2] - px) > AXIS_TRACE_RESERVOIR_R)
                    continue;
                kept_v++;
            }
            AxisSlab_add_vertex(&slabs[k], &v[j * 3]);
        }
        for (size_t j = 0; j < nf; j++) {
            const float *a = &v[(size_t)f[j * 3 + 0] * 3];
            const float *b = &v[(size_t)f[j * 3 + 1] * 3];
            const float *c = &v[(size_t)f[j * 3 + 2] * 3];
            double zc = ((double)a[0] + b[0] + c[0]) / 3.0;
            long k = (long)floor((zc - zmin) / slab_vox);
            if (k < 0 || (size_t)k >= nslab) continue;
            if (have_prior) {
                double py = 0.0, px = 0.0;
                double yc = ((double)a[1] + b[1] + c[1]) / 3.0;
                double xc = ((double)a[2] + b[2] + c[2]) / 3.0;
                AxisWarp_eval(&prior, zc, &py, &px);
                if (hypot(yc - py, xc - px) > AXIS_TRACE_RESERVOIR_R) continue;
                kept_f++;
            }
            AxisSlab_add_face(&slabs[k], a, b, c);
        }
        total_v += nv; total_f += nf;
        Arena_dispose(&ar);
        if ((i + 1) % 200 == 0 || i + 1 == npile)
            fprintf(stderr, "  streamed %zu/%zu entries (%zu verts, %zu faces) %.0fs\n",
                    i + 1, npile, total_v, total_f,
                    (double)(clock() - t0) / CLOCKS_PER_SEC);
    }

    /* debug: the polar raster of every slab about a given centre */
    if (dump_polar) {
        for (size_t k = 0; k < nslab; k++) {
            snprintf(path, sizeof path, "%s/%s_polar_z%05.0f.png", out_dir, tag, rows[k].z);
            AxisTrack_polar_debug(&slabs[k], dump_y, dump_x, 600, path);
        }
    }

    /* 3. seed on the lowest slab with material */
    double seed_y = 0.0, seed_x = 0.0;
    double cl_y[SAT_MAX_CLUSTERS], cl_x[SAT_MAX_CLUSTERS], cl_d1[SAT_MAX_CLUSTERS];
    size_t cl_m[SAT_MAX_CLUSTERS], cl_r[SAT_MAX_CLUSTERS], ncl = 0, seed_slab = 0;
    int seeded = 0;
    if (have_prior)
        fprintf(stderr, "  pass 2 reservoir kept %zu verts, %zu faces near the prior curve\n", kept_v, kept_f);
    for (size_t k = 0; k < nslab && !seeded; k++) {
        double wy = 0.0, wx = 0.0, wr = 0.0;
        if (slabs[k].nv < (have_prior ? 5000u : 20000u)) continue;
        if (have_prior) { AxisWarp_eval(&prior, rows[k].z, &wy, &wx); wr = 1.5 * AXIS_TRACE_SEED_GRID; }
        if (AxisTrack_seed_in(&slabs[k], wy, wx, wr, &seed_y, &seed_x, cl_y, cl_x, cl_m, cl_d1, cl_r,
                              SAT_MAX_CLUSTERS, &ncl) == 0) { seeded = 1; seed_slab = k; }
    }
    if (!seeded) {
        fprintf(stderr, "scroll_axis_track: no slab yielded a seed (abstain everywhere)\n");
        return 2;
    }
    fprintf(stderr, "  seed: slab %zu (z %.0f) -> (y %.1f, x %.1f) from %zu cluster(s)\n",
            seed_slab, rows[seed_slab].z, seed_y, seed_x, ncl);
    for (size_t c = 0; c < ncl; c++)
        fprintf(stderr, "    cluster %zu: (%.1f, %.1f) members=%zu void=%.2f resid_p50=%.1f%s\n",
                c, cl_y[c], cl_x[c], cl_m[c], (double)cl_r[c] / 1000.0, cl_d1[c],
                (double)cl_r[c] / 1000.0 >= AXIS_TRACE_VOID_MAX ? " (no void: not the umbilicus)" : "");
    if (have_manifest)
        fprintf(stderr, "  manifest point (%.0f, %.0f) is %.1f vox from the derived seed%s\n",
                man_y, man_x, hypot(man_y - seed_y, man_x - seed_x),
                hypot(man_y - seed_y, man_x - seed_x) > AXIS_TABLE_MANIFEST_RETIRE
                    ? " -> RETIRE it" : "");

    /* 4. forward pass from the seed slab upward, then downward */
    {
        double py = seed_y, px = seed_x;
        double hist_z[3], hist_y[3], hist_x[3];
        size_t nh = 0;
        for (int dir = 1; dir >= -1; dir -= 2) {
            py = seed_y; px = seed_x; nh = 0;
            for (long k = (long)seed_slab; k >= 0 && (size_t)k < nslab; k += dir) {
                AxisRow *r = &rows[k];
                double prior_y = py, prior_x = px;
                if (dir == -1 && (size_t)k == seed_slab) continue;
                if (nh >= 2) {   /* linear trend of the last locked slabs */
                    double dz = hist_z[0] - hist_z[nh - 1];
                    if (fabs(dz) > 1e-9) {
                        double sy = (hist_y[0] - hist_y[nh - 1]) / dz;
                        double sx = (hist_x[0] - hist_x[nh - 1]) / dz;
                        prior_y = hist_y[0] + sy * (r->z - hist_z[0]);
                        prior_x = hist_x[0] + sx * (r->z - hist_z[0]);
                    }
                }
                AxisTrack_primary(&slabs[k], prior_y, prior_x, r);
                AxisTrack_secondary(&slabs[k], r->primary_ok ? r->y1 : prior_y,
                                    r->primary_ok ? r->x1 : prior_x, r);
                fprintf(stderr, "  z %6.0f: primary (%.1f, %.1f) n=%zu cond=%.2f sectors=%d void=%.2f resid=%.1f jump=%.1f %s | "
                        "secondary (%.1f, %.1f) ridges=%zu d1=%.1f pitch=%+.2f rwall=%.0f %s\n",
                        r->z, r->y1, r->x1, r->n_in, r->cond, r->sectors, r->void_ratio, r->resid_p50, r->jump,
                        r->primary_ok ? "ok" : "ABSTAIN", r->y2, r->x2, r->ridges, r->d1,
                        r->pitch, r->r_wall, r->secondary_ok ? "ok" : "ABSTAIN");
                if (r->primary_ok) {
                    for (size_t h = (nh < 3 ? nh : 2); h > 0; h--) {
                        hist_z[h] = hist_z[h - 1]; hist_y[h] = hist_y[h - 1]; hist_x[h] = hist_x[h - 1];
                    }
                    hist_z[0] = r->z; hist_y[0] = r->y1; hist_x[0] = r->x1;
                    if (nh < 3) nh++;
                    py = r->y1; px = r->x1;
                } else { py = prior_y; px = prior_x; }
            }
        }
    }
    /* 5. backward consistency: re-run from the top with the top's centre as
     * prior; a slab whose two passes disagree abstains */
    {
        size_t top = nslab;
        double py = 0.0, px = 0.0;
        int have = 0;
        for (size_t k = nslab; k-- > 0;)
            if (rows[k].primary_ok) { top = k; py = rows[k].y1; px = rows[k].x1; have = 1; break; }
        if (have) {
            size_t disagreements = 0;
            for (size_t k = top + 1; k-- > 0;) {
                AxisRow *b = &rows_b[k];
                memset(b, 0, sizeof *b);
                b->z = rows[k].z;
                AxisTrack_primary(&slabs[k], py, px, b);
                /* Backward consistency compares the position authority only.
                 * Re-running the large polar workspace here supplied no gate
                 * evidence (b->secondary is never consumed) and made dense
                 * pass-2 sweeps needlessly fragile. */
                if (b->primary_ok) { py = b->y1; px = b->x1; }
                if (rows[k].primary_ok && b->primary_ok &&
                    hypot(rows[k].y1 - b->y1, rows[k].x1 - b->x1) > AXIS_TRACE_FWD_BWD_AGREE) {
                    fprintf(stderr, "  z %6.0f: forward (%.1f, %.1f) vs backward (%.1f, %.1f) disagree -> ABSTAIN\n",
                            rows[k].z, rows[k].y1, rows[k].x1, b->y1, b->x1);
                    rows[k].primary_ok = 0;
                    disagreements++;
                }
            }
            fprintf(stderr, "  backward pass: %zu disagreement(s) over %zu slabs\n", disagreements, top + 1);
        }
    }

    /* 6. finish + sense */
    size_t nlock = AxisTrack_finish(rows, nslab);
    double sense_agree = 0.0;
    int sense = AxisTrack_sense(rows, nslab, &sense_agree);
    size_t n_secondary = 0;
    for (size_t k = 0; k < nslab; k++) if (rows[k].secondary_ok) n_secondary++;
    {
        double dmax = 0.0, smax = 0.0;
        size_t nint = 0, next = 0, nabs = 0;
        for (size_t k = 0; k < nslab; k++) {
            if (rows[k].status == AXIS_ROW_INTERP) nint++;
            else if (rows[k].status == AXIS_ROW_EXTRAP) next++;
            else if (rows[k].status == AXIS_ROW_ABSTAIN) nabs++;
            if (rows[k].status == AXIS_ROW_LOCK) {
                double d = hypot(rows[k].y - seed_y, rows[k].x - seed_x);
                if (d > dmax) dmax = d;
                for (size_t j = k + 1; j < nslab; j++)
                    if (rows[j].status == AXIS_ROW_LOCK) {
                        double sl = hypot(rows[j].y - rows[k].y, rows[j].x - rows[k].x) / (rows[j].z - rows[k].z);
                        if (sl > smax) smax = sl;
                        break;
                    }
            }
        }
        fprintf(stderr, "scroll_axis_track: %zu lock, %zu interp, %zu extrap, %zu abstain of %zu slabs; "
                "drift from seed max %.0f vox, max slope %.3f; winding sense %+d (agree %.2f)\n",
                nlock, nint, next, nabs, nslab, dmax, smax, sense, sense_agree);
    }

    /* 7. outputs */
    {
        char header[4096];
        int hl = 0;
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "tool: scroll_axis_track (src/whole/axis_track.c); tag %s\n", tag);
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "input: %s (%zu entries, %zu verts, %zu faces)\n", input, npile, total_v, total_f);
        if (have_prior)
            hl += snprintf(header + hl, sizeof header - (size_t)hl, "pass 2: reservoir within %.0f vox of the prior curve %s (%zu rows); kept %zu verts, %zu faces\n", AXIS_TRACE_RESERVOIR_R, prior_path, prior.n, kept_v, kept_f);
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "slab %.0f vox; z [%.0f, %.0f); reservoir %d faces + %d verts per slab\n", slab_vox, zmin, zmax, AXIS_TRACE_RESERVOIR_FACES, AXIS_TRACE_RESERVOIR_VERTS);
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "seed: slab z %.0f -> (y %.1f, x %.1f) from %zu cluster(s); the manifest point is never the seed\n", rows[seed_slab].z, seed_y, seed_x, ncl);
        for (size_t c = 0; c < ncl && hl < (int)sizeof header - 200; c++)
            hl += snprintf(header + hl, sizeof header - (size_t)hl, "seed cluster %zu: (%.1f, %.1f) members=%zu void=%.2f resid_p50=%.1f%s\n", c, cl_y[c], cl_x[c], cl_m[c], (double)cl_r[c] / 1000.0, cl_d1[c], (double)cl_r[c] / 1000.0 >= AXIS_TRACE_VOID_MAX ? " (no void)" : "");
        if (have_manifest)
            hl += snprintf(header + hl, sizeof header - (size_t)hl, "manifest point (%.0f, %.0f): %.1f vox from the derived seed (%s)\n", man_y, man_x, hypot(man_y - seed_y, man_x - seed_x), hypot(man_y - seed_y, man_x - seed_x) > AXIS_TABLE_MANIFEST_RETIRE ? "RETIRED" : "consistent");
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "constants: local_r %.0f trim %.2f irls %d min_samples %d min_cond %.2f max_resid %.0f max_jump %.0f fwd_bwd %.0f ridge_arc %.0f max_slope %.2f gap_rows %d extrap %.0f\n",
                       AXIS_TRACE_LOCAL_R, AXIS_TRACE_TRIM_FRAC, AXIS_TRACE_IRLS_ROUNDS, AXIS_TRACE_MIN_SAMPLES, AXIS_TRACE_MIN_COND, AXIS_TRACE_MAX_RESID_P50, AXIS_TRACE_MAX_JUMP, AXIS_TRACE_FWD_BWD_AGREE, AXIS_TRACE_RIDGE_MIN_ARC, AXIS_TABLE_MAX_SLOPE, AXIS_TABLE_MAX_GAP_ROWS, AXIS_TABLE_EXTRAP_VOX);
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "position authority: PRIMARY normal-line estimator (the polar harmonic check locked on %zu of %zu slabs; wraps are far from circular on this scroll)\n", n_secondary, nslab);
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "winding sense: %+d (agreement %.2f over slabs where the polar check locked; sign of the ridge slope; 0 = no evidence)\n", sense, sense_agree);
        hl += snprintf(header + hl, sizeof header - (size_t)hl, "locked %zu of %zu slabs; status per row below (lock|interp|extrap); abstained rows listed next", nlock, nslab);
        snprintf(path, sizeof path, "%s/%s.axis.csv", out_dir, tag);
        if (AxisTrack_write_csv(path, rows, nslab, header) != 0) {
            fprintf(stderr, "scroll_axis_track: cannot write %s\n", path);
            return 1;
        }
        fprintf(stderr, "  wrote %s\n", path);
    }
    {
        FILE *jf = NULL;
        snprintf(path, sizeof path, "%s/%s.axis.json", out_dir, tag);
        jf = fopen(path, "wb");
        if (jf != NULL) {
            fprintf(jf, "{\n  \"tag\": \"%s\",\n  \"slab_vox\": %.0f,\n  \"z_min\": %.0f,\n  \"z_max\": %.0f,\n"
                    "  \"seed\": {\"z\": %.0f, \"y\": %.2f, \"x\": %.2f, \"clusters\": %zu},\n"
                    "  \"winding_sense\": %d,\n  \"winding_sense_agree\": %.3f,\n  \"locked\": %zu,\n  \"rows\": [\n",
                    tag, slab_vox, zmin, zmax, rows[seed_slab].z, seed_y, seed_x, ncl, sense, sense_agree, nlock);
            for (size_t k = 0; k < nslab; k++) {
                const AxisRow *r = &rows[k];
                fprintf(jf, "    {\"z\": %.1f, \"status\": \"%s\", \"y\": %.2f, \"x\": %.2f, \"y1\": %.2f, \"x1\": %.2f, "
                        "\"y2\": %.2f, \"x2\": %.2f, \"n\": %zu, \"cond\": %.3f, \"sectors\": %d, \"void_ratio\": %.3f, "
                        "\"resid_p50\": %.2f, \"jump\": %.2f, "
                        "\"d1\": %.2f, \"agree\": %.2f, \"ridges\": %zu, \"pitch\": %.3f, \"r_wall\": %.1f}%s\n",
                        r->z, sat_status(r->status), r->y, r->x, r->y1, r->x1, r->y2, r->x2, r->n_in, r->cond,
                        r->sectors, r->void_ratio, r->resid_p50, r->jump, r->d1, r->agree, r->ridges, r->pitch, r->r_wall,
                        k + 1 < nslab ? "," : "");
            }
            fprintf(jf, "  ]\n}\n");
            fclose(jf);
            fprintf(stderr, "  wrote %s\n", path);
        }
    }
    /* contact sheet + per-slab panels */
    {
        double ylo = 1e30, yhi = -1e30, xlo = 1e30, xhi = -1e30, cell = 1.0;
        int S = 0, cols = 8, rowsn = 0, cellw = 0, cellh = 0, W = 0, H = 0;
        uint8_t *rgb = NULL;
        for (size_t k = 0; k < nslab; k++)
            for (size_t i = 0; i < slabs[k].nv; i++) {
                if (slabs[k].vy[i] < ylo) ylo = slabs[k].vy[i];
                if (slabs[k].vy[i] > yhi) yhi = slabs[k].vy[i];
                if (slabs[k].vx[i] < xlo) xlo = slabs[k].vx[i];
                if (slabs[k].vx[i] > xhi) xhi = slabs[k].vx[i];
            }
        if (yhi > ylo && xhi > xlo) {
            double span = fmax(yhi - ylo, xhi - xlo);
            cell = fmax(AXIS_TRACE_RASTER_VOX, span / 512.0);
            S = (int)ceil(span / cell) + 1;
            if (S < SAT_POLAR_W) S = SAT_POLAR_W;
            cellw = S + 8; cellh = S + 4 + SAT_POLAR_H + 2 + 3 + 8;
            if ((int)nslab < cols) cols = (int)nslab;
            rowsn = ((int)nslab + cols - 1) / cols;
            W = cols * cellw; H = rowsn * cellh;
            rgb = (uint8_t *)calloc((size_t)W * (size_t)H * 3, 1);
            if (rgb != NULL) {
                uint8_t *one = (uint8_t *)calloc((size_t)cellw * (size_t)cellh * 3, 1);
                for (size_t k = 0; k < nslab; k++) {
                    int cx0 = (int)(k % (size_t)cols) * cellw, cy0 = (int)(k / (size_t)cols) * cellh;
                    sat_render_slab(&slabs[k], &rows[k], k > 0 ? &rows[k - 1] : NULL, ylo, xlo,
                                    cell, S, have_manifest, man_y, man_x, rgb, W, H, cx0, cy0);
                    if (one != NULL) {
                        memset(one, 0, (size_t)cellw * (size_t)cellh * 3);
                        sat_render_slab(&slabs[k], &rows[k], k > 0 ? &rows[k - 1] : NULL, ylo, xlo,
                                        cell, S, have_manifest, man_y, man_x, one, cellw, cellh, 0, 0);
                        snprintf(path, sizeof path, "%s/%s_slab_z%05.0f.png", out_dir, tag, rows[k].z);
                        PngWrite_rgb(path, cellw, cellh, one);
                    }
                }
                snprintf(path, sizeof path, "%s/%s_contact.png", out_dir, tag);
                if (PngWrite_rgb(path, W, H, rgb) == 0)
                    fprintf(stderr, "  wrote %s (%dx%d, %d cols; occupancy panel %.1f vox/px over y[%.0f,%.0f] x[%.0f,%.0f]; polar panel 2 vox/row)\n",
                            path, W, H, cols, cell, ylo, yhi, xlo, xhi);
                free(one);
            }
            free(rgb);
        }
    }
    for (size_t k = 0; k < nslab; k++) AxisSlab_dispose(&slabs[k]);
    free(slabs); free(rows); free(rows_b); free(pile);
    AxisWarp_dispose(&prior);
    fprintf(stderr, "scroll_axis_track: done in %.0fs\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
    return nlock > 0 ? 0 : 2;
}
