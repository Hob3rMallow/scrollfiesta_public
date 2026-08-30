/* uv_compact.c -- monotone u compaction of a parameterized atlas.
 *
 * The unwrapped atlas is a spiral laid flat, so material legitimately spans the
 * whole u range; the empty space is not a misplaced-chart problem a winding
 * search can fix.  Measured on a PHerc0139 10x10x10: 4,111 of 20,963 columns
 * are under 10% filled, scattered over 35 contiguous spans, and they hold 0.9%
 * of the material.  They are arc positions where the sheet is genuinely sparse.
 *
 * This squeezes those spans by a fixed factor with a MONOTONE piecewise-linear
 * map, so u order -- and therefore the sheet's continuity -- is preserved
 * exactly and no discontinuity can be introduced.  Dense columns are untouched,
 * so the distortion is confined to the sparse spans: at threshold 0.05 and
 * factor 2 that is 0.40% of the material, well under the p99 distortion cut.
 *
 * Coverage is measured on a column grid of the given du, from the mesh's own
 * UVs, so no baked raster is needed.
 *
 *   uv_compact <in.vmesh|obj> <out.obj> [--du F=2] [--dv F=1]
 *              [--thresh F=0.05] [--factor F=2]
 *   uv_compact --selftest
 * Exit: 0 ok / selftest pass, 1 IO, 2 usage, 3 selftest fail.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/mesh_bin.h"
#include "../common/obj_io.h"

/* Build the monotone map: dense columns keep unit width, sparse ones shrink to
 * 1/factor.  Returns the new u for a given old u by prefix-summing widths. */
static void uc_build_map(const uint8_t *dense, size_t ncol, double du,
                         double factor, double u0, double *edge)
{
    double at = 0.0;
    for (size_t i = 0; i < ncol; i++) {
        edge[i] = at;
        at += dense[i] ? du : du / factor;
    }
    edge[ncol] = at;
    (void)u0;
}

static double uc_map_u(double u, double u0, double du, size_t ncol,
                       const double *edge)
{
    double t = (u - u0) / du;
    if (t <= 0.0) return edge[0] + t * du;          /* linear outside, monotone */
    if (t >= (double)ncol) return edge[ncol] + (t - (double)ncol) * du;
    size_t i = (size_t)t;
    double f = t - (double)i;
    return edge[i] + f * (edge[i + 1] - edge[i]);
}

static int uc_run(const char *in, const char *out, double du, double dv,
                  double thresh, double factor)
{
    MeshBinData m;
    char binpath[2048];
    if (MeshBin_companion_path(in, binpath, sizeof binpath) != 0 ||
        MeshBin_read_malloc(binpath, &m) != 0) {
        fprintf(stderr, "uv_compact: cannot read %s\n", binpath);
        return 1;
    }
    if (m.uv == NULL) {
        fprintf(stderr, "uv_compact: %s has no UVs\n", binpath);
        MeshBin_dispose(&m);
        return 1;
    }
    double u0 = m.uv[0], u1 = m.uv[0], v0 = m.uv[1], v1 = m.uv[1];
    for (size_t i = 1; i < m.nv; i++) {
        double u = m.uv[i * 2], v = m.uv[i * 2 + 1];
        if (u < u0) u0 = u;
        if (u > u1) u1 = u;
        if (v < v0) v0 = v;
        if (v > v1) v1 = v;
    }
    size_t ncol = (size_t)((u1 - u0) / du) + 1;
    size_t nrow = (size_t)((v1 - v0) / dv) + 1;
    if (ncol == 0 || nrow == 0) { MeshBin_dispose(&m); return 1; }
    /* Column occupancy from the vertices themselves: a column counts a row as
     * covered if any vertex lands in that cell.  Coarse, but the threshold is
     * a fraction of the FULL v height, so only near-empty columns are caught. */
    uint8_t *seen = (uint8_t *)calloc(ncol * nrow, 1);
    if (seen == NULL) { MeshBin_dispose(&m); return 1; }
    for (size_t i = 0; i < m.nv; i++) {
        size_t ci = (size_t)((m.uv[i * 2] - u0) / du);
        size_t ri = (size_t)((m.uv[i * 2 + 1] - v0) / dv);
        if (ci < ncol && ri < nrow) seen[ci * nrow + ri] = 1;
    }
    uint8_t *dense = (uint8_t *)calloc(ncol, 1);
    double *edge = (double *)calloc(ncol + 1, sizeof(double));
    if (dense == NULL || edge == NULL) {
        free(seen); free(dense); free(edge); MeshBin_dispose(&m); return 1;
    }
    size_t nsparse = 0;
    for (size_t c = 0; c < ncol; c++) {
        size_t n = 0;
        for (size_t r = 0; r < nrow; r++) n += seen[c * nrow + r];
        dense[c] = (double)n >= thresh * (double)nrow;
        if (!dense[c]) nsparse++;
    }
    free(seen);
    uc_build_map(dense, ncol, du, factor, u0, edge);
    for (size_t i = 0; i < m.nv; i++)
        m.uv[i * 2] = (float)(u0 + uc_map_u(m.uv[i * 2], u0, du, ncol, edge));
    double neww = edge[ncol];
    fprintf(stderr,
            "uv_compact: %zu columns, %zu sparse (<%.0f%% of v) squeezed %.1fx; "
            "u span %.0f -> %.0f (%.1f%%)\n",
            ncol, nsparse, thresh * 100.0, factor, (u1 - u0), neww,
            100.0 * neww / (u1 - u0));
    char outbin[2048];
    int rc = 0;
    if (ObjIO_write_uv(out, m.verts, m.nv, m.faces, m.nf, m.uv) != 0 ||
        MeshBin_companion_path(out, outbin, sizeof outbin) != 0 ||
        MeshBin_write(outbin, m.verts, m.nv, m.faces, m.nf, m.uv) != 0) {
        fprintf(stderr, "uv_compact: cannot write %s and its VMESH companion\n",
                out);
        rc = 1;
    }
    free(dense); free(edge);
    MeshBin_dispose(&m);
    return rc;
}

static int uc_selftest(void)
{
    int fail = 0;
    /* three columns, middle one sparse: the map must stay strictly increasing
     * and shrink only the middle. */
    uint8_t dense[3] = { 1, 0, 1 };
    double edge[4];
    uc_build_map(dense, 3, 2.0, 2.0, 0.0, edge);
    fail |= !(edge[0] == 0.0 && edge[1] == 2.0 && edge[2] == 3.0 &&
              edge[3] == 5.0);
    for (int i = 0; i < 3; i++) fail |= !(edge[i] < edge[i + 1]);
    double prev = -1e300;
    for (double u = -1.0; u <= 7.0; u += 0.25) {
        double y = uc_map_u(u, 0.0, 2.0, 3, edge);
        fail |= !(y > prev);           /* strictly monotone everywhere */
        prev = y;
    }
    fprintf(stderr, "uv_compact selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}

int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--selftest")) return uc_selftest() ? 3 : 0;
    if (argc < 3) {
        fprintf(stderr, "usage: %s <in.vmesh|obj> <out.obj> [--du F] [--dv F]"
                        " [--thresh F] [--factor F]\n       %s --selftest\n",
                argv[0], argv[0]);
        return 2;
    }
    double du = 2.0, dv = 1.0, thresh = 0.05, factor = 2.0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--du") && i + 1 < argc) du = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dv") && i + 1 < argc) dv = atof(argv[++i]);
        else if (!strcmp(argv[i], "--thresh") && i + 1 < argc) thresh = atof(argv[++i]);
        else if (!strcmp(argv[i], "--factor") && i + 1 < argc) factor = atof(argv[++i]);
        else { fprintf(stderr, "uv_compact: unknown arg %s\n", argv[i]); return 2; }
    }
    if (!(du > 0.0) || !(dv > 0.0) || !(factor >= 1.0) ||
        !(thresh >= 0.0 && thresh < 1.0)) {
        fprintf(stderr, "uv_compact: bad parameters\n");
        return 2;
    }
    return uc_run(argv[1], argv[2], du, dv, thresh, factor);
}
