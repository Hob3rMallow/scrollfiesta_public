/* swirl_slim.c -- recover a marbled umbilicus region by (1) cutting radius-gated
 * inter-wrap bridges, then (2) SLIM-refining (symmetric-Dirichlet) the UV of
 * every freed component.  The compressed centre -- where the spiral's turns are
 * crammed into a few U columns because the circumference vanishes -- is the high-
 * distortion region SLIM expands, so each turn gets its own columns and stops
 * sampling across neighbours (the marble).  Components pack largest-first; the
 * small leftovers land at the tail as kibble.
 *
 * usage:
 *   swirl_slim <in.vmesh> <out.vmesh> [options]
 *     --umbilicus Y X   source-space axis point (y,x); radius = |(y,x)-(Y,X)|
 *     --cut-gate F      drop faces whose vertex-radius spread exceeds F vox
 *                       (default 4.75 = wrap_spacing/2; 0 disables the cut)
 *     --iterations N    SLIM local/global iterations (default 40)
 *     --pad F           component packing gutter, vox (default 20)
 *     --verbose
 *   swirl_slim --selftest
 */
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../flatten/slim_refine.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int rc = SlimRefine_selftest() + MeshBin_selftest();
        fprintf(stderr, "=== swirl_slim selftest (%d failure%s) ===\n",
                rc, rc == 1 ? "" : "s");
        return rc ? 1 : 0;
    }
    if (argc < 3) {
        fprintf(stderr, "usage: swirl_slim <in.vmesh> <out.vmesh> "
                "[--umbilicus Y X] [--cut-gate F] [--iterations N] [--pad F] "
                "[--verbose]\n");
        return 1;
    }
    const char *in_path = argv[1], *out_path = argv[2];
    double uy = 3405.0, ux = 2878.0, cut_gate = 4.75, pad = 20.0;
    int iterations = 40, verbose = 0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--umbilicus") && i + 2 < argc) {
            uy = atof(argv[++i]); ux = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--cut-gate") && i + 1 < argc) {
            cut_gate = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--iterations") && i + 1 < argc) {
            iterations = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--pad") && i + 1 < argc) {
            pad = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--verbose")) {
            verbose = 1;
        } else {
            fprintf(stderr, "swirl_slim: unknown arg %s\n", argv[i]);
            return 1;
        }
    }

    char resolved[4096];
    MeshBinData mesh;
    if (MeshBin_companion_path(in_path, resolved, sizeof resolved) != 0 ||
        MeshBin_read_malloc(resolved, &mesh) != 0) {
        fprintf(stderr, "swirl_slim: cannot read %s\n", in_path);
        return 1;
    }
    if (mesh.uv == NULL) {
        fprintf(stderr, "swirl_slim: input has no UV (SLIM needs a start)\n");
        MeshBin_dispose(&mesh);
        return 1;
    }
    size_t nv = mesh.nv, nf = mesh.nf;

    /* per-vertex radius from the umbilicus (axis along Z; verts are (z,y,x)) */
    double *rad = (double *)malloc(nv * sizeof *rad);
    int32_t *kf = (int32_t *)malloc(nf * 3 * sizeof *kf);
    if (rad == NULL || kf == NULL) {
        free(rad); free(kf); MeshBin_dispose(&mesh); return 1;
    }
    for (size_t v = 0; v < nv; v++) {
        double dy = (double)mesh.verts[v * 3 + 1] - uy;
        double dx = (double)mesh.verts[v * 3 + 2] - ux;
        rad[v] = sqrt(dy * dy + dx * dx);
    }
    /* keep faces whose three vertex radii span at most cut_gate */
    size_t nkf = 0, ncut = 0;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = mesh.faces[f * 3 + 0];
        int32_t b = mesh.faces[f * 3 + 1];
        int32_t c = mesh.faces[f * 3 + 2];
        double r0 = rad[a], r1 = rad[b], r2 = rad[c];
        double lo = r0 < r1 ? (r0 < r2 ? r0 : r2) : (r1 < r2 ? r1 : r2);
        double hi = r0 > r1 ? (r0 > r2 ? r0 : r2) : (r1 > r2 ? r1 : r2);
        if (cut_gate > 0.0 && hi - lo > cut_gate) { ncut++; continue; }
        kf[nkf * 3 + 0] = a; kf[nkf * 3 + 1] = b; kf[nkf * 3 + 2] = c; nkf++;
    }
    free(rad);
    fprintf(stderr, "[swirl_slim] nv=%zu nf=%zu; cut %zu bridges (gate %.2f) -> "
            "%zu faces\n", nv, nf, ncut, cut_gate, nkf);
    if (nkf < 1) { free(kf); MeshBin_dispose(&mesh); return 1; }

    double *uv_in = (double *)malloc(nv * 2 * sizeof *uv_in);
    double *uv_out = (double *)malloc(nv * 2 * sizeof *uv_out);
    if (uv_in == NULL || uv_out == NULL) {
        free(kf); free(uv_in); free(uv_out); MeshBin_dispose(&mesh); return 1;
    }
    for (size_t i = 0; i < nv * 2; i++) uv_in[i] = (double)mesh.uv[i];

    SlimRefineOpts opts;
    SlimRefine_defaults(&opts);
    opts.iterations = iterations;
    opts.strip_pack = 1;      /* pack components largest-first -> kibble tail */
    opts.padding = pad;
    opts.guard_boundary = 1;  /* reject self-intersecting chart boundaries */
    opts.verbose = verbose;

    Arena_T arena = Arena_new();
    SlimRefineStats st;
    memset(&st, 0, sizeof st);
    int rc = (arena == NULL) ? -1
           : SlimRefine_run_double(arena, mesh.verts, nv, kf, nkf, uv_in,
                                   &opts, uv_out, &st);
    if (rc != 0) {
        fprintf(stderr, "swirl_slim: SlimRefine_run_double failed (rc=%d)\n", rc);
        Arena_dispose(&arena);
        free(kf); free(uv_in); free(uv_out); MeshBin_dispose(&mesh);
        return 1;
    }
    fprintf(stderr,
            "[swirl_slim] components=%zu  flips %zu->%zu  energy %.5g->%.5g  "
            "min_det %.4g->%.4g  qc_max %.3f->%.3f  atlas %.0fx%.0f  %.1fs  %s\n",
            st.components, st.flips_before, st.flips_after,
            st.energy_before, st.energy_after, st.min_det_before, st.min_det_after,
            st.qc_max_before, st.qc_max_after, st.atlas_width, st.atlas_height,
            st.total_seconds,
            st.converged ? "converged" : (st.reverted ? "REVERTED" : "stopped"));

    float *uvf = (float *)malloc(nv * 2 * sizeof *uvf);
    if (uvf == NULL) {
        Arena_dispose(&arena);
        free(kf); free(uv_in); free(uv_out); MeshBin_dispose(&mesh); return 1;
    }
    for (size_t i = 0; i < nv * 2; i++) uvf[i] = (float)uv_out[i];
    int wrc = MeshBin_write(out_path, mesh.verts, nv, kf, nkf, uvf);
    if (wrc != 0)
        fprintf(stderr, "swirl_slim: cannot write %s (rc=%d)\n", out_path, wrc);
    else
        fprintf(stderr, "[swirl_slim] wrote %s (%zu verts, %zu faces)\n",
                out_path, nv, nkf);

    free(uvf); Arena_dispose(&arena);
    free(kf); free(uv_in); free(uv_out); MeshBin_dispose(&mesh);
    return wrc ? 1 : 0;
}
