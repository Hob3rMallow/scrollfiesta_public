/* ============================================================================
 * ribbon_relax_probe.c -- fast step4-only harness for tuning the ribbon-relax
 * optimizer WITHOUT re-running the ~10-min winding/ribbon/overlap/snap stages.
 *
 * Loads a snapped UV OBJ (v + vt + f, as written by scroll_ribbon's
 * <id>_snapped.obj) and the step3 rawtex TIF, recovers the fiber field, runs
 * RibbonRelax with CLI-tunable parameters, writes <out>.obj and prints the
 * before/after stretch, fiber axis-error, and vertex displacement. The qc-reject
 * gate drops the overlap-relocated/degenerate faces, so no face_skip mask is
 * needed here. Bake the result with obj_bake_raw to view it.
 * ==========================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/obj_io.h"
#include "../flatten/fiber_field.h"
#include "../flatten/ribbon_relax.h"

/* Parse "v z y x", "vt u v", "f a/a b/b c/c" (or "f a b c"). vt is 1:1 with v
 * (scroll_ribbon writes a/a). malloc's outputs; caller frees. */
static int load_uv_obj(const char *path, float **pv, float **puv, int32_t **pf,
                       size_t *pnv, size_t *pnf)
{
    FILE *fp = fopen(path, "r");
    char line[1024];
    size_t nv = 0, nvt = 0, nf = 0, vi = 0, vti = 0, fi = 0;
    float *V, *UV; int32_t *F;
    if (fp == NULL) return -1;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'v' && line[1] == ' ') nv++;
        else if (line[0] == 'v' && line[1] == 't') nvt++;
        else if (line[0] == 'f' && line[1] == ' ') nf++;
    }
    if (nv == 0 || nf == 0) { fclose(fp); return -1; }
    V = (float *)malloc(nv * 3 * sizeof *V);
    UV = (float *)malloc((nvt ? nvt : nv) * 2 * sizeof *UV);
    F = (int32_t *)malloc(nf * 3 * sizeof *F);
    if (V == NULL || UV == NULL || F == NULL) { free(V); free(UV); free(F); fclose(fp); return -1; }
    rewind(fp);
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'v' && line[1] == ' ') {
            double z, y, x;
            if (sscanf(line + 2, "%lf %lf %lf", &z, &y, &x) == 3 && vi < nv) {
                V[vi*3] = (float)z; V[vi*3+1] = (float)y; V[vi*3+2] = (float)x; vi++;
            }
        } else if (line[0] == 'v' && line[1] == 't') {
            double u, v;
            if (sscanf(line + 3, "%lf %lf", &u, &v) == 2 && vti < (nvt ? nvt : nv)) {
                UV[vti*2] = (float)u; UV[vti*2+1] = (float)v; vti++;
            }
        } else if (line[0] == 'f' && line[1] == ' ') {
            int a, b, c;
            if ((sscanf(line + 2, "%d/%*d %d/%*d %d/%*d", &a, &b, &c) == 3 ||
                 sscanf(line + 2, "%d %d %d", &a, &b, &c) == 3) && fi < nf) {
                F[fi*3] = a-1; F[fi*3+1] = b-1; F[fi*3+2] = c-1; fi++;
            }
        }
    }
    fclose(fp);
    if (vti != vi) { fprintf(stderr, "  WARN: %zu vt != %zu v (using identity for missing)\n", vti, vi); }
    *pv = V; *puv = UV; *pf = F; *pnv = vi; *pnf = fi;
    return 0;
}

int main(int argc, char **argv)
{
    const char *obj = NULL, *tex = NULL, *out = NULL;
    double grad_sigma = 1.5, tensor_sigma = 6.0;
    RibbonRelaxOpts o;
    RibbonRelaxStats st;
    float *V = NULL, *UV = NULL, *OUT = NULL; int32_t *F = NULL;
    size_t nv = 0, nf = 0; int i;
    uint8_t *img = NULL; int td = 0, tw = 0, th = 0;
    Arena_T arena;
    FiberField fib;

    RibbonRelax_defaults(&o); o.verbose = 1;
    for (i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "--lambda") == 0 && i+1 < argc)        o.lambda_align = atof(argv[++i]);
        else if (strcmp(argv[i], "--sweeps") == 0 && i+1 < argc)        o.sweeps = atoi(argv[++i]);
        else if (strcmp(argv[i], "--coh-gate") == 0 && i+1 < argc)      o.coh_gate = atof(argv[++i]);
        else if (strcmp(argv[i], "--max-disp") == 0 && i+1 < argc)      o.max_disp = atof(argv[++i]);
        else if (strcmp(argv[i], "--qc-reject") == 0 && i+1 < argc)     o.qc_reject = atof(argv[++i]);
        else if (strcmp(argv[i], "--free-boundary") == 0)               o.fix_boundary = 0;
        else if (strcmp(argv[i], "--grad-sigma") == 0 && i+1 < argc)    grad_sigma = atof(argv[++i]);
        else if (strcmp(argv[i], "--tensor-sigma") == 0 && i+1 < argc)  tensor_sigma = atof(argv[++i]);
        else if (obj == NULL) obj = argv[i];
        else if (tex == NULL) tex = argv[i];
        else if (out == NULL) out = argv[i];
    }
    if (obj == NULL || tex == NULL || out == NULL) {
        fprintf(stderr, "usage: ribbon_relax_probe <snapped_uv.obj> <step3_rawtex.tif> <out.obj>\n"
                "   [--lambda f] [--sweeps n] [--coh-gate f] [--free-boundary]\n"
                "   [--max-disp f] [--qc-reject f] [--grad-sigma f] [--tensor-sigma f]\n");
        return 2;
    }
    setvbuf(stderr, NULL, _IONBF, 0);

    if (load_uv_obj(obj, &V, &UV, &F, &nv, &nf) != 0) { fprintf(stderr, "load fail: %s\n", obj); return 1; }
    fprintf(stderr, "[relax_probe] %s: %zu verts, %zu faces\n", obj, nv, nf);

    arena = Arena_new();
    if (TiffIO_load(arena, tex, &img, &td, &th, &tw) != 0) { fprintf(stderr, "tex load fail: %s\n", tex); return 1; }
    if (FiberField_compute(arena, img, tw, th, grad_sigma, tensor_sigma, 15.0, 4, &fib) != 0) {
        fprintf(stderr, "fiber compute fail\n"); return 1;
    }
    fprintf(stderr, "[relax_probe] fiber: %dx%d mean_coh=%.3f axis_err=%.2f\n", tw, th, fib.mean_coh, fib.mean_axis_err_deg);

    OUT = (float *)malloc(nv * 2 * sizeof *OUT);
    fprintf(stderr, "[relax_probe] opts: lambda=%.2f sweeps=%d coh_gate=%.2f max_disp=%.1f qc_reject=%.0f fix_boundary=%d\n",
            o.lambda_align, o.sweeps, o.coh_gate, o.max_disp, o.qc_reject, o.fix_boundary);
    if (RibbonRelax_run(arena, V, nv, F, nf, UV, NULL, &fib, 1.0, 1.0, &o, OUT, &st) != 0) {
        fprintf(stderr, "relax fail\n"); return 1;
    }
    fprintf(stderr,
        "[relax_probe] RESULT  stretch %.4f -> %.4f (max %.2f)  axis-err %.3f -> %.3f deg\n"
        "              disp mean=%.3f max=%.3f vox  interior=%zu moved=%zu reject=%zu fiber=%zu flips %d->%d%s\n",
        st.stretch_mean_before, st.stretch_mean_after, st.stretch_max_after,
        st.axis_err_before, st.axis_err_after,
        st.mean_disp, st.max_disp, st.n_interior, st.n_moved, st.n_reject, st.n_fiber_faces,
        st.flips_before, st.flips_after, st.reverted ? "  (REVERTED)" : "");

    if (ObjIO_write_uv(out, V, nv, F, nf, OUT) != 0) fprintf(stderr, "  WARN: write %s failed\n", out);
    else fprintf(stderr, "[relax_probe] wrote %s\n", out);

    free(V); free(UV); free(F); free(OUT);
    Arena_dispose(&arena);
    return 0;
}
