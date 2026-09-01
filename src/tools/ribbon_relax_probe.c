/* ============================================================================
 * ribbon_relax_probe.c -- fast UV-only Sander/symmetric-Dirichlet harness.
 *
 * Historical mode loads a snapped UV OBJ plus its raw-texture TIF and tunes
 * the optional fiber-alignment term.  The VMESH mode is the topology-preserving
 * post-process used by the probabilistic-winding pipeline: it can run pure
 * geometric isometry, retain the accepted winding/layout as a soft guide, and
 * write an identically indexed VMESH plus a machine-readable report.
 * ==========================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/obj_io.h"
#include "../common/tiff_io.h"
/* Avoid ves_platform.h's annotation macros here: this CLI retains the
 * conventional local name `out`.  The implementation is already linked. */
int ves_ensure_parent_dir(const char *filepath);
#include "../flatten/fiber_field.h"
#include "../flatten/ribbon_relax.h"

/* Parse "v z y x", "vt u v", "f a/a b/b c/c" (or "f a b c"). */
static int load_uv_obj(const char *path, float **pv, float **puv, int32_t **pf,
                       size_t *pnv, size_t *pnf)
{
    FILE *fp = fopen(path, "r");
    char line[1024];
    size_t nv = 0, nvt = 0, nf = 0, vi = 0, vti = 0, fi = 0;
    float *V, *UV;
    int32_t *F;
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
    if (V == NULL || UV == NULL || F == NULL) {
        free(V); free(UV); free(F); fclose(fp); return -1;
    }
    rewind(fp);
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == 'v' && line[1] == ' ') {
            double z, y, x;
            if (sscanf(line + 2, "%lf %lf %lf", &z, &y, &x) == 3 && vi < nv) {
                V[vi*3] = (float)z;
                V[vi*3+1] = (float)y;
                V[vi*3+2] = (float)x;
                vi++;
            }
        } else if (line[0] == 'v' && line[1] == 't') {
            double u, v;
            if (sscanf(line + 3, "%lf %lf", &u, &v) == 2 &&
                vti < (nvt ? nvt : nv)) {
                UV[vti*2] = (float)u;
                UV[vti*2+1] = (float)v;
                vti++;
            }
        } else if (line[0] == 'f' && line[1] == ' ') {
            int a, b, c;
            if ((sscanf(line + 2, "%d/%*d %d/%*d %d/%*d", &a, &b, &c) == 3 ||
                 sscanf(line + 2, "%d %d %d", &a, &b, &c) == 3) && fi < nf) {
                F[fi*3] = a-1;
                F[fi*3+1] = b-1;
                F[fi*3+2] = c-1;
                fi++;
            }
        }
    }
    fclose(fp);
    if (vti != vi)
        fprintf(stderr, "  WARN: %zu vt != %zu v\n", vti, vi);
    *pv = V; *puv = UV; *pf = F; *pnv = vi; *pnf = fi;
    return 0;
}

static void json_string(FILE *fp, const char *text)
{
    const unsigned char *p = (const unsigned char *)text;
    fputc('"', fp);
    for (; *p; p++) {
        if (*p == '"' || *p == '\\') fprintf(fp, "\\%c", *p);
        else if (*p == '\n') fputs("\\n", fp);
        else if (*p == '\r') fputs("\\r", fp);
        else if (*p == '\t') fputs("\\t", fp);
        else if (*p < 0x20) fprintf(fp, "\\u%04x", (unsigned)*p);
        else fputc(*p, fp);
    }
    fputc('"', fp);
}

static int write_relax_report(const char *out_path, const char *input_path,
                              const RibbonRelaxOpts *o,
                              const RibbonRelaxStats *st,
                              size_t nv, size_t nf)
{
    char path[4096];
    int n = snprintf(path, sizeof path, "%s.relax.json", out_path);
    FILE *fp;
    if (n < 0 || (size_t)n >= sizeof path) return -1;
    fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    fputs("{\n  \"schema\": \"vmesh-sander-relax-v1\",\n  \"input\": ", fp);
    json_string(fp, input_path);
    fputs(",\n  \"output\": ", fp);
    json_string(fp, out_path);
    fprintf(fp,
        ",\n  \"vertices\": %zu, \"faces\": %zu,\n"
        "  \"policy\": {\"orientation_preserving\": true, "
        "\"topology_preserving\": true, \"lambda_align\": %.9g, "
        "\"guide_weight\": %.9g, \"max_displacement\": %.9g, "
        "\"max_stretch_growth\": %.9g, \"stretch_guard_floor\": %.9g, "
        "\"sweeps\": %d},\n"
        "  \"result\": {\"iterations\": %d, \"converged\": %s, "
        "\"reverted\": %s, \"stretch_mean_before\": %.9g, "
        "\"stretch_mean_after\": %.9g, \"stretch_max_before\": %.9g, "
        "\"stretch_max_after\": %.9g, \"energy_before\": %.17g, "
        "\"energy_after\": %.17g, \"mean_displacement\": %.9g, "
        "\"max_displacement\": %.9g, \"flips_before\": %d, "
        "\"flips_after\": %d, \"moved_updates\": %zu, "
        "\"rejected_faces\": %zu, \"stretch_guard_rejections\": %zu}\n}\n",
        nv, nf, o->lambda_align, o->guide_weight, o->max_disp,
        o->max_stretch_growth, o->stretch_guard_floor, o->sweeps,
        st->iterations, st->converged ? "true" : "false",
        st->reverted ? "true" : "false",
        st->stretch_mean_before, st->stretch_mean_after,
        st->stretch_max_before, st->stretch_max_after,
        st->energy_before, st->energy_after, st->mean_disp, st->max_disp,
        st->flips_before, st->flips_after, st->n_moved, st->n_reject,
        st->n_stretch_guard_reject);
    if (fclose(fp) != 0) return -1;
    fprintf(stderr, "[relax_probe] wrote %s\n", path);
    return 0;
}

static void usage(void)
{
    fprintf(stderr,
        "usage:\n"
        "  ribbon_relax_probe --vmesh <in.vmesh> <out.vmesh> "
        "--isometry-only [options]\n"
        "  ribbon_relax_probe <snapped_uv.obj> <rawtex.tif> <out.obj> "
        "[options]\n"
        "options:\n"
        "  --lambda F --guide-weight F --sweeps N --movement-tolerance F\n"
        "  --max-disp F --qc-reject F --line-search N\n"
        "  --max-stretch-growth F --stretch-guard-floor F\n"
        "  --free-boundary | --fix-boundary | --global-injective\n"
        "  --texture rawtex.tif --grad-sigma F --tensor-sigma F\n");
}

int main(int argc, char **argv)
{
    const char *obj = NULL, *tex = NULL, *out = NULL;
    const char *vmesh_in = NULL, *vmesh_out = NULL;
    double grad_sigma = 1.5, tensor_sigma = 6.0;
    int isometry_only = 0;
    RibbonRelaxOpts o;
    RibbonRelaxStats st;
    float *V = NULL, *UV = NULL, *OUT = NULL;
    int32_t *F = NULL;
    size_t nv = 0, nf = 0;
    uint8_t *img = NULL;
    int td = 0, tw = 0, th = 0;
    Arena_T arena = NULL;
    FiberField fib;
    int i;

    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return RibbonRelax_selftest();

    RibbonRelax_defaults(&o);
    o.verbose = 1;
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--vmesh") == 0 && i + 2 < argc) {
            vmesh_in = argv[++i];
            vmesh_out = argv[++i];
        } else if (strcmp(argv[i], "--texture") == 0 && i + 1 < argc) {
            tex = argv[++i];
        } else if (strcmp(argv[i], "--isometry-only") == 0) {
            isometry_only = 1;
        } else if (strcmp(argv[i], "--guide-weight") == 0 && i + 1 < argc) {
            o.guide_weight = atof(argv[++i]);
        } else if (strcmp(argv[i], "--movement-tolerance") == 0 && i + 1 < argc) {
            o.movement_tolerance = atof(argv[++i]);
        } else if (strcmp(argv[i], "--line-search") == 0 && i + 1 < argc) {
            o.line_search = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--lambda") == 0 && i + 1 < argc) {
            o.lambda_align = atof(argv[++i]);
        } else if (strcmp(argv[i], "--sweeps") == 0 && i + 1 < argc) {
            o.sweeps = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--coh-gate") == 0 && i + 1 < argc) {
            o.coh_gate = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-disp") == 0 && i + 1 < argc) {
            o.max_disp = atof(argv[++i]);
        } else if (strcmp(argv[i], "--qc-reject") == 0 && i + 1 < argc) {
            o.qc_reject = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-stretch-growth") == 0 && i + 1 < argc) {
            o.max_stretch_growth = atof(argv[++i]);
        } else if (strcmp(argv[i], "--stretch-guard-floor") == 0 && i + 1 < argc) {
            o.stretch_guard_floor = atof(argv[++i]);
        } else if (strcmp(argv[i], "--free-boundary") == 0) {
            o.fix_boundary = 0;
        } else if (strcmp(argv[i], "--fix-boundary") == 0) {
            o.fix_boundary = 1;
        } else if (strcmp(argv[i], "--global-injective") == 0) {
            o.global_injective = 1;
        } else if (strcmp(argv[i], "--grad-sigma") == 0 && i + 1 < argc) {
            grad_sigma = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tensor-sigma") == 0 && i + 1 < argc) {
            tensor_sigma = atof(argv[++i]);
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "unknown option: %s\n", argv[i]);
            usage();
            return 2;
        } else if (obj == NULL) obj = argv[i];
        else if (tex == NULL) tex = argv[i];
        else if (out == NULL) out = argv[i];
        else {
            fprintf(stderr, "unexpected positional argument: %s\n", argv[i]);
            return 2;
        }
    }
    if (isometry_only) o.lambda_align = 0.0;
    if (o.sweeps < 0 || o.line_search < 1 || o.guide_weight < 0.0 ||
        o.max_disp < 0.0 || o.movement_tolerance < 0.0 ||
        (o.max_stretch_growth > 0.0 && o.max_stretch_growth < 1.0) ||
        o.stretch_guard_floor < 1.0) {
        fprintf(stderr, "invalid relaxation option\n");
        return 2;
    }
    setvbuf(stderr, NULL, _IONBF, 0);

    if (vmesh_in != NULL || vmesh_out != NULL) {
        MeshBinData mesh;
        double *din = NULL, *dout = NULL;
        float *fout = NULL;
        const FiberField *fibp = NULL;
        char resolved[4096];
        int rc = 1;
        memset(&mesh, 0, sizeof mesh);
        if (vmesh_in == NULL || vmesh_out == NULL) {
            fprintf(stderr, "--vmesh requires input and output paths\n");
            return 2;
        }
        if (MeshBin_companion_path(vmesh_in, resolved, sizeof resolved) != 0 ||
            MeshBin_read_malloc(resolved, &mesh) != 0 || mesh.uv == NULL) {
            fprintf(stderr, "vmesh load fail: %s\n", vmesh_in);
            return 1;
        }
        fprintf(stderr, "[relax_probe] %s: %zu verts, %zu faces\n",
                resolved, mesh.nv, mesh.nf);
        arena = Arena_new();
        if (arena == NULL) goto vmesh_done;
        if (o.lambda_align > 0.0) {
            if (tex == NULL) {
                fprintf(stderr,
                        "fiber alignment requires --texture; use "
                        "--isometry-only for geometric relaxation\n");
                goto vmesh_done;
            }
            if (TiffIO_load(arena, tex, &img, &td, &th, &tw) != 0 ||
                FiberField_compute(arena, img, tw, th, grad_sigma,
                                   tensor_sigma, 15.0, 4, &fib) != 0) {
                fprintf(stderr, "texture/fiber load fail: %s\n", tex);
                goto vmesh_done;
            }
            fibp = &fib;
        }
        din = (double *)malloc(mesh.nv * 2 * sizeof(*din));
        dout = (double *)malloc(mesh.nv * 2 * sizeof(*dout));
        fout = (float *)malloc(mesh.nv * 2 * sizeof(*fout));
        if (din == NULL || dout == NULL || fout == NULL) goto vmesh_done;
        for (size_t k = 0; k < mesh.nv * 2; k++) din[k] = mesh.uv[k];
        o.reference_uv = din;
        if (o.guide_weight > 0.0) o.guide_uv = din;
        fprintf(stderr,
                "[relax_probe] vmesh opts: lambda=%.3g guide=%.3g "
                "sweeps=%d tolerance=%.3g max_disp=%.3g qc_reject=%.3g "
                "stretch_guard=(%.3g, %.3g) fix_boundary=%d\n",
                o.lambda_align, o.guide_weight, o.sweeps,
                o.movement_tolerance, o.max_disp, o.qc_reject,
                o.max_stretch_growth, o.stretch_guard_floor,
                o.fix_boundary);
        if (RibbonRelax_run_double(arena, mesh.verts, mesh.nv, mesh.faces,
                                   mesh.nf, din, NULL, fibp, 1.0, 1.0,
                                   &o, dout, &st) != 0) {
            fprintf(stderr, "relax fail\n");
            goto vmesh_done;
        }
        fprintf(stderr,
            "[relax_probe] RESULT stretch %.6f -> %.6f (max %.3g -> %.3g) "
            "energy %.9g -> %.9g\n"
            "              disp mean=%.4f max=%.4f vox interior=%zu "
            "moved=%zu reject=%zu stretch_guard_reject=%zu flips %d->%d%s\n",
            st.stretch_mean_before, st.stretch_mean_after,
            st.stretch_max_before, st.stretch_max_after,
            st.energy_before, st.energy_after, st.mean_disp, st.max_disp,
            st.n_interior, st.n_moved, st.n_reject,
            st.n_stretch_guard_reject,
            st.flips_before, st.flips_after,
            st.reverted ? " (REVERTED)" : "");
        for (size_t k = 0; k < mesh.nv * 2; k++) fout[k] = (float)dout[k];
        if (ves_ensure_parent_dir(vmesh_out) != 0 ||
            MeshBin_write(vmesh_out, mesh.verts, mesh.nv, mesh.faces,
                          mesh.nf, fout) != 0) {
            fprintf(stderr, "vmesh write fail: %s\n", vmesh_out);
            goto vmesh_done;
        }
        fprintf(stderr, "[relax_probe] wrote %s\n", vmesh_out);
        if (write_relax_report(vmesh_out, resolved, &o, &st,
                               mesh.nv, mesh.nf) != 0)
            fprintf(stderr, "  WARN: could not write relaxation report\n");
        rc = 0;
vmesh_done:
        free(din); free(dout); free(fout);
        if (arena != NULL) Arena_dispose(&arena);
        MeshBin_dispose(&mesh);
        return rc;
    }

    if (obj == NULL || tex == NULL || out == NULL) {
        usage();
        return 2;
    }
    if (load_uv_obj(obj, &V, &UV, &F, &nv, &nf) != 0) {
        fprintf(stderr, "load fail: %s\n", obj);
        return 1;
    }
    fprintf(stderr, "[relax_probe] %s: %zu verts, %zu faces\n", obj, nv, nf);
    arena = Arena_new();
    if (arena == NULL ||
        TiffIO_load(arena, tex, &img, &td, &th, &tw) != 0 ||
        FiberField_compute(arena, img, tw, th, grad_sigma, tensor_sigma,
                           15.0, 4, &fib) != 0) {
        fprintf(stderr, "texture/fiber load fail: %s\n", tex);
        free(V); free(UV); free(F);
        if (arena != NULL) Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr,
            "[relax_probe] fiber: %dx%d mean_coh=%.3f axis_err=%.2f\n",
            tw, th, fib.mean_coh, fib.mean_axis_err_deg);
    OUT = (float *)malloc(nv * 2 * sizeof *OUT);
    if (OUT == NULL ||
        RibbonRelax_run(arena, V, nv, F, nf, UV, NULL, &fib, 1.0, 1.0,
                        &o, OUT, &st) != 0) {
        fprintf(stderr, "relax fail\n");
        free(V); free(UV); free(F); free(OUT); Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr,
        "[relax_probe] RESULT stretch %.4f -> %.4f (max %.2f) "
        "axis-err %.3f -> %.3f deg\n"
        "              disp mean=%.3f max=%.3f vox interior=%zu "
        "moved=%zu reject=%zu fiber=%zu flips %d->%d%s\n",
        st.stretch_mean_before, st.stretch_mean_after, st.stretch_max_after,
        st.axis_err_before, st.axis_err_after, st.mean_disp, st.max_disp,
        st.n_interior, st.n_moved, st.n_reject, st.n_fiber_faces,
        st.flips_before, st.flips_after, st.reverted ? " (REVERTED)" : "");
    if (ObjIO_write_uv(out, V, nv, F, nf, OUT) != 0)
        fprintf(stderr, "  WARN: write %s failed\n", out);
    else
        fprintf(stderr, "[relax_probe] wrote %s\n", out);
    free(V); free(UV); free(F); free(OUT);
    Arena_dispose(&arena);
    return 0;
}
