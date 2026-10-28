/* Native whole-sheet rigid alignment; physical fibre arms are measured by
 * ink_closeup --measure-fibers and remain in original face order. */
#include "../flatten/fiber_rigid.h"
#include "../common/mesh_bin.h"
#include "../common/ves_platform.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int rigid_cli_selftest(void);

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return FiberRigid_selftest() || rigid_cli_selftest();
    if (argc != 5) {
        fprintf(stderr, "usage: vmesh_fiber_rigid input.vmesh arms.f64 output.vmesh report.json\n");
        return 2;
    }
    Arena_T arena = Arena_new();
    MeshBinData mesh = {0};
    FiberRigidStats stats = {0};
    FILE *fp = NULL;
    int rc = 1;
    double started = ves_clock_sec();
    /* Never overwrite either an input or an earlier result. */
    for (int k = 3; k < 5; k++) {
        fp = fopen(argv[k], "rb");
        if (fp) { fclose(fp); fp = NULL; fprintf(stderr, "ERROR: output exists: %s\n", argv[k]); goto done; }
    }
    if (MeshBin_read_precise_arena(arena, argv[1], &mesh) || (!mesh.uv && !mesh.uv64) ||
            mesh.nf > SIZE_MAX/sizeof(SlimRefineDirection)) goto done;
    double *uv = mesh.uv64;
    if (!uv) {
        uv = ARENA_ALLOC(arena, 2*mesh.nv*sizeof *uv);
        for (size_t k = 0; k < 2*mesh.nv; k++) uv[k] = mesh.uv[k];
    }
    SlimRefineDirection *arms = ARENA_ALLOC(arena, mesh.nf*sizeof *arms);
    fp = fopen(argv[2], "rb");
    if (!fp || fread(arms, sizeof *arms, mesh.nf, fp) != mesh.nf || fgetc(fp) != EOF || ferror(fp)) goto done;
    fclose(fp); fp = NULL;
    if (FiberRigid_fit(arena, mesh.verts, mesh.nv, mesh.faces, mesh.nf, uv, arms, &stats)) goto done;
    double *out = ARENA_ALLOC(arena, 2*mesh.nv*sizeof *out);
    if (FiberRigid_apply(uv, mesh.nv, &stats, out)) goto done;
    double magnitude = 1, max_vector = 0, max_length = 0;
    double co = cos(stats.radians), si = sin(stats.radians);
    for (size_t k = 0; k < 2*mesh.nv; k++) {
        if (!isfinite(out[k])) goto done;
        magnitude = fmax(magnitude, fmax(fabs(uv[k]), fabs(out[k])));
    }
    for (size_t f = 0; f < mesh.nf; f++) for (size_t k = 0; k < 3; k++) {
        size_t a = (size_t)mesh.faces[3*f+k], b = (size_t)mesh.faces[3*f+(k+1)%3];
        double x = uv[2*b]-uv[2*a], y = uv[2*b+1]-uv[2*a+1];
        double u = out[2*b]-out[2*a], v = out[2*b+1]-out[2*a+1];
        max_vector = fmax(max_vector, hypot(u-(co*x-si*y), v-(si*x+co*y)));
        max_length = fmax(max_length, fabs(hypot(u, v)-hypot(x, y)));
    }
    double roundoff_bound = 64*DBL_EPSILON*magnitude;
    if (max_vector > roundoff_bound || max_length > roundoff_bound) goto done;
    ves_ensure_parent_dir(argv[3]);
    if (MeshBin_write_uv64(argv[3], mesh.verts, mesh.nv, mesh.faces, mesh.nf, out)) goto done;
    ves_ensure_parent_dir(argv[4]);
    fp = fopen(argv[4], "wb");
    if (!fp) goto done;
    fprintf(fp, "{\n \"schema\":\"whole-sheet-fiber-rigid-v1\",\"passed\":true,\n"
            " \"vertices\":%zu,\"faces\":%zu,\"supported_faces\":%zu,\n"
            " \"area\":%.17g,\"supported_area\":%.17g,\"supported_area_fraction\":%.17g,\n"
            " \"evidence_mass\":%.17g,\"rotation_radians\":%.17g,\"rotation_degrees\":%.17g,\n"
            " \"center\":[%.17g,%.17g],\"rms_before_degrees\":%.17g,\"rms_after_degrees\":%.17g,\n"
            " \"signed_mean_before_degrees\":%.17g,\"signed_mean_after_degrees\":%.17g,\n"
            " \"max_edge_vector_error\":%.17g,\"max_edge_length_error\":%.17g,\"roundoff_bound\":%.17g,\n"
            " \"one_common_rotation\":true,\"xyz_faces_unchanged\":true,\"iterative_solves\":0,\n"
            " \"seconds\":%.6f,\"production_qualified\":false\n}\n",
            mesh.nv, mesh.nf, stats.supported_faces, stats.area, stats.supported_area, stats.supported_area/stats.area,
            stats.evidence_mass, stats.radians, stats.radians*57.29577951308232,
            stats.center[0], stats.center[1], stats.rms_before_degrees, stats.rms_after_degrees,
            stats.signed_mean_before_degrees, stats.signed_mean_after_degrees,
            max_vector, max_length, roundoff_bound, ves_clock_sec()-started);
    rc = fclose(fp) != 0; fp = NULL;
    fprintf(stderr, "[fiber rigid] turn %.9g degrees; RMS %.9g -> %.9g; %.3fs\n",
            stats.radians*57.29577951308232, stats.rms_before_degrees, stats.rms_after_degrees, ves_clock_sec()-started);
done:
    if (fp) fclose(fp);
    Arena_dispose(&arena);
    if (rc) fprintf(stderr, "ERROR: rigid alignment failed; no qualification\n");
    return rc;
}

static int rigid_cli_selftest(void)
{
    char dir[512], input[600], arms_path[600], output[600], report[600];
    if (ves_selftest_dir(dir, sizeof dir, "vmesh_fiber_rigid")) return 1;
    snprintf(input, sizeof input, "%s/input.vmesh", dir);
    snprintf(arms_path, sizeof arms_path, "%s/arms.f64", dir);
    snprintf(output, sizeof output, "%s/output.vmesh", dir);
    snprintf(report, sizeof report, "%s/report.json", dir);
    const float xyz[9] = {0,0,0, 2,0,0, 0,3,0};
    const int32_t faces[3] = {0,1,2};
    const double co = cos(.2), si = sin(.2);
    const double uv[6] = {100,200, 100+2*co,200+2*si, 100-3*si,200+3*co};
    const SlimRefineDirection arms = {{1,0,0}, {0,1,0}, 1};
    int failed = MeshBin_write_uv64(input, xyz, 3, faces, 1, uv) != 0;
    FILE *fp = fopen(arms_path, "wb");
    if (!fp) return 1;
    if (fwrite(&arms, sizeof arms, 1, fp) != 1) failed = 1;
    if (fclose(fp)) failed = 1;
    char *args[5] = {"vmesh_fiber_rigid", input, arms_path, output, report};
    if (!failed && main(5, args)) failed = 1;
    MeshBinData result = {0};
    if (!failed && MeshBin_read_precise_malloc(output, &result)) failed = 1;
    if (!failed && (result.nv != 3 || result.nf != 1 || !result.uv64 ||
            memcmp(result.verts, xyz, sizeof xyz) || memcmp(result.faces, faces, sizeof faces) ||
            fabs(result.uv64[3]-result.uv64[1]) > 1e-10 ||
            fabs(result.uv64[4]-result.uv64[0]) > 1e-10)) failed = 1;
    MeshBin_dispose(&result);
    if (!failed && main(5, args) == 0) failed = 1; /* immutable output guard */
    if (!failed) (void)ves_selftest_clean(dir);
    fprintf(stderr, "FiberRigid CLI roundtrip and overwrite guard: %s\n", failed ? "FAIL" : "PASS");
    return failed;
}
