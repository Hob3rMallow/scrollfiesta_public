/* qslim_obj.c -- decimate an OBJ mesh to a target face fraction with QEM
 * (Ozaki/Kanai LME-MDD, whole-mesh, quadric-error-driven; boundary verts are
 * pinned by QEM so open edges survive). OBJ in -> OBJ out. The face budget is
 * allocated across the whole mesh by error, so flat regions decimate hard and
 * detailed/curved regions keep faces -- better than forcing every component to
 * the same fraction.
 *
 * Usage:
 *   qslim_obj <in.obj> <out.obj> <keep_ratio>   keep_ratio in (0,1), e.g. 0.10
 *   qslim_obj --selftest
 *
 * Note: QEM collapses can locally mis-wind faces; the downstream seam-weld
 * (grid_weld) runs OrientMesh, so a decimated mesh fed straight into the weld
 * gets re-oriented there. (qslim_obj itself does not re-orient.)
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/qem.h"
#include "../remesh/isotropic_remesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* half-edge for the winding-consistency check in run_selftest */
typedef struct { int32_t lo, hi; int8_t fwd; } HE;
static int he_cmp(const void *pa, const void *pb)
{
    const HE *a = (const HE *)pa, *b = (const HE *)pb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}

/* Count "same_dir" manifold edges: an edge shared by exactly two faces whose
 * two half-edges traverse it the SAME way -> the faces are wound opposite (one
 * is back-face culled). 0 = globally consistent winding. Scratch on arena a. */
static long count_same_dir_he(Arena_T a, const int32_t *f, size_t nf)
{
    if (nf == 0) return 0;
    HE *he = (HE *)ARENA_ALLOC(a, (long)(nf * 3 * sizeof(HE)));
    size_t n = 0;
    for (size_t fi = 0; fi < nf; fi++) {
        for (int e = 0; e < 3; e++) {
            int32_t p = f[fi*3 + (size_t)e], q = f[fi*3 + (size_t)((e+1)%3)];
            he[n].lo = p < q ? p : q;
            he[n].hi = p < q ? q : p;
            he[n].fwd = (int8_t)(p < q ? 1 : 0);
            n++;
        }
    }
    qsort(he, n, sizeof(HE), he_cmp);
    long same_dir = 0;
    for (size_t i = 0; i + 1 < n; ) {
        size_t j = i + 1;
        while (j < n && he[j].lo == he[i].lo && he[j].hi == he[i].hi) j++;
        if (j - i == 2 && he[i].fwd == he[i+1].fwd) same_dir++;
        i = j;
    }
    return same_dir;
}

static int run_selftest(void)
{
    /* A flat W x W quad grid (2*W*W tris). QEM collapses the coplanar interior
     * freely (zero quadric error) while pinning the boundary, so decimation to
     * 1/4 faces is reachable and must reduce the mesh. */
    enum { W = 20 };
    const size_t side = W + 1;
    const size_t nv = side * side;
    const size_t nf = (size_t)(2 * W * W);
    Arena_T a = Arena_new();
    float *verts = (float *)ARENA_ALLOC(a, (long)(nv * 3 * sizeof(float)));
    int32_t *faces = (int32_t *)ARENA_ALLOC(a, (long)(nf * 3 * sizeof(int32_t)));
    for (size_t j = 0; j < side; j++) {
        for (size_t i = 0; i < side; i++) {
            size_t v = j * side + i;
            verts[v * 3 + 0] = (float)i;
            verts[v * 3 + 1] = (float)j;
            verts[v * 3 + 2] = 0.0f;
        }
    }
    size_t fi = 0;
    for (size_t j = 0; j < W; j++) {
        for (size_t i = 0; i < W; i++) {
            int32_t v00 = (int32_t)(j * side + i);
            int32_t v10 = (int32_t)(j * side + i + 1);
            int32_t v01 = (int32_t)((j + 1) * side + i);
            int32_t v11 = (int32_t)((j + 1) * side + i + 1);
            faces[fi * 3 + 0] = v00; faces[fi * 3 + 1] = v10; faces[fi * 3 + 2] = v11; fi++;
            faces[fi * 3 + 0] = v00; faces[fi * 3 + 1] = v11; faces[fi * 3 + 2] = v01; fi++;
        }
    }
    size_t target_nf = nf / 4;
    float *ov = NULL; int32_t *of = NULL; size_t onv = 0, onf = 0;
    int rc = QEM_simplify(a, verts, nv, faces, nf, target_nf, &ov, &onv, &of, &onf);
    int fails = 0;
    if (rc != 0)            { fprintf(stderr, "[selftest] QEM_simplify rc=%d -> FAIL\n", rc); fails++; }
    if (onf == 0)           { fprintf(stderr, "[selftest] empty output -> FAIL\n"); fails++; }
    if (onf > target_nf)    { fprintf(stderr, "[selftest] onf=%zu > target=%zu -> FAIL\n", onf, target_nf); fails++; }
    if (onf >= nf)          { fprintf(stderr, "[selftest] no reduction (onf=%zu nf=%zu) -> FAIL\n", onf, nf); fails++; }
    if (!fails)
        fprintf(stderr, "[selftest] grid %zu->%zu f, %zu->%zu v (target %zu) -> ok\n",
                nf, onf, nv, onv, target_nf);

    /* Case 2: a CONSISTENTLY-WOUND open cylinder (a wrapped surface). The
     * default QEM winding-repair flips faces against a single global reference
     * normal, which on a wrapped surface opposes ~half the faces -> whole
     * regions flip -> the "winding flips around each wrap" bug. With
     * skip_winding_repair the input's consistent winding must survive: count
     * "same_dir" manifold edges (both incident faces traverse the shared edge
     * the same way = a flip) and require ~0. */
    {
        enum { NU = 48, NH = 24 };   /* open cylinder: seam NOT wrapped */
        const size_t cnv = (size_t)NU * NH;
        const size_t cnf = (size_t)(NU - 1) * (NH - 1) * 2;
        float *cv = (float *)ARENA_ALLOC(a, (long)(cnv * 3 * sizeof(float)));
        int32_t *cf = (int32_t *)ARENA_ALLOC(a, (long)(cnf * 3 * sizeof(int32_t)));
        for (int jh = 0; jh < NH; jh++) {
            for (int iu = 0; iu < NU; iu++) {
                double ang = 2.0 * 3.14159265358979323846 * (double)iu / (double)NU;
                size_t v = (size_t)jh * NU + (size_t)iu;
                cv[v*3+0] = (float)(4.0 * (double)jh);      /* z = height */
                cv[v*3+1] = (float)(50.0 * sin(ang));       /* y */
                cv[v*3+2] = (float)(50.0 * cos(ang));       /* x */
            }
        }
        size_t ci = 0;
        for (int jh = 0; jh < NH - 1; jh++) {
            for (int iu = 0; iu < NU - 1; iu++) {           /* consistent CCW winding */
                int32_t a0 = (int32_t)((size_t)jh*NU + (size_t)iu);
                int32_t b0 = a0 + 1;
                int32_t c0 = (int32_t)((size_t)(jh+1)*NU + (size_t)iu);
                int32_t d0 = c0 + 1;
                cf[ci*3+0]=a0; cf[ci*3+1]=b0; cf[ci*3+2]=d0; ci++;
                cf[ci*3+0]=a0; cf[ci*3+1]=d0; cf[ci*3+2]=c0; ci++;
            }
        }
        QemOpts qo; memset(&qo, 0, sizeof qo);
        qo.free_boundary = 1;
        qo.skip_winding_repair = 1;
        float *cov = NULL; int32_t *cof = NULL; size_t conv = 0, conf = 0;
        int crc = QEM_simplify_opts(a, cv, cnv, cf, cnf, cnf / 4, &qo,
                                    &cov, &conv, &cof, &conf, NULL);
        /* same_dir manifold edges: both half-edges traverse the shared edge the
         * same way -> the two incident faces are wound opposite. 0 = consistent. */
        long same_dir = (crc == 0 && conf > 0) ? count_same_dir_he(a, cof, conf) : -1;
        int c2fail = 0;
        if (crc != 0)               { fprintf(stderr, "[selftest2] rc=%d -> FAIL\n", crc); c2fail++; }
        if (conf == 0 || conf >= cnf) { fprintf(stderr, "[selftest2] no reduction conf=%zu -> FAIL\n", conf); c2fail++; }
        if (same_dir > (long)(conf / 100 + 4)) {
            fprintf(stderr, "[selftest2] winding NOT preserved: %ld same_dir edges of %zu faces -> FAIL\n",
                    same_dir, conf); c2fail++;
        }
        if (!c2fail)
            fprintf(stderr, "[selftest2] wound cylinder %zu->%zu f, same_dir=%ld -> ok\n",
                    cnf, conf, same_dir);
        fails += c2fail;
    }

    /* Case 3: the --winding-repair path (skip_winding_repair=0) on a genuinely
     * PLANAR patch (reuse case 1's flat grid). Here the global-reference repair
     * is correct -- one outward normal, all faces agree -- so it must reduce the
     * mesh AND leave winding consistent (same_dir==0). Keeps the opt-in repair
     * branch covered; the repair is only wrong on WRAPPED input (case 2). */
    {
        QemOpts qo; memset(&qo, 0, sizeof qo);
        qo.skip_winding_repair = 0;      /* repair ON: the --winding-repair path */
        float *rov = NULL; int32_t *rof = NULL; size_t ronv = 0, ronf = 0;
        int rrc = QEM_simplify_opts(a, verts, nv, faces, nf, nf / 4, &qo,
                                    &rov, &ronv, &rof, &ronf, NULL);
        long same_dir = (rrc == 0 && ronf > 0) ? count_same_dir_he(a, rof, ronf) : -1;
        int c3fail = 0;
        if (rrc != 0)                { fprintf(stderr, "[selftest3] rc=%d -> FAIL\n", rrc); c3fail++; }
        if (ronf == 0 || ronf >= nf) { fprintf(stderr, "[selftest3] no reduction ronf=%zu -> FAIL\n", ronf); c3fail++; }
        if (same_dir != 0)           { fprintf(stderr, "[selftest3] winding-repair on planar left %ld same_dir -> FAIL\n", same_dir); c3fail++; }
        if (!c3fail)
            fprintf(stderr, "[selftest3] planar winding-repair %zu->%zu f, same_dir=%ld -> ok\n",
                    nf, ronf, same_dir);
        fails += c3fail;
    }

    Arena_dispose(&a);
    fprintf(stderr, "=== qslim_obj selftest %s (%d failure%s) ===\n",
            fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return run_selftest();
    if (argc < 4) {
        fprintf(stderr,
            "Usage: %s <in.obj> <out.obj> <keep_ratio>   (keep_ratio in (0,1), e.g. 0.10)\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    const char *in_path = argv[1];
    const char *out_path = argv[2];
    double ratio = atof(argv[3]);
    int do_remesh = 0, use_prob = 1, remesh_mls = 0, winding_repair = 0, ai;
    if (!(ratio > 0.0 && ratio < 1.0)) {
        fprintf(stderr, "qslim_obj: keep_ratio must be in (0,1), got %s\n", argv[3]);
        return 1;
    }
    /* Optional flags. Probabilistic quadrics default ON for this tool (the
     * hierarchical LOD path wants them); --no-prob restores deterministic QEM.
     * --remesh runs the post-decimation isotropic remesh quality pass.
     * Winding: QEM's global-reference winding-repair is SKIPPED by default. It
     * flips every face opposing one per-component "outward" reference normal --
     * on a surface wrapped around an axis (the scroll) ~half the faces oppose
     * ANY single reference, so the "repair" reverses whole regions and
     * manufactures same_dir edges. This tool's inputs are already consistently
     * wound (the weld runs OrientMesh first), so we preserve that winding.
     * --winding-repair re-enables the repair for a genuinely planar, mis-wound
     * input. */
    for (ai = 4; ai < argc; ai++) {
        if      (!strcmp(argv[ai], "--remesh"))         do_remesh = 1;
        else if (!strcmp(argv[ai], "--no-prob"))        use_prob = 0;
        else if (!strcmp(argv[ai], "--mls"))            remesh_mls = 1;
        else if (!strcmp(argv[ai], "--winding-repair")) winding_repair = 1;
        else { fprintf(stderr, "qslim_obj: unknown arg %s\n", argv[ai]); return 1; }
    }

    Arena_T arena = Arena_new();
    float *verts = NULL; int32_t *faces = NULL; size_t nv = 0, nf = 0;
    if (ObjIO_read(arena, in_path, &verts, &nv, &faces, &nf) != 0) {
        fprintf(stderr, "qslim_obj: cannot read %s\n", in_path);
        Arena_dispose(&arena); return 1;
    }
    size_t target_nf = (size_t)((double)nf * ratio + 0.5);
    if (target_nf < 4) target_nf = 4;   /* never collapse below a tetra's worth */

    /* Decimate. Request the compacted post-decimation pin mask (the auto-pinned
     * open boundary) so the remesh pass can freeze exactly that boundary. */
    QemOpts qo; memset(&qo, 0, sizeof qo);
    qo.prob_quadrics = use_prob;
    qo.skip_winding_repair = !winding_repair;   /* default: preserve input winding (see flags) */
    float *ov = NULL; int32_t *of = NULL; size_t onv = 0, onf = 0;
    uint8_t *opin = NULL;
    int rc = QEM_simplify_opts(arena, verts, nv, faces, nf, target_nf, &qo,
                               &ov, &onv, &of, &onf, do_remesh ? &opin : NULL);
    if (rc != 0) {
        fprintf(stderr, "qslim_obj: QEM_simplify failed (rc=%d) on %s\n", rc, in_path);
        Arena_dispose(&arena); return 1;
    }

    /* Isotropic remesh (fail-closed: falls back to the decimated mesh verbatim
     * if it would degrade topology, so onf/ov/of are always valid). */
    int remesh_rc = -1;
    if (do_remesh && onf > 0) {
        RemeshOpts ro; Remesh_default_opts(&ro);
        if (remesh_mls) ro.reproject = REMESH_REPROJECT_MLS;
        float *rv = NULL; int32_t *rf = NULL; size_t rnv = 0, rnf = 0;
        remesh_rc = Remesh_isotropic(arena, ov, onv, of, onf, opin, &ro,
                                     &rv, &rnv, &rf, &rnf);
        ov = rv; onv = rnv; of = rf; onf = rnf;
    }

    ves_ensure_parent_dir(out_path);
    if (ObjIO_write(out_path, ov, onv, of, onf) != 0) {
        fprintf(stderr, "qslim_obj: cannot write %s\n", out_path);
        Arena_dispose(&arena); return 1;
    }
    fprintf(stderr,
        "qslim_obj: %s  v %zu->%zu  f %zu->%zu (%.2f%% of faces, target %zu)%s%s%s -> %s\n",
        in_path, nv, onv, nf, onf,
        nf ? 100.0 * (double)onf / (double)nf : 0.0, target_nf,
        use_prob ? " prob" : "",
        winding_repair ? " wind-repair" : " wind-preserve",
        do_remesh ? (remesh_rc == 0 ? " remesh" : " remesh(fallback)") : "",
        out_path);
    Arena_dispose(&arena);
    return 0;
}
