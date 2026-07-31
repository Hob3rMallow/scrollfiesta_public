/*
 * cvt_simplify -- CVT-decimate a welded block, per connected component,
 * boundary-preserving. The quality-preserving replacement for qslim in the
 * hierarchical LOD pyramid (hierarchical_weld --cvt-simplify).
 *
 *   cvt_simplify <in.obj> <out.obj> --keep-ratio R
 *                [--iters N] [--min-faces M] [--seed S]
 *
 * qslim decimates fast but re-scars every LOD tier with slivers, undoing the
 * band-CVT weld quality. Instead: split the welded block into connected
 * components (scroll wraps), and uniformly CVT-remesh each to keep_ratio of its
 * vertices. CVT's boundary seeds sit ON each component's boundary loop and are
 * projected back to it every iteration, so the block's OUTER boundary survives
 * (resampled but geometrically on the same polyline) as a weldable seam -- the
 * next pyramid level's seam_refine re-refines it and the bridge re-welds it,
 * exactly the leaf-weld path. Per-component isolation makes inter-wrap fusion
 * impossible (CVT never sees two wraps at once). ManifoldGuard cleans the
 * occasional RVD-dual bowtie. Fail-closed: a component whose CVT errors or
 * returns empty is kept at full resolution (never dropped).
 *
 * Because every level stays coarse, grid_weld's interior hole-fill never faces
 * the million-loop meshes that made it hang at undecimated upper levels -- it
 * stays functional at every level, which is the point.
 */
#include "../remesh/cvt_remesh.h"
#include "../remesh/manifold_guard.h"
#include "../common/obj_io.h"
#include "../common/arena.h"
#include "../common/except.h"
#include "../common/mesh_types.h"
#include "../common/pipeline_constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* union-find with path halving */
static int32_t uf_find(int32_t *par, int32_t x)
{
    while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
    return x;
}

/* (component root, face id) grouped by root so a component's faces are a
 * contiguous run after sorting. */
typedef struct { int32_t root, face; } FR;
static int fr_cmp(const void *pa, const void *pb)
{
    const FR *a = (const FR *)pa, *b = (const FR *)pb;
    if (a->root != b->root) return a->root < b->root ? -1 : 1;
    return a->face < b->face ? -1 : (a->face > b->face ? 1 : 0);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr,
            "usage: %s <in.obj> <out.obj> [--keep-ratio R=0.25] [--iters N=12]\n"
            "          [--min-faces M=64] [--seed S=1]\n", argv[0]);
        return 2;
    }
    const char *in_path = argv[1], *out_path = argv[2];
    double keep_ratio = 0.25;
    int    iters = CVT_PIPELINE_ITERS;
    size_t min_faces = 64;
    uint32_t seed = 1;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--keep-ratio") && i + 1 < argc) keep_ratio = atof(argv[++i]);
        else if (!strcmp(argv[i], "--iters") && i + 1 < argc) iters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min-faces") && i + 1 < argc) min_faces = (size_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = (uint32_t)atoll(argv[++i]);
        else { fprintf(stderr, "cvt_simplify: unknown arg %s\n", argv[i]); return 2; }
    }
    if (!(keep_ratio > 0.0 && keep_ratio <= 1.0)) { fprintf(stderr, "keep-ratio in (0,1]\n"); return 2; }

    Arena_T arena = Arena_new();
    int rc_exit = 1;

    TRY
        float *V = NULL; int32_t *F = NULL; size_t nv = 0, nf = 0;
        if (ObjIO_read(arena, in_path, &V, &nv, &F, &nf) != 0 || nv < 3 || nf < 1) {
            fprintf(stderr, "cvt_simplify: read failed / empty: %s\n", in_path);
            RAISE(IO_Failed);
        }

        /* connected components over shared verts */
        int32_t *par = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
        for (size_t v = 0; v < nv; v++) par[v] = (int32_t)v;
        for (size_t f = 0; f < nf; f++) {
            int32_t a = F[f*3+0], b = F[f*3+1], c = F[f*3+2];
            int32_t ra = uf_find(par, a), rb = uf_find(par, b);
            if (rb != ra) par[rb] = ra;
            int32_t rc = uf_find(par, c);
            ra = uf_find(par, a);
            if (rc != ra) par[rc] = ra;
        }
        /* face -> component root, sorted so each component is a contiguous run */
        FR *fr = (FR *)ARENA_ALLOC(arena, (long)(nf * sizeof(FR)));
        for (size_t f = 0; f < nf; f++) { fr[f].root = uf_find(par, F[f*3+0]); fr[f].face = (int32_t)f; }
        qsort(fr, nf, sizeof(FR), fr_cmp);
        size_t nroot = 0;
        for (size_t f = 0; f < nf; f++) if (f == 0 || fr[f].root != fr[f-1].root) nroot++;

        /* per component: gather local mesh, CVT-decimate (fail-closed). Store
         * results as a ComponentMesh array for the guard, then concatenate. */
        ComponentMesh *cms = (ComponentMesh *)ARENA_CALLOC(arena, (long)nroot, (long)sizeof(ComponentMesh));
        size_t ncm = 0;
        int32_t *g2l = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
        for (size_t v = 0; v < nv; v++) g2l[v] = -1;

        size_t n_cvt = 0, n_kept = 0;
        size_t pi = 0;
        while (pi < nf) {
            size_t pj = pi + 1;
            while (pj < nf && fr[pj].root == fr[pi].root) pj++;
            size_t cf = pj - pi;
            int32_t *lf = (int32_t *)ARENA_ALLOC(arena, (long)(cf * 3 * sizeof(int32_t)));
            int32_t *l2g = (int32_t *)ARENA_ALLOC(arena, (long)(cf * 3 * sizeof(int32_t)));
            size_t lnv = 0;
            for (size_t k = 0; k < cf; k++) {
                for (int e = 0; e < 3; e++) {
                    int32_t gv = F[(size_t)fr[pi+k].face*3+(size_t)e];
                    if (g2l[gv] < 0) { g2l[gv] = (int32_t)lnv; l2g[lnv] = gv; lnv++; }
                    lf[k*3+(size_t)e] = g2l[gv];
                }
            }
            float *lv = (float *)ARENA_ALLOC(arena, (long)(lnv * 3 * sizeof(float)));
            for (size_t v = 0; v < lnv; v++) {
                lv[v*3+0] = V[(size_t)l2g[v]*3+0];
                lv[v*3+1] = V[(size_t)l2g[v]*3+1];
                lv[v*3+2] = V[(size_t)l2g[v]*3+2];
            }
            for (size_t v = 0; v < lnv; v++) g2l[l2g[v]] = -1;   /* reset */

            float *ov = lv; int32_t *of = lf; size_t onv = lnv, onf = cf;
            if (cf >= min_faces) {
                size_t target = (size_t)(keep_ratio * (double)lnv + 0.5);
                if (target < CVT_MIN_SITES) target = CVT_MIN_SITES;
                if (target < lnv) {                 /* only if it actually coarsens */
                    CvtOpts co; CVT_default_opts(&co);
                    co.n_iters = iters; co.seed = seed + (uint32_t)ncm;
                    float *nv2 = NULL; int32_t *nf2 = NULL; size_t nnv = 0, nnf = 0;
                    int crc = CVT_remesh(arena, lv, lnv, lf, cf, target, &co, NULL,
                                         &nv2, &nnv, &nf2, &nnf);
                    if (crc == 0 && nnf > 0) { ov = nv2; of = nf2; onv = nnv; onf = nnf; n_cvt++; }
                    else n_kept++;
                } else n_kept++;
            } else n_kept++;

            cms[ncm].verts = ov; cms[ncm].faces = of; cms[ncm].nv = onv; cms[ncm].nf = onf;
            cms[ncm].comp_id = (int)ncm + 1; cms[ncm].self = &cms[ncm];
            ncm++;
            pi = pj;
        }

        /* ManifoldGuard each component (resolve RVD-dual bowties / NM edges,
         * reorient). Per-component => no cross-wrap interaction. */
        for (size_t i = 0; i < ncm; i++) {
            if (cms[i].nf == 0) continue;
            ManifoldGuardStats mg;
            ManifoldGuard_process(arena, &cms[i], 1, 1 /*reorient*/, &mg);
        }

        /* concatenate components into one mesh */
        size_t out_nv = 0, out_nf = 0;
        for (size_t i = 0; i < ncm; i++) { out_nv += cms[i].nv; out_nf += cms[i].nf; }
        float *OV = (float *)ARENA_ALLOC(arena, (long)(out_nv * 3 * sizeof(float)));
        int32_t *OF = (int32_t *)ARENA_ALLOC(arena, (long)(out_nf * 3 * sizeof(int32_t)));
        size_t voff = 0, foff = 0;
        for (size_t i = 0; i < ncm; i++) {
            for (size_t v = 0; v < cms[i].nv; v++) {
                OV[(voff+v)*3+0] = cms[i].verts[v*3+0];
                OV[(voff+v)*3+1] = cms[i].verts[v*3+1];
                OV[(voff+v)*3+2] = cms[i].verts[v*3+2];
            }
            for (size_t f = 0; f < cms[i].nf; f++) {
                OF[(foff+f)*3+0] = cms[i].faces[f*3+0] + (int32_t)voff;
                OF[(foff+f)*3+1] = cms[i].faces[f*3+1] + (int32_t)voff;
                OF[(foff+f)*3+2] = cms[i].faces[f*3+2] + (int32_t)voff;
            }
            voff += cms[i].nv; foff += cms[i].nf;
        }

        if (ObjIO_write(out_path, OV, out_nv, OF, out_nf) != 0) {
            fprintf(stderr, "cvt_simplify: write failed: %s\n", out_path);
            RAISE(IO_Failed);
        }
        fprintf(stderr,
            "cvt_simplify: %zu comps (%zu cvt, %zu kept), %zu -> %zu faces "
            "(%.0f%%), %zu -> %zu verts, keep=%.3f\n",
            ncm, n_cvt, n_kept, nf, out_nf,
            nf ? 100.0 * (double)out_nf / (double)nf : 0.0, nv, out_nv, keep_ratio);
        rc_exit = 0;
    EXCEPT(IO_Failed)
        rc_exit = 1;
    EXCEPT(Arena_Failed)
        fprintf(stderr, "cvt_simplify: OOM\n");
        rc_exit = 1;
    END_TRY;

    Arena_dispose(&arena);
    return rc_exit;
}
