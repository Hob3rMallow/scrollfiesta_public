/* obj_components.c -- connected-component breakdown + inter-component proximity.
 *
 * Diagnoses a BPA surface that fragmented into multiple connected components
 * (triangle islands sharing no edge) when the underlying point cloud is clearly
 * ONE sheet -- i.e. BPA failed to pivot a ball across a sub-voxel gap and left a
 * stitch open. obj_topology.ps1 reports only the component COUNT; this tool
 * reports, per component, its size and bounding box, and -- the crux -- the
 * minimum vertex-to-vertex GAP between every pair of components plus how many
 * vertex pairs sit within a "should-have-welded" threshold. A pair separated by
 * a fraction of a voxel over many vertex pairs is one sheet BPA dropped; a pair
 * separated by >= a scroll-wrap gap (>=7 vox) is correctly distinct.
 *
 * Components are labelled 1..C by descending face count (pipeline convention).
 * Coordinate convention matches the pipeline OBJs: "v c0 c1 c2 [r g b]" with
 * c0=z, c1=y, c2=x.
 *
 * Usage:
 *   obj_components <mesh.obj> [--gap <vox=1.5>] [--cell <vox=4.0>] [--dump <out.obj>]
 *   obj_components --selftest
 *
 *   --gap   near-pair threshold: vertex pairs closer than this are "weld-distance".
 *   --cell  proximity search window (= grid cell). Gaps are measured EXACTLY up to
 *           this distance; pairs farther apart are reported ">cell (separate)".
 *   --dump  write the mesh recoloured per component (distinct palette), with
 *           near-seam verts (within --gap of another component) brightened, for
 *           MeshLab inspection.
 *
 * Exit: 0 = ok / selftest pass, 1 = IO error, 2 = usage error, 3 = selftest fail.
 *
 * Standalone C99; no project deps (own tiny vcxproj). Modeled on pinhole_verdict.c.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define INF 1e300

/* ---- growable arrays (flattened triples), as in pinhole_verdict.c ---- */
typedef struct { double *p; size_t n, cap; } DArr;
typedef struct { int    *p; size_t n, cap; } IArr;
static void darr_push3(DArr *d, double x, double y, double z) {
    if (d->n + 3 > d->cap) { d->cap = d->cap ? d->cap * 2 : 3072;
        d->p = realloc(d->p, d->cap * sizeof(double)); }
    d->p[d->n++] = x; d->p[d->n++] = y; d->p[d->n++] = z;
}
static void iarr_push3(IArr *a, int x, int y, int z) {
    if (a->n + 3 > a->cap) { a->cap = a->cap ? a->cap * 2 : 3072;
        a->p = realloc(a->p, a->cap * sizeof(int)); }
    a->p[a->n++] = x; a->p[a->n++] = y; a->p[a->n++] = z;
}

/* ---- union-find with path halving ---- */
static int uf_find(int *uf, int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

/* ---- per-component record (sorted by descending face count) ---- */
typedef struct {
    int    id;            /* 1-based label */
    int    root;          /* union-find root vertex index */
    long   nf;            /* triangle count */
    long   nv;            /* referenced-vertex count */
    double bbmin[3], bbmax[3], cen[3];   /* (z,y,x) */
} Comp;
static int comp_cmp_nf(const void *a, const void *b) {
    const Comp *ca = a, *cb = b;
    if (ca->nf != cb->nf) return ca->nf < cb->nf ? 1 : -1;   /* descending */
    return 0;
}

/* ---- analysis result (for selftest assertions + --brief one-liner) ---- */
enum { K_SINGLE = 0, K_SEPARATE = 1, K_SEAM = 2, K_PARTIAL = 3, K_DOUBLE = 4 };
static const char *KLASS_TAG[5] = { "SINGLE", "SEPARATE", "SEAM-SPLIT", "PARTIAL-DOUBLE", "DOUBLE" };
typedef struct {
    int  comps;
    double closest_min;   /* smallest inter-comp gap found (vox); INF if all separate */
    int  closest_a, closest_b;   /* 1-based ids of the closest pair, -1 if none in window */
    long closest_near;           /* #vertex pairs within --gap for that pair */
    double cover_a, cover_b;     /* %% of each closest-pair island coincident with the other */
    int  klass;                  /* K_* classification of the closest pair */
    double crease_deg;           /* mean ACUTE angle between facing normals across the seam
                                  * (0 = coplanar continuation, big = a fold/crease); <0 = not computed */
} Result;

/* ===================================================================
 * Core analysis: union-find components, per-comp geometry, and a uniform
 * spatial grid that measures the exact min gap (<= cell) between every pair of
 * components and counts vertex pairs within --gap.
 * =================================================================== */
static Result analyze(const double *V, size_t nv, const int *F, size_t nf,
                      double gap, double cell, const char *dumppath, int verbose, int want_crease)
{
    Result R = { 0, INF, -1, -1, 0, 0.0, 0.0, K_SINGLE, -1.0 };
    if (nv == 0 || nf == 0) { if (verbose) printf("  (empty mesh: nv=%zu nf=%zu)\n", nv, nf); return R; }

    /* connected components via union-find over face edges */
    int *uf = malloc(nv * sizeof(int));
    for (size_t i = 0; i < nv; i++) uf[i] = (int)i;
    for (size_t fi = 0; fi < nf; fi++) {
        int a = F[fi*3], b = F[fi*3+1], c = F[fi*3+2];
        int ra = uf_find(uf,a), rb = uf_find(uf,b);
        if (ra != rb) uf[rb] = ra;
        ra = uf_find(uf,a); int rc = uf_find(uf,c);
        if (ra != rc) uf[rc] = ra;
    }

    /* referenced verts + per-root counts */
    char *ref = calloc(nv, 1);
    for (size_t i = 0; i < nf*3; i++) ref[F[i]] = 1;
    long *rootnf = calloc(nv, sizeof(long));
    long *rootnv = calloc(nv, sizeof(long));
    for (size_t fi = 0; fi < nf; fi++) rootnf[uf_find(uf, F[fi*3])]++;
    for (size_t i = 0; i < nv; i++) if (ref[i]) rootnv[uf_find(uf,(int)i)]++;

    /* collect + sort components by face count */
    int C = 0;
    for (size_t i = 0; i < nv; i++) if (ref[i] && uf_find(uf,(int)i) == (int)i) C++;
    Comp *comp = malloc((size_t)C * sizeof(Comp));
    int ci = 0;
    for (size_t i = 0; i < nv; i++) if (ref[i] && uf_find(uf,(int)i) == (int)i) {
        comp[ci].root = (int)i; comp[ci].nf = rootnf[i]; comp[ci].nv = rootnv[i];
        for (int k = 0; k < 3; k++) { comp[ci].bbmin[k] = INF; comp[ci].bbmax[k] = -INF; comp[ci].cen[k] = 0; }
        ci++;
    }
    qsort(comp, (size_t)C, sizeof(Comp), comp_cmp_nf);
    for (int k = 0; k < C; k++) comp[k].id = k + 1;

    /* root -> 1-based id, then vertex -> id */
    int *root2id = calloc(nv, sizeof(int));
    for (int k = 0; k < C; k++) root2id[comp[k].root] = comp[k].id;
    int *cov = calloc(nv, sizeof(int));   /* component-of-vertex; 0 = unreferenced */
    for (size_t i = 0; i < nv; i++) if (ref[i]) cov[i] = root2id[uf_find(uf,(int)i)];

    /* per-component bbox + centroid */
    for (size_t i = 0; i < nv; i++) if (ref[i]) {
        Comp *c = &comp[cov[i]-1];
        for (int k = 0; k < 3; k++) { double x = V[i*3+k];
            if (x < c->bbmin[k]) c->bbmin[k] = x;
            if (x > c->bbmax[k]) c->bbmax[k] = x;
            c->cen[k] += x; }
    }
    for (int k = 0; k < C; k++) for (int j = 0; j < 3; j++)
        comp[k].cen[j] /= (double)comp[k].nv;

    if (verbose) {
        printf("  components: %d   (labelled 1..%d by descending face count)\n", C, C);
        printf("    id   faces    verts      bbox z[lo,hi]        y[lo,hi]        x[lo,hi]    extent  centroid(z,y,x)\n");
        for (int k = 0; k < C; k++) {
            Comp *c = &comp[k];
            double ez = c->bbmax[0]-c->bbmin[0], ey = c->bbmax[1]-c->bbmin[1], ex = c->bbmax[2]-c->bbmin[2];
            double emax = ez > ey ? (ez > ex ? ez : ex) : (ey > ex ? ey : ex);
            printf("   %3d %7ld %8ld   [%7.1f,%7.1f] [%7.1f,%7.1f] [%7.1f,%7.1f] %6.1f  (%.1f,%.1f,%.1f)\n",
                   c->id, c->nf, c->nv,
                   c->bbmin[0], c->bbmax[0], c->bbmin[1], c->bbmax[1], c->bbmin[2], c->bbmax[2],
                   emax, c->cen[0], c->cen[1], c->cen[2]);
        }
    }
    R.comps = C;
    if (C < 2) {
        R.klass = K_SINGLE;
        if (verbose) printf("  single component -- nothing to compare.\n");
        free(uf); free(ref); free(rootnf); free(rootnv); free(comp); free(root2id); free(cov);
        return R;
    }

    /* ---- uniform spatial grid over referenced verts (cell = search window) ---- */
    double gmin[3] = { INF, INF, INF }, gmax[3] = { -INF, -INF, -INF };
    for (size_t i = 0; i < nv; i++) if (ref[i]) for (int k = 0; k < 3; k++) {
        double x = V[i*3+k]; if (x < gmin[k]) gmin[k] = x; if (x > gmax[k]) gmax[k] = x; }
    double cs = cell > 1e-6 ? cell : 1.0;
    int dim[3];
    size_t ncell;
    for (;;) {
        for (int k = 0; k < 3; k++) { dim[k] = (int)((gmax[k]-gmin[k])/cs) + 1; if (dim[k] < 1) dim[k] = 1; }
        ncell = (size_t)dim[0]*(size_t)dim[1]*(size_t)dim[2];
        if (ncell <= 50000000u) break;
        cs *= 2.0;   /* keep the grid bounded if cell was set tiny */
    }
    int  *cidx = malloc(nv * sizeof(int));        /* cell index per vert, -1 if unreferenced */
    long *coff = calloc(ncell + 1, sizeof(long)); /* CSR offsets */
    for (size_t i = 0; i < nv; i++) {
        if (!ref[i]) { cidx[i] = -1; continue; }
        int g[3];
        for (int k = 0; k < 3; k++) { int q = (int)((V[i*3+k]-gmin[k])/cs);
            if (q < 0) q = 0; if (q >= dim[k]) q = dim[k]-1; g[k] = q; }
        size_t ic = ((size_t)g[0]*(size_t)dim[1] + (size_t)g[1])*(size_t)dim[2] + (size_t)g[2];
        cidx[i] = (int)ic; coff[ic+1]++;
    }
    for (size_t i = 0; i < ncell; i++) coff[i+1] += coff[i];
    int  *cellv = malloc((size_t)coff[ncell] * sizeof(int));
    long *cur = malloc(ncell * sizeof(long));
    for (size_t i = 0; i < ncell; i++) cur[i] = coff[i];
    for (size_t i = 0; i < nv; i++) if (cidx[i] >= 0) cellv[cur[cidx[i]]++] = (int)i;
    free(cur);

    /* ---- proximity scan: per-vertex nearest cross-comp neighbour + per-pair stats ---- */
    double *pmin2 = malloc((size_t)C*(size_t)C * sizeof(double));
    long   *pcnt  = calloc((size_t)C*(size_t)C, sizeof(long));
    for (size_t i = 0; i < (size_t)C*(size_t)C; i++) pmin2[i] = INF;
    double *nn2 = malloc(nv * sizeof(double));
    int    *nnc = calloc(nv, sizeof(int));
    int    *nn_idx = malloc(nv * sizeof(int));   /* nearest cross-comp vertex index */
    for (size_t i = 0; i < nv; i++) { nn2[i] = INF; nn_idx[i] = -1; }
    double gap2 = gap*gap;
    int planeYZ = dim[1]*dim[2];
    for (size_t i = 0; i < nv; i++) {
        if (cidx[i] < 0) continue;
        int ai = cov[i];
        double xi = V[i*3], yi = V[i*3+1], zi = V[i*3+2];
        int g0 = cidx[i] / planeYZ, rem = cidx[i] % planeYZ, g1 = rem / dim[2], g2 = rem % dim[2];
        for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
            int n0 = g0+dz, n1 = g1+dy, n2 = g2+dx;
            if (n0 < 0 || n1 < 0 || n2 < 0 || n0 >= dim[0] || n1 >= dim[1] || n2 >= dim[2]) continue;
            size_t nc = ((size_t)n0*(size_t)dim[1] + (size_t)n1)*(size_t)dim[2] + (size_t)n2;
            for (long t = coff[nc]; t < coff[nc+1]; t++) {
                int w = cellv[t];
                int aw = cov[w];
                if (aw == ai) continue;                 /* same component */
                double ddx = V[w*3]-xi, ddy = V[w*3+1]-yi, ddz = V[w*3+2]-zi;
                double d2 = ddx*ddx + ddy*ddy + ddz*ddz;
                if (d2 < nn2[i]) { nn2[i] = d2; nnc[i] = aw; nn_idx[i] = w; }
                if ((size_t)w > i) {                     /* canonical: count each pair once */
                    int lo = ai < aw ? ai : aw, hi = ai < aw ? aw : ai;
                    size_t pij = (size_t)(lo-1)*(size_t)C + (size_t)(hi-1);
                    if (d2 < pmin2[pij]) pmin2[pij] = d2;
                    if (d2 <= gap2) pcnt[pij]++;
                }
            }
        }
    }

    /* ---- report inter-component gaps ---- */
    if (verbose)
        printf("\n  inter-component gaps (exact within search window cell=%.1f vox; near-pair threshold gap=%.2f vox):\n", cs, gap);
    double best = INF; int ba = -1, bb = -1; long bnear = 0;
    for (int a = 0; a < C; a++) for (int b = a+1; b < C; b++) {
        size_t pij = (size_t)a*(size_t)C + (size_t)b;
        if (pmin2[pij] < INF/2) {
            double mn = sqrt(pmin2[pij]);
            if (verbose)
                printf("    comp %d <-> comp %d : min gap %.4f vox,  %ld vertex pair%s within %.2f vox\n",
                       a+1, b+1, mn, pcnt[pij], pcnt[pij]==1?"":"s", gap);
            if (mn < best) { best = mn; ba = a+1; bb = b+1; bnear = pcnt[pij]; }
        } else if (verbose) {
            printf("    comp %d <-> comp %d : gap > %.1f vox (no vertices within window -- SEPARATE)\n", a+1, b+1, cs);
        }
    }
    R.closest_min = best; R.closest_a = ba; R.closest_b = bb; R.closest_near = bnear;
    if (ba < 0) R.klass = K_SEPARATE;

    /* ---- coincidence coverage for the closest pair ----
     * What fraction of each island has a partner-island vertex within --gap?
     * HIGH on both sides => the islands occupy the SAME surface region (a doubled
     * / shingled overlap -- a z-fighting double the merge must collapse). LOW
     * (only a thin band) => they tile complementary area and meet along a SEAM
     * (one layer BPA's front split; the fix is to bridge, not collapse). */
    if (ba > 0) {
        long cov_a = 0, cov_b = 0; double sda = 0, sdb = 0;
        for (size_t i = 0; i < nv; i++) {
            if (!ref[i] || nn2[i] > gap2) continue;
            if (cov[i] == ba && nnc[i] == bb) { cov_a++; sda += sqrt(nn2[i]); }
            else if (cov[i] == bb && nnc[i] == ba) { cov_b++; sdb += sqrt(nn2[i]); }
        }
        long na = comp[ba-1].nv, nb = comp[bb-1].nv;
        double fa = na ? 100.0*(double)cov_a/(double)na : 0.0;
        double fb = nb ? 100.0*(double)cov_b/(double)nb : 0.0;
        R.cover_a = fa; R.cover_b = fb;
        R.klass = (fa > 50.0 && fb > 50.0) ? K_DOUBLE
                : (fa > 50.0 || fb > 50.0) ? K_PARTIAL : K_SEAM;
        if (verbose) {
            printf("\n  closest-pair coincidence (within %.2f vox):\n", gap);
            printf("    comp %d: %ld/%ld verts (%.0f%%) coincide with comp %d  (mean %.3f vox)\n",
                   ba, cov_a, na, fa, bb, cov_a ? sda/(double)cov_a : 0.0);
            printf("    comp %d: %ld/%ld verts (%.0f%%) coincide with comp %d  (mean %.3f vox)\n",
                   bb, cov_b, nb, fb, ba, cov_b ? sdb/(double)cov_b : 0.0);
            printf("    -> %s\n", R.klass == K_DOUBLE
                   ? "high both sides: islands OVERLAP the same region (doubled/shingled sheet)"
                   : R.klass == K_PARTIAL
                     ? "high one side: one island sits atop part of the other (partial double)"
                     : "thin band only: they meet along a SEAM (one split layer)");
        }
    }

    /* ---- crease / dihedral test for the closest pair ----
     * Is the seam a sharp FOLD (the two islands are the two faces of a crease,
     * with surface normals swinging across the line) or a flat CONTINUATION (the
     * normals stay parallel and the split is purely a sampling/empty-ball stall)?
     * Normals come from face WINDING (no MLS needed). We measure, per near-seam
     * vertex pair, the angle between the two facing vertex normals, folded to the
     * ACUTE range [0,90] so it is independent of the two fronts' winding sign:
     * ~0 = coplanar continuation, large = a real fold. */
    if (want_crease && ba > 0) {
        double *vn = (double *)calloc(nv * 3, sizeof(double));   /* area-weighted vertex normals */
        for (size_t fi = 0; fi < nf; fi++) {
            int i0 = F[fi*3], i1 = F[fi*3+1], i2 = F[fi*3+2];
            double e1[3], e2[3], fnv[3];
            for (int k = 0; k < 3; k++) { e1[k] = V[i1*3+k]-V[i0*3+k]; e2[k] = V[i2*3+k]-V[i0*3+k]; }
            fnv[0] = e1[1]*e2[2]-e1[2]*e2[1];
            fnv[1] = e1[2]*e2[0]-e1[0]*e2[2];
            fnv[2] = e1[0]*e2[1]-e1[1]*e2[0];           /* magnitude = 2*area (area weight) */
            for (int k = 0; k < 3; k++) { vn[i0*3+k]+=fnv[k]; vn[i1*3+k]+=fnv[k]; vn[i2*3+k]+=fnv[k]; }
        }
        for (size_t i = 0; i < nv; i++) {
            double *p = &vn[i*3], L = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
            if (L > 1e-12) { p[0]/=L; p[1]/=L; p[2]/=L; }
        }
        static const double DEG = 57.29577951308232;
        static const double ab_edge[7] = { 5,15,30,45,60,75,90 };
        long hist[7] = {0,0,0,0,0,0,0}, ncr = 0;
        double sum_acute = 0, sum_raw = 0;
        double mna[3] = {0,0,0}, mnb[3] = {0,0,0};   /* seam-band mean normals, each side */
        for (size_t i = 0; i < nv; i++) {
            if (!ref[i] || nn2[i] > gap2 || nn_idx[i] < 0) continue;
            int vci = cov[i], vcj = cov[nn_idx[i]];
            int ab = (vci == ba && vcj == bb), ba_ = (vci == bb && vcj == ba);
            if (!ab && !ba_) continue;
            double *ni = &vn[i*3], *nj = &vn[nn_idx[i]*3];
            double dot = ni[0]*nj[0]+ni[1]*nj[1]+ni[2]*nj[2];
            if (dot > 1) dot = 1; if (dot < -1) dot = -1;
            double adot = dot < 0 ? -dot : dot;
            double acute = acos(adot)*DEG, raw = acos(dot)*DEG;
            sum_acute += acute; sum_raw += raw; ncr++;
            int bk = 0; while (bk < 6 && acute >= ab_edge[bk]) bk++; hist[bk]++;
            if (ab) for (int k = 0; k < 3; k++) mna[k] += ni[k];
            else    for (int k = 0; k < 3; k++) mnb[k] += nj[k];
        }
        double mean_acute = ncr ? sum_acute/(double)ncr : -1.0;
        R.crease_deg = mean_acute;
        /* bulk angle between the two seam-band mean normals (acute) */
        double la = sqrt(mna[0]*mna[0]+mna[1]*mna[1]+mna[2]*mna[2]);
        double lb = sqrt(mnb[0]*mnb[0]+mnb[1]*mnb[1]+mnb[2]*mnb[2]);
        double bulk = -1.0;
        if (la > 1e-9 && lb > 1e-9) {
            double d = (mna[0]*mnb[0]+mna[1]*mnb[1]+mna[2]*mnb[2])/(la*lb);
            if (d < 0) d = -d; if (d > 1) d = 1;
            bulk = acos(d)*DEG;
        }
        if (verbose) {
            printf("\n  crease/dihedral test (closest pair comp %d <-> comp %d, face-winding normals):\n", ba, bb);
            printf("    facing-normal angle over %ld seam vertex pairs: mean ACUTE = %.1f deg (raw mean = %.1f deg)\n",
                   ncr, mean_acute, ncr ? sum_raw/(double)ncr : 0.0);
            printf("    acute hist [0-5 5-15 15-30 30-45 45-60 60-75 75-90]: %ld %ld %ld %ld %ld %ld %ld\n",
                   hist[0],hist[1],hist[2],hist[3],hist[4],hist[5],hist[6]);
            printf("    seam-band mean-normal angle: %.1f deg (acute)\n", bulk);
            printf("    -> %s\n", mean_acute < 10.0
                   ? "COPLANAR continuation (< 10 deg) -- NOT a crease; the split is a sampling/empty-ball stall"
                   : mean_acute < 25.0
                     ? "mild bend (10-25 deg) -- borderline"
                     : "CREASE / fold CONFIRMED (>= 25 deg) -- the seam runs along a sharp surface bend");
        }
        free(vn);
    }

    /* ---- verdict ---- */
    if (verbose) {
        printf("\n  VERDICT:\n");
        if (ba < 0) {
            printf("    all %d components are > %.1f vox apart -- genuinely separate sheets.\n", C, cs);
        } else {
            Comp *ca = &comp[ba-1], *cb = &comp[bb-1];
            if (best <= gap) {
                printf("    comp %d (%ld f) and comp %d (%ld f) meet at %.4f vox over %ld vertex pairs\n"
                       "    (<= weld distance %.2f vox). They are ONE sheet that BPA failed to stitch:\n"
                       "    the points coincide but no triangle bridges the islands.\n",
                       ba, ca->nf, bb, cb->nf, best, bnear, gap);
            } else if (best <= 2.5) {
                printf("    closest pair comp %d (%ld f) <-> comp %d (%ld f) is %.4f vox apart (just beyond\n"
                       "    weld %.2f vox) -- almost certainly one sheet; the BPA ball (~rho) couldn't bridge it.\n",
                       ba, ca->nf, bb, cb->nf, best, gap);
            } else {
                printf("    closest pair comp %d <-> comp %d is %.4f vox apart.\n", ba, bb, best);
            }
        }
    }

    /* ---- optional per-component colored dump ---- */
    if (dumppath) {
        static const double PAL[6][3] = {
            {0.75,0.75,0.78}, {0.15,0.85,0.25}, {0.90,0.15,0.15},
            {0.20,0.45,0.95}, {0.95,0.75,0.10}, {0.65,0.25,0.85} };
        FILE *of = fopen(dumppath, "w");
        if (!of) { printf("  ERROR: cannot write dump %s\n", dumppath); }
        else {
            for (size_t i = 0; i < nv; i++) {
                double col[3] = { 0.40, 0.40, 0.40 };
                if (ref[i]) {
                    const double *pc = PAL[(cov[i]-1) % 6];
                    col[0] = pc[0]; col[1] = pc[1]; col[2] = pc[2];
                    if (nn2[i] <= gap2) { col[0] = 1.0; col[1] = 1.0; col[2] = 0.15; }  /* near-seam */
                }
                fprintf(of, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                        V[i*3], V[i*3+1], V[i*3+2], col[0], col[1], col[2]);
            }
            for (size_t fi = 0; fi < nf; fi++)
                fprintf(of, "f %d %d %d\n", F[fi*3]+1, F[fi*3+1]+1, F[fi*3+2]+1);
            fclose(of);
            if (verbose)
                printf("\n  wrote per-component colored OBJ -> %s\n"
                       "   (comp1 gray, comp2 green, comp3 red, comp4 blue, ...; near-seam verts within %.2f vox = bright yellow)\n",
                       dumppath, gap);
        }
    }

    free(uf); free(ref); free(rootnf); free(rootnv); free(comp); free(root2id); free(cov);
    free(cidx); free(coff); free(cellv); free(pmin2); free(pcnt); free(nn2); free(nnc); free(nn_idx);
    return R;
}

/* ===================================================================
 * Self-test: synthetic meshes with KNOWN component structure and gaps.
 * =================================================================== */
static int selftest(void)
{
    int fails = 0;
    printf("=== obj_components selftest ===\n");

    /* case 1: two triangles sharing an edge -> ONE connected component */
    {
        double V[] = { 0,0,0,  0,0,1,  0,1,0,  0,1,1 };
        int    F[] = { 0,1,2,  1,3,2 };
        printf("[case1] connected quad (expect comps=1)\n");
        Result r = analyze(V, 4, F, 2, 1.5, 4.0, NULL, 1, 0);
        if (r.comps != 1) { printf("  FAIL: comps=%d != 1\n", r.comps); fails++; }
        else printf("  ok: comps=1\n");
    }
    /* case 2: two triangles 10 vox apart -> 2 comps, beyond window = SEPARATE */
    {
        double V[] = { 0,0,0,  0,0,1,  0,1,0,   10,0,0,  10,0,1,  10,1,0 };
        int    F[] = { 0,1,2,  3,4,5 };
        printf("[case2] two far triangles (expect comps=2, separate)\n");
        Result r = analyze(V, 6, F, 2, 1.5, 4.0, NULL, 1, 0);
        if (!(r.comps == 2 && r.closest_a < 0)) {
            printf("  FAIL: comps=%d closest_a=%d (expect 2, -1)\n", r.comps, r.closest_a); fails++;
        } else printf("  ok: comps=2, correctly separate\n");
    }
    /* case 3: two triangles 0.5 vox apart, disjoint verts -> 2 comps, SHOULD-JOIN */
    {
        double V[] = { 0,0,0,  0,0,1,  0,1,0,   0.5,0,0,  0.5,0,1,  0.5,1,0 };
        int    F[] = { 0,1,2,  3,4,5 };
        printf("[case3] two near triangles @0.5 vox (expect comps=2, gap~0.5, near>0)\n");
        Result r = analyze(V, 6, F, 2, 1.5, 4.0, NULL, 1, 0);
        if (!(r.comps == 2 && r.closest_a > 0 && r.closest_min < 1.5 && r.closest_near > 0)) {
            printf("  FAIL: comps=%d closest_min=%.3f near=%ld\n", r.comps, r.closest_min, r.closest_near); fails++;
        } else printf("  ok: comps=2 gap=%.3f near=%ld (BPA-miss signature)\n", r.closest_min, r.closest_near);
    }
    /* case 4: crease -- A in z=0 plane (normal z), B in x=0.3 plane (normal x),
     * B within gap of A -> facing normals 90 deg apart. Expect crease_deg ~90. */
    {
        double V[] = { 0,0,0,  0,2,0,  0,0,2,    1,0,0.3,  1,2,0.3,  3,0,0.3 };
        int    F[] = { 0,1,2,  3,4,5 };
        printf("[case4] 90-deg crease (expect crease_deg ~90)\n");
        Result r = analyze(V, 6, F, 2, 1.5, 4.0, NULL, 1, 1);
        if (!(r.comps == 2 && r.crease_deg > 80.0 && r.crease_deg <= 90.0001)) {
            printf("  FAIL: comps=%d crease_deg=%.1f (expect ~90)\n", r.comps, r.crease_deg); fails++;
        } else printf("  ok: crease_deg=%.1f (fold confirmed)\n", r.crease_deg);
    }
    /* case 5: coplanar continuation -- A and B both in z=0 plane, B within gap of
     * A -> facing normals parallel. Expect crease_deg ~0. */
    {
        double V[] = { 0,0,0,  0,2,0,  0,0,2,    0,0,2.3,  0,2,2.3,  0,0,4 };
        int    F[] = { 0,1,2,  3,4,5 };
        printf("[case5] coplanar continuation (expect crease_deg ~0)\n");
        Result r = analyze(V, 6, F, 2, 1.5, 4.0, NULL, 1, 1);
        if (!(r.comps == 2 && r.crease_deg >= 0.0 && r.crease_deg < 10.0)) {
            printf("  FAIL: comps=%d crease_deg=%.1f (expect ~0)\n", r.comps, r.crease_deg); fails++;
        } else printf("  ok: crease_deg=%.1f (coplanar, not a crease)\n", r.crease_deg);
    }

    printf("=== selftest %s (%d failure%s) ===\n",
           fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s <mesh.obj> [--gap <vox=1.5>] [--cell <vox=4.0>] [--dump <out.obj>] [--brief] [--crease]\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    const char *path = argv[1];
    double gap = 1.5, cell = 4.0;
    const char *dump = NULL;
    int brief = 0, want_crease = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--gap")  && i+1 < argc) gap  = atof(argv[++i]);
        else if (!strcmp(argv[i], "--cell") && i+1 < argc) cell = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dump") && i+1 < argc) dump = argv[++i];
        else if (!strcmp(argv[i], "--brief")) brief = 1;
        else if (!strcmp(argv[i], "--crease")) want_crease = 1;
        else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 2; }
    }
    if (gap > cell) { fprintf(stderr, "note: --gap (%.2f) > --cell (%.2f); raising cell to gap\n", gap, cell); cell = gap; }

    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", path); return 1; }
    DArr V = {0}; IArr Fc = {0};
    char line[8192];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            double x = 0, y = 0, z = 0;
            if (sscanf(line+2, "%lf %lf %lf", &x, &y, &z) == 3) darr_push3(&V, x, y, z);
        } else if (line[0] == 'f' && line[1] == ' ') {
            int idx[64]; int ni = 0;
            const char *p = line + 1;
            while (*p && ni < 64) {
                while (*p == ' ' || *p == '\t') p++;
                if (*p == '\0' || *p == '\n' || *p == '\r') break;
                int vi = (int)strtol(p, (char**)&p, 10);
                while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;  /* skip /vt/vn */
                if (vi < 0) vi = (int)(V.n/3) + vi; else vi -= 1;
                idx[ni++] = vi;
            }
            for (int i = 1; i+1 < ni; i++) iarr_push3(&Fc, idx[0], idx[i], idx[i+1]);
        }
    }
    fclose(f);

    size_t nv = V.n/3, nf = Fc.n/3;
    if (!brief) {
        printf("=== obj_components: %s ===\n", path);
        printf("  V=%zu  F=%zu  gap=%.2f vox  cell=%.2f vox\n", nv, nf, gap, cell);
    }
    Result r = analyze(V.p, nv, Fc.p, nf, gap, cell, dump, !brief, want_crease);

    if (brief) {
        /* one line per file, for scanning a corpus. basename only. */
        const char *base = path;
        for (const char *s = path; *s; s++) if (*s == '/' || *s == '\\') base = s + 1;
        char crease[32] = "";
        if (want_crease && r.crease_deg >= 0.0) snprintf(crease, sizeof(crease), " crease=%.0fdeg", r.crease_deg);
        if (r.comps < 2) {
            printf("%-52s comps=%d  %s\n", base, r.comps, KLASS_TAG[K_SINGLE]);
        } else if (r.closest_a < 0) {
            printf("%-52s comps=%d  %s\n", base, r.comps, KLASS_TAG[r.klass]);
        } else {
            printf("%-52s comps=%d  c%d-c%d gap=%.3f near=%ld cover=%.0f/%.0f%%%s  %s\n",
                   base, r.comps, r.closest_a, r.closest_b, r.closest_min,
                   r.closest_near, r.cover_a, r.cover_b, crease, KLASS_TAG[r.klass]);
        }
    }

    free(V.p); free(Fc.p);
    return 0;
}
