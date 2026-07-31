/* crossseam_fill.c -- find, rank, and dump the interior hole-fills that the grid
 * weld added ACROSS a cube seam.
 *
 * Background. grid_weld concatenates the per-cube charts at world origins, bridges
 * the ~1-vox seam gap with BPA, then runs HoleFill_process_ex(interior_only) over
 * the merged mesh. Most filled holes sit inside a single cube; the interesting ones
 * straddle a seam plane -- their boundary loop has rim vertices owned by two adjacent
 * cubes and the fill triangulation bridges the seam. `--dump-stages` writes the weld
 * pipeline OBJs, of which two matter here:
 *     ..._06_cleanup.obj   (input to interior hole-fill)
 *     ..._07_holefill.obj  (output of interior hole-fill)
 * Holefill APPENDS Steiner vertices but REORDERS faces, so the added triangles are
 * recovered by a set-difference on sorted vertex triples, NOT by taking a tail.
 *
 * What it does.
 *   1. PRE faces -> sorted-triple set (qsort + bsearch, exact, no hash collisions).
 *   2. POST faces not in that set  == the triangles hole-fill added (the "fill faces").
 *   3. union-find the fill faces by shared vertex  == one patch per filled hole.
 *   4. per patch: rim (boundary) verts, Steiner verts, bbox, single-loop + manifold
 *      check, and -- for each candidate seam plane -- whether the patch straddles it.
 *   5. rank patches by a "clean cross-seam fill" score and print the leaders, naming
 *      the two cubes the seam divides.
 *   6. --dump <id> writes that patch (red) plus the surrounding POST mesh (gray) as a
 *      vertex-coloured OBJ you can drop straight into a viewer.
 *
 * Coordinate convention: pipeline OBJ vertices are (Z, Y, X) world voxels -- axis
 * 0 = Z, 1 = Y, 2 = X. Cube origins are multiples of 128, so a vertex's owning cube
 * origin on an axis is floor(coord/128)*128, and the internal seam planes are the
 * multiples of 128 strictly inside the mesh bbox.
 *
 * Usage:
 *   crossseam_fill <pre.obj> <post.obj> [opts]
 *       --seam <Z|Y|X>=<coord>   seam plane (repeatable; default = auto, every
 *                                multiple of 128 strictly inside the bbox)
 *       --band <vox>             rim verts within this of a plane count as "on" it
 *                                when measuring the seam gap        (default 4.0)
 *       --min <n> --max <n>      rim-loop size window for "good"     (default 4..400)
 *       --top <k>                how many ranked patches to print    (default 20)
 *       --dump <out.obj>         write a coloured crop ...
 *       --patch <id>             ... of this patch id (1-based from the ranking)
 *       --radius <vox>           crop neighbourhood radius           (default 14.0)
 *       --all                    rank ALL fill patches, not just cross-seam ones
 *   crossseam_fill --selftest
 *
 * Standalone C99; no project deps so it builds as its own tiny vcxproj.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/* growable arrays (flattened), as in pinhole_verdict.c / obj_components.c */
typedef struct { double *p; size_t n, cap; } DArr;
typedef struct { int    *p; size_t n, cap; } IArr;
static void darr_push3(DArr *d, double x, double y, double z) {
    if (d->n + 3 > d->cap) { d->cap = d->cap ? d->cap * 2 : 4096;
        d->p = realloc(d->p, d->cap * sizeof(double)); }
    d->p[d->n++] = x; d->p[d->n++] = y; d->p[d->n++] = z;
}
static void iarr_push3(IArr *a, int x, int y, int z) {
    if (a->n + 3 > a->cap) { a->cap = a->cap ? a->cap * 2 : 4096;
        a->p = realloc(a->p, a->cap * sizeof(int)); }
    a->p[a->n++] = x; a->p[a->n++] = y; a->p[a->n++] = z;
}

/* ------------------------------------------------------------------ */
/* OBJ loader. Reads vertex xyz (ignores trailing rgb), triangulates polygon
 * fans, resolves negative indices. Returns 0 on success. */
static int load_obj(const char *path, DArr *V, IArr *F) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", path); return 1; }
    char line[8192];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            double x = 0, y = 0, z = 0;
            if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3) darr_push3(V, x, y, z);
        } else if (line[0] == 'f' && line[1] == ' ') {
            int idx[64]; int ni = 0;
            const char *p = line + 1;
            while (*p && ni < 64) {
                while (*p == ' ' || *p == '\t') p++;
                if (*p == '\0' || *p == '\n' || *p == '\r') break;
                int vi = (int)strtol(p, (char **)&p, 10);
                while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++; /* skip /vt/vn */
                if (vi < 0) vi = (int)(V->n / 3) + vi; else vi -= 1;
                idx[ni++] = vi;
            }
            for (int i = 1; i + 1 < ni; i++) iarr_push3(F, idx[0], idx[i], idx[i + 1]);
        }
    }
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ */
/* sorted-triple set for the PRE faces */
static void sort3(int *t) {        /* ascending in place */
    int a = t[0], b = t[1], c = t[2], tmp;
    if (a > b) { tmp = a; a = b; b = tmp; }
    if (b > c) { tmp = b; b = c; c = tmp; }
    if (a > b) { tmp = a; a = b; b = tmp; }
    t[0] = a; t[1] = b; t[2] = c;
}
static int tri_cmp(const void *pa, const void *pb) {
    const int *a = pa, *b = pb;
    if (a[0] != b[0]) return a[0] < b[0] ? -1 : 1;
    if (a[1] != b[1]) return a[1] < b[1] ? -1 : 1;
    if (a[2] != b[2]) return a[2] < b[2] ? -1 : 1;
    return 0;
}
/* bsearch over the sorted PRE triple table */
static int tri_in_set(const int *set, size_t nset, const int *key) {
    size_t lo = 0, hi = nset;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int c = tri_cmp(set + mid * 3, key);
        if (c == 0) return 1;
        if (c < 0) lo = mid + 1; else hi = mid;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
static int uf_find(int *uf, int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }

/* undirected edge, packed lo<<32|hi, for per-patch manifold/loop checks */
static int u64cmp(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}
static uint64_t ekey(int a, int b) {
    uint32_t lo = (uint32_t)(a < b ? a : b), hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | hi;
}

/* ------------------------------------------------------------------ */
typedef struct { int axis; double coord; } Seam;

typedef struct {
    int    id;                 /* 1-based, assigned after ranking */
    int    root;               /* uf root */
    long   nf;                 /* fill triangles */
    int    nrim;               /* rim (pre-existing boundary) verts in patch */
    int    nsteiner;           /* steiner verts in patch */
    double bbmin[3], bbmax[3], cen[3];
    int    n_loops;            /* boundary loops over the fill faces (1 = clean disk) */
    int    max_edge_mult;      /* max undirected-edge multiplicity (2 = manifold) */
    /* seam straddle (best plane) */
    int    seam_axis;          /* -1 if none */
    double seam_coord;
    int    n_side_lo, n_side_hi;   /* rim verts on each side, within band of the plane */
    double seam_gap;               /* coord span of band rim verts across the plane */
    int    n_planes_crossed;       /* how many candidate planes it straddles */
    double score;
} Patch;

static int patch_cmp_score(const void *a, const void *b) {
    const Patch *pa = a, *pb = b;
    if (pa->score != pb->score) return pa->score < pb->score ? 1 : -1; /* descending */
    return 0;
}

static const char *AXIS_NAME = "ZYX";

/* owning cube origin on an axis for a world coord (origins are multiples of 128) */
static long cube_origin(double coord) { return (long)floor(coord / 128.0) * 128; }

/* ================================================================== */
/* core analysis. returns patch array (caller frees), count in *out_n. */
static Patch *analyze(const DArr *Vp, const IArr *Fp /*post*/,
                      const int *preset, size_t npreset,
                      size_t nv_pre, const Seam *seams, int nseam,
                      double band, size_t *out_n, IArr *fill_face_idx /*opt*/) {
    const double *V = Vp->p;
    size_t nv = Vp->n / 3, nf = Fp->n / 3;
    const int *F = Fp->p;

    /* ---- 1. fill faces = POST faces not in PRE set ---- */
    int *uf = malloc(nv * sizeof(int));
    for (size_t i = 0; i < nv; i++) uf[i] = (int)i;
    size_t n_fill = 0;
    char *is_fill = calloc(nf, 1);
    for (size_t fi = 0; fi < nf; fi++) {
        int key[3] = { F[fi*3], F[fi*3+1], F[fi*3+2] };
        sort3(key);
        if (!tri_in_set(preset, npreset, key)) {
            is_fill[fi] = 1; n_fill++;
            int a = F[fi*3], b = F[fi*3+1], c = F[fi*3+2];
            int ra = uf_find(uf,a), rb = uf_find(uf,b); if (ra!=rb) uf[rb]=ra;
            ra = uf_find(uf,a); int rc = uf_find(uf,c); if (ra!=rc) uf[rc]=ra;
            if (fill_face_idx) { /* keep for caller */
                if (fill_face_idx->n + 1 > fill_face_idx->cap) {
                    fill_face_idx->cap = fill_face_idx->cap ? fill_face_idx->cap*2 : 1024;
                    fill_face_idx->p = realloc(fill_face_idx->p, fill_face_idx->cap*sizeof(int)); }
                fill_face_idx->p[fill_face_idx->n++] = (int)fi;
            }
        }
    }

    /* ---- 2. roots among fill faces -> patches ---- */
    /* count fill faces per root */
    long *rootnf = calloc(nv, sizeof(long));
    for (size_t fi = 0; fi < nf; fi++) if (is_fill[fi]) rootnf[uf_find(uf, F[fi*3])]++;
    int P = 0;
    for (size_t i = 0; i < nv; i++) if (rootnf[i] > 0) P++;
    Patch *patch = calloc((size_t)(P > 0 ? P : 1), sizeof(Patch));
    int *root2idx = malloc(nv * sizeof(int));
    for (size_t i = 0; i < nv; i++) root2idx[i] = -1;
    int pi = 0;
    for (size_t i = 0; i < nv; i++) if (rootnf[i] > 0) {
        root2idx[i] = pi;
        patch[pi].root = (int)i; patch[pi].nf = rootnf[i];
        patch[pi].seam_axis = -1;
        for (int k = 0; k < 3; k++) { patch[pi].bbmin[k] = 1e300; patch[pi].bbmax[k] = -1e300; }
        pi++;
    }

    /* ---- 3. per-patch vertex stats (rim/steiner, bbox, centroid) ---- */
    /* a vertex belongs to a patch if it is touched by a fill face of that patch */
    char *vseen = calloc(nv, 1);
    long *vcount = calloc((size_t)P ? (size_t)P : 1, sizeof(long)); /* verts per patch (for centroid avg) */
    for (size_t fi = 0; fi < nf; fi++) {
        if (!is_fill[fi]) continue;
        int p = root2idx[uf_find(uf, F[fi*3])];
        for (int k = 0; k < 3; k++) {
            int v = F[fi*3+k];
            if (vseen[v]) continue;        /* count each vert once per whole mesh; patches are disjoint */
            vseen[v] = 1;
            if ((size_t)v >= nv_pre) patch[p].nsteiner++; else patch[p].nrim++;
            const double *c = V + (size_t)v*3;
            for (int a = 0; a < 3; a++) {
                if (c[a] < patch[p].bbmin[a]) patch[p].bbmin[a] = c[a];
                if (c[a] > patch[p].bbmax[a]) patch[p].bbmax[a] = c[a];
                patch[p].cen[a] += c[a];
            }
            vcount[p]++;
        }
    }
    for (int p = 0; p < P; p++) if (vcount[p] > 0)
        for (int a = 0; a < 3; a++) patch[p].cen[a] /= (double)vcount[p];

    /* ---- 4. per-patch topology: boundary-loop count + edge multiplicity ---- */
    /* gather undirected edges of each patch's fill faces, sort, count runs */
    for (int p = 0; p < P; p++) {
        long m = patch[p].nf;
        uint64_t *ed = malloc((size_t)m * 3 * sizeof(uint64_t));
        size_t ne = 0;
        for (size_t fi = 0; fi < nf; fi++) {
            if (!is_fill[fi]) continue;
            if (root2idx[uf_find(uf, F[fi*3])] != p) continue;
            int a = F[fi*3], b = F[fi*3+1], c = F[fi*3+2];
            ed[ne++] = ekey(a,b); ed[ne++] = ekey(b,c); ed[ne++] = ekey(c,a);
        }
        qsort(ed, ne, sizeof(uint64_t), u64cmp);
        int maxmult = 0, nbound = 0;
        size_t i = 0;
        while (i < ne) {
            size_t j = i; while (j < ne && ed[j] == ed[i]) j++;
            int mult = (int)(j - i);
            if (mult > maxmult) maxmult = mult;
            if (mult == 1) nbound++;        /* boundary edge */
            i = j;
        }
        patch[p].max_edge_mult = maxmult;
        /* a clean disk has exactly nrim boundary edges forming ONE cycle; we report
         * loop count as boundary-edges / (avg per loop) -- but the robust, cheap proxy
         * is: nbound boundary edges, and for a single closed loop nbound == #distinct
         * boundary verts.  Track #loops via degree walk below. */
        patch[p].n_loops = nbound; /* placeholder, refined next */
        free(ed);
    }
    /* refine n_loops: count connected boundary cycles via uf over boundary edges */
    for (int p = 0; p < P; p++) {
        long m = patch[p].nf;
        uint64_t *ed = malloc((size_t)m * 3 * sizeof(uint64_t));
        size_t ne = 0;
        for (size_t fi = 0; fi < nf; fi++) {
            if (!is_fill[fi]) continue;
            if (root2idx[uf_find(uf, F[fi*3])] != p) continue;
            int a = F[fi*3], b = F[fi*3+1], c = F[fi*3+2];
            ed[ne++] = ekey(a,b); ed[ne++] = ekey(b,c); ed[ne++] = ekey(c,a);
        }
        qsort(ed, ne, sizeof(uint64_t), u64cmp);
        /* boundary edges only */
        IArr be = {0};
        size_t i = 0;
        while (i < ne) {
            size_t j = i; while (j < ne && ed[j] == ed[i]) j++;
            if (j - i == 1) {
                int lo = (int)(ed[i] >> 32), hi = (int)(ed[i] & 0xffffffffu);
                iarr_push3(&be, lo, hi, 0);
            }
            i = j;
        }
        /* count boundary cycles: uf over boundary-edge endpoints, keyed by raw
         * vertex id over the whole mesh (sparse, but cheap and unambiguous) */
        if (be.n == 0) { patch[p].n_loops = 0; free(ed); free(be.p); continue; }
        size_t neb = be.n / 3;
        int *luf = malloc(nv * sizeof(int));
        for (size_t k = 0; k < nv; k++) luf[k] = (int)k;
        for (size_t e = 0; e < neb; e++) {
            int a = be.p[e*3], b = be.p[e*3+1];
            int ra = uf_find(luf,a), rb = uf_find(luf,b); if (ra!=rb) luf[rb]=ra;
        }
        /* count distinct roots among boundary verts */
        int loops = 0;
        /* mark roots seen */
        char *seenr = calloc(nv,1);
        for (size_t e = 0; e < neb; e++) {
            int r = uf_find(luf, be.p[e*3]);
            if (!seenr[r]) { seenr[r] = 1; loops++; }
        }
        patch[p].n_loops = loops;
        free(seenr); free(luf); free(be.p); free(ed);
    }

    /* ---- 5. seam straddle per patch ---- */
    for (int p = 0; p < P; p++) {
        int crossed = 0, best = -1; double best_balance = -1;
        for (int s = 0; s < nseam; s++) {
            int a = seams[s].axis; double c = seams[s].coord;
            if (!(patch[p].bbmin[a] < c - 0.05 && patch[p].bbmax[a] > c + 0.05)) continue;
            crossed++;
            /* rim verts within band on each side */
            int lo = 0, hi = 0;
            for (size_t fi = 0; fi < nf; fi++) {
                if (!is_fill[fi]) continue;
                if (root2idx[uf_find(uf, F[fi*3])] != p) continue;
                for (int k = 0; k < 3; k++) {
                    int v = F[fi*3+k];
                    if ((size_t)v >= nv_pre) continue; /* rim only */
                    double d = V[(size_t)v*3 + a] - c;
                    if (d < 0 && -d <= band) lo++;
                    else if (d > 0 && d <= band) hi++;
                }
            }
            int bal = lo < hi ? lo : hi;   /* balance: both sides represented */
            if (bal > best_balance) { best_balance = bal; best = s;
                patch[p].n_side_lo = lo; patch[p].n_side_hi = hi;
                patch[p].seam_axis = a; patch[p].seam_coord = c;
                patch[p].seam_gap = patch[p].bbmax[a] - patch[p].bbmin[a];
            }
        }
        patch[p].n_planes_crossed = crossed;
    }

    /* ---- 6. score: clean single-seam disk that genuinely bridges the gap ---- */
    for (int p = 0; p < P; p++) {
        Patch *q = &patch[p];
        double sc = 0;
        if (q->seam_axis >= 0 && q->n_side_lo > 0 && q->n_side_hi > 0) {
            int a = q->seam_axis;
            double ext_n = q->bbmax[a] - q->bbmin[a];               /* thickness across seam */
            double ext_o = 0;                                       /* span along seam */
            for (int k = 0; k < 3; k++) if (k != a) {
                double e = q->bbmax[k] - q->bbmin[k]; if (e > ext_o) ext_o = e;
            }
            int bal = q->n_side_lo < q->n_side_hi ? q->n_side_lo : q->n_side_hi;
            sc  = 100.0;
            sc += 8.0 * bal;                       /* well-anchored on both cubes */
            sc += (q->n_loops == 1) ? 40.0 : -30.0 * (q->n_loops - 1); /* clean disk */
            sc += (q->max_edge_mult <= 2) ? 25.0 : -40.0;             /* manifold */
            sc += (q->n_planes_crossed == 1) ? 20.0 : -15.0;          /* not a corner */
            sc -= 6.0 * ext_n;                     /* thin crack, not a gaping hole */
            sc += 3.0 * (ext_o > ext_n ? 1.0 : 0.0);
            if (q->nrim < 4) sc -= 50.0;
            if (q->nrim > 200) sc -= 0.2 * (q->nrim - 200);  /* mega-holes are not "clean" */
        }
        q->score = sc;
    }

    free(uf); free(is_fill); free(rootnf); free(root2idx);
    free(vseen); free(vcount);
    *out_n = (size_t)P;
    return patch;
}

/* ------------------------------------------------------------------ */
/* dump one patch (red) + surrounding POST faces within radius (gray) */
static int dump_patch(const char *out, const DArr *Vp, const IArr *Fp,
                      const int *preset, size_t npreset, size_t nv_pre,
                      const Seam *seams, int nseam, double band,
                      int want_id, double radius) {
    size_t np = 0; IArr fillfaces = {0};
    Patch *patch = analyze(Vp, Fp, preset, npreset, nv_pre, seams, nseam, band, &np, NULL);
    /* rank to map want_id -> root (same order as printed) */
    qsort(patch, np, sizeof(Patch), patch_cmp_score);
    if (want_id < 1 || (size_t)want_id > np) {
        fprintf(stderr, "dump: patch id %d out of range (1..%zu)\n", want_id, np);
        free(patch); free(fillfaces.p); return 1;
    }
    Patch *q = &patch[want_id - 1];
    double cen[3] = { q->cen[0], q->cen[1], q->cen[2] };

    /* recompute is_fill + roots to identify this patch's faces (re-run uf cheaply) */
    const double *V = Vp->p; size_t nv = Vp->n/3, nf = Fp->n/3; const int *F = Fp->p;
    int *uf = malloc(nv*sizeof(int)); for (size_t i=0;i<nv;i++) uf[i]=(int)i;
    char *is_fill = calloc(nf,1);
    for (size_t fi=0; fi<nf; fi++){ int key[3]={F[fi*3],F[fi*3+1],F[fi*3+2]}; sort3(key);
        if(!tri_in_set(preset,npreset,key)){ is_fill[fi]=1;
            int a=F[fi*3],b=F[fi*3+1],c=F[fi*3+2];
            int ra=uf_find(uf,a),rb=uf_find(uf,b); if(ra!=rb)uf[rb]=ra;
            ra=uf_find(uf,a); int rc=uf_find(uf,c); if(ra!=rc)uf[rc]=ra; } }
    int target_root = uf_find(uf, q->root);

    /* select faces: this patch's fill faces, plus any POST face fully within radius */
    double r2 = radius * radius;
    char *vsel  = calloc(nv, 1);   /* selected into the crop */
    char *vfill = calloc(nv, 1);   /* touched by a fill triangle of THIS patch -> red */
    IArr keepf = {0};
    char *fis_fill = calloc(nf,1);
    for (size_t fi=0; fi<nf; fi++) {
        int a=F[fi*3],b=F[fi*3+1],c=F[fi*3+2];
        int patchface = is_fill[fi] && uf_find(uf,a)==target_root;
        int near = 0;
        if (!patchface) {
            int allnear = 1;
            for (int k=0;k<3;k++){ const double *cc=V+(size_t)F[fi*3+k]*3;
                double dz=cc[0]-cen[0],dy=cc[1]-cen[1],dx=cc[2]-cen[2];
                if (dz*dz+dy*dy+dx*dx > r2){ allnear=0; break; } }
            near = allnear;
        }
        if (patchface || near) {
            iarr_push3(&keepf, a, b, c);
            fis_fill[keepf.n/3 - 1] = (char)(patchface ? 1 : 0);
            vsel[a]=vsel[b]=vsel[c]=1;
            if (patchface) { vfill[a]=vfill[b]=vfill[c]=1; }
        }
    }
    /* compact verts */
    int *remap = malloc(nv*sizeof(int)); for (size_t i=0;i<nv;i++) remap[i]=-1;
    FILE *fo = fopen(out, "w");
    if (!fo) { fprintf(stderr,"dump: cannot write %s\n", out); free(patch); return 1; }
    char seamtag[24];
    if (q->seam_axis >= 0) snprintf(seamtag, sizeof seamtag, "%c=%.1f", AXIS_NAME[q->seam_axis], q->seam_coord);
    else snprintf(seamtag, sizeof seamtag, "interior");
    fprintf(fo, "# crossseam_fill dump: patch #%d  seam %s  centroid(Z,Y,X)=%.1f,%.1f,%.1f\n",
            want_id, seamtag, cen[0], cen[1], cen[2]);
    fprintf(fo, "# red = hole-fill triangles across the seam ; gray = surrounding welded mesh\n");
    int next = 1;
    for (size_t i=0;i<nv;i++) if (vsel[i]) {
        remap[i] = next++;
        const double *cc = V + i*3;
        /* red = vertex on the hole rim or a Steiner point of this fill; gray = context */
        int isfillv = vfill[i];
        fprintf(fo, "v %.4f %.4f %.4f %.3f %.3f %.3f\n", cc[0], cc[1], cc[2],
                isfillv ? 0.95 : 0.62, isfillv ? 0.15 : 0.62, isfillv ? 0.15 : 0.62);
    }
    /* faces, grouped */
    fprintf(fo, "g context\n");
    for (size_t fi=0; fi<keepf.n/3; fi++) if (!fis_fill[fi])
        fprintf(fo, "f %d %d %d\n", remap[keepf.p[fi*3]], remap[keepf.p[fi*3+1]], remap[keepf.p[fi*3+2]]);
    fprintf(fo, "g fill\n");
    for (size_t fi=0; fi<keepf.n/3; fi++) if (fis_fill[fi])
        fprintf(fo, "f %d %d %d\n", remap[keepf.p[fi*3]], remap[keepf.p[fi*3+1]], remap[keepf.p[fi*3+2]]);
    fclose(fo);

    long nctx=0, nfil=0; for (size_t fi=0; fi<keepf.n/3; fi++) { if (fis_fill[fi]) nfil++; else nctx++; }
    printf("dumped patch #%d -> %s  (%ld fill tris red, %ld context tris gray, %d verts)\n",
           want_id, out, nfil, nctx, next-1);

    free(patch); free(uf); free(is_fill); free(vsel); free(vfill); free(keepf.p);
    free(fis_fill); free(remap); free(fillfaces.p);
    return 0;
}

/* ================================================================== */
static void print_patch_line(const Patch *q, int rank) {
    char seam[32] = "interior";
    if (q->seam_axis >= 0)
        snprintf(seam, sizeof seam, "%c=%.0f", AXIS_NAME[q->seam_axis], q->seam_coord);
    char topo[48];
    snprintf(topo, sizeof topo, "%s%s",
             q->n_loops == 1 ? "disk" : "MULTI-LOOP",
             q->max_edge_mult > 2 ? ",NONMANIFOLD" : "");
    printf("  #%-3d score=%-6.0f %-9s rim=%-4d steiner=%-3d tris=%-4ld %-18s "
           "sides(%d|%d) thick=%.2f  @Z,Y,X=%.0f,%.0f,%.0f\n",
           rank, q->score, seam, q->nrim, q->nsteiner, q->nf, topo,
           q->n_side_lo, q->n_side_hi, q->seam_gap,
           q->cen[0], q->cen[1], q->cen[2]);
}

/* ------------------------------------------------------------------ */
static int run_selftest(void) {
    /* PRE: two separated strips either side of plane X=10, each an open quad rim.
     *   left strip  x in {8,9}, right strip x in {11,12}, z,y span {0,1}
     * POST: PRE + two fill triangles bridging the 9->11 gap into one quad. */
    DArr Vpre = {0}, Vpost = {0}; IArr Fpre = {0}, Fpost = {0};
    /* verts (Z,Y,X): index order matters */
    /* 0: 0,0,8  1: 1,0,8  2: 0,0,9  3: 1,0,9   (left strip)
       4: 0,0,11 5: 1,0,11 6: 0,0,12 7: 1,0,12  (right strip) */
    double pts[8][3] = {
        {0,0,8},{1,0,8},{0,0,9},{1,0,9},
        {0,0,11},{1,0,11},{0,0,12},{1,0,12}
    };
    for (int i=0;i<8;i++){ darr_push3(&Vpre, pts[i][0],pts[i][1],pts[i][2]);
                           darr_push3(&Vpost,pts[i][0],pts[i][1],pts[i][2]); }
    /* left quad (0,1,3,2) and right quad (4,5,7,6) as two tris each -> PRE */
    iarr_push3(&Fpre,0,1,3); iarr_push3(&Fpre,0,3,2);
    iarr_push3(&Fpre,4,5,7); iarr_push3(&Fpre,4,7,6);
    /* POST keeps those (reordered to mimic holefill) + bridge quad (2,3,5,4) */
    iarr_push3(&Fpost,4,5,7); iarr_push3(&Fpost,0,1,3);
    iarr_push3(&Fpost,4,7,6); iarr_push3(&Fpost,0,3,2);
    iarr_push3(&Fpost,2,3,5); iarr_push3(&Fpost,2,5,4);    /* the cross-seam fill */

    /* build PRE set */
    size_t npre = Fpre.n/3;
    int *preset = malloc(npre*3*sizeof(int));
    for (size_t i=0;i<npre;i++){ int t[3]={Fpre.p[i*3],Fpre.p[i*3+1],Fpre.p[i*3+2]};
        sort3(t); preset[i*3]=t[0];preset[i*3+1]=t[1];preset[i*3+2]=t[2]; }
    qsort(preset, npre, 3*sizeof(int), tri_cmp);

    Seam seams[1] = { { 2, 10.0 } };   /* axis 2 = X, plane X=10 */
    size_t np = 0;
    Patch *patch = analyze(&Vpost,&Fpost, preset,npre, 8 /*nv_pre*/, seams,1, 4.0, &np, NULL);

    int ok = 1;
    if (np != 1) { fprintf(stderr,"selftest FAIL: expected 1 patch, got %zu\n", np); ok=0; }
    if (ok && patch[0].nf != 2) { fprintf(stderr,"selftest FAIL: expected 2 fill tris, got %ld\n", patch[0].nf); ok=0; }
    if (ok && patch[0].seam_axis != 2) { fprintf(stderr,"selftest FAIL: seam axis %d != 2(X)\n", patch[0].seam_axis); ok=0; }
    if (ok && !(patch[0].n_side_lo>0 && patch[0].n_side_hi>0)) { fprintf(stderr,"selftest FAIL: not straddling (lo=%d hi=%d)\n", patch[0].n_side_lo, patch[0].n_side_hi); ok=0; }
    if (ok && patch[0].nrim != 4) { fprintf(stderr,"selftest FAIL: expected rim=4, got %d\n", patch[0].nrim); ok=0; }
    if (ok && patch[0].n_loops != 1) { fprintf(stderr,"selftest FAIL: expected 1 loop, got %d\n", patch[0].n_loops); ok=0; }
    if (ok && patch[0].max_edge_mult != 2) { fprintf(stderr,"selftest FAIL: expected manifold (mult=2), got %d\n", patch[0].max_edge_mult); ok=0; }

    /* second test: a purely interior fill (no seam) must NOT be flagged cross-seam */
    {
        DArr V2={0}; IArr Fpre2={0},Fpost2={0};
        double q[4][3]={{0,0,0},{1,0,0},{0,1,0},{1,1,0}};
        for(int i=0;i<4;i++) darr_push3(&V2,q[i][0],q[i][1],q[i][2]);
        /* PRE: one tri; POST: + a second tri completing the quad, far from X=10 */
        iarr_push3(&Fpre2,0,1,3);
        iarr_push3(&Fpost2,0,1,3); iarr_push3(&Fpost2,0,3,2);
        size_t npre2=Fpre2.n/3; int *ps2=malloc(npre2*3*sizeof(int));
        for(size_t i=0;i<npre2;i++){int t[3]={Fpre2.p[i*3],Fpre2.p[i*3+1],Fpre2.p[i*3+2]};sort3(t);
            ps2[i*3]=t[0];ps2[i*3+1]=t[1];ps2[i*3+2]=t[2];}
        qsort(ps2,npre2,3*sizeof(int),tri_cmp);
        size_t np2=0; Patch *pt2=analyze(&V2,&Fpost2,ps2,npre2,4,seams,1,4.0,&np2,NULL);
        if (np2!=1){fprintf(stderr,"selftest FAIL(2): expected 1 patch got %zu\n",np2);ok=0;}
        if (ok && pt2[0].seam_axis!=-1){fprintf(stderr,"selftest FAIL(2): interior fill wrongly tagged seam axis %d\n",pt2[0].seam_axis);ok=0;}
        free(pt2);free(ps2);free(V2.p);free(Fpre2.p);free(Fpost2.p);
    }

    free(patch); free(preset);
    free(Vpre.p);free(Vpost.p);free(Fpre.p);free(Fpost.p);
    if (ok) printf("selftest OK\n");
    return ok ? 0 : 1;
}

/* ================================================================== */
int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return run_selftest();
    if (argc < 3) {
        fprintf(stderr,
          "usage: %s <pre.obj> <post.obj> [--seam Z=4480 ...] [--band v] [--min n]\n"
          "          [--max n] [--top k] [--all] [--check] [--dump out.obj --patch id [--radius v]]\n"
          "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    const char *pre = argv[1], *post = argv[2];
    Seam seams[64]; int nseam = 0;
    double band = 4.0, radius = 14.0;
    int minrim = 4, maxrim = 400, top = 20, want_all = 0, want_check = 0;
    const char *dumpout = NULL; int dump_id = 0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i],"--seam") && i+1<argc) {
            char ax; double cv;
            if (sscanf(argv[++i], "%c=%lf", &ax, &cv)==2) {
                int a = ax=='Z'||ax=='z'?0 : ax=='Y'||ax=='y'?1 : ax=='X'||ax=='x'?2 : -1;
                if (a<0){ fprintf(stderr,"bad --seam axis '%c'\n",ax); return 2; }
                if (nseam<64){ seams[nseam].axis=a; seams[nseam].coord=cv; nseam++; }
            }
        }
        else if (!strcmp(argv[i],"--band")  && i+1<argc) band   = atof(argv[++i]);
        else if (!strcmp(argv[i],"--radius")&& i+1<argc) radius = atof(argv[++i]);
        else if (!strcmp(argv[i],"--min")   && i+1<argc) minrim = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--max")   && i+1<argc) maxrim = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--top")   && i+1<argc) top    = atoi(argv[++i]);
        else if (!strcmp(argv[i],"--dump")  && i+1<argc) dumpout= argv[++i];
        else if (!strcmp(argv[i],"--patch") && i+1<argc) dump_id= atoi(argv[++i]);
        else if (!strcmp(argv[i],"--all")) want_all = 1;
        else if (!strcmp(argv[i],"--check")) want_check = 1;
        else { fprintf(stderr,"unknown arg: %s\n", argv[i]); return 2; }
    }

    DArr Vpre={0},Vpost={0}; IArr Fpre={0},Fpost={0};
    if (load_obj(pre,&Vpre,&Fpre)) return 1;
    if (load_obj(post,&Vpost,&Fpost)) return 1;
    size_t nv_pre = Vpre.n/3, nv_post = Vpost.n/3;
    printf("=== crossseam_fill ===\n");
    printf("  PRE  %s : %zu v, %zu f\n", pre, nv_pre, Fpre.n/3);
    printf("  POST %s : %zu v, %zu f  (+%zu v, +%zu f)\n", post, nv_post, Fpost.n/3,
           nv_post - nv_pre, Fpost.n/3 - Fpre.n/3);

    /* auto-infer seams = multiples of 128 strictly inside the POST bbox, per axis */
    if (nseam == 0) {
        double mn[3]={1e300,1e300,1e300}, mx[3]={-1e300,-1e300,-1e300};
        for (size_t i=0;i<nv_post;i++) for(int a=0;a<3;a++){
            double c=Vpost.p[i*3+a]; if(c<mn[a])mn[a]=c; if(c>mx[a])mx[a]=c; }
        for (int a=0;a<3;a++) {
            long k0 = (long)ceil((mn[a]+1.0)/128.0), k1=(long)floor((mx[a]-1.0)/128.0);
            for (long k=k0;k<=k1;k++) if (nseam<64){ seams[nseam].axis=a; seams[nseam].coord=(double)(k*128); nseam++; }
        }
        printf("  seams (auto): ");
        for (int s=0;s<nseam;s++) printf("%c=%.0f ", AXIS_NAME[seams[s].axis], seams[s].coord);
        printf("\n");
    }

    /* build PRE sorted-triple set */
    size_t npre = Fpre.n/3;
    int *preset = malloc((npre?npre:1)*3*sizeof(int));
    for (size_t i=0;i<npre;i++){ int t[3]={Fpre.p[i*3],Fpre.p[i*3+1],Fpre.p[i*3+2]};
        sort3(t); preset[i*3]=t[0];preset[i*3+1]=t[1];preset[i*3+2]=t[2]; }
    qsort(preset, npre, 3*sizeof(int), tri_cmp);

    size_t np=0; IArr fillidx={0};
    Patch *patch = analyze(&Vpost,&Fpost, preset,npre, nv_pre, seams,nseam, band, &np, &fillidx);
    qsort(patch, np, sizeof(Patch), patch_cmp_score);

    long cross=0; for (size_t i=0;i<np;i++) if (patch[i].seam_axis>=0 && patch[i].n_side_lo>0 && patch[i].n_side_hi>0) cross++;
    printf("  fill faces: %zu ; fill patches (holes filled): %zu ; cross-seam: %ld\n\n",
           fillidx.n, np, cross);

    if (want_check) {
        const double *V = Vpost.p; const int *F = Fpost.p;
        double min_area = 1e300, min_edge = 1e300; long ndegen = 0, nshort = 0;
        long deg_rr = 0, deg_rs = 0, deg_ss = 0;   /* degenerate-tri shortest-edge class */
        for (size_t t = 0; t < fillidx.n; t++) {
            size_t fi = (size_t)fillidx.p[t];
            int ia = F[fi*3+0], ib = F[fi*3+1], ic = F[fi*3+2];
            const double *A = V + (size_t)ia*3, *B = V + (size_t)ib*3, *C = V + (size_t)ic*3;
            double e1[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
            double e2[3] = { C[0]-A[0], C[1]-A[1], C[2]-A[2] };
            double cx = e1[1]*e2[2]-e1[2]*e2[1], cy = e1[2]*e2[0]-e1[0]*e2[2], cz = e1[0]*e2[1]-e1[1]*e2[0];
            double area = 0.5 * sqrt(cx*cx+cy*cy+cz*cz);
            if (area < min_area) min_area = area;
            double l1 = sqrt((B[0]-A[0])*(B[0]-A[0])+(B[1]-A[1])*(B[1]-A[1])+(B[2]-A[2])*(B[2]-A[2]));
            double l2 = sqrt((C[0]-B[0])*(C[0]-B[0])+(C[1]-B[1])*(C[1]-B[1])+(C[2]-B[2])*(C[2]-B[2]));
            double l3 = sqrt((A[0]-C[0])*(A[0]-C[0])+(A[1]-C[1])*(A[1]-C[1])+(A[2]-C[2])*(A[2]-C[2]));
            double mn = l1 < l2 ? (l1 < l3 ? l1 : l3) : (l2 < l3 ? l2 : l3);
            if (mn < min_edge) min_edge = mn;
            if (mn < 0.1) nshort++;
            if (area < 1e-4) {
                ndegen++;
                /* classify the shortest edge's endpoints: rim (<nv_pre) vs Steiner */
                int ea, eb;
                if (l1 <= l2 && l1 <= l3) { ea = ia; eb = ib; }
                else if (l2 <= l3)        { ea = ib; eb = ic; }
                else                      { ea = ic; eb = ia; }
                int sa = (size_t)ea >= nv_pre, sb = (size_t)eb >= nv_pre;
                if (sa && sb) deg_ss++; else if (sa || sb) deg_rs++; else deg_rr++;
            }
        }
        printf("  [check] fill tris=%zu  degenerate(area<1e-4)=%ld  short-edge(<0.1vox)=%ld  "
               "min_area=%.5f  min_edge=%.4f\n", fillidx.n, ndegen, nshort, min_area, min_edge);
        printf("  [check] degenerate by shortest-edge class: rim-rim=%ld rim-steiner=%ld steiner-steiner=%ld\n\n",
               deg_rr, deg_rs, deg_ss);
    }

    printf("ranked %s fill patches (rim in [%d,%d]):\n",
           want_all ? "ALL" : "cross-seam", minrim, maxrim);
    int shown=0;
    for (size_t i=0;i<np && shown<top;i++) {
        Patch *q=&patch[i];
        int is_cross = q->seam_axis>=0 && q->n_side_lo>0 && q->n_side_hi>0;
        if (!want_all && !is_cross) continue;
        if (q->nrim<minrim || q->nrim>maxrim) continue;
        print_patch_line(q, (int)i+1);
        shown++;
    }
    if (shown==0) printf("  (none matched)\n");

    /* name the two cubes the top cross-seam patch divides */
    for (size_t i=0;i<np;i++){ Patch *q=&patch[i];
        if (q->seam_axis>=0 && q->n_side_lo>0 && q->n_side_hi>0){
            int a=q->seam_axis;
            long o[3]; for(int k=0;k<3;k++) o[k]=cube_origin(q->cen[k]);
            long lo_o[3]={o[0],o[1],o[2]}, hi_o[3]={o[0],o[1],o[2]};
            lo_o[a]=cube_origin(q->seam_coord-1.0); hi_o[a]=cube_origin(q->seam_coord+1.0);
            printf("\nbest cross-seam fill = patch #%zu : seam %c=%.0f divides cubes\n", i+1,
                   AXIS_NAME[a], q->seam_coord);
            printf("    z%05ld_y%05ld_x%05ld  <->  z%05ld_y%05ld_x%05ld\n",
                   lo_o[0],lo_o[1],lo_o[2], hi_o[0],hi_o[1],hi_o[2]);
            printf("    %d-vertex rim, %ld fill tris (+%d Steiner), %s, seam gap %.2f vox\n",
                   q->nrim, q->nf, q->nsteiner,
                   q->n_loops==1 ? (q->max_edge_mult<=2?"clean manifold disk":"single loop, non-manifold")
                                 : "multi-loop", q->seam_gap);
            break;
        }
    }

    if (dumpout) {
        if (dump_id<1) { fprintf(stderr,"--dump needs --patch <id>\n"); }
        else dump_patch(dumpout, &Vpost,&Fpost, preset,npre, nv_pre, seams,nseam, band, dump_id, radius);
    }

    free(patch); free(preset); free(fillidx.p);
    free(Vpre.p);free(Vpost.p);free(Fpre.p);free(Fpost.p);
    return 0;
}
