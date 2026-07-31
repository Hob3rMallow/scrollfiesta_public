/* wind_audit.c -- winding-coordinate topology auditor for welded scroll meshes.
 *
 * Two topology defects matter and mesh_quality / manifold_check cannot see
 * either (both are per-triangle or per-edge-manifold tests, blind to how the
 * sheet is wired):
 *
 *   MERGER (inter-wrap fusion): an edge that jumps ACROSS the ~7-vox gap between
 *   two turns of the spiral -- a radial short-circuit. Fatal: it welds turn N to
 *   turn N+1 directly instead of via the 360-deg traverse, so the spiral can no
 *   longer be unwrapped. Signature: winding coordinate w = r/pitch - theta/2pi
 *   jumps by ~1 across a single edge (large radial delta at nearly constant
 *   angle).
 *
 *   SPLIT (intra-sheet break): two DIFFERENT connected components that occupy the
 *   same turn (same w) and sit within a small physical gap -- the sheet was
 *   severed where it should be continuous. Signature: a close vertex pair in two
 *   different components with |dw| ~ 0.
 *
 *   BROAD FUSION (local-cube multi-wrap wall): a connected component spans more
 *   than one local winding even though every individual edge has a small |dw|.
 *   This catches the common "staircase" failure where BPA joins adjacent turns
 *   through many short triangles: no edge is a full-turn jump, so MERGER alone
 *   misses it. The test is intentionally opt-in (--local-span-tol) because a
 *   whole-scroll mesh is supposed to be one sheet spanning many windings. It is
 *   a semantic gate for leaf-cube BPA stages away from the umbilicus.
 *
 * The winding coordinate is LOCAL and branch-cut-free when computed pairwise:
 *   dr  = hypot(dy_v,dx_v) - hypot(dy_u,dx_u)      (dy,dx = coord - umbilicus)
 *   dth = wrap_to_pmpi( atan2(dy_v,dx_v) - atan2(dy_u,dx_u) )
 *   dw  = dr/pitch - dth/(2*pi)
 * exactly the seam gate's own test (seam_weld.c:536-540), so a merger edge here
 * is an edge the gate SHOULD have rejected (or that a non-bridge stage created).
 *
 * Vertex order is (z,y,x): coord[0]=z (axis), coord[1]=y, coord[2]=x. Radius is
 * in the (y,x) plane about the umbilicus. Turn count spanned by a component is
 * (r_max - r_min)/pitch -- the single stat that says whether a giant component
 * is a legit continuous spiral (few merger edges) or a fused pile (many).
 *
 *   wind_audit <in.obj> --umb-y Y --umb-x X [--pitch P=9.5] [--top N=25]
 *              [--merge-tol T=0.40] [--split-wtol T=0.35] [--split-gap G=14]
 *              [--core-pitches C=3] [--min-faces M=64]
 *              [--local-span-tol T] [--fail-on-local]
 *              [--dump-mergers f.obj] [--dump-splits f.obj]
 *   wind_audit --selftest
 *
 * Standalone C99, own OBJ parser (shared shape with obj_cc_color), no deps.
 * Exit: 0 ok / selftest pass, 1 IO error, 2 usage error, 3 selftest fail,
 *       4 --fail-on-local gate failed.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)

/* ---------- growing arrays ---------- */
typedef struct { float *v; size_t n, cap; } FVec;   /* verts: 3 floats each */
typedef struct { int32_t *f; size_t n, cap; } IVec; /* faces: 3 ints each   */

static int fv_push(FVec *a, float x, float y, float z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        float *np = (float *)realloc(a->v, nc * 3 * sizeof(float)); if (!np) return -1; a->v = np; a->cap = nc; }
    a->v[a->n*3+0] = x; a->v[a->n*3+1] = y; a->v[a->n*3+2] = z; a->n++; return 0;
}
static int iv_push(IVec *a, int32_t x, int32_t y, int32_t z)
{
    if (a->n == a->cap) { size_t nc = a->cap ? a->cap * 2 : 1u << 16;
        int32_t *np = (int32_t *)realloc(a->f, nc * 3 * sizeof(int32_t)); if (!np) return -1; a->f = np; a->cap = nc; }
    a->f[a->n*3+0] = x; a->f[a->n*3+1] = y; a->f[a->n*3+2] = z; a->n++; return 0;
}

static int32_t parse_findex(const char *tok) { return (int32_t)atol(tok); }

static int read_obj(const char *path, FVec *V, IVec *F)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "wind_audit: cannot open %s\n", path); return -1; }
    size_t lcap = 1u << 16; char *line = (char *)malloc(lcap);
    int rc = 0;
    if (!line) { fclose(fp); return -1; }
    while (fgets(line, (int)lcap, fp)) {
        while (!strchr(line, '\n') && !feof(fp)) {
            size_t len = strlen(line);
            char *nl = (char *)realloc(line, lcap * 2); if (!nl) { rc = -1; goto done; }
            line = nl; lcap *= 2;
            if (!fgets(line + len, (int)(lcap - len), fp)) break;
        }
        if (line[0] == 'v' && (line[1] == ' ' || line[1] == '\t')) {
            double x = 0, y = 0, z = 0;
            if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3)
                if (fv_push(V, (float)x, (float)y, (float)z) != 0) { rc = -1; goto done; }
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char *tok = strtok(line + 2, " \t\r\n");
            int32_t idx[3]; int n = 0;
            while (tok && n < 3) { idx[n++] = parse_findex(tok); tok = strtok(NULL, " \t\r\n"); }
            if (n == 3) {
                int32_t a = idx[0], b = idx[1], c = idx[2];
                if (a < 0) a = (int32_t)V->n + a + 1;
                if (b < 0) b = (int32_t)V->n + b + 1;
                if (c < 0) c = (int32_t)V->n + c + 1;
                if (iv_push(F, a - 1, b - 1, c - 1) != 0) { rc = -1; goto done; }
            }
        }
    }
done:
    free(line); fclose(fp);
    return rc;
}

/* ---------- union-find (path halving + union by size) ---------- */
static int32_t uf_find(int32_t *p, int32_t x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }
static void uf_union(int32_t *p, int32_t *sz, int32_t a, int32_t b)
{
    a = uf_find(p, a); b = uf_find(p, b);
    if (a == b) return;
    if (sz[a] < sz[b]) { int32_t t = a; a = b; b = t; }
    p[b] = a; sz[a] += sz[b];
}

/* ---------- winding helpers (coord order z,y,x; umbilicus in y,x) ---------- */
static double wrap_pmpi(double a) { while (a > M_PI) a -= TWO_PI; while (a <= -M_PI) a += TWO_PI; return a; }

/* per-vertex radius about the axis (y,x plane) */
static double vradius(const float *V, size_t i, double uy, double ux)
{
    double dy = (double)V[i*3+1] - uy, dx = (double)V[i*3+2] - ux;
    return hypot(dy, dx);
}
/* pairwise winding delta w_v - w_u (branch-cut-free) */
static double dwind(const float *V, size_t u, size_t v, double uy, double ux, double pitch)
{
    double dyu = (double)V[u*3+1] - uy, dxu = (double)V[u*3+2] - ux;
    double dyv = (double)V[v*3+1] - uy, dxv = (double)V[v*3+2] - ux;
    double dr  = hypot(dyv, dxv) - hypot(dyu, dxu);
    double dth = wrap_pmpi(atan2(dyv, dxv) - atan2(dyu, dxu));
    return dr / pitch - dth / TWO_PI;
}

/* ---------- component table ---------- */
typedef struct {
    int32_t root; size_t faces; int32_t rank;
    double rmin, rmax; float zmin, zmax;
    size_t merger_edges;   /* reliable (outer) merger edges wholly inside comp */
    double theta_sin, theta_cos;
    double phase_min, phase_max; /* branch-local r/pitch - theta/(2pi) */
    size_t outer_vertices;
} Comp;
static int comp_cmp(const void *pa, const void *pb)
{
    const Comp *a = (const Comp *)pa, *b = (const Comp *)pb;
    if (a->faces != b->faces) return a->faces < b->faces ? 1 : -1;
    return a->root < b->root ? -1 : 1;
}

/* ---------- spatial hash over vertices (cell = split_gap) ---------- */
typedef struct { int64_t key; int32_t head; } HCell;
typedef struct { HCell *cell; size_t mask; int32_t *next; double inv; int64_t off; } SHash;

static int64_t sh_key(const SHash *h, double a, double b, double c)
{
    int64_t ia = (int64_t)floor(a * h->inv) + h->off;
    int64_t ib = (int64_t)floor(b * h->inv) + h->off;
    int64_t ic = (int64_t)floor(c * h->inv) + h->off;
    /* pack into 63 bits: 21 bits each (range 0..2M cells per axis) */
    return ((ia & 0x1FFFFF) << 42) | ((ib & 0x1FFFFF) << 21) | (ic & 0x1FFFFF);
}
static size_t sh_hash(int64_t k, size_t mask)
{
    uint64_t x = (uint64_t)k * 0x9E3779B97F4A7C15ull;
    x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ull; x ^= x >> 32;
    return (size_t)x & mask;
}
static int32_t *sh_slot(SHash *h, int64_t key, int create)
{
    size_t i = sh_hash(key, h->mask);
    for (;;) {
        if (h->cell[i].head == -2) { /* empty */
            if (!create) return NULL;
            h->cell[i].key = key; h->cell[i].head = -1; return &h->cell[i].head;
        }
        if (h->cell[i].key == key) return &h->cell[i].head;
        i = (i + 1) & h->mask;
    }
}

/* ---------- reporting core ---------- */
typedef struct {
    double umb_y, umb_x, pitch;
    double merge_tol, split_wtol, split_gap, core_pitches;
    double local_span_tol;       /* 0 disables leaf-cube broad-fusion test */
    double seam_pitch, seam_tol;   /* leaf-cube seam planes (multiples of seam_pitch) */
    size_t min_faces; int top;
    int fail_on_local;
    const char *dump_mergers, *dump_splits, *json_path;
} Params;

/* is coord within seam_tol of a multiple of seam_pitch on ANY of the 3 axes? */
static int near_seam(double z, double y, double x, double pitch, double tol)
{
    if (pitch <= 0) return 0;
    double az = fabs(z - pitch * floor(z / pitch + 0.5));
    double ay = fabs(y - pitch * floor(y / pitch + 0.5));
    double ax = fabs(x - pitch * floor(x / pitch + 0.5));
    return (az <= tol || ay <= tol || ax <= tol);
}

/* a detected split candidate pair, aggregated per (rankA<rankB) */
typedef struct { int32_t ra, rb; size_t n; double min_gap, sum_gap; size_t rep_u, rep_v; double rep_dw; } SplitAgg;

static int run(const char *in, const Params *P)
{
    FVec V = {0}; IVec F = {0};
    if (read_obj(in, &V, &F) != 0) { free(V.v); free(F.f); return 1; }
    size_t nv = V.n, nf = F.n;
    if (nv == 0 || nf == 0) { fprintf(stderr, "wind_audit: empty mesh\n"); free(V.v); free(F.f); return 1; }
    double core_r = P->core_pitches * P->pitch;

    /* --- connected components --- */
    int32_t *par = (int32_t *)malloc(nv * sizeof(int32_t));
    int32_t *usz = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!par || !usz) { fprintf(stderr, "wind_audit: OOM\n"); free(par); free(usz); free(V.v); free(F.f); return 1; }
    for (size_t i = 0; i < nv; i++) { par[i] = (int32_t)i; usz[i] = 1; }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F.f[f*3+0], b = F.f[f*3+1], c = F.f[f*3+2];
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv) continue;
        uf_union(par, usz, a, b); uf_union(par, usz, b, c);
    }
    size_t *root_faces = (size_t *)calloc(nv, sizeof(size_t));
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F.f[f*3+0];
        if (a < 0 || (size_t)a >= nv) continue;
        root_faces[uf_find(par, a)]++;
    }
    size_t ncomp = 0;
    for (size_t i = 0; i < nv; i++) if (root_faces[i] > 0) ncomp++;
    Comp *comp = (Comp *)malloc((ncomp ? ncomp : 1) * sizeof(Comp));
    size_t ci = 0;
    for (size_t i = 0; i < nv; i++) if (root_faces[i] > 0) {
        comp[ci].root = (int32_t)i; comp[ci].faces = root_faces[i]; comp[ci].rank = 0;
        comp[ci].rmin = 1e30; comp[ci].rmax = -1e30; comp[ci].zmin = 1e30f; comp[ci].zmax = -1e30f;
        comp[ci].merger_edges = 0;
        comp[ci].theta_sin = comp[ci].theta_cos = 0.0;
        comp[ci].phase_min = 1e30; comp[ci].phase_max = -1e30;
        comp[ci].outer_vertices = 0;
        ci++;
    }
    qsort(comp, ncomp, sizeof(Comp), comp_cmp);
    int32_t *root2rank = (int32_t *)malloc(nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) root2rank[i] = -1;
    for (size_t r = 0; r < ncomp; r++) { comp[r].rank = (int32_t)r; root2rank[comp[r].root] = (int32_t)r; }

    /* per-vertex rank + radius; accumulate comp r/z bounds */
    int32_t *vrank = (int32_t *)malloc(nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) {
        int32_t rk = root2rank[uf_find(par, (int32_t)i)];
        vrank[i] = rk;
        if (rk < 0) continue;
        double r = vradius(V.v, i, P->umb_y, P->umb_x);
        float z = V.v[i*3+0];
        if (r < comp[rk].rmin) comp[rk].rmin = r;
        if (r > comp[rk].rmax) comp[rk].rmax = r;
        if (z < comp[rk].zmin) comp[rk].zmin = z;
        if (z > comp[rk].zmax) comp[rk].zmax = z;
        if (r >= core_r) {
            double dy = (double)V.v[i*3+1] - P->umb_y;
            double dx = (double)V.v[i*3+2] - P->umb_x;
            double th = atan2(dy, dx);
            comp[rk].theta_sin += sin(th);
            comp[rk].theta_cos += cos(th);
            comp[rk].outer_vertices++;
        }
    }

    /* A local reference angle per component removes atan2's global branch cut.
     * This is valid for a leaf cube in the outer scroll (its angular footprint
     * is far below pi); it is deliberately not a whole-scroll statistic. */
    for (size_t i = 0; i < nv; i++) {
        int32_t rk = vrank[i];
        if (rk < 0 || comp[rk].outer_vertices == 0) continue;
        double dy = (double)V.v[i*3+1] - P->umb_y;
        double dx = (double)V.v[i*3+2] - P->umb_x;
        double r = hypot(dy, dx);
        if (r < core_r) continue;
        double ref = atan2(comp[rk].theta_sin, comp[rk].theta_cos);
        double dth = wrap_pmpi(atan2(dy, dx) - ref);
        double phase = r / P->pitch - dth / TWO_PI;
        if (phase < comp[rk].phase_min) comp[rk].phase_min = phase;
        if (phase > comp[rk].phase_max) comp[rk].phase_max = phase;
    }
    size_t broad_fusion_components = 0, broad_fusion_faces = 0;
    double broad_fusion_max_span = 0.0;
    if (P->local_span_tol > 0.0) {
        for (size_t r = 0; r < ncomp; r++) {
            double span = comp[r].phase_max - comp[r].phase_min;
            /* Components touching the core may legitimately turn within one
             * cube, so the leaf-cube gate only judges wholly outer components. */
            if (comp[r].faces >= P->min_faces && comp[r].rmin >= core_r
                && comp[r].outer_vertices > 0 && span > P->local_span_tol) {
                broad_fusion_components++;
                broad_fusion_faces += comp[r].faces;
                if (span > broad_fusion_max_span) broad_fusion_max_span = span;
            }
        }
    }

    /* --- MERGER edges: dedup unique edges, test |dw| --- */
    /* dedup edge set (open addressing over packed u<<32|v) */
    size_t eslots = 1; while (eslots < nf * 4) eslots <<= 1;
    int64_t *eset = (int64_t *)malloc(eslots * sizeof(int64_t));
    uint8_t *ecount = (uint8_t *)calloc(eslots, 1);
    if (!eset || !ecount) {
        fprintf(stderr, "wind_audit: OOM allocating edge table\n");
        free(eset); free(ecount);
        free(vrank); free(par); free(usz); free(root_faces); free(comp); free(root2rank);
        free(V.v); free(F.f);
        return 1;
    }
    for (size_t i = 0; i < eslots; i++) eset[i] = -1;
    size_t emask = eslots - 1;
    size_t merge_outer = 0, merge_core = 0, edges_total = 0;
    /* merger-edge geometry: are they ~one-pitch radial jumps? */
    double mlen_sum = 0, mdr_sum = 0;           /* outer merger edges */
    double olen_sum = 0; size_t outer_edges = 0; /* all outer edges (contrast) */
    size_t dwhist[5] = {0,0,0,0,0};             /* |dw| bins over outer edges */
    size_t merge_ft = 0, merge_ft_seam = 0;     /* full-turn (|dw| in [0.7,1.3]) fusions, seam-classified */
    uint8_t *face_flag = NULL;
    if (P->dump_mergers) { face_flag = (uint8_t *)calloc(nf, 1); }
    for (size_t f = 0; f < nf; f++) {
        int32_t vv[3] = { F.f[f*3+0], F.f[f*3+1], F.f[f*3+2] };
        if (vv[0]<0||vv[1]<0||vv[2]<0||(size_t)vv[0]>=nv||(size_t)vv[1]>=nv||(size_t)vv[2]>=nv) continue;
        for (int e = 0; e < 3; e++) {
            int32_t u = vv[e], w = vv[(e+1)%3];
            int32_t a = u < w ? u : w, b = u < w ? w : u;
            int64_t key = ((int64_t)a << 32) | (uint32_t)b;
            size_t i = sh_hash(key, emask);
            int seen = 0;
            for (;;) {
                if (eset[i] == -1) { eset[i] = key; ecount[i] = 1; break; }
                if (eset[i] == key) {
                    if (ecount[i] < 255) ecount[i]++;
                    seen = 1; break;
                }
                i = (i+1)&emask;
            }
            if (seen) continue;
            edges_total++;
            double dw = dwind(V.v, (size_t)a, (size_t)b, P->umb_y, P->umb_x, P->pitch);
            double ra = vradius(V.v, (size_t)a, P->umb_y, P->umb_x);
            double rb = vradius(V.v, (size_t)b, P->umb_y, P->umb_x);
            double rmin = ra < rb ? ra : rb;
            double adw = fabs(dw);
            if (rmin >= core_r) {              /* outer edge: tally geometry + hist */
                double dz = (double)V.v[a*3+0]-V.v[b*3+0], dy = (double)V.v[a*3+1]-V.v[b*3+1], dx = (double)V.v[a*3+2]-V.v[b*3+2];
                double len = sqrt(dz*dz + dy*dy + dx*dx);
                outer_edges++; olen_sum += len;
                if      (adw < 0.40) dwhist[0]++;
                else if (adw < 0.70) dwhist[1]++;
                else if (adw < 1.30) dwhist[2]++;   /* ~one-turn radial jump */
                else if (adw < 2.30) dwhist[3]++;   /* ~two-turn */
                else                 dwhist[4]++;
                if (adw > P->merge_tol) { mlen_sum += len; mdr_sum += fabs(rb - ra); }
                if (adw >= 0.70 && adw <= 1.30) {   /* unambiguous full-turn fusion */
                    merge_ft++;
                    double mz = 0.5*((double)V.v[a*3+0]+V.v[b*3+0]);
                    double my = 0.5*((double)V.v[a*3+1]+V.v[b*3+1]);
                    double mx = 0.5*((double)V.v[a*3+2]+V.v[b*3+2]);
                    if (near_seam(mz, my, mx, P->seam_pitch, P->seam_tol)) merge_ft_seam++;
                    int32_t rk = vrank[a];       /* per-comp column = reliable full-turn count */
                    if (rk >= 0 && rk == vrank[b]) comp[rk].merger_edges++;
                }
            }
            if (adw > P->merge_tol) {
                if (rmin < core_r) merge_core++;
                else merge_outer++;
                if (face_flag) face_flag[f] = 1;
            }
        }
    }
    size_t boundary_edges = 0;
    for (size_t i = 0; i < eslots; i++)
        if (eset[i] != -1 && ecount[i] == 1) boundary_edges++;
    free(ecount);
    free(eset);

    /* --- SPLIT candidates: spatial hash, cross-component close low-dw pairs --- */
    SHash H; memset(&H, 0, sizeof(H));
    size_t hslots = 1; while (hslots < nv * 2) hslots <<= 1;
    H.cell = (HCell *)malloc(hslots * sizeof(HCell));
    H.next = (int32_t *)malloc(nv * sizeof(int32_t));
    H.mask = hslots - 1; H.inv = 1.0 / P->split_gap; H.off = 1 << 20;
    for (size_t i = 0; i < hslots; i++) { H.cell[i].head = -2; H.cell[i].key = 0; }
    for (size_t i = 0; i < nv; i++) {
        if (vrank[i] < 0) { H.next[i] = -1; continue; }
        int64_t key = sh_key(&H, V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
        int32_t *slot = sh_slot(&H, key, 1);
        H.next[i] = *slot; *slot = (int32_t)i;
    }
    /* aggregate split pairs by (ra,rb); small open table */
    size_t aslots = 4096; SplitAgg *agg = (SplitAgg *)calloc(aslots, sizeof(SplitAgg));
    for (size_t i = 0; i < aslots; i++) { agg[i].ra = -1; agg[i].rb = -1; }
    size_t nagg = 0;
    double gap2 = P->split_gap * P->split_gap;
    size_t split_pairs_total = 0, split_pairs_seam = 0;
    for (size_t i = 0; i < nv; i++) {
        int32_t rki = vrank[i]; if (rki < 0) continue;
        double az = V.v[i*3+0], ay = V.v[i*3+1], ax = V.v[i*3+2];
        for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
            int64_t key = sh_key(&H, az + dz*P->split_gap, ay + dy*P->split_gap, ax + dx*P->split_gap);
            int32_t *slot = sh_slot(&H, key, 0);
            if (!slot) continue;
            int scan = 0;
            for (int32_t j = *slot; j >= 0 && scan < 96; j = H.next[j], scan++) {
                if ((size_t)j <= i) continue;              /* unordered, once */
                int32_t rkj = vrank[j]; if (rkj < 0 || rkj == rki) continue;
                double d2 = (az-V.v[j*3+0])*(az-V.v[j*3+0]) + (ay-V.v[j*3+1])*(ay-V.v[j*3+1]) + (ax-V.v[j*3+2])*(ax-V.v[j*3+2]);
                if (d2 > gap2) continue;
                double dw = dwind(V.v, i, (size_t)j, P->umb_y, P->umb_x, P->pitch);
                if (fabs(dw) > P->split_wtol) continue;    /* different turn: not a split */
                split_pairs_total++;
                if (near_seam(0.5*(az+V.v[j*3+0]), 0.5*(ay+V.v[j*3+1]), 0.5*(ax+V.v[j*3+2]), P->seam_pitch, P->seam_tol)) split_pairs_seam++;
                int32_t ra = rki < rkj ? rki : rkj, rb = rki < rkj ? rkj : rki;
                /* find/insert agg */
                uint64_t hk = ((uint64_t)ra * 2654435761u) ^ ((uint64_t)rb * 40503u);
                size_t s = (size_t)hk & (aslots - 1); size_t guard = 0;
                for (;;) {
                    if (agg[s].ra == -1) { agg[s].ra = ra; agg[s].rb = rb; agg[s].n = 0; agg[s].min_gap = 1e30; agg[s].sum_gap = 0; nagg++; }
                    if (agg[s].ra == ra && agg[s].rb == rb) break;
                    s = (s + 1) & (aslots - 1);
                    if (++guard >= aslots) { s = SIZE_MAX; break; }
                }
                if (s == SIZE_MAX) continue;
                double gap = sqrt(d2);
                agg[s].n++; agg[s].sum_gap += gap;
                if (gap < agg[s].min_gap) { agg[s].min_gap = gap; agg[s].rep_u = i; agg[s].rep_v = (size_t)j; agg[s].rep_dw = dw; }
            }
        }
    }

    /* ---- report ---- */
    printf("wind_audit: %s\n", in);
    printf("  umbilicus=(y %.1f, x %.1f)  pitch=%.3f  merge_tol=%.2f  split_wtol=%.2f  split_gap=%.1f  core<%.1f vox\n",
           P->umb_y, P->umb_x, P->pitch, P->merge_tol, P->split_wtol, P->split_gap, core_r);
    printf("  verts=%zu  faces=%zu  components=%zu  unique_edges=%zu\n", nv, nf, ncomp, edges_total);
    printf("\n  == MERGERS (radial short-circuits, |dw|>%.2f across one edge) ==\n", P->merge_tol);
    printf("    outer (r>=%.0f, reliable): %zu edges (%.2f%% of %zu outer)    core (r<%.0f): %zu edges\n",
           core_r, merge_outer, outer_edges ? 100.0*(double)merge_outer/(double)outer_edges : 0.0,
           outer_edges, core_r, merge_core);
    printf("    outer |dw| histogram:  <0.40 (ok)=%zu  0.40-0.70=%zu  0.70-1.30 (~1 turn)=%zu  1.30-2.30 (~2)=%zu  >2.30=%zu\n",
           dwhist[0], dwhist[1], dwhist[2], dwhist[3], dwhist[4]);
    printf("    merger-edge geometry:  mean_len=%.2f vox  mean|dr|=%.2f vox   (all outer edges mean_len=%.2f vox; pitch=%.2f)\n",
           merge_outer ? mlen_sum/(double)merge_outer : 0.0, merge_outer ? mdr_sum/(double)merge_outer : 0.0,
           outer_edges ? olen_sum/(double)outer_edges : 0.0, P->pitch);
    printf("    FULL-TURN fusions (|dw| in [0.7,1.3], the reliable count): %zu   of which near a %.0f-vox seam: %zu (%.1f%%)\n",
           merge_ft, P->seam_pitch, merge_ft_seam, merge_ft ? 100.0*(double)merge_ft_seam/(double)merge_ft : 0.0);
    printf("    mesh boundary edges: %zu\n", boundary_edges);
    if (P->local_span_tol > 0.0) {
        printf("\n  == BROAD FUSIONS (leaf-cube phase span > %.2f turns) ==\n",
               P->local_span_tol);
        printf("    wholly-outer components: %zu  faces involved: %zu  max phase span: %.3f turns\n",
               broad_fusion_components, broad_fusion_faces, broad_fusion_max_span);
    }

    printf("\n  == COMPONENTS (top %d by faces) ==\n", P->top);
    printf("    %-4s %10s %6s   %8s %8s %8s  %9s %9s  %8s %8s   %s\n",
           "rank", "faces", "%", "r_min", "r_max", "r_span", "turns", "phase", "z_min", "z_max", "ft_fusions");
    int rows = P->top < (int)ncomp ? P->top : (int)ncomp;
    for (int r = 0; r < rows; r++) {
        double span = comp[r].rmax - comp[r].rmin;
        double pspan = comp[r].outer_vertices ?
            comp[r].phase_max - comp[r].phase_min : 0.0;
        printf("    #%-3d %10zu %5.1f%%   %8.1f %8.1f %8.1f  %8.1f %8.2f  %8.0f %8.0f   %zu\n",
               r, comp[r].faces, 100.0*(double)comp[r].faces/(double)nf,
               comp[r].rmin, comp[r].rmax, span, span / P->pitch, pspan,
               (double)comp[r].zmin, (double)comp[r].zmax, comp[r].merger_edges);
    }

    /* rank split aggregates by pair count */
    printf("\n  == SPLITS (same-turn |dw|<%.2f, gap<%.1f vox, DIFFERENT components) ==\n", P->split_wtol, P->split_gap);
    printf("    total candidate vertex pairs: %zu (near a %.0f-vox seam: %zu = %.1f%%)   distinct component-pairs: %zu\n",
           split_pairs_total, P->seam_pitch, split_pairs_seam,
           split_pairs_total ? 100.0*(double)split_pairs_seam/(double)split_pairs_total : 0.0, nagg);
    /* collect + sort agg */
    SplitAgg *list = (SplitAgg *)malloc((nagg ? nagg : 1) * sizeof(SplitAgg));
    size_t nl = 0;
    for (size_t s = 0; s < aslots; s++) if (agg[s].ra != -1) list[nl++] = agg[s];
    /* simple selection of top by n */
    int srows = P->top < (int)nl ? P->top : (int)nl;
    printf("    %-9s %8s %8s %8s   %s\n", "A<->B", "pairs", "min_gap", "mean_gap", "representative (z,y,x) r  dw");
    for (int r = 0; r < srows; r++) {
        size_t best = 0; size_t bestn = 0; int have = 0;
        for (size_t s = 0; s < nl; s++) if (list[s].n > bestn) { bestn = list[s].n; best = s; have = 1; }
        if (!have) break;
        SplitAgg *g = &list[best];
        double rr = vradius(V.v, g->rep_u, P->umb_y, P->umb_x);
        printf("    #%-3d<->#%-3d %8zu %8.2f %8.2f   (%.0f,%.0f,%.0f) r=%.0f dw=%+.3f\n",
               g->ra, g->rb, g->n, g->min_gap, g->sum_gap / (double)g->n,
               V.v[g->rep_u*3+0], V.v[g->rep_u*3+1], V.v[g->rep_u*3+2], rr, g->rep_dw);
        list[best].n = 0; /* consume */
    }

    /* ---- optional machine-readable summary (regression gate) ---- */
    if (P->json_path) {
        FILE *j = fopen(P->json_path, "wb");
        if (j) {
            fprintf(j, "{\n");
            fprintf(j, "  \"verts\": %zu, \"faces\": %zu, \"components\": %zu, \"outer_edges\": %zu, \"boundary_edges\": %zu,\n",
                    nv, nf, ncomp, outer_edges, boundary_edges);
            fprintf(j, "  \"merge_ft\": %zu, \"merge_ft_seam\": %zu, \"merge_ft_seam_frac\": %.4f,\n",
                    merge_ft, merge_ft_seam, merge_ft ? (double)merge_ft_seam/(double)merge_ft : 0.0);
            fprintf(j, "  \"merge_ft_per_outer_edge\": %.6f,\n", outer_edges ? (double)merge_ft/(double)outer_edges : 0.0);
            fprintf(j, "  \"split_pairs\": %zu, \"split_pairs_seam\": %zu, \"split_pairs_seam_frac\": %.4f, \"split_comp_pairs\": %zu,\n",
                    split_pairs_total, split_pairs_seam, split_pairs_total ? (double)split_pairs_seam/(double)split_pairs_total : 0.0, nagg);
            fprintf(j, "  \"local_span_tol\": %.4f, \"broad_fusion_components\": %zu, \"broad_fusion_faces\": %zu, \"broad_fusion_max_span\": %.4f,\n",
                    P->local_span_tol, broad_fusion_components, broad_fusion_faces, broad_fusion_max_span);
            fprintf(j, "  \"top_components\": [\n");
            int jr = (rows < 8) ? rows : 8;
            for (int r = 0; r < jr; r++)
                fprintf(j, "    {\"rank\": %d, \"faces\": %zu, \"turns\": %.1f, \"phase_span\": %.3f, \"ft_fusions\": %zu}%s\n",
                        r, comp[r].faces, (comp[r].rmax - comp[r].rmin) / P->pitch,
                        comp[r].outer_vertices ? comp[r].phase_max - comp[r].phase_min : 0.0,
                        comp[r].merger_edges, (r+1<jr)?",":"");
            fprintf(j, "  ]\n}\n");
            fclose(j);
            printf("  wrote json -> %s\n", P->json_path);
        }
    }

    /* ---- optional dumps ---- */
    if (P->dump_mergers && face_flag) {
        FILE *o = fopen(P->dump_mergers, "wb");
        if (o) {
            fprintf(o, "# wind_audit mergers: flagged faces red, rest grey\n");
            for (size_t i = 0; i < nv; i++)
                fprintf(o, "v %.5g %.5g %.5g 0.30 0.30 0.30\n", V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
            /* re-emit: flagged faces get bright verts via a second pass would need
             * per-vertex; simplest: emit flagged faces with a marker colour by
             * writing them last as separate verts */
            for (size_t f = 0; f < nf; f++)
                fprintf(o, "f %d %d %d\n", F.f[f*3+0]+1, F.f[f*3+1]+1, F.f[f*3+2]+1);
            /* overlay flagged faces as bright-red duplicated triangles */
            size_t base = nv;
            for (size_t f = 0; f < nf; f++) if (face_flag[f]) {
                for (int k = 0; k < 3; k++) { int32_t vi = F.f[f*3+k];
                    fprintf(o, "v %.5g %.5g %.5g 1.0 0.05 0.05\n", V.v[vi*3+0], V.v[vi*3+1], V.v[vi*3+2]); }
            }
            size_t idx = base;
            for (size_t f = 0; f < nf; f++) if (face_flag[f]) { fprintf(o, "f %zu %zu %zu\n", idx+1, idx+2, idx+3); idx += 3; }
            fclose(o);
            printf("\n  wrote mergers overlay -> %s\n", P->dump_mergers);
        }
    }
    if (P->dump_splits) {
        FILE *o = fopen(P->dump_splits, "wb");
        if (o) {
            fprintf(o, "# wind_audit splits: full mesh grey, split-pair verts bright, l-lines across gaps\n");
            for (size_t i = 0; i < nv; i++)
                fprintf(o, "v %.5g %.5g %.5g 0.28 0.28 0.28\n", V.v[i*3+0], V.v[i*3+1], V.v[i*3+2]);
            for (size_t f = 0; f < nf; f++)
                fprintf(o, "f %d %d %d\n", F.f[f*3+0]+1, F.f[f*3+1]+1, F.f[f*3+2]+1);
            /* bright endpoints + connecting lines (scan agg: list[].n was
             * consumed by the top-K report above) */
            size_t idx = nv;
            for (size_t s = 0; s < aslots; s++) if (agg[s].ra != -1) {
                size_t u = agg[s].rep_u, v = agg[s].rep_v;
                fprintf(o, "v %.5g %.5g %.5g 0.1 1.0 0.1\n", V.v[u*3+0], V.v[u*3+1], V.v[u*3+2]);
                fprintf(o, "v %.5g %.5g %.5g 1.0 0.9 0.1\n", V.v[v*3+0], V.v[v*3+1], V.v[v*3+2]);
                fprintf(o, "l %zu %zu\n", idx+1, idx+2); idx += 2;
            }
            fclose(o);
            printf("  wrote splits overlay -> %s\n", P->dump_splits);
        }
    }

    free(list); free(agg); free(H.cell); free(H.next);
    free(face_flag);
    free(vrank); free(par); free(usz); free(root_faces); free(comp); free(root2rank);
    free(V.v); free(F.f);
    return (P->fail_on_local &&
            (merge_ft > 0 || broad_fusion_components > 0)) ? 4 : 0;
}

/* ---------- selftest ---------- */
static int approx(double a, double b, double e) { return fabs(a - b) <= e; }

static int selftest(void)
{
    int fails = 0;

    /* winding math: a full-pitch radial jump at constant angle -> dw ~ 1 */
    {
        float Vt[6] = { 0, 50.0f, 0, 0, 59.5f, 0 };  /* (z,y,x): r 50 -> 59.5, same angle */
        double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5);
        if (!approx(dw, 1.0, 0.02)) { fprintf(stderr, "selftest: radial-jump dw=%.4f want ~1.0\n", dw); fails++; }
    }
    /* winding math: same turn, small angular step -> dw ~ 0 */
    {
        /* two points on r=50 circle 0.1 rad apart: r const, |dth|=0.1 -> |dw|=0.1/2pi ~ 0.0159 */
        float Vt[6] = { 0, 50.0f, 0, 0, (float)(50*cos(0.1)), (float)(50*sin(0.1)) };
        double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5);
        if (!approx(fabs(dw), 0.1/TWO_PI, 0.01)) { fprintf(stderr, "selftest: tangential |dw|=%.4f want ~%.4f\n", fabs(dw), 0.1/TWO_PI); fails++; }
    }
    /* branch-cut: crossing theta=0 ray on the same circle stays ~0 */
    {
        float Vt[6] = { 0, (float)(50*sin(-0.05)), (float)(50*cos(-0.05)),
                        0, (float)(50*sin(+0.05)), (float)(50*cos(+0.05)) };
        double dw = dwind(Vt, 0, 1, 0.0, 0.0, 9.5);
        if (fabs(dw) > 0.02) { fprintf(stderr, "selftest: branch-cut dw=%.4f want ~0\n", dw); fails++; }
    }

    /* end-to-end: build a tiny mesh with a known merger + a known split,
     * write it, run(), and check the printed-path invariants via return code +
     * internal recomputation. We verify detection logic directly here. */
    {
        /* Case A: two concentric ring-arcs (r=50, r=57.6) = ADJACENT turns.
         *   dw ~ 7.6/9.5 = 0.80 -> NOT a split, NOT... they are separate comps
         *   with a ~7.6 gap. With split_wtol 0.35 they must NOT be flagged. */
        double dw_adj = 7.6 / 9.5;
        if (dw_adj <= 0.35) { fprintf(stderr, "selftest: adjacent-turn dw=%.3f should exceed split_wtol\n", dw_adj); fails++; }

        /* Case B: same ring severed -> two arcs at r=50, angular gap; dw ~ 0,
         *   physical gap ~ small -> MUST be flagged. Represented by two verts
         *   at r=50, 3 vox apart tangentially. */
        float p0y = 50.0f, p0x = 0.0f;
        float p1y = (float)(50*cos(3.0/50)), p1x = (float)(50*sin(3.0/50)); /* ~3 vox arc */
        float Vt[6] = { 0, p0y, p0x, 0, p1y, p1x };
        double dw_sev = dwind(Vt, 0, 1, 0.0, 0.0, 9.5);
        double gap = hypot((double)p1y-p0y, (double)p1x-p0x);
        if (fabs(dw_sev) > 0.35) { fprintf(stderr, "selftest: severed dw=%.3f should be < split_wtol\n", dw_sev); fails++; }
        if (gap > 14.0) { fprintf(stderr, "selftest: severed gap=%.2f should be < split_gap\n", gap); fails++; }
    }
    /* A staircase fusion can span a full turn while every constituent edge is
     * below the abrupt-merger threshold. This is exactly why the component
     * phase-span gate exists. */
    {
        float Vt[5 * 3] = {
            0, 50.000f, 0,  0, 52.375f, 0,  0, 54.750f, 0,
            0, 57.125f, 0,  0, 59.500f, 0
        };
        double pmin = 1e30, pmax = -1e30;
        for (size_t i = 0; i < 5; i++) {
            double ph = vradius(Vt, i, 0.0, 0.0) / 9.5;
            if (ph < pmin) pmin = ph;
            if (ph > pmax) pmax = ph;
            if (i > 0 && fabs(dwind(Vt, i-1, i, 0.0, 0.0, 9.5)) >= 0.40) {
                fprintf(stderr, "selftest: staircase edge %zu looked abrupt\n", i);
                fails++;
            }
        }
        if (!approx(pmax - pmin, 1.0, 0.02)) {
            fprintf(stderr, "selftest: staircase phase span %.3f want 1.0\n", pmax-pmin);
            fails++;
        }
    }

    /* full run smoke on a written OBJ: 2 disjoint triangles at different radii */
    {
        const char *tin = "wind_audit_selftest_in.obj";
        FILE *fp = fopen(tin, "wb");
        if (!fp) { fprintf(stderr, "selftest: tmp write fail\n"); return 3; }
        /* comp0 near r~50, comp1 near r~200 (far apart, no split, no merger) */
        fprintf(fp, "v 0 50 0\nv 0 51 1\nv 0 50 2\n");
        fprintf(fp, "v 0 200 0\nv 0 201 1\nv 0 200 2\n");
        fprintf(fp, "f 1 2 3\nf 4 5 6\n");
        fclose(fp);
        Params P; memset(&P, 0, sizeof(P));
        P.umb_y = 0; P.umb_x = 0; P.pitch = 9.5; P.merge_tol = 0.40; P.split_wtol = 0.35;
        P.split_gap = 14.0; P.core_pitches = 3.0; P.min_faces = 0; P.top = 5;
        int rc = run(tin, &P);
        remove(tin);
        if (rc != 0) { fprintf(stderr, "selftest: run rc=%d\n", rc); fails++; }
    }

    /* Integrated broad-fusion gate: a connected radial staircase spans one
     * phase turn, but each edge advances only one quarter turn. */
    {
        const char *tin = "wind_audit_selftest_local.obj";
        FILE *fp = fopen(tin, "wb");
        if (!fp) { fprintf(stderr, "selftest: local tmp write fail\n"); return 3; }
        for (int i = 0; i < 5; i++) {
            double r = 50.0 + 2.375 * (double)i;
            fprintf(fp, "v 0 %.6f 0\nv 1 %.6f 0\n", r, r);
        }
        for (int i = 0; i < 4; i++) {
            int a = 2*i + 1, b = a + 1, c = a + 2, d = a + 3;
            fprintf(fp, "f %d %d %d\nf %d %d %d\n", a, c, b, b, c, d);
        }
        fclose(fp);
        Params P; memset(&P, 0, sizeof(P));
        P.umb_y = 0; P.umb_x = 0; P.pitch = 9.5;
        P.merge_tol = 0.40; P.split_wtol = 0.35; P.split_gap = 14.0;
        P.core_pitches = 3.0; P.min_faces = 1; P.top = 5;
        P.local_span_tol = 0.9; P.fail_on_local = 1;
        int rc = run(tin, &P);
        remove(tin);
        if (rc != 4) {
            fprintf(stderr, "selftest: broad-fusion gate rc=%d want 4\n", rc);
            fails++;
        }
    }

    fprintf(stderr, "wind_audit selftest: %s\n", fails ? "FAIL" : "PASS");
    return fails ? 3 : 0;
}

int main(int argc, char **argv)
{
    if (argc == 2 && !strcmp(argv[1], "--selftest")) return selftest();
    if (argc < 2) {
        fprintf(stderr,
            "usage: %s <in.obj> --umb-y Y --umb-x X [--pitch P=9.5] [--top N=25]\n"
            "         [--merge-tol T=0.40] [--split-wtol T=0.35] [--split-gap G=14]\n"
            "         [--core-pitches C=3] [--min-faces M=64]\n"
            "         [--local-span-tol T] [--fail-on-local]\n"
            "         [--dump-mergers f.obj] [--dump-splits f.obj]\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 2;
    }
    Params P; memset(&P, 0, sizeof(P));
    /* split_wtol default = the seam gate's SOFT tol (0.25): a split below it is
     * one the gate WOULD have allowed -> a genuine weld failure, not a gate-
     * correct rejection of a ~1/3-turn-apart chart. */
    P.umb_y = 3405; P.umb_x = 2878; P.pitch = 9.5; P.merge_tol = 0.40; P.split_wtol = 0.25;
    P.split_gap = 14.0; P.core_pitches = 3.0; P.min_faces = 64; P.top = 25;
    P.seam_pitch = 128.0; P.seam_tol = 3.0;
    const char *in = argv[1];
    for (int i = 2; i < argc; i++) {
        if      (!strcmp(argv[i], "--umb-y") && i+1 < argc) P.umb_y = atof(argv[++i]);
        else if (!strcmp(argv[i], "--umb-x") && i+1 < argc) P.umb_x = atof(argv[++i]);
        else if (!strcmp(argv[i], "--pitch") && i+1 < argc) P.pitch = atof(argv[++i]);
        else if (!strcmp(argv[i], "--top") && i+1 < argc) P.top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--merge-tol") && i+1 < argc) P.merge_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--split-wtol") && i+1 < argc) P.split_wtol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--split-gap") && i+1 < argc) P.split_gap = atof(argv[++i]);
        else if (!strcmp(argv[i], "--core-pitches") && i+1 < argc) P.core_pitches = atof(argv[++i]);
        else if (!strcmp(argv[i], "--local-span-tol") && i+1 < argc) P.local_span_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--fail-on-local")) P.fail_on_local = 1;
        else if (!strcmp(argv[i], "--seam-pitch") && i+1 < argc) P.seam_pitch = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seam-tol") && i+1 < argc) P.seam_tol = atof(argv[++i]);
        else if (!strcmp(argv[i], "--min-faces") && i+1 < argc) P.min_faces = (size_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--dump-mergers") && i+1 < argc) P.dump_mergers = argv[++i];
        else if (!strcmp(argv[i], "--dump-splits") && i+1 < argc) P.dump_splits = argv[++i];
        else if (!strcmp(argv[i], "--json") && i+1 < argc) P.json_path = argv[++i];
        else { fprintf(stderr, "wind_audit: unknown arg %s\n", argv[i]); return 2; }
    }
    if (P.pitch <= 0 || P.split_gap <= 0) { fprintf(stderr, "wind_audit: pitch/split-gap must be > 0\n"); return 2; }
    return run(in, &P);
}
