#include "scaffold.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ calib */

void Scaffold_calib_default(ScaffoldCalib *c)
{
    assert(c);
    c->spiral_a = 0.85;
    c->spiral_b = -9.5;
    c->pitch    = 9.5;
    c->sense    = -1;
    c->axis_point[0] = 0.0f; c->axis_point[1] = 3405.0f; c->axis_point[2] = 2878.0f;
    c->axis_dir[0]   = 1.0f; c->axis_dir[1]   = 0.0f;    c->axis_dir[2]   = 0.0f;
}

/* crude key scan (same convention as seam_own's so_jfind). */
static const char *sc_jfind(const char *s, const char *key)
{
    char pat[128];
    snprintf(pat, sizeof(pat), "\"%s\"", key);
    const char *p = strstr(s, pat);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(pat), ':');
    return p != NULL ? p + 1 : NULL;
}

int Scaffold_read_calib(const char *placed_dir, ScaffoldCalib *c)
{
    assert(placed_dir && c);
    Scaffold_calib_default(c);

    char path[1024];
    snprintf(path, sizeof(path), "%s/placed_index.json", placed_dir);
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';

    const char *p = NULL;
    double a = 0, b = 0, cc = 0;
    if ((p = sc_jfind(buf, "axis_point_zyx")) != NULL
        && sscanf(p, " [ %lf , %lf , %lf", &a, &b, &cc) == 3) {
        c->axis_point[0] = (float)a; c->axis_point[1] = (float)b; c->axis_point[2] = (float)cc;
    }
    if ((p = sc_jfind(buf, "axis_dir_zyx")) != NULL
        && sscanf(p, " [ %lf , %lf , %lf", &a, &b, &cc) == 3) {
        c->axis_dir[0] = (float)a; c->axis_dir[1] = (float)b; c->axis_dir[2] = (float)cc;
    }
    if ((p = sc_jfind(buf, "pitch")) != NULL && sscanf(p, " %lf", &a) == 1) c->pitch = a;
    if ((p = sc_jfind(buf, "spiral_a")) != NULL && sscanf(p, " %lf", &a) == 1) c->spiral_a = a;
    if ((p = sc_jfind(buf, "spiral_b")) != NULL && sscanf(p, " %lf", &a) == 1) c->spiral_b = a;
    /* sense lives in the calibration sub-block: "sense": -1 (an int). */
    { int s = 0; if ((p = sc_jfind(buf, "sense")) != NULL && sscanf(p, " %d", &s) == 1
                  && (s == 1 || s == -1)) c->sense = s; }
    return 0;
}

/* ------------------------------------------------------------------ kstats */

static int cmp_double(const void *pa, const void *pb)
{
    double a = *(const double *)pa, b = *(const double *)pb;
    return a < b ? -1 : (a > b ? 1 : 0);
}

/* Open-addressing (cube,gid) -> [w_min, w_max] accumulator. */
typedef struct { uint64_t *key; double *wmin; double *wmax; size_t cap, mask, n; } GroupMap;

static void gm_init(Arena_T ar, GroupMap *g, size_t hint)
{
    size_t cap = 64; while (cap < hint * 2) cap <<= 1;
    g->cap = cap; g->mask = cap - 1; g->n = 0;
    g->key  = (uint64_t *)ARENA_ALLOC(ar, (size_t)(cap * sizeof(uint64_t)));
    g->wmin = (double *)ARENA_ALLOC(ar, (size_t)(cap * sizeof(double)));
    g->wmax = (double *)ARENA_ALLOC(ar, (size_t)(cap * sizeof(double)));
    for (size_t i = 0; i < cap; i++) g->key[i] = UINT64_MAX;
}

static void gm_add(GroupMap *g, int32_t cube, int32_t gid, double w)
{
    uint64_t k = ((uint64_t)(uint32_t)cube << 20) | (uint64_t)((uint32_t)gid & 0xFFFFFu);
    uint64_t h = (k * 1099511628211ull) & g->mask;
    for (;;) {
        if (g->key[h] == k) {
            if (w < g->wmin[h]) g->wmin[h] = w;
            if (w > g->wmax[h]) g->wmax[h] = w;
            return;
        }
        if (g->key[h] == UINT64_MAX) {
            g->key[h] = k; g->wmin[h] = w; g->wmax[h] = w; g->n++;
            return;
        }
        h = (h + 1) & g->mask;
    }
}

int Scaffold_kstats(Arena_T arena, const PieceSet *ps, const ScaffoldCalib *c,
                    ScaffoldKStats *out, const char *hist_path)
{
    assert(arena && ps && c && out);
    memset(out, 0, sizeof(*out));
    if (ps->nf == 0 || ps->nv == 0 || ps->phi == NULL) return -1;

    Arena_Mark mark = Arena_save(arena);

    const double PI = SCAFFOLD_2PI * 0.5;

    /* vert -> cube (verts load contiguous per cube). */
    int32_t *vcube = (int32_t *)ARENA_ALLOC(arena, (size_t)(ps->nv * sizeof(int32_t)));
    for (size_t ci = 0; ci < ps->n_cubes; ci++)
        for (size_t v = ps->cube_voff[ci]; v < ps->cube_voff[ci + 1]; v++)
            vcube[v] = (int32_t)ci;

    /* unwound turn coordinate per vert (only meaningful where phi is set). */
    double *w = (double *)ARENA_ALLOC(arena, (size_t)(ps->nv * sizeof(double)));
    for (size_t v = 0; v < ps->nv; v++) w[v] = Scaffold_unwound(c, (double)ps->phi[v]);

    uint8_t *vused = (uint8_t *)ARENA_CALLOC(arena, (size_t)ps->nv, 1L);
    double *spread = (double *)ARENA_ALLOC(arena, (size_t)(ps->nf * sizeof(double)));
    size_t nsp = 0, face_ge_pi = 0;

    for (size_t fidx = 0; fidx < ps->nf; fidx++) {
        int32_t a = ps->faces[fidx * 3 + 0];
        int32_t b = ps->faces[fidx * 3 + 1];
        int32_t cvt = ps->faces[fidx * 3 + 2];
        double wa = w[a], wb = w[b], wc = w[cvt];
        double lo = wa < wb ? (wa < wc ? wa : wc) : (wb < wc ? wb : wc);
        double hi = wa > wb ? (wa > wc ? wa : wc) : (wb > wc ? wb : wc);
        double sp = hi - lo;
        spread[nsp++] = sp;
        if (sp >= PI) face_ge_pi++;
        vused[a] = vused[b] = vused[cvt] = 1;
    }

    /* k coverage over used verts + group accumulation. */
    int kmin = 0x7fffffff, kmax = -0x7fffffff;
    size_t nv_used = 0;
    GroupMap gm; gm_init(arena, &gm, ps->n_cubes * 32 + 64);
    for (size_t v = 0; v < ps->nv; v++) {
        if (!vused[v]) continue;
        nv_used++;
        int k = (int)floor(w[v] / SCAFFOLD_2PI);
        if (k < kmin) kmin = k;
        if (k > kmax) kmax = k;
        int32_t gid = ps->gid != NULL ? ps->gid[v] : -1;
        if (gid >= 0) gm_add(&gm, vcube[v], gid, w[v]);
    }

    size_t grp_ge_pi = 0;
    for (size_t i = 0; i < gm.cap; i++)
        if (gm.key[i] != UINT64_MAX && (gm.wmax[i] - gm.wmin[i]) >= PI) grp_ge_pi++;

    qsort(spread, nsp, sizeof(double), cmp_double);

    out->nv_used = nv_used;
    out->nf_used = nsp;
    out->k_min = kmin; out->k_max = kmax;
    out->k_extent = (nv_used > 0) ? (kmax - kmin + 1) : 0;
    out->face_spread_ge_pi = face_ge_pi;
    out->face_spread_frac = nsp > 0 ? (double)face_ge_pi / (double)nsp : 0.0;
    out->p50_face_spread = nsp > 0 ? spread[nsp / 2] : 0.0;
    out->p99_face_spread = nsp > 0 ? spread[(size_t)((double)nsp * 0.99)] : 0.0;
    out->n_groups = gm.n;
    out->group_spread_ge_pi = grp_ge_pi;

    /* optional k x cube-z histogram (COVERED face counts). */
    if (hist_path != NULL && nv_used > 0) {
        int nk = kmax - kmin + 1;
        /* distinct cube-z origins */
        long zvals[512]; int nz = 0;
        for (size_t ci = 0; ci < ps->n_cubes && nz < 512; ci++) {
            long oz = ps->cube_org[ci][0];
            int seen = 0; for (int t = 0; t < nz; t++) if (zvals[t] == oz) { seen = 1; break; }
            if (!seen) zvals[nz++] = oz;
        }
        /* sort zvals ascending (small nz) */
        for (int i = 1; i < nz; i++) { long t = zvals[i]; int j = i - 1;
            while (j >= 0 && zvals[j] > t) { zvals[j+1] = zvals[j]; j--; } zvals[j+1] = t; }
        size_t *H = (size_t *)ARENA_CALLOC(arena, (size_t)((size_t)nk * (size_t)(nz > 0 ? nz : 1)),
                                           (size_t)sizeof(size_t));
        for (size_t fidx = 0; fidx < ps->nf; fidx++) {
            int32_t cf = ps->face_cube[fidx];
            long oz = ps->cube_org[cf][0];
            int zi = 0; for (; zi < nz; zi++) if (zvals[zi] == oz) break;
            double wc = (w[ps->faces[fidx*3+0]] + w[ps->faces[fidx*3+1]]
                         + w[ps->faces[fidx*3+2]]) / 3.0;
            int k = (int)floor(wc / SCAFFOLD_2PI) - kmin;
            if (k >= 0 && k < nk && zi < nz) H[(size_t)k * (size_t)nz + (size_t)zi]++;
        }
        FILE *hf = fopen(hist_path, "w");
        if (hf != NULL) {
            fprintf(hf, "{\n  \"k_min\": %d, \"k_max\": %d,\n  \"cube_z\": [", kmin, kmax);
            for (int zi = 0; zi < nz; zi++) fprintf(hf, "%s%ld", zi ? ", " : "", zvals[zi]);
            fprintf(hf, "],\n  \"counts_k_by_z\": [\n");
            for (int k = 0; k < nk; k++) {
                fprintf(hf, "    [");
                for (int zi = 0; zi < nz; zi++)
                    fprintf(hf, "%s%zu", zi ? "," : "", H[(size_t)k * (size_t)nz + (size_t)zi]);
                fprintf(hf, "]%s\n", k + 1 < nk ? "," : "");
            }
            fprintf(hf, "  ]\n}\n");
            fclose(hf);
        }
    }

    Arena_restore(arena, mark);
    return 0;
}

/* ------------------------------------------------------------------ selftest */

/* Build a tiny synthetic PieceSet: one cube, a strip of quads laid on an
 * analytic spiral so k is exactly known. Verts (z,y,x); phi = -theta - 2pi*k
 * with sense -1 so unwound w = theta + 2pi*k increases outward. */
int Scaffold_selftest(void)
{
    int fails = 0;
    Arena_T ar = Arena_new();

    ScaffoldCalib c; Scaffold_calib_default(&c);   /* sense -1, pitch 9.5 */

    /* Wrap-index math: sense -1, phi = -(theta + 2pi*k). */
    {
        double phi = -(0.3 + SCAFFOLD_2PI * 7.0);
        int k = Scaffold_wrap_index(&c, phi);
        if (k != 7) { fprintf(stderr, "[scaffold] wrap_index got %d want 7\n", k); fails++; }
        double phi2 = -(SCAFFOLD_2PI * 0.0 + 0.0);
        if (Scaffold_wrap_index(&c, phi2) != 0) { fprintf(stderr, "[scaffold] wrap_index 0 fail\n"); fails++; }
    }

    /* Synthetic piece set: a SPIRAL RIBBON of 4 turns. Faces run ALONG the
     * spiral (consecutive arc samples, small d-phi) and across a narrow z
     * width -- exactly like a real scroll sheet, where radially-adjacent turns
     * are NOT mesh-connected. k = floor(t) spans 0..3; gid = floor(t) (one
     * group per turn, so groups legitimately span ~2pi). */
    {
        int NT = 200, WROWS = 2;               /* 200 arc samples, 2 z rows */
        double turns = 3.95;                    /* just under 4 => k in {0,1,2,3} */
        size_t nv = (size_t)NT * (size_t)WROWS;
        size_t nf = (size_t)(NT - 1) * (size_t)(WROWS - 1) * 2;
        PieceSet ps; memset(&ps, 0, sizeof ps);
        ps.verts = (float *)ARENA_ALLOC(ar, (size_t)(nv * 3 * sizeof(float)));
        ps.uv    = (float *)ARENA_ALLOC(ar, (size_t)(nv * 2 * sizeof(float)));
        ps.phi   = (float *)ARENA_ALLOC(ar, (size_t)(nv * sizeof(float)));
        ps.normals = (float *)ARENA_CALLOC(ar, (size_t)(nv * 3), (size_t)sizeof(float));
        ps.gid   = (int32_t *)ARENA_ALLOC(ar, (size_t)(nv * sizeof(int32_t)));
        ps.faces = (int32_t *)ARENA_ALLOC(ar, (size_t)(nf * 3 * sizeof(int32_t)));
        ps.face_cube = (int32_t *)ARENA_CALLOC(ar, (size_t)nf, (size_t)sizeof(int32_t));
        ps.ids = (char (*)[48])ARENA_ALLOC(ar, (size_t)(1 * 48));
        memcpy(ps.ids[0], "z00000_y00000_x00000", 21);
        ps.cube_voff = (size_t *)ARENA_ALLOC(ar, (size_t)(2 * sizeof(size_t)));
        ps.cube_voff[0] = 0; ps.cube_voff[1] = nv;
        ps.cube_org = (long (*)[3])ARENA_ALLOC(ar, (size_t)(3 * sizeof(long)));
        ps.cube_org[0][0] = 0; ps.cube_org[0][1] = 0; ps.cube_org[0][2] = 0;
        ps.n_cubes = 1; ps.nv = nv; ps.nf = nf;

        for (int i = 0; i < NT; i++) {
            double t = (double)i / (double)(NT - 1) * turns;   /* turns along spiral */
            double theta = SCAFFOLD_2PI * t;
            double r = c.spiral_a + (-c.spiral_b) * t;
            for (int row = 0; row < WROWS; row++) {
                size_t v = (size_t)i * (size_t)WROWS + (size_t)row;
                ps.verts[v*3+0] = (float)row;                       /* z (width) */
                ps.verts[v*3+1] = (float)(3405.0 + r * sin(theta)); /* y */
                ps.verts[v*3+2] = (float)(2878.0 + r * cos(theta)); /* x */
                ps.phi[v] = (float)(-theta);        /* sense -1 => w = theta = 2pi*t */
                ps.uv[v*2+0] = (float)(r * theta); ps.uv[v*2+1] = (float)row;
                ps.gid[v] = (int32_t)floor(t);      /* one group per turn */
            }
        }
        size_t fi = 0;
        for (int i = 0; i < NT - 1; i++)
            for (int row = 0; row < WROWS - 1; row++) {
                int32_t a = (int32_t)(i * WROWS + row);
                int32_t b = (int32_t)((i + 1) * WROWS + row);
                int32_t cc = (int32_t)(i * WROWS + row + 1);
                int32_t d = (int32_t)((i + 1) * WROWS + row + 1);
                ps.faces[fi*3+0]=a; ps.faces[fi*3+1]=b; ps.faces[fi*3+2]=cc; fi++;
                ps.faces[fi*3+0]=b; ps.faces[fi*3+1]=d; ps.faces[fi*3+2]=cc; fi++;
            }
        ps.u_min = 0; ps.u_max = 1; ps.v_min = 0; ps.v_max = WROWS;

        ScaffoldKStats st;
        if (Scaffold_kstats(ar, &ps, &c, &st, NULL) != 0) {
            fprintf(stderr, "[scaffold] kstats failed\n"); fails++;
        } else {
            if (st.k_extent != 4) { fprintf(stderr, "[scaffold] k_extent %d want 4\n", st.k_extent); fails++; }
            if (st.k_min != 0 || st.k_max != 3) { fprintf(stderr, "[scaffold] k range [%d,%d] want [0,3]\n", st.k_min, st.k_max); fails++; }
            /* faces span ~1/200 of 4 turns in phi => spread << pi => 0 over-pi. */
            if (st.face_spread_ge_pi != 0) { fprintf(stderr, "[scaffold] %zu faces >= pi (want 0)\n", st.face_spread_ge_pi); fails++; }
            if (st.n_groups != 4) { fprintf(stderr, "[scaffold] n_groups %zu want 4\n", st.n_groups); fails++; }
        }
    }

    Arena_dispose(&ar);
    if (fails == 0) fprintf(stderr, "[scaffold] selftest PASSED\n");
    else fprintf(stderr, "[scaffold] selftest FAILED (%d)\n", fails);
    return fails;
}
