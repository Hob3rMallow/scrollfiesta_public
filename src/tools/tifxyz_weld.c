/*
 * tifxyz_weld.c -- weld z-slab tifxyz strips into full-height per-wrap segments.
 *
 * The full-scroll pipeline processes the band as independent z-slabs (1x21x21
 * strips), each exported to its own tifxyz atlas from ONE shared global solve
 * (scroll_unroll --z-range). Because the strips share the u lattice and cover
 * disjoint v (world z) bands, welding is a vertical stack: for each wrap k this
 * tool places every strip's k-segment at its absolute (u,v) into one
 * full-height canvas (first cover wins any halo overlap) and writes a unified
 * atlas. The C form of python/scripts/weld_slab_atlas.py -- reassembly is the
 * scaling mechanism for the 241k-cube whole scroll.
 *
 * usage: tifxyz_weld <out_root> <strip_root...>   (each root holds seg/<name>/)
 *        tifxyz_weld --selftest
 *
 * out_root gets seg/welded_w###/{x,y,z,mask,provenance}.tif + meta.json and
 * atlas.json. Segments carry absolute origin_uv; a wrap's members are all the
 * "*_w<k>_*" segments across the input roots.
 */
#include "../common/ves_platform.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
#endif

#include "../common/arena.h"
#include "../common/tiff_io.h"

/* crude json key scan (same convention as the placed_index / tifxyz_render
 * readers): find "key" then the ':' after it. */
static const char *jfind(const char *s, const char *key)
{
    char pat[128];
    snprintf(pat, sizeof(pat), "\"%s\"", key);
    const char *p = strstr(s, pat);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(pat), ':');
    return p != NULL ? p + 1 : NULL;
}

typedef struct {
    char   dir[1024];        /* full segment dir (…/seg/<name>) */
    int    k;                /* wrap index */
    double ou, ov, du, dv;   /* absolute origin (vox) + grid step */
    int    W, H;
} SegRef;

/* wrap k from a segment name "…_w<digits>_…" (first such run). -1 if none. */
static int parse_wrap_k(const char *name)
{
    const char *p = name;
    while ((p = strstr(p, "_w")) != NULL) {
        if (p[2] >= '0' && p[2] <= '9') return atoi(p + 2);
        p += 2;
    }
    return -1;
}

/* read origin_uv/du/dv/flip from meta.json + W/H from mask.tif. 0 on success. */
static int read_ref(Arena_T scr, const char *seg_dir, SegRef *r)
{
    char path[1200];
    snprintf(path, sizeof(path), "%s/meta.json", seg_dir);
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';
    const char *p;
    r->du = r->dv = 1.0; r->ou = r->ov = 0.0;
    int flip_u = 0, flip_v = 0;
    if ((p = jfind(buf, "du")) != NULL) sscanf(p, " %lf", &r->du);
    if ((p = jfind(buf, "dv")) != NULL) sscanf(p, " %lf", &r->dv);
    if ((p = jfind(buf, "origin_uv")) != NULL) sscanf(p, " [ %lf , %lf", &r->ou, &r->ov);
    if ((p = jfind(buf, "flip_u")) != NULL) sscanf(p, " %d", &flip_u);
    if ((p = jfind(buf, "flip_v")) != NULL) sscanf(p, " %d", &flip_v);
    if (flip_u || flip_v) {
        fprintf(stderr, "tifxyz_weld: flip_u/flip_v unsupported (%s)\n", seg_dir);
        return -1;
    }
    snprintf(path, sizeof(path), "%s/mask.tif", seg_dir);
    uint8_t *m = NULL; int D = 0, H = 0, W = 0;
    if (TiffIO_load(scr, path, &m, &D, &H, &W) != 0) return -1;
    r->W = W; r->H = H;
    snprintf(r->dir, sizeof(r->dir), "%s", seg_dir);
    return 0;
}

/* enumerate <root>/seg/<name>/ dirs, append matching SegRefs. */
static void scan_segs(Arena_T scr, const char *root,
                      SegRef **refs, size_t *n, size_t *cap)
{
    char segroot[1100];
    snprintf(segroot, sizeof(segroot), "%s/seg", root);
#ifdef _WIN32
    char glob[1200]; snprintf(glob, sizeof(glob), "%s\\*", segroot);
    WIN32_FIND_DATAA fd; HANDLE h = FindFirstFileA(glob, &fd);
    if (h == INVALID_HANDLE_VALUE) return;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) continue;
        if (fd.cFileName[0] == '.') continue;
        int k = parse_wrap_k(fd.cFileName);
        if (k < 0) continue;
        char sd[1200]; snprintf(sd, sizeof(sd), "%s/%s", segroot, fd.cFileName);
        Arena_free(scr);
        SegRef r; r.k = k;
        if (read_ref(scr, sd, &r) != 0) continue;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 256;
            *refs = (SegRef *)realloc(*refs, *cap * sizeof(SegRef)); }
        (*refs)[(*n)++] = r;
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(segroot);
    if (d == NULL) return;
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (e->d_name[0] == '.') continue;
        int k = parse_wrap_k(e->d_name);
        if (k < 0) continue;
        char sd[1200]; snprintf(sd, sizeof(sd), "%s/%s", segroot, e->d_name);
        Arena_free(scr);
        SegRef r; r.k = k;
        if (read_ref(scr, sd, &r) != 0) continue;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 256;
            *refs = (SegRef *)realloc(*refs, *cap * sizeof(SegRef)); }
        (*refs)[(*n)++] = r;
    }
    closedir(d);
#endif
}

typedef struct { int k, W, H; double ou, ov; long valid, contested; int nslab; } Piece;

/* weld all wraps found under the strip roots into out_root. 0 on success. */
static int weld_run(const char *out_root, char **roots, int nroot,
                    size_t *out_nwrap, long *out_valid)
{
    Arena_T scr = Arena_new();
    SegRef *refs = NULL; size_t nseg = 0, cap = 0;
    for (int i = 0; i < nroot; i++) scan_segs(scr, roots[i], &refs, &nseg, &cap);
    if (nseg == 0) {
        fprintf(stderr, "tifxyz_weld: no segments under the %d root(s)\n", nroot);
        free(refs); Arena_free(scr); return 1;
    }
    int kmin = INT_MAX, kmax = INT_MIN;
    for (size_t i = 0; i < nseg; i++) {
        if (refs[i].k < kmin) kmin = refs[i].k;
        if (refs[i].k > kmax) kmax = refs[i].k;
    }
    fprintf(stderr, "tifxyz_weld: %zu strip segments, wraps k[%d,%d]\n",
            nseg, kmin, kmax);

    Piece *pieces = (Piece *)malloc((size_t)(kmax - kmin + 1) * sizeof(Piece));
    size_t npiece = 0; long total_valid = 0;

    for (int k = kmin; k <= kmax; k++) {
        /* members = strip segments for this wrap */
        double u0 = 1e300, u1 = -1e300, v0 = 1e300, v1 = -1e300, du = 1.0, dv = 1.0;
        int nmem = 0;
        for (size_t i = 0; i < nseg; i++) {
            if (refs[i].k != k) continue;
            SegRef *s = &refs[i];
            du = s->du; dv = s->dv;
            if (s->ou < u0) u0 = s->ou;
            if (s->ou + s->W * s->du > u1) u1 = s->ou + s->W * s->du;
            if (s->ov < v0) v0 = s->ov;
            if (s->ov + s->H * s->dv > v1) v1 = s->ov + s->H * s->dv;
            nmem++;
        }
        if (nmem == 0) continue;
        long Wk = (long)floor((u1 - u0) / du + 0.5);
        long Hk = (long)floor((v1 - v0) / dv + 0.5);
        if (Wk < 1 || Hk < 1 || (double)Wk * (double)Hk > 950000000.0) {
            fprintf(stderr, "  w%03d: skip (%ldx%ld px)\n", k, Wk, Hk);
            continue;
        }
        size_t np = (size_t)Wk * (size_t)Hk;
        float *X = (float *)malloc(np * sizeof(float));
        float *Y = (float *)malloc(np * sizeof(float));
        float *Z = (float *)malloc(np * sizeof(float));
        uint8_t *M = (uint8_t *)calloc(np, 1);
        uint8_t *P = (uint8_t *)calloc(np, 1);
        if (!X || !Y || !Z || !M || !P) { free(X);free(Y);free(Z);free(M);free(P);
            free(pieces); free(refs); Arena_free(scr); return 1; }
        for (size_t i = 0; i < np; i++) { X[i] = -1.0f; Y[i] = -1.0f; Z[i] = -1.0f; }

        for (size_t i = 0; i < nseg; i++) {
            if (refs[i].k != k) continue;
            SegRef *s = &refs[i];
            Arena_free(scr);
            char p[1300]; float *sx=NULL,*sy=NULL,*sz=NULL; int sw=0,sh=0,w2,h2;
            snprintf(p, sizeof p, "%s/x.tif", s->dir);
            if (TiffIO_load_float2d(scr, p, &sx, &sw, &sh) != 0) continue;
            snprintf(p, sizeof p, "%s/y.tif", s->dir);
            if (TiffIO_load_float2d(scr, p, &sy, &w2, &h2) != 0 || w2!=sw || h2!=sh) continue;
            snprintf(p, sizeof p, "%s/z.tif", s->dir);
            if (TiffIO_load_float2d(scr, p, &sz, &w2, &h2) != 0 || w2!=sw || h2!=sh) continue;
            uint8_t *sm=NULL,*sp=NULL; int D3=0,hh=0,ww=0;
            snprintf(p, sizeof p, "%s/mask.tif", s->dir);
            if (TiffIO_load(scr, p, &sm, &D3, &hh, &ww) != 0 || ww!=sw || hh!=sh) continue;
            snprintf(p, sizeof p, "%s/provenance.tif", s->dir);
            if (TiffIO_load(scr, p, &sp, &D3, &hh, &ww) != 0 || ww!=sw || hh!=sh) sp = NULL;
            long col0 = (long)floor((s->ou - u0) / du + 0.5);
            long row0 = (long)floor((s->ov - v0) / dv + 0.5);
            for (int r = 0; r < sh; r++) {
                long dr = r + row0;
                if (dr < 0 || dr >= Hk) continue;
                for (int c = 0; c < sw; c++) {
                    size_t si = (size_t)r * (size_t)sw + (size_t)c;
                    if (sm[si] < 255) continue;
                    long dc = c + col0;
                    if (dc < 0 || dc >= Wk) continue;
                    size_t di = (size_t)dr * (size_t)Wk + (size_t)dc;
                    if (M[di] != 0) continue;               /* first cover wins */
                    X[di] = sx[si]; Y[di] = sy[si]; Z[di] = sz[si];
                    M[di] = 255; P[di] = sp ? sp[si] : 1;
                }
            }
        }

        long nv = 0, nc = 0;
        double blo[3] = { 1e300,1e300,1e300 }, bhi[3] = { -1e300,-1e300,-1e300 };
        for (size_t i = 0; i < np; i++) {
            if (M[i] < 255) continue;
            nv++;
            if (P[i] == 3) nc++;
            if (X[i] < blo[0]) blo[0] = X[i]; if (X[i] > bhi[0]) bhi[0] = X[i];
            if (Y[i] < blo[1]) blo[1] = Y[i]; if (Y[i] > bhi[1]) bhi[1] = Y[i];
            if (Z[i] < blo[2]) blo[2] = Z[i]; if (Z[i] > bhi[2]) bhi[2] = Z[i];
        }
        if (nv == 0) { free(X);free(Y);free(Z);free(M);free(P); continue; }

        char od[1200], p[1300];
        snprintf(od, sizeof od, "%s/seg/welded_w%03d", out_root, k);
        snprintf(p, sizeof p, "%s/x.tif", od);
        if (ves_ensure_parent_dir(p) != 0) { fprintf(stderr, "mkdir %s failed\n", od); }
        int wr = 0;
        wr |= TiffIO_save_float2d(p, X, (int)Wk, (int)Hk);
        snprintf(p, sizeof p, "%s/y.tif", od); wr |= TiffIO_save_float2d(p, Y, (int)Wk, (int)Hk);
        snprintf(p, sizeof p, "%s/z.tif", od); wr |= TiffIO_save_float2d(p, Z, (int)Wk, (int)Hk);
        snprintf(p, sizeof p, "%s/mask.tif", od); wr |= TiffIO_save(p, M, 1, (int)Hk, (int)Wk);
        snprintf(p, sizeof p, "%s/provenance.tif", od); wr |= TiffIO_save(p, P, 1, (int)Hk, (int)Wk);
        snprintf(p, sizeof p, "%s/meta.json", od);
        FILE *mf = fopen(p, "w");
        if (mf) {
            fprintf(mf,
                "{\"scale\": [%.8f, %.8f], \"type\": \"seg\", \"format\": \"tifxyz\", "
                "\"uuid\": \"welded_w%03d\", \"source\": \"scrollfiesta-welded\", "
                "\"scrollfiesta\": {\"du\": %.6f, \"dv\": %.6f, \"flip_u\": 0, \"flip_v\": 0, "
                "\"origin_uv\": [%.3f, %.3f], \"wrap_k\": %d, \"n_slabs\": %d, "
                "\"valid_px\": %ld, \"contested_px\": %ld, "
                "\"provenance\": \"provenance.tif: 0=empty 1=single 2=overlap 3=contested\"}, "
                "\"bbox\": [[%.3f, %.3f, %.3f], [%.3f, %.3f, %.3f]]}\n",
                1.0/du, 1.0/dv, k, du, dv, u0, v0, k, nmem, nv, nc,
                blo[0], blo[1], blo[2], bhi[0], bhi[1], bhi[2]);
            fclose(mf);
        }
        if (wr) fprintf(stderr, "  w%03d: WARN a band write failed\n", k);

        pieces[npiece].k = k; pieces[npiece].W = (int)Wk; pieces[npiece].H = (int)Hk;
        pieces[npiece].ou = u0; pieces[npiece].ov = v0;
        pieces[npiece].valid = nv; pieces[npiece].contested = nc; pieces[npiece].nslab = nmem;
        npiece++; total_valid += nv;
        free(X); free(Y); free(Z); free(M); free(P);
    }

    /* atlas.json */
    char ap[1200]; snprintf(ap, sizeof ap, "%s/atlas.json", out_root);
    if (ves_ensure_parent_dir(ap) != 0) { /* out_root exists from seg writes */ }
    FILE *af = fopen(ap, "w");
    if (af) {
        double vlo = 1e300, vhi = -1e300;
        for (size_t i = 0; i < npiece; i++) {
            if (pieces[i].ov < vlo) vlo = pieces[i].ov;
            if (pieces[i].ov + pieces[i].H > vhi) vhi = pieces[i].ov + pieces[i].H;
        }
        fprintf(af, "{\n  \"format\": \"scrollfiesta-atlas-welded\", \"version\": 1,\n"
                    "  \"prefix\": \"welded\", \"n_pieces\": %zu,\n"
                    "  \"v_full\": [%.1f, %.1f], \"total_valid_px\": %ld,\n  \"pieces\": [\n",
                npiece, npiece ? vlo : 0.0, npiece ? vhi : 0.0, total_valid);
        for (size_t i = 0; i < npiece; i++)
            fprintf(af, "    {\"uuid\": \"welded_w%03d\", \"k\": %d, \"n_slabs\": %d, "
                        "\"W\": %d, \"H\": %d, \"origin_uv\": [%.3f, %.3f], "
                        "\"valid_px\": %ld, \"contested_px\": %ld}%s\n",
                    pieces[i].k, pieces[i].k, pieces[i].nslab, pieces[i].W, pieces[i].H,
                    pieces[i].ou, pieces[i].ov, pieces[i].valid, pieces[i].contested,
                    i + 1 < npiece ? "," : "");
        fprintf(af, "  ]\n}\n");
        fclose(af);
    }
    fprintf(stderr, "tifxyz_weld: %zu welded wraps, %ld valid px -> %s\n",
            npiece, total_valid, out_root);
    if (out_nwrap) *out_nwrap = npiece;
    if (out_valid) *out_valid = total_valid;
    free(pieces); free(refs); Arena_free(scr);
    return 0;
}

/* ============================================================================
 * Self-test: two synthetic z-strips (disjoint v bands, shared u) each with one
 * w005 segment -> weld -> one full-height welded_w005 with exact pixel parity.
 * ==========================================================================*/
static int tw_write_strip(const char *root, double ou, double ov, int W, int H,
                          float xbase)
{
    char od[512], p[640];
    snprintf(od, sizeof od, "%s/seg/tw_w005_z%05d", root, (int)ov);
    snprintf(p, sizeof p, "%s/x.tif", od);
    if (ves_ensure_parent_dir(p) != 0) return -1;
    size_t np = (size_t)W * (size_t)H;
    float *X = (float *)malloc(np*sizeof(float)), *Y = (float *)malloc(np*sizeof(float)),
          *Z = (float *)malloc(np*sizeof(float));
    uint8_t *M = (uint8_t *)malloc(np), *P = (uint8_t *)malloc(np);
    for (size_t i = 0; i < np; i++) {
        X[i] = xbase + (float)i; Y[i] = 20.0f; Z[i] = (float)(ov) + (float)(i / (size_t)W);
        M[i] = 255; P[i] = 1;
    }
    int rc = 0;
    rc |= TiffIO_save_float2d(p, X, W, H);
    snprintf(p, sizeof p, "%s/y.tif", od); rc |= TiffIO_save_float2d(p, Y, W, H);
    snprintf(p, sizeof p, "%s/z.tif", od); rc |= TiffIO_save_float2d(p, Z, W, H);
    snprintf(p, sizeof p, "%s/mask.tif", od); rc |= TiffIO_save(p, M, 1, H, W);
    snprintf(p, sizeof p, "%s/provenance.tif", od); rc |= TiffIO_save(p, P, 1, H, W);
    snprintf(p, sizeof p, "%s/meta.json", od);
    FILE *mf = fopen(p, "w");
    if (mf) { fprintf(mf, "{\"scale\":[1,1],\"scrollfiesta\":{\"du\":1.0,\"dv\":1.0,"
                          "\"flip_u\":0,\"flip_v\":0,\"origin_uv\":[%.1f,%.1f]}}\n", ou, ov);
              fclose(mf); } else rc = -1;
    free(X);free(Y);free(Z);free(M);free(P);
    return rc;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) {
        const char *base = "output/_selftest_tifxyz_weld";
        char s0[256], s1[256], out[256];
        snprintf(s0, sizeof s0, "%s/strip0", base);
        snprintf(s1, sizeof s1, "%s/strip1", base);
        snprintf(out, sizeof out, "%s/welded", base);
        int fails = 0;
        /* strip0: u[100,104) v[0,3);  strip1: same u, v[10,13) */
        if (tw_write_strip(s0, 100.0, 0.0, 4, 3, 1000.0f) != 0) { fprintf(stderr,"FAIL write s0\n"); fails++; }
        if (tw_write_strip(s1, 100.0, 10.0, 4, 3, 2000.0f) != 0) { fprintf(stderr,"FAIL write s1\n"); fails++; }
        char *roots[2] = { s0, s1 };
        size_t nw = 0; long nv = 0;
        if (weld_run(out, roots, 2, &nw, &nv) != 0) { fprintf(stderr,"FAIL weld rc\n"); fails++; }
        if (nw != 1) { fprintf(stderr, "FAIL: expected 1 welded wrap, got %zu\n", nw); fails++; }
        if (nv != 24) { fprintf(stderr, "FAIL: expected 24 valid px (2x12), got %ld\n", nv); fails++; }
        /* readback: welded_w005 must be 4 wide x 13 tall (v 0..12), strip0 rows
         * 0-2, gap rows 3-9 empty, strip1 rows 10-12. */
        Arena_T ar = Arena_new();
        char p[512]; uint8_t *M = NULL; int D=0,H=0,W=0;
        snprintf(p, sizeof p, "%s/seg/welded_w005/mask.tif", out);
        if (TiffIO_load(ar, p, &M, &D, &H, &W) != 0) { fprintf(stderr,"FAIL readback mask\n"); fails++; }
        else {
            if (W != 4 || H != 13) { fprintf(stderr, "FAIL: welded %dx%d, want 4x13\n", W, H); fails++; }
            int gap_empty = 1, s0_full = 1, s1_full = 1;
            for (int r = 0; r < H && r < 13; r++)
                for (int c = 0; c < 4; c++) {
                    uint8_t v = M[r*W+c];
                    if (r <= 2 && v != 255) s0_full = 0;
                    else if (r >= 3 && r <= 9 && v != 0) gap_empty = 0;
                    else if (r >= 10 && v != 255) s1_full = 0;
                }
            if (!s0_full) { fprintf(stderr, "FAIL: strip0 rows not full\n"); fails++; }
            if (!gap_empty) { fprintf(stderr, "FAIL: gap rows not empty\n"); fails++; }
            if (!s1_full) { fprintf(stderr, "FAIL: strip1 rows not full\n"); fails++; }
        }
        /* x value at strip1 row10 col0 must be strip1's xbase (2000), proving
         * absolute-v placement (not overwritten / mis-rowed). */
        float *X = NULL;
        snprintf(p, sizeof p, "%s/seg/welded_w005/x.tif", out);
        if (TiffIO_load_float2d(ar, p, &X, &W, &H) == 0 && H == 13) {
            if (fabs((double)X[0] - 1000.0) > 1e-3) { fprintf(stderr,"FAIL: row0 x!=1000\n"); fails++; }
            if (fabs((double)X[10*W] - 2000.0) > 1e-3) { fprintf(stderr,"FAIL: row10 x!=2000\n"); fails++; }
        } else { fprintf(stderr, "FAIL readback x\n"); fails++; }
        Arena_dispose(&ar);
        fprintf(stderr, "[tifxyz_weld selftest] %s (%d failures)\n",
                fails ? "FAILED" : "PASSED", fails);
        return fails ? 1 : 0;
    }
    if (argc < 3) {
        fprintf(stderr, "usage: tifxyz_weld <out_root> <strip_root...>\n"
                        "       tifxyz_weld --selftest\n"
                        "  each strip_root holds seg/<name>/{x,y,z,mask}.tif + meta.json;\n"
                        "  same-wrap segments (*_w<k>_*) are stacked in v onto the shared\n"
                        "  u lattice into out_root/seg/welded_w###/ + atlas.json\n");
        return 2;
    }
    return weld_run(argv[1], &argv[2], argc - 2, NULL, NULL);
}
