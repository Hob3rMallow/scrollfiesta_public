/* obj_biggest.c -- extract the largest connected component(s) from an OBJ.
 *
 * A welded grid surface is one big recto sheet plus a scatter of small "phantom"
 * components: islands BPA split off, verso flaps that papyrus separation left
 * disconnected from the recto, and stray specks. Downstream unwrapping wants ONE
 * sheet, so this tool keeps only the largest connected component(s) (triangles
 * sharing an edge) and drops the rest.
 *
 * Connectivity is pure topology -- vertex COORDINATES are never parsed. Kept
 * vertex lines are copied through byte-for-byte (so per-vertex colors, precision,
 * and the "v z y x r g b" layout survive exactly); only face indices are
 * renumbered to the compacted vertex set. Components are ranked by triangle
 * count, matching the pipeline's 1..C descending-face-count convention.
 *
 * Usage:
 *   obj_biggest <in.obj>                      # dry run: print component breakdown
 *   obj_biggest <in.obj> <out.obj> [opts]     # write the kept component(s)
 *   obj_biggest --selftest
 *
 *   --top N         keep the N largest components (default 1 = just the biggest).
 *   --min-faces F   instead keep EVERY component with >= F faces (overrides --top).
 *   --quiet         suppress the per-component table (keep the one-line summary).
 *
 * Faces are emitted as triangles (input polygons, if any, are fan-triangulated);
 * pipeline OBJs are already triangle soup so this is a verbatim pass for them.
 *
 * Exit: 0 = ok / selftest pass, 1 = IO or alloc error, 2 = usage, 3 = selftest fail.
 *
 * Standalone C99; no project deps (own tiny vcxproj). Modeled on obj_components.c.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/* growable arrays                                                      */
/* ------------------------------------------------------------------ */
typedef struct { int32_t *p; size_t n, cap; } I32Arr;
static void i32_push(I32Arr *a, int32_t v) {
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : (size_t)1 << 20;
        a->p = (int32_t *)realloc(a->p, a->cap * sizeof(int32_t));
        if (!a->p) { fprintf(stderr, "ERROR: out of memory (faces)\n"); exit(1); }
    }
    a->p[a->n++] = v;
}

/* a vertex line = a (pointer,length) slice into the file buffer, copied verbatim */
typedef struct { char *p; uint32_t len; } VSlice;
typedef struct { VSlice *p; size_t n, cap; } VArr;
static void v_push(VArr *a, char *ptr, uint32_t len) {
    if (a->n == a->cap) {
        a->cap = a->cap ? a->cap * 2 : (size_t)1 << 20;
        a->p = (VSlice *)realloc(a->p, a->cap * sizeof(VSlice));
        if (!a->p) { fprintf(stderr, "ERROR: out of memory (verts)\n"); exit(1); }
    }
    a->p[a->n].p = ptr; a->p[a->n].len = len; a->n++;
}

/* ------------------------------------------------------------------ */
/* union-find with path halving                                        */
/* ------------------------------------------------------------------ */
static int32_t uf_find(int32_t *uf, int32_t x) {
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
}

/* ------------------------------------------------------------------ */
/* connected components, labelled 1..C by descending face count        */
/* ------------------------------------------------------------------ */
typedef struct {
    int32_t *vcomp;    /* [nv]  1-based component id per vertex, 0 if unreferenced */
    long    *comp_nf;  /* [C]   triangle count of component (id-1) */
    long    *comp_nv;  /* [C]   referenced-vertex count of component (id-1) */
    int      C;        /* number of components */
} CC;

typedef struct { int32_t root; long nf, nv; } CompTmp;
static int comptmp_cmp(const void *A, const void *B) {
    const CompTmp *a = A, *b = B;
    if (a->nf != b->nf) return a->nf < b->nf ? 1 : -1;   /* descending faces */
    if (a->nv != b->nv) return a->nv < b->nv ? 1 : -1;   /* tie: descending verts */
    return a->root < b->root ? -1 : (a->root > b->root ? 1 : 0);  /* stable */
}

static void cc_free(CC *cc) {
    free(cc->vcomp); free(cc->comp_nf); free(cc->comp_nv);
    cc->vcomp = NULL; cc->comp_nf = cc->comp_nv = NULL; cc->C = 0;
}

/* Build components over triangles `tris` (0-based, 3 indices each, all in [0,nv)).
 * Returns 0 on success, -1 on alloc failure. On success `cc` owns three arrays. */
static int cc_compute(size_t nv, const int32_t *tris, size_t ntri, CC *cc)
{
    cc->vcomp = NULL; cc->comp_nf = cc->comp_nv = NULL; cc->C = 0;
    if (nv == 0) return 0;

    int32_t *uf = (int32_t *)malloc(nv * sizeof(int32_t));
    char    *ref = (char *)calloc(nv, 1);
    long *rootnf = (long *)calloc(nv, sizeof(long));
    long *rootnv = (long *)calloc(nv, sizeof(long));
    if (!uf || !ref || !rootnf || !rootnv) { free(uf); free(ref); free(rootnf); free(rootnv); return -1; }
    for (size_t i = 0; i < nv; i++) uf[i] = (int32_t)i;

    for (size_t t = 0; t < ntri; t++) {
        int32_t a = tris[3*t], b = tris[3*t+1], c = tris[3*t+2];
        int32_t ra = uf_find(uf, a), rb = uf_find(uf, b);
        if (ra != rb) uf[rb] = ra;
        ra = uf_find(uf, a); int32_t rc = uf_find(uf, c);
        if (ra != rc) uf[rc] = ra;
    }
    for (size_t i = 0; i < 3*ntri; i++) ref[tris[i]] = 1;
    for (size_t t = 0; t < ntri; t++) rootnf[uf_find(uf, tris[3*t])]++;
    for (size_t i = 0; i < nv; i++) if (ref[i]) rootnv[uf_find(uf, (int32_t)i)]++;

    int C = 0;
    for (size_t i = 0; i < nv; i++) if (ref[i] && uf_find(uf, (int32_t)i) == (int32_t)i) C++;

    CompTmp *tmp = (CompTmp *)malloc((size_t)(C ? C : 1) * sizeof(CompTmp));
    int32_t *root2id = (int32_t *)calloc(nv, sizeof(int32_t));
    int32_t *vcomp = (int32_t *)calloc(nv, sizeof(int32_t));
    long *comp_nf = (long *)malloc((size_t)(C ? C : 1) * sizeof(long));
    long *comp_nv = (long *)malloc((size_t)(C ? C : 1) * sizeof(long));
    if (!tmp || !root2id || !vcomp || !comp_nf || !comp_nv) {
        free(uf); free(ref); free(rootnf); free(rootnv);
        free(tmp); free(root2id); free(vcomp); free(comp_nf); free(comp_nv);
        return -1;
    }
    int ci = 0;
    for (size_t i = 0; i < nv; i++) if (ref[i] && uf_find(uf, (int32_t)i) == (int32_t)i) {
        tmp[ci].root = (int32_t)i; tmp[ci].nf = rootnf[i]; tmp[ci].nv = rootnv[i]; ci++;
    }
    qsort(tmp, (size_t)C, sizeof(CompTmp), comptmp_cmp);
    for (int k = 0; k < C; k++) { root2id[tmp[k].root] = k + 1; comp_nf[k] = tmp[k].nf; comp_nv[k] = tmp[k].nv; }
    for (size_t i = 0; i < nv; i++) if (ref[i]) vcomp[i] = root2id[uf_find(uf, (int32_t)i)];

    free(uf); free(ref); free(rootnf); free(rootnv); free(tmp); free(root2id);
    cc->vcomp = vcomp; cc->comp_nf = comp_nf; cc->comp_nv = comp_nv; cc->C = C;
    return 0;
}

/* ------------------------------------------------------------------ */
/* parse an OBJ into vertex-line slices + resolved 0-based triangles   */
/* ------------------------------------------------------------------ */
typedef struct {
    char   *buf;     /* owns the whole file; VSlice pointers alias into it */
    VArr    verts;   /* one slice per "v " line, in file order */
    I32Arr  tris;    /* 3 ids per triangle, 0-based, all valid in [0,nv) */
    long    n_other; /* vt/vn/vp lines skipped */
    long    n_badface; /* triangles dropped for out-of-range indices */
} ObjData;

static int obj_parse(const char *path, ObjData *od)
{
    memset(od, 0, sizeof(*od));
    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "ERROR: cannot open %s\n", path); return 1; }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz < 0) { fclose(f); fprintf(stderr, "ERROR: ftell failed on %s\n", path); return 1; }
    char *buf = (char *)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); fprintf(stderr, "ERROR: out of memory reading %s\n", path); return 1; }
    size_t rd = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    buf[rd] = '\0';
    od->buf = buf;

    /* signed original OBJ indices, resolved after nv is final (handles negatives) */
    I32Arr raw = {0};
    char *p = buf, *end = buf + rd;
    while (p < end) {
        char *nl = (char *)memchr(p, '\n', (size_t)(end - p));
        char *ls = p, *le = nl ? nl : end;           /* line = [ls,le) exclusive */
        char *te = le;
        if (te > ls && te[-1] == '\r') te--;          /* trim CR */
        size_t llen = (size_t)(te - ls);

        if (llen >= 2 && ls[0] == 'v' && (ls[1] == ' ' || ls[1] == '\t')) {
            v_push(&od->verts, ls, (uint32_t)llen);
        } else if (llen >= 2 && ls[0] == 'f' && (ls[1] == ' ' || ls[1] == '\t')) {
            int32_t idx[256]; int cnt = 0;
            char *q = ls + 1;
            while (q < te && cnt < 256) {
                while (q < te && (*q == ' ' || *q == '\t')) q++;
                if (q >= te) break;
                char *ep = q;
                long v = strtol(q, &ep, 10);          /* stops at '/', space, EOL */
                if (ep == q) break;                   /* not a number */
                idx[cnt++] = (int32_t)v;
                q = ep;
                while (q < te && *q != ' ' && *q != '\t') q++;   /* skip /vt/vn */
            }
            for (int i = 1; i + 1 < cnt; i++) {       /* fan triangulate */
                i32_push(&raw, idx[0]); i32_push(&raw, idx[i]); i32_push(&raw, idx[i+1]);
            }
        } else if (llen >= 2 && ls[0] == 'v' &&
                   (ls[1] == 't' || ls[1] == 'n' || ls[1] == 'p')) {
            od->n_other++;
        }
        p = nl ? nl + 1 : end;
    }

    /* resolve signed -> 0-based; drop triangles with any out-of-range corner */
    size_t nv = od->verts.n;
    size_t nraw = raw.n / 3;
    for (size_t t = 0; t < nraw; t++) {
        int32_t r[3]; int ok = 1;
        for (int k = 0; k < 3; k++) {
            long s = raw.p[3*t + k];
            long v = s > 0 ? s - 1 : (long)nv + s;     /* s<0 = relative; s==0 invalid */
            if (v < 0 || (size_t)v >= nv) { ok = 0; break; }
            r[k] = (int32_t)v;
        }
        if (ok) { i32_push(&od->tris, r[0]); i32_push(&od->tris, r[1]); i32_push(&od->tris, r[2]); }
        else od->n_badface++;
    }
    free(raw.p);
    return 0;
}

static void obj_free(ObjData *od) {
    free(od->buf); free(od->verts.p); free(od->tris.p);
    memset(od, 0, sizeof(*od));
}

/* ------------------------------------------------------------------ */
/* run: parse -> components -> keep top-N (or >= min-faces) -> write    */
/* ------------------------------------------------------------------ */
typedef struct {
    int  C;            /* total components */
    int  kept_comps;   /* components written */
    long total_nv, total_nf;
    long kept_nv, kept_nf;
    long largest_nf;   /* faces in the single biggest component (0 if none) */
} Stats;

static const char *base_name(const char *path) {
    const char *b = path;
    for (const char *s = path; *s; s++) if (*s == '/' || *s == '\\') b = s + 1;
    return b;
}

/* out == NULL => dry run (analyse + report only, no file written). */
static int run_file(const char *in, const char *out, int top, long min_faces,
                    int verbose, Stats *st)
{
    memset(st, 0, sizeof(*st));
    ObjData od;
    int rc = obj_parse(in, &od);
    if (rc) return rc;

    size_t nv = od.verts.n, ntri = od.tris.n / 3;
    if (nv > (size_t)INT32_MAX) { fprintf(stderr, "ERROR: %zu verts exceeds int32 index range\n", nv); obj_free(&od); return 1; }

    CC cc;
    if (cc_compute(nv, od.tris.p, ntri, &cc) != 0) {
        fprintf(stderr, "ERROR: out of memory during component labelling\n");
        obj_free(&od); return 1;
    }
    st->C = cc.C;
    st->total_nv = (long)nv;
    st->total_nf = (long)ntri;
    st->largest_nf = cc.C ? cc.comp_nf[0] : 0;

    if (verbose) {
        printf("=== obj_biggest: %s ===\n", base_name(in));
        printf("  V=%zu  F=%zu  components=%d", nv, ntri, cc.C);
        if (od.n_other)   printf("  (skipped %ld vt/vn/vp)", od.n_other);
        if (od.n_badface) printf("  (dropped %ld out-of-range face(s))", od.n_badface);
        printf("\n");
        int show = cc.C < 12 ? cc.C : 12;
        printf("    rank     faces     verts   %%faces\n");
        for (int k = 0; k < show; k++)
            printf("   %5d %9ld %9ld   %5.1f%%\n", k + 1, cc.comp_nf[k], cc.comp_nv[k],
                   ntri ? 100.0 * (double)cc.comp_nf[k] / (double)ntri : 0.0);
        if (cc.C > show) printf("    ... %d smaller component(s)\n", cc.C - show);
    }

    /* selection: min_faces>0 keeps every comp with >= min_faces; else the top N */
    int kept_comps = 0; long kept_nv = 0, kept_nf = 0;
    char *keep = (char *)calloc((size_t)(cc.C ? cc.C : 1), 1);
    if (!keep) { cc_free(&cc); obj_free(&od); fprintf(stderr, "ERROR: out of memory\n"); return 1; }
    for (int k = 0; k < cc.C; k++) {
        int sel = (min_faces > 0) ? (cc.comp_nf[k] >= min_faces) : (k < top);
        if (sel) { keep[k] = 1; kept_comps++; kept_nv += cc.comp_nv[k]; kept_nf += cc.comp_nf[k]; }
    }
    st->kept_comps = kept_comps;
    st->kept_nv = kept_nv;
    st->kept_nf = kept_nf;

    if (verbose) {
        if (min_faces > 0)
            printf("  keeping %d component(s) with >= %ld faces:", kept_comps, min_faces);
        else
            printf("  keeping top %d component(s):", kept_comps);
        printf(" %ld / %zu faces (%.1f%%), %ld / %zu verts\n",
               kept_nf, ntri, ntri ? 100.0 * (double)kept_nf / (double)ntri : 0.0, kept_nv, nv);
    }

    if (!out) { free(keep); cc_free(&cc); obj_free(&od); return 0; }

    FILE *of = fopen(out, "wb");
    if (!of) { fprintf(stderr, "ERROR: cannot write %s\n", out); free(keep); cc_free(&cc); obj_free(&od); return 1; }
    static char iobuf[1 << 24];              /* 16 MB write buffer */
    setvbuf(of, iobuf, _IOFBF, sizeof(iobuf));

    fprintf(of, "# obj_biggest: kept %d of %d component(s) from %s\n", kept_comps, cc.C, base_name(in));
    if (min_faces > 0) fprintf(of, "# selection: components with >= %ld faces\n", min_faces);
    else               fprintf(of, "# selection: top %d component(s) by face count\n", top);
    fprintf(of, "# %ld verts, %ld faces\n", kept_nv, kept_nf);

    /* compact vertices: new 1-based id in original order for kept components */
    int32_t *newid = (int32_t *)malloc(nv * sizeof(int32_t));
    if (!newid) { fclose(of); free(keep); cc_free(&cc); obj_free(&od); fprintf(stderr, "ERROR: out of memory\n"); return 1; }
    long out_nv = 0;
    for (size_t i = 0; i < nv; i++) {
        int32_t comp = cc.vcomp[i];
        if (comp && keep[comp-1]) {
            newid[i] = (int32_t)(++out_nv);
            fwrite(od.verts.p[i].p, 1, od.verts.p[i].len, of);
            fputc('\n', of);
        } else newid[i] = 0;
    }
    long out_nf = 0;
    for (size_t t = 0; t < ntri; t++) {
        int32_t a = od.tris.p[3*t];
        int32_t comp = cc.vcomp[a];
        if (!comp || !keep[comp-1]) continue;          /* all 3 share a component */
        fprintf(of, "f %d %d %d\n",
                (int)newid[a], (int)newid[od.tris.p[3*t+1]], (int)newid[od.tris.p[3*t+2]]);
        out_nf++;
    }
    if (fclose(of) != 0) { fprintf(stderr, "ERROR: write failed on %s\n", out); free(newid); free(keep); cc_free(&cc); obj_free(&od); return 1; }

    if (verbose) printf("  wrote %s  (%ld verts, %ld faces)\n", out, out_nv, out_nf);
    st->kept_nv = out_nv; st->kept_nf = out_nf;

    free(newid); free(keep); cc_free(&cc); obj_free(&od);
    return 0;
}

/* ================================================================== */
/* self-test                                                           */
/* ================================================================== */
static int write_tmp(const char *path, const char *content) {
    FILE *f = fopen(path, "wb");
    if (!f) return 1;
    fputs(content, f);
    fclose(f);
    return 0;
}
static int count_lines(const char *path, long *nv, long *nf) {
    *nv = *nf = 0;
    FILE *f = fopen(path, "rb");
    if (!f) return 1;
    char ln[512];
    while (fgets(ln, sizeof(ln), f)) {
        if (ln[0]=='v' && ln[1]==' ') (*nv)++;
        else if (ln[0]=='f' && ln[1]==' ') (*nf)++;
    }
    fclose(f);
    return 0;
}

static int selftest(void)
{
    int fails = 0;
    printf("=== obj_biggest selftest ===\n");

    /* case 1: pure cc_compute -- two triangles sharing an edge => 1 component */
    {
        int32_t T[] = { 0,1,2,  1,3,2 };
        CC cc; int rc = cc_compute(4, T, 2, &cc);
        printf("[case1] shared-edge quad (expect C=1, nf=2)\n");
        if (rc || cc.C != 1 || cc.comp_nf[0] != 2) { printf("  FAIL: rc=%d C=%d nf=%ld\n", rc, cc.C, cc.C?cc.comp_nf[0]:-1); fails++; }
        else printf("  ok: C=1 nf=2\n");
        cc_free(&cc);
    }

    /* case 2: cc_compute ranking -- a 3-tri comp and a 1-tri comp, verts disjoint */
    {
        /* comp A: verts 0..3, 3 tris (a strip); comp B: verts 4..6, 1 tri */
        int32_t T[] = { 0,1,2,  1,2,3,  0,2,3,   4,5,6 };
        CC cc; int rc = cc_compute(7, T, 4, &cc);
        printf("[case2] 3-tri + 1-tri (expect C=2, comp1 nf=3, comp2 nf=1)\n");
        if (rc || cc.C != 2 || cc.comp_nf[0] != 3 || cc.comp_nf[1] != 1) {
            printf("  FAIL: rc=%d C=%d nf=[%ld,%ld]\n", rc, cc.C, cc.C>0?cc.comp_nf[0]:-1, cc.C>1?cc.comp_nf[1]:-1); fails++;
        } else printf("  ok: C=2 ranked 3,1\n");
        cc_free(&cc);
    }

    /* case 3: end-to-end -- big comp (3 tris) + small comp (1 tri, colored vert),
     * keep top 1 => small comp's faces gone; big comp's 4 verts + 3 faces remain */
    {
        const char *in  = "obj_biggest_selftest_in.obj";
        const char *out = "obj_biggest_selftest_out.obj";
        const char *content =
            "# test\n"
            "v 0 0 0\n" "v 1 0 0\n" "v 0 1 0\n" "v 1 1 0\n"   /* big comp verts 1..4 */
            "v 9 9 9 0.5 0.6 0.7\n" "v 9 9 8\n" "v 9 8 9\n"    /* small comp verts 5..7 */
            "f 1 2 3\n" "f 2 4 3\n" "f 1 3 4\n"                /* big: 3 tris */
            "f 5 6 7\n";                                        /* small: 1 tri */
        printf("[case3] end-to-end top-1 (expect out nv=4 nf=3)\n");
        Stats st;
        if (write_tmp(in, content)) { printf("  FAIL: cannot write tmp\n"); fails++; }
        else if (run_file(in, out, 1, 0, 0, &st)) { printf("  FAIL: run_file error\n"); fails++; }
        else {
            long nv=0, nf=0; count_lines(out, &nv, &nf);
            if (!(st.C==2 && nv==4 && nf==3 && st.kept_nf==3)) {
                printf("  FAIL: C=%d out nv=%ld nf=%ld (kept_nf=%ld)\n", st.C, nv, nf, st.kept_nf); fails++;
            } else printf("  ok: C=2, kept biggest -> nv=4 nf=3\n");
        }
        remove(in); remove(out);
    }

    /* case 4: --top 2 keeps both comps (nv=7 nf=4) */
    {
        const char *in  = "obj_biggest_selftest_in.obj";
        const char *out = "obj_biggest_selftest_out.obj";
        const char *content =
            "v 0 0 0\n" "v 1 0 0\n" "v 0 1 0\n" "v 1 1 0\n"
            "v 9 9 9\n" "v 9 9 8\n" "v 9 8 9\n"
            "f 1 2 3\n" "f 2 4 3\n" "f 1 3 4\n" "f 5 6 7\n";
        printf("[case4] --top 2 (expect out nv=7 nf=4)\n");
        Stats st;
        if (write_tmp(in, content) || run_file(in, out, 2, 0, 0, &st)) { printf("  FAIL: run error\n"); fails++; }
        else {
            long nv=0, nf=0; count_lines(out, &nv, &nf);
            if (!(nv==7 && nf==4)) { printf("  FAIL: nv=%ld nf=%ld\n", nv, nf); fails++; }
            else printf("  ok: both kept nv=7 nf=4\n");
        }
        remove(in); remove(out);
    }

    /* case 5: color + coordinates on the kept vertex survive verbatim */
    {
        const char *in  = "obj_biggest_selftest_in.obj";
        const char *out = "obj_biggest_selftest_out.obj";
        const char *content =
            "v 4357.412598 3079.225830 2629.845703 0.9880 0.5530 0.3840\n"
            "v 1 0 0\n" "v 0 1 0\n"
            "v 9 9 9\n" "v 9 9 8\n" "v 9 8 9\n"
            "f 1 2 3\n"                       /* big comp (rank1: 1 tri, lower id) */
            "f 4 5 6\n";                      /* other comp (1 tri) */
        printf("[case5] verbatim color passthrough on kept vert\n");
        Stats st;
        int ok = 0;
        if (!write_tmp(in, content) && !run_file(in, out, 1, 0, 0, &st)) {
            FILE *f = fopen(out, "rb");
            if (f) {
                char ln[512]; int found = 0;
                while (fgets(ln, sizeof(ln), f))
                    if (strstr(ln, "4357.412598 3079.225830 2629.845703 0.9880 0.5530 0.3840")) found = 1;
                fclose(f);
                ok = found;
            }
        }
        if (!ok) { printf("  FAIL: colored vertex line not preserved verbatim\n"); fails++; }
        else printf("  ok: colored coords preserved exactly\n");
        remove(in); remove(out);
    }

    /* case 6: negative (relative) face indices resolve correctly */
    {
        const char *in  = "obj_biggest_selftest_in.obj";
        const char *out = "obj_biggest_selftest_out.obj";
        const char *content =
            "v 0 0 0\n" "v 1 0 0\n" "v 0 1 0\n"
            "f -3 -2 -1\n";                   /* == f 1 2 3 */
        printf("[case6] negative face indices (expect nv=3 nf=1)\n");
        Stats st;
        if (write_tmp(in, content) || run_file(in, out, 1, 0, 0, &st)) { printf("  FAIL: run error\n"); fails++; }
        else {
            long nv=0, nf=0; count_lines(out, &nv, &nf);
            if (!(nv==3 && nf==1 && st.C==1)) { printf("  FAIL: C=%d nv=%ld nf=%ld\n", st.C, nv, nf); fails++; }
            else printf("  ok: relative indices resolved, nv=3 nf=1\n");
        }
        remove(in); remove(out);
    }

    /* case 7: --min-faces threshold keeps only comps at/above the bar */
    {
        int32_t T[] = { 0,1,2,  1,2,3,  0,2,3,   4,5,6 };  /* nf 3 and 1 */
        CC cc; cc_compute(7, T, 4, &cc);
        printf("[case7] min-faces logic sanity (comp nf = 3,1)\n");
        int keep_at_2 = 0; for (int k=0;k<cc.C;k++) if (cc.comp_nf[k] >= 2) keep_at_2++;
        if (keep_at_2 != 1) { printf("  FAIL: %d comps >=2 faces (expect 1)\n", keep_at_2); fails++; }
        else printf("  ok: exactly 1 comp has >=2 faces\n");
        cc_free(&cc);
    }

    /* case 8: empty / no-faces input => C=0, graceful */
    {
        CC cc; int rc = cc_compute(0, NULL, 0, &cc);
        printf("[case8] empty mesh (expect C=0, no crash)\n");
        if (rc || cc.C != 0) { printf("  FAIL: rc=%d C=%d\n", rc, cc.C); fails++; }
        else printf("  ok: C=0\n");
        cc_free(&cc);
    }

    printf("=== selftest %s (%d failure%s) ===\n",
           fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}

/* ================================================================== */
int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s <in.obj> [out.obj] [--top N] [--min-faces F] [--quiet]\n"
            "       %s --selftest\n"
            "  With no out.obj: prints the component breakdown (dry run).\n"
            "  --top N        keep the N largest components (default 1).\n"
            "  --min-faces F  keep every component with >= F faces (overrides --top).\n",
            argv[0], argv[0]);
        return 2;
    }
    const char *in = NULL, *out = NULL;
    int top = 1, quiet = 0;
    long min_faces = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--top") && i+1 < argc) top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min-faces") && i+1 < argc) min_faces = atol(argv[++i]);
        else if (!strcmp(argv[i], "--quiet")) quiet = 1;
        else if (argv[i][0] == '-' && argv[i][1] == '-') { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 2; }
        else if (!in) in = argv[i];
        else if (!out) out = argv[i];
        else { fprintf(stderr, "unexpected arg: %s\n", argv[i]); return 2; }
    }
    if (!in) { fprintf(stderr, "ERROR: no input OBJ\n"); return 2; }
    if (top < 1) top = 1;

    clock_t t0 = clock();
    Stats st;
    int rc = run_file(in, out, top, min_faces, !quiet, &st);
    if (rc) return rc;
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    if (!quiet) printf("  done in %.2f s\n", secs);
    return 0;
}
