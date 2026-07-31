/* manifold_check.c -- standalone 2-manifold audit for a triangle-soup OBJ.
 *
 * Answers one question per mesh: is it a 2-manifold (with boundary)? Two
 * independent topology tests, both purely combinatorial -- only face indices
 * matter, so vertex coordinates are not even parsed for geometry:
 *
 *   EDGE test (lifted from grid_weld.c::manifold_audit): build directed
 *   half-edges, sort by undirected key, bucket each undirected edge by how
 *   many half-edges (= incident faces) it collects:
 *       run==1   -> boundary edge      (open hole/border; NOT a defect)
 *       run==2   -> interior edge: opposite winding = manifold pair,
 *                   SAME winding = "same_dir" (a doubled or flipped face)
 *       run>=3   -> NON-MANIFOLD edge  (3+ faces share one edge)
 *
 *   VERTEX test (catches what the edge test cannot -- the bowtie/pinch): two
 *   triangle fans meeting at a single vertex but sharing NO edge are
 *   edge-manifold yet vertex-NON-manifold. For each vertex, union its
 *   incident triangles whenever they share a vertex-incident edge; >1
 *   resulting fan => non-manifold vertex.
 *
 * VERDICT = MANIFOLD iff nm_edges==0 AND nm_verts==0. "same_dir" pairs are
 * reported separately as an orientation/duplication defect: the edge still
 * has exactly two faces (edge-manifold), but the surface is not consistently
 * oriented there (or a triangle is doubled).
 *
 * Usage:
 *   manifold_check [--brief] <mesh.obj> [<mesh.obj> ...]
 *   manifold_check [--brief] --list <paths.txt>      (one obj path per line)
 *   manifold_check --selftest
 *
 *   --brief  one parseable line per mesh, plus a corpus summary footer.
 *            Auto-enabled for --list or when more than one file is given.
 *
 * Exit: 0 = analysis ran (check the printed VERDICT -- a non-manifold mesh is
 *           a finding, not a tool error), 1 = IO error, 2 = usage error,
 *           3 = selftest fail.
 *
 * Standalone C99; no project deps (own tiny vcxproj). Modeled on
 * obj_components.c / pinhole_verdict.c.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ---- growable int array (flattened face triples), as in obj_components.c ---- */
typedef struct { int *p; size_t n, cap; } IArr;
static void iarr_push3(IArr *a, int x, int y, int z) {
    if (a->n + 3 > a->cap) { a->cap = a->cap ? a->cap * 2 : 3072;
        a->p = realloc(a->p, a->cap * sizeof(int)); }
    a->p[a->n++] = x; a->p[a->n++] = y; a->p[a->n++] = z;
}

/* ---- per-mesh audit result ---- */
typedef struct {
    long nv, nf, ndegen, nskip;
    long boundary, manifold_pairs, same_dir, nm_edges;  /* edge test */
    long nm_verts;                                       /* vertex test */
} Audit;

static int audit_is_manifold(const Audit *a) {
    return a->nm_edges == 0 && a->nm_verts == 0;
}

/* ===================================================================
 * Edge test: directed half-edges, sorted by undirected key, run lengths.
 * (cmp + walk lifted from grid_weld.c::manifold_audit.)
 * =================================================================== */
typedef struct { int32_t src, dst; } DHE;

static int cmp_dhe_undirected(const void *pa, const void *pb) {
    const DHE *a = pa, *b = pb;
    int32_t a0 = (a->src < a->dst) ? a->src : a->dst;
    int32_t a1 = (a->src < a->dst) ? a->dst : a->src;
    int32_t b0 = (b->src < b->dst) ? b->src : b->dst;
    int32_t b1 = (b->src < b->dst) ? b->dst : b->src;
    if (a0 != b0) return (a0 < b0) ? -1 : 1;
    if (a1 != b1) return (a1 < b1) ? -1 : 1;
    return 0;
}

static void edge_audit(const int *F, long nf, Audit *out) {
    if (nf == 0) return;
    size_t hn = (size_t)nf * 3;
    DHE *he = malloc(hn * sizeof(DHE));
    for (long f = 0; f < nf; f++) {
        int v0 = F[f*3+0], v1 = F[f*3+1], v2 = F[f*3+2];
        he[f*3+0].src = v0; he[f*3+0].dst = v1;
        he[f*3+1].src = v1; he[f*3+1].dst = v2;
        he[f*3+2].src = v2; he[f*3+2].dst = v0;
    }
    qsort(he, hn, sizeof(DHE), cmp_dhe_undirected);
    size_t i = 0;
    while (i < hn) {
        size_t j = i + 1;
        int32_t a0 = (he[i].src < he[i].dst) ? he[i].src : he[i].dst;
        int32_t a1 = (he[i].src < he[i].dst) ? he[i].dst : he[i].src;
        while (j < hn) {
            int32_t b0 = (he[j].src < he[j].dst) ? he[j].src : he[j].dst;
            int32_t b1 = (he[j].src < he[j].dst) ? he[j].dst : he[j].src;
            if (b0 != a0 || b1 != a1) break;
            j++;
        }
        size_t run = j - i;
        if (run == 1) { out->boundary++; }
        else if (run > 2) { out->nm_edges++; }
        else {
            int dir0 = (he[i].src   < he[i].dst)   ? 0 : 1;
            int dir1 = (he[i+1].src < he[i+1].dst) ? 0 : 1;
            if (dir0 == dir1) out->same_dir++;
            else out->manifold_pairs++;
        }
        i = j;
    }
    free(he);
}

/* ===================================================================
 * Vertex test: count triangle fans around each vertex. Two incident
 * triangles belong to the same fan iff they share a vertex-incident edge,
 * i.e. they share one of the vertex's two triangle-neighbours. >1 fan => a
 * non-manifold (bowtie/pinch) vertex.
 *
 * Implemented with a CSR of (vertex -> incident corners), where each corner
 * carries the two OTHER vertices of its triangle. Per vertex we union its
 * corners that name a common neighbour.
 * =================================================================== */
typedef struct { int nbr; int loc; } NbrPair;
static int cmp_nbrpair(const void *a, const void *b) {
    const NbrPair *x = a, *y = b;
    return (x->nbr < y->nbr) ? -1 : (x->nbr > y->nbr ? 1 : 0);
}
static int lf_find(int *p, int x) { while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; } return x; }

static void vertex_audit(long nv, const int *F, long nf, Audit *out) {
    if (nv == 0 || nf == 0) return;
    long *deg = calloc((size_t)nv, sizeof(long));
    for (long f = 0; f < nf; f++)
        for (int k = 0; k < 3; k++) deg[F[f*3+k]]++;

    /* CSR offsets over corners (3 per face) */
    long *off = calloc((size_t)nv + 1, sizeof(long));
    long maxdeg = 0;
    for (long v = 0; v < nv; v++) { off[v+1] = off[v] + deg[v]; if (deg[v] > maxdeg) maxdeg = deg[v]; }
    long total = off[nv];                       /* == 3*nf minus degenerate self-refs */
    int *cn0 = malloc((size_t)total * sizeof(int));   /* corner's first  neighbour */
    int *cn1 = malloc((size_t)total * sizeof(int));   /* corner's second neighbour */
    long *cur = malloc(((size_t)nv) * sizeof(long));
    for (long v = 0; v < nv; v++) cur[v] = off[v];
    for (long f = 0; f < nf; f++) {
        int a = F[f*3+0], b = F[f*3+1], c = F[f*3+2];
        long pa = cur[a]++; cn0[pa] = b; cn1[pa] = c;   /* at a, neighbours are b,c */
        long pb = cur[b]++; cn0[pb] = a; cn1[pb] = c;   /* at b, neighbours are a,c */
        long pc = cur[c]++; cn0[pc] = a; cn1[pc] = b;   /* at c, neighbours are a,b */
    }

    /* per-vertex fan count */
    int     *parent = malloc((size_t)(maxdeg > 0 ? maxdeg : 1) * sizeof(int));
    NbrPair *pairs  = malloc((size_t)(maxdeg > 0 ? maxdeg : 1) * 2 * sizeof(NbrPair));
    for (long v = 0; v < nv; v++) {
        long s = off[v], e = off[v+1];
        long k = e - s;
        if (k <= 1) continue;                    /* 0 or 1 incident face -> trivially one fan */
        for (long t = 0; t < k; t++) parent[t] = (int)t;
        /* collect (neighbour, local-corner) pairs, skipping self-refs (degenerate faces) */
        long np = 0;
        for (long t = 0; t < k; t++) {
            int n0 = cn0[s+t], n1 = cn1[s+t];
            if (n0 != v) { pairs[np].nbr = n0; pairs[np].loc = (int)t; np++; }
            if (n1 != v) { pairs[np].nbr = n1; pairs[np].loc = (int)t; np++; }
        }
        qsort(pairs, (size_t)np, sizeof(NbrPair), cmp_nbrpair);
        for (long t = 1; t < np; t++) {
            if (pairs[t].nbr == pairs[t-1].nbr) {
                int ra = lf_find(parent, pairs[t].loc);
                int rb = lf_find(parent, pairs[t-1].loc);
                if (ra != rb) parent[rb] = ra;
            }
        }
        long fans = 0;
        for (long t = 0; t < k; t++) if (lf_find(parent, (int)t) == (int)t) fans++;
        if (fans > 1) out->nm_verts++;
    }

    free(deg); free(off); free(cn0); free(cn1); free(cur); free(parent); free(pairs);
}

/* ===================================================================
 * OBJ loader -- topology only. Counts 'v' lines (for index bounds) and
 * fan-triangulates each 'f' polygon. Handles "f a/b/c" slash forms and
 * negative (relative) indices, exactly as obj_components.c. Faces with any
 * out-of-range index are skipped (counted in nskip); faces with a repeated
 * index are counted in ndegen but still audited (they are harmless to both
 * tests and self-refs are filtered in the vertex test).
 * =================================================================== */
static int load_obj(const char *path, long *out_nv, int **out_F, long *out_nf,
                    long *out_skip, long *out_degen) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    long nv = 0;
    IArr Fc = {0};
    char line[8192];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            nv++;                                /* coordinates irrelevant to topology */
        } else if (line[0] == 'f' && line[1] == ' ') {
            int idx[64]; int ni = 0;
            const char *p = line + 1;
            while (*p && ni < 64) {
                while (*p == ' ' || *p == '\t') p++;
                if (*p == '\0' || *p == '\n' || *p == '\r') break;
                int vi = (int)strtol(p, (char **)&p, 10);
                while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++; /* skip /vt/vn */
                if (vi < 0) vi = (int)nv + vi; else vi -= 1;   /* OBJ 1-based; neg = relative */
                idx[ni++] = vi;
            }
            for (int i = 1; i + 1 < ni; i++) iarr_push3(&Fc, idx[0], idx[i], idx[i+1]);
        }
    }
    fclose(f);

    /* validate / filter face indices in [0, nv) */
    int *F = Fc.p;
    long nf_raw = (long)(Fc.n / 3);
    long nf = 0, nskip = 0, ndegen = 0;
    for (long t = 0; t < nf_raw; t++) {
        int a = F[t*3+0], b = F[t*3+1], c = F[t*3+2];
        if (a < 0 || a >= nv || b < 0 || b >= nv || c < 0 || c >= nv) { nskip++; continue; }
        if (a == b || b == c || a == c) ndegen++;
        F[nf*3+0] = a; F[nf*3+1] = b; F[nf*3+2] = c;   /* compact in place */
        nf++;
    }
    *out_nv = nv; *out_F = F; *out_nf = nf; *out_skip = nskip; *out_degen = ndegen;
    return 0;
}

static const char *basename_of(const char *path) {
    const char *base = path;
    for (const char *s = path; *s; s++) if (*s == '/' || *s == '\\') base = s + 1;
    return base;
}

/* analyse one file; fill *out. returns 0 ok, -1 IO error. */
static int analyze_file(const char *path, Audit *out, int brief) {
    memset(out, 0, sizeof *out);
    int *F = NULL; long nv = 0, nf = 0, nskip = 0, ndegen = 0;
    if (load_obj(path, &nv, &F, &nf, &nskip, &ndegen) != 0) {
        fprintf(stderr, "ERROR: cannot open %s\n", path);
        return -1;
    }
    out->nv = nv; out->nf = nf; out->nskip = nskip; out->ndegen = ndegen;
    edge_audit(F, nf, out);
    vertex_audit(nv, F, nf, out);
    free(F);

    int man = audit_is_manifold(out);
    if (brief) {
        printf("%-44s V=%-8ld F=%-8ld bnd=%-6ld nm_edge=%-5ld nm_vert=%-5ld same_dir=%-5ld VERDICT=%s\n",
               basename_of(path), out->nv, out->nf, out->boundary,
               out->nm_edges, out->nm_verts, out->same_dir,
               man ? "MANIFOLD" : "NONMANIFOLD");
    } else {
        printf("=== manifold_check: %s ===\n", path);
        printf("  V=%ld  F=%ld", out->nv, out->nf);
        if (ndegen) printf("  (degenerate faces: %ld)", ndegen);
        if (nskip)  printf("  (skipped out-of-range faces: %ld)", nskip);
        printf("\n");
        printf("  edges:    boundary=%ld  manifold=%ld  same_dir=%ld  non_manifold=%ld\n",
               out->boundary, out->manifold_pairs, out->same_dir, out->nm_edges);
        printf("  vertices: non_manifold(bowtie/pinch)=%ld\n", out->nm_verts);
        if (man) {
            printf("  VERDICT: MANIFOLD%s\n", out->boundary ? " (with boundary)" : " (closed)");
        } else {
            printf("  VERDICT: NON-MANIFOLD  [nm_edge=%ld nm_vert=%ld]\n",
                   out->nm_edges, out->nm_verts);
        }
        if (out->same_dir)
            printf("  NOTE: %ld same_dir edge(s) -- doubled or flipped face(s); "
                   "edge-manifold but not consistently oriented.\n", out->same_dir);
    }
    return 0;
}

/* ===================================================================
 * Selftest
 * =================================================================== */
static int selftest(void) {
    int fails = 0;
    printf("=== manifold_check selftest ===\n");

    /* case 1: consistently-oriented closed tetrahedron -> fully manifold, closed */
    {
        int F[] = { 0,2,1,  0,1,3,  0,3,2,  1,2,3 };
        Audit a; memset(&a, 0, sizeof a);
        edge_audit(F, 4, &a); vertex_audit(4, F, 4, &a);
        int ok = a.boundary==0 && a.manifold_pairs==6 && a.same_dir==0 &&
                 a.nm_edges==0 && a.nm_verts==0;
        printf("[case1] closed tetra: bnd=%ld man=%ld sd=%ld nmE=%ld nmV=%ld -> %s\n",
               a.boundary, a.manifold_pairs, a.same_dir, a.nm_edges, a.nm_verts,
               ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }
    /* case 2: two triangles sharing one edge -> manifold WITH boundary */
    {
        int F[] = { 0,1,2,  1,3,2 };
        Audit a; memset(&a, 0, sizeof a);
        edge_audit(F, 2, &a); vertex_audit(4, F, 2, &a);
        int ok = a.boundary==4 && a.manifold_pairs==1 && a.nm_edges==0 &&
                 a.same_dir==0 && a.nm_verts==0;
        printf("[case2] open quad: bnd=%ld man=%ld nmE=%ld nmV=%ld -> %s\n",
               a.boundary, a.manifold_pairs, a.nm_edges, a.nm_verts, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }
    /* case 3: NON-MANIFOLD EDGE -- 3 triangles sharing edge {0,1} */
    {
        int F[] = { 0,1,2,  0,1,3,  0,1,4 };
        Audit a; memset(&a, 0, sizeof a);
        edge_audit(F, 3, &a); vertex_audit(5, F, 3, &a);
        int ok = a.nm_edges==1 && !audit_is_manifold(&a);
        printf("[case3] 3-fan edge: nmE=%ld (expect 1) -> %s\n", a.nm_edges, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }
    /* case 4: NON-MANIFOLD VERTEX (bowtie) -- two tris sharing only vertex 0.
     * Edge test sees only boundary edges; the vertex test must catch it. */
    {
        int F[] = { 0,1,2,  0,3,4 };
        Audit a; memset(&a, 0, sizeof a);
        edge_audit(F, 2, &a); vertex_audit(5, F, 2, &a);
        int ok = a.nm_edges==0 && a.nm_verts==1 && !audit_is_manifold(&a);
        printf("[case4] bowtie vertex: nmE=%ld nmV=%ld (expect 0,1) -> %s\n",
               a.nm_edges, a.nm_verts, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }
    /* case 5: doubled face -> 3 same_dir edges, edge/vertex-manifold */
    {
        int F[] = { 0,1,2,  0,1,2 };
        Audit a; memset(&a, 0, sizeof a);
        edge_audit(F, 2, &a); vertex_audit(3, F, 2, &a);
        int ok = a.same_dir==3 && a.nm_edges==0 && a.nm_verts==0;
        printf("[case5] doubled face: sd=%ld nmE=%ld nmV=%ld (expect 3,0,0) -> %s\n",
               a.same_dir, a.nm_edges, a.nm_verts, ok ? "ok" : "FAIL");
        if (!ok) fails++;
    }

    printf("=== selftest %s (%d failure%s) ===\n",
           fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}

/* ===================================================================
 * main
 * =================================================================== */
typedef struct { char name[256]; Audit a; } Row;

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) {
        fprintf(stderr,
            "Usage: %s [--brief] <mesh.obj> [<mesh.obj> ...]\n"
            "       %s [--brief] --list <paths.txt>\n"
            "       %s --selftest\n", argv[0], argv[0], argv[0]);
        return 2;
    }

    int brief = 0;
    const char *listfile = NULL;
    const char **files = malloc((size_t)argc * sizeof(char *));
    int nfiles = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--brief")) brief = 1;
        else if (!strcmp(argv[i], "--list") && i + 1 < argc) listfile = argv[++i];
        else if (argv[i][0] == '-' && argv[i][1] == '-') {
            fprintf(stderr, "unknown arg: %s\n", argv[i]); free(files); return 2;
        } else files[nfiles++] = argv[i];
    }

    /* read --list into the file set */
    char *listbuf = NULL;
    if (listfile) {
        FILE *lf = fopen(listfile, "r");
        if (!lf) { fprintf(stderr, "ERROR: cannot open list %s\n", listfile); free(files); return 1; }
        /* slurp + split lines in place */
        fseek(lf, 0, SEEK_END); long sz = ftell(lf); fseek(lf, 0, SEEK_SET);
        listbuf = malloc((size_t)sz + 1);
        size_t got = fread(listbuf, 1, (size_t)sz, lf); listbuf[got] = '\0'; fclose(lf);
        files = realloc(files, ((size_t)argc + got / 2 + 8) * sizeof(char *));
        char *p = listbuf;
        while (*p) {
            while (*p == '\r' || *p == '\n' || *p == ' ' || *p == '\t') p++;
            if (!*p) break;
            char *start = p;
            while (*p && *p != '\r' && *p != '\n') p++;
            if (*p) *p++ = '\0';
            /* trim trailing whitespace */
            char *end = start + strlen(start);
            while (end > start && (end[-1] == ' ' || end[-1] == '\t')) *--end = '\0';
            if (*start) files[nfiles++] = start;
        }
        brief = 1;
    }

    if (nfiles == 0) { fprintf(stderr, "ERROR: no input files\n"); free(files); free(listbuf); return 2; }
    if (nfiles > 1) brief = 1;

    Row *rows = malloc((size_t)nfiles * sizeof(Row));
    int nrows = 0, io_errors = 0;
    for (int i = 0; i < nfiles; i++) {
        Audit a;
        if (analyze_file(files[i], &a, brief) != 0) { io_errors++; continue; }
        snprintf(rows[nrows].name, sizeof(rows[nrows].name), "%s", basename_of(files[i]));
        rows[nrows].a = a;
        nrows++;
    }

    int rc = 0;
    if (brief) {
        long nman = 0, nnm = 0, nwind = 0;
        for (int i = 0; i < nrows; i++) {
            if (audit_is_manifold(&rows[i].a)) {
                nman++;
                if (rows[i].a.same_dir) nwind++;
            } else nnm++;
        }
        printf("\n=== summary: %d mesh(es) analysed", nrows);
        if (io_errors) printf(", %d IO error(s)", io_errors);
        printf(" ===\n");
        printf("  manifold:            %ld\n", nman);
        printf("  NON-MANIFOLD:        %ld   (nm_edge>0 or nm_vert>0)\n", nnm);
        printf("  (of the manifold, %ld also have same_dir winding/doubled defects)\n", nwind);
        if (nnm) {
            printf("\nNON-MANIFOLD meshes:\n");
            for (int i = 0; i < nrows; i++) if (!audit_is_manifold(&rows[i].a))
                printf("  %-44s nm_edge=%-5ld nm_vert=%-5ld same_dir=%-5ld\n",
                       rows[i].name, rows[i].a.nm_edges, rows[i].a.nm_verts, rows[i].a.same_dir);
            rc = 0;   /* a finding, not a tool error */
        } else {
            printf("\nAll %d mesh(es) are 2-manifold.\n", nrows);
        }
        if (nwind) {
            printf("\nWinding/doubled-face defects (same_dir>0, still edge+vertex manifold):\n");
            for (int i = 0; i < nrows; i++)
                if (audit_is_manifold(&rows[i].a) && rows[i].a.same_dir)
                    printf("  %-44s same_dir=%ld\n", rows[i].name, rows[i].a.same_dir);
        }
    }
    if (io_errors && nrows == 0) rc = 1;

    free(rows); free(files); free(listbuf);
    return rc;
}
