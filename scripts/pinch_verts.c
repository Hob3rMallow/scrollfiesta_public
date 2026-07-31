/* Detect non-manifold ("pinch"/"bowtie") vertices in an OBJ mesh.
 *
 * A vertex is manifold iff its incident triangles form a single fan
 * (connected via shared edges through that vertex). If the incident
 * faces split into >=2 fan components, the vertex is a pinch/bowtie:
 * two sheet corners touching at one point. Every edge can still be
 * <=2 faces (edge-manifold) while the vertex is non-manifold, so the
 * usual edge-based audit misses these.
 *
 * Also reports boundary-edge degree per vertex (a manifold boundary
 * vertex has exactly 2 incident boundary edges; >2 => pinch).
 *
 * Plain C99, single file, reads "v z y x" and "f a b c" (1-based). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct { float z, y, x; } V3;

static int *uf;
static int uf_find(int a) { while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; } return a; }
static void uf_union(int a, int b) { a = uf_find(a); b = uf_find(b); if (a != b) uf[a] = b; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s mesh.obj [max_report=10]\n", argv[0]); return 1; }
    int max_report = (argc >= 3) ? atoi(argv[2]) : 10;
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }

    int vcap = 1 << 16, nv = 0, fcap = 1 << 16, nf = 0;
    V3 *V = malloc((size_t)vcap * sizeof(V3));
    int *F = malloc((size_t)fcap * 3 * sizeof(int));
    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            V3 v;
            if (sscanf(line, "v %f %f %f", &v.z, &v.y, &v.x) == 3) {
                if (nv >= vcap) { vcap *= 2; V = realloc(V, (size_t)vcap * sizeof(V3)); }
                V[nv++] = v;
            }
        } else if (line[0] == 'f' && line[1] == ' ') {
            int a, b, c;
            if (sscanf(line, "f %d %d %d", &a, &b, &c) == 3) {
                if (nf >= fcap) { fcap *= 2; F = realloc(F, (size_t)fcap * 3 * sizeof(int)); }
                F[nf*3+0] = a-1; F[nf*3+1] = b-1; F[nf*3+2] = c-1; nf++;
            }
        }
    }
    fclose(f);
    fprintf(stderr, "loaded %d v %d f\n", nv, nf);

    /* vertex -> incident faces (CSR) */
    int *deg = calloc((size_t)nv, sizeof(int));
    for (int i = 0; i < nf; i++)
        for (int k = 0; k < 3; k++) deg[F[i*3+k]]++;
    int *off = malloc((size_t)(nv+1) * sizeof(int));
    off[0] = 0;
    for (int i = 0; i < nv; i++) off[i+1] = off[i] + deg[i];
    int total = off[nv];
    int *inc = malloc((size_t)total * sizeof(int)); /* incident face ids */
    int *cur = malloc((size_t)nv * sizeof(int));
    memcpy(cur, off, (size_t)nv * sizeof(int));
    for (int i = 0; i < nf; i++)
        for (int k = 0; k < 3; k++) { int v = F[i*3+k]; inc[cur[v]++] = i; }

    /* boundary-edge degree per vertex: build edge->face count via hash */
    #define HSZ (1u<<23)
    int *bk = malloc(HSZ * sizeof(int));
    for (size_t i = 0; i < HSZ; i++) bk[i] = -1;
    typedef struct { int a, b, cnt, next; } E;
    E *edges = malloc((size_t)nf * 3 * sizeof(E));
    int en = 0;
    for (int i = 0; i < nf; i++) {
        int v[3] = { F[i*3], F[i*3+1], F[i*3+2] };
        for (int e = 0; e < 3; e++) {
            int a = v[e], b = v[(e+1)%3];
            if (a > b) { int t=a; a=b; b=t; }
            uint64_t h = (uint64_t)(uint32_t)a * 0x9E3779B97F4A7C15ULL + (uint64_t)(uint32_t)b;
            h ^= h>>33; h *= 0xff51afd7ed558ccdULL; h ^= h>>33;
            size_t bb = h & (HSZ-1);
            int idx = bk[bb];
            while (idx >= 0) { if (edges[idx].a==a && edges[idx].b==b) break; idx = edges[idx].next; }
            if (idx < 0) { idx = en++; edges[idx].a=a; edges[idx].b=b; edges[idx].cnt=0; edges[idx].next=bk[bb]; bk[bb]=idx; }
            edges[idx].cnt++;
        }
    }
    int *bdeg = calloc((size_t)nv, sizeof(int)); /* incident boundary-edge count */
    for (int i = 0; i < en; i++) {
        if (edges[i].cnt == 1) { bdeg[edges[i].a]++; bdeg[edges[i].b]++; }
    }

    /* fan-component count per vertex via union-find over incident faces.
     * Two incident faces are in the same fan iff they share an edge that
     * contains v. */
    uf = malloc((size_t)nf * sizeof(int));
    int n_pinch = 0, n_bdry = 0, n_bdry_pinch = 0, reported = 0;
    int hist_fans[16] = {0};
    for (int v = 0; v < nv; v++) {
        int s = off[v], e = off[v+1], d = e - s;
        if (d <= 0) continue;
        for (int j = s; j < e; j++) uf[inc[j]] = inc[j];
        /* For each incident face, its two edges through v go to the other
         * two verts. Union faces sharing the same "other endpoint". */
        for (int j = s; j < e; j++) {
            int fi = inc[j];
            int o1 = -1, o2 = -1;
            for (int k = 0; k < 3; k++) { int w = F[fi*3+k]; if (w != v) { if (o1<0) o1=w; else o2=w; } }
            for (int j2 = j+1; j2 < e; j2++) {
                int gj = inc[j2];
                int p1=-1,p2=-1;
                for (int k=0;k<3;k++){ int w=F[gj*3+k]; if(w!=v){ if(p1<0)p1=w; else p2=w; } }
                if (o1==p1||o1==p2||o2==p1||o2==p2) uf_union(fi, gj);
            }
        }
        int fans = 0;
        for (int j = s; j < e; j++) if (uf_find(inc[j]) == inc[j]) fans++;
        if (fans < 16) hist_fans[fans]++;
        int is_bdry = (bdeg[v] > 0);
        if (is_bdry) n_bdry++;
        if (fans >= 2) {
            n_pinch++;
            if (is_bdry) n_bdry_pinch++;
            if (reported < max_report) {
                printf("PINCH v%d at (%.3f,%.3f,%.3f): %d incident faces, %d fans, bdeg=%d\n",
                       v, (double)V[v].z, (double)V[v].y, (double)V[v].x, d, fans, bdeg[v]);
                reported++;
            }
        }
    }

    /* boundary-degree histogram */
    int bd_hist[16] = {0}, bd_over = 0;
    for (int v = 0; v < nv; v++) {
        int bd = bdeg[v];
        if (bd == 0) continue;
        if (bd < 16) bd_hist[bd]++; else bd_over++;
    }

    printf("\n=== fan-component histogram (incident-face fans per vertex) ===\n");
    for (int i = 0; i < 16; i++) if (hist_fans[i]) printf("  %2d fans: %d verts\n", i, hist_fans[i]);
    printf("\n=== boundary-edge-degree histogram (boundary verts only) ===\n");
    for (int i = 1; i < 16; i++) if (bd_hist[i]) printf("  bdeg %2d: %d verts%s\n", i, bd_hist[i], (i==2)?"  (manifold boundary)":(i>2?"  (PINCH)":""));
    if (bd_over) printf("  bdeg>=16: %d verts\n", bd_over);

    printf("\n=== SUMMARY ===\n");
    printf("  vertices:                 %d\n", nv);
    printf("  boundary vertices:        %d\n", n_bdry);
    printf("  NON-MANIFOLD vertices:    %d  (fans>=2; pinch/bowtie)\n", n_pinch);
    printf("    of which on boundary:   %d\n", n_bdry_pinch);
    return 0;
}
