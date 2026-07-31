/* Print details about non-manifold edges in a mesh. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define HASH_SZ (1u<<22)
typedef struct { int va, vb, count, faces[8], next; } Edge;
typedef struct { float z, y, x; } V3;

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s mesh.obj [max=20]\n", argv[0]); return 1; }
    int max = (argc >= 3) ? atoi(argv[2]) : 20;
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    int vcap = 1<<16, nv = 0, fcap = 1<<16, nf = 0;
    V3 *V = malloc(vcap*sizeof(V3));
    int *F = malloc(fcap*3*sizeof(int));
    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (line[0]=='v' && line[1]==' ') {
            V3 v;
            if (sscanf(line, "v %f %f %f", &v.z, &v.y, &v.x) == 3) {
                if (nv>=vcap) { vcap*=2; V = realloc(V, vcap*sizeof(V3)); }
                V[nv++] = v;
            }
        } else if (line[0]=='f' && line[1]==' ') {
            int a, b, c;
            if (sscanf(line, "f %d %d %d", &a, &b, &c) == 3) {
                if (nf>=fcap) { fcap*=2; F = realloc(F, fcap*3*sizeof(int)); }
                F[nf*3+0]=a-1; F[nf*3+1]=b-1; F[nf*3+2]=c-1; nf++;
            }
        }
    }
    fclose(f);
    fprintf(stderr, "loaded %d v %d f\n", nv, nf);

    int *bk = malloc(HASH_SZ*sizeof(int));
    for (size_t i=0; i<HASH_SZ; i++) bk[i] = -1;
    Edge *E = malloc(nf*3*sizeof(Edge));
    int en = 0;
    for (int i=0; i<nf; i++) {
        int v[3] = { F[i*3], F[i*3+1], F[i*3+2] };
        for (int e=0; e<3; e++) {
            int va = v[e], vb = v[(e+1)%3];
            if (va > vb) { int t=va; va=vb; vb=t; }
            uint64_t h = ((uint64_t)(uint32_t)va) * 0x9E3779B97F4A7C15ULL + (uint64_t)(uint32_t)vb;
            h ^= h>>33; h *= 0xff51afd7ed558ccdULL; h ^= h>>33;
            size_t b = h & (HASH_SZ-1);
            int idx = bk[b];
            while (idx >= 0) {
                if (E[idx].va==va && E[idx].vb==vb) break;
                idx = E[idx].next;
            }
            if (idx < 0) {
                idx = en++;
                E[idx].va=va; E[idx].vb=vb; E[idx].count=0;
                E[idx].next = bk[b]; bk[b]=idx;
            }
            if (E[idx].count < 8) E[idx].faces[E[idx].count] = i;
            E[idx].count++;
        }
    }
    int n_nm = 0;
    for (int i=0; i<en; i++) {
        if (E[i].count < 3) continue;
        n_nm++;
        if (n_nm > max) continue;
        printf("\n=== NM edge (%d,%d) count=%d ===\n", E[i].va, E[i].vb, E[i].count);
        printf("  va: %.4f %.4f %.4f\n", V[E[i].va].z, V[E[i].va].y, V[E[i].va].x);
        printf("  vb: %.4f %.4f %.4f\n", V[E[i].vb].z, V[E[i].vb].y, V[E[i].vb].x);
        for (int k=0; k<E[i].count && k<8; k++) {
            int fi = E[i].faces[k];
            int a=F[fi*3+0], b=F[fi*3+1], c=F[fi*3+2];
            /* find third vert */
            int third = -1;
            if (a != E[i].va && a != E[i].vb) third = a;
            else if (b != E[i].va && b != E[i].vb) third = b;
            else third = c;
            printf("  f%d=(%d,%d,%d) third=%d at (%.4f,%.4f,%.4f)\n",
                fi, a, b, c, third,
                V[third].z, V[third].y, V[third].x);
        }
    }
    printf("\nTotal non-manifold edges: %d\n", n_nm);
    return 0;
}
